#ifndef GF11_3D_ABS_SYM_HPP
#define GF11_3D_ABS_SYM_HPP
#include "find_root.hpp"
#include "throw_exception.hpp"
#include "tolerance.hpp"
#include <boost/math/constants/constants.hpp>
#include <boost/math/special_functions/erf.hpp>
#include <boost/math/tools/roots.hpp>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <cmath>
#include <cstdint>

namespace gf11
{

class GreensFunction3DAbsSym
{
  public:
    using real_type = double;

  public:
    GreensFunction3DAbsSym(const real_type D, const real_type a) noexcept
        : D_(D), a_(a)
    {}
    ~GreensFunction3DAbsSym() = default;

    real_type drawTime(const real_type rnd) const;
    real_type drawR   (const real_type rnd, const real_type t) const;

    real_type a() const noexcept {return a_;}
    real_type D() const noexcept {return D_;}

    static const char* name() noexcept {return "GreensFunction3DAbsSym";}

  private:

    real_type p_survival  (const real_type t) const;
    real_type p_r_fourier (const real_type r, const real_type t) const;
    real_type p_int_r     (const real_type r, const real_type t) const noexcept;
    real_type p_int_r_free(const real_type r, const real_type t) const noexcept;

    static real_type ellipticTheta4Zero(const real_type q);

  private:

    static constexpr std::uintmax_t max_iteration_count() noexcept {return 100;}

    static constexpr real_type CUTOFF()     noexcept {return 1e-10;}
    static constexpr real_type CUTOFF_H()   noexcept {return 6.0;}
    static           real_type log_CUTOFF() noexcept {return std::log(CUTOFF());}

    real_type D_;
    real_type a_;
};

inline GreensFunction3DAbsSym::real_type
GreensFunction3DAbsSym::p_survival(const real_type t) const
{
    const real_type a2(a_ * a_);
    const real_type pi2(boost::math::constants::pi_sqr<real_type>());
    const real_type q(-D_ * pi2 * t / a2);
    return real_type(1) - ellipticTheta4Zero(std::exp(q));
}

inline GreensFunction3DAbsSym::real_type
GreensFunction3DAbsSym::p_int_r_free(
        const real_type r, const real_type t) const noexcept
{
    const real_type Dt(this->D_ * t);
    const real_type sqrtDt(std::sqrt(Dt));
    const real_type sqrtPI(boost::math::constants::root_pi<real_type>());

    return boost::math::erf(r / (sqrtDt + sqrtDt)) -
           r * std::exp(- r * r / (4.0 * Dt)) / (sqrtPI * sqrtDt);
}

inline GreensFunction3DAbsSym::real_type
GreensFunction3DAbsSym::p_int_r(
        const real_type r, const real_type t) const noexcept
{
    if(std::abs(this->p_int_r_free(r, t)) < CUTOFF())
    {// p_int_r is always smaller than p_free.
        return 0.0;
    }

    real_type value(0.0);

    const real_type inva  = 1. / a_;
    const real_type PI2   = boost::math::constants::pi_sqr<real_type>();
    const real_type PIr   = boost::math::constants::pi<real_type>() * r;
    const real_type PIr_a = PIr * inva;
    const real_type Dt    = D_ * t;
    const real_type DtPIsq_asq = D_ * t * PI2 * inva * inva;

    const real_type maxn = (a_ / boost::math::constants::pi<real_type>()) *
                            std::sqrt((DtPIsq_asq - log_CUTOFF())/ Dt);
//  original: std::sqrt(std::log(std::exp(DtPIsq_asq) / CUTOFF) / Dt);

    const std::int64_t N_MAX = 10000;
    const std::int64_t N = (std::min)(
            static_cast<std::int64_t>(std::ceil(maxn) + 1), N_MAX);

    if (N == N_MAX)
    {
        std::cerr << "p_int_r: didn't converge\n";
    }

    for(std::int64_t n = 1; n <= N; ++n)
    {
        const real_type term1(std::exp(- n * n * DtPIsq_asq));
        const real_type angle_n(n * PIr_a);
        const real_type sin_n = std::sin(angle_n); // sincos
        const real_type cos_n = std::cos(angle_n);
        const real_type term2(a_ * sin_n);
        const real_type term3(n  * PIr * cos_n);

        const real_type term(term1 * (term2 - term3) / n);
        value += term;
    }

    const real_type factor(2 * inva / boost::math::constants::pi<real_type>());
    return value * factor;
}

inline GreensFunction3DAbsSym::real_type
GreensFunction3DAbsSym::p_r_fourier(
        const real_type r, const real_type t) const
{
    const real_type a2  = a_ * a_;
    const real_type r2  = r  * r;
    const real_type Dt  = D_ * t;
    const real_type pi2 = boost::math::constants::pi_sqr<real_type>();
    const real_type two_pi = boost::math::constants::two_pi<real_type>();

    real_type        value = 0.0;
    const std::int64_t N = 100;
    std::int64_t       n = 1;
    while(true)
    {
        const real_type term1(
                std::exp( -(pi2 * r2 + a2 * n*n) / (4.0 * pi2 * Dt)));

        // XXX the original version
//         const real_type term2(boost::math::constants::pi<real_type>() * r *
//                 std::exp(std::log(std::cosh(a_ * r * n / (two_pi * Dt)))));
//         const real_type term3(a_ * n *
//                 std::exp(std::log(std::sinh(a_ * r * n / (two_pi * Dt)))));

        const real_type term2(boost::math::constants::pi<real_type>() * r *
                std::cosh(a_ * r * n / (two_pi * Dt)));
        const real_type term3(a_ * n *
                std::sinh(a_ * r * n / (two_pi * Dt)));

        const real_type term(term1 * r * (term2 - term3));
        value += term;

        if(std::abs(value) * 1e-8 > std::abs(term)){break;}

        if(n > N)
        {
            std::cerr << "p_r_fourier: didn't converge\n";
            break;
        }
        ++n;
    }
    const real_type factor((boost::math::constants::root_two<real_type>() *
                            pi2 * std::pow(Dt, 1.5)));
    return value / factor;
}

inline GreensFunction3DAbsSym::real_type
GreensFunction3DAbsSym::drawTime(const real_type rnd) const
{
    if(rnd >= 1.0 || rnd < 0.0)
    {
        throw_exception<std::invalid_argument>(
                "GreensFunction3DAbsSym: 0.0 <= (rnd=", rnd, ") < 1.0");
    }
    if(D_ == 0.0 || a_ == std::numeric_limits<real_type>::infinity())
    {
        return std::numeric_limits<real_type>::infinity();
    }
    if (a_ == 0.0)
    {
        return 0.0;
    }

    const real_type t_guess(a_ * a_ / (6. * D_));
    real_type low  = t_guess;
    real_type high = t_guess;

    auto is_surviving = [rnd, this](const real_type t) -> real_type {
        return rnd - this->p_survival(t);
    };

    const real_type value(is_surviving(t_guess));
    real_type low_value  = value;
    real_type high_value = value;

    if(value < 0.0)
    {
        high *= 10.0;

        while(true)
        {
            high_value = is_surviving(high);
            if(high_value >= 0.0)
            {
                break;
            }
            if(std::abs(high) >= t_guess * 1e6)
            {
                throw std::runtime_error(
                        "GreensFunction3DAbsSym: couldn't adjust high");
            }
            high *= 10.0;
        }
    }
    else
    {
        real_type low_value_prev(value);
        low *= 0.1;

        while(true)
        {
            low_value = is_surviving(low);
            if(low_value <= 0.0)
            {
                break;
            }
            if(std::abs(low) <= t_guess * 1e-6 ||
               std::abs(low_value - low_value_prev) < CUTOFF())
            {
                std::cerr << "couldn't adjust low\n";
                return low;
            }
            low_value_prev = low_value;
            low *= .1;
        }
    }
    tolerance<real_type> tol(/*absolute =*/1e-18, /*relative =*/1e-12);
    return find_root(is_surviving, low, high, low_value, high_value, tol,
                 max_iteration_count(), "GreensFunction3DAbsSym::drawTime()");
}

inline GreensFunction3DAbsSym::real_type
GreensFunction3DAbsSym::drawR(
        const real_type rnd, const real_type t) const
{
    if (rnd >= 1.0 || rnd < 0.0)
    {
        throw_exception<std::invalid_argument>(
                "GreensFunction3DAbsSym: 0.0 <= (rnd=", rnd,") < 1.0", rnd);
    }
    if (t < 0.0)
    {
        throw_exception<std::invalid_argument>(
                "GreensFunction3DAbsSym: 0.0 < (t = ", t, ")");
    }
    if (a_ == 0.0 || t == 0.0 || D_ == 0.0)
    {
        return 0.0;
    }

    const real_type thresholdDistance(CUTOFF_H() * std::sqrt(6.0 * D_ * t));

    const real_type low(0.0);
    const real_type high(a_);
    tolerance<real_type> tol(/*absolute =*/1e-18, /*relative =*/1e-12);

    if (a_ <= thresholdDistance) // use r.
    {
        const real_type psurv = this->p_int_r(a_, t);
        if(psurv == 0.0)
        {
            return a_;
        }
        if(psurv < 0.0)
        {
            throw_exception<std::logic_error>(
                    "survival probability(", psurv, ") is negative");
        }
        const real_type target = psurv * rnd;
        auto diffuse_r = [target, t, this](const real_type r) -> real_type {
            return this->p_int_r(r, t) - target;
        };
        return find_root(diffuse_r, low, high, tol, max_iteration_count(),
                     "GreensFunction3DAbsSym::drawR()");
    }
    else // use r_free.
    {
        if (p_int_r_free(a_, t) < rnd)
        {
            std::cerr << "p_int_r_free(a, t) < rnd, returning a\n";
            return a_;
        }
        auto diffuse_r_free = [rnd, t, this](const real_type r) -> real_type {
            return this->p_int_r_free(r, t) - rnd;
        };

        return find_root(diffuse_r_free, low, high, tol, max_iteration_count(),
                     "GreensFunction3DAbsSym::drawR()");
    }
}

inline GreensFunction3DAbsSym::real_type
GreensFunction3DAbsSym::ellipticTheta4Zero(const real_type q)
{
    if (std::abs(q) > 1.0)
    {
        throw_exception<std::invalid_argument>(
            "GreensFunction3DAbsSym::ellipticTheta4Zero: "
            "abs(q) = (", q, ") <= 1.0");
    }

    const std::size_t N(1000);
    real_type value(1.0), q_n(q), q_2n(1.0);

    for(std::size_t n = 1; n <= N; ++n)
    {
        const real_type term2(1.0 - q_2n * q);  // q^(2n-1) = (q^(n-1))^2 * q

        q_2n = q_n * q_n;

        const real_type term1(1.0 - q_2n);
        const real_type term(term1 * term2 * term2);
        const real_type value_prev(value);
        value *= term;

        if (std::abs(value - value_prev) < 1e-18) {return value;}
        q_n *= q;
    }
    return value;
}

template<typename charT, typename traitsT>
inline std::basic_ostream<charT, traitsT>&
operator<<(std::basic_ostream<charT, traitsT>& os, const GreensFunction3DAbsSym& gf)
{
    os << "GreensFunction3DAbsSym("
       << "D=" << gf.D() << ", " << "a=" << gf.a() << ")";
    return os;
}

}// gf11
#endif//GF11_3D_ABS_SYM_HPP

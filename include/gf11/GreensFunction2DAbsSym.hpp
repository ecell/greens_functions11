#ifndef GF11_2D_ABS_SYM_HPP
#define GF11_2D_ABS_SYM_HPP
#include <boost/math/constants/constants.hpp>
#include <boost/math/special_functions/bessel.hpp>
#include <boost/math/special_functions/erf.hpp>
#include "throw_exception.hpp"
#include "find_root.hpp"
#include "tolerance.hpp"
#include <vector>
#include <stdexcept>
#include <iosfwd>
#include <limits>
#include <cstdint>
#include <cassert>
#include <cmath>

namespace gf11
{

template<typename realT>
class GreensFunction2DAbsSym
{
  public:
    using real_type = realT;

    GreensFunction2DAbsSym(const real_type D, const real_type a) noexcept
        : D_(D), a_(a)
    {}
    ~GreensFunction2DAbsSym() = default;

    real_type drawTime(const real_type rnd) const;
    real_type drawR   (const real_type rnd, const real_type t) const;

    real_type a() const noexcept {return a_;}
    real_type D() const noexcept {return D_;}

    static const char* name() noexcept {return "GreensFunction2DAbsSym";}

  private:

    real_type p_survival  (const real_type t)                    const;
    real_type p_int_r     (const real_type r, const real_type t) const;

  private:

    // H = 4.0: ~3e-5, 4.26: ~1e-6, 5.0: ~3e-7, 5.2: ~1e-7,
    //     5.6: ~1e-8, 6.0: ~1e-9
    static constexpr real_type CUTOFF   = 1e-10;
    static constexpr real_type CUTOFF_H = 6.0;

    real_type D_;
    real_type a_;
};

template<typename realT>
constexpr typename GreensFunction2DAbsSym<realT>::real_type
GreensFunction2DAbsSym<realT>::CUTOFF;
template<typename realT>
constexpr typename GreensFunction2DAbsSym<realT>::real_type
GreensFunction2DAbsSym<realT>::CUTOFF_H;

template<typename realT>
typename GreensFunction2DAbsSym<realT>::real_type
GreensFunction2DAbsSym<realT>::p_survival(const real_type t) const
{
    constexpr int N_max = 100;

    const real_type Dt = this->D_ * t;
    const real_type ra = real_type(1) / this->a_;

    real_type sum = 0;
    for(int n=1; n<=N_max; ++n)
    {
        const real_type aAn    = boost::math::cyl_bessel_j_zero<real_type>(0, n);
        const real_type J1_aAn = boost::math::cyl_bessel_j(1, aAn);
        const real_type An     = aAn * ra;
        const real_type term   = std::exp(-Dt * An * An) / (An * J1_aAn);
        sum += term;
        if(std::abs(term / sum) < CUTOFF)
        {
            break;
        }
    }
    return 2 * ra * sum;
}

template<typename realT>
typename GreensFunction2DAbsSym<realT>::real_type
GreensFunction2DAbsSym<realT>::p_int_r(const real_type r, const real_type t) const
{
    constexpr int N_max = 10000;
    const real_type Dt = this->D_ * t;
    const real_type ra = real_type(1) / this->a_;

    real_type sum = 0.0;
    for(int n=1; n <= N_max; ++n)
    {
        const real_type aAn = boost::math::cyl_bessel_j_zero<real_type>(0, n);
        const real_type An  = aAn * ra;
        const real_type rAn = r * An;
        const real_type J1_aAn = boost::math::cyl_bessel_j(1, aAn); // CylBesselGen
        const real_type J1_rAn = boost::math::cyl_bessel_j(1, rAn);
        const real_type term = (std::exp(-Dt * An * An) * r * J1_rAn) /
                               (An * J1_aAn * J1_aAn);
        sum += term;
        if(std::abs(term / sum) <= CUTOFF)
        {
            break;
        }
    }
    return 2.0 * ra * ra * sum;
}

template<typename realT>
typename GreensFunction2DAbsSym<realT>::real_type
GreensFunction2DAbsSym<realT>::drawTime(const real_type rnd) const
{
    if(!(0.0 <= rnd && rnd < 1.0))
    {
        throw_exception<std::invalid_argument>("GreensFunction2DRadAbs::drawTime: "
            "rnd(", rnd, ") must be in the range [0,1)");
    }
    if(this->D_ == 0 || this->a_ == std::numeric_limits<real_type>::infinity())
    {
        return std::numeric_limits<real_type>::infinity();
    }
    if(this->a_ == 0)
    {
        return 0;
    }

    // the equation to solve
    const auto p_surv_eq = [=](const real_type t) -> real_type {
        return 1 - this->p_survival(t) - rnd;
    };

    const real_type t_guess = this->a_ * this->a_ / (4 * this->D_);

    real_type value  = p_surv_eq(t_guess);
    real_type low_v  = value;
    real_type high_v = value;
    real_type low    = t_guess;
    real_type high   = t_guess;

    if(value < 0)
    {
        while(true)
        {
            high *= 10;
            value = p_surv_eq(high);

            if(std::abs(high) >= t_guess * 1e+6)
            {
                throw_exception<std::invalid_argument>("GreensFunction2DRadAbs::"
                    "drawTime: unable to find upper limit of time");
            }
            high_v = value;
            if(value > 0.0) {break;}
        }
    }
    else
    {
        while(true)
        {
            low *= 0.1;
            value = p_surv_eq(low);
            if(std::abs(low) <= t_guess * 1e-6 ||
               std::abs(value - low_v) < CUTOFF)
            {
                // couldn't adjust lower limit. it's almost zero.
                // XXX: here, low_v contains the previous value of `value`.
                //      so `value - low_v` means the difference from prev
                return low;
            }
            low_v = value;

            if(value < 0.0) {break;}
        }
    }
    const tolerance<real_type> tol(/*absolute =*/1e-18, /*relative =*/1e-12);
    return find_root(p_surv_eq, low, high, low_v, high_v, tol, 100,
                     "GreensFunction2DRadAbs::drawTime");
}

template<typename realT>
typename GreensFunction2DAbsSym<realT>::real_type
GreensFunction2DAbsSym<realT>::drawR(const real_type rnd, const real_type t) const
{
    if(!(0.0 <= rnd && rnd < 1.0))
    {
        throw_exception<std::invalid_argument>("GreensFunction2DRadAbs::drawR: "
            "rnd(", rnd, ") must be in the range [0,1)", rnd);
    }
    if (t < 0)
    {
        throw_exception<std::invalid_argument>("GreensFunction2DAbsSym::drawR: "
            "time must be positive: 0.0 < (t = ", t, ")");
    }
    if (this->a_ == 0 || t == 0 || this->D_ == 0)
    {
        return 0;
    }

    const real_type psurv  = this->p_survival(t);
    const real_type target = psurv * rnd;
    assert(0 < psurv);

    const auto p_r_eq = [=](const real_type r) -> real_type {
        return this->p_int_r(r, t) - target;
    };
    const tolerance<real_type> tol(/*absolute =*/1e-18, /*relative =*/1e-12);

    return find_root(p_r_eq, real_type(0), this->a_, tol, 100,
                     "GreensFunction2DAbsSym::drawR");
}

template<typename charT, typename traitsT, typename realT>
std::basic_ostream<charT, traitsT>&
operator<<(std::basic_ostream<charT, traitsT>& os,
           const GreensFunction2DAbsSym<realT>& gf)
{
    os << gf.name() << "(D=" << gf.D() << ",a=" << gf.a() << ")";
    return os;
}
} // gf11
#endif// GF11_2D_ABS_SYM_HPP

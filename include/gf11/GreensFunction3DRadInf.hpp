#ifndef GF11_3D_RAD_INF_H
#define GF11_3D_RAD_INF_H

#include "throw_exception.hpp"
#include "factorial.hpp"
#include "tolerance.hpp"
#include "sumup.hpp"
#include "find_root.hpp"
#include "gf_math.hpp"
#include "SphericalBesselGenerator.hpp"

#include <boost/math/tools/roots.hpp>
#include <boost/math/special_functions/expm1.hpp>
#include <boost/math/special_functions/erf.hpp>
#include <boost/format.hpp>

#include <gsl/gsl_sf_legendre.h>
#include <gsl/gsl_integration.h>

#include <array>
#include <iostream>
#include <limits>
#include <memory>
#include <stdexcept>
#include <utility>
#include <vector>

#include <cstdint>
#include <cmath>

namespace gf11
{

class GreensFunction3DRadInf
{
  public:
    using real_type = double;

  public:
    GreensFunction3DRadInf(const real_type D,  const real_type kf,
                           const real_type r0, const real_type sigma)
        : D_(D), kf_(kf), r0_(r0), sigma_(sigma),
          kD_(4.0 * boost::math::constants::pi<real_type>() * sigma * D),
          alpha_((1.0 + (kf / kD_)) * (std::sqrt(D) / sigma))
    {
        if(sigma > r0)
        {
            throw_exception<std::invalid_argument>("GreensFunction3DRadInf: "
                "r0(", r0, ") must larger than sigma(", sigma, ")");
        }
    }
    ~GreensFunction3DRadInf() = default;

    real_type D()     const noexcept {return this->D_;}
    real_type kf()    const noexcept {return this->kf_;}
    real_type kD()    const noexcept {return this->kD_;}
    real_type r0()    const noexcept {return this->r0_;}
    real_type sigma() const noexcept {return this->sigma_;}

    real_type drawTime (const real_type rnd) const;
    real_type drawR    (const real_type rnd, const real_type t) const;
    real_type drawTheta(const real_type rnd, const real_type r,
                        const real_type t  ) const;

    static const char* name() noexcept {return "GreensFunction3DRadInf";}

  private:

    real_type p_reaction(const real_type t) const;
    real_type p_survival(const real_type t) const
    {
        return 1.0 - this->p_reaction(t);
    }

    real_type p_int_r(const real_type r, const real_type t) const;

    real_type ip_free(const real_type theta, const real_type r,
                      const real_type t) const;
    real_type ip_corr(const real_type theta, const real_type r,
                      const std::vector<real_type>& Rn_table) const;

    // `t` is used to generate Rn_table, in greens_functions...
    real_type ip_theta(const real_type r, const real_type theta,
        const real_type t, const std::vector<real_type>& Rn_table) const;

    real_type p_free_max(const real_type r, const real_type t) const;
    real_type p_irr(const real_type r, const real_type t) const;

    real_type p_corr_R(const real_type alpha, const unsigned int n,
                       const real_type r, const real_type t) const;

    std::vector<real_type>
    make_Rn_table(const real_type r, const real_type t) const;

    struct p_corr_R_params // for GSL
    {
        const GreensFunction3DRadInf* gf;
        unsigned int n;
        real_type r;
        real_type t;
    };

  private:

    static constexpr real_type   TOLERANCE()       noexcept {return 1e-8;}
    static constexpr real_type   THETA_TOLERANCE() noexcept {return 1e-5;}
    static constexpr real_type   MIN_T()           noexcept {return 1e-12;}
    static constexpr real_type   H()               noexcept {return 4.0;}
    static constexpr std::size_t MAX_ORDER()       noexcept {return 70;}

    real_type D_;
    real_type kf_;
    real_type r0_;
    real_type sigma_;
    real_type kD_;
    real_type alpha_;
};

inline GreensFunction3DRadInf::real_type
GreensFunction3DRadInf::p_reaction(const real_type t) const
{
    const auto sqrt_t = std::sqrt(t);
    const auto sqrt_D = std::sqrt(D_);

    const auto r0_m_sigma_over_sqrt4D_t((r0_ - sigma_) / (2 * sqrt_D * sqrt_t));

    const auto Wf     = W(r0_m_sigma_over_sqrt4D_t, alpha_ * sqrt_t);
    const auto factor = sigma_ * kf_ / (r0_ * (kf_ + kD_));

    return factor * (boost::math::erfc(r0_m_sigma_over_sqrt4D_t) - Wf);
}

inline GreensFunction3DRadInf::real_type
GreensFunction3DRadInf::p_int_r(const real_type r, const real_type t) const
{
    const real_type Dt    = D_ * t;
    const real_type kf_kD = kf_ + kD_;
    const real_type Dt4   = 4.0 * Dt;

    const real_type rcp_sqrtDt4  = 1.0 / std::sqrt(Dt4);
    const real_type ksigma2      = 2.0 * kf_ * sigma_;
    const real_type alpha_sqrt_t = alpha_ * std::sqrt(t);

    const real_type r_r0_2s_sqrtDt4 = (r - 2.0 * sigma_ + r0_) * rcp_sqrtDt4;
    const real_type r_r0_sqrtDt4    = (r - r0_)                * rcp_sqrtDt4;
    const real_type r0_s_sqrtDt4    = (r0_ - sigma_)           * rcp_sqrtDt4;

    const real_type term1 =
        (boost::math::expm1(-r_r0_2s_sqrtDt4 * r_r0_2s_sqrtDt4) -
         boost::math::expm1(-r_r0_sqrtDt4    * r_r0_sqrtDt4)) *
        std::sqrt(Dt * 2 * boost::math::constants::one_div_two_pi<real_type>());

    const real_type erf_r_r0_2s_sqrtDt4 = boost::math::erf(r_r0_2s_sqrtDt4);

    const real_type term2 = kf_kD * r0_ * boost::math::erf(r_r0_sqrtDt4) +
        kf_kD * r0_ * erf_r_r0_2s_sqrtDt4 +
        ksigma2 * (boost::math::erf(r0_s_sqrtDt4) - erf_r_r0_2s_sqrtDt4);

    const real_type term3 =
         kf_ * sigma_ * W(r0_s_sqrtDt4, alpha_sqrt_t) -
        (kf_ * r + kD_ * (r - sigma_)) * W(r_r0_2s_sqrtDt4, alpha_sqrt_t);

    return (term1 + (1.0 / kf_kD)) * (0.5 * term2 + term3) / r0_;
}

inline GreensFunction3DRadInf::real_type
GreensFunction3DRadInf::ip_free(
        const real_type theta, const real_type r, const real_type t) const
{
    return ip_theta_free(theta, r, this->r0_, t, this->D_); // gf_math
}

inline GreensFunction3DRadInf::real_type
GreensFunction3DRadInf::ip_corr(const real_type theta, const real_type r,
        const std::vector<real_type>& Rn_table) const
{
    constexpr real_type pi = boost::math::constants::pi<real_type>();
    if(Rn_table.size() == 0)
    {
        return 0.0;
    }

    const real_type cos_theta = std::cos(theta);

    std::vector<real_type> legendre_table(Rn_table.size() + 2, 0.0);
    legendre_table[0] = 1.0; // n = -1
    gsl_sf_legendre_Pl_array(Rn_table.size(), cos_theta, std::addressof(legendre_table[1]));
    //XXX: not Rn_table.size() + 1? It looks weird.

    const real_type p = sumup_all(
        [this, &Rn_table, &legendre_table](const std::size_t n) -> real_type
        {
            const auto lgnd_n_m1(legendre_table[n]);   // n-1
            const auto lgnd_n_p1(legendre_table[n+2]); // n+1
            return Rn_table[n] * (lgnd_n_m1 - lgnd_n_p1);// / (1.0 + 2.0 * n);
        },
        Rn_table.size());

    return -p / (4.0 * pi * std::sqrt(r * r0_));
}

inline GreensFunction3DRadInf::real_type
GreensFunction3DRadInf::ip_theta(const real_type r,
    const real_type theta, const real_type t,
    const std::vector<real_type>& Rn_table) const
{
    return this->ip_free(theta, r, t) + this->ip_corr(theta, r, Rn_table);
}

inline GreensFunction3DRadInf::real_type
GreensFunction3DRadInf::p_irr(const real_type r, const real_type t) const
{
    constexpr real_type pi = boost::math::constants::pi<real_type>();
    const real_type sqrtD = std::sqrt(D_);

    const real_type alpha = (1.0 + (kf_ / kD_)) * (std::sqrt(D_) / sigma_);

    const real_type Dt4 = 4.0 * D_ * t;
    const real_type r_plus_r0_minus_2sigma = r + r0_ - 2.0 * sigma_;

    const auto square = [](const real_type x) noexcept -> real_type {return x*x;};

    const real_type num1 = std::exp(-square(r - r0_) / Dt4);
    const real_type num2 = std::exp(-square(r_plus_r0_minus_2sigma) / Dt4);
    const real_type num3 = W(r_plus_r0_minus_2sigma / std::sqrt(Dt4), alpha * std::sqrt(t));

    const real_type num      = (num1 + num2) / std::sqrt(4.0 * pi * t) - alpha * num3;
    const real_type den      = 4.0 * pi * r * r0_ * sqrtD;
    const real_type jacobian = 4.0 * pi * r * r;

    return num / den * jacobian;
}

inline GreensFunction3DRadInf::real_type
GreensFunction3DRadInf::p_free_max(const real_type r, const real_type t) const
{
    constexpr real_type pi = boost::math::constants::pi<real_type>();
    const real_type Dt4   = 4.0 * D_ * t;
    const real_type Dt4Pi = Dt4 * pi;
    const real_type dr_sq = (r - r0_) * (r - r0_);
    return std::exp(-dr_sq / Dt4) / std::sqrt(Dt4Pi * Dt4Pi * Dt4Pi);
}

inline GreensFunction3DRadInf::real_type
GreensFunction3DRadInf::p_corr_R(const real_type alpha, const unsigned int n,
                                 const real_type r, const real_type t) const
{
    constexpr real_type pi = boost::math::constants::pi<real_type>();

    const real_type ks     = kf_ * sigma_;
    const real_type realn  = static_cast<real_type>(n);
    const real_type ks_m_n = ks - realn;

    const real_type alphasq = alpha * alpha;

    const real_type term1 = std::exp(-D_ * t * alphasq);

    const real_type sAlpha  = sigma_ * alpha;
    const real_type rAlpha  = r      * alpha;
    const real_type r0Alpha = r0_    * alpha;

    SphericalBesselGenerator<real_type> s;

    const real_type js  = s.j(n,   sAlpha);
    const real_type ys  = s.y(n,   sAlpha);
    const real_type js1 = s.j(n+1, sAlpha);
    const real_type ys1 = s.y(n+1, sAlpha);
    const real_type jr  = s.j(n,   rAlpha);
    const real_type yr  = s.y(n,   rAlpha);
    const real_type jr0 = s.j(n,   r0Alpha);
    const real_type yr0 = s.y(n,   r0Alpha);

    const real_type R1 = (ks_m_n * js + sAlpha * js1);
    const real_type R2 = (ks_m_n * ys + sAlpha * ys1);

    const real_type F1R1 = R1 * jr * jr0 - R1 * yr * yr0;
    const real_type F2   = jr0 * yr + jr * yr0;

    const real_type num = 2.0 * std::sqrt(r * r0_) * alphasq * R1 * (F1R1 + F2 * R2);
    const real_type den = pi * (R1 * R1 + R2 * R2);

    const real_type result = term1 * num / den;

    assert(std::isfinite(result));

    return result;
}

inline std::vector<GreensFunction3DRadInf::real_type>
GreensFunction3DRadInf::make_Rn_table(const real_type r, const real_type t) const
{
    constexpr real_type pi = boost::math::constants::pi<real_type>();
    std::vector<real_type> Rn_table;
    {
        // First, estimate the size of p_corr, and if it's small enough,
        // we don't need to calculate it in the first place.

        const real_type pirr       = p_irr(r, t);
        const real_type ipfree_max = ip_free(pi, r, t) * 2 * pi * r * r;

        if(std::abs((pirr - ipfree_max) / ipfree_max) < 1e-8)
        {
            return Rn_table; // empty!
        }
    }

    const real_type pfreemax = p_free_max(r, t);

    std::unique_ptr<
        gsl_integration_workspace, decltype(&gsl_integration_workspace_free)
        > workspace(gsl_integration_workspace_alloc(2000),
                   &gsl_integration_workspace_free);

    real_type Rn_prev(0.0);
    const real_type Rn_factor = 1.0 / (4.0 * pi * std::sqrt(r * r0_));

    const real_type integration_tolerance = pfreemax / Rn_factor * THETA_TOLERANCE();
    const real_type truncation_tolerance  = pfreemax * THETA_TOLERANCE() * 1e-1;

    // function to be integrated
    auto p_corr_R_F = [](double alpha, void* p) noexcept -> double {
        const p_corr_R_params* para = (const p_corr_R_params*)p;
        return para->gf->p_corr_R(alpha, para->n, para->r, para->t);
    };

    for(unsigned int n=0; n<=MAX_ORDER(); ++n)
    {
        p_corr_R_params params = {this, n, r, t};

        gsl_function F = {p_corr_R_F, std::addressof(params)};

        const real_type umax = std::sqrt(40.0 / (this->D_ * t));
        real_type Rn;
        real_type error;

        gsl_integration_qag(&F, 0.0, umax, integration_tolerance, THETA_TOLERANCE(),
                            2000, GSL_INTEG_GAUSS61, workspace.get(),
                            std::addressof(Rn), std::addressof(error));

        Rn_table.push_back(Rn);

        // truncate when converged enough.
        const real_type abs_Rn = std::abs(Rn);
        if(abs_Rn * Rn_factor < truncation_tolerance && abs_Rn < Rn_prev)
        {
            break;
        }

        Rn_prev = abs_Rn;
    }
    return Rn_table;
}

inline GreensFunction3DRadInf::real_type
GreensFunction3DRadInf::drawTime(const real_type rnd) const
{
    if(rnd >= 1.0 || rnd < 0.0)
    {
        throw_exception<std::invalid_argument>(
                "GreensFunction3DRadInf: 0.0 <= (rnd = ", rnd, ") < 1.0");
    }
    if(rnd >= this->p_reaction(std::numeric_limits<real_type>::infinity()))
    {
        return std::numeric_limits<real_type>::infinity();
    }

    const tolerance<real_type> tol(/*absolute =*/1e-18, /*relative =*/1e-12);
    const std::uintmax_t max_iteration_count = 100;

    const real_type low  = std::numeric_limits<real_type>::min();
    const real_type high = 100.0; //TODO ??? I think we need to be a bit more careful.

    // because abs boundary locates far away, particle never escape from it.
    // it reacts, or not.
    const auto is_surviving = [rnd, this](const real_type t) -> real_type {
        return rnd - this->p_reaction(t);
    };

    return find_root(is_surviving, low, high, tol, max_iteration_count,
                     "GreensFunction3DRadInf::drawTime()");
}

inline GreensFunction3DRadInf::real_type
GreensFunction3DRadInf::drawR(const real_type rnd, const real_type t) const
{
    if (!(rnd < 1.0 && rnd >= 0.0))
    {
        throw_exception<std::invalid_argument>(
                "GreensFunction3DRadInf: 0.0 <= (rnd = ", rnd, ") < 1.0");
    }
    if (t < 0.0)
    {
        throw_exception<std::invalid_argument>(
                "GreensFunction3DRadInf: 0.0 <= (t = ", t, ")");
    }
    if(t == 0.0)
    {
        return r0_;
    }

    const real_type target = rnd * this->p_survival(t); // normalize
    const auto p_int_r_eq = [t, target, this](const real_type r) -> real_type {
        return this->p_int_r(r, t) - target;
    };

    // adjust low and high starting from r0.
    // this is necessary to avoid root finding in the long tails where
    // numerics can be unstable.

    real_type low  = r0_;
    real_type high = r0_;
    const real_type val_at_r0 = p_int_r_eq(r0_);

    real_type low_val  = val_at_r0;
    real_type high_val = val_at_r0;

    const real_type sqrt6Dt(std::sqrt(6.0 * D_ * t));
    if(val_at_r0 < 0.0) // update high
    {
        std::uint32_t h = 3;
        while(true)
        {
            high     = r0_ + h * sqrt6Dt;
            high_val = p_int_r_eq(high);
            if(high_val > 0.0)
            {
                break;
            }

            ++h;
            if(h > 20)
            {
                throw_exception<std::runtime_error>(
                    "GreensFunction3DRadInf::drawR: "
                    "h > 20 while adjusting upper bound of r");
            }
        }
    }
    else
    {
        std::uint32_t h = 3;
        while(true)
        {
            low = r0_ - h * sqrt6Dt;
            if(low < sigma_)
            {
                if(p_int_r_eq(sigma_) > 0.0)
                {
                    // WARN: maybe because of the numerical error,
                    //       p_int_r - target still positive at sigma.
                    //       but the result should be in [sigma, infty),
                    //       here it returns sigma
                    return sigma_;
                }
                low = sigma_;
                break;
            }
            low_val = p_int_r_eq(low);
            if(low_val < 0.0)
            {
                break;
            }
            ++h;
        }
    }

    const tolerance<real_type> tol(/* abs = */1e-15, /* rel = */TOLERANCE());
    const std::uintmax_t max_iteration_count = 100;

    return find_root(p_int_r_eq, low, high, low_val, high_val, tol,
                     max_iteration_count, "GreensFunction3DRadInf::drawR()");
}

inline GreensFunction3DRadInf::real_type
GreensFunction3DRadInf::drawTheta(const real_type rnd,
        const real_type r, const real_type t) const
{
    if (!(rnd < 1.0 && rnd >= 0.0))
    {
        throw_exception<std::invalid_argument>(
            "GreensFunction3DRadInf::drawTheta: 0 <= (rnd = ", rnd, ") < 1");
    }
    if (r < sigma_)
    {
        throw_exception<std::invalid_argument>(
            "GreensFunction3DRadInf::drawTheta: sigma <= (r = ", r, ")");
    }
    if (t < 0.0)
    {
        throw_exception<std::invalid_argument>(
            "GreensFunction3DRadInf::drawTheta: 0 <= (t = ", t, ")");
    }
    if(t == 0.0)
    {
        return 0.0;
    }

    constexpr real_type pi = boost::math::constants::pi<real_type>();

    const auto Rn_table = this->make_Rn_table(r, t);

    // value at pi
    const real_type target = this->ip_theta(pi, r, t, Rn_table) * rnd;

    const auto ip_theta_eq =
        [&Rn_table, r, t, target, this](const real_type theta) -> real_type {
        return this->ip_theta(theta, r, t, Rn_table);
    };

    const tolerance<real_type> tol(/*abs =*/1e-15, /*rel =*/THETA_TOLERANCE());
    const std::uintmax_t max_iteration_count = 100;

    return find_root(ip_theta_eq, 0.0, pi, tol, max_iteration_count,
                     "GreensFunction3DRadInf::drawTheta()");
}

template<typename charT, typename traitsT>
std::basic_ostream<charT, traitsT>&
operator<<(std::basic_ostream<charT, traitsT>& os, const GreensFunction3DRadInf& gf)
{
    os << boost::format("%1%(D=%2%,r0=%3%,sigma=%4%,kf=%5%,kD=%6%)") %
          gf.name() % gf.D() % gf.r0() % gf.sigma() % gf.kf() % gf.kD();
    return os;
}

}// gf11
#endif// GF11_3D_RAD_INF

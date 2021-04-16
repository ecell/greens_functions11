#ifndef GF11_3D_RAD_INF_HPP
#define GF11_3D_RAD_INF_HPP
#include "throw_exception.hpp"
#include <boost/math/constants/constants.hpp>
#include <boost/format.hpp>
#include <iosfwd>
#include <stdexcept>
#include <vector>
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

template<typename charT, typename traitsT>
std::basic_ostream<charT, traitsT>&
operator<<(std::basic_ostream<charT, traitsT>& os, const GreensFunction3DRadInf& gf)
{
    os << boost::format("%1%(D=%2%,r0=%3%,sigma=%4%,kf=%5%,kD=%6%)") %
          gf.name() % gf.D() % gf.r0() % gf.sigma() % gf.kf() % gf.kD();
    return os;
}

}// gf11

#if defined(GF11_HEADER_ONLY)
#include "GreensFunction3DRadInf.cpp"
#endif

#endif// GF11_3D_RAD_INF_HPP

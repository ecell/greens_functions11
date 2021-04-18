#ifndef GF11_3D_RAD_ABS_HPP
#define GF11_3D_RAD_ABS_HPP

#include "throw_exception.hpp"
#include "ipv_event_kind.hpp"
#include "factorial.hpp"

#include <boost/format.hpp>
#include <boost/math/constants/constants.hpp>

#include <array>
#include <iosfwd>
#include <stdexcept>
#include <utility>
#include <vector>

#include <cstdint>

namespace gf11
{

class GreensFunction3DRadAbs
{
  public:
    using real_type = double;

    // XXX This name of enum is too confusing to use, but is kept for backward
    //     compatibility (originally, it was a normal enum defined in a base
    //     class, `GreensFunction`.)
    using event_type = GreensFunction;

  public:
    GreensFunction3DRadAbs(
            const real_type D,     const real_type kf, const real_type r0,
            const real_type sigma, const real_type a)
        : D_(D), kf_(kf), r0_(r0), sigma_(sigma), a_(a),
          h_(kf / (4 * boost::math::constants::pi<real_type>() * sigma * sigma * D)),
          h_sigma_plus_1(h_ * sigma + 1.0)
    {
        if(a < sigma)
        {
            throw_exception<std::invalid_argument>("GreensFunction3DRadAbs: "
                "a(", a, ") must be larger than sigma(", sigma, ")");
        }
        if (!(sigma <= r0 && r0 <= a))
        {
            throw_exception<std::invalid_argument>("GreensFunction3DRadAbs: "
                "r0(", r0, ") must be in the range [sigma(", sigma, "), a(", a, ")]");
        }

        alpha_offset_table.fill(0xFFFFFFFF);
        alpha_offset_table[0] = 0;
        alpha_table.fill(std::vector<real_type>{});
    }
    ~GreensFunction3DRadAbs() = default;

    real_type D()     const noexcept {return this->D_;}
    real_type kf()    const noexcept {return this->kf_;}
    real_type r0()    const noexcept {return this->r0_;}
    real_type a()     const noexcept {return this->a_;}
    real_type h()     const noexcept {return this->h_;}
    real_type sigma() const noexcept {return this->sigma_;}

    real_type  drawTime     (const real_type rnd) const;
    event_type drawEventType(const real_type rnd, const real_type t) const;
    real_type  drawR        (const real_type rnd, const real_type t) const;
    real_type  drawTheta    (const real_type rnd, const real_type r, const real_type t) const;

    static const char* name() noexcept {return "GreensFunction3DRadAbs";}

  private: // member functions

    real_type f_alpha0_aux(const real_type alpha) const;
    real_type alpha0_i    (const std::uint32_t i) const;
    real_type get_alpha0  (const std::uint32_t i) const;

    real_type num_r0(const real_type alpha) const;

    real_type leave_s_i(const real_type alpha) const;
    real_type leave_a_i(const real_type alpha) const;
    real_type leave_s  (const real_type t)     const;
    real_type leave_a  (const real_type t)     const;

    std::uint32_t guess_max_i(const real_type t) const;
    real_type p_survival_irr         (const real_type t) const;
    real_type p_survival_nocollision (const real_type t) const;
    void      create_p_survival_table(std::vector<real_type>& tab) const;
    real_type p_survival_i           (const real_type t) const;
    real_type p_survival(const real_type t, std::vector<real_type>& tab) const;

    real_type p_int_r_i(
            const real_type r, const real_type t, const real_type num_r0) const;
    real_type p_int_r(const real_type r, const real_type t) const;

    std::uint32_t alpha_offset(const std::size_t n) const;

    real_type f_alpha_aux(const real_type alpha, const std::uint32_t n) const;
    real_type f_alpha    (const real_type alpha, const std::uint32_t n) const;
    real_type alpha_i    (const std::uint32_t i, const std::uint32_t n) const;
    real_type get_alpha  (const std::size_t   n, const std::uint32_t i) const;

    real_type p_n_alpha(const std::uint32_t i, const std::uint32_t n,
                        const real_type r,     const real_type t) const;
    real_type p_n      (const std::size_t n, const real_type r,
                        const real_type t,   const real_type max_alpha) const;
    std::vector<real_type>
    make_p_n_table(const real_type r, const real_type t) const;

    real_type dp_n_alpha_at_a(const std::uint32_t i, const std::uint32_t n,
                              const real_type t) const;
    real_type dp_n_at_a(const std::uint32_t n, const real_type t,
                        const real_type max_alpha) const;
    std::vector<real_type>
    make_dp_n_at_a_table(const real_type r) const;

    real_type ip_theta(const real_type theta, const real_type r,
                       const std::vector<real_type>& p_n_table) const;

    static real_type G(const std::uint32_t n, const std::uint32_t k)
    {
        return factorial<real_type>(n + k) *
              (factorial_r<real_type>(k) * factorial_r<real_type>(n - k));
    }

    static real_type P(const std::uint32_t n, const real_type x)
    {
        const real_type x2sq_r = 1.0 / ((x + x) * (x + x));

        real_type    sx2    = 1.0;
        std::int32_t sign   = 1;
        real_type    result = 0.0;
        for(std::uint32_t m = 0, maxm = n / 2; m <= maxm; ++m)
        {
            result += sign * sx2 * G(n, 2 * m);

            sign *= -1;
            sx2  *= x2sq_r;
        }

        return result;
    }

    static std::pair<real_type, real_type>
    P2(const std::uint32_t n, const real_type x)
    {
        const real_type     x2sq_r = 1.0 / ((x + x) * (x + x));
        const std::uint32_t np1    = n + 1;

        real_type    sx2  = 1.0;
        std::int32_t sign = 1;

        real_type result  = 0.0;
        real_type resultp = 0.0;
        for(std::uint32_t m = 0, maxm = n / 2; m <= maxm; ++m)
        {
            const real_type     sx2p = sign * sx2;
            const std::uint32_t m2   = 2 * m;

            result  += sx2p * G(n, m2);
            resultp += sx2p * G(np1, m2);

            sign *= -1;
            sx2  *= x2sq_r;
        }
        if(n % 2)
        {
            resultp += sign * sx2 * G(np1, np1);
        }
        return std::make_pair(result, resultp);
    }

    static real_type Q(const std::uint32_t n, const real_type x)
    {
        real_type       sx2  = 1.0 / (x + x);
        const real_type x2sq = sx2 * sx2;

        // sum_(0)^((n-1)/2)
        real_type    result = 0.0;
        std::int32_t sign   = 1;
        for(std::uint32_t m = 0, maxm = (n + 1) / 2; m < maxm; ++m)
        {
            result += sign * sx2 * G(n, 2 * m + 1);
            sign *= -1;
            sx2  *= x2sq;
        }
        return result;
    }

    static std::pair<real_type, real_type>
    Q2(const std::uint32_t n, const real_type x)
    {
        real_type           sx2  = 1.0 / (x + x);
        const real_type     x2sq = sx2 * sx2;
        const std::uint32_t np1  = n + 1;

        std::int32_t sign = 1;
        real_type result  = 0.0;
        real_type resultp = 0.0;
        for(std::uint32_t m = 0, maxm = (n + 1) / 2; m < maxm; ++m)
        {
            const real_type     sx2p = sign * sx2;
            const std::uint32_t m2p1 = 2 * m + 1;

            result  += sx2p * G(n,   m2p1);
            resultp += sx2p * G(np1, m2p1);

            sign *= -1;
            sx2  *= x2sq;
        }
        if (!(n % 2))
        {
            resultp += sign * sx2 * G(np1, np1);
        }
        return std::make_pair(result, resultp);
    }

  private: // member variables

    // Error tolerance used by default.
    static constexpr real_type   TOLERANCE    () noexcept {return 1e-8;}
    static constexpr std::size_t MAX_ITERATION() noexcept {return 100;}
    // SphericalBesselGenerator's accuracy,
    // used by some theta-related calculations.
    static constexpr real_type THETA_TOLERANCE() noexcept {return 1e-5;}
    static constexpr real_type MIN_T_FACTOR   () noexcept {return 1e-8;}
    // for lookup-table
    static constexpr std::size_t MAX_ORDER    () noexcept {return 50;}
    static constexpr std::size_t MAX_ALPHA_SEQ() noexcept {return 2000;}

    real_type D_;
    real_type kf_;
    real_type r0_;
    real_type sigma_;
    real_type a_;
    real_type h_;
    real_type h_sigma_plus_1;

    // lookup tables in range [0, MAX_ORDER]
    mutable std::array<std::uint32_t,          51> alpha_offset_table;
    mutable std::array<std::vector<real_type>, 51> alpha_table;
};

template<typename charT, typename traitsT>
std::basic_ostream<charT, traitsT>&
operator<<(std::basic_ostream<charT, traitsT>& os, const GreensFunction3DRadAbs& gf)
{
    os << boost::format("%1%(D=%2%,r0=%3%,sigma=%4%,a=%5%,kf=%6%,h=%7%)") %
          gf.name() % gf.D() % gf.r0() % gf.sigma() % gf.a() % gf.kf() % gf.h();
    return os;
}

}// gf11

#if defined(GF11_HEADER_ONLY)
#include "GreensFunction3DRadAbs.cpp"
#endif


#endif// GF11_3D_RAD_ABS_HPP

#ifndef GF11_2D_RAD_ABS_HPP
#define GF11_2D_RAD_ABS_HPP
#include "ipv_event_kind.hpp"
#include "throw_exception.hpp"

#include <boost/math/constants/constants.hpp>

#include <algorithm>
#include <array>
#include <iosfwd>
#include <stdexcept>
#include <tuple>
#include <utility>
#include <vector>

#include <cstdint>

namespace gf11
{

class GreensFunction2DRadAbs
{
  public:
    using real_type = double;

    // XXX This name of enum is too confusing to use, but is kept for backward
    //     compatibility (originally, it was a normal enum defined in a base
    //     class, `GreensFunction`.)
    using event_type = GreensFunction;

  public:

    GreensFunction2DRadAbs(const real_type D,  const real_type kf,
                           const real_type r0, const real_type sigma,
                           const real_type a) noexcept
        : D_(D), kf_(kf), r0_(r0), sigma_(sigma), a_(a),
          h_(kf / (D * 2 * boost::math::constants::pi<real_type>() * sigma)),
          estimated_alpha_root_distance_(
                  boost::math::constants::pi<real_type>() / (a - sigma))
    {
        if(a < sigma)
        {
            throw_exception<std::invalid_argument>("GreensFunction2DRadAbs: "
                "a(", a, ") must be larger than sigma(", sigma, ")");
        }
        if (!(sigma <= r0 && r0 <= a))
        {
            throw_exception<std::invalid_argument>("GreensFunction2DRadAbs: "
                "r0(", r0, ") must be in the range [sigma(", sigma, "), a(", a, ")]");
        }

        // Clears all vectors in the alpha_table
        for(auto& vec : this->alpha_table_)
        {
            vec.clear();
        }

        // Sets all values of the alpha_x_scan_table_ to zero.
        std::fill(this->alpha_x_scan_table_.begin(),
                  this->alpha_x_scan_table_.end(),
                  SCAN_START() * estimated_alpha_root_distance_);

        // Sets all values of the alpha_correctly_estimated_ table to zero.
        std::fill(this->alpha_correctly_estimated_.begin(),
                  this->alpha_correctly_estimated_.end(),
                  0);
    }
    ~GreensFunction2DRadAbs() = default;

    real_type  drawTime     (const real_type rnd) const;
    event_type drawEventType(const real_type rnd, const real_type t) const;
    real_type  drawR        (const real_type rnd, const real_type t) const;
    real_type  drawTheta    (const real_type rnd, const real_type r, const real_type t) const;

    real_type D()     const noexcept {return D_;}
    real_type kf()    const noexcept {return kf_;}
    real_type r0()    const noexcept {return r0_;}
    real_type sigma() const noexcept {return sigma_;}
    real_type a()     const noexcept {return a_;}
    real_type h()     const noexcept {return h_;}

    static const char* name() noexcept {return "GreensFunction2DRadAbs";}

  private:

    // ------------------------------------------------------------------------

    real_type leaves      (const real_type t) const;
    real_type leaves_i_exp(const std::uint32_t i, const real_type t) const;
    real_type leaves_i    (const real_type alpha) const;

    real_type leavea      (const real_type t) const;
    real_type leavea_i_exp(const std::uint32_t i, const real_type t) const;

    real_type calc_A_i_0(const real_type alpha) const;

    // ------------------------------------------------------------------------

    real_type f_alpha0(const real_type alpha) const;
    real_type f_alpha (const real_type alpha, const std::uint32_t n) const;

    real_type get_alpha(const std::size_t n, const std::size_t i) const;

    real_type get_alpha_root_0(const real_type low, const real_type high) const;
    real_type get_alpha_root_N(const real_type low, const real_type high, const std::size_t n) const;
    real_type get_alpha_root  (const real_type low, const real_type high, const std::size_t n) const;

    std::pair<real_type, real_type> give_root_interval       (const std::size_t n) const;
    std::pair<real_type, real_type> give_root_interval_simple(const std::size_t n, const std::size_t j) const;

    void decide_on_method2(const std::size_t n, const std::size_t j) const;

    // ------------------------------------------------------------------------
    // drawTime related funcs

    real_type p_survival(const real_type t) const;
    real_type p_survival_table(const real_type t, std::vector<real_type>& table) const;
    real_type p_survival_i(const real_type alpha) const;
    real_type p_survival_i_exp_table(const std::uint32_t i, const real_type t,
                                     std::vector<real_type>& psurv_table) const;
    void create_psurv_table(std::vector<real_type>& table) const;

    std::size_t guess_maxi(const real_type t) const;

    // ------------------------------------------------------------------------
    // drawR related funcs

    void create_Y0J0_tables(std::vector<real_type>& Y0_table,
                            std::vector<real_type>& J0_table,
                            std::vector<real_type>& Y0J1J0Y1_table,
                            const real_type t) const;

    std::tuple<real_type, real_type, real_type>
    Y0J0J1_constants(const real_type alpha, const real_type t) const;

    real_type p_int_r_table(const real_type r,
                const std::vector<real_type>& Y0_aAn_table,
                const std::vector<real_type>& J0_aAn_table,
                const std::vector<real_type>& Y0J1J0Y1_table) const;

    real_type p_int_r_i_exp_table(const std::uint32_t i, const real_type r,
                const std::vector<real_type>& Y0_aAn_table,
                const std::vector<real_type>& J0_aAn_table,
                const std::vector<real_type>& Y0J1J0Y1_table) const;

    // ------------------------------------------------------------------------
    // drawTheta related funcs

    real_type ip_theta_table(const real_type theta,
                             const std::vector<real_type>& p_n_table) const;
    real_type ip_theta_n(const std::uint32_t n, const real_type theta,
                         const std::vector<real_type>& p_n_table) const;

    real_type p_m_alpha(const std::uint32_t n, const std::uint32_t m,
                        const real_type r, const real_type t) const;
    real_type p_m(const std::uint32_t m,
                  const real_type r, const real_type t) const;

    void make_p_m_table(std::vector<real_type>& p_m_table,
                        const real_type r, const real_type t) const;

    real_type dp_m_at_a(const std::uint32_t m, const real_type t) const;
    real_type dp_m_alpha_at_a(const std::uint32_t n, const std::uint32_t m,
                              const real_type t) const;

    void make_dp_m_at_a_table(std::vector<real_type>& p_m_table,
                              const real_type t) const;

  private:

    real_type D_;
    real_type kf_;
    real_type r0_;
    real_type sigma_;
    real_type a_;
    real_type h_;
    real_type estimated_alpha_root_distance_;

    // MAX_ORDER+1
    mutable std::array<std::vector<real_type>, 31> alpha_table_;
    mutable std::array<real_type,              31> alpha_x_scan_table_;
    mutable std::array<std::uint32_t,          31> alpha_correctly_estimated_;

  private:
    // Error tolerance used by default.
    static constexpr real_type TOLERANCE() {return 1e-8;}

    // SphericalBesselGenerator's accuracy, used by some
    // theta-related calculations.

    static constexpr real_type MIN_T_FACTOR() {return 1e-8;}

    static constexpr real_type L_TYPICAL() {return 1e-7;} // typical length scale
    static constexpr real_type T_TYPICAL() {return 1e-5;} // typical time scale
    static constexpr real_type EPSILON()   {return 1e-12;} // relative numeric error

    static constexpr unsigned int MAX_ORDER()     {return 30;}  // The maximum number of m terms
    static constexpr unsigned int MAX_ALPHA_SEQ() {return 500;} // The maximum number of n terms

    // Parameters for alpha-root finding
    // ======
    // See getAlpha() in cpp file for more information.
    //
    // Parameters for scanning method
    // Left boundary of 1st search interval 1st root
    static constexpr real_type SCAN_START() {return 0.001;}
    // Length of the scanning interval relative to estimated interval
    static constexpr real_type FRACTION_SCAN_INTERVAL() {return 0.5;} // TODO CHANGED THIS FROM .5 to .2

    // Other paramters
    // After CONVERGENCE_ASSUMED subsequent roots that lay within +/-
    // INTERVAL_MARGIN from the distance to which the distance is known to
    // converge, it is assumed all following roots have a distances inbetween
    // that don't deviate for more than INTERVAL_MARGIN from the distance to
    // which the roots are known to converge (Pi/(a-sigma)).
    static constexpr real_type CONVERGENCE_ASSUMED() {return 25;}
    static constexpr real_type INTERVAL_MARGIN()     {return 0.33;}
};

template<typename charT, typename traitsT>
std::basic_ostream<charT, traitsT>&
operator<<(std::basic_ostream<charT, traitsT>& os, const GreensFunction2DRadAbs& gf)
{
    os << gf.name()
       << ": D     = " << gf.D()
       << ", kf    = " << gf.kf()
       << ", r0    = " << gf.r0()
       << ", sigma = " << gf.sigma()
       << ", a     = " << gf.a()
       << ", h     = " << gf.h()
       ;
    return os;
}

} // gf11

#if defined(GF11_HEADER_ONLY)
#include "GreensFunction2DRadAbs.cpp"
#endif

#endif// GF11_2D_RAD_ABS_HPP

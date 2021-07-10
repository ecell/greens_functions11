#ifndef GF11_HEADER_ONLY
#include "GreensFunction2DRadAbs.hpp"
#endif

#include "config.hpp"
#include "find_root.hpp"
#include "ipv_event_kind.hpp"
#include "sumup.hpp"
#include "tolerance.hpp"

#include <boost/math/special_functions/bessel.hpp>
#include <boost/math/special_functions/lambert_w.hpp>
#include <boost/format.hpp>

#include <limits>
#include <cassert>
#include <cmath>

namespace gf11
{
GF11_INLINE GreensFunction2DRadAbs::real_type
GreensFunction2DRadAbs::leaves(const real_type t) const
{
    // calculates the flux leaving through the inner interface at a given moment
    // FIXME: This is inaccurate for small t's!!

    constexpr real_type pi_sq = boost::math::constants::pi_sqr<real_type>();

    const real_type p = sumup_until_convergence(
            [this, t](const std::uint32_t i) noexcept -> real_type {
                return this->leaves_i_exp(i, t);
            }, MAX_ALPHA_SEQ(), 1e-7);

    // The minus is not there because the flux is in the negative r
    // direction, and the direction is already accounted for in the derivative of B0,n(r)
    // See also leaves_i
    return (pi_sq / 2) * D_ * sigma_ * p;
}

GF11_INLINE GreensFunction2DRadAbs::real_type
GreensFunction2DRadAbs::leaves_i_exp(const std::uint32_t i, const real_type t) const
{
    // adds the exponential with the time to the sum.
    // Needed for the inner interface (reaction)

    const real_type alpha = this->get_alpha(0, i);

    return std::exp(-D_ * t * alpha * alpha) * this->leaves_i(alpha);
}

GF11_INLINE GreensFunction2DRadAbs::real_type
GreensFunction2DRadAbs::leaves_i(const real_type alpha) const
{
    namespace policies = boost::math::policies;
    constexpr auto policy = policies::make_policy(policies::promote_double<false>());

    // Calculates the n-th term of the summation for calculating the flux through
    // the inner interface (reaction)

    const real_type s_An = sigma_ * alpha;
    const real_type a_An = a_     * alpha;

    const real_type J1_sAn = boost::math::cyl_bessel_j(1, s_An, policy);
    const real_type Y1_sAn = boost::math::cyl_neumann (1, s_An, policy);
    const real_type J0_aAn = boost::math::cyl_bessel_j(0, a_An, policy);
    const real_type Y0_aAn = boost::math::cyl_neumann (0, a_An, policy);

    // the original comment says this is either An,0 or A0,n; TODO: write a doc
    const real_type A_i_0  = this->calc_A_i_0(alpha);

    // dBn,0(sigma)/dr
    const real_type dB_n_0dr = -alpha * (J1_sAn * Y0_aAn - Y1_sAn * J0_aAn);

    return A_i_0 * dB_n_0dr;
}

GF11_INLINE GreensFunction2DRadAbs::real_type
GreensFunction2DRadAbs::leavea(const real_type t) const
{
    // calculates the flux leaving through the outer interface at a given moment

    constexpr real_type pi = boost::math::constants::pi<real_type>();

    const real_type p = sumup_until_convergence(
            [this, t](const std::uint32_t i) noexcept -> real_type {
                return this->leavea_i_exp(i, t);
            }, MAX_ALPHA_SEQ(), 1e-7);

    return pi * (this->D_) * p;
}

GF11_INLINE GreensFunction2DRadAbs::real_type
GreensFunction2DRadAbs::leavea_i_exp(const std::uint32_t i, const real_type t) const
{
    // adds the exponential with the time to the sum.
    // Needed for the calculation of the flux throught the outer interface

    const real_type alpha = this->get_alpha(0, i);
    return std::exp(-this->D_ * t * alpha * alpha) * this->calc_A_i_0(alpha);
}

GF11_INLINE GreensFunction2DRadAbs::real_type
GreensFunction2DRadAbs::calc_A_i_0(const real_type alpha) const
{
    namespace policies = boost::math::policies;
    constexpr auto policy = policies::make_policy(policies::promote_double<false>());
    // Calculates the factor An,0 for (for example) determination of the flux
    // through the outer interface

    const real_type s_An = sigma_ * alpha;
    const real_type a_An = a_     * alpha;
    const real_type r0An = r0_    * alpha;


    const real_type J0_sAn  = boost::math::cyl_bessel_j(0, s_An, policy);
    const real_type J1_sAn  = boost::math::cyl_bessel_j(1, s_An, policy);

    const real_type J0_aAn  = boost::math::cyl_bessel_j(0, a_An, policy);
    const real_type Y0_aAn  = boost::math::cyl_neumann (0, a_An, policy);

    const real_type J0_r0An = boost::math::cyl_bessel_j(0, r0An, policy);
    const real_type Y0_r0An = boost::math::cyl_neumann (0, r0An, policy);


    const real_type alpha_sq = alpha * alpha;

    const real_type rho    = h_ * J0_sAn + alpha * J1_sAn;
    const real_type rho_sq = rho * rho;

    const real_type B_n_0  = J0_r0An * Y0_aAn - Y0_r0An * J0_aAn;

    const real_type A_i_0  = (alpha_sq * rho_sq * B_n_0) /
                             (rho_sq - (J0_aAn * J0_aAn) * (h_ * h_ + alpha_sq));
    return A_i_0;
}

// -----------------------------------------------------------------------------

GF11_INLINE GreensFunction2DRadAbs::real_type
GreensFunction2DRadAbs::f_alpha0(const real_type alpha) const
{
    namespace policies = boost::math::policies;
    constexpr auto policy = policies::make_policy(policies::promote_double<false>());
    // The method evaluates the equation for finding the alphas for given alpha.
    // This is needed to find the alpha's at which the expression is zero ->
    // alpha is the root.
    const real_type s_An = this->sigma_ * alpha;
    const real_type a_An = this->a_     * alpha;

    const real_type J0_s_An = boost::math::cyl_bessel_j(0, s_An, policy);
    const real_type J1_s_An = boost::math::cyl_bessel_j(1, s_An, policy);
    const real_type J0_a_An = boost::math::cyl_bessel_j(0, a_An, policy);

    const real_type Y0_s_An = boost::math::cyl_neumann (0, s_An, policy);
    const real_type Y1_s_An = boost::math::cyl_neumann (1, s_An, policy);
    const real_type Y0_a_An = boost::math::cyl_neumann (0, a_An, policy);

    const real_type rho1 = (this->h_ * J0_s_An + alpha * J1_s_An) * Y0_a_An;
    const real_type rho2 = (this->h_ * Y0_s_An + alpha * Y1_s_An) * J0_a_An;

    return rho1 - rho2;
}

GF11_INLINE GreensFunction2DRadAbs::real_type
GreensFunction2DRadAbs::f_alpha(
        const real_type alpha, const std::uint32_t n) const
{
    namespace policies = boost::math::policies;
    constexpr auto policy = policies::make_policy(policies::promote_double<false>());

    // The roots (y=0) of this function are constants in the Green's Functions.
    const real_type s_An    = sigma_ * alpha;
    const real_type a_An    = a_     * alpha;
    const real_type h_sigma = h_ * sigma_;

    const real_type Jn_s_An  = boost::math::cyl_bessel_j(n,   s_An, policy);
    const real_type Jn1_s_An = boost::math::cyl_bessel_j(n+1, s_An, policy);
    const real_type Jn_a_An  = boost::math::cyl_bessel_j(n,   a_An, policy);

    const real_type Yn_s_An  = boost::math::cyl_neumann (n,   s_An, policy);
    const real_type Yn1_s_An = boost::math::cyl_neumann (n+1, s_An, policy);
    const real_type Yn_a_An  = boost::math::cyl_neumann (n,   a_An, policy);

    const real_type rho1 = (h_sigma * Jn_s_An + s_An * Jn1_s_An - n * Jn_s_An) * Yn_a_An;
    const real_type rho2 = (h_sigma * Yn_s_An + s_An * Yn1_s_An - n * Yn_s_An) * Jn_a_An;

    return rho1 - rho2;
}

GF11_INLINE GreensFunction2DRadAbs::real_type
GreensFunction2DRadAbs::get_alpha(const std::size_t n, /* order of the Bessel fn */
                                  const std::size_t i  /* ith root */) const
{
    // This function searches for roots (y=0) of the so-called alpha-function
    // (::f_alpha()). It either moves a small search interval along the x-axis
    // to check for sign-change (which would indicate a root), and calls the
    // root finder, or directly calls the root finder if the spacing between
    // roots is found to be converged.

    real_type low, high;

    // n is the order of the Bessel function.
    assert(n < this->alpha_table_.size());
    auto& alpha_table = this->alpha_table_[n];

    // # Expansion of root table

    // If doesn't contain requested value, expand table until value
    const auto old_size = alpha_table.size();
    if(old_size <= i)
    {
        // Expand the table, temporarily fill with zeroes
        alpha_table.resize(i+1, 0.0);

        // # Calculating the necessary roots to expand the table
        for(std::size_t j = old_size; j <= i; j++)
        {
            if (alpha_table[j] != 0) // the root is already found
            {
                std::cerr << "tried accessing root that's not 0. Didn't search..\n";
                std::cerr << boost::format("    i = %1%, old_size = %2%, j = %3%\n") %
                                           i % old_size % j;
            }
            else
            {
                // Method 1. SCANNING. If the roots are not expected to lie close enough
                // to the estimate, use the "scanning" procedure to find an interval
                // that contains a root. (More robust method.)
                //      If it is established that we can use method 2,
                // alpha_x_scan_table_[n] will contain a value < 0.
                if (0 <= alpha_x_scan_table_[n])
                {
                    // ### Gets estimate of interval by sign-change-searching
                    //      high and low are the return-values of GiveRootInterval.
                    std::tie(low, high) = this->give_root_interval(n);

                    // ### Finds the root
                    alpha_table[j] = this->get_alpha_root(low, high, n);

                    // Check if we can use method 2 for next roots
                    this->decide_on_method2(n, j);
                }
                // Method 2. ASSUMING ROOTS AT ~FIXED INTERVAL. If next root is expected
                // to lie at distance close enough to estimated distance the root finder
                // can be called using a simple estimated interval.
                else
                {
                    // ### Get interval by simple extrapolation
                    std::tie(low, high) = this->give_root_interval_simple(n, j);

                    // ### Finds the root
                    alpha_table[j] = this->get_alpha_root(low, high, n);
                }
            }
        }
    }
    return alpha_table[i];
}

// =============================================================================
// Root finding algorithm for alpha function.
// Roots (y=0) are necessary to calculate
// Modified by wehrens@amolf.nl. nov, 2011
//
// Note mar 2012:
// The modified root finding algorithm can only work if the roots are calculated
// in a "chronological" sequence (i.e. i = 0, i = 1, i = 2 etc). Therefore roots
// are now all calculated immediately when the roots-table is expanded.
// =============================================================================

GF11_INLINE std::pair<GreensFunction2DRadAbs::real_type, GreensFunction2DRadAbs::real_type>
GreensFunction2DRadAbs::give_root_interval(const std::size_t n /* Order of Bessel functions */) const
{
    // Scans for next interval where a sign change is observed, and thus where a
    // root is expected.

    // # Get/calculate boundaries and determine width of interval to check for
    // sign change.
    const real_type interval( FRACTION_SCAN_INTERVAL() * estimated_alpha_root_distance_);

    // If order (n) is zero, the offset is zero, otherwhise take n-1 offset
    // value as first estimate.
    //      (This can be done because the offsets only get bigger, the roots
    // shift to the right with the order of the Bessel functions n.)
    if (alpha_x_scan_table_[n] == 0) // which implies i == 0
    {
        if (0 < n)
        {
            alpha_x_scan_table_[n] = this->alpha_table_[n-1][0];
        }
    }

    // Define new interval as between x1=offset ("low") and x2=(offset+interval)
    // ("high"). Make sure "low" is > 0.
    real_type low  = alpha_x_scan_table_[n];
    real_type high = alpha_x_scan_table_[n] + interval;
    if (low <= 0) // TODO this check should be redundant
    {
        //low = EPSILON/L_TYPICAL;
        throw_exception<std::runtime_error>("Left alpha search interval boundary < 0.");
    }

    // # Look for the sign change:

    // Get the values of the function at x "low" and x "high".
    // Different for n=0 because this results in a simpler function.
    //      (Note: This code could be optimized by duplicating all involved
    // functions and removing this if-statement.

    // Variables for function values @ resp. left and right boundary interval.
    real_type f_low, f_high;

    if (n == 0)
    {
        f_low  = this->f_alpha0(low);
        f_high = this->f_alpha0(high);
    }
    else
    {
        f_low  = this->f_alpha(low,n);
        f_high = this->f_alpha(high,n);
    }

    // Continue shifting the search interval until a sign change is detected.
    while(0.0 < f_low * f_high)
    {
        low   = high;
        f_low = f_high;

        high  += interval;
        f_high = this->f_alpha(high, n);
    }

    // When above loop has finished, low and high have values inbetween which
    // a root of the alpha function should lie have been found. Make sure that
    // scanning for the next root starts at the end of the domain [low, high]
    // found here.
    alpha_x_scan_table_[n] = high;

    return std::make_pair(low, high);
}

GF11_INLINE std::pair<GreensFunction2DRadAbs::real_type, GreensFunction2DRadAbs::real_type>
GreensFunction2DRadAbs::give_root_interval_simple(
            const std::size_t n /* Order of Bessel functions */,
            const std::size_t i /* i-th root*/) const
{
    // Simply returns an interval based upon previous root, estimated interval
    // inbetween roots and INTERVAL_MARGIN (see .hpp).

    // Offset is simply based on previous root, the interval in which the first
    // root (i=0) lies is never calculated with this function.
    const real_type previous_root = this->get_alpha(n, i-1);

    // Calculates interval [low, high] where root is expected based on the
    // assumption where in a converging regime, where the deviation from this
    // estimate is not more than INTERVAL_MARGIN.
    const real_type low  = previous_root +
                estimated_alpha_root_distance_ * (1 - INTERVAL_MARGIN());
    const real_type high = previous_root +
                estimated_alpha_root_distance_ * (1 + INTERVAL_MARGIN());

    return std::make_pair(low, high);
}

GF11_INLINE GreensFunction2DRadAbs::real_type
GreensFunction2DRadAbs::get_alpha_root(const real_type low, const real_type high,
                                       const std::size_t n /* order of bessel */) const
{
    if (n == 0)
    {
        return get_alpha_root_0(low, high);
    }
    else
    {
        return get_alpha_root_N(low, high, n);
    }
}

GF11_INLINE GreensFunction2DRadAbs::real_type
GreensFunction2DRadAbs::get_alpha_root_0(const real_type low, const real_type high) const
{
    const auto rel_tol = EPSILON();
    const auto abs_tol = rel_tol / L_TYPICAL();

    return find_root<real_type>(
            [this](const real_type alpha) -> real_type {
                return this->f_alpha0(alpha);
            }, low, high, tolerance<real_type>{abs_tol, rel_tol}, 100,
            "GreensFunction2DRadAbs::get_alpha_root_0");
}

GF11_INLINE GreensFunction2DRadAbs::real_type
GreensFunction2DRadAbs::get_alpha_root_N(const real_type low, const real_type high,
                                         const std::size_t n /* order of bessel */) const
{
    const auto rel_tol = EPSILON();
    const auto abs_tol = rel_tol / L_TYPICAL();

    return find_root<real_type>(
            [this, n](const real_type alpha) -> real_type {
                return this->f_alpha(alpha, n);
            }, low, high, tolerance<real_type>{abs_tol, rel_tol}, 100,
            "GreensFunction2DRadAbs::get_alpha_root_N");
}

GF11_INLINE void GreensFunction2DRadAbs::decide_on_method2(
        const std::size_t n, const std::size_t i) const
{
    // The function just counts the number of times the root lies in the interval
    // that can be used to estimate the next root. If this happens enough, then it
    // edits alpha_correctly_estimated_ to a negative value to signify that it can
    // be assumed that for all next roots the distance inbetween roots will
    // equal the interval +/- margin.
    // TODO: This could be more sophisticated.

    // Since the function can only decide with two alpha's already calculated,
    // it can't do anything if i = 0.
    if (i > 0)
    {
        // note the recursiveness!
        const real_type dx = this->get_alpha(n, i) - this->get_alpha(n, i-1);

        // If the relative deviation from the expected difference is smaller
        // than the expected margin, increase number of would-be correct
        // guesses.
        if (std::abs(1.0 - dx / estimated_alpha_root_distance_) < INTERVAL_MARGIN())
        {
            ++alpha_correctly_estimated_[n];
        }
        else
        {
            alpha_correctly_estimated_[n] = 0;
        }

        // If guessing would have worked for last CONVERGENCE_ASSUMED roots,
        // assume it will for all following roots.
        if (alpha_correctly_estimated_[n] > CONVERGENCE_ASSUMED())
        {
            alpha_x_scan_table_[n] = -2; // permanently switch
        }
    }
    return;
}

// -----------------------------------------------------------------------------

GF11_INLINE GreensFunction2DRadAbs::real_type
GreensFunction2DRadAbs::p_survival(const real_type t) const
{
    // Calculates the survival probability at a given time.
    std::vector<real_type> psurv_table;
    return this->p_survival_table(t, psurv_table);
}

GF11_INLINE GreensFunction2DRadAbs::real_type
GreensFunction2DRadAbs::p_survival_table(const real_type t,
                                         std::vector<real_type>& psurv_table) const
{
    // This actually calculates the Survival probability at time t given the
    // particle was at r0 at time 0.
    // It uses the p_surv_table for efficiency (so you don't have to calculate
    // all the constant factors all the time)

    constexpr real_type pi_sq = boost::math::constants::pi_sqr<real_type>();

    const std::size_t maxi = this->guess_maxi(t); // guess the maximum number of iterations required
//  const std::size_t maxi( 500 );             // THIS LEADS TO BIZARRE RESULTS

    // If the estimated # terms needed for convergence is bigger than number
    // of terms summed over (MAX_ALPHA_SEQ), give error.
    if(maxi == MAX_ALPHA_SEQ())
    {
        std::cerr << boost::format("GreensFunction2DRadAbs::p_survival_table (used by drawTime) "
                "couldn't converge; max terms reached: %1%\n") % maxi;
    }

    if(psurv_table.size() < maxi + 1)
    {
        this->get_alpha(0, maxi);              // updates the table of roots
        this->create_psurv_table(psurv_table); // fill data
    }

    // TODO: A sum over terms is performed, where convergence is assumed.
    //       It is not clear if this is a just assumption.
    const real_type p = sumup_all(
            [this, t, &psurv_table](const std::uint32_t i) -> real_type {
                return this->p_survival_i_exp_table(i, t, psurv_table);
            }, maxi);

    return p * pi_sq / 2;
}

GF11_INLINE std::size_t
GreensFunction2DRadAbs::guess_maxi(const real_type t) const
{
    // This tries to guess the maximum number of n iterations it needs for
    // calculating the survival probability

    constexpr real_type pi = boost::math::constants::pi<real_type>();
    constexpr std::uint32_t safety = 2;

    if(std::numeric_limits<real_type>::infinity() <= t)
    {
        return safety;
    }

    const real_type alpha0 = this->get_alpha(0, 0);
    const real_type Dt     = D_ * t;
    const real_type thr    = (std::exp(-Dt * alpha0 * alpha0 ) / alpha0) * EPSILON() * 1e-1;
    const real_type thrsq  = thr * thr;

    if(thrsq <= 0.0)
    {
        return MAX_ALPHA_SEQ();
    }

    const real_type max_alpha = 1.0 / std::sqrt(std::exp(boost::math::lambert_w0(2 * Dt / thrsq)) * thrsq);
    const std::size_t    maxi = safety + static_cast<std::size_t>(max_alpha * (a_ - sigma_) / pi);

    return (std::min<std::size_t>)(maxi, MAX_ALPHA_SEQ());
}

GF11_INLINE void
GreensFunction2DRadAbs::create_psurv_table(std::vector<real_type>& table) const
{
    // get the roots for the survival probability
    const auto& alpha_table_0 = this->alpha_table_[0];

    table.clear();
    table.reserve(alpha_table_0.size());

    for(const auto& alpha : alpha_table_0)
    {
        table.push_back(this->p_survival_i(alpha));
    }
    return ;
}

GF11_INLINE GreensFunction2DRadAbs::real_type
GreensFunction2DRadAbs::p_survival_i(const real_type alpha) const
{
    namespace policies = boost::math::policies;
    constexpr auto policy = policies::make_policy(policies::promote_double<false>());

    // calculates the constant part of the i-th term for the survival probability
    constexpr real_type pi = boost::math::constants::pi<real_type>();

    const real_type s_An = sigma_ * alpha;
    const real_type a_An = a_     * alpha;
    const real_type alpha_sq = alpha * alpha;

    const real_type J0_aAn = boost::math::cyl_bessel_j(0, a_An, policy);
    const real_type J1_sAn = boost::math::cyl_bessel_j(1, s_An, policy);
    const real_type Y0_aAn = boost::math::cyl_neumann (0, a_An, policy);
    const real_type Y1_sAn = boost::math::cyl_neumann (1, s_An, policy);

    // calculate C0,n
    const real_type C_i_0 = this->calc_A_i_0(alpha);

    // calculate the integral over Bn,0
    const real_type dB_n_0dr = J1_sAn * Y0_aAn - Y1_sAn * J0_aAn;

    // this is only the part without alpha of dB0,n(sigma)/dr
    const real_type B_n_0_int = real_type(2.0) / (pi * alpha_sq) - (sigma_ / alpha) * dB_n_0dr;

    return C_i_0 * B_n_0_int;
}

GF11_INLINE GreensFunction2DRadAbs::real_type
GreensFunction2DRadAbs::p_survival_i_exp_table(
        const std::uint32_t i, const real_type t,
        std::vector<real_type>& psurv_table) const
{
    // calculates the ith term with exponent and time for the survival probability
    const real_type alpha = this->get_alpha(0, i);
    return std::exp(-D_ * t * alpha * alpha) * psurv_table[i];
}

// ----------------------------------------------------------------------------

GF11_INLINE void
GreensFunction2DRadAbs::create_Y0J0_tables(std::vector<real_type>& Y0_table,
                                           std::vector<real_type>& J0_table,
                                           std::vector<real_type>& Y0J1J0Y1_table,
                                           const real_type t) const
{
    // Creates the tables with various Bessel functions used in drawR,
    // the table is used to speed things up

    const auto& alpha_table_0 =  this->alpha_table_[0];

    Y0_table      .clear();
    J0_table      .clear();
    Y0J1J0Y1_table.clear();

    Y0_table      .reserve(alpha_table_0.size());
    J0_table      .reserve(alpha_table_0.size());
    Y0J1J0Y1_table.reserve(alpha_table_0.size());

    real_type Y0, J0, Y0J1_J0Y1;
    for(std::size_t i=0; i < alpha_table_0.size(); i++)
    {
        std::tie(Y0, J0, Y0J1_J0Y1) = this->Y0J0J1_constants(alpha_table_0[i], t);

        Y0_table      .push_back(Y0);
        J0_table      .push_back(J0);
        Y0J1J0Y1_table.push_back(Y0J1_J0Y1);
    }
    return;
}

GF11_INLINE std::tuple<GreensFunction2DRadAbs::real_type,
                  GreensFunction2DRadAbs::real_type,
                  GreensFunction2DRadAbs::real_type>
GreensFunction2DRadAbs::Y0J0J1_constants(const real_type alpha, const real_type t) const
{
    namespace policies = boost::math::policies;
    constexpr auto policy = policies::make_policy(policies::promote_double<false>());

    const real_type s_An     = sigma_ * alpha;
    const real_type a_An     = a_     * alpha;
    const real_type alpha_sq = alpha  * alpha;

    const real_type J0_aAn = boost::math::cyl_bessel_j(0, a_An, policy);
    const real_type Y0_aAn = boost::math::cyl_neumann (0, a_An, policy);
    const real_type J1_sAn = boost::math::cyl_bessel_j(1, s_An, policy);
    const real_type Y1_sAn = boost::math::cyl_neumann (1, s_An, policy);

    // calculate An,0
    const real_type A_i_0 = this->calc_A_i_0(alpha);
    // _sq * rho_sq * B_n_0)/( rho_sq - J0_bAn*J0_bAn*(h*h + alpha_sq)));

    // calculate the exponent with the time and product
    const real_type expT     = std::exp(-D_ * alpha_sq * t);
    const real_type Ai0_expT = A_i_0 * expT / alpha;

    // calculate the large constant term in the intergral of Bn,0
    const real_type Y0J1_J0Y1 = Y0_aAn * sigma_ * J1_sAn - J0_aAn * sigma_ * Y1_sAn;

    return std::make_tuple(Ai0_expT * Y0_aAn, Ai0_expT * J0_aAn, Ai0_expT * Y0J1_J0Y1);
}

GF11_INLINE GreensFunction2DRadAbs::real_type
GreensFunction2DRadAbs::p_int_r_table(const real_type r,
                                      const std::vector<real_type>& Y0_aAn_table,
                                      const std::vector<real_type>& J0_aAn_table,
                                      const std::vector<real_type>& Y0J1J0Y1_table) const
{
    // calculates the sum of the sequence for drawR based upon the values in the tables and r

    constexpr real_type pi_sq = boost::math::constants::pi_sqr<real_type>();

    const real_type p = sumup_until_convergence(
        [this, r, &Y0_aAn_table, &J0_aAn_table, &Y0J1J0Y1_table]
        (const std::uint32_t i) -> real_type {
            return p_int_r_i_exp_table(i, r, Y0_aAn_table, J0_aAn_table, Y0J1J0Y1_table);
        }, Y0_aAn_table.size(), 1e-7);

    return p * pi_sq / 2;
}

GF11_INLINE GreensFunction2DRadAbs::real_type
GreensFunction2DRadAbs::p_int_r_i_exp_table(const std::uint32_t i, const real_type r,
                                            const std::vector<real_type>& Y0_aAn_table,
                                            const std::vector<real_type>& J0_aAn_table,
                                            const std::vector<real_type>& Y0J1J0Y1_table) const
{
    namespace policies = boost::math::policies;
    constexpr auto policy = policies::make_policy(policies::promote_double<false>());

    const real_type alpha = this->get_alpha(0, i);
    const real_type r_An  = r * alpha;

    const real_type J1_rAn = boost::math::cyl_bessel_j(1, r_An, policy);
    const real_type Y1_rAn = boost::math::cyl_neumann (1, r_An, policy);

    return Y0_aAn_table[i] * r * J1_rAn - J0_aAn_table[i] * r * Y1_rAn - Y0J1J0Y1_table[i];
}

// ----------------------------------------------------------------------------
// drawTheta related stuff

GF11_INLINE GreensFunction2DRadAbs::real_type
GreensFunction2DRadAbs::ip_theta_n(const std::uint32_t m,
                                   const real_type theta,
                                   const std::vector<real_type>& p_n_table) const
{
    // This calculates the m-th term of the summation for the drawTheta calculation
    // Note that m here starts at 0 and in the equations the sum starts at 1!

    // artificial increase of m to make sure m starts at 1
    const auto m_p1 = m + 1;
    return p_n_table[m_p1] * std::sin(m_p1 * theta) / m_p1;
}

GF11_INLINE GreensFunction2DRadAbs::real_type
GreensFunction2DRadAbs::ip_theta_table(const real_type theta,
                                       const std::vector<real_type>& p_n_table) const
{
    // calculates the cummulative probability of finding the particle at a certain theta
    // It is used by the drawTheta method
    // It uses the p_n_table for it to speed things up

    // get the length of the sum
    // it is shifted one because the first entry should
    // be used (m=0)
    const std::size_t maxm = p_n_table.size() - 1;

    return sumup_all([this, theta, &p_n_table](const std::size_t i) -> real_type {
            return this->ip_theta_n(i, theta, p_n_table);
        }, maxm);
}

GF11_INLINE GreensFunction2DRadAbs::real_type
GreensFunction2DRadAbs::p_m_alpha(const std::uint32_t n, // n-th root
                                  const std::uint32_t m, // bessel order
                                  const real_type r,
                                  const real_type t) const
{
    namespace policies = boost::math::policies;
    constexpr auto policy = policies::make_policy(policies::promote_double<false>());

    // The calculates constant factor m,n for the drawing of theta. These factors are summed later.

    // Gets the n-th root using the bessel functions of order m.
    const real_type alpha = this->get_alpha(m, n);

    const real_type alpha_sq = alpha * alpha;
    const real_type msq = real_type(m) * real_type(m);
    const real_type ssq = sigma_ * sigma_;

    const real_type s_Anm = sigma_ * alpha;
    const real_type a_Anm = a_     * alpha;
    const real_type r0Anm = r0_    * alpha;
    const real_type r_Anm = r      * alpha;

    // calculate the needed bessel functions
    const real_type Jm_sAnm   = boost::math::cyl_bessel_j(m,   s_Anm, policy);
    const real_type Jmp1_sAnm = boost::math::cyl_bessel_j(m+1, s_Anm, policy);    // prime

    const real_type Jm_aAnm   = boost::math::cyl_bessel_j(m, a_Anm, policy);
    const real_type Ym_aAnm   = boost::math::cyl_neumann (m, a_Anm, policy);

    const real_type Jm_r0Anm  = boost::math::cyl_bessel_j(m, r0Anm, policy);
    const real_type Ym_r0Anm  = boost::math::cyl_neumann (m, r0Anm, policy);

    const real_type Jm_rAnm   = boost::math::cyl_bessel_j(m, r_Anm, policy);
    const real_type Ym_rAnm   = boost::math::cyl_neumann (m, r_Anm, policy);

    // calculating An,m
    const real_type h_ma   = h_ - m / sigma_;
    const real_type rho    = h_ma * Jm_sAnm + alpha * Jmp1_sAnm;
    const real_type rho_sq = rho * rho;

    // calculating Bn,m(r')
    const real_type B_n_m_r0 = Jm_r0Anm * Ym_aAnm - Ym_r0Anm * Jm_aAnm;

    const real_type A_n_m = (alpha_sq * rho_sq * B_n_m_r0) /
        (rho_sq - (Jm_aAnm * Jm_aAnm) * (h_ * h_ + alpha_sq - msq / ssq));

    // calculating Bn,m(r*)
    const real_type B_n_m_r = Jm_rAnm * Ym_aAnm - Ym_rAnm * Jm_aAnm;

    return A_n_m * B_n_m_r * std::exp(-D_ * alpha_sq * t);
}

GF11_INLINE GreensFunction2DRadAbs::real_type
GreensFunction2DRadAbs::p_m(const std::uint32_t m,
                            const real_type r,
                            const real_type t) const
{
    // This calculates the m-th constant factor for the drawTheta method.
    // The m-th factor is a summation over n

    return sumup_until_convergence(
            [this, m, r, t](const std::uint32_t i) -> real_type {
                return this->p_m_alpha(i, m, r, t);
            }, MAX_ALPHA_SEQ(), EPSILON());
}

GF11_INLINE void
GreensFunction2DRadAbs::make_p_m_table(std::vector<real_type>& p_m_table,
                                       const real_type r, const real_type t) const
{
    // this should make the table of constants used in the iteration for finding the root for drawTheta
    // The index of the array is consistent with the index of the summation
    p_m_table.clear();

    // This is the p_m where m is 0, for the denominator
    const real_type p_0 = this->p_m(0, r, t);
    const real_type p_1 = this->p_m(1, r, t) / p_0;
    p_m_table.push_back(p_0);
    p_m_table.push_back(p_1);

    if(p_1 == 0)
    {
        return; // all the terms are zero? We are finished
    }

    // get a measure for the allowed error, is this correct?
    const real_type threshold = std::abs(EPSILON() * p_1);

    real_type p_m_abs = std::abs(p_1);
    real_type p_m_prev_abs;
    std::uint32_t m = 1;

    do
    {
        m++;
        if(MAX_ORDER() <= m) // If the number of terms is too large
        {
            std::cerr << boost::format("p_m didn't converge "
                "(m=%1%, t=%2%, r0=%3%, r=%4%, t_est=%5%, continuing...") %
                m % t % r0_ % r % ((r - r0_) * (r - r0_) / D_);
            std::cerr << *this << std::endl;
            break;
        }

        p_m_prev_abs = p_m_abs;
        const real_type p_m = this->p_m(m, r, t) / p_0; // get the next term

        if( ! std::isfinite(p_m))
        {
            std::cerr << boost::format(
                    "makep_m_table: invalid value (p_m = %1%, m=%2%)") % p_m % m;
            break;
        }

        p_m_table.push_back(p_m);
        p_m_abs = std::abs(p_m);
    }
    while (p_m_abs >= threshold || p_m_prev_abs >= threshold || p_m_abs >= p_m_prev_abs );

    return ;
}

GF11_INLINE GreensFunction2DRadAbs::real_type
GreensFunction2DRadAbs::dp_m_alpha_at_a(const std::uint32_t n, // n-th root
                                        const std::uint32_t m, // order of bessel
                                        const real_type t) const
{
    namespace policies = boost::math::policies;
    constexpr auto policy = policies::make_policy(policies::promote_double<false>());
    // This method calculates the constants for the drawTheta method when the particle is at the boundary

    // get the n-th root using the besselfunctions of order m
    const real_type alpha = this->get_alpha(m, n);

    const real_type alpha_sq = alpha * alpha;
    const real_type msq = real_type(m) * real_type(m);
    const real_type ssq = sigma_ * sigma_;

    const real_type s_Anm = sigma_ * alpha;
    const real_type a_Anm = a_     * alpha;
    const real_type r0Anm = r0_    * alpha;

    const real_type Jm_sAnm   = boost::math::cyl_bessel_j(m,   s_Anm, policy);
    const real_type Jmp1_sAnm = boost::math::cyl_bessel_j(m+1, s_Anm, policy);
    const real_type Jm_aAnm   = boost::math::cyl_bessel_j(m,   a_Anm, policy);
    const real_type Ym_aAnm   = boost::math::cyl_neumann (m,   a_Anm, policy);
    const real_type Jm_r0Anm  = boost::math::cyl_bessel_j(m,   r0Anm, policy);
    const real_type Ym_r0Anm  = boost::math::cyl_neumann (m,   r0Anm, policy);

    // calculating An,m
    const real_type h_ma   = h_ - m / sigma_;
    const real_type rho    = h_ma * Jm_sAnm + alpha * Jmp1_sAnm;
    const real_type rho_sq = rho * rho;

    // calculating Bn,m(r')
    const real_type B_n_m_r0 = Jm_r0Anm * Ym_aAnm - Ym_r0Anm * Jm_aAnm;

    const real_type A_n_m = (alpha_sq * rho_sq * B_n_m_r0) /
        (rho_sq - (Jm_aAnm * Jm_aAnm) * (h_ * h_ + alpha_sq - msq / ssq));

    return A_n_m * std::exp(-D_ * alpha_sq * t);
}

GF11_INLINE GreensFunction2DRadAbs::real_type
GreensFunction2DRadAbs::dp_m_at_a(const std::uint32_t m, const real_type t) const
{
    // Makes the sum over n for order m for the constants for the drawtheta Method

    return sumup_until_convergence(
            [this, m, t](const std::uint32_t i) -> real_type {
                return this->dp_m_alpha_at_a(i, m, t);
            },
            MAX_ALPHA_SEQ(), EPSILON());
}

GF11_INLINE void
GreensFunction2DRadAbs::make_dp_m_at_a_table(std::vector<real_type>& p_m_table,
                                             const real_type  t) const
{
    p_m_table.clear();

    // This is the p_m where m is 0, for the denominator
    const real_type p_0 = this->dp_m_at_a(0, t);
    const real_type p_1 = this->dp_m_at_a(1, t) / p_0;
    p_m_table.push_back(p_0);
    p_m_table.push_back(p_1);

    if(p_1 == 0.0)
    {
        return; // all the terms are zero? We are finished
    }

    // get a measure for the allowed error
    const real_type threshold = std::abs(EPSILON() * p_1);

    real_type p_m_abs = std::abs(p_1);
    real_type p_m_prev_abs;
    std::uint32_t m = 1;

    do
    {
        m++;
        if (MAX_ORDER() <= m) // If the number of terms is too large
        {
            std::cerr << boost::format("dp_m didn't converge (m=%1%), continuing...") % m;
            break;
        }

        p_m_prev_abs = p_m_abs;
        const real_type p_m = this->dp_m_at_a(m, t) / p_0; // get the next term

        // DEBUG (something to check in the future?)
        if (p_m_abs == 0.0)
        {
           std::cerr << "Zero valued term found, but convergence is:" <<
               p_m_table[p_m_table.size()-1-1] / p_m_table[p_m_table.size()-2-1];
        }
        // END DEBUG

        if( ! std::isfinite(p_m))
        {
            std::cerr << "make_dp_m_at_a_table: invalid value "
                      << boost::format("(p_m=%1%, m=%2%, t=%3%, p_0=%4%)") %
                      p_m % m % t % p_0;
            break;
        }

        p_m_table.push_back(p_m);
        p_m_abs = std::abs(p_m);
    }
    while (p_m_abs >= threshold || p_m_prev_abs >= threshold || p_m_abs >= p_m_prev_abs);
    // truncate when converged enough.
    // if the current term is smaller than threshold
    // AND the previous term is also smaller than threshold
    // AND the current term is smaller than the previous
}


// ----------------------------------------------------------------------------

GF11_INLINE GreensFunction2DRadAbs::real_type
GreensFunction2DRadAbs::drawTime(const real_type rnd) const
{
    // Draws a first passage time, this could be an escape (through the outer
    // boundary) or a reaction (through the inner boundary)
    constexpr real_type pi = boost::math::constants::pi<real_type>();

    if(!(0.0 <= rnd && rnd < 1.0))
    {
        throw_exception<std::invalid_argument>("GreensFunction2DRadAbs::drawTime: "
            "rnd(", rnd, ") must be in the range [0,1)");
    }
    if (r0_ == a_ || a_ == sigma_)
    {
        return 0.0;
    }

    real_type t_guess;

    const auto pow2 = [](const real_type x) noexcept -> real_type {return x*x;};

    // get some initial guess for the time, dr=sqrt(2dDt) with d
    // the dimensionality (2 in this case)
    const real_type t_Abs = pow2(a_ - r0_) / ( 4.0 * D_);
    if ( kf_ == 0.0 ) // if there was only one absorbing boundary
    {
        t_guess = t_Abs;
    }
    else
    {
        const real_type t_Rad = D_ / pow2(kf_ / (2 * pi * a_)) +
                                     pow2(r0_ - sigma_) / D_;
        // take the shortest time to a boundary
        t_guess = (std::min)( t_Abs, t_Rad );
    }

    t_guess *= 0.1;

    // something with determining the lowest possible t
    const real_type minT = (std::min)(sigma_ * sigma_ / D_ * MIN_T_FACTOR(), t_guess * 1e-7);

    std::vector<real_type> p_surv_table;
    const auto p_surv_eq =
        [this, rnd, &p_surv_table](const real_type t) -> real_type {
            return rnd - this->p_survival_table(t, p_surv_table);
        };

    // put in a upper and lower limit (the picked time cannot be infinite!)
    real_type low  = t_guess;
    real_type high = t_guess;

    const real_type rel_tol = EPSILON();
    const real_type abs_tol = rel_tol * T_TYPICAL();

    // adjust high and low to make sure that f( low ) and f( high ) straddle.
    real_type value = p_surv_eq(t_guess);

    if (value < 0.0)  // if the function is below zero at the guess the upper
    {                 // boundary should be moved (passed the zero point)
        do
        {
            high *= 10;
            value = p_surv_eq(high);

            if(1e10 <= std::abs(high))
            {
                throw_exception<std::runtime_error>("GreensFunction2DRadAbs::drawTime: "
                    "Couldn't adjust high. F(", high, ") = ", value, "; ", (*this));
            }
        }
        while (value < 0.0);
    }
    else // if the function is over zero (or at zero!) then the lower
    {    // boundary should be moved
        real_type value_prev = value;
        do
        {
            low *= 0.1;
            value = p_surv_eq(low);

            if(std::abs(low) <= minT || std::abs(value - value_prev) < abs_tol)
            {
                std::cerr << boost::format("Couldn't adjust low. F(%1%) = %2%") % low % value;
                return low;
            }
            value_prev = value;
        }
        while ( value >= 0.0 );
    }

    return find_root(p_surv_eq, low, high, tolerance<real_type>{abs_tol, rel_tol}, 100,
                    "GreensFunction2DRadAbs::drawTime");
}

GF11_INLINE GreensFunction2DRadAbs::event_type
GreensFunction2DRadAbs::drawEventType(const real_type rnd, const real_type t) const
{
    if(!(0.0 <= rnd && rnd < 1.0))
    {
        throw_exception<std::invalid_argument>("GreensFunction2DRadAbs::drawEventType: "
            "rnd(", rnd, ") must be in the range [0,1)", rnd);
    }
    if (t < 0)
    {
        throw_exception<std::invalid_argument>("GreensFunction2DRadAbs::drawEventTime: "
            "time must be positive: 0.0 < (t = ", t, ")");
    }

    // if there cannot be any flow through the radiating boundary it is always an escape
    if( kf_ == 0.0 )
    {
        return GreensFunction::IV_ESCAPE;
    }

    // First, check if r0 is close only either to a or sigma relative
    // to Dt.  In such cases, the event type is always ESCAPE or REACTION,
    // respectively.   This avoids numerical instability in calculating
    // leavea() and/or leaves().

    // Here, use a rather large threshold for safety.
    const std::uint32_t H = 6; // 6 times the msd travelled as threshold

    const real_type max_dist = H * std::sqrt( 4.0 * D_ * t );
    const real_type a_dist   = a_  - r0_;
    const real_type s_dist   = r0_ - sigma_;

    if (max_dist < a_dist)
    {
        if (s_dist < max_dist)
        {
            return GreensFunction::IV_REACTION;
        }
    }
    else // a_dist < max_dist
    {
        if (max_dist < s_dist)
        {
            return GreensFunction::IV_ESCAPE;
        }
    }

    const real_type reaction = this->leaves(t); // flux through rad boundary
    const real_type escape   = this->leavea(t); // flux through abs boundary
    const real_type value    = reaction / ( reaction + escape );

    if (rnd <= value)
    {
        return GreensFunction::IV_REACTION;
    }
    else
    {
        return GreensFunction::IV_ESCAPE;
    }
}

GF11_INLINE GreensFunction2DRadAbs::real_type
GreensFunction2DRadAbs::drawR(const real_type rnd, const real_type t) const
{
    // This draws a radius R at a given time, provided that the particle was at r0 at t=0

    if(!(0.0 <= rnd && rnd < 1.0))
    {
        throw_exception<std::invalid_argument>("GreensFunction2DRadAbs::drawR: "
            "rnd(", rnd, ") must be in the range [0,1)", rnd);
    }
    if (t < 0)
    {
        throw_exception<std::invalid_argument>("GreensFunction2DRadAbs::drawR: "
            "time must be positive: 0.0 < (t = ", t, ")");
    }

    if( t == 0.0 ) // if no time has passed
    {
        return r0_;
    }

    std::vector<real_type> Y0_aAn_table;
    std::vector<real_type> J0_aAn_table;
    std::vector<real_type> Y0J1J0Y1_table;
    this->create_Y0J0_tables(Y0_aAn_table, J0_aAn_table, Y0J1J0Y1_table, t);

    // adjust low and high starting from r0.
    // this is necessary to avoid root finding in the long tails where
    // numerics can be unstable.
    real_type low   = r0_;
    real_type high  = r0_;
    real_type value = 0;
    std::uint32_t H = 3;

    const real_type rel_tol = EPSILON();
    const real_type abs_tol = rel_tol * L_TYPICAL();

    const real_type psurv       = this->p_survival(t);
    const real_type p_threshold = rnd * psurv;
    const auto p_int_r_eq = [this, t, &Y0_aAn_table, &J0_aAn_table, Y0J1J0Y1_table, p_threshold]
        (const real_type r) -> real_type {
            return this->p_int_r_table(r, Y0_aAn_table, J0_aAn_table, Y0J1J0Y1_table) - p_threshold;
        };

    const real_type msd = std::sqrt(4.0 * D_ * t);

    if(p_int_r_eq(r0_) < 0.0)
    {
        do
        {
            high = r0_ + H * msd;
            if(a_ < high)
            {
                if(p_int_r_eq(a_) < 0.0)
                {
                    // something is very wrong, this should never happen
                    std::cerr << "drawR: p_int_r_table(a) < 0.0. Returning a.\n";
                    return a_;
                }
                high = a_;
                break;
            }
            value = p_int_r_eq(high);
            ++H;
        }
        while (value < 0.0);
    }
    else
    {
        do
        {
            low = r0_ - H * msd;
            if(low < sigma_)
            {
                if(0.0 < p_int_r_eq(sigma_))
                {
                    // something is very wrong, this should never happen
                    std::cerr << "drawR: p_int_r_table(sigma) > 0.0. "
                                 "Returning sigma.\n";
                    return sigma_;
                }
                low = sigma_;
                break;
            }
            value = p_int_r_eq(low);
            ++H;
        }
        while (0.0 < value);
    }

    return find_root(p_int_r_eq, low, high, tolerance<real_type>{abs_tol, rel_tol}, 100,
                     "GreensFunction2DRadAbs::drawR");
}

GF11_INLINE GreensFunction2DRadAbs::real_type
GreensFunction2DRadAbs::drawTheta(const real_type rnd,
                                  const real_type r,
                                  const real_type t) const
{
    // This method draws a theta given a certain r and time (and intial condition of course)

    if(!(0.0 <= rnd && rnd < 1.0))
    {
        throw_exception<std::invalid_argument>("GreensFunction2DAbsSym::drawR: "
            "rnd(", rnd, ") must be in the range [0,1)");
    }
    if (t <= 0)
    {
        throw_exception<std::invalid_argument>("GreensFunction2DAbsSym::drawR: "
            "time must be positive: 0.0 < (t = ", t, ")");
    }
    if(!(sigma_ <= r && r < a_))
    {
        throw_exception<std::invalid_argument>("GreensFunction2DAbsSym::drawR: "
            "r(", r, ") must be in the range [sigma(", sigma_, "), a(", a_, "))");
    }

    // t == 0 means no move.
    if (t <= T_TYPICAL() * EPSILON()                  ||
        D_ == 0                                       ||
        std::abs(r0_ - a_) <= EPSILON() * L_TYPICAL() ||
        rnd <= EPSILON())
    {
        return 0.0;
    }
    else if (r == sigma_) // a reaction has occured, the angle is irrelevant
    {
        return 0.0;
    }

    // making the tables with constants

    std::vector<real_type> p_m_table;
    if(std::abs(r - a_) <= EPSILON() * L_TYPICAL()) // If the r is at the outer boundary
    {
        // making the table if particle on the outer boundary
        this->make_dp_m_at_a_table(p_m_table, t);
    }
    else
    {
        // making the table of constants for the regular case
        this->make_p_m_table(p_m_table, r, t);
    }

    const auto ip_theta_eq = [this, r, t, &p_m_table, rnd]
        (const real_type theta) -> real_type {
            constexpr auto one_div_two_pi = boost::math::constants::one_div_two_pi<real_type>();
            constexpr auto one_div_pi     = 2 * one_div_two_pi;
            return theta * one_div_two_pi +
                   this->ip_theta_table(theta, p_m_table) * one_div_pi -
                   rnd * 0.5;
        };

    constexpr auto pi = boost::math::constants::pi<real_type>();

    return find_root<real_type>(ip_theta_eq, 0.0, pi,
            tolerance<real_type>{EPSILON(), EPSILON()}, 100,
            "GreensFunction2DRadAbs::drawTheta");
}

} // gf11

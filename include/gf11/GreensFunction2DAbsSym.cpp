#ifndef GF11_HEADER_ONLY
#include "GreensFunction2DAbsSym.hpp"
#endif

#include "config.hpp"
#include "find_root.hpp"
#include "throw_exception.hpp"
#include "tolerance.hpp"

#include <boost/math/constants/constants.hpp>
#include <boost/math/special_functions/bessel.hpp>
#include <boost/math/special_functions/erf.hpp>

#include <gsl/gsl_sf_bessel.h>

#include <limits>
#include <stdexcept>
#include <vector>
#include <cstdint>
#include <cassert>
#include <cmath>

namespace gf11
{

GF11_INLINE GreensFunction2DAbsSym::real_type
GreensFunction2DAbsSym::p_survival(const real_type t) const
{
    namespace policies = boost::math::policies;

    constexpr int N_max = 100;

    const real_type Dt = this->D_ * t;
    const real_type ra = real_type(1) / this->a_;

    real_type sum = 0;
    for(int n=1; n<=N_max; ++n)
    {
        // XXX Here, it uses both boost and GSL for the sake of speed.

        // Normal boost is slower than GSL, because of double promotion.
//         const real_type aAn    = boost::math::cyl_bessel_j_zero<real_type>(0, n);
//         const real_type J1_aAn = boost::math::cyl_bessel_j(1, aAn);

        // GSL is ~1.25x faster than the normal boost
//         const real_type aAn    = gsl_sf_bessel_zero_J0(n);
//         const real_type J1_aAn = gsl_sf_bessel_J1(aAn);

        // The following combination is the fastest. This is ~3.0x faster than
        // the normal boost. `promote_double<false>` diables calculation using
        // `long double`. Yes, by default, boost calculates the special functions
        // using `long double` internally.
        // By defining -DBOOST_MATH_PROMOTE_DOUBLE_POLICY=false, boost only uses
        // `double` to calculate the stuff. Since the relative difference between
        // Boost implementation and the original gsl implementation in
        // GreensFunctions is normally 1e-16~1e-15, I think the differences are
        // ignorable...

        const real_type aAn    = gsl_sf_bessel_zero_J0(n);
        const real_type J1_aAn = boost::math::cyl_bessel_j(1, aAn,
            policies::make_policy(policies::promote_double<false>()));

        const real_type An     = aAn * ra;
        const real_type term   = std::exp(-Dt * An * An) / (An * J1_aAn);
        sum += term;
        if(std::abs(term / sum) < CUTOFF())
        {
            break;
        }
    }
    return 2 * ra * sum;
}

GF11_INLINE GreensFunction2DAbsSym::real_type
GreensFunction2DAbsSym::p_int_r(const real_type r, const real_type t) const
{
    namespace policies = boost::math::policies;
    if(r == real_type(0.0))
    {
        // if r == 0, then all the terms are zero because `term` contains `r`
        // in its numerator.
        return 0.0;
    }

    constexpr int N_max = 10000;
    const real_type Dt = this->D_ * t;
    const real_type ra = real_type(1) / this->a_;

    real_type sum = 0.0;
    for(int n=1; n <= N_max; ++n)
    {
        // gsl_sf_bessel_zero_J0 is faster than boost::cyl_bessel_j_zero.
//         const real_type aAn = boost::math::cyl_bessel_j_zero<real_type>(0, n);
        const real_type aAn = gsl_sf_bessel_zero_J0(n);
        const real_type An  = aAn * ra;
        const real_type rAn = r * An;

        // boost::cyl_bessel_j is faster than gsl_sf_bessel_J1.
        const real_type J1_aAn = boost::math::cyl_bessel_j(1, aAn,
            policies::make_policy(policies::promote_double<false>()));
        const real_type J1_rAn = boost::math::cyl_bessel_j(1, rAn,
            policies::make_policy(policies::promote_double<false>()));
        const real_type term = (std::exp(-Dt * An * An) * r * J1_rAn) /
                               (An * J1_aAn * J1_aAn);
        sum += term;
        if(std::abs(term / sum) <= CUTOFF())
        {
            break;
        }
    }
    return 2.0 * ra * ra * sum;
}

GF11_INLINE GreensFunction2DAbsSym::real_type
GreensFunction2DAbsSym::drawTime(const real_type rnd) const
{
    if(!(0.0 <= rnd && rnd < 1.0))
    {
        throw_exception<std::invalid_argument>("GreensFunction2DAbsSym::drawTime: "
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
                throw_exception<std::invalid_argument>("GreensFunction2DAbsSym::"
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
               std::abs(value - low_v) < CUTOFF())
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
                     "GreensFunction2DAbsSym::drawTime");
}

GF11_INLINE GreensFunction2DAbsSym::real_type
GreensFunction2DAbsSym::drawR(const real_type rnd, const real_type t) const
{
    if(!(0.0 <= rnd && rnd < 1.0))
    {
        throw_exception<std::invalid_argument>("GreensFunction2DAbsSym::drawR: "
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


} // gf11

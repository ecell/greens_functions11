#ifndef GF11_HEADER_ONLY
#include "GreensFunction3DRadAbs.hpp"
#endif

#include "config.hpp"
#include "tolerance.hpp"
#include "sumup.hpp"
#include "find_root.hpp"
#include "gf_math.hpp"
#include "SphericalBesselGenerator.hpp"

#include <boost/math/tools/roots.hpp>
#include <boost/math/special_functions/legendre.hpp>

#include <iostream>
#include <limits>
#include <utility>

#include <cmath>

namespace gf11
{

GF11_INLINE GreensFunction3DRadAbs::real_type
GreensFunction3DRadAbs::get_alpha0(const std::uint32_t i) const
{
    // update alpha_table[0]
    auto& alpha0_table = this->alpha_table[0];
    std::size_t current_size = alpha0_table.size();
    if(current_size <= i)
    {
        alpha0_table.resize(i+1);
        for(std::size_t m = current_size; m <= i; ++m)
        {
            alpha0_table[m] = this->alpha0_i(m);
        }
    }
    return alpha0_table[i];
}

GF11_INLINE GreensFunction3DRadAbs::real_type
GreensFunction3DRadAbs::alpha0_i(const std::uint32_t i) const
{
    constexpr real_type pi      = boost::math::constants::pi<real_type>();
    constexpr real_type half_pi = boost::math::constants::half_pi<real_type>();

    const real_type target = i * pi + half_pi;

    // We know the range of the solution from - Pi/2 <= atan <= Pi/2.
    const real_type interval = pi / (a_ - sigma_);

    const real_type low  = i * interval + std::numeric_limits<real_type>::epsilon();
    const real_type high = (i+1) * interval;

    const auto f_alpha0_aux_eq = [target, this](const real_type x) -> real_type {
        return this->f_alpha0_aux(x) - target;
    };

    const tolerance<real_type> tol(/*absolute =*/0.0, /*relative =*/1e-15);
    constexpr std::uintmax_t iter_limit = 100;

    return find_root(f_alpha0_aux_eq, low, high, tol, iter_limit,
                     "GreensFunction3DRadAbs::alpha0_i()");
}

GF11_INLINE GreensFunction3DRadAbs::real_type
GreensFunction3DRadAbs::f_alpha0_aux(const real_type alpha) const
{
    const real_type term1 = (a_ - sigma_) * alpha;
    const real_type term2 = std::atan(this->h_sigma_plus_1 / (sigma_ * alpha));
    return term1 - term2;
}

GF11_INLINE GreensFunction3DRadAbs::real_type
GreensFunction3DRadAbs::num_r0(const real_type alpha) const
{
    const real_type angle_r0 = alpha * (r0_ - sigma_);
    const real_type sin_r0   = std::sin(angle_r0);
    const real_type cos_r0   = std::cos(angle_r0);
    return alpha * sigma_ * cos_r0 + h_sigma_plus_1 * sin_r0;
}

GF11_INLINE GreensFunction3DRadAbs::real_type
GreensFunction3DRadAbs::leave_s(const real_type t) const
{
    return sumup_until_convergence([t, this](const std::uint32_t i) -> real_type
        {
            const real_type alpha = this->get_alpha0(i);
            return std::exp(-alpha * alpha * this->D_ * t) *
                   this->leave_s_i(alpha);
        },
        MAX_ALPHA_SEQ(), TOLERANCE());
}

GF11_INLINE GreensFunction3DRadAbs::real_type
GreensFunction3DRadAbs::leave_a(const real_type t) const
{
    return sumup_until_convergence([t, this](const std::uint32_t i) -> real_type
        {
            const real_type alpha = this->get_alpha0(i);
            return std::exp(-alpha * alpha * this->D_ * t) *
                   this->leave_a_i(alpha);
        },
        MAX_ALPHA_SEQ(), TOLERANCE());
}

GF11_INLINE GreensFunction3DRadAbs::real_type
GreensFunction3DRadAbs::leave_s_i(const real_type alpha) const
{
    const real_type alpha_sq = alpha  * alpha;
    const real_type sigma_sq = sigma_ * sigma_;

    const real_type numer = h_ * alpha * this->num_r0(alpha);
    const real_type denom = 2 * boost::math::constants::pi<real_type>() *
        r0_ * ((a_ - sigma_) * sigma_sq * alpha_sq + h_sigma_plus_1 *
               (a_ + a_ * h_ * sigma_ - h_ * sigma_sq));
    return -D_ * numer / denom;
}

GF11_INLINE GreensFunction3DRadAbs::real_type
GreensFunction3DRadAbs::leave_a_i(const real_type alpha) const
{
    const real_type sigma_sq = sigma_ * sigma_;
    const real_type alpha_sq = alpha  * alpha;
    const real_type angle_a  = alpha * (a_ - sigma_);

    const real_type numer1 = alpha * std::cos(angle_a) * (
            h_sigma_plus_1 * h_sigma_plus_1 + sigma_sq * alpha_sq);

    const real_type numer2 = this->num_r0(alpha);

    const real_type denom = 2 * a_ * boost::math::constants::pi<real_type>() *
            r0_ * h_sigma_plus_1 * (h_sigma_plus_1 *
            (a_ + a_ * h_ * sigma_ - h_ * sigma_sq) +
            (a_ - sigma_) * sigma_sq * alpha_sq);

    return D_ * numer1 * numer2 / denom;
}

GF11_INLINE std::uint32_t
GreensFunction3DRadAbs::guess_max_i(const real_type t) const
{
    constexpr std::uint32_t safety = 2;
    if (t >= std::numeric_limits<real_type>::infinity())
    {
        return safety;
    }

    const real_type alpha0 = this->get_alpha0(0);
    const real_type Dt     = D_ * t;
    const real_type thr    = std::exp(-Dt * alpha0 * alpha0) * TOLERANCE() * 0.1;
    if (thr <= 0.0)
    {
        return MAX_ALPHA_SEQ();
    }

    const real_type max_alpha = std::sqrt(alpha0 * alpha0 - std::log(thr) / Dt);

    const std::uint32_t maxi = safety + static_cast<std::uint32_t>(
        max_alpha * (a_ - sigma_) / boost::math::constants::pi<real_type>());

    return (std::min)(maxi, static_cast<std::uint32_t>(MAX_ALPHA_SEQ()));
}

GF11_INLINE GreensFunction3DRadAbs::real_type
GreensFunction3DRadAbs::p_survival_irr(const real_type t) const
{
    const real_type kD    = 4 * boost::math::constants::pi<real_type>() * sigma_ * D_;
    const real_type alpha = (1.0 + (kf_ / kD)) * (std::sqrt(D_) / sigma_);

    const real_type sqrt_t = std::sqrt(t);
    const real_type sqrt_D = std::sqrt(D_);

    const real_type r0_minus_sigma_over_sqrt4D_t =
        (r0_ - sigma_) / ((sqrt_D + sqrt_D) * sqrt_t);

    const real_type Wf = ::gf11::W(r0_minus_sigma_over_sqrt4D_t, alpha * sqrt_t);
    const real_type factor = sigma_ * kf_ / (r0_ * (kf_ + kD));

    return 1.0 - factor * (std::erfc(r0_minus_sigma_over_sqrt4D_t) - Wf);
}

GF11_INLINE GreensFunction3DRadAbs::real_type
GreensFunction3DRadAbs::p_survival_nocollision(const real_type t) const
{
    constexpr real_type pi    = boost::math::constants::pi    <real_type>();
    constexpr real_type pi_sq = boost::math::constants::pi_sqr<real_type>();
    constexpr real_type one_over_pi  =
        2 * boost::math::constants::one_div_two_pi<real_type>();

    const real_type Dt    = D_  * t;
    const real_type asq   = a_  * a_;
    const real_type a_r   = 1.0 / a_;
    const real_type asq_r = a_r * a_r;
    const real_type PIr0  = pi  * r0_;

    const real_type angle_factor = PIr0 * a_r;
    const real_type exp_factor   = -Dt * pi_sq * asq_r;

    const std::uint32_t i_max = (std::max)(static_cast<std::uint32_t>(std::ceil(
        std::sqrt(pi_sq - asq * std::log(TOLERANCE()) / Dt) * one_over_pi)), 2u);

    real_type p    = 0.0;
    real_type sign = 1.0;
    for(std::size_t i=1; i<=i_max; ++i)
    {
        p += sign * std::exp(exp_factor * i * i) * std::sin(angle_factor * i)/i;
        sign = -sign;
    }
    return p * (a_ + a_) / PIr0;
}

GF11_INLINE GreensFunction3DRadAbs::real_type
GreensFunction3DRadAbs::p_survival_i(const real_type alpha) const
{
    const real_type sigma_sq = sigma_ * sigma_;
    const real_type alpha_sq = alpha  * alpha;

    const real_type angle_a = alpha * (a_ - sigma_);
    const real_type cos_a   = std::cos(angle_a);

    const real_type num1 = h_ * sigma_sq * h_sigma_plus_1 - a_ *
        (h_sigma_plus_1 * h_sigma_plus_1 + sigma_sq * alpha_sq) * cos_a;
    const real_type num2 = this->num_r0(alpha);

    const real_type den = r0_ * h_sigma_plus_1 * alpha *
        (-h_sigma_plus_1 * (a_ + a_ * h_ * sigma_ - h_ * sigma_sq) +
         (sigma_ - a_) * sigma_sq * alpha_sq);

    return -2.0 * num1 * num2 / den;
}

GF11_INLINE void GreensFunction3DRadAbs::create_p_survival_table(
        std::vector<real_type>& table) const
{
    const auto& alpha_table_0 = this->alpha_table[0];

    table.clear();
    table.resize(alpha_table_0.size());

    std::transform(alpha_table_0.begin(), alpha_table_0.end(), table.begin(),
        [this](const real_type alpha) -> real_type {
            return this->p_survival_i(alpha);
        });
    return;
}

GF11_INLINE GreensFunction3DRadAbs::real_type
GreensFunction3DRadAbs::p_survival(
        const real_type t, std::vector<real_type>& table) const
{
    const real_type dist_to_a = a_  - r0_;
    const real_type dist_to_s = r0_ - sigma_;

    constexpr real_type H = 6.0; // a fairly strict criterion for safety.
    const real_type max_dist = H * std::sqrt(6.0 * D_ * t);

    if(dist_to_a > max_dist)
    {
        if(dist_to_s > max_dist) // far from anything; it'll survive.
        {
            return 1.0;
        }
        else // close only to s, ignore a
        {
            return this->p_survival_irr(t);
        }
    }
    else
    {
        if (dist_to_s > max_dist)  // close only to a, ignore s.
        {
            return this->p_survival_nocollision(t);
        }
        else  // close to both boundaries.  do the normal calculation.
        {
            const std::uint32_t maxi = this->guess_max_i(t);
            if(table.size() < maxi + 1)
            {
                this->get_alpha0(maxi); // updates alpha_table[0]
                this->create_p_survival_table(table);
            }
            return sumup_all(
                [&table, t, this](const std::uint32_t i) -> real_type {
                    const real_type alpha = this->get_alpha0(i);
                    return std::exp(-(this->D_) * t * alpha * alpha) * table[i];
                }, maxi);
        }
    }
}

GF11_INLINE GreensFunction3DRadAbs::real_type
GreensFunction3DRadAbs::p_int_r_i(
        const real_type r, const real_type alpha, const real_type num_r0) const
{
    const real_type angle_r = alpha * (r - this->sigma_);
    const real_type sin_r = std::sin(angle_r);
    const real_type cos_r = std::cos(angle_r);
    // we can expect the compiler uses __sincos here.

    const real_type sigma_sq = sigma_ * sigma_;
    const real_type alpha_sq = alpha  * alpha;
    const real_type h_sigma  = h_ * sigma_;

    const real_type numer = alpha *
        (h_sigma * sigma_ - h_sigma * r * cos_r - (r - sigma_) * cos_r) +
        (h_sigma_plus_1 + r * sigma_ * alpha_sq) * sin_r;

    const real_type denom = r0_ * alpha_sq *
        ((a_ - sigma_) * sigma_sq * alpha_sq +
         h_sigma_plus_1 * (a_ + a_ * h_ * sigma_ - h_ * sigma_sq));

    return 2.0 * num_r0 * numer / denom;
}

GF11_INLINE GreensFunction3DRadAbs::real_type
GreensFunction3DRadAbs::p_int_r(
        const real_type r, const real_type t) const
{
    // do not use sumup_all. series acceleration is indispensable.
    return sumup_until_convergence(
            [r, t, this](const std::uint32_t i) -> real_type {
                const real_type alpha = this->get_alpha0(i);
                return std::exp(-(this->D_ * t * alpha * alpha)) *
                       this->p_int_r_i(r, alpha, this->num_r0(alpha));
            }, MAX_ALPHA_SEQ(), TOLERANCE());
}

GF11_INLINE GreensFunction3DRadAbs::real_type
GreensFunction3DRadAbs::ip_theta(const real_type theta,
        const real_type r, const std::vector<real_type>& p_n_table) const
{
    const real_type cos_theta = std::cos(theta);

    //XXX legendre_table is offset by 1 to incorporate the n=-1 case.
    //    legendre_table[0] contains legendre_p(-1, x),
    //    legendre_table[1] contains legendre_p( 0, x) ...
    std::vector<real_type> legendre_table(p_n_table.size() + 2);
    legendre_table[0] = 1.0; // n = -1
    legendre_table[1] = boost::math::legendre_p(0, cos_theta); // n = 0
    legendre_table[2] = boost::math::legendre_p(1, cos_theta); // n = 1
    for(std::size_t i=3; i<legendre_table.size(); ++i)
    {
        // legendre_next(l, Pl(x), Pl-1(x)) calculates Pl+1(x) by using
        // recurrence relationship.
        legendre_table[i] = boost::math::legendre_next(
                i-2, cos_theta, legendre_table[i-1], legendre_table[i-2]);
    }

    return sumup_all([&p_n_table, &legendre_table, this]
        (const std::uint32_t n) -> real_type {
            const real_type Pn_m1 = legendre_table[n];  // n-1
            const real_type Pn_p1 = legendre_table[n+2];// n+1
            return p_n_table[n] * (Pn_m1 - Pn_p1);
        }, p_n_table.size());
}

GF11_INLINE GreensFunction3DRadAbs::real_type
GreensFunction3DRadAbs::f_alpha(
        const real_type alpha, const std::uint32_t n) const
{
    const real_type a_alpha         = this->a_     * alpha;
    const real_type sigma_alpha     = this->sigma_ * alpha;
    const real_type h_sigma         = this->h_     * this->sigma_;
    const real_type real_n          = static_cast<real_type>(n);
    const real_type h_sigma_minus_n = h_sigma - real_n;

    const SphericalBesselGenerator<real_type> s;

    const real_type js1 = s.j(n,   sigma_alpha);
    const real_type ys1 = s.y(n,   sigma_alpha);
    const real_type js2 = s.j(n+1, sigma_alpha);
    const real_type ys2 = s.y(n+1, sigma_alpha);
    const real_type ja  = s.j(n,   a_alpha);
    const real_type ya  = s.y(n,   a_alpha);

    constexpr real_type one_div_pi =
        2.0 * boost::math::constants::one_div_two_pi<real_type>();

    const real_type term1 = (h_sigma_minus_n * js1 + sigma_alpha * js2) * ya;
    const real_type term2 = (h_sigma_minus_n * ys1 + sigma_alpha * ys2) * ja;
    const real_type factor = 2.0 * alpha * std::sqrt(a_ * sigma_) * one_div_pi;

    return (term1 - term2) * factor;
}

GF11_INLINE std::uint32_t
GreensFunction3DRadAbs::alpha_offset(const std::size_t n) const
{
    if(this->alpha_offset_table[n] != 0xFFFFFFFF)
    {
        return this->alpha_offset_table[n];
    }
    assert(this->alpha_offset_table.size() >= n);

    constexpr real_type pi      = boost::math::constants::pi<real_type>();
    constexpr real_type half_pi = boost::math::constants::half_pi<real_type>();

    std::uint32_t   offset = this->alpha_offset_table[n-1];
    const real_type factor = 1.0 / (this->a_ - this->sigma_);
    real_type       target = offset * pi + half_pi;

    // We know the range of the solution from - Pi/2 <= atan <= Pi/2.
    const real_type alpha_mid        = target  * factor;
    const real_type alpha_half_range = half_pi * factor;

    real_type low  = alpha_mid - alpha_half_range * (1.0 - 1e-3); // avoid zero.
    real_type high = alpha_mid + alpha_half_range;

    // Here we find the interval where the first positive root is in.
    // We find the first pair of alpha
    // (Pi * offset + Pi/2) +- Pi/2 / (a - sigma)
    // where the values of f_alpha() straddle.
    // The assumption is the interval between roots is not much
    // smaller than Pi / (a - sigma).

    real_type low_value  = this->f_alpha(low,  n);
    real_type high_value = this->f_alpha(high, n);

    // this can be much faster if better initial guess is given.
    while(low_value * high_value >= 0)
    {
        ++offset;
        target = pi * offset + half_pi;
        low    = (target - half_pi) * factor;
        high   = (target + half_pi) * factor;

        low_value  = high_value;
        high_value = this->f_alpha(high, n);
    }
    this->alpha_offset_table[n] = offset;
    return offset;
}

GF11_INLINE GreensFunction3DRadAbs::real_type
GreensFunction3DRadAbs::f_alpha_aux(
        const real_type alpha, const std::uint32_t n) const
{
    if (alpha == 0.0)
    {
        return -1.0;
    }

    const real_type a_alpha         = a_     * alpha;
    const real_type sigma_alpha     = sigma_ * alpha;
    const real_type n_minus_h_sigma = n - h_ * sigma_;

    /*(a - s) u -
      ArcTan[(P[n, a u] ((-n + h s) P[n, s u] + s u Q[1 + n, s u]) -
               Q[n, a u] (s u P[1 + n, s u] + (n - h s) Q[n, s u]))/
             (Q[n, a u] ((-n + h s) P[n, s u] + s u Q[1 + n, s u]) + 
               P[n, a u] (s u P[1 + n, s u] + (n - h s) Q[n, s u]))]
    */

    const real_type Pa = P(n, a_alpha);
    const real_type Qa = Q(n, a_alpha);

    real_type Ps, Psp;
    std::tie(Ps, Psp) = P2(n, sigma_alpha);

    real_type Qs, Qsp;
    std::tie(Qs, Qsp) = Q2(n, sigma_alpha);

    const real_type n_minus_h_sigma_Ps = n_minus_h_sigma * Ps;
    const real_type n_minus_h_sigma_Qs = n_minus_h_sigma * Qs;
    const real_type sigma_alpha_Psp    = sigma_alpha * Psp;
    const real_type sigma_alpha_Qsp    = sigma_alpha * Qsp;

    const real_type Qa_Pa = Qa / Pa;

    const real_type A = sigma_alpha_Qsp - n_minus_h_sigma_Ps;
    const real_type B = sigma_alpha_Psp + n_minus_h_sigma_Qs;

    // this form, dividing all terms by Pa, prevents overflow.
    const real_type angle = (A - Qa_Pa * B) / (Qa_Pa * A + B);
    const real_type term1 = (a_ - sigma_) * alpha;
    const real_type term2 = std::atan(angle);
    return term1 - term2;
}


GF11_INLINE GreensFunction3DRadAbs::real_type
GreensFunction3DRadAbs::alpha_i(
        const std::uint32_t i, const std::uint32_t n) const
{
    constexpr real_type pi      = boost::math::constants::pi<real_type>();
    constexpr real_type half_pi = boost::math::constants::half_pi<real_type>();

    const real_type target = pi * i + half_pi;
    const real_type factor = 1.0 / (this->a_ - this->sigma_);
    const real_type low    = (target - half_pi) * factor;
    const real_type high   = (target + half_pi) * factor;

    const auto f_alpha_aux_eq =
        [n, target, this](const real_type alpha) -> real_type {
        return this->f_alpha_aux(alpha, n) - target;
    };
    const tolerance<real_type> tol(/*abs = */1e-6, /*rel = */1e-15);
    return find_root(f_alpha_aux_eq, low, high, tol, MAX_ITERATION(),
                     "GreensFunction3DRadAbs::alpha_i()");
}

GF11_INLINE GreensFunction3DRadAbs::real_type
GreensFunction3DRadAbs::get_alpha(
        const std::size_t n, const std::uint32_t i) const
{
    auto& alpha_table = this->alpha_table[n];
    const auto current_size = alpha_table.size();

    if(current_size <= i)
    {
        alpha_table.resize(i+1);
        const std::size_t offset = this->alpha_offset(n);

        for(std::size_t m = current_size; m <= i; ++m)
        {
            alpha_table[m] = this->alpha_i(m + offset, n);
        }
    }
    return alpha_table[i];
}

GF11_INLINE GreensFunction3DRadAbs::real_type
GreensFunction3DRadAbs::p_n_alpha(
        const std::uint32_t i, const std::uint32_t n,
        const real_type r, const real_type t) const
{
    const real_type mDt             = -(this->D_) * t;
    const real_type alpha           = this->get_alpha(n, i);
    const real_type alpha_sq        = alpha * alpha;
    const real_type a_alpha         = this->a_ * alpha;
    const real_type sigma_alpha     = this->sigma_ * alpha;
    const real_type h_sigma         = this->h_ * this->sigma_;
    const real_type real_n          = static_cast<real_type>(n);
    const real_type h_sigma_minus_n = h_sigma - real_n;

    const real_type term1 = alpha_sq * alpha_sq * std::exp(mDt * alpha_sq);

    const SphericalBesselGenerator<real_type> s;

    const real_type js1 = s.j(n,   sigma_alpha);
    const real_type js2 = s.j(n+1, sigma_alpha);
    const real_type ja  = s.j(n,   a_alpha);
    const real_type ya  = s.y(n,   a_alpha);
    const real_type jr  = s.j(n,   r * alpha);
    const real_type yr  = s.y(n,   r * alpha);
    const real_type jr0 = s.j(n,   r0_ * alpha);
    const real_type yr0 = s.y(n,   r0_ * alpha);

    const real_type J   = h_sigma_minus_n * js1 + sigma_alpha * js2;
    const real_type Jsq = J * J;
    const real_type JY1 = ja * yr  - ya * jr;
    const real_type JY2 = ja * yr0 - ya * jr0;

    const real_type numer = Jsq * JY1 * JY2;
    const real_type denom = (a_ * (real_n + real_n * real_n -
        sigma_ * (h_ + h_ * h_ * sigma_ + sigma_ * alpha_sq)) * ja * ja) +
        sigma_ * Jsq;
    return term1 * numer / denom;
}

GF11_INLINE GreensFunction3DRadAbs::real_type
GreensFunction3DRadAbs::p_n(
        const std::size_t n, const real_type r, const real_type t,
        const real_type max_alpha) const
{
    constexpr std::uint32_t min_i = 2;

    real_type p = 0.0;
    for(std::uint32_t i = 0; i<=MAX_ALPHA_SEQ(); ++i)
    {
        const real_type alpha = get_alpha(n, i);
        const real_type p_i   = this->p_n_alpha(i, n, r, t);
        p += p_i;

        if(alpha >= max_alpha && i >= min_i)
        {
            break;
        }
    }
    return p;
}

GF11_INLINE std::vector<GreensFunction3DRadAbs::real_type>
GreensFunction3DRadAbs::make_p_n_table(
        const real_type r, const real_type t) const
{
    constexpr real_type one_div_two_pi =
        boost::math::constants::one_div_two_pi<real_type>();

    const real_type factor  = this->a_ * this->sigma_ * one_div_two_pi;
    const real_type Dt      = this->D_ * t;
    const real_type alpha00 = this->get_alpha(0, 0);

    const real_type max_alpha = std::sqrt(
        alpha00 * alpha00 - std::log(THETA_TOLERANCE() * 0.1) / Dt);

    const real_type p_0 = this->p_n(0, r, t, max_alpha) * factor;

    std::vector<real_type> p_n_table;
    p_n_table.reserve(MAX_ORDER() + 1);

    p_n_table.push_back(p_0);
    if(p_0 == 0)
    {
        return p_n_table;
    }

    const real_type threshold = std::abs(THETA_TOLERANCE() * p_0);
    real_type    p_n_prev_abs = std::abs(p_0);

    // we have the lookup table for $n \in [0, MAX_ORDER]$. n can be MAX_ORDER.
    for(std::uint32_t n=1; n<=MAX_ORDER(); ++n)
    {
        if(this->get_alpha(n, 0) >= max_alpha)
        {
            return p_n_table;
        }

        const real_type p_n_ = this->p_n(n, r, t, max_alpha) * factor;

        p_n_table.push_back(p_n_);
        const real_type p_n_abs = std::abs(p_n_);

        // truncate when converged enough.
        if(p_n_abs <  threshold && p_n_prev_abs < threshold &&
           p_n_abs <= p_n_prev_abs)
        {
            return p_n_table;
        }
        p_n_prev_abs = p_n_abs;
    }
    return p_n_table;
}

GF11_INLINE GreensFunction3DRadAbs::real_type
GreensFunction3DRadAbs::dp_n_alpha_at_a(
    const std::uint32_t i, const std::uint32_t n, const real_type t) const
{
    const real_type Dt              = this->D_ * t;
    const real_type alpha           = this->get_alpha(n, i);
    const real_type alpha_sq        = alpha * alpha;
    const real_type a_alpha         = this->a_     * alpha;
    const real_type sigma_alpha     = this->sigma_ * alpha;
    const real_type r0_alpha        = this->r0_    * alpha;
    const real_type h_sigma         = this->h_     * this->sigma_;
    const real_type real_n          = static_cast<real_type>(n);
    const real_type h_sigma_minus_n = h_sigma - real_n;

    const SphericalBesselGenerator<real_type> s;

    const real_type js1 = s.j(n,   sigma_alpha);
    const real_type js2 = s.j(n+1, sigma_alpha);
    const real_type ja  = s.j(n,   a_alpha);
    const real_type ya  = s.y(n,   a_alpha);
    const real_type jr0 = s.j(n,   r0_alpha);
    const real_type yr0 = s.y(n,   r0_alpha);

    const real_type J   = h_sigma_minus_n * js1 + sigma_alpha * js2;
    const real_type Jsq = J * J;
    const real_type JY  = -jr0 * ya + ja * yr0;

    const real_type numer1 = alpha_sq * alpha * std::exp(-Dt * alpha_sq);
    const real_type numer2 = Jsq * JY;

    const real_type denom1 = a_ * (real_n + real_n * real_n -
            sigma_ * (h_ + h_ * h_ * sigma_ + sigma_ * alpha_sq)) * ja * ja;
    const real_type denom2 = sigma_ * Jsq;

    return numer1 * numer2 / (denom1 + denom2);
}

GF11_INLINE GreensFunction3DRadAbs::real_type
GreensFunction3DRadAbs::dp_n_at_a(
    const std::uint32_t n, const real_type t, const real_type max_alpha) const
{
    constexpr std::uint32_t min_i = 2u;

    real_type p = 0.0;
    for(std::uint32_t i = 0; i <= MAX_ALPHA_SEQ(); ++i)
    {
        p += this->dp_n_alpha_at_a(i, n, t);

        if(this->get_alpha(n, i) >= max_alpha && i >= min_i)
        {
            break;
        }
    }
    return p;
}

GF11_INLINE std::vector<GreensFunction3DRadAbs::real_type>
GreensFunction3DRadAbs::make_dp_n_at_a_table(const real_type t) const
{
    constexpr real_type two_pi = boost::math::constants::two_pi<real_type>();

    const real_type factor    = this->D_ * this->sigma_ / (this->a_ * two_pi);
    const real_type Dt        = this->D_ * t;
    const real_type alpha00   = this->get_alpha(0, 0);
    const real_type max_alpha = std::sqrt(Dt * alpha00 * alpha00 -
                                          std::log(THETA_TOLERANCE() * 0.1) / Dt);
    const real_type p_0       = this->dp_n_at_a(0, t, max_alpha) * factor;

    if(p_0 == 0)
    {
        return {p_0};
    }

    std::vector<real_type> p_n_table;
    p_n_table.reserve(MAX_ORDER());
    p_n_table.push_back(p_0);

    const real_type  threshold = std::abs(THETA_TOLERANCE() * p_0);
    real_type     p_n_prev_abs = std::abs(p_0);
    for(std::uint32_t n=1; n<=MAX_ORDER(); ++n)
    {
        if(this->get_alpha(n, 0) >= max_alpha)
        {
            return p_n_table;
        }

        const real_type p_n     = this->dp_n_at_a(n, t, max_alpha) * factor;
        const real_type p_n_abs = std::abs(p_n);

        p_n_table.push_back(p_n);

        // truncate when converged enough.
        if(p_n_abs      <  threshold &&
           p_n_prev_abs <  threshold &&
           p_n_abs      <= p_n_prev_abs)
        {
            return p_n_table;
        }
        p_n_prev_abs = p_n_abs;
    }
    return p_n_table;
}

GF11_INLINE GreensFunction3DRadAbs::real_type
GreensFunction3DRadAbs::drawTime(const real_type rnd) const
{
    if(!(0.0 <= rnd && rnd < 1.0))
    {
        throw_exception<std::invalid_argument>("GreensFunction3DRadAbs::drawTime: "
            "rnd(", rnd, ") must be in the range [0,1)");
    }
    if (r0_ == a_ || a_ == sigma_)
    {
        return 0.0;
    }

    const real_type dist = (kf_ == 0) ? (a_-r0_) : (std::min)(a_-r0_, r0_-sigma_);
    const real_type t_guess = (dist * dist / (6.0 * D_)) * 0.1;
    const real_type min_t = (std::min)(sigma_ * sigma_ / D_ * MIN_T_FACTOR(),
                                       t_guess * 1e-6);

    std::vector<real_type> table;
    auto p_surv_eq = [&table, rnd, this](const real_type t) -> real_type {
        return rnd - this->p_survival(t, table);
    };
    real_type low  = t_guess;
    real_type high = t_guess;

    const real_type value = p_surv_eq(t_guess);
    real_type low_value  = value;
    real_type high_value = value;
    if(value < 0.0)
    {
        while(true)
        {
            high *= 10.0;
            high_value = p_surv_eq(high);

            if(0.0    <= high_value) {break;}
            if(1.0e10 <= std::abs(high))
            {
                throw_exception<std::runtime_error>("GreensFunction3DRadAbs: "
                    "couldn't adjust high. p_surv_eq(", high, ")=", high_value,
                    ", rnd=", rnd, ", gf=", *this);
            }
        }
    }
    else // 0.0 < p_surv(t_guess)
    {
        real_type low_value_prev = value;
        while(true)
        {
            low *= 0.1;
            low_value = p_surv_eq(low);

            if(low_value <= 0.0) {break;}
            if(std::abs(low) <= min_t ||
               std::abs(low_value - low_value_prev) < TOLERANCE())
            {
                // converged!?
                return low;
            }
            low_value_prev = low_value;
        }
    }

    const tolerance<real_type> tol(/*absolute =*/0.0, /*relative =*/TOLERANCE());
    return find_root(p_surv_eq, low, high, low_value, high_value, tol, MAX_ITERATION(),
                     "GreensFunction3DRadAbs::drawTime()");
}

GF11_INLINE GreensFunction
GreensFunction3DRadAbs::drawEventType(
        const real_type rnd, const real_type t) const
{
    if(!(0.0 <= rnd && rnd < 1.0))
    {
        throw_exception<std::invalid_argument>(
            "GreensFunction3DRadAbs::drawEventType: "
            "rnd(", rnd, ") must be in the range [0,1)");
    }
    if(t <= 0.0)
    {
        throw_exception<std::invalid_argument>(
            "GreensFunction3DRadAbs::drawEventType: t(", t, ") must be positive");
    }
    if (kf_ == 0.0)
    {
        return GreensFunction::IV_ESCAPE;
    }

    // First, check if r0 is close only either to a or sigma relative
    // to Dt.  In such cases, the event type is always IV_ESCAPE or
    // IV_REACTION, respectively. This avoids numerical instability in
    // calculating leavea() and/or leaves().

    // Here, use a rather large threshold for safety.
    constexpr std::uint32_t H = 6u;
    const real_type max_dist = H * std::sqrt(6.0 * D_ * t);
    const real_type a_dist = a_  - r0_;
    const real_type s_dist = r0_ - sigma_;

    if (a_dist > max_dist)
    {
        if (s_dist < max_dist)
        {
            return GreensFunction::IV_REACTION;
        }
    }
    else // a_dist < max_dist
    {
        if (s_dist > max_dist)
        {
            return GreensFunction::IV_ESCAPE;
        }
    }

    constexpr real_type pi = boost::math::constants::pi<real_type>();
    const real_type reaction = this->leave_s(t) * 4.0 * pi * sigma_ * sigma_;
    const real_type escape   = this->leave_a(t) * 4.0 * pi * a_ * a_;
    const real_type value    = reaction / (reaction + escape);

    if (rnd <= value)
    {
        return GreensFunction::IV_REACTION; // leave_s
    }
    else
    {
        return GreensFunction::IV_ESCAPE;   // leave_a
    }
}

GF11_INLINE GreensFunction3DRadAbs::real_type
GreensFunction3DRadAbs::drawR(
        const real_type rnd, const real_type t) const
{
    if(!(0.0 <= rnd && rnd < 1.0))
    {
        throw_exception<std::invalid_argument>("GreensFunction3DRadAbs::drawR: "
            "rnd(", rnd, ") must be in the range [0,1)", rnd);
    }
    if (t == 0.0)
    {
        return r0_;
    }

    std::vector<real_type> table; //TODO check: can be stored as a member val?
    const real_type target = rnd * this->p_survival(t, table);

    const auto p_int_r_eq = [t, target, this](const real_type r) -> real_type {
        return this->p_int_r(r, t) - target;
    };

    // adjust low and high starting from r0.
    // this is necessary to avoid root finding in the long tails where
    // numerics can be unstable.

    real_type low  = this->r0_;
    real_type high = this->r0_;
    const real_type sqrt6Dt = std::sqrt(6.0 * this->D_ * t);

    const real_type value = p_int_r_eq(this->r0_);
    real_type low_value  = value;
    real_type high_value = value;
    if(value < 0.0)
    {
        // low = r0
        std::uint32_t H = 3;
        while(true)
        {
            high = this->r0_ + H * sqrt6Dt;
            if(high > this->a_)
            {
                if(p_int_r_eq(this->a_) < 0.0)
                {
                    return this->a_;
                }
                high = this->a_;
                high_value = p_int_r_eq(high);
                break;
            }

            high_value = p_int_r_eq(high);
            if(high_value > 0.0)
            {
                break;
            }
            ++H;
        }
    }
    else
    {
        // high = r0
        std::uint32_t H = 3;
        while(true)
        {
            low = this->r0_ - H * sqrt6Dt;
            if (low < this->sigma_)
            {
                if(p_int_r_eq(this->sigma_) > 0.0)
                {
                    return this->sigma_;
                }
                low = this->sigma_;
                low_value = p_int_r_eq(low);
                break;
            }

            low_value = p_int_r_eq(low);
            if(low_value < 0.0)
            {
                break;
            }
            ++H;
        }
    }
    const tolerance<real_type> tol(/*absolute =*/1e-15, /*relative =*/TOLERANCE());
    return find_root(p_int_r_eq, low, high, low_value, high_value, tol, MAX_ITERATION(),
                     "GreensFunction3DRadAbs::drawR");
}

GF11_INLINE GreensFunction3DRadAbs::real_type
GreensFunction3DRadAbs::drawTheta(
        const real_type rnd, const real_type r, const real_type t) const
{
    if(!(rnd < 1.0 && rnd >= 0.0))
    {
        throw_exception<std::invalid_argument>("GreensFunction3DRadAbs::drawTheta:"
            " rnd(", rnd, ") must be in the range [0,1)");
    }
    if(this->sigma_ > r)
    {
        throw_exception<std::invalid_argument>("GreensFunction3DRadAbs::drawTheta:"
            " r(", r, ") must be larger than or equal to sigma(", sigma_, ")");
    }
    if(t < 0.0)
    {
        throw_exception<std::invalid_argument>("GreensFunction3DRadAbs::drawTheta:"
            " t(", t, ") must be positive");
    }

    if(t == 0.0)
    {
        return 0.0;
    }

    constexpr real_type pi = boost::math::constants::pi<real_type>();

    const std::vector<real_type> p_n_table = (r >= this->a_) ?
        this->make_dp_n_at_a_table(t) : this->make_p_n_table(r, t);

    const real_type target = rnd * this->ip_theta(pi, r, p_n_table);

    const auto ip_theta_eq =
        [&p_n_table, r, target, this](const real_type theta) -> real_type {
        return this->ip_theta(theta, r, p_n_table) - target;
    };

    const tolerance<real_type> tol(
            /*absolute =*/1e-11, /*relative =*/THETA_TOLERANCE());

    return find_root(ip_theta_eq, 0.0, pi, tol, MAX_ITERATION(),
                     "GreensFunction3DRadAbs::drawTheta()");
}

} // gf11

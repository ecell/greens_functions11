#ifndef GF11_MATHEMATICAL_UTILITY_H
#define GF11_MATHEMATICAL_UTILITY_H
#include <boost/math/special_functions/erf.hpp>
#include <boost/math/constants/constants.hpp>
#include <cmath>

namespace gf11
{

/**
   Calculates std::exp(x^2) * erfc(x)

   See asymptotic expansion here:
   http://en.wikipedia.org/wiki/Error_function
*/
template<typename realT>
inline realT expxsq_erfc(const realT x) noexcept
{
    const auto xsq = x * x;
    if(x > realT(26.0))
    {
        constexpr auto one_div_root_pi =
            boost::math::constants::one_div_root_pi<realT>();
        const auto x2sq_r = realT(1) / (realT(2) * xsq); // 2 / (2 x)^2
        /*
          up to second term in the expansion.
          abs err ~= 9e-8 at x == 20, 3e-8 at x == 25

          the third term
          - (8 / (x2sq * x2sq * x2sq))
          and beyond doesn't have a major contribution for large x.
        */
        return (one_div_root_pi / x) * (realT(1) - x2sq_r + x2sq_r * x2sq_r);
    }
    return std::exp(xsq) * boost::math::erfc(x);
}

// W(a, b) := std::exp(2 a b + b^2) erfc(a + b)
template<typename realT>
inline realT W(const realT a, const realT b) noexcept
{
    // std::exp(2 a b + b^2) erfc(a + b) ==
    //               std::exp(- a^2) std::exp((a + b)^2) erfc(a + b)
    return std::exp(- a * a) * expxsq_erfc(a + b);
}

} // gf11
#endif //GF11_MATHEMATICAL_UTILITY_H

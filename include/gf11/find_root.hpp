#ifndef GF11_FIND_ROOT_HPP
#define GF11_FIND_ROOT_HPP
#include "throw_exception.hpp"
#include "tolerance.hpp"
#include "type_traits.hpp"
#include <boost/math/tools/roots.hpp>
#include <utility>

namespace gf11
{

template<typename realT, typename F>
realT find_root(F&& f, realT low, realT high, tolerance<realT> tol,
                const std::uintmax_t iter_limit, const char* fname)
{
    static_assert(std::is_same<realT, return_type_of_t<F, realT>>::value, "");

    auto iteration = iter_limit;
    const auto result = boost::math::tools::toms748_solve(
        std::forward<F>(f), low, high, tol, iteration);

    if(iteration == iter_limit)
    {
        throw_exception<std::runtime_error>(fname, ": failed to find root "
            "(iteration=", iter_limit, ", low=", result.first, ", high=", result.second, ")");
    }
    return (result.first + result.second) / 2;
}

template<typename realT, typename F>
realT find_root(F&& f, realT low, realT high, realT low_value, realT high_value,
                tolerance<realT> tol, const std::uintmax_t iter_limit,
                const char* fname)
{
    static_assert(std::is_same<realT, return_type_of_t<F, realT>>::value, "");

    auto iteration = iter_limit;
    const auto result = boost::math::tools::toms748_solve(
        std::forward<F>(f), low, high, low_value, high_value, tol, iteration);

    if(iteration == iter_limit)
    {
        throw_exception<std::runtime_error>(fname, ": failed to find root "
            "(iteration=", iter_limit, ", low=", result.first, ", high=", result.second, ")");
    }
    return (result.first + result.second) / 2;
}

} // gf11
#endif// GF11_FIND_ROOT_HPP

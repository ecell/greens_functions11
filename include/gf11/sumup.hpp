#ifndef GF11_SUMUP_HPP
#define GF11_SUMUP_HPP
#include "type_traits.hpp"
#include <gsl/gsl_sum.h> // series acceleration method
#include <vector>
#include <memory>
#include <cstdint>
#include <cmath>

namespace gf11
{

// XXX: sumup_until_convergence depends on GSL, so the return type is fixed as double.

// F should take an unsigned integer and return a double.
template<typename F>
double sumup_until_convergence(F&& f, std::size_t max_terms, double tol)
{
    using function_type = typename std::remove_const<typename std::remove_reference<F>::type>::type;
    static_assert(std::is_same<return_type_of_t<function_type, std::uint32_t>, double>::value, "");

    // if(relative difference < tolerance) for 4 (by default) concecutive times,
    // it's considered as converged.
    constexpr std::uint32_t convergence_check = 4;

    double sum = f(0u);
    if(sum == 0.0)
    {
        return 0; // this assumption is only available with the currently-used greens functions.
    }

    std::vector<double> table(max_terms);
    table[0] = sum; // assign f(0) to table.at(0).

    std::uint32_t convergence_count = 0;
    for(std::size_t i=1; i<max_terms; ++i)
    {
        const double fi = f(i);
        table[i] = fi;
        sum += fi;

        if(std::abs(sum) * tol >= std::abs(fi)) // never omit `=`
        {
            ++convergence_count;
        }
        else
        {
            convergence_count = 0;
        }

        if(convergence_count >= convergence_check)
        {
            return sum; // converged. return right now.
        }
    }

    // didn't converge in [0, max_terms).
    // use series acceleration method provided by GSL.

    double error;

    std::unique_ptr<
        gsl_sum_levin_utrunc_workspace, decltype(&gsl_sum_levin_utrunc_free)
        > workspace(gsl_sum_levin_utrunc_alloc(max_terms), &gsl_sum_levin_utrunc_free);

    gsl_sum_levin_utrunc_accel(table.data(), table.size(), workspace.get(),
                               std::addressof(sum), std::addressof(error));

    if (std::abs(error) >= std::abs(sum * tol * 10.0))
    {
        throw_exception<std::runtime_error>("gf11::sumup: didn't converge(i<max)."
                " sum = ", sum, ", error = ", error, ", tolereance = ", tol,
                ", max_terms = ", max_terms);
    }
    return sum;
}

template<typename F>
double sumup_all(F&& f, std::size_t max_terms) noexcept
{
    using function_type = typename std::remove_const<typename std::remove_reference<F>::type>::type;
    static_assert(std::is_same<return_type_of_t<function_type, std::uint32_t>, double>::value, "");

    double sum = f(0u);
    if(sum == 0.0)
    {
        return 0;
    }
    for(std::size_t i=1; i<max_terms; ++i)
    {
        sum += f(i);
    }
    return sum;
}

} // gf11
#endif// GF11_SUMUP_HPP

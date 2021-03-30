#ifndef GF11_FACTORIAL_H
#define GF11_FACTORIAL_H
#include "index_sequence.hpp"
#include <limits>

namespace gf11
{

inline constexpr std::size_t max_factorial() noexcept {return 171;} // for double

namespace detail
{

template<typename valueT, std::size_t N>
struct compiletime_array
{
    static_assert(N > 0, "compiletime_array must have size");
    valueT value[N];
};

template<typename realT>
inline constexpr realT factorial_impl_rec(
    const realT current, const std::size_t i, const std::size_t maximum) noexcept
{
    return (i > maximum) ? current : // for max == 0 case
            factorial_impl_rec<realT>(current * static_cast<realT>(i), i+1, maximum);
}

template<typename realT>
inline constexpr realT factorial_impl(std::size_t i) noexcept
{
    return factorial_impl_rec<realT>(realT(1.0), 1, i);
}

template<typename realT, std::size_t ...Ns>
inline constexpr compiletime_array<realT, sizeof...(Ns)>
make_factorial_table(index_sequence<Ns...>) noexcept
{
    return {{factorial_impl<realT>(Ns)...}};
}

template<typename realT, std::size_t N>
struct factorial_table
{
    using value_type = realT;
    static constexpr std::size_t size = N;
    static constexpr compiletime_array<realT, N> table =
        make_factorial_table<realT>(make_index_sequence<N>{});

    static inline constexpr value_type get(std::size_t n) noexcept
    {
        return (n < N) ? factorial_table<realT, N>::table.value[n] :
                         std::numeric_limits<value_type>::infinity();
    }
};

template<typename realT, std::size_t N>
constexpr std::size_t factorial_table<realT, N>::size;
template<typename realT, std::size_t N>
constexpr compiletime_array<realT, N> factorial_table<realT, N>::table;

// factorial_r

template<typename realT>
inline constexpr realT factorial_r_impl(std::size_t i) noexcept
{
    return realT(1.) / factorial_table<realT, max_factorial()>::get(i);
}

template<typename realT, std::size_t ...Ns>
inline constexpr compiletime_array<realT, sizeof...(Ns)>
make_factorial_r_table(index_sequence<Ns...>) noexcept
{
    return {{factorial_r_impl<realT>(Ns)...}};
}

template<typename realT, std::size_t N>
struct factorial_r_table
{
    using value_type = realT;
    static constexpr std::size_t size = N;
    static constexpr compiletime_array<realT, N> table =
        make_factorial_r_table<realT>(make_index_sequence<N>{});

    static inline constexpr value_type get(std::size_t n) noexcept
    {
        return (n < N) ? factorial_r_table<realT, N>::table.value[n] : value_type(0.0);
    }
};
template<typename realT, std::size_t N>
constexpr std::size_t factorial_r_table<realT, N>::size;
template<typename realT, std::size_t N>
constexpr compiletime_array<realT, N> factorial_r_table<realT, N>::table;

} // detail

template<typename realT>
inline constexpr realT factorial(std::size_t n) noexcept
{
    return detail::factorial_table<realT, max_factorial()>::get(n);
}

template<typename realT>
inline constexpr realT factorial_r(std::size_t n) noexcept
{
    return detail::factorial_r_table<realT, max_factorial()>::get(n);
}

static_assert(factorial<double>(0) == 1.0, "");
static_assert(factorial<double>(1) == 1.0, "");
static_assert(factorial<double>(2) == 2.0, "");
static_assert(factorial<double>(3) == 2.0 * 3.0, "");
static_assert(factorial<double>(4) == 2.0 * 3.0 * 4.0, "");
static_assert(factorial<double>(5) == 2.0 * 3.0 * 4.0 * 5.0, "");

} // gf11
#endif// GF11_FACTORIAL_H

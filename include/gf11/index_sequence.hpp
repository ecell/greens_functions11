#ifndef GF11_INDEX_SEQUENCE_HPP
#define GF11_INDEX_SEQUENCE_HPP
#include <type_traits>
#include <cstddef>

namespace gf11
{

template<std::size_t ... Ns>
struct index_sequence
{
    using value_type = std::size_t;
    static constexpr std::size_t size = sizeof...(Ns);
};
template<std::size_t ... Ns>
constexpr std::size_t index_sequence<Ns...>::size;

namespace detail
{
template<typename T1, typename T2>
struct index_sequence_concatenator;

template<std::size_t ... v1, std::size_t ... v2>
struct index_sequence_concatenator<
    ::gf11::index_sequence<v1...>, ::gf11::index_sequence<v2...>>
{
    using type = ::gf11::index_sequence<v1..., v2...>;
};

// [First, Last]: both ends are included.
template<std::size_t First, std::size_t Last>
struct index_sequence_range
{
    static_assert(First < Last, "");

    // {0, 1, ..., N} := {0, 1, ..., N/2} ++ {N/2+1, N/2+2, ... N}
    using type = typename index_sequence_concatenator<
            typename index_sequence_range<First, (First + Last) / 2>::type,
            typename index_sequence_range<(First + Last) / 2 + 1, Last>::type
        >::type;
};
template<std::size_t I>
struct index_sequence_range<I, I>
{
    using type = index_sequence<I>;
};

// [First, Last): Last is not included.
template<std::size_t N>
struct index_sequence_generator
{
    using type = typename index_sequence_range<0, N-1>::type;
};
template<>
struct index_sequence_generator<0>
{
    using type = ::gf11::index_sequence<0>;
};
} // detail

// make {0, 1, ..., N-1}
template<std::size_t N>
using make_index_sequence = typename detail::index_sequence_generator<N>::type;

// tests
static_assert(std::is_same<make_index_sequence< 0>, index_sequence<0>>::value, "");
static_assert(std::is_same<make_index_sequence< 1>, index_sequence<0>>::value, "");
static_assert(std::is_same<make_index_sequence< 2>, index_sequence<0,1>>::value, "");
static_assert(std::is_same<make_index_sequence< 4>, index_sequence<0,1,2,3>>::value, "");
static_assert(std::is_same<make_index_sequence< 8>, index_sequence<0,1,2,3,4,5,6,7>>::value, "");
static_assert(std::is_same<make_index_sequence<16>, index_sequence<0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15>>::value, "");

} // gf11
#endif //GF11_INDEX_SEQUENCE_HPP

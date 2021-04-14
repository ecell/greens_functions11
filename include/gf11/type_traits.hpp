#ifndef GF11_TYPE_TRAITS_HPP
#define GF11_TYPE_TRAITS_HPP
#include <type_traits>

namespace gf11
{

// result_of is deprecated after C++17, use invoke_result instead.
#if __cplusplus >= 201703L
    template<typename F, typename ... Args>
    using return_type_of_t = std::invoke_result_t<F, Args...>;
#else
    template<typename F, typename ... Args>
    using return_type_of_t = typename std::result_of<F(Args...)>::type;
#endif

template<typename T>
using remove_cvref_t = typename std::remove_cv<
        typename std::remove_reference<T>::type
    >::type;

} // gf11
#endif// GF11_TYPE_TRAITS_HPP

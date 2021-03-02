#ifndef GF11_THROW_EXCEPTION_H
#define GF11_THROW_EXCEPTION_H
#include <sstream>
#include <string>
#include <utility>

namespace gf11
{

namespace detail
{

inline std::string concat_arguments_to_string_impl(std::ostringstream& oss)
{
    return oss.str();
}
template<typename T, typename ... Ts>
std::string concat_arguments_to_string_impl(
        std::ostringstream& oss, T&& head, Ts&& ... tail)
{
    oss << std::forward<T>(head);
    return concat_arguments_to_string_impl(oss, std::forward<Ts>(tail)...);
}
template<typename ... Ts>
std::string concat_arguments_to_string(Ts&& ... args)
{
    std::ostringstream oss;
    return concat_arguments_to_string_impl(oss, std::forward<Ts>(args)...);
}

} // detail

template<class Exception, typename ... Ts>
[[noreturn]] void throw_exception(Ts&& ... args)
{
    throw Exception(detail::concat_arguments_to_string(std::forward<Ts>(args)...));
}

} // gf11
#endif// GF11_THROW_EXCEPTION_H

#ifndef GF11_IPV_EVENT_KIND_HPP
#define GF11_IPV_EVENT_KIND_HPP
#include <cstdint>
#include <ostream>

namespace gf11
{

// This name is a bit confusing, but is kept for the sake of backward compatibility.
enum class GreensFunction : std::uint8_t
{
    IV_ESCAPE,
    IV_REACTION
};

inline std::ostream& operator<<(std::ostream& os, GreensFunction kind)
{
    if(kind == GreensFunction::IV_ESCAPE)
    {
        os << "iv_escape";
    }
    else if(kind == GreensFunction::IV_REACTION)
    {
        os << "iv_reaction";
    }
    else
    {
        os << "unknown";
    }
    return os;
}

} // gf11
#endif// GF11_TAGS_HPP

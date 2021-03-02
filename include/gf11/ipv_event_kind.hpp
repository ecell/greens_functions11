#ifndef GF11_IPV_EVENT_KIND_H
#define GF11_IPV_EVENT_KIND_H
#include <cstdint>
#include <ostream>

namespace gf11
{

enum class IpvEventKind : std::uint8_t
{
    ipv_escape,
    ipv_reaction,
};

inline std::ostream& operator<<(std::ostream& os, IpvEventKind k)
{
    if(k == IpvEventKind::ipv_escape)
    {
        os << "ipv_escape";
    }
    else if(k == IpvEventKind::ipv_reaction)
    {
        os << "ipv_reaction";
    }
    else
    {
        os << "unknown";
    }
    return os;
}

} // gf11
#endif// GF11_TAGS_H

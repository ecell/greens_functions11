#ifndef GF11_TOLERANCE_HPP
#define GF11_TOLERANCE_HPP
#include <cmath>

namespace gf11
{

template<typename realT>
struct tolerance
{
    using real_type = realT;

    tolerance(const real_type a, const real_type r) noexcept
        : abs_(a), rel_(r)
    {}

    bool operator()(const real_type lw, const real_type up) const noexcept
    {
        return (std::abs(up - lw) < abs_) || (std::abs(up / lw - real_type(1.)) < rel_);
    }

  private:
    real_type abs_, rel_;
};

} // gf11
#endif// GF11_TOLERANCE_HPP

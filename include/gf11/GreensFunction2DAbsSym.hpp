#ifndef GF11_2D_ABS_SYM_HPP
#define GF11_2D_ABS_SYM_HPP
#include <iosfwd>

namespace gf11
{

class GreensFunction2DAbsSym
{
  public:
    using real_type = double;

    GreensFunction2DAbsSym(const real_type D, const real_type a) noexcept
        : D_(D), a_(a)
    {}
    ~GreensFunction2DAbsSym() = default;

    real_type drawTime(const real_type rnd) const;
    real_type drawR   (const real_type rnd, const real_type t) const;

    real_type a() const noexcept {return a_;}
    real_type D() const noexcept {return D_;}

    static const char* name() noexcept {return "GreensFunction2DAbsSym";}

  private:

    real_type p_survival  (const real_type t)                    const;
    real_type p_int_r     (const real_type r, const real_type t) const;

  private:

    // H = 4.0: ~3e-5, 4.26: ~1e-6, 5.0: ~3e-7, 5.2: ~1e-7,
    //     5.6: ~1e-8, 6.0: ~1e-9
    static constexpr real_type CUTOFF()   noexcept {return 1e-10;}
    static constexpr real_type CUTOFF_H() noexcept {return 6.0;}

    real_type D_;
    real_type a_;
};


template<typename charT, typename traitsT>
std::basic_ostream<charT, traitsT>&
operator<<(std::basic_ostream<charT, traitsT>& os, const GreensFunction2DAbsSym& gf)
{
    os << gf.name() << "(D=" << gf.D() << ",a=" << gf.a() << ")";
    return os;
}

} // gf11

#if defined(GF11_HEADER_ONLY)
#include "GreensFunction2DAbsSym.cpp"
#endif

#endif// GF11_2D_ABS_SYM_HPP

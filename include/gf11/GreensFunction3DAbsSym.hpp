#ifndef GF11_3D_ABS_SYM_HPP
#define GF11_3D_ABS_SYM_HPP
#include <iosfwd>
#include <cmath>
#include <cstddef>

namespace gf11
{

class GreensFunction3DAbsSym
{
  public:
    using real_type = double;

  public:
    GreensFunction3DAbsSym(const real_type D, const real_type a) noexcept
        : D_(D), a_(a)
    {}
    ~GreensFunction3DAbsSym() = default;

    real_type drawTime(const real_type rnd) const;
    real_type drawR   (const real_type rnd, const real_type t) const;

    real_type a() const noexcept {return a_;}
    real_type D() const noexcept {return D_;}

    static const char* name() noexcept {return "GreensFunction3DAbsSym";}

  private:

    real_type p_survival  (const real_type t) const;
    real_type p_r_fourier (const real_type r, const real_type t) const;
    real_type p_int_r     (const real_type r, const real_type t) const noexcept;
    real_type p_int_r_free(const real_type r, const real_type t) const noexcept;

    static real_type ellipticTheta4Zero(const real_type q);

  private:

    static constexpr std::uintmax_t max_iteration_count() noexcept {return 100;}

    static constexpr real_type CUTOFF()     noexcept {return 1e-10;}
    static constexpr real_type CUTOFF_H()   noexcept {return 6.0;}
    static           real_type log_CUTOFF() noexcept {return std::log(CUTOFF());}

    real_type D_;
    real_type a_;
};

template<typename charT, typename traitsT>
std::basic_ostream<charT, traitsT>&
operator<<(std::basic_ostream<charT, traitsT>& os, const GreensFunction3DAbsSym& gf)
{
    os << "GreensFunction3DAbsSym("
       << "D=" << gf.D() << ", " << "a=" << gf.a() << ")";
    return os;
}

}// gf11

#if defined(GF11_HEADER_ONLY)
#include "GreensFunction3DAbsSym.cpp"
#endif

#endif//GF11_3D_ABS_SYM_HPP

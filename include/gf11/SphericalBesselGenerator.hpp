#ifndef GF11_SPHERICAL_BESSEL_GENERATOR_HPP
#define GF11_SPHERICAL_BESSEL_GENERATOR_HPP
#include <boost/math/special_functions/bessel.hpp>
#include <algorithm>
#include <vector>
#include <cmath>

namespace gf11
{
namespace detail
{

template<typename realT>
struct Table
{
    using real_type = realT;

    std::size_t N;
    real_type   x_start;
    real_type   delta_x;
    std::vector<real_type> y;
};

// We can start table interpolation from zero because there is
// no singularity in bessel_j for z>=0.

template<typename realT>
inline constexpr realT minz_j(const std::size_t /*n*/) noexcept {return 0.0;}
template<typename realT>
inline constexpr realT minz_y(const std::size_t /*n*/) noexcept {return 0.5;}

template<typename realT>
inline realT maxz_j(const std::size_t n) noexcept
{
    const realT z = (n * n + n + 1) / realT(1.221e-4);
    if(z >= realT(1000))
    {
        return (std::max<realT>)(1000, n * n);
    }
    return z;
}
template<typename realT>
inline realT maxz_y(const std::size_t n) noexcept
{
    // from gsl/special/bessel_y.c:
    //  else if(GSL_ROOT3_DBL_EPSILON * x > (l*l + l + 1.0)) {
    //     int status = gsl_sf_bessel_Ynu_asympx_e(l + 0.5, x, result);
    //     ...
    const realT z = (n * n + n + 1) / realT(6.06e-6);
    // ... but this is usually too big.
    if(z >= 2000)
    {
        return (std::max<realT>)(2000, n * n);
    }
    return z;
}

template<typename realT>
std::vector<realT> gen_z(const realT from, const realT to, const realT delta)
{
    const std::size_t N = std::ceil((to - from) / delta);

    std::vector<realT> table(N);
    realT z_val = from;
    for(std::size_t i=0; i<N; ++i)
    {
        table[i] = z_val;
        z_val += delta;
    }
    return table;
}

template<typename realT>
realT interpolate(const realT x, const Table<realT>& table)
{
    // hermite interpolation
    const realT h_inv = 1.0 / table.delta_x;

    realT i_real;
    const realT       x_lo  = std::modf((x - table.x_start) * h_inv, &i_real);
    const realT       x_hi  = 1.0 - x_lo;
    const std::size_t index = static_cast<size_t>(i_real * 2);

    const realT y_lo    = table.y[index    ];
    const realT ydot_lo = table.y[index + 1] * table.delta_x;
    const realT y_hi    = table.y[index + 2];
    const realT ydot_hi = table.y[index + 3] * table.delta_x;

    return x_hi * x_hi * (y_lo + x_lo * (2 * y_lo + ydot_lo)) +
           x_lo * x_lo * (y_hi + x_hi * (2 * y_hi - ydot_hi));
}

} // detail

template<typename realT>
class SphericalBesselGenerator
{
  public:
    using real_type  = realT;
    using table_type = detail::Table<real_type>;

    static constexpr std::size_t sj_table_min =  4;
    static constexpr std::size_t sj_table_max = 51;
    static constexpr std::size_t sy_table_min =  3;
    static constexpr std::size_t sy_table_max = 40;
    static constexpr std::size_t sjy_table_resolution = 35;

    static constexpr real_type delta =
        boost::math::constants::pi<real_type>() / sjy_table_resolution;

  public:

    SphericalBesselGenerator()  = default;
    ~SphericalBesselGenerator() = default;

    real_type j(const std::size_t n, const real_type z) const noexcept;
    real_type y(const std::size_t n, const real_type z) const noexcept;

    static constexpr std::size_t get_min_n_j() noexcept {return sj_table_min;}
    static constexpr std::size_t get_min_n_y() noexcept {return sy_table_min;}
    static constexpr std::size_t get_max_n_j() noexcept {return sj_table_max;}
    static constexpr std::size_t get_max_n_y() noexcept {return sy_table_max;}

  private:

    static std::vector<table_type> gen_sj_table();
    static std::vector<table_type> gen_sy_table();

    static real_type calc_j(const std::size_t n, const real_type z) noexcept
    {
        return boost::math::sph_bessel(n, z);
    }
    static real_type calc_y(const std::size_t n, const real_type z) noexcept
    {
        return boost::math::sph_neumann(n, z);
    }

    static real_type j_for_small_n(const std::size_t n, const real_type z) noexcept;
    static real_type y_for_small_n(const std::size_t n, const real_type z) noexcept;

    static table_type const& get_sj_table(const std::size_t n) noexcept {return sj_table_[n];}
    static table_type const& get_sy_table(const std::size_t n) noexcept {return sy_table_[n];}

  private:

    static std::vector<table_type> sj_table_;
    static std::vector<table_type> sy_table_;
};

template<typename realT>
constexpr std::size_t SphericalBesselGenerator<realT>::sj_table_min;
template<typename realT>
constexpr std::size_t SphericalBesselGenerator<realT>::sj_table_max;
template<typename realT>
constexpr std::size_t SphericalBesselGenerator<realT>::sy_table_min;
template<typename realT>
constexpr std::size_t SphericalBesselGenerator<realT>::sy_table_max;

template<typename realT>
std::vector<typename SphericalBesselGenerator<realT>::table_type>
SphericalBesselGenerator<realT>::sj_table_ =
    SphericalBesselGenerator<realT>::gen_sj_table();

template<typename realT>
std::vector<typename SphericalBesselGenerator<realT>::table_type>
SphericalBesselGenerator<realT>::sy_table_ =
    SphericalBesselGenerator<realT>::gen_sy_table();

template<typename realT>
typename SphericalBesselGenerator<realT>::real_type
SphericalBesselGenerator<realT>::j(
        const std::size_t n, const real_type z) const noexcept
{
    if     (n < sj_table_min) {return j_for_small_n(n, z);}
    else if(sj_table_max < n) {return calc_j(n, z);}

    const auto& table = get_sj_table(n);
    const real_type min_z = table.x_start + table.delta_x * 3;
    const real_type max_z = table.x_start + table.delta_x * (table.N - 3);

    if(min_z <= z && z < max_z)
    {
        return detail::interpolate(z, table);
    }
    return calc_j(n, z);
}

template<typename realT>
typename SphericalBesselGenerator<realT>::real_type
SphericalBesselGenerator<realT>::y(
        const std::size_t n, const real_type z) const noexcept
{
    if     (n < sy_table_min) {return y_for_small_n(n, z);}
    else if(sy_table_max < n) {return calc_y(n, z);}

    const auto& table = get_sy_table(n);
    const real_type min_z = table.x_start + table.delta_x * 3;
    const real_type max_z = table.x_start + table.delta_x * (table.N - 3);

    if(min_z <= z && z < max_z)
    {
        return detail::interpolate(z, table);
    }
    return calc_y(n, z);
}

template<typename realT>
typename SphericalBesselGenerator<realT>::real_type
SphericalBesselGenerator<realT>::j_for_small_n(
        const std::size_t n, const real_type z) noexcept
{
    assert(0 <= n && n <= 3);
    if(z == 0.0)
    {
        if(n == 0) {return 1.0;} else {return 0.0;}
    }
    if(n == 0)
    {
        return std::sin(z) / z;
    }

    const real_type sin_z = std::sin(z);
    const real_type cos_z = std::cos(z);
    const real_type z_r   = 1.0 / z;

    switch(n)
    {
        case 1:
        {
            return (sin_z * z_r - cos_z) * z_r;
        }
        case 2:
        {
            const real_type x3_zsq = 3.0 * z_r * z_r;
            return (x3_zsq - 1.0) * sin_z * z_r - x3_zsq * cos_z;
        }
        case 3:
        {
            const real_type x15_zsq = 15.0 * z_r * z_r;
            return ((x15_zsq - 6.0) * sin_z * z_r - (x15_zsq - 1.0) * cos_z) * z_r;
        }
        default: return std::numeric_limits<real_type>::quiet_NaN();
    }
}

template<typename realT>
typename SphericalBesselGenerator<realT>::real_type
SphericalBesselGenerator<realT>::y_for_small_n(
        const std::size_t n, const real_type z) noexcept
{
    assert(0 <= n && n <= 2);
    if(n == 0)
    {
        return -std::cos(z) / z;
    }

    const real_type sin_z = std::sin(z);
    const real_type cos_z = std::cos(z);
    const real_type z_r   = 1.0 / z;

    if(n == 1)
    {
        return -(cos_z * z_r + sin_z) * z_r;
    }
    else
    {
        const realT x3_zsq = 3.0 * z_r * z_r;
        return (1.0 - x3_zsq) * cos_z * z_r - x3_zsq * sin_z;
    }
}

template<typename realT>
std::vector<typename SphericalBesselGenerator<realT>::table_type>
SphericalBesselGenerator<realT>::gen_sj_table()
{
    const auto z_table = detail::gen_z(
            detail::minz_j<real_type>(sj_table_max),
            detail::maxz_j<real_type>(sj_table_max), delta);

    std::vector<table_type> sj_table(sj_table_max + 1);
    for(const auto z : z_table)
    {
        sj_table[sj_table_min-1].y.push_back(j_for_small_n(sj_table_min-1, z));
        sj_table[sj_table_min-1].y.push_back(0.0); // dummy (will not be used)
    }

    for(std::size_t n = sj_table_min; n <= sj_table_max; ++n)
    {
        sj_table[n].N       = z_table.size();
        sj_table[n].x_start = z_table.front();
        sj_table[n].delta_x = delta;

        sj_table[n].y.clear();
        sj_table[n].y.reserve(2 * z_table.size());
        for(std::size_t i=0; i<z_table.size(); ++i)
        {
            const real_type z = z_table[i];
            const real_type jn_1_z = sj_table[n-1].y[2*i];
            const real_type jn_z   = boost::math::sph_bessel(n, z);
            sj_table[n].y.push_back(jn_z);
            sj_table[n].y.push_back(jn_1_z - ((n + 1) / z) * jn_z);
        }
    }
    return sj_table;
}

template<typename realT>
std::vector<typename SphericalBesselGenerator<realT>::table_type>
SphericalBesselGenerator<realT>::gen_sy_table()
{
    const auto z_table = detail::gen_z(
            detail::minz_y<real_type>(sy_table_max),
            detail::maxz_y<real_type>(sy_table_max), delta);

    std::vector<table_type> sy_table(sy_table_max + 1);
    for(const auto z : z_table)
    {
        sy_table[sy_table_min-1].y.push_back(y_for_small_n(sy_table_min-1, z));
        sy_table[sy_table_min-1].y.push_back(0.0); // dummy (will not be used)
    }

    for(std::size_t n = sy_table_min; n <= sy_table_max; ++n)
    {
        sy_table[n].N       = z_table.size();
        sy_table[n].x_start = z_table.front();
        sy_table[n].delta_x = delta;

        sy_table[n].y.clear();
        sy_table[n].y.reserve(2 * z_table.size());
        for(std::size_t i=0; i<z_table.size(); ++i)
        {
            const real_type z = z_table[i];
            const real_type yn_1_z = sy_table[n-1].y[2*i];
            const real_type yn_z   = boost::math::sph_neumann(n, z);
            sy_table[n].y.push_back(yn_z);
            sy_table[n].y.push_back(yn_1_z - ((n + 1) / z) * yn_z);
        }
    }
    return sy_table;
}


} // gf11
#endif // GF11_SPHERICAL_BESSEL_GENERATOR_HPP

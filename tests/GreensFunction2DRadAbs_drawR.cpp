#include <gf11/GreensFunction2DRadAbs.hpp>
#include <greens_functions/GreensFunction2DRadAbs.hpp>
#include <iostream>
#include <chrono>

template<typename Rep, typename Period>
double format_as_seconds(const std::chrono::duration<Rep, Period>& dur)
{
    return std::chrono::duration_cast<std::chrono::microseconds>(dur).count() * 1e-6;
}

int main(int argc, char **argv)
{
    namespace gf98 = greens_functions;

    constexpr std::size_t N = 10000;
    constexpr double   drnd = 1.0 / N;

    const double D     = 1e-4;
    const double kf    = 1e-5;
    const double r0    = 0.020;
    const double sigma = 0.010;
    const double a     = 0.050;

    gf11::GreensFunction2DRadAbs gf_11(D, kf, r0, sigma, a);
    gf98::GreensFunction2DRadAbs gf_98(D, kf, r0, sigma, a);

    // -----------------------------------------------------------------------
    // drawR

    const auto t_example = gf_11.drawTime(0.5);

    std::vector<double> gf98_drawR(N);
    const auto gf98_drawR_start = std::chrono::system_clock::now();
    for(std::size_t i=0; i<N; ++i)
    {
        const double rnd = drnd * i;
        gf98_drawR[i] = gf_98.drawR(rnd, t_example);
    }
    const auto gf98_drawR_dur = std::chrono::system_clock::now() - gf98_drawR_start;

    std::vector<double> gf11_drawR(N);
    const auto gf11_drawR_start = std::chrono::system_clock::now();
    for(std::size_t i=0; i<N; ++i)
    {
        const double rnd = drnd * i;
        gf11_drawR[i] = gf_11.drawR(rnd, t_example);
    }
    const auto gf11_drawR_dur = std::chrono::system_clock::now() - gf11_drawR_start;

    std::cerr << "# drawR x" << N << "\n";
    std::cerr << "gf98: " << format_as_seconds(gf98_drawR_dur) << " [sec]\n";
    std::cerr << "gf11: " << format_as_seconds(gf11_drawR_dur) << " [sec]\n";

    std::cout << "# drawR\n";
    std::cout << "# rnd   gf11   gf98    abserr (gf11 - gf98)   relerr (|gf11/gf98 - 1|)\n";
    for(std::size_t i=0; i<N; ++i)
    {
        const double rnd = drnd * i;
        const double t11 = gf11_drawR[i];
        const double t98 = gf98_drawR[i];
        std::cout << rnd << ' ' << t11 << ' ' << t98 << ' '
                  << t11 - t98 << ' ' << std::abs(t11 / t98 - 1.0) << '\n';
    }
    return 0;
}

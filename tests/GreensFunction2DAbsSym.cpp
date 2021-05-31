#include <gf11/GreensFunction2DAbsSym.hpp>
#include <greens_functions/GreensFunction2DAbsSym.hpp>
#include <iostream>
#include <chrono>

int main(int argc, char **argv)
{
    namespace gf98 = greens_functions;

    constexpr std::size_t N = 100000;
    constexpr double   drnd = 1.0 / N;

    gf11::GreensFunction2DAbsSym gf_11(/* D = */1.0, /* a = */1.0);
    gf98::GreensFunction2DAbsSym gf_98(/* D = */1.0, /* a = */1.0);

    // -----------------------------------------------------------------------
    // drawTime

    std::vector<double> gf98_drawTime(N);
    const auto gf98_drawTime_start = std::chrono::system_clock::now();
    for(std::size_t i=0; i<N; ++i)
    {
        const double rnd = drnd * i;
        gf98_drawTime[i] = gf_98.drawTime(rnd);
    }
    const auto gf98_drawTime_dur = std::chrono::system_clock::now() - gf98_drawTime_start;

    std::vector<double> gf11_drawTime(N);
    const auto gf11_drawTime_start = std::chrono::system_clock::now();
    for(std::size_t i=0; i<N; ++i)
    {
        const double rnd = drnd * i;
        gf11_drawTime[i] = gf_11.drawTime(rnd);
    }
    const auto gf11_drawTime_dur = std::chrono::system_clock::now() - gf11_drawTime_start;

    std::cerr << "# drawTime x" << N << "\n";
    std::cerr << "gf98: " << std::chrono::duration_cast<std::chrono::microseconds>(gf98_drawTime_dur).count() / 1.0e6 << " [sec]\n";
    std::cerr << "gf11: " << std::chrono::duration_cast<std::chrono::microseconds>(gf11_drawTime_dur).count() / 1.0e6 << " [sec]\n";

    std::cout << "# drawTime\n";
    std::cout << "# rnd   gf11   gf98    abserr (gf11 - gf98)   relerr (|gf11/gf98 - 1|)\n";
    for(std::size_t i=0; i<N; ++i)
    {
        const double rnd = drnd * i;
        const double t11 = gf11_drawTime[i];
        const double t98 = gf98_drawTime[i];
        std::cout << rnd << ' ' << t11 << ' ' << t98 << ' '
                  << t11 - t98 << ' ' << std::abs(t11 / t98 - 1.0) << '\n';
    }

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
    std::cerr << "gf98: " << std::chrono::duration_cast<std::chrono::microseconds>(gf98_drawR_dur).count() / 1.0e6 << " [sec]\n";
    std::cerr << "gf11: " << std::chrono::duration_cast<std::chrono::microseconds>(gf11_drawR_dur).count() / 1.0e6 << " [sec]\n";

    std::cout << "# drawR\n";
    std::cout << "# rnd   gf11   gf98    abserr (gf11 - gf98)   relerr (|gf11/gf98 - 1|)\n";
    for(std::size_t i=0; i<N; ++i)
    {
        const double rnd = drnd * i;
        const double R11 = gf11_drawR[i];
        const double R98 = gf98_drawR[i];
        std::cout << rnd << ' ' << R11 << ' ' << R98 << ' '
                  << R11 - R98 << ' ' << std::abs(R11 / R98 - 1.0) << '\n';
    }
    return 0;
}

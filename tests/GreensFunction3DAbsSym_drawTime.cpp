#include <gf11/GreensFunction3DAbsSym.hpp>
#include <greens_functions/GreensFunction3DAbsSym.hpp>
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

    constexpr std::size_t N = 100000;
    constexpr double   drnd = 1.0 / N;

    gf11::GreensFunction3DAbsSym gf_11(/* D = */1.0, /* a = */1.0);
    gf98::GreensFunction3DAbsSym gf_98(/* D = */1.0, /* a = */1.0);

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
    std::cerr << "gf98: " << format_as_seconds(gf98_drawTime_dur) << " [sec]\n";
    std::cerr << "gf11: " << format_as_seconds(gf11_drawTime_dur) << " [sec]\n";

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

    return 0;
}

#include <gf11/GreensFunction3DRadInf.hpp>
#include <greens_functions/GreensFunction3DRadInf.hpp>
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

    const double D     = 1e-12;
    const double kf    = 1e-8;
    const double r0    = 5e-8;
    const double sigma = 1e-8;

    gf11::GreensFunction3DRadInf gf_11(D, kf, r0, sigma);
    gf98::GreensFunction3DRadInf gf_98(D, kf, r0, sigma);

    // -----------------------------------------------------------------------
    // drawR

    const auto t_example = gf_11.drawTime(0.5);

    std::vector<double> gf98_drawR(N);
    const auto gf98_drawR_start = std::chrono::system_clock::now();
    for(std::size_t i=1; i<N+1; ++i)
    {
        const double rnd = drnd * i;
        gf98_drawR[i] = gf_98.drawR(rnd, t_example);
    }
    const auto gf98_drawR_dur = std::chrono::system_clock::now() - gf98_drawR_start;

    std::vector<double> gf11_drawR(N);
    const auto gf11_drawR_start = std::chrono::system_clock::now();
    for(std::size_t i=1; i<N+1; ++i)
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
        const double R11 = gf11_drawR[i];
        const double R98 = gf98_drawR[i];
        std::cout << rnd << ' ' << R11 << ' ' << R98 << ' '
                  << R11 - R98 << ' ' << std::abs(R11 / R98 - 1.0) << '\n';
    }
    return 0;
}

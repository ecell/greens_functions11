#include <gf11/GreensFunction3DRadAbs.hpp>
#include <greens_functions/GreensFunction3DRadAbs.hpp>
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
    const double a     = 1e-7;

    gf11::GreensFunction3DRadAbs gf_11(D, kf, r0, sigma, a);
    gf98::GreensFunction3DRadAbs gf_98(D, kf, r0, sigma, a);

    // -----------------------------------------------------------------------
    // drawTheta

    const auto t_example = gf_11.drawTime(0.5);
    const auto r_example = gf_11.drawR(0.5, t_example);

    std::vector<double> gf98_drawTheta(N);
    const auto gf98_drawTheta_start = std::chrono::system_clock::now();
    for(std::size_t i=0; i<N; ++i)
    {
        const double rnd = drnd * i;
        gf98_drawTheta[i] = gf_98.drawTheta(rnd, t_example, r_example);
    }
    const auto gf98_drawTheta_dur = std::chrono::system_clock::now() - gf98_drawTheta_start;

    std::vector<double> gf11_drawTheta(N);
    const auto gf11_drawTheta_start = std::chrono::system_clock::now();
    for(std::size_t i=0; i<N; ++i)
    {
        const double rnd = drnd * i;
        gf11_drawTheta[i] = gf_11.drawTheta(rnd, t_example, r_example);
    }
    const auto gf11_drawTheta_dur = std::chrono::system_clock::now() - gf11_drawTheta_start;

    std::cerr << "# drawTheta x" << N << "\n";
    std::cerr << "gf98: " << format_as_seconds(gf98_drawTheta_dur) << " [sec]\n";
    std::cerr << "gf11: " << format_as_seconds(gf11_drawTheta_dur) << " [sec]\n";

    std::cout << "# drawTheta\n";
    std::cout << "# rnd   gf11   gf98    abserr (gf11 - gf98)   relerr (|gf11/gf98 - 1|)\n";
    for(std::size_t i=0; i<N; ++i)
    {
        const double rnd = drnd * i;
        const double R11 = gf11_drawTheta[i];
        const double R98 = gf98_drawTheta[i];
        std::cout << rnd << ' ' << R11 << ' ' << R98 << ' '
                  << R11 - R98 << ' ' << std::abs(R11 / R98 - 1.0) << '\n';
    }
    return 0;
}

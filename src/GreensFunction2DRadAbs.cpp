#include <gf11/GreensFunction2DRadAbs.hpp>
#include <greens_functions/GreensFunction2DRadAbs.hpp>
#include <iostream>

int main()
{
    namespace gf98 = greens_functions;

    const double D     = 1e-4;
    const double kf    = 1e-5;
    const double r0    = 0.020;
    const double sigma = 0.010;
    const double a     = 0.050;

    gf11::GreensFunction2DRadAbs gf_11(D, kf, r0, sigma, a);
    gf98::GreensFunction2DRadAbs gf_98(D, kf, r0, sigma, a);
    std::cout << "# rnd   gf11   gf98    abserr (gf11 - gf98)   relerr (|gf11/gf98 - 1|)" << std::endl;
    for(std::size_t i=1; i<10000; ++i)
    {
        const double rnd = 0.0001 * i;
        const double t98 = gf_98.drawTime(rnd);
        std::cout << rnd << ' ' << t98 << std::flush;
        const double t11 = gf_11.drawTime(rnd);
        std::cout << ' ' << t11 << ' '
                  << t11 - t98 << ' ' << std::abs(t11 / t98 - 1.0) << std::endl;
    }

    return 0;
}

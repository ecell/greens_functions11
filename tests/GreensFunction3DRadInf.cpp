#include <gf11/GreensFunction3DRadInf.hpp>
#include <greens_functions/GreensFunction3DRadInf.hpp>
#include <iostream>

int main()
{
    namespace gf98 = greens_functions;

    const double D     = 1e-12;
    const double kf    = 1e-8;
    const double r0    = 5e-8;
    const double sigma = 1e-8;

    gf11::GreensFunction3DRadInf gf_11(D, kf, r0, sigma);
    gf98::GreensFunction3DRadInf gf_98(D, kf, r0, sigma);
    std::cout << "# rnd   gf11   gf98    abserr (gf11 - gf98)   relerr (|gf11/gf98 - 1|)\n";
    for(std::size_t i=0; i<10000; ++i)
    {
        const double rnd = 0.0001 * i;
        const double t11 = gf_11.drawTime(rnd);
        const double t98 = gf_98.drawTime(rnd);
//         std::cout << rnd << ' ' << t98 << '\n';
        std::cout << rnd << ' ' << t11 << ' ' << t98 << ' '
                  << t11 - t98 << ' ' << std::abs(t11 / t98 - 1.0) << std::endl;
    }
    return 0;
}

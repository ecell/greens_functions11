#include <gf11/GreensFunction3DAbsSym.hpp>
#include <greens_functions/GreensFunction3DAbsSym.hpp>
#include <iostream>

int main()
{
    gf11::GreensFunction3DAbsSym<double>     gf_11(1.0, 1.0);
    greens_functions::GreensFunction3DAbsSym gf_98(1.0, 1.0);
    std::cout << "# rnd   gf11   gf98    abserr (gf11 - gf98)   relerr (|gf11/gf98 - 1|)\n";
    for(std::size_t i=0; i<10000; ++i)
    {
        const double rnd = 0.0001 * i;
        const double t11 = gf_11.drawTime(rnd);
        const double t98 = gf_98.drawTime(rnd);
        std::cout << rnd << ' ' << t11 << ' ' << t98 << ' '
                  << t11 - t98 << ' ' << std::abs(t11 / t98 - 1.0) << '\n';
    }
    return 0;
}

# gf11

Implementations of Green's functions for the enhanced Greens Functions Reaction Dynamics (eGFRD) in C++11.

- [x] 2DAbsSym
- [x] 2DAbsSym
- [ ] 2DRefWedgeAbs
- [x] 3DAbsSym
- [x] 3DRadInf
- [x] 3DRadAbs

## Build

### Prerequisites

- C++11 compliant compiler
- GSL
- Boost

### Usage as a header-only library

Define `GF11_HEADER_ONLY`. Then include it.

### Usage as a pre-built library

Compile the library via CMake giving `-DGF11_BUILD_LIBRARY=ON`.

### Testing

In `tests/`, there is a code to check the relative and absolute difference between
`greens_functions` and `greens_functions11`. The results are put on `stdout`.

It also compares the time took. The results are put on `stderr`.
To check correctly, it is recommended to turn `-O3 -march=native -mtune=native`
in `greens_functions`.

## Benchmark Results

The following micro benchmarks does benchmark and also error checking.
So the durations include the time took to store the data into vector and the
real efficiency might differ.

| function                     | gf98           | gf11                         |
|:-----------------------------|:---------------|:-----------------------------|
| 2DAbsSym.drawTime() x100'000 | 1.54186  [sec] | 0.367498 [sec]: ~4.2x faster |
| 2DAbsSym.drawR()    x100'000 | 1.56546  [sec] | 0.37903  [sec]: ~4.1x faster |
| 3DAbsSym.drawTime() x100'000 | 0.325472 [sec] | 0.118635 [sec]: ~2.7x faster |
| 3DAbsSym.drawR()    x100'000 | 0.667269 [sec] | 0.393544 [sec]: ~1.7x faster |

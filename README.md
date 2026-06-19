# Sturm

C++ implementation of Sturm sequences to bracket real roots of polynomials.

This repository provides:
- Library implementing Sturm sequences (sturm.h / sturm.cpp)
- CLI (sturm_app) and unit tests (Catch2)
- CMake build system and GitHub Actions CI

## Quick start

Build and run tests with CMake:

    mkdir -p build && cmake -S . -B build && cmake --build build
    ctest --test-dir build --output-on-failure

## Usage

Run the CLI and enter polynomial degree and coefficients when prompted:

    ./build/sturm_app

## Example:

```
Please enter order of polynomial: 6

Please enter coefficient number 6: -1
Please enter coefficient number 5: 1
Please enter coefficient number 4: 2
Please enter coefficient number 3: -2
Please enter coefficient number 2: 1
Please enter coefficient number 1: 2
Please enter coefficient number 0: -1

Sturm sequence for:
-1.000000 1.000000 2.000000 -2.000000 1.000000 2.000000 -1.000000 

-1.000000 1.000000 2.000000 -2.000000 1.000000 2.000000 -1.000000 
-1.000000 0.833333 1.333333 -1.000000 0.333333 0.333333 
-1.000000 0.965517 -0.620690 -2.137931 1.172414 
-1.000000 -0.667969 0.304687 -0.097656 
1.000000 0.764977 -0.502304 
1.000000 0.180000 
0.607600 

4 distinct real roots for x: -1.246980 -1.000000 0.445042 1.801938 
```

## References

- Original [C implementation](https://github.com/erich666/GraphicsGems/tree/master/gems/Sturm) (Graphics Gems) — D.G. Hook and P. R. McAree.
- Further reading: ["Roots of Polynomials"](docs/polyroots.pdf) (PDF) — Yan-Bin Jia

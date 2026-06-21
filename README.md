# Sturm

[![CI](https://github.com/viniciusjavs/Sturm/actions/workflows/ci.yml/badge.svg)](https://github.com/viniciusjavs/Sturm/actions/workflows/ci.yml)

Modern C++ implementation of **Sturm sequences** for bracketing the distinct real
roots of a real polynomial.

This repository provides:

- A small, header-clean library (`include/sturm/sturm.hpp`, `src/sturm.cpp`)
- A command-line front end (`sturm_app`)
- Unit tests (Catch2 v3) and a CMake build with cross-platform CI

## Highlights

- **Value semantics, no shared state** — a `Sturm` object is immutable after
  construction, so it is safe to reuse and to query concurrently from multiple
  threads.
- **No fixed degree limit** — coefficients live in `std::vector`, so there is no
  hard-coded maximum order and no fixed-buffer overflow risk.
- **Safe inputs** — trailing zero coefficients are trimmed to the true degree,
  and an identically-zero polynomial throws `sturm::invalid_polynomial`.
- **No I/O in the library** — all console interaction lives in the CLI.

## Quick start

Build and run the tests with CMake:

    cmake -S . -B build && cmake --build build
    ctest --test-dir build --output-on-failure

## CLI usage

Run the CLI and enter the polynomial degree and coefficients when prompted
(highest order first):

    ./build/sturm_app

### Example

Input for `-x^6 + x^5 + 2x^4 - 2x^3 + x^2 + 2x - 1` (degree `6`, then
coefficients `-1 1 2 -2 1 2 -1`):

```
Sturm sequence:
-1.000000 1.000000 2.000000 -2.000000 1.000000 2.000000 -1.000000
-1.000000 0.833333 1.333333 -1.000000 0.333333 0.333333
-1.000000 0.965517 -0.620690 -2.137931 1.172414
-1.000000 -0.667969 0.304687 -0.097656
1.000000 0.764977 -0.502304
1.000000 0.180000
0.607600

4 distinct real roots for x: -1.246980 -1.000000 0.445042 1.801938
```

## Library usage

```cpp
#include <iostream>

#include "sturm/sturm.hpp"

int main() {
    // p(x) = x^2 - 1, coefficients in ascending order: c0 + c1*x + c2*x^2
    sturm::Sturm poly(std::vector<double>{-1.0, 0.0, 1.0});

    std::cout << "degree: " << poly.degree() << '\n';
    for (double r : poly.real_roots()) std::cout << "root: " << r << '\n';

    // Or supply coefficients highest-order first:
    auto q = sturm::Sturm::from_descending(std::vector<double>{1.0, 0.0, -1.0});
}
```

Link against the `Sturm::sturm` target from CMake.

## Notes & limitations

- Roots are **distinct**: a repeated root (e.g. `(x - 1)^2`) is reported once,
  as Sturm's theorem counts distinct real roots, not multiplicities.
- Bisection isolates each root inside a bracket. In the rare case where a root
  sits exactly on the midpoint of a symmetric bracket (for example the root `0`
  of `x^3 - 25x`), the isolation can miss that root — a known limitation
  inherited from the original Graphics Gems routine.

## References

- Original [C implementation](https://github.com/erich666/GraphicsGems/tree/master/gems/Sturm)
  (Graphics Gems) — D. G. Hook and P. R. McAree.
- Further reading: ["Roots of Polynomials"](docs/polyroots.pdf) (PDF) — Yan-Bin Jia.

## License

The numerical method derives from *Graphics Gems* and is used under the Graphics
Gems license; the C++ port and additions are MIT licensed. See [LICENSE](LICENSE).

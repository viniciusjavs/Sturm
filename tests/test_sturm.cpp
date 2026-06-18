#define CATCH_CONFIG_MAIN
#include <catch2/catch.hpp>
#include <vector>

#include "sturm.h"

TEST_CASE("Quadratic x^2 - 1 has roots -1 and 1") {
    Sturm s;
    s.set_coeffs(2, std::vector<double>{-1.0, 0.0, 1.0});  // coef[0] = -1, coef[2] = 1
    auto roots = s.get_real_roots();
    REQUIRE(roots.size() == 2);
    // roots should be approximately -1 and 1 (order unspecified)
    REQUIRE(((fabs(roots[0] - (-1.0)) < 1e-9 && fabs(roots[1] - 1.0) < 1e-9) ||
             (fabs(roots[1] - (-1.0)) < 1e-9 && fabs(roots[0] - 1.0) < 1e-9)));
}

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

TEST_CASE("No real roots: x^2 + 1") {
    Sturm s;
    s.set_coeffs(2, std::vector<double>{1.0, 0.0, 1.0});  // x^2 + 1
    auto roots = s.get_real_roots();
    REQUIRE(roots.empty());
}

TEST_CASE("Repeated root: (x-1)^2") {
    Sturm s;
    s.set_coeffs(2, std::vector<double>{1.0, -2.0, 1.0});  // x^2 - 2x + 1
    auto roots = s.get_real_roots();
    REQUIRE(roots.size() == 1);
    REQUIRE(std::fabs(roots[0] - 1.0) < 1e-8);
}

TEST_CASE("Linear polynomial: x - 2") {
    Sturm s;
    s.set_coeffs(1, std::vector<double>{-2.0, 1.0});  // x - 2
    auto roots = s.get_real_roots();
    REQUIRE(roots.size() == 1);
    REQUIRE(std::fabs(roots[0] - 2.0) < 1e-8);
}

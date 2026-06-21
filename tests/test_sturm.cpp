#include <algorithm>
#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <cmath>
#include <thread>
#include <vector>

#include "sturm/sturm.hpp"

using sturm::Sturm;

namespace {

// Sort roots and compare against expected values (order-independent).
void check_roots(std::vector<double> got, std::vector<double> expected, double tol = 1e-7) {
    std::sort(got.begin(), got.end());
    std::sort(expected.begin(), expected.end());
    REQUIRE(got.size() == expected.size());
    for (std::size_t i = 0; i < expected.size(); ++i)
        REQUIRE(got[i] == Catch::Approx(expected[i]).margin(tol));
}

}  // namespace

TEST_CASE("Quadratic x^2 - 1 has roots -1 and 1") {
    Sturm s(std::vector<double>{-1.0, 0.0, 1.0});
    check_roots(s.real_roots(), {-1.0, 1.0});
}

TEST_CASE("No real roots: x^2 + 1") {
    Sturm s(std::vector<double>{1.0, 0.0, 1.0});
    REQUIRE(s.real_roots().empty());
}

TEST_CASE("Repeated root collapses to a single distinct root: (x-1)^2") {
    Sturm s(std::vector<double>{1.0, -2.0, 1.0});
    check_roots(s.real_roots(), {1.0}, 1e-7);
}

TEST_CASE("Linear polynomial: x - 2") {
    Sturm s(std::vector<double>{-2.0, 1.0});
    check_roots(s.real_roots(), {2.0});
}

TEST_CASE("Degree-6 example from the README has four distinct roots") {
    // -1 + 2x + x^2 - 2x^3 + 2x^4 + x^5 - x^6  (ascending coefficients)
    Sturm s(std::vector<double>{-1.0, 2.0, 1.0, -2.0, 2.0, 1.0, -1.0});
    auto roots = s.real_roots();
    REQUIRE(roots.size() == 4);
    for (double r : roots) REQUIRE(s(r) == Catch::Approx(0.0).margin(1e-6));
}

TEST_CASE("from_descending matches the ascending constructor") {
    // x^2 - 1: ascending {-1, 0, 1}, descending {1, 0, -1}.
    Sturm asc(std::vector<double>{-1.0, 0.0, 1.0});
    Sturm desc = Sturm::from_descending(std::vector<double>{1.0, 0.0, -1.0});
    check_roots(desc.real_roots(), asc.real_roots());
}

TEST_CASE("Evaluation uses Horner's method") {
    Sturm s(std::vector<double>{-6.0, 11.0, -6.0, 1.0});  // (x-1)(x-2)(x-3)
    REQUIRE(s(0.0) == Catch::Approx(-6.0));
    REQUIRE(s(2.0) == Catch::Approx(0.0).margin(1e-12));
    check_roots(s.real_roots(), {1.0, 2.0, 3.0});
}

TEST_CASE("Degenerate and edge-case inputs are handled safely") {
    SECTION("identically-zero polynomial throws") {
        REQUIRE_THROWS_AS(Sturm(std::vector<double>{0.0, 0.0, 0.0}), sturm::invalid_polynomial);
    }
    SECTION("empty coefficient list throws") {
        REQUIRE_THROWS_AS(Sturm(std::vector<double>{}), sturm::invalid_polynomial);
    }
    SECTION("non-zero constant has no roots") {
        Sturm s(std::vector<double>{5.0});
        REQUIRE(s.degree() == 0);
        REQUIRE(s.real_roots().empty());
    }
    SECTION("trailing zero coefficients are trimmed to the true degree") {
        Sturm s(std::vector<double>{2.0, 3.0, 0.0, 0.0});  // really 2 + 3x
        REQUIRE(s.degree() == 1);
        check_roots(s.real_roots(), {-2.0 / 3.0});
    }
}

TEST_CASE("Roots far outside [-1, 1] are bracketed correctly") {
    // (x + 4)(x - 2)(x - 3) = x^3 - x^2 - 14x + 24
    Sturm s(std::vector<double>{24.0, -14.0, -1.0, 1.0});
    check_roots(s.real_roots(), {-4.0, 2.0, 3.0}, 1e-6);
}

TEST_CASE("Repeated calls on the same object are stable (reentrant)") {
    // Roots outside [-1, 1] force the bracket to widen; a buggy implementation
    // that kept the widened bracket as state would drift on the second call.
    Sturm s(std::vector<double>{24.0, -14.0, -1.0, 1.0});
    auto first = s.real_roots();
    auto second = s.real_roots();
    REQUIRE(first == second);
}

TEST_CASE("Concurrent queries are thread-safe") {
    const Sturm s(std::vector<double>{24.0, -14.0, -1.0, 1.0});
    const auto expected = s.real_roots();

    std::vector<std::thread> threads;
    std::vector<std::vector<double>> results(8);
    for (std::size_t i = 0; i < results.size(); ++i)
        threads.emplace_back([&, i] { results[i] = s.real_roots(); });
    for (auto& t : threads) t.join();

    for (const auto& r : results) REQUIRE(r == expected);
}

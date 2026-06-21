#ifndef STURM_STURM_HPP
#define STURM_STURM_HPP

// Sturm sequences for bracketing the distinct real roots of a real polynomial.
//
// Modern C++ port of the routine from "Graphics Gems" (Academic Press, 1990):
//   "Using Sturm Sequences to Bracket Real Roots of Polynomial Equations"
//   by D. G. Hook and P. R. McAree.
// See LICENSE for attribution and licensing.

#include <span>
#include <stdexcept>
#include <vector>

namespace sturm {

// Thrown when a polynomial cannot be processed (e.g. it is identically zero).
class invalid_polynomial : public std::invalid_argument {
   public:
    using std::invalid_argument::invalid_argument;
};

// An immutable real polynomial.
//
// Coefficients are stored in ascending order:
//     p(x) = c[0] + c[1]*x + c[2]*x^2 + ... + c[n]*x^n
// Trailing (highest-order) zero coefficients are trimmed on construction so
// that degree() always reflects the true degree.
//
// The type is a value type: copyable, movable, and free of shared state, so
// instances are safe to reuse and to query concurrently from multiple threads.
class Sturm {
   public:
    // Construct from coefficients in ascending order (c[0] is the constant
    // term). Throws invalid_polynomial if every coefficient is zero.
    explicit Sturm(std::vector<double> ascending_coeffs);

    // Construct from coefficients in descending order (c[n] first, c[0] last),
    // matching the common "write the polynomial left to right" convention.
    [[nodiscard]] static Sturm from_descending(std::span<const double> descending_coeffs);

    // Degree of the polynomial (>= 0; a non-zero constant has degree 0).
    [[nodiscard]] int degree() const noexcept;

    // Evaluate p(x) using Horner's method.
    [[nodiscard]] double operator()(double x) const noexcept;

    // The Sturm sequence p0, p1, ..., pm. Each row holds one polynomial's
    // coefficients in descending order; row 0 is the polynomial itself.
    [[nodiscard]] std::vector<std::vector<double>> sequence() const;

    // The distinct real roots, in ascending order. Repeated roots are reported
    // once (Sturm's theorem counts distinct roots). Returns an empty vector when
    // there are no real roots.
    [[nodiscard]] std::vector<double> real_roots() const;

   private:
    std::vector<double> coeffs_;  // ascending order, trailing zeros trimmed
};

}  // namespace sturm

#endif  // STURM_STURM_HPP

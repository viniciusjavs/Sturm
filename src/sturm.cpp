// Sturm sequence root bracketing — modern C++ port.
//
// Numerical method from "Graphics Gems" (Academic Press, 1990):
//   "Using Sturm Sequences to Bracket Real Roots of Polynomial Equations"
//   by D. G. Hook and P. R. McAree.
//
// The numerical kernel below preserves the original algorithm exactly; the
// modernization is structural: std::vector instead of fixed C arrays (no degree
// cap, no out-of-bounds risk), no shared mutable state (reentrant and
// thread-safe), const-correctness, and no I/O in the library.

#include "sturm/sturm.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>

namespace sturm {
namespace {

constexpr double rel_error = 1.0e-14;     // smallest relative error we want
constexpr int max_pow = 32;               // max power of 10 we search to
constexpr int max_it = 800;               // max number of iterations
constexpr double small_enough = 1.0e-12;  // coefficients below this are zero

// One polynomial in the chain, coefficients ascending: coef[0..ord].
struct Poly {
    int ord = 0;
    std::vector<double> coef;
};

// Evaluate a chain polynomial at x via Horner's method.
double eval(const Poly& p, double x) noexcept {
    double f = p.coef[p.ord];
    for (int i = p.ord - 1; i >= 0; --i) f = x * f + p.coef[i];
    return f;
}

// r = remainder of u / v (Graphics Gems convention). Assumes |lead(v)| == 1.
// Returns the degree of r (0 when r is constant).
int poly_mod(const Poly& u, const Poly& v, Poly& r) {
    r.coef.assign(u.coef.begin(), u.coef.begin() + (u.ord + 1));

    if (v.coef[v.ord] < 0.0) {
        for (int k = u.ord - v.ord - 1; k >= 0; k -= 2) r.coef[k] = -r.coef[k];
        for (int k = u.ord - v.ord; k >= 0; --k)
            for (int j = v.ord + k - 1; j >= k; --j)
                r.coef[j] = -r.coef[j] - r.coef[v.ord + k] * v.coef[j - k];
    } else {
        for (int k = u.ord - v.ord; k >= 0; --k)
            for (int j = v.ord + k - 1; j >= k; --j) r.coef[j] -= r.coef[v.ord + k] * v.coef[j - k];
    }

    int k = v.ord - 1;
    while (k >= 0 && std::fabs(r.coef[k]) < small_enough) {
        r.coef[k] = 0.0;
        --k;
    }
    r.ord = (k < 0) ? 0 : k;
    return r.ord;
}

// Build the Sturm chain for a polynomial of degree >= 1.
std::vector<Poly> build_chain(const std::vector<double>& coeffs) {
    const int order = static_cast<int>(coeffs.size()) - 1;

    std::vector<Poly> seq;
    seq.reserve(coeffs.size());

    Poly p0;
    p0.ord = order;
    p0.coef = coeffs;
    seq.push_back(std::move(p0));

    // p1 = p' normalised so that |leading coefficient| == 1.
    Poly p1;
    p1.ord = order - 1;
    p1.coef.resize(order);
    const double lead = std::fabs(seq[0].coef[order] * order);
    for (int i = 1; i <= order; ++i) p1.coef[i - 1] = seq[0].coef[i] * i / lead;
    seq.push_back(std::move(p1));

    for (;;) {
        Poly r;
        const Poly& u = seq[seq.size() - 2];
        const Poly& v = seq[seq.size() - 1];
        if (poly_mod(u, v, r) > 0) {
            // Reverse the sign and normalise (leading coefficient -> -1).
            const double g = -std::fabs(r.coef[r.ord]);
            for (int i = r.ord; i >= 0; --i) r.coef[i] /= g;
            seq.push_back(std::move(r));
        } else {
            r.coef[0] = -r.coef[0];  // reverse the sign of the final constant
            seq.push_back(std::move(r));
            break;
        }
    }
    return seq;
}

// Sign changes among the chain's leading terms as x -> +/- infinity.
// neg_inf == false gives the count at +infinity, true at -infinity.
int changes_at_infinity(const std::vector<Poly>& seq, bool neg_inf) {
    const auto lead = [neg_inf](const Poly& p) {
        const double l = p.coef[p.ord];
        return (neg_inf && (p.ord & 1)) ? -l : l;
    };
    int changes = 0;
    double lf = lead(seq[0]);
    for (std::size_t i = 1; i < seq.size(); ++i) {
        const double f = lead(seq[i]);
        if (lf == 0.0 || lf * f < 0) ++changes;
        lf = f;
    }
    return changes;
}

// Number of sign changes in the chain evaluated at a.
int num_changes(const std::vector<Poly>& seq, double a) {
    int changes = 0;
    double lf = eval(seq[0], a);
    for (std::size_t i = 1; i < seq.size(); ++i) {
        const double f = eval(seq[i], a);
        if (lf == 0.0 || lf * f < 0) ++changes;
        lf = f;
    }
    return changes;
}

// Modified regula-falsi on the single root of p in [a, b]. Returns true and
// writes the root to val on success, false if it cannot converge.
bool mod_rf(const Poly& p, double a, double b, double& val) {
    const int ord = p.ord;
    const std::vector<double>& coef = p.coef;

    double fa = coef[ord];
    double fb = fa;
    for (int i = ord - 1; i >= 0; --i) {
        fa = a * fa + coef[i];
        fb = b * fb + coef[i];
    }

    if (fa * fb > 0.0) return false;
    if (std::fabs(fa) < rel_error) {
        val = a;
        return true;
    }
    if (std::fabs(fb) < rel_error) {
        val = b;
        return true;
    }

    double lfx = fa;
    for (int its = 0; its < max_it; ++its) {
        const double x = (fb * a - fa * b) / (fb - fa);
        double fx = coef[ord];
        for (int i = ord - 1; i >= 0; --i) fx = x * fx + coef[i];

        if (std::fabs(x) > rel_error) {
            if (std::fabs(fx / x) < rel_error) {
                val = x;
                return true;
            }
        } else if (std::fabs(fx) < rel_error) {
            val = x;
            return true;
        }

        if (fa * fx < 0) {
            b = x;
            fb = fx;
            if (lfx * fx > 0) fa /= 2;
        } else {
            a = x;
            fa = fx;
            if (lfx * fx > 0) fb /= 2;
        }
        lfx = fx;
    }
    return false;
}

// Isolate roots in [min, max] into the (exactly sized) span via bisection.
void bisect(const std::vector<Poly>& seq, double min, double max, int atmin, int atmax,
            std::span<double> roots) {
    const int nroot = atmin - atmax;
    double mid = 0.0;
    int its = 0;

    if (nroot == 1) {
        if (mod_rf(seq[0], min, max, roots[0])) return;

        for (its = 0; its < max_it; ++its) {
            mid = (min + max) / 2;
            const int atmid = num_changes(seq, mid);
            if (std::fabs(mid) > rel_error) {
                if (std::fabs((max - min) / mid) < rel_error) {
                    roots[0] = mid;
                    return;
                }
            } else if (std::fabs(max - min) < rel_error) {
                roots[0] = mid;
                return;
            }
            if ((atmin - atmid) == 0)
                min = mid;
            else
                max = mid;
        }
        roots[0] = mid;  // failed to converge: best estimate
        return;
    }

    for (its = 0; its < max_it; ++its) {
        mid = (min + max) / 2;
        const int atmid = num_changes(seq, mid);
        const int n1 = atmin - atmid;
        const int n2 = atmid - atmax;
        if (n1 != 0 && n2 != 0) {
            const auto split =
                std::clamp<std::size_t>(static_cast<std::size_t>(n1), 0, roots.size());
            bisect(seq, min, mid, atmin, atmid, roots.subspan(0, split));
            bisect(seq, mid, max, atmid, atmax, roots.subspan(split));
            return;
        }
        if (n1 == 0)
            min = mid;
        else
            max = mid;
    }

    // Roots too close together: report the shared midpoint for the interval.
    for (double& r : roots) r = mid;
}

}  // namespace

Sturm::Sturm(std::vector<double> ascending_coeffs) : coeffs_(std::move(ascending_coeffs)) {
    while (!coeffs_.empty() && coeffs_.back() == 0.0) coeffs_.pop_back();
    if (coeffs_.empty()) throw invalid_polynomial("Sturm: polynomial is identically zero");
}

Sturm Sturm::from_descending(std::span<const double> descending_coeffs) {
    return Sturm(std::vector<double>(descending_coeffs.rbegin(), descending_coeffs.rend()));
}

int Sturm::degree() const noexcept { return static_cast<int>(coeffs_.size()) - 1; }

double Sturm::operator()(double x) const noexcept {
    double f = 0.0;
    for (auto it = coeffs_.rbegin(); it != coeffs_.rend(); ++it) f = x * f + *it;
    return f;
}

std::vector<std::vector<double>> Sturm::sequence() const {
    std::vector<std::vector<double>> rows;
    if (degree() < 1) {
        rows.push_back({coeffs_.rbegin(), coeffs_.rend()});
        return rows;
    }

    const std::vector<Poly> seq = build_chain(coeffs_);
    rows.reserve(seq.size());
    for (const Poly& p : seq) {
        std::vector<double> row;
        row.reserve(static_cast<std::size_t>(p.ord) + 1);
        for (int j = p.ord; j >= 0; --j) row.push_back(p.coef[j]);
        rows.push_back(std::move(row));
    }
    return rows;
}

std::vector<double> Sturm::real_roots() const {
    if (degree() < 1) return {};  // a non-zero constant has no roots

    const std::vector<Poly> seq = build_chain(coeffs_);

    int atmax = changes_at_infinity(seq, /*neg_inf=*/false);
    int atmin = changes_at_infinity(seq, /*neg_inf=*/true);
    if (atmin - atmax <= 0) return {};

    // Widen the bracket [-1, 1] until it captures every root. These are locals,
    // so repeated calls (and concurrent calls) never interfere.
    double min = -1.0;
    double max = 1.0;
    int nchanges = num_changes(seq, min);
    for (int i = 0; nchanges != atmin && i != max_pow; ++i) {
        min *= 10.0;
        nchanges = num_changes(seq, min);
    }
    if (nchanges != atmin) atmin = nchanges;

    nchanges = num_changes(seq, max);
    for (int i = 0; nchanges != atmax && i != max_pow; ++i) {
        max *= 10.0;
        nchanges = num_changes(seq, max);
    }
    if (nchanges != atmax) atmax = nchanges;

    const int nroots = atmin - atmax;
    if (nroots <= 0) return {};

    std::vector<double> roots(static_cast<std::size_t>(nroots), 0.0);
    bisect(seq, min, max, atmin, atmax, roots);
    return roots;
}

}  // namespace sturm

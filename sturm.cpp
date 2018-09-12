// A class based C++14 program fork of the one located at:
// https://webdocs.cs.ualberta.ca/~graphics/books/GraphicsGems/gems/Sturm/
// Author: @vjavs
// Timestamp: 14 Mar 2017

/*
  Using Sturm Sequences to Bracket Real Roots of Polynomial Equations by D.G.
  Hook and P.R. McAree from "Graphics Gems", Academic Press, 1990
 */

#include <iostream>
#include <cmath>
#include <vector>

#include "sturm.h"

// Some useful constants
constexpr double rel_error {1.0e-14};    // smallest relative error we want
constexpr int max_pow = 32;              // max power of 10 we wish to search to
constexpr int max_it = 800;              // max number of iterations
constexpr double small_enough {1.0e-12}; // a coefficient smaller than
                                         // small_enough is considered to be
                                         // zero (0.0).

void Sturm::set_parameters()
/*
  Input the parameters
 */
{
    std::cout << "Please enter order of polynomial: ";
    std::cin >> order;
    std::cout << '\n';
    for (auto i = order; i >= 0; --i) {
        std::cout << "Please enter coefficient number " << i << ": ";
        std::cin >> sturm_seq[0].coef[i];
    }
    std::cout << '\n';
}

std::vector<std::vector<double>> Sturm::get_sturm_sequence()
/*
  Build the Sturm sequence and returns a vector of coefficients
 */
{
    std::vector<std::vector<double>> seq;
    std::vector<double> coefs;
    
    num_poly = build_sturm();
    for (auto i = order; i >= 0; --i)
	coefs.push_back(sturm_seq[0].coef[i]);
    seq.push_back(coefs);
    for (auto i = 0; i <= num_poly; ++i) {
	coefs.clear();
        for (auto j = sturm_seq[i].ord; j >= 0; --j)
	    coefs.push_back(sturm_seq[i].coef[j]);
	seq.push_back(coefs);
    }
    return seq;
}

void Sturm::show_sturm_sequence(const std::vector<std::vector<double>> &seq)
{
    std::cout << "Sturm sequence for:\n";
    std::cout << std::fixed;
    auto first = true;
    for (const auto &poly : seq) {
	for(const auto &coef : poly)
	    std::cout << coef << ' ';
	if (first) { std::cout << "\n"; first = false;}
        std::cout << "\n";
    }
}

int Sturm::build_sturm()
/*
  Build up a sturm sequence for a polynomial in smat, returning the number of
  polynomials in the sequence
 */
{
    double f, *fp, *fc;
    Poly *sp;

    sturm_seq[0].ord = order;
    sturm_seq[1].ord = order - 1;

    // Calculate the derivative and normalise the leading coefficient.
    f = fabs(sturm_seq[0].coef[order] * order);
    fp = sturm_seq[1].coef;
    fc = sturm_seq[0].coef + 1;
    for (auto i = 1; i <= order; i++)
        *fp++ = *fc++ *i / f;

    // Construct the rest of the Sturm sequence
    for (sp = sturm_seq + 2; modp(sp - 2, sp - 1, sp); sp++) {
        // Reverse the sign and normalise
        f = -fabs(sp->coef[sp->ord]);
        for (fp = &sp->coef[sp->ord]; fp >= sp->coef; fp--)
            *fp /= f;
    }

    sp->coef[0] = -sp->coef[0]; // reverse the sign

    return sp - sturm_seq;
}

int Sturm::modp(Poly *u, Poly *v, Poly *r)
/*
  Calculates the modulus of u(x) / v(x) leaving it in r, it  returns 0 if r(x)
  is a constant. Note: this function assumes the leading coefficient of v is 1
  or -1
 */
{
    int k, j;
    double *nr, *end, *uc;

    nr = r->coef;
    end = &u->coef[u->ord];

    uc = u->coef;
    while (uc <= end)
        *nr++ = *uc++;

    if (v->coef[v->ord] < 0.0) {
        for (k = u->ord - v->ord - 1; k >= 0; k -= 2)
            r->coef[k] = -r->coef[k];
        for (k = u->ord - v->ord; k >= 0; k--)
            for (j = v->ord + k - 1; j >= k; j--)
                r->coef[j] = -r->coef[j] - r->coef[v->ord + k] * v->coef[j - k];
    }
    else {
        for (k = u->ord - v->ord; k >= 0; k--)
            for (j = v->ord + k - 1; j >= k; j--)
                r->coef[j] -= r->coef[v->ord + k] * v->coef[j - k];
    }

    k = v->ord - 1;
    while (k >= 0 && fabs(r->coef[k]) < small_enough) {
        r->coef[k] = 0.0;
        k--;
    }

    r->ord = (k < 0) ? 0 : k;

    return r->ord;
}

std::vector<double> Sturm::get_real_roots()
{
    // Get the number of real roots
    nroots = num_roots();
    if (nroots == 0) {
	std::cout << "solve: no real roots\n";
        exit (0); ///exception?
    }
    
    // Calculate the bracket that the roots live in
    nchanges = num_changes(min);
    for (auto i = 0; nchanges != atmin && i != max_pow; ++i) {
        min *= 10.0;
        nchanges = num_changes(min);
    }
    if (nchanges != atmin) {
	std::cout << "solve: unable to bracket all negative roots\n";
        atmin = nchanges;
    }

    nchanges = num_changes(max);
    for (auto i = 0; nchanges != atmax && i != max_pow; ++i) {
        max *= 10.0;
        nchanges = num_changes(max);
    }
    if (nchanges != atmax) {
	std::cout << "solve: unable to bracket all positive roots\n";
        atmax = nchanges;
    }
    nroots = atmin - atmax;

    // Perform the bisection.
    bisect(min, max, atmin, atmax, roots);

    // fills vector of roots
    std::vector<double> roots_vec;
    for (auto i = 0; i != nroots; ++i)
        roots_vec.push_back(roots[i]);

    return roots_vec;
}

int Sturm::num_roots()
/*
 Return the number of distinct real roots of the polynomial* described in sseq.
 */
{
    int atposinf = 0, atneginf = 0;
    Poly *s;
    double f, lf;

    // Changes at positive infinity
    lf = sturm_seq[0].coef[sturm_seq[0].ord];

    for (s = sturm_seq + 1; s <= sturm_seq + num_poly; ++s) {
        f = s->coef[s->ord];
        if (lf == 0.0 || lf * f < 0) ++atposinf;
        lf = f;
    }

    // Changes at negative infinity
    if (sturm_seq[0].ord & 1) lf = -sturm_seq[0].coef[sturm_seq[0].ord];
    else lf = sturm_seq[0].coef[sturm_seq[0].ord];

    for (s = sturm_seq + 1; s <= sturm_seq + num_poly; ++s) {
        if (s->ord & 1) f = -s->coef[s->ord];
        else f = s->coef[s->ord];
        if (lf == 0.0 || lf * f < 0) ++atneginf;
        lf = f;
    }

    atmin = atneginf;
    atmax = atposinf;

    return atneginf - atposinf;
}

int Sturm::num_changes(const double& a)
/*
  Return the number of sign changes in the Sturm sequence in  sseq at the value 
  a.
 */
{
    double f, lf;
    Poly *s;
    int changes = 0;

    lf = eval_poly(sturm_seq[0].ord, sturm_seq[0].coef, a);

    for (s = sturm_seq + 1; s <= sturm_seq + num_poly; s++) {
        f = eval_poly(s->ord, s->coef, a);
        if (lf == 0.0 || lf * f < 0) ++changes;
        lf = f;
    }

    return changes;
}

double Sturm::eval_poly(int ord, double *coef, double x)
{
/*
  Evaluate polynomial defined in coef returning its value.
 */
    double *fp, f;
    fp = &coef[ord];
    f = *fp;
    for (fp--; fp >= coef; fp--)
        f = x * f + *fp;
    return f;
}

void Sturm::bisect(double min, double  max, const int &atmin, const int &atmax, double *roots)
/*
  Uses a bisection based on the sturm sequence for the polynomial described in
  sseq to isolate intervals in which roots occur, the roots are returned in the
  roots array in order of magnitude.
 */
{
    double mid;
    int n1 = 0, n2 = 0, its, atmid, nroot;

    if ((nroot = atmin - atmax) == 1) {

        // First try a less expensive technique.
        if (mod_rf(sturm_seq->ord, sturm_seq->coef, min, max, &roots[0]))
            return;

        // if we get here we have to evaluate the root the hard way by using the
        // Sturm sequence.
        for (its = 0; its < max_it; its++) {
            mid = (min + max) / 2;
            atmid = num_changes(mid);
            if (fabs(mid) > rel_error) {
                if (fabs((max - min) / mid) < rel_error) {
                    roots[0] = mid;
                    return;
                }
            }
            else if (fabs(max - min) < rel_error) {
                roots[0] = mid;
                return;
            }
            if ((atmin - atmid) == 0) min = mid;
            else max = mid;
        }
	if (its == max_it) {
            std::cerr << "bisect: overflow min " << min << " max " << max
                      << " diff " << max - min << " nroot " << nroot << " n1 "
                      << n1 << " n2 " << n2 << "\n";
            roots[0] = mid;
        }
        return;
    }

    // More than one root in the interval, we have to bisect...
    for (its = 0; its < max_it; its++) {
        mid = (min + max) / 2;
        atmid = num_changes(mid);
        n1 = atmin - atmid;
        n2 = atmid - atmax;
        if (n1 != 0 && n2 != 0) {
            bisect(min, mid, atmin, atmid, roots);
            bisect(mid, max, atmid, atmax, &roots[n1]);
            break;
        }

        if (n1 == 0) min = mid;
        else max = mid;
    }

    if (its == max_it) {
	std::cerr << "bisect: roots too close together\n";
        std::cerr << "bisect: overflow min " << min << " max " << max
                  << " diff " << max - min << " nroot " << nroot << " n1 " << n1
                  << " n2 " << n2 << "\n";
        for (n1 = atmax; n1 < atmin; n1++) roots[n1 - atmax] = mid;
    }
}

int Sturm::mod_rf(int ord, double *coef, double a, double b, double *val)
{
/*
  Uses the modified regula-falsi method to evaluate the root in interval
  [a,b] of the polynomial described in coef. The root is returned is
  returned in *val. The routine returns zero  if it can't converge.
 */
    int its;
    double fa, fb, x, fx, lfx;
    double *fp, *scoef, *ecoef;

    scoef = coef;
    ecoef = &coef[ord];

    fb = fa = *ecoef;
    for (fp = ecoef - 1; fp >= scoef; fp--) {
        fa = a * fa + *fp;
        fb = b * fb + *fp;
    }

    // if there is no sign difference the method won't work
    if (fa * fb > 0.0)
        return 0;
    if (fabs(fa) < rel_error) {
        *val = a;
        return 1;
    }
    if (fabs(fb) < rel_error) {
        *val = b;
        return 1;
    }

    lfx = fa;

    for (its = 0; its < max_it; its++) {
        x = (fb * a - fa * b) / (fb - fa);
        fx = *ecoef;
        for (fp = ecoef - 1; fp >= scoef; fp--)
            fx = x * fx + *fp;
        if (fabs(x) > rel_error) {
            if (fabs(fx / x) < rel_error) {
                *val = x;
                return 1;
            }
        }
        else if (fabs(fx) < rel_error) {
            *val = x;
            return 1;
        }
        if ((fa * fx) < 0) {
            b = x;
            fb = fx;
            if ((lfx * fx) > 0)
                fa /= 2;
        }
        else {
            a = x;
            fa = fx;
            if ((lfx * fx) > 0)
                fb /= 2;
        }

        lfx = fx;
    }

    std::cerr << "mod_rf overflow " << a << " " << b << " " << fx << "\n";

    return 0;
}

void Sturm::show_roots(const std::vector<double> &roots)
/*
  Write out the roots...
*/
{
    if (roots.size() == 1) {
        std::cout << "\n1 distinct real root at x = " << roots.front() << "\n";
    }
    else {
        std::cout << "\n" << roots.size() << " distinct real roots for x: ";
        for (const auto& root : roots) std::cout << root << ' ';
        std::cout << '\n';
    }
}

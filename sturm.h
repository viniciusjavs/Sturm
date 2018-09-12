#ifndef STURM_H
#define STURM_H

#include <iostream>
#include <cmath>
#include <vector>

constexpr int max_order = 12; // maximum order for a polynomial

class Sturm {
    
    // Defined type for representing a polynomial
    typedef struct Polynomial {
      public:
        int ord;
        double coef[max_order];
    } Poly;

    Poly sturm_seq[max_order];
    int order, num_poly, nroots, nchanges, atmin, atmax;
    double min = -1.0, max = 1.0, roots[max_order];
    
    int build_sturm();
    int modp(Poly *, Poly *, Poly *);
    int num_roots();
    int num_changes(const double &);
    void bisect(double, double, const int &, const int &, double *);
    double eval_poly(int, double *, double);
    int mod_rf(int, double *, double, double, double *);

  public:
    void set_parameters();
    std::vector<std::vector<double>> get_sturm_sequence();
    void show_sturm_sequence(const std::vector<std::vector<double>> &);
    std::vector<double> get_real_roots();
    void show_roots(const std::vector<double> &);
};

#endif


//
// Numerical Integration (Adaptive Gauss--Lobatto formula)
//
// Description:
//   Gauss--Lobatto formula is a numerical integrator
//   that is exact for polynomials of degree <= 2n+1.
//   Adaptive Gauss--Lobatto recursively decomposes the
//   domain and computes integral by using G-L formula.
//
// Algorithm:
//   Above.
//
// Complexity:
//   O(#pieces) for a piecewise polynomials.
//   In general, it converges in O(1/n^6) for smooth functions.
//   For (possibly) non-smooth functions, this is the best integrator.
//
// Verified: 
//   AOJ 2034
//
// References:
//   W. Gander and W. Gautschi (2000):
//   Adaptive quadrature - revisited.
//   BIT Numerical Mathematics, vol.40, no.1, pp.84--101.
//
#include <iostream>
#include <vector>
#include <limits>
#include <cstdio>
#include <algorithm>
#include <functional>

using namespace std;

#define fst first
#define snd second
#define all(c) ((c).begin()), ((c).end())

template <class F>
double integrate(F f, double lo, double hi, double eps = 1e-8) {
  const double th = eps / 1e-14; // (= eps / machine_epsilon)
  function<double (double,double,double,double,int)> rec =
    [&](double x0, double x6, double y0, double y6, int d) {
    const double a = sqrt(2.0/3.0), b = 1.0 / sqrt(5.0);
    double x3 = (x0 + x6)/2, y3 = f(x3), h = (x6 - x0)/2;
    double x1 = x3-a*h, x2 = x3-b*h, x4 = x3+b*h, x5 = x3+a*h;
    double y1 = f(x1),  y2 = f(x2),  y4 = f(x4),  y5 = f(x5);
    double I1 = (y0+y6 + 5*(y2+y4)) * (h/6);
    double I2 = (77*(y0+y6) + 432*(y1+y5) + 625*(y2+y4) + 672*y3) * (h/1470);
    if (x3 + h == x3 || d > 50) return 0.0; 
    if (d > 4 && th + (I1-I2) == th) return I2; // avoid degeneracy
    return (double)(rec(x0, x1, y0, y1, d+1) + rec(x1, x2, y1, y2, d+1) 
                  + rec(x2, x3, y2, y3, d+1) + rec(x3, x4, y3, y4, d+1) 
                  + rec(x4, x5, y4, y5, d+1) + rec(x5, x6, y5, y6, d+1));
  };
  return rec(lo, hi, f(lo), f(hi), 0);
}

int main() {
  for (int n; scanf("%d", &n); ) {
    if (n == 0) break;
    double r; 
    scanf("%lf\n", &r);

    vector<pair<double,double>> p(n);
    for (int i = 0; i < n; ++i) 
      scanf("%lf %lf", &p[i].fst, &p[i].snd);
    function<double (double)> f = [&](double x) {
      if (x < p[0].fst) return 0.0;
      if (x >= p.back().fst) return 0.0;
      for (int i = 0; i+1 < n; ++i) 
        if (x < p[i+1].fst) 
          return p[i].snd + (p[i+1].snd-p[i].snd)/(p[i+1].fst-p[i].fst)*(x-p[i].fst);
    };
    function<double (double)> g = [&](double x) {
      return f(x) * f(x + r);
    };
    printf("%.12f\n", integrate(g, p[0].fst-100, p.back().fst+100));
  }
}


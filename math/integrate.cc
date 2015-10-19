//
// Numerical Integration (double exponential integration)
//
// Description:
//   It computes int_a^b f(x) dx by using transformation
//     x = tanh(pi/2 sinh(t))
//   with the trapezoital rule. 
//   For analytic functions, this is the most accurate
//   formula in theory.
//
// Complexity:
//   Accuracy is exp(-O(N/log N)) with N function evaluations.
//
// Remark:
//   In general, if f is smooth upto k-th derivative,
//   O(1/N^(k+1)) is best possible.
//   If you know non-smooth points of f, you have to
//   split f by these points.
//
//
#include <iostream>
#include <vector>
#include <cmath>
#include <cstdio>
#include <algorithm>
#include <functional>

using namespace std;

#define fst first
#define snd second
#define all(c) ((c).begin()), ((c).end())

// doubly exponential integration; int_a^b f(x) dx; a = -infty or b = infty is allowed.
template <class F>
double integrate(F f, double a, double b, double eps = 1e-8) { 
  const double C = asin(1.0), T = 10.0; 
  function<double (double)> x, g; // g = dx/dt
  if (!isinf(a) && !isinf(b)) { // [a, b]
    x = [&](double t) { return (b+a)/2 + (b-a)/2*tanh(C*sinh(t)); };
    g = [&](double t) { return C*(b-a)/2*cosh(t)/pow(cosh(C*sinh(t)), 2); };
  } else if (!isinf(a)) {      // [a, infty]
    x = [&](double t) { return a + exp(C*sinh(t)); };
    g = [&](double t) { return C*cosh(t)*exp(C*sinh(t)); };
  } else if (!isinf(b)) {      // [-infty, b]
    x = [&](double t) { return b - exp(C*sinh(-t)); };
    g = [&](double t) { return C*cosh(t)*exp(C*sinh(-t)); };
  } else {                     // [-infty, infty]
    x = [&](double t) { return sinh(C*sinh(t)); };
    g = [&](double t) { return C*cosh(t)*cosh(C*sinh(t)); };
  }
  for (double S = 0, I, h = 1.0/C; ; S = I, h /= 2) {
    I = 0;
    for (double t = 0; t < T; t > 0 ? t = -t : t = -t+h) {
      double dI = f(x(t)) * g(t) * h;
      //printf("[%lf %lf %lf %lf %lf\n", t, I, dI, x(t), g(t));
      if (!isnan(dI)) I += dI;
    }
    //printf("{{{ %f }}}\n", I);
    if (fabs(I - S) < fabs(I) * eps) return I;
  }
}



// for comparison: Simpson integration
template <class F>
double integrate_S(F f, double a, double b, double eps = 1e-8) {
  auto simpson = [&](double w, double a, double b, double c) {
    return (b + a + 4 * c) * w / 6;
  };
  function<double (double,double,double,double)> rec = [&](double a, double b, double eps, double A) {
    double c = (a + b) / 2;
    double fa = f(a), fb = f(b), fc = f(c);
    double L = simpson(c-a, fa, fb, f((c+a)/2));
    double R = simpson(b-c, fc, fb, f((b+c)/2));
    if (fabs(L+R-A) <= 15*eps) return L+R+(L+R-A)/15;
    return rec(a, c, eps/2, L) + rec(c, b, eps/2, R);
  };
  return rec(a, b, eps, simpson(b - a, f(a), f(b), f((b + a) / 2)));
}


int main() {
  //auto f = [&](double x) { return exp(-x); };
  if (1) {
    int eval = 0;
    auto f = [&](double x) { ++eval; return exp(x); };
    double a = -10.0, b = 10.0;
    printf("EXACT = %.12lf\n", exp(b) - exp(a));
    printf("DE    = %.12lf; ", integrate(f, a, b, 1e-8));
    printf("%d function calls\n", eval);
    /*
    eval = 0;
    printf("SIMP  = %.12lf; ", integrate_S(f, a, b, 1e-8));
    printf("%d function calls\n", eval);
    */
  }
  if (1) {
    int eval = 0;
    auto f = [&](double x) { ++eval; return exp(-x); };
    double a = 0.0, b = 1.0/0.0;
    printf("EXACT = %.12lf\n", 1.0);
    printf("DE    = %.12lf; ", integrate(f, a, b, 1e-8));
    printf("%d function calls\n", eval);
  }
  if (1) {
    int eval = 0;
    auto f = [&](double x) { ++eval; return exp(x); };
    double a = -1.0/0.0, b = 0;
    printf("EXACT = %.12lf\n", 1.0);
    printf("DE    = %.12lf; ", integrate(f, a, b, 1e-8));
    printf("%d function calls\n", eval);
  }
  if (1) {
    int eval = 0;
    auto f = [&](double x) { ++eval; return exp(-x*x); };
    double a = -1.0/0.0, b = 1.0/0.0;
    printf("EXACT = %.12lf\n", sqrt(4.0*atan(1.0)));
    printf("DE    = %.12lf; ", integrate(f, a, b, 1e-8));
    printf("%d function calls\n", eval);
  }
}

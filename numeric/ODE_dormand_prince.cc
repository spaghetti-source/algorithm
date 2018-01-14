// 
// Dormand-Prince ODE solver
//
// Description:
//   It numerically solves an ordinary differential equation
//     dx/dt = f(t, x)
//
// Algorithm:
//   Dormand and Prince's adaptive step-size Runge-Kutta,
//   which is used in Matlab as "ode45" function.
//   This performs 4th and 5th-order Runge-Kutta methods
//   to estimate the optimal step size.
//
// References:
//   J. R. Dormand and P. J. Prince (1980):
//   A family of embedded Runge-Kutta formulae.
//   Journal of Computational and Applied Mathematics, vol.6, no.1, pp.13--26.
//
#include <iostream>
#include <vector>
#include <cstdio>
#include <algorithm>
#include <cmath>
#include <functional>

using namespace std;

#define fst first
#define snd second
#define all(c) ((c).begin()), ((c).end())
#define TEST(s) if (!(s)) { cout << __LINE__ << " " << #s << endl; exit(-1); }

template <class F>
double dormand_prince(F f, double t, double tend, double x) {
  double c[] = {0, 1./5, 3./10, 4./5, 8./9, 1, 1};
  double a[7][7] = {
    {0},
    {1./5},
    {3./40, 9./40},
    {44./45, -56./15, 32./9},
    {19372./6561, -25360./2187, 64448./6561, -212./729},
    {9017./3168, -355./33, 46732./5247, 49./176, -5103./18656},
    {35./384, 0, 500./1113, 125./192, -2187./6784, 11./84}
  };
  double b[] = {5179./57600, 0, 7571./16695, 393./640, -92097./339200, 187./2100, 1./40};
  double e[] = {71./57600, 0, -71./16695, 71./1920, -17253./339200, 22./525, -1./40};
  
  const double EPS = 1e-5;
  double h = EPS;
  while (t < tend) {
    if (t + h >= tend) h = tend - t;
    double k[7];
    for (int i = 0; i < 7; ++i) {
      double u = 0;
      for (int j = 0; j < i; ++j) 
        u += a[i][j] * k[j];
      k[i] = h * f(t + c[i] * h, x + u);
    }
    double err = 0;
    for (int i = 0; i < 7; ++i) x   += b[i] * k[i];
    t += h; // (t, x)

    for (int i = 0; i < 7; ++i) err += e[i] * k[i];
    double s = pow(EPS * h / (2 * abs(err)), 1./5);
    s = min(max(s, 1./4), 4.);
    h = s * h;
  }
  return x;
}

// for comparison
template <class F>
double runge_kutta(F f, double t, double tend, double x) {
  const double EPS = 1e-5;
  for (double h = EPS; t < tend; ) {
    if (t + h >= tend) h = tend - t;
    double k1 = h * f(t      , x       );
    double k2 = h * f(t + h/2, x + k1/2);
    double k3 = h * f(t + h/2, x + k2/2);
    double k4 = h * f(t + h  , x + k3  );
    x += (k1 + 2 * k2 + 2 * k3 + k4) / 6;
    t += h; // (t, x)
  }
  return x;
}

// for comparison
template <class F>
double euler(F f, double t, double tend, double x) {
  const double EPS = 1e-5;
  for (double h = EPS; t < tend; ) {
    if (t + h >= tend) h = tend - t;
    x += h * f(t, x);
    t += h;
  }
  return x;
}

int main() {
  auto f = [](double t, double x) {
    return t * x;
  };
  printf("%f\n", dormand_prince(f, 0, 1, 1));
  printf("%f\n", runge_kutta(f, 0, 1, 1));
  printf("%f\n", euler(f, 0, 1, 1));
  printf("%f\n", exp(1.0/2.0));
}

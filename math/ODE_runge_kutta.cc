// 
// The 4th order Runge-Kutta ODE solver
//
// Description:
//   It numerically solves an ordinary differential equation
//     dx/dt = f(t, x)
//
// Algorithm:
//   The 4th order Runge-Kutta algorithm. Let
//     k1 = h f(t,     x)
//     k2 = h f(t+h/2, x+k1/2)
//     k3 = h f(t+h/2, x+k2/2)
//     k4 = h f(t+h,   x+k3)
//   Then
//     x(t+h) = x(t) + (k1 + 2 k2 + 2 k3 + k4) / 6.
//
//   This is the most commonly used ODE solver, 
//   and referred as *THE* Runge Kutta method or RK4.
//   RK4 with sufficiently small step size gives usually 
//   sufficient result in solving non-stiff ODEs.
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
  printf("%f\n", runge_kutta(f, 0, 1, 1));
  printf("%f\n", euler(f, 0, 1, 1));
  printf("%f\n", exp(1.0/2.0));
}


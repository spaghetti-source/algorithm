// 
// Numerical Derivative by Ridder's method.
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
#define TEST(s) if (!(s)) { cout << __LINE__ << " " << #s << endl; exit(-1); }

template <class F>
double differentiate(F f, double x, double eps = 1e-8) {
  const int n = 10;
  const double alpha = 1.4;
  double h = 1e-2, a[n][n], ans = 1.0/0.0, err = 1.0/0.0;

  a[0][0] = (f(x + h) - f(x - h)) / (2 * h);
  for (int i = 1; i < n; ++i) {
    h /= alpha;
    a[0][i] = (f(x + h) - f(x - h))/(2 * h);
    double fac = alpha * alpha;
    for (int j = 1; j <= i; ++j) {
      a[j][i] = (a[j-1][i] * fac - a[j-1][i-1])/(fac - 1.0);
      fac *= alpha * alpha;
      double errt = max(fabs(a[j][i] - a[j-1][i]), fabs(a[j][i] - a[j-1][i-1]));
      if (errt <= err) {
        err = errt;
        ans = a[j][i];
        if (err < eps) return ans;
      }
    }
    if (fabs(a[i][i] - a[i-1][i-1]) >= 2 * err) break;
  }
  return ans;
}

double f(double x) {
  return exp(-x*x);
}
double df(double x) {
  return -2 * x * exp(-x*x);
}

int main() {
  for (int i = 0; i < 10; ++i) {
    double x = rand() / (1.0 + RAND_MAX);
    cout << differentiate(f, x) - df(x) << endl;
  }
}


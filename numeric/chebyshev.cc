//
// Chebyshev approximation of smooth function
//
// Description:
//
//   Given a smooth function f: [0,1] -> R. It approximates f by
//     f(x) \sim p(x) := sum_{k=0}^{n-1} ck Tk(x)
//   where Tk(x) is the Chebyshev polynomial of the first kind:
//     Tk(cos(u)) = cos(ku).
//   This almost minimizes the maximum error:
//     sup_x | f(x) - p(x) |
//   This also provides accurate dp/dx and \int_{a^x} p(u) du.
//
//   We can construct the Chebyshev approximation in O(n^2) time
//   by using the orthogonality with Clenshow's relation.
//
// Complexity:
//
//   O(n^2), where n is usually 10--30.
//
// References:
//
//   C. W. Clenshaw (1995): 
//   "A note on the summation of Chebyshev series." 
//   Mathematical Tables and Other Aids to Computation.
//   vol.9, no.51, pp.118--120.
//

#include <iostream>
#include <vector>
#include <cstdio>
#include <cmath>
#include <algorithm>
#include <functional>

using namespace std;

#define fst first
#define snd second
#define all(c) ((c).begin()), ((c).end())
#define TEST(s) if (!(s)) { cout << __LINE__ << " " << #s << endl; exit(-1); }

using Real = long double;
struct Chebyshev {
  static constexpr Real PI = acos(-1.0);
  int n;
  Real a, b;
  vector<Real> c;
  Chebyshev(Real a, Real b, int n) : n(n), a(a), b(b), c(n) { }

  template <class F>
  Chebyshev(F f, Real a, Real b, int n = 20) : n(n), a(a), b(b), c(n) {
    vector<Real> h(n);
    for (int k = 0; k < n; ++k) {
      Real y = cos(PI*(k+0.5)/n);
      h[k] = f((b-a)/2*y + (b+a)/2);
    }
    for (int j = 0; j < n; ++j) {
      for (int k = 0; k < n; ++k)
        c[j] += h[k] * cos(PI*j*(k+0.5)/n);
      c[j] *= 2.0/n;
    }
  }
  Real operator()(Real x) const {
    Real y = (2*x - a-b)/(b-a), u = 0, v = 0;
    for (int j = n-1; j >= 1; --j) {
      Real w = 2*y*u - v + c[j];
      v = u; u = w;
    }
    return y*u - v + 0.5*c[0];
  }
};
Chebyshev differentiate(Chebyshev f) {
  Chebyshev g = f;
  g.c[f.n-2] = 2 * (f.n-1) * f.c[f.n-1];
  for (int j = f.n-3; j >= 0; --j) 
    g.c[j] = g.c[j+2] + 2 * (j+1) * f.c[j+1];
  for (int j = 0; j < g.n; ++j)
    g.c[j] *= 2.0/(g.b - g.a);
  return g;
}
Chebyshev integrate(Chebyshev f) {
  Chebyshev g = f;
  Real sum = 0, coef = (f.b-f.a)/4, sign = 1.0;
  for (int j = 1; j <= g.n-2; ++j) {
    g.c[j] = coef * (f.c[j-1] - f.c[j+1]) / j;
    sum += sign * g.c[j];
    sign = -sign;
  }
  g.c[f.n-1] = coef * f.c[f.n-2] / (f.n-1);
  sum += sign * g.c[f.n-1];
  g.c[0] = 2 * sum;
  return g;
}

int main() {
  auto f = [&](Real x) {
    return sqrt(1 + x*x);
  };
  auto df = [&](Real x) {
    return x / sqrt(1 + x*x);
  };
  auto F = [&](Real x) {
    return (x*sqrt(1 + x*x) + asinh(x))/2;
  };
  Chebyshev g(f, 0, 1, 10);
  Chebyshev dg = differentiate(g);
  Chebyshev G = integrate(g);

  cout << "---" << endl;
  for (int i = 0; i < 10; ++i) {
    Real u = i / 10.0;
    cout << f(u) - g(u) << endl;
  }
  cout << "---" << endl;
  for (int i = 0; i < 10; ++i) {
    Real u = i / 10.0;
    cout << df(u) - dg(u) << endl;
  }
  cout << "---" << endl;
  for (int i = 0; i < 10; ++i) {
    Real u = i / 10.0;
    cout << F(u) - G(u) << endl;
  }
}

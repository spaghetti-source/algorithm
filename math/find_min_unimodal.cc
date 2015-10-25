//
// Minimum of unimodal function (goldsection search)
//
// Description:
//   Unimodal function is a function that has unique peak.
//
#include <iostream>
#include <vector>
#include <cstdio>
#include <algorithm>
#include <functional>

using namespace std;

#define fst first
#define snd second
#define all(c) ((c).begin()), ((c).end())
#define TEST(s) if (!(s)) { cout << __LINE__ << " " << #s << endl; exit(-1); }

template <class F>
double find_min(F f, double a, double d, double eps = 1e-8) { 
  const double r = 2 / (3 + sqrt(5.));
  double b = a + r*(d-a), c = d - r*(d-a), fb = f(b), fc = f(c);
  while (d - a > eps) {
    if (fb > fc) { // '<': maximum, '>': minimum
      a = b; b = c; c = d - r * (d - a);
      fb = fc; fc = f(c);
    } else {
      d = c; c = b; b = a + r * (d - a);
      fc = fb; fb = f(b);
    }
  }
  return c;
}

int main() {
  cout << find_min([](double x) { return x*x; }, -1, 1) << endl;
}

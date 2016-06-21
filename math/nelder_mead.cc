//
// Nelder Mead method (aka. Downhill Simplex Method)
//
// Description:
//   Nelder Mead method is a first-order optimization method
//   that only requires function evaluation oracle.
//   Typically, it performs well on a function on a small 
//   dimensional space.
//
// Algorithm:
//   Consider a simplex on the domain.
//   Evaluate function values on each vertex,
//   and replace the vertices according to some rules.
//   See the reference.
//   
// Reference:
//   Fuchang Gao and Lixing Han (2010):
//   Implementing the Nelder-Mead simplex algorithm with adaptive parameters.
//   Computational Optimization and Applications, vol.51, no.1, pp.259--277.
//
#include <iostream>
#include <vector>
#include <cstdio>
#include <algorithm>
#include <functional>
#include <cmath>

using namespace std;

#define fst first
#define snd second
#define all(c) ((c).begin()), ((c).end())
#define TEST(s) if (!(s)) { cout << __LINE__ << " " << #s << endl; exit(-1); }


#include <vector>
using namespace std;

template <class T>
vector<T> operator+(vector<T> x, vector<T> y) {
  for (int i = 0; i < y.size(); ++i) x[i] += y[i];
  return x;
}
template <class T>
vector<T> operator-(vector<T> x, vector<T> y) {
  for (int i = 0; i < y.size(); ++i) x[i] -= y[i];
  return x;
}
template <class T>
vector<T> operator*(T a, vector<T> x) {
  for (int i = 0; i < x.size(); ++i) x[i] *= a;
  return x;
}


template <class T>
ostream &operator<<(ostream &os, const vector<T> &v) {
  os << "[";
  for (int i = 0; i < v.size(); os << v[i++]) 
    if (i > 0) os << " ";
  os << "]";
  return os;
}
template <class T>
ostream &operator<<(ostream &os, const vector<vector<T>> &v) {
  os << "[";
  for (int i = 0; i < v.size(); os << v[i++]) 
    if (i > 0) os << endl << " ";
  os << "]";
  return os;
}

const double EPS = 1e-8;
template <class F>
pair<double, vector<double>> minimize(F f, vector<double> x0, double R = 100) {
  const int n = x0.size();
  const double alpha = 1.0, beta = 1.0 + 2.0/n, 
               gamma = 0.75 - 1.0/(2.0 * n), delta = 1.0 - 1.0/n;

  vector<vector<double>> xs = {x0};
  vector<double> fs = {f(x0)};
  vector<int> ord(n+1);
  for (int i = 1; i <= n; ++i) {
    xs.push_back(x0);
    xs[i][i-1] += (x0[i-1] == 0.0 ? 0.00025 : 0.05);
    fs.push_back(f(xs[i]));
    ord[i] = i;
  }
  sort(all(ord), [&](int i, int j) { return fs[i] < fs[j]; });

  auto centroid = [&]() {
    vector<double> x(n);
    for (int i = 0; i <= n; ++i) 
      if (i != ord[n]) x = x + xs[i];
    return (1.0/n)*x;
  };
  vector<double> c = centroid();

  for (int iter = 0; iter < 10000; ++iter) {
    auto replace = [&](vector<double> x, double v) {
      xs[ord[n]] = x; fs[ord[n]] = v;
      int j = n, k = ord[n];
      for (; j > 0 && fs[ord[j-1]] > fs[k]; --j) 
        ord[j] = ord[j-1];
      ord[j] = k;
      c = c + (1.0 / n) * (x - xs[ord[n]]);
    };
    auto shrink = [&]() {
      for (int i = 1; i <= n; ++i) {
        xs[i] = xs[0] + delta * (xs[i] - xs[0]);
        fs[i] = f(xs[i]);
      }
      sort(all(ord), [&](int i, int j) { return fs[i] < fs[j]; });
      c = centroid();
    };
    auto run = [&]() {
      auto xr = c + alpha * (c - xs[ord[n]]);
      auto fr = f(xr);
      if (fr < fs[ord[0]]) {
        auto xe = c + beta * (xr - c);
        auto fe = f(xe);
        if (fe < fr) return replace(xe, fe);
        else         return replace(xr, fr);
      }
      if (fr <= fs[ord[n-1]]) return replace(xr, fr);
      auto xc = c + (fr > fs[ord[n]] ? -gamma : gamma) * (xr - c);
      auto fc = f(xc);
      if (fc <= fr) return replace(xc, fc);
      else          return shrink();
    }; run();
    auto conv = [&]() { // conversion check
      for (int i = 0; i <= n; ++i) {
        if (fs[i] >= fs[ord[0]] + EPS) return false;
        for (int j = 0; j < n; ++j) 
          if (fabs(xs[i][j] - xs[ord[0]][j]) >= EPS) return false;
      }
      return true;
    }; if (conv()) break;
  }
  return {fs[ord[0]], xs[ord[0]]};
}

double rosenbrock(vector<double> x){ 
  return (1.0 - x[0])*(1.0 - x[0])
    + 100 * (x[1] - x[0]*x[0]) *  (x[1] - x[0]*x[0]);
}
double beale(vector<double> x) {
  return pow(1.5-x[0]+x[0]*x[1],2) + pow(2.25-x[0]+x[0]*x[1]*x[1],2) + pow(2.625-x[0]+x[0]*x[1]*x[1]*x[1], 2);
}
double bukin(vector<double> x) {
  return 100 * sqrt(fabs(x[1] - 0.01 * x[0]*x[0])) + 0.01 * fabs(x[0] + 10);
}
double egg_holder(vector<double> x) {
  return -(x[1]+47)*sin(sqrt(fabs(x[0]/2+(x[1]+47)))) - x[0]*sin(sqrt(fabs(x[0] - (x[1]+47))));
}
double holder_table(vector<double> x) {
  return -fabs(sin(x[0])*cos(x[1])*exp(fabs(1 - sqrt(x[0]*x[0]+x[1]*x[1])/M_PI)));
}

int main() {
  vector<double> x = {0, 0};
  auto ans = minimize(rosenbrock, x);
  cout << ans.fst << " at " << ans.snd << endl;

  ans = minimize(beale, x);
  cout << ans.fst << " at " << ans.snd << endl;

  x = {-11, 2};
  ans = minimize(bukin, x);
  cout << ans.fst << " at " << ans.snd << endl;

  ans = minimize(egg_holder, x);
  cout << ans.fst << " at " << ans.snd << endl;

  x = {8, 9};
  ans = minimize(holder_table, x);
  cout << ans.fst << " at " << ans.snd << endl;
}

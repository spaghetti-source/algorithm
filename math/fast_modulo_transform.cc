//
// Fast Modulo Transform and Fast Convolution in any Modulo
//
// Description:
//   For a given array x of size n (= 2^k), it computes
//     \hat x[i] = \sum_j x[i] h^j mod M
//   in O(n log n) time, where h^n == 1 (mod M).
//   For some nice M, we can efficiently compute h.
//
//   By using fmt, we can compute a convolution in arbitrary modulo
//   by using Chinese Remainder Theorem.
//
// Algorithm:
//   Fast modulo transform, which is a Z/Z_p version of fast fourier transform.
//
// Note:
//   We assume n is a power of 2 and n < 2^23 (>= 8*10^6)
//

#include <iostream>
#include <vector>
#include <cstdio>
#include <algorithm>
#include <functional>
#include <cassert>

using namespace std;

#define fst first
#define snd second
#define all(c) ((c).begin()), ((c).end())

typedef long long ll;
ll inv(ll b, ll M) {
  ll u = 1, x = 0, s = b, t = M; 
  while (s) {
    ll q = t / s;
    swap(x -= u * q, u);
    swap(t -= s * q, s);
  }
  return (x %= M) >= 0 ? x : x + M;
}
ll pow(ll a, ll b, ll M) {
  ll x = 1;
  for (; b > 0; b >>= 1) {
    if (b & 1) x = (a * x) % M;
    a = (a * a) % M;
  }
  return x;
}
// fast modulo transform
//   (1) n = 2^k < 2^23 
//   (2) only predetermined mod can be used
void fmt(vector<ll> &x, ll mod, int sign = +1) {
  int n = x.size();
  ll h = pow(3, (mod-1)/n, mod);
  if (sign < 0) h = inv(h, mod);
  for (int i = 0, j = 1; j < n-1; ++j) {
    for (int k = n >> 1; k > (i ^= k); k >>= 1);
    if (j < i) swap(x[i], x[j]);
  }
  for (int m = 1; m < n; m *= 2) {
    ll w = 1, wk = pow(h, n/(2*m), mod);
    for (int i = 0; i < m; ++i) {
      for (int j = i; j < n; j += 2*m) {
        ll u = x[j], d = x[j+m] * w % mod;
        if ((x[j]   = u + d) >= mod) x[j]   -= mod;
        if ((x[j+m] = u - d) <    0) x[j+m] += mod;
      }
      w = w * wk % mod;
    }
  }
  if (sign < 0) {
		ll n_inv = inv(n, mod);
		for (auto& a : x) a = (a * n_inv) % mod;
  }
}
// convolution via fast modulo transform
vector<ll> conv(vector<ll> x, vector<ll> y, ll mod) {
  fmt(x, mod, +1); fmt(y, mod, +1);
  for (int i = 0; i < x.size(); ++i) 
    x[i] = (x[i] * y[i]) % mod;
  fmt(x, mod, -1);
  return x;
}
// general convolution by using fmts with chinese remainder thm.
vector<ll> convolution(vector<ll> x, vector<ll> y, ll mod) {
  for (auto &a: x) a %= mod;
  for (auto &b: y) b %= mod;
  int n = x.size() + y.size() - 1, size = n - 1;
  for (int s: {1,2,4,8,16}) size |= (size >> s);
  size += 1;
  x.resize(size);
  y.resize(size);
  ll A = 167772161, B = 469762049, C = 1224736769, D = (A*B%mod);
  vector<ll> z(n), a = conv(x,y,A), b = conv(x,y,B), c = conv(x,y,C);
  for (int i = 0; i < n; ++i) {
    z[i] =  A * (104391568 * (b[i] - a[i]) % B);
    z[i] += D * (721017874 * (c[i] - (a[i] + z[i]) % C) % C);
    if ((z[i] = (z[i] + a[i]) % mod) < 0) z[i] += mod;
  }
  return z;
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

// naive algorithm for comparison
void fmt_n(vector<ll> &a, ll mod) {
  int n = a.size();
  vector<ll> b(n);
  ll g = pow(3, (mod-1)/n, mod), h = 1;
  for (int i = 0; i < n; ++i) {
    ll f = 1;
    for (int j = 0; j < n; ++j) {
      b[i] = (b[i] + (a[j] * f) % mod) % mod;
      f = (f * h) % mod;
    }
    h = (h * g) % mod;
  }
  swap(a, b);
}
vector<ll> convolution_n(vector<ll> a, vector<ll> b, ll mod) {
  vector<ll> c(a.size() + b.size() - 1);
  for (int i = 0; i < a.size(); ++i)
    for (int j = 0; j < b.size(); ++j)
      c[i+j] = (c[i+j] + a[i] * b[j]) % mod;
  return c;
}


// === tick a time ===
#include <ctime>
double tick() {
  static clock_t oldtick;
  clock_t newtick = clock();
  double diff = 1.0*(newtick - oldtick) / CLOCKS_PER_SEC;
  oldtick = newtick;
  return diff;
}

int main() {
  for (int iter = 2; iter < 100; ++iter) {
    cout << iter << endl;
    srand( iter );
    vector<ll> a(rand() + 10), b(rand() + 10);
    ll mod = 10000 + rand();
    for (auto &u: a) u = rand() % mod;
    for (auto &v: b) v = rand() % mod;
//    cout << a << endl;
//    cout << b << endl; 
    vector<ll> c = convolution(a, b, mod);
    vector<ll> d = convolution_n(a, b, mod);
//    cout << c << endl;
//    cout << d << endl;
    if (c != d) break;
  }
}

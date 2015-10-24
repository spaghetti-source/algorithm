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
// Verify:
//   SPOJ 235: Very Fast Multiplication (fmt is 3 times slower than fft)

#include <iostream>
#include <vector>
#include <cstdio>
#include <algorithm>
#include <cmath>
#include <complex>
#include <string>
#include <cstring>
#include <functional>
#include <cassert>

using namespace std;

#define fst first
#define snd second
#define all(c) ((c).begin()), ((c).end())

typedef long long ll;
ll mul(ll a, ll b, ll M) { // a * b (mod M)
  ll r = a*b - (ll)((long double)(a)*b/M+.5)*M;
  return r < 0 ? r + M: r;
}
ll div(ll a, ll b, ll M) { // solve b x == a (mod M)
  ll u = 1, x = 0, s = b, t = M; 
  while (s) { // extgcd for b x + M s = t
    ll q = t / s;
    swap(x -= u * q, u);
    swap(t -= s * q, s);
  }
  if (a % t) return -1; // infeasible
  return mul(x < 0 ? x + M : x, a / t, M); // b (xa/t) == a (mod M)
}
ll pow(ll a, ll b, ll M) {
  ll x = 1;
  for (; b > 0; b >>= 1) {
    if (b & 1) x = (a * x) % M;
    a = (a * a) % M;
  }
  return x;
}
// assume: size of a/b is power of two, mod is predetermined
template <int mod, int sign>
void fmt(vector<ll>& x) {
  const int n = x.size();
  int h = pow(3, (mod-1)/n, mod);
  if (sign < 0) h = div(1, h, mod);
  for (int i = 0, j = 1; j < n-1; ++j) {
    for (int k = n >> 1; k > (i ^= k); k >>= 1);
    if (j < i) swap(x[i], x[j]);
  }
  for (int m = 1; m < n; m *= 2) {
    ll w = 1, wk = pow(h, n / (2*m), mod);
    for (int i = 0; i < m; ++i) {
      for (int s = i; s < n; s += 2*m) {
        ll u = x[s], d = x[s + m] * w % mod;
        if ((x[s] = u + d) >= mod) x[s] -= mod;
        if ((x[s + m] = u - d) < 0) x[s + m] += mod;
      }
      w = w * wk % mod;
    }
  }
  if (sign < 0) { 
    ll inv = div(1, n, mod);
    for (auto &a: x) 
      a = a * inv % mod;
  }
}
// assume: size of a/b is power of two, mod is predetermined
template <int mod>
vector<ll> conv(vector<ll> a, vector<ll> b){
  fmt<mod,+1>(a); fmt<mod,+1>(b);
  for (int i = 0; i < a.size(); ++i) 
    a[i] = a[i] * b[i] % mod;
  fmt<mod,-1>(a);
  return a;
}
// general convolution where mod < 2^31.
vector<ll> conv(vector<ll> a, vector<ll> b, ll mod){
  int n = a.size() + b.size() - 1;
  for (int k: {1,2,4,8,16}) n |= (n >> k); ++n;
  a.resize(n); b.resize(n);
  const int A = 167772161, B = 469762049, C = 1224736769, D = (ll)(A) * B % mod;
  vector<ll> x = conv<A>(a,b), y = conv<B>(a,b), z = conv<C>(a,b);
  for (int i = 0; i < x.size(); ++i) {
    ll X = (y[i] - x[i]) * 104391568;
    if ((X %= B) < 0) X += B;
    ll Y = (z[i] - (x[i] + A * X) % C) * 721017874;
    if ((Y %= C) < 0) Y += C;
    x[i] += A * X + D * Y;
    if ((x[i] %= mod) < 0) x[i] += mod;
  }
  x.resize(n);
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
vector<ll> conv_n(vector<ll> a, vector<ll> b, ll mod) {
  vector<ll> c(a.size() + b.size() - 1);
  for (int i = 0; i < a.size(); ++i)
    for (int j = 0; j < b.size(); ++j)
      c[i+j] = (c[i+j] + a[i] * b[j]) % mod;
  return c;
}

const int WIDTH = 5;
const long long RADIX = 100000; // = 10^WIDTH

vector<ll> parse(const char s[]) {
  int n = strlen(s);
  int m = (n + WIDTH-1) / WIDTH;
  vector<ll> v(m);
  for (int i = 0; i < m; ++i) {
    int b = n - WIDTH * i, x = 0;
    for (int a = max(0, b - WIDTH); a < b; ++a)
      x = x * 10 + s[a] - '0';
    v[i] = x;
  }
  v.push_back(0);
  return v;
}

void print(const vector<ll> &v) {
  int i, N = v.size();
  vector<long long> digits(N + 1, 0);

  for (i = 0; i < N; i++) {
    digits[i] = v[i];
  }
  long long c = 0;
  for (i = 0; i < N; i++) {
    c += digits[i];
    digits[i] = c % RADIX;
    c /= RADIX;
  }
  for (i = N-1; i > 0 && digits[i] == 0; i--);
  printf("%lld", digits[i]);
  for (i--; i >= 0; i--)
    printf("%.*lld", WIDTH, digits[i]);
  printf("\n");
}

int main() {
  static char a_str[310000], b_str[310000];

  int T; scanf("%d", &T);
  while (T--) {
    scanf("%s %s", a_str, b_str);
    vector<ll> A = parse(a_str);
    vector<ll> B = parse(b_str);
    vector<ll> C = conv(A, B, 9000000000000000000ll);
    print(C);
  }
}

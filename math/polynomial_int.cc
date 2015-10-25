//
// Polynomial in Z/MZ
//
// Implemented routines:
//   1) addition 
//   2) subtraction
//   3) multiplication (naive, Karatsuba, FMT)
//   4) division (naive, Newton)
//   5) gcd
//
// TODO:
//   6) multipoint evaluation
//   7) interpolation
//   8) polynomial shift
//
#include <iostream>
#include <unordered_set> 
#include <vector>
#include <cassert>
#include <cstdio>
#include <algorithm>
#include <functional>

using namespace std;

#define fst first
#define snd second
#define all(c) ((c).begin()), ((c).end())
#define TEST(s) if (!(s)) { cout << __LINE__ << " " << #s << endl; exit(-1); }

typedef long long ll;
ll add(ll a, ll b, ll M) { // a + b (mod M)
  return (a += b) >= M ? a - M : a; 
}
ll sub(ll a, ll b, ll M) { // a - b (mod M)
  return (a -= b) < 0 ? a + M : a; 
}
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

// p(x) = p[0] + p[1] x + ... + p[n-1] x^n-1
// assertion: p.back() != 0
typedef vector<ll> poly;
ostream& operator<<(ostream &os, const poly &p) {
  const double EPS = 1e-4;
  bool head = true;
  for (int i = 0; i < p.size(); ++i) {
    if (p[i] == 0) continue;
    if (!head) os << " + ";
    os << p[i];
    head = false;
    if (i >= 1) os << " x";
    if (i >= 2) os << "^" << i;
  }
  return os;
}


// value of p(x)
ll eval(poly p, ll x, ll M) {
  ll ans = 0;
  for (int i = p.size()-1; i >= 0; --i)
    ans = add(mul(ans, x, M), p[i], M);
  return ans;
};

// q(x) = p(x + a)
poly shift(poly p, ll a, ll M) {
  poly q(p.size());
  for (int i = p.size()-1; i >= 0; --i) {
    for (int j = p.size()-i-1; j >= 1; --j) 
      q[j] = add(mul(q[j], a, M), q[j-1], M);
    q[0] = add(mul(q[0], a, M), p[i], M);
  }
  return q;
}

poly add(poly p, const poly &q, ll M) {
  if (p.size() < q.size()) p.resize(q.size());
  for (int i = 0; i < q.size(); ++i)
    p[i] = add(p[i], q[i], M);
  while (!p.empty() && !p.back()) p.pop_back();
  return p;
}
poly sub(poly p, const poly &q, ll M) {
  if (p.size() < q.size()) p.resize(q.size());
  for (int i = 0; i < q.size(); ++i)
    p[i] = sub(p[i], q[i], M);
  while (!p.empty() && !p.back()) p.pop_back();
  return p;
}

// naive multiplication in O(n^2)
poly mul_n(const poly &p, const poly &q, ll M) {
  if (p.empty() || q.empty()) return {};
  poly r(p.size() + q.size() - 1);
  for (int i = 0; i < p.size(); ++i)
    for (int j = 0; j < q.size(); ++j)
      r[i+j] = add(r[i+j], mul(p[i], q[j], M), M);
  while (!r.empty() && !r.back()) r.pop_back();
  return r;
}
// naive division (long division) in O(n^2)
pair<poly,poly> divmod_n(poly p, poly q, ll M) {
  poly u(p.size() - q.size() + 1); 
  ll inv = div(1, q.back(), M);
  for (int i = u.size()-1; i >= 0; --i) {
    u[i] = mul(p.back(), inv, M);
    for (int j = 0; j < q.size(); ++j) 
      p[j+p.size()-q.size()] = sub(p[j+p.size()-q.size()], mul(q[j], u[i], M), M);
    p.pop_back();
  }
  return {u, p};
}
// Karatsuba multiplication; this works correctly for M in [long long]
poly mul_k(poly p, poly q, ll M) {
  int n = max(p.size(), q.size()), m = p.size() + q.size() - 1;
  for (int k: {1,2,4,8,16}) n |= (n >> k); ++n; // n is power of two
  p.resize(n); q.resize(n);
  poly r(6*n);
  function<void(ll*, ll*, int, ll*)> rec = [&](ll *p0, ll *q0, int n, ll *r0) {
    if (n <= 4) { // 4 is the best threshold
      fill_n(r0, 2*n, 0);
      for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j)
          r0[i+j] = add(r0[i+j], mul(p0[i], q0[j], M), M);
      return;
    }
    ll *p1=p0+n/2,*q1=q0+n/2,*r1=r0+n/2,*r2=r0+n,*u=r0+5*n,*v=u+n/2,*w=r0+2*n;
    for (int i = 0; i < n/2; ++i) {
      u[i] = add(p0[i], p1[i], M);
      v[i] = add(q0[i], q1[i], M);
    }
    rec(p0, q0, n/2, r0);
    rec(p1, q1, n/2, r2);
    rec( u,  v, n/2, w);
    for (int i = 0; i < n; ++i) w[i] = sub(w[i], add(r0[i], r2[i], M), M);
    for (int i = 0; i < n; ++i) r1[i] = add(r1[i], w[i], M);
  }; rec(&p[0], &q[0], n, &r[0]);
  r.resize(m);
  return r;
}

// FFT-based multiplication: this works correctly for M in [int]
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
poly mul(poly p, poly q, ll M) { return conv(p, q, M); }

// Newton division: O(M(n) log n); M is the complexity of multiplication
// fast when FFT multiplication is used
pair<poly,poly> divmod(poly p, poly q, ll M) {
  poly s = p;
  reverse(all(p)); reverse(all(q));
  poly t = {div(1, q[0], M)}; 
  if (p.size() < q.size() || t[0] < 0) return { {}, {} };
  for (int k = 1; k <= 2*(p.size()-q.size()+1); k *= 2) {
    poly s = mul(mul(t, q, M), t, M);
    t.resize(k);
    for (int i = 0; i < k; ++i)
      t[i] = sub(2*t[i], s[i], M);
  }
  t.resize(p.size() - q.size() + 1);
  t = mul(t, p, M);
  t.resize(p.size() - q.size() + 1); 
  reverse(all(t)); reverse(all(p)); reverse(all(q));
  while (!t.empty() && !t.back()) t.pop_back();
  return {t, sub(p, mul(q, t, M), M) };
}
// polynomial GCD: O(D(n) log n); D is the complexity of division
poly gcd(poly p, poly q, ll M) {
  for (; !p.empty(); swap(p, q = divmod(q, p, M).snd));
  return p;
}

/* WIP
vector<ll> evaluate(poly p, vector<ll> x, ll M) {
  vector<ll> y(x.size());

  function<void(int,int,poly&,poly&)> rec = [&](int i, int j, poly p, poly &q) {
    if (i+1 == j) {
      y[i] = (p.empty() ? 0 : p[0]);
      q = {y[i], 1};
      return;
    }
    rec(i, (i+j)/2, divmod(p, pim).snd);
    rec((i+j)/2, j, divmod(p, pmj).snd);
  };
  return y;
}


// solve
// p(x) mod (x - x[i].fst) == x[i].snd for i = 0, ..., n-1
// ==> 
poly interpolate(vector<pair<ll,ll>> x, ll M) {
  // 
}
*/


// find p s.t. p(x[i].fst) = x[i].snd
// assert: all x[i] must be distinct
// O(n^2)
poly interpolate_n(vector<pair<ll,ll>> x, ll M) {
  int n = x.size();
  vector<ll> dp(n+1);
  dp[0] = 1;
  for (int i = 0; i < n; ++i) {
    for (int j = i; j >= 0; --j) {
      dp[j+1] = add(dp[j+1], dp[j], M);
      dp[j] = mul(dp[j], M - x[i].fst, M);
    }
  }
  poly r(n);
  for (int i = 0; i < n; ++i) {
    ll den = 1, res = 0;
    for (int j = 0; j < n; ++j) 
      if (i != j) den = mul(den, sub(x[i].fst, x[j].fst, M), M);
    den = div(1, den, M);
    
    for (int j = n-1; j >= 0; --j) {
      res = add(dp[j+1], mul(res, x[i].fst, M), M);
      r[j] = add(r[j], mul(res, mul(den, x[i].snd, M), M), M);
    }
  }
  while (!r.empty() && !r.back()) r.pop_back();
  return r;
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
bool TEST_DIVMOD() {
  for (int seed = 0; seed < 1; ++seed) {
    srand(seed);
    const ll M = 100000007;
    for (int n = 2; n < 1000; ++n) {
      int m = 1 + rand() % n;
      poly p(n), q(m);
      for (int i = 0; i < p.size(); ++i) p[i] = rand() % M;
      for (int i = 0; i < q.size(); ++i) q[i] = rand() % M;
      while (p.back() == 0) p.back() = rand() % M;
      while (q.back() == 0) q.back() = rand() % M;
      auto pq1 = divmod(p, q, M);
      TEST(sub(p, add(mul(pq1.fst, q, M), pq1.snd, M), M).empty()); 
      auto pq2 = divmod_n(p, q, M);
      TEST(sub(p, add(mul(pq2.fst, q, M), pq2.snd, M), M).empty()); 
      TEST(sub(pq1.fst, pq2.fst, M).empty() && sub(pq1.snd, pq2.snd, M).empty());
    }
  }
  return true;
}
bool TEST_MUL() {
  for (int seed = 0; seed < 1; ++seed) {
    srand(seed);
    const ll M = 100000007; 
    for (int n = 2; n < 1000; ++n) {
      poly p(n), q(n);
      for (int i = 0; i < p.size(); ++i) p[i] = rand() % M;
      for (int i = 0; i < q.size(); ++i) q[i] = rand() % M;
      while (p.back() == 0) p.back() = rand() % M;
      while (q.back() == 0) q.back() = rand() % M;
      auto pq1 = mul(p, q, M);
      auto pq2 = mul_n(p, q, M);
      TEST(sub(pq1, pq2, M).empty());
    }
  }
  return true;
}
bool TEST_INTERPOLATE() {
  for (int seed = 0; seed < 1; ++seed) {
    srand(seed);
    const ll M = 100000007; 
    for (int n = 2; n < 1000; ++n) {
      unordered_set<ll> xi;
      while (xi.size() < n) xi.insert(rand() % M);
      vector<pair<ll,ll>> x;
      for (auto &a: xi) x.push_back({a, rand() % M});
      poly q = interpolate_n(x, M);
      for (int i = 0; i < x.size(); ++i) 
        TEST(sub(eval(q, x[i].fst, M), x[i].snd, M) == 0);
    }
  }
  return true;
}

int main() {
  TEST(TEST_INTERPOLATE());
  TEST(TEST_DIVMOD());
  TEST(TEST_MUL());
}

//
// Polynomial with integer coefficient (mod M)
//
// Implemented routines:
//   1) addition 
//   2) subtraction
//   3) multiplication (naive O(n^2), Karatsuba O(n^1.5..), FFT O(n log n))
//   4) division (naive O(n^2), Newton O(M(n)))
//   5) gcd
//   6) multipoint evaluation (divide conquer: O(M(n) log |X|))
//   7) interpolation (naive O(n^2), divide conquer O(M(n) log n))
//   8) polynomial shift (naive, fast)
//
//   *) n! mod M in O(n^{1/2} log n) time
//
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
poly mul(poly p, poly q, ll M) { 
  poly pq = conv(p, q, M);
  pq.resize(p.size() + q.size() - 1);
  while (!pq.empty() && !pq.back()) pq.pop_back();
  return pq;
}

// Newton division: O(M(n)); M is the complexity of multiplication
// fast when FFT multiplication is used
//
// Note: complexity = M(n) + M(n/2) + M(n/4) + ... <= 2 M(n).
pair<poly,poly> divmod(poly p, poly q, ll M) {
  if (p.size() < q.size()) return { {}, p };
  reverse(all(p)); reverse(all(q));
  poly t = {div(1, q[0], M)}; 
  if (t[0] < 0) return { {}, {} }; // infeasible
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
// polynomial GCD: O(M(n) log n); 
poly gcd(poly p, poly q, ll M) {
  for (; !p.empty(); swap(p, q = divmod(q, p, M).snd));
  return p;
}

// value of p(x)
ll eval(poly p, ll x, ll M) {
  ll ans = 0;
  for (int i = p.size()-1; i >= 0; --i)
    ans = add(mul(ans, x, M), p[i], M);
  return ans;
};
//
// faster multipoint evaluation
// fast if |x| >= 10000.
//
// algo:
//   evaluate(p, {x[0], ..., x[n-1]})
//   = evaluate(p mod (X-x[0])...(X-x[n/2-1]), {x[0], ..., x[n/2-1]}),
//   + evaluate(p mod (X-x[n/2])...(X-x[n-1]), {x[n/2], ..., x[n-1]}),
//
// f(n) = 2 f(n/2) + M(n) ==> O(M(n) log n)
//
vector<ll> evaluate(poly p, vector<ll> x, ll M) {
  vector<poly> prod(8*x.size()); // segment tree
  function<poly(int,int,int)> run = [&](int i, int j, int k) {
    if (i   == j) return prod[k] = (poly){1};
    if (i+1 == j) return prod[k] = (poly){M-x[i], 1}; 
    return prod[k] = mul(run(i,(i+j)/2,2*k+1), run((i+j)/2,j,2*k+2), M);
  }; run(0, x.size(), 0);
  vector<ll> y(x.size());
  function<void(int,int,int,poly)> rec = [&](int i, int j, int k, poly p) {
    if (j - i <= 8) {
      for (; i < j; ++i) y[i] = eval(p, x[i], M);
    } else {
      rec(i, (i+j)/2, 2*k+1, divmod(p, prod[2*k+1], M).snd);
      rec((i+j)/2, j, 2*k+2, divmod(p, prod[2*k+2], M).snd);
    }
  }; rec(0, x.size(), 0, p);
  return y;
}

poly interpolate_n(vector<ll> x, vector<ll> y, ll M) {
  int n = x.size();
  vector<ll> dp(n+1);
  dp[0] = 1;
  for (int i = 0; i < n; ++i) {
    for (int j = i; j >= 0; --j) {
      dp[j+1] = add(dp[j+1], dp[j], M);
      dp[j] = mul(dp[j], M - x[i], M);
    }
  }
  poly r(n);
  for (int i = 0; i < n; ++i) {
    ll den = 1, res = 0;
    for (int j = 0; j < n; ++j) 
      if (i != j) den = mul(den, sub(x[i], x[j], M), M);
    den = div(1, den, M);
    
    for (int j = n-1; j >= 0; --j) {
      res = add(dp[j+1], mul(res, x[i], M), M);
      r[j] = add(r[j], mul(res, mul(den, y[i], M), M), M);
    }
  }
  while (!r.empty() && !r.back()) r.pop_back();
  return r;
}
// 
// faster algo to find a poly p such that
//   p(x[i]) = y[i]  for each i
//
// see http://people.mpi-inf.mpg.de/~csaha/lectures/lec6.pdf
//
poly interpolate(vector<ll> x, vector<ll> y, ll M) {
  vector<poly> prod(8*x.size()); // segment tree
  function<poly(int,int,int)> run = [&](int i, int j, int k) {
    if (i   == j) return prod[k] = (poly){1};
    if (i+1 == j) return prod[k] = (poly){M-x[i], 1}; 
    return prod[k] = mul(run(i,(i+j)/2,2*k+1), run((i+j)/2,j,2*k+2), M);
  }; run(0, x.size(), 0); // preprocessing in O(n log n) time

  poly H = prod[0]; // newton polynomial
  for (int i = 1; i < H.size(); ++i) H[i-1] = mul(H[i], i, M);
  do H.pop_back(); while (!H.empty() && !H.back());

  vector<ll> u(x.size());
  function<void(int,int,int,poly)> rec = [&](int i, int j, int k, poly p) {
    if (j - i <= 8) {
      for (; i < j; ++i) u[i] = eval(p, x[i], M);
    } else {
      rec(i, (i+j)/2, 2*k+1, divmod(p, prod[2*k+1], M).snd);
      rec((i+j)/2, j, 2*k+2, divmod(p, prod[2*k+2], M).snd);
    }
  }; rec(0, x.size(), 0, H); // multipoint evaluation

  for (int i = 0; i < x.size(); ++i) u[i] = div(y[i], u[i], M); 

  function<poly(int,int,int)> f = [&](int i, int j, int k) {
    if (i   >= j) return poly();
    if (i+1 == j) return (poly){u[i]};
    return add(mul(f(i,(i+j)/2,2*k+1), prod[2*k+2], M),
               mul(f((i+j)/2,j,2*k+2), prod[2*k+1], M), M);
  };
  return f(0, x.size(), 0);
}


//
// return p(x+a) 
// 
poly shift_n(poly p, ll a, ll M) {
  poly q(p.size());
  for (int i = p.size()-1; i >= 0; --i) {
    for (int j = p.size()-i-1; j >= 1; --j) 
      q[j] = add(mul(q[j], a, M), q[j-1], M);
    q[0] = add(mul(q[0], a, M), p[i], M);
  }
  return q;
}

//
// faster algorithm for computing p(x + a) 
//
// fast if n >= 4096
// algo: p(x+a) = p_h(x) (x+a)^m + q_h(x)
// cplx: preproc: O(M(n))
//       div-con: O(M(n) log n)
//
poly shift(poly p, ll a, ll M) {
  vector<poly> pow(p.size()); 
  pow[0] = {1}; pow[1] = {a,1};
  int m = 2; 
  for (; m < p.size(); m *= 2) 
    pow[m] = mul(pow[m/2], pow[m/2], M);
  function<poly(poly,int)> rec = [&](poly p, int m) {
    if (p.size() <= 1) return p;
    while (m >= p.size()) m /= 2;
    poly q(p.begin() + m, p.end());
    p.resize(m);
    return add(mul(rec(q, m), pow[m], M), rec(p, m), M);
  };
  return rec(p, m);
}


//
// overpeform when n >= 134217728 lol
// 
ll factmod(ll n, ll M) {
  if (n <= 1) return 1;
  ll m = sqrt(n); 
  function<poly(int,int)> get = [&](int i, int j) {
    if (i   == j) return poly();
    if (i+1 == j) return (poly){i,1};
    return mul(get(i, (i+j)/2), get((i+j)/2, j), M);
  };
  poly p = get(0, m); // = x (x+1) (x+2) ... (x+(m-1))
  vector<ll> x(m);
  for (int i = 0; i < m; ++i) x[i] = 1 + i * m;
  vector<ll> y = evaluate(p, x, M);
  ll fac = 1;
  for (int i = 0; i < m; ++i)
    fac = mul(fac, y[i], M);
  for (ll i = m*m+1; i <= n; ++i)
    fac = mul(fac, i, M);
  return fac;
}
ll factmod_n(ll n, ll M) {
  ll fac = 1;
  for (ll k = 1; k <= n; ++k)
    fac = mul(k, fac, M);
  return fac;
}
ll factmod_p(ll n, ll M) { // only works for prime M
  ll fac = 1;
  for (; n > 1; n /= M) {
    fac = mul(fac, (n / M) % 2 ? M - 1 : 1, M);
    for (ll i = 2; i <= n % M; ++i)
      fac = mul(fac, i, M);
  }
  return fac;
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
bool TEST_EVALUATE() {
  const ll M = 100000007;
  for (int seed = 0; seed < 1; ++seed) {
    srand(seed);
    for (int n = 2; n < 1000; ++n) {
      poly p(n);
      for (int i = 0; i < n; ++i) p[i] = rand() % M;
      while (p.back() == 0) p.back() = rand() % M;
      int m = 1 + rand() % (2 * n);
      vector<ll> x(m);
      for (int i = 0; i < m; ++i) x[i] = rand() % M;
      tick();
      vector<ll> y = evaluate(p, x, M);
      //cout << n << " " << tick() << " ";
      vector<ll> z(m);
      for (int i = 0; i < m; ++i) z[i] = eval(p, x[i], M);
      //cout << tick() << endl;
      for (int i = 0; i < x.size(); ++i) 
        TEST(sub(y[i], eval(p, x[i], M), M) == 0);
    }
  }
  return true;
}
bool TEST_INTERPOLATE() {
  for (int seed = 0; seed < 1; ++seed) {
    srand(seed);
    const ll M = 100000007; 
    for (int n = 2; n < 100000; n*=2) {
      unordered_set<ll> xi;
      while (xi.size() < n) xi.insert(rand() % M);
      vector<ll> x(all(xi)), y(x.size());
      for (int i = 0; i < y.size(); ++i) y[i] = rand() % M;
      cout << n << " ";
      poly q = interpolate(x, y, M);
      cout << tick() << " ";
      for (int i = 0; i < x.size(); ++i) 
        TEST(sub(eval(q, x[i], M), y[i], M) == 0);
      poly r = interpolate_n(x, y, M);
      cout << tick() << endl;
      for (int i = 0; i < x.size(); ++i) 
        TEST(sub(eval(r, x[i], M), y[i], M) == 0);
      TEST(sub(q, r, M).empty());
    }
  }
  return true;
}
bool TEST_SHIFT() {
  for (int seed = 0; seed < 1; ++seed) {
    srand(seed);
    const ll M = 1000000007; 
    for (int n = 2; n < 100000; n *= 2) {
      poly p(n);
      for (int i = 0; i < n; ++i) p[i] = rand() % M;
      while (p.back() == 0) p.back() = rand() % M;
      ll a = rand() % M;
      cout << n << " ";
      tick();
      poly q = shift(p, a, M);
      cout << tick() << " ";
      poly r = shift_n(p, a, M);
      cout << tick() << endl;
      TEST(sub(q, r, M).empty());
      ll b = rand() % M;
      TEST(sub(eval(p, add(a,b,M), M), eval(q, b, M), M) == 0);
    }
  }
}

int main() {

  TEST_SHIFT();

  const ll M = 1000000007;

  poly p = {0, 0, 0, 0, 0, 1};
  poly q = shift(p, 1, M);
  cout << q << endl;

  return 0;
  for (int i = 1; i < M; i *= 2) {
    cout << i << " ";
    cout << factmod_p(i, M) << " ";
    cout << tick() << "  ";
    cout << factmod_n(i, M) << " ";
    cout << tick() << endl;
  }
  return 0;
  TEST(TEST_INTERPOLATE());
  TEST(TEST_EVALUATE());
  TEST(TEST_DIVMOD());
  TEST(TEST_MUL());
  return 0;
}

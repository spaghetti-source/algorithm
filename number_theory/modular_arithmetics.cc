//
// Modular arithmetics (long long)
//
// Note:
//   int       < 2^31 < 10^9
//   long long < 2^63 < 10^18
//
// feasible for M < 2^62 (10^18 < 2^62 < 10^19)
//
//
// Verified:
//   SPOJ 11409: Fibonacci With a Twist
//   SPOJ 9832: Matrix Inverse
//
//

#include <cassert>
#include <cstring>
#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <vector>
#include <ctime>
#include <functional>
#include <algorithm>

using namespace std;


#define All(c) c.begin(), c.end()
#define FOR(i,c) for(typeof(c.begin())i=c.begin();i!=c.end();++i)
#define REP(i,n) for(int i=0;i<n;++i)
#define fst first
#define snd second

typedef long long ll;
typedef vector<ll> vec;
typedef vector<vec> mat;

ll add(ll a, ll b, ll M) {
  a += b;
  if (a >= M) a -= M;
  return a;
}
ll sub(ll a, ll b, ll M) { 
  if (a < b) a += M;
  return a - b;
}

// Correctness of mul
//   ab = floor(ab / M) * M + (ab % M)
//   -> (ab % M) = ab - floor(ab / M) * M
ll mul(ll a, ll b, ll M) {
  ll r = a*b - floor(1.0*a*b/M)*M;
  return r < 0 ? r + M : r >= M ? r - M : r;
}
ll pow(ll a, ll b, ll M) {
  ll x = 1;
  for (; b > 0; b >>= 1) {
    if (b & 1) x = mul(x, a, M);
    a = mul(a, a, M);
  }
  return x;
}
ll div(ll a, ll b, ll M) {
  ll u = 1, x = 0, s = b, t = M; 
  while (s) {
    ll q = t / s;
    swap(x -= u * q, u);
    swap(t -= s * q, s);
  }
  if (a % t) return -1; // infeasible
  return mul(x < 0 ? x + M : x, a / t, M);
}


// Modular Matrix
mat eye(int n) {
  mat I(n, vec(n));
  REP(i, n) I[i][i] = 1;
  return I;
}
mat zeros(int n) {
  return mat(n, vec(n));
}
mat mul(mat A, mat B, ll M) {
  int l = A.size(), m = B.size(),  n = B[0].size();
  mat C(l, vec(n));
  REP(i,l) REP(k,m) REP(j,n) 
    C[i][j] = add(C[i][j], mul(A[i][k], B[k][j], M), M);
  return C;
}
mat pow(mat A, ll b, ll M) {
  mat X = eye(A.size());
  for (; b > 0; b >>= 1) {
    if (b & 1) X = mul(X, A, M);
    A = mul(A, A, M);
  }
  return X;
}
// assume: M is prime (singular ==> 
// verify: SPOJ9832
mat inv(mat A, ll M) {
  int n = A.size();
  mat B(n, vec(n));
  for (int i = 0; i < n; ++i) 
    B[i][i] = 1;

  for (int i = 0; i < n; ++i) {
    int j = i;
    while (j < n && A[j][i] == 0) ++j;
    if (j == n) return {};
    swap(A[i], A[j]); 
    swap(B[i], B[j]);

    ll inv = div(1, A[i][i], M);
    for (int k = i; k < n; ++k)
      A[i][k] = mul(A[i][k], inv, M);
    for (int k = 0; k < n; ++k)
      B[i][k] = mul(B[i][k], inv, M);
    for (int j = 0; j < n; ++j) {
      if (i == j || A[j][i] == 0) continue;
      ll cor = A[j][i];
      for (int k = i; k < n; ++k)
        A[j][k] = sub(A[j][k], mul(cor, A[i][k], M), M);
      for (int k = 0; k < n; ++k)
        B[j][k] = sub(B[j][k], mul(cor, B[i][k], M), M);
    }
  }
  return B;
}



void disp(mat A) {
  cout << "[";
  REP(i, A.size()) {
    if (i != 0) cout << " ";
    REP(j, A[i].size()) {
      cout << A[i][j];
      if (j != A[i].size()-1) cout << ", ";
      else                    cout << "; ";
    }
    cout << endl;
  }
  cout << endl;
}


ll binomial(ll n, ll k, ll M) {
  ll num = 1, den = 1;
  while (n > 0 || k > 0) {
    ll m = n % M, l = k % M;
    if (m < l) return 0;
    if (l > m - l) l = m - l;
    while (l > 0) {
      num = mul(num, m--, M);
      den = mul(den, l--, M);
    }
    n /= M; k /= M;
  }
  return div(num, den, M);
}


// tick a time 
double tick() {
  static clock_t oldtick;
  clock_t newtick = clock();
  double diff = 1.0*(newtick - oldtick) / CLOCKS_PER_SEC;
  oldtick = newtick;
  return diff;
}


const int N = 10000000;
ll x[N];
int main() {
}

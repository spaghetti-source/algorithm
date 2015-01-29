// 
// Carmichael Lambda (Universal Totient Function)
//
// Description:
//   lambda(n) is a smallest number that satisfies
//     a^lambda(n) = 1 (mod n)
//   for all a coprime with n. 
//   This is also known as an universal totien tunction psi(n).
//
// Complexity:
//   carmichael_lambda(n):     O(sqrt(n)) by trial division.
//   carmichael_lambda(lo,hi): O((hi-lo) loglog(hi)) by prime sieve.
//
// Verified:
//   Upto 10000

#include <iostream>
#include <vector>
#include <cstdio>
#include <algorithm>
#include <unordered_map>
#include <functional>

using namespace std;

#define fst first
#define snd second
#define all(c) ((c).begin()), ((c).end())

typedef long long ll;


ll gcd(ll a, ll b) {
  for (; a; swap(b %= a, a));
  return b;
}
ll lcm(ll a, ll b) {
  return a * (b / gcd(a, b));
}

ll carmichael_lambda(ll n) {
  ll lambda = 1;
  if (n % 8 == 0) n /= 2;
  for (ll d = 2; d*d <= n; ++d) {
    if (n % d == 0) {
      n /= d;
      ll y = d - 1;
      while (n % d == 0) {
        n /= d;
        y *= d;
      }
      lambda = lcm(lambda, y);
    }
  }
  if (n > 1) lambda = lcm(lambda, n-1);
  return lambda;
}

vector<ll> primes(ll lo, ll hi) { // primes in [lo, hi)
  const ll M = 1 << 14, SQR = 1 << 16;
  vector<bool> composite(M), small_composite(SQR);

  vector<pair<ll,ll>> sieve; 
  for (ll i = 3; i < SQR; i+=2) {
    if (!small_composite[i]) {
      ll k = i*i + 2*i*max(0.0, ceil((lo - i*i)/(2.0*i)));
      sieve.push_back({2*i, k});
      for (ll j = i*i; j < SQR; j += 2*i) 
        small_composite[j] = 1;
    }
  }
  vector<ll> ps; 
  if (lo <= 2) { ps.push_back(2); lo = 3; }
  for (ll k = lo|1, low = lo; low < hi; low += M) {
    ll high = min(low + M, hi);
    fill(all(composite), 0);
    for (auto &z: sieve) 
      for (; z.snd < high; z.snd += z.fst)
        composite[z.snd - low] = 1;
    for (; k < high; k+=2) 
      if (!composite[k - low]) ps.push_back(k);
  }
  return ps;
}
vector<ll> primes(ll n) { // primes in [0,n)
  return primes(0,n);
}
vector<ll> carmichael_lambda(ll lo, ll hi) { // lambda(n) for all n in [lo, hi)
  vector<ll> ps = primes(sqrt(hi)+1);
  vector<ll> res(hi-lo), lambda(hi-lo, 1);
  iota(all(res), lo);

  for (ll k = ((lo+7)/8)*8; k < hi; k += 8) res[k-lo] /= 2;
  for (ll p: ps) {
    for (ll k = ((lo+(p-1))/p)*p; k < hi; k += p) {
      if (res[k-lo] < p) continue;
      ll t = p - 1;
      res[k-lo] /= p;
      while (res[k-lo] > 1 && res[k-lo] % p == 0) {
        t *= p;
        res[k-lo] /= p; 
      }
      lambda[k-lo] = lcm(lambda[k-lo], t);
    }
  }
  for (ll k = lo; k < hi; ++k) {
    if (res[k-lo] > 1) 
      lambda[k-lo] = lcm(lambda[k-lo], res[k-lo]-1);
  }
  return lambda; // lambda[k-lo] = lambda(k)
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
  int n = 1000000;
  auto lmbd = carmichael_lambda(0, n);
  for (int i = 0; i < n; ++i) {
    if (lmbd[i] != carmichael_lambda(i)) {
      cout << i << " " << lmbd[i] << " " << carmichael_lambda(i) << endl;
    }
  }

/*
  for (int iter = 0; iter < 100; ++iter) {
    int lo = rand(), hi = lo + rand();
    auto x = euler_phi(lo, hi);
    for (int n = lo; n < hi; ++n)
      if (x[n-lo] != euler_phi(n)) cout << "!!" << endl;
  }
  */
}

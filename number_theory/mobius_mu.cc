// 
// Mobius Mu 
//
// Description:
//     mu(n) =  1 if n is square-free, even number of prime factors
//             -1 if ...               odd  ...
//              0 if n has a squared prime factor.
//   equivalently, multiplicative with
//     mu(p^k) =  1 if k = 0
//               -1 if k = 1
//                0 if k > 1.
//
// Complexity:
//   mobius_mu(n):     O(sqrt(n)) by trial division.
//   mobius_mu(lo,hi): O((hi-lo) loglog(hi)) by prime sieve.
//
// Verified:
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

typedef long long ll;

ll mobius_mu(ll n) {
  if (n == 0) return 0;
  ll mu = 1;
  for (ll x = 2; x*x <= n; ++x) {
    if (n % x == 0) {
      mu = -mu;
      n /= x;
      if (n % x == 0) return 0;
    }
  }
  return n > 1 ? -mu : mu;
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
vector<ll> mobius_mu(ll lo, ll hi) { // phi(n) for all n in [lo, hi)
  vector<ll> ps = primes(sqrt(hi)+1);
  vector<ll> res(hi-lo), mu(hi-lo, 1);
  iota(all(res), lo);

  for (ll p: ps) {
    for (ll k = ((lo+(p-1))/p)*p; k < hi; k += p) {
      mu[k-lo] = -mu[k-lo];
      if (res[k-lo] % p == 0) {
        res[k-lo] /= p;
        if (res[k-lo] % p == 0) {
          mu[k-lo] = 0;
          res[k-lo] = 1;
        }
      }
    }
  }
  for (ll k = lo; k < hi; ++k) {
    if (res[k-lo] > 1) 
      mu[k-lo] = -mu[k-lo];
  }
  return mu; // mu[k-lo] = mu(k)
}

int main() {
  for (int iter = 0; iter < 1000; ++iter) {
    int lo = rand(), hi = lo + rand();
    auto x = mobius_mu(lo, hi);
    for (int i = lo; i < hi; ++i) {
      if (x[i-lo] != mobius_mu(i)) {
        cout << lo << " "  << hi << "  " << mobius_mu(i) << " " << x[i-lo] << endl;
      }
    }
  }
}

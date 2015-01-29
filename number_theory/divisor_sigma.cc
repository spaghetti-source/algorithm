//
// Divisor function
//
// Description:
//     sigma(n) = sum[n % d == 0] d
//   equivalently,
//     sigma(p^k) = 1 + p + p^2 + ... + p^k
//   with multiplicative.
//
// Complexity:
//   divisor_sigma(n):     O(sqrt(n)) by trial division.
//   divisor_sigma(lo,hi): O((hi-lo) loglog(hi)) by prime sieve.
//

#include <iostream>
#include <cstdio>
#include <vector>
#include <algorithm>
 
using namespace std;

#define fst first
#define snd second
#define all(c) ((c).begin()), ((c).end())
 
typedef long long ll;
ll divisor_sigma(ll n) {
  ll sigma = 0, d = 1;
  for (; d*d < n; ++d)
    if (n % d == 0) sigma += d + n/d;
  if (d*d == n) sigma += d;
  return sigma;
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
vector<ll> divisor_sigma(ll lo, ll hi) { // sigma(n) for all n in [lo, hi)
  vector<ll> ps = primes(sqrt(hi)+1);
  vector<ll> res(hi-lo), sigma(hi-lo, 1);
  iota(all(res), lo);

  for (ll p: ps) {
    for (ll k = ((lo+(p-1))/p)*p; k < hi; k += p) {
      ll b = 1;
      while (res[k-lo] > 1 && res[k-lo] % p == 0) {
        res[k-lo] /= p; 
        b = 1 + b * p;
      }
      sigma[k-lo] *= b;
    }
  }
  for (ll k = lo; k < hi; ++k) 
    if (res[k-lo] > 1) 
      sigma[k-lo] *= (1 + res[k-lo]);
  return sigma; // sigma[k-lo] = sigma(k)
}

int main() {
  for (int i = 0; i < 17; ++i)
    cout << divisor_sigma(i) << " ";
  cout << endl;

  auto x = divisor_sigma(0, 17);
  for (int i = 0; i < 17; ++i)
    cout << x[i] << " ";
  cout << endl;

  for (int iter = 0; iter < 100; ++iter) {
    int lo = rand(), hi = lo + rand();
    auto x = divisor_sigma(lo, hi);
    for (int n = lo; n < hi; ++n)
      if (x[n-lo] != divisor_sigma(n)) cout << "!!" << endl;
  }
}

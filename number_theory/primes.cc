// 
// Generate primes 
//
// Description:
//   Generating primes from 1 to n.
//
// Algorithm:
//   Segmented sieve. It first enumerates small primes (upto sqrt(n))
//   Then it sieve the rest numbers.
//
// Complexity: 
//   O(n log log n)
//
// Verified:
//   SPOJ2, SPOJ503
//
#include <iostream>
#include <vector>
#include <cstdio>
#include <cstdlib>
#include <map>
#include <cmath>
#include <cstring>
#include <functional>
#include <algorithm>

using namespace std;

#define fst first
#define snd second
#define all(c) c.begin(), c.end()

typedef long long ll;
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

// SPOJ503
int main() { 
  int n; scanf("%d", &n);
  for (int i = 0; i < n; ++i) {
    ll lo, hi;
    scanf("%lld %lld", &lo, &hi);
    auto x = primes(lo, hi+1);
    for (ll p: x) printf("%lld\n", p);
  }
}
/* 
// SPOJ2 
int main() { 
  int n; scanf("%d", &n);
  for (int i = 0; i < n; ++i) {
    if (i > 0) printf("\n");
    ll lo, hi;
    scanf("%lld %lld", &lo, &hi);
    auto x = primes(lo, hi+1);
    for (ll p: x) printf("%lld\n", p);
  }
}
*/

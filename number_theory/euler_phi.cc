// 
// Euler Phi (Totient Function)
//
// Description:
//     phi(n) = #{ k <= n : k is coprime to n }
//            = n (1 - 1/p1) ... (1 - 1/pm).
//   or equivalently
//     phi(p^k) = (p-1) p^{k-1}.
//   with multiplicative.
//
// Complexity:
//   euler_phi(n):     O(sqrt(n)) by trial division.
//   euler_phi(lo,hi): O((hi-lo) loglog(hi)) by prime sieve.
//
// Verified:
//   SPOJ 22268
//
// Note:
//   Complexity of sieve ver. equals to the sum of exponents in hi!/lo!.
//   This is known to be O(hi loglog hi - lo loglog lo).

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

ll euler_phi(ll n) {
  if (n == 0) return 0;
  ll ans = n;
  for (ll x = 2; x*x <= n; ++x) {
    if (n % x == 0) {
      ans -= ans / x;
      while (n % x == 0) n /= x;
    }
  }
  if (n > 1) ans -= ans / n;
  return ans;
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
vector<ll> euler_phi(ll lo, ll hi) { // phi(n) for all n in [lo, hi)
  vector<ll> ps = primes(sqrt(hi)+1);
  vector<ll> res(hi-lo), phi(hi-lo, 1);
  iota(all(res), lo);

  for (ll p: ps) {
    for (ll k = ((lo+(p-1))/p)*p; k < hi; k += p) {
      if (res[k-lo] < p) continue;
      phi[k-lo] *= (p - 1);
      res[k-lo] /= p;
      while (res[k-lo] > 1 && res[k-lo] % p == 0) {
        phi[k-lo] *= p;
        res[k-lo] /= p; 
      }
    }
  }
  for (ll k = lo; k < hi; ++k) {
    if (res[k-lo] > 1) 
      phi[k-lo] *= (res[k-lo]-1);
  }
  return phi; // phi[k-lo] = phi(k)
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

  for (int iter = 0; iter < 100; ++iter) {
    int lo = rand(), hi = lo + rand();
    auto x = euler_phi(lo, hi);
    for (int n = lo; n < hi; ++n)
      if (x[n-lo] != euler_phi(n)) cout << "!!" << endl;
  }
}

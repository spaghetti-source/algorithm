// 
// Solving a quadratic equation ax^2 + bx + c == 0
//
// Description
//   The solution is given by x = (-b \pm sqrt(D)) / 2a, where D = b^2 - 4ac.
//   To avoid the numerical errors, it should be computed by
//     if b > 0 then x1 = (-b - sqrt(D))/2a
//     otherwise     x1 = (-b + sqrt(D))/2a
//   and
//     x2 = c/(a1*x1)
//

#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath> 
#include <cstdio>
#include <cstring>
 
using namespace std;
#define fst first
#define snd second
#define all(c) ((c).begin()), ((c).end())

double quad_eqn(double a, double b, double c) { // assuming a != 0
  double D = b*b - 4*a*c, x1, x2;
  if (b > 0) x1 = (-b - sqrt(D))/(2*a);
  else       x1 = (-b + sqrt(D))/(2*a);
  x2 = c / (a * x1);
  return max(x1, x2);
}


// verify: SPOJ22329

typedef long long ll;
// number of eggs in x days for the hens with ability d
ll egg(ll d, ll x) {
  ll k = quad_eqn(1, 2*d-1, -2*x);
  while (k*k + (2*d-1)*k <= 2*x) ++k;
  while (k*k + (2*d-1)*k  > 2*x) --k;
  return k;
}
ll eggs(vector<ll> &ds, ll x) {
  ll ans = 0;
  for (auto d: ds) ans += egg(d, x);
  return ans;
}

int main() {
  int ncase; scanf("%d", &ncase);
  for (int icase = 0; icase < ncase; ++icase) {
    ll n, k; scanf("%lld %lld", &n, &k);
    vector<ll> ds(n);
    for (int i = 0; i < n; ++i) 
      scanf("%lld", &ds[i]);

    ll lo = 0, hi = 1; // x <= lo ==> eggs < k; x >= hi ==> eggs >= k
    while (eggs(ds, hi) < k) {
      lo = hi; 
      hi *= 2;
    }
    while (lo+1 < hi) {
      ll mi = (lo + hi) / 2;
      if (eggs(ds, mi) < k) lo = mi;
      else                  hi = mi;
    }
    printf("%lld\n", hi);
  }
}

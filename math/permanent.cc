// 
// Permanent
//
// Description:
//   Permanent of n x n matrix A is defined by
//     perm(A) := \sum_{sigma: permutation} \sum_{i=1}^n A_{i sigma(i)}.
//   Note that the determinant is defined by
//     det(A) := \sum_{sigma: permutation} sgn(sigma) \sum_{i=1}^n A_{i sigma(i)}.
//   Thus, these differs only sig(sigma) factor.
//   Computing permanent is known to be #P-hard (even for mod 3)
//
// Algorithm:
//   Ryser's inclusion exclusion princple. Let V = {1, ..., n}.
//     perm(A) = sum_{S subseteq V} (-1)^{|V - S|} a(S),
//   where
//     a(S) = prod_{j in V} sum_{i in S} A_{i j}. 
//
// Complexity:
//   O(n 2^n).
// 
// Verified:
//   SPOJ423
//
// References:
//   H. J. Ryser (1963):
//   Combinatorial Mathematics.
//   The Mathematical Association of America.

#include <iostream>
#include <vector>
#include <cstdio>
#include <algorithm>
#include <functional>
#include <unordered_map>

using namespace std;

#define fst first
#define snd second
#define all(c) ((c).begin()), ((c).end())

typedef long long ll;
typedef vector<ll> vec;
typedef vector<vec> mat;
ll permanent(mat A) {
  int n = A.size();
  vector<ll> a(n);   // row sum
  vector<int> io(n); // included or not 
  ll perm = 0;
  for (int i = 1; i < (1<<n); ++i) { 
    int k = __builtin_ffs(i) - 1; // = least significant bit
    ll dir = (io[k] ^= 1) ? +1 : -1;
    ll term = ((n - i) & 1 ? -1 : 1);
    for (int j = 0; j < n; ++j) {
      a[j] += dir * A[k][j];
      term *= a[j];
    }
    perm += term;
  }
  return perm;
}

int main() {
  int ncase = 0; scanf("%d", &ncase);
  for (int icase = 0; icase < ncase; ++icase) {
    int n; scanf("%d", &n);
    mat A(n, vec(n));
    for (int i = 0; i < n; ++i) 
      for (int j = 0; j < n; ++j) 
        scanf("%lld", &A[i][j]);
    printf("%lld\n", permanent(A));
  }
}

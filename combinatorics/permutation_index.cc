// 
// Permutation Index
//
// Description:
//   index_perm computes the lexicographical index of the permutation x of [0,n),
//   i.e., it gives k such that x = next_permutation^k [0,n).
//
//   unindex_perm is the inverse function of index_perm.
//
// Algorithm:
//   We represent index in the factorial number system:
//     index = d_0 (n-1)! + d_1 (n-2)! + ... + d_{n-1} 0!
//
//   index_perm:
//     for i = 0, ..., n-1: 
//       d[i] = the number of active numbers j < x[i] 
//       deactivate x[i]
//
//   unindex_perm:
//     for i = 0, ..., n-1:
//       x[i] = d[i]-th active element
//       deactivate x[i]
//   
//   Both can be computed in O(n log n) time by using Fenwick tree.
// 
// Complexity:
//   O(n log n) time, O(n) space.
//
// Note:
//   It correctly works for n <= 20 since 20! < 2^64.
//   
// References:
//   B. Bonet (2008):
//   Efficient algorithms to rank and unrank permutations in lexicographic order.
//   In Proceedings of the AAAI Workshop on Search in AI and Robotics, pp.18--23.

#include <iostream>
#include <vector>
#include <cstdio>
#include <algorithm>
#include <numeric>
#include <functional>

using namespace std;

#define fst first
#define snd second
#define all(c) ((c).begin()), ((c).end())

typedef long long ll;

ll index_perm(vector<int> x) {
  ll r = 0;
  vector<int> a(x.size()+1);
  for (int i = 0; i < x.size(); ++i) {
    int s = x[i];
    for (int k = x[i]; k > 0; k &= k-1) s -= a[k];
    r = (x.size() - i) * r + s;
    for (int k = x[i]+1; k < x.size(); k += k&-k) ++a[k];
  }
  return r;
}
vector<int> unindex_perm(ll r, int n) {
  vector<int> d(n), x(n, n);
  for (int i = n-1; i >= 0; --i) {
    d[i] = r % (n - i); r /= (n - i);
  }
  vector<int> a(n+1);
  for (int k = 1; k <= n; ++k) a[k] = k & -k;
  for (int i = 0; i < n; ++i) {
    for (int s: {1,2,4,8,16}) x[i] |= (x[i] >> s); 
    for (int p = ++x[i]; p > 0; p >>= 1, x[i] |= p)
      if (x[i] <= n && a[x[i]] <= d[i]) d[i] -= a[x[i]]; else x[i] ^= p;
    for (int k = x[i]+1; k < x.size(); k += k&-k) --a[k];
  }
  return x;
}

int main() {
  int n = 20;
  vector<int> x(n);
  iota(all(x), 0);
  for (int i = 0; i < 100; ++i) {
    random_shuffle(all(x));
    ll r = index_perm(x);
    auto a = unindex_perm(r, n);
    for (int i = 0; i < n; ++i) {
      if (a[i] != x[i]) {
        printf("wrong\n");
      }
    }
  }
}

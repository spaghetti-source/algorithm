// 
// Permutation Hash
//
// Description:
//   hash_perm gives one-to-one correspondence between 
//   permutations over [0,n) and the integer less than n!.
//   unhash_perm is the inverse function of hash_perm.
//
// Algorithm:
//   The idea is based on the Fisher-Yates shuffle algorithm:
//     while n > 1:
//       swap(x[n-1], x[rand() % n]);
//       --n;
//   For an integer given by a factorial number system:
//     hash = d_0 (n-1)! + d_1 (n-2)! + ... + d_{n-1} 0!
//   The algorithm computes 
//     while n > 1:
//       swap(x[n-1], x[d_{n-1}]
//       --n;
// 
// Complexity:
//   O(n) time, O(n) space.
//   
// Verification:
//   self. 


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

vector<int> unhash_perm(ll r, int n) {
  vector<int> x(n);
  iota(all(x), 0);
  for (; n > 0; --n) {
    swap(x[n-1], x[r % n]);
    r /= n;
  }
  return x;
}
ll hash_perm(vector<int> x) {
  int n = x.size();
  vector<int> y(n);
  for (int i = 0; i < n; ++i) y[x[i]] = i;
  ll c = 0, fac = 1;
  for (; n > 1; --n) {
    c += fac * x[n-1]; fac *= n;
    swap(x[n-1], x[y[n-1]]);
    swap(y[n-1], y[x[y[n-1]]]);
  }
  return c;
}

int main() {
  int n = 9;
  vector<int> x(n);
  iota(all(x), 0);
  do {
    ll r = hash_perm(x);
    cout << r << ": ";
    auto a = unhash_perm(r, n);
    for (int i = 0; i < n; ++i) {
      cout << a[i] << " ";
      if (a[i] != x[i]) exit(-1);
    }
    cout << endl;
  } while (next_permutation(all(x)));
  return 0;
}

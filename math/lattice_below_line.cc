// 
// Number of lattice points below a line 
//
// Description:
//
//   Let a, b, n, m be nonnegative integers. The task is to compute 
//     sum_{i in [0,n)} floor((a + ib)/m).
//
//   We compute this quantity in two directions alternately.
//   First, let 
//     a = (a/m) m + (a%m), 
//     b = (b/m) m + (b%m).
//   Then the quantity is
//     sum [(a/m)+i*(b/m)] + floor(((a%m) + i(b%m))/m)
//   Here, the first term is analytically evaluated.
//   The second term is zero if b%m == 0. Otherwise, the task is 
//   reduced to compute
//     sum_{i in [0,n)} floor((a + ib)/m)
//   where a < m, b < m. By changing the axes, this quantity is 
//     sum_{i in [0,n')} floor((a' + ib')/m')
//   where
//     n' = (a + b n) / m,
//     a' = (a + b n) % m,
//     b' = m,
//     m' = b.
//
//   We evaluate the number of iterations. Since the computation 
//   between b and m is the same as the one of the Euclidean 
//   algorithm. Thus it terminates in O(log m) time.
//
// Complexity:
//
//   O(log m).
//
// Verified:
// 
//   Somewhere

#include <iostream>
#include <vector>
#include <cstdio>
#include <algorithm>
#include <functional>

using namespace std;

#define fst first
#define snd second
#define all(c) ((c).begin()), ((c).end())

//
// sum_{0<=i<n} floor(a + b i)/m; assume n, m >= 0, a >= 0, b >= 0
//
// 
using Int = long long;
Int latticeBelowLine(Int n, Int a, Int b, Int m) {
  Int ans = 0;
  while (m) {
    ans += (n-1)*n/2*(b/m) + n*(a/m); 
    a %= m;
    b %= m;
    auto z = (a+b*n);
    a = z%m;
    n = z/m;
    swap(b, m);
  }
  return ans;
}

int main() {
  srand(time(0));
  Int a = rand(), b = rand(), n = rand(), m = rand();
  cout << latticeBelowLine(n, a, b, m) << endl;
}

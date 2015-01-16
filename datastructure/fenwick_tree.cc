//
// Fenwick Tree (aka. binary indexed tree)
//
// Description:
//   A data structure that allows 
//     add(k,a): x[k] = x[k] '+' a
//     sum(i,j): x[i] '+' x[i+1] '+' ... '+' x[j]
//   where '+' is an associative operator.
//
// Algorithm:
//     y        x
//     0 = sum [0]
//     1 = sum [0,1]
//     2 = sum [2]
//     3 = sum [0,1,2,3]
//     4 = sum [4]
//     5 = sum [4,5]
//     6 = sum [6]
//     7 = sum [0,1,2,3,4,5,6,7]
//
// Complexity:
//   O(log n) access, n space.
//
// Verified:
//   SPOJ3377

#include <iostream>
#include <vector>
#include <cstdio>
#include <algorithm>

using namespace std;

#define fst first
#define snd second
#define all(c) ((c).begin()), ((c).end())

template <class T>
struct fenwick_tree {
  vector<T> x;
  fenwick_tree(int n, T a = T(0)) : x(n, a) { }
  void add(int k, int a) { // x[k] += a
    for (; k < x.size(); k |= k+1) x[k] += a;
  }
  T sum(int k) { // return x[0] + ... + x[k]
    T s = 0;
    for (; k >= 0; k = (k&(k+1))-1) s += x[k];
    return s;
  }
};

int main() {
  int cases;
  scanf("%d", &cases);
  for (int icases = 0; icases < cases; ++icases) {
    int n, u;
    scanf("%d %d", &n, &u);

    fenwick_tree<int> FT(n+1);
    for (int k = 0; k < u; ++k) {
      int i, j, a;
      scanf("%d %d %d", &i, &j, &a);
      FT.add(i,    a);
      FT.add(j+1, -a);
    }
    int q;
    scanf("%d", &q);
    for (int k = 0; k < q; ++k) {
      int i;
      scanf("%d", &i);
      printf("%d\n", FT.get(i));
    }
  }
}

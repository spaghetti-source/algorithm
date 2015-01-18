//
// 2D Fenwick Tree 
//
// Description:
//   A data structure that allows 
//     add(k,a): x[k] += a
//     sum(k):   sum { x[i] : i <= k }
//   where k denotes a point in [0,n)x[0,m).
//
// Algorithm:
//   Fenwick tree of Fenwick tree.
//
// Complexity:
//   O((log n)^2) time, O(n^2) space.
//
// Verified:
//   SPOJ1029

#include <iostream>
#include <vector>
#include <cstdio>
#include <algorithm>

using namespace std;

#define fst first
#define snd second
#define all(c) ((c).begin()), ((c).end())

template <class T>
struct fenwick_tree_2d {
  vector<vector<T>> x;
  fenwick_tree_2d(int n, int m) : x(n, vector<T>(m)) { }
  void add(int k1, int k2, int a) { // x[k] += a
    for (; k1 < x.size(); k1 |= k1+1) 
      for (int k=k2; k < x[k1].size(); k |= k+1) x[k1][k] += a;
  }
  T sum(int k1, int k2) { // return x[0] + ... + x[k]
    T s = 0;
    for (; k1 >= 0; k1 = (k1&(k1+1))-1) 
      for (int k=k2; k >= 0; k = (k&(k+1))-1) s += x[k1][k];
    return s;
  }
};

void doit() {
  int n; scanf("%d", &n);
  fenwick_tree_2d<int> FT(n,n);
  while (1) {
    char cmd[5];
    scanf("%s", cmd);
    switch (cmd[1]) {
      case 'E': {
        int i, j, a;
        scanf("%d %d %d", &i, &j, &a);
        int b = FT.sum(i, j) + FT.sum(i-1,j-1) - FT.sum(i-1,j) - FT.sum(i,j-1);
        FT.add(i, j, a - b);
      }
      break;
      case 'U': {
        int i, j, k, l;
        scanf("%d %d %d %d", &i, &j, &k, &l);
        int b = FT.sum(k, l) + FT.sum(i-1,j-1) - FT.sum(i-1,l) - FT.sum(k,j-1);
        printf("%d\n", b);
        break;
      }
      case 'N': {
        printf("\n");
        return;
      }
    }
  }
}
int main() {
  int cases;
  scanf("%d", &cases);
  for (int icases = 0; icases < cases; ++icases) doit();
}

//
// Knapsack Problem (branch-and-bound with expanding core)
//
// Description:
//   We are given a set of items with profit p_i and weight w_i.
//   The problem is to find a subset of items that maximizes
//   the total profit under the total weight less than some capacity c.
//
//   1) c is small     ==> capacity DP
//   2) p is small     ==> price DP
//   3) both are large ==> branch and bound.
// 
// Algorithm:
//   Branch and bound with expanding core.
//   We first sort the items by p_i/w_i in descending order.
//   Then, the fractional solution is given by selecting
//   integral {1 ... b-1} and fractional b.
//   Branch-and-bound method recursively finds a solution for b = 0 and 1.
//
//   For an efficient implementation, we maintain a interval [s,t];
//   which means that all items i <= s is selected, and all items j >= t
//   is not selected. The algorithm recursively select s or discard t.
//   
//   Intuitively, the algorithm enumerates all possibilities in [s,t].
//   The set of items in [s,t] is called "core." 
//
// Complexity:
//   O(2^c), where c is the size of core. Basically, c is not so large.
//
// Verified:
//   SPOJ3321.
//
// References:
//   H. Kellerer, U. Pferschy, and D. Pisinger (2004):
//   Knapsack problems.
//   Springer Science & Business Media.

#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <vector>
#include <ctime>
#include <functional>
#include <queue>
#include <algorithm>
#include <numeric>

using namespace std;

#define all(c) c.begin(),c.end()

template <class T>
struct knapsack {
  T c; 
  struct item { T p, w; }; // price/weight
  vector<item> is;
  void add_item(T p, T w) {
    is.push_back({p, w});
  }
  T det(T a, T b, T c, T d) {
    return a * d - b * c;
  }
  T z;
  void expbranch(T p, T w, int s, int t) {
    if (w <= c) {
      if (p >= z) z = p;
      for (; t < is.size(); ++t) {
        if (det(p - z - 1, w - c, is[t].p, is[t].w) < 0) return;
        expbranch(p + is[t].p, w + is[t].w, s, t + 1);
      }
    } else {
      for (; s >= 0; --s) {
        if (det(p - z - 1, w - c, is[s].p, is[s].w) < 0) return;
        expbranch(p - is[s].p, w - is[s].w, s - 1, t);
      }
    }
  }
  T solve() {
    sort(all(is), [](const item &a, const item &b) { 
      return a.p * b.w > a.w * b.p;
    });
    T p = 0, w = 0;
    z = 0;
    int b = 0; 
    for (; b < is.size() && w <= c; ++b) {
      p += is[b].p;
      w += is[b].w;
    }
    expbranch(p, w, b-1, b);
    return z;
  }
};

int main() {
  int s, n;
  scanf("%d %d", &s, &n);
  knapsack<int> solver;
  solver.c = s;
  for (int i = 0; i < n; ++i) {
    int v, w; 
    scanf("%d %d", &w, &v);
    solver.add_item(v, w);
  }
  printf("%d\n", solver.solve());
}

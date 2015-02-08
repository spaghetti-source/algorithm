// 
// Interval Scheduling (weighted)
//
// Description:
//   We are given a set of intervals [b_i, e_i] with weight w_i, i = 1, ..., n.
//   Find a set of disjoint intervals with maximum weight.
//
// Algorithm:
//   Dynamic programming. Suppose e_1 <= e_2 <= ... <= e_n.
//   Let f(k) be the optimal score for intervals {1, ..., k}.
//   Then, we have
//     f(k) = max( f(p(k)), f(k-1) ),
//   where 
//     p(k) = max { j : e_j <= b_k }.
//   p(k) is computed in O(log n). Thereore the above DP is 
//   performed in O(n log n).
//     
// Complexity:
//   O(n log n).
//
// Verified:
//   SPOJ11515.

#include <iostream>
#include <vector>
#include <cstdio>
#include <algorithm>
#include <functional>

using namespace std;

#define fst first
#define snd second
#define all(c) ((c).begin()), ((c).end())

template <class T>
struct weighted_interval_scheduling {
  struct interval { T b, e, w; };
  vector<interval> is;
  void add_interval(T b, T e, T w) {
    is.push_back({b, e, w});
  }
  T max_scheduling() {
    int n = is.size();
    sort(all(is), [](interval x, interval y) { return x.e < y.e; });
    vector<int> p(n);
    for (int i = 0; i < n; ++i) {
      auto cond = [](T key, interval x) { return key < x.e; };
      p[i] = upper_bound(all(is), is[i].b, cond) - is.begin();
    }
    vector<T> f(n+1);
    for (int i = 0; i < n; ++i) 
      f[i+1] = max(f[p[i]] + is[i].w, f[i]);
    return f[n];
  }
};

int main() {
  int ncase; scanf("%d", &ncase);
  for (int icase = 0; icase < ncase; ++icase) {
    int n; scanf("%d", &n);
    weighted_interval_scheduling<int> is;
    for (int i = 0; i < n; ++i) {
      int u, v;
      scanf("%d %d", &u, &v);
      is.add_interval(u, v, 1);
    }
    printf("%d\n", is.max_scheduling());
  }
}

// 
// Interval Scheduling (unweighted)
//
// Description:
//   We are given a set of intervals [b_i, e_i], i = 1, ..., n.
//   Find a set of disjoint intervals with maximum cardinality.
//
// Algorithm:
//   Greedy. Sort by e_i and then take the intervals greedily.
// 
//   Consider an optimal scheduling and I = [b, e] be the first interval.
//   If there exists an interval I' = [b', e'] with e' < e, then
//   we can replace I by I' without increasing the cardinality.
//   This shows that there exists a solution that contains the interval
//   with the earliest end point.
//     
// Complexity:
//   O(n log n).
//

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
struct unweighted_interval_scheduling {
  struct interval { T b, e; };
  vector<interval> is;
  void add_interval(T b, T e) {
    is.push_back({b, e});
  }
  T max_scheduling() {
    sort(all(is), [](interval x, interval y) {
      return x.e < y.e;
    });
    int score = 0;
    T sweep = is[0].b - 1;
    for (int i = 0; i < is.size(); ++i) {
      if (is[i].b >= sweep) {
        ++score;
        sweep = is[i].e;
      }
    }
    return score;
  }
};

int main() {
  int ncase; scanf("%d", &ncase);
  for (int icase = 0; icase < ncase; ++icase) {
    int n; scanf("%d", &n);
    unweighted_interval_scheduling<int> is;
    for (int i = 0; i < n; ++i) {
      int u, v;
      scanf("%d %d", &u, &v);
      is.add_interval(u, v);
    }
    printf("%d\n", is.max_scheduling());
  }
}

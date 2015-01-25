//
// Covered Range
//
// Description:
//   Given a set of intervals [l_j, r_j).
//   Find a measure of union of the intervals.
//
// Algorithm:
//   Plane sweep.
//
// Complexity:
//   O(n log n).
//
// Verified:
//   SPOJ18531
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

// measure of union of x[j].
template <class T>
T covered_range(vector<pair<T, T>> x) {
  typedef pair<T, int> event;
  vector<event> es;
  for (int i = 0; i < x.size(); ++i) {
    es.push_back({x[i].fst, i});
    es.push_back({x[i].snd,~i});
  }
  sort(all(es));
  int c = 0;
  T a = es[0].fst, ans = 0;
  for (auto e: es) {
    if (c > 0) ans += e.fst - a;
    if (e.snd >= 0) ++c;
    else            --c;
    a = e.fst;
  }
  return ans;
}

int main() {
  vector<pair<int,int>> x;
  x.push_back({0,3});
  x.push_back({2,5});
  x.push_back({4,6});
  cout << covered_range(x) << endl;
}

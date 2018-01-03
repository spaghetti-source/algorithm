// 
// Partially Persistent Union Find
//
// Description:
//   It is a persistent version of union find data structure.
//   It allows us to apply "find" to the all versions and 
//   "unite" to the latest version.
//
// Complexity:
//   O(log n) for each query.
//
// Verified:
//   CODE THANKS FESTIVAL 2017, Union Set: 
//     https://code-thanks-festival-2017-open.contest.atcoder.jp/tasks/code_thanks_festival_2017_h
//

#include <bits/stdc++.h>
using namespace std;

#define fst first
#define snd second
#define all(c) ((c).begin()), ((c).end())
using namespace std;

struct PartiallyPersistentUnionFind {
  vector<vector<pair<int, int>>> parent; // (parent index, modified time)
  int now = 0; // time = 0 is the initial state
  PartiallyPersistentUnionFind(int n) : parent(n, {{-1,0}}) { }
  bool unite(int u, int v) { 
    ++now;
    u = root(u, now); v = root(v, now);
    if (u == v) return false;
    if (parent[u].back().fst > parent[v].back().fst) swap(u, v);
    parent[u].push_back({parent[u].back().fst+parent[v].back().fst, now});
    parent[v].push_back({u, now});
    return true;
  }
  bool find(int u, int v, int t) { return root(u, t) == root(v, t); }
  int root(int u, int t) { 
    if (parent[u].back().fst >= 0 && parent[u].back().snd <= t) 
      return root(parent[u].back().fst, t);
    return u;
  }
  int size(int u, int t) { 
    u = root(u, t);
    int lo = 0, hi = parent[u].size();
    while (lo + 1 < hi) {
      int mi = (lo + hi) / 2;
      if (parent[u][mi].snd <= t) lo = mi;
      else                        hi = mi;
    }
    return -parent[u][lo].fst;
  }
};

int main() {
  int n, m, q, a, b;
  scanf("%d %d", &n, &m);

  PartiallyPersistentUnionFind uf(n);
  for (int i = 0; i < m; ++i) {
    scanf("%d %d", &a, &b); --a; --b;
    uf.unite(a, b);
  }
  scanf("%d", &q);
  for (int i = 0; i < q; ++i) {
    scanf("%d %d", &a, &b); --a; --b;
    int lo = 0, hi = uf.now+1;
    while (lo+1 < hi) {
      int mi = (lo + hi) / 2;
      if (uf.find(a, b, mi)) hi = mi;
      else                   lo = mi;
    }
    if (hi == uf.now+1) printf("-1\n");
    else                printf("%d\n", hi);
  }
}

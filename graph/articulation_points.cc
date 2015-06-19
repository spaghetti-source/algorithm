// 
// Articulation points / Biconnected components
//
// Description:
//   Let G = (V, E). If G-v is disconnected, v in V is said to
//   be an articulation point. If G has no articulation points,
//   it is said to be biconnected.
//   A biconnected component is a maximal biconnected subgraph.
//   The algorithm finds all articulation points and biconnected 
//   components.
//
// Algorithm:
//   Hopcroft-Tarjan's DFS based algorithm.
//
// Complexity:
//   O(n + m).
//
// Verified:
//   SPOJ 14956: Submerging Island (articulation point)
//   POJ 2942: Knights of the Round Table (biconnected components)
//
// References:
//   J. Hopcroft and R. E. Tarjan (1973):
//   Efficient algorithms for graph manipulation.
//   Communications of the ACM, vol.16, no.6, pp.372-378.
//

#include <iostream>
#include <vector>
#include <cstdio>
#include <algorithm>
#include <functional>
#include <unordered_set>

using namespace std;

#define fst first
#define snd second
#define all(c) ((c).begin()), ((c).end())

struct graph {
  int n;
  vector<vector<int>> adj;
  graph(int n) : n(n), adj(n) { }
  void add_edge(int src, int dst) {
    adj[src].push_back(dst);
    adj[dst].push_back(src);
  }
};

void biconnected_components(graph g) {
  vector<int> arts(g.n), num(g.n), low(g.n), S;
  vector<vector<int>> comps;
  function<void(int,int,int&)> dfs = [&](int p, int i, int &t) { 
    num[i] = low[i] = ++t; 
    S.push_back(i);
    for (int j: g.adj[i]) {
      if (j == p) continue;
      if (num[j] == 0) {
        dfs(i, j, t);
        low[i] = min(low[i], low[j]);
        if (num[i] <= low[j]) {
          if (num[i] != 1 || low[j] > 2) arts[i] = true;
          for (comps.push_back({i}); comps.back().back() != j; S.pop_back()) 
            comps.back().push_back(S.back()); // 
        }
      } else low[i] = min(low[i], num[j]);
    }
  };
  for (int i = 0, t; i < g.n; ++i)
    if (num[i] == 0) dfs(-1, i, t = 0);

  // SPOJ SUBMERGE
  int count = 0;
  for (int i = 0; i < g.n; ++i) 
    count += arts[i];
  printf("%d\n", count);
}

int main() {
  for (int n, m; ~scanf("%d %d", &n, &m) && n; ) {
    graph g(n);
    for (int i = 0; i < m; ++i) {
      int u, v; scanf("%d %d", &u, &v);
      g.add_edge(u-1, v-1);
    }
    biconnected_components(g);
  }
}

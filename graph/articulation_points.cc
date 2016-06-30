// 
// Articulation points / Biconnected components
//
// Description:
//   Let G = (V, E). If G-v is disconnected, v in V is said to
//   be an articulation point. If G has no articulation points,
//   it is said to be biconnected.
//
//   A biconnected component is a maximal biconnected subgraph.
//   The algorithm finds all articulation points and biconnected 
//   components.
//
//   The most important fact is that by contracting biconnected
//   components we obtain a tree, which is called the block tree.
//
//
// Algorithm:
//   Hopcroft-Tarjan's DFS based algorithm.
//
//   Single DFS finds a block tree rooted from the component
//   that contains the specified root.
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

  void biconnected_components() {
    vector<int> num(n), low(n), S;
    unordered_set<int> arts;

    function<void(int,int,int&)> dfs = [&](int p, int u, int &t) { 
      num[u] = low[u] = ++t; 
      S.push_back(u);
      for (int v: adj[u]) {
        if (v == p) continue;
        if (num[v] == 0) {
          dfs(u, v, t);
          low[u] = min(low[u], low[v]);
          if (num[u] <= low[v]) {
            if (num[u] != 1 || num[v] > 2) {
              // here, u is an articulation point if
              //   (a). u is non-root
              //   (b). u is root with two more children 
            }
            vector<int> C = {u}; // biconnected component
            while (C.back() != v) {
              C.push_back(S.back());
              S.pop_back();
            }
          }
        } else low[u] = min(low[u], num[v]);
      }
    };
    for (int u = 0, t; u < n; ++u) 
      if (!num[u]) dfs(-1, u, t = 0);
    cout << arts.size() << endl;
  }
};

int main() {
  for (int n, m; ~scanf("%d %d", &n, &m) && n; ) {
    graph g(n);
    for (int i = 0; i < m; ++i) {
      int u, v; scanf("%d %d", &u, &v);
      g.add_edge(u-1, v-1);
    }
    g.biconnected_components();
  }
}

// 
// Bipartite graph recognition
//
// Description:
//   A graph is bipartite if there is a 2-partition 
//   such that there are no edges in the same components.
//
// Algorithm:
//   Depth first search.
//
// Verified:
//   SPOJ3337
//
#include <iostream>
#include <cstdio>
#include <string>
#include <vector>
#include <iterator>
#include <functional>
#include <algorithm>

using namespace std;

struct graph {
  int n;
  vector<vector<int>> adj;
  graph(int n = 0) : n(n), adj(n) { }
  void add_edge(int src, int dst) {
    n = max(n, max(src, dst)+1);
    adj.resize(n);
    adj[src].push_back(dst);
    adj[dst].push_back(src);
  }
};
bool is_bipartite(graph g) {
  vector<int> color(g.n, -1);
  for (int u = 0; u < g.n; ++u) {
    if (color[u] != -1) continue;
    color[u] = 0;
    for (vector<int> S = {u}; !S.empty(); ) {
      int v = S.back(); S.pop_back();
      for (auto w: g.adj[v]) {
        if (color[w] == color[v]) return false;
        if (color[w] == -1) {
          color[w] = !color[v];
          S.push_back(w);
        }
      }
    }
  }
  return true;
}


int main() {
  int ncase;
  scanf("%d", &ncase);
  for (int icase = 0; icase < ncase; ++icase) {
    printf("Scenario #%d:\n", icase+1);
    int n, m;
    scanf("%d %d", &n, &m);
    graph g(n);
    for (int i = 0; i < m; ++i) {
      int u, v;
      scanf("%d %d", &u, &v);
      g.add_edge(u-1, v-1);
    }
    if (is_bipartite(g)) printf("No suspicious bugs found!\n");
    else                 printf("Suspicious bugs found!\n");
  }
}

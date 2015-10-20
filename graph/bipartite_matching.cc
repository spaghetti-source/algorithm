//
// Ford-Fulkerson' maximum bipartite matching
//
// Description:
//   Compute the maximum cardinality matching for bipartite graph.
//
// Algorithm:
//   Ford-Fulkerson type (DFS-based) augmentaing path algorithm.
//
// Complexity:
//   O(m n) time
//
// Verified:
//   AOJ Matching
// 
// Note: 
//   TLE in SPOJ 4206: Fast Maximum Matching
//

#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <vector>
#include <queue>
#include <cstring>
#include <functional>
#include <algorithm>

using namespace std;
#define all(c) (c).begin(), (c).end()

struct graph {
  int L, R;
  vector<vector<int>> adj;
  graph(int L, int R) : L(L), R(R), adj(L+R) { }
  void add_edge(int u, int v) {
    adj[u].push_back(v+L);
    adj[v+L].push_back(u);
  }
  int maximum_matching() {
    vector<int> visited(L), mate(L+R, -1);
    function<bool(int)> augment = [&](int u) { // DFS
      if (visited[u]) return false;
      visited[u] = true;
      for (int w: adj[u]) {
        int v = mate[w];
        if (v < 0 || augment(v)) {
          mate[u] = w;
          mate[w] = u;
          return true;
        }
      }
      return false;
    };
    int match = 0;
    for (int u = 0; u < L; ++u) {
      fill(all(visited), 0);
      if (augment(u)) ++match;
    }
    return match;
  }
};

int main() {
  int L, R, m; 
  scanf("%d %d %d", &L, &R, &m);
  graph g(L, R);
  for (int i = 0; i < m; ++i) {
    int u, v;
    scanf("%d %d", &u, &v);
    g.add_edge(u, v);
  }
  printf("%d\n", g.maximum_matching());
}

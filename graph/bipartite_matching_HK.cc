//
// Hopcroft-Karp's maximum cardinality bipartite matching
//
// Description:
//   Compute the maximum cardinality matching for bipartite graph.
//
// Algorithm:
//   The algorithm iterates following procedures:
//     (1) BFS from the source to get the distance to the sink.
//         If not reachable, there are no augment path hence break.
//     (2) Find vertex disjoint shortest augment paths by DFS.
//   It can be shown that the outer-loop is atmost O(\sqrt{n}) times
//   therefore the whole complexity is O(m \sqrt{n}).
//   Note that this is a specialzation of Dinic's maximum flow.
//
//
// Complexity:
//   O(m \sqrt{n}) time
//
// Verified:
//   SPOJ 4206: Fast Maximum Matching
//
// References:
//   J. E. Hopcroft and R. M. Karp (1973):
//   An n^5/2 algorithm for maximum matchings in bipartite graphs.
//   SIAM Journal on Computing, vol.2, no.4, pp.225-231.
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

struct graph {
  int L, R;
  vector<vector<int>> adj;
  graph(int L, int R) : L(L), R(R), adj(L+R) { }
  void add_edge(int u, int v) {
    adj[u].push_back(v+L);
    adj[v+L].push_back(u);
  }
  int maximum_matching() {
    vector<int> level(L), mate(L+R, -1);

    function<bool(void)> levelize = [&]() { // BFS
      queue<int> Q;
      for (int u = 0; u < L; ++u) {
        level[u] = -1;
        if (mate[u] < 0) {
          level[u] = 0;
          Q.push(u); 
        }
      }
      while (!Q.empty()) {
        int u = Q.front(); Q.pop();
        for (int w: adj[u]) {
          int v = mate[w];
          if (v < 0) return true;
          if (level[v] < 0) {
            level[v] = level[u] + 1;
            Q.push(v); 
          }
        }
      }
      return false;
    };
    function<bool(int)> augment = [&](int u) { // DFS
      for (int w: adj[u]) {
        int v = mate[w];
        if (v < 0 || (level[v] > level[u] && augment(v))) {
          mate[u] = w;
          mate[w] = u;
          return true;
        }
      }
      return false;
    };
    int match = 0;
    while (levelize()) 
      for (int u = 0; u < L; ++u) 
        if (mate[u] < 0 && augment(u)) 
          ++match;
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

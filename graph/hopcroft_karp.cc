//
// Hopcroft-Karp's maximum cardinality bipartite matching
//
// Description:
//   Compute the maximum cardinality matching for bipartite graph.
//
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
//
// Verified:
//   SPOJ 4206: Fast Maximum Matching
//
//
// References:
//   J. E. Hopcroft and R. M. Karp (1973):
//   An n^5/2 algorithm for maximum matchings in bipartite graphs.
//   SIAM Journal on Computing, vol. 2, no. 4, pp. 225-231.
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

struct edge {
  int src, dst;
};
int L, R;
vector<vector<edge>> adj;
void add_edge(int u, int v) {
  v += L;
  adj[u].push_back({u, v});
  adj[v].push_back({v, u});
}
vector<int> mate, level;
bool levelize() { 
  level.assign(L, -1);
  queue<int> Q;
  for (int u = 0; u < L; ++u) {
    if (mate[u] == -1) {
      Q.push(u); 
      level[u] = 0;
    }
  }
  while (!Q.empty()) {
    int u = Q.front(); Q.pop();
    for (auto e: adj[u]) {
      int v = mate[e.dst];
      if (v < 0) return true;
      if (level[v] < 0) {
        Q.push(v); 
        level[v] = level[u] + 1;
      }
    }
  }
  return false;
}
bool augment(int u) {
  for (auto e: adj[u]) {
    int v = mate[e.dst];
    if (v < 0 || (level[v] > level[u] && augment(v))) {
      mate[e.src] = e.dst;
      mate[e.dst] = e.src;
      return true;
    }
  }
  return false;
}
int maximum_matching() {
  mate.assign(L+R, -1);
  int match = 0;
  while (levelize()) 
    for (int u = 0; u < L; ++u) 
      if (mate[u] == -1 && augment(u)) ++match;
  return match;
}

int main() {
  int L, R, m; 
  scanf("%d %d %d", &L, &R, &m);
  for (int i = 0; i < m; ++i) {
    int u, v;
    scanf("%d %d", &u, &v);
    add_edge(u-1, v-1);
  }
  printf("%d\n", maximum_matching());
}

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

#define ALL(c) c.begin(), c.end()
#define FOR(i,c) for(typeof(c.begin())i=c.begin();i!=c.end();++i)
#define REP(i,n) for(int i=0;i<n;++i)
#define fst first
#define snd second

struct Edge {
  int src, dst;
  Edge(int src, int dst) : src(src), dst(dst) { }
};
struct BipartiteGraph {
  int L, R, n;
  vector< vector<Edge> > adj;
  BipartiteGraph(int L, int R) : L(L), R(R), n(L+R), adj(n) { }
  void addEdge(int u, int v) {
    v += L;
    adj[u].push_back(Edge(u, v));
    adj[v].push_back(Edge(v, u));
  }

  vector<int> visited, mate, dist;
  bool levelize() { // true iff reachable
    dist.assign(L, -1);
    queue<int> Q;
    REP(u, L) if (mate[u] == -1) {
      Q.push(u); 
      dist[u] = 0;
    }
    while (!Q.empty()) {
      int u = Q.front(); Q.pop();
      FOR(e, adj[u]) {
        int v = mate[e->dst];
        if (v < 0) return true;
        if (dist[v] < 0) {
          Q.push(v); 
          dist[v] = dist[u] + 1;
        }
      }
    }
    return false;
  }
  bool augment(int u) {
    if (visited[u]++) return false;
    FOR(e, adj[u]) {
      int v = mate[e->dst];
      if (v < 0 || (dist[v] > dist[u] && augment(v))) {
        mate[e->src] = e->dst;
        mate[e->dst] = e->src;
        return true;
      }
    }
    return false;
  }
  int maximumMatching() {
    mate.assign(n, -1);
    int match = 0;
    while (levelize()) {
      visited.assign(L, 0);
      REP(u, L) if (mate[u] == -1 && augment(u)) ++match;
    } 
    return match;
  }
};

int main() {
  int L, R, m; 
  scanf("%d %d %d", &L, &R, &m);
  BipartiteGraph G(L, R);
  REP(i, m) {
    int u, v;
    scanf("%d %d", &u, &v);
    G.addEdge(u-1, v-1);
  }
  printf("%d\n", G.maximumMatching());
}

//
// Maximum Flow (Dinic)
// 
// Description:
//   Given a directed network G = (V, E) with edge capacity c: E->R.
//   The algorithm finds a maximum flow. 
//
// Algorithm:
//   Dinic's blocking flow algorithm.
//
// Complexity:
//   O(n^2 m), but very fast in practice.
// 
// Verified:
//   SPOJ FASTFLOW
//
// Reference:
//   E. A. Dinic (1970):
//   Algorithm for solution of a problem of maximum flow in networks with power estimation.
//   Soviet Mathematics Doklady, vol. 11, pp. 1277-1280.
//
//   B. H. Korte and Jens Vygen (2008):
//   Combinatorial Optimization: Theory and Algorithms.
//   Springer Berlin Heidelberg. 
//

#include <iostream>
#include <vector>
#include <queue>
#include <cstdio>
#include <algorithm>

using namespace std;

typedef long long flow_type;
flow_type INF = (1ULL<<50);
struct edge {
  int src, dst;
  flow_type capacity;
  int rev;
  flow_type residue;
};
struct graph {
  int n;
  vector<vector<edge>> adj;
  graph(int n = 0) : n(n), adj(n) { }
  void add_edge(int src, int dst, flow_type capacity) {
    n = max(n, max(src, dst)+1);
    adj.resize(n);
    adj[src].push_back({src, dst, capacity, (int)adj[dst].size()});
    adj[dst].push_back({dst, src, capacity, (int)adj[src].size()-1});
  }
  vector<int> level, iter;
  flow_type augment(int u, int t, flow_type cur) {
    if (u == t) return cur;
    for (int &i = iter[u]; i < adj[u].size(); ++i) {
      edge &e = adj[u][i];
      if (e.residue > 0 && level[u] < level[e.dst]) {
        flow_type f = augment(e.dst, t, min(cur, e.residue));
        if (f > 0) {
          e.residue -= f;
          adj[e.dst][e.rev].residue += f;
          return f;
        }
      }
    }
    return 0;
  }
  int bfs(int s, int t) {
    level.assign(n, -1); level[s] = 0;
    queue<int> Q; Q.push(s);
    while (!Q.empty()) {
      int u = Q.front(); Q.pop();
      if (u == t) break;
      for (auto &e: adj[u]) {
        if (e.residue > 0 && level[e.dst] < 0) {
          Q.push(e.dst);
          level[e.dst] = level[u] + 1;
        }
      }
    }
    return level[t];
  }
  flow_type max_flow(int s, int t) {
    for (int u = 0; u < n; ++u) // initialize
      for (auto &e: adj[u]) e.residue = e.capacity;
    flow_type flow = 0;
    int itera = 0;
    while (bfs(s, t) >= 0) {
      iter.assign(n, 0);
      for (flow_type f; (f = augment(s, t, INF)) > 0; )
        flow += f;
    } // level[u] == -1 ==> t-side
    return flow;
  }
};

int main() {
  for (int n, m; scanf("%d %d", &n, &m) == 2; ) {
    graph g(n);
    for (int i = 0; i < m; ++i) {
      int u, v, w;
      scanf("%d %d %d", &u, &v, &w);
      g.add_edge(u-1, v-1, w);
    }
    printf("%lld\n", g.max_flow(0, n-1));
  }
}

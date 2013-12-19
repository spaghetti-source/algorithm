//
// Dinic's maximum flow
//
//
// Description:
//   Compute the maximum flow.
//
//
// Algorithm:
//   The algorithm iterates following procedures:
//     (1) BFS from the source to get the distance to the sink.
//         If not reachable, there are no augment path hence break.
//     (2) Find vertex disjoint shortest augment paths by DFS.
//   Since the length of augmenting path increases at least 1, for each outer-loop,
//   the number of outer-loop is O(n). For each inner-loop, it takes O(nm) because
//   it consists of one BFS and O(n) DFSs. Therefore the whole complexity is O(n^2 m).
//
//
// Complexity:
//   O(n^2 m) time
//
//
// Verified:
//   PKU 1459, Power Network
//   SPOJ 4410, Fast Maximum Flow
//
// References:
//   E. A. Dinic (1970):
//   Algorithm for solution of a problem of maximum flow in networks with power estimation.
//   Soviet Mathematics Doklady, vol. 11, pp. 1277-1280.

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
  flow_type residue;
  int rev;
};
struct maximum_flow {
  vector<edge> edges;
  void add_edge(int src, int dst, flow_type capacity) {
    edges.push_back({src, dst, capacity});
  }
  int n;
  vector<vector<edge>> adj;
  void make_graph(int n_ = 0) {
    n = n_;
    for (auto e: edges) 
      n = max(n, max(e.src, e.dst)+1);
    adj.resize(n);
    for (auto e: edges) {
      edge r = {e.dst, e.src, 0}; 
      adj[e.src].push_back(e);
      adj[r.src].push_back(r);
      adj[e.src].back().rev = adj[e.dst].size()-1;
      adj[r.src].back().rev = adj[r.dst].size()-1;
    }
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
    level.assign(n, -1);
    level[s] = 0;
    queue<int> Q;
    Q.push(s);
    while (!Q.empty()) {
      int u = Q.front(); Q.pop();
      if (u == t) break;
      for (int i = 0; i < adj[u].size(); ++i) {
        edge &e = adj[u][i];
        if (e.residue > 0 && level[e.dst] < 0) {
          Q.push(e.dst);
          level[e.dst] = level[u] + 1;
        }
      }
    }
    return level[t];
  }
  flow_type solve(int s, int t) {
    flow_type flow = 0;
    while (bfs(s, t) >= 0) {
      iter.assign(n, 0);
      for (flow_type f; (f = augment(s, t, INF)) > 0; )
        flow += f;
    }
    return flow;
  }
};

int main() {
  int n, m; 
  scanf("%d %d", &n, &m);
  maximum_flow mf;
  for (int i = 0; i < m; ++i) {
    int u, v, w;
    scanf("%d %d %d", &u, &v, &w);
    mf.add_edge(u-1, v-1, w);
  }
  mf.make_graph(n);
  printf("%lld\n", mf.solve(0, n-1));
}

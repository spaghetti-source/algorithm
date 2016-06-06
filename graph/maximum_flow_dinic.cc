//
// Dinic's maximum flow
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
//   In particular, for a unit capacity graph, 
//   it runs in O(m min{m^{1/2}, n^{2/3}}).
// 
// Verified:
//   SPOJ FASTFLOW
//
// Reference:
//   E. A. Dinic (1970):
//   Algorithm for solution of a problem of maximum flow in networks with power estimation.
//   Soviet Mathematics Doklady, vol. 11, pp. 1277-1280.
//
//   B. H. Korte and J. Vygen (2008):
//   Combinatorial Optimization: Theory and Algorithms.
//   Springer Berlin Heidelberg. 
//

#include <iostream>
#include <vector>
#include <cstdio>
#include <queue>
#include <algorithm>
#include <functional>

using namespace std;

#define fst first
#define snd second
#define all(c) ((c).begin()), ((c).end())

const long long INF = (1ll << 50);
struct graph {
  int n;
  struct edge {
    int src, dst;
    long long capacity, residue;
    size_t rev;
  };
  vector<vector<edge>> adj;
  graph(int n) : n(n), adj(n) { }
  void add_edge(int src, int dst, long long capacity) {
    adj[src].push_back({src, dst, capacity, 0, adj[dst].size()});
    adj[dst].push_back({dst, src, capacity, 0, adj[src].size()-1}); // bidirectional edge
  }
  long long max_flow(int s, int t) {
    vector<int> level(n), iter(n);
    function<int(void)> levelize = [&]() {
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
    };
    function<long long(int, long long)> augment = [&](int u, long long cur) {
      if (u == t) return cur;
      for (int &i = iter[u]; i < adj[u].size(); ++i) {
        edge &e = adj[u][i];
        if (e.residue > 0 && level[u] < level[e.dst]) {
          long long f = augment(e.dst, min(cur, e.residue));
          if (f > 0) {
            e.residue -= f;
            adj[e.dst][e.rev].residue += f;
            return f;
          }
        }
      }
      return 0ll;
    };
    for (int u = 0; u < n; ++u) // initialize
      for (auto &e: adj[u]) e.residue = e.capacity;
    long long total = 0;
    while (levelize() >= 0) {
      fill(all(iter), 0);
      for (long long f; (f = augment(s, INF)) > 0; )
        total += f;
    } // level[u] == -1 ==> t-side
    return total;
  }
};

int main() {
  for (int n, m; scanf("%d %d", &n, &m) == 2; ) {
    graph g(n);
    for (int i = 0; i < m; ++i) {
      int u, v, w;
      scanf("%d %d %d", &u, &v, &w);
      //g.add_edge(u, v, w);
      g.add_edge(u-1, v-1, w);
    }
    printf("%lld\n", g.max_flow(0, n-1));
  }
}

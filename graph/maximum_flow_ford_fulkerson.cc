//
// Ford-Fulkerson's maximum flow 
// 
// Description:
//   Given a directed network G = (V, E) with edge capacity c: E->R.
//   The algorithm finds a maximum flow. 
//
// Algorithm:
//   Ford-Fulkerson's augmenting path algorithm
//
// Complexity:
//   O(m F), where F is the maximum flow value.
// 
// Verified:
//   AOJ GRL_6_A: Maximum Flow
//
// Reference:
//   B. H. Korte and J. Vygen (2008):
//   Combinatorial Optimization: Theory and Algorithms.
//   Springer Berlin Heidelberg. 
//

#include <iostream>
#include <vector>
#include <cstdio>
#include <algorithm>
#include <functional>

using namespace std;

#define fst first
#define snd second
#define all(c) ((c).begin()), ((c).end())

const int INF = 1 << 30;
struct graph {
  int n;
  struct edge {
    int src, dst;
    int capacity, residue;
    size_t rev;
  };
  vector<vector<edge>> adj;
  graph(int n) : n(n), adj(n) { }
  void add_edge(int src, int dst, int capacity) {
    adj[src].push_back({src, dst, capacity, 0, adj[dst].size()});
    adj[dst].push_back({dst, src, 0, 0, adj[src].size()-1});
  }
  int max_flow(int s, int t) {
    vector<bool> visited(n);
    function<int (int,int)> augment = [&](int u, int f = INF) {
      if (u == t) return f;
      visited[u] = true;
      for (auto &e: adj[u]) {
        if (!visited[e.dst] && e.residue > 0) {
          int d = augment(e.dst, min(e.residue, f));
          if (d > 0) {
            e.residue -= d;
            adj[e.dst][e.rev].residue += d;
            return d;
          }
        }
      }
      return 0;
    };
    for (int u = 0; u < n; ++u)
      for (auto &e: adj[u]) e.residue = e.capacity;

    int total = 0;
    for (int f = 1; f; ) {
      fill(all(visited), false);
      total += (f = augment(s, INF));
    } // { u : visited[u] == true } is s-side
    return total;
  }
};

int main() {
  for (int n, m; scanf("%d %d", &n, &m) == 2; ) {
    graph g(n);
    for (int i = 0; i < m; ++i) {
      int u, v, w;
      scanf("%d %d %d", &u, &v, &w);
      g.add_edge(u, v, w);
    }
    printf("%d\n", g.max_flow(0, n-1));
  }
}

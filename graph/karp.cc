// 
// Minimum Mean Cycle (Karp)
//
//
// Description:
//   Given a directed graph G = (V, E) with edge length w.
//   Find a minimum mean cycle C, i.e., min w(C)/|C|.
//
//
// Algorithm:
//   Karp's algorithm. Fix a vertex s. 
//   By a dynamic programming, we can compute
//   the shortest path from s to v, with "exactly" k edges.
//   We write d(s,u;k) for this value.
//   Then, we can show that
//     min_u max_k [d(s,u;n) - d(s,u;k)]/(n-k)
//   is the length of minimum mean cycle.
//
//   (prf. Note that d(s,u;n) consists of cycle + path.
//    Subtract the path, we obtain a length of cycle.)
//
// 
// Complexity:
//   O(nm) time, O(n^2) space.
//
//
// Remark:
//   The algorithm can be applied to a "directed" network.
//   For an undirected network, MMC problem can be solved by
//   b-matching/T-join. See Korte and Vygen, Ch. 12.
//
//
// References
//   R. M. Karp (1978):
//   A characterization of the minimum cycle mean in a digraph.
//   Discrete Mathematics, vol. 23, pp. 309-311.
//
#include <iostream>
#include <vector>
#include <cstdio>
#include <cstdlib>
#include <map>
#include <cmath>
#include <cstring>
#include <functional>
#include <algorithm>
#include <unordered_map>

using namespace std;

typedef int weight_type;
const weight_type INF = 99999999;
struct edge {
  int src, dst;
  weight_type weight;
};
struct minimum_mean_cycle {
  vector<edge> edges;
  void add_edge(int src, int dst, weight_type weight) {
    // must be directed
    edges.push_back({src, dst, weight});
  }
  int n;
  vector<vector<edge>> adj;
  void make_graph(int n_ = 0) {
    n = n_;
    for (auto e: edges) 
      n = max(n, max(e.src, e.dst)+1);
    adj.resize(n);
    for (auto e: edges) 
      adj[e.src].push_back(e);
  }

  double solve() {
    vector<vector<weight_type>> dist(n+1, vector<weight_type>(n, INF));
    dist[0][0] = 0;
    for (int k = 0; k < n; ++k) 
      for (int u = 0; u < n; ++u) 
        for (auto e: adj[u]) 
          dist[k+1][e.dst] = min(dist[k+1][e.dst], dist[k][e.src] + e.weight);
    weight_type num = 1;
    int den = 0;
    for (int k = 0; k < n; ++k) 
      for (int u = 0; u < n; ++u) 
        if (dist[k][u] < INF)
          if (num * (n-k) > (dist[n][u]-dist[k][u]) * den) {
            num = (dist[n][u] - dist[k][u]);
            den = n-k;
          }
    return 1.0*num/den;
  }
};



int main() {
  minimum_mean_cycle solver;
  int n = 10;
  for (int i = 0; i < n; ++i) 
    for (int j = 0; j < n; ++j) 
      solver.add_edge(i, j, (rand() % (2*n)) - n);
  solver.make_graph(n);

  cout << solver.solve() << endl;
}

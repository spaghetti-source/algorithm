// 
// Minimum Mean Cycle (Karp)
//
// Description:
//   Given a directed graph G = (V, E) with edge length w.
//   Find a minimum mean cycle C, i.e., min w(C)/|C|.
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
//   (prf. Note that d(s,u;k) consists of cycle + path.
//    Subtract the path, we obtain a length of cycle.)
// 
// Complexity:
//   O(nm) time, O(n^2) space.
//
// Remark:
//   The algorithm can be applied to a "directed" network.
//   For an undirected network, MMC problem can be solved by
//   b-matching/T-join. See Korte and Vygen, Ch. 12.
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
#define fst first
#define snd second
#define all(c) ((c).begin()), ((c).end())

typedef int weight_type;
const weight_type INF = 99999999;
struct edge {
  int src, dst;
  weight_type weight;
};
struct graph {
  int n;
  vector<vector<edge>> adj;
  graph(int n = 0) : n(n) { }
  void add_edge(int src, int dst, weight_type weight) {
    // must be directed
    // negative edges are allowed
    n = max(n, max(src, dst)+1);
    adj.resize(n);
    adj[src].push_back({src, dst, weight}); 
  }
  pair<weight_type, weight_type> min_mean_cycle() {
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
    return {num, den}; // ratio = num/den
  }
};



int main() {
  graph solver;
  int n = 10;
  for (int i = 0; i < n; ++i) 
    for (int j = 0; j < n; ++j) 
      solver.add_edge(i, j, (rand() % (2*n)) - n);

  auto p = solver.min_mean_cycle();
  cout << p.fst << "/" << p.snd << endl;
}

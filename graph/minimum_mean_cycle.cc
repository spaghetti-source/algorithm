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

struct graph {
  typedef int weight_type;
  const weight_type INF = 99999999;
  struct edge {
    int src, dst;
    weight_type weight;
  };
  int n;
  vector<vector<edge>> adj;
  graph(int n) : n(n), adj(n) { }
  void add_edge(int src, int dst, weight_type weight) {
    adj[src].push_back({src, dst, weight}); 
  }
  typedef pair<weight_type, int> fraction;
  fraction min_mean_cycle() {
    vector<vector<weight_type>> dist(n+1, vector<weight_type>(n));
    vector<vector<int>> prev(n+1, vector<int>(n, -1));
    fill(all(prev[0]), 0);

    for (int k = 0; k < n; ++k) {
      for (int u = 0; u < n; ++u) {
        if (prev[k][u] < 0) continue;
        for (auto e: adj[u]) {
          if (prev[k+1][e.dst] < 0 ||
              dist[k+1][e.dst] > dist[k][e.src] + e.weight) {
            dist[k+1][e.dst] = dist[k][e.src] + e.weight;
            prev[k+1][e.dst] = e.src;
          }
        }
      }
    }
    int v = -1;
    fraction opt = {1, 0}; // +infty
    for (int u = 0; u < n; ++u) {
      fraction f = {-1, 0}; // -infty
      for (int k = n-1; k >= 0; --k) {
        if (prev[k][u] < 0) continue;
        fraction g = {dist[n][u] - dist[k][u], n - k};
        if (f.fst * g.snd < f.snd * g.fst) f = g;
      }
      if (opt.fst * f.snd > f.fst * opt.snd) { opt = f; v = u; }
    }
    if (v >= 0) { // found a loop
      vector<int> p; // path
      for (int k = n; p.size() < 2 || p[0] != p.back(); v = prev[k--][v]) 
        p.push_back(v);
      reverse(all(p));
    }
    return opt; 
  }
};



int main() {
  srand( 4 );
  int n = 4;
  graph solver(n);
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      if (i == j) continue;
      int c = (rand() % (2*n)) - n;
      solver.add_edge(i, j, c);
    }
  }

  auto p = solver.min_mean_cycle();
  cout << p.fst << "/" << p.snd << endl;
}

//
// Betweenness centrality of undirected unweighted graph (Brandes)
//
// Description:
// 
//   Compute betweenness centrality, defined by
//     f(u) := \sum_{u,t \eq v} |s-t shortest paths that contains v|/|s-t shortest paths|
//
// Algorithm:
//
//   Brandes's algorithm, O(nm) time, O(m) space.
//
// References:
//
//   U. Brandes (2001): A faster algorithm for betweenness centrality.
//   Journal of Mathematical Sociology, vol.25, pp.163–177.

#include <iostream>
#include <vector>
#include <cstdio>
#include <cstdlib>
#include <queue>
#include <map>
#include <cmath>
#include <cstring>
#include <functional>
#include <algorithm>
#include <unordered_map>
#include <unordered_set>

using namespace std;

#define fst first
#define snd second
#define all(c) ((c).begin()), ((c).end())

struct edge {
  size_t src, dst;
};
struct graph {
  vector<edge> edges;
  void add_edge(size_t src, size_t dst) {
    edges.push_back({src, dst});
  }
  size_t n;
  vector<vector<edge>> adj;
  void make_graph(int n_ = 0) {
    n = n_;
    for (auto e: edges) 
      n = max(n, max(e.src, e.dst)+1);
    adj.resize(n);
    for (auto e: edges) {
      adj[e.src].push_back(e);
      swap(e.src, e.dst);
      adj[e.src].push_back(e);
    }
  }

  vector<double> betweeness_centrality() {
    vector<double> centrality(n);

    for (size_t s = 0; s < n; ++s) {
      vector<size_t> S;
      vector<double> sigma(n); sigma[s] = 1;
      vector<int> dist(n, -1); dist[s]  = 0;
      queue<size_t> que;       que.push(s);
      while (!que.empty()) {
        size_t u = que.front();
        S.push_back(u);
        que.pop();
        for (auto e: adj[u]) {
          if (dist[e.dst] < 0) {
            dist[e.dst] = dist[e.src] + 1;
            que.push(e.dst);
          }
          if (dist[e.dst] == dist[e.src] + 1) {
            sigma[e.dst] += sigma[e.src];
          }
        }
      }
      vector<double> delta(n);
      while (!S.empty()) {
        size_t u = S.back();
        S.pop_back();
        for (auto e: adj[u]) {
          if (dist[e.dst] == dist[e.src] + 1) {
            delta[e.src] += sigma[e.src] / sigma[e.dst] * (1 + delta[e.dst]);
          }
        }
        if (u != s) centrality[u] += delta[u];
      }
    }
    return centrality;
  }
};

int main() {
  graph g;
  g.add_edge( 0, 1 );
  g.add_edge( 0, 2 );
  g.add_edge( 0, 3 );

  g.add_edge( 1, 4 );
  g.add_edge( 1, 5 );
  g.add_edge( 4, 5 );
  g.add_edge( 4, 6 );
  g.add_edge( 5, 6 );

  g.add_edge( 2, 7 );
  g.add_edge( 2, 8 );
  g.add_edge( 7, 8 );
  g.add_edge( 7, 9 );
  g.add_edge( 8, 9 );

  g.add_edge( 3,10 );
  g.add_edge( 3,11 );
  g.add_edge(10,11 );
  g.add_edge(10,12 );
  g.add_edge(11,12 );

  g.make_graph();
  g.betweeness_centrality();
}


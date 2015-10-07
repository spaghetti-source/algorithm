// 
// Undirected Eulerian Path (Hierholzer's algorithm)
// 
// Description:
//   A path (walk) P is Eulerian if P visits all edges exactly once.
//   A graph admits Eulerian path if and only if there are exactly two
//   odd-degree vertices (they are initial and terminal of a path).
//
// Algorithm:
//   Hierholzer's algorithm performs DFS from the initial vertex.
//
// Complexity:
//   O(n + m)
//
// Verified:
//   UVA 10054 The Necklace (tour)
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

struct graph {
  int n;
  struct edge { int src, dst, rev; };
  vector<vector<edge>> adj;
  graph() : n(0) { }
  //graph(int n) : n(n), adj(n) { }

  void add_edge(int src, int dst) {
    n = max(n, max(src, dst)+1);
    if (adj.size() < n) adj.resize(n);
    adj[src].push_back({src, dst, (int)adj[dst].size()});
    adj[dst].push_back({dst, src, (int)adj[src].size()-1});
  }

  // destructive
  vector<int> path;
  void visit(int u) {
    while (!adj[u].empty()) {
      auto e = adj[u].back();
      adj[u].pop_back();
      if (e.src >= 0) {
        adj[e.dst][e.rev].src = -1;
        visit(e.dst);
      }
    }
    path.push_back(u);
  }
  vector<int> eulerian_path() {
    int m = 0, s = -1;
    for (int u = 0; u < n; ++u) {
      m += adj[u].size();
      if (adj[u].size() % 2 == 1) s = u;
    }
    path.clear(); if (s >= 0) visit(s);
    if (path.size() != m/2 + 1) return {};
    return path;
  }
  vector<int> eulerian_tour() {
    int m = 0, s = 0;
    for (int u = 0; u < n; ++u) {
      m += adj[u].size();
      if (adj[u].size() > 0) s = u;
    }
    path.clear(); visit(s);
    if (path.size() != m/2 + 1 || path[0] != path.back()) return {};
    return path;
  }
};


int main() {
  int ncase;
  scanf("%d", &ncase);
  for (int icase = 0; icase < ncase; ++icase) {
    if (icase > 0) printf("\n");
    printf("Case #%d\n", icase+1);

    int m;
    scanf("%d", &m);
    graph g;
    for (int i = 0; i < m; ++i) {
      int s, t;
      scanf("%d %d", &s, &t);
      g.add_edge(s-1, t-1);
    }
    auto path = g.eulerian_tour();
    if (path.empty()) {
      printf("some beads may be lost\n");
    } else {
      for (int i = 0; i+1 < path.size(); ++i) 
        printf("%d %d\n", path[i]+1, path[i+1]+1);
    }
  }
}

// 
// Transitive Reduction of DAG
//
// Description:
//   A transitive reduction of a graph G = (V, E) is a
//   a graph H = (V, F) such that transitive closures of 
//   H and G are the same. 
//   There are possibly many transitive reductions with 
//   the fewest edges, and finding one of them is NP-hard.
//   On he other hand, if a graph is directed acyclic, 
//   its transitive reduction uniquely exists and can be
//   found in O(nm) time.
//   Note that transitive closure and reduction have the 
//   same time complexity on DAG.
//
// Algorithm:
//   For each vertex u, compute longest path distance from u
//   to v in adj[u]. Then, remove all edges (u,v) with d(u,v) > 1.
//
// Complexity:
//   O(nm). Usually the coefficient is not so large.
//
#include <iostream>
#include <vector>
#include <cstdio>
#include <set>
#include <algorithm>
#include <functional>

using namespace std;

#define fst first
#define snd second
#define all(c) ((c).begin()), ((c).end())
#define TEST(s) if (!(s)) { cout << __LINE__ << " " << #s << endl; exit(-1); }

struct graph { // DAG
  int n;
  vector<vector<int>> adj;
  graph(int n) : n(n), adj(n) { }
  void add_edge(int i, int j) { 
    adj[i].push_back(j);
  }

  void transitive_reduction() {
    vector<int> ord, d(n, -1);
    function<void(int)> rec = [&](int u) {
      d[u] = 0;
      for (int v: adj[u]) 
        if (d[v] < 0) rec(v);
      ord.push_back(u); 
      for (int v: ord) d[v] = 0;
      for (int i = ord.size()-1; i >= 0; --i) 
        for (int w: adj[ord[i]]) 
          d[w] = max(d[w], d[ord[i]] + 1);

      adj[u].erase(
          remove_if(all(adj[u]), [&](int v) { return d[v] > 1; }),
          adj[u].end()
      );
    };
    for (int u = 0; u < n; ++u)
      if (d[u] < 0) rec(u);
  }
};


// === tick a time ===
#include <ctime>
double tick() {
  static clock_t oldtick;
  clock_t newtick = clock();
  double diff = 1.0*(newtick - oldtick) / CLOCKS_PER_SEC;
  oldtick = newtick;
  return diff;
}

int main() {
  int n = 10000, m = 100 * n;
  set<pair<int, int>> edges;
  while (edges.size() < m) {
    int u = rand() % n, v = rand() % n;
    if (u == v) continue;
    if (u > v) swap(u, v);
    edges.insert({u, v});
  }
  graph g(n);
  for (auto p: edges) 
    g.add_edge(p.fst, p.snd);


  tick();
  g.transitive_reduction();
  cout << tick() << endl;

  int mm = 0;
  for (int u = 0; u < n; ++u)
    mm += g.adj[u].size();
  cout << m << " ==> " << mm << " ( " << mm/double(m) << ")" << endl;
}

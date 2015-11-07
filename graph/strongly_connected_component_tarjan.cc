// 
// Tarjan's strongly connected component
//
// Description:
//   For a graph G = (V, E), u and v are strongly connected if
//   there are paths u -> v and v -> u. This defines an equivalent
//   relation, and its equivalent class is called a strongly 
//   connected component.
//
// Algorithm:
//   Tarjan's single stack algorithm. 
//   This does not need to keep reverse edges;
//   thus it is memory efficient than Kosaraju's algorithm.
//
// Complexity:
//   O(n + m)
//
// Verified:
//   SPOJ 6818
//
// References: 
//   R. E. Tarjan (1972):
//   Depth-first search and linear graph algorithms.
//   SIAM Journal on Computing, vol.1, no.2, pp.146--160.
//
#include <iostream>
#include <vector>
#include <cstdio>
#include <cstdlib>
#include <map>
#include <set>
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

struct graph {
  int n;
  vector<vector<int>> adj;
  graph(int n) : n(n), adj(n) { }
  void add_edge(int src, int dst) {
    adj[src].push_back(dst);
  }

  vector<vector<int>> strongly_connected_components() {
    vector<int> open, id(n);
    vector<vector<int>> scc;
    int t = -n-1;
    auto argmin = [&](int u, int v) { return id[u] < id[v] ? u : v; };
    function<int(int)> dfs = [&](int u) {
      open.push_back(u);
      id[u] = t++;
      int w = u;
      for (int v: adj[u]) {
        if      (id[v] == 0) w = argmin(w, dfs(v));
        else if (id[v]  < 0) w = argmin(w, v);
      }
      if (w == u) {
        scc.push_back({});
        while (1) {
          int v = open.back();
          open.pop_back();
          id[v] = scc.size();
          scc.back().push_back(v);
          if (u == v) break;
        }
      }
      return w;
    };
    for (int u = 0; u < n; ++u) 
      if (id[u] == 0) dfs(u);
    return scc;
  }
};

int main() {
  int n, m;
  scanf("%d %d", &n, &m);
  graph g(n);
  for (int k = 0; k < m; ++k) {
    int i, j;
    scanf("%d %d", &i, &j);
    g.add_edge(i-1, j-1);
  }

  vector<vector<int>> scc = g.strongly_connected_components();
  vector<int> outdeg(scc.size());
  vector<int> id(n);
  for (int i = 0; i < scc.size(); ++i)
    for (int u: scc[i]) id[u] = i;
  for (int u = 0; u < n; ++u) 
    for (int v: g.adj[u]) 
      if (id[u] != id[v]) ++outdeg[id[u]];

  if (count(all(outdeg), 0) != 1) {
    printf("0\n");
  } else {
    int i = find(all(outdeg), 0) - outdeg.begin();
    sort(all(scc[i]));
    printf("%d\n%d", scc[i].size(), scc[i][0]+1);
    for (int j = 1; j < scc[i].size(); ++j) 
      printf(" %d", scc[i][j]+1);
    printf("\n");
  }
}

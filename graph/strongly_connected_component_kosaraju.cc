// 
// Kosaraju's strongly connected component
//
// Description:
//   For a graph G = (V, E), u and v are strongly connected if
//   there are paths u -> v and v -> u. This defines an equivalent
//   relation, and its equivalent class is called a strongly 
//   connected component.
//
// Algorithm:
//   Kosaraju's algorithm performs DFS on G and rev(G).
//   First DFS finds topological ordering of SCCs, and
//   the second DFS extracts components.
//
// Complexity:
//   O(n + m)
//
// Verified:
//   SPOJ 6818
//
// References:
//   A. V. Aho, J. E. Hopcroft, and J. D. Ullman (1983):
//   Data Structures and Algorithms,
//   Addison-Wesley.
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
  vector<vector<int>> adj, rdj;
  graph(int n) : n(n), adj(n), rdj(n) { }
  void add_edge(int src, int dst) {
    adj[src].push_back(dst);
    rdj[dst].push_back(src);
  }

  vector<vector<int>> strongly_connected_components() { // kosaraju
    vector<int> ord, visited(n);
    vector<vector<int>> scc;
    function<void(int,vector<vector<int>>&, vector<int>&)> dfs 
      = [&](int u, vector<vector<int>> &adj, vector<int> &out) {
      visited[u] = true;
      for (int v: adj[u]) 
        if (!visited[v]) dfs(v, adj, out);
      out.push_back(u);
    };
    for (int u = 0; u < n; ++u)
      if (!visited[u]) dfs(u, adj, ord);
    fill(all(visited), false);
    for (int i = n-1; i >= 0; --i) 
      if (!visited[ord[i]]) 
        scc.push_back({}), dfs(ord[i], rdj, scc.back()); 
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

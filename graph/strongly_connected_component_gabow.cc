//
// Gabow's strongly connected component
//
// Description:
//   For a graph G = (V, E), u and v are strongly connected if
//   there are paths u -> v and v -> u. This defines an equivalent
//   relation, and its equivalent class is called a strongly 
//   connected component.
//
//
// Algorithm:
//   Gabow's double stack algorithm.
//   This is simpler than Tarjan's algorithm.
//   In my opinion, use Gabow or Kosaraju for SCC problems.
//
// Complexity:
//   O(n + m)
//
// Verified:
//   SPOJ 6818
//
// References: 
//   H. N. Gabow (2000):
//   Path-based depth first search strong and biconnected components.
//   Information Processing Letters, vol.74 no.3-4, pp.107-114.
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
    vector<vector<int>> scc;
    vector<int> S, B, I(n);
    function<void(int)> dfs = [&](int u) {
      B.push_back(I[u] = S.size());
      S.push_back(u);
      for (int v: adj[u]) {
        if (!I[v]) dfs(v);
        else while (I[v] < B.back()) B.pop_back();
      }
      if (I[u] == B.back()) {
        scc.push_back({});
        B.pop_back();
        for (; I[u] < S.size(); S.pop_back()) {
          scc.back().push_back(S.back());
          I[S.back()] = n + scc.size();
        }
      }
    };
    for (int u = 0; u < n; ++u)
      if (!I[u]) dfs(u);
    return scc; // I[u] - n is the index of u
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

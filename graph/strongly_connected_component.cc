//
// Strongly Connected Components
//
// Description:
//   Compute the strongly connected components(SCCs) decomposition.
//   A set of vertices S of a digraph is a SCC iff 
//   for each u, v in S, there is a path from u to v.
//   Any graph can be uniquely decomposed into SCCs.
//
// Algorithm:
//   Tarjan's single DFS / single stack algorithm.
//
// Complexity:
//   O(n + m) time and space. 
//   Note that this is much faster than Kosaraju's two DFS algorithm,
//   and almost as fast as Gabow's single DFS / two stacks algorithm.
//
// Verified:
//   SPOJ6818
//
// References: 
// - R. E. Tarjan (1972):
//   Depth-first search and linear graph algorithms.
//   SIAM Journal on Computing, vol.1, no.2, pp.146–160.
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

using namespace std;

#define fst first
#define snd second
#define all(c) ((c).begin()), ((c).end())

struct edge {
  int src, dst;
};
struct graph {
  int n;
  vector<vector<edge>> adj;
  void add_edge(int src, int dst) {
    adj[src].push_back({src, dst});
  }
  graph(int n) : n(n), adj(n) { }

  vector<int> open, id;
  vector<vector<int>> scc;
  int visit(int i, int &t) {
    int w = i;
    open.push_back(i);
    id[i] = t++;
    for (auto e: adj[i]) {
      if (id[e.dst] == 0) {
        int k = visit(e.dst, t);
        if (id[w] > id[k]) w = k;
      } else if (id[e.dst] < 0) {
        if (id[w] > id[e.dst]) w = e.dst;
      }
    }
    if (w == i) {
      scc.push_back({});
      while (1) {
        int j = open.back();
        open.pop_back();
        id[j] = scc.size();
        scc.back().push_back(j);
        if (i == j) break;
      }
    }
    return w;
  }
  int strongly_connected_components() {
    scc.clear();
    open.clear();
    id.assign(n, 0);
    for (int i = 0, t = -n-1; i < n; ++i) {
      if (id[i] == 0) {
        visit(i, t);
      }
    }
    return scc.size();
  }

  void solve() {
    strongly_connected_components();
    vector<int> outdeg(scc.size());
    for (int i = 0; i < n; ++i) {
      for (auto e: adj[i]) {
        if (id[e.src] != id[e.dst]) ++outdeg[ id[e.src]-1 ];
      }
    }
    int c = find(all(outdeg), 0) - outdeg.begin();
    if (c == outdeg.size()) {
      printf("0\n");
      return;
    }
    vector<int> component;
    for (int i = 0; i < n; ++i) {
      if (id[i]-1 == c) component.push_back(i);
    }
    printf("%d\n%d", component.size(), component[0]+1);
    for (int i = 1; i < component.size(); ++i) 
      printf(" %d", component[i]+1);
    printf("\n");
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
  g.solve();
}

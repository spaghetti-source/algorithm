//
// Dominator Tree (Lengauer-Tarjan)
//
// Description:
//   Let G = (V, E) be a directed graph and fix r in V.
//   v is a dominator of u if any paths from r to u through v.
//   The set of dominators of u forms a total order, and 
//   the closest dominator is called the immediate dominator.
//   The set { (u,v) : v is the immediate dominator } forms a tree, 
//   which is called the dominator tree.
//
// Algorithm:
//   Lengauer and Tarjan's DFS-based algorithm.
//
// Complexity:
//   O(m log n)
//
// Verified:
//   SPOJ964
// 
// References:
//   T. Lengauer and R. Tarjan (1979):
//   A fast algorithm for finding dominators in a flowgraph.
//   ACM Transactions on Programming Languages and Systems,
//   vol.1, no.1, pp.121-141.
//


#include <iostream>
#include <vector>
#include <cstdio>
#include <cstdlib>
#include <functional>
#include <algorithm>

using namespace std;

#define fst first
#define snd second
#define all(c) ((c).begin()), ((c).end())

struct edge { int src, dst; };
struct graph {
  int n;
  vector<vector<edge>> adj, rdj;
  graph(int n = 0) : n(0) { } 
  void add_edge(int src, int dst) {
    n = max(n, max(src, dst)+1);
    adj.resize(n); rdj.resize(n);
    adj[src].push_back({src, dst});
    rdj[dst].push_back({dst, src});
  }

  vector<int> rank, semi, low, anc;
  int eval(int v) { 
    if (anc[v] < n && anc[anc[v]] < n) {
      int x = eval(anc[v]);
      if (rank[semi[low[v]]] > rank[semi[x]]) low[v] = x;
      anc[v] = anc[anc[v]];
    }
    return low[v];
  }
  vector<int> prev, ord;
  void dfs(int u) {
    rank[u] = ord.size();
    ord.push_back(u);
    for (auto e: adj[u]) {
      if (rank[e.dst] < n) continue;
      dfs(e.dst);
      prev[e.dst] = u;
    }
  }
  vector<int> idom; // idom[u] is an immediate dominator of u
  void dominator_tree(int r) {
    idom.assign(n, n); prev = rank = anc = idom;
    semi.resize(n); iota(all(semi), 0); low = semi;
    ord.clear(); dfs(r);

    vector<vector<int>> dom(n);
    for (int i = ord.size()-1; i >= 1; --i) {
      int w = ord[i];
      for (auto e: rdj[w]) {
        int u = eval(e.dst);
        if (rank[semi[w]] > rank[semi[u]]) semi[w] = semi[u];
      }
      dom[semi[w]].push_back(w);
      anc[w] = prev[w];
      for (int v: dom[prev[w]]) {
        int u = eval(v);
        idom[v] = (rank[prev[w]] > rank[semi[u]] ? u : prev[w]);
      }
      dom[prev[w]].clear();
    }
    for (int i = 1; i < ord.size(); ++i) {
      int w = ord[i];
      if (idom[w] != semi[w]) idom[w] = idom[idom[w]];
    }
  }
  vector<int> dominators(int u) {
    vector<int> S;
    for (; u < n; u = idom[u]) S.push_back(u);
    return S;
  }
};

int main() {
  int t; scanf("%d", &t);
  while (t--) {
    int n, m; scanf("%d %d", &n, &m);
    graph g(n);
    for (int i = 0; i < m; ++i) {
      int u, v; scanf("%d %d", &u, &v);
      g.add_edge(u-1, v-1);
    }
    g.dominator_tree(0);
    auto S = g.dominators(n-1);
    cout << S[S.size()-2]+1 << endl;
  }
}

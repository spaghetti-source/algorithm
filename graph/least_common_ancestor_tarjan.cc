//
// Offline least common ancestor
//
// Description
//   For a rooted tree T, LCA(u,v) is a vertex u
//   that is the deepest node that is a common ancestor of u and v.
//   It computes all lcas of (u_j, v_j) for v = 1, ..., q. 
//
// Algorithm
//   Tarjan's dfs and union-find.
//
// Complexity:
//   O((m+q) a(n)), where a(n) is the inverse Ackermann function.
//
// Verified:
//   SPOJ14932

#include <iostream>
#include <vector>
#include <cstdio>
#include <unordered_map>
#include <algorithm>

using namespace std;

#define fst first
#define snd second
#define all(c) ((c).begin()), ((c).end())

struct graph {
  int n;
  vector<vector<int>> adj;
  graph(int n = 0) : n(n), adj(n) { }
  void add_edge(int src, int dst) {
    n = max(n, max(src, dst)+1);
    adj.resize(n);
    adj[src].push_back(dst);
  }
  struct query { int u, v, a; };
  struct union_find {
    vector<int> p; 
    union_find(int n) : p(n, -1) { };
    bool unite(int u, int v) { 
      if ((u = root(u)) == (v = root(v))) return false;
      if (p[u] > p[v]) swap(u, v);
      p[u] += p[v]; p[v] = u;
      return true;
    }
    int root(int u) { return p[u] < 0 ? u : p[u] = root(p[u]); }
  };
  void lca(vector<query> &queries) {
    vector<vector<query*>> Q(n);
    for (auto &q: queries) {
      Q[q.u].push_back(&q);
      Q[q.v].push_back(&q);
    }
    union_find uf(n);
    vector<int> anc(n), color(n);
    iota(all(anc), 0);
    function<void (int)> rec = [&](int u) {
      for (auto v: adj[u]) {
        rec(v);
        uf.unite(u, v);
        anc[uf.root(u)] = u;
      }
      color[u] = 1;
      for (auto it: Q[u]) {
        if (it->u != u) swap(it->u, it->v);
        if (color[it->v] == 1) it->a = anc[uf.root(it->v)];
      }
    };
    vector<int> deg(n);
    for (int u = 0; u < n; ++u)
      for (auto v: adj[u])
        ++deg[v];
    for (int u = 0; u < n; ++u)
      if (deg[u] == 0) rec(u);
  }
};

int main() {
  int ncase; scanf("%d", &ncase);
  for (int icase = 0; icase < ncase; ++icase) {
    printf("Case %d:\n", icase+1);
    int n; scanf("%d", &n);
    graph g(n);
    for (int i = 0; i < n; ++i) {
      int k; scanf("%d", &k);
      for (int j = 0; j < k; ++j) {
        int l; scanf("%d", &l);
        g.add_edge(i, l-1);
      }
    }
    int q; scanf("%d", &q);
    vector<graph::query> queries;
    for (int i = 0; i < q; ++i) {
      int u, v; scanf("%d %d", &u, &v);
      queries.push_back({u-1, v-1, -1});
    }
    g.lca(queries);
    for (auto q: queries) 
      printf("%d\n", q.a+1);
  }
}

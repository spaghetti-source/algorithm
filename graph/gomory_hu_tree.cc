//
// Gomory-Hu tree
//
// Description:
//   A cut tree T = (V, F) of a undirected graph G = (V, E) 
//   is a tree such that
//     mincut_T(u,v) = mincut_G(u,v)  for all u,v in V.
//   Such tree always exists, and is called a Gomory-Hu tree.
//  
// Algorithm:
//   Gusfield's simplified version of Gomory-Hu algorithm.
//   
// Complexity:
//   O(n-1) max-flow call.
//
// Verified:
//   SPOJ3900
//
// References:
//   R. E. Gomory and T. C. Hu (1961):
//   Multi-terminal network flows. 
//   Journal of the Society for Industrial and Applied Mathematics,
//   vol. 9, 1961.
//
//   D. Gusfield (1990):
//   Very Simple Methods for All Pairs Network Flow Analysis.
//   SIAM Journal of Computing, vol.19, no.1, pp.143-155.
//
#include <iostream>
#include <cstdio>
#include <queue>
#include <map>
#include <algorithm>

using namespace std;
#define fst first
#define snd second
#define all(c) ((c).begin()), ((c).end())

typedef long long flow_type;
flow_type INF = (1ULL<<50);
struct edge {
  int src, dst;
  flow_type capacity;
  int rev;
  flow_type residue;
};
struct graph {
  int n;
  vector<vector<edge>> adj;
  graph(int n = 0) : n(n), adj(n) { }
  void add_edge(int src, int dst, flow_type capacity) {
    n = max(n, max(src, dst)+1);
    adj.resize(n);
    adj[src].push_back({src, dst, capacity, (int)adj[dst].size()});
    adj[dst].push_back({dst, src, capacity, (int)adj[src].size()-1});
  }
  vector<int> level, iter;
  flow_type augment(int u, int t, flow_type cur) {
    if (u == t) return cur;
    for (int &i = iter[u]; i < adj[u].size(); ++i) {
      edge &e = adj[u][i];
      if (e.residue > 0 && level[u] < level[e.dst]) {
        flow_type f = augment(e.dst, t, min(cur, e.residue));
        if (f > 0) {
          e.residue -= f;
          adj[e.dst][e.rev].residue += f;
          return f;
        }
      }
    }
    return 0;
  }
  int bfs(int s, int t) {
    level.assign(n, -1); level[s] = 0;
    queue<int> Q; Q.push(s);
    while (!Q.empty()) {
      int u = Q.front(); Q.pop();
      if (u == t) break;
      for (auto &e: adj[u]) {
        if (e.residue > 0 && level[e.dst] < 0) {
          Q.push(e.dst);
          level[e.dst] = level[u] + 1;
        }
      }
    }
    return level[t];
  }
  flow_type max_flow(int s, int t) {
    for (int u = 0; u < n; ++u) // initialize
      for (auto &e: adj[u]) e.residue = e.capacity;
    flow_type flow = 0;
    int itera = 0;
    while (bfs(s, t) >= 0) {
      iter.assign(n, 0);
      for (flow_type f; (f = augment(s, t, INF)) > 0; )
        flow += f;
    } // level[u] == -1 ==> t-side
    return flow;
  }
  vector<edge> tree;
  void gomory_hu() {
    tree.clear();
    vector<int> parent(n);
    for (int u = 1; u < n; ++u) {
      tree.push_back({u, parent[u], max_flow(u, parent[u])});
      for (int v = u+1; v < n; ++v)
        if (level[v] >= 0 && parent[v] == parent[u])
          parent[v] = u;
    }
  }
};



struct union_find {
  vector<int> p; 
  union_find(int n) : p(n, -1) { };
  bool unite(int u, int v) { 
    if ((u = root(u)) == (v = root(v))) return false;
    if (p[u] > p[v]) swap(u, v);
    p[u] += p[v]; p[v] = u;
    return true;
  }
  bool find(int u, int v) { return root(u) == root(v); }
  int root(int u) { return p[u] < 0 ? u : p[u] = root(p[u]); }
  int size(int u) { return -p[root(u)]; }
};

int main() {
  int ncase;
  scanf("%d", &ncase);
  for (int icase = 0; icase < ncase; ++icase) {
    if (icase > 0) printf("\n");
    int n, m;
    scanf("%d %d", &n, &m);
    graph g(n);
    for (int i = 0; i < m; ++i) {
      int u, v, c; scanf("%d %d %d", &u, &v, &c);
      g.add_edge(u-1, v-1, c);
    }
    g.gomory_hu();

    map<flow_type, vector<edge>> edges;
    for (auto e: g.tree)
      edges[-e.capacity].push_back(e);

    union_find uf(n);
    vector<flow_type> fs;
    vector<int> cs;
    int sum = n*(n-1)/2;
    fs.push_back(INF);
    cs.push_back(sum);
    for (auto &p: edges) {
      fs.push_back(-p.fst);
      cs.push_back(sum);
      for (auto e: p.snd) {
        if (!uf.find(e.src, e.dst)) {
          sum -= uf.size(e.src) * uf.size(e.dst);
          uf.unite(e.src, e.dst);
        }
      }
    }
    fs.push_back(-1);
    cs.push_back(0);
    reverse(all(fs));
    reverse(all(cs));
    int q; scanf("%d", &q);
    for (int i = 0; i < q; ++i) {
      int a; scanf("%d", &a);
      int k = distance(fs.begin(), lower_bound(all(fs), a));
      while (fs[k] > a) --k;
      printf("%d\n", cs[k]);
    }
  }
}

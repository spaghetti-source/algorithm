//
// Dominator Tree (Lengauer-Tarjan)
//
// Description:
//   Let G = (V, E) be a directed graph and fix r in V.
//   v is a dominator of u if any path between r and u throughs v.
//   The set of dominators of u forms a total order, and 
//   the minimum dominator is called the immediate dominator.
//   The dominator tree is { (u,v) : v is the immediate dominator }.
//
// Algorithm:
//   Lengauer and Tarjan.
//
// Complexity:
//   O(m log n)
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
  void add_edge(int src, int dst) {
    n = max(n, max(src, dst)+1);
    adj.resize(n); rdj.resize(n);
    adj[src].push_back({src, dst});
    rdj[dst].push_back({dst, src});
  }

  vector<int> rank, semi, label, ancestor;
  void compress(int v) {
    int a = ancestor[v];
    if (ancestor[a] < 0) return;
    compress(a);
    if (rank[semi[label[v]]] > rank[semi[label[a]]]) label[v] = label[a];
    ancestor[v] = ancestor[a];
  }
  int eval(int v) {
    if (ancestor[v] < 0) return v;
    compress(v); 
    return label[v];
  }
  vector<int> parent, ord;
  void dfs(int u) {
    rank[u] = ord.size();
    ord.push_back(u);
    for (auto e: adj[u]) {
      if (rank[e.dst] >= 0) continue;
      parent[e.dst] = u;
      dfs(e.dst);
    }
  }
  vector<int> dom;
  void dominator_tree(int r) {
    rank.assign(n, -1); 
    parent.resize(n);
    ord.clear();
    dfs(r);

    ancestor.assign(n, -1);
    semi.resize(n); iota(all(semi), 0);
    label = semi;

    dom.assign(n, -1);
    vector<vector<int>> bucket(n);
    for (int i = ord.size()-1; i >= 1; --i) {
      int w = ord[i];
      for (auto e: rdj[w]) {
        int u = eval(e.dst);
        if (rank[semi[u]] < rank[semi[w]]) semi[w] = semi[u];
      }
      bucket[semi[w]].push_back(w);
      ancestor[w] = parent[w];
      for (int v: bucket[parent[w]]) {
        int u = eval(v);
        dom[v] = (rank[semi[u]] < rank[semi[v]] ? u : parent[w]);
      }
      bucket[parent[w]].clear();
    }
    for (int i = 1; i < ord.size(); ++i) {
      int w = ord[i];
      if (dom[w] != semi[w]) dom[w] = dom[dom[w]];
    } // i's immediate dominator = dom[i]
  }
  vector<int> dominators(int u) {
    vector<int> S;
    for (; u >= 0; u = dom[u]) S.push_back(u);
    return S;
  }
};


int main() {
  graph g;
  g.add_edge(0, 1);
  g.add_edge(1, 2);
  g.add_edge(2, 3);
  g.add_edge(3, 4);
  g.add_edge(4, 5);
  g.add_edge(0, 4);

  g.dominator_tree(0);
  for (auto a: g.dominators(5))
    cout << a << endl;
}

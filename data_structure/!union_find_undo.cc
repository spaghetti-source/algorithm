//
// Undoable Union Find
//
// It can undo each unite operation. 
//
// Naming:
//   Data structures with undo operation is weaker than
//   partially-persistence (accessible all older vers) and
//   slightly stronger than semi-persistence (backtrackable).
//
#include <iostream>
#include <set>
#include <map>
#include <vector>
#include <algorithm>
#include <iterator>

using namespace std;
#define fst first
#define snd second
#define all(c) begin(c), end(c)


struct UndoableUnionFind {
  vector<int> parent; 
  vector<tuple<int,int,int>> history;
  UndoableUnionFind(int n) : parent(n, -1) { };
  bool unite(int u, int v) { 
    u = root(u); v = root(v);
    if (u == v) return false;
    if (parent[u] > parent[v]) swap(u, v);
    history.push_back(make_tuple(u, v, parent[v]));
    parent[u] += parent[v]; parent[v] = u;
    return true;
  }
  void undo() {
    int u, v, w; 
    tie(u, v, w) = history.back();
    history.pop_back();
    parent[v] = w;
    parent[u] -= parent[v];
  }
  bool find(int u, int v) { return root(u) == root(v); }
  int root(int u) { return parent[u] < 0 ? u : parent[u] = root(parent[u]); }
  int size(int u) { return -parent[root(u)]; }
};


struct OfflineDynamicConnectivity {
  int n;
  UndoableUnionFind uf;
  OfflineDynamicConnectivity(int n) : n(n), uf(n) { }

  typedef pair<int,int> Edge;
  vector<Edge> query;
  vector<int> ans;
  map<Edge, vector<pair<int,int>>> appear;

  void addEdge(int u, int v) {
    if (u > v) swap(u, v);
    appear[{u,v}].push_back({query.size(), 1<<30});
  }
  void eraseEdge(int u, int v) {
    if (u > v) swap(u, v);
    appear[{u,v}].back().snd = query.size();
  }
  int isConnected(int u, int v) {
    query.push_back({u,v});
    return query.size()-1;
  }
  vector<set<Edge>> edges;
  void insert(int l, int r, int k, int s, int t, Edge e) {
    s = max(s, l); t = min(r, t);
    if (s >= t) return;
    if (s == l && t == r) {
      edges[k].insert(e);
    } else {
      insert(l, (l+r)/2, 2*k+1, s, t, e);
      insert((l+r)/2, r, 2*k+2, s, t, e);
      if (edges[2*k+1].count(e) && edges[2*k+2].count(e)) {
        edges[2*k+1].erase(e);
        edges[2*k+2].erase(e);
        edges[k].insert(e);
      }
    }
  }
  void rec(int l, int r, int k) {
    if (l >= r) return;
    for (Edge e: edges[k]) uf.unite(e.fst, e.snd);
    if (l+1 == r) {
      ans[l] = uf.find(query[l].fst, query[l].snd);
    } else {
      rec(l, (l+r)/2, 2*k+1);
      rec((l+r)/2, r, 2*k+2);
    }
    for (Edge e: edges[k]) uf.undo();
  }
  void solve() {
    int q = query.size();
    edges.resize(4*q);
    for (auto a: appear) 
      for (auto b: a.snd) 
        insert(0, q, 0, b.fst, b.snd, a.fst);
    ans.assign(q, 0);
    rec(0, q, 0);
  }
};

int main() {
  OfflineDynamicConnectivity solver(3);
  solver.isConnected(0,1);
  solver.addEdge(0,1);
  solver.isConnected(0,1);
  solver.isConnected(1,2);
  solver.addEdge(1,2);
  solver.isConnected(0,2);
  solver.eraseEdge(0,1);
  solver.isConnected(0,2);
  solver.solve();
  for (int i = 0; i < solver.query.size(); ++i) {
    cout << solver.query[i].fst << " " << solver.query[i].snd << " " << solver.ans[i] << endl;
  }
}

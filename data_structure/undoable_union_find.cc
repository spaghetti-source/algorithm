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

struct undoable_union_find {
  vector<int> p;
  undoable_union_find(int n) : p(n, -1) { };
  vector<tuple<int,int,int>> hist;
  bool unite(int u, int v) { 
    if ((u = root(u)) == (v = root(v))) return false;
    if (p[u] > p[v]) swap(u, v);
    hist.push_back(make_tuple(u, v, p[v]));
    p[u] += p[v]; p[v] = u;
    return true;
  }
  void undo() {
    int u, v, w; tie(u, v, w) = hist.back();
    hist.pop_back();
    p[v] = w;
    p[u] -= p[v];
  }
  bool find(int u, int v) { return root(u) == root(v); }
  int root(int u) { for (; p[u] >= 0; u = p[u]); return u; }
  int size(int u) { return -p[root(u)]; }
};

struct offline_dynamic_connectivity {
  int n;
  undoable_union_find uf;
  offline_dynamic_connectivity(int n) : n(n), uf(n) { }

  typedef pair<int,int> edge;
  vector<edge> query;
  vector<int> ans;
  map<edge, vector<pair<int,int>>> appear;

  void add_edge(int u, int v) {
    if (u > v) swap(u, v);
    appear[{u,v}].push_back({query.size(), 1<<30});
  }
  void erase_edge(int u, int v) {
    if (u > v) swap(u, v);
    appear[{u,v}].back().snd = query.size();
  }
  int is_connected(int u, int v) {
    query.push_back({u,v});
    return query.size()-1;
  }

  vector<set<edge>> es;
  void insert(int l, int r, int k, int s, int t, edge e) {
    s = max(s, l); t = min(r, t);
    if (s >= t) return;
    if (s == l && t == r) {
      es[k].insert(e);
    } else {
      insert(l, (l+r)/2, 2*k+1, s, t, e);
      insert((l+r)/2, r, 2*k+2, s, t, e);
      if (es[2*k+1].count(e) && es[2*k+2].count(e)) {
        es[2*k+1].erase(e);
        es[2*k+2].erase(e);
        es[k].insert(e);
      }
    }
  }
  void rec(int l, int r, int k) {
    if (l >= r) return;
    for (edge e: es[k]) uf.unite(e.fst, e.snd);
    if (l+1 == r) {
      ans[l] = uf.find(query[l].fst, query[l].snd);
    } else {
      rec(l, (l+r)/2, 2*k+1);
      rec((l+r)/2, r, 2*k+2);
    }
    for (edge e: es[k]) uf.undo();
  }
  void solve() {
    int q = query.size();
    es.resize(4*q);
    for (auto a: appear) 
      for (auto b: a.snd) 
        insert(0, q, 0, b.fst, b.snd, a.fst);
    ans.assign(q, 0);
    rec(0, q, 0);
  }
};

int main() {
  offline_dynamic_connectivity solver(3);
  solver.is_connected(0,1);
  solver.add_edge(0,1);
  solver.is_connected(0,1);
  solver.is_connected(1,2);
  solver.add_edge(1,2);
  solver.is_connected(0,2);
  solver.erase_edge(0,1);
  solver.is_connected(0,2);
  solver.solve();
  for (int i = 0; i < solver.query.size(); ++i) {
    cout << solver.query[i].fst << " " << solver.query[i].snd << " " << solver.ans[i] << endl;
  }
}

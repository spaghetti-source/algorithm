// 
// Strongly Connected Component Comparison
//
// Gabow = Tarjan < Kosaraju.
// This difference can be ignorable in most cases.
//
// Conclusion:
//   Use Gabow (fast and memory efficient) or Kosaraju (short).
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

#include "tick.hh"

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

  vector<vector<int>> scc_kosaraju() {
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

  vector<vector<int>> scc_gabow() {
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

  vector<vector<int>> scc_tarjan() {
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
          scc.back().push_back(u);
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


template <class T>
ostream &operator<<(ostream &os, const vector<T> &v) {
  os << "[";
  for (int i = 0; i < v.size(); os << v[i++]) 
    if (i > 0) os << " ";
  os << "]";
  return os;
}
template <class T>
ostream &operator<<(ostream &os, const vector<vector<T>> &v) {
  os << "[";
  for (int i = 0; i < v.size(); os << v[i++]) 
    if (i > 0) os << endl << " ";
  os << "]";
  return os;
}

int main() {
  srand( time(0) );

  double ta = 0, tb = 0, tc = 0;
  for (int iter = 0; iter < 10; ++iter) {
    int n = 50000, m = n * 10;
    graph g(n);
    set<pair<int,int>> edges;
    while (edges.size() < m) {
      int u = rand() % n, v = rand() % n;
      if (u == v) continue;
      edges.insert({u, v});
    }
    for (auto p: edges) 
      g.add_edge(p.fst, p.snd);

    vector<int> order = {0,1,2};
    random_shuffle(all(order));
    int a, b, c;
    for (int i = 0; i < 3; ++i) {
      if (order[i] == 0) {
        tick();
        a = g.scc_kosaraju().size();
        ta += tick();
      } else if (order[i] == 1) {
        tick();
        b = g.scc_tarjan().size();
        tb += tick();
      } else {
        tick();
        c = g.scc_gabow().size();
        tc += tick();
      }
    }
    if (a != b || a != c || b != c) {
      cout << "***" << endl;
    }
  }
  cout << "kosaraju: " << ta << endl 
       << "tarjan:   " << tb << endl
       << "gabow:    " << tc << endl;

  /*
  int n, m;
  scanf("%d %d", &n, &m);
  graph g(n);
  for (int k = 0; k < m; ++k) {
    int i, j;
    scanf("%d %d", &i, &j);
    g.add_edge(i-1, j-1);
  }
  g.solve();
  */
}

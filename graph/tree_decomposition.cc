// 
// Tree Decomposition
//
// Description:
//   Let G = (V, E) be a graph. (T, F) is a tree-decomposition of G if
//   - X in T is a subsets of V, and \cup_{X in T} X = V.
//   - each e in E is contained in some X in T.
//   - for each v in V, { X in T : v in X } is connected.
//   A treewidth is max_{X in T} |X|-1. 
//
//   Simply say, if G is close to tree, it has a small treewidth.
//   - G is tree            <=> tw(G) = 1
//   - G is series-parallel <=> tw(G) = 2
//   - G is planar          ==> tw(G) = O(sqrt(n))
//   - G is complete        ==> tw(G) = n-1.
//
//   Many NP-hard problems can be solved in polynomial time
//   for constant treewidth graphs.
//
//
// Algorithm:
//   Minimum degree heuristics.
//   for i = 1 ... n:
//     Choose v with the smallest degree
//     Make v's neighbor to a clique, and eliminate v.
//   A star-representation is adopted to reduce the complexity. 
//
// Complexity:
//    it has no approximation guarantee, but practically works very well.
//    O(TW n), where TW is an obtained treewidth.

#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <vector>
#include <cstring>
#include <cassert>
#include <queue>
#include <functional>
#include <algorithm>

using namespace std;

#define fst first
#define snd second


struct graph {
  vector<pair<int, int>> edges;
  void add_edge(int u, int v) { 
    edges.push_back({u, v});
  }
  int n;
  vector<vector<int>> adj;
  void make_graph() {
    n = 0;
    for (auto e: edges) 
      n = max(n, max(e.fst, e.snd)+1);
    adj.resize(n);
    for (auto e: edges) {
      adj[e.fst].push_back(e.snd);
      adj[e.snd].push_back(e.fst);
    }
    for (auto &nbh: adj) {
      sort(nbh.begin(), nbh.end());
      nbh.erase(unique(nbh.begin(), nbh.end()), nbh.end());
    }
  }

  vector<int> parent; 
  int root(int v) { // union-find data structure
    if (parent[v] == v || parent[v] == -1) return v;
    return parent[v] = root(parent[v]);
  }
  void normalize(vector<int> &S) {
    for (auto &v: S) 
      v = root(v);
    sort(S.begin(), S.end());
    S.erase(unique(S.begin(), S.end()), S.end());
  }
  vector<int> neighbor(int u) {
    vector<int> nbh;
    normalize(adj[u]);
    for (auto v: adj[u]) {
      if (parent[v] == v) {
        nbh.push_back(v);
      } else {
        normalize(adj[v]);
        for (auto w: adj[v]) {
          if (parent[w] == w) {
            nbh.push_back(w);
          }
        }
      }
    }
    normalize(nbh);
    return nbh;
  }
  void contract(int u) {
    vector<int> live, dead;
    for (auto v: adj[u]) {
      if (parent[v] == v) live.push_back(v);
      else                dead.push_back(v);
    }
    parent[u] = -1;
    adj[u].swap(live);
    for (auto v: dead) {
      normalize(adj[v]);
      adj[u].insert(adj[u].end(), adj[v].begin(), adj[v].end());
      adj[v].clear();
      parent[v] = u;
    }
  }

  int solve() {
    typedef pair<int, int> node; // (deg, vertex)
    int tree_width = 0;
    parent.resize(n);
    priority_queue<node, vector<node>, greater<node>> Q;
    for (int u = 0; u < n; ++u) {
      parent[u] = u;
      Q.push(node(adj[u].size(), u));
    }
    while (!Q.empty()) {
      int deg = Q.top().fst; 
      int u = Q.top().snd;
      Q.pop();

      vector<int> nbh = neighbor(u);
      if (nbh.size() > deg) {
        Q.push({nbh.size(), u});
        continue;
      }
      tree_width = max(tree_width, nbh.size());
      vector<int> intersect;
      set_intersection(bags.back().begin(), bags.back(),end(),
                       nhb.begin(), nbh.end(), intersect);
      if (intersect.size() != nbh.size())
        bags.push_back(nbh);
      contract(u);
    }
    return tree_width;
  }
};


graph g;
void read_edges(char *file) {
  FILE *fp = fopen(file, "r");
  vector< pair<int, int> > edge;
  for (char buf[256]; fgets(buf, 256, fp); ) {
    int u, v;
    sscanf(buf, "%d %d", &u, &v);
    g.add_edge(u, v);
  }
  g.make_graph();
}

int main(int argc, char *argv[]) {
  read_edges(argv[1]);
  g.solve();
}

// 
// Least common ancestor by heavy-light decomposition
//
// Description:
//   For a rooted tree T, LCA(u,v) is a vertex u
//   that is the deepest node that is a common ancestor of u and v.
//
// Algorithm:
//   A heavy-light decomposition finds a partition of a tree into
//   a set of paths that is recursively defined as follows:
//   1) root is contained in a "heavy" path that spans root-to-child
//   2) other childs are linked by "light" edges and the subtrees are
//      recursively decomposed by the heavy-light decomposition.
//      Here, all subtrees have fewer nodes than the "heavy" path.
//
//   If two nodes u and v are located in the same heavy-path, LCA is
//   immediately obtained (shallower node is the LCA). Otherwise,
//   climbing light links until these are located in the same heavy-path.
//   By construction, only O(log n) light-link-climb is enough to reach
//   the root of the tree; thus it gives O(log n) algorithm.
//
// Complexity:
//   O(n) for preprocessing,
//   O(log n) for query
//
// Verified:
//   AOJ GRL_5C
//
#include <iostream>
#include <vector>
#include <cstdio>
#include <algorithm>
#include <functional>

using namespace std;

#define fst first
#define snd second
#define all(c) ((c).begin()), ((c).end())
#define TEST(s) if (!(s)) { cout << __LINE__ << " " << #s << endl; exit(-1); }


struct tree {
  int n;
  vector<vector<int>> adj;
  tree(int n) : n(n), adj(n) { }
  void add_edge(int s, int t) {
    adj[s].push_back(t);
    adj[t].push_back(s);
  }
  vector<int> size, depth, head, parent;
  void rootify(int r) {
    size = depth = parent = head = vector<int>(n, -1);
    function<int (int,int)> dfs = [&](int u, int p) {
      parent[u] = p;
      depth[u] = depth[p]+1;
      for (int v: adj[u]) 
        if (v != p) size[u] += dfs(v, u);
      return ++size[u];
    }; dfs(r, r);
    function<void (int,int)> dec = [&](int u, int s) {
      head[u] = s;
      int z = -1;
      for (int v: adj[u]) 
        if (head[v] < 0 && (z < 0 || size[z] < size[v])) z = v;
      for (int v: adj[u]) 
        if (head[v] < 0) dec(v, v == z ? s : v);
    }; dec(r, r);
  }
  int lca(int u, int v) {
    while (head[u] != head[v]) 
      if (depth[head[u]] < depth[head[v]]) v = parent[head[v]];
      else                                 u = parent[head[u]];
    return depth[u] < depth[v] ? u : v;
  }
};

int main() {
  int n;
  scanf("%d", &n);
  tree T(n);
  for (int u = 0; u < n ; ++u) {
    int k; 
    scanf("%d", &k);
    for (int j = 0; j < k; ++j) {
      int v;
      scanf("%d", &v);
      T.add_edge(u, v);
    }
  }
  T.rootify(0);
  int q;
  scanf("%d", &q);
  for (int i = 0; i < q; ++i) {
    int u, v;
    scanf("%d %d", &u, &v);
    printf("%d\n", T.lca(u, v));
  }
}

// 
// Least common ancestor by euler tour + sparse table
//
// Description:
//   For a rooted tree T, LCA(u,v) is a vertex u
//   that is the deepest node that is a common ancestor of u and v.
//
// Algorithm:
//   This first finds an euler tour of the tree.
//   Then, RMQ(pos[u], pos[v]) = LCA(u, v), where
//   RMQ is the range minimum query between i and j,
//   where the weight is defined by the depth of the node.
//
// Complexity:
//   O(n log n) for preprocessing,
//   O(1) for query
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
  vector<int> pos, tour, depth;
  vector<vector<int>> table;
  int argmin(int i, int j) { return depth[i] < depth[j] ? i : j; }
  void rootify(int r) {
    pos.resize(n);
    function<void (int,int,int)> dfs = [&](int u, int p, int d) {
      pos[u] = depth.size();
      tour.push_back(u);
      depth.push_back(d);
      for (int v: adj[u]) {
        if (v != p) {
          dfs(v, u, d+1);
          tour.push_back(u);
          depth.push_back(d);
        }
      }
    }; dfs(r, r, 0);
    int logn = sizeof(int)*__CHAR_BIT__-1-__builtin_clz(tour.size()); // log2
    table.resize(logn+1, vector<int>(tour.size()));
    iota(all(table[0]), 0);
    for (int h = 0; h < logn; ++h) 
      for (int i = 0; i+(1<<h) < tour.size(); ++i)
        table[h+1][i] = argmin(table[h][i], table[h][i+(1<<h)]);
  }
  int lca(int u, int v) {
    int i = pos[u], j = pos[v]; if (i > j) swap(i, j);
    int h = sizeof(int)*__CHAR_BIT__-1-__builtin_clz(j-i); // = log2
    return i == j ? u : tour[argmin(table[h][i], table[h][j-(1<<h)])];
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

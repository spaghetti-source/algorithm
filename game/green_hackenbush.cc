//
// Green Hackenbush
//
// Description:
//   Consider a two player game on a graph with a specified vertex (root).
//   In each turn, a player eliminates one edge. 
//   Then, if a subgraph that is disconnected from the root, it is removed. 
//   If a player cannot select an edge (i.e., the graph is singleton), 
//   he will lose.
//
//   Compute the Grundy number of the given graph.
//
// Algorithm:
//   We use two principles:
//     1. Colon Principle: Grundy number of a tree is the xor of 
//        Grundy number of child subtrees.
//        (Proof: easy).
//
//     2. Fusion Principle: Consider a pair of adjacent vertices u, v
//        that has another path (i.e., they are in a cycle). Then,
//        we can contract u and v without changing Grundy number.
//        (Proof: difficult)
//
//   We first decompose graph into two-edge connected components.
//   Then, by contracting each components by using Fusion Principle,
//   we obtain a tree (and many self loops) that has the same Grundy
//   number to the original graph. By using Colon Principle, we can
//   compute the Grundy number.
//
// Complexity:
//   O(m + n).
//
// Verified:
//   SPOJ 1477: Play with a Tree
//   IPSC 2003 G: Got Root?
//
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

struct hackenbush {
  int n;
  vector<vector<int>> adj;

  hackenbush(int n) : n(n), adj(n) { }
  void add_edge(int u, int v) {
    adj[u].push_back(v);
    if (u != v) adj[v].push_back(u);
  }

  // r is the only root connecting to the ground
  int grundy(int r) {
    vector<int> num(n), low(n);
    int t = 0;
    function<int(int,int)> dfs = [&](int p, int u) { 
      num[u] = low[u] = ++t;
      int ans = 0;
      for (int v: adj[u]) {
        if (v == p) { p += 2*n; continue; } 
        if (num[v] == 0) {
          int res = dfs(u, v);
          low[u] = min(low[u], low[v]);
          if (low[v] > num[u]) ans ^= (1 + res) ^ 1; // bridge
          else                 ans ^= res;           // non bridge
        } else low[u] = min(low[u], num[v]);
      }
      if (p > n) p -= 2*n;
      for (int v: adj[u]) 
        if (v != p && num[u] <= num[v]) ans ^= 1;
      return ans;
    };
    return dfs(-1, r);
  }
};

int main() {
  int cases; scanf("%d", &cases);
  for (int icase = 0; icase < cases; ++icase) {
    int n; scanf("%d", &n);
    vector<int> ground(n);
    int r;
    for (int i = 0; i < n; ++i) {
      scanf("%d", &ground[i]);
      if (ground[i] == 1) r = i;
    }
    int ans = 0;
    hackenbush g(n);
    for (int i = 0; i < n-1; ++i) {
      int u, v;
      scanf("%d %d", &u, &v);
      --u; --v;
      if (ground[u]) u = r;
      if (ground[v]) v = r;
      if (u == v) ans ^= 1;
      else g.add_edge(u, v);
    }
    int res = ans ^ g.grundy(r);
    printf("%d\n", res != 0);
  }
}

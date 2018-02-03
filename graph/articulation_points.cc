// 
// Block-Cut Tree (Articulation points / Biconnected components)
//
// Description:
//   Let G = (V, E). If G-v is disconnected, v in V is said to
//   be an articulation point. If G has no articulation points,
//   it is said to be biconnected.
//
//   A biconnected component is a maximal biconnected subgraph.
//   The algorithm finds all articulation points and biconnected
//   components. It can be obtained by the Hopcroft-Tarjan DFS.
//
//   We maintain the biconnected component decomposition by the
//   block tree whose vertices are the blocks and the articulation
//   points. By contracting the graph by the articulation points,
//   we obtain the intersection graph of the blocks.
//
// Complexity:
//   O(n + m).
//
// Verified:
//   AOJ_GRL_3_A (articulation point)
//
// References:
//   J. Hopcroft and R. E. Tarjan (1973):
//   Efficient algorithms for graph manipulation.
//   Communications of the ACM, vol.16, no.6, pp.372-378.
//

// g++ -std=c++17 -O3 -fmax-errors=1 -fsanitize=undefined
#include <bits/stdc++.h>

using namespace std;

#define fst first
#define snd second
#define all(c) ((c).begin()), ((c).end())
#define TEST(s) if (!(s)) { cout << __LINE__ << " " << #s << endl; exit(-1); }

struct Graph {
  int n;
  vector<vector<int>> adj;
  Graph(int n) : n(n), adj(n) { }
  void addEdge(int u, int v) {
    adj[u].push_back(v);
    adj[v].push_back(u);
  }
};

//
// bcc.n-bcc.block.size() is the number of articulation points
// index[u] is the corresponding node in the block-cut tree.
// Here, if u in block[k] then
//   if u is an articulation point, k in adj[index[u]] 
//   otherwise, index[u] = k
//
struct BiconnectedComponents : Graph {
  vector<int> is_articulation, index;
  vector<vector<int>> block;

  BiconnectedComponents(Graph g) : Graph(0) {
    vector<int> low(g.n), num(g.n), cur(g.n), par(g.n, -1), path;
    is_articulation.resize(g.n);
    for (int s = 0; s < g.n; ++s) {
      if (num[s]) continue;
      int time = 0;
      vector<int> stack = {s};
      while (!stack.empty()) {
        int u = stack.back();
        if (cur[u] == 0) {
          low[u] = num[u] = ++time;
          path.push_back(u);
        }
        if (cur[u] == g.adj[u].size()) {
          stack.pop_back();
        } else if (cur[u] >= 0) {
          int v = g.adj[u][cur[u]++];
          if (num[v] == 0) {
            cur[u] = ~cur[u];
            stack.push_back(v);
          } else if (v != par[u]) {
            low[u] = min(low[u], num[v]);
          }
        } else {
          cur[u] = ~cur[u];
          int v = g.adj[u][cur[u]-1];
          low[u] = min(low[u], low[v]);
          if (num[u] <= low[v]) {
            is_articulation[u] = (num[u] > 1 || num[v] > 2);
            block.push_back({u});
            while (block.back().back() != v) {
              block.back().push_back(path.back());
              path.pop_back();
            }
          }
        }
      }
    }
    index.resize(g.n); // make a block tree
    n = block.size();
    for (int u = 0; u < g.n; ++u) 
      if (is_articulation[u]) index[u] = n++;
    adj.resize(n);
    for (int k = 0; k < block.size(); ++k) {
      for (int u: block[k]) {
        if (!is_articulation[u]) index[u] = k;
        else addEdge(k, index[u]);
      }
    }
  }
};

void AOJ_GRL_3_A() {
  int n, m; 
  scanf("%d %d", &n, &m);
  Graph g(n);
  for (int i = 0; i < m; ++i) {
    int u, v; scanf("%d %d", &u, &v);
    g.addEdge(u, v);
  }
  BiconnectedComponents bcc(g);
  for (int u = 0; u < g.n; ++u) 
    if (bcc.is_articulation[u]) cout << u << endl;
}

int main() {
  AOJ_GRL_3_A();
  //SPOJ_SUBMERGE();
  //test();
  /*
  */
}


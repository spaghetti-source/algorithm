// 
// Bridge-Block Tree (Bridge / Two-edge connected component)
//
// Description:
//   Let G = (V, E). e in E is said to be a cut edge if G-e is 
//   disconnected If G has no cut edges, it is said to be two-edge 
//   connected.
//   A two-edge connected component is a maximal two-edge connected
//   subgraph. The algorithm finds all bridges with the two-edge
//   connected components.
//
//   We maintain the two-edge connected component decomposition by
//   a bridge-block tree whose nodes are the two-edge connected 
//   components.
//
//
// Complexity:
//   O(n + m).
//
// Verified:
//
// References:
//

#include <iostream>
#include <vector>
#include <cstdio>
#include <algorithm>
#include <functional>
#include <unordered_set>

using namespace std;

#define fst first
#define snd second
#define all(c) ((c).begin()), ((c).end())

struct Graph {
  int n;
  vector<vector<int>> adj;
  Graph(int n) : n(n), adj(n) { }
  void addEdge(int src, int dst) {
    adj[src].push_back(dst);
    adj[dst].push_back(src);
  }
};

struct BridgeBlockTree : Graph {
  vector<int> index;         // index[u] is the block index containing u
  vector<vector<int>> block; // u in block[k] <=> index[u] == k

  BridgeBlockTree(Graph g) : Graph(0) {
    index.assign(g.n, -1);
    vector<int> num(g.n), par(g.n,-1), cur(g.n);

    for (int s = 0; s < g.n; ++s) {
      if (num[s]) continue;
      int time = 0;
      vector<int> snum, path, stack = {s};
      while (!stack.empty()) {
        int u = stack.back();
        if (cur[u] == 0) {
          num[u] = ++time;
          path.push_back(u); 
          snum.push_back(num[u]);
        }
        if (cur[u] == g.adj[u].size()) {
          if (num[u] == snum.back()) {
            snum.pop_back();
            block.push_back({});
            while (1) {
              int w = path.back(); path.pop_back(); 
              block.back().push_back(w);
              index[w] = block.size()-1;
              if (u == w) break;
            }
          }
          stack.pop_back();
        } else {
          int v = g.adj[u][cur[u]++];
          if (!num[v]) {
            par[v] = u;
            stack.push_back(v);
          } else if (v != par[u] && index[v] < 0) {
            while (snum.back() > num[v]) snum.pop_back();
          }
        }
      }
    }
    n = block.size();
    adj.resize(n);
    for (int u = 0; u < g.n; ++u)
      if (par[u] >= 0 && index[u] != index[par[u]])
        addEdge(index[u], index[par[u]]);
  }
};


// === tick a time ===
#include <ctime>
double tick() {
  static clock_t oldtick;
  clock_t newtick = clock();
  double diff = 1.0*(newtick - oldtick) / CLOCKS_PER_SEC;
  oldtick = newtick;
  return diff;
}

int main() {
  Graph g(6);
  g.addEdge(0, 1);
  g.addEdge(1, 2);
  g.addEdge(2, 0);
  g.addEdge(2, 3);
  g.addEdge(3, 4);
  g.addEdge(4, 5);
  g.addEdge(5, 3);
  BridgeBlockTree t(g);

  cout << t.n << endl;
  for (int u = 0; u < t.n; ++u) {
    for (int v: t.adj[u]) {
      cout << u << " " << v << endl;
    }
  }

  for (auto B: t.block) {
    for (int u: B) {
      cout << u << " ";
    }
    cout << endl;
  }
  for (int u = 0; u < 6; ++u) {
    cout << t.index[u] << " ";
  }
  cout << endl;

  //g.bridgeless_component();
  /*
  for (int n, m; ~scanf("%d %d", &n, &m) && n; ) {
    graph g(n);
    for (int i = 0; i < m; ++i) {
      int u, v; scanf("%d %d", &u, &v);
      g.add_edge(u-1, v-1);
    }
    g.biconnected_components();
  }
  */
}

// 
// k-Core Decomposition
//
// Description:
//   This finds a k-Core decomposition of given graph G. 
//   Here, k-core decomposition is a layered decomposition 
//     V \supset C_1 \supset C_2 \supset ... \supset C_m
//   such that each C_k is a k-connected subgraph.
//   The largest k is is known as the degeneracy of graph.
//
// Algorithm:
//   Greedy. Pick smallest connected vertex u and 
//   remove u and the adjacent edges.
//
// Complexity:
//   O(m log n).
//
// References:
//   S. B. Seidman (1983):
//   Network structure and minimum degree.
//   Social Networks, vol. 5, 269-287.
//
#include <iostream>
#include <vector>
#include <cstdio>
#include <cstdlib>
#include <queue>
#include <algorithm>

using namespace std;

struct edge {
  int src, dst;
};
struct k_core_decomposition {
  vector<edge> edges;
  void add_edge(int src, int dst) {
    edges.push_back({src, dst});
  }
  int n;
  vector<vector<edge>> adj;
  void make_graph(int n_ = 0) {
    n = n_;
    for (auto e: edges) 
      n = max(n, max(e.src, e.dst)+1);
    adj.resize(n);
    for (auto e: edges) 
      adj[e.src].push_back(e);
  }
  // A subgraph C_k := { v : k[v] >= k } is a maximal k-connected subgraph
  vector<int> kindex;
  int solve() {
    typedef pair<int, int> node;
    priority_queue<node, vector<node>, greater<node>> Q;
    kindex.assign(n, -1);
    vector<int> degree(n);
    for (int u = 0; u < n; ++u) 
      Q.push({degree[u] = adj[u].size(), u});
    while (!Q.empty()) {
      auto p = Q.top(); Q.pop();
      if (degree[p.second] < p.first) continue;
      kindex[p.second] = degree[p.second];
      for (edge e: adj[p.second]) 
        if (kindex[e.dst] < 0)
          Q.push({--degree[e.dst], e.dst});
    }
    return *max_element(kindex.begin(), kindex.end());
  }
};

int main() {
  cout << plus<int>()(2,3) << endl;
}

//
// Topological Sort
//
//
// Description:
//
//   Let G = (V, E) be a graph. An ordering ord: [n] -> V is a topological 
//   ordering if i > j then there is no edge from ord[i] to ord[j].
//   G has a topological ordering if and only if G is DAG.
//
//   A topological order can be obtained in O(n + m) time by using
//   an iterative method (Kahn's algorithm) or a recursive method 
//   (by Tarjan's algorithm). The following implementation is a 
//   Kahn's algorithm.
//
//   Note that if you want to find the all topological orders,
//
//
// Complexity:
//
//   O(n + m)
//
//
// Verified:
//
//   AOJ GPL_4_B 
//
// References:
//
//   Arthur B. Kahn (1962):
//   "Topological sorting of large networks".
//   Communications of the ACM, 5 (11): 558--562.
//
#include <bits/stdc++.h>
using namespace std;

#define fst first
#define snd second
#define all(c) ((c).begin()), ((c).end())

struct Graph {
  int n;
  vector<vector<int>> adj;
  Graph(int n) : n(n), adj(n) { }
  void addEdge(int u, int v) {
    adj[u].push_back(v);
  }
};

// return empty list if g has no topological order
vector<int> topologicalSort(Graph g) {
  vector<int> deg(g.n);
  for (int u = 0; u < g.n; ++u)
    for (int v: g.adj[u]) ++deg[v];
  vector<int> stack; 
  for (int u = 0; u < g.n; ++u)
    if (!deg[u]) stack.push_back(u);

  vector<int> order;
  while (!stack.empty()) {
    int u = stack.back(); stack.pop_back();
    order.push_back(u);
    for (int v: g.adj[u]) 
      if (!--deg[v]) stack.push_back(v);
  }
  return order.size() == g.n ? order : vector<int>();
}

int main() {
  int n, m; cin >> n >> m;
  Graph g(n);
  for (int i = 0; i < m; ++i) {
    int u, v; cin >> u >> v;
    g.addEdge(u, v);
  }
  auto ord = topologicalSort(g);
  for (int i = 0; i < ord.size(); ++i) {
    if (i > 0) cout << " ";
    cout << ord[i];
  }
}

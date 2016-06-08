//
// Maximum Flow (Goldberg-Tarjan, aka. Push-Relabel, Preflow-Push)
// 
// Description:
//   Given a directed network G = (V, E) with edge capacity c: E->R.
//   The algorithm finds a maximum flow. 
//
// Algorithm:
//   Goldberg-Tarjan's push-relabel algorithm with gap-heuristics.
//
// Complexity:
//   O(n^3)
// 
// Verified:
//   SPOJ FASTFLOW
//
// Reference:
//   B. H. Korte and Jens Vygen (2008):
//   Combinatorial Optimization: Theory and Algorithms.
//   Springer Berlin Heidelberg. 
//

#include <iostream>
#include <vector>
#include <cstdio>
#include <queue>
#include <algorithm>
#include <functional>

using namespace std;

#define fst first
#define snd second
#define all(c) ((c).begin()), ((c).end())

const long long INF = (1ll << 50);
struct graph {
  typedef long long flow_type;
  struct edge {
    int src, dst;
    flow_type capacity, flow;
    size_t rev;
  };
  int n;
  vector<vector<edge>> adj;
  graph(int n) : n(n), adj(n) { }

  void add_edge(int src, int dst, int capacity) {
    adj[src].push_back({src, dst, capacity, 0, adj[dst].size()});
    adj[dst].push_back({dst, src, 0, 0, adj[src].size() - 1});
  }

  flow_type max_flow(int s, int t) {
    vector<flow_type> excess(n);
    vector<int> dist(n), active(n), count(2*n);
    queue<int> Q;
    auto enqueue = [&](int v) { 
      if (!active[v] && excess[v] > 0) { active[v] = true; Q.push(v); } 
    };
    auto push = [&](edge &e) {
      flow_type f = min(excess[e.src], e.capacity - e.flow);
      if (dist[e.src] <= dist[e.dst] || f == 0) return;
      e.flow += f;
      adj[e.dst][e.rev].flow -= f;
      excess[e.dst] += f;    
      excess[e.src] -= f;
      enqueue(e.dst);
    };

    dist[s] = n; active[s] = active[t] = true;
    count[0] = n-1; count[n] = 1;
    for (int u = 0; u < n; ++u)
      for (auto &e: adj[u]) e.flow = 0;
    for (auto &e: adj[s]) {
      excess[s] += e.capacity;
      push(e);
    }
    while (!Q.empty()) {
      int u = Q.front(); Q.pop();
      active[u] = false;

      for (auto &e: adj[u]) push(e);
      if (excess[u] > 0) {
        if (count[dist[u]] == 1) {
          int k = dist[u]; // Gap Heuristics
          for (int v = 0; v < n; v++) {
            if (dist[v] < k) continue;
            count[dist[v]]--;
            dist[v] = max(dist[v], n+1);
            count[dist[v]]++;
            enqueue(v);
          }
        } else {
          count[dist[u]]--; // Relabel
          dist[u] = 2*n;
          for (auto &e: adj[u]) 
            if (e.capacity > e.flow)
              dist[u] = min(dist[u], dist[e.dst] + 1);
          count[dist[u]]++;
          enqueue(u);
        }
      }
    }

    flow_type flow = 0;
    for (auto e: adj[s]) flow += e.flow;
    return flow;
  }
};

int main() {
  for (int n, m; scanf("%d %d", &n, &m) == 2; ) {
    graph g(n);
    for (int i = 0; i < m; ++i) {
      int u, v, w;
      scanf("%d %d %d", &u, &v, &w);
      g.add_edge(u, v, w);
    }
    printf("%d\n", g.max_flow(0, n-1));
  }
}

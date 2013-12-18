//
// Minimum Cost Flow (Tomizawa, Edmonds-Karp)
//
// 
// Description:
//   Given a directed graph G = (V,E) with nonnegative capacity c and cost w.
//   The algorithm find a maximum s-t flow of G with minimum cost.
// 
//
// Algorithm:
//   Tomizawa (1971), and Edmonds and Karp (1972)'s 
//   successive shortest path algorithm,
//   which is also known as the primal-dual method.
//
//
// Complexity:
//   O(F m log n), where F is the amount of maximum flow.
//
// 
// References:
//   N. Tomizawa (1971):
//   On some techniques useful for solution of transportation network problems.
//   Networks, vol. pp. 173-194.
//
//   J. Edmonds and R.M. Karp (1972):
//   Theoretical improvements in algorithmic efficiency for network flow problems.
//   Journal of ACM, vol. 19, pp. 248-264.
//
//
// Historical Note:
//   The successive shortest path type algorithm was developped
//   independently by Jewell (1958), Iri (1960), and Busacker and Gowen (1961).
//   Later, Tomizawa (1971), and Edmonds and Karp (1972) independently 
//   suggested to use vertex potential in shortest path algorithm.
// 


#include <iostream>
#include <vector>
#include <queue>
#include <cstdio>
#include <cstdlib>
#include <map>
#include <cmath>
#include <cstring>
#include <functional>
#include <algorithm>
#include <unordered_map>

using namespace std;


typedef long long cost_type;
typedef long long flow_type;
const flow_type INF = 99999999;
struct min_cost_max_flow {
  struct edge {
    int src, dst;
    flow_type residue; 
    cost_type cost;
    int rev;
  };
  vector<edge> edges;
  void add_edge(int src, int dst, flow_type cap, cost_type cost) {
    edges.push_back({src, dst, cap, cost});
  }
  int n;
  vector<vector<edge>> adj;
  void make_graph(int n_ = 0) {
    n = n_;
    for (auto e: edges) 
      n = max(n, max(e.src, e.dst)+1);
    adj.resize(n);
    for (auto e: edges) {
      edge r = {e.dst, e.src, 0, -e.cost};
      adj[e.src].push_back(e);
      adj[r.src].push_back(r);
      adj[e.src].back().rev = adj[r.src].size()-1;
      adj[r.src].back().rev = adj[e.src].size()-1;
    }
  }
  vector<cost_type> potential;
  cost_type rcost(const edge &e) {
    return e.cost + potential[e.src] - potential[e.dst];
  }
  void bellman_ford(int s) {
    for (int k = 0; k < n; ++k) 
      for (int u = 0; u < n; ++u)
        for (auto e: adj[u]) 
          if (e.residue > 0 && rcost(e) < 0)
            potential[e.dst] += rcost(e);
  }
  vector<cost_type> dist;
  vector<edge*> back;
  flow_type dijkstra(int s, int t) {
    fill(dist.begin(), dist.end(), INF);
    dist[s] = 0;
    typedef pair<cost_type, int> node;
    priority_queue<node, vector<node>, greater<node>> Q;
    Q.push({0, s});
    while (!Q.empty()) {
      auto p = Q.top(); Q.pop();
      if (dist[p.second] < p.first) continue;
      if (p.second == t) break;
      for (edge &e: adj[p.second]) {
        if (e.residue <= 0) continue;
        if (dist[e.dst] > dist[e.src] + rcost(e)) {
          dist[e.dst] = dist[e.src] + rcost(e);
          back[e.dst] = &e;
          Q.push({dist[e.dst], e.dst});
        }
      }
    }
    return dist[t];
  }
  pair<flow_type, cost_type> solve(int s, int t) {
    flow_type flow = 0;
    cost_type cost = 0;

    potential.assign(n, 0);
    dist.resize(n);
    back.resize(n);
    //bellman_ford(s); // remove negative costs

    while (dijkstra(s, t) < INF) {
      for (int u = 0; u < n; ++u)
        if (dist[u] < dist[t]) 
          potential[u] += dist[u] - dist[t];

      flow_type f = INF;
      for (int u = t; u != s; u = back[u]->src) 
        f = min(f, back[u]->residue);
      for (int u = t; u != s; u = back[u]->src) {
        back[u]->residue -= f;
        adj[back[u]->dst][back[u]->rev].residue += f;
      }
      flow += f;
      cost += f * (potential[t] - potential[s]);
    }
    return {flow, cost};
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
  for (int seed = 0; seed < 100; ++seed) {
    cout << "--------------------" << endl;

  cout << "seed = " << seed << endl;
  srand(seed);

  int n = 500;
  min_cost_max_flow mcmf;

  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < i; ++j) {
      if (rand() % 2) mcmf.add_edge(i, j, 1+rand()%100, 1+rand()%100);
      if (rand() % 2) mcmf.add_edge(j, i, 1+rand()%100, 1+rand()%100);
    }
  }
  mcmf.make_graph(n);
  int s = 0, t = n-1;
  tick();
  auto a = mcmf.solve(s, t);
  cout << a.first << " " << a.second << " " << tick() << endl;
  }
}

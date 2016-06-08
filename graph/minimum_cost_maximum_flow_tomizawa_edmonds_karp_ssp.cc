//
// Minimum Cost Maximum Flow (Tomizawa, Edmonds-Karp's successive shortest path)
//
// Description:
//   Given a directed graph G = (V,E) with nonnegative capacity c and cost w.
//   The algorithm find a maximum s-t flow of G with minimum cost.
// 
// Algorithm:
//   Tomizawa (1971), and Edmonds and Karp (1972)'s 
//   successive shortest path algorithm,
//   which is also known as the primal-dual method.
//
// Complexity:
//   O(F m log n), where F is the amount of maximum flow.
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

#define fst first
#define snd second
#define all(c) ((c).begin()), ((c).end())
#define TEST(s) if (!(s)) { cout << __LINE__ << " " << #s << endl; exit(-1); }

const long long INF = 99999999;
struct graph {
  typedef int flow_type;
  typedef int cost_type;
  struct edge {
    int src, dst;
    flow_type capacity, flow;
    cost_type cost;
    size_t rev;
  };
  vector<edge> edges;
  void add_edge(int src, int dst, flow_type cap, cost_type cost) {
    adj[src].push_back({src, dst, cap, 0, cost, adj[dst].size()});
    adj[dst].push_back({dst, src, 0, 0, -cost, adj[src].size()-1});
  }
  int n;
  vector<vector<edge>> adj;
  graph(int n) : n(n), adj(n) { }

  pair<flow_type, cost_type> min_cost_max_flow(int s, int t) {
    flow_type flow = 0;
    cost_type cost = 0;

    for (int u = 0; u < n; ++u) // initialize
      for (auto &e: adj[u]) e.flow = 0;

    vector<cost_type> p(n, 0);

    auto rcost = [&](edge e) { return e.cost + p[e.src] - p[e.dst]; };
    for (int iter = 0; ; ++iter) {
      vector<int> prev(n, -1); prev[s] = 0;
      vector<cost_type> dist(n, INF); dist[s] = 0;
      if (iter == 0) { // use Bellman-Ford to remove negative cost edges
        vector<int> count(n); count[s] = 1;
        queue<int> que; 
        for (que.push(s); !que.empty(); ) {
          int u = que.front(); que.pop();
          count[u] = -count[u];
          for (auto &e: adj[u]) {
            if (e.capacity > e.flow && dist[e.dst] > dist[e.src] + rcost(e)) {
              dist[e.dst] = dist[e.src] + rcost(e);
              prev[e.dst] = e.rev;
              if (count[e.dst] <= 0) {
                count[e.dst] = -count[e.dst] + 1;
                que.push(e.dst);
              }
            }
          }
        }
      } else { // use Dijkstra 
        typedef pair<cost_type, int> node;
        priority_queue<node, vector<node>, greater<node>> que;
        que.push({0, s});
        while (!que.empty()) {
          node a = que.top(); que.pop();
          if (a.snd == t) break;
          if (dist[a.snd] > a.fst) continue;
          for (auto e: adj[a.snd]) {
            if (e.capacity > e.flow && dist[e.dst] > a.fst + rcost(e)) {
              dist[e.dst] = dist[e.src] + rcost(e);
              prev[e.dst] = e.rev;
              que.push({dist[e.dst], e.dst});
            }
          }
        }
      }
      if (prev[t] == -1) break;

      for (int u = 0; u < n; ++u) 
        if (dist[u] < dist[t]) p[u] += dist[u] - dist[t];

      function<flow_type(int,flow_type)> augment = [&](int u, flow_type cur) {
        if (u == s) return cur;
        edge &r = adj[u][prev[u]], &e = adj[r.dst][r.rev];
        flow_type f = augment(e.src, min(e.capacity - e.flow, cur));
        e.flow += f; r.flow -= f;
        return f;
      };
      flow_type f = augment(t, INF);
      flow += f;
      cost += f * (p[t] - p[s]);
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
  graph g(n);

  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < i; ++j) {
      if (rand() % 2) g.add_edge(i, j, 1+rand()%100, 1+rand()%100);
      if (rand() % 2) g.add_edge(j, i, 1+rand()%100, 1+rand()%100);
    }
  }
  int s = 0, t = n-1;
  tick();
  auto a = g.min_cost_max_flow(s, t);
  cout << a.fst << " " << a.snd << " " << tick() << endl;
}
}

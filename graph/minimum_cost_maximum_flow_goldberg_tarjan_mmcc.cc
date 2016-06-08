//
// Minimum Cost Maximum Flow (Goldberg-Tarjan's minimum mean cycle canceling)
//
// Description:
//   Given a directed graph G = (V,E) with nonnegative capacity c and cost w.
//   The algorithm find a maximum s-t flow of G with minimum cost.
// 
// Algorithm:
//   Goldberg-Tarjan (Tomizawa (1971)'s algorithm.
//   It first finds a feasible (i.e., maximum) flow. Then it successively
//   finds negative cycles, and cancels them.
//   By finding minimum mean cycle, it finishes in strongly polynomial time.
//
// Complexity:
//   O(n^2 m^2 log(nC), where C is the maximum capacity. 
//   Practically, this is much slower than other algorithms.
// 
// References:
//   A. Goldberg and R. E. Tarjan (1989):
//   Finding minimum-cost circulations by canceling negative cycles.
//   Journal of the ACM, vol. 36, no. 4, pp. 873--886.
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

#define DEBUG 0

const long long INF = 99999999;
namespace NegativeCycleCanceling {
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

  flow_type max_flow(int s, int t) {
    vector<int> level(n), iter(n);
    function<int(void)> levelize = [&]() { // foward levelize
      level.assign(n, -1); level[s] = 0;
      queue<int> Q; Q.push(s);
      while (!Q.empty()) {
        int u = Q.front(); Q.pop();
        if (u == t) break;
        for (auto &e: adj[u]) {
          if (e.capacity > e.flow && level[e.dst] < 0) {
            Q.push(e.dst);
            level[e.dst] = level[u] + 1;
          }
        }
      }
      return level[t];
    };
    function<flow_type(int, flow_type)> augment = [&](int u, flow_type cur) {
      if (u == t) return cur;
      for (int &i = iter[u]; i < adj[u].size(); ++i) {
        edge &e = adj[u][i], &r = adj[e.dst][e.rev];
        if (e.capacity > e.flow && level[u] < level[e.dst]) {
          flow_type f = augment(e.dst, min(cur, e.capacity - e.flow));
          if (f > 0) {
            e.flow += f;
            r.flow -= f;
            return f;
          }
        }
      }
      return flow_type(0);
    };
    for (int u = 0; u < n; ++u) // initialize
      for (auto &e: adj[u]) e.flow = 0;

    flow_type flow = 0;
    while (levelize() >= 0) {
      fill(all(iter), 0);
      for (flow_type f; (f = augment(s, INF)) > 0; )
        flow += f;
    }
    return flow;
  }

  void disp() {
    if (!DEBUG) return;
    cout << "--------" << endl;
    for (int u = 0; u < n; ++u)
      for (auto e: adj[u]) {
        cout << e.src << " " << e.dst << " " << e.capacity << " " << e.flow << " " << e.cost << endl;
      }
    cost_type cost = 0;
    for (int u = 0; u < n; ++u) 
      for (auto &e: adj[u])
        if (e.flow > 0) cost += e.flow * e.cost;
    cout << "cost = " << cost << endl;
  }

  bool minimum_mean_cycle_cancel() {
    vector<vector<flow_type>> dist(n+1, vector<flow_type>(n));
    vector<vector<int>> prev(n+1, vector<int>(n, -1));
    fill(all(prev[0]), 0);

    for (int k = 0; k < n; ++k) {
      for (int u = 0; u < n; ++u) {
        if (prev[k][u] < 0) continue;
        for (auto e: adj[u]) {
          if (e.capacity <= e.flow) continue;
          if (prev[k+1][e.dst] < 0 || dist[k+1][e.dst] > dist[k][u] + e.cost) {
            dist[k+1][e.dst] = dist[k][u] + e.cost;
            prev[k+1][e.dst] = e.rev;
          }
        }
      }
    }
    flow_type num = INF;
    int v, den = 1;
    for (int u = 0; u < n; ++u) {
      flow_type num_u = -INF;
      int den_u = 0;
      for (int k = 0; k < n; ++k) {
        if (prev[k][u] < 0) continue;
        if (num_u * (n-k) < (dist[n][u]-dist[k][u]) * den_u) {
            num_u = (dist[n][u] - dist[k][u]); den_u = n-k;
        }
      }
      if (den_u > 0 && num * den_u > num_u * den) {
        num = num_u; den = den_u; v = u;
      }
    }
    if (num >= 0) return false;
    vector<int> back(n, -1);
    for (int k = n; back[v] < 0; --k) {
      back[v] = prev[k][v];
      edge &r = adj[v][back[v]];
      v = r.dst;
    }
    function<flow_type(int,flow_type)> augment = [&](int u, flow_type cur) {
      if (cur < INF && u == v) return cur;
      edge &r = adj[u][back[u]], &e = adj[r.dst][r.rev];
      flow_type f = augment(r.dst, min(e.capacity - e.flow, cur));
      e.flow += f;
      r.flow -= f;
      return f;
    }; 
    augment(v, INF);
    return true;
  }
  pair<flow_type, cost_type> min_cost_max_flow(int s, int t) {
    flow_type flow = max_flow(s, t);
    while (minimum_mean_cycle_cancel());
    cost_type cost = 0;
    for (int u = 0; u < n; ++u) 
      for (auto &e: adj[u])
        if (e.flow > 0) cost += e.flow * e.cost;
    return {flow, cost};
  }
};
}

// for comparison
namespace Tomizawa_Edmonds_Karp {
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
}

  

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
  for (int seed = 0; seed < 10000; ++seed) {
    cout << "--------------------" << endl;

    cout << "seed = " << seed << endl;
    srand(seed);

    int n = 200;
    Tomizawa_Edmonds_Karp::graph g1(n);
    NegativeCycleCanceling::graph g2(n);

    /*
    g2.add_edge(0, 1, 1, 1);
    g2.add_edge(1, 2, 1, 1);
    g2.add_edge(2, 3, 1, 1);
    g2.add_edge(1, 3, 1, 4);
    */
    /*
    g2.add_edge(0, 1, 1, 1);
    g2.add_edge(1, 2, 1, 1);
    g2.add_edge(2, 3, 1, 1);
    g2.add_edge(1, 3, 1, 1);
    */

    for (int i = 0; i < n; ++i) {
      for (int j = 0; j < n; ++j) {
        if (i == j) continue;
        if (rand() % 2 == 0) {
          int a = 1 + rand() % 10; 
          int b = 1 + rand() % 10;
          g1.add_edge(i, j, a, b);
          g2.add_edge(i, j, a, b);
        } 
      }
    }

    int s = 0, t = n-1;
    tick();
    auto a = g1.min_cost_max_flow(s, t);
    cout << a.fst << " " << a.second << " " << tick() << endl;
    auto b = g2.min_cost_max_flow(s, t);
    cout << b.fst << " " << b.second << " " << tick() << endl;
    if (a.snd != b.snd) break;
  }
}

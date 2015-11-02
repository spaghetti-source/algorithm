//
// Single pair shortest path (Bidirectional dijkstra)
//
// Description:
//   For a pair of vertices, it finds a shortest path between
//   these two vertices.
//
// Algorithm:
//   Bidirectional dijkstra algorithm that performs Dijkstra 
//   algorithm from s and t simultaneously.
//   Usually, it is much faster than standard Dijkstra.
//
// Verified:
//   SPOJ SHPATH
// 
#include <iostream>
#include <queue>
#include <unordered_set>
#include <vector>
#include <string>
#include <map>
#include <algorithm>

using namespace std;

#define fst first
#define snd second

const int INF = 99999999;
struct graph {
  int n;
  struct edge { int src, dst, weight; };
  vector<vector<edge>> adj, rdj;
  graph(int n) : n(n), adj(n), rdj(n) { }
  void add_edge(int u, int v, int w) {
    adj[u].push_back({u, v, w});
    rdj[v].push_back({v, u, w});
  }

  int shortest_path(int s, int t) {
    if (s == t) return 0;
    vector<int> ds(n, INF), dt(n, INF);
    typedef pair<int, int> node;
    priority_queue<node, vector<node>, greater<node>> Qs, Qt;
    Qs.push({ds[s] = 0, s});
    Qt.push({dt[t] = 0, t});
    int mu = INF;
    while (!Qs.empty() && !Qt.empty()) {
      if (Qs.top().fst + Qt.top().fst >= mu) break;
      if (Qs.top().fst <= Qt.top().fst) {
        node x = Qs.top(); Qs.pop();
        if (ds[x.snd] > x.fst) continue;
        for (edge &e: adj[x.snd]) {
          if (ds[e.src] + e.weight < ds[e.dst]) {
            mu = min(mu, ds[e.src] + e.weight + dt[e.dst]);
            Qs.push({ds[e.dst] = ds[e.src] + e.weight, e.dst});
          }
        }
      } else {
        node x = Qt.top(); Qt.pop();
        if (dt[x.snd] > x.fst) continue;
        for (edge &e: rdj[x.snd]) {
          if (dt[e.src] + e.weight < dt[e.dst]) {
            mu = min(mu, dt[e.src] + e.weight + ds[e.dst]);
            Qt.push({dt[e.dst] = dt[e.src] + e.weight, e.dst});
          }
        }
      }
    }
    return mu;
  }
};


void solve() {
  int n; 
  scanf("%d", &n);
  graph g(n);
  map<string, int> id;
  for (int u = 0; u < n; ++u) {
    char name[1024];
    scanf("%s", name);
    id[name] = u;
    int k;
    scanf("%d", &k);
    for (int j = 0; j < k; ++j) {
      int v, d;
      scanf("%d %d", &v, &d);
      g.add_edge(u, v-1, d);
    }
  }
  int k;
  scanf("%d", &k);
  for (int i = 0; i < k; ++i) {
    char s[1024], t[1024];
    scanf("%s %s", s, t);
    printf("%d\n", g.shortest_path(id[s], id[t]));
  }
}

int main() {
  int ncase; scanf("%d", &ncase);
  for (int icase = 0; icase < ncase; ++icase) solve();
}

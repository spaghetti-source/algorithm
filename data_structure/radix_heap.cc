//
// Radix Heap
//
// Description:
//
//   A radix heap is a monotonic heap for keys of type unsigned int.
//   Here, a monotonic heap is a heap that allows push(key) only for 
//   key >= top(). It is useful in the Dijkstra shortest path algorithm
//   with integral edge weights.
//   see: https://en.wikipedia.org/wiki/Radix_heap
//
// Complexity:
//
//   O(|bit|) for all operations. Moreover,
//   O(1), amortized, for pop.
//
// Verified:
// 
//   SPOJ_SHPATH
//
// References:
// 
//   B. V. Cherkassky, A. V. Goldberg, C. Silverstein (1997):
//   "Buckets, Heaps, Lists and Monotone Priority Queues."
//   in Proceedings of the 8th Annual ACM-SIAM Symposium on 
//   Discrete Algorithms (SODA'97), pp. 83-92.
//   
#include <bits/stdc++.h>

using namespace std;

#define fst first
#define snd second
#define all(c) ((c).begin()), ((c).end())
#define TEST(s) if (!(s)) { cout << __LINE__ << " " << #s << endl; exit(-1); }

// Min Heap
template <class T>
struct RadixHeap {
  using uint = uint32_t; // fix bit size
  static int bsr(uint a) { return a ? 31 - __builtin_clz(a) : -1; }
  uint size, last;
  vector<pair<uint, T>> v[33];
  RadixHeap() : size(0), last(0) { }

  bool empty() const { return size == 0; }
  void aux(pair<uint, T> p) { v[bsr(p.fst^last)+1].push_back(p); }
  pair<uint, T> top() {
    if (v[0].empty()) {
      int i = 1;
      while (v[i].empty()) ++i;
      last = min_element(all(v[i]))->fst;
      for (auto p: v[i]) aux(p);
      v[i].clear();
    }
    return v[0].back();
  }
  void push(uint key, T value) { ++size; aux({key, value}); }
  void pop() { --size; top(); v[0].pop_back(); }
};

template <class Weight>
struct WeightedGraph {
  struct Edge {
    int to;
    Weight weight;
  };
  int n;
  vector<vector<Edge>> adj, rdj;
  WeightedGraph(int n) : n(n), adj(n), rdj(n) { }
  void addEdge(int u, int v, Weight w) {
    adj[u].push_back({v, w});
    rdj[v].push_back({u, w}); // can be omitted for the standard dijkstra
  }
};

template <class Weight>
struct ShortestPath { }; // omit

template <>
struct ShortestPath<int> {
  using Weight = int;
  WeightedGraph<Weight> g;

  vector<Weight> dist;
  vector<int> prev;
  ShortestPath(WeightedGraph<Weight> g) : g(g), dist(g.n), prev(g.n) { }

  void solve(int s) {
    fill(all(prev), -1);
    fill(all(dist), -1); dist[s] = 0;

    RadixHeap<int> que;
    que.push(0,s);
    while (!que.empty()) {
      auto d = que.top().fst;
      auto u = que.top().snd;
      que.pop();
      if (dist[u] < d) continue;
      for (auto e: g.adj[u]) {
        auto v = e.to;
        auto w = e.weight;
        if (dist[v] >= 0 && dist[v] <= dist[u]+w) continue;
        dist[v] = dist[u] + w;
        prev[v] = u;
        que.push(dist[v], v);
      }
    }
  }
  // Bidirectional Dijkstra
  int solve(int s, int t) {
    if (s == t) return dist[s] = 0;

    fill(all(dist), -1); dist[s] = 0;
    vector<Weight> drev(g.n, -1); drev[t] = 0;
    RadixHeap<int> qs, qt;
    qs.push(0,s); qt.push(0,t);
    int mu = -1;
    while (!qs.empty() && !qt.empty()) {
      if (mu >= 0 && qs.top().fst + qt.top().fst >= mu) break;
      if (qs.top().fst <= qt.top().fst) {
        auto d = qs.top().fst;
        auto u = qs.top().snd;
        qs.pop();
        if (dist[u] > d) continue;
        for (auto e: g.adj[u]) {
          auto v = e.to;
          auto w = e.weight;
          if (dist[v] >= 0 && dist[v] <= dist[u] + w) continue;
          dist[v] = dist[u] + w;
          qs.push(dist[v], v);
          if (drev[v] >= 0) {
            auto nu = dist[v] + drev[v];
            if (mu < 0 || mu > nu) mu = nu;
          }
        }
      } else {
        auto d = qt.top().fst;
        auto u = qt.top().snd;
        qt.pop();
        if (drev[u] > d) continue;
        for (auto e: g.rdj[u]) {
          auto v = e.to;
          auto w = e.weight;
          if (drev[v] >= 0 && drev[v] <= drev[u] + w) continue;
          drev[v] = drev[u] + w;
          qt.push(drev[v], v);
          if (dist[v] >= 0) {
            auto nu = dist[v] + drev[v];
            if (mu < 0 || mu > nu) mu = nu;
          }
        }
      }
    }
    return mu;
  }
};

void SPOJ_SHPATH() {
  int ncase; scanf("%d", &ncase);
  for (int icase = 0; icase < ncase; ++icase) {
    int n; 
    scanf("%d", &n);
    WeightedGraph<int> g(n);
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
        g.addEdge(u, v-1, d);
      }
    }
    ShortestPath<int> solver(g);
    int k;
    scanf("%d", &k);
    for (int i = 0; i < k; ++i) {
      char s[1024], t[1024];
      scanf("%s %s", s, t);
      printf("%d\n", solver.solve(id[s], id[t]));
    }
  }
}
int main() {
  SPOJ_SHPATH();
}

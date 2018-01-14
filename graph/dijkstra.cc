//
// Dijkstra's Single Source Shortest Path
//
// Description:
//
//   Dijkstra algorithm finds a single source shortest path on
//   a nonnegative weighted graph. 
//   It implements two algorithms with two data structures.
//   (1) standard dijkstra with standard heap
//   (2) bidirectional dijkstra with standard heap
//   (3) standard dijkstra with radix heap
//   (4) bidirectional dijkstra with radix heap
//   
//   For a simple test (see the code), we observe that 
//     Binomial,Unidirectional  14.74[s]
//     Binomial Bidirectional    0.52[s]
//     Radix,Unidirectional      4.90[s]
//     Radix,Bidirectional       0.31[s]
//
// Complexity:
//
//   O(m log n) for binomial heap
//   O(m) for radix heap
//
// Verify:
// 
//   SPOJ_SHPATH  (bidirectional search)
//
//
#include <bits/stdc++.h>
using namespace std;

#define fst first
#define snd second
#define all(c) ((c).begin()), ((c).end())

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
struct ShortestPath {
  WeightedGraph<Weight> g;
  ShortestPath(WeightedGraph<Weight> g) : g(g), dist(g.n), prev(g.n) { }

  vector<Weight> dist;
  vector<int> prev;
  void solve(int s) {
    prev.assign(g.n,-1);
    dist.assign(g.n,-1); dist[s] = 0;

    using Node = pair<Weight, int>;
    priority_queue<Node, vector<Node>, greater<Node>> que;
    que.push({0,s});
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
        que.push({dist[v], v});
      }
    }
  }
  int solve(int s, int t) {
    if (s == t) return dist[s] = 0;
    fill(all(dist), -1); dist[s] = 0;
    vector<Weight> drev(g.n, -1); drev[t] = 0;

    using Node = pair<Weight, int>;
    priority_queue<Node, vector<Node>, greater<Node>> qs, qt;
    qs.push({0,s}); qt.push({0,t});
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
          qs.push({dist[v], v});
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
          qt.push({drev[v], v});
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


//
// Dijkstra with Radix Heap
//
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

// === tick a time ===
#include <ctime>
double tick() {
  static clock_t oldtick;
  clock_t newtick = clock();
  double diff = 1.0*(newtick - oldtick) / CLOCKS_PER_SEC;
  oldtick = newtick;
  return diff;
}
int test() {
  int n = 5000; 
  WeightedGraph<int> g(n);
  for (int u = 0; u < n; ++u) {
    int d = 1 + rand() % 30;
    unordered_set<int> N;
    while (N.size() < d) {
      int v = rand() % n;
      if (u != v) N.insert(v);
    }
    for (auto v: N) {
      int w = 1 + rand() % 10;
      g.addEdge(u, v, w);
    }
  }
  ShortestPath<int> solver(g);
  double t = 0;
  for (int u = 0; u < n; ++u) {
    int v = rand() % n;
    tick();
    solver.solve(u);
    int b = solver.dist[v];
    t += tick();
  }
  cout << t << endl;
  return true;
}

int main() {
  test();
//  SPOJ_SHPATH();
  /*
  WeightedGraph<int> g(3);
  g.addEdge(0,1,10);
  g.addEdge(0,2,20);

  ShortestPath solver(g);
  solver.solve(0);
  for (int u = 0; u < g.n; ++u) {
    cout << solver.dist[u] << " ";
  }
  cout << endl;
  */
}

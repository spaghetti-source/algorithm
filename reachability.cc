#include <iostream>
#include <vector>
#include <cstdio>
#include <algorithm>
#include <queue>
#include <functional>
#include <unordered_map>
#include <random>
#include <map>

using namespace std;

#define fst first
#define snd second
#define all(c) ((c).begin()), ((c).end())
#define TEST(s) if (!(s)) { cout << __LINE__ << " " << #s << endl; exit(-1); }

const int INF = 99999999;
struct graph {
  int n;
  vector<vector<int>> adj, rdj;
  graph(int n) : n(n), adj(n), rdj(n) { }
  void add_edge(int u, int v) {
    adj[u].push_back(v);
  }
  typedef pair<int,int> interval;
  vector<interval> merge(vector<interval> x, vector<interval> y) {
    vector<interval> z;
    for (int i = 0, j = 0; i < x.size() || j < y.size(); ) {
      interval a = {INF, INF};
      if (i < x.size()) a.fst = min(a.fst, x[i].fst);
      if (j < y.size()) a.fst = min(a.fst, y[j].fst);
      a.snd = a.fst;
      while (1) {
        if      (i < x.size() && x[i].fst <= a.snd) a.snd = max(a.snd, x[i++].snd);
        else if (j < y.size() && y[j].fst <= a.snd) a.snd = max(a.snd, y[j++].snd);
        else break;
      }
      z.push_back(a); 
    }
    return z;
  }
  bool contain(vector<interval> x, int a) {
    if (a < x[0].fst) return false;
    int lo = 0, hi = x.size();
    while (hi - lo >= 2) { 
      int mi = (hi + lo) / 2;
      if (x[mi].fst <= a) lo = mi;
      else                hi = mi;
    }
    return a < x[lo].snd;
  }
  vector<int> rank, ord;
  vector<vector<pair<int,int>>> is;
  void build() {
    rank.assign(n, -1);
    ord.clear();
    function<void(int)> dfs = [&](int u) {
      rank[u] = 0;
      for (int v: adj[u]) 
        if (rank[v] < 0) dfs(v);
      rank[u] = ord.size();
      ord.push_back(u);
    }; 
    for (int u = 0; u < n; ++u) 
      if (rank[u] < 0) dfs(u);

    is.resize(n);
    for (int u = 0; u < n; ++u) 
      is[u] = { {rank[u], rank[u]+1} };
    
    while (1) {
      bool updated = false;
      for (int u: ord) {
        vector<interval> next = is[u];
        for (int v: adj[u]) 
          next = merge(next, is[v]);
        if (is[u] != next) {
          is[u] = next;
          updated = true;
        }
      }
      if (!updated) break;
    }
  }
  bool is_reachable(int u, int v) {
    return contain(is[u], rank[v]);
  }

  // for verification
  bool is_reachable_n(int s, int t) {
    vector<int> seen(n, 0);
    function<bool(int)> dfs = [&](int u) {
      seen[u] = true;
      if (u == t) return true;
      for (int v: adj[u]) 
        if (!seen[v] && dfs(v)) return true;
      return false;
    };
    return dfs(s);
  }
};
graph random_graph(int n, int d) {
  graph g(n);
  for (int i = 0; i < n; ++i) {
    for (int k = 0; k < d; ++k) {
      int j = rand() % n;
      if (i > j) g.add_edge(i, j);
      else       g.add_edge(j, i);
    }
  }
  return g;
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
  int n = 100000;
  graph g = random_graph(n, 1);;
  g.build();
  cout << n << endl;
  double t1 = 0, t2 = 0;
  for (int i = 0; i < n; ++i) {
    int u = rand() % n;
    int v = rand() % n;
    tick();
    int a = g.is_reachable_n(u, v);
    t1 += tick();
    int b = g.is_reachable(u, v);
    t2 += tick();
    if (a != b) cout << "*" << endl;
  }
  cout << t1 << " " << t2 << endl;
}

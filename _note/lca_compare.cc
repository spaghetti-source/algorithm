// 
// Least Common Ancestor (Heavy-light decomposition)
//
// 1) Heavy-Light decomposition
// 2) Doubling
// 3) Euler-tour + RMQ
//
// Use 1 or 3.
//

#include <iostream>
#include <vector>
#include <cstdio>
#include <algorithm>
#include <functional>
#include <numeric>

using namespace std;

#define fst first
#define snd second
#define all(c) ((c).begin()), ((c).end())


// === tick a time ===
#include <ctime>
double tick() {
  static clock_t oldtick;
  clock_t newtick = clock();
  double diff = 1.0*(newtick - oldtick) / CLOCKS_PER_SEC;
  oldtick = newtick;
  return diff;
}


namespace HEAVY_LIGHT {
struct tree {
  int n;
  vector<vector<int>> adj;
  tree(int n) : n(n), adj(n) { }
  void add_edge(int s, int t) {
    adj[s].push_back(t);
    adj[t].push_back(s);
  }
  vector<int> size, depth, head, parent;
  void rootify(int r) {
    size = depth = parent = head = vector<int>(n, -1);
    function<int (int,int)> dfs = [&](int u, int p) {
      parent[u] = p;
      depth[u] = depth[p]+1;
      for (int v: adj[u]) 
        if (v != p) size[u] += dfs(v, u);
      return ++size[u];
    }; dfs(r, r);
    function<void (int,int)> dec = [&](int u, int s) {
      head[u] = s;
      int z = -1;
      for (int v: adj[u]) 
        if (head[v] < 0 && (z < 0 || size[z] < size[v])) z = v;
      for (int v: adj[u]) 
        if (head[v] < 0) dec(v, v == z ? s : v);
    }; dec(r, r);
  }
  int lca(int u, int v) {
    while (head[u] != head[v]) 
      if (depth[head[u]] < depth[head[v]]) v = parent[head[v]];
      else                                 u = parent[head[u]];
    return depth[u] < depth[v] ? u : v;
  }
};
}

// Doubling (memoise 2^h-ancestors)
// Not fast. No need to use
namespace DOUBLING {
struct tree {
  int n, logn;
  vector<vector<int>> adj;
  tree(int n) : n(n), adj(n) { 
    logn = 1;
    for (int k = 1; k < n; k *= 2) ++logn;
  }
  void add_edge(int s, int t) {
    adj[s].push_back(t);
    adj[t].push_back(s);
  }
  vector<vector<int>> parent;
  vector<int> rank, depth;
  void rootify(int r) {
    parent.assign(logn, vector<int>(n, r));
    rank.resize(n);
    depth.assign(n, 0);
    int id = 0;
    function<void (int,int)> dfs = [&](int u, int p) {
      rank[u] = id++;
      parent[0][u] = p;
      for (int i = 0; i+1 < logn; ++i) 
        parent[i+1][u] = parent[i][parent[i][u]];
      for (int v: adj[u]) 
        if (v != p) { depth[v] = depth[u]+1; dfs(v, u); }
    }; dfs(r, r);
  }

  int lca(int u, int v) {
    if (depth[u] > depth[v]) swap(u, v);
    for (int i = depth[v]-depth[u], k = 0; i; i /= 2) {
      if (i & 1) v = parent[k][v];
      ++k;
    }
    if (u == v) return u;
    for (int i = logn-1; i >= 0; --i) {
      if (parent[i][u] != parent[i][v]) {
        u = parent[i][u];
        v = parent[i][v];
      }
    }
    return parent[0][u];
  }
};
}

namespace EULER_TOUR_SPARSE_TABLE {
struct tree {
  int n;
  vector<vector<int>> adj;
  tree(int n) : n(n), adj(n) { }
  void add_edge(int s, int t) {
    adj[s].push_back(t);
    adj[t].push_back(s);
  }
  vector<int> pos, tour, depth;
  vector<vector<int>> table;
  int argmin(int i, int j) { return depth[i] < depth[j] ? i : j; }
  void rootify(int r) {
    pos.resize(n);
    tick();
    function<void (int,int,int)> dfs = [&](int u, int p, int d) {
      pos[u] = depth.size();
      tour.push_back(u);
      depth.push_back(d);
      for (int v: adj[u]) {
        if (v != p) {
          dfs(v, u, d+1);
          tour.push_back(u);
          depth.push_back(d);
        }
      }
    }; dfs(r, r, 0);
    int logn = sizeof(int)*__CHAR_BIT__-1-__builtin_clz(tour.size()); // log2
    table.resize(logn+1, vector<int>(tour.size()));
    iota(all(table[0]), 0);
    for (int h = 0; h < logn; ++h) 
      for (int i = 0; i+(1<<h) < tour.size(); ++i)
        table[h+1][i] = argmin(table[h][i], table[h][i+(1<<h)]);
  }
  int lca(int u, int v) {
    int i = pos[u], j = pos[v]; if (i > j) swap(i, j);
    int h = sizeof(int)*__CHAR_BIT__-1-__builtin_clz(j-i); // = log2
    return i == j ? u : tour[argmin(table[h][i], table[h][j-(1<<h)])];
  }
};
};

int main() {
  int n = 1000000;
  HEAVY_LIGHT::tree T1(n);
  EULER_TOUR_SPARSE_TABLE::tree T2(n);
  for (int i = 1; i < n; ++i) {
    int u = rand() % i;
    T1.add_edge(u, i);
    T2.add_edge(u, i);
  }
  int count = 0;
  tick();
  T1.rootify(0);
  cout << tick() << " ";
  for (int i = 0; i < n; ++i) {
    int u = rand() % n;
    int v = rand() % n;
    count += T1.lca(u,v);
  }
  cout << tick() << endl;

  tick();
  T2.rootify(0);
  cout << tick() << " ";
  for (int i = 0; i < n; ++i) {
    int u = rand() % n;
    int v = rand() % n;
    count += T2.lca(u,v);
  }
  cout << tick() << endl;

  if (count == 0) cout << "lucky" << endl;
  for (int i = 0; i < n; ++i) {
    int u = rand() % n;
    int v = rand() % n;
    if (T1.lca(u,v) != T2.lca(u,v)) {
      cout << "wrong" << endl;
      return -1;
    }
  }
}

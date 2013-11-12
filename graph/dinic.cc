//
// Maximum Flow (Dinic)
//
#include <iostream>
#include <vector>
#include <queue>
#include <cstdio>
#include <algorithm>

using namespace std;

int n;
struct edge {
  int src, dst;
  double residue;
  edge *rev;
};
vector<vector<edge*>> adj;
void add_edge(int u, int v, double c) {
  adj[u].push_back(new edge{u, v, c});
  adj[v].push_back(new edge{v, u, 0}); // for directed
  // adj[v].push_back(new edge{v, u, c}); // for undirected
  adj[u].back()->rev = adj[v].back();
  adj[v].back()->rev = adj[u].back();
}
vector<int> level, visiting;
double augment(int u, int t, double cur) {
  if (u == t || cur <= 0) return cur;
  if (visiting[u]) return 0;
  visiting[u] = true;
  for (auto e: adj[u]) {
    if (level[e->dst] > level[u]) {
      double f = augment(e->dst, t, min(cur, e->residue));
      if (f > 0) {
        e->residue -= f;
        e->rev->residue += f;
        visiting[u] = false;
        return f;
      }
    }
  }
  return 0;
}
double maximum_flow(int s, int t) {
  double total = 0;
  while (1) {
    level.assign(n, -1);
    level[s] = 0;
    queue<int> Q;
    Q.push(s);
    
    while (!Q.empty()) {
      int u = Q.front();
      if (u == t) break;
      Q.pop();
      for (auto e: adj[u]) {
        if (e->residue > 0 && level[e->dst] == -1) {
          Q.push(e->dst);
          level[e->dst] = level[e->src] + 1;
        }
      }
    }
    double prev = total;
    visiting.assign(n, 0);
    while (1) {
      double f = augment(s, t, 1.0/0.0);
      if (f == 0) break;
      total += f;
    }
    if (prev == total) break;
  }
  return total;
}

// PKU, Power Network
int main() {
  for (int np, nc, m; scanf("%d %d %d %d", &n, &np, &nc, &m) == 4; ) {
    int s = n, t = n+1;
    n += 2;
    adj.clear();
    adj.resize(n);
    for (int i = 0; i < m; ++i){
      int src, dst, cap;
      scanf(" (%d ,%d )%d", &src, &dst, &cap);
      add_edge(src, dst, cap);
    }
    for (int i = 0; i < np; ++i){
      int dst, cap;
      scanf(" ( %d )%d", &dst, &cap);
      add_edge(s, dst, cap);
    }
    for(int i = 0; i < nc; ++i){
      int src, cap;
      scanf(" ( %d )%d ", &src, &cap);
      add_edge(src, t, cap);
    }
    double w = maximum_flow(s, t);
    printf("%d\n", (int)w);
  }
}

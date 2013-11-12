//
// Dinic's maximum flow
//
//
// Description:
//   Compute the maximum flow.
//
//
// Algorithm:
//   The algorithm iterates following procedures:
//     (1) BFS from the source to get the distance to the sink.
//         If not reachable, there are no augment path hence break.
//     (2) Find vertex disjoint shortest augment paths by DFS.
//   Since the length of augmenting path increases at least 1, for each outer-loop,
//   the number of outer-loop is O(n). For each inner-loop, it takes O(nm) because
//   it consists of one BFS and O(n) DFSs. Therefore the whole complexity is O(n^2 m).
//
//
// Complexity:
//   O(n^2 m) time
//
//
// Verified:
//   PKU 1459, Power Network
//
//
// References:
//   E. A. Dinic (1970):
//   Algorithm for solution of a problem of maximum flow in networks with power estimation.
//   Soviet Mathematics Doklady, vol. 11, pp. 1277-1280.


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
  int rev;
};
vector< vector<edge> > adj;
void add_edge(int u, int v, double c) {
  adj[u].push_back((edge){u, v, c, adj[v].size()});
  adj[v].push_back((edge){v, u, 0, adj[u].size()-1});
}
vector<int> level, iter;
double augment(int u, int t, double cur) {
  if (u == t) return cur;
  for (int &i = iter[u]; i < adj[u].size(); ++i) {
    edge &e = adj[u][i];
    if (e.residue > 0 && level[u] < level[e.dst]) {
      double f = augment(e.dst, t, min(cur, e.residue));
      if (f > 0) {
        e.residue -= f;
        adj[e.dst][e.rev].residue += f;
        return f;
      }
    }
  }
  return 0;
}
int bfs(int s, int t) {
  level.assign(n, -1);
  level[s] = 0;
  queue<int> Q;
  Q.push(s);
  while (!Q.empty()) {
    int u = Q.front(); Q.pop();
    for (int i = 0; i < adj[u].size(); ++i) {
      edge &e = adj[u][i];
      if (e.residue > 0 && level[e.dst] < 0) {
        Q.push(e.dst);
        level[e.dst] = level[u] + 1;
      }
    }
  }
  return level[t];
}
double maximum_flow(int s, int t) {
  double flow = 0;
  while (1) {
    if (bfs(s, t) < 0) return flow;
    iter.assign(n, 0);
    for (double f; (f = augment(s, t, 1.0/0.0)) > 0; ) 
      flow += f;
  }
}

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

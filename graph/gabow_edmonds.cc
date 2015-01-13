//
// General Graph Matching (Gabow-Edmonds)
//
// Description:
//   It computes a maximum matching in a general graph.
//
// Algorithm:
//   Gabow's simplified version of Edmonds' blossom algorithm.
// 
// Comlexity:
//   O(n^3)
//
// Verified:
//   LA3820, LA4130
//
// References:
//   H.Gabow (1976):
//   An efficient implementation of Edmonds' algorithm for maximum matching on graphs.
//   Journal of the ACM, vol.23, no.2, pp.221-234.
//

#include <cstring>
#include <vector>
#include <queue>
#include <iostream>
#include <cstdlib>
#include <cstdio>
#include <algorithm>
#include <set>

using namespace std;

struct graph {
  int n;
  vector<vector<int>> adj;
  graph(int n) : n(n), adj(n) { };
  void add_edge(int x, int y) {
    adj[x].push_back(y);
    adj[y].push_back(x);
  }
  queue<int> Q;
  vector<int> label, mate, cycle;
  void rematch(int x, int y) {
    int m = mate[x]; mate[x] = y;
    if (mate[m] == x) {
      if (label[x] < n) {
        rematch(mate[m] = label[x], m);
      } else {
        int s = (label[x]-n)/n, t = (label[x]-n)%n;
        rematch(s, t); rematch(t, s);
      }
    }
  }
  void traverse(int x) {
    vector<int> save = mate;
    rematch(x, x);
    for (int u = 0; u < n; ++u) 
      if (mate[u] != save[u]) cycle[u] ^= 1;
    save.swap(mate);
  }
  void relabel(int x, int y) {
    cycle = vector<int>(n, 0); 
    traverse(x);
    traverse(y);
    for (int u = 0; u < n; ++u) {
      if (!cycle[u] || label[u] >= 0) continue;
      label[u] = n+x+y*n;
      Q.push(u);
    }
  }
  int augment(int r) {
    label.assign(n, -2);
    label[r] = -1;
    Q = queue<int>(); Q.push(r);
    while (!Q.empty()) {
      int x = Q.front(); Q.pop();
      for (int y: adj[x]) {
        if (mate[y] < 0 && r != y) {
          rematch(mate[y] = x, y); return 1;
        } else if (label[y] >= -1) {
          relabel(x, y);
        } else if (label[mate[y]] < -1) {
          label[mate[y]] = x;
          Q.push(mate[y]);
        }
      }
    }
    return 0;
  }
  int maximum_matching() {
    mate.assign(n, -2);
    int matching = 0;
    for (int u = 0; u < n; ++u) 
      if (mate[u] < 0) matching += augment(u);
    return matching;
  }
};

int doit() {
  int n, m; 
  scanf("%d %d", &n, &m);

  vector<int> v(n);
  for (int i = 0; i < n; ++i) {
    scanf("%d", &v[i]);
  }
  set<int> S;
  for (int i = 0; i < m; ++i) {
    int x;
    scanf("%d", &x);
    S.insert(x);
  }

  graph g(n);
  for (int i = 0; i < n; ++i) {
    for (int j = i+1; j < n; ++j) {
      if (S.count(v[i] + v[j])) 
        g.add_edge(i, j);
    }
  }
  return g.maximum_matching();
}
int main() {
  doit();
}

// 
// Maximum Cut (approximation)
//
//
// Description:
//   Given a weighted undirected graph G = (V, E), w: E -> R+.
//   Find a maximum cut, i.e., a subset S such that w(S,V-S) is maximum.
//   This problem is known to be NP-hard. 
//
//
// Algorithm:
//   Sahni and Gonzalez's greedy algorithm.
//   The algorithm has only 1/2 approximation ratio, 
//   but this has one of the best performance in practice.
//
//
// Complexity:
//   O(m log n).
//   
//
// Reference
//   S. Kahruman, E. Kolotoglu, S. Butenko, and I. V. Hicks (2007):
//   On greedy construction heuristics for the MAX-CUT problem.
//   International Journal on Computational Science and Engineering,
//   vol. 3, no. 3, pp. 211-218.
//
//   S. Sahni and T. Gonzales (1976):
//   P-complete approximation problems.
//   Journal of the ACM, 
//   vol. 23, no. 3, pp. 555-565.
//

#include <iostream>
#include <vector>
#include <cstdio>
#include <cstdlib>
#include <queue>
#include <algorithm>

using namespace std;

struct maximum_cut {
  int n;
  struct edge {
    int src, dst;
    double weight;
  };
  vector<edge> edges;
  vector<vector<edge>> adj;

  void add_edge(int s, int t, double w = 1) {
    edges.push_back({s, t, w});
  }
  void make_graph() {
    n = 0;
    for (auto e: edges) 
      n = max(n, max(e.src, e.dst)+1);
    adj.resize(n);
    for (auto e: edges) {
      adj[e.src].push_back(e);
      swap(e.src, e.dst);
      adj[e.src].push_back(e);
    }
  }

  priority_queue<pair<double, int>> que; // (score,vertex)
  vector<double> wS, wT;
  vector<bool> S, T;
  double score(int u) { return abs(wS[u] - wT[u]); }
  void insert(int u, vector<bool> &U, vector<double> &wU) {
    U[u] = 1;
    for (auto e: adj[u]) {
      int v = e.dst;
      if (S[v] || T[v]) continue;
      wU[v] += e.weight;
      que.push({score(v), v});
    }
  }
  double solve() {
    while (!que.empty()) que.pop();
    for (int i = 0; i < n; ++i) 
      que.push( {0, i} );

    S.assign(n, 0);  T.assign(n, 0);
    wS.assign(n, 0); wT.assign(n, 0);

    edge emax = *max_element(edges.begin(), edges.end(),
        [](edge a, edge b) { return a.weight < b.weight; });

    double cut = emax.weight;
    insert(emax.src, S, wS);
    insert(emax.dst, T, wT);

    while (!que.empty()) {
      int u = que.top().second;
      double d = que.top().first;
      que.pop();
      if (S[u] || T[u]) continue; // already fixed
      if (abs(score(u) - d) > 1e-8) continue; // already score modified
      cut += max(wS[u], wT[u]);
      if (wS[u] < wT[u]) insert(u, S, wS);
      else               insert(u, T, wT);
    }
    return cut;
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
  maximum_cut solver;

  int n = 100000, m = 1000000;
  for (int i = 0; i < m; ++i) {
    int s, t;
    do {
      s = rand() % n;
      t = rand() % n;
    } while (s == t);
    solver.add_edge(s, t, rand() * n / (1.0 + RAND_MAX));
  }

  solver.make_graph();

  tick();
  printf("%f", solver.solve());
  printf(", %f[s]\n", tick());
}

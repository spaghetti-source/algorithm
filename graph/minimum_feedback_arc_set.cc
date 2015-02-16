//
// Minimum Feedback Arc Set
//
// Description:
//   We are given an undirected weighted graph G = (V, E; w).
//   A set of edges F is said to be a feedback arc set if
//   (V, E-F) is has no cycles.
//   The problem is to find a minimum weight feedback arc set.
//
// Algorithm:
//   Min FAS problem is equivalent to find an ordering p 
//   that minimizes the total weight of inconsistence edges:
//     \sum { w(u,v) : u <_p v }.
//   This gives a dynamic programming solution.
//   For a subset S of V, let f(S) be the value of MFAS for G(S).
//   Then we have
//     f(S) = min { f(S-u) + sum { w(u,v) : v in S } }.
//   This means that the optimal solution for S is obtained by
//   adding largest element u to the optimal solution S-u.
//   
// Complexity:
//   O(m 2^n).
//
// References:
//   V. Raman and S. Saurabh (2007):
//   Improved fixed parameter tractable algorithms for two edge
//   problems: MAXCUT and MAXDAG.
//   Information Processing Letters, vol. 104, pp. 65--72.

#include <iostream>
#include <vector>
#include <cstdio>
#include <algorithm>
#include <functional>

using namespace std;

#define fst first
#define snd second
#define all(c) ((c).begin()), ((c).end())


template <class T>
struct graph {
  struct edge {
    int src, dst;
    T weight;
  };
  int n;
  vector<vector<edge>> adj;
  T inf;
  graph(int n) : n(n), adj(n), inf(0) { }
  void add_edge(int src, int dst, T weight) {
    adj[src].push_back({src, dst, weight});
    inf += e.weight;
  }

  T min_feedback_arc_set() {
    vector<T> f(1<<n, inf);
    f[0] = 0;
    for (int S = 0; S < (1<<n); ++S) {
      for (int u = 0; u < n; ++u) {
        if (S & (1<<u)) continue;
        T w = 0;
        for (edge e: adj[u]) 
          if (S & (1<<e.dst)) 
            w += e.weight;
        f[S|(1<<u)] = min(f[S|(1<<u)], f[S] + w);
      }
    }
    return f[(1<<n)-1];
  }


  // for verification
  T min_feedback_arc_set_naive() {
    T opt = inf;
    vector<int> pi(n); iota(all(pi), 0);
    do {
      T ans = 0;
      for (int u = 0; u < n; ++u)
        for (edge e: adj[u]) 
          if (pi[e.src] < pi[e.dst]) 
            ans += e.weight;
      opt = min(opt, ans);
    } while (next_permutation(all(pi)));
    return opt;
  }
};

int main() {
  for (int seed = 0; seed < 1000; ++seed) {
    srand(seed);
    int n = 7;
    graph<int> g(n);
    for (int i = 0; i < n; ++i) {
      for (int j = 0; j < n; ++j) {
        if (i == j) continue;
        int w = rand() % 10;
        g.add_edge(i, j, w);
//        cout << i << " " << j << " " << w << endl;
      }
    }
    int a = g.min_feedback_arc_set();
    int b = g.min_feedback_arc_set_naive();
    if (a != b) {
      cout << seed << endl;
      cout << "DP = " << a << endl;
      cout << "naive = " << b << endl;
      break;
    }
  }
}

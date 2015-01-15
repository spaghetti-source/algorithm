// Minimum Arborescence (Chu-Liu/Edmonds)
//
// Description:
//   Let G = (V, E) be a weighted directed graph.
//   For a vertex r, an edge-set T is called r-arborescense if
//     (1) T is a spanning tree (with forgetting directions),
//     (2) for each u in V, indeg_T(u) <= 1, indeg_T(r) = 0.
//   The program finds the minimum weight of r-arborescence.
//
//
// Algorithm:
//   Chu-Liu/Edmonds' recursive shrinking.
//   At first, it finds a minimum incomming edge for each v in V.
//   Then, if it forms a arborescence, it is a solution,
//   and otherwise, it contracts a cycle and iterates the procedure.
//
//
// Complexity: 
//   O(mn)
//
//
// Remark:
//   More efficient (but long) version is also available.
//
//
// References: 
//   Y. J. Chu and T. H. Liu (1965): 
//   On the shortest arborescence of a directed graph,
//   Science Sinica, vol. 14, pp. 1396--1400.
//
//   J. Edmonds (1967): 
//   Optimum branchings
//   Journal on Research of the National Bureau of Standards, 71B,
//   pp. 233--240.
//

#include <iostream>
#include <vector>
#include <cstdio>
#include <cstdlib>
#include <map>
#include <cmath>
#include <cstring>
#include <functional>
#include <algorithm>
#include <unordered_map>

using namespace std;

#define fst first
#define snd second
#define all(c) (c).begin(), (c).end()

struct graph {
  int n;
  graph(int n) : n(n) { }
  struct edge {
    int src, dst;
    int weight;
  };
  vector<edge> edges;
  void add_edge(int u, int v, int w) {
    edges.push_back({u, v, w});
  }
  int arborescence(int r) {
    int N = n;
    for (int res = 0; ;) {
      vector<edge> in(N, {-1,-1,(int)INF});
      vector<int> C(N, -1);
      for (auto e: edges) // cheapest comming edges
        if (in[e.dst].weight > e.weight) 
          in[e.dst] = e;
      in[r] = {r, r, 0};

      for (int u = 0; u < N; ++u) { // no comming edge ==> no aborescense
        if (in[u].src < 0) return -1;
        res += in[u].weight;
      }
      vector<int> mark(N, -1); // contract cycles
      int index = 0;
      for (int i = 0; i < N; ++i) {
        if (mark[i] != -1) continue;
        int u = i;
        while (mark[u] == -1) {
          mark[u] = i;
          u = in[u].src;
        }
        if (mark[u] != i || u == r) continue;
        for (int v = in[u].src; u != v; v = in[v].src) C[v] = index;
        C[u] = index++;
      }
      if (index == 0) return res; // found arborescence
      for (int i = 0; i < N; ++i) // contract
        if (C[i] == -1) C[i] = index++;

      vector<edge> next;
      for (auto &e: edges) 
        if (C[e.src] != C[e.dst] && C[e.dst] != C[r]) 
          next.push_back({C[e.src], C[e.dst], e.weight - in[e.dst].weight});
      edges.swap(next);
      N = index; r = C[r];
    }
  }
};


#include <ctime>


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

  for (int q = 0; q <= 10; ++q) {
  srand( q );

  int n = 1000;
  graph g(n);

  for (int i = 0; i < n; ++i)
    for (int k = 0; k < 10; ++k) {
      int j;
      do {
        j = rand() % n;
      } while (i == j);
      int w = 1 + rand() % 100;
      g.add_edge(i, j, w);
    }

  tick();
  double b = g.arborescence(0);
  cout << b << " " << tick() << endl;

  cout << abs(b - a) << endl;

  }
}

//
// Enumerating maximal cliques (Bron and Kerbosch)
//
// Description:
//   A clique is a fully-connected subset of vertices.
//   This algorithm finds all maximal cliques.
//   It can also be used for enumerating maximal independent sets.
// 
// Algorithm:
//   Bron and Kerbosch's exhaustive search.
//
// Complexity:
//   O(3^{n/3}) without output, which is practical for n <= 64.
//
// References: 
//   C. Bron and J. Kerbosch (1973):
//   Algorithm 457: finding all cliques of an undirected graph.
//   Communication of the ACM, vol. 16, no. 9, pp. 575-577.
//
// Comment:
//   We have also implemented the Tomita's pivotting strategy
//   and/or degeneracy ordering; 
//   however, this version is the most efficient.
//

#include <iostream>
#include <vector>
#include <cstdio>
#include <algorithm>
#include <functional>
#include <unordered_map>

using namespace std;

#define fst first
#define snd second
#define all(c) ((c).begin()), ((c).end())

struct small_graph { // assume: |V| <= 64
  typedef unsigned long long state;
  int n;
  unordered_map<state, state> N; // neighbor of i
  small_graph(int n = 0) : n(n) { }
  void add_edge(int src, int dst) {
    n = max(n, max(src, dst)+1);
    N[1ull<<src] |= (1ull<<dst);
    N[1ull<<dst] |= (1ull<<src);
  }
  int maximal_cliques() {
    // enumerate cliques that include R, exclude X, arbitrary P.
    function<int (state,state,state)> rec = [&](state R, state P, state X) {
      if (!(P | X)) { // R is a bitset of maximal clique
        for (int i = 0; i < n; ++i, R >>= 1) 
          if (R & 1) cout << i << " ";
        cout << endl;
        return 1;
      }
      int count = 0;
      state u = (P | X); u = u & ~(u-1);
      while (P & ~N[u]) {
        state v = P & ~(P-1);
        count += rec(R | v, P & N[v], X & N[v]);
        P ^= v; X |= v;
      }
      return count;
    };
    return rec(0, (1ull<<n)-1, 0);
  }
  // independent set = clique in the complement graph
  int maximal_independent_set() {
    function<int (state,state,state)> rec = [&](state R, state P, state X) {
      if (!(P | X)) { // R is a bitset of maximal independent set
        for (int i = 0; i < n; ++i, R >>= 1) 
          if (R & 1) cout << i << " ";
        cout << endl;
        return 1; 
      }
      int count = 0;
      state u = (P | X); u = u & ~(u-1);
      while (P & (N[u]|u)) {
        state v = P & ~(P-1); 
        count += rec(R | v, P & ~(N[v]|v), X & ~(N[v]|v));
        P ^= v; X |= v;
      }
      return count;
    };
    return rec(0, (1ull<<n)-1, 0);
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
  small_graph g;
  g.add_edge(0,1);
  g.add_edge(1,2);
  g.add_edge(2,3);
  g.add_edge(3,0);
  g.add_edge(3,4);
  g.add_edge(4,5);
  g.add_edge(5,0);
  g.maximal_independent_set();
  /*
  int n = 62;
  small_graph g;
  for (int i = 0; i < n; ++i) {
    for (int k = 0; k < n/2; ++k) {
      while (1) {
        int j = rand() % n;
        if (i != j) {
          g.add_edge(i, j);
          break;
        }
      }
    }
  }
  tick();
  cout << g.maximum_clique() << endl;
  cout << tick() << endl;
  */
}

// 
// Hamilton Cycle for graphs with Ore condition
//
// Description:
//   When an undirected graph satisfies the Ore condition,
//     d(u) + d(v) >= n, (u,v) not in E.
//   it has a Hamiltonian cycle.
//
// Algorithm:
//   Palmer's construction.
//   Arrange vertices on n-cycle (regardless of the edges).
//   For two vertices, i and i+1, with (i,i+1) not in E,
//   we find vertex j such that (i,j) and (j+1,i+1) in E.
//   We rearrange cycles by 
//     i -> j -> j-1 -> ... -> i -> j+1,
//   which increases the length of a span.
//
// Complexity:
//   O(n^2)
//
// References:
//   E. M. Palmer (1997): The hidden algorithm of Ore's theorem on Hamiltonian cycles.
//   Computers & Mathematics with Applications, vol.34, no.11, pp.113--119.
//   

#include <iostream>
#include <vector>
#include <cstdio>
#include <algorithm>
#include <unordered_set>
#include <numeric>
#include <functional>
#include <random>

using namespace std;

#define fst first
#define snd second
#define all(c) ((c).begin()), ((c).end())

struct graph {
  int n;
  vector<unordered_set<int>> adj;
  graph(int n) : n(n), adj(n) { }
  void add_edge(int src, int dst) {
    adj[src].insert(dst);
    adj[dst].insert(src);
  }
  vector<int> hamilton_cycle() {
    auto X = [&](int i) { return i < n ? i : i - n; }; // faster than mod
    vector<int> cycle(n);
    iota(all(cycle), 0);
    while (1) { 
      bool updated = false;
      for (int i = 0; i < n; ++i) {
        if (adj[cycle[i]].count(cycle[X(i+1)])) continue;
        for (int j = i+2; j < i+n; ++j) {
          if (adj[cycle[i]].count(cycle[X(j)]) &&
              adj[cycle[X(i+1)]].count(cycle[X(j+1)])) {
            for (int k = i+1, l = j; k < l; ++k, --l)
              swap(cycle[X(k)], cycle[X(l)]);
            updated = true;
            break;
          }
        }
      }
      if (!updated) break;
    }
    return cycle;
  }
};


template <class T>
ostream &operator<<(ostream &os, const vector<T> &v) {
  os << "[";
  for (int i = 0; i < v.size(); os << v[i++]) 
    if (i > 0) os << " ";
  os << "]";
  return os;
}
template <class T>
ostream &operator<<(ostream &os, const vector<vector<T>> &v) {
  os << "[";
  for (int i = 0; i < v.size(); os << v[i++]) 
    if (i > 0) os << endl << " ";
  os << "]";
  return os;
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
  for (int seed = 0; seed < 100; ++seed) {
    cout << "seed = " << seed << endl;
    int n = 500;
    graph g(n);
    vector<int> ord(n);
    iota(all(ord), 0);
    //shuffle(all(ord), default_random_engine(time(0)));
    auto engine = default_random_engine(seed);
    shuffle(all(ord), engine);
    for (int i: ord) for (int j: ord) {
      if (i == j) continue;
      if (engine() % (n/2)) g.add_edge(i, j);
    }
    while (1) {
      bool finish = true;
      for (int i: ord) for (int j: ord) {
        if (i == j) continue;
        if (g.adj[i].count(j)) continue;
        if (g.adj[i].size() + g.adj[j].size() < n) {
          g.add_edge(i, j);
          finish = false;
        }
      }
      if (finish) break;
    }

    tick();
    auto cycle = g.hamilton_cycle();
    cout << tick() << endl;
   // cout << cycle << endl;
    for (int i = 0; i < n; ++i) {
      if (!g.adj[cycle[i]].count(cycle[(i+1)%n])) {
        cout << "seed = " << seed << ": something wrong\n";
      }
    }
  }
}

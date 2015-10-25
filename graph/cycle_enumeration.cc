//
// Hawick and James' cycle enumeration
//
// Description:
//   For a directed graph, it enumerates all cycles
//
// Algorithm:
//   Hawick and James's implementation of Johnson algorithm.
//
// Complexity
//   O(n+m) average time for each cycles, 
//   O(n+m) space.
//
// References:
//   K. A. Hawick and H. A. James (2007): 
//   Enumerating circuits and loops in graphs with self-arcs and multiple-arcs.
//   in Proceedings of International Conference on Foundations of Computer Science,
//   pp.14--20.
//
//   D. B. Johnson (1975): 
//   Finding all the elementary circuits of a directed graph.
//   SIAM Journal on Computing, vol.4, pp.77--84.
//
#include <iostream>
#include <vector>
#include <cstdio>
#include <algorithm>
#include <unordered_set>
#include <functional>

using namespace std;

#define fst first
#define snd second
#define all(c) ((c).begin()), ((c).end())
#define TEST(s) if (!(s)) { cout << __LINE__ << " " << #s << endl; exit(-1); }

struct graph {
  int n;
  vector<vector<int>> adj;
  graph(int n) : n(n), adj(n) { }
  void add_edge(int u, int v) {
    adj[u].push_back(v);
  }

  void all_cycles() {
    vector<int> S;
    for (int s = 0; s < n; ++s) {
      vector<bool> blocked(n);
      vector<unordered_set<int>> B(n);

      function<void(int)> unblock = [&](int u) {
        blocked[u] = false;
        for (int v: B[u]) 
          if (blocked[v]) unblock(v);
        B[u].clear();
      };
      function<bool(int)> rec = [&](int u) {
        bool is_cycle = false;
        blocked[u] = true;
        S.push_back(u);
        for (int v: adj[u]) {
          if (v <  s) continue;
          if (v == s) {
            is_cycle = true; // S forms a cycle
            cout << "found cycle" << endl;
            for (auto w: S) cout << " " << w;
            cout << endl;
          } else if (!blocked[v] && rec(v)) is_cycle = true;
          if (is_cycle) unblock(u);
          else {
            for (int v: adj[u]) {
              if (v < s) continue;
              if (!B[v].count(u)) B[v].insert(u);
            }
          }
        }
        S.pop_back();
        return is_cycle;
      }; rec(s);
    }
  }
};

int main() {
  int n = 4;
  graph g(n);
  g.add_edge(0,1);
  g.add_edge(1,2);
  g.add_edge(2,0);
  g.add_edge(1,3);
  g.add_edge(3,0);
  g.add_edge(2,3);

  g.all_cycles();
}

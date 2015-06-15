//
// Claw-free graph recognition
//
// Description:
//   A graph is claw-free if it contains no claws, i.e., 
//     o--o--o
//        | 
//        o
//   as a subgraph. In other words, in a claw-free graph,
//   any neighbors N(v) in contains no triangles.
//
// Algorithm:
//   We test the triangle-freeness in each N(v) of G~.
//   Here, we can use that |N(v)| <= 2 sqrt(m) due to the 
//   Turan's theorem (hand-shaking type theorem).
//
// Complexity:
//   O(n d^3). Here, d <= 2 sqrt(m).
//   Note that d^3 can be reduced to d^\omega by using
//   fast matrix multiplication.
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

struct graph {
  int n;
  vector<vector<int>> adj;
  graph(int n) : n(n), adj(n) { }
  void add_edge(int src, int dst) {
    adj[src].push_back(dst);
    adj[dst].push_back(src);
  }
};
bool is_claw_free(graph g) {
  int threshold = 0;
  vector<unordered_set<int>> N(g.n);
  for (int s = 0; s < g.n; ++s) {
    threshold += g.adj[s].size();
    for (int v: g.adj[s]) N[s].insert(v);
  }
  threshold = 2 * sqrt(threshold);
  for (int s = 0; s < g.n; ++s) {
    vector<int> &nbh = g.adj[s];
    if (nbh.size() > threshold) return false; // Turan's theorem
    for (int i = 0; i < nbh.size(); ++i) 
      for (int j = i+1; j < nbh.size(); ++j) 
        if (!N[nbh[j]].count(nbh[i])) 
          for (int k = i+2; k < nbh.size(); ++k) 
            if (!N[nbh[k]].count(nbh[i]) && !N[nbh[k]].count(nbh[j]))
              return false;
  }
  return true;
}

int main() {

}

//
// Cograph Recognition
//
// Description:
//   Cographs are recursively defined as follows:
//     (1) A singleton is a cograph
//     (2) A disjoint union of cographs is a cograph
//     (3) A complement of a cograph is a cograph
//   Cographs are also characterized by P_4-free graphs.
//
// Algorithm:
//    By using the above constructive characterization, we have
//      G is a cograph iff all connected components of G~ are cographs.
//    This yields a polynomial time algorithm.
//    Note that there are O(n + m) algorithm for this problem
//    by using lexicographic BFS.
//
// Complexity:
//   O(n^3) in worst case.
// 
// Verify:
//   POJ3236: Michelle's Evaluation
//
#include <iostream> 
#include <vector>
#include <list>
#include <cstdio>
#include <algorithm>
#include <functional>
#include <map>
#include <set>

using namespace std;

#define fst first
#define snd second
#define all(c) ((c).begin()), ((c).end())

struct graph {
  int n;
  vector<vector<int> > adj;
  graph(int n) : n(n), adj(n) { }
  void add_edge(int src, int dst) {
    adj[src].push_back(dst);
    adj[dst].push_back(src);
  }
};

bool is_cograph(graph g, int depth = 0) {
  if (g.n <= 3) return true;
  for (auto &nbh: g.adj) sort(all(nbh));
  vector<int> index(g.n, -1);
  for (int i = 0; i < g.n; ++i) {
    if (index[i] >= 0) continue;
    int size = 0;
    index[i] = size++;
    vector<int> S(1,i), comps;
    while (!S.empty()) {
      int j = S.back(); S.pop_back();
      comps.push_back(j);
      for (int k: g.adj[j]) {
        if (index[k] < 0) {
          index[k] = size++;
          S.push_back(k);
        }
      }
    }
    graph h(size);
    for (int j = 0; j < comps.size(); ++j) 
      for (int k = j+1; k < comps.size(); ++k) 
        if (!binary_search(all(g.adj[comps[j]]), comps[k]))
          h.add_edge(index[comps[j]], index[comps[k]]);
    if (depth > 0 && g.n == h.n) return false; // both G and G' are connected 
    if (!is_cograph(h, depth+1)) return false; // some component is not a cograph
  }
  return true;
}

int main() {
  int n, m;
  scanf("%d %d", &n, &m);

  graph g(n);
  for (int i = 0; i < m; ++i) {
    int u, v;
    scanf("%d %d", &u, &v);
    g.add_edge(u-1, v-1);
  }
  printf("%s\n", (is_cograph_n(g) ? "Yes" : "No"));
}

//
// Minimum Spanning Tree (Kruskal)
//
// 
// Description
//   Finding minimum spanning tree.
//

#include <iostream>
#include <vector>
#include <algorithm>

using namespace std;

#define fst first
#define snd second
#define all(c) ((c).begin()), ((c).end())

struct edge {
  int src, dst;
  int weight;
};
struct graph {
  int n;
  vector<edge> edges;
  graph(int n = 0) : n(n) { }
  void add_edge(int src, int dst, int weight) {
    n = max(n, max(src, dst)+1);
    edges.push_back({src, dst, weight});
  }
  vector<int> p; 
  int root(int i) { 
    return p[i] < 0 ? i : p[i] = root(p[i]); 
  }
  bool unite(int i, int j) {
    if ((i = root(i)) == (j = root(j))) return false;
    if (p[i] > p[j]) swap(i, j);
    p[i] += p[j]; p[j] = i;
    return true;
  }
  int kruskal() {
    p.assign(n, -1);
    sort(all(edges), [](edge x, edge y) {
        return x.weight < y.weight; 
    });
    int result = 0;
    for (auto e: edges) 
      if (unite(e.src, e.dst)) 
        result += e.weight;
    return result;
  }
};

graph random_graph(int n, int d) {
  graph g(n);
  for (int i = 0; i < n; ++i) {
    for (int k = 0; k < d; ++k) {
      int j = rand() % n;
      g.add_edge(i, j, rand() % n);
    }
  }
  return g;
}
int main() {
  auto g = random_graph(100000, 100);
  cout << "[kruskal] " << g.kruskal() << endl;
}

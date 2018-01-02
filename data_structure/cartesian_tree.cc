// 
// Cartesian Tree
//
// Description:
//   For a given sequence xs, the Cartesian tree is recursively
//   defined as follows:
//   - the root is the minimum in xs.
//   - the left child is the Cartesian tree for the the left segment,
//     and the right child is the Cartesian tree for the right segment.
//   It is constructed in O(n) time by the left-to-right traversal.
//
//   By using the Cartesian tree, we can solve
//   - Sorting in O(n log k) time, where k is the number of consecutive
//     sorted subsegments.
//   - LCA(i,j) = RMQ(i,j). Thus, by using the Tarjan's Offline LCA,
//     we can solve m RMQs in O(m a(n)) time.
// 
// Verified:
//   SPOJ RPLN: http://www.spoj.com/problems/RPLN/
// 
// References:
//   C. Levcopoulos and O. Petersson (1989):
//   Heapsort - Adapted for Presorted Files.
//   in Proceedings of the Workshop on Algorithms and Data Structures, 
//   pp. 499-509.

#include <iostream>
#include <vector>
#include <cstdio>
#include <queue>
#include <algorithm>
#include <functional>

using namespace std;

#define fst first
#define snd second
#define all(c) ((c).begin()), ((c).end())

template <class T>
struct CartesianTree {
  int n, root;
  vector<T> xs;
  struct Node {
    int left = -1, right = -1, parent = -1;
  };
  vector<Node> node;
  CartesianTree(const vector<T> &xs) : n(xs.size()), xs(xs), node(n) {
    root = 0;
    for (int i = 1; i < n; ++i) {
      int j = i-1;
      while (node[j].parent >=0 && xs[i] < xs[j]) 
        j = node[j].parent;
      if (xs[i] < xs[j]) {
        node[j].parent = i;
        node[i].left = j;
        root = i;
      } else {
        if (node[j].right >= 0) node[node[j].right].parent = i;
        node[i].left = node[j].right;
        node[i].parent = j;
        node[j].right = i;
      } 
    }
  }
  // In-order traverse gives an original sequence
  void traverse(int t, int tab = 0) {
    if (t < 0) return;
    traverse(node[t].left, tab + 2);
    for (int i = 0; i < tab; ++i) cout << " ";
    cout << xs[t] << endl;
    traverse(node[t].right, tab + 2);
  }
  void traverse() {
    traverse(root);
  }
  // Sorting in O(n log k) by Levcopoulos-Petersson 
  void sort() {
    auto comp = [&](int i, int j) { return xs[i] > xs[j]; };
    priority_queue<int, vector<int>, decltype(comp)> que(comp);
    que.push(root);
    while (!que.empty()) {
      int t = que.top(); que.pop();
      cout << xs[t] << " ";
      if (node[t].left >= 0) que.push(node[t].left);
      if (node[t].right >= 0) que.push(node[t].right);
    }
    cout << endl;
  }

  struct UnionFind {
    vector<int> parent; // parent[root] is the negative of the size.
    UnionFind(int n) : parent(n, -1) { };
    bool unite(int u, int v) { 
      u = root(u); v = root(v);
      if (u == v) return false;
      if (parent[u] > parent[v]) swap(u, v);
      parent[u] += parent[v]; parent[v] = u;
      return true;
    }
    bool find(int u, int v) { return root(u) == root(v); }
    int root(int u) { return parent[u] < 0 ? u : parent[u] = root(parent[u]); }
    int size(int u) { return -parent[root(u)]; }
  };

  struct query { int u, v, a; };
  void rangeMinQueries(vector<query> &queries) {
    vector<vector<query*>> Q(n);
    for (auto &q: queries) {
      Q[q.u].push_back(&q);
      Q[q.v].push_back(&q);
    }
    UnionFind uf(n);
    vector<int> anc(n), color(n);
    iota(all(anc), 0);
    function<void (int)> rec = [&](int u) {
      for (int c: {node[u].left, node[u].right}) {
        if (c < 0) continue;
        rec(c);
        uf.unite(c, u);
        anc[uf.root(u)] = u;
      }
      color[u] = 1;
      for (auto it: Q[u]) {
        if (it->u != u) swap(it->u, it->v);
        if (color[it->v] == 1) it->a = anc[uf.root(it->v)];
      }
    };
    rec(root);
  }
};


int main() {
  int ncase; scanf("%d", &ncase);
  for (int icase = 0; icase < ncase; ++icase) {
    printf("Scenario #%d:\n", icase+1);

    int n, q; scanf("%d %d", &n, &q);
    vector<int> x(n);
    for (int i = 0; i < n; ++i) 
      scanf("%d", &x[i]);
    CartesianTree<int> T(x);

    vector<CartesianTree<int>::query> qs(q);
    for (int i = 0; i < q; ++i) {
      scanf("%d %d", &qs[i].u, &qs[i].v);
      --qs[i].u; --qs[i].v;
    }
    T.rangeMinQueries(qs);
    for (auto q: qs) 
      printf("%d\n", x[q.a]);
  }
}

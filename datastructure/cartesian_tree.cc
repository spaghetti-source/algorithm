// 
// Cartesian Tree
//
// Description:
//   For a given sequence x, the Cartesian tree is recursively
//   defined as follows:
//   - root is the minimum in x.
//   - the left child is the Cartesian tree for the the left of the min.,
//     and the right child is the Cartesian tree for the right of the min.
// 
// Algorithm:
//   Left-to-right construction. 
//
// Applications:
//   In order traversal gives the original sequence. 
//   Sorting can be performed in O(n log k)
//   LCA of two element gives RMQ in original sequence.
//
// Complexity:
//   construction: O(n)
//   sort: O(n log k) by using priority queue.
//   offline RMQ: O(n + q) by Tarjan's offline LCA.
//
// Verified:
//   SPOJ11772, SPOJ1005604 for RMQ.
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
struct cartesian_tree {
  int n, root;
  vector<T> x;
  vector<int> l, r, p;
  cartesian_tree(const vector<T> &x) 
    : n(x.size()), x(x), l(n,-1), r(n,-1), p(n,-1) {
    root = 0;
    for (int i = 1; i < n; ++i) {
      int j = i-1;
      while (p[j] >=0 && x[i] < x[j]) j = p[j];
      if (x[i] < x[j]) {
        p[j] = i;
        l[i] = j;
        root = i;
      } else {
        if (r[j] >= 0) p[r[j]] = i;
        l[i] = r[j];
        p[i] = j;
        r[j] = i;
      } 
    }
  }
  // In-order traverse gives an original sequence
  void traverse(int t, int tab = 0) {
    if (t < 0) return;
    traverse(l[t], tab + 2);
    for (int i = 0; i < tab; ++i) cout << " ";
    cout << x[t] << endl;
    traverse(r[t], tab + 2);
  }
  void traverse() {
    traverse(root);
  }
  // Sorting in O(n log k) by Levcopoulos-Petersson 
  void sort() {
    auto comp = [&](int i, int j) { return x[i] > x[j]; };
    priority_queue<int, vector<int>, decltype(comp)> que(comp);
    que.push(root);
    while (!que.empty()) {
      int t = que.top(); que.pop();
      cout << x[t] << " ";
      if (l[t] >= 0) que.push(l[t]);
      if (r[t] >= 0) que.push(r[t]);
    }
    cout << endl;
  }

  struct union_find {
    vector<int> p; 
    union_find(int n) : p(n, -1) { };
    bool unite(int u, int v) { 
      if ((u = root(u)) == (v = root(v))) return false;
      if (p[u] > p[v]) swap(u, v);
      p[u] += p[v]; p[v] = u;
      return true;
    }
    int root(int u) { return p[u] < 0 ? u : p[u] = root(p[u]); }
  };
  struct query { int u, v, a; };
  void range_min_queries(vector<query> &queries) {
    vector<vector<query*>> Q(n);
    for (auto &q: queries) {
      Q[q.u].push_back(&q);
      Q[q.v].push_back(&q);
    }
    union_find uf(n);
    vector<int> anc(n), color(n);
    iota(all(anc), 0);
    function<void (int)> rec = [&](int u) {
      for (int c: {l[u], r[u]}) {
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
    cartesian_tree<int> T(x);

    vector<cartesian_tree<int>::query> qs(q);
    for (int i = 0; i < q; ++i) {
      scanf("%d %d", &qs[i].u, &qs[i].v);
      --qs[i].u; --qs[i].v;
    }
    T.range_min_queries(qs);
    for (auto q: qs) 
      printf("%d\n", x[q.a]);
  }
}

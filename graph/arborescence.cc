// Minimum Arborescence 
//
// Description:
//   Let G = (V, E) be a weighted directed graph.
//   For a vertex r, an edge-set T is called r-arborerscense if
//     (1) T is a spanning tree (with forgetting directions),
//     (2) for each u in V, indeg_T(u) <= 1, indeg_T(r) = 0.
//   The program finds the minimum weight of r-arborescence.
//
// Algorithm:
//   Gabow et al.'s efficient implementation of Chu-Liu/Edmonds.
//   Use disjoint set and mergeable heap data structure.
//
// Complexity: 
//   O(m log n)
//
// Remark:
//   To recover the edges, it requires additional works. 
//   See Camerini et al.
//
// References: 
//   H. N. Gabow, Z. Galil, T. Spencer, and R. E. Tarjan (1986):
//   Efficient algorithms for finding minimum spanning trees in undirected and directed graphs.
//   Combinatorica, vol 6, pp. 109--122.
//
//   Y. J. Chu and T. H. Liu (1965): 
//   On the shortest arborescence of a directed graph.
//   Science Sinica, vol. 14, pp. 1396--1400.
//
//   J. Edmonds (1967): 
//   Optimum branchings.
//   Journal on Research of the National Bureau of Standards, 71B,
//   pp. 233--240.
//
//   P. M. Camerini, L. Fratta, and F. Maffioli (1979):
//   A note on finding optimum branchings,
//   Networks, vol. 9, pp. 309--312.
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

struct edge { 
  int src, dst; 
  double weight;
};
struct union_find {
  vector<int> p; 
  union_find(int n) : p(n, -1) { };
  bool unite(int u, int v) { 
    if ((u = root(u)) == (v = root(v))) return false;
    if (p[u] > p[v]) swap(u, v);
    p[u] += p[v]; p[v] = u;
    return true;
  }
  bool find(int u, int v) { return root(u) == root(v); }
  int root(int u) { return p[u] < 0 ? u : p[u] = root(p[u]); }
  int size(int u) { return -p[root(u)]; }
};

struct skew_heap {
  struct node {
    node *ch[2];
    edge key;
    double delta;
  } *root;
  skew_heap() : root(0) { }
  void propagate(node *a) {
    a->key.weight += a->delta;
    if (a->ch[0]) a->ch[0]->delta += a->delta;
    if (a->ch[1]) a->ch[1]->delta += a->delta;
    a->delta = 0;
  }
  node *merge(node *a, node *b) {
    if (!a || !b) return a ? a : b;
    propagate(a); propagate(b);
    if (a->key.weight > b->key.weight) swap(a, b);
    a->ch[1] = merge(b, a->ch[1]);
    swap(a->ch[0], a->ch[1]);
    return a;
  }
  void push(edge key) {
    node *n = new node();
    n->ch[0] = n->ch[1] = 0;
    n->key = key; n->delta = 0;
    root = merge(root, n);
  }
  void pop() {
    propagate(root);
    node *temp = root;
    root = merge(root->ch[0], root->ch[1]);
  }
  edge top() {
    propagate(root);
    return root->key;
  }
  bool empty() { 
    return !root;
  }
  void add(double delta) {
    root->delta += delta;
  }
  void merge(skew_heap x) {
    root = merge(root, x.root);
  }
};

const double INF = 9999;
struct minimum_aborescense {
  vector<edge> edges;
  void add_edge(int src, int dst, double weight) {
    edges.push_back({src, dst, weight});
  }
  int n;
  void make_graph() {
    n = 0;
    for (auto e: edges) 
      n = max(n, max(e.src, e.dst)+1);
  }

  double solve(int r) {
    union_find uf(n);
    vector<skew_heap> heap(n);
    for (auto e: edges) 
      heap[e.dst].push(e);

    double score = 0;
    vector<int> seen(n, -1);
    seen[r] = r;
    for (int s = 0; s < n; ++s) {
      vector<int> path;
      for (int u = s; seen[u] < 0;) {
        path.push_back(u);
        seen[u] = s;
        if (heap[u].empty()) return INF; 

        edge min_e = heap[u].top(); 
        score += min_e.weight;
        heap[u].add(-min_e.weight);
        heap[u].pop();

        int v = uf.root(min_e.src);
        if (seen[v] == s) {
          skew_heap new_heap;
          while (1) {
            int w = path.back();
            path.pop_back();
            new_heap.merge(heap[w]);
            if (!uf.unite(v, w)) break;
          }
          heap[uf.root(v)] = new_heap;
          seen[uf.root(v)] = -1;
        }
        u = uf.root(v);
      }
    }
    return score;
  }
};



// ---------------test routine---------------

double bestscore;
vector<edge> bestsolution;

void brute_force(vector<edge> edges, vector<edge> &solution, vector<int> &deg, union_find &uf, int k, double score) {
  if (k == edges.size()) {
    if (solution.size()+1 < uf.p.size()) score = 9999;
    if (bestscore > score) {
      bestscore = score;
      bestsolution = solution;
    }

    return;
  }

  edge e = edges[k];
  if (deg[e.dst] == 0 && uf.root(e.src) != uf.root(e.dst)) {
    union_find uf2 = uf;
    uf2.unite(e.src, e.dst);
    solution.push_back(e);
    deg[e.dst] += 1;
    brute_force(edges, solution, deg, uf2, k+1, score + e.weight);
    deg[e.dst] -= 1;
    solution.pop_back();
  }
  brute_force(edges, solution, deg, uf, k+1, score);
}

double brute_force(vector<edge> edges, int r) {
  int n = 0;
  for (auto e: edges) { n = max(n, max(e.src, e.dst)+1); }

  vector<int> deg(n);
  deg[r] = 1;

  vector<edge> solution;
  union_find uf(n); 
  bestsolution.clear();
  bestscore = 9999;
  brute_force(edges, solution, deg, uf, 0, 0);

  cout << "brute force solution = " << bestscore << endl;
  for (auto e: bestsolution) {
    cout << "(" << e.src << " " << e.dst << ") "; 
  }
  cout << endl;

  return bestscore;
}



#include <ctime>

int main() {

  for (int q = 0; q <= 100; ++q) {
  srand( q );

  minimum_aborescense solver;

  for (int i = 0; i < 100; ++i)
    for (int j = 0; j < i; ++j) {
      if (rand() % 5 == 0) solver.add_edge(i, j, 1 + rand() % 100);
      if (rand() % 5 == 0) solver.add_edge(j, i, 1 + rand() % 100);
    }

  /*
  cout << "input "<< endl;
  for (auto e: solver.edges) {
    cout << e.src << " " << e.dst << " " << e.weight << endl;
  }
  cout << endl;
  */


  solver.make_graph();
  double a = solver.solve(0);
  cout << a << endl;
  /*
  double b = brute_force(solver.edges, 0);
  cout << "seed = " << q << ":  " << a << " " << b <<  endl;
  if (a != b) return 0;
  */
  }
  cout << "*** NO PROBLEM ***" << endl;
}

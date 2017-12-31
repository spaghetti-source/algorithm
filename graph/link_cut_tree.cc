// 
// Link Cut Tree (Slator-Tarjan)
//
// Decription:
//   It maintains rooted arborescences with the following operations
//     link(u,v)     : add link from u to v
//     cut(u)        : cut link from u (to the root direction)
//     lca(u,v)      : least common ancestor of u and v
//
// Algorithm:
//   Classify links into solid and dashed.
//   Here, each vertex has at most one solid link.
//   Then, maintain solid paths by splay trees,
//   which are sorted in the depth of the vertices.
//   see: http://planarity.org/Klein_splay_trees_and_link-cut_trees.pdf
// 
// Complexity:
//   O(log n), amortized.
//
// References:
//   D. D. Sleator and R. E. Tarjan (1983):
//   A Data Structure for Dynamic Trees.
//   Journal oF Computer and System Sciences, vol. 26, no. 3, pp. 362-391.
//
// Verified: AOJ Spaceships


#include <iostream>
#include <cstdio>
#include <vector>
#include <cassert>

using namespace std;

struct link_cut_tree {
  struct node { 
    node *child[2], *parent;
  };
  bool is_root(node *x) {
    return !x->parent || (x->parent->child[0] != x 
                      &&  x->parent->child[1] != x);
  }
  int dir(node *x) { return x->parent && x->parent->child[1] == x; }
  void rot(node* t) {
    node *p = t->parent;
    int d = dir(t);
    p->child[d] = t->child[!d];
    if (p->child[d]) p->child[d]->parent = p;
    if (!is_root(p)) p->parent->child[dir(p)] = t;
    t->parent = p->parent;
    t->child[!d] = p;
    p->parent = t;
  }
  void splay(node *x) {
    while (!is_root(x)) {
      if (!is_root(x->parent)) {
        if (dir(x) == dir(x->parent)) rot(x->parent); 
        else                          rot(x);
      }
      rot(x);
    }
  }
  node *expose(node *x) {
    node *r = 0;
    for (node *p = x; p; p = p->parent) {
      splay(p);
      p->child[1] = r;
      r = p;
    }
    splay(x);
    return r;
  }

  vector<node> ns;
  link_cut_tree(int n) : ns(n) {
    for (int i = 0; i < n; ++i) 
      ns[i].child[0] = ns[i].child[1] = ns[i].parent = 0;
  }
  void link(int x, int y) {
    expose(&ns[x]);
    expose(&ns[y]);
    ns[y].child[1] = &ns[x];
    ns[x].parent = &ns[y];
  }
  void cut(int x) {
    expose(&ns[x]); 
    node *y = ns[x].child[0];
    ns[x].child[0] = y->parent = 0;
  }
  int lca(int x, int y) {
    expose(&ns[x]);
    node *u = expose(&ns[y]);
    return ns[x].parent ? u - &ns[0] : -1;
  }
};

int main() {
  int n, q, t, a, b;
  scanf("%d %d", &n, &q);

  link_cut_tree LCT(n);
  for (int i = 0; i < q; ++i) {
    scanf("%d", &t);
    if (t == 1) {
      scanf("%d %d", &a, &b);
      LCT.link(a-1, b-1);
    } else if (t == 2) {
      scanf("%d", &a);
      LCT.cut(a-1);
    } else {
      scanf("%d %d", &a, &b);
      int c = LCT.lca(a-1, b-1);
      printf("%d\n", c < 0 ? c : c+1);
    }
  }
}

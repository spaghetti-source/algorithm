// 
// Link Cut Tree
//
// Decription:
//   It maintains rooted forests with link/cut operations
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

#include <iostream>
#include <cstdio>
#include <vector>

using namespace std;

struct link_cut_tree {
  struct node { 
    int x, s; // value and sum
    node *ch[2], *p;
  };
  int sum(node *t) { return t ? t->s : 0; }
  node *update(node *t) {
    if (t) t->s = t->x + sum(t->ch[0]) + sum(t->ch[1]);
    return t;
  }
  node *make_node(int x) { return new node({x, x, 0, 0, 0}); }

  int dir(node *t) { return t != t->p->ch[0]; }
  bool is_root(node *t) {
    return !t->p || (t->p->ch[0] != t && t->p->ch[1] != t);
  }
  void connect(node *p, node *t, int d) {
    p->ch[d] = t; if (t) t->p = p;
    update(p);
  }
  void rot(node *t) {
    node *p = t->p;
    int d = dir(t); 
    if (!is_root(p)) connect(p->p, t, dir(p));
    else             t->p = p->p;
    connect(p, t->ch[!d], d);
    connect(t, p, !d);
  }
  void splay(node *t) {
    for (; !is_root(t); rot(t))
      if (!is_root(t->p)) rot(dir(t) == dir(t->p) ? t->p : t);
  }
  node *expose(node *t) {
    node *l = 0;
    for (node *s = t; s; s = s->p) {
      splay(s);
      connect(s, l, 1);
      l = s;
    }
    splay(t);
    return l;
  }
  void link(node *t, node *p) { 
    expose(t);
    expose(p);
    t->p = p;
  }
  void cut(node *t) {
    expose(t); 
    t->ch[0] = t->ch[0]->p = 0;
  }
  node *lca(node *s, node *t) {
    expose(s);
    node *u = expose(t);
    return !s->p ? 0 : u;
  }
  int sum_to_root(node *t) {
    expose(t);
    return sum(t->ch[0]) + t->x;
  }
};

int main() {
  link_cut_tree LCT;
  int n, m; cin >> n >> m;

  vector<link_cut_tree::node*> a(n);
  for (int i = 0; i < n; ++i) {
    a[i] = LCT.make_node(i+1);
  }

  int u, v;
  link_cut_tree::node *p;
  for (int k = 0; k < m; ++k) {
    int t; cin >> t;
    switch (t) {
      case 1:
        cin >> u >> v; --u; --v;
        LCT.link(a[u], a[v]);
        break;
      case 2:
        cin >> u; --u;
        LCT.cut(a[u]);
        break;
      case 3:
        cin >> u >> v; --u; --v;
        p = LCT.lca(a[u], a[v]);
        if (!p) cout << "-1" << endl;
        else    cout << p->x << endl;
      break;
    }
  }
}

//
// Euler Tour Tree
//
// Description:
//   Maintain dynamic unrooted tree with supporting
//   - make_node(x)       : return singleton with value x
//   - link(u,v)          : add link between u and v
//   - cut(uv)            : remove edge uv
//   - sum_in_component(u): return sum of all values in the component
//
// Algorithm:
//   Maintain euler tours by splay trees.
//   This data structure is originally proposed by Miltersen et al,
//   and independently by Fredman and Henzinger.
//
// Complexity:
//   O(log n)
//
// References:
//   M. L. Fredman and M. R. Henzinger (1998):
//   Lower bounds for fully dynamic connectivity problems in graphs.
//   Algorithmica, vol. 22, no. 3, pp. 351–362.
//
//   P. B. Miltersen, S. Subramanian, J. S. Vitter, and R. Tamassia (1994):
//   Complexity models for incremental computation.
//   Theoretical Computer Science, vol. 130. no. 1, pp. 203–236.

\bibitem{old2}

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
#include <unordered_set>

using namespace std;

#define fst first
#define snd second
#define all(c) ((c).begin()), ((c).end())

struct euler_tour_tree {
  struct node {
    int x, s; // value, sum
    node *ch[2], *p, *r;
  };
  int sum(node *t) { return t ? t->s : 0; }
  node *update(node *t) {
    if (t) t->s = t->x + sum(t->ch[0]) + sum(t->ch[1]);
    return t;
  }

  int dir(node *t) { return t != t->p->ch[0]; }
  void connect(node *p, node *t, int d) {
    p->ch[d] = t; if (t) t->p = p;
    update(p);
  }
  node *disconnect(node *t, int d) {
    node *c = t->ch[d]; t->ch[d] = 0; if (c) c->p = 0; 
    update(t);
    return c;
  }
  void rot(node *t) {
    node *p = t->p;
    int d = dir(t); 
    if (p->p) connect(p->p, t, dir(p));
    else      t->p = p->p;
    connect(p, t->ch[!d], d);
    connect(t, p, !d);
  }
  void splay(node *t) {
    for (; t->p; rot(t))
      if (t->p->p) rot(dir(t) == dir(t->p) ? t->p : t);
  }
  void join(node *s, node *t) {
    if (!s || !t) return;
    for (; s->ch[1]; s = s->ch[1]); splay(s);
    for (; t->ch[0]; t = t->ch[0]); connect(t, s, 0);
  }
  node *make_node(int x, node *l = 0, node *r = 0) {
    node *t = new node({x});
    connect(t, l, 0); connect(t, r, 1);
    return t;
  }
  node *link(node *u, node *v, int x = 0) {
    splay(u); splay(v); 
    node *uv = make_node(x, u, disconnect(v,1));
    node *vu = make_node(0, v, disconnect(u,1));
    uv->r = vu; vu->r = uv;
    join(uv, vu);
    return uv;
  }
  void cut(node *uv) {
    splay(uv); disconnect(uv,1); splay(uv->r); 
    join(disconnect(uv,0), disconnect(uv->r,1));
    delete uv, uv->r;
  }
  int sum_in_component(node *u) {
    splay(u);
    return u->s;
  }
};

int main() {
  euler_tour_tree ETT;

  auto a = ETT.make_node(3);
  auto b = ETT.make_node(1);
  auto c = ETT.make_node(4);
  auto d = ETT.make_node(1);
  auto e = ETT.make_node(5);
  auto f = ETT.make_node(9);
  auto g = ETT.make_node(2);

  //        3
  //      1   4
  //     1 5 9 2

  auto ab = ETT.link(a, b);
  ETT.link(a, c);
  ETT.link(b, d);
  ETT.link(b, e);
  auto cf = ETT.link(c, f);
  ETT.link(c, g);

  cout << ETT.sum_in_component(a) << endl;
  cout << ETT.sum_in_component(b) << endl;
  cout << ETT.sum_in_component(c) << endl;
  cout << ETT.sum_in_component(d) << endl;
  cout << endl;

  ETT.cut(ab);
  cout << ETT.sum_in_component(a) << endl;
  cout << ETT.sum_in_component(b) << endl;
  cout << ETT.sum_in_component(c) << endl;
  cout << ETT.sum_in_component(d) << endl;
  cout << endl;

  ETT.cut(cf);
  cout << ETT.sum_in_component(a) << endl;
  cout << ETT.sum_in_component(b) << endl;
  cout << ETT.sum_in_component(c) << endl;
  cout << ETT.sum_in_component(d) << endl;
  cout << endl;

  ETT.link(a,b);
  cout << ETT.sum_in_component(a) << endl;
  cout << ETT.sum_in_component(b) << endl;
  cout << ETT.sum_in_component(c) << endl;
  cout << ETT.sum_in_component(d) << endl;
  cout << endl;
}

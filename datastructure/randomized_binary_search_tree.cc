// 
// Randomized Binary Search Tree (merge-split)
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

using namespace std;

#define fst first
#define snd second
#define all(c) ((c).begin()), ((c).end())

template <class T>
struct randomized_binary_search_tree {
  struct node {
    int s;
    T x;
    node *l, *r;
  } *root;
  randomized_binary_search_tree() : root(0) { }
  node *make_node(const T&x, node *l = 0, node *r = 0) {
    node *t = new node({1, x, l, r});
    if (l) t->s += l->s;
    if (r) t->s += r->s;
    return t;
  }
  node *merge(node *a, node *b) {
    if (!a || !b) return a ? a : b;
    if (rand() % (a->s + b->s) < a->s) {
      a->r = merge(a->r, b);
      return a;
    } else {
      b->l = merge(a, b->l);
      return b;
    }
  }
  pair<node*, node*> split(node *a, const T &x) {
    if (!a) return {0, 0};
    if (a->x < x) {
      auto p = split(a->r, x);
      a->r = p.fst;
      return {a, p.snd};
    } else {
      auto p = split(a->l, x);
      a->l = p.snd;
      return {p.fst, a};
    }
  }
  void insert(const T& x) {
    auto p = split(root, x);
    root = merge(merge(p.fst, make_node(x)), p.snd);
  }
  node *remove_min(node *t) {
    if (!t->l) return t->r;
    t->l = remove_min(t->l);
    return t;
  }
  void remove(const T& x) {
    auto p = split(root, x);
    return merge(remove_min(p.fst), p.snd);
  }
};

int main() {
  randomized_binary_search_tree<int> T;
}

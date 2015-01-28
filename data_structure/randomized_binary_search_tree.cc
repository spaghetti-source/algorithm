// 
// Randomized Binary Search Tree (merge-split)
//

#include <iostream>
#include <vector>
#include <cstdio>
#include <cstdlib>
#include <map>
#include <set>
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
    T x;
    node *l, *r;
    int s;
  } *root;
  randomized_binary_search_tree() : root(0) { }

  node *update(node *t) {
    if (!t) return t;
    t->s = 1;
    if (t->l) t->s += t->l->s;
    if (t->r) t->s += t->r->s;
    return t;
  }
  node *make_node(const T&x, node *l = 0, node *r = 0) {
    return update(new node({x, l, r}));
  }
  node *merge(node *a, node *b) {
    if (!a || !b) return a ? a : b;
    if (rand() % (a->s + b->s) < a->s) {
      a->r = merge(a->r, b);
      return update(a);
    } else {
      b->l = merge(a, b->l);
      return update(b);
    }
  }
  pair<node*, node*> split(node *a, const T &x) {
    if (!a) return {0, 0};
    if (a->x < x) {
      auto p = split(a->r, x);
      a->r = p.fst;
      return {update(a), p.snd};
    } else {
      auto p = split(a->l, x);
      a->l = p.snd;
      return {p.fst, update(a)};
    }
  }
  void insert(const T& x) {
    auto p = split(root, x);
    root = merge(merge(p.fst, make_node(x)), p.snd);
  }
  node *remove(node *t, const T &x) {
    if (!t) return t;
    if (t->x == x) return merge(t->l, t->r);
    if (t->x < x) t->r = remove(t->r, x);
    else          t->l = remove(t->l, x);
    return update(t);
  }
  void remove(const T& x) {
    root = remove(root, x);
  }
  node *find(const T &x) {
    node *t = root;
    while (t && t->x != x) 
      t = (t->x < x ? t->r : t->l);
    return t;
  }
};

#include <ctime>
double tick() {
  static clock_t oldtick;
  clock_t newtick = clock();
  double diff = 1.0*(newtick - oldtick) / CLOCKS_PER_SEC;
  oldtick = newtick;
  return diff;
}

int main() {
  const int n = 10000;
  tick();
  multiset<int> S;
  randomized_binary_search_tree<int> T;
  for (int i = 0; i < n; ++i) {
    int x = rand() % 100;
    T.insert(x);
    S.insert(x);
  }
  for (int i = 0; i < n; ++i) {
    int x = rand() % 100;
    T.remove(x);
    auto it = S.find(x);
    if (it != S.end()) S.erase(it);
  }
  for (int i = 0; i < n; ++i) {
    int x = rand() % 100;
    if ((T.find(x) != 0) != (S.find(x) != S.end())) cout << "!" << endl;
  }
  cout << "tick: " << tick() << endl;
}

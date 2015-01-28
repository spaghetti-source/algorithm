// 
// Splay Tree
//
// Decription:
//   A self-organizing binary search tree.
//
// Algorithm:
//   splay brings a node to a root by a sequence of rotations.
// 
// Complexity:
//   O(log n), amortized.
//

#include <iostream>
#include <cstdio>
#include <vector>
#include <algorithm>

using namespace std;
#define fst first
#define snd second

template <class T>
struct splay_tree {
  struct node {
    T x;
    node *ch[2];
    int size;
  } *root;
  splay_tree() : root(0) { }
  int size(node *t) { return t ? t->size : 0; }
  node *update(node *t) {
    if (t) t->size = size(t->ch[0]) + size(t->ch[1]) + 1;
    return t;
  }
  node *rot(node *t, int d) { // t->ch[d] is returned
    node *s = t->ch[d];
    t->ch[d] = s->ch[!d]; update(t);
    s->ch[!d] = t;        update(s);
    return s;
  }
  node *splay(node *t, int k) { // splay k-th element
    if (!t) return t;
    if (size(t->ch[0]) == k) return t;
    int d = size(t->ch[0]) < k;
    if (d) k = k - size(t->ch[0]) - 1;
    if (size(t->ch[d]->ch[0]) != k) {
      int c = size(t->ch[d]->ch[0]) < k;
      if (c) k = k - size(t->ch[d]->ch[0]) - 1;
      t->ch[d]->ch[c] = splay(t->ch[d]->ch[c], k);
      if (c == d) t = rot(t, d);
      else        t->ch[d] = rot(t->ch[d], !d);
    }
    return rot(t, d);
  }
  T &operator[](int k) {
    root = splay(root, k);
    return root->x;
  }
  void insert(int k, T x) {
    if (size(root) == k) {
      root = update(new node({x, {root, 0}}));
    } else {
      node *t = splay(root, k);
      root = new node({x, {t->ch[0], t}});
      root->ch[1]->ch[0] = 0;
      update(root->ch[1]);
      update(root);
    }
  }
  void erase(int k) {
    root = splay(root, k);
    root->ch[1] = splay(root->ch[1], 0);
    if (root->ch[1]) {
      root->ch[1]->ch[0] = root->ch[0];
      root = update(root->ch[1]);
    } else root = root->ch[0];
  }
};

// for comparison
template <class T>
struct randomized_binary_search_tree {
  struct node {
    T x;
    node *l, *r;
    int size;
  } *root;
  randomized_binary_search_tree() : root(0) { }

  int size(node *t) { return t ? t->size : 0; }
  node *update(node *t) {
    if (t) t->size = 1 + size(t->l) + size(t->r);
    return t;
  }
  node *make_node(const T&x, node *l = 0, node *r = 0) {
    return update(new node({x, l, r}));
  }
  node *merge(node *a, node *b) {
    if (!a || !b) return a ? a : b;
    if (rand() % (size(a) + size(b)) < size(a)) {
      a->r = merge(a->r, b);
      return update(a);
    } else {
      b->l = merge(a, b->l);
      return update(b);
    }
  }
  pair<node*, node*> split(node *t, int k) {
    if (!t) return {0, 0};
    if (size(t->l) < k) {
      auto p = split(t->r, k - size(t->l) - 1);
      t->r = p.fst;
      return {update(t), p.snd};
    } else {
      auto p = split(t->l, k);
      t->l = p.snd;
      return {p.fst, update(t)};
    }
  }
  T &operator[](int k) {
    node *t = root;
    for (node *t = root; size(t->l) != k; ) {
      if (size(t->l) < k) {
        k = k - size(t->l) - 1;
        t = t->l;
      } else t = t->r;
    }
    return t->x;
  }
  void insert(int k, T x) {
    auto p = split(root, k);
    root = merge(merge(p.fst, new node({x, 0, 0, 1})), p.snd);
  }
  node *erase(node *t, int k) {
    if (!t) return t;
    if (size(t->l) == k) return merge(t->l, t->r);
    if (size(t->l) < k) t->r = erase(t->r, k - size(t->l) - 1);
    else                t->l = erase(t->l, k);
    return update(t);
  }
  void erase(int k) {
    root = erase(root, k);
  }
};


// === tick a time ===
#include <ctime>
double tick() {
  static clock_t oldtick;
  clock_t newtick = clock();
  double diff = 1.0*(newtick - oldtick) / CLOCKS_PER_SEC;
  oldtick = newtick;
  return diff;
}

int main() {
  int n = 1000000;
  vector<int> a(n), b(n);
  for (int i = 0; i < n; ++i) {
    a[i] = rand();
    b[i] = rand() % (min(i+1, 100));
  }
  { // splay tree
    tick();
    splay_tree<int> T;
    for (int i = 0; i < n; ++i) {
      T.insert(b[i], a[i]);
    }
    for (int i = n-1; i >= 0; --i) {
      T.erase(b[i]);
    }
    cout << tick() << endl;
  }
  { // randomized bst
    tick();
    randomized_binary_search_tree<int> T;
    for (int i = 0; i < n; ++i) {
      T.insert(b[i], a[i]);
    }
    for (int i = n-1; i >= 0; --i) {
      T.erase(b[i]);
    }
    cout << tick() << endl;
  }
}

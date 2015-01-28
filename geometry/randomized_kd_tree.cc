//
// Randomized KD Tree (for d = 2)
//
// Description
//   Randomized KD tree is a binary space partition tree.
//   Each tree node has a point, direction, and two childs.
//   The left descendants are located to the `left' of the point p
//   and the right descendants are located to the `right' of the point p,
//   where q is `left' of p means q[dir] < p[dir].
//
//   By randomizing the direction, we can easily implement
//   insertion and deletion.
//
// Complexity
//   O(log n) insertion and deletion.
//   O(log n) search for a random points.
//

#include <iostream>
#include <cstdio>
#include <complex>
#include <algorithm>
#include <ctime>

using namespace std;

#define fst first
#define snd second
#define all(c) ((c).begin()), ((c).end())

double tick() {
  static clock_t oldtick;
  clock_t newtick = clock();
  double diff = 1.0*(newtick - oldtick) / CLOCKS_PER_SEC;
  oldtick = newtick;
  return diff;
}

typedef complex<double> point;
struct randomized_kd_tree {
  struct node {
    point p;
    int d, s;
    node *l, *r;
    bool is_left_of(node *x) {
      if (x->d) return real(p) < real(x->p);
      else      return imag(p) < imag(x->p);
    }
  } *root;
  randomized_kd_tree() : root(0) { }
  int size(node *t) { return t ? t->s : 0; }
  node *update(node *t) {
    t->s = 1 + size(t->l) + size(t->r);
    return t;
  }
  pair<node*, node*> split(node *t, node *x) {
    if (!t) return {0, 0};
    if (t->d == x->d) {
      if (t->is_left_of(x)) {
        auto p = split(t->r, x);
        t->r = p.fst; 
        return {update(t), p.snd};
      } else {
        auto p = split(t->l, x);
        t->l = p.snd; 
        return {p.fst, update(t)};
      }
    } else {
      auto l = split(t->l, x);
      auto r = split(t->r, x);
      if (t->is_left_of(x)) {
        t->l = l.fst;
        t->r = r.fst;
        return {update(t), join(l.snd, r.snd, t->d)};
      } else {
        t->l = l.snd;
        t->r = r.snd;
        return {join(l.fst, r.fst, t->d), update(t)};
     }
    }
  }
  node *join(node *l, node *r, int d) {
    if (!l) return r; 
    if (!r) return l;
    if (rand() % (size(l) + size(r)) < size(l)) {
      if (l->d == d) { 
        l->r = join(l->r, r, d);
        return update(l);
      } else {
        auto p = split(r, l);
        l->l = join(l->l, p.fst, d);
        l->r = join(l->r, p.snd, d);
        return update(l);
      }
    } else {
      if (r->d == d) {
        r->l = join(l, r->l, d);
        return update(r);
      } else {
        auto p = split(l, r);
        r->l = join(p.fst, r->l, d);
        r->r = join(p.snd, r->r, d);
        return update(r);
      }
    }
  }
  node *insert(node *t, node *x) {
    if (rand() % (size(t)+1) == 0) {
      auto p = split(t, x);
      x->l = p.fst;
      x->r = p.snd;
      return update(x);
    } else {
      if (x->is_left_of(t)) t->l = insert(t->l, x);
      else                  t->r = insert(t->r, x);
      return update(t);
    }
  }
  void insert(point p) {
    root = insert(root, new node({p, rand()%2}));
  }
  node *remove(node *t, node *x) {
    if (!t) return t;
    if (t->p == x->p) return join(t->l, t->r, t->d);
    if (x->is_left_of(t)) t->l = remove(t->l, x);
    else                  t->r = remove(t->r, x);
    return update(t);
  }
  void remove(point p) {
    node n = {p};
    root = remove(root, &n);
  }
  void closest(node *t, point p, pair<double, node*> &ub) {
    if (!t) return;
    double r = norm(t->p - p);
    if (r < ub.fst) ub = {r, t};
    node *fst = t->r, *snd = t->l;
    double w = t->d ? real(p - t->p) : imag(p - t->p);
    if (w < 0) swap(fst, snd);
    closest(fst, p, ub);
    if (ub.fst > w*w) closest(snd, p, ub);
  }
  point closest(point p) {
    pair<double, node*> ub(1.0/0.0, 0);
    closest(root, p, ub);
    return ub.snd->p;
  }

  // verification
  int height(node *n) {
    return n ? 1 + max(height(n->l), height(n->r)) : 0;
  }
  int height() { 
    return height(root);
  }
  int size_rec(node *n) {
    return n ? 1 + size_rec(n->l) + size_rec(n->r) : 0;
  }
  int size_rec() {
    return size_rec(root);
  }
  void display(node *n, int tab = 0) {
    if (!n) return;
    display(n->l, tab+2);
    for (int i = 0; i < tab; ++i) cout << " ";
    cout << n->p << " (" << n->d << ")" << endl;
    display(n->r, tab+2);
  }
  void display() {
    display(root);
  }
};


int main() {
  srand( 0xdeadbeef );

  int n = 100000;
  vector<point> ps;
  for (int i = 0; i < n; ++i) 
    ps.push_back(point(rand()%n, rand()%n));

  randomized_kd_tree T;

  // sequential insertion
  random_shuffle(all(ps));
  tick();
  for (int i = 0; i < n; ++i) {
    T.insert(ps[i]);
  }
  cout << "insert " << n << " points: "  << tick() << "[s]" << endl;
  cout << "height of " << n << " points: "  << T.height() << endl;

  // search
  random_shuffle(all(ps));
  tick();
  for (int i = 0; i < n; ++i) {
    T.closest(ps[i]);
  }
  cout << "search " << n << " points: "  << tick() << "[s]" << endl;

  // sequential removal
  random_shuffle(all(ps));
  tick();
  for (int i = 0; i < n; ++i) {
    T.remove(ps[i]);
  }
  cout << "remove " << n << " points: "  << tick() << "[s]" << endl;

  // verification
  n = 1000;

  if (T.size_rec() != 0) {
    cout << "ERROR_1" << endl;
    return 0;
  }
  for (int i = 0; i < n; ++i) 
    T.insert(ps[i]);
  for (int i = 1; i < n; i+=2) 
    T.remove(ps[i]);
  if (T.size_rec() != T.root->s) {
    cout << "ERROR_2" << endl;
    return 0;
  }
  vector<point> qs;
  for (int i = 0; i < n; i+=2) 
    qs.push_back(ps[i]);

  for (int i = 0; i < n; ++i) {
    point p(rand(), rand());
    point Tp = T.closest(p);

    point Tq = qs[0];
    for (auto q: qs) 
      if (norm(p - Tq) > norm(p - q)) Tq = q;
    
    if (abs(norm(Tp - p) - norm(Tq - p)) > 1e-8) {
      cout << norm(Tp - p) << endl;
      cout << norm(Tq - p) << endl;
      cout << "ERROR_3" << endl;
      return 0;
    }
  }
  cout << "verification passed" << endl;
}

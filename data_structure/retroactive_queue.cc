//
// Fully Retroactive Queue
//
// Description:
//   It maintains a list of actions ("push" and "pop").
//   We can insert/erase the actions and ask status ("front") in 
//   any position of the list.
//
//   To implement this structure, we maintain the list by a BST 
//   and keep track the number of pushs/pops in each position
//   by using the range addition and range minimum.
//
// Complexity:
//   O(n log n).
//
// References:
//   E. Demaine, J. Iacono, and S. Langerman (2007):
//   "Retroactive data structures". ACM Transactions of Algorithms,
//   vol.3, no.2, pp.1--21.
//

#include <bits/stdc++.h>
using namespace std;

#define fst first
#define snd second
#define all(c) ((c).begin()), ((c).end())

template <class T>
struct RetroactiveQueue {
  struct Node { // Splay Tree
    int type; // 0: none, 1: push, 2: pop
    T value;
    Node *child[2], *parent;
    int a_push, d_push, a_pop, d_pop, min_remain;
  } *head;
  int remain(Node *x) { 
    return !x ? 0 : x->a_push + x->d_push - x->a_pop - x->d_pop;
  }
  RetroactiveQueue() : head(new Node({0})) { }
  Node *update(Node *x) {
    if (!x) return x;
    x->min_remain = remain(x);
    for (int i: {0, 1}) 
      if (x->child[i]) 
        x->min_remain = min(x->min_remain, x->child[i]->min_remain);
    return x; 
  }
  Node *pushdown(Node *x) {
    if (!x) return x;
    for (int i: {0, 1}) {
      if (x->child[i]) {
        x->child[i]->d_push += x->d_push;
        x->child[i]->d_pop += x->d_pop;
      }
    }
    x->a_push += x->d_push;
    x->a_pop += x->d_pop;
    x->d_push = x->d_pop = 0;
    return x;
  }
  int dir(Node *x) { return x->parent && x->parent->child[1] == x; }
  void link(Node *x, Node *y, int d) {
    if (x) x->child[d] = y;
    if (y) y->parent = x;
  }
  void rot(Node *x) {
    int d = dir(x);
    Node *p = x->parent;
    pushdown(p->parent); pushdown(p); pushdown(x);
    link(p->parent, x, dir(p));
    link(p, x->child[!d], d);
    link(x, p, !d);
    update(p); update(x);
  }
  void splay(Node *x) {
    if (!x) return;
    while (x->parent) {
      if (x->parent->parent) {
        if (dir(x) == dir(x->parent)) rot(x->parent); 
        else                          rot(x);
      }
      rot(x);
    }
    pushdown(x);
  }
  Node *insert(Node *x, Node *y) {
    splay(x);
    y->child[0] = x;
    x->parent = y;
    y->child[1] = x->child[1];
    x->child[1] = 0;
    if (y->child[1]) {
      y->child[1]->parent = y;
      if      (y->type == 1) y->child[1]->d_push += 1;
      else if (y->type == 2) y->child[1]->d_pop += 1;
    }
    y->a_push = x->a_push + x->d_push + (y->type & 1);
    y->a_pop  = x->a_pop + x->d_pop + (y->type >> 1);
    update(y->child[0]);
    update(y->child[1]);
    update(y);
    return y;
  }
  Node *insert_push(Node *x, T a) { return insert(x, new Node({1, a})); }
  Node *insert_pop(Node *x) { return insert(x, new Node({2})); }
  Node *erase(Node *x) {
    splay(x);
    Node *y = x->child[1];
    if (!y) {
      x = x->child[0];
      x->parent = 0;
      return x;
    }
    if      (x->type == 1) y->d_push -= 1;
    else if (x->type == 2) y->d_pop -= 1;
    y->parent = 0;
    update(y);
    while (y->child[0]) y = y->child[0];
    splay(y);
    y->child[0] = x->child[0];
    if (y->child[0]) y->child[0]->parent = y;
    update(y->child[0]);
    update(y);
    return y;
  }
  bool valid(Node *x) { splay(x); return x->min_remain >= 0; }
  T front(Node *x) {
    splay(x); 
    int k = x->a_pop + x->d_pop;
    for (Node *y = x; y; ) {
      pushdown(y);
      if (y->a_push > k) {
        x = y;
        y = y->child[0];
      } else {
        y = y->child[1];
      }
    }
    return x->value;
  }


  void display() {
    splay(head);
    display(head, 0);
  }
  void display(Node *x, int tab = 0) {
    if (!x) return;
    display(x->child[0], tab+2);
    cout << string(tab, ' ') << x->type << " " << x->value << " " << x->a_push << " " << x->d_push << " " << x->a_pop << " " << x->d_pop << " " << x->min_remain << endl;
    display(x->child[1], tab+2);
  }
};

int main() {
  RetroactiveQueue<int> que;
  auto i1 = que.insert_push(que.head, 10);
  auto i2 = que.insert_push(i1, 20);
  auto i3 = que.insert_push(i2, 30);
  auto i4 = que.insert_pop(i3);
  //auto i5 = que.insert_pop(i2);
  que.erase(i1);
  que.display();
  cout << que.front(i1) << endl;
  cout << que.front(i2) << endl;
  cout << que.front(i3) << endl;
  cout << que.front(i4) << endl;
}

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
//   Left-to-right method. 
//
// Complexity:
//   O(n).
//
// Applications:
//   Sorting can be performed in O(n log k)
//   LCA of two element gives RMQ in original sequence.
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

using namespace std;

#define fst first
#define snd second
#define all(c) ((c).begin()), ((c).end())

template <class T>
struct cartesian_tree {
  struct node {
    T a;
    node *l, *r, *p;
  } *root;
  vector<node*> ns;
  cartesian_tree(const vector<T> &x) {
    for (int i = 0; i < x.size(); ++i) 
      ns.push_back(new node({x[i]}));
    root = ns[0];
    for (int i = 1; i < x.size(); ++i) {
      node *s = ns[i-1];
      while (s->p && ns[i]->a < s->a) s = s->p;
      if (ns[i]->a < s->a) {
        s->p = ns[i];
        ns[i]->l = s;
        root = ns[i];
      } else {
        if (s->r) s->r->p = ns[i];
        ns[i]->l = s->r;
        ns[i]->p = s;
        s->r = ns[i];
      } 
    }
  }
  // In-order traverse gives an original sequence
  void traverse(node *t, int tab = 0) {
    if (!t) return;
    traverse(t->l, tab + 2);
    for (int i = 0; i < tab; ++i) cout << " ";
    cout << t->a << endl;
    traverse(t->r, tab + 2);
  }
  void traverse() {
    traverse(root);
  }
  // Sorting in O(n log k) by Levcopoulos-Petersson 
  void sort() {
    function<bool (node *, node*)> comp = 
      [&](node *s, node *t) { return s->a > t->a; };
    priority_queue<node*, vector<node*>, decltype(comp)> que(comp);
    que.push(root);
    while (!que.empty()) {
      node *t = que.top(); que.pop();
      cout << t->a << " ";
      if (t->l) que.push(t->l);
      if (t->r) que.push(t->r);
    }
    cout << endl;
  }
};


int main() {
  vector<int> x = {3,1,4,1,5,9,2};
  cartesian_tree<int> T(x);
  T.traverse();
  T.sort();
}

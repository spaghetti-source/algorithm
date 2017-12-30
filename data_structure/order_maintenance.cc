//
// Order Maintenance
//
// - create_node(): return new node x
// - insert(x, y): insert node y after node x
// - erase(x): erase node x from the list
// - order(x, y): return true if x is before y
//
// Running Time:
//   worst case O(1) for create_node, erase, and order.
//   amortized O(log n) for insert; very small constant.
//
// Reference:
//   P. Dietz and D. Sleator (1987): 
//   "Two algorithms for maintaining order in a list".
//   STOC.
//
#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <functional>
#include <vector>
#include <bitset>
#include <cassert>
#include <algorithm>
#include <ctime>
#include <cmath>

#define fst first 
#define snd second
#define all(c) (c).begin(), (c).end()

using namespace std;


// === tick a time ===
#include <ctime>
double tick() {
  static clock_t oldtick;
  clock_t newtick = clock();
  double diff = 1.0*(newtick - oldtick) / CLOCKS_PER_SEC;
  oldtick = newtick;
  return diff;
}

struct order_maintenance {
  using label_type = unsigned long long;
  using signed_label_type = long long;
  const label_type M = ~((~label_type(0))>>1);
  struct node {
    node *prev = 0, *next = 0;
    label_type label; // allows to contain at most sqrt(M) elements
  } *head;

  order_maintenance() {
    head = new node();
    head->next = head->prev = head;
  }
  label_type width(node *x, node *y) {
    label_type ret = y->label - x->label;
    if (ret - 1 >= M) ret += M;
    return ret;
  }
  void insert(node *x, node *u) {
    label_type label = x->label;
    if (width(x, x->next) <= 1) {
      node *mid = x->next, *end = mid->next;
      label_type required = 3;
      while (width(x, end) <= 4 * width(x, mid) && end != x) {
        ++required;
        end = end->next;
        if (end == x) break;
        ++required;
        end = end->next;
        mid = mid->next;
      }
      label_type gap = (x == end ? M : width(x, end)) / required;
      label_type val = end->label;
      while (1) {
        if (end == head) val += M;
        end = end->prev;
        if (end == x) break;
        val -= gap;
        end->label = val;
      }
    }
    u->label = label + width(x, x->next)/2;
    u->next = x->next;
    u->prev = x; 
    u->next->prev = u;
    u->prev->next = u;
  }
  void erase(node *u) { 
    u->prev->next = u->next->prev;
    u->next->prev = u->prev->next;
  }
  node *create_node() {
    return new node;
  }
  bool order(node *x, node *y) {
    return x->label < y->label;
  }
};

int main() {
  order_maintenance T;
  int n = 1000000;
  vector<order_maintenance::node*> nodes(n);
  vector<int> order(n);
  for (int i = 0; i < n; ++i) {
    nodes[i] = T.create_node();
    order[i] = i;
  }
  random_shuffle(all(order));
  tick();
  T.insert(T.head, nodes[order[0]]);
  for (int i = 1; i < n; ++i) {
    //cout << i << endl;
    T.insert(nodes[order[i-1]], nodes[order[i]]);
  }
  cout << tick() << endl;
}

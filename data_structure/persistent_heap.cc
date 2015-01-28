//
// Persistent Heap (based on randomized meldable heap)
//
// Description:
//   Meldable heap with O(1) copy.
//
// Algorithm:
//   It is a persistence version of randomized meldable heap.
//
// Complexity:
//   O(log n) time/space for each operations.
//
// 

#include <iostream>
#include <vector>
#include <cstdio>
#include <queue>
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

template <class T>
struct persistent_heap {
  struct node {
    T x;
    node *l, *r;
  } *root;
  persistent_heap() : root(0) { }
  node *merge(node *a, node *b) {
    if (!a || !b) return a ? a : b;
    if (a->x > b->x) swap(a, b);
    if (rand() % 2) return new node({a->x, merge(a->l, b), a->r});
    else            return new node({a->x, a->l, merge(a->r, b)});
  }
  void merge(persistent_heap b) { root = merge(root, b.root); }
  void push(T x) { root = merge(root, new node({x})); }
  void pop() { root = merge(root->l, root->r); }
  T top() const { return root->x; }
};


int main() {
  priority_queue<int, vector<int>, greater<int>> que;
  persistent_heap<int> heap;
  int n = 10000;
  for (int i = 0; i < n; ++i) {
    int x = rand();
    heap.push(x);
    que.push(x);
  }
  auto tmp_que = que;
  auto tmp_heap = heap;
  while (!que.empty()) {
    if (heap.top() != que.top()) cout << "***" << endl;
    que.pop();
    heap.pop();
  }
  heap = tmp_heap;
  que = tmp_que;

  while (!que.empty()) {
    if (heap.top() != que.top()) cout << "***" << endl;
    que.pop();
    heap.pop();
  }
}

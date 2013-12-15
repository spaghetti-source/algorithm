// Skew Heap
//
// Description:
//   Heap data structure with the following operations.
//
//   1. push a value O(log n)
//   2. pop the smallest value O(log n)
//   3. merge two heaps O(log n + log m)
//   4. add a value to all elements O(1)
//

struct skew_heap {
  struct node {
    node *ch[2];
    int key;
    int delta;
  } *root;
  skew_heap() : root(0) { }
  void propagate(node *a) {
    a->key += a->delta;
    if (a->ch[0]) a->ch[0]->delta += a->delta;
    if (a->ch[1]) a->ch[1]->delta += a->delta;
    a->delta = 0;
  }
  node *merge(node *a, node *b) {
    if (!a || !b) return a ? a : b;
    propagate(a); propagate(b);
    if (a->key > b->key) swap(a, b); // min heap
    a->ch[1] = merge(b, a->ch[1]);
    swap(a->ch[0], a->ch[1]);
    return a;
  }
  void push(int key) {
    node *n = new node();
    n->ch[0] = n->ch[1] = 0;
    n->key = key; n->delta = 0;
    root = merge(root, n);
  }
  void pop() {
    propagate(root);
    node *temp = root;
    root = merge(root->ch[0], root->ch[1]);
  }
  int top() {
    propagate(root);
    return root->key;
  }
  bool empty() { 
    return !root;
  }
  void add(int delta) {
    if (root) root->delta += delta;
  }
  void merge(skew_heap x) { // destroy x
    root = merge(root, x.root);
  }
};

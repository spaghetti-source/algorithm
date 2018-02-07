//
// Leftist Heap
//
// Description:
//
//   Leftist heap is a heap data structure that allows
//   the meld (merge) operation in O(log n) time.
//   Use this for persistent heaps.
//
// Complexity:
//
//   O(1) for top, O(log n) for push/pop/meld
//
// g++ -std=c++17 -O3 -fmax-errors=1 -fsanitize=undefined
#include <bits/stdc++.h>

using namespace std;

#define fst first
#define snd second
#define all(c) ((c).begin()), ((c).end())
#define TEST(s) if (!(s)) { cout << __LINE__ << " " << #s << endl; exit(-1); }

template <class T>
struct LeftistHeap {
  struct Node {
    T key;
    Node *left = 0, *right = 0;
    int dist = 0;
  } *root = 0;
  static Node *merge(Node *x, Node *y) {
    if (!x) return y;
    if (!y) return x;
    if (x->key > y->key) swap(x, y);
    x->right = merge(x->right, y);
    if (!x->left || x->left->dist < x->dist) swap(x->left, x->right);
    x->dist = (x->right ? x->right->dist : 0) + 1;
    return x;
  }
  void push(T key) { root = merge(root, new Node({key})); }
  void pop() { root = merge(root->left, root->right); }
  T top() { return root->key; }
};

//
// Persistent Implementaiton. (allow copy)
//
template <class T>
struct PersistentLeftistHeap {
  struct Node {
    T key;
    Node *left = 0, *right = 0;
    int dist = 0;
  } *root = 0;
  static Node *merge(Node *x, Node *y) {
    if (!x) return y;
    if (!y) return x;
    if (x->key > y->key) swap(x, y);
    x = new Node(*x);
    x->right = merge(x->right, y);
    if (!x->left || x->left->dist < x->dist) swap(x->left, x->right);
    x->dist = (x->right ? x->right->dist : 0) + 1;
    return x;
  }
  void push(T key) { root = merge(root, new Node({key})); }
  void pop() { root = merge(root->left, root->right); }
  T top() { return root->key; }
};

int main() {
  PersistentLeftistHeap<int> heap;
  heap.push(3);
  heap.push(1);
  heap.push(4);
  heap.push(1);
  heap.push(5);
  cout << heap.top() << endl; heap.pop();
  cout << heap.top() << endl; heap.pop();
  auto temp = heap;
  cout << heap.top() << endl; heap.pop();
  cout << heap.top() << endl; heap.pop();
  cout << temp.top() << endl; temp.pop();
  cout << temp.top() << endl; temp.pop();
}

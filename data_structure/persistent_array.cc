// 
// Persistent Array (sqrt decomposition)
//
// Description:
//   An array with O(sqrt(n)) operations (inc. copy)
//
// Algorithm:
//   Store base aray and operation sequence.
//   If the length of operation sequence exceeds sqrt(n),
//   update base array and clear the operation sequence.
//
// Complexity:
//   Copy O(sqrt(n))
//   Get O(sqrt(n))
//   Set O(sqrt(n)) time and space per operation
// 
// Comment:
//   This implementation is much faster than the implementation
//   based on a complete binary tree representation, 
//   which runs in O(log n) time / extra space per operation.
//  
//

#include <iostream>
#include <vector>
#include <cstdio>
#include <functional>
#include <algorithm>

using namespace std;

#define fst first
#define snd second

template <class T>
struct persistent_array {
  const int n;
  T *arr;
  vector<pair<int,T>> op;
  persistent_array(int n, T x = T(0)) : n(n) {
    arr = new T[n];
    fill(arr, arr+n, x);
  }
  const T& get(int k) {
    for (int i = op.size()-1; i >= 0; --i) 
      if (op[i].fst == k) return op[i].snd;
    return arr[k];
  }
  const T& set(int k, const T &x) {
    op.push_back({k, x});
    if (op.size()*op.size() > n) {
      T *new_arr = new T[n];
      copy(arr, arr+n, new_arr);
      arr = new_arr;
      for (int i = 0; i < op.size(); ++i)
        arr[op[i].fst] = op[i].snd;
      op.clear();
    }
    return x;
  }
};


namespace slow {
// theoretically fast, but practically slow due to indirect memory access
template <class T>
struct persistent_array {
  struct node {
    int k;
    T x;
    node *l, *r;
  } *root;
  persistent_array(int n, T x = T(0)) {
    function<node* (int, int)> r = [&](int i, int j) {
      if (i  > j) return (node *)0;
      if (i == j) return new node({i, x});
      int m = (i + j) / 2;
      return new node({m, x, r(i,m-1), r(m+1,j)});
    };
    root = r(0, n-1);
  }
  const T& get(int k) {
    for (node *t = root; ;) {
      if      (t->k < k) t = t->r;
      else if (t->k > k) t = t->l;
      else return t->x;
    }
  }
  const T& set(int k, const T &x) {
    function<node* (node*)> r = [&](node *t) {
      if (!t || t->k == k) return new node({k, x, t->l, t->r});
      if (t->k < k) return new node({t->k, t->x, t->l, r(t->r)});
      else          return new node({t->k, t->x, r(t->l), t->r});
    };
    root = r(root);
    return x;
  }
};
}

#include <ctime>
double tick() {
  static clock_t oldtick;
  clock_t newtick = clock();
  double diff = 1.0*(newtick - oldtick) / CLOCKS_PER_SEC;
  oldtick = newtick;
  return diff;
}

int main() {
  tick();
  vector< persistent_array<int> > arr;

  const int n = 10000; 
  const int m = 100;
  arr.push_back( persistent_array<int>(n));
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < m; ++j) {
      int k = rand() % n;
      arr.back().set(k, rand());
    }
    arr.push_back(arr.back());
  }
  /*
  for (int i = 0; i < arr.size(); ++i) {
    for (int j = 0; j < n; ++j) {
      cout << arr[i].get(j) << " ";
    }
    cout << endl;
  }
  */
  cout << tick() << endl;
}

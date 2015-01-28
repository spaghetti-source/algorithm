// 
// Persistence Union Find
//
// Description:
//   Use persistent array instead of standard array in union find data structure
//
// Complexity:
//   O(a* T(n)), where T(n) is a complexity of persistent array
//

#include <iostream>
#include <vector>
#include <cstdio>
#include <cstdlib>
#include <map>
#include <cmath>
#include <cstring>
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

struct persistent_union_find {
  persistent_array<int> p;
  persistent_union_find(int n) : p(n, -1) { }
  bool unite(int u, int v) { 
    if ((u = root(u)) == (v = root(v))) return false;
    if (p.get(u) > p.get(v)) swap(u, v);
    p.set(u, p.get(u) + p.get(v)); p.set(v, u);
    return true;
  }
  bool find(int u, int v) { return root(u) == root(v); }
  int root(int u) { return p.get(u) < 0 ? u : p.set(u, root(p.get(u))); }
  int size(int u) { return -p.get(root(u)); }
};

int main() {
  persistent_union_find uf(8);

  uf.unite(0,1);
  uf.unite(1,2);
  uf.unite(2,3);

  uf.unite(4,5);
  uf.unite(5,6);
  uf.unite(6,7);

  persistent_union_find tmp = uf;

  cout << uf.find(0,7) << endl;
  uf.unite(3,4);
  cout << uf.find(0,7) << endl;

  cout << tmp.find(0,7) << endl;
}

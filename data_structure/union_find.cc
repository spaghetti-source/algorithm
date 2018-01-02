// 
// Union Find Data Structure
// 
// Description:
//   An union-find data structure (aka. disjoint set data structure) 
//   maintains a disjoint sets and supports the following operations.
//   - unite(u, v): merge sets containing u and v.
//   - find(u, v) : return true if u and v are in the same set
//   - size(u)    : size of the set containing u.
//
//   The weighted version additionally maintains the values for the
//   elements and supports the following operations:
//   - add(u, a)   : value[u] += a
//   - addSet(u, a): value[v] += a for v in set(u)
//   - get(u)      : return value[u]
//   - getSet(u)   : return sum(value[v] for v in set(u))
//
// Complexity:
//   Amortized O(a(n)) for all operations.
//   Here, a(n) is the inverse Ackermann function, which is
//   less than five for a realistic size input.
//
// Verified:
//   AOJ1330 http://judge.u-aizu.ac.jp/onlinejudge/description.jsp?id=1330
//   and other many problems.
//

#include <bits/stdc++.h>
using namespace std;

#define fst first
#define snd second
#define all(c) ((c).begin()), ((c).end())

struct UnionFind {
  vector<int> parent; // parent[root] is the negative of the size.
  UnionFind(int n) : parent(n, -1) { };
  bool unite(int u, int v) { 
    u = root(u); v = root(v);
    if (u == v) return false;
    if (parent[u] > parent[v]) swap(u, v);
    parent[u] += parent[v]; parent[v] += u;
    return true;
  }
  bool find(int u, int v) { return root(u) == root(v); }
  int root(int u) { return parent[u] < 0 ? u : parent[u] = root(parent[u]); }
  int size(int u) { return -parent[root(u)]; }
};

template <class T>
struct WeightedUnionFind {
  struct Data {
    int parent = -1;
    T value = 0, delta = 0, total = 0;
  };
  vector<Data> data;
  WeightedUnionFind(int n) : data(n) { }

  void add(int u, T a) {
    data[u].value += a;
    data[root(u)].total += a;
  }
  void addSet(int u, T a) {
    data[root(u)].delta += a;
    data[root(u)].total -= data[root(u)].parent * a;
  }
  T get(int u) {
    return data[u].value + root_(u).snd;
  }
  T getSet(int u) {
    return data[root(u)].total;
  }
  bool unite(int u, int v) { 
    u = root(u); v = root(v);
    if (u == v) return false;
    if (data[u].parent > data[v].parent) swap(u, v);
    data[u].parent += data[v].parent; 
    data[v].parent = u;
    data[v].delta -= data[u].delta;
    data[u].total += data[v].total;
    return true;
  }
  pair<int, T> root_(int u) {
    if (data[u].parent < 0) return {u, data[u].delta};
    auto p = root_(data[u].parent);
    p.snd += data[u].delta;
    data[u].parent = p.fst;
    data[u].delta = p.snd - data[p.fst].delta;
    return p;
  }
  int root(int u) { return root_(u).fst; } 
  bool find(int u, int v) { 
    return root(u) == root(v); 
  }
  int size(int u) { 
    return -data[root(u)].parent;
  }
};

int test() {
  int n = 5;
  WeightedUnionFind<int> uf(n);
  for (int i = 0; i < n; ++i) {
    uf.add(i, i);
  }
  uf.unite(2, 4);
  for (int i = 0; i < n; ++i) {
    cout << uf.get(i) << " " << uf.getSet(i) << endl;
  }
  cout << endl;
  uf.unite(3, 4);
  for (int i = 0; i < n; ++i) {
    cout << uf.get(i) << " " << uf.getSet(i) << endl;
  }
  cout << endl;
  uf.add(3,10);
  for (int i = 0; i < n; ++i) {
    cout << uf.get(i) << " " << uf.getSet(i) << endl;
  }
  cout << endl;
  uf.addSet(4,5);
  for (int i = 0; i < n; ++i) {
    cout << uf.get(i) << " " << uf.getSet(i) << endl;
  }
}

void AOJ1330() {
  for (int n, m; cin >> n >> m && n; ) {
    WeightedUnionFind<long long> uf(n);
    for (int i = 0; i < m; ++i) {
      char c[2];
      int a, b;
      scanf("%s %d %d", c, &a, &b);
      --a; --b;
      if (c[0] == '!') {
        long long x;
        scanf("%Ld", &x);
        uf.addSet(b, uf.get(a)-uf.get(b)-x);
        uf.unite(a, b);
      } else {
        if (uf.find(a, b)) {
          printf("%Ld\n", uf.get(a)-uf.get(b));
        } else {
          printf("UNKNOWN\n");
        }
      }
    }
  }
}
int main() { AOJ1330(); }

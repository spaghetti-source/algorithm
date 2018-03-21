// 
// Disjoint Sparse Table
//
// Description:
//
//   Let `otimes` be a binary associative operator.
//   The disjoint sparse table is a data structure for a 
//   sequence xs that admits a query
//     prod(i,j) = xs[i] `otimes` ... `otimes` xs[j-1]
//   in time O(1).
//
//   The structure is a segment tree whose node maintains 
//   prod(i,m) and prod(m,j) for all i, j in the segment.
//   Then prod(i,j) is evaluated by finding the node that 
//   splits [i,j) and returning prod(i,m)*prod(m,j).
//
// Complexity:
//
//   preprocessing O(n log n)
//   query O(1)
//
#include <bits/stdc++.h>

using namespace std;

#define fst first
#define snd second
#define all(c) ((c).begin()), ((c).end())
#define TEST(s) if (!(s)) { cout << __LINE__ << " " << #s << endl; exit(-1); }

template <class T, class Op>
struct DisjointSparseTable {
  vector<vector<T>> ys;
  Op otimes;
  DisjointSparseTable(vector<T> xs, Op otimes_) : otimes(otimes_) {
    int n = 1;
    while (n <= xs.size()) n *= 2;
    xs.resize(n);
    ys.push_back(xs);
    for (int h = 1; ; ++h) {
      int range = (2 << h), half = (range /= 2);
      if (range > n) break;
      ys.push_back(xs);
      for (int i = half; i < n; i += range) {
        for (int j = i-2; j >= i-half; --j) 
          ys[h][j] = otimes(ys[h][j], ys[h][j+1]);
        for (int j = i+1; j < min(n, i+half); ++j) 
          ys[h][j] = otimes(ys[h][j-1], ys[h][j]);
      }
    }
  }
  T prod(int i, int j) { // [i, j) query
    --j;
    int h = sizeof(int)*__CHAR_BIT__-1-__builtin_clz(i ^ j);
    return otimes(ys[h][i], ys[h][j]);
  }
};
template <class T, class Op>
auto makeDisjointSparseTable(vector<T> xs, Op op) {
  return DisjointSparseTable<T, Op>(xs, op);
}

int main() {
  vector<int> xs = {3,1,4,1,5,1};
  int n = xs.size();
  auto otimes = [](int a, int b) { return max(a, b); };
  auto dst = makeDisjointSparseTable(xs, otimes);

  for (int i = 0; i < n; ++i) {
    for (int j = i+1; j <= n; ++j) {
      cout << i << " " << j << " " << dst.prod(i, j) << " ";
      int a = xs[i];
      for (int k = i+1; k < j; ++k)
        a = otimes(a, xs[k]);
      cout << a << endl;
    }
  }
}

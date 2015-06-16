//
// Fenwick Tree (aka. binary indexed tree)
//
// Description:
//   A data structure that allows 
//     add(k,a):       b[k] += a
//     sum(k):         b[0] + ... + b[k-1]
//     lower_bound(a): min { k : sum(k) >= a }
//
// Algorithm:
//  [                    1000                     ]
//  [         100         ] [                     ]
//  [   010   ] [         ] [   110   ] [         ]
//  [001] [   ] [011] [   ] [101] [   ] [111] [   ]
//
//  - x[k] maintains the segment b[*,k).
//  - k + (k & (k+1)) is the immediate ancestor of k
//  - k - (k & (k-1)) is the rightmost left segment of k
//
// Complexity:
//   O(log n) access, n space.
//
// Verified:
//   SPOJ3266, SPOJ3267, SPOJ3377

#include <iostream>
#include <vector>
#include <cstdio>
#include <algorithm>

using namespace std;

#define fst first
#define snd second
#define all(c) ((c).begin()), ((c).end())

template <class T>
struct fenwick_tree {
  vector<T> x;
  fenwick_tree(int n) : x(n+1) { } 
  void add(int k, T a) { 
    for (++k; k < x.size(); k += k&-k) x[k] += a;
  }
  T sum(int k) {
    T s = 0;
    for (; k > 0; k &= k-1) s += x[k];
    return s;
  }
  int lower_bound(T a) {
    if (a <= 0) return 0;
    int k = x.size()-1; 
    for (int s: {1,2,4,8,16}) k |= (k >> s); 
    for (int p = ++k; p > 0; p >>= 1, k |= p)
      if (k < x.size() && x[k] < a) a -= x[k]; else k ^= p;
    return k+1;
  }
};

namespace full {

template <class T>
struct fenwick_tree {
  vector<T> x;
  fenwick_tree(int n) : x(n+1) { } 
  // initialize by a constant
  fenwick_tree(int n, T a) : x(n+1, a) { 
    x[0] = 0;
    for (int k = 1; k+(k&-k) <= n; ++k) x[k+(k&-k)] += x[k];
  } 
  // initialize by a vector
  fenwick_tree(vector<T> y) : x(y.size()+1) {
    for (int k = 0; k < y.size(); ++k) x[k+1] = y[k];
    for (int k = 1; k+(k&-k) < x.size(); ++k) x[k+(k&-k)] += x[k];
  }
  // b[k] += a
  void add(int k, T a) { 
    for (++k; k < x.size(); k += k&-k) x[k] += a;
  }
  // sum b[0,k)
  T sum(int k) {
    T s = 0;
    for (; k > 0; k &= k-1) s += x[k];
    return s;
  }
  // min { k : sum(k) >= a }; it requires b[k] >= 0
  int lower_bound(T a) {
    if (a <= 0) return 0;
    int k = x.size()-1; 
    for (int s: {1,2,4,8,16}) k |= (k >> s); 
    for (int p = ++k; p > 0; p >>= 1, k |= p)
      if (k < x.size() && x[k] < a) a -= x[k]; else k ^= p;
    return k+1;
  }
  // max { k : sum(k) <= a }; it requires b[k] >= 0
  int upper_bound(T a) {
    int k = x.size()-1; 
    for (int s: {1,2,4,8,16}) k |= (k >> s); 
    for (int p = ++k; p > 0; p >>= 1, k |= p)
      if (k < x.size() && x[k] <= a) a -= x[k]; else k ^= p;
    return k;
  }
};

} // full

int main() {
  full::fenwick_tree<int> T(6, 1);

  cout << "lower_bound" << endl;
  cout << T.lower_bound(0) << endl;
  cout << T.lower_bound(1) << endl;
  cout << T.lower_bound(2) << endl;
  cout << T.lower_bound(3) << endl;
  cout << T.lower_bound(4) << endl;

  cout << "upper_bound" << endl;
  cout << T.upper_bound(0) << endl;
  cout << T.upper_bound(1) << endl;
  cout << T.upper_bound(2) << endl;
  cout << T.upper_bound(3) << endl;
  cout << T.upper_bound(4) << endl;

  return 0;
  int cases;
  scanf("%d", &cases);
  for (int icases = 0; icases < cases; ++icases) {
    int n, u;
    scanf("%d %d", &n, &u);

    fenwick_tree<int> FT(n+1);
    for (int k = 0; k < u; ++k) {
      int i, j, a;
      scanf("%d %d %d", &i, &j, &a);
      FT.add(i,    a);
      FT.add(j+1, -a);
    }
    int q;
    scanf("%d", &q);
    for (int k = 0; k < q; ++k) {
      int i;
      scanf("%d", &i);
      printf("%d\n", FT.sum(i));
    }
  }
}

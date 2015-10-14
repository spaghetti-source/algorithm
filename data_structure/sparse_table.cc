//
// Sparse Table for Range Minimum Query
//
// Description:
//   The sparse table stores 
//     table[h][i] = min(x[i], ..., x[i+2^h-1])
//   to solve
//     RMQ(i,j) = min { x[i], ..., x[j-1] }.
//  
// Algorithm:
//   table[h+1][i] = min(table[h][i], table[h][i+2^h]).
//   RMQ(i,j) = min(table[h][i], table[h][j-2^h-1].
//
// Complexity:
//   O(n log n) for construction,
//   O(1) for query.
//

#include <iostream>
#include <vector>
#include <cstdio>
#include <numeric>
#include <algorithm>
#include <functional>

using namespace std;

#define fst first
#define snd second
#define all(c) ((c).begin()), ((c).end())

template <class T>
struct sparse_table {
  const vector<T> &x;
  vector<vector<int>> table;
  int argmin(int i, int j) { return x[i] < x[j] ? i : j; }
  sparse_table(const vector<T> &x) : x(x) {
    int logn = sizeof(int)*__CHAR_BIT__-1-__builtin_clz(x.size());
    table.assign(logn+1, vector<int>(x.size()));
    iota(all(table[0]), 0);
    for (int h = 0; h+1 <= logn; ++h) 
      for (int i = 0; i+(1<<h) < x.size(); ++i)
        table[h+1][i] = argmin(table[h][i], table[h][i+(1<<h)]);
  }
  T range_min(int i, int j) { // = min x[i,j)
    int h = sizeof(int)*__CHAR_BIT__-1-__builtin_clz(j-i); // = log2
    return x[argmin(table[h][i], table[h][j-(1<<h)])];
  }
};

int main() {
  vector<int> a = {5,3,8,2,1,5,6};
  int n = a.size();
  sparse_table<int> ST(a);

  for (int i = 0; i < n; ++i) {
    for (int j = i+1; j <= n; ++j) {
      cout << ST.range_min(i, j) << " ";
    }
    cout << endl;
  }
}

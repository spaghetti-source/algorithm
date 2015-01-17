//
// Longest increasing subsequence
//
// Description:
//   We are given a sequence x[1,n]. 
//
// Algorithm:
//   During iterations for k = 1, ..., n, 
//   we maintain two arrays, length and tail:
//     length[k] = length of LIS ending a[k],
//     tail[l]   = last element in LIS with length l,
//   where 
//     tail[length[k]] = k.
//   Since tail is a decreasing sequence, length[k] can be
//   computed in O(log n) time by binary search.
//   (using lower_bound, it finds a strict LIS, and
//    using upper_bound, it finds a weak LIS.)
//
// Complexity:
//   O(n log n).
//
// Comment:
//   The algorithm finds a last element (in dictionary order of index).
//   For other order of the sequence, it is useful to work with the bucket
//     bucket[l] = { k : length[k] == l }.

#include <iostream>
#include <vector>
#include <cstdio>
#include <algorithm>

using namespace std;

#define fst first
#define snd second
#define all(c) ((c).begin()), ((c).end())

template <class T>
vector<T> longest_increasing_subsequence(const vector<T> &x) {
  int n = x.size();
  vector<int> length(n);
  vector<int> tail;
  for (int k = 0; k < n; ++k) {
    length[k] = distance(tail.begin(), upper_bound(all(tail), k, 
      [&](int i, int j) { return x[i] < x[j]; }
    ));
    if (length[k] == tail.size()) tail.push_back(k);
    else                          tail[length[k]] = k;
  }
  int m = *max_element(all(length));
  vector<T> y(m+1);
  for (int i = n-1; i >= 0; --i) 
    if (length[i] == m) y[m--] = x[i];
  return y;
}


int main() {
  vector<int> x = {3,1,4,2};
  auto y = lis(x);
  for (auto a: y) cout << a << " "; cout << endl;
}

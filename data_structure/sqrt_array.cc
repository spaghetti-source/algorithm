//
// SQRT Array
//
// Description:
//   An array with o(n) deletion and insertion
//
// Algorithm:
//   Decompose array into O(sqrt(n)) subarrays.
//   Then the all operation is performed in O(sqrt(n)).
//
// Complexity:
//   O(sqrt(n)); however, due to the cheap constant factor,
//   it is comparable with binary search trees.
//   If only deletion is required, it is better choice.
//
// Verified:
//   erase: SPOJ16016
//

#include <iostream>
#include <vector>
#include <cstdio>
#include <algorithm>

using namespace std;

#define fst first
#define snd second
#define all(c) ((c).begin()), ((c).end())

template <class T>
struct sqrt_array {
  int n;
  vector<vector<T>> x;
  sqrt_array(int n) : n(n) {
    int sqrtn = sqrt(n);
    for (int t = n; t > 0; t -= sqrtn) 
      x.push_back(vector<T>(min(t, sqrtn)));
  }
  void erase(int k) { 
    --n;
    int i = 0;
    for (; k >= x[i].size();  k -= x[i++].size());
    x[i].erase(x[i].begin()+k);
    if (x[i].empty()) x.erase(x.begin()+i);
  }
  void insert(int k, T a = T()) { 
    if (n++ == 0) x.push_back({});
    int i = 0;
    for (; i < x.size() && k >= x[i].size(); k -= x[i++].size());
    if (i == x.size()) x[--i].push_back(a);
    else               x[i].insert(x[i].begin()+k, a);
    int sqrtn = sqrt(n);
    if (x[i].size() > 2*sqrtn) {
      vector<T> y(x[i].begin()+sqrtn, x[i].end());
      x[i].resize(sqrtn);
      x.insert(x.begin()+i+1, y);
    }
  }
  T &operator[](int k) {
    int i = 0;
    for (; k >= x[i].size();  k -= x[i++].size());
    return x[i][k];
  }
  int size() const { return n; }
};

int main() {
  int n = 100;
  sqrt_array<int> x(0);
  cout << "here" << endl;
  for (int i = 0; i < n; ++i) 
    x.insert(rand() % max(x.size(), 1), i);
  for (int i = 0; i < n; ++i) 
    x[i] = i;
  for (int i = 0; i < n; ++i) 
    x.erase(rand() % x.size());
}

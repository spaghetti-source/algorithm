// 
// Horowitz and Sahni's subset sum in O(n 2^{n/2}) time/space
//
// Description:
//   We are given a set of numbers a[0,n) and a number S.
//   The task is to determine there is a subset I such that
//   a[I] == S.
//
// Algorithm:
//   Split a[0,n) into a[0,m) and a[m,n).
//   Let L = { all subset sums in a[0,m) }. 
//   and then check whether L contains (sum - s) 
//   for all subset sums s in a[m,n).
//
// Complexity:
//   O(n 2^{n/2}) = O(n 1.41^n), which is practical for n <= 45.
//   If all values are bounded, it also runs in O(n^2 amax).
//
// Remark:
//   Schroeppel and Shamir improved the space complexity in O(2^{n/4}).
//   It generates L/R in the increasing/decreasing order by using O(2^{n/4}) space:
//   Split a into four component, generate all subset sums, and use priority queue.
//   However, in my implementation (see below), it is much slower than the H-S algorithm.
//
// References:
//   E. Horowitz and S. Sahni (1974):
//   "Computing Partitions with Applications to the Knapsack Problem",
//   Journal of the ACM., vol. 21, no. 2, pp. 277-292.
//
//   R. Schroeppel and A. Shamir (1981):
//   " A T = O(2n=2); S = O(2n=4) algorithm for certain NP-Complete Problems",
//   SIAM Journal of Computing, vol. 10, no. 3, pp. 456-464.
//
#include <iostream>
#include <vector>
#include <cstdio>
#include <algorithm>
#include <functional>
#include <unordered_set>
#include <unordered_map>
#include <queue>

using namespace std;

#define fst first
#define snd second
#define all(c) ((c).begin()), ((c).end())

template <class T>
bool subset_sum(const vector<T>& x, const T &sum) {
  int n = x.size(), m = n / 2;
  vector<T> L(1, 0), R(1, 0);
  for (int i = 0; i < m; ++i) {
    int size = L.size();
    for (int j = 0; j < size; ++j) {
      const T a = L[j] + x[i];
      if (a == sum) return true;
      if (!binary_search(L.begin(), L.begin()+size, a)) 
        L.push_back(a);
    }
    inplace_merge(L.begin(), L.begin()+size, L.end());
  }
  for (int i = m; i < n; ++i) {
    int size = R.size();
    for (int j = 0; j < size; ++j) {
      const T a = R[j] + x[i];
      if (binary_search(L.begin(), L.end(), sum - a)) return true;
      if (!binary_search(R.begin(), R.begin()+size, a))
        R.push_back(a);
    }
    inplace_merge(R.begin(), R.begin()+size, R.end());
  }
  return false;
}

// count the number of solutions
template <class T>
long long subset_sum_count(const vector<T>& x, const T &sum) {
  int n = x.size(), m = n / 2;
  unordered_map<T, long long> L, R;
  L[0] = R[0] = 1;

  for (int i = 0; i < m; ++i) {
    auto X = L;
    for (auto &p: X) 
      L[p.fst + x[i]] += p.snd;
  }
  for (int i = m; i < n; ++i) {
    auto X = R;
    for (auto &p: X) 
      R[p.fst + x[i]] += p.snd;
  }
  long long ans = 0;
  for (auto &p: L) 
    if (R.count(sum - p.fst)) 
      ans += p.snd * R[sum - p.fst];
  return ans;
}

// low memory usage but slow :(
template <class T>
bool subset_sum_lowmemory(const vector<T>& x, const T &sum) {

  vector<T> X = {0}, Y = {0}, Z = {0}, W = {0};
  function<bool(int,int,vector<T>&)> rec = 
    [&](int i, int j, vector<T> &X) {
    for (; i < j; ++i) {
      int n = X.size();
      for (int k = 0; k < n; ++k) {
        T a = x[i] + X[k];
        if (a == sum) return true;
        if (!binary_search(all(X), a)) 
          X.push_back(a);
      }
      inplace_merge(X.begin(), X.begin()+n, X.end());
    }
    return false;
  };
  int n = x.size();
  if (rec(0*n/4, 1*n/4, X) || rec(1*n/4, 2*n/4, Y) || 
      rec(2*n/4, 3*n/4, Z) || rec(3*n/4, 4*n/4, W)) return true;

  typedef tuple<T,int,int> node;
  priority_queue<node, vector<node>, less<node>> P;
  priority_queue<node, vector<node>, greater<node>> Q;
  for (int i = 0; i < X.size(); ++i) 
    P.push(make_tuple(X[i]+Y.back(), i, Y.size()-1));
  for (int i = 0; i < Z.size(); ++i)
    Q.push(make_tuple(Z[i]+W[0], i, 0));

  while (!P.empty() && !Q.empty()) {
    T a = get<0>(P.top()), b = get<0>(Q.top());
    if (a + b == sum) return true;
    if (a + b > sum) { 
      int i = get<1>(P.top()), j = get<2>(P.top()); P.pop();
      if (--j >= 0) P.push(make_tuple(X[i]+Y[j], i, j));
    } else {
      int i = get<1>(Q.top()), j = get<2>(Q.top()); Q.pop();
      if (++j < W.size()) Q.push(make_tuple(Z[i]+W[j], i, j));
    }
  }
  return false;
}



// === tick a time ===
#include <ctime>
double tick() {
  static clock_t oldtick;
  clock_t newtick = clock();
  double diff = 1.0*(newtick - oldtick) / CLOCKS_PER_SEC;
  oldtick = newtick;
  return diff;
}



int main() {

  int n = 50, k = 5, mod = 1000000;
  vector<int> a(n), b(k);
  for (int i = 0; i < n; ++i) 
    a[i] = ((rand() << 16) ^ (rand())) % mod; 
  for (int i = 0; i < k; ++i) 
    b[i] = ((rand() << 16) ^ (rand())) % (n * mod); 

  b[0] = 0;
  for (int i = 0; i < n; ++i) 
    b[0] += a[i];
  if (b[0] % 2) { a[0] += 1; b[0] += 1; }
  b[0] /= 2;


  if (1) {
    tick();
    for (int i = 0; i < k; ++i) 
      cout << subset_sum(a, b[i]) << endl;
    cout << tick()/k << endl;
  }
  cout << "---" << endl;
  {
    tick();
    for (int i = 0; i < k; ++i) 
      cout << subset_sum_count(a, b[i]) << endl;
    cout << tick()/k << endl;
  }
}

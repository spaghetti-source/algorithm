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
// References:
//   E. Horowitz and S. Sahni (1974):
//   "Computing Partitions with Applications to the Knapsack Problem",
//   Journal of the ACM., vol. 21, no. 2, pp. 277-292.
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

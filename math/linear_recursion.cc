// 
// Linear Recurrence Solver
//
// Description:
//   Consider 
//     x[i+n] = a[0] x[i] + a[1] x[i+1] + ... + a[n-1] x[i+n-1].
//   with initial solution x[0], x[1], ..., x[n-1].
//   We compute k-th term of x in O(n^2 log k) time.
//
// Algorithm:
//   Since x[k] is linear in x[0], ..., x[n-1],
//   there exists function f: Z -> R^n such that
//     x[k] = f(k)[0] x[0] + ... + f(k)[n-1] x[n-1].
//   Here, f satisfies the following identities:
//     x[2k] = f(k)[0] x[k] + ... + f(k)[n-1] x[k+n-1]
//           = f(k)[0] (f(k)[0] x[0] + ... + f(k)[n-1] x[n-1])
//           + f(k)[1] (f(k)[0] x[1] + ... + f(k)[n-1] x[n-1+1])
//           + ...
//           = sum{0<=i<n, 0<0j<n} f(k)[i] f(k)[j] x[i+j].
//           = t[0] x[0] + ... + t[2*n-1] x[2*n-1].
//   Also, we have
//     x[n+k] = a[0] x[n] + ... + a[n+k-1] x[n+k-1], 
//   thus
//     t[0] x[0] + ... + t[2*n-1] x[2*n-1]
//     = t[0] x[0] + ... + t[2*n-1] (a[0] x[n] + ... + a[n+k-1] x[n+k-1])
//     = t'[0] x[0] + ... + t[2*n-2]' x[2*n-2].
//     ...
//     = t''[0] x[0] + ... + t''[n-1] x[n-1].
//   This means, we can compute f(2*k) from f(k) in O(n^2) time.      
//   
// Complexity:
//   O(n^2 log k) time, O(n log k) space.
//
#include <iostream>
#include <vector>
#include <cstdio>
#include <algorithm>
#include <functional>

using namespace std;

#define fst first
#define snd second
#define all(c) ((c).begin()), ((c).end())

int linear_recurrence(vector<int> a, vector<int> x, int k) {
  int n = a.size();
  vector<int> t(2*n+1);

  function<vector<int> (int)> rec = [&](int k) {
    vector<int> c(n);
    if (k < n) c[k] = 1;
    else {
      vector<int> b = rec(k / 2);
      fill(all(t), 0);
      for (int i = 0; i < n; ++i) 
        for (int j = 0; j < n; ++j)
          t[i+j+(k&1)] += b[i]*b[j];
      for (int i = 2*n-1; i >= n; --i) 
        for (int j = 0; j < n; j++) 
          t[i-n+j] += a[j]*t[i];
      for (int i = 0; i < n; ++i) 
        c[i] = t[i];
    }
    return c;
  };
  vector<int> c = rec(k);
  int ans = 0;
  for (int i = 0; i < x.size(); ++i) 
    ans += c[i]*x[i];
  return ans;
}


int main() {
  // x[n+k] = x[n] + 2*x[n+1] + 3*x[n+2];
  // x[0] = 6, x[1] = 5, x[2] = 4.
  // 10-th term = 220696
  cout << linear_recurrence({1,2,3}, {6,5,4}, 10) << endl;
}

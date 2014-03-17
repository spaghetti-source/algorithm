//
// Fast linear recursion computation
//
//
// Description
//   Compute m-th term of the sequence
//   x(n+k) = c(0) x(n) + c(1) x(n+1) + ... + c(k-1) x(n+k-1).
//
//
// Algorithm
//   Compute A^m x(0), where A is a companion matrix.
//   Note that we can compute a 2-power of companion matrix in O(k^2 log n) time;
//   let B := A^n then C = B^2 is obtained by
//     C(1) = B(1) * B, 
//     C(2) = C(1) * A,
//     ..., 
//     C(k) = C(k-1) A.
//
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
#include <unordered_map>

using namespace std;

vector<int> operator*(const vector<int> &x, const vector<vector<int>> &A) {
  vector<int> y(A[0].size());
  for (size_t j = 0; j < y.size(); ++j) 
    for (size_t i = 0; i < x.size(); ++i) 
      y[j] += x[i] * A[i][j];
  return y;
}
vector<int> operator*(const vector<vector<int>> &A, const vector<int> &x) {
  vector<int> y(A[0].size());
  for (size_t i = 0; i < y.size(); ++i) 
    for (size_t j = 0; j < x.size(); ++j) 
      y[i] += A[i][j] * x[j];
  return y;
}
// x[m+k] = c[0] x[m] + c[1] x[m+1] + ... + c[k-1] x[m+k-1]
int linear_recurrence(vector<int> c, vector<int> x, int m) {
  size_t k = c.size();
  vector<vector<int>> B(k, vector<int>(k));
  for (size_t i = 0; i+1 < k; ++i) 
    B[i][i+1] = 1;
  B[k-1] = c;
  for (; m > 0; m /= 2) {
    if (m % 2) x = B * x;
    B[0] = B[0] * B;
    for (size_t i = 1; i < k; ++i) {
      for (size_t j = 0; j < k; ++j) { // B[i] = B[i-1] * A
        B[i][j] = B[i-1][k-1] * c[j]; 
        if (j > 0) B[i][j] += B[i-1][j-1];
      }
    }
  }
  return x[0];
}

int main() {
  vector<int> c(3); c[0] = 4; c[1] = 5; c[2] = 6;
  vector<int> x(3); x[0] = 0; x[1] = 1; x[2] = 3;
  for (int m; cin >> m; ) {
    cout << linear_recurrence(c, x, m) << endl;
  }
}

//
// 0/1 Knapsack Problem (Dynamic Programming)
//
// Description:
//   We are given a set of items with profit p_i and weight w_i.
//   The problem is to find a subset of items that maximizes
//   the total profit under the total weight less than some capacity c.
//
//   1) c is small     ==> weight DP
//   2) p is small     ==> price DP
//   3) both are large ==> branch and bound.
//
// Algorithm:
//   weight DP:
//   Let F[a,i] be the max profit for weight <= a using items 1 ... i.
//   Then we have
//     F[a, i] = max(F[a, i-1], F[a-w[i], i-1] + p[i]).
//   The solution is F[c,n].
//
//   Profit DP:
//   Let F[a,i] be the min weight for profit >= a using items 1 ... i.
//   Then we have
//     F[a, i] = min(F[a, i-1], F[a-p[i], i-1] + w[i]).
//   The solution is max { a : F[a,n] <= c }
//
// Complexity:
//   O(n c) for weight DP.
//   O(n (sum p)) for profit DP.
//
// Verified:
//   SPOJ3321.

#include <iostream>
#include <vector>
#include <cstdio>
#include <algorithm>

using namespace std;

// weight DP
// Complexity: O(nc)
//
// F[a] := maximum profit for weight >= a
//
int knapsackW(vector<int> p, vector<int> w, int c) {
  int n = w.size();
  vector<int> F(c+1);
  for (int i = 0; i < n; ++i)
    for (int a = c; a >= w[i]; --a)
      F[a] = max(F[a], F[a-w[i]] + p[i]);
  return F[c];
}

// Profit DP
// Complexity: O(n sum p)
//
// F[a] := minimum weight for profit a
//
int knapsackP(vector<int> p, vector<int> w, int c) {
  int n = p.size(), P = accumulate(all(p), 0);
  vector<int> F(P+1, c+1); F[0] = 0;
  for (int i = 0; i < n; ++i)
    for (int a = P; a >= p[i]; --a)
      F[a] = min(F[a], F[a-p[i]] + w[i]);
  for (int a = P; a >= 0; --a)
    if (F[a] <= c) return a;
}


int main() {
  vector<int> p = {3,1,4,1,5,9};
  vector<int> w = {2,6,5,3,5,8};
  int c = 10;

  cout << knapsackW(p, w, c) << endl;
  cout << knapsackP(p, w, c) << endl;
}

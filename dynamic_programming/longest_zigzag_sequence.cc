//
// Longest ZigZag Subsequence
//
// Description:
// 
//   A sequence xs is zigzag if x[i] < x[i+1], x[i+1] > x[i+2], for all i
//   (initial direction can be arbitrary). The maximum length zigzag 
//   subsequence is computed in O(n) time by a greedy method.
//
//   First, we contract contiguous same numbers. Then, the number of
//   peaks corresponds to the longest zig-zag subsequence.
//
#include <bits/stdc++.h>

using namespace std;

#define fst first
#define snd second
#define all(c) ((c).begin()), ((c).end())
#define TEST(s) if (!(s)) { cout << __LINE__ << " " << #s << endl; exit(-1); }

template <class T>
int longestZigZagSubsequence(vector<T> xs) {
  int n = xs.size(), len = 1, prev = -1; 
  for (int i = 0, j; i < n; i = j) {
    for (j = i+1; j < n && xs[i] == xs[j]; ++j);
    if (j < n) {
      int sign = (xs[i] < xs[j]);
      if (prev != sign) ++len;
      prev = sign;
    }
  }
  return len;
}

// DP for verification
template <class T>
int longestZigZagSubsequenceN(vector<T> A) {
  int n = A.size();
  vector<vector<int>> Z(n, vector<int>(2));
  Z[0][0] = 1;
  Z[0][1] = 1;
  int best = 1;
  for(int i = 1; i < n; i++){
    for(int j = i-1; j>= 0; j--){
      if(A[j] < A[i]) Z[i][0] = max(Z[j][1]+1, Z[i][0]);
      if(A[j] > A[i]) Z[i][1] = max(Z[j][0]+1, Z[i][1]);
    }
    best = max(best, max(Z[i][0],Z[i][1]));
  }
  return best;
}

int main() {
  for (int seed = 0; seed < 10000; ++seed) {
    srand(seed);
    int n = 100;
    vector<int> a(n);
    for (int i = 0; i < n; ++i) {
      a[i] = rand() % n;
    }
    assert(longestZigZagSubsequence(a) == longestZigZagSubsequenceN(a));
  }
}

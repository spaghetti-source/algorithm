//
// Non-succinct Wavelet matrix
//
// Description:
//   Wavelet matrix is a data structure for a list of integers
//   in [0, 2^BIT) that admits the following queries in O(BIT):
//     1) rank(a,k)    = |{ i in [0,k) : x[i] == a }|
//     2) select(a,k)  = min { i in N : rank(a,i) >= k }
//   Usually, wavelet matrix is a succinct data structure;
//   however, the following implementation uses much memory
//   to simplify and accelerate operations.
//
// Algorithm:
//   See Claude, Navarro, and Ordonez.
// 
// Complexity:
//   O(log(BIT) n) for construction,
//   O(log(BIT)) for rank and select.
//
// References:
//   Francisco Claude, Gonzalo Navarro, and Alberto Ordonez (2015):
//   The Wavelet Matrix: An Efficient Wavelet Tree for Large Alphabets.
//   Information Systems, vol. 47, pp. 15--32.
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


struct wavelet_matrix {
  static const int BIT = 16; // 0 <= x[i] < 2^BIT
  int n;
  vector<int> z, l; // number of zeros, and last occurrence
  vector<vector<vector<int>>> r, s; // rank/select vector
  wavelet_matrix(vector<int> x) : n(x.size()), z(n), l(1<<BIT, -1),
    r(BIT, vector<vector<int>>(n+1, vector<int>(2, 0))),
    s(BIT, vector<vector<int>>(n+1, vector<int>(2,-1))) {
    vector<int> y(n);
    for (int h = BIT-1; h >= 0; --h) {
      int j = 0, k = 0, c[2] = {0, 0};
      for (int i = 0; i < n; ++i) {
        int a = (x[i] >> h) & 1;
        if (a) y[j++] = x[i]; else x[k++] = x[i];
        for (int p = 0; p < 2; ++p)
          r[h][i+1][p] = r[h][i][p] + (p == a);
        s[h][++c[a]][a] = i+1;
      }
      z[h] = k;
      for (int i = 0; i < j; ++i) x[k++] = y[i];
    }
    for (int i = n-1; i >= 0; --i) l[x[i]] = i;
  }
  int rank(int a, int k) { // count a in [0,k) 
    if (l[a] < 0) return 0;
    for (int h = BIT-1; h >= 0; --h) {
      int t = (a >> h) & 1;
      k = r[h][k][t];
      if (t) k += z[h];
    }
    return k - l[a];
  }
  int select(int a, int k) { // k-th position+1 of a
    if (k <= 0 || k > rank(a,n)) return -1;
    k += l[a];
    for (int h = 0; h < BIT; ++h) {
      int t = (a >> h) & 1;
      if (t) k -= z[h];
      k = s[h][k][t];
    }
    return k;
  }
};




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
  int n = 1000000;
  vector<int> x(n);
  for (int i = 0; i < n; ++i) {
    x[i] = rand() & ( (1 << wavelet_matrix::BIT) - 1 );
  }
  tick();
  wavelet_matrix W(x);
  printf("construction: %lf[s]\n", tick());

  int count = 0;
  for (int a = 0; a < (1 << W.BIT); ++a) {
    for (int k = 1; k <= W.rank(a,n); ++k) {
      ++count;
      if (W.rank(a, W.select(a,k)) != k) {
        printf("select(%d,%d) = %d\n", a, k, W.select(a, k));
        printf("rank(%d,%d) = %d\n", a, W.select(a, k), W.rank(a, W.select(a, k)));
        return -1;
      }
    }
  }
  printf("%lf [s/op]\n", tick() / count);
  return 0;

  /*
  //                0 1  2 3 4 5 6  7  8 9  0  1  2 3 4  5 6  7 8 9 0  1 2 3  4  5 6 7 8 9  0 1
  vector<int> x = {11,0,15,6,5,2,7,12,11,0,12,12,13,4,6,13,1,11,6,1,7,10,2,7,14,11,1,7,5,4,14,6};
  printf("%d\n", W.select(0,2));
  printf("%d\n", W.rank(13, 12));
  printf("%d\n", W.rank(13, 13));
  printf("%d\n", W.rank(13, 14));
  printf("%d\n", W.rank(13, 15));
  printf("%d\n", W.rank(13, 16));
  */
}

//
// Bradley-Terry model for pairwise comparison
//
// Description:
//   Consider pairwise comparisons between n players.
//   This model assumes that each player i has a strength w_i,
//   and player i beats player j with probability w_i/(w_i + w_j).
//   The algorithm estimates the strengths from a comparison data.
//
// Algorithm:
//   We maximize the log-likelihood:
//     L(w) := sum ( log(w_i) - log(w_i + w_j) ).
//   We derive an iterative algorithm. Let w be a current solution.
//   We observe that the second term is an obstraction for maximizatino
//   because log(w_i + w_j) is concave. we Wpproximate this term by
//     log(w_i' + w_j') <= log(w_i + w_j) + (w_i' + w_j')/(w_i + w_j).
//   Then, we have a minorization
//     L(w') >= sum ( log(w_i') - (w_i' + w_j')/(w_i + w_j) ) + const.
//   Here, the right hand side has a unique maximum that is given by
//     w_i = n_i ( sum 1/(w_i' + w_j') )^{-1}.
//   By iterating this process, w converges to a local optimal solution.
//   Under mild assumptions, the likelihood has a unique optimal solution;
//   therefore this algorithm converges to a global optimal solution.
//
// Complexity:
//   O(n^2) per iteration. The number of iterations are usually small.
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

struct bradley_terry {
  int n;
  vector<double> w;
  vector<vector<int>> a;
  bradley_terry(int n) : n(n), w(n,1) { regularize(); }

  // reguralization avoids no-match pairs
  void regularize() {
    a.assign(n, vector<int>(n, 1));
    for (int i = 0; i < n; ++i)
      a[i][i] = n-1;
  }

  // win beats lose num times
  void add_match(int win, int lose, int num = 1) {
    a[win][lose] += num;
    a[win][win]  += num;
  }

  // estimate the strengths
  void learning() {
    for (int iter = 0; iter < 100; ++iter) {
      double norm = 0;
      vector<double> z(n);
      for (int i = 0; i < n; ++i) {
        double sum = 0;
        for (int j = 0; j < n; ++j) 
          if (i != j) sum += (a[i][j] + a[j][i]) / (w[i] + w[j]);
        z[i] = a[i][i] / sum;
        norm += z[i];
      }
      double err = 0;
      for (int i = 0; i < n; ++i)  {
        err += abs(w[i] - z[i] / norm);
        w[i] = z[i] / norm;
      }
      if (err < 1e-6) break;
    }
  }

};


// 2014, NPB (Nippon Professional Baseball) Central League
int data[6][6] = {
  { 0, 13, 13, 16, 11, 13}, // Giants
  {11,  0, 14, 12, 16, 13}, // Tigers
  {10, 10,  0, 14, 15, 16}, // Carp
  { 8, 11, 10,  0, 14, 11}, // Dragons
  {13,  8,  8,  9,  0, 16}, // BayStars
  {11, 11,  8, 12,  8,  0}, // Swallows
};

int main() {
  int n = 6;
  bradley_terry BT(n);
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      BT.add_match(i, j, data[i][j]);
    }
  }
  BT.learning();

  for (int i = 0; i < n; ++i) {
    cout << BT.w[i] << endl;
  }
}

// 
// Minimum assignment (simplified Jonker-Volgenant)
//
// Description:
//   We are given a cost table of size n times m with n <= m.
//   It finds a minimum cost assignment, i.e.,
//     min sum_{ij} c(i,j) x(i,j)
//     where sum_{i in [n]} x(i,j)  = 1,
//           sum_{j in [m]} x(i,j) <= 1.
//
// Algorithm:
//   Simplified version of Jonker-Volgenant algorithm,
//   which omits a heuristical initialization step.
//
// Complexity:
//   O(n^3).
//   Much faster than the Kuhn-Munkres algorithm.
//
// Note:
//   It finds minimum cost maximal matching.
//   To find the minimum cost non-maximal matching,
//   we add n dummy vertices to the right side.
//
// Verified: 
//   SPOJ 286: Selfish City
//
// References:
//   R. Jonker and A. Volgenant (1987):
//   A shortest augmenting path algorithm for dense and sparse linear assignment problems.
//   Computing, vol.38, no.4, pp.325-340.
//
//   A. Volgenant (1996): 
//   Linear and Semi Assignment Problems: a core oriented approach.
//   Computers and Operations Research, vol.23, pp.917-932.
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
#include <unordered_set>

using namespace std;

#define fst first
#define snd second
#define all(c) ((c).begin()), ((c).end())

typedef int value_type;
const value_type inf = 99999999;

value_type min_assignment(const vector<vector<value_type>> &c) {
  const int n = c.size(), m = c[0].size(); // assert(n <= m);
  vector<value_type> v(m), dist(m);        // v: potential
  vector<int> matchL(n,-1), matchR(m,-1);  // matching pairs
  vector<int> index(m), prev(m);
  iota(all(index), 0);

  auto residue = [&](int i, int j) { return c[i][j] - v[j]; };
  for (int f = 0; f < n; ++f) {
    for (int j = 0; j < m; ++j) {
      dist[j] = residue(f, j);
      prev[j] = f;
    }
    value_type w;
    int j, l;
    for (int s = 0, t = 0;;) {
      if (s == t) {
        l = s; w = dist[index[t++]]; 
        for (int k = t; k < m; ++k) {
          j = index[k];
          value_type h = dist[j];
          if (h <= w) {
            if (h < w) { t = s; w = h; }
            index[k] = index[t]; index[t++] = j;
          }
        }
        for (int k = s; k < t; ++k) {
          j = index[k];
          if (matchR[j] < 0) goto aug;
        }
      }
      int q = index[s++], i = matchR[q];
      for (int k = t; k < m; ++k) {
        j = index[k];
        value_type h = residue(i,j) - residue(i,q) + w;
        if (h < dist[j]) { 
          dist[j] = h; prev[j] = i;
          if (h == w) {
            if (matchR[j] < 0) goto aug;
            index[k] = index[t]; index[t++] = j;
          }
        }
      }
    }
aug:for(int k = 0; k < l; ++k) 
      v[index[k]] += dist[index[k]] - w;
    int i;
    do {
      matchR[j] = i = prev[j]; 
      swap(j, matchL[i]);
    } while (i != f);
  }
  value_type opt = 0;
  for (int i = 0; i < n; ++i) 
    opt += c[i][matchL[i]]; // (i, matchL[i]) is a solution
  return opt;
}

int main() {
  int ncase;
  scanf("%d", &ncase);
  for (int icase = 0; icase < ncase; ++icase) {
    int n, m; scanf("%d %d", &n, &m);
    vector<vector<int>> c(n, vector<int>(m+n, 0));
    for (int s, t, u; scanf("%d %d %d", &s, &t, &u) && s; ) {
      c[s-1][t-1] = -u;
    }
    printf("%d\n", -min_assignment(c));
  }
}

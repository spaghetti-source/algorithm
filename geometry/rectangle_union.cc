//
// Area of Union of Rectangles (Bentley)
//
// Description:
//   For a given set of rectangles, it gives the area of the union.
//   This problem is sometines called the Klee's measure problem [Klee'77].
//
// Algorithm:
//   Bentley's plane-sweep algorithm [Bentley'77].
//   We first apply the coordinate compression technique.
//   Then the y-structure, which is called measure tree, is simply implemented 
//   by using segment tree data structure.
//
// Complexity:
//   O(n log n) time and O(n) space.
//
// Verify:
//   LightOJ 1120: Rectangle Union
//
// References:
//
//   V. Klee (1977): 
//   Can the measure of \cup[a_i, b_i] be computed in less than O(n \log n) steps?
//   American Mathematical Monthly, vol.84, pp. 284--285.
//
//   J. L. Bentley (1977): 
//   Algorithms for Klee's rectangle problems.
//   Unpublished notes, Computer Science Department, Carnegie Mellon University.
//

#include <iostream>
#include <vector>
#include <map>
#include <cstdio>
#include <algorithm>
#include <functional>

using namespace std;

#define fst first
#define snd second
#define all(c) ((c).begin()), ((c).end())

struct rectangle { int xl, yl, xh, yh; };
long long rectangle_area(vector<rectangle> rs) {
  vector<int> ys; // coordinate compression
  for (int i = 0; i < rs.size(); ++i) {
    ys.push_back(rs[i].yl);
    ys.push_back(rs[i].yh);
  }
  sort(all(ys)); ys.erase(unique(all(ys)), ys.end());

  int n = ys.size(); // measure tree
  vector<int> C(8*n), A(8*n);
  function<void (int,int,int,int,int,int)> aux = 
  [&](int a, int b, int c, int l, int r, int k) {
    if ((a = max(a,l)) >= (b = min(b,r))) return;
    if (a == l && b == r) C[k] += c;
    else {
      aux(a, b, c, l, (l+r)/2, 2*k+1);
      aux(a, b, c, (l+r)/2, r, 2*k+2);
    }
    if (C[k]) A[k] = ys[r] - ys[l];
    else      A[k] = A[2*k+1] + A[2*k+2];
  };

  struct event { int x, l, h, c; }; // plane sweep
  vector<event> es;
  for (auto r: rs) {
    int l = distance(ys.begin(), lower_bound(all(ys), r.yl));
    int h = distance(ys.begin(), lower_bound(all(ys), r.yh));
    es.push_back({r.xl, l, h, +1});
    es.push_back({r.xh, l, h, -1});
  }
  sort(all(es), [](event a, event b) { return a.x != b.x ? a.x < b.x : a.c > b.c; });
  long long area = 0, prev = 0;
  for (auto &e: es) {
    area += (e.x - prev) * A[0];
    prev = e.x;
    aux(e.l,e.h,e.c,0,n,0);
  }
  return area;
}


int main() {
  int ncase; scanf("%d", &ncase);
  for (int icase = 0; icase < ncase; ++icase) {
    int n; scanf("%d", &n);
    vector<rectangle> rs(n);
    for (int i = 0; i < n; ++i) 
      scanf("%d %d %d %d", &rs[i].xl, &rs[i].yl, &rs[i].xh, &rs[i].yh); 
    printf("Case %d: %lld\n", icase+1, rectangle_area(rs));
  }
}

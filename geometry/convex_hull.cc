//
// Convex hull of 2D points
//
// Description:
//   Find a convex hull of point sets.
//
// Algorithm:
//   Andrew's monotone chain.
//
// References:
//   A. M. Andrew (1979): 
//     Another efficient algorithm for convex hulls in two dimensions.
//     Information Processing Letters, vol.9, pp.216-219.
//

#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <map>
#include <vector>
#include <cstring>
#include <functional>
#include <algorithm>
#include <cmath>
#include <complex>

using namespace std;

#define ALL(c) c.begin(), c.end()
#define FOR(i,c) for(typeof(c.begin())i=c.begin();i!=c.end();++i)
#define REP(i,n) for(int i=0;i<n;++i)
#define fst first
#define snd second

typedef long long Value;
typedef complex<Value> Point;
#define X real()
#define Y imag()
Value dot(Point a, Point b)   { return real(conj(a)*b); }
Value cross(Point a, Point b) { return imag(conj(a)*b); }
Value dist2(Point a, Point b) { return dot(a-b, a-b); }

int ccw(Point a, Point b, Point c) {
  b -= a; c -= a;
  if (cross(b,c) > 0)      return +1; // counter clockwise
  if (cross(b,c) < 0)      return -1; // clockwise
  if (dot(b,c) < 0)        return +2; // c--a--b on line
  if (dot(b,b) < dot(c,c)) return -2; // a--b--c on lne
  return 0;
}

// Convex Hull
//
// Algorithm:
//   Andrew's monotone chain
namespace std { 
  bool operator < (Point a, Point b) { // bottom-left
    return a.Y != b.Y ? a.Y < b.Y : a.X < b.X; 
  }
}
vector<Point> convexHull(vector<Point> p) {
  int n = p.size(), k = 0;
  vector<Point> h(2*n);
  sort(ALL(p));
  for (int i = 0; i < n; h[k++] = p[i++]) 
    while (k >= 2 && ccw(h[k-2], h[k-1], p[i]) <= 0) --k;
  for (int i = n-2, t = k+1; i >= 0; h[k++] = p[i--]) 
    while (k >= t && ccw(h[k-2], h[k-1], p[i]) <= 0) --k;
  return vector<Point>(h.begin(), h.begin() + k - (k > 1));
}


// SPOJ 26: Build the Fence
#define prev(p, i) ((i)-1>=0       ? p[(i)-1]: p[(i)-1+p.size()])
#define curr(p, i) ((i) < p.size() ? p[i]    : p[(i) - p.size()])
#define next(p, i) ((i)+1<p.size() ? p[(i)+1]: p[(i)+1-p.size()])
map<Point, int> dic;
void solve() {
  static int _count;
  if (_count++ > 0) printf("\n");

  dic.clear();
  int n; scanf("%d", &n);
  vector<Point> p(n);
  REP(i,n) {
    scanf("%lld %lld", &p[i].X, &p[i].Y);
    if (dic[p[i]] == 0) dic[p[i]] = i+1;
  }
  vector<Point> ch = convexHull(p);
  double len = 0;
  vector<int> out;
  REP(i, ch.size()) len += sqrt(dist2(curr(ch,i), next(ch,i)));
  printf("%.2lf\n", len);
  REP(i, ch.size()) {
    if (i > 0) printf(" ");
    printf("%d", dic[ch[i]]);
  }
  printf("\n");
}

int main() {
  int T; scanf("%d", &T);
  while (T--) {
    solve();
  }
}

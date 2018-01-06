// 
// 3D Coordinate-Wise Domination
//
// Description:
//
//   Point (x,y,z) dominates (x',y',z') if
//     x < x', y < y', and z < z'
//   holds. Kung-Luccio-Preparata proposed an algorithm to compute
//   the all set of dominating points in O(n log n) time.
//
//   It maintains a data structure to check the domination in (y,z) plane,
//   and proceed the points in the decreasing order of x.
//   By the processing order, the new point is never dominated by the latter
//   points and is dominated by the previous points if its (y,z) coordinate
//   is dominated by them.
//
// Complexity:
// 
//   O(n log n). By using this method recursively,
//   we can solve d-dimensional domination in O(n log^{d-2} n).
//
// Reference:
//   
//   Hsiang-Tsung Kung, Fabrizio Luccio, Franco P. Preparata (1975):
//   "On finding the maxima of a set of vectors." Journal of the ACM,
//   vol.22, no.4, pp.469-476.
//

#include <bits/stdc++.h>
using namespace std;

#define fst first
#define snd second
#define all(c) ((c).begin()), ((c).end())

using Real = int;
struct Point {
  Real x, y, z;
  bool operator < (Point p) const {
    if (x != p.x) return x < p.x;
    if (y != p.y) return y < p.y;
    return z < p.z;
  }
};
int domination(vector<Point> xs) {
  multiset<Point> frontier;
  auto bad = [&](multiset<Point>::iterator it) {
    auto jt = next(it);
    if (jt == frontier.end()) return false;
    return it->y < jt->y && it->z < jt->z; 
  };
  int n = xs.size();
  sort(all(xs));
  int count = 0;
  for (int i = n-1; i >= 0; --i) {
    auto proc = [&] {
      if (i < n-1 && xs[i].x == xs[i+1].x) return true;
      Point p(xs[i]); p.x = 0;
      auto it = frontier.insert(p);
      if (bad(it)) {
        frontier.erase(it);
        return false;
      } else {
        while (it != frontier.begin() && bad(prev(it))) 
          frontier.erase(prev(it));
        return true;
      }
    };
    if (proc()) ++count;
  }
  return count;
}


int main() {
  int ncase; scanf("%d", &ncase);
  for (int icase = 0; icase < ncase; ++icase) {
    int n; scanf("%d", &n);
    vector<Point> ps(n);
    for (int i = 0; i < n; ++i) {
      scanf("%d %d %d", &ps[i].x, &ps[i].y, &ps[i].z);
      ps[i].x = -ps[i].x;
      ps[i].y = -ps[i].y;
      ps[i].z = -ps[i].z;
    }
    printf("%d\n", domination(ps));
  }
}

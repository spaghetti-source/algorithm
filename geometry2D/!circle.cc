//
// Basic Predicates for circles
//
#include <iostream>
#include <vector>
#include <cstdio>
#include <cstdlib>
#include <map>
#include <cmath>
#include <cstring>
#include <complex>
#include <functional>
#include <algorithm>
#include <unordered_map>
#include <unordered_set>

using namespace std;

#define fst first
#define snd second
#define all(c) ((c).begin()), ((c).end())

typedef complex<double> point;
namespace std {
bool operator < (point p, point q) {
  if (real(p) != real(q)) return real(p) < real(q);
  return imag(p) < imag(q);
}
};
double dot(point p, point q)   { return real(conj(p) * q); }
double cross(point p, point q) { return imag(conj(p) * q); }
double EPS = 1e-8;
int sign(double x) {
  if (x < -EPS) return -1;
  if (x > +EPS) return +1;
  return 0;
}

struct circle { point p; double r; };
struct line { point p, q; };

vector<point> intersect(circle C, circle D) {
  double d = abs(C.p - D.p);
  if (sign(d - C.r + D.r) > 0) return {};       // too far
  if (sign(d - abs(C.r - D.r)) <= 0) return {}; // too close
  double a = (C.r*C.r - D.r*D.r + d*d)/(2*d);
  double h = sqrt(C.r*C.r - a*a);
  point v = (C.p - D.p) / d;
  if (sign(h) == 0) return {C.p + v*a};         // touch
  return {C.p + v*a + point(0,1)*v*h,           // intersect
          C.p + v*a - point(0,1)*v*h};
}
//
// Intersection of line L and circle C.
//
// Let p(t) = L.p + t (L.q - L.p). We solve
//   |p(t) - C.p|^2 = C.r^2
// By letting u = L.p - L.q, v = L.p - C.p, the above is
//   t^2 (u,u) + 2 t (u,v) + (v,v) = C.r^2
//       ~~a~~       ~~b~~   ~~c~~
// Thus
//   det = b^2 - ac,
//   t in { (b + sqrt(det))/a, c/(b + sqrt(det)) }
//
vector<point> intersect(line L, circle C) {
  point  u = L.p - L.q, v = L.p - C.p;
  double a = norm(u), b = dot(u,v), c = norm(v),
         det = b*b - a*c, r = b + sqrt(max(det, 0.0));
  if (sign(det) <  0) return {};              // no solution
  if (sign(det) == 0) return {L.p - b/a * u}; // touch
  return {L.p - (b + sqrt(det))/a*u,
          L.p - c/(b + sqrt(det))*u};
}

//
// Tangent point(s) of point p and circle C
// 
// Let q be a tangent point.
// The angle between q-p-c.p is 
//   sin(t) = r/|p - c.p|.
// and the solution is
//   p + (c.p - p) * exp(\pm it).
//
// Verified: SPOJ18531
//
vector<point> tangent(point p, circle c) {
  double sin2 = c.r*c.r/norm(p - c.p);
  if (sign(1-sin2)  < 0) return {};
  if (sign(1-sin2) == 0) return {p};
  point z(sqrt(1-sin2), sqrt(sin2));
  return {p+(c.p-p)*conj(z), p+(c.p-p)*z};
}

int main() {
}

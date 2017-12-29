//
//
// http://cs.smith.edu/classwiki/images/c/c8/Mitchell.MinPerim.pdf
// https://www.ics.uci.edu/~eppstein/pubs/EppOveRot-DCG-92.pdf
// https://en.wikipedia.org/wiki/Well-separated_pair_decomposition
// http://www.algorithmic-solutions.info/leda_manual/geo_alg.html
//
// future work:
//   boolean operations (polygonal region overlay; 人間には実装無理では)
//   ham sandwitch cut https://www.ti.inf.ethz.ch/ew/lehre/CG09/materials/v11.pdf

#include <list>
#include <iostream>
#include <vector>
#include <fstream>
#include <cstdio>
#include <algorithm>
#include <cmath>
#include <map>
#include <queue>
#include <set>
#include <functional>
#include <ctime>
#include <numeric>
#include <unordered_set>
#include <unordered_map>
#include <random>
#include <list>

using namespace std;

#define fst first
#define snd second
#define all(c) ((c).begin()), ((c).end())
#define TEST(x,a) { auto y=(x); if (sign(y-a)) { cout << "line " << __LINE__  << #x << " = " << y << " != " << a << endl; exit(-1); } }
double urand() { return rand() / (1.0 + RAND_MAX); }

// === tick a time ===
#include <ctime>
double tick() {
  static clock_t oldtick;
  clock_t newtick = clock();
  double diff = 1.0*(newtick - oldtick) / CLOCKS_PER_SEC;
  oldtick = newtick;
  return diff;
}

template <class T>
ostream &operator<<(ostream &os, const vector<T> &v) {
  os << "[";
  for (int i = 0; i < v.size(); os << v[i++])
    if (i > 0) os << " ";
  os << "]";
  return os;
}
template <class T>
ostream &operator<<(ostream &os, const vector<vector<T>> &v) {
  os << "[";
  for (int i = 0; i < v.size(); os << v[i++])
    if (i > 0) os << endl << " ";
  os << "]";
  return os;
}

const double PI = acos(-1.0);

// implementation note: use EPS only in this function
// usage note: check sign(x) < 0, sign(x) > 0, or sign(x) == 0
const double EPS = 1e-8;
int sign(double x) {
  if (x < -EPS) return -1;
  if (x > +EPS) return +1;
  return 0;
}

using real = long double;
struct point {
  real x, y;
  point &operator+=(point p) { x += p.x; y += p.y; return *this; }
  point &operator-=(point p) { x -= p.x; y -= p.y; return *this; }
  point &operator*=(real a)     { x *= a;   y *= a;   return *this; }
  point &operator/=(real a)     { return *this *= (1.0/a); }
  point operator-() const    { return {-x, -y}; }
  bool operator<(point p) const {
    int s = sign(x - p.x);
    return s ? s < 0 : sign(y - p.y) < 0;
  }
};
bool operator==(point p, point q) { return !(p < q) && !(q < p); }
bool operator!=(point p, point q) { return p < q || q < p; }
bool operator<=(point p, point q) { return !(q < p); }
point operator+(point p, point q) { return p += q; }
point operator-(point p, point q) { return p -= q; }
point operator*(real a, point p) { return p *= a; }
point operator*(point p, real a) { return p *= a; }
point operator/(point p, real a) { return p /= a; }
real dot(point p, point q) { return p.x*q.x+p.y*q.y; }
real cross(point p, point q) { return p.x*q.y-p.y*q.x; } // left turn > 0
real norm2(point p) { return dot(p,p); }
point orth(point p) { return {-p.y, p.x}; }
real norm(point p) { return sqrt(dot(p,p)); }
real arg(point p) { return atan2(p.y, p.x); }
real arg(point p, point q){ return atan2(cross(p,q), dot(p,q)); }

istream &operator>>(istream &is, point &p) { is>>p.x>>p.y;return is; }
ostream &operator<<(ostream &os, const point &p) { os<<"("<<p.x<<","<<p.y<<")"; return os; }


// exact comparison by polar angle
// usage: sort(all(ps), polar_angle(origin, direction));
struct polar_angle {
  const point o;
  const int s; // +1 for ccw, -1 for cw
  polar_angle(point p = {0,0}, int s = +1) : o(p), s(s) { }
  int quad(point p) const {
    for (int i = 1; i <= 4; ++i, swap(p.x = -p.x, p.y))
      if (p.x > 0 && p.y >= 0) return i;
    return 0;
  }
  bool operator()(point p, point q) const {
    p = p - o; q = q - o;
    if (quad(p) != quad(q)) return s*quad(p) < s*quad(q);
    if (cross(p, q)) return s*cross(p, q) > 0;
    return norm2(p) < norm2(q); // closer first
  }
};

struct line { point p, q; };
bool operator==(line l, line m) {
  return !sign(cross(l.p-l.q,m.p-m.q)) && !sign(cross(l.p-l.q,m.p-l.p));
}

struct segment { point p, q; };
bool operator==(segment l, line m) {
  return (l.p==m.p && l.q==m.q) || (l.p==m.q && l.q==m.p); // do not consider the direction
}
struct circle { point p; real r; };
bool operator==(circle c, circle d) { return c.p == d.p && !sign(c.r - d.r); }

// inner tangent polygon.
// sometimes, circle problem can be reduced to discretized problem
vector<point> discretize(circle c, int n = 50) {
  vector<point> ps;
  for (int i = 0; i < n; ++i) {
    double x = cos(2*PI*i/n), y = sqrt(1 - x*x);
    ps.push_back(c.p + c.r * point({x,y}));
  }
  return ps;
}

typedef vector<point> polygon;


//-----------------------------------------------------------------------------
// visualizer
//-----------------------------------------------------------------------------
struct visualizer {
  real minx, maxx, miny, maxy, scale;
  ofstream ofs;
  visualizer(string s = "data.js") : ofs(s) {
    minx = miny =  1.0/0.0;
    maxx = maxy = -1.0/0.0;
  }
  void update(point p) {
    minx = min(minx, p.x); miny = min(miny, p.y);
    maxx = max(maxx, p.x); maxy = max(maxy, p.y);
    scale = 480/(max(maxx-minx, max(maxy-miny,1.0l)));
  }
  double X(point p) { return scale * (p.x-minx) + 10; }
  double Y(point p) { return 490 - scale * (p.y-miny); }
  vector<point> ps;
  vector<segment> ss;
  vector<circle> cs;
  visualizer &operator<<(circle c) {
    cs.push_back(c);
    update(c.p + point({ c.r, c.r}));
    update(c.p + point({-c.r,-c.r}));
    return *this;
  }
  visualizer &operator<<(point p) {
    ps.push_back(p);
    update(p);
    return *this;
  }
  visualizer &operator<<(segment s) {
    ss.push_back(s);
    update(s.p);
    update(s.q);
    return *this;
  }
  ~visualizer() {
    for (point p: ps)
      ofs << "circle(" << X(p) << "," << Y(p) << ",3.0)" << endl;
    for (segment s: ss)
      ofs << "line(" << X(s.p) << "," << Y(s.p) << ","
                     << X(s.q) << "," << Y(s.q) << ")" << endl;
    for (circle c: cs)
      ofs << "circle(" << X(c.p) << "," << Y(c.p) << "," << scale*c.r << ")" << endl;
  }
};


//-----------------------------------------------------------------------------
// point p is on line l
//-----------------------------------------------------------------------------
vector<point> intersect(line l, point p) {
  if (sign(cross(l.p-p, l.q-p)) != 0) return {};
  return {p};
}
vector<point> intersect(point p, line l) {
  return intersect(l, p);
}

//-----------------------------------------------------------------------------
// point p is on segment s
//-----------------------------------------------------------------------------
vector<point> intersect(segment s, point p) {
  if (sign(cross(s.p-p, s.q-p)) != 0) return {};
  if (sign(  dot(s.p-p, s.q-p))  > 0) return {}; // >  for strict intersect
  return {p};                                    // >= for weak intersect
}
vector<point> intersect(point p, segment s) {
  return intersect(s, p);
}

//-----------------------------------------------------------------------------
// intersection of two lines
//
// Derivation: The crosspoint equation is
//   l.p + (l.q - l.p) t == m.p + (m.q - m.p) s ... (CP.1)
// Taking a cross product with (m.q - m.p), we have
//   (m.q - m.p) x (l.q - l.p) t == (m.q - m.p) x (m.p - l.p)
//   ~~~~~~~~~~~~a~~~~~~~~~~~~      ~~~~~~~~~~~~b~~~~~~~~~~~~
//
// If a != 0, these two lines have a proper crosspoint.
// If a == 0 and b == 0, these two lines are the same line.
// Otherwise, these two lines are parallel.
//-----------------------------------------------------------------------------
vector<point> intersect(line l, line m) {
  auto a = cross(m.q - m.p, l.q - l.p);
  auto b = cross(l.p - m.p, l.q - l.p);
  if ( sign(a)) return {m.p + b/a*(m.q - m.p)}; // properly crossing
  if (!sign(b)) return {m.p, m.q};              // same line
  return {};                                    // disjoint parallel
}

//-----------------------------------------------------------------------------
// intersection of line and segment
//
// In CP.1, t must be in [0, 1]. Thus 0 <= b/a <= 1 is required.
// By assuming a > 0, it is equivalent to b >= 0 and b-a <= 0
//-----------------------------------------------------------------------------
vector<point> intersect(line l, segment s) {
  auto a = cross(s.q - s.p, l.q - l.p);
  auto b = cross(l.p - s.p, l.q - l.p);
  if (a < 0) { a *= -1; b *= -1; }
  if (sign(b) < 0 || sign(a-b) < 0) return {};      // no intersect
  if (sign(a) != 0) return {s.p + b/a*(s.q - s.p)}; // properly crossing
  if (sign(b) == 0) return {s.p, s.q};              // same line
  return {};                                        // disjoint parallel
}

//-----------------------------------------------------------------------------
// intersection of two segments
//
// Solve two CP.1s simultaneously.
// If two segments are overlapping, it is bit difficult.
//-----------------------------------------------------------------------------
vector<point> intersect(segment s, segment t) {
  auto a = cross(s.q - s.p, t.q - t.p);
  auto b = cross(t.p - s.p, t.q - t.p);
  auto c = cross(s.q - s.p, s.p - t.p);
  if (a < 0) { a = -a; b = -b; c = -c; }
  if (sign(b) < 0 || sign(a-b) < 0 ||
      sign(c) < 0 || sign(a-c) < 0) return {};      // disjoint
  if (sign(a) != 0) return {s.p + b/a*(s.q - s.p)}; // properly crossing
  vector<point> ps;                                 // same line
  auto insert_if_possible = [&](point p) {
    for (auto q: ps) if (sign(dot(p-q, p-q)) == 0) return;
    ps.push_back(p);
  };
  if (sign(dot(s.p-t.p, s.q-t.p)) <= 0) insert_if_possible(t.p);
  if (sign(dot(s.p-t.q, s.q-t.q)) <= 0) insert_if_possible(t.q);
  if (sign(dot(t.p-s.p, t.q-s.p)) <= 0) insert_if_possible(s.p);
  if (sign(dot(t.p-s.q, t.q-s.q)) <= 0) insert_if_possible(s.q);
  return ps;
}

//-----------------------------------------------------------------------------
// reflected point with respect to l
//-----------------------------------------------------------------------------
point reflection(point p, line l) {
  auto a = dot(l.p-l.q, l.p-l.q);
  auto b = dot(l.p-p,   l.p-l.q);
  return 2 * (l.p + a/b*(l.q - l.p)) - p;
}

//-----------------------------------------------------------------------------
// closest point on l
//-----------------------------------------------------------------------------
point projection(point p, line l) {
  auto a = dot(l.p-l.q, l.p-l.q);
  auto b = dot(l.p-p,   l.p-l.q);
  return l.p + a/b*(l.q - l.p);
}

//-----------------------------------------------------------------------------
// closest point on s
//-----------------------------------------------------------------------------
point projection(point p, segment s) {
  auto a = dot(s.p-s.q, s.p-s.q);
  auto b = dot(s.p - p, s.p-s.q);
  if (sign(b)   < 0) return s.p;
  if (sign(a-b) < 0) return s.q;
  return s.p + b/a*(s.q - s.p);
}

//-----------------------------------------------------------------------------
// closest point on c
//-----------------------------------------------------------------------------
point projection(point p, circle c) {
  point v = p - c.p;
  return c.p + c.r * v / norm(v);
}
//-----------------------------------------------------------------------------
// distances
//-----------------------------------------------------------------------------
real dist(point p, point q) {
  return norm(p-q);
}
real dist(point p, line l) {
  return dist(p, projection(p, l));
}
real dist(line l, point p) {
  return dist(p, l);
}
real dist(line l, line m) {
  if (sign(cross(l.p-l.q,m.p-m.q))) return 0; // cross
  return dist(l.p, m);
}
real dist(point p, segment s) {
  return dist(p, projection(p, s));
}
real dist(line l, segment s) {
  if (intersect(l, s).size()) return 0;
  return min(dist(l, s.p), dist(l, s.q));
}
real dist(segment s, line l) {
  return dist(l, s);
}
real dist(segment s, segment T) {
  if (intersect(s, T).size()) return 0;
  return min({dist(s.p,T), dist(s.q,T), dist(T.p,s), dist(T.q,s)});
}


//-----------------------------------------------------------------------------
// intersection points of two circles
//-----------------------------------------------------------------------------
vector<point> intersect(circle c, circle d) {
  if (c.r < d.r) swap(c, d);
  auto pow2 = [](real a) { return a*a; };
  auto g = dot(c.p-d.p, c.p-d.p);
  if (sign(g) == 0) {
    if (sign(c.r-d.r) != 0) return {};
    return {{c.p.x+c.r, c.p.y}, {c.p.x, c.p.y+c.r}, {c.p.x-c.r, c.p.y}};
  }
  int inner = sign(g-pow2(c.r-d.r));
  int outer = sign(g-pow2(c.r+d.r));
  if (inner < 0 || outer > 0) return {};
  if (inner == 0) return {(c.r*d.p-d.r*c.p)/(c.r-d.r)};
  if (outer == 0) return {(c.r*d.p+d.r*c.p)/(c.r+d.r)};
  auto eta = (pow2(c.r) - pow2(d.r) + g)/(2*g);
  auto mu = sqrt(pow2(c.r)/g - pow2(eta));
  point q = c.p + eta*(d.p-c.p), v = mu*orth(d.p - c.p);
  return {q + v, q - v};
}

//-----------------------------------------------------------------------------
// intersection of line and circle
//-----------------------------------------------------------------------------
vector<point> intersect(line l, circle c) {
  point u = l.q - l.p, v = l.p - c.p;
  auto a = dot(u,u), b = dot(u,v)/a, t = (dot(v,v) - c.r*c.r)/a;
  auto det = b*b - t;
  if (sign(det) <  0) return {};          // no solution
  if (sign(det) == 0) return {l.p - b*u}; // touch inner/outer
  return {l.p - (b + sqrt(det))*u,        // properly intersect
          l.p - (b - sqrt(det))*u};
}
vector<point> intersect(circle c, line l) {
  return intersect(l, c);
}

//-----------------------------------------------------------------------------
// intersection of segment and circle
// number of points:
//   0 ==> no intersect
//   1 ==> touch
//   2 ==> contained
//   3 ==> penetrating single side
//   4 ==> penetrating both sides
// sorted from s.p to s.q; usually, one would use ans[0] and ans.back().
//-----------------------------------------------------------------------------
vector<point> intersect(circle c, segment s) {
  point u = s.q - s.p, v = s.p - c.p;
  auto a = dot(u,u), b = dot(u,v)/a, t = (dot(v,v) - c.r*c.r)/a;
  auto det = b*b - t;
  if (sign(det) <  0) return {};          // no solution

  auto t1 = -b - sqrt(det), t2 = -b + sqrt(det);
  vector<point> ps;
  auto insert_if_possible = [&](point p) {
    for (auto q: ps) if (sign(dot(p-q, p-q)) == 0) return;
    ps.push_back(p);
  };
  if (sign(c.r - norm(s.p-c.p))   >= 0) insert_if_possible(s.p);
  if (sign(t1) >= 0 && sign(1-t1) >= 0) insert_if_possible(s.p+t1*u);
  if (sign(t2) >= 0 && sign(1-t2) >= 0) insert_if_possible(s.p+t2*u);
  if (sign(c.r - norm(s.q-c.p))   >= 0) insert_if_possible(s.q);
  return ps;
}
vector<point> intersect(segment s, circle c) {
  return intersect(c, s);
}
bool contains(circle c, point p) {
  return sign(dot(c.p-p, c.p-p) - c.r * c.r) <= 0;
}


//-----------------------------------------------------------------------------
// polygon contains a point
// half-line crossing method
//-----------------------------------------------------------------------------
int contains(polygon ps, point p) {
  bool in = false;
  for (int i = 0; i < ps.size(); ++i) {
    int j = (i+1 == ps.size() ? 0 : i+1);
    point a = ps[i] - p, b = ps[j] - p;
    if (a.y > b.y) swap(a, b);
    if (a.y <= 0 && 0 < b.y && cross(a, b) < 0) in = !in;
    if (!sign(cross(a, b)) && sign(dot(a, b)) <= 0)
      return 1; // point on the edge
  }
  return in ? 2 : 0; // point in:out the polygon
}

//-----------------------------------------------------------------------------
// convex hull with the given points
//
// Andrew's monotone chain
//
// [verified]
//
// A. M. Andrew (1979):
// Another efficient algorithm for convex hulls in two dimensions.
// Information Processing Letters 9 (5): 216-219.
//-----------------------------------------------------------------------------
polygon convex_hull(vector<point> ps) {
  int n = ps.size(), k = 0;
  sort(all(ps));
  vector<point> ch(2*n);
  auto cond = [&](point p, point q, point o) {
    int a = sign(cross(p-o, q-o)); // return a < 0 if no three points on line
    return a ? a < 0 : sign(dot(p-o, q-o)) >= 0;
  };
  for (int i = 0; i < n; ch[k++] = ps[i++]) // lower-hull
    for (; k >= 2 && cond(ch[k-2], ch[k-1], ps[i]); --k);
  for (int i = n-2, t = k+1; i >= 0; ch[k++] = ps[i--]) // upper-hull
    for (; k >= t && cond(ch[k-2], ch[k-1], ps[i]); --k);
  ch.resize(k-1);
  return ch;
}


//-----------------------------------------------------------------------------
// farthest pair of points
//
// is located on the convex hull. perform rotating caliper.
//
// [verified ACAC002]
//-----------------------------------------------------------------------------
pair<point,point> farthest_pair(vector<point> ps) {
  vector<point> ch = convex_hull(ps);
  int n = ch.size();
  int u = 0, v = 1;
  real best = -1;
  for (int i = 0, j = 1; i < n; ++i) {
    while (1) {
      int k = (j+1 < n ? j+1 : 0);
      real len = norm2(ch[j] - ch[i]);
      if (sign(len - norm2(ch[k] - ch[i])) <= 0) j = k;
      else {
        if (best < len) { best = len; u = i; v = j; }
        break;
      }
    }
  }
  return make_pair(ch[u], ch[v]);
}


//-----------------------------------------------------------------------------
// signed area of arbitrary polygon
//
// [verified]
//-----------------------------------------------------------------------------
real area(polygon ps) {
  if (ps.size() <= 2) return 0;
  auto a = cross(ps.back(), ps[0]);
  for (int i = 0; i+1 < ps.size(); ++i)
    a += cross(ps[i], ps[i+1]);
  return a/2;
}

//-----------------------------------------------------------------------------
// left side of a convex polygon with respect to a line
//
// [verified]
//-----------------------------------------------------------------------------
polygon convex_cut(polygon ps, line l) {
  vector<point> qs;
  for (int i = 0; i < ps.size(); ++i) {
    int j = (i+1 == ps.size() ? 0 : i+1);
    if (sign(cross(l.p - ps[i], l.q - ps[i])) >= 0) qs.push_back(ps[i]);
    if (sign(cross(l.p - ps[i], l.q - ps[i])) *
        sign(cross(l.p - ps[j], l.q - ps[j])) < 0) {
      auto a = cross(ps[j] - ps[i], l.q - l.p);
      auto b = cross(l.p - ps[i], l.q - l.p);
      qs.push_back(ps[i] + b/a*(ps[j] - ps[i]));
    }
  }
  return qs;
}

//-----------------------------------------------------------------------------
// the circle through three points
//-----------------------------------------------------------------------------
circle three_point_circle(point p, point q, point r) {
  point u = orth(q - p), v = r - p;
  point o = (p+q + u*dot(r-q,v)/dot(u,v))/2;
  return {o, norm(p-o)};
}
//-----------------------------------------------------------------------------
// area of quadrilateral; 
// not verified
//-----------------------------------------------------------------------------
real quadrilateral_area(real a, real b, real c, real d) {
  real s = a+b+c+d;
  return sqrt((s-a)*(s-b)*(s-c)*(s-d))/4;
}
real triangle_area(real a, real b, real c) {
  return quadrilateral_area(a, b, c, 0);
}
real excircle_radius(real a, real b, real c) {
  return a*b*c/4/triangle_area(a, b, c);
}
real incircle_radius(real a, real b, real c) {
  return 2*triangle_area(a,b,c)/(a+b+c);
}

//-----------------------------------------------------------------------------
// tangent line of a circle at the nearest point from p
//
// Let q be a tangent point.
// The angle between q -- p -- c.p is
//   sin(t) = r/|p - c.p|.
// and the solution is
//   p + (c.p - p) * exp(\pm it).
//
// [verified]
//-----------------------------------------------------------------------------
vector<line> tangent(point p, circle c) {
  point u = p - c.p, v = orth(u);
  auto g = dot(u,u) - c.r*c.r;
  if (sign(g) < 0) return {};
  if (sign(g) == 0) return {{p, p + v}};
  u = u * (c.r*c.r/dot(u,u));
  v = v * (c.r*sqrt(g)/dot(v,v));
  return {{p, c.p + u - v}, {p, c.p + u + v}};
}
vector<line> tangent(circle c, point p) {
  return tangent(p, c);
}

//-----------------------------------------------------------------------------
// common tangents of two circles
//
// [verified]
//-----------------------------------------------------------------------------
vector<line> tangent(circle c, circle d) {
  if (c.r < d.r) swap(c, d);
  auto g = dot(c.p-d.p, c.p-d.p);
  if (sign(g) == 0) return {}; // same origin
  point u = (d.p-c.p)/sqrt(g), v = orth(u); // coordinate systems
  vector<line> ls;
  for (int s = +1; s >= -1; s -= 2) {
    auto h = (c.r+s*d.r)/sqrt(g); // = cos(theta)
    if (sign(1 - h*h) == 0) { // touch inner/outer
      ls.push_back({c.p+c.r*u, c.p+c.r*(u+v)});
    } else if (sign(1 - h*h) > 0) { // properly intersect
      point uu = h*u, vv = sqrt(1-h*h)*v;
      ls.push_back({c.p+c.r*(uu+vv), d.p-d.r*(uu+vv)*s});
      ls.push_back({c.p+c.r*(uu-vv), d.p-d.r*(uu-vv)*s});
    }
  }
  return ls;
}

//-----------------------------------------------------------------------------
// common tangent circles of two lines of radius r
//
//   o|o
//  --+--
//   o|o
//
// number of solutions:
//   0: parallel and distance < r
//   1: parallel and distance = r
//   4: non-parallel 
//
// [non-verified]
//-----------------------------------------------------------------------------
vector<circle> tangent_circles(line l, line m, real r) {
  vector<circle> cs;
  real a = cross(l.p-m.p, l.q-l.p), b = cross(m.q-m.p, l.q-l.p);
  if (!sign(b)) { // parallel case
    if (l.q < l.p) swap(l.p, l.q);
    if (m.q < m.p) swap(m.p, m.q);
    if (sign(cross(m.p-l.p, m.q-l.p)) >= 0) swap(l, m); // l is lower
    point v = orth(m.q - m.p); v /= norm(v);
    real d = a / cross(l.q - l.p, v);
    if (sign(d - 2*r) == 0) cs.push_back({l.p + r * v, r});
  } else {
    point o = m.p + a/b * (m.q - m.p), u = l.q - l.p, v = m.q - m.p;
    u /= norm(u); v /= norm(v);
    for (int i = 0; i < 4; ++i) {
      point w = u + v; w /= norm(w);
      real t = r / sqrt(1 - dot(v,w)*dot(v,w));
      cs.push_back({o + t * w, r});
      u *= -1; swap(u, v);
    }
  }
  return cs;
}

//-----------------------------------------------------------------------------
// find a line that passes the maximum number of points
//
// [non-verified]
//-----------------------------------------------------------------------------
int maximum_points_line(vector<point> ps) {
  sort(all(ps)); // assume dx >= 0; if dx == 0 then dy >= 0
  int max_count = 0;
  for (int i = 0; i < ps.size(); ++i) {
    vector<point> qs;
    int base = 1;
    for (int j = 0; j < i; ++j)
      if (ps[j] == ps[i]) ++base;
      else qs.push_back(ps[j] - ps[i]);
    sort(all(qs), polar_angle());
    for (int j = 0, k; j < qs.size(); j += k) {
      for (k = 1; j+k < qs.size() && sign(cross(qs[j], qs[j+k])) == 0; ++k);
      max_count = max(max_count, base + k);
    }
  }
  return max_count;
}

// for comparison
int maximum_points_line_n(vector<point> ps) {
  sort(all(ps));
  int ans = 0;
  for (int i = 0; i < ps.size(); ++i) {
    for (int j = i+1; j < ps.size(); ++j) {
      if (ps[i] == ps[j]) continue;
      int count = 0;
      for (int k = 0; k < ps.size(); ++k) {
        if (intersect((line){ps[i], ps[j]}, ps[k]).size() > 0) ++count;
      }
      ans = max(ans, count);
    }
  }
  return ans;
}
void verify_maximum_points_line() {
  int n = 100;
  vector<point> ps(n);
  for (int i = 0; i < n; ++i) {
    ps[i].x = rand() % 20;
    ps[i].y = rand() % 20;
  }
  cout << maximum_points_line(ps) << endl;
  cout << maximum_points_line_n(ps) << endl;
}


//-----------------------------------------------------------------------------
// triangulate simple polygon in O(n) time.
//
// [non-verified]; future work for polygonal overlay
//-----------------------------------------------------------------------------
real triangulate(vector<point> ps) {
  int n = ps.size();
  vector<int> next(n);
  for (int i = 0; i < n-1; ++i) next[i] = i+1;
  auto is_ear = [&](int i, int j, int k) {
    if (sign(cross(ps[j]-ps[i], ps[k]-ps[i])) <= 0) return false;
    for (int l = next[k]; l != i; l = next[l])
      if (sign(cross(ps[i]-ps[l], ps[j]-ps[l])) >= 0
       && sign(cross(ps[j]-ps[l], ps[k]-ps[l])) >= 0
       && sign(cross(ps[k]-ps[l], ps[i]-ps[l])) >= 0) return false;
    return true;
  };
  real area = 0;
  for (int i = 0; next[next[i]] != i; ) {
    if (is_ear(i, next[i], next[next[i]])) {
      area  += abs(cross(ps[next[i]]-ps[i], ps[next[next[i]]] - ps[i])) / 2;
      next[i] = next[next[i]];
    } else i = next[i];
  }
  return area;
}

//-----------------------------------------------------------------------------
// area of intersection of two circles
// [verified]
//-----------------------------------------------------------------------------
real intersection_area(circle c, circle d) {
  if (c.r < d.r) swap(c, d);
  auto A = [&](real r, real h) {
    return r*r*acos(h/r)-h*sqrt(r*r-h*h);
  };
  auto l = norm(c.p - d.p), a = (l*l + c.r*c.r - d.r*d.r)/(2*l);
  if (sign(l - c.r - d.r) >= 0) return 0; // far away
  if (sign(l - c.r + d.r) <= 0) return PI*d.r*d.r;
  if (sign(l - c.r) >= 0) return A(c.r, a) + A(d.r, l-a);
  else return A(c.r, a) + PI*d.r*d.r - A(d.r, a-l);
}
//-----------------------------------------------------------------------------
// circle-polygon intersection area
// [verified]
//-----------------------------------------------------------------------------
real intersection_area(vector<point> ps, circle c) {
  auto tri = [&](point p, point q){
    point d = q - p;
    auto a = dot(d,p)/dot(d,d), b = (dot(p,p)-c.r*c.r)/dot(d,d);
    auto det = a*a - b;
    if (det <= 0) return arg(p,q) * c.r*c.r / 2;
    auto s = max(0.l, -a-sqrt(det)), t = min(1.l, -a+sqrt(det));
    if (t < 0 || 1 <= s) return c.r*c.r*arg(p,q)/2;
    point u = p + s*d, v = p + t*d;
    return arg(p,u)*c.r*c.r/2 + cross(u,v)/2 + arg(v,q)*c.r*c.r/2;
  };
  auto sum = 0.0;
  for (int i = 0; i < ps.size(); ++i)
    sum += tri(ps[i] - c.p, ps[(i+1)%ps.size()] - c.p);
  return sum;
}
real intersection_area(circle c, vector<point> ps) {
  return intersection_area(ps, c);
}

//-----------------------------------------------------------------------------
// find the closest pair of points by sweepline
// [verified]
//-----------------------------------------------------------------------------
pair<point,point> closest_pair(vector<point> ps) {
  sort(all(ps), [](point p, point q) { return p.y < q.y; });
  auto u = ps[0], v = ps[1];
  auto best = dot(u-v, u-v);
  auto update = [&](point p, point q) {
    auto dist = dot(p-q, p-q);
    if (best > dist) { best = dist; u = p; v = q; }
  };
  set<point> S; S.insert(u); S.insert(v);
  for (int l = 0, r = 2; r < ps.size(); ++r) {
    if (S.count(ps[r])) return {ps[r], ps[r]};
    if ((ps[l].y-ps[r].y)*(ps[l].y-ps[r].y) > best) S.erase(ps[l++]);
    auto i = S.insert(ps[r]).fst;
    for (auto j = i; ; ++j) {
      if (j == S.end() || (i->x-j->x)*(i->x-j->x) > best) break;
      if (i != j) update(*i, *j);
    }
    for (auto j = i; ; --j) {
      if (i != j) update(*i, *j);
      if (j == S.begin() || (i->x-j->x)*(i->x-j->x) > best) break;
    }
  }
  return {u, v};
}
//-----------------------------------------------------------------------------
// find the closest pair of points by divide and conquer
// (slower than the sweepline)
//-----------------------------------------------------------------------------
pair<point,point> closest_pair2(vector<point> ps) {
  sort(all(ps), [](point p, point q) { return p.y < q.y; });
  auto u = ps[0], v = ps[1];
  auto best = dot(u-v, u-v);
  auto update = [&](point p, point q) {
    auto dist = dot(p-q, p-q);
    if (best > dist) { best = dist; u = p; v = q; }
  };
  function<void(int,int)> rec = [&](int l, int r) {
    if (r - l <= 1) {
      for (int i = l; i < r; ++i)
        for (int j = i+1; j < r; ++j)
          update(ps[i], ps[j]);
      stable_sort(&ps[l], &ps[r]);
    } else {
      int m = (l + r) / 2;
      auto y = ps[m].y;
      rec(l, m); rec(m, r);
      inplace_merge(&ps[l], &ps[m], &ps[r]);
      vector<point> qs;
      for (int i = l; i < r; ++i) {
        if ((ps[i].y-y)*(ps[i].y-y) >= best) continue;
        for (int j = (int)qs.size()-1; j >= 0; --j) {
          if ((qs[j].x-ps[i].x)*(qs[j].x-ps[i].x) >= best) break;
          update(qs[j], ps[i]);
        }
        qs.push_back(ps[i]);
      }
    }
  };
  rec(0, ps.size());
  return {u, v};
}


//-----------------------------------------------------------------------------
// half-plane intersection in O(n log n)
// [un-verified; buggy?]
//-----------------------------------------------------------------------------
vector<point> half_plane_intersection(vector<line> ls) {
  int n = ls.size();
  sort(all(ls), [&](line l, line m) { return arg(l.p-l.q) < arg(m.p-m.q); });
  int L = 0, R = 0;
  vector<line> bd(2*n);
  vector<point> ps(2*n);
  bd[R] = ls[0];
  auto left = [&](point p, line l) { return sign(cross(l.p-p, l.q-p)) >= 0; };
  for (int i = 1; i < n; ++i) {
    if (!sign(cross(bd[R].p-bd[R].q, ls[i].p-ls[i].q))) {
      if (left(ls[i].p, bd[R])) bd[R] = ls[i];
    } else {
      while (L < R && !left(ps[R-1], ls[i])) --R;
      while (L < R && !left(ps[L]  , ls[i])) ++L;
      bd[++R] = ls[i];
    }
    if (R > L) ps[R-1] = intersect(bd[R-1], bd[R])[0];
  }
  while (L < R && !left(ps[R-1], bd[L])) --R;
  if (R - L <= 1) return {};
  if (R > L) ps[R] = intersect(bd[L], bd[R])[0];
  vector<point> ch = {ps[L]};
  for (int i = L+1; i <= R; ++i)
    if (!(ch.back() == ps[i])) ch.push_back(ps[i]);
  if (ch.size() > 1 && ch.back() == ch[0]) ch.pop_back();
  return ch;
}



//-----------------------------------------------------------------------------
// geometric median, i.e., p in argmin sum |ps[i] - p|_2
//
// Weiszfield's iterative reweighted method least squares
// with avoiding equal points
// [verified: LA5102]
//-----------------------------------------------------------------------------
point geometric_median(vector<point> ps) {
  point g = {0,0};
  real eta = 1000.0;
  for (int iter = 0; iter < 1000000; ++iter, eta /= 2) {
    real w = 0;
    point h = {0,0};
    for (point p: ps) {
      real a = eta + norm(p - g);
      h = h + p/a; w = w + 1/a;
    }
    h = h / w; swap(g, h);
    if (norm(g - h) < EPS) break;
  }
  return g;
}


//-----------------------------------------------------------------------------
// find a circle of radius r that contains many points as possible
//
// quad-tree search (this is faster than the next sweepline solution)
// [verified: AOJ 1132]
//-----------------------------------------------------------------------------
int maximum_circle_cover(vector<point> ps, double r) {
  const double dx[] = {1,-1,-1,1}, dy[] = {1,1,-1,-1};
  point best_p;
  int best = 0;
  function<void(point,double,vector<point>)>
    rec = [&](point p, double w, vector<point> ps) {
    w /= 2;
    point qs[4];
    vector<point> pss[4];
    for (int i = 0; i < 4; ++i) {
      qs[i] = p + w * point({dx[i], dy[i]});
      int lo = 0;
      for (point q: ps) {
        auto d = dist(qs[i], q);
        if (sign(d - r) <= 0) ++lo;
        if (sign(d - w*sqrt(2) - r) <= 0) pss[i].push_back(q);
      }
      if (lo > best) { best = lo; best_p = qs[i]; }
    }
    for (int i = 0, j; i < 4; ++i) {
      for (int j = i+1; j < 4; ++j)
        if (pss[i].size() < pss[j].size())
          swap(pss[i], pss[j]), swap(qs[i], qs[j]);
      if (pss[i].size() <= best) break;
      rec(qs[i], w, pss[i]);
    }
  };
  real w = 0;
  for (point p: ps) w = max(w, max(abs(p.x), abs(p.y)));
  rec({0,0}, w, ps);
  return best; //best_p;
}

//-----------------------------------------------------------------------------
// find a circle of radius r that contains many points as possible
//
// sweepline O(n^2 log n).
// [verified: AOJ 1132]
//-----------------------------------------------------------------------------
int maximum_circle_cover2(vector<point> ps, double r) {
  int best = 0;
  for (point p: ps) {
    int count = 0;
    vector<pair<double,int>> aux;
    for (point q: ps) {
      auto d = dist(p, q);
      if (sign(d) == 0) ++count;
      else if (sign(d - 2*r) <= 0) {
        auto theta = arg(q-p);
        auto phi   = acos(min(1.l, d/(2*r)));
        aux.push_back({theta-phi, -1});
        aux.push_back({theta+phi, +1});
      }
    }
    sort(all(aux));
    best = max(best, count);
    for (auto a: aux)
      best = max(best, count -= a.snd);
  }
  return best;
}

//-----------------------------------------------------------------------------
// area of union of rectangles
// Bentley's sweepline with segment tree.
//
// [accepted, LightOJ 1120 Rectangle Union]
//-----------------------------------------------------------------------------
struct rectangle { point p, q; }; // lower-left and upper-right
real rectangle_union(vector<rectangle> rs) {
  vector<real> ys; // plane sweep with coordinate compression
  struct event {
    real x, ymin, ymax;
    int add;
    bool operator<(event e) const { return x < e.x; }
  };
  vector<event> es;
  for (auto r: rs) {
    ys.push_back(r.p.y);
    ys.push_back(r.q.y);
    es.push_back({r.p.x, r.p.y, r.q.y, +1});
    es.push_back({r.q.x, r.p.y, r.q.y, -1});
  }
  sort(all(es));
  sort(all(ys));
  ys.erase(unique(all(ys)), ys.end());

  vector<real> len(4 * ys.size()); // segment tree on sweepline
  vector<int> sum(4 * ys.size());
  function<void(real, real, int, int,int,int)> update
    = [&](real ymin, real ymax, int add, int i, int j, int k) {
    ymin = max(ymin, ys[i]); ymax = min(ymax, ys[j]);
    if (ymin >= ymax) return;
    if (ys[i] == ymin && ys[j] == ymax) sum[k] += add;
    else {
      update(ymin, ymax, add, i, (i+j)/2, 2*k+1);
      update(ymin, ymax, add, (i+j)/2, j, 2*k+2);
    }
    if (sum[k]) len[k] = ys[j] - ys[i];
    else        len[k] = len[2*k+1] + len[2*k+2];
  };
  real area = 0;
  for (int i = 0; i+1 < es.size(); ++i) {
    update(es[i].ymin, es[i].ymax, es[i].add, 0, ys.size()-1, 0);
    area += len[0] * (es[i+1].x - es[i].x);
  }
  return area;
}

//-----------------------------------------------------------------------------
// number of points below ps[i]--ps[j]
//
// TODO: めっちゃバグる・バグってるなう
//
// [non-verified]
//-----------------------------------------------------------------------------
struct points_counter {
  int n;
  vector<point> ps;
  vector<vector<int>> low; // low[i][j] is the number of points below (ps[i], ps[j])
  points_counter(vector<point> ps_) : n(ps_.size()), ps(ps_), low(n, vector<int>(n)) {
    double sint = 1e-1, cost = sqrt(1 - sint*sint); // perturbation
    for (point &p: ps) p = {p.x*cost - p.y*sint, p.x*sint + p.y*cost};

    vector<int> is(n); iota(all(is), 0);
    sort(all(is), [&](int i, int j) { return ps[i] < ps[j]; });

    vector<int> left;
    vector<vector<int>> iss(n); // iss[i] = points smaller than i, sorted by the polar angle
    for (int o: is) {
      sort(all(left), [&](int j, int k) { // sort in clock-wise order
        point u = ps[j] - ps[o], v = ps[k] - ps[o];
        int s = sign(cross(u, v));
        return s ? s < 0 : dot(u,u) > dot(v,v);
      });
      iss[o] = left;
      left.push_back(o);
    }
    for (int o: is) {
      cout << "processing " << ps[o] << endl;
      cout << iss[o].size() << endl;
      for (int i = 0; i+1 < iss[o].size(); ++i) {
        int a = iss[o][i], b = iss[o][i+1];
        cout << "  comparing " << ps[a] << " " << ps[b] << endl;
        if (ps[b] < ps[a]) low[b][o] = low[a][o] + low[a][b] + 1;
        else               low[b][o] = low[a][o] - low[a][b];
        low[o][b] = low[b][o];
      }
    }
    cout << low << endl;


    vector<vector<int>> low2(n, vector<int>(n));
    for (int i = 0; i < n; ++i) {
      for (int j = i+1; j < n; ++j) {
        for (int k = 0; k < n; ++k) {
          if (k == i || k == j) continue;
          point p = ps[i], q = ps[j], r = ps[k];
          if (q < p) swap(p, q); // assume p < q
          if (r.x < p.x || r.x >= q.x) continue;
          if (sign(cross(q-r, p-r)) > 0) low2[i][j]++;
        }
        low2[j][i] = low2[i][j];
      }
    }
    cout << low2 << endl;
  }
  int count(int a, int b, int o) { // i -> j -> k ccw order
    if (ps[b] < ps[a]) swap(a, b);
    if (ps[o] < ps[b]) swap(b, o);
    if (ps[b] < ps[a]) swap(a, b);
    if (ps[b].y < ps[a].y) swap(a, b);
    cout << "comparing " << ps[a] << " " << ps[b] << " " << ps[o] << endl;
    cout << "values " << low[o][a] << " " << low[a][b] << " " << low[o][b] << endl;
    return low[o][a] + low[a][b] - low[o][b] + (ps[a].x > ps[b].x);
  }
};


//-----------------------------------------------------------------------------
// segment arrangement
//
// graph whose vertices are the points and edges are the segments.
// complexity is O(m^2 log m).
//
// [accepted: AOJ 1226]
// [accepted: AOJ 2068]
// [accepted: AOJ 0273]
//-----------------------------------------------------------------------------
namespace arrangement_slow {
struct arrangement {
  struct edge {
    int src, dst;
    real weight;
    size_t id, rev;
  };
  int n;
  vector<point> ps;
  map<point,int> id;
  vector<vector<edge>> adj;

  arrangement(vector<segment> ss) : n(0) {
    vector<vector<pair<real, int>>> group(ss.size());
    auto append = [&](int i, point p) {
      if (!id.count(p)) { id[p] = n++; ps.push_back(p); }
      group[i].push_back({norm(ss[i].p - p), id[p]});
    };
    for (int i = 0; i < ss.size(); ++i) {
      append(i, ss[i].p); append(i, ss[i].q);
      for (int j = 0; j < i; ++j) {
        for (point p: intersect(ss[i], ss[j])) {
          append(i, p); append(j, p);
        }
      }
    }
    adj.resize(n);
    for (auto &vs: group) {
      sort(all(vs));
      for (int i = 0; i+1 < vs.size(); ++i) {
        auto u = vs[i].snd, v = vs[i+1].snd;
        if (u == v) continue;
        auto len = vs[i+1].fst - vs[i].fst;
        adj[u].push_back({u, v, len});
        adj[v].push_back({v, u, len});
      }
    }
    // remove duplicates and orient edges
    vector<unordered_map<int, int>> idx(n);
    for (int u = 0; u < n; ++u) {
      auto eq = [&](edge e, edge f) { return e.dst == f.dst; };
      auto lt = [&](edge e, edge f) { return e.dst <  f.dst; };
      sort(all(adj[u]), lt);
      adj[u].erase(unique(all(adj[u]), eq), adj[u].end());
      sort(all(adj[u]), [&](edge e, edge f) {
        return arg(ps[e.dst] - ps[e.src]) > arg(ps[f.dst] - ps[f.src]);
      });
      for (int i = 0; i < adj[u].size(); ++i) {
        adj[u][i].id = i;
        int v = adj[u][i].dst;
        idx[u][v] = i;
        if (idx[v].count(u)) {
          int j = idx[v][u];
          adj[u][i].rev = j;
          adj[v][j].rev = i;
        }
      }
    }
  }
  // traverse utilities for DCEL (planar graph)
  // traversing edges surrounding the left face
  edge twin(edge e) const { return adj[e.dst][e.rev]; }
  edge next(edge e) const {
    int j = adj[e.dst][e.rev].id + 1;
    if (j >= adj[e.dst].size()) j = 0;
    return adj[e.dst][j];
  }
  edge outer() const {
    int s = 0; // leftmost, outerface edge
    for (int i = 1; i < n; ++i) if (ps[i] < ps[s]) s = i;
    for (int i = 0; i < adj[s].size(); ++i) {
      edge e1 = adj[s][i], e2 = adj[s][(i+1)%adj[s].size()];
      if (cross(ps[e1.dst]-ps[s], ps[e2.dst]-ps[s]) <= 0) return e1;
    }
  }

  // [accepted: AOJ 2068 Hide-and-Seek]
  void shortest_path(point sp) {
    int s = id[sp];
    vector<real> dist(n, 1.0/0.0);
    vector<int> prev(n, -1);
    typedef pair<real, int> node;
    priority_queue<node, vector<node>, greater<node>> Q;
    Q.push(node(dist[s] = 0, s));
    while (!Q.empty()) {
      node z = Q.top(); Q.pop();
      if (dist[z.snd] <= z.fst) {
        for (auto e: adj[z.snd]) {
          if (dist[e.dst] > dist[e.src] + e.weight) {
            dist[e.dst] = dist[e.src] + e.weight;
            prev[e.dst] = e.src;
            Q.push({dist[e.dst], e.dst});
          }
        }
      }
    }
    real ans = 0;
    for (int u = 0; u < n; ++u) {
      for (edge e: adj[u]) {
        real s = (dist[e.dst] - dist[e.src] + e.weight)/2;
        ans = max(ans, dist[e.src] + s);
      }
    }
    printf("%.12lf\n", ans);
  }

  // [accepted: AOJ 1226 Fishnet]
  template <class F>
  void traverse(edge e, F func) {
    edge s = e;
    do {
      func(e);
      e = next(e);
    } while (e.src != s.src || e.dst != s.dst);
  }
  void traverse_all_faces() {
    vector<unordered_set<int>> visited(n);
    real ans = 0;
    for (int u = 0; u < n; ++u) {
      for (edge e: adj[u]) {
        if (!visited[e.src].count(e.dst)) {
          real area = 0;
          traverse(e, [&](edge e) {
            visited[e.src].insert(e.dst);
            area += cross(ps[e.src], ps[e.dst]);
          });
          ans = max(ans, area);
        }
      }
    }
    printf("%.6f\n", ans/2);
  }
  void traverse_all_faces2() {
    vector<unordered_set<int>> visited(n);
    real ans = 0;
    for (int u = 0; u < n; ++u) {
      for (edge e: adj[u]) {
        if (!visited[e.src].count(e.dst)) {
          real area = 0;
          traverse(e, [&](edge e) {
            visited[e.src].insert(e.dst);
            area += cross(ps[e.src], ps[e.dst]);
          });
          if (area > 0) ans += area;
        }
      }
    }
    printf("%.12f\n", ans/2);
  }
  edge contains() {
    vector<unordered_set<int>> visited(n);
    auto contains = [&](edge e, point p) {
      edge s = e;
      bool in = false;
      do {
        point a = ps[e.src] - p, b = ps[e.dst] - p;
        if (a.y > b.y) swap(a, b);
        if (a.y <= 0 && 0 < b.y && cross(a, b) < 0) in = !in;
        if (!sign(cross(a, b)) && sign(dot(a, b)) <= 0)
          return 1; // point on the edge
        visited[e.src].count(e.dst);
        e = next(e);
      } while (e.src != s.src || e.dst != s.dst);
      return in ? 2 : 0; // point in:out the polygon
    };
    for (int u = 0; u < n; ++u)
      for (edge e: adj[u])
        if (!visited[e.src].count(e.dst))
          if (contains(e, (point){0,0})) return e;
    return {-1};
  }
};

void AOJ0273() {
  for (int c, w; ~scanf("%d %d", &c, &w) && c; ) {
    vector<point> ps(c);
    for (int i = 0; i < c; ++i) {
      scanf("%lf %lf", &ps[i].x, &ps[i].y);
    }
    vector<segment> ss;
    for (int i = 0; i < w; ++i) {
      int u, v; scanf("%d %d", &u, &v);
      ss.push_back({ps[u-1], ps[v-1]});
    }
    arrangement arr(ss);
    typedef arrangement::edge edge;
    edge e = arr.outer(); // outer-face edge

    vector<unordered_map<int,int>> step(arr.n);
    queue<edge> que;
    auto proceed = [&](edge s, int value) {
      if (step[s.src].count(s.dst)) return;
      edge e = s;
      do {
        step[e.src][e.dst] = value;
        que.push(arr.twin(e));
        e = arr.next(e);
      } while (e.src != s.src || e.dst != s.dst);
    };
    proceed(e, 0);
    int ans = 0;
    while (!que.empty()) {
      edge e = que.front(); que.pop();
      ans = step[e.dst][e.src];
      proceed(e, ans + 1);
    }
    printf("%d\n", ans);
  }
}



void AOJ1226() {
  for (int n; ~scanf("%d", &n) && n; ) {
    vector<vector<double>> a(4, vector<double>(n));
    for (int k = 0; k < 4; ++k)
      for (int i = 0; i < n; ++i)
        scanf("%lf", &a[k][i]);

    vector<segment> ss = {
      {{0,0},{0,1}},
      {{0,1},{1,1}},
      {{1,1},{1,0}},
      {{1,0},{0,0}}
    };
    for (int i = 0; i < n; ++i) {
      ss.push_back({{a[0][i],0},{a[1][i],1}});
      ss.push_back({{0,a[2][i]},{1,a[3][i]}});
    }
    arrangement arr(ss);
    arr.traverse_all_faces();
  }
}

void AOJ2448() {
  int n; scanf("%d", &n);
  vector<point> ps(n);
  for (int i = 0; i < n; ++i)
    scanf("%lf %lf", &ps[i].x, &ps[i].y);
  vector<segment> ss;
  for (int i = 0; i+1 < n; ++i)
    ss.push_back({ps[i], ps[i+1]});
  arrangement arr(ss);
  arr.traverse_all_faces2();
}

void AOJ1247() {
  for (int n; ~scanf("%d", &n) && n; ) {
    vector<segment> ss;
    for (int i = 0; i < n; ++i) {
      double x1, y1, x2, y2;
      scanf("%lf %lf %lf %lf", &x1, &y1, &x2, &y2);
      ss.push_back({{x1,y1},{x2,y2}});
    }
    arrangement arr(ss);
    if (arr.contains().src >= 0) printf("yes\n");
    else                         printf("no\n");
  }
}
}
// ---------------------------------------------------------------------------
// segment arrangement (Bentley-Ottman's plane sweep)
//
// graph whose vertices are the points and edges are the segments.
// complexity is O((n + k) log n); where k is the number of intersections.
//
// Verified: AOJ1226, AOJ2448
// ---------------------------------------------------------------------------
struct arrangement {
  struct Vertex; struct Edge; // Doubly Connected Edge List
  struct Vertex {
    point p;    // 
    Edge *edge; // any edge incident to this vertex
  };
  struct Edge { 
    Vertex *vertex;    // origin vertex of the edge
    Edge *twin;        // reverse edge of e
    Edge *prev, *next; // surrounding edges of face (CCW)
  };
  map<Vertex*, map<Vertex*, Edge*>> adj;
  vector<Vertex*> vertices;
  vector<Edge*> edges;

  vector<segment> segs;

  struct node { // Sweep-Line Structure (RBST)
    int index;
    int size = 1;
    node *left = 0, *right = 0;
  } *root = 0;
  vector<node> ns;

  node *update(node *x) {
    if (x) {
      x->size = 1;
      if (x->left)  x->size += x->left->size;
      if (x->right) x->size += x->right->size;
    }
    return x;
  }
  node *merge(node *x, node *y) {
    if (!x) return y;
    if (!y) return x;
    if (rand() % (x->size + y->size) < x->size) {
      x->right = merge(x->right, y);
      return update(x);
    } else {
      y->left = merge(x, y->left);
      return update(y);
    }
  }
  template <class C> // 3-way split: cond(x) < 0, cond(x) == 0, cond(x) > 0
  tuple<node*, node*, node*> split(node *x, C cond) {
    if (!x) return make_tuple(x,x,x);
    if (cond(x) == 0) {
      auto a = split(x->left, cond);
      auto b = split(x->right, cond);
      x->left = x->right = 0; update(x);
      get<1>(a) = merge(merge(get<1>(a), x), get<1>(b));
      get<2>(a) = get<2>(b);
      return a;
    }
    if (cond(x) < 0) {
      auto a = split(x->right, cond);
      x->right = 0; update(x);
      get<0>(a) = merge(x, get<0>(a));
      return a;
    }
    if (cond(x) > 0) {
      auto a = split(x->left, cond);
      x->left = 0; update(x);
      get<2>(a) = merge(get<2>(a), x);
      return a;
    }
  }
  node *leftmost(node *x) { while (x && x->left) x = x->left; return x; }
  node *rightmost(node *x) { while (x && x->right) x = x->right; return x; }
  template <class F>
  void process(node *x, F func) {
    if (!x) return;
    process(x->left, func);
    func(x);
    process(x->right, func);
  }

  arrangement(vector<segment> segs_) : segs(segs_) {
    ns.resize(segs.size());

    set<point> events;
    map<point, set<int>> L, R;
    for (int i = 0; i < segs.size(); ++i) {
      if (segs[i].q < segs[i].p) swap(segs[i].p, segs[i].q);
      events.insert(segs[i].p);
      events.insert(segs[i].q);
      L[segs[i].p].insert(i);
      R[segs[i].q].insert(i);
      ns[i].index = i;
    }
    vector<Vertex*> last(segs.size());

    while (!events.empty()) {
      const point p = *events.begin();
      events.erase(events.begin());

      Vertex *u = new Vertex({p});
      vertices.push_back(u);
      
      auto cond = [&](node *x) {
        const segment &s = segs[x->index];
        if (sign(s.q.x - s.p.x) == 0) {
          if (sign(p.y - s.p.y) < 0) return -1;
          if (sign(s.q.y - p.y) < 0) return +1;
          return 0;
        }
        return -sign(cross(s.p - p, s.q - p));
      };
      auto z = split(root, cond);
      vector<node*> inserter;
      process(get<1>(z), [&](node *x) {
        Vertex *v = last[x->index];
        if (!adj[u][v]) {
          adj[u][v] = new Edge({u});
          adj[v][u] = new Edge({v});
          adj[u][v]->twin = adj[v][u];
          adj[v][u]->twin = adj[u][v];
          edges.push_back(adj[u][v]);
          edges.push_back(adj[v][u]);
        }
        if (!R[p].count(x->index)) 
          inserter.push_back(x);
      });
      for (int i: L[p]) 
        if (!R[p].count(i))
          inserter.push_back(&ns[i]);

      sort(all(inserter), [&](node *x, node *y) {
        const segment &s = segs[x->index], &t = segs[y->index];
        return sign(cross(s.q - s.p, t.q - t.p)) >= 0;
      });
      auto add_event = [&](node *x, node *y) {
        if (!x || !y) return;
        vector<point> ps = intersect(segs[x->index], segs[y->index]);
        for (point q: ps) 
          if (p < q) events.insert(q);
      };
      if (inserter.empty()) {
        add_event(rightmost(get<0>(z)), leftmost(get<2>(z)));
      } else {
        add_event(rightmost(get<0>(z)), inserter[0]);
        add_event(leftmost(get<2>(z)), inserter.back());
      }
      root = 0;
      for (int i = 0; i < inserter.size(); ++i) {
        last[inserter[i]->index] = u;
        inserter[i]->left = inserter[i]->right = 0;
        root = merge(root, update(inserter[i]));
      }
      root = merge(merge(get<0>(z), root), get<2>(z));
    }
    for (auto &pp: adj) {
      Vertex *u = pp.fst;
      vector<Edge*> es;
      for (auto z: pp.snd) es.push_back(z.snd);
      sort(all(es), [&](Edge *e, Edge *f) {
        auto quad = [](point p) {
          for (int i = 1; i <= 4; ++i, swap(p.x = -p.x, p.y))
            if (p.x > 0 && p.y >= 0) return i;
          return 0;
        };
        const point p = e->twin->vertex->p - e->vertex->p;
        const point q = f->twin->vertex->p - f->vertex->p;
        if (quad(p) != quad(q)) return quad(p) < quad(q);
        return sign(cross(p, q)) > 0;
      });
      u->edge = es.back();
      for (Edge *e: es) {
        u->edge->next = e;
        u->edge->next->prev = u->edge;
        u->edge = u->edge->next;
      }
    }
  }
};
void AOJ1226() {
  for (int n; ~scanf("%d", &n) && n; ) {
    //cout << n << endl;
    vector<vector<double>> a(4, vector<double>(n));
    for (int k = 0; k < 4; ++k)
      for (int i = 0; i < n; ++i)
        scanf("%lf", &a[k][i]);

    vector<segment> ss = {
      {{0,0},{0,1}},
      {{0,1},{1,1}},
      {{1,1},{1,0}},
      {{1,0},{0,0}}
    };
    for (int i = 0; i < n; ++i) {
      ss.push_back({{a[0][i],0},{a[1][i],1}});
      ss.push_back({{0,a[2][i]},{1,a[3][i]}});
    }
    arrangement arr(ss);

    double result = 0;
    set<arrangement::Edge*> seen;
    for (auto *e: arr.edges) {
      if (seen.count(e)) continue;
      auto *f = e;
      double area = 0;
      do {
        seen.insert(f);
        area += cross(f->vertex->p, f->twin->vertex->p);
        f = f->twin->prev;
      } while (f != e);
      result = max(result, area);
    }
    printf("%.6f\n", result/2);
  }
}
void AOJ2448() {
  int n; scanf("%d", &n);
  vector<point> ps(n);
  for (int i = 0; i < n; ++i)
    scanf("%lf %lf", &ps[i].x, &ps[i].y);
  vector<segment> ss;
  for (int i = 0; i+1 < n; ++i)
    ss.push_back({ps[i], ps[i+1]});
  arrangement arr(ss);

  double result = 0;
  set<arrangement::Edge*> seen;
  for (auto *e: arr.edges) {
    if (seen.count(e)) continue;
    auto *f = e;
    double area = 0;
    do {
      seen.insert(f);
      area += cross(f->vertex->p, f->twin->vertex->p);
      f = f->twin->prev;
    } while (f != e);
    if (area > 0) result += area;
  }
  printf("%.12f\n", result/2);
}


//-----------------------------------------------------------------------------
// visibility graph in O(|points|^2 |obstacles|)
//
// Here, the visibility of p and q is:
//   1. p, q in [a,b]          ==> visible
//   2. a in (p,q)             ==> invisible
//   3. (p,q) intersects (a,b) ==> invisible
//   4. m=(p+q)/2 in int(obs)  ==> invisible
//
// [non-verified]; 
//-----------------------------------------------------------------------------
struct visibility_graph {
  struct edge {
    int src, dst;
    real weight;
  };
  int n;
  vector<point> ps;
  vector<vector<edge>> adj;
  visibility_graph(vector<point> ps, vector<polygon> os) :
    n(ps.size()), ps(ps), adj(n) {
    auto blocked = [&](point s, point t, polygon &obs) {
      bool in = false, on = false;
      point m = (s + t) / 2;
      for (int i = 0; i < obs.size(); ++i) {
        int j = (i+1 == obs.size() ? 0 : i+1);
        point a = obs[i], b = obs[j];
        if (!sign(cross(a-s,b-s)) && sign(dot(a-s,b-s)) <= 0 && // s is on [a,b] and t is on [a,b]
            !sign(cross(a-t,b-t)) && sign(dot(a-t,b-t)) <= 0) return false;
        if (!sign(cross(s-a,t-a)) && sign(dot(s-a,t-a)) < 0) return false; // (s,t) contains a
        if (sign(cross(a-s,b-s))*sign(cross(a-t,b-t)) < 0 && // (s,t) intersects (a,b)
            sign(cross(s-a,t-a))*sign(cross(s-b,t-b)) < 0) return true;
        if (b.y < a.y) swap(a, b); // a is lower
        if (a.y <= m.y && m.y < b.y)
          if (sign(cross(a-m, b-m)) < 0) in = !in;
      }
      return in; // midpoint is on obs
    };
    for (int i = 0; i < ps.size(); ++i) {
      for (int j = i+1, k; j < ps.size(); ++j) {
        for (k = 0; k < os.size(); ++k)
          if (blocked(ps[i], ps[j], os[k])) break;
        if (k == os.size()) {
          auto len = norm(ps[i] - ps[j]);
          adj[i].push_back({i, j, len});
          adj[j].push_back({j, i, len});
        }
      }
    }
  }
};
//
// weighted geometric shortest path
//
// find a shortest path on 2D plane,
// where each obstacle is assigned a weight (inverse velocity).
// if weight = INF, it reduced to the standard geometric shortest path problem.
//
// basic algorithm: A* search on implicit visibility graph.
// It performs O(n^2) A* search with dist(u,t) as a heuristics.
//
// on a weighted problem, divide edges to approximate solutions.
//
// reference:
// L. Aleksandrov, A. Maheshwari, and J. -R. Sack (1998):
// Determining approximate shortest paths on weighted polyhedral surfaces.
//
// see aso: http://arxiv.org/pdf/0904.2550.pdf
//
struct geometric_shortest_path {
  vector<point> ps;
  vector<int> prev, next;
  vector<vector<int>> region;
  vector<real> weight;
  void add_region(vector<point> poly, real w = 1.0/0.0) {
    int n = ps.size(), k = poly.size();
    vector<int> obj;
    for (int i = 0; i < k; ++i) {
      ps.push_back(poly[i]);
      obj.push_back(n+i);
      next.push_back(n+i+1);
      prev.push_back(n+i-1);
    }
    next[n+k-1] = n;
    prev[n] = n+k-1;
    region.push_back(obj);
    weight.push_back(w);
  }
  void refine(int v) {
    auto divide = [&](int u, int w) {
      int v = ps.size();
      ps.push_back((ps[u] + ps[w])/2);
      prev.push_back(u);
      next.push_back(w);
      prev[next[v]] = v;
      next[prev[v]] = v;
    };
    divide(prev[v], v);
    divide(v, next[v]);
  }
  enum { INTERSECT, CONTAINED, ONLINE, VISIBLE };
  int visibility(point s, point t, vector<int> &reg) {
    bool in = false;
    point m = (s + t) / 2;
    for (int i = 0; i < reg.size(); ++i) {
      point a = ps[reg[i]], b = ps[reg[(i+1)%reg.size()]];
      if (!sign(cross(a-s,b-s)) && sign(dot(a-s,b-s)) <= 0 && // s is on [a,b] and t is on [a,b]
          !sign(cross(a-t,b-t)) && sign(dot(a-t,b-t)) <= 0) return ONLINE;
      if (sign(cross(a-s,b-s))*sign(cross(a-t,b-t)) < 0 && // (s,t) intersects (a,b)
          sign(cross(s-a,t-a))*sign(cross(s-b,t-b)) < 0) return INTERSECT;
      if (b.y < a.y) swap(a, b);
      if (a.y <= m.y && m.y < b.y && sign(cross(a-m, b-m)) < 0) in = !in;
    }
    return in ? CONTAINED : VISIBLE;
  }
  real length(int u, int v) {
    point s = ps[u], t = ps[v];
    auto len = norm(s - t);
    for (int i = 0; i < region.size(); ++i) {
      int a = visibility(s, t, region[i]);
      if (a == ONLINE) return len;
      if (a == CONTAINED) return weight[i] * len;
      if (a == INTERSECT) return 1.0 / 0.0;
    }
    return len;
  }
  pair<real, vector<int>> shortest_path(point p, point q) {
    ps.push_back(p); ps.push_back(q);
    int n = ps.size(), s = n-2, t = n-1;
    vector<real> dist(n, 1.0/0.0), h(n); dist[s] = 0;
    for (int u = 0; u < n; ++u) h[u] = norm(ps[u]-ps[t]); // modify here if weight < 1.0
    vector<int> last(n, ~s); // negative -> unfixed
    while (last[t] < 0) {
      int u = t;
      for (int v = 0; v < n; ++v)
        if (last[v] < 0 && dist[v] + h[v] < dist[u] + h[u]) u = v;
      last[u] = ~last[u];
      for (int v = 0; v < n; ++v) {
        if (last[v] >= 0) continue;
        auto len = dist[u] + length(u, v);
        if (sign(dist[v] - len) > 0) {
          dist[v] = len;
          last[v] = ~u;
        }
      }
    }
    auto cost = dist[t];
    vector<int> path;
    if (cost < 1.0 / 0.0) {
      while (t != s) {
        path.push_back(t);
        t = last[t];
      }
      path.push_back(s);
    }
    ps.pop_back(); ps.pop_back();
    return {cost, path};
  }
};

//-----------------------------------------------------------------------------
// merge segments in O(n log n)
// [non-verified]
//-----------------------------------------------------------------------------
vector<segment> merge_segments(vector<segment> ss) {
  auto compare = [&](segment s, segment t) {
    int a = sign(cross(s.q-s.p, t.p-t.q));
    return a ? a < 0 : sign(cross(t.p-s.p, t.q-s.p)) < 0;
  };
  for (segment &s: ss) if (s.q < s.p) swap(s.p, s.q);
  sort(all(ss), compare);
  vector<segment> res;
  for (int i = 0, j; i < ss.size(); i = j) {
    for (j = i+1; j < ss.size(); ++j)
      if (compare(ss[i], ss[j])) break;
    point o = ss[i].p, v = ss[i].q - ss[i].p;
    sort(ss.begin()+i, ss.begin()+j, [&](segment s, segment &t) {
      auto a = dot(s.p - t.p, v);
      if (a) return a < 0;
      return dot(s.q - t.q, v) < 0;
    });
    for (int a = i, b; a < j; a = b) {
      for (b = a+1; b < j; ++b) {
        if (sign(dot(ss[a].q - ss[b].p, v)) < 0) break;
        if (dot(ss[a].q - ss[b].q, v) < 0) ss[a].q = ss[b].q;
      }
      res.push_back(ss[a]);
    }
  }
  return res;
}

//-----------------------------------------------------------------------------
// relative neighborhood graph
//   there is an edge between p and q iff
//     d(p,q) < d(p,r), and d(p,q) < d(q,r)
//   for all other point r.
//
// it contains euclidean mst.
//
// O(n^2) by Katajainen and Nevalainen
//
// [non-verified]
//-----------------------------------------------------------------------------
struct relative_neighborhood_graph {
  struct edge {
    int src, dst;
    real len;
  };
  int n;
  vector<point> ps;
  vector<vector<edge>> adj;
  relative_neighborhood_graph(vector<point> ps) :
    n(ps.size()), ps(ps), adj(n) {

    auto quadrant = [&](point p) {
      for (int i = 0; i < 8; ++i) {
        if (p.x >= 0 && p.y >= 0 && p.y < p.x) return i;
        p = (point){p.x + p.y, -p.x + p.y};
      }
      return -1;
    };
    for (int i = 0; i < n; ++i) {
      vector<vector<int>> cand(8);
      vector<real> dist(8, 1.0/0.0);
      for (int j = 0; j < n; ++j) {
        if (i == j) continue;
        int k = quadrant(ps[j] - ps[i]);
        real d = norm2(ps[j] - ps[i]);
        if      (sign(d - dist[k])  < 0) cand[k] = {j}, dist[k] = d;
        else if (sign(d - dist[k]) == 0) cand[k].push_back(j);
      }
      for (int k = 0, l; k < 8; ++k) {
        //cout << dist[k] << endl;
        for (int j: cand[k]) {
          if (j >= i) continue;
          for (l = 0; l < n; ++l) {
            if (l == i || l == j) continue;
            if (sign(dist[k] - norm2(ps[l] - ps[i])) > 0) break;
          }
          if (l == n) {
            adj[i].push_back({i, j, sqrt(dist[k])});
            adj[j].push_back({j, i, sqrt(dist[k])});
          }
        }
      }
    }
  }
};


// ---------------------------------------------------------------------------
// maximum area of empty convex k-gon
//
// David P. Dovkin, Harbert Edelsbrunner, and Mark H. Overmars (1990):
// Searching for Empty Convex Polygons.
// Algorithmica, vol.5, no.1, pp.561--571.
//
// Complexity: O(n^3)
//
// [non-verified]
// ---------------------------------------------------------------------------
real maximum_area_empty_convex_polygon(vector<point> ps) {
  int n = ps.size();
  real ans = 0;
  for (int o = 0; o < n; ++o) {
    vector<point> vs;
    for (point p: ps)
      if (p < ps[o]) vs.push_back(p - ps[o]);
    int m = vs.size();
    sort(all(vs), [&](point p, point q) {
      int s = sign(cross(p, q));
      return s ? s > 0 : sign(dot(q,q) - dot(p,p)) > 0;
    });
    vector<vector<real>> dp(m, vector<real>(m));
    for (int k = 1; k < m; ++k) {
      int i = k-1;
      while (i >= 0 && !sign(cross(vs[i], vs[k]))) --i;
      vector<int> seq;
      for (int j; i > 0; i = j) {
        seq.push_back(i);
        for (j = i-1; j > 0 && sign(cross(vs[i]-vs[k], vs[j]-vs[k])) > 0; --j);
        auto v = cross(vs[i], vs[k]);
        if (j >= 0) v += dp[i][j];
        dp[k][i] = v;
        ans = max(ans, v);
      }
      for (int i = (int)seq.size()-2; i >= 0; --i)
        dp[k][seq[i]] = max(dp[k][seq[i]], dp[k][seq[i+1]]);
    }
  }
  return ans;
}

//-----------------------------------------------------------------------------
// delaunay triangulation
// O(n^2) in the worst case, O(n log n) -- O(n sqrt(n)) in pratice.
//
// [verified: RUPC 2013 F]
//-----------------------------------------------------------------------------
struct delaunay {
  struct edge {
    int src, dst;
  };
  int n;
  vector<point> ps;
  vector<vector<edge>> adj; // optional
  vector<int> inner;
  int incircle(int a, int b, int c, int p) {
    point u = ps[a]-ps[p], v = ps[b]-ps[p], w = ps[c]-ps[p];
    return sign(norm2(u)*cross(v,w)
               +norm2(v)*cross(w,u)
               +norm2(w)*cross(u,v)) > 0;
  }
  bool orient(int a, int b, int p) {
    point u = ps[a]-ps[b], v = ps[p]-ps[b];
    int s = sign(cross(u, v));
    return s ? s > 0 : sign(dot(u, v)) > 0;
  }
  delaunay(vector<point> ps) : n(ps.size()), ps(ps), adj(n), inner(n) {
    if (n <= 1) return;
    vector<unordered_map<int,int>> ccw(n); // ccw[u][v] is the third pt for (u,v)
    auto make_triangle = [&](int a, int b, int c) {
      ccw[a][b] = c; ccw[b][c] = a; ccw[c][a] = b;
    };
    vector<int> is(n); iota(all(is), 0);
    sort(all(is), [&](int i, int j) { return ps[i] < ps[j]; });
    // delaunay flips
    function<void(int,int)> rec = [&](int a, int b) {
      if (!ccw[a].count(b) || !ccw[b].count(a)) return;
      int c = ccw[a][b], d = ccw[b][a];
      if (incircle(a, b, c, d) > 0) {
        ccw[a].erase(b); ccw[b].erase(a);
        make_triangle(d, c, a);
        make_triangle(c, d, b);
        rec(a, d); rec(d, b); rec(b, c); rec(c, a);
      }
    };
    // lexicographic triangulation
    vector<int> next(n,-1), prev(n,-1);
    next[is[0]] = prev[is[0]] = is[1];
    next[is[1]] = prev[is[1]] = is[0];
    for (int i = 2; i < n; ++i) {
      int h = is[i], u = is[i-1], v = u;
      while ( orient(u, next[u], h)) u = next[u];
      while (!orient(v, prev[v], h)) v = prev[v];
      for (int w = v; w != u; w = next[w])
        if (sign(cross(ps[next[w]]-ps[h], ps[w]-ps[h])) > 0)
          make_triangle(w, h, next[w]);
      next[h] = u; prev[u] = h;
      prev[h] = v; next[v] = h;
    }
    for (int u: is) {
      auto nbh = ccw[u]; // hardcopy
      for (auto z: nbh) rec(z.fst, z.snd); // flip
    }
    // complete graph structure
    for (int u: is) {
      int v = ccw[u].begin()->fst, s = v;
      while (ccw[s].count(u)) {
        s = ccw[s][u];
        if (s == v) break;
      } // v != s means that u is on the outer face
      if (v != s) { inner[u] = false; v = s; }
      do {
        adj[u].push_back({u, v});
        if (!ccw[u].count(v)) break;
        v = ccw[u][v];
      } while (v != s);
    }
  }
};
//-----------------------------------------------------------------------------
// voronoi diagram (from delaunay triangulation)
// O(sum deg^2) = O(n^2) in worst case, O(n) for random input.
//
// [verified: RUPC 2013 F]
//-----------------------------------------------------------------------------
struct voronoi {
  struct edge {
    int src, dst;
    real len;
  };
  int n, m;
  vector<point> ps, qs; // qs is the voronoi vertices
  map<point,int> id;
  vector<vector<int>> cell;
  vector<vector<edge>> adj;
  void add_edge(int u, int v) {
    auto len = norm(qs[u] - qs[v]);
    adj[u].push_back({u, v, len});
    adj[v].push_back({v, u, len});
  }
  int node(point p) {
    if (!id.count(p)) { id[p] = m++; qs.push_back(p); adj.push_back({}); }
    return id[p];
  }
  voronoi(delaunay DT, vector<point> domain) :
    n(DT.n), m(0), ps(DT.ps), cell(n) {
    for (int u = 0; u < n; ++u) {
      vector<point> region = domain;
      for (auto e: DT.adj[u]) {
        point s = (ps[e.src]+ps[e.dst])/2, d = orth(ps[e.dst]-ps[e.src]);
        region = convex_cut(region, {s, s+d});
      }
      for (int i = 0; i < region.size(); ++i) {
        add_edge(node(region[i]), node(region[(i+1)%region.size()]));
        cell[u].push_back(node(region[i]));
      }
    }
  }
  voronoi(vector<point> ps, vector<point> domain) :
    voronoi(delaunay(ps), domain) { }
};

/* 
DCEL にしたい

struct Delaunay {
  struct Vertex; struct Edge;
  struct Vertex {
    point p;
    Edge *edge; // any adjacent edge
  };
  struct Edge {
    Vertex *vertex; // origin
    Edge *twin;
    Edge *prev, *next;
  };
  int incircle(Vertex *a, Vertex *b, Vertex *c, Vertex *p) {
    point u = a->p - p->p, v = b->p - p->p, w = c->p - p->p;
    return sign(norm2(u)*cross(v,w)
               +norm2(v)*cross(w,u)
               +norm2(w)*cross(u,v)) > 0;
  }
  bool orient(Vertex *a, Vertex *b, Vertex *p) {
    point u = a->p - b->p, v = p->p - b->p;
    int s = sign(cross(u, v));
    return s ? s > 0 : sign(dot(u, v)) > 0;
  }
  Delaunay(vector<point> ps) {
    
  }
      z   d   w            d  
        /   \            / | \
    枝 a -e-> b   ==>   a  |  b
        \   /            \ v /
      y   c   x            c

    を逆転させる手続き
    if incircle:
      r = e->twin
      x = e->next y = e->prev
      z = r->next w = r->prev
      x->vertex->edge = x
      z->vertex->edge = z
      w->next = x x->next = r r->next = w
      w->prev = r x->prev = w r->prev = x
      y->next = z z->next = e e->next = y
      y->prev = e z->prev = y e->prev = z
      e->vertex = w->vertex
      r->vertex = y->vertex
      rec(x); rec(y); rec(z); rec(w);
    
    Vertex *u :: rightmost
    Edge *e = u->edge
    while  orient: e = e->prev
    
    span(e->next)
    e->twin->next で CH を traverse できる

    while !orient: 
      r = new Edge, s = new Edge, t = new Edge
      r->next = s s->next = t t->next = r
      r->prev = b 
    


      p  ---- v
       \     /    
         q /

    auto make_triangle = [&](int a, int b, int c) {
      ccw[a][b] = c; ccw[b][c] = a; ccw[c][a] = b;
    };
    vector<int> is(n); iota(all(is), 0);
    sort(all(is), [&](int i, int j) { return ps[i] < ps[j]; });
    // delaunay flips
    function<void(int,int)> rec = [&](int a, int b) {
      if (!ccw[a].count(b) || !ccw[b].count(a)) return;
      int c = ccw[a][b], d = ccw[b][a];
      if (incircle(a, b, c, d) > 0) {
        ccw[a].erase(b); ccw[b].erase(a);
        make_triangle(d, c, a);
        make_triangle(c, d, b);
        rec(a, d); rec(d, b); rec(b, c); rec(c, a);
      }
    };

};
*/


//-----------------------------------------------------------------------------
// [non-verified]
//-----------------------------------------------------------------------------
bool is_congruence(polygon ps, polygon qs) {
  if (ps.size() != qs.size()) return false;
  int n = ps.size();
  vector<point> as, bs;
  for (int i = 0; i < n; ++i) {
    int j = (i+1 < n ? i+1 : 0), k = (j+1 < n ? j+1 : 0);
    as.push_back({dot(ps[i]-ps[j], ps[k]-ps[j]),
                cross(ps[i]-ps[j], ps[k]-ps[j])});
    bs.push_back({dot(qs[i]-qs[j], qs[k]-qs[j]),
                cross(qs[i]-qs[j], qs[k]-qs[j])});
  }
  vector<int> fail(n+1, -1); // knuth morris pratt
  for (int i = 1, j = -1; i <= n; ++i) {
    while (j >= 0 && !(as[j] == as[i-1])) j = fail[j];
    fail[i] = ++j;
  }
  for (int i = 0, k = 0; i < 2*n; ++i) {
    while (k >= 0 && !(bs[i%n] == as[k])) k = fail[k];
    if (++k == n) return true; // ps[0,n) == qs[i-n+1, *)
  }
  return false;
}
//-----------------------------------------------------------------------------
// [non-verified]
//-----------------------------------------------------------------------------
bool is_similar(polygon ps, polygon qs) {
  if (ps.size() != qs.size()) return false;
  int n = ps.size();
  vector<point> as, bs;
  for (int i = 0; i < n; ++i) {
    int j = (i+1 < n ? i+1 : 0), k = (j+1 < n ? j+1 : 0);
    as.push_back({dot(ps[i]-ps[j], ps[k]-ps[j]),
                cross(ps[i]-ps[j], ps[k]-ps[j])});
    bs.push_back({dot(qs[i]-qs[j], qs[k]-qs[j]),
                cross(qs[i]-qs[j], qs[k]-qs[j])});
  }
  real maxa = as[0].x, maxb = bs[0].y;
  for (int i = 1; i < n; ++i)
    maxa = max(maxa, as[i].x), maxb = max(maxb, bs[i].x);
  vector<int> fail(n+1, -1); // knuth morris pratt
  for (int i = 1, j = -1; i <= n; ++i) {
    while (j >= 0 && !(as[j] == as[i-1])) j = fail[j];
    fail[i] = ++j;
  }
  for (int i = 0, k = 0; i < 2*n; ++i) {
    while (k >= 0 && !(maxa*bs[i%n] == maxb*as[k])) k = fail[k];
    if (++k == n) return true; // maxb * ps[0,n) == maxa * qs[i-n+1, *)
  }
  return false;
}
//-----------------------------------------------------------------------------
// [non-verified]
//-----------------------------------------------------------------------------
polygon simplify(polygon ps) {
  int n = ps.size();
  polygon qs;
  for (int i = 0; i < n; ++i) {
    int j = (i+1 < n ? i+1 : 0), k = (j+1 < n ? j+1 : 0);
    if (sign(dot(ps[j]-ps[i], ps[k]-ps[i])) < 0 ||
        sign(cross(ps[j]-ps[i], ps[k]-ps[i])))
      qs.push_back(ps[j]);
  }
  return qs;
}


//-----------------------------------------------------------------------------
// VP tree (metric tree)
//
// [non-verified]
//-----------------------------------------------------------------------------
struct vantage_point_tree {
  int n;
  vector<point> ps;
  vector<pair<real, int>> aux;
  vantage_point_tree(vector<point> ps) : n(ps.size()), ps(ps), aux(n) {
    for (int i = 0; i < n; ++i) aux[i].snd = i;
    function<void(int,int)> rec = [&](int l, int r) {
      if (l == r) return;
      swap(aux[l], aux[l+rand()%(r-l)]);
      for (int i = l+1; i < r; ++i)
        aux[i].fst = norm(ps[aux[i].snd] - ps[aux[l].snd]);
      int m = (l+1 + r) / 2;
      nth_element(aux.begin()+l+1, aux.begin()+m, aux.begin()+r);
      aux[l].fst = aux[m].fst;
      rec(l+1, m); rec(m, r);
    };
    rec(0, n);
  }
  vector<int> closest(point p, int k = 1) {
    priority_queue<pair<real, int>> topk;
    function<void(int,int)> rec = [&](int l, int r) {
      int m = (l+1 + r) / 2;
      if (l == r) return;
      auto d = norm(p - ps[aux[l].snd]);
      if (topk.size() < k) topk.push({d, aux[l].snd});
      else if (topk.top().fst > d) {
        topk.pop();
        topk.push({d, aux[l].snd});
      }
      if (d <= aux[l].fst) {
        rec(l+1, m);
        if (aux[l].fst - d < topk.top().fst) rec(m, r);
      } else {
        rec(m, r);
        if (d - aux[l].fst < topk.top().fst) rec(l+1, m);
      }
    };
    rec(0, n);
    vector<int> ans;
    for (; !topk.empty(); topk.pop()) ans.push_back(topk.top().snd);
    reverse(all(ans));
    return ans;
  }
};

//-----------------------------------------------------------------------------
// random ball over
//   O(n f(n)) construction, O(n/f(n)) query
//   If there are T queries, set f(n) = O(sqrt(T)).
//
// [non-verified]
//-----------------------------------------------------------------------------
struct random_ball_cover {
  int n, m;
  vector<point> ps, cs;
  vector<vector<pair<real, int>>> vs;
  random_ball_cover(vector<point> ps) : n(ps.size()), ps(ps), cs(ps) {
    for (int k = m = 1; m < n && k < n; k *= 2, ++m); // m = O(log(n))
    cs.resize(m); vs.resize(m);
    for (int u = 0; u < n; ++u) {
      int j = 0;
      real best = norm(cs[j] - ps[u]);
      for (int k = 1; k < m; ++k) {
        real len = norm(cs[k] - ps[u]);
        if (best > len) { best = len; j = k; }
      }
      vs[j].push_back({best, u});
    }
    for (int j = 0; j < m; ++j) {
      sort(all(vs[j])); reverse(all(vs[j]));
    }
  }
  vector<int> closest(point p, int k = 1) {
    vector<pair<real, int>> order(m);
    for (int j = 0; j < m; ++j)
      order[j] = {norm(cs[j] - p), j};
    sort(all(order));
    priority_queue<pair<real, int>> topk;
    for (auto ord: order) {
      real rad = ord.fst;
      for (auto z: vs[ord.snd]) {
        if (topk.size() < k || topk.top().fst > ord.fst - z.fst) {
          topk.push({norm(ps[z.snd] - p), z.snd});
          if (topk.size() > k) topk.pop();
        } else break;
      }
    }
    vector<int> ans;
    for (; !topk.empty(); topk.pop()) ans.push_back(topk.top().snd);
    reverse(all(ans));
    return ans;
  }
};


//-----------------------------------------------------------------------------
// Dynamic Convex Hull (Overmars-Leeuwen)
//
// Memoise divide-and-conquer procedure in segment tree.
// Each tree nodes stores convex hulls by binary search tree.
//
// [verified] CODECHEF MGCHGEOM, ACM/ICPC 2017-D
//-----------------------------------------------------------------------------

template <int turn>
struct half_hull {
  struct node {
    point p;
    node *child[2], *parent;
    node *next[2];

    real area;
    double perimeter;
  };
  node *update(node *t) {
    if (t) {
      t->area = 0;
      if (t->child[0]) t->area += t->child[0]->area + cross(t->p, t->next[0]->p);
      if (t->child[1]) t->area += t->child[1]->area - cross(t->p, t->next[1]->p);

      t->perimeter = 0;
      if (t->child[0]) t->perimeter += t->child[0]->perimeter + norm(t->p - t->next[0]->p);
      if (t->child[1]) t->perimeter += t->child[1]->perimeter + norm(t->p - t->next[1]->p);
    }
    return t;
  }
  int dir(node *x, node *y) { return x && x->child[0] == y ? 0 : 1; }
  void link(node *x, node *y, int d) {
    if (x) x->child[d] = y;
    if (y) y->parent = x;
  }
  void rot(node *x, int d) {
    node *y = x->child[d], *z = x->parent;
    link(x, y->child[!d], d);
    link(y, x, !d);
    link(z, y, dir(z, x));
    update(x); update(y);
  }
  void splay(node *x) {
    if (!x) return;
    while (x->parent) {
      node *y = x->parent, *z = y->parent;
      int dy = dir(y, x), dz = dir(z, y);
      if (!z)            { rot(y, dy); }
      else if (dy == dz) { rot(z, dz), rot(y, dy); }
      else               { rot(y, dy), rot(z, dz); }
    }
  }
  int n;
  vector<node> vs;
  vector<node*> ch, sh; 
  vector<point> ps;
  half_hull(vector<point> ps_) : ps(ps_) {
    sort(all(ps));
    ps.erase(unique(all(ps)), ps.end());
    n = ps.size();
    vs.resize(n); 
    for (int i = 0; i < n; ++i) vs[i] = {ps[i]};
    ch.resize(4*n); sh.resize(4*n);
  }
  pair<node*, node*> bridge(node *lch, node *rch) {
    if (lch && rch) {
      while (1) {
        int it = 0;
        for (splay(rch); ; ++it) {
          if (rch->next[0] && turn*sign(cross(rch->next[0]->p - lch->p, rch->p - lch->p)) <= 0) {
            rch = rch->child[0];
          } else if (rch->next[1] && turn*sign(cross(rch->p - lch->p, rch->next[1]->p - lch->p)) > 0) {
            rch = rch->child[1];
          } else break;
        }
        for (splay(lch); ; ++it) {
          if (lch->next[1] && turn*sign(cross(lch->p - rch->p, lch->next[1]->p - rch->p)) <= 0) {
            lch = lch->child[1];
          } else if (lch->next[0] && turn*sign(cross(lch->next[0]->p - rch->p, lch->p - rch->p)) > 0) {
            lch = lch->child[0];
          } else break;
        }
        if (it == 0) break;
      }
    } else if (lch) {
      for (splay(lch); lch->child[1]; lch = lch->child[1]);
    } else if (rch) {
      for (splay(rch); rch->child[0]; rch = rch->child[0]);
    }
    return {lch, rch};
  }
  node *apart(node *x, int d) {
    if (!x) return 0;
    node *y = x->child[d];
    if (!y) return 0;
    y->parent = x->child[d] = 0;
    for (; y->child[!d]; y = y->child[!d]);
    splay(y);
    x->next[d] = y->next[!d] = 0;
    update(x);
    return y;
  }
  node *join(node *x, node *y) {
    if (!x) return y;
    if (!y) return x;
    x->child[1] = x->next[1] = y;
    y->parent = y->next[0] = x;
    update(x);
    return x;
  }
  template <int I>
  void build(int l, int r, int k, int i) {
    if (l+1 == r) { ch[k] = I ? &vs[i] : 0; return; }
    splay(ch[2*k+1]); 
    apart(ch[2*k+1], 1);
    ch[2*k+1] = join(ch[2*k+1], sh[2*k+1]);
    splay(ch[2*k+2]);
    ch[2*k+2] = join(sh[2*k+2], ch[2*k+2]);
    if (i < (l+r)/2) build<I>(l, (l+r)/2, 2*k+1, i);
    else             build<I>((l+r)/2, r, 2*k+2, i);
    tie(ch[2*k+1], ch[2*k+2]) = bridge(ch[2*k+1], ch[2*k+2]);
    splay(ch[2*k+1]);
    splay(ch[2*k+2]);
    sh[2*k+1] = apart(ch[2*k+1], 1);
    sh[2*k+2] = apart(ch[2*k+2], 0);
    ch[k] = join(ch[2*k+1], ch[2*k+2]);
  }
  int index(point p) const { return lower_bound(all(ps), p) - ps.begin(); }
  void insert(point p) { build<1>(0, n, 0, index(p)); }
  void erase(point p)  { build<0>(0, n, 0, index(p)); }

  node *find(point p) {
    node *t = ch[0];
    while (t && t->p != p) 
      t = t->child[t->p < p];
    return t;
  }
};
struct dynamic_convex_hull {
  half_hull<+1> upper;
  half_hull<-1> lower;

  dynamic_convex_hull(vector<point> ps) : upper(ps), lower(ps) { }
  void insert(point p) { upper.insert(p); lower.insert(p); }
  void erase(point p) { upper.erase(p); lower.erase(p); }

  real area() {
    real area = 0;
    if (upper.ch[0]) area += upper.ch[0]->area; 
    if (lower.ch[0]) area -= lower.ch[0]->area; 
    return area;
  }
  double perimeter() {
    double perimeter = 0;
    if (upper.ch[0]) perimeter += upper.ch[0]->perimeter; 
    if (lower.ch[0]) perimeter += lower.ch[0]->perimeter; 
    return perimeter;
  }

  struct hull_iterator {
    half_hull<+1>::node *upper_ptr;
    half_hull<-1>::node *lower_ptr;
    hull_iterator(half_hull<+1>::node *upper_ptr, half_hull<-1>::node *lower_ptr) : 
		  upper_ptr(upper_ptr), lower_ptr(lower_ptr) { }

    bool operator==(hull_iterator it) {
			auto is_end = [&](hull_iterator it) { return !it.upper_ptr && !it.lower_ptr; };
			if (is_end(*this) && is_end(it)) return true;
			if (is_end(*this) || is_end(it)) return false;
			return *(*this) == *it;
    }
    bool operator!=(hull_iterator it) { return !operator==(it); }
    void operator++() { 
      if (upper_ptr) upper_ptr = upper_ptr->next[1];
      if (!upper_ptr) {
        lower_ptr = lower_ptr->next[0];
        if (lower_ptr && !lower_ptr->next[0]) lower_ptr = 0;
      }
    }
    point operator*() { 
			return upper_ptr ? upper_ptr->p : lower_ptr->p;
    }
  };
  hull_iterator begin() {
    if (!upper.ch[0]) return hull_iterator(0, 0);
    auto upper_ptr = upper.ch[0];
    for (; upper_ptr->child[0]; upper_ptr = upper_ptr->child[0]);
    auto lower_ptr = lower.ch[0];
    for (; lower_ptr->child[1]; lower_ptr = lower_ptr->child[1]);
    return hull_iterator(upper_ptr, lower_ptr);
  }
  hull_iterator end() { return hull_iterator(0, 0); }
  hull_iterator find(point p) { 
		hull_iterator it = begin();
		it.upper_ptr = upper.find(p);
    if (!it.upper_ptr) it.lower_ptr = lower.find(p);
		return it;
  }
	hull_iterator next(hull_iterator it) {
		++it;
		if (it == end()) it = begin();
		return it;
	}
};

int ACMICPC_JP_2017D() {
  int n;
  scanf("%d", &n);
  vector<point> ps;
  for (int i = 0; i < n; ++i) {
    double x, y;
    scanf("%lf %lf", &x, &y);
    ps.push_back({x, y});
  }
  dynamic_convex_hull hull(ps);
  for (int i = 0; i < n; ++i) 
		hull.insert(ps[i]);

  double per = hull.perimeter();
  double best1 = 0;
  { // best non-adjacent
    vector<pair<double, point>> top;
		for (auto it = hull.begin(); it != hull.end(); ++it) {
      hull.erase(*it);
      top.push_back({per - hull.perimeter(), *it});
      hull.insert(*it);
    }
    sort(all(top)); reverse(all(top));
    for (int i = 0; i < min(4, (int)top.size()); ++i) {
      point p = top[i].snd;
      for (int j = i+1; j < min(4, (int)top.size()); ++j) {
        point q = top[j].snd;
        auto it = hull.next(hull.find(p));
        if (*it == q) continue;
        it = hull.next(hull.find(q));
        if (*it == p) continue;
        best1 = max(best1, top[i].fst + top[j].fst);
      }
    }
  }
	double best2 = 0;
  { // best consecurive two
		for (auto it = hull.begin(); it != hull.end(); ++it) {
			auto jt = hull.next(it);
      hull.erase(*it);
      hull.erase(*jt);
      best2 = max(best2, per - hull.perimeter());
      hull.insert(*it);
      hull.insert(*jt);
    }
  }
  double best3 = 0;
  { // best successive two
		for (auto it = hull.begin(); it != hull.end(); ++it) {
			auto jt = hull.next(it), kt = hull.next(jt);
      hull.erase(*jt);
			for (auto lt = it; lt != kt; lt = hull.next(lt)) {
        hull.erase(*lt);
        best3 = max(best3, per - hull.perimeter());
        hull.insert(*lt);
			}
      hull.insert(*jt);
		}
  }
  printf("%.5f\n", max({best1, best2, best3}));
}

//-----------------------------------------------------------------------------
// kd-tree with distance biseparator.
//
// if we construct naively, it tooks O(n^2).
// using candle burning method (by Bentley), it can construct in O(n log n).
// The following implementation is lazy so O(n log^2 n).
//-----------------------------------------------------------------------------
struct split_tree {
  struct node {
    node *l, *r;
    int d;
    point h;
    int id;
  } *root;
  int n;
  vector<point> ps;
  pair<point,int> ref(pair<point,int> p) {
    return {{p.fst.y, p.fst.x}, p.snd};
  }
  split_tree(vector<point> ps) : n(ps.size()), ps(ps) {
    set<pair<point,int>> s[2];
    for (int i = 0; i < n; ++i) {
      pair<point,int> z = {ps[i], i};
      s[0].insert(z);
      s[1].insert(ref(z));
    }
    function<node*(set<pair<point,int>>*)> rec = [&](set<pair<point,int>> s[]) {
      if (s[0].empty()) return (node*)0;
      if (s[0].size() == 1) {
        int id = s[0].begin()->snd;
        return new node({0, 0, -1, ps[id], id});
      }
      int d = 0;
      if (prev(s[0].end())->fst - s[0].begin()->fst <
          prev(s[1].end())->fst - s[1].begin()->fst) d = 1;
      auto it = s[d].begin(), jt = prev(s[d].end());
      point h = (jt->fst + it->fst) / 2; // left if smaller than h
      while (1) {
        if (it == jt || h <= it->fst) break; ++it;
        if (it == jt || jt->fst < h) break; --jt;
      }
      node *l, *r;
      set<pair<point,int>> ss[2];
      if (it->fst < h) { // violate jt side
        ++jt;
        do {
          ss[!d].insert(ref(*jt));
          s[!d].erase(ref(*jt));
          ss[d].insert(*jt);
          jt = s[d].erase(jt);
        } while (jt != s[d].end());
        r = rec(ss);
        l = rec(s);
      } else { // violate it side
        do {
          --it;
          ss[!d].insert(ref(*it));
          s[!d].erase(ref(*it));
          ss[d].insert(*it);
          it = s[d].erase(it);
        } while (it != s[d].begin());
        l = rec(ss);
        r = rec(s);
      }
      return new node({l, r, d, h, -1});
    };
    root = rec(s);
  }

  int depth(node *u) {
    if (!u) return 0;
    return 1+max(depth(u->l), depth(u->r));
  }
  void disp(node *u, int tab = 0) {
    if (!u) return;
    if (u->d < 0) {
      cout << string(tab, ' ') << u->h << endl;
      return;
    } else {
      cout << string(tab, ' ') << "split at " << (u->h) << " (" << u->d << ")" << endl;
      disp(u->l, tab+2);
      disp(u->r, tab+2);
    }
  }
};

// naive. O(n^2) in worst case. or O(n log(largest gap/smallest gap))
struct split_tree_n {
  struct node {
    node *l, *r;
    int d;
    point h;
    int id;
  } *root;

  int n;
  vector<point> ps;
  pair<point,int> ref(pair<point,int> p) {
    return {{p.fst.y, p.fst.x}, p.snd};
  }
  split_tree_n(vector<point> ps) : n(ps.size()), ps(ps) {
    function<node*(vector<int>)> rec = [&](vector<int> id) {
      if (id.size() == 0) return (node*)0;
      if (id.size() == 1) {
        return new node({0, 0, -1, ps[id[0]], id[0]});
      }
      vector<double> xs, ys;
      for (int i: id) {
        xs.push_back(ps[i].x);
        ys.push_back(ps[i].y);
      }
      sort(all(xs));
      sort(all(ys));
      if (xs.back() - xs[0] > ys.back() - ys[0]) {
        real d = (xs.back()+xs[0])/2;
        vector<int> ls, rs;
        for (int i: id) {
          if (ps[i].x < d) ls.push_back(i);
          else             rs.push_back(i);
        }
        return new node({rec(ls), rec(rs), 0, {d,0}, -1});
      } else {
        real d = (ys.back()+ys[0])/2;
        vector<int> ls, rs;
        for (int i: id) {
          if (ps[i].y < d) ls.push_back(i);
          else             rs.push_back(i);
        }
        return new node({rec(ls), rec(rs), 1, {0,d}, -1});
      }
    };
    vector<int> id(n); iota(all(id), 0);
    root = rec(id);
  }

  int depth(node *u) {
    if (!u) return 0;
    return 1+max(depth(u->l), depth(u->r));
  }
  void disp(node *u, int tab = 0) {
    if (!u) return;
    if (u->d < 0) {
      cout << string(tab, ' ') << u->h << endl;
      return;
    } else {
      cout << string(tab, ' ') << "split at " << (u->h) << " (" << u->d << ")" << endl;
      disp(u->l, tab+2);
      disp(u->r, tab+2);
    }
  }
};

// TODO
// https://graphics.stanford.edu/courses/cs468-06-fall/Papers/02%20vladlen%20notes.pdf
// epsilon net
// Let (X, F) be a set system. N \subset X is called epsilon net if
// N cat S \neq empty for all S in F with u(S) > epsilon.
// (if S in F has area greater than epsilon, S intersect with N).
//
// Theorem. there exists e-net of size O(d (1/e) log(1/e))
//



// vor(i) is the polygon such that for all p in vor(i)
// ps[i] is the farthest point from p.
//
// clearly: vor(i) is non-empty if ps[i] is on the convex hull.
//
// TODO
struct farthest_voronoi {
  int n;
  vector<point> ps, qs;
  map<point,int> id;
  /*
  int node(point p) {
    if (!id.count(p)) { id[p] = m++; qs.push_back(p); adj.push_back({}); }
    return id[p];
  }
  */
  vector<vector<point>> vor;

  farthest_voronoi(vector<point> ps, vector<point> domain) :
    n(ps.size()), ps(ps) {
    vector<int> id(n); iota(all(id), 0);
    sort(all(id), [&](int i, int j) { return ps[i] < ps[j]; });

    vector<int> ch(2*n);
    auto cond = [&](int i, int j, int k) {
      return sign(cross(ps[i]-ps[k], ps[j]-ps[k])) < 0;
    };
    int k = 0;
    for (int i = 0; i < n; ch[k++] = i++)
      for (; k >= 2 && cond(ch[k-2], ch[k-1], i); --k);
    for (int i = n-2, t = k+1; i >= 0; ch[k++] = i--)
      for (; k >= t && cond(ch[k-2], ch[k-1], i); --k);
    ch.resize(k-1);

    cout << ch << endl;

    vector<int> prev(n, -1), next(n, -1);
    for (int i = 0; i < ch.size(); ++i) {
      int j = i + 1; if (j >= ch.size()) j = 0;
      next[ch[i]] = ch[j];
      prev[ch[j]] = ch[i];
    }

    auto disp = [&](int s) {
      int u = s;
      do {
        cout << u << " ";
        u = next[u];
      } while (u != s);
      cout << endl;
    };
    shuffle(all(ch), mt19937());
    vor.resize(n);
    for (int i = (int)ch.size()-1; i >= 1; --i) {
      disp(ch[0]);
      next[prev[ch[i]]] = next[ch[i]];
      prev[next[ch[i]]] = prev[ch[i]];
    }
    // つらい
    // 基本 DCEL がらみはしんどい
    vor[ch[0]] = domain;
    cout << "---" << endl;
    for (int i = 1; i < ch.size(); ++i) {
      next[prev[ch[i]]] = ch[i];
      prev[next[ch[i]]] = ch[i];
      disp(ch[0]);
      // ch[i] と next[ch[i]] の bisector で切る
      // 領域を hyperplane で管理するのが良いのか
      /*
      for (int j = 0; j < i; ++j) {
        point s = (ps[e.src]+ps[e.dst])/2, d = orth(ps[e.dst]-ps[e.src]);
        region = convex_cut(region, {s, s+d});
      }
      */
    }

    // ランダムに挿入しながら構築する．
    // 点 i を挿入するとき．next[i] の voronoi region を i, next[i] の bisector で切る
    // bisector があたる領域を j の voronoi region として i, j の bisector で切る
    // これを繰り返す
    //
    // データ構造処理がデス辛い
    //
    // 必要な演算
    // voronoi region の辺集合

  }

};

/*
struct doubly_connected_edge_list {
  struct Vertex;
  struct Edge;
  struct Face;
  struct vertex {
    Edge *edge; // any incident edge
    point p; // 
  };
  struct Edge {
    Vertex *vertex; // origin
    Edge *prev, *next; // surrounding edges of face
    Edge *twin; // reverse edge of e
    Face *face; // left face
  };
  struct Face {
    Edge *edge; // Any incident edge
  };

  // process edges incident to v in CCW order
  void incident(Vertex *v) {
    Edge *e = vertex->edge;
    do {
      // process e
      e = e->prev->twin;
    } while (e != v->edge)
    v->edge;
  }
  // process vertices incident to e in SRC-DST order
  void incident(Edge *e) {
    // process e->vertex;
    // process e->twin->vertex;
  }
  // process edges incident to f in CCW order
  void incident(Face *f) {
    Edge *e = f->edge;
    do {
      // process e
      e = e->next;
    } while (e != f->edge);
  }
};
*/

// 
// Polyhedral Region Overlay
//
struct DCEL {
  struct Vertex;
  struct Edge;
  struct Vertex {
    point p;
    Edge *edge;
  };
  struct Edge {
    Vertex *vertex;
    Edge *twin;
    Edge *prev, *next;
    // Face *face;
  };

  DCEL copy() {

  }
};
void overlay(DCEL *a, DCEL *b) {
  // 
  // a と b で segment intersection をとく
  // 各 intersection について新しい点をおいて edge を分割
  // その後，face を assign する
  // 
};


// Convex Hull of circles; めっちゃバグってる
//
// divide and conquer？？
//
// ２つの凸包から，それらの union の凸包を作ればいい
//
// 最も右端にある円を選び，真上への接直線を引く
// dominant なほうを追加する
//
// TODO
//  http://ac.els-cdn.com/092577219290015K/1-s2.0-092577219290015K-main.pdf?_tid=f2321600-62e3-11e6-aa7a-00000aacb35e&acdnat=1471264364_f7d8c0796b081532ec0e84e22028f111
//
// [non-verified; in progress]
vector<circle> convex_hull(vector<circle> cs) {
  int n = cs.size();
  typedef vector<int> hull;
  function<hull(hull,hull)> merge = [&](hull P, hull Q) {
    cout << "---" << endl;
    cout << "merge ";
    for (auto i: P) cout << "[" << cs[i].p << "," << cs[i].r << "] ";
    cout << "and";
    for (auto j: Q) cout << " [" << cs[j].p << "," << cs[j].r << "]";
    cout << endl;
    hull S;
    unordered_set<int> elem;
    auto add = [&](int k) {
      if (elem.count(k)) return;
      cout << "add [" << cs[k].p << "," << cs[k].r << "]" << endl;
      S.push_back(k);
      elem.insert(k);
    };
    auto supp = [&](int k, point v) {
      circle c = cs[k];
      return c.p - c.r * orth(v);
    };
    auto dom = [&](int i, int j, point v) {
      point p = supp(P[i], v), q = supp(Q[j], v);
      cout << "p = " << p << ", q = " << q << endl;
      int s = sign(cross(p - q, v));
      cout << "s = " << s << endl;
      return s ? s > 0 : sign(dot(q - p, v)) > 0;
    };
    auto tangent = [&](int k, int l) {
      if (k == l) return point({0,0});
      auto u = cs[l].p - cs[k].p;
      auto b = cs[k].r - cs[l].r, g = norm(u);
      u /= g;
      auto h = b / g;
      cout << h << endl;
      if (sign(1 - h*h) < 0) return point({0,0});
      return orth(u*h - orth(u)*sqrt(max(0.0l, 1 - h*h)));
    };
    auto compare = [&](point a, point b, point v) {
      cout << "compare " << a << " " << b << " " << v << endl;
      if (sign(dot(a,a)) == 0) return false;
      if (sign(dot(b,b)) == 0) return true;
      // TODO: v と a, b は同じ方向ではいけない（大丈夫？）
      // a はゼロかもしれない
      int s = sign(cross(v, a)), t = sign(cross(v, b));
      cout << "s = " << s << ", t = " << t << endl;
      if (s == 0 && dot(v, a) > 0) return true;
      if (t == 0 && dot(v, b) > 0) return false;
      return s != t ? s > t : sign(cross(a, b)) > 0;
    };
    auto advance = [&](int &i, int &j, point &v) {
      int I = (i+1) % P.size(), J = (j+1) % Q.size();
      point a = tangent(P[i], Q[j]), b = tangent(P[i], P[I]),
            c = tangent(Q[j], Q[J]), d = tangent(Q[j], P[i]);
      cout << a << " " << b << " " << c << " " << d << endl;
      if (compare(a, b, v) && compare(a, c, v)) { cout << a << " is the first turn" << endl; add(Q[j]); }
      if (compare(d, b, v) && compare(d, c, v)) { cout << d << " is the first turn" << endl; add(P[i]); }
      if (compare(b, c, v)) { v = b; i = I; }
      else                  { v = c; j = J; }
    };
    int i = 0, j = 0;
    point v = {1, 0};
    for (int iter = 0; iter < 3+P.size()+Q.size(); ++iter) {
      if (sign(norm(v)) == 0) break;
      cout << "currently " << cs[P[i]].p << " and " << cs[Q[j]].p << endl;
      if (dom(i, j, v)) {
        cout << cs[P[i]].p << " is the lowermost" << endl;
        add(P[i]);
      } else {
        cout << cs[Q[j]].p << " is the lowermost" << endl;
        add(Q[j]);
      }
      advance(i, j, v);
    }
    cout << "==> ";
    for (auto i: S) cout << "[" << cs[i].p << "," << cs[i].r << "] ";
    cout << endl;
    return S;
  };
  function<hull(int,int)> rec = [&](int l, int r) {
    if (l+1 == r) return (hull){l};
    auto P = rec(l, (l+r)/2), Q = rec((l+r)/2, r);
    return merge(P, Q);
  };
  auto ch = rec(0, n);
  vector<circle> res;
  for (int i: ch) res.push_back(cs[i]);
  return res;
}

void verify_convex_hull_discs() {
  int n = 20;
  vector<circle> cs;
  for (int i = 0; i < n; ++i)
    cs.push_back({ {urand(), urand()}, urand()/3 });
  auto ch = convex_hull(cs);

  visualizer vis;
  for (auto c: cs) vis << c;
  point shift = point({vis.maxx+1, 0});
  for (auto c: ch) vis << circle({c.p + shift, c.r});
}

void verify_farthest_voronoi() {

  vector<point> ps = {
    {0,0},
    {2,0},
    {1,1},
    {2,2},
    {0,2},
  };
  vector<point> region = {
    {-10,-10},
    { 10,-10},
    { 10, 10},
    {-10, 10}
  };
  farthest_voronoi V(ps, region);
}


// verify
//
void verify_delaunay() {
  {
  vector<point> ps = {
    //{0,0}, {1,0}, {1,1},
    /*
    {0,0}, {0,1}, {0,2},
    {1,0}, {1,1},
    */
    /*
    {0,0}, {1,0}, {2,0},
    {0,1}, {1,1}, {2,1},
    {0,2}, {1,2}, {2,2},
    */
  };
  }
  int n = 100000;
  vector<point> ps(n);
  for (int i = 0; i < n; ++i)
    ps[i].x = urand(), ps[i].y = urand();
  tick();
  delaunay DT(ps);
  cout << tick() << endl;

  visualizer vis;
  for (point p: ps) vis << p;
  for (int u = 0; u < DT.n; ++u) {
    for (auto e: DT.adj[u])
      vis << segment({ps[e.src], ps[e.dst]});
  }
  /*
  for (int n = 5; n < 7; n *= 1.5) {
    vector<point> ps(n);
    for (int i = 0; i < n; ++i)
      ps[i].x = 100 * urand(), ps[i].y = 100 * urand();
    delaunay DT(ps);
  }
  */
}
void verify_voronoi() {
  /*
  vector<point> ps = {
    {0,0}, {1,0}, {0,1}
  };
  */
  for (int n = 5; n < 100000; n *= 1.5) {
    vector<point> ps(n);
    for (int i = 0; i < n; ++i)
      ps[i].x = urand(), ps[i].y = urand();

    vector<point> region = {
      {-0.1,-0.1},
      { 1.1,-0.1},
      { 1.1, 1.1},
      {-0.1, 1.1}
    };
    cout << "n = " << n << " ";
    voronoi Vor(ps, region);

    visualizer vis;
    for (point p: ps) vis << p;
    for (int u = 0; u < Vor.m; ++u) {
      //vis << Vor.qs[u];
      for (auto e: Vor.adj[u])
        vis << segment({Vor.qs[e.src], Vor.qs[e.dst]});
    }
  }
}

void POJ2235() {
  for (int n; ~scanf("%d", &n); ) {
    vector<point> ps(n);
    for (int i = 0; i < n; ++i)
      scanf("%lf %lf", &ps[i].x, &ps[i].y);
    delaunay DT(ps);
  }
}
void AOJ1514() {
  for (int n, m; scanf("%d %d", &n, &m); ) {
    if (n == 0) break;
    vector<point> ps(n);
    for (int i = 0; i < n; ++i)
      scanf("%lf %lf", &ps[i].x, &ps[i].y);
    delaunay DT(ps);
    vector<double> score(n);
    double ans = 0.0;
    for (int j = 0; j < m; ++j) {
      point c, d;
      double s;
      scanf("%lf %lf %lf %lf %lf", &c.x, &c.y, &d.x, &d.y, &s);
      vector<point> domain = {
        c + point({-d.x, -d.y}),
        c + point({+d.x, -d.y}),
        c + point({+d.x, +d.y}),
        c + point({-d.x, +d.y}),
      };
      voronoi V(DT, domain);
      visualizer vis;
      for (point p: ps) vis << p;
      for (int u = 0; u < V.m; ++u) {
        for (auto e: V.adj[u])
          vis << segment({V.qs[e.src], V.qs[e.dst]});
      }
      vector<point> &qs = V.qs;
      for (int i = 0; i < n; ++i) {
        real area = 0;
        int K = V.cell[i].size();
        for (int k = 0; k < K; ++k)
          area += cross(qs[k], qs[(k+1)%K]);
        score[i] += s * area / 2 / (4 * d.x * d.y);
        ans = max(ans, score[i]);
      }
    }
    printf("%.12lf\n", ans);
  }
}

void verify_rectangle_union() {
  vector<rectangle> rs = {
    { {0,0},{2,2} },
    { {1,1},{4,2} },
    { {3,0},{5,2} },
  };
  cout << rectangle_union(rs) << endl;
}

void verify_relative_neighborhood_graph() {
  vector<point> ps = {
    {0,0},
    {1,0},
    {0,1},
    {1,1}
  };
  relative_neighborhood_graph rng(ps);
  for (int u = 0; u < rng.n; ++u) {
    for (auto e: rng.adj[u]) {
      cout << ps[e.src] << " " << ps[e.dst] << endl;
    }
  }
}


void verify_nearest_neighbor_structure() {
  int n = 1000000, m = 100, k = 10;
  vector<point> ps(n);
  for (int i = 0; i < n; ++i) {
    ps[i].x = urand();
    ps[i].y = urand();
  }
  tick();
  vantage_point_tree VPT(ps);
  cout << "vantage point tree: " << tick() << endl;
  random_ball_cover RBC(ps);
  cout << "random ball cover: " << tick() << endl;

  double t1 = 0, t2 = 0;
  for (int i = 0; i < m; ++i) {
    point p;
    p.x = urand();
    p.y = urand();

    cout << "---" << endl;
    tick();
    cout << VPT.closest(p, k) << endl;
    t1 += tick();
    cout << RBC.closest(p, k) << endl;
    t2 += tick();
  }
  cout << "vantage point tree: " << t1 << endl;
  cout << "random ball over: " << t2 << endl;
}



void verify_points_counter() {
  vector<point> ps = {
    {0,1},
    {1,0},
    {1,1},
    {2,1},
  };
  points_counter T(ps);
}









void verify_intersectCC() {
  auto th = 3.141592653589*rand()/(1.0+RAND_MAX);
  auto dx = 1 - 2*rand()/(1.0+RAND_MAX);
  auto dy = 1 - 2*rand()/(1.0+RAND_MAX);
  auto S = [&](point p) {
    real c = cos(th), s = sin(th);
    p.x -= dx; p.y -= dy;
    return (point){c * p.x + s * p.y, -s * p.x + c*p.y};
  };
  auto T = [&](circle C) {
    return (circle){ S(C.p), C.r };
  };
  {
    circle C = { {0,0}, 1 };
    circle D = { {2,0}, 1 };
    TEST(intersect(T(C), T(D)).size(), 1); // touch outer
  }
  {
    circle C = { {0,0}, 1 };
    circle D = { {1,0}, 1 };
    TEST(intersect(T(C), T(D)).size(), 2); // properly intersect inner
  }
  {
    circle C = { {0,0}, 2 };
    circle D = { {1,0}, 1 };
    TEST(intersect(T(C), T(D)).size(), 1); // intersect inner
  }
  {
    circle C = { {0,0}, 2 };
    circle D = { {3,0}, 2 };
    TEST(intersect(T(C), T(D)).size(), 2); // properly intersect outer
  }
  {
    circle C = { {0,0}, 3 };
    circle D = { {1,0}, 1 };
    TEST(intersect(T(C), T(D)).size(), 0); // too close
  }
  {
    circle C = { {0,0}, 1 };
    circle D = { {3,0}, 1 };
    TEST(intersect(T(C), T(D)).size(), 0); // too far
  }
}

void verify_three_point_circle() {
  for (int iter = 0; iter < 100; ++iter) {
    point p = {urand(), urand()};
    point q = {urand(), urand()};
    point r = {urand(), urand()};
    circle c = three_point_circle(p, q, r);
    TEST(dot(c.p - p, c.p - p), dot(c.p - q, c.p - q));
    TEST(dot(c.p - p, c.p - p), dot(c.p - r, c.p - r));
  }
}

void verify_triangulate() {
  vector<point> ps;
  ps.push_back({0, 0});
  ps.push_back({2, 0});
  ps.push_back({2, 2});
  ps.push_back({1, 1});
  ps.push_back({0, 2});
  TEST(triangulate(ps), area(ps));
}

void verify_convex_cut() {
  int n;
  scanf("%d", &n);
  vector<point> ps(n);
  for (int i = 0; i < n; ++i) {
    scanf("%lf %lf", &ps[i].x, &ps[i].y);
  }
  int m;
  scanf("%d", &m);
  for (int i = 0; i < m; ++i) {
    line L;
    scanf("%lf %lf %lf %lf", &L.p.x, &L.p.y, &L.q.x, &L.q.y);
    vector<point> qs = convex_cut(ps, L);
    printf("%f\n", area(qs));
  }
}
void verify_convex_cut2() {
  vector<point> ps = {
    {0.,0.},
    {1.,0.},
    {2.,0.},
    {2.,1.},
    {2.,2.},
    {1.,2.},
    {0.,2.},
    {0.,1.}
  };
  line L = {{1,0},{1.5,0}};
  vector<point> qs = convex_cut(ps, L);
  for (auto q: qs) cout << q << " "; cout << endl;
}
void verify_tangent() {
  point p;
  circle c;
  scanf("%lf %lf", &p.x, &p.y);
  scanf("%lf %lf %lf", &c.p.x, &c.p.y, &c.r);
  vector<point> ps;
  for (auto L: tangent(c, p)) {
    if (L.p == p) ps.push_back(L.q);
    else          ps.push_back(L.p);
  }
  sort(all(ps));
  for (auto p: ps)
    printf("%.12f %.12f\n", p.x + EPS, p.y + EPS);
}
void verify_tangentCC() {
  TEST(tangent({{-2.0, 0.0}, 2.0}, {{-2.0, 0.0}, 1.0}).size(), 0);
  TEST(tangent({{-2.0, 0.0}, 2.0}, {{-1.0, 0.0}, 1.0}).size(), 1);
  TEST(tangent({{-2.0, 0.0}, 2.0}, {{-0.5, 0.0}, 1.0}).size(), 2);
  TEST(tangent({{-2.0, 0.0}, 2.0}, {{+1.0, 0.0}, 1.0}).size(), 3);
  TEST(tangent({{-2.0, 0.0}, 2.0}, {{+2.0, 0.0}, 1.0}).size(), 4);
  TEST(tangent({{+2.0, 0.0}, 2.0}, {{+2.0, 0.0}, 1.0}).size(), 0);
  TEST(tangent({{+2.0, 0.0}, 2.0}, {{+1.0, 0.0}, 1.0}).size(), 1);
  TEST(tangent({{+2.0, 0.0}, 2.0}, {{+0.5, 0.0}, 1.0}).size(), 2);
  TEST(tangent({{+2.0, 0.0}, 2.0}, {{-1.0, 0.0}, 1.0}).size(), 3);
  TEST(tangent({{+2.0, 0.0}, 2.0}, {{-2.0, 0.0}, 1.0}).size(), 4);

  for (int iter = 0; iter < 100; ++iter) {
    circle c = {{urand(), urand()}, urand()};
    circle d = {{urand(), urand()}, urand()};
    for (line l: tangent(c, d)) {
      TEST(intersect(l, c).size(), 1);
      TEST(intersect(l, d).size(), 1);
    }
  }
}
void verify_tangentCC2() {
  circle c, d;
  cin >> c.p.x >> c.p.y >> c.r;
  cin >> d.p.x >> d.p.y >> d.r;
  vector<line> ls = tangent(c, d);
  vector<point> ps;
  for (line l: ls) {
    ps.push_back(intersect(l, c)[0]);
  }
  sort(all(ps));
  for (point p: ps) {
    double x = p.x; if (sign(x) == 0) x = 0;
    double y = p.y; if (sign(y) == 0) y = 0;
    printf("%.12f %.12f\n", x, y);
  }
}
void verify_intersect_area() {
  vector<point> ps = {
    {0.,0.},
    {3.,0.},
    {0.,3.},
  };
  circle c = {{0.0,0.0},2.8};
  cout << intersection_area(ps, c) << endl;
}

void verify_closest_pair() {
  double t1 = 0, t2 = 0;
  for (int seed = 3; seed < 100; ++seed) {
    srand(seed);
    int n = 20000;
    vector<point> ps;
    for (int i = 0; i < n; ++i) {
      ps.push_back({urand(), urand()});
    }
    tick();
    auto uv = closest_pair(ps);
    t1 += tick();
    auto wx = closest_pair2(ps);
    t2 += tick();
    auto d1 = dist(uv.fst, uv.snd);
    auto d2 = dist(wx.fst, wx.snd);
    cout << uv.fst << " " << uv.snd << " " << dist(uv.fst, uv.snd) << endl;
    cout << wx.fst << " " << wx.snd << " " << dist(wx.fst, wx.snd) << endl;
    cout << endl;
    if (sign(d1 - d2) != 0) {
      cout << seed << endl;
      break;
    }
  }
  cout << t1 << " " << t2 << endl;
}

vector<point> half_plane_intersection_n(vector<line> ls) {
  real s = 99999999;
  vector<point> ps = {{-s,-s},{s,-s},{s,s},{-s,s}};
  for (int i = 0; i < ls.size(); ++i)
    ps = convex_cut(ps, ls[i]);
  return ps;
}
void verify_half_plane_intersection() {
  for (int seed = 0; seed < 1000; ++seed) {
    srand(seed);
    int n = 3 + rand() % 100;
    vector<line> ls;
    for (int i = 0; i < n; ++i) {
      double th1 = 2 * PI * i / n;
      double th2 = 2 * PI * (i+1) / n;
      double r1 = 10 + 100 * urand();
      double r2 = 10 + 100 * urand();
      ls.push_back({{r1*cos(th1), r1*sin(th1)},
                    {r2*cos(th2), r2*sin(th2)}});
    }
    auto a1 = area(half_plane_intersection(ls));
    auto a2 = area(half_plane_intersection_n(ls));
    if (sign(1 - a2/a1)) { cout << a1 << " " << a2 << " " << 1 - a2/a1 << " " << sign(1-a2/a1) << seed << endl; exit(-1); }
  }
}

// TEST
void angular_sort_n(vector<point> &ps) {
  vector<tuple<real,real,point>> qs;
  for (point p: ps) {
    qs.push_back(make_tuple(arg(p), norm2(p), p));
    if (get<0>(qs.back()) < 0) get<0>(qs.back()) += 2 * PI;
  }
  sort(all(qs));
  for (int i = 0; i < ps.size(); ++i)
    ps[i] = get<2>(qs[i]);
}
void verify_polar_angle() {
  double t1 = 0, t2 = 0;
  int n = 1000000;
  vector<point> ps;
  for (int i = 0; i < n; ++i)
    ps.push_back({(double)(rand()%2001)-1000, (double)(rand()%2001)-1000});
  vector<point> qs = ps;
  tick();
  sort(all(ps), polar_angle());
  t1 += tick();
  angular_sort_n(qs);
  t2 += tick();
  for (int i = 0; i < ps.size(); ++i) {
    if (!(ps[i] == qs[i])) { cout << "*"; }
    //cout << ps[i] << " " << qs[i] << endl;
  }
  cout << t1 << " " << t2 << endl;
}

void verify_maximum_circle_cover() {
  int n = 10000;
  vector<point> ps;
  for (int i = 0; i < n; ++i) {
    double x = rand() % 1000, y = rand() % 1000;
    ps.push_back({x, y});
  }
  tick();
  cout << maximum_circle_cover(ps, 20) << endl;
  cout << tick() << endl;
  cout << maximum_circle_cover2(ps, 20) << endl;
  cout << tick() << endl;
}
void verify_maximum_circle_cover2() {
  for (int n; scanf("%d", &n) && n; ) {
    vector<point> ps(n);
    for (int i = 0; i < n; ++i)
      scanf("%lf %lf", &ps[i].x, &ps[i].y);
    printf("%d\n", maximum_circle_cover2(ps, 1.0));
  }
}


void verify_circle_intersection_area() {
  circle c = {{0,0}, 2};
  circle d = {{2.8,0}, 1};
  cout << intersection_area(c, d) << endl;

  // Monte-Carlo integration
  real w = 10.0, vol = w*w; // w times w box centered at (0,0)
  int count = 0, maxcount = 100000000;
  for (int i = 0; i < maxcount; ++i) {
    point p = {w*urand() - w/2, w*urand() - w/2};
    if (contains(c, p) && contains(d, p)) ++count;
  }
  real r = 1.0*count/maxcount;
  real mu = vol * r, sigma = vol*sqrt((r-r*r)/maxcount);
  cout << "95\% confidence: [" << mu-3*sigma << ", " << mu+3*sigma << "]" << endl;
}

void verify_geometric_median() {
  int n = 10;
  vector<point> ps(n);
  for (int i = 0; i < n; ++i) {
    ps[i].x = urand();
    ps[i].y = urand();
  }
  cout << geometric_median(ps) << endl;
}
void verify_geometric_median2() {
  while (1) {
    vector<point> ps(4);
    for (int i = 0; i < 4; ++i)
      cin >> ps[i].x >> ps[i].y;
    if (ps[0].x < 0) break;
    point g = geometric_median(ps);
    //cout << g << endl;
    double w = 0;
    for (int i = 0; i < 4; ++i)
      w += norm(ps[i] - g);
    printf("%.4f\n", w);
  }
}
void verify_tangent_circles() {
  line l = { {0,0},{1,0} };
  line m = { {0,2},{1,2} };
  for (circle c: tangent_circles(l, m, 1.0)) {
    cout << c.p << " " << c.r << endl;
    cout << "   "; for (point p: intersect(c, l)) { cout << p << " "; } cout << endl;
    cout << "   "; for (point p: intersect(c, m)) { cout << p << " "; } cout << endl;
  }
}


void verify_arrangement() {
  vector<segment> ss = {
    {{0,0},{2,0}},
    {{3,0},{5,0}},
    {{1,0},{4,0}}
  };
  arrangement arr(ss);
}

void verify_is_congruence() {
  int n = 10;
  vector<point> ps(n);
  for (int i = 0; i < n; ++i) {
    ps[i] = {urand(), urand()};
  }
  double r = urand();
  point u = {urand(), urand()};
  double c = urand(), s = sqrt(1 - c*c);
  vector<point> qs;
  for (int i = 0; i < ps.size(); ++i) {
    qs.push_back(r * (point){c*ps[i].x+s*ps[i].y, -s*ps[i].x+c*ps[i].y} + u);
  }
  int k = rand() % qs.size();
  rotate(qs.begin(), qs.begin()+k, qs.end());
  cout << is_congruence(ps, qs) << endl;
  cout << is_similar(ps, qs) << endl;
}


void verify_merge_segment() {
  vector<segment> ss = {
    {{0,0},{2,0}},
    {{0,0},{0,2}},
    {{4,0},{7,0}},
    {{5,0},{6,0}},
  };
  merge_segments(ss);
}
void verify_farthest_pair() {
  int n = 30;
  vector<point> ps(n);
  for (int i = 0; i < n; ++i) {
    ps[i].x = rand() % 10;
    ps[i].y = rand() % 10;
  }
  auto ans = farthest_pair(ps);

  real best = -1;
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      best = max(best, norm(ps[i]-ps[j]));
    }
  }

  cout << norm(ans.fst-ans.snd) << " " << best << endl;
}

void ACAC002() {
  int n; cin >> n;
  vector<point> ps(n);
  for (int i = 0; i < n; ++i)
    cin >> ps[i];
  auto ans = farthest_pair(ps);
  printf("%.12f\n", norm(ans.snd-ans.fst));
}

void verify_split_tree() {
  for (int n = 1; n < 1000000; n *= 2) {
    vector<point> ps;
    double x = 1, y = 0;
    for (int i = 0; i < n; ++i) {
      x *= 1.001;
      y = 0; //rand() / (RAND_MAX +1.0);
      ps.push_back({x,y});
    }
    sort(all(ps)); ps.erase(unique(all(ps)), ps.end());
    tick();
    split_tree T(ps);
    cout << T.depth(T.root) << endl;
    cout << n << " " << tick() << endl;
  }
}


int main() {
  //verify_intersectCC();
  //verify_three_point_circle();
  //verify_triangulate();
  //verify_convex_cut2();
  //verify_tangent();
  //verify_tangentCC();
  //verify_tangentCC2();
  //verify_intersect_area();
  //verify_closest_pair();
  //verify_half_plane_intersection();
  //verify_polar_angle();
  //verify_maximum_circle_cover();
  //verify_maximum_circle_cover2();
  //verify_rectangle_union();
  //verify_points_counter();
  //verify_circle_intersection_area();
  //verify_geometric_median();
  //verify_tangent_circles();
  //verify_maximum_points_line();
  //verify_arrangement();
  //verify_merge_segment();
  //AOJ1226();
  //AOJ2448();
  //AOJ1247();
  //verify_relative_neighborhood_graph();
  //verify_nearest_neighbor_structure();
  //verify_delaunay();
  //POJ2235();
  // verify_voronoi();
  //verify_is_congruence();
  //AOJ1514();
  //AOJ0273();
  //verify_farthest_pair();
  //ACAC002();
  //verify_compressed_quad_tree();
  //verify_split_tree();
  //verify_farthest_voronoi();
  //verify_convex_hull_discs();
}

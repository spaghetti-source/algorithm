//
// Random Ball Cover 
//
// Description:
// 
//   Random ball cover is a data structure for metric nearest neighbor 
//   search. It is useful if the only distance function is available, 
//   and/or dimension is moderately large (10 to 1000). 
//
//   It first select O(sqrt{n})-size random sample from the points as 
//   representatives. Then, it assigns other points to the nearest
//   representatives. The search procedure uses the triangle inequality
//   to prune unnecessary searching.
//
//
// Complexity:
//
//   Construction: O(n log n).
//   Search: O(sqrt(n)) in a random distribution.
//
//
// Reference:
//
//   Lawrence Cayton (2012):
//   "Accelerating nearest neighbor search on manycore systems."
//   Parallel & Distributed Processing Symposium (IPDPS), 
//   in IPDPS, pp.402-413.
//
#include <bits/stdc++.h>

using namespace std;

#define fst first
#define snd second
#define all(c) ((c).begin()), ((c).end())
#define TEST(s) if (!(s)) { cout << __LINE__ << " " << #s << endl; exit(-1); }


// === tick a time ===
#include <ctime>
double tick() {
  static clock_t oldtick;
  clock_t newtick = clock();
  double diff = 1.0*(newtick - oldtick) / CLOCKS_PER_SEC;
  oldtick = newtick;
  return diff;
}


using Real = long double;
const Real PI = acos(-1.0);
const Real EPS = 1e-8;
int sign(Real x) { return (x > EPS) - (x < -EPS); }
struct Point {
  Real x, y;
  Point &operator+=(Point p) { x += p.x; y += p.y; return *this; }
  Point &operator*=(Real a)  { x *= a;   y *= a;   return *this; }
  Point operator+() const    { return {+x, +y}; }
  Point operator-() const    { return {-x, -y}; }

  Point &operator-=(Point p) { return *this += -p; }
  Point &operator/=(Real a)  { return *this *= 1/a; }
};
Point operator+(Point p, Point q) { return p += q; }
Point operator-(Point p, Point q) { return p -= q; }
Point operator*(Real a, Point p) { return p *= a; }
Point operator*(Point p, Real a) { return p *= a; }
Point operator/(Point p, Real a) { return p /= a; }

int compare(Point p, Point q) { 
  int s = sign(p.x - q.x);
  return s ? s : sign(p.y - q.y);
}
bool operator==(Point p, Point q) { return compare(p,q)==0; }
bool operator!=(Point p, Point q) { return compare(p,q)!=0; }
bool operator<=(Point p, Point q) { return compare(p,q)<=0; }
bool operator>=(Point p, Point q) { return compare(p,q)>=0; }
bool operator<(Point p, Point q) { return compare(p,q)<0; }
bool operator>(Point p, Point q) { return compare(p,q)>0; }

Real dot(Point p, Point q) { return p.x*q.x+p.y*q.y; }
Real cross(Point p, Point q) { return p.x*q.y-p.y*q.x; } // left turn > 0
Real norm2(Point p) { return dot(p,p); }
Point orth(Point p) { return {-p.y, p.x}; }
Real norm(Point p) { return sqrt(dot(p,p)); }
Real arg(Point p) { return atan2(p.y, p.x); }
Real arg(Point p, Point q){ return atan2(cross(p,q), dot(p,q)); }

istream &operator>>(istream &is, Point &p) { is>>p.x>>p.y;return is; }
ostream &operator<<(ostream &os, const Point &p) { os<<"("<<p.x<<","<<p.y<<")"; return os; }

struct RandomBallCover {
  // must satisfy the triangle inequality
  Real dist(Point x, Point y) {
    return sqrt(dot(x-y,x-y));
  }

  vector<Point> ps;
  vector<vector<pair<Real,int>>> list;
  RandomBallCover(vector<Point> ps) : ps(ps) {
    int n = ps.size(), sqrtn = sqrt(n + 0.5);
    vector<int> idx(n);
    iota(all(idx), 0);
    random_shuffle(all(idx));
    for (int i = 0; i < sqrtn; ++i) 
      list.push_back({make_pair(0.0, idx[i])});
    for (int i = sqrtn; i < n; ++i) {
      Real nearest = 1.0/0.0;
      int id;
      for (int j = 0; j < sqrtn; ++j) {
        Real d = dist(ps[list[j][0].snd], ps[idx[i]]);
        if (nearest > d) {
          nearest = d;
          id = j;
        }
      }
      list[id].push_back({nearest, idx[i]});
    }
    for (int i = 0; i < list.size(); ++i)
      sort(all(list[i]));
  }
  int nearestNeighbor(Point p) { 
    vector<Real> dis(list.size()), rem(list.size());
    vector<int> cand(list.size());
    Real best = 1.0/0.0;
    int id;
    for (int i = 0; i < list.size(); ++i) {
      dis[i] = dist(p, ps[list[i][0].snd]);
      if (dis[i] < best) {
        best = dis[i];
        id = list[i][0].snd;
      }
      cand[i] = i;
    }
    sort(all(cand), [&](int i, int j) { return dis[i] < dis[j]; });
    for (int i: cand) {
      for (int k = list[i].size()-1; k >= 0; --k) {
        if (best <= dis[i] - list[i][k].fst) break;
        int j = list[i][k].snd;
        Real d = dist(p, ps[j]);
        if (d < best) {
          best = d;
          id = j;
        }
      }
    }
    return id;
  }
  int nearestNeighborNaive(Point p) {
    int id = 0;
    for (int i = 0; i < ps.size(); ++i) 
      if (dist(p, ps[id]) > dist(p, ps[i])) id = i;
    return id;
  }
};


// for comparison
// Kd-Tree for k = 2
struct KdTree {
  struct Node {
    int id;
    Node *left, *right;
  } *root = 0;
  vector<Point> ps;
  KdTree(vector<Point> ps) : ps(ps) {
    int idx[ps.size()];
    iota(idx, idx+ps.size(), 0);
    root = build<0>(idx, idx+ps.size());
  }
  template <int d>
  Node *build(int *l, int *r) {
    if (l - r >= 0) return 0;
    auto comp = [&](int i, int j) {
      if (d == 0) return ps[i].x < ps[j].x;
      if (d == 1) return ps[i].y < ps[j].y;
    };
    int *m = l + (r-l)/2;
    nth_element(l, m, r, comp);
    return new Node({*m, build<!d>(l,m), build<!d>(m+1,r)});
  }
  template <int d>
  void nearestNeighborSearchRec(Node *x, Point p, pair<Real,int> &ub) {
    if (!x) return;
    Point q = p - ps[x->id];
    Real r = norm(q), w;
    if (r < ub.fst) ub = {r, x->id};
    if (d == 0) w = q.x;
    else        w = q.y;
    Node *fst = x->left, *snd = x->right;
    if (w > 0) swap(fst, snd);
    nearestNeighborSearchRec<!d>(fst, p, ub);
    if (ub.fst > abs(w)) nearestNeighborSearchRec<!d>(snd, p, ub);
  }
  int nearestNeighbor(Point p) {
    pair<Real,int> ub(1.0/0.0, -1);
    nearestNeighborSearchRec<0>(root, p, ub);
    return ub.snd;
  }
};

int main() {
  int n = 10000;
  vector<Point> ps;
  for (int i = 0; i < n; ++i) {
    Real x = rand() % n;
    Real y = rand() % n;
    ps.push_back({x, y});
  }
  RandomBallCover X(ps);
  KdTree T(ps);
  double t1 = 0, t2 = 0;
  for (int i = 0; i < n; ++i) {
    Real x = rand() % n;
    Real y = rand() % n;
    Point p = {x, y};
    tick();
    int j = X.nearestNeighbor(p);
    t1 += tick();
    int k = T.nearestNeighbor(p);
    t2 += tick();
    assert(sign(norm(p - ps[j]) - norm(p - ps[k])) == 0);
  }
  cout << t1 << " " << t2 << endl;
}

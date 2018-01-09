//
// Kd-Tree 
//
// Description:
//
//   Kd-Tree is a spatial partitioning tree defined as follows[1].
//   Each node has a point, called a pivot, and its left
//   and right childs store the points that are lefter and righter
//   on the pivot. Here, "left" and "right" are defined by the 
//   coordinate that switches cyclically over the depth.
//
//   The Kd-Tree has the depth O(log n), and if the input is random,
//   the searching has complexity O(log n) in expectation.
//   However, in an adversarial input, the worst case complexity 
//   becomes O(n^{1-1/d}), which tends to a naive bound as d to infty[2].
//
//   Note: Computing NNs for n points requires O(n^{2-1/d}) time.
//   Basically, it is practical and has advantage over the naive search 
//   if n = 10^5.
//
// Complexity:
//
//   Construction: O(n log n)
//   Nearest Neighbor: O(n^{1-1/d}) in worst case, 
//                     O(log n) in random case.
//
// References:
//
//   [1] Jon Louis Bentley (1975):
//   "Multidimensional binary search trees used for associative searching."
//   Communications of the ACM, vol. 18, no. 9, pp. 509--517.
//
//   [2] Der-Tsai Lee and C. K. Wong (1977):
//   "Worst-case analysis for region and partial region searches in 
//    multidimensional binary search trees and balanced quad trees."
//   Acta Informatica, vol.9, no.1, pp. 23--29.
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


  // for verification
  int nearestNeighborNaive(Point p) {
    int id = 0;
    for (int i = 1; i < ps.size(); ++i)
      if (norm(ps[i] - p) < norm(ps[id] - p)) id = i;
    return id;
  }
  void display(Node *x, int tab=0) {
    if (!x) return;
    display(x->left, tab+2);
    cout << string(tab, ' ') << x->id << ": " << ps[x->id] << endl;
    display(x->right,tab+2);
  }
};
int main() {
  int n = 10000;
  vector<Point> ps(n);
  for (int i = 0; i < n; ++i) {
    ps[i].x = rand() % n;
    ps[i].y = rand() % n;
  }
  KdTree T(ps);
  //T.display(T.root);

  double t1 = 0, t2 = 0;
  for (int k = 0; k < n; ++k) {
    Point p = {rand() % n, rand() % n};
    tick();
    Point q = ps[T.nearestNeighbor(p)];
    t1 += tick();
    Point r = ps[T.nearestNeighborNaive(p)];
    t2 += tick();
    TEST(sign(norm(q-p) - norm(r-p)) == 0);
    //cout << p << " " << norm(q - p) << " " << norm(r - p) << endl;
  }
  cout << t1 << " [s] by Kd search" << endl << 
          t2 << " [s] by naive search" << endl;
}

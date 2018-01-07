// 
// Segment Arrangement (Bentley-Ottman's Plane-Sweep)
//
// Description:
//   Given a set of segments, it finds the all intersections of the
//   segments and construct the graph structure. By the plane-sweep
//   with the balanced binary search tree, it runs in O(k log n)
//   time, where k is the size of output.
//
// Complexity:
//   O(k log n), where k is the size of output.
//
// Verified:
//   AOJ 1226
//   AOJ 2448
//
// References:
//   CGAA
//

#include <bits/stdc++.h>
using namespace std;

#define fst first
#define snd second
#define all(c) ((c).begin()), ((c).end())

using Real = double;
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

struct Segment { Point p, q; };

vector<Point> intersect(Segment s, Segment t) {
  auto a = cross(s.q - s.p, t.q - t.p);
  auto b = cross(t.p - s.p, t.q - t.p);
  auto c = cross(s.q - s.p, s.p - t.p);
  if (a < 0) { a = -a; b = -b; c = -c; }
  if (sign(b) < 0 || sign(a-b) < 0 ||
      sign(c) < 0 || sign(a-c) < 0) return {};      // disjoint
  if (sign(a) != 0) return {s.p + b/a*(s.q - s.p)}; // properly crossing
  vector<Point> ps;                                 // same line
  auto insert_if_possible = [&](Point p) {
    for (auto q: ps) if (sign(dot(p-q, p-q)) == 0) return;
    ps.push_back(p);
  };
  if (sign(dot(s.p-t.p, s.q-t.p)) <= 0) insert_if_possible(t.p);
  if (sign(dot(s.p-t.q, s.q-t.q)) <= 0) insert_if_possible(t.q);
  if (sign(dot(t.p-s.p, t.q-s.p)) <= 0) insert_if_possible(s.p);
  if (sign(dot(t.p-s.q, t.q-s.q)) <= 0) insert_if_possible(s.q);
  return ps;
}


// 
// each vertex has one incident edge
// each edge has the origin vertex, twin edge, and the left face.
// two edges (prev, next) representing a list of edges surrounding a face. 
//
struct DoublyConnectedEdgeList {
  vector<int> incident_edge; // vertex
  vector<int> origin, twin, prev, next, incident_face; // edge
  vector<int> component; // face
  int edges()    const { return origin.size(); }
  int vertices() const { return incident_edge.size(); }
  int faces()    const { return component.size(); }

  vector<Point> point;
  int newVertex(Point p, int e = -1) {
    point.push_back(p);
    incident_edge.push_back(e);
    return vertices()-1;
  }
  int newEdge(int o = -1) {
    origin.push_back(o);
    twin.push_back(-1);
    prev.push_back(-1);
    next.push_back(-1);
    incident_face.push_back(-1);
    return edges()-1;
  }
  int newFace(int e = -1) {
    component.push_back(e);
    return component.size()-1;
  }
  void completeFaces() {
    component.clear();
    fill(all(incident_face), -1);
    for (int e = 0; e < edges(); ++e) {
      if (incident_face[e] >= 0) continue;
      int f = newFace(e), x = e;
      do {
        incident_face[x] = f;
        x = next[x];
      } while (x != e);
    }
  }
};

struct Arrangement : DoublyConnectedEdgeList {
  vector<Segment> segs;

  unordered_map<int, unordered_map<int, int>> adj; // (Vertex, Vertex) -> Edge

  struct Node { // Sweep-Line Structure (RBST)
    int index, size = 1;
    Node *left = 0, *right = 0;
  } *root = 0;
  vector<Node> ns;
  Node *update(Node *x) {
    if (x) {
      x->size = 1;
      if (x->left)  x->size += x->left->size;
      if (x->right) x->size += x->right->size;
    }
    return x;
  }
  Node *merge(Node *x, Node *y) {
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
  tuple<Node*, Node*, Node*> split(Node *x, C cond) {
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
  Node *leftmost(Node *x) { while (x && x->left) x = x->left; return x; }
  Node *rightmost(Node *x) { while (x && x->right) x = x->right; return x; }
  template <class F>
  void process(Node *x, F func) {
    if (!x) return;
    process(x->left, func);
    func(x);
    process(x->right, func);
  }
  Arrangement(vector<Segment> segs_) : segs(segs_) {
    ns.resize(segs.size());
    set<Point> events;
    map<Point, set<int>> L, R;

    for (int i = 0; i < segs.size(); ++i) {
      if (segs[i].q < segs[i].p) swap(segs[i].p, segs[i].q);
      events.insert(segs[i].p);
      events.insert(segs[i].q);
      L[segs[i].p].insert(i);
      R[segs[i].q].insert(i);
      ns[i].index = i;
    }
    vector<int> last(segs.size(), -1);

    while (!events.empty()) {
      const Point p = *events.begin();
      events.erase(events.begin());
      int u = newVertex(p);
      
      auto cond = [&](Node *x) {
        const Segment &s = segs[x->index];
        if (sign(s.q.x - s.p.x) == 0) {
          if (sign(p.y - s.p.y) < 0) return -1;
          if (sign(s.q.y - p.y) < 0) return +1;
          return 0;
        }
        return -sign(cross(s.p - p, s.q - p));
      };
      auto z = split(root, cond); 
      vector<Node*> inserter;
      process(get<1>(z), [&](Node *x) {
        int v = last[x->index];
        if (!adj[u].count(v)) {
          int e = newEdge(u), f = newEdge(v);
          adj[u][v] = e;
          adj[v][u] = f;
          twin[e] = f;
          twin[f] = e;
        }
        if (!R[p].count(x->index)) 
          inserter.push_back(x);
      });
      for (int i: L[p]) 
        if (!R[p].count(i))
          inserter.push_back(&ns[i]);

      sort(all(inserter), [&](Node *x, Node *y) {
        const Segment &s = segs[x->index], &t = segs[y->index];
        return sign(cross(s.q - s.p, t.q - t.p)) >= 0;
      });
      auto addEvent = [&](Node *x, Node *y) {
        if (!x || !y) return;
        vector<Point> ps = intersect(segs[x->index], segs[y->index]);
        for (Point q: ps) 
          if (p < q) events.insert(q);
      };
      if (inserter.empty()) {
        addEvent(rightmost(get<0>(z)), leftmost(get<2>(z)));
      } else {
        addEvent(rightmost(get<0>(z)), inserter[0]);
        addEvent(leftmost(get<2>(z)), inserter.back());
      }
      root = 0;
      for (int i = 0; i < inserter.size(); ++i) {
        last[inserter[i]->index] = u;
        inserter[i]->left = inserter[i]->right = 0;
        root = merge(root, update(inserter[i]));
      }
      root = merge(merge(get<0>(z), root), get<2>(z));
    }
    vector<int> next_inc(next.size()), prev_inc(prev.size());
    for (auto &pp: adj) {
      int u = pp.fst;
      vector<int> es;
      for (auto z: pp.snd) es.push_back(z.snd);
      sort(all(es), [&](int e, int f) { // polar sort
        auto quad = [](Point p) {
          for (int i = 1; i <= 4; ++i, swap(p.x = -p.x, p.y))
            if (p.x > 0 && p.y >= 0) return i;
          return 0;
        };
        const Point p = point[origin[twin[e]]] - point[origin[e]];
        const Point q = point[origin[twin[f]]] - point[origin[f]];
        if (quad(p) != quad(q)) return quad(p) < quad(q);
        return sign(cross(p, q)) > 0;
      });
      incident_edge[u] = es[0];
      for (int i = 0; i < es.size(); ++i) {
        int j = (i + 1) % es.size();
        next_inc[es[i]] = es[j];
        prev_inc[es[j]] = es[i];
      }
    }
    for (int e = 0; e < edges(); ++e) {
      next[e] = prev_inc[twin[e]];
      prev[e] = twin[next_inc[e]];
    }
  }
};
void AOJ1226() {
  for (int n; ~scanf("%d", &n) && n; ) {
    vector<vector<double>> a(4, vector<double>(n));
    for (int k = 0; k < 4; ++k)
      for (int i = 0; i < n; ++i)
        scanf("%lf", &a[k][i]);

    vector<Segment> ss = {
      {{0,0},{0,1}},
      {{0,1},{1,1}},
      {{1,1},{1,0}},
      {{1,0},{0,0}}
    };
    for (int i = 0; i < n; ++i) {
      ss.push_back({{a[0][i],0},{a[1][i],1}});
      ss.push_back({{0,a[2][i]},{1,a[3][i]}});
    }
    Arrangement arr(ss);
    arr.completeFaces();

    double result = 0;
    for (int f = 0; f < arr.faces(); ++f) {
      double area = 0;
      int e = arr.component[f];
      do {
        area += cross(arr.point[arr.origin[e]],
                      arr.point[arr.origin[arr.twin[e]]]),
        e = arr.next[e];
      } while (e != arr.component[f]);
      result = max(result, area);
    }
    printf("%.6f\n", result/2);
  }
}
void AOJ2448() {
  int n; scanf("%d", &n);
  vector<Point> ps(n);
  for (int i = 0; i < n; ++i)
    scanf("%lf %lf", &ps[i].x, &ps[i].y);
  vector<Segment> ss;
  for (int i = 0; i+1 < n; ++i)
    ss.push_back({ps[i], ps[i+1]});
  Arrangement arr(ss);
  arr.completeFaces();

  double result = 0;
  for (int f = 0; f < arr.faces(); ++f) {
    double area = 0;
    int e = arr.component[f];
    do {
      area += cross(arr.point[arr.origin[e]],
                    arr.point[arr.origin[arr.twin[e]]]),
      e = arr.next[e];
    } while (e != arr.component[f]);
    if (area > 0) result += area;
  }
  printf("%.12f\n", result/2);
}

int main() {
  AOJ2448();
  //AOJ1226();
}

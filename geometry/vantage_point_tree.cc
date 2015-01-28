//
// Vantage Point Tree (vp tree)
//
// Description
//   Vantage point tree is a metric tree.
//   Each tree node has a point, radius, and two childs.
//   The points of left descendants are contained in the ball B(p,r)
//   and the points of right descendants are exluded from the ball.
//
//   We can find k-nearest neighbors of a given point p efficiently
//   by pruning search.
//
//   The data structure is independently proposed by J. Uhlmann and 
//   P. N. Yianilos.
//
// Complexity:
//   Construction: O(n log n).
//   Search: O(log n)
//
//   In my implementation, its construction is few times slower than kd tree
//   and its search is bit faster than kd tree.
//
// References
//   J. Uhlmann (1991): 
//   Satisfying General Proximity/Similarity Queries with Metric Trees.
//   Information Processing Letters, vol. 40, no. 4, pp. 175--179.
//
//   Peter N. Yianilos (1993): 
//   Data structures and algorithms for nearest neighbor search in general metric spaces.
//   in Proceedings of the 4th Annual ACM-SIAM Symposium on Discrete algorithms,
//   Society for Industrial and Applied Mathematics Philadelphia, PA, USA. pp. 311--321. 
//
#include <iostream>
#include <vector>
#include <queue>
#include <complex>
#include <algorithm>

using namespace std;

#include <ctime>
double tick() {
  static clock_t oldtick;
  clock_t newtick = clock();
  double diff = 1.0*(newtick - oldtick) / CLOCKS_PER_SEC;
  oldtick = newtick;
  return diff;
}

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
struct vantage_point_tree {
  struct node {
    point p;
    double th;
    node *l, *r; 
  } *root;
  vector<pair<double, point>> aux;
  vantage_point_tree(vector<point> ps) {
    for (int i = 0; i < ps.size(); ++i)
      aux.push_back({0, ps[i]});
    root = build(0, ps.size());
  }
  node *build(int l, int r) {
    if (l == r) return 0;
    swap(aux[l], aux[l + rand() % (r - l)]);
    point p = aux[l++].snd;
    if (l == r) return new node({p});
    for (int i = l; i < r; ++i) 
      aux[i].fst = norm(p - aux[i].snd);
    int m = (l + r) / 2;
    nth_element(aux.begin()+l, aux.begin()+m, aux.begin()+r);
    return new node({p, sqrt(aux[m].fst), build(l, m), build(m, r)});
  }
  priority_queue<pair<double, node*>> que;
  void k_nn(node *t, point p, int k) {
    if (!t) return;
    double d = abs(p - t->p);
    if (que.size() < k) que.push({d, t});
    else if (que.top().fst > d) {
      que.pop();
      que.push({d, t});
    }
    if (!t->l && !t->r) return;
    if (d < t->th) {
      k_nn(t->l, p, k);
      if (t->th - d <= que.top().fst) k_nn(t->r, p, k);
    } else {
      k_nn(t->r, p, k);
      if (d - t->th <= que.top().fst) k_nn(t->l, p, k);
    }
  }
  vector<point> k_nn(point p, int k) {
    k_nn(root, p, k);
    vector<point> ans;
    for (; !que.empty(); que.pop())
      ans.push_back(que.top().snd->p);
    reverse(all(ans));
    return ans;
  }
};

int main() {
  srand( 0xdeadbeef );

  int n = 100000;
  vector<point> ps;
  for (int i = 0; i < n; ++i) 
    ps.push_back(point(rand()%n, rand()%n));

  tick();
  vantage_point_tree T(ps);
  cout << "construct " << n << " points: "  << tick() << "[s]" << endl;

  // search
  tick();
  for (int i = 0; i < n; ++i) 
    T.k_nn(ps[i], 1);
  cout << "search " << n << " points: "  << tick() << "[s]" << endl;

  // verify
  for (int i = 0; i < 100; ++i) {
    point p(rand(), rand());
    point Tp = T.k_nn(p, 1)[0];

    point Tq = ps[0];
    for (auto q: ps) 
      if (norm(p - Tq) > norm(p - q)) Tq = q;
    
    if (abs(norm(Tp - p) - norm(Tq - p)) > 1e-8) {
      cout << norm(Tp - p) << endl;
      cout << norm(Tq - p) << endl;
      cout << "ERROR" << endl;
      return 0;
    }
  }
  cout << "verification passed" << endl;
}

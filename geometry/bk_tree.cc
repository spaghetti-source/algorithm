//
// Burkhard-Keller Tree (metric tree)
//
// Description:
//   Let V be a (finite) set, and d: V x V -> R be a metric.
//   BK tree supports the following operations:
//     - insert(p): insert a point p, O((log n)^2)
//     - traverse(p,d): enumerate all q with d(p,q) <= d
//
// Remark:
//   To delete elements and/or rebalance the tree,
//   we can use the same technique as the scapegoat tree.
// 
// Reference
//   W. Burkhard and R. Keller (1973):
//   Some approaches to best-match file searching, 
//   Communications of the ACM, vol. 16, issue. 4, pp. 230--236.
//

#include <map>
#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <vector>
#include <cstring>
#include <functional>
#include <unordered_map>
#include <algorithm>

using namespace std;

#define fst first
#define snd second
typedef pair<int,int> PII;
int dist(PII a, PII b) { return max(abs(a.fst-b.fst),abs(a.snd-b.snd)); }
void process(PII a) { printf("%d %d\n", a.fst,a.snd); }

template <class T>
struct bk_tree {
  typedef int dist_type;
  struct node {
    T p;
    unordered_map<dist_type, node*> ch;
  } *root;
  bk_tree() : root(0) { }

  node *insert(node *n, T p) {
    if (!n) { n = new node(); n->p = p; return n; }
    dist_type d = dist(n->p, p);
    n->ch[d] = insert(n->ch[d], p);
    return n;
  }
  void insert(T p) { root = insert(root, p); }
  void traverse(node *n, T p, dist_type dmax) {
    if (!n) return;
    dist_type d = dist(n->p, p);
    if (d < dmax) {
      process(n->p); // write your process
    }
    for (auto i: n->ch) 
      if (-dmax <= i.fst - d && i.fst - d <= dmax) 
        traverse(i.snd, p, dmax);
  }
  void traverse(T p, dist_type dmax) { traverse(root, p, dmax); }
};

int main() {
 
  bk_tree<PII> B;
  for (int i = 0; i < 10; ++i)
    for (int j = 0; j < 10; ++j)
      B.insert(PII(i,j));

  B.traverse( PII(0,0), 2 );

  cout << endl;

  B.traverse( PII(3,3), 2 );

  cout << endl;
}

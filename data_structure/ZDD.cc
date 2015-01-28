// 
// Zero-suppressed binary decision diagram with family algebra operations
//
// Description: 
//   ZDD maintains a set of sets of integers. Possible operations are
//     - cup(A, B): union of A and B              { 1, 12, 123 } cup { 1, 23 } = { 1, 12, 123, 23 } 
//     - cap(A, B): intersection of A and B       { 1, 12, 123 } cap { 12 } = { 12 }
//     - sub(A, B): difference of A and B         { 1, 12, 123 } sub { 12 } = { 1, 123 }
//     - mul(A, B): Kronecker product of A and B, { 1, 12, 123 } mul { 34 } = { 134, 1234 }
//     - div(A, B): quotient of A and B,          { 1, 12, 123 } div { 12 } = { {}, 3 }
//     - mod(A, B): reminder of A and B           { 1, 12, 123 } mod { 12 } = { 1 }
//     - count(A) : the number of sets in A
//
// Complexity:
//   O( #essentially disjoint subsets ).
//
// References: 
//   S. Minato (1993): 
//     Zero-suppressed BDDs for set manipulation in combinatorial problems.
//     Proceedings of the 30st annual Design Automation Conference, pp. 272-277.
//   S. Minato (1994):
//     Calculation of unate cube set algebra using zero-suppressed BDDs.
//     Proceedings of the 31st annual Design Automation Conference, pp. 420-424.
// 
#include <iostream>
#include <vector>
#include <cstdio>
#include <cstdlib>
#include <map>
#include <cmath>
#include <set>
#include <cstring>
#include <functional>
#include <algorithm>

using namespace std;
typedef long long LL;

#define ALL(c) c.begin(), c.end()
#define FOR(i,c) for(typeof(c.begin())i=c.begin();i!=c.end();++i)
#define REP(i,n) for(int i=0;i<n;++i)
#define fst first
#define snd second


// ZDD (Family Algebra)
struct ZDD {
  struct Node { 
    int v, lo, hi; 
    bool operator<(Node y) const { 
      return v != y.v ? v < y.v : lo != y.lo ? lo < y.lo : hi < y.hi;
    }
  };
  vector<Node> node;
  static const int zero = 0, unit = 1;
  ZDD() : node(2, (Node){-1,0,0}) { } // 0 = \emptyset, 1 = { \emptyset }

  int getnode(int v, int lo, int hi) {
    static map<Node, int> H;
    if (!hi) return lo;
    Node t = {v, lo, hi};
    if (H.count(t)) return H[t];
    node.push_back(t);
    return H[t] = node.size()-1;
  }
  int var(int v) { return getnode(v, 0, 1); }
  int cup(int x, int y) {
    static map<pair<int, int>, int> H;
    if (x > y) swap(x, y);
    if (!x || x == y) return y;
    if (H.count(make_pair(x,y))) return H[make_pair(x,y)];
    return H[make_pair(x,y)] = 
      node[x].v > node[y].v ? getnode(node[x].v, cup(node[x].lo, y), node[x].hi) :
      node[x].v < node[y].v ? getnode(node[y].v, cup(node[y].lo, x), node[y].hi) :
      getnode(node[x].v, cup(node[x].lo, node[y].lo), cup(node[x].hi, node[y].hi));
  }
  int cap(int x, int y) {
    static map<pair<int, int>, int> H;
    if (x > y) swap(x, y);
    if (!x || x == y) return x;
    if (H.count(make_pair(x,y))) return H[make_pair(x,y)];
    return H[make_pair(x,y)] = 
      node[x].v > node[y].v ? cap(node[x].lo, y) :
      node[x].v < node[y].v ? cap(node[y].lo, x) :
      getnode(node[x].v, cap(node[x].lo, node[y].lo), cap(node[x].hi, node[y].hi));
  }
  int sub(int x, int y) { 
    static map<pair<int, int>, int> H;
    if (!x || !y) return x;
    if (x == y) return 0;
    if (H.count(make_pair(x,y))) return H[make_pair(x,y)];
    return H[make_pair(x,y)] = 
      node[x].v > node[y].v ? getnode(node[x].v, sub(node[x].lo, y), node[x].hi) :
      node[x].v < node[y].v ? sub(x, node[y].lo) :
      getnode(node[x].v, sub(node[x].lo, node[y].lo), sub(node[x].hi, node[y].hi));
  }
  int mul(int x, int y) { 
    static map<pair<int, int>, int> H;
    if (x > y) swap(x, y);
    if (!x || y == 1) return x;
    if (!y || x == 1) return y;
    if (H.count(make_pair(x,y))) return H[make_pair(x,y)];
    return H[make_pair(x,y)] = 
      node[x].v > node[y].v ? getnode(node[x].v, mul(node[x].lo, y), mul(node[x].hi, y)) :
      node[x].v < node[y].v ? getnode(node[y].v, mul(node[y].lo, x), mul(node[y].hi, x)) :
      getnode(node[x].v, mul(node[x].lo, node[y].lo), cup(cup(
          mul(node[x].hi, node[y].hi), mul(node[x].hi, node[y].lo)), mul(node[x].lo, node[y].hi)));
  }
  int div(int x, int y) {
    static map<pair<int, int>, int> H;
    if (y == 1) return x;
    if (x <= 1) return 0;
    if (x == y) return 1;
    if (node[x].v < node[y].v) return 0; // node[y].v does not occur in x
    if (H.count(make_pair(x,y))) return H[make_pair(x,y)];
    if (node[x].v != node[y].v) 
      return H[make_pair(x,y)] = getnode(node[x].v, div(node[x].lo, y), div(node[x].hi, y));
    int z = div(node[x].hi, node[y].hi);
    return H[make_pair(x,y)] = z && node[y].lo ? cap(z, div(node[x].lo, node[y].lo)) : z;
  }
  int mod(int x, int y) { return sub(x, mul(y,div(x,y))); }

  LL count(int x) {
    static map<int, LL> H;
    if (x <= 1) return x;
    if (H.count(x)) return H[x];
    return H[x] = count(node[x].lo) + count(node[x].hi);
  }

  void disp(int x, vector<int> p = vector<int>()) {
    if (x == 1) {
      for (int i = 0; i < p.size(); ++i) cout << p[i] << ".";
      cout << "_ ";
    } else if (x != 0) { 
      p.push_back(node[x].v);
      disp(node[x].hi, p);
      p.pop_back();
      disp(node[x].lo, p);
    }
  }
};


int main() {
  ZDD zdd;

  int x = zdd.unit;
  for (int i = 0; i < 20; ++i) 
    x = zdd.cup(x, zdd.var(i));

  int y = zdd.unit;
  for (int i = 10; i < 30; ++i) 
    y = zdd.cup(y, zdd.var(i));

  int a = zdd.unit;
  for (int k = 0; k < 5; ++k) a = zdd.mul(a, x);
  for (int k = 0; k < 5; ++k) a = zdd.mul(a, y);

  cout << zdd.count(a) << endl;
}

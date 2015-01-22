// 
// Two-Sat
//
// Description:
//   In a SAT problem, if all clauses contains at most two literals,
//   this problem is called a two-sat, and solved in linear time.
//
// Algorithm:
//   For simplicity, we assume that all clauses contain exactly
//   two different literals (otherwise, we need a unit propagation).
//   We construct a graph G = (V, E), where V is a set of all literals 
//   and the negatives (i.e., |V| = 2n).
//   For a clause (u, v), we add a link (~u -> v) and (~v -> u).
//   Compute a strong connected component in this graph,
//   and if there ie u such that u and ~u are in the same component, 
//   the CNF cannot satisfied.
//
// Complexity:
//   O(n+m)
//   Here, n is the number of literals and is the number of clauses.
//
// References:
//   B. Aspvall, M. Plass, R. E. Tarjan (1979):
//   A linear-time algorithm for testing the truth of certain quantified boolean formulas.
//   Information Processing Letters, vol. 8, no. 3, pp. 121-123.

#include <iostream>
#include <vector>
#include <cstdio>
#include <unordered_set>
#include <algorithm>
#include <cstring>
#include <queue>

using namespace std;

#define fst first
#define snd second
#define all(c) ((c).begin()), ((c).end())

#define dout if(1){}else cout


// positive literal  x in [0,n), 
// negative literal ~x in [-n,0)
struct satisfiability_twosat {
  int n;
  vector<int> x, num, col, ord;
  vector<vector<int>> imp, rmp;
  satisfiability_twosat(int n) : 
    n(n), x(2*n), num(2*n), col(2*n), imp(2*n), rmp(2*n) { } 
  void add_clause(int u, int v) { // assert: u != v 
    if      (u == ~v) return;   // ignore trivial clause
    else if (u ==  v) x[u+n] = 1; // pure clause ==> set
    else {
      imp[~u+n].push_back(v);
      imp[~v+n].push_back(u);
      rmp[u+n].push_back(~v);
      rmp[v+n].push_back(~u);
    }
  }
  void visit(int u, bool b) {
    if (num[u+n]++ > 0) return;
    x[u+n] |= b;
    for (auto v: imp[u+n]) visit(v, x[u+n]);
    ord.push_back(u);
  }
  void rvisit(int u, int k) {
    if (num[u+n]++ > 0) return;
    col[u+n] = k;
    for (auto v: rmp[u+n]) rvisit(v, k);
  }
  bool solve() {
    for (int u = -n; u < n; ++u) 
      if (x[u+n]) visit(u, 0);
    for (int u = -n; u < n; ++u) 
      if (x[u+n] && x[~u+n]) return false;
      else visit(u, 0);
    num = col; 
    reverse(all(ord)); 
    for (int v: ord) rvisit(v, v);
    for (int u = 0; u < n; ++u) 
      if (col[u+n] == col[~u+n]) return false;
    return true;
  }
};



// for verification
struct satisfiability {
  int n;
  vector<int> x, occ, pos, neg;
  vector<vector<int>> adj, lit;
  vector<vector<int>> decision_stack;
  vector<int> unit_stack, pure_stack;
  satisfiability(int n) : 
    n(n), x(2*n), adj(2*n), occ(2*n), decision_stack(1) { }
  void add_clause(vector<int> c) { 
    sort(all(c)); c.erase(unique(all(c)), c.end());
    for (int i = 0; i < c.size(); ++i)
      if (binary_search(all(c), ~c[i])) return;
    for (auto u: c) {
      adj[u+n].push_back(lit.size());
      occ[u+n] += 1;
    }
    lit.push_back(c);
    pos.push_back(0);
    neg.push_back(0);
  } 
  void push(int u) {
    x[u+n] = 1;
    decision_stack.back().push_back(u);
    for (auto i: adj[ u+n])
      if (pos[i]++ == 0) 
        for (auto u: lit[i])
          --occ[u+n];
    for (auto i: adj[~u+n]) { 
      ++neg[i];
      if (pos[i] == 0) unit_stack.push_back(i); 
    }
  }
  void pop() {
    int u = decision_stack.back().back();
    decision_stack.back().pop_back();
    x[u+n] = 0;
    for (auto i: adj[ u+n])
      if (--pos[i] == 0)
        for (auto u: lit[i])
          ++occ[u+n];
    for (auto i: adj[~u+n]) --neg[i];
  }
  bool reduction() {
    while (!unit_stack.empty() || !pure_stack.empty()) {
      if (!pure_stack.empty()) { // pure literal elimination
        int u = pure_stack.back();
        pure_stack.pop_back();
        if (occ[u+n] == 1 && occ[~u+n] == 0) push(u);
      } else {                   // unit propagation
        int i = unit_stack.back();
        unit_stack.pop_back();
        if (pos[i] > 0) continue;
        if (neg[i]     == lit[i].size()) return false;
        if (neg[i] + 1 == lit[i].size()) {
          int w = n;
          for (int u: lit[i]) if (!x[u+n] && !x[~u+n]) w = u;
          if (x[~w+n]) return false;
          push(w);
        }
      }
    }
    return true;
  }
  bool solve() {
    while (1) { 
      if (reduction()) {
        int s = 0;
        for (int u = 0; u < n; ++u)
          if (occ[s+n]+occ[~s+n] < occ[u+n]+occ[~u+n]) s = u;
        if (occ[s+n] + occ[~s+n] == 0) return true;
        decision_stack.push_back({});
        push(s);
      } else {
        int s = decision_stack.back()[0];
        while (!decision_stack.back().empty()) pop();
        decision_stack.pop_back();
        if (decision_stack.empty()) return false;
        push(~s);
      }
    }
  }
};



int main() {
  // check: 2-sat
  for (int iter = 0; iter < 10000; ++iter) {
    srand(iter);
    int n = 100, m = 150;
    satisfiability_twosat twosat(n);
    satisfiability gen(n);
    for (int i = 0; i < m; ++i) {
      int u = rand() % n;
      int v = rand() % n;
      if (rand() % 2) u = ~u;
      if (rand() % 2) v = ~v;
      twosat.add_clause(u, v);
      gen.add_clause({u, v});
      if (u >= 0) dout << u << " or ";
      else        dout << " ~" << ~u << " or ";
      if (v >= 0) dout << v << endl;
      else        dout << "~" << ~v << endl;
    }

    dout << "seed = " << iter << endl;
    dout << "begin two" << endl;
    int a = twosat.solve();
    dout << "begin gen" << endl;
    int b = gen.solve();
    //cout << a << " " << b << endl;

    if (a != b) {
      cout << "seed = " << iter << endl;
      cout << "2sat = " << a << " / gen = " << b << endl;
      break;
    }
  }
}

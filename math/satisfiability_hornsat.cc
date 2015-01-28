// 
// Horn-SAT
//
// Description:
//   In a SAT problem, if all clauses contains at most one positive literals,
//   this problem is called a Horn-SAT, and solved in linear time. 
//
// Algorithm:
//   If there is no singleton clause, by assigning false for all literals,
//   we obtain a solution. Otherwise, let {u} be a singleton clause.
//   Then, u must be positive. We eliminate all clauses that contains u
//   and remove ~u from all clauses containing ~u. We repeat this procedure.
//   Note that this procedure is implemented as similar as the BFS.
//
// Complexity:
//   O(n+m)
//   Here, n is the number of literals and is the number of clauses.
//
// References:
//   W. Dowling and J. Gallier (1984) 
//   Linear-time algorithms for testing the satisfiability of propositional Horn formulae.
//   Journal of Logic Programming, vol. 3, pp. 267-284.
 
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

struct satisfiability_hornsat {
  int n, T, F;
  struct clause { int num, u; }; // num of neg vars, pos var.
  vector<clause> cs;
  vector<vector<int>> adj;
  satisfiability_hornsat(int n) : n(n), adj(n) { }
  void add_clause(vector<int> c) {
    sort(all(c)); c.erase(unique(all(c)), c.end());
    for (int i = 0; i < c.size(); ++i)
      if (binary_search(all(c), ~c[i])) return;
    int u = n;
    for (auto a: c) {
      if (a >= 0) adj[u=a].push_back(cs.size());
      else        adj[ ~a].push_back(cs.size());
    }
    cs.push_back({(int)c.size() - (u != n), u});
  }
  bool solve() {
    vector<int> front, x(n);
    for (int i = 0; i < cs.size(); ++i) 
      if (cs[i].num == 0)
        front.push_back(i);
    while (!front.empty()) {
      int u = cs[front.back()].u; front.pop_back();
      if (x[u]) continue;
      if (u == n) return false;
      x[u] = true;
      for (int i: adj[u]) 
        if (--cs[i].num == 0) 
          front.push_back(i);
    }
    return true;
  }
};

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
  // check: horn-sat
  for (int iter = 0; iter < 10; ++iter) {
    srand(iter);
    int n = 50, m = 200;
    satisfiability_hornsat horn(n);
    satisfiability gen(n);
    
    vector<vector<int>> in;
    for (int i = 0; i < m; ++i) {
      int k = 0;
      for (int j = 0; j < 10; ++j)
        k = k + rand() % 2;
      vector<int> c;
      for (int j = 0; j < k; ++j) {
        int u = rand() % n;
        c.push_back(~u);
      }
      if (rand() % 2) {
        int u = rand() % n;
        c.push_back(u);
      }
      if (c.size() > 0) 
        in.push_back(c);
    }

    for (auto v: in) {
      horn.add_clause(v);
      gen.add_clause(v);
    }

    for (auto v: in) {
      for (int i = 0; i < v.size(); ++i) {
        if (i > 0) dout << " or ";
        int a = v[i];
        if (a >= 0) dout << a;
        else        dout << "~" << ~a;
      }
      dout << endl;
    }

    dout << "seed = " << iter << endl;
    dout << "begin horn" << endl;
    int a = horn.solve();
    dout << "begin gen" << endl;
    int b = gen.solve();

    if (a != b) {
      cout << "!!!" << endl;
      cout << "seed = " << iter << endl;
      cout << "horn = " << a << " / gen = " << b << endl;
      for (auto v: in) {
        for (int i = 0; i < v.size(); ++i) {
          if (i > 0) cout << " or ";
          int a = v[i];
          if (a >= 0) cout << a;
          else        cout << "~" << ~a;
        }
        cout << endl;
      }
      break;
    }
  }
}

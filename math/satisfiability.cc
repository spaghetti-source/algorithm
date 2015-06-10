// 
// SAT solver (DPLL method)
//
// Description
//   We are given a CNF, e.g., 
//     \phi(x) = (x_1 or ~x_2) and (x_3 or ~x_4 or ~x_5) and ... .
//   SAT finds an assignment x for phi(x) = true.
//
// Algorithm:
//   Davis-Putnum-Logemann-Loveland's algorithm that performs 
//     1) unit propagation, 
//     2) pure literal elimination,
//     3) branch and search.
//  
// Complexity:
//   O(2^n) in worst case. 
//   This implementation is practical for n = 100.
//
// Verified:
//   Uniform 3-SAT instances from http://www.cs.ubc.ca/~hoos/SATLIB/benchm.html


#include <iostream>
#include <vector>
#include <cstdio>
#include <algorithm>
#include <cstring>
#include <functional>

using namespace std;

#define fst first
#define snd second
#define all(c) ((c).begin()), ((c).end())
#define dout if(0){}else cout

// positive literal  x in [0,n), 
// negative literal ~x in [-n,0)
//
struct satisfiability {
  int n;
  vector<int> occ, pos, neg;
  vector<vector<int>> adj, lit;
  satisfiability(int n) : n(n), adj(2*n), occ(2*n) { }
  void add_clause(const vector<int> &c) { 
    for (auto u: c) {
      adj[u+n].push_back(lit.size());
      occ[u+n] += 1;
    }
    lit.push_back(c);
  } 
  vector<bool> x;
  vector<vector<int>> decision_stack;
  vector<int> unit_stack, pure_stack;
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
    x.assign(2*n,0);
    pos = neg = vector<int>(lit.size());
    decision_stack.assign(1,{});
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



void moin(int seed) {
#if 1
  vector<vector<int>> input;
  int n;
  for (char s[1024]; fgets(s, sizeof(s), stdin); ) {
    if (s[0] == 'c') continue;
    if (s[0] == 'p') {
      sscanf(s, "p cnf %d %*d", &n);
    } else if (s[0] == '%') {
      break;
    } else {
      int u, v, w;
      sscanf(s, "%d %d %d", &u, &v, &w);
      input.push_back({u, v, w});
    }
  }
  satisfiability solver(n);
  for (auto &c: input) {
    for (auto &a: c) 
      if (a >= 0) a -= 1;
    solver.add_clause(c);
  }
  int x = solver.solve();
  cout << x << endl;
#else
  srand(seed);
  satisfiability solver(9);
  solver.add_clause({0,1,5});
  solver.add_clause({0,1,6});
  solver.add_clause({0,1,7});
  solver.add_clause({0,1,8});
  cout << solver.solve() << endl;
#endif
}

int main() {
  moin(0);
}

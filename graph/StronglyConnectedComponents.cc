//
// Strongly Connected Components
//
// Description:
//   Compute the strongly connected components(SCCs) decomposition.
//   A set of vertices S of a digraph is a SCC iff 
//   for each u, v in S, there is a path from u to v.
//   Any digraph can be uniquely split into SCCs.
//
//
// Algorithm:
//   Gabow's one-pass / two-stacks algorithm.
//
//
// Complexity:
//   O(n + m) time and space.
//
//
// Verified:
//   SPOJ 6818: Capital City
//
//
// References: 
//
// - H. N. Gabow (2000):
//   Path-based depth first search strong and biconnected components.
//   Information Processing Letters, vol.74, no.3-4, pp.107-114.
//
#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <vector>
#include <cstring>
#include <functional>
#include <algorithm>
#include <stack>

using namespace std;

#define ALL(c) c.begin(), c.end()
#define FOR(i,c) for(typeof(c.begin())i=c.begin();i!=c.end();++i)
#define REP(i,n) for(int i=0;i<n;++i)
#define fst first
#define snd second

struct Graph {
  int n;
  vector< vector<int> > adj;
  stack<int> S, B;
  vector<int> I;
  Graph(int n) : n(n), adj(n), I(n) { }
  void addEdge(int u, int v) {
    adj[u].push_back(v); // directed
  }
  void visit(int v) {
    B.push(I[v] = S.size());
    S.push(v);
    FOR(e, adj[v]) { 
      int w = *e;
      if (!I[w]) visit(w);
      else while (I[w] < B.top()) B.pop();
    }
    if (I[v] == B.top()) {
      B.pop();
      for (; I[v] < S.size(); S.pop()) 
        cout << S.top() << " ";
      cout << endl;
    }
  }
  void stronglyConnectedComponents() {
    REP(v, n) if (!I[v]) visit(v);
  }
};

int main() { }

// usage
namespace SPOJ6818 {
struct Graph {
  int n;
  vector< vector<int> > adj;
  stack<int> S, B;
  vector<int> I;
  int C;
  vector<int> id; 
  Graph(int n) : n(n), adj(n), I(n), C(0), id(n) { }
  void addEdge(int u, int v) {
    adj[u].push_back(v); // directed
  }
  void visit(int v) {
    B.push(I[v] = S.size());
    S.push(v);
    FOR(e, adj[v]) { 
      int w = *e;
      if (!I[w]) visit(w);
      else while (I[w] < B.top()) B.pop();
    }
    if (I[v] == B.top()) {
      B.pop();
      for (; I[v] < S.size(); S.pop()) 
        id[ S.top() ] = C;
      ++C;
    }
  }
  void stronglyConnectedComponents() {
    REP(v, n) if (!I[v]) visit(v);
  }

  // SPOJ 6818
  void solve() {
    stronglyConnectedComponents();
    vector<int> outdeg(C);
    REP(u, n) FOR(e, adj[u]) 
      if (id[u] != id[*e]) ++outdeg[ id[u] ];
    int c = find(ALL(outdeg), 0) - outdeg.begin();
    if (c == outdeg.size()) {
      printf("0\n");
      return;
    }
    vector<int> component;
    REP(u, n) if (id[u] == c) component.push_back(u);
    printf("%d\n%d", component.size(), component[0]+1);
    for (int i = 1; i < component.size(); ++i) 
      printf(" %d", component[i]+1);
    printf("\n");
  }
};

int main() {
  int n, m;
  scanf("%d %d", &n, &m);
  Graph G(n);
  REP(i, m) {
    int src, dst;
    scanf("%d %d", &src, &dst);
    G.addEdge(src-1, dst-1);
  }
  G.solve();
}
}

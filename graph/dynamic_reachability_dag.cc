// 
// Italiano's dynamic reachability data structure for DAG
//
// Description:
//   It is a data structure that admits the following operations:
//     add_edge(s, t):     insert edge (s,t) to the network if
//                         it does not make a cycle
//
//     is_reachable(s, t): return true iff there is a path s --> t
//
// Algorithm:
//   We maintain reachability trees T(u) for all u in V.
//   Then is_reachable(s, t) is solved by checking "t in T(u)".
//   For add_edge(s, t), if is_reachable(s, t) or is_reachable(t, s) then
//   no update is performed. Otherwise, we meld T(s) and T(t).
//
// Complexity:
//   amortized O(n) per update
//
// Verified:
//   SPOJ 9458: Ghosts having fun
//
// References:
//   Giuseppe F. Italiano (1988):
//   Finding paths and deleting edges in directed acyclic graphs.
//   Information Processing Letters, vol. 28, no. 1, pp. 5--11.
//   

#include <iostream>
#include <vector>
#include <cstdio>
#include <algorithm>
#include <functional>

using namespace std;

#define fst first
#define snd second
#define all(c) ((c).begin()), ((c).end())

struct dag_reachability {
  int n;
  vector<vector<int>> parent;
  vector<vector<vector<int>>> child;
  dag_reachability(int n) : n(n), parent(n, vector<int>(n, -1)), child(n, vector<vector<int>>(n)) { }

  bool is_reachable(int src, int dst) {
    return src == dst || parent[src][dst] >= 0;
  }
  bool add_edge(int src, int dst) {
    if (is_reachable(dst, src)) return false; // break DAG condition
    if (is_reachable(src, dst)) return true;  // no-modification performed
    for (int p = 0; p < n; ++p) 
      if (is_reachable(p, src) && !is_reachable(p, dst)) 
        meld(p, dst, src, dst);
    return true;
  }
  void meld(int root, int sub, int u, int v) {
    parent[root][v] = u;
    child[root][u].push_back(v);
    for (int c: child[sub][v]) 
      if (!is_reachable(root, c)) 
        meld(root, sub, v, c);
  }
};
// SPOJ 9458: Ghosts having fun
int main() {
  int n, m;
  scanf("%d %d", &n, &m);
  dag_reachability g(n);
  for (int i = 0; i < m; ++i) {
    int s, t;
    scanf("%d %d", &s, &t);
    if (!g.add_edge(s-1, t-1))
      printf("%d %d\n", s, t);
  }
  printf("0 0\n");
}



// test routine
struct dag_reachability_n {
  int n;
  vector<vector<int>> adj;
  dag_reachability_n(int n) : n(n), adj(n) { }

  bool is_reachable(int src, int dst) {
    vector<int> visited(n); 
    function<bool (int)> visit = [&](int u) {
      visited[u] = true;
      if (u == dst) return true;
      for (int v: adj[u]) 
        if (!visited[v]) 
          if (visit(v)) return true;
      return false;
    };
    return visit(src);
  }
  bool add_edge(int src, int dst) {
    if (is_reachable(dst, src)) return false; // break DAG condition
    if (is_reachable(src, dst)) return true;  // no-modification performed
    adj[src].push_back(dst);
    return true;
  }
};

void verify() {
  int n = 50;
  dag_reachability g(n);
  dag_reachability_n h(n);
  for (int i = 0; i < 4000; ++i) {
    int s, t;
    while (1) {
      s = rand() % n;
      t = rand() % n;
      if (s != t) break;
    }
    int a = g.add_edge(s, t);
    int b = h.add_edge(s, t);
    if (a != b) {
      cout << "*****" << endl;
      exit(0);
    }
  }
}

// === tick a time ===
#include <ctime>
double tick() {
  static clock_t oldtick;
  clock_t newtick = clock();
  double diff = 1.0*(newtick - oldtick) / CLOCKS_PER_SEC;
  oldtick = newtick;
  return diff;
}
void speed() {
  int n = 1000, m = 300000;
  int a = 0;
  vector<pair<int,int>> query;
  for (int i = 0; i < m; ++i) {
    int s, t;
    while (1) {
      s = rand() % n;
      t = rand() % n;
      if (s != t) break;
    }
    query.push_back({s, t});
  }
  {
    tick();
    dag_reachability g(n);
    for (int i = 0; i < m; ++i) {
      int s = query[i].fst, t = query[i].snd;
      a += g.add_edge(s, t);
    }
  }
  cout << "italiano " << tick() << endl;
  cout << a << endl;
  {
    dag_reachability_n g(n);
    for (int i = 0; i < m; ++i) {
      int s = query[i].fst, t = query[i].snd;
      a -= g.add_edge(s, t);
    }
  }
  cout << "naive " << tick() << endl;
  cout << a << endl;
}
void test() {
  verify();
  speed();
}

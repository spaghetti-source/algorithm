//
// Chordal graph recognition
//
// Description:
//   A graph is chordal if any cycle C with |C| >= 4 has a chord. 
//   It is also characterized by a perfect elimination ordering (PEO):
//   An ordering pi is a PEO if, for all u in V, 
//     B(u) := {u} \cup { v in N(u) : pi(v) > pi(u) } 
//   forms a clique. A graph is chordal if and only if it admits a PEO.
//
//   Many problems on a chordal graph can be solved in P by using PEO.
// 
// Algorithm:
//   We find a PEO to regocnize a chordal graph. 
//   The maximum cardinality search (MCS), which iterates
//     select u in V with largest |N(u) \cap selected|,
//   gives a PEO [Rose and Tarjan'75]; thus we perform MCS and 
//   then verify whether it is PEO or not.
//
// Complexity:
//   O(n+m) time and space
//
// References
//   D. J. Rose and R. Endre Tarjan (1975):
//   Algorithmic aspects of vertex elimination.
//   In Proceedings of the 7th Annual ACM Symposium on Theory of Computing (STOC'75),
//   pp. 245--254.
//   

#include <iostream>
#include <vector>
#include <list>
#include <cstdio>
#include <algorithm>
#include <functional>
#include <unordered_set>
#include <unordered_map>
#include <map>
#include <set>

using namespace std;

#define fst first
#define snd second
#define all(c) ((c).begin()), ((c).end())

struct graph {
  int n;
  vector<vector<int>> adj;
  graph(int n) : n(n), adj(n) { }
  void add_edge(int src, int dst) {
    adj[src].push_back(dst);
    adj[dst].push_back(src);
  }
};
bool is_chordal(graph g) {
  int top = 0;
  vector<int> rank(g.n, -1), score(g.n);
  vector<unordered_set<int>> R(g.n);
  vector<vector<int>> bucket(g.n); // bucket Dijkstra
  for (int i = 0; i < g.n; ++i)
    bucket[0].push_back(i);
  for (int i = 0; i < g.n; ) {
    while (bucket[top].empty()) --top;
    int u = bucket[top].back(); bucket[top].pop_back();
    if (rank[u] >= 0) continue;
    rank[u] = i++;
    int p = -1;
    for (int v: g.adj[u]) {
      if (rank[v] >= 0) {
        R[u].insert(v);
        if (p == -1 || rank[p] < rank[v]) p = v;
      } else {
        bucket[++score[v]].push_back(v);
        top = max(top, score[v]);
      }
    }
    if (p >= 0) 
      for (int v: R[u]) 
        if (v != p && !R[p].count(v)) return false;
  }
  return true;
}

int maximum_clique_chordal(graph g) {
  int top = 0;
  vector<int> pi(g.n), rank(g.n, -1), score(g.n);
  vector<vector<int>> bucket(g.n); // bucket dijkstra
  for (int i = 0; i < g.n; ++i)
    bucket[0].push_back(i);
  for (int i = 0; i < g.n; ) {
    while (bucket[top].empty()) --top;
    int u = bucket[top].back(); bucket[top].pop_back();
    if (rank[u] >= 0) continue;
    pi[i] = u; rank[u] = i++;
    for (int v: g.adj[u]) {
      if (rank[v] >= 0) continue;
      bucket[++score[v]].push_back(v);
      top = max(top, score[v]);
    }
  }
  int ans = 0;
  unordered_set<int> processed;
  vector<unordered_set<int>> R(g.n);
  for (int u: pi) {
    int p = -1;
    for (int v: g.adj[u]) {
      if (!processed.count(v)) continue;
      R[u].insert(v);
      if (p == -1 || rank[p] < rank[v]) p = v;
    }
    ans = max(ans, (int)R[u].size());
    processed.insert(u);
  }
  return ans;
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

bool test1() {
  graph g(3);
  g.add_edge(0, 1);
  g.add_edge(1, 2);
  g.add_edge(2, 0);
  return is_chordal(g);
}
bool test2() {
  graph g(4);
  g.add_edge(0, 1);
  g.add_edge(1, 2);
  g.add_edge(2, 3);
  g.add_edge(3, 0);
  return !is_chordal(g);
}
bool test3() {
  graph g(4);
  g.add_edge(0, 1);
  g.add_edge(1, 2);
  g.add_edge(2, 3);
  g.add_edge(3, 0);
  g.add_edge(3, 1);
  return is_chordal(g);
}
bool test4() {
  graph g(5);
  g.add_edge(0, 1);
  g.add_edge(1, 2);
  g.add_edge(2, 3);
  g.add_edge(3, 0);
  g.add_edge(4, 0);
  g.add_edge(4, 1);
  g.add_edge(4, 2);
  g.add_edge(4, 3);
  return !is_chordal(g);
}
bool test5() {
  graph g(5);
  g.add_edge(0, 1);
  g.add_edge(1, 2);
  g.add_edge(2, 3);
  g.add_edge(3, 4);
  g.add_edge(4, 0);
  g.add_edge(4, 1);
  g.add_edge(4, 2);
  return is_chordal(g);
}
bool test6() {
  int n = 1000;
  graph g(n);
  for (int i = 0; i < n; ++i)
    for (int j = i+1; j < n; ++j)
      if (rand() % 100) g.add_edge(i, j);
  tick();
  is_chordal(g);
  cout << tick() << endl;
  return true;
}

int main() {
  if      (!test1()) cout << "test1 failed" << endl;
  else if (!test2()) cout << "test2 failed" << endl;
  else if (!test3()) cout << "test3 failed" << endl;
  else if (!test4()) cout << "test4 failed" << endl;
  else if (!test5()) cout << "test5 failed" << endl;
  else if (!test6()) cout << "test6 failed" << endl;
  else cout << "all tests passed" << endl;
}

//
// Tree Isomorphism
//
// Description:
//   Rooted trees (S,s) and (T,t) are isomorphic
//   iff there is a bijection between childs of s and childs of t
//   such that s.child[i] and t.child[pi[i]] is isomorphic.
//
//   Two trees are isomorphic iff these are isomorphic with some roots.
//
// Algorithm:
//   Aho-Hopcroft-Ullmann's algorithm for rooted isomorphism.
//   It simultaneously scans the vertices in the two trees from bottom to root,
//   and assigns unique id (among the level) for each subtrees.
//
//   For unrooted isomorphism, it first finds centers of tree,
//   and try all possibility of rooted tree isomorphism.
//   Since the number of centers in a tree is at most two,
//   it can be performed as the same complexity as rooted isomorphism.
//
// Complexity:
//   O(n log n). 
//
// Verified:
//   SPOJ7826.
//
// References:
//   A. V. Aho, J. E. Hopcroft, and J. D. Ullman (1974):
//   The Design and Analysis of Computer Algorithms.
//   Addison-Wesley.


#include <iostream>
#include <algorithm>
#include <vector>
#include <queue>
#include <map>
#include <cstdio>
#include <cstring>

using namespace std;
#define fst first
#define snd second
#define all(c) ((c).begin()), ((c).end())

struct tree {
  int n;
  vector<vector<int>> adj;
  tree(int n) : n(n), adj(n) { }
  void add_edge(int src, int dst) {
    adj[src].push_back(dst);
    adj[dst].push_back(src);
  }
  vector<int> centers() {
    vector<int> prev;
    int u = 0;
    for (int k = 0; k < 2; ++k) { // double sweep
      queue<int> que;
      prev.assign(n, -1);
      que.push(prev[u] = u); 
      while (!que.empty()) {
        u = que.front(); que.pop();
        for (auto v: adj[u]) {
          if (prev[v] >= 0) continue;
          que.push(v);
          prev[v] = u;
        }
      }
    }
    vector<int> path = {u}; // median on a path
    while (u != prev[u]) 
      path.push_back(u = prev[u]);
    int m = path.size(); 
    if (m % 2 == 0) return {path[m/2-1], path[m/2]};
    else            return {path[m/2]};
  }

  vector<vector<int>> layer;
  vector<int> prev;
  int levelize(int r) { // split vertices into levels
    prev.assign(n,-1); prev[r] = n;
    layer = {{r}};
    while (1) {
      vector<int> next;
      for (int u: layer.back()) {
        for (int v: adj[u]) {
          if (prev[v] >= 0) continue;
          prev[v] = u;
          next.push_back(v);
        }
      }
      if (next.empty()) break;
      layer.push_back(next);
    }
    return layer.size();
  }
};

bool isomorphic(tree S, int s, tree T, int t) {
  if (S.n != T.n) return false;
  if (S.levelize(s) != T.levelize(t)) return false;
  
  vector<vector<int>> longcodeS(S.n+1), longcodeT(T.n+1);
  vector<int> codeS(S.n), codeT(T.n);
  for (int h = S.layer.size()-1; h >= 0; --h) {
    map<vector<int>, int> bucket;
    for (int u: S.layer[h]) {
      sort(all(longcodeS[u]));
      bucket[longcodeS[u]] = 0;
    }
    for (int u: T.layer[h]) {
      sort(all(longcodeT[u]));
      bucket[longcodeT[u]] = 0;
    }
    int id = 0;
    for (auto &p: bucket) p.snd = id++;
    for (int u: S.layer[h]) {
      codeS[u] = bucket[longcodeS[u]];
      longcodeS[S.prev[u]].push_back(codeS[u]);
    }
    for (int u: T.layer[h]) {
      codeT[u] = bucket[longcodeT[u]];
      longcodeT[T.prev[u]].push_back(codeT[u]);
    }
  }
  return codeS[s] == codeT[t];
}
bool isomorphic(tree S, tree T) {
  auto x = S.centers(), y = T.centers();
  if (x.size() != y.size()) return false;
  if (isomorphic(S, x[0], T, y[0])) return true;
  return x.size() > 1 && isomorphic(S, x[1], T, y[0]);
}

void doit() {
  int n; scanf("%d", &n);
  tree S(n), T(n);
  for (int i = 0; i < n-1; ++i) {
    int u, v; scanf("%d %d", &u, &v);
    S.add_edge(u-1, v-1);
  }
  for (int i = 0; i < n-1; ++i) {
    int u, v; scanf("%d %d", &u, &v);
    T.add_edge(u-1, v-1);
  }
  if (isomorphic(S, T)) printf("YES\n");
  else                  printf("NO\n");
}
int main() {
  int ncase; scanf("%d", &ncase);
  for (int icase = 0; icase < ncase; ++icase) doit();
}


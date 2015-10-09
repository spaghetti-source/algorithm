//
// Edit Distance of two NFAs
//   
// Description:
//   Let L, M be languages. The edit distance of L and M
//   is defined by d(L,LM) = min { d(s,t) : s in L, t in M }.
//   Here, we assume that L and M is described by NFAs.
//
// Algorithm:
//   Consider a graph G = (V,E) such that
//     V = product of states in NFAs L and M
//     Suppose 
//       a: u --> v
//       b: s --> t
//     Then
//       (u,s) --> (v,s) with cost(a,empty)
//       (u,s) --> (u,t) with cost(empty,b)
//       (u,s) --> (v,t) with cost(a, b)
//   Then a path from (Lbegin,Mbegin) to (Lend,Mend)
//   corresponds to the correct edit operations.
//   We emulate the Dijkstra algorithm on this graph.
//
// Verified:
//   SPOJ 8056: Regex Edit Distance
//
// References:
//   Horst Bunke (1996): Edit Distance of Regular Languages,
//   in Proceedings of the 5th Annual Symposium on Document Analysis 
//   and Information Retrieval, pp.113--124.
//

#include <iostream>
#include <queue>
#include <vector>
#include <unordered_map> 
#include <unordered_set> 
#include <map>
#include <cstring>
#include <set>
#include <cstdio>
#include <bitset>

#include <algorithm>
#include <functional>

using namespace std;

#define fst first
#define snd second
#define all(c) ((c).begin()), ((c).end())
#define TEST(s) if (!s) { cout << __LINE__ << " " << #s << endl; exit(-1); }

const int NFA_STATE = 256, DFA_STATE = 500, ALPHA = 'c';
typedef bitset<NFA_STATE> subset;
struct NFA {
  static int size;
  static vector<int> next[NFA_STATE][ALPHA];
  static int new_node() {
    for (int a = 0; a < ALPHA; ++a) 
      next[size][a].clear();
    return size++;
  }
  static NFA symbol(char a) {
    int begin = new_node(), end = new_node();
    next[begin][a].push_back(end);
    return {begin, end};
  }
  static NFA unite(NFA x, NFA y) {
    int begin = new_node(), end = new_node();
    next[begin][0].push_back(x.begin);
    next[begin][0].push_back(y.begin);
    next[x.end][0].push_back(end);
    next[y.end][0].push_back(end);
    return {begin, end};
  }
  static NFA concat(NFA x, NFA y) {
    next[x.end][0].push_back(y.begin);
    return {x.begin, y.end};
  }
  static NFA star(NFA x) {
    int begin = new_node(), end = new_node();
    next[begin][0].push_back(x.begin);
    next[begin][0].push_back(end);
    next[x.end][0].push_back(x.begin);
    next[x.end][0].push_back(end);
    return {begin, end};
  }
  int begin, end;

  void closure(int u, subset &x) {
    x[u] = 1;
    for (int v: next[u][0]) 
      if (!x[v])
        closure(v, x);
  }
  bool run(const char *s) {
    subset x;
    closure(begin, x);
    for (; *s; ++s) {
      subset y;
      for (int u = 0; u < size; ++u) 
        if (x[u]) 
          for (int v: next[u][*s])
            closure(v, y);
      x = y;
    }
    return x[end];
  }
};
int NFA::size;
vector<int> NFA::next[NFA_STATE][ALPHA];

NFA parse(const char *s) {
  function<NFA ()> regex, factor, term;
  regex = [&]() {
    NFA a = factor();
    if (*s == '|') { ++s; a = NFA::unite(a, regex()); }
    return a;
  };
  factor = [&]() {
    NFA a = term();
    if (*s == '*') { a = NFA::star(a); ++s; }
    if (*s && *s != '|' && *s != ')') a = NFA::concat(a, factor());
    return a;
  };
  term = [&]() {
    if (*s == '(') { ++s; NFA a = regex(); ++s; return a; } 
    else { NFA a = NFA::symbol(*s); ++s; return a; }
  };
  return regex();
}


int edit_distance(NFA x, NFA y) {
  int dist[NFA::size][NFA::size];
  for (int i = 0; i < NFA::size; ++i)
    for (int j = 0; j < NFA::size; ++j)
      dist[i][j] = 99999999;
  struct node {
    int cost, u, v;
    bool operator<(node x) const { return cost > x.cost; }
  };
  priority_queue<node> que;
  que.push({0, x.begin, y.begin});
  dist[x.begin][y.begin] = 0;
  vector<int> A = {0, 'a', 'b'};
  while (!que.empty()) {
    node t = que.top(); que.pop();
    if (dist[t.u][t.v] < t.cost) continue;
    if (t.u == x.end && t.v == y.end) return t.cost;

    for (int a: A) {
      int len = (a == 0 ? 0 : 1); // insert a
      for (int u: x.next[t.u][a]) {
        if (dist[u][t.v] > t.cost + len) {
          dist[u][t.v] = t.cost + len;
          que.push({dist[u][t.v], u, t.v});
        }
      }
    }
    for (int b: A) {
      int len = (b == 0 ? 0 : 1); // delete b
      for (int v: y.next[t.v][b]) {
        if (dist[t.u][v] > t.cost + len) {
          dist[t.u][v] = t.cost + len;
          que.push({dist[t.u][v], t.u, v});
        }
      }
    }
    for (int a: A) {
      for (int b: A) {
        int len = (a == b ? 0 : 1);
        for (int u: x.next[t.u][a]) {
          for (int v: y.next[t.v][b]) {
            if (dist[u][v] > t.cost + len) {
              dist[u][v] = t.cost + len;
              que.push({dist[u][v], u, v});
            }
          }
        }
      }
    }
  }
  return -1; // something wrong
}


int main() {
  int ncase; scanf("%d", &ncase);
  for (int icase = 0; icase < ncase; ++icase) {
    NFA::size = 0;
    char s[1024], t[1024];
    scanf("%s %s", s, t);
    printf("%d\n", edit_distance(parse(s), parse(t)));
  }
}

//
// Suffix Tree (Ukkonen's algorithm)
//
// Description:
//   A suffix tree of a given string s is a 
//   patricia (i.e., compressed trie) of all suffixes of s.
//
// Algorithm:
//   Ukkonen's left-to-right scan algorithm.
//   See the original paper.
//
// Complexity:
//   O(n).
//   Due to a large constant factor, it is useful in
//   n <= 10000 in programming contest.
//   Use suffix array, instead of suffix tree.
//
// References:
//   E. Ukkonen (1995):
//   On-line Construction of Suffix Trees.
//   Algorithmica, vol.14, no.3, pp.249--260.

#include <iostream>
#include <vector>
#include <cstdio>
#include <cstring>
#include <algorithm>
#include <functional>

using namespace std;

#define fst first
#define snd second
#define all(c) ((c).begin()), ((c).end())


struct suffix_tree {
  const char *str;
  const int n;

  int root;
  struct edge { int to, a, b; };
  vector<vector<edge>> adj;
  vector<int> suf;
  int add_node() {
    adj.push_back(vector<edge>(0x100, {-1}));
    suf.push_back(0);
    return suf.size()-1;
  }
  int add_edge(int r, const edge &e) {
    adj[r][str[e.a]] = e;
  }
  int test_and_split(int s, int k, int p, char c) {
    if (k > p) return adj[s][c].to == -1 ? s : 0;
    edge e = adj[s][str[k]];
    if (c == str[e.a + p - k + 1]) return 0;
    int r = add_node();
    add_edge(r, {e.to, e.a + p - k + 1, e.b});
    add_edge(s, {r, e.a, e.a + p - k});
    return r;
  }
  int canonize(int s, int &k, int p) {
    if (p < k) return s;
    if (s == 0) { s = root; ++k; if (p < k) return s; }
    edge e = adj[s][str[k]];
    while (e.b - e.a <= p - k) {
      k = k + e.b - e.a + 1;
      s = e.to;
      if (k <= p) e = adj[s][str[k]];
    }
    return s;
  }
  int update(int s, int &k, int i) {
    int oldr = root, r = test_and_split(s, k, i-1, str[i]);
    while (r) {
      int rp = add_node();
      add_edge(r, {rp, i, n-1});
      if (oldr != root) suf[oldr] = r;
      oldr = r;
      s = canonize(suf[s], k, i-1);
      r = test_and_split(s, k, i-1, str[i]);
    }
    if (oldr != root) suf[oldr] = s;
    return s;
  }
  suffix_tree(const char *str) : str(str), n(strlen(str)) {
    add_node(); // bottom = 0
    root = add_node();
    int s = root, k = 0;
    for (int i = 0; i < n; ++i) {
      s = update(s, k, i);
      s = canonize(s, k, i);
    }
  }

  void display(int s, int tab = 0) {
    if (s == -1 || s == 0) return;
    for (int c = 0; c < 0x100; ++c) {
      if (adj[s][c].to != -1) {
        cout << string(tab, ' ');
        for (int i = adj[s][c].a; i <= adj[s][c].b; ++i)
          cout << str[i];
        cout << endl;
        display(adj[s][c].to, tab+2);
      }
    }
  }
  void display() { display(root); }
};


// === tick a time ===
#include <ctime>
double tick() {
  static clock_t oldtick;
  clock_t newtick = clock();
  double diff = 1.0*(newtick - oldtick) / CLOCKS_PER_SEC;
  oldtick = newtick;
  return diff;
}

int main() {

  //suffix_tree T("abracadabra");
  //suffix_tree T("abracadabra");
  //suffix_tree T("AAAAABCDEAAAAA#");
  const int n = 10000;
  char s[n+1];
  for (int i = 0; i < n; ++i)
    s[i] = "abc"[rand() % 3];

  tick();
  suffix_tree T(s);
  cout << "construction: " << tick() << endl;
  //T.display();

}

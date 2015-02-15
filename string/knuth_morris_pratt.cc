// 
// Knuth-Morris-Pratt string matching
//
// Description:
//   It processes a pattern string to find
//   all occurrence of a given text.
//
// Algorithm:
//   It constructs an automaton of accepting a pattern,
//   and then gives a text into the automaton..
//
// Complexity:
//   O(n + |occur|) with O(m) preprocessing.
//   In a random case, it reduced to O(n/m + |occur|).
//
// Verified:
//   SPOJ 21524
// 
// Comment:
//   KMP is often considered slower than BM.
//   However, in the programming contest setting, 
//   these are equally fast.
//
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

struct knuth_morris_pratt {
  int m;
  const char *p;
  vector<int> fail;
  knuth_morris_pratt(const char *p) : p(p), m(strlen(p)) {
    fail.resize(m+1, -1);
    for (int i = 1, j = -1; i <= m; ++i) {
      while (j >= 0 && p[j] != p[i-1]) j = fail[j];
      fail[i] = ++j;
    }
  }
  vector<int> match(const char *s) {
    int n = strlen(s);
    vector<int> occur;
    for (int i = 0, k = 0; i < n; ++i) {
      while (k >= 0 && s[i] != p[k]) k = fail[k];
      if (++k == m) { 
        /* match at s[i-m+1 ... i] */
        occur.push_back(i-m+1);
      }
    }
    return occur;
  }
};


// for comparison: boyer moore
struct boyer_moore {
  int m;
  const char *p;
  vector<int> skip, next;
  boyer_moore(const char *p) : p(p), m(strlen(p)) {
    skip.resize(0x100); // bad char heuristics
    for (int i = 0; i < m; ++i) 
      skip[p[i]] = m - i - 1;

    vector<int> g(m, m); // good suffix heuristics
    next.resize(m);
    for (int i = 0; i < m; ++i)
      next[i] = 2*m-i-1;
    for (int i = m-1, j = m; i >= 0; --i, --j) {
      g[i] = j;
      while (j < m && p[j] != p[i]) {
        next[j] = min(next[j], m-i-1);
        j = g[j];
      }
    }
  }
  vector<int> match(const char s[]) {
    int n = strlen(s);
    vector<int> occur;
    for (int i = m-1; i < n; ) {
      int j = m-1;
      while (j >= 0 && s[i] == p[j]) --i, --j;
      if (j < 0) { 
        /* match at s[i+1, ..., i+m] */
        occur.push_back(i+1);
        i += m + 1;
      } else i += max(skip[s[i]], next[j]);
    }
    return occur;
  }
};

int main() {
  int ncase; scanf("%d", &ncase);
  for (int icase = 0; icase < ncase; ++icase) {
    if (icase > 0) printf("\n");
    char s[1000010], p[1000010];
    scanf("%s %s", s, p);
    knuth_morris_pratt M(p);
    auto v = M.match(s);
    if (v.empty()) {
      printf("Not Found\n");
    } else {
      printf("%d\n", v.size());
      for (int i = 0; i < v.size(); ++i) {
        if (i > 0) printf(" ");
        printf("%d", v[i]+1);
      }
      printf("\n");
    }
  }
}


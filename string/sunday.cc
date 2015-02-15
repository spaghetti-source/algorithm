// 
// Sunday string matching
//
// Description:
//   It processes a pattern string to find
//   all occurrence of a given text.
//
// Algorithm:
//   It simplifies a boyer-moore algorithm.
//
// Complexity:
//   O(nm) in worst case; but fastest in a random case.
//
// Verified:
//   SPOJ 21524
// 
// Comment:
//   We recommend you to use KMP or BM in programming contest
//   because, there may be some "worst case" instances.
//   You can use this algorithm in more "practical" use.
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

struct sunday {
  int m;
  const char *p;
  vector<int> skip;
  sunday(const char *p) : p(p), m(strlen(p)) {
    skip.assign(0x100, m+1);
    for (int i = 0; i < m; ++i)
      skip[p[i]] = m - i;
  }
  vector<int> match(const char s[]) {
    int n = strlen(s);
    vector<int> occur;
    for (int i = 0; i <= n - m; ) {
      if (memcmp(p, s + i, m) == 0) {
        /* match at s[i, ..., i+m-1] */
        occur.push_back(i);
      }
      i += skip[s[i + m]];
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
    sunday M(p);
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


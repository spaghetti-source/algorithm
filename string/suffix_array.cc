//
// Suffix Array (Manbar and Myers' O(n (log n)^2))
//
// Description:
//   For a string s, tts suffix array is a lexicographically sorted
//   list of suffixes of s. For example, for s = "abbab", its SA is
//    0 ab
//    1 abbab
//    2 b
//    3 bab
//    4 bbab
// 
// Algorithm:
//   Manbar and Myers' doubling algorithm.
//   Suppose that suffixes are sorted by its first h characters.
//   Then, the comparison of first 2h characters is computed by
//     suf(i) <_2h suf(j) == if (suf(i) !=_h suf(j)) suf(i) <_h suf(j)
//                           else                    suf(i+h) <_h suf(j+h)
// 
// Complexity:
//   O(n (log n)^2). 
//   If we use radix sort instead of standard sort,
//   we obtain O(n log n) algorithm. However, it does not improve
//   practical performance so much.
//
// Verify:
//   SPOJ 6409: SARRAY (80 pt)
//
#include <iostream>
#include <vector>
#include <cstdio>
#include <cstring>
#include <algorithm>
#include <numeric>
#include <functional>

using namespace std;

#define fst first
#define snd second
#define all(c) ((c).begin()), ((c).end())

struct suffix_array {
  int n;
  vector<int> x;
  suffix_array(const char *s) : n(strlen(s)), x(n) { 
    vector<int> r(n), t(n); 
    for (int i = 0; i < n; ++i) r[x[i] = i] = s[i];
    for (int h = 1; t[n-1] != n-1; h *= 2) {
      auto cmp = [&](int i, int j) {
        if (r[i] != r[j]) return r[i] < r[j];
        return i+h < n && j+h < n ? r[i+h] < r[j+h] : i > j;
      };
      sort(all(x), cmp);
      for (int i = 0; i+1 < n; ++i) t[i+1] = t[i] + cmp(x[i], x[i+1]);
      for (int i = 0; i < n; ++i)   r[x[i]] = t[i];
    }
  }
  int operator[](int i) const { return x[i]; }
};

int main() {
  char s[100010];
  scanf("%s", s);
  suffix_array sary(s);
  for (int i = 0; i < sary.n; ++i) 
    printf("%d\n", sary[i]);
}

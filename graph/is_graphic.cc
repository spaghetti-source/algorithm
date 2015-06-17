//
// Grahpic sequence recognition
//
// Description
//   Erdos-Gallai theorem states that a sequence 
//   d = [d_1, ..., d_n] is a degree sequence of some
//   simple graph if and only if it satisfies
//     1. d_1+...+d_n is even,
//   and
//     2. d_1+...+d_k <= k(k-1)+min(k,d_k+1)+...+min(k, d_n).
//   for all k = 1, ..., n.
//
// Algorithm:
//   The right hand side of the 2nd condition is rewritten as
//     k(k-1)+k+...+k+d_p+1+...+d_n = k(p-1)+d_p+...+d_n.
//   where k > d_p >= ... >= d_n. Such p is obtained by a 
//   binary search; therefore it runs in O(log n) time.
//   Therefore the total complexity is O(n log n).
//
// Verified:
//   UVA 11414: Dream
//   UVA 10720: Graph Construction

#include <iostream>
#include <vector>
#include <cstdio>
#include <algorithm>
#include <functional>

using namespace std;

#define fst first
#define snd second
#define all(c) ((c).begin()), ((c).end())

bool is_graphic(vector<int> d) {
  int n = d.size();
  sort(all(d), greater<int>());
  vector<int> s(n+1);
  for (int i = 0; i < n; ++i) s[i+1] = s[i] + d[i];
  if (s[n] % 2) return false;
  for (int k = 1; k <= n; ++k) {
    int p = distance(d.begin(), 
        lower_bound(d.begin()+k, d.end(), k, greater<int>()));
    if (s[k] > k * (p-1) + s[n] - s[p]) return false;
  }
  return true;
}

// UVA 10720: Graph Construction
/*
int main() {
  int ncase; scanf("%d", &ncase);
  for (int icase = 0; icase < ncase; ++icase) {
    int n; scanf("%d", &n);
    vector<int> d(n);
    for (int i = 0; i < n; ++i)
      scanf("%d", &d[i]);
    printf("%s\n", (is_graphic(d) ? "Yes" : "No"));
  }
}
*/
int main() {
  for (int n; scanf("%d", &n); ) {
    if (n == 0) break;
    vector<int> d(n);
    for (int i = 0; i < n; ++i)
      scanf("%d", &d[i]);
    printf("%s\n", (is_graphic(d) ? "Possible" : "Not possible"));
  }
}

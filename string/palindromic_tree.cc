// 
// Palindromic Tree
//
// Description:
//   It is a tre data structure for a string s such that
//   1) each node corresponds to a palindromic substring,
//   2) next[p][c] is a link from "p" to "cpc"
//   3) suf[p] is a maximal palindromic suffix of p
// 
// Algorithm:
//   Manacher-like construction
//
// Complexity:
//   O(n)
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

struct palindromic_tree {
  vector<vector<int>> next;
  vector<int> suf, len;
  int new_node() {
    next.push_back(vector<int>(256,-1));
    suf.push_back(0);
    len.push_back(0);
    return next.size() - 1;
  }
  palindromic_tree(char *s) {
    len[new_node()] = -1;
    len[new_node()] =  0;
    int t = 1; 
    for (int i = 0; s[i]; ++i) {
      int p = t;
      for (; i-1-len[p] < 0 || s[i-1-len[p]] != s[i]; p = suf[p]);
      if ((t = next[p][s[i]]) >= 0) continue;
      t = new_node();
      len[t] = len[p] + 2;
      next[p][s[i]] = t;
      if (len[t] == 1) { 
        suf[t] = 1; // EMPTY
      } else {
        p = suf[p];
        for (; i-1-len[p] < 0 || s[i-1-len[p]] != s[i]; p = suf[p]);
        suf[t] = next[p][s[i]];
      }
    }
  }
  void display() {
    vector<char> buf;
    function<void (int)> rec = [&](int p) {
      if (len[p] > 0) {
        for (int i = buf.size()-1; i >=  0; --i) cout << buf[i];
        for (int i = len[p] % 2; i < buf.size(); ++i) cout << buf[i];
        cout << endl;
      }
      for (int a = 0; a < 256; ++a) {
        if (next[p][a] >= 0) {
          buf.push_back(a);
          rec(next[p][a]);
          buf.pop_back();
        }
      }
    };
    rec(0); rec(1);
  }
};

char s[30010];
int main() {
  int k; 
  scanf("%d %s", &k, s);
  palindromic_tree T(s);
  T.display();
}

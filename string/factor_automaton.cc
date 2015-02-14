// 
// Factor Automaton (aka. Directed Acyclig Word Graph)
//
// Description:
//   For a given string s, the factor automaton is an
//   automaton that accepts all substrings s[i,j).
//   The automaton has O(n) state.
//
// Algorithm:
//   Blumer et al's online construction.
//   See Section 6.3 of Crochemore-Rytter.
//
// Complexity:
//   O(n) time and space.
//
// Verified: 
//   SPOJ 22531
// 
// References:
//   A. Blumer, J. Blumer, D. Haussler, A. Ehrenfeucht, 
//   M. T. Chen, and J. Seiferas (1985):
//   The smallest automation recognizing the subwords of a text.
//   Theoretical Computer Science, vol. 40, pp. 31-55.
//
//   M. Crochemore, and W. Rytter (1994):
//   Text Algorithms.
//   Oxford University Press.
//
#include <iostream>
#include <vector>
#include <cstdio>
#include <queue>
#include <algorithm>
#include <functional>

using namespace std;

#define fst first
#define snd second
#define all(c) ((c).begin()), ((c).end())

struct factor_automaton {
  vector<vector<int>> link;
  vector<vector<bool>> bold;
  vector<int> suf;
  int root;
  int add_node() {
    link.push_back(vector<int>(0x100));
    bold.push_back(vector<bool>(0x100));
    suf.push_back(0);
    return link.size()-1;
  }
  factor_automaton(const char s[]) {
    add_node(); // = nil
    int sink = root = add_node();
    suf[root] = 0;
    for (int i = 0; s[i]; ++i) {
      char a = s[i];
      int newsink = add_node();
      link[sink][a] = newsink;
      bold[sink][a] = true;
      int w = suf[sink];
      while (w != 0 && link[w][a] == 0) {
        link[w][a] = newsink;
        bold[w][a] = false;
        w = suf[w];
      }
      int v = link[w][a];
      if (w == 0) suf[newsink] = root;
      else if (bold[w][a]) suf[newsink] = v;
      else {
        int newnode = add_node();
        link[newnode] = link[v];
        link[w][a] = newnode;
        bold[w][a] = true;
        suf[newsink] = newnode;
        suf[newnode] = suf[v]; 
        suf[v] = newnode;
        w = suf[w];
        while (w != 0 && link[w][a] == v && bold[w][a] == false) {
          link[w][a] = newnode;
          w = suf[w];
        }
      }
      sink = newsink;
    }
  }

  void disp() {
    cout << "--- display ---" << endl;
    for (int i = 1; i < link.size(); ++i) {
      cout << "  state " << i << endl;
      for (int c = 0; c < 0x100; ++c) {
        if (link[i][c] > 0) {
          if (bold[i][c]) cout << "    ==" << (char)c << "==> " << link[i][c] << endl;
          else            cout << "    --" << (char)c << "--> " << link[i][c] << endl;
        }
      }
    }

  }
};


int main() {
  for (char text[100000]; ~scanf("%s", text); ) {
    factor_automaton FA(text);

    vector<pair<int,char>> prev(FA.link.size());
    queue<int> que;
    que.push(FA.root);
    while (prev[0].fst == 0) {
      int s = que.front();
      que.pop();
      for (int c = 'A'; c <= 'Z'; ++c) {
      //for (int c = 'A'; c <= 'B'; ++c) { // for test
        if (prev[FA.link[s][c]].fst == 0) {
          prev[FA.link[s][c]] = {s, c};
          que.push(FA.link[s][c]);
          //cout << s << " --" << (char)c << "--> " << FA.link[s][c] << endl;
        }
      }
    }
    int s = 0;
    vector<char> result;
    do {
      result.push_back(prev[s].snd);
      s = prev[s].fst;
    } while (s > 1);
    for (int i = result.size()-1; i >= 0; --i)
      printf("%c", result[i]);
    printf("\n");
  }
}

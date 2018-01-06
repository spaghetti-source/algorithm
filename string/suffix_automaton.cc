//
// Suffix Automaton
//
// Description:
//
//   It is an automaton that accepts all suffixes
//   of a given string.
//
// Verified:
//   
//   SPOJ_SUBST1
//   SPOJ_SUBLEX
//
#include <bits/stdc++.h>

using namespace std;

#define fst first
#define snd second
#define all(c) ((c).begin()), ((c).end())
#define TEST(s) if (!(s)) { cout << __LINE__ << " " << #s << endl; exit(-1); }

struct SuffixAutomaton {
  vector<int> length = {0}, parent = {-1};
  vector<map<int,int>> next = {map<int,int>()};
  int size() const { return next.size(); }

  int extend(int p, int c) {
    int head = next.size();
    next.push_back(map<int,int>());
    length.push_back(length[p]+1);
    parent.push_back(0);
    while (p >= 0 && !next[p].count(c)) {
      next[p][c] = head;
      p = parent[p];
    }
    if (p >= 0) {
      int q = next[p][c];
      if (length[p]+1 == length[q]) {
        parent[head] = q;
      } else { 
        int clone = next.size(); // clone of q
        length.push_back(length[p]+1);
        parent.push_back(parent[q]); 
        next.push_back(next[q]);
        parent[q] = parent[head] = clone;
        while (p >= 0 && next[p][c] == q) {
          next[p][c] = clone;
          p = parent[p];
        }
      }
    }
    return head;
  }
  vector<int> topological_order;
  SuffixAutomaton(const char s[]) { // can be constructed online
    int p = 0;
    for(int i = 0; s[i]; ++i) 
      p = extend(p, s[i]);
    vector<int> mark(size()); // topological sort
    for (int i = 0; i < size(); ++i) 
      for (auto z: next[i]) ++mark[z.snd];
    topological_order.push_back(0);
    for (int i = 0; i < size(); ++i)
      for (auto z: next[i]) 
        if (--mark[z.snd] == 0) topological_order.push_back(z.snd);
  }
  template <class F> // topological order
  void process(F func) { 
    for (int i = 0; i < topological_order.size(); ++i) func(i);
  }
  template <class F>
  void processRev(F func) { // reverse topological order
    for (int i = topological_order.size()-1; i >= 0; --i) func(i);
  }
};
int countDistinctSubstrings(const char s[]) {
  SuffixAutomaton M(s);
  vector<int> dp(M.size());
  int ans = 0;
  dp[0] = 1;
  M.process([&](int i) {
    ans += dp[i];
    for (auto z: M.next[i])
      dp[z.snd] += dp[i];
  });
  return ans;
}

string kthSubstring(const char s[], int k) {
  SuffixAutomaton M(s);
  vector<int> num_terminal(M.size()); // preprocess that can be reused
  M.processRev([&](int i) { 
    num_terminal[i] = 1;
    for (auto z: M.next[i])
      num_terminal[i] += num_terminal[z.snd];
  });
  string ret; // answer to query
  for (int p = 0; k-- > 0; ) {
    for (auto z: M.next[p]) {
      if (num_terminal[z.snd] > k) {
        ret.push_back(z.fst);
        p = z.snd; break;
      } else k -= num_terminal[z.snd];
    }
  }
  return ret;
}

void SPOJ_SUBST1() {
  int ncase;
  scanf("%d", &ncase);
  for (int icase = 0; icase < ncase; ++icase) {
    char s[100000];
    scanf("%s", s);
    printf("%d\n", countDistinctSubstrings(s)-1);
  }
}
void SPOJ_SUBLEX() {
  char s[100000];
  scanf("%s", s);
  SuffixAutomaton M(s);

  vector<int> num_terminal(M.size());
  M.processRev([&](int i) { 
    num_terminal[i] = 1;
    for (auto z: M.next[i])
      num_terminal[i] += num_terminal[z.snd];
  });
  int ncase; scanf("%d", &ncase);
  for (int icase = 0; icase < ncase; ++icase) {
    int k; scanf("%d", &k);
    string ret;
    for (int p = 0; k-- > 0; ) {
      for (auto z: M.next[p]) {
        if (num_terminal[z.snd] > k) {
          ret.push_back(z.fst);
          p = z.snd; break;
        } else k -= num_terminal[z.snd];
      }
    }
    printf("%s\n", ret.c_str());
  }
}

int main() {
  //SPOJ_SUBST1();
  SPOJ_SUBLEX();
}

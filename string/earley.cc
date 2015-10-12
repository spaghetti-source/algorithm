//
// Earley Parser
//
// Description:
//   We are given CFG, i.e.,
//     A -> B
//     A -> aAa|bAb
//     B -> aa|bb
//     B -> a|b
//  It determines that a given string is matched by the CFG.
//
// Algorithm:
//   Earley algorithm. It generates all states with memoisation.
//   Here, state is given by (rule, pos-in-rule, pos-in-text).
//
// Complexity:
//   O(|G|^2 n^3) in the worst case.
//   If a grammar is simple, it usually reduced to O(|G|^2 n^2).
//
// Remark:
//   Because of simplicity, This implementation does not allow the 
//   epsilon rule. Please expand epsilon rule by hand.
//   (TODO!)
//   
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

struct earley_parser {
  vector<int> terminal;
  vector<vector<vector<int>>> grammar; 
  int add_symbol(char c = 0) {
    terminal.push_back(c);
    grammar.push_back({});
    return grammar.size()-1;
  }
  void add_grammar(int A, vector<int> As) {
    As.push_back(0);
    grammar[A].push_back(As);
  }
  earley_parser() { add_symbol(); add_symbol(); } 
  bool parse(const char s[], int init) {
    int n = strlen(s);
    struct state { int a, k, p, i; }; 
    vector<vector<vector<state>>> chart(n+1, vector<vector<state>>(grammar.size()));
    auto enqueue = [&](vector<state> &curr, const state &S) {
      for (auto &T: curr) 
        if (T.a == S.a && T.k == S.k && T.p == S.p && T.i == S.i) return;
      curr.push_back(S);
    };
    auto symbol = [&](const state &S) { return grammar[S.a][S.k][S.p]; };
    grammar[1] = { {init, 0} };
    vector<state> curr = {{1, 0, 0, 0}}, next;
    for (int k = 0; k <= n; ++k) {
      for (int i = 0; i < curr.size(); ++i) {
        state S = curr[i];
        int B = symbol(S);
        if (B) {
          if (!terminal[B]) {
            for (int j = 0; j < grammar[B].size(); ++j)
              enqueue(curr, {B, j, 0, k});
          } else if (terminal[B] == s[k]) {
            enqueue(next, {S.a, S.k, S.p+1, S.i});
          }
        } else {
          for (auto &T: chart[S.i][S.a]) 
            enqueue(curr, {T.a, T.k, T.p+1, T.i});
        }
      }
      for (auto &T: curr) 
        chart[k][symbol(T)].push_back(T);
      curr.swap(next);
      next.clear();
    }
    for (auto &T: chart[n][0]) 
      if (T.a == 1) return true;
    return false;
  }
};

int main() {
  earley_parser parser;
//     A -> B
//     A -> aAa|bAb
//     B -> aa|bb
//     B -> a|b
  int A = parser.add_symbol();
  int B = parser.add_symbol();
  int a = parser.add_symbol('a');
  int b = parser.add_symbol('b');
  parser.add_grammar(A, {B});
  parser.add_grammar(A, {a,A,a});
  parser.add_grammar(A, {b,A,b});
  parser.add_grammar(B, {a});
  parser.add_grammar(B, {b});
  parser.add_grammar(B, {a,a});
  parser.add_grammar(B, {b,b});
  for (char s[1024]; cin >> s; ) {
    cout << parser.parse(s, A) << endl;
  }
}

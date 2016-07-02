//
// Alpha-Beta Pruning with Transposition Table (Tic-Tac-Toe)
//
// Description:
//   Alpha-Beta pruning finds the game value of the state 
//   under the condition of "in [alpha, beta)".
//   Transposition table  memoises states in the search tree.
//
// Algorithm
//   The basic code for fail-soft alpha-beta pruning is given as
//
//   def alphabeta(s, alpha, beta):
//     if s is a terminate node:
//       return s.score
//     low = -inf
//     for each t: child of s
//       low = max(low, -alphabeta(t, -beta, -max(alpha,low)))
//       if low >= beta: break
//     return low
//   
//   Here, fail-soft means it always returns a realizable value.
//   This implementation is suitable for memoisation.
//
#include <iostream>
#include <vector>
#include <cstdio>
#include <algorithm>
#include <functional>
#include <unordered_map>

using namespace std;

#define fst first
#define snd second
#define all(c) ((c).begin()), ((c).end())
#define TEST(s) if (!(s)) { cout << __LINE__ << " " << #s << endl; exit(-1); }

const int w = 3, n = 9;
struct state {
  int p, T[n];

  state() { 
    p = +1;
    for (int i = 0; i < 9; ++i)
      T[i] = 0;
  }

  // must implement.
  // return the board status (finished or not)
  // and set board score into score variable
  // (player p win  ==> positive
  //           lose ==> negative )
  int score;
  bool finished() {
    const vector<vector<int>> lines = {
      {0,1,2},{3,4,5},{6,7,8},
      {0,3,6},{1,4,7},{2,5,8},
      {0,4,8},{2,4,6}
    };
    score = n;
    for (int i = 0; i < n; ++i) 
      if (T[i] != 0) --score;
    if (score == 0) return true;
    for (int i = 0; i < lines.size(); ++i) {
      bool same = true;
      for (int j = 1; j < lines[i].size(); ++j) 
        if (T[lines[i][0]] != T[lines[i][j]]) same = false;
      if (same && T[lines[i][0]]) {
        score += 100; if (T[lines[i][0]] != p) score *= -1;
        return true;
      }
    }
    return false;
  }
  void disp() {
    for (int i = 0; i < 3; ++i) {
      for (int j = 0; j < 3; ++j) {
        printf("%c ", T[3*i+j] > 0 ? 'o' : T[3*i+j] < 0 ? 'x' : '.');
      }
      printf("\n");
    }
  }
  // must implement
  bool operator == (const state &s) const {
    for (int i = 0; i < 9; ++i) 
      if (T[i] != s.T[i]) return false;
    return true;
  }
};
// must implement
namespace std {
  template <>
  struct hash<state> {
    size_t operator()(const state &s) const {
      size_t h = 0;
      for (int i = 0; i < 9; ++i) 
        h = 4 * h + (s.T[i] + 1);
      return h;
    }
  };
};

const int INF = 99999999;
unordered_map<state, int> cache;
pair<int, int> alphabeta(state s, int alpha, int beta) {
  int move = -1;
  if (s.finished()) return {s.score, move};
  int low = (cache.count(s) ? cache[s] : -INF);

  for (int k = 0; k < 9; ++k) {
    if (s.T[k]) continue;
    s.T[k] = s.p; s.p = -s.p;
    auto ans = alphabeta(s, -beta, -max(alpha,low));
    s.T[k] = 0; s.p = -s.p;
    if (low <= -ans.fst) {
      low = -ans.fst;
      move = k;
    }
    if (low >= beta) break;
  }
  return {cache[s] = low, move};
}


int main() {
  state s;
  alphabeta(s, -INF, INF);
  cout << cache.size() << endl;

  int player = +1;
  while (1) {
    s.disp();
    if (s.finished()) {
      if (s.score == 0) cout << "game is draw" << endl;
      else {
        if (s.p > 0) cout << "player +1 " << (s.score > 0 ? "win" : "lose") << endl;
        else         cout << "player -1 " << (s.score < 0 ? "lose" : "win") << endl;
      }
      break;
    }
    if (s.p == player) {
      int k; 
      while (1) {
        cin >> k;
        if (s.T[k] == 0) break;
      }
      s.T[k] = s.p; s.p = -s.p;
    } else {
      auto ans = alphabeta(s, -INF, INF);
      cout << ans.fst << endl;
      int k = ans.snd;
      s.T[k] = s.p; s.p = -s.p;
    }
  }
}

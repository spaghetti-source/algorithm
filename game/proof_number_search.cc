// 
// Proof Number Search
//
// Description:
//   Proof number search is search algorithm used in two-player endgame.
//   It is a best-first search that expands most provable nodes.
//
// Note:
//   It is slow if deciding win/lose is difficult until the playout.
//   (in this case, naive search may outperform) 
// 
//

#include <iostream>
#include <vector>
#include <cstdio>
#include <algorithm>
#include <unordered_map>
#include <unordered_set>
#include <functional>

using namespace std;

#define fst first
#define snd second
#define all(c) ((c).begin()), ((c).end())
#define TEST(s) if (!(s)) { cout << __LINE__ << " " << #s << endl; exit(-1); }

//  0  1  2  3
//  4  5  6  7
//  8  9 10 11
// 12 13 14 15

struct state {
  int p, T[16], status;
  state() { 
    p = +1;
    for (int i = 0; i < 16; ++i)
      T[i] = 0;
  }

  // must implement
  //   terminal: +1 (alice win) or -1 (bob win)
  //   non-terminal: 0
  int eval() {
    const vector<vector<int>> lines = {
      {0,1,2,3},{4,5,6,7},{8,9,10,11},{12,13,14,15},
      {0,4,8,12},{1,5,9,13},{2,6,10,14},{3,7,11,15},
      {0,5,10,15},{3,6,9,12}
      /*
      {0,1,2}, {1,2,3}, {4,5,6}, {5,6,7},
      {8,9,10}, {9,10,11}, {12,13,14}, {13,14,15},
      {0,4,8}, {4,8,12}, {1,5,9}, {5,9,13},
      {2,6,10}, {6,10,14}, {3,7,11}, {7,11,15},
      {0,5,10}, {1,6,11}, {2,5,8}, {3,6,9},
      {4,9,14}, {5,10,15}, {6,9,12}, {7,10,13}
      */
    };
    for (int i = 0; i < lines.size(); ++i) {
      bool same = true;
      for (int j = 1; j < lines[i].size(); ++j) 
        if (T[lines[i][0]] != T[lines[i][j]]) same = false;
      if (same && T[lines[i][0]]) return T[lines[i][0]];
    }
    int free = 0;
    for (int i = 0; i < 16; ++i)
      free += T[i] == 0;
    if (free == 0) return -1;
    return 0;
  }

  // must implement
  vector<state> next() const {
    vector<state> ss;
    for (int i = 0; i < 16; ++i) {
      if (T[i] == 0) {
        state t(*this);
        t.T[i] = t.p; t.p = -t.p;
        ss.push_back(t);
      }
    }
    return ss;
  }
  // must implement
  bool operator == (const state &s) const {
    for (int i = 0; i < 16; ++i) 
      if (T[i] != s.T[i]) return false;
    return true;
  }

  void disp() {
    for (int i = 0; i < 4; ++i) {
      for (int j = 0; j < 4; ++j) {
        printf("%c ", T[4*i+j] > 0 ? 'o' : T[4*i+j] < 0 ? 'x' : '.');
      }
      printf("\n");
    }
  }

};

// must implement
namespace std {
  template <>
  struct hash<state> {
    size_t operator()(const state &s) const {
      size_t h = 0;
      for (int i = 0; i < 16; ++i) 
        h = 4 * h + (s.T[i] + 1);
      return h;
    }
  };
};

int cached_search(state s) {
  unordered_map<state, int> cache;
  function<int(state&)> rec = [&](state &s) { 
    if (cache.count(s)) return cache[s];
    int eval = s.eval();
    if (eval) return eval;
    for (state t: s.next()) {
      if (rec(t) == s.p) return cache[s] = s.p;
    }
    return cache[s] = -s.p;
  };
  int ret = rec(s);
  cout << ret << " / opened nodes = " << cache.size() << endl;
  return 0;
}

const int INF = 99999999;
bool proof_number_search(state root) {
  vector<state> states;

  unordered_map<state, int> pn, dn;
  unordered_set<state> expanded;

  function<void(state)> augment = [&](state s) {
    // foward expand
    vector<state> next = s.next();
    if (!expanded.count(s)) {
      expanded.insert(s);
      for (state t: next) {
        if (pn.count(t)) continue;
        auto eval = t.eval();
        if (eval) {
          expanded.insert(t);
          if (eval > 0) { pn[t] = 0; dn[t] = INF; }
          else          { pn[t] = INF; dn[t] = 0; }
        } else {
          pn[t] = dn[t] = 1;
        }
      }
    } else {
      if (!next.empty()) {
        state t = next[0];
        for (int i = 1; i < next.size(); ++i) 
          if      (s.p == +1 && pn[t] > pn[next[i]]) t = next[i];
          else if (s.p == -1 && dn[t] > dn[next[i]]) t = next[i];
        augment(t);
      }
    }
    // backward update
    if (s.p == +1) {
      pn[s] = INF; dn[s] = 0;
      for (state t: next) {
        pn[s]  = min(pn[s], pn[t]);
        dn[s] += dn[t];
      }
    } else {
      pn[s] = 0; dn[s] = INF;
      for (state t: next) {
        pn[s] += pn[t];
        dn[s]  = min(dn[s], dn[t]);
      }
    }
  };
  pn[root] = dn[root] = 1;
  int iter = 0;
  while (pn[root] > 0 && dn[root] > 0) augment(root);
  //cout << "number of expanded nodes: " << expanded.size() << endl;
  return pn[root] == 0;
}

int main() {
  state s;
  cout << cached_search(s) << endl;
  cout << proof_number_search(s) << endl;
}

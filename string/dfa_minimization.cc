//
// Deterministic Finite Automata Minimization
//
// Description:
//
//   Hopcroft minimization algorithm. See Wikipedia.
//
// Complexity:
//
//   O(|A| n log n).
//
#include <bits/stdc++.h>

using namespace std;

#define fst first
#define snd second
#define all(c) ((c).begin()), ((c).end())
#define TEST(s) if (!(s)) { cout << __LINE__ << " " << #s << endl; exit(-1); }

struct Automaton {
  vector<vector<int>> trans;
  vector<bool> is_accept;
  int init = 0;
  int next(int state, int a) { return trans[state][a]; }
  bool accept(int state) { return is_accept[state]; }
  int size() { return trans.size(); }
};
template <class AM>
Automaton minimizeAutomaton(AM A) {
  // remove unreachables
  vector<int> seen(A.size());
  seen[A.init] = 1;
  vector<int> partition = {A.init}, pos, label;
  for (int i = 0; i < partition.size(); ++i) {
    pos.push_back(i);
    label.push_back(0);
    int state = partition[i];
    for (int a = 0; a <= 9; ++a) {
      int state_ = A.next(state, a);
      if (!seen[state_]) {
        seen[state_] = 1;
        partition.push_back(state_);
      }
    }
  }
  // make inverse mapping
  vector<int> inverse[partition.size()][10];
  for (int i = 0; i < partition.size(); ++i) {
    int state = partition[i];
    for (int a = 0; a <= 9; ++a) {
      inverse[A.next(state, a)][a].push_back(state);
    }
  }
  // Hopcroft minimization
  vector<int> begin = {0}, mid = {0}, end = {partition.size()};
  auto mark = [&](int state) {
    int x = label[state];
    if (pos[state] < mid[x]) return;
    int state_ = partition[mid[x]++];
    swap(pos[state], pos[state_]);
    swap(partition[pos[state]], partition[pos[state_]]);
  };
  auto refine = [&](int x) {
    if (mid[x] == begin[x]) return -1;
    if (mid[x] == end[x]) { mid[x] = begin[x]; return -1; }
    int y = begin.size();
    if (mid[x] - begin[x] < end[x] - mid[x]) {
      begin.push_back(begin[x]);
      end.push_back(mid[x]);
      mid.push_back(begin[x]);
      begin[x] = mid[x];
    } else {
      begin.push_back(mid[x]);
      end.push_back(end[x]);
      mid.push_back(mid[x]);
      end[x] = mid[x];
      mid[x] = begin[x];
    }
    for (int i = begin.back(); i < end.back(); ++i)
      label[partition[i]] = y;
    return y;
  };
  for (int state: partition) 
    if (A.accept(state)) mark(state);

  if (refine(0) >= 0) { // do Hopcroft minimization
    fill(all(seen), 0);
    seen[0] = seen[1] = 1;
    vector<int> process = {0, 1};
    vector<int> is_suspect(partition.size());
    while (!process.empty()) {
      int x = process.back(); process.pop_back();
      seen[x] = 0;
      for (int a = 0; a <= 9; ++a) {
        vector<int> suspect;
        for (int i = begin[x]; i < end[x]; ++i) {
          int u = partition[i];
          for (int state: inverse[u][a]) {
            int y = label[state];
            mark(state);
            if (!is_suspect[y]) {
              is_suspect[y] = 1;
              suspect.push_back(y);
            }
          }
        }
        for (int y: suspect) {
          is_suspect[y] = 0;
          int z = refine(y);
          if (z < 0) continue;
          if (seen[y]) {
            process.push_back(z);
            seen[z] = 1;
          } else {
            if (end[y] - begin[y] < end[z] - begin[z]) {
              process.push_back(y);
              seen[y] = 1;
            } else {
              process.push_back(z);
              seen[z] = 1;
            }
          }
        }
      }
    }
  }
  Automaton M; // completion
  M.trans.assign(begin.size(), vector<int>(10));
  M.is_accept.resize(begin.size());
  for (int x = 0; x < begin.size(); ++x) {
    int state = partition[begin[x]];
    M.is_accept[x] = A.accept(state);
    for (int a = 0; a <= 9; ++a) {
      int y = label[A.next(state, a)];
      M.trans[x][a] = label[A.next(state, a)];
    }
  }
  M.init = label[A.init];
  return M;
}


// state =  x : x == n % mod
struct ModuloAutomaton {
  int mod;
  ModuloAutomaton(int mod) : mod(mod) { }
  int init = 0;
  int size() { return mod; }
  int next(int state, int a) { return (10 * state + a) % mod; }
  bool accept(int state) { return state == 0; }
};
// state =  0        : empty
//          1        : fail
//          2 ... 10 : singleton and last number is state-1
//         11 ... 19 : increased and last number is state-10
//         20 ... 28 : decreased and last number is state-20
struct ZigZagAutomaton {
  int init = 0;
  int size() { return 29; }
  int next(int state, int a) {
    if (state == 0) return a == 0 ? 0 : a + 1;
    if (state == 1) return 1; 
    if (state <= 10) {
      int last = state - 1;
      if      (a > last) return a + 10;
      else if (a < last) return a + 20;
    } else if (state <= 19) {
      int last = state - 10;
      if (a < last) return a + 20;
    } else if (state <= 28) {
      int last = state - 20;
      if (a > last) return a + 10;
    }
    return 1;
  }
  bool accept(int state) { return state != 1; }
};
template <class Automaton1, class Automaton2>
Automaton intersectionAutomaton(Automaton1 A, Automaton2 B) {
  Automaton M;
  vector<vector<int>> table(A.size(), vector<int>(B.size(), -1));
  vector<int> x = {A.init}, y = {B.init};
  table[x[0]][y[0]] = 0;
  for (int i = 0; i < x.size(); ++i) {
    M.trans.push_back(vector<int>(10, -1));
    M.is_accept.push_back(A.accept(x[i]) && B.accept(y[i]));
    for (int a = 0; a <= 9; ++a) {
      int u = A.next(x[i], a), v = B.next(y[i], a);
      if (table[u][v] == -1) {
        table[u][v] = x.size();
        x.push_back(u);
        y.push_back(v);
      }
      M.trans[i][a] = table[u][v];
    }
  }
  return M;
}

template <class Automaton>
int digitDP(string num, Automaton A, int eq = 1) {
  int n = num.size();
  vector<vector<vector<int>>> dp(n+1);

  dp[0] = vector<vector<int>>(2, vector<int>(A.size()));
  dp[0][1][A.init] = 1;
  auto addTo = [&](int &x, int y) {
    if ((x += y) >= 10000) x -= 10000;
  };
  for (int i = 0; i < n; ++i) {
    dp[i+1] = vector<vector<int>>(2, vector<int>(A.size()));
    for (int tight = 0; tight <= 1; ++tight) {
      for (int state = 0; state < A.size(); ++state) {
        if (dp[i][tight][state] == 0) continue;
        int lim = (tight ? num[i] - '0' : 9);
        for (int d = 0; d <= lim; ++d) {
          int tight_ = tight && d == lim;
          int state_ = A.next(state, d);
          addTo(dp[i+1][tight_][state_], dp[i][tight][state]);
        }
      }
    }
    dp[i].clear();
  }
  int ans = 0;
  for (int tight = 0; tight <= eq; ++tight) 
    for (int state = 0; state < A.size(); ++state) 
      if (A.accept(state)) addTo(ans, dp[n][tight][state]);
  return ans;
}

template <class Automaton>
int debug(string num, Automaton A) {
  function<void(int,int,int,string)> rec 
    = [&](int i, int tight, int state, string s) {
    if (i == num.size()) {
      if (A.accept(state)) cout << s << endl;
      return;
    }
    int lim = (tight ? num[i] - '0' : 9);
    for (int d = 0; d <= lim; ++d) {
      int tight_ = tight && d == lim;
      int state_ = A.next(state, d);
      s.push_back('0' + d);
      rec(i+1, tight_, state_, s);
      s.pop_back();
    }
  };
  rec(0, 1, A.init, "");
}

void AOJ_ZIGZAG() {
  char A[1000], B[1000];
  int M;
  scanf("%s %s %d", A, B, &M);
  ZigZagAutomaton zigzag;
  ModuloAutomaton modulo(M);
  auto IM = minimizeAutomaton(intersectionAutomaton(zigzag, modulo));
  int a = digitDP(A, IM, 0);
  int b = digitDP(B, IM, 1);
  cout << (b + (10000 - a)) % 10000 << endl;
}
int main() { 
  AOJ_ZIGZAG();
  /*
  Automaton A;
  A.trans.assign(7, vector<int>(10, 6));
  A.trans[0][0] = 1;
  A.trans[0][1] = 2;
  A.trans[1][0] = 0;
  A.trans[1][1] = 3;
  A.trans[2][0] = 4;
  A.trans[2][1] = 5;
  A.trans[3][0] = 4;
  A.trans[3][1] = 5;
  A.trans[4][0] = 4;
  A.trans[4][1] = 5;
  A.trans[5][0] = 5;
  A.trans[5][1] = 5;
  A.is_accept.assign(7, 0);
  A.is_accept[2] = 1;
  A.is_accept[3] = 1;
  A.is_accept[4] = 1;
  auto B = minimizeAutomaton(A);
  cout << B.size() << endl;
  */
  //minimizeAutomaton(ZigZagAutomaton());
}

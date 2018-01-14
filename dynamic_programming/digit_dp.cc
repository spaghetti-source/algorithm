//
// Digit DP
//
// Description:
// 
//   Digit DP is a framework to solve problems of counting
//   the numbers less than equal to a given number and 
//   whose digits satisfy some constraint.
//
//   More generally, it can compute the sum-product
//     sum { prod(x) : 0 <= x <= z }
//   where
//     prod(x) = (((e * x[0]) * x[1])...) * x[n-1].
//
//   The sum operator + is required to be commutative, and 
//   right-distributive with respect to * as
//     (u + v) * d = (u * d + v * d)
//
//   The constraint of digits should be represented by a 
//   finite automaton that reads digits from left to right.
//   The DP has the table
//     dp[digit][tight][status]
//   with the following DP
//
//     dp[0][1][M.init] = e
//     for (int i = 0; i < n; ++i) {
//       for (int tight = 0; tight <= 1; ++tight) {
//         for (int state = 0; state < M.size(); ++state) {
//           int lim = tight ? 'z'-0 : 9;
//           for (int d = 0; d <= lim; ++d) {
//             oplusTo(dp[i+1][tight&&d==lim][next(state,d)], 
//               otimes(dp[i][tight][state], d));
//           }
//         }
//       }
//     }
//     for (int tight = 0; tight <= 1; ++tight)
//       for (int state = 0; state < M.size(); ++state)
//         if (M.accept(state)) oplusTo(ans, dp[n][tight][state];
//     return ans;
// 
//   Here, tight means the left digits are tight so that the current
//   digit cannot run over 0 to 9, and state means the state of the
//   automaton, which compresses the 10^len states to |Automaton| states.
//
// Verified:
//   SPOJ_CPCRC1C
//   AOJ_ZIGZAG
//   ABC_007D

#include <bits/stdc++.h>
using namespace std;

#define fst first
#define snd second
#define all(c) ((c).begin()), ((c).end())


// struct Value {
//   Value &operator+(Value y) 
//   Value &operator*(int d)
// };
// struct Automaton {
//   int init
//   int size()
//   int next(int state, int d)
//   bool accept(int state)
// };
template <class Value, class Automaton>
Value digitDP(string z, Value e, Automaton M, bool eq = 1) {
  struct Maybe {
    Value value;
    bool undefined = true;
  };
  auto oplusTo = [&](Maybe &x, Maybe y) {
    if (x.undefined) x = y;
    else if (!y.undefined) x.value += y.value;
  };
  auto otimes = [&](Maybe x, int d) {
    x.value *= d;
    return x;
  };
  int n = z.size();
  vector<vector<Maybe>> curr(2, vector<Maybe>(M.size()));
  curr[1][M.init] = {e, false};
  for (int i = 0; i < n; ++i) {
    vector<vector<Maybe>> next(2, vector<Maybe>(M.size()));
    for (int tight = 0; tight <= 1; ++tight) {
      for (int state = 0; state < M.size(); ++state) {
        if (curr[tight][state].undefined) continue;
        int lim = (tight ? z[i] - '0' : 9);
        for (int d = 0; d <= lim; ++d) {
          int tight_ = tight && d == lim;
          int state_ = M.next(state, d);
          oplusTo(next[tight_][state_], otimes(curr[tight][state], d));
        }
      }
    }
    curr = next;
  }
  Maybe ans;
  for (int tight = 0; tight <= eq; ++tight) 
    for (int state = 0; state < M.size(); ++state) 
      if (M.accept(state)) oplusTo(ans, curr[tight][state]);
  return ans.value;
}

template <class T>
string toString(T x) { 
  stringstream ss;
  ss << x;
  return ss.str();
}

// 
// Sum of digits.
// Since sum is not distributive, (u + v) + d != (u + d + v + d), 
// we need to augment the number to contain the number of numbers.
// Let + and * be defined by 
//   (u,a) + (v,b) = (u+v,a+b),
//   (u,a) * d = (u+a*d,a).
// Then, they are right-distributive as
//   ((u,a) + (v,b)) * c = (u+v,a+b) * c = (u+v+ac+bc,a+b),
//   (u,a) * c + (v,b) * c = (u+ac,a) + (v+bc,b) = (u+v+ac+bc,a+b).
//
using Int = long long;
Int sumOfDigits(string z, bool eq = true) {
  struct Value {
    Int count, sum;
    Value &operator+=(Value y) { count+=y.count; sum+=y.sum; return *this; }
    Value &operator*=(int d) { sum+=count*d; return *this; }
  };
  struct Automaton {
    int init = 0;
    int size() { return 1; }
    int next(int s, int d) { return 0; }
    int accept(int s) { return true; }
  };
  return digitDP(z, (Value){1,0}, Automaton(), eq).sum;
}
void SPOJ_CPCRC1C() {
  for (long long a, b; cin >> a >> b; ) {
    if (a < 0 && b < 0) break;
    cout << sumOfDigits(toString(b), true) 
          - sumOfDigits(toString(a), false) << endl;
  }
}

struct Automaton {
  vector<vector<int>> trans;
  vector<bool> is_accept;
  int init = 0;
  int next(int state, int a) { return trans[state][a]; }
  bool accept(int state) { return is_accept[state]; }
  int size() { return trans.size(); }
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

//
// Count the zigzag numbers that is a multiple of M.
// Here, a number is zigzag if its digits are alternatively 
// increasing and decreasing, like 14283415...
// Since there are multiple conditions, we use automaton 
// composition to simplify the approach.
//
void AOJ_ZIGZAG() {
  char A[1000], B[1000];
  int M;
  scanf("%s %s %d", A, B, &M);

  struct Value {
    int value = 0;
    Value &operator+=(Value x) {
      if ((value += x.value) >= 10000) value -= 10000;
      return *this;
    }
    Value &operator*=(int d) {
      return *this;
    }
  } e = (Value){1};

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
  } zigzag;

  // state =  x : x == n % mod
  struct ModuloAutomaton {
    int mod;
    ModuloAutomaton(int mod) : mod(mod) { }
    int init = 0;
    int size() { return mod; }
    int next(int state, int a) { return (10 * state + a) % mod; }
    bool accept(int state) { return state == 0; }
  } modulo(M);

  auto IM = intersectionAutomaton(zigzag, modulo);
  int a = digitDP(A, e, IM, 0).value;
  int b = digitDP(B, e, IM, 1).value;
  cout << (b + (10000 - a)) % 10000 << endl;
}

//
// Count the numbers that does not contain 4 and 7 in each digit.
//
void ABC007D() {
  string a, b;
  cin >> a >> b;

  struct ForbiddenNumber {
    int init = 0;
    int size() { return 2; }
    int next(int state, int a) { 
      if (state == 1) return 1;
      if (a == 4 || a == 9) return 1;
      return 0;
    }
    bool accept(int state) { return state == 1; }
  };
  struct Counter {
    long long value = 0;
    Counter &operator+=(Counter x) {
      value += x.value;
      return *this;
    }
    Counter &operator*=(int d) {
      return *this;
    }
  };
  cout << digitDP(b, (Counter){1}, ForbiddenNumber(), true).value
        - digitDP(a, (Counter){1}, ForbiddenNumber(), false).value << endl;
}

int main() {
  ABC007D();
  //SPOJ_CPCRC1C();
  //AOJ_ZIGZAG();
}

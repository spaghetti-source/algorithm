//
// Regex --> NFA --> DFA
//
// Description:
//   We are given a regex, eg, "a(b|c)d*e".
//   It construct NFA by recursive descent parsing,
//   and then construct DFA by powerset construction.
//
// Algorithm:
//   recursive descent parsing.
//   Robin-Scott's powerset construction.
//
// Complexity:
//   recursive descent parsiong: O(n) 
//   powerset construction: O(2^n) in worst case.
//
// Verified:
//   SPOJ 10354: Count Strings
//

#include <iostream>
#include <queue>
#include <vector>
#include <unordered_map> 
#include <unordered_set> 
#include <map>
#include <cstring>
#include <set>
#include <cstdio>
#include <bitset>

#include <algorithm>
#include <functional>

using namespace std;

#define fst first
#define snd second
#define all(c) ((c).begin()), ((c).end())
#define TEST(s) if (!s) { cout << __LINE__ << " " << #s << endl; exit(-1); }



const int NFA_STATE = 128, DFA_STATE = 1000, ALPHA = 256;
typedef bitset<NFA_STATE> subset;
struct NFA {
  static int size;
  static vector<int> next[NFA_STATE][ALPHA];
  static int new_node() {
    for (int a = 0; a < ALPHA; ++a) 
      next[size][a].clear();
    return size++;
  }
  static NFA symbol(char a) {
    int begin = new_node(), end = new_node();
    next[begin][a].push_back(end);
    return {begin, end};
  }
  static NFA unite(NFA x, NFA y) {
    int begin = new_node(), end = new_node();
    next[begin][0].push_back(x.begin);
    next[begin][0].push_back(y.begin);
    next[x.end][0].push_back(end);
    next[y.end][0].push_back(end);
    return {begin, end};
  }
  static NFA concat(NFA x, NFA y) {
    next[x.end][0].push_back(y.begin);
    return {x.begin, y.end};
  }
  static NFA star(NFA x) {
    int begin = new_node(), end = new_node();
    next[begin][0].push_back(x.begin);
    next[begin][0].push_back(end);
    next[x.end][0].push_back(x.begin);
    next[x.end][0].push_back(end);
    return {begin, end};
  }
  int begin, end;

  void closure(int u, subset &x) {
    x[u] = 1;
    for (int v: next[u][0]) 
      if (!x[v])
        closure(v, x);
  }
  bool run(const char *s) {
    subset x;
    closure(begin, x);
    for (; *s; ++s) {
      subset y;
      for (int u = 0; u < size; ++u) 
        if (x[u]) 
          for (int v: next[u][*s])
            closure(v, y);
      x = y;
    }
    return x[end];
  }
};
int NFA::size;
vector<int> NFA::next[NFA_STATE][ALPHA];

NFA parse(const char *s) {
  function<NFA ()> regex, factor, term;
  regex = [&]() {
    NFA a = factor();
    if (*s == '|') { ++s; a = NFA::unite(a, regex()); }
    return a;
  };
  factor = [&]() {
    NFA a = term();
    if (*s == '*') { a = NFA::star(a); ++s; }
    if (*s && *s != '|' && *s != ')') a = NFA::concat(a, factor());
    return a;
  };
  term = [&]() {
    if (*s == '(') { ++s; NFA a = regex(); ++s; return a; } 
    else { NFA a = NFA::symbol(*s); ++s; return a; }
  };
  return regex();
}


struct DFA {
  static int size, next[DFA_STATE][ALPHA];
  static int new_node() {
    memset(next[size], -1, sizeof(next[size]));
    return size++;
  }
  int begin, end[DFA_STATE];

  bool run(const char *s) {
    int u = begin;
    for (; *s; ++s) {
      u = next[u][*s];
      if (u < 0) return false;
    }
    return end[u];
  }
};
int DFA::size = 0, DFA::next[DFA_STATE][ALPHA];


DFA convert(NFA x) {
  DFA z;
  unordered_map<subset, int> states;
  vector<subset> process(1);
  x.closure(x.begin, process[0]);
  states[process[0]] = z.begin = z.new_node();
  while (!process.empty()) {
    auto S = process.back();
    process.pop_back();
    for (int a = 1; a < ALPHA; ++a) {
      subset T;
      for (int u = 0; u < x.size; ++u)
        if (S[u])
          for (int v: x.next[u][a])
            x.closure(v, T);
      if (T == 0) continue;
      if (!states.count(T)) {
        states[T] = z.new_node();
        z.end[states[T]] = T[x.end];
        process.push_back(T);
      }
      DFA::next[states[S]][a] = states[T];
    }
  }
  return z;
}


typedef long long ll;
void mul(int n, ll A[], ll B[], ll M) {
  ll C[n*n]; 
  for (int i = 0; i < n*n; ++i) C[i] = 0;
  for (int i = 0; i < n; ++i) 
    for (int k = 0; k < n; ++k)
      for (int j = 0; j < n; ++j) 
        C[n*i+j] = (C[n*i+j] + A[n*i+k] * B[n*k+j]) % M;
  for (int i = 0; i < n*n; ++i) A[i] = C[i];
}
void mulvec(int n, ll A[], ll x[], ll M) {
  ll y[n];
  for (int i = 0; i < n; ++i) {
    y[i] = 0;
    for (int j = 0; j < n; ++j) 
      y[i] = (y[i] + A[n*i+j] * x[j]) % M;
  }
  for (int i = 0; i < n; ++i) x[i] = y[i];
}
void powmul(int n, ll A[], ll k, ll x[], ll M) {
  for (; k > 0; k >>= 1) {
    if (k & 1) mulvec(n, A, x, M);
    mul(n, A, A, M);
  }
}

const ll M = 1000000007;
int main() {
  int ncase;
  scanf("%d", &ncase);
  for (int icase = 0; icase < ncase; ++icase) {
    DFA::size = NFA::size = 0;
    char s[1024];
    int l;
    scanf("%s %d", s, &l);
    DFA x = convert(parse(s));

    ll A[x.size*x.size];
    for (int i = 0; i < x.size*x.size; ++i) A[i] = 0;

    for (int i = 0; i < x.size; ++i) 
      for (int a = 1; a <= ALPHA; ++a) 
        if (x.next[i][a] >= 0)
          A[x.size*x.next[i][a]+i] = 1;

    ll a[x.size];
    for (int i = 0; i < x.size; ++i) a[i] =  0;
    a[0] = 1;
    powmul(x.size, A, l, a, M);
    ll ans = 0;
    for (int i = 0; i < x.size; ++i) if (x.end[i]) {
      ans += a[i];
      if (ans >= M) ans -= M;
    }
    printf("%lld\n", ans);
  }
}


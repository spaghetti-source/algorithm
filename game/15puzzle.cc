// 
// Fifteen Puzzle (IDA* with disjoint pattern database heuristics)
//

#include <iostream>
#include <vector>
#include <cstdio>
#include <queue>
#include <algorithm>
#include <unordered_map>
#include <functional>
#include <unordered_set>

using namespace std;

#define fst first
#define snd second
const int INF = 9999999;

// 16 is for blank tile, 0 is for abstracted tiles
//
//  1  2  3  4
//  5  6  7  8
//  9 10 11 12
// 13 14 15 16
//
struct state {
  int p, T[16]; 
  long long hash_value, pow[16];
  const long long P = 17, M = 1e9+7;
  state(vector<int> T_) { 
    copy(T_.begin(), T_.end(), T);
    for (p = 0; T[p] != 16; ++p);

    pow[0] = 1;
    for (int i = 0; i < 15; ++i) 
      pow[i+1] = (pow[i] * P) % M;
    hash_value = 0;
    for (int i = 0; i < 16; ++i) 
      hash_value = (hash_value + pow[i] * T[i]) % M;
  }

  bool move(int d) { // (d + 2) % 4 is reverse direction of d
    int q;
    if (d == 0) {
      if (p % 4 == 0) return false;
      q = p - 1;
    } else if (d == 1) {
      if (p < 4) return false;
      q = p - 4;
    } else if (d == 2) {
      if (p % 4 == 3) return false;
      q = p + 1;
    } else {
      if (p >= 12) return false;
      q = p + 4;
    } 
    hash_value -= (pow[p]*T[p] + pow[q]*T[q]);
    swap(T[p], T[q]);
    hash_value += (pow[p]*T[p] + pow[q]*T[q]);
    while (hash_value < 0) hash_value += M;
    while (hash_value >= M) hash_value -= M;
    p = q;
    return true;
  }

  bool operator == (const state &s) const {
    for (int i = 0; i < 16; ++i)
      if (T[i] != s.T[i]) return false;
    return true;
  }

  void disp() {
    for (int i = 0; i < 16; ++i) {
      printf(" %2d", T[i]);
      if (i % 4 == 3) printf("\n");
    }
  }
};
namespace std {
  template <>
  struct hash<state> { 
    size_t operator()(const state &s) const {
      return s.hash_value; 
    } 
  };
}
// 16 is for blank tile, 0 is for abstracted tiles
//
//  1  2  3  4
//  5  6  7  8
//  9 10 11 12
// 13 14 15 16
//

vector<vector<int>> pat = {
  /*
  {0,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,1},
  {0,0,0,0,0,0,1,1,1,1,1,0,0,0,0,0,1},
  {0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1},
  */
  /*
  {0,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,1},
  {0,0,0,0,1,1,1,1,0,0,0,0,0,0,0,0,1},
  {0,0,0,0,0,0,0,0,1,1,1,1,0,0,0,0,1},
  {0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1}
  */
  /*
  {0,1,1,0,0,1,0,0,0,0,0,0,0,0,0,0,1},
  {0,0,0,1,1,0,0,0,1,0,0,0,0,0,0,0,1},
  {0,0,0,0,0,0,1,1,0,0,1,0,0,0,0,0,1},
  {0,0,0,0,0,0,0,0,0,1,0,0,0,1,1,0,1},
  {0,0,0,0,0,0,0,0,0,0,0,1,1,0,0,1,1}
  */
  {0,1,1,0,0,1,1,0,0,1,0,0,0,0,0,0,1},
  {0,0,0,1,1,0,0,1,1,0,0,0,1,0,0,0,1},
  {0,0,0,0,0,0,0,0,0,0,1,1,0,1,1,1,1},
};
vector<unordered_map<state, int>> ph(pat.size());

void pattern_database(vector<int> pat, unordered_map<state, int> &ph) {
  fprintf(stderr, "construting pattern database\n");
  vector<int> goal = {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16};
  for (int i = 0; i < 16; ++i) 
    if (!pat[goal[i]]) goal[i] = 0;
  queue<state> que;
  que.push(goal);
  ph[goal] = 0;
  
  while (!que.empty()) {
    state s = que.front(); que.pop();
    int c = ph[s];
    for (int d = 0; d < 4; ++d) {
      int p = s.p;
      if (s.move(d)) {
        if (!ph.count(s)) {
          ph[s] = c + (pat[s.T[p]] != 0);
          que.push(s);
        }
        s.move((d+2)%4);
      }
    }
  }
  fprintf(stderr, "done constructing pattern database, size %zd\n", ph.size());
}


int manhattan_heuristics(state s) {
  // for comparison
  int d = 0;
  for (int k = 0; k < 16; ++k) {
    int a = s.T[k] - 1;
    if (s.T[k] != 16) d += abs(a%4 - k%4) + abs(a/4 - k/4);
  }
  return d;
}

int pattern_database_heuristics(state s) {
  int d = 0;
  for (int k = 0; k < pat.size(); ++k) {
    vector<int> T(s.T, s.T+16);
    for (int i = 0; i < 16; ++i)
      if (!pat[k][T[i]]) T[i] = 0;
    state t(T);
    d += ph[k][t];
  }
  return d;
}

int ida_star(state s) {
  int bound, opened = 0;
  function<bool(state&,int,int)> rec = [&](state &s, int depth, int prev) {
    ++opened;
    //int h = manhattan_heuristics(s);
    int h = pattern_database_heuristics(s);
    if (h == 0) return true;
    if (depth < h) return false;
    for (int d = 0; d < 4; ++d) {
      if ((d + 2) % 4 == prev) continue;
      if (s.move(d)) {
        if (rec(s, depth-1, d)) return true;
        s.move((d + 2) % 4);
      }
    }
    return false;
  };
  for (bound = 0; !rec(s, bound, 100); ++bound) {
    printf("%d\n", bound);
  }
  printf("%d\n", bound);
  printf("opened nodes = %d\n", opened);
  return 0;
}

// === tick a time ===
#include <ctime>
double tick() {
  static clock_t oldtick;
  clock_t newtick = clock();
  double diff = 1.0*(newtick - oldtick) / CLOCKS_PER_SEC;
  oldtick = newtick;
  return diff;
}

int main() {

  for (int k = 0; k < pat.size(); ++k)
    pattern_database(pat[k], ph[k]);

  vector<int> idx = {2,3,1,4,5,16,11,8,9,7,10,12,13,14,6,15};
  //vector<int> idx = {1,5,9,13,2,6,10,14,3,7,11,15,4,8,12,16}; // <- very difficult 

  state s(idx);
  tick();
  ida_star(s);
  // manhattan:         opened nodes = 101622
  // 3-4-4-4 pattern:   opened nodes = 9946 
  // 3-3-3-3-3 pattern: opened nodes = 7003
  // 5-5-5   pattern:   opened nodes = 1691
  printf("%f\n", tick());
}

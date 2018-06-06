//
// Segment Recognizer 
//
// Description:
//   Let M be an automaton and x be a sequence of alphabets.
//   The segment recognizer computes the transitioned state 
//   starting from s and reading x[i,j) in O(|M|) time.
//   The preprocessing requires O(|M| |x|) time and space.
//
//   The same method is implemented by the segment tree, 
//   where the time complexity is O(log n) and the space 
//   complexity is O(n log n). Thus, the segment recognizer 
//   is efficient if |M| is small.
//
// Algorithm:
//   Basically, it stores all the runs from all initial 
//   position i and initial state s. To reduce the space,
//   it merges two runs if they yields the same state.
//
// Reference
//   Mikola Bojanczyk (2009): "Factorization forests", 
//   International Conference on Developments in Language Theory,
//   pp. 1--17.
//
#include <bits/stdc++.h>

using namespace std;

#define fst first
#define snd second
#define all(c) ((c).begin()), ((c).end())
#define TEST(s) if (!(s)) { cout << __LINE__ << " " << #s << endl; exit(-1); }


// === tick a time ===
#include <ctime>
double tick() {
  static clock_t oldtick;
  clock_t newtick = clock();
  double diff = 1.0*(newtick - oldtick) / CLOCKS_PER_SEC;
  oldtick = newtick;
  return diff;
}

template <int MOD>
struct ModuloAutomaton {
  const int init = 0;
  int size() const { return MOD; }
  int next(int s, int d) const { return (s+d)%MOD; }
  int accept(int s) const { return s==0; }
};

// 0: free
// 1: selected
// 2: bottom
struct IndependenceAutomaton {
  const int init = 0;
  int size() const { return 3; }
  int next(int s, int d) const { 
    if (s == 0) return d;
    if (s == 1) return 2*d;
    if (s == 2) return s;
  }
  int accept(int s) const { return s!=2; }
};

template <class Automaton>
struct SegmentRecognizer {
  Automaton M;
  vector<int> x;

  struct Tape {
    int begin;
    vector<int> sequence;
  };
  vector<vector<int>> index;
  vector<Tape> tapes;

  SegmentRecognizer(Automaton M, vector<int> x) : M(M), x(x) { 
    index.assign(x.size()+1, vector<int>(M.size()));
    vector<int> stripe;
    for (int r = 0; r < M.size(); ++r) {
      stripe.push_back(r);
      index[0][r] = stripe[r];
      tapes.push_back({0, {r}});
    }
    for (int i = 0; i < x.size(); ++i) {
      unordered_set<int> available;
      for (int s = 0; s < M.size(); ++s)
        available.insert(s);
      vector<int> reallocate;
      for (int r = 0; r < M.size(); ++r) {
        int next = M.next(tapes[stripe[r]].sequence.back(), x[i]);
        if (available.count(next)) {
          available.erase(next);
          index[i+1][next] = stripe[r];
          tapes[stripe[r]].sequence.push_back(next);
        } else {
          reallocate.push_back(r);
        }
      }
      for (int r: reallocate) {
        int s = *available.begin();
        stripe[r] = tapes.size();
        index[i+1][s] = stripe[r];
        tapes.push_back({i+1, {s}});
        available.erase(s);
      }
    }
  }

  int getState(int i, int s, int j) {
    while (1) {
      auto &tape = tapes[index[i][s]];
      if (j - tape.begin < tape.sequence.size()) {
        return tape.sequence[j - tape.begin];
      } else {
        i = tape.begin + tape.sequence.size();
        s = M.next(tape.sequence.back(), x[i-1]);
      }
    }
  }
};
template <class Automaton>
SegmentRecognizer<Automaton> makeSegmentRecognizer(Automaton M, vector<int> s) {
  return SegmentRecognizer<Automaton>(M, s);
}

int main() {
  IndependenceAutomaton M;

  for (int n = 2; n < (1<<24); n*=2) {
    vector<int> x(n);
    for (int i = 0; i < n; ++i) {
      x[i] = (rand() % 10 == 0);
    }
    auto recognizer = makeSegmentRecognizer(M, x);

    tick();
    int count = 0;
    for (int iter = 0; iter < n; ++iter) {
      int v = (rand() % n) + 1;
      int u = rand() % v;
      count += recognizer.getState(u, 0, v);
    }
    double t = tick();
    cout << n << " " << t / n << endl;
  }
}

//
// Cube container
//
// Descrption:
//   It has a cube container which has six faces on
//   front, up, down, left, right, and bottom.
//   It admits the following three rotations
//     rotX: front -> up    -> back  -> down
//     rotY: left  -> up    -> right -> down
//     rotZ: left  -> front -> right -> down
//
// Algorithm:
//   Trivial.
//
// Verified:
//   SPOJ 21526
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

template <class T>
struct cube {
  T F, B, U, D, L, R;
  void rotX() { T x = D; D = B; B = U; U = F; F = x; } // FUBD -> DFUB
  void rotY() { T x = D; D = R; R = U; U = L; L = x; } // LURD -> DLUR
  void rotZ() { T x = B; B = R; R = F; F = L; L = x; } // LFRB -> BLFR

  // usage:
  // for (int i = 0; i < 24; ++i) {
  //   /* do something */
  //   c.next_roll();
  // }
  //
  // or 
  //
  // do {
  //   /* do something */
  // } while (c.next_roll());
  bool next_roll() { 
    static int it = 0;
    rotZ();
    if (it % 8 == 3) rotX();
    if (it % 8 == 7) rotY();
    return it = (it == 23 ? 0 : it+1);
  }
  bool operator==(cube c) const {
    for (int k = 0; k < 6; ++k) {
      if (U != c.U || F != c.F) continue;
      for (int i = 0; i < 4; ++i) {
        if (L == c.L && F == c.F && R == c.R && L == c.R) return true;
        c.rotZ();
      }
      if (k % 2) c.rotY(); else c.rotZ();
    }
    return false;
  }
};
template <class T>
ostream &operator<<(ostream &ofs, cube<T> c) {
  return (ofs << c.F << c.B << c.U << c.D << c.L << c.R);
}

void test() {
  cube<int> c = {1,2,3,4,5,6};
  do {
    cout << c << endl;
  } while (c.next_roll());
}

int main() {
  test();
  return 0;

  int ncase; scanf("%d", &ncase);
  for (int icase = 0; icase < ncase; ++icase) {
    cube<int> c;
    c.F = 0; c.U = 1; c.D = 2; c.L = 3; c.R = 4; c.B = 5;
    char s[6][1024];
    for (int i = 0; i < 6; ++i) 
      scanf("%s", s[i]);
    int q; scanf("%d", &q);
    for (int i = 0; i < q; ++i) {
      char t[1024];
      int k;
      scanf("%s %d", t, &k);
      k %= 4;
      if (t[0] == 'X') {
        while (k--) c.rotX();
      } else if (t[0] == 'Y') {
        while (k--) c.rotY();
      } else {
        while (k--) c.rotZ();
      }
    }
    printf("%s %s %s %s %s %s\n",
        s[c.F], s[c.U], s[c.D], s[c.L], s[c.R], s[c.B]);
  }
}

//
// Cube data structure
//
// Descrption:
//   A data structure for cube has six faces on
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
  T F, U, D, L, R, B;
  void rotX() { swap(D, B); swap(B, U); swap(U, F); } // FUBD -> DFUB
  void rotY() { swap(D, R); swap(R, U); swap(U, L); } // LURD -> DLUR
  void rotZ() { swap(B, R); swap(R, F); swap(F, L); } // LFRB -> BLFR
};

int main() {
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

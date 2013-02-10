//
// Karatsuba multiplication of polynomial
//
// Description:
// 
//   Given A(x) := a[0] + a[1] x + ... + a[an-1] x^{an-1} and
//         B(x) := b[0] + b[1] x + ... + b[bn-1] x^{bn-1},
//   compute C(x) = A(x) B(x).
//
//   In other words, compute the convolution of a and b:
//     c[k] = (a * b)[k] := sum_{i+j == k} a[i] b[j]
//
// Algorithm:
//
//   Karatsuba's divide-and-conquer.
//
//   Let n := an/2. Divide A, B as follows:
//     A(x) = Ah(x) x^n + Al(x),
//     B(x) = Bh(x) x^n + Bl(x).
//   Then C(x) = Ah(x) Bh(x) x^{2n} 
//             + (Ah(x) Bl(x) + Al(x) Bh(x)) x^n
//             + Al(x) Bl(x).
//   This formula requires 4 multiplications.
//
//   To reduce the number of multiplications, 
//   we introduce the auxilialy term:
//     D(x) = (Ah(x) + Al(x)) (Bh(x) + Bl(x)),
//   then we have
//     (Ah(x) Bl(x) + Al(x) Bh(x)) = D(x) - Ah(x) Al(x) - Bh(x) Bl(x).
//   This formula requires only 3 multiplications.
//
//
// Complexity
//
//   O(n^log 3).
//   In practive, it is faster than naive method for n >= 1000.
//
//
// Verified
//   SPOJ 31: Fast Multiplication
//   (Remark: TLE for VFMUL)
//
//
// References
//
// - A. Karatsuba and Y. Ofman (1962): 
//   Multiplication of Many-Digital Numbers by Automatic Computers,
//   Proceedings of the USSR Academy of Sciences, vol.145, pp.293-294

#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <vector>
#include <cstring>
#include <functional>
#include <algorithm>

using namespace std;

#define REP(i,n) for(int i=0;i<n;++i)

typedef long long LL;
void mul(LL a[], int an, LL b[], int bn, LL c[], int &cn) {
  if (an < bn) { swap(an, bn); swap(a, b); }
  if (bn == 0) { cn = 0; return; }
  if (bn == 1) { cn = an; REP(i, cn) c[i] = a[i] * b[0]; return; }
  cn = an + bn - 1;
  memset(c, 0, sizeof(c[0])*cn);

  int m = (an+1)/2, n = min(m, bn), tn1, tn2, tn3;
  LL *al = a, *ah = a+m, *bl = b, *bh = b+m;
  LL tmp1[2*m], tmp2[2*m], tmp3[2*m];
  REP(i, m) tmp1[i] = al[i] + (i < an-m ? ah[i] : 0);
  REP(i, n) tmp2[i] = bl[i] + (i < bn-n ? bh[i] : 0);
  mul(tmp1, m, tmp2, n, tmp3, tn3);     // = (al + ah)(bl + bh)
  mul(al, m, bl, n, tmp1, tn1);         // = al bl
  mul(ah, an-m, bh, bn-n, tmp2, tn2);   // = ah bh

  REP(i, tn1) { c[i] += tmp1[i]; c[i+m] -= tmp1[i]; }
  REP(i, tn2) { c[i+2*m] += tmp2[i]; c[i+m] -= tmp2[i]; }
  REP(i, tn3) c[i+m] += tmp3[i];
}


// a naive method
void mul0(LL a[], int an, LL b[], int bn, LL c[], int &cn) {
  cn = an + bn - 1;
  memset(c, 0, sizeof(c[0])*cn);
  REP(i, an) REP(j, bn) c[i+j] += a[i]*b[j];
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

const int N = 700;
LL a[N], b[N], c[2*N];
int main() {
  int an = N, bn = N-1, cn;
  REP(i, an) a[i] = rand();
  REP(i, bn) b[i] = rand();

  LL s;
  tick();
  s = 0;
  mul0(a, an, b, bn, c, cn);
  REP(i, cn) s += c[i];
  cout << s << " " << tick() << endl;

  s = 0;
  mul(a, an, b, bn, c, cn);
  REP(i, cn) s += c[i];
  cout << s << " " << tick() << endl;
}

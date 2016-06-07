// Karatsuba multiplication
//
// Given A(x) := a[0] + a[1] x + ... + a[an-1] x^{an-1} and
//       B(x) := b[0] + b[1] x + ... + b[bn-1] x^{bn-1},
// Compute C(x) = A(x) B(x)
//
// i.e., c[k] = (a * b)[k] := sum_{i+j == k} a[i] b[j]  (convolution)
//
//
// Complexity
//   O(n^log 3), always faster than naive method 
//
// Verified
//   SPOJ 31: Fast Multiplication
//   SPOJ 235: Very Fast Multiplication

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
  cn = an + bn - 1;
  memset(c, 0, sizeof(c[0])*cn);
  if (bn <= 32) {
    REP(i, an) REP(j, bn) c[i+j] += a[i]*b[j];
    return;
  }
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

/*
const int N = 20000;
LL a[N], b[N], c[2*N];
int main() {
  for (int n = 2; n < 10000; n *= 2) {
    printf("%d ", n);
    int an = n, bn = n, cn;
    REP(i, an) a[i] = rand();
    REP(i, bn) b[i] = rand();
    tick();
    mul0(a, an, b, bn, c, cn);
    printf("%f ", tick());
    LL s;
    s = 0;
    REP(i, cn) s += c[i];

    mul(a, an, b, bn, c, cn);
    printf("%f ", tick());
    REP(i, cn) s -= c[i];
    printf("%d\n", s);
    //cout << s << " " << tick() << endl;
  } 
}
*/

char s[300010];
int d = 6;
LL B = 1000000; // = 10^d
void read(LL a[], int &an) {
  scanf("%s", s);
  int n = strlen(s); 
  reverse(s, s+n);
  REP(i,2*d) s[n+i] = '0';

  an = 0;
  for (int i = 0; i < n; i += d) {
    a[an] = 0;
    for (int k = d-1; k >= 0; --k)
      a[an] = a[an] * 10 + s[i+k] - '0';
    ++an;
  }
}

const int N = 300010;
LL a[N], b[N], c[2*N];
int main() {
  int T; scanf("%d", &T);
  while (T--) {
    memset(a, 0, sizeof(a));
    memset(b, 0, sizeof(b));
    memset(c, 0, sizeof(c));
    int an, bn, cn;
    read(a, an);
    read(b, bn);
    mul(a, an, b, bn, c, cn);
    REP(i, cn) {
      if (c[i] >= B) {
        if (i == cn-1) ++cn;
        c[i+1] += c[i] / B;
        c[i] %= B;
      }
    }
    while (cn > 1 && c[cn-1] == 0) --cn;
    printf("%d", c[--cn]);
    while (cn > 0) printf("%0*d", d, c[--cn]);
    printf("\n");
  }
}

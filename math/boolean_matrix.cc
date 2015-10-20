// 
// Boolean matrix
//
// Description:
//   This admits very fast operations for boolean matrices.
//
// Algorithm:
//   Block matrix decomposition technique:
//   For a matrix A of size n x n, we split A as the block matrix
//   each block is of size n/W x n/W. Here, computation of each 
//   W x W block is performed by bit operations;
//
// Complexity: (in practice)
//   O(n^3); however, for practical n, it is basically
//   20--100 times faster than the naive implementation.
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

namespace bitmatrix {

typedef unsigned long long ull;
struct mat {
  int n, m;
  vector<vector<ull>> x;
  mat(int n, int m) : n(n), m(m), x(1+n/8, vector<ull>(1+m/8)) { }
  bool get(int i, int j) { 
    return x[i/8][j/8] & (1ull << (8*(i%8)+(j%8)));
  }
  void set(int i, int j, int b) {
    if (b) x[i/8][j/8] |=  (1ull << (8*(i%8)+(j%8))); 
    else   x[i/8][j/8] &= ~(1ull << (8*(i%8)+(j%8)));
  } 
};
void disp(ull a) {
  for (int i = 0; i < 8; ++i) {
    for (int j = 0; j < 8; ++j) {
      printf("%d", !!(a & 1));
      a >>= 1;
    }
    printf("\n");
  }
}
mat eye(int n) {
  mat I(n, n);
  for (int i = 0; i < I.x.size(); ++i) 
    I.x[i][i] = 0x8040201008040201;
  return I;
}
mat add(mat A, const mat &B) { 
  for (int i = 0; i < A.x.size(); ++i)
    for (int j = 0; j < A.x[0].size(); ++j) 
      A.x[i][j] |= B.x[i][j];
  return A;
}
// 64 x 64 matrix product in 160 operations. 
ull mul(ull a, ull b) { // C[i][j] |= A[i][k] & B[k][j]
  ull x, y, c = 0, u = 0x101010101010101, v = 0xff;
  x = a&u; u<<=1; x |= (x<<1); x |= (x<< 2); x |= (x << 4); 
  y = b&v; v<<=8; y |= (y<<8); y |= (y<<16); y |= (y <<32); c |= (x&y);
  x = a&u; u<<=1; x |= (x>>1); x |= (x<< 2); x |= (x << 4);
  y = b&v; v<<=8; y |= (y>>8); y |= (y<<16); y |= (y <<32); c |= (x&y);
  x = a&u; u<<=1; x |= (x<<1); x |= (x>> 2); x |= (x << 4);
  y = b&v; v<<=8; y |= (y<<8); y |= (y>>16); y |= (y <<32); c |= (x&y);
  x = a&u; u<<=1; x |= (x>>1); x |= (x>> 2); x |= (x << 4);
  y = b&v; v<<=8; y |= (y>>8); y |= (y>>16); y |= (y <<32); c |= (x&y);
  x = a&u; u<<=1; x |= (x<<1); x |= (x<< 2); x |= (x >> 4);
  y = b&v; v<<=8; y |= (y<<8); y |= (y<<16); y |= (y >>32); c |= (x&y);
  x = a&u; u<<=1; x |= (x>>1); x |= (x<< 2); x |= (x >> 4);
  y = b&v; v<<=8; y |= (y>>8); y |= (y<<16); y |= (y >>32); c |= (x&y);
  x = a&u; u<<=1; x |= (x<<1); x |= (x>> 2); x |= (x >> 4);
  y = b&v; v<<=8; y |= (y<<8); y |= (y>>16); y |= (y >>32); c |= (x&y);
  x = a&u; u<<=1; x |= (x>>1); x |= (x>> 2); x |= (x >> 4);
  y = b&v; v<<=8; y |= (y>>8); y |= (y>>16); y |= (y >>32); c |= (x&y);
  return c;
}
mat mul(mat A, mat B) {
  mat C(A.n, B.m);
  for (int i = 0; i < A.x.size(); ++i) 
    for (int k = 0; k < B.x.size(); ++k) 
      for (int j = 0; j < B.x[0].size(); ++j) 
        C.x[i][j] |= mul(A.x[i][k], B.x[k][j]);
  return C;
}
mat pow(mat A, int k) {
  mat X = eye(A.n);
  for (; k > 0; k >>= 1) {
    if (k & 1) X = mul(X, A);
    A = mul(A, A);
  }
  return X;
}
}

namespace vector_bool {
typedef vector<bool> vec;
typedef vector<vec> mat;
mat eye(int n) {
  mat I(n, vec(n));
  for (int i = 0; i < n; ++i)
    I[i][i] = 1;
  return I;
}
mat add(mat A, const mat &B) {
  for (int i = 0; i < A.size(); ++i)
    for (int j = 0; j < A[0].size(); ++j)
      A[i][j] = (A[i][j] | B[i][j]);
  return A;
}
mat mul(mat A, const mat &B) {
  for (int i = 0; i < A.size(); ++i) {
    vec x(A[0].size());
    for (int k = 0; k < B.size(); ++k) 
      for (int j = 0; j < B[0].size(); ++j) 
        x[j] = (x[j] | (A[i][k] & B[k][j]));
    A[i].swap(x);
  }
  return A;
}
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
  int n = 2000;
  {
    using namespace bitmatrix;
    tick();
    mat A(n, n);
    for (int i = 0; i < n; ++i) {
      for (int j = 0; j < n; ++j) {
        A.set(i, j, rand() % 2);
      }
    }
    A = pow(add(A, eye(n)), n); // (A + I)^n is a transitive closure
    printf("%f\n", tick());
  }

  {
    using namespace vector_bool;
    tick();
    mat A(n, vec(n));
    for (int i = 0; i < n; ++i) {
      for (int j = 0; j < n; ++j) {
        A[i][j] = rand() % 2;
      }
      A[i][i] = 1;
    }
    for (int k = 0; k < n; ++k)
      for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j)
          A[i][j] = (A[i][j] | (A[i][k] & A[k][j])) ;
    printf("%f\n", tick());
  }


  /*
  {
    tick();
    ull a = 0x123402010444107ull, b = 0x140455113142107ull;
    for (int iter = 0; iter < 10000000; ++iter) 
      a = mul_n(a, b);
    cout << tick() << endl;
  }
  */
}

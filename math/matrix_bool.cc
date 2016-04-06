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
//   50--60 times faster than the naive implementation.
//

#include <iostream>
#include <vector>
#include <cstdio>
#include <algorithm>
#include <functional>
#include <ctime>

using namespace std;

namespace bitmatrix {
typedef unsigned long long ull;
struct mat {
  int n, m;
  vector<vector<ull>> x;
  mat(int m, int n) : m(m), n(n), x(1+m/8, vector<ull>(1+n/8)) { }
  bool get(int i, int j) const { 
    return x[i/8][j/8] & (1ull << (8*(i%8)+(j%8)));
  }
  void set(int i, int j, int b) {
    if (b) x[i/8][j/8] |=  (1ull << (8*(i%8)+(j%8))); 
    else   x[i/8][j/8] &= ~(1ull << (8*(i%8)+(j%8)));
  } 
};
ostream &operator<<(ostream &os, const mat &A) {
  for (int i = 0; i < A.m; ++i) {
    for (int j = 0; j < A.n; ++j) 
      os << A.get(i, j);
    os << endl;
  }
  return os;
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

void disp(ull a) {
  for (int i = 0; i < 8; ++i) {
    for (int j = 0; j < 8; ++j) {
      printf("%d", !!(a & 1));
      a >>= 1;
    }
    printf("\n");
  }
}

ull mul(ull a, ull b) { // C[i][j] |= A[i][k] & B[k][j]
  const ull u = 0xff, v = 0x101010101010101;
  ull c = 0;
  for (;a && b; a >>= 1, b >>= 8)
    c |= (((a & v) * u) & ((b & u) * v));
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

ull transpose(ull a) {
  ull t = (a ^ (a >> 7)) & 0x00aa00aa00aa00aa;
  a = a ^ t ^ (t << 7);
  t = (a ^ (a >> 14)) & 0x0000cccc0000cccc;
  a = a ^ t ^ (t << 14);
  t = (a ^ (a >> 28)) & 0x00000000f0f0f0f0;
  a = a ^ t ^ (t << 28);
  return a;
}
mat transpose(mat A) {
  mat B(A.m, A.n);
  for (int i = 0; i < A.x.size(); ++i) 
    for (int j = 0; j <= A.x[0].size(); ++j) 
      B.x[j][i] = transpose(A.x[i][j]);
  return B;
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
mat pow(mat A, int k) {
  mat X = eye(A.size());
  for (; k > 0; k >>= 1) {
    if (k & 1) X = mul(X, A);
    A = mul(A, A);
  }
  return X;
}
}

int main() {
  for (int n = 1; n < 10000; n *= 2) {
    printf("%d\t", n);
    {
      using namespace bitmatrix;
      mat A(n, n);
      for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
          A.set(i, j, rand() % 2);
        }
      }
      double _time = clock();
      A = pow(add(A, eye(n)), n); // (A + I)^n is a transitive closure
      _time = clock() - _time;
      printf("%f\t", _time / CLOCKS_PER_SEC);
    }
    printf("\n");
    continue;
    {
      using namespace vector_bool;
      mat A(n, vec(n));
      for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
          A[i][j] = rand() % 2;
        }
        A[i][i] = 1;
      }
      double _time = clock();
      A = pow(add(A, eye(n)), n); // (A + I)^n is a transitive closure
      _time = clock() - _time;
      printf("%f\n", _time / CLOCKS_PER_SEC);
    }
  }
}

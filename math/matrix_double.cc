//
// Matrix Computation Algorithms (double)
//

#include <iostream>
#include <vector>
#include <cmath>
#include <cstdio>
#include <algorithm>
#include <functional>

using namespace std;

#define fst first
#define snd second
#define all(c) ((c).begin()), ((c).end())

template <class T>
ostream &operator<<(ostream &os, const vector<T> &v) {
  os << "[";
  for (int i = 0; i < v.size(); os << v[i++]) 
    if (i > 0) os << " ";
  os << "]";
  return os;
}
template <class T>
ostream &operator<<(ostream &os, const vector<vector<T>> &v) {
  os << "[";
  for (int i = 0; i < v.size(); os << v[i++]) 
    if (i > 0) os << endl << " ";
  os << "]";
  return os;
}


typedef vector<double> vec;
typedef vector<vec> mat;
int sign(double x) { return x < 0 ? -1 : 1; }
mat eye(int n) {
  mat I(n, vec(n));
  for (int i = 0; i < n; ++i)
    I[i][i] = 1;
  return I;
}
mat add(mat A, const mat &B) {
  for (int i = 0; i < A.size(); ++i)
    for (int j = 0; j < A[0].size(); ++j)
      A[i][j] += B[i][j];
  return A;
}
mat mul(mat A, const mat &B) {
  for (int i = 0; i < A.size(); ++i) {
    vec x(A[0].size());
    for (int k = 0; k < B.size(); ++k) 
      for (int j = 0; j < B[0].size(); ++j) 
        x[j] += A[i][k] * B[k][j];
    A[i].swap(x);
  }
  return A;
}
mat pow(mat A, int k) {
  mat X = eye(A.size());
  for (; k > 0; k /= 2) {
    if (k & 1) X = mul(X, A);
    A = mul(A, A);
  }
  return X;
}
double diff(vec a, vec b) {
  double S = 0;
  for (int i = 0; i < a.size(); ++i)
    S += (a[i] - b[i]) * (a[i] - b[i]);
  return sqrt(S);
}
double diff(mat A, mat B) {
  double S = 0;
  for (int i = 0; i < A.size(); ++i)
    for (int j = 0; j < A[0].size(); ++j)
      S += (A[i][j]-B[i][j])*(A[i][j]-B[i][j]);
  return sqrt(S);
}
vec mul(mat A, vec b) {
  vec x(A.size());
  for (int i = 0; i < A.size(); ++i)
    for (int j = 0; j < A[0].size(); ++j)
      x[i] += A[i][j] * b[j];
  return x;
}
mat transpose(mat A) {
  for (int i = 0; i < A.size(); ++i)
    for (int j = 0; j < i; ++j)
      swap(A[i][j], A[j][i]);
  return A;
}
double det(mat A) {
  double D = 1;
  for (int i = 0; i < A.size(); ++i) { 
    int p = i;
    for (int j = i+1; j < A.size(); ++j) 
      if (fabs(A[p][i]) < fabs(A[j][i])) p = j;
    swap(A[p], A[i]); 
    for (int j = i+1; j < A.size(); ++j) 
      for (int k = i+1; k < A.size(); ++k) 
        A[j][k] -= A[i][k] * A[j][i] / A[i][i];
    D *= A[i][i];
    if (p != i) D = -D;
  }
  return D;
}
vec solve(mat A, vec b) { // assume: A is non-singular
  for (int i = 0; i < A.size(); ++i) {
    int p = i;
    for (int j = i+1; j < A.size(); ++j) 
      if (fabs(A[p][i]) < fabs(A[j][i])) p = j;
    swap(A[p], A[i]); swap(b[p], b[i]); 
    for (int j = i+1; j < A.size(); ++j) { 
      for (int k = i+1; k < A.size(); ++k) 
        A[j][k] -= A[i][k] * A[j][i] / A[i][i];
      b[j] -= b[i] * A[j][i] / A[i][i];
    }
  }
  for (int i = A.size()-1; i >= 0; --i) {
    for (int j = i+1; j < A.size(); ++j)
      b[i] -= A[i][j] * b[j];
    b[i] /= A[i][i];
  }
  return b;
}

// TODO: verify
mat solve(mat A, mat B) { // A^{-1} B
  for (int i = 0; i < A.size(); ++i) { // forward elimination
    int p = i;
    for (int j = i+1; j < A.size(); ++j) 
      if (fabs(A[p][i]) < fabs(A[j][i])) p = j;
    swap(A[p], A[i]); swap(B[p], B[i]);
    for (int j = i+1; j < A.size(); ++j) { 
      double coef = A[j][i] / A[i][i];
      for (int k = i; k < A.size(); ++k) 
        A[j][k] -= A[i][k] * coef;
      for (int k = 0; k < B[0].size(); ++k) 
        B[j][k] -= B[i][k] * coef;
    }
  }
  for (int i = A.size()-1; i >= 0; --i) { // backward substitution
    for (int j = i+1; j < A.size(); ++j)
      for (int k = 0; k < 0; ++k) 
        B[i][k] -= A[i][j] * B[j][k];
    for (int k = 0; k < B[0].size(); ++k) 
      B[i][k] /= A[i][i];
  }
  return B;
}

// LU factorization
struct lu_data { mat A; vector<int> pi; };
lu_data lu(mat A) {
  vector<int> pi;
  for (int i = 0; i < A.size(); ++i) {
    int p = i;
    for (int j = i+1; j < A.size(); ++j) 
      if (fabs(A[p][i]) < fabs(A[j][i])) p = j;
    pi.push_back(p);
    swap(A[p], A[i]);
    for (int j = i+1; j < A.size(); ++j) {
      for (int k = i+1; k < A.size(); ++k)
        A[j][k] -= A[i][k] * A[j][i] / A[i][i];
      A[j][i] /= A[i][i];
    }
  }
  return {A, pi};
}
vec solve(lu_data LU, vec b) {
  mat &A = LU.A;
  vector<int> &pi = LU.pi;
  for (int i = 0; i < pi.size(); ++i)
    swap(b[i], b[pi[i]]);
  for (int i = 0; i < A.size(); ++i)
    for (int j = 0; j < i; ++j)
      b[i] -= A[i][j] * b[j];
  for (int i = A.size()-1; i >= 0; --i) {
    for (int j = i+1; j < A.size(); ++j)
      b[i] -= A[i][j] * b[j];
    b[i] /= A[i][i];
  }
  return b;
}

// QR factorization (when will we use?)
const double EPS = 1e-8;
struct qr_data { mat Q, R; };
qr_data qr(mat A) {
  const int n = A.size();
  mat Q = eye(A.size());
  for (int i = 0; i < A.size(); ++i) {
    for (int j = i+1; j < A.size(); ++j) {
      double c = A[i][i], s = A[j][i], r = sqrt(c*c + s*s);
      if (r <= EPS) continue;
      c /= r; s /= r;
      auto rot = [&](double &x, double &y) {
        tie(x, y) = make_pair(c*x+s*y, -s*x+c*y);
      };
      for (int k = i; k < A.size(); ++k) rot(A[i][k], A[j][k]);
      for (int k = 0; k < A.size(); ++k) rot(Q[k][i], Q[k][j]);
    }
  }
  return {Q, A};
}
vec solve(qr_data QR, vec b) {
  mat &Q = QR.Q, &R = QR.R; 
  vec x(Q.size());
  for (int i = 0; i < Q.size(); ++i)
    for (int j = 0; j < Q.size(); ++j)
      x[i] += Q[j][i] * b[j];
  for (int i = Q.size()-1; i >= 0; --i) {
    for (int j = i+1; j < Q.size(); ++j)
      x[i] -= R[i][j] * x[j];
    x[i] /= R[i][i];
  }
  return x;
}

// Cholesky factorization
//
// A = L L' for positive definite A
// (upper part of L is garbage)
struct cholesky_data { mat L; };
cholesky_data chol(mat A) {
  for (int i = 0; i < A.size(); ++i) {
    for (int j = i; j < A.size(); ++j) {
      for (int k = i-1; k >= 0; --k) 
        A[j][i] -= A[i][k] * A[j][k];
      if (i == j) A[i][i] = sqrt(A[i][i]);
      else        A[j][i] /= A[i][i];
    }
  }
  return {A};
}
vec solve(cholesky_data C, vec b) {
  mat &L = C.L;
  for (int i = 0; i < L.size(); ++i) {
    for (int j = 0; j < i; ++j) 
      b[i] -= L[i][j] * b[j];
    b[i] /= L[i][i];
  }
  for (int i = L.size()-1; i >= 0; --i) {
    for (int j = i+1; j < L.size(); ++j)
      b[i] -= L[j][i] * b[j];
    b[i] /= L[i][i];
  }
  return b;
}

// Eigenvalues of symmetric matrix
//
// Perform Hessenberg QR for tridiagonalization,
// and then perform implicit LR for diagonalization.
//
// O(n^3 + T n^2); where T is the number of iterations.
//
struct tridiag_data { vec d, e; };
tridiag_data tridiag(mat A) {
  vec u(A.size()), v(A.size());
  for(int k = 1; k < A.size(); ++k) {
    double norm = 0;
    for (int i = k; i < A.size(); ++i) {
      u[i] = A[k-1][i];
      norm += u[i] * u[i];
    }
    u[k] += sign(u[k]) * sqrt(norm);
    norm += u[k]*u[k] - A[k-1][k]*A[k-1][k];

    for (int i = 0; i < A.size(); ++i) {
      v[i] = 0;
      for (int j = k; j < A.size(); ++j)
        v[i] += A[i][j] * u[j];
      v[i] *= 2 / norm;
    }
    double gamma = 0;
    for (int j = k; j < A.size(); ++j)
      gamma += u[j] * v[j];
    gamma /= norm;
    for (int j = k; j < A.size(); ++j)
      v[j] -= gamma * u[j];
    for (int i = 0; i < A.size(); ++i) {
      for (int j = k; j < A.size(); ++j) {
        A[i][j] -= v[i] * u[j];
        A[j][i] -= v[i] * u[j];
      }
    }
  }
  for (int i = 0; i   < A.size(); ++i) u[i] = A[i][i];
  for (int i = 0; i+1 < A.size(); ++i) v[i] = A[i][i+1];
  return {u, v};
}
vec trieig(tridiag_data tri) {
  vec &d = tri.d, &e = tri.e;
  for (int i, m, l = 0; l < d.size(); ) {
    for (m = l; m < d.size()-1; ++m)
      if (fabs(e[m]) < EPS) break;
    if (m == l) { ++l; continue; }
    double s = 1, c = 1, p = 0;
    double g = (d[l+1]-d[l])/(2*e[l]);
    if (g > 0) g = d[m]-d[l] + e[l]/(g + sign(g) * sqrt(g*g+1));
    for (i = m-1; i >= l; --i) {
      double f = s*e[i], b = c*e[i], r = sqrt(f*f+g*g);
      if ((e[i+1] = r) < EPS) break;
      s = f/r; c = g/r;
      r = (d[i]-d[i+1]+p)*s + 2*c*b;
      double dp = s*r - p;
      d[i+1] += dp; p += dp; g = c*r - b;
    }
    d[i+1] -= p; e[i+1] = g; e[m] = 0;
  }
  return d;
}
vec symeig(mat A) { return trieig(tridiag(A)); }



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
  int n = 10;
  mat A(n, vec(n));
  for (int i = 0; i < n; ++i)
    for (int j = 0; j < n; ++j)
      A[i][j] = A[j][i] = rand() % 10;

  for (int i = 0; i < n; ++i) {
    A[i][i] = 1;
    for (int j = 0; j < n; ++j) 
      if (i != j) 
        A[i][i] += fabs(A[i][j]);
  }
  
  vec b(n);
  for (int i = 0; i < n; ++i)
    b[i] = rand() % 100;

  tick();
  vec x = solve(A, b);
  //cout << x << endl;
  cout << diff(b, mul(A, x)) << endl;
  cout << tick() << endl;

  auto LU = lu(A);
  vec y = solve(LU, b);
  //cout << y << endl;
  cout << diff(b, mul(A, y)) << endl;

  tick();
  auto QR = qr(A);
  vec z = solve(QR, b);
  //cout << z << endl;
  cout << diff(b, mul(A, z)) << endl;

  tick();
  auto C = chol(A);
  vec w = solve(C, b);
  cout << tick() << endl;
  cout << diff(b, mul(A, w)) << endl;

  //cout << symeig(A) << endl;
}

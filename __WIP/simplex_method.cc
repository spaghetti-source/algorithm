//
// Simplex Method for Linear Programming
//
// 
//
#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <vector>
#include <cmath>
#include <cstring>
#include <functional>
#include <algorithm>

using namespace std;

#define ALL(c) c.begin(), c.end()
#define FOR(i,c) for(typeof(c.begin())i=c.begin();i!=c.end();++i)
#define REP(i,n) for(int i=0;i<n;++i)
#define fst first
#define snd second

const double EPS = 1e-6, INF = 1./0.;

#define LT(a,b) (a+EPS < b)
#define LE(a,b) (a-EPS < b)

//
// Primal Simplex Method
//
// min. c x 
// s.t. A x <= b
//        x >= 0
// (for b >= 0)
//
// <=> 
//
// min. c x 
// s.t. A x + s == b
//         x, s >= 0
//
// There is a primal feasible solution (0, b).
// (therefore this LP is always feasible)
//
//
// n variables, m constraints
// (c: n-vector, b: m-bector, A: nxm-array)
//
// Remark:
//   You can apply this for c >= 0 by taking LP-dual.
//
double simplexMethodP(double c[], int n, double b[], int m, double A[]) {
  double T[m+1][n+m+1];
  memset(T, 0, sizeof(T));
  REP(j, m) {
    REP(i, n) T[j][i] = A[j*n+i];
    T[j][n+j] = 1;
    T[j][n+m] = b[j];
  }
  REP(i, n) T[m][i] = c[i];
  while (1) {
    int p = 0, q = 0;
    REP(i, n+m) if (T[m][i] <= T[m][p]) p = i;
    if (LE(0, T[m][p])) return T[m][n+m];
    REP(j, m) if (LE(0, T[j][p]))
      if (T[j][p]*fabs(T[q][n+m]) >= T[q][p]*fabs(T[j][n+m])) q = j;
    if (LE(T[q][p], 0)) return INF; // primal infeasible
    REP(j, m+1) if (j != q) {
      double alpha = T[j][p] / T[q][p];
      REP(i, n+m+1) T[j][i] -= T[q][i] * alpha;
    }
  }
}

// min. c x
// s.t. A x <= b
//        x >= 0
//
// assume c >= 0
//
double simplexMethodD(double c[], int n, double b[], int m, double A[]) {
  double T[m+1][n+m+1];
  memset(T, 0, sizeof(T));
  REP(j, m) {
    REP(i, n) T[j][i] = A[j*n+i];
    T[j][n+j] = 1;
    T[j][n+m] = b[j];
  }
  REP(i, n) T[m][i] = c[i];
  while (1) {
    int p = 0, q = 0;
    REP(j, m) if (T[j][n+m] <= T[q][n+m]) q = j;
    if (LE(0, T[q][n+m])) return -T[m][n+m];
    REP(i, n+m+1) T[q][i] *= -1;
    REP(i, n+m) if (LE(0, T[q][i]))
      if (T[q][i]*fabs(T[m][p]) >= T[q][p]*fabs(T[m][i])) p = i;
    if (LE(T[q][p], 0)) return -INF; // dual infeasible
    REP(j, m+1) if (j != q) {
      double alpha = T[j][p] / T[q][p];
      REP(i, n+m+1) T[j][i] -= T[q][i] * alpha;
    }
  }
}

//
// Parametric Self-Dual simplex method
//
// min. c x
// s.t. A x <= b
//        x >= 0
//
// for general c and b.
//   
double simplexMethodPD(double c[], int n, double b[], int m, double A[]) {
  double T[m+1][n+m+1];
  memset(T, 0, sizeof(T));
  REP(j, m) {
    REP(i, n) T[j][i] = A[j*n+i];
    T[j][n+j] = 1;
    T[j][n+m] = b[j];
  }
  REP(i, n) T[m][i] = c[i];
  while (1) { 
    int p = 0, q = 0;
    REP(i, n+m) if (T[m][i] <= T[m][p]) p = i;
    REP(j, m) if (T[j][n+m] <= T[q][n+m]) q = j;
    double t = min(T[m][p], T[q][n+m]);
    if (t >= -EPS) return -T[m][n+m]; // optimal

    if (t < T[q][n+m]) { // tight on c -> primal update
      REP(j, m) if (T[j][p] >= EPS) 
        if (T[j][p]*fabs(T[q][n+m]-t) >= T[q][p]*fabs(T[j][n+m]-t)) q = j;
      if (T[q][p] <= EPS) return INF; // primal infeasible
    } else { // tight on b -> dual update
      REP(i, n+m+1) T[q][i] *= -1;
      REP(i, n+m) if (T[q][i] >= EPS) 
        if (T[q][i]*fabs(T[m][p]-t) >= T[q][p]*fabs(T[m][i]-t)) p = i;
      if (T[q][p] <= EPS) return -INF; // dual infeasible
    }
    REP(j, m+1) if (j != q) T[j][p] /= T[q][p]; T[q][p] = 1; // pivot(q,p)
    REP(j, m+1) if (j != q) {
      double alpha = T[j][p] / T[q][p];
      REP(i, n+m+1) T[j][i] -= T[q][i] * alpha;
    }
  }
}

struct lp_solve {
  int n, m;
  vector<double> b, c;
  vector<vector<double>> A;

  vector<int> B;
  vector<vector<double>> T;
  void disp() {
    for (int i = 0; i < T.size(); ++i) {
      printf("%2d: ", B[i]);
      for (int j = 0; j < T[i].size(); ++j) {
        printf("  %.1f", T[i][j]);
      }
      printf("\n");
    }
    printf("\n");
  }

  void init() {
    B.assign(m+1, 0);
    n = c.size(); m = b.size();
    T.assign(m+1, vector<double>(n+m+1));
    for (int j = 1; j <= m; ++j) {
      for (int i = 1; i <= n; ++i) T[j][i] = A[j-1][i-1];
      T[j][n+j] = 1;
      T[j][0] = b[j-1];
      B[j] = n+j;
    }
    for (int i = 1; i <= n; ++i) T[0][i] = c[i-1];
  }
  
  double optimize() {
    while (1) { 
      disp();
      int p = 1, q = 1;

      for (int i = 1; i <= n+m; ++i) if (T[0][i] <= T[0][p]) p = i;
      for (int j = 1; j <= m; ++j) if (T[j][0] <= T[q][0]) q = j;
      double t = min(T[0][p], T[q][0]);

      if (t >= -EPS) return -T[0][0]; // optimal
      if (t < T[q][0]) { // tight on c -> primal update
        for (int j = 1; j <= m; ++j) if (T[j][p] >= EPS) 
          if (T[j][p]*fabs(T[q][0]-t) >= T[q][p]*fabs(T[j][0]-t)) q = j;
        if (T[q][p] <= EPS) return INF; // primal infeasible
      } else { // tight on b -> dual update
        for (int i = 0; i <= n+m; ++i) T[q][i] *= -1;
        for (int i = 1; i <= n+m; ++i) 
          if (T[q][i] >= EPS) 
            if (T[q][i]*fabs(T[0][p]-t) >= T[q][p]*fabs(T[0][i]-t)) p = i;
        if (T[q][p] <= EPS) return -INF; // dual infeasible
      }
      B[q] = p;
 
      double alpha = T[q][p];
      for (int i = 0; i <= n+m; ++i) 
        T[q][i] /= alpha;
      for (int j = 0; j <= m; ++j) if (j != q) {
        double alpha = T[j][p];
        for (int i = 0; i <= n+m; ++i) 
          T[j][i] -= T[q][i] * alpha;
      }
    }
  }
};

#include <ctime>

int doit(int seed) {

  //cout << "seed = " << seed << endl;
  srand(seed);

  int n = 3, m = 3;

  double c[n];
  double A[m*n];
  double b[m];
  REP(i, n) c[i] = (rand() % 21) - 10;
  REP(i, m) b[i] = (rand() % 21) - 5;
  REP(i, n) REP(j, m) A[j*n+i] = (rand() % 21) - 5;

  /*
  cout << "c = [";
  REP(i,n) cout << c[i] << (i == n-1 ? "];" : ","); cout << endl;
  cout << "b = [";
  REP(i,m) cout << b[i] << (i == m-1 ? "];" : ";"); cout << endl;

  cout << "A = [";
  REP(j, m) {
    REP(i, n) cout << A[j*n+i] << (i == n-1 ? j == m-1 ? "];" : ";" : ",");
  }
  cout << endl;

  cout << "[sol,opt] = linprog(c,[A;-eye(size(c,2))],[b;zeros(size(c,2),1)])" << endl;
  cout << endl;
  */

  double prev = simplexMethodPD(c, n, b, m, A);

  lp_solve solver;
  solver.c.resize(n);
  for (int i = 0; i < n; ++i) solver.c[i] = c[i];
  solver.b.resize(m);
  for (int i = 0; i < m; ++i) solver.b[i] = b[i];
  solver.A.resize(m, vector<double>(n));
  for (int i = 0; i < m; ++i)
    for (int j = 0; j < n; ++j) solver.A[i][j] = A[n*i+j];
  solver.init();

  double curr = solver.optimize();

  if (fabs(prev - curr) > 1e-4) {
    cout << seed << endl;
  }
}


void moin() {
  lp_solve solver;
  solver.c = {-1, -1, -1};
  solver.b = {1, 1, 1};
  solver.A = {
    {1, 1, 0},
    {1, 0, 1},
    {0, 1, 1},
  };

  solver.init();
  double z = solver.optimize();
  cout << z << endl;
}
int main() {
  moin();
  return 0;

  //doit( time(0) );

  int seed = time(0);
  for (int i = 0; i < 100; ++i) {
    doit(seed + i);
    break;
  }

//  doit( 1364741712 );
}

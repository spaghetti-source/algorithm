//
// Parametric Self-Dual Simplex method
//
// Description:
//   Solve a canonical LP:
//      min. c x
//      s.t. A x <= b
//           x >= 0
// 
//
// Algorithm:
//   Parametric self-dual simplex method.
//
//
// Complexity:
//   O(n+m) iterations on average.
// 
//
// References:
//
// - G. B. Dantzig (1963):
//   Linear Programming and Extensions.
//   Princeton University Press.
//
// - R. J. Vanderbei (2007):
//   Linear programming: Foundations and Extensions.
//   3rd eds., Springer.
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
        if (T[j][p]*(T[q][n+m]-t) >= T[q][p]*(T[j][n+m]-t)) q = j;
      if (T[q][p] <= EPS) return INF; // primal infeasible
    } else { // tight on b -> dual update
      REP(i, n+m+1) T[q][i] *= -1;
      REP(i, n+m) if (T[q][i] >= EPS) 
        if (T[q][i]*(T[m][p]-t) >= T[q][p]*(T[m][i]-t)) p = i;
      if (T[q][p] <= EPS) return -INF; // dual infeasible
    }
    REP(i, m+n+1) if (i != p) T[q][i] /= T[q][p]; T[q][p] = 1; // pivot(q,p)
    REP(j, m+1) if (j != q) {
      double alpha = T[j][p];
      REP(i, n+m+1) T[j][i] -= T[q][i] * alpha;
    }
  }
}

#include <ctime>
int doit(int seed) {
  // verify with matlab
  cout << "seed = " << seed << endl;
  srand(seed);

  int n = 7, m = 3;
  double c[n];
  double A[m*n];
  double b[m];
  REP(i, n) c[i] = (rand() % 21) - 5;
  REP(i, m) b[i] = (rand() % 21) - 5;
  REP(i, n) REP(j, m) A[j*n+i] = (rand() % 21) - 5;

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

  cout << simplexMethodPD(c, n, b, m, A) << endl;
}

int main() { doit( time(0) ); }

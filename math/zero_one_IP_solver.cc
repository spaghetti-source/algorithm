//
// Zero-One IP Solver
//
// Description:
//   Balas's branch-and-bound 0-1 IP sovler.
//   It seeks an integer solution to min c'x s.t. Ax <= b, x in {0,1}^n.
//   The algorithm enumerates x with some branch-and-bound techniques
//   without solving the corresponding LP. 
//
//   It scales about n <= 30.
//
//
// References:
//   E. Balas (1965): 
//   An additive algorithm for solving linear programs with zero-one variables.
//   Operations Research, vol. 13, no. 4, pp. 517--546.
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
#define TEST(s) if (!(s)) { cout << __LINE__ << " " << #s << endl; exit(-1); }

template <class It>
bool next_radix(It begin, It end, int base) {
  for (It cur = begin; cur != end; ++cur) {
    if ((*cur += 1) >= base) *cur = 0;
    else return true;
  }
  return false;
}


const double INF = 1.0 / 0.0;

struct ip_solve {
  int m, n;
  vector<double> c, b;
  vector<vector<double>> A;

  double gauge;
  vector<int> neg;
  ip_solve(vector<double> c, vector<vector<double>> A, vector<double> b) :
    m(b.size()), n(c.size()), c(c), b(b), A(A), gauge(0), neg(n) { 
    // normalize c[j] >= 0
    for (int j = 0; j < n; ++j) {
      if (c[j] >= 0) continue;
      neg[j] = 1;
      gauge += c[j];
      c[j] = -c[j];
      for (int i = 0; i < m; ++i) {
        b[i] -= A[i][j];
        A[i][j] = -A[i][j];
      }
    }
  }

  double obj;
  vector<int> x, sol;
  void explore(double cost) {
    vector<int> history;

    int negatives = 0;
    for (int i = 0; i < m; ++i) {
      if (b[i] >= 0) continue;
      ++negatives;
      double residual = b[i];
      double bound = cost;
      double except = -INF;

      for (int j = 0; j < n; ++j) {
        if (x[j] >= 0) continue;
        if (c[j] + cost >= obj) {
          history.push_back(j);
          x[j] = 0;
        } else if (A[i][j] < 0) {
          residual -= A[i][j];
          bound += c[j];
          except = max(except, A[i][j]);
        }
      }
      if (residual < 0) { negatives = -1; break; }
      if (residual + except < 0) {
        if (bound >= obj) { negatives = -1; break; }
        for (int j = 0; j < n; ++j) {
          if (x[j] < 0 && A[i][j] < 0) {
            history.push_back(j);
            x[j] = 1;
            cost += c[j];
            for (int i = 0; i < m; ++i) 
              b[i] -= A[i][j];
          }
        }
      }
    }
    if (negatives == 0) {
      obj = cost;
      for (int j = 0; j < n; ++j)
        sol[j] = max(0, x[j]);
    } else if (negatives > 0) {
      int p = -1;
      double min_infeasibility = -INF, min_cost = INF;
      for (int j = 0; j < n; ++j) {
        if (x[j] >= 0) continue;
        double infeasibility = 0;
        for (int i = 0; i < m; ++i) 
          infeasibility += min(0.0, b[i] - A[i][j]);
        if (infeasibility > min_infeasibility || 
            (infeasibility == min_infeasibility && c[j] < min_cost)) {
          min_cost = c[j];
          min_infeasibility = infeasibility;
          p = j;
        }
      }
      if (p >= 0) {
        x[p] = 1;
        for (int i = 0; i < m; ++i) b[i] -= A[i][p];
        explore(cost + c[p]);
        x[p] = 0;
        for (int i = 0; i < m; ++i) b[i] += A[i][p];
        explore(cost);
        x[p] = -1;
      }
    }
    for (int j: history) {
      if (x[j] == 1) {
        for (int i = 0; i < m; ++i) 
          b[i] += A[i][j];
      }
      x[j] = -1;
    }
  }

  void solve() {
    obj = INF;
    x = sol= vector<int>(n, -1);
    explore(0);
    obj += gauge;
    for (int i = 0; i < n; ++i) 
      sol[i] ^= neg[i];
  }
};

int main() {
  int n = 50, m = 20;
  vector<double> c(n, 1);

  for (int i = 0; i < n; ++i) 
    c[i] = rand() % 11 - 5;

  vector<double> b(m, -1);

  for (int i = 0; i < m; ++i)
    b[i] = rand() % 11 - 5;

  vector<vector<double>> A(m, vector<double>(n));
  for (int i = 0; i < m; ++i) {
    for (int j = 0; j < n; ++j)
      A[i][j] = rand() % 11 - 5;
  }

  ip_solve solver(c, A, b);
  solver.solve();
  printf("%.0f\n", solver.obj);
}

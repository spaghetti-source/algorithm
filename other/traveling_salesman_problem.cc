// 
// Traveling Salesman Problem (Branch-and-Bound)
//
// Description:
//   It finds TSP by branch-and-bound algorithm that computes
//   a lower bound from the reduced cost matrix.
//   This procedure is similar to the Kuhn-Munkres algorithm.
//
//   In practice, it scales to networks with |V| <= 40.
//
//   
// References:
//   J. D. C. Little, K. G. Murty, D. W. Sweeney, and C. Karel (1963):
//   An algorithm for traveling salesman problem. 
//   Operations Research, vol. 11, pp. 972--989.
//
// Verified:
//   AOJ DPL_2_A: Traveling Salesman Problem

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

const double INF = 1.0 / 0.0;
struct traveling_salesman_problem {
  int n;
  vector<vector<double>> w; // set solver.w[i][j] = infty, etc
  traveling_salesman_problem(int n) : n(n), w(n, vector<double>(n, INF)) { }

  double dynamic_programming() {
    vector<vector<double>> c(n, vector<double>(1<<n, INF));
    vector<vector<int>> p(n, vector<int>(1<<n,-1));
    c[0][1] = 0; // 0 start
    for (size_t S = 0; S < (1<<n); ++S) {
      for (int i = 0; i < n; ++i) {
        if (!(S & (1 << i))) continue;
        for (int j = 0; j < n; ++j) {
          if (S & (1 << j)) continue;
          if (c[j][S|(1<<j)] > c[i][S] + w[i][j]) {
            c[j][S|(1<<j)] = c[i][S] + w[i][j];
            p[j][S|(1<<j)] = i;
          }
        }
      }
    }
    int last;
    double ans = INF;
    for (int i = 0; i < n; ++i) {
      if (ans > c[i][(1<<n)-1] + w[i][0]) {
        ans = c[i][(1<<n)-1] + w[i][0];
        last = i;
      }
    }
    vector<int> path(n);
    int S = (1 << n) - 1;
    path[n-1] = last;
    for (int i = n-1; i >= 1; --i) {
      path[i-1] = p[path[i]][S];
      S ^= (1 << path[i]);
    }
    return ans;
  }
  double brute_force() {
    vector<int> tour;
    for (int i = 0; i < n; ++i)
      tour.push_back(i);
    double ans = INF;
    do {
      double cost = 0;
      for (int i = 0; i < tour.size(); ++i) {
        cost += w[tour[i]][tour[(i+1)%n]];
      }
      ans = min(ans, cost);
    } while (next_permutation(tour.begin()+1, tour.end()));
    return ans;
  }

  vector<int> next, prev, best;
  double upper_bound;
  void explore(int edges, double cost, vector<int> row, vector<int> col) {
    int size = row.size();
    vector<double> rowred(size, INF), colred(size, INF);
    for (int i = 0; i < size; ++i) {
      for (int j = 0; j < size; ++j) 
        rowred[i] = min(rowred[i], w[row[i]][col[j]]);
      for (int j = 0; j < size; ++j) 
        if (w[row[i]][col[j]] < INF) 
          w[row[i]][col[j]] -= rowred[i];
      cost += rowred[i];
    }
    for (int j = 0; j < size; ++j) {
      for (int i = 0; i < size; ++i) 
        colred[j] = min(colred[j], w[row[i]][col[j]]);
      for (int i = 0; i < size; ++i) 
        if (w[row[i]][col[j]] < INF) 
          w[row[i]][col[j]] -= colred[j];
      cost += colred[j];
    }
    if (cost < upper_bound) {
      if (edges == n - 2) {
        upper_bound = cost;
        best = next;
        best[row[0]] = col[0];
        best[row[1]] = col[1];
        if (w[row[0]][col[0]] >= INF) 
          swap(best[row[0]], best[row[1]]);
      } else {
        int r, c;
        double most = -INF; 
        for (int i = 0; i < size; ++i) {
          for (int j = 0; j < size; ++j) {
            if (w[row[i]][col[j]] != 0) continue;
            double minrow = INF, mincol = INF;
            for (int k = 0; k < size; ++k) {
              if (i != k) minrow = min(minrow, w[row[k]][col[j]]);
              if (j != k) mincol = min(mincol, w[row[i]][col[k]]);
            }
            if (most < minrow + mincol) {
              most = minrow + mincol;
              r = i;
              c = j;
            }
          }
        }
        next[row[r]] = col[c];
        prev[col[c]] = row[r];

        int last = col[c], first = row[r];
        while (next[last] >= 0) last = next[last];
        while (prev[first] >= 0) first = prev[first];
        double colrowval = w[last][first];
        w[last][first] = INF;
        vector<int> newrow = row, newcol = col;
        newrow.erase(newrow.begin() + r);
        newcol.erase(newcol.begin() + c);
        explore(edges + 1, cost, newrow, newcol);
        w[last][first] = colrowval;

        prev[col[c]] = -1;
        next[row[r]] = -1;

        if (cost + most < upper_bound) {
          w[row[r]][col[c]] = INF;
          explore(edges, cost, row, col);
          w[row[r]][col[c]] = 0;
        }
      }
    }
    for (int i = 0; i < size; ++i)
      for (int j = 0; j < size; ++j)
        w[row[i]][col[j]] += rowred[i] + colred[j];
  }

  double solve() {
    vector<int> row(n), col(n);
    for (int i = 0; i < n; ++i) row[i] = col[i] = i;
    upper_bound = INF;
    prev = next = vector<int>(n, -1);
    explore(0, 0, row, col);

    if (upper_bound < INF) {
      vector<int> path(n); // solution
      path[0] = 0;
      for (int i = 1; i < n; ++i) 
        path[i] = best[path[i-1]];
    }
    return upper_bound;
  }
};


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
int main() {
  for (int iter = 0; iter < 1000; ++iter) {
    srand(iter);
    int n = 40;
    traveling_salesman_problem tsp(n);
    for (int i = 0; i < n; ++i)
      for (int j = 0; j < n; ++j)
        if (i != j) tsp.w[i][j] = rand() % 1000 + 1;

    tick();
    double a = 0;// tsp.brute_force();
    double ta = tick();
    double b = 0;// tsp.dynamic_programming();
    double tb = tick();
    double c = tsp.solve();
    double tc = tick();

    a = b =c;
    printf("%f %f\n%f %f\n%f %f\n\n", a, ta, b, tb, c, tc);
    if (a != b || b != c) {
      printf("seed = %d\n", iter);
      exit(-1);
    }
  }
}
*/

// AOJ: DPL_2_A
int main() {
  int n, m;
  scanf("%d %d", &n, &m);
  traveling_salesman_problem tsp(n);
  for (int i = 0; i < m; ++i) {
    int u, v, w;
    scanf("%d %d %d", &u, &v, &w);
    tsp.w[u][v] = w;
  }
  double ans = tsp.solve();
  if (ans < INF) printf("%.0f\n", ans);
  else           printf("-1\n");
}

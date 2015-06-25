//
// Exact Cover 
//
// Description:
//   We are given a family of sets F on [0,n).
//   The exact cover problem is to find a subfamily of F
//   such that each k in [0,n) is covered exactly once.
//   For example, if F consists from
//     {1,2}, {1,2,3}, {3,4}, {5,6}, {3,5,6},
//   then the exact cover is 
//     {1,2}, {3,4}, {5,6}.
//     
// Algorithm:
//   Knuth's algorithm X is the following recursive algorithm:
//
//     select some k in [0,n)
//     for each subset S that covers k
//       select S and remove all conflicting sets
//       recursion
//    
//    To implement this algorithm efficiently, we can use
//    a data structure, which is called dancing links. 
//   
// Verified: 
//   SPOJ 1428: EASUDOKU
//   SPOJ 1110: SUDOKU
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

vector<int> exact_cover(vector<vector<int>> sets) { 
  int m = 0, M = 10;
  for (auto &v: sets) {
    m = max(m, *max_element(all(v))+1);
  }
  M += (1 + sets.size()) * m;
  vector<int> L(M), R(M), U(M), D(M), S(M), C(M), A(M);
  for (int i = 0; i <= m; ++i) {
    L[i] = i-1; R[i] = i+1;
    D[i] = U[i] = C[i] = i;
  }
  L[0] = m; R[m] = 0;
  int p = m+1;
  for (int row = 0; row < sets.size(); ++row) { // add sets
    for (int i = 0; i < sets[row].size(); ++i) {
      int col = sets[row][i];
      C[p] = col; A[p] = row; ++S[col];
      D[p] = D[col]; U[p] = col; D[col] = U[D[p]] = p;
      if (i == 0) { L[p] = R[p] = p; }
      else { L[p] = p-1; R[p] = R[p-1]; R[p-1] = L[R[p]] = p; }
      ++p;
    }
  }
  auto remove = [&](int x) {
    L[R[x]] = L[x]; R[L[x]] = R[x];
    for (int i = D[x]; i != x; i = D[i]) 
      for (int j = R[i]; j != i; j = R[j]) 
      { U[D[j]] = U[j], D[U[j]] = D[j], --S[C[j]]; }
  };
  auto resume = [&](int x) {
    for (int i = U[x]; i != x; i = U[i]) 
      for (int j = L[i]; j != i; j = L[j]) 
      { U[D[j]] = j, D[U[j]] = j, ++S[C[j]]; }
    L[R[x]] = x; R[L[x]] = x;
  };

  vector<int> solution;
  function<bool(void)> rec = [&]() {
    if (R[m] == m) return true; // found
    int col = R[m];
    for (int i = R[m]; i != m; i = R[i]) 
      if (S[i] < S[col]) col = i;
    if (S[col] == 0) return false;
    remove(col);
    for (int i = D[col]; i != col; i = D[i]) {
      solution.push_back(A[i]);
      for (int j = R[i]; j != i; j = R[j]) remove(C[j]);
      if (rec()) return true;
      for (int j = L[i]; j != i; j = L[j]) resume(C[j]);
      solution.pop_back();
    }
    resume(col);
    return false;
  };
  rec();
  return solution;
}



int w = 3;
bool sudoku(vector<vector<int>> &b) {
  vector<vector<int>> sets;
  vector<tuple<int,int,int>> ns;

  auto id = [](int a, int b, int c) { return w*w*w*w*a + w*w*b + c; };
  auto add_set =[&](int i, int j, int k) {
    sets.push_back({});
    sets.back().push_back(id(0, i, j));
    sets.back().push_back(id(1, i, k));
    sets.back().push_back(id(2, j, k));
    sets.back().push_back(id(3, w*(i/w)+(j/w), k));
    ns.push_back(make_tuple(i, j, k));
  };
  for (int i = 0; i < w*w; ++i) {
    for (int j = 0; j < w*w; ++j) {
      if (b[i][j] == 0) {
        for (int k = 0; k < w*w; ++k) 
          add_set(i, j, k);
      } else {
        add_set(i, j, b[i][j]-1);
      }
    }
  }
  auto x = exact_cover(sets);
  if (x.empty()) return false;
  for (auto a: x) {
    int i = get<0>(ns[a]);
    int j = get<1>(ns[a]);
    int k = get<2>(ns[a]);
    b[i][j] = k;
  }
  return true;
}

void SPOJ_EASUDOKU() {
  w = 3;
  int ncase; scanf("%d", &ncase);
  for (int icase = 0; icase < ncase; ++icase) {
    vector<vector<int>> b(w*w, vector<int>(w*w));
    for (int i = 0; i < w*w; ++i) 
      for (int j = 0; j < w*w; ++j) 
        scanf("%d", &b[i][j]);
    if (sudoku(b)) {
      for (int i = 0; i < w*w; ++i) {
        for (int j = 0; j < w*w; ++j) {
          if (j > 0) printf(" ");
          printf("%d", b[i][j]+1);
        }
        printf("\n");
      }
      printf("\n");
    } else {
      printf("No solution\n");
    }
  }
}
void SPOJ_SUDOKU() {
  w = 4;
  int ncase; scanf("%d", &ncase);
  for (int icase = 0; icase < ncase; ++icase) {
    if (icase > 0) printf("\n");
    vector<vector<int>> b(w*w, vector<int>(w*w));
    for (int i = 0; i < w*w; ++i) {
      char s[1024];
      scanf("%s", s);
      for (int j = 0; j < w*w; ++j) {
        char c = s[j];
        if (c == '-') b[i][j] = 0;
        else          b[i][j] = c - 'A' + 1;
      }
    }
    if (sudoku(b)) {
      for (int i = 0; i < w*w; ++i) {
        for (int j = 0; j < w*w; ++j) {
          printf("%c", b[i][j]+'A');
        }
        printf("\n");
      }
    }
  }

}

int main() {
  // SPOJ_EASUDOKU();
   SPOJ_SUDOKU();
}

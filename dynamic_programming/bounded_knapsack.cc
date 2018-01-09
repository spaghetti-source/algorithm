//
// Bounded Knapsack Problem
//
// Description:
// 
//   There are n kinds of items of profit pi, weight bi, and 
//   amount mi (i in {0...n-1}). For a given W, we want to select items 
//   to maximize the total profit subject to the total weight
//   is at most W and the number of each item is at most mi.
//
//   This problem can be solved by the following DP.
//     E[j][w] = max {E[j-1][w], E[j-1][w-bj]+pj, E[j-1][w-2bj]+2pj. ...}
//   A naive implementation requires O(nmW) time. However, we can reduce
//   this to O(nW) as follows. 
//
//   We compute E[j][s], E[j][s+bj]. E[s+2bj] ... for each s in [0,bj).
//   For simplicity, we consider s = 0. Then, we have
//     E[j][w+bj]  = max {E[j-1][w+bj], E[j-1][w]+pj, E[j-1][w-bj]+2pj. ...}
//   By comparing this with the original formula, 
//   - E[j][w+bj] contains E[j-1][w+bj] term
//   - E[j][w+bj] does not contain E[j-1][w-mjbj] term
//   - The all terms have been added pj
//   Thus, by using a data structure that supports these operations, 
//   we can perform the DP efficiently. The data structure is implemented 
//   by a maximum queue with one accumulation parameter.
//
// Complexity:
//   
//   O(n W)
//
#include <bits/stdc++.h>

using namespace std;

#define fst first
#define snd second
#define all(c) ((c).begin()), ((c).end())
#define TEST(s) if (!(s)) { cout << __LINE__ << " " << #s << endl; exit(-1); }

int boundedKnapsackDP(vector<int> ps,
                       vector<int> ws,
                       vector<int> ms,
                       int W) {
  int n = ps.size();
  vector<vector<int>> dp(n+1, vector<int>(W+1));
  for (int i = 0; i < n; ++i) {
    for (int s = 0; s < ws[i]; ++s) {
      int alpha = 0;
      queue<int> que;
      deque<int> peek;
      for (int w = s; w <= W; w += ws[i]) {
        alpha += ps[i];
        int a = dp[i][w]-alpha;
        que.push(a);
        while (!peek.empty() && peek.back() < a) peek.pop_back();
        peek.push_back(a);
        while (que.size() > ms[i]+1) {
          if (que.front() == peek.front()) peek.pop_front();
          que.pop();
        }
        dp[i+1][w] = peek.front()+alpha;
      }
    }
  }
  int ans = 0;
  for (int w = 0; w <= W; ++w) 
    ans = max(ans, dp[n][w]);
  return ans;
}

int boundedKnapsackDPNaive(vector<int> ps,
                           vector<int> ws,
                           vector<int> ms,
                           int W) {
  int n = ps.size();
  vector<vector<int>> dp(n+1, vector<int>(W+1));
  for (int i = 0; i < n; ++i) {
    for (int w = 0; w <= W; ++w) {
      dp[i+1][w] = dp[i][w];
      for (int j = 1; j <= ms[i]; ++j) {
        if (w - j*ws[i] < 0) break;
        dp[i+1][w] = max(dp[i+1][w], dp[i][w-j*ws[i]]+j*ps[i]);
      }
    }
  }
  int ans = 0;
  for (int w = 0; w <= W; ++w) 
    ans = max(ans, dp[n][w]);
  return ans;
}

int main() {
  int seed; 
  seed = 20;
  //cin >> seed;
  seed = time(0);
  srand(seed);
  int n = 100;
  int W = 10000;
  vector<int> ps(n), ws(n), ms(n);
  for (int i = 0; i < n; ++i) {
    ps[i] = rand() % n + 1;
    ws[i] = rand() % n + 1;
    ms[i] = rand() % n + 1;
    //cout << ps[i] << " " << ws[i] << " " << ms[i] << endl;
  }

  cout << boundedKnapsackDP(ps, ws, ms, W) << endl;
  cout << boundedKnapsackDPNaive(ps, ws, ms, W) << endl;
}

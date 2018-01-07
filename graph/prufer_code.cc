// 
// Prufer Code in Linear Time
//
// Description:
//   Prufer code gives one-to-one correspondence
//   between labeled trees and integers of length n-2.
//   This immediately shows the Cayley theorem:
//   the number of labeled trees is n^{n-2}.
//   The algorithm computes the Prufer code in linear
//   time. 
//
// Complexity:
//   O(n)
//
// Reference:
//   Xiaodong Wang, Lei Wang, and Yingjie Wui (2009):
//   "An Optimal Algorithm for Prufer Codes", 
//   Journal of Software Engineering and Applications, 
//   vol.2, 111--115.
//

#include <bits/stdc++.h>

using namespace std;

#define fst first
#define snd second
#define all(c) ((c).begin()), ((c).end())
#define TEST(s) if (!(s)) { cout << __LINE__ << " " << #s << endl; exit(-1); }

struct Tree {
  int n;
  vector<vector<int>> adj;
  Tree(int n) : n(n), adj(n) { }
  void addEdge(int u, int v) {
    adj[u].push_back(v);
    adj[v].push_back(u);
  }
};

vector<int> labeledTreeToCode(Tree T) {
  vector<int> deg(T.n), parent(T.n, -1), code;
  function<void(int)> dfs = [&](int u) {
    deg[u] = T.adj[u].size();
    for (int v: T.adj[u]) {
      if (v != parent[u]) {
        parent[v] = u;
        dfs(v);
      }
    }
  }; dfs(T.n-1);

  int index = -1;
  while (deg[++index] != 1);
	for (int u = index, i = 0; i < T.n-2; ++i) {
    int v = parent[u];
		code.push_back(v);
    if (--deg[v] == 1 && v < index) {
      u = v;
    } else {
      while (deg[++index] != 1);
      u = index;
    }
	}
	return code;
}

Tree codeToLabeledTree(vector<int> code) {
	int n = code.size() + 2;
  Tree T(n);
	vector<int> deg(n, 1);
  for (int i = 0; i < n-2; ++i)
    ++deg[code[i]];

  int index = -1;
  while (deg[++index] != 1);
  for (int u = index, i = 0; i < n-2; ++i) {
    int v = code[i];
		T.addEdge(u, v);
		--deg[u]; --deg[v];
		if (deg[v] == 1 && v < index) {
			u = v;
		} else {
			while (deg[++index] != 1);
			u = index;
		}
	}
	for (int u = 0; u < n-1; ++u) 
		if (deg[u] == 1) T.addEdge(u, n-1);
	return T;
}

int main() {
  Tree T(6);
  T.addEdge(0, 3);
  T.addEdge(1, 3);
  T.addEdge(2, 3);
  T.addEdge(3, 4);
  T.addEdge(4, 5);
  auto code = labeledTreeToCode(T);
  for (int u: code) {
    cout << u << " ";
  }
  cout << endl;

  Tree G = codeToLabeledTree(code);
  for (int u = 0; u < G.adj.size(); ++u) {
    for (int v: G.adj[u]) {
      if (u < v) cout << u << " " << v << endl;
    }
  }
}

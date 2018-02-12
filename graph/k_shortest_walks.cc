//
// K Shortest Walks (Simplified Eppstein)
//
// Description
//
//   We are given a weighted graph. The k-shortest walks problem
//   seeks k different s-t walks (paths allowing repeated vertices)
//   in the increasing order of the lengths.
//
//   If we maintain each walks explicitly, it must costs O(k^2 m) time.
//   To avoid this complexity, we maintain the walks in a compact format.
//   Let us fix a reverse shortest path tree from t. A deviation is an 
//   edge that is not on the tree. Any walk is represented by a concatenation
//   of deviations and paths on the tree. We enumerate all possible 
//   deviations and use the best-first search to find the k-th solution.
//
//   The Eppstein's algorithm maintains the set of deviations by 
//   the augmented persistent heaps and emurates the best-first search. 
//   Here, we implemented a simplified version of the Eppstein's algorithm,
//   which uses the simple persistent heaps instead of the augmented ones.
//   It increases the space from O(m + n log n) to O(m log n).
//   
// Complexity:
//
//   O(m log m) construction
//   O(k log k) for k-th search
//
// Verified:
//
//   UTPC2013_10 J K-th Cycle
//
// References:
//
//   David Eppstein (1998): 
//   "Finding the k shortest paths",
//   SIAM Journal on computing, vol.28, no.2, pp.652--673.
//
#include <bits/stdc++.h>

using namespace std;

#define fst first
#define snd second
#define all(c) ((c).begin()), ((c).end())
#define TEST(s) if (!(s)) { cout << __LINE__ << " " << #s << endl; exit(-1); }

struct Graph {
  int n, m = 0;
  vector<int> head;                 // Vertex
  vector<int> src, dst, next, prev; // Edge

  using Weight = long long;
  vector<Weight> weight;
  Graph(int n) : n(n), head(n, -1) { }
  int addEdge(int u, int v, Weight w) {
    next.push_back(head[u]);
    src.push_back(u);
    dst.push_back(v);
    weight.push_back(w);
    return head[u] = m++;
  }
};
constexpr Graph::Weight INF = 1e15;

struct KShortestWalks {
  Graph g;
  vector<Graph::Weight> dist;
  vector<int> tree, order;
  void reverseDijkstra(int t) {
    vector<vector<int>> adj(g.n);
    for (int u = 0; u < g.n; ++u) 
      for (int e = g.head[u]; e >= 0; e = g.next[e]) 
        adj[g.dst[e]].push_back(e);
    dist.assign(g.n, INF);
    tree.assign(g.n, ~g.m);
    using Node = tuple<Graph::Weight,int>;
    priority_queue<Node, vector<Node>, greater<Node>> que;
    que.push(make_tuple(0, t));
    dist[t] = 0;
    while (!que.empty()) {
      int u = get<1>(que.top()); que.pop();
      if (tree[u] >= 0) continue;
      tree[u] = ~tree[u];
      order.push_back(u);
      for (int e: adj[u]) {
        int v = g.src[e];
        if (dist[v] > dist[u] + g.weight[e]) {
          tree[v] = ~e;
          dist[v] = dist[u] + g.weight[e];
          que.push(Node(dist[v], v));
        }
      }
    }
  }
  struct Node { // Persistent Heap (Leftist Heap)
    int e;
    Graph::Weight delta;
    Node *left = 0, *right = 0;
    int rnk = 0;
  } *root = 0;
  static Node *merge(Node *x, Node *y) {
    if (!x) return y;
    if (!y) return x;
    if (x->delta > y->delta) swap(x, y);
    x = new Node(*x);
    x->right = merge(x->right, y);
    if (!x->left || x->left->rnk < x->rnk) swap(x->left, x->right);
    x->rnk = (x->right ? x->right->rnk : 0) + 1;
    return x;
  }
  vector<Node*> deviation;
  void buildHeap() {
    deviation.resize(g.n);
    for (int u: order) {
      int v = -1;
      for (int e = g.head[u]; e >= 0; e = g.next[e]) {
        if (tree[u] == e) v = g.dst[e]; 
        else if (dist[g.dst[e]] < INF) {
          auto delta = g.weight[e] - dist[g.src[e]] + dist[g.dst[e]];
          deviation[u] = merge(deviation[u], new Node({e, delta}));
        }
      }
      if (v >= 0) deviation[u] = merge(deviation[u], deviation[v]);
    }
  }
  KShortestPaths(Graph g_, int t) : g(g_) {
    reverseDijkstra(t);
    buildHeap();
  }
  void enumerate(int s, int kth) {
    int k = 0;
    Node *x = deviation[s];
    Graph::Weight len = dist[s];
    ++k;
    using SearchNode = tuple<Node*, Graph::Weight>;
    auto comp = [](SearchNode x, SearchNode y) { return get<1>(x) > get<1>(y); };
    priority_queue<SearchNode, vector<SearchNode>, decltype(comp)> que(comp);
    if (x) que.push(SearchNode(x, len + x->delta));
    while (!que.empty() && k < kth) {
      tie(x, len) = que.top(); que.pop();
      int e = x->e, u = g.src[e], v = g.dst[e];
      cout << len << endl; ++k;
      if (deviation[v]) que.push(SearchNode(deviation[v], len+deviation[v]->delta));
      for (Node *y: {x->left, x->right}) 
        if (y) que.push(SearchNode(y, len + y->delta-x->delta));
    }
    while (k < kth) { cout << -1 << endl; ++k; }
  }
};

void KSH_test() {
  int n = 4;
  Graph g(n);
  g.addEdge(0, 1, 2);
  g.addEdge(0, 2, 2);
  g.addEdge(1, 3, 4);
  g.addEdge(2, 3, 2);
  g.addEdge(1, 2, 1);
  g.addEdge(2, 1, 1);
  KShortestPaths ksh(g, 3);
  ksh.enumerate(0, 10);
}

void UTPC2013_10() {
  int n, m, k;
  scanf("%d %d %d", &n, &m, &k);
  Graph g(n);
  for (int i = 0; i < m; ++i) {
    int u, v;
    long long w;
    scanf("%d %d %lld", &u, &v, &w);
    g.addEdge(u, v, w);
  }
  KShortestPaths ksh(g, 0);
  ksh.enumerate(0, k+1);
}

int main() {
  UTPC2013_10();
  //KSH_test();
}

//
// Graph Isomorphism
//
// Description:
//   Directed graphs G and H are isomprhic
//   iff there is a bijection between these vertices such that
//     (u,v) in E(G) iff (pi(u),pi(v)) in E(H).
//   It is open that testing graph isomorphism is NP-complete.
// 
// Algorithm:
//   Cordella-Foggia-Sansone-Vento's algorithm (aka, VF2 algorithm)
//   with large degree heuristics.
//
// Complexity:
//   O(d n!) time, O(n) space.
//
// Verify:
//   SPOJ401
//
// References:
//   L. P. Cordella, P. Foggia, C. Sansone, and M. Vento (2004):
//   A (sub)graph isomorphism algorithm for matching large graphs.
//   IEEE Transactions on Pattern Analysis and Machine Intelligence,
//   vol.28, no.10, pp.1367--1372.

#include <iostream>
#include <algorithm>
#include <vector>
#include <queue>
#include <set>
#include <map>
#include <cstdio>
#include <cstring>
#include <unordered_set>

using namespace std;
#define fst first
#define snd second
#define all(c) ((c).begin()), ((c).end())

struct graph {
  int n;
  vector<unordered_set<int>> adj, rdj;
  graph(int n) : n(n), adj(n), rdj(n) { } 
  void add_edge(int src, int dst) {
    adj[src].insert(dst);
    rdj[dst].insert(src);
  }
};


bool isomorphic(graph G, graph H) {
  if (G.n != H.n) return false;
  int n = G.n;
  vector<int> GtoH(n, -1), HtoG(n, -1); 
  vector<int> inG(n), outG(n), inH(n), outH(n); 

  function<bool (int)> match = [&](int k) {
    if (k == n) return true;

    int s = -1; // large degree heuristics
    for (int t = 0; t < n; ++t) {
      if (GtoH[t] >= 0) continue;
      if (s == -1 || inG[s] + outG[s] < inG[t] + outG[t]) s = t;
    }
    for (int u = 0; u < n; ++u) {
      auto check = [&](int s, int u) { // VF2 feasibility function
        if (HtoG[u] >= 0) return false;
        if (inG[s] != inH[u]) return false;
        if (outG[s] != outH[u]) return false;
        int TinG = 0, ToutG = 0, TnewG = 0, TinH = 0, ToutH = 0, TnewH = 0;
        for (int t: G.adj[s]) {
          if (GtoH[t] >= 0) { 
            if (!H.adj[u].count(GtoH[t])) return false;
          } else {
            if (inG[t]) ++TinG;
            if (outG[t]) ++ToutG;
            if (!inG[t] && !outG[t]) ++TnewG;
          }
        }
        for (int v: H.adj[u]) {
          if (HtoG[v] >= 0) {
            if (!G.adj[s].count(HtoG[v])) return false;
          } else {
            if (inH[v]) ++TinH;
            if (outH[v]) ++ToutH;
            if (!inH[v] && !outH[v]) ++TnewH;
          }
        }
        if (TinG != TinH || ToutG != ToutH || TnewG != TnewH) return false;
        for (int t: G.rdj[s]) {
          if (GtoH[t] >= 0) {
            if (!H.rdj[u].count(GtoH[t])) return false;
          } else {
            if (inG[t]) ++TinG;
            if (outG[t]) ++ToutG;
            if (!inG[t] && !outG[t]) ++TnewG;
          }
        }
        for (int v: H.rdj[u]) {
          if (HtoG[v] >= 0) {
            if (!G.rdj[s].count(HtoG[v])) return false;
          } else {
            if (inH[v]) ++TinH;
            if (outH[v]) ++ToutH;
            if (!inH[v] && !outH[v]) ++TnewH;
          }
        }
        if (TinG != TinH || ToutG != ToutH || TnewG != TnewH) return false;
        return true;
      };
      if (!check(s, u)) continue;
      ++inG[s]; ++outG[s];
      for (int t: G.rdj[s]) ++inG[t];
      for (int t: G.adj[s]) ++outG[t];
      ++inH[u]; ++outH[u];
      for (int v: H.rdj[u]) ++inH[v];
      for (int v: H.adj[u]) ++outH[v];

      GtoH[s] = u; HtoG[u] = s;
      if (match(k+1)) return true;
      GtoH[s] = -1; HtoG[u] = -1;

      --inG[s]; --outG[s];
      for (int t: G.rdj[s]) --inG[t];
      for (int t: G.adj[s]) --outG[t];
      --inH[u]; --outH[u];
      for (int v: H.rdj[u]) --inH[v];
      for (int v: H.adj[u]) --outH[v];
    }
    return false;
  };
  return match(0);
}


// --- tester ---

// slow; but correct
bool isomorphic_naive(graph G, graph H) {
  if (G.n != H.n) return false;
  int n = G.n;
  vector<int> GtoH(n, -1), HtoG(n, -1); // matched pair

  function<bool (int)> match = [&](int s) {
    if (s == n) {
      for (int t = 0; t < n; ++t) {
        vector<bool> nbhG(n), nbhH(n);
        for (int p: G.adj[t]) nbhG[GtoH[p]] = 1;
        for (int q: H.adj[GtoH[t]]) nbhH[q] = 1;
        if (nbhG != nbhH) return false;
      }
      return true;
    } else {
      for (int u = 0; u < n; ++u) {
        if (HtoG[u] >= 0) continue;
        GtoH[s] = u; HtoG[u] = s;
        if (match(s+1)) return true;
        GtoH[s] = -1; HtoG[u] = -1;
      }
      return false;
    }
  };
  return match(0);
}




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
  // correctness
  if (0) {
    for (int seed = 0; seed < 10000; ++seed) {
      int n = 5, r = 2;
      graph g(n), h(n);
      vector<int> deg(n);

      for (int u = 0; u < n; ++u) {
        set<int> adj;
        while (adj.size() < r) {
          int v = rand() % n;
          if (v != u) adj.insert(v);
        }
        for (auto v: adj)
          h.add_edge(u, v);
      }
      for (int u = 0; u < n; ++u) {
        set<int> adj;
        while (adj.size() < r) {
          int v = rand() % n;
          if (v != u) adj.insert(v);
        }
        for (auto v: adj)
          g.add_edge(u, v);
      }
      int naive = isomorphic_naive(g, h);
      int algo  = isomorphic(g, h);
      if (naive != algo) {
        cout << "incorrect" << endl;
        exit(-1);
      }
    }
    cout << "correctness passed" << endl;
  }
  // efficiency

  if (1) {
    int n = 1000, r = 10;
    graph g(n), h(n);
    vector<int> deg(n);

    for (int u = 0; u < n; ++u) {
      set<int> adj;
      while (adj.size() < r) {
        int v = rand() % n;
        if (v != u) adj.insert(v);
      }
      for (auto v: adj)
        h.add_edge(u, v);
    }
    for (int u = 0; u < n; ++u) {
      set<int> adj;
      while (adj.size() < r) {
        int v = rand() % n;
        if (v != u) adj.insert(v);
      }
      for (auto v: adj)
        g.add_edge(u, v);
    }
    tick();
    int algo  = isomorphic(g, h);
    cout << algo << endl;
    cout << "efficiency: " << tick() << endl;
  }
  if (1) {
    int n = 1000, r = 10;
    vector<int> pi(n); iota(all(pi), 0);
    random_shuffle(all(pi));

    graph g(n), h(n);
    vector<int> deg(n);

    for (int u = 0; u < n; ++u) {
      set<int> adj;
      while (adj.size() < r) {
        int v = rand() % n;
        if (v != u) adj.insert(v);
      }
      for (auto v: adj) {
        h.add_edge(u, v);
        g.add_edge(pi[u], pi[v]);
      }
    }
    tick();
    int algo  = isomorphic(g, h);
    cout << algo << endl;
    cout << "efficiency: " << tick() << endl;
  }
}

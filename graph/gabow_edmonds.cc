//
// General Graph Matching (Gabow-Edmonds)
//
// Description:
//   
//   For a graph G = (V, E), a matching M is a set of edges
//   such that any vertex is contained in M at most once.
//   The matching with maximum cardinality is computed by
//   the Edmonds blossom algorithm.
//
//   This implementation is the Gabow's simplified version
//   with the lazy update technique to improve the complexity
//   in sparse graphs.
//
//
// Complexity:
//
//   O(n m log n)
//
//
// Verified:
//
//   SPOJ ADABLOOM
//
//
// References:
//   H.Gabow (1976):
//   An efficient implementation of Edmonds' algorithm for maximum matching on graphs.
//   Journal of the ACM, vol.23, no.2, pp.221-234.
//
//   http://min-25.hatenablog.com/entry/2016/11/21/222625
//
//

// g++ -std=c++17 -O3 -fmax-errors=1 -fsanitize=undefined
#include <bits/stdc++.h>
using namespace std;

#define fst first
#define snd second
#define all(c) ((c).begin()), ((c).end())

struct Graph {
  int n;
  vector< vector<int> > adj;
  Graph(int n) : n(n), adj(n) { };
  void addEdge(int u, int v) {
    adj[u].push_back(v);
    adj[v].push_back(u);
  }

  vector<int> mate;
  int maximumMatching() {
    mate.assign(n+1, n);
    vector<int> first(n+1, n), que(n);
    vector<pair<int,int>> label(n+1, make_pair(-1,-1));
    int head = 0, tail = 0;

    function<void(int,int)> rematch = [&](int v, int w) {
      int t = mate[v]; mate[v] = w;
      if (mate[t] != v) return;
      if (label[v].snd == -1) {
        mate[t] = label[v].fst;
        rematch(mate[t], t);
      } else {
        int x, y; tie(x, y) = label[v];
        rematch(x, y); rematch(y, x);
      }
    };
    auto relabel = [&](int x, int y) {
      function<int(int)> findFirst = [&](int u) {
        return label[first[u]].fst < 0 ? first[u] :
               first[u] = findFirst(first[u]);
      };
      int r = findFirst(x), s = findFirst(y);
      if (r == s) return; 
      auto h = make_pair(~x, y);
      label[r] = label[s] = h;
      int join;
      while (1) {
        if (s != n) swap(r, s);
        r = findFirst(label[mate[r]].fst);
        if (label[r] == h) {
          join = r;
          break;
        } else {
          label[r] = h;
        }
      }
      for (int v: {first[x], first[y]}) {
        for (; v != join; v = first[label[mate[v]].fst]) {
          label[v] = make_pair(x, y);
          first[v] = join;
          que[tail++] = v;
        }
      }
    };
    auto augment = [&](int u) {
      label[u] = make_pair(n, -1);
      first[u] = n;
      head = tail = 0;
      for (que[tail++] = u; head < tail;) {
        int x = que[head++];
        for (int y: adj[x]) {
          if (mate[y] == n && y != u) {
            mate[y] = x;
            rematch(x, y);
            return true; 
          } else if (label[y].fst >= 0) {
            relabel(x, y);
          } else if (label[mate[y]].fst == -1) {
            label[mate[y]].fst = x;
            first[mate[y]] = y;
            que[tail++] = mate[y];
          } 
        }
      }
      return false;
    };
    int matching = 0;
    for (int u = 0; u < n; ++u) {
      if (mate[u] < n || !augment(u)) continue;
      ++matching;
      for (int i = 0; i < tail; ++i) 
        label[que[i]] = label[mate[que[i]]] = make_pair(-1,-1);
      label[n] = make_pair(-1, -1);
    }
    return matching;
  }
};


void LA3820() {
  int ncase; scanf("%d", &ncase);
  for (int icase = 0; icase < ncase; ++icase) {
    int n, m; 
    scanf("%d %d", &n, &m);
    vector<int> v(n);
    for (int i = 0; i < n; ++i) {
      scanf("%d", &v[i]);
    }
    set<int> S;
    for (int i = 0; i < m; ++i) {
      int x;
      scanf("%d", &x);
      S.insert(x);
    }
    Graph G(n);
    for (int i = 0; i < n; ++i) {
      for (int j = i+1; j < n; ++j) {
        if (S.count(v[i] + v[j])) 
          G.addEdge(i, j);
      }
    }
    printf("%d\n", G.maximumMatching());
  }
}

void UOJ79() {
  int n, m; cin >> n >> m;
  Graph g(n);
  for (int i = 0; i < m; ++i) {
    int u, v;
    cin >> u >> v;
    g.addEdge(u-1, v-1);
  }
  cout << g.maximumMatching() << endl;
  for (int u = 0; u < n; ++u) {
    if (u > 0) cout << " ";
    if (g.mate[u] >= n) cout << 0;
    else                cout << g.mate[u]+1;
  }
  cout << endl;
}
void SPOJ_ADABLOOM() {
  int ncase; scanf("%d", &ncase);
  for (int icase; icase < ncase; ++icase) {
    int n; scanf("%d", &n);
    vector<long long> a(n);
    for (int i = 0; i < n; ++i) 
      scanf("%lld", &a[i]);
    random_shuffle(all(a));
    Graph g(n);
    for (int i = 0; i < n; ++i) { 
      for (int j = 0; j < n; ++j) {
        if (a[i] < (a[i] ^ a[j]) && (a[i] ^ a[j]) < a[j]) g.addEdge(i, j);
      }
    }
    cout << g.maximumMatching() << endl;
  }
}

int main() {
  SPOJ_ADABLOOM();
  //UOJ79();
  //LA3820();
}

//
// Counting Perfect Matchings in Plane Graph 
// (Fisher, Kasteleyn, and Temperley)
//
// Description:
//
//   Pfaffian Orientation; see https://en.wikipedia.org/wiki/FKT_algorithm
//
// Complexity:
//
//   O(n^3).
//
// g++ -std=c++17 -O3 -fmax-errors=1 -fsanitize=undefined
#include <bits/stdc++.h>

using namespace std;

#define fst first
#define snd second
#define all(c) ((c).begin()), ((c).end())
#define TEST(s) if (!(s)) { cout << __LINE__ << " " << #s << endl; exit(-1); }

using Real = long double;
struct Point { 
  Real x, y;
};

struct PlaneGraph {
  vector<int> incident_edge; // vertex record
  vector<int> origin, twin, prev, next, incident_face; // edge record
  vector<int> component; // face record
  int edges()    const { return origin.size(); }
  int vertices() const { return incident_edge.size(); }
  int faces()    const { return component.size(); }

  vector<Point> point;
  int newVertex(Point p, int e = -1) {
    point.push_back(p);
    incident_edge.push_back(e);
    return vertices()-1;
  }
  int newEdge(int o = -1) {
    origin.push_back(o);
    twin.push_back(-1);
    prev.push_back(-1);
    next.push_back(-1);
    incident_face.push_back(-1);
    return edges()-1;
  }
  int newFace(int e = -1) {
    component.push_back(e);
    return component.size()-1;
  }
  void completeFaces() {
    component.clear();
    fill(all(incident_face), -1);
    for (int e = 0; e < edges(); ++e) {
      if (incident_face[e] >= 0) continue;
      int f = newFace(e), x = e;
      do {
        incident_face[x] = f;
        x = next[x];
      } while (x != e);
    }
  }

  // assume connected
  vector<int> pfaffianOrientation() {
    // take any spanning tree T
    vector<int> dir(edges(), -2), seen(vertices());
    function<void(int)> dfs1 = [&](int u) {
      seen[u] = 1;
      int e = incident_edge[u];
      do {
        int v = origin[twin[e]];
        if (!seen[v]) {
          dir[e] = 1;
          dir[twin[e]] = -dir[e];
          dfs1(v);
        }
        e = next[twin[e]];
      } while (e != incident_edge[u]);
    };
    for (int u = 0; u < vertices(); ++u)
      if (!seen[u]) dfs1(u);

    // take any dual spanning tree that does not cross T
    seen = vector<int>(faces());
    vector<int> come(faces(), -1);
    function<void(int,int)> dfs2 = [&](int f, int p) {
      int parity = 0;
      int free_edge = -1;
      seen[f] = 1;
      int e = component[f];
      do {
        int g = incident_face[twin[e]];
        if (dir[e] == -2 && !seen[g]) dfs2(g, twin[e]);
        if (dir[e] == -2) { assert(free_edge == -1); free_edge = e; }
        else if (dir[e] == 1) ++parity;
        e = next[e];
      } while (e != component[f]);
      if (free_edge != -1) {
        dir[free_edge] = -(parity % 2 == 0 ? -1 : +1);
        dir[twin[free_edge]] = -dir[free_edge];
      }
    };
    dfs2(0, -1);
    return dir;
  }
  using Int = __int128_t;
  Int countPerfectMatching() {
    vector<int> dir = pfaffianOrientation();
    int n = vertices();
    vector<vector<Int>> A(n, vector<Int>(n));
    for (int e = 0; e < edges(); ++e) 
      A[origin[e]][origin[twin[e]]] = dir[e];

    // compute determinant
    Int det = 1;
    for (int j = 0; j < n; ++j) {
      for (int i = j+1; i < n; ++i) {
        while (A[i][j]) { 
          det = -det;
          Int t = A[j][j] / A[i][j];
          for (int k = j; k < n; ++k) 
            swap(A[i][k], A[j][k] -= t * A[i][k]);
        }
      }
      det *= A[j][j]; // % mod
    }
    return sqrt(det + 0.1);
  }
};

using Int = long long;
Int dominoCount(vector<vector<char>> table) {
  int m = table.size(), n = table[0].size();
  vector<vector<int>> index(m, vector<int>(n, -1));
  PlaneGraph g;
  unordered_map<int,unordered_map<int,int>> adj;
  for (int i = 0; i < m; ++i) {
    for (int j = 0; j < n; ++j) {
      if (table[i][j] == '.') {
        index[i][j] = g.newVertex(Point({i,j}));
      }
    }
  }
  unordered_map<int,int> next_inc, prev_inc;
  for (int i = 0; i < m; ++i) {
    for (int j = 0; j < n; ++j) {
      vector<int> inc;
      int x = index[i][j];
      int dx[] = {1,0,-1,0}, dy[] = {0,1,0,-1};
      for (int p = 0; p < 4; ++p) {
        int k = i+dx[p], l = j+dy[p];
        if (k < 0 || l < 0) continue;
        if (k >= table.size() || l >= table[k].size()) continue;
        if (table[i][j] != '.' || table[k][l] != '.') continue;
        int y = index[k][l];
        if (!adj[x].count(y)) adj[x][y] = g.newEdge(x);
        if (!adj[y].count(x)) adj[y][x] = g.newEdge(y);
        g.twin[adj[x][y]] = adj[y][x];
        g.twin[adj[y][x]] = adj[x][y];
        g.incident_edge[x] = adj[x][y];
        g.incident_edge[y] = adj[y][x];
        inc.push_back(adj[x][y]);
      }
      for (int i = 0; i < inc.size(); ++i) {
        int j = (i == inc.size()-1 ? 0 : i+1);
        next_inc[inc[i]] = inc[j];
        prev_inc[inc[j]] = inc[i];
      }
    }
  }
  for (int e = 0; e < g.edges(); ++e) {
    g.next[e] = prev_inc[g.twin[e]];
    g.prev[e] = g.twin[next_inc[e]];
  }
  g.completeFaces();
  return g.countPerfectMatching();
}

void SPOJ_GNY07H() {
  int ncase;
  cin >> ncase;
  for (int icase = 0; icase < ncase; ++icase) {
    int w;
    cin >> w;
    vector<vector<char>> table(4, vector<char>(w, '.')) ;
    cout << icase+1 << " " << dominoCount(table) << endl;
  }
}

int main() {
  SPOJ_GNY07H();
}

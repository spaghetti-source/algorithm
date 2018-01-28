//
// Modular Arithmetics
//
//  long long: 10^18 < 2^63-1 < 10^19 (strict inequality)
// __int128_t: 10^38 < 2^127-1 < 10^39
//
// g++ -std=c++17 -O3 -fmax-errors=1 -fsanitize=undefined
#include <bits/stdc++.h>
#pragma GCC optimize ("O3")

using namespace std;
#define fst first
#define snd second
#define all(c) begin(c), end(c)


// === tick a time ===
#include <ctime>
double tick() {
  static clock_t oldtick;
  clock_t newtick = clock();
  double diff = 1.0*(newtick - oldtick) / CLOCKS_PER_SEC;
  oldtick = newtick;
  return diff;
}


template <class T>
ostream &operator<<(ostream &os, const vector<T> &v) {
  os << "[";
  for (int i = 0; i < v.size(); os << v[i++]) 
    if (i > 0) os << " ";
  os << "]";
  return os;
}
template <class T>
ostream &operator<<(ostream &os, const vector<vector<T>> &v) {
  os << "[";
  for (int i = 0; i < v.size(); os << v[i++]) 
    if (i > 0) os << endl << " ";
  os << "]";
  return os;
}
using Int = long long;
struct ModInt {
  static Int mod; // Set this value
  Int val;
  ModInt(Int v=0) : val(v%mod) { }
  ModInt operator-() const { return ModInt(val?mod-val:val); }
  ModInt &operator+=(ModInt a) { 
    if ((val += a.val) >= mod) val -= mod;
    return *this;
  }
  ModInt &operator-=(ModInt a) { 
    if ((val -= a.val) < 0) val += mod;
    return *this;
  }
  ModInt &operator*=(ModInt a) { 
    val = (__uint128_t(val) * a.val) % mod;
    return *this;
  }
  ModInt &operator/=(ModInt a) { 
    Int u = 1, v = a.val, s = 0, t = mod;
    while (v) {
      Int q = t / v;
      swap(s -= u * q, u);
      swap(t -= v * q, v);
    }
    a.val = (s < 0 ? s + mod : s);
    val /= t;
    return (*this) *= a;
  }
  ModInt inv() const { return ModInt(1) /= (*this); }
  bool operator<(ModInt x) const { return val < x.val; }
};
Int ModInt::mod = 1e9+7;

ostream &operator<<(ostream &os, ModInt a) { os << a.val; return os; }
ModInt operator+(ModInt a, ModInt b) { return a += b; }
ModInt operator-(ModInt a, ModInt b) { return a -= b; }
ModInt operator*(ModInt a, ModInt b) { return a *= b; }
ModInt operator/(ModInt a, ModInt b) { return a /= b; }
ModInt pow(ModInt a, Int e) { 
  ModInt x(1);
  for (; e > 0; e /= 2) {
    if (e % 2 == 1) x *= a;
    a *= a;
  }
  return x;
}
ModInt stringToModInt(string s) {
  Int val = 0;
  for (int i = 0; i < s.size(); ++i) 
    val = (val*10 + (s[i]-'0')) % ModInt::mod;
  return ModInt(val);
}

// compute inv[1], inv[2], ..., inv[mod-1] in O(n) time
vector<ModInt> inverse() {
  Int mod = ModInt::mod;
  vector<ModInt> inv(mod);
  inv[1].val = 1;
  for (Int a = 2; a < mod; ++a) 
    inv[a] = inv[mod % a] * ModInt(mod - mod/a);
  return inv;
}


//
// Solve x^2 = n; mod should be a prime
// 
// Verified: Code Forces Quadratic Equations
//
bool isQuadraticResidue(ModInt n) {
  return n.val == 0 || n.mod == 2 || pow(n, (n.mod-1)/2).val == 1;
}
ModInt sqrt(ModInt n) { 
  const Int mod = ModInt::mod;
  if (n.val == 0 || mod == 2) return n;
  int M = __builtin_ctz(mod-1), Q = (mod-1)>>M;
  ModInt z(2);
  while (isQuadraticResidue(z)) ++z.val;
  ModInt c = pow(z, Q);
  ModInt t = pow(n, Q);
  ModInt R = pow(n, (Q+1)/2);
  while (t.val != 1) {
    int i = 0;
    for (ModInt s = t; s.val != 1; s *= s) ++i;
    if (M == i) exit(0);
    ModInt b = pow(c, 1<<(M-i-1));
    M = i;
    c = b*b;
    t *= c;
    R *= b;
  }
  return R;
}
vector<ModInt> quadraticEquation(ModInt a, ModInt b, ModInt c) {
  if (ModInt::mod == 2) {
    vector<ModInt> ans;
    if (c.val == 0) ans.push_back(c);
    if ((a + b + c).val == 0) ans.push_back(ModInt(1));
    return ans;
  } else {
    b /= (a+a); c /= a; 
    ModInt D = b*b - c;
    if (!isQuadraticResidue(D)) return {};
    ModInt s = sqrt(D), x = -b+s, y = -b-s;
    return (x.val < y.val) ? vector<ModInt>({x, y}) : 
           (x.val > y.val) ? vector<ModInt>({y, x}) : 
                             vector<ModInt>({x});
  }
}

// Discrete Logarithm by Shanks' Baby-Step Giant-Step
// Find k such that a^k == b 
Int log(ModInt a, ModInt b) {
  Int h = ceil(sqrt(ModInt::mod+1e-9));
  unordered_map<Int,Int> hash;
  ModInt x(1);
  for (Int i = 0; i < h; ++i) {
    if (!hash.count(x.val)) hash[x.val] = i;
    x *= a;
  }
  x = x.inv(); 
  ModInt y = b;
  for (int i = 0; i < h; ++i) {
    if (hash.count(y.val)) return i*h+hash[y.val];
    y *= x;
  }
  return -1;
}

struct ModMatrix {
  int m, n; // m times n matrix
  vector<vector<ModInt>> val;
  ModInt &operator()(int i, int j) { return val[i][j]; }
  ModMatrix(int m, int n) : 
    m(m), n(n), val(m, vector<ModInt>(n)) { }
  ModMatrix operator-() const {
    ModMatrix A(m, n);
    for (int i = 0; i < m; ++i) 
      for (int j = 0; j < n; ++j) 
        A.val[i][j] = -val[i][j];
    return A;
  }
  ModMatrix &operator+=(ModMatrix A) {
    for (int i = 0; i < m; ++i) 
      for (int j = 0; j < n; ++j) 
        val[i][j] += A.val[i][j];
    return *this;
  }
  ModMatrix &operator-=(ModMatrix A) {
    for (int i = 0; i < m; ++i) 
      for (int j = 0; j < n; ++j) 
        val[i][j] -= A.val[i][j];
    return *this;
  }
  ModMatrix &operator*=(ModMatrix A) {
    for (int i = 0; i < m; ++i) {
      vector<ModInt> row(A.n);
      for (int j = 0; j < A.n; ++j) {
        for (int k = 0; k < A.m; ++k) 
          row[j] += val[i][k] * A.val[k][j];
      }
      val[i] = row;
    }
    return *this;
  }
  static ModMatrix eye(int n) {
    ModMatrix I(n, n);
    for (int i = 0; i < n; ++i) I.val[i][i].val = 1;
    return I;
  }
  static ModMatrix zero(int n) {
    return ModMatrix(n, n);
  }
  // mod should be prime
  ModMatrix inv() const { 
    ModMatrix B = eye(n);
    vector<vector<ModInt>>  a = val;
    vector<vector<ModInt>> &b = B.val;
    for (int i = 0, j, k; i < n; ++i) {
      for (j = i; j < n && a[j][i].val == 0; ++j);
      if (j == n) return ModMatrix(0,0); // regularity is checked by m = 0
      swap(a[i], a[j]); 
      swap(b[i], b[j]);
      ModInt inv = a[i][i].inv();
      for (k = i; k < n; ++k) a[i][k] *= inv;
      for (k = 0; k < n; ++k) b[i][k] *= inv;
      for (j = 0; j < n; ++j) {
        if (i == j || a[j][i].val == 0) continue;
        ModInt c = a[j][i];
        for (k = i; k < n; ++k) a[j][k] -= c * a[i][k];
        for (k = 0; k < n; ++k) b[j][k] -= c * b[i][k];
      }
    }
    return B;
  }
  // It can be used for any composite modulo.
  ModInt det() const {
    vector<vector<ModInt>> a = val; 
    ModInt D(1);
    for (int j = 0; j < n; ++j) {
      for (int i = j+1; i < n; ++i) {
        while (a[i][j].val) { 
          D = -D;
          ModInt t(a[j][j].val/a[i][j].val);
          for (int k = j; k < n; ++k) 
            swap(a[i][k], a[j][k] -= t * a[i][k]);
        }
      }
      D *= a[j][j];
    }
    return D;
  }
};
ModMatrix operator+(ModMatrix A, ModMatrix B) { return A += B; }
ModMatrix operator-(ModMatrix A, ModMatrix B) { return A -= B; }
ModMatrix operator*(ModMatrix A, ModMatrix B) { return A *= B; }
ModMatrix pow(ModMatrix A, int k) {
  ModMatrix X = ModMatrix::eye(A.n);
  for (; k > 0; k /= 2) {
    if (k % 2 == 1) X *= A;
    A *= A;
  }
  return X;
}
ModInt dot(ModMatrix A, ModMatrix B) {
  ModInt val;
  for (int i = 0; i < A.m; ++i)
    for (int j = 0; j < A.n; ++j)
      val += A.val[i][j] * B.val[i][j];
  return val;
}
using ModVector = vector<ModInt>;
ModVector operator*(ModMatrix A, ModVector x) {
  vector<ModInt> y(A.m);
  for (int i = 0; i < A.m; ++i) 
    for (int j = 0; j < A.n; ++j)
      y[i] += A.val[i][j] * x[j];
  return y;
}
ModInt dot(ModVector a, ModVector b) {
  ModInt val;
  for (int i = 0; i < a.size(); ++i)
    val += a[i] * b[i];
  return val;
}

//
// Only available for prime modulos.
// If you want to compute multiple inverses, 
// use LU decomposition instead of computing the inverse.
//
struct LUDecomposition {
  int n;
  vector<int> pi;
  vector<vector<ModInt>> val;
  LUDecomposition(ModMatrix A) : n(A.n), val(A.val) {
    pi.resize(n+1);
    iota(all(pi), 0);
    for (int i = 0, j, k; i < n; ++i) {
      for (k = i; k < n; ++k) 
        if (val[k][i].val) break; 
      if (k == n) { pi[n] = -1; return; } // NG
      if (k != i) {
        swap(pi[i], pi[k]);
        swap(val[i], val[k]);
        ++pi[n];
      }
      for (j = i+1; j < n; ++j) {
        if (val[j][i].val == 0) continue;
        val[j][i]/= val[i][i];
        for (k = i+1; k < n; ++k) 
          val[j][k] -= val[j][i] * val[i][k];
      }
    }
  }
  bool isRegular() const { return pi[n] >= 0; }
  ModVector solve(ModVector b) {
    vector<ModInt> x(b.size());
    for (int i = 0; i < n; ++i) {
      x[i] = b[pi[i]];
      for (int k = 0; k < i; ++k) 
        x[i] -= val[i][k] * x[k];
    }
    for (int i = n-1; i >= 0; --i) {
      for (int k = i+1; k < n; ++k)
        x[i] -= val[i][k] * x[k];
      x[i] /= val[i][i];
    }
    return x;
  }
  ModMatrix inverse() { // do not compute the inverse
    ModMatrix B(n, n);
    for (int j = 0; j < n; ++j) {
      for (int i = 0; i < n; ++i) {
        if (pi[i] == j) B.val[i][j].val = 1;
        for (int k = 0; k < i; k++)
          B.val[i][j] -= val[i][k] * B.val[k][j];
      }
      for (int i = n-1; i >= 0; --i) {
        for (int k = i+1; k < n; ++k) 
          B.val[i][j] -= val[i][k] * B.val[k][j];
        B.val[i][j] /= val[i][i];
      }
    }
    return B;
  }
  ModInt det() {
    ModInt D = val[0][0];
    for (int i = 1; i < n; i++)
      D *= val[i][i];
    return ((pi[n] - n) % 2 != 0) ? -D : D;
  }
};

void mulTest() {
  ModInt::mod = 1e9+7;
  int m = 4, n = 4;
  ModMatrix A(m,n);
  ModMatrix B(m,n);
  vector<ModInt> b(n);
  for (int i = 0; i < m; ++i) {
    for (int j = 0; j < n; ++j) {
      A(i,j).val = rand() % ModInt::mod;
      B(i,j).val = rand() % ModInt::mod;
    }
    b[i].val = rand() % ModInt::mod;
  }
  ModMatrix C = A * B;

  for (int i = 0; i < C.m; ++i) {
    for (int j = 0; j < C.n; ++j) {
      cout << C(i,j) << " ";
    }
    cout << endl;
  }
  cout << C.det() << endl;

  LUDecomposition LU(C);
  cout << LU.det() << endl;

  // A^{-1} b = x
  auto x = LU.solve(b);
  cout << b << " " << (C * x) << endl;
  cout << "end" << endl;
}

void CF_QUADRATIC_EQUATIONS() {
  int ncase; cin >> ncase;
  for (int icase = 0; icase < ncase; ++icase) {
    Int a, b, c, p;
    cin >> a >> b >> c >> p;
    ModInt::mod = p;
    vector<ModInt> ans = quadraticEquation(
      ModInt(a), ModInt(b), ModInt(c));
    cout << ans.size();
    for (int i = 0; i < ans.size(); ++i) 
      cout << " " << ans[i].val;
    cout << endl;
  }
}

void CF_DISCLOG() {
  ModInt::mod = 999998999999;
  ModInt a(21309), b(696969);
  cout << log(a, b) << endl;
}

int SPOJ_MIFF() {
  for (int icase = 0; ; ++icase) {
    int n, p; scanf("%d %d", &n, &p);
    if (n == 0 && p == 0) break;
    if (icase > 0) printf("\n");
    ModInt::mod = p;
    ModMatrix A(n, n);
    for (int i = 0; i < n; ++i)
      for (int j = 0; j < n; ++j)
        scanf("%d", &A(i,j));

    ModMatrix B = A.inv();
    if (B.m > 0) {
      for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
          printf("%d ", B(i,j));
        }
        printf("\n");
      }
    } else {
      printf("singular\n");
    }
  }
}

int main() {
  SPOJ_MIFF();
  //CF_DISCLOG();
  //CF_QUADRATIC_EQUATIONS();
  //mulTest();
}

// 
// Fast Fourier Transformation
//
// Description:
//   Given a complex sequence a[0,n), where n is a power of two.
//   Compute
//     A[k] = sum_k a[k] E^k 
//   where E = exp(2 pi i / n).
//
// Algorithm:
//   Cooley-Turkey's algorithm.
//
// Complexity:
//   O(n log n).
//
// Verified:
//   SPOJ235.

#include <iostream>
#include <vector>
#include <string>
#include <cmath>
#include <algorithm>
#include <complex>
#include <cstring>

using namespace std;

typedef complex<double> C;
void fft(vector<C> &a, int sign = 1) {
  int n = a.size(); // n should be a power of two
  double theta = 8 * sign * atan(1.0) / n; 
  for (int i = 0, j = 1; j < n - 1; ++j) {
    for (int k = n >> 1; k > (i ^= k); k >>= 1);
    if (j < i) swap(a[i], a[j]);
  }
  for (int m, mh = 1; (m = mh << 1) <= n; mh = m) {
    int irev = 0;
    for (int i = 0; i < n; i += m) {
      C w = exp(C(0, theta*irev));
      for (int k = n >> 2; k > (irev ^= k); k >>= 1);
      for (int j = i; j < mh + i; ++j) {
        int k = j + mh;
        C x = a[j] - a[k];
        a[j] += a[k];
        a[k] = w * x;
      }
    }
  }
}


const int WIDTH = 5;
const long long RADIX = 100000; // = 10^WIDTH

vector<C> parse(const char s[]) {
  int n = strlen(s);
  int m = (n + WIDTH-1) / WIDTH;
  vector<C> v(m);
  for (int i = 0; i < m; ++i) {
    int b = n - WIDTH * i, x = 0;
    for (int a = max(0, b - WIDTH); a < b; ++a)
      x = x * 10 + s[a] - '0';
    v[i] = x;
  }
  return v;
}

void print(const vector<C> &v) {
  int i, N = v.size();
  vector<long long> digits(N + 1, 0);
  long double err = 0;

  for (i = 0; i < N; i++) {
    digits[i] = (long long)(v[i].real() + 0.5);
  }
  long long c = 0;
  for (i = 0; i < N; i++) {
    c += digits[i];
    digits[i] = c % RADIX;
    c /= RADIX;
  }
  for (i = N-1; i > 0 && digits[i] == 0; i--);
  printf("%lld", digits[i]);
  for (i--; i >= 0; i--)
    printf("%.*lld", WIDTH, digits[i]);
  printf("\n");
}

char a[310000], b[310000];
int main() {

  int T; scanf("%d", &T);
  while (T--) {
    scanf("%s %s", a, b);
    vector<C> A = parse(a);
    vector<C> B = parse(b);

    int N = 1;
    while (N < max(A.size(), B.size())) N *= 2;
    N *= 2;
    A.resize(N);
    B.resize(N);

    fft(A, +1);
    fft(B, +1);
    for (int i = 0; i < N; i++) A[i] *= B[i];
    fft(A, -1);
    for (int i = 0; i < N; i++) A[i] /= N;

    print(A);
  }
}

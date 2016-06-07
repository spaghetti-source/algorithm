//
// Coordinate Compression
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
#define TEST(s) if (!(s)) { cout << __LINE__ << " " << #s << endl; exit(-1); }

template <class T>
struct coordinate_compression {
  vector<T> xs, ys; // O(n log n)
  coordinate_compression(vector<T> xs) : xs(xs), ys(xs) {
    sort(all(ys));
    ys.erase(unique(all(ys)), ys.end());
  }
  int index(T a) const { // O(log n)
    auto it = lower_bound(all(ys), a);
    if (it == ys.end() || *it != a) return -1;
    return distance(ys.begin(), it);
  }
  int value(int k) const { return ys[k]; } // k in [0, ys.size())
  int size() const { return ys.size(); }
};


int main() {
  vector<int> x = {3,1,4,1,5,9};
  coordinate_compression<int> compressor(x);
  for (int i = 0; i <x.size(); ++i)
    cout << compressor.index(x[i]) << " ";
  cout << endl;
}

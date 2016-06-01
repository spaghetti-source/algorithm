//
// 0 0 0
// 1 0 0
// 2 0 0
// 0 1 0
// 1 1 0
// 2 1 0
// ...
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

template <class It>
bool next_radix(It begin, It end, int base) {
  for (It cur = begin; cur != end; ++cur) {
    if ((*cur += 1) >= base) *cur = 0;
    else return true;
  }
  return false;
}

int main() {
  vector<int> a(3, 0);
  do {
    for (int i = 0; i < 3; ++i) cout << a[i];
    cout << endl;
  } while (next_radix(all(a), 2));
}

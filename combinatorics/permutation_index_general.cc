// 
// Lexicographical index of a general permutation
//
// Description:
//   index_perm computes the lexicographical index of a permutation x,
//   i.e., it gives k such that x = next_permutation^k (sort(x)).
//   Note that x can have duplicated elements.
//
//   unindex_perm is the inverse function of index_perm.
//
// Algorithm:
//   index(x) = 
//     id = 0
//     for each b < head(x)
//       id += number of permutations start from b
//     return id + index(tail(x))
// 
// Complexity:
//   O(n |A|) time, O(n) space,
//   where |A| is the number of distinct elements in x
//   If A has no distinct elements, use O(n log n) algo.
//

#include <iostream>
#include <vector>
#include <map>
#include <cstdio>
#include <algorithm>
#include <functional>

using namespace std;

#define fst first
#define snd second
#define all(c) ((c).begin()), ((c).end())
#define TEST(s) if (!(s)) { cout << __LINE__ << " " << #s << endl; exit(-1); }

typedef unsigned long long ull;
template <class T>
ull index_perm(vector<T> x) {
  ull suff = 1; 
  map<T, ull> count;
  for (int i = 0; i < x.size(); ++i) {
    suff *= (i + 1);
    suff /= (++count[x[i]]);
  }
  ull id = 0;
  for (int i = 0; i < x.size(); ++i) {
    for (auto &a: count) {
      suff = suff * a.snd-- / (x.size() - i);
      if (a.fst == x[i]) {
        if (a.snd == 0) count.erase(a.fst);
        break;
      }
      id += suff;
      suff = suff * (x.size() - i) / ++a.snd;
    }
  }
  return id;
}
template <class T>
vector<T> unindex_perm(ull id, vector<T> x) {
  ull suff = 1; 
  map<T, ull> count;
  for (int i = 0; i < x.size(); ++i) {
    suff *= (i + 1);
    suff /= (++count[x[i]]);
  }
  for (int i = 0; i < x.size(); ++i) {
    for (auto &a: count) {
      suff = suff * a.snd-- / (x.size() - i);
      if (id < suff) {
        x[i] = a.fst;
        if (a.snd == 0) count.erase(a.fst);
        break;
      }
      id -= suff;
      suff = suff * (x.size() - i) / ++a.snd;
    }
  }
  return x;
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

bool TEST_index_perm() {
  for (int iter = 0; iter < 1000; ++iter) { 
    int n = 10;
    vector<int> x(n), c = {0,1,9};
    for (int i = 0; i < n; ++i)
      x[i] = c[rand() % c.size()];
    int idx = index_perm(x);
    auto y = x;
    sort(all(y));
    int it = 0;
    do {
      if (it++ == idx && x != y) return false;
    } while (next_permutation(all(y)));
  }
  return true;
}

bool TEST_unindex_perm() {
  for (int iter = 0; iter < 1000; ++iter) { 
    int n = 10;
    vector<int> x(n), c = {0,1,9};
    for (int i = 0; i < n; ++i)
      x[i] = c[rand() % c.size()];
    int idx = index_perm(x);
    auto shuffled_x = x;
    random_shuffle(all(shuffled_x));
    auto y = unindex_perm(idx, shuffled_x);
    if (x != y) return false;
  }
  return true;
}

int main() {
  TEST(TEST_index_perm());
  TEST(TEST_unindex_perm());
}

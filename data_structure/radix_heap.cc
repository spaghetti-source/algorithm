#include <iostream>
#include <cstring>
#include <vector>
#include <cstdio>
#include <algorithm>
#include <functional>

using namespace std;

#define fst first
#define snd second
#define all(c) ((c).begin()), ((c).end())
#define TEST(s) if (!(s)) { cout << __LINE__ << " " << #s << endl; exit(-1); }

// min heap
template <class T>
struct radix_heap {
  typedef unsigned int uint;
  static int bsr(uint a) { return a ? 31 - __builtin_clz(a) : -1; }
  uint size, last;
  vector<pair<uint, T>> v[33];
  radix_heap() : size(0), last(0) { }

  bool empty() const { return size == 0; }
  void aux(const pair<uint, T> &p) { v[bsr(p.fst^last)+1].push_back(p); }
  pair<uint, T> top() {
    if (v[0].empty()) {
      int i = 1;
      while (v[i].empty()) ++i;
      last = min_element(all(v[i]))->fst;
      for (auto &p: v[i]) aux(p);
      v[i].clear();
    }
    return v[0].back();
  }
  void push(uint key, T value) { ++size; aux({key, value}); }
  void pop() { --size; top(); v[0].pop_back(); }
};

int main() {
  radix_heap<string> heap;
  heap.push(strlen("test"), "test");
  heap.push(strlen("a"), "a");
  heap.push(strlen("ab"), "ab");
  heap.push(strlen("aaa"), "aaa");
  heap.push(strlen("xyzz"), "xyzz");
  while (!heap.empty()) {
    cout << heap.top().snd << endl;
    heap.pop();
  }
}

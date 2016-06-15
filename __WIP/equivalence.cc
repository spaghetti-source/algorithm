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

// specify operator == 
template <class T, class F>
void equivalence(vector<T> x, F eq) {
  vector<int> parent(x.size());
  for (int i = ; i < x.size(); ++i) {
    parent[i] = i;
    for (int j = 0; j < i; ++j) {
      parent[j] = parent[parent[j]];
      if (eq(x[i], x[j])) parent[parent[parent[j]]] = i;
    }
  }
  for (int i = 0; i < x.size(); ++i) parent[i] = parent[parent[i]];

  for (int i = 0; i < parent.size(); ++i) 
    printf("%d %d \n", x[i], parent[i]);
}

int main() {
  vector<int> x = {3,1,4,1,5,9,2,6,5,3,5,8,9};
  equivalence(x, [&](int a, int b) { return a == b; });
}

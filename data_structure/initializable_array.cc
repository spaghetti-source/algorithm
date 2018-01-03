//
// Initializable Array
//
// Description:
//   It allows the following operations in O(1) time.
//     - init(a):    initialize xs[i] = a for all i
//     - xs[i]:      return xs[i]
//     - set(i, a):  set xs[i] = a
//   The important operation is "init", which usually
//   requires O(n) time. By maintaining timestamps,
//   we can "emulate" the initialization.
//
// Complexity:
//   O(1).
//
// References:
//   J. Bentley (1986): Programming pearls. Addison-Wesley.
//
#include <bits/stdc++.h>
using namespace std;

#define fst first
#define snd second
#define all(c) ((c).begin()), ((c).end())

template <class T>
struct InitializableArray {
  T initv, *value;
  size_t b, *from, *to;
  InitializableArray(int n) {
    value = new T[n];
    from  = new size_t[n];
    to    = new size_t[n];
  }
  bool chain(int i) {
    int j = from[i];
    return j < b && to[j] == i;
  }
  void init(T a) { 
    initv = a;
    b = 0;
  }
  T operator[](int i) { 
    return chain(i) ? value[i] : initv;
  }
  void set(int i, T a) { 
    if (!chain(i)) {
      from[i] = b;
      to[b++] = i;
    }
    value[i] = a;
  }
};

int main() {
  InitializableArray<int> a(3);
  cout << a.value[0] << endl;
  a.init(0);
  a.init(2);
  a.set(1, 5);
  cout << a[0] << endl;
  cout << a[1] << endl;
  cout << a[2] << endl;
  a.set(2, 3);
  a.init(0);
}

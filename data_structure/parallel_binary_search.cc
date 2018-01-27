// 
// Parallel Binary Search
//
// Description:
//
//   Let W(t) be a data structure depending on time t.
//   Suppose that there are n agents. Each agent j wants to
//   find the smallest time t(j) such that cond(j, W(t(j))) == true
//   where cond(j, W(t)) is monotone in t.
//
//   If we perform n binary searches independently, we will construct 
//   W(t) multiple times. Thus, we avoid the redundant construction
//   by performing n binary searches in parallel. Imagine a binary 
//   search tree on t. For each node, we first construct W(t). 
//   Then we process multiple agents in parallel. Then, the total 
//   number of constructions is O(log T), which is independent to 
//   the number of agents.
//
// Complexity:
//
//   Suppose that W(t) is constructed from W(t') in time M(|t-t'|),
//   and the condition cond(j,W(t)) is evaluated in time Q.
//   Then, it runs in O(M(T log T) + n Q log T) time.
//
//   Even if W does not have decremental operations, i.e., W(t) 
//   cannot be constructed from W(t') with t' > t, we can still use
//   the parallel binary search that runs in O(M(T) log T + n Q log T);
//   time. The constant factor is twice worse than the above.
//   Similar result holds if W does not have incremental operations.
//
#include <bits/stdc++.h>

using namespace std;

#define fst first
#define snd second
#define all(c) ((c).begin()), ((c).end())
#define TEST(s) if (!(s)) { cout << __LINE__ << " " << #s << endl; exit(-1); }

using Int = __int128_t;

// Point Query, Range Update 
template <class T>
struct FenwickTree {
  vector<T> x;
  FenwickTree(int n) : x(n) { }
  void add(int k, T a) { // aux
    for (; k < x.size(); k |= k+1) x[k] += a;
  }
  void add(int i, int j, T a) { // add x[k] += a for all k in [i,j]
    add(i, a);
    if (j+1 < x.size()) add(j+1, -a);
  }
  T get(int k) { // return x[k]
    T sum = 0;
    for (; k >= 0; k = (k&(k+1))-1) sum += x[k];
    return sum;
  }
};

template <class Update, class Cond>
vector<int> parallelBinarySearch(
  int n, int lo, int hi, Update update, Cond cond) {
  using It = vector<int>::iterator;
  vector<int> agents(n), solution(n, lo);
  iota(all(agents), 0);

  It begin = agents.begin(), end = agents.end();
  deque<tuple<int,int,It,It>> stack = {make_tuple(lo, hi, begin, end)};
  while (!stack.empty()) {
    // invariant: elems in [begin, end) satisfy "!cond(lo) and cond(hi)"
    tie(lo, hi, begin, end) = stack.back();
    stack.pop_back();

    if (begin == end) continue;
    if (lo+1 == hi) {
      for_each(begin, end, [&](int k) { solution[k] = hi; });
      continue;
    }
    int mi = (lo + hi) / 2;
    update(mi);
    It mid = partition(begin, end, [&](int k) { return cond(k); });
    stack.push_back(make_tuple(mi, hi, mid, end));
    stack.push_back(make_tuple(lo, mi, begin, mid));
  }
  return solution;
}

void SPOJ_METEORS() {
  int n, m, k;
  scanf("%d %d", &n, &m);
  vector<vector<int>> S(n); 
  vector<Int> p(n);
  for (int i = 0; i < m; ++i) {
    int j; scanf("%d", &j);
    S[j-1].push_back(i);
  }
  for (int i = 0; i < n; ++i)
    scanf("%lld", &p[i]);
  scanf("%d", &k);
  vector<int> l(k), r(k);
  vector<Int> a(k);
  for (int i = 0; i < k; ++i) {
    scanf("%d %d %lld", &l[i], &r[i], &a[i]);
    --l[i]; --r[i];
  }

  FenwickTree<Int> FT(m);
  int curr = -1; 
  auto update = [&](int time) {
    while (curr < time) {
      ++curr;
      if (l[curr] <= r[curr]) {
        FT.add(l[curr], r[curr], a[curr]);
      } else {
        FT.add(l[curr], m-1, a[curr]);
        FT.add(0, r[curr], a[curr]);
      }
    }
    while (curr > time) {
      if (l[curr] <= r[curr]) {
        FT.add(l[curr], r[curr], -a[curr]);
      } else {
        FT.add(l[curr], m-1, -a[curr]);
        FT.add(0, r[curr], -a[curr]);
      }
      --curr;
    }
  };
  // minimum time such that cond = true.
  auto cond = [&](int j) {
    Int total = 0;
    for (int i: S[j]) {
      total += FT.get(i);
    }
    return total >= p[j];
  };

  auto solution = parallelBinarySearch(n, -1, k, update, cond);

  for (int i = 0; i < n; ++i) {
    if (solution[i] >= k) cout << "NIE" << endl;
    else cout << 1+solution[i] << endl;
  }
}

int main() {
  SPOJ_METEORS();
}

// 
// Min-Map Heap (four heaps technique)
//
// Description:
//   A data structure for push/min/max/popmin/popmax.
//
// Algorithm:
//   Maintain two min heaps (minh, minp) an two max heaps (maxh, maxp).
//   The actual elements in the heap is (minh \cup maxh) \setminus (\minp \cup maxp).
//   
// Complexity:
//   Amortized O(1) for push and top, and O(log n) for pop. 
//

template <class T, class L = less<T>, class G = greater<T>>
struct minmax_heap {
  priority_queue<T, vector<T>, G> minh, minp;
  priority_queue<T, vector<T>, L> maxh, maxp;
  void normalize() {
    while (!minp.empty() && minp.top() == minh.top()) {
      minp.pop();
      minh.pop();
    }
    while (!maxp.empty() && maxp.top() == maxh.top()) {
      maxp.pop();
      maxh.pop();
    }
  }
  void push(T x) { minh.push(x); maxh.push(x); }
  T min() { normalize(); return minh.top(); }
  T max() { normalize(); return maxh.top(); }
  void pop_min() { normalize(); maxp.push(minh.top()); minh.pop(); }
  void pop_max() { normalize(); minp.push(maxh.top()); maxh.pop(); }
};

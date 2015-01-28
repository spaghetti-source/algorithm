struct union_find {
  int components;
  vector<int> root, next, size;
  void clear(int n) {
    components = n;
    root.resize(n);
    iota(root.begin(), root.end(), 0);
    next.assign(n, -1);
    size.assign(n,  1);
  }
  bool find(int x, int y) { return root[x] == root[y]; }
  bool unite(int x, int y) {
    if ((x = root[x]) == (y = root[y])) return false;
    if (size[x] < size[y]) swap(x, y);
    size[x] += size[y];
    for (int z; y >= 0; y = z) {
      z = next[y];
      next[y] = next[x];
      next[x] = y;
      root[y] = x;
    }
    --components;
    return true;
  }
};

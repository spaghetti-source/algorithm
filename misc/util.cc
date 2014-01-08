template <class T>
struct eictionary {
  unordered_map<T, size_t> dict;
  vector<T> idict;
  size_t id(T s) { 
    if (!dict.count(s)) {
      dict[s] = idict.size();
      idict.push_back(s);
    }
    return dict[s];
  }
  T value(size_t id) {
    return idict[id];
  }
};


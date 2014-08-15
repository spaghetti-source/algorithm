#pragma once
#include <sstream>
#include <string>
#include <unordered_map>
#include <cstdio>

namespace parse_args {
  std::unordered_map<std::string, std::string> args;
  void init(int argc, char *argv[]) {
    for (int i = 1; argv[i]; ++i) {
      char cmd[256], param[256];
      std::sscanf(argv[i], "%[^=]=%s", cmd, param);
      args[cmd] = param;
    }
  }
  bool has(std::string s) {
    return args.count(s);
  }
  template <class T>
  T get(std::string s, T x = T()) {
    std::stringstream ss(args[s]);
    T value;
    ss >> value;
    return value;
  }
};

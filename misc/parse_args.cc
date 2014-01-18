struct parse_args {
  unordered_map<string, string> args;
  parse_args(int argc, char *argv[]) {
    for (int i = 1; argv[i]; ++i) {
      char cmd[256], param[256];
      sscanf(argv[i], "%[^=]=%s", cmd, param);
      args[string(cmd)] = string(param);
    }
  }
  bool has(string s) {
    return args.count(s);
  }
  template <class T>
  T get(string s) {
    stringstream ss(args[s]);
    T value;
    ss >> value;
    return value;
  }
};

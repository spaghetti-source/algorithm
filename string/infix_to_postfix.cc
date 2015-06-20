//
// Convert infix notation to postfix notation
//
// Description:
//   1 * 2 + 3 * (4 + 5) <- infix notation
//   1 2 * 3 4 5 + * +   <- postfix notation
//   Postfix notation is easy to evaluate.
//
// Algorithm:
//   Shunting-yard algorithm by Dijkstra.
//   
// Verified:
//   SPOJ 4: ONP - Transform the Expression
//

#include <iostream>
#include <stack>
#include <vector>
#include <cstdio>
#include <algorithm>
#include <functional>
#include <sstream>

using namespace std;

string infix_to_postfix(string s) {
  s += '_'; // terminal symbol
  stringstream ss;
  vector<char> op = {'_'};
  auto rank = [](char c) { return string("(^/*-+)").find(c); };
  for (char c: s) {
    if (isalnum(c)) ss << c; 
    else {
      for (; op.back() != '('; op.pop_back()) {
        if (rank(op.back()) >= rank(c)) break;
        ss << op.back();
      }
      if (c == ')') op.pop_back();
      else          op.push_back(c);
    }
  }
  return ss.str();
}

int main() {
  int ncase; scanf("%d", &ncase);
  for (int icase = 0; icase < ncase; ++icase) {
    char s[1024]; scanf("%s", s);
    printf("%s\n", infix_to_postfix(s).c_str());
  }
}

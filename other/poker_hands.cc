// 
// Poker Hands 
//
// Description:
//   It determines the hand of the poker.
//
// Algorithm:
//   Naive.
//
// Verified:
//   UVA 10315	Poker Hands (
//
#include <iostream>
#include <string>
#include <vector>
#include <cstdio>
#include <algorithm>
#include <functional>

using namespace std;

#define fst first
#define snd second
#define all(c) ((c).begin()), ((c).end())

struct card {
  int rank, suit;
  card(int rank, int suit) : rank(rank), suit(suit) { }
  card(char s[]) { // s is, i.e., "7H"
    for (rank = 0; s[0] != "_A23456789TJQK"[rank]; ++rank);
    for (suit = 0; s[1] != "CDHS"[suit]; ++suit);
  }
};

enum {
  HIGHEST_CARD, ONE_PAIR, TWO_PAIR, THREE_OF_A_KIND,
  STRAIGHT, FLUSH, FULL_HOUSE, FOUR_OF_A_KIND, STRAIGHT_FLUSH };
pair<int,vector<int>> poker_hand(vector<card> cs) {
  vector<int> code, freq(30);
  for (card c: cs) freq[c.rank+13] = freq[c.rank] += 1;
  for (int i = 14; i >= 2; --i) 
    if (freq[i]) code.push_back(i);
  stable_sort(all(code), [&](int i, int j) { return freq[i] > freq[j]; });

  bool straight = false, flush = true;
  for (int i = 1, j; i <= 10; ++i) { // beginning of straight (10 == 'X')
    for (j = 0; j < 5 && freq[i+j]; ++j);
    if (j == 5) straight = true;
  }
  for (int i = 1; i < 5; ++i)
    if (cs[0].suit != cs[i].suit) flush = false;

  if (straight && flush) return {STRAIGHT_FLUSH, code};
  if (freq[code[0]] == 4) return {FOUR_OF_A_KIND, code};
  if (freq[code[0]] == 3 && freq[code[1]] == 2) return {FULL_HOUSE, code};
  if (flush) return {FLUSH, code};
  if (straight) return {STRAIGHT, code};
  if (freq[code[0]] == 3) return {THREE_OF_A_KIND, code};
  if (freq[code[0]] == 2 && freq[code[1]] == 2) return {TWO_PAIR, code};
  if (freq[code[0]] == 2) return {ONE_PAIR, code};
  return {HIGHEST_CARD, code};
};

int main() {
  while (1) {
    vector<card> cs[2];
    pair<int, vector<int>> res[2];
    for (int k = 0; k < 2; ++k) {
      for (int i = 0; i < 5; ++i) {
        char s[120];
        if (scanf("%s", s) != 1) return 0;;
        cs[k].push_back(card(s));
      }
      res[k] = poker_hand(cs[k]);
    }
    if      (res[0] < res[1]) printf("White wins.\n");
    else if (res[0] > res[1]) printf("Black wins.\n");
    else                      printf("Tie.\n");
  }
}

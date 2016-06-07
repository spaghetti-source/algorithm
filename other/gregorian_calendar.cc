// 
// Calendar
//
// Description:
//   A code for converting gregorian date <=> julian day number.
//

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

struct gregorian_date { 
  int y, m, d; 
  gregorian_date(int y, int m, int d) : y(y), m(m), d(d) { }

  // from julian day number to gregorian date
  gregorian_date(int x) {
    int e = 4*x + 4*((((4*x+274277)/146097)*3)/4) + 5455;
    int h = 5*((e% 1461) / 4) + 2;
    d = (h%153)/5 + 1;
    m = (h/153+2)%12 + 1;
    y = (e/1461) + (14-m)/12 - 4716;
  }
  // number of days from the epoch
  int julian_day_number() const {
    int a = (14-m)/12, Y = y+4800-a, M = m+12*a-3;
    return d + (153*M+2)/5 + 365*Y + (Y/4) - (Y/100) + (Y/400) - 32045;
  }
  // week of the day
  int day_of_week() const {
    return (julian_day_number() + 1) % 7;
  }
  bool is_leap() const {
    return (y % 4 == 0 && y % 100 != 0) || y % 400 == 0;
  }
};
ostream &operator<<(ostream &os, const gregorian_date &g) {
  os << g.y << "/" << g.m << "/" << g.d;
  return os;
}

int main() {
  int JDN = gregorian_date(2000, 1, 1).julian_day_number();
  cout << JDN << endl;
  auto g = gregorian_date(JDN);
  cout << g << endl;
}

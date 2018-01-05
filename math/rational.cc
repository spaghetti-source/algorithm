//
// Rational Type
//
#include <bits/stdc++.h>
using namespace std;

#define fst first
#define snd second
#define all(c) ((c).begin()), ((c).end())

using Int = long long;
struct Rational {
  Int num, den; // x = num/den

  static Int gcd(Int a, Int b) {
    for (; a; swap(a, b %= a));
    return b;
  }
  Rational(Int a = 0) : num(a), den(1) { }
  Rational(Int a, Int b, bool do_normalize = false) : num(a), den(b) { 
    if (do_normalize) {
      auto g = gcd(num, den);
      num /= g; den /= g;
      if (den < 0) { num = -num; den = -den; }
    }
  }
  static Rational inf() { return Rational(1,0); }

  Rational inv() const { 
    if (num < 0) return {-den, -num};
    else         return { den,  num};
  }
  Rational operator+() const { return *this; }
  Rational operator-() const { return Rational(-num,den); }
  Rational &operator+=(Rational x) { 
    auto g = gcd(den, x.den);
    num = num*(x.den/g) + (den/g)*x.num;
    den = den*(x.den/g);
    return *this;
  }
  Rational &operator*=(Rational x) {
    auto g = gcd(num, x.den), h = gcd(den, x.num);
    num = (num/g)*(x.num/h);
    den = (den/h)*(x.den/g);
    if (den < 0) { num = -num; den = -den; }
    return *this;
  }
  Rational &operator-=(Rational x) { return *this += -x; }
  Rational &operator/=(Rational x) { return *this *= x.inv(); }
};
Rational operator+(Rational x, Rational y) { return x += y; }
Rational operator-(Rational x, Rational y) { return x -= y; }
Rational operator*(Rational x, Rational y) { return x *= y; }
Rational operator/(Rational x, Rational y) { return x /= y; }
int compare(Rational x, Rational y) { // sign(x-y)
  if (x.num == 0) return (y.num < 0) - (y.num > 0);
  if (x.num <  0) return y.num >= 0 ? -1 : compare(-y, -x);
  if (y.num <= 0) return 1;
  while (x.den != 0 && y.den != 0) {
    auto a = x.num/x.den, b = y.num/y.den;
    if (a != b) return a - b;
    swap(x.num -= a*x.den, x.den);
    swap(y.num -= b*y.den, y.den);
    swap(x, y);
  }
  if (x.den != 0) return -y.num;
  if (y.den != 0) return  x.num;
  return x.num - y.num;
}
bool operator==(Rational x, Rational y) { return compare(x,y)==0; }
bool operator!=(Rational x, Rational y) { return compare(x,y)!=0; }
bool operator<=(Rational x, Rational y) { return compare(x,y)<=0; }
bool operator>=(Rational x, Rational y) { return compare(x,y)>=0; }
bool operator<(Rational x, Rational y) { return compare(x,y)<0; }
bool operator>(Rational x, Rational y) { return compare(x,y)>0; }

ostream &operator<<(ostream &os, Rational x) { os<<x.num<<"/"<<x.den; return os; }

// Do not provide a cast operator since it makes debug difficult
using Real = double;
Real toReal(Rational x) { return (Real)x.num/x.den; }
Rational toRational(Real a) {
  Int num = 1, den = 0, as[20] = {0};
  auto b = a;
  for (int k = 0; k < 20; ++k) {
    auto num_ = num, den_ = den;
    num = 1; den = 0;
    as[k] = round(b);
    for (int i = k; i >= 0; --i) 
      swap(num, den += as[i] * num);
    if (abs(num) >= 1e9 || abs(den) >= 1e9) {
      num = num_; den = den_;
      break; // overflow; backup and forgive
    }
    auto error = fabs(Real(num)/den - a);
    b -= as[k];
    if (b == 0 || error < 1e-16) break;
    b = 1 / b;
  }
  if (den < 0) { num = -num; den = -den; }
  return Rational(num, den, false);
}


/*
Rational approx(Rational::Real a) { 
  Rational::Int num = 1, den = 0;
  Rational::Int as[30] = {0};
  auto b = a;
  for (int k = 0; k < 30; ++k) {
    as[k] = round(b);
    auto numk = num, denk = den;
    for (int i = k; i >= 0; --i) 
      swap(numk, denk += as[i] * numk);
    if (abs(numk) >= 1e9 || abs(denk) >= 1e9) break; // overflow 
    auto error = fabs(Real(numk)/denk - a);
    if (error < 1e-9) {
      num = numk; den = denk;
      break;
    }
    b -= as[k];
    if (b == 0) break;
    b = 1 / b;
  }
  return Rational(num, den);
}
*/
void verify_compare() {
  srand( time(0) );
  for (int iter = 0; iter < 100000; ++iter) {
    cout << iter << endl;
    Rational x(rand()%100, 1 + rand()%100);
    Rational y(rand()%100, 1 + rand()%100);
    assert( (x < y) == (-y < -x) );
    assert(  x ==  x );
    assert( -x == -x );
    assert( (x == y) == (-x == -y) );
    assert( (x != y) == (-x != -y) );
    if (x != Rational(0,1)) assert( (y + x) != (y - x) );
    if (x > Rational(0,1)) assert(  (y-x) <  (y+x) );
    if (x > Rational(0,1)) assert( (-y-x) < (-y+x) );
  }
}

int main() {
  // 0/1 < 1/1
  /*
  Rational x(2,1), y(3,1);
  cout << x + y << endl;
  */
  /*
  Rational x = atan(sqrt(Rational(1,2)));
  cout << x << endl;
  printf("%.12f\n%.12f\n", (double)x, atan(sqrt(1.0/2.0)));
  */
  verify_compare();
}

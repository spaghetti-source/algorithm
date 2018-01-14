// 
// Automatic Differentiation by Dual Numbers
//
// Description:
//
//   Dual number is an extended real number of the form a + epsilon b,
//   where epsilon is the first order infinitesimal (i.e, epsilon^2 = 0).
//   By overloading each function for dual numbers, as in the code, 
//   we can obtain the derivative of f at a by evaluating f(a + epsilon).
//
// Complexity:
// 
//   Linear in composition depth.
//
#include <bits/stdc++.h>

using namespace std;

#define fst first
#define snd second
#define all(c) ((c).begin()), ((c).end())
#define TEST(s) if (!(s)) { cout << __LINE__ << " " << #s << endl; exit(-1); }

using Real = double;
struct DualNumber {
  Real a, b; // a + epsilon b
  DualNumber(Real a = 0, Real b = 0) : a(a), b(b) { }
  DualNumber &operator+=(DualNumber x) { b+=x.b; a+=x.a; return *this; }
  DualNumber &operator*=(DualNumber x) { b=b*x.a+a*x.b; a*=x.a; return *this; }
  DualNumber operator+() const { return *this; }
  DualNumber operator-() const { return {-a, -b}; }
  DualNumber inv() const { return {1.0/a, -b/(a*a)}; }
  DualNumber &operator-=(DualNumber x) { return *this += -x; }
  DualNumber &operator/=(DualNumber x) { return *this *= x.inv(); }
};
DualNumber operator+(DualNumber x, DualNumber y) { return x += y; }
DualNumber operator-(DualNumber x, DualNumber y) { return x -= y; }
DualNumber operator*(DualNumber x, DualNumber y) { return x *= y; }
DualNumber operator/(DualNumber x, DualNumber y) { return x /= y; }

// define functions with its derivative
DualNumber pow(DualNumber x, Real e) { return {pow(x.a,e),x.b*pow(x.a,e-1)}; }
DualNumber sqrt(DualNumber x) { return pow(x,0.5); }
DualNumber exp(DualNumber x) { return {exp(x.a),x.b*exp(x.a)}; }
DualNumber cos(DualNumber x) { return {cos(x.a),-x.b*sin(x.a)}; }
DualNumber sin(DualNumber x) { return {sin(x.a),x.b*cos(x.a)}; }
DualNumber tan(DualNumber x) { return sin(x)/cos(x); }
DualNumber log(DualNumber x) { return {log(x.a),x.b/x.a}; }

int main() {
  DualNumber x = 3, y = 4;
  auto f = [&](DualNumber x) {
    return sin(x*x) + cos(exp(x)) + tan(x);
  };
  x.b = 1; // set infinitesimal part
  cout << f(x).b << endl;
}


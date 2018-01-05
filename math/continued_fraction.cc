//
// Continued Fraction 
//
// Description:
//   see https://perl.plover.com/yak/cftalk/
//
// Gosper's algorithm
//
#include <bits/stdc++.h>
using namespace std;

#define fst first
#define snd second
#define all(c) ((c).begin()), ((c).end())

using Int = long long;
struct Generator {
  using Pointer = shared_ptr<Generator>;
  virtual Int value() = 0;
  virtual Pointer next() = 0;
};

// (a x + b)/(c x + d)
struct HoloGen : Generator {
  Int a, b, c, d;
  Pointer x;
  HoloGen(Pointer x, Int a, Int b, Int c, Int d) : x(x), a(a), b(b), c(c), d(d) { }
  void normalize() {
    while (1) {
      if (c == 0 && d == 0) break;
      double b_ac = (double)a/c, b_bd = (double)b/d;
      if ((Int)b_ac == (Int)b_bd) break;
      if (!x) {
        tie(a,b,c,d) = make_tuple(a,a,c,c);
      } else {
        Int p = x->value(); x = x->next();
        tie(a,b,c,d) = make_tuple(a*p+b,a,c*p+d,c);
      }
    }
  }
  virtual Int value() { normalize(); return a/c; }
  virtual Pointer next() {
    Int r = value();
    if (a-c*r == 0 && b-d*r == 0) return 0;
    return Pointer(new HoloGen(x,c,d,a-c*r,b-d*r));
  }
};
// (a x y + b x + c y + d)/(e x y + f x + g y + h)
struct ArithGen : Generator {
  Int a, b, c, d, e, f, g, h;
  Pointer x, y; 
  ArithGen(Pointer x, Pointer y, Int a, Int b, Int c, Int d, Int e, Int f, Int g, Int h) : x(x), y(y), a(a), b(b), c(c), d(d), e(e), f(f), g(g), h(h) { }
  void normalize() {
    while (1) {
      if (e == 0 && f == 0 && g == 0 && h == 0) break;
      double b_ae = (double)a/e, b_cg = (double)c/g,
             b_bf = (double)b/f, b_dh = (double)d/h;
      if ((Int)b_ae == (Int)b_cg && 
          (Int)b_cg == (Int)b_bf && 
          (Int)b_bf == (Int)b_dh) break;
      auto idiff = [&](double a, double b) {
        if (isinf(a) && isinf(b)) return 0.0;
        if (isinf(a) || isinf(b)) return 1.0/0.0;
        return fabs(a - b);
      };
      if (max(idiff(b_ae, b_cg), idiff(b_bf, b_dh))
        > max(idiff(b_ae, b_bf), idiff(b_cg, b_dh))) {
        if (!x) {
          tie(a,b,c,d,e,f,g,h) = make_tuple(a,b,a,b,e,f,e,f);
        } else {
          Int p = x->value(); x = x->next();
          tie(a,b,c,d,e,f,g,h) = make_tuple(a*p+c,b*p+d,a,b,e*p+g,f*p+h,e,f);
        }
      } else { 
        if (!y) {
          tie(a,b,c,d,e,f,g,h) = make_tuple(a,a,c,c,e,e,g,g);
        } else {
          Int p = y->value(); y = y->next();
          tie(a,b,c,d,e,f,g,h) = make_tuple(a*p+b,a,c*p+d,c,e*p+f,e,g*p+h,g);
        }
      }
    }
  }
  virtual Int value() { normalize(); return a/e; }
  virtual Pointer next() {
    Int r = value();
    if (a-e*r == 0 && b-f*r == 0 && c-g*r == 0 && d-h*r == 0) return 0;
    return Pointer(new ArithGen(x,y,e,f,g,h,a-e*r,b-f*r,c-g*r,d-h*r));
  }
};
struct ContinuedFraction {
  using Pointer = shared_ptr<Generator>;
  Pointer head;
  ContinuedFraction(Pointer head) : head(head) { }
  double toReal() {
    Pointer run = head;
    function<double(Pointer, int)> eval = [&](Pointer run, int depth) {
      if (!run) return 1.0/0.0;
      return run->value() + 1.0 / eval(run->next(), depth-1);
    };
    return eval(head, 20);
  }
};
ContinuedFraction operator+(ContinuedFraction x, ContinuedFraction y) {
  return ContinuedFraction(ContinuedFraction::Pointer(
    new ArithGen(x.head, y.head, 0, 1, 1, 0, 0, 0, 0, 1)));
}
ContinuedFraction operator-(ContinuedFraction x, ContinuedFraction y) {
  return ContinuedFraction(ContinuedFraction::Pointer(
    new ArithGen(x.head, y.head, 0, 1, -1, 0, 0, 0, 0, 1)));
}
ContinuedFraction operator*(ContinuedFraction x, ContinuedFraction y) {
  return ContinuedFraction(ContinuedFraction::Pointer(
    new ArithGen(x.head, y.head, 1, 0, 0, 0, 0, 0, 0, 1)));
}
ContinuedFraction operator/(ContinuedFraction x, ContinuedFraction y) {
  return ContinuedFraction(ContinuedFraction::Pointer(
    new ArithGen(x.head, y.head, 0, 1, 0, 0, 0, 0, 1, 0)));
}
ContinuedFraction operator*(ContinuedFraction x, Int a) {
  return ContinuedFraction(ContinuedFraction::Pointer(
    new HoloGen(x.head, a, 0, 0, 1)));
}

int compare(ContinuedFraction x, ContinuedFraction y) {
  auto i = x.head, j = y.head;
  while (1) {
    if (!i && !j) return 0;
    if (!j || i->value() < j->value()) return -1;
    if (!i || i->value() > j->value()) return +1;
    tie(i, j) = make_tuple(j->next(), i->next());
  }
}
bool operator==(ContinuedFraction x, ContinuedFraction y) { return compare(x, y) == 0; }
bool operator!=(ContinuedFraction x, ContinuedFraction y) { return compare(x, y) != 0; }
bool operator<=(ContinuedFraction x, ContinuedFraction y) { return compare(x, y) <= 0; }
bool operator>=(ContinuedFraction x, ContinuedFraction y) { return compare(x, y) >= 0; }
bool operator<(ContinuedFraction x, ContinuedFraction y) { return compare(x, y) < 0; }
bool operator>(ContinuedFraction x, ContinuedFraction y) { return compare(x, y) > 0; }

struct Rational : ContinuedFraction {
  struct RationalGen : Generator {
    Int num, den;
    RationalGen(Int num, Int den) : num(num), den(den) { }
    virtual Int value() { 
      return num / den - (num < 0); 
    }
    virtual Pointer next() { 
      if (num - value() * den == 0) return 0;
      return Pointer(new RationalGen(den, num - value() * den));
    }
  };
  Rational(Int a, Int b) : 
    ContinuedFraction(Pointer(new RationalGen(a, b))) { }
};
struct NapierGenerator : Generator {
  int i;
  NapierGenerator(int i = 0) : i(i) { }
  virtual Int value() {
    if (i == 0) return 2;
    if (i % 3 == 2) return 2 * (i/3) + 2;
    return 1;
  }
  virtual Pointer next() {
    return Pointer(new NapierGenerator(value()));
  }
};
ContinuedFraction e(ContinuedFraction::Pointer(new NapierGenerator(0)));



int main() {
  ContinuedFraction x = Rational(-8, 11), y = Rational(-1, 2);
  cout << compare(x, y) << endl;
  cout << x.toReal() << endl;
}

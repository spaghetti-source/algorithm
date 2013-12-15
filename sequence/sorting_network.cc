// Optimal Sorting (Sorting Network)
//   Twice faster than std::sort
// 
// Reference:
//   http://www.angelfire.com/blog/ronz/Articles/sn2-13_16_horz_2.gif

#define SW(a,b) if (a > b) swap(a, b);
#define SORT3(a,b,c) SW(b,c)SW(a,b)SW(b,c)
#define SORT4(a,b,c,d) SW(a,b)SW(c,d)SW(b,d)SW(a,c)SW(b,c)
#define SORT5(a,b,c,d,e) \
  SW(b,c)SW(d,e)SW(b,d)SW(a,c)SW(c,e)SW(a,d)SW(a,b)SW(c,d)SW(b,c)
#define SORT6(a,b,c,d,e,f) SW(a,b)SW(a,b)SW(c,d)SW(e,f)SW(a,c)SW(d,f)\
  SW(b,e)SW(a,b)SW(c,d)SW(e,f)SW(b,c)SW(d,e)SW(c,d)

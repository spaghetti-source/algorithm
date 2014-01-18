#include <ctime>
#include <sys/time.h>
double tick() {
  static double old;

  struct timeval tv;
  gettimeofday(&tv, NULL);
  double now = tv.tv_sec + tv.tv_usec * 1e-6, ret = now;

  ret -= old;
  old = now;
  return ret;
}

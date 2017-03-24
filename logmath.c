#include <math.h>
#include <stdio.h>
#include "logmath.h"

double _logadd(const double p, const double q) {
  return p + log1p(exp(q - p));
}

double logadd(const double p, const double q) {
  return (p > q) ? _logadd(p, q) : _logadd(q, p);
}

double logCheckAdd(const double p, const double q) {
  if (p == log(0.0)){
    return q;
  }
  else{
    return logadd(p,q);
  }
}
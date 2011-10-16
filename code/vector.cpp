#include "vector.h"
#include <cmath>

double v3_norm(const listdouble& v3) {
  double norm = sqrt(pow(v3[0], 2.0) + pow(v3[1], 2.0) + pow(v3[2], 2.0));
  return norm;
}

listdouble v3_cross(const listdouble& v1, const listdouble& v2) {
  listdouble v3(3);

  v3[0] = v1[1] * v2[2] - v1[2] * v2[1];
  v3[1] = v1[2] * v2[0] - v1[0] * v2[2];
  v3[2] = v1[0] * v2[1] - v1[1] * v2[0];

  return v3;
}

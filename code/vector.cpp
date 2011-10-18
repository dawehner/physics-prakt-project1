#include "vector.h"
#include <cmath>

double v3_norm(const listdouble& v3) {
  double norm = sqrt(pow(v3[0], 2.0) + pow(v3[1], 2.0) + pow(v3[2], 2.0));
  return norm;
}

listdouble v3_slice(const listdouble& vx, const int start) {
  listdouble v3(3);
  for (int i = 0; i < 3; i++) {
    v3[i] = vx[start + i];
  }

  return v3;
}

listdouble v3_cross(const listdouble& v1, const listdouble& v2) {
  listdouble v3(3);

  v3[0] = v1[1] * v2[2] - v1[2] * v2[1];
  v3[1] = v1[2] * v2[0] - v1[0] * v2[2];
  v3[2] = v1[0] * v2[1] - v1[1] * v2[0];

  return v3;
}

double v3_scalar(const listdouble& v1, const listdouble& v2) {
  return v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2];
}


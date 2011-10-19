#include "main.h"
#include "vector.h"
#include <cmath>

void calc_accel_multiple(listdouble& a, const listdouble& y, const listdouble& m) {
  const int count_bodies = m.size();

  listdouble r_diff(3);
  double r_diff_norm = 0.0;

  for (int i = 0; i < count_bodies; i++) {
    a[i*3 + 2] = a[i*3 + 1] = a[i*3 + 0] = 0.0;

    for (int j = 0; j < count_bodies; j++) {
      if (i == j) {
        continue;
      }

      // Calc the connection vector between the bodies.
      // r_ij
      for (int k = 0; k < 3; k++) {
        r_diff[k] = y[j*3 + k] - y[i*3 + k];
      }
      r_diff_norm = v3_norm(r_diff);

      double r_diff_norm_3 = pow(r_diff_norm, 3.0);
      for (int k = 0; k < 3; k++) {
        a[i*3 + k] = a[i*3 + k] + m[j] * r_diff[k] / r_diff_norm_3;
      }
    }
  }
}

void calc_accel_change_multiple(listdouble& da, const listdouble& y, const listdouble& m) {
  const int count_bodies = m.size();
  const int size = count_bodies * 3;

  listdouble r_diff(3);
  double r_diff_norm = 0.0;
  listdouble v_diff(3);
  double v_diff_norm = 0.0;

  for (int i = 0; i < count_bodies; i++) {
    da[i*3 + 2] = da[i*3 + 1] = da[i*3 + 0] = 0.0;

    for (int j = 0; j < count_bodies; j++) {
      if (i == j) {
        continue;
      }

      // Calc the connection and v_connection.
      // r_ij, v_ij
      for (int k = 0; k < 3; k++) {
        r_diff[k] = y[j*3 + k] - y[i*1 + k];
        v_diff[k] = y[size + j*3 + k] - y[size + i*1 + k];
      }
      r_diff_norm = v3_norm(r_diff);
      v_diff_norm = v3_norm(v_diff);
      double r_diff_norm_3 = pow(r_diff_norm, 3.0);
      double r_diff_norm_5 = pow(r_diff_norm, 5.0);
      double r_diff_v_diff = v3_scalar(r_diff, v_diff);

      for (int k = 0; k < 3; k++) {
        da[i*3 + k] = da[i*3 + k] + m[j] *
        ( (v_diff[k] / r_diff_norm_3)
          - 3 * r_diff_v_diff * r_diff[k] / r_diff_norm_5);
      }
    }
  }
}

void calc_dydx(listdouble& dydx, const listdouble& y, const listdouble& m) {
  const int size = y.size();
  const int size_2 = size / 2;

  // Calc speed.
  for (int i = 0; i < size_2; i++) {
    dydx[i] = y[size_2 + i];
  }

  // Calc accelleration.
  listdouble a(size_2);
  calc_accel_multiple(a, y, m);
  for (int i = 0; i < size_2; i++) {
    // start on the v's.
    dydx[size_2 + i] = a[i];
  }
}

void calc_dydx2(listdouble& dydx2, const listdouble& y, const listdouble& m) {
  const int size = y.size();
  const int size_2 = size/2;

  listdouble a(size);
  listdouble da(size);

  calc_accel_multiple(a, y, m);
  calc_accel_change_multiple(da, y, m);

  for (int i = 0; i < size_2; i++) {
    dydx2[i] = a[i];
  }
  for (int i = size_2; i < size; i++) {
    dydx2[i] = da[i - size];
  }

  
}


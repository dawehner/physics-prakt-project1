#include "integration.h"
#include "vector.h"
#include <cmath>

void integration_euler(listdouble& y, const listdouble& dydx, const listdouble& m, const double& dt) {
  int size = y.size();
  for (int i = 0; i < size; i++) {
    y[i] = y[i] + dt * dydx[i];
  }
}

void integration_heun(listdouble& y, const listdouble& dydx, const listdouble& m, const double& dt) {}
void integration_verlet(listdouble& y, const listdouble& dydx, const listdouble& m, const double& dt) {}
void integration_leapfrog(listdouble& y, const listdouble& dydx, const listdouble& m, const double& dt) {}
void integration_rk4(listdouble& y, const listdouble& dydx, const listdouble& m, const double& dt) {
  int size = y.size();
  listdouble y1 = y;
  listdouble y2(size);
  listdouble y3(size);
  listdouble y4(size);

  listdouble dydx2(size);
  listdouble dydx3(size);
  listdouble dydx4(size);

  // K2
  for (int i = 0; i < size; i++) {
    y2[i] = y[i] + 0.5 * dt * dydx[i];
  }
  calc_dydx(dydx2, y2, m);

  // K3
  for (int i = 0; i < size; i++) {
    y3[i] = y[i] + 0.5 * dt * dydx2[i];
  }
  calc_dydx(dydx3, y3, m);

  // K4
  for (int i = 0; i < size; i++) {
    y4[i] = y[i] + 0.5 * dt * dydx3[i];
  }
  calc_dydx(dydx4, y4, m);

  //y_{k+1} = y_k = h/6 (K1 + 2K2 + 2K3 + k4)
  for (int i = 0; i < size; i++) {
    y[i] = y[i] + (dt/6.0) * (dydx[i] + 2.0 * dydx2[i] + 2.0 * dydx3[i] + dydx4[i]);
  }
}

void calc_accel_multiple(listdouble& a, const listdouble& y, const listdouble& m) {
  const int count_bodies = m.size();
  const int size = count_bodies * 3;

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
  const int size = m.size() * 3;

  // Calc speed.
  for (int i = 0; i < size; i++) {
    dydx[i] = y[size + i];
  }

  // Calc accelleration.
  listdouble a(size);
  calc_accel_multiple(a, y, m);
  for (int i = 0; i < size; i++) {
    // start on the v's.
    dydx[size + i] = a[i];
  }
}


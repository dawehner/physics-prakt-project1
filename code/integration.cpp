#include "integration.h"
#include "vector.h"

void integration_euler(listdouble& y, const listdouble& dydx, const listdouble& m, const double& dt) {
  int size = y.size();
  for (int i = 0; i < size; i++) {
    y[i] = y[i] + dt * dydx[i];
  }
}

void integration_heun(listdouble& y, const listdouble& dydx, const listdouble& m, const double& dt) {}
void integration_verlet(listdouble& y, const listdouble& dydx, const listdouble& m, const double& dt) {}
void integration_leapfrog(listdouble& y, const listdouble& dydx, const listdouble& m, const double& dt) {}
void integration_rk4(listdouble& y, const listdouble& dydx, const listdouble& m, const double& dt) {}

void calc_accel_multiple(listdouble& a, const listdouble& y, const listdouble& m) {
  const int count_bodies = m.size();
  const int size = count_bodies * 3;

  listdouble r_diff(3);
  double r_diff_norm = 0.0;
  double r_diff_helper = 0.0;

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

      r_diff_helper = m[j] / pow(r_diff_norm, 3.0);
      for (int k = 0; k < 3; k++) {
        a[i*3 + k] = a[i*3 + k] + m[j] * r_diff[k] / r_diff_helper;
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

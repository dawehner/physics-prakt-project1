#include "integration.h"
#include "vector.h"
#include <cmath>
#include "dydx.h"

void integration_euler(listdouble& y, const listdouble& dydx, const listdouble& m, double& dt) {
  const int size = y.size();
  for (int i = 0; i < size; i++) {
    y[i] = y[i] + dt * dydx[i];
  }
}


void integration_euler_cromer(listdouble& y, const listdouble& dydx, const listdouble& m, double& dt) {
  const int size = y.size();
  const int size_2 = size/2;
  listdouble dydx_copy = dydx;

  // First calculate the new v's
  for (int i = 0; i < size_2; i++) {
    y[size_2 + i] = dydx_copy[i] = y[size_2 + i] + dt * dydx[size_2 + i];
  }

  // Then the r's
  for (int i = 0; i < size_2; i++) {
    y[i] = y[i] + dt * dydx[i];
  }
}


void integration_heun(listdouble& y, const listdouble& dydx, const listdouble& m, double& dt) {
  const int size = y.size();

  listdouble y2(size);
  listdouble dydx2(size);

  for (int i = 0; i < size; i++) {
    y2[i] = y[i] + dt * dydx[i];
  }
  calc_dydx(dydx2, y2, m);

  for (int i = 0; i < size; i++) {
    y[i] = y[i] + 0.5 * dt * (dydx[i] + dydx2[i]);
  }
}

/**
 * @todo: this calcs the a two times.
 */
void integration_verlet(listdouble& y, const listdouble& dydx, const listdouble& m, double& dt) {
  const int size = y.size();
  const int size_2 = size/2;

  listdouble v_s(size_2);
  // Calc v^~
  for (int i = 0; i < size_2; i++) {
    v_s[i] = y[size_2 + i] + (1.0/2.0) * dt * dydx[size_2 + i];
  }

  // Calc r
  for (int i = 0; i < size_2; i++) {
    y[i] = y[i] + dt * v_s[i];
  }

  listdouble a_n(size_2);
  calc_accel_multiple(a_n, y, m);
  // Calc v
  for (int i = 0; i < size_2; i++) {
    y[size_2 + i] = v_s[i] + (1.0/2.0) * dt * a_n[i];
  }
}

void integration_leapfrog(listdouble& y, const listdouble& dydx, const listdouble& m, double& dt) { }

void integration_rk4(listdouble& y, const listdouble& dydx, const listdouble& m, double& dt) {
  const int size = y.size();
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
    y4[i] = y[i] + dt * dydx3[i];
  }
  calc_dydx(dydx4, y4, m);

  //y_{k+1} = y_k = h/6 (K1 + 2K2 + 2K3 + k4)
  for (int i = 0; i < size; i++) {
    y[i] = y[i] + (dt/6.0) * (dydx[i] + 2.0 * dydx2[i] + 2.0 * dydx3[i] + dydx4[i]);
  }
}

void integration_hermit(listdouble& y, const listdouble& dydx, const listdouble& m, double& dt) {
  const int size = y.size();
  const int size_2 = size/2;
  listdouble da(size_2);
  listdouble y_p(size);

  calc_accel_change_multiple(da, y, m);

  // Now calculate the predicted r,v.
  for (int i = 0; i < size; i++) {
    // r
    if (i < size_2) {
      // r = r + v + a + a'.
      y_p[i] = y[i] + dydx[i] * dt + (1.0/2.0) * dydx[size_2 + i] * pow(dt, 2.0) + (1.0/6.0) * da[i] * pow(dt, 3.0);
    }
    // v
    else {
      y_p[i] = y[i] + dydx[i] * dt + (1.0/2.0) * da[i - size_2] * pow(dt, 2.0);
    }
  }

  // Calculate the predicted a, da.
  listdouble ap(size_2);
  listdouble dap(size_2);

  calc_accel_multiple(ap, y_p, m);
  calc_accel_change_multiple(dap, y_p, m);

  // Calculate dda, ddda.
  listdouble dda(size_2);
  listdouble ddda(size_2);
  for (int i = 0; i < size_2; i++) {
    dda[i] = (- 3.0 * (dydx[size_2 + i] - ap[i]) / (pow(dt, 2.0)) - (2.0 * da[i] + dap[i]) / dt ) * 2.0;
    ddda[i] = (2.0 * ((dydx[size_2 + i] - ap[i]) / pow(dt, 3.0)) + (da[i] + dap[i]) / (pow(dt, 2.0)) ) * 6.0;
  }

  // Calculate the corrected values.
  for (int i = 0; i < size; i++) {
    // r_c
    if (i < size_2) {
      y[i] = y_p[i] + (1.0/24.0) * dda[i] * pow(dt, 4.0) + (1.0/120.0) * ddda[i] * pow(dt, 5.0);
    }
    // v_c
    else {
      y[i] = y_p[i] + (1.0/6.0) * dda[i - size_2] * pow(dt, 3.0) + (1.0/24.0) * ddda[i - size_2] * pow(dt, 4.0);
    }
  }
}

void integration_hermit_iter(listdouble& y, const listdouble& dydx, const listdouble& m, double& dt) {
  const int size = y.size();
  const int size_2 = size/2;
  listdouble da(size_2);
  const int hermit_iteration_count = 5;
  listdouble y_p(size);

  // Copy y
  listdouble y_orig = y;

  calc_accel_change_multiple(da, y, m);

  // Now calculate the predicted r,v.
  for (int i = 0; i < size; i++) {
    // r
    if (i < size_2) {
      // r = r + v + a + a'.
      y_p[i] = y[i] + dydx[i] * dt + (1.0/2.0) * dydx[size_2 + i] * pow(dt, 2.0) + (1.0/6.0) * da[i] * pow(dt, 3.0);
    }
    // v
    else {
      y_p[i] = y[i] + dydx[i] * dt + (1.0/2.0) * da[i - size_2] * pow(dt, 2.0);
    }
  }

  for (int j = 0; j < hermit_iteration_count; j++) {
    // Calculate the predicted a, da.
    listdouble ap(size_2);
    listdouble dap(size_2);

    calc_accel_multiple(ap, y, m);
    calc_accel_change_multiple(dap, y, m);

    // v
    for (int i = 0; i < size_2; i++) {
      y[size_2 + i] = y_orig[size_2 + i] + (1.0/2.0) * (ap[i] + dydx[size_2 + i]) * dt + (1.0/12.0) * (dap[i] - da[i]) * pow(dt, 2.0);
    }

    // r
    for (int i = 0; i < size_2; i++) {
      y[i] = y_orig[i] + (1.0/2.0) * (y[size_2 + i] + y_orig[size_2 + i]) * dt + (1.0/12.0) * (ap[i] - dydx[size_2 + i]) * pow(dt, 2.0);
    }
  }
}

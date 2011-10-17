#include "integration.h"
#include "vector.h"
#include <cmath>
#include "dydx.h"

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

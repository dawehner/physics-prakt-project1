#ifndef INTEGRATION_H
#define INTEGRATION_H

#include "main.h"
#include "integration.cpp"

const int INTEGRATION_METHOD_EULER = 0;
const int INTEGRATION_METHOD_HEUN = 1;
const int INTEGRATION_METHOD_VERLET = 2;
const int INTEGRATION_METHOD_LEAPFROG = 3;
const int INTEGRATION_METHOD_RK4 = 4;
const int INTEGRATION_METHOD_HERMIT = 5;
const int INTEGRATION_METHOD_HERMIT_ITER = 6;
const int INTEGRATION_METHOD_EULER_CROMER = 7;

void integration_euler(listdouble& y, const listdouble& dydx, const listdouble& m, double& dt);
void integration_euler_cromer(listdouble& y, const listdouble& dydx, const listdouble& m, double& dt);
void integration_heun(listdouble& y, const listdouble& dydx, const listdouble& m, double& dt);
void integration_verlet(listdouble& y, const listdouble& dydx, const listdouble& m, double& dt);
void integration_leapfrog(listdouble& y, const listdouble& dydx, const listdouble& m, double& dt);
void integration_rk4(listdouble& y, const listdouble& dydx, const listdouble& m, double& dt);
void integration_hermit(listdouble& y, const listdouble& dydx, const listdouble& m, double& dt);
void integration_hermit_iter(listdouble& y, const listdouble& dydx, const listdouble& m, double& dt);

#endif
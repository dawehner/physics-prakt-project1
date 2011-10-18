#ifndef DYDX_H
#define DYDX_H

#include "math.h"
#include "dydx.cpp"

void calc_accel_multiple(listdouble& a, const listdouble& y, const listdouble& m);
void calc_accel_change_multiple(listdouble& da, const listdouble& y, const listdouble& m);
void calc_dydx(listdouble& dydx, const listdouble& y, const listdouble& m);
void calc_dydx2(listdouble& dydx2, const listdouble& y, const listdouble& m);

#endif

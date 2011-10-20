#ifndef QUANTITIES_H
#define QUANTITIES_H

#include "main.h"
#include "quantities.cpp"
#include <list>

double calc_energy(const listdouble& y, const listdouble& m, const int i);
double calc_total_momentum(const listdouble& y, const listdouble& m);
double calc_total_angular_momentum(const listdouble& y, const listdouble& m);
double calc_total_energy(const listdouble& y, const listdouble& m);
listdouble calc_r_center_mass(const listdouble& y, const listdouble& m);
double calc_total_mass(const listdouble& m);
void calc_2body_values_inital(const listdouble& y, const listdouble& m, listdouble& list_total_mass, listdouble& list_start_energy, listdouble& list_start_great_axis, listdouble& list_start_excentric);
void _calc_2body_values(const listdouble& y, const listdouble& m, const int i, listdouble& j_spec, listdouble& e_runge, double& excentric, double& great_half_axis, double& energy);

#endif
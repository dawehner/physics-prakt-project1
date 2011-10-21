#ifndef MAIN_H
#define MAIN_H

#include <vector>
#include <string>

using namespace std;

const int G = 1;

/**
 * Shortcut to get a list of doubles.
 */
typedef vector< long double> listdouble;

void nbody_load_from_file(string& filename, listdouble& y, listdouble& dydx, listdouble& m, double& eta, double& t_max);

void nbody_adapt_timestamp(double& dt, const double& eta, const listdouble& dydx, const listdouble& da);

/**
 * Initialize the problem.
 * - Normalize the problem to M(total_mass) = 1
 * - move the problem to the center of mass.
 */
void nbody_init_problem(listdouble& y, listdouble& m);

/**
 * Add the current positions to the position file.
 */
void nbody_write_pos(ofstream &file_pos, listdouble& y, listdouble& dydx, listdouble& m);

/**
 * Add the current conserved quantities to the conserved file.
 */
void nbody_write_conservered(ofstream& conserved_pos, const double t, const double energy, const double start_energy, const double total_momentum, const double start_total_momentum, const double total_angular_momentum, const double start_total_angular_momentum, const listdouble& r_cm);

void nbody_2body_values_write(ofstream& conserved_2body_file, const double t, const listdouble& y, const listdouble& m, const listdouble& list_total_mass, const listdouble& list_start_energy, const listdouble& list_start_great_axis, const listdouble& list_start_excentric);

#endif
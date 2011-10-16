#ifndef MAIN_H
#define MAIN_H

#include <vector>
#include <string>

using namespace std;

const int G = 1;

/**
 * Shortcut to get a list of doubles.
 */
typedef vector< double> listdouble;

void nbody_load_from_file(string& filename, listdouble& y, listdouble& dydx, listdouble& m, double& eta, double& t_max);

void nbody_adapt_timestamp(double& dt, const double& eta);

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

#endif
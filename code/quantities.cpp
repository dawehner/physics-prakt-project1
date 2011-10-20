#include "quantities.h"
#include "vector.h"
#include <cmath>
#include <bits/algorithmfwd.h>
#include "main.h"
#include "quantities.h"

double calc_energy(const listdouble& y, const listdouble& m, const int j) {
  double energy_kinetic = 0.0;
  const int count_bodies = m.size();
  const int size = count_bodies * 3;
  listdouble v3 = v3_slice(y, size + j * 3);
  const listdouble r3 = v3_slice(y, j*3);

  energy_kinetic += 0.5 * m[j] * pow(v3_norm(v3), 2.0);

  // v_i = sum_j G mj * mi / r_ij
  listdouble r_diff(3);
  double r_diff_norm = 0.0;
  double energy_potential = 0.0;
  for (int i = 0; i < count_bodies; i++) {
    if (i == j) {
      continue;
    }
    for (int k = 0; k < 3; k++) {
      r_diff[k] = r3[k] - y[i*3 + k];
    }
    r_diff_norm = v3_norm(r_diff);
    energy_potential += - m[i] * m[j] / r_diff_norm;
  }

  return energy_potential + energy_kinetic;
}

double calc_total_energy(const listdouble& y, const listdouble& m) {
  const int count_bodies = m.size();

  double total_energy = 0.0;
  for (int i = 0; i < count_bodies; i++) {
    total_energy += calc_energy(y, m, i);
  }

  return total_energy;
}


double calc_total_momentum(const listdouble& y, const listdouble& m) {
  const int count_bodies = m.size();
  const int size_2 = y.size() / 2.0;
  double total_momentum = 0.0;

  listdouble v3(3);
  for (int i = 0; i < count_bodies; i++) {
    v3 = v3_slice(y, size_2 + i * 3);
    total_momentum += m[i] * v3_norm(v3);
  }

  return total_momentum;
}

double calc_total_angular_momentum(const listdouble& y, const listdouble& m) {
  const int count_bodies = m.size();
  const int size_2 = y.size() / 2.0;
  double total_angular_momentum = 0.0;

  listdouble r3(3);
  listdouble p3(3);
  for (int i = 0; i < count_bodies; i++) {
    r3 = v3_slice(y, i * 3);
    p3 = v3_slice(y, i * 3 + size_2);
    for (int k = 0; k < 3; k++) {
      p3[k] = m[i] * p3[k];
    }
    total_angular_momentum += v3_norm(v3_cross(r3, p3));
  }

  return total_angular_momentum;
}

listdouble calc_r_center_mass(const listdouble& y, const listdouble& m) {
  const int count_bodies = m.size();
  listdouble r_cm(3);
  // Iterate over x,y,z
  for (int i = 0; i < 3; i++) {
    // then over each body.
    for (int j = 0; j < count_bodies; j++) {
      // And sum up the values of a certain component of each body.
      r_cm[i] = r_cm[i] + m[j] * y[j * 3 + i];
    }
  }

  return r_cm;
}

double calc_total_mass(const listdouble& m) {
  const int size = m.size();
  double total_mass = 0.0;
  for (int i = 0; i < size; i++) {
    total_mass += m[i];
  }
  return total_mass;
}

void _calc_2body_values(const listdouble& y, const listdouble& m, const int i,
  listdouble& j_spec, listdouble& e_runge, double& excentric, double& great_half_axis, double& energy) {
  const int size_2 = y.size()/2;
  const double total_mass = m[0] + m[i];
  listdouble r_rel(3);
  listdouble v_rel(3);
  listdouble v_cross_j(3);
  double r_rel_norm = 0.0;
  for (int k = 0; k < 3; k++) {
    r_rel[k] = y[k] - y[i * 3 + k];
    v_rel[k] = y[size_2 + k] - y[size_2 + i * 3 + k];
  }
  j_spec = v3_cross(r_rel, v_rel);
  v_cross_j = v3_cross(v_rel, j_spec);
  r_rel_norm = v3_norm(r_rel);
  for (int k = 0; k < 3; k++) {
    e_runge[k] = (v_cross_j[k] / total_mass) - r_rel[k] / r_rel_norm;
  }
  excentric = v3_norm(e_runge);
  great_half_axis = (pow(v3_norm(j_spec), 2.0) / total_mass) / (1.0 - pow(excentric, 2.0));
  energy = calc_energy(y, m, i);
}

void calc_2body_values_inital(const listdouble& y, const listdouble& m, listdouble& list_total_mass, listdouble& list_start_energy, listdouble& list_start_great_axis, listdouble& list_start_excentric) {
  const int count_bodies = m.size();
  for (int i = 1; i < count_bodies; i++) {
    listdouble r_rel(3);
    listdouble v_rel(3);
    listdouble j_spec(3);
    listdouble e_runge(3);
    listdouble v_cross_j(3);
    double excentric = 0.0;
    double great_half_axis = 0.0;
    double energy = 0.0;
    _calc_2body_values(y, m, i, j_spec, e_runge, excentric, great_half_axis, energy);
    list_start_energy.push_back(energy);
    list_start_great_axis.push_back(great_half_axis);
    list_start_excentric.push_back(great_half_axis);
  }
}


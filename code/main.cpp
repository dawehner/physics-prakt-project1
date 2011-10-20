#include <iostream>
#include <fstream>

#include "main.h"
#include "integration.h"
#include "dydx.h"
#include "quantities.h"

#include <string>
#include <vector>

// Boost awesomeness
#include <boost/program_options.hpp>
#include <boost/tokenizer.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/algorithm/string.hpp>

using namespace boost;
namespace po = boost::program_options;

ofstream emptystream;

int main(int argc, char **argv) {
  // Setup variables needed to run the programm

  // Programm options:
  bool adapt_timestamp = true;
  int integration_method = INTEGRATION_METHOD_EULER;
  string output_filename_prefix = "prefix";
  bool write_to_files = true;
  bool break_closed_encounter = false;
  string input_filename = "";
  bool calc_2body_values = false;


  // Stores variables for the time.
  double t_max = 0.0;
  double dt = 0.0;
  double eta = 0.0;
  double t = 0.0;

  // Actual integration variables.
  // Stores all positions and velocities
  listdouble y;
  listdouble dydx;
  listdouble m;

  // Stores a pointer to the current used integration method
  void (*integration_function) (listdouble& y, const listdouble& dydx, const listdouble& m, double& dt) = integration_euler;

  // Load the options to the function.
// Describe the available programm options.
  po::options_description desc("Allowed options", 1024);
  desc.add_options()
    ("help", "Produce help message")
    ("timestamp-adaption,t", po::value<bool>(&adapt_timestamp)->default_value(false), "Adapt timestamp")
    ("integration-method,i", po::value<int>(&integration_method)->default_value(0), "Integration method")
    ("output,o", po::value<string>(&output_filename_prefix)->default_value("output"), "The output file prefix")
    ("write-to-files,w", po::value<bool>(&write_to_files)->default_value(true), "Write to files at all")
    // Currently unsupported as t_max is in the input file.
    //     ("period-counts,c", po::value<int>(&P_count)->default_value(10), "How many orbits should be calculated")
    //     ("steps-per-orbit,s", po::value<int>(&steps_per_orbit)->default_value(100), "The initial amount of steps per orbit")
    ("input,f", po::value<string>(&input_filename)->default_value(""), "Specify the file which has the initial parameters")
    ("break-closed-encounters,e", po::value<bool>(&break_closed_encounter)->default_value(false), "Should the programm be stoped if a closed encounter is detected")
    ("calc-2body-values", po::value<bool>(&calc_2body_values)->default_value(false), "Should the two body values(specific impuls, great half axis, excentric) be calculated.")
    ;
  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, desc), vm);
  po::notify(vm);

  if (vm.count("help")) {
    cout << desc << endl;
    return EXIT_FAILURE;
  }

  switch (integration_method) {
    case INTEGRATION_METHOD_EULER:
      integration_function = integration_euler;
      break;
    case INTEGRATION_METHOD_EULER_CROMER:
      integration_function = integration_euler_cromer;
      break;
    case INTEGRATION_METHOD_HEUN:
      integration_function = integration_heun;
      break;
    case INTEGRATION_METHOD_VERLET:
      integration_function = integration_verlet;
      break;
    case INTEGRATION_METHOD_RK4:
      integration_function = integration_rk4;
      break;
    case INTEGRATION_METHOD_HERMIT:
      integration_function = integration_hermit;
      break;
    case INTEGRATION_METHOD_HERMIT_ITER:
      integration_function = integration_hermit_iter;
      break;
//     case INTEGRATION_ANALYTIC:
//       integration_method = integration_analytic;
//       break;
//     case INTEGRATION_METHOD_LEAPFROG:
//       integration_function = integration_leap_frog;
//       break;
  }

  ofstream pos_file;
  string pos_file_name = output_filename_prefix + ".dat";
  if (write_to_files) {
    pos_file.open(pos_file_name.c_str());
  }
  ofstream conserved_file;
  string conserved_file_name = output_filename_prefix + "-conserved.dat";
  if (write_to_files) {
    conserved_file.open(conserved_file_name.c_str());
  }
  ofstream conserved_2body_file;
  string conserved_2body_file_name = output_filename_prefix + "-conserved-2body.dat";
  if (write_to_files) {
    conserved_2body_file.open(conserved_2body_file_name.c_str());
  }

  // Load the initial parameters
  nbody_load_from_file(input_filename, y, dydx, m, eta, t_max);
  dt = eta;
  const int count_bodies = m.size();
  listdouble da(count_bodies * 3);

  // Setup the initial values
  nbody_init_problem(y, m);

  // Calculate the initial conserved quantities.
  // Additional calculate things like the specific momentum, if needed.
  listdouble list_total_mass(1), list_start_energy(1), list_start_great_axis(1), list_start_excentric(1);
  double start_total_energy, start_total_momentum, start_total_angular_momentum = 0.0;
  if (calc_2body_values) {
    calc_2body_values_inital(y, m, list_total_mass, list_start_energy, list_start_great_axis, list_start_excentric);
  }
  start_total_energy = calc_total_energy(y, m);
  start_total_momentum = calc_total_momentum(y, m);
  start_total_angular_momentum = calc_total_angular_momentum(y, m);

  // Store initial conserved quantities

  // The main loop
  while (t < t_max) {
    // @TODO Calculate the next positions/velos etc.
    calc_dydx(dydx, y, m);
    integration_function(y, dydx, m, dt);
    // @TODO Adapt the timestamp and update the current time.
    calc_accel_change_multiple(da, y, m);

    if (adapt_timestamp) {
      nbody_adapt_timestamp(dt, eta, dydx, da);
    }

    if (write_to_files) {
      nbody_write_pos(pos_file, y, dydx, m);
    }
    if (calc_2body_values) {
      nbody_2body_values_write(conserved_2body_file, t, y, m, list_start_energy, list_start_energy, list_start_great_axis, list_start_excentric);
    }

    double total_energy, total_momentum, total_angular_momentum = 0.0;
    listdouble r_cm(3);
    total_momentum = calc_total_momentum(y, m);
    total_angular_momentum = calc_total_angular_momentum(y, m);
    r_cm = calc_r_center_mass(y, m);
    total_energy = calc_total_energy(y, m);
    nbody_write_conservered(conserved_file, t, abs(total_energy), start_total_energy,
                            total_momentum, start_total_momentum,
                            total_angular_momentum, start_total_angular_momentum, r_cm);

    t += dt;
  }

  // Stop the programm.
  if (write_to_files) {
    pos_file.close();
    conserved_file.close();
  }

  return 0;
}

void nbody_load_from_file(string& filename, listdouble& y, listdouble& dydx, listdouble& m, double& eta, double& t_max) {
  ifstream file(filename.c_str());
  string input_str;

  // On which line of the file are we at the moment.
  int line_count = 0;

  // Variables we store.
  // How many bodies do we integrate.
  int count_bodies = 0;

  typedef boost::tokenizer<boost::char_separator<char> > tokenizer;
  char_separator<char> sep(" ");
  while (getline(file, input_str)) {
    trim(input_str);
    // drop_empty_tokens removes for example "".
    tokenizer tokens(input_str, sep);

    // Calculate the initial steps.
    if (line_count == 0) {
      int count = 0;
      for (tokenizer::iterator beg = tokens.begin(); beg != tokens.end(); ++beg) {
        switch (count++) {
          case 0:
            count_bodies = lexical_cast<int>(*beg);
            // Resize to the amount of bodies we have.
            y.resize(count_bodies * 3 * 2);
            dydx.resize(count_bodies * 3 * 2);
            m.resize(count_bodies);

            break;
          case 1:
            t_max = lexical_cast<double>(*beg);
            break;
          case 2:
            eta = lexical_cast<double>(*beg);
            break;
        }
      }
    }

    // Here comes the masses
    // For example 2 bodies:
    // first line: 0
    // m1: 1
    // m2: 2
    else if (line_count <= count_bodies) {
      m[line_count - 1] = lexical_cast<double>(*tokens.begin());
    }

    // Here comes the initial positions and velocities.
    // In the file the r and v are stored in seperate lines.
    else {
      int line_rv = (line_count - (1 + count_bodies));
//       int current_body = line_rv / 2;

      int count = 0;
      for (tokenizer::iterator beg = tokens.begin(); beg != tokens.end(); ++beg) {
        double value = lexical_cast<double>(*beg);
        y[line_rv * 3 + count] = value;
        count++;
      }
    }
    line_count++;
  }

//   cout << "Loaded initial data:" << "N: " << count_bodies << " eta: " << eta << " t_max: " << t_max << endl;
}


// void nbody_adapt_timestamp(double& dt, const double& eta) { }

void nbody_init_problem(listdouble& y, listdouble& m) {
  // Normalize the masses.
  double total_mass = 0.0;

  int size = m.size();
  for (int i = 0; i < size; i++) {
    total_mass += m[i];
  }
  for (int i = 0; i < size; i++) {
    m[i] = m[i] / total_mass;
  }
  total_mass = 1;

  // Find the center of mass.
  // \vec{r}_cm = \frac{1}{M} \cdot \sum_i^N m_i \cdot \vec{r}_i
  listdouble r_cm(3);

  r_cm = calc_r_center_mass(y, m);

  // Now we have the center of mass.

  // Let's move all components according to the center of mass.
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < size; j++) {
      y[j*3 + i] = y[j*3 + i] - r_cm[i];
    }
  }

  // Now do the same for the velocities as well.
  listdouble v_cm(3);
  int vel_offset = 3 * size;

  // Iterate over x,y,z
  for (int i = 0; i < 3; i++) {
    // then over each body.
    for (int j = 0; j < size; j++) {
      // And sum up the values of a certain component of each body.
      v_cm[i] = v_cm[i] + m[j] * y[vel_offset + j * 3 + i];
    }
    v_cm[i] = v_cm[i] / total_mass;
  }

  // Now we have the center of mass.

  // Let's move all components according to the center of mass.
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < size; j++) {
      y[vel_offset + j*3 + i] = y[vel_offset + j*3 + i] - v_cm[i];
    }
  }
}

void nbody_adapt_timestamp(double& dt, const double& eta, const listdouble& dydx, const listdouble& da) {
  double min = numeric_limits<double>::max();
  double min_help = 0.0;
  const int size_2 = dydx.size() / 2;
  const int count_bodies = dydx.size() / 6;

  listdouble a3(3);
  listdouble da3(3);
  double a3_norm = 0.0;
  double da3_norm = 0.0;
  for (int i = 0; i < count_bodies; i++) {
    for (int k = 0; k < 3; k++) {
      a3[k] = dydx[size_2 + i*3 + k];
      da3[k] = da3[i*3 + k];
    }
    da3_norm = v3_norm(da3);
    a3_norm = v3_norm(a3);
    min_help = a3_norm / da3_norm;
    min = min_help < min ? min_help : min;
  }
  dt = eta * min;
}

void nbody_write_pos(ofstream &pos_file, listdouble& y, listdouble& dydx, listdouble& m) {
  if (pos_file.is_open()) {
    int count_bodies = m.size();
    for (int i = 0; i < count_bodies; i++) {
      // output r, then v
      for (int j = 0; j < 6; j++) {
        pos_file << scientific << y[i * 6 + j] << " ";
      }
      // output a
      for (int j = 3; j < 6; j++) {
        pos_file << scientific << dydx[j] << " ";
      }
      pos_file << m[i];
      pos_file << endl;
    }
  }
}

void nbody_write_conservered(ofstream &conserved_pos, const double t,
  const double energy, const double start_energy,
  const double total_momentum, const double start_total_momentum,
  const double total_angular_momentum , const double start_total_angular_momentum,
  const listdouble& r_cm) {
  if (conserved_pos.is_open()) {
    conserved_pos
     << t << " "
     << abs(start_energy - energy) << " "
     << abs(start_total_momentum - total_momentum) << " "
     << abs(start_total_angular_momentum - total_angular_momentum) << " "
     << r_cm[0] << " "
     << r_cm[1] << " "
     << r_cm[2] << " "
     << endl;
  }
}


void nbody_2body_values_write(ofstream &conserved_2body_file, const double t, const listdouble& y, const listdouble& m, const listdouble& list_total_mass, const listdouble& list_start_energy, const listdouble& list_start_great_axis, const listdouble& list_start_excentric) {
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
    if (conserved_2body_file.is_open()) {
      conserved_2body_file << t << " "
      << abs(energy - list_start_energy[i])<< " "
      << v3_norm(j_spec) << " "
      << abs(excentric - list_start_excentric[i]) << " "
      << great_half_axis << " "
      << endl;
    }
  }
}

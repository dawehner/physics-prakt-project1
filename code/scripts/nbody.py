# @file
# General helper functions for the nbody code.

from math import *


# Create a init file by calculating r, v with a,e and mass.
# @param filename
#   The filename into which the initial data should be written.
# @param values
#   An array of arrays which have
#    - 0: mass
#    - 1: Excentric
#    - 2: Great half axis
#    - TODO: 3: variable to move the body to another place, for example \pi
def nbody_generate_init_file_ea(filename, values, t_max = 2, eta = 0.1):
  output = ""
  m1 = 1.0
  counter = 0
  count_bodies = len(values)

  # The first line with the count of bodies, eta and t_max.
  output += "{0} {1} {2}\n".format(count_bodies, t_max, eta);

  # Now N lines with the masses
  line = ""
  for value in values:
    line += "{0}\n".format(value[0])
  output += line

  # Now the r,v.
  for value in values:
    if counter == 0:
      output += "{0} {1} {2} {3} {4} {5} \n".format(0.0, 0.0, 0.0, 0.0, 0.0, 0.0)

    else:
      m = float(value[0])
      e = float(value[1])
      a = float(value[2])
      v = float(value[3])
      x = a * (1.0 + e) * v
      y = 0.0
      vx = 0.0
      vy = sqrt(((1.0 * (m1 + m)) / a) * (1.0 - e) / (1.0 + e)) * v
      # set vz = z = 0
      output += "{0} {1} {2}\n{3} {4} {5} \n".format(x, y, 0.0, vx, vy, 0.0);

    counter = counter + 1

  input_file = open(filename, 'w')
  input_file.write(output)
  input_file.close()

# @file
# General helper functions for the nbody code.

from math import *
#import subprocess
import os;
import Gnuplot;


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
  m1 = 1.0
  counter = 0
  count_bodies = len(values)
  input_data = []

  # The first line with the count of bodies, eta and t_max.
  input_data.append([count_bodies, t_max, eta])

  # Now N lines with the masses
  for value in values:
    input_data.append([value[0]])

  # Now the r,v.
  for value in values:
    if counter == 0:
      input_data.append([0.0, 0.0, 0.0, 0.0, 0.0, 0.0])

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
      input_data.append([x, y, 0.0, vx, vy, 0.0]);

    counter = counter + 1

  nbody_generate_init_file(filename, input_data)

def nbody_generate_init_file(filename, lines):
  output = ""
  # Convert the number of prticles to a int.
  lines[0][0] = int(lines[0][0])
  for line in lines:
    output += " ".join(map(str, line))
    output += "\n"

  input_file = open(filename, 'w')
  input_file.write(output);
  input_file.close()


def nbody_load_init_file(filename):
  input_file = open(filename, "r")
  values = []

  for line in input_file:
    value = []
    line = line.strip()
    line_values = line.rsplit(" ")
    line_values = filter(lambda x: bool(x), line_values)
    for line_value in line_values:
      value.append(float(line_value))
    values.append(value)
  return values

# Run the nbody programm with a certain input params
# and move the data to a folder.
def nbody_run_collect(params, folder):
  default_params = {
    "--integration-method": 4,
    "--write-to-files": 1,
    "--input": "'in2\.txt'",
  }
  default_params.update(params)

  input = ["./build/nbody"]
  for key, value in default_params.iteritems():
    input.append(key)
    input.append(str(value))

  program = " ".join(input)
  os.system(program)

  if default_params['--write-to-files']:
    os.system("rm %s -Rf" %(folder))
    os.mkdir(folder)
    os.system("mv output*.dat {0}/".format(folder))
  #subprocess.call(input, shell=True)

def nbody_output_gnuplot_abstract(name = "", rows = [], custom_settings = {}):
  # Provide some default settings
  settings = {
    "terminal": "pngcairo font 'Droid Sans,14' linewidth 1 rounded size 1280,1024",
    "size": "square",
    "output": "'{0}.png'".format(name.replace('.', ''))
  }

  settings.update(custom_settings)

  filename = "output{}.dat".format(name)
  directory = name + "-result"

  plot = Gnuplot.Gnuplot()

  # Apply the settings
  for key,value in settings.iteritems():
    if key in ['size', 'terminal', 'log']:
      plot("set {} {}".format(key, value))
    else:
      plot("set {} '{}'".format(key, value))


  plots = []

  for value in rows:
    args = {'using':'', 'title':'', 'every':'','style':''}
    if value.has_key('using'):
      args['using'] = "using {}".format(value['using'])
    if value.has_key('title'):
      args['title'] = "title '{}'".format(value['title'])
    if value.has_key('every'):
      args['every'] = "every {}".format(value['every'])
    if value.has_key('style'):
      args['style'] = "with {}".format(value['style'])
    if value.has_key('filename'):
      args['filename'] = "{}".format(value['filename'])
    else:
      args['filename'] = "{}/{}".format(directory, filename)

    plots.append("'{filename}' {using} {every} {title} {style}".format(**args))

  plot("plot {}".format(", ".join(plots)))

  plot("quit")

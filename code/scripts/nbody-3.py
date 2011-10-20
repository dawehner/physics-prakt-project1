from nbody import *
import time

input_files = {100: 'pl.100.txt', 1000: 'pl.1k.txt'}


rows = []
for count_bodies, data_filename in input_files.iteritems():
  input_data = nbody_load_init_file(data_filename)
  input_data[0][2] = 0.005
  data_filename = data_filename + "-adapted"
  nbody_generate_init_file(data_filename, input_data)
  t_0 = time.time()
  nbody_run_collect({"--input": data_filename, "--integration-method": 2, "--write-to-files" : 1}, "output/3/nbody-{0}".format(count_bodies));
  print "{0} bodies needed {1} seconds".format(count_bodies, time.time() - t_0)

  rows.append({
    "title" : "count_bodies: {0}".format(count_bodies),
    "using" : "1:2",
    "filename": "output/3/nbody-{0}/output-conserved.dat".format(count_bodies)
  })

# Plot energies
settings = {"title": "Energy: {0} bodies".format(count_bodies), "output": "output/3/energy.png", "log" : "y"}
nbody_output_gnuplot_abstract("energy-3", rows, settings)
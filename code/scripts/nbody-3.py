from nbody import *
import time

input_files = {100: 'pl.100.txt', 1000: 'pl.1k.txt'}

for count_bodies, data_filename in input_files.iteritems():
  t_0 = time.time()
  nbody_run_collect({"--input": data_filename, "--integration-method": 4, "--write-to-files" : 0},
      "output/3/nbody");
  print "{0} bodies needed {1} seconds".format(count_bodies, time.time() - t_0)

from nbody import *

eta_list = [0.5, 0.1, 0.05, 0.001, 0.0005, 0.0001]
data = orig_data = nbody_load_init_file("./in2.txt")
used_integrators = {6: "Hermit_Iter."}
#used_integrators = {0: "Euler", 1: "Heun", 2: "Verlet", 4: "RK4", 5: "Hermit", 6: "Hermit_Iter."}

# 1) First part: different eta's.
for eta in eta_list:
  data[0][2] = eta
  # Generate a input file, run the programm and plot E, e, a foreach integration method.
  data_filename = "nbody-2-eta-{0}_input".format(str(eta).replace('.', ''))
  nbody_generate_init_file(data_filename, data)

  for method, name in used_integrators.iteritems():
    print "run eta: {0} method: {1}".format(eta, name)
    nbody_run_collect({"--input": data_filename, "--integration-method": method, "--calc-2body-values" : 1},
    "output/2/nbody-2-eta-{0}-method-{1}".format(eta, method))

print "Start plotting"
# Generate one plot for the energy/e, a forech integration_method but with different eta.
for method, name in used_integrators.iteritems():

  rows = []
  settings = {"title": "Energy: {0}".format(name), "output": "output/2/energy-{0}.png".format(name), "log" : "y"}
  for eta in eta_list:
    rows.append({"title" : "eta = {0}".format(eta), "using" : "1:2",
    "filename": "output/2/nbody-2-eta-{0}-method-{1}/output-conserved-2body.dat".format(eta, method)})

  nbody_output_gnuplot_abstract(name, rows, settings)

  rows = []
  settings = {"title": "Specific Momentum: {0}".format(name), "output": "output/2/momentum-{0}.png".format(name), "log" : "y"}
  for eta in eta_list:
    rows.append({"title" : "eta = {0}".format(eta), "using" : "1:3",
    "filename": "output/2/nbody-2-eta-{0}-method-{1}/output-conserved-2body.dat".format(eta, method)})

  nbody_output_gnuplot_abstract(name, rows, settings)

  rows = []
  settings = {"title": "Excentric: {0}".format(name), "output": "output/2/excentric-{0}.png".format(name), "log" : "y"}
  for eta in eta_list:
    rows.append({"title" : "eta = {0}".format(eta), "using" : "1:4",
    "filename": "output/2/nbody-2-eta-{0}-method-{1}/output-conserved-2body.dat".format(eta, method)})

  nbody_output_gnuplot_abstract(name, rows, settings)

  rows = []
  settings = {"title": "Great Half axis: {0}".format(name), "output": "output/2/axis-{0}.png".format(name), "log" : "y"}
  for eta in eta_list:
    rows.append({"title" : "eta = {0}".format(eta), "using" : "1:5",
    "filename": "output/2/nbody-2-eta-{0}-method-{1}/output-conserved-2body.dat".format(eta, method)})

  nbody_output_gnuplot_abstract(name, rows, settings)

# 2) Second part: Different excentric start values foreach integration method.

from nbody import *

#eta_list = [0.5, 0.1, 0.05, 0.01, 0.005, 0.001, 0.0005]
#data = orig_data = nbody_load_init_file("./in2.txt")
#used_integrators = {6: "Hermit_Iter."}
#used_integrators = {0: "Euler", 1: "Heun", 2: "Verlet", 4: "RK4", 5: "Hermit", 8: "Hermit_Iter1", 9: "Hermit_Iter2", 10: "Hermit_Iter5"}

## 1) First part: different eta's.
#for eta in eta_list:
  #data[0][2] = eta
  ## Generate a input file, run the programm and plot E, e, a foreach integration method.
  #data_filename = "nbody-2-eta-{0}_input".format(str(eta).replace('.', ''))

  #for method, name in used_integrators.iteritems():
    #print "run eta: {0} method: {1}".format(eta, name)

    ## Plot more for verlet to show the constance of the energy over long time.
    #if (method == 2):
      #data[0][1] = 2 * pi * 20
    #else:
      #data[0][1] = 10
    #nbody_generate_init_file(data_filename, data)
    #nbody_run_collect({"--input": data_filename, "--integration-method": method, "--calc-2body-values" : 1},
    #"output/2/nbody-2-eta-{0}-method-{1}".format(eta, method))

#print "Start plotting"
## Generate one plot for the energy/e, a forech integration_method but with different eta.
#for method, name in used_integrators.iteritems():

  #rows = []
  #settings = {"title": "Energy: {0}".format(name), "output": "output/2/energy-{0}.png".format(name), "log" : "y"}
  #for eta in eta_list:
    #rows.append({"title" : "eta = {0}".format(eta), "using" : "1:2",
    #"filename": "output/2/nbody-2-eta-{0}-method-{1}/output-conserved-2body.dat".format(eta, method)})

  #nbody_output_gnuplot_abstract(name, rows, settings)

  #rows = []
  #settings = {"title": "Specific Momentum: {0}".format(name), "output": "output/2/momentum-{0}.png".format(name), "log" : "y"}
  #for eta in eta_list:
    #rows.append({"title" : "eta = {0}".format(eta), "using" : "1:3",
    #"filename": "output/2/nbody-2-eta-{0}-method-{1}/output-conserved-2body.dat".format(eta, method)})

  #nbody_output_gnuplot_abstract(name, rows, settings)

  #rows = []
  #settings = {"title": "Excentric: {0}".format(name), "output": "output/2/excentric-{0}.png".format(name), "log" : "y"}
  #for eta in eta_list:
    #rows.append({"title" : "eta = {0}".format(eta), "using" : "1:4",
    #"filename": "output/2/nbody-2-eta-{0}-method-{1}/output-conserved-2body.dat".format(eta, method)})

  #nbody_output_gnuplot_abstract(name, rows, settings)

  #rows = []
  #settings = {"title": "Great Half axis: {0}".format(name), "output": "output/2/axis-{0}.png".format(name), "log" : "y"}
  #for eta in eta_list:
    #rows.append({"title" : "eta = {0}".format(eta), "using" : "1:5",
    #"filename": "output/2/nbody-2-eta-{0}-method-{1}/output-conserved-2body.dat".format(eta, method)})

  #nbody_output_gnuplot_abstract(name, rows, settings)

# Easy plotting of rk4 output for the tex file.
rows = []
settings = {"title": "Position of the two bodies", "output": "output/2/positions.png"}
rows.append({"title": "Body1", "using": "1:2", "every": "2::1", "filename": "output/2/nbody-2-eta-0.001-method-4/output.dat", "style": "dots"})
rows.append({"title": "Body2", "using": "1:2", "every": "2::2", "filename": "output/2/nbody-2-eta-0.001-method-4/output.dat", "style": "dots"})
nbody_output_gnuplot_abstract("", rows, settings)

# 2) Second part: Different excentric start values foreach integration method.

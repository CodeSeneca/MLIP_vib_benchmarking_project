#! /apps/python/3.12-conda/bin/python3

###############################################################################
#####
##### Author: Maximilian Bechtel <maxi.bechtel@fau.de>
#####         Julien Steffen     <julien.steffen@fau.de>
#####
###############################################################################

import sys
from input_output import read_command_line_arguments

if __name__ == "__main__":
  ##############################################################################
  ##### MAIN PROGRAM
  ##############################################################################

  # Process command line arguments
  mode, input_file = read_command_line_arguments()

  if not mode:
    print("")
    print("Usage: dynamic.py [options]\n")
    print("The options are:")
    print("-input=[input_file]:  Run a calculation specified in [input_file]")
    print("-average:             Determine average values from a previous MD run")
    print("                      stored in md.log")
    print("")
    sys.exit(-1)
  elif mode == "average":
    from calc import calc_averages
    calc_averages()
  elif mode == "calculate" and not input_file:
    print("")
    print("No suitable input file given.")
    print("")
    sys.exit(-1)
  elif mode == "calculate":

    from calc import set_pes
    from input_output import read_input_file, read_atoms

    # Read in the input file
    pes_method, mace_mlip_type, mace_mlip_file, uma_model, uma_task, \
    gptff_file, ocp_model, ocp_cache, ensemble, thermostat, T_init, pressure, pfactor, \
    num_steps, dt, num_freq, smass, tchain, seed, dispersion, stationary, \
    zero_rotation, device, pdamp, pchain, \
    npt_method, mattersim_model, cueq = read_input_file(input_file)

    # Print information about the input parameters given
    print("")
    print("Used method to calculate energies and gradients:", pes_method)
    if pes_method == "mace":
      print("Used model type:", mace_mlip_type)
      print("Used MACE model:", mace_mlip_file)
      if cueq == True:
        print("cuEquivariance acceleration will be used.")
    elif pes_method == "uma" and uma_model == "small":
      print("Used uma model: uma-s-1p1")
      print("Performed task:", uma_task)
    elif pes_method == "uma" and uma_model == "medium":
      print("Used uma model: uma-m-1p1")
      print("Performed task:", uma_task)
    elif pes_method == "orb":
      print("Currently only the orb_v3_conservative_inf_omat model is available")
    elif pes_method == "gptff":
      print("Used GPTFF model:", gptff_file)
    elif pes_method == "ocpcalc":
      print("The OCPCalculator is used for a MLIP")
      print("Used model file:", ocp_model)
      print("Local cache: ", ocp_cache)
    elif pes_method == "mattersim":
      print("Used model type:", mattersim_model)
    elif pes_method == "sevennet":
      print("Currently only the SevenNet-MF-ompa-mpa model is available")

    else:
      print("No suitable model given. Aborting with exit code -5 ...")
    # if uma_task not "oc20" or uma_task not "omat" or uma_task not "omol" or uma_task not "odac" or uma_task not "omc":
    #   print("No suitable UMA task given. Aborting with exit code -6 ...")
    if dispersion:
      print("The Grimme D3 dispersion correction will be applied.")
    if device == "cuda":
      print("Cuda acceleration will be used.")
    elif device == "cpu":
      print("The calculation will be performed on the CPU(s).")
    else:
      print("No suitable acceleration method given (cpu or cuda). Aborting with exit code -3 ...")
      sys.exit(-3)
    print("")

    # Read in the initial geometry and set the calculator
    atoms_object = read_atoms("POSCAR")
    set_pes(atoms_object, pes_method, mace_mlip_type, mace_mlip_file, \
    uma_model, uma_task, gptff_file, ocp_model, ocp_cache, dispersion, \
    device, mattersim_model, cueq)

    # Initialize the MD simulation and create a dynamics object for it
    # but only if ensemble as Master Keyword is set
    if ensemble == "nvt" or "npt":

      from md import init_md, run_md

      dynamics_object = init_md(atoms_object, ensemble, thermostat, T_init, \
      pressure, pfactor, seed, stationary, zero_rotation, dt, smass, tchain, \
      pdamp, pchain, npt_method)

      if not dynamics_object:
        print("The MD dynamics object could not be initialized. Aborting with exit code -4 ...")
        sys.exit(-4)

      if ensemble == "nvt":
        print("A MD simulation in the NVT ensemble will be performed.")
        print("Thermostat:             ", thermostat)
        if thermostat == "nose-hoover":
          print("Number of chains:       ", tchain)
        print("Temperature (K):        ", T_init)
        print("Smass (fs):             ", smass)
        print("Tchain:                 ", tchain)
        print("Time step dt:           ", dt)
        print("Total number of steps:  ", num_steps)
        print("Write out frequency:    ", num_freq)
      elif ensemble == "npt" and npt_method == "parrinello_rahman":
        print("A MD simulation in the NpT ensemble will be performed.")
        print("Temperature (K):        ", T_init)
        print("External pressure (bar):", pressure)
        print("Smass (fs):             ", smass)
        print("Pfactor (GPa*fs^2):     ", pfactor)
        print("Time step dt:           ", dt)
        print("Total number of steps:  ", num_steps)
        print("Write out frequency:    ", num_freq)
      elif ensemble == "npt" and npt_method == "isotropic_mtk":
        print("A MD simulation in the NpT ensemble will be performed.")
        print("Integrator:              Isotropic Martyna-Tobias-Klein (MTK)")
        print("Temperature (K):        ", T_init)
        print("External pressure (bar):", pressure)
        print("External pressure (ev/A^3):", pressure*6.241509074e-7)
        print("Smass (fs):             ", smass)
        print("Pdamp (fs):             ", pdamp)
        print("Tchain:                 ", tchain)
        print("Pchain:                 ", pchain)
        print("Time step dt:           ", dt)
        print("Total number of steps:  ", num_steps)
        print("Write out frequency:    ", num_freq)
      elif ensemble == "npt" and npt_method == "mtk":
        print("A MD simulation in the NpT ensemble will be performed.")
        print("Integrator:              Martyna-Tobias-Klein (MTK)")
        print("Temperature (K):        ", T_init)
        print("External pressure (bar):", pressure)
        print("External pressure (ev/A^3):", pressure*6.241509074e-7)
        print("Smass (fs):             ", smass)
        print("Pdamp (fs):             ", pdamp)
        print("Tchain:                 ", tchain)
        print("Pchain:                 ", pchain)
        print("Time step dt:           ", dt)
        print("Total number of steps:  ", num_steps)
        print("Write out frequency:    ", num_freq)

      # Run the MD simulation
      print("")
      print("Entering MD loop ...")
      run_md(atoms_object, dynamics_object, num_steps, ensemble, num_freq)
    else:
      print("")
      print("No suitable calculation type has been given.")
      print("")

  print("")
  print("********** dynamic.py EXITED NORMALLY **********")
  print("")


#!/usr/bin/env python3
###############################################################################
#####
##### Author: Maximilian Bechtel <maxi.bechtel@fau.de>
#####
###############################################################################


import sys
from calc import set_pes, calc_averages
from input_output import read_command_line_arguments, read_input_file, read_atoms

if __name__ == "__main__":
  ##############################################################################
  ##### MAIN PROGRAM
  ##############################################################################

  # Process command line arguments
  mode, input_file = read_command_line_arguments()

  if mode == "average":
    calc_averages()

  if mode == "calculate" and not input_file:
    print("")
    print("No suitable input file given.")
    print("Usage: dynamic.py -input=[input_filename]")
    print("")
    sys.exit(-1)

<<<<<<< HEAD
  if mode == "calculate":
    # Read in the input file
    pes_method, mace_mlip_type, mace_mlip_file, uma_model, uma_task, ensemble, \
    thermostat, T_init, pressure, pfactor, num_steps, dt, num_freq, smass, \
    tchain, seed, dispersion, stationary, zero_rotation, \
    device, pdamp, pchain, npt_method = read_input_file(input_file)

    # Print information about the input parameters given
    print("")
    print("Used method to calculate energies and gradients:", pes_method)
    if pes_method == "mace":
      print("Used model type:", mace_mlip_type)
      print("Used MACE model:", mace_mlip_file)
    elif pes_method == "uma" and uma_model == "small":
      print("Used uma model: uma-s-1p1")
      print("Performed task:", uma_task)
    elif pes_method == "uma" and uma_model == "medium":
      print("Used uma model: uma-m-1p1")
      print("Performed task:", uma_task)
    elif pes_method == "orb":
      print("Currently only the orb_v3_conservative_inf_omat model is available")
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
=======
  # Read in the input file
  pes_method, mace_mlip_type, mace_mlip_file, uma_model, uma_task, gptff_file, ensemble, \
  thermostat, T_init, pressure, pfactor, num_steps, dt, num_freq, smass, \
  tchain, seed, dispersion, stationary, zero_rotation, \
  device, pdamp, pchain, npt_method = read_input_file(input_file)
  # Set number of chains in the Nose-Hoover chain to three
  num_chains=3

  # Print information about the input parameters given
  print("")
  print("Used method to calculate energies and gradients:", pes_method)
  if pes_method == "mace":
    print("Used model type:", mace_mlip_type)
    print("Used MACE model:", mace_mlip_file)
  elif pes_method == "uma" and uma_model == "small":
    print("Used uma model: uma-s-1p1")
    print("Performed task:", uma_task)
  elif pes_method == "uma" and uma_model == "medium":
    print("Used uma model: uma-m-1p1")
    print("Performed task:", uma_task)
  elif pes_method == "orb":
    print("Currently only the orb_v3_conservative_inf_omat model is available")
  elif pes_method == "gptff":
    print("Used GPTFF model:",gptff_file)
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
  uma_model, uma_task, gptff_file, dispersion, device)

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
        print("Number of chains:       ", num_chains)
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
>>>>>>> 5ece41147f7a2cc36c689fab1f19dbf62733849c
    print("")

    # Read in the initial geometry and set the calculator
    atoms_object = read_atoms("POSCAR")
    set_pes(atoms_object, pes_method, mace_mlip_type, mace_mlip_file, \
    uma_model, uma_task, dispersion, device)

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
        print("External pressure (eV/A^3):", pressure*6.241509074e-7)
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
      print("")

      # Run the MD simulation
      print("Entering MD loop ...")
      print("")
      run_md(atoms_object, dynamics_object, num_steps, num_freq)
    else:
      print("")
      print("No suitable calculation type has been given.")
      print("")

  print("")
  print("********** dynamic.py EXITED NORMALLY **********")
  print("")


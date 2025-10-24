#! /apps/python/3.12-conda/bin/python

###############################################################################
#####
##### Author: Maximilian Bechtel <maxi.bechtel@fau.de>
#####
###############################################################################

import sys
from calc import set_pes
from input_output import read_command_line_arguments, read_input_file, read_atoms

if __name__ == "__main__":
  ##############################################################################
  ##### MAIN PROGRAM
  ##############################################################################

  # Process command line arguments
  input_file = read_command_line_arguments()
  if not input_file:
    print("")
    print("No suitable input file given.")
    print("Usage: dynamic.py -input=[input_filename]")
    print("")
    sys.exit(-1)

  # Read in the input file
  pes_method, mace_mlip_type, mace_mlip_file, uma_model, uma_task, ensemble, \
  thermostat, T_init, num_steps, dt, num_freq, smass, num_chains, seed, \
  dispersion, stationary, zero_rotation, device = read_input_file(input_file)

  # Print information about the input parameters given
  print("")
  print("Used method to calculate energies and gradients:", pes_method)
  if pes_method == "mace":
    print("Used model type:", mace_mlip_type)
    print("Used MACE model:", mace_mlip_file)
  if pes_method == "uma" and uma_model == "small":
    print("Used uma model: uma-s-1p1")
    print("Performed task:", uma_task)
  elif pes_method == "uma" and uma_model == "medium":
    print("Used uma model: uma-m-1p1")
    print("Performed task:", uma_task)
  else:
    print("No suitable UMA model given. Aborting with exit code -5 ...")
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
  uma_model, uma_task, dispersion, device)

  # Initialize the MD simulation and create a dynamics object for it
  # but only if ensemble as Master Keyword is set
  if ensemble == "nvt":

    from md import init_md, run_md

    print("")
    print("A MD simulation in the NVT ensemble will be performed.")

    dynamics_object = init_md(atoms_object, ensemble, thermostat, T_init, \
    seed, stationary, zero_rotation, dt, smass, num_chains)
    if pes_method == "mace" and mace_mlip_type == "omol" and dispersion==True:
      print("WARNING: The MACE OMOL models have been selected. Currently MACE OMOL does not support D3 corrections natively!")
      print("D3 corrections will be added manually by the TorchDFTD3Calculator, which is still experiemtally.")

    if not dynamics_object:
      print("The MD dynamics object could not be initialized. Aborting with exit code -4 ...")
      print("")
      sys.exit(-4)

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


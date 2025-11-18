// Copyright 2025
// Author: Maximilian Bechtel, <maxi.bechtel@fau.de>
//
// Program: eval_aimd
// Evaluate the results from a VASP AIMD calculation
//

#include <iostream>
#include "./io.h"

// Main function to drive everything
int main(int argc, char** argv) {
  // Object to handle input/output
  IO io(argc, argv);
  // Exit code of various functions
  int exit_code = 0;

  // Print information about usage of the script
  io.print_information();
  // Process command line arguments
  exit_code = io.process_cla();
  if (exit_code != 0) {
    return exit_code;
  }

  // Read in the POSCAR file
  exit_code = io.read_poscar();
  if (exit_code != 0) {
    return exit_code;
  }
  // Read in the OUTCAR file
  exit_code = io.read_outcar();
  if (exit_code != 0) {
    return exit_code;
  }

  // Write the md.log file
  io.write_md_log();
  // Write the MACE trainset file
  io.write_mace();

  return exit_code;
}


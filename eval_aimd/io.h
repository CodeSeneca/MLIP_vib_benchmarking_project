// Copyright 2025
// Author: Maximilian Bechtel, <maxi.bechtel@fau.de>
//
// Module io: process input and output
//

#ifndef IO_H_
#define IO_H_

#include <string>
#include <vector>
#include "./data.h"

// Class to handle input and output related tasks
class IO {
 private:
  // Number of command line arguments
  int _argc;
  // Array of the actual command line arguments
  char** _argv;
  // Path to the POSCAR file
  std::string _poscar_path;
  // Path to the OUTCAR file
  std::string _outcar_path;
  // Name of the MACE trainset file
  static const char* mace_file;
  // Frequency in which MD steps should be read
  size_t _read_freq;
  // Determine if the MACE trainset file should be written
  bool _write_mace;

  // Total number of atoms
  size_t _num_atoms;
  // Vector of all element symbols
  std::vector<std::string> _element_symbols;
  // Vector of all element numbers
  std::vector<size_t> _element_numbers;

  // Vector of all electronic energies
  std::vector<double> _energies_el;
  // Vector of all unit cells
  std::vector<Cell> _cells;
  // Vector of all molecules
  std::vector<std::vector<Atom>> _molecules;

 public:
  // Constructor to get argc and argv
  IO(int argc, char** argv);

  // Process command line arguments
  int process_cla();

  // Print information about usage of the script
  void print_information();

  // Read in the POSCAR file to get all element numbers
  int read_poscar();

  // Read in the OUTCAR file
  int read_outcar();

  // Write a MACE trainset file
  void write_mace();
};

#endif  // IO_H_


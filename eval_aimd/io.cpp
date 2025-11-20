// Copyright 2025
// Author: Maximilian Bechtel, <maxi.bechtel@fau.de>
//
// Module io: process input and output
//

#include <stdlib.h>
#include <unistd.h>
#include <stdio.h>
#include <getopt.h>

#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include <vector>
#include <iomanip>

#include "./io.h"
#include "./data.h"

const char* IO::mace_file = "trainset_mace.xyz";
const char* IO::md_log_file = "md.log";

// _____________________________________________________________________________
IO::IO(int argc, char** argv) {
  _argc = argc;
  _argv = argv;

  _num_atoms = 0;
  _read_freq = 1;
  _poscar_path = "POSCAR";
  _outcar_path = "OUTCAR";
  _write_mace = false;
}

// _____________________________________________________________________________
int IO::process_cla() {
  // Definition of all options
  struct option options[] = {
    { "write_mace", 0, NULL, 0 },
    { "read_freq", 1, NULL, 'r' },  // optional
    { "outcar_file", 1, NULL, 'o' },  // optional
    { "poscar_file", 1, NULL, 'p' },  // optional
    { "help", 0, NULL, 'h' },  // optional
    { NULL, 0, NULL, 0 }
  };

  // Process all command line arguments
  optind = 1;
  int option_index = 0;
  while (true) {
    char c = getopt_long(_argc, _argv, "r:o:p:h", options, NULL);
    if (c == -1) break;  // Leave loop if all arguments have been processed

    switch (c) {
      case 'r':
        _read_freq = atoi(optarg);
        break;
      case 'o':
        _outcar_path = optarg;
        break;
      case 'p':
        _poscar_path = optarg;
        break;
      case 'h':
        printf("Usage: ./eval_aimd [options]\n"
               "\n"
               "--write_mace                  write a MACE training set file\n"
               "--read_freq, -r <n>:          process only each nth MD step\n"
               "                              default: 1\n"
               "--outcar_file, -o <filename>: path to OUTCAR file\n"
               "                              default: OUTCAR\n"
               "--poscar_file, -p <filename>: path to POSCAR file\n"
               "                              default: POSCAR\n");
        return 1;
      case 0:
        if ( std::string(options[option_index].name) == "write_mace" ) {
          _write_mace = true;
          break;
        }
    }
  }

  return 0;
}

// _____________________________________________________________________________
void IO::print_information() {
  std::cout << "PROGRAM eval_aimd: Evaluate the results from a VASP AIMD" <<
  " calculation" << std::endl << std:: endl;
}

// _____________________________________________________________________________
int IO::read_poscar() {
  // Open the POSCAR file
  std::string line;
  std::stringstream line_stream;
  std::string symbol;
  int number;

  std::cout << "Reading in the POSCAR file ..." << std::endl;
  std::cout << "-------------------------------------" << std::endl;

  std::ifstream poscar_file(_poscar_path);

  if (!poscar_file.is_open()) {
    std::cerr << "Error while opening file " << _poscar_path << std::endl;
    return 1;
  }

  // Skip the first 5 lines
  for (int i = 1; i <= 5; i++) {
    getline(poscar_file, line);
  }
  // Read in all atomic symbols
  getline(poscar_file, line);
  line_stream.str(line);
  // Split the string into vector _element_symbols along whitespaces
  while (line_stream >> symbol) {
    _element_symbols.push_back(symbol);
  }
  std::cout << "Found elements: ";
  for (const std::string& element : _element_symbols) {
    std::cout << element << " ";
  }
  std::cout << std::endl;

  // Read in the number of atoms for each atomic symbol
  getline(poscar_file, line);
  line_stream.clear();
  line_stream.str(line);
  while (line_stream >> number) {
    _element_numbers.push_back(number);
  }
  std::cout << "Found numbers: ";
  for (int atomic_number : _element_numbers) {
    std::cout << atomic_number << " ";
  }
  std::cout << std::endl;

  // Add up all numbers
  for (int atomic_number : _element_numbers) {
    _num_atoms += atomic_number;
  }
  std::cout << "Total number of atoms: " << _num_atoms << std::endl;
  std::cout << "-------------------------------------" << std::endl
  << std::endl;

  // Close the POSCAR file
  poscar_file.close();

  return 0;
}

// _____________________________________________________________________________
int IO::read_outcar() {
  // One line of the OUTCAR file
  std::string line;
  std::stringstream line_stream;
  // Dummy string for reading
  std::string dummy;
  // Unit cell
  Cell cell;
  // One atom in one MD step
  Atom atom;
  // Electronic energy of one MD step
  double energy_el;
  // Potential, kinetic and total energy of one MD step
  double epot, ekin, etot;
  // Vector of atoms representing the molecule of one MD step
  std::vector<Atom> molecule;

  std::cout << "Reading in the OUTCAR file ..." << std::endl;

  // Open the OUTCAR file
  std::ifstream outcar_file(_outcar_path);

  if (!outcar_file.is_open()) {
    std::cerr << "Error while opening file " << _outcar_path << std::endl;
    return 1;
  }
  // Read in the OUTCAR file until EOF reached
  while (getline(outcar_file, line)) {
    if (line.find("VOLUME and BASIS-vectors are now :") != std::string::npos) {
      // Skip the next two lines
      getline(outcar_file, line);
      getline(outcar_file, line);
      // Read in the current volume
      getline(outcar_file, line);
      line_stream.clear();
      line_stream.str(line);
      line_stream >> dummy >> dummy >> dummy >> dummy >> cell.volume;

      // Skip the next line
      getline(outcar_file, line);
      // Vector a1
      getline(outcar_file, line);
      line_stream.clear();
      line_stream.str(line);
      line_stream >> cell.a1[0] >> cell.a1[1] >> cell.a1[2];

      // Vector a2
      getline(outcar_file, line);
      line_stream.clear();
      line_stream.str(line);
      line_stream >> cell.a2[0] >> cell.a2[1] >> cell.a2[2];

      // Vector a3
      getline(outcar_file, line);
      line_stream.clear();
      line_stream.str(line);
      line_stream >> cell.a3[0] >> cell.a3[1] >> cell.a3[2];

      _cells.push_back(cell);

    } else if (line.find("TOTAL-FORCE") != std::string::npos) {
      // Skip the next line
      getline(outcar_file, line);
      // Clear the old molecule
      molecule = std::vector<Atom>();
      // Read in atomic positions and gradients
      while (getline(outcar_file, line)) {
        if (line.find("---") != std::string::npos) {
          break;
        }

        line_stream.clear();
        line_stream.str(line);
        line_stream >> atom.pos[0] >> atom.pos[1] >> atom.pos[2]
        >> atom.grad[0] >> atom.grad[1] >> atom.grad[2];

        molecule.push_back(atom);
      }

      _molecules.push_back(molecule);
    } else if (line.find("FREE ENERGIE OF THE ION-ELECTRON SYSTEM")
               != std::string::npos) {
      // Skip the next three lines
      getline(outcar_file, line);
      getline(outcar_file, line);
      getline(outcar_file, line);
      // Read in the electronic energy (sigma->0)
      getline(outcar_file, line);
      line_stream.clear();
      line_stream.str(line);
      line_stream >> dummy >> dummy >> dummy >> dummy >> dummy >> dummy
      >> energy_el;

      _energies_el.push_back(energy_el);
    } else if (line.find("ENERGY OF THE ELECTRON-ION-THERMOSTAT SYSTEM (eV)")
      != std::string::npos) {
      // Skip the next line
      getline(outcar_file, line);
      // Read in the energy ion-electron TOTEN = Epot
      getline(outcar_file, line);
      line_stream.clear();
      line_stream.str(line);
      line_stream >> dummy >> dummy >> dummy >> dummy >> epot >> dummy >> dummy;
      // std::cout << epot << std::endl;
      _epot.push_back(epot);
      // Read in Ekin
      getline(outcar_file, line);
      line_stream.clear();
      line_stream.str(line);
      line_stream >> dummy >> dummy >> dummy >> dummy >> ekin;
      _ekin.push_back(ekin);
      // Caluclate Etot = Ekin + Epot
      etot = ekin + epot;
      _etot.push_back(etot);
      }
  }

  // Close the OUTCAR file
  outcar_file.close();

  return 0;
}

// _____________________________________________________________________________
void IO::write_md_log() {
  // The current MD step that is written to md.log
  size_t step = 1;

  std::cout << "Writing the md.log file ..." << std::endl;
  if (_read_freq == 1) {
    std::cout << "Writing out every step" << std::endl;
  } else if (_read_freq == 2) {
    std::cout << "Writing out every " << _read_freq << "nd step" << std::endl;
  } else if (_read_freq == 2) {
    std::cout << "Writing out every " << _read_freq << "rd step" << std::endl;
  } else {
    std::cout << "Writing out every " << _read_freq << "nd step" << std::endl;
  }

  // Open file md.log for writing
  std::ofstream md_log(md_log_file);

  // Write the header of the md.log file
  md_log << "# MD step     free energy TOTEN (eV)     kinetic energy (eV)";
  md_log << "     Epot+Ekin (eV)    volume of cell (A^3)"
  << std::endl << std::endl;

  for (size_t i = 0; i < _cells.size(); i++) {
    if (step%_read_freq == 0) {
      md_log << "    " << step << "            " << std::fixed <<
      std::setprecision(5) <<
      _epot[i] << "                  " << _ekin[i] << "           "
      << _etot[i] << "              " <<
      std::fixed << std::setprecision(2) << _cells[i].volume << "     "
      << std::endl;
    }
    step++;
  }

  // Close the md.log file
  md_log.close();
}

// _____________________________________________________________________________
void IO::write_mace() {
  if (_write_mace) {
    // The current MD step that is written to the MACE trainset file
    size_t step = 1;

    std::cout << "-------------------------------------" << std::endl;
    std::cout << "Number of MD steps: " << _energies_el.size() << std::endl;
    std::cout << "Write out frequency: " << _read_freq << std::endl;
    std::cout << _energies_el.size()/_read_freq <<
    " MD steps will be evaluated" << std::endl;
    std::cout << "Writing MD steps from " << _read_freq << " to " <<
    _energies_el.size()/_read_freq * _read_freq << std::endl;
    std::cout << "-------------------------------------" << std::endl;
    std::cout << "Writing the MACE trainset file ..." << std::endl;

    // Open the MACE trainset file for writing
    std::ofstream trainset_mace(mace_file);

    for (size_t i = 0; i < _energies_el.size(); i++) {
      if (step%_read_freq == 0) {
        trainset_mace << "     " << _num_atoms << std::endl;

        trainset_mace << "\"Lattice=    " << std::fixed << std::setprecision(8)
        << _cells[i].a1[0] << "   " << _cells[i].a1[1] << "  "<< _cells[i].a1[2]
        << "  "
        << _cells[i].a2[0] << "   " << _cells[i].a2[1] << "  "<< _cells[i].a2[2]
        << "  "
        << _cells[i].a3[0] << "   " << _cells[i].a3[1] << "  "<< _cells[i].a3[2]
        << " Properties=species:S:1:pos:R:3:molID:I:1:REF_forces:R:3"
        << " Nmols=1 REF_energy=    " << _energies_el[i]<< " pbc=\"T T T"<< "\""
        << std::endl;

        // The current element that is written out
        size_t current_element = 0;
        // How many lines of the current element were already written
        size_t written_for_element = 0;
        // First element symbol
        for (size_t l = 0; l < _molecules[i].size(); l++) {
          // Check if the current element was completely written
          if (written_for_element >= _element_numbers[current_element]) {
            current_element++;
            written_for_element = 0;
          }

          trainset_mace << _element_symbols[current_element] << "    " <<
          std::fixed << std::setprecision(8) <<
          _molecules[i][l].pos[0] << "   " << _molecules[i][l].pos[1] << "   "<<
          _molecules[i][l].pos[2] << "   " << "0        " <<
          _molecules[i][l].grad[0] << "   " << _molecules[i][l].grad[1] << "   "
          << _molecules[i][l].grad[2] << std::endl;

          written_for_element++;
        }

        // The rest of the element symbols
        // for (size_t k = 1; k < _element_symbols.size(); k++) {
        //   for (size_t j = _element_numbers[k-1]; j <
        //        _element_numbers[k-1]+_element_numbers[k]; j++) {
        //     trainset_mace << _element_symbols[k] << "    " <<
        //     std::fixed << std::setprecision(8) <<
        //     _molecules[i][j].pos[0] << "   " << _molecules[i][j].pos[1]
        //     << "   " <<
        //     _molecules[i][j].pos[2] << "   " << "0        " <<
        //     _molecules[i][j].grad[0]<< "   "<< _molecules[i][j].grad[1]<<" "
        //     << _molecules[i][j].grad[2] << std::endl;
        //   }
        // }
      }
      step++;
    }

    // Close the MACE trainset file
    trainset_mace.close();
  }
}

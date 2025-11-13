// Copyright 2025
// Author: Maximilian Bechtel, <maxi.bechtel@fau.de>
//
// Module data: Definition of mathematical objects and molecule related data
//

#ifndef DATA_H_
#define DATA_H_

#include <string>

// Unit cell of a MD simulation
struct Cell {
  // Lattice vectors of the unit cell
  double a1[3];
  double a2[3];
  double a3[3];
  // Volume of the unit cell
  double volume;
};

// Chemical atom
struct Atom {
  // Element symbol
  std::string symbol;
  // Coordinates
  double pos[3];
  // Gradient
  double grad[3];
};

#endif  // DATA_H_

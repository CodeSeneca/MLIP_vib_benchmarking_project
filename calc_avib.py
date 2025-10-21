#! /apps/python/3.12-conda/bin/python

import sys
import numpy as np
from scipy.integrate import simpson
from scipy.constants import N_A
from scipy.constants import k as kB
from scipy.constants import h as h_planck
import matplotlib.pyplot as plt

# Process command line arguments
def process_cla():
  # Number of command line argument
  argc = len(sys.argv)
  # File with the VDOS data
  filename = "power_spectrum_global.csv"
  # Check whether the normalized VDOS spectrum should be shown
  show = False
  # Temperature
  T = 298.15
  # Vibrational degrees of freedom
  fvib = 9

  for arg in sys.argv:
    arg_split = arg.split("=")

    if arg_split[0] == "-file":
      filename = arg_split[1]
    elif arg_split[0] == "-show":
      show = True
    elif arg_split[0] == "-T":
      T = float(arg_split[1])
    elif arg_split[0] == "-fvib":
      fvib = int(arg_split[1])

  return filename, show, T, fvib

# Plot the VDOS
def plot_vdos(x, y):
  fig, ax = plt.subplots(figsize=(10,7))
  ax.set_xlabel(r"wavenumber ($cm^{-1}$)", fontsize=12)
  ax.set_ylabel(r"$D(\nu)$", fontsize=12)

  ax.plot(x, y, label="DSI", linewidth=1.3)
  plt.legend(fontsize=12)

  plt.show()

# Integrate a function given gy sample points (xi, yi)
def integrate(x, y):
  nf = simpson(y, x=x)
  return nf

# Quantum mechanical expression for vibrational frequencies
def W_QM(v, T):
  beta = 1.0/(kB*T)

  w_qm = ( 1.0 - np.exp(-beta*h_planck*v) ) / np.exp(-0.5*beta*h_planck*v)
  w_qm = np.log(w_qm)
  return w_qm

# Calculate the vibrational part of the free energy (in J/mol)
def calc_avib(vibs, vdos_normalized, T):
  # Calculate the QM expression
  dim = vibs.size
  w_qm = np.empty(shape=dim)

  for i in range(0, dim):
    w_qm[i] = W_QM(vibs[i], T)

  # Calculcate the free energy
  prod = vdos_normalized * w_qm
  ln_Qvib = integrate(vibs, prod)

  Avib_J = kB*T * ln_Qvib        # in J
  Avib_J_mol = Avib_J * N_A      # in J/mol

  return Avib_J_mol

if __name__ == "__main__":
  # Process the command line arguments given
  filename, show, T, fvib = process_cla()

  # Read in the VDOS data
  data = np.genfromtxt(filename, delimiter=';', comments='#')
  wavenumbers = data[:,0]
  vdos = data[:,1]

  # Normalize the VDOS spectrum to fvib
  I = integrate(wavenumbers, vdos)
  norm = fvib/I
  vdos_norm = norm * vdos
  I_norm = integrate(wavenumbers, vdos_norm)

  print(f"f_vib before normalization: {I:.5f}")
  print(f"f_vib after normalization: {I_norm:.5f}")

  # Remove the first entry if it is 0.00 cm-1
  if wavenumbers[0] == 0.00:
    wavenumbers = wavenumbers[1:]
    vdos_norm = vdos_norm[1:]

  Avib_J_mol = calc_avib(wavenumbers, vdos_norm, T)    # in J/mol
  Avib_kJ_mol = Avib_J_mol/1000.0                      # in kJ/mol
  Avib_eV = Avib_kJ_mol*0.0104                         # in eV

  print("")
  print(f"Avib = {Avib_kJ_mol:.5f} kJ/mol")
  print(f"Avib = {Avib_eV:.5f} eV")
  print("")

  # Plot the normalized VDOS spectrum
  if show:
    plot_vdos(wavenumbers, vdos_norm)


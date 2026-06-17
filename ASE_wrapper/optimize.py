#! /home/bechtel/progs/envs/mace/bin/python3

import os
import sys
import shutil
from ase.io import read, write
from ase.optimize import BFGS

def process_args():
  """ Process command line arguments """

  # Total number of command line arguments
  num_args = len(sys.argv)
  # File to process: POSCAR or complete XDATCAR with MD trajectory are possible
  filename = None
  # CUDA or CPU
  device = "cpu"
  # Convergence criterion: max_a(|Fa|) < fmax (in ev/A)
  fmax = 0.025
  # Start and end frame to optimize
  start = 0
  end = -1
  # Process only each Nth MD step
  step = 1
  # Extract the N frames with lowest energy
  lowest = None

  # Path to the MACE MLIP model file
  path = None
  # Used head of the multihead MH-1 model
  head = "omol"


  for arg in sys.argv:
    if arg == "-h" or arg == "--help":
      print("Usage: optimize.py [options]\n")
      print("The options are:")
      print("--file=[filename]")
      print("--fmax=[maximum force norm]")
      print("--start=[start frame to optimize]")
      print("--end=[end frame to optimize]")
      print("--step=[process each Nth frame]")
      print("--lowest=[the N frames with lowest energies to extract]")
      print("--device=[cpu or cuda]")
      print("--path=[path to MACE model]")
      print("--head=[head for the multihead MACE MH1 model]")
      print("")
      sys.exit(0)
    else:
      arg_split = arg.split("=")
      if arg_split[0] == "--fmax":
        fmax = float(arg_split[1])
      if arg_split[0] == "--file":
        filename = arg_split[1]
      if arg_split[0] == "--device":
        device = arg_split[1]
      if arg_split[0] == "--start":
        start = int(arg_split[1])
      if arg_split[0] == "--end":
        end = int(arg_split[1])
      if arg_split[0] == "--step":
        step = int(arg_split[1])
      if arg_split[0] == "--lowest":
        lowest = int(arg_split[1])

  return filename, path, head, device, fmax, start, end, step, lowest

def read_molecules(filename, start, end, step):
  """ Read in the molecules from filename """

  atoms = read(filename, index=f"{start}:{end}:{step}")
  return atoms

def set_pes(atoms, path, head, disp=False, dtype="float64", device="cpu", cueq=False):
  """ Set the calculator for the Potential Energy Surface (PES) """

  from mace.calculators import mace_mp

  calc = mace_mp(model=path, head=head, dispersion=disp, default_dtype=dtype, device=device, enable_cueq=cueq)
  atoms.calc = calc

def convert_traj(traj_file):
  """ Convert an ASE .traj file into a VASP XDATCAR file """

  frames = read(traj_file, index=":")
  last_frame = read(traj_file, index=-1)

  # write("XDATCAR", frames, format="vasp-xdatcar")
  write("trajectory.xyz", frames, format="xyz")
  write("CONTCAR", last_frame, format="vasp")

def init_geometry_optimization(atoms):
  """ Perform a geometry optimization for the given atoms object """

  dyn = BFGS(atoms=atoms, logfile="-", trajectory="opt.traj")
  return dyn

def main():
  # List of all minimum energies
  energies = []

  filename, path, head, device, fmax, start, end, step, lowest = process_args()
  if not filename:
    print("No suitable file for geometry optimization given. Possible are a POSCAR or a complete XDATCAR from a previous MD run.")
    sys.exit(-1)

  print(f"\nReading in the molecular structures from {filename} ...")
  print(f"From frame {start} to frame {end} every {step}th step will be processed.")
  molecules = read_molecules(filename=filename, start=start, end=end, step=step)

  print("\n=================================")
  print("      Geometry optimization    ")
  print("=================================")

  print("Algorithm: BFGS")
  print(f"Setting the MACE Calculator to {model} ...")
  if device == "cpu":
    print("The calculations will be performed on the CPU(s).")

  elif device == "cuda":
    print("The calculations will be performed on the GPU.")
  else:
    print("No suitable device given! The options are cpu or cuda.")
    sys.exit(1)

  print(f"Maximum norm of the force vectors fmax: {fmax}")
  print("The trajectory will be written to opt.traj. Open it with ase gui opt.traj ...\n")

  # Perform a geometry optimization for each structure in filename
  with open("energies.dat", 'w') as energy_file:
    energy_file.write("# Frame no.     Energy (eV)     rel. Energy (kJ/mol)\n")

    for i, mol in enumerate(molecules):
      set_pes(mol, path=path, head=head, device=device)
      dyn = init_geometry_optimization(mol)

      os.makedirs(f"frame{i*step + start}", exist_ok=True)
      os.chdir(f"frame{i*step + start}")

      dyn.run(fmax=fmax)
      e_min = mol.get_potential_energy()
      energies.append(e_min)
      print(f"Minimum Energy: {e_min:.5f}\n")
      energy_file.write(f"{(i*step + start):11}     {e_min:15.5f}     {((e_min - energies[0])*96.485):12.2f}\n")
      energy_file.flush()

      # Convert the .traj file written by ASE to a VASP XDATCAR
      print("Converting opt.traj into a xyz file ...")
      print("The last frame will be written to CONTCAR ...\n")
      convert_traj("opt.traj")

      os.chdir("..")

    # Get the 20 frames with the lowest energy
    if lowest:
      lowest_energies = sorted( enumerate(energies), key=lambda x: x[1] )[0:lowest]
      energy_file.write("\n# The frames with lowest energies are:\n")

      for pair in lowest_energies:
        if pair[0]%20 == 0 and pair[0] != 0:
          energy_file.write(f"{pair[0]*step}\n")
        else:
          energy_file.write(f"{pair[0]*step}  ")

      energy_file.write("\n")

      # Copy the frames with lowest energies into the folder lowest_energies
      print("Copying the frames with lowest energies into the folder lowest_energies ...")
      os.makedirs("lowest_energies", exist_ok=True)
      for pair in lowest_energies:
        src = f"frame{pair[0] * step}"
        dst = os.path.join( "lowest_energies", os.path.basename(src) )
        shutil.copytree(src, dst)

  print("\n***** optimize.py exited normally *****")

if __name__ == "__main__":
  main()


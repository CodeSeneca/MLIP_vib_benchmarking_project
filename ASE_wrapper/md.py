import numpy as np

from ase import units
from ase.io import write
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution as Maxwell
from ase.md.velocitydistribution import Stationary, ZeroRotation
from ase.md.nose_hoover_chain import NoseHooverChainNVT as NoseHoover
from ase.md.nose_hoover_chain import IsotropicMTKNPT, MTKNPT
from ase.md.npt import NPT

##############################################################################
##### GLOBAL VARIABLES
##############################################################################

# Average temperature after the whole MD simulation
temp_aver = 0
# Average pressure after the whole MD simulation
pressure_aver = 0
# Average volume after the whole MD simulation
volume_aver = 0
# Average density after the whole MD simulation
density_aver = 0
# Average total energy after the whole MD simulation
etot_aver = 0
# Logfile for energy output in each MD step
logfile_name = "md.log"
# Configuration files for each MD step
traj_filename = "XDATCAR"
conf_filename = "CONTCAR"

##############################################################################
##### FUNCTIONS
##############################################################################

def print_md_step(logfile, traj_file, md_step, atoms, num_freq):
  """Print Epot, Ekin, Etot, temperature T, volume V and density of one MD time step"""

  # Use temp_aver and pressure_aver as global variables here
  global temp_aver
  global pressure_aver
  global volume_aver
  global density_aver
  global etot_aver

  # Flush the RAM buffer to actually write the MD step
  logfile.flush()
  traj_file.flush()

  # Calculate all desired values
  epot = atoms.get_potential_energy()
  ekin = atoms.get_kinetic_energy()
  etot = ekin + epot
  etot_aver += etot

  num_atoms = len(atoms)
  temp = ekin / (1.5 * num_atoms * units.kB)
  temp_aver += temp
  volume = atoms.get_volume()
  volume_aver += volume
  mass_total = sum(atoms.get_masses())
  density = mass_total/volume * 1.66053906660
  density_aver += density

  stress = atoms.get_stress(include_ideal_gas=True)/units.GPa
  # Negative stress means positive pressure -> compression
  pressure_trace = -(stress[0] + stress[1] + stress[2])/3.0
  pressure_aver += pressure_trace

  # Print all calculated values
  md_step = md_step * num_freq
  logfile.write("    %i        %.5f        %.5f       %.5f       %.5f       %.5f        %.5f        %.5f\n" % (md_step, epot, ekin, etot, temp, volume, density, pressure_trace))

  # Print the current configuration to XDATCAR and CONTCAR
  # The CONTCAR file is overwritten in each MD step
  write(traj_file, atoms, format="vasp-xdatcar")
  with open(conf_filename, "w") as conf_file:
    write(conf_file, atoms, format="vasp")

def init_md(atoms, ensemble, thermostat, T_init, pressure, pfactor, seed, \
stationary, zero_rotation, dt, smass, tchain, pdamp, pchain, npt_method):
  """Initialize the MD simulation"""

  # Initialize the velocities randomly via a Maxwell-Boltzmann distribution
  # at the desired initial temperature T_init
  # if seed = None chose the seed randomly
  if not seed:
    print("The seed for the Maxwell Boltzmann distribution will be set randomly.")
    Maxwell(atoms=atoms, temperature_K=T_init)
  else:
    print(f"The seed for the Maxwell Boltzmann distribution was set to {seed}.")
    Maxwell(atoms=atoms, temperature_K=T_init, rng=np.random.RandomState(seed))
  print("")

  # Set the center-of-mass momentum to zero
  if stationary:
    Stationary(atoms=atoms)
    print("The center-of-mass momentum was set to zero.")
  # Set the total angular momentum to zero by counteracting rigid rotations
  if zero_rotation:
    ZeroRotation(atoms=atoms)
    print("The total angular momentum was set to zero.")
  print("")

  # Initialize a NVT dynamics object
  dyn = None

  # Nose-Hoover Chain Thermostat
  if ensemble == "nvt" and thermostat == "nose-hoover":
    dyn = NoseHoover(atoms=atoms,
                     timestep=dt*units.fs,
                     temperature_K=T_init,
                     tdamp=smass*units.fs,
                     tchain=tchain)

  # Combined Nose-Hoover and Parrinello-Rahman dynamics, suggested by Melchionna et al.
  elif ensemble == "npt" and npt_method == "parrinello_rahman":
    dyn = NPT(atoms,
              dt*units.fs,
              temperature_K=T_init,
              externalstress=pressure*units.bar,  # isotropic pressure equals (-p, -p, -p, 0, 0, 0), positive p means compression
              ttime=smass*units.fs,
              pfactor=pfactor*units.GPa*(units.fs**2),
              mask=[[1,0,0],[0,1,0],[0,0,1]])  # Change the volume V but not its shape during the NpT simulation (isotropic changes)

  # Nose-Hoover Chain Thermostat with its pressure extension by Martyna, Tobias and Klein
  elif ensemble == "npt" and npt_method == "isotropic_mtk":
    dyn = IsotropicMTKNPT(atoms=atoms,
                          timestep=dt*units.fs,
                          temperature_K=T_init,
                          pressure_au=pressure*6.241509074e-7,  # in ev/A^3
                          tdamp=smass*units.fs,
                          pdamp=pdamp*units.fs,
                          tchain=tchain,
                          pchain=pchain)

  elif ensemble == "npt" and npt_method == "mtk":
    dyn = MTKNPT(atoms=atoms,
                 timestep=dt*units.fs,
                 temperature_K=T_init,
                 pressure_au=pressure*6.241509074e-7,  # in eV/A^3
                 tdamp=smass*units.fs,
                 pdamp=pdamp*units.fs,
                 tchain=tchain,
                 pchain=pchain)

  return dyn

def run_md(atoms_object, dynamics_object, num_steps, num_freq):
  "Run the MD simulation"

  with open(logfile_name, "w") as logfile, open(traj_filename, "w") as traj_file:

    logfile.write("# Step     pot.energy (eV)    kin.energy (eV)   tot.energy (eV)   temperature (K)  volume(A^3)  density(g/cm^3)  pressure (GaPa)\n")
    logfile.write("\n")


    # Print the initial configuration at t = 0.0 fs
    print_md_step(logfile, traj_file, 0, atoms_object, num_freq)

    for i in range(1, int(num_steps/num_freq) + 1):
      # Move the dynamics one step further
      dynamics_object.run(num_freq)
      print_md_step(logfile, traj_file, i, atoms_object, num_freq)

    # Caluclate the average temperature after the performed MD run
    calc_aver(num_steps, num_freq, logfile)

def calc_aver(num_steps, num_freq, logfile):
  """Calculate the average temperature and pressure from a pervious MD run"""

  # Use temp_aver, pressure_aver as global varibles here
  global temp_aver
  global pressure_aver
  global volume_aver
  global density_aver
  global etot_aver

  temp_aver = temp_aver / (num_steps/num_freq + 1)
  pressure_aver = pressure_aver / (num_steps/num_freq + 1)
  volume_aver = volume_aver / (num_steps/num_freq + 1)
  density_aver = density_aver / (num_steps/num_freq + 1)
  etot_aver = etot_aver / (num_steps/num_freq + 1)
  logfile.write("\n")
  logfile.write(f"# The average temperature of the MD run is:                                {temp_aver:.5f} K\n")
  logfile.write(f"# The average external pressure of the MD run is:                          {pressure_aver:.5f} GaPa\n")
  logfile.write(f"# The average volume of the MD run is:                                     {volume_aver:.5f} A^3\n")
  logfile.write(f"# The average density of the MD run is:                                    {density_aver:.5f} g/cm^3\n")
  logfile.write(f"# The average total energy of the MD run is (i.e. the internal energy U):  {etot_aver:.5f} eV\n")


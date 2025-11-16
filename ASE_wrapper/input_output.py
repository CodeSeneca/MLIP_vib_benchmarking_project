import sys

def read_command_line_arguments():
  """Process the command line arguments given by the user"""

  # Filename of the input file
  mode = None
  input_file = None

  for arg in sys.argv:
    if arg[0:7] == "-input=":
      mode = "calculate"
      input_file = arg[7:]
    if arg[0:8] == "-average":
      mode = "average"

  return mode, input_file

def read_input_file(input_filename):
  """Read in the input file and set all given input parameters"""

  ##############################################################################
  ##### GENERAL INPUT VARIABLES
  ##### HERE DEFAULT VALUES FOR ALL INPUT PARAMETERS ARE GIVEN
  ##############################################################################
  
  # The method that should be used to calculate energies and gradients
  pes_method = "mace"
  # Type of the MACE model
  # mp:   foundation models trained on Materials Project data (for crystals)
  # omol: foundation models trained on OMOL data              (for molecules)
  mace_mlip_type = "mp"
  # Path to the desired MACE model
  mace_mlip_file = None
  # Type of the UMA model
  # uma-s-1p1 (small model) or uma-m-1p1 (medium model)
  uma_model = "small"
  # The task that should be performed with the UMA model
  # Available are:
  #   oc20 :  catalysis
  #   omat :  inorganic materials
  #   omol :  molecules
  #   odac :  metal organic frameworks (MOFs)
  #   omc  :  molecular crystals
  uma_task = "omol"
  # For GPTFF: the file with the model weigths (currently only gptff_v1.pth 
  #   and gptff_v2.pth available)
  gptff_file = None
  # For the OCPCalculator: enables a couple of MLIPs, mostly older ones
  

  device="cpu"

  # Master Keyword to activate the MD routine
  # available are: NVT, NpT
  ensemble = None
  # The desired NVT method
  thermostat = "nose-hoover"
  # The desired NpT method
  npt_method = "isotropic_mtk"
  # The desired initial temperature in K
  T_init = 300.0
  # Desired external pressure in bar
  pressure = 1.0
  # Constant in the Parrinello-Rahman Barostat differential equation (in GPa*fs**2)
  pfactor = 2e6
  # The desired total number of MD steps
  num_steps = 100
  # The desired MD time step in fs
  dt = 1.0
  # The write out frequency (every N steps, energy and structure are written)
  num_freq = 1
  # The fictious mass of the Nose-Hoover Chain Thermostat (as multiples of time steps)
  smass = 40.0
  # The characteristic time scale for the barostat in ASE time units
  pdamp = 1000
  # The chain length of the Nose-Hoover Chain Thermostat
  # NB: MDALGO = 2 IN VASP USES THE CLASSICAL NOSE-HOOVER THERMOSTAT WITHOUT CHAINS
  tchain = 3 # default is 3 in ASE
  # The chain length of the MTK barostat
  pchain = 3
  # Random seed for initialization of the velocities
  seed = None
  # Determine if the Grimme D3 dispersion correction should be used
  dispersion = False
  # Determine if translation of the center-of-mass and total angular momentum
  # should be removed
  stationary = False
  zero_rotation = False
  
  ##############################################################################

  try:
    with open(input_filename, "r") as input_file:
      # Read in the input file line by line
      for line in input_file:
        # Split the line into a Python list
        line_list = line.rstrip().split()
        # Skip the line if it is empty
        if not line_list:
          continue

        if line_list[0] == "ensemble":
          ensemble = line_list[1]
        elif line_list[0] == "steps":
          num_steps = int(line_list[1])
        elif line_list[0] == "dt":
          dt = float(line_list[1])
        elif line_list[0] == "steps_freq":
          num_freq = int(line_list[1])
        elif line_list[0] == "stationary":
          stationary = True
        elif line_list[0] == "zero_rotation":
          zero_rotation = True
        elif line_list[0] == "dispersion":
          dispersion = True
        elif line_list[0] == "seed":
          seed = int(line_list[1])
        elif line_list[0] == "pes":
          pes_method = line_list[1]

        # If the Master Keyword ensemble was set to NVT check if the
        # NVT object was correctly initialized in the input file
        # form: nvt {
        # temperature ...
        # thermostat ...
        # smass ...
        # tchain ...
        # }
        elif ensemble == "nvt" and line_list[0] == "nvt" and line_list[1] == "{":
          next_line = "xxxx"
          while next_line != "}":
            next_line = input_file.readline().rstrip()
            next_line_split = next_line.split()
            if next_line_split[0] == "temperature":
              T_init = float(next_line_split[1])
            if next_line_split[0] == "thermostat":
              thermostat = next_line_split[1]
            if next_line_split[0] == "smass":
              smass = int(next_line_split[1])
            if next_line_split[0] == "tchain":
              tchain = int(next_line_split[1])

        elif pes_method == "mace" and line_list[0] == "mace" and line_list[1] == "{":
          next_line = "xxxx"
          while next_line != "}":
            next_line = input_file.readline().rstrip()
            next_line_split = next_line.split()
            if next_line_split[0] == "mlip_file":
              mace_mlip_file = next_line_split[1]
            if next_line_split[0] == "mlip_type":
              mace_mlip_type = next_line_split[1]
            if next_line_split[0] == "device":
              device = next_line_split[1]

        elif pes_method == "uma" and line_list[0] == "uma" and line_list[1] == "{":
          next_line = "xxxx"
          while next_line != "}":
            next_line = input_file.readline().rstrip()
            next_line_split = next_line.split()
            if next_line_split[0] == "model":
              uma_model = next_line_split[1]
            if next_line_split[0] == "task":
              uma_task = next_line_split[1]
            if next_line_split[0] == "device":
              device = next_line_split[1]

        elif pes_method == "orb" and line_list[0] == "orb" and line_list[1] == "{":
          next_line = "xxxx"
          while next_line != "}":
            next_line = input_file.readline().rstrip()
            next_line_split = next_line.split()
            if next_line_split[0] == "device":
              device = next_line_split[1]

        elif pes_method == "gptff" and line_list[0] == "gptff" and line_list[1] == "{":
          next_line = "xxxx"
          while next_line != "}":
            next_line = input_file.readline().rstrip()
            next_line_split = next_line.split()
            if next_line_split[0] == "mlip_file":
              gptff_file = next_line_split[1]
            if next_line_split[0] == "device":
              device = next_line_split[1]

        elif pes_method == "ocpcalc" and line_list[0] == "ocpcalc" and line_list[1] == "{":
          next_line = "xxxx"
          while next_line != "}":
            next_line = input_file.readline().rstrip()
            next_line_split = next_line.split()
            if next_line_split[0] == "ocp_model":
              ocp_model = next_line_split[1]
            if next_line_split[0] == "local_cache":
              ocp_cache = next_line_split[1]
            if next_line_split[0] == "device":
              device = next_line_split[1]


        elif ensemble == "npt" and line_list[0] == "npt" and line_list[1] == "{":
          next_line = "xxxx"
          while next_line != "}":
            next_line = input_file.readline().rstrip()
            next_line_split = next_line.split()
            if next_line_split[0] == "method":
              npt_method = next_line_split[1]
            if next_line_split[0] == "temperature":
              T_init = float(next_line_split[1])
            if next_line_split[0] == "smass":
              smass = int(next_line_split[1])
            if next_line_split[0] == "pressure":
              pressure = float(next_line_split[1])
            if next_line_split[0] == "pfactor":
              pfactor = float(next_line_split[1])
            if next_line_split[0] == "pdamp":
              pdamp = int(next_line_split[1])
            if next_line_split[0] == "tchain":
              tchain = int(next_line_split[1])
            if next_line_split[0] == "pchain":
              pchain = int(next_line_split[1])

  except FileNotFoundError:
    print("")
    print("Input file", input_filename, "was not found. Aborting with exit code -2 ...")
    print("")
    sys.exit(-2)

  return pes_method, mace_mlip_type, mace_mlip_file, uma_model, uma_task, \
  gptff_file, ocp_model, ocp_cache, ensemble, thermostat, T_init, pressure, pfactor, num_steps, dt, num_freq, \
  smass, tchain, seed, dispersion, stationary, zero_rotation, device, pdamp, \
  pchain, npt_method

def read_atoms(file_type):
  """Read in the geometry from the POSCAR file and return an atoms object
  from it"""

  from ase.io import read

#  if pes_method == "gptff":
#     atoms = adp.get_atoms(file_type)
#  else:
  atoms = read(file_type)

  return atoms


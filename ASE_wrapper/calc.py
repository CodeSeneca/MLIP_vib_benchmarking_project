import numpy as np

def set_pes(atoms, pes_method, mace_mlip_type, mace_mlip_file, uma_model, \
uma_task, gptff_file, ocp_model, ocp_cache, dispersion, device):
  """Set the Calculator used for all calculations"""

  # Define a MACE model with possible D3 dispersion correction
  if pes_method == "mace" and mace_mlip_type == "mp":

    from mace.calculators import mace_mp

    mace_mlip = mace_mp(model=mace_mlip_file, \
    dispersion=dispersion, default_dtype="float32", device=device)
    atoms.calc = mace_mlip
  elif pes_method == "mace" and mace_mlip_type == "omol":

    from mace.calculators import mace_omol
    # from ase.calculators.mixing import SumCalculator
    # from torch_dftd.torch_dftd3_calculator import TorchDFTD3Calculator as DFTD3

    mace_mlip = mace_omol(model=mace_mlip_file, \
    default_dtype="float32", device=device)
    atoms.calc = mace_mlip

    if dispersion:
      print("WARNING: The MACE OMOL models with the Grimme D3 dispersion correction was chosen.")
      print("This is not necessary since the OMOL25 training set is based on the wB97M-VV10 functional,")
      print("which already includes non-local dispersions.")
      dipsersion=False
      print("Dispersion was set to False.")

      # dft_d3_calc = DFTD3(atoms=atoms, device=device, damping="bj")
      # print("")
      # print("Pytorch implementation of DFTD3 initialized.")
      # combined_calc = SumCalculator([mace_mlip, dft_d3_calc])
      # atoms.calc = combined_calc
    # else:
      # atoms.calc = mace_mlip

  # Define an UMA model by the FAIR Chemistry team
  if pes_method == "uma":
    from fairchem.core import pretrained_mlip, FAIRChemCalculator

    # Set the UMA model: uma-s-1p1 (small model) or uma-m-1p1 (medium model)
    if uma_model == "small":
      predictor = pretrained_mlip.get_predict_unit("uma-s-1p1", device=device)
    elif uma_model == "medium":
      predictor = pretrained_mlip.get_predict_unit("uma-m-1p1", device=device)

    # Set the task that should be performed with the UMA model
    # Available are:
    #   oc20 :  catalysis
    #   omat :  inorganic materials
    #   omol :  molecules
    #   odac :  metal organic frameworks (MOFs)
    #   omc  :  molecular crystals
    calc = FAIRChemCalculator(predictor, task_name=uma_task)

    # Add the Grimme D3 correction manually
    if dispersion:
      from ase.calculators.mixing import SumCalculator
      from torch_dftd.torch_dftd3_calculator import TorchDFTD3Calculator as DFTD3

      dft_d3_calc = DFTD3(atoms=atoms, device=device, damping="bj")
      print("")
      print("Pytorch implementation of DFTD3 initialized.")
      combined_calc = SumCalculator([calc, dft_d3_calc])
      atoms.calc = combined_calc
    else:
      atoms.calc = calc

  if pes_method == "orb":
    from orb_models.forcefield import pretrained
    from orb_models.forcefield.calculator import ORBCalculator

    orbff = pretrained.orb_v3_conservative_inf_omat(device=device, precision="float32-high")
    calc = ORBCalculator(orbff, device=device)

    if dispersion:
      from ase.calculators.mixing import SumCalculator
      from torch_dftd.torch_dftd3_calculator import TorchDFTD3Calculator as DFTD3

      dft_d3_calc = DFTD3(atoms=atoms, device=device, damping="bj")
      print("")
      print("Pytorch implementation of DFTD3 initialized.")
      combined_calc = SumCalculator([calc, dft_d3_calc])
      atoms.calc = combined_calc
    else:
      atoms.calc = calc

  if pes_method == "gptff":
    from gptff.model.mpredict import ASECalculator
    from pymatgen.core import Structure
    from pymatgen.io.ase import AseAtomsAdaptor

    model_weight = gptff_file
    p = ASECalculator(model_weight, device)
    adp = AseAtomsAdaptor()
    struc = Structure.from_file('POSCAR')
    atoms = adp.get_atoms(struc)
#    atoms.set_calculator(p)
    atoms.calc=p
    #if dispersion:
    #  from ase.calculators.mixing import SumCalculator
    #  from torch_dftd.torch_dftd3_calculator import TorchDFTD3Calculator as DFTD3
#
#      dft_d3_calc = DFTD3(atoms=atoms, device=device, damping="bj")
#      print("")
#      print("Pytorch implementation of DFTD3 initialized.")
#      combined_calc = SumCalculator([calc, dft_d3_calc])
#      atoms.calc = combined_calc
#    else:
#      atoms.calc = calc

  if pes_method == "ocpcalc":
    from fairchem.core import OCPCalculator
    if device == "cpu":
      cpu=True
    elif device == "cuda":
      cpu=False

    calc = OCPCalculator(
      model_name=ocp_model,
      local_cache=ocp_cache,
      cpu=cpu
    )
    atoms.calc = calc

def calc_averages():
  """Calculate average quantaties from md.log file written by a MD simulation"""

  data = np.genfromtxt("md.log")

  # Etot
  etot = data[:,3]
  etot_mean = np.mean(etot)
  etot_std = np.std(etot)
  # Temperature
  temperature = data[:,4]
  temperature_mean = np.mean(temperature)
  temperature_std = np.std(temperature)
  # Volume
  volume = data[:,5]
  volume_mean = np.mean(volume)
  volume_std = np.std(volume)
  # Density
  density = data[:,6]
  density_mean = np.mean(density)
  density_std = np.std(density)
  # Pressure
  pressure = data[:,7]
  pressure_mean = np.mean(pressure)
  pressure_std = np.std(pressure)

  print("")
  print(f"Performing statistics over {len(etot)} MD time steps")
  print("======================================================================================================")
  print(f"The average temperature of the MD run is:                                {temperature_mean:.5f} K")
  print(f"The average pressure of the MD run is:                                   {pressure_mean*10000:.5f} bar")
  print(f"The average volume of the MD run is:                                     {volume_mean:.5f} A^3")
  print(f"The average density of the MD run is:                                    {density_mean:.5f} g/cm^3")
  print(f"The average total energy of the MD run is (i.e. the internal energy U):  {etot_mean:.5f} eV")
  print("")


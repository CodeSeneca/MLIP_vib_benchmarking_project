def set_pes(atoms, pes_method, mace_mlip_type, mace_mlip_file, dispersion, device):
  """Set the Calculator used for all calculations"""

  # Define a MACE model with possible D3 dispersion correction
  if pes_method == "mace" and mace_mlip_type == "mp":

    from mace.calculators import mace_mp

    mace_mlip = mace_mp(model=mace_mlip_file, \
    dispersion=dispersion, default_dtype="float32", device=device)
    atoms.calc = mace_mlip
  elif pes_method == "mace" and mace_mlip_type == "omol":

    from ase.calculators.mixing import SumCalculator
    from mace.calculators import mace_omol
    from torch_dftd.torch_dftd3_calculator import TorchDFTD3Calculator as DFTD3

    mace_mlip = mace_omol(model=mace_mlip_file, \
    default_dtype="float32", device=device)

    if dispersion:
      dft_d3_calc = DFTD3(atoms=atoms, device=device, damping="bj")
      print("")
      print("Pytorch implementation of DFTD3 initialized.")
      combined_calc = SumCalculator([mace_mlip, dft_d3_calc])
      atoms.calc = combined_calc
    else:
      atoms.calc = mace_mlip


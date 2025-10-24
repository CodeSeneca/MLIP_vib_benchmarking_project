def set_pes(atoms, pes_method, mace_mlip_type, mace_mlip_file, uma_model, \
uma_task, dispersion, device):
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
    atoms.calc = calc


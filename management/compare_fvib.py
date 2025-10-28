#! /apps/python/3.12-conda/bin/python

import sys
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats

################################################################################
##### FUNCTIONS
################################################################################

def process_command_line_arguments():
  """Process the command line arguments given by the user

  the arguments have to be:
  arg1    the first fvib_results.dat file
  arg2    the second fvib_results.dat file
  """

  # The number of arguments given
  num_arguments = len(sys.argv) - 1
  # File 1 and fil2
  arg1 = arg2 = None

  if num_arguments != 2:
    return arg1, arg2

  arg1 = sys.argv[1]
  arg2 = sys.argv[2]

  return arg1, arg2

def get_fvib(filename1, filename2):
  """Get the average Fvib values for each run from both files"""

  fvib_1 = []
  fvib_2 = []

  # Read in Fvib from the two files
  with open(filename1, 'r') as file1, open(filename2, 'r') as file2:
    print(f"Reading average Fvib values from {filename1} ...")
    for line in file1:
      line = line.rstrip().split()
      if not line or line[0] == '#' or len(line) != 1:
        continue
      else:
        value = float(line[0])

      # Check if the Fvib entry is empty
      if value == int(value):
        value = "undefined"
        fvib_1.append(value)
        # Skip the next two lines
        file1.readline()
        file1.readline()
      else:
        fvib_1.append(value)

    print(f"Reading average Fvib values from {filename2} ...")
    for line in file2:
      line = line.rstrip().split()
      if not line or line[0] == '#' or len(line) != 1:
        continue
      else:
        value = float(line[0])

      # Check if the Fvib entry is empty
      if value == int(value):
        value = "undefined"
        fvib_2.append(value)
        # Skip the next two lines
        file2.readline()
        file2.readline()
      else:
        fvib_2.append(value)

  return fvib_1, fvib_2

def get_diff(values1, values2):
  """Calculate the Fvib difference with respect to values1"""

  diff = []

  for i in range(0, len(values1)):
    if values1[i] == "undefined" or values2[i] == "undefined":
      continue
    else:
      values_diff = abs(values2[i] - values1[i])
      diff.append(values_diff)

  return diff

def plot_histogram(diff):
  """Plot the errors as a histogram plot"""

  fig, ax = plt.subplots(figsize=(10,7))
  ax.set_xlim(0.0, 2.5)
  ax.set_xlabel(r"error in $F_{vib}$ $(eV)$", fontsize=12)
  ax.set_ylabel("relative frequency", fontsize=12)

  counts, bins = np.histogram(diff, bins=100)
  frequencies = counts / float(counts.sum())
  ax.bar(bins[:-1], frequencies, width=np.diff(bins), align="edge", color="darkred", label="MACE-MP-agnesi-medium model")
  ax.legend(fontsize=12)

  plt.savefig(fname="error_histogram.png", format="png")
  plt.show()

def plot_scatter(fvib_1, fvib_2, diff, filename1, filename2):
  """Plot the two data lists against each other"""

  # Remove all undefined values from fvib_1
  undefined_indeces = []
  for i in range(0, len(fvib_1)):
    if fvib_1[i] == "undefined":
      undefined_indeces.append(i)

  print(f"{len(undefined_indeces)} undefined values were detected in {filename1}. They will be removed in the scatter plot.")
  for i in range(0, len(undefined_indeces)):
    fvib_1.pop(undefined_indeces[i]-i)
    fvib_2.pop(undefined_indeces[i]-i)

  # Remove all undefined values from fvib_2
  undefined_indeces = []
  for i in range(0, len(fvib_2)):
    if fvib_2[i] == "undefined":
      undefined_indeces.append(i)

  print(f"{len(undefined_indeces)} undefined values were detected in {filename2}. They will be removed in the scatter plot.")
  for i in range(0, len(undefined_indeces)):
    fvib_2.pop(undefined_indeces[i]-i)
    fvib_1.pop(undefined_indeces[i]-i)

  fig, ax = plt.subplots(figsize=(10,7))
  ax.set_xlabel(r"$F_{vib}$ VASP PBE D3-BJ $(eV)$", fontsize=12)
  ax.set_ylabel(r"$F_{vib}$ MLIP model $(eV)$", fontsize=12)

  # Transform fvib_1 and fvib_2 to numpy arrays
  fvib_1 = np.array(fvib_1)
  fvib_2 = np.array(fvib_2)

  ax.scatter(fvib_1, fvib_2, color="blue", label=r"$F_{vib}$ data points")

  # Perform a linear regression for the two datasets
  print("")
  print("Performing a linear regression for the two datasets ...")
  res = stats.linregress(fvib_1, fvib_2)
  print(f"R-squared: {res.rvalue**2:.5f}")
  print(f"f(x) = {res.slope:.5f}*x + {res.intercept:.5f}")
  # Print the linear regression function
  ax.plot(fvib_1, res.intercept + res.slope*fvib_1, color='r', label="linear fit")

  # Calculate the mean absolute error (MAE) for the two datasets
  diff = np.array(diff)
  mae = diff.sum()/len(diff)
  print(f"MAE: {mae:.5f}")

  # Show the R^2 and MAE value in the plot
  ax.text(
    0.80, 0.10,
    f"$R^2 = {res.rvalue**2:.5f}$\nMAE = {mae:.5f}",
    transform=ax.transAxes,          # relative position with respect to x,y axis
    fontsize=12,
    verticalalignment='top',
    bbox=dict(boxstyle="round", facecolor="white", alpha=0.7)  # box around the text
  )

  ax.legend(fontsize=12)
  plt.savefig(fname="fvib_scatter.png", format="png")
  plt.show()

################################################################################
##### MAIN PROGRAM
################################################################################

if __name__ == "__main__":
  # Print information about basic usage of the script
  print("Script compare_fvib.py: Compare the results from different MLIPs")
  print("                        written by benchmarking.sh -make_vDOS")
  print("")
  print("Usage: compare_fvib.py [fvib_results_1.dat] [fvib_results_2.dat]")
  print("")

  # Get the two files to compare
  filename1, filename2 = process_command_line_arguments()

  # Exit if no files were given
  if not filename1 or not filename2:
    print("No suitable files to compare given ...")
    print("")
    sys.exit(-1)

  # Print information about the files given
  print("Files to compare given:")
  print("Reference file:", filename1)
  print("Second file:", filename2)
  print("")

  # Get the average Fvib values
  fvib_1, fvib_2 = get_fvib(filename1, filename2)

  print("")
  print(f"{filename1} has {len(fvib_1)} entries")
  print(f"{filename2} has {len(fvib_2)} entries")
  print("")

  # Calculate the Fvib difference for all runs
  fvib_diff = get_diff(fvib_1, fvib_2)

  # Plot the error per atom as histogram plot
  print("Generating the histogram plot ...")
  plot_histogram(fvib_diff)
  # Plot the values against each other
  print("Generating the scatter plot ...")
  plot_scatter(fvib_1, fvib_2, fvib_diff, filename1, filename2)


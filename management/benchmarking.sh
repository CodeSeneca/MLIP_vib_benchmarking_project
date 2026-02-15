#######################################################################
#####
##### MASTER SCRIPT TO MANAGE THE MLIP VIB BENCHMARKING PROJECT
##### Author: Maximilian Bechtel <maxi.bechtel@fau.de>
#####
#######################################################################

#######################################################################
##### GLOBAL VARIABLES
#######################################################################

# The mode of the program: -make_vdos or -copy_input
mode="undefined"
# The number of run folders
num_runs="not set"
# The number of traj$j folders per run$i folder
num_trajs="not set"
# Path to the POSCAR files, where to copy from
path="not set"
# Start run folder to be submitted
start_run="not set"
# End run folder to be submitted
end_run="not set"
# Format of the trajectory file: XDATCAR or trajectory.xyz
format="not set"

# Temperature in Kelvin
temp="not set"
# Output file
out_file="fvib_results.dat"

# The copy command via rsync
CP="rsync -auvz"

# The location at which the bash script is running
script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

#######################################################################
##### FUNCTIONS
#######################################################################

# Parse command line arguments
parse_command_line_args() {
  # The number of command line arguments given
  local num_arguments=$#

  if [ $num_arguments -eq 0 ]; then
    echo ""
    echo "Usage: $0 [options] "
    echo ""
    echo "***** The options are: *****"
    echo "-copy_input [number of run folders] [number of traj folders per run folder] [path where to copy POSCAR files from]"
    echo "-make_vDOS [start] [end] [temperature] [fileformat]"
    echo "-submit_jobs [start] [end]"
    echo "-render_molecules [start] [end]"
    echo "-stability [start] [end]"
    echo ""
    exit -1
  fi

  if [ "$1" == "-copy_input" ]; then
    mode="copy_input"
    num_runs=$2
    num_trajs=$3
    path="$4"
  elif [ "$1" == "-make_vDOS" ]; then
    mode="make_vDOS"
    start_run=$2
    end_run=$3
    temp=$4
    format=$5
  elif [ "$1" == "-submit_jobs" ]; then
    mode="submit_jobs"
    start_run=$2
    end_run=$3
  elif [ "$1" == "-render_molecules" ]; then
    mode="render_molecules"
    start_run=$2
    end_run=$3
  elif [ "$1" == "-stability" ]; then
    mode="stability"
    start_run=$2
    end_run=$3
  fi
}

# Copy all POSCAR files from the OMOL25 database subset
# all POSCAR files will be stored in folders of the form run$i, each run$i
# folder consists of folders of the form traj$j
# the number of traj$j folders per run$i folder has to be given as argument
copy_input() {
  if [ -z $num_runs ] || [ -z $num_trajs ] || [ -z $path ]; then
    echo "No suitable value for number of run folders, traj folders or the path was given. Aborting ..."
    exit -3
  fi

  echo "Number of molecules (run folders): $num_runs"
  echo "Number of individual MD runs (traj folders) per run folder: $num_trajs"
  echo "Path, where to copy POSCAR files from: $path"
  echo ""

  echo "Creating all run folders ..."
  # Create all run$i folders with traj$j folders each
  # if they already exist, all run$i folders will be overwritten
  if [ -d "run1" ]; then
      echo "The run folders already exist. They will be overwritten ..."
  fi

  for i in $(seq 1 $num_runs); do
    if [ -d "run$i" ]; then
      rm -r "run$i"
    fi
    mkdir "run$i"
  done

  # Create each traj$j folder and copy POSCAR to it
  for i in $(seq 1 $num_runs); do
    cd "run$i"
    echo "Copying POSCAR to run$i ..."

    $CP "$path/run$i/traj1/POSCAR" . > copy.log

      for j in $(seq 1 $num_trajs); do
        mkdir "traj$j"
        cp POSCAR "traj$j"
      done
      rm POSCAR
      rm copy.log

    cd ..
  done
  
}

# Determine the number of atoms in a POSCAR file
# $1: filename (get_number_of_atoms $1)
get_number_of_atoms() {
  local filename="$1"
  local num_atoms_array=($(head -n 7 $filename | tail -n 1))

  # Add up all numbers in num_atoms_array
  local sum=0
  for atom_number in "${num_atoms_array[@]}"
  do
    sum=$(($sum+$atom_number))
  done

  echo "$sum"
}

# Go into each run$i folder to create a vDOS spectrum and calculate
# the vibrational free energy F_vib from it
make_vDOS() {
  if [ -z $start_run ] || [ -z $end_run ] || [ -z $temp ] || [ -z $format ]; then
    echo "No suitable values were given. Aborting ..."
    echo ""
    exit -4
  fi

  echo "A vDOS calculation for the following run folders will be done:"
  echo "Trajectory format: $format"
  echo "Start run folder: $start_run"
  echo "End run folder: $end_run"
  echo "Temperature: $temp K"
  echo ""

  # Delete output file from previous run
  if [ -f $out_file ]; then
    rm $out_file
  fi
  
  for i in $(seq $start_run $end_run)
  do
    # Average vibrational free energy
    fvib_aver=0
  
    echo "Creating vDOS spectra in run$i ..."
    echo "# Molecule type $i" >> $out_file
  
    cd run$i

    if [[ "$format" == "trajectory.xyz" ]]; then
      xyz2poscar.py traj1/start.xyz traj1/POSCAR
    fi

    elements=$(head -n 6 traj1/POSCAR | tail -n 1)
    echo "# Elements: $elements" >> "../fvib_results.dat"
  
    # Determine the number of traj$i folder
    traj_number=$(find . -maxdepth 1 -type d -name "traj*" | wc -l)
    for j in $(seq 1 $traj_number)
    do
      cd traj$j
  
        # Remove old travis.log and power_spectrum_global.csv file
        if [ -f "travis.log" ]; then
          rm "travis.log"
        fi
        if [ -f "power_spectrum_global.csv" ]; then
          rm "power_spectrum_global.csv"
        fi
  
        if [[ "$format" == "XDATCAR" ]]; then
          # Create a new XDATCAR file in the NpT format
          modify_xdatcar -print_npt > modify_xdatcar.log
          rm modify_xdatcar.log
          # Rename old XDATCAR to XDATCAR_NVT
          mv XDATCAR XDATCAR_NVT
          # Rename XDATCAR (in NpT format) to XDATCAR for subsequent analysis
          # with TRAVIS
          mv XDATCAR_mod XDATCAR
        fi

        ###########################
        # HERE ANALYSIS WITH TRAVIS
        ##########################
        if [[ "$format" == "XDATCAR" ]]; then
          travis -p XDATCAR -i ../../travis_vDOS_input.txt > power_spectrum.log
        elif [[ "$format" == "trajectory.xyz" ]]
        then
          travis -p trajectory.xyz -i ../../travis_vDOS_input.txt > power_spectrum.log
        fi

        rm power_spectrum.log

        if [[ "$format" == "XDATCAR" ]]; then
          # Undo the renaming
          mv XDATCAR XDATCAR_mod
          mv XDATCAR_NVT XDATCAR
          rm xdat_length
        fi
  
        # Determine the number of atoms in the POSCAR file
        if [[ "$format" == "trajectory.xyz" ]]; then
          xyz2poscar.py start.xyz POSCAR
        fi
        num_atoms=$(get_number_of_atoms "POSCAR")
        # Determine the number of vibrational degrees of freedom fvib = 3N - 6
        num_fvib=$(echo "3*$num_atoms-6" | bc -l)
        # Calculate the vibrational free energy F_vib
        python $script_dir/calc_avib.py -T=$temp -fvib=$num_fvib > fvib.log
        avib_out=$(tail -n 2 fvib.log | head -n 1)
        avib=$(echo "$avib_out" | awk '{print$3}')
        echo "$num_fvib    $avib" >> "../../fvib_results.dat"
  
      fvib_aver=$(echo "scale=5; $fvib_aver+$avib" | bc -l)
      cd ..
    done
  
    fvib_aver=$(echo "scale=5; $fvib_aver/$traj_number" | bc -l)
    cd ..
    echo "# --------------------" >> $out_file
    echo "     $fvib_aver" >> $out_file
    echo "" >> $out_file
  done
}

# Go into each run$i/traj$j folder, copy input files and submit job
submit_jobs() {
  if [ -z $start_run ] || [ -z $end_run ]; then
    echo "No suitable values were given. Aborting ..."
    echo ""
    exit -4
  fi

  echo "The following runs will be submitted:"
  echo "Start run folder: $start_run"
  echo "End run folder: $end_run"

  # Go into each run$i folder and copy input files to each traj$j folder
  echo ""
  echo "Copying slurm_script to all given run folders ..."
  for i in $(seq $start_run $end_run); do
    cd "run$i"

    # Determine the number of traj$i folder
    traj_number=$(find . -maxdepth 1 -type d -name "traj*" | wc -l)

    for j in $(seq 1 $traj_number); do
      cp ../slurm_script traj$j

      # Rename the job name in the slurm_script to frame$i.$j
      cd traj$j
      sed -i "5a\\#SBATCH --job-name=run$i.$j" slurm_script
      sed -i '5d' slurm_script
      cd ..
    done

    cd ..
  done

  # Go into each run$i/traj$j folder and submit job
  for i in $(seq $start_run $end_run); do
    echo "Submitting jobs in run$i ..."
    cd "run$i"
    for j in $(seq 1 $traj_number); do
      cd "traj$j"
      sbatch slurm_script
      sleep 0.5
      cd ..
    done
    cd ..
  done
}

# Go into each run$i folder and render pictures of the POSCAR present there
render_molecules() {
  echo "Pictures for the following run folders will be rendered:"
  echo "Start run folder: $start_run"
  echo "End run folder: $end_run"
  echo ""

  if [ -d "rendered_pictures" ]; then
      echo "The rendered_pictures folder already exists. New rendered pictures will be appended but some pictures might be overwritten!"
  else
      mkdir rendered_pictures
  fi


  # Go into each run$i folder
  for i in $(seq $start_run $end_run)
  do
    cd run$i
    echo "Rendering pictures for run$i ..."

    cd traj1
      # Convert the POSCAR in .xyz format -> gives the file output.xyz
      python3 $script_dir/poscar2xyz.py >> poscar2xyz.tmp
      rm poscar2xyz.tmp

      # Rotate the molecule in output.xyz to standard rotation
      python3 $script_dir/rotate_mol.py output.xyz >> rotate_mol.tmp
      rm rotate_mol.tmp

      # The actual rendering with VMD
      $script_dir/render_mol.sh struc_rotated.xyz >> render_mol.tmp
      rm render_mol.tmp
      rm tmp.vmd
      cp struc_rotated_top.png ../../rendered_pictures/molecule$i.png
      mv output.xyz molecule$i.xyz
    cd ..

    cd ..
  done

  if [ -z $start_run ] || [ -z $end_run ]; then
    echo "No suitable values were given. Aborting ..."
    echo ""
    exit -4
  fi
}

determine_stability() {
  if [ -z $start_run ] || [ -z $end_run ]; then
    echo "No suitable values were given. Aborting ..."
    echo ""
    exit -4
  fi

  echo "The stability of all MD runs will be evaluated:"
  echo "Start run folder: $start_run"
  echo "End run folder: $end_run"
  echo ""

  # Flag if MD trajectory is stable or not
  # 0 = stable, 1 = unstable
  local stable=0
  # Name of the stability output file
  local outfile="stability.dat"

  # Remove file stability.dat if it exists
  if [ -f "$outfile" ]; then
    echo "The file $outfile already exists. It will be overwritten ..."
    rm $outfile
  fi

  echo "# A boolean mask is applied with" >> stability.dat
  echo "# 0 = stable" >> stability.dat
  echo "# 1 = unstable" >> stability.dat
  echo "" >> stability.dat

  for i in $(seq $start_run $end_run)
  do
    cd run$i
    echo "Evaluating stability in run$i ..."
    echo "# Molecule type $i:" >> ../stability.dat

    # Determine the number of traj$i folder
    traj_number=$(find . -maxdepth 1 -type d -name "traj*" | wc -l)
    for j in $(seq 1 $traj_number)
    do
      cd traj$j

      # Determine the stability of all bonds in the MD trajectory
      analyze_md -stability > stability.log
      sleep 1
      rm xdat_length

      # Check if the folder stability is empty
      # yes: no bonds have broken
      if [ -z "$(ls -A "analysis")" ]; then
        echo "MD trajectory is stable."
        stable=0
      else
        echo "MD trajectory is not stable."
       stable=1
      fi
      echo "$stable" >> ../../stability.dat

      cd ..
    done

    cd ..
  done
}

#######################################################################
##### MAIN PROGRAM
#######################################################################

# Parse the command line arguments given by the user
parse_command_line_args "$@"
echo ""
echo "The following mode was chosen: $mode"

if [ "$mode" == "copy_input" ]; then
  copy_input
elif  [ "$mode" == "make_vDOS" ]; then
  make_vDOS
elif  [ "$mode" == "submit_jobs" ]; then
  submit_jobs
elif  [ "$mode" == "render_molecules" ]; then
  render_molecules
elif  [ "$mode" == "stability" ]; then
  determine_stability
else
  echo "No suitable mode was given. Aborting benchmarking.sh ..."
  echo ""
  exit -2
fi

echo ""
echo "********** benchmarking.sh exited normally **********"
echo ""


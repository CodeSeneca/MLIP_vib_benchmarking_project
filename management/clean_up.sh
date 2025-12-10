#! /usr/bin/bash

#############################################################
##### GO INTO EACH RUN FOLDER AND CLEAN UP UNNECESSARY FILES
#############################################################

for i in {1..1000}
do
  echo "Cleaning up run$i ..."
  cd run$i

  for j in {1..3}
  do 
    cd traj$j

    rm run*
    rm slurm_script
    rm CONTCAR
    rm XDATCAR_mod
    rm fvib.log
    rm travis.log

    cd ..
  done

  cd ..
done

#!/bin/bash

#script to change text inside numbered script to send to SLURM

for i in 21 25 33 # OH density of the SAM 
do
  for j in 216 1000 2000 3000 4000 5000 6500 7000 8000 9000 10000  # number of water molecules
  do
  	#cd /gfs1/work/beclaila/s${i}_w${j}
  	cd /Volumes/BACKUPS/MyTimeMachine/HLRN/s${i}_w${j}/
  	
  	cp ${i}pc_${j}.top ${i}pc_${j}_cuda.top 
  	
    file=${i}pc_${j}_cuda.top
    sed -i '' 's/ffG53a6\.itp/\/home\/leixeres\/programs\/gromacs\/share\/gromacs\/top\/gromos53a6\.ff\/forcefield\.itp/g' $file
    sed -i '' 's/spce\.itp/\/home\/leixeres\/programs\/gromacs\/share\/gromacs\/top\/gromos53a6\.ff\/spce\.itp/g' $file
    
    
   done
done
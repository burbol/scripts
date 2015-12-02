#!/bin/bash

source /usr/local/gromacs/bin/GMXRC

for i in 0 5 11 17 # OH density of the SAM 
do
  for j in 1000 2000 3000 4000  # number of water molecules
  do
  
   echo 0 20000 |/home/shavkat/GMX/bin/trjcat -f NVT_sam${i}_water${j}.xtc NVT2_sam${i}_water${j}.xtc -o NVT_sam${i}_water${j}_MERGED.xtc -settime 
  

	done        
  done
#!/bin/bash

# Script for using trjconv recursively 

source /usr/local/gromacs/bin/GMXRC

cd /Users/burbol/Desktop/Testing/files_for_laila/sam0
  for i in 3000 4000 # number of water molecules
  do

    for j in 0 2 4 6 8 10 12 14 16 18  # segments of time in ns
    do

        k=$((j+2))
        
        echo 1 4 4 | g_densmap -f NVT_surf_water${i}_0.xtc -s NVT_surf_water${i}_0.tpr -n index0_${i}.ndx -amax 5 -rmax 5 -bin 0.05 -od densmap_0pc_${i}_${j}ns_${k}ns.dat -b ${j}000 -e ${k}000
        
       

done
done

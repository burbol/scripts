#!/bin/bash

# Script for using trjconv recursively 

source /usr/local/gromacs/bin/GMXRC

for i in 66 # OH density of the SAM 
do
  for j in 4000  # number of water molecules
  do
  mkdir sam${i}_water${j}
  mkdir densmaps_s${i}_w${j}
    for k in 2 4 6 8  # segments of time in ns
    do
        start=$((k*1000))
        ending=$((start+2000))
        l=$((k+2))
        
        cp ../NVT_sam${i}_water${j}.tpr .
        
        echo "Creating a dummy.gro file to visualize results using the dummy.xtc file" >> output${k}.dat
        
        echo 0 | trjconv -f dummy_${k}ns_${l}ns.xtc -s NVT_sam${i}_water${j}.tpr -trans ${x} ${y} 0 -o dummy_${k}ns_${l}ns.gro -pbc atom -b ${start} -e ${ending}

        cp dummy_${k}ns_${l}ns.gro ./sam${i}_water${j}
        rm dummy_${k}ns_${l}ns.gro

done
done
done

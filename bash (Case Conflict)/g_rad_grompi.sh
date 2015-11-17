#!/bin/bash

# Script for using g_density recursively and g_rad_density in time intervals of 500 ps

#source /usr/local/gromacs/bin/GMXRC

n=4
cd /home/eixeres/files_for_laila
mkdir ./global_density_maps
mkdir ./g_rad_densmaps

for i in 5 11 17 66 # OH density of the SAM 
do
  for j in 1000 2000 3000 4000 # number of water molecules
  do
  
  if [$i -eq 0 ]
  then 
    n=3
  else 
    n=4
  fi
  
#    /home/shavkat/GMX/bin/make_ndx -f NVT_sam${i}_water${j}.gro -o index${i}_${j}.ndx

#    /home/eixeres/gromacs_tpi_compiled/bin/grompp -f NVT.mdp -c NVT_sam${i}_water${j}.gro -p ${i}pc1_${j}.top -n index${i}_${j}.ndx -o g_rad_NVT_sam${i}_water${j}.tpr -maxwarn 1

echo ${n} | /home/shavkat/GMX/bin/g_density -f NVT_sam${i}_water${j}.xtc -s NVT_sam${i}_water${j}.tpr -o g_density_NVT_sam${i}_water${j}.xvg -sl 1000

	    
    for k in {0..195..5}  # segments of time to be multiplied by 100 so that we get [ps]
    do
        nanosecs = k/10
        start=$((k*100))
        ending=$((start+500))

        echo "Creating densmap g_rad_dmap_$i pc_w$j _$start ps_$ending ps.dat"
        
        echo ${n} ${n} | /home/eixeres/g_rad_density -f NVT_sam${i}_water${j}.xtc -s NVT_sam${i}_water${j}.tpr -n index${i}_${j}.ndx -max 5 -min -5 -sz 200 -o g_rad_dmap_${i}pc_w${j}_${start}ps_${ending}ps.xvg -b ${start} -e ${ending}

        
        cp g_rad_dmap_${i}pc_w${j}_${start}ps_${ending}ps.xvg ./g_rad_densmaps/
        rm g_rad_dmap_${i}pc_w${j}_${start}ps_${ending}ps.xvg
	
	cp g_density_NVT_sam${i}_water${j}.xvg ./g_density/
	rm g_density_NVT_sam${i}_water${j}.xvg
        
    done
  done
done

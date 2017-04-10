#!/bin/bash

#################################
# For testing:

#for i in 0 
#do
#	for j in 1000
#    do
#################################

cd /Volumes/Version_v2

for i in 11 22 33 37 44 50
do
   for j in 1000 2000 3000 4000 5000 6500 7000 8000 9000 10000 
   do
   
   cd  s${i}_w${j}/
  
  	# copy to sheldon-ng
  	#scp s${i}_w${j}/* eixeres@sheldon-ng.physik.fu-berlin.de:/scratch/eixeres/Version_v2/scripts/s${i}_w${j}/
  	
  	# copy to yoshi
  	#scp s${i}_w${j}/pairtypes.itp eixeres@yoshi.physik.fu-berlin.de:/scratch/eixeres/Version_v2/scripts/s${i}_w${j}/
  	#cp {${i}pc_${j}_old.top,index${i}_${j}.ndx} /scratch/eixeres/Version_v2/s${i}_w${j}/
  	
   # copy to iMac
   scp {eixeres@sheldon-ng.physik.fu-berlin.de:/net/data03/eixeres/scratch_yoshi/Version_v2/s${i}_w${j}/*dil.xtc,eixeres@sheldon-ng.physik.fu-berlin.de:/net/data03/eixeres/scratch_yoshi/Version_v2/s${i}_w${j}/NVT*.gro} .
  	
  	cd ..

  done
done
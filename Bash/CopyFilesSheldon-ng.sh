#!/bin/bash

#cd /Users/burbol/Desktop/scripts/Python/SCRIPT_CREATION/NewVersion1/
#cd /Users/burbol/Desktop/scripts/SAM_CREATION/SAMs/NEW/drop_placement/NewVersion1/
#cd /Users/burbol/Desktop/scripts/SAM_CREATION/SAMs/NEW/drop_placement/NewVersion1/ITPnew
#cd /Users/burbol/Desktop/scripts/Python/SCRIPT_CREATION/scheldon-ng/
#cd /Users/burbol/MEGAsync/scripts/Python/SCRIPT_CREATION/scheldon/

for i in 25 # OH density of the SAM 
do
	#for j in 2000 3000 4000 5000 6500 7000 8000 9000  # number of water molecules
	for j in 4000
	#for j in 6500 7000 8000 9000 10000   # cuda
	#for j in 216 1000 4000 3000 4000 5000 # sheldon-ng
  do 
  	#mkdir s${i}_w${j}
  	#cd /Users/burbol/MEGAsync/scripts/SAM_CREATION/SAMs/NEW/drop_placement/NewVersion3/s${i}_w${j}
  	
  	#cd /Users/burbol2/Dropbox/scripts/Python/SCRIPT_CREATION/SLURM/s${i}_w${j}
  	
  	#cd /home/eixeres/Downloads/s21_w1000
  	#scp * eixeres@sheldon-ng.physik.fu-berlin.de:/scratch/eixeres/NewVersion3/s${i}_w${j}/
  	#scp * eixeres@yoshi.physik.fu-berlin.de:/scratch/eixeres/NewVersion3/s${i}_w${j}/
  	mkdir /Volumes/UNI/Densmaps_NewVersion4/s${i}_w${j}/
  	scp eixeres@yoshi.physik.fu-berlin.de:/net/data/eixeres/NewVersion4/FINISHED/s${i}_w${j}/*.xvg /Volumes/UNI/Densmaps_NewVersion4/s${i}_w${j}/
  
  	#cp *.itp ../s${i}_w${j}/
  	#cp ./soroban/s${i}_w${j}/s${i}_w${j}* ./NewVersion1/s${i}_w${j}/
  	
  	#scp burbol@burbol.ddns.net:/Volumes/BACKUPS/MyTimeMachine/HLRN/s${i}_w${j}/*.top ./s${i}_w${j}
  	#scp burbol@burbol.ddns.net:/Volumes/BACKUPS/MyTimeMachine/HLRN/s${i}_w${j}/*.itp ./s${i}_w${j}
  	#scp burbol@burbol.ddns.net:/Volumes/BACKUPS/MyTimeMachine/HLRN/s${i}_w${j}/sam*.gro ./s${i}_w${j}
  	#scp burbol@burbol.ddns.net:/Volumes/BACKUPS/MyTimeMachine/HLRN/s${i}_w${j}/NVT* ./s${i}_w${j}
  	#scp burbol@burbol.ddns.net:/Volumes/BACKUPS/MyTimeMachine/HLRN/s${i}_w${j}/NVT* ./s${i}_w${j}
  	
  	#scp burbol@burbol.ddns.net:/Users/burbol/Desktop/scripts/Python/SCRIPT_CREATION/scheldon-ng/s${i}_w${j}/s${i}_w${j}_* ./s${i}_w${j}
  	#scp sam${i}_water${j}.gro eixeres@bert.physik.fu-berlin.de:/home/eixeres/Downloads/
  	

  done
done
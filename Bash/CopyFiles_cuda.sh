#!/bin/bash

for i in 0 # OH density of the SAM 
do
	#for j in 216 1000 2000 3000 4000 5000 6500 7000 8000 9000 10000  # number of water molecules
	#for j in 6500 7000 8000 9000 10000   # cuda
	for j in 6500 

  do 
  	#mkdir s${i}_w${j}
  
  	#scp burbol@burbol.ddns.net:/Volumes/BACKUPS/MyTimeMachine/HLRN/s${i}_w${j}/*.top ./s${i}_w${j}
  	#scp burbol@burbol.ddns.net:/Volumes/BACKUPS/MyTimeMachine/HLRN/s${i}_w${j}/*.itp ./s${i}_w${j}
  	#scp burbol@burbol.ddns.net:/Volumes/BACKUPS/MyTimeMachine/HLRN/s${i}_w${j}/sam*.gro ./s${i}_w${j}
  	
  	#scp burbol@burbol.ddns.net:/Volumes/BACKUPS/MyTimeMachine/HLRN/s${i}_w${j}/NVT* ./s${i}_w${j}
  	
  	#scp burbol@burbol.ddns.net:/Users/burbol/Desktop/NVT_20ns_cuda.mdp ./s${i}_w${j}
  	
  	#scp -r /Users/burbol/MEGAsync/scripts/Python/SCRIPT_CREATION/NewVersionBIG/s${i}_w${j} leixeres@login.megware.de:/home/leixeres/
  	#scp /Users/burbol/MEGAsync/scripts/Python/SCRIPT_CREATION/NewVersionBIG/s${i}_w${j}/s${i}_w${j}_cuda* leixeres@login.megware.de:/home/leixeres/s${i}_w${j}/
  	#scp -r /Users/burbol/MEGAsync/scripts/SAM_CREATION/SAMs/NEW/drop_placement/NewVersion3/s${i}_w${j} eixeres@sheldon-ng.physik.fu-berlin.de:/scratch/eixeres/NewVersionBIG/
  	#scp /Users/burbol/MEGAsync/scripts/Python/SCRIPT_CREATION/sheldon-ng/s${i}_w${j}/s${i}_w${j}* eixeres@sheldon-ng.physik.fu-berlin.de:/scratch/eixeres/NewVersionBIG/s${i}_w${j}/
  	scp /Users/burbol/MEGAsync/scripts/Python/SCRIPT_CREATION/sheldon/s${i}_w${j}_next/s${i}_w${j}_next* eixeres@sheldon.physik.fu-berlin.de:/home/eixeres/files_for_laila/s${i}_w${j}/
  	#scp /Users/burbol/MEGAsync/scripts/SAM_CREATION/SAMs/NEW/drop_placement/NewVersion3/s${i}_w${j}/*.top eixeres@sheldon-ng.physik.fu-berlin.de:/scratch/eixeres/NewVersionBIG/s${i}_w${j}
  	#scp /Users/burbol/MEGAsync/scripts/SAM_CREATION/SAMs/NEW/drop_placement/NewVersion3/s${i}_w${j}/*.gro eixeres@sheldon-ng.physik.fu-berlin.de:/scratch/eixeres/NewVersionBIG/s${i}_w${j}
  	

  done
done
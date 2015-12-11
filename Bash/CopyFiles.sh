#!/bin/bash

#cd /Users/burbol/Desktop/scripts/Python/SCRIPT_CREATION/NewVersion1/
#cd /Users/burbol/Desktop/scripts/SAM_CREATION/SAMs/NEW/drop_placement/NewVersion1/
#cd /Users/burbol/Desktop/scripts/SAM_CREATION/SAMs/NEW/drop_placement/NewVersion1/ITPnew
#cd /Users/burbol/Desktop/scripts/Python/SCRIPT_CREATION/scheldon-ng/
cd /Volumes/Backup/YosemiteFiles/MEGAsync/scripts/SAM_CREATION/SAMs/NEW/drop_placement/NewVersion3/

#{216:60, 1000:80, 2000:80, 3000:80, 4000:80, 5000:100, 6500:100, 7000:100, 8000:120, 9000:120, 10000:120}

for i in 21 25 # OH density of the SAM 
do
	for j in 1000 2000 3000 4000 5000 6500 7000 8000 9000 10000  # number of water molecules
	#for j in 6500 7000 8000 9000 10000   # cuda
	#for j in 1000 4000 3000 4000 5000 # sheldon-ng
	#for j in 216 2000 4000  #
  do 
  	#mkdir ${i}_w${j}
  	#cp s${i}/${i}pc_${j}.top s${i}_w${j}/
  
  	#cp *.itp ../s${i}_w${j}/
  	#cp ./soroban/s${i}_w${j}/s${i}_w${j}* ./NewVersion1/s${i}_w${j}/
#   	scp burbol@burbol.ddns.net:/Volumes/BACKUPS/MyTimeMachine/HLRN/s${i}_w${j}/*.top ./s${i}_w${j}
#   	scp burbol@burbol.ddns.net:/Volumes/BACKUPS/MyTimeMachine/HLRN/s${i}_w${j}/*.itp ./s${i}_w${j}
#   	scp burbol@burbol.ddns.net:/Volumes/BACKUPS/MyTimeMachine/HLRN/s${i}_w${j}/sam*.gro ./s${i}_w${j}
#   	scp burbol@burbol.ddns.net:/Volumes/BACKUPS/MyTimeMachine/HLRN/s${i}_w${j}/NVT* ./s${i}_w${j}
#   	scp burbol@burbol.ddns.net:/Volumes/BACKUPS/MyTimeMachine/HLRN/s${i}_w${j}/NVT* ./s${i}_w${j}
#   	
#   	scp burbol@burbol.ddns.net:/Users/burbol/Desktop/scripts/Python/SCRIPT_CREATION/scheldon-ng/s${i}_w${j}/s${i}_w${j}_* ./s${i}_w${j}	

	echo "hi"
  done
done
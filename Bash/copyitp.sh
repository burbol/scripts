#!/bin/bash

## OLD
## cd /Users/burbol/MEGAsync/scripts/SAM_CREATION/SAMs/NEW/drop_placement/ #folder where the gromacs files are saved together
#cd /Users/burbol/MEGAsync/scripts/Python/SCRIPT_CREATION/sheldon
cd /Users/burbol/MEGAsync/scripts/SAM_CREATION/SAMs/NEW/drop_placement/NewVersion3/NewVersionBIG_Backup
cd /Users/burbol/MEGAsync/scripts/SAM_CREATION/SAMs/NEW/drop_placement/NewVersion3

# folders where the sripts are should already exist!!!
# here: /Users/burbol/MEGAsync/scripts/Python/SCRIPT_CREATION


for i in 21 25 # OH density of the SAM 
do
  #for j in 216 1000 2000 3000 4000 5000 6500 7000 8000 9000 10000  # number of water molecules
  for j in 6500 7000 8000 9000 10000 # For testing!!!!
  do
	
	#mkdir s${i}_w${j}
	#cp ./Wrong_ITP/s${i}_w${j}/*.itp ./s${i}_w${j}/
	#cp ./Wrong_ITP/s${i}_w${j}/*.mdp ./s${i}_w${j}/
	#cp ./Wrong_ITP/s${i}_w${j}/*.gro ./s${i}_w${j}/
	#cp ./Wrong_ITP/s${i}_w${j}/*.top ./s${i}_w${j}/
	
	#cp /Users/burbol/Downloads/*.itp  s${i}_w${j}/
	
	#cp sam${i}_water${j}.gro ../s${i}_w${j}/
	#cp ${i}pc_${j}.top ../s${i}_w${j}/
	
	mv sam${i}_water${j}.gro s${i}/sam${i}_water${j}.gro
	mv ${i}pc_${j}.top s${i}/${i}pc_${j}.top

	
	
	## copy .itp files 
	##cp NewVersion3/s${i}_w${j}/*_cuda.top /Users/burbol/MEGAsync/scripts/Python/SCRIPT_CREATION/NewVersionBIG/s${i}_w${j}/
	##cp NewVersion3/s${i}_w${j}/*_cuda.top /Users/burbol/MEGAsync/scripts/SAM_CREATION/SAMs/NEW/drop_placement/NewVersion3/s${i}_w${j}/
	
	## cp Mini_cuda.mdp /Users/burbol/MEGAsync/scripts/Python/SCRIPT_CREATION/NewVersionBIG/s${i}_w${j}/
	##	 cp Mini_cuda.mdp /Users/burbol/MEGAsync/scripts/SAM_CREATION/SAMs/NEW/drop_placement/NewVersion3/s${i}_w${j}/
	

  done
done
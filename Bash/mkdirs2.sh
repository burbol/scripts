#!/bin/bash

#cd /Volumes/Backup/YosemiteFiles/MEGAsync/scripts/SAM_CREATION/SAMs/NEW/drop_placement/DropsGroTop
cd /Volumes/Backup/YosemiteFiles/MEGAsync/scripts/SAM_CREATION/SAMs/NEW/drop_placement/NewVersion4

for i in 5 21 25 41 # OH density of the SAM 
do
  for j in 216 1000 2000 3000 4000 5000 6500 7000 8000 9000 10000  # number of water molecules
  do 
	
	mkdir ./s${i}_w${j}
	#cp sam${i}_water${j}.gro s${i}_w${j}/
	#cp *.itp s${i}_w${j}/
	#cp *.mdp s${i}_w${j}/
  
  done
done
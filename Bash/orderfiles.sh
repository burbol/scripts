#!/bin/bash

# THIS SCRIPT PUTS THE FILES NEEDED TO START A GROMACS SIMULATION IN ORDERED FOLDERS
# AND THEN COPIES THEM TO THE FOLDER,
# WHERE THE CORRESPONDING SUBMISSION SCRIPTS ARE 

cd /Users/burbol/MEGAsync/scripts/SAM_CREATION/SAMs/NEW/drop_placement/NewVersion3  #folder where the gromacs files are saved together

# folders where the sripts are should already exist!!!
# here: /Users/burbol/MEGAsync/scripts/Python/SCRIPT_CREATION

#mkdir /Users/burbol/MEGAsync/scripts/Python/SCRIPT_CREATION/NewVersionBIG

for i in 21 25 # OH density of the SAM 
do
  #for j in 216 1000 2000 3000 4000 5000 6500 7000 8000 9000 10000  # number of water molecules
  for j in 1000 2000 3000 4000 5000 # For testing!!!!
  do
  	#cd /Users/burbol/MEGAsync/scripts/SAM_CREATION/SAMs/NEW/drop_placement/NewVersionBIG/s${i}  #folder where the gromacs files are saved together
	mkdir s${i}_w${j}
	#mkdir /Users/burbol/MEGAsync/scripts/Python/SCRIPT_CREATION/NewVersionBIG/s${i}_w${j}/
  
	# copy .gro files
	#cp sam${i}_water${j}.gro /Users/burbol/MEGAsync/scripts/Python/SCRIPT_CREATION/NewVersionBIG/s${i}_w${j}/
	cp sam${i}_water${j}.gro ./s${i}_w${j}/
  
	# copy .top files 
	#cp ${i}pc_${j}.top /Users/burbol/MEGAsync/scripts/Python/SCRIPT_CREATION/NewVersionBIG/s${i}_w${j}/
	cp ${i}pc_${j}.top ./s${i}_w${j}/
	
	# copy .itp files 
	#cp *.itp /Users/burbol/MEGAsync/scripts/Python/SCRIPT_CREATION/NewVersionBIG/s${i}_w${j}/
	#cp *.itp ./s${i}_w${j}/
	
	# copy .mdp files 
	#cp *.mdp /Users/burbol/MEGAsync/scripts/Python/SCRIPT_CREATION/NewVersionBIG/s${i}_w${j}/
	#cp *.mdp ./s${i}_w${j}/
  
  done
done
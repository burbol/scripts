#!/bin/bash

# THIS SCRIPT PUTS THE FILES NEEDED TO START A GROMACS SIMULATION IN ORDERED FOLDERS
# AND THEN COPIES THEM TO THE FOLDER,
# WHERE THE CORRESPONDING SUBMISSION SCRIPTS ARE 

#cd /Users/burbol/Downloads/small_sams2  #folder where the gromacs files are saved together
#cd /Users/burbol/MEGAsync/scripts/Python/SCRIPT_CREATION/sheldon/ptensor/
#cd /Users/burbol/MEGAsync/scripts/Python/SCRIPT_CREATION/soroban/ptensor/

# folders where the sripts are should already exist!!!
# here: /Users/burbol/MEGAsync/scripts/Python/SCRIPT_CREATION

for i in 0 5 11 17 21 25 33 50 66 # OH density of the SAM 
do
	#mkdir w${i}
	#cd /Users/burbol/Downloads/small_sams2/w${i}
	#cd w${i}
	#cd /Users/burbol/MEGAsync/scripts/Python/SCRIPT_CREATION/soroban/ptensor/w${i}
	
	#cd /Users/burbol/MEGAsync/scripts/Python/SCRIPT_CREATION/sheldon/ptensor/w${i}
	#pwd
	
	#mkdir sheldon
	# copy .gro files
	#cp water${i}_ptensor.gro ./w${i}/
	#cp ./w${i}/water${i}_ptensor.gro ..
	#cp sam${i}_water_ptensor.gro ./w${i}/
	#cp ./w${i}/sam${i}_water_ptensor.gro ..
  
	# copy .top files 
	#cp water${i}_ptensor.top ./w${i}/
	#cp sam${i}_water_ptensor.top ./w${i}/
	#rm ./w${i}/water${i}_ptensor.gro
	#rm ./w${i}/water${i}_ptensor.top
	#rm ./w${i}/w${i}_ptensor*
	#rm ./w${i}/${i}pc_ptensor.top
	#cp s_w${i}_ptensor* sheldon/
	#rm s_w${i}_ptensor*
	
	#cp w${i}/s_w${i}_ptensor*  w${i}/sheldon/
	#rm w${i}/s_w${i}_ptensor*
	
	#cp w${i}/w${i}_ptensor*  w${i}/sheldon/
	#cp /Users/burbol/Downloads/small_sams2/w${i}/sheldon/w${i}_ptensor* ./w${i}/sheldon/
	
	#cp /Users/burbol/Downloads/small_sams2/w${i}/soroban/s_w${i}_ptensor* .
	#$cp /Users/burbol/Downloads/small_sams2/w${i}/sheldon/s_w${i}_ptensor* .
	#rm w${i}/w${i}_ptensor*
	
	# copy .itp files 
	#cp *.itp ./w${i}/
	#cp *.itp /Users/burbol/MEGAsync/scripts/Python/SCRIPT_CREATION/sheldon/ptensor/w${i}
	#cp /Users/burbol/Downloads/small_sams2/w${i}/*.itp /Users/burbol/Downloads/small_sams2/WaterSamSeparated/w${i}/
	cp /Users/burbol/Downloads/small_sams2/WaterSamSeparated/*.mdp /Users/burbol/Downloads/small_sams2/WaterSamSeparated/w${i}/
	#rm /Users/burbol/Downloads/small_sams2/WaterSamSeparated/w${i}/\#*
	#Mini.mdp		
	#NPT_PR1.mdp	
	#NPT_PR2.mdp
	
	# copy .mdp files 
	#cp *.mdp ./w${i}/
	#cp *.mdp /Users/burbol/MEGAsync/scripts/Python/SCRIPT_CREATION/sheldon/ptensor/w${i}
	
	#copy scripts
	#cp /Users/burbol/MEGAsync/scripts/Python/SCRIPT_CREATION/sheldon/ptensor/w${i}/s_w${i}_ptensor* ./w${i}/
	#cp /Users/burbol/MEGAsync/scripts/Python/SCRIPT_CREATION/soroban/ptensor/w${i}/s_w${i}_ptensor* ./w${i}/

done
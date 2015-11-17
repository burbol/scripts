#!/bin/bash

#cd /Users/burbol/MEGAsync/scripts/Python/SCRIPT_CREATION/NewVersion1
cd /Users/burbol/Downloads/small_sams2

#for i in 0 5 11 17 21 25 33 50 66 # OH density of the SAM 
for i in 0
do
  #for j in 216 1000 2000 3000 4000 5000 6500 7000 8000 9000 10000  # number of water molecules
  #for j in 216 2000 3000 3000 10000  # number of water molecules s21
  #for j in 216 1000 4000  # number of water molecules s25
  #for j in 216 2000 4000  # number of water molecules s33
  #do 
  
  	#HLRN
  	#scp /Users/burbol/MEGAsync/scripts/SAM_CREATION/SAMs/NEW/drop_placement/*.itp beclaila@hlogin2.hlrn.de:/gfs1/work/beclaila/s${i}_w${j}/
	#scp /Users/burbol/MEGAsync/scripts/Python/SCRIPT_CREATION/s${i}_w${j}/s${i}_w${j}_119 beclaila@hlogin2.hlrn.de:/gfs1/work/beclaila/s${i}_w${j}/
	#scp /Users/burbol/MEGAsync/scripts/Python/SCRIPT_CREATION/s${i}_w${j}/s${i}_w${j}_159 beclaila@hlogin2.hlrn.de:/gfs1/work/beclaila/s${i}_w${j}/
	#scp /Users/burbol/MEGAsync/scripts/Python/SCRIPT_CREATION/s${i}_w${j}/s${i}_w${j}_199 beclaila@hlogin2.hlrn.de:/gfs1/work/beclaila/s${i}_w${j}/
	#scp /Users/burbol/MEGAsync/scripts/Python/SCRIPT_CREATION/s${i}_w${j}/s${i}_w${j}_239 beclaila@hlogin2.hlrn.de:/gfs1/work/beclaila/s${i}_w${j}/

	#SOROBAN	
	#scp -r /Users/burbol/MEGAsync/scripts/Python/SCRIPT_CREATION/soroban/NEWVersionBIG/s${i}_w${j} eixeres@soroban.zedat.fu-berlin.de:/scratch/eixeres/NewVersionBIG/
	#scp /Users/burbol/MEGAsync/scripts/Python/SCRIPT_CREATION/soroban/NEWVersionBIG/s${i}_w${j}/s${i}_w${j}* eixeres@soroban.zedat.fu-berlin.de:/scratch/eixeres/s${i}_w${j}/
	#cp /Users/burbol/MEGAsync/scripts/Python/SCRIPT_CREATION/soroban/NEWVersionBIG/s${i}_w${j}/s${i}_w${j}* /Volumes/BACKUPS/MyTimeMachine/HLRN/s${i}_w${j}/
	#scp /Users/burbol/MEGAsync/scripts/Python/SCRIPT_CREATION/soroban/s${i}_w${j}/s${i}_w${j}*  eixeres@soroban.zedat.fu-berlin.de:/scratch/eixeres/s${i}_w${j}/
	#scp /Users/burbol/MEGAsync/scripts/Python/SCRIPT_CREATION/soroban/NEWVersionBIG/s${i}_w${j}/s${i}_w${j}*  eixeres@soroban.zedat.fu-berlin.de:/scratch/eixeres/NewVersionBIG/s${i}_w${j}/
	#scp /Users/burbol/MEGAsync/scripts/Python/SCRIPT_CREATION/soroban/s${i}_w${j}/s${i}_w${j}*  eixeres@soroban.zedat.fu-berlin.de:/scratch/eixeres/s${i}_w${j}/
  	#scp burbol@burbol.ddns.net:/Users/burbol/MEGAsync/scripts/Python/SCRIPT_CREATION/soroban/s${i}_w${j}/s${i}_w${j}*  s${i}_w${j}/*
  	#scp /Users/burbol/MEGAsync/scripts/Python/SCRIPT_CREATION/s${i}_w${j}/*.itp  eixeres@soroban.zedat.fu-berlin.de:/scratch/eixeres/s${i}_w${j}/

 #./w${i}/sam${i}_water_ptensor.top  eixeres@soroban.zedat.fu-berlin.de:/scratch/eixeres/ptensor/w${i}/
 #./w${i}/sam${i}_water_ptensor.gro  eixeres@soroban.zedat.fu-berlin.de:/scratch/eixeres/ptensor/w${i}/
  scp /Users/burbol/Downloads/small_sams2/w${i}/soroban/s_w${i}_ptensor* eixeres@soroban.zedat.fu-berlin.de:/scratch/eixeres/w${i}/
 
 #./w${i}/sam${i}_water_ptensor.top eixeres@sheldon.physik.fu-berlin.de:/home/eixeres/files_for_laila/ptensor/w${i}
 #./w${i}/sam${i}_water_ptensor.gro eixeres@sheldon.physik.fu-berlin.de:/home/eixeres/files_for_laila/ptensor/w${i}
  scp /Users/burbol/Downloads/small_sams2/w${i}/sheldon/s_w${i}_ptensor* eixeres@sheldon.physik.fu-berlin.de:/home/eixeres/files_for_laila/ptensor/w${i}
  	
  	echo "copied s_w${i}_ptensor"

done
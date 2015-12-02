#!/bin/bash

cd /Users/burbol/Desktop/scripts/Python/SCRIPT_CREATION/NewVersion1

#{216:60, 1000:80, 2000:80, 3000:80, 4000:80, 5000:100, 6500:100, 7000:100, 8000:120, 9000:120, 10000:120}

for i in 33 # OH density of the SAM 
do
  #for j in 216 1000 2000 3000 4000 5000 6500 7000 8000 9000 10000  # number of water molecules
  for j in 216 2000 4000  # number of water molecules
  do 
  
  	#scp /Users/burbol/Desktop/scripts/SAM_CREATION/SAMs/NEW/drop_placement/*.itp beclaila@hlogin2.hlrn.de:/gfs1/work/beclaila/s${i}_w${j}/
  	scp /Users/burbol/Desktop/scripts/Python/SCRIPT_CREATION/soroban/s${i}_w${j}/s${i}_w${j}*  eixeres@soroban.zedat.fu-berlin.de:/scratch/eixeres/s${i}_w${j}/
	#scp /Users/burbol/Desktop/scripts/Python/SCRIPT_CREATION/s${i}_w${j}/s${i}_w${j}_119 beclaila@hlogin2.hlrn.de:/gfs1/work/beclaila/s${i}_w${j}/
	#scp /Users/burbol/Desktop/scripts/Python/SCRIPT_CREATION/s${i}_w${j}/s${i}_w${j}_159 beclaila@hlogin2.hlrn.de:/gfs1/work/beclaila/s${i}_w${j}/
	#scp /Users/burbol/Desktop/scripts/Python/SCRIPT_CREATION/s${i}_w${j}/s${i}_w${j}_199 beclaila@hlogin2.hlrn.de:/gfs1/work/beclaila/s${i}_w${j}/
	#scp /Users/burbol/Desktop/scripts/Python/SCRIPT_CREATION/s${i}_w${j}/s${i}_w${j}_239 beclaila@hlogin2.hlrn.de:/gfs1/work/beclaila/s${i}_w${j}/
  
  done
done

# ojo
#  s25_w4000
#  s25_w1000
#  s25_w216

# s33_w4000
#  s33_w216
#  s33_w2000

# s21_w216
# s21_w2000
# s21_w3000
# s21_w10000
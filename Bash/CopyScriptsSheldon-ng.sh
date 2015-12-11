#!/bin/bash

cd /Users/burbol/MEGAsync/scripts/Python/SCRIPT_CREATION/sheldon-ng/

#for i in 33 50 # OH density of the SAM 
for i in 25 
do
 for j in 1000 2000 3000 4000 5000 6500 7000 8000 9000 10000
 #for j in 216 3000 4000 5000 6500 7000 8000 9000 10000
 #for j in 7000 8000 9000 10000


 do

  cd s${i}_w${j}
  #scp s${i}_w${j}_* eixeres@sheldon-ng.physik.fu-berlin.de:/scratch/eixeres/NewVersionBIG/s${i}_w${j}/
  scp *s${i}_w${j}_new_* eixeres@sheldon-ng.physik.fu-berlin.de:/scratch/eixeres/Dec14_Last_Sims/s${i}_w${j}/

  cd ..
  done
done
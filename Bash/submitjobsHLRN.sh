#!/bin/bash

for i in 33 # OH density of the SAM 
#for i in 21 # For testing!!!!
do
  for j in 216 1000 2000 3000 4000 5000 6500 7000 8000 9000 10000  # number of water molecules
  do 
  
  cd s${i}_w${j}
  msub s${i}_w${j}_0
  cd ..
  
  done
done
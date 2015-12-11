#!/bin/bash

for i in 21 25 # OH density of the SAM 
do
  for j in 216 1000 2000 3000 4000 5000 6500 7000 8000 9000 10000  # number of water molecules
  do 
	
	mkdir /scratch/eixeres/HLRN/s${i}_w${j}
  
  done
done
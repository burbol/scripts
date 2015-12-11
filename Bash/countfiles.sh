#!/bin/bash

# Script for counting number of output files (there should be 9 for each simulation)


for i in 0 5 11 17 # OH density of the SAM 
do
  for j in 1000 2000 3000 4000 5000 6500 7000 8000 9000 10000 # number of water molecules
  do
      
        echo -e "Number of  output files for simulation NVT_sam${i}_water${j}_next  \c " && ls NVT_sam${i}_water${j}_next.* | wc -l
	#echo -e "A \c " && echo "B"
	# printf "%s %s\n" "String 1" "String 2"

  done
done


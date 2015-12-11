#!/bin/bash

source /usr/local/gromacs/bin/GMXRC

#script to move droplets on top of sams -> they still need to be copied & .top files changed

cd /Users/burbol/Desktop/scripts/SAM_CREATION/SAMs/NEW/drop_placement

#y-sizes
#sam21=$(echo "12.725" | bc)
#sam21_half=$(echo "sam21/2" | bc)
sam21="12"".""725"
for i in 21
do
	samname='sam'$i
	#echo "value $sam21"
	echo "name: $samname" #->**
	#echo "value through name: $(($samname))"
	#value=$(echo "$((sam$i))" | bc)
	echo "$((sam$i))"
	#echo "value : $value"

done
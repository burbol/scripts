#!/bin/bash

# This script changes the names of the density profiles produced with g_rad_grompi_all.sh
# The value of shift in that script should have been zero. Thus, all the numbers in the 
# names of the produced xvg files are wrong (shifted, in this case by 42.5ns).
# This script corrects those file names . For example the file 
# g_rad_dmap_5pc_w2000_42.5ns_43.0ns.xvg  will be renamed to g_rad_dmap_5pc_w2000_0ns_0.5ns.xvg

# the round function:
round()
{
echo $(printf %.$2f $(echo "scale=$2;(((10^$2)*$1)+0.5)/(10^$2)" | bc))
};



for i in 5 21 25 41   # OH density of the SAM 
do
  #for j in 2000
  for j in 3000 4000 5000 6500 7000 8000 9000  # number of water molecules
  do
    cd /net/data/eixeres/NewVersion4/FINISHED/s${i}_w${j}
    
    k1=0
    shift=425
    while [[ $k1 -le 995 ]]  # segments of time to be multiplied by 100 so that we get [ps]
      do    
        shift=425
        k2=$((k1+5))
        nanosecs1=$(echo $(round $k1/10+$shift 1)) # the sum is computed before the division!! (why?)
        nanosecs2=$(echo $(round $k2/10+$shift 1)) # the sum is computed before the division!!
        start=$((k1*100))
        ending=$((start+500))

        shift=0
        nanosecs1a=$(echo $(round $k1/10+$shift 1)) # the sum is computed before the division!! (why?)
        nanosecs2a=$(echo $(round $k2/10+$shift 1)) # the sum is computed before the division!!

        FILE=g_rad_dmap_${i}pc_w${j}_${nanosecs1}ns_${nanosecs2}ns.xvg
        NEWFILE=g_rad_dmap_${i}pc_w${j}_${nanosecs1a}ns_${nanosecs2a}ns.xvg
        #echo ${FILE}
        #echo ${NEWFILE}
        mv $FILE $NEWFILE ;
        
        #output to see what the script is doing
        #echo "mv" $FILE $NEWFILE ;
        echo $i $j
        
        k1=$(($k1+5))
      done



done
  done
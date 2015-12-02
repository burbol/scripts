#!/bin/bash

# Script for using recursively g_rad_density in time intervals of 500 ps from 20ns to 40ns to produce radial (cilindrical) density maps

#source /usr/local/gromacs/bin/GMXRC


# the round function:
round()
{
echo $(printf %.$2f $(echo "scale=$2;(((10^$2)*$1)+0.5)/(10^$2)" | bc))
};


n=4
cd /home/eixeres/files_for_laila
#mkdir ./global_density_maps
mkdir ./g_rad_densmaps_s11_repeated

#FIRST ROUND FOR BIG DROPLETS
#for i in 0 5 17 66 # OH density of the SAM 
#do
  #for j in 5000 6500 7000 8000 9000 10000 # number of water molecules
  #do 
  
  #if [ $i -eq 0 ]
  #then 
    #n=3
  #else 
    #n=4
  #fi

    #k=200
    #while [[ $k -le 495 ]]  # segments of time to be multiplied by 100 so that we get [ps]
    #do    
        #k2=$((k+5))
        #nanosecs1=$(echo $(round $k/10 1))
        #nanosecs2=$(echo $(round $k2/10 1))
        #start=$((k*100))
        #ending=$((start+500))

        #echo "Creating densmap from $start ps to $ending ps"
        
        #echo ${n} ${n} | /home/eixeres/g_rad_density -f NVT_sam${i}_water${j}.xtc -s g_rad_NVT_sam${i}_water${j}.tpr -n index${i}_${j}.ndx -sz 200 -o g_rad_dmap_${i}pc_w${j}_${nanosecs1}ns_${nanosecs2}ns.xvg -b ${start} -e ${ending}

        
        #cp g_rad_dmap_${i}pc_w${j}_${nanosecs1}ns_${nanosecs2}ns.xvg ./g_rad_densmaps_p2/
        #rm g_rad_dmap_${i}pc_w${j}_${nanosecs1}ns_${nanosecs2}ns.xvg
	
	    #k=$(($k+5))
	#done        
  #done
#done

#SECOND ROUND FOR SAM 11%
#for i in 11 # OH density of the SAM 
#do
  #for j in 5000 7000 8000 10000 # number of water molecules
  #do
  
  #if [ $i -eq 0 ]
  #then 
    #n=3
  #else 
    #n=4
  #fi

    #k=200
    #while [[ $k -le 495 ]]  # segments of time to be multiplied by 100 so that we get [ps]
    #do    
        #k2=$((k+5))
        #nanosecs1=$(echo $(round $k/10 1))
        #nanosecs2=$(echo $(round $k2/10 1))
        #start=$((k*100))
        #ending=$((start+500))

        #echo "Creating densmap from $start ps to $ending ps"
        
        #echo ${n} ${n} | /home/eixeres/g_rad_density -f NVT_sam${i}_water${j}.xtc -s g_rad_NVT_sam${i}_water${j}.tpr -n index${i}_${j}.ndx -sz 200 -o g_rad_dmap_${i}pc_w${j}_${nanosecs1}ns_${nanosecs2}ns.xvg -b ${start} -e ${ending}

        
        #cp g_rad_dmap_${i}pc_w${j}_${nanosecs1}ns_${nanosecs2}ns.xvg ./g_rad_densmaps_p2/
        #rm g_rad_dmap_${i}pc_w${j}_${nanosecs1}ns_${nanosecs2}ns.xvg
	
	    #k=$(($k+5))
	#done        
  #done
#done

#THIRD ROUND FOR THE SMALL DROPLETS
for i in 11 # OH density of the SAM 
do
  for j in 6500 9000  # number of water molecules
  do
  
  /home/eixeres/gromacs_tpi_compiled/bin/grompp -f NVT.mdp -c NVT1_sam${i}_water${j}.gro -p ${i}pc1_${j}.top -n index${i}_${j}.ndx -o g_rad_NVT_sam${i}_water${j}.tpr -maxwarn 1
   
  
  if [ $i -eq 0 ]
  then 
    n=3
  else 
    n=4
  fi

    k=0
    while [[ $k -le 395 ]]  # segments of time to be multiplied by 100 so that we get [ps]
    do    
        k2=$((k+5))
        nanosecs1=$(echo $(round $k/10 1))
        nanosecs2=$(echo $(round $k2/10 1))
        start=$((k*100))
        ending=$((start+500))

        echo "Creating densmap from $start ps to $ending ps"
        
        echo ${n} ${n} | /home/eixeres/g_rad_density -f NVT1_sam${i}_water${j}.xtc -s g_rad_NVT_sam${i}_water${j}.tpr -n index${i}_${j}.ndx -sz 200 -o g_rad_dmap_${i}pc_w${j}_${nanosecs1}ns_${nanosecs2}ns.xvg -b ${start} -e ${ending}

        
        cp g_rad_dmap_${i}pc_w${j}_${nanosecs1}ns_${nanosecs2}ns.xvg ./g_rad_densmaps_s11_repeated/
        rm g_rad_dmap_${i}pc_w${j}_${nanosecs1}ns_${nanosecs2}ns.xvg
	
	    k=$(($k+5))
	    
	    cp g_rad_dmap_${i}pc_w${j}_${nanosecs1}ns_${nanosecs2}ns.xvg ./g_rad_densmaps_s11_repeated/
        rm g_rad_dmap_${i}pc_w${j}_${nanosecs1}ns_${nanosecs2}ns.xvg
        rm \#mdout*
	done        
  done
done
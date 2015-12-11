#!/bin/bash

# Script for using g_density recursively and g_rad_density in time intervals of 500 ps
# Run script from slurm server

#source /usr/local/gromacs/bin/GMXRC
#######source /Volumes/UNI/SHELDON/CLOSEDsheldon/eixeres/gromacs_tpi_compiled/bin/GMXRC
module load gromacs/single/openmpi1.4.5/4.6.5


# the round function:
round()
{
echo $(printf %.$2f $(echo "scale=$2;(((10^$2)*$1)+0.5)/(10^$2)" | bc))
};


n=4
#cd  /media/SHELDON/files_for_laila/
#cd /media/SHELDON/CLOSEDsheldon/eixeres/files_for_laila/s11_w2000_merged
#cd  /Volumes/UNI/SHELDON/files_for_laila/s0_w1000/
#cd /net/clusterhome/eixeres/files_for_laila/Simulations_to_delete/s0_w1000
#cd /Volumes/Backup/YosemiteFiles/Simulations
#cd /media/SHELDON/Last_Simulations_Dec14/s33_w9000
#mkdir ./global_density_maps
#mkdir ./g_rad_densmaps_next

for i in 41 # OH density of the SAM 
#for i in 5 21 25 41 # OH density of the SAM 
do
  for j in 8000
  #for j in 3000 4000 5000 6500 7000 8000 9000
  #for j in 1000 2000 3000 4000 5000 6500 7000 8000 9000 10000  # number of water molecules
  do
  
  #cd /home/eixeres/files_for_laila/s${i}_w${j}_merged
  #cd /media/SHELDON/CLOSEDsheldon/eixeres/files_for_laila/s${i}_w${j}_merged
  cd /net/data/eixeres/NewVersion4/FINISHED/s${i}_w${j}
  #cd s${i}_w${j}
  
  if [ $i -eq 0 ]
  then 
    n=3
  else 
    n=4
  fi
  

#######################################################################
########## THESE LINES WERE ADDED TO CREATE THE .gro FILES OF SIMULATIONS NOT ENDED NORMALLY ##########
    temp=0
    if [ $j -eq 1000 ]
    then temp=19000
    elif [ $j -eq 2000 ]
    then temp=43000
    elif [ $j -eq 3000 ]
    then temp=51000
    elif [ $j -eq 4000 ]
    then temp=47000
    elif [ $j -eq 5000 ]
    then temp=48000
    elif [ $j -eq 6500 ]
    then temp=44000
    elif [ $j -eq 7000 ]
    then temp=93000
    elif [ $j -eq 8000 ]
    then temp=63000
    elif [ $j -eq 9000 ]
    then temp=19000
    elif [ $j -eq 10000 ]
    then temp=26000
    fi
    
    #echo "0"|/home/shavkat/GMX/bin/trjconv -f NVT_sam${i}_water${j}.xtc # watch output for last frame determination
    #/home/shavkat/GMX/bin/grompp -f NVT.mdp -p ${i}pc_${j}.top -t NVT_sam${i}_water${j}.trr -c Mini_sam${i}_water${j}.gro -o NVT_sam${i}_water${j}_temp.tpr -maxwarn 2
    #echo "0"|/home/shavkat/GMX/bin/trjconv -f NVT_sam${i}_water${j}.xtc  -s NVT_sam${i}_water${j}_temp.tpr -o NVT_sam${i}_water${j}_temp.gro -b $temp -sep
    #Afterwards look for the last .gro file, rename and delete the rest
    
########## INSERTION ENDED  ##########
#######################################################################
    
    echo "0" "q"|make_ndx -f sam${i}_water${j}.gro -o index${i}_${j}.ndx 
 
    /net/data/eixeres/sheldon-old/gromacs_tpi_compiled/bin/grompp -f NVT.mdp -c NVT_sam${i}_water${j}.gro -p ${i}pc_${j}.top -n index${i}_${j}.ndx -o g_rad_NVT_sam${i}_water${j}.tpr -maxwarn 1
    
    echo ${n} | g_density -f NVT_sam${i}_water${j}.xtc -s g_rad_NVT_sam${i}_water${j}.tpr -o g_density_NVT_sam${i}_water${j}.xvg -sl 1000
    echo "6" | g_density -f NVT_sam${i}_water${j}.xtc -s g_rad_NVT_sam${i}_water${j}.tpr -o g_density_SAM_sam${i}_water${j}.xvg -sl 1000 # for SAMs density maps
    

    # atom number global density maps
    
    #grompp -f NVT.mdp -c NVT_sam${i}_water${j}.gro -p ${i}pc_${j}.top -n index${i}_${j}.ndx -o g_rad_NVT_sam${i}_water${j}.tpr -maxwarn 1
    echo ${n} | g_density -dens number -f NVT_sam${i}_water${j}.xtc -s g_rad_NVT_sam${i}_water${j}.tpr -o ng_density_NVT_sam${i}_water${j}.xvg -sl 1000
    #cp    ng_density_NVT_sam${i}_water${j}.xvg /home/eixeres/Dropbox/Apps/Computable/December/StepByStepDropletMethod/
    
    echo "5" |g_density -dens number -f NVT_sam${i}_water${j}.xtc -s g_rad_NVT_sam${i}_water${j}.tpr -o ng_density_SAM_sam${i}_water${j}.xvg -sl 1000
    #cp ng_density_SAM_sam${i}_water${j}.xvg /home/eixeres/Dropbox/Apps/Computable/December/StepByStepDropletMethod/

#cp g_density_SAM_sam${i}_water${j}.xvg /home/eixeres/files_for_laila/global_SAMS_densmaps/
#rm g_density_NVT_sam${i}_water${j}.xvg	

    #k=200
    k=0
    #while [[ $k -le 395 ]]  # segments of time to be multiplied by 100 so that we get [ps]
    while [[ $k -le 995 ]]  # segments of time to be multiplied by 100 so that we get [ps]
    do    
        k2=$((k+5))
        nanosecs1=$(echo $(round $k/10 1))
        nanosecs2=$(echo $(round $k2/10 1))
        start=$((k*100))
        ending=$((start+500))

        #echo "Creating densmap from $start ps to $ending ps"

        echo ${n} ${n} | /net/data/eixeres/sheldon-old/g_rad_density -f NVT_sam${i}_water${j}.xtc -s g_rad_NVT_sam${i}_water${j}.tpr -n index${i}_${j}.ndx -sz 200 -o g_rad_dmap_${i}pc_w${j}_${nanosecs1}ns_${nanosecs2}ns.xvg -b ${start} -e ${ending}
        #echo ${n} ${n} | /home/eixeres/g_rad_density -f NVT_sam${i}_water${j}_next.xtc -s g_rad_NVT_sam${i}_water${j}_next.tpr -n index${i}_${j}.ndx -sz 200 -o g_rad_dmap_${i}pc_w${j}_${nanosecs1}ns_${nanosecs2}ns.xvg -b ${start} -e ${ending}
        #echo ${n} ${n} | /home/eixeres/g_rad_density -f NVT_sam${i}_water${j}_40ns.xtc -s g_rad_NVT_sam${i}_water${j}_40ns.tpr -n index${i}_${j}.ndx -sz 200 -o g_rad_dmap_${i}pc_w${j}_${nanosecs1}ns_${nanosecs2}ns.xvg -b ${start} -e ${ending}

        
        #cp g_rad_dmap_${i}pc_w${j}_${nanosecs1}ns_${nanosecs2}ns.xvg /home/eixeres/files_for_laila/g_rad_densmaps_all/
        
        #cp g_rad_dmap_${i}pc_w${j}_${nanosecs1}ns_${nanosecs2}ns.xvg /home/eixeres/files_for_laila/g_rad_densmaps_next/
        #rm g_rad_dmap_${i}pc_w${j}_${nanosecs1}ns_${nanosecs2}ns.xvg
	
	    k=$(($k+5))
	done    
  done
done

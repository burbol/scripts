#!/bin/bash

# Script for equilibrating water boxes with only water moelcules

source /usr/local/gromacs/bin/GMXRC
cd /Users/burbol/Downloads/small_sams2/WaterSamSeparated

#for i in 0 5 11 17 21 25 33 50 66 # OH density of the SAM 
for i in 0 # For testing
do
        
		cd w${i}

        echo "grompp -f minim.mdp -c sam${i}_water_ptensor.gro -p sam${i}_water_ptensor.top -o Mini_sam${i}_water_ptensor.tpr -maxwarn 1"
        #grompp -f minim.mdp -c sam${i}_water_ptensor.gro -p sam${i}_water_ptensor.top -o Mini_sam${i}_water_ptensor.tpr -maxwarn 1
        echo "mdrun -testverlet -deffnm Mini_sam${i}_water_ptensor"
        #mdrun -deffnm Mini_sam${i}_water_ptensor
        echo "grompp -f nvt.mdp -c Mini_sam${i}_water_ptensor.gro -p sam${i}_water_ptensor.top -o NVT_sam${i}_water_ptensor.tpr -maxwarn 1"
        #grompp -f nvt.mdp -c Mini_sam${i}_water_ptensor.gro -p sam${i}_water_ptensor.top -o NVT_sam${i}_water_ptensor.tpr -maxwarn 1
        echo "mdrun -testverlet -deffnm NVT_sam${i}_water_ptensor"
        #mdrun -deffnm NVT_sam${i}_water_ptensor
        echo "grompp -f npt.mdp -c NVT_sam${i}_water_ptensor.gro -p sam${i}_water_ptensor.top -o NPT_sam${i}_water_ptensor.tpr -maxwarn 1"
        #grompp -f npt.mdp -c NVT_sam${i}_water_ptensor.gro -p sam${i}_water_ptensor.top -o NPT_sam${i}_water_ptensor.tpr -maxwarn 1
        mdrun -deffnm NPT_sam${i}_water_ptensor
        
        
        #echo "grompp -f npt.mdp -c NVT_sam${i}_water_ptensor.gro -p sam${i}_water_ptensor.top -o NPT_PR1_sam${i}_water_ptensor.tpr -maxwarn 1"
        #grompp -f npt.mdp -c NVT_sam${i}_water_ptensor.gro -p sam${i}_water_ptensor.top -o NPT_PR1_sam${i}_water_ptensor.tpr -maxwarn 1
        #echo "mdrun -testverlet -deffnm NPT_PR1_sam${i}_water_ptensor"
        #mdrun -testverlet -deffnm NPT_PR1_sam${i}_water_ptensor
        #echo "grompp -f NPT_PR1.mdp -c NVT_sam${i}_water_ptensor.gro -p sam${i}_water_ptensor.top -o NPT_PR1_sam${i}_water_ptensor.tpr -maxwarn 1"
        #grompp -f NPT_PR1.mdp -c NVT_sam${i}_water_ptensor.gro -p sam${i}_water_ptensor.top -o NPT_PR1_sam${i}_water_ptensor.tpr -maxwarn 1
        #echo "mdrun -testverlet -deffnm NPT_PR1_sam${i}_water_ptensor"
        #echo "grompp -f NPT_PR2.mdp -c NPT_PR1_sam${i}_water_ptensor.gro -p sam${i}_water_ptensor.top -o NPT_PR2_sam${i}_water_ptensor.tpr -maxwarn 1"
        #grompp -f NPT_PR2.mdp -c NPT_PR1_sam${i}_water_ptensor.gro -p sam${i}_water_ptensor.top -o NPT_PR2_sam${i}_water_ptensor.tpr -maxwarn 1
        #mdrun -testverlet -deffnm NPT_PR2_sam${i}_water_ptensor"
        #grompp -f NPT.mdp -c NPT_PR2_sam${i}_water_ptensor.gro -p sam${i}_water_ptensor.top -o NPT_sam${i}_water_ptensor.tpr -maxwarn 1"
        #mdrun -testverlet -deffnm NPT_sam${i}_water_ptensor"
        cd ..
done
#rm \#*

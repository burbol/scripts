#### THIS SCRIPT DOESN'T WORK, because arrays can't handle floating point numbers
####!!! USE placedrop.py instad. It's a python implementation of this program

#!/bin/bash

source /usr/local/gromacs/bin/GMXRC

#script to move droplets on top of sams -> they still need to be copied & .top files changed

cd /Users/burbol/Desktop/scripts/SAM_CREATION/SAMs/NEW/drop_placement

#y-sizes divided by 2
sam[21]=$(echo "12.725" | bc)
sam[25]=$(echo "12.725" | bc)
sam[33]=$(echo "13.585" | bc)

xbox[21]=$(echo "10.000" | bc)
ybox[21]=$(echo "13.191" | bc)
zbox[21]=$(echo "12.000" | bc)

xbox[25]=$(echo "10.000" | bc)
ybox[25]=$(echo "13.191" | bc)
zbox[25]=$(echo "12.000" | bc)

xbox[33]=$(echo "13.500" | bc)
ybox[33]=$(echo "13.856" | bc)
zbox[33]=$(echo "12.000" | bc)

water[1000]=$(echo "3.325" | bc)
water[2000]=$(echo "4.122" | bc)
water[3000]=$(echo "4.635" | bc)
water[4000]=$(echo "5.123" | bc)
water[5000]=$(echo "5.501" | bc)
water[6500]=$(echo "6.005" | bc)
water[7000]=$(echo "6.121" | bc)
water[8000]=$(echo "6.395" | bc)
water[9000]=$(echo "6.660" | bc)
water[10000]=$(echo "6.902" | bc)


for j in 1000
#for j in 1000 2000 3000 4000 5000 6500 7000 8000 9000 10000
do
	for i in 21 
	#for i in 21 25 33
	do

    #file names
	#startwaterfile='NPT_water'$j'.gro'
	#startsamfile='start'$i'.gro'
	waterfile='NPT_water'$j'_c.gro'
	samfile='start'$i'_c.gro'
	
	newsystfile='NPT_sam'$i'_water'$j'_c.gro'
	
	dist=$(echo "$waterhalf+$samhalf" | bc)
	
	#we change the water box size too the box size of the sams
	#editconf -f $startwaterfile -o $waterfile -box $((xbox[$i])) $((ybox[$i])) $((zbox[$i]))
	
	#we center every system (water and sam!)
	#editconf -f $samfile -o $samfile -c
	#editconf -f $waterfile  -o $waterfile -c
	
	#we move the drops up and the sams down
	#editconf -f $samfile -o $samfile -translate 0 $samhalf 0 
	#editconf -f $waterfile -o $waterfile -translate 0 $waterhalf 0 
	
	##### TEST #######
	echo "water= $water"
	echo "sam= $sam"
	echo "waterhalf = $waterhalf"
	echo "samhalf = $samhalf"
	echo "-f $water -o $water -translate 0 $waterhalf 0 "
	echo "editconf -f $sam -o $sam -translate 0 $samhalf 0 "
	
	#we copy the drops .gro file at the end to the sam's .gro file (the rest has to be done by hand)
	#cat $sam>>$newsyst
	#cat $water>>$newsyst
	done
done


#####  echo $((grofile[$i]))  #####
#####  echo $((water$i)) #####
#####  ydist=$(echo "$6" | bc)  #####
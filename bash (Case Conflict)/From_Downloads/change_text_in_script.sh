#!/bin/bash

#script to change text inside numbered script to send to SLURM

for j in 17
do
   #for i in {0..79}
   #do
   
    i=4
    cd /Users/burbol/ownCloud/scripts/sheldon-ng/s17_w1000/
    job='s'$j'_w1000_'$i
    sed -i -e 's/nodes\=3/nodes\=1/g' $job
    #sed -i -e 's/nodes\=1/nodes\=1\ \n\ \#SBATCH\ --exclusive/' $job
    
    #To repair:
    #sed -i -e 's/\#SBATCH\ --exclusive/\ /' $job
    
    #OLD
	#job='sam'$j'_water216_continue'
    #sed -i -e 's/job = \'sam17_water216_continue\' + str\(i\)/job='sam17_water216_continue'$i/g' $job
    #sed -i -e 's/ + str(i)/$i/g' $job
   #done
done
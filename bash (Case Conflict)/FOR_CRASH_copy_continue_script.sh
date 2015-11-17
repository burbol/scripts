#!/bin/bash

#script to create N numbered copies of a script to send to SLURM

for i in 1 2 3 4 5 6 7 8 9
do
j=$((i+1))
job='sam0_water216_continue'$i
newjob='sam17_water216_continue'$j
cp $job $newjob
done
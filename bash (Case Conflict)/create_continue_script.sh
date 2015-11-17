#!/bin/bash

#script to create N numbered copies of a script to send to SLURM

for i in 1 2 3 4 5 6 7 8 9
do
job='sam0_water216_continue'$i
newjob='sam0_water216_continue'$i
cp $job $newjob
done
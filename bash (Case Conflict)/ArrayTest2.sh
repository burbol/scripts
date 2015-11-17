#!/bin/bash

i=0;
for j in `cat gro_files.txt` 
do
   grofile[$i]=$j 
   echo $i
   echo $grofile[i]
    i=$(($i+1)) 
   
done 
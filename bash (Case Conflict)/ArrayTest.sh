#!/bin/bash
i=1
read grofile < gro_files.txt
for i in 8 
 do
 
  echo "i is $i"
  echo $grofile[$i]
done

#!/bin/bash


cd /Users/burbol/Desktop/scripts/SAM_CREATION/SAMs/NEW/drop_placement
grep -c "SAM" NPT_sam21_water2000.gro
grep -c "OAM" NPT_sam21_water2000.gro

#echo "SAM num= $i, OAM num= $j
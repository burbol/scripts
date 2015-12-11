#!/bin/bash

cd /Users/burbol/Downloads/SAMs

for i in 0 5 11 17 # OH density of the SAM 
do

awk '{print $6}' C19HeadGroups${i}.dat > C19xCoordsam${i}.txt

awk '{print $7}' C19HeadGroups${i}.dat > C19yCoordsam${i}.txt

awk '{print $8}' C19HeadGroups${i}.dat > C19zCoord${i}.txt

awk '{print $6}' O1HeadGroups${i}.dat > O1xCoord${i}.txt

awk '{print $7}' O1HeadGroups${i}.dat > O1yCoord${i}.txt

awk '{print $8}' O1HeadGroups${i}.dat > O1zCoord${i}.txt

done
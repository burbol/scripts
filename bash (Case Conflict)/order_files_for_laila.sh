#!/bin/bash

#cd /Volumes/UNI/SHELDON/files_for_laila/
cd /Volumes/UNI/SHELDON/files_for_laila/Simulations_to_delete/

cd /net/clusterhome/eixeres/files_for_laila/Simulations_to_delete

for i in 0 5 11 17 66 # OH density of the SAM 
#for i in 0
do
  for j in 10000 1000 2000 3000 4000 5000 6500 7000 8000 9000  # number of water molecules
 #for j in 1000
  do
     mkdir s${i}_w${j}
     cp *sam${i}_water${j}* s${i}_w${j}/ 
     cp *_${j}_${i} s${i}_w${j}/
     cp ${i}pc_${j}.top s${i}_w${j}/
     cp ${i}pc1_${j}.top s${i}_w${j}/
    
      cp index${i}_${j}.ndx  s${i}_w${j}/
      cp *_s${i}_w${j} s${i}_w${j}/
      cp *_s${i}_w${j}_next s${i}_w${j}/

    # echo "ls files_for_laila/*sam${i}_water${j}*" 
    # echo "ls files_for_laila/s${i}_w${j}/" 
    # echo "ls files_for_laila/*_${j}_${i} ls "
    # echo "ls files_for_laila/s${i}_w${j}/"
    # echo "ls files_for_laila/${i}pc_${j}.top "
    # echo "ls files_for_laila/s${i}_w${j}/"
     
    # cp *.itp s${i}_w${j}/
    # cp *.mdp s${i}_w${j}/
    
      rm index${i}_${j}.ndx  
      rm *_s${i}_w${j} 
      rm *_s${i}_w${j}_next 
 
     rm *sam${i}_water${j}*  
     rm *_${j}_${i}
     rm ${i}pc_${j}.top
     rm ${i}pc1_${j}.top
  
  done
done
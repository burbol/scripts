#!/bin/bash

for i in 21 25 # OH density of the SAM 
do
  for j in 1000 2000 3000 4000 5000 6500 7000 8000 9000 10000  # number of water molecules
  do 
  
# First we uncomment the next line to copy the files to BERT
#scp -r /Volumes/Backup/YosemiteFiles/MEGAsync/scripts/SAM_CREATION/SAMs/NEW/drop_placement/NewVersion3/s${i}_w${j} eixeres@bert.physik.fu-berlin.de:/home/eixeres/Downloads/
#scp -r /Volumes/Backup/YosemiteFiles/MEGAsync/scripts/SAM_CREATION/SAMs/NEW/drop_placement/NewVersion3/s${i}_w${j}/*cuda.top eixeres@bert.physik.fu-berlin.de:/home/eixeres/Downloads/s${i}_w${j}/
#scp /Volumes/Backup/YosemiteFiles/MEGAsync/scripts/SAM_CREATION/SAMs/NEW/drop_placement/NewVersion3/s50_w4000/*cuda.mdp eixeres@bert.physik.fu-berlin.de:/home/eixeres/Downloads/s${i}_w${j}/
##scp /Users/burbol2/Dropbox/scripts/Python/SCRIPT_CREATION/HLRN/s${i}_w${j}/*_cuda_* eixeres@bert.physik.fu-berlin.de:/home/eixeres/Downloads/s${i}_w${j}/
#scp /Volumes/Backup/YosemiteFiles/MEGAsync/scripts/SAM_CREATION/SAMs/NEW/drop_placement/NVT_*ns.mdp eixeres@bert.physik.fu-berlin.de:/home/eixeres/Downloads/s${i}_w${j}/
#scp /Volumes/Backup/YosemiteFiles/MEGAsync/scripts/SAM_CREATION/SAMs/NEW/drop_placement/NewVersion3/s${i}_w${j}/${i}pc_${j}.top eixeres@bert.physik.fu-berlin.de:/home/eixeres/Downloads/s${i}_w${j}/


### Then we copy the next line manually to the shell to connect to BERT
###ssh eixeres@bert.physik.fu-berlin.de


# Then we uncomment the next line to copy the folders to the work directory
#scp -r /home/eixeres/Downloads/s${i}_w${j} beclaila@blogin.hlrn.de:/gfs1/work/beclaila/
#scp /home/eixeres/Downloads/s${i}_w${j}/*cuda.mdp beclaila@blogin.hlrn.de:/gfs1/work/beclaila/s${i}_w${j}/
#scp /home/eixeres/Downloads/s${i}_w${j}/*cuda.top beclaila@blogin.hlrn.de:/gfs1/work/beclaila/s${i}_w${j}/
scp /home/eixeres//Dropbox/scripts/Python/SCRIPT_CREATION/HLRN/s${i}_w${j}/s${i}_w${j}_cuda* beclaila@blogin.hlrn.de:/gfs1/work/beclaila/s${i}_w${j}/
#scp /home/eixeres/Downloads/s${i}_w${j}/NVT_*ns.mdp beclaila@blogin.hlrn.de:/gfs1/work/beclaila/s${i}_w${j}/
#scp /home/eixeres/Downloads/s${i}_w${j}/${i}pc_${j}.top beclaila@blogin.hlrn.de:/gfs1/work/beclaila/s${i}_w${j}/


###Then we connect MANUALLY to HLRN (through BERT) 
###ssh beclaila@blogin.hlrn.de

### We copy this script to HLRN and run it from there: 
#we uncomment the next line to copy the missing .itp and .mdp files to each folder
#cp /gfs1/work/beclaila/s21_w1000/*.itp /gfs1/work/beclaila/s${i}_w${j}/
#cp /gfs1/work/beclaila/s21_w1000/*.mdp /gfs1/work/beclaila/s${i}_w${j}/

  done
done


# #test scripts copied manually:
# scp /Users/burbol2/Dropbox/scripts/Python/SCRIPT_CREATION/HLRN/s21_w1000/s21_w1000_testcuda eixeres@bert.physik.fu-berlin.de:/home/eixeres/Downloads/s21_w1000/
# scp /Users/burbol2/Dropbox/scripts/Python/SCRIPT_CREATION/HLRN/s25_w1000/s25_w1000_testcuda eixeres@bert.physik.fu-berlin.de:/home/eixeres/Downloads/s25_w1000/
# scp /home/eixeres/Downloads/s21_w1000/s21_w1000_testcuda beclaila@blogin.hlrn.de:/gfs1/work/beclaila/s21_w1000/
# scp /home/eixeres/Downloads/s25_w1000/s25_w1000_testcuda beclaila@blogin.hlrn.de:/gfs1/work/beclaila/s25_w1000/


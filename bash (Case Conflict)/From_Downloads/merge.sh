
source /usr/local/gromacs/bin/GMXRC

for i in 0 5 11 17 # OH density of the SAM 
do
  for j in 1000 2000 3000 4000  # number of water molecules
  do
  
  echo 0.000 2000.000/home/shavkat/trjcat -f NVT_sam${i}_water${j}.xtc NVT2_sam${i}_water${j}.xtc -n index${i}_${j}.ndx -o NVT_sam${i}_water${j}_MERGED.tpr -settime
  
  rm \#mdout*
        
	done        
  done
done
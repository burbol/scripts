#!/bin/bash

# Different water droplets situated on top of various SAMs have been simulated with GROMACS
# This bash script translates the x- and y-coordinates of the center of each droplet each 2ns
# to exactly the center of the top of the surface.
# 


# To control see the outputk.dat and output2k.dat files with k being the time period in ns
# ATTENTION! After each run type rm ./#average* 

source /usr/local/gromacs/bin/GMXRC

for i in 66 # OH density of the SAM 
do
  for j in 2000  # number of water molecules
  do
  mkdir sam${i}_water${j}
  mkdir densmaps_s${i}_w${j}
    for k in 6 8 10 12 14 16 18   # segments of time in ns
    do
        start=$((k*1000))
        ending=$((start+2000))
        l=$((k+2))

# First we produce a coordinate file ".xvg" of the water droplet
        echo water | g_traj -f NVT_sam${i}_water${j}.xtc -s  NVT_sam${i}_water${j}.tpr -ox  coord${k}_sam${i}_water${j}.xvg -com -b ${start} -e ${ending}

# Now we analize the ".xvg" file to determine the x- and y-coordinates of the average center of mass of the water droplet 
       g_analyze -f coord${k}_sam${i}_water${j}.xvg -av | grep SS1 > av${k}_sam${i}_water${j}
	   rm average.xvg
       g_analyze -f coord${k}_sam${i}_water${j}.xvg -av | grep SS2 >> av${k}_sam${i}_water${j}
	   rm average.xvg
 
# Now we determine the simulation's box size and then we divide the x- and y-coordinates by 2 to get the center of the box, which is also the center of the surface  
        awk 'END {lastline = $NR; print '$lastline'}' NVT_sam${i}_water${j}.gro > boxsize_sam${i}_water${j}.dat
        
        
        oldx=$(awk '/SS1/ {print 1*$2}' av${k}_sam${i}_water${j})

        newx=$(awk 'END {print 0.5*$1}' boxsize_sam${i}_water${j}.dat)
        
        oldy=$(awk '/SS2/ {print 1*$2}' av${k}_sam${i}_water${j})
      
        newy=$(awk 'END {print 0.5*$2}' boxsize_sam${i}_water${j}.dat)
           
        echo "Moving sam${i}_water${j} in the simulation period from $k to $l ns " > output${k}.dat
        
        echo "The okd x is $oldx" >> output${k}.dat
        echo "The new x is $newx" >> output${k}.dat
        echo "The old y is $oldy" >> output${k}.dat
        echo "The new y is $newy" >> output${k}.dat
        
        x=$(echo "$newx - $oldx" | bc)
        y=$(echo "$newy - $oldy" | bc)   
        echo "Moving $x in x direction" >> output${k}.dat
        echo "Moving $y in y direction" >> output${k}.dat
        
      # Now we have to check out if the distance we have to move the y-coordinate is negative. If so we change the format so tha GROMACS can read it     
         if [ $(echo "$oldy < $newy"|bc) -eq 1 ]
           then echo "Distance in Y-diretion is negative: $y" >> output${k}.dat
            z=$(echo "-1 * $y" | bc)
            y=$(echo "-0$z") 
            echo "Format has been changed to: $y">> output${k}.dat           
         fi   
 
 # Here we translate the water droplet to the calculated position of the surface's center
        echo 0 | trjconv -f NVT_sam${i}_water${j}.xtc -s NVT_sam${i}_water${j}.tpr -trans ${x} ${y} 0 -o dummy_${k}ns_${l}ns.xtc -pbc atom -b ${start} -e ${ending} 

        echo "Checking out new position" >> output${k}.dat # We analyze again the coordinates of the center of the water droplet and check if now it is in the correct position
        
        echo water | g_traj -f dummy_${k}ns_${l}ns.xtc -s  NVT_sam${i}_water${j}.tpr -ox  newcoord${k}_sam${i}_water${j}.xvg -com -b ${start} -e ${ending}

       g_analyze -f newcoord${k}_sam${i}_water${j}.xvg -av | grep SS1 > newav${k}_sam${i}_water${j}
	   rm average.xvg
       g_analyze -f newcoord${k}_sam${i}_water${j}.xvg -av | grep SS2 >> newav${k}_sam${i}_water${j}
	   rm average.xvg
        

        movedx=$(awk '/SS1/ {print 1*$2}' newav${k}_sam${i}_water${j})
        movedy=$(awk '/SS2/ {print 1*$2}' newav${k}_sam${i}_water${j})
    
        
        diffx=$(echo "$movedx - $newx"|bc)
        diffy=$(echo "$movedy - $newy"|bc)
        
# diffx and diffy are the distance difference between the center of the surface and the center of the droplet. 
# We need to change the format of these numbers and put them in absolute value in order to check if they are small enough 

 # First for x
        if [ $(echo "$diffx <1"|bc) = 1 ]        
           then if [ $(echo "$diffx>0 "|bc) = 1 ] ## Checking if number is between 0 and 1
           then z=$diffx
           diffx=$(echo "0$z")
           echo "This is our new formatted diffx: $diffx"
           fi
       fi
       
       if [ $(echo "$diffx <0"|bc) = 1 ]        
           then if [ $(echo "$diffx>-1 "|bc) = 1 ] ## Checking if number is between -1 and 0    
           then z=$(echo "-1 * $diffx" |bc) 
           diffx=$(echo "0$z")
              echo "This is our new formatted diffx in absolute value: $diffx"
           fi
       fi    
 #Now for y
 
    if [ $(echo "$diffy < 1"|bc) = 1 ]        
       then    if [ $(echo "$diffy>0 "|bc) = 1 ] ## Checking if number is between 0 and 1
          then  z=$diffy
           diffy=$(echo "0$z")
           echo "This is our new formatted diffy: $diffy"
           fi
    fi
    
      if [ $(echo "$diffy < 0"|bc) = 1 ]        
           then if [ $(echo "$diffy>-1 "|bc) = 1 ] ## Checking if number is between -1 and 0    
           then z=$(echo "-1 * $diffy" |bc) # Changing the sign because we want the absolute value
           diffy=$(echo "0$z")
              echo "This is our new formatted diffx in absolute value: $diffy"
           fi
      fi
           
        
        
         if [ $(echo "$diffx < 0.001"|bc) = 1 ]
           then
           echo "This is the value of diffx inside the if statement: $diffx"
           echo "New x coordinate $movedx is equal to the surface center at $newx" 
           else echo "Error! New x coordinate $movedx is not equal to the surface center at $newx" 
         fi
         if [ $(echo "$diffy < 0.001"|bc) = 1 ]
           then 
           echo "This is the value of diffy inside the if statement: $diffy"
           echo "New y coordinate $movedy is equal to the surface center at $newy"
           else echo "Error! New y coordinate $movedy is not equal to the surface center at $newy" 
         fi
         
         
 #       echo "Creating a dummy.gro file to visualize results using the dummy.xtc file" >> output${k}.dat
        
 #       echo 0 | trjconv -f NVT_sam${i}_water${j}.xtc -s NVT_sam${i}_water${j}.tpr -trans ${x} ${y} 0 -o dummy_${k}ns_${l}ns.gro -pbc atom -b ${start} -e ${ending}


        echo "Producing density map" >> output${k}.dat # Since we have everything that's needed, we produce the density map for this time segment 
        
        echo 1 4 4 | g_densmap -f dummy_${k}ns_${l}ns.xtc -s NVT_sam${i}_water${j}.tpr  -n index${i}_${j}.ndx -amax 5 -rmax 5 -bin 0.05 -od densmap_${i}pc_${j}_${k}ns_${l}ns.dat 
        
        echo "Copying files to specific folder" >> output${k}.dat
        
        cp dummy_${k}ns_${l}ns.xtc ./sam${i}_water${j}
        cp dummy_${k}ns_${l}ns.gro ./sam${i}_water${j}
        cp coord${k}_sam${i}_water${j}.xvg ./sam${i}_water${j}
        cp av${k}_sam${i}_water${j} ./sam${i}_water${j}
        cp boxsize_sam${i}_water${j}.dat ./sam${i}_water${j}
        cp newcoord${k}_sam${i}_water${j}.xvg ./sam${i}_water${j}
        cp newav${k}_sam${i}_water${j} ./sam${i}_water${j}
        cp output${k}.dat ./sam${i}_water${j}
 
        cp densmap_${i}pc_${j}_${k}ns_${l}ns.dat ./densmaps_s${i}_w${j}
                
        rm dummy_${k}ns_${l}ns.xtc
        rm dummy_${k}ns_${l}ns.gro
        rm coord${k}_sam${i}_water${j}.xvg 
        rm av${k}_sam${i}_water${j} 
        rm boxsize_sam${i}_water${j}.dat 
        rm newcoord${k}_sam${i}_water${j}.xvg 
        rm newav${k}_sam${i}_water${j} 
        rm output${k}.dat
        
        rm densmap_${i}pc_${j}_${k}ns_${l}ns.dat
   
        
done
done
done

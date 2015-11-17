  #!/bin/bash

#These are the values for period from 10 to 12 ns of sam66_water3000  
        movedx=6.6671
        newx=6.75
        movedy=6.98666
        newy=6.928
 
 
 #These are invented values
#         movedy=6.749999
 #        newy=6.75
  #       newx=6.546
   #      movedx=6.55555

        diffx=$(echo "$movedx - $newx"|bc)
        diffy=$(echo "$movedy - $newy"|bc)
        
        echo "this is diffx: $diffx"
        echo "this is diffy: $diffy"

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

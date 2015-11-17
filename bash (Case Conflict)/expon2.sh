#!/bin/sh


        oldx=$(awk '/SS1/ {print 1*$2}' av_sam5_water1000)

        newx=$(awk 'END {print 0.5*$1}' boxsize.dat)
        
        oldy=$(awk '/SS2/ {print 1*$2}' av_sam5_water1000)
      
        newy=$(awk 'END {print 0.5*$2}' boxsize.dat)
   
        echo "The old x is $oldx"
        echo "The new x is $newx"
        echo "The old y is $oldy"
        echo "The new y is $newy"
    
        x=$(echo "$newx - $oldx" | bc)
        y=$(echo "$newy - $oldy" | bc)   
        echo "move $x in x direction"
        echo "move $y in y direction"
        
        
        a="1.21231"
b="2.22454" 
c=$(echo "$a < $b" | bc)
if [ $c == '1' ]; then 
    echo 'a is smaller than b'
    echo "this is c $c"
else 
    echo 'a is larger than b'
      echo "this is c $c"
fi
        
        uno=$(echo 6.675 | bc)
        dos=$(echo 6.674 | bc) 
        lvar=5 
        tres=$(echo "$uno - $dos" |bc)
         if [ $(echo "$tres <= 1.5"|bc) -eq 1 ]
         then echo "$tres is smaller or equal then 1.5"
       
          lvar=6
          echo "$lvar" > 'lvar.txt'
         else echo "$tres is greater then 1.5"
         lvar=6
         echo "$lvar" > 'lvar.txt'
         
        fi 

echo "is this a 6?"

lvar=$(cat lvar.txt)
echo "$lvar"
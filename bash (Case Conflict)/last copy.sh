#!/bin/sh

awk 'END {lastline = $NR; print '$lastline'}' NVT_sam5_water1000.gro > boxsize.dat


awk -F, ' {
        for(i=1;i<=NF;i++) {
               if(($i~/E\+/)||($i~/E-/))
                        $i=sprintf("%8.f",$i);            
        }
}1' OFS=, av_sam5_water1000 > av2_sam5_water1000

        oldx=$(awk '/SS1/ {print $2}' av2_sam5_water1000)
        
        printf "%7.4f\n" "the number oldx is $oldx"
 
        oldy=$(awk '/SS2/ {print $2}' av2_sam5_water1000)
        
        echo "the number oldy is $oldy"
        
        newx=$(awk 'END {print 0.5*$1}' boxsize.dat)

        echo "the number new2x is $newx"
        
        newy=$(awk 'END {print 0.5*$2}' boxsize.dat)
        
        echo "the number new2y is $newy"
              

        
        echo "The old x is $oldx"
        echo "The new x is $newx"
        echo "The old y is $oldy"
        echo "The new y is $newy"
        
        xdiff=$oldx-$newx
        ydiff=$oldy-$newy
        
        echo "The x difference is $xdiff"
        echo "The y difference is $ydiff"
        
        echo |awk -v ox="$oldx" -v nx="$newx" 'END {print $nx "   " $ox}' 
        
        awk -v a="1.234e23" -v b="9.876e14" 'BEGIN{print (a * b)}'
#!/bin/bash


for FILE in NVT3_sam11_water2000.* ;
 do NEWFILE=`echo ${FILE} | sed 's/NVT3/NVT2/g'`;
   echo ${FILE}
   echo ${NEWFILE}
   mv $FILE $NEWFILE ;
done


       # x=$(echo "$newx - $oldx" | bc)
       #NEWFILE=`echo $FILE | sed 's/ /_/g'`
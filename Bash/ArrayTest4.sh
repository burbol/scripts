#!/bin/bash
# define as first and second parameter for script:
total=0
declare -a sum[10]
limit=1.234
echo "Limit is equal to $limit"
ara=5.6
for (( i=1; i<=$limit; i++ ))
do
sum[$i]=`echo $ara | cut -d " " -f $i`
echo ${sum[$i]}
if [ $total -eq 0 ]
then
total=${sum[$i]}
else
total=`echo $total + echo ${sum[$i]} | bc`
fi
done
echo "Total Hard Disk Size is $total"
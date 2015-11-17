#!/bin/bash

echo "Hello World!" 
n1=21
n2=28
k=25
k2=$(($k+5))
echo "This is n1",  $n1, "This is n2",  $n2
n3=$(echo "$n1/$n2"|bc)
echo $n3
echo "This is k2", $k2
# the round function:
round()
{
echo $(printf %.$2f $(echo "scale=$2;(((10^$2)*$1)+0.5)/(10^$2)" | bc))
};

echo $(round $n1/$n2 2);

nanosecs1=$(echo $(round $k/10 1))
nanosecs2=$(echo $(round $k2/10 1))

echo $nanosecs1, $nanosecs2

#echo "$n1/$n2" | python -c "print round(float(raw_input()))"

while [[ $count -le 95 ]]
do
   echo "$count"
   count=$(( count + 5 ))
done

i=0
  if [ $i -eq 0 ]
  then 
    n=3    
  else 
    n=4
  fi
echo $n
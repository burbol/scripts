#!/bin/bash

#####################################################################
# Default scale used by float functions.

float_scale=5


#####################################################################
# Evaluate a floating point number expression.

function float_eval()
{
    local stat=0
    local result=0.0
    if [[ $# -gt 0 ]]; then
        result=$(echo "scale=$float_scale; $*" | bc -q 2>/dev/null)
        stat=$?
        if [[ $stat -eq 0  &&  -z "$result" ]]; then stat=1; fi
    fi
    echo $result
    return $stat
}


##EXAMPLE of usage:
  #float_eval '12.0 / 3.0'
  #a=12.0
  #b=3.0
  #c=$(float_eval "$a / $b")
 
 
###################################################################### 
a=1.234
j=3;
#for i in 1 2 3 4
#do
   #grofile[$i]=$(float_eval "$a")
   #echo $a
   #grofile[$i]=$a
   #echo $(float_eval "$((grofile[$i]))")
   #new[$i]=$(echo "$((grofile[$i]))/2" |bc)
   #echo $((new[$i]))
#done 

#######################################################################
###### From other source######
#!/bin/ksh

for i in 1 2 3 4 5 6 7 8
do
  a[$i]=$(bc <<-%
  sc=6
  float_eval "1 / $i"
  %)
  printf "%d %f\n" $i ${a[$i]}
done

#######################################################################
###### From other source######
#sum=0
#for i in "${thearray[@]}"; do
#  sum=$(echo $sum + $i | bc -l);
#done
#echo "Sum = ${sum}"
#average=$(echo $sum / ${#thearray[@]} | bc -l)
#echo "Average = ${average}"

#!/bin/bash

n1 = 0
n2 = 5
n3 = 10
n4 = 15

printf -v v1 "%.1f" "$n1"
printf -v v2 "%.1f" "$n2"
printf -v v3 "%.1f" "$n3"

echo $v1, $v2, $v3 
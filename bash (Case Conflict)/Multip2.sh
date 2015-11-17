#!/bin/bash

y=$(echo "-.2314"| bc)
z=$(echo "-1 * $y" |bc) 
echo "This is z: $z"
y=$(echo "-0$z")
echo " This is y: $y"
#!/bin/bash

  x=$(echo "-.123" | bc)

  y=$(echo "-1 * $x" | bc)
  z=$(echo "-0$y") 
  echo "$z"
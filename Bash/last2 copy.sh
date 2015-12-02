#!/bin/sh

awk '/SS1/ {print "x "$2}
  END  {print "hola"}
' av_sam5_water1000

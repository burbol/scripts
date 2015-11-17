#!/bin/sh

awk 'END {lastline = $NR; print '$lastline'}' NVT_sam5_water1000.gro > boxsize.dat

awk '{print '$2' }' av_sam${i}_water${j}.xvg > boxsize.dat

awk 'END {lastline = $NR; print '$lastline'}' NVT_sam5_water1000.gro > boxsize.dat
awk 'END {print '$2'}' boxsize.dat > boxsize.dat


SS1  5.55  bla bla
SS2  6.00  bla bla
3 4 7


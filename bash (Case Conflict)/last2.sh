#!/bin/sh

awk 'END {print '$2' }' av_sam5_water1000 | read coordx coordy
echo "$coordx $coordy"
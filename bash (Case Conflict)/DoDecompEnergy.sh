#!/bin/sh 


rm DecompositionAll.dat
rm DecompositionAll*.dat

for i in 35 37 39 41 43 45 47 49 51 53 55; do

cd ../static;

mkdir DecomposEnergy
echo 15 | trjconv -f Static${i}.trr -s Static${i}.tpr -n ../StartDecomposition/index -o DecomposEnergy/Static${i}_PS.trr
echo 15 | trjconv -f Static${i}.gro -s Static${i}.tpr -n ../StartDecomposition/index -o DecomposEnergy/Static${i}_PS.gro

echo 16 | trjconv -f Static${i}.trr -s Static${i}.tpr -n ../StartDecomposition/index -o DecomposEnergy/Static${i}_SW.trr
echo 16 | trjconv -f Static${i}.gro -s Static${i}.tpr -n ../StartDecomposition/index -o DecomposEnergy/Static${i}_SW.gro

echo 17 | trjconv -f Static${i}.trr -s Static${i}.tpr -n ../StartDecomposition/index -o DecomposEnergy/Static${i}_PW.trr
echo 17 | trjconv -f Static${i}.gro -s Static${i}.tpr -n ../StartDecomposition/index -o DecomposEnergy/Static${i}_PW.gro

echo 1  | trjconv -f Static${i}.trr -s Static${i}.tpr -n ../StartDecomposition/index -o DecomposEnergy/Static${i}_PP.trr
echo 1  | trjconv -f Static${i}.gro -s Static${i}.tpr -n ../StartDecomposition/index -o DecomposEnergy/Static${i}_PP.gro

echo 13 | trjconv -f Static${i}.trr -s Static${i}.tpr -n ../StartDecomposition/index -o DecomposEnergy/Static${i}_WW.trr
echo 13 | trjconv -f Static${i}.gro -s Static${i}.tpr -n ../StartDecomposition/index -o DecomposEnergy/Static${i}_WW.gro


cd DecomposEnergy

grompp -c Static${i}_PS.gro -f ../../StartDecomposition/PS.mdp -p ../../StartDecomposition/PS.top -o PS${i}.tpr
mdrun -rerun Static${i}_PS.trr -s PS${i}.tpr -deffnm Static${i}_PS

grompp -c Static${i}_SW.gro -f ../../StartDecomposition/SW.mdp -p ../../StartDecomposition/SW.top -o SW${i}.tpr
mdrun -rerun Static${i}_SW.trr -s SW${i}.tpr -deffnm Static${i}_SW

grompp -c Static${i}_PW.gro -f ../../StartDecomposition/PW.mdp -p ../../StartDecomposition/PW.top -o PW${i}.tpr
mdrun -rerun Static${i}_PW.trr -s PW${i}.tpr -deffnm Static${i}_PW

grompp -c Static${i}_PP.gro -f ../../StartDecomposition/PP.mdp -p ../../StartDecomposition/PP.top -o PP${i}.tpr
mdrun -rerun Static${i}_PP.trr -s PP${i}.tpr -deffnm Static${i}_PP

grompp -c Static${i}_WW.gro -f ../../StartDecomposition/WW.mdp -p ../../StartDecomposition/WW.top -o WW${i}.tpr
mdrun -rerun Static${i}_WW.trr -s WW${i}.tpr -deffnm Static${i}_WW

echo Potential  Kinetic-En.  Total-Energy | g_energy -f Static${i}_PP.edr -b 5000 -o Energy${i}_PP
echo Potential  Kinetic-En.  Total-Energy | g_energy -f Static${i}_PS.edr -b 5000 -o Energy${i}_PS
echo Potential  Kinetic-En.  Total-Energy | g_energy -f Static${i}_PW.edr -b 5000 -o Energy${i}_PW
echo Potential  Kinetic-En.  Total-Energy | g_energy -f Static${i}_SW.edr -b 5000 -o Energy${i}_SW
echo Potential  Kinetic-En.  Total-Energy | g_energy -f Static${i}_WW.edr -b 5000 -o Energy${i}_WW


g_analyze -f Energy${i}_PP.xvg -av | grep SS1 > EnergyAllPP.dat
g_analyze -f Energy${i}_PS.xvg -av | grep SS1 > EnergyAllPS.dat
g_analyze -f Energy${i}_PW.xvg -av | grep SS1 > EnergyAllPW.dat
g_analyze -f Energy${i}_SW.xvg -av | grep SS1 > EnergyAllSW.dat
g_analyze -f Energy${i}_WW.xvg -av | grep SS1 > EnergyAllWW.dat


rm *.*#
rm *.*~

cd ../
rm *.*#
rm *.*~
cd ../EnergyDecomp

awk '{print "POS", $2, $3, $4}' ../static/DecomposEnergy/EnergyAllPP.dat >> DecompositionAllPP.dat;
sed s/POS/${i}/g DecompositionAllPP.dat > src;
mv src DecompositionAllPP.dat;

awk '{print "POS", $2, $3, $4}' ../static/DecomposEnergy/EnergyAllPW.dat >> DecompositionAllPW.dat;
sed s/POS/${i}/g DecompositionAllPW.dat > src;
mv src DecompositionAllPW.dat;

awk '{print "POS", $2, $3, $4}' ../static/DecomposEnergy/EnergyAllWW.dat >> DecompositionAllWW.dat;
sed s/POS/${i}/g DecompositionAllWW.dat > src;
mv src DecompositionAllWW.dat;

awk '{print "POS", $2, $3, $4}' ../static/DecomposEnergy/EnergyAllPS.dat >> DecompositionAllPS.dat;
sed s/POS/${i}/g DecompositionAllPS.dat > src;
mv src DecompositionAllPS.dat;

awk '{print "POS", $2, $3, $4}' ../static/DecomposEnergy/EnergyAllSW.dat >> DecompositionAllSW.dat;
sed s/POS/${i}/g DecompositionAllSW.dat > src;
mv src DecompositionAllSW.dat;

done
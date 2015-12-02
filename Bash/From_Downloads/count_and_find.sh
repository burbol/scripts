# Find value and display lines with occurrences

cat start11.gro |grep "C " HeadGroups11.dat

# Find value and copy lines with occurrences to file

grep "C " start11.gro >> HeadGroups11.dat
 
grep "C1 " start.gro >> C1HeadGroups0.dat
 
grep "O1 " start17.gro >> O1HeadGroups17.dat

# Count occurrences and display number

grep -c "H42 " start11.gro 

grep -c "O1 " start11.gro 
 
grep -c "H41 " start.gro 
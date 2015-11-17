#!/usr/bin/env python

import os
import math
import fnmatch # for filtering the file list to get the top,gro,mdp,ndx files.

Dir = './temp/'
Filename = 'jobscript_'

N = 50

walltime = 20
Email = 'bla'
NumberOfProcessorsPerNode = 5


ProductionFilename = 'topol'
TrajFilename = 'traj'
EnergyFilename = 'ener'
LogFilename = 'md'
MaxH = 20
TopFile = 'topology.top'
CurGro = 'input_gro.gro'

for i in range(N):
    JobOut = open(Dir + '/' + Filename + str(i),'w')
    JobOut.write('#!/bin/bash\n')
    JobOut.write('#SBATCH --job-name=' + str(i) + '\n')
    JobOut.write('#SBATCH --mail-type=end\n')
    JobOut.write('#SBATCH --mem=2048\n')
    JobOut.write('#SBATCH --time=' + str(walltime)  + '\n')
    JobOut.write('#SBATCH --mail-user=' + Email + '\n')
    JobOut.write('#SBATCH --ntasks=1\n')
    JobOut.write('#SBATCH --cpus-per-task=' + str(NumberOfProcessorsPerNode) + '\n')
    JobOut.write('\n')
    # first jobscript should contain grompp:
    if i == 0:
        JobOut.write('grompp -f production.mdp -n index.ndx -p ' + TopFile + ' -c ' + CurGro + 
    	    ' -o ' + ProductionFilename + '.tpr -maxwarn 1\n')
        JobOut.write('\n')
    # first jobscript does not have '-cpi state.cpt' in mdrun call:
    if i == 0:
        JobOut.write('mdrun -nt 8 -s ' + ProductionFilename + '.tpr -o ' + 
    	    TrajFilename + ' -x ' + TrajFilename + ' -e ' + 
    	    EnergyFilename + ' -g ' + LogFilename + ' -px -pf -maxh ' +
    	    str(MaxH) + '\n')
    else:
    	JobOut.write('mdrun -nt 8 -s ' + ProductionFilename + '.tpr -o ' + 
    		    TrajFilename + ' -x ' + TrajFilename + ' -e ' + 
    		    EnergyFilename + ' -g ' + LogFilename + ' -px -pf -maxh ' +
    		    str(MaxH) + ' -cpi state.cpt\n')	
    JobOut.write('\n')
    # first 49 jobscripts check if confout.gro exists and - if not - submit the next jobscript
    if i < N-1: # note that we count from 0 to N-1, so the 50th jobscript has i = 49
        JobOut.write('NotDone=false\n')
        JobOut.write('# check if the file "confout.gro" is there\n')
        JobOut.write('if [ ! -f confout.gro ]; then\n')
        JobOut.write('  NotDone=true\n')
        JobOut.write('  echo "confout.gro not found in current directory."\n')
        JobOut.write('fi\n')
        JobOut.write('\n')
        JobOut.write('\n')
        JobOut.write('if [ "$NotDone" = true ]; then\n')
        JobOut.write('sbatch ' + str(Filename) + str(i+1) + '\n')
        JobOut.write('fi\n')
        JobOut.write('\n')
        JobOut.close()


#!/usr/bin/env python

# !!IMPORTANT!! When not running tests, last PRINTED line should NOT be a comment (it's a "cp" command to copy the whole simulations folder)

import os
import math
import fnmatch # for filtering the file list to get the top,gro,mdp,ndx files.

Dir = '.' # directory to save script 
SimDir = '/scratch/eixeres/s17_w1000/'   #directory from where to run simulation
BackupDir = '/home/eixeres/simulations/finished/'   # directory to copy (backup) simulation output
Filename = 'NVT_sam17_water1000_next'  # name of files
Scriptname = 's17_w1000_8test'  # name of script (when NOT testing, set equal to "Filename")
Jobname = '8tests17w1000' # Name of job

N = 1  # NUMBER OF SCRIPTS TO SUBMIT = total time running will be this number times the Walltime

Walltime = '0:30:00'  # Walltime in format 'h:mm:ss'
SimTime = 0.5   # maximum run length of simulation in hours, it can be also a fraction of an hour
Email = 'laila.e@fu-berlin.de'

NodesNum = 1  # soroban = 12 cpus/node ; sheldon-ng = 8 cpus/node
CpuNum = 8
Memory = 1024
SimLen = 20  # total simulation time length in ps
Partition = 'test' # use 'test' for testing and 'main' otherwise

#next 3 lines only useful when runnig also grompp
#IndexFile = ' '
#FirstGro = 'sam5_water216'
#FirstMdp = 'Mini.mdp'



for i in range(N):
    JobOut = open(Dir + '/' + Scriptname + str(i),'w')
    JobOut.write('#!/bin/bash\n')
    JobOut.write('\n')
    JobOut.write('#SBATCH -p '+ Partition +'\n')
    JobOut.write('\n')
    JobOut.write('#SBATCH --mem=' + str(Memory) +'\n')
    JobOut.write('#SBATCH --job-name=' + str(i) + Jobname + '\n')
    JobOut.write('#SBATCH --output=' + Scriptname + '.out\n')
    JobOut.write('\n')
    JobOut.write('#SBATCH --mail-user=' + Email + '\n')
    JobOut.write('#SBATCH --mail-type=end\n')
    JobOut.write('#SBATCH --mail-type=fail\n')
    JobOut.write('\n')
    JobOut.write('#SBATCH --ntasks=' + str(CpuNum)+'\n')
    JobOut.write('#SBATCH --nodes=' + str(NodesNum) +'\n')
    JobOut.write('\n')
    JobOut.write('#SBATCH --time=' + Walltime  + '\n')
    JobOut.write('\n')
    JobOut.write('module load slurm \n')
    JobOut.write('\n')
    JobOut.write('module load gromacs/openmpi/gcc/64/4.5.4 \n')
    JobOut.write('\n')
    JobOut.write('STARTTIME=$(date +%s)\n' + '\n'+ '#use sleep or testing... \n')
    JobOut.write('\n')
    JobOut.write('cd ' + SimDir + '\n')
    JobOut.write('\n')
    
    # first jobscript changes "mdp" options with tpbconv to extend simulation
    if i == 0:
        JobOut.write('tpbconv -s ' + Filename + '.tpr  -until ' + str(SimLen) + ' -o ' + Filename + '_next.tpr \n')
        JobOut.write('\n')
        JobOut.write('mpirun -np ' + str(CpuNum) + ' mdrun -s ' + Filename + '_next.tpr -cpi  '+ Filename + '.cpt -maxh ' +
    	    str(SimTime) + '\n')
    else:
        JobOut.write('mpirun -np ' + str(CpuNum) + ' mdrun -s ' + Filename + '.tpr -cpi  '+ Filename + '.cpt -maxh ' +
    	    str(SimTime) + '\n')
    	    
    JobOut.write('\n')
    
    # first N-1 jobscripts check if runtime is less then 15 sec, if not, submit next script 
    if i < N-1: # note that we count from 0 to N-1, so the 50th jobscript has i = 49
        JobOut.write( 'RUNTIME=$(($(date +%s)-$STARTTIME))\n' + '\n'+ 'echo \"the job took $RUNTIME seconds...\"\n' + '\n' + 'if [[ $RUNTIME -lt 10 ]]; then\n' +
                    '   echo "job took less than 10 seconds to run, aborting."\n' + '   exit\n' + 'else\n' + '   echo "everything fine..."\n' + '   sbatch ' + Scriptname + str(i+1) + '\n' + '   cp -r ' + SimDir + ' ' + BackupDir + '\n') 
        JobOut.write('   fi\n')
        JobOut.write('\n')
        
    # after the last run, we also want to backup the simulation files
    elif i == N-1:
        JobOut.write('#cp -r ' + SimDir + ' ' + BackupDir + '\n')
    JobOut.close()


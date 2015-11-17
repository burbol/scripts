#!/usr/bin/env python

# !!IMPORTANT!! If files don't need to be copied, last PRINTED line should be a comment (it's a "cp" command to copy the whole simulations folder)

import os
import math
import fnmatch # for filtering the file list to get the top,gro,mdp,ndx files.


Dir = '.' # directory to save script 
SimDir = '/home/eixeres/files_for_laila/s5_w9000'   #directory from where to run simulation
BackupDir = '/home/eixeres/files_for_laila/'   # directory to copy (backup) simulation output
Filename = 'NVT_sam5_water9000_next'  # name of files
Scriptname = 's5_w9000_next'  # name of script (when NOT testing, set equal to "Filename")
Jobname = 's5_w9000_next' # Name of job

N = 20  # NUMBER OF SCRIPTS TO SUBMIT = total time running will be this number times the Walltime
Walltime = '24:00:00'  # Walltime in format 'h:mm:ss'
SimTime = 24   # maximum run length of simulation in hours, it can be also a fraction of an hour

Email = 'laila.e@fu-berlin.de'

NodesNum = 2  # soroban = 12 cpus/node ; sheldon/sheldon-ng = 8 cpus/node
CpuNum = 12  # ratio cpus/node!!!
SimLen = 60000  # total simulation time length in ps

#next 3 lines only useful when runnig also grompp
#IndexFile = ' '
#FirstGro = 'sam5_water216'
#FirstMdp = 'Mini.mdp'


for i in range(N):

    JobOut = open(Dir + '/' + Scriptname + str(i),'w')

    JobOut.write('#!/bin/bash\n')

    JobOut.write('#PBS -N ' + Jobname  + str(i) +  '\n')

    JobOut.write('#PBS -o ' + Scriptname + '.stdout\n')
    
    JobOut.write('#PBS -e ' + Scriptname + '.stderr\n')

    JobOut.write('#PBS -M ' + Email + '\n')

    JobOut.write('\n')

    JobOut.write('#PBS-l nodes=' + str(NodesNum) + ':ppn=' + str(CpuNum)+'\n')

    JobOut.write('#PBS -l walltime=' + Walltime  + '\n')

    JobOut.write('\n')

    JobOut.write('cd $PBS_O_WORKDIR \n')

    JobOut.write('NPROCS=`wc -l < $PBS_NODEFILE`\n')

    JobOut.write('\n')

    

    # first jobscript changes "mdp" options with tpbconv to extend simulation

    if i == 0:

        JobOut.write('tpbconv -s ' + Filename + '.tpr  -until ' + str(SimLen) + ' -o ' + Filename + '_next.tpr \n')

        JobOut.write('mpirun -np ' + str(CpuNum) + ' mdrun -s ' + Filename + '_next.tpr -cpi  '+ Filename + '.cpt -maxh ' +

    	    str(SimTime) + '\n')

    else:

        JobOut.write('mpirun -np ' + str(CpuNum) + ' mdrun -s ' + Filename + '.tpr -cpi  '+ Filename + '.cpt -maxh ' +

    	    str(SimTime) + '\n')

    	    
    JobOut.write('\n')
    
    if i < N-1:
    
    	JobOut.write('qsub ' + Scriptname + str(i+1) + '\n')
    
    elif i == N-1:

        JobOut.write('#cp -r ' + SimDir + ' ' + BackupDir + '\n')
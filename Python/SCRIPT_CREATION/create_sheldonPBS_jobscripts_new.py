#!/usr/bin/env python

# !!IMPORTANT!! If files don't need to be copied, last PRINTED line should be a comment (it's a "cp" command to copy the whole simulations folder)

import os
import math
import fnmatch # for filtering the file list to get the top,gro,mdp,ndx files.

N = 40  # NUMBER OF SCRIPTS TO SUBMIT = total time running will be this number times the Walltime
Walltime = '24:00:00'  # Walltime in format 'h:mm:ss'
SimTime = 24   # maximum run length of simulation in hours, it can be also a fraction of an hour

Email = 'laila.e@fu-berlin.de'

CpuNum = 8  # ratio cpus/node!!! # soroban = 12 cpus/node ; sheldon/sheldon-ng = 8 cpus/node

#next 3 lines only useful when runnig also grompp
#IndexFile = ' '
#FirstGro = 'sam5_water216'
#FirstMdp = 'Mini.mdp'

#pc = [0,5,17]
#molec = [1000, 6500,9000]
pc = [17]
molec = [1000,6500]

Simns = {1000:60, 6500:80, 9000:80, 4000:40, 5000:40, 7000:80 , 8000:80}  # total simulation time length in ns
#Nodes = {1000:2, 6500:3, 9000:3} 
Nodes = {1000:1, 6500:1, 9000:1, 4000:1, 5000:1, 10000:1,7000:1, 8000:1} 

for i in pc:
	for j in molec: 
		NodesNum = Nodes[j]
		SimLen = Simns[j]*1000  # total simulation time length in ps
		
		#Dir = './sheldon/s'+str(i)+'_w'+str(j)+'_new' # directory to save script  ##for s17 w1000 & w6500
		Dir = './sheldon/s'+str(i)+'_w'+str(j) # directory to save script 
		SimDir = '/home/eixeres/files_for_laila/s'+str(i)+'_w'+str(j)   #directory from where to run simulation
		BackupDir = '/home/eixeres/files_for_laila/'   # directory to copy (backup) simulation output
		Filename = 'NVT_sam'+str(i)+'_water'+str(j)+'_next'  # name of files
		Scriptname = 's'+str(i)+'_w'+str(j)+'_new'  # name of script (when NOT testing, set equal to "Filename")

		os.system("mkdir "+ Dir)
		for k in range(N):

			Jobname = 's'+str(i)+'_w'+str(j)+'_next'+ str(k)  # Name of job

			JobOut = open(Dir + '/' + Scriptname + str(k),'w')

			JobOut.write('#!/bin/bash\n')

			JobOut.write('#PBS -N ' + Jobname  + str(k) +  '\n')

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
			
			JobOut.write('source /net/opt/gromacs/openmpi1.4.5/4.5.4/bin/GMXRC\n')
			#JobOut.write('source /net/opt/gromacs/single/thread/4.6.5/bin/GMXRC\n')
			
			JobOut.write('\n')

	

			# first jobscript changes "mdp" options with tpbconv to extend simulation

			if k == 0:

				JobOut.write('tpbconv -s ' + Filename + '.tpr  -until ' + str(SimLen) + ' -o ' + Filename + '_new.tpr \n')

				#JobOut.write('mpirun -np ' + str(CpuNum*NodesNum) + ' /home/shavkat/GMX/bin/g_tune_pme_mpi -s ' + Filename + '_next.tpr -cpo  '+ Filename + '.cpt -so ' + Filename + '.tpr -launch\n')
				JobOut.write('mdrun -nt '+ str(CpuNum*NodesNum)+' -cpi '+ Filename + '.cpt -s ' + Filename + '_new.tpr -maxh ' +str(SimTime) + ' -v -testverlet\n')
				#JobOut.write('mdrun -nt '+ str(CpuNum*NodesNum)+' -cpi '+ Filename + '.cpt -s ' + Filename + '.tpr -maxh ' +str(SimTime) + ' -v -testverlet\n')

			else:

				JobOut.write('mdrun -nt '+ str(CpuNum*NodesNum)+' -cpi '+ Filename + '.cpt -s ' + Filename + '.tpr -maxh ' +str(SimTime) + ' -v -testverlet\n')

				JobOut.write( 'RUNTIME=$(($(date +%s)-$STARTTIME))\n' + '\n'+ 'echo \"the job took $RUNTIME seconds...\"\n' + '\n' + 'if [[ $RUNTIME -lt 10 ]]; then\n' +

							'   echo "job took less than 10 seconds to run, aborting."\n' + '   exit\n' + 'else\n' + '   echo "everything fine..."\n' + '   qsub ' + Scriptname + str(k+1) + '\n') 
				
				JobOut.write('   fi\n')

				JobOut.write('\n')


			
			JobOut.write('\n')
	
			if k < N-1:
	
				JobOut.write('qsub ' + Scriptname + str(k+1) + '\n')
	
			elif k == N-1:

				JobOut.write('#cp -r ' + SimDir + ' ' + BackupDir + '\n')
		JobOut.close()
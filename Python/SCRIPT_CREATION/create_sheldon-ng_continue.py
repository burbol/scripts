#!/usr/bin/env python

# !!IMPORTANT!! When not running tests, last PRINTED line should NOT be a comment (it's a "cp" command to copy the whole simulations folder)

import os
import math
import fnmatch # for filtering the file list to get the top,gro,mdp,ndx files.


#N = 40  # NUMBER OF SCRIPTS TO SUBMIT = total time running will be this number times the Walltime

Days=12
Walltime = str*(Days)+'-0:00:00'  # Walltime in format 'h:mm:ss'
SimTime = 24*Days   # maximum run length of simulation in hours, it can be also a fraction of an hour

Email = 'laila.e@fu-berlin.de'

Cpus = 8  # soroban = 12 cpus/node ; sheldon-ng = 8 cpus/node
Memory = 2048
SimLen = 60000  # total simulation time length in ps
Partition = 'main' # use 'test' for testing and 'main' otherwise 

#next 3 lines only useful when runnig also grompp
#IndexFile = ' '
#FirstGro = 'sam17_water216'
#FirstMdp = 'Mini.mdp'


Simns = {216:60, 1000:80, 2000:80, 3000:80, 4000:80, 5000:100, 6500:100, 7000:100, 8000:120, 9000:120, 10000:120}  # total simulation time length in ns
Nodes = {216:2, 1000:3, 2000:3, 3000:3, 4000:3, 5000:3, 6500:3, 7000:5, 8000:6, 9000:6, 10000:6} 
#Nodes = {1000:1, 2000:1, 3000:1, 5000:1} # FOR TESTING!!

pc = [21]
#molec = [1000, 2000, 3000, 5000] # FOR TESTING!!
molec = [1000, 2000, 3000, 4000, 5000, 6500, 7000, 8000, 9000, 10000]


BackupDir = '/home/eixeres/NewVersion3/'   # directory to copy (backup) simulation output

for i in pc:
	for j in molec:
		
		Dir = '/Users/burbol/MEGAsync/scripts/Python/SCRIPT_CREATION/sheldon-ng/s'+str(i)+'_w'+str(j)
		SimDir = '/scratch/eixeres/NewVersion3/s'+str(i)+'_w'+str(j)   #directory from where to run simulation
		Filename = 'NVT_sam'+str(i)+'_water'+str(j)  # name of files
		Scriptname = 's'+str(i)+'_w'+str(j)  # name of script (when NOT testing, set equal to "Filename")
		
		Startfile='sam'+str(i)+'_water'+str(j)
		Minifile='Mini_sam'+str(i)+'_water'+str(j)
		NVTfile = 'NVT_sam'+str(i)+'_water'+str(j)
		topfile= str(i)+'pc_'+str(j)+'.top'
		Minimdp = 'Mini.mdp'
		NVTmdp = 'NVT_'+str(Simns[j])+'ns.mdp'
		
		SimLen = Simns[j]*1000  # total simulation time length in ps
		NodesNum = Nodes[j]
		CpuNum = NodesNum*Cpus
		N = Simns[j]  # NUMBER OF SCRIPTS TO SUBMIT: 
	
		os.system("mkdir "+Dir)
		for k in range(1,N):
		
			#os.system("mkdir " + Dir)

			Jobname = 's'+str(i)+'_w'+str(j)+ str(k)  # Name of job
			JobOut = open(Dir + '/' + Scriptname +'_'+ str(k),'w')

			JobOut.write('#!/bin/bash\n')

			JobOut.write('\n')

			JobOut.write('#SBATCH -p '+ Partition +'\n')

			JobOut.write('\n')

			JobOut.write('#SBATCH --mem=' + str(Memory) +'\n')

			JobOut.write('#SBATCH --job-name=' + str(k) + Jobname + '\n')

			JobOut.write('#SBATCH --output=' + Scriptname + '.out\n')

			JobOut.write('\n')

			JobOut.write('#SBATCH --mail-user=' + Email + '\n')

			JobOut.write('#SBATCH --mail-type=end\n')

			JobOut.write('#SBATCH --mail-type=fail\n')

			JobOut.write('\n')

			#JobOut.write('#SBATCH --ntasks=' + str(CpuNum)+'\n')
			JobOut.write('#SBATCH --tasks-per-node=' + str(Cpus)+'\n')

			JobOut.write('#SBATCH --nodes=' + str(NodesNum) +'\n')

			JobOut.write('\n')

			JobOut.write('#SBATCH --time=' + str(Walltime)  + '\n')

			JobOut.write('\n')

			JobOut.write('module load slurm \n')

			JobOut.write('\n')
			
			#JobOut.write('module load gromacs/single/openmpi1.4.5/4.6.5\n')
			JobOut.write('module load gromacs/single/thread/4.6.5\n') # creo q este es mejor!!

			JobOut.write('\n')

			JobOut.write('STARTTIME=$(date +%s)\n' + '\n'+ '#use sleep or testing... \n')

			JobOut.write('\n')

			JobOut.write('cd ' + SimDir + '\n')

			JobOut.write('\n')

	
	

						# first jobscript changes "mdp" options with tpbconv to extend simulation

			if k == 0:

				#JobOut.write('tpbconv -s ' + Filename + '.tpr  -until ' + str(SimLen) + ' -o ' + Filename + '.tpr \n')
				JobOut.write('mpirun -np ' + str(CpuNum) + ' mdrun -cpi '+ Filename + '.cpt -s ' + Filename + '.tpr -deffnm   '+ Filename + ' -maxh ' +str(SimTime) + ' -v -testverlet\n')

			else:

				JobOut.write('mpirun -np ' + str(CpuNum) + ' mdrun -cpi '+ Filename + '.cpt -s ' + Filename + '.tpr -deffnm   '+ Filename + ' -maxh ' +str(SimTime) + ' -v -testverlet\n')
				print('mpirun -np ' + str(CpuNum) + ' mdrun -cpi '+ Filename + '.cpt -s ' + Filename + '.tpr -deffnm   '+ Filename + ' -maxh ' +str(SimTime) + ' -v -testverlet\n')
			

			JobOut.write('\n')

	

			# first N-1 jobscripts check if runtime is less then 15 sec, if not, submit next script 

			if k < N-1: # note that we count from 0 to N-1, so the 50th jobscript has i = 49

				JobOut.write( 'RUNTIME=$(($(date +%s)-$STARTTIME))\n' + '\n'+ 'echo \"the job took $RUNTIME seconds...\"\n' + '\n' + 'if [[ $RUNTIME -lt 10 ]]; then\n' +

							'   echo "job took less than 10 seconds to run, aborting."\n' + '   exit\n' + 'else\n' + '   echo "everything fine..."\n' + '   sbatch ' + Scriptname +'_'+ str(k+1) + '\n') 

				JobOut.write('   fi\n')

				JobOut.write('\n')

		

			# after the last run, we also want to backup the simulation files

			elif k == N-1:

				JobOut.write("mkdir " + BackupDir + '\n')
				JobOut.write('cp -r ' + SimDir + ' ' + BackupDir + '\n')

			JobOut.close()




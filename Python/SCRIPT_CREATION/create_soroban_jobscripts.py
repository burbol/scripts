#!/usr/bin/env python

# !!IMPORTANT!! When not running tests, last PRINTED line should NOT be a comment (it's a "cp" command to copy the whole simulations folder)

import os
import math
import fnmatch # for filtering the file list to get the top,gro,mdp,ndx files.


#N = 40  # NUMBER OF SCRIPTS TO SUBMIT = total time running will be this number times the Walltime

Walltime = '24:00:00'  # Walltime in format 'h:mm:ss'
SimTime = 24   # maximum run length of simulation in hours, it can be also a fraction of an hour

Email = 'laila.e@fu-berlin.de'

NodesNum =  2 # soroban = 12 cpus/node ; sheldon-ng = 8 cpus/node
Cpus = 12
Memory = 1024
SimLen = 60000  # total simulation time length in ps
Partition = 'main' # use 'test' for testing and 'main' otherwise 

#next 3 lines only useful when runnig also grompp
#IndexFile = ' '
#FirstGro = 'sam17_water216'
#FirstMdp = 'Mini.mdp'


Simns = {216:40, 1000:80, 2000:80, 3000:80, 4000:80, 5000:100, 6500:100, 7000:100, 8000:120, 9000:120, 10000:120}  # total simulation time length in ns
Nodes = {216:2, 1000:2, 2000:2, 3000:2, 4000:2, 5000:2, 6500:3, 7000:3, 8000:3, 9000:3, 10000:3} 

#pc = [0,5,11,17,33,50]
pc = [21, 25]
molec = [216, 1000, 2000, 3000, 4000, 5000, 6500, 7000, 8000, 9000, 10000]

for i in pc:
	for j in molec:
		
		
		Dir = './soroban/NEWVersionBIG/s'+str(i)+'_w'+str(j) # directory to save script
		SimDir = '/scratch/eixeres/s'+str(i)+'_w'+str(j)   #directory from where to run simulation
		BackupDir = '/home/eixeres/simulations/finished/'   # directory to copy (backup) simulation output
		Filename = 'NVT_sam'+str(i)+'_water'+str(j)  # name of files
		Scriptname = 's'+str(i)+'_w'+str(j)  # name of script (when NOT testing, set equal to "Filename")
		
		Startfile='sam'+str(i)+'_water'+str(j)
		Minifile='Mini_sam'+str(i)+'_water'+str(j)
		NVTfile = 'NVT_sam'+str(i)+'_water'+str(j)
		topfile= str(i)+'pc_'+str(j)+'.top'
		Minimdp = 'Mini.mdp'
		NVTmdp = 'NVT_'+str(Simns[j])+'ns.mdp'
		
		SimLen = Simns[j]*1000  # total simulation time length in p
		NodesNum = Nodes[j]
		CpuNum = NodesNum*Cpus
		N = Simns[j]  # NUMBER OF SCRIPTS TO SUBMIT: 2 scripts each ns 
		os.system("mkdir " + Dir)
		for k in range(N):

			Jobname = 's'+str(i)+'_w'+str(j)+ str(k)  # Name of job

			JobOut = open(Dir + '/' + Scriptname + str(k),'w')

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

			JobOut.write('#SBATCH --ntasks=' + str(CpuNum)+'\n')

			JobOut.write('#SBATCH --nodes=' + str(NodesNum) +'\n')

			JobOut.write('\n')

			JobOut.write('#SBATCH --time=' + Walltime  + '\n')

			JobOut.write('\n')

			JobOut.write('module load slurm \n')

			JobOut.write('\n')

			#JobOut.write('module load gromacs/openmpi/gcc/64/4.5.4\n')
			JobOut.write('module load gromacs/intelmpi/intel/64/4.6\n')

			JobOut.write('\n')

			JobOut.write('STARTTIME=$(date +%s)\n' + '\n'+ '#use sleep or testing... \n')

			JobOut.write('\n')

			JobOut.write('cd ' + SimDir + '\n')

			JobOut.write('\n')

	

			# first jobscript changes "mdp" options with tpbconv to extend simulation

			if k == 0:

				
				JobOut.write('grompp -f ' + Minimdp + ' -c ' + Startfile + '.gro -p ' + topfile + ' -o '+ Minifile +'.tpr -maxwarn 1\n')

				JobOut.write('mpirun -np ' + str(CpuNum) + ' mdrun -deffnm ' + Minifile + ' -maxh ' + str(SimTime) + '\n')
				
				JobOut.write('grompp -f ' + NVTmdp + ' -c ' + Minifile + '.gro -p ' + topfile + ' -o '+ NVTfile +'.tpr -maxwarn 1\n')
				
				
				#JobOut.write('tpbconv -s ' + Filename + '.tpr  -until ' + str(SimLen) + ' -o ' + Filename + '_next.tpr \n')

				JobOut.write('\n')

				JobOut.write('mpirun -np ' + str(CpuNum) + ' mdrun -s ' + Filename + '.tpr -deffnm  '+ Filename + ' -maxh ' +

					str(SimTime) + '\n')
				print('mpirun -np ' + str(CpuNum) + ' mdrun -s ' + Filename + '.tpr -deffnm  '+ Filename + ' -maxh ' +

					str(SimTime) + '\n')

			else:

				JobOut.write('mpirun -np ' + str(CpuNum) + ' mdrun -cpi '+ Filename + '.cpt -s ' + Filename + '.tpr -deffnm  '+ Filename + ' -maxh ' +

					str(SimTime) + '\n')

			

			JobOut.write('\n')

	

			# first N-1 jobscripts check if runtime is less then 15 sec, if not, submit next script 

			if k < N-1: # note that we count from 0 to N-1, so the 50th jobscript has i = 49

				JobOut.write( 'RUNTIME=$(($(date +%s)-$STARTTIME))\n' + '\n'+ 'echo \"the job took $RUNTIME seconds...\"\n' + '\n' + 'if [[ $RUNTIME -lt 10 ]]; then\n' +

							'   echo "job took less than 10 seconds to run, aborting."\n' + '   exit\n' + 'else\n' + '   echo "everything fine..."\n' + '   qsub ' + Scriptname + str(k+1) + '\n') 

				JobOut.write('\n')

		

			# after the last run, we also want to backup the simulation files

			elif k == N-1:

				JobOut.write('cp -r ' + SimDir + ' ' + BackupDir + '\n')

	JobOut.close()




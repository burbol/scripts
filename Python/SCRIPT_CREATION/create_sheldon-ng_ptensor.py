#!/usr/bin/env python

# !!IMPORTANT!! If files don't need to be copied, last PRINTED line should be a comment (it's a "cp" command to copy the whole simulations folder)

import os
import math
import fnmatch # for filtering the file list to get the top,gro,mdp,ndx files.

Partition = 'main' # use 'test' for testing and 'main' otherwise 
Memory = 1024
N = 3  # NUMBER OF SCRIPTS TO SUBMIT = total time running will be this number times the Walltime
Walltime = '48:00:00'  # Walltime in format 'h:mm:ss'
SimTime = 48   # maximum run length of simulation in hours, it can be also a fraction of an hour
NodesNum = 2
Cpus = 8  # ratio cpus/node!!! # soroban = 12 cpus/node ; sheldon/sheldon-ng = 8 cpus/node
CpuNum = NodesNum*Cpus

Email = 'laila.e@fu-berlin.de'

pc = [0,11,17,33,50,66,21,25,5]

for i in pc:

	Dir = './sheldon-ng/ptensor/w'+str(i) # directory to save script 
	Filename = 'sam'+str(i)+'_water_ptensor'  # name of files
	Scriptname = 'w'+str(i)+'_ptensor'

	Startfile=Filename
	Minifile='Mini_'+str(Filename)
	NVTfile = 'NVT_'+str(Filename)
	NPTfilePR1 = 'NPT_PR1_'+str(Filename)
	NPTfilePR2 = 'NPT_PR2_'+str(Filename)
	NPTfile = 'NPT_'+str(Filename)
	topfile= str(Filename)+'.top'
	Minimdp = 'Mini.mdp'
	NVTmdp = 'NVT.mdp'
	NPTmdpPR1 = 'NPT.mdp'
	NPTmdpPR2 = 'NPT.mdp'
	NPTmdp = 'NPT.mdp'
	
	#os.system("mkdir "+Dir)
	for k in range(N):

		Jobname = str(Scriptname)+ str(k)  # Name of job

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
		#JobOut.write('#SBATCH --ntasks=' + str(CpuNum)+'\n')
		JobOut.write('#SBATCH --tasks-per-node=' + str(Cpus)+'\n')
		JobOut.write('#SBATCH --nodes=' + str(NodesNum) +'\n')			
		#JobOut.write('#SBATCH --exclusive')  # we want to run on only 1 node
		JobOut.write('\n')
		JobOut.write('#SBATCH --time=' + Walltime  + '\n')
		JobOut.write('\n')
		JobOut.write('module load slurm \n')
		JobOut.write('\n')			
		JobOut.write('module load gromacs/single/openmpi1.4.5/4.6.5\n')
		JobOut.write('\n')

		if k == 0:
			JobOut.write('grompp -f '+ str(Minimdp) +' -c ' + Startfile+ '.gro -p '+str(topfile)+' -o '+str(Minifile)+'.tpr -maxwarn 1\n')
			JobOut.write('mpirun -np '+ str(CpuNum)+' mdrun -deffnm '+str(Minifile)+' -v -testverlet\n')

			JobOut.write('grompp -f '+ str(NVTmdp) + ' -c '+str(Minifile)+'.gro -p '+str(topfile)+' -o '+str(NVTfile)+'.tpr -maxwarn 1\n')
			JobOut.write('mpirun -np '+ str(CpuNum)+' mdrun -deffnm '+str(NVTfile)+' -v -testverlet\n')
			
			JobOut.write('grompp -f '+ str(NPTmdpPR1) + ' -c '+str(NVTfile)+'.gro -p '+str(topfile)+' -o '+str(NPTfilePR1)+'.tpr -maxwarn 1\n')
			JobOut.write('mpirun -np '+ str(CpuNum)+' mdrun -deffnm '+str(NPTfilePR1)+' -v -testverlet\n')
			
			JobOut.write('grompp -f '+ str(NPTmdpPR2) + ' -c '+str(NPTfilePR1)+'.gro -p '+str(topfile)+' -o '+str(NPTfilePR2)+'.tpr -maxwarn 1\n')
			JobOut.write('mpirun -np '+ str(CpuNum)+' mdrun -deffnm '+str(NPTfilePR2)+' -v -testverlet\n')

			JobOut.write('grompp -f '+ str(NPTmdp) + ' -c '+str(NPTfilePR2)+'.gro -p '+str(topfile)+' -o '+str(NPTfile)+'.tpr -maxwarn 1\n')
			JobOut.write('mpirun -np '+ str(CpuNum)+' mdrun -deffnm '+str(NPTfile)+' -maxh ' +str(SimTime) + '-v -testverlet\n')
		else:
			JobOut.write('mpirun -np '+ str(CpuNum)+' mdrun -cpi '+ str(NPTfile) + '.cpt -s ' + str(NPTfile) + '.tpr -deffnm '+str(NPTfile)+' -maxh ' +str(SimTime) + '-v -testverlet\n')	
			JobOut.write('\n')
		if k < N-1:
			JobOut.write('sbatch ' + Scriptname +'_'+ str(k+1) + '\n') 
	JobOut.close()
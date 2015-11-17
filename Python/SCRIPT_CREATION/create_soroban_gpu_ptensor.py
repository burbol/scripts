#!/usr/bin/env python

# !!IMPORTANT!! When not running tests, last PRINTED line should NOT be a comment (it's a "cp" command to copy the whole simulations folder)

import os
import math
import fnmatch # for filtering the file list to get the top,gro,mdp,ndx files.

Email = 'laila.e@fu-berlin.de'

Memory = 1024
Partition = 'gpu' # use 'test' for testing and 'main' otherwise 

N = 3  # NUMBER OF SCRIPTS TO SUBMIT = total time running will be this number times the Walltime
Walltime = '72:00:00'  # Walltime in format 'h:mm:ss'
SimTime = 72   # maximum run length of simulation in hours, it can be also a fraction of an hour
NodesNum = 1
Cpus = 4  # ratio cpus/node!!! # soroban = 12 cpus/node ; sheldon/sheldon-ng = 8 cpus/node
CpuNum = NodesNum*Cpus

BackupDir = '/home/eixeres/simulations/finished/'   # directory to copy (backup) simulation output

pc = [0,11,17,33,50,66,21,25,5]

for i in pc:

	#Dir = './soroban/ptensor/w'+str(i) # directory to save script 
	Dir = '/Users/burbol/Downloads/small_sams2/w'+str(i)+'/soroban'
	SimDir = '/home/eixeres/files_for_laila/ptensor/w'+str(i) #directory from where to run simulation
	Filename = 'sam'+str(i)+'_water_ptensor'  # name of files
	Scriptname = 's_w'+str(i)+'_ptensor'

	Startfile=Filename
	Minifile='Mini_'+str(Filename)
	NVTfile = 'NVT_'+str(Filename)
	NPTfilePR1 = 'NPT_PR1_'+str(Filename)
	NPTfilePR2 = 'NPT_PR2_'+str(Filename)
	NPTfile = 'NPT_'+str(Filename)
	topfile= str(Filename)+'.top'
	Minimdp = 'Mini.mdp'
	NVTmdp = 'NVT.mdp'
	NPTmdpPR1 = 'NPT_PR1.mdp'
	NPTmdpPR2 = 'NPT_PR2.mdp'
	NPTmdp = 'NPT.mdp'
	
	#os.system("mkdir " + Dir)
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
		JobOut.write('#SBATCH --ntasks=' + str(CpuNum)+'\n')
		JobOut.write('#SBATCH --nodes=' + str(NodesNum) +'\n')			
		#JobOut.write('#SBATCH --exclusive\n')
		JobOut.write('\n')
		JobOut.write('#SBATCH --time=' + Walltime  + '\n')
		JobOut.write('\n')
		JobOut.write('module load slurm \n')
		JobOut.write('\n')
		JobOut.write('module unload gromacs\n')
		# testing these modules
		JobOut.write('module load openmpi/intel/64/1.6.5\n')             
		JobOut.write('module load intel/mkl/64/11.0/2013.5.192\n')    
		JobOut.write('module load slurm/2.4.5\n')                    
		JobOut.write('module load intel/compiler/64/13.1/2013.5.192\n')  
		JobOut.write('module load intel-mpi/64/4.0.3/008\n')
		## testing ended
		#JobOut.write('module load gromacs/openmpi/gcc/64/4.5.4 \n')
		JobOut.write('module load gromacs/intelmpi/intel/64/4.6 \n')
		#JobOut.write('module load openmpi/intel/64/1.6.5 \n')
		JobOut.write('\n')
		JobOut.write('STARTTIME=$(date +%s)\n' + '\n'+ '#use sleep or testing... \n')
		JobOut.write('\n')
		JobOut.write('cd ' + SimDir + '\n')
		JobOut.write('\n')
		# first jobscript changes "mdp" options with tpbconv to extend simulation
		if k == 0:
			JobOut.write('grompp -f '+ str(Minimdp) +' -c ' + Startfile+ '.gro -p '+str(topfile)+' -o '+str(Minifile)+'.tpr -maxwarn 1\n')
			JobOut.write('mpirun -np '+ str(CpuNum)+' mdrun -deffnm '+str(Minifile)+' \n')

			JobOut.write('grompp -f '+ str(NVTmdp) + ' -c '+str(Minifile)+'.gro -p '+str(topfile)+' -o '+str(NVTfile)+'.tpr -maxwarn 1\n')
			JobOut.write('mpirun -np '+ str(CpuNum)+' mdrun -deffnm '+str(NVTfile)+' \n')
		
			JobOut.write('grompp -f '+ str(NPTmdpPR1) + ' -c '+str(NVTfile)+'.gro -p '+str(topfile)+' -o '+str(NPTfilePR1)+'.tpr -maxwarn 1\n')
			JobOut.write('mpirun -np '+ str(CpuNum)+' mdrun -deffnm '+str(NPTfilePR1)+' \n')
		
			JobOut.write('grompp -f '+ str(NPTmdpPR2) + ' -c '+str(NPTfilePR1)+'.gro -p '+str(topfile)+' -o '+str(NPTfilePR2)+'.tpr -maxwarn 1\n')
			JobOut.write('mpirun -np '+ str(CpuNum)+' mdrun -deffnm '+str(NPTfilePR2)+' \n')

			JobOut.write('grompp -f '+ str(NPTmdp) + ' -c '+str(NPTfilePR2)+'.gro -p '+str(topfile)+' -o '+str(NPTfile)+'.tpr -maxwarn 1\n')
			JobOut.write('mpirun -np '+ str(CpuNum)+' mdrun -deffnm '+str(NPTfile)+' -maxh ' +str(SimTime) + '\n')
		else:
			JobOut.write('mpirun -np '+ str(CpuNum)+' mdrun -cpi '+ str(NPTfile) + '.cpt -s ' + str(NPTfile) + '.tpr -deffnm '+str(NPTfile)+' -maxh ' +str(SimTime) + '\n')	
			JobOut.write('\n')
		
		JobOut.write('\n')
	
		# first N-1 jobscripts check if runtime is less then 15 sec, if not, submit next script 
		if k < N-1: # note that we count from 0 to N-1, so the 50th jobscript has i = 49
			JobOut.write( 'RUNTIME=$(($(date +%s)-$STARTTIME))\n' + '\n'+ 'echo \"the job took $RUNTIME seconds...\"\n' + '\n' + 'if [[ $RUNTIME -lt 10 ]]; then\n' +
						'   echo "job took less than 10 seconds to run, aborting."\n' + '   exit\n' + 'else\n' + '   echo "everything fine..."\n' + '   sbatch ' + Scriptname + str(k+1) + '\n') 
			JobOut.write('   fi\n')
			JobOut.write('\n')
		# after the last run, we also want to backup the simulation files
		elif k == N-1:
			JobOut.write('cp -r ' + SimDir + ' ' + BackupDir + '\n')
			print("finished"+Scriptname) 

	JobOut.close()




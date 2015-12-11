#!/usr/bin/env python

# !!IMPORTANT!! If files don't need to be copied, last PRINTED line should be a comment (it's a "cp" command to copy the whole simulations folder)

import os
import math
import fnmatch # for filtering the file list to get the top,gro,mdp,ndx files.

N = 3  # NUMBER OF SCRIPTS TO SUBMIT = total time running will be this number times the Walltime
Walltime = '72:00:00'  # Walltime in format 'h:mm:ss'
SimTime = 72   # maximum run length of simulation in hours, it can be also a fraction of an hour
NodesNum = 1
Cpus = 4  # ratio cpus/node!!! # soroban = 12 cpus/node ; sheldon/sheldon-ng = 8 cpus/node
CpuNum = NodesNum*Cpus

Email = 'laila.e@fu-berlin.de'

#pc = [17,21,33,50,66]
#pc = [5,11]
#pc = [0]
pc = [25]

for i in pc:
	#old folders:
	#Dir = './sheldon/ptensor/w'+str(i) # directory to save script 
	#Dir = '/Users/burbol/Downloads/small_sams2/w'+str(i)+'/sheldon'
	
	#new folders:
	Dir = '/Users/burbol/Downloads/small_sams2/WaterSamSeparated/w'+str(i)
	#Dir = '/Users/burbol/Downloads/small_sams2/WaterSamSeparated'  # 0% is in the upper folder WaterSamSeparated
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
	#NVTmdp = 'NVT.mdp'
	#NPTmdpPR1 = 'NPT_PR1.mdp'
	#NPTmdpPR2 = 'NPT_PR2.mdp'
	NVTmdp = 'NVT_100ps.mdp'
	NPTmdpPR1 = 'NPT_double_PR1.mdp'
	NPTmdpPR2 = 'NPT_double_PR2.mdp'
	NPTmdp = 'NPT.mdp'
	
	#os.system("mkdir "+Dir)
	for k in range(N):

		Jobname = str(Scriptname)+ str(k)  # Name of job

		JobOut = open(Dir + '/' + Scriptname + str(k),'w')
		print("openning"+Dir + '/' + Scriptname + str(k))

		JobOut.write('#!/bin/bash\n')

		JobOut.write('#PBS -N ' + Jobname + '\n')

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
		
		JobOut.write('\n')

	

		# first jobscript changes "mdp" options with tpbconv to extend simulation

		if k == 0:
			JobOut.write('grompp -f '+ str(Minimdp) +' -c ' + Startfile+ '.gro -p '+str(topfile)+' -o '+str(Minifile)+'.tpr -maxwarn 1\n')
			JobOut.write('mdrun -nt '+ str(CpuNum)+' -deffnm '+str(Minifile)+' \n')

			JobOut.write('grompp -f '+ str(NVTmdp) + ' -c '+str(Minifile)+'.gro -p '+str(topfile)+' -o '+str(NVTfile)+'.tpr -maxwarn 1\n')
			JobOut.write('mdrun -nt '+ str(CpuNum)+' -deffnm '+str(NVTfile)+' \n')
			
			JobOut.write('grompp -f '+ str(NPTmdpPR1) + ' -c '+str(NVTfile)+'.gro -p '+str(topfile)+' -o '+str(NPTfilePR1)+'.tpr -maxwarn 1\n')
			JobOut.write('mdrun -nt '+ str(CpuNum)+' -deffnm '+str(NPTfilePR1)+' \n')
			
			JobOut.write('grompp -f '+ str(NPTmdpPR2) + ' -c '+str(NPTfilePR1)+'.gro -p '+str(topfile)+' -o '+str(NPTfilePR2)+'.tpr -maxwarn 1\n')
			JobOut.write('mdrun -nt '+ str(CpuNum)+' -deffnm '+str(NPTfilePR2)+' \n')

			JobOut.write('grompp -f '+ str(NPTmdp) + ' -c '+str(NPTfilePR2)+'.gro -p '+str(topfile)+' -o '+str(NPTfile)+'.tpr -maxwarn 1\n')
			JobOut.write('mdrun -nt '+ str(CpuNum)+' -deffnm '+str(NPTfile)+' -maxh ' +str(SimTime) + '\n')
		else:
			JobOut.write('mdrun -nt '+ str(CpuNum)+' -cpi '+ str(NPTfile) + '.cpt -s ' + str(NPTfile) + '.tpr -deffnm '+str(NPTfile)+' -maxh ' +str(SimTime) + '\n')	
			JobOut.write('\n')
		if k < N-1:
			JobOut.write('qsub ' + Scriptname + str(k+1) + '\n')
		#elif k == N-1:
			#JobOut.write('qsub ' + Scriptname + str(0) + '\n')
		JobOut.close()

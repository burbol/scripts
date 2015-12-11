#!/usr/bin/env python

# !!IMPORTANT!! If files don't need to be copied, last PRINTED line should be a comment (it's a "cp" command to copy the whole simulations folder)

#/sw/chem/gromacs/4.6.7/mpp1/GNU/bin  gromacs4 folder
#/sw/chem/gromacs/5.0.4/mpp1/GNU/bin  gromacs5 folder
# working command: aprun -n 2 mdrun -deffnm NVT_sam21_water216 -maxh 0.1
#msub -I -l feature=mpp1 -l walltime=00:20:00 -l nodes=4:ppn=24


import os
import math
import fnmatch # for filtering the file list to get the top,gro,mdp,ndx files.

Dir = '/Users/burbol2/Dropbox/scripts/Python/SCRIPT_CREATION/HLRN/' # directory to save scripts in folders 

#N = 40  # NUMBER OF SCRIPTS TO SUBMIT = total time running will be this number times the Walltime
Walltime = '12:00:00'  # Walltime in format 'h:mm:ss'
SimTime = 12   # maximum run length of simulation in hours, it can be also a fraction of an hour

#SimLen=Simns*100000  # total simulation time length in ps

Email = 'laila.e@fu-berlin.de'

#NodesNum = 6  # hlrn = 24 ppn & min nodes=4 !!; soroban = 12 cpus/node ; sheldon/sheldon-ng = 8 cpus/node
CpuNum = 24  # processors per node
NodeType = 'mpp1' #mpp1, smp1, data, prepost1, test, ccm1 
#=> for gromacs5 use mpp1, which is only in Berlin!
#(check https://www.hlrn.de/home/view/System3/BatchSystem#Requesting_Compute_Resources)

pc = [21, 25]
molec = [1000, 2000, 3000, 4000, 5000, 6500, 7000, 8000, 9000, 10000]

#Simns = {216:60, 1000:80, 2000:80, 3000:80, 4000:80, 5000:100, 6500:100, 7000:100, 8000:120, 9000:120, 10000:120}  # total simulation time length in ns
Simns = {1000:40, 2000:40, 3000:40, 4000:40, 5000:40, 6500:40, 7000:40, 8000:40, 9000:40, 10000:40}  # total simulation time length in ns
#NodesNum = {216:4, 1000:4, 2000:4, 3000:4, 4000:4, 5000:6, 6500:6, 7000:6, 8000:8, 9000:8, 10000:8}
NodesNum = {1000:4, 2000:4, 3000:4, 4000:4, 5000:4, 6500:4, 7000:4, 8000:4, 9000:4, 10000:4} 



#pc = [21]  # For testing!!!!
#molec = [216]  # For testing!!!!

os.chdir(Dir)
for i in pc:
	for j in molec: 

		SimDir = '/gfs1/work/beclaila/s'+str(i)+'_w'+str(j)   #directory from where to run simulation
		#BackupDir = '/qfs1/perm/beclaila/NewVersion1/s'+str(i)+'_w'+str(j)   # directory to copy (backup) simulation output
		BackupDir = '/home/b/beclaila/s'+str(i)+'_w'+str(j) 
		Startfile='sam'+str(i)+'_water'+str(j)
		Minifile='Mini_sam'+str(i)+'_water'+str(j)
		NVTfile = 'NVT_sam'+str(i)+'_water'+str(j)
		topfile= str(i)+'pc_'+str(j)+'_cuda.top'
		Minimdp = 'Mini_cuda.mdp'
		NVTmdp = 'NVT_' + str(Simns[j]) +'ns_cuda.mdp'
		
		Scriptname = 's'+str(i)+'_w'+str(j)+'_cuda' 
		Jobname = 's'+str(i)+'_w'+str(j) # Name of job
		Indexfile = 'index'+str(i)+'_'+str(j)
		
		os.system('mkdir '+Jobname )
		
		N = Simns[j]  # NUMBER OF SCRIPTS TO SUBMIT: 1 scripts each ns 
		Totalcpus=NodesNum[j]*CpuNum
		#N=3 # For testing!!!!

		for k in range(N):

			f1 = open(Dir + '/' + Jobname +'/'+ Scriptname + '_' + str(k),'w+')

			f1.write('#!/bin/bash\n')

			f1.write('#PBS -N ' + Jobname  +'_'+ str(k) +  '\n')

			f1.write('#PBS -o ' + Scriptname + '.stdout\n')

			f1.write('#PBS -e ' + Scriptname + '.stderr\n')

			f1.write('#PBS -M ' + Email + '\n')

			f1.write('\n')
			
			f1.write('#PBS -l feature='+ NodeType +'\n')

			f1.write('#PBS -l nodes=' + str(NodesNum[j]) + ':ppn=' + str(CpuNum)+'\n')

			f1.write('#PBS -l walltime=' + Walltime  + '\n')

			f1.write('\n')

			f1.write('cd $PBS_O_WORKDIR\n')

			f1.write('\n')

			#f1.write('module load gromacs/4.6.7\n')
			
			f1.write('module unload PrgEnv-cray/5.2.40\n')
			f1.write('module load PrgEnv-gnu\n')
			f1.write('module load gromacs/5.0.4\n')

			f1.write('\n')

			f1.write('STARTTIME=$(date +%s)\n' + '\n')

			f1.write('\n')

			if k == 0:

				f1.write('aprun -n 1 grompp -f ' + Minimdp + ' -c ' + Startfile + '.gro -p ' + topfile + ' -o '+ Minifile +'.tpr -maxwarn 1\n')

				f1.write('aprun -n ' + str(Totalcpus) + ' mdrun -deffnm ' + Minifile + ' -maxh ' + str(SimTime) + '\n')
				
				f1.write('aprun -n 1 grompp -f ' + NVTmdp + ' -c ' + Minifile + '.gro -p ' + topfile + ' -o '+ NVTfile +'.tpr -maxwarn 1\n')
				
				f1.write('aprun -n ' + str(Totalcpus) + ' mdrun -deffnm ' + NVTfile +  ' -tunepme -v -testverlet -cpo '+ NVTfile + '.cpt  -maxh ' +str(SimTime) + '\n')
			else:
				
				f1.write('aprun -n ' + str(Totalcpus) + ' mdrun -deffnm '+ NVTfile + ' -s ' + NVTfile + '.tpr -tunepme -v -testverlet -cpi  '+ NVTfile + '.cpt -maxh ' + str(SimTime) + '\n')

			f1.write('\n')
			
			if k < N-1:
				f1.write( 'RUNTIME=$(($(date +%s)-$STARTTIME))\n' + '\n'+ 'echo \"the job took $RUNTIME seconds...\"\n' +  \
							'\n' + 'if [[ $RUNTIME -lt 10 ]]; then\n' +'   echo "job took less than 10 seconds to run, aborting."\n' +
							'   exit\n' + 'else\n' +'   echo "everything fine..."\n'+ '   msub ' + Scriptname + '_' + str(k+1) + '\n')

				f1.write('fi\n')
				f1.write('\n')

			# after the last run, we also want to backup the simulation files
			elif k == N-1:
				f1.write('cp -r ' + SimDir + ' ' + BackupDir + '\n')
		
				#f1.write('msub ' + Scriptname + '_' + str(k+1) + '\n')
				
			f1.close()

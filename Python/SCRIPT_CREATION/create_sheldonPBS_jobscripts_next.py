#!/usr/bin/env python

# !!IMPORTANT!! If files don't need to be copied, last PRINTED line should be a comment (it's a "cp" command to copy the whole simulations folder)

import os
import math
import fnmatch # for filtering the file list to get the top,gro,mdp,ndx files.

N = 5  # NUMBER OF SCRIPTS TO SUBMIT = total time running will be this number times the Walltime
Walltime = '300:00:00'  # Walltime in format 'h:mm:ss'
SimTime = 300   # maximum run length of simulation in hours, it can be also a fraction of an hour

Email = 'laila.e@fu-berlin.de'

Cpus = 12  # ratio cpus/node!!! # soroban = 12 cpus/node ; sheldon/sheldon-ng = 8 cpus/node

#next 3 lines only useful when runnig also grompp
#IndexFile = ' '
#FirstGro = 'sam5_water216'
#FirstMdp = 'Mini.mdp'

#pc = [0,5,17]
#molec = [1000, 6500,9000]

pc = [25]
molec = [1000, 2000, 3000, 4000, 5000, 6500, 7000, 8000, 9000, 10000]

Simns = {216:60, 1000:60, 2000:60, 3000:60, 4000:60, 5000:60, 6500:60, 7000:80, 8000:80, 9000:80, 10000:100}  # total simulation time length in ns
Nodes = {216:1, 1000:1, 2000:1, 3000:1, 4000:1, 5000:1, 6500:2, 7000:2, 8000:2, 9000:2, 10000:2} 

for i in pc:
	for j in molec: 

		SimDir = '/scratch/eixeres/Dec14_Last_Sims/s'+str(i)+'_w'+str(j)   #directory from where to run simulation
		
		Startfile='sam'+str(i)+'_water'+str(j)
		Minifile='Mini_sam'+str(i)+'_water'+str(j)
		NVTfile = 'NVT_sam'+str(i)+'_water'+str(j)
		topfile= str(i)+'pc_'+str(j)+'.top'
		Minimdp = 'Mini.mdp'
		NVTmdp = 'NVT_'+str(Simns[j])+'ns.mdp'
		#N = Simns[j]  # NUMBER OF SCRIPTS TO SUBMIT: 
		
	
	
		NodesNum = Nodes[j]
		CpuNum = NodesNum*Cpus
		SimLen = Simns[j]*1000  # total simulation time length in ps
		
		Dir = './sheldon/s'+str(i)+'_w'+str(j) # directory to save script 
		Filename = 'NVT_sam'+str(i)+'_water'+str(j)  # name of files
		Scriptname = 's'+str(i)+'_w'+str(j)  # name of script (when NOT testing, set equal to "Filename")

		os.system("mkdir "+ Dir)
		for k in range(2,N):

			Jobname = Scriptname +'_'+ str(k)  # Name of job

			JobOut = open(Dir + '/' + Jobname,'w')

			JobOut.write('#!/bin/bash\n')

			JobOut.write('#PBS -N ' + Jobname +  '\n')

			JobOut.write('#PBS -o ' + Scriptname + '.stdout\n')
	
			JobOut.write('#PBS -e ' + Scriptname + '.stderr\n')

			JobOut.write('#PBS -M ' + Email + '\n')

			JobOut.write('\n')

			JobOut.write('#PBS-l nodes=' + str(NodesNum) + ':ppn=' + str(Cpus)+'\n')

			JobOut.write('#PBS -l walltime=' + Walltime  + '\n')

			JobOut.write('\n')

			JobOut.write('cd $PBS_O_WORKDIR \n')

			JobOut.write('NPROCS=`wc -l < $PBS_NODEFILE`\n')

			JobOut.write('\n')
			
			#JobOut.write('source /net/opt/gromacs/openmpi1.4.5/4.5.4/bin/GMXRC\n')
			#JobOut.write('source /net/opt/gromacs/single/thread/4.6.5/bin/GMXRC\n')
			
			JobOut.write('\n')

	

			# first jobscript changes "mdp" options with tpbconv to extend simulation

			if k == 0:
			
				JobOut.write('/home/shavkat/GMX/bin/grompp -f ' + Minimdp + ' -c ' + Startfile + '.gro -p ' + topfile + ' -o '+ Minifile +'.tpr -maxwarn 1\n')
				JobOut.write('mpirun -np ' + str(CpuNum) + ' /home/shavkat/GMX/bin/mdrun_mpi -deffnm ' + Minifile + ' -maxh ' + str(SimTime) + ' -v -testverlet\n')				
				JobOut.write('/home/shavkat/GMX/bin/grompp -f ' + NVTmdp + ' -c ' + Minifile + '.gro -p ' + topfile + ' -o '+ NVTfile +'.tpr -maxwarn 1\n')								
				JobOut.write('\n')
			
			
				#JobOut.write('/home/shavkat/GMX/bin/tpbconv -s ' + Filename + '.tpr  -until ' + str(SimLen) + ' -o ' + Filename + '.tpr \n')
				JobOut.write('mpirun -np ' + str(CpuNum) + ' /home/shavkat/GMX/bin/mdrun_mpi -deffnm ' + Filename + ' -maxh ' +str(SimTime) + ' -v -testverlet\n')

                # OLDER commands:
				#JobOut.write('tpbconv -s ' + Filename + '.tpr  -until ' + str(SimLen) + ' -o ' + Filename + '.tpr \n')
				#JobOut.write('mpirun -np ' + str(CpuNum) + ' /home/shavkat/GMX/bin/g_tune_pme_mpi -s ' + Filename + '.tpr -cpo  '+ Filename + '.cpt -so ' + Filename + '.tpr -launch\n')
				##JobOut.write('mdrun -nt '+ str(CpuNum)+' -cpi '+ Filename + '.cpt -s ' + Filename + '.tpr -maxh ' +str(SimTime) + ' -v -testverlet\n')
				#JobOut.write('mdrun -nt '+ str(CpuNum)+' -cpi '+ Filename + '.cpt -s ' + Filename + '.tpr -maxh ' +str(SimTime) + ' -v -testverlet\n')

			else:

				#JobOut.write('mdrun -nt '+ str(CpuNum)+' -cpi '+ Filename + '.cpt -s ' + Filename + '.tpr -maxh ' +str(SimTime) + ' -v -testverlet\n')
				#print('mdrun -nt '+ str(CpuNum)+' -cpi '+ Filename + '.cpt -s ' + Filename + '.tpr -maxh ' +str(SimTime) + ' -v -testverlet\n')
				JobOut.write('mpirun -np ' + str(CpuNum) + ' /home/shavkat/GMX/bin/mdrun_mpi -cpi '+ Filename + '.cpt -s ' + Filename + '.tpr -deffnm ' + Filename + ' -maxh ' +str(SimTime) + ' -v -testverlet\n')

				JobOut.write('\n')
				JobOut.write( 'RUNTIME=$(($(date +%s)-$STARTTIME))\n' + '\n'+ 'echo \"the job took $RUNTIME seconds...\"\n' + '\n' + 'if [[ $RUNTIME -lt 10 ]]; then\n' +

							'   echo "job took less than 10 seconds to run, aborting."\n' + '   exit\n' + 'else\n' + '   echo "everything fine..."\n' + '   qsub ' + Scriptname + '_'+str(k+1) + '\n') 

				JobOut.write('   fi\n')
				
				JobOut.write('   fi\n')

				JobOut.write('\n')


			
			JobOut.write('\n')
	
			if k < N-1:
	
				JobOut.write('qsub ' + Scriptname +'_'+ str(k+1) + '\n')
	
		JobOut.close()
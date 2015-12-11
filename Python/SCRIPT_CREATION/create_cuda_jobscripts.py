#!/usr/bin/env python

# !!IMPORTANT!! When not running tests, last PRINTED line should NOT be a comment (it's a "cp" command to copy the whole simulations folder)

import os
import math
import fnmatch # for filtering the file list to get the top,gro,mdp,ndx files.


N = 5  # NUMBER OF SCRIPTS TO SUBMIT = total time running will be this number times the Walltime

Walltime = '48:00:00'  # Walltime in format 'h:mm:ss'
SimTime = 48   # maximum run length of simulation in hours, it can be also a fraction of an hour

Email = 'laila.e@fu-berlin.de'

Simns = {216:60, 1000:80, 2000:80, 3000:80, 4000:80, 5000:100, 6500:100, 7000:100, 8000:120, 9000:120, 10000:120}  # total simulation time length in ns

pc = [33,50]
molec = [6500, 7000, 8000, 9000, 10000]

for i in pc:
	for j in molec:
		
		Dir = '/Users/burbol/MEGAsync/scripts/Python/SCRIPT_CREATION/NewVersionBIG/s'+str(i)+'_w'+str(j) # directory to save script
		SimDir = '/home/leixeres/s'+str(i)+'_w'+str(j)   #directory from where to run simulation
		Filename = 'NVT_sam'+str(i)+'_water'+str(j)  # name of files
		Scriptname = 's'+str(i)+'_w'+str(j)+'_cuda'  # name of script (when NOT testing, set equal to "Filename")
		
		Startfile='sam'+str(i)+'_water'+str(j)
		Minifile='Mini_sam'+str(i)+'_water'+str(j)
		NVTfile = 'NVT_sam'+str(i)+'_water'+str(j)
		topfile= str(i)+'pc_'+str(j)+'_cuda.top'
		Minimdp = 'Mini_cuda.mdp'
		NVTmdp = 'NVT_'+str(Simns[j])+'ns_cuda.mdp'
		
		SimLen = Simns[j]*1000  # total simulation time length in p
	
		os.chdir(Dir)
		for k in range(N):

			Jobname = 's'+str(i)+'_w'+str(j)+ str(k)  # Name of job

			JobOut = open(Dir + '/' + Scriptname + str(k),'w')

			JobOut.write('#!/bin/bash\n')

			JobOut.write('\n')

			JobOut.write('#SBATCH --job-name=' + str(k) + Jobname + '\n')

			JobOut.write('\n')

			JobOut.write('#SBATCH --exclusive\n')

			JobOut.write('#SBATCH --time=' + Walltime  + '\n')

			JobOut.write('\n')
			JobOut.write('module add intel-studio-2013-sp1\n')
			JobOut.write('module add nvidia/cuda/5.5\n')
			JobOut.write('module add cuda/sdk/2.3\n')
			JobOut.write('module add comp/python/2.7.8-gcc\n')
			JobOut.write('module add nvidia/arrayfire/1.9\n')
			JobOut.write('module add nvidia/cuda/5.5\n')
			JobOut.write('module add mpi/openmpi/1.7.4/gcc_cuda\n')
			JobOut.write('module add mpi/openmpi/1.8.1/gcc\n')

			JobOut.write('\n')
			
			JobOut.write('export PKG_CONFIG_PATH="/home/leixeres/my_fftw/lib/pkgconfig"\n')
			JobOut.write('export CPPFLAGS="-I/home/leixeres/my_fftw/include"\n')
			JobOut.write('export LDFLAGS="-L/home/leixeres/my_fftw/lib"\n')
			JobOut.write('\n')
			
			JobOut.write('export PATH="/cluster/mpi/gcc/openmpi-1.7.4_gcc-4.8.2-cuda/bin":${PATH}\n')
			JobOut.write('export LD_LIBRARY_PATH="/cluster/mpi/gcc/openmpi-1.7.4_gcc-4.8.2-cuda/lib/":${LD_LIBRARY_PATH}\n')
			JobOut.write('export INCLUDE="/cluster/mpi/gcc/openmpi-1.7.4_gcc-4.8.2-cuda/include/":${INCLUDE}\n')
			JobOut.write('export MANPATH="/cluster/mpi/gcc/openmpi-1.7.4_gcc-4.8.2-cuda/share/man/":${MANPATH}\n')
			JobOut.write('\n')
			
			JobOut.write('export MPI_C_LIBRARIES="/cluster/mpi/gcc/openmpi-1.7.4_gcc-4.8.2-cuda/lib/"\n')
			JobOut.write('export MPI_CXX_LIBRARIES="/cluster/mpi/gcc/openmpi-1.7.4_gcc-4.8.2-cuda/lib/"\n')
			JobOut.write('export CUDA_NVCC_HOST_COMPILER="/cluster/mpi/gcc/openmpi-1.7.4_gcc-4.8.2-cuda/bin/mpiCC"\n')
			JobOut.write('\n')
			
			JobOut.write('source /home/leixeres/programs/gromacs/bin/GMXRC\n')
			JobOut.write('\n')

			JobOut.write('STARTTIME=$(date +%s)\n' + '\n')
			JobOut.write('\n')

			# first jobscript changes "mdp" options with tpbconv to extend simulation

			if k == 0:

				JobOut.write('/home/leixeres/programs/gromacs/bin/grompp_mpi -f ' + Minimdp + ' -c ' + Startfile + '.gro -p ' + topfile + ' -o '+ Minifile +'.tpr -maxwarn 1\n')

				JobOut.write('/home/leixeres/programs/gromacs/bin/mdrun_mpi -deffnm ' + Minifile + ' -maxh ' + str(SimTime) + ' -nb gpu_cpu\n')
				
				JobOut.write('/home/leixeres/programs/gromacs/bin/grompp_mpi -f ' + NVTmdp + ' -c ' + Minifile + '.gro -p ' + topfile + ' -o '+ NVTfile +'.tpr -maxwarn 1\n')

				JobOut.write('/home/leixeres/programs/gromacs/bin/mdrun_mpi -s -deffnm  '+ NVTfile + ' -maxh ' + str(SimTime) + ' -nb gpu_cpu\n')
				print('/home/leixeres/programs/gromacs/bin/mdrun_mpi -s -deffnm  '+ NVTfile + ' -maxh ' + str(SimTime) + ' -nb gpu_cpu\n')
				JobOut.write('\n')
			else:

				JobOut.write('/home/leixeres/programs/gromacs/bin/mdrun_mpi -cpi '+ Filename + '.cpt -s ' + NVTfile + '.tpr -deffnm  '+ NVTfile + ' -maxh ' +str(SimTime) + ' -nb gpu_cpu\n')

				JobOut.write('\n')

	

			# first N-1 jobscripts check if runtime is less then 15 sec, if not, submit next script

			if k < N-1: # note that we count from 0 to N-1, so the 50th jobscript has i = 49

				JobOut.write( 'RUNTIME=$(($(date +%s)-$STARTTIME))\n' + '\n'+ 'echo \"the job took $RUNTIME seconds...\"\n' + '\n' + 'if [[ $RUNTIME -lt 10 ]]; then\n' +

							'   echo "job took less than 10 seconds to run, aborting."\n' + '   exit\n' + 'else\n' + '   echo "everything fine..."\n' + '   sbatch -N 1 -p TESLAK20C -w "sandybridge08" ' + Scriptname + str(k+1) + '\n') 

				JobOut.write('   fi\n')

				JobOut.write('\n')

			JobOut.close()




#!/usr/bin/env python

import os
import math
import fnmatch # for filtering the file list to get the top,gro,mdp,ndx files.

#	- subdir with starting files (except for .gro!)
#	- subdir with .gro files
#	- list of names of simulation, simulation times, velocities, spring constants, nstxtcouts (or if it should be trr?!), number of times the simulation should be performed
# do:
#	- create subdirs according to the names + '_' + number
#	- copy starting files into dir
#	- copy one of the .gro files into dir
#	- adapt production.mdp (sim time, pullcode, nstxtcout, maybe even trr?)
#	- write a jobscript file
# at the end: create one bash script to submit.

StartingFilesDir='starting_files'
GroFilesDir='gro_files'

MolTypes='ii' # isoleucine
ListOfSimulationNames=['hola','quiero','trabajar','mas',
		'efectivo','v0.01','v0.05','v0.05']
ListOfSimulationLengths=[10110,1060,440,160,
			1060,1060,600,600] # in ns
ListOfTimesSimulationsShouldBePerformed=[1,1,1,1,
					1,1,1,1]

# this all goes into the changing of the mdp file
ListOfPullingVelocities=[0.005,0.01,0.1,1.0,
			0.01,0.01,0.05,0.05] # in m/s. 0 = no pull
#ListOfSpringConstants=[150,200,250,300] # in ??
ListOfSpringConstants=[150,150,200,250,
		      150,130,150,150]
ListOfPullGroups=['a_1-19','a_1-19','a_1-19','a_1-19',
		 'a_1-19','a_1-19','a_1-19','a_1-19']
ListOfRefGroups=['a_191-380','a_191-380','a_191-380','a_191-380',
		'a_191-380','a_191-380','a_191-380','a_191-380']
UseTrr=[False,False,False,False,
       False,False,False,False]
OutputFrequencies=[2000,2000,2000,20000,
		   250,200,50,50] # in simulation steps
OutputOnlyPeptide=[False,False,False,False,
		  True,True,True,True]

ProductionFilename = 'topol'
TrajFilename = 'traj'
EnergyFilename = 'ener'
LogFilename = 'md'

# jobscript stuff
JobscriptStyle='slurm'

NameOfJobscriptFile = 'jobscript'
NumberOfNodes = '1'
NumberOfProcessorsPerNode = '8'
ListOfWalltimes = ['48:00:00','48:00:00','48:00:00','48:00:00',
		  '48:00:00','48:00:00','48:00:00','48:00:00']
ListOfFilesizes = ['9000M','9000M','2000M',
		  '9000M','9000M','2000M','2000M','4000M']
Email = 'laila@physik.fu-berlin.de'
ListOfMaxH = ['48','48','48','48',
	     '48','48','48','48']

AddRestartPartForSlurm = True


def ActivatePull(Fname):
	File=open(Fname)
	Lines=File.readlines()
	File.close()
	out=open(Fname,'w')
	for Line in Lines:
		if Line[:5] == ";pull":
			out.write(Line[1:])
		else:
			out.write(Line)
	out.close()

def GetListOfFilesWithExtension(Directory,Extension):
	ListOfFiles = os.listdir(Directory)
	pattern = '*.' + str(Extension)
	return fnmatch.filter(ListOfFiles, pattern)

####
# 1. Go into StartingFilesDir and find the relevant files (.mdp, .top, .ndx)
TempFilelist = GetListOfFilesWithExtension(StartingFilesDir,'mdp')
if len(TempFilelist) != 1: sys.exit("RuntimeError:\tThere should be exactly one mdp file in " + StartingFilesDir + ", but I found " + str(len(TempFilelist)) + ".")
MdpFile = TempFilelist[0]

TempFilelist = GetListOfFilesWithExtension(StartingFilesDir,'top')
if len(TempFilelist) != 1: sys.exit("RuntimeError:\tThere should be exactly one top file in " + StartingFilesDir + ", but I found " + str(len(TempFilelist)) + ".")
TopFile = TempFilelist[0]

TempFilelist = GetListOfFilesWithExtension(StartingFilesDir,'ndx')
if len(TempFilelist) != 1: sys.exit("RuntimeError:\tThere should be exactly one ndx file in " + StartingFilesDir + ", but I found " + str(len(TempFilelist)) + ".")
NdxFile = TempFilelist[0]

# Get list of gro files in GroFilesDir
GroFilesList = GetListOfFilesWithExtension(GroFilesDir,'gro')
if len(GroFilesList) == 0: sys.exit("RuntimeError:\tNo gro files in " + str(GroFilesDir))
RunningIndex=0

OverallListOfSimulationDirs=[]
#
for i in range(len(ListOfSimulationNames)):
	CurSimName=ListOfSimulationNames[i]
	CurSimLen=ListOfSimulationLengths[i]
	CurOutFreq=OutputFrequencies[i]
# Create subdirs
	if ListOfTimesSimulationsShouldBePerformed[i] != 1:
		for j in range(ListOfTimesSimulationsShouldBePerformed[i]):
			CurDir=ListOfSimulationNames[i] + '_' + str(j)
			OverallListOfSimulationDirs.append(CurDir)
			CurMdp=CurDir+'/'+MdpFile
			os.system("mkdir -p " + CurDir)
			os.system("cp -r ./" + StartingFilesDir + "/* ./" + CurDir)
	else:
		CurDir=ListOfSimulationNames[i] + '_k' + str(ListOfSpringConstants[i])
		OverallListOfSimulationDirs.append(CurDir)
		CurMdp=CurDir+'/'+MdpFile
		os.system("mkdir -p " + CurDir)
		os.system("cp -r ./" + StartingFilesDir + "/* ./" + CurDir)
# Apply changes to mdp file
	# pull 
		if ListOfPullingVelocities[i] != 0:
			ActivatePull(CurMdp)
			File=open(CurMdp)
			Lines=File.readlines()
			File.close()
			out=open(CurMdp,'w')
			for Line in Lines:
				if Line[:10] == 'pull_rate1':
					out.write("pull_rate1      = " + str(ListOfPullingVelocities[i]/1000.0) + ' ; nm/ps = ' + str(ListOfPullingVelocities[i]) + ' m/s\n')
				elif Line[:7] == 'pull_k1':
					out.write("pull_k1         = " + str(ListOfSpringConstants[i]) + ' ; kJ mol^-1 nm^-2\n')
				elif Line[:11] == 'pull_group0':
					out.write("pull_group0     = " + str(ListOfRefGroups[i]) + '\n')
				elif Line[:11] == 'pull_group1':
					out.write("pull_group1     = " + str(ListOfPullGroups[i]) + '\n')
				else:
					out.write(Line)
			out.close()
	# simulation length and output
		File=open(CurMdp)
		Lines=File.readlines()
		File.close()
		out=open(CurMdp,'w')
		for Line in Lines:
			if Line[:2] == 'dt':
				dt=float(Line.split()[2]) # this is in ps
				out.write(Line)
			elif Line[:6] == 'nsteps':
				nsteps=int(math.ceil(ListOfSimulationLengths[i]*1000.0/dt))
				out.write("nsteps                   = " + str(nsteps) + ' ; ' + str(ListOfSimulationLengths[i]) + ' ns\n')
			elif Line[:7] == 'nstxout':
				out.write("nstxout                  = ")
				if UseTrr[i] == True:
					out.write(str(OutputFrequencies[i]) + ' ; ' + str(dt*OutputFrequencies[i]) + ' ps\n')
				else:
					out.write('0\n')
			elif Line[:9] == 'nstxtcout':
				out.write("nstxtcout                = ")
				if UseTrr[i] == True:
					out.write('0\n')
				else:
					out.write(str(OutputFrequencies[i]) + ' ; ' + str(dt*OutputFrequencies[i]) + ' ps\n')
			elif Line[:8] == 'xtc-grps':
				out.write("xtc-grps                 = ")
				if OutputOnlyPeptide[i]:
					out.write("Protein\n")
				else:
					out.write("\n")
			else:
				out.write(Line)
		out.close()
#	copy gro file
		CurGro = GroFilesList[RunningIndex%len(GroFilesList)]
		RunningIndex += 1
		os.system("cp ./" + GroFilesDir + '/' + CurGro + ' ./' + CurDir + '/' + CurGro)

## Create jobscript
		if JobscriptStyle == 'torque':
			JobOut = open('./' + CurDir + '/' + NameOfJobscriptFile,'w')
			JobOut.write('#!/bin/bash\n')
			JobOut.write('#PBS -N ' + MolTypes + '_' + CurDir + '\n')
			JobOut.write('#PBS -l nodes=' + NumberOfNodes + ':ppn=' + NumberOfProcessorsPerNode + '\n')
			JobOut.write('#PBS -l walltime=' + ListOfWalltimes[i] + '\n')
			JobOut.write('#PBS -l file=' + ListOfFilesizes[i] + '\n')
			JobOut.write('#PBS -m ea -M ' + Email + '\n')
			JobOut.write('\n')
			JobOut.write('cd $PBS_O_WORKDIR\n')
			JobOut.write('\n')
			JobOut.write('workdir=/local_scratch/$PBS_JOBID')
			JobOut.write('\n')
			JobOut.write('cd $PBS_O_WORKDIR\n')
			JobOut.write('cp -r $PBS_O_WORKDIR/* $workdir\n')
			JobOut.write('cd $workdir\n')
			JobOut.write('\n')
			JobOut.write('grompp -f production.mdp -n index.ndx -p ' + TopFile + ' -c ' + CurGro + 
				    ' -o ' + ListOfSimulationNames[i] + '.tpr\n')
			JobOut.write('\n')
			JobOut.write('mdrun -nt 8 -s ' + ListOfSimulationNames[i] + '.tpr -o ' + 
				    ListOfSimulationNames[i] + ' -x ' + ListOfSimulationNames[i] + ' -e ' +
				    ListOfSimulationNames[i] + ' -g ' + ListOfSimulationNames[i] + ' -px -pf -maxh ' +
				    ListOfMaxH[i] + '\n')
			JobOut.write('\n')
			JobOut.write('cp -r * $PBS_O_WORKDIR')
			JobOut.close()
		elif JobscriptStyle == 'slurm':
			JobOut = open('./' + CurDir + '/' + NameOfJobscriptFile,'w')
			JobOut.write('#!/bin/bash\n')
			JobOut.write('#SBATCH --job-name=' + MolTypes + '_' + CurDir + '\n')
			JobOut.write('#SBATCH --mail-type=end\n')
			JobOut.write('#SBATCH --mem=2048\n')
			JobOut.write('#SBATCH --time=' + ListOfWalltimes[i] + '\n')
			JobOut.write('#SBATCH --mail-user=' + Email + '\n')
			JobOut.write('#SBATCH --ntasks=1\n')
			JobOut.write('#SBATCH --cpus-per-task=' + NumberOfProcessorsPerNode + '\n')
			JobOut.write('\n')
			JobOut.write('grompp -f production.mdp -n index.ndx -p ' + TopFile + ' -c ' + CurGro + 
				    ' -o ' + ProductionFilename + '.tpr -maxwarn 1\n')
			JobOut.write('\n')
			JobOut.write('mdrun -nt 8 -s ' + ProductionFilename + '.tpr -o ' + 
				    TrajFilename + ' -x ' + TrajFilename + ' -e ' + 
				    EnergyFilename + ' -g ' + LogFilename + ' -px -pf -maxh ' +
				    ListOfMaxH[i] + '\n')
			JobOut.write('\n')
			if AddRestartPartForSlurm:
				JobOut.write('sbatch ' + NameOfJobscriptFile + '_continue')
			JobOut.close()
			if AddRestartPartForSlurm:
				JobOut = open('./' + CurDir + '/' + NameOfJobscriptFile + '_continue','w')
				JobOut.write('#!/bin/bash\n')
				JobOut.write('#SBATCH --job-name=' + MolTypes + '_' + CurDir + '\n')
				JobOut.write('#SBATCH --mail-type=end\n')
				JobOut.write('#SBATCH --mem=2048\n')
				JobOut.write('#SBATCH --time=' + ListOfWalltimes[i] + '\n')
				JobOut.write('#SBATCH --mail-user=' + Email + '\n')
				JobOut.write('#SBATCH --ntasks=1\n')
				JobOut.write('#SBATCH --cpus-per-task=' + NumberOfProcessorsPerNode + '\n')
				JobOut.write('\n')
				JobOut.write('#grompp -f production.mdp -n index.ndx -p ' + TopFile + ' -c ' + CurGro + 
					    ' -o ' + ProductionFilename + '.tpr -maxwarn 1\n')
				JobOut.write('\n')
				JobOut.write('mdrun -nt 8 -s ' + ProductionFilename + '.tpr -o ' + 
					    TrajFilename + ' -x ' + TrajFilename + ' -e ' + 
					    EnergyFilename + ' -g ' + LogFilename + ' -px -pf -maxh ' +
					    ListOfMaxH[i] + ' -cpi state.cpt\n')
				JobOut.write('\n')
				JobOut.write('FirstTest=false\n')
				JobOut.write('SecondTest=false\n')
				JobOut.write('\n')
				JobOut.write('# check if the file "confout.gro" is there\n')
				JobOut.write('if [ ! -f confout.gro ]; then\n')
				JobOut.write('  FirstTest=true\n')
				JobOut.write('  echo "confout.gro not found in current directory."\n')
				JobOut.write('fi\n')
				JobOut.write('\n')
				JobOut.write('# check if there are more than 50 "slurm*" files\n')
				JobOut.write('count=`ls -l | grep slurm -c`\n')
				JobOut.write('if [ $count -lt 50 ]; then\n')
				JobOut.write('  SecondTest=true\n')
				JobOut.write('  echo "less than 50 simulations have been started from this directory"\n')
				JobOut.write('fi\n')
				JobOut.write('\n')
				JobOut.write('if [ "$FirstTest" = true -a "$SecondTest" = true ]; then\n')
				JobOut.write('sbatch jobscript_continue\n')
				JobOut.write('fi\n')
				JobOut.write('\n')
				JobOut.close()

# create submit.bash
Out = open('submit.sh','w')
Out.write('#!/bin/bash\n\n')
if JobscriptStyle == 'torque':
	for Sim in OverallListOfSimulationDirs:
		Out.write('cd ' + Sim + '\n')
		Out.write('qsub ' + NameOfJobscriptFile + '\n')
		Out.write('cd ..\n\n')
elif JobscriptStyle == 'slurm':
	for Sim in OverallListOfSimulationDirs:
		Out.write('cd ' + Sim + '\n')
		Out.write('sbatch ' + NameOfJobscriptFile + '\n')
		Out.write('cd ..\n\n')
Out.close()
os.system("chmod 755 submit.sh")


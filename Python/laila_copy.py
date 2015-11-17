#!/usr/bin/env python

import numpy as np
import math
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import sys

Nx = 21 # number of gridpoints in x-direction
Ny = 22 # number of gridpoints in y-direction

# V = (Vx,Vy) first basis vector 
Vx = 0.123123
Vy = 0.3468346412

# W = (Wx,Wy) second basis vector
Wx = 0.91274816745
Wy = 0.1111

zPos = 0.5

# percentage to replace
PercentageToReplace = int(sys.argv[1]) #40

# create arrays that hold particle positions and types
xPos = np.zeros([Nx,Ny],dtype=float) # x positions of particles
yPos = np.zeros([Nx,Ny],dtype=float) # y positions of particles
pType = np.zeros([Nx,Ny],dtype=str) # particle type
# xPos[i,j] = x coordinate of the particle at position [i,j]
# 0 <= i <= Nx
# 0 <= j <= Ny


def MinDist(iR,jR,iRand,jRand):
  MinimalDistance = Nx*Ny
  for l in range(Nrep):
    CurDist = np.sqrt((iR[l]-iRand)**2+(jR[l]-jRand)**2)
    if CurDist < MinimalDistance:
      MinimalDistance = CurDist
  return MinimalDistance



# set standard particle type
for i in range(Nx):
  for j in range(Ny):
    pType[i,j] = 'C' 

# set particle positions
for i in range(Nx):
  for j in range(Ny):
    xPos[i,j] = i*Vx + j*Wx
    yPos[i,j] = i*Vy + j*Wy

# exchange PercentageToReplace particles
#    How many do we need to replace?
Nrep = int(math.ceil(0.01*PercentageToReplace*Nx*Ny)) # number of particles to replace
#    Get indices of particles to replace:
iRep = np.ones(Nrep,dtype=int)*(-2)
jRep = np.ones(Nrep,dtype=int)*(-2)
for k in range(Nrep):
  iRandom = np.random.randint(Nx)
  jRandom = np.random.randint(Ny)
  # draw again if particles have already been changed
  # or are next to a particle that has already been
  # changed:
  while (MinDist(iRep,jRep,iRandom,jRandom) <= 1.2): # if you replace the 1.2 by 1.5,
                # there will also be no diagonal neighbors (highest percentage for this
		# seems to be 21)
    iRandom = np.random.randint(Nx)
    jRandom = np.random.randint(Ny)
  pType[iRandom,jRandom] = 'O'
  iRep[k] = iRandom
  jRep[k] = jRandom
  

# save the particle types and positions
with open('output_' + str(PercentageToReplace) + '.txt','w') as f:
  for i in range(Nx):
    for j in range(Ny):
      f.write(pType[i,j] + '\t' +
	    str(xPos[i,j]) + '\t' + 
	    str(yPos[i,j]) + '\t' + 
	    str(zPos) + '\n')

# plot the particle types and positions
#   create array with 0s where there is O
#   and ones where there is C:
Z = np.zeros([Nx,Ny],dtype=int)
for i in range(Nx):
  for j in range(Ny):
    if pType[i,j] == 'C':
      Z[i,j] = 1
#   plot
#    first we plot the heatmap
fig, ax = plt.subplots()
ax.imshow(Z, cmap=plt.cm.winter, interpolation='nearest')
fig.savefig('output_' + str(PercentageToReplace) + '_heatmap.pdf',format='pdf')
#    then we plot the actual positions
SizeOfDots = 15
fig, ax = plt.subplots()
for i in range(Nx):
  for j in range(Ny):
    if pType[i,j] == 'O':
      ax.plot([xPos[i,j]],[yPos[i,j]],marker='o',markersize=SizeOfDots,color='blue')
    else:
      ax.plot([xPos[i,j]],[yPos[i,j]],marker='o',markersize=SizeOfDots,color='green')
ax.set_xlim(xPos[0,0]-0.5,xPos[-1,-1]+0.5)
ax.set_ylim(yPos[0,0]-0.5,yPos[-1,-1]+0.5)
fig.savefig('output_' + str(PercentageToReplace) + '_positions.pdf',format='pdf')


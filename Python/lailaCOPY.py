#!/usr/bin/env python

import numpy as np
import math
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import sys

# percentage to replace
#PercentageToReplace = int(sys.argv[1]) #40

Nx = 10 # number of gridpoints in x-direction (carbons)
Ny = 8 # number of gridpoints in y-direction (carbons)

Nxoxy = 5 # number of oxygens in x-direction 
Nyoxy = 4 # number of oxygens in y-direction 



a0 = 5.00
# V = (Vx,Vy) first basis vector 
Vx = a0*np.sin(np.pi/6)
Vy = a0*np.cos(np.pi/6)

# W = (Wx,Wy) second basis vector
Wx = a0*np.sin(np.pi/6)
Wy = -a0*np.cos(np.pi/6)

zPos = 21.680


# create arrays that hold particle positions and types
xPos = np.zeros([Nx,Ny],dtype=float) # x positions of particles
yPos = np.zeros([Nx,Ny],dtype=float) # y positions of particles
pType = np.zeros([Nx,Ny],dtype=str) # particle type
# xPos[i,j] = x coordinate of the particle at position [i,j]
# 0 <= i <= Nx
# 0 <= j <= Ny


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
#    Get indices of particles to replace:
xinterval = Nx/Nxoxy
yinterval = Ny/Nyoxy
for k in range(0,Nx,xinterval):
    for l in range(0,Ny,yinterval):
      pType[k,l] = 'O' 

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
plt.show()
#fig.savefig('output_' + str(PercentageToReplace) + '_positions.pdf',format='pdf')
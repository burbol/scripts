# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <codecell>

%pylab inline

# <codecell>

#!/usr/bin/env python

import numpy as np
import math
##import matplotlib as mpl
##mpl.use('Agg')
import matplotlib.pyplot as plt
import sys

# <codecell>

###################    SET THE PERCENTAGE MANUALLY     #################### 
####  FIRST, WE DETERMINE THE NUMBER OF CARBONS AND OXYGENS NEEDED  ####
percentage = 15
percentage = percentage/100.0
zPos=21.680
# Here we define the functions and variables that we will use to find the number of carbons and oxygens.
dummy, Nx, Ny, NoxyTotal, Exactpc = False, 0, 0, 0, 0
def goodratio(x):
    if (abs(round(x,2)-round(x)))<=0.01 and round(x**2)>=7:
        return True
    else:
        return False  
def searchratio(percentage):
    for n in range(2,100):
        for m in range(2,100):
            NoxyInt = sqrt((n*m*8.*percentage)/(1-percentage))
            if goodratio(NoxyInt):
                NoxyTotal = round(NoxyInt**2) #total number of oxygens
                Nx = m  # Number of carbons in x-direction = Nx
                Ny = n*8  # Number of carbons in y-direction = Ny
                Exactpc = round(NoxyInt**2,6)/(n*m*8 + round(NoxyInt**2,6)) # Exact percentage
                #return True
                return Nx, Ny, NoxyTotal, Exactpc
    print "Desired ratio not found. Check returned values!!"
    return Nx, Ny, NoxyTotal, Exactpc 

# This line is where the number of oxygens and carbons is set!!!
Nx, Ny, NoxyTotal, Exactpc = searchratio(percentage)   
Noxy = int(sqrt(NoxyTotal)) # number of oxygens in each direction --> Nox**2 = NoxTotal


#Here we check if the particle numbers found really produce the desired percentage
#print dummy, Nx, Ny, NoxyTotal, Exactpc
#print "total number of oxygens = ", NoxyTotal
#print "total number of carbons = ",Nx*Ny
#print "Number of carbons in x-direction = Nx = ", Nx
#print "Number of carbons in y-direction = Ny = ", Ny
#print "Exact percentage = ", Exactpc
#Roundpc = round(Exactpc,3) # Rounded percentage
#print "Rounded percentage = ", Roundpc

#If we need to set the number of oxygens and carbons MANUALLY we can uncomment the next 3 lines:
Nyoxy =int(percentage**2*150)
Ny = int(float(Nyoxy)/percentage) # number of carbons in y-direction
Nx = 5 # number of carbons in x-direction
Nxoxy = Nx
carbons = Nx*Ny
NoxyTotal = Nxoxy*Nyoxy # total number of oxygens
print "!!!!!!!!!!!Nyoxy=", Nyoxy, "carbons=", carbons, " NoxyTotal",NoxyTotal

##########################   HERE STARTS THE LATTICE CONSTRUNCTION   ##########################
# Arrays that will hold carbon particle positions
xPos = np.zeros([Nx,Ny],dtype=float) # x positions of particles
yPos = np.zeros([Nx,Ny],dtype=float) # y positions of particles

# Arrays that will hold oxygen particle positions
xPosOXY = np.zeros([Nxoxy,Nyoxy],dtype=float) # x positions of particles
yPosOXY = np.zeros([Nxoxy,Nyoxy],dtype=float) # y positions of particles


##########################   CARBONS: HEXAGONAL LATTICE   #################################

# WE SET THE CELL PARAMETER a0 
a0= 5.000

# We define a "unit" cell in the orthogonal standard basis. 
# V0, W0 span the unit cell (-->these "0's" are zeros NOT the "O-letter" from Oxygen!),
# which is a square with side length = a0 and area =a0**2!!!

# V0 = (V0x,V0y) first basis vector 
V0x = a0
V0y = 0
# W0 = (W0x,W0y) second basis vector
W0x = 0
W0y = a0

#We define the standard basis vectors V, W of the hexagonal lattice: we will need them to transform the unit cell 
#from the standard coordinate system to the coordinate system of the hexagonal lattice

# V = (Vx,Vy) first basis vector 
Vx = 1.0
Vy = 0.0
Vabs = sqrt(Vx**2+Vy**2)  # absolute value of basis vector should be =1
# W = (Wx,Wy) second basis vector
Wx = np.sin(np.pi/3.) #sin(60deg)
Wy = np.cos(np.pi/3.) #cos(60deg)
Wabs = sqrt((Wx**2)+(Wy**2)) # absolute value of basis vector should be =1

print 'Vabs=', Vabs, "Wabs=", Wabs # checking absolute value of basis vectors

# This are equivalent basis vectors V,W of the hexagonal lattice (uncomment if needed)
# V = (Vx,Vy) first basis vector 
#Vx = (3/2.)
#Vy = sqrt(3)/2.
# W = (Wx,Wy) second basis vector
#Wx = -(3/2.)
#Wy = sqrt(3)/2.

# Now we transform the "unit" cell from the standard coordinate system to the coord. system of the hex. lattice
# The transformed "unit" cell will be the area spanned by the vetors Vcarb, Wcarb.
# Vcarb and Wcarb are the vectors V0 and W0 muliplied by a transformations matrix created by the coordinates of Vx and Vy

# Vcarb = (Vxcarb,Vycarb) first "unit" cell vector 
Vxcarb = Vx*V0x + Wx*V0y
Vycarb = Vy*V0x + Wy*V0y
# Wcarb = (Wxcarb,Wycarb) second "unit" cell vector
Wxcarb = Vx*W0x + Wx*W0y
Wycarb = Vy*W0x + Wy*W0y

# We create the finite vector space spanned by the basis vectors of the "unit" cell 
# of our hexagonal lattice: Vcarb, Wcarb.
# This finite vector space is our hexagonal lattice (or grid) and has a total number of vectors = Nx*Ny.
# The positions of the carbon particles of the hexagonal lattice will be stored in 2D-arrays xPos and yPos.
for i in range(Nx):
  for j in range(Ny):
        xPos[i,j] = i*Vxcarb + j*Wxcarb
        yPos[i,j] = i*Vycarb + j*Wycarb
        
#print 'Checking if CARBON grid has correct size. These numbers should be equal:'
#print 'Side x = Max(xPos)+Vxcarb+Wxcarb = ',np.max(xPos)+Vxcarb+Wxcarb, ", Nx*Vxcarb + Ny*Wxcarb= ", float(Nx*Vxcarb + Ny*Wxcarb)
#print 'Side y = Max(yPos)+Vycarb+Wycarb = ',np.max(yPos)+Vycarb+Wycarb, ", Nx*Vycarb + Ny*Wycarb= ", float(Nx*Vycarb + Ny*Wycarb)
#print 'Carbon lattice total AREA = Side x * Side y  = ', (np.max(xPos)+Vxcarb+Wxcarb)*(np.max(yPos)+Vycarb+Wycarb), ' = ', \
#(Nx*Vxcarb + Ny*Wxcarb)*(Nx*Vycarb + Ny*Wycarb)
#print ' '

##########################   OXYGENS   #################################

# WE SET THE CELL PARAMETERS of the oxygen's "unit" cell. Here we need 2 values,
# because each side of the cell can have a different length. So we set a0x and a0y as the x- and y-coordinates 
a0x = (Nx*a0)/Nxoxy
a0y = (Ny*a0)/Nyoxy
print "(Nx*a0)/Nxoxy=",(Nx*a0)/Nxoxy
print "(Ny*a0)/Nyoxy=",(Ny*a0)/Nyoxy

# We define the basis vectors V0oxy, W0oxy of the oxygen's "unit" cell in the orthogonal unitary (standard) basis.
# It should define a rectangle with side lengths = a0x and a0y; and an area =a0x*a0y !!!

# V0oxy = (V0xoxy,V0yoxy) first basis vector 
V0xoxy = a0x
V0yoxy = 0
# W0oxy = (W0xoxy,W0yoxy) second basis vector
W0xoxy = 0
W0yoxy = a0y

# Now we transform the oxygen "unit" cell from the standard coordinate system to the coord. system of the hex. lattice
# The transformed "unit" cell will be the area spanned by the vetors Voxy, Woxy.
# Voxy and Woxy are the vectors V0oxy and W0oxy muliplied by a transformations matrix created before with the coordinates of Vx and Vy.

# Voxy = (Vxoxy,Vyoxy) first "unit" cell vector 
Vxoxy = Vx*V0xoxy + Wx*V0yoxy
Vyoxy = Vy*V0xoxy + Wy*V0yoxy
# Woxy = (Wxoxy,Wyoxy) second "unit" cell vector
Wxoxy = Vx*W0xoxy + Wx*W0yoxy
Wyoxy = Vy*W0xoxy + Wy*W0yoxy

# We create the finite vector space spanned by the basis vectors of the oxygen's "unit" cell Voxy and Woxy.
# This finite vector space is the oxygen's lattice (or grid) and has a total number of vectors = NoxyTotal = Nxoxy**Nyoxy.
# The positions of the oxygen particles will be stored in the 2D-arrays xPosOXY and yPosOXY.
for i in range(Nxoxy):
  for j in range(Nyoxy):
        xPosOXY[i,j] = i*Vxoxy + j*Wxoxy
        yPosOXY[i,j] = i*Vyoxy + j*Wyoxy

#print 'Checking if OXYGEN grid has correct size. These numbers should be equal:'
#print 'Side x = Max(xPosOXY)=',np.max(xPosOXY)+Vxoxy+Wxoxy, "Nxoxy*(Vxoxy+Wxoxy)= ", float(Nxoxy*(Vxoxy+Wxoxy)) 
#print 'Side y = Max(yPosOXY)=',np.max(yPosOXY)+Vyoxy+Wyoxy, "Nyoxy*(Vyoxy+Wyoxy)= ", float(Nyoxy*(Vyoxy+Wyoxy))
#print 'Oxygen lattice total AREA = Side x * Side y  = ', float((np.max(xPosOXY)+Vxoxy+Wxoxy)*(np.max(yPosOXY)+Vyoxy+Wyoxy)), ' = ', \
#float(Nxoxy*Nyoxy*(Vxoxy+Wxoxy)*(Vyoxy+Wyoxy))
##########################  PLOT OF THE RESULTING LATTICE  ##########################
SizeOfDotsC = 8
SizeOfDotsO = 3
fig, ax = plt.subplots()
ax.plot(xPos,yPos,marker='o',markersize=SizeOfDotsC,color='yellow')
ax.plot(xPosOXY,yPosOXY,marker='o',markersize=SizeOfDotsO,color='blue')
#ax.plot(xOPosFirst,yOPosFirst,marker='o',markersize=SizeOfDots,color='blue')
#ax.set_xlim(xPos[0,0]-0.5,xPos[-1,-1]+0.5)
#ax.set_xlim(120,200)
#ax.set_ylim(yPos[0,0]-0.5,yPos[-1,-1]+0.5)
##fig.savefig('output_' + str(PercentageToReplace) + '_positions.pdf',format='pdf')
plt.show()
#for i in range(Nx):
#    print xPos[i,0]
#for i in range(Noxy):
#    print xPosOXY[i,0]

# <codecell>

# save the particle types and positions
with open('output_' + str(PercentageToReplace) + '.txt','w') as f:
  for i in range(Nx):
    for j in range(Ny):
      f.write( '\t' + str(xPos[i,j]) + '\t' + str(yPos[i,j]) + '\t' + '\n')

# <codecell>

i=10
j=13

chainNum=80
totalpos=13
atompos=19
atomtype='C'
chaintype='SAM'
title = 'sam 0% OH-coverage'
atomnum=str(atomtype)+str(atompos)
occupancy=1.00
temp=0.00
print 'TITLE     ' + title
print 'REMARK    THIS IS A SIMULATION BOX'
print 'CRYST1  150.000  138.564  120.000  90.00  90.00  90.00 P 1           1'
print 'MODEL        1'
for i in range(2,5):
    #print '{0:6}{1:5}{2:4} {3} {4:5}  {5:10.3f} {6:7.3f} {7:7.3f}  {8:5} {9:5}'.format("ATOM ", l, 'C19', 'SAM',s, 
                                                                                    #round(xPos[i,j],3),round(yPos[i,j],3),
                                                                                    #zPos,'1.00','0.00')
    x= round(xPos[i,j],3)
    y= round(yPos[i,j],3)
    z= zPos
    print "%-6s%5d %4s%1s%3s %1s%4d%1s   %8.3f%8.3f%8.3f%6.2f%6.2f          %2s%2s"%("ATOM ", totalpos, atomnum, '',chaintype,'',chainNum,'', 
                                                                                    x,y,z,occupancy,temp,'','')	
    j = j+1
    totalpos=59516
    localpos=3

# <codecell>

#############################################  PDB FILE WRITTING    ###################################################

############## SAM & OAM CHAINS ##############
xLengthCC =
yLengthCC =
zLengthCC =

xLengthCH1 =
yLengthCH1 =
zLengthCH1 =

xLengthCH2 =
yLengthCH2 =
zLengthCH2 =

################# Head Groups  ##################

####  SAMs Head Groups  ####
xLengthHeadCH1 =
yLengthHeadCH1 =
zLengthHeadCH1 =

xLengthHeadCH2 =
yLengthHeadCH2 =
zLengthHeadCH2 =

xLengthHeadCH3 =
yLengthHeadCH3 =
zLengthHeadCH3 =

####  OAM Head Groups  ####

xLengthCO =
yLengthCO =
zLengthCO =

xLengthOH1 =
yLengthOH1 =
zLengthOH1 =

xLengthOH2 =
yLengthOH2 =
zLengthOH2 =

xLengthOH3 =
yLengthOH3 =
zLengthOH3 =


atomtype='C'
chaintype='SAM'
title = 'sam 0% OH-coverage'
atomnum=str(atomtype)+str(atompos)
occupancy=1.00
temp=0.00
print 'TITLE     ' + title
print 'REMARK    THIS IS A SIMULATION BOX'
print 'CRYST1  150.000  138.564  120.000  90.00  90.00  90.00 P 1           1'
print 'MODEL        1'

# We start all the counters
chainlength = 62
totalpos=1
atompos=1
chainNum = 1

# x-coord of the first atom of the file
xfirst = 0.000  
yfirst = 0.000
zfirst = 1.500


# First we write the top Head Group
i = 0
j = 0
atomtype='C'
xold= round(xPos[i,j],3)
yold= round(yPos[i,j],3)
zold= zPos
print "%-6s%5d %4s%1s%3s %1s%4d%1s   %8.3f%8.3f%8.3f%6.2f%6.2f          %2s%2s"%("ATOM ", totalpos, atomnum, '',chaintype,'',chainNum,'', xold, yold, zold, occupancy, temp,'','')
for i in range(Nx):
  for j in range(Ny): 
        # At the beginning of each chain we write a Head Group 
        if totalpos = (chainlength*(chainNum-1))+1:
            xnew = xold + xLengthHeadCH1
            ynew = yold + yLengthHeadCH1
            znew = zold + zLengthHeadCH1
            print "%-6s%5d %4s%1s%3s %1s%4d%1s   %8.3f%8.3f%8.3f%6.2f%6.2f          %2s%2s"%("ATOM ", totalpos, atomnum, '',chaintype,'',chainNum,'', x, y, z, occupancy, temp,'','')
        # At the end of each chain we also write a Head Group
        elif totalpos = (chainlength*chainNum)-3:
            xnew = xold + xLengthHeadCH1
            ynew = yold + yLengthHeadCH1
            znew = zold + zLengthHeadCH1
            print "%-6s%5d %4s%1s%3s %1s%4d%1s   %8.3f%8.3f%8.3f%6.2f%6.2f          %2s%2s"%("ATOM ", totalpos, atomnum, '',chaintype,'',chainNum,'', x, y, z, occupancy, temp,'','')
        else:
            xnew = xold + xLengthCH1
            ynew = yold + yLengthCH1
            znew = zold + zLengthCH1
        print "%-6s%5d %4s%1s%3s %1s%4d%1s   %8.3f%8.3f%8.3f%6.2f%6.2f          %2s%2s"%("ATOM ", totalpos, atomnum, '',chaintype,'',chainNum,'', xnew, ynew, znew, occupancy, temp,'','')
        
        xold=xnew
        yold=ynew
        zold=znew
        totalpos = totalpos + 1
        

# <codecell>

## Testing how the number of oxygens and carbons is searched"
percentage = 33
percentage = percentage/100.0
zPos=21.680
# Here we define the functions and variables that we will use to find the number of carbons and oxygens.
dummy, Nx, Ny, NoxyTotal, Exactpc = False, 0, 0, 0, 0
def goodratio(x):
    if (abs(round(x,2)-round(x)))<=0.01 and round(x**2)>=7:
        return True
    else:
        return False  
def searchratio(percentage):
    for n in range(2,100):
        for m in range(2,100):
            NoxyInt = sqrt((n*m*8.*percentage)/(1-percentage))
            if goodratio(NoxyInt):
                NoxyTotal = round(NoxyInt**2) #total number of oxygens
                Nx = m  # Number of carbons in x-direction = Nx
                Ny = n*8  # Number of carbons in y-direction = Ny
                Exactpc = round(NoxyInt**2,6)/(n*m*8 + round(NoxyInt**2,6)) # Exact percentage
                #return True
                return Nx, Ny, NoxyTotal, Exactpc
    print "Desired ratio not found. Check returned values!!"
    return Nx, Ny, NoxyTotal, Exactpc 

# <codecell>

# This line is where the number of oxygens and carbons is set!!!
Nx, Ny, NoxyTotal, Exactpc = searchratio(percentage)   
Noxy = int(sqrt(NoxyTotal)) # number of oxygens in each direction --> Nox**2 = NoxTotal


#Here we check if the particle numbers found really produce the desired percentage
print dummy, Nx, Ny, NoxyTotal, Exactpc
print "total number of oxygens = ", NoxyTotal
print "total number of carbons = ",Nx*Ny
print "Number of carbons in x-direction = Nx = ", Nx
print "Number of carbons in y-direction = Ny = ", Ny
print "Exact percentage = ", Exactpc
Roundpc = round(Exactpc,3) # Rounded percentage
print "Rounded percentage = ", Roundpc

# <codecell>



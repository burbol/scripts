#!/usr/bin/python

import numpy as np
import math
import matplotlib.pyplot as plt
import sys

zPos = 21.680
####  FIRST, WE DETERMINE THE NUMBER OF CARBONS AND OXYGENS NEEDED  ####
############    SET THE PERCENTAGE MANUALLY     ############# 
percentage = 0.015

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
print dummy, Nx, Ny, NoxyTotal, Exactpc
print "total number of oxygens = ", NoxyTotal
print "total number of carbons = ",Nx*Ny
print "Number of carbons in x-direction = Nx = ", Nx
print "Number of carbons in y-direction = Ny = ", Ny
print "Exact percentage = ", Exactpc
Roundpc = round(Exactpc,3) # Rounded percentage
print "Rounded percentage = ", Roundpc

#If we need to set the number of oxygens and carbons MANUALLY we can uncomment the next 3 lines:
#Nx = 24 # number of carbons in x-direction
#Ny = 16 # number of carbons in y-direction
#NoxyTotal = 9 # total number of oxygens

################   HERE STARTS THE LATTICE CONSTRUNCTION   ################

# Arrays that will hold carbon particle positions
xPos = np.zeros([Nx,Ny],dtype=float) # x positions of particles
yPos = np.zeros([Nx,Ny],dtype=float) # y positions of particles

# Arrays that will hold oxygen particle positions
xPosOXY = np.zeros([Noxy,Noxy],dtype=float) # x positions of particles
yPosOXY = np.zeros([Noxy,Noxy],dtype=float) # y positions of particles


##################   CARBONS: HEXAGONAL LATTICE   #########################

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
        
print 'Checking if CARBON grid has correct size. These numbers should be equal:'
print 'Side x = Max(xPos)+Vxcarb+Wxcarb = ',np.max(xPos)+Vxcarb+Wxcarb, ", Nx*Vxcarb + Ny*Wxcarb= ", float(Nx*Vxcarb + Ny*Wxcarb)
print 'Side y = Max(yPos)+Vycarb+Wycarb = ',np.max(yPos)+Vycarb+Wycarb, ", Nx*Vycarb + Ny*Wycarb= ", float(Nx*Vycarb + Ny*Wycarb)
print 'Carbon lattice total AREA = Side x * Side y  = ', (np.max(xPos)+Vxcarb+Wxcarb)*(np.max(yPos)+Vycarb+Wycarb), ' = ', \
(Nx*Vxcarb + Ny*Wxcarb)*(Nx*Vycarb + Ny*Wycarb)
print ' '

##########################   OXYGENS   #################################

# WE SET THE CELL PARAMETERS of the oxygen's "unit" cell. Here we need 2 values,
# because each side of the cell can have a different length. So we set a0x and a0y as the x- and y-coordinates 
a0x = (Nx*a0)/Noxy
a0y = (Ny*a0)/Noxy

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
# This finite vector space is the oxygen's lattice (or grid) and has a total number of vectors = NoxyTotal = Noxy**2.
# The positions of the oxygen particles will be stored in the 2D-arrays xPosOXY and yPosOXY.
for i in range(Noxy):
  for j in range(Noxy):
        xPosOXY[i,j] = i*Vxoxy + j*Wxoxy
        yPosOXY[i,j] = i*Vyoxy + j*Wyoxy

print 'Checking if OXYGEN grid has correct size. These numbers should be equal:'
print 'Side x = Max(xPosOXY)=',np.max(xPosOXY)+Vxoxy+Wxoxy, "Noxy*(Vxoxy+Wxoxy)= ", float(Noxy*(Vxoxy+Wxoxy)) 
print 'Side y = Max(yPosOXY)=',np.max(yPosOXY)+Vyoxy+Wyoxy, "Noxy*(Vyoxy+Wyoxy)= ", float(Noxy*(Vyoxy+Wyoxy))
print 'Oxygen lattice total AREA = Side x * Side y  = ', float((np.max(xPosOXY)+Vxoxy+Wxoxy)*(np.max(yPosOXY)+Vyoxy+Wyoxy)), ' = ', \
float(Noxy**2*(Vxoxy+Wxoxy)*(Vyoxy+Wyoxy))


########################## WRITE .PDB FILE ##########################
PercentageToReplace = int(percentage*1000)
i = 1
with open('start' + str(PercentageToReplace) + '.pdb','w') as f:
  for i in range(Nx):
    for j in range(Ny):
      f.write('ATOM      ' + str(i) +'C19' + '   ' +
        ' H1  OAM     1       ' + 
	    str(xPos[i,j]) + '   ' + 
	    str(yPos[i,j]) + '   ' + 
	    str(zPos) + '\n'+
	    '  1.00  0.00')
     print '{0:10} ==> {1:10d}'.format(i, j,round(xPos[i,j],3),round(yPos[i,j],3),zPos)	    
print 'ATOM      {}  C19  SAM     {}       0.100   0.000   0.510  1.00  0.00'    

###################  PLOT OF THE RESULTING LATTICE  ######################
SizeOfDots = 5
fig, ax = plt.subplots()
ax.plot(xPos,yPos,marker='o',markersize=SizeOfDots,color='yellow')
ax.plot(xPosOXY,yPosOXY,marker='o',markersize=SizeOfDots,color='blue')

#ax.set_xlim(xPos[0,0]-0.5,xPos[-1,-1]+0.5)
#ax.set_ylim(yPos[0,0]-0.5,yPos[-1,-1]+0.5)
fig.savefig('start' + str(PercentageToReplace) + '.pdf',format='pdf')
plt.show()
#for i in range(Nx):
#    print xPos[i,0]
#for i in range(Noxy):
#    print xPosOXY[i,0]
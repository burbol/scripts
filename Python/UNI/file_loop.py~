#!/usr/bin/python
import os 
import numpy as np

for i in [5]:
    for j in [2]:
        for k in range(0, 20, 2):
            l = k + 2
            filename = 'densmap_%dpc_%d000_%d_%d.dat'%(i,j,k,l, )
            foldername = 'densmaps_s%d_w%d000'%(i,j, )

os.chdir("/home/eixeres/Downloads/Berendsen")

os.chdir(foldername)

f = open('fileloop.txt', 'w')
print >> f, 'Filename:', filename  # or f.write('...\n')
print >> f, "Forldename:", foldername
print >> f, "Current directory:", os.getcwd()



print >> f, '{0}  {1}  {2}'.format('Case', 'Contact Angle', 'Base Radius')
z = np.zeros(11)
for x in range(1,11):
    z[x] = x**2 -5
    print >> f, '{0}  {1}  {2}'.format(z[x], z[x]**2, z[x]**3)

f.close()
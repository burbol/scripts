# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <codecell>

%pylab inline

# <codecell>

import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit

# <codecell>

%cd /Users/burbol/Desktop/Testing/files_for_laila/densmaps_berendsen

# <codecell>

%cd densmaps_s5_w1000/

# <codecell>

input = np.loadtxt("densmap_5pc_1000_12ns_14ns.dat")

# <codecell>

x = input[0]
slice = 90
y = input[slice]

# <codecell>

len(np.nonzero(y))

# <codecell>

x[np.nonzero(y)],  np.nonzero(y), len(np.nonzero(y))

# <codecell>

def checkzero(s):
        if len(np.nonzero(s)) > 1: 
            return True
        else: 
            return False

# <codecell>

zeroslice = 0
for slice in range(70,150):
    y = input[slice]
    if checkzero(y): 
       zeroslice = slice
       print zeroslice
       break
    else: continue

# <codecell>

print zeroslice
y = input[zeroslice]

# <codecell>

pyplot.plot(x,y)
pyplot.axis([0.05, 5, 0, 150])
pyplot.show()

# <codecell>

def func(r, ro, R, d):
    return (ro/2)*(1-np.tanh(2*(r-R)/d))

# <codecell>

popt, pcov = curve_fit(func, x, y,[100, 3, 1])

# <codecell>

import pylab as pl
pl.plot(x, func(x, *popt), 'r-')
pl.plot(x,y,'k.')

# <codecell>

z = y[0]
print z
shift = 4.7

# <codecell>

slicemin = 110
slicemax = 158
dims = slicemax - slicemin + 1
i=0
z2 = np.zeros(dims)
R = np.zeros(dims)
for a in range(slicemin, slicemax):
    y2 = input[a]
    params, pcov2 = curve_fit(func, x, y2,[100, 3, 1])
 #   z2[i] = y2[0] This is the correct z-value of each slice!
    z2[i] = 0.05*a-4.25
    R[i] = params[1]
    i = i + 1
pl.plot(z2, R, 'k.')

# <codecell>

def circle(z, R, m):
    return np.sqrt(R**2 - (z - m)**2)

# <codecell>

params2, pcov3 = curve_fit(circle, z2, R, [2.6, 2])

# <codecell>

z3 = np.zeros(49)

# <codecell>

for b in range(1,49):
    z3[b] = z3[b-1] + (b * 0.8)

# <codecell>

z3.shape

# <codecell>

pl.plot(z3, circle(z3, *params2), 'r-')
pyplot.plot(z2, R,'k.')

# <codecell>

print params2

# <codecell>

R.shape

# <codecell>

np.where?

# <codecell>



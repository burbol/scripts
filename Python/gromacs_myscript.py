# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <codecell>

# type first on the shell "source /usr/local/gromacs/bin/GMXRC"
import gromacs
gromacs.config.setup()

# <codecell>

%pylab inline

# <codecell>

import matplotlib 
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy import optimize

# <codecell>

>>> import numpy as np
>>> import gromacs.formats
>>> X = np.linspace(-10,10,50000)
>>> yerr = np.random.randn(len(X))*0.05
>>> data = np.vstack((X, np.sin(X) + yerr, np.random.randn(len(X))*0.05))
>>> xvg = gromacs.formats.XVG(array=data)

# <codecell>

%cd /Users/burbol/downloads/Testing

# <codecell>

xvg.read("g_rad_densmap_s0_w1000_8ns_10ns.xvg")

# <codecell>

xvg.plot()

# <codecell>

#colclass = matplotlib.colors.Colormap(jet, N=256)
#matplotlib.cm.register_cmap(name='jet', cmap='colclass')
#matplotlib.colors.ListedColormap([0,0,0], name='jet', N=None)
#mycolor = matplotlib.cm.get_cmap(name='jet', lut=None)
xvg.plot(columns=[0,80], maxpoints=None, alpha=1)

# <codecell>

xvg.errorbar(columns=[0,80,81], maxpoints=1000, color="red")

# <codecell>

input = xvg.array

# <codecell>

x = input[0]
print len(input), len(x)
slice = 56
y = input[slice]
print x[slice], y[0]

# <codecell>

ticki=np.zeros(len(range(-10,70,1)))
for lm in range(-10,70,1):    
    xm = 1.1
    xm = float(lm)
    xm /=10
    ticki[lm + 10] = xm
#    print type(xm)
#    print ticki[lm + 10]

# <codecell>

slicedens = np.zeros(len(input[0]))
dens = np.zeros(len(input)-1)
zcoord = np.zeros(len(input)-1)
for i in range(1,len(input)): 
    slicedens = input[i]
    dens[i-1] = sum(slicedens)
    zcoord[i-1] = 0.060 *(i-1)
    #print dens[j]
#print len(dens)
#print len(x)
fig = plt.figure()
ax = fig.add_subplot(111)
ax.plot(zcoord,dens)
ax.plot(zcoord,dens,'k.')
ax.axis([-1, 7, 0, 250000])
ax.set_xticks(ticki)
ax.set_yticks(range(0,250000,5000))
ax.grid(True)
plt.setp(ax.get_xticklabels(), visible=False)
plt.setp(ax.get_yticklabels(), visible=False)
plt.show()

# <codecell>

from matplotlib.ticker import MultipleLocator, FormatStrFormatter

majorLocator   = MultipleLocator(1)
majorFormatter = FormatStrFormatter('%d')
minorLocator   = MultipleLocator(0.2)


fig, ax = plt.subplots()
#plt.plot(t,s)

#fig = plt.figure()
#ax = fig.add_subplot(111)
plt.plot(zcoord,dens)
plt.plot(zcoord,dens,'k.')

ax.xaxis.set_major_locator(majorLocator)
ax.xaxis.set_major_formatter(majorFormatter)

#for the minor ticks, use no labels; default NullFormatter
ax.xaxis.set_minor_locator(minorLocator)

ax.xaxis.grid(True, which='minor')
ax.yaxis.grid(True, which='major')

#ax.grid(True)
plt.show()

# <codecell>

dens = numpy.array(dens)
#print max(dens)
#print x[dens.argmax()]
max_y = max(dens)  # Find the maximum density value
max_x = zcoord[dens.argmax()]  # Find the x value corresponding to the maximum density value
print max_x, max_y
zeroslice = nonzero(zcoord == max_x)[0][0]
print zeroslice
shift = max_x

# <codecell>

# DELETE THIS CELL!!
# we set a different shift value (for g_rad_densmap_s0_w1000_8ns_10ns.xvg) 
# from the one determined above because it's is not the global max
shift = 2.52
print zcoord[42]
print nonzero(zcoord == 2.52)[0][0]

# <codecell>

# THESE ARE ONLY TESTS TO FIND THE FIRST LOCAL MAX
maxdensi = xvg.max
plt.plot(xvg.max)
max_y = max(maxdensi)  # Find the maximum density value
max_x = zcoord[maxdensi.argmax()]  # Find the x value corresponding to the maximum density value
print max_x, max_y
zeroslice = nonzero(zcoord == max_x)[0][0]
print zeroslice
shift = max_x

# <codecell>

#xvg.errorbar(columns=[0,80,81])

# <codecell>

def func(r, ro, R, d):
    return (ro/2)*(1-np.tanh(2*(r-R)/d))

# <codecell>

y = input[47]
popt, pcov = curve_fit(func, x, y,[1000, 3, 1])
print popt

# <codecell>

plt.plot(x, func(x, *popt), 'r-')
plt.plot(x,y,'k.')

# <codecell>

k = func(0,*popt)
print k

# <codecell>

def checkslice(s):
    popt, pcov = curve_fit(func, x, s,[100, 3, 1])
    k = func(0,*popt)
    if k > 885:
        return True
    else: 
        return False

# <codecell>

slicemin = 0
for slice in range(0,150):
    y = input[slice]
    if checkslice(y): 
       slicemin = slice
       print slicemin
       break
    else: continue

# <codecell>

slicemax = 0
for slice in range(200,0,-1):
    y = input[slice]
    if checkslice(y): 
       slicemax = slice
       print slicemax
       break
    else: continue

# <codecell>

# USE THIS CELL TO DETERMINE BULK VALUE
densmean = np.zeros(200)
i = 0
for slice in range(1,200):
    ymean = input[slice]
    popt, pcov = curve_fit(func, x, ymean,[1000, 3, 1])
    k = func(0,*popt)
    densmean[slice]=k
#    print k
    if k > 800:
        start = slice - i
        i = i + 1

end = start + i
print "start =", start, "end =", end
print zcoord[start]
densmean2 = np. zeros(i + 1)
i2 = 0
for slice2 in range(0,200):
    ymean = input[slice2]
    popt, pcov = curve_fit(func, x, ymean,[1000, 3, 1])
    k = func(0,*popt)
    if k > 800:
        densmean2[i2] = k
        i2 = i2 + 1

        
print 'length of densmap2', len(densmean2)
plt.plot(zcoord,densmean,'k.')
import statistics
print mean(densmean2)
print median(densmean2)
print statistics.median_low(densmean2)
print statistics.median_high(densmean2)
#print statistics.mode(densmean2)

# <codecell>

## DELETE THIS CELL, IT GIVES ONLY HOW MANY TIMES DOES EACH VALUE OF densmean APPEAR -> NOT USEFUL
from collections import Counter
c = Counter(densmean)  
#print densmean
#list(c)
c += Counter()   
#print c.most_common()       # n least common elements 
#print c[-3.2638837653955745e-17]
#sum(c.values())  

# <codecell>

#slicemin2 = slicemin + 16
slicemin2 = slicemin 
slicemax = slicemax 
print zeroslice, slicemin, slicemin2, slicemax
#y = input[zeroslice]
#shift = y[0]
print shift

# <codecell>

def from_slice_circle(slicemax, slicemin2, input, zcoord):
    
    dims = slicemax - slicemin2 + 1
    maxdim = slicemax + 1
    i=0
    z = np.zeros(dims)
    R = np.zeros(dims)
    w = np.zeros(dims)
    
    def func(r, ro, R, d):
        return (ro/2)*(1-np.tanh(2*(r-R)/d))
    
    for a in range(slicemin2, maxdim):
        y2 = input[a]
        i = a - slicemin2
        z[i] = zcoord[a] - 2.52
    #    w[i] = 0.05*(a - 1) - 5.0 - shift -> Ahother way to determine z
        params, pcov2 = curve_fit(func, x, y2,[1000, 3, 1])
        R[i] = params[1]
#    plt.plot(z, R, 'k.')

#      == LEAST SQUARES CIRCLE ==
#     Advanced usage of optimize.leastsq with jacobian

    method  = "leastsq with jacobian"
#    from scipy    import optimize

    # coordinates of the barycenter

    x_m = mean(z)
    y_m = mean(R)

    def calc_R(xc, yc):
        """ calculate the distance of each data points from the center (xc, yc) """
        return sqrt((z-xc)**2 + (R-yc)**2)

    def f(c):
        """ calculate the algebraic distance between the 2D points and the mean circle centered at c=(xc, yc) """
        Ri = calc_R(*c)
        return Ri - Ri.mean()

    def Df(c):
        """ Jacobian of f
        The axis corresponding to derivatives must be coherent with the col_deriv option of leastsq"""
        xc, yc     = c
        df_dc    = empty((len(c), z.size))

        Ri = calc_R(xc, yc)
        df_dc[ 0] = (xc - z)/Ri                   # dR/dxc
        df_dc[ 1] = (yc - R)/Ri                   # dR/dyc
        df_dc       = df_dc - df_dc.mean(axis=1)[:, newaxis]

        return df_dc

    center_estimate = x_m, y_m
    center, ier = optimize.leastsq(f, center_estimate, Dfun=Df, col_deriv=True)

    xc, yc = center
    Ri        = calc_R(xc, yc)
    R_2       = Ri.mean()
    residu    = sum((Ri - R_2)**2)
    residu2   = sum((Ri**2-R_2**2)**2)
    
    #Here we move the circle down, so that the y-coordinate of its center yc = 0.
    s = yc*np.ones(len(R))
    R = R - s
    yc = yc - yc
    return xc, yc, R_2, residu, residu2, z, R

# <codecell>

def distance(z, R, xc, yc):
    a = array([z,R])
    b = array([xc,yc])
    d = np.linalg.norm(a - b)
    return d

# <codecell>

slicemin = 47
slicemax = 92
epsilon = 0.02
xc, yc, R_2, residu, residu2, z, R= from_slice_circle(slicemax, slicemin, input, zcoord)
print xc, yc, R_2, z[0], R[0], residu, residu2
print "shift is", shift
slicesnr = slicemax - slicemin
print len(z), len(R), slicesnr
# For slicemin
d = distance(z[0], R[0], xc, yc)-R_2
print d
#if d <= epsilon:

# For slicemax
d = distance(z[slicesnr], R[slicesnr], xc, yc)-R_2
print d

# <codecell>

theta_fit = linspace(-pi, pi, 180)
x_fit2 = xc + R_2*cos(theta_fit)
y_fit2 = yc + R_2*sin(theta_fit)
plt.plot(x_fit2, y_fit2, color="blue", label='method', lw=2)
# plot(X, C, color="blue", linewidth=2.5, linestyle="-")
# Set x limits
xlim(-3,3)
# Set y limits
ylim(0,4)
plt.plot(z, R, 'k.')

# <codecell>

-1.79326905749 0.0898083934024 3.00014070197 -0.81 0.0227674754812 0.807332782233

The Contact Angle is   53.292626884
The base radius is   2.40520899696

-2.31124195341 0.0310632265358 3.06574017087 -1.05 0.00371061172299 0.139538471801


The Contact Angle is   41.0712096432
The base radius is   2.01418058477


-2.31600176047 0.0184026567677 3.0787098635 -1.14 0.00341249773831 0.129385137185


The Contact Angle is   41.2131486412
The base radius is   2.02844528375

# <codecell>

def angle(R, m):
    t = m / np.sqrt(R**2 - m**2)
    t2 = math.degrees(np.arctan(t)) + 90
    return t2

# <codecell>

def radius(R, m):
    r = np.sqrt(R**2 - m**2)
    return r

# <codecell>

theta = angle(R_2, xc)
base_r = radius(R_2, xc)
print "The Contact Angle is  ", theta
print "The base radius is  ", base_r

# <codecell>

a = array([0,0])
b = array([2,2])
print np.linalg.norm(a - b)

# <codecell>

a = [1, 2, 3]
b = 4*np.ones(3)
c = a - b
print c

# <codecell>



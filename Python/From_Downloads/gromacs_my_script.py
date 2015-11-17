#!/usr/bin/python

# type first on the shell "source /usr/local/gromacs/bin/GMXRC"
import gromacs
gromacs.config.setup()
import matplotlib 
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy import optimize
import numpy as np
import gromacs.formats

os.chdir('/Users/burbol/downloads/Testing')

peaks = np.loadtxt('WaterPeaksFinal.txt', skiprows=1)

xvg.read("g_rad_densmap_s0_w1000_8ns_10ns.xvg")

xvg.plot()

xvg.plot(columns=[0,80], maxpoints=None, alpha=1)

input = xvg.array

x = input[0]
y = input[slice]

slicedens = np.zeros(len(input[0]))
dens = np.zeros(len(input)-1)
zcoord = np.zeros(len(input)-1)
for i in range(1,len(input)): 
    zcoord[i-1] = 0.060 *(i-1)

def func(r, ro, R, d):
    return (ro/2)*(1-np.tanh(2*(r-R)/d))

k = func(0,*popt)

def checkslice(s):
    popt, pcov = curve_fit(func, x, s,[100, 3, 1])
    k = func(0,*popt)
    if k > 885:
        return True
    else: 
        return False

slicemin = 0
for slice in range(0,150):
    y = input[slice]
    if checkslice(y): 
       slicemin = slice
       print slicemin
       break
    else: continue

slicemax = 0
for slice in range(200,0,-1):
    y = input[slice]
    if checkslice(y): 
       slicemax = slice
       print slicemax
       break
    else: continue

# FUNCTION TO DETERMINE BULK VALUE OF WATER DENSITY
def waterbulk(input, zcoord, func):
    x = input[0]
    densmean = np.zeros(200)
    i = 0
    for slice in range(1,200):
        ymean = input[slice]
        popt, pcov = curve_fit(func, x, ymean,[1000, 3, 1])
        k = func(0,*popt)
        densmean[slice]=k
        if k > 800:
            start = slice - i
            i = i + 1

    end = start + i
    densmean2 = np. zeros(i + 1)
    i2 = 0
    for slice2 in range(0,200):
        ymean = input[slice2]
        popt, pcov = curve_fit(func, x, ymean,[1000, 3, 1])
        k = func(0,*popt)
        if k > 800:
            densmean2[i2] = k
            i2 = i2 + 1
    return median(densmean2)

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

def distance(z, R, xc, yc):
    a = array([z,R])
    b = array([xc,yc])
    d = np.linalg.norm(a - b)
    return d

slicesnr = slicemax - slicemin

# For slicemin
d = distance(z[0], R[0], xc, 0)-R_2
if d <= epsilon:
   while (d <= epsilon):
        oldslicemin = slicemin
        slicemin = oldslicemin + 1
        xc, yc, R_2, residu, residu2, z, R= from_slice_circle(slicemax, slicemin, input, zcoord)
        d = distance(z[0], R[0], xc, 0)-R_2
   if d > epsilon:
      slicemin = oldslicemin
elif d > epsilon:
   while (d > epsilon):
        oldslicemin = slicemin
        slicemin = oldslicemin + 1
        xc, yc, R_2, residu, residu2, z, R= from_slice_circle(slicemax, slicemin, input, zcoord)
        d = distance(z[0], R[0], xc, 0)-R_2
#   if d <= epsilon:
#      break
print "slicemin is", slicemin, "and its d is", d

# For slicemax
slicesnr = slicemax - slicemin -1
d = distance(z[slicesnr], R[slicesnr], xc, 0)-R_2
if d <= epsilon:
   while (d <= epsilon):   
       oldslicemax = slicemax
       slicemax = oldslicemax + 1
       xc, yc, R_2, residu, residu2, z, R= from_slice_circle(slicemax, slicemin, input, zcoord)
       slicesnr = slicemax - slicemin
       d = distance(z[slicesnr], R[slicesnr], xc, 0)-R_2
   if d > epsilon:
        slicemax = oldslicemax
elif d > epsilon:
   while (d > epsilon):
       oldslicemax = slicemax
       slicemax = oldslicemax + 1
       xc, yc, R_2, residu, residu2, z, R= from_slice_circle(slicemax, slicemin, input, zcoord)
       slicesnr = slicemax - slicemin -1
       d = distance(z[slicesnr], R[slicesnr], xc, 0)-R_2
#   if d <= epsilon:
#      break
            
print "slicemax is", slicemax, "and its d is", d
    
     
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



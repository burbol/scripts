#!/usr/bin/python

# THIS PROGRAM CALCULATES,FROM DENSITY FILES CREATED WITH GROMACS,THE CONTACT ANGLE AND 
# THE RADIUS AT THE BASE OF A WATER DROPLET WHICH HAS BEEN SIMULATED BEEING ON TOP OF A SURFACE.

import gromacs
gromacs.config.setup()
import matplotlib 
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
import os
import math
import pylab as pl
import gromacs.formats
from matplotlib.backends.backend_pdf import PdfPages

# Main directory where all the density files are stored, and where the results will be stored too
os.chdir("/home/eixeres/Downloads/Berendsen")

#Read file with first water peak positions and save in variable "peaks"
peaks = np.loadtxt('WaterPeaksFinal.txt', skiprows=1)

# Name of the file where the Contact Angles will be saved
myfile = open('Contact_Angles.txt', 'w')
myfile2 = open('Contact_Angles2.txt', 'w')

# First line of the output file: the titles of the 3 columns of data
print >> myfile, '{0}  {1}  {2}'.format('File', 'Contact Angle', 'Base Radius')
print >> myfile2, '{0}  {1}'.format('Contact Angle', 'Base Radius')

# Path to the pdf-file, where the plots of the fitted circles will be stored
pp = PdfPages('/home/eixeres/Downloads/Berendsen/Results.pdf')

# This function allows to use ranges, just as the built-in function range(), but with float arguments.
# Taken from "http://code.activestate.com/recipes/66472-frange-a-range-function-with-float-increments/"
def frange(start, end=None, inc=None):
    "A range function, that does accept float increments..."

    if end == None:
        end = start + 0.0
        start = 0.0

    if inc == None:
        inc = 1.0

    L = []
    while 1:
        next = start + len(L) * inc
        if inc > 0 and next >= end:
            break
        elif inc < 0 and next <= end:
            break
        L.append(next)
        
    return L

# "peaknum" will be the position of the array "peaks" 
peaknum = 0

#Here we read from a file the positions of the first water peaks of each system and save them in the array "peak"
peak = np.loadtxt('WaterPeaksFinal.txt', skiprows=1)

# We create the dictionary "systems" which relates the systems with the positions of the array "peakcount". 
# For example: "systems[(11, 3000)]" would give back the position 22, which can be used to call "peakcount[22]" 
# and get back the shift value of the system with a SAM  with 11% OH- coverage and a water droplet of 3000 molecules.
peakpos = 0
for b in [0, 5, 11, 17]:
    for c in [1000, 2000, 3000, 4000, 5000, 6500, 7000, 8000, 9000, 10000]:
        z = b, c
        systems[z] = peakpos
        peakpos = peakpos + 1
#print systems[(11, 3000)]

# Loop through all the density files in the main directory. In this case "b" stands for the OH-percentage of the SAM,
# and "c" stands for the number of molecules (2 for 2000, 3 for 3000, etc)
peakcount = 0
for b in [0 5, 11, 17]:
    for c in [1, 2, 3, 4]:
        os.chdir('/home/eixeres/Downloads/Berendsen')
        foldername = 'densmaps_s%d_w%d000'%(b,c, )
        os.chdir(foldername)
        print >> myfile, '             '
        peakcount = systems[(b, c)]
        shift = peak[peakcount]
        for d in frange(0, 20, 0.5):
            l = d + 0.5
            l = str(l)
            d = str(d)            
            filename = 'g_rad_densmap_s%d_w%d_%sns_%sns.xvg'%(b,c,d,l, )
            xvg.read(filename)
            input = xvg.array

#Here we create the array "zcoord" with the z-coordinates of the slices taking into account that they are 0.060nm appart
            zcoord = np.zeros(len(input)-1)
            for i in range(1,len(input)): 
                zcoord[i-1] = 0.060 *(i-1)

# The function "checkzero" checks, given a specific slice "s", if there is some density value different then zero. 
# In that case it returns "True", if all the values are zero, it returns "False" (it ignores the first value of the input array, 
# because it only stores the slice position and not any density information).
            x = input[0]
            def checkzero(s):
                i = len(np.transpose(np.nonzero(s)))
                if i > 1: 
                   return True
                else: 
                   return False

# Here we use the function "checkzero" to find the first slice where there is some point of density different
# then zero. Then we set that slice number as the variable "zeroslice". We can change the slices being checked (from 70 to 150). 
            zeroslice = 0
            for slice in range(70,150):
                y = input[slice]
                if checkzero(y): 
                   zeroslice = slice
                   break
                else: continue

# The function "func" is the sigmoid fit of the density data points.
            def func(r, ro, R, d):
                return (ro/2)*(1-np.tanh(2*(r-R)/d))
                
# Searching for initian guess for the first parameter of the fit
            def guess(y):
                m = max(y)
                m2 = min(y)
                third = (m - m2)/3
                bottom = m - third
                tot = 0
                count = 0
                for i in range(0, len(y)):
                    if y[i] >= bottom:
                       tot = tot + y[i]
                       count = count + 1
                At = tot / count
                return At
 
 
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
                
# The function "checkslice" has 2 arguments: "s" is the slice number, and "t" is the initial guess 
# of the first parameter of the fit. It makes a fit of the density data, and the checks if the value of 
# the fitted function at R=0 is greater then 90. If so, it returns "True". If not, it returns "False".
# "maxfev" is the allowed number of funtion calls when doing the fit. The default number was 800.
            def checkslice(s, t):
                popt, pcov = curve_fit(func, x, s,[t, 2, 1], maxfev=10000)
                k = func(0,*popt)
                bulk = waterbulk(input, zcoord, func)
                if k > bulk:
                   return True
                else: 
                   return False
                   
# Here we use "checkslice" to find the first slice with density (at R=0) greater the 90. We call that slice "slicemin"
            slicemin = 0
            for slice in range(70,150):
                y = input[slice]
                At = guess(y)
                if checkslice(y, At):
		           slicemin = slice
                   break
                else: continue               

# Here we use "checkslice" to find the last slice with density (at R=0) greater the 90. We call that slice "slicemax"
            slicemax = 0
            for slice in range(179,130,-1):
                y = input[slice]
                At = guess(y)
                if checkslice(y, At):
                   slicemax = slice
                   break
                else: continue

# "slicemin2" is the first slice that we will actually analyze. 
            slicemin2 = zeroslice + 16
            if slicemin2 < slicemin:
	           slicemin2 = slicemin
	    
# Some output on the screen to follow what the program is doing
	    print "Checking ", filename
	    print "slices found: zero, min, zero+16, max"
	    print zeroslice, slicemin, slicemin2, slicemax
	    
# Now we determine R and z. 
# R comes from the fit, and z is a vector with the shifted positions of the slices.
            dims = slicemax - slicemin2 + 1
            maxdim = slicemax + 1
            i=0
            z = np.zeros(dims)
            R = np.zeros(dims)
            w = np.zeros(dims)
            for a in range(slicemin2, maxdim):
                y2 = input[a]
                i = a - slicemin2
                z[i] = y2[0] - shift
          #    w[i] = 0.05*(a - 1) - 5.0 - shift -> Ahother way to determine z
                params, pcov2 = curve_fit(func, x, y2,[100, 3, 1])
                R[i] = params[1]

# Here we use part of the source code of a website, to fit points to a circle.--> "http://wiki.scipy.org/Cookbook/Least_Squares_Circle"

#  == LEAST SQUARES CIRCLE ==
# Advanced usage of optimize.leastsq with jacobian

            method  = "leastsq with jacobian"
            from scipy    import optimize

             # coordinates of the barycenter
            x_m = np.mean(z)
            y_m = np.mean(R)

            def calc_R(xc, yc):
#    """ calculate the distance of each data points from the center (xc, yc) """
                return np.sqrt((z-xc)**2 + (R-yc)**2)

            def f(c):
#    """ calculate the algebraic distance between the 2D points and the mean circle centered at c=(xc, yc) """
                Ri = calc_R(*c)
                return Ri - Ri.mean()

            def Df(c):
#    """ Jacobian of f
#    The axis corresponding to derivatives must be coherent with the col_deriv option of leastsq"""
                xc, yc     = c
                df_dc    = np.empty((len(c), z.size))

                Ri = calc_R(xc, yc)
                df_dc[ 0] = (xc - z)/Ri                   # dR/dxc
                df_dc[ 1] = (yc - R)/Ri                   # dR/dyc
                df_dc       = df_dc - df_dc.mean(axis=1)[:, np.newaxis]
                return df_dc

            center_estimate = x_m, y_m
            center, ier = optimize.leastsq(f, center_estimate, Dfun=Df, col_deriv=True)

# The results from the fitted circle, with the coordinates of the circle center: xc, yc, and Circle radius: R_2.
            xc, yc = center
            Ri        = calc_R(xc, yc)
            R_2       = Ri.mean()
            residu    = sum((Ri - R_2)**2)

# The function "plot_circ" makes a plot of the fitted circle together with the points that have been fitted. 
# The arguments: z is the shifted position of each slice, R the calculated radius corresponding to each slice,
# xc and yc are the coordinates of the center of the fitted circle, and R_2 is the radius of the fitted circle.
# Other options regarding the aspect of the plots should be added/changed here inside de function.
            def plot_circ(z, R, xc, yc, R_2):                    
                fig = pl.figure()           
                theta_fit = np.linspace(-np.pi, np.pi, 180)
                x_fit2 = xc + R_2*np.cos(theta_fit)
                y_fit2 = yc + R_2*np.sin(theta_fit)               
                pl.plot(x_fit2, y_fit2, 'b-', z, R, 'k.')
                fig.suptitle(filename)
                pl.xlim(0,5)
                pl.ylim(0,5)
                pl.xlabel('z [nm]')
                pl.ylabel('R [nm]')
                pl.grid(True)
                return fig
    

# The function "angle" calculates the contact angle.
            def angle(R, m):
                t = m / np.sqrt(R**2 - m**2)
                t2 = math.degrees(np.arctan(t)) + 90
                return t2

# The function "radius" calculates the base radius of the droplet.
            def radius(R, m):
                r = np.sqrt(R**2 - m**2)
                return r

# Here are our final results, with "theta" as the contact angle, and "base_r" as the base radius.
            theta = angle(R_2, xc)
            base_r = radius(R_2, xc)

# The results are printed to the opened file.
            print >> myfile, '{0}  {1}  {2}'.format(filename, theta, base_r)
            print >> myfile2, '{0}'.format('File')
	    print >> myfile2, '{0}  {1}'.format(theta, base_r)

# We call the function "plot_circ" and save each plot on a single page inside the same file.            
            plot1 =  plot_circ(z, R, xc, yc, R_2)
            pp.savefig(plot1)
            

myfile.close()
pp.close()



#!/usr/bin/python

import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
import os
import math
import pylab as pl

os.chdir("/home/eixeres/Downloads/Berendsen")
myfile = open('Contact_Angles.txt', 'w')
print >> myfile, '{0}  {1}  {2}'.format('File', 'Contact Angle', 'Base Radius')

for b in [5]:
    for c in [2]:
        os.chdir("/home/eixeres/Downloads/Berendsen")
        foldername = 'densmaps_s%d_w%d000'%(b,c, )
        os.chdir(foldername)
        print >> myfile, '             '
        for d in range(0, 20, 2):
            l = d + 2
            filename = 'densmap_%dpc_%d000_%dns_%dns.dat'%(b,c,d,l, )
            



# os.chdir('/Users/burbol/Desktop/Testing/files_for_laila/densmaps_berendsen')
# os.chdir('densmaps_s5_w2000/')

# input = np.loadtxt("densmap_5pc_2000_12ns_14ns.dat")

            input = np.loadtxt(filename)

            x = input[0]
            def checkzero(s):
                i = len(np.transpose(np.nonzero(s)))
                if i > 1: 
                   return True
                else: 
                   return False

            zeroslice = 0
            for slice in range(70,150):
                y = input[slice]
                if checkzero(y): 
                   zeroslice = slice
                   break
                else: continue

            def func(r, ro, R, d):
                return (ro/2)*(1-np.tanh(2*(r-R)/d))
                
            #searching for initian guess for the first parameter of the fit
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
                
# ATTENTION: change the parameters of curve_fit if necessary!!
            def checkslice(s, t):
                popt, pcov = curve_fit(func, x, s,[t, 2, 1], maxfev=10000)
                k = func(0,*popt)
                if k > 90:
                   return True
                else: 
                   return False
                   

            slicemin = 0
            for slice in range(70,150):
                y = input[slice]
                #print "slicemin: I'm checking slice number:", slice, "of file:", filename
                At = guess(y)
                if checkslice(y, At):
		   slicemin = slice
                   break
                else: continue               
        
            slicemax = 0
            for slice in range(179,130,-1):
                y = input[slice]
                #print "slicemax: I'm checking slice number:", slice, "of file:", filename
                At = guess(y)
                if checkslice(y, At):
                   slicemax = slice
                   break
                else: continue

            slicemin2 = zeroslice + 16
            if slicemin2 < slicemin:
	      slicemin2 = slicemin
	    
            print filename
	    print "slices: zero, min, zero+16, max", zeroslice, slicemin, slicemin2, slicemax
            y = input[zeroslice]
            shift = y[0]

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

            xc, yc = center
            Ri        = calc_R(xc, yc)
            R_2       = Ri.mean()
            residu    = sum((Ri - R_2)**2)

# print "Circle center:" ,xc, yc
# print "Circle radius:", R_2

            def plot_circ(z, R, xc, yc, R_2):
                theta_fit = np.linspace(-np.pi, np.pi, 180)
                x_fit2 = xc + R_2*np.cos(theta_fit)
                y_fit2 = yc + R_2*np.sin(theta_fit)
                d2 = d/2
                l2 = l/2
                ax = pl.subplot(10,1,d2)
                pl.plot(x_fit2, y_fit2, color="blue", label=method, lw=2)
                ax.set_xlim(0,4.0)
                ax.set_ylim(0,3)
                ax.plot(z, R, 'k.')
                fig = plt.gcf()
                fig.set_size_inches(18.5,10.5)
    

            def angle(R, m):
                t = m / np.sqrt(R**2 - m**2)
                t2 = math.degrees(np.arctan(t)) + 90
                return t2

            def radius(R, m):
                r = np.sqrt(R**2 - m**2)
                return r

            theta = angle(R_2, xc)
            base_r = radius(R_2, xc)
#            print "The Contact Angle is  ", theta
#            print "The base radius is  ", base_r
            

            print >> myfile, '{0}  {1}  {2}'.format(filename, theta, base_r)
            plot_circ(z, R, xc, yc, R_2)



plt.savefig("/home/eixeres/Downloads/Berendsen/Results.pdf", dpi = 100)
myfile.close()




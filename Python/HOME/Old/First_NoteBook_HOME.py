#!/usr/bin/python

import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
import os
import math
import pylab as pl

os.chdir('/Users/burbol/Desktop/Testing/files_for_laila/densmaps_berendsen')
os.chdir('densmaps_s5_w2000/')

input = np.loadtxt("densmap_5pc_2000_12ns_14ns.dat")

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
popt, pcov = curve_fit(func, x, y,[100, 3, 1])


def checkslice(s):
    popt, pcov = curve_fit(func, x, s,[100, 3, 1])
    k = func(0,*popt)
    if k > 90:
        return True
    else: 
        return False

slicemin = 0
for slice in range(70,150):
    y = input[slice]
    if checkslice(y): 
       slicemin = slice
       break
    else: continue

slicemax = 0
for slice in range(200,130,-1):
    y = input[slice]
    if checkslice(y): 
       slicemax = slice
       break
    else: continue

slicemin2 = slicemin + 16
print "slices: zero, min, min+16, max", zeroslice, slicemin, slicemin2, slicemax
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
    """ calculate the distance of each data points from the center (xc, yc) """
    return np.sqrt((z-xc)**2 + (R-yc)**2)

def f(c):
    """ calculate the algebraic distance between the 2D points and the mean circle centered at c=(xc, yc) """
    Ri = calc_R(*c)
    return Ri - Ri.mean()

def Df(c):
    """ Jacobian of f
    The axis corresponding to derivatives must be coherent with the col_deriv option of leastsq"""
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
    theta_fit = linspace(-pi, pi, 180)
    x_fit2 = xc + R_2*cos(theta_fit)
    y_fit2 = yc + R_2*sin(theta_fit)
    pl.plot(x_fit2, y_fit2, color="blue", label=method, lw=2)
    xlim(0,4.0)
    ylim(0,3)
    pl.plot(z, R, 'k.')


def angle(R, m):
    t = m / np.sqrt(R**2 - m**2)
    t2 = math.degrees(np.arctan(t)) + 90
    return t2

def radius(R, m):
    r = np.sqrt(R**2 - m**2)
    return r

theta = angle(R_2, xc)
base_r = radius(R_2, xc)
print "The Contact Angle is  ", theta
print "The base radius is  ", base_r




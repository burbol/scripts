# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <codecell>

%pylab inline

# <codecell>

from scipy    import optimize

# <codecell>

#! /usr/bin/env python

"""
http://www.scipy.org/Cookbook/Least_Squares_Circle
"""
from numpy import *

# Coordinates of the 2D points

x = r_[ 1.25, 1.3, 1.35, 1.4, 1.45, 1.5, 1.55, 1.6, 1.65, 1.7, 1.75, 1.8, 1.85, 1.9, 1.95, 2., 2.05, 2.1, 2.15, 2.2, 2.25, 2.3, 2.35, 2.4, 2.45, 2.5, 2.55, 2.6, 2.65, 2.7, 2.75, 2.8, 2.85, 2.9, 2.95, 3., 3.05, 3.1, 3.15, 3.2, 3.25, 3.3, 3.35, 3.4, 3.45, 3.5, 3.55, 3.6, 0.]
y = r_[ 1.93981182, 2.00310173, 1.99274388, 2.0116961, 2.0416055, 2.03737162, 2.03590569, 2.04440828, 2.05315168, 2.05268199, 2.04403347, 2.04894743, 2.04726235, 2.04443174, 2.04967703, 2.04804559, 2.02013272, 2.0184885, 2.01394277, 2.00961606, 1.996566, 1.97486746, 1.95327747, 1.96275451, 1.93143096, 1.91394371, 1.91355069, 1.87852594, 1.84949912, 1.83748102, 1.80284304, 1.7891118, 1.76267059, 1.72186055, 1.68523042, 1.66175036, 1.646192, 1.58001539, 1.53496118, 1.49672587, 1.42635437, 1.41212638, 1.34340227, 1.27809465, 1.19550884, 1.1508165, 1.07445865, 0.94906878, 0.]
basename = 'circle'

# x = r_[36, 36, 19, 18, 33, 26]
# y = r_[14, 10, 28, 31, 18, 26]
# basename = 'arc'



# == METHOD 1 ==
method_1 = 'algebraic'

# coordinates of the barycenter
x_m = mean(x)
y_m = mean(y)

# calculation of the reduced coordinates
u = x - x_m
v = y - y_m

# linear system defining the center in reduced coordinates (uc, vc):
#    Suu * uc +  Suv * vc = (Suuu + Suvv)/2
#    Suv * uc +  Svv * vc = (Suuv + Svvv)/2
Suv  = sum(u*v)
Suu  = sum(u**2)
Svv  = sum(v**2)
Suuv = sum(u**2 * v)
Suvv = sum(u * v**2)
Suuu = sum(u**3)
Svvv = sum(v**3)

# Solving the linear system
A = array([ [ Suu, Suv ], [Suv, Svv]])
B = array([ Suuu + Suvv, Svvv + Suuv ])/2.0
uc, vc = linalg.solve(A, B)

xc_1 = x_m + uc
yc_1 = y_m + vc

# Calculation of all distances from the center (xc_1, yc_1)
Ri_1      = sqrt((x-xc_1)**2 + (y-yc_1)**2)
R_1       = mean(Ri_1)
residu_1  = sum((Ri_1-R_1)**2)
residu2_1 = sum((Ri_1**2-R_1**2)**2)

# Decorator to count functions calls
import functools
def countcalls(fn):
    "decorator function count function calls "

    @functools.wraps(fn)
    def wrapped(*args):
        wrapped.ncalls +=1
        return fn(*args)

    wrapped.ncalls = 0
    return wrapped




# == METHOD 2b ==
# Advanced usage of least squares, with jacobian
method_2b  = "leastsq with jacobian"

def calc_R(xc, yc):
    """ calculate the distance of each 2D points from the center c=(xc, yc) """
    return sqrt((x-xc)**2 + (y-yc)**2)

@countcalls
def f_2b(c):
    """ calculate the algebraic distance between the 2D points and the mean circle centered at c=(xc, yc) """
    Ri = calc_R(*c)
    return Ri - Ri.mean()

@countcalls
def Df_2b(c):
    """ Jacobian of f_2b
    The axis corresponding to derivatives must be coherent with the col_deriv option of leastsq"""
    xc, yc     = c
    df2b_dc    = empty((len(c), x.size))

    Ri = calc_R(xc, yc)
    df2b_dc[ 0] = (xc - x)/Ri                   # dR/dxc
    df2b_dc[ 1] = (yc - y)/Ri                   # dR/dyc
    df2b_dc       = df2b_dc - df2b_dc.mean(axis=1)[:, newaxis]

    return df2b_dc

center_estimate = x_m, y_m
center_2b, ier = optimize.leastsq(f_2b, center_estimate, Dfun=Df_2b, col_deriv=True)

xc_2b, yc_2b = center_2b
Ri_2b        = calc_R(xc_2b, yc_2b)
R_2b         = Ri_2b.mean()
residu_2b    = sum((Ri_2b - R_2b)**2)
residu2_2b   = sum((Ri_2b**2-R_2b**2)**2)
ncalls_2b    = f_2b.ncalls

print """
Method 2b :
print "Functions calls : f_2b=%d Df_2b=%d""" % ( f_2b.ncalls, Df_2b.ncalls)

# == METHOD 3 ==
# Basic usage of odr with an implicit function definition
from scipy      import  odr

method_3  = "odr"

@countcalls
def f_3(beta, x):
    """ implicit definition of the circle """
    return (x[0]-beta[0])**2 + (x[1]-beta[1])**2 -beta[2]**2

# initial guess for parameters
R_m = calc_R(x_m, y_m).mean()
beta0 = [ x_m, y_m, R_m]

# for implicit function :
#       data.x contains both coordinates of the points
#       data.y is the dimensionality of the response
lsc_data   = odr.Data(row_stack([x, y]), y=1)
lsc_model  = odr.Model(f_3, implicit=True)
lsc_odr    = odr.ODR(lsc_data, lsc_model, beta0)
lsc_out    = lsc_odr.run()

xc_3, yc_3, R_3 = lsc_out.beta
Ri_3       = calc_R(xc_3, yc_3)
residu_3   = sum((Ri_3 - R_3)**2)
residu2_3  = sum((Ri_3**2-R_3**2)**2)
ncalls_3   = f_3.ncalls


# Summary
fmt = '%-22s %10.5f %10.5f %10.5f %10d %10.6f %10.6f %10.2f'
print ('\n%-22s' +' %10s'*7) % tuple('METHOD Xc Yc Rc nb_calls std(Ri) residu residu2'.split())
print '-'*(22 +7*(10+1))
print  fmt % (method_1 , xc_1 , yc_1 , R_1 ,        1 , Ri_1.std() , residu_1 , residu2_1 )
# print  fmt % (method_2 , xc_2 , yc_2 , R_2 , ncalls_2 , Ri_2.std() , residu_2 , residu2_2 )
print  fmt % (method_2b, xc_2b, yc_2b, R_2b, ncalls_2b, Ri_2b.std(), residu_2b, residu2_2b)
# print  fmt % (method_3 , xc_3 , yc_3 , R_3 , ncalls_3 , Ri_3.std() , residu_3 , residu2_3 )
# print  fmt % (method_3b, xc_3b, yc_3b, R_3b, ncalls_3b, Ri_3b.std(), residu_3b, residu2_3b)

# plotting functions
from matplotlib                 import pyplot as p, cm, colors
p.close('all')

def plot_all(residu2=False):
    """ Draw data points, best fit circles and center for the three methods,
    and adds the iso contours corresponding to the fiel residu or residu2
    """

    f = p.figure( facecolor='white')  #figsize=(7, 5.4), dpi=72,
    p.axis('equal')

    theta_fit = linspace(-pi, pi, 180)

    x_fit1 = xc_1 + R_1*cos(theta_fit)
    y_fit1 = yc_1 + R_1*sin(theta_fit)
    p.plot(x_fit1, y_fit1, 'b-' , label=method_1, lw=2)

#    x_fit2 = xc_2b + R_2b*cos(theta_fit)
#    y_fit2 = yc_2b + R_2b*sin(theta_fit)
#    p.plot(x_fit2, y_fit2, 'k--', label=method_2, lw=2)

#    x_fit3 = xc_3 + R_3*cos(theta_fit)
#    y_fit3 = yc_3 + R_3*sin(theta_fit)
#    p.plot(x_fit3, y_fit3, 'r-.', label=method_3, lw=2)

    p.plot([xc_1], [yc_1], 'bD', mec='y', mew=1)
    p.plot([xc_2b], [yc_2b], 'gD', mec='r', mew=1)
#    p.plot([xc_3], [yc_3], 'kD', mec='w', mew=1)

    # draw
    p.xlabel('x')
    p.ylabel('y')



    # plot data
    p.plot(x, y, 'ro', label='data', ms=8, mec='b', mew=1)
    p.legend(loc='best',labelspacing=0.1 )

    p.xlim(xmin=vmin, xmax=vmax)
    p.ylim(ymin=vmin, ymax=vmax)

    p.grid()
    p.title('Least Squares Circle')
    p.savefig('%s_residu%d.png' % (basename, 2 if residu2 else 1))

plot_all(residu2=False)
plot_all(residu2=True )

p.show()
# vim: set et sts=4 sw=4:

# <codecell>



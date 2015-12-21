#!/usr/bin/python

# THIS PYTHON PROGRAM USES THE (.xvg) DENSITY FILES CREATED WITH g_rad_density, A MODIFIED VERSION OF THE GROMACS PROGRAM g_densmap.
# IT CALCULATES THE CONTACT ANGLE BETWEEN DIFFERENT SIZED WATER DROPLETS AND VARIOUS SAMS WITH DIFFERENT PERCETAGES OF OH-COVERAGE.
# THE RADIUS OF THE CIRCULAR SURFACE AT THE BASE OF THE WATER DROPLETS IS ALSO DETERMINED.
# THE RESULTS ARE SAVED IN TWO .TXT FILES AND ONE PLOT OF EACH SYSTEM IN A DIFFERENT PAGE OF A .PDF FILE.
# (THE SIMULATIONS OF THE SYSTEMS OF DROPLET + SAM HAVE BEEN DONE WITH GROMACS.)

import gromacs
gromacs.config.setup()
import matplotlib 
import matplotlib.pyplot as plt
import numpy as np
from numpy import linalg
from scipy.optimize import curve_fit
from scipy import optimize
from scipy import odr
import os
import math
import pylab as pl
from matplotlib.backends.backend_pdf import PdfPages
from gromacs.formats import XVG

# This function allows to use ranges, just as the built-in function range(), but with float arguments.
# Taken from "http://code.activestate.com/recipes/66472-frange-a-range-function-with-float-increments/"
def frange(start, end=None, inc=None):
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

#SAMs=[33,50]
#Waters=[1000, 2000, 3000, 4000, 5000, 6500, 7000, 8000, 9000, 10000]
SAMs=[5]
Waters=[2000]

start_time=0
end_time=100

# Main directory where all the density files are stored, and where the results will be stored too
densfolder="/Volumes/UNI/SHELDON/g_rad_densmaps_s33s50"
#densfolder="/media/SHELDON/g_rad_densmaps_s33s50"
#densfolder="/media/SHELDON/Last_Simulations_Dec14/FINISHED_1_Jan_2015/g_rad_densmaps/"

#File with the positions of z=0
#peaksinput='MiddlePoint_s33s50.txt' #z=0 is the middle point between the first water peak  and the last SAM peak of each system
#peaksinput='SAMPeaks_s33s50.txt' #z=0 is the last SAM peak of each system
#peaksinput='WaterPeaks_s33s50.txt' #z=0 is the first water peak of each system
peaksinput='/Users/burbol2/Dropbox/Apps/Computable/Definitions_Interface/All_systems_without66/GDS_positions_CORRECTED.txt'

# Files where the Contact Angles will be saved. WE USE 2 FILES: txtouput has
# a third column with the names of the density maps files.txtoutput2 is better for analysis
#txtoutput='Contact_Angles_MiddlePoint_s'+str(SAMs[0])+'.txt'
#txtoutput2='Contact_Angles2_MiddlePoint_s'+str(SAMs[0])+'.txt'
#txtoutput='Contact_Angles_SAMPeak_s'+str(SAMs[0])+'_10000.txt'
#txtoutput2='Contact_Angles2_SAMPeak_s'+str(SAMs[0])+'_10000.txt'
#txtoutput='Contact_Angles_WaterPeak_s'+str(SAMs[0])+'.txt'
#txtoutput2='Contact_Angles2_WaterPeak_s'+str(SAMs[0])+'.txt'
txtoutput='Contact_Angles_GDS_s'+str(SAMs[0])+'.txt'
txtoutput2='Contact_Angles2_GDS_s'+str(SAMs[0])+'.txt'

# Pdf-file, where the plots of the fitted circles will be stored
#plotsOutput='MiddlePeakCircle_plots_s'+str(SAMs[0])+'.pdf'
#plotsOutput='SAMPeakCircle_plots_s'+str(SAMs[0])+'.pdf'
#plotsOutput='WaterPeakCircle_plots_s'+str(SAMs[0])+'_10000.pdf'
plotsOutput='GDSCircle_plots_all_s'+str(SAMs[0])+'.pdf'

# ZLENGTH is the system total length in the z-direction Change this value if necessary.
ZLENGTH = 20.0 # *(for the "older" systems use ZLENGTH =12.0)

# systems that appear in the peaks txt file peaksinput
SAMpeaks=[0,5,11,17,33,50]
Waterpeaks=[1000, 2000, 3000, 4000, 5000, 6500, 7000, 8000, 9000, 10000]

os.chdir(densfolder)

#Here we read the z=0 positions and save them in the array "peak"
peak = np.loadtxt(peaksinput, skiprows=1)

# We create the dictionary "systems" which relates the systems with the positions of the array "peakcount". 
# For example: "systems[(11, 3000)]" would give back the position 22, which can be used to call "peakcount[22]" 
# and get back the shift value of the system with a SAM	 with 11% OH- coverage and a water droplet of 3000 molecules.
# This is done for the case of using this script with only a subgroup of systems
peakpos = 0
systems={}
for b in SAMpeaks:
	for c in Waterpeaks:
		z = b, c
		systems[z] = peakpos
		peakpos = peakpos + 1
#print systems[(11, 3000)]

# Create the PdfPages object to which we will save the pages:
# The with statement makes sure that the PdfPages object is closed properly at
# the end of the block, even if an Exception occurs.
with PdfPages(plotsOutput) as pp, open(txtoutput, 'w') as myfile, open(txtoutput2, 'w') as myfile2:

# these 3 lines were needed before inserting the "with" statement
#pp = PdfPages(plotsOutput) 
#myfile = open(txtoutput, 'a')
#myfile2 = open(txtoutput2, 'a')

	# We write the first line of the output file: the titles of the 3 columns of data
	print >> myfile, '{0}  {1}	{2}'.format('File', 'Contact Angle', 'Base Radius')
	print >> myfile2, '{0}	{1}'.format('Contact Angle', 'Base Radius')

	# Loop through all the density files in the main directory. In this case "b" stands for the OH-percentage of the SAM,
	# and "c" stands for the number of molecules (2 for 2000, 3 for 3000, etc)
	# We save each time the data of the density file in the variable "input" and
	# the radial coordinates in the variable "x"
	peakcount = 0
	for b in SAMs:
		for c in Waters:
			print >> myfile, '             '
			peakcount = systems[(b, c)]
			shift = peak[peakcount]
			print >> myfile2, '{0}	{1} {2}'.format('#File:', b, c)
			for d in frange(start_time, end_time, 0.5):
				l = d + 0.5
				l = str(l)
				d = str(d)
				xvg = XVG()
				filename = 'g_rad_dmap_%dpc_w%d_%sns_%sns.xvg'%(b,c,d,l, )
				if not os.path.isfile(filename):
					print >> myfile, '{0}  {1}  {2}'.format(filename, 'nan', 'nan')
					print >> myfile2, '{0}  {1}'.format('nan', 'nan')
					continue 
				xvg.read(filename)
				input = xvg.array
				x = input[0]
	# Here we create the array "zcoord" with the z-coordinates of the slices subtracting the shift value.
				slicestotal = len(input) - 1
				S = ZLENGTH/ slicestotal
				zcoord = np.zeros(len(input)-1)
				for i in range(1,len(input)): 
					zcoord[i-1] = (S *(i-1))-shift

	# We set the number of the first slice as "zeroslice", which will be used later for leaving the first 8nm above the surface out of the fitted data.
				zeroslice = 0
				cnt = 0
				for slicepos in zcoord:
					if slicepos >= 0.0: 
						zeroslice = cnt
						break
					cnt = cnt + 1

	# The function "func" is the sigmoid fit of the density data points.
				def func(r, ro, R, d):
					return (ro/2)*(1-np.tanh(2*(r-R)/d))
				
	# Searching for initial guess for the first parameter of the fit
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
	# It takes 3 arguments: 
	# input: all the data of the density map of one system
	# zcoord: the z-coordinates of the slices of the system
	# func: sigmoidal function used to fit the density data of one slice
				def waterbulk(input, zcoord, func):
					x = input[0]
					densmean = np.zeros(200)
					i = 0
					start = 0
					for slice in range(1,200):
						ymean = input[slice]
						At = guess(ymean)					
						popt, pcov = curve_fit(func, x, ymean,[At, 2, 1], maxfev=10000)
						k = func(0,*popt)
						densmean[slice]=k
						if k > 800:
							start = slice - i
							i = i + 1

					end = start + i
					densmean2 = np.zeros(i + 1)
					i2 = 0
					for slice2 in range(0,200):
						ymean = input[slice2]
						At = guess(ymean)
						popt, pcov = curve_fit(func, x, ymean,[At, 2, 1], maxfev=10000)
						k = func(0,*popt)
						if k > 800:
							densmean2[i2] = k
							i2 = i2 + 1
					return np.median(densmean2)
 
	# We use te function waterbulk, and set the variable "bulk" as the returned value.
				bulk = waterbulk(input, zcoord, func)
				print "bulk value=", bulk

	# The function "checkslice" makes a fit of the density data, and the checks if the value of 
	# the fitted function at R=0 is greater then 80% of the bulk value. If so, it returns "True". If not, it returns "False".
	# "maxfev" is the allowed number of funtion calls when doing the fit. The default number was 800.
	# It takes 3 arguments as input:
	# s: (type: list) density data of a specific slice i, in our program it would be "input[i]"
	# t: (type: float) the initial guess of the first parameter of the fit.
				def checkslice(s, t, bulk):
					popt, pcov = curve_fit(func, x, s,[t, 2, 1], maxfev=10000)
					k = func(0,*popt)
					bulklimit = 0.9*bulk
					if k > bulklimit:
						return True
					else: 
						return False
				
	# Here we use "checkslice" to find the first slice with density (at R=0) greater then 80% of the bulk value. We call that slice "slicemin"
				slicemin = 0
				for slice in range(0,150):
					y = input[slice]
					At = guess(y)
					if checkslice(y, At, bulk):
						slicemin = slice
						break
					else: 
						continue 
				#print "SLICEMIN DETERMINED", slicemin
	# We leave the first 8nm above the surface out. If slicemin is already 0.8nm or more above the surface, we leave slicemin without change.
	#			slicemin2 = zeroslice + 13
	#			if slicemin2 > slicemin:
	#				slicemin = slicemin2

	# Here we use "checkslice" to find the last slice with density (at R=0) greater then 80% of the bulk value. We call that slice "slicemax"
				slicemax = 0
				for slice in range(200,0,-1):
					y = input[slice]
					At = guess(y)
					if checkslice(y, At, bulk):
						slicemax = slice
						break
					else:
						continue
				#print "SLICEMAX DETERMINED -->",slicemax

	#We leave the slices below zero out
				if zeroslice > slicemin:
					slicemin = zeroslice 

	# Some output on the screen to follow what the program is doing
				print "Checking ", filename
				print "slices changed to: zero, min, max"
				print zeroslice, slicemin, slicemax


	# FUNCTION "from_slice_circle" makes a circle fit of some consecutive slices
	# slicemin: number of the first slice
	# slicemax: number of the last slice
	# input: array containing all the data of the density map of all the slices 
	# zcoord: list with the z-coordinates of all the slices
	# returns 7 variables: 
	# xc, yc, R_2: the coordinates of the circle center and its radius (xc corresponds to the z-position and yc to the radial distance to the central axis of the drop). 
	# But yc won't be used for later calculations and has to be set to yc=0).
	# residu, residu2: error and quadratic error of the circle fit 
	# z, R: coordinates of the resulting points of the fitted radial density profile of each slice between slicemin and slicemax
	# We use part of an imported source code to fit points to a circle. Source code can be found here--> "http://wiki.scipy.org/Cookbook/Least_Squares_Circle"

# The functions calc_R and f_3 will be used in the function from_slice_circle
				def calc_R(xc, yc):
					""" calculate the distance of each 2D points from the center c=(xc, yc) """
					return np.sqrt((x-xc)**2 + (y-yc)**2)

				def f_3(beta, x):
					""" implicit definition of the circle """
					return (x[0]-beta[0])**2 + (x[1])**2 -beta[2]**2
					
				def from_slice_circle(slicemax, slicemin, input, zcoord):
				# Now we determine R and z. 
				# R comes from the fit, and z is a vector with the shifted positions of the slices.
					dims = slicemax - slicemin + 1
					x = input[0]
					i=0
					z = np.zeros(dims)
					R = np.zeros(dims)
					for a in range(slicemin, slicemax +1):
						y2 = input[a]
						i = a - slicemin
						z[i] = zcoord[a-1] 
						At = guess(y2)
						params, pcov2 = curve_fit(func, x, y2,[At, 2, 1], maxfev=10000)
						R[i] = params[1]
						
					# coordinates of the barycenter
					x_m = np.mean(z)
					y_m = np.mean(R)
					x = z
					y = R
					
					# initial guess for parameters
					R_m = calc_R(x_m, 0).mean()
					beta0 = [ x_m, 0, R_m]
					
					# for implicit function :
					# data.x contains both coordinates of the points
					# data.y is the dimensionality of the response
					lsc_data= odr.Data(np.row_stack([x, y]), y=1)
					lsc_model  = odr.Model(f_3, implicit=True)
					lsc_odr	   = odr.ODR(lsc_data, lsc_model, beta0)
					lsc_out	   = lsc_odr.run()
					
					xc_3, yc_3, R_3 = lsc_out.beta
					Ri_3 = calc_R(xc_3, yc_3)
					residu_3 = sum((Ri_3 - R_3)**2)
					residu2_3 = sum((Ri_3**2-R_3**2)**2)

					return xc_3, yc_3, R_3, residu_3, residu2_3, z, R


	# The function "distance" calculates the distance between the fitted circle and a point. It takes 5 argumets:
	# x, y are the coordinates of the point; xc, yc are the coordinates of the circle center; R is its radius.
				def distance(x, y, xc, yc,R):
					a = np.asarray([x,y])
					b = np.asarray([xc,yc])
					d = np.linalg.norm(a - b)- R
					d = math.fabs(d)
					return d 
 
	 # Now we search for the "best" slices to make the circle fit and calculate later the contact angle and base radius.
	 # We look for the first and last fitted LG-interface points (related to the first and last droplet slices: slicemin and slicemax) with a distance to the circle smaller then epsilon, which we set to 0.029nm.
						
				epsilon = 0.03
				epsilon2 = 0.031
	# For slicemin
				#print "Starting the search for the best slicemin. We start with slicemin =", slicemin
				xc, yc, R_2, residu, residu2, z, R= from_slice_circle(slicemax, slicemin, input, zcoord)
				d = distance(z[0], R[0], xc, yc, R_2)
				if d <= epsilon:
					#print "d is smaller then epsilon: d=", d
					while d <= epsilon:
						oldslicemin = slicemin
						slicemin = oldslicemin - 1
						#print "subtracted 1 from slicemin. Now slicemin=", slicemin
						old_d = d
						xc, yc, R_2, residu, residu2, z, R= from_slice_circle(slicemax, slicemin, input, zcoord)
						d = distance(z[0], R[0], xc, yc, R_2)
						#print "New d=", d
				elif d > epsilon:
					#print "d is bigger then epsilon: d=", d
					while d > epsilon:
						oldslicemin = slicemin
						slicemin = oldslicemin + 1
						#print "added 1 to slicemin. Now slicemin=", slicemin
						old_d = d
						xc, yc, R_2, residu, residu2, z, R= from_slice_circle(slicemax, slicemin, input, zcoord)
						d = distance(z[0], R[0], xc, yc, R_2)
						#print "New d=", d
				if d > epsilon:
					#print "We subtracted one too much. d has become bigger then epsilon ->d=", d
					slicemin = oldslicemin
					d = old_d
					#print "We add one slice again. Now slicemin=", slicemin, "And d=", d

				#print "Out of search loop: slicemin=", slicemin, "and its d=", d 

	# For slicemax
				slicesnr = slicemax - slicemin
				#print "Starting the search for the best slicemax. We start with slicemax =", slicemax, "And slicesnr=", slicesnr
				xc, yc, R_2, residu, residu2, z, R= from_slice_circle(slicemax, slicemin, input, zcoord)
				d = distance(z[-1], R[-1], xc, yc, R_2)
				if d <= epsilon2:
					#print "d is smaller then epsilon2: d=", d
					while d <= epsilon2:	  
						oldslicemax = slicemax
						slicemax = oldslicemax + 1
						slicesnr = slicemax - slicemin +1
						#print "added 1 to slicemax. Now slicemax=", slicemax, "And slicesnr=", slicesnr
						old_d = d
						xc, yc, R_2, residu, residu2, z, R= from_slice_circle(slicemax, slicemin, input, zcoord)
						slicesnr = slicemax - slicemin - 1
						d = distance(z[-1], R[-1], xc, yc, R_2)
						#print "New d=", d
				elif d > epsilon2:
					#print "d is bigger then epsilon2: d=", d
					while d > epsilon2:
						oldslicemax = slicemax
						slicemax = oldslicemax - 1
						slicesnr = slicemax - slicemin +1
						#print "subtracted 1 from slicemax. Now slicemax=", slicemax, "And slicesnr=", slicesnr
						old_d = d
						xc, yc, R_2, residu, residu2, z, R= from_slice_circle(slicemax, slicemin, input, zcoord)
						slicesnr = slicemax - slicemin -1
						d = distance(z[-1], R[-1], xc, yc, R_2)
						#print "New d=", d
				if d > epsilon2:
					#print "We added one too much. d has become bigger then epsilon ->d=", d
					slicemax = oldslicemax
					d = old_d
					#print "We subtract one slice again. Now slicemax=", slicemax, "And d=", d

				#print "Out of search loop: slicemax=s", slicemax, "and its d is", d

	#We leave (again) the slices below zero out of the fit
				if zeroslice > slicemin:
					slicemin = zeroslice 

	# Once we have the "correct" values of slicemin and slicemax we make the circle fit again
				xc, yc, R_2, residu, residu2, z, R= from_slice_circle(slicemax, slicemin, input, zcoord)
				print "xc, yc, R_2, residu, residu2 =", xc, yc, R_2, residu, residu2

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
				print "TESTING: theta, base_r=", theta, base_r

	# The results are printed to the opened file.
				print >> myfile, '{0}  {1}	{2}'.format(filename, theta, base_r)
				print >> myfile2, '{0}	{1}'.format(theta, base_r)
		
	# The function "plot_circ" makes a plot of the fitted circle together with the points that have been fitted. 
	# The arguments: z is the shifted position of each slice, R the calculated radius corresponding to each slice,
	# xc and yc are the coordinates of the center of the fitted circle, and R_2 is the radius of the fitted circle.
	# It returns the figure object "fig".
	# Other options regarding the aspect of the plots should be added/changed here inside de function.
				def plot_circ(z, R, xc, yc, R_2):
					fig = pl.figure()
					theta_fit = np.linspace(-np.pi, np.pi, 180)
					x_fit2 = xc + R_2*np.cos(theta_fit)
					y_fit2 = yc + R_2*np.sin(theta_fit)
					pl.plot(y_fit2, x_fit2, 'b-', R, z, 'k.')
					fig.suptitle(filename)
					pl.xlim(0,base_r+ 1.0)
					pl.ylim(0,base_r+ 1.0)
					#pl.xlim(0,3*base_r) # for big drops
					#pl.ylim(0,3*base_r) # for big drops
					pl.xlabel('R [nm]')
					pl.ylabel('z [nm]')
					pl.grid(True)
					return fig
	
	# We call the function "plot_circ" and save each plot on a single page inside the same file.			
				plot1 =	 plot_circ(z, R, xc, yc, R_2)
				pp.savefig(plot1)
				pl.close() # Is this line necessary after using the "with" statement??
			

			#myfile.close() With the "with" statement, there is no need to close the files anymore!
			#pp.close()


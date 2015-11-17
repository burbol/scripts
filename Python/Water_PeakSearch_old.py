#!/usr/bin/python

import gromacs
gromacs.config.setup()
import matplotlib 
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import numpy as np
import os
import math
import pylab as pl
import gromacs.formats
from matplotlib.backends.backend_pdf import PdfPages

os.chdir('/Users/burbol/Desktop/scripts/Python')
import peakdetect
pk = peakdetect.peakdetect
os.chdir('/Users/burbol/Downloads/global_density_maps')
xvg = gromacs.formats.XVG()


# Name of the file where the Contact Angles will be saved
myfile1 = open('WaterPeaks1.txt', 'w')
print >> myfile1, '{0}  {1}'.format('File', 'Water Peaks')
myfile = open('WaterPeaks.txt', 'w')
print >> myfile, '{0}'.format('Water Peaks Position')

pp = PdfPages('/Users/burbol/Downloads/global_density_maps/densplots.pdf')

# Loop through all the density files in the main directory
for b in [0, 5, 11, 17]:
    for c in [1000, 2000, 3000, 4000, 5000, 6500, 7000, 8000, 9000, 10000]:
        os.chdir('/Users/burbol/Downloads/global_density_maps')
        print >> myfile1, '             '
        print >> myfile, '             '
        filename = 'g_density_NVT_sam%d_water%d.xvg'%(b,c, ) 
        xvg.read(filename)
        input = xvg.array
        x = input[0]
        y = input[1]
        [MAXTAB, MINTAB] = pk(y, x, 1, 8)

# The results are printed to the opened file.
        print >> myfile1, '{0}  {1}'.format(filename, MAXTAB[0])
        print >> myfile, '{0}'.format(filename, MAXTAB[0])            
        plot =  xvg.plot()
        pp.savefig(plot)

myfile.close()
pp.close()



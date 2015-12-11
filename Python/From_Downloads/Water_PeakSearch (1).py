#!/usr/bin/python

#source /usr/local/gromacs/bin/GMXRC

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
from matplotlib.ticker import MultipleLocator, FormatStrFormatter

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

ticki=np.zeros(24)
for lm in range(0,24): 
    xm = lm*0.5   
    ticki[lm] = xm 

def plot_dens(x, y, filename, point):                   

    majorLocator   = MultipleLocator(1)
    majorFormatter = FormatStrFormatter('%d')
    minorLocator   = MultipleLocator(0.2)
  
#    fig = pl.figure()   
    fig, ax = plt.subplots()                        
    pl.plot(x, y, 'b-', x, y, 'k.')
    fig.suptitle(filename)
    # pl.plot(z, R, 'k.')
    #fig.set_xlim(0,4.0)
    #fig.set_ylim(0,3)                
    # fig = plt.gcf()
    plt.xticks(ticki)
#    plt.xticks(range(0,d_max_y,50))
    pl.xlim(0,12)
    pl.ylim(0,max(y))
    pl.xlabel('z [nm]')
    pl.ylabel('Density [kg/nm3]')
#    pl.grid(b=True, which='both', axis='both', color='r', linestyle='-', linewidth=0.5)
    # fig.set_size_inches(18.5,10.5)
    pl.annotate('first peak', xy=(point[0], point[1]),  xycoords='data', xytext=(-50, 30), textcoords='offset points', arrowprops=dict(arrowstyle="->"))

    ax.xaxis.set_major_locator(majorLocator)
    ax.xaxis.set_major_formatter(majorFormatter)

    #for the minor ticks, use no labels; default NullFormatter
    ax.xaxis.set_minor_locator(minorLocator)

    ax.xaxis.grid(True, which='minor')
    ax.yaxis.grid(True, which='major')

    return fig
    

# Loop through all the density files in the main directory
for b in [0, 5, 11, 17]:
    for c in [1000, 2000, 3000, 4000, 5000, 6500, 7000, 8000, 9000, 10000]:
        os.chdir('/Users/burbol/Downloads/global_density_maps')
        print >> myfile1, '             '
#        print >> myfile, '             '
        filename = 'g_density_NVT_sam%d_water%d.xvg'%(b,c, ) 
        xvg.read(filename)
        input = xvg.array
        x = input[0]
        y = input[1]
        [MAXTAB, MINTAB] = pk(y, x, 1, 8)
        point = MAXTAB[0]

# The results are printed to the opened file.
        print >> myfile1, '{0}  {1}'.format(filename, MAXTAB[0])
        print >> myfile, '{0}'.format(point[0]) 

        plot1 =  plot_dens(x, y, filename, point)
        pp.savefig(plot1)


pp.close()
myfile.close()



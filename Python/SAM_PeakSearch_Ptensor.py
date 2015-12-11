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

#os.chdir('/Users/burbol/Desktop/scripts/Python')
import peakdetect
pk = peakdetect.peakdetect

# The function "peakdetect" is included in the module "peakdetect.py", 
# which can be found at https://gist.github.com/sixtenbe/1178136.
# This is a a copy of its original dlastpeakescription:
# peakdetect(y_axis, x_axis = None, lookahead = 300, delta=0):

"""
    Converted from/based on a MATLAB script at: 
    http://billauer.co.il/peakdet.html
    
    function for detecting local maximas and minmias in a signal.
    Discovers peaks by searching for values which are surrounded by lower
    or larger values for maximas and minimas respectively
    
    keyword arguments:
    y_axis -- A list containg the signal over which to find peaks
    x_axis -- (optional) A x-axis whose values correspond to the y_axis list
        and is used in the return to specify the postion of the peaks. If
        omitted an index of the y_axis is used. (default: None)
    lookahead -- (optional) distance to look ahead from a peak candidate to
        determine if it is the actual peak (default: 200) 
        '(sample / period) / f' where '4 >= f >= 1.25' might be a good value
    delta -- (optional) this specifies a minimum difference between a peak and
        the following points, before a peak may be considered a peak. Useful
        to hinder the function from picking up false peaks towards to end of
        the signal. To work well delta should be set to delta >= RMSnoise * 5.
        (default: 0)
            delta function causes a 20% decrease in speed, when omitted
            Correctly used it can double the speed of the function
    
    return -- two lists [max_peaks, min_peaks] containing the positive and
        negative peaks respectively. Each cell of the lists contains a tupple
        of: (position, peak_value) 
        to get the average peak value do: np.mean(max_peaks, 0)[1] on the
        results to unpack one of the lists into x, y coordinates do: 
        x, y = zip(*tab)
"""

#os.chdir('/Users/burbol/Downloads/global_SAMS_densmaps')
#os.chdir('/Volumes/Backup/YosemiteFiles/MEGAsync/Density_Profiles/ptensor/ptensor_density_profiles')
os.chdir('/media/SHELDON/CLOSEDsheldon/eixeres/files_for_laila/ptensor/density_profiles')
xvg = gromacs.formats.XVG()

SAMs=[0, 5, 11, 17,33,50,66]
mindens=500.0 # change this value to get different peaks!!


# Name of the file where the Contact Angles will be saved
myfile1 = open('SAMPeaks_Ptensor1.txt', 'w')
print >> myfile1, '{0}  {1}'.format('File', 'Last SAM\'s Peak')
myfile = open('SAMPeaks_Ptensor.txt', 'w')
print >> myfile, '{0}'.format('Last SAM\'s Peak Position')

pp = PdfPages('SAMs_densplots_Ptensor.pdf')

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
    pl.xlim(0,3)
    pl.ylim(0,10000)
    pl.xlabel(r'$z [nm]$')
    pl.ylabel(r'$\rho (z) [\frac{kg}{{nm}^{3}}]$')
#    pl.grid(b=True, which='both', axis='both', color='r', linestyle='-', linewidth=0.5)
    # fig.set_size_inches(18.5,10.5)
    pl.annotate('last peak', xy=(point[0], point[1]),  xycoords='data', xytext=(-50, 30), textcoords='offset points', arrowprops=dict(arrowstyle="->"))

    ax.xaxis.set_major_locator(majorLocator)
    ax.xaxis.set_major_formatter(majorFormatter)

    #for the minor ticks, use no labels; default NullFormatter
    ax.xaxis.set_minor_locator(minorLocator)

    ax.xaxis.grid(True, which='minor')
    ax.yaxis.grid(True, which='major')

    return fig
    

# Loop through all the density files in the main directory
for b in SAMs:
    print >> myfile1, '             '
#    print >> myfile, '             '
    filename = 'dens2_NPT_sam%d_water_ptensor.xvg'%(b, )
    xvg.read(filename)
    input = xvg.array
    x = input[0]
    y = input[1]
# Here we search for all the peaks. In this case, all the points (density values) different then zero are already density peaks
# The lists "peaks_all" and "peaks_all_pos" store, respectively, the density peak values and their positions.
    peaks_all = []
    peaks_all_pos = []
    peak_number = 0
    for m in range(0,999):
        check = y[m]
        check_pos =x[m]
        if check > mindens: 
            peak_number = peak_number +1
            peaks_all.append(check)
            peaks_all_pos.append(check_pos)

# Now we check if the last peak is bigger then the mean density value, because we don't want peaks that are too small.
# When the "correct" peak is selected, we store its value in the list variable "peak".
    lastpeak = peaks_all[-1]
    lastpeak_pos = peaks_all_pos[-1]
    j=1
    while lastpeak < np.mean(peaks_all):
        lastpeak = peaks_all[-1-j]
        lastpeak_pos = peaks_all_pos[-1-j]
        j = j + 1
    peak = [lastpeak_pos, lastpeak]     

# The results are printed to the opened file.
    print >> myfile1, '{0}  {1}'.format(filename, lastpeak_pos)
    print >> myfile, '{0}'.format(lastpeak_pos) 

    plot1 =  plot_dens(x, y, filename, peak)
    pp.savefig(plot1)


pp.close()
myfile.close()



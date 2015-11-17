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

# The function "peakdetect" is included in the module "peakdetect.py", 
# which can be found at https://gist.github.com/sixtenbe/1178136.
# This is a a copy of its original description:
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

    majorLoca
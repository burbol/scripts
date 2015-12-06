#!/usr/bin/python

'''
This Python script computes the contact angle between water 
droplets and SAMs simulated with GROMACS. The radius of the 
circular surface at the base of the droplets is also determined. 
The script is designed to loop through multiple input density 
profiles.

Input:
Files with format ".xvg" produced with "g_rad_density"

Output: 
        -- Two files with format ".txt" (contact angle and base 
        radius in two columns. One of them has a third column
        with the name of the input file. Each line corresponds
        to one input file).
        -- One file with format ".pdf" (circle fit plots. Each 
        page corresponds to one input file).
        Dependencies: C_AnglesModule.py, gromacs wrapper, scipy
        
(g_rad_density is a modification of the GROMACS program g_densmap
created in the research group AG Netz at the FU Berlin. The
difference is that g_rad_density uses the same coordinates as 
the simulation files while g_densmaps sets the center at the com).

It uses the self-developed module C_AnglesModule.py
Also the module GromacsWrapper from:
Oliver Beckstein, http://github.com/orbeckst/GromacsWrapper
'''
import C_AnglesModule
import numpy as np
import os
import math
from matplotlib.backends.backend_pdf import PdfPages
import gromacs
gromacs.config.setup()
from gromacs.formats import XVG

#parentfolder = "/Volumes/UNI/SHELDON/g_rad_densmaps/"
#peaksfolder = "/Users/burbol2/Dropbox/Apps/Computable/Definitions_Interface/NewNames/"
#outputfolder = "/Users/burbol2/Desktop/"

#parentfolder = "/net/data/eixeres/NewVersion4/FINISHED/"
parentfolder = "/Volumes/UNI/Densmaps_NewVersion4/"
peaksfolder = parentfolder
outputfolder = parentfolder

# Systems to analyze:

#SAM percentage of -OH coverage (polarity)
SAMs=[25]
#Waters=[3000]

#Number of water molecules in the droplet
Waters=[2000, 3000, 4000, 5000, 6500, 7000, 8000, 9000]


#Length of longest interval analysed. 
#For shorter simulations 'nan' is returned.
start_time = 0
end_time = 100
timestep = 0.5


# Simulation box length in the z-direction:
boxlength={5: 20.0, 21: 20.0, 25: 20.0, 41: 20.0}


# Select position of solid-liquid interface (z = 0).
# a) for Water Peak
# b) for GDS
# c) for SAM Peak
# d) for Middle Point
# e) for highest carbon in SAM

option= 'e'
if option == 'a':
    interface='WaterPeaks'
elif option == 'b':
    interface='GDS'
elif option == 'c':
    interface='SAMPeaks'
elif option == 'd':
    interface='MiddlePoint'
elif option == 'e':
    interface='SAM_last_C_atom'

# File with interface positions:
file_pos_interface = peaksfolder + interface + '.txt'

# systems that appear in file with interface positions:
SAMpeaks = [5,21,25,41]
Waterpeaks = Waters
#SAMpeaks=[0,5,11,17,21,25,33,41,50,66]
#Waterpeaks=[2000, 3000, 4000, 5000, 6500, 7000, 8000, 9000]
            
#Read file with interface positions.
interface_data= np.loadtxt(file_pos_interface, skiprows=1)

# Save interface position of each system in dictionary for easier
# access. The keys are the polarity and molecule number of the system.
# For example: interface_pos[(11, 1000)]
i = 0
interface_pos={}
for pc in SAMpeaks:
    for molecs in Waterpeaks:
        z = pc, molecs
        interface_pos[z] = interface_data[i]
        i = i + 1


os.chdir(parentfolder)
for pc in SAMs:

    txtoutput2 = outputfolder + 'Contact_Angles2_' + interface +'s'+str(pc)+'.txt'
    txtoutput = outputfolder + 'Contact_Angles_' + interface +'_s'+str(pc)+'.txt'
    plotsOutput = outputfolder + interface + 'Circle_plots_s'+str(pc)+'.pdf'

    #with PdfPages(plotsOutput) as pp, open(txtoutput, 'w') as myfile, open(txtoutput2, 'w') as myfile2:
    # the former line was changed for the next two because of "AttributeError: __exit__" when using 
    # PdfPages with a with-statement, probably due to an older python version on the server
    with open(txtoutput, 'w') as myfile, open(txtoutput2, 'w') as myfile2:
        pp = PdfPages(plotsOutput)
        
        # Write titles for the 3 columns of data
        
        print >> myfile, '{0}  {1}  {2}'.format(
        'File', 'Contact Angle', 'Base Radius')
        print >> myfile2, '{0}  {1}'.format(
        'Contact Angle', 'Base Radius')
        
        
        # Loop through input files and save in density in variable
        # "myinput" and radial coordinate in variable "x".
        
        for molecs in Waters:
        
            simulation_folder = parentfolder + 's'+str(pc) + '_w' + str(molecs)
            
            # If simulation folder doesn't exist jump to the next time step
            if not os.path.isdir(simulation_folder):
                continue
                
            os.chdir(simulation_folder)
            print >> myfile, '             '
            print >> myfile2, '{0}  {1} {2}'.format('#File:', pc, molecs)
            
            for start in C_AnglesModule.frange(
            start_time, end_time, timestep):
                
                # Set input file name.
                end = start + timestep
                filename = 'g_rad_dmap_%dpc_w%d_%sns_%sns.xvg'%(
                pc, molecs, str(start), str(end), )
                
            
                # If file doesn't exist print "nan" and jump to the 
                # next time step
                
                if not os.path.isfile(filename): 
                   print >> myfile, '{0}  {1}   {2}'.format(
                   filename, 'nan', 'nan')
                   print >> myfile2, '{0}   {1}'.format('nan', 'nan')
                   continue
                                  
                # Read density files and save in multidimensional 
                # array (matrix). The first element of array (first 
                # column of matrix) are the z-coordinates.
                
                xvg = XVG()   
                xvg.read(filename)
                myinput = xvg.array
                x = myinput[0] 
                
                
                # Create array "zcoord" subtracting the position of 
                # the solid-liquid interface (shift) from z-coordinates.
                
                slicestotal = len(myinput) - 1
                sliceheight = boxlength[pc]/ slicestotal
                zcoord = np.zeros(len(myinput)-1)
                
                shift = interface_pos[(pc, molecs)]
                
                for i in range(1,len(myinput)): 
                    zcoord[i-1] = (sliceheight *(i-1))-shift
        
        
                # Set slice with zcoord = 0.0 as "zeroslice". Needed
                # to discard first 8 nm above the surface from fit.
                
                zeroslice = 0
                cnt = 0
                for slicepos in zcoord:
                    if slicepos >= 0.0: 
                       zeroslice = cnt
                       break
                    cnt = cnt + 1
                
                # Arrays with arameters of sigmoid fit for all slices:
                
                r0, R, d = C_AnglesModule.sigmoidfit(myinput)
                
                
                # Bulk density.
                
                bulk = C_AnglesModule.waterbulk(myinput, r0, R, d )
                
                
                # Fiirst estimate of last slice "good" for fit.
                
                slicemax = 0
                for slice in range(slicestotal,0,-1):
                    
                    if C_AnglesModule.checkslice(
                    r0[slice], R[slice], d[slice], bulk):
                        
                        slicemax = slice
                        break
                    else:
                        continue
                        
                        
                # First estimate of first slice "good" for fit (slicemin).
                
                slicemin = 0
                middleslice = int(0.75*slicestotal)
                
                for slice in range(0,middleslice):
                    
                    if C_AnglesModule.checkslice(
                    r0[slice], R[slice], d[slice], bulk):
                        slicemin = slice
                        break
                    else: 
                        continue 
                        
  
                # Shift "slicemin" if too close to surface.
                # 0.8nm above "zeroslice" (if needed).
                
                discard = 0.8 # use nm!
                slices_out =  int(discard/sliceheight)
                slicemin2 = zeroslice + 13
                
                if slicemin < slicemin2:
                    slicemin = slicemin2

        
                #Leave slices below zero out (if needed).
                
                if slicemin < zeroslice :
                    slicemin = zeroslice 
                
                               
                if ((slicemin==None) and (slicemax==None)) or (slicemin >= slicemax - 1):
                    print >> myfile, '{0}  {1}  {2}'.format(
                    filename,'nan', 'nan')
                    print >> myfile2, '{0}  {1}'.format('nan', 'nan')
                    print "ERROR 2"
                    continue
                
                # Boundaries for fit:
                
                slicemin, slicemax = C_AnglesModule.find_best_slices(
                slicemax, slicemin, zcoord, R)                   
        
                #Leave slices below zero out.
                
                if zeroslice > slicemin:
                    slicemin = zeroslice 
                
                
                # Make the circle fit again.
                xc, yc, Rc, z_circle, R_circle  = \
                    C_AnglesModule.from_slice_circle(
                    slicemax, slicemin,zcoord, R)
                    
                    
                if math.isnan(xc) or math.isnan(Rc):
                    print >> myfile, '{0}  {1}  {2}'.format(
                    filename,'nan', 'nan')
                    print >> myfile2, '{0}  {1}'.format('nan', 'nan')
                    print "ERROR 3"
                    continue
    
                # Final results: 

                theta = C_AnglesModule.angle(Rc, xc)  #contact angle
                base_r = C_AnglesModule.radius(Rc, xc)  #base radius
                
                
                # Print results on screen.
                print "File analyzed ", filename
                print "Slices used: zero slice=",  zeroslice, 
                "min=",  slicemin, "max=",  slicemax
                print "\n Results: theta=", theta, "base_r=", base_r
        
        
                # Save results in output files.
                
                print >> myfile, '{0}  {1}  {2}'.format(
                filename, theta, base_r)
                print >> myfile2, '{0}  {1}'.format(theta, base_r)
        
        
                # Plot circle with fitted points.
                
                circle_fit_plot =  C_AnglesModule.plot_circ_points(
                z_circle, R_circle, xc, yc, Rc, filename)
                
                
                #Save one plot per page.
                
                pp.savefig(circle_fit_plot)
                
                
                #pp.close()
    pp.close()

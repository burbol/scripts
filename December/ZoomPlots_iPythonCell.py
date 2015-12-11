b=0

minblocksize=20 # The minimum block size is 10ns=20 points
beg = 10 # We leave out the first 5ns=10 points
blocksNum=3 #the function best_start uses also 3 blocks for averaging

Nrows = 5
Ncolumns = 2
fig, ax = plt.subplots(Ncolumns, Nrows,figsize=(8,4),dpi=1000)


#fig.subplots_adjust(bottom=0.5,top=2.5, left=0.0, right=1.5, hspace = 0.5)
fig.subplots_adjust(bottom=0.0,top=4.5, left=0.0, right=1.5, hspace = 0.35)
matplotlib.rcParams['legend.handlelength'] = 0
matplotlib.rcParams['legend.markerscale'] = 0
matplotlib.rcParams['legend.handletextpad'] = 0
matplotlib.rcParams['legend.markerscale'] = 2 

#matplotlib.rcParams['font.family'] = 'Times New Roman'    
#matplotlib.rcParams['font.family'] = 'Times New Roman Bold' 

errortheta_w[b]=[]
theta_w[b]=[]
start_theta_w[b]=[]
end_theta_w[b]=[]
equil_theta_w[b]=[]

i = 0
for c in Waters:
#for c in [1000]:
    n=Waters.index(c)

    theta2 = array(angles_w[(b, c)])

    end=endpoint(theta2)
    start=best_start(theta2,t,beg,minblocksize)
    slope, intercept, delete1, delete2, delete3 = stats.linregress(t[start:end],theta2[start:end])
    thetamean, errrorbar =blockAverage((theta2[start:end]),blocksNum)
    print "for SAM ",b,"% and ", c, " molecules:"
    if start == None:
        print " Not equilibrated!"
        sampling = 0
        start = beg
        end = endpoint(theta2)
        symb='<'
    else:
        sampling = (end-start)/2.
        print "start point=", start, "length of sampling =", sampling, "ns"
        print "slope=", slope, "intercept=",intercept
        print "average=",thetamean, "error=",errrorbar
        symb='>'

    ax = plt.subplot(Nrows, Ncolumns, i+1)

    m=str(round(abs(slope)*(end-start)/2.,3))
    
    line1, = ax.plot(t[start:end], func(t[start:end], slope, intercept),'r', label=str(sampling)+r'ns: $\sigma_{t} =$ '+str(round(2*errrorbar,3))+'$^\circ $'+symb+'  '+m+'$^\circ=$'+' '+'$ m_{t}$')
    ax.plot(t[start:end], theta2[start:end],'-',linewidth=0.5)   
    ax.plot(t[start:end], func(t[start:end], 0, thetamean))
    ax.plot(t[start:end], func(t[start:end], 0, thetamean+errrorbar),'g-')
    ax.plot(t[start:end], func(t[start:end], 0, thetamean-errrorbar),'g-')   
    
    matplotlib.rcParams['font.family'] = 'Arial'
    
    # Set common labels
    ax.set_xlabel(r'$t$ [ns]',fontsize=16)
    plt.ylabel(r'$\theta_{mic}$ [deg]',fontsize=17)
    #ax.set_ylabel(r'R$\ _{BASE} $ [nm]')
    
    #titles for each subplot:
    matplotlib.rcParams['font.family'] = 'Times New Roman Bold' 
    rcParams['mathtext.default']=u'it'
    #####ax.set_title(str(c)+' water molec: $\Delta(t)=$'+str(sampling)+'ns')
    ax.set_title(r'$\Phi=$'+str(b)+', $n_{W}=$'+str(c),fontsize=19)

    matplotlib.rcParams['font.family'] = 'Arial' 
    # Create a legend for the first line.
    #first_legend = plt.legend(handles=[line1],bbox_to_anchor=(0.2, 0.9), loc=3,borderaxespad=0.,borderpad=0.2,fontsize=11)
    first_legend = plt.legend(loc=1,borderaxespad=0.,borderpad=0.2,fontsize=14)
    # Add the legend manually to the current Axes.
    ax = plt.gca().add_artist(first_legend)

    
    # SAVING RESULTS
    if sampling != 0:
        (errortheta_w[b]).append(errrorbar)
        (theta_w[b]).append(thetamean)
        (start_theta_w[b]).append(start)
        (end_theta_w[b]).append(end)
        (equil_theta_w[b]).append(c)
    
    i = i+1
 
matplotlib.rcParams['font.family'] = 'Times New Roman' 
#mytitle = plt.suptitle(r'Contact Angle $\theta_{mic}$ vs. Time for SAM'+str(b)+'%', fontsize=19, fontweight='bold',x=0.7, y=4.65)
###mytitle = plt.suptitle(r'Contact Angle $\theta_{mic}$ vs. Time for SAM'+str(b)+'%', fontsize=17,x=0.7, y=4.65)
plt.show()
fig.savefig('equil_theta_t_s'+str(b)+'.jpg',bbox_extra_artists=(first_legend,mytitle,), bbox_inches='tight')
# check for better quality:
# plt.savefig('destination_path.eps', format='eps', dpi=1000)
b=0

minblocksize=20 # The minimum block size is 10ns=20 points
beg = 10 # We leave out the first 5ns=10 points
blocksNum=3 #the function best_start uses also 3 blocks for averaging

Nrows = 5
Ncolumns = 2
fig, ax = plt.subplots(Ncolumns, Nrows)
#fig.subplots_adjust(bottom=0.5,top=2.5, left=0.0, right=1.5, hspace = 0.5)
fig.subplots_adjust(bottom=0.0,top=4.5, left=0.0, right=1.5, hspace = 0.35)
matplotlib.rcParams['legend.handlelength'] = 0
matplotlib.rcParams['legend.markerscale'] = 0
matplotlib.rcParams['legend.handletextpad'] = 0
matplotlib.rcParams['legend.markerscale'] = 2 

#matplotlib.rcParams['font.family'] = 'Times New Roman'    
#matplotlib.rcParams['font.family'] = 'Times New Roman Bold' 

errortheta_m[b]=[]
theta_m[b]=[]
start_theta_m[b]=[]
end_theta_m[b]=[]
equil_theta_m[b]=[]

errortheta_s[b]=[]
theta_s[b]=[]
start_theta_s[b]=[]
end_theta_s[b]=[]
equil_theta_s[b]=[]

i = 0
for c in Waters:
#for c in [1000]:
    print c," molec:"
    n=Waters.index(c)

    theta2_w = array(angles_w[(b, c)])
    theta2_m = array(angles_m[(b, c)])
    theta2_s = array(angles_s[(b, c)])
    
    end=endpoint(theta2_m)
    start=best_start(theta2_m,t,beg,minblocksize)
    thetamean, errrorbar =blockAverage((theta2_m[start:end]),blocksNum)
    # SAVING RESULTS
    if start != None:
        (errortheta_m[b]).append(errrorbar)
        (theta_m[b]).append(thetamean)
        (start_theta_m[b]).append(start)
        (end_theta_m[b]).append(end)
        (equil_theta_m[b]).append(c)
        
    end=endpoint(theta2_s)
    start=best_start(theta2_s,t,beg,minblocksize)
    thetamean, errrorbar =blockAverage((theta2_s[start:end]),blocksNum)
    # SAVING RESULTS
    if start != None:
        (errortheta_s[b]).append(errrorbar)
        (theta_s[b]).append(thetamean)
        (start_theta_s[b]).append(start)
        (end_theta_s[b]).append(end)
        (equil_theta_s[b]).append(c)
    
    ax = plt.subplot(Nrows, Ncolumns, i+1)
    
    ####line1, = ax.plot(t, theta2_w,'-',linewidth=0.5)
    #ax.plot(t, theta2_w,'-',linewidth=0.5)
    #ax.plot(t, theta2_m,'-',linewidth=0.5)
    #ax.plot(t, theta2_s,'-',linewidth=0.5)
    
    mylines=ax.plot(t, theta2_w,'o', t, theta2_m, 'o', t,theta2_s,'o')
    plt.setp(mylines, markersize=3.0)
    
    matplotlib.rcParams['font.family'] = 'Arial'
    
    # Set common labels
    ax.set_xlabel('time [ns]',fontsize=13)
    plt.ylabel(r'$\theta_{mic}$ [deg]',fontsize=14)
    #ax.set_ylabel(r'R$\ _{BASE} $ [nm]')
    
    #titles for each subplot:
    matplotlib.rcParams['font.family'] = 'Times New Roman Bold' 
    #####ax.set_title(str(c)+' water molec: $\Delta(t)=$'+str(sampling)+'ns')
    ax.set_title(str(c)+' water molec.',fontweight='bold',fontsize=15)

    matplotlib.rcParams['font.family'] = 'Arial' 
    
    i = i+1
 
matplotlib.rcParams['font.family'] = 'Times New Roman' 
mytitle = plt.suptitle(r'Contact Angle $\theta_{mic}$ vs. Time for SAM'+str(b)+'%', fontsize=19, fontweight='bold',x=0.7, y=4.65)
###mytitle = plt.suptitle(r'Contact Angle $\theta_{mic}$ vs. Time for SAM'+str(b)+'%', fontsize=17,x=0.7, y=4.65)
plt.show()
fig.savefig('theta_t_s'+str(b)+'.jpg',bbox_extra_artists=(first_legend,mytitle,), bbox_inches='tight')
# check for better quality:
# plt.savefig('destination_path.eps', format='eps', dpi=1000)
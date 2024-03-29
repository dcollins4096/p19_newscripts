# UPDATED: AlphaBAnalysis.py. 
# use this version for ideas scattered throughout file
'''
 PART OF P19_NEWSCRIPTS REPOSITORY:
This script must be placed on: ~/p19_newscripts
With: starter1.py, starter2.py on same directory.
Started as: ~/p19_newscripts/tools_tracks/density_radius.py
In conjuction with: profiles66_0#.py

 notes: 
 for debug purposes, long_list = long_list[:3]
'''

# - - - - - - - - - - - - - - - - - - - - 
from starter2 import *
import data_locations as dl
import davetools
reload(davetools)

import scipy
from scipy import stats
from scipy.optimize import curve_fit
import matplotlib.gridspec as gridspec
import matplotlib.ticker as mtick
# - - - - - - - - - - - - - - - - - - - - 
# MOST LIKELY "MAIN"

plt.close('all')

# RUN ONCE with ipython -i file_name.py, then this if statement should save you time
# RUN SECOND time with run -i file_name.py in python shell

if 'this_simname' not in dir():
    this_simname = 'u203' 

def testforEleven(loop,location):
    if 11 in loop.tr.frames:
        print(':------(',location)

if 'this_looper' not in dir():
    file_list=glob.glob(dl.every_ten[this_simname])#[:2] 
    #file_list = ['/archive1/luzlourdes/u202/u202_all_primitives_primitives_c0106_nXXX0.h5']
    this_looper=looper.core_looper(directory=dl.sims[this_simname])
    for nfile,fname in enumerate(file_list):
        this_looper.load_loop(fname)
        #testforEleven(this_looper,fname)  #EXAMPLE TESTforELEVEN 
        print( "Reading file %d of %d"%(nfile,len(file_list)))
    thtr = this_looper.tr
    thtr.sort_time()
    all_cores = np.unique(thtr.core_ids)

print("HELLO")
# - - - - - - - - - - - - - - - - - - - - 
# PICK YOUR CORES
#core_list=all_cores
core_list = looper.get_all_nonzero(dl.n_particles[this_simname])  #CAREFUL, check looper.py!!
#core_list = [21,70,85,165,275,297]  #u05
#core_list = [70,165,275,297]  #u05
#core_list = [75,105,277,362]  #u202


# TRIED THINGS...
#cl, p, anz = looper.get_all_nonzero()  #change looper.py for profiles.py, UO5!
#core_part = parts[128:,0]
#particle_part = parts[128:,1]
frame_list = dl.frame_list[this_simname]


# GLOBAL TFF 
G = 1620/(4*np.pi)
rho_mean = 1
t_ff = np.sqrt(3*np.pi/(32*G*rho_mean)) 

# TFF PERCENTAGE:
tff_p = thtr.times[:2]/t_ff  #MUST be editted for suitable time steps 
tff_labels = ['%.2f'%s for s in tff_p]


# - - - - - - - - - - - - - - - - - - - - 
rm = rainbow_map(len(all_cores)) 
if 'rho_extents' not in dir():
    rho_extents=davetools.extents()
    magfield_extents = davetools.extents()
    vorticity_extents = davetools.extents() #ADDED
    r_extents=davetools.extents()

    for nc,core_id in enumerate(all_cores): 
        ms = trackage.mini_scrubber(thtr,core_id)
        if ms.nparticles == 1:
            continue
        density = thtr.c([core_id],'density')
        zeroes = density[0][0]
        if zeroes == 0:
            print(':------(',location)
        magfield = thtr.c([core_id],'magnetic_field_strength')
        magfield_x = thtr.c([core_id],'magnetic_field_z') #ADDED 
        rho_extents(density)
        magfield_extents(magfield)
        r_extents(ms.r)


# - - - - - - - - - - - - - - - - - - - - 
# ARRAYS FOR HISTOGRAMS
betarr = np.empty([0],dtype=float)
beta_neg_c = np.empty([0],dtype=float)
beta_neg_p = np.empty([0],dtype=float)
time_stamp = np.empty([0],dtype=float)

beta_curvarr = np.empty([0],dtype=float)
pearsonr = np.empty([0],dtype=float)
pearson_lmr = np.empty([0],dtype=float)
pearson_mr = np.empty([0],dtype=float)
spearmanr = np.empty([0],dtype=float)


# SHOULD CREATE A DICTIONARY...
pear0 = np.empty([0],dtype=float)
pear1 = np.empty([0],dtype=float)
pear2 = np.empty([0],dtype=float)
pear3 = np.empty([0],dtype=float)
pear4 = np.empty([0],dtype=float)
pear5 = np.empty([0],dtype=float)
pear6 = np.empty([0],dtype=float)
pear7 = np.empty([0],dtype=float)
pear8 = np.empty([0],dtype=float)
pear9 = np.empty([0],dtype=float)
pear10 = np.empty([0],dtype=float)
pear11 = np.empty([0],dtype=float)
pear12 = np.empty([0],dtype=float)
pear13 = np.empty([0],dtype=float)
pear14 = np.empty([0],dtype=float)

bear0 = np.empty([0],dtype=float)
bear1 = np.empty([0],dtype=float)
bear2 = np.empty([0],dtype=float)
bear3 = np.empty([0],dtype=float)
bear4 = np.empty([0],dtype=float)
bear5 = np.empty([0],dtype=float)
bear6 = np.empty([0],dtype=float)
bear7 = np.empty([0],dtype=float)
bear8 = np.empty([0],dtype=float)
bear9 = np.empty([0],dtype=float)
bear10 = np.empty([0],dtype=float)
bear11 = np.empty([0],dtype=float)
bear12 = np.empty([0],dtype=float)
bear13 = np.empty([0],dtype=float)
bear14 = np.empty([0],dtype=float)

# - - - - - - - - - - - - - - - - - - - - 
# LABELS FOR ALL PLOTS
def labelled(ax,xscale=None,yscale=None,xlabel=None,ylabel=None,\
             xlim=None,ylim=None,title=None,linthreshx=0.1,linthreshy=0.1):
    ax.set_xscale(xscale)
    ax.set_yscale(yscale)
    if xscale == 'symlog':
        ax.set_xscale(xscale,linthreshx=linthreshx)
    if yscale == 'symlog':
        ax.set_yscale(yscale,linthreshy=linthreshy)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)
    ax.set_title(title)


#~~~~~~~~~~ CORELOOP
for nc,core_id in enumerate(core_list):  #WATCH    
    #if core_id == 269:  #FOR U202
    #    continue
    #if core_id == 211:  #FOR U202
    #    continue
    #if core_id == 166:  #FOR U201
    #    continue 
    #if core_id == 34:  #FOR U201
    #    continue
    #if core_id == 14:  #FOR U201
    #    continue

    # miniscrubber computes distance, r^2, several other quantities 
    ms = trackage.mini_scrubber(thtr,core_id)
    tmap=rainbow_map(ms.ntimes)
    if ms.nparticles == 1:
        continue

    asort =  np.argsort(thtr.times)  #EDIT
    density = thtr.c([core_id],'density')
    magfield = thtr.c([core_id],'magnetic_field_strength')
    cellvolume = thtr.c([core_id],'cell_volume')
    if (asort != sorted(asort)).any():
        print("Warning: times not sorted.") 

    # SET-UP FOR PLOTTING!
    # FOR ALL TIMES, turn this "on"
    if 0 :
        fig, ax1=plt.subplots(1,1) 
    # FOR TWO PANEL SUBPLOTS, turn this "on"
    if 0:
        fig, (ax1,ax2)=plt.subplots(1,2)
    # FOR MULTIPLE PANELS, turn this "on" 
    if 0:
        fig = plt.figure()
        # exploring a common axis with fig.text OR ax = fig.add_subplot(111) 
        fig.text(0.365,0.03,r'$\rho/\rho_{o}$')#, ha='center', va='center') 

        ax1 = plt.subplot(331)
        ax2 = plt.subplot(332)
        ax4 = plt.subplot(334)
        ax5 = plt.subplot(335)
        ax7 = plt.subplot(337)
        ax8 = plt.subplot(338)

        ax3 = plt.subplot(133)  #`ax6, ax9` 
        #lplots = [0,ax1,0,ax2,0,ax4,0,0,ax5,0,ax7,0,ax8,0]  #for u05  
        #lplots = [0,ax1,0,ax2,0,ax4,0,ax5,0,ax7,0,ax8,0]  #for u202  
        lplots = [0,ax1,0,ax2,0,ax4,ax5,0,ax7,0,ax8,0]  #for u203 
        fig.subplots_adjust(wspace=0, hspace=0)


#~~~~~~~~~~ TIMELOOP
    for n_count,n_time in enumerate(asort):  
        #print("TIME begin",n_time)
        time=thtr.times[n_time]  #these next three lines for ALL TIME only
        if time == 0:
            continue
        c=tmap(n_count,ms.nparticles) # previously as n_count
      
        this_r=ms.r[:,n_time]+0 
        this_B = magfield[:,n_time]+0  #no mini_scrubber needed
        this_rho=ms.density[:,n_time]+0 
        this_cv=ms.cell_volume[:,n_time]+0 

        r_un = nar(sorted(np.unique(this_r)))

        # FOR ALL TIMES do field only, choose respective fig, and change lplots[n_time] for ax1
        # FOR EACH TIME FRAME, add [:,n_time]
 
        X = np.log10(density).flatten()   
        Y = np.log10(magfield).flatten()  
        #Z = np.log10(cellvolume[:,n_time]).flatten()  #ADDED 

        XX = 10 ** X 
        X2 = np.linspace(X.min()+2,X.max()-3,num=len(X))  # FOR ALL TIME, ONE PANEL PLOTS 
        #X2 = np.linspace(X.min(),X.max(),num=len(X)) # FOR PER FRAME, MULTIPLE PANEL PLOTS 
        XX2 = 10 ** X2 
   
 
 
        # NOTE: the "else" statements are a holder for the pearsonR's that are nans;
        # these need to be accounted for rightfully. suggestion: reverse the time,core loops
        if 0:
        #~~~~~~~~~~ PEARSONr: SCIPY vs MANUAL 
            xs = np.std(X)
            mxs = np.std(this_rho)
            ys = np.std(Y)
            mys = np.std(this_B)

        # SCIPY PEARSONr: scipy logged
            if xs != 0 and ys != 0:
                pearX,pearY = scipy.stats.pearsonr(X,Y)  
                #print("TIME",n_time)
            else:
                print("A zero encountered!!",xs,ys)
                pearX = 0  
            pearsonr = np.append(pearsonr,pearX)  # REMIND SELF HOW TO FIX THIS OR HOW THIS MAY BE OK

        # SCIPY PEARSONr: scipy not logged
            #if mxs != 0  and mys != 0:
            #    pear0,pear1 = scipy.stats.pearsonr(this_rho,this_B)  
            #else:
            #    pear0 = 0  
            #pearsonr = np.append(pearsonr,pear0)  #otherwise pear0 is taking the same value twice
     
        # MANUAL PEARSONr: not logged 
            #avB = np.sum(this_B*this_cv)/np.sum(this_cv)
            #varB = np.sum((this_B-avB)**2*this_cv)/np.sum(this_cv)
            #avRho = np.sum(this_rho*this_cv)/np.sum(this_cv)
            #varRho = np.sum((this_rho-avRho)**2*this_cv)/np.sum(this_cv)
            #covar = np.sum((this_B-avB)*(this_rho-avRho)*this_cv)/np.sum(this_cv)
            #if varB != 0 and varRho != 0:
            #    pearson_r = covar/(np.sqrt(varB)*np.sqrt(varRho)) 
            #else:
            #    pearson_r = 0
            #pearson_mr = np.append(pearson_mr,pearson_r) 

        # MANUAL PEARSONr: logged parameters
            #l_avB = np.sum(Y * Z) / np.sum(Z)
            #l_varB = np.sum((Y - l_avB)**2 * Z) / np.sum(Z)
            #l_avRho = np.sum(X * Z) / np.sum(Z)
            #l_varRho = np.sum((X - l_avRho)**2 * Z) / np.sum(Z)
            #l_covar = np.sum((Y-l_avB) * (X-l_avRho) * Z) / np.sum(Z)
            ##if l_varB <=0  or l_varRho <= 0:
            #if l_varB != 0 and l_varRho != 0:  # TEST
            #    pearsonR = l_covar / (np.sqrt(l_varB) * np.sqrt(l_varRho))  
            #    #if np.isnan(pearsonR):  # or if a != a ...true = problem
            #    #    pdb.set_trace()
            #else:
            #    pearsonR = 0 
            #pearson_lmr = np.append(pearson_lmr,pearsonR)


        # NOTES 
        # B = B0 rho^beta --> lnB = ln(B_0*rho^beta) -->
        # lnB = beta*ln(rho) + lnB_0

        # CURVE-FIT:  -revise!
        #def testing(x, a, b):
        #    return a*x+b 
        #param, param_cov = curve_fit(testing, X, Y)
        #Y_B = 10 ** (param[0]*X + param[1])
        #beta_curve = param[0]

        # POLY-FITTING:
        #print('TESTING_2') 
        pfit = np.polyfit(X,Y,1) 
        other = np.poly1d(pfit)
        beta = pfit[0]
        B_o = pfit[1]
       

        if n_time == asort[-1]:  #FOR ALL TIME ONLY
            betarr = np.append(betarr,beta)  #unindent if per frame
        #print('TESTING_3') 
        #print(betarr) 
        #beta_curvarr = np.append(beta_curvarr, beta_curve)

        # THE NEGATIVE OUTLIERS
        if beta < 0:
            #comment out when using profiles.py
            #beta_neg_p = np.append(beta_neg_p,particle_part[nc])  
            #beta_neg_c = np.append(beta_neg_c,core_part[nc])  
            time_stamp = np.append(time_stamp,n_time)

#1,3,5,8,10,12: u05
#1,3,5,7,9,11:  u202
        YY = 10 ** (pfit[0]*X2 + pfit[1]) 
        if 0:  
            #if n_time in {1,3,5,6,8,10}:  # EDIT for diff sims 
                # FOR ALL TIME indent 4 left, comment if above, change lplots[n_time] to ax1, & comment out labelled to ax8
            ax1.scatter(density[:,n_time],magfield[:,n_time],c=c,label=thtr.times[n_time],s=0.1)  #EDIT thtr.times          
            ax1.plot(XX2,YY,c='k',linewidth=1.0) #c=c[0] for colors
            print("plotted")             

                #labelled(lplots[n_time],xscale='log',yscale='log',xlabel=None,ylabel=None,
                #         xlim=rho_extents.minmax, ylim=magfield_extents.minmax)
                #ax2.tick_params(axis='y',labelleft=False)
                #ax4.set_ylabel(r'$\mid B \mid (\mu G)$')
                #ax5.tick_params(axis='y',labelleft=False)
                #ax8.tick_params(axis='y',labelleft=False)
                
            # ADD LINEAR FITS FOR CONTRAST
            if 0: 
                #YY_four = 10 ** ((2/5)*X2 + pfit[1])
                #YY_five = 10 ** ((1/2)*X2 + pfit[1])
                #YY_six = 10 ** ((2/3)*X2 + pfit[1]) 
          
                #ax1.plot(XX2,YY_four,c='r',linewidth=0.5)
                #ax1.plot(XX2,YY_five,c='b',linewidth=0.5)
                #ax1.plot(XX2,YY_six,c='g',linewidth=0.5) 
    
                labelled(ax1,xscale='log',yscale='log',xlabel=r'$\rho/\rho_{o}$',ylabel=r'$\mid B \mid (\mu G)$',
                         xlim=rho_extents.minmax, ylim=magfield_extents.minmax,title=r'$\alpha_{B\rho} = %.3f$'%beta)

            outname = 'BrhoTff_c%04d'%(core_id) 
            if 0:  #WHEN should I turn this on?
                if n_time == asort[-1]:
                    plt.savefig(outname)
                    print("saved "+outname)
                    plt.close(fig) 
 
        # FOR HISTOGRAMS PER FRAME, SCATTER BETA, AND BOX PLOTS
        if 0:
            if n_time == asort[-1]:              
                #pear0 = np.append(pear0,pearsonr[0])
                #pear1 = np.append(pear1,pearsonr[1]) 
                #pear2 = np.append(pear2,pearsonr[2]) 
                #pear3 = np.append(pear3,pearsonr[3]) 
                #pear4 = np.append(pear4,pearsonr[4]) 
                #pear5 = np.append(pear5,pearsonr[5]) 
                #pear6 = np.append(pear6,pearsonr[6])  
                #pear7 = np.append(pear7,pearsonr[7]) 
                #pear8 = np.append(pear8,pearsonr[8]) 
                #pear9 = np.append(pear9,pearsonr[9]) 
                #pear10 = np.append(pear10,pearsonr[10]) 
                #pear11 = np.append(pear11,pearsonr[11])  #final u203 pear contains nans 
                #pear12 = np.append(pear12,pearsonr[12])  #final u202 pear, but want to discard final frame   
                #pear13 = np.append(pear13,pearsonr[13]) 
                #pear14 = np.append(pear14,pearsonr[14]) 
              
                bear0 = np.append(bear0,betarr[0])
                bear1 = np.append(bear1,betarr[1]) 
                bear2 = np.append(bear2,betarr[2]) 
                bear3 = np.append(bear3,betarr[3]) 
                bear4 = np.append(bear4,betarr[4]) 
                bear5 = np.append(bear5,betarr[5]) 
                bear6 = np.append(bear6,betarr[6])  
                bear7 = np.append(bear7,betarr[7]) 
                bear8 = np.append(bear8,betarr[8]) 
                bear9 = np.append(bear9,betarr[9])
                bear10 = np.append(bear10,betarr[10])  #for u203
                bear11 = np.append(bear11,betarr[11])  #for u202 
                #bear12 = np.append(bear12,betarr[12])    
                #bear13 = np.append(bear13,betarr[13])
                #bear14 = np.append(bear14,betarr[14])

                # COMMENT OFF/ON AS NEEDED
                #pearson_lmr = np.empty([0],dtype=float)
                #pearsonr = np.empty([0],dtype=float) 
                #betarr = np.empty([0],dtype=float)
                #beta_curvarr = np.empty([0],dtype=float)

                if 0: 
                    #plt.clf()
                    tmap2 = rainbow_map(len(betarr)) 
                    c2 = [tmap2(n) for n in range(len(betarr))] 
                    #(np.arange(0,1,0.075)  #u05
                    #(np.arange(0.075,1,0.075)  #u202
                    ax3.scatter(np.arange(0.075,0.9,0.075),betarr,c=c2)  
                    ax3.plot(np.arange(0.075,0.9,0.075),betarr,c='k',linewidth=1.0)  
                    ax3.plot([0,1],[0,0],c=[0.5]*4)
                 
                    tff_mod = [0.0,0.0,0.3,0.6,0.9]  #this should match MultipleLocator
                    ax3.set_xticklabels(tff_mod)   
                    ax3.xaxis.set_major_locator(plt.MultipleLocator(0.3))

                    ax3.yaxis.tick_right()
                    ax3.set_xlabel(r'$t_{\rm{ff}}$')  #\rm{  },to avoid italization 
                    ax3.set_yscale('linear')
                    ax3.set_ylim(-0.5,0.5)     
                    ax3.set_title(r"$\alpha_{B}$")#, P=%.3f, S=%.3f"%(p1,s1)) 
          
                    #spearmanr = np.empty([0],dtype=float)
                    betarr = np.empty([0],dtype=float)

    if 0:
        #plt.tight_layout()
        fig.savefig(outname)
        print("saved "+outname)
        plt.close(fig)

pears = [pear0,pear1,pear2,pear3,pear4,pear5,pear6,\
         pear7,pear8,pear9,pear10,pear11,pear12,pear13,pear14]
bears = [bear0,bear1,bear2,bear3,bear4,bear5,bear6,\
         bear7,bear8,bear9,bear10,bear11,bear12,bear13,bear14]

# - - - - - SCATTER_PLOT FOR ALL FRAMES - revise
if 0:
    newbetas = np.empty([0],dtype=float)
    for i in range(13):
        newbetas = np.append(newbetas,betas[i])
    plt.plot([0,13],[0,0],c=[0.5]*4)
    # 175 cores * 14 frames = 2450. 14/2450 = 0.00571429
    plt.scatter(np.arange(0,13,0.00571429),newbetas)
    plt.ylim(-1.0,1.0)
    plt.xlabel('Time')
    plt.ylabel('Beta')
    plt.savefig('AllFramesBetas')

# - - - - - BOXPLOT FOR ALL FRAMES
# NOTE: for PearsonR, nan values are not accepted
print("BEFORE BOXVIOPLOT")
if 0:  
    count = 0 #temp
    Pears = {}
    index = []
    data = []
    #tsorted = tsorted[2:] 
    for i,num in enumerate(tff_labels):         
        Pears[tff_labels[i]] = pears[i]  #TEST, previously "frames" 
    for j, (key, val) in enumerate(Pears.items()):
    #    for one, two in enumerate(val):  # DOUBLE CHECK if this does what I want, and if necessary, maybe 0 vals are already ignored
    #        print(two)
    #        count = count + 1
    #        if two == 0:
    #            continue
        #if val == 0:  #EDIT: val is an entire array, sort through each value
        #gives ValueError: The truth value of an array with more than one element is ambiguous. Use a.any() or a.all() 
        index.append(key)
        data.append(val)

    fig, ax1=plt.subplots(1,1)    

    #To compare boxplot or violinplots with a zero, however I think this is affecting the position of the plots
    ax1.plot([1,14],[0,0],c=[0.5]*4)  #U05  
    #ax1.plot([1,11],[0,0],c=[0.5]*4)  #U203 
    #ax1.plot([1,12],[0,0],c=[0.5]*4)  #U202

    # trying to find how to use 'position' suitingly as to match time...
    if 0:
        bparts = ax1.boxplot(data,showfliers=True,showmeans=True,meanline=True)  #bp is now a dictionary     
    if 1:
        vparts = ax1.violinplot(data,showmeans=True, showextrema=True, showmedians=True)
        vparts['cmeans'].set_color('r')   
    tff_mod = [0.0,0.0,0.01,0.07,0.15,0.23,0.30,0.38,0.45,0.53,0.61,0.68,0.75,0.82,0.89]#,0.95]  #U05 add a 0.0 
    #tff_mod = [0.0,0.0,0.07,0.15,0.23,0.30,0.38,0.45,0.53,0.61,0.68,0.75]  #U201 
    #tff_mod = [0.0,0.0,0.08,0.17,0.26,0.35,0.44,0.52,0.60,0.69,0.77,0.86]#,0.91]  #U203
    #tff_mod = [0.0,0.0,0.08,0.16,0.25,0.33,0.41,0.50,0.58,0.66,0.74,0.82,0.90]#,0.97]  #U202
    
    #ax1.xaxis.set_major_locator(plt.MultipleLocator(2))  #U05 order swapped with set_xticklabels
    ax1.xaxis.set_major_locator(plt.MultipleLocator(1))  #U10,U11

    ax1.set_xticklabels(tff_mod)     
    ax1.set_xlabel(r'$t_{\rm{ff}}$')
    #ax1.set_ylim(-1.10,1.10)  #EDIT: pearsonR
    ax1.set_ylim(-1.0,1.0)
    #ax1.set_ylim(-1.25,1.25) 
    #ax1.set_ylabel(r'$\alpha_B$')
    ax1.set_ylabel(r'$r_{B\rho}$')
    ax1.set_title(r'$\beta = 0.2$')  # EDIT
    fig.savefig('PearsBoxplotTen_%s'%this_simname)
    print("SAVED")
    plt.close(fig)

# - - - - - HISTOGRAM FOR ALL CORES FOR ALL TIME
print("before histos!!")  #TURN ON 4 SWTICHES FOR ALL CORE TIMES
if 0:
    #BETA ENTRIES, AVERAGE & STD 
    betavg = np.mean(betarr)
    betastd = np.std(betarr)
    entries = len(betarr) 
  
    fig, ax1 = plt.subplots()

    ax1.hist(betarr, 50, density=False, histtype='step', color='r')  #change color for diff sims
    ax1.set_xlabel(r'$\alpha_{B\rho}$') 
    ax1.set_ylabel('PDF')

    y_vals = ax1.get_yticks()
    ax1.set_yticklabels(['{:.3f}'.format(x/len(betarr)) for x in y_vals])

    name_save = "BetaHistogramTff_u202" 

    t = 'Cores: %d\n'%entries +  r'Mean $\beta = %.3f$'%betavg + '\n' + r'$\sigma = %.3f$'%betastd
    ax1.text(0.56, 15, t, color='green', bbox=dict(facecolor='none', edgecolor='black', boxstyle='round,pad=1'))  #u05 
    #ax1.text(-0.3, 3, t, color='green', bbox=dict(facecolor='none', edgecolor='black', boxstyle='round,pad=1'))  #u101 
    #ax1.text(0.47, 13.5, t, color='blue', bbox=dict(facecolor='none', edgecolor='black', boxstyle='round,pad=1'))  #u10 
    #ax1.text(0.43, 17.5, t, color='blue', bbox=dict(facecolor='none', edgecolor='black', boxstyle='round,pad=1'))  #u102 
    #ax1.text(0.50, 12.5, t, color='m', bbox=dict(facecolor='none', edgecolor='black', boxstyle='round,pad=1'))  #u11 
    #ax1.text(0.48,11, t, color='m', bbox=dict(facecolor='none', edgecolor='black', boxstyle='round,pad=1'))  #u103 
    plt.savefig(name_save)
    print('saved '+name_save)
   
# - - - - - HISTOGRAM FOR ALL CORES PER TIME FRAME   
if 0: 
    for i,vals in enumerate(bears): 
        i = i+2 
        fig, ax1 = plt.subplots()
        
        ax1.hist(vals, 50, density=False, histtype='step', color='g')
        ax1.set_xlim(-1.5,1.5)  #best to find max and min of betas.
        ax1.set_xlabel(r'$\beta$')
        
        ax1.set_ylabel('Ocurrence')
        y_vals = ax1.get_yticks()
        ax1.set_yticklabels(['{:.3f}'.format(x/len(vals)) for x in y_vals])
        
        ax1.set_title('Frame %d'%i)   
        outname = 'u05cut_Bear_histogram%d'%i
        plt.savefig(outname)
        print("saved "+outname) 
        plt.clf()  #OR? plt.close(fig)


# - - - - - SCATTTER PLOT FOR PROFILES66.PY 
def plot_particles(axes,nframe):  
    for core_id in all_cores:
        ms = trackage.mini_scrubber(thtr,core_id)
        tmap=rainbow_map(ms.ntimes)
        if ms.nparticles == 1:
            continue
        asort =  np.argsort(thtr.times[1:])
        density = thtr.c([core_id],'density')
        magfield = thtr.c([core_id],'magnetic_field_strength')
        vorticity = thtr.c([core_id],'vorticity_magnitude')  #ADDED

        rho_extents(density)  #ADDED
        magfield_extents(magfield)  #ADDED
        vorticity_extents(vorticity)  #ADDED 

        vortRho = vorticity/density
        magRho = magfield/density

        if (asort != sorted(asort)).any():
            print("Warning: times not sorted.")   
       
        for n_count,n_time in enumerate(asort): 
            #time=thtr.times[n_time]
            #if time == 0:
            #    continue
            c=tmap(n_count,ms.nparticles)
            this_r=ms.r[:,n_time]+0
            r_un = nar(sorted(np.unique(this_r)))
        
            X = np.log10(density[:,n_time]).flatten()   
            Y = np.log10(magfield[:,n_time]).flatten()  
            
            if n_time == nframe: 
                axes.scatter(density[:,n_time],magfield[:,n_time],c='k',label=thtr.times[n_time],s=0.1)  #EDIT possibly thtr.times           
    return

# - - - - - SCATTER PLOT FOR <beta> vs time: EDIT
if 0:
    def avgparams(bears):
        return np.mean(bears)
    averaged = map(avgparams, bears) 

# maybe try this later....
#    def squared(bears):
#        for i in len(bears):
#            squared = map(square, bears[i])
#        return np.mean(squared)

    plt.plot([0,13],[0,0],c=[0.5]*4)
    plt.scatter(tff_labels,list(averaged))
    plt.ylim(-0.5,0.5)
    plt.xlabel('Time')
    plt.ylabel('Beta Average')
    plt.savefig('BetaAvg_Time')


# - - - - - FROM MASS_TIME.PY, Q: ARE WE DOUBLE COUNTING ZONES?
dx = 1./2048
nx = 1./dx
'''
for nc,core_id in enumerate(core_list):

    for nf,frame in enumerate(thtr.frames):
        x =np.floor(thtr.c([core_id],'x')/dx)[:,nf]
        y =np.floor(thtr.c([core_id],'y')/dx)[:,nf]
        z =np.floor(thtr.c([core_id],'z')/dx)[:,nf]
        density = thtr.c([core_id],'density')[:,nf]
        cell_volume = thtr.c([core_id],'cell_volume')[:,nf]
        index = x + nx*(y * nx*z)
        ar = np.argsort(index)
        rs = np.argsort(ar)
        isorted=index[ar]
        mask = np.ones_like(density,dtype='bool')
        mask[1:] = isorted[1:]-isorted[:-1] != 0
        mask2 = mask[ rs]
        mass = (density[mask2]*cell_volume[mask2]).sum()
'''

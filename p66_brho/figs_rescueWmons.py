
'''
needs to be in the same folder as tsung_spheres
'''

# edit
from starter2 import *
import davetools
reload(davetools)
import p49_fields
reload(p49_fields)
import math
import pcolormesh_helper as pch
from matplotlib.ticker import PercentFormatter
np.set_printoptions(threshold=sys.maxsize)

import track_loader as TL
import cfpack as cfp
from cfpack import stop,print
import lmfit
import latex
import matplotlib as mpl
mpl.rcParams.update(mpl.rcParamsDefault)
import seaborn as sns

import tsing
reload(tsing)
import tsung_spheres
reload(tsung_spheres)
import monster
reload(monster)
# --- --- --- --- --- --- ---


class withspheres(): 
    def __init__(self,the_loop):
        self.this_looper = the_loop

        self.bmag_sph = defaultdict(list)
        self.rhoave_sph = defaultdict(list)
        self.bmag_sph_rinf = defaultdict(list)
        self.rhoave_sph_rinf = defaultdict(list)

        self.bmag_parts = defaultdict(list)
        self.rhoave_parts = defaultdict(list)

        self.ncolumn_sph = defaultdict(list)
        self.blos_sph = defaultdict(list)
        self.ncolumn_sph_rinf = defaultdict(list)
        self.blos_sph_rinf = defaultdict(list)
        
        self.the_rinfs = defaultdict(list)

    def framescores(self, sim, typeofstudy='withspheres', core_list=None): 
        print('inside!!')
        thtr = self.this_looper.tr

        # CORES
        all_cores = np.unique(thtr.core_ids)
        if core_list is None:
            core_list = all_cores  #or debug

        # EVERY TEN FRAMES
        for nf,frame in enumerate(thtr.frames[1:]):  #or debug 

            # CORE-LOOP
            for nc,core_id in enumerate(core_list):
                ms = trackage.mini_scrubber(thtr,core_id,do_velocity=True)
                monster.load([sim])
                the_monster = monster.closet[sim]

                the_radius_rinf = the_monster.get_r_inflection(core_id,frame)
                self.the_rinfs[nf].append(the_radius_rinf)
                if the_radius_rinf < 0.0078125:  #which is 1/128
                    print('radius smaller than the_radius_rinf',the_radius_rinf)
                #    continue
                the_sphere_rinf = the_monster.get_sphere(core_id,frame,'rinf')
                the_sphere = the_monster.get_sphere(core_id,frame,'r1') 

                # get and store the avg fields
                # WITH SPHERES
                self.bmag_sph[nf].append((the_sphere['density']* the_sphere['magnetic_field_strength'] * the_sphere['cell_volume']).sum()/the_sphere['cell_mass'].sum()) 
                self.rhoave_sph[nf].append((the_sphere['density'] * the_sphere['cell_volume']).sum()/ the_sphere['cell_volume'].sum())
                self.bmag_sph_rinf[nf].append((the_sphere_rinf['density']* the_sphere_rinf['magnetic_field_strength'] * \
                                               the_sphere_rinf['cell_volume']).sum()/the_sphere_rinf['cell_mass'].sum()) 
                self.rhoave_sph_rinf[nf].append((the_sphere_rinf['density'] * the_sphere_rinf['cell_volume']).sum()/ the_sphere_rinf['cell_volume'].sum())

                # WITH PARTICLES 
                if typeofstudy == 'particles':
                    mask = ms.compute_unique_mask(core_id, dx=1./2048,frame=nf)
                    # get the fields.
                    density = thtr.c([core_id],'density')[mask,nf]
                    cell_volume = thtr.c([core_id],'cell_volume')[mask,nf]
                    cell_mass = density * cell_volume  
                    bx = thtr.c([core_id],'magnetic_field_x')[mask,nf]
                    by = thtr.c([core_id],'magnetic_field_y')[mask,nf]
                    bz = thtr.c([core_id],'magnetic_field_z')[mask,nf]                
                    bb = np.sqrt(bx*bx+by*by+bz*bz) 
                    # store the avgs
                    self.bmag_parts[nf].append((bb * density * cell_volume).sum()/cell_mass.sum())  
                    self.rhoave_parts[nf].append((density * cell_volume).sum()/(cell_volume.sum()))  
     
        if typeofstudy == 'particles':
            data_bmagparts = [*self.bmag_parts.values()] 
            data_rhoparts = [*self.rhoave_parts.values()]
        data_bmagsph = [*self.bmag_sph.values()] 
        data_rhosph = [*self.rhoave_sph.values()]
        data_bmagsph_rinf = [*self.bmag_sph_rinf.values()] 
        data_rhosph_rinf = [*self.rhoave_sph_rinf.values()]
        data_the_rinfs = [*self.the_rinfs.values()]

        print('before saving to h5 files!!')
        if 1: 
            hfivename = 'p66_brho/brho_sph_r1rinf_%s.h5'%(sim)
            Fptr = h5py.File(hfivename,'w')
            Fptr['bfield_sph'] = data_bmagsph 
            Fptr['rhoavg_sph'] = data_rhosph
            Fptr['bfield_sph_rinf'] = data_bmagsph_rinf 
            Fptr['rhoavg_sph_rinf'] = data_rhosph_rinf
            Fptr['the_rinfs'] = data_the_rinfs 
            if typeofstudy == 'particles':
                Fptr['bfield_parts'] = data_bmagparts 
                Fptr['rhoavg_parst'] = data_rhoparts
            Fptr.close()
            print('h5 file written. closing.')


    def syntheticobs(self,sim,core_list=None): 
        print('inside synthetic observations!')
        thtr = self.this_looper.tr

        # CORES
        all_cores = np.unique(thtr.core_ids)
        if core_list is None:
            core_list = all_cores #or debug

        # EVERY TEN FRAMES 
        for nf,frame in enumerate(thtr.frames[1:]): 

            # CORE-LOOP
            for nc,core_id in enumerate(core_list):
                monster.load([sim])
                the_monster = monster.closet[sim]
                the_sphere = the_monster.get_sphere(core_id,frame,'r1') 
                the_area = np.pi * (1/128)**2
                the_radius_rinf = the_monster.get_r_inflection(core_id,frame)
                the_sphere_rinf = the_monster.get_sphere(core_id,frame,'rinf') 
                the_area_rinf = np.pi * the_radius_rinf**2
                
                B = ['magnetic_field_x','magnetic_field_y','magnetic_field_z']
                for j in range(3): 
                    self.ncolumn_sph[nf].append((the_sphere['density'] * the_sphere['cell_volume']).sum()/the_area)
                    self.blos_sph[nf].append((the_sphere['density'] * the_sphere[B[j]] * the_sphere['cell_volume']).sum()/the_sphere['gas','cell_mass'].sum())
 
                    self.ncolumn_sph_rinf[nf].append((the_sphere_rinf['density'] * the_sphere_rinf['cell_volume']).sum()/the_area_rinf)
                    self.blos_sph_rinf[nf].append((the_sphere_rinf['density'] * the_sphere_rinf[B[j]] * the_sphere_rinf['cell_volume']).sum()/the_sphere_rinf['gas','cell_mass'].sum())
 
        data_blossph = [*self.blos_sph.values()] 
        data_ncolumnsph = [*self.ncolumn_sph.values()]
        data_blossph_rinf = [*self.blos_sph_rinf.values()] 
        data_ncolumnsph_rinf = [*self.ncolumn_sph_rinf.values()]
        print('before saving to h5 files')
        if 1: 
            hfivename = 'p66_brho/blosncol_sph_r1rinf_%s.h5'%(sim)
            Fptr = h5py.File(hfivename,'w')
            #pdb.set_trace()
            Fptr['blos_sph'] = data_blossph 
            Fptr['ncolumn_sph'] = data_ncolumnsph
            Fptr['blos_sph_rinf'] = data_blossph_rinf 
            Fptr['ncolumn_sph_rinf'] = data_ncolumnsph_rinf
            Fptr.close()
            print('h5 file written. closing.')




# YOU ENTER HERE
# TO GET DATA AND STORE
# COMPARING SPHERES WITH PARTICLES
# note: take a look at p19_play/psedo.py for a brief example with monster.py
if 0:
    sims=['u601']#, 'u602','u603']
    TL.load_tracks(sims)
    for sim in sims:
        core_list=None
        running = withspheres(TL.loops[sim])
        if 1:
            running.framescores(sim)
        if 1:
            running.syntheticobs(sim)
# SPHERES SYNCED TO TSUNG
G = 1620/(4*np.pi)
rho_mean = 1
t_ff = np.sqrt(3*np.pi/(32*G*rho_mean))
if 0:
    sims=['u501', 'u502', 'u503']
    import three_loopers_u500 as TL   #EDIT THIS!! and put pdb.set_trace() back in this file
    #TL.load_tracks(sims)
    if 'tsing_tool' not in dir():
        tsing_tool={}
        for ns,sim in enumerate(sims):
            obj=tsing.te_tc(TL.loops[sim])
            tsing_tool[sim]=obj
            tsing_tool[sim].run()
    if 'mp' not in dir():
        for sim in sims:
            all_cores=np.unique(TL.loops[sim].tr.core_ids)
            core_list=list(all_cores)
            #core_list=core_list[:1]  #debug

            mp=tsung_spheres.tsungspheres(TL.loops[sim])   #EDIT! make tsung part of this class, then split THIS file into two respectively
            mp.run(core_list=core_list, tsing=tsing_tool[sim], obs=True, whichdataset='tsungtimes')



# TO READ DATA FROM STORAGE
if 1:   
    figtype = 'tsungtff'  #kappaperframe, kappatff, rinfradii: histos, tsungtff 
    individually = 'no'  #if kappaperframe, yes: one panel per frame, no: frame time series
    parts_or_spheres ='synthetic' #parts, spheres, sphparts, sph_tsung, synthetic  
    kappadyn = 'yes' #no: do ratios, yes: do kappa dynamical; EDIT: outnames should reflect this too
    whichradius = None#'rinf' 
    
    series = '500'  #500 or 600
    # for tsung plots
    if series == '500':  
        sims=['u503']#, 'u502', 'u503']  #EDIT: need a loop!!
        if figtype!= 'tsungtff':
            hfivename = 'p66_brho/h5files/brho_sphtsung_%s.h5'%(sims[0])  
            hfivename_synth = 'p66_brho/h5files/blosncol_sphtsung_%s.h5'%(sims[0])  
            Fptr = h5py.File(hfivename,'r')
            Fptr_synth = h5py.File(hfivename_synth,'r')

            #if parts_or_spheres == 'spheres':
            b_sph = Fptr['bmag_sph'][()] 
            rho_sph = Fptr['rho_sph'][()]
            if parts_or_spheres == 'synthetic':
                b_sph_synth = Fptr_synth['bmag_sph'][()] 
                rho_sph_synth = Fptr_synth['rho_sph'][()]
        if figtype == 'tsungtff':
            hfivename = 'p66_brho/h5files/coreids_tsungtff_%s.h5'%(sims[0])  
            Fptr = h5py.File(hfivename,'r')
            tsung_tff = Fptr['tsungtff'][()] 


    # for tff plots
    if series == '600':
        sims=['u601']#, 'u602', 'u603']  #EDIT: need a loop!!
        if whichradius == 'rinf':
            hfivename = 'p66_brho/h5files/brho_sph_r1rinf_%s.h5'%(sims[0])
            hfivename_synth = 'p66_brho/h5files/blosncol_sph_r1rinf_%s.h5'%(sims[0])
        else:
            hfivename = 'p66_brho/h5files/brho_sphparts_%s.h5'%(sims[0])  
            hfivename_synth = 'p66_brho/h5files/blosncol_sph_%s.h5'%(sims[0])  
        Fptr = h5py.File(hfivename,'r')
        Fptr_synth = h5py.File(hfivename_synth,'r')

        if figtype == 'kappaperframe' or figtype == 'kappatff':
            b_sph = Fptr['bfield_sph'][()] 
            rho_sph = Fptr['rhoavg_sph'][()] 
            if whichradius == 'rinf':
                b_sph_rinf = Fptr['bfield_sph_rinf'][()]  #added 
                rho_sph_rinf = Fptr['rhoavg_sph_rinf'][()]#added 
            if parts_or_spheres == 'parts':
                b_parts = Fptr['bfield_parts'][()] 
                rho_parts = Fptr['rhoavg_parst'][()] 
            if parts_or_spheres == 'synthetic':
                b_sph_synth = Fptr_synth['blos_sph'][()] 
                rho_sph_synth = Fptr_synth['ncolumn_sph'][()] 
                if whichradius == 'rinf':
                    b_sph_synth_rinf = Fptr_synth['blos_sph_rinf'][()]     #added 
                    rho_sph_synth_rinf = Fptr_synth['ncolumn_sph_rinf'][()]#added 
        if figtype == 'rinfradii':
            rinf_radii = Fptr['the_rinfs'][()]  
 


    def afunct(x, a, b):
        y = a * x + b
        return y

    if figtype == 'kappaperframe':
        kappas_sph = []
        kappas_sph_synth = []
        kappas_parts = []
        kappas_sph_rinf = []
        kappas_sph_rinf_synth = []
        for i in range(len(rho_sph)):   
            # relics of ignoring the 0's!
            #pre_rhosph = rho_sph[i][rho_sph[i] != 0]
            #rhosph_log = np.log10(pre_rhosph)  
            rhosph_log = np.log10(rho_sph[i])
            bsph_log = np.log10(abs(b_sph[i]))
            rets_sph = cfp.fit(afunct, rhosph_log, bsph_log)  
            kappas_sph = np.append(kappas_sph, rets_sph.popt[0])
            if whichradius == 'rinf':
                rhosph_rinf_log = np.log10(rho_sph_rinf[i])
                bsph_rinf_log = np.log10(b_sph_rinf[i])
                rets_sph_rinf = cfp.fit(afunct, rhosph_rinf_log, bsph_rinf_log)  
                kappas_sph_rinf = np.append(kappas_sph_rinf, rets_sph_rinf.popt[0])

            if parts_or_spheres == 'synthetic':
                rhosph_synth_log = np.log10(rho_sph_synth[i])
                bsph_synth_log = np.log10(abs(b_sph_synth[i]))
                rets_sph_synth = cfp.fit(afunct, rhosph_synth_log, bsph_synth_log)  
                kappas_sph_synth = np.append(kappas_sph_synth, rets_sph_synth.popt[0])
                if whichradius == 'rinf':
                    rhosph_synth_rinf_log = np.log10(rho_sph_synth_rinf[i])
                    bsph_synth_rinf_log = np.log10(abs(b_sph_synth_rinf[i]))
                    rets_sph_synth_rinf = cfp.fit(afunct, rhosph_synth_rinf_log, bsph_synth_rinf_log)  
                    kappas_sph_rinf_synth = np.append(kappas_sph_rinf_synth, rets_sph_synth_rinf.popt[0])

            if parts_or_spheres == 'parts':
                bparts_log = np.log10(b_parts[i]) 
                rhoparts_log = np.log10(rho_parts[i])  
                rets_parts = cfp.fit(afunct, rhoparts_log, bparts_log)  
                kappas_parts = np.append(kappas_parts, rets_parts.popt[0]) 
            '''
            # if we wanted to plot the fit
            the_xrecipe = np.linspace(rhosph_log.min(),rhosph_log.max(),num=500)
            the_xten =10**the_xrecipe
            the_yrecipe_ = afunct(the_xrecipe, *rets.popt)
            the_yten = 10**the_yrecipe
            ''' 
            # one panel for each time frame; edit x,y labels and saving names as necessary
            if individually == 'yes':
                fig,ax = plt.subplots(1,1)
                ax.scatter(rho_sph[i], b_sph[i], color='b', label='spheres', alpha=0.5)  
                if parts_or_spheres == 'parts':
                    ax.scatter(rho_parts[i], b_parts[i], color='r', label='particles', alpha=0.5)  
                ax.legend()
                ax.set(xlabel=r'$\rho_{ave}$', ylabel=r'$|B|$', xscale='log', yscale='log', xlim=(1e-2,1e8), ylim=(1e0,1e4))
                outname = 'p66_brho/sphtsung_frame%d_%s'%(i,sims[0])
                plt.savefig(outname)
                print('figure saved!')
                plt.clf()      

        # time frames in time series; edit x,y labels and saving names as necessary
        if individually == 'no':
            fig,ax = plt.subplots(1,1)
            the_x = np.linspace(1,len(rho_sph),len(rho_sph))  
            if kappadyn == 'no':
                kappas_ratio = kappas_sph/kappas_sph_synth
                ax.scatter(the_x, kappas_ratio, color='orange', label='ratio')  
                outname = 'p66_brho/sphwsynth_kappadynratio_scatter_%s'%sims[0] #tsung or tff
            if kappadyn == 'yes':
                #ax.scatter(the_x, kappas_sph, color='b', label='spheres_r1')  
                if whichradius == 'rinf':
                    ax.scatter(the_x, kappas_sph_rinf, color='b', label='spheres_rinf')  
                if parts_or_spheres == 'synthetic':
                    #ax.scatter(the_x, kappas_sph_synth, color='g', label='synthetic_r1')  
                    if whichradius == 'rinf':
                        ax.scatter(the_x, kappas_sph_rinf_synth, color='g', label='synthetic_rinf')  
                if parts_or_spheres == 'parts':
                    ax.plot(the_x, kappas_parts, color='r', label='particles')  
                outname = 'p66_brho/sphsynth_rinf_kappadyn_scatter_%s'%sims[0]  #tsung or tff
            ax.legend()
            xlabels = [r'$t_{tsung,dummy}$', r'$t_{tff,dummy}$'] 
            ylabels = [r'$\kappa/\kappa_synth$', r'$\kappa$' ]  
            ax.set(xlabel=xlabels[1], ylabel=ylabels[1])#, ylim=(0,4)) #,\kappa_{parts}$', ylim=(0,0.95))
            plt.savefig(outname)
            print('figure saved!')
            plt.clf()    


    if figtype == 'kappatff':
        rhotff_sph = []
        btff_sph = []
        rhotff_parts = []
        btff_parts = []
        fig,ax = plt.subplots(1,1)
        for i in range(len(rho_sph)):   
            rhotff_sph = np.append(rhotff_sph, rho_sph[i]) 
            btff_sph = np.append(btff_sph, b_sph[i]) 
            if parts_or_spheres == 'parts':
                rhotff_parts = np.append(rhotff_parts, rho_parts[i]) 
                btff_parts = np.append(btff_parts, b_parts[i])   

        rhotff_sphlog = np.log10(rhotff_sph)
        btff_sphlog = np.log10(btff_sph) 
        rets_sphlog = cfp.fit(afunct, rhotff_sphlog, btff_sphlog)  
        if parts_or_spheres == 'parts':
            rhotff_partslog = np.log10(rhotff_parts)
            btff_partslog = np.log10(btff_parts)
            rets_partslog = cfp.fit(afunct, rhotff_partslog, btff_partslog)  

        tmap = rainbow_map(len(rhotff_sph)) 
        ctr = [tmap(n) for n in range(len(rhotff_sph))]
        color_opt = [ctr,'b','r']       
        if 1:
            ax.scatter(rhotff_sph, btff_sph, color=color_opt[0], alpha=0.4, label='spheres')
        if 0:
            ax.scatter(rhotff_parts, btff_parts, color=color_opt[0], alpha=0.4, label='particles')

        title=[r'$\kappa_{sph_tsung}=%f$'%rets_sphlog.popt[0]]#, r'$\kappa_{sph}=%f$, $\kappa_{parts}=%f$'%(rets_sphlog.popt[0],rets_partslog.popt[0]), 'particles', 'spheres']
        ax.set(xlabel=r'$\rho_{ave}$', ylabel=r'$|B|$', xscale='log', yscale='log', xlim=(10e-3,10e3), ylim=(10e-2,10e3), \
               title=title[0])  
        ax.legend(loc='best')

        outname = 'p66_brho/%s_kappatff_scatter_%s'%(parts_or_spheres,sims[0])
        plt.savefig(outname)
        print('figure saved!')
        plt.clf()    


    # to make histograms of the rinfs of cores per frame
    if figtype == 'rinfradii':
        typeofhistos = 'original'
        for i in range(len(rinf_radii)):   
            if typeofhistos == 'freedman':
                IQR = np.percentile(rinf_radii[i], 75) - np.percentile(rinf_radii[i], 25)
                fd_bin_width = 2 * IQR / (len(rinf_radii[i]))**(1/3)
                fd_bins = int(np.ceil((np.max(rinf_radii[i]) - np.min(rinf_radii[i])) / fd_bin_width))
                sns.histplot(rinf_radii[i], bins=fd_bins, kde=False)
                print('saving frame %d'%i)
                plt.savefig('p66_brho/plotexplore/sns_rinfs_frame_%d'%i)
                plt.clf()    

            if typeofhistos == 'pdfform':
                fig,ax = plt.subplots(1,1)
                the_min = rinf_radii[i].min()
                the_max = rinf_radii[i].max()
                the_bins = np.linspace(the_min, the_max, num=64)
                the_PDF, xbins = np.histogram(rinf_radii[i], bins=the_bins, weights = None, density=True)
                bin_centers = 0.5*(xbins[1:]+xbins[:-1])
                plot_x = bin_centers
                plot_y = the_PDF 
                ax.plot(plot_x,plot_y, c='k')
                print('saving frame %d'%i)
                plt.savefig('p66_brho/plotexplore/rinfs_frame_%d'%i)
                plt.clf()    

            if typeofhistos == 'original':
                fig,ax = plt.subplots(1,1)
                the_min = rinf_radii[i].min()
                the_max = rinf_radii[i].max()
                the_bins = np.linspace(the_min, the_max, num=64)  #there was a comment made at the group meet about this.
                ax.hist(rinf_radii[i], bins=the_bins, density=True, histtype='step', color='k')
                print('saving frame %d'%i)
                plt.savefig('p66_brho/plotexplore/og_rinfs_frame_%d'%i)

    if figtype == 'tsungtff':
        fig,ax = plt.subplots(1,1)
        tff_ratio = tsung_tff/t_ff
        the_min = tff_ratio.min()
        the_max = tff_ratio.max()
        the_bins = np.linspace(the_min, the_max, num=64)  #there was a comment made at the group meet about this.
        ax.hist(tff_ratio, bins=the_bins, density=True, histtype='step', color='k')
        tsungtffratio = tff_ratio.mean()
        ax.set(title='mean tsung/tff ratio: %f'%tsungtffratio, xlabel=r'$t_{\mathrm{sung}}/t_{\mathrm{ff}}$', ylabel='count')
        print('saving histos')
        plt.savefig('p66_brho/plotexplore/tsungratiotff_%s'%sims[0])

    print('closing fig and h5 file!')
    Fptr.close()






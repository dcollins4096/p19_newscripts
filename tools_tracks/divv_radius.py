from starter2 import *
import data_locations as dl
from collections import defaultdict

import davetools
reload(davetools)

plt.close('all')

def dcolormesh(ax,xbins, ybins, hist):

    nx = len(xbins) ; ny=len(ybins)
    TheX = np.r_[(ny)*[xbins]].transpose()
    TheY = np.r_[(nx)*[ybins]]
    cmap = copy.copy(mpl.cm.get_cmap("viridis"))
    cmap.set_under('w')
    minmin = hist[hist>0].min()
    norm = mpl.colors.LogNorm(vmin=minmin,vmax=hist.max())
    ploot=ax.pcolormesh(TheX, TheY, hist, cmap=cmap,norm=norm,shading='nearest')
    return ploot

class div_r_tool():
    def __init__(self,this_looper):
        self.this_looper=this_looper
        self.unique_mass=defaultdict(list)
        self.divv=defaultdict(list)
        self.cores_used=[]
    def run(self,core_list=None, divmin=None, nybins=100):
        dx=1./2048
        nx = 1./dx
        thtr = self.this_looper.tr
        all_cores = np.unique(thtr.core_ids)
        if core_list is None:
            core_list = all_cores

        thtr.sort_time()

        tsorted = thtr.times
        self.core_list=core_list
        for core_id in core_list:
            ms = trackage.mini_scrubber(thtr,core_id)
            if ms.nparticles < 10:
                continue
            print('go ', core_id)
            self.cores_used.append(core_id)
            G=1620./(4*np.pi)
            tff_global = np.sqrt(3*np.pi/(32*G*1))

            self.times = thtr.times/tff_global


            divv = thtr.c([core_id],'velocity_divergence')
            cell_volume = thtr.c([core_id],'cell_volume')
            the_r =  ms.r
            the_r[ the_r < 1/2048] = 1/2048

            maxmax = np.abs(divv).max()
            maxmax = 1000

            fig,ax = plt.subplots(1,1)
            ax.clear()
            rm = rainbow_map(len(self.times))
            for nt, time in enumerate(self.times):
                ax.scatter(the_r[:,nt], divv[:,nt],c=rm(nt, ms.nparticles),s=0.1)
                ax.set_yscale('symlog',linthresh=1000)
                ax.set_xscale('log')
                ax.set_ylim(-10000,1200)
                ax.set_xlim(1/2048,0.5)

                outname = 'plots_to_sort/divv_radius_%s_c%04d_n%04d.png'%(self.this_looper.out_prefix, core_id,nt)
                fig.savefig(outname)
                print(outname)

            
            ax.grid(True)
            fig.savefig('plots_to_sort/divv_radius_%s_c%04d.png'%(self.this_looper.out_prefix, core_id))

            plt.close(fig)


            #self.h = np.histogram()


import three_loopers_1tff as tl
if 'clobber' not in dir():
    clobber=True
if 'div_r_tool1' not in dir() or clobber:
    div_r_tool1=div_r_tool(tl.looper1)
    core_list = [31]
    #core_list = None
    div_r_tool1.run(core_list=core_list, divmin=None,nybins=100)


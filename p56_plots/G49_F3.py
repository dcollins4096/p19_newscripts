from starter2 import *

import convex_hull_tools as CHT
import convex_hull_plot2d as CHP
reload(CHP)
import matplotlib.colors as colors
import figure_sublabel
reload(figure_sublabel)

reload(CHT)
import hair_dryer
reload(hair_dryer)
import close_tool
#import three_loopers_tenfour as TL4
import three_loopers_six as TL6
MOD = TL6
sim_list=['u601','u602','u603']
sim_list=['u603']
#sim_list=['u402']
import supersets
reload(supersets)
if 'ht' not in dir() :
    ht = {}
    for this_simname in sim_list:
        ht[this_simname] = CHT.hull_tool(MOD.loops[this_simname])
        ht[this_simname].make_hulls()
        ht[this_simname].make_overlaps()

#if 'ct' not in dir():
#    ct = {}
#    for this_simname in sim_list:
#        ct[this_simname] = close_tool.close_tool( MOD.loops[this_simname])
#        ct[this_simname].make_distance()

if 'st' not in dir():
    reload(supersets)
    st={}
    for this_simname in sim_list:
        st[this_simname] = supersets.superset( MOD.loops[this_simname], ht[this_simname])
        st[this_simname].find()
import colors
if 1:
    for sim in sim_list:
        stool=st[this_simname]
        htool=ht[this_simname]
        fig,ax=plt.subplots(2,2, figsize=(12,8))
        delta = 0.07
        fig.subplots_adjust(wspace=0,hspace=0, left=delta, right=1-delta)
        #for aaa in ax.flatten():
        #    aaa.set_aspect('equal')

        for nset,this_superset in enumerate(stool.supersets):
            if nset != 1:
                continue
            for los in [0,1,2]:
                core_list=np.sort(list(this_superset))
                color_dict = colors.make_core_cmap(core_list, cmap = 'tab20', seed = -1)
                #color_dict = colors.make_core_cmap(core_list)
                if 1:
                    CHP.plot_2d(htool,core_list=core_list,frames=[0],
                                accumulate=True,label_cores=[-1], 
                                prefix = "S%02d"%nset, color_dict=color_dict,
                                axis_to_plot=[0,1,2],
                                plot_square=False,
                                plotstyle = "Figure8",
                                center_image=False,
                                external_axis=ax.flatten()
                               )
            reload(hair_dryer)
            hd = hair_dryer.hair_time(MOD.loops[sim])
            hd.run(core_list=core_list, color_dict=color_dict,
                   output_prefix="%s_S%02d"%(sim,nset),
                  los_list=[1], external_axis=ax[1][1])
            ax[1][0].set_xlim(nar( ax[0][0].get_xlim()))
            ax[1][1].yaxis.tick_right()
            ax[0][1].yaxis.tick_right()
            ax[1][1].yaxis.set_label_position('right')
            ax[0][1].yaxis.set_label_position('right')
            ax[0][0].xaxis.tick_top()
            ax[0][1].xaxis.tick_top()
            ax[0][0].xaxis.set_label_position('top')
            ax[0][1].xaxis.set_label_position('top')

            ax[1][1].set_ylim( ax[1][0].get_ylim())
            ax[1][1].set_ylabel( ax[1][0].get_ylabel())
            ax[1][1].set_xlabel(r'$t/t_{\rm{ff}}$')

            figure_sublabel.labs(0.68,-0.08  ,r'$3a$')(ax[0][0])
            figure_sublabel.labs(0.758,-0.08 ,r'$3b$')(ax[0][1])
            figure_sublabel.labs(0.68,1.15   ,r'$3c$')(ax[1][0])
            figure_sublabel.labs(0.0,1.15    ,r'$3d$')(ax[1][1])




            fig.savefig('plots_to_sort/overlap_hair_%s_S%02d.pdf'%(sim,nset))


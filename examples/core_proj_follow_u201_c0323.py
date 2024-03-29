from starter2 import *
import xtra_energy
import data_locations as dl
import loop_tools
reload(dl)

reload(looper)
reload(loop_apps)
reload(loop_tools)

#import three_loopers_1tff as tl
#this_looper = tl.looper1
import three_loopers_mountain_top as TLM
this_simname = 'u301'
this_looper = TLM.loops[this_simname]


if 'this_looper' not in dir():
    output_base = "%s_cores"%this_simname
    core_list =  np.array([323])
    frame_list = list(range(10))
    fields = ['x','y','z','density']
    derived = []
    core_323_fname = 'u201_c0323_T1.h5'
    #for making
    this_looper = looper.core_looper(directory= dl.sims[this_simname],
                                     sim_name = this_simname,
                                     out_prefix = this_simname,
                                     target_frame = dl.target_frames[this_simname],
                                     frame_list = frame_list,
                                     core_list =  core_list,
                                     fields_from_grid=fields,
                                     derived = derived,
                                     do_shift=False
                                  )
    this_looper.plot_directory = "/home/dccollins/PigPen"
    this_looper.get_target_indices(h5_name=dl.peak_list[this_simname])
    #                                 bad_particle_list=dl.bad_particles.get(this_simname,None))
    this_looper.get_tracks()
    this_looper.save(core_323_fname)


#loop_apps.core_proj_follow(this_looper,axis_list=[0], field='PotentialField',
#                           zoom=False, grids=False, particles=False, moving_center=True) #, frame_list=[1])
#slab = {'zmin':zmin, 'zmax':zmax}
#import shift_snaps
#shift_snaps.shift_snaps(this_looper)

#core_cmap = colors.make_core_cmap( this_looper.core_list)

loop_apps.core_proj_follow_b(this_looper,axis_list=[0], field='density', core_list=[323], slab=False, zoom=True, only_sphere=True, center_on_sphere=True, annotate=False, frame_list = [0,40,100,125], p_size=7)
                           #zoom=False, grids=False, particles=False, moving_center=True) #, frame_list=[1])
#loop_apps.core_proj_multiple(this_looper,axis_list=[0], field='density', color_dict=color_dict, particles=True, fields=False, core_list )#, frame_list=[31])

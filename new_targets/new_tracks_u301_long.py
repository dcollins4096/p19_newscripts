#
# run get_mountain_tops first
#
from starter2 import *
import xtra_energy
import data_locations as dl
reload(dl)
reload(looper)

if 'this_simname' not in dir():
    this_simname = 'u301'

mountain_top_fname = "datasets_small/%s_mountain_tops_take_9.h5"%this_simname
outname = 'u301_long_pos_only_TEST.h5'

if 1:
    """this set of parameters extracts all primitive quantities"""
    target_frame = dl.target_frames[this_simname]
    frame_list = list(range(0,target_frame+1,1))
    frame_list = range(0,8)
    fields = ['x','y','z','density']#,'velocity_magnitude','magnetic_field_strength', 'velocity_divergence']
    #fields += ['velocity_x','velocity_y','velocity_z']
    #fields += ['magnetic_field_%s'%s for s in 'xyz']
    #fields += ['PotentialField']
    fields += ['particle_pos_x', 'particle_pos_y', 'particle_pos_z', 'particle_index']
    derived=[]

print("G1 make")
disk_or_new = 'new'
if ('new_looper' not in dir() and disk_or_new) or True:
    new_looper = looper.core_looper(directory= dl.sims[this_simname],
                                     sim_name = this_simname,
                                     out_prefix = this_simname,
                                     target_frame = dl.target_frames[this_simname],
                                     frame_list = frame_list,
                                     core_list =  None,
                                     fields_from_grid=fields,
                                     derived = derived,
                                    do_shift=False
                                  )
    new_looper.plot_directory = "./plots_to_sort"
    new_looper.read_targets(mountain_top_fname)


print("G2 verify")
if 0:
    #verify and remove bad particles.
    #Don't need to do this every time, it takes a minute.
    #need to have the bad_particles file 
    for frame in frame_list:
        new_looper.verify_all_particles(frame)
    new_looper.save_bad_particles('datasets_small/%s_bad_particles_TEST.h5'%this_simname)

print("G3 get new")
if 1:
    new_looper.read_bad_particles('datasets_small/%s_bad_particles_TEST.h5'%this_simname)
    new_looper.remove_bad_particles()
    new_looper.get_tracks()


if 1:
    import tracks_read_write
    tracks_read_write.save_loop_trackage_only( new_looper, outname)


if 0:
    print("ZEROS: %d"%( (new_looper.tr.track_dict['density']==0).sum() ))
if 'new_looper' not in dir() and disk_or_new == 'disk':
    print("reload")
    file_list=['TEMP.h5']
    new_looper=looper.core_looper(directory=dl.sims[this_simname])
    new_looper.plot_directory = "/home/dccollins/PigPen"
    for nfile,fname in enumerate(file_list):
        new_looper.load_loop(fname)
        print( "File %d of %d"%(nfile,len(file_list)))
    new_looper.tr.sort_time()
    #import loop_tools
    #loop_tools.re_shift_snaps(new_looper)


if 0:
    import colors
    reload(colors)
    core_cmap = colors.make_core_cmap( new_looper.core_list)
    reload(loop_apps)
    loop_apps.core_proj_multiple(new_looper,#core_list=[258],
                                 #frame_list=[0],
                                 axis_list=[1], color_dict=core_cmap, particles=True,
                                 only_sphere=False,zoom=False,
                                 center_on_sphere=False, annotate=True,tracker_positions=True, shifted_tracker=False)

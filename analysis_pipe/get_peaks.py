from starter2 import *
import loop_tools
#directory = 'u05-r4-l4-128-Beta0.2  u10_r4_l4_128-Beta2  u11_r4_l4_128-Beta20'

def get_peaks(trackname):
    this_track = track_info.tracks[trackname]
    frame = this_track.target_frame
    simdir = this_track.sim_directory
    h5name = this_track.peak_fname
    if os.path.exists(h5name):
        print("File exists, exiting.", h5name)
        return 0


    ds = yt.load("%s/DD%04d/data%04d"%(simdir,frame,frame))

    #get the clumps
    master_clump = loop_tools.get_leaf_clumps(ds,small_test=this_track.clump_parameters.get('small_test'),
                                              c_min=this_track.clump_parameters.get('c_min',None),
                                              step=this_track.clump_parameters.get('step',None))
    #get the peaks from the clump.  Save in h5name
    loop_tools.get_peak_indices(master_clump, ds, h5name) #,h5_name="NEW_PEAKS.h5" )
    print("made peaks in %s"%h5name)





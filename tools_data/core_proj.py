from starter2 import *
import annotate_particles_4
reload(annotate_particles_4)


def core_proj_multiple(looper,  field='density', zoom=1,axis_list=[0,1,2], color_dict={}
                    core_list=None,frame_list=None, 
                       annotate_grids=True, plot_particles=True, annotate=False, 
                      annotate_fields=False, annotate_gravity=False, annotate_velocity=False,annotate_velocity_streamlines=False, annotate_lic=False, 
                       annotate_particles=False, annotate_core_ids=True,
                       code_length=True, 
                      tracker_positions=True, shifted_tracker=True, float_positions=False,
                      marker_size=1, verbose=False,
                       zlim=None,
                       cmap='Greys', weight_field=None):
    """
    Plots an collection of cores in a smooth manner.
    FIRST loop over frames,
        Draws boxes around the particles,  
    THEN smooths the path of the edges.
    THEN make all the plots
    """

    tr = looper.tr
    if core_list is None:
        core_list = np.unique(tr.core_ids)
    if frame_list is None:
        frame_list = looper.frame_list
    if monotonic is True:
        #monotonic  ={'zlim':MonotoneEnforcer()}#, 'left':MonotoneEnforcer2(nvalues=6),'right':MonotoneEnforcer2(nvalues=6)}
        monotonic  ={'zlim':MonotoneColorbar()}#, 'left':MonotoneEnforcer2(nvalues=6),'right':MonotoneEnforcer2(nvalues=6)}
        #monotonic  ={'zlim':MonotoneEnforcer(), 'left':MonotoneEnforcer(),'right':MonotoneEnforcer()}
    tracker_index =  [np.where(looper.tr.frames == frame)[0][0] for frame in frame_list]
    times=nar(looper.tr.times[ tracker_index] )
    all_times=looper.tr.times


    #
    #get all the miniscrubbers at once.
    #We should speed this code up.
    #

    if 
    if verbose:
        print("Mini scrubbers")
    mini_scrubbers = {}
    if annotate_particles:
        for core_id in core_list:
            if mean_velocity:
                do_velocity=True
            else:
                do_velocity=False
            mini_scrubbers[core_id]=  trackage.mini_scrubber(looper.tr,core_id, do_velocity=do_velocity)

    for frame in frame_list:

        # Check to see if the image was made already,
        # and skips it if it has.
        if len(core_list) == 1:
            suffix = "c%04d"%core_list[0]
        else:
            suffix = 'multi'
        ds = looper.load(frame)

        #
        # main plot loop
        #
        for ax in axis_list:
            proj = ds.proj(field,ax,center=center, weight_field=weight_field)
            pw = proj.to_pw(center = center, origin='domain')
            
            pw.zoom(zoom)
            if zlim is not None:
                pw.set_zlim(field,zlim[0],zlim[1])

            pw.set_cmap(field,cmap)

            Hcoord=ds.coordinates.x_axis[ax]
            Vcoord=ds.coordinates.y_axis[ax]
            if annotate_lic:
                pw.annotate_line_integral_convolution('magnetic_field_x','magnetic_field_y', lim=(0.5,0.65))
            if annotate_fields:
                #pw.annotate_magnetic_field(plot_args={'color':'b'})
                pw.annotate_streamlines("magnetic_field_x","magnetic_field_y",plot_args={'color':'b'})
            if annotate_velocity:
                pw.annotate_velocity()
            if annotate_velocity_streamlines:
                GH = [YT_velocity_x, YT_velocity_y, YT_velocity_z][Hcoord]
                GV = [YT_velocity_x, YT_velocity_y, YT_velocity_z][Vcoord]
                pw.annotate_streamlines( GH, GV, plot_args={'color':'k'})
            if annotate_gravity:
                GH = [YT_acceleration_x, YT_acceleration_y, YT_acceleration_z][Hcoord]
                GV = [YT_acceleration_x, YT_acceleration_y, YT_acceleration_z][Vcoord]
                pw.annotate_streamlines( GH, GV, plot_args={'color':'g'})


            if annotate_grids:
                pw.annotate_grids()
            if code_length:
                pw.set_axes_unit('code_length')
            if field in ['MagVelAngle']:
                pw.set_cmap('MagVelAngle','hsv')
                pw.set_zlim('MagVelAngle',0,180)
            looper.pw = pw
            if verbose:
                print('save')
            if verbose:
                print('annotate')
            for core_id in core_list:
                positions = position_dict[core_id]
                color=color_dict.get(core_id, 'k')
                #color = [1.0,0.0,0.0,0.8]
                #alpha=0.1
                alpha=0.4
                marker_size=1
                if annotate_core_ids:
                    #pw.annotate_sphere(snapshot.R_centroid,Rmax, circle_args={'color':color} ) #R_mag.max())
                    centroid=positions.mean(axis=0)
                    pw.annotate_text(centroid,
                                     "%d"%core_id,text_args={'color':color}, 
                                     inset_box_args={'visible':False},
                                     coord_system='data')
                if annotate_particles:
                    pw.annotate_these_particles4(1.0, col=[color]*positions.shape[0], positions=positions, 
                                                 p_size=marker_size, alpha=alpha)





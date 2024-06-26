from starter2 import *
import xtra_energy
import three_loopers_six as TL
import mountain_top
import xtra_energy
import pcolormesh_helper as pch
reload(xtra_energy)

def energy_plots(this_looper, core_list=None, r_inflection=None, r_mass=None):

    sim_name = this_looper.sim_name
    thtr=this_looper.tr
    if core_list is None:
        core_list=np.unique(this_looper.tr.core_ids)

    frame=this_looper.target_frame
    ds = this_looper.load(frame)
    xtra_energy.add_energies(ds)
    xtra_energy.add_gravity(ds)
    reload(mountain_top)
    #radius=1e-2
    for core_id in core_list:
        print("Proj energies %s c%04d"%(sim_name, core_id))
        ms = trackage.mini_scrubber(thtr,core_id,do_velocity=True) 
        ms.get_central_at_once(core_id)
        ms.particle_pos(core_id)
        ms.compute_ge(core_id)
        ms.compute_ke(core_id)

        if 0:
            r_infl = r_inflection[core_id]
            radius=max([1/128, 1.2*r_infl])
        if 0:
            radius = max([2/128, ms.r[:,-1].max()])
        if 0:
            radius = 2/128
        if 1:
            radius = 1/128

        radius = ds.arr(radius,'code_length')

        peak = this_looper.targets[core_id].peak_location

        proj_axis=0
        if proj_axis == 1:
            print("ERROR: plotting messed up for axis 1")
            raise
        left = peak - radius.v
        right = peak + radius.v
        dx = 1/2048
        nzones = np.floor((right-left)/dx).astype('int')

        cg = ds.covering_grid(4,left,nzones, num_ghost_zones=1)
        #cg = ds.covering_grid(0,[0.0]*3,[128]*3, num_ghost_zones=1)
        def from_cg(field):
            ccc = cg[field]#.swapaxes(0,2)
            return ccc.v
        def proj(arr):
            pj = arr.sum(axis=proj_axis)
            avg_pj = pj/arr.shape[proj_axis]
            return avg_pj
        def proj_m(arr,weight):
            pj = (arr*weight).sum(axis=proj_axis)
            wj = (weight).sum(axis=proj_axis)
            return pj/wj


        h_axis = ds.coordinates.x_axis[proj_axis]
        v_axis = ds.coordinates.y_axis[proj_axis]
        x_h = proj( from_cg( ('gas','xyz'[v_axis]) ) )#/nzones[h_axis]
        x_v = proj( from_cg( ('gas','xyz'[h_axis]) ) )#/nzones[v_axis]

        TheX, TheY = np.meshgrid( np.unique( x_h), np.unique( x_v))
        KEP = 'ge_ke'
        field_list=[YT_grav_energy, YT_kinetic_energy, YT_ge_ke]


        center = 0.5*(left+right)
        reg = ds.region(center,left,right)
        #YT_ge_ke = ('gas','ge_ke')
        #def kege(field, data):
        #    out = data.ds.arr(np.zeros_like( data[YT_grav_energy].v), 'dimensionless')
        #    ok = np.abs(data[YT_grav_energy])>0
        #    out[ok]=data[YT_grav_energy][ok].v/data[YT_kinetic_energy][ok].v
        #    return out.v
        #ds.add_field(YT_ge_ke, function=kege, sampling_type='cell')
        fig,ax=plt.subplots(1,1)

        density = from_cg(YT_density)
        #field=YT_grav_energy
        field = YT_ge_ke
        #field = YT_grad_rho
        if type(field) is tuple:
            field_name = field[1]
        else:
            field_name=field
        if field_name not in ['keprime_ge']:
            rho= np.abs(from_cg(field))
            dv = from_cg(YT_cell_volume)
            PPP = proj_m(rho,density)
        else:
            raise
            vx = from_cg(YT_velocity_x) - ms.mean_vx[-1]
            vy = from_cg(YT_velocity_y) - ms.mean_vy[-1]
            vz = from_cg(YT_velocity_z) - ms.mean_vz[-1]
            rho= from_cg(YT_density)
            GE = from_cg(YT_grav_energy)
            p = 0.5*rho*(vx*vx+vy*vy+vz*vz)
            PPP = proj(GE/p)

        vmin,vmax = 1e-3,1e3
        cmap = 'seismic'
        norm = mpl.colors.LogNorm( vmin=vmin, vmax=vmax)

        X_to_plot = TheY.transpose()
        Y_to_plot = TheX.transpose()
        PPP_to_plot = PPP.transpose()
        ploot=ax.pcolormesh( X_to_plot, Y_to_plot, PPP_to_plot, norm=norm, shading='nearest', cmap=cmap)
        fig.colorbar(ploot, ax=ax)

        if 0:
            circle_x = peak[v_axis]
            circle_y = peak[h_axis]
            circle = plt.Circle( (circle_x, circle_y), r_infl, edgecolor='r', facecolor='None')
            axlist[nf].add_patch(circle)
    if 0:
        #fig7,ax7=plt.subplots(1,1)
        xyz = [ms.particle_x[:,-1], ms.particle_y[:,-1], ms.particle_z[:,-1]]
        axlist[0].scatter( xyz[h_axis], xyz[v_axis], c='k')
        axlist[1].scatter( xyz[h_axis], xyz[v_axis], c='k')
        axlist[2].scatter( xyz[h_axis], xyz[v_axis], c='k')
        #xyz = [ms.this_x[:,-1], ms.this_y[:,-1], ms.this_z[:,-1]]
        #xyz = [ms.particle_x[:,-1], ms.particle_y[:,-1], ms.particle_z[:,-1]]
        #axlist[1].scatter( xyz[h_axis], xyz[v_axis])
        #print(ms.nparticles)
        #fig7.savefig('plots_to_sort/wtf.png')

    outname='plots_to_sort/EG_EK_%s_%s_c%04d_n%04d.png'%('xyz'[proj_axis],sim_name, core_id, frame)
    fig.savefig(outname)
    print(outname)
    return cg


if 0:
    fig,ax=plt.subplots(1,3,figsize=(12,4))

    for nsim,sim in enumerate(TL.loops):
        this_looper = TL.loops[sim]
        core_list=np.unique(this_looper.tr.core_ids)
        mountain_top_fname = "datasets_small/u30%d_mountain_tops_take_9.h5"%(nsim+1)
        this_looper.read_targets_only(mountain_top_fname)
        for core_id in core_list:
            ax[nsim].scatter( this_looper.targets[core_id].min_density,  this_looper.targets[core_id].peak_density)
        axbonk(ax[nsim],xscale='log',yscale='log', xlabel='min density',ylabel='peak density')

    fig.savefig('plots_to_sort/test.png')
sim_list=['u601','u602','u603']
#sim_list=['u601']
#sim_list=['u602','u603']


if 1:
    sim_list=sim_list[0:1]
    for nsim,sim in enumerate(sim_list):
        this_looper = TL.loops[sim]
        core_list=np.unique(this_looper.tr.core_ids)
        mountain_top_fname = "datasets_small/u30%d_mountain_tops_take_9.h5"%(nsim+1)
        if this_looper.targets is None:
            this_looper.read_targets_only(mountain_top_fname)
        infl=None
        massedge=None
        #infl=inflection[sim].rinflection
        #massedge=mass_edge[sim].edge
        #core_list=core_list[13:17]
        core_list=[323]
        
        output=energy_plots(this_looper,core_list=core_list,
                          r_inflection=infl)
                          #r_mass=massedge)


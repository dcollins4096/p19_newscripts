import matplotlib
matplotlib.use('Agg')
import math
import matplotlib.pyplot as plt
import numpy as np
nar = np.array
import time
from importlib import reload
import h5py
import copy 
import pdb
from importlib import reload
from collections import defaultdict
import weakref
import colors
if 0:
    import loop_tools
    reload(loop_tools)
"""
trackage.
the tool package to get history for individual particles
"""

def is_sorted(array):
    zero = np.abs(array-np.sort(array)).sum()
    if zero > 0:
        print("False")
        pdb.set_trace()
        return
    print("True")

class track_manager():
    """Manages tracks.
    track_manager[field] returns a list of [particle, time]
    so 
    for i in track_manager[field]:
        plt.plot(track_manager.time,i)
    plots the field vs. time.
    self.ingest(snapshot) populates the master list
    """

    def __init__(self,my_loop=None, h5ptr=None):
        #self.my_loop = None
        #if my_loop is not None:
        #    self.my_loop = weakref.proxy(my_loop)
        self.track_dict={}
        if h5ptr is None:
            self.particle_ids=nar([],dtype='int64')
            self.core_ids = nar([],dtype='int64')
            self.frames = nar([],dtype='int64')
            self.times = nar([],dtype='float64')
            self.V_rad = nar([],dtype='float64')
            self.V_rel = nar([],dtype='float64')
            self.shape=(0,0) #particles, frames
        else:
            #self.shape=(0,0) 
            self.read(fptr=h5ptr)

    def sort_time(self):
        """fix up the time ordering of the members"""
        #pdb.set_trace()
        asort =  np.argsort(self.times)
        if (asort != sorted(asort)).any():
            for k in self.track_dict:
                self.track_dict[k]=self.track_dict[k][:,asort]
            self.times = self.times[asort]
            self.frames = self.frames[asort]

    def write(self,fname="TRACKAGE_SAVE.ht", fptr=None):
        create_file_here = False
        if fptr is None:
            create_file_here=True
            fptr = h5py.File(fname,'w')
        try:
            fptr['particle_ids'] = self.particle_ids
            fptr['core_ids'] = self.core_ids
            fptr['frames'] = self.frames
            fptr['times'] = self.times
            for k in self.track_dict:
                fptr[k] = self.track_dict[k]
        except:
            raise
        finally:
            if create_file_here:
                fptr.close()

    def read(self,fname = None, fptr=None):
        create_file_here = False
        if fptr is None:
            create_file_here=True
            fptr = h5py.File(fname,'r')
        try:
            for key in fptr:
                if str(key) in ['particle_ids', 'core_ids', 'frames', 'times']:
                    self.__dict__[key]=fptr[key][:]
                else:
                    self.track_dict[key] = fptr[key][:] 
        except:
            raise
        finally:
            if create_file_here:
                fptr.close()

    def merge(self,fname):
        fptr = h5py.File(fname,'r')
        temp_dict = {}
        temp_track_dict = {}
        try:
            for key in fptr:
                if str(key) in ['particle_ids', 'core_ids', 'frames', 'times']:
                    temp_dict[key]=fptr[key][:]
                else:
                    temp_track_dict[key] = fptr[key][:] 
        except:
            raise
        finally:
            fptr.close()
        if len(temp_dict['particle_ids']) != len(self.particle_ids):
            print("Error with merging: wrong number of particles")
        elif np.abs(temp_dict['particle_ids'] - self.particle_ids).sum() > 0:
            print("error with merging: wrong particles")
        elif np.abs( temp_dict['core_ids']-self.core_ids).sum() > 0:
            print("error with merging: core ids")
        else:
            frames = np.concatenate([self.frames,temp_dict['frames']])
            args= np.argsort(frames)
            self.frames=frames[args]
            self.times=np.concatenate([self.times,temp_dict['times']])[args]
            for key in self.track_dict:
                oot = np.concatenate( [self.track_dict[key], temp_track_dict[key]],axis=1)[:,args]
                self.track_dict[key]=oot

        return {'td':temp_dict,'ttd':temp_track_dict}

    def ingest(self,snapshot):
        #pdb.set_trace()
        particle_ids = copy.copy(snapshot.ind)

        if hasattr(snapshot, 'core_id'):
            """old style looper"""
            if snapshot.core_id not in self.core_ids:
                #this might not be the best place for the parent step.
                core_ids = np.ones_like(particle_ids) * snapshot.core_id
                if hasattr(core_ids,'v'):
                    core_ids = core_ids.v #need it to not have units.
                self.core_ids = np.append(self.core_ids, core_ids)
                self.particle_ids = np.append(self.particle_ids, particle_ids)

            particle_start = np.where(self.particle_ids==particle_ids[0])[0][0]
            particle_end=particle_start+particle_ids.size
        else:
            """looper2"""

            if len(self.particle_ids) == 0:
                self.particle_ids = copy.copy(particle_ids)
                self.core_ids = copy.copy(snapshot.core_ids)
            if snapshot.core_ids[0] not in self.core_ids:
                self.particle_ids = np.concatenate([self.particle_ids, particle_ids])
                self.core_ids = np.concatenate([self.core_ids,snapshot.core_ids])
            #particle_start = np.where(self.particle_ids==particle_ids[0])[0][0]
            #particle_end=particle_start+particle_ids.size
            particle_start = np.where(self.core_ids == snapshot.core_ids[0])[0][0]
            particle_end = particle_start + particle_ids.size
            
            
            self_cores_set=set(self.core_ids)
            snap_cores_set=set(snapshot.core_ids)
            if not self_cores_set.issuperset(snap_cores_set):
                if not self_cores_set.isdisjoint(snap_cores_set):
                    print(" Don't know what to do with only some of the cores being present.")
                    pdb.set_trace()
                self.core_ids = np.concatenate([self.core_ids, snapshot.core_ids])
                self.particle_ids = np.concatenate([self.particle_ids, particle_ids])
                particle_start = np.where( self.core_ids ==  snapshot.core_ids[0])[0][0]
                particle_end = particle_start + particle_ids.size


        #check that the particles we're inserting
        #are in the same order as the previous ones
        check_particles = np.abs(self.particle_ids[particle_start:particle_end] - particle_ids).max()

        if check_particles:
            print("Fatal Error: particle order.  Email collins.")
            pdb.set_trace()

        #check if we were handed a list or an int
        frames = snapshot.frame
        try:
            test = iter(frames)
        except:
            frames = [frames]
        times = snapshot.time
        try:
            test = iter(snapshot.time)
        except:
            times = [snapshot.time]

        for nf,frame in enumerate(frames):
            if frame not in self.frames:
                self.frames=np.append(self.frames,frame)

                self.times=np.append(self.times,times[nf]) #ds['InitialTime'])

        frame_id = np.where(self.frames == frames[0])[0][0]


        for yt_field_name in snapshot.field_values:
            if type(yt_field_name) is tuple:
                field = yt_field_name[1]
            else:
                field=yt_field_name
            current_shape = self[field].shape
            new_shape = [self.particle_ids.size,
                         self.frames.size]
            temp_frame = np.zeros(new_shape,dtype='float64')
            old_slice = (slice(None,current_shape[0]),
                         slice(None,current_shape[1]))
            temp_frame[old_slice]= self[field]
            new_slice = (slice(particle_start,particle_end),
                         slice(frame_id,frame_id+len(frames)))
            nuggle=np.array(snapshot.field_values[yt_field_name])
            nuggle.shape=(particle_ids.size,len(frames)) #does this break something?
            temp_frame[new_slice]=nuggle
            self[field]=temp_frame
    
    def p(self,particle_list,field):
        output = None
        for particle in particle_list:
            loc = np.where(self.particle_ids == particle)
            parts = self[field][loc,:]
            if output is None:
                output = parts
            else:
                output= np.append(output, parts,axis=0)
        return output

    def c(self,core_list,field):
        if field in ['particle_id']:
            for ncore, core_id in enumerate(core_list):
                these_pids = self.particle_ids[ self.core_ids == core_id]
                if ncore == 0:
                    output = these_pids
                else:
                    output = np.append(output, these_pids)
        else:
            if type(core_list) is int:
                core_list = [core_list]
            nf = self.frames.size
            output = None
            for core in core_list:
                loc = self.core_ids == core
                core_values = self[field][loc,:]
                if output is None:
                    output = core_values
                else:
                    output = np.append(output,core_values,axis=0)
        return output

    def __getitem__(self,item):
        output = None
        if item in ['particles', 'particle_ids']:
            return self.particle_ids
        if item in ['core_ids','cores']:
            return self.core_ids
        if item not in self.track_dict:
            self.track_dict[item]=np.zeros(self.shape,dtype='float64')
        return self.track_dict[item]

    def __setitem__(self,item,value):
        self.track_dict[item]=value
    #can I remove this?
    #def setup_fields(self,field_list=None):
    #    for frame in self.my_loop.field_list:
    #        for core_id in self.my_loop.core_list:
    #            this_snapshot = looper.make_snapshot(frame,core_id)
    #            if field_list is None:
    #                local_field_list = this_snapshot.field_values.keys()
    #            else:
    #                local_field_list = field_list
    #            for field in local_field_list:
    #                self[field].ingest(this_snapshot)


def shift_down(pos):
    #shift based on the time history: if a particle jumps more than half the box,
    #move it.
    # NOT A GOOD IDEA.  Later realized that this is buggy
    sign = -1
    out = copy.copy(pos)
    mean_pos = np.median(out[:,-1])
    distance_from_final =        out- mean_pos
    ft = np.abs(distance_from_final) > 0.5
    out[ft] -=  1*np.sign(distance_from_final[ft])
    return out#,delta

def shift_6(pos):
    #Shift based on large shifts in the trajectory.
    #Assume we want to be on the same side as the end.
    #Find all the times when the particle jumps by more than half the box,
    #and shift it back.
    out = copy.copy(pos)
    delta = out[1:] - out[:-1]
    index =  np.where(np.abs(delta) > 0.5)[0]
    sign = list(np.sign( delta[index]))
    index=list(index)
    if len(index) % 2:
        index = [-1]+index
        sign = [0]+sign
    while len(index):
        Ie = index.pop(-1)+1
        Is = index.pop(-1)+1
        Sgn = sign.pop(-1)
        othersgn = sign.pop(-1)
        out[Is:Ie]  += Sgn
    return out

def shift_4_minus(arr):
    #Loop over particles in the array and call shift_6
    #first version, kill later
    out = np.zeros_like(arr)
    for n,p in enumerate(arr):
        out[n,:]=shift_6(arr[n,:])
    return out  

def shift_4(arr):
    #Loop over particles in the array and call shift_6
    #first shift based on each particles endpoints
    out = np.zeros_like(arr)
    for n,p in enumerate(arr):
        out[n,:]=shift_6(arr[n,:])
    #edit: this code is obsolote if we enforce 
    #      having the last point and enforce the last point is in [0,1]
    #then shift entire tracks that are on the wrong side
    #centroid = out[:,-1].mean()
    #delta = out[:,-1]-centroid
    #if ( delta > 0.5).any():
    #    more_shift = np.where( delta>0.5)[0]
    #    shift = np.sign(delta[more_shift])
    #    shift.shape = shift.size,1
    #    out[more_shift,:] -= shift
    return out  


        
class mini_scrubber():
    def __init__(self,trk,core_id,do_velocity=True, do_magnetic=False, do_means=True, get_help=True, do_central=False,
                do_ge=False,do_ke=False):
        self.trk=trk
        self.scrub(core_id,do_velocity=do_velocity, do_magnetic=do_magnetic, do_means=do_means, get_help=get_help)
        if do_central:
            self.get_central_at_once(core_id)
        if do_ge:
            self.compute_ge(core_id)
        if do_ke:
            self.compute_ke_rel(core_id)
        self.axis=0
                
    def compute_unique_mask(self,core_id, dx,frame):
        """Computes the unique particle mask.
        Code note: we should get this tool to figure out dx by itself."""
        nx = int(1/dx)
        ix =np.floor(self.trk.c([core_id],'x')/dx)[:,frame]
        iy =np.floor(self.trk.c([core_id],'y')/dx)[:,frame]
        iz =np.floor(self.trk.c([core_id],'z')/dx)[:,frame]
        index = ix + nx*(iy * nx*iz)
        ar = np.argsort(index)
        rs = np.argsort(ar)
        isorted=index[ar]
        mask = np.ones_like(ix,dtype='bool')
        mask[1:] = isorted[1:]-isorted[:-1] != 0
        mask2 = mask[ rs]
        return mask2

    def scrub(self,core_id, axis=0, do_velocity=True, do_magnetic=False, do_means=True, get_help=True):
        if core_id not in self.trk.core_ids:
            print("Core %d not found in looper"%core_id)
            print("  (also please write a better error handler)")
            raise
        self.particle_id = self.trk.c([core_id],'particle_id')
        self.raw_x = self.trk.c([core_id],'x')
        self.raw_y = self.trk.c([core_id],'y')
        self.raw_z = self.trk.c([core_id],'z')

        if 1:
            #do the shift
            self.this_x = shift_4(self.raw_x)
            self.this_y = shift_4(self.raw_y)
            self.this_z = shift_4(self.raw_z)
        else:
            #don't actually shift, for testing purposes
            self.this_x = self.raw_x+0
            self.this_y = self.raw_y+0
            self.this_z = self.raw_z+0

        if get_help:
            self.density = self.trk.c([core_id],'density')
            self.cell_volume = self.trk.c([core_id],'cell_volume')
            self.mass = self.density*self.cell_volume
            self.mass_total=self.mass.sum(axis=0)
            self.density_tot = self.density.sum(axis=0)
            self.particle_ids = self.trk.c([core_id],'particle_id')


        if do_means:

            if 0: 
                #Temp debugging code.  Kill later.
                if 'shift_x' not in self.trk.track_dict or True:
                    self.trk.track_dict['shift_x']=np.zeros_like( self.trk.track_dict['density'])
                    self.trk.track_dict['shift_y']=np.zeros_like( self.trk.track_dict['density'])
                    self.trk.track_dict['shift_z']=np.zeros_like( self.trk.track_dict['density'])
                    self.trk.track_dict['test_rho']=np.zeros_like( self.trk.track_dict['density'])
                    self.trk.track_dict['test_z']=np.zeros_like( self.trk.track_dict['density'])
                core_mask = self.trk.core_ids == core_id
                self.shift_x = self.this_x - self.raw_x
                self.shift_y = self.this_y - self.raw_y
                self.shift_z = self.this_z - self.raw_z
                self.trk.track_dict['shift_x'][core_mask,:]= self.shift_x
                self.trk.track_dict['shift_y'][core_mask,:]= self.shift_y
                self.trk.track_dict['shift_z'][core_mask,:]= self.shift_z
                self.trk.track_dict['test_rho'][core_mask,:]= self.raw_z
                self.trk.track_dict['test_z'][core_mask,:]= self.this_z

            #print("kludge: raw mean")
            self.mean_x = np.mean(self.this_x,axis=0)
            self.mean_y = np.mean(self.this_y,axis=0)
            self.mean_z = np.mean(self.this_z,axis=0)
            self.mean_center=nar([self.mean_x,self.mean_y,self.mean_z])
            #self.mean_x = np.sum(self.this_x*self.mass,axis=0)/self.mass_total
            #self.mean_y = np.sum(self.this_y*self.mass,axis=0)/self.mass_total
            #self.mean_z = np.sum(self.this_z*self.mass,axis=0)/self.mass_total
            self.mean_xc = np.sum(self.this_x*self.density,axis=0)/self.density_tot
            self.mean_yc = np.sum(self.this_y*self.density,axis=0)/self.density_tot
            self.mean_zc = np.sum(self.this_z*self.density,axis=0)/self.density_tot
            self.mean_center_density=nar([self.mean_xc,self.mean_yc,self.mean_zc])


            self.nparticles,self.ntimes=self.this_x.shape
            self.meanx2 = np.tile(self.mean_x,(self.raw_x.shape[0],1))
            self.meany2 = np.tile(self.mean_y,(self.raw_x.shape[0],1))
            self.meanz2 = np.tile(self.mean_z,(self.raw_z.shape[0],1))

            self.meanx2c = np.tile(self.mean_xc,(self.raw_x.shape[0],1))
            self.meany2c = np.tile(self.mean_yc,(self.raw_x.shape[0],1))
            self.meanz2c = np.tile(self.mean_zc,(self.raw_z.shape[0],1))


            self.rx_rel=self.this_x-self.meanx2
            self.ry_rel=self.this_y-self.meany2
            self.rz_rel=self.this_z-self.meanz2

            self.rx_relc=self.this_x-self.meanx2c
            self.ry_relc=self.this_y-self.meany2c
            self.rz_relc=self.this_z-self.meanz2c

            self.mean_xm = np.sum(self.this_x*self.density*self.cell_volume,axis=0)/self.mass_total
            self.mean_ym = np.sum(self.this_y*self.density*self.cell_volume,axis=0)/self.mass_total
            self.mean_zm = np.sum(self.this_z*self.density*self.cell_volume,axis=0)/self.mass_total
            self.mean_xmt = np.tile(self.mean_xm,(self.raw_x.shape[0],1))
            self.mean_ymt = np.tile(self.mean_ym,(self.raw_x.shape[0],1))
            self.mean_zmt = np.tile(self.mean_zm,(self.raw_z.shape[0],1))
            self.rx_relm=self.this_x-self.mean_xmt
            self.ry_relm=self.this_y-self.mean_ymt
            self.rz_relm=self.this_z-self.mean_zmt
            self.r2m = self.rx_relm**2+self.ry_relm**2+self.rz_relm**2
            self.rm = np.sqrt(self.r2m)

            self.r2 = self.rx_rel**2+self.ry_rel**2+self.rz_rel**2
            self.r2c = self.rx_relc**2+self.ry_relc**2+self.rz_relc**2

            self.r=np.sqrt(self.r2)
            self.rc=np.sqrt(self.r2c)
            self.rmax = np.max(self.r,axis=0)

            self.rms = np.sqrt( np.mean(self.r2,axis=0))

        if do_magnetic:
            self.bx = self.trk.c([core_id],'magnetic_field_x')
            self.by = self.trk.c([core_id],'magnetic_field_y')
            self.bz = self.trk.c([core_id],'magnetic_field_z')
            self.b2 = self.bx*self.bx+self.by*self.by+self.bz*self.bz

        if do_velocity:
            self.raw_vx = self.trk.c([core_id],'velocity_x')
            self.raw_vy = self.trk.c([core_id],'velocity_y')
            self.raw_vz = self.trk.c([core_id],'velocity_z')

            self.sqr_vx = np.sum(self.raw_vx**2,axis=0)
            self.sqr_vy = np.sum(self.raw_vy**2,axis=0)
            self.sqr_vz = np.sum(self.raw_vz**2,axis=0)

            self.raw_v2 = self.raw_vx**2+self.raw_vy**2+self.raw_vz**2
            #print("KLUDGE: using raw mean for velocity")
            self.mean_vx = np.mean(self.raw_vx,axis=0)
            self.mean_vy = np.mean(self.raw_vy,axis=0)
            self.mean_vz = np.mean(self.raw_vz,axis=0)

            self.mass_mean_vx = np.sum(self.raw_vx*self.mass,axis=0)/self.mass_total
            self.mass_mean_vy = np.sum(self.raw_vy*self.mass,axis=0)/self.mass_total
            self.mass_mean_vz = np.sum(self.raw_vz*self.mass,axis=0)/self.mass_total

            self.rel_vx = self.raw_vx-self.mean_vx
            self.rel_vy = self.raw_vy-self.mean_vy
            self.rel_vz = self.raw_vz-self.mean_vz
            self.rel_vmag = (self.rel_vx**2+self.rel_vy**2+self.rel_vz**2)**(0.5)

            self.rx_hat = np.zeros_like(self.r)
            self.ry_hat = np.zeros_like(self.r)
            self.rz_hat = np.zeros_like(self.r)
            ok = self.r > 0
            self.rx_hat[ok]= self.rx_rel[ok]/self.r[ok]
            self.ry_hat[ok]= self.ry_rel[ok]/self.r[ok]
            self.rz_hat[ok]= self.rz_rel[ok]/self.r[ok]

            self.norm_r = (self.rx_hat**2+self.ry_hat**2+self.rz_hat**2)**(0.5)

            self.vr_raw = self.rx_hat*self.raw_vx+\
                      self.ry_hat*self.raw_vy+\
                      self.rz_hat*self.raw_vz
            self.vr_rel = self.rx_hat*self.rel_vx+\
                      self.ry_hat*self.rel_vy+\
                      self.rz_hat*self.rel_vz
            self.vr_x = self.vr_rel*self.rx_hat
            self.vr_y = self.vr_rel*self.ry_hat
            self.vr_z = self.vr_rel*self.rz_hat
            self.vt2_raw = (self.raw_vx-self.vr_raw*self.rx_hat)**2+\
                       (self.raw_vy-self.vr_raw*self.ry_hat)**2+\
                       (self.raw_vz-self.vr_raw*self.rz_hat)**2
            self.vt2_rel = (self.rel_vx-self.vr_rel*self.rx_hat)**2+\
                       (self.rel_vy-self.vr_rel*self.ry_hat)**2+\
                       (self.rel_vz-self.vr_rel*self.rz_hat)**2
            self.vt_x = self.rel_vx-self.vr_rel*self.rx_hat
            self.vt_y = self.rel_vy-self.vr_rel*self.ry_hat
            self.vt_z = self.rel_vz-self.vr_rel*self.rz_hat


        self.axis = axis
        if self.axis == 0:
            self.this_h = self.this_y
            self.h_label='y'
            self.this_v = self.this_z
            self.v_label='z'
        elif self.axis == 1:
            self.this_h = self.this_z
            self.h_label='z'
            self.this_v = self.this_x
            self.v_label='x'
        if self.axis == 2:
            self.this_h = self.this_x
            self.h_label='x'
            self.this_v = self.this_y
            self.v_label='y'

    def get_central_velocity(self,core_id,nt):

        this_r = self.r[:,nt]
        vx = self.trk.c([core_id],'velocity_x')[:,nt]
        vy = self.trk.c([core_id],'velocity_y')[:,nt]
        vz = self.trk.c([core_id],'velocity_z')[:,nt]

        asort = np.argsort(this_r)
        ind_min = asort[0]

        vxc = vx[ind_min]
        vyc = vy[ind_min]
        vzc = vz[ind_min]

        self.cen_vx = self.raw_vx-vxc
        self.cen_vy = self.raw_vy-vyc
        self.cen_vz = self.raw_vz-vzc
        self.cen_vmag = (self.cen_vx**2+self.cen_vy**2+self.cen_vz**2)**(0.5)

        self.rc_vmag = self.rx_hat*self.cen_vx+\
                       self.ry_hat*self.cen_vy+\
                       self.rz_hat*self.cen_vz

    def get_central_velocity2(self,core_id,nt):

        #relative to the density weighted center.
        this_r = self.rc[:,nt]
        vx = self.trk.c([core_id],'velocity_x')[:,nt]
        vy = self.trk.c([core_id],'velocity_y')[:,nt]
        vz = self.trk.c([core_id],'velocity_z')[:,nt]

        asort = np.argsort(this_r)
        ind_min = asort[0]

        vxc = vx[ind_min]
        vyc = vy[ind_min]
        vzc = vz[ind_min]

        self.d_cen_vx = self.raw_vx-vxc
        self.d_cen_vy = self.raw_vy-vyc
        self.d_cen_vz = self.raw_vz-vzc
        self.d_cen_vmag = (self.d_cen_vx**2+self.d_cen_vy**2+self.d_cen_vz**2)**(0.5)

        self.rcd_vmag = self.rx_hat*self.d_cen_vx+\
                       self.ry_hat*self.d_cen_vy+\
                       self.rz_hat*self.d_cen_vz
    def get_central_at_once(self,core_id):

        this_r = self.r
        vx = self.trk.c([core_id],'velocity_x')
        vy = self.trk.c([core_id],'velocity_y')
        vz = self.trk.c([core_id],'velocity_z')

        argmin = np.argmin(this_r,axis=0)
        ind = np.arange(this_r.shape[1])
        take = np.ravel_multi_index( nar([argmin,ind]), this_r.shape)

        self.vxc = vx.flatten()[take]
        self.vyc = vy.flatten()[take]
        self.vzc = vz.flatten()[take]
        self.vcentral=np.stack([self.vxc,self.vyc,self.vzc])

        self.cen_vx = self.raw_vx-self.vxc
        self.cen_vy = self.raw_vy-self.vyc
        self.cen_vz = self.raw_vz-self.vzc
        self.cen_vmag = (self.cen_vx**2+self.cen_vy**2+self.cen_vz**2)**(0.5)

        self.rc_vmag = self.rx_hat*self.cen_vx+\
                       self.ry_hat*self.cen_vy+\
                       self.rz_hat*self.cen_vz
    def Moment(self):
        self.moment_of_inertia_z = (self.mass*(self.rx_rel**2+self.ry_rel**2)).sum(axis=0)
        self.moment_of_inertia_x = (self.mass*(self.rz_rel**2+self.ry_rel**2)).sum(axis=0)
        self.moment_of_inertia_y = (self.mass*(self.rz_rel**2+self.rx_rel**2)).sum(axis=0)
        self.moment_of_inertia_xy = self.moment_of_inertia_yx = -(self.mass*(self.rx_rel*self.ry_rel)).sum(axis = 0)
        self.moment_of_inertia_yz = self.moment_of_inertia_zy = -(self.mass*(self.ry_rel*self.rz_rel)).sum(axis = 0)
        self.moment_of_inertia_xz = self.moment_of_inertia_zx = - (self.mass*(self.rx_rel*self.rz_rel)).sum(axis = 0) 
        self.moment_of_inertia_zii = (self.mass*(self.rx_rel**2+self.ry_rel**2)).sum()
        self.moment_of_inertia_xii = (self.mass*(self.rz_rel**2+self.ry_rel**2)).sum()
        self.moment_of_inertia_yii = (self.mass*(self.rz_rel**2+self.rx_rel**2)).sum()
        self.moment_of_inertia_xyii = self.moment_of_inertia_yxii = -(self.mass*(self.rx_rel*self.ry_rel)).sum()
        self.moment_of_inertia_yzii = self.moment_of_inertia_zyii = -(self.mass*(self.ry_rel*self.rz_rel)).sum()
        self.moment_of_inertia_xzii = self.moment_of_inertia_zxii = - (self.mass*(self.rx_rel*self.rz_rel)).sum() 

        self.I_ii = self.mass*(self.r**2)

    def momenta(self):
        self.angular_v_x = ((self.ry_rel*self.rel_vz-self.rz_rel*self.rel_vy)/self.r**2)
        self.angular_v_y = ((self.rz_rel*self.rel_vx - self.rx_rel*self.rel_vz)/self.r**2)
        self.angular_v_z = ((self.rx_rel*self.rel_vy - self.ry_rel*self.rel_vx)/self.r**2)
        self.linear_momentum_rel_x = self.mass*(self.rel_vx)
        self.linear_momentum_rel_y = self.mass*(self.rel_vy)
        self.linear_momentum_rel_z = self.mass*(self.rel_vz)
        self.angular_momentum_rel_x = self.ry_rel*self.linear_momentum_rel_z-self.rz_rel*self.linear_momentum_rel_y
        self.angular_momentum_rel_y = self.rz_rel*self.linear_momentum_rel_x-self.rx_rel*self.linear_momentum_rel_z
        self.angular_momentum_rel_z = self.rx_rel*self.linear_momentum_rel_y-self.ry_rel*self.linear_momentum_rel_x
        self.angular_momentum_mag = (self.angular_momentum_rel_x**2+self.angular_momentum_rel_y**2+self.angular_momentum_rel_z**2)**0.5
        self.j_hat_x = self.angular_momentum_rel_x/self.angular_momentum_mag
        self.j_hat_y = self.angular_momentum_rel_y/self.angular_momentum_mag
        self.j_hat_z = self.angular_momentum_rel_z/self.angular_momentum_mag
        self.j_hat_r = self.angular_momentum_mag
        self.j_hat_theta = np.arctan2(self.j_hat_y,self.j_hat_x)
        self.j_hat_phi = np.arccos(self.j_hat_z)
        self.r_dot_angular_moment = self.rx_rel*self.angular_momentum_rel_x + self.ry_rel*self.angular_momentum_rel_y + self.rz_rel*self.angular_momentum_rel_z

    def particle_pos(self,core_id):
        shift_x = self.this_x - self.raw_x
        shift_y = self.this_y - self.raw_y
        shift_z = self.this_z - self.raw_z
        name_to_use_x= 'particle_pos_x'
        name_to_use_y= 'particle_pos_y'
        name_to_use_z= 'particle_pos_z'
        if name_to_use_x not in self.trk.track_dict:
            name_to_use_x= 'particle_position_x'
            name_to_use_y= 'particle_position_y'
            name_to_use_z= 'particle_position_z'

        self.particle_x = self.trk.c([core_id],name_to_use_x) + shift_x
        self.particle_y = self.trk.c([core_id],name_to_use_y) + shift_y
        self.particle_z = self.trk.c([core_id],name_to_use_z) + shift_z

    def make_floats(self, core_id):
        name_to_use_x= 'particle_pos_x'
        name_to_use_y= 'particle_pos_y'
        name_to_use_z= 'particle_pos_z'
        if name_to_use_x not in self.trk.track_dict:
            name_to_use_x= 'particle_position_x'
            name_to_use_y= 'particle_position_y'
            name_to_use_z= 'particle_position_z'
        self.float_x = self.trk.c([core_id],name_to_use_x)
        self.float_y = self.trk.c([core_id],name_to_use_y)
        self.float_z = self.trk.c([core_id],name_to_use_z)
        shift_x = self.this_x - self.raw_x
        shift_y = self.this_y - self.raw_y
        shift_z = self.this_z - self.raw_z
        self.float_x += shift_x
        self.float_y += shift_y
        self.float_z += shift_z
    def compute_ge(self,core_id):
        self.gx = self.trk.c([core_id],'grav_x')
        self.gy = self.trk.c([core_id],'grav_y')
        self.gz = self.trk.c([core_id],'grav_z')


        self.ge = -1/(np.pi*8*colors.G)*(self.gx**2+self.gy**2+self.gz**2)
    def compute_ke(self,core_id):
        self.ke = 0.5*self.density*(self.raw_vx**2+self.raw_vy**2+self.raw_vz**2)
    def compute_ke_rel(self,core_id):
        self.ke_rel = 0.5*self.density*(self.rel_vx**2+self.rel_vy**2+self.rel_vz**2)
    """
    def compute_Rb(self):
        core_id=self.core_id
            bx   =extract(thtr.c([core_id],'magnetic_field_x'))
            by   =extract(thtr.c([core_id],'magnetic_field_y'))
            bz   =extract(thtr.c([core_id],'magnetic_field_z'))
            vx   =extract(thtr.c([core_id],'velocity_x'))
            vy   =extract(thtr.c([core_id],'velocity_y'))
            vz   =extract(thtr.c([core_id],'velocity_z'))
            dxvx =extract(thtr.c([core_id],'dxvx'))
            dxvy =extract(thtr.c([core_id],'dxvy'))
            dxvz =extract(thtr.c([core_id],'dxvz'))
            dyvx =extract(thtr.c([core_id],'dyvx'))
            dyvy =extract(thtr.c([core_id],'dyvy'))
            dyvz =extract(thtr.c([core_id],'dyvz'))
            dzvx =extract(thtr.c([core_id],'dzvx'))
            dzvy =extract(thtr.c([core_id],'dzvy'))
            dzvz =extract(thtr.c([core_id],'dzvz'))

            if OOM:
                self.dxvx=dxvx
                self.dxvy=dxvy
                self.dxvz=dxvz
                self.dyvx=dyvx
                self.dyvy=dyvy
                self.dyvz=dyvz
                self.dzvx=dzvx
                self.dzvy=dzvy
                self.dzvz=dzvz
                self.dixj=[[self.dxvx,self.dxvy,self.dxvz],
                           [self.dyvx,self.dyvy,self.dyvz],
                           [self.dzvx,self.dzvy,self.dzvz]]
            Sx = bx*dxvx+by*dyvx+bz*dzvx
            Sy = bx*dxvy+by*dyvy+bz*dzvy
            Sz = bx*dxvz+by*dyvz+bz*dzvz
            Stretch= bx*Sx+by*Sy+bz*Sz

            R = Stretch/(B2*divv)
            """




from starter2 import *



class camera_1():
    def __init__(self, looper, method):
        self.looper=looper
        self.all_center={}
        self.all_left={}
        self.all_right={}
        self.all_positions={}
        self.method=method
        self.times=[]
    def run(self, core_list, frame_list, mini_scrubbers):
        if self.method == 'domain':
            self.run_domain(core_list, frame_list,mini_scrubbers)
        if self.method.startswith('tight'):
            self.run_tight(core_list,frame_list, mini_scrubbers )
        if self.method == 'smooth_zoom':
            self.run_smooth_zoom(core_list,frame_list,mini_scrubbers)
        else:
            print("ill defined camera")
            raise
            
    def run_smooth_zoom(self,core_list,frame_list, mini_scrubbers):
        looper=self.looper
        ds = looper.load( frame_list[0])
        particle_left=[]
        particle_right=[]
        for frame in frame_list:
            frame_ind = np.where(looper.tr.frames == frame)[0]
            self.times.append(looper.tr.times[frame_ind])

            #get extents and bounding region
            #it's backwards because we're looking for extrema

            position_dict={}
            left =  ds.domain_right_edge.v
            right = ds.domain_left_edge.v
            for core_id in core_list:
                ms = mini_scrubbers[core_id]
                
                ms.make_floats(core_id)
                this_x = ms.float_x[:,frame_ind]
                this_y = ms.float_y[:,frame_ind]
                this_z = ms.float_z[:,frame_ind]

                positions = np.column_stack([this_x,this_y,this_z])
                position_dict[core_id] = positions

                this_left =  positions.min(axis=0)
                this_right = positions.max(axis=0)
                left = np.row_stack([this_left,left]).min(axis=0)
                right = np.row_stack([this_right,right]).max(axis=0)

            self.all_positions[frame]=position_dict

            center = 0.5*(left+right)
            particle_left.append( left)
            particle_right.append(right)
        particle_left = np.stack(particle_left)
        particle_right = np.stack(particle_right)
        particle_center=0.5*(particle_left+particle_right)
        first_center = particle_center[0,:]
        last_center = particle_center[-1,:]
        first_time = self.times[0]
        last_time = self.times[-1]
        first_left = first_center - 0.25
        first_right= first_center + 0.25
        last_left = last_center - 4./128
        last_right= last_center + 4./128
        dL = last_left - first_left
        dR = last_right - first_right
        dt = last_time - first_time


        if 1:
            for nf,frame in enumerate(frame_list):
                frame_ind = np.where(looper.tr.frames == frame)[0]
                time=looper.tr.times[frame_ind]-first_time
                self.all_left[frame] = first_left + dL/dt * time
                self.all_right[frame] = first_right + dR/dt * time
                self.all_center[frame]=0.5*(self.all_left[frame]+self.all_right[frame])


        if 0:
            for nf,frame in enumerate(frame_list):
                frame_ind = np.where(looper.tr.frames == frame)[0]
                time=looper.tr.times[frame_ind]
                Rt = 0.25 + (4./128-0.25)*time/colors.tff
                self.all_left[frame] = particle_center[nf] - Rt
                self.all_right[frame] = particle_center[nf] + Rt
                self.all_center[frame] = particle_center[nf]


    def run_domain(self,core_list,frame_list, mini_scrubbers):
        looper=self.looper
        ds = looper.load( frame_list[0])
        for frame in frame_list:
            frame_ind = np.where(looper.tr.frames == frame)[0]

            #get extents and bounding region
            #it's backwards because we're looking for extrema

            position_dict={}
            for core_id in core_list:
                ms = mini_scrubbers[core_id]
                
                ms.make_floats(core_id)
                this_x = ms.float_x[:,frame_ind]
                this_y = ms.float_y[:,frame_ind]
                this_z = ms.float_z[:,frame_ind]

                positions = np.column_stack([this_x,this_y,this_z])
                position_dict[core_id] = positions

            self.all_positions[frame]=position_dict

            left =  ds.domain_left_edge.v
            right = ds.domain_right_edge.v
            center = 0.5*(left+right)
            self.all_left[frame] = left
            self.all_right[frame]= right
            self.all_center[frame]=center

    def run_tight(self,core_list,frame_list, mini_scrubbers, derived=None):
        looper=self.looper
        for frame in frame_list:
            ds = looper.load(frame)
            frame_ind = np.where(looper.tr.frames == frame)[0]
            if derived is not None:
                for ddd in derived:
                    ddd(ds)

            #get extents and bounding region
            #it's backwards because we're looking for extrema
            left =  ds.domain_right_edge.v
            right = ds.domain_left_edge.v

            #
            # Find the extents of all cores.
            # Fill position dict
            #
            mean_pos=None
            position_dict={}
            for core_id in core_list:
                ms = mini_scrubbers[core_id]
                
                if self.method == 'tight_this':
                    this_x=ms.this_x[:,frame_ind]
                    this_y=ms.this_y[:,frame_ind]
                    this_z=ms.this_z[:,frame_ind]
                elif self.method == 'tight_raw':
                    this_x=ms.raw_x[:,frame_ind]
                    this_y=ms.raw_y[:,frame_ind]
                    this_z=ms.raw_z[:,frame_ind]
                elif self.method == 'tight_float':
                    ms.make_floats(core_id)
                    this_x = ms.float_x[:,frame_ind]
                    this_y = ms.float_y[:,frame_ind]
                    this_z = ms.float_z[:,frame_ind]

                positions = np.column_stack([this_x,this_y,this_z])
                position_dict[core_id] = positions
                if mean_pos is None:
                    mean_pos=copy.copy(positions)
                else:
                    mean_pos = np.row_stack([mean_pos,positions])

                this_left =  positions.min(axis=0)
                this_right = positions.max(axis=0)
                left = np.row_stack([this_left,left]).min(axis=0)
                right = np.row_stack([this_right,right]).max(axis=0)
            self.all_positions[frame]=position_dict

            center = 0.5*(left+right)
            left = np.row_stack([left,center - 1/128]).min(axis=0)
            right = np.row_stack([right,center + 1/128]).max(axis=0)

            #Make it a square
            if 1:
                centroid = mean_pos.mean(axis=0)
                centroid.shape = 1,centroid.size
                distance=(mean_pos-centroid)
                extent = np.sqrt( ((distance)**2).sum(axis=1).max() )
                size = max([extent,1/128])
                centroid.shape = centroid.size
                square_left  = centroid - size
                square_right = centroid + size

                #left  = np.row_stack([square_left,left]).min(axis=0)
                #right = np.row_stack([square_right,right]).max(axis=0)
                left = square_left
                right = square_right
                if ( left >= right).any():
                    print('wtf')
                    pdb.set_trace()



            self.all_left[frame] = ds.arr(left,'code_length')
            self.all_right[frame]= ds.arr(right,'code_length')
            self.all_center[frame]=ds.arr(center,'code_length')
            if (self.all_left[frame] >= self.all_right[frame]).any():
                print("Should not have left bigger than right")
                pdb.set_trace()


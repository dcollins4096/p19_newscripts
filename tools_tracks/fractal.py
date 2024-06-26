


from starter2 import *
import matplotlib.image as mpimg
from scipy.spatial import ConvexHull
import matplotlib.patches as patches
import data_locations as dl
reload(dl)
plt.close('all')

def make_fractals(looper, framelist, core_list=None):
    
    name =  looper.sim_name
    fractals={}
    for nframe, frame in enumerate(framelist):
        fractals[frame] = fractal_tool( looper )
        fractals[frame].run(frame=frame, do_plots=False, core_list=core_list)
    return fractals

#from https://gist.github.com/viveksck/1110dfca01e4ec2c608515f0d5a5b1d1
def fractal_dimension(Z, threshold=0.9, do_plots=False, plotname='plot'):

    # Only for 3d image
    assert(len(Z.shape) == 3)

    # From https://github.com/rougier/numpy-100 (#87)
    def boxcount(Z, k):
        S=np.add.reduceat(Z, np.arange(0, Z.shape[0], k), axis=0)
        S=np.add.reduceat(S, np.arange(0, Z.shape[1], k), axis=1)
        S=np.add.reduceat(S, np.arange(0, Z.shape[2], k), axis=2)
        #S = np.add.reduceat(
        #    np.add.reduceat(
        #    np.add.reduceat(Z, np.arange(0, Z.shape[0], k), axis=0),
        #                       np.arange(0, Z.shape[1], k), axis=1),
        #                       np.arange(0, Z.shape[1], k), axis=2)

        ## the original did this.  We want full boxes as well.
        ## We count non-empty (0) and non-full boxes (k*k)
        #return len(np.where((S > 0) & (S < k*k))[0])
        return S



    # Transform Z into a binary array
    Z = (Z > threshold)


    # Minimal dimension of image
    p = min(Z.shape)

    # Greatest power of 2 less than or equal to p
    n = 2**np.floor(np.log(p)/np.log(2))
    #pdb.set_trace()

    # Extract the exponent
    n = int(np.log(n)/np.log(2))

    # Build successive box sizes (from 2**n down to 2**1)
    #sizes = 2**np.arange(n, 1, -1)
    sizes = 2**np.arange(n, -1, -1)

    if 0:
        #maybe it makes more sense to keep everything?
        #don't repeat sizes that are too big.
        image_length= np.max( nar(np.where(Z)))
        image_length = 2**np.round(np.log2(image_length))
        sizes = sizes[ sizes <= image_length]

    # Actual box counting with decreasing size
    counts = []
    boxes = []
    for ns,size in enumerate(sizes):
        boxes.append( boxcount(Z,size))
        counts.append( (boxes[-1]>0).sum() )

    if do_plots:
        fig, axes = plt.subplots(int(len(boxes)/2),2, figsize=(8,12))#, figsize=(4*len(sizes),4))
        ax=axes.flatten()
        fig.subplots_adjust(wspace=0, hspace=0)
        for ns,box in enumerate(boxes[1:]):
            img = box.sum(axis=0)
            img = img[:int(img.shape[0]//2),:int(img.shape[1]//2)]
            norm=mpl.colors.Normalize(vmin=1,vmax=img.max())
            cmap=copy.copy(mpl.cm.get_cmap("viridis"))
            cmap.set_under('w')
            ax[ns].imshow(img.transpose(),norm=norm,cmap=cmap, interpolation='nearest',origin='lower')
        for aaa in ax:
            aaa.set_yticks([])
            aaa.set_xticks([])
        fig.savefig(plotname)
        plt.close(fig)

    # Fit the successive log(sizes) with log (counts)

    sizes=nar(sizes)
    counts=nar(counts)
    #ok = counts > 1
    ok = slice(None)
    if True: #ok.sum():
        coeffs = np.polyfit(np.log(sizes[ok]), np.log(counts[ok]), 1)
    else:
        coeffs = [0,0]
    #plt.clf()
    #plt.plot(np.log(sizes), np.log(counts),c='k',marker='x')
    #plt.plot(np.log(sizes), coeffs[0]*np.log(sizes)+coeffs[1])
    #nfractal=len(glob.glob("plots_to_sort/nfractal*")) #plt.savefig("plots_to_sort/nfractal_%s.png"%nfractal)
    return -coeffs[0], {'counts':counts,'sizes':sizes, 'coeffs':coeffs}

dx=1./128
class fractal_tool():
    def __init__(self,this_looper):
        self.this_looper=this_looper
        self.cores_used=[]
        self.fractal_dim=[]

    def run(self,core_list=None,frame=0, do_plots=False, plot_prefix='FRACTAL'):
        thtr = self.this_looper.tr
        all_cores = np.unique(thtr.core_ids)
        rm = rainbow_map(len(all_cores))
        if core_list is None:
            core_list = all_cores
        all_frames=thtr.frames
        nf = np.where(all_frames==frame)[0][0]
        print(nf)
        for core_id in core_list:
            ms = trackage.mini_scrubber(thtr,core_id, do_velocity=True, do_means=True)
            if ms.this_x.shape[0] <= 4:
                continue
            #print("fractal dim on ", core_id)
            self.cores_used.append(core_id)

            #x =np.floor(thtr.c([core_id],'x')/dx)[:,nf]#or whatever the number of zones is
            #y =np.floor(thtr.c([core_id],'y')/dx)[:,nf]
            #z =np.floor(thtr.c([core_id],'z')/dx)[:,nf]
            x = (np.floor(ms.this_x/dx)[:,nf] ).astype('int')
            y = (np.floor(ms.this_y/dx)[:,nf] ).astype('int')
            z = (np.floor(ms.this_z/dx)[:,nf] ).astype('int')
            #x,y,z = np.mgrid[0:128:1, 0:128:1, 37:38:1].astype('int')
            density = thtr.c([core_id],'density')[:,nf]
            cell_volume = thtr.c([core_id],'cell_volume')[:,nf]

            x = x-x.min()
            y = y-y.min()
            z = z-z.min()

            image_size=[x.max()+1,y.max()+1,z.max()+1]
            image = np.zeros([128,128,128])
            mask = (x < 128)*(y<128)*(z<128)
            x = x[ mask]
            y = y[ mask]
            z = z[ mask]
            image[(x,y,z)]=1
                        
                
            #print("Total points", x.size, image.sum())
            plotname = 'plots_to_sort/fractal_covering_%s_c%04d'%(plot_prefix,core_id)
            #print(plotname)
            stuff = fractal_dimension(image,do_plots=do_plots, plotname=plotname)
            dimension = stuff[0]
            self.fractal_dim.append(dimension)
            frame = self.this_looper.tr.frames[nf]
            print("fractal dim c%04d n%04d %0.2f"%(core_id,frame,self.fractal_dim[-1]))
            if 0:
                #image.
                plt.clf()
                cmap = copy.copy(mpl.cm.get_cmap("viridis"))
                cmap.set_under('w')
                image_2d = image.sum(axis=0)
                norm = mpl.colors.Normalize(vmin=1,vmax=image_2d.max())
                plt.imshow( image_2d,cmap=cmap,norm=norm)

                plt.savefig('plots_to_sort/%s_fractal_c%04d.png'%(self.this_looper.out_prefix,core_id))

            if do_plots:
                counts = stuff[1]['counts']
                sizes  = stuff[1]['sizes']
                coeffs  = stuff[1]['coeffs']
                plt.clf()
                plt.plot(np.log10(sizes), np.log(counts),c='k',marker='x')
                plt.plot(np.log10(sizes), coeffs[0]*np.log(sizes)+coeffs[1])
                plt.xlabel('Size'); plt.ylabel('Count')
                plt.title('c%04d D=%0.2f'%(core_id, dimension))
                plt.savefig('plots_to_sort/%s_counts_c%04d.png'%(self.this_looper.out_prefix,core_id))


#import three_loopers_1tff as tl
#import three_loopers_mountain_top as TLM
#import three_loopers_tenfour as TL4
import track_loader as TL
sims = ['u601','u602','u603']
sims = ['u602']
#track_loader.load(sims)


if 'toolshed' not in dir():
    toolshed = {}



if 0:
    framelist=[0]#,1,2,3,4,5,6,7,8]
    looper_list = TL4.loops
    
    for simname in looper_list:
        looper=looper_list[simname]
        name =  looper.sim_name

        if name not in toolshed:
            toolshed[ name ] = {}
            toolshed[ name ]['looper'] = looper
        for nframe, frame in enumerate(framelist):
            if nframe not in framelist:
                continue
            if frame in toolshed[ name]:
                continue
            toolshed[name][frame] = fractal_tool( looper )
            toolshed[name][frame].run(nf=nframe, do_plots=False)

if 0:
    for simname in ['u401']: #['u401','u402','u403']:
        tool = fractal_tool( TL4.loops[simname])
        core_list = None #np.unique(TL4.loops[simname].tr.core_ids)[:5]
        core_list=[323]
        tool.run(core_list=core_list, do_plots=True, plot_prefix=simname)

if 0:
    fig,ax=plt.subplots(1,1)
    for nsim,simname in enumerate(toolshed):
        print(nsim,simname)
        ax.hist(toolshed[simname][0].fractal_dim, histtype='step',color=colors.color[simname])
    outname="%s/fractal_dims.png"%"plots_to_sort"
    fig.savefig(outname)

if 0:
    for nf,ft in enumerate([ft1,ft2,ft3]):
        name=ft.this_looper.out_prefix

        plt.hist(ft.fractal_dim,histtype='step',bins=10,color=colors.color[name],label=name)
        #plt.savefig("plots_to_sort/wtf_%d.png"%nf)

    plt.legend(loc=0)
    plt.savefig('plots_to_sort/fractal_dist.png')

if 0:
    rm = rainbow_map(10)
    for name in toolshed:
        dims = []
        #mylooper=looper_list[0]
        mylooper=toolshed[name]['looper']
        for nframe,frame in enumerate(toolshed[name]):
            tool = toolshed[name][frame]
            lab=ft.this_looper.out_prefix + " %s"%frame
            plt.hist(tool.fractal_dim,histtype='step',bins=10,color=rm(nframe),label=lab)
            dims.append(tool.fractal_dim)
        plt.savefig("plots_to_sort/wtf_%d.png"%nf)
        plt.savefig('plots_to_sort/fractal_dist.png')
        fig,ax=plt.subplots(1,1)
        ax.violinplot(dims)
        fig.savefig('plots_to_sort/violin_test.png')
        plt.close(fig)
        fig,axes=plt.subplots(1,2)
        dims=nar(dims)

        cores_used = tool.cores_used
        OverAchievers = []
        for nb,b in enumerate(dims.transpose()):
            core_id = cores_used[nb]
            if (b > 2.5).any():
                OverAchievers.append(core_id)
            nparticles=mylooper.snaps[0][core_id].target_indices.size
            if nparticles < 10:
                continue

            axes[0].plot(mylooper.tr.times[framelist], b/b[0:3].mean())
            axes[1].plot(mylooper.tr.times[framelist], b)
        fig.savefig('plots_to_sort/streams.png')

#FFF=fractal_dimension(image)

#Z=image
#k=4
#b=np.add.reduceat(Z, np.arange(0, Z.shape[0], k), axis=0)
#b=np.add.reduceat(b, np.arange(0, Z.shape[1], k), axis=1)
#b=np.add.reduceat(b, np.arange(0, Z.shape[2], k), axis=2)

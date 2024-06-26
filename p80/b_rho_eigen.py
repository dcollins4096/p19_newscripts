from starter2 import *

import three_loopers_u900 as TL

from scipy import ndimage

import movie_frames
reload(movie_frames)
import heat_map
reload(heat_map)
plt.close('all')
import pcolormesh_helper as pch
reload(pch)
import progressbar
import time
import eigen
reload(eigen)
def splat(array, tcenter, ax, title,bins):
    smooth= ndimage.gaussian_filter1d(array, 2, 0)
    smooth = 0.5*(smooth[1:,:]+smooth[:-1,:])
    ds_x,ds_y,ds_h,ds_dv,ds_p=heat_map.heat_map( smooth.transpose(), tcenter, bins=bins, ax=ax)
    ax.set_yscale('symlog',linthresh=100)
    ax.set_title(title)
    return ds_x,ds_y,ds_h,ds_dv,ds_p

def plot_monster4(TTT):

    core_id=TTT.core_id
    which = np.where( nar(TTT.collector['cores_used'])==core_id)[0][0]
    fig,ax=plt.subplots(1,2)
    ax0=ax[0];ax1=ax[1]
    ybins = np.geomspace(TTT.rho.min(),TTT.rho.max(),64)/TTT.collector['rho_0'][which]
    Q = TTT.bmag**2/(TTT.collector['B_0'][which])
    #ybins = np.geomspace((TTT.bmag**2).min(),(TTT.bmag**2).max(),64)/TTT.collector['B_0'][which]
    xbins = np.geomspace(Q.min(),Q.max(),64)

    y=TTT.rho.flatten()
    x=Q.flatten()

    pch.simple_phase( x, y, bins=[xbins,ybins],ax=ax0)
    p0=TTT.collector['pfit0'][which]
    p1=TTT.collector['pfit1'][which]
    ax0.set(xscale='log',yscale='log')
    ax0.plot(xbins,xbins, 'k')
    ax0.plot( xbins, 10**(p0*np.log10(xbins)+p1), 'k--')


    y=TTT.rho.flatten()**(1-TTT.collector['RB_mean'][which])
    x=Q.flatten()

    THISAX=ax1  
    pch.simple_phase( x, y, bins=[xbins,ybins],ax=THISAX)
    p0=TTT.collector['pfit0'][which]
    p1=TTT.collector['pfit1'][which]
    THISAX.set(xscale='log',yscale='log', xlabel='B',ylabel='rho',title='1-RB %0.1e'%(1-TTT.collector['RB_mean'][which]))
    THISAX.plot(xbins,xbins, 'k')
    THISAX.plot( xbins, 10**(p0*np.log10(xbins)+p1), 'k--')


    fig.savefig('plots_to_sort/more_test_%s_c%04d'%(TTT.this_looper.sim_name,core_id))




def plot_monster3(TTT):
    print('burp')
    core_id=TTT.core_id

    fig,ax=plt.subplots(2,3)
    T0 = (TTT.A0*TTT.Bp0*TTT.Br0)
    T1 = (TTT.A1*TTT.Bp1*TTT.Br1)
    T2 = (TTT.A2*TTT.Bp2*TTT.Br2)
    R = TTT.R

    ax0=ax[0][0]
    ax1=ax[0][1]
    ax2=ax[0][2]
    ax3=ax[1][0]
    ax4=ax[1][1]
    ax5=ax[1][2]

    ext = extents()
    ext( np.abs(T0));ext(np.abs(T1));ext(np.abs(T2))
    big_bins = np.geomspace(ext.minmax[0],ext.minmax[1],64)
    bins = np.linspace(-2,2,128)
    if 1:
        THIS_AX=ax0
        x = R.flatten()
        y = T0.flatten()
        pch.simple_phase(x,y,ax=THIS_AX,bins=[bins,bins])
        THIS_AX.set(title = 'T0 vs R')
    if 1:
        THIS_AX=ax1
        x = R.flatten()
        y = T1.flatten()
        pch.simple_phase(x,y,ax=THIS_AX,bins=[bins,bins])
        THIS_AX.set(title = 'T1 vs R')
        #THIS_AX.set_xscale('symlog',linthresh=10)
        #THIS_AX.set_yscale('symlog',linthresh=10)
    if 1:
        THIS_AX=ax2
        x = R.flatten()
        y = T2.flatten()
        pch.simple_phase(x,y,ax=THIS_AX,bins=[bins,bins])
        THIS_AX.set(title = 'T2 vs R')
        #THIS_AX.set_xscale('symlog',linthresh=10)
        #THIS_AX.set_yscale('symlog',linthresh=10)

    TheMin=1e-2
    TheMax=1e3
    binsa = np.geomspace(TheMin,TheMax,16)
    if 1:
        THIS_AX=ax3
        x = TTT.rho.flatten()
        y = T0.flatten()
        binsx = np.geomspace(x.min(),x.max(),128)
        print(binsa)
        print(T0[T0>0].min())
        binsy = np.concatenate([-binsa[::-1],binsa])
        pch.simple_phase(x,y,ax=THIS_AX,bins=[binsx,binsy])
        MEAN_T0=(TTT.rho*T0*TTT.dv).sum()/(TTT.rho*TTT.dv).sum()
        THIS_AX.axhline(MEAN_T0)
        THIS_AX.set(title = 'T0 vs rho %0.1f'%MEAN_T0, xscale='log')
        THIS_AX.set_yscale('symlog', linthresh=0.1)
        print(MEAN_T0)
    if 1:
        THIS_AX=ax4
        x = TTT.rho.flatten()
        y = T1.flatten()
        binsx = np.geomspace(x.min(),x.max(),128)
        print(binsa)
        print(T0[T0>0].min())
        binsy = np.concatenate([-binsa[::-1],binsa])
        pch.simple_phase(x,y,ax=THIS_AX,bins=[binsx,binsy])
        MEAN_T1=(np.abs(TTT.rho*T1*TTT.dv)).sum()/(TTT.rho*TTT.dv).sum()
        THIS_AX.axhline(MEAN_T1)
        THIS_AX.set(title = 'T1 vs rho %0.1e'%MEAN_T1, xscale='log')
        THIS_AX.set_yscale('symlog', linthresh=0.1)
        print(MEAN_T1)
    if 1:
        THIS_AX=ax5
        x = TTT.rho.flatten()
        y = T2.flatten()
        binsx = np.geomspace(x.min(),x.max(),128)
        print(binsa)
        print(T0[T0>0].min())
        binsy = np.concatenate([-binsa[::-1],binsa])
        pch.simple_phase(x,y,ax=THIS_AX,bins=[binsx,binsy])
        MEAN_T2=(TTT.rho*T2*TTT.dv).sum()/(TTT.rho*TTT.dv).sum()
        THIS_AX.axhline(MEAN_T2)
        print(MEAN_T2)
        THIS_AX.set(title = 'T1 vs rho %0.1e'%MEAN_T2, xscale='log')
        THIS_AX.set_yscale('symlog', linthresh=0.1)

    if 0:
        THIS_AX=ax0
        srt=np.abs(T0.flatten())
        srt.sort()
        Nbins=32
        MinBin = srt[int(srt.size/Nbins*8)]
        dBin = (np.log10(T0.max())-np.log10(MinBin))/Nbins
        binsa = 10**(np.arange( np.log10(MinBin), np.log10(T0.max()),dBin))
        binsb = -10**(np.arange( np.log10(MinBin), np.log10(-T0.min()),dBin)[::-1])
        bins=np.concatenate([binsb,binsa])
        x = T0.flatten()
        y = T1.flatten()
        pch.simple_phase(x,y,ax=THIS_AX,bins=[bins,bins])
        THIS_AX.set_xscale('symlog')
        THIS_AX.set_yscale('symlog')
    if 0:
        THIS_AX=ax4
        bins = np.linspace(-1,1,128)
        x = (T0+T1+T2).flatten()
        y = T0.flatten()
        pch.simple_phase(x,y,ax=THIS_AX,bins=[bins,bins])
        #THIS_AX.set_xscale('symlog',linthresh=10)
        #THIS_AX.set_yscale('symlog',linthresh=10)
    if 0:
        THIS_AX=ax0
        x=np.abs(T0).flatten()
        y=np.abs(T1).flatten()
        binsx=np.geomspace(x.min(),x.max(),64)
        binsy=np.geomspace(y.min(),y.max(),64)
        pch.simple_phase(x,y,ax=THIS_AX,bins=[binsx,big_bins])
        THIS_AX.plot(binsx,binsx)
        THIS_AX.set(xscale='log',yscale='log',title='T0 vs T1')
    if 0:
        THIS_AX=ax0
        x=np.abs(T0).flatten()
        y=np.abs(T0+T1+T2).flatten()
        binsx=np.geomspace(x.min(),x.max(),64)
        binsy=np.geomspace(y.min(),y.max(),64)
        pch.simple_phase(x,y,ax=THIS_AX,bins=[binsx,big_bins])
        THIS_AX.plot(binsx,binsx)
        THIS_AX.set(xscale='log',yscale='log',title='R vs T0')
    if 0:
        THIS_AX=ax1
        x=np.abs(T1).flatten()
        y=np.abs(T0+T1+T2).flatten()
        binsx=np.geomspace(x.min(),x.max(),64)
        binsy=np.geomspace(y.min(),y.max(),64)
        pch.simple_phase(x,y,ax=THIS_AX,bins=[binsx,big_bins])
        THIS_AX.plot(binsx,binsx)
        THIS_AX.set(xscale='log',yscale='log',title='R vs T1')
    fig.savefig('plots_to_sort/the_T_%s_c%04d'%(TTT.this_looper.sim_name,core_id))

def plot_monster2(TTT):
    print('burp')
    core_id=TTT.core_id
    if 1:
        fig, ax=plt.subplots(4,6,figsize=(20,12))
        ax0=ax[0][0]
        ax1=ax[0][1]
        ax2=ax[0][2]
        ax3=ax[0][3]
        ax4=ax[0][4]
        ax5=ax[0][5]

        ax6=ax[1][0]
        ax7=ax[1][1]
        ax8=ax[1][2]
        ax9=ax[1][3]
        ax10=ax[1][4]
        ax11=ax[1][5]

        ax12=ax[2][0]
        ax13=ax[2][1]
        ax14=ax[2][2]
        ax15=ax[2][3]
        ax16=ax[2][4]
        ax17=ax[2][5]

        ax18=ax[3][0]
        ax19=ax[3][1]
        ax20=ax[3][2]
        ax21=ax[3][3]
        ax22=ax[3][4]
        ax23=ax[3][5]
        def DODO(V0,V1,V2,a0,a1,a2,a3, log=False):
            if log:
                binner=np.geomspace
            else:
                binner=np.linspace
            if a0:
                X = V0.flatten()
                bins=binner(X.min(),X.max(),32)
                a0.hist(V0.flatten(), histtype='step',bins=bins, cumulative=False, density=False, label='V0')
                a0.hist(V1.flatten(), histtype='step',bins=bins, cumulative=False, density=False, label='V1')
                a0.hist(V2.flatten(), histtype='step',bins=bins, cumulative=False, density=False, label='V2')
                #a0.legend(loc=0)
                #a0.set(xscale='log',yscale='log',title='P(A0)')
            if a1:
                #P(A1A0)
                X =(V1/V0).flatten()
                bins=binner(X.min(),X.max(),32)
                a1.hist(X, histtype='step',bins=bins, cumulative=True, density=True,label='V1/V0')
                X =(V2/V0).flatten()
                bins=binner(X.min(),X.max(),32)
                a1.hist(X, histtype='step',bins=bins, cumulative=True, density=True,label='V2/V0')
                #ax1.legend(loc=0)
                #ax1.set(xscale='log',yscale='linear')#,title='P(A1/A0)')
            if a2:
                #A0 vs A1
                X =(V0.flatten())
                Y =(V1.flatten())
                THIS_AX=a2
                ext=extents(Y)
                pch.simple_phase(X,Y,ax=THIS_AX,log=log)
                #THIS_AX.set(xscale='log',yscale='log', title='A1 vs A0')
                THIS_AX.plot(ext.minmax,ext.minmax,c=[0.5]*3)
            if a3:
                #A0 vs A2
                X =(V0.flatten())
                Y =(V2.flatten())
                THIS_AX=a3
                ext=extents(Y)
                pch.simple_phase(X,Y,ax=THIS_AX,log=log)
                #THIS_AX.set(xscale='log',yscale='log', title='A2 vs A0')
                THIS_AX.plot(ext.minmax,ext.minmax,c=[0.5]*3)
        DODO( TTT.A0, TTT.A1, TTT.A2, ax0, ax1, ax2, ax3, log=True)
        ax0.set(xscale='log',yscale='log',title=r'$P(a_i)$')
        ax0.legend(loc=0)
        ax1.set(xscale='log',yscale='linear',title=r'$P(a_{1,2}/a_0)$')
        ax2.set(xscale='log',yscale='log',title='a1 vs a0')
        ax3.set(xscale='log',yscale='log',title='a2 vs a0')
        DODO( TTT.Bp0, TTT.Bp1, TTT.Bp2, ax4, None, None, None, log=False)
        ax4.set(title=r'P(b_i) new rot')
        ax4.legend(loc=0)
        DODO( TTT.Br0, TTT.Br1, TTT.Br2, ax5, None, None, None, log=False)
        ax5.set(title=r'P(b_i) new')
        ax5.legend(loc=0)
        T0 = (TTT.A0*TTT.Bp0*TTT.Br0)
        T1 = (TTT.A1*TTT.Bp1*TTT.Br1)
        T2 = (TTT.A2*TTT.Bp2*TTT.Br2)

        def GOO(A0,Bp0,Br0,n,nbins, ax6,ax7,ax8,ax9,ax10,ax11):

            THIS_AX=ax6
            x = A0.flatten()
            y = Bp0.flatten()
            binsx = np.geomspace( x.min(),x.max(),64)
            binsy = np.linspace(y.min(),y.max(),64)
            pch.simple_phase(x,y,ax=THIS_AX,bins=[binsx,binsy])
            THIS_AX.set(xscale='log',title='a0 vs bp')
        
            THIS_AX=ax7
            x = A0.flatten()
            y = Br0.flatten()
            binsx = np.geomspace( x.min(),x.max(),64)
            binsy = np.linspace(y.min(),y.max(),64)
            pch.simple_phase(x,y,ax=THIS_AX,bins=[binsx,binsy])
            THIS_AX.set(xscale='log',title='a0 vs br')
        
            THIS_AX=ax8
            x = Bp0.flatten()
            y = Br0.flatten()
            binsx = np.linspace( x.min(),x.max(),64)
            binsy = np.linspace(y.min(),y.max(),64)
            pch.simple_phase(x,y,ax=THIS_AX,bins=[binsx,binsy])
            THIS_AX.set(xscale='linear',title='bp vs br')
        
            THIS_AX=ax9
            x = A0.flatten()
            y = Br0.flatten()*Bp0.flatten()
            binsx = np.geomspace( x.min(),x.max(),64)
            binsy = np.linspace(y.min(),y.max(),64)
            pch.simple_phase(x,y,ax=THIS_AX,bins=[binsx,binsy])
            THIS_AX.set(xscale='log',title='Br0*Bp0 vs A0')
        
            THIS_AX=ax10
            #x = A0.flatten()
            y = Br0.flatten()*Bp0.flatten()
            #binsx = np.geomspace( x.min(),x.max(),64)
            binsy = np.linspace(y.min(),y.max(),64)
            THIS_AX.hist(y, histtype='step',bins=binsy)
            THIS_AX.set(xscale='linear',title='Br0*Bp0')
            if 1:
                #T0 vs Br0*Bp0
                THIS_AX=ax11
                x = Br0.flatten()*Bp0.flatten()
                y = (A0*Bp0*Br0).flatten()
                ok = slice(None)
                x = x #np.abs(x[ok])
                y =np.abs(y[ok])
                if (y<0).sum():
                    pdb.set_trace()
                binsx = np.linspace( x.min(),x.max(),64)
                #binsy = np.geomspace(y.min(),y.max(),64)
                binsy = np.geomspace(5e-6,1e3,64)
                pch.simple_phase(x,y,ax=THIS_AX,bins=[binsx,binsy])
                #THIS_AX.plot(x,x)
                THIS_AX.set(xscale='linear',yscale='log',title='T0 vs Br0*Bp0')

            if 0:
                #Br0 vs Br0*Bp0
                THIS_AX=ax11
                x = Br0.flatten()
                y = Br0.flatten()*TTT.Bp0.flatten()
                binsx = np.linspace( x.min(),x.max(),64)
                binsy = np.linspace(y.min(),y.max(),64)
                pch.simple_phase(x,y,ax=THIS_AX,bins=[binsx,binsy])
                THIS_AX.set(xscale='linear',title='Br0*Bp0 vs A0')
            if 0:
                #Product of two uniforms is -logz
                x = np.random.random(1000)
                y = np.random.random(1000)
                ax11.hist(x,histtype='step',density=True,bins=32)
                ax11.hist(y,histtype='step',density=True,bins=32)
                ax11.hist(x*y,histtype='step',density=True,bins=32)
                xx=np.linspace(1e-3,1,100)
                ax11.plot(xx,-np.log(xx))
                #ax11.plot(xx,xx-xx*np.log(xx))
        GOO(TTT.A0,TTT.Bp0,TTT.Br0,0,64,ax6,ax7,ax8,ax9,ax10,ax11)
        GOO(TTT.A1,TTT.Bp1,TTT.Br1,0,64,ax12,ax13,ax14,ax15,ax16,ax17)
        GOO(TTT.A2,TTT.Bp2,TTT.Br2,0,64,ax18,ax19,ax20,ax21,ax22,ax23)




        #pdb.set_trace()
        #DODO( T0, T1, T2, ax12, ax13, ax14, ax15, log=False)

        if 0:
            MinBin=np.max([T0[T0>0].min(),-T0[T0<0].max()])
            srt=np.abs(T0.flatten())
            srt.sort()
            Nbins=16
            MinBin = srt[int(srt.size/Nbins*4)]
            dBin = (np.log10(T0.max())-np.log10(MinBin))/Nbins
            binsa = 10**(np.arange( np.log10(MinBin), np.log10(T0.max()),dBin))
            binsb = -10**(np.arange( np.log10(MinBin), np.log10(-T0.min()),dBin)[::-1])
            bins=np.concatenate([binsb,binsa])
            print(bins)
            ax12.hist(T0.flatten(),histtype='step',bins=bins,density=True)
            ax12.hist(T1.flatten(),histtype='step',bins=bins,density=True)
            ax12.hist(T2.flatten(),histtype='step',bins=bins,density=True)
            ax12.set_xscale('symlog', linthresh=1)
            ax12.set(yscale='log')

        if 0:
            #A1 vs A2
            X =(TTT.A1.flatten())
            Y =(TTT.A2.flatten())
            ext=extents(X)
            pch.simple_phase(X,Y,ax=ax10,log=True)
            ax10.set(xscale='log',yscale='log',xlabel='A1',ylabel='A2')
            ax10.plot(ext.minmax,ext.minmax,c=[0.5]*3)
        if 0:
            T0 = (TTT.A0*TTT.Bp0*TTT.Br0)
            T1 = (TTT.A1*TTT.Bp1*TTT.Br1)
            T2 = (TTT.A2*TTT.Bp2*TTT.Br2)
            #X =np.abs((TTT.R.flatten()))
            X =np.abs((T0.flatten()))
            Y =np.abs((T2).flatten())
            print(Y.shape)
            print(X.shape)
            print(np.isnan(X).any())
            print(np.isnan(Y).any())
            ext=extents(Y)
            pch.simple_phase(X,Y,ax=ax11,log=True)
            ax11.set(xscale='log',yscale='log',xlabel='A0',ylabel='A1',title='R vs the first')
            ax11.plot(ext.minmax,ext.minmax,c=[0.5]*3)
        fig.savefig('plots_to_sort/the_ratios_%s_c%04d'%(TTT.this_looper.sim_name,core_id))


def plot_monster1(TTT):
    core_id=TTT.core_id
    if 1:

        fig, ax=plt.subplots(5,3,figsize=(20,12))
        ax0=ax[0][0]
        ax1=ax[0][1]
        ax2=ax[0][2]
        ax3=ax[1][0]
        ax4=ax[1][1]
        ax5=ax[1][2]
        ax6=ax[2][0]
        ax7=ax[2][1]
        ax8=ax[2][2]

        ax9= ax[3][0]
        ax10=ax[3][1]
        ax11=ax[3][2]
        ax12=ax[4][0]
        ax13=ax[4][1]
        ax14=ax[4][2]
        ext=extents()
        #ext(A0);ext(A1);ext(A2)
        ext(TTT.A0);ext(TTT.A1);ext(TTT.A2)
        #print(A0.min())
        #bins2=np.geomspace(ext.minmax[0],ext.minmax[1],64)
        if ext.minmax[0]<0:
            binsa = np.geomspace( np.abs(TTT.A0)[np.abs(TTT.A0)>0].min(), TTT.A0.max(), 32)
            bins = np.concatenate([-binsa[::-1],binsa])
        else:
            bins=np.geomspace(ext.minmax[0],ext.minmax[1],64)


        if 1:
            splat(TTT.A0, TTT.tcenter, ax0, 'A0',bins=bins)
            splat(TTT.A1, TTT.tcenter, ax1, 'A1',bins=bins)
            splat(TTT.A2, TTT.tcenter, ax2, 'A2',bins=bins)
            #ax0.set(yscale='log')#,ylim=ext.minmax)
            #ax1.set(yscale='log')#,ylim=ext.minmax)
            #ax2.set(yscale='log')#,ylim=ext.minmax)
            ax0.set_yscale('symlog')
            ax1.set_yscale('symlog')
            ax2.set_yscale('symlog')

        if 1:

            ext=extents()
            ext(TTT.Br0);ext(TTT.Br1);ext(TTT.Br2)
            bins=np.linspace(ext.minmax[0],ext.minmax[1],64)

            splat(TTT.Bp0, TTT.tcenter, ax3, 'B_r_n 0',bins=bins)
            splat(TTT.Bp1, TTT.tcenter, ax4, 'B_r_n 1',bins=bins)
            splat(TTT.Bp2, TTT.tcenter, ax5, 'B_r_n 2',bins=bins)
            #ax3.set(yscale='linear',ylim=ext.minmax)
            #ax3.set_yscale('symlog',linthresh=1e-2)
            ax3.set(yscale='linear',ylim=ext.minmax)
            ax4.set(yscale='linear',ylim=ext.minmax)
            ax5.set(yscale='linear',ylim=ext.minmax)

            splat(TTT.Br0, TTT.tcenter, ax6, 'B_n 0',bins=bins)
            splat(TTT.Br1, TTT.tcenter, ax7, 'B_n 1',bins=bins)
            splat(TTT.Br2, TTT.tcenter, ax8, 'B_n 2',bins=bins)
            ax6.set(yscale='linear',ylim=ext.minmax)
            ax7.set(yscale='linear',ylim=ext.minmax)
            ax8.set(yscale='linear',ylim=ext.minmax)


        if 1:
            bins=np.linspace(-100,100)
            bins=np.linspace(-1,1,128)
            T0 = (TTT.A0*TTT.Bp0*TTT.Br0)
            splat(T0,TTT.tcenter, ax9,'T0',bins=bins)
            ax9.set(yscale='linear')
            T1 = (TTT.A1*TTT.Bp1*TTT.Br1)
            splat(T1,TTT.tcenter, ax10,'T1',bins=bins)
            ax10.set(yscale='linear')
            T2 = (TTT.A2*TTT.Bp2*TTT.Br2)
            splat(T2,TTT.tcenter, ax11,'T2',bins=bins)
            ax11.set(yscale='linear')
        if 1:
            binsa = np.geomspace( T0[T0>0].min(),np.abs(T0).max(),16)
            bins = np.concatenate([-binsa[::-1],binsa])
            splat(T0,TTT.tcenter, ax12,'T0',bins=bins)
            #ax12.set(yscale='linear')
            ax12.set_yscale('symlog')
            splat(T1,TTT.tcenter, ax13,'T1',bins=bins)
            ax13.set_yscale('symlog')
            splat(T2,TTT.tcenter, ax14,'T2',bins=bins)
            ax14.set_yscale('symlog')

        if 0:
            X =np.abs((TTT.R.flatten()))
            Y =np.abs((T0+T1+T2).flatten())
            print(Y.shape)
            print(X.shape)
            print(np.isnan(X).any())
            print(np.isnan(Y).any())
            ext=extents(Y)
            pch.simple_phase(X,Y,ax=ax11,log=True)
            ax11.set(xscale='log',yscale='log',xlabel='A0',ylabel='A1',title='word')
            ax11.plot(ext.minmax,ext.minmax,c=[0.5]*3)




        fig.savefig('plots_to_sort/the_stuff_%s_c%04d'%(TTT.this_looper.sim_name,core_id))





sim_list=['u902']
new_ddd=False
if 'clobber' not in dir():
    clobber = False
if 'ddd' not in dir() or clobber:
    for sim in sim_list:
        ddd = eigen.dq_dt2(TL.loops[sim])
        core_list=None
        #core_list=[7]
        core_list=[74]#, 112]
        #core_list=[112]
        #core_list=[74, 112]
        #core_list=[112]
        core_list=[74]
        #core_list=TL.loops[sim].core_by_mode['Alone']
        #core_list=core_list[8:9]
        P=ddd.run(core_list=core_list)
        new_ddd=True

if not new_ddd:
    for sim in sim_list:
        TTT = ddd
        if 1:
            if 1:
                #plot_monster1(TTT)
                #plot_monster2(TTT)
                plot_monster3(TTT)
                #plot_monster4(TTT)
                #plot_monster5(TTT)


from seabreezediag.configdir import Config, Meta
from netCDF4 import Dataset as nc,num2date,date2num
import numpy as np
import os
def get_meta(fn, vn, M):
    '''Return the sturcture of data that is collected'''
    with nc(fn) as f:
        tn=f.variables[vn].dimensions[0] #The name of the time variable
        times = num2date(f.variables[tn][0:2],f.variables[tn].units)
    dt = (times[1]-times[0]).seconds
    secperday=86400.
    ts=int(secperday/dt)
    D = np.zeros([ts,len(M.lat),len(M.lon)])
    data={}
    for v in ('sb_con','thc',vn):
        data[v]={}
        for s in ('DJF','MAM','JJA','SON'):
            data[v][s]=np.zeros_like(D)
    return ts,data

def get_data(fn,vn):
    '''Return the data in a netcdf array'''

    with nc(fn) as f:
        tn=f.variables[vn].dimensions[0] #The name of the time variable
        times = num2date(f.variables[tn][:],f.variables[tn].units)
        dt = int((times[1]-times[0]).seconds)
        #Get the number of days in the array
        nday = int(dt*len(times)/86400)
        shape=f.variables[vn].shape
        data = np.mean(f.variables[vn][:].reshape(int(nday),int(86400/dt),shape[-2],shape[-1]),axis=0)
    return data

def main(**kwargs):
    config = kwargs['config']
    Cfg = Config(kwargs['config'])
    M = Meta(Cfg)
    M.get_dates()
    tt = {}
    ncout = os.path.join(M.datadir,'%ssb_con.nc'%M.prefix)
    mon2seas={'01':'DJF','02':'DJF','03':'MAM','04':'MAM','05':'MAM','06':'JJA','07':'JJA','08':'JJA','09':'SON','10':'SON','11':'SON','12':'DJF'}
    seas=dict(DJF=1,MAM=4,JJA=6,SON=10)
    if not os.path.isfile(ncout):
        for nn,ts in enumerate(M.dates):
            year,mon=tuple(ts.split('_'))
            tstring = os.path.join(M.datadir,year,'%sXXX_%s.nc'%(M.prefix,ts))
            f_sb=tstring.replace('XXX','sb')
            sys.stdout.flush()
            sys.stdout.write('Adding information from %s/%s to %s... '%(mon,year,mon2seas[mon]))
            sys.stdout.flush()
            vars={'sb_con':f_sb,'thc':f_sb}
            vars[Cfg['vtheta']]=tstring.replace('XXX',Cfg['vtheta'])
            if nn == 0:
                nt,data=get_meta(vars[Cfg['vtheta']],Cfg['vtheta'],M)
                times = np.arange(0,24,24/int(nt))*60**2
                for s in seas.keys():
                    tt[s]=1
            for v in vars.keys():
                data[v][mon2seas[mon]] += get_data(vars[v],v)
            tt[mon2seas[mon]] += 1
            sys.stdout.write('ok\n')
        vars=list(vars.keys())
        for v in vars:
            shape=data[v][list(seas.keys())[0]].shape
            D = np.zeros([len(seas.keys())*shape[0],shape[1],shape[2]])
            #T=np.zeros([leddn(seas.keys())*len(times)])
            T=[]
            dd = 0
            for s in seas.keys():
                if v == Cfg.vtheta:
                    vn = 'temp'
                else:
                    vn = v
                D[dd*len(times):(dd*len(times))+len(times)]=data[v][s]/tt[s]
                for ii in times:
                    T.append(num2date(ii,'Seconds since 2017-%02i-15 00:00:00'%seas[s]))
                #T[dd*len(times):(dd*len(times))+len(times)]=num2date(times,'Seconds since 2017-%02i-15 00:00:00'%seas[s])
                dd += 1
            M.create_nc(D,ncout,vn,T)
    import matplotlib.pyplot as plt
    from mpl_toolkits.basemap import Basemap, cm
    cmap = cm.GMT_polar
    cmap2= cm.sstanom
    mon2seas={'01':'DJF','04':'MAM','06':'JJA','10':'SON'}
    with nc(ncout) as f:
        tt = 0
        sys.stdout.flush()
        sys.stdout.write('Createing map ... ')
        m = Basemap(lon_0 = 10, resolution = 'i', area_thresh = 10000)
        sys.stdout.write('ok\n')
        fig, axes = plt.subplots(2, 1,figsize=(14.22,8.875))
        plot = os.path.join(os.path.expanduser('~/PhD/UM/'),'test_run_ERAI-temp.png')
        nf = '/home/wilfred/Data/ERAI/netcdf2/1987/Erai_sb_1987_08.nc'
        g1 = nc(nf)
        g2 = nc(nf.replace('_sb_','_t2m_'))
        axes[0].set_title("2M Temperature ")
        axes[1].set_title("Sea-Level Temperature (Moist Adiabatic Descent)")
        im = m.pcolormesh(M.lon,M.lat,g2.variables['t2m'][0]-273.15,latlon=True,vmin=-37.,vmax=37.,ax=axes[0],cmap=cmap2)
        im = m.pcolormesh(M.lon,M.lat,g1.variables['thc'][0]-273.15,latlon=True,vmin=-37.,vmax=37.,ax=axes[1],cmap=cmap2)
        m.drawcoastlines(ax=axes[0],linewidth=0.25)
        m.drawcoastlines(ax=axes[1],linewidth=0.25)
        plt.subplots_adjust(bottom=0.1, right=0.98, top=0.98,left = 0.01,hspace=0.10,wspace=0.10)
        #cax = plt.axes([0.01, 0.05, 0.99, 0.15])
        #m.colorbar(cax=cax,location='bottom',ax=axes[0,0])
        fig.colorbar(im,ax=axes.ravel().tolist(),shrink=0.95,pad=0.01, aspect=80,label='Temperature ($^{\\circ}$C)')
        fig.savefig(plot,dpi=150,bbox_inches='tight')
        g1.close(),g2.close()

        for mon,seas in mon2seas.items():
            fig, axes = plt.subplots(2, 2,figsize=(14.22,8.875))
            plot = os.path.join(os.path.expanduser('~/PhD/UM/'),'test_run_ERAI-%s.png'%(seas))
            sys.stdout.flush()
            sys.stdout.write('Plotting data for %s ... '%seas)
            hours=range(0,24,6)
            h=0
            for ii in range(0,2):
                for jj in range(0,2):
                    data = f.variables['sb_con'][tt]
                    axes[ii,jj].set_title("Subgrid Sea-Breeze Convergence at %02i UTC (%s)"%(hours[h],seas))
                    im = m.pcolormesh(M.lon,M.lat,f.variables['sb_con'][tt],latlon=True,vmin=-5.,vmax=5.,ax=axes[ii,jj],cmap=cmap)
                    m.drawcoastlines(ax=axes[ii,jj],linewidth=0.25)
                    h += 1
                    tt += 1
            plt.subplots_adjust(bottom=0.2, right=0.98, top=0.98,left = 0.01,hspace=0.00,wspace=0.00)
            #cax = plt.axes([0.01, 0.05, 0.99, 0.15])
            #m.colorbar(cax=cax,location='bottom',ax=axes[0,0])
            fig.colorbar(im,ax=axes.ravel().tolist(),shrink=0.95,pad=0.01, aspect=80)
            fig.savefig(plot,dpi=150,bbox_inches='tight')
            sys.stdout.write('ok\n')
            fig.clf()



def kw(**kwargs):
    """
    Method that returns the default parameter keywordarguments
    """
    path=os.path.join(os.path.dirname(os.path.abspath(sys.argv[0])))
    KW=dict(config=os.path.join(path,'run.conf'))
    for key,value in kwargs.items():
        KW[key]=value
    return KW



if __name__ == "__main__":
    import sys
    kwargs=kw()
    helpstring="""
    Module to run the seabreeze detection.

    Usage:
        python %s --option=value
        options are given below

    Options:
        config           : The configuration file that contains all information
                           to run the sea-breeze detection %s
            """%(
            sys.argv[0],
            str(kwargs['config'])\
        )
    try:
        for arg in sys.argv[1:]:
            try:
                key,value = arg.strip('--').split('=')
            except ValueError:
                sys.exit(helpstring)
            if key.lower()=='help':
                sys.exit(helpstring)
            #Everything else shold go here
            else :
                try: #The key is an integer set it
                    kwargs[key.lower()] = float(value)
                except ValueError:
                    try: #Do we have a list of integers seperated by kommas?
                        kwargs[key.lower()] = [float(v) for v in value.split(',')]
                    except ValueError: #If not, probably a string
                        kwargs[key.lower()] = value

    except IndexError:
        pass
    kwargs = kw(**kwargs)
    main(**kwargs)


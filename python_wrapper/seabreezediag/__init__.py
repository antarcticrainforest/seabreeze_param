import numpy as np
import warnings
from seabreeze import diag as dg, get_edges, get_dist
"""
Python wrapper for seabreeze diag fortran routine

Methods
-------
    diag : call the sea-breeze diagnosis fortran subroutine
    c2f  : convert column/row major arrays
    read : read data from a netcdf-file
"""
def c2f(array):
    '''
    Convert from column/row major to row/column major
    
    Parameters
    ----------
    array : ND - array
            row or column major array
    Returns
    -------
    ND - array row or column major array

    Example
    -------
    Create a c style array and make it a fortran style one
    >>> import numpy as np
    >>> c_style = np.random_sample([2,3])
    >>> c_style
    ... array([[ 0.3542732 ,  0.80694807,  0.56701353],
               [ 0.31814638,  0.95021108,  0.42089208]])
    >>> f_style = c2f(c_style)
    >>> f_style
    ... array([[ 0.3542732 ,  0.31814638],
               [ 0.80694807,  0.95021108],
               [ 0.56701353,  0.42089208]])
    '''
    sh = array.shape[::-1]
    return array.ravel().reshape(sh,order='F')

def check(kw,key,value):
    ''' Check for the presence in a dictionary and return the value, and if not
        present return a default value
    '''
    try:
        value = kw[key]
        del kw[key]
    except KeyError:
        pass
    return value

def read_nc(fnv, fnu, fntheta, fnci,vv='v', vu='u', vtheta='t2m', vci='ci', 
            vpres='pres', vtime='time'):
    """
      Function to read meta data from various netcdf-sources

      Variables:
        fnv (str)     : name of the v-wind netcdf file
        fnu (str)     : name of the u-wind netcdf file
        fntheta (str) : name of the temp netcdf file
        fnci (str)    : name of the sea-ice frac. netcdf file
        vv (str)      : name of the v-wind netcdf variable (default : v)
        vu (str)      : name of the u-wind netcdf variable (default : u)
        vtheta (str)  : name of the of the surf temp var. (default : t2m)
        vpres (str)   : name of the of the presure var. (default : pres)
        vci (str)     : name of the sea-ice variable (default : ci)

       Returns named tuple:
    """
    from netCDF4 import Dataset as nc, num2date
    from collections import namedtuple
    import os
    data = dict( v=vv, u=vu, theta=vtheta, ci=vci, pres=vpres)
    meta = namedtuple('meta','time dt nc '+' '.join(data.keys()))
    get_meta = True
    usr = os.path.expanduser
    meta.nc = {'v':nc(usr(fnv)),'u':nc(usr(fnu)),'theta':nc(usr(fntheta))
               ,'ci':nc(usr(fnci))}

    for v, ncf in meta.nc.items():
        setattr(meta,v,ncf.variables[data[v]])

    dims = meta.nc['v'].dimensions
    meta.time = num2date(meta.nc['v'].variables[vtime][:],
                         meta.nc['v'].variables[vtime].units)
    meta.pres = meta.nc['v'].variables[vpres][:]
    meta.dt = (meta.time[1]-meta.time[0]).seconds/60.
    return meta

def diag(tt, lsm, z, std, lon, lat, pres, *args, **kwargs):
    """
    Calculate potential strength of sea-breeze convergence within a predefined
    coastal area (see Bergemann et al. 2017 doi:10.1002/2017MS001048 for detail)

    Parameters
    ----------
    tt: int
        Number of timestep since start of application, used to determine the
        beginning of simulation and and the application of the algorithm in the
        considered interval
    lsm: 2-D array of type-float
        Land-sea mask of the model (lat,lon)
    z  : 2-D array of type-float
        Surface elevation data (lat,lon)
    std : 2-D array of type-float
        Standard deviation of subgrid orography
    lon: 1-D array of type-float
        the longitude vector of the model
    lat: 1-D array of type-float
        the latitude vector of the model
    pres: 1-D array of type-float
        Pressure levels in Pa stored in a 1D array
    u: N-D array of type-float
        U-component of the wind  ([time],pres,lat,lon)
        time is optional
    v: N-D array of type-float
        V-component of the wind stored ([time],pres,lat,lon)
        time is optional
    t: N-D array of type-float
        Surface temperature array  ([time],lat,lon)
        time is optional
    ci: N-D array of type-float
        Fraction of sea-ice cover [0..1] ([time],lat,lon)
        time is optional. Ci can be set to None, in which case it won't be
        taken into account for calculation.

    Keywords
    --------
    windspeed: 2-D array of type-float
        Wind speed, passed through the application period (lat,lon)
        Default: zero array
    winddir: 2-D array of type-float
        Wind direction, is passed through the application period (lat,lon)
        Default: zero array
    thc: 2-D array of type-float
        Thermal heating contrast, passed through the application period (lat,lon)
        Default: zero array
    target_plev : float
        The pressure level in hPa for wind threshold application Default: 700.0
    thresh_wind :  float
        Threshold of wind speed in m/s Default: 11
    thresh_winddir :  float
        Threshold of change in wind direction in deg Default: 90.0
    thresh_windch :  float
        Threshold of change in wind speed in m/s Default: 5.0
    thresh_thc :  float
        Threshold in of thermal heating contrast in K Default: 0.75
    target_time :  float
        Time period for the application of the thresholds in h Default: 6
    maxdist :  float
        Distance from the coast that is supposed to be influenced by
        sea-breeezes in km Default: 180
    timestep :  float
        Time stepping of the input data in minutes Default: 24

    Returns
    -------
        timestep number: integer indicating the timestep since the start of
                        the application of diagnostic
        sea-breeze-strength,
        thermal-heating-contrast, windspeep wind-direction

    Example
    -------
    >>> from netCDF4 import Dataset as nc, num2date
    >>> import numpy as np
    >>> from datetime import datetime, timedelta
    >>> import seabreezediag as sb
    >>> tt = 1
    >>> THC = np.zeros([T.shape[-3],T.shape[-2]])
    >>> WS,WD = np.zeros_like(thc),np.zeros_like(thc)
    >>> for fn in  ('input_netcdf-month1.nc',input_netcdf-month2.nc'):
            test_nc = nc(fn,'a') #Open the data netcdf-file
            time = num2date(test_nc.variables['time'][:],
            ...                 test_nc.variables['time'].units)
            >>> P = test_nc.variables['pres'][:]
            >>> U = test_nc.variables['uwind'][:]
            >>> V = test_nc.variables['vwind'][:]
            >>> T = test_nc.variables['tsurf'][:]
            >>> LSM = test_nc.variables['lsm'][:]
            >>> CI = test_nc.variables['s_ice'][:]
            >>> H = test_nc.variables['z'][:]
            >>> dt = (time[1] - time[0]).minutes
            >>> tt,breeze,THC,WS,WD = sb.diag(tt,LSM,H,CI,P,U,V,T,CI,WS,WD,
            ...                                    THC,timestep=dt)
            >>> test_nc.variables['breeze'][:]=breeze
            >>> test_nc.close()
    """
    import os,sys
    ws = check(kwargs,'ws',None)
    wd = check(kwargs,'wd',None)
    thc = check(kwargs,'thc',None)
    meta= check(kwargs,'meta',None)
    if type(meta) == type(None):
        u, v, t, ci = args
    else:
        u, v, t = meta.u, meta.v, meta.theta
        try:
            ci = meta.ci
        except AttributeError:
            ci = None
    tt = max(1,tt)
    if type(ws) == type(None):
        if tt > 1 :
            warnings.warn('Windspeed should be given from previous timestep')
        ws = np.zeros_like(lsm)
    if type(wd) == type(None):
        if tt > 1 :
            warnings.warn('Wind direction should be given from previous timestep')
        wd = np.zeros_like(lsm)
    if type(thc) == type(None):
        if tt > 1 :
            warnings.warn('Heating contrast should be given from previous timestep')
        thc = np.zeros_like(lsm)
    if type(ci) == type(None):
        dist = c2f(get_dist(get_edges(c2f(lsm),c2f(np.zeros_like(c))),
            c2f(lsm),lon,lat))
    if len(v.shape) > 3:
        #Time is also contained
        sb_con=np.zeros([4,len(v),t.shape[-2],t.shape[-1]])
        for ts in range(len(v)):
            if type(ci) != type(None):
                try:
                    ic = ci[ts].filled(0)
                except AttributeError:
                    ic = ci[ts]
                dist = c2f(get_dist(get_edges(c2f(lsm),c2f(ic)),c2f(lsm),lon,lat))
            out = c2f(dg(tt,
                         c2f(pres),
                         c2f(z),
                         c2f(std),
                         c2f(t[ts]),
                         c2f(v[ts]),
                         c2f(u[ts]),
                         c2f(dist),
                         c2f(ws),
                         c2f(wd),
                         c2f(thc),**kwargs))
            sb_con[0,ts] = out[0]
            sb_con[1,ts] = out[1]
            sb_con[2,ts] = out[2]
            sb_con[3,ts] = out[3]
            thc,ws,wd = out[1],out[2],out[3]
            tt += 1
    else:
        if type(ci) != type(None):
            try:
                ic = ci.filled(0)
            except AttributeError:
                ic = ci
            dist = c2f(get_dist(get_edges(c2f(lsm),c2f(ic)),c2f(lsm),lon,lat))
        sb_con=np.zeros([4,1,t.shape[-2],t.shape[-1]])
        out = c2f(dg(tt,c2f(p),c2f(z),c2f(std),c2f(t[:]),c2f(v[:]),c2f(u[:]),
            c2f(dist),c2f(ws),c2f(wd),c2f(thc),**kwargs))
        sb_con[0] = out[0]
        sb_con[1] = out[1]
        sb_con[2] = out[2]
        sb_con[3] = out[3]
        thc,ws,wd = out[1],out[2],out[3]
        tt += 1
    del out
    return tt,sb_con[0],thc,ws,wd


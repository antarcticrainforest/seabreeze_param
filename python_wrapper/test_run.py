from netCDF4 import Dataset as nc, date2num,num2date
from datetime import timedelta, datetime
import seabreezediag as sbd
from seabreezediag.configdir import Config, Meta
import numpy as np
import os,glob

def main(**kwargs):
    config = kwargs['config']
    Cfg = Config(kwargs['config']) #Read the config file
    M = Meta(Cfg) #Create the meta object
    #Initialize the heating contr., windspeed, wind dir. and output arrays
    thc = np.zeros([len(M.lat),len(M.lon)])
    windspeed = np.zeros_like(thc)
    winddir = np.zeros_like(thc)
    sb_con = np.zeros_like(thc)
    tt = 1 #Define the first timestep
    for ts in M.dates:
        #Ts is of form of 'YYYY_MM' or 'YYYY_MM_DD'
        year=ts.split('_')[0]
        tstring = os.path.join(M.datadir,year,'%sXXX_%s.nc'%(M.prefix,ts))
        f_sb=tstring.replace('XXX', 'sb') #Define the output netcdf-file
        sys.stdout.flush()
        sys.stdout.write('Creating sea-breeze data for %s ... '%os.path.basename(f_sb))
        sys.stdout.flush()
        vars={}
        for v in ('vu', 'vv', 'vtheta', 'vci'):
            vars[v]=tstring.replace('XXX', Cfg[v])
        # Read the netcdf-file containing the input data
        nc_data = sbd.read_nc(vars['vv'],
                              vars['vu'],
                              vars['vtheta'],
                              vars['vci'],
                              vv=Cfg.vv,
                              vu=Cfg.vu,
                              vpres=Cfg.vpres,
                              vtime=Cfg.vtime,
                              vtheta=Cfg.vtheta)
        tt, sb_con, thc, windspeed, winddir = sbd.diag(tt,
                                                       M.landfrac,
                                                       M.z,
                                                       M.std,
                                                       M.lon,
                                                       M.lat,
                                                       nc_data.pres,
                                                       meta=nc_data,
                                                       ws=windspeed,
                                                       wd=winddir,
                                                       thc=thc)

        M.create_nc(sb_con,f_sb,'sb_con', nc_data.time)
        #Close all opened netcdf files
        for f in nc_data.nc.values():
            f.close()
        del sb_con
        sys.stdout.flush()
        sys.stdout.write('ok\n')
        sys.stdout.flush()
def kw(**kwargs):
    """
    Method that returns the default parameter keywordarguments
    """
    KW=dict(
        config=os.path.join(os.path.dirname(os.path.abspath(sys.argv[0])),'run.conf')\
        )
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


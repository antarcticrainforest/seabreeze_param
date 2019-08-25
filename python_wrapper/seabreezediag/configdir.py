'''
Created on August 13, 2013

@author: Martin Bergemann

@institution: Monash University School of Mathematical

@description: Collect information from a given config file and return the
              content as a dict-object
'''
import os,sys,string,glob
from datetime import timedelta,datetime
from netCDF4 import Dataset as nc, date2num,num2date
import numpy as np
class Meta(object):
    '''
    This class reads important (static)-meta data from netcdf-files;
    '''
    def __init__(self,C):
        '''
        Variables :
            C (Config ojcet) : A Config oject containing the configuration
        
        Instances :
            lon = The longitude vector (numpy array)
            lat = The latitude vector (numpy array)
            start = the first date of the considered period (datetime object)
            end = the last date of the considered period (datetime object)
            dates = list of netcdf filenames with data between start and end (list object)
            datadir = the path to the data that should be read (str)
            prefix = the prefix of the netcdf-filenames containing the data (like ERAI_ , str)
            vp = the pressure vector, 1-D  (numpy array)
            vv = the v-wind field, 3-D (numpy array)
            vu = the u-wind field, 3-D (numpy array)
         Methods:
             read     : Read as variable from an opened netcdf file
             creat_nc : Output the seabreeze data into a netcdf-file
        '''

        #Open the meta-data of the file containing the distance form the 
        # coast information 
        for fn, vn in ((C.landfracfile, 'landfrac'),
                       (C.topofile,'z'),
                       (C.orofile,'std')):
            with nc(os.path.expanduser(fn)) as f:
                try:
                    setattr(self,vn,f.variables[C['v%s'%vn]][0,0,:,:])
                except (IndexError,ValueError):
                    try:
                        setattr(self,vn,f.variables[C['v%s'%vn]][0,:,:])
                    except (IndexError,ValueError):
                        setattr(self,vn,f.variables[C['v%s'%vn]][:])
                self.lon = f.variables[C.vlon][:]
                self.lat = f.variables[C.vlat][:]
        self.start=datetime.strptime(C.start,'%Y-%m-%d_%H:%M')
        self.end=datetime.strptime(C.end,'%Y-%m-%d_%H:%M')
        for i in ('vtheta', 'prefix', 'vpres', 'datadir', 'vu', 'vv'):
            if i == 'datadir':
                setattr(self, i, os.path.expanduser(C[i]))
            else:
                setattr(self, i, C[i])
        #Get all netcdf files between start and end
        self.__get_dates()

    def read(self, f, varname, timestep=None):
        '''
        This method reads the a netcdf variable from an opened netcdf file
        Variables:
            f (netcdf-ojbect) : object of the netcdf
            varname (str) : variable name to be read
            timestep (int) : the time step that should be read, if None 
                             the whole array is read (default: None)
        Returns:
            ND-array : data 
        '''
        if type(timestep) == type(None):
            data = f.variables[varname][:]
        else :
            data = f.variables[varname][timestep]
        return data

    def __get_dates(self):
        ''' Return all netcdf-files containing data for the considered period
        Variables : None
        Returns : None
        '''
        ts = self.start
        dt = timedelta(days=1)
        dates = []
        fn = os.path.join(self.datadir,'%04i'%ts.year,self.prefix+'*'+\
                    '%s'%ts.strftime('%Y_??_??.nc'))
        if len(glob.glob(fn))>0:
            daily=True
        else:
            fn = os.path.join(self.datadir,'%04i'%ts.year,self.prefix+'*'\
                    +'%s'%ts.strftime('%Y_??.nc'))
            if len(glob.glob(fn))>0:
                daily=False
            else:
                raise ValueError('Only daily or monthly file-format is supported\n')

        while ts < self.end:
            if daily:
                tstring = ts.strftime('%Y_%m_%d')
            else:
                tstring = ts.strftime('%Y_%m')
            fn = os.path.join(self.datadir,'%04i'%ts.year,self.prefix+'*_'\
                    +'%s.nc'%tstring)
            all_present = True
            for v in (self.vv, self.vu,self.vtheta):
                F = fn.replace('*',v)
                if not os.path.isfile(F):
                    all_present=False

            if all_present and tstring not in dates:
                dates.append(tstring)
            ts += dt

        self.dates = dates

    def create_nc(self,data,fname,varname,times,add=''):
        '''Output the seabreeze data into a netcdf-file
        Variables:
            data (ND-array) : The data that should be written to a netcdf file
            fname (str) : The netcdf-filename
            varname (str) : The name of the netcdf variable
            times (1D-array) : The time vector
            add : suffix for additional information to be added to the netcdf.long_name attribute 
                  (dfault : None)
        '''
        if os.path.isfile(fname):
            mode = 'a'
        else:
            mode = 'w'
        lookup = dict(\
                thc=dict(name='Thermal Heating Contrast Between Land and Ocean',units='K'),
                sb_con=dict(name='Subgrid Sea-Breeze Convergence',units=' '),
                windspeed=dict(name='Coastal Windspeed',units='m/s'),
                winddir=dict(name='Coastal Wind Direction',units='deg'),
                temp=dict(name='2M Temperture', units='degC'))


        with nc(fname,mode) as f:
            for i in ('lat','lon','time'):
                try:
                    if i == 'time':
                        ll = None
                        typ='i'
                    else:
                        ll=len(getattr(self,i))
                        typ='f'
                    f.createDimension(i,ll)
                    f.createVariable(i,typ,(i,))
                except (OSError,RuntimeError):
                    pass
            f.variables['lon'].units='degrees_east'
            f.variables['lat'].units='degrees_north'
            f.variables['lon'].axis='X'
            f.variables['lat'].axis='Y'
            f.variables['lat'].long_name='Latitude'
            f.variables['lon'].long_name='Longitude'
            f.variables['time'].units='Seconds since 1970-01-01 00:00:00'
            f.variables['time'].long_name='Time'
            f.variables['time'].axis='T'

            f.variables['lon'][:]=self.lon
            f.variables['lat'][:]=self.lat
            f.variables['time'][:]=date2num(times,f.variables['time'].units)

            try:
                f.createVariable(varname,'f',('time','lat','lon'))
            except (OSError,RuntimeError):
                pass


            f.variables[varname][:]=data
            f.variables[varname].long_name=lookup[varname]['name']+add
            f.variables[varname].units=lookup[varname]['units']
            f.variables[varname].grid = 'lonlat'
            f.variables[varname].missing_value=np.float32(2.0e20)



class Config(dict):
    """
    Config-class : Config(filename,**kwargs)
    Read a configuration file and return it's content stored in a dictionary 
    object.

    Parameters
    ----------
    filename : (str-obj) 
         the name of the configuration file to be read
    Keywords
    --------
    maketuple : (bool) Default: True
        if maketuple is True the method tries to interprete values with , 
        as seperaters for different tuple values
    skipwhitespace: (bool) Default: True
        whitspaces wont be considered if set True
    split : (str) Default =
        the seperator to seperate key and value in the configuration file

    Example
    -------
    Consider the following simple configuration file:
    #Filename of the test data
    filename = 'foo.nc' #
    variable = bar # The variable to be considered
     x1 = 9.0 # First index
    x2 =10  # Last index
    update = true
    times = 1,2,3 # Some time steps

    >>> C = Config('files.conf')
    >>> print(C)
    >>> Keys     | Values
    >>> -------------------------
    >>> update   | True
    >>> times    | (1.0, 2.0, 3.0)
    >>> x2       | 10
    >>> filename | foo.nc
    >>> variable | bar
    >>> x1       | 9.0

    """
    def __setattr__(self,k,v):
        if k in self.keys():
            self[k] = v
        elif not hasattr(self,k):
            self[k] = v
        else:
            raise AttributeError( "Cannot set '%s', class attribute already \
                    exists" % (k, ))
    def __getattr__(self, k):
        if k in self.keys():
            return self[k]
        raise AttributeError("Attribute '%s', deos not exist, available\
                attirbutes are: %s" %(k,self.keys().__repr__().strip(']')\
                .strip('[')))
    def __repr__(self):
        a=8
        b=8
        for v,k in self.items():
            if len(str(v)) > a:
                a = len(str(v))
            if len(str(k)) > b:
                b = len(str(k))
        a = a + 1
        b = b + 1
        s="Keys"+(a-4)*' '+"| Values\n"
        s=s+(a+b)*'-'+'\n'
        for v,k in self.items():
            s=s+str(v)+(a-len(str(v)))*' '+'| '+str(k)+'\n'
        return s
    def __init__(self,filename,maketuple=True,skipwhitespace=True,split='='):
        """ Try to open the file filename, read it and create a class instance 
            for every key entry 
            NOTE: Every key is an instance of Config but Config itself is of type dict
                  therefore the intances can be accesses both ways
        Example:
        -------
        from seabreezediag.configdir import Config
        C = Config('Path/to/run.conf')
        >>> T = C['times']
        >>> t = C.times
        >>> T == t
        >>> Ture


        """

        try:
            if isinstance(filename,str):
                f=open(filename)
            file = f.readlines()
            f.close()
        except IOError:
            print( "Could not open %s: No such file or directory" %filename)
            return
        self.__get(file,maketuple,skipwhitespace,split=split)
        for key,value in self.items():
            try:
                if "$" in value:
                    var=value.split('/')[0]
                    try:
                        path=os.environ[var.replace('$','')]
                        self[key]=value.replace(var,path)
                    except KeyError:
                        raise KeyError('Environmet variable %s not set'%var)
            except TypeError:
                pass
                


    def __get(self,file,maketuple,skipwhitespace,split='='):
        """Explanation: This function takes a variable-name and extracts
        the according value form the config-file
        varname : the name of the variable that should be looked up"""
        
        #Define a blacklist of characters
        blacklist=[']','[','{','}','@','#','"',"'"]
        strip={'\t':'=','=':'\t'}[split]
        for i,line in enumerate(file):
            try:
                if not line.strip('\n')[0] in blacklist and split in line:
                    var=line.strip('\n').strip(strip).strip()\
                            .split(split)
                    if skipwhitespace :
                        value=var[1].replace(' ','')
                    else:
                        value=var[1]
                    for pos,v in enumerate(value):
                        if v.startswith('#'):
                            value=value[:pos+1]
                            break

                    for b in blacklist:
                        value=value.replace(b,'')
                    try:
                        value = int(value)
                    except ValueError:
                        try :
                            value = float(value)
                        except ValueError:
                            if value.lower() == 'false':
                                value = False
                            elif value.lower() == 'true':
                                value = True
                            elif value.lower() == 'none':
                                value = None
                    if isinstance(value,str):
                            #value=value.strip()
                        if ',' in value and maketuple:
                            value2=[]
                            for v in value.split(','):
                                try:
                                    value2.append(float(v.strip(')').strip('(')))
                                except ValueError:
                                    value2.append(v.strip(')').strip('('))
                            value=tuple(value2)

                    self[var[0].replace(' ','').strip()]=value
            except IndexError:
                pass

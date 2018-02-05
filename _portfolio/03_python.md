---
title: Python Wrapper
permalink: /python
---
This documentation describes the python wrapper to test the parametrization and 
to create some offline climatologies. The code is located in the 
[python_wrapper](https://github.com/antarcticrainforest/seabreeze_param/python_wrapper)
sub folder. 
## Prerequisites
A working python 2 or python 3 version is required to run the code. You'll also 
need a working [numpy](http://www.numpy.org) python package. It is also recommended 
to have the [netCDF4](https://pypi.python.org/pypi/netCDF4) package available to your 
python distribution. You will also need the gfrotran and gcc compilers in order 
to build the python module from the fortran source files. The compiling 
process involves the fortran to C wrapper [f2py](http://www.f2py.com) which comes 
with the numpy's distutils package.

## Building and Installing
Building and installing is most simple:
```bash
$: python setup.py build
$: python setup.py intall
```
The latter command will install the built python package into the root file system; 
usually ```/usr/lib/pythonX.X/site-packages```. If you don't have write acces to 
```/usr/lib/``` install the packages with the ```--user``` option:
```bash
$: python setup.py install --user
```
This will install the packages in ```$HOME/.local/lib/pythonX.X/site-packages/```. 
With X.X being the python version (e.g 3.6). If this directory doesn't already exsist create it with:
```bash
$: mkdir -p $HOME/.local/lib/pythonX.X/site-packages/
```

You should now be able to import the seabreeze package:
```python
>>> import seabreezediag as sbd
```
## The Classes
The python package contains a couple of classes which are described here.
### The ```Config``` class
Read a configuration file and return it's content stored in a dictionary 
object.
<ul style="list-style-type:none"><li>
Variables:
<pre><code class="language-b"><font size='2'><b>filename</b> (str) : name of the config file to be read
<b></b> (2D-array) : Land-sea mask of the model (lat,lon)
</font></code></pre></li></ul>
<ul style="list-style-type:none"><li>
Keywords:
<pre><code class="language-b"><font size='2'><b>maketuple</b> (bool)    : if maketuple is True the method tries 
                         to interprete values with , as seperators for tuple values 
                         (<em>defautl</em> : True)
<b>skipwhitespace</b> (bool) : whitespaces wont be considered if set True
                         (<em>defautl</em> : True)
<b>split</b> (str) : the seperator to seperate key and value in the config file
                      (<em>defautl</em> : =)
</font></code></pre></li></ul>

Example:

Consider the following simple configuration file:
```bash
$: cat test.conf
#Filename of the test data
filename = 'foo.nc' #
variable = bar # The variable to be considered
 x1 = 9.0 # First index
x2 =10  # Last index
update = true
times = 1,2,3 # Some time steps
```

```python
>>> from seabreezediag.configdir import Config
>>> C = Config('test.conf')
print(C)
Keys     | Values
-------------------------
update   | True
times    | (1.0, 2.0, 3.0)
x2       | 10
filename | foo.nc
variable | bar
x1       | 9.0
```

The class tries to open the file filename, read it and create a class instance 
for every key entry

**Note**: Every key is an instance of Conf but Conf itself is of type dict
          therefore the intances can be accesses both ways 
{: .notice--info}
Example:
```python
>>> T = C['times']
>>> t = C.times
>>> T == t
Ture
```


### The ```Meta``` class
The ```Meta``` class is not a meta class. It can be used to read important (static) 
meta data from netcdf-files; like Longitude, Latitude vectors and other data.

The object is created by calling :
```python
from seabreezediag.configdir import Meta,Config
C = Config('path/to/config.conf')
M = Meta(C)
```
The only variable that ```Meta``` is created with has to be of type ```Config```.

```Meta``` has the following *instances*:

<ul style="list-style-type:none"><li>
<pre><code class="language-b"><font size='2'><b>lon</b> (1D-array)   : the longitude vector
<b>lat</b> (1D-array)   : the latitude vector
<b>start</b> (datetime) : the first date of the considered period
<b>end</b> (datetime)   : the last date of the considered period
<b>dates</b> (list)     : list of netcdf filenames with data between start and end
<b>datadir</b> (str)    : the path to the data that should be read
<b>prefix</b> (str)     : the prefix of the netcdf-filenames containing the data (like ERAI)
<b>vp</b> (1D-array)    : the pressure vector
<b>vv</b> (ND-array)    : the v-wind field
<b>vu</b> (ND-array)    : the u-wind field
<b>vtheta</b> (ND-array): the surface temperature field
</font></code></pre></li></ul>
The following *methods* are available:

#### ```Meta.create_nc(fname,varname,times,add=None)```:

<ul style="list-style-type:none"><li>
Output the seabreeze data into a netcdf-file
<ul style="list-style-type:none"><li>
Variables:
<pre><code class="language-b"><font size='2'><b>data</b> (ND-array)  : the data that should be written to a netcdf file
<b>fname</b> (str)      : the netcdf-filename
<b>varname</b> (str)    : the name of the netcdf variable
<b>times</b> (1D-array) : the time vector
<b>add</b>              : suffix for additional information to be added to the netcdf.long_name attribute (dfault : None)
</font></code></pre></li></ul>
</li></ul>
#### ```Meta.read(varname,timestep=None):```
<ul style="list-style-type:none"><li>
This method reads a variable from an opened netcdf file
<ul style="list-style-type:none"><li>
Variables:
<pre><code class="language-b"><font size='2'><b>f</b> (netcdf-ojbect)  : object of the netcdf
<b>varname</b> (str)      : variable name to be read
<b>timestep</b> (int)     : the time step that should be read, if None the whole array is read (default: None)
</font></code></pre></li></ul>
<ul style="list-style-type:none"><li>
Returns:
<pre><code class="language-b"><font size='2'><b>ND-array</b></font></code></pre>
</li></ul>
</li></ul>
<ul style="list-style-type:none">
<li><b>Note</b> For Meta to work Config needs the following entries:</li>
<li><b>landfracfile</b> : Filename of land area fraction (land-sea mask) data</li>
<li><b> topofile</b>     : Filename of orography data</li>
<li><b>orofile</b>      : Filename of std of sub-grid orography</li>
<li><b>vlon</b>         : Variable name for longitude vector</li>
<li><b>vlat</b>         : Variable name for latitude vector</li>
<li><b>start</b>        : Start date of the considered period</li>
<li><b>last</b>         : End date of the considered period</li>
<li><b>datadir</b>      : Parent directory of the netcdf data</li>
<li><b>prefix</b>       : Prefix of the netcdf-file names, like ERAI</li>
<li><b>vp</b>           : Name of the pressure variable</li>
<li><b>vu</b>           : Name of the u-wind variable</li>
<li><b>vv</b>           : Name of the v-wind variable</li>
<li><b>vtheta</b>       : Name of the surf. temp. variable</li></ul>
{: .notice--info}
<ul class="notice--info">
<p><b>Note</b>: The structure of the netcdf files containing the data to be read 
should have the following format:</p>
<pre><code class="language-b">datadir/YYYY/prefix_YYYY_MM_DD.nc</code></pre>
<p>for daily data or:</p>
<pre><code class="language-b">datadir/YYYY/prefix_YYYY_MM.nc</code></pre>
<p>for monthly data.</p></ul>

### The ```seabreezediag``` class
This is the class that does the actual work. It is simply imported by:
```python
>>> import seabreezediag as sbd
```
The class offers the following methods:
#### ```sbd.f2c(array)```
<ul style="list-style-type:none"><li>
Convert from column/row major to row/column major. This function has been added, 
to convert 
<ul style="list-style-type:none"><li>
Variables:
<pre><code class="language-b"><font size='2'><b>f</b> (netcdf-ojbect)  : object of the netcdf
<b>array</b> (ND-array)  : array of rank N and any typ to be converted
</font></code></pre></li></ul>
<ul style="list-style-type:none"><li>
Returns:
<pre><code class="language-b">
<font size='2'><b>ND-array</b>         : Converted array</font></code></pre></li></ul>
</li></ul>

#### ```sbd.read_nc(fnv, fnu, fntheta, fnci, vv='v', vu='u', vtheta='t2m', vci='ci')```
<ul style="list-style-type:none"><li>
Method to read relevant data from various netcdf-sources
<ul style="list-style-type:none"><li>
Variables:
<pre><code class="language-b"><font size='2'><b>f</b> (netcdf-ojbect)  : object of the netcdf
<b>fnv</b> (str)     : name of the v-wind netcdf file
<b>fnu</b> (str)     : name of the u-wind netcdf file
<b>fntheta</b> (str) : name of the temp netcdf file
<b>fnci</b> (str)    : name of the sea-ice frac. netcdf file
<b>vv</b> (str)      : name of the v-wind netcdf variable (default : v)
<b>vu</b> (str)      : name of the u-wind netcdf variable (default : u)
<b>vtheta</b> (str)  : name of the of the surf temp var. (default : t2m)
<b>vci</b> (str)     : name of the sea-ice variable (default : ci)
</font></code></pre></li></ul>
<ul style="list-style-type:none"><li>
Returns:
<pre><code class="language-b"><font size='2'><b>namedtuple</b>         : Named tuple containing all information</font></code></pre>
<p>Example:</p>
<div class="language-python highlighter-rouge"><div class="highlight"><pre class="highlight"><code><span class="o">&gt;&gt;&gt;</span> <span class="n">e</span> <span class="o">=</span> <span class="s">'~/Data/ERAI/netcdf/1988/Erai_'</span>
<span class="o">&gt;&gt;&gt;</span><span> nc_data </span> <span class="o">=</span> <span class="n">sbd</span><span class="o">.</span><span class="n">read_nc</span><span class="p">(</span><span class="n">e</span><span class="o">+</span><span class="s">'v_1988_09.nc'</span><span class="p">,</span><span class="n">e</span><span class="o">+</span><span class="s">'u_1988_09.nc'</span><span class="p">,</span>
        <span class="n">e</span><span class="o">+</span><span class="s">'t2m_1988_09.nc'</span><span class="p">,</span> <span class="n">e</span><span class="o">+</span><span class="s">'ci_1988_09.nc'</span><span class="p">)</span>
</code></pre></div></div>
<p>The returned <em>namedtuple</em> has the following <em>field names</em>:</p>
<div class="language-python highlighter-rouge"><div class="highlight"><pre class="highlight"><code><span class="o">&gt;&gt;&gt;</span> <span class="k">print</span><span class="p">(</span><span class="n">nc_data</span><span class="o">.</span><span class="n">_fields</span><span class="p">)</span>
<span class="o">&gt;&gt;&gt;</span> <span class="p">(</span><span class="s">'time'</span><span class="p">,</span> <span class="s">'pres'</span><span class="p">,</span> <span class="s">'dt'</span><span class="p">,</span> <span class="s">'nc'</span><span class="p">,</span> <span class="s">'v'</span><span class="p">,</span> <span class="s">'u'</span><span class="p">,</span> <span class="s">'theta'</span><span class="p">,</span> <span class="s">'ci'</span><span class="p">)</span>
</code></pre></div></div>

<pre><code class="language-b"><font size='2'>
<b>time</b>  : time vector (1d-array)
<b>pres</b>  : pressure variable object (netcdf-object)
<b>dt</b>    : model timestep (int)
<b>nc</b>    : dictionary with all opened netcdf-file objects (dict)
<b>v</b>     : v-wind variable (netcdf-object)
<b>u</b>     : u-wind variable object (netcdf-object)
<b>theta</b> : surf-temp variable object (netcdf-objcet)
<b>ci</b>    : Sea-ice fract netcdf-variable object (netcdf-object)
</font></code></pre>
</li></ul>
</li></ul>
<br>
#### ```sbd.diag(tt, lsm, z, std, lon, lat, *args, **kwargs)```
<ul style="list-style-type:none"><li>
Method to calculate potential strength of sea-breeze convergence within a predefined
coastal area.
<ul style="list-style-type:none"><li>
Variables:
<pre><code class="language-b"><font size='2'><b>tt</b> (int)       : Number of timestep since start of application,
                 used to determine the beginning of simulation and and the
                 application of the algorithm in the considered interval
<b>lsm</b> (2D-array) : Land-sea mask of the model (lat,lon)
<b>z</b> (2D-array)   : Surface elevation data (lat,lon)
<b>std</b> (2D-array) : Standard deviation of subgrid orography
<b>lon</b> (1D-array) : The longitude vector of the model
<b>lat</b> (1D-array) : The latitude vector of the model
<b>p</b> (1D-array)   : Pressure levels in Pa stored in a 1D array
<b>u</b> (ND-array)   : U-component of the wind stored ([time],pres,lat,lon) ,time is optional
<b>v</b> (ND-array)   : V-component of the wind stored ([time],pres,lat,lon) ,time is optional
<b>t</b> (ND-array)   : Surface temperature array  ([time],lat,lon), time is optional
<b>ci</b> (ND-array)  : Fraction of sea-ice cover [0..1] ([time],lat,lon)
                 time is optional. Ci can be set to None, in which case 
                 it won't be taken into account for calculation.
</font></code></pre></li></ul>
<ul style="list-style-type:none"><li>
Keywords:
<pre><code class="language-b" ><font size='2'><b>windspeeed</b> (2D-array)   : Wind speed, passed through the application period (lat,lon)
                          (<em>default</em>: zero array)
<b>winddir</b> (2D-array)      : Wind direction, passed through the application period (lat,lon)
                          (<em>default</em>: zero array)
<b>thc</b> (2D-array)          : Thermal heating contrast, passed through the application period (lat,lon)
                          (<em>default</em>: zero array)
<b>target_plev</b> (float)     : P-level (hPa) for wind thresh application 
                          (<em>default</em>: 700.)
<b>thresh_wind</b> (float)     : Thresh. of wind speed  (m/s)
                          (<em>default</em>: 11.)
<b>thresh_winddir</b> (float)  : Thresh. of wind speed  change (m/s)
                          (<em>default</em>: 11.)
<b>thresh_winddir</b> (float)  : Thresh. of wind direction  change (deg.)
                          (<em>default</em>: 90.)
<b>thresh_windch</b> (float)   : Thresh. of wind speed change (m/s)
                          (<em>default</em>: 5.)
<b>thresh_thc</b> (float)      : Thresh. of thermal heating contrast (K)
                          (<em>default</em>: 0.75)
<b>target_time</b> (float)     : Time perid for application of thresholds (h)
                       (<em>default</em>: 6.)
<b>maxdist</b> (float)         : Area of influence by the sea-breeze on-/offshore (km)
                       (<em>default</em>: 180.)
<b>timestep</b> (float)      : Time stepping of the model (mins)
                       (<em>default</em>: 24.)
<b>meta</b> (float)          : If meta is of type <em>namedtuple</em> the arguments 
                       for the u-array, v-array, t-array and ci-array are skipped and
                       the data is taken from meta (<em>default</em>: None)
</font></code></pre></li></ul>


<ul style="list-style-type:none"><li>
Returns:
<pre><code class="language-b"><font size='2'><b>timestep</b> (int)   : timestep after the application
<b>breeze</b> (ND-array): sub-grid scale see-breeze convergence ([time],lat,lon)
<b>thc</b> (2D-array)   : thermal heating contrast of last timestep as initial cond. for next timestep
<b>ws</b> (2D-array)    : windspeed of last timestep as init. cond. for next timestep
<b>wd</b> (2D-array)    : wind direction of last timestep as init. for next timestep
</font></code></pre>
<p>Example:</p>
</li></ul>
</li></ul>
```python
>>> from netCDF4 import Dataset as nc, num2date
>>> import numpy as np
>>> from datetime import datetime, timedelta
>>> import seabreezediag as sbd
>>> tt = 1 #Set the first timestep to one
>>> THC = np.zeros([T.shape[-3],T.shape[-2]])
>>> WS,WD = np.zeros_like(thc),np.zeros_like(thc)
>>> for fn in  ('input_netcdf-month1.nc','input_netcdf-month2.nc'):
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
        >>> tt,breeze,THC,WS,WD = sbd.diag(tt,LSM,H,CI,P,U,V,T,CI,WS,WD,
        ...                                    THC,timestep=dt)
        >>> test_nc.variables['breeze'][:]=breeze
        >>> test_nc.close()
```


The file ```test_run.py``` in the  ```python_wrapper``` directory provides an 
example application of the module.


If you are considering implementing this routine take a look at the  [about](/zz_about) 
section to learn more on how to contribute and improve this project. You are encouraged 
to get in touch via [GitHub](https://github.com/antarcticrainforest/seabreeze_param). 
Bugs should be reported either on the GitHub [issues](https://github.com/antarcticrainforest/seabreeze_param/issues) 
pages or by sending an [email](mailto:martin.bergemann@monash.edu) to the author of this page.

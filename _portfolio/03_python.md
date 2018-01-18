---
title: Python Wrapper
permalink: /03_python
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
to compile to bake the python module from the fortran source files. The compiling 
process involves a fortran to C wrapper [f2py](http://www.f2py.com) which comes 
with the numpy's distutils package.

## Building and Installing
Building and installing is most simple:
```bash
$: python setup.py build
$: python setup.py intall
```
Building 
The last command will install the builded python package into the root file system; 
usually ```/usr/lib/pythonX.X/site-packages```. If you don't have write acces to 
```/usr/lib/``` install the packages with the ```--user``` option:
```bash
$: python setup.py install --user
```
This will install the packages in ```$HOME/.local/lib/python3.6/site-packages/```. 
If this directory doesn't already exsist create it with:
```bash
$: mkdir -p $HOME/.local/lib/pythonX.X/site-packages/
```

You should now be able to import the seabreeze package:
```python
>>> import seabreezediag as sbd
```
## Usage
The python package contains a couple of classes which are described here.
### The ```Config``` class
Read a configuration file and return it's content stored in a dictionary 
object.

Parameters:

filename : (str-obj) 
     the name of the configuration file to be read
Keywords:

maketuple : (bool) Default: True
    if maketuple is True the method tries to interprete values with , 
    as seperaters for different tuple values
skipwhitespace: (bool) Default: True
    whitspaces wont be considered if set True
split : (str) Default =
    the seperator to seperate key and value in the configuration file

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

```Meta``` has the following instances:

* lon (1D-array)
  * The longitude vector
* lat (1D-array)
  * The latitude vector
* start (datetime)
  * the first date of the considered period
* end (datetime)
  * the last date of the considered period
* dates (list)
  * list of netcdf filenames with data between start and end
* datadir (str)
  * the path to the data that should be read
* prefix (str)
  * the prefix of the netcdf-filenames containing the data (like ERAI)
* vp (1D-array)
  * the pressure vector
* vv (ND-array)
  * the v-wind field
* vu (ND-array)
  * the u-wind field
* vtheta (ND-array)
  * the surface temperature field

The following methods are available:
#### ```create_nc```:
Output the seabreeze data into a netcdf-file

Variables:

* data (ND-array) 
  * The data that should be written to a netcdf file
* fname (str) 
  * The netcdf-filename
* varname (str) 
  * The name of the netcdf variable
* times (1D-array) : 
  * The time vector
* add 
  * suffix for additional information to be added to the netcdf.long_name attribute
              (dfault : None)

#### ```read```
This method reads a variable from an opened netcdf file

Variables:

* f (netcdf-ojbect) 
   * object of the netcdf
* varname (str) 
   * variable name to be read
* timestep (int) 
   * the time step that should be read, if None the whole array is read (default: None)

Returns:

*   ND-array

<ul style="list-style-type:none">
<li><b>Note</b>: For Meta to work Config needs the following entries:</li>
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

If you are considering implementing this routine take a look at the  [about](/zz_about) 
section to learn more on how to contribute and improve this project. You are encouraged 
to get in touch via [GitHub](https://github.com/antarcticrainforest/seabreeze_param). 
Bugs should be reported either on the GitHub [issues](https://github.com/antarcticrainforest/seabreeze_param/issues) 
pages or by sending an [email](mailto:martin.bergemann@monash.edu) to the author of this page.

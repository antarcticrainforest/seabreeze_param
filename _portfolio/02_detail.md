---
title: Code Structure
permalink: /generic
---
This page describes the generic structure of the code. It is not model specific 
and servers only as an outline to understand the purpose of the routines. You will 
not be able to apply this source (everything in the ```model/generic``` folder) 
out of the box. The aim of this page is to explain how the source code is organized 
and how the routines can be called. The preparation varies from model to model, 
therfore it is essantail that you pay attention to the **Note** boxes.

<h2 id='the-sea_breeze_diag-module'>The <code class='highlighter-rouge'><font color='#6E6E6E'>sea_breeze_diag</font></code> module</h2>
The code for the sea-breeze parametrization is located in ```sea_breeze_diag.F90```
and contains a module with all subroutines that are need to apply the [previously](/) 
introduced algorithm. This section should explain how the module can be loaded
and how its contained subroutines should be called.

If you're using a model other than the UK MetOffice's UM make sure that the source 
file ```sea_breeze_diag.F90``` is located in a meaningful subdirectory of the model 
and that code gets compiled and linked. If you feel adding a new docu on how this can 
be done for your model please refer to the [about](/zz_about) section on this site.

Once correctly compiled and linked the sea-breeze module can be loaded by:
```fortran
use sea_breeze_diag_mod, only: &
  sea_breeze_diag, get_edges, get_dist
```
although module contains more the subroutines only ```sea_breeze_diag```, ```get_dist```,
and ```get_edges``` have to be imported to apply the parametrization.

## Subroutines
The ```sea_breeze_diag_mod``` contains the following subroutines:
<li> <code class='highlighter-rouge'><font color='#6E6E6E'>sea_breeze_diag</font></code>
  <ul style="list-style-type:none">
  <li><font size='3'>Calculates the sub-grid scale sea-breeze strength (see 
  <a href='/#the-algorithm'>The Algorithm</a>)</font></li>
  </ul>
</li>
<li> <code class='highlighter-rouge'><font color='#6E6E6E'>get_edges</font></code>
 <ul style="list-style-type:none">
  <li><font size='3'>apply sobel edge detection to calculate the area of influence 
  (see <a href='/#filtering-of-coastal-areas'>coastal filtering</a>)</font></li>
 </ul>
</li>
<li> <code class='highlighter-rouge'>get_dist</code>
 <ul style="list-style-type:none">
  <li><font size='3'>calculates the euclidian distance to the closest points, 
  that where deemed to be coastal output of <code class='highlighter-rouge'>get_edges</code>)</font></li>
 </ul>
</li>
<li> <code class='highlighter-rouge'>sigmoid</code>
 <ul style="list-style-type:none">
  <li><font size='3'>only called by <code class='highlighter-rouge'>sea_breeze_diag</code>
  and applies a sigmoid function to the sub-grid orography field in 
  <a href='/#steep-terrain'>steep terrain</a>.</font></li>
 </ul>
</li>
The following sub-section provide more details on how to call the modules 
subroutines.

### ```get_edges``` and ```get_dist```
This subroutine creates a coastal mask from the land-area and the 
sea-ice fraction fields. Alternativly the land-sea mask instead of land-area 
fraction can be use. Both arrays (land-area fraction/land-sea mask and sea-ice 
fraction) are excpected to be of 2D shape and of type real. After the coastal mask 
has been created the ```get_dist``` routine has to be called to create a field 
with the euclidian distance to the next coast point 
(defining the area of influence by the sea-breeze circulation). 
Bcause ```get_dist``` makes use of information from neigboring grid points 
a *boundary swapping* of the coastal mask has to be performed between the calls. 
This is a model specific task (but yet important), therefore ```get_dist``` *can't* 
be called from within ```get_edges```.
#### Variables:
```get_edges``` requires following variables:
* ```mask, intent(in)```: 
   
   real valued array of rank 2. This can be either the sum of 
the land-area fraction and sea-ice fraction fields, or the land-sea mask. A sum 
of the land-sea mask and the sea-ice fraction is also possible. 
* ```coast, intent(out)```: 
   
   the routine returns real valued array of rank 2. The 
field contains a mask of the coastline.
```get_dist``` requires :
* ```coast, intent(in)``` : 
   
   the coastline mask (returned by ```get_edges```)
* ```mask, intent(in)``` : 
   the same mask that was passed into ```get_edges``` 
(e.g. land_area_frac + sea_ice_frac)
* ```dist, intent(out)``` : 
   
   euclidian distance to the closest coastal point, in ```mask```

The above variables are mandatory. Other variables might be subject to the models
implementations and infrastructure.

The following code snipset serves as an minimal example of how to call the routines:
```fortran

subroutine just_a_test( landfrac,icefrac,dist,lon, &
                        lat,n,nlons,nlats )
  use sea_breeze_diag_mod, only:                   &
    sea_breeze_diag, get_edges, get_dist
  use a_swapping_boundaries_mod,                   &
    only swap ! made up - example only ;)

  implicit none

  integer, intent(in) :: nlons, nlats,n
  real, intent(in), dimension(nlons,nlats)::       &
    landfrac, icefrac !from somewhere
  real, intent(in), dimension(nlons) ::            &
    lon ! lon vector
  real, intent(in), dimension(nlats) ::            &
    lat ! lat vector
  real, intent(out), dimension(nlons+n,nlats+n) :: &
    dist ! coastal distance
  real, dimension( nlons+n,nlats+n ) ::            &
    coast, & ! The coastal points mask
    mask
  
  !Crate the mask that is passed into get_edges
  mask( 1+n/2:nlons+n/2,1+n/2:nlats ) = landfrac + icefrac
  
  !Make a perform a boundary swapping (model specific)
  swap(mask,n,nlons,nlats)
  get_edges(mask,coast,n,nlons,nlats)
  !Swap the boundaries for the coastal mask
  swap(coast,n,nlons,nlats)
  get_edges(coast,mask,lon,lat,dist,n,nlons,nlats)
end subroutine just_a_test
```

**Note:** The ```a_swapping_boundaries_mod``` module and the ```swap``` have to be
replaced by the appropriate modules and calls according to the models boundary 
swapping routine.
{: .notice--info}


**Caution:** Since boundary swapping is necessary to apply the above routines it is 
strongly suggested that they are called *before* any atmospheric physics is called.
Calling the routine from within any atmospheric physics routine will (depending 
on the implementation of the model) result in a communication overflow and hence 
reduce the performence of the model significantly.
{: .notice--warning} 



### ```sea_breeze_diag```

```sea_breeze_diag``` is the routine that calcultates the actual sub-grid sea-breeze 
strength. This routine can be called from within any of the atmospheric physics 
routines, but perferably *before* the cumulus scheme is involved. 

**Caution:** The routine makes use of the surface temperature field. If this 
field hasn't been subject to boundary swapping it is adviced to do boundary_swapping 
of the surface temp field *before* ```sea_breeze_diag``` is called.
{: .notice--warning}


#### Variables:
```sea_breeze_diag``` requires following variables:
* ```timestep, intent(in), float```:

   the timestep of the model (in seconds)
* ```timestep_number, intent(in), integer```: 

   scalar for the timestep number since start of the integration (N+1 per model timestep)

* ```p, intent(in), float```:

   array of rank 3 that contains the pressure field on model levels.
* ```u, intent(in), float```:

   array of rank 3 that contains the u component of the wind field on model levels.
* ```v, intent(in), float```:

   array of rank 3 that contains the v component of the wind field on model levels.

* ```theta, intent(in), float```:

   array of rank 2 that cointains the surface temperature field.
* ```cdist, intent(in), float```:

   array of rank 2 containing the euclidian distance to the closest coastal point
* ```windspeed, intent(in), float```:

  array of rank 2 containing the wind speed at 600 hPa (can be initialized with 0 at 
  the beginning of the modeling period)
* ```winddir, intent(in), float```:

  array of rank 2 containing the wind direction at 600 hPa (can be initialized with 0 at 
  the beginning of the modeling period)
* ```thc, intent(in), float```:

  array of rank 2 containing the thermal heating contrast (can be initialized with 0 at 
  the beginning of the modeling period)
* ```sb_con, intent(out), float```:
  array of rank 2. The output to the routine.


**Note:** the input arrays ```windspeed```, ```winddir``` and ```thc``` have to be 
primary variables. That means that they should not go out of scope when the physics routine 
that that calls ```sea_breeze_diag``` goes out of scope. The simplist way of preserving 
the content of the array is involving a ```save``` statement. But many other ways 
are possible, like pointers, depending on the models ifrastructure. The above variables 
are mandatory. Other variables might be subject to the models implementations 
and infrastructure. 
{: .notice--info}

The following code snipset should serve as a minimal example of how to call the routine:
```fortran
module some_example
  !This routine should initialize some fields
  implicit none

  real, allocatable, save :: & 
    windspeed(:,:), winddir(:,:), thc(:,:)
  contains
  subroutine init(nlons,nlats,data)
    implicit none
    integer, intent(in) :: nlons,nlats
    real, allocatable, save, intent(in) :: &
      data(:,:)
    if ( .not. allocated(data) ) then
      allocate(data(nlons,nlats))
      data(:,:) = 0.0
    endif
  end subroutine
end module some_example_init

  subroutine yet_another_test(timestep,timestep_number,&
     theta, p,u,v,sb_con,nlons,nlats,nk)
    use some_example
    use sea_breeze_diag_mod, only: sea_breeze_diag

    integer, intent(in) :: & 
      timestep_number,nlons,nlats,nk
    real, intent(in),dimension(nlons,nlats) :: &
      theta
    real, intent(in),dimension(nlons,nlats,nk) :: & 
      p , u , v
    real, intent(out),dimension(nlons,nlats) :: &
      sb_con

    call init(windspeed)
    call init(wind_dir)
    call init(thc)

    call sea_breeze_diag(timestep,timestep_number, theta, &
      p,u,v,windspeed,winddir,thc,sb_con,nlon,nlats,nk)


  end subroutine yet_another_test
```

If you are considering implementing this routine take a look at the  [about](/zz_about) 
section to learn more on how to contribute and improve this project. You are encouraged 
to get in touch via [GitHub](https://github.com/antarcticrainforest/seabreeze_param). 
Bugs should be reported either on the GitHub [issues](https://github.com/antarcticrainforest/seabreeze_param/issues) 
pages or by sending an [email](mailto:martin.bergemann@monash.edu) to the author of this page.

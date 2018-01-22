---
title: UM Setup
permalink: /04_install_um
---
This documentation should explain how to include the sea-breeze diagnosis into 
the UK Met Office Unified Model (UM). 

## The easy way
If you want to start with a *fresh* copy of the source code and the configuration 
it is recommendet to make a copy of the source code and the rose configuration.
### Copy the config
To copy a working rose configuration simply type 
```bash
$: mosrs-auth
$: rosie copy u-an760/trunk
```
The ```mosrs-auth``` command is used to sign in to the met office repository server. The 
```rosie``` command creates a new copy of the config in ```~/rose/```.
### Create a new branch of the source
Create a copy of the branch with the sea-breeze parametrization:
```bash
$: fcm branch-create --branch-of-branch vn10.7_anybranchname fcm:um.x/branches/dev/martinbergemann/vn10.7_SeaBreeze
```
this will create a copy of the source code in 
```config
 fcm:um.x/branches/dev/yourusername/revision#_vn10.7_anybranname
```
Now you can checkout the source for the freshly created branch:
```bash
$: fcm checkout fcm:um.x/branches/dev/yourusername/revision#_vn10.7_anybranchname vn10.7_anybranchname
```
You can find the revision number with :
```bash
$: cd vn10.7_anybranchname
$: fcm bls
fcm:um.x/branches/dev/yourusername/r48828_vn10.7_test@49179
```
The number after the ```@``` sign is the assigned revision number, in this case 49179.
Now navigate to the rose config directory and edit ```um_sources``` entry in ```/app/fcm_make_um/rose-app.conf```
```config
um_sources=branches/dev/yourusername/r48828_vn10.7_test@49179
```

The revision number has to be updated anytime the source code is changed and the changes have 
been commited.

## The manual way
### Adding variables to the stash
The sea breeze parametrization needs following new variables:

* 2D array of land area fraction
* thermal heating contrast between land and ocean (2D array)
* windspeed (2D array)
* wind direction (2D array)
* surface temperature, with halo (2D array)
* coastal mask, with extendet halo (2D array)
* sub grid sea-breeze strength (2D output array)

All variables except for the output array have to be *primary* variables. The 
land area fraction array is needed because the models original land area fraction 
field is compressed onto land points only. This field is the only one that needs 
an ancillary file for initialisation. All other fields can be initialized with one. 
To add the fields copy the content in the ```top_of_wokring_dir/rose-meta/um-atmos/HEAD/etc/stash/STASHmaster/```
directory to ```~/roses/roesID/app/um-app/file/```.

If the latter directory doesn't excist, create it. 
Copy the following lines into the STASHmaster_A file in the ```~/roses/roesID/app/um/file/``` directory:
```config
1|    1 |    0 |  301 |2D array for land area fraction     |
2|    2 |    0 |    1 |    1 |    5 |   -1 |   -1 |    0 |    0 |    0 |    0 |
3| 000000000000000000000000000000 | 00000000000000000001 |    3 |
4|    1 |    2 |  -3  -10   -3   -3  -12   21   -3  -99  -99  -99 |
5|    0 |  395 |    0 |  129 |    0 |    0 |    0 | 9999 |   74 |
#
1|    1 |    0 |  594 |Thermal heating contrast land-ocean |
2|    2 |    0 |    1 |    1 |    5 |   -1 |   -1 |    0 |    0 |    0 |    0 |
3| 000000000000000000000000000000 | 00000000000000000001 |    3 |
4|    1 |    2 | -99  -99  -99  -99  -99  -99  -99  -99  -99  -99 |
5|    0 |   16 |    0 |    0 |    0 |    0 |    0 |    0 |    0 |
#
1|    1 |    0 |  595 |Windspeed at sea-breeze routine call|
2|    2 |    0 |    1 |    1 |    5 |   -1 |   -1 |    0 |    0 |    0 |    0 |
3| 000000000000000000000000000000 | 00000000000000000001 |    3 |
4|    1 |    2 | -99  -99  -99  -99  -99  -99  -99  -99  -99  -99 |
5|    0 |   56 |    0 |    0 |    0 |    0 |    0 |    0 |    0 |
#
1|    1 |    0 |  596 |Wind direction at sb routine call   |
2|    2 |    0 |    1 |    1 |    5 |   -1 |   -1 |    0 |    0 |    0 |    0 |
3| 000000000000000000000000000000 | 00000000000000000001 |    3 |
4|    1 |    2 | -99  -99  -99  -99  -99  -99  -99  -99  -99  -99 |
5|    0 |   56 |    0 |    0 |    0 |    0 |    0 |    0 |    0 |
#
1|    1 |    0 |  597 |Surface temperature with halo       |
2|    2 |    0 |    1 |    1 |    5 |   -1 |   -1 |    0 |    0 |    0 |    0 |
3| 000000000000000000000000000000 | 00000000000000000001 |    1 |
4|    1 |    2 | -99  -99  -99  -99  -99  -99  -99  -99  -99  -99 |
5|    0 |   16 |    0 |    0 |    0 |    0 |    0 |    0 |    0 |
#
1|    1 |    0 |  598 |Coastal mask with extended halo     |
2|    2 |    0 |    1 |    1 |    5 |   -1 |   -1 |    0 |    0 |    0 |    0 |
3| 000000000000000000000000000000 | 00000000000000000001 |    2 |
4|    1 |    2 | -99  -99  -99  -99  -99  -99  -99  -99  -99  -99 |
5|    0 |   16 |    0 |    0 |    0 |    0 |    0 |    0 |    0 |
#
1|    1 |    5 |  594 |Subgrid Sea-Breeze Strength         |
2|    0 |    0 |    1 |    1 |    5 |   -1 |   -1 |    0 |    0 |    0 |    0 |
3| 000000000000000000000000000000 | 00000000000001100000 |    3 |
4|    1 |    2 | -99  -99  -99  -99  -99  -99  -99  -99  -99  -99 |
5|    0 |   16 |    0 |    0 |    0 |    0 |    0 |    0 |    0 |
```
It is also recommendet to add information about the variables into the ```STASHmaster-meta.conf```.

You now need to edit your UM app to point to this STASHmaster. This is done by setting the STASHMASTER environment variable to the current directory where the FORTRAN executable is being run (the task working directory).
```config
[env]
STASHMASTER=.
```
Either manually add this variable to the [env] section of the app using a text editor, or open the app in rose edit, navigate to env &rarr; Runtime controls. From the Page drop-down menu select Add  &rarr;  Add latent variable &rarr; STASHMASTER
The ancillary field for the land-area fraction is readable and located in ```/short/public/mb6059/UM_ancils/lfrac```. 
The other primary fields can be initialized with 0. 
For more information on how to use rose visit the official [documentation site](https://code.metoffice.gov.uk/doc/um/latest/um-training/index.html).


### Updating the code
To add new variables new pointers and fields have to be set to the code. 
#### Add new pointers
The first step is to define the new pointers in ```top_level/control/atm_d1_indices_mod.F90```:
```fortran
INTEGER :: jcoastal        ! Information about coastal points
INTEGER :: jdlfrac         ! Land fraction in grid box (2D field)
INTEGER :: jthc            ! Thermal heating contrast land-ocean
INTEGER :: jwindspeed      ! Windspeed for sea-breeze detection
INTEGER :: jwind_dir       ! Wind direction for sea-breeze detection
INTEGER :: jtsurf          ! Surface temperature
```
Then add the new fields in ```control/top_level/atm_fields_real_mod.F90```:
```fortran
REAL, POINTER :: coastal(:,:)          ! coastal mask
REAL, POINTER :: dlfrac(:,:)           ! 2D Land area fraction
REAL, POINTER :: thc(:,:)              ! Thermal heating contrast land-ocean
REAL, POINTER :: tsurf(:,:)            ! surface temperature
REAL, POINTER :: windspeed(:,:)        ! Windspeed for sea-breeze diag
REAL, POINTER :: wind_dir(:,:)         ! Wind direction for sea-breeze diag
```
Then point the new pointer to the appropriate section the in the D1 data array. 
This is done in the ```control/top_level/set_atm_fields.F90``` file:
```fortran
use atm_d1_indices_mod, only: jthc, jtsurf, jwindspeed, jwind_dir,jdlfrac
! Sea breeze diagnostics                                                        
thc(tdims%i_start:tdims%i_end, tdims%j_start:tdims%j_end)                 &         
      => d1(jthc : jthc+field_length(theta_points,no_halo,1) -1)                
tsurf(tdims_s%i_start:tdims_s%i_end, tdims_s%j_start:tdims_s%j_end)     &           
      => d1(jtsurf  : jtsurf +field_length(theta_points,single_halo,1) -1)          
windspeed(tdims%i_start:tdims%i_end, tdims%j_start:tdims%j_end)           &         
      => d1(jwindspeed :jwindspeed+field_length(theta_points,no_halo,1) -1)         
wind_dir(tdims%i_start:tdims%i_end, tdims%j_start:tdims%j_end)           &          
      => d1(jwind_dir : jwind_dir+field_length(theta_points,no_halo,1) -1)          
coastal(tdims_l%i_start:tdims_l%i_end, tdims_l%j_start:tdims_l%j_end)     &         
      => d1(jcoastal  : jcoastal +field_length(theta_points,extended_halo,1) -1)
dlfrac(tdims%i_start:tdims%i_end, tdims%j_start:tdims%j_end)                 &  
      => d1(jdlfrac : jdlfrac+field_length(theta_points,no_halo,1) -1)
```
Then associate the pointers with the stach items, this is done in ```control/top_level/set_atm_pointers.F90```:
```fortran
use atm_d1_indices_mod, only: jthc, jtsurf, jwindspeed, jwind_dir,jdlfrac
! Sea breeze detection variables                                                
jdlfrac                          = si(301, Sect_No, im_index)                   
jthc                             = si(594, Sect_No, im_index)                   
jwindspeed                       = si(595, Sect_No, im_index)                   
jwind_dir                        = si(596, Sect_No, im_index)                   
jtsurf                           = si(597, Sect_No, im_index)                   
jcoastal                         = si(598, Sect_No, im_index)
```
#### Add the output field
The output field is a diagnostic field, which not preserved across timesteps. 
Diagnostic fields can be added in ```control/top_level/atmos_physics2_alloc.F90```.
To add the output field of the sea-breeze detection add the following lines to ```control/top_level/atmos_physics2_alloc.F90```.
```fortran
real, allocatable :: sb_con(:,:)
allocate(sb_con(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end))
```

#### Edit atm_step_4A.f90
The creation of the mask that defines coastal areas involves the call of the ```swap_bounds``` routine 
to comunicate across processing nodes. ```swap_bounds``` should be involved on a high level stage, like in ```atm_step_4a.f90```, 
which performes the integration of the atmospheric model. Add the following code to ```control/top_level/atm_step_4A.F90```.
```fortran
!$omp  parallel default(none)                                          &    
!$omp& shared( wdims,tsurf, tstar, r_u_p2, r_v_p2_n )                  &    
!$omp& private( i, j)                                                       
!$omp do schedule( static )                                                 
!Copy the tstar to tsurf                                                    
do j = wdims%j_start, wdims%j_end                                           
  do i = wdims%i_start, wdims%i_end                                         
    tsurf(i,j) = tstar(i,j)                                                 
  endd o                                                                    
enddo                                                                       
!$omp enddo                                                                
!$omp end parallel                                                           
! Call swap bounds for the surfce temperature                               
                                                                            
call swap_bounds(tsurf,                                              &          
               tdims_s%i_len - 2*tdims_s%halo_i,                     &          
               tdims_s%j_len - 2*tdims_s%halo_j,1,                   &          
               tdims_s%halo_i, tdims_s%halo_j, fld_type_w,.FALSE.)          
! Call the get_edges routine to get information about the coastline         
call get_edges(coastal,ice_fraction,dlfrac)
```
Add this code-snipset *before* the call of ```atmos_physics2```. The call 
agruments of atmos_physics2 also have to be changed. Add the following lines to the call of ```atmos_physics2:```
```fortran
! Sea-breeze detection
thc, tsurf, windspeed, wind_dir,coastal
```
#### Edit ```atmos_physics2.f90```
You also have to change the argument list in ```control/top_level/atmos_physics2.F90```.
Add the following lines to ```control/top_level/atmos_phycics2.F90``` in order to declare 
the new variables:
```fortran
use sea_breeze_diag_mod, only : seabreeze_diag
.
.
.
real intent(in)::                                                      &           
    thc(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),          &           
                    ! Thermal heating contrast land-ocean                             
    tsurf(tdims_s%i_start:tdims_s%i_end,tdims_s%j_start:tdims_s%j_end),&           
                    ! Surface temperature with halo                                   
    windspeed(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),    &           
                    ! Windspeed at sea-breeze timestep                                
    wind_dir(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)                     
                    ! Wind direction at sea-breeze timestep
real, intent(inout) ::                                                 &       
           coastal (tdims_l%i_start:tdims_l%i_end,                     &       
             tdims_l%j_start:tdims_l%j_end)
```
Then call the sea-breeze parametrization *after* ```v_to_p``` and ```u_to_p``` and *before* 
the convection scheme (```ni_conv_ctl```) gets called:
```fortran
!$omp  parallel default(none)                                 &
!$omp shared(timestep,timestep_number,pdims,tdims,p,u_p,v_p,tdims_s,tdims_l, &
!$omp        tsurf,coastal, windspeed, wind_dir, thc, sb_con, Error_code)
  ! Before calling the diganosis for moist convection call the sea-breeze
  ! diagnosis:
  call seabreeze_diag( timestep, timestep_number, &
                      p( pdims%i_start:pdims%i_end, pdims%j_start:pdims%j_end, &
                         pdims%k_start:pdims%k_end ), &
                      u_p( pdims%i_start:pdims%i_end, pdims%j_start:pdims%j_end, &
                           pdims%k_start:pdims%k_end ),  &
                      v_p( pdims%i_start:pdims%i_end, pdims%j_start:pdims%j_end, &
                           pdims%k_start:pdims%k_end ),  &
                      tsurf(tdims_s%i_start:tdims_s%i_end, &
                            tdims_s%j_start:tdims_s%j_end), &
                      coastal(tdims_l%i_start:tdims_l%i_end, &
                                    tdims_l%j_start:tdims_l%j_end), &
                      windspeed(tdims%i_start:tdims%i_end, &
                                tdims%j_start:tdims%j_end), &
                      wind_dir(tdims%i_start:tdims%i_end, &
                               tdims%j_start:tdims%j_end), &
                      thc(tdims%i_start:tdims%i_end, &
                          tdims%j_start:tdims%j_end), &
                      sb_con(pdims%i_start:pdims%i_end, &
                             pdims%j_start:pdims%j_end), Error_code )
!$omp end parallel
```
## Copy the source file
You should copy the ```sea_breeze_diag.F90``` in the ```UM/vn10.7``` sub-folder 
on the [git-hub](https://github.com/antarcticrainforest/seabreeze_param) repository into
any appropriate directory in the UM source code folder structure. The new source code 
will be picked up automatically by the compiling process.

Before submitting any jobs you have to add and commit the changes by:
```bash
$: fcm add sea_breeze_diag.F90
$: fcm commit
```

This should conclude the adaptation of the UM code to implement the sea-breeze 
parametrization into the UM framework.

If you are considering implementing this routine take a look to the  [about](/zz_about) 
section to learn more on who to contribute and improve this project. You are encouraged 
to get in touch via [GitHub](https://github.com/antarcticrainforest/seabreeze_param). 
Bugs should be reported either on the GitHub [issues](https://github.com/antarcticrainforest/seabreeze_param/issues) 
pages or by sending an [email](mailto:martin.bergemann@monash.edu) to the author of this page.

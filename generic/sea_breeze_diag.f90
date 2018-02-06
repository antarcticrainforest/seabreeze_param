!--------------------------------------------------------------------
! Sea Breeze triggering module
!--------------------------------------------------------------------
! MODULE : sea_breeze_diag
!
!> \author{Martin Bergemann, Monash University, Australia}
!
!> Description:
!! This module should detect potential sea-breeze conditions, soley
!! based on large-scale atmosphereric contions (wind, thermal heating
!! contrast).
!
!
! Revision History:
! 2017-08-28 - Initial Version
!
!--------------------------------------------------------------------

module sea_breeze_diag_mod
    implicit none
    !use UM_ParParams
    character(len=*), parameter, private :: ModuleName='sea_breeze_diag_mod'


contains
   !---------------------------------------------------------------------------
   !> @author
   !> Martin Bergemann, Monash University, Australia
   !
   !> Description:
   !> It applies serveral thresholds for change in large-scale
   !! windspeed and direction, overall windspeed and temperature
   !! difference between land and ocean to determine favorable
   !! sea-breeze conditions.
   !
   !> @brief
   !> sea-breeze conditions if : \f$ f = |\Delta \vec{V}| < a \&
   !!      |\Delta \alpha| < b \& |\vec{V}| < c \& |\Delta T| > d \f$
   !> the binary output \f$ f \f$ is then scaled by
   !!      \f$ f_{sb} = f \cdot (|\vec{V}| \cdot |\Delta T|)^{-1} \f$
   !> to get continous vlues.
   !
   ! Revision History:
   ! 2017-08-28 - Initial Version
   !
   !> @param[in] row_lenght, row : processed field dims and subset
   !> @param[in] bl_levels, land_points, p , P_theta_lev, exner_rho
   !! rho_only, rho_theta, z_full, z_half : vertical model grid
   !> @param[in] u,v : windsppeeds
   !> @param[in] theta: potential temperature
   !> @param[out] sb_con : scaled sea-breeze conditions
   !
   !---------------------------------------------------------------------------

subroutine seabreeze_diag(timestep,timestep_number, &
    p, u, v, theta, mask, z, sigma, windspeed, winddir, thc, sb_con)
  ! Definitions of prognostic variable array sizes
  ! Timestep number
  ! use timestep_mod, only: timestep_number

  implicit none

  !----------------------------------------------------------------------------
  ! Subroutine args with intent in
  !----------------------------------------------------------------------------
   integer, intent(in) ::                &
      timestep_number                      !> @var the number of the timestep

  real, intent(in) ::                   &
      p( : ,  : , : ), &    !> @var P  on rho levels (Pa)
      u( : , : , : ), &     !> @var u on rho levels (Pa)
      v( : , : , : ), &     !> @var v on rho levels (Pa)
      theta( :, :),&   !> @var Surface theta (Kelvin)
      timestep                             !> @var timestep (seconds)
    real, intent(in)     ::             &
      mask( :, : ), &    !> var distance form the coast
      z(: , :) , &       !> var surface height
      sigma(: , :)       !> var std of surface height

   !---------------------------------------------------------------------------
   ! Subroutine args with intent inout
   !---------------------------------------------------------------------------
   real, intent(inout)  ::                &
   sb_con( :, :), &     !> @var Sea breeze strength
   windspeed( :, :), &  !> @var Windspeed
   winddir( :, : ), &   !> @var Wind direction in deg
   thc( :, : )          !> @var Thermal heating contrast

   !---------------------------------------------------------------------------
   ! Local variables
   !---------------------------------------------------------------------------
   character (len=*), parameter ::  RoutineName = 'sea_breeze_diag'
   real, dimension(size(theta,1),size(theta,2)) :: t0 !Surface temp. (tmp array)
   integer ::                                      &
   nlats,nlons,model_lev,     & !> @var local number of lons and lats,model_lev
   elat,slat,elon,slon,       & !> @var first and last lon/lat
   k_start,k_end,             & !> @var start and end point of z-vector
   i,j,iu,ju,il,jl,           & !> @var local Loop counter (horiz. field index).
   ii, jj,nn,                 & !> @var local compressed array counter.
   k                            !> @var local Loop counter (vert. level index).
   real ,dimension(size(thc,1),size(thc,2))::   &
   windspeed_abs,             & !> @var local wind speep (abs value)
   thc_abs,                   & !> @var abs value of  thermal heating contrast
   n_windspeed,               & !> @var updated var. of windspeed
   n_thc,                     & !> @var updated var. of thermal heating contrast
   n_winddir,                 & !> @var updated var of the winddirection
   dwinddir,                  & !> @var change in wind direction
   dwindspeed,                & !> @var change in wind speed
   smod_sigma,                & !> @var sigmoid of std orography
   mwindspeed                   !> @var mean wind speed
   real :: mul,               & !> @var multiplyer for coastal land,sea points
     scale_wind,scale_thc       !> @var scaling factors for wind and thc
   integer           :: p_lev         ! target pres. level for wind
   real ::         &
   p_lev_diff,                &  !> @var tempor. variable to get the
                                 !! pressure level of the wind
   T_l,                       &  !> @var avg. Temp over land points
   T_s,                       &  !> @var avg. Temp over sea points
   n_l,                       &  !> @var number of coast. land points
   n_s                           !> @var number of coast. sea points

   !---------------------------------------------------------------------------

   !---------------------------------------------------------------------------
   ! Model constants aka tuning parameters:
   !---------------------------------------------------------------------------
   real, parameter    ::  &
   rad2deg =  57.2957,         & ! Real to Degrees
   target_plev = 100 * 700.,   & ! The wind pres. level in Pa
   thresh_wind = 11.,          & ! windspeed threshold
   thresh_winddir = 90.,       & ! delta wind direct.  thresh.
   thresh_windch = 5.,         & ! delta windspeed threshold
   thresh_thc    = 0.75,       & ! thermal heating con. thresh.
   target_time   = 6.*60**2      ! considered time period [sec]
   real, parameter ::  &
   maxdist = 180.                ! maximum distance that should be influenced by
                                 ! sea-breeze conditions [km]
   real, parameter :: gmma = -0.0060956     !K/m moist adiabatic lapse-rate
   logical :: found                         !switch to indicate if coast was found


   !---------------------------------------------------------------------------
   ! 1.0 Initialisation
   !---------------------------------------------------------------------------
   k_start = 1
   k_end = size(p,3)
   nlats = size(p,2)
   nlons = size(p,1)
   model_lev    = k_end - k_start + 1
   slon = 1
   elon = size(p,1)
   slat = 1
   elat = size(p,2)
   ! only call compute_chunk_size if compiling with OMP
   ! The procedure call is protected by the optional compile sentinel
   !$ CALL compute_chunk_size(1,model_lev,k_start,k_end)
   !$ CALL compute_chunk_size(pdims%i_start,pdims%i_end,slon,elon)
   !$ CALL compute_chunk_size(pdims%j_start,pdims%j_end,slat,elat)
    !--------------------------------------------------------------------------
    ! 1.1 Calculate a theoretical temperature at sea level from moist adiabatic
    !     descent. This is useful for calculating the temperature contrast in
    !     coastal areas with steep terrain,
    !    The calculation is derived from the definiton of the moist
    !    adiabatic laps rate gamma = dT/dz --> T0 = T1 - gamma * z1
    !--------------------------------------------------------------------------
    call sigmoid(sigma,smod_sigma)
    t0 = theta - ( gmma * z * smod_sigma )


  do i = slat,elat
    do j = slon,elon
      !n_s = count(mask(j-1:j+1,i-1:i+1),kind=jprb)/real(9,kind=jprb)
      !n_s = sum(mask(j-1:j+1,i-1:i+1))
      if ( abs(mask(j,i)) > maxdist ) then
        ! Coast is too far away
        sb_con(j,i) = 0.0
      else
        ! -----------------------------------------------------------------
        ! 2.1 Calculate thermal heating contrast
        !------------------------------------------------------------------

        if ( mask(j,i) >= 0.0 ) then ! Land point
          mul = 1
        else
          mul = -1
        endif
      
      nn = 1 !The search radius to find coastal land AND ocean points
      !Increase the search Radius until we have found coasta land AND
      !ocean points (false =- .true.)
      do while (found .eqv. .false.)
        n_l = 0 ! Set all lil' helper values to zero
        n_s = 0
        T_l = 0
        T_s = 0
        ! Lop through the neighborhood of the point, check if land or ocean
        ! and calculate the average land and ocen temp T_{l,o}/n_{l,o}
        do ii=i-nn,i+nn
          do jj=j-nn,j+nn
            if (mask(jj,ii) >= 0.0) then !Land point
              T_l = T_l + t0(jj,ii)
              n_l = n_l + 1
            else
              T_s = T_s + t0(jj,ii)
              n_s = n_s +1
            endif
          enddo
        enddo
        if (n_s > 0 .and. n_l > 0 ) then ! land & ocean points are found
          found = .true.                 ! no need to increas search rad.
        else
          nn = nn + 1                    !land or ocean points missing
        end if
      end do ! while loop
      ! calculate thermal heating contrast
      n_thc(j,i) = mul * ((T_l/n_l) - (T_s/n_s))

        !------------------------------------------------------------------
        ! 2.1 Calculate the wind speed and direction
        !------------------------------------------------------------------

        !Get model level for the target pressure
        p_lev = int(minloc(abs(p(j,i,:) - target_plev),1))
        ! Calculate wind speed
        n_windspeed(j,i) = sqrt(u(j,i,p_lev)**2 + v(j,i,p_lev)**2)
        ! And the direction of the wind in a 2D plane
        n_winddir(j,i) = atan2(-1*u(j,i,p_lev),-1*v(j,i,p_lev)) * rad2deg

        !------------------------------------------------------------------
        ! 2.2 Check if seabreeze conditions should be updated
        !     only the first timestep or every target_time seconds
        !------------------------------------------------------------------

        ! If it is the first timestep, neglect wind change and direction
        if ( timestep_number < 2 ) then
          thc(j,i) = n_thc(j,i)
          winddir(j,i) = n_winddir(j,i)
          windspeed(j,i) = n_windspeed(j,i)
        endif
          
        ! Calculate avg. thermal heating contrast, and wind conditions
        thc_abs(j,i) = abs(n_thc(j,i))
        mwindspeed(j,i) = (windspeed(j,i) + n_windspeed(j,i)) / 2.
        dwindspeed(j,i) = abs(windspeed(j,i) - n_windspeed(j,i))
        ! And get the change in wind direction
        dwinddir(j,i) = abs(modulo((&
                        (winddir(j,i) - n_winddir(j,i)) + 180.),360.) - 180.)
        ! Now apply the threshold rules
        if (  dwinddir(j,i)   < thresh_winddir .and. &
              dwindspeed(j,i) < thresh_windch  .and. &
              mwindspeed(j,i) < thresh_wind    .and. &
              thc_abs(j,i)    > thresh_thc ) then ! All categories are met
          ! Apply the scaling
          scale_wind = (thresh_wind-mwindspeed(j,i)) / max(real(1),mwindspeed(j,i))
          scale_thc =  (thc_abs(j,i) - thresh_thc) / n_thc(j,i)
          sb_con(j,i) = scale_thc*scale_wind
        else
          sb_con(j,i) = 0.0 !Something was missing, no sea breeze
        endif
        !update thermal heating contrast and windspeed
        windspeed(j,i) = n_windspeed(j,i)
        thc(j,i) = n_thc(j,i)
        
        if ( modulo(real(timestep_number)*timestep,target_time)<0.0001 ) then
          winddir(j,i) = n_winddir(j,i) !Update the windspeed only every six hours
        endif
      endif ! End of Check for coastal points if
    end do ! End of row_length loop
  end do ! End of rows loop

 end subroutine seabreeze_diag

  subroutine get_edges(mask,icefrac,landfrac,halo_size)
  !> Description 
  !!: get_edges applies a modified binary version of a sobel operator 
  !! to extract the coastline from a given land sea mask
    use halo_exchange_mod, only : swap_bounds

    implicit none
    character (len=*), parameter ::  RoutineName = 'get_edges'
    integer :: nlons, nlats
    integer,intent(in) :: halo_size
    !----------------------------------------------------------------------
    ! Variables with intent out
    !----------------------------------------------------------------------
    real , intent(out), dimension(:,:) ::  &
      mask( : , : )     !> var land-sea mask + sea ice
    !----------------------------------------------------------------------
    ! Variables with intent in 
    !----------------------------------------------------------------------
    real , intent(in), dimension(:,:) :: &
      landfrac,      &! land-area fraction
      icefrac         ! sea-ice fraction

    !----------------------------------------------------------------------
    ! Local Variables
    !----------------------------------------------------------------------
    real, allocatable :: &
      coast( :, :)          ! helper array for coast area
    integer :: i,j,x,y,slon,elon,slat,elat
    integer , dimension(3,3) :: weight   !Sobel operator
    real :: px, & !Gradient in x
            py, & !Gradient in y
            p     ! sqrt(px**2+py**2)
    ! Define the sobel operator
    weight = reshape ((/-1,-2,-1, &
                         0, 0, 0, &
                         1, 2, 1 /),shape(weight))
   nlons = size(landfrac,1)
   nlats = size(landfrac,2)
   allocate( coast(nlons,nlats) )
   !----------------------------------------------------------------------
   ! Get multi thread info from openmp
   !----------------------------------------------------------------------
   ! call compute_chunk_size(tdims%i_start,tdims%i_end,slon,elon)
   ! call compute_chunk_size(tdims%j_start,tdims%j_end,slat,elat)

   !----------------------------------------------------------------------
   ! Create a new land-sea mask from landfrac and ice fract
   !----------------------------------------------------------------------
   !$omp parallel default(shared) private(x,y)
   !$omp do schedule(static)
   do y = 1,nlats
    do x = 1,nlons
      if ( icefrac(x,y) <= 0.2 ) then
        if ( landfrac(x,y) >= 0.5 ) then
          mask(x+halo_size,y+halo_size) = 1.
        else
          mask(x+halo_size,y+halo_size) = 0.
        endif
      else
        if (landfrac(x,y) + icefrac(x,y) >= 0.5 ) then
          mask(x+halo_size,y+halo_size) = 1.
        else
          mask(x+halo_size,y+halo_size) = 0.
        endif
      endif
    enddo
   enddo
   !$omp end do nowait
   !$omp end parallel
    call swap_bounds(mask, halo_size)
   !----------------------------------------------------------------------
   ! Create the coastal mask by applying sobels algorithm
   !----------------------------------------------------------------------
   !$omp parallel shared(tdims,weight,mask,coast) private (x,y,px,py,i,j,p)
   !$omp do schedule(static)
   do y = 1,nlats
    do x = 1,nlons
      px = 0.0
      py = 0.0
      do j=-1,1
        do i=-1,1
          px = px + weight(i + 2,j + 2) * mask(x+j,y+i) ! Apply a sobel operator
          py = py + weight(j + 2,i + 2) * mask(x+j,y+i)
        enddo
      enddo
      p = sqrt(px**2+py**2)
      !Sobel operator is gray scale make it binary
      if (p == 0) then
        coast(x,y) = 0.
      else
        coast(x,y) = 1
      endif
    enddo
   enddo
   !$omp end do nowait
   !$omp end parallel
   mask( 1+halo_size:nlons+halo_size,1+halo_size:nlats+halo_size ) = coast(:,:)

   call swap_bounds(mask, halo_size)

  end subroutine get_edges

  subroutine get_dist(coast,landfrac,lon,lat,maxdist,cdist,halo_size)
  !> Description
  !! get_dist calculates the distance to the next
  !! coastal grid point within a given search radius around the this 
  !! this coastal grid point

  !! f2py declarations
    !f2py integer       , intent(in)        :: nlons,nlats
    !f2py real          , intent(in)        :: coast, landfrac
    !f2py real          , intent(in)        :: lon,lat
    !f2py real optional , intent(in)        :: maxdist = 180
    !f2py real          , intent(out)       :: cdist


  implicit none
    integer, intent(in)  :: halo_size
    real,  intent(in), dimension(:) :: lon
    real,  intent(in), dimension(:) :: lat
    real,  intent(in), dimension(:,:) :: coast ! Coastal data
    real,  intent(in), dimension(:,:) :: landfrac ! Land-area fraction

    real,intent(out), dimension(size(coast,1),size(coast,2)) :: & 
      cdist !Shortest distance to a coastal point (output)

    integer :: nlats, nlons

    integer :: jj,ii,i,j,xx,yy,tlat      !Helper variables
    real :: dlam,dphi,a,c,dx,maxdist,l1,l2 !more helper variables
    real, dimension(size(lon,1)) :: lam1         ! lon in radians
    real, dimension(size(lat,1)) :: phi1         ! lat in radians
    real, parameter :: R = 6370.9989       ! Earth's radius in km
    real, parameter :: pi = 3.1415926      ! pi
    real, parameter :: r2d    = 180.0/pi   ! radian to degree
    real, parameter :: d2r    = pi/180.0   ! degree to radian
    
    nlons = size(lon,1)
    nlats = size(lat,1)
    !Initialize the distance with 12000 km
    do i=1,nlats
      do j=1,nlons
        cdist(j,i) = 12000.
      enddo
    enddo

    do i=1,nlats
      do j=1,nlons
        if ( coast(j,i) > 0.) then ! Check for coastal point
          do ii=-halo_size,halo_size ! Loop through neighborhood ii
              yy = i+ii !Index for lats is bounded
              dphi = phi1(i) - phi1(yy)
            do jj=-halo_size,halo_size ! Loop through neighborhood jj
             ! Periodic boundaries for lats
              xx = j+jj
              if (lon(j) > 180) then !Correct lon
                l1 = d2r * ( lon(j) - 360. )
              else
                l1 = d2r * (lon(j))
              endif ! end-if correct lon
              if (lon(xx) > 180) then ! Correct lon
                l2 = d2r * ( lon(xx) - 360. )
              else
                l2 = d2r * lon(xx)
              endif ! end-if correct lon
              dlam = l1 - l2
              a = sin(dphi/2)**2+(cos(phi1(i))*(cos(phi1(yy))*sin(dlam/2)**2))
              c = R*2*atan2(sqrt(a),sqrt(1-a))+0.5
              if ( c < abs(cdist(xx,yy)) ) then ! Check for min coastal dist
                if (landfrac(xx,yy) > 0.0) then ! Check sea or land point
                  cdist(xx,yy) = c
                else
                  cdist(xx,yy) = -c
                endif ! end-if sea or land point
              endif ! end-if min coastal dist
            enddo ! end-do jj neighborhood
          enddo !end-do ii neighborhood
        endif !end-if coastal point check
        if ( abs(cdist(j,i)) > 2*maxdist) cdist(j,i) = 12000.
      enddo !loop lon
  enddo ! loop lat

  end subroutine get_dist
  
  subroutine sigmoid(ary,sm)
    implicit none
    real, dimension(:,:) ,intent(in) :: ary
    real, dimension(size(ary,1),size(ary,2)) ,intent(out) :: sm
    integer  :: nlons,nlats
    real   :: mean,std,var,r
    integer :: i,j
    nlons = size(ary,1)
    nlats = size(ary,2)
    mean = sum(ary)/(nlons*nlats)
    var = 0

    !$omp parallel default(shared) private(j,i)
    !$omp do schedule(static)
    do i=1,nlats
      do j=1,nlons
        var = var + (ary(j,i)-mean)**2
      enddo
    enddo
    !$omp end do nowait
    !$omp end parallel
    std = 2/sqrt( var / (nlons*nlats) )
    r = ( maxval(ary) - minval(ary) ) / 4.
    sm = 1 / ( 1 + exp(-std*(ary-r) ) )
  end subroutine sigmoid
end module sea_breeze_diag_mod

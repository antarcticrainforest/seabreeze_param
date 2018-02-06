!-----------------------------------------------------------------------------
! Sea Breeze triggering module
!------------------------------------------------------------------------------
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
!------------------------------------------------------------------------------


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

subroutine diag(timestep_number, &
    p, z , std, theta, v, u, cdist, windspeed, winddir, thc,  &
    target_plev,thresh_wind,thresh_winddir,thresh_windch, &
    thresh_thc,target_time,maxdist,timestep,nps,nlons,nlats,output)

  !$ use omp_lib
  implicit none


  integer :: nps
  integer :: nlats
  integer :: nlons
  real ::                                &
    target_plev,                         & ! target pressure level
    thresh_wind,                         & ! wind speed threshold
    thresh_winddir,                      & ! wind direction change threshold
    thresh_windch,                       & ! wind speed change threshold
    thresh_thc,                          & ! thermal heating contrast
    target_time,                         & ! considered time period
    timestep,                            & ! time step of the input data
    maxdist                                ! maximum distance for sea-breezes
  real, dimension (nlons,nlats) :: t0      ! dry adiabatic temp at sea-level
  !----------------------------------------------------------------------------
  ! args with intent in
  !----------------------------------------------------------------------------
   integer ::                            &
      timestep_number                      ! the number of the timestep
    real, dimension(nps) ::              &
      p                                    ! P on rho levels (Pa)
    real, dimension(nlons,nlats,nps) ::  &
      v                                    ! v on rho levels (Pa)
    real, dimension(nlons,nlats,nps) ::  &
      u                                    ! u on rho levels (Pa)
    real, dimension(nlons,nlats) ::      &
      cdist,                             & ! eclidian dist. to the next coast
      theta,                             & ! surface temperature
      z,                                 & ! surface elevation
      std,                               & ! std of subgrid orography
      sb_con,                            & ! Sea breeze strength
      windspeed,                         & ! Windspeed
      winddir,                           & ! Winddirection in deg
      thc                                  ! Thermal heating contrast

   !---------------------------------------------------------------------------
   ! args with intent out
   !---------------------------------------------------------------------------
    real, dimension(nlons,nlats,4) ::    &
      output                              !f2py supports only one output array
                                          !this is where that should be returned
                                          !is stored

   !---------------------------------------------------------------------------
   ! Local variables
   !---------------------------------------------------------------------------
   integer ::                             &
   i,j,                                   & !Loop counter (horiz. field).
   ii, jj, nn, ki, kj                       !Local counters.
   real, dimension(nlons,nlats)::         &
     thc_abs,                             & !abs value of thermal heating contr.
     n_windspeed,                         & !update of windspeed
     n_thc,                               & !update of thermal heating contrast
     n_winddir,                           & !update of the winddirection
     dwinddir,                            & !change in wind direction
     dwindspeed,                          & !change in wind speed
     smod_std,                            & !sigmoid of subgrid orography
     mwindspeed                             !avg. wind speed
   integer ::                             &
     p_lev                                  !target pres. level for wind
   real ::                                &
     T_l,                                 & !avg. Temp over land points
     T_s,                                 & !avg. Temp over sea points
     n_l,                                 & !number of coast. land points
     n_s,                                 & !number of coast. sea points
     mul,                                 & !multipl. for coastal land,sea point
     scale_wind,                          & !scaling factor for wind
     scale_thc                              !scaling factor for therm. heat. con
   real, parameter :: rad2deg = 57.2957
   real, parameter :: gmma = -0.0060956     !K/m moist adiabatic lapse-rate
   logical :: found                         !switch to indicate if coast was found

    !! f2py declarations
    !f2py integer       , intent(in)        :: timestep_number
    !f2py integer       , intent(in)        :: nps
    !f2py integer       , intent(in)        :: nlats
    !f2py integer       , intent(in)        :: nlons
    !f2py real          , intent(in)        :: p, u, v, theta, z, std
    !f2py real          , intent(in)        :: cdist
    !f2py real          , intent(in)        :: thc, windspeed, winddir
    !f2py real optional , intent(in)        :: target_plev=700.,thresh_wind=11,
    !f2py real optional , intent(in)        :: thresh_winddir=90.,thresh_windch=5.
    !f2py real optional , intent(in)        :: thresh_thc=0.75,target_time=6
    !f2py real optional , intent(in)        :: timestep = 24
    !f2py real optional , intent(in)        :: maxdist = 180
    !f2py real          , intent(out)       :: output
    !--------------------------------------------------------------------------
    ! 1.0 Convert timestep and target_time to seconds and target_plev to Pa
    !--------------------------------------------------------------------------
    timestep = timestep * 60. ! minutes to sec
    target_time = target_time * 60.**2 !hours to sec
    target_plev = target_plev * 100.

    !--------------------------------------------------------------------------
    ! 1.1 Calculate a theoretical temperature at sea level from moist adiabatic
    !     descent. This is useful for calculating the temperature contrast in
    !     coastal areas with steep terrain,
    !    The calculation is derived from the definiton of the moist
    !    adiabatic laps rate gamma = dT/dz --> T0 = T1 - gamma * z1
    !--------------------------------------------------------------------------
    call sigmoid(std,nlons,nlats,smod_std)
    t0 = theta - ( gmma * z * smod_std )

    !$omp parallel default(private) shared(nlats,nlons,sb_con,cdist,n_thc, t0, & 
    !$omp target_plev,target_time,u,v,p,timestep_number,thresh_winddir,thc_abs, &
    !$omp thresh_windch,thresh_wind,thresh_thc,winddir,n_winddir,windspeed,n_windspeed, &
    !$omp thc, output, mwindspeed, dwindspeed,timestep,maxdist)
    !$omp do schedule(static)
    do i=1,nlats-1
      do j=1,nlons
      found = .false.
      !------------------------------------------------------------------------
      ! 2.0 Check for coastal point
      !------------------------------------------------------------------------
        if ( abs(cdist(j,i)) > maxdist ) then
          ! Coast is too far away, set missing value
          sb_con(j,i) = 2.0E20
        else
          ! There is a coastal point in the vicinity, lets do the calculation

          ! -----------------------------------------------------------------
          ! 2.1 Calculate thermal heating contrast
          !------------------------------------------------------------------

          if ( cdist(j,i) >= 0.0 ) then ! Coastal land point
            mul = 1
          else
            mul = -1
          endif !End land/sea check

          nn = 1 !The search radius to find coastal land AND ocean points

          !Increase the search Radius until we have found coasta land AND
          !ocean points (false =- .true.)
          do while (found .eqv. .false.)
            n_l = 0 ! Set all lil' helper values to zero
            n_s = 0
            T_l = 0
            T_s = 0
            ! Loop through the neighborhood of the point, check if land or ocean
            ! and calculate the average land and ocen temp T_{l,o}/n_{l,o}
            do ii=i-nn,i+nn
             do jj=j-nn,j+nn

                ki = min(max(1,ii),nlats) !Index for lats is bounded
                kj = max(1,modulo(jj,nlons)) ! Periodic boundaries for lons

                if (cdist(kj,ki) >= 0.0) then !Coastal land point
                   T_l = T_l + t0(kj,ki)
                   n_l = n_l + 1
                else
                   T_s = T_s + t0(kj,ki)
                   n_s = n_s + 1
                endif !End land/sea point check
             enddo
            enddo
            if (n_s > 0 .and.  n_l > 0) then ! land and ocean points are found
              found = .true.                 ! no need to increas the search
            else                             ! radius
              nn = nn + 1                    ! land or ocean point missing
            endif                            ! increase the search radius
          enddo !while loop

          ! Calculate the avg thermal heating contrast
          n_thc(j,i) = mul * ((T_l/n_l) - (T_s/n_s))

          !--------------------------------------------------------------------
          ! 2.2 Calculate the wind speed and direction
          !--------------------------------------------------------------------

          !Get model level for the target pressure
          p_lev = int(minloc(abs(p - target_plev),1))

          ! Calculate wind speed
          n_windspeed(j,i) = sqrt(u(j,i,p_lev)**2 + v(j,i,p_lev)**2)
          ! And the direction of the wind in a 2D plane
          n_winddir(j,i) = atan2(-1*u(j,i,p_lev),-1*v(j,i,p_lev)) * rad2deg

          ! If it is the first timestep, neglect wind change and direction
          if ( timestep_number < 2 ) then
            thc(j,i) = n_thc(j,i)
            winddir(j,i) = n_winddir(j,i)
            windspeed(j,i) = n_windspeed(j,i)
          endif !End of first timestep check
          !--------------------------------------------------------------------
          ! 2.3 Calculate avg. thermal heating contrast, and wind conditions
          !--------------------------------------------------------------------
          thc_abs(j,i) = abs(n_thc(j,i))
          mwindspeed(j,i) = (windspeed(j,i) + n_windspeed(j,i)) / 2.
          dwindspeed(j,i) = abs(windspeed(j,i) - n_windspeed(j,i))
          ! And get the change in wind direction
          dwinddir(j,i) = abs(modulo((&
                         (winddir(j,i) - n_winddir(j,i)) + 180.),360.) - 180.)
          !--------------------------------------------------------------------
          ! 2.4 Now apply the threshold rules and calculate the scaling factor
          !--------------------------------------------------------------------
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
          endif ! End threshold test

          !--------------------------------------------------------------------
          ! 2.5 update thermal heating contrast windspeed and direction
          !--------------------------------------------------------------------
          thc(j,i) = n_thc(j,i)

          ! Update the windspeed only every target_time hours
          if ( modulo(real(timestep_number)*timestep,target_time)<0.0001 ) then
            windspeed(j,i) = n_windspeed(j,i)
            winddir(j,i) = n_winddir(j,i)
          endif !End of time check

        endif ! End of Check for coastal points if
        output(j,i,1) = sb_con(j,i)
        output(j,i,2) = t0(j,i)
        output(j,i,3) = windspeed(j,i)
        output(j,i,4) = winddir(j,i)
      end do ! End of row_length loop
    end do ! End of rows loop
    !$omp end do nowait
    !$omp end parallel
end subroutine diag

subroutine sigmoid(ary,nlons,nlats,sm)
  !$ use omp_lib
  implicit none
  integer , intent(in) :: nlons,nlats
  real, dimension(nlons,nlats) ,intent(in) :: ary
  real, dimension(nlons,nlats) ,intent(out) :: sm
  real   :: mean,std,var,r
  integer :: i,j

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

      !-------------------------------------------------------------------------
      ! This code is not intended to wrok, it should only outline how
      ! the sea-breeze parametrization should be implemented to a new
      ! model.
      !
      ! Suppose this program perfroms a model time step and calls all
      ! physics routines. This programm should demostrate how to call
      ! the routines that are involved when the sub-grid sea-breeze
      ! strength is calculated.
      !
      ! The calculation steps are the following:
      ! 1   call get_edges to create the models coastline
      ! 2   perform an exchange of coastal field across boundaries
      !       (boundary swapping)
      ! 3   call get_dist to calculate the distances to the closest
      !      coastal points
      ! 4   perform an exchange of data across boundaries (boundary
      !      swappgin)
      ! 5   call physics
      ! 5.1 from physics call the sea-breeze diagnosis
      !
      !------------------------------------------------------------------------

      program dummy_model
        implicit none
        contains
          subroutine atmos_step()
            use sea_breeze_diag_mod, only : get_edges
            use get_all_fields_mod
            implicit none
            real, dimension(nx,ny) :: cdist
            call get_edges(mask, ice_frac, land_frac, halo_size)
            call get_dist(mask, land_frac, lon, lat, 180, cdist, &
              halo_size)
            call physics(p, u, v, theta, z, sigma ,mask, windspeed, & 
              winddir, thc, sb_con,timestep_number, timestep)
          end subroutine atmos_step

          subroutine physics(p,u,v, theta, z ,sigma, mask, windspeed,& 
              winddir, thc, sb_con, timestep_number, timestep)
            use sea_breeze_diag_mod, only : seabreeze_diag
            implicit none
            real, intent(in), dimension(:,:,:) :: & 
              p, u, v
            real, intent(in), dimension(:,:) :: &
              mask, theta, z, sigma
            real, intent(inout), dimension(:,:) :: &
              thc, windspeed, winddir, sb_con
            real, intent(in) :: timestep
            integer, intent(in) :: timestep_number

            call seabreeze_diag(timestep,timestep_number, &
              p, u, v, theta, z, sigma, mask, windspeed, winddir, thc, & 
              sb_con)
          end subroutine physics
      end program dummy_model

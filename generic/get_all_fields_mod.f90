      !------------------------------------------------------------------------
      ! This dummy module should provide some empty fields. Fields that
      ! are neeede to call the sea-breeze routine. The module is only
      ! for demonstrative purpose. In reality the fields come from
      ! ceratain routine calls performed by the climate model
      module get_all_fields_mod
        implicit none

        integer, parameter ::  nx = 128 , ny = 96 , nz = 56 
        integer, parameter :: halo_size = 4
        real, parameter :: timestep = 24/60.
        integer :: timestep_number
        real, dimension(nx) :: lon
        real, dimension(ny) :: lat
        real, dimension(nx,ny,nz) :: &
          p, u , v
        real, dimension(nx,ny) :: &
          sb_con, land_frac, ice_frac, windspeed, winddir, thc, z, sigma
        real, dimension(nx+halo_size,ny+halo_size) ::  mask, theta
      end module get_all_fields_mod

      !------------------------------------------------------------------------
      ! This module is for demonstrative purpose only, it stands for
      ! various routines about communications across processors to
      ! exchange information about split domain boundaries
      !
      ! In the present case it doesn't do anything
      !------------------------------------------------------------------------
      module halo_exchange_mod
        implicit none

        contains
          subroutine swap_bounds(field,halo_size)
            implicit none

            integer, intent(in) :: halo_size
            real, intent(inout) :: field(:,:)
        end subroutine swap_bounds
      end module halo_exchange_mod

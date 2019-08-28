      !------------------------------------------------------------------------
      ! MODULE : sobel
      !
      !> \author{Martin Bergemann, Monash University, Australia}
      !
      !> Description:
      !! This module contains two subroutines (get_edges and get_dist)
      !
      !> Soubroutines:
      !! get_edges: it applies a modified binary version of a sobel operator 
      !! to extract the coastline from a given land sea mask
      !!
      !! get_dist: this routine calculates the distance to the next
      !! coastal grid point within a given search radius around the this 
      !! this coastal grid point
      !
      !------------------------------------------------------------------------

      subroutine get_edges(lsm,ci,nlons,nlats,coast)
        !! f2py declarations
        !f2py integer       , intent(in)        :: nlats
        !f2py integer       , intent(in)        :: nlons
        !f2py real          , intent(in)        :: lsm,ci
        !f2py real          , intent(out)       :: coast
        
        !$ use omp_lib
        implicit none
        !----------------------------------------------------------------------
        ! Variables with intent in
        !----------------------------------------------------------------------
        integer :: nlats, nlons
        real , dimension(nlons,nlats) ::  &
          lsm,                            & ! Land-sea mask
          ci                                ! Sea-ice data
        !----------------------------------------------------------------------
        ! Variables with inten out
        !----------------------------------------------------------------------
        real , dimension(nlons,nlats) :: coast ! coastal data

        !----------------------------------------------------------------------
        ! Local variables
        !----------------------------------------------------------------------
        integer :: i,j,x,y,xx,yy
        integer , dimension(3,3) :: weight   !Sobel operator
        real :: px, & !Gradient in x
                py, & !Gradient in y
                p,  & ! sqrt(px**2+py**2)
                mm    ! Mask for land or sea point

        real, dimension(nlons,nlats) :: mask
        mask = lsm + ci  !Define a new mask from land-sea mask and sea-ice

        ! Define the sobel operator
        weight = reshape ((/-1,-2,-1, &
                             0, 0, 0, &
                             1, 2, 1 /),shape(weight))
        ! Initialize the output array
        !$omp parallel default(shared) private (x,y,j,i,px,py,yy,xx,mm,p)
        !$omp do schedule(static)
        do y=1,nlats
          do x=1,nlons
            px = 0.0
            py = 0.0
            do j=-1,1
              do i=-1,1
                !Take care about the boundaries
                yy = min(max(1,y+i),nlats) ! Upper and lower bounds a limited
                xx = max(1,modulo(x+j,nlons)) ! East-west bounds are periodic
                if ( mask(xx,yy) > 0.4 ) then ! Land-or sea points
                  mm = 1 ! Land point
                else
                  mm = 0 ! sea point
                endif
                px = px + weight(i + 2,j + 2) * mm ! Apply a sobel operator
                py = py + weight(j + 2,i + 2) * mm
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
      end subroutine get_edges

      subroutine get_dist(coast,mask,lon,lat,nlons,nlats,maxdist,cdist)
      !! f2py declarations
        !f2py integer       , intent(in)        :: nlons,nlats
        !f2py real          , intent(in)        :: coast, mask
        !f2py real          , intent(in)        :: lon,lat
        !f2py real optional , intent(in)        :: maxdist = 180
        !f2py real          , intent(out)       :: cdist

      !$ use omp_lib

      implicit none
        integer  :: nlons,nlats
        real,  dimension(nlons) :: lon
        real,  dimension(nlats) :: lat
        real,  dimension(nlons,nlats) :: coast ! Coastal data
        real,  dimension(nlons,nlats) :: mask  ! Land-sea mask

        real,  dimension(nlons,nlats) :: cdist !Shortest distance to a coastal
                                               ! point (output)

        integer :: k,jj,ii,i,j,xx,yy,tlat      !Helper variables
        real :: dlam,dphi,a,c,dx,maxdist,l1,l2 !more helper variables
        real, dimension(nlons) :: lam1         ! lon in radians
        real, dimension(nlats) :: phi1         ! lat in radians
        real, parameter :: R = 6370.9989       ! Earth's radius in km
        real, parameter :: pi = 3.1415926      ! pi
        real, parameter :: r2d    = 180.0/pi   ! radian to degree
        real, parameter :: d2r    = pi/180.0   ! degree to radian
        integer :: nt
        !The maximum search radius in given is km, this has to be
        !translated into degree, to apply the search radius to gridded
        !data close to the poles a search radius of x-km would become
        !unnecessary large and would make the algorithm inefficient.
        !Therefore the search radius in degrees is bounded at a limit of
        !70 deg N/S. This means the the radius for points north/south of
        !70 becomes slightly smaller, but the algorithm is sufficiently
        !fast.

        tlat = int(minloc(abs(70 - lat),1)) !Get the index of the 70 deg lat
        phi1 = d2r * lat ! lat in rad
        lam1 = d2r * lon ! lon in rad
        dphi = phi1(tlat+1) - phi1(tlat) !differenz at 70 deg
        dlam = lam1(2) - lam1(1) 
        !Calculate the distance in km
        a = sin(dphi/2)**2+(cos(phi1(tlat+1))*(cos(phi1(tlat))*sin(dlam/2)**2))
        dx = R*2*atan2(sqrt(a),sqrt(1-a)) 
        k = int(maxdist / dx) ! set the max number of grid points for
        !the search radius k

        !Initialize the distance with 12000 km
        !$omp parallel default(shared) private (j,i)
        !$omp do schedule(static)
        do i=1,nlats
          do j=1,nlons
            cdist(j,i) = 12000.
          enddo
        enddo
        !$omp end do nowait
        !$omp end parallel

        nt = 0
        !$omp parallel default(shared) private(j,i,ii,jj,xx,yy,l1,l2,dlam,dphi,a,c,nt)
        !$omp do schedule(static)
        do i=1,nlats
          do j=1,nlons
            !$ nt = omp_get_thread_num()
            if ( coast(j,i) > 0.) then ! Check for coastal point
              do ii=-k,k ! Loop through neighborhood ii
                  yy = min(max(1,ii+i),nlats) !Index for lats is bounded
                  dphi = phi1(i) - phi1(yy)
                do jj=-k,k ! Loop through neighborhood jj
                 ! Periodic boundaries for lats
                  xx = modulo(j+jj,nlons)
                  if ( xx == 0 ) xx = nlons
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
                    if (mask(xx,yy) > 0.0) then ! Check sea or land point
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
        !$omp end do nowait
        !$omp end parallel
      end subroutine get_dist

      subroutine get_threads(nt)
        
        !! f2py declarations
        !f2py integer       , intent(out)        :: nt
        !$ use omp_lib
        implicit none
        integer :: nt
        nt = 0
        !$ nt = omp_get_max_threads()


      end subroutine get_threads

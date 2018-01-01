

module sweepsize
!=======================================================================
! Dimension of 1D sweeps.  maxsweep must be as long as the longest of the 
! 3D arrays PLUS the ghost zones: 
! ie, must have maxsweep >= max(imax,jmax,kmax) + 12
!----------------------------------------------------------------------

integer, parameter :: maxsweep=1036 

end module sweepsize

module sweeps      
!=======================================================================
! Data structures used in 1D sweeps, dimensioned maxsweep  (set in sweepsize.mod)
!----------------------------------------------------------------------

use sweepsize

character(len=1) :: sweep                                    ! direction of sweep: x,y,z
integer :: nmin, nmax, ngeom, nleft, nright                  ! number of first and last real zone  
real, dimension(maxsweep) :: r, p, e, q, u, v, w             ! fluid variables
real, dimension(maxsweep) :: xa, xa0, dx, dx0, dvol          ! coordinate values
real, dimension(maxsweep) :: f, flat                         ! flattening parameter
real, dimension(maxsweep,5) :: para                          ! parabolic interpolation coefficients
real :: radius, theta, stheta

end module sweeps


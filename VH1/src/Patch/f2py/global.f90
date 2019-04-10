module global
!=======================================================================
! global variables accessible by (most) anything
!-----------------------------------------------------------------------

 character(len=50) :: prefix              ! prefix for output filenames

 integer :: ncycle, ncycp, ncycm, ncycd  ! cycle number
 integer :: nfile                        ! output file numbers
 integer :: ndim

 real(kind=8) :: time, dt, timem, timep, svel, vdtext
 real(kind=8) :: gam, gamm

 real(kind=8), parameter :: courant = 0.5           ! timestep fraction of courant limit
 real(kind=8), parameter :: pi = 3.1415926535897931 ! shouldn't computers know this?
 real(kind=8), parameter :: xwig = 0.00             ! fraction of a zone to wiggle grid for dissipation
 real(kind=8), parameter :: smallp = 1.0e-30        ! Set small values to prevent divide by zero
 real(kind=8), parameter :: smallr = 1.0e-30
 real(kind=8), parameter :: small  = 1.0e-30

 real(kind=8) :: uinflo, dinflo, vinflo, winflo, pinflo, einflo 
 real(kind=8) :: uotflo, dotflo, votflo, wotflo, potflo, eotflo
      
end module global

module sweepsize
!=======================================================================
! Dimension of 1D sweeps.  maxsweep must be as long as the longest of the 
! 3D arrays PLUS the ghost zones: 
! ie, must have maxsweep >= max(imax,jmax,kmax) + 12
!----------------------------------------------------------------------

integer, parameter :: maxsweep=1048588 

end module sweepsize

module sweeps      
!=======================================================================
! Data structures used in 1D sweeps, dimensioned maxsweep  (set in sweepsize.mod)
!----------------------------------------------------------------------

use sweepsize

character(len=1) :: sweep                                    ! direction of sweep: x,y,z
integer :: nmin, nmax, ngeom, nleft, nright                  ! number of first and last real zone  
real(kind=8), dimension(maxsweep) :: r, p, e, q, u, v, w, g          ! fluid variables
real(kind=8), dimension(maxsweep) :: xa, xa0, dx, dx0, dvol          ! coordinate values
real(kind=8), dimension(maxsweep) :: f, flat                         ! flattening parameter
real(kind=8), dimension(maxsweep,5) :: para                          ! parabolic interpolation coefficients
real(kind=8) :: radius, theta, stheta

end module sweeps

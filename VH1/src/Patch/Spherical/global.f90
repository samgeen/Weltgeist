module global
!=======================================================================
! global variables accessible by (most) anything
!-----------------------------------------------------------------------

 character(len=50) :: prefix              ! prefix for output filenames

 integer :: ncycle, ncycp, ncycm, ncycd  ! cycle number
 integer :: nfile                        ! output file numbers
 integer :: ndim

 real :: time, dt, timem, timep, svel 
 real :: gam, gamm

 real, parameter :: courant = 0.5           ! timestep fraction of courant limit
 real, parameter :: pi = 3.1415926535897931 ! shouldn't computers know this?
 real, parameter :: xwig = 0.00             ! fraction of a zone to wiggle grid for dissipation
 real, parameter :: smallp = 1.0e-15        ! Set small values to prevent divide by zero
 real, parameter :: smallr = 1.0e-15
 real, parameter :: small  = 1.0e-15

 real :: uinflo, dinflo, vinflo, winflo, pinflo, einflo 
 real :: uotflo, dotflo, votflo, wotflo, potflo, eotflo
      
end module global

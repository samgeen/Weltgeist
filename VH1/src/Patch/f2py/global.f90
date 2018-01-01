module global
!=======================================================================
! global variables accessible by (most) anything
!-----------------------------------------------------------------------

 character(len=50) :: prefix              ! prefix for output filenames

 integer :: ncycle, ncycp, ncycm, ncycd  ! cycle number
 integer :: nfile                        ! output file numbers
 integer :: ndim

 real(kind=8) :: time, dt, timem, timep, svel 
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

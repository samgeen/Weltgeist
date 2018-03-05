! Hacked VH1 entry point

module data

! Allow external access to global variables
use global
use sweeps
use sweepsize
use zone

implicit none

! Calculate this every timestep for use by the code at large
!real(kind=8), dimension(maxsweep) :: vols

! Dimensions of grid in x, y and z dimensions
!f2py INTEGER :: imax, jmax, kmax   ! Memory dimensions

! number of first and last real zone
!f2py INTEGER :: ngeomx, ngeomy, ngeomz       ! XYZ Geometry flag
!f2py INTEGER :: nleftx, nlefty, nleftz       ! XYZ Lower Boundary Condition
!f2py INTEGER :: nrightx,nrighty,nrightz      ! XYZ Upper Boundary Condition
! This should be large enough... (shouldn't need to change it, probably...)
!f2py   integer,parameter :: maxsweep=1048588
! fluid variables
!f2py   real(kind=8), dimension(maxsweep) :: r, p, e, q, u, v, w
! coordinate values 
!f2py REAL(kind=8), ALLOCATABLE,DIMENSION(:) :: zxa, zdx, zxc
!f2py REAL(kind=8), ALLOCATABLE,DIMENSION(:) :: zya, zdy, zyc
!f2py REAL(kind=8), ALLOCATABLE,DIMENSION(:) :: zza, zdz, zzc
!f2py REAL(kind=8), ALLOCATABLE,DIMENSION(:,:,:) :: zro
!f2py REAL(kind=8), ALLOCATABLE,DIMENSION(:,:,:) :: zpr
!f2py REAL(kind=8), ALLOCATABLE,DIMENSION(:,:,:) :: zux
!f2py REAL(kind=8), ALLOCATABLE,DIMENSION(:,:,:) :: zuy
!f2py REAL(kind=8), ALLOCATABLE,DIMENSION(:,:,:) :: zuz
!f2py REAL(kind=8), ALLOCATABLE,DIMENSION(:,:,:) :: zfl
 
!f2py   real(kind=8) :: xmin, xmax, ymin, ymax, zmin, zmax
!f2py   real(kind=8) :: time, dt, timem, timep, svel, vdtext 
!f2py   real(kind=8) :: gam

contains

subroutine setup
  implicit none
  call init

end subroutine setup

! A single hydro step
subroutine step
  implicit none
  REAL(kind=8) :: olddt
  
  call dtcon   ! Check constraints on the timestep


  olddt  = dt
  svel   = 0.
  vdtext = 0.

! Alternate sweeps to approximate 2nd order operator splitting
! SAM GEEN - intending only to use sweepx but might as well keep this

                call sweepx
  if(ndim > 1)  call sweepy
  if(ndim == 3) call sweepz

  time  = time  + dt
  timep = timep + dt
  timem = timem + dt

  ! if(ndim == 3) call sweepz
  ! if(ndim > 1)  call sweepy 
  !               call sweepx

  ! time  = time  + dt
  ! timep = timep + dt
  ! timem = timem + dt
  dt = olddt

end subroutine step

end module data

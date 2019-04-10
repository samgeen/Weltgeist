subroutine init

! Sod shock tube problem (a whimpy test) in 1, 2, or 3 dimensions
! 24jan92 blondin
!=======================================================================
! GLOBALS
use global
use zone
use sweepsize

IMPLICIT NONE

! LOCALS
INTEGER :: i, j, k
REAL(kind=8) :: ridt, xvel, yvel, zvel, width, widthz, widthy
REAL(kind=8) :: dleft, pleft, dright, pright, plane

!--------------------------------------------------------------------------------
! Set up geometry and boundary conditions of grid
!
! Boundary condition flags : nleft, nright
!   = 0  :  reflecting boundary condition
!   = 1  :  inflow/outflow boundary condition (zero gradients)
!   = 2  :  fixed inflow boundary condition (values set by dinflo, pinflo, etc.)
!   = 3  :  periodic (nmax+1 = nmin; nmin-1 = nmax)
!
! Geometry flag : ngeom                         |  Cartesian:
!   = 0  :  planar                              |    gx = 0, gy = 0, gz = 0   (x,y,z)
!   = 1  :  cylindrical radial                  |  Cylindrical:
!   = 2  :  spherical   radial             3D= {     gx = 1, gy = 3, gz = 0   (s,phi,z)
!   = 3  :  cylindrical angle                   |
!   = 4  :  spherical polar angle (theta)       |  Spherical:
!   = 5  :  spherical azimu angle (phi)         |    gx = 2, gy = 4, gz = 5   (r,theta,phi)

! If any dimension is angular, multiply coordinates by pi...
if(ngeomy >= 3) then
   ymin = ymin * pi
   ymax = ymax * pi
endif
if(ngeomz >= 3) then
   zmin = zmin * pi
   zmax = zmax * pi
endif

!======================================================================
! Allocate hydro variables in zonemod.f90 (module zone)
allocate(zro(imax,jmax,kmax))
allocate(zpr(imax,jmax,kmax))
allocate(zux(imax,jmax,kmax))
allocate(zuy(imax,jmax,kmax))
allocate(zuz(imax,jmax,kmax))
allocate(zfl(imax,jmax,kmax))
allocate(zgr(imax,jmax,kmax))

allocate(zxa(imax))
allocate(zdx(imax))
allocate(zxc(imax))

allocate(zya(jmax))
allocate(zdy(jmax))
allocate(zyc(jmax))

allocate(zza(kmax))
allocate(zdz(kmax))
allocate(zzc(kmax))

!======================================================================
! Set up the number of dimensions

! Check that arrays are large enough for desired number of physical zones
if (max(imax,jmax,kmax)+12 > maxsweep) then
  write(*,*) 'maxsweep too small, recompile with a larger value!'
  stop
endif

! Set the number of dimensions based on array sizes
if (jmax*kmax==1) then
  ndim = 1
else if (kmax==1) then
  ndim = 2
else
  ndim = 3
endif

!======================================================================
! Set up thermal properties

gamm = gam - 1.0

!=======================================================================
! set time and cycle counters

time   = 0.0
timep  = 0.0
timem  = 0.0
ncycle = 0
ncycp  = 0
ncycd  = 0
ncycm  = 0
nfile  = 0
vdtext = 0.

! Set up grid coordinates

call grid(imax,xmin,xmax,zxa,zxc,zdx)
call grid(jmax,ymin,ymax,zya,zyc,zdy)
call grid(kmax,zmin,zmax,zza,zzc,zdz)

if (ndim <= 2) zzc(1) = 0.0
if (ndim == 1) zyc(1) = 0.0

!======================================================================
! Log parameters of problem in history file
! HACK - Comment this out to allow more generic initial conditions
!!$
!!$write (8,"('Oblique Sod shock tube in ',i1,' dimensions.')") ndim
!!$if (ndim==1) then
!!$ write (8,"('Grid dimensions: ',i4)") imax
!!$else if (ndim==2) then
!!$ write (8,"('Grid dimensions: ',i4,' x ',i4)") imax, jmax
!!$else
!!$ write (8,"('Grid dimensions: ',i4,' x ',i4,' x ',i4)") imax, jmax, kmax
!!$endif
!!$write (8,*) 
!!$write (8,"('Ratio of specific heats = ',f5.3)") gam
!!$write (8,"('Pressure ratio is ',f5.3)") pright
!!$write (8,"('Density ratio is ',f5.3)") dright
!!$write (8,*) 

! initialize grid to zero (make density and pressure 1 to prevent errors)
zro = 1d0
zpr = 1d0
zux = 0d0
zuy = 0d0
zuz = 0d0
zfl = 0d0
zgr = 0d0

return
end

!#########################################################################

subroutine grid( nzones, xmin, xmax, xa, xc, dx )

! Create grid to cover physical size from xmin to xmax
! number of physical grid zones is nzones
!
! xa(1) is left boundary location - at xmin
! xa(nzones+1) is right boundary location - at xmax
!----------------------------------------------------------------------
IMPLICIT NONE

! LOCALS
integer :: nzones, n
real(kind=8), dimension(nzones) :: xa, dx, xc
real(kind=8) :: dxfac, xmin, xmax

!=======================================================================

dxfac = (xmax - xmin) / float(nzones)
do n = 1, nzones
  xa(n) = xmin + (n-1)*dxfac
  dx(n) = dxfac
  xc(n) = xa(n) + 0.5*dx(n)
enddo

return
end

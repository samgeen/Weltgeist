subroutine dtcon

! set the timestep using various constraints.
! NB - additional constraints may be required depending on the problem!
!-----------------------------------------------------------------------

! GLOBALS
use global
use zone

IMPLICIT NONE

! LOCALS
INTEGER :: i, j, k
REAL(kind=8) ::  ridt, dtx, dt3, xvel, yvel, zvel
REAL(kind=8)::   widthy, widthz, width

!------------------------------------------------------------------------
!   Hydro constraint on timestep.  Use R*d(theta) if y geometry is angular
ridt = 0.


if(ndim==1) then
  do i = 1, imax
    svel = sqrt(gam*zpr(i,1,1)/zro(i,1,1))/zdx(i)
    xvel = abs(zux(i,1,1)) / zdx(i)
    ridt = max(xvel,ridt,svel)
  enddo
else if(ndim==2) then
  do j = 1, jmax
   do i = 1, imax
     widthy = zdy(j)
     if(ngeomy > 2) widthy = widthy*zxc(i)
     width  = min(zdx(i),widthy)
     svel = sqrt(gam*zpr(i,j,1)/zro(i,j,1))/width
     xvel = abs(zux(i,j,1)) / zdx(i)
     yvel = abs(zuy(i,j,1)) / widthy
     ridt = max(xvel,yvel,svel,ridt)
   enddo
  enddo
else if(ndim==3) then 
  do k = 1, kmax
   do j = 1, jmax
    do i = 1, imax
      widthy = zdy(j)
      widthz = zdz(k)
      if(ngeomy >  2) widthy = widthy*zxc(i)
      if(ngeomz >  2) widthz = widthz*zxc(i)
      if(ngeomz == 5) widthz = widthz*sin(zyc(j))
      width  = min(zdx(i),widthy,widthz)
      svel = sqrt(gam*zpr(i,j,k)/zro(i,j,k))/width
      xvel = abs(zux(i,j,k)) / zdx(i)
      yvel = abs(zuy(i,j,k)) / widthy
      zvel = abs(zuz(i,j,k)) / widthz
      ridt = max(xvel,yvel,zvel,ridt)
    enddo
   enddo
  enddo
endif

ridt = max(ridt,vdtext)
dtx  = courant / ridt     ! global time constraint for given courant parameter

if (dt .gt. 0d0) then
   dt3  = 1.1d0 * dt        ! limiting constraint on rate of increase of dt
   dt   = min( dt3, dtx ) ! use smallest required timestep
else
   dt = dtx
endif
      
if (time/dt > 1.e20) then   ! if timestep becomes too small, stop the program!
  write(*,*) 'Timestep has become too small: dt = ',dt
  write(*,*) '                             time = ',time
  call prin('ABORT')
  stop
endif

return
end

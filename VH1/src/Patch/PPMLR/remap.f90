subroutine remap

! Remap mass, momentum, and energy from the updated lagrangian grid
! to the fixed Eulerian grid, using piecewise parabolic functions.
!-----------------------------------------------------------------------
! GLOBALS
use global
use sweeps

IMPLICIT NONE

! LOCALS
INTEGER :: n, nn, nad
REAL(kind=8), DIMENSION(maxsweep) :: du, ul, u6, dv, vl, v6, dw, wl, w6, de, el, e6
REAL(kind=8), DIMENSION(maxsweep) :: dq, ql, q6, dm, dm0, delta, dvol0
REAL(kind=8), DIMENSION(maxsweep) :: fluxu, fluxv, fluxw, fluxe, fluxq
REAL(kind=8), DIMENSION(maxsweep,nadvect+1) :: dr, r6, rl, fluxr
REAL(kind=8) :: fractn, fractn2, ekin, deltx

REAL(kind=8), PARAMETER :: third  = 1.0 / 3.0
REAL(kind=8), PARAMETER :: twothd = 2.0 / 3.0
REAL(kind=8), PARAMETER :: fourthd= 4.0 / 3.0

!---------------------------------------------------------------------------
! Generate interpolation functions, saving da, al for 
! constructing left and right total energy states. 

call paraset (nmin-1, nmax+1, para, dx, xa)
do nad = 1, nadvect+1
  call parabola(nmin-1, nmax+1, para, r(:,nad), dr(:,nad), r6(:,nad), rl(:,nad), flat)
enddo
call parabola(nmin-1, nmax+1, para, u, du, u6, ul, flat)
call parabola(nmin-1, nmax+1, para, v, dv, v6, vl, flat)
call parabola(nmin-1, nmax+1, para, w, dw, w6, wl, flat)
call parabola(nmin-1, nmax+1, para, q, dq, q6, ql, flat)
call parabola(nmin-1, nmax+1, para, e, de, e6, el, flat)

call volume (nmin, nmax, ngeom, radius, xa0, dx0, dvol0)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Calculate the volume of the overlapping subshells (delta)

select case (ngeom)
 case (0) 
  do n = nmin, nmax+1
    delta(n) = xa(n) - xa0(n)
  enddo
 case (1) 
  do n = nmin, nmax+1
    delta(n) = xa(n) - xa0(n)
    delta(n) = delta(n)*(xa0(n)+.5*delta(n))
  enddo
 case (2) 
  do n = nmin, nmax+1
    delta(n) = xa(n) - xa0(n)
    delta(n) = delta(n)*(xa0(n)*(xa0(n) + delta(n)) + delta(n)**2*third)
  enddo
 case (3) 
  do n = nmin, nmax+1
    delta(n) = (xa(n) - xa0(n)) * radius
  enddo
 case (4) 
  do n = nmin, nmax+1
    delta(n) = (cos(xa0(n)) - cos(xa(n))) * radius
  enddo
 case (5) 
  do n = nmin, nmax+1
    delta(n) = (xa(n) - xa0(n)) * radius 
  enddo
end select

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Calculate the total mass (fluxr), momentum (fluxu), and energy (fluxe)
! in the subshell created by the overlap of the Lagrangian and Eulerican grids.
! If the zone face has moved to the left (deltx > 0), use the integral from the
! left side of zone n (fluxrr).  If the zone face has moved to the right 
! (deltx < 0), use the integral from the right side of zone nn=n-1 (fluxrl).

fluxr = 0.0  ! zero out fluxes to ensure nmin-1 and nmax+2 are zero
fluxu = 0.0
fluxv = 0.0
fluxw = 0.0
fluxe = 0.0
fluxq = 0.0

do n = nmin, nmax + 1
  deltx = xa(n) - xa0(n)
  if(deltx >= 0.0) then
    nn = n - 1
    fractn  = 0.5*deltx/dx(nn)
    fractn2 = 1. - fourthd*fractn
    do nad = 1, nadvect+1
      fluxr (n,nad) = (rl(nn,nad) + dr(nn,nad) - fractn*(dr(nn,nad) - fractn2*r6(nn,nad)))*delta(n)
    enddo
    fluxu (n) = (ul(nn) + du(nn) - fractn*(du(nn) - fractn2*u6(nn)))*fluxr(n,1)
    fluxv (n) = (vl(nn) + dv(nn) - fractn*(dv(nn) - fractn2*v6(nn)))*fluxr(n,1)
    fluxw (n) = (wl(nn) + dw(nn) - fractn*(dw(nn) - fractn2*w6(nn)))*fluxr(n,1)
    fluxe (n) = (el(nn) + de(nn) - fractn*(de(nn) - fractn2*e6(nn)))*fluxr(n,1)
    fluxq (n) = (ql(nn) + dq(nn) - fractn*(dq(nn) - fractn2*q6(nn)))*fluxr(n,1)
  else
    fractn   = 0.5*deltx/dx(n)
    fractn2  = 1. + fourthd*fractn
    do nad = 1, nadvect+1
      fluxr(n,nad) = (rl(n,nad) - fractn*(dr(n,nad) + fractn2*r6(n,nad)))*delta(n)
    enddo
    fluxu(n) = (ul(n) - fractn*(du(n) + fractn2*u6(n)))*fluxr(n,1)
    fluxv(n) = (vl(n) - fractn*(dv(n) + fractn2*v6(n)))*fluxr(n,1)
    fluxw(n) = (wl(n) - fractn*(dw(n) + fractn2*w6(n)))*fluxr(n,1)
    fluxe(n) = (el(n) - fractn*(de(n) + fractn2*e6(n)))*fluxr(n,1)
    fluxq(n) = (ql(n) - fractn*(dq(n) + fractn2*q6(n)))*fluxr(n,1)
  endif
enddo

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Advect mass, momentum, and energy by moving the subshell quantities 
! into the appropriate Eulerian zone. 

do n = nmin-1, nmax+1  ! must update nmin-1, nmax+1 for possible second remap
  do nad = 1,nadvect+1 ! treat all advected quantities including total mass density
    dm (n) = r(n,nad) * dvol(n)
    dm0(n) = (dm(n) + fluxr(n,nad) - fluxr(n+1,nad))
    r  (n,nad) = dm0(n)/dvol0(n)
    r  (n,nad) = max(smallr,r(n,nad))
  enddo
  dm0(n) = 1./(r(n,1)*dvol0(n))
  u  (n) = (u(n)*dm(n) + fluxu(n)-fluxu(n+1))*dm0(n)
  v  (n) = (v(n)*dm(n) + fluxv(n)-fluxv(n+1))*dm0(n)
  w  (n) = (w(n)*dm(n) + fluxw(n)-fluxw(n+1))*dm0(n)
  e  (n) = (e(n)*dm(n) + fluxe(n)-fluxe(n+1))*dm0(n)
  q  (n) = (q(n)*dm(n) + fluxq(n)-fluxq(n+1))*dm0(n)
enddo
         
! If flow is highly supersonic remap on internal energy, else on total E
do n = nmin, nmax
  ekin = 0.5*(u(n)**2+v(n)**2+w(n)**2)
  if(ekin/q(n) < 100.0) q(n) = e(n) - ekin
  p(n) = gamm*r(n,1)*q(n)
  p(n) = max(smallp,p(n))
enddo

return
end


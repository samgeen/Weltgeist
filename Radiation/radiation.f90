! Very simple radiation tracing module
! Assumes all neutral hydrogen absorbs all photons
! Sam Geen, March 2018

subroutine trace_radiation(Sphotons,x,nH,xHII,T,ncell)
! Trace a ray through the spherical grid and ionise everything in the way
! Use very simple instant ionisation/recombination model
! Sphotons = d(nphotons)/dt = integrate(4*pi*r^2*nH(r)^2*alpha_B).dr
implicit none

! Interface variables
integer::ncell
real(kind=8)::Sphotons
real(kind=8),dimension(1:ncell)::x,nH,xHII
! Internal variables
integer::icell


! Loop through 


end subroutine trace_radiation
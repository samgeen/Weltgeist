! Very simple radiation tracing module
! Assumes all neutral hydrogen absorbs all photons with a simple dust model
! Sam Geen, March 2018

subroutine trace_radiation(dr,Qion,sigmaDust,nH,drecombinationsdr,ncell)
! Trace a ray through the spherical grid and ionise everything in the way
! Use very simple instant ionisation/recombination model
! This part is largely to make iterating through the ray faster, other parts done in numpy in radiation.py
implicit none

! Interface variables
real(kind=8)::dr
integer::ncell
real(kind=8),dimension(1:ncell)::Qion,sigmaDust,nH,drecombinationsdr
! Internal variables
integer::icell

! Loop through cells
do icell=2,ncell
    ! See Draine (2011) equation 2
    Qion(icell) = Qion(icell-1) - (drecombinationsdr(icell) + nH(icell) * sigmaDust(icell) * Qion(icell-1))*dr
    ! If no photons left, set the remaining bins to zero and exit the loop
    if (Qion(icell) .lt. 0.0) then
        Qion(icell:ncell) = 0.0
        exit
    endif 
end do


! Loop through 


end subroutine trace_radiation
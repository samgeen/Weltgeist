! Wrapper around singlestar functions that makes it easier to call from f2py
! Sam Geen, January 2018

SUBROUTINE star_setup(tableloc)
use singlestar_module
implicit none

character(len=200), intent(in)::tableloc
! f2py character, intent(in),dimension(200)::tableloc

! For simplicity, assume the table is in this folder
call ssm_setup(tableloc)

END SUBROUTINE star_setup

SUBROUTINE star_lifetime(mass_ini,lifetime)
  ! Return the lifetime of a star in seconds
  ! mass_ini - initial stellar mass in Msun
  ! RETURNS
  ! lifetime - lifetime of the star in seconds
  use singlestar_module
  implicit none

  integer,parameter::dp=kind(1.0D0) ! real*8  

  real(dp),intent(in)::mass_ini
  real(dp),intent(out)::lifetime
  ! f2py real(dp),intent(in)::mass_ini
  ! f2py real(dp),intent(out)::lifetime
  ! Interpolate over table to find lifetime for this mass
  call ssm_lifetime(mass_ini,lifetime)
  
END SUBROUTINE star_lifetime


SUBROUTINE star_winds(mass_ini,t,dt,energy,massloss)
  ! Wind energy, mass loss between t and dt (from birth of star)
  ! mass_ini - initial stellar mass in Msun
  ! t - age of star in seconds
  ! dt - timestep length in seconds
  ! RETURNS
  ! energy - energy emitted in ergs between t and dt
  ! massloss - mass lost in g between t and dt
  use singlestar_module
  implicit none

  integer,parameter::dp=kind(1.0D0) ! real*8 

  real(dp),intent(in)::mass_ini
  real(dp),intent(in)::t
  real(dp),intent(in)::dt
  real(dp),intent(out)::energy
  real(dp),intent(out)::massloss
  ! f2py real(dp),intent(in)::mass_ini
  ! f2py real(dp),intent(in)::t
  ! f2py real(dp),intent(in)::dt
  ! f2py real(dp),intent(out)::energy
  ! f2py real(dp),intent(out)::massloss

  call ssm_winds(mass_ini,t,dt,energy,massloss)

END SUBROUTINE star_winds


SUBROUTINE star_radiation(mass_ini,t,dt,nphotons)
  ! Number of photons to emit in each group from a star
  ! mass_ini - initial stellar mass in Msun
  ! t - age of star in seconds
  ! dt - timestep length in seconds
  ! RETURNS
  ! nphotons - number of photons emitted between t and t+dt
  use singlestar_module
  implicit none

  integer,parameter::dp=kind(1.0D0) ! real*8  
  integer,parameter::ngroups=5

  real(dp),intent(in)::mass_ini
  real(dp),intent(in)::t
  real(dp),intent(in)::dt
  real(dp),dimension(ngroups),intent(out)::nphotons
  ! f2py real(dp),intent(in)::mass_ini
  ! f2py real(dp),intent(in)::t
  ! f2py real(dp),intent(in)::dt
  ! f2py real(dp),dimension(ngroups),intent(out)::nphotons

  call ssm_radiation(mass_ini,t,dt,nphotons)

END SUBROUTINE star_radiation

SUBROUTINE star_bandenergies(mass_ini,t,dt,energies)
  ! Energies in photon bands emitted between t and t+dt
  ! mass_ini - initial stellar mass in Msun
  ! t - age of star in seconds
  ! dt - timestep length in seconds
  ! RETURNS
  ! energies - photon energy emitted in this band between t and t+dt
  use singlestar_module
  implicit none

  integer,parameter::dp=kind(1.0D0) ! real*8  
  integer,parameter::nbands=4

  real(dp),intent(in)::mass_ini
  real(dp),intent(in)::t
  real(dp),intent(in)::dt
  real(dp),dimension(nbands),intent(out)::energies
  ! f2py real(dp),intent(in)::mass_ini
  ! f2py real(dp),intent(in)::t
  ! f2py real(dp),intent(in)::dt
  ! f2py real(dp),dimension(nbands),intent(out)::energies

  call ssm_bandenergies(mass_ini,t,dt,energies)

END SUBROUTINE star_bandenergies

SUBROUTINE star_effectivetemperature(mass_ini,t,Teff)
  ! Stellar effective temperature
  ! mass_ini - initial stellar mass in Msun
  ! t - age of star in seconds
  ! RETURNS
  ! Teff - effective temperature of star in Kelvin at time t
  use singlestar_module
  implicit none

  integer,parameter::dp=kind(1.0D0) ! real*8  

  real(dp),intent(in)::mass_ini
  real(dp),intent(in)::t
  real(dp),intent(out)::Teff
  ! f2py real(dp),intent(in)::mass_ini
  ! f2py real(dp),intent(in)::t
  ! f2py real(dp),intent(out)::Teff
  
  call ssm_Teffective(mass_ini,t,Teff)

END SUBROUTINE star_effectivetemperature


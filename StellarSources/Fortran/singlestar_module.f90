! Output single star properties from Starburst99/Geneva tables
! Sam Geen, April 2017

MODULE singlestar_module
  
  use table_1d_module
  
  implicit none

  logical::ssm_is_setup=.false.

  integer,parameter::dp=kind(1.0D0) ! real*8  

  ! Number of photon groups the spectrum is divided into
  ! HACK TO MAKE THIS WORK WITHOUT RADIATION MODULE
  integer,parameter::ngroups=5
  ! Additional bands tracking photon luminosities
  integer,parameter::nbands=4
  
  private
  
  ! Public functions (i.e. useable outside module) are:
  ! ssm_setup(tableloc_in)
  ! ssm_lifetime(mass_ini,lifetime)
  ! ssm_winds(mass_ini,t,dt,energy,massloss)
  ! ssm_radiation(mass_ini,t,dt,nphotons)
  public::ssm_setup, ssm_lifetime, ssm_winds, ssm_radiation, ssm_bandenergies, ssm_Teffective, ssm_supernovae

  ! Private functions (used only inside the module) are:
  private::ssm_mtoi, ssm_itom, ssm_filename, ssm_interpolate, ssm_title
  
  ! Number of tables to read
  integer::ssm_numtables
  ! Location of tables to read
  character(len=200)::ssm_tableloc
  ! Lifetimes of sources by mass
  type(lookup_table)::ssm_lifetimes
  ! Supernovae (by mass)
  type(lookup_table)::ssm_sn_masslosses
  type(lookup_table)::ssm_sn_energies
  type(lookup_table)::ssm_sn_yields
  ! Winds
  type(lookup_table),dimension(:), allocatable::ssm_energies
  type(lookup_table),dimension(:), allocatable::ssm_masslosses
  ! Radiation
  type(lookup_table),dimension(:,:), allocatable::ssm_rads
  type(lookup_table),dimension(:,:), allocatable::ssm_bands
  type(lookup_table),dimension(:), allocatable::ssm_Teff
CONTAINS
  
! PUBLIC / INTERFACE FUNCTIONS
  
SUBROUTINE ssm_setup(tableloc_in)
  ! Set up the Single Star Module
  ! IMPORTANT! MUST BE RUN BEFORE YOU USE THE MOULE
  ! tableloc_in - string giving location of tables to read

  character(len=*),intent(in)::tableloc_in
  character(len=200)::filename
  integer::it,ip
  ssm_tableloc = tableloc_in
  ! Read the lifetimes tables
  ! NOTE - this is a lifetime per star
  ! These stars MUST match the stars in the other tables
  ! We assume stars are 5,10,15,...115,120 Msun
  if (ssm_is_setup) return
  filename = TRIM(ssm_tableloc)//"lifetimes.dat"
  call setup_table(ssm_lifetimes,filename)
  ssm_numtables = ssm_lifetimes%size

  ! Read the supernovae tables
  filename = TRIM(ssm_tableloc)//"sn_masslosses.dat"
  call setup_table(ssm_sn_masslosses,filename)
  filename = TRIM(ssm_tableloc)//"sn_energies.dat"
  call setup_table(ssm_sn_energies,filename)
  filename = TRIM(ssm_tableloc)//"sn_yields.dat"
  call setup_table(ssm_sn_yields,filename)

  ! Allocate wind and radiation tables
  allocate(ssm_energies(ssm_numtables))
  allocate(ssm_masslosses(ssm_numtables))
  allocate(ssm_rads(ssm_numtables,ngroups))
  allocate(ssm_bands(ssm_numtables,nbands))
  allocate(ssm_Teff(ssm_numtables))
  ! Run through tables and set them up
  do it=1,ssm_numtables
     ! Winds
     call ssm_filename(it,"cumulenergy",filename)
     call setup_table(ssm_energies(it),filename)
     call ssm_filename(it,"cumulmassloss",filename)
     call setup_table(ssm_masslosses(it),filename)
     ! Radiation spectrum
     ! TODO: Make more generic for, e.g., 5 groups?
     ip = 0
     if (ngroups.gt.3) then
        ip = 2
        call ssm_filename(it,"cumulIR",filename)
        call setup_table(ssm_rads(it,1),filename)
        call ssm_filename(it,"cumulOpt",filename)
        call setup_table(ssm_rads(it,2),filename)
     endif
     call ssm_filename(it,"cumulHII",filename)
     call setup_table(ssm_rads(it,1+ip),filename)
     call ssm_filename(it,"cumulHeII",filename)
     call setup_table(ssm_rads(it,2+ip),filename)
     call ssm_filename(it,"cumulHeIII",filename)
     call setup_table(ssm_rads(it,3+ip),filename)

     ! Radiation energy bands
     call ssm_filename(it,"cumulEKband",filename)
     call setup_table(ssm_bands(it,1),filename)
     call ssm_filename(it,"cumulEVband",filename)
     call setup_table(ssm_bands(it,2),filename)
     call ssm_filename(it,"cumulLbol",filename)
     call setup_table(ssm_bands(it,3),filename)
     call ssm_filename(it,"cumulEion",filename)
     call setup_table(ssm_bands(it,4),filename)
     ! Stellar surface temperature
     call ssm_filename(it,"Teff",filename)
     call setup_table(ssm_Teff(it),filename)
  enddo
  ssm_is_setup=.true.
  ! TODO: Should I also deallocate later???
  
END SUBROUTINE ssm_setup

SUBROUTINE ssm_lifetime(mass_ini,lifetime)
  ! Return the lifetime of a star in seconds
  ! mass_ini - initial stellar mass in Msun
  ! RETURNS
  ! lifetime - lifetime of the star in seconds

  real(dp),intent(in)::mass_ini
  real(dp),intent(out)::lifetime
  ! Interpolate over table to find lifetime for this mass
  call find_value(ssm_lifetimes,mass_ini,lifetime)
  
END SUBROUTINE ssm_lifetime

SUBROUTINE ssm_supernovae(mass_ini,energy,massloss,yield)
  ! Return the lifetime of a star in seconds
  ! mass_ini - initial stellar mass in Msun
  ! RETURNS
  ! energy - energy of the supernova in ergs
  ! massloss - mass injected by the supernova in g
  ! yield - absolute metallicity of ejecta (solar = 0.014)

  real(dp),intent(in)::mass_ini
  real(dp),intent(out)::energy
  real(dp),intent(out)::massloss
  real(dp),intent(out)::yield
  ! Interpolate over table to find energy, mass loss and yield for this mass
  call find_value(ssm_sn_energies,mass_ini,energy)
  call find_value(ssm_sn_masslosses,mass_ini,massloss)
  call find_value(ssm_sn_yields,mass_ini,yield)
  
END SUBROUTINE ssm_supernovae

SUBROUTINE ssm_winds(mass_ini,t,dt,energy,massloss)
  ! Wind energy, mass loss between t and dt (from birth of star)
  ! mass_ini - initial stellar mass in Msun
  ! t - age of star in seconds
  ! dt - timestep length in seconds
  ! RETURNS
  ! energy - energy emitted in ergs between t and dt
  ! massloss - mass lost in g between t and dt

  real(dp),intent(in)::mass_ini
  real(dp),intent(in)::t
  real(dp),intent(in)::dt
  real(dp),intent(out)::energy
  real(dp),intent(out)::massloss
  real(dp)::v1,v2
  ! Get energies and mass losses for t and dt
  call ssm_interpolate(ssm_energies, mass_ini, t, v1)
  call ssm_interpolate(ssm_energies, mass_ini, t+dt, v2)
  energy = v2 - v1
  call ssm_interpolate(ssm_masslosses, mass_ini, t, v1)
  call ssm_interpolate(ssm_masslosses, mass_ini, t+dt, v2)
  massloss = v2 - v1

END SUBROUTINE ssm_winds

SUBROUTINE ssm_radiation(mass_ini,t,dt,nphotons)
  ! Number of photons to emit in each group from a star
  ! mass_ini - initial stellar mass in Msun
  ! t - age of star in seconds
  ! dt - timestep length in seconds
  ! RETURNS
  ! nphotons - number of photons emitted between t and t+dt

  real(dp),intent(in)::mass_ini
  real(dp),intent(in)::t
  real(dp),intent(in)::dt
  real(dp),dimension(ngroups),intent(out)::nphotons
  real(dp)::v1,v2
  integer::ig

  do ig=1,ngroups
     call ssm_interpolate(ssm_rads(:,ig),mass_ini,t,v1)
     call ssm_interpolate(ssm_rads(:,ig),mass_ini,t+dt,v2)
     nphotons(ig) = v2 - v1
  enddo

END SUBROUTINE ssm_radiation

SUBROUTINE ssm_bandenergies(mass_ini,t,dt,energies)
  ! Energies in photon bands emitted between t and t+dt
  ! mass_ini - initial stellar mass in Msun
  ! t - age of star in seconds
  ! dt - timestep length in seconds
  ! RETURNS
  ! energies - photon energy emitted in this band between t and t+dt

  real(dp),intent(in)::mass_ini
  real(dp),intent(in)::t
  real(dp),intent(in)::dt
  real(dp),dimension(nbands),intent(out)::energies
  real(dp)::v1,v2
  integer::ib

  do ib=1,nbands
     call ssm_interpolate(ssm_bands(:,ib),mass_ini,t,v1)
     call ssm_interpolate(ssm_bands(:,ib),mass_ini,t+dt,v2)
     energies(ib) = v2 - v1
  enddo

END SUBROUTINE ssm_bandenergies

SUBROUTINE ssm_Teffective(mass_ini,t,Teff)
  ! Ionised gas temperature
  ! mass_ini - initial stellar mass in Msun
  ! t - age of star in seconds
  ! RETURNS
  ! Teff - effective temperature of star in Kelvin at time t

  real(dp),intent(in)::mass_ini
  real(dp),intent(in)::t
  real(dp),intent(out)::Teff
  real(dp)::v1
  ! Get energies and mass losses for t and dt
  call ssm_interpolate(ssm_Teff, mass_ini, t, v1)
  Teff = v1

END SUBROUTINE ssm_Teffective

! PRIVATE FUNCTIONS
! CANNOT BE CALLED OUTSIDE THE MODULE

SUBROUTINE ssm_mtoi(mass_ini,index)
  ! Convert mass to index (1 index every 5 Msun)
  
  real(dp),intent(in)::mass_ini
  integer,intent(out)::index

  index = int(mass_ini/5d0)
  
END SUBROUTINE ssm_mtoi

SUBROUTINE ssm_itom(index,mass_ini)
  ! Convert mass to index (1 index every 5 Msun)

  integer,intent(in)::index
  real(dp),intent(out)::mass_ini

  mass_ini = 5d0*index
  
END SUBROUTINE ssm_itom

SUBROUTINE ssm_filename(index,prop,filename)
  ! Returns a filename for each of the tables

  integer,intent(in)::index
  character(len=*),intent(in)::prop
  character(len=*),intent(out)::filename
  character(LEN=3)::nmass
  ! Index is 1 to 24 for masses 5 to 120
  call ssm_title(index*5,nmass)
  ! Put together filename
  filename = TRIM(ssm_tableloc)//"_m"//TRIM(nmass)// &
       & prop//".dat"
  ! TODO: Add error checking for file existence
END SUBROUTINE ssm_filename

SUBROUTINE ssm_interpolate(tables, mass_ini, time, output)
  ! Find interpolated values between two masses
  ! This works by assuming that each stellar track is "stretched"
  !  so that a point at x=[0,1]*lifetime is similar between tracks
  type(lookup_table),dimension(:),intent(in)::tables
  real(dp),intent(in)::mass_ini,time
  real(dp),intent(out)::output
  integer::i1,i2
  real(dp)::m1,m2,x1,x2,x,life,life1,life2,o1,o2
  type(lookup_table)::table1,table2
  ! Find tables to interpolate between
  call ssm_mtoi(mass_ini,i1)
  i2 = i1+1
  ! Table overflow: return value from largest star
  if (i2.gt.ssm_numtables) then
     table2 = tables(ssm_numtables)
     call find_value(table2,time,output)
  ! Underflow: assume small stars do nothing
  else if (i1.lt.1) then
     output = 0d0
  ! Interpolate
  else
     ! Get masses and lifetimes of bounding stellar tracks
     call ssm_itom(i1,m1)
     call ssm_itom(i2,m2)
     call ssm_lifetime(mass_ini,life)
     call ssm_lifetime(m1,life1)
     call ssm_lifetime(m2,life2)
     table1 = tables(i1)
     table2 = tables(i2)
     ! Scale time to a scale-free value
     x = time / life
     ! Find value in each track scaled to their lifetimes
     call find_value(table1,x*life1,o1)
     call find_value(table2,x*life2,o2)
     ! Interpolate between these values in mass
     !write(*,*) "OUTPUT", o1, o2, m1, m2, mass_ini, life1, life2, x
     output = (o2 * (mass_ini - m1) + o1 * (m2 - mass_ini)) / (m2 - m1)
  end if
  
END SUBROUTINE ssm_interpolate

SUBROUTINE ssm_title(n,nchar)
  ! Returns integer as a 3-character zero-filled string
  ! Used to generate filenames for tables

  implicit none
  integer,intent(in)::n
  character(LEN=3),intent(out)::nchar

  character(LEN=1)::nchar1
  character(LEN=2)::nchar2
  character(LEN=3)::nchar3

  if(n.ge.100)then
     write(nchar3,'(i3)') n
     nchar = nchar3
  elseif(n.ge.10)then
     write(nchar2,'(i2)') n
     nchar = '0'//nchar2
  else
     write(nchar1,'(i1)') n
     nchar = '00'//nchar1
  endif

END SUBROUTINE ssm_title

END MODULE singlestar_module

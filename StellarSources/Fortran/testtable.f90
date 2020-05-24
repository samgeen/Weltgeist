program testtable
  use table_1d_module
  use singlestar_module
  implicit none
  integer,parameter::dp=kind(1.0D0) ! real*8  
  real(dp),parameter::Myrins = 3.1556926d13
  integer::i
  real(dp)::e1,e2,mstar,t,dt
  real(dp),dimension(5)::nphotons

  ! Test the basic table setup
  type(lookup_table)::energy_table
  call setup_table(energy_table, "../../Compressed/singlestar_z0.014_m020cumulenergy.dat")
  call find_value(energy_table,Myrins,e1)
  call find_value(energy_table,1.00001*Myrins,e2)
  write(*,*) "Test reading tables"
  write(*,*) "ENERGY AT 1 MYR, 1.00001 MYR", e1, e2
  write(*,*) "DIFFERENCE:", e2 - e1
  call find_value(energy_table,100*Myrins,e1)
  call find_value(energy_table,200*Myrins,e2)
  write(*,*) "ENERGY AT 100 MYR, 200 MYR (SHOULD BE THE SAME)", e1, e2
  call find_value(energy_table,0d0,e1)
  call find_value(energy_table,-100*Myrins,e2)
  write(*,*) "ENERGY AT 0 MYR, -100 MYR (SHOULD BE 0)", e1, e2
  call find_value(energy_table,1d-10*Myrins,e1)
  call find_value(energy_table,1d-7*Myrins,e2)
  write(*,*) "ENERGY AT 1e-10,1e-7 MYR", e1, e2

  ! Test single star module
  write(*,*)"-------------------"
  mstar = 8.29
  write(*,*)"Testing single star module with star mass", mstar, "Msun"
  !debug_lkup = .true.
  !call ssm_setup("../../Outputs/singlestar_z0.0140")
  call ssm_setup("../../Compressed/singlestar_z0.014")
  call ssm_lifetime(mstar,e1)
  write(*,*) "Lifetime:", e1/Myrins
  t = 2d0*Myrins
  dt = 1d-2*Myrins
  call ssm_winds(mstar,t,dt,e1,e2)
  write(*,*) "Energy, mass loss at 2 Myr for 0.01 Myr", e1, e2
  call ssm_radiation(mstar,t,dt,nphotons)
  write(*,*) "Number of photons in each group for 2+0.01Myr", nphotons
  ! Output a sample of tracks
  !do i=5,120,5
  !   call output_track(real(i,dp))
  !end do
  !call output_track(32.5d0)
  call output_track(5d0)
  !call output_track(8.29d0)
  !call output_track(9.99d0)
  !call output_track(62.5d0)
  !call output_track(1d0)
  !call output_track(200d0)

end program testtable


SUBROUTINE Replace_Text (s,text,rep,out)
implicit none
CHARACTER(*):: text,rep
CHARACTER(1)::c
CHARACTER(20):: s, out
INTEGER:: i, nt, nr

out = s
DO i=1,LEN(s)
   c = s(i:i)
   if (c.eq.text) then
      out(i:i) = rep
   endif
END DO
END SUBROUTINE Replace_Text

SUBROUTINE output_track(mstar)

  use singlestar_module
  implicit none
  integer,parameter::dp=kind(1.0D0) ! real*8  

  real(dp),intent(in)::mstar
  real(dp)::life,t,dt,e,ml
  real(dp),dimension(5)::np
  integer::i
  integer,parameter::ilun=101
  integer,parameter::nt=1000
  character(len=20)::mstr,tmp
  character(len=200)::filename
  real(dp),dimension(nt)::ts,es,mls,np1s,np2s,np3s,np4s,np5s
  ! Set up file
  write(mstr,'(F10.3)') mstar
  call replace_text(mstr," ","0",tmp)
  call replace_text(tmp,".","p",mstr)
  filename="/home/samgeen/Data/Simulations/StellarSources/test/track"//TRIM(mstr)//".dat"
  open(ilun,file=filename,form="unformatted",status='REPLACE')
  ! Find lifetime (gives number of points)
  call ssm_lifetime(mstar,life)
  dt = life/real(nt,dp)
  ! Read data
  write(*,*) mstar, life
  do i=1,nt
     t = (i*life)/real(nt,dp)
     call ssm_winds(mstar,t,dt,e,ml)
     call ssm_radiation(mstar,t,dt,np)
     ts(i) = t
     es(i) = e
     mls(i) = ml
     np1s(i) = np(1)
     np2s(i) = np(2)
     np3s(i) = np(3)
     np4s(i) = np(4)
     np5s(i) = np(5)
  enddo
  ! Make cumulative values
  !do i=2,nt
  !   es(i) = es(i)+es(i-1)
  !   mls(i) = mls(i)+mls(i-1)
  !   np1s(i) = np1s(i)+np1s(i-1)
  !   np2s(i) = np2s(i)+np2s(i-1)
  !   np3s(i) = np3s(i)+np3s(i-1)
  !enddo
  ! Write
  write(ilun) ts
  write(ilun) es
  write(ilun) mls
  write(ilun) np1s
  write(ilun) np2s
  write(ilun) np3s
  write(ilun) np4s
  write(ilun) np5s
  close(ilun)
  
END SUBROUTINE output_track

module cooling_module

! Define double precision as in RAMSES
integer,parameter::dp=kind(1.0D0) ! real*8

! Define gamma as in RAMSES
real(dp)::gamma=1.4d0

! Define mu as in Patrick Hennebelle's original module
real(dp)::mu=1.4d0

! HACK - set rt by hand for now to keep RAMSES-like interface
logical,parameter::rt=.false.

! HARD-CODED PARAMETERS FOR SAFE BOUNDING OF PHYSICAL VALUES
real(dp),parameter::smallr=1d-10

contains

!=======================================================================
subroutine solve_cooling(nH,T2,xHII,zsolar,dt,deltaT2,ncell)
!=======================================================================
  implicit none
  ! Input/output variables to this function
  ! nH - hydrogen number density in PHYSICAL units
  ! T2 - temperature / mu in PHYSICAL units
  ! zsolar - Metallicity in solar units (Zphys / 0.02)
  ! dt - cooling timestep in seconds
  ! deltaT2 - temperature change in K
  ! ncell - number of elements in the vector
  integer::ncell
  real(kind=8)::dt
  real(kind=8),dimension(1:ncell)::nH,T2,deltaT2,zsolar,xHII
  ! Temporary variables
  integer::i
  real(kind=8)::TT_ini, mu
  ! Units
  ! TODO - FIGURE THIS UNITS STUFF OUT
  !real(kind=8) :: scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2
  !call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)
  do i=1,ncell
     TT_ini = T2(i)
     call calc_temp(nH(i),T2(i),xHII(i),dt)
     deltaT2(i) = (T2(i) - TT_ini)
  end do
end subroutine solve_cooling

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine  calc_temp(NN,TT,XX,dt_tot)
    implicit none

    integer :: n,i,j,k,idim, iter, itermax,ii

    real(dp) :: dt, dt_tot, temps, dt_max, itermoy
    real(dp) :: rho,temp,dt_tot_unicode

    real(dp) :: mm,uma, kb, alpha,kb_mm
    real(dp) :: NN,TT, TTold, ref,dRefdT, eps, vardt,varrel, dTemp, XX
    real(dp) :: rhoutot2
    !real(dp) :: scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2
    !
    ! Cette routine fonctionne en cgs
    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    kb  =  1.38062d-16   ! erg/degre
    !  uma =  1.660531e-24  ! gramme
    !  mu  =  1.4
    !  mm = mu*uma
    !  kb_mm = kb / mm
    !  TT = TT  / kb  !/ kb_mm

    !call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)
    !scale_T2 = scale_T2 * mu

    if( TT .le. 0.) then
        TT = 50!. / scale_T2
        return
    endif

    !if( TT*scale_T2 .gt. 50.) then
    !TT = 50. / scale_T2
    !return
    !endif

    vardt = 10.**(1./10.); varrel = 0.2

    !  nn = (rho/(gramme/cm3)) /mm

    itermax = 0 ; itermoy = 0.

    if (NN .le. smallr) then
        if( NN .le. 0)  write(*,*) 'prob dens',NN
        NN = smallr
    endif

    alpha = NN*kb/(gamma-1.)

    iter  = 0 ; temps = 0.
    do while ( temps < dt_tot)
        if (TT .lt.0) then
            write(*,*) 'prob Temp',TT, NN
            !         write(*,*) 'repair assuming isobariticity'
            NN = max(NN,smallr)
            TT = min(4000./NN,8000.)
        endif


        TTold = TT

        ! Calculate cooling rates
        call calc_cooling_rates(TT,NN,XX,ref,dRefdT)

        ! Calculate for first iteration
        if (iter == 0) then
            if (dRefDT .ne. 0.) then
                dt = abs(1.0E-1 * alpha/dRefDT)
            else
                dt = 1.0E-1 * dt_tot
            endif
            dt_max = dt_tot - temps
            if (dt > 0.7*dt_max) dt = dt_max*(1.+1.0E-12)
        endif

        ! Calculate change in temperature
        dTemp = ref/(alpha/dt - dRefdT)

        eps = abs(dTemp/TT)
        if (eps > 0.2) dTemp = 0.2*TTold*dTemp/abs(dTemp)

        TT = TTold + dTemp
        if (TT < 0.) then
            write(*,*) 'Temperature negative !!!'
            write(*,*) 'TTold,TT   = ',TTold,TT
            write(*,*) 'rho   = ',rho
            TT = 100.  !*kelvin
        endif


        ! Set up for next iteration
        iter = iter + 1

        temps = temps + dt

        dt = vardt*varrel*dt/Max(vardt*eps, varrel)

        dt_max = dt_tot - temps
        if (dt > 0.7*dt_max) dt = dt_max*(1.+1.0E-12)
    enddo

    return
end subroutine calc_temp

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


subroutine calc_cooling_rates(T,n,xin,ref,drefdT)
  implicit none

  real(dp) :: T,n,xin,ref,T2,ref2,drefdT
  real(dp),parameter::eps=1.0d5

!------------------------------------------------------

  T2 = T*(1.+eps)
  ! Get cooling rates for T and T+dT
  call cool_rate(T, n,xin,ref )
  call cool_rate(T2,n,xin,ref2)

  ! Get local gradient in cooling rate
  dRefdT = (ref2-ref)/(T2-T)

  return

end subroutine calc_cooling_rates

subroutine cool_rate(T,n,xin,ref)

  implicit none

  real(dp) :: T,n,xin,ref

!------------------------------------------------------

  ! Get single rate
  if (T < 10035.d0) then
     call cool_rate_lowT(T,n,ref)
  else
     call cool_rate_highT(T,n,ref)
  end if

  ! TODO: include photoheating/cooling of metals !!!

end subroutine cool_rate

subroutine cool_rate_highT(T,n,ref)
    implicit none

    real(dp) :: T,P,N,x,ne                         !x est le taux d'ionisation
    real(dp) :: T2, ref2,dRefDT
    real(dp) :: froid,chaud,ref,nc, froidc
    real(dp) :: froidcII, froido, froidh
    real(dp) :: froidc_m,froidcII_m,froido_m
    real(dp) :: param, G0, epsilon,k1,k2,bet,froidrec
    real(dp) :: eps

    real(dp) :: logT, intcst, logT2
    real(dp) :: ion, neut
    real(dp) :: refion

    !taux de refroidissement base sur Dopita et Sutherland

    logT=log10(T)

    if (logT .LT. 4.0) then
        froid=0.1343*logT**3-1.3906*logT**2+5.1554*logT-31.967
    else if (logT .LT. 4.25) then
        froid=12.64*logT-75.56
    else if (logT .LT. 4.35) then
        froid=-0.3*logT-20.565
    else if (logT .LT. 4.9) then
        froid=1.745*logT-29.463
    else if (logT .LT. 5.4) then
        froid=-20.9125
    else if (logT .LT. 5.9) then
        froid=-1.795*logT-11.219
    else if (logT .LT. 6.2) then
        froid=-21.8095
    else if (logT .LT. 6.7) then
        froid=-1.261*logT-13.991
    else
        froid=-22.44
    endif

    froid=-1.0*10.0**(froid)
    ref = froid*(n**2)


end subroutine cool_rate_highT


subroutine cool_rate_lowT(T,n,ref)


  implicit none

  real(dp) :: T,P,N,x,ne                      !x est le taux d'ionisation
  real(dp) :: dRefDT, x_ana
  real(dp) :: froid,chaud,ref,nc, froidc
  real(dp) :: froidcII, froido, froidh
  real(dp) :: froidc_m,froidcII_m,froido_m
  real(dp) :: param, G0, epsilon,k1,k2,bet,froidrec


  !fonction de chauffage et de refroidissement calculee a partir des
  !refroidissements des differents elements



  !abondance de carbone 3.5 10-4, depletion 0.4

!!! on calcule l'ionisation
!!! on suppose que si x est superieure a 1.d-4 elle est domine par
!!! l'hydrogene et que la valeur inf est donne par le carbone
!!! et vaut 3.5 1.d-4 * depletion * densite

!!! Pour les electrons dus a l'hydrogene on prend la
!!! formule donnee par Wolfire et al. 2003 appendice C2.
!!! on prend un taux d'ionisation de 1.d-16 G0'=GO/1.7
!!! Z'd = 1 et phipah=0.5



  ne = 2.4E-3*((T/100.)**0.25)/0.5 !formule C15 Wolfire et al. 2003

  ! Analytic ionisation in absence of photoionisation
  x_ana = ne / N   ! ionisation
  x_ana = min(x_ana,0.1)
  x_ana = max(x_ana,3.5E-4*0.4)

  !transition hyperfine a basse temperature: carbone et oxygene
  !chiffre pris dans la these de Karl Joulain

  !refroidissement par le carbone




  !      froidcII = ( 2.2d-23                     !excitation par H
  !     c          + 5.5d-20 * 2 / sqrt(T) * x ) !excitation par e
  !     c              * 3.5d-4 * 0.4d0 * exp(-92.d0 / T)


 ! NOTE - HERE WE USE THE NON-PHOTOIONISED RATES AS THIS MIGHT 
 !        BE TOO HIGH AT x=1
  froidcII =  92. * 1.38E-16 * 2. * (2.8E-7* ((T/100.)**(-0.5))*x_ana + 8.E-10*((T/100.)**(0.07))) &
       * 3.5E-4 * 0.4 * exp(-92./ T)
  !     c               3.d-4  * exp(-92. / T)


  !refroidissement par l'oxygene
  !abondance 8.6 10-4 depletion 0.8

  froido = 1.E-26 * sqrt(T) * (24. * exp(-228./ T) + 7. * exp(-326./ T) )


  !      froido =  230.d0*1.38d-16 * (
  !     c            1.4d-8*x + 9.2d-11 *(T /100.d0)**(0.67) )
  !     c     * exp(-230.d0 / T)

  !      froido = froido +  330.d0*1.38d-16 *(
  !     c            1.4d-8*x + 4.3d-11 *(T /100.d0)**(0.8) )
  !     c     * exp(-330.d0 / T)

  !      froido = froido +  98.d0*1.38d-16 * (
  !     c            5.d-9 *x + 1.1d-10* (T /100.d0)**(0.44) )
  !    c      * exp(-98.d0 / T)


  !       froido = 2.5d-27 * (T/100)**0.4 * exp(-228.d0 / T)


  !on tient compte de l'abondance du
  froido = froido * 4.5E-4


  !refroidissement par l'hydrogene
  ! formule de Spitzer 1978
! STG HACK - ONLY RUN THIS IF WE DON'T HAVE RT
  froidh = 0d0
  if (.not. rt) then
     froidh = 7.3E-19 * x * exp(-118400./ T )
  endif

  !refroidissement par les raies metastables des metaux
  !chiffre pris dans Hollenbach and McKee 1989 (ApJ 342, 306)





  !carbone une fois ionise ,1 transitions 2P 4P
  ! le poids est 1
  ! 2P->4P :
  ! les expressions des coefficients d'excitation ont une dependance
  !en la temperature differente au dela de 10000K
  !abondance 3.5 d-4 depletion 0.4

         froidcII_m = 6.2d4 * 1.38d-16 * 1.d0 * &    !transition 2P->4P
        ( 2.3d-8* (T/10000.)**(-0.5) * x + 1.d-12 ) *exp(-6.2d4 / T) &
           * 3.5d-4 * 0.4




         if ( T .le. 1.d4 ) then
         froido_m = 2.3d4 * 1.38d-16 / 3.d0 * &
        ( 5.1d-9 * (T/10000.)**(0.57) * x + 1.d-12) *exp(-2.3d4/T)
  
         froido_m = froido_m + &
              4.9d4 * 1.38d-16 / 3.d0  * &
        ( 2.5d-9 * (T/10000.)**(0.57) * x + 1.d-12) *exp(-4.9d4/T)
  

         froido_m = froido_m + &
              2.6d4 * 1.38d-16 * 1.d0  * &
        ( 5.2d-9 * (T/10000.)**(0.57) * x + 1.d-12) *exp(-2.6d4/T)

         else

         froido_m = 2.3d4 * 1.38d-16 / 3.d0 * &
        ( 5.1d-9 * (T/10000.)**(0.17) * x + 1.d-12) *exp(-2.3d4/T)
  
         froido_m = froido_m + &
              4.9d4 * 1.38d-16 / 3.d0  * &
        ( 2.5d-9 * (T/10000.)**(0.13) * x + 1.d-12) *exp(-4.9d4/T)


         froido_m = froido_m + &
              2.6d4 * 1.38d-16 * 1.d0  * &
        ( 5.2d-9 * (T/10000.)**(0.15) * x + 1.d-12) *exp(-2.6d4/T)


         endif

  !! abondance de l'oxygene
         froido_m = froido_m *   4.5d-4



!!! on somme les refroidissements
  froid = froidcII  + froidh  + froido  + froido_m +  froidcII_m


  !      froid=froid*1.d-13    !conversion en MKS


  !refroidissement par le carbone neutre. On suppose l'equilibre
  ! de la reaction C + hv <=> C+ + e-
  ! les taux de reactions et de refroidissement sont pris dans
  !la these de Karl Joulain.

  ! abondance du carbone relative a n (MKS)


  !    C+ + e- => C
  !       k1 = 4.4d-12 * (T/300.)**(-0.61) !s^-1 cm^-3

  !       k1 = k1


  !    C => C+ + e-
  !       k2 = 2.2d-10


  ! on a : [C] = k1/k2 [C+] * [e-]
  ! on suppose que tout le carbone est sous forme C+
  ! et que [e-] = [C+]

  ! l'abondance relative de carbone
  !      nc = k1/k2 * (3.5d-4*0.4)**2 * n


  !      froidc =  1.0d-24 * ( 1.4d0 * exp( -23.d0 / T ) &
  !                     + 3.8d0 * exp( -62.d0 / T )   )

  !      froidc = froidc * nc !(nc est l'abondance relative du carbone)


  !       n=exp(log(10.d0)*logn) !ici pas besoin de log

  !       valeur utilisees au 23/08/98
  !       chaud=4.d0*exp(-24.5d0*log(10.d0))*1.d-7  !un peu empirique ....



!!!! on calcule le chauffage
!!! on prend en compte le chauffage sur les grains
!!! formules 1 et 2  de Wolfire et al. 1995

!!!! G0 est le flux UV par rapport au flux defini par Habing et
!!!! Draine

  G0 = 1./1.7

  param = G0 * sqrt(T)/(n*x)
  epsilon = 4.9E-2 / (1. + (param/1925.)**0.73)
  epsilon  = epsilon + 3.7E-2 * (T/1.E4)**0.7 / (1. + (param/5.E3) )


  chaud = 1.E-24 * epsilon


  ! pour un flux de rayonnement G0/1.7
  chaud = chaud * G0

  !refroidissement recombinaison sur les grains charges positivement
  bet = 0.74/(T**0.068)
  froidrec = 4.65E-30*(T**0.94)*(param**bet)*x



  !! chaud=1.d-32 !un peu empirique ....

  !      froidc=0.d0


  ref= chaud*n - (n**2)*(froid + froidrec) !!!+ froidc)

  return

end subroutine cool_rate_lowT

end module cooling_module

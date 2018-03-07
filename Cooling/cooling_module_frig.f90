!! comes from frig version and then all the specific development done by Valeska
!! to take into account extinction have been moved there
!! PH 19/01/2017
!=======================================================================
subroutine solve_cooling_frig(nH,T2,zsolar,dt,gamma,ncell,deltaT2)
!=======================================================================
  implicit none
  ! BRIDGE FUNCTION WITH SAME INTERFACE AS SOLVE_COOLING 
  ! Input/output variables to this function
  ! nH - hydrogen number density in PHYSICAL units
  ! T2 - temperature / mu in PHYSICAL units
  ! zsolar - Metallicity in solar units (Zphys / 0.02)
  ! dt - cooling timestep in seconds
  ! ncell - number of elements in the vector
  ! deltaT2 - temperature change in K/mu (??)
  integer,intent(in)::ncell
  real(kind=8),intent(in)::dt,gamma
  real(kind=8),dimension(1:ncell),intent(in)::nH,T2,zsolar
  real(kind=8),dimension(1:ncell),intent(out)::deltaT2
  ! Input/output variables to analytic function calc_temp 
  real(kind=8)::NN,TT, dt_tot
  ! Temporary variables
  integer::i
  real(kind=8)::TT_ini, mu
  ! Loop over cells
  dt_tot = dt
  do i=1,ncell
     NN = nH(i) ! NOTE!! THE CODE BELOW ASSUMES scale_nH=1 !!
                ! SO WE LEAVE THIS AS IT IS TO KEEP UNITS CONSISTENCY
     TT = T2(i)
     TT_ini = TT
     call calc_temp(NN,TT,dt_tot,gamma)
     deltaT2(i) = (TT - TT_ini)
  end do
end subroutine solve_cooling_frig

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine  calc_temp(NN,TT,dt_tot,gamma)
    !use amr_parameters
    !use hydro_commons

    implicit none

    integer :: n,i,j,k,idim, iter, itermax,ii

    real(kind=8) :: dt, dt_tot, temps, dt_max, itermoy
    real(kind=8) :: rho,temp

    !alpha replaced by alpha_ct because of conflict with another alpha by PH 19/01/2017
    real(kind=8) :: mm,uma, kb, alpha_ct,mu,kb_mm
    real(kind=8) :: NN,TT, TTold, ref,ref2,dRefdT, eps, vardt,varrel, dTemp,dummy
    real(kind=8) :: rhoutot2,gamma
    ! HARD-CODED mu TO MAKE TEMPERATURE AGREE WITH HENNEBELLE CODE
    mu = 1.4d0
    !
    ! Cette routine fonctionne en cgs
    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    kb  =  1.38062d-16   ! erg/degre
    !  uma =  1.660531e-24  ! gramme
    !  mu  =  1.4
    !  mm = mu*uma
    !  kb_mm = kb / mm
    !  TT = TT  / kb  !/ kb_mm

    if( TT .le. 0.) then
        TT = 10.
        return
    endif

    vardt = 10.**(1./10.); varrel = 0.2

    !  nn = (rho/(gramme/cm3)) /mm

    itermax = 0 ; itermoy = 0.



    ! if (NN .le. smallr) then
    !     if( NN .le. 0)  write(*,*) 'prob dens',NN
    !     NN = smallr  !max(NN,smallr)
    ! endif


    ! alpha_ct = NN*kb_mm/(gamma-1.)
    alpha_ct = NN*kb/(gamma-1.)

    ! eps - a small offset of T to find gradient in T
    eps = 1d-5

    iter  = 0 ; temps = 0.
    do while ( temps < dt_tot)
        if (TT .lt.0) then
            write(*,*) 'prob Temp',TT, NN
            !         write(*,*) 'repair assuming isobariticity'
            !NN = max(NN,smallr)
            TT = min(4000./NN,8000.)  !2.*4000. / NN
        endif


        TTold = TT

        ! Calculate cooling rate
        !NN is assumed to be in cc and TT in Kelvin
        if (TT < 10035.d0) then
            call cooling_low(TT,NN,ref,dummy)
            call cooling_low(TT*(1d0+eps),NN,ref2,dummy)
        else
            call cooling_high(TT,NN,ref)
            call cooling_high(TT*(1d0+eps),NN,ref2)
        end if
        
        ! dT = T*(1+eps)-T = eps*T
        dRefdT = (ref2-ref)/(TT*eps)

        ! TODO STG - COPY THIS FUNCTION UP TO HERE, USE ref, drefdT TO 
        !            REPLACE rt_cmp_metals SOMEHOW


        !       write(*,*) 'check',TTold, TT,NN,ref,dRefdT,iter


        if (iter == 0) then
            if (dRefDT .ne. 0.) then
                dt = abs(1.0E-1 * alpha_ct/dRefDT)
            else
                dt = 1.0E-1 * dt_tot
            endif
            dt_max = dt_tot - temps
            if (dt > 0.7*dt_max) dt = dt_max*(1.+1.0E-12)
        endif

        dTemp = ref/(alpha_ct/dt - dRefdT)

        eps = abs(dTemp/TT)
        if (eps > 0.2) dTemp = 0.2*TTold*dTemp/abs(dTemp)

        TT = TTold + dTemp
        if (TT < 0.) then
            write(*,*) 'Temperature negative !!!'
            write(*,*) 'TTold,TT   = ',TTold,TT
            write(*,*) 'rho   = ',rho
            TT = 10.  !*kelvin
        endif


        iter = iter + 1

        temps = temps + dt

        dt = vardt*varrel*dt/Max(vardt*eps, varrel)

        dt_max = dt_tot - temps
        if (dt > 0.7*dt_max) dt = dt_max*(1.+1.0E-12)
        !        write(*,987) temps, TT
        !987     format(E10.3,2x,E10.3)
        !        read (*,*)
    enddo

    return
end subroutine calc_temp

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine cooling_high(T,n,ref)
    !use amr_parameters
    implicit none

    real(kind=8) :: T,P,N,x,ne                    !x est le taux d'ionisation
    real(kind=8) :: T2, ref2
    real(kind=8) :: froid,chaud,ref,nc, froidc
    real(kind=8) :: froidcII, froido, froidh
    real(kind=8) :: froidc_m,froidcII_m,froido_m
    real(kind=8) :: param, G0, epsilon,k1,k2,bet,froidrec
    real(kind=8) :: eps

    real(kind=8) :: logT, intcst, logT2
    real(kind=8) :: ion, neut
    real(kind=8) :: refion

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

!    chaud = 1.E-25
    chaud = 0.
    ref= chaud*n + (n**2)*(froid)

end subroutine cooling_high

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine cooling_low(T,n,ref,chaud)
! Note: (et oui c'est en anglais maintenant, desole)
! T,n - temperature and density input
! ref - *total* cooling rate with cooling as -ve (i.e. heating - cooling)
! chaud - heating rate *only* 
!    (used for mixed cooling functions in, e.g., rt_metal_cooling_module)

  !use amr_parameters

  implicit none

  real(kind=8) :: T,P,N,x,ne,x_ana                 !x est le taux d'ionisation
  real(kind=8) :: froid,chaud,ref,nc, froidc
  real(kind=8) :: froidcII, froido, froidh
  real(kind=8) :: froidc_m,froidcII_m,froido_m
  real(kind=8) :: param, G0, epsilon,k1,k2,bet,froidrec


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



  ne = 2.4d-3*((T/100d0)**0.25d0)/0.5d0 !formule C15 Wolfire et al. 2003

  ! Analytic ionisation in absence of photoionisation
  x_ana = ne / N   ! ionisation
  x_ana = min(x_ana,0.1d0)
  x_ana = max(x_ana,3.5d-4*0.4d0)
  x = x_ana ! (Different variables in case we need photoionised values later)

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
  froidh = 7.3E-19 * x * exp(-118400./ T )

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

  chaud = chaud*n
  ref= chaud - (n**2)*(froid + froidrec) !!!+ froidc)

  return

end subroutine cooling_low

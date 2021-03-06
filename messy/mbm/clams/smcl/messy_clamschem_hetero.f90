Module messy_clamschem_hetero

contains

subroutine hetero
!
!   version of subroutine hetero that calculates the sum of heterogeneous
!   reaction rates on liquid, NAT and ice aerosol.
!   In this version chindex has been changed to an 2-dimensional array.
!   for "no available reaction" chindex=numhet+1 is used
!
!   J.U. Grooss, 4.11.1998
!    
!   05.07.1999 Fortran 90 version of this routine, intermediate routine rufhet
!              is not needed anymore
!   17.03.2003 workaround for HCl = 0 that caused some artificial high reaction rates
!   05.01.2005 Implemented new structure for sedimenation, ensure
!              that only 90% sediments (c.lemmen)
!
 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
  USE messy_clams_global,        ONLY: prec
  use messy_clamschem_global,    only: ntraj, missing_index, &
                                       parth, wt, ar, con, densaero, aerh2so4, &
                                       sedinucl, &
                                       shindex, chindex

  use messy_clamschem_hetero_shi,only: hetero_shi

  use messy_clamschem_globalhet, only: numhet,itr0, &
                                       densnat, densice, denssat, &
                                       cnatinit, ciceinit, aer_h2so4_default, &
                                       satmelting, liquids, &
                                       teold, prold, kstate, astate, &
                                       laststate, lstate, &
                                       told, pressold, densnat_old, densice_old, &
                                       denssat_old, cnat_old, cice_old, mixedpop, &
                                       log_state,  natcore, t_lt_tnat, vliq_save

  use messy_clamschem_asad_mod,       only: f, p, t, tnd, rk, nhrkx
  use messy_clamschem_asad_mod_clams, only: hetdt, jphk

  implicit none 
  
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
  integer ::  ij,  itr 
  real(PREC)  :: het1(numhet+1) 
  real(PREC) :: ltemp, lpress, chcl, chocl, ch2o, cclno3, cbrno3, chbr, &
       chobr,  chno3, parthcl, parthcl0, parthno3, parthno30, parth2o, &
       parth2o0, parthbr, parthbr0, aliq, asat, anat, aice, asatliq, wts, &
       wtn, wtcl, wtbr, wthocl, wthobr, xp, vliq, snat 

  real(PREC) :: terminal_velocity(2),mean_free_path,mol_gas_const&
       &,mol_weight_air,viscosity, denschange
  real(PREC), parameter :: pi = 3.1415926 
  logical :: tincrease
!-----------------------------------------------
!
!
!******************************************
  ! loop over air parcels
  do itr = 1, ntraj 
    if (missing_index(itr)) cycle

    ! change of density that applies to ftr in cdrive must also
    ! be considered for con !!!   (jug, 05.04.05)
    denschange= TeOLD(itr)/t(itr) * 0.01* p(itr)/prold(itr)
    tincrease=(t(itr) > TeOLD(itr))

    ! add condensed phase mixing ratios to gasphase 
    ! HCl in molec/ccm
    f(itr,shindex(1)) = f(itr,shindex(1)) + con(itr,1) * denschange
    ! HNO3 in molec/ccm
    f(itr,shindex(2)) = f(itr,shindex(2)) + con(itr,2) * denschange
    ! H2O  in molec/ccm
    f(itr,shindex(3)) = f(itr,shindex(3)) + con(itr,3) * denschange
    ! HBr  in molec/ccm
    f(itr,shindex(4)) = f(itr,shindex(4)) + con(itr,4) * denschange



    ! TRANSFER INPUT VARIABLES TO INTERNAL PARAMETERS
    ltemp   = t(itr) 
    lpress  = p(itr)/100.          ! Pressure in hPa
    chcl   = f(itr,shindex(1))    ! HCl in molec/ccm
    chno3  = f(itr,shindex(2))    ! HNO3 in molec/ccm
    ch2o   = f(itr,shindex(3))    ! H2O  in molec/ccm
    chbr   = f(itr,shindex(4))    ! HBr  in molec/ccm
    chocl  = f(itr,shindex(5))    ! HOCl in molec/ccm
    cclno3 = f(itr,shindex(6))    ! ClONO2 in molec/ccm
    chobr  = f(itr,shindex(7))    ! HOBr in molec/ccm
    cbrno3 = f(itr,shindex(8))    ! BrONO2 in molec/ccm

    parthcl  = parth(itr,1,1)     ! (old) gas fraction of HCl (0<parthcl<1)
    parthno3 = parth(itr,1,2)     ! (old) gas fraction of HNO3 (0<part<1)
    parth2o  = parth(itr,1,3)     ! (old) gas fraction of H2O (0<part<1)
    parthbr  = parth(itr,1,4)     ! (old) gas fraction of HBr (0<part<1)
    xp       = densaero(itr)      ! (old) liquid aerosol density

    kstate        = astate(itr)
    laststate     = lstate(itr)
    TOLD          = TeOLD(itr)
    pressold      = Prold(itr)
    densnat       = densnat_old(itr)
    densice       = densice_old(itr)
    denssat       = denssat_old(itr)
    cnatinit      = cnat_old(itr)
    ciceinit      = cice_old(itr)
    aer_h2so4_default = aerh2so4(itr)
    mixedpop      = log_state(1,itr)
    natcore       = log_state(2,itr)
    satmelting    = log_state(3,itr)
    liquids       = log_state(4,itr)


  ! STORE OLD VALUES OF HET AND OUTPUT VARIABLES
    parthcl0 = parthcl 
    parthno30 = parthno3 
    parth2o0 = parth2o 
    parthbr0 = parthbr 

    ! Berechnung der mittleren freien Weglaenge (Seinfeld, S. 455)
    ! Angabe des Drucks in Pa, Temperatur in K, alles andere s. u.
    ! Das Ergebnis hat die Dimension m, darum wird mit dem Faktor 100 in cm
    ! fuer die Uebergabe an hetero_shi multipliziert.

    viscosity = 1.72e-5      ! [kg/(m s)]
    mol_weight_air = 28.97   ! [kg/kmol]
    mol_gas_const = 8314.    ! [J/(kmol K)]

    mean_free_path=100.*(2.*viscosity)/(p(itr)*sqrt(8.*mol_weight_air/&
         &(pi*mol_gas_const*t(itr))))

    ! Festlegen des HCl-Mischungsverhaeltnis fuer die Rechung auf mind. 1ppt
!!!!! tnd (in messy_clamschem_asad_mod) vom Typ REAL, chcl vom Typ REAL(RP) 
    !     -> evtl. unterschiedl. -> Warning beim Compilieren!
    chcl=max(chcl,tnd(itr)*1E-12)

    itr0= itr
    call hetero_shi (lpress, ltemp, chno3, chcl, ch2o, chocl, chbr, chobr, &
         cclno3, cbrno3, parthno3, parthcl, parthbr, parth2o, wts, wtn, wtcl, &
         wtbr, wthocl, wthobr, aice, anat, aliq, asat, asatliq, xp, het1, &
         & terminal_velocity,mean_free_path,vliq,snat)


  ! ===============================================================
  ! TRANSFER RESULTS OF OUTPUT ARRAY FOR TRANSFER TO ASAD:

    parth(itr,1,1) = parthcl 
    parth(itr,1,2) = parthno3 
    parth(itr,1,3) = parth2o 
    parth(itr,1,4) = parthbr 
    parth(itr,2,1) = parthcl0 
    parth(itr,2,2) = parthno30 
    parth(itr,2,3) = parth2o0 
    parth(itr,2,4) = parthbr0 

    wt(itr,1)      = wts      ! Weight fraction of Sulphur
    wt(itr,2)      = wtn      ! Weight fraction of Nitrogen
    wt(itr,3)      = wtcl     ! Weight fraction of Chlorine
    wt(itr,4)      = wtbr     ! Weight fraction of Bromine
    wt(itr,5)      = wthocl   ! Weight fraction of HOCl
    wt(itr,6)      = wthobr   ! Weight fraction of HOBr
    ar(itr,1)      = aice     ! surface area ice
    ar(itr,2)      = anat     ! surface area NAT 
    ar(itr,3)      = aliq     ! surface area liq
    ar(itr,4)      = asat     ! surface area SAT
    ar(itr,5)      = asatliq  ! surface area SATLIQ
    densaero(itr)  = xp       ! areosol number density

!    sedinucl(itr)= (astate(itr) <=2 .and. kstate > 2 )
!!$    sedinucl(itr)  = t_lt_tnat
    astate(itr)    =  kstate 
    lstate(itr)    =  laststate 
    TeOLD(itr)     =  t(itr) ! was TOLD
    Prold(itr)     =  p(itr)/100. ! was pressold
    densnat_old(itr) = densnat
    densice_old(itr) = densice
    denssat_old(itr) = denssat      
    cnat_old(itr)    = cnatinit
    cice_old(itr)    = ciceinit
    log_state(1,itr) = mixedpop
    log_state(2,itr) = natcore
    log_state(3,itr) = satmelting
    log_state(4,itr) = liquids
    vliq_save(itr) = vliq


!!$    ! jug/ 2/2012 For increasing temperatures, calculated snat may be below 1
!!$    ! while still NAT exists (as it is calculated using hno3*parthno3) 
!!$    if (tincrease .and. parthno3 < 0.9999) snat = max(snat,snatmax(itr))
!!$
!!$    ! jug, 2/2012 save snatmax and tmin for sedi...
!!$    if (snat > 1. .and. snat > snatmax(itr))  then
!!$       snatmax(itr)=snat
!!$       tmin(itr)=t(itr)
!!$    endif
!!$
!!$    ! jug, 2/2012 reset snatmax and tmin, if T > T_NAT
!!$    if (snat < 0.95 )  then 
!!$       snatmax(itr)=0.
!!$       tmin(itr)=t(itr)
!!$    endif
!!$
!!$
!!$    ! jug, 7/2013 on the first day of NAT particles
!!$    !(snat0 < 1, snatmax gt 1) set snatmax to 1.0
!!$    if (snat0(itr) < 1.0 .and. snat > 1. ) snatmax(itr) = 1.0001


    ! re-partitioning of gasphas and condensed phase
    con(itr,1) = f(itr,shindex(1))*(1.0 - parth(itr,1,1)) 
    con(itr,2) = f(itr,shindex(2))*(1.0 - parth(itr,1,2)) 
    con(itr,3) = f(itr,shindex(3))*(1.0 - parth(itr,1,3)) 
    con(itr,4) = f(itr,shindex(4))*(1.0 - parth(itr,1,4)) 
 

    ! HCl in molec/ccm
    f(itr,shindex(1)) = f(itr,shindex(1)) - con(itr,1)
    ! HNO3 in molec/ccm
    f(itr,shindex(2)) = f(itr,shindex(2)) - con(itr,2) 
    ! H2O  in molec/ccm
    f(itr,shindex(3)) = f(itr,shindex(3)) - con(itr,3)
    ! HBr  in molec/ccm
    f(itr,shindex(4)) = f(itr,shindex(4)) - con(itr,4) 

    ! set reaction rate coefficients:
    do ij = 1, jphk 
       ! sum up reaction rates on all 4 substrates
       rk(itr,nhrkx(ij)) = sum(het1(chindex(ij,:))) 
    end do
  end do
  return  
end subroutine hetero 

end Module messy_clamschem_hetero

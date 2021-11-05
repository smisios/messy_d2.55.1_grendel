MODULE messy_ddep

  ! MESSy
  USE messy_main_constants_mem, ONLY: DP, SP

  IMPLICIT NONE
  PRIVATE
  SAVE

  CHARACTER(len=*), PARAMETER,  PUBLIC :: MODSTR = 'ddep'
  CHARACTER(len=*), PARAMETER,  PUBLIC :: modver = '2.1'
  
  PUBLIC :: dp, sp

  INTRINSIC  ABS, ATAN, EXP, LOG, MAX, MIN, SQRT , TINY

  ! CTRL-NAMELIST PARAMETER
  ! switch for considering role of whitecaps
  LOGICAL, PUBLIC :: l_whitecap = .false.
  ! switch for considering relative humidity effect
  LOGICAL, PUBLIC :: l_rh = .false.
  ! switch for considering wetness at leaf surface
  LOGICAL, PUBLIC :: l_ganzeori = .true. !ju_te_20180622
  
  ! SO2
  REAL(dp), PARAMETER :: diffrb_so2=1.6_dp
  REAL(dp), PARAMETER :: rsoil_so2=250._dp
  REAL(dp), PARAMETER :: rwater_so2=1._dp
  REAL(dp), PARAMETER :: rws_so2=100._dp
  REAL(dp), PARAMETER :: rsnow_so2=1._dp
  REAL(dp), PARAMETER :: rmes_so2=1._dp
  REAL(dp), PARAMETER :: rcut_so2=1.e5_dp
  REAL(dp), PARAMETER :: diff_so2=1.9_dp

  ! O3
  REAL(dp), PARAMETER :: diffrb_o3=1.2_dp
  REAL(dp), PARAMETER :: rsoil_o3=400._dp
  REAL(dp), PARAMETER :: rwater_o3=2000._dp
  REAL(dp), PARAMETER :: rws_o3=2000._dp
  REAL(dp), PARAMETER :: rsnow_o3=2000._dp
  REAL(dp), PARAMETER :: rmes_o3=1._dp
  REAL(dp), PARAMETER :: rcut_o3=1.e5_dp
  REAL(dp), PARAMETER :: diff_o3=1.6_dp
  !

  ! declaration of the the maximum number of vegetation
  ! layers. This should resemble the number of layers being distinguished
  ! in the input vegetation files (default 2 layers, must be extended
  ! to 4 layers for hight resolution calculations of radiation profiles.
  ! This is required due to the sensitivity of the isoprene emission
  ! calculations for the vertical discretazation, October 2001).

  ! variables exchanged between drydep_vdbl_parameter and  drydep_vdbl
  REAL(dp), ALLOCATABLE ::  rahcan(:)
  REAL(dp), ALLOCATABLE ::  rsnowhno3so2(:)
  REAL(dp), ALLOCATABLE ::  rsoilso2(:)

  ! SUBROUTINES
  PUBLIC :: ddep_read_nml_ctrl
  PUBLIC :: drydep_calc_rs                 ! subroutine
  PUBLIC :: drydep_calc_ra                 ! subroutine
  PUBLIC :: drydep_vdbl_parameter
  PUBLIC :: drydep_vdbl
  PUBLIC :: drydep_vdbl_dealloc
  
  PUBLIC :: drydep_vdaer_parameter
  PUBLIC :: drydep_vdaer
  
  PUBLIC :: drydep_calc_abswind           ! function
  PUBLIC :: drydep_calc_radius_mass       ! function
  PUBLIC :: drydep_calc_radius_num        ! function
  PUBLIC :: drydep_calc_landtypefractions ! subroutine
  PUBLIC :: drydep_calc_vegfrac           ! subroutine
  PUBLIC :: drydep_calc_layerthickness    ! function
  PUBLIC :: drydep_calc_vd_eff            ! function
  PUBLIC :: drydep_posfinit               ! function
  PUBLIC :: drydep_calc_rcut              ! subroutine ! ju_te_20180622
CONTAINS

  !===========================================================================
  
  SUBROUTINE ddep_read_nml_ctrl(status, iou)

    ! DRYDEP MODULE ROUTINE (CORE)
    !
    ! READ DRYDEP NAMELIST, CHECK IT, AND INITIALIZE GLOBAL VARIABLES
    !
    ! Author: Patrick Joeckel, MPICH, Feb 2002
    ! Modified: Laurens Ganzeveld, MPICH, 14-11-2002

    USE messy_main_tools, ONLY: read_nml_open, read_nml_check, read_nml_close

    IMPLICIT NONE

    ! I/O
    INTEGER, INTENT(OUT) :: status
    INTEGER, INTENT(IN)  :: iou   ! logical I/O unit

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER       :: substr='drydep_read_nml_ctrl'
    LOGICAL                           :: lex          ! file exists ?
    INTEGER                           :: fstat        ! file status

    NAMELIST /CTRL/  l_whitecap, l_rh, l_ganzeori ! ju_te_20180622

    status = 1 ! ERROR ON RETURN

    CALL read_nml_open(lex, substr, iou, 'CTRL', modstr)
    IF (.not.lex) RETURN    ! <modstr>.nml does not exist

    READ(iou, NML=CTRL, IOSTAT=fstat)
    CALL read_nml_check(fstat, substr, iou, 'CTRL', modstr)
    IF (fstat /= 0) RETURN  ! error while reading namelist

    WRITE(*,*) 'Switch wetness dependent calculation of the cuticular resistance ... '
    IF (.NOT.l_ganzeori) THEN
      WRITE(*,*) '... consider wetness (Zhang 2002)'
    ELSE 
      WRITE(*,*) '... no wetness dependence - old scheme (Wesely 1998)'
    END IF

    ! DIAGNOSE NAMELIST AND SET GLOBAL SWITCHES
    CALL read_nml_close(substr, iou, modstr)

    status = 0

  END SUBROUTINE ddep_read_nml_ctrl

  !========================================================================

  SUBROUTINE drydep_calc_rs(name, tracnum, number,henry, molweight, dryreac &
                         , diff, diffrb, rmes, rsoil, rwater, rsnow, rws, rcut)
 
    IMPLICIT NONE

    CHARACTER(LEN=*), INTENT(IN)  :: name
    INTEGER, INTENT(IN)  :: number, tracnum
    REAL(dp), INTENT(IN) :: henry
    REAL(dp), INTENT(IN) :: molweight
    REAL(dp), INTENT(IN) :: dryreac
    
    REAL(dp), INTENT(INOUT), DIMENSION(tracnum) :: &
         diff, diffrb, rmes, rsoil, rwater, rsnow  ! ju_te_20180625 
    REAL(dp), INTENT(INOUT), OPTIONAL,  DIMENSION(tracnum) :: rcut, rws
    INTRINSIC SQRT, TRIM
    LOGICAL :: l_dry =.FALSE.    ! ju_te_20180629

    !--------------------------------------------------------------------------
    ! subroutine drydep_calc_rs, to calculate the values of the uptake 
    ! resistances required to calculate the trace gas dry deposition velocity. 
    ! this routine is based on an approach by wesely, 1989, in which the uptake
    ! resistances of trace gases, for which the dry deposition velocities have
    ! not been observed, are estimated based on the henry coefficient and a 
    ! reactivity coefficient and the uptake resistances of so2 and o3, of the
    ! "big leaf" dry deposition scheme
    ! by ganzeveld and j. lelieveld j. geophys. res., 100, 20,999-21,012,1995,
    ! ganzeveld et al.,j. geophys. res., 103, 5679-5694, 1998 and 
    ! ganzeveld et al,  submitted to j. geophys. res., 2001. 
    ! for more information of the wesely
    ! approach see atmospheric environment vol 23, no 6, 1293-1304.
    !
    ! the program needs as input data the molecular mass of the defined trace
    ! gases, the henry coefficient [m atm-1] and an estimated reactivity
    ! coefficient which has 3 distinct values: 0 for non-reactive species,
    ! (e.g, so2, acetaldehyde), 0.1 for moderately reactive species (e.g. pan),
    ! and 1 for reactive species (e.g., o3, hno3).
    !-------------------------------------------------------------------------

    ! definition of the different resistances for the 6
    ! species of the original dry deposition scheme:
    !     ------------o3, hno3, no, no2, so2 and so4--------------
  
    ! ju_te_20180629+
    ! flag for use of the optional variables
     IF (PRESENT (rcut) .AND. PRESENT(rws)) THEN
       l_dry = .TRUE.   
     ELSE 
       l_dry = .FALSE.
     END IF
    ! ju_te_20180629-

     SELECT CASE(TRIM(name))

     CASE('SO2', 'H2SO4')

       diffrb(number)=diffrb_so2
       rsoil(number)=rsoil_so2
       rwater(number)=rwater_so2
       rsnow(number)=rsnow_so2
       rmes(number)=rmes_so2
       diff(number)=diff_so2
       IF(l_dry) THEN  ! ju_te_20180629
         rcut(number)=rcut_so2
         rws(number)=rws_so2
       END IF ! ju_te_20180629

     CASE('O3')

       diffrb(number)=diffrb_o3
       rsoil(number)=rsoil_o3
       rwater(number)=rwater_o3
       rsnow(number)=rsnow_o3
       rmes(number)=rmes_o3
       diff(number)=diff_o3
       IF(l_dry) THEN ! ju_te_20180629
         rcut(number)=rcut_o3
         rws(number)=rws_o3
       END IF         ! ju_te_20180629
      
    CASE('SO4')

       diffrb(number)=1.8_dp
       rsoil(number)=1.e5_dp
       rwater(number)=1.e5_dp
       rsnow(number)=1.e5_dp
       rmes(number)=1.e5_dp
       diff(number)=2.7_dp
       IF(l_dry) THEN          ! ju_te_20180629
         rcut(number)=1.e5_dp
         rws(number)=1.e5_dp
       END IF                  ! ju_te_20180629

     CASE('HNO3')

       diffrb(number)=1.4_dp
       rsoil(number)=1._dp
       rwater(number)=1._dp
       rsnow(number)=1._dp
       rmes(number)=1._dp
       diff(number)=1.9_dp
       IF(l_dry) THEN        ! ju_te_20180629
         rcut(number)=1._dp
         rws(number)=1._dp
       END IF                ! ju_te_20180629

     CASE('NO')

       diffrb(number)=1.1_dp
       rsoil(number)=1.e5_dp
       rwater(number)=1.e5_dp
       rsnow(number)=1.e5_dp
       rmes(number)=500._dp
       diff(number)=1.3_dp
       IF(l_dry) THEN         ! ju_te_20180629
         rcut(number)=1.e5_dp
         rws(number)=1.e5_dp
       END IF                 ! ju_te_20180629

     CASE('NO2')

       diffrb(number)=1.2_dp
       rsoil(number)=600._dp
       rwater(number)=1.e5_dp
       rsnow(number)=1.e5_dp
       rmes(number)=1._dp
       diff(number)=1.6_dp
       IF(l_dry) THEN          ! ju_te_20180629
          rcut(number)=1.e5_dp
          rws(number)=1.e5_dp
       END IF                  ! ju_te_20180629
       
    CASE DEFAULT

        IF (name/='SO2'  .AND. name/='O3') THEN

           ! calculation of term which is used to correct the stomatal 
           ! resistance for differences in the diffusitivy 
           
           diff(number)=SQRT(molweight/18._dp)

           ! calculation of the term to correct for differences in diffusivity
           ! between water vapor and the trace gas. it is calculated from:
           ! diff bl=(v/dx)**2/3, with v=0.189 sm-1 and dx= dh2o/sqrt(mh2o/mx),
           ! with dh2o=0.212

           diffrb(number)=(0.189_dp/(0.212_dp/diff(number)))**(2._dp/3._dp)

           ! calculation of rmx, the mesophyll resistance

           rmes(number)=1._dp/(henry/3000._dp+100._dp*dryreac)

          ! calculation of rlux, the cuticular resistance, equation 7 of 
          ! Wesely's  paper

           IF(l_dry) THEN   ! ju_te_20180629
             rcut(number)=1._dp/(1.e-5_dp*henry+dryreac)*rcut_o3
             rws(number)=1._dp/(1._dp/(3._dp*rws_so2)+1.e-7_dp*henry+ &
                     dryreac/rws_o3)
           END IF           ! ju_te_20180629

           ! calculation of rgsx, the soil resistance, equation 9 of wesely's
           ! paper

           rsoil(number)=1._dp/(henry/(1.e5_dp*rsoil_so2)+ &
                dryreac/rsoil_o3)

           ! the snow resistance is similar as the soil resistance
           rsnow(number)=1._dp/(henry/(1.e5_dp*rsnow_so2)+ &
                dryreac/rsnow_o3)

           ! calculation of sea uptake resistance, using equation 9 of wesely's
           ! paper

           rwater(number)=1._dp/(henry/(1.e5_dp*rwater_so2)+ &
                dryreac/rwater_o3)

        ENDIF
        
     END SELECT

  END SUBROUTINE drydep_calc_rs
!=============================================================================

!=============================================================================
! ju_te_20180629+
  SUBROUTINE drydep_calc_rcut(kproma, prh_2m, pustveg,     &
         plai, dryreac, henry, rsfl, rsfc, rcut_2d, rws_2d) 

    !=======================================================
    ! Author: Tamara Emmerichs, FZ Juelich, 26-06-2018
    ! t.emmerichs@juelich.de
    ! drydep_calc_rcut to calculate the cuticular resistance
    ! over dry and wet canopies at canopy scale according to
    ! Zhang et al., Atm. Environment 36 (2002) 4787-4799
    ! new -  dependency on friction velocity, LAI, RH
    ! for all tracer - scaling with henry constant,
    ! reactivity of each tracer (approach by Wesely 1989)
    ! Zhang et al., Atmos. Chem. Phys., 3, 2067-2082, 2003
    !========================================================
    
    IMPLICIT NONE
 
    ! Formalparameter
    INTEGER, INTENT(IN) :: kproma
    REAL(dp), DIMENSION(1:kproma), INTENT(IN)  :: plai,    &
                               pustveg, prh_2m, rsfl, rsfc
    REAL(dp), INTENT(IN)  :: dryreac, henry

    !LOCAL
    INTEGER :: jl, rcutw0_so2
    REAL(dp)  :: rcutw_o3, rcutd_o3, rcutw_so2,  &
                 ftrac
    REAL(dp)  :: pplai, ppustveg
    ! input parameters in s/m see Zhang et al. 2002 (Table 2)
    INTEGER, PARAMETER :: rcutd0=5000._dp, rcutw0=300._dp,   &
                          rcutd0_so2=2000._dp
    REAL(dp), DIMENSION(1:kproma), INTENT(OUT)  :: rcut_2d, &
               rws_2d
   
    ! tracer-dependent factor for the resistance over dry canopy
    ftrac = 1.e-5_dp*henry+dryreac                                 
    ! ftrac_wet= 1._dp/(3._dp*rws_so2)+1.e-7_dp*henry 
    ! tracer-dependent factor for the resistance over wet canopy

    DO jl=1,kproma
       ! threshold avoid negative LAI over the ocean
      pplai = MAX(plai(jl),1.e-5_dp)                               
      ppustveg = MAX(0.1_dp,pustveg(jl))
      ! cuticular resistance of ozone above dry canopy (s/m), upper threshold 
      ! from L.Zhang, J.R.Brook and R.Vet (2003)
      !rcutd_so2 = rcutd0_so2/MAX(1.e-10_dp,(exp(0.03_dp*prh_2m(jl) & 
      !              *100._dp)*pplai**(0.25_dp)*ppustveg))
      IF ((rsfl(jl)+rsfc(jl)) > 0._dp) THEN
        rcutw0_so2 = 50._dp                                        ! rain
      ELSE
        rcutw0_so2 = 100._dp                                       ! dew
      END IF
      ! cuticular resistance of sulphur dioxide above wet canopy (s/m)
      rcutw_so2 = rcutw0_so2/(pplai**(0.5_dp) *ppustveg)           

      ! cuticular resistance of ozone above dry canopy (s/m)
      ! upper threshold from Kerkweg et al. 2006
      rcutd_o3 = MIN(1.e5_dp,(rcutd0/(exp(0.03_dp*prh_2m(jl)&      
                    *100._dp)*pplai**(0.25_dp)*ppustveg)))         
      ! cuticular resistance of ozone above wet canopy (s/m)
      rcutw_o3 = rcutw0/MAX(1.e-10_dp,(pplai**(0.5_dp)      &      
                 *ppustveg))
      ! tracer-dependent cuticular resistance above dry canopy (s/m)
      rcut_2d(jl) = rcutd_o3/MAX(1.e-15_dp,ftrac)                  
      ! tracer-dependent cuticular resistance above wet canopy (s/m)
      rws_2d(jl) = 1._dp/MAX(1.e-10_dp,(1._dp/(3._dp*rcutw_so2) &  
                   +1.e-7_dp*henry +dryreac/rcutw_o3)) 
    END DO

  END SUBROUTINE drydep_calc_rcut
! ju_te_20180629-
!=============================================================================

!=============================================================================

SUBROUTINE drydep_calc_ra(kproma,  loland, &
       prahl,   prahw,    prahi,   prahveg, prahslsn,       &
       pustarl, pustarw,  pustari, pustveg, pustslsn,       &
       psurf, pcfml, pcfncl,pcdnl, sqrtwind, ptvir, ptvl,   &
       pril, paz0, pz0m, pcdnw, pcfmw, pcfncw, priw, ptvw , &
       pcdni, pcfmi, pcfnci, prii, ptvi )

  USE messy_main_constants_mem, ONLY: g

    IMPLICIT NONE

    INTEGER :: kproma

    ! external parameters

    LOGICAL, DIMENSION(kproma), INTENT(in)  :: loland

    ! internal parameters

    REAL(dp), DIMENSION(kproma), INTENT(in)  :: ptvl
    REAL(dp), DIMENSION(kproma), INTENT(in)  :: ptvir
    REAL(dp), DIMENSION(kproma), INTENT(in)  :: ptvw
    REAL(dp), DIMENSION(kproma), INTENT(in)  :: ptvi
    REAL(dp), DIMENSION(kproma), INTENT(in)  :: pcfml
    REAL(dp), DIMENSION(kproma), INTENT(in)  :: pcfncl
    REAL(dp), DIMENSION(kproma), INTENT(in)  :: pcdnl
    REAL(dp), DIMENSION(kproma), INTENT(in)  :: pcdnw
    REAL(dp), DIMENSION(kproma), INTENT(in)  :: pcfmw
    REAL(dp), DIMENSION(kproma), INTENT(in)  :: pcfncw
    REAL(dp), DIMENSION(kproma), INTENT(in)  :: pcdni
    REAL(dp), DIMENSION(kproma), INTENT(in)  :: pcfmi
    REAL(dp), DIMENSION(kproma), INTENT(in)  :: pcfnci
    REAL(dp) ::  zcml(1:kproma),                  &
         cmw(1:kproma),cmi(1:kproma),cmveg(1:kproma),&
         cmslsn(1:kproma)
    REAL(dp), DIMENSION(kproma) :: pustarl,pustarw,pustari,pustveg,pustslsn
    REAL(dp), DIMENSION(kproma) :: prahl,prahw,prahi,prahveg,prahslsn
    REAL(dp) ::  zcdnveg(1:kproma),zcdnslsn(1:kproma), psih(1:kproma)

    REAL(dp) :: zoverl(1:kproma),zxzsurf(1:kproma),&
         zxzref(1:kproma)

    REAL(dp), DIMENSION(kproma) :: psurf
    REAL(dp)                    :: zmonin(1:kproma)
    REAL(dp), DIMENSION(kproma) :: sqrtwind
    REAL(dp), DIMENSION(kproma) :: pril, priw, prii
    REAL(dp), DIMENSION(kproma) :: paz0
    REAL(dp), DIMENSION(kproma) :: pz0m

    REAL(dp), PARAMETER :: zkmkh=0.74_dp 
    REAL(sp), PARAMETER :: ckap = 0.4 ! von Karman constant

    ! calculation of drag coefficient and u*. A change with
    ! the previous calculations in echam4 is that the exchange
    ! coefficients, calculated in echam5 over land, ice and water,
    ! is being used to explicitly calculate the aerodynamic
    ! resistances over these surface

    ! ==================================================================
    ! OVER LAND
    ! ==================================================================

    WHERE (ABS(pcfncl(:)) >  TINY(pcfncl(1)))
       zcml(:)=pcdnl(:)*pcfml(:)/pcfncl(:)
    ELSEWHERE
       zcml(:)=0.
    ENDWHERE

    pustarl(:)=SQRT(zcml(:))*sqrtwind(:)

    WHERE (pril(:).gt.0._dp)

       ! calculating the Monin-Obukhov lenght directly applying the
       ! formula given by Stull, 9.7.5k, page 386
       
       zmonin(:)=(pustarl(:)*((ptvir(:)+ptvl(:))/2._dp)*  &
            sqrtwind(:))/(ckap*g*(ptvir(:)- ptvl(:))) 
       zoverl(:)=psurf(:)/zmonin(:)
       psih(:)=-4.7_dp*zoverl(:)

    ELSEWHERE

       zmonin(:)=psurf(:)/pril(:)
       zoverl(:)=psurf(:)/zmonin(:)
       zxzsurf(:)=zkmkh*(1._dp-9.*(zoverl(:)))**(0.5)
       zxzref(:)=zkmkh
       psih(:)= &
            (2.*LOG((1._dp+zxzsurf(:))/2.)+LOG((1._dp+zxzsurf(:)**2.)/2.)-  &
            2.*ATAN(zxzsurf(:)))-  & ! primitive function value for z
            (2.*LOG((1._dp+zxzref(:))/2.)+LOG((1._dp+zxzref(:)**2.)/2.)- &
            2.*ATAN(zxzref(:)))      ! primitive function value for zz
       
    END WHERE
    
    prahl(:)= MAX(1._dp,(1./(pustarl(:)*ckap))* &
         (LOG(psurf(:)/paz0(:))-psih(:)))

    WHERE (pz0m(:).gt.0._dp .or.loland(:))
       zcdnveg(:)=(ckap/LOG(1._dp+psurf(:)/(MAX(0.02_dp,pz0m(:)))))**2
       zcdnslsn(:)=(ckap/LOG(1._dp+psurf(:)/(0.005)))**2
       cmveg(:)=zcdnveg(:)*pcfml(:)/pcfncl(:)
       cmslsn(:)=zcdnslsn(:)*pcfml(:)/pcfncl(:)
       pustveg(:)=SQRT(cmveg(:))*sqrtwind(:)
       pustslsn(:)=SQRT(cmslsn(:))*sqrtwind(:)
       prahveg(:)=MAX(1._dp,(1./(pustveg(:)*ckap))* &
               (LOG((psurf(:))/ MAX(0.02_dp,pz0m(:)))-psih(:)))
       prahslsn(:)=MAX(1._dp,(1./(pustslsn(:)*ckap))* &
               (LOG(psurf(:)/ 0.005)-psih(:)))
    ELSEWHERE
       cmveg(:)=zcml(:)
       cmslsn(:)=zcml(:)
       pustveg(:)=pustarl(:)
       pustslsn(:)=pustarl(:)
       prahveg(:)=prahl(:)
       prahslsn(:)=prahl(:)
     END WHERE

    ! ==================================================================
    ! OVER WATER
    ! ==================================================================
     
     WHERE (ABS(pcfncw(:))> TINY(pcfncw(1)))
        cmw(:)=pcdnw(:)*pcfmw(:)/pcfncw(:)
     ELSEWHERE
        cmw(:)=0.
     ENDWHERE

     pustarw(:)=SQRT(cmw(:))*sqrtwind(:)

     ! Computation of stability correction term, The stability
     ! correction functions are taken from Stull (page 383-385)
     ! (08-11-98) and are slightly different from those by Williams
     ! and Hicks et al., which were originaly being used in the dry
     ! deposition scheme.

     WHERE (priw(:).gt.0._dp) 
        
        ! calculating the Monin-Obukhov lenght directly applying the
        ! formula given by Stull, 9.7.5k, page 386

          zmonin(:)=(pustarw(:)*((ptvir(:)+ ptvw(:))/2.)*             &
               sqrtwind(:))/(ckap*g*(ptvir(:)-ptvw(:)))
          zoverl(:)=psurf(:)/zmonin(:)
          psih(:)=-4.7*zoverl(:)
       ELSEWHERE
          zmonin(:)=psurf(:)/priw(:)
          zoverl(:)=psurf(:)/zmonin(:)
          zxzsurf(:)=zkmkh*(1._dp-9.*(zoverl(:)))**(0.5)
          zxzref(:)=zkmkh
          psih(:)= &
               (2.*LOG((1._dp+zxzsurf(:))/2.)+LOG((1._dp+zxzsurf(:)**2.)/2.)-  &
               2.*ATAN(zxzsurf(:)))-  & ! primitive function value for z
               (2.*LOG((1._dp+zxzref(:))/2.)+LOG((1._dp+zxzref(:)**2.)/2.)- &
               2.*ATAN(zxzref(:)))      ! primitive function value for zz
       ENDWHERE

       prahw(:)=MAX(1._dp,(1./(pustarw(:)*ckap))* &
            (LOG((psurf(:))/paz0(:))-psih(:)))

       ! ==================================================================
       ! OVER ICE
       ! ==================================================================
       

       WHERE (ABS(pcfnci(:)) > TINY(pcfnci(1))) 
          cmi(:)=pcdni(:)*pcfmi(:)/pcfnci(:)
       ELSEWHERE
          cmi(:) = 0._dp
       ENDWHERE

       pustari(:)=SQRT(cmi(:))*sqrtwind(:)

       ! Computation of stability correction term, The stability
       ! correction functions are taken from Stull (page 383-385)
       ! (08-11-98) and are slightly different from those by Williams
       ! and Hicks et al., which were originaly being used in the dry
       ! deposition scheme.

       WHERE (prii(:).gt.0._dp) 

          ! calculating the Monin-Obukhov lenght directly applying the
          ! formula given by Stull, 9.7.5k, page 386
          
          zmonin(:)=(pustari(:)*((ptvir(:)+ptvi(:))/2.)*             &
               sqrtwind(:))/(ckap*g*(ptvir(:)-ptvi(:)))
          zoverl(:)=psurf(:)/zmonin(:)
          psih(:)=-4.7*zoverl(:)

       ELSEWHERE

          zmonin(:)=psurf(:)/prii(:)
          zoverl(:)=psurf(:)/zmonin(:)
          zxzsurf(:)=zkmkh*(1._dp-9.*(zoverl(:)))**(0.5)
          zxzref(:)=zkmkh
          psih(:)= &
               (2.*LOG((1._dp+zxzsurf(:))/2.)+LOG((1._dp+zxzsurf(:)**2.)/2.)-  &
               2.*ATAN(zxzsurf(:)))-  & ! primitive function value for z
               (2.*LOG((1._dp+zxzref(:))/2.)+LOG((1._dp+zxzref(:)**2.)/2.)- &
               2.*ATAN(zxzref(:)))      ! primitive function value for zz
       ENDWHERE

       prahi(:)=MAX(1._dp,(1./(pustari(:)*ckap))* &
            (LOG((psurf(:))/paz0(:))-psih(:)))


END SUBROUTINE drydep_calc_ra

!=============================================================================


  SUBROUTINE drydep_vdbl_parameter(kproma,  &
        ustveg, ptslm1, prh_2m, plai, phc, psoilph)
    !=====================================================================
    !    Program to calculate the dry deposition velocity for
    !    a selection of trace gases and the sulfate aerosol
    !    according to the "big leaf" approach" using the ECHAM
    !    surface cover characterization and surface parameters
    !    See for more information about the code the papers by
    !    Ganzeveld and Lelieveld, JGR 100, 1995,
    !    Ganzeveld et al., JGR 103, 1998
    !    Ganzeveld et al., 2001
    !=====================================================================

    ! code converted to f90 by Rolf Sander, April 2001, further
    ! modified by Laurens Ganzeveld for implementation in ECHAM5
    ! October 2001

    IMPLICIT NONE

    INTEGER :: kproma

    REAL(dp), DIMENSION(:)   :: ptslm1, prh_2m
    REAL(dp), DIMENSION(:)   :: phc, plai
    REAL(dp), DIMENSION(:,:) :: psoilph 
    REAL(dp), DIMENSION(:)   :: ustveg

    ALLOCATE(rahcan(1:kproma))
    ALLOCATE(rsnowhno3so2(1:kproma))
    ALLOCATE(rsoilso2(1:kproma))

    ! calculation of standard deviation of wind direction, used
    ! for the calculation of the SO4-- deposition velocity, this should
    ! be replaced once by using the Richardson number as criteria to
    ! select between two stability classes relevant to the SO4-- deposition
    ! calculations.

    ! from now on using LAI for calculation of within 
    ! canopy aerodynamic resistance

    rahcan(:)=MAX(1._dp,14._dp*plai(:)*phc(:)/ ustveg(:))
    
    ! calculation of the HNO3 and SO2 snow resistance as a
    ! function of the snow temperature

    rsnowhno3so2(:)=MIN(MAX(10._dp,10._dp**(-0.09*(ptslm1(:)-273._dp) &
         +2.4_dp)),1.e5_dp)

    ! calculation of soil resistance of SO2 as a function of the RH
    ! and the soil pH. There are five different soil pH classes and
    ! each class has a assigned surface resistance, which is subsequently
    ! corrected for the RH

    !   the soil pH classes are:
    !   SOILPH(JL,1,JR) - soil pH <5.5
    !   SOILPH(JL,2,JR) - soil 5.5 <pH <7.3
    !   SOILPH(JL,3,JR) - soil 7.3< pH <8.5
    !   SOILPH(JL,4,JR) - soil 8.5 <pH
    !   SOILPH(JL,5,JR) - soil 4 < pH <8.5
    
       rsoilso2(:)=MAX(50._dp, &
            psoilph(:,1)*115._dp+ &
            psoilph(:,2)*65._dp+ &
            psoilph(:,3)*25._dp+ &
            psoilph(:,4)*25._dp+ &
            psoilph(:,5)*70._dp+ &
            MIN(MAX(0._dp,1000._dp*EXP(MIN(10._dp,269._dp-ptslm1(:)))), &
            1.e5_dp))
       WHERE (prh_2m(:).LT.0.6_dp.AND.prh_2m(:).GT.0.4_dp) ! semi-arid regions
          rsoilso2(:)=MAX(50._dp,(rsoilso2(:)*3.41_dp-85._dp)+ &
               MIN(MAX(0._dp,1000.*EXP(MIN(10._dp,269._dp-ptslm1(:)))), &
               1.e5_dp))
       ENDWHERE
       WHERE (prh_2m(:).LE.0.4_dp)       ! arid regions
          rsoilso2(:)=MAX(50._dp,MIN(1000._dp,(rsoilso2(:)*3.41-85._dp+ &
               ((0.4_dp-prh_2m(:))/0.4_dp)*1.e5_dp)+ &
               MIN(MAX(0._dp,1000.*EXP(MIN(10._dp,269._dp-ptslm1(:)))), &
               1.e5_dp)))
       ENDWHERE
       
     END SUBROUTINE drydep_vdbl_parameter

!============================================================================

 SUBROUTINE drydep_vdbl(kproma,  tracnum, laland ,      &
      pvdrydep, trname, itr,                            &
       pfri,   pcvs,  pcvw,    pvgrat,    rahw,   rahi, &
       rahveg,  rahslsn,   ustarw, ustari,              &
       ustveg,  ustslsn, plai,prco_leaf, pfws,          &
       diff, diffrb, rmes, rsoil, rwater, rsnow, rws,   & ! ju_te_20180625 
       rcut, rcut_2d, rws_2d)                             ! ju_te_20180625 

!=====================================================================
   !    Program to calculate the dry deposition velocity for
   !    a selection of trace gases and the sulfate aerosol
   !    according to the "big leaf" approach" using the ECHAM
   !    surface cover characterization and surface parameters
   !    See for more information about the code the papers by
   !    Ganzeveld and Lelieveld, JGR 100, 1995,
   !    Ganzeveld et al., JGR 103, 1998
   !    Ganzeveld et al., 2001
   !    Tamara Emmerichs, FZ juelich, 2018-06-26:
   !    calculation of cuticular resistance
   !    depending on canopy wetness, LAI, friction velocity, RH
   !    can be switch on with l_ganzeori=.FALSE. (in ddep.nml)
   !    For more information about the calculation see
   !    Zhang et al., Atm. Environment 36 (2002) 4787-4799;
   !    Zhang et al., Atmos. Chem. Phys., 3, 2067-2082, 2003
   !=====================================================================
  
   ! code converted to f90 by Rolf Sander, April 2001, further
   ! modified by Laurens Ganzeveld for implementation in ECHAM5
   ! October 2001
   
   IMPLICIT NONE
   
   INTEGER,          INTENT(IN) :: kproma
   CHARACTER(len=*), INTENT(IN) :: trname
   ! number of actual tracer index in interface
   INTEGER,          INTENT(IN) :: itr   

   INTEGER,  INTENT(IN)                        :: tracnum
   REAL(dp), INTENT(INOUT), DIMENSION(tracnum) :: &
         diff, diffrb, rmes, rsoil, rwater, rsnow
 
   REAL(dp), DIMENSION(1:kproma),  INTENT(IN)  :: plai
   REAL(dp), DIMENSION(1:kproma),  INTENT(IN)  :: prco_leaf, pfws

   LOGICAL,  DIMENSION(1:kproma),  INTENT(IN)  :: laland

   REAL(dp), DIMENSION(1:kproma),  INTENT(IN)  :: pcvs,pcvw,pvgrat,pfri

   REAL(dp), DIMENSION(1:kproma),  INTENT(OUT) :: pvdrydep
   ! ju_te_20180625+
   REAL(dp), DIMENSION(1:kproma), OPTIONAL, INTENT(IN) :: rws_2d, rcut_2d
   REAL(dp), INTENT(IN), OPTIONAL, DIMENSION(tracnum)  :: rcut, rws
   ! ju_te_20180629-

   ! LOCAL
   REAL(dp) :: rbw(1:kproma),rbi(1:kproma),          &
         rbveg(1:kproma),     rbslsn(1:kproma),      &
         rc0x(1:kproma),      rsveg(1:kproma),       &
         rleaf(1:kproma),     vdveg(1:kproma),       &
         vdsoil(1:kproma),    vdsn(1:kproma),        &
         vdland(1:kproma),    vdwat(1:kproma),       &
         vdws(1:kproma)
   
   REAL(dp), DIMENSION(:) :: rahw,rahi,rahveg,rahslsn, &
        ustarw,ustari,ustveg,ustslsn

   REAL(dp)             :: zrmes(1:kproma)
   REAL(dp)             :: fsic(1:kproma),fsicm1(1:kproma)
   
   REAL(dp), DIMENSION(1:kproma) :: rc0x_o3

    INTRINSIC TRIM, MAX, MIN

    ! correcting for soil moisture
    rc0x_o3(:) = diff_o3 *prco_leaf(:)/ MAX(1.e-5_dp,pfws(:)) 
   
    ! calculation of quasi-laminar boundary layer resistances
    rbveg(:)=(2./(ustveg(:)*0.40))*diffrb(itr)
    rbslsn(:)=(1./(ustslsn(:)*0.40))*diffrb(itr)
    rbw(:)=(1./(ustarw(:)*0.40))*diffrb(itr)
    rbi(:)=(1./(ustari(:)*0.40))*diffrb(itr)
    
    ! Stomatal resistance computation for each species, by
    ! multiplying the leaf stomatal resistance with the
    ! diffusivity term
    
    ! correcting for soil moisture
    rc0x(:)=diff(itr)*prco_leaf(:)/MAX(1.e-5_dp,pfws(:)) 
    
    ! definition of the NO and NO2 mesophyllic resistance as a function
    ! of the ozone mesophyllic resistances (see Ganzeveld and Lelieveld, '95)
    
    SELECT CASE(TRIM(trname))
    CASE('NO')
       zrmes(:)=5._dp*rc0x_o3(:)
    CASE('NO2')
       zrmes(:)=0.5_dp*rc0x_o3(:)
    CASE DEFAULT
       zrmes(:) = rmes(itr)
    END SELECT

    IF (.NOT.l_ganzeori) THEN  ! ju_te_20180625+
       WHERE(rc0x(:)+zrmes(:) .gt. 0.0_dp .and. rcut_2d(:) .gt. 0.0_dp)
          rleaf(:)=(1._dp/((1._dp/rcut_2d(:))+(1._dp/(rc0x(:)+zrmes(:)))))
       ELSEWHERE
          rleaf(:) = 0._dp
       END WHERE
    ELSE                       ! um_ak_20130926-
       WHERE(rc0x(:)+zrmes(:) > 0.0_dp)
          rleaf(:)=(1./((1./rcut(itr))+(1./(rc0x(:)+zrmes(:)))))
       ELSEWHERE
          rleaf(:) = 0._dp
       END WHERE
    END IF ! ju_te_20180625

    IF (l_ganzeori) THEN ! ju_te_20190503
       WHERE (rleaf(:) > 0._dp)
          rsveg(:)=(1./((1./(rahcan(:)+ &
               (1./(ustveg(:)*0.40))*diffrb(itr)+rsoil(itr)))+ &
               (1./(rleaf(:)/MAX(1.e-5_dp,plai(:))))))
       ELSEWHERE
          rsveg(:) = 0._dp
       END WHERE
    ELSE
       WHERE (rleaf(:) > 0._dp)
          ! 1./(ustveg(:)*0.40))*diffrb(itr) - rbveg for individiual leasves 
          ! is considered by R_cut through parametrizing it as a function of 
          ! friction velocity (Zhang et al. 2003), rleaf given at canopy scale
          rsveg(:)=(1./((1./(rahcan(:)+rsoil(itr)))+ (1./(rleaf(:))))) 
       ELSEWHERE
          rsveg(:)=0
       END WHERE
    END IF

    SELECT CASE(TRIM(trname))
    CASE('HNO3', 'N2O5')
      rsveg(:)=MAX(10._dp,rsveg(:))
      vdveg(:)=(1./(rahveg(:)+rbveg(:)+rsveg(:)))
      vdsoil(:)=(1./(rahslsn(:)+rbslsn(:)+rsoil(itr)))
      IF (.NOT.l_ganzeori) THEN ! ju_te_20180629+
         vdws(:)=1._dp/(rahveg(:)+rbveg(:)+rws_2d(:))
      ELSE                      ! ju_te_20180629-
         vdws(:)=(1./(rahveg(:)+rbveg(:)+rws(itr)))
      END IF   ! ju_te_20180629
      vdsn(:)= &
        (1./(rahslsn(:)+rbslsn(:)+rsnowhno3so2(:)))

    CASE('SO2','H2SO4')
      vdveg(:)=(1./(rahveg(:)+rbveg(:)+rsveg(:)))
      vdsoil(:)=(1./(rahslsn(:)+rbslsn(:) &
        +rsoilso2(:)))
      IF (.NOT.l_ganzeori) THEN ! ju_te_20180629+
         vdws(:)=1._dp/(rahveg(:)+rbveg(:)+rws_2d(:))
      ELSE                      ! ju_te_20180629-
         vdws(:)=(1./(rahveg(:)+rbveg(:)+rws(itr)))
      END IF ! ju_te_20180629
      vdsn(:)= &
        (1./(rahslsn(:)+rbslsn(:)+rsnowhno3so2(:)))
    CASE DEFAULT 
      vdveg(:)=(1./(rahveg(:)+rbveg(:)+rsveg(:)))
      vdsoil(:)=(1./(rahslsn(:)+rbslsn(:)+rsoil(itr)))
      IF (.NOT.l_ganzeori) THEN ! ju_te_20180629+
         vdws(:)=1._dp/(rahveg(:)+rbveg(:)+rws_2d(:))
      ELSE                      ! ju_te_20180629-
         vdws(:)=(1./(rahveg(:)+rbveg(:)+rws(itr)))
      END IF   ! ju_te_20180629
      vdsn(:)=(1./(rahslsn(:)+rbslsn(:)+rsnow(itr)))
    END SELECT

    ! calculation of the over-land dry deposition velocity, considering
    ! the surface cover fractions
    
    vdland(:)=pcvs(:)*vdsn(:)+ &
         (1._dp-pcvs(:))*(1._dp-pcvw(:))*pvgrat(:)*vdveg(:)+ &
         (1._dp-pcvs(:))*(1._dp-pvgrat(:))*(1._dp-pcvw(:))*vdsoil(:)+ &
         (1._dp-pcvs(:))*pcvw(:)*vdws(:)
    
    ! There is a modification in echam5 compared to echam4
    ! with pfri > 0 over land for grids with frozen lake surfaces.
    ! The value is not contineous but just 0 or 1 over land (to avoid
    ! the introduction of serious complications in some of the physical
    ! routines, personal communications, Uli Schlese, 13-02-2002).
    ! This change requires a modification of the use of the parameter
    ! seaice for use in the surface trace gas and aerosol exchange
    ! processes.
    
    fsic(:)=MIN(1._dp,pfri(:))
    fsicm1(:)=1._dp-fsic(:)

    SELECT CASE(TRIM(trname))
    CASE('HNO3', 'N2O5')
       vdwat(:)=fsic(:)*(1./(rahi(:) &        ! different formula for hno3
            +rbi(:)+rsnowhno3so2(:)))+ &
            fsicm1(:)*(1./(rahw(:)+rbw(:)+rwater(itr)))
    CASE('SO2', 'H2SO4')
       vdwat(:)=fsic(:)*(1./(rahi(:) &        ! different formula for so2
            +rbi(:)+rsnowhno3so2(:)))+ &
            fsicm1(:)*(1./(rahw(:)+rbw(:)+rwater(itr)))
    CASE DEFAULT
       vdwat(:)=fsic(:)*(1./(rahi(:)+rbi(:)+rsnow(itr)))+ &
            fsicm1(:)*(1./(rahw(:)+rbw(:)+rwater(itr)))
    END SELECT

    ! finally, calculation of the grid deposition velocity from the
    ! over-land and ocean dry deposition velocities, the deposition
    ! velocities in pvdrydep are given in [m s-1]!
    
    WHERE(laland(:))
       pvdrydep(:)=vdland(:)
    ELSEWHERE
       pvdrydep(:)=vdwat(:)
    END WHERE
    pvdrydep(:)=MAX(1.e-5_dp,pvdrydep(:))
    
  END SUBROUTINE drydep_vdbl

!==========================================================================

  SUBROUTINE drydep_vdbl_dealloc
    
    DEALLOCATE(rahcan)
    DEALLOCATE(rsnowhno3so2)
    DEALLOCATE(rsoilso2)

  END SUBROUTINE drydep_vdbl_dealloc
!=============================================================================

!=============================================================================

  SUBROUTINE drydep_vdaer_parameter( kproma,   &
       loland,            pustarw, &
       abswind, abswind10, prh_2m, &
       palpha, palphae, pbubble, pbeta)

    !=====================================================================
    !    Program to calculate the dry deposition velocity for
    !    aerosols as a function of the particle radius considering the
    !    impaction and difussion. Sedimentation is also calculated but not
    !    included in the calculated deposition velocity since this process
    !    is considered in a separate routine
    !=====================================================================

    !   Author:
    !   -------
    !   Laurens Ganzeveld, MPI Mainz                          04/2001
    !
    !   Modifications:
    !   --------------
    !   Philip Stier,      MPI Hamburg (adoption to ECHAM/M7) 12/2001
    !   Swen Metzger,      MPI Mainz                          02/2004
    !                     (generalization - independance of M7)
    !   Astrid Kerkweg,    MPI Mainz, modularization,         04/2004 
    !                                 MESSy-conformity   
    !
    !   Methodology:
    !
    !     This program calculates the integrated deposition velocity from
    !     the mass size distribution for a log normal aerosol distribution
    !     defined by the radius and the log sigma. The model calculates the
    !     Vd over land and over sea considering the diffusion,impaction and
    !     sedimentation. Over sea the effect op particle growth due to the
    !     large relative humidity is accounted for and the effect of
    !     bubble bursting is also considered. The bubble bursting enhances
    !     the dry deposition since it causes the breakdown of the laminar
    !     boundary layer and the scavenging of the particles by the sea spray.
    !     This model version does not contain yet a parameterization which
    !     specifically considers the deposition to vegetated surfaces. For
    !     these surfaces, a surface resistance as a function of canopy
    !     structure should be incorporated. This might possible be implemented
    !     in the future (Laurens Ganzeveld, 29-01-2002)
    !
    !     The model requires as input parameters:
    !     -------------------------------------------------------------------
    !     - UM, windspeed at reference height
    !     - PAZ0M, surface roughness, over sea this term is calculated from UM
    !     - PTSLM1, Air or surface temperature
    !     - LOLAND, land-sea mask
    !     - N, aerosol number [cm -3]
    !     - R, aerosol radius [m]
    !     - LSIGMA, LOG sigma of the log normal distribution
    !     - RHOA_AEROS, the density of the aerosol   
    !          @@@ Replaced by: densaer
    ! =======================================================================

    USE messy_main_constants_mem, ONLY: api => pi

    IMPLICIT NONE

    ! LOCAL
    !CHARACTER(LEN=*),  PARAMETER :: substr='drydep_vdaer'

    INTEGER :: kproma

    LOGICAL :: loland(:)

    REAL(dp)    :: pustarw(:)

    REAL(dp) :: abswind(:), abswind10(:)
    REAL(dp) :: prh_2m(:)        
    ! LOCAL

    REAL(dp) :: palpha(:) , pbubble(:), pbeta(:), palphae(:)
    
    REAL(dp) :: s(1:kproma), phi(1:kproma)
    REAL(dp) :: alpha1(1:kproma), alpharat(1:kproma)
    REAL(dp) :: vk1(1:kproma), vk2(1:kproma)
    
    REAL(dp) :: eff, rdrop, zdrop
    REAL(dp) :: qdrop(1:kproma)
    REAL(dp) :: eps


    INTRINSIC MAX, MIN, EXP
    !--- 0) Initializations:

    !--- 1) Calculate correction factors:

    ! calculation of some parameters required for the deposition calculations
    ! bubble bursting effect,see Hummelshoj, equation 10
    ! relationship by Wu (1988), note that Hummelshoj has not
    ! considered the cunningham factor which yields a different
    ! vb curve, with smaller values for small particles
    
    palpha(:)=MAX(1.e-10_dp,1.7e-6*abswind10(:)**3.75)   
    eff=0.5_dp
    rdrop=0.005_dp      ! cm
    zdrop=10.0_dp       ! cm
    qdrop=5.*(100.*palpha(:))
    IF (l_whitecap) THEN
       pbubble(:)=((100.*pustarw(:))**2.)/(100.*abswind(:))+eff* & 
            !! mz_ak_20041025 ustarl-> ustarw
            (2.*api*rdrop**2.)*(2.*zdrop)*(qdrop(:)/palpha(:))
    ELSE 
       pbubble(:)=0._dp
       palpha(:) =0._dp
    ENDIF
    
    !--- Correction of particle radius for particle growth close to the
    !    surface according to Fitzgerald, 1975, the relative humidity over
    !    the ocean is restricted to 98% (0.98) due to the salinity

    s(:)        = MIN(0.98_dp,prh_2m(:))
    eps         = 0.6_dp
    pbeta(:)    = EXP((0.00077_dp*s(:))/(1.009_dp-s(:)))
    phi(:)      = 1.058_dp-((0.0155_dp*(s(:)-0.97_dp))/(1.02_dp-s(:)**1.4))
    alpha1(:)   = 1.2_dp*EXP((0.066_dp*s(:))/(phi(:)-s(:)))
    vk1(:)      = 10.2_dp-23.7*s(:)+14.5*s(:)**2.
    vk2(:)      = -6.7_dp+15.5*s(:)-9.2*s(:)**2.
    alpharat(:) = 1._dp-vk1(:)*(1._dp-eps)-vk2(:)*(1._dp-eps**2)
    palphae(:)  = alpharat(:)*alpha1(:)
    
    !--- Over land no correction for large humidity close to the surface:

    WHERE (loland(:).OR..NOT.l_rh)
       palphae(:)=1._dp
       pbeta(:)=1._dp
    ENDWHERE
    
  END SUBROUTINE drydep_vdaer_parameter

!============================================================================

  SUBROUTINE drydep_vdaer( pvdrydep,  kproma,             &
         pvgrat,   pcvs,      pcvw,     pfri,     pcvbs,  &
         pcvwat,   pustarw,   pustveg,  pustslsn, prahw,   &
         prahveg,  prahslsn,  pevap,    paz0m,            &
         abswind,  pradius,   pdensaer, ptslm1,           &
         palpha,   palphae,   pbubble,  pbeta             )

    !=====================================================================
    !    Program to calculate the dry deposition velocity for
    !    aerosols as a function of the particle radius considering the
    !    impaction and diffusion.
    !=====================================================================

    !   Author:
    !   -------
    !   Laurens Ganzeveld, MPI Mainz                          04/2001
    !
    !   Modifications:
    !   --------------
    !   Philip Stier,      MPI Hamburg (adoption to ECHAM/M7) 12/2001
    !   Swen Metzger,      MPI Mainz                          02/2004
    !                     (generalization - independance of M7)
    !   Astrid Kerkweg,    MPI Mainz, modularization,         04/2004 
    !                                 MESSy-conformity   
    !
    !   Methodology:
    !
    !     This program calculates the integrated deposition velocity from
    !     the mass size distribution for a log normal aerosol distribution
    !     defined by the radius and the log sigma. The model calculates the
    !     Vd over land and over sea considering the diffusion,impaction and
    !     sedimentation. Over sea the effect of particle growth due to the
    !     large relative humidity is accounted for and the effect of
    !     bubble bursting is also considered. The bubble bursting enhances
    !     the dry deposition since it causes the breakdown of the laminar
    !     boundary layer and the scavenging of the particles by the sea spray.
    !     This model version does not contain yet a parameterization which
    !     specifically considers the deposition to vegetated surfaces.
    !
    !     The model requires as input parameters:
    !     -------------------------------------------------------------------
    !     - UM, windspeed at reference height
    !     - PAZ0M, surface roughness, over sea this term is calculated from UM
    !     - PTSLM1, Air or surface temperature
    !     - LOLAND, land-sea mask
    !     - N, aerosol number [cm -3]
    !     - R, aerosol radius [m]
    !     - LSIGMA, LOG sigma of the log normal distribution
    !     - RHOA_AEROS, the density of the aerosol   
    !          @@@ Replaced by: densaer
    !========================================================================

    USE messy_main_constants_mem, ONLY:  pi, g

    IMPLICIT NONE

    !--- input parameters from echam, the parameter pevap has still to be
    !    checked more carefully, if this is a used constant or actually the
    !    echam evaporation

    ! LOCAL
    !CHARACTER(LEN=*),  PARAMETER :: substr='drydep_vdaer_parameter'

    INTEGER :: kproma

    REAL(dp), INTENT(out) :: pvdrydep(:) 
    
    REAL(dp)    :: pustarw(:),  pustveg(:),  pustslsn(:),  &
               prahw(:),     prahslsn(:), pevap(:), prahveg(:)
    REAL(dp)   ::  pvgrat(:),  pcvs(:),  pcvw(:),   &
               pcvbs(:),     pcvwat(:),  pfri(:)   

    REAL(dp) :: abswind(:)
    REAL(dp) ::  paz0m(:)
    REAL(dp) :: ptslm1(:), pdensaer(:)

    ! LOCAL

    REAL(dp) :: palpha(:) , pbubble(:), palphae(:), pbeta(:)
    
    REAL(dp) :: pradius(:)
    REAL(dp) :: cunning(1:kproma), dc(1:kproma), relax(1:kproma)
    REAL(dp) :: sc(1:kproma), st_veg(1:kproma)
    REAL(dp) :: st_slsn(1:kproma), st_wat(1:kproma)
    REAL(dp) ::  vb_veg(1:kproma),vb_slsn(1:kproma)
    REAL(dp) :: vim_veg(1:kproma), vim_slsn(1:kproma)
    REAL(dp) :: vkd_veg(1:kproma), vin_veg(1:kproma)
    REAL(dp) :: zrebound(1:kproma), vkd_slsn(1:kproma)
    REAL(dp) :: re(1:kproma) , vkd_wat(1:kproma)
    REAL(dp) :: vb_sea(1:kproma), vim_sea(1:kproma)
    REAL(dp) :: vkdaccsea(1:kproma)
    REAL(dp) :: vkc_veg(1:kproma), vkc_slsn(1:kproma)
    REAL(dp) :: vkc_wat(1:kproma), vdpart_veg(1:kproma)
    REAL(dp) :: vdpart_slsn(1:kproma), vdpart_wat(1:kproma)
    

    REAL, PARAMETER :: daccm=0.         ! factor which corrects for 
                                        !  evaporation (see paper slinn)
    REAL, PARAMETER :: visc=0.15        ! molecular viscocity [cm2 s-1]
    REAL, PARAMETER :: dynvisc=1.789e-4 ! g cm-1 s-1
    REAL, PARAMETER :: cl=0.066*1e-4    ! mean free path [cm] (particle size 
                                        ! also in cm)
    REAL, PARAMETER :: bc= 1.38e-16     ! boltzman constant [g cm2 s-2 K-1] 
                                        ! (1.38e-23 J deg-1)
    REAL, PARAMETER :: kappa=1.         ! shapefactor

    REAL, PARAMETER :: ZAS    = 10.E-6*1.E2  ! um -> CM, see paper
                                             !Gallagher and Slinn, 1982
    REAL, PARAMETER :: zAL    = 1.E-3*1.E2   ! mm -> CM, see paper
                                             ! Gallagher and Slinn, 1982

    REAL(dp), PARAMETER :: vkar=0.40_dp       ! von karman constant

   REAL(dp), PARAMETER :: crmin=0.01E-6_dp ! smallest radius for which s dry 
                                    ! deposition velocity is calculated

   INTEGER :: jl

    !---  Calculation of the dry deposition velocity utilizing 
    !     the log-normal properties:

    !--- Do calculations only for particles larger than 0.0001 um :

    DO jl =1, kproma
       pevap (jl) = 0._dp  ! makes no sense but where does the value come from ?

       IF (pradius(jl) > crmin) THEN

          !--- Cunningham factor:
       
          IF ( ABS(palphae(jl)) > TINY(palphae(1))) THEN
          
             cunning(jl)=1._dp+(cl/(palphae(jl)*pradius(jl))**pbeta(jl))*      &
                  (1.257_dp+0.400*EXP(-1.1_dp*(pradius(jl))/cl))
             
             !--- Diffusivity:
             
             dc(jl)=(bc*ptslm1(jl)*cunning(jl))/(3.*pi*dynvisc*  &
                  (palphae(jl)*pradius(jl))**pbeta(jl))

       ELSE
          cunning(jl) = 0.
          dc(jl) = 0.
       ENDIF

       ! Relaxation:

       relax(jl)=(pdensaer(jl)*1.e-3_dp*&
            (((palphae(jl)*pradius(jl))**pbeta(jl))**2.)* &
            cunning(jl))/(18._dp*dynvisc*kappa)

       ! Sedimentation is calculated operator split in the 
       ! amodule handeling sedimentation:

       ! Calculation of schmidt and stokes number

       IF (ABS(dc(jl)) > TINY(dc(jl))) THEN
          sc(jl)      =visc/dc(jl)
       ELSE
          sc(jl) = 0.
       ENDIF
       !  modified calculation of the Stokes number
       !  over vegetated surfaces, see paper by Gallagher et al., 
       !  JGR 2002. zAL is a characteristic radius for the
       !  largest collectors comprising the surface
       st_veg(jl)  =MAX(1.e-1_dp,(relax(jl)*(100.*pustveg(jl))**2.)/&
                        (g*1.e2*zAL))
       st_slsn(jl) =MAX(1.e-1_dp,(relax(jl)*(100.*pustslsn(jl))**2.)/visc )
       st_wat(jl)  =MAX(1.e-1_dp,(relax(jl)*(100.*pustarw(jl))**2.)/visc)

       !--- Calculation of the dry deposition velocity
       !    See paper slinn and slinn, 1980, vd is related to d**2/3
       !    over land, whereas over sea there is accounted for slip
       !    vb represents the contribution in vd of the brownian diffusion
       !    and vi represents the impaction.
       !
       !--- Over land, the vegetation and snow and bare soil fractions
       !    are considered:
       
       vb_veg(jl)  =(1./vkar)*((pustveg(jl)/abswind(jl))**2)&
                     *100.*abswind(jl)*(sc(jl)**(-2._dp/3._dp))
       vb_slsn(jl)  =(1./vkar)*((pustslsn(jl)/abswind(jl))**2) &
                     *100.*abswind(jl)*(sc(jl)**(-2._dp/3._dp))

       !  modified calculation of the impaction over
       !  vegetated surfaces, see paper by Gallagher et al., JGR 2002. 
       !  We have applied here the parameterization by Slinn [1982] 
       !  over vegetated surfaces
       vim_veg(jl)   =(1./vkar)*((pustveg(jl)/abswind(jl))**2)&
            *100.*abswind(jl)* (st_veg(jl)**2/(1._dp+st_veg(jl)**2))

       vim_slsn(jl)  =(1./vkar)*((pustslsn(jl)/abswind(jl))**2) &
                      *100.*abswind(jl)*(10.**(-3./st_slsn(jl)))

       !  the term evap has not been defined yet (30-01-2002) and the
       !  term daccm is set to zero anyhow. It still must be checked if
       !  this term should be included and how it relates to the factors
       !  that correct for particle growth close to the surface, according
       !  to Fitzjarald, 1975 (see ALPHAE and BETA)

       !!qqq delete ?
       vkd_veg(jl)  =pevap(jl)*daccm+(vb_veg(jl)+vim_veg(jl))
                      
       ! modified calculation of the surface resistance over 
       ! vegetated surfaces, see paper by Gallagher et al., JGR 2002. 
       ! The calculation includes the interception collection efficiency 
       !  vim and a rebound correction factor R

       vin_veg(jl)  =(1./vkar)*((pustveg(jl)/abswind(jl))**2)&
            *100.*abswind(jl)* (1._dp/2._dp)*(pradius(jl)/zAS)**2   
       zrebound(jl) =exp(-st_veg(jl)**0.5)       
       vkd_veg(jl) =pevap(jl)*daccm+zrebound(jl)&
            *(vb_veg(jl)+vim_veg(jl)+vin_veg(jl))
       vkd_slsn(jl) =pevap(jl)*daccm+(vb_slsn(jl)+vim_slsn(jl))

       !--- Over sea:
       !    Brownian diffusion for rough elements, see Hummelshoj
       !    re is the reynolds stress:
       re(jl)        =(100.*pustarw(jl)*100.*paz0m(jl))/visc
       vb_sea(jl)    =(1._dp/3._dp)*100.*pustarw(jl)*((sc(jl)**(-0.5))&
                       *re(jl)**(-0.5))
       vim_sea(jl)   =100.*pustarw(jl)*10.**(-3._dp/st_wat(jl))
       vkdaccsea(jl) =vb_sea(jl)+vim_sea(jl)
       vkd_wat(jl)   =(1._dp-palpha(jl))*vkdaccsea(jl)+palpha(jl)*pbubble(jl)

       ! Slinn and Slinn parameterization, calculation of vd:
       ! ====================================================
       ! without considering the role of sedimentation !
       ! This process is being calculated in a separate
       ! subroutine in echam5!
       ! ====================================================

       vkc_veg(jl)     =(100./prahveg(jl))
       vkc_slsn(jl)     =(100./prahslsn(jl))
       vkc_wat(jl)     =(100./prahw(jl))
       IF (ABS(vkd_veg(jl)) > TINY(vkd_veg(jl))) THEN
          vdpart_veg(jl)  =1./((1./vkc_veg(jl))+(1./vkd_veg(jl)))
       ELSE
          vdpart_veg(jl)  =vkc_veg(jl)
       ENDIF
       IF (ABS(vkd_slsn(jl)) > TINY(vkd_slsn(jl))) THEN
          vdpart_slsn(jl) =1./((1./vkc_slsn(jl))+(1./vkd_slsn(jl)))
       ELSE
          vdpart_slsn(jl) =vkc_slsn(jl)
       ENDIF
       IF (ABS(vkd_wat(jl)) > TINY(vkd_wat(jl))) THEN
          vdpart_wat(jl)  =1./((1./vkc_wat(jl))+(1./vkd_wat(jl)))
       ELSE
          vdpart_wat(jl)  = vkc_wat(jl)
       ENDIF
       !--- Calculate the dry deposition velocity weighted according 
       ! to the surface types:
       
       pvdrydep(jl) = pcvs(jl)    *vdpart_slsn(jl)  + &   ! snow fraction
                     pcvbs(jl)    *vdpart_slsn(jl)  + &   ! bare soil fraction
                     (1._dp-pcvs(jl))*(1._dp-pcvw(jl))*pvgrat(jl)* &
                           vdpart_veg(jl)           + &   ! vegetation fraction
                     (1._dp-pcvs(jl))*pcvw(jl)      * &
                           vdpart_veg(jl)           + &   ! wet skin fraction
                      pfri(jl)    *vdpart_slsn(jl)  + &   ! sea ice fraction
                      pcvwat(jl)  *vdpart_wat(jl)         ! water fraction

    ELSE

       pvdrydep(jl)=0.

    ENDIF
    ENDDO !jl

END SUBROUTINE drydep_vdaer

!===========================================================================

!=============================================================================

elemental real(dp) function drydep_calc_abswind(lowlimit, u  , v)

  real(dp), intent(in) :: lowlimit
  real(dp), intent(in) :: u, v

  INTRINSIC sqrt, MAX

  drydep_calc_abswind = sqrt(MAX(lowlimit, u**2+ v**2))

end function drydep_calc_abswind
!===========================================================================
elemental real(dp) function drydep_calc_radius_mass(pradius, sigma)

REAL(dp), intent(in) :: pradius
REAL(dp), INTENT(in) :: sigma

INTRINSIC EXP, LOG

! consistency with submodels SEDI and SCAV
drydep_calc_radius_mass = pradius *1.E2*EXP(3.*(LOG(sigma))**2)

END function drydep_calc_radius_mass

!===========================================================================
!===========================================================================
elemental real(dp) function drydep_calc_radius_num(pradius)

REAL(dp), intent(in) :: pradius

drydep_calc_radius_num = 1.E2 * pradius 

END function drydep_calc_radius_num

!===========================================================================

SUBROUTINE drydep_calc_landtypefractions(frac_water, frac_ice, &
                                          l_land, sealandmask, seaice )
  REAL(dp), INTENT(OUT) :: frac_water(:), frac_ice(:)
  LOGICAL,  INTENT(OUT) :: l_land(:)
  REAL(dp), INTENT(IN)  :: sealandmask(:), seaice(:)

  l_land(:) = .false.

  WHERE (sealandmask(:) >= 0.5_dp)
     l_land(:) = .true. 
  ENDWHERE

  frac_water(:) = (1._dp-sealandmask(:))*(1._dp-seaice(:))
  frac_ice(:)   = 1._dp - sealandmask(:)-frac_water(:)

END SUBROUTINE drydep_calc_landtypefractions

!===========================================================================

SUBROUTINE drydep_calc_vegfrac(pvgrat, pcvbs, &
                               pfvgrat, pcvs, pcvw, l_land)

REAL(dp), DIMENSION(:), INTENT(OUT) :: pvgrat, pcvbs
REAL(dp), DIMENSION(:), INTENT(IN)  :: pfvgrat, pcvs, pcvw
LOGICAL,  DIMENSION(:), INTENT(IN)  :: l_land

    WHERE (l_land(:))
       ! vegetation fraction
       pvgrat(:)=pfvgrat(:)                              
       ! bare soil fraction
       pcvbs(:)=(1._dp-pcvs(:))*(1._dp-pvgrat(:))*(1._dp-pcvw(:))  
    ELSEWHERE
       pvgrat(:) = 0.
       pcvbs(:)  = 0.
    END WHERE


END SUBROUTINE drydep_calc_vegfrac

!===========================================================================
 elemental real(dp) function drydep_calc_layerthickness(geopot)

USE messy_main_constants_mem, ONLY: g 

real(dp), INTENT(IN) :: geopot

drydep_calc_layerthickness = geopot/g

end function drydep_calc_layerthickness

!===========================================================================
elemental real(dp) function drydep_calc_vd_eff(vd,deltat,deltaz)

  real(dp), INTENT(in) :: vd, deltat, deltaz

  real(dp)  :: zqvd, vd_int
  
  INTRINSIC MAX, EXP

  zqvd=1./MAX(1.e-20_dp, vd)
  vd_int = -1./(zqvd*deltaz/deltat)

  drydep_calc_vd_eff = deltaz/deltat*(1._dp-EXP(vd_int)) 

end function drydep_calc_vd_eff

!===========================================================================
elemental real(dp) function drydep_posfinit(conctr,zddep,deltat,fac)

  real(dp), INTENT(in) :: conctr, zddep, deltat, fac
  real(dp) :: tmp1, tmp2

  tmp1 = conctr - zddep/fac*deltat
  IF (tmp1 < 0.0_dp) THEN
     tmp2 = conctr/deltat*fac
  ELSE
     tmp2 = zddep
  ENDIF

  drydep_posfinit = tmp2

end function drydep_posfinit
!===========================================================================

END MODULE messy_ddep

! *********************************************************************************
MODULE  messy_rad_fubrad_sr_str
  !
  ! PURPOSE:
  !    Provide a parameterization of the short wave heating rates due to 
  !    absorption of UV radiation in the Schumann-Runge continuum () and
  !    bands (175 - 205 nm) by oxygen.
  !
  ! REFERENCES: 
  !    Strobel, D. F., Parameterization of the atmospheric heating rate from 
  !    15 to 120 km due to O2 and O3 absorption of solar radiation,
  !    J. Geophys. Res., 83, 6225-6230, 1978.
  !
  ! CODE HISTORY:
  !  
  !
  !    Markus Kunze, FU-Berlin, 09/2015: 
  !
  !********************************************************************************
  
  USE messy_main_constants_mem, ONLY: dp, g, pi, N_A, M_O2, h_Planck, c_light
  USE messy_rad_fubrad_mem
  
  IMPLICIT NONE
  PRIVATE
  SAVE
  
  PUBLIC :: schumann_runge           &
          , src_hr_o2                &
          , srb_hr_o2                &
          , src_eff                  &
          , srb_hr_zhu               &
          , schuru_hr_o2             &
          , schuru_hr_o2_srb         &
          , schuru_flx_o2            &
          , src_flx_o2

  ! Interface of public subroutines
  
  INTERFACE schumann_runge
     MODULE PROCEDURE schumann_runge
  END INTERFACE
  
  INTERFACE srb_hr_zhu
     MODULE PROCEDURE srb_hr_zhu
  END INTERFACE
  
  INTERFACE schuru_hr_o2
     MODULE PROCEDURE schuru_hr_o2
  END INTERFACE
  
  INTERFACE schuru_hr_o2_srb
     MODULE PROCEDURE schuru_hr_o2_srb
  END INTERFACE
  
  INTERFACE src_hr_o2
     MODULE PROCEDURE src_hr_o2
  END INTERFACE
  
  INTERFACE srb_hr_o2
     MODULE PROCEDURE srb_hr_o2
  END INTERFACE
  
  INTERFACE schuru_flx_o2
     MODULE PROCEDURE schuru_flx_o2
  END INTERFACE
  
  INTERFACE src_flx_o2
     MODULE PROCEDURE src_flx_o2
  END INTERFACE
  
CONTAINS
  !
  !=============================================================================*
  !
  SUBROUTINE schumann_runge (KBDIM, KPROMA, KFLEV, KLEV, psrbht, psrcht, PO2, PPF)
    !                                                                         
    !*** *SCHUMANN_RUNGE* - MIDDLE ATMOSPHERIC SOLAR HEATING RATES DUE TO O2
    !
    !     S. PAWSON,  U. LANGEMATZ      BERLIN       8/92, 2/93
    !     K. NISSEN, BERLIN, 8/2005, CHANGES FOR MESSY
    !                        1/2006, Multiplication by Efficiency Factor
    !     K. MATTHES,BERLIN, 2/2008, changes in QSCH1, QSCH2, QSCHB
    !     M. KUNZE,  BERLIN, 2/2015, changes for variable oxygen concentration
    !
    !     PURPOSE
    !     -------
    !
    !     SOLAR HEATING RATES DUE TO O2 ABSORPTION OF SOLAR RADIATION IN THE
    !     MIDDLE ATMOSPHERE ARE CALCULATED USING THE MODEL OF *STROBEL (1978).
    !
    !**   METHOD
    !     ------
    !
    !     RADIATIVE HEATING RATES ARE CALCULATED USING A MODIFIED VERSION OF THE
    !     *STROBEL (1978) RADIATION SCHEME. ECMWF SOLAR GEOMETRICAL TERMS ARE
    !     USED. THE RESULTS ARE CALCULATED AT FULL MODEL LEVELS (HENCE PFLP).
    !
    !**   INTERFACE
    !     ---------
    !
    !     *SCHUMANN_RUNGE IS CALLED FROM *RADHEAT
    !
    !**   EXTERNALS
    !     ---------
    !
    !     NONE
    !
    !**   REFERENCES
    !     ----------
    !
    !     STROBEL (1978), J. GEOPHYS. RES., VOL 83, NO. C12, 6225-6230
    !
    
    ! 0.1 DEFINE VARIABLES
    ! --------------------
    
    !ARGUMENTS
    !
    ! INPUT:
    ! ------
    INTEGER, INTENT(in) :: KBDIM   !first dimension of 2-d arrays
    INTEGER, INTENT(in) :: KPROMA  !Number of longitudes 
    INTEGER, INTENT(in) :: KFLEV   !Number of full short wave levels
    INTEGER, INTENT(in) :: KLEV    !Number of full levels
    !
    REAL(dp), DIMENSION(KBDIM,KLEV), INTENT(in) :: ppf  ! Full Level Pressure
    REAL(dp), DIMENSION(kbdim,klev), INTENT(in) :: po2  ! O2 mass mixing ratio
    !
    ! OUTPUT:
    ! -------
    REAL(dp), DIMENSION(kbdim,kflev), INTENT(out):: psrbht ! hr schumann-runge bands
    REAL(dp), DIMENSION(kbdim,kflev), INTENT(out):: psrcht ! hr schumann-runge continuum
    
    !LOCAL
    INTEGER :: JK  ! Loop counter
    
    REAL(dp), DIMENSION(KBDIM) :: ZXL2  ! O2 column density along radiation path in cm^-2
    REAL(dp), DIMENSION(KBDIM) :: ZXN2  ! O2 number density along radiation path in cm^-3
    !                                 
    !  2.  INITIALISATION OF SOME DATA                                  
    !  -------------------------------
    
    !*    2.2  MAGNIFICATION FACTOR                                                
    !
    !ZSEC(1:KPROMA) = 35._dp/(DSQRT(1224._dp*(PMU0(1:KPROMA)*PMU0(1:KPROMA))+1._dp))
    !
    !  3.  CALCULATE HEATING RATES                                      
    !  ---------------------------                                     
    
    DO JK=1,KFLEV
       !
       ZXL2(1:KPROMA) = OXNDC_fct (po2(1:KPROMA,jk)) * &
                                   PRZSEC(1:KPROMA) * PPF(1:KPROMA,JK)
       ZXN2(1:KPROMA) = OXCNST_fct(po2(1:KPROMA,jk))
       !
       ! PSRCHT: Flux Divergence Schumann-Runge Continuum in  W/(m2 cm)
       ! PSRBHT: Flux Divergence Schumann-Runge Band      in  W/(m2 cm)
       ! SRC_EFF (REAL FUNCTION) EFFICIENCY FACTOR 
       !
       psrbht(1:kproma,jk) = QSCHB(ZXN2(1:kproma),ZXL2(1:kproma)) * oxfac
       psrcht(1:kproma,jk) = (QSCH1(zxn2(1:kproma),ZXL2(1:kproma))  + &
                              QSCH2(zxn2(1:kproma),ZXL2(1:kproma))) * &
                              SRC_EFF(ppf(1:kproma,jk)/100._dp) * oxfac
    END DO
    !
    RETURN
  END SUBROUTINE schumann_runge
  !
  SUBROUTINE src_hr_o2 (kbdim, kproma, kflev, klev, psrcht, po2, ppf)
    !
    ! INPUT:
    ! ------
    INTEGER, INTENT(in) :: kbdim   !first dimension of 2-d arrays
    INTEGER, INTENT(in) :: kproma  !Number of longitudes 
    INTEGER, INTENT(in) :: kflev   !Number of full short wave levels
    INTEGER, INTENT(in) :: klev    !Number of full levels
    !
    REAL(dp), DIMENSION(kbdim,klev), INTENT(in) :: ppf  ! Full Level Pressure
    REAL(dp), DIMENSION(kbdim,klev), INTENT(in) :: po2  ! O2 mass mixing ratio
    !
    ! OUTPUT:
    ! -------
    REAL(dp), DIMENSION(kbdim,kflev), INTENT(out):: psrcht ! hr shumann-runge continuum
    
    !LOCAL
    INTEGER :: jk  ! Loop counter
    
    REAL(dp), DIMENSION(kbdim) :: zxl2  ! O2 column density along radiation path in cm^-2
    REAL(dp), DIMENSION(kbdim) :: zxn2  ! O2 number density along radiation path in cm^-3
    !
    !  3.  CALCULATE HEATING RATES
    !  ---------------------------
    
    DO jk = 1, kflev
       !
       zxl2(1:kproma) = OXNDC_fct (po2(1:kproma,jk)) * &
                                   przsec(1:kproma) * ppf(1:kproma,jk)
       zxn2(1:kproma) = OXCNST_fct(po2(1:kproma,jk))
       !
       ! PSRCHT: Flux Divergence Schumann-Runge Continuum in  W/(m2 cm)
       ! SRC_EFF (REAL FUNCTION) EFFICIENCY FACTOR 
       !
       psrcht(1:kproma,jk) = (QSCH1(zxn2(1:kproma),zxl2(1:kproma))  + &
                              QSCH2(zxn2(1:kproma),zxl2(1:kproma))) * &
                              SRC_EFF(ppf(1:kproma,jk)/100._dp) * oxfac
    END DO
    !
    RETURN
  END SUBROUTINE src_hr_o2
  ! -
  SUBROUTINE srb_hr_o2(KBDIM, KPROMA, KFLEV, KLEV, psrbht, PO2, PPF)
    !
    ! INPUT:
    ! ------
    INTEGER, INTENT(in) :: KBDIM   !first dimension of 2-d arrays
    INTEGER, INTENT(in) :: KPROMA  !Number of longitudes 
    INTEGER, INTENT(in) :: KFLEV   !Number of full short wave levels
    INTEGER, INTENT(in) :: KLEV    !Number of full levels
    !
    REAL(dp), DIMENSION(KBDIM,KLEV), INTENT(in) :: ppf  ! Full Level Pressure
    REAL(dp), DIMENSION(kbdim,klev), INTENT(in) :: po2  ! O2 mass mixing ratio
    !
    ! OUTPUT:
    ! -------
    REAL(dp), DIMENSION(kbdim,kflev), INTENT(out):: psrbht ! hr schumann-runge bands
    
    !LOCAL
    INTEGER :: JK  ! Loop counter
    
    REAL(dp), DIMENSION(KBDIM) :: ZXL2  ! O2 column density along radiation path in cm^-2
    REAL(dp), DIMENSION(KBDIM) :: ZXN2  ! O2 number density along radiation path in cm^-3
    !
    !  3.  CALCULATE HEATING RATES                                      
    !  ---------------------------                                     
    
    DO JK=1,KFLEV
       !
       ZXL2(1:KPROMA) = OXNDC_fct (po2(1:KPROMA,jk)) * &
                                   PRZSEC(1:KPROMA) * PPF(1:KPROMA,JK)
       ZXN2(1:KPROMA) = OXCNST_fct(po2(1:KPROMA,jk))
       !
       ! PSRBHT: Flux Divergence Schumann-Runge Band      in  W/(m2 cm)
       ! SRC_EFF (REAL FUNCTION) EFFICIENCY FACTOR 
       !
       psrbht(1:kproma,jk) = QSCHB(ZXN2(1:kproma),ZXL2(1:kproma)) * oxfac

    END DO
    !
    RETURN
  END SUBROUTINE srb_hr_o2
  !
  ELEMENTAL FUNCTION QSCHB(OXCNST, ZL2) RESULT(qsrb)
    !
    ! Purpose: QSCHB: Schumann-Runge Bands (175-205nm)
    ! --------        a=0.67, b=3.44E9
    ! 
    ! Input:  zl2    - O2 column density along radiation path in cm^-2
    ! ------  OXCNST - O2 number density along radiation path in cm^-3
    !
    REAL(dp), INTENT(in) :: OXCNST
    REAL(dp), INTENT(in) :: ZL2
    REAL(dp)             :: qsrb
    !
    !TOTAL
    IF (ZL2 > 1.e18) THEN
       qsrb = cgs2si * soscale * OXCNST / (0.143_dp * ZL2 + 9.64E8_dp * DSQRT(ZL2))
    ELSE
       qsrb = cgs2si * soscale * OXCNST * 9.03e-19_dp 
    ENDIF

    RETURN
  END FUNCTION QSCHB
  !
  ELEMENTAL FUNCTION QSCH1(OXCNST, ZL2) RESULT(qsrc1)
    !
    ! Purpose: Schuman-Runge Continuum (125-175nm); QSCH1: 125-152nm
    ! --------
    !
    ! QSCH1 = value * oxcnst * exp(average cross section * O2 column)
    ! value = zFsch1 (flux in erg cm-2 s-1) * 1.E-3 (convert cgs to si units) 
    !         * 1.E-17 (cross section)
    !
    ! Input:
    ! ------
    REAL(dp), INTENT(in) :: OXCNST
    REAL(dp), INTENT(in) :: ZL2
    REAL(dp)             :: qsrc1
    !
    !
    ! TOTAL HEATING RATE
    
    qsrc1 = cgs2si * zFsch1 * 1.E-17_dp * OXCNST * EXP(-1.0E-17_dp * ZL2)
    
    RETURN
  END FUNCTION QSCH1
  !
  ELEMENTAL FUNCTION QSCH2(OXCNST,ZL2) RESULT(qsrc2)

    ! Purpose: Schuman-Runge Continuum (125-175nm); QSCH2: 152-166nm, 166-175nm
    ! QSCH2 = 1.E-3 (convert cgs to si units) * oxcnst * 
    !               (zFl *exp(cross section *O2 column)) +
    !               (zFd *exp(cross section * O2 column)) - 
    !               (zFs *exp(cross section * O2 column))
    !
    ! Input:
    ! ------
    REAL(dp), INTENT(in) :: OXCNST
    REAL(dp), INTENT(in) :: ZL2
    REAL(dp)             :: qsrc2
    !
    !TOTAL HEATING RATE
    
    qsrc2 = cgs2si * OXCNST *(zFl *EXP(-2.9E-19_dp * ZL2) &
         &    + zFd * EXP(-1.54E-18_dp * ZL2)     &
         &    - zFs * EXP(-1.1E-17_dp * ZL2)) / ZL2
    !
    RETURN
  END FUNCTION QSCH2
  !-
  ELEMENTAL FUNCTION src_eff(pres) RESULT(srceff)
    !
    ! PURPOSE:
    ! --------
    ! CALCULATE EFFICIENCY FACTOR FOR THE SCHUMANN-RUNGE CONTINUUM
    ! THIS FOLLOWS Mlynczak and Solomon, JGR, vol 98, p 10517, 1993   
    !
    ! AUTHOR: 
    ! -------
    ! K. Nissen, 16.1.2006, Free University of Berlin 
    !
    REAL(dp), INTENT(in) :: PRES    !pressure in hPa
    REAL(dp)             :: srceff
    !
    REAL(dp), PARAMETER, DIMENSION(4) :: c =  &
                      (/0.75349_dp, 0.0036_dp, 0.059468_dp, -0.022795_dp/)
    REAL(dp) :: X
    
    IF (pres > 1.e-2_dp) THEN
       srceff = 0.7938_dp
    ELSE IF (pres > 1.e-4_dp) THEN
       x = LOG10(pres)+3._dp
       srceff = c(1) + x*c(2) + x*x*c(3) + x*x*x*c(4)
    ELSE
       srceff = 0.8320_dp
    END IF
    !
    RETURN
  END FUNCTION src_eff
  !
  !==============================================================================
  !==============================================================================
  !
  ! ------------------------------------------------------------------
  SUBROUTINE srb_hr_zhu(nbdim, nproma, nswlev, nlev, po2, ppf, psrbht, psrcht)
    ! INPUT:
    ! ------
    INTEGER,                         INTENT(in) :: nbdim   ! first dimension of 2-d arrays
    INTEGER,                         INTENT(in) :: nproma  ! Number of longitudes 
    INTEGER,                         INTENT(in) :: nswlev  ! Number of full FUBRAD levels
    INTEGER,                         INTENT(in) :: nlev    ! Number of full levels
    REAL(dp), DIMENSION(nbdim,nlev), INTENT(in) :: po2     ! O2 mass mixing ratio
    REAL(dp), DIMENSION(nbdim,nlev), INTENT(in) :: ppf     ! pressure on full levels
    ! OUTPUT:
    ! -------
    REAL(dp), DIMENSION(nbdim,nswlev), INTENT(out):: psrbht ! hr shumann-runge bands
    REAL(dp), DIMENSION(nbdim,nswlev), INTENT(out):: psrcht ! hr shumann-runge continuum
    
    !LOCAL
    INTEGER :: jk  ! Loop counter
    
    !
    ! srb_eff  - to achieve a realistic heating rate with the provided integrated
    !            flux in the Schumann-Runge bands, it is necessary to scale the 
    !            parametrized heating rate with a scaling factor 0.333.
    !            But according to Mlynczak and Solomon (1993) the heating efficiency
    !            in the Schumann-Runge bands is 1., as the photons are not sufficiently
    !            energetic to create exited O_1D photolysis products.
    REAL(dp), PARAMETER        :: srb_eff = 1._dp/3._dp
    REAL(dp), DIMENSION(nbdim) :: zxl2_top ! O2 column density (uppermost level)
    REAL(dp), DIMENSION(nbdim) :: zxl2  ! O2 column density along radiation path in cm^-2
    REAL(dp), DIMENSION(nbdim) :: zxn2  ! O2 number density along radiation path in cm^-3
    
    zxl2_top(1:nproma) = OXNDC_fct ( po2(1:nproma,1) ) * &
                                  przsec(1:nproma) * ppf(1:nproma,1)
    DO jk = 1, nswlev
       !
       zxl2(1:nproma) = OXNDC_fct ( po2(1:nproma,jk) ) * &
                                 przsec(1:nproma) * ppf(1:nproma,jk)
       zxn2(1:nproma) = OXCNST_fct( po2(1:nproma,jk) )
       !
       ! PSRBHT: Flux Divergence Schumann-Runge Band      in  W/(m2 cm)
       ! SRC_EFF (REAL FUNCTION) EFFICIENCY FACTOR 
       !
       psrbht(1:nproma,jk) = QSCHB_ZHU(zxn2(1:nproma),zxl2(1:nproma), &
                                   zxl2_top(1:nproma)) * &
                                   oxfac
                             !srb_eff * oxfac
       psrcht(1:nproma,jk) = (QSCH1(zxn2(1:nproma),zxl2(1:nproma))  + &
                              QSCH2(zxn2(1:nproma),zxl2(1:nproma))) * &
                             SRC_EFF(ppf(1:nproma,jk)/100._dp) * oxfac
    END DO
    RETURN 
  END SUBROUTINE srb_hr_zhu
  !
  ELEMENTAL FUNCTION QSCHB_ZHU(pxn2, pxl2, pxl2_top) RESULT(qsrb)
    !
    ! Purpose: Calculate heating rates in Schumann-Runge Bands (175-205nm)
    ! -------- according to Zhu (1994).
    ! 
    ! Input:  pxl2 - O2 column density along radiation path in cm^-2
    ! ------  pxn2 - O2 number density along radiation path in cm^-3
    !
    REAL(dp), INTENT(in) :: pxn2
    REAL(dp), INTENT(in) :: pxl2
    REAL(dp), INTENT(in) :: pxl2_top
    REAL(dp)             :: qsrb
    REAL(dp), PARAMETER  :: zsigsrb = 2.07E-20_dp ! O2 cross section in cm^2
    REAL(dp), PARAMETER  :: ysrb = 0.0152_dp
    REAL(dp)             :: xsrb                  ! efficiency factor
    REAL(dp)             :: fact1, fact2
    !
    xsrb = (pxl2_top/pxl2)**0.3_dp
    !
    fact1 = DSQRT(1._dp + 4._dp * zsigsrb * pxl2 / (pi * ysrb))
    fact2 = pi * ysrb * 0.5_dp
    !
    qsrb = pxn2 * zFsrb(1) * xsrb * zsigsrb / fact1 * EXP( -fact2 * (fact1 - 1._dp) )
    !
    RETURN
  END FUNCTION QSCHB_ZHU
  !                                                                         
  !*** *schuru_hr_o2* - MIDDLE ATMOSPHERIC SOLAR HEATING RATES DUE TO O2
  !
  !     S. PAWSON,  U. LANGEMATZ      BERLIN       8/92, 2/93
  !     K. NISSEN,                    BERLIN, 8/2005, CHANGES FOR MESSY 
  !                                           1/2006, Multiplication by Efficiency Factor
  !
  !     PURPOSE
  !     -------
  !
  !     SOLAR HEATING RATES DUE TO O2 ABSORPTION OF SOLAR RADIATION IN THE 
  !     MIDDLE ATMOSPHERE ARE CALCULATED USING THE MODEL OF *STROBEL (1978).
  ! 
  !**   METHOD
  !     ------
  !     
  !     RADIATIVE HEATING RATES ARE CALCULATED USING A MODIFIED VERSION OF THE 
  !     *STROBEL (1978) RADIATION SCHEME. ECMWF SOLAR GEOMETRICAL TERMS ARE
  !     USED. THE RESULTS ARE CALCULATED AT FULL MODEL LEVELS (HENCE PFLP).
  !       
  !**   REFERENCES                                   
  !     ----------                                             
  !
  !     STROBEL (1978), J. GEOPHYS. RES., VOL 83, NO. C12, 6225-6230
  !                                         
  ! --------------------
  SUBROUTINE schuru_hr_o2(kbdim, kproma, kswlev, klev, psrb_hr, psrc_hr, po2, ppf)
    !
    ! 0.1 DEFINE VARIABLES
    ! --------------------

    !ARGUMENTS
    INTEGER,  INTENT(in) :: kbdim           ! first dimension of 2-d arrays
    INTEGER,  INTENT(in) :: kproma          ! Number of longitudes 
    INTEGER,  INTENT(in) :: kswlev          ! Number of full short wave levels
    INTEGER,  INTENT(in) :: klev            ! Number of full levels
    !
    REAL(dp), DIMENSION(kbdim,klev), INTENT(in) :: ppf  ! Full Level Pressure
    REAL(dp), DIMENSION(kbdim,klev), INTENT(in) :: po2  ! O2 mass mixing ratio
    !
    REAL(dp), DIMENSION(kbdim,kswlev), INTENT(out):: psrb_hr ! hr shumann-runge bands
    REAL(dp), DIMENSION(kbdim,kswlev), INTENT(out):: psrc_hr ! hr shumann-runge continuum

    !LOCAL
    INTEGER                    :: jk, i  ! Loop counter
    REAL(dp), DIMENSION(kbdim) :: zxn2   ! O2 number density along radiation path in cm^-3
    REAL(dp), DIMENSION(kbdim) :: zxl2   ! O2 Column
    REAL(dp), DIMENSION(kbdim) :: qsrb   ! heating rates Schumann-Runge bands
    REAL(dp), DIMENSION(kbdim) :: qsrc   ! heating rates Schumann-Runge continuum
    !
    !  3.  CALCULATE HEATING RATES                                      
    !  ---------------------------                                     
    DO jk = 1, kswlev
       zxl2(1:kproma) = OXNDC_fct (po2(1:kproma,jk)) * przsec(1:kproma) * ppf(1:kproma,jk)
       zxn2(1:kproma) = OXCNST_fct(po2(1:kproma,jk))
       qsrb(:) = 0._dp
       qsrc(:) = 0._dp
       !
       !CDIR unroll=nsrb_st
       DO i = 1, nsrb_st
          qsrb(1:kproma) = qsrb(1:kproma) + &
                           SRBsig(zxl2(1:kproma),i) * zxn2(1:kproma) * &
                           FSRB(zFsrb(i), zxl2(1:kproma),i)
       END DO
       !CDIR unroll=nsrc
       DO i = 1, nsrc
          qsrc(1:kproma) = qsrc(1:kproma) +  &
                           zsigsrc(i) * zxn2(1:kproma) * &
                           zFsrc(i) * EXP(-zsigsrc(i) * zxl2(1:kproma))
       END DO
       !
       !   psrc_hr: Flux Divergence Schumann-Runge Continuum in  W/(m2 cm)
       !   psrb_hr: Flux Divergence Schumann-Runge Band      in  W/(m2 cm)
       !   SRC_EFF (REAL FUNCTION) EFFICIENCY FACTOR 
       !
       psrc_hr(1:kproma,jk) = qsrc(1:kproma)* oxfac * SRC_EFF(ppf(1:kproma,jk)/100._dp)
       psrb_hr(1:kproma,jk) = qsrb(1:kproma)* oxfac
    END DO
    RETURN                                                                    
  END SUBROUTINE schuru_hr_o2
  SUBROUTINE schuru_hr_o2_srb(kbdim, kproma, kswlev, klev, psrb_hr, po2, ppf)
    !
    ! 0.1 DEFINE VARIABLES
    ! --------------------

    !ARGUMENTS
    INTEGER,  INTENT(in) :: kbdim           ! first dimension of 2-d arrays
    INTEGER,  INTENT(in) :: kproma          ! Number of longitudes 
    INTEGER,  INTENT(in) :: kswlev          ! Number of full short wave levels
    INTEGER,  INTENT(in) :: klev            ! Number of full levels
    !
    REAL(dp), DIMENSION(kbdim,klev), INTENT(in) :: ppf  ! Full Level Pressure
    REAL(dp), DIMENSION(kbdim,klev), INTENT(in) :: po2  ! O2 mass mixing ratio
    !
    REAL(dp), DIMENSION(kbdim,kswlev), INTENT(out):: psrb_hr ! hr shumann-runge bands

    !LOCAL
    INTEGER                    :: jk, i  ! Loop counter
    REAL(dp), DIMENSION(kbdim) :: zxn2   ! O2 number density along radiation path in cm^-3
    REAL(dp), DIMENSION(kbdim) :: zxl2   ! O2 Column
    REAL(dp), DIMENSION(kbdim) :: qsrb   ! heating rates Schumann-Runge bands
    !
    !  3.  CALCULATE HEATING RATES                                      
    !  ---------------------------                                     
    DO jk = 1, kswlev
       zxl2(1:kproma) = OXNDC_fct (po2(1:kproma,jk)) * przsec(1:kproma) * ppf(1:kproma,jk)
       zxn2(1:kproma) = OXCNST_fct(po2(1:kproma,jk))
       qsrb(:) = 0._dp
       !
       !CDIR unroll=nsrb_st
       DO i = 1, nsrb_st
          qsrb(1:kproma) = qsrb(1:kproma) + &
                           SRBsig(zxl2(1:kproma),i) * zxn2(1:kproma) * &
                           FSRB(zFsrb(i), zxl2(1:kproma),i)
       END DO
       !
       !   psrb_hr: Flux Divergence Schumann-Runge Band      in  W/(m2 cm)
       !
       psrb_hr(1:kproma,jk) = qsrb(1:kproma)* oxfac
    END DO
    RETURN                                                                    
  END SUBROUTINE schuru_hr_o2_srb
  !
  ELEMENTAL FUNCTION OXCNST_fct(o2_mmr) RESULT(x)
    !  
    ! Purpose:
    ! --------
    ! constant to derive number density in cm-3 from MMR:
    !
    ! oxcnst = O2_MMR * Air density * (N_Avogadro/M_O2) * 1E-6
    !
    ! 1e-6 converts from cm-3 to m-3.
    ! The density factor is ignored as it cancels when the 
    ! heating rate in K s-1 is computed from that in W m-3
    !
    !REAL(DP), PARAMETER :: oxcnst=4.3549E18_dp ! original value
    !
    ! AUTHOR: 
    ! -------
    ! M. Kunze, Freie Universitaet Berlin, 2010. 
    !
    REAL(dp), INTENT(in) :: o2_mmr ! O2 mass mixing ratio
    ! OUTPUT
    REAL(dp) :: x                  ! number density in molecules cm^-3
    ! 
    ! As M_O2 is in g/mol use a factor of 1.E-3 to convert to kg/mol.
    !
    !
    ! The factor 1.E-6_dp is used to convert from m^-3 to cm^-3.
    !
    x = o2_mmr * N_A/(M_O2*1.E-3_dp) * 1.E-6_dp
    !
    RETURN
  END FUNCTION OXCNST_fct
  ! - -----------------------------------------------------------------------
  ELEMENTAL FUNCTION OXNDC_fct(o2_mmr) RESULT(x)
    !
    ! Purpose: Derive oxygen number density in cm-3 from MMR:
    ! --------
    ! Calculate oxndc (almost) oxygen number density
    ! used as a factor to determine o2 column.
    ! 
    ! oxndc = O2_MMR * (N_Avogadro/M_O2) * (1E-4/g)

    !
    ! Equation to determine O2_column:
    ! 
    !     O2_VMR*(M_O2/M_air)*Air density*(N_Avogadro/M_O2)*1e-6
    ! =>  O2_MMR             *Air density*(N_Avogadro/M_O2)*1e-6
    ! 
    ! substitute air density by using hydrostatic balance and
    ! integrate from the model top to layer
    !
    !REAL(DP), PARAMETER :: oxndc=4.442E19_dp ! original value
    !
    ! To get the original value, o2_vmr = 0.20955775081737177.
    ! With M_air in g/mol.
    ! o2_mmr * N_A/(M_O2*1.E-3) * 1.E-4/g = 4.442E19
    !
    ! AUTHOR: 
    ! -------
    ! M. Kunze, June 2010, Free University of Berlin 
    !
    REAL(dp), INTENT(in) :: o2_mmr
    REAL(dp) :: x
    ! 
    ! As M_O2 is in g/mol use a factor of 1.E-3 to convert to kg/mol.
    !
    x = o2_mmr * N_A/(M_O2*1.E-3_dp) * 1.E-4_dp/g
    ! 
    RETURN
  END FUNCTION OXNDC_fct
  ! - -----------------------------------------------------------------------
  ! - -----------------------------------------------------------------------
  SUBROUTINE schuru_flx_o2(kbdim, kproma, kswlevp1, klevp1, psrb_fl, psrc_fl, po2c)
    !                                                                         
    !*** *schuru_flx_o2* - MIDDLE ATMOSPHERIC FLUXES IN THE SCHUMANN-RUNGE
    !                     CONTINUUM AND BANDS
    !     PURPOSE
    !     -------
    !     MIDDLE ATMOSPHERIC FLUXES IN THE SCHUMANN-RUNGE CONTINUUM AND BANDS ARE 
    !     CALCULATED USING THE DETAILED INTERVALS OF STROBEL (1978).
    ! 
    !**   METHOD
    !     ------
    !     MIDDLE ATMOSPHERIC FLUXES IN THE SCHUMANN-RUNGE CONTINUUM AND BANDS ARE
    !     CALCULATED USING THE DETAILED INTERVALS OF THE STROBEL (1978) RADIATION SCHEME. 
    !     ECMWF SOLAR GEOMETRICAL TERMS ARE USED. THE FLUXES ARE CALCULATED AT 
    !     HALF MODEL LEVELS (HENCE PFLP).
    !     SOLAR HEATING RATES DUE TO O2 ABSORPTION OF SOLAR RADIATION IN THE 
    !     MIDDLE ATMOSPHERE ARE CALCULATED USING THE MODEL OF *STROBEL (1978).
    !       
    !**   REFERENCES                                   
    !     ----------                                             
    !     STROBEL (1978), J. GEOPHYS. RES., VOL 83, NO. C12, 6225-6230
    ! --------------------
    !ARGUMENTS
    INTEGER,  INTENT(in) :: kbdim        ! first dimension of 2-d arrays
    INTEGER,  INTENT(in) :: kproma       ! Number of longitudes 
    INTEGER,  INTENT(in) :: kswlevp1     ! Number of half short wave levels
    INTEGER,  INTENT(in) :: klevp1       ! Number of half levels
    !
    REAL(dp), DIMENSION(kbdim,klevp1), INTENT(in)   :: po2c  ! O2 column in path length (atm-cm)
    !
    REAL(dp), DIMENSION(kbdim,kswlevp1), INTENT(out):: psrb_fl ! flux shumann-runge bands
    REAL(dp), DIMENSION(kbdim,kswlevp1), INTENT(out):: psrc_fl ! flux shumann-runge continuum

    !LOCAL
    INTEGER                    :: jk, i  ! Loop counter
    REAL(dp), DIMENSION(kbdim) :: zxl2   ! effective O2 Column
    !
    !  1.  CALCULATE FLUXES
    !  --------------------  
    psrc_fl(:,:) = 0._dp
    psrb_fl(:,:) = 0._dp
    DO jk = 1, kswlevp1
       !
       zxl2(1:kproma) = po2c(1:kproma,jk) * przsec(1:kproma)
       !
       !CDIR unroll=nsrb_st
       DO i = 1, nsrb_st
          psrb_fl(1:kproma,jk) = psrb_fl(1:kproma,jk) + &
                                 FSRB(zFsrb(i), zxl2(1:kproma),i)
       END DO
       !CDIR unroll=nsrc
       DO i = 1, nsrc
          psrc_fl(1:kproma,jk) = psrc_fl(1:kproma,jk) +  &
                                 zFsrc(i) * EXP(-zsigsrc(i) * zxl2(1:kproma))
       END DO
    END DO
    RETURN                                                                    
  END SUBROUTINE schuru_flx_o2
  !
  !   ------------------------------------------------------------------
  ELEMENTAL FUNCTION FSRB(pFsrb, zl2, ib) RESULT(x)
    !
    ! FSRB: Calculate the flux in the Schumann-Runge bands (175-205nm)
    ! -----
    !
    ! Reference: Strobel, 1978: Parametrization of atmospheric heating rates
    ! ----------                from 15 to 120 km due to O2 and O3 absorption
    !                           of solar radiation, JGR, Vol. 83, No. C12.
    !
    ! Author: Markus Kunze, FU-Berlin, 2014
    ! ------- 
    !
    REAL(dp), INTENT(in) :: pFsrb ! solar flux in band ib of Schumann-Runge bands in W m-2
    REAL(dp), INTENT(in) :: zl2   ! column density of oxygen
    INTEGER,  INTENT(in) :: ib    ! index of actual band
    ! OUTPUT:
    REAL(dp) :: x
    !
    REAL(dp), DIMENSION(nsrb_st), PARAMETER ::  &
       gamma = (/ 5.2699E-20_dp, 2.1935E-20_dp, 8.7466E-21_dp, 3.0367E-21_dp, &
                  2.1219E-21_dp, 1.1067E-21_dp, 4.6920E-22_dp, 2.1977E-22_dp, &
                  1.0391E-22_dp, 7.7757E-23_dp, 4.9293E-23_dp, 1.5313E-23_dp, &
                  1.8859E-23_dp, 1.4519E-23_dp, 1.2204E-23_dp, 1.2899E-23_dp, &
                  1.3033E-23_dp, 1.2662E-23_dp, 1.1983E-23_dp &
                 /)
    REAL(dp), DIMENSION(nsrb_st), PARAMETER ::  &
       delta = (/ 3.1771E-10_dp, 2.4514E-10_dp, 1.8609E-10_dp, 1.0362E-10_dp, &
                  6.1440E-11_dp, 6.7403E-11_dp, 9.0615E-11_dp, 5.0541E-11_dp, &
                  3.5539E-11_dp, 3.1147E-11_dp, 1.9898E-11_dp, 1.2232E-11_dp, &
                  8.3986E-12_dp, 5.4926E-12_dp, 2.6838E-12_dp, 1.7876E-12_dp, &
                  3.0529E-13_dp, 9.6723E-14_dp, 6.4880E-14_dp &
                /)
    !
    x = pFsrb * EXP( -(gamma(ib)*zl2 + delta(ib)*DSQRT(zl2)) )
    !
    RETURN
  END FUNCTION FSRB
  ! -
  ELEMENTAL FUNCTION SRBsig(zl2, iband) RESULT (sig)
    ! 
    ! SRBsig: Calculate the absobtion cross-sections for the 
    ! ------- Schumann-Runge bands (175-205nm).
    !
    ! Reference: Strobel, 1978: Parametrization of atmospheric heating rates
    ! ----------                from 15 to 120 km due to O2 and O3 absorption
    !                           of solar radiation, JGR, Vol. 83, No. C12.
    !
    ! Author: Markus Kunze, FU-Berlin, 2014
    ! ------- 
    !
    REAL(dp), INTENT(in) :: zl2   ! column density of oxygen
    INTEGER,  INTENT(in) :: iband ! index of actual band
    ! -
    REAL(dp) :: sig ! absorption cross section
    ! -
    REAL(dp), DIMENSION(nsrb_st), PARAMETER ::  &
         alpha = (/  3.3967E18_dp, 5.0849E18_dp, 5.5205E18_dp, 8.2969E18_dp, &
                     4.0394E18_dp, 1.2199E18_dp, 9.0621E18_dp, 7.2957E17_dp, &
                     1.9521E19_dp, 5.9996E19_dp, 6.8431E19_dp, 1.4269E20_dp, &
                     4.6314E20_dp, 1.2321E21_dp, 1.4252E22_dp, 1.8757E22_dp, &
                     7.3112E22_dp, 8.0106E22_dp, 8.2595E22_dp &
                  /)
    REAL(dp), DIMENSION(nsrb_st), PARAMETER ::  &
         beta  = (/  1.4959E09_dp, 2.5956E09_dp, 5.2439E09_dp, 1.1984E10_dp, &
                     2.4928E10_dp, 2.4265E10_dp, 1.5782E10_dp, 3.1870E10_dp, &
                     3.9664E10_dp, 4.6831E10_dp, 7.8005E10_dp, 1.3011E11_dp, &
                     1.7798E11_dp, 2.5293E11_dp, 1.0586E11_dp, 1.0081E11_dp, &
                     2.2040E09_dp,-3.5546E09_dp, 6.9850E08_dp  &
                    /)           
    ! -
    sig = 1._dp/(alpha(iband) + beta(iband) * DSQRT(zl2))
    ! -
    RETURN
  END FUNCTION SRBsig
  !
  SUBROUTINE src_flx_o2(kbdim, kproma, kswlevp1, klevp1, psrc_fl, po2c)
    !                                                                         
    !*** *src_flx_o2* - MIDDLE ATMOSPHERIC FLUXES IN THE SCHUMANN-RUNGE
    !                   CONTINUUM
    !     PURPOSE
    !     -------
    !     MIDDLE ATMOSPHERIC FLUXES IN THE SCHUMANN-RUNGE CONTINUUM ARE CALCULATED 
    !     USING THE DETAILED INTERVALS OF STROBEL (1978).
    ! 
    !**   METHOD
    !     ------
    !     MIDDLE ATMOSPHERIC FLUXES IN THE SCHUMANN-RUNGE CONTINUUM CALCULATED 
    !     USING A MODIFIED VERSION OF THE STROBEL (1978) RADIATION SCHEME. 
    !     ECMWF SOLAR GEOMETRICAL TERMS ARE USED. THE FLUXES ARE CALCULATED AT 
    !     HALF MODEL LEVELS (HENCE PFLP).
    !       
    !**   REFERENCES                                   
    !     ----------                                             
    !     STROBEL (1978), J. GEOPHYS. RES., VOL 83, NO. C12, 6225-6230
    ! --------------------
    !ARGUMENTS
    INTEGER,  INTENT(in) :: kbdim        ! first dimension of 2-d arrays
    INTEGER,  INTENT(in) :: kproma       ! Number of longitudes 
    INTEGER,  INTENT(in) :: kswlevp1     ! Number of half short wave levels
    INTEGER,  INTENT(in) :: klevp1       ! Number of half levels
    !
    REAL(dp), DIMENSION(kbdim,klevp1), INTENT(in)   :: po2c  ! O2 column in path length (atm-cm)
    !
    REAL(dp), DIMENSION(kbdim,kswlevp1), INTENT(out):: psrc_fl ! flux shumann-runge continuum

    !LOCAL
    INTEGER                    :: jk, i  ! Loop counter
    REAL(dp), DIMENSION(kbdim) :: zxl2   ! effective O2 Column
    !
    !  1.  CALCULATE FLUXES                                  
    !  ---------------------  
    psrc_fl(:,:) = 0._dp
    DO jk = 1, kswlevp1
       !
       zxl2(1:kproma) = po2c(1:kproma,jk) * przsec(1:kproma)
       !
       !CDIR unroll=nsrc
       DO i = 1, nsrc
          psrc_fl(1:kproma,jk) = psrc_fl(1:kproma,jk) +  &
                                 zFsrc(i) * EXP(-zsigsrc(i) * zxl2(1:kproma))
       END DO
    END DO
    RETURN                                                                    
  END SUBROUTINE src_flx_o2
  !
  FUNCTION eff_srb(iband) RESULT(eff)
    INTEGER,  INTENT(in) :: iband ! index of actual band
    
    REAL(dp) :: x1, x2
    REAL(dp), DIMENSION(nsrb_st), PARAMETER ::  &
         lambda = (/ 175.00_dp, 176.32_dp, 176.86_dp, 177.46_dp, 178.26_dp, &
                     179.26_dp, 180.36_dp, 181.64_dp, 183.06_dp, 184.62_dp, &
                     186.34_dp, 188.22_dp, 190.24_dp, 192.40_dp, 194.70_dp, &
                     197.18_dp, 198.50_dp, 200.00_dp, 202.50_dp  &
                   /)
    REAL(dp) :: eff
    LOGICAL, SAVE :: lfirst = .TRUE.
    
    x1 = h_Planck * c_light/(lambda(iband)*1.E-9_dp)
    x2 = x1 - 494.E3_dp
    IF (lfirst) PRINT *,':X1:',x1,':x2:',x2
    eff = x2 / x1
    lfirst = .FALSE.
    RETURN
  END FUNCTION
  
END MODULE  messy_rad_fubrad_sr_str

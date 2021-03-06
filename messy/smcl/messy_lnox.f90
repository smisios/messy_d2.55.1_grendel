! ***********************************************************************
MODULE messy_lnox
! ***********************************************************************

  ! Lightning parameterizations 
  ! (Price and Rind / Grewe / Allen and Pickering (2) / Finney et al. / Dahl)
  !
  ! Authors: V. Grewe, DLR Oberpfaffenhofen, (volker.grewe@dlr.de)
  !             - Written for ECHAM4.L39(DLR)/CHEM
  !          P. Joeckel, MPI-CH Mainz, (joeckel@mpch-mainz.mpg.de), Aug 2003
  !             - f77 -> f90, IMPICIT NONE
  !             - ECHAM4 independent MODULE
  !             - messy_lnox_e5 (later _si): ECHAM5/MESSy INTERFACE for LNOx
  !             - partly separation of resolution dependent parts
  !               => MESSy-SMCL FOR SUBMODEL LNOX
  !               => LNOX = LIGHTNING NOx EMISSION PARAMETERIZATION (LNOx)
  !          H. Tost, MPICH, January 2007
  !             - additional parameterisations implemented 
  !               (following Allen and Pickering)
  !          A. Kerkweg, UNI-MZ, June 2012
  !             - additional parameterisation implemented
  !               (Dahl, 2010)
  !          P. Joeckel, DLR, November 2012
  !             - revision, l_mode_scal introduced to distinguish
  !               production mode from "readjust scaling mode"
  !             - claned up
  !          P. Joeckel, DLR, September 2014
  !             - implementation of ice flux parameterisation
  !             - submodel completely revised
  !          Francisco JavierPerez-Invernon, DLR, Aug 2020:
  !             - implementation of extended Finney parameterisation
  !               (ice flux through external iso-surface)
  !            
  ! Literature:
  ! V. Grewe, D. Brunner, M. Dameris, J.L. Grenfell, R. Hein, D. Shindell,
  !    and J. Staehelin, Origin and variability of upper tropospheric
  !    nitrogen oxides and ozone at northern mid-latitudes, Atm. Env., 35,
  !    3421-3433, 2001.
  ! C. Kurz and V. Grewe, Meteorol. Z., 11, 379-393, 2002.
  ! C. Price, J. Penner, and M. Prather, NOx from lightning 1. Global
  !    distribution based on lightning physics, JGR, 102(D5), 5929-5941,
  !    1997.
  ! K. E. Pickering, Y. Wang, W.-K. Tao, C. Price, and J. F. Muller, 
  !    Vertical distributions of lightning NOx for use in regional and global
  !    chemical transport models, JGR, 103(D23), 31203-31216, 1998.
  ! C. Price, and D. Rind, Modeling Global Distributions in a General
  !    Circulation Model, Month. Wea. Rev., 122, 1930-1939, 1994.
  ! C. Price, and D. Rind, A Simple Parameterization for Calculating
  !    Global Lightning Distributions, JGR, 97(D9), 9919-9933, 1992.
  ! D. J. Allen and K. E. Pickering, Evaluation of lightning flash 
  !    parameterizations for use in a global chemical transport model, 
  !    JGR, 107(D23), doi:10.129/2002JD002066, 2002
  ! Finney, D. L., Doherty, R. M., Wild, O., Huntrieser, H., Pumphrey, H. C.,
  !    and Blyth, A. M.: Using cloud ice flux to parametrise large-scale
  !    lightning, Atmos. Chem. Phys., 14, 12665-12682, 
  !    doi:10.5194/acp-14-12665-2014, 2014.  
  !
  ! TODO:
  !    see TODO list in messy_lnox_si.f90

  USE messy_main_constants_mem, ONLY: DP

  IMPLICIT NONE
  PUBLIC
  SAVE

  CHARACTER(LEN=*), PARAMETER :: modstr = 'lnox'
  CHARACTER(LEN=*), PARAMETER :: modver = '3.0'

  ! PARAMETERISATIONS FOR FLASH FREQUENCY CALCULATION
  INTEGER, PARAMETER :: IPARAM_PaR_T = 1  ! 1: Price & Rind
  INTEGER, PARAMETER :: IPARAM_Grewe = 2  ! 2: Grewe
  INTEGER, PARAMETER :: IPARAM_AaP_M = 3  ! 3: Allen & Pickering: massflux
  INTEGER, PARAMETER :: IPARAM_AaP_P = 4  ! 4: Allen & Pickering: precipitation
  INTEGER, PARAMETER :: IPARAM_FinIF = 5  ! 5: Finney et al.: cloud ice flux
  INTEGER, PARAMETER :: IPARAM_extIF = 6  ! 6: Finney et al. extended: ice flux
  INTEGER, PARAMETER :: IPARAM_DahlC = 7  ! 7: Dahl (mmust be last in list!)
  INTEGER, PARAMETER :: NPARAM = 7        ! number of parameterisations
  INTEGER, PARAMETER :: IPARAM_ALL   = 9  ! run all diagnostically

  !
  CHARACTER(LEN=6), DIMENSION(NPARAM) :: paramnames = (/ &
       '_PaR_T', '_Grewe', '_AaP_M', '_AaP_P', '_FinIF', '_extIF', '_DahlC'/)

  ! see funcion lflash below
  REAL(DP), PARAMETER :: fpm2s_min = 1.433E-14_dp ! [1/(m^2 s)]

  ! CTRL NAMELIST PARAMETERS
  LOGICAL         :: l_mode_scal = .FALSE. ! TRUE: diagnostic simulation to
  !                                               re-adjust scaling
  !                                        ! FALSE: production simulation
  ! NOTES: 
  !     - To readjust the scaling parameter(s) r_scal_ff in the namelist,
  !       perform a simulation (no need to switch on "chemistry"!!!) with
  !       l_mode_scal = .TRUE. and (preferably) r_scal_ff(.) = 1.0.
  !       From the resulting output (make sure that output is T in 
  !       channel.nml!!!) calculate the required scaling factor r_scal_ff:
  !     - The global flash frequency should be on the order of 50-100 flashes/s
  !       (LIS OTD global annual average (1995-2005): 46.6 flashes/s.
  !     - The cloud-to-ground flash rate should be about 1/10 the in-cloud 
  !       flash rate, i.e., ffcg/ffic = pg / (1-pg) ~ 0.1.
  !     - Flash frequencies below a certain cut-off will be ignored,
  !       except for in the scaling mode; see function lflash.
  !     - The average ligthning NOx production in Tg(N)/year is:
  !       (sum over boxes) flag * r_scal_ff * (ffcg + r_eff*ffic) * r_noxpf ...
  !       ... * 3600 s/h * 24 h/d * 365 d/y / 1e9 kg/Tg
  ! 
  INTEGER         :: i_ffcalc = 0   ! calculation of flash frequ. 
  !                                 ! (IPARAM_*, see above)
  !
  INTEGER         :: i_iccg   = 0   ! flash handling
  !                                 ! 1: CG only
  !                                 ! 2: CG + IC
  INTEGER         :: i_shape  = 0   ! vertical distribution
  !                                 ! 1: const. mix. ratio
  !                                 ! 2: C-shape according to Pickering
  !
  ! scaling of flash frequency required for 'coupling' of different
  ! convection schemes or base models or model resolutions or ...:
  REAL(DP), DIMENSION(NPARAM) :: r_scal_ff = 1.0_dp
  !
  ! average NOx production per CG flash [kgN/flash]
  ! (P&R: 6.7E26 molec/flash -> ~ 15.6 kgN/flash):
  REAL(DP), DIMENSION(NPARAM) :: r_noxpf  = 15.6_dp
  !
  ! NOx production efficiency ratio (IC/CG):
  REAL(DP), DIMENSION(NPARAM) :: r_eff    = 0.1_dp 
  !

  ! DEFAULT PARAMETERS FOR IPARAM_Grewe
  ! (Grewe-parameterisation, for Tiedtke convection scheme) 
  ! ff =  A * (sqrt(D)*MF)**B
  REAL(DP) :: r_Grewe_A = 1.54e-5_dp
  REAL(DP) :: r_Grewe_B = 4.9_dp

  INTERFACE search_level_index
     MODULE PROCEDURE search_level_index_sigma
     MODULE PROCEDURE search_level_index_press
  END INTERFACE

CONTAINS

  ! =======================================================================
  ELEMENTAL SUBROUTINE ff_PaR_cth(ff, cth, lfr)

    IMPLICIT NONE
    
    ! I/O
    REAL(DP), INTENT(OUT) :: ff   ! total flash frequency  [1/s]
    REAL(DP), INTENT(IN)  :: cth  ! cloud top height [m]
    REAL(DP), INTENT(IN)  :: lfr  ! land fraction    [0...1]

    ! LOCAL
    REAL(DP) :: ff_land, ff_ocean, cth_km
    
    ! convert cth from m to km
    cth_km = cth * 1.e-3_dp

    ff_land  = 3.44e-5_dp * (cth_km**4.90_dp) 
    ff_ocean = 6.40e-4_dp * (cth_km**1.73_dp)

    ff = ( lfr * ff_land + (1._dp - lfr) * ff_ocean ) 

    ! NOTE: ff = ff*C with C = 0.97241*exp(0.048203*R)
    !                 and  R = delta-lat * delta-lon [deg^2]
    ! for resolution independene is omitted here

    ! convert from 1/min to 1/s
    ff = ( ff / 60.0_dp ) * r_scal_ff(IPARAM_PaR_T)

  END SUBROUTINE ff_PaR_cth
  ! =======================================================================

  ! =======================================================================
  ELEMENTAL SUBROUTINE ff_Grewe_muv(ff, cth, cbh, muv)

    IMPLICIT NONE
    INTRINSIC :: SQRT

    ! I/O
    REAL(DP), INTENT(OUT) :: ff   ! flash frequency  [1/s]
    REAL(DP), INTENT(IN)  :: cth  ! cloud top height [m]
    REAL(DP), INTENT(IN)  :: cbh  ! cloud bottom height [m]
    REAL(DP), INTENT(IN)  :: muv  ! mean updraft velocity [m/s]
!!$    REAL(DP), INTENT(IN)  :: A    ! correlation ...
!!$    REAL(DP), INTENT(IN)  :: B    ! ... parameter

    ff = ( r_Grewe_A * (muv * SQRT(cth-cbh))**r_Grewe_B )
    ! cpnvert from 1/min to 1/s
    ff = ( ff / 60.0_dp ) * r_scal_ff(IPARAM_Grewe)

  END SUBROUTINE ff_Grewe_muv
  ! =======================================================================

  ! =======================================================================
  ELEMENTAL SUBROUTINE ffcg_AaP_umf(ff, umf, area)

    IMPLICIT NONE
    INTRINSIC :: MIN

    ! I/O
    REAL(DP), INTENT(OUT) :: ff   ! flash frequency (cloud-to-ground)  [1/s]
    REAL(DP), INTENT(IN)  :: umf  ! (convective) updraft mass flux [kg/m^2/s]
    !                             ! at sigma = p/ps = 0.44 (or at p=440 hPa ?)
    REAL(DP), INTENT(IN)  :: area ! grid box area [m^2]

    ! LOCAL
    REAL(DP) :: mf

    ! kg/m^2/min ?, limit to 9.6
    mf = MIN(umf * 60._dp, 9.6_dp)

    ff = -3.71e-2_dp * mf
    ff = (ff + 5.23e-1_dp) * mf
    ff = (ff - 7.19e-1_dp) * mf
    ff = (ff + 3.08e-1_dp) * mf
    ff = (ff - 2.34e-2_dp) 

    ! rescaling with reference grid box area:
    ! reference: a box of 2? x 2.5? centered at 30?N: 5.353e10_dp m^2
    ff = ff * area / 5.353e10_dp

    ! 1/min -> 1/s
    ff = (ff / 60.0_dp) * r_scal_ff(IPARAM_AaP_M)

  END SUBROUTINE ffcg_AaP_umf
  ! =======================================================================

  ! =======================================================================
  ELEMENTAL SUBROUTINE ffcg_AaP_precip(ff, cprec, lfr, area)

    IMPLICIT NONE
    INTRINSIC :: MAX, MIN

    ! I/O
    REAL(DP), INTENT(OUT) :: ff     ! flash frequency  [1/s]
    REAL(DP), INTENT(IN)  :: cprec  ! convective precipitation [mm/s]
    REAL(DP), INTENT(IN)  :: lfr    ! land fraction    [0...1]
    REAL(DP), INTENT(IN)  :: area   ! grid box area [m^2]

    ! LOCAL
    REAL(DP) :: pr, ff_land, ff_ocean

    ! convert from mm/s to mm/day
    pr = cprec * 86400._dp

    ! limit to 7 ... 90 mm/day
    IF (pr < 7.0_dp) THEN 
       ff = 0.0_dp
    ELSE
       pr = MIN(pr, 90.0_dp)

       ! over land
       ff_land = -2.93e-6_dp * pr
       ff_land = (ff_land + 3.21e-4_dp) * pr
       ff_land = (ff_land + 5.41e-3_dp) * pr
       ff_land = (ff_land - 4.76e-2_dp) * pr
       ff_land = (ff_land - 3.75e-2_dp)

       ! over ocean and ice (!!!)
       ff_ocean = -2.42e-7_dp * pr
       ff_ocean = (ff_ocean + 3.68e-5_dp) * pr
       ff_ocean = (ff_ocean + 5.45e-3_dp) * pr
       ff_ocean = (ff_ocean - 4.80e-2_dp) * pr
       ff_ocean = (ff_ocean - 5.23e-2_dp)

       ff = ( lfr * ff_land + (1._dp - lfr) * ff_ocean )

       ! rescaling with reference grid box area:
       ! reference: a box of 2? x 2.5? centered at 30?N: 5.353e10_dp m^2
       ff = ff * area / 5.353e10_dp

    END IF

    ! 1/min -> 1/s
    ff = (ff / 60.0_dp) * r_scal_ff(IPARAM_AaP_P)

  END SUBROUTINE ffcg_AaP_precip
  ! =======================================================================

  ! =======================================================================
  ELEMENTAL SUBROUTINE ff_Finney_cif(ff, cif, lfr, area, which)

    IMPLICIT NONE

    ! I/O
    REAL(DP), INTENT(OUT) :: ff     ! flash frequency  [1/s]
    REAL(DP), INTENT(IN)  :: cif    ! convective cloud ice flux [kg/m^2/s]
    !                               ! at 440 hPa
    REAL(DP), INTENT(IN)  :: lfr    ! land fraction    [0...1]
    REAL(DP), INTENT(IN)  :: area   ! grid box area [m^2]
    INTEGER,  INTENT(IN)  :: which  ! IPARAM_FinIF or IPARAM_extIF

    ! LOCAL
    REAL(DP) :: ff_land, ff_ocean

    ff_land  = 6.58E-07_dp * cif ! flash density [1/m^2/s]
    ff_ocean = 9.08E-08_dp * cif ! flash density [1/m^2/s]

    ! * area: flash rate [1/s]
    ff = ( lfr * ff_land + (1._dp - lfr) * ff_ocean ) * area

    ff = ff * r_scal_ff(which)

  END SUBROUTINE ff_Finney_cif
  ! =======================================================================

  ! =======================================================================
  ELEMENTAL SUBROUTINE cg_fraction(pg, cdh)

    IMPLICIT NONE

    ! I/O
    REAL(DP), INTENT(OUT) :: pg     ! fraction of cloud-to-ground (CG) flashes
    REAL(DP), INTENT(IN)  :: cdh    ! cloud depth above 0 degC [m]

    ! LOCAL
    REAL(DP) :: zcdh_km ! cloud height below 0 degC [km]

    SELECT CASE(i_iccg)

    CASE(1)
       ! CG only
       pg = 1.0_dp

    CASE(2)
       ! IC + CG

       IF (cdh < 5500._dp) THEN
          ! only CG if 0 degC level is below 5.5 km
          pg=1._dp
       ELSE
          IF (cdh < 14000._dp) THEN
             ! 0 degC level is below 14 km
             zcdh_km = cdh * 1.e-3_dp
             pg = 0.021_dp * zcdh_km - 0.648_dp
             pg = pg * zcdh_km +  7.49_dp
             pg = pg * zcdh_km - 36.54_dp
             pg = pg * zcdh_km + 64.09_dp
             pg = 1._dp / pg
          ELSE
             ! constant fraction of 2% CG in clouds where
             ! 0 degC level is above 14 km
             pg = 0.02_dp
          END IF
       END IF

    END SELECT

  END SUBROUTINE cg_fraction
  ! =======================================================================

  ! =======================================================================
  ELEMENTAL SUBROUTINE ffs(IPAR, pg, ff, ffcg, ffic)

    IMPLICIT NONE
    
    ! I/O
    INTEGER,  INTENT(IN)    :: IPAR ! selected parameterisation
    REAL(DP), INTENT(IN)    :: pg   ! fraction of cloud-to-ground (CG) flashes
    REAL(DP), INTENT(INOUT) :: ff   ! total flash frequency [1/s]
    REAL(DP), INTENT(INOUT) :: ffcg ! cloud-to-ground flash frequency [1/s]
    REAL(DP), INTENT(INOUT) :: ffic ! intra-cloud flash frequency [1/s]

    SELECT CASE(IPAR)
    CASE(IPARAM_PaR_T, IPARAM_Grewe, IPARAM_FinIF, IPARAM_extIF)
       ffcg = ff * pg
       ffic = ff * (1._dp - pg)
    CASE(IPARAM_AaP_M, IPARAM_AaP_P)
       IF (pg > 0.0_dp) THEN
          ffic = ffcg / pg * (1._dp - pg)
       ELSE
          ! If the fraction of coud-to-ground flashes is zero, then
          ! ffcg is zero and thus, ffic must also be zero ...
          ffic = 0.0_dp
       ENDIF
       ff   = ffcg + ffic
    CASE DEFAULT
       !
    END SELECT    

  END SUBROUTINE ffs
  ! =======================================================================

  ! =======================================================================
  ELEMENTAL FUNCTION lflash(ff, area, coslat)

    ! CUT-OFF CRITERION:
    ! Flash frequencies below a certain cut-off will be ignored,
    ! except for in the scaling-mode!
    !
    !   NOTE: Originally, the criterion was 
    !         "minimum 1 flash per time step in grid box".
    !         This is, however, largely scale dependent and does not work
    !         for small scale models. Therefore, the criterion is replaced by
    !         "minimum fpm2s_min flashes per square-meter and second". 
    !
    ! An attempt to make this scale-independent is to normalise this condition
    ! per time and area. The gridbox area, however, depends on the latitude, 
    ! therefore it needs to be normalised to be approx. constant 
    ! (except near the poles).
    ! 
    ! ff * dtime > 1    !  / (grid-box-area / cos(lat)) / dtime
    ! ff / (grid-box-area / cos(lat)) > 1 / (grid-box-area / cos(lat)) / dtime
    ! ff / (grid-box-area / cos(lat)) > fpm2s_min 
    ! with
    ! fpm2s_min = Reference [1 / (grid-box-area / cos(lat)) / dtime]
    ! Reference: T42 (~ 2.8?x2.8?); dtime = 12min = 720s
    ! Ferret:
    !   let a_norm = gboxarea[l=1,x=0]/cos(Y[gy=gboxarea]/180*3.1415)
    !   list 1/a_norm[y=@ave]/720
    !   => fpm2s_min = 1.433E-14 [1/(m^2 s)]

    ! I/O
    LOGICAL :: lflash
    REAL(DP), INTENT(IN) :: ff     ! flash frequency [1/s]
    REAL(DP), INTENT(IN) :: area   ! surface area [m^2]
    REAL(DP), INTENT(IN) :: coslat ! cos(latitude)

    lflash = (ff / (area/coslat) >= fpm2s_min ) .OR. l_mode_scal

  END FUNCTION lflash
  ! =======================================================================

  ! =======================================================================
  SUBROUTINE cloud_heights(kvec, klev, deltaz, alti, temp, ibot, itop &
       , cth, cbh, czh, izero, ldeep, rthick)

    USE messy_main_constants_mem, ONLY: g

    IMPLICIT NONE

    ! I/O
    INTEGER,                           INTENT(IN)  :: kvec  ! vector length
    INTEGER,                           INTENT(IN)  :: klev  ! number of levels
    ! layer thinkness in [m]
    REAL(DP), DIMENSION(kvec, klev),   INTENT(IN)  :: deltaz
    ! altitude above ground at vertical box mid [m] 
    ! (Note: this is only required for backward compatibility, see below!)
    REAL(DP), DIMENSION(kvec, klev),   INTENT(IN)  :: alti
    ! temperature [K]
    REAL(DP), DIMENSION(kvec, klev),   INTENT(IN)  :: temp
    INTEGER,  DIMENSION(kvec),         INTENT(IN)  :: ibot  ! bottom level
    INTEGER,  DIMENSION(kvec),         INTENT(IN)  :: itop  ! top level
    !
    REAL(DP), DIMENSION(kvec),         INTENT(OUT) :: cth   ! top height [m]
    REAL(DP), DIMENSION(kvec),         INTENT(OUT) :: cbh   ! bottom height [m]
    ! depth of cloud "above" 0 degC level (i.e., T < 0 degC) [m]
    REAL(DP), DIMENSION(kvec),         INTENT(OUT) :: czh
    ! 0 degC [level]
    INTEGER, DIMENSION(kvec),          INTENT(OUT) :: izero
    ! flag for "deep" clouds
    LOGICAL, DIMENSION(kvec),          INTENT(OUT) :: ldeep
    ! thickness definition of "deep" clouds
    REAL(DP),                          INTENT(IN), OPTIONAL :: rthick

    ! LOCAL
    INTEGER  :: jv, jk
    REAL(DP) :: grheight(kvec, klev)   ! height of box [m]
    REAL(DP) :: zrthick

    IF (PRESENT(rthick)) THEN
       zrthick = rthick
    ELSE
       zrthick = 3000.0_dp  ! default: 3000 m
    END IF

    ! INIT
    cth(:)      = 0._dp
    cbh(:)      = 0._dp
    czh(:)      = 0._dp
    ldeep(:)    = .FALSE.

    grheight(:,1) = 0.0_dp ! in reality, this would be infinity
    DO jk=2, klev
       grheight(:,jk) = deltaz(:,jk)
    END DO

    DO jv=1, kvec
       ! top height: integrate from ground (klev) upwards (-1) to top level
       IF (itop(jv) > 0) THEN
          DO jk=klev,itop(jv),-1
             cth(jv) = cth(jv) + grheight(jv,jk)
          END DO
       END IF

       ! bottom height: integrate from ground (klev) upwards (-1) to 
       ! bottom level
       IF (ibot(jv) > 0) THEN
          DO jk=klev,ibot(jv),-1
             cbh(jv) = cbh(jv) + grheight(jv,jk)
          END DO
       END IF

       ! level of zero degrees
       izero(jv) = itop(jv)   
       IF ((itop(jv) > 0) .AND. (ibot(jv) > 0)) THEN
          DO jk=itop(jv),ibot(jv)
             ! height above 0 degC level
             IF (temp(jv,jk) <= 273.15_dp) THEN
                czh(jv) = czh(jv) + grheight(jv,jk)
                izero(jv) = jk ! 0 degC level
             END IF
          END DO
       ELSE
          czh(jv)      = 0._dp 
          izero(jv)    = 0
       END IF

    END DO

! NOTE: This is the consistent solution where only geopoti is used!
!!$    ldeep(:) = &
!!$         (itop(:) > 0) .AND.             &
!!$         (itop(:) < ibot(:)) .AND.       &
!!$         ( (cth(:) - cbh(:)) > zrthick)

! However, ...
    ! Note: In order to achieve identical results as in the old code,
    !       we use here the geopotential at box mid and approximate
    !       g =~ 10 m/s^2
    !  changed to altitude
    DO jv=1, kvec
       ldeep(jv) = &
            (itop(jv) > 0) .AND.             &
            (itop(jv) < ibot(jv))
       
       ! Note: without this separation, some compilers will
       !       access array elements out of bounds, even if the
       !       first and second condition above are false
       IF (ldeep(jv)) THEN
          ldeep(jv) = ldeep(jv) .AND. &
               ( (alti(jv,itop(jv)) - alti(jv,ibot(jv))) > zrthick )
      END IF
    END DO

  END SUBROUTINE cloud_heights
  ! =======================================================================

  ! =======================================================================
  SUBROUTINE search_level_index_sigma(k, kvec, klev, val, press, ps)

    IMPLICIT NONE

    ! I/O
    INTEGER,                         INTENT(IN)  :: kvec  ! vector length
    INTEGER,                         INTENT(IN)  :: klev  ! number of levels
    INTEGER,  DIMENSION(kvec),       INTENT(OUT) :: k     ! level index
    REAL(DP),                        INTENT(IN)  :: val   ! sigma value (p/ps)
    REAL(DP), DIMENSION(kvec, klev), INTENT(IN)  :: press ! pressure [Pa, hPa]
    ! surface pressure [Pa, hPa]
    REAL(DP), DIMENSION(kvec),       INTENT(IN)  :: ps

    ! LOCAL
    INTEGER  :: jv, jk

    DO jv=1, kvec
       k(jv) = klev
       DO jk=klev-1,1,-1
          IF ( press(jv,jk)/ps(jv) < val ) THEN
             k(jv) = jk+1
             EXIT
          END IF
       END DO
    END DO

  END SUBROUTINE search_level_index_sigma
  ! =======================================================================

  ! =======================================================================
  SUBROUTINE search_level_index_press(k, kvec, klev, val, press)

    IMPLICIT NONE

    ! I/O
    INTEGER,                        INTENT(IN)  :: kvec  ! vector length
    INTEGER,                        INTENT(IN)  :: klev  ! number of levels
    INTEGER,  DIMENSION(kvec),      INTENT(OUT) :: k     ! level index
    ! pressure value [Pa, hPa]
    REAL(DP),                       INTENT(IN)  :: val
    REAL(DP), DIMENSION(kvec,klev), INTENT(IN)  :: press ! pressure [Pa, hPa]

    ! LOCAL
    INTEGER  :: jv, jk

    DO jv=1, kvec
       k(jv) = klev
       DO jk=klev-1,1,-1
          IF ( press(jv,jk) < val ) THEN
             k(jv) = jk+1
             EXIT
          END IF
       END DO
    END DO

  END SUBROUTINE search_level_index_press
  ! =======================================================================

  ! =======================================================================
  SUBROUTINE mean_updraft_vel(kvec, klev, deltaz, umflx, rho, ibot, itop &
       , muv)

    USE messy_main_constants_mem, ONLY: g

    IMPLICIT NONE

    ! I/O
    INTEGER                                        :: kvec ! vector lenght
    INTEGER                                        :: klev ! number of levels
    ! layerthickness of vertical box [m] (see below)
    REAL(DP), DIMENSION(kvec, klev),   INTENT(IN)  :: deltaz
    ! updraft mass flux [kg m-^2 s^-1]
    REAL(DP), DIMENSION(kvec, klev),   INTENT(IN)  :: umflx
    ! air density [kg/m^3]
    REAL(DP), DIMENSION(kvec, klev),   INTENT(IN)  :: rho
    INTEGER,  DIMENSION(kvec),         INTENT(IN)  :: ibot  ! bottom level
    INTEGER,  DIMENSION(kvec),         INTENT(IN)  :: itop  ! top level
    ! mean updraft velocity [m/s]
    REAL(DP), DIMENSION(kvec),         INTENT(OUT) :: muv
    
    ! LOCAL
    INTEGER  :: jv, jk
    REAL(DP) :: grheight(kvec, klev)   ! height of box [m]    
    REAL(DP) :: cth, cbh

    grheight(:,1) = 0.0_dp ! in reality, this would be infinity
    DO jk=2, klev
       grheight(:,jk) = deltaz(:,jk)
    END DO

    muv(:)   = 0.0_dp

    DO jv = 1, kvec

       cth = 0.0_dp
       cbh = 0.0_dp

       IF ( (itop(jv) > 0) .AND. (ibot(jv) > 0) ) THEN
          DO jk = klev, itop(jv), -1
             cth = cth + grheight(jv,jk)
          END DO

          DO jk = klev, ibot(jv), -1
             cbh = cbh + grheight(jv,jk)
          END DO

          DO jk = itop(jv), ibot(jv)
             muv(jv)   = muv(jv) + ( umflx(jv,jk) / &
                  rho(jv,jk) ) * grheight(jv,jk)
          END DO
          IF ( (cth-cbh) > 0.0_dp ) THEN
             muv(jv) = muv(jv) / (cth-cbh)
          ELSE
             muv(jv) = 0.0_dp
          ENDIF
       END IF
    END DO

  END SUBROUTINE mean_updraft_vel
  ! =======================================================================

  ! =======================================================================
  ELEMENTAL SUBROUTINE nox_prod(IPAR, dtime, ffcg, ffic, NOx_cg, NOx_ic)

    IMPLICIT NONE

    ! I/O
    INTEGER,  INTENT(IN)    :: IPAR     ! selected parameterisation
    REAL(DP), INTENT(IN)    :: dtime    ! time interval [s]
    REAL(DP), INTENT(IN)    :: ffcg     ! cloud-to-ground flash frequency [1/s]
    REAL(DP), INTENT(IN)    :: ffic     ! intra-cloud flash frequency [1/s]
    REAL(DP), INTENT(OUT)   :: NOx_cg   ! NOx production [kg(N)]
    REAL(DP), INTENT(OUT)   :: NOx_ic   ! NOx production [kg(N)]

    NOx_cg = r_noxpf(IPAR) * dtime * ffcg
    NOx_ic = r_noxpf(IPAR) * dtime * ffic * r_eff(IPAR)

  END SUBROUTINE nox_prod
  ! =======================================================================

  ! =======================================================================
  SUBROUTINE shape_flat(kvec, klev, grmass, itop, ibot, izero, lland &
       , wcg, wic)

    USE messy_main_constants_mem, ONLY: M_air, MN, N_A

    IMPLICIT NONE

    INTRINSIC :: SUM

    ! I/O
    INTEGER,                         INTENT(IN)  :: kvec  ! vector length
    INTEGER,                         INTENT(IN)  :: klev  ! number of levels
    ! (dry) air mass in grid-box [kg]
    REAL(DP), PARAMETER                          :: scal = M_air/MN
    REAL(DP), DIMENSION(kvec, klev), INTENT(IN)  :: grmass ! box mass [kg]
    INTEGER,  DIMENSION(kvec),       INTENT(IN)  :: itop   ! top level
    INTEGER,  DIMENSION(kvec),       INTENT(IN)  :: ibot   ! bottom level
    INTEGER,  DIMENSION(kvec),       INTENT(IN)  :: izero  ! level of 0 degC 
    LOGICAL,                         INTENT(IN)  :: lland  ! land or sea/ice
    REAL(DP), DIMENSION(kvec, klev), INTENT(OUT) :: wcg    ! see below
    REAL(DP), DIMENSION(kvec, klev), INTENT(OUT) :: wic    ! ...

    ! LOCAL
    INTEGER  :: jv, jk, ibase
    REAL(dp) :: mass   ! air mass in convective column [kg]

    ! INIT
    wcg(:,:) = 0.0_dp
    wic(:,:) = 0.0_dp

    DO jv=1, kvec

       ! cloud-to-ground (top level of cloud to ground level) over land
       !                 (top level of cloud to bottom level of cloud
       !                  over sea and ice ?)
       IF (lland) THEN
          ibase = klev
       ELSE
          ibase = ibot(jv)
       ENDIF

       IF ( (itop(jv) > 0) .AND. (ibase > 0) ) THEN
          mass = SUM(grmass(jv,itop(jv):ibase))

          ! wcg is the "weight" in [mol(N)/mol(air) / kg(N)], i.e.
          ! multiplied with kg(N)/dtime it yields the tendency in [mol/mol/s] 
          IF (mass > 0.0_dp) &
               wcg(jv,itop(jv):ibase) = 1.0_dp/mass * scal
       ELSE
          wcg(jv,:) = 0.0_dp
       END IF

       ! intra-cloud (top level of cloud to level of 0 degC in cloud)
       ibase = izero(jv)

       IF ( (itop(jv) > 0) .AND. (ibase > 0) ) THEN
          mass = SUM(grmass(jv,itop(jv):ibase))
          ! wic ...
          IF (mass > 0.0_dp) &
               wic(jv,itop(jv):ibase) = 1.0_dp/mass * scal
       ELSE
          wic(jv,:) = 0.0_dp
       END IF

    END DO

  END SUBROUTINE shape_flat
  ! =======================================================================

  ! =======================================================================
  SUBROUTINE shape_C(kvec, klev, grmass, grvol, itop, ibot, lland, w)

    USE messy_main_constants_mem, ONLY: M_air, MN, N_A

    IMPLICIT NONE

    ! I/O
    INTEGER,                         INTENT(IN)  :: kvec  ! vector length
    INTEGER,                         INTENT(IN)  :: klev  ! number of levels
    ! (dry) air mass in grid-box [kg]
    REAL(DP), PARAMETER                          :: scal = M_air/MN
    REAL(DP), DIMENSION(kvec, klev), INTENT(IN)  :: grmass ! box mass [kg]
    REAL(DP), DIMENSION(kvec, klev), INTENT(IN)  :: grvol  ! box volume [m^3]
    INTEGER,  DIMENSION(kvec),       INTENT(IN)  :: itop   ! top level
    INTEGER,  DIMENSION(kvec),       INTENT(IN)  :: ibot   ! bottom level
    LOGICAL,                         INTENT(IN)  :: lland  ! land or sea/ice
    REAL(DP), DIMENSION(kvec, klev), INTENT(OUT) :: w      ! see below

    ! LOCAL
    INTEGER  :: jv, jk, ibase
    REAL(dp) :: mass   ! air mass in convective column [kg]
    REAL(dp) :: vol    ! volume of convective column [m^3]
    REAL(dp) :: rfa, rfb, rfc
    INTEGER  :: lbase, ltop
    REAL(dp), DIMENSION(klev) :: zw
    REAL(dp) :: zlevt

    ! INIT
    w(:,:) = 0.0_dp

    IF (lland) THEN ! land
       rfa =  0.12_dp
       rfb =  1.20_dp - 16._dp*rfa
       rfc = -2.80_dp + 39._dp*rfa
       lbase =  1
       ltop  = 13
    ELSE ! sea or ice
       rfa = 14.5_dp/81._dp
       rfb = (14.5_dp - 117._dp*rfa)/9._dp
       rfc = (198._dp*rfa - 11._dp)/9._dp
       lbase =  2
       ltop  = 11
    END IF

    DO jv=1, kvec

       IF (lland) THEN
          ibase = klev
       ELSE
          ibase = ibot(jv)
       ENDIF

       zw(:) = 0.0_dp
       vol = 0.0_dp
       IF ( (itop(jv) > 0) .AND. (ibase > 0) ) THEN
          one_box: IF (itop(jv) == ibase) THEN
             ! cloud is within one box
             zw(ibase) = 1.0_dp
             vol = grvol(jv,ibase)
          ELSE
             DO jk=itop(jv),ibase
                ! Note: linearization of model levels using the
                !       level index (!), such that always
                !       zbase <= zlevt <= ztop [km]
                ! 
                zlevt=REAL(ltop-lbase,dp)/REAL(itop(jv)-ibase,dp) &
                     *REAL(jk-ibase,dp)+REAL(lbase,dp) ! [km]
                ! 
                ! parabola parameterization (C-shape) using the
                ! parameters rfa, rfb, rfc of Pickering
                zw(jk) = (rfa*zlevt*zlevt+rfb*zlevt+rfc) ! [1]
                vol = vol + zw(jk)*grvol(jv,jk)
             END DO
          END IF one_box
       
          ! w is the "weight" in [mol(N)/mol(air) / kg(N)], i.e.
          ! multiplied with kg(N)/dtime it yields the tendency in [mol/mol/s]
          IF (vol > 0.0_dp) THEN
             DO jk=itop(jv),ibase
                w(jv,jk) = (zw(jk) * grvol(jv,jk) / vol) / grmass(jv,jk) * scal
             END DO
          END IF
       ELSE
          w(jv,:) = 0.0_dp
       ENDIF

   END DO

  END SUBROUTINE shape_C
  ! =======================================================================

  ! =========================================================================
  SUBROUTINE lnox_read_nml_ctrl(status, iou)

    USE messy_main_tools, ONLY: read_nml_open, read_nml_check, read_nml_close

    IMPLICIT NONE

    ! I/O
    INTEGER, INTENT(OUT) :: status ! error status
    INTEGER, INTENT(IN)  :: iou    ! logical I/O unit

    NAMELIST /CTRL/ l_mode_scal &
         , i_ffcalc, i_iccg, i_shape, r_scal_ff, r_noxpf, r_eff &
         , r_Grewe_A, r_Grewe_B

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER       :: substr='lnox_read_nml_ctrl'
    LOGICAL                           :: lex          ! file exists ?
    INTEGER                           :: fstat        ! file status

    ! INITIALIZE
    status = 1 ! ERROR

    CALL read_nml_open(lex, substr, iou, 'CTRL', modstr)
    IF (.not.lex) RETURN    ! <modstr>.nml does not exist

    READ(iou, NML=CTRL, IOSTAT=fstat)
    CALL read_nml_check(fstat, substr, iou, 'CTRL', modstr)
    IF (fstat /= 0) RETURN  ! error while reading namelist

    ! DIAGNOSE NAMELIST AND SET GLOBAL SWITCHES
    !
    IF (l_mode_scal) THEN
       WRITE(*,*) 'WARNING: LNOX-Mode: diagnostic to assess scaling parameters'
    ELSE
       WRITE(*,*) 'LNOX-Mode: production simulation'
    ENDIF
    !
    WRITE(*,*) 'Lightning-NOx emission parameterization:'
    !
    WRITE(*,*) 'Calculation of flash frequency ...'
    SELECT CASE(i_ffcalc)
    CASE(IPARAM_PaR_T)
       WRITE(*,*) '... Price and Rind, 1992'
    CASE(IPARAM_Grewe)
       WRITE(*,*) '... Grewe et al., 2001'
    CASE(IPARAM_AaP_M)
       WRITE(*,*) '... Allen and Pickering (Massflux), 2002'
    CASE(IPARAM_AaP_P)
       WRITE(*,*) '... Allen and Pickering (Precipitation), 2002'
    CASE(IPARAM_FinIF)
       WRITE(*,*) '... Finney et al. (Cloud ice flux), 2014'
    CASE(IPARAM_extIF)
       WRITE(*,*) '... Finney et al. (Cloud ice flux), 2014, extended'
    CASE(IPARAM_DahlC)
       WRITE(*,*) '... Dahl, 2010 (COSMO only!)'
    CASE(IPARAM_ALL )
       WRITE(*,*) '... All parameterisations are used diagnostically:'
       WRITE(*,*) '... (1) Price and Rind, 1992'
       WRITE(*,*) '... (2) Grewe et al., 2001'
       WRITE(*,*) '... (3) Allen and Pickering (Massflux), 2002'
       WRITE(*,*) '... (4) Allen and Pickering (Precipitation), 2002'
       WRITE(*,*) '... (5) Finney et al. (Cloud ice flux), 2014'
       WRITE(*,*) '... (6) Dahl, 2010 (COSMO only!)'
       WRITE(*,*) '... (7) Finney et al. (Cloud ice flux), 2014, extended'
    CASE DEFAULT
       WRITE(*,*) 'ERROR: Unknown i_ffcalc!'
       RETURN
    END SELECT
    
    IF (i_ffcalc /= IPARAM_ALL ) THEN
       WRITE(*,*) 'Rescaling factor for flash frequency ...'
       WRITE(*,*) '... ', r_scal_ff(i_ffcalc)
       WRITE(*,*) 'Average NOx production per CG flash ...'
       WRITE(*,*) '... ',r_noxpf(i_ffcalc),'[kgN/flash]'
    ELSE
       WRITE(*,*) 'Rescaling factors for flash frequency ...'
       WRITE(*,*) '... ', r_scal_ff(:)
       WRITE(*,*) 'Average NOx production per CG flash ...'
       WRITE(*,*) '... ',r_noxpf(:),'[kgN/flash]'
    END IF

    WRITE(*,*) 'Flash types ...'
    SELECT CASE(i_iccg)
    CASE(1)
       WRITE(*,*) '... cloud-to-ground flashes only'
    CASE(2)
       IF (i_ffcalc /= 9) THEN
          WRITE(*,*) '... cloud-to-ground and intra-cloud flashes'
          WRITE(*,*) 'NOx production efficiency ratio (IC/CG) ...'
          WRITE(*,*) '... ',r_eff(i_ffcalc)
       ELSE
          WRITE(*,*) '... cloud-to-ground and intra-cloud flashes'
          WRITE(*,*) 'NOx production efficiency ratios (IC/CG) ...'
          WRITE(*,*) '... ',r_eff(:)
       END IF
    CASE DEFAULT
       WRITE(*,*) 'ERROR: Unknown i_iccg!'
       RETURN
    END SELECT
    !
    WRITE(*,*) 'Vertical NOx distribution ...'
    SELECT CASE(i_shape)
    CASE(1)
       WRITE(*,*) '... constant mixing ratio'
    CASE(2)
       WRITE(*,*) '... C-shape, Pickering et al., 1998'
    CASE DEFAULT
       WRITE(*,*) 'ERROR: Unknown i_shape!'
       RETURN
    END SELECT

    CALL read_nml_close(substr, iou, modstr)
    status = 0 ! NO ERROR

  END SUBROUTINE lnox_read_nml_ctrl
  ! =========================================================================

! ***********************************************************************
END MODULE messy_lnox
! ***********************************************************************

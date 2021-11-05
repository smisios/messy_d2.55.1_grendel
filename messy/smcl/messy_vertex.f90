MODULE MESSY_VERTEX
  ! AUTHOR:
  !  H.G. Ouwersloot, MPIC, April 2016

  USE messy_main_constants_mem, ONLY: dp,  ckap=>c_vKar, cvdifts, rd, rv

  IMPLICIT NONE

  PRIVATE
  CHARACTER(len=*), PUBLIC, PARAMETER :: MODSTR='vertex'
  CHARACTER(LEN=*), PUBLIC, PARAMETER :: modver='1.0'

  SAVE

  !------------------------------------------------------------------------------
  ! Constants previously defined in E5VDIFF
  !------------------------------------------------------------------------------
  REAL(dp) :: clam    = 150._dp        !  *asymptotic mixing length for momentum.
  REAL(dp) :: cc      = 5._dp          !  *stability parameter for unstable cases.
  REAL(dp) :: cchar   = 0.018_dp       !  *charnock constant.
  REAL(dp) :: cfreec  = 0.001_dp       !  *free convection parameter
  REAL(dp) :: cgam    = 1.25_dp        !  *free convection parameter
  REAL(dp) :: cz0ice  = 0.001_dp       !  *roughness over sea-ice
  REAL(dp) :: csncri  = 5.85036E-3_dp  !  *critical snow depth for soil computations

 ! Constants used for computation of saturation mixing ratio
 ! over liquid water (*c_les*) or ice(*c_ies*)
 ! Based on Murray (Journal of Applied Meteorology, 1967)
  REAL(dp), PARAMETER :: c1es  = 610.78_dp
  REAL(dp), PARAMETER :: c2es  = c1es*rd/rv
  REAL(dp), PARAMETER :: c3les = 17.269_dp
  REAL(dp), PARAMETER :: c3ies = 21.875_dp
  REAL(dp), PARAMETER :: c4les = 35.86_dp
  REAL(dp), PARAMETER :: c4ies = 7.66_dp

  ! Variables for computation of evapotranspiration
  REAL(dp), PARAMETER :: cva   = 5000._dp
  REAL(dp), PARAMETER :: cvb   = 10._dp
  REAL(dp), PARAMETER :: cvc   = 100._dp
  REAL(dp), PARAMETER :: cvbc  = cvb*cvc
  REAL(dp), PARAMETER :: cvabc = (cva+cvbc)/cvc
  REAL(dp), PARAMETER :: cvk   = .9_dp
  REAL(dp), PARAMETER :: cvrad = 0.55_dp
  !------------------------------------------------------------------------------

  REAL(dp) :: zm1  = 0.92_dp                         ! Constant to calculate S_N_M - A1
  REAL(dp) :: zm2  = 16.6_dp                         ! Constant to calculate S_N_M - B1
  REAL(dp) :: zm4  = 0.08_dp                         ! Constant to calculate S_N_M - C1
  REAL(dp) :: zh1  = 0.74_dp                         ! Constant to calculate S_N_H - A2
  REAL(dp) :: zh2                                    ! Constant to calculate S_N_H - gamma_1
  REAL(dp) :: zshn                                   ! Neutral exchange coefficient for heat, S_N_H, used to calculate S_H
  REAL(dp) :: zsmn                                   ! Neutral exchange coefficient for momentum, S_N_M, used to calculate S_M
  REAL(dp) :: zplmin  = 0.35_dp                      ! Fraction of field capacity for wilting point, below which the soil moisture stress correction function = 0
  INTEGER :: irstom = 0 ! ju_te_20180625
  INTEGER :: izwet  = 0 ! ju_te_20180912
  INTEGER :: ifws   = 0 ! ju_te_20190920

  PUBLIC   :: dp

  PUBLIC   :: cc                                     ! Closure constant - 5
  PUBLIC   :: cchar                                  ! Charnock constant, used to calculate z0m over sea - 0.018
  PUBLIC   :: cfreec                                 ! Factor (beta) to calculate term for (C_R) for convective C_H over sea (Miller et al., 1992) - 0.001
  PUBLIC   :: cgam                                   ! Power (gamma) to calculate convective C_H over sea (Miller et al., 1992) - 1.25
  PUBLIC   :: ckap                                   ! Von Karman constant - 0.4
  PUBLIC   :: clam                                   ! Standard mixing length in the ABL - 150 m
  PUBLIC   :: csncri                                 ! Critical snow depth for soil computations - 5.85036 mm
  PUBLIC   :: cva                                    ! Constant to define stomatal resistance - 5000 J/m3
  PUBLIC   :: cvb                                    ! Constant to define stomatal resistance - 10 W/m2
  PUBLIC   :: cvabc                                  ! (cva+cvb*cvc)/cvc
  PUBLIC   :: cvc                                    ! Minum stomatal resistance - 100 s/m
  PUBLIC   :: cvk                                    ! Extinction coefficient - 0.9
  PUBLIC   :: cvdifts                                ! Factor for timestep weighting, alpha - 1.5 for the applied over-implicit scheme
  PUBLIC   :: cvrad                                  ! Fraction of the net SW radiation contributing to PAR - 0.55
  PUBLIC   :: cz0ice                                 ! Roughness over sea-ice - 0.001 m
  PUBLIC   :: c2es                                   ! Constant to compute saturation mixing ratio over liquid water or ice - 610.78*rd/rv
  PUBLIC   :: c3ies                                  ! Constant to compute saturation mixing ratio over ice - 21.875
  PUBLIC   :: c3les                                  ! Constant to compute saturation mixing ratio over liquid water - 17.269
  PUBLIC   :: c4ies                                  ! Constant to compute saturation mixing ratio over ice - 7.66
  PUBLIC   :: c4les                                  ! Constant to compute saturation mixing ratio over liquid water - 35.86
! op_pj_20161107: see above

  PUBLIC   :: zm1
  PUBLIC   :: zm2
  PUBLIC   :: zm4
  PUBLIC   :: zh1
  PUBLIC   :: zh2
  PUBLIC   :: zshn
  PUBLIC   :: zsmn
  PUBLIC   :: zplmin
  PUBLIC   :: irstom ! ju_te_20180625
  PUBLIC   :: izwet  ! ju_te_20180912
  PUBLIC   :: ifws   ! ju_te_20190920

  PUBLIC   :: pblheightRi

  ! ju_te_20180713+
  ! Temperature and drought stress for rstom (functions)
  PUBLIC :: fT, fD
  
  ! Subroutines
  PUBLIC :: vertex_read_nml_ctrl
  PUBLIC :: vertex_calc_rstom
  ! ju_te_20180713-
  PUBLIC :: vertex_surftemp

  CONTAINS

! ju_te_20180713+
!=========================================================================
   SUBROUTINE vertex_read_nml_ctrl(status, iou) 

    ! VERTEX MODULE ROUTINE (CORE)
    !
    ! READ VERTEX NAMELIST, CHECK IT, AND INITIALIZE GLOBAL VARIABLES
    !
    ! Author: Tamara Emmerichs, FZ Juelich, 26-06-2018 (ju_te_20180625)
    
    USE messy_main_tools, ONLY: read_nml_open, read_nml_check, read_nml_close
    IMPLICIT NONE
     
    ! I/O
    INTEGER, INTENT(OUT) :: status
    INTEGER, INTENT(IN)  :: iou   ! logical I/O unit
    
    ! LOCAL
    CHARACTER(LEN=*), PARAMETER       :: substr='vertex_read_nml_ctrl'
    LOGICAL                           :: lex          ! file exists ?
    INTEGER                           :: fstat        ! file status

    NAMELIST /CTRL/  irstom, izwet, ifws

    status = 1 ! ERROR ON RETURN
    CALL read_nml_open(lex, substr, iou, 'CTRL', modstr)
    IF (.not.lex) RETURN    ! <modstr>.nml does not exist
    READ(iou, NML=CTRL, IOSTAT=fstat)
    CALL read_nml_check(fstat, substr, iou, 'CTRL', modstr)
    IF (fstat /= 0) RETURN  ! error while reading namelist
    
    WRITE(*,*) 'Switch application of reduction factors to evapotranspiration on'
    SELECT CASE (izwet)
    CASE(0)
      WRITE(*,*) '... DEFAULT'
    CASE(1)
      WRITE(*,*) '... with stress factors - Jarvis 1976/Katul 2009'
    END SELECT

    WRITE(*,*) 'Choose method for rstom calculation '
    SELECT CASE (irstom)
    CASE(0)
       WRITE(*,*) '... DEFAULT (Ganzeveld et al., 1995)'
    CASE(1)
       WRITE(*,*) 'leaf scale resistance with variable LAI (Ganzeveld et al., 1995)'
    CASE(2)
       WRITE(*,*) 'canopy scale resistance with variable LAI (Ganzeveld et al., 1998)'
    CASE(3)
       WRITE(*,*) 'canopy scale resistance with temperature stress factor (Jarvis, 1976)'
    CASE(4)
       WRITE(*,*) 'canopy scale resistance with drought stress factor (Katul et al., 2009)'
    CASE(5)
       WRITE(*,*) 'canopy scale resistance with both stress factors (Jarvis, 1976/Katul et al., 2009)'
    END SELECT

    WRITE(*,*) 'Method used for fws (soil moisture stress function for DDEP and ET) calculation'

    SELECT CASE (ifws)
    CASE(0)
       WRITE(*,*) '... DEFAULT'
    CASE(1)
       WRITE(*,*) '... without wilting point according to the original parametrization of Delworth and Manabe (1988) '
    END SELECT

    ! DIAGNOSE NAMELIST AND SET GLOBAL SWITCHES
    CALL read_nml_close(substr, iou, modstr)

    status = 0

  END SUBROUTINE vertex_read_nml_ctrl
!=========================================================================
! ju_te_20180713-

!=========================================================================
  SUBROUTINE pblheightRi(pblh, pblh_i, pblh_f,                            &
                         kproma, klev,                                    &
                         theta, q, geopot, u, v,                          &
                         buoy, ustar, ricralt                            )

  USE messy_main_constants_mem, ONLY: g, zcrdq=>vtmpc1
  USE messy_main_tools,         ONLY: iso2ind

  IMPLICIT NONE

  ! ======================================================================
  !
  ! Routine to calculate boundary-layer height.
  ! Based on routine provided by Bert Holtslag, March 1990,
  ! (KNMI, de Bilt, the Netherlands)
  !
  ! Adjustments by Bert Holtslag (July '91) for countergradient terms,
  ! and temperature excesses. Mechanical mixing depth removed.
  !
  !    * updated to EC_PBLHGHT by E. van Meijgaard, KNMI, Januari  '93
  !    * updated to EC_PBLHGHT by E. van Meijgaard, KNMI, December '95
  !      in order to be callable by the ECHAM3
  !      local vertical diffusion scheme EC_VDIFF /EC4_VDIFF
  !    * F77 to F90 for ECHAM5 done by Michael Traub, MPCH, Nov 2003
  !
  !    Modifications by D. Vogelezang, May 1993:
  !    1. 'Parcel temperature' calculation at first model level
  !    2. 'Delta-U method' in  bound. lay. height formulation
  !    Modifications by D. Vogelezang, Fall 1994:
  !    3. 'Ustar-shear production' term added
  !
  !    Based on implementation in TROPOP routine
  !    Cleaned by removing calculations whose results did not impact output
  !    Code was optimized for vector machine - now uses IF-statements

  ! Output
  REAL(DP), DIMENSION(:),   INTENT(OUT)  :: pblh    ! Boundary-layer height [m]
  REAL(DP), DIMENSION(:),   INTENT(OUT)  :: pblh_i  ! Boundary-layer height index
  REAL(DP), DIMENSION(:),   INTENT(OUT)  :: pblh_f  ! Fraction of grid box at pblh_i within boundary layer

  ! Input
  INTEGER                 , INTENT(IN)   :: kproma  ! Vector length of current row
  INTEGER                 , INTENT(IN)   :: klev    ! Amount of vertical levels
  REAL(DP), DIMENSION(:,:), INTENT(IN)   :: theta   ! Potential temperature [K]
  REAL(DP), DIMENSION(:,:), INTENT(IN)   :: q       ! Specific humidity [kg/kg]
  REAL(DP), DIMENSION(:,:), INTENT(IN)   :: geopot  ! Geopotential height [m]
  REAL(DP), DIMENSION(:,:), INTENT(IN)   :: u       ! Windspeed Eastward direction [m/s]
  REAL(DP), DIMENSION(:,:), INTENT(IN)   :: v       ! Windspeed Northward direction [m/s]
  REAL(DP), DIMENSION(:),   INTENT(IN)   :: buoy    ! Surface buoyancy flux [Km/s]
  REAL(DP), DIMENSION(:),   INTENT(IN)   :: ustar   ! Surface friction velocity [m/s]
  REAL(DP), OPTIONAL,       INTENT(IN)   :: ricralt ! Critical Richardson number - if deviating from standard 0.25

  ! Local parameters
  REAL(DP), PARAMETER :: zacb  = 100._dp
  REAL(DP), PARAMETER :: onet  = 0.33333333_dp
  REAL(DP), PARAMETER :: c1    = 0.6_dp ! = 0.4 x 1.5 = 1.5 times von Karman
  REAL(DP), PARAMETER :: fak   = 8.5_dp

  ! Local variables
  REAL(DP) :: ricr = 0.25_dp
  INTEGER  :: nlevp, npbl, nlevc
  INTEGER  :: i, k
  REAL(DP) :: vvk, tkv, zrino
  REAL(DP) :: wsc

  LOGICAL,  DIMENSION(kproma) :: jwork
  INTEGER,  DIMENSION(kproma) :: pblhi
  REAL(DP), DIMENSION(kproma) :: tav
  REAL(DP), DIMENSION(kproma) :: tlv
  REAL(DP), DIMENSION(kproma) :: ustr
  REAL(DP), DIMENSION(kproma) :: wrino

  ! Init
  IF (PRESENT(ricralt)) THEN
     ricr = ricralt
  END IF
  nlevp = klev+1
  npbl  = klev-1

  ! Compute bottom level virtual potential temperature and Obukhov length

  DO  i=1,kproma
     tav(i)    = theta(i,klev)*(1._dp + zcrdq*q(i,klev))                         ! Theta_v at surface
     ustr(i)   = max(ustar(i),0.01_dp)                                           ! u* at surface
  END DO

  !  calculation of boundary layer height

  jwork  = .TRUE.                                                                ! Switch - if TRUE, keep looking for boundary-layer height; still within boundary layer
  wrino  = 0._dp                                                                 ! Richardson number at the previously evaluated level
  pblh   = 0._dp                                                                 ! Boundary-layer height

  DO k=2,npbl
     nlevc = nlevp-k                                                             ! Current level
     DO i=1,kproma
        IF (jwork(i)) THEN
           vvk   = (u(i,nlevc)-u(i,klev))**2 + (v(i,nlevc)-v(i,klev))**2 +     & ! Corrected squared difference in U between current level and lowest level
                   zacb*ustr(i)*ustr(i)
           tkv   = theta(i,nlevc)*(1. + zcrdq*q(i,nlevc))                        ! Theta_v at current level
           zrino = (tkv-tav(i))*(geopot(i,nlevc)-geopot(i,klev))/(tav(i)*vvk)    ! Bulk Richardson number between current level and lowest level
           IF (zrino .GE. ricr) THEN                                             ! First level for which Ri >= Ri{critical}
              pblh(i) = (geopot(i,nlevc+1) + ((ricr-wrino(i))/(zrino-wrino(i)))& ! Boundary-layer height determined by interpolation to height at which Ri = Ri{critical}
                        * (geopot(i,nlevc)-geopot(i,nlevc+1)))/g                 !   Based on height and Ri of previously evaluated level (underlying adjacent level)
              jwork(i) = .FALSE.                                                 ! Set switch to .FALSE.; from here on located above the boundary layer
           ELSE
              wrino(i) = zrino                                                   ! Save Richardson number; to be used if next evaluated level is located above the boundary layer
           END IF
        END IF
     END DO
  END DO

  DO i=1,kproma
     IF (jwork(i)) THEN
        pblh(i)   = geopot(i,2) / g                                              ! If not found, boundary-layer height is at second highest grid center
     END IF
     IF (buoy(i) .GT. 0._dp) THEN                                                ! Unstable case; account for thermal excess
        jwork(i)  = .TRUE.                                                       ! For this position, re-evaluate the boundary-layer height
        wsc       = (ustr(i)**3 + c1 * (g/tav(i)) * pblh(i) * buoy(i))**onet     ! Turbulent velocity scale: w_m = (u*^3 + w*^3)^(1/3) (Holtslag et al., 1995)
        tlv(i)    = tav(i) + fak * buoy(i) / wsc                                 ! Virtual potential temperature at surface adapted for thermal excess
     ENDIF
  END DO

  wrino = 0._dp                                                                  ! Reset Richardson number at the previously evaluated level

  !  Improve phbl estimate under convective conditions using convective temperature excess

  DO k=2,npbl
     nlevc = nlevp-k                                                             ! Current level
     DO i=1,kproma
        IF (jwork(i)) THEN
           vvk   = (u(i,nlevc)-u(i,klev))**2 + (v(i,nlevc)-v(i,klev))**2 +     & ! Corrected squared difference in U between current level and lowest level
                   zacb*ustr(i)*ustr(i)
           tkv   = theta(i,nlevc)*(1. + zcrdq*q(i,nlevc))                        ! Theta_v at current level
           zrino = (tkv-tlv(i))*(geopot(i,nlevc)-geopot(i,klev))/(tlv(i)*vvk)    ! Bulk Richardson number between current level and lowest level
           IF (zrino .GE. ricr) THEN                                             ! First level for which Ri >= Ri{critical}
              pblh(i) = (geopot(i,nlevc+1) + ((ricr-wrino(i))/(zrino-wrino(i)))& ! Boundary-layer height determined by interpolation to height at which Ri = Ri{critical}
                        * (geopot(i,nlevc)-geopot(i,nlevc+1)))/g                 !   Based on height and Ri of previously evaluated level (underlying adjacent level)
              jwork(i) = .FALSE.                                                 ! Set switch to .FALSE.; from here on located above the boundary layer
           ELSE
              wrino(i) = zrino                                                   ! Save Richardson number; to be used if next evaluated level is located above the boundary layer
           END IF
        END IF
     END DO
  END DO
  !
  DO i=1,kproma
     IF (jwork(i)) THEN
        pblh(i)   = geopot(i,2) / g                                              ! If not found, boundary-layer height is at second highest grid center
     END IF
  END DO

  CALL iso2ind(kproma, geopot/g, pblh, pblhi, f=pblh_f, lrev=.TRUE.)            ! Determine pblh_i and pblh_f
  pblh_i(1:kproma) = REAL(pblhi,dp)

  END SUBROUTINE pblheightRi
!=========================================================================

! ju_te_20180713+
!=========================================================================
! Temperature-stress factor for stomatal conductance according to Jarvis (1976)
  ELEMENTAL REAL(dp) FUNCTION fT(zptslm1)
    implicit none
    REAL(dp), INTENT(IN) :: zptslm1
    REAL(dp), PARAMETER  :: tlow = 268.15, to =   293.15, thigh =   318.15 !K        ! Jarvis (1976) Table 2
    REAL(dp)             :: b3, b4                                                   ! leaf adjustment according to Jarvis(1976)
    REAL(dp)             :: pTs                                                      ! temporary temperature variable

    ! define temporary temperature variable according to the parametrization by Jarvis (1976)
    pTs = MIN((thigh-0.15_dp),MAX(tlow+0.15_dp,zptslm1))

! leaf adjustment according to Jarvis(1976)
  b4 =   (thigh-to)/(thigh-tlow)
  b3 = 1._dp/((to-tlow)*(thigh-to)**b4)

! Temperature stress factors for the conductance
  fT = MIN(1._dp,MAX(0.1_dp,(b3*(pTs-tlow)*(thigh-pTs)**b4)))

  END FUNCTION fT

! Drought-stress factor for stomatal conductance according to  Katul et al. (2009)
  ELEMENTAL REAL(dp) FUNCTION fD(zptslm1,zrh_2m)
    implicit none
    REAL(dp), INTENT(IN) :: zptslm1, zrh_2m
    REAL(dp) :: D                                                                    ! calculated according to Kraus (2004)

    D = (1._dp-MIN(0.99_dp,zrh_2m)) * &
        (0.61078_dp*EXP((17.1_dp*(zptslm1-273.15_dp))/(235._dp+zptslm1-273.15_dp)))  ! D=(1-RH) * esat(T) (Vapour pressure deficit)
    fD = (MAX(0.1_dp,D))**(-0.5_dp)
  END FUNCTION fD

! ju_te_20180613-

  SUBROUTINE vertex_calc_rstom(kproma,zwet, fws, rleaf, pws, pwsmx, &
                               zsrfll, zlai,ptslm1,                 &
                               rh_2m)

  !==================================================================
  ! Author: Tamara Emmerichs, FZ Juelich, 28-06-2018 
  ! t.emmerichs@fz-juelich.de
  ! vertex_calc_rstom to calculate the stomatal resistance 
  ! 0)the default: according to Ganzeveld et al., 
  !   JOURNAL OF GEOPHYSICAL RESEARCH, VOL. 100, NO. D10, 
  !   PAGES 20,999-21,012, OCTOBER 20, 1995
  ! 1)with variable LAI, at leaf scale
  ! 2)with variable LAI, at canopy scale: according to Ganzeveld et al.,
  !   JOURNAL OF GEOPHYSICAL RESEARCH, VOL. 103, NO. D5, 
  !   PAGES 5679-5694, MARCH 20, 1998
  ! 3)CASE 2+ temperature stress factor: according to Jarvis
  !   Phil. Trans. R. Soc. Lond. B. 273, 593-610 (1976)
  ! 4)CASE 2+ drought stress factor: according to Jarvis/Katul
  !   Plant, Cell and Environment (2009) 32, 968-979
  ! 5)all CASES together
  ! =================================================================

  ! Output
  REAL(dp), DIMENSION(:), INTENT(OUT)  :: zwet                                   ! Stomatal resistance [s/m]
  REAL(dp), DIMENSION(:), INTENT(OUT)  :: fws                                    ! soil moisture stress function
  REAL(dp), DIMENSION(:), INTENT(OUT)  :: rleaf                                  ! r_stom*fws

  ! Input
  REAL(dp), DIMENSION(:), INTENT(IN)  :: pws                                     ! Surface soil wetness
  REAL(dp), DIMENSION(:), INTENT(IN)  :: pwsmx                                   ! Field capacity of the soil
  REAL(dp), DIMENSION(:), INTENT(IN)  :: zsrfll                                  ! Net solar radiation over land
  REAL(dp), DIMENSION(:), INTENT(IN)  :: zlai                                    ! Leaf Area Index
  INTEGER,                INTENT(IN)  :: kproma                                

  REAL(dp), OPTIONAL,  DIMENSION(:)  :: rh_2m, ptslm1
  
  ! Parameters
  REAL(dp), PARAMETER :: zplmax  = 0.75_dp                                       ! Fraction of field capacity above with fws = 1
  REAL(dp), PARAMETER :: zepsr = 1e-10_dp                                        ! Minimum radiation 10^-10 W/m2
  REAL(dp), PARAMETER :: zepevap = 1e-10_dp                                      ! Minimum conductivity of leaf stomata 10^-10 m/s

  ! Local
  INTEGER             :: jl
  REAL(dp)            :: zabcs, zrsi, zsrfl, zwcrit, zwpwp
  REAL(dp)            :: zln1, zln2
  REAL(dp)            :: tmp_rleaf  ! leaf adjustment according to Jarvis(1976)
  REAL(dp)            :: zzlai
  
  DO jl=1,kproma                                                                    ! Calculating stomatal resistance, r_stom, according to Eq. (8) of Ganzev        eld and Lelieveld (1995)
    zwcrit=zplmax*pwsmx(jl)                                                         ! Soil moisture content above which the soil moisture stress correction f        unction = 1
    zwpwp=zplmin*pwsmx(jl)                                                          ! Wilting point
    ! soil moisture stress function  ! ju_te_20191217+
    SELECT CASE(ifws)
    CASE(0) ! Default
       fws(jl)=MAX(0._dp,MIN(1._dp,(pws(jl)-zwpwp)/(zwcrit-zwpwp))) 
    CASE(1) ! original formular of fws independent of the wilting point according to Delworth, Thomas L., and Syukuro Manabe. JOURNAL OF CLIMATE 1.5 (1988): 523-547.
       IF (pws(jl).GT.zwcrit) THEN
          fws(jl)=1._dp
       ELSE IF (pws(jl).LE.zwcrit) THEN
          fws(jl)=pws(jl)/zwcrit
       END IF
    END SELECT
    ! ju_te_20191217-
   
    zsrfl=MAX(zepsr,zsrfll(jl)*cvrad)                                               ! Net PAR over land
    zabcs=(cva+cvb*cvc)/(cvc*zsrfl)                                                 ! Coefficient d * PAR
    zln1=LOG((zabcs*EXP(cvk*zlai(jl))+1._dp)/(zabcs+1._dp))                         ! First logarithm in Eq. (8)
    zln2=LOG((zabcs+EXP(-cvk*zlai(jl)))/(zabcs+1._dp))                              ! Second logarithm in Eq. (8)
    zrsi=(cvb*zln1/cvabc-zln2)/(cvk*cvc)                                            ! 1 / (r_stom * F(ws) )
   
    zzlai=MAX(1.e-5_dp,zlai(jl)) ! avoid negative values over the ocean
    ! reduction factors for evapotranspiration
    SELECT CASE (izwet)
    CASE(0) ! DEFAULT
       ! canopy scale resistance
       zwet(jl)=1._dp/(zrsi*fws(jl)+zepevap)
     CASE(1)
        ! adjustment according to Jarvis (1976)
        zwet(jl)=1._dp/(zrsi*fws(jl)+zepevap)*1._dp/fT(ptslm1(jl))* &
             1._dp/fD(ptslm1(jl),rh_2m(jl))
     END SELECT

    SELECT CASE (irstom)
      CASE(0) ! DEFAULT
        ! calculation of stomatal canopy scale resistance using the echam
        ! equation applying an LAI of 1, included by Laurens Ganzeveld, 18/10/01
        rleaf(jl) = cvk*cvc/(cvb*LOG((zabcs*EXP(cvk)+1._dp)                        & 
                 /(zabcs+1._dp))/cvabc-LOG((zabcs+EXP(-cvk))/(zabcs+1._dp)))        
      CASE(1)
        ! leaf scale resistance       
        rleaf(jl)=cvk*cvc/(cvb*LOG((zabcs*EXP(cvk*zzlai)+1._dp)                 &! the correct formula for the canopy-scale rco_leaf is given but then it is
            /(zabcs+1._dp))/cvabc-LOG((zabcs+EXP(-cvk*zzlai))/(zabcs+1._dp)))   &! scaled to leaf-level stomatal resistance. rcut and rc0x in DDEP 
             *zzlai                                                              ! are assumed to be leaf-level resistances which there divided by LAI
                  
      CASE(2)
        ! canopy scale resistance -use also canopy scale in ddep (l_gannzeori=F)
        rleaf(jl)=cvk*cvc/(cvb*LOG((zabcs*EXP(cvk*zzlai)+1._dp)            &
            /(zabcs+1._dp))/cvabc-LOG((zabcs+EXP(-cvk*zzlai))/(zabcs+1._dp)))  
      
      CASE(3)
        ! canopy scale + adjustment according to Jarvis (1976) - Temperature stress factor
        ! rco_leaf(PAR,LAI) *1/f(T)
        tmp_rleaf=cvk*cvc/(cvb*LOG((zabcs*EXP(cvk*zzlai)+1._dp)                 &
            /(zabcs+1._dp))/cvabc-LOG((zabcs+EXP(-cvk*zzlai))/(zabcs+1._dp)))         
        rleaf(jl) = tmp_rleaf   * 1._dp/fT(ptslm1(jl))


      CASE(4)
        ! canopy scale + adjustment according to Katul (2009) - Drought stress factor
        ! rco_leaf(PAR,LAI) *1/f(D)
        tmp_rleaf=cvk*cvc/(cvb*LOG((zabcs*EXP(cvk*zzlai)+1._dp) &
            /(zabcs+1._dp))/cvabc-LOG((zabcs+EXP(-cvk*zzlai))/(zabcs+1._dp)))
        rleaf(jl) = tmp_rleaf   * 1._dp/fD(ptslm1(jl),rh_2m(jl))

      CASE(5)
        ! all changes - rco_leaf(PAR,LAI) *f(T)*f(D)                        
        tmp_rleaf=cvk*cvc/(cvb*LOG((zabcs*EXP(cvk*zzlai)+1._dp)                 &
            /(zabcs+1._dp))/cvabc-LOG((zabcs+EXP(-cvk*zzlai))/(zabcs+1._dp)))   
        rleaf(jl) = tmp_rleaf * 1._dp/fT(ptslm1(jl)) * 1._dp/fD(ptslm1(jl),rh_2m(jl))
     END SELECT
  
  END DO
  
  END SUBROUTINE vertex_calc_rstom    
!=========================================================================
!=========================================================================
  SUBROUTINE vertex_surftemp(klon, pdt,                                   &
! CONSTANTS
          pemi, pboltz, pcp, pc16, platev, platsu,                        &
! COEFFICIENTS FROM THE ELIMINATION
          pfscoe, pescoe, pfqcoe, peqcoe,                                 &
! OLD VALUES AT THE SURFACE
          psold, pqsold, pdqsold,                                         &
! OTHER FLUXES
          pnetrad, pgrdfl,                                                &
! DIFFUSION COEFFICIENTS, CAIR AND CSAT FOR EVAP AND SOIL HEAT CAPACITY
          pcfh, pcair, pcsat, pfracsu, pgrdcap,                           &
! Logical land mask
          lpland,                                                         &
! OUTPUT
          psnew, pqsnew)
!
! COMPUTES THE ENERGY BALANCE AT THE SURFACE WITH AN IMPLICIT SCHEME
! THAT IS CONNECTED TO THE RICHTMYER AND MORTON ALGORITHM OF THE PBL.
!
! INPUT
! -----
! KLON     : LENGTH OF ARRAYS TO BE USED
! PDT      : LEAP-FROG TIMESTEP IN SECONDS TIMES ALPHA (1.5)
!            ALPHA USED TO COMPENSATE FOR EVALUATING X* INSTEAD OF
!            X^(t+1); X* = ALPHA X^(t+1) + ( 1 - ALPHA ) X^(t-1)
!
! PEMI     : SURFACE EMISSIVITY
! PBOLTZ   : STEFAN-BOLTZMANN CONSTANT
! PCP      : SPECIFIC HEAT OF AIR
! PC16     : CPD*VTMPC2=CPD*(DELTA-1) (CF. VDIFF),FOR SENS.HEAT FL.
! PLATEV   : LATENT HEAT OF EVAPORATION
! PLATSU   : LATENT HEAT OF SUBLIMATION
!
! PFSCOE, PESCOE : COEFFICIENTS OF THE RICHTMYER AND MORTON SCHEME
!                  FOR DRY STATIC ENERGY
! PFQCOE, PEQCOE : AS ABOVE BUT FOR SPECIFIC HUMIDITY
!
! PSOLD   : OLD SURFACE DRY STATIC ENERGY (TS * CP)
! PQSOLD  : SATURATED  SPECIFIC HUMIDITY FOR OLD TEMPERATURE
! PDQSOLD : DERIVATIVE OF SATURATED  SPECIFIC HUMIDITY AT THE
!           OLD TEMPERATURE
!
! PNETRAD : NET RADIATION AT THE SURFACE (UPWARD LONGWAVE IS
!           INCLUDED BUT FOR THE OLD SURFACE TEMPERATURE)
! PGRDFL  : GROUND HEAT FLUX
!
! PCFH    : DIFFUSION COEFFICIENT FOR STATIC ENERGY AND MOISTURE
! PCAIR   : COEFFICIENT IN LATENT HEAT FLUX FORMULA (SEE VDIFF)
! PCSAT   : COEFFICIENT IN LATENT HEAT FLUX FORMULA (SEE VDIFF)
! PFRACSU : FRACTION OF SURFACE FOR SUBLIMATION
! PGRDCAP : SURFACE HEAT CAPACITY
!
! OUTPUT
! ------
! PSNEW   : NEW SURFACE STATIC ENERGY
! PQSNEW  : NEW SATURATED SURFACE AIR MOISTURE
!
!
! AUTHOR.
! -------
!
! J. POLCHER  *LMD*  AND  J.-P. SCHULZ  *MPI*,  MAY 1995
! H.G. OUWERSLOOT *MPIC*, 2016
!
!
! MODIFICATIONS.
! --------------
!
! J.-P. SCHULZ  *MPI*,  OCTOBER 1997:
!    MODIFY ACCORDING TO LATENT HEAT FLUX FORMULATION IN VDIFF
!    USING ZCAIR AND ZCSAT COEFFICIENTS.
!
! J.-P. SCHULZ  *MPI*,  AUGUST 1998:
!    MODIFY ACCORDING TO SENSIBLE HEAT FLUX FORMULATION IN VDIFF.
!
! PF.COE IS ADAPTED TO EVALUATE X INSTEAD OF X/ALPHA IN EQUATIONS
!
! OUWERSLOOT *MPIC*, 2016:
! TRANSFER FROM ECHAM TO MESSY AND CLARIFICATIONS
!
    IMPLICIT NONE

!   ARGUMENTS
    INTEGER :: klon
    REAL(dp):: pdt, pemi, pboltz, pcp(klon), pc16, platev, platsu
    REAL(dp):: pfscoe(klon), pescoe(klon), pfqcoe(klon), peqcoe(klon)
    REAL(dp):: psold(klon), pqsold(klon), pdqsold(klon)
    REAL(dp):: pnetrad(klon), pgrdfl(klon)
    REAL(dp):: pcfh(klon), pcair(klon), pcsat(klon), pfracsu(klon)
    REAL(dp):: pgrdcap(klon)
    REAL(dp):: psnew(klon), pqsnew(klon)
    LOGICAL :: lpland(klon)

    INTEGER :: jl
    REAL(dp):: zcolin, zcohfl, zcoind, zicp, zca, zcs

!   MAIN ROUTINE
    DO jl = 1,klon
      IF (lpland(jl)) THEN
      ! Defining zcolin, zcohfl and zcoind such that at the surface
      ! s^{new} = ( zcolin s^{old} + zcoind) / ( zcolin + zcohfl )
      ! (Cs/cp)*(s^{new}-s^{old})/(alpha*2*Delta{T} =
      !   R_n^{impl.} + LE^{impl.} + H^{impl.} + G + X, (fluxes pointing to the surface)
      ! See Eq. (24) of Schulz et al. (2001); extra term X
      ! X = - ( LE^{impl.} / L ) * (cpv - cpd) T_S^{old}
      ! Compensates the use of dry static energy flux, H, (based on s_S and s_air) instead of
      ! directly using SH itself (based on T_S and T_air):
      !   w's' = w'(cp T)' => rho cp w'T' = rho w's' - rho T w'cp'
      !   w'cp' = w'q' (cpv - cpd)
      !   => SH = H - LE / L * T * (cpv - cpd)

        zicp = 1._dp/pcp(jl)                                                ! 1 / cp

        zca    = platsu*pfracsu(jl) +  platev*(pcair(jl) - pfracsu(jl))     ! L * beta
        zcs    = platsu*pfracsu(jl) +  platev*(pcsat(jl) - pfracsu(jl))     ! L * beta * h

        zcolin = pgrdcap(jl)*zicp +                                       & ! Cs / cp + 2 alpha Delta{t} * ( 1/cp * 4*epsilon*sigma*(s_old/c_p)^3
                          pdt*(zicp*4._dp*pemi*pboltz*                    & !        - rho C_H |U| ( L beta E_q - L beta h
                          ((zicp*psold(jl))**3._dp) -                     & !             - 1/cp * (cpv - cpd) s^{old} ( beta E_q - beta h) )
                          pcfh(jl)*(zca*peqcoe(jl) - zcs -                & !        * 1/cp d{q_s}/d{T} )
                                    zicp*pc16*psold(jl)*                  &
                                    (pcair(jl)*peqcoe(jl) - pcsat(jl)))*  &
                               zicp*pdqsold(jl))

        zcohfl = -pdt*pcfh(jl)*(pescoe(jl)-1._dp)                           ! - 2 alpha Delta{t} rho C_H |U| ( E_s - 1)

        zcoind = pdt * (pnetrad(jl) + pcfh(jl)*pfscoe(jl) +  pcfh(jl)*    & ! 2 alpha Delta{t} * (R_n^{expl.} + rho C_H |U| F_s +
                          ((zca*peqcoe(jl)-zcs)*pqsold(jl) +              & !     rho C_H |U| ( ( L beta E_q - L beta h ) q_s^{old}
                          zca*pfqcoe(jl) - zicp*pc16*psold(jl)*           & !         + L beta F_q - 1/cp * (cpv - cpd) s^{old} *
                         ((pcair(jl)*peqcoe(jl) - pcsat(jl))*pqsold(jl) + & !             ( ( beta E_q - beta h) q_s^{old} + beta F_q ) )
                             pcair(jl)*pfqcoe(jl))) + pgrdfl(jl))           !     + G )

        psnew(jl) = (zcolin * psold(jl) + zcoind) / (zcolin + zcohfl)
        pqsnew(jl) = pqsold(jl) + zicp * pdqsold(jl) * (psnew(jl) -       & ! q_s,new = q_s,old + (1 / cp) * d{q_s}/d{T} * (s_new - s_old)
                                                              psold(jl))
      ELSE
        psnew(jl) = psold(jl)
        pqsnew(jl) = pqsold(jl)
      END IF
    END DO
    RETURN
  END SUBROUTINE vertex_surftemp
!=========================================================================

END MODULE MESSY_VERTEX

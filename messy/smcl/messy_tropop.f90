! *************************************************************************
MODULE messy_tropop
! *************************************************************************

  ! TROPOPAUSE DIAGNOSTICS
  !
  ! Author: Patrick Joeckel, MPICH, Mainz, Sep 2003

  USE messy_main_constants_mem, ONLY: dp

  IMPLICIT NONE
  PRIVATE

  CHARACTER(LEN=*), PARAMETER, PUBLIC :: modstr = 'tropop'
  CHARACTER(LEN=*), PARAMETER, PUBLIC :: modver = '2.1'

  ! a - b*cos^2(latitude) [hPa] -> climatological tropopopause
  REAL(DP), PUBLIC :: r_climtp(2) = (/300._dp, 215._dp/)
  LOGICAL, PUBLIC  :: l_wmo_clim_corr = .false.  ! correct WMO tropopause with
                                                ! climatolog. average ?
  REAL(DP), PUBLIC :: r_lat      = 30._dp ! [deg] interface-|latitude| in PV-WMO
  REAL(DP), PUBLIC :: r_dyntp_PV = 3.5_dp ! [PVU] |PV| at dynamical tropopause
  ! look for PV-isoline in this pressure interval [Pa]
  REAL(DP), PUBLIC :: r_press_range_PV(2) = (/0.0_dp, 110000._dp/)
  LOGICAL, PUBLIC  :: l_o3_PV     = .false.      ! calculate O3(PV) ?
  LOGICAL, PUBLIC  :: l_n2o       = .false.      ! calculate N2O(O3)
  LOGICAL, PUBLIC  :: l_noy       = .false.      ! calculate NOy(N2O)
  LOGICAL, PUBLIC  :: l_pblh      = .false.      ! calculate planet. bnd. layer
  LOGICAL, PUBLIC  :: l_slp       = .false.      ! calculate sea-level pressure

  LOGICAL, PUBLIC  :: l_TI       = .false.       ! turbulence index 
  LOGICAL, PUBLIC  :: l_N2       = .false.       ! Brunt-Vaeisaelae frequency

  LOGICAL, PUBLIC  :: l_cpt       = .false.      ! coldpoint diagnostics

  LOGICAL, PUBLIC  :: l_windspeed =.false.       ! calculate wind speed
  
  INTERFACE calc_PV
     MODULE PROCEDURE calc_PV_1d
     MODULE PROCEDURE calc_PV_2d
  END INTERFACE

  INTERFACE pblindex
     MODULE PROCEDURE pblindex_1d
     MODULE PROCEDURE pblindex_2d
  END INTERFACE

  PUBLIC :: tropop_read_nml_ctrl ! read namelist
  PUBLIC :: calc_PV     ! calculation of potential vorticity (1 column)
  PUBLIC :: if_PVtropop ! level index and trop. box fraction of PV-tropopause
  PUBLIC :: WMOtropop   ! calculate WMO tropopause pressure
  PUBLIC :: pv_to_o3    ! Ozone parameterization using PV
  PUBLIC :: o3_to_n2o   ! N2O parameterization using O3
  PUBLIC :: n2o_to_noy  ! NOy parameterization using N2O

  ! ROUTINES FOR PBLH CALCULATION
  PUBLIC :: pblheight   ! calculation of planetary boundary layer height
  PUBLIC :: pblindex    ! calculation of pblh index of ECHAM5 vertical grid

  PUBLIC :: calc_bulk_richardson
  PUBLIC :: calc_pbl_brn

  PUBLIC :: calc_slp
  PUBLIC :: cptdiag

  PUBLIC :: calc_TI
  PUBLIC :: calc_N2

  PUBLIC :: calc_windspeed

CONTAINS

! ----------------------------------------------------------------------
  SUBROUTINE tropop_read_nml_ctrl(status, iou)

    ! TROPOP MODULE ROUTINE (CORE)
    !
    ! READ NAMELIST, CHECK IT, AND INITIALIZE GLOBAL VARIABLES
    !
    ! Author: Patrick Joeckel, MPICH, Sep 2003

    USE messy_main_tools, ONLY: read_nml_open, read_nml_check, read_nml_close

    IMPLICIT NONE

    ! I/O
    INTEGER, INTENT(OUT) :: status ! error status
    INTEGER, INTENT(IN)  :: iou    ! logical I/O unit

    NAMELIST /CTRL/ r_climtp  &
         , l_wmo_clim_corr    &
         , r_lat            &
         , r_dyntp_PV       &
         , r_press_range_PV &
         , l_o3_PV          &
         , l_n2o            &
         , l_noy            &
         , l_pblh           &
         , l_slp            &    ! um_ak_20120103
         , l_TI             &    ! um_ch_20150127
         , l_N2             &    ! um_ch_20150127
         , l_cpt            &    ! op_sd_20121214
         , l_windspeed

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'tropop_read_nml_ctrl'
    LOGICAL                     :: lex          ! file exists ?
    INTEGER                     :: fstat        ! file status

    ! INITIALIZE
    status = 1 ! ERROR

    ! INITIALIZE GLOBAL CONTROL VARIABLES
    ! -> DEFAULT VALUES ARE SET AT DECLARATION ABOVE

    CALL read_nml_open(lex, substr, iou, 'CTRL', modstr)
    IF (.not.lex) RETURN    ! <modstr>.nml does not exist

    READ(iou, NML=CTRL, IOSTAT=fstat)
    CALL read_nml_check(fstat, substr, iou, 'CTRL', modstr)
    IF (fstat /= 0) RETURN  ! error while reading namelist

    ! DIAGNOSE NAMELIST AND SET GLOBAL SWITCHES
    ! CHECK TRACER INITIALIZATION
    WRITE(*,*) ' ... climatological tropopause at'
    WRITE(*,*) '     ',r_climtp(1),' - ',r_climtp(2),' * cos^2(latitude) hPa'

    IF (l_wmo_clim_corr) THEN
       WRITE(*,*) ' ... WMO gaps filled with climatol. average: ON'
    ELSE
       WRITE(*,*) ' ... WMO gaps filled with climatol. average: OFF'
    END IF

    WRITE(*,*) ' ... dynamical tropopause at PV = ',r_dyntp_PV,' PVU'
    WRITE(*,*) '     in pressure interval ['  &
         , r_press_range_PV(1),',',r_press_range_PV(2),'] [Pa]'

    IF ((r_lat < 0.) .OR. (r_lat > 90.)) THEN
       WRITE(*,*) ' ERROR: PV/WMO-interface-|latitude| (r_lat) out of range!'
       RETURN
    ELSE
       WRITE(*,*) ' ... PV/WMO-interface-|latitude|: ',r_lat,' deg'
    END IF

    IF (l_NOy) THEN
       WRITE(*,*) ' ... calculate NOy(N2O) ...'
       IF (.NOT.l_N2O) THEN
          WRITE(*,*) ' ... -> N2O required ; switching ON'
          l_N2O = .true.
       END IF
    ELSE
       WRITE(*,*) ' ... skip calculation of NOy(N2O) ...'
    END IF

    IF (l_N2O) THEN
       WRITE(*,*) ' ... calculate N2O(O3) ...'
       IF (.NOT.l_o3_PV) THEN
          WRITE(*,*) ' ... -> O3 required ; switching ON'
          l_o3_PV = .true.
       END IF
    ELSE
       WRITE(*,*) ' ... skip calculation of N2O(O3) ...'
    END IF

    IF (l_o3_PV) THEN
       WRITE(*,*) ' ... calculate O3(PV, month) ...'
    ELSE
       WRITE(*,*) ' ... skip calculation of O3(PV, month)'
    END IF

    IF (l_pblh) THEN
       WRITE(*,*) ' ... calculate planetary boundary layer height ...'
    ELSE
       WRITE(*,*) ' ... skip calc. of planetary boundary layer height ...'
    END IF

    IF (l_slp) THEN
       WRITE(*,*) ' ... calculate sea-level pressure ...'
    ELSE
       WRITE(*,*) ' ... skip calc. of  sea-level pressure...'
    END IF

    IF (l_cpt) THEN
       WRITE(*,*) ' ... calculate coldpoint diagnostics ...'
    ELSE
       WRITE(*,*) ' ... skip calculation of coldpoint diagnostics ...'
    END IF

    IF (l_windspeed) THEN
       WRITE(*,*) ' ... calculate 3D windspeed ...'
    ELSE
       WRITE(*,*) ' ... skip calculation of 3d windspeed ...'
    END IF

    CALL read_nml_close(substr, iou, modstr)

    status = 0  ! no ERROR

  END SUBROUTINE tropop_read_nml_ctrl
! ----------------------------------------------------------------------

! -------------------------------------------------------------------------
SUBROUTINE calc_pv_1d(pvort, pmid, temp, vort, corio, rd, cpd, g)

  IMPLICIT NONE

  ! I/O
  REAL(DP), DIMENSION(:), INTENT(OUT) :: pvort  ! potential vorticity [PVU]
  REAL(DP), DIMENSION(:), INTENT(IN)  :: pmid   ! pressure [Pa]
  REAL(DP), DIMENSION(:), INTENT(IN)  :: temp   ! temperature [K]
  REAL(DP), DIMENSION(:), INTENT(IN)  :: vort   ! ??
  REAL(DP),               INTENT(IN)  :: corio  ! coriolis parameter
  REAL(DP),               INTENT(IN), OPTIONAL :: rd  ! gas const. of dry air
                                                  ! [J/K/kg]
  REAL(DP),               INTENT(IN), OPTIONAL :: cpd ! spec. heat of dry air
                                                  ! at constant pressure
                                                  ! [J/K/kg]
  REAL(DP),               INTENT(IN), OPTIONAL :: g ! gravity accel. [m/s2]

  ! LCOAL
  INTEGER                             :: nk
  REAL(DP), DIMENSION(:), ALLOCATABLE :: zth
  REAL(DP)                            :: zrd
  REAL(DP)                            :: zg
  REAL(DP)                            :: zcpd
  INTEGER                             :: jk
  REAL(DP)                            :: zwf, zdtdp

  IF (PRESENT(rd)) THEN
     zrd = rd
  ELSE
     zrd = 287.05_dp ! [J/K/kg]
  END IF

  IF (PRESENT(cpd)) THEN
     zcpd = cpd
  ELSE
     zcpd = 1005.46_dp ! [J/K/kg]
  END IF

  IF (PRESENT(g)) THEN
     zg = g
  ELSE
     zg = 9.80665_dp ! [m/s2]
  END IF

  nk = SIZE(pmid)
  ALLOCATE(zth(nk))

  pvort(:) = 0.0_dp
  zth(:) = temp(:) * (1.e5_dp / pmid(:) ) ** (zrd/zcpd)

  DO jk = 2, nk-2

     zwf  = (pmid(jk)-pmid(jk-1))    &
          / (pmid(jk+1)-pmid(jk-1))
     
     zdtdp = zwf * (zth(jk+1)-zth(jk))     &
          / (pmid(jk+1)-pmid(jk))          &
          + (1.-zwf) * (zth(jk)-zth(jk-1)) &
          / (pmid(jk)-pmid(jk-1))
     
     ! PV calculated in PV units [1e-6 K m-2 kg-1 s-1]
     pvort(jk) = -1.e6_dp * zdtdp * (vort(jk) + corio) * zg

  END DO

  DEALLOCATE(zth)

END SUBROUTINE calc_pv_1d
! -------------------------------------------------------------------------

! -------------------------------------------------------------------------
SUBROUTINE calc_PV_2d(kproma, pvort, pmid, temp, vort, corio, rd, cpd, g)

  IMPLICIT NONE

  ! I/O
  INTEGER,              INTENT(IN)  :: kproma
  REAL(DP), DIMENSION(:,:), INTENT(OUT) :: pvort  ! potential vorticity [PVU]
  REAL(DP), DIMENSION(:,:), INTENT(IN)  :: pmid   ! pressure [Pa]
  REAL(DP), DIMENSION(:,:), INTENT(IN)  :: temp   ! temperature [K]
  REAL(DP), DIMENSION(:,:), INTENT(IN)  :: vort   ! relative vorticity [1/s]
  REAL(DP), DIMENSION(:),   INTENT(IN)  :: corio  ! coriolis parameter [1/s]
  REAL(DP),           INTENT(IN), OPTIONAL :: rd  ! gas const. of dry air
                                                  ! [J/K/kg]
  REAL(DP),           INTENT(IN), OPTIONAL :: cpd ! specific heat of dry air
                                                  ! at constant pressure
                                                  ! [J/K/kg]
  REAL(DP),           INTENT(IN), OPTIONAL :: g   ! gravity acceleration [m/s2]

  ! LCOAL
  INTEGER                             :: nk
  REAL(DP), DIMENSION(:,:), ALLOCATABLE :: zth
  REAL(DP)                            :: zrd
  REAL(DP)                            :: zg
  REAL(DP)                            :: zcpd
  INTEGER                             :: jl, jk
  REAL(DP)                            :: zwf, zdtdp

  IF (PRESENT(rd)) THEN
     zrd = rd
  ELSE
     zrd = 287.05_dp ! [J/K/kg]
  END IF

  IF (PRESENT(cpd)) THEN
     zcpd = cpd
  ELSE
     zcpd = 1005.46_dp ! [J/K/kg]
  END IF

  IF (PRESENT(g)) THEN
     zg = g
  ELSE
     zg = 9.80665_dp ! [m/s2]
  END IF

  nk = SIZE(pmid,2)
  ALLOCATE(zth(kproma,nk))

  pvort(:,:) = 0.0_dp
  zth = temp * (1.e5_dp / pmid ) ** (zrd/zcpd)

  DO jk = 2, nk-2
     DO jl = 1, kproma

        zwf  = (pmid(jl,jk)-pmid(jl,jk-1))    &
             / (pmid(jl,jk+1)-pmid(jl,jk-1))
     
        zdtdp = zwf * (zth(jl,jk+1)-zth(jl,jk))     &
             / (pmid(jl,jk+1)-pmid(jl,jk))          &
             + (1.-zwf) * (zth(jl,jk)-zth(jl,jk-1)) &
             / (pmid(jl,jk)-pmid(jl,jk-1))
     
        ! PV calculated in PV units [1e-6 K m-2 kg-1 s-1]
        pvort(jl,jk) = -1.e6_dp * zdtdp * (vort(jl,jk) + corio(jl)) * zg

     END DO
  END DO

  DEALLOCATE(zth)

END SUBROUTINE calc_PV_2d
! -------------------------------------------------------------------------

! -------------------------------------------------------------------------
SUBROUTINE if_PVtropop(pvort, kproma, klev, k, ft, pv)

  ! index and (troposphere) box fraction of PV tropopause

  IMPLICIT NONE

  ! I/O
  REAL(DP), DIMENSION(:,:), INTENT(IN)  :: pvort  ! potential vorticity [PVU]
  INTEGER,                  INTENT(IN)  :: kproma ! vector length (e.g. nlon)
  INTEGER,                  INTENT(IN)  :: klev   ! number of levels
  INTEGER,  DIMENSION(:),   INTENT(OUT) :: k      ! PV tropopause level
  REAL(DP), DIMENSION(:),   INTENT(OUT), OPTIONAL :: ft ! fraction in trpsph.
  REAL(DP),                 INTENT(IN),  OPTIONAL :: pv ! PV-value at tropopause

  ! LOCAL
  REAL(DP) :: zpv
  INTEGER  :: jk, jp
  REAL(DP) :: zft, dk
  INTEGER  :: zk

  IF (PRESENT(pv)) THEN
     zpv = pv
  ELSE
     zpv = 3.5_dp
  END IF

  k(:) = 0

  ! FIND HIGHEST TROPOSPHERIC GRIDBOX
  DO jp=1, kproma

     DO jk = 2, klev-1
        IF (ABS(pvort(jp, jk)) < zpv) THEN
           k(jp) = jk
           EXIT      ! (EXIT LEVEL LOOP) LOOPS MUST NOT BE CHANGED !!!
        ENDIF
     ENDDO

     ! ADJUST INDEX
     ! CALCULATE FRACTION OF BOX IN TROPOSPHERE
     !
     ! METHOD: LINEAR INTERPOLATION
     !
     zk = k(jp)
     ! THE FOLLOWING CONDITION MUST ALWAYS BE .TRUE.,
     ! SINCE THE FIRST LEVEL WITH PV < PV_TP FROM THE TOP IS SEARCHED
     ! (SEE LOOP ABOVE !)
     !  -> PV(k-1) > PV(k)  FOR PV(k) < PV_TP
     IF ( ABS( (ABS(pvort(jp,zk)) - ABS(pvort(jp, zk-1)) ) ) &
          > TINY(0.0_dp) ) THEN
        zft = (zpv-ABS(pvort(jp, zk-1))) &
             /(ABS(pvort(jp, zk))-ABS(pvort(jp, zk-1)))  ! e [0,1)
     ELSE
        zft = 0.5  ! SHOULD BE NEVER REACHED !!!
     END IF
     ! zft must be e [0,1]; however at low latitudes PV goes to infinity
     ! and therefore the index of the PV-tropopause can go -> 0
     zft = MAX(0.0_dp, zft)
     zft = MIN(1.0_dp, zft)

     dk = INT(zft+0.5_dp)
     ! zft e [0,0.5] -> dk = 0 -> ONE LEVEL ABOVE (-1)
     ! zft e (0.5,1) -> dk = 1 -> KEEP LEVEL
     k(jp) = k(jp) - 1 + dk

     ! CALCULATE FRACTION OF BOX IN TROPOSPHERE
     ! ONE LEVEL ABOVE (dk = 0) -> zft e [0, 0.5]
     !                          -> ft e [0.5, 0]
     !       EXAMPLE: zft = 0   -> TP AT BOX MID -> trop. fract. = 0.5
     !                zft = 0.5 -> TP AT LOWER INTERFACE -> trop. fract. = 0.0
     ! KEEP LEVEL      (dk = 1) -> zft e (0.5, 1)
     !                          -> ft e (1, 0.5)
     !       EXAMPLE: zft = 0.5 -> TP AT UPPER INTERFACE -> trop. fract. = 1.0
     !                zft = 1   -> TP AT BOX MID -> trop. fract. = 0.5

     !                              FOR  dk = 1          FOR  dk=0
     IF (PRESENT(ft)) ft(jp) = (1.5_dp-zft)*REAL(dk,dp) &
                                + (0.5_dp-zft)*REAL(1-dk,dp)

  END DO

END SUBROUTINE if_PVtropop
! -------------------------------------------------------------------------

! -------------------------------------------------------------------------
SUBROUTINE WMOtropop(ptropo, kproma, klev &
     , temp, pmid, p0, g, rgas, cp, phi) !um_ak_20110627 p0 added

!   original: XTTPHWMO
!   programmed by Dieter Nodorp                 V1.0  March 95
!   modified by Thomas Reichler                 V1.2  May 95
!   Rewritten for use in ECHAM  Christine Land  V2.0  Jan 96
!   Modified for use in MATCH   Patrick Joeckel V2.1  May 98
!   Modified for use in ECHAM5  Patrick Joeckel V2.2  Sep 2003
!            - renamed to WMOtropop
!            - conversion to F90     (STEP 1)
!            - automatic choice of top/bottom level for search
!
!     purpose
!     -------
!     xttphwmo (WMOtropop) computes the tropopause height following the
!     definition of the height of the tropopause as postulated
!     by the WMO (1957).
!
!     interface
!     ---------
!     call xttphwmo(phi,t,pmid,ptropo)
!
!          phi             : gaussian latitude in rad
!          t(plond,plev)   : temperature                 --> ptm1
!          pmid(plond,plev): pressure at layer midpoints --> papm1
!          ptropo(plon)    : tropopause pressure
!
!     call WMOtropop(ptropo, kproma, klev &
!                    , temp, pmid, g, rgas, cp, phi)
!
!          OUTPUT:
!            ptropo(:)       : tropopause pressure
!
!          INPUT:
!            kproma          : vector length (e.g. nlon)
!            klev            : number of levels
!            temp(:,:)       : temperature
!            pmid(:,:)       : pressure at layer midpoints
!          INPUT (OPTIONAL):
!            g               : gravity acc. [m/s2]
!            rgas            : gas constant for dry air [J/K/kg]
!            cp              : specific heat of dry air
!                              at constant pressure [J/K/kg]
!            phi(:)          : gaussian latitude in rad
!
!     method
!     ------
!   - listing of errors
!   - p**k interpolation for the temperature Tm and
!     the temperature gradient Gamma of each model layer
!   - p**kappa interpolation for p-wmo
!
!     references
!     ----------
!
!     WMO (1992): International meteorological vocabulary, Genf, 784pp.
!
!   " 1. The first tropopause is defined as the lowest level at which
!        the lapse rate decreases to 2 deg C per kilometer or less,
!        provided also the average lapse rate between this level and
!        all higher levels within 2 kilometers does not exceed 2 deg C
!        per kilometer."
!
!
!     Reichler (1995): Eine globale Klimatologie der Tropopausenhoehe
!                 auf der Basis von ECMWF-Analysen, Diplomarbeit
!                 Universitaet Augsburg, 145pp.
!
!     Reichler, T., M. Dameris, R. Sausen (2003): Determining the
!                 tropopause height from gridded data,
!                 Geophys. Res. L., 30, No. 20, 2042
!
!     Dameris, M., D. Nodorp and R.Sausen, 1995: Correlation between
!                 Tropopause Height Pressure and TOMS-Data for the
!                 EASOE-Winter 1991/1992, Beitr. Phys. Atmosph., 68,
!                 227-232.
!
!     externals
!     ---------
!           none
!

  IMPLICIT NONE

  ! I/O
  REAL(DP), DIMENSION(:),   INTENT(OUT) :: ptropo ! tropopause pressure [Pa]
  INTEGER,                  INTENT(IN)  :: kproma ! vector length (e.g. nlon)
  INTEGER,                  INTENT(IN)  :: klev   ! number of levels
  REAL(DP), DIMENSION(:,:), INTENT(IN)  :: temp   ! temperature [K]
  REAL(DP), DIMENSION(:,:), INTENT(IN)  :: pmid   ! pressure [Pa]
  REAL(DP), DIMENSION(:),   INTENT(IN)  :: p0     ! stand. press. [Pa]
  REAL(DP),               INTENT(IN), OPTIONAL  :: g    ! gravity acc. [m/s2]
  REAL(DP),               INTENT(IN), OPTIONAL  :: rgas ! gas constant
                                                        ! for dry air [J/K/kg]
  REAL(DP),               INTENT(IN), OPTIONAL :: cp ! spec. heat of dry air
                                                     ! at constant pressure
                                                     ! [J/K/kg]
  REAL(DP), DIMENSION(:), INTENT(IN), OPTIONAL :: phi ! gaussian latitude [rad]

  ! LOCAL
  REAL(DP), PARAMETER         :: zgwmo = -0.002     ! WMO tropopause ...
  REAL(DP), PARAMETER         :: zdeltaz=2000.0     ! ... definitions
  REAL(DP)                    :: zg, zrgas, zcp
  INTEGER                     :: iplimb  ! bottom level
  INTEGER                     :: iplimt  ! top level
  REAL(DP), DIMENSION(kproma) :: zplimb  ! pressure at bottom level
  REAL(DP), DIMENSION(kproma) :: zplimt  ! pressure at top level
  INTEGER                     :: jp, jk, jj, kcount
  REAL(DP), DIMENSION(kproma, klev) :: zpmk, zpm, za, zb, ztm, zdtdz
  REAL(DP)                    :: zptph, zp2km, zkappa, zzkap
  REAL(DP)                    :: zfaktor, zasum, zag, zbg, zaquer, zptf
  LOGICAL                     :: ldtdz               ! dt/dz is valid

  ! INIT
  IF (PRESENT(g)) THEN
     zg = g
  ELSE
     zg = 9.80665_dp ! m/s2
  END IF

  IF (PRESENT(rgas)) THEN
     zrgas = rgas
  ELSE
     zrgas = 287.05_dp ! [J/K/kg]
  END IF
  IF (PRESENT(cp)) THEN
     zcp = cp
  ELSE
     zcp = 1005.46_dp  ! [J/K/kg]
  END IF

  zkappa=zrgas/zcp
  zzkap=1._dp/zkappa
  zfaktor=(-1._dp)*zg/zrgas

  ! INIT
  zpmk(:,:)  = 0.0_dp 
  zpm(:,:)   = 0.0_dp 
  za(:,:)    = 0.0_dp 
  zb(:,:)    = 0.0_dp 
  ztm(:,:)   = 0.0_dp 
  zdtdz(:,:) = 0.0_dp 

  DO iplimt=3, klev
     IF ( p0(iplimt) >= 5000._dp ) EXIT
  END DO
  ! now: 3 <= iplimt <= klev 

  DO iplimb=klev-1, 1, -1
     IF ( p0(iplimb) <= 60000._dp ) EXIT
  END DO
  ! now: 1 <= iplimb <= klev-1 

  DO jp=1,kproma
     ptropo(jp) = -999.0_dp      ! undefined tropopause pressure
     zplimb(jp)=pmid(jp,iplimb)  ! pressure at bottom level
     zplimt(jp)=pmid(jp,iplimt)  ! pressure at top level
  END DO

  ! 2.1 compute dt/dz
  ! -----------------
! Note: 2 <= iplimt <= klev AND 1 <=iplimb <= klev-1 (see loops above)
!       Since in the next loop below, jk+1 is used, one more element
!       must be calculated here.      
  DO jk=iplimb+1,iplimt-1,-1          ! loop over levels
     DO jp=1,kproma                ! vector loop
        ! ztm    lineare Interpolation in p**kappa
        ! gamma  dt/dp = a * kappa + pmid(jp,jk)**(kappa-1.)
        zpmk(jp,jk)=0.5_dp* (pmid(jp,jk-1)**zkappa + pmid(jp,jk)**zkappa)
        zpm(jp,jk)=zpmk(jp,jk)**zzkap                   ! p mitte
        za(jp,jk)=(temp(jp,jk-1)-temp(jp,jk))/ &
             (pmid(jp,jk-1)**zkappa-pmid(jp,jk)**zkappa)
        zb(jp,jk) = temp(jp,jk)-(za(jp,jk)*pmid(jp,jk)**zkappa)
        ztm(jp,jk)=za(jp,jk)*zpmk(jp,jk)+zb(jp,jk)      ! T mitte

        zdtdz(jp,jk)=zfaktor*zkappa*za(jp,jk) &
             * zpm(jp,jk)**zkappa/ztm(jp,jk)
     END DO
  END DO

  DO 2000 jp=1,kproma                     ! loop over longitudes
     ! Note: This loop must be limited to iplimt ... iplimb to avoid
     !       array access out of bounds.

     DO 1000 jk=iplimb,iplimt,-1          ! loop over levels

        ! 2.2 First test: valid dt/dz ?
        ! -----------------------------

        IF ( (zdtdz(jp,jk) > zgwmo) .AND.    &    ! dt/dz > -2K/km
             zpm(jp,jk) <= zplimb(jp)) THEN        ! zpm not too low

           ! 2.3 dtdz is valid > berechne p-wmo durch
           ! Interpolation zwischen aktueller und
           ! darunterliegender Schicht
           ! ----------------------------------------

           ! 1.lineare in p^kappa (= Dieters neue Methode)

           zag = (zdtdz(jp,jk)-zdtdz(jp,jk+1))/ &
                &        (zpmk(jp,jk)-zpmk(jp,jk+1))    ! a-gamma
           zbg = zdtdz(jp,jk+1) - zag*zpmk(jp,jk+1)     ! b-gamma

           IF (((zgwmo-zbg)/zag) < 0._dp) THEN
              zptf = 0._dp
           ELSE
              zptf = 1._dp
           END IF

           zptph = zptf*ABS((zgwmo-zbg)/zag)**zzkap

           ldtdz = (zdtdz(jp,jk+1) < zgwmo)

           IF ( .NOT.ldtdz ) THEN
              zptph=zpm(jp,jk)
           END IF

           IF (zptph < zplimt(jp)) THEN
              ptropo(jp)=-999._dp
              GOTO 1500
           END IF

           ! 2.4 2nd test: dt/dz above 2km must not be lower than -2K/km
           ! -----------------------------------------------------------
           zp2km = zptph + zdeltaz*zpm(jp,jk) &
                / ztm(jp,jk)*zfaktor ! p at ptph + 2km
           zasum = 0.0_dp                          ! zdtdz above
           kcount = 0                              ! number of levels above

           ! 2.5 Test until pm < p2km
           ! --------------------------

           DO 265 jj=jk,iplimt-1,-1
              if(zpm(jp,jj) > zptph)  GOTO 265               ! doesn't happen
              if(zpm(jp,jj).lt.zp2km) GOTO 888               ! ptropo valid

              zasum = zasum+zdtdz(jp,jj)
              kcount = kcount+1
              zaquer = zasum/float(kcount)          ! dt/dz mean

              ! 2.6 discard tropo ?
              ! --------------------

              IF (zaquer <= zgwmo) GOTO 1000 ! dt/dz above < 2K/1000
                                             ! discard it

265           CONTINUE                       ! test next level

              zptph=-999._dp
888           CONTINUE
              ! zptph is valid
              ptropo(jp) = zptph
              GOTO 2000

              ! at current level no valid ptph found
              ! continue search at next higher level
           ELSE
              zptph=-999._dp
           END IF
1000       CONTINUE  ! next level
1500       CONTINUE  ! next level

           ! all levels tested, no valid zptph found
           ! zptph=-999.

2000       CONTINUE


  ! FILL OPTIONALLY  WITH CLIMATOLOGICAL TROPOPAUSE, IN CASE
  ! WHERE NO TROPOPAUSE LEVEL WAS FOUND
  IF (PRESENT(phi)) THEN
     DO jp=1, kproma
        IF (ptropo(jp) < 0.0_dp) ptropo(jp) =   &
            r_climtp(1)*100._dp - (r_climtp(2)*100._dp)*cos(phi(jp))*cos(phi(jp))
!             30000.-21500.*cos(phi(jp))*cos(phi(jp))
     END DO
  END IF

END SUBROUTINE WMOtropop
! -------------------------------------------------------------------------

! -------------------------------------------------------------------------
ELEMENTAL SUBROUTINE pv_to_o3(o3, pv, month)

  ! References:
  !
  !  Roelofs GJ, and Lelieveld, J.,
  !    MODEL ANALYSIS OF STRATOSPHERE-TROPOSPHERE EXCHANGE OF OZONE
  !    AND ITS ROLE IN THE TROPOSPHERIC OZONE BUDGET,
  !    in C.S. Zerofs et al. (eds.), Chemistry and Radiation Changes
  !    in the Ozone Layer, 25-43, Kluwer Academic Publishers, 2000.
  !
  !  Roelofs GJ, and Lelieveld J.,
  !    MODEL STUDY OF THE INFLUENCE OF CROSS-TROPOPAUSE
  !    O-3 TRANSPORTS ON TROPOSPHERIC O-3 LEVELS,
  !    Tellus Series B-Chemical & Physical Meteorology.,
  !    49(1):38-55, 1997.

  ! I/O
  REAL(DP),    INTENT(OUT)      :: o3       ! OZONE [mol/mol]
  REAL(DP),    INTENT(IN)       :: pv       ! potential vorticity [PVU]
  INTEGER, INTENT(IN)           :: month    ! month of the year

  ! LOCAL
  REAL(DP), DIMENSION(12), PARAMETER :: slope = &
       (/ 53.75_dp, 73.43_dp, 75.92_dp, 72.82_dp, 71.99_dp, 65.12_dp, &
          47.70_dp, 42.91_dp, 35.76_dp, 35.80_dp, 36.97_dp, 54.02_dp  &
          /) * 1.e-09_dp   ! mol/mol/PVU

  REAL(DP), DIMENSION(12), PARAMETER :: intersect = &
       (/ 73.12_dp, 145.61_dp, 135.64_dp, 108.17_dp, 94.86_dp,  84.38_dp, &
          36.83_dp,  46.48_dp,  34.93_dp,  47.78_dp, 41.71_dp, 114.63_dp  &
          /) * 1.e-09_dp ! mol/mol

  o3 = ( MAX( ABS(pv) * slope(month) - intersect(month), 0.0_dp) )

END SUBROUTINE pv_to_o3
! -------------------------------------------------------------------------

! -------------------------------------------------------------------------
ELEMENTAL SUBROUTINE o3_to_n2o(n2o, o3, lat)

  ! References:
  !
  !  D.M. Murphy and D.W. Fahey, 
  !    AN ESTIMATE OF THE FLUX OF STRATOSPHERIC REACTIVE NITROGEN
  !    AND OZONE INTO THE TROPOSPHERE,
  !    J. Geophys. Res, Vol. 99(D3), p 5325-5332, 1994
  !
  !    Original correlation:
  !
  !       45N ... 62N : [O3] = 6449 - 20.6*[N2O]
  !       20N ... 40N : [O3] = 6691 - 23.5*[N2O]
  !       20S ... 40S : [O3] = 6846 - 22.4*[N2O]
  !       45S ... 62S : [O3] = 5006 - 16.1*[N2O]
  !
  !    Intervals extrapolated:
  !       45N ... 62N  -> 42.5N ... 90N
  !       20N ... 40N  ->  0.0N ... 42.5N
  !       20S ... 40S  ->  0.0S ... 42.5S
  !       45S ... 62S  -> 42.5S ... 90S
  !

  ! I/O
  REAL(DP), INTENT(OUT) :: n2o    ! N2O [mol/mol]
  REAL(DP), INTENT(IN)  :: o3     ! O3  [mol/mol]
  REAL(DP), INTENT(IN)  :: lat    ! latitude [degrees_north]

  ! LOCAL
  INTEGER :: imh, il, is, in

  in = INT(SIGN(-0.5_dp,lat)+1._dp)      ! northern hemisphere lat >= 0
  is = 1-in                        ! southern hemisphere
  
  ! MID OR HIGH LATITUDE: |lat| > 42.5
  imh = INT(SIGN(-0.5_dp,ABS(lat)-42.5_dp)+1._dp) 
  il  = 1-imh                            ! LOW LATITUDE
  
  n2o = ( (6449.0e-09_dp - o3)/20.6_dp ) * REAL(in*imh,dp) &
      + ( (6691.0e-09_dp - o3)/23.5_dp ) * REAL(in*il,dp)  &
      + ( (6846.0e-09_dp - o3)/22.4_dp ) * REAL(is*il,dp)  &
      + ( (5006.0e-09_dp - o3)/16.1_dp ) * REAL(is*imh,dp)

  ! CUT OFF NEGATIVES
  n2o = MAX(n2o, 0.0_dp)

END SUBROUTINE o3_to_n2o
! -------------------------------------------------------------------------

! -------------------------------------------------------------------------
ELEMENTAL SUBROUTINE n2o_to_noy(noy, n2o, lat)

  ! References:
  !
  !  D.M. Murphy and D.W. Fahey, 
  !    AN ESTIMATE OF THE FLUX OF STRATOSPHERIC REACTIVE NITROGEN
  !    AND OZONE INTO THE TROPOSPHERE,
  !    J. Geophys. Res, Vol. 99(D3), p 5325-5332, 1994
  !
  !    Original correlation:
  !
  !       45N ... 62N : [NOy] = 19.80 - 0.0623*[N2O]
  !       20N ... 40N : [NOy] = 20.85 - 0.0715*[N2O]
  !       20S ... 40S : [NOy] = 21.82 - 0.0739*[N2O]
  !       45S ... 62S : [NOy] = 18.80 - 0.0618*[N2O]
  !
  !    Intervals extrapolated:
  !       45N ... 62N  -> 42.5N ... 90N
  !       20N ... 40N  ->  0.0N ... 42.5N
  !       20S ... 40S  ->  0.0S ... 42.5S
  !       45S ... 62S  -> 42.5S ... 90S
  !

  ! I/O
  REAL(DP), INTENT(OUT) :: noy    ! N2O [mol/mol]
  REAL(DP), INTENT(IN)  :: n2o    ! N2O [mol/mol]
  REAL(DP), INTENT(IN)  :: lat    ! latitude [degrees_north]

  ! LOCAL
  INTEGER :: imh, il, is, in

  in = INT(SIGN(-0.5_dp,lat)+1._dp)      ! northern hemisphere lat >= 0
  is = 1-in                        ! southern hemisphere
  
  ! MID OR HIGH LATITUDE: |lat| > 42.5
  imh = INT(SIGN(-0.5_dp,ABS(lat)-42.5_dp)+1._dp)
  il  = 1-imh                            ! LOW LATITUDE

  noy = (19.80e-09_dp - 0.0623_dp*n2o) * REAL(in*imh,dp) &
      + (20.85e-09_dp - 0.0715_dp*n2o) * REAL(in*il,dp)  &
      + (21.82e-09_dp - 0.0739_dp*n2o) * REAL(is*il,dp)  &
      + (18.80e-09_dp - 0.0618_dp*n2o) * REAL(is*imh,dp)

  ! CUT OFF NEGATIVES
  noy = MAX(noy, 0.0_dp)

END SUBROUTINE n2o_to_noy
! -------------------------------------------------------------------------

! -------------------------------------------------------------------------
SUBROUTINE pblheight (pblh, kproma, nlev, th, q, phi   &   ! ,obklen,buoypr
     !, u, v, ths, heat, qflx, ustar, cn ,ch ,cm        &
     , u, v,      heat, qflx, ustar                    &
     , g )                                                 ! OPTIONAL

  ! ATTENTION: obklen [ obukhov length [m]]
  !            buoypr [ buoyancy production [m^2/s^3]] 
  !            are also diagnostic output of this routine.
  !            Write out to channel ?????

  IMPLICIT NONE

  ! ======================================================================
  !
  ! model for boundary layer turbulent mixing .
  ! routine provided by bert holtslag, march 1990.
  ! (knmi, po box 201, 3730 ae de bilt, the netherlands)
  !
  !    made to work in a general circulation model context by byron bovill
  !    modified by jim rosinski to compile efficiently on a vector machine
  !    further modified by j. j. hack to correct counter-gradient formulat
  !
  !    adjustments by bert holtslag (july 91) for countergradient terms,
  !    and temperature excesses. mechanical mixing depth removed.
  !
  !    * updated to EC_PBLHGHT by E. van Meijgaard, KNMI, Januari  '93
  !    * updated to EC_PBLHGHT by E. van Meijgaard, KNMI, December '95
  !      in order to be callable by the ECHAM3
  !      local vertical diffusion scheme EC_VDIFF /EC4_VDIFF
  !    * F77 to F90 for ECHAM5 done by Michael Traub, MPCH, Nov 2003  
  !  
  !    
  !    Modifications by D. Vogelezang, May 1993:
  !    1. 'Parcel temperature' calculation at first model level
  !    2. 'Delta-U method' in  bound. lay. height formulation 
  !    Modifications by D. Vogelezang, Fall 1994:
  !    3. 'Ustar-shear production' term added
  !
    
  !I/O
  REAL(DP), DIMENSION(:),   INTENT(OUT)  :: pblh  ! boundary-layer height [m] 
  INTEGER                 , INTENT(IN)   :: kproma
  INTEGER                 , INTENT(IN)   :: nlev
  REAL(DP), DIMENSION(:,:), INTENT(IN)   :: th    ! potential temperature [k]
  REAL(DP), DIMENSION(:,:), INTENT(IN)   :: q     ! specific humidity [kg/kg]
  REAL(DP), DIMENSION(:,:), INTENT(IN)   :: phi   ! geopotential height [m]
  REAL(DP), DIMENSION(:,:), INTENT(IN)   :: u     ! windspeed x-direction [m/s]
  REAL(DP), DIMENSION(:,:), INTENT(IN)   :: v     ! windspeed y-direction [m/s]
! REAL(DP), DIMENSION(:),   INTENT(IN)   :: ths  ! surface (potential) temp. [k]
  REAL(DP), DIMENSION(:),   INTENT(IN)   :: heat ! surface kin. heat flux [Km/s]
  REAL(DP), DIMENSION(:),   INTENT(IN)   :: qflx ! surf. kin. moist. flux [m/s]
  REAL(DP), DIMENSION(:),   INTENT(IN)   :: ustar ! surface frict. vel. [m/s] 
! REAL(DP), DIMENSION(:),   INTENT(IN)   :: cn   ! sfc neutral echange coeff.
! REAL(DP), DIMENSION(:),   INTENT(IN)   :: ch   ! sfc exchange coeff. for heat
! REAL(DP), DIMENSION(:),   INTENT(IN)   :: cm   ! sfc exchange coeff. ...
                                                 !  ... for momentum
  REAL(DP),                 INTENT(IN), OPTIONAL  :: g ! gravity [m/s2]

  ! LOCAL
  REAL(DP), DIMENSION(kproma) :: jwork
  REAL(DP), DIMENSION(kproma) :: wrino
  REAL(DP), DIMENSION(kproma) :: tav
  REAL(DP), DIMENSION(kproma) :: ustr
  REAL(DP), DIMENSION(kproma) :: vv2
  REAL(DP), DIMENSION(kproma) :: heatv
  REAL(DP), DIMENSION(kproma) :: obklen
  !REAL(DP), DIMENSION(kproma) :: buoypr
  !REAL(DP), DIMENSION(kproma) :: wtseff
  REAL(DP), DIMENSION(kproma) :: tlv
  !
  REAL(DP) :: vvk,fmt,wsc,tkv,zrino,ztherm
  !REAL(DP) :: zbn,zbh
  REAL(DP) :: zrw
  REAL(DP) :: zg
  INTEGER  :: i,k, km
  INTEGER  :: nlevp, npbl
  INTEGER  :: jq
  
  !  pbl constants
  REAL(DP), PARAMETER :: fak   = 8.5_dp
  REAL(DP), PARAMETER :: onet  = 0.33333333_dp

! correction according to Vogelezang and Holtslag (1996), Eq. (3)

  REAL(DP), PARAMETER :: ricr  = 0.25_dp
  REAL(DP), PARAMETER :: vk    = 0.4_dp
  REAL(DP), PARAMETER :: binm  = 1.5_dp
  REAL(DP), PARAMETER :: zcrdq = 0.61_dp
  REAL(DP), PARAMETER :: zacb  = 100._dp 
  REAL(DP), PARAMETER :: tiny  = 1.e-9_dp  ! bound wind**2 to guard divide by 0
    
  ! INIT
  IF (PRESENT(g)) THEN
     zg = g
  ELSE
     zg = 9.80665_dp  ! [m/s2]
  END IF
  nlevp = nlev+1
  npbl  = nlev-1
      
  !  set up local arrays: compute bottom level virtual temperature
  !  and heat flux

  DO  i=1,kproma
     tav(i) = th(i,nlev)*(1._dp + zcrdq*q(i,nlev))
     ustr(i) = max(ustar(i),0.01_dp)
     vv2(i) = u(i,nlev)*u(i,nlev) + v(i,nlev)*v(i,nlev)  &
          + zacb*ustr(i)*ustr(i)
     vv2(i) = max(vv2(i),tiny)
     heatv(i) = heat(i) + zcrdq*th(i,nlev)*qflx(i)
     obklen(i) = -tav(i)*ustr(i)**3/                              &
          (zg*vk*(heatv(i) + sign(1.e-10_dp,heatv(i))))
     !buoypr(i) = zg*heatv(i)/tav(i)
     
     !  interpolate ts and qs from at a fixed height of 10m in the surf
     !  layer using the surface fluxes and the surface and first level
     !  compute the corresponding virtual temperature.
     
     !zbn = vk/sqrt(cn(i))
     !zbh = vk*sqrt(cm(i)) / ch(i)
     !
     !wtseff(i) = tav(i)
     
  END DO
  
  !  calculation of boundary layer height
  
  DO i=1,kproma
     jwork(i) = 1
     vvk = u(i,nlev)*u(i,nlev) + v(i,nlev)*v(i,nlev) &
          + zacb*ustr(i)*ustr(i)
     vvk = max(vvk,tiny)
     wrino(i) = (tav(i)-tav(i))*phi(i,nlev)/(tav(i)*vvk)
     jq= jqif(wrino(i),ricr)
     pblh(i) = jq*phi(i,nlev)/zg+(1-jq)*0._dp
     jwork(i) = (1-jq)
  END DO
  
  DO k=2,npbl
     km = k - 1
     DO i=1,kproma
        vvk = (u(i,nlevp-k)-u(i,nlev))**2 +         &
             (v(i,nlevp-k)-v(i,nlev))**2            &
             + zacb*ustr(i)*ustr(i)
        vvk = max(vvk,tiny)
        tkv = th(i,nlevp-k)*(1. + zcrdq*q(i,nlevp-k))
        zrino = (tkv-tav(i))*(phi(i,nlevp-k)-phi(i,nlev))     &
             /(tav(i)*vvk)
        zrino=jwork(i)*zrino
        !
        ! jq=0 if zrino.ge.ricr
        !
        jq=jqif(ricr,zrino)
        zrw=(zrino-wrino(i)) + sign(1.e-5_dp,zrino-wrino(i))
        pblh(i) = jq*pblh(i)+ (1-jq)*                 &
             (  phi(i,nlevp-km) + (ricr-wrino(i))/    &   
             zrw*(phi(i,nlevp-k)-phi(i,nlevp-km)) )/zg
        !
        jwork(i) = jq*jwork(i)
        
        ! if pblh is found already avoid dividing by zero next
        ! level by a faked value
        ! of wrino(i) (jwork=0 already)
        !
        wrino(i)=zrino+(1-jwork(i))*0.1_dp
     END DO
  END DO
 
  DO i=1,kproma
     pblh(i) = (1-jwork(i))*pblh(i)+jwork(i)*phi(i,nlevp-npbl)/zg
     jq=jqif(heatv(i),0._dp)
     
     !  unstable case; compute first estimate of velocity scale,
     !  and thermal temperature excess
     
     fmt = (jq*(1. - binm*pblh(i)/obklen(i))+(1-jq))**onet
     wsc = ustr(i)*fmt
     !ztherm = (heat(i) + zcrdq*th(i,nlev)*qflx(i))*fak/wsc
     ztherm =  heatv(i)*fak/wsc
     !  improve phbl estimate under convective conditions using
     !  convective temperature excess (ztherm)
     
     !dv   Error found 2/2/93. Wrong calculation of initial wrino(i) in
     !dv   loop 507. wrino(i) now divided by vv2(i) instead of (i) !!!!
     !  
     wrino(i)= jq*(-ztherm*(phi(i,nlev)-phi(i,nlev)) &
          /(tav(i)*vv2(i)))+(1-jq)*0.1_dp
     tlv(i) = jq*( tav(i) + ztherm)
     jwork(i)=jq
  END DO
 
  DO k=2,npbl
     km = k - 1
     
     DO i=1,kproma
        vvk = (u(i,nlevp-k)-u(i,nlev))**2 +    &
             (v(i,nlevp-k)-v(i,nlev))**2       &
             + zacb*ustr(i)*ustr(i)
        vvk=max(vvk,tiny)
        tkv = th(i,nlevp-k)*(1. + zcrdq*q(i,nlevp-k))
        zrino = (tkv-tlv(i))*(phi(i,nlevp-k)-phi(i,nlev))    &
             /(vvk*tav(i))
        zrino=zrino*jwork(i)
        !
        ! jq=0 if zrino.ge.ricr
        !
        jq=jqif(ricr,zrino)
        zrw=(zrino-wrino(i)) + sign(1.e-5_dp,zrino-wrino(i))
        pblh(i) = jq*pblh(i)+(1-jq)*                      &
             (  phi(i,nlevp-km) + (ricr-wrino(i))/        &
             zrw*(phi(i,nlevp-k)-phi(i,nlevp-km)) )/zg
        jwork(i) = jq*jwork(i)
        !
        ! if pblh is found already avoid dividing by zero
        ! next level by a faked value
        ! of wrino(i) (jwork=0 already)
        !
        wrino(i)=zrino+(1-jwork(i))*0.1_dp
     END DO
  END DO
  !
  DO i=1,kproma
     pblh(i) = (1-jwork(i))*pblh(i) + jwork(i)*phi(i,nlevp-npbl)/zg
  END DO
  !+
  ! pbl height is now available in pblh(i)
  !

CONTAINS

  ! -----------------------------------------------------------------------
  INTEGER FUNCTION jqif (val1,val2)
    
    IMPLICIT NONE
    
    ! jqif=1 if xxxz>yyyz else jqif=0
    
    REAL(DP), INTENT(IN) :: val1, val2

    ! ORIGINAL VERSION DOES NOT WORK CORRECTLY ON FORTRAN SYSTEMS
    ! WHICH DISTINGUISH BETWEEN -0 AND +0
    jqif = nint(0.5_dp + sign(0.5_dp , (val1-val2))) 

    RETURN
    
  END FUNCTION jqif
  ! -----------------------------------------------------------------------
   
END SUBROUTINE pblheight
! -------------------------------------------------------------------------

! -------------------------------------------------------------------------
SUBROUTINE pblindex_1d(k, pblh, altitude)

  ! 
  ! CALCULATE THE INDEX OF THE PLANETARY BOUNDARY LAYER HEIGHT 
  !
  ! Author: Michael Traub, MPICH, Nov 2003

  IMPLICIT NONE

  !I/O
  INTEGER,                INTENT(OUT)           :: k      ! index of BL-HEIGHT
  REAL(DP),               INTENT(IN)            :: pblh   ! BL-HEIGHT
  REAL(DP), DIMENSION(:), INTENT(IN)            :: altitude ! height above ground at vertical interfaces
  
  ! LOCAL
  INTEGER  :: nk

  ! INIT
  !
  nk   = SIZE(altitude)

  DO k=nk-1,1,-1
     IF (pblh < altitude(k)) EXIT     
  END DO

END SUBROUTINE pblindex_1d
! -------------------------------------------------------------------------

! -------------------------------------------------------------------------
SUBROUTINE pblindex_2d(kproma, k, pblh, altitude)

  ! as pblindex_1d, however for vectors ... 

  IMPLICIT NONE

  !I/O  
  INTEGER,                  INTENT(IN)            :: kproma
  INTEGER,  DIMENSION(:),   INTENT(OUT)           :: k      ! index of BL-HEIGHT
  REAL(DP), DIMENSION(:),   INTENT(IN)            :: pblh   ! BL-HEIGHT
  REAL(DP), DIMENSION(:,:), INTENT(IN)            :: altitude ! height above ground at vertical interfaces
 
  ! LOCAL
  INTEGER  :: jl, jk, nk
  INTEGER, DIMENSION(kproma) :: ilfound

  ! INIT
  !
  nk   = SIZE(altitude,2)
  ilfound(:) = 0

  DO jk=nk-1,1,-1
     DO jl = 1, kproma
        IF ( ilfound(jl) == 1 ) CYCLE
        IF (pblh(jl) < altitude(jl,jk)) THEN
           k(jl) = jk
           ilfound(jl) = 1
        END IF
     END DO
  END DO

END SUBROUTINE pblindex_2d
! -------------------------------------------------------------------------

! The following two subroutines are adopted from the COSMO module file 
! pp_utilities
!==============================================================================

SUBROUTINE calc_bulk_richardson(idim,jdim, kdim,jrow, brn,te,qve,ue,ve,prs &
                                ,hsurf, prs_surf, t_surf, qv_surf, hhl)
  !----------------------------------------------------------------------------
  ! Description:
  !
  ! The bulk Richardson Number (BRN) is a dimensionless number defined as the 
  ! bulk ratio between the buoyant consumption term and the mechanical production
  ! term in the TKE equation (Stull, 1988).
  !
  ! BRN(z)=g* (z-topo)*(theta_v(z)-theta_v_ground)/(theta_v_mean * (u**2+v**2))
  ! 
  ! theta_v_mean is the virtual potential temperature in the layer comprised
  ! between the ground and the height z.
  !
  ! We here consider u=0 m/s and v=0 m/s at the ground (some authors consider
  ! u and v at 2meters). Some others consider theta_v_ground instead of 
  ! theta_v_mean.
  !
  ! References:
  ! ----------- 
  ! - Stull R.,B., 1988: An Introduction to Boundary Layer Meteorology. 
  !   Kluwer Academic Publishers. ISBN 90-277-2768-6
  !
  ! - Sorensen J.H., Rasmussen A., Svensmark H.,1996: Forecast of Atmospheric 
  !   Boundary-Layer Height Utilised for ETEX Real-Time Dispersion Modelling. 
  !   Phys.Chem. Earth, Vol. 21, No. 5-6, pp. 435-439, 1997 Elsevier Science Ltd.
  !
  ! - Jacobson M.Z., 1999: Fundamentals of atmospheric modeling.
  !   Cambridge University Press. ISBN 0-521-63143-2
  !----------------------------------------------------------------------------

  USE messy_main_constants_mem,  ONLY: cp_air, rd, vtmpc1, g

  IMPLICIT NONE

  ! Input data
  !----------- 
  INTEGER, INTENT (IN) ::  &
       jrow,                   & ! current 2nd hor. dim index
       idim ,                  & ! array dimension in "zonal" direction
       jdim ,                  & ! array dimension in "meridional" direction
       kdim                      ! array dimension in vertical direction
  
  REAL(dp),    INTENT (IN) ::  &
       te  (idim,kdim),   & ! environment temperature
       qve (idim,kdim),   & ! environment humidity
       ue  (idim,jdim,kdim),   & ! environment zonal wind 
       ve  (idim,jdim,kdim),   & ! environment meridional wind
       prs (idim,kdim),   & ! full level pressure
       hsurf(idim),       & ! topography
       prs_surf(idim) ,   & ! surface pressure
       t_surf(idim)   ,   & ! surface temperature
       qv_surf(idim)  ,   & ! specific water vapor content on the surface
                            ! (kg/kg)
       hhl (idim,kdim+1)    ! height of half levels
    
  ! Output data
  !------------ 
  REAL (dp), INTENT (OUT) :: brn(idim,kdim) ! Bulk Richardson Number [-]
  
  
  ! Local scalars and automatic arrays
  !-----------------------------------
  REAL (dp) ::       &      
    u_avg(idim,kdim),   & ! averaged zonal wind at the center of the cells
    v_avg(idim,kdim),   & ! aver. meridional wind at the center of the cells
    hfl(idim,kdim),     & ! height of full levels
    theta_v(kdim),           & ! virtual potential temperature on each level
    pasl,                    & ! sea level pressure (ref)
    theta_s,                 & ! potential temperature at the surface 
    theta_v_s,               & ! virtual potential temperature at the surface 
    theta,                   & ! potential temperature on each level
    cumul,                   & ! cumul of virtual potential temperature
    theta_v_cum,             & ! cumulated virtual potential temperature
    theta_v_avg                ! mean virtual potential temperature between 
                               ! k and kdim
  
  INTEGER :: i, k              ! Indices of input/output fields
  
  !-----------------------------------------------------------------------------
  
  ! Reference pressure for potential temperature calculation
  pasl = 100000._dp    
  
  ! Compute height on full LM COARSE levels
  DO k = 1, kdim
     hfl(:,k) = 0.5_dp * hhl(:,k+1) + 0.5_dp * hhl(:,k)
  ENDDO
  
  
  !-------------------------------------------------------------------------
  ! Compute averaged wind velocities at the center of the cells
  !-------------------------------------------------------------------------
  ! zonal wind
  
  DO k=1, kdim
     DO i=1, idim-1
        u_avg(i,k) = (ue(i,jrow,k) + ue(i+1,jrow,k)) / 2       
     ENDDO
     u_avg(idim,k) = u_avg(idim-1,k)
  ENDDO
  
  ! meridional wind
  IF (jrow /= jdim) THEN
     DO k=1, kdim
        DO i=1, idim
           v_avg(i,k) = (ve(i,jrow,k) + ve(i,jrow+1,k)) / 2             
        ENDDO
     ENDDO
  ELSE
     v_avg(1:idim,1:kdim) =  &
          (ve(1:idim,jdim-1,1:kdim) + ve(1:idim,jdim,1:kdim)) / 2  
  ENDIF
  
  !Grid points loop
  
     DO i=1, idim  !x loop
        
        !-----------------------------------------------------------------------
        ! Calculation of all components of the bulk richardson number and 
        ! of the bulk Richardson number on each full level
        !-----------------------------------------------------------------------
        
        ! Potential temperature at the surface 
        theta_s = t_surf(i) * (pasl / prs_surf(i))**(rd/cp_air)
        
        ! Virtual potential temperature at the surface           
        theta_v_s = theta_s * (1. + vtmpc1 * qv_surf(i))
        
        ! Initial cumulated virtual potential temperature
        cumul = 0.0_dp
        
        DO k = kdim, 1, -1  ! Level loop 
           
           ! Virtual potential temperature at each level
           theta = te(i,k) * (pasl / prs(i,k))**(rd/cp_air)
           theta_v (k)= theta * (1 + vtmpc1 * qve(i,k)) 
           
           ! Mean virtual potential temperature between k and kdim
           theta_v_cum  = cumul + theta_v(k)
           theta_v_avg  = theta_v_cum / (kdim - k + 1)
           
           !Bulk Richardson number
           !**********************
           
           ! Limitation of the wind: 1cm/s minimum for each component
           
           IF((u_avg(i,k)**2 + v_avg(i,k)**2) .LE. 2.0E-4_dp)THEN
              brn(i,k) = g *(hfl(i,k)-hsurf(i))*(theta_v(k)-theta_v_s)/ &
                   (theta_v_avg*2.0E-4_dp)
              
           ELSE
              brn(i,k) = g *(hfl(i,k)-hsurf(i))*(theta_v(k)-theta_v_s)/ &
                   (theta_v_avg*(u_avg(i,k)**2+v_avg(i,k)**2))
           ENDIF
           
           
           !Update cumulated virtual potential temperature
           cumul = theta_v_cum
        ENDDO ! Level loop 
     ENDDO ! x
  
END  SUBROUTINE calc_bulk_richardson
!===============================================================================

SUBROUTINE calc_pbl_brn(idim, kdim, te, qve, prs, hhl, hsurf, brn, hpbl, k_pbl)

  !-----------------------------------------------------------------------------
  ! Description:
  !
  ! The PBL height is the estimated planetary boundary layer height in meters 
  ! above the ground. There are several methods to estimate it. In general they
  ! perform well for convective situations. The stable situations are more 
  ! problematic.
  ! 
  ! We choose the bulk richardson method, which is a widely used and robust 
  ! method. 
  ! The height of the PBL is given by the height at which the bulk Richardson
  ! number reaches a prescribed value, the critical Richardson number.
  ! In the literature, one finds critical values between say 0.2 and 1.
  ! We here consider a value of 0.33 under stable conditions (Wetzel, 1982)and
  ! of 0.22 under convective conditions (Vogelezang and Holtslag, 1996).
  !
  ! References:
  ! -----------
  ! - Stull R.,B., 1988: An Introduction to Boundary Layer Meteorology.
  !   Kluwer Academic Publishers. ISBN 90-277-2768-6
  !
  ! - Sorensen J.H., Rasmussen A., Svensmark H.,1996: Forecast of Atmospheric 
  !   Boundary-Layer Height Utilised for ETEX Real-Time Dispersion Modelling. 
  !   Phys.Chem. Earth, Vol. 21, No. 5-6, pp. 435-439, 1997 Elsevier Science Ltd.
  !
  ! - Seibert P. and all, 1999: Review and intercomparison of operational methods
  !   for the determination of the mixing height.
  !   Atmospheric Environment 34 (2000) 1001-1027, 2000 Elsevier Science Ltd.  
  !
  ! - Vogelezang D.H.P., and Holtslag A.A.M., 1996: Evaluation and model impacts
  !   of alternative boundary-layer height formulations.
  !   Boundary-Layer Meteorol. 81, 245-269.
  !
  ! - Wetzel P.J., 1982: Toward parametrization of the stable boundary layer.
  !   J. Appl. Meteorol. 21, 7-13.
  !-----------------------------------------------------------------------------

  USE messy_main_constants_mem,  ONLY: cp_air, rd, vtmpc1

  IMPLICIT NONE

  ! Input data
  !----------- 
  INTEGER, INTENT (IN) ::  &
       idim ,                  & ! array dimension in zonal direction
  !     jdim ,                  & ! array dimension meridional direction
       kdim                      ! array dimension in vertical direction
  
  REAL(dp),    INTENT (IN) ::  &
       te  (idim,kdim),   & ! environment temperature
       qve (idim,kdim),   & ! environment humidity
       prs (idim,kdim),   & ! full level pressure
       hhl (idim,kdim+1), & ! height of half levels
       hsurf(idim),       & ! topo
       brn (idim,kdim)      ! Bulk Richardson Number [-]

  ! Output data
  !------------ 
  REAL(dp), INTENT (OUT) :: hpbl(idim)   ! PBL height [m above ground]
  ! k-index corresponding to the top of the PBL
  REAL(dp), INTENT(OUT)  :: k_pbl(idim)                       
  
  
  ! Local scalars and automatic arrays
  !-----------------------------------
  REAL(dp)               ::  &      
       hfl(idim,kdim),            & ! height of full levels
       zh(idim),zh2(idim),        & ! working arrays for the linear regression
                                    ! zh=sum(hfl); zh2=sum(hfl**2)
       zttav(idim),zhttav(idim),  & ! working arrays for the linear regression
                                    ! zttav=sum(theta_v); zhttav=sum(hfl*theta_v)
       beta(idim),                & ! coefficient of linear regression
       theta_v(kdim),             & ! virtual potential temperature on 4 levels
       theta,                     & ! potential temperature on each level
       pasl,                      & ! sea level pressure (ref) 
       brn_cr                       ! bulk richardson critical value
  
  
  INTEGER :: i, k, kpbl(idim)      ! working indexes
  
  !------------------------------------------------------------------------
  ! Initialisations
  !------------------------------------------------------------------------
  ! Sums for linear regression
  zh(:)    = 0.0_dp
  zttav(:) = 0.0_dp
  zh2(:)   = 0.0_dp
  zhttav(:)= 0.0_dp
  ! Index of the top of the PBL
  kpbl(:)  = -1
  ! ref pressure
  pasl     = 100000.0_dp
  
  ! Compute height on full LM COARSE levels
  DO k = 1, kdim
     hfl(:,k) = 0.5_dp * hhl(:,k+1) + 0.5_dp * hhl(:,k)
  ENDDO
  
  
  ! Determine if stable or unstable (neutral) conditions
  !------------------------------------------------------------------------
  ! calculation of the linear regression coefficient in order to approximate
  ! the temperature gradient in the first 4 layers and evaluate the current
  ! conditions.
  ! if delta_theta_v/delta_z > 0: stable (beta >0)
  ! if delta_theta_v/delta_z < 0: unstable (or neutral) (beta <=0)
  
     DO i=1,idim
        
        DO k=kdim-3,kdim
           ! potential temperature
           theta = te(i,k) * (pasl / prs(i,k))**(rd/cp_air) 
           ! virtual potentiel temperatur
           theta_v (k)= theta * (1 + vtmpc1 * qve(i,k))     
           
           ! sum of the heights        
           zh (i) = zh (i) + hfl(i,k)                   
           ! sum of the (heights)**2   
           zh2(i) = zh2(i) + hfl(i,k)* hfl(i,k)       
           ! sum of the theta_v
           zttav (i) = zttav (i) + theta_v(k)             
           ! sum of the (height*theta_v)
           zhttav(i) = zhttav(i) + hfl(i,k)* theta_v(k) 
        ENDDO
        
        ! coefficient of linear regression
        beta (i)=(zhttav(i)-zh(i)*zttav(i)/4.0_dp)/(zh2(i)-zh(i)*zh(i)/4.0_dp)  
        
        ! set the critical bulk richardson number according to the current
        ! conditions 
        
        IF(beta (i) .GT. 0.0_dp) THEN     ! Stable conditions 
           brn_cr = 0.33_dp                        ! Wetzel value    
        ELSE                                    ! Unstable or neutral conditions
           brn_cr = 0.22_dp                        ! Vogelezang value   
        ENDIF
        
        !----------------------------------------------------------------------
        ! Find k index corresponding to the upper limit of the BL and compute
        ! PBL height
        ! ---------------------------------------------------------------------
        
        ! Scan the levels from the top to the ground and find
        ! the first level where BRN < brn_cr
        
        DO k= 1, kdim
           IF(brn(i,k) .GE. brn_cr) THEN
              kpbl(i) = k
           ENDIF
        ENDDO
        
        ! Compute the corresponding height. If no level found
        ! where BRN < brn_cr, assign a missing_value for the PBL
        ! height.
        
        IF (kpbl(i).LE.0) THEN                     ! no level found
           hpbl(i) = 1.e-37
        ELSE                                       ! algorithm ok
           hpbl(i) =  hfl(i,kpbl(i))-hsurf(i)
        ENDIF
        
     ENDDO !x loop

     k_pbl(:) = REAL(kpbl(:),dp)
  
END  SUBROUTINE calc_pbl_brn

!==============================================================================

!==============================================================================
SUBROUTINE calc_slp (pmsl, ps, t, rho0, dp0, hsurf, ie, g, r_d )
! The subroutine was adopted from the subroutine pp_utilities of the COSMO model
!------------------------------------------------------------------------------
!
! Description:
!   This routine computes the mean sea-level pressure by extrapolating the
!   surface pressure at height z = hsurf hydrostatically to z = 0 whereever
!   the surface height is greater than 0.01 m. For hsurf < 0.01 m, the
!   mean sea-level pressure is simply set to ps.
!   
!   The fields are passed twodimensional, that means the calling procedure has
!   to choose the proper level and time level.
!
! Method:
!   Formula to be described in the documentation
!
!------------------------------------------------------------------------------
  IMPLICIT NONE
  ! Parameter list:
  
  INTEGER, INTENT (IN)    ::    &
       ie              ! horizontal dimensions of the fields
  
  REAL(dp), INTENT (IN)          ::    &
       ps   (ie)    , & ! suface pressure      
       t    (ie)    , & ! temperature at lowest model level
       dp0  (ie)    , & ! full level pressure thickness of lowest layer
       rho0 (ie)    , & ! reference density at lowest full model level
       hsurf(ie)        ! geometrical height of surface

  REAL(dp), INTENT (OUT)         ::    &
       pmsl (ie)        ! mean sea-level pressure
  
  REAL(dp), INTENT (IN)          ::    &
       g,                & ! gravity acceleration
       r_d                 ! gas constant for dry air
  
  ! Local Variables
  
  REAL(dp)         ::    &
       zlapse,           & ! lapse rate ( 6.5K / km )
       zrg,              & ! 1/g                      
       zr3,              & ! 1/3                      
       zdz,              & ! height difference: lowest level - surface 
       ztstar,           & ! temperature at surface (extrapolated from lowest level)
       ztmsl,            & ! temperature at mean sea level (extrapolated from tstar)
       zalph,            & ! zlapse*r_d/g
       zprt,             & ! 
       zprtal              !
  
  INTEGER          ::    &
       i                 ! Loop indicees
  
!------------------------------------------------------------------------------

! Begin subroutine calpmsl

! Initializations
  zlapse = 0.0065_dp
  zrg    = 1.0_dp/g
  zr3    = 1.0_dp/3.0_dp

! Set mean sea-level pressure to ps
  pmsl (:) = ps (:)

  
  DO i = 1, ie
     
     ! Calculate surface temperature "tstar" and mean sea-level 
     ! temperature "tmsl" and correct lapse rate "zalph"
     zdz    = 0.5_dp*dp0(i)/ ( g*rho0(i) )
     ztstar = t(i) + zlapse*zdz
     IF( ztstar < 255.0_dp) THEN
        ztstar = 0.5*( 255.0_dp + ztstar )
     ENDIF
     ztmsl = ztstar + zlapse*hsurf(i)
     zalph = r_d*zlapse*zrg
     IF( ztstar > 290.5_dp ) THEN
        ztstar = 0.5*( 290.5_dp + ztstar )
        ztmsl  = ztstar
        zalph  = 0.0_dp
     ENDIF
     IF( ztstar <= 290.5_dp .AND. ztmsl > 290.5_dp) THEN
        ztmsl  = 290.5_dp 
        zalph  = r_d*(ztmsl-ztstar)*zrg/MAX(hsurf(i), 0.01_dp)
     ENDIF
     
     ! Calculate mean sea level pressure for elevated terrain
     zprt   = g*hsurf(i)/(r_d*ztstar)
     zprtal = zprt*zalph 
     pmsl(i) = ps(i)*EXP &
          ( zprt*(1.0_dp - zprtal*(0.5_dp - zr3*zprtal)) )
  ENDDO
   
END SUBROUTINE calc_slp
!==============================================================================

!==============================================================================
SUBROUTINE cptdiag(temp, p, tpres, PCPT, JLEVCPT)
!

  IMPLICIT NONE

  !I/O
  REAL(DP), DIMENSION(:), INTENT(IN)  :: temp     ! temperature [K]
  REAL(DP), DIMENSION(:), INTENT(IN)  :: p        ! pressure [Pa]
  REAL(DP),               INTENT(IN)  :: tpres    ! tropopause pressure [Pa]
  REAL(DP),               INTENT(OUT) :: PCPT     ! press. coldpoint [Pa]
  INTEGER,                INTENT(OUT) :: JLEVCPT  ! index coldpoint

  ! LOCAL
  INTEGER :: nlev
  INTEGER :: jk, JLEVST

  nlev = SIZE(temp)

  JLEVCPT = nlev
  JLEVST=0
  DO jk=1, nlev
     IF ( (temp(jk) < temp(JLEVCPT)) .AND. &
          (JLEVST == 0) .AND. &
          p(jk) > 5000._dp ) THEN
        JLEVCPT = jk
     ENDIF

     IF ( (p(jk) > tpres) .AND. (JLEVST == 0) ) THEN
        JLEVST = jk
     ENDIF
  ENDDO

  ! pressure heights at cold-point
  PCPT = p(JLEVCPT)

END SUBROUTINE cptdiag
!==============================================================================

!==============================================================================
SUBROUTINE calc_TI (idim,jdim, kdim,jrow, TI, U, V, altitude, dlon, dlat, crlat)
!Calculation of turbulence index (Ellrod-Knapp)
!------------------------------------------------------------------------------

  USE messy_main_constants_mem, ONLY: r_earth => radius_earth, pi

  IMPLICIT NONE
  ! Input data
  !----------- 
  INTEGER, INTENT (IN) ::  &
       jrow , &
       idim ,                  & ! array dimension in "zonal" direction
       jdim ,                  & ! array dimension in "meridional" direction
       kdim                      ! array dimension in vertical direction
       
  REAL(dp), INTENT (IN)          ::    &
       U(idim,jdim,kdim)       , & 
       V(idim,jdim,kdim)       , & 
       altitude(idim,jdim,kdim), &
       dlat ,                    &
       dlon ,                    &
       crlat

  REAL(dp), INTENT (OUT)         ::  TI(idim,kdim) ! turbulent index
     
  ! Local scalars and automatic arrays
  !-----------------------------------
  REAL (dp) ::       &   
       dudx(idim,kdim) ,&
       dvdy(idim,kdim) ,&
       dudy(idim,kdim) ,&
       dvdx(idim,kdim), &
       VWS(idim,kdim)
  REAL(dp) :: dx, dy    
  INTEGER          ::    &
       i,k                 ! Loop indicees

!------------------------------------------------------------------------------
! Begin subroutine calc_TI
  dx = r_earth * (pi/180.0_dp) * dlon * crlat
  dy = r_earth * (pi/180.0_dp) * dlat

  TI   = 0._dp
  dudx = 0._dp
  dvdy = 0._dp
  dudy = 0._dp
  dvdx = 0._dp
  VWS  = 0._dp
 
  IF (jrow /= 1  .and. jrow /= jdim) THEN 
     DO k = 2, kdim-1  
        DO i = 2, idim-1
           dudx(i,k) = (U(i,jrow,k)- U(i-1,jrow,k))/dx
           dvdy(i,k) = (V(i,jrow,k)- V(i,jrow-1,k))/dy
           dudy(i,k) = ((U(i-1,jrow+1,k)+ U(i,jrow+1,k))*0.5_dp   &
                - (U(i-1,jrow-1,k)+ U(i,jrow-1,k))*0.5_dp)/(2._dp*dy)
           dvdx(i,k) = ((V(i+1,jrow,k)+ V(i+1,jrow-1,k))*0.5_dp   &
                - (V(i-1,jrow,k)+ V(i-1,jrow-1,k))*0.5_dp)/(2.0_dp*dy)
           
           VWS(i,k) = SQRT(((((U(i,jrow,k-1)-U(i,jrow,k+1))    &
                +(U(i-1,jrow,k-1)-U(i-1,jrow,k+1)))*0.5_dp)**2)   &
                +(((V(i,jrow,k-1)-V(i,jrow,k+1))               &
                +(V(i,jrow-1,k-1)-V(i,jrow-1,k+1)))*0.5_dp)**2)   &
                /((altitude(i,jrow,k-1)-altitude(i,jrow,k+1)))
           
           TI(i,k) = SQRT((dudx(i,k)-dvdy(i,k))**2                 &
                +(dvdx(i,k)+dudy(i,k))**2)*VWS(i,k)*REAL(10**7,dp)
           
        END DO
     END DO
  END IF
 
END SUBROUTINE calc_TI
!==============================================================================

!==============================================================================
SUBROUTINE calc_N2 (idim,kdim,N2,tpot,altitude)
!Calculation of the Brunt-Vaeisaelae frequency
!------------------------------------------------------------------------------

  USE messy_main_constants_mem, ONLY: g 

  IMPLICIT NONE
  ! Input data
  !----------- 
  INTEGER,  INTENT (IN) ::  &
       idim ,               & ! array dimension in "zonal" direction
       kdim                   ! array dimension in vertical direction
       
  REAL(dp), INTENT (IN) ::  &
       tpot(idim,kdim)    , & 
       altitude(idim,kdim)  

  REAL(dp), INTENT (OUT) ::  N2(idim,kdim) ! turbulent index
  
  ! Local scalars and automatic arrays
  !-----------------------------------
  INTEGER          ::    &
       i,k                 ! Loop indicees
!------------------------------------------------------------------------------
! Begin subroutine calc_N2
  N2      = 0._dp
    
  DO i = 1, idim
     DO k = 2, kdim-1  
        N2(i,k)=(g/tpot(i,k))*((tpot(i,k-1)-tpot(i,k+1))/   &
             ((altitude(i,k-1)-altitude(i,k+1))))
     END DO
  END DO
  
END SUBROUTINE calc_N2
!==============================================================================

elemental real(dp) function calc_windspeed(u,v)

  REAL(dp), INTENT(IN) :: u,v

  REAL(dp) :: speed2

  speed2 = u*u + v*v

  calc_windspeed = SQRT(speed2)

END function calc_windspeed
! *************************************************************************

END MODULE messy_tropop
! *************************************************************************

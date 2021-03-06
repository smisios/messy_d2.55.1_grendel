#include "messy_main_ppd_bi.inc"

!
#ifdef ECHAM5
#define _TFCORR
#endif
! Note: Time Filter Correction (_TFCORR) is only implemented for
!       "classic emssion pints", since emission time series do not
!       have well-defined start and end times.
!       Use with caution!
!
! **********************************************************************
MODULE messy_trexp_si
! **********************************************************************

  ! Tracer Release EXPeriment
  !
  ! INTERFACE FOR ECHAM5 (MESSy/SMIL)
  !
  ! Author: Patrick Joeckel, MPICH, December 2004
  !
  ! References:
  !
  ! Jöckel, P., Kerkweg, A., Pozzer, A., Sander, R., Tost, H., Riede, H.,
  ! Baumgaertner, A., Gromov, S., & Kern, B.: Development cycle 2 of
  ! the Modular Earth Submodel System (MESSy2), Geoscientific Model
  ! Development, 3, 717–752, doi: 10.5194/gmd-3-717-2010,
  ! URL http://www.geosci-model-dev.net/3/717/2010/ (2010)
  !
  ! Patrick Joekel, Carl A. M. Brenninkmeijer, Hanwant B. Singh, and Paul J.
  ! Crutzen, Investigation of the Global Atmospheric Oxidation Efficiency and
  ! Its Trends: A Proposal to Initiate IGAC-GHOST (Global HO Systematic Tests),
  ! IGACtivities Newsletter, Issue No. 28, 2-5, May 2003.
  ! http://www.mpch-mainz.mpg.de/~joeckel/docs/IGAC_Newsletter_May03.pdf
  !
  ! Patrick Joekel, Carl A. M. Brenninkmeijer, and Paul J. Crutzen,
  ! A discussion on the determination of atmospheric OH and its trends,
  ! Atmos. Chem. Phys., 3, 107-118, 2003.
  ! http://www.copernicus.org/EGU/acp/acp/3/107/acp-3-107.pdf
  !
  ! 2011/05/25: Volker Grewe, DLR:
  !    Oscillation due to leap-frog time filtering
  !    avoided for first emission timestep.
  !    Tested for Lagrangian approach only, not finalized
  !    for grid tracers.
  ! 2011/06/29: Patrick Joeckel, DLR:
  !    New feature to use time series emissione (via import_ts)
  ! 2012/07/24: Sabine Brinkop, DLR:
  !    Correction of index calculation for non-equidistant grid
  !    nn_index -> nl_index
  !

  USE messy_main_blather_bi,    ONLY: start_message_bi, end_message_bi
#ifdef MESSYTENDENCY
  USE messy_main_tendency_bi,   ONLY: mtend_get_handle,       &
                                      mtend_get_start_l,      &
                                      mtend_add_l,            &
                                      mtend_register,           &
                                      mtend_id_tracer,        &
                                      mtend_id_t, mtend_id_q
#endif
  USE messy_main_constants_mem, ONLY: STRLEN_MEDIUM, STRLEN_LONG
  USE messy_main_channel,       ONLY: STRLEN_CHANNEL, STRLEN_OBJECT
  USE messy_main_tools,         ONLY: PTR_3D_ARRAY
  USE messy_main_tracer,        ONLY: STRLEN_TRSET
  USE messy_trexp

  IMPLICIT NONE
  PRIVATE
  SAVE

  INTRINSIC :: NULL

  ! NAMELIST LIMITS
  INTEGER, PARAMETER :: NMAXTRAC  = 300 ! NUMBER OF TRACERS
  INTEGER, PARAMETER :: NMAXPOINT = 300 ! NUMBER OF RELEASE POINTS
  INTEGER, PARAMETER :: MAXSTRLEN = 200 ! max. string length

  ! HEIGHT UNITS
  INTEGER, PARAMETER :: HU_undef = 0 ! undefined
  INTEGER, PARAMETER :: HU_hPa   = 1 ! hPa
  INTEGER, PARAMETER :: HU_Pa    = 2 ! Pa
  INTEGER, PARAMETER :: HU_mASL  = 3 ! metre above sea level
  INTEGER, PARAMETER :: HU_mAGL  = 4 ! metre above ground level
  CHARACTER(LEN=5), PARAMETER :: HU_STRING(0:4) = (/ 'undef', 'hPa  ', 'Pa   ', 'm ASL', 'm AGL' /)

  ! EMISSION UNITS
  INTEGER, PARAMETER :: EU_undef  = 0 ! undefined
  INTEGER, PARAMETER :: EU_kg     = 1 ! kg (total mass during emission interval)
  INTEGER, PARAMETER :: EU_kgs    = 2 ! kg/s
  INTEGER, PARAMETER :: EU_mols   = 3 ! mol/s
  INTEGER, PARAMETER :: EU_molecs = 4 ! molecules/s
  CHARACTER(LEN=11), PARAMETER :: EU_STRING(0:4) = (/ 'undef      ', 'kg         ', 'kg/s       ' &
     &                                              , 'mol/s      ', 'molecules/s' /)

  ! USER INTERFACE
  ! ... define tracers, decay and / or reaction partners ...
  TYPE T_RELEASED_TRACER_IO
     CHARACTER(LEN=10*STRLEN_TRSET+10) :: trset  = '' ! TRACER SET NAME(S)
     CHARACTER(LEN=2*STRLEN_MEDIUM+1)  :: trname = '' ! TRACER NAME[_SUBNAME]
     INTEGER                       :: order   = 0      ! ORDER OF REACTION
     REAL(DP)                      :: ka      = 0.0_DP ! Arrhenius-A OR
     !                                                 ! Decay Const.
     REAL(DP)                      :: Ta      = 0.0_DP ! Activation Temperature
     CHARACTER(LEN=STRLEN_CHANNEL) :: channel = '' ! REACTION PARTNER (channel)
     CHARACTER(LEN=STRLEN_OBJECT)  :: object  = '' ! REACTION PARTNER (object)
  END TYPE T_RELEASED_TRACER_IO

  ! ... release points in space and time
  TYPE T_RELEASE_POINT_IO
     !
     INTEGER  :: type   =    0       ! 0: off
     !                               ! 1: classic,
     !                               ! 2: external time series
     !                               ! 3: as 1, but create channel objects with
     !                                          (sum of) emissions
     !                               ! 4: as 2, but create channel objects with
     !                                          (sum of) emissions
     !
     ! shared
     REAL(DP) :: lon    =    0.0_DP  ! degrees east
     REAL(DP) :: lat    =    0.0_DP  ! degreas north
     !
     ! 1: classic
     REAL(DP) :: height =     1000.0_DP  ! hPa (default) / Pa / m ASL / m AGL
     CHARACTER(LEN=STRLEN_MEDIUM) ::   hunit = ''
     REAL(DP) :: emission =      0.0_DP  ! kg / kg/s / mol/s / molec/s
     CHARACTER(LEN=STRLEN_MEDIUM) ::   eunit = ''
     INTEGER, DIMENSION(6) :: ds = (/ 0,0,0,0,0,0 /) ! START: YY, MM, DD, ...
     INTEGER, DIMENSION(6) :: de = (/ 0,0,0,0,0,0 /) ! END:   ... HR, MI, SE
     !
     ! 2: external time series
     CHARACTER(LEN=STRLEN_CHANNEL) :: channel      = ''
     CHARACTER(LEN=STRLEN_OBJECT)  :: obj_height   = ''
     CHARACTER(LEN=STRLEN_OBJECT)  :: obj_flux     = ''
     REAL(DP)                      :: scal         =    1.0_DP
     !
     ! shared
     CHARACTER(LEN=MAXSTRLEN) :: trset  = ''
     CHARACTER(LEN=MAXSTRLEN) :: trlist = ''
     !
  END TYPE T_RELEASE_POINT_IO

  ! WORKSPACE
  ! ... reaction partner, unit conversion
  TYPE T_RELEASED_TRACER
     TYPE(T_RELEASED_TRACER_IO)          :: io
     REAL(DP), DIMENSION(:,:,:), POINTER :: rp   => NULL() ! REACTION PARTNER
     REAL(DP), DIMENSION(:,:,:), POINTER :: rpte => NULL() ! TENDENCY
     LOGICAL                             :: lv2c = .FALSE. ! UNIT CONVERSION ?
  END TYPE T_RELEASED_TRACER

  ! ... cross reference to release points
  TYPE T_RELEASED_TRACER_XREF
     INTEGER                             :: idt    ! TRACER ID
     INTEGER                             :: np = 0 ! NUMBER OF RELEAESE POINTS
     INTEGER, DIMENSION(:), POINTER      :: ixp => NULL() ! REL.POINT INDICES
  END TYPE T_RELEASED_TRACER_XREF

  ! ... release points in decomposition; timing
  TYPE T_RELEASE_POINT
     TYPE(T_RELEASE_POINT_IO) :: io
     !
     ! shared
     INTEGER                :: jgx    = 0   ! longitude index (global)
     INTEGER                :: jgy    = 0   ! latitude index (global)
     !
     INTEGER                :: pe     = 0   ! on which CPU [0...p_ncpus]
     INTEGER                :: jp     = 0   ! which column [1...kproma]
     INTEGER                :: jrow   = 0   ! which row    [1...jrow]
     INTEGER                :: ierr   = 0   ! status
     !
     ! 1: classic
     REAL(DP)               :: dt     = 0   ! emission time in seconds
     LOGICAL                :: lnow   = .FALSE. ! EMISSION NOW ?
     REAL(DP), POINTER      :: rnow   => NULL() ! -> channel object
     REAL(DP)               :: estart = 0._dp   ! time of emission start (jd)
     REAL(DP)               :: estop  = 0._dp   ! time when emission ends (jd)
#ifdef _TFCORR
     LOGICAL                :: lfirst = .FALSE.   ! first time step of emission
     LOGICAL                :: l2nd   = .FALSE.  ! 2nd time step of emission
#endif
     !
     INTEGER                :: height_unit   = HU_undef
     REAL(DP)               :: height_scal   = 100._dp ! standard hPa -> Pa
     INTEGER                :: emission_unit = EU_undef
     !
     ! 2: external time series
     REAL(DP), DIMENSION(:), POINTER :: height => NULL() ! [Pa], [m]
     REAL(DP), DIMENSION(:), POINTER :: flux   => NULL() ! [molec/s]
     !
  END TYPE T_RELEASE_POINT

  ! ... cross reference to tracers
  TYPE T_RELEASE_POINT_XREF
     INTEGER                         :: ntrac = 0      ! NUMBER OF ASSOC. TRAC.
     INTEGER,  DIMENSION(:), POINTER :: idt  => NULL() ! TRACER IDs
     REAL(DP), DIMENSION(:), POINTER :: flx => NULL()  ! [mol/s] emission flux
     REAL(DP), DIMENSION(:), POINTER :: molarmass => NULL()
     LOGICAL,  DIMENSION(:), POINTER :: ltr => NULL()  ! LOCAL TRACER ASSOC. ?
     TYPE(PTR_3D_ARRAY), DIMENSION(:), POINTER :: emis => NULL()
  END TYPE T_RELEASE_POINT_XREF

  ! WORKSPACE
  ! ... representation independent
  TYPE(T_RELEASED_TRACER_IO),    DIMENSION(NMAXTRAC)  :: TR
  TYPE(T_RELEASED_TRACER), ALLOCATABLE,  DIMENSION(:) :: XTR
  TYPE(T_RELEASE_POINT_IO),     DIMENSION(NMAXPOINT)  :: POINT
  TYPE(T_RELEASE_POINT),   ALLOCATABLE,  DIMENSION(:) :: XPOINT
  ! allow point emissions also for non-TREXP-tracers ?
  LOGICAL :: l_force_emis = .FALSE.
  ! allow reaction / decay also for non-TREX-tracers ?
  LOGICAL :: l_allow_ext = .FALSE.
  !
  LOGICAL :: l_tf_corr = .FALSE. ! time filter correction for classic point
  !                              ! sources
  !
  TYPE(T_RELEASE_POINT_XREF),   DIMENSION(:), ALLOCATABLE :: XXREF
  TYPE(T_RELEASED_TRACER_XREF), DIMENSION(:), ALLOCATABLE :: XTRAC
  !
  ! ACTUAL NUMBERS
  INTEGER :: NTRAC  = 0  ! number of tracers
  INTEGER :: NPOINT = 0  ! number of release points

  ! REPRESENTATION SWITCHES
  LOGICAL :: L_GP = .TRUE.  ! grid-point
  LOGICAL :: L_LG = .TRUE.  ! Lagrangian ATTILA
  LOGICAL :: L_CL = .TRUE.  ! Lagrangian CLaMS
  !

#ifdef MESSYTENDENCY
  INTEGER :: my_handle
#endif

#ifdef ECHAM5
!!#D attila +
  ! SPECIAL: LAGRANGIAN WORKSPACE ATTILA
  REAL(DP), DIMENSION(:),   POINTER :: iplon  => NULL() ! index
  REAL(DP), DIMENSION(:),   POINTER :: iplat  => NULL() ! index
  REAL(DP), DIMENSION(:),   POINTER :: iplev  => NULL() ! index
  REAL(DP),                 POINTER :: amcell => NULL() ! kg
!!#D attila -

!!#D clams +
  ! SPECIAL: LAGRANGIAN WORKSPACE CLaMS
  REAL(DP), DIMENSION(:,:), POINTER :: POS       => NULL() ! position index
  ! number of parcels per echam grid box in parallel decomposition
  REAL(DP), DIMENSION(:,:,:),   POINTER :: spc_g => NULL()
  REAL(DP), DIMENSION(:,:,:),   POINTER :: pc_g  => NULL()   ! global field
!!#D clams -
#endif
  REAL(DP), DIMENSION(:,:,:),   POINTER :: temp_3d => NULL() ! temp. [K]
  REAL(DP), DIMENSION(:,:,:),   POINTER :: cair_3d => NULL() ! 1/cm3

  PUBLIC :: trexp_initialize
  PUBLIC :: trexp_new_tracer
  PUBLIC :: trexp_init_memory
  PUBLIC :: trexp_init_coupling
  !         ! -> init_trac_cpl
  !         ! -> init_emis_cpl
  PUBLIC :: trexp_global_start
  PUBLIC :: trexp_physc
  PUBLIC :: trexp_global_end
  PUBLIC :: trexp_free_memory

  !PRIVATE: trexp_read_nml_cpl

CONTAINS

  ! -------------------------------------------------------------------
  SUBROUTINE trexp_initialize

    ! trexp MODULE ROUTINE (ECHAM-5 INTERFACE)
    !
    ! INITIALIZATION OF GLOBAL VARIABLES FROM NAMELIST
    ! IN PARALLEL ENVIRONMENT
    !
    ! Author: Patrick Joeckel, MPICH, Dec 2004

    USE messy_main_mpi_bi,        ONLY: p_parallel_io, p_io, p_bcast
    USE messy_main_blather_bi,    ONLY: error_bi
    USE messy_main_tools,         ONLY: find_next_free_unit  &
                                      , strcrack
    USE messy_main_timer,         ONLY: time_span_d, gregor2julian
#ifdef ICON
    USE messy_main_channel_bi,    ONLY: n_dom
#endif

    IMPLICIT NONE

    INTRINSIC :: TRIM, TINY, ABS
#ifdef ICON
    INTRINSIC :: INDEX
#endif
    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'trexp_initialize'
    INTEGER                     :: iou    ! I/O unit
    INTEGER                     :: status ! error status
    INTEGER                     :: k, jt, i, num
#ifdef ICON
    INTEGER                     :: idx, n
#endif
    ! TRACER SETS
    CHARACTER(LEN=STRLEN_TRSET), DIMENSION(:), POINTER :: trs
    CHARACTER(LEN=STRLEN_TRSET)                        :: trsname
    !
    INTEGER                     :: noe
    INTEGER                     :: ITRAC
    INTEGER                     :: IPOINT

    NULLIFY(trs)

    ! INITIALIZE CPL
    IF (p_parallel_io) THEN
       iou = find_next_free_unit(100,200)
       CALL trexp_read_nml_cpl(status, iou)
       IF (status /= 0) CALL error_bi(' ',substr)
    END IF
    CALL p_bcast(L_GP, p_io)
    CALL p_bcast(L_LG, p_io)
    CALL p_bcast(L_CL, p_io)
    CALL p_bcast(l_force_emis, p_io)
    CALL p_bcast(l_allow_ext, p_io)
    CALL p_bcast(l_tf_corr, p_io)

    CALL start_message_bi(modstr,'INITIALISATION ',substr)

    IF (p_parallel_io) THEN
       WRITE(*,*) '======================================================'
       WRITE(*,*) '                   TRACERS                            '
       WRITE(*,*) '======================================================'
    END IF

    IF (p_parallel_io) THEN
       ! ############### TRACERS #################################
       ! GET NUMBER OF NEW TRACERS
       NTRAC = 0
       DO jt=1, NMAXTRAC

          IF (TRIM(TR(jt)%trset) == '') CYCLE
          IF (TRIM(TR(jt)%trname) == '') CYCLE

          ! CREATE INDIVIDUAL TRACER ENTRY FOR EACH TRACER SET IN LIST
          CALL strcrack(TRIM(TR(jt)%trset), ';', trs, num)
#if defined(ICON)
!! qqq do this only for GPTRSTR???
          ! check for specific tracer sets (domain specific)
          noe = 0
          DO k=1,  num
             trsname = trs(k)
             idx = INDEX(trsname, '_D')
             IF (idx == 0) THEN
                ! if domain is unspecified, use it as wildcard for all domains;
                ! -> add n_dom entries
                noe = noe + n_dom
             ELSE
                ! only one specific domain
                noe = noe + 1
             END IF
          END DO
#else
          noe = num
#endif
          NTRAC = NTRAC + noe
          IF (ASSOCIATED(trs)) THEN
             DEALLOCATE(trs); NULLIFY(trs)
          END IF
          WRITE(*,*) 'NO. OF TRACERS: ',noe
          WRITE(*,*) '-----------------------------------------------'
       END DO
       WRITE(*,*) '-----------------------------------------------'
       WRITE(*,*) 'TOTAL NUMBER OF TRACERS: ',NTRAC
       WRITE(*,*) '-----------------------------------------------'
    END IF
    CALL p_bcast(NTRAC, p_io)

    ALLOCATE(XTR(NTRAC))
    ALLOCATE(XTRAC(NTRAC))

    IF (p_parallel_io) THEN
       ITRAC = 0
       ! COPY DATA AND PARSE STR
       tracer_loop: DO jt=1, NMAXTRAC

          IF (TRIM(TR(jt)%trset) == '') CYCLE
          IF (TRIM(TR(jt)%trname) == '') CYCLE

          ! CREATE INDIVIDUAL TRACER ENTRY FOR EACH TRACER SET IN LIST
          CALL strcrack(TRIM(TR(jt)%trset), ';', trs, num)

#if defined(ICON)
!! qqq do this only for GPTRSTR???
          ! check for specific tracer sets (domain specific)
          DO k=1,  num
             trsname = trs(k)
             idx = INDEX(trsname, '_D')
             IF (idx == 0) THEN
                ! if domain is unspecified, use it as wildcard for all domains;
                ! -> add n_dom entries
                DO n=1, n_dom
                   ITRAC = ITRAC+1
                   IF (ITRAC > NTRAC) &
                        CALL error_bi('tracer counter out of bounds',substr)
                   CALL copy2x_tr(status,XTR(ITRAC),TR(jt),trsname,n)
                   IF (status /= 0) &
                        CALL error_bi('copy2x_tr reported an error ',substr)
                END DO
             ELSE
                ITRAC = ITRAC+1
                IF (ITRAC > NTRAC) &
                     CALL error_bi('nudging counter out of bounds',substr)
                CALL copy2x_tr(status, XTR(ITRAC),TR(jt),trsname)
                IF (status /= 0) &
                     CALL error_bi('copy2x_tr reported an error ',substr)
             END IF
          END DO
#else
          DO k=1, num
             trsname = trs(k)
             ITRAC = ITRAC+1
             IF (ITRAC > NTRAC) &
                  CALL error_bi('tracer counter out of bounds',substr)
             CALL copy2x_tr(status, XTR(ITRAC),TR(jt),trsname)
             IF (status /= 0) &
                  CALL error_bi('copy2x_tr reported an error ',substr)
          END DO
#endif
          IF (ASSOCIATED(trs)) THEN
             DEALLOCATE(trs); NULLIFY(trs)
          END IF

       END DO tracer_loop
    END IF

    ! BROADCAST RESULTS
    DO jt=1, NTRAC
       CALL p_bcast(XTR(jt)%io%trset,   p_io)
       CALL p_bcast(XTR(jt)%io%trname,  p_io)
       CALL p_bcast(XTR(jt)%io%order,   p_io)
       CALL p_bcast(XTR(jt)%io%ka,      p_io)
       CALL p_bcast(XTR(jt)%io%Ta,      p_io)
       CALL p_bcast(XTR(jt)%io%channel, p_io)
       CALL p_bcast(XTR(jt)%io%object,  p_io)
       ! idt             -> set in trexp_new_tacer
       ! np, ixp         -> set in trexp_init_coupling
       ! rp, rpte, lv2c  -> set in trexp_init_coupling
    END DO

    IF (p_parallel_io) THEN
       WRITE(*,*) ' ---> ',NTRAC,' TRACER(S) SPECIFIED !'
    END IF

    IF (p_parallel_io) THEN
       WRITE(*,*) '======================================================'
       WRITE(*,*) '               RELEASE POINTS                         '
       WRITE(*,*) '======================================================'
    END IF

    IF (p_parallel_io) THEN
       ! ############### RELEASE POINTS #################################
       ! GET NUMBER OF RELEASE POINTS
       NPOINT = 0
       ! COPY DATA AND PARSE STR
       point_loop1: DO i=1, NMAXPOINT
          ! TYPE
          IF (POINT(i)%type == 0) CYCLE  ! OFF
          ! TRACER SET and TRACER
          IF (TRIM(POINT(i)%trset)  == '') CYCLE
          IF (TRIM(POINT(i)%trlist) == '') CYCLE
          ! CREATE ONE INDIVIDUAL ENTRY FOR EACH TRACER SET,
          ! HOWEVER, KEEP THE LIST OF TRACERS PER TRACER SET
          CALL strcrack(TRIM(POINT(i)%trset), ';', trs, num)
#if defined(ICON)
!! qqq do this only for GPTRSTR???
          ! check for specific domains
          noe = 0
          DO k=1,  num
             trsname = trs(k)
             idx = INDEX(trsname, '_D')
             IF (idx == 0) THEN
                ! if domain is unspecified, use it as wildcard for all domains;
                ! -> add n_dom entries
                noe = noe + n_dom
             ELSE
                ! only one specific domain
                noe = noe + 1
             END IF
          END DO
#else
          noe = num
#endif
          IF (ASSOCIATED(trs)) THEN
             DEALLOCATE(trs); NULLIFY(trs)
          END IF
          !
          NPOINT = NPOINT + noe
       END DO point_loop1

    END IF
    CALL p_bcast(NPOINT, p_io)

    ALLOCATE(XPOINT(NPOINT))
    ALLOCATE(XXREF(NPOINT))

    IF (p_parallel_io) THEN
       ! ############### RELEASE POINTS #################################
       ! GET NUMBER OF RELEASE POINTS
       IPOINT = 0
       ! COPY DATA AND PARSE STR
       point_loop: DO i=1, NMAXPOINT
          !
          ! SHARED
          !
          ! TYPE
          IF (POINT(i)%type == 0) CYCLE  ! OFF

          WRITE(*,*) 'NUMBER / TYPE           : ' &
                  , '(',i,') / ', POINT(i)%type

          IF (TRIM(POINT(i)%trset) == '') THEN
             WRITE(*,*) '    ... empty tracer set list ... skipping ...'
             CYCLE
          END IF

          IF (TRIM(POINT(i)%trlist) == '') THEN
             WRITE(*,*) '    ... empty tracer list ... skipping ...'
             CYCLE
          END IF
          !
          ! CREATE ONE SPECIFIC ENTRY FOR EACH TRACER SET,
          ! HOWEVER, KEEP THE LIST OF TRACERS PER TRACER SET
          CALL strcrack(TRIM(POINT(i)%trset), ';', trs, num)
#if defined(ICON)
!! qqq do this only for GPTRSTR???
          ! check for specific domains
          DO k=1,  num
             trsname = trs(k)
             idx = INDEX(trsname, '_D')
             IF (idx == 0) THEN
                DO n=1, n_dom
                   IPOINT = IPOINT + 1
                   IF (IPOINT > NPOINT) &
                        CALL error_bi('point counter out of bounds',substr)
                   CALL copy2x_pt(status &
                        , XPOINT(IPOINT),POINT(i),trsname,n)
                   IF (status /= 0) &
                        CALL error_bi('copy2x_pt reported an error ',substr)
                END DO
             ELSE
                IPOINT = IPOINT + 1
                IF (IPOINT > NPOINT) &
                     CALL error_bi('point counter out of bounds',substr)
                CALL copy2x_pt(status &
                     , XPOINT(IPOINT),POINT(i),trsname)
                IF (status /= 0) &
                     CALL error_bi('copy2x_pt reported an error ',substr)
             END IF
          END DO
#else
          DO k=1, num
             trsname = trs(k)
             IPOINT = IPOINT + 1
             IF (IPOINT > NPOINT) &
                  CALL error_bi('point counter out of bounds',substr)
             CALL copy2x_pt(status &
                  , XPOINT(IPOINT),POINT(i),trsname)
             IF (status /= 0) &
                  CALL error_bi('copy2x_pt reported an error ',substr)
          END DO
#endif
          IF (ASSOCIATED(trs)) THEN
             DEALLOCATE(trs); NULLIFY(trs)
          END IF

       END DO point_loop
    END IF

    ! BROADCAST RESULTS
    DO i=1, NPOINT
       ! SHARED
       CALL p_bcast(XPOINT(i)%io%type,    p_io)
       CALL p_bcast(XPOINT(i)%io%lon,     p_io)
       CALL p_bcast(XPOINT(i)%io%lat,     p_io)
       CALL p_bcast(XPOINT(i)%io%trset,   p_io)
       CALL p_bcast(XPOINT(i)%io%trlist,  p_io)
       ! jgx, jgy                          -> set in trexp_init_coupling
       ! pe, jp, jrow, ierr                -> set in trexp_init_coupling
       !
       ! CLASSIC
       CALL p_bcast(XPOINT(i)%io%height,   p_io)
       CALL p_bcast(XPOINT(i)%io%hunit,    p_io)
       CALL p_bcast(XPOINT(i)%io%emission, p_io)
       CALL p_bcast(XPOINT(i)%io%eunit,    p_io)
       CALL p_bcast(XPOINT(i)%io%ds(:),    p_io)
       CALL p_bcast(XPOINT(i)%io%de(:),    p_io)
       CALL p_bcast(XPOINT(i)%dt,          p_io)
       CALL p_bcast(XPOINT(i)%estop,       p_io)
       CALL p_bcast(XPOINT(i)%estart,      p_io)
       ! ntrac, idt, ltr, flx              -> set in trexp_init_coupling
       ! lnow                              -> set in trexp_global_start
       !
       ! UNITS (until now only point sources of type 1,3)
       CALL p_bcast(XPOINT(i)%height_unit,   p_io)
       CALL p_bcast(XPOINT(i)%height_scal,   p_io)
       CALL p_bcast(XPOINT(i)%emission_unit, p_io)
       ! EXTERNAL TIME SERIES
       CALL p_bcast(XPOINT(i)%io%channel,    p_io)
       CALL p_bcast(XPOINT(i)%io%obj_height, p_io)
       CALL p_bcast(XPOINT(i)%io%obj_flux,   p_io)
       CALL p_bcast(XPOINT(i)%io%scal,       p_io)
       ! height, flux                      -> set in trexp_init_coupling
       !
    END DO

    IF (p_parallel_io) THEN
       WRITE(*,*) ' ---> ',NPOINT,' RELEASE POINT(S) SPECIFIED !'
       WRITE(*,*) '------------------------------------------------------'
    END IF

    CALL end_message_bi(modstr,'INITIALISATION ',substr)

  CONTAINS

    SUBROUTINE copy2x_tr(status, X, T, trs, dom)

      USE messy_main_tools,       ONLY: int2str, to_upper

      IMPLICIT NONE

      ! I/O
      INTEGER,                     INTENT(OUT) :: status
      TYPE(T_RELEASED_TRACER),     INTENT(OUT) :: X
      TYPE(T_RELEASED_TRACER_IO),  INTENT(IN)  :: T
      CHARACTER(LEN=STRLEN_TRSET), INTENT(IN)  :: trs
      INTEGER, OPTIONAL,           INTENT(IN)  :: dom

      ! LOCAL
      CHARACTER(LEN=2) :: dstr

      status = 0

      ! copy everything
      X%io = T

      ! set tracer set (depending on domain)
      IF (PRESENT(dom)) THEN
         CALL int2str(dstr, dom)
         X%io%trset   = TRIM(trs)//'_D'//dstr
      ELSE
         X%io%trset   = TRIM(trs)
      END IF

      WRITE(*,*) '-----------------------------------------------'
      WRITE(*,*) 'TRACER SET     : '//TRIM(X%io%trset)
      WRITE(*,*) 'TRACER         : '//TRIM(X%io%trname)

      SELECT CASE(X%io%order)
      CASE(0)
         WRITE(*,*) '... 0-ORDER REACTION (DECAY) WITH ...'
         WRITE(*,*) '    ... ka = ',X%io%ka,' 1/s'
         ! Ta, channel, object -> DEFAULT
      CASE(1)
         WRITE(*,*) '... 1st-ORDER REACTION (k=ka*exp(-Ta/T)) WITH ...'
         WRITE(*,*) '    ... ka = ',X%io%ka,' cm^3/s'
         WRITE(*,*) '    ... Ta = ',X%io%Ta,' K'
         !
         IF (TRIM(X%io%channel) == '') THEN
            WRITE(*,*) ' *** ERROR *** undefined reaction partner (channel)'
            status = 1
            RETURN
         END IF
         !
         IF (TRIM(X%io%object) == '') THEN
            WRITE(*,*) ' *** ERROR *** undefined reaction partner (object)'
            status = 2
            RETURN
         END IF
         WRITE(*,*) '    ... Y  = ''',TRIM(X%io%channel),''' ' &
              ,'''',TRIM(X%io%object),''''
      CASE(-1)
         WRITE(*,*) '... 1st-ORDER REACTION (k=ka+c_air*kb) WITH ...'
         WRITE(*,*) '    ... ka = ',X%io%ka,' cm^3/s'
         WRITE(*,*) '    ... kb = ',X%io%Ta,' cm^6/s'
         !
         IF (TRIM(X%io%channel) == '') THEN
            WRITE(*,*) ' *** ERROR *** undefined reaction partner (channel)'
            status = 3
            RETURN
         END IF
         !
         IF (TRIM(X%io%object) == '') THEN
            WRITE(*,*) ' *** ERROR *** undefined reaction partner (object)'
            status = 4
            RETURN
         END IF
         !
         WRITE(*,*) '    ... Y  = ''',TRIM(X%io%channel),''' ' &
              ,'''',TRIM(X%io%object),''''
      CASE(-2)
         WRITE(*,*) '... special case (CH4 + OH)'
         IF (TRIM(X%io%channel) == '') THEN
            WRITE(*,*) ' *** ERROR *** undefined reaction partner (channel)'
            status = 5
            RETURN
         END IF
         !
         IF (TRIM(X%io%object) == '') THEN
            WRITE(*,*) ' *** ERROR *** undefined reaction partner (object)'
            status = 6
            RETURN
         END IF
      CASE(-3)
         WRITE(*,*) '... special case (SO2 + OH)'
         IF (TRIM(X%io%channel) == '') THEN
            WRITE(*,*) ' *** ERROR *** undefined reaction partner (channel)'
            status = 7
            RETURN
         END IF
         IF (TRIM(X%io%object) == '') THEN
            WRITE(*,*) ' *** ERROR *** undefined reaction partner (object)'
            status = 8
            RETURN
         END IF
      CASE DEFAULT
         WRITE(*,*) ' *** ERROR *** unknown ORDER OF REACTION (='  &
                     , X%io%order,')'
         status = 9
         RETURN
      END SELECT

    END SUBROUTINE copy2x_tr

    SUBROUTINE copy2x_pt(status, X, P, trs, dom)

      USE messy_main_tools,       ONLY: int2str, to_upper

      IMPLICIT NONE

      ! I/O
      INTEGER,                     INTENT(OUT) :: status
      TYPE(T_RELEASE_POINT),       INTENT(OUT) :: X
      TYPE(T_RELEASE_POINT_IO),    INTENT(IN)  :: P
      CHARACTER(LEN=STRLEN_TRSET), INTENT(IN)  :: trs
      INTEGER, OPTIONAL,           INTENT(IN)  :: dom

      ! LOCAL
      CHARACTER(LEN=2) :: dstr

      status = 0

      ! copy everything
      X%io = P

      ! set tracer set (depending on domain)
      IF (PRESENT(dom)) THEN
         CALL int2str(dstr, dom)
         X%io%trset   = TRIM(trs)//'_D'//dstr
      ELSE
         X%io%trset   = TRIM(trs)
      END IF

      WRITE(*,*) '-----------------------------------------------'
      WRITE(*,*) 'TRACER SET     : ',TRIM(X%io%trset)
      WRITE(*,*) 'TRACER LIST    : ',TRIM(X%io%trlist)
      WRITE(*,*) 'TYPE           : ',X%io%type
      WRITE(*,*) 'LONGITUDE      : ',X%io%lon
      IF ((X%io%lon < -180.0_DP) .OR. (POINT(i)%lon > 180.0_DP)) THEN
         WRITE(*,*) ' *** ERROR *** LONGITUDE out of range [-180, 180]'
         status = 1
         RETURN
      END IF
      WRITE(*,*) 'LATITUDE       : ',X%io%lat
      IF ((X%io%lat < -90.0_DP) .OR. (POINT(i)%lat > 90.0_DP)) THEN
         WRITE(*,*) ' *** ERROR *** LATITUDE out of range [-90, 90]'
         status = 2
         RETURN
      END IF
             !
      SELECT CASE(X%io%type)
      CASE (0)
         ! OFF (see above)
      CASE (1, 3) ! op_cf_20150414: type 3
         ! CLASSIC
         ! evaluate emission height units and emission units
         SELECT CASE(TRIM(to_upper(X%io%hunit)))
         CASE ('HPA')
            X%height_unit = HU_hPa
            X%height_scal  = 100._dp
         CASE ('PA')
            X%height_unit = HU_Pa
            X%height_scal  = 1._dp
         CASE ('M', 'M ASL', 'MASL', 'M A.S.L')
            X%height_unit = HU_mASL
            X%height_scal  = 1._dp
         CASE ('M AGL', 'MAGL', 'M A.G.L')
            X%height_unit = HU_mAGL
            X%height_scal  = 1._dp
         CASE DEFAULT
            WRITE(*,*) ' *** ERROR *** Unknown unit for emission height: '//TRIM(X%io%hunit)
            status = 3
            RETURN
         END SELECT
         SELECT CASE(TRIM(to_upper(X%io%eunit)))
         CASE ('KG')
            X%emission_unit = EU_kg
         CASE ('KG/S', 'KG S-1')
            X%emission_unit = EU_kgs
         CASE ('MOL/S', 'MOL S-1')
            X%emission_unit = EU_mols
         CASE ('MOLEC/S', 'MOLECULES/S', 'MOLEC S-1')
            X%emission_unit = EU_molecs
         CASE DEFAULT
            WRITE(*,*) ' *** ERROR *** Unknown unit for emission: '//TRIM(X%io%eunit)
            status = 4
            RETURN
         END SELECT
         ! write information on point emission
         WRITE(*,*) 'HEIGHT         : ', X%io%height &
            &      ,' ['//TRIM(HU_STRING(X%height_unit))//']'
         WRITE(*,*) 'EMISSION       : ', X%io%emission &
            &      ,' ['//TRIM(EU_STRING(X%emission_unit))//']'
         WRITE(*,*) 'START          : '  &
              , X%io%ds(1), '-', X%io%ds(2), '-', X%io%ds(3) ,' ' &
              , X%io%ds(4), ':', X%io%ds(5), ':', X%io%ds(6)
         X%estart = gregor2julian( &
              X%io%ds(1), X%io%ds(2), X%io%ds(3), X%io%ds(4), &
              X%io%ds(5), X%io%ds(6))
         WRITE(*,*) 'END            : ' &
              , X%io%de(1), '-', X%io%de(2) ,'-', X%io%de(3) ,' ' &
              , X%io%de(4), ':', X%io%de(5) ,':', X%io%de(6)
         X%estop = gregor2julian( &
              X%io%de(1), X%io%de(2), X%io%de(3), X%io%de(4), &
              X%io%de(5), X%io%de(6))

         ! CALCULATE TIME INTERVAL IN SECONDS
         CALL time_span_d(X%dt, &
              X%io%ds(1), X%io%ds(2), X%io%ds(3) &
              , X%io%ds(4), X%io%ds(5), X%io%ds(6) &
              , X%io%de(1), X%io%de(2), X%io%de(3) &
              , X%io%de(4), X%io%de(5), X%io%de(6) &
              )
         ! days -> seconds
         X%dt = X%dt * 86400.0_dp
         WRITE(*,*) 'dT             : ', X%dt,' [s]'
         IF (X%dt < 0.0_DP) THEN
            WRITE(*,*) ' *** ERROR *** TIME INTERVAL < 0 : ', X%dt
            status = 5
            RETURN
         END IF
         IF (ABS(X%dt) < TINY(0.0_DP)) THEN
            WRITE(*,*) ' *** ERROR *** TIME INTERVAL ~ 0 : ', X%dt
            status = 6
            RETURN
         END IF

      CASE (2, 4) ! op_cf_20150414: type 4
         ! EXTERNAL TIME SERIES
         !
         ! SCALING
         WRITE(*,*) 'SCALING        : ', X%io%scal
         !
         ! CHANNEL (WITH TIME SERIES) AND OBJECTS
         IF (TRIM(X%io%channel) == '') THEN
            WRITE(*,*) ' *** ERROR *** empty CHANNEL'
            status = 7
            RETURN
         ELSE
            WRITE(*,*) 'CHANNEL        : ', TRIM(X%io%channel)
         END IF
         !
         IF (TRIM(X%io%obj_height) == '') THEN
            WRITE(*,*) ' *** ERROR *** empty OBJECT (height)'
            status = 8
            RETURN
         ELSE
            WRITE(*,*) 'OBJECT (height)  : ', TRIM(X%io%obj_height)
         END IF
         !
         IF (TRIM(X%io%obj_flux) == '') THEN
            WRITE(*,*) ' *** ERROR *** empty OBJECT (flux)'
            status = 9
            RETURN
         ELSE
            WRITE(*,*) 'OBJECT (flux)  : ', TRIM(X%io%obj_flux)
         END IF
                !
      CASE DEFAULT
         WRITE(*,*) ' *** ERROR *** unknown type ('&
              , X%io%type
         status = 10
         RETURN
      END SELECT
      WRITE(*,*) '-----------------------------------------------'

    END SUBROUTINE copy2x_pt

  END SUBROUTINE trexp_initialize
  ! -------------------------------------------------------------------

  ! -------------------------------------------------------------------
  SUBROUTINE trexp_new_tracer

    ! trexp MODULE ROUTINE (ECHAM-5 INTERFACE, PRIVATE)
    !
    ! defines specific tracers
    !
    ! Author: Patrick Joeckel, MPICH, December 2004

    USE messy_main_mpi_bi,          ONLY: p_parallel_io
    USE messy_main_tracer_mem_bi,   ONLY: GPTRSTR, LGTRSTR, CLTRSTR
    USE messy_main_tracer_tools_bi, ONLY: tracer_halt
    USE messy_main_blather_bi,      ONLY: error_bi, warning_bi
    USE messy_main_tracer,          ONLY: new_tracer, set_tracer, R_molarmass &
                                        , get_tracer, full2base_sub

    IMPLICIT NONE

    INTRINSIC :: TRIM

    ! LOCAL
    INTEGER :: status
    CHARACTER(LEN=*), PARAMETER :: substr = 'trexp_new_tracer'
    INTEGER :: jt
    CHARACTER(LEN=STRLEN_MEDIUM) :: basename    = '' ! name of tracer
    CHARACTER(LEN=STRLEN_MEDIUM) :: subname     = '' ! OPTIONAL subname
    REAL(dp)                     :: zmolarmass

    CALL start_message_bi(modstr, 'TRACER REQUEST', substr)

#ifdef MESSYTENDENCY
    my_handle = mtend_get_handle(modstr)
#endif

    DO jt=1, NTRAC

#ifdef ICON
       IF (TRIM(XTR(jt)%io%trset) /= GPTRSTR) CYCLE
#endif
       IF ((.NOT. L_GP) .AND. (TRIM(XTR(jt)%io%trset) == TRIM(GPTRSTR))) &
            CYCLE
       IF ((.NOT. L_LG) .AND. (TRIM(XTR(jt)%io%trset) == TRIM(LGTRSTR))) &
            CYCLE
       IF ((.NOT. L_CL) .AND. (TRIM(XTR(jt)%io%trset) == TRIM(CLTRSTR))) &
            CYCLE

       IF (p_parallel_io) THEN
          WRITE(*,*) 'TRACER SET : ', TRIM(XTR(jt)%io%trset)
          WRITE(*,*) 'TRACER     : ', TRIM(XTR(jt)%io%trname)
       END IF

       ! does tracer already exist (e.g., defined in ptrac)?
       CALL full2base_sub(status, TRIM(XTR(jt)%io%trname) &
            , basename, subname)
       CALL tracer_halt(substr, status)

       CALL get_tracer(status, TRIM(XTR(jt)%io%trset), basename &
            , subname=subname, idx = XTRAC(jt)%idt)
       IF (status /= 0) THEN
          ! NEW TRACER
          CALL new_tracer(status, TRIM(XTR(jt)%io%trset) &
               , basename, modstr       &
               , subname = subname      &
               , idx = XTRAC(jt)%idt    &
               , unit='mol/mol')
          CALL tracer_halt(substr, status)
          CALL get_tracer(status, TRIM(XTR(jt)%io%trset), XTRAC(jt)%idt &
               , R_molarmass, zmolarmass)
          CALL tracer_halt(substr, status)
          IF (zmolarmass <= 0.0_dp) THEN
             CALL warning_bi('setting molar mass of tracer '//&
                  &TRIM(XTR(jt)%io%trname)//' in tracer set '//&
                  &TRIM(XTR(jt)%io%trset)//' to 1.0 g/mol' , substr)
             CALL set_tracer(status, TRIM(XTR(jt)%io%trset) &
                  , XTRAC(jt)%idt          &
                  , R_molarmass, 1.0_dp)
             CALL tracer_halt(substr, status)
          ENDIF
       ELSE
          ! TRACER EXISTS
          IF (l_allow_ext) THEN
             CALL warning_bi('tracer '//TRIM(XTR(jt)%io%trname)&
                  &//' modified by TREXP', substr)
          ELSE
             CALL error_bi('tracer '//TRIM(XTR(jt)%io%trname)&
                  &//' modified by TREXP', substr)
          ENDIF
       END IF

    END DO

    CALL end_message_bi(modstr, 'TRACER REQUEST', substr)

  END SUBROUTINE trexp_new_tracer
  ! -------------------------------------------------------------------

  ! -------------------------------------------------------------------
  SUBROUTINE trexp_init_memory

!!#D attila +
#ifdef ECHAM5
    USE messy_main_grid_def_mem_bi, ONLY: nproma, nlev, ngpblks
#endif
!!#D attila -

    IMPLICIT NONE

    ! LOCAL
    ! CHARACTER(LEN=*), PARAMETER :: substr = 'trexp_init_memory'

!!#D attila +
#ifdef ECHAM5
    IF (L_LG) THEN
       IF (.NOT. ASSOCIATED(temp_3d)) ALLOCATE(temp_3d(nproma, nlev, ngpblks))
       IF (.NOT. ASSOCIATED(cair_3d)) ALLOCATE(cair_3d(nproma, nlev, ngpblks))
    END IF
#endif
!!#D attila -

!!#D clams +
#ifdef ECHAM5
    IF (L_CL) THEN
       IF (.NOT. ASSOCIATED(temp_3d)) ALLOCATE(temp_3d(nproma, nlev, ngpblks))
       IF (.NOT. ASSOCIATED(cair_3d)) ALLOCATE(cair_3d(nproma, nlev, ngpblks))
    END IF
#endif
!!#D clams -

#ifdef MESSYTENDENCY
    CALL mtend_register(my_handle, mtend_id_tracer)
#endif

  END SUBROUTINE trexp_init_memory
  ! -------------------------------------------------------------------

  ! -------------------------------------------------------------------
  SUBROUTINE trexp_init_coupling

    USE messy_main_mpi_bi,           ONLY: p_parallel_io, p_io, p_bcast
    USE messy_main_blather_bi,       ONLY: error_bi, info_bi
    USE messy_main_tracer_mem_bi,    ONLY: GPTRSTR, LGTRSTR, CLTRSTR &
                                         , gp_channel, lg_channel    &
                                         , cl_channel                &
                                         , xt, xtte
    USE messy_main_tracer_tools_bi,  ONLY: tracer_halt
    USE messy_main_transform_bi,     ONLY: locate_in_decomp
    USE messy_main_channel_error_bi, ONLY: channel_halt
    USE messy_main_tracer,           ONLY: get_tracer              &
                                         , t_trinfo, full2base_sub &
                                         , STRLEN_FNAME
    USE messy_main_tools,            ONLY: strcrack, to_upper
    USE messy_main_constants_mem,    ONLY: STRLEN_ULONG
    USE messy_main_channel,          ONLY: get_channel_object, get_attribute
    USE messy_main_channel,          ONLY: new_channel_object, new_channel &
                                         , new_attribute, STRLEN_OBJECT &
                                         , STRLEN_CHANNEL, get_channel_info
    USE messy_main_channel_error_bi, ONLY: channel_halt
    USE messy_main_channel_bi,       ONLY: SCALAR, GP_3D_MID, LG_ATTILA &
                                         , REPR_LG_CLAMS, REPR_UNDEF
    USE messy_main_tools,            ONLY: int2str

    IMPLICIT NONE

    INTRINSIC :: TRIM

    ! LOCAL
    INTEGER                          :: status
    CHARACTER(LEN=*), PARAMETER      :: substr = 'trexp_init_coupling'
    INTEGER                          :: jt, i
    CHARACTER(LEN=2*STRLEN_MEDIUM+1) :: fullname
    CHARACTER(LEN=STRLEN_MEDIUM)     :: basename, subname
    TYPE(t_trinfo)                   :: trinfo
    INTEGER                          :: ierr
    INTEGER                          :: idx
    !
    CHARACTER(LEN=STRLEN_MEDIUM)     :: unit
    CHARACTER(LEN=STRLEN_ULONG)      :: att_unit
    REAL(DP)                         :: molarmass
    !
    CHARACTER(LEN=4)                 :: istr = ''  ! 0001-9999
    INTEGER                          :: j, idt
    CHARACTER(LEN=STRLEN_FNAME)      :: trfname = ''
    CHARACTER(LEN=STRLEN_OBJECT)     :: objname = ''
    CHARACTER(LEN=STRLEN_CHANNEL)    :: chaname = ''
    INTEGER                          :: reprid

    CALL start_message_bi(modstr, 'COUPLING INITIALIZATION', substr)

    ! INITIALIZE RELEASE POINTS, IF TRACER(S) EXISTS
    IF (p_parallel_io) THEN
       WRITE(*,*) '------------------------------------------------------'
       WRITE(*,*) '                RELEASE POINTS                        '
       WRITE(*,*) '------------------------------------------------------'
    END IF

    point_loop: DO i=1, NPOINT

#ifdef ICON
       IF (TRIM(XPOINT(i)%io%trset) /= TRIM(GPTRSTR)) CYCLE
#endif
       IF ((.NOT. L_GP) .AND. (TRIM(XPOINT(i)%io%trset) == TRIM(GPTRSTR))) &
            CYCLE
       IF ((.NOT. L_LG) .AND. (TRIM(XPOINT(i)%io%trset) == TRIM(LGTRSTR))) &
            CYCLE
       IF ((.NOT. L_CL) .AND. (TRIM(XPOINT(i)%io%trset) == TRIM(CLTRSTR))) &
            CYCLE

       ! LOCATE POSITION
#ifndef ICON
       CALL locate_in_decomp(XPOINT(i)%ierr      &
            , XPOINT(i)%io%lon, XPOINT(i)%io%lat &
            , XPOINT(i)%pe, XPOINT(i)%jp, XPOINT(i)%jrow &
            , XPOINT(i)%jgx, XPOINT(i)%jgy)
#else
       CALL locate_in_decomp(XPOINT(i)%ierr      &
            , XPOINT(i)%io%lon, XPOINT(i)%io%lat &
            , XPOINT(i)%pe, XPOINT(i)%jp, XPOINT(i)%jrow)
#endif
       IF (XPOINT(i)%ierr /= 0) THEN
          IF (p_parallel_io) THEN
             SELECT CASE(XPOINT(i)%ierr)
             CASE(1)
                WRITE(*,*) ' longitude out of range: ', XPOINT(i)%io%lon, &
                     '; skipping POINT ', i
             CASE(2)
                WRITE(*,*) ' latitude out of range: ',  XPOINT(i)%io%lat, &
                     '; skipping POINT ', i
             CASE(555)
                WRITE(*,*) ' POINT not located in model domain: ' &
                     , XPOINT(i)%io%lon, XPOINT(i)%io%lat,         &
                     '; skipping POINT ', i
             END SELECT
          END IF
          CYCLE
       END IF

       ! DIAGNOSTIC OUTPUT
       IF (p_parallel_io) THEN
          WRITE(*,'(1a,i3.3,1x,a5,f9.4,1x,a5,f9.4,1x,a5,i5,1x,a5,i4,'//&
               &'1x,a5,i4,1x,a5,i4,1x,a5,i4)')  &
               '#',i                                                  &
               , ' LAT=',XPOINT(i)%io%lat, ' LON=',XPOINT(i)%io%lon   &
               , '  PE=',XPOINT(i)%pe,     '  JP=', XPOINT(i)%jp      &
               , 'JROW=',XPOINT(i)%jrow                               &
               , ' JGX=',XPOINT(i)%jgx,    ' JGY=',XPOINT(i)%jgy
       ENDIF

       SELECT CASE(XPOINT(i)%io%type)
       CASE(1, 3)
          ! CLASSIC
          IF (p_parallel_io) THEN
             WRITE(*,'(6x,2(a6,i4,a1,i2.2,a1,i2.2,1x,i2.2,a1,i2.2,a1,i2.2,1x)'//&
                  &',a4,e12.4,a4)')  &
                  'START='   &
                  , XPOINT(i)%io%ds(1),'-', XPOINT(i)%io%ds(2) &
                  , '-', XPOINT(i)%io%ds(3) &
                  , XPOINT(i)%io%ds(4),':', XPOINT(i)%io%ds(5) &
                  , ':', XPOINT(i)%io%ds(6) &
                  , '  END=' &
                  , XPOINT(i)%io%de(1),'-', XPOINT(i)%io%de(2) &
                  , '-', XPOINT(i)%io%de(3) &
                  , XPOINT(i)%io%de(4),':', XPOINT(i)%io%de(5) &
                  , ':', XPOINT(i)%io%de(6) &
                  , ' Dt=',XPOINT(i)%dt,' [s]'
          END IF
       CASE(2, 4)
          ! EXTERNAL TIME SERIES
          IF (p_parallel_io) WRITE(*,'(6x,a,a,a,a)') 'height: ' &
               , TRIM(XPOINT(i)%io%channel),' / ', TRIM(XPOINT(i)%io%obj_height)
          CALL get_channel_object(status &
               , TRIM(XPOINT(i)%io%channel) &
               , TRIM(XPOINT(i)%io%obj_height) &
               , p1=XPOINT(i)%height )
          CALL channel_halt(substr, status)
          IF (p_parallel_io) THEN
             CALL get_attribute(status &
                  , TRIM(XPOINT(i)%io%channel), TRIM(XPOINT(i)%io%obj_height) &
                  , 'units',  c=att_unit )
             CALL channel_halt(substr, status)
             SELECT CASE(TRIM(to_upper(att_unit)))
             CASE('HPA')
                XPOINT(i)%height_unit = HU_hPa
                XPOINT(i)%height_scal = 100._dp
             CASE('PA')
                XPOINT(i)%height_unit = HU_Pa
                XPOINT(i)%height_scal = 1._dp
             CASE('M', 'M ASL', 'MASL', 'M A.S.L.')
                XPOINT(i)%height_unit = HU_mASL
                XPOINT(i)%height_scal = 1._dp
             CASE('M AGL', 'MAGL', 'M A.G.L.')
                XPOINT(i)%height_unit = HU_mAGL
                XPOINT(i)%height_scal = 1._dp
             CASE DEFAULT
                CALL error_bi('Unknown unit of height channel object: '//TRIM(att_unit) &
                     , substr)
             END SELECT
          END IF
          CALL p_bcast(XPOINT(i)%height_unit, p_io)
          CALL p_bcast(XPOINT(i)%height_scal, p_io)

          IF (p_parallel_io) WRITE(*,'(6x,a,a,a,a)') 'Flux : ' &
               , TRIM(XPOINT(i)%io%channel),' / ', TRIM(XPOINT(i)%io%obj_flux)
          CALL get_channel_object(status &
               , TRIM(XPOINT(i)%io%channel) &
               , TRIM(XPOINT(i)%io%obj_flux) &
               , p1=XPOINT(i)%flux )
          CALL channel_halt(substr, status)
          !
          IF (p_parallel_io) THEN
             ! Test concistency of flux and height object obtained
             ! there has to be one more of the height bonudaries than
             ! the vertical emission boxes of the flux
             IF (SIZE(XPOINT(i)%flux, 1) + 1 /= SIZE(XPOINT(i)%height, 1)) THEN
                CALL error_bi("Vertical axis of '"//TRIM(XPOINT(i)%io%obj_height) &
                   //"' has to be one element larger than vertical axis of '"      &
                   //TRIM(XPOINT(i)%io%obj_flux)//"' on channel '"                  &
                   //TRIM(XPOINT(i)%io%channel)//"'" &
                     , substr)
             END IF
             CALL get_attribute(status &
                  , TRIM(XPOINT(i)%io%channel), TRIM(XPOINT(i)%io%obj_flux) &
                  , 'units',  c=att_unit )
             CALL channel_halt(substr, status)
             SELECT CASE(TRIM(to_upper(att_unit)))
             CASE('KG/S', 'KG S-1')
                XPOINT(i)%emission_unit = EU_kgs
             CASE('MOL/S', 'MOL S-1')
                XPOINT(i)%emission_unit = EU_mols
             CASE('MOLEC/S', 'MOLECULES/S', 'MOLEC S-1')
                XPOINT(i)%emission_unit = EU_molecs
             CASE DEFAULT
                CALL error_bi('Unknown unit of emission flux channel object: '//TRIM(att_unit) &
                     , substr)
             END SELECT
          END IF
          CALL p_bcast(XPOINT(i)%emission_unit, p_io)
       END SELECT

    END DO point_loop

    ! DEFINE NEW CHANNEL
!qqq in CASE IF ICON with > 1 domain, this will be called multiple time
!    resulting in a error ???
    CALL new_channel(status, modstr)
    CALL channel_halt(substr, status)

    point_loop2: DO i=1, NPOINT

#ifdef ICON
       IF (TRIM(XPOINT(i)%io%trset) /= TRIM(GPTRSTR)) CYCLE
#endif
       IF ((.NOT. L_GP) .AND. (TRIM(XPOINT(i)%io%trset) == TRIM(GPTRSTR))) &
            CYCLE
       IF ((.NOT. L_LG) .AND. (TRIM(XPOINT(i)%io%trset) == TRIM(LGTRSTR))) &
            CYCLE
       IF ((.NOT. L_CL) .AND. (TRIM(XPOINT(i)%io%trset) == TRIM(CLTRSTR))) &
            CYCLE

       IF (XPOINT(i)%ierr /= 0) CYCLE
       CALL int2str(istr, i, '0', 'x')

       CALL new_channel_object(status, modstr, 'rnow_'//istr &
            , p0=XPOINT(i)%rnow, reprid=SCALAR)
       CALL channel_halt(substr, status)
       CALL new_attribute(status, modstr, 'rnow_'//istr &
            , 'long_name' &
            , c='1 if emission occurs at this timestep, 0 else' )
       CALL channel_halt(substr, status)

    END DO point_loop2

    IF (L_GP) THEN
       IF (p_parallel_io) THEN
          WRITE(*,*) '------------------------------------------------------'
          WRITE(*,*) '                GP TRACER COUPLING                    '
          WRITE(*,*) '------------------------------------------------------'
       END IF
       CALL init_trac_cpl(XXREF, GPTRSTR)
    END IF

!!#D attila +
#ifdef ECHAM5
    IF (L_LG) THEN
       IF (p_parallel_io) THEN
          WRITE(*,*) '------------------------------------------------------'
          WRITE(*,*) '                LG TRACER COUPLING                    '
          WRITE(*,*) '------------------------------------------------------'
       END IF
       CALL init_trac_cpl(XXREF, LGTRSTR)
    END IF
#endif
!!#D attila -

!!#D clams +
#ifdef ECHAM5
    IF (L_CL) THEN
       IF (p_parallel_io) THEN
          WRITE(*,*) '------------------------------------------------------'
          WRITE(*,*) '                CL TRACER COUPLING                    '
          WRITE(*,*) '------------------------------------------------------'
       END IF
       CALL init_trac_cpl(XXREF, CLTRSTR)
    END IF
#endif
!!#D clams -

    ! INITIALIZE POINTERS TO REACTION PARTNERS
    IF (p_parallel_io) THEN
       WRITE(*,*) '------------------------------------------------------'
       WRITE(*,*) '                    REACTIONS                         '
       WRITE(*,*) '------------------------------------------------------'
    END IF

    tracer_loop: DO jt=1, NTRAC
#ifdef ICON
       IF (TRIM(XTR(jt)%io%trset) /= TRIM(GPTRSTR)) CYCLE
#endif
       IF ((.NOT. L_GP) .AND. (TRIM(XTR(jt)%io%trset) == TRIM(GPTRSTR))) &
            CYCLE
       IF ((.NOT. L_LG) .AND. (TRIM(XTR(jt)%io%trset) == TRIM(LGTRSTR))) &
            CYCLE
       IF ((.NOT. L_CL) .AND. (TRIM(XTR(jt)%io%trset) == TRIM(CLTRSTR))) &
            CYCLE

       fullname=TRIM(XTR(jt)%io%trname)

       SELECT CASE(XTR(jt)%io%order)
       CASE(0)
          ! --- 0-ORDER REACTION (DECAY -> NO REACTION PARTNER) ---------
          IF (p_parallel_io) THEN
             WRITE(*,*) ' -> ',TRIM(fullname),' (DECAY)'
          END IF
          XTR(jt)%lv2c = .FALSE.
          ! -------------------------------------------------------------
       CASE(1,-1,-2,-3)
          ! --- 1st-ORDER RACTION ---------------------------------------
          IF (p_parallel_io) THEN
             WRITE(*,*) ' -> ',TRIM(fullname),' + ' &
                  , TRIM(XTR(jt)%io%object),' (',TRIM(XTR(jt)%io%channel),')'
          END IF
          SELECT CASE(TRIM(XTR(jt)%io%channel))
          CASE(gp_channel)
             ! ... GP-TRACER ............................................
             CALL full2base_sub(status, TRIM(XTR(jt)%io%object) &
                  , basename, subname)
             CALL tracer_halt(substr, status)
             !
             CALL get_tracer(ierr, GPTRSTR,  basename &
                  , subname=subname                   &
                  , idx=idx                           &
                  , unit = unit)
             CALL tracer_halt(substr, ierr)
             XTR(jt)%rp   => xt  (_RI_XYZN_(:,:,:,idx))
             XTR(jt)%rpte => xtte(_RI_XYZN_(:,:,:,idx))
             ! ..........................................................
          CASE(lg_channel)
             ! ... LG- TRACER ...........................................
             CALL error_bi(&
                  'LAGRANGIAN TRACER (LG) AS REACTION-PARTNER NOT POSSIBLE' &
                  , substr)
             ! ..........................................................
          CASE(cl_channel)
             ! ... CL- TRACER ...........................................
             CALL error_bi(&
                  'LAGRANGIAN TRACER (CL) AS REACTION-PARTNER NOT POSSIBLE' &
                  , substr)
             ! ..........................................................
          CASE default
             ! ... NON-TRACER CHANNEL OBJECT ............................
             CALL get_channel_object(status &
                  , TRIM(XTR(jt)%io%channel), TRIM(XTR(jt)%io%object) &
                  , p3=XTR(jt)%rp )
             CALL channel_halt(substr, status)
             CALL get_attribute(status &
                  , TRIM(XTR(jt)%io%channel), TRIM(XTR(jt)%io%object) &
                  , 'units',  c=att_unit )
             unit = TRIM(att_unit)
             ! ..........................................................
          END SELECT
          ! CECK UNIT
          SELECT CASE(TRIM(unit))
          CASE('mol/mol','vmr')
             XTR(jt)%lv2c = .TRUE.
             IF (p_parallel_io) THEN
                WRITE(*,*) '    -> CONVERSION MR -> CONC: YES'
             END IF
          CASE('cm^-3','cm-3','cm^(-3)')
             XTR(jt)%lv2c = .FALSE.
             IF (p_parallel_io) THEN
                WRITE(*,*) '    -> CONVERSION MR -> CONC: NO'
             END IF
          CASE('1/s','s-1','s^(-1)')
             XTR(jt)%lv2c = .FALSE.
             IF (p_parallel_io) THEN
                WRITE(*,*) '    UNIT OF REACTION PARTNER IS '''//&
                     &TRIM(unit)//''' !!!'
                WRITE(*,*) '    ASSUMING THAT ka, Ta ARE SET CORRECTLY ...'
                WRITE(*,*) '    -> CONVERSION MR -> CONC: NO'
             END IF
             ! op_fr_20180727: unknown added for SCALC
          CASE('<unspecified>','unknown')
             XTR(jt)%lv2c = .FALSE.
             IF (p_parallel_io) THEN
                WRITE(*,*) '!!! WARNING: UNIT IS '''//TRIM(unit)//''' !!!'
                WRITE(*,*) '    ASSUMING ''cm^(-3)'' ...'
                WRITE(*,*) '    -> CONVERSION MR -> CONC: NO'
             END IF
          CASE DEFAULT
             CALL error_bi('UNKNOWN UNIT OF REACTION PARTNER: '//TRIM(unit) &
                  , substr)
          END SELECT
          ! -------------------------------------------------------------
       CASE DEFAULT
          ! -------------------------------------------------------------
          ! ORDER NOT 0 OR 1 OR -1 OR -2
          ! ... SHOULD NEVER BE REACHED
          ! -------------------------------------------------------------
       END SELECT

    END DO tracer_loop

    ! INITIALIZE CROSS REFERENCE INFORMATION
    IF (p_parallel_io) THEN
       WRITE(*,*) '------------------------------------------------------'
       WRITE(*,*) '                TRACER EMISSIONS                      '
       WRITE(*,*) '------------------------------------------------------'
    END IF
    CALL init_emis_cpl(XXREF, XTRAC)

    IF (p_parallel_io) THEN
       WRITE(*,*) '------------------------------------------------------'
    END IF

!!#D attila +
#ifdef ECHAM5
    IF (L_LG) THEN
       CALL get_channel_object(status, 'attila', 'IPLAT', p1=iplat)
       CALL channel_halt(substr, status)
       CALL get_channel_object(status, 'attila', 'IPLON', p1=iplon)
       CALL channel_halt(substr, status)
       CALL get_channel_object(status, 'attila', 'IPLEV', p1=iplev)
       CALL channel_halt(substr, status)
       CALL get_channel_object(status, 'attila', 'AMCELL', p0=amcell)
       CALL channel_halt(substr, status)
    END IF
#endif
!!#D attila -

!!#D clams +
#ifdef ECHAM5
    IF (L_CL) THEN
       CALL get_channel_object(status, 'clams', 'POS', p2=pos)
       CALL channel_halt(substr, status)
       CALL get_channel_object(status, 'clams', 'SPC_G', p3=spc_g)
       CALL channel_halt(substr, status)
       IF (p_parallel_io) THEN
          WRITE(*,*) '------------------------------------------------------'
          WRITE(*,*) '                POS and SPC_G from clams coupled      '
          WRITE(*,*) '------------------------------------------------------'
       END IF
    END IF
#endif
!!#D clams -

    IF (p_parallel_io) THEN
       WRITE(*,*) '------------------------------------------------------'
       WRITE(*,*) '         EMISSION DIAGNOSTICS (TYPE 3 and 4)          '
       WRITE(*,*) '------------------------------------------------------'
    END IF

    point_loop3: DO i=1, NPOINT

#ifdef ICON
       IF (TRIM(XPOINT(i)%io%trset) /= GPTRSTR) CYCLE
#endif
       IF ((.NOT. L_GP) .AND. (TRIM(XPOINT(i)%io%trset) == TRIM(GPTRSTR))) &
            CYCLE
       IF ((.NOT. L_LG) .AND. (TRIM(XPOINT(i)%io%trset) == TRIM(LGTRSTR))) &
            CYCLE
       IF ((.NOT. L_CL) .AND. (TRIM(XPOINT(i)%io%trset) == TRIM(CLTRSTR))) &
            CYCLE

       IF (XPOINT(i)%ierr /= 0) CYCLE

       chaname = modstr//'_'//TRIM(XPOINT(i)%io%trset)

       reprid = REPR_UNDEF
       IF (TRIM(XPOINT(i)%io%trset) == TRIM(GPTRSTR)) reprid = GP_3D_MID
       IF (TRIM(XPOINT(i)%io%trset) == TRIM(LGTRSTR)) reprid = LG_ATTILA
       IF (TRIM(XPOINT(i)%io%trset) == TRIM(CLTRSTR)) reprid = REPR_LG_CLAMS

       IF (reprid == REPR_UNDEF) &
            CALL error_bi('representation could not be determined',substr)

       ALLOCATE(XXREF(i)%emis(XXREF(i)%ntrac))
       ! ONLY FOR TYPE 3 OR TYPE 4 POINTS
       IF ((XPOINT(i)%io%type /= 3) .AND. (XPOINT(i)%io%type /= 4)) CYCLE

       trac_loop3: DO j=1, XXREF(i)%ntrac
          ! GET TRACER ID
          idt   = XXREF(i)%idt(j)
          IF (idt == 0) CYCLE        ! tracer does not exist

          ! SET TRACER ID, FIND NAME AND CREATE CHANNEL OBJECT NAME
          CALL get_tracer(status, TRIM(XPOINT(i)%io%trset) &
               , idt, fullname=trfname)
          CALL tracer_halt(substr, status)
          objname = 'emis_'//TRIM(trfname)

          ! CREATE NEW CHANNEL, IF IT DOES NOT ALREADY EXIST
          CALL get_channel_info(status, TRIM(chaname))
          IF (status == 3003) THEN
             CALL info_bi('creating new channel '//TRIM(chaname) ,substr)
             CALL new_channel(status, TRIM(chaname))
          END IF
          CALL channel_halt(substr, status)

          ! CREATE CHANNEL OBJECT, IF IT DOES NOT ALREADY EXIST
          ! (FROM ANOTHER EMISSION POINT!)
          CALL get_channel_object(status, TRIM(chaname), TRIM(objname), &
               p3=XXREF(i)%emis(j)%ptr)
          IF (status == 0) THEN
             CALL info_bi('object '//TRIM(objname)//' in channel '//&
                  &TRIM(chaname)//' exists already', substr)
          ELSE
             CALL info_bi('creating object '//TRIM(objname)//&
                  &' in channel '//TRIM(chaname), substr)
             CALL new_channel_object(status, TRIM(chaname) &
                  , TRIM(objname) &
                  , p3=XXREF(i)%emis(j)%ptr, reprid=reprid)
             CALL channel_halt(substr, status)
             CALL new_attribute(status, TRIM(chaname) &
                  , TRIM(objname) &
                  , 'long_name' &
                  , c='sum of point emissions (only type 3 or 4)'//&
                  &' of tracer '//TRIM(trfname) )
             CALL channel_halt(substr, status)
             CALL new_attribute(status, TRIM(chaname) &
                  , TRIM(objname) &
                  , 'units', c='molec/m^3/s')
             CALL channel_halt(substr, status)
          END IF

       END DO trac_loop3
    END DO point_loop3

    IF (p_parallel_io) THEN
       WRITE(*,*) '------------------------------------------------------'
    END IF

    CALL end_message_bi(modstr, 'COUPLING INITIALIZATION', substr)

  CONTAINS

    ! -------------------------------------------------------------------
    SUBROUTINE init_trac_cpl(XXREF, TRSTR)

      USE messy_main_tracer,        ONLY: R_molarmass
      USE messy_main_constants_mem, ONLY: I8, N_A
      USE messy_main_timer,         ONLY: time_step_len

      IMPLICIT NONE
      INTRINSIC :: REAL

      ! I/O
      ! x referene table
      TYPE(T_RELEASE_POINT_XREF), DIMENSION(:), INTENT(INOUT) :: XXREF
      ! tracer set
      CHARACTER(LEN=*),                         INTENT(IN)    :: TRSTR

      ! LOCAL
      INTEGER :: i, jt
      CHARACTER(LEN=STRLEN_LONG), DIMENSION(:), POINTER :: tracs => NULL()
      INTEGER :: zntrac
      ! correction factor for flux to achieve exact total mass
      REAL(DP)    :: fcorr
      REAL(DP)    :: iets  ! emission time in seconds
      REAL(DP)    :: idts  ! time step length in seconds
      REAL(DP)    :: irest
      INTEGER(I8) :: nsteps

      idts = time_step_len

      point_loop: DO i=1, NPOINT

         IF (TRIM(XPOINT(i)%io%trset) /= TRIM(TRSTR)) CYCLE

         IF (XPOINT(i)%ierr /= 0) CYCLE

         ! PARSE TRACER LIST
         CALL strcrack(XPOINT(i)%io%trlist, ';', tracs, zntrac)

         ! INITIALIZE
         XXREF(i)%ntrac = zntrac
         !
         ALLOCATE(XXREF(i)%idt(zntrac))
         XXREF(i)%idt(:) = 0

         ! CLASSIC ONLY
         IF ((XPOINT(i)%io%type == 1) .OR. (XPOINT(i)%io%type == 3)) THEN
            ALLOCATE(XXREF(i)%flx(zntrac))
            XXREF(i)%flx(:) = 0.0_DP
         END IF
         ALLOCATE(XXREF(i)%ltr(zntrac))
         XXREF(i)%ltr(:) = .FALSE.
         ALLOCATE(XXREF(i)%molarmass(zntrac))
         XXREF(i)%molarmass(:) = 1.0_dp
         !
         trac_loop: DO jt=1, zntrac
            !
            IF (p_parallel_io) THEN
               WRITE(*,*) 'LOOKING FOR TRACER '//TRIM(tracs(jt))//' IN '//&
                    &'TRACER SET '//TRIM(TRSTR)
            END IF
            fullname = TRIM(tracs(jt))
            CALL full2base_sub(status, TRIM(fullname) &
                 , basename, subname)
            CALL tracer_halt(substr, status)

            CALL get_tracer(status, TRSTR, basename       &
                 , subname=subname, idx=XXREF(i)%idt(jt), trinfo=trinfo)
            CALL tracer_halt(substr, status)

            molarmass = trinfo%meta%cask_r(R_molarmass)
            IF (molarmass <= ZERO_EPS) THEN   ! ZERO
               XXREF(i)%molarmass(jt) = 1.0_dp
            ELSE
               XXREF(i)%molarmass(jt) = molarmass
            END IF

            IF ((XPOINT(i)%io%type == 1) .OR. (XPOINT(i)%io%type == 3)) THEN
               IF (molarmass <= ZERO_EPS) THEN   ! ZERO
                  XXREF(i)%flx(jt) = 0.0_DP
               ELSE
                  SELECT CASE(XPOINT(i)%emission_unit)
                  CASE (EU_kg)
                     XXREF(i)%flx(jt) = (XPOINT(i)%io%emission      &
                          * (1.E3_DP/molarmass) / XPOINT(i)%dt)
                     ! FLUX CORRECTION FOR TOTAL MASS
                     iets = XPOINT(i)%dt
                     irest = iets / idts - AINT(iets / idts)
                     IF (irest <= ZERO_EPS ) THEN
                        nsteps = INT(iets / idts, I8)
                        fcorr = 1.0_dp
                     ELSE
                        nsteps = INT(iets / idts, I8) + 1_i8
                        fcorr = XPOINT(i)%dt / &
                             REAL(nsteps * idts, dp)
                     ENDIF
                     XXREF(i)%flx(jt) = XXREF(i)%flx(jt) * fcorr
                  CASE (EU_kgs)
                     XXREF(i)%flx(jt) = XPOINT(i)%io%emission      &
                          * (1.E3_DP/molarmass)
                  CASE (EU_mols)
                     XXREF(i)%flx(jt) = XPOINT(i)%io%emission
                  CASE (EU_molecs)
                     XXREF(i)%flx(jt) = XPOINT(i)%io%emission / N_A
                  END SELECT
               END IF
            END IF

            IF (p_parallel_io) THEN
               SELECT CASE(XPOINT(i)%io%type)
               CASE(1, 3)
                  WRITE(*,*) '     -> ',TRIM(fullname),           &
                       ' (id = ', XXREF(i)%idt(jt),';',            &
                       ' molar mass = ',molarmass,' [g/mol];'
                  WRITE(*,*) '           flux = ',XXREF(i)%flx(jt) &
                       ,'[mol/s])'
                  WRITE(*,*) '           (corrected with factor ',fcorr,';'
                  WRITE(*,*) '            DT_emis = ',iets,' s ;' &
                       , ' DT = ',idts,' s ;'
                  WRITE(*,*) '            #steps = ',nsteps,')'
               CASE(2, 4)
                  WRITE(*,*) '     -> ',TRIM(fullname),             &
                       ' (id = ', XXREF(i)%idt(jt),';',             &
                       ' molar mass = ',molarmass,' [g/mol])'
               END SELECT
            END IF

         END DO trac_loop
         !
         IF (ASSOCIATED(tracs)) THEN
            DEALLOCATE(tracs)
            NULLIFY(tracs)
         END IF
         !
      END DO point_loop

    END SUBROUTINE init_trac_cpl
    ! -------------------------------------------------------------------

    ! -------------------------------------------------------------------
    SUBROUTINE init_emis_cpl(XXREF, XTRAC)

      IMPLICIT NONE

      ! I/O
      TYPE(T_RELEASE_POINT_XREF), DIMENSION(:),   INTENT(INOUT) :: XXREF
      TYPE(T_RELEASED_TRACER_XREF), DIMENSION(:), INTENT(INOUT) :: XTRAC

      ! LOCAL
      INTEGER :: i, jt, j
      INTEGER, DIMENSION(NMAXPOINT)    :: icrp

      tracer_loop: DO jt=1, NTRAC
         !
#ifdef ICON
         IF (TRIM(XTR(jt)%io%trset) /= TRIM(GPTRSTR)) CYCLE
#endif
         IF ((.NOT. L_GP) .AND. (TRIM(XTR(jt)%io%trset) == TRIM(GPTRSTR))) &
              CYCLE
         IF ((.NOT. L_LG) .AND. (TRIM(XTR(jt)%io%trset) == TRIM(LGTRSTR))) &
              CYCLE
         IF ((.NOT. L_CL) .AND. (TRIM(XTR(jt)%io%trset) == TRIM(CLTRSTR))) &
              CYCLE
         !
         icrp(:) = 0
         !
         point_loop: DO i=1, NPOINT
            !
            IF (XPOINT(i)%ierr /= 0) CYCLE

            IF (TRIM(XPOINT(i)%io%trset) /= TRIM(XTR(jt)%io%trset)) CYCLE

            trac_loop: DO j=1, XXREF(i)%ntrac
               !
               IF (XXREF(i)%idt(j) == 0) CYCLE          ! NO ASSOCIATED TRACER
               IF (XXREF(i)%idt(j) == XTRAC(jt)%idt) THEN
                  ! THIS TRACER IS TREXP - TRACER
                  XTRAC(jt)%np = XTRAC(jt)%np + 1
                  icrp(XTRAC(jt)%np) = i
                  XXREF(i)%ltr(j) = .TRUE.
               END IF
               !
            END DO trac_loop
            !
         END DO point_loop

         ! SAVE RESULT
         ALLOCATE(XTRAC(jt)%ixp(XTRAC(jt)%np))
         XTRAC(jt)%ixp(:) = icrp(:XTRAC(jt)%np)

         ! DIAGNOSTIC OUTPUT
         fullname=TRIM(XTR(jt)%io%trname)
         IF (p_parallel_io) THEN
            WRITE(*,*) TRIM(fullname),' -> ',XTRAC(jt)%ixp(:)
         END IF

      END DO tracer_loop

    END SUBROUTINE init_emis_cpl
    ! -------------------------------------------------------------------

  END SUBROUTINE trexp_init_coupling
  ! -------------------------------------------------------------------

  ! -------------------------------------------------------------------
  SUBROUTINE trexp_global_start

    USE messy_main_timer,         ONLY: YEAR, MONTH, DAY, HOUR, MINUTE, SECOND
    USE messy_main_timer,         ONLY: gregor2julian
    USE messy_main_tracer_mem_bi, ONLY: GPTRSTR, LGTRSTR, CLTRSTR
#ifdef _TFCORR
    USE messy_main_constants_mem, ONLY: OneDay
    USE messy_main_timer,         ONLY: time_step_len
#endif

    IMPLICIT NONE

    ! LOCAL
    INTEGER  :: i
    REAL(dp) :: now
#ifdef _TFCORR
    REAL(dp)            :: half_time_step_len_in_days

    half_time_step_len_in_days=time_step_len/2._dp/OneDay
#endif

    ! SET TIME-FLAG FOR EMISSION
    point_loop: DO i=1, NPOINT

       IF (XPOINT(i)%ierr /= 0)  CYCLE

#ifdef ICON
       IF (TRIM(XPOINT(i)%io%trset) /= TRIM(GPTRSTR)) CYCLE
#endif
       IF ((.NOT. L_GP) .AND. (TRIM(XPOINT(i)%io%trset) == TRIM(GPTRSTR))) &
            CYCLE
       IF ((.NOT. L_LG) .AND. (TRIM(XPOINT(i)%io%trset) == TRIM(LGTRSTR))) &
            CYCLE
       IF ((.NOT. L_CL) .AND. (TRIM(XPOINT(i)%io%trset) == TRIM(CLTRSTR))) &
            CYCLE

       ! ONLY FOR CLASSIC POINTS
       IF ((XPOINT(i)%io%type /= 1) .AND. (XPOINT(i)%io%type /= 3)) CYCLE

       now = gregor2julian(YEAR, MONTH, DAY, HOUR, MINUTE, SECOND)
       XPOINT(i)%lnow = (now >= XPOINT(i)%estart) .AND. (now <  XPOINT(i)%estop)

       IF (XPOINT(i)%lnow) THEN
          XPOINT(i)%rnow = 1.0_dp
       ELSE
          XPOINT(i)%rnow = 0.0_dp
       ENDIF

#ifdef _TFCORR
       IF (l_tf_corr) THEN
          !     lfirst     l2nd
          !    |----------.----------|
          !    |--x-------.--]-------|
          XPOINT(i)%lfirst = &
               (now >= XPOINT(i)%estart) .AND. &
               (now < ( XPOINT(i)%estart + half_time_step_len_in_days ))

          ! 2nd time step is shifted by half_time_step_len_in_days
          XPOINT(i)%l2nd = &
               ( (now - half_time_step_len_in_days) >= XPOINT(i)%estart) .AND. &
               ( (now - half_time_step_len_in_days) <  &
               (XPOINT(i)%estart + half_time_step_len_in_days ))
       END IF
#endif

    END DO point_loop

  END SUBROUTINE trexp_global_start
  ! -------------------------------------------------------------------

  ! -------------------------------------------------------------------
  SUBROUTINE trexp_physc

    USE messy_main_mpi_bi,         ONLY: p_pe
#ifndef MESSYTENDENCY
#ifndef ICON
    USE messy_main_data_bi,        ONLY: tm1_3d, tte_3d  &
                                       , qm1 => qm1_3d   &
                                       , qte => qte_3d
#else
    USE messy_main_data_bi,        ONLY: temp_3d, qv_3d
#endif
#endif
    USE messy_main_tracer_mem_bi,  ONLY: ntrac_gp
#ifdef _TFCORR
    USE messy_main_tracer_mem_bi,  ONLY: qxtf, qxt
    USE messy_main_data_bi,        ONLY: eps
#endif
#if ! defined(MESSYTENDENCY) || defined(_TFCORR)
    USE messy_main_tracer_mem_bi,  ONLY: qxtm1, qxtte
#endif
    USE messy_main_grid_def_mem_bi,ONLY: jrow, kproma, nlev
    USE messy_main_grid_def_bi,    ONLY: grmass=>grmassdry &
                                       , zbound=>altitudei_msl &
                                       , grvol
    USE messy_main_data_bi,        ONLY: pmid=>press_3d    &
                                       , pbound=>pressi_3d
    USE messy_main_timer,          ONLY: ztmst=>time_step_len
    USE messy_main_tools,          ONLY: nl_index
    USE messy_main_constants_mem,  ONLY: M_air, N_A, R_gas, M_H2O
    USE messy_main_tracer_mem_bi,  ONLY: GPTRSTR

    IMPLICIT NONE

    INTRINSIC :: EXP, ASSOCIATED, MIN, MAX

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER         :: substr='trexp_physc'
    REAL(dp),                 PARAMETER :: vtmpc1 = M_air / M_H2O - 1._dp
    INTEGER                             :: jt, i, j, k, ii
    INTEGER                             :: jp, jk, jk1, jk2, kmin, kmax
    INTEGER                             :: idt
    ! emission [mol/mol/s]
    REAL(DP), DIMENSION(:,:,:), POINTER :: efld => NULL() ! EMISSION FIELD
    REAL(DP), DIMENSION(:,:),   POINTER :: k0   => NULL() ! REACTION PARTNER
    REAL(DP), DIMENSION(:,:),   POINTER :: temp => NULL() ! temperature [K]
    REAL(DP), DIMENSION(:,:),   POINTER :: cair => NULL() ! 1/cm3
    REAL(DP), DIMENSION(:,:),   POINTER :: sphum => NULL() ! specific humidity
    !                                                != m(H2O)/m(air) [kg/kg]
    REAL(DP)                            :: height ! emission 'height' [Pa], [m]
    REAL(DP)                            :: hmin, hmax, uifc, lifc, frac
    !
    REAL(DP), DIMENSION(kproma,nlev)    :: zxt     ! TRACER AT t=0
    REAL(DP), DIMENSION(kproma,nlev)    :: zxtte   ! TRACER TENDENCY (Q+decay)
    REAL(DP), DIMENSION(kproma,nlev)    :: zloxtte ! LOCAL TRACER TENDENCY
    ! LOCL TRACER TENDENCIES
    REAL(DP), DIMENSION(_RI_X_ZN_(kproma,nlev,ntrac_gp)), TARGET :: zfxtte
!!$#ifdef MESSYTENDENCY
!!$    REAL(DP), DIMENSION(:,:),   POINTER :: szfxtte => NULL() ! TRACER TENDENCY
!!$#endif

    ! INIT

    ALLOCATE(temp(kproma, nlev))
#ifndef MESSYTENDENCY
#ifndef ICON
    temp(:,:) = tm1_3d(_RI_XYZ__(1:kproma,jrow,:)) + &
                tte_3d(_RI_XYZ__(1:kproma,jrow,:))*ztmst
#else
   temp(:,:) = temp_3d(_RI_XYZ__(1:kproma,jrow,:))
#endif
#else
    CALL mtend_get_start_l(mtend_id_t, v0=temp)
#endif

!!#D attila +
#ifdef ECHAM5
    IF (L_LG) THEN
       temp_3d(1:kproma,:,jrow) = temp(1:kproma,:)
    END IF
#endif
!!#D attila -
!!#D clams +
#ifdef ECHAM5
    IF (L_CL) THEN
       temp_3d(1:kproma,:,jrow) = temp(1:kproma,:)
    END IF
#endif
!!#D clams -

    ! ... SPECIFIC HUMIDITY
    ALLOCATE(sphum(kproma,nlev))

#ifndef MESSYTENDENCY
#ifndef ICON
    sphum(1:kproma,:) = MAX( &
         qm1(_RI_XYZ__(1:kproma,jrow,:)) + &
         qte(_RI_XYZ__(1:kproma,jrow,:)) * ztmst,  &
         0._dp)
#else
    sphum(:,:) = MAX(qv_3d(_RI_XYZ__(1:kproma,jrow,:)), 0._dp)
#endif
#else
    CALL mtend_get_start_l(mtend_id_q, v0=sphum)
    sphum(:,:) = MAX(sphum(:,:), 0.0_dp)
#endif

    ! ... CONCENTRATION OF AIR
    ALLOCATE(cair(kproma,SIZE(sphum,2)))
    cair(1:kproma,:)  = (N_A/1.E6_dp) * pmid(_RI_XYZ__(1:kproma,jrow,:)) / &
         (R_gas*temp(1:kproma,:)*(1.0_dp+vtmpc1*sphum(1:kproma,:)))
!!#D attila +
#ifdef ECHAM5
    IF (L_LG) THEN
       cair_3d(1:kproma,:,jrow) = cair(1:kproma,:)
    END IF
#endif
!!#D attila -
!!#D clams +
#ifdef ECHAM5
    IF (L_CL) THEN
       cair_3d(1:kproma,:,jrow) = cair(1:kproma,:)
    END IF
#endif
!!#D clams -

    ! GP representation only
    IF (.NOT. L_GP) THEN
       IF (ASSOCIATED(temp)) THEN
          DEALLOCATE(temp); NULLIFY(temp)
       ENDIF
       IF (ASSOCIATED(sphum)) THEN
          DEALLOCATE(sphum); NULLIFY(sphum)
       ENDIF
       IF (ASSOCIATED(cair)) THEN
          DEALLOCATE(cair); NULLIFY(cair)
       ENDIF
       RETURN
    ENDIF

    ALLOCATE(efld(kproma,nlev, NTRAC))
    efld(:,:,:) = 0.0_DP
    !
    ALLOCATE(k0(kproma, nlev))
    k0(:,:) = 0.0_DP

    ! INITIAISE EMISSION DIAGNOSTIC
    DO i = 1, NPOINT
       IF (TRIM(XPOINT(i)%io%trset) /= TRIM(GPTRSTR)) CYCLE
       DO jt=1, XXREF(i)%ntrac
          IF (ASSOCIATED(XXREF(i)%emis(jt)%ptr)) THEN
             XXREF(i)%emis(jt)%ptr(_RI_XYZ__(:,jrow,:)) = 0.0_dp
          END IF
       END DO
    END DO

    ! 1st: FILL EMISSION FIELD FOR ALL TREXP-TRACERS
    trexp_tracer_loop: DO jt=1, NTRAC
       !
       IF (TRIM(XTR(jt)%io%trset) /= TRIM(GPTRSTR)) CYCLE
       !
       ! SET TRACER ID
       idt = XTRAC(jt)%idt
       !
       point_loop: DO i=1, XTRAC(jt)%np
          !
          ii = XTRAC(jt)%ixp(i)
          !
          IF (TRIM(XPOINT(ii)%io%trset) /= TRIM(GPTRSTR)) CYCLE
          IF (XPOINT(ii)%ierr /= 0)    CYCLE
          !
          IF (XPOINT(ii)%pe   /= p_pe) CYCLE   ! LOCATION IS NOT ON THIS pe
          IF (XPOINT(ii)%jrow /= jrow) CYCLE   ! ... IS NOT IN CURRENT ROW
          ! GET VECTOR INDEX
          jp = XPOINT(ii)%jp

          SELECT CASE(XPOINT(ii)%io%type)
          CASE (1, 3)
             !
             IF (.NOT. XPOINT(ii)%lnow)   CYCLE   ! ... IS NOT ACTIVE NOW
             !
             ! GET LEVEL INDEX
             SELECT CASE (XPOINT(ii)%height_unit)
             CASE(HU_hPa, HU_Pa)
                ! Pressure coordinates, scale if necessary (hPa -> Pa)
                height = XPOINT(ii)%io%height * XPOINT(ii)%height_scal
                CALL nl_index(pbound(_RI_XYZ__(jp,jrow,1:nlev)), height, jk)
             CASE(HU_mASL)
                height = XPOINT(ii)%io%height  ! m ASL !!!
                CALL nl_index(zbound(_RI_XYZ__(jp,jrow,1:nlev)), height, jk)
             CASE(HU_mAGL)
                height = zbound(_RI_XYZ__(jp,jrow,nlev+1))                    &
                   &   + XPOINT(ii)%io%height ! m AGL -> m ASL !!!
                CALL nl_index(zbound(_RI_XYZ__(jp,jrow,1:nlev)), height, jk)
             END SELECT

             trac_loop_cl: DO j=1, XXREF(ii)%ntrac
                !
                ! RELEASED TRACER IS NON-TREXP
                IF (.NOT. XXREF(ii)%ltr(j))  CYCLE
                ! ... IS NOT CURRENT TRACER
                IF (XXREF(ii)%idt(j) /= idt) CYCLE
                !
                ! (mol/s) / (kg/(kg/mol)) -> mol/mol/s
                efld(jp, jk, jt) = efld(jp, jk, jt) + XXREF(ii)%flx(j)  &
                     / ( grmass(_RI_XYZ__(jp,jrow,jk)) * ( 1.E3_DP/M_air ) )
                !
                ! (mol/s) * Na / m^3 -> molec/m^3/s
                IF (XPOINT(ii)%io%type == 3) THEN
                   XXREF(ii)%emis(j)%ptr(_RI_XYZ__(jp,jrow,jk)) = &
                        XXREF(ii)%emis(j)%ptr(_RI_XYZ__(jp,jrow,jk)) + &
                        XXREF(ii)%flx(j) * N_A / grvol(_RI_XYZ__(jp,jrow,jk))
                ENDIF
                !
             END DO trac_loop_cl

          CASE (2, 4)
             height_loop: DO k=1, SIZE(XPOINT(ii)%height, 1)-1
                ! GET LEVEL INDICES
                SELECT CASE (XPOINT(ii)%height_unit)
                CASE(HU_Pa, HU_hPa)
                   hmin = MIN(XPOINT(ii)%height(k), XPOINT(ii)%height(k+1)) &
                      &   * XPOINT(ii)%height_scal
                   hmax = MAX(XPOINT(ii)%height(k), XPOINT(ii)%height(k+1)) &
                      &   * XPOINT(ii)%height_scal
                   CALL nl_index(pbound(_RI_XYZ__(jp,jrow,1:nlev)), hmin, jk1)
                   CALL nl_index(pbound(_RI_XYZ__(jp,jrow,1:nlev)), hmax, jk2)
                CASE(HU_mASL)
                   hmin = MIN(XPOINT(ii)%height(k), XPOINT(ii)%height(k+1))
                   hmax = MAX(XPOINT(ii)%height(k), XPOINT(ii)%height(k+1))
                   CALL nl_index(zbound(_RI_XYZ__(jp,jrow,1:nlev)), hmin, jk1)
                   CALL nl_index(zbound(_RI_XYZ__(jp,jrow,1:nlev)), hmax, jk2)
                CASE(HU_mAGL)
                   hmin = MIN(XPOINT(ii)%height(k), XPOINT(ii)%height(k+1)) &
                      &   + zbound(_RI_XYZ__(jp,jrow,nlev+1))
                   hmax = MAX(XPOINT(ii)%height(k), XPOINT(ii)%height(k+1)) &
                      &   + zbound(_RI_XYZ__(jp,jrow,nlev+1))
                   CALL nl_index(zbound(_RI_XYZ__(jp,jrow,1:nlev)), hmin, jk1)
                   CALL nl_index(zbound(_RI_XYZ__(jp,jrow,1:nlev)), hmax, jk2)
                END SELECT
                kmin = MIN(jk1, jk2)
                kmax = MAX(jk1, jk2)

                trac_loop_ts: DO j=1, XXREF(ii)%ntrac
                   !
                   ! RELEASED TRACER IS NON-TREXP
                   IF (.NOT. XXREF(ii)%ltr(j))  CYCLE
                   ! ... IS NOT CURRENT TRACER
                   IF (XXREF(ii)%idt(j) /= idt) CYCLE
                   !
                   DO jk=kmin, kmax
                      SELECT CASE (XPOINT(ii)%height_unit)
                      CASE (HU_hPa, HU_Pa)
                         lifc = MIN(pbound(_RI_XYZ__(jp,jrow,jk)) &
                            &     , pbound(_RI_XYZ__(jp,jrow,jk+1)))
                         uifc = MAX(pbound(_RI_XYZ__(jp,jrow,jk)) &
                            &     , pbound(_RI_XYZ__(jp,jrow,jk+1)))
                      CASE (HU_mASL, HU_mAGL)
                         lifc = MIN(zbound(_RI_XYZ__(jp,jrow,jk)) &
                            &     , zbound(_RI_XYZ__(jp,jrow,jk+1)))
                         uifc = MAX(zbound(_RI_XYZ__(jp,jrow,jk)) &
                            &     , zbound(_RI_XYZ__(jp,jrow,jk+1)))
                      END SELECT
                      frac = (MIN(uifc, hmax) - MAX(lifc, hmin)) / (hmax - hmin)
                      SELECT CASE (XPOINT(ii)%emission_unit)
                      CASE (EU_kgs)
                         ! (kg/s) / (mol/mol) / kg -> mol/mol/s
                         efld(jp, jk, jt) = efld(jp, jk, jt)      &
                            + frac * XPOINT(ii)%flux(k)           &
                            * (M_air / XXREF(ii)%molarmass(j))    &
                            / grmass(_RI_XYZ__(jp,jrow,jk))
                         ! (kg/s) * (mol/kg) * Na / m^3 -> molec/m^3/s
                         IF (XPOINT(ii)%io%type == 4) THEN
                            XXREF(ii)%emis(j)%ptr(_RI_XYZ__(jp,jrow,jk)) =    &
                                 XXREF(ii)%emis(j)%ptr(_RI_XYZ__(jp,jrow,jk)) &
                                 + frac * XPOINT(ii)%flux(k)                  &
                                 * 1.E3_DP / XXREF(ii)%molarmass(j) * N_A     &
                                 / grvol(_RI_XYZ__(jp,jrow,jk))
                         ENDIF
                      CASE (EU_mols)
                         ! (mol/s) / (kg/(kg/mol) -> mol/mol/s
                         efld(jp, jk, jt) = efld(jp, jk, jt)                  &
                                          + frac * XPOINT(ii)%flux(k)         &
                                          / ( grmass(_RI_XYZ__(jp,jrow,jk))   &
                                          * ( 1.E3_DP/M_air ) )
                         ! (mol/s) * Na / m^3 -> molec/m^3/s
                         IF (XPOINT(ii)%io%type == 4) THEN
                            XXREF(ii)%emis(j)%ptr(_RI_XYZ__(jp,jrow,jk)) =    &
                                 XXREF(ii)%emis(j)%ptr(_RI_XYZ__(jp,jrow,jk)) &
                                 + frac * XPOINT(ii)%flux(k) * N_A            &
                                 / grvol(_RI_XYZ__(jp,jrow,jk))
                         ENDIF
                      CASE (EU_molecs)
                         ! (molec/s) / Na / (kg/(kg/mol) -> mol/mol/s
                         efld(jp, jk, jt) = efld(jp, jk, jt)                  &
                                      + frac * XPOINT(ii)%flux(k)             &
                                      / N_A / ( grmass(_RI_XYZ__(jp,jrow,jk)) &
                                                * ( 1.E3_DP/M_air ) )
                         ! (molec/s) / m^3 -> molec/m^3/s
                         IF (XPOINT(ii)%io%type == 4) THEN
                            XXREF(ii)%emis(j)%ptr(_RI_XYZ__(jp,jrow,jk)) =    &
                                 XXREF(ii)%emis(j)%ptr(_RI_XYZ__(jp,jrow,jk)) &
                                 + frac * XPOINT(ii)%flux(k)                  &
                                 / grvol(_RI_XYZ__(jp,jrow,jk))
                         ENDIF
                      END SELECT
                   END DO
                END DO trac_loop_ts
             END DO height_loop

          END SELECT
             !
       END DO point_loop
       !
    END DO trexp_tracer_loop

    ! 2nd: TRACER TENDENCY: EMISSION + DECAY / REACTION
    trexp_tracer_loop2: DO jt=1, NTRAC
       !
       IF (TRIM(XTR(jt)%io%trset) /= TRIM(GPTRSTR)) CYCLE
       !
       zloxtte(:,:) = 0.0_dp
       !
       ! SET TRACER ID
       idt = XTRAC(jt)%idt
       !
       ! REACTION PARTNER
       ! ... start value
       SELECT CASE(XTR(jt)%io%order)
       CASE(0)
          ! DECAY
          k0(:,:) = XTR(jt)%io%ka
          !
        CASE(1,-1,-2,-3)
          ! REACTION PARTNER
          k0(:,:) = XTR(jt)%rp(_RI_XYZ__(1:kproma,jrow,:))
          IF (ASSOCIATED(XTR(jt)%rpte)) THEN
             ! TRACER + TENDENCY
             k0(:,:) = k0(:,:) + XTR(jt)%rpte(_RI_XYZ__(1:kproma,jrow,:))*ztmst
          END IF
          ! ... unit conversion -> cm^(-3)
          IF (XTR(jt)%lv2c) THEN
             CALL vmr2conc(k0(:,:), pmid(_RI_XYZ__(1:kproma,jrow,:)), temp(:,:))
          ENDIF
          !
          SELECT CASE(XTR(jt)%io%order)
          CASE(1)
             ! -> 1/s
             k0(:,:) = k0(:,:) * XTR(jt)%io%ka &
                  * exp(- XTR(jt)%io%Ta / temp(:,:))
          CASE(-1)
             ! -> 1/s
             k0(:,:) = k0(:,:) * (XTR(jt)%io%ka + cair(:,:) * XTR(jt)%io%Ta)
          CASE(-2)
             ! special case CH4 + OH -> 1/s
             k0(:,:) = k0(:,:) * ( &
                  1.85E-20_dp * exp(2.82_dp*log(temp(:,:)) - 987._dp/temp(:,:))&
                  )
          CASE(-3)
             ! special case SO2 + OH
             ! reaction rate from mecca.eqn
             k0(:,:) = k0(:,:) * k_3rd(temp(:,:),cair(:,:),3.3E-31_dp, &
                       4.3_dp,1.6E-12_dp,0._dp,0.6_dp)
          END SELECT
          !
       CASE DEFAULT
          !
          ! SHOULD NEVER BE REACHED
          !
       END SELECT
       !
       ! START VALUE
#ifndef MESSYTENDENCY
       zxt(:,:) = qxtm1(_RI_X_ZN_(1:kproma,:,idt)) &
            + qxtte(_RI_X_ZN_(1:kproma,:,idt)) * ztmst
#else
       CALL mtend_get_start_l(idt, v0=zxt)
#endif
       !
       ! INTEGRATE
       CALL solve(zxtte(:,:), zxt(:,:), efld(:,:,jt), k0(:,:), ztmst)
       !
#ifdef _TFCORR
       ! The "2nd emission step" correction needs to be applied WITHOUT
       ! the current emission tendency, i.e., before qxtte = qxtte + ...
       ! The "2nd emission step" correction can only be applied once
       ! per tracer ... (e.g. for emissions at multiple points into one tracer).
       ! Note: The correction is applied at ANY "2nd emssion step".
       IF (l_tf_corr) THEN
          point_loop_2nd: DO i=1, XTRAC(jt)%np
             !
             ii = XTRAC(jt)%ixp(i)
             !
             IF (XPOINT(ii)%ierr /= 0)      CYCLE
             IF (XPOINT(ii)%pe   /= p_pe)   CYCLE  ! LOCATION IS NOT ON THIS pe
             IF (XPOINT(ii)%jrow /= jrow)   CYCLE  ! ... IS NOT IN CURRENT ROW
             IF (.NOT. XPOINT(ii)%lnow)     CYCLE  ! ... IS NOT ACTIVE NOW
             !
             ! ONLY FOR CLASSIC POINTS
             IF ((XPOINT(ii)%io%type /= 1) .AND. (XPOINT(ii)%io%type /= 3)) &
                  CYCLE
             !
             jp = XPOINT(ii)%jp
             SELECT CASE (XPOINT(ii)%height_unit)
             CASE(HU_hPa , HU_Pa)
                height = XPOINT(ii)%io%height * XPOINT(ii)%height_scal
                CALL nl_index(pbound(_RI_XYZ__(jp,jrow,1:nlev)), height, jk)
             CASE(HU_mASL)
                height = XPOINT(ii)%io%height  ! m ASL !!!
                CALL nl_index(zbound(_RI_XYZ__(jp,jrow,1:nlev)), height, jk)
             CASE(HU_mAGL)
                height = zbound(_RI_XYZ__(jp,jrow,nlev+1)) &
                   &   + XPOINT(ii)%io%height ! m AGL -> m ASL !!!
                CALL nl_index(zbound(_RI_XYZ__(jp,jrow,1:nlev)), height, jk)
             END SELECT
             !
             IF (XPOINT(ii)%l2nd) THEN
                zloxtte(jp,jk) = - eps * qxt(_RI_X_ZN_(jp,jk,idt)) / ztmst
             END IF
          END DO point_loop_2nd
       END IF
#endif
       !
       ! ADD TENDENCY
       zloxtte(:,:) = zloxtte(:,:) + zxtte(:,:)
#ifndef MESSYTENDENCY
       qxtte(_RI_X_ZN_(1:kproma,:,idt)) = qxtte(_RI_X_ZN_(1:kproma,:,idt)) &
            + zloxtte(:,:)
#else
       CALL mtend_add_l(my_handle, idt, px=zloxtte)
#endif
       !
#ifdef _TFCORR
       ! The "1st emission step" correction needs to be applied WITH
       ! the current emission tendency, i.e., after qxtte = qxtte + ...
       ! The "1st emission step" correction is only be applied once
       ! per tracer ... (e.g. for emissions at multiple points into one tracer).
       ! Note: The correction is applied at ANY "1st emssion step".
       IF (l_tf_corr) THEN
          point_loop_1st: DO i=1, XTRAC(jt)%np
             !
             ii = XTRAC(jt)%ixp(i)
             !
             IF (XPOINT(ii)%ierr /= 0)    CYCLE
             IF (XPOINT(ii)%pe   /= p_pe) CYCLE    ! LOCATION IS NOT ON THIS pe
             IF (XPOINT(ii)%jrow /= jrow) CYCLE    ! ... IS NOT IN CURRENT ROW
             IF (.NOT. XPOINT(ii)%lnow)   CYCLE    ! ... IS NOT ACTIVE NOW
             !
             ! ONLY FOR CLASSIC POINTS
             IF ((XPOINT(ii)%io%type /= 1) .AND. (XPOINT(ii)%io%type /= 3)) CYCLE
             !
             jp = XPOINT(ii)%jp
             SELECT CASE (XPOINT(ii)%height_unit)
             CASE(HU_hPa, HU_Pa)
                height = XPOINT(ii)%io%height * XPOINT(ii)%height_scal
                CALL nl_index(pbound(_RI_XYZ__(jp,jrow,1:nlev)), height, jk)
             CASE(HU_mASL)
                height = XPOINT(ii)%io%height  ! m ASL !!!
                CALL nl_index(zbound(_RI_XYZ__(jp,jrow,1:nlev)), height, jk)
             CASE(HU_mAGL)
                height = zbound(_RI_XYZ__(jp,jrow,nlev+1)) &
                   &   + XPOINT(ii)%io%height ! m AGL -> m ASL !!!
                CALL nl_index(zbound(_RI_XYZ__(jp,jrow,1:nlev)), height, jk)
             END SELECT
             !
             IF (XPOINT(ii)%lfirst) THEN
                qxtf(_RI_X_ZN_(jp,jk,idt)) = qxtm1(_RI_X_ZN_(jp,jk,idt)) + &
                     qxtte(_RI_X_ZN_(jp,jk,idt)) * ztmst
             END IF
          END DO point_loop_1st
       END IF
#endif
       !
    END DO trexp_tracer_loop2

    ! 3rd: IF REQUIRED: ADD EMISSION TENDENCY FOR ALL NON-TREXP-TRACERS
    force: IF (l_force_emis) THEN
       !
       zfxtte(:,:,:) = 0.0_dp
       !
       point_loop2: DO i=1, NPOINT
          !
          IF (TRIM(XPOINT(i)%io%trset) /= TRIM(GPTRSTR)) CYCLE
          IF (XPOINT(i)%ierr /= 0)     CYCLE
          IF (XPOINT(i)%pe   /= p_pe)  CYCLE   ! LOCATION IS NOT ON THIS pe
          IF (XPOINT(i)%jrow /= jrow)  CYCLE   ! LOCATION IS NOT IN CURRENT ROW
          ! GET VECTOR INDEX
          jp    = XPOINT(i)%jp

          SELECT CASE(XPOINT(i)%io%type)

          CASE (1, 3)
             IF (.NOT. XPOINT(i)%lnow)    CYCLE   ! EMISSION NOT AT THIS TIME
             !
             ! GET LEVEL INDEX
             SELECT CASE (XPOINT(i)%height_unit)
             CASE(HU_hPa, HU_Pa)
                height = XPOINT(i)%io%height * XPOINT(i)%height_scal
                CALL nl_index(pbound(_RI_XYZ__(jp,jrow,1:nlev)), height, jk)
             CASE(HU_mASL)
                height = XPOINT(i)%io%height  ! m ASL !!!
                CALL nl_index(zbound(_RI_XYZ__(jp,jrow,1:nlev)), height, jk)
             CASE(HU_mAGL)
                height = zbound(_RI_XYZ__(jp,jrow,nlev+1)) &
                   &   + XPOINT(i)%io%height ! m AGL -> m ASL !!!
                CALL nl_index(zbound(_RI_XYZ__(jp,jrow,1:nlev)), height, jk)
             END SELECT
             !
             trac_loop_cl2: DO j=1, XXREF(i)%ntrac
                !
                IF (XXREF(i)%ltr(j)) CYCLE  ! TRACER IS TREXP-TRACER
                !
                ! GET TRACER ID
                idt   = XXREF(i)%idt(j)
                IF (idt == 0) CYCLE        ! tracer does not exist
                !
#ifdef _TFCORR
                IF (l_tf_corr) THEN
                   IF (XPOINT(i)%l2nd) THEN
                      zfxtte(_RI_X_ZN_(jp,jk,idt)) = &
                           - eps * qxt(_RI_X_ZN_(jp,jk,idt)) / ztmst
                   END IF
                END IF
#endif
                ! ADD EMISSION TO TRACER TENDENCY
                ! (mol/s) / (kg/(kg/mol)) -> mol/mol/s
                zfxtte(_RI_X_ZN_(jp,jk,idt)) = zfxtte(_RI_X_ZN_(jp,jk,idt)) &
                     + XXREF(i)%flx(j)                                      &
                     / ( grmass(_RI_XYZ__(jp,jrow,jk)) * ( 1.E3_DP/M_air ) )
                ! (mol/s) * Na / m^3 -> molec/m^3/s
                IF (XPOINT(i)%io%type == 3) THEN
                   XXREF(i)%emis(j)%ptr(_RI_XYZ__(jp,jrow,jk)) =              &
                        XXREF(i)%emis(j)%ptr(_RI_XYZ__(jp,jrow,jk))           &
                        + XXREF(i)%flx(j) * N_A / grvol(_RI_XYZ__(jp,jrow,jk))
                END IF
                !
#ifdef _TFCORR
                IF (l_tf_corr) THEN
                   IF (XPOINT(i)%lfirst) THEN
                      qxtf(_RI_X_ZN_(jp,jk,idt)) = &
                           qxtm1(_RI_X_ZN_(jp,jk,idt)) +  &
                           ( qxtte(_RI_X_ZN_(jp,jk,idt)) + &
                             zfxtte(_RI_X_ZN_(jp,jk,idt)) ) * ztmst
                   END IF
                END IF
#endif
                !
             END DO trac_loop_cl2

          CASE(2, 4)
             height_loop2: DO k=1, SIZE(XPOINT(i)%height, 1)-1
                ! GET LEVEL INDICES
                SELECT CASE (XPOINT(i)%height_unit)
                CASE(HU_Pa, HU_hPa)
                   hmin = MIN(XPOINT(i)%height(k), XPOINT(i)%height(k+1))   &
                      &   * XPOINT(i)%height_scal
                   hmax = MAX(XPOINT(i)%height(k), XPOINT(i)%height(k+1))   &
                      &   * XPOINT(i)%height_scal
                   CALL nl_index(pbound(_RI_XYZ__(jp,jrow,1:nlev)), hmin, jk1)
                   CALL nl_index(pbound(_RI_XYZ__(jp,jrow,1:nlev)), hmax, jk2)
                CASE(HU_mASL)
                   hmin = MIN(XPOINT(i)%height(k), XPOINT(i)%height(k+1))
                   hmax = MAX(XPOINT(i)%height(k), XPOINT(i)%height(k+1))
                   CALL nl_index(zbound(_RI_XYZ__(jp,jrow,1:nlev)), hmin, jk1)
                   CALL nl_index(zbound(_RI_XYZ__(jp,jrow,1:nlev)), hmax, jk2)
                CASE(HU_mAGL)
                   hmin = MIN(XPOINT(i)%height(k), XPOINT(i)%height(k+1))   &
                      &   + zbound(_RI_XYZ__(jp,jrow,nlev+1))
                   hmax = MAX(XPOINT(i)%height(k), XPOINT(i)%height(k+1))   &
                      &   + zbound(_RI_XYZ__(jp,jrow,nlev+1))
                   CALL nl_index(zbound(_RI_XYZ__(jp,jrow,1:nlev)), hmin, jk1)
                   CALL nl_index(zbound(_RI_XYZ__(jp,jrow,1:nlev)), hmax, jk2)
                END SELECT
                kmin = MIN(jk1, jk2)
                kmax = MAX(jk1, jk2)

                trac_loop_ts2: DO j=1, XXREF(i)%ntrac
                   !
                   IF (XXREF(i)%ltr(j)) CYCLE  ! TRACER IS TREXP-TRACER
                   !
                   ! GET TRACER ID
                   idt   = XXREF(i)%idt(j)
                   IF (idt == 0) CYCLE        ! TRACER DOES NOT EXIST
                   !
                   ! ADD EMISSION TO TRACER TENDENCY
                   DO jk=kmin, kmax
                      SELECT CASE (XPOINT(i)%height_unit)
                      CASE (HU_hPa, HU_Pa)
                         lifc = MIN(pbound(_RI_XYZ__(jp,jrow,jk)) &
                            &     , pbound(_RI_XYZ__(jp,jrow,jk+1)))
                         uifc = MAX(pbound(_RI_XYZ__(jp,jrow,jk)) &
                            &     , pbound(_RI_XYZ__(jp,jrow,jk+1)))
                      CASE (HU_mASL, HU_mAGL)
                         lifc = MIN(zbound(_RI_XYZ__(jp,jrow,jk)) &
                            &     , zbound(_RI_XYZ__(jp,jrow,jk+1)))
                         uifc = MAX(zbound(_RI_XYZ__(jp,jrow,jk)) &
                            &     , zbound(_RI_XYZ__(jp,jrow,jk+1)))
                      END SELECT
                      frac = (MIN(uifc, hmax) - MAX(lifc, hmin)) / (hmax - hmin)
                      SELECT CASE (XPOINT(i)%emission_unit)
                      CASE (EU_kgs)
                         ! (kg/s) / ((kg/mol)/(kg/mol)) -> mol/mol/s
                         zfxtte(_RI_X_ZN_(jp,jk,idt)) =          &
                              zfxtte(_RI_X_ZN_(jp,jk,idt))       &
                              + frac * XPOINT(i)%flux(k)         &
                              * (M_air / XXREF(ii)%molarmass(j)) &
                              / grmass(_RI_XYZ__(jp,jrow,jk))
                         ! (kg/s) / (kg/mol) * Na / m^3 -> molec/m^3/s
                         IF (XPOINT(i)%io%type == 3) THEN
                            XXREF(i)%emis(j)%ptr(_RI_XYZ__(jp,jrow,jk)) =     &
                                 XXREF(i)%emis(j)%ptr(_RI_XYZ__(jp,jrow,jk))  &
                                 + frac * XPOINT(i)%flux(k)                   &
                                 * 1.3_DP / XXREF(ii)%molarmass(j) * N_A      &
                                 / grvol(_RI_XYZ__(jp,jrow,jk))
                         ENDIF
                      CASE (EU_mols)
                         ! (mol/s) / (kg/(kg/mol)) -> mol/mol/s
                         zfxtte(_RI_X_ZN_(jp,jk,idt)) =         &
                              zfxtte(_RI_X_ZN_(jp,jk,idt))      &
                              + frac * XPOINT(i)%flux(k)        &
                              / ( grmass(_RI_XYZ__(jp,jrow,jk)) &
                              * ( 1.E3_DP/M_air ) )
                         ! (mol/s) * Na / m^3 -> molec/m^3/s
                         IF (XPOINT(i)%io%type == 3) THEN
                            XXREF(i)%emis(j)%ptr(_RI_XYZ__(jp,jrow,jk)) =    &
                                 XXREF(i)%emis(j)%ptr(_RI_XYZ__(jp,jrow,jk)) &
                                 + frac * XPOINT(i)%flux(k) * N_A            &
                                 / grvol(_RI_XYZ__(jp,jrow,jk))
                         END IF
                      CASE (EU_molecs)
                         ! (molec/s) / Na / (kg/(kg/mol)) -> mol/mol/s
                         zfxtte(_RI_X_ZN_(jp,jk,idt)) =         &
                              zfxtte(_RI_X_ZN_(jp,jk,idt))      &
                              + frac * XPOINT(i)%flux(k) / N_A  &
                              / ( grmass(_RI_XYZ__(jp,jrow,jk)) &
                              * ( 1.E3_DP/M_air ) )
                         ! (molec/s) / m^3 -> molec/m^3/s
                         IF (XPOINT(i)%io%type == 3) THEN
                            XXREF(i)%emis(j)%ptr(_RI_XYZ__(jp,jrow,jk)) =     &
                                 XXREF(i)%emis(j)%ptr(_RI_XYZ__(jp,jrow,jk))  &
                                 + frac * XPOINT(i)%flux(k)                   &
                                 / grvol(_RI_XYZ__(jp,jrow,jk))
                         END IF
                      END SELECT
                   END DO
                END DO trac_loop_ts2
             END DO height_loop2
          END SELECT
          !
       END DO point_loop2
       !
#ifndef MESSYTENDENCY
       qxtte(_RI_X_ZN_(1:kproma,:,:)) = &
            qxtte(_RI_X_ZN_(1:kproma,:,:)) + &
            zfxtte(_RI_X_ZN_(1:kproma,:,:))
#else
       CALL mtend_add_l(my_handle, mtend_id_tracer, pxt=zfxtte)
#endif
!!$       DO idt=1, ntrac_gp
!!$#ifndef MESSYTENDENCY
!!$          qxtte(_RI_X_ZN_(1:kproma,:,idt)) = &
!!$               qxtte(_RI_X_ZN_(1:kproma,:,idt)) + &
!!$               zfxtte(_RI_X_ZN_(1:kproma,:,idt))
!!$#else
!!$          szfxtte => zfxtte(_RI_X_ZN_(1:kproma,:,idt))
!!$          CALL mtend_add_l(my_handle, idt, px=szfxtte)
!!$#endif
!!$       END DO
!!$#ifdef MESSYTENDENCY
!!$       NULLIFY(szfxtte)
!!$#endif
       !
    END IF force

    ! CLEAN UP
    DEALLOCATE(efld)
    NULLIFY(efld)
    DEALLOCATE(k0)
    NULLIFY(k0)
    DEALLOCATE(temp)
    NULLIFY(temp)
    DEALLOCATE(sphum)
    NULLIFY(sphum)
    DEALLOCATE(cair)
    NULLIFY(cair)

  END SUBROUTINE trexp_physc
  ! -------------------------------------------------------------------

  ! -------------------------------------------------------------------
  SUBROUTINE trexp_global_end

#ifdef ECHAM5

!!#D clams +
    USE messy_main_tracer_mem_bi,  ONLY: qxt_c, qxtte_c
    USE messy_clams_tools_e5,      ONLY: gp2cl_e5
    USE messy_main_tracer_mem_bi,  ONLY: dnparts, dnparts_max
    USE messy_clams_global,        ONLY: cmcell
    USE messy_main_tracer_mem_bi,  ONLY: CLTRSTR
!!#D clams -

!!#D attila +
    USE messy_main_tracer_mem_bi,  ONLY: qxtm1_a, qxtte_a
    USE messy_main_tracer_mem_bi,  ONLY: NCELL
    USE messy_attila_tools_e5,     ONLY: gp2lg_e5, GNCB
    USE messy_main_tracer_mem_bi,  ONLY: LGTRSTR
#ifdef _TFCORR
    USE messy_main_tracer_mem_bi,  ONLY: qxtf_a , qxt_a
    USE messy_main_data_bi,        ONLY: eps
#endif
!!#D attila -

    USE messy_main_blather_bi,     ONLY: error_bi
    USE messy_main_transform_bi,   ONLY: trp_gpdc_gpgl
    USE messy_main_grid_def_mem_bi, ONLY: nlev, nproma, npromz, ngpblks
    USE messy_main_grid_def_bi,     ONLY: grvol
    USE messy_main_data_bi,        ONLY: pmid=>press_3d                &
                                       , pbound=>pressi_3d
    USE messy_main_timer,          ONLY: ztmst=>time_step_len
    USE messy_main_tools,          ONLY: nl_index
    USE messy_main_constants_mem,  ONLY: M_air, N_A

    IMPLICIT NONE

    INTRINSIC :: EXP, ASSOCIATED, SUM, NINT, MIN, MAX

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER         :: substr='trexp_global_end'
    INTEGER                             :: jt, i, j, k, ii, jc, nc
    INTEGER                             :: jgx, jgy, jk
    INTEGER                             :: jk1, jk2, kmin, kmax
    INTEGER                             :: jgxlg, jgylg, jklg
    INTEGER                             :: idt
    REAL(DP), DIMENSION(:,:,:), POINTER :: press_gl => NULL() ! pressure [Pa]
    REAL(DP)                            :: hmin, hmax, uifc, lifc, frac
    ! emission [mol/mol/s]
    ! ATTILA
    REAL(DP), DIMENSION(:,:),   POINTER :: efld_lg => NULL()  ! EMISSION FIELD
    REAL(DP), DIMENSION(:,:,:), POINTER :: k0_gp => NULL() ! REACT. PARTNER
    REAL(DP), DIMENSION(:),     POINTER :: k0_lg => NULL() ! REACT. PARTNER
    REAL(DP), DIMENSION(:),     POINTER :: zxt  => NULL() ! TRACER AT t=0
    REAL(DP), DIMENSION(:),     POINTER :: zxtte=> NULL() ! TRACER TENDENCY
    ! CLaMS
    REAL(DP), DIMENSION(:,:),   POINTER :: efld_cl => NULL()  ! EMISSION FIELD
    REAL(DP), DIMENSION(:),     POINTER :: k0_cl => NULL() ! REACT. PARTNER
    REAL(DP), DIMENSION(:),     POINTER :: zxt_c  => NULL() ! TRACER AT t=0
    REAL(DP), DIMENSION(:),     POINTER :: zxtte_c=> NULL() ! TRACER TENDENCY
    REAL(DP)                            :: height ! emission 'height' [Pa]
    REAL(DP), DIMENSION(:,:,:), POINTER :: grvol_gl => NULL() ! volume [m^3]

    ! LG representation only
    IF (.NOT. L_LG .AND. .NOT. L_CL) RETURN

!!#D attila +
    atilla: IF (L_LG) THEN
       ! INIT
       ! set rest of last row to 1 for division below ...
       if (npromz .lt. nproma) then
          temp_3d(npromz+1:,:,ngpblks) = 1.0_DP
       endif
       ! globalize pressure (for position information) ...
       CALL trp_gpdc_gpgl(1, pbound, press_gl)
       CALL trp_gpdc_gpgl(1, grvol, grvol_gl)
       !
       ALLOCATE(efld_lg(NCELL, NTRAC))
       efld_lg(:,:) = 0.0_DP
       !
       ALLOCATE(k0_lg(NCELL))
       k0_lg(:) = 0.0_DP
       !
       ALLOCATE(zxt(NCELL))
       ALLOCATE(zxtte(NCELL))

       ! INITIALISE EMISSION DIAGNOSTIC
       DO i=1, NPOINT
          IF (TRIM(XPOINT(i)%io%trset) /= TRIM(LGTRSTR)) CYCLE
          DO jt=1, XXREF(i)%ntrac
             IF (ASSOCIATED(XXREF(i)%emis(jt)%ptr)) THEN
                XXREF(i)%emis(jt)%ptr(:,:,:) = 0.0_dp
             END IF
          END DO
       END DO

       ! 1st: FILL EMISSION FIELD FOR ALL TREXP-TRACERS
       tracer_loop: DO jt=1, NTRAC
          !
          IF (TRIM(XTR(jt)%io%trset) /= TRIM(LGTRSTR)) CYCLE
          !
          ! SET TRACER ID
          idt = XTRAC(jt)%idt
          !
          point_loop: DO i=1, XTRAC(jt)%np
             !
             ii = XTRAC(jt)%ixp(i)
             IF (XPOINT(ii)%ierr /= 0)    CYCLE
             !
             ! GET INDICES IN GLOBAL FIELD
             ! ... lon, lat
             jgx = XPOINT(ii)%jgx  ! lon index
             jgy = XPOINT(ii)%jgy  ! lat index

             SELECT CASE(XPOINT(ii)%io%type)
             CASE (1, 3)
                IF (.NOT. XPOINT(ii)%lnow)   CYCLE   ! ... IS NOT ACTIVE NOW
                ! ... lev
                SELECT CASE (XPOINT(ii)%height_unit)
                CASE (HU_hPa, HU_Pa)
                   height = XPOINT(ii)%io%height * XPOINT(ii)%height_scal
                CASE DEFAULT
                   CALL error_bi('Emission height in metre not implemented for ATTILA!',substr)
                END SELECT
                CALL nl_index(press_gl(jgx,1:nlev,jgy), height, jk)

                trac_loop_cl: DO j=1, XXREF(ii)%ntrac
                   !
                   ! RELEASED TRACER IS NON-TREXP
                   IF (.NOT. XXREF(ii)%ltr(j))  CYCLE
                   ! ... IS NOT CURRENT TRACER
                   IF (XXREF(ii)%idt(j) /= idt) CYCLE
                   !
                   DO jc=1, NCELL
                      jgxlg = NINT(iplon(jc))
                      jgylg = NINT(iplat(jc))
                      jklg  = NINT(iplev(jc))
                      IF (jgxlg /= jgx) CYCLE   ! wrong longitude
                      IF (jgylg /= jgy) CYCLE   ! wrong latitude
                      IF (jklg > jk)    CYCLE   ! don't emit below
                      ! emission at this level; if empty then closest above ...
                      ! ... jklg == jk -> nc = GNCB(jgx,jk,jgy)
                      ! ... jklg <  jk -> nc >= GNCB(jgx,jk,jgy)
                      nc = NINT(SUM(GNCB(jgx,jklg:jk,jgy)))
                      ! ... there are still boxes below ...
                      IF (nc > NINT(GNCB(jgx,jklg,jgy))) CYCLE
                      !
                      ! (mol/s) / (kg/(kg/mol)) -> mol/mol/s
                      efld_lg(jc, jt) = efld_lg(jc, jt) + XXREF(ii)%flx(j)  &
                           / ( AMCELL * ( 1.E3_DP/M_air ) ) &
                           / GNCB(jgx,jklg,jgy)
                      ! (mol/s) * Na / m^3 -> molec/m^3/s
                      IF (XPOINT(ii)%io%type == 3) THEN
                         XXREF(ii)%emis(j)%ptr(jc,1,1) = &
                              XXREF(ii)%emis(j)%ptr(jc,1,1) + &
                              XXREF(ii)%flx(j) &
                              * N_A / grvol_gl(jgxlg,jklg,jgylg)
                      ENDIF
                   END DO
                END DO trac_loop_cl

             CASE (2, 4)
                height_loop3: DO k=1, SIZE(XPOINT(ii)%height, 1)-1
                   ! ... levels
                   SELECT CASE (XPOINT(ii)%height_unit)
                   CASE(HU_hPa, HU_Pa)
                      hmin = MIN(XPOINT(ii)%height(k), XPOINT(ii)%height(k+1)) &
                         &   * XPOINT(ii)%height_scal
                      hmax = MAX(XPOINT(ii)%height(k), XPOINT(ii)%height(k+1)) &
                         &   * XPOINT(ii)%height_scal
                      CALL nl_index(press_gl(jgx,1:nlev,jgy), hmin, jk1)
                      CALL nl_index(press_gl(jgx,1:nlev,jgy), hmax, jk2)
                   CASE DEFAULT
                      CALL error_bi('Emission height in metre not implemented for ATTILA!',substr)
                   END SELECT
                   kmin = MIN(jk1, jk2)
                   kmax = MAX(jk1, jk2)

                   trac_loop_ts: DO j=1, XXREF(ii)%ntrac
                      !
                      ! RELEASED TRACER IS NON-TREXP
                      IF (.NOT. XXREF(ii)%ltr(j))  CYCLE
                      ! ... IS NOT CURRENT TRACER
                      IF (XXREF(ii)%idt(j) /= idt) CYCLE
                      !
                      DO jc=1, NCELL
                         jgxlg = NINT(iplon(jc))
                         jgylg = NINT(iplat(jc))
                         jklg  = NINT(iplev(jc))
                         IF (jgxlg /= jgx) CYCLE   ! wrong longitude
                         IF (jgylg /= jgy) CYCLE   ! wrong latitude
                         IF (jklg > kmax)    CYCLE   ! don't emit below
                         IF (jklg < kmin)    CYCLE   ! don't emit above
                         !
                         SELECT CASE (XPOINT(ii)%height_unit)
                         CASE (HU_hPa, HU_Pa)
                            lifc = MIN(pbound(_RI_XYZ__(jgxlg,jgylg,jklg)), &
                               &       pbound(_RI_XYZ__(jgxlg,jgylg,jklg+1)))
                            uifc = MAX(pbound(_RI_XYZ__(jgxlg,jgylg,jklg)), &
                               &       pbound(_RI_XYZ__(jgxlg,jgylg,jklg+1)))
                         CASE DEFAULT
                            CALL error_bi('Emission height in metre not implemented for ATTILA!',substr)
                         END SELECT
                         frac = (MIN(uifc, hmax) - MAX(lifc, hmin)) / (hmax - hmin)
                         SELECT CASE (XPOINT(ii)%emission_unit)
                         CASE (EU_kgs)
                            ! (kg/s) / (kg/mol)/(kg/mol)) -> mol/mol/s
                            efld_lg(jc, jt) = efld_lg(jc, jt)                &
                                 + frac * XPOINT(ii)%flux(k)                 &
                                 * (M_air / XXREF(ii)%molarmass(j)) / AMCELL &
                                 / GNCB(jgx,jklg,jgy)
                            ! (kg/s) / (kg/mol) * Na / m^3 -> molec/m^3/s
                            IF (XPOINT(ii)%io%type == 3) THEN
                               XXREF(ii)%emis(j)%ptr(jc,1,1) =          &
                                    XXREF(ii)%emis(j)%ptr(jc,1,1)       &
                                    + frac * XPOINT(ii)%flux(k)           &
                                    * (1000._DP / XXREF(ii)%molarmass(j)) &
                                    * N_A / grvol_gl(jgxlg,jklg,jgylg)
                            ENDIF
                         CASE (EU_mols)
                            ! (mol/s) / (kg/(kg/mol)) -> mol/mol/s
                            efld_lg(jc, jt) = efld_lg(jc, jt)     &
                                 + frac * XPOINT(ii)%flux(k)      &
                                 / ( AMCELL * ( 1.E3_DP/M_air ) ) &
                                 / GNCB(jgx,jklg,jgy)
                            ! (mol/s) * Na / m^3 -> molec/m^3/s
                            IF (XPOINT(ii)%io%type == 3) THEN
                               XXREF(ii)%emis(j)%ptr(jc,1,1) =    &
                                    XXREF(ii)%emis(j)%ptr(jc,1,1) &
                                    + frac * XPOINT(ii)%flux(k)     &
                                    * N_A / grvol_gl(jgxlg,jklg,jgylg)
                            ENDIF
                         CASE (EU_molecs)
                            ! (molec/s) / Na / (kg/(kg/mol)) -> mol/mol/s
                            efld_lg(jc, jt) = efld_lg(jc, jt)      &
                                 + frac * XPOINT(ii)%flux(k) / N_A &
                                 / ( AMCELL * ( 1.E3_DP/M_air ) )  &
                                 / GNCB(jgx,jklg,jgy)
                            ! (molec/s) / m^3 -> molec/m^3/s
                            IF (XPOINT(ii)%io%type == 3) THEN
                               XXREF(ii)%emis(j)%ptr(jc,1,1) =    &
                                    XXREF(ii)%emis(j)%ptr(jc,1,1) &
                                    + frac * XPOINT(ii)%flux(k)     &
                                    / grvol_gl(jgxlg,jklg,jgylg)
                            ENDIF
                         END SELECT
                      END DO
                   END DO trac_loop_ts
                END DO height_loop3

             END SELECT
             !
          END DO point_loop
          !
       END DO tracer_loop

       ! 2nd: TRACER TENDENCY: EMISSION + DECAY / REACTION
       tracer_loop2: DO jt=1, NTRAC
          !
          IF (TRIM(XTR(jt)%io%trset) /= TRIM(LGTRSTR)) CYCLE
          !
          ! SET TRACER ID
          idt = XTRAC(jt)%idt
          !
          ! REACTION PARTNER
          ! ... start value
          SELECT CASE(XTR(jt)%io%order)
          CASE(0)
             ! DECAY
             k0_lg(:) = XTR(jt)%io%ka
             !
          CASE(1,-1,-2,-3)
             ! REACTION PARTNER
             ! NOTE: in the current implementation, the reaction partner can
             !       only be in GP representation ...
             ALLOCATE(k0_gp(nproma, nlev, ngpblks))
             !
             k0_gp(:,:,:) = XTR(jt)%rp(:,:,:)
             IF (ASSOCIATED(XTR(jt)%rpte)) THEN
                ! TRACER + TENDENCY
                k0_gp(:,:,:) = k0_gp(:,:,:) + XTR(jt)%rpte(:,:,:)*ztmst
             END IF
             ! ... unit conversion -> cm^(-3)
             IF (XTR(jt)%lv2c) THEN
                CALL vmr2conc(k0_gp(:,:,:), pmid(:,:,:), temp_3d(:,:,:))
             ENDIF
             !
             SELECT CASE(XTR(jt)%io%order)
             CASE(1)
                ! -> 1/s
                k0_gp(:,:,:) = k0_gp(:,:,:) * XTR(jt)%io%ka &
                     * exp(- XTR(jt)%io%Ta / temp_3d(:,:,:))
             CASE(-1)
                ! -> 1/s
                k0_gp(:,:,:) = k0_gp(:,:,:) * (XTR(jt)%io%ka + &
                     cair_3d(:,:,:) * XTR(jt)%io%Ta)
             CASE(-2)
                ! special case CH4 + OH -> 1/s
                k0_gp(:,:,:) = k0_gp(:,:,:) * ( &
                     1.85E-20_dp * exp(2.82_dp*log(temp_3d(:,:,:)) &
                     - 987._dp/temp_3d(:,:,:)) )
             CASE(-3)
                ! special case SO2 + OH
                ! reaction rate from mecca.eqn
                k0_gp(:,:,:) = k0_gp(:,:,:) * k_3rd(temp_3d(:,:,:), &
                     cair_3d(:,:,:),3.3E-31_dp, &
                     4.3_dp,1.6E-12_dp,0._dp,0.6_dp)
             END SELECT
             !
             CALL gp2lg_e5(k0_gp, k0_lg)
             !
             DEALLOCATE(k0_gp)
             NULLIFY(k0_gp)
             !
          CASE DEFAULT
             !
             ! SHOULD NEVER BE REACHED
             !
          END SELECT
          !
          ! START VALUE
          zxt(:) = qxtm1_a(:,idt) + qxtte_a(:,idt) * ztmst
          !
          ! INTEGRATE
          CALL solve(zxtte(:), zxt(:), efld_lg(:,jt), k0_lg(:), ztmst)
          !
#ifdef _TFCORR
          ! The "2nd emission step" correction needs to be applied WITHOUT
          ! the current emission tendency, i.e., before qxtte = qxtte + ...
          ! The "2nd emission step" correction can only be applied once
          ! per tracer ... (e.g. for emissions at multiple points into one tracer).
          ! Note: The correction is applied at ANY "2nd emssion step".
          IF (l_tf_corr) THEN
             point_loop_2nd: DO i=1, XTRAC(jt)%np
                !
                ii = XTRAC(jt)%ixp(i)
                !
                IF (XPOINT(ii)%ierr /= 0)    CYCLE
                IF (.NOT. XPOINT(ii)%lnow)   CYCLE   ! ... IS NOT ACTIVE NOW
                !
                ! ONLY FOR CLASSIC POINTS
                IF ((XPOINT(ii)%io%type /= 1) .AND. (XPOINT(ii)%io%type /= 3)) CYCLE
                !
                jgx = XPOINT(ii)%jgx                  ! lon index
                jgy = XPOINT(ii)%jgy                  ! lat index
                SELECT CASE (XPOINT(ii)%height_unit)
                CASE(HU_hPa, HU_Pa)
                   height = XPOINT(ii)%io%height * XPOINT(ii)%height_scal
                CASE DEFAULT
                      CALL error_bi('Emission height in metre not implemented for ATTILA!',substr)
                END SELECT
                CALL nl_index(press_gl(jgx,1:nlev,jgy), height, jk)
                !
                IF (XPOINT(ii)%l2nd) THEN
                   DO jc=1, NCELL
                      jgxlg = NINT(iplon(jc))
                      jgylg = NINT(iplat(jc))
                      jklg  = NINT(iplev(jc))
                      IF (jgxlg /= jgx) CYCLE   ! wrong longitude
                      IF (jgylg /= jgy) CYCLE   ! wrong latitude
                      IF (jklg > jk)    CYCLE   ! don't emit below
                      ! emission at this level; if empty then closest above ...
                      ! ... jklg == jk -> nc = GNCB(jgx,jk,jgy)
                      ! ... jklg <  jk -> nc >= GNCB(jgx,jk,jgy)
                      nc = NINT(SUM(GNCB(jgx,jklg:jk,jgy)))
                      ! ... there are still boxes below ...
                      IF (nc > NINT(GNCB(jgx,jklg,jgy))) CYCLE
                      !
                      qxtte_a(jc, idt) = qxtte_a(jc, idt) &
                           - eps * qxt_a(jc,idt) / ztmst
                   END DO
                END IF
             END DO point_loop_2nd
          END IF
#endif
          !
          ! ADD TENDENCY
          qxtte_a(:,idt) = qxtte_a(:,idt) + zxtte(:)
          !
#ifdef _TFCORR
          ! The "1st emission step" correction needs to be applied WITH
          ! the current emission tendency, i.e., after qxtte = qxtte + ...
          ! The "1st emission step" correction is only be applied once
          ! per tracer ... (e.g. for emissions at multiple points into one tracer).
          ! Note: The correction is applied at ANY "1st emssion step".
          IF (l_tf_corr) THEN
             point_loop_1st: DO i=1, XTRAC(jt)%np
                !
                ii = XTRAC(jt)%ixp(i)
                !
                IF (XPOINT(ii)%ierr /= 0)    CYCLE
                IF (.NOT. XPOINT(ii)%lnow)   CYCLE   ! ... IS NOT ACTIVE NOW
                !
                ! ONLY FOR CLASSIC POINTS
                IF ((XPOINT(ii)%io%type /= 1) .AND. (XPOINT(ii)%io%type /= 3)) CYCLE
                !
                jgx = XPOINT(ii)%jgx                  ! lon index
                jgy = XPOINT(ii)%jgy                  ! lat index
                SELECT CASE (XPOINT(ii)%height_unit)
                CASE(HU_hPa, HU_Pa)
                   height = XPOINT(ii)%io%height * XPOINT(ii)%height_scal
                CASE DEFAULT
                      CALL error_bi('Emission height in metre not implemented for ATTILA!',substr)
                END SELECT
                CALL nl_index(press_gl(jgx,1:nlev,jgy), height, jk)
                !
                IF (XPOINT(ii)%lfirst) THEN
                   DO jc=1, NCELL
                      jgxlg = NINT(iplon(jc))
                      jgylg = NINT(iplat(jc))
                      jklg  = NINT(iplev(jc))
                      IF (jgxlg /= jgx) CYCLE   ! wrong longitude
                      IF (jgylg /= jgy) CYCLE   ! wrong latitude
                      IF (jklg > jk)    CYCLE   ! don't emit below
                      ! emission at this level; if empty then closest above ...
                      ! ... jklg == jk -> nc = GNCB(jgx,jk,jgy)
                      ! ... jklg <  jk -> nc >= GNCB(jgx,jk,jgy)
                      nc = NINT(SUM(GNCB(jgx,jklg:jk,jgy)))
                      ! ... there are still boxes below ...
                      IF (nc > NINT(GNCB(jgx,jklg,jgy))) CYCLE
                      !
                      qxtf_a(jc,idt)=qxtm1_a(jc,idt)+ qxtte_a(jc,idt) * ztmst
                   END DO
                END IF
                !
             END DO point_loop_1st
          END IF
#endif
          !
       END DO tracer_loop2

       ! 3rd: IF REQUIRED: ADD EMISSION TENDENCY FOR ALL NON-TREXP-TRACERS
       force: IF (l_force_emis) THEN
          !
          point_loop2: DO i=1, NPOINT
             !
             IF (TRIM(XPOINT(i)%io%trset) /= TRIM(LGTRSTR)) CYCLE
             !
             IF (XPOINT(i)%ierr /= 0)     CYCLE
             !
             ! GET GLOBAL INDICES
             jgx  = XPOINT(i)%jgx
             jgy  = XPOINT(i)%jgy

             SELECT CASE(XPOINT(i)%io%type)

             CASE (1, 3)

                IF (.NOT. XPOINT(i)%lnow)    CYCLE   ! EMISSION NOT AT THIS TIME

                ! GET LEVEL INDEX
                SELECT CASE (XPOINT(i)%height_unit)
                CASE(HU_hPa, HU_Pa)
                   height = XPOINT(i)%io%height * XPOINT(i)%height_scal
                CASE DEFAULT
                      CALL error_bi('Emission height in metre not implemented for ATTILA!',substr)
                END SELECT
                CALL nl_index(press_gl(jgx,1:nlev,jgy), height, jk)
                !
                trac_loop_cl2: DO j=1, XXREF(i)%ntrac
                   !
                   IF (XXREF(i)%ltr(j)) CYCLE  ! TRACER IS TREXP-TRACER
                   !
                   ! GET TRACER ID
                   idt   = XXREF(i)%idt(j)
                   IF (idt == 0) CYCLE        ! tracer does not exist
                   !
                   ! ADD EMISSION TO TRACER TENDENCY
                   ! (mol/s) / (kg/(kg/mol)) -> mol/mol/s
                   DO jc=1, NCELL
                      jgxlg = NINT(iplon(jc))
                      jgylg = NINT(iplat(jc))
                      jklg  = NINT(iplev(jc))
                      IF (jgxlg /= jgx) CYCLE   ! wrong longitude
                      IF (jgylg /= jgy) CYCLE   ! wrong latitude
                      IF (jklg > jk)    CYCLE   ! don't emit below
                      ! emission at this level; if empty then closest above ...
                      ! ... iplev(jc) == jk -> nc = GNCB(jgx,jk,jgy)
                      ! ... iplev(jc) <  jk -> nc >= GNCB(jgx,jk,jgy)
                      nc = NINT(SUM(GNCB(jgx,jklg:jk,jgy)))
                      ! ... there are still boxes below ...
                      IF (nc > NINT(GNCB(jgx,jklg,jgy))) CYCLE
                      !
#ifdef _TFCORR
                      IF (l_tf_corr) THEN
                         IF (XPOINT(i)%l2nd) THEN
                            qxtte_a(jc, idt) = qxtte_a(jc, idt) &
                                 - eps * qxt_a(jc,idt) / ztmst
                         END IF
                      END IF
#endif
                      !
                      ! (mol/s) / (kg/(kg/mol)) -> mol/mol/s
                      qxtte_a(jc, idt) = qxtte_a(jc, idt)   &
                           + XXREF(i)%flx(j)             &
                           / ( amcell * ( 1.E3_DP/M_air ) ) &
                           / GNCB(jgx,jklg,jgy)
                      ! (mol/s) * Na / m^3 -> molec/m^3/s
                      IF (XPOINT(i)%io%type == 3) THEN
                         XXREF(i)%emis(j)%ptr(jc,1,1) = &
                              XXREF(i)%emis(j)%ptr(jc,1,1) + &
                              XXREF(i)%flx(j) &
                              * N_A / grvol_gl(jgxlg,jklg,jgylg)
                      ENDIF
                      !
#ifdef _TFCORR
                      IF (l_tf_corr) THEN
                         IF (XPOINT(i)%lfirst) THEN
                            qxtf_a(jc,idt)=qxtm1_a(jc,idt)+ qxtte_a(jc,idt) * ztmst
                         END IF
                      END IF
#endif
                      !
                   END DO
                END DO trac_loop_cl2

             CASE (2, 4)
                height_loop4: DO k=1, SIZE(XPOINT(i)%height, 1)-1
                   ! GET LEVEL INDICES
                   SELECT CASE (XPOINT(i)%height_unit)
                   CASE(HU_hPa, HU_Pa)
                      hmin = MIN(XPOINT(i)%height(k), XPOINT(i)%height(k+1)) &
                         &   * XPOINT(i)%height_scal
                      hmax = MAX(XPOINT(i)%height(k), XPOINT(i)%height(k+1)) &
                         &   * XPOINT(i)%height_scal
                      CALL nl_index(press_gl(jgx,1:nlev,jgy), hmin, jk1)
                      CALL nl_index(press_gl(jgx,1:nlev,jgy), hmax, jk2)
                   CASE DEFAULT
                      CALL error_bi('Emission height in metre not implemented for ATTILA!',substr)
                   END SELECT
                   kmin = MIN(jk1, jk2)
                   kmax = MAX(jk1, jk2)

                   trac_loop_ts2: DO j=1, XXREF(i)%ntrac
                      !
                      IF (XXREF(i)%ltr(j)) CYCLE  ! TRACER IS TREXP-TRACER
                      !
                      ! GET TRACER ID
                      idt   = XXREF(i)%idt(j)
                      IF (idt == 0) CYCLE        ! tracer does not exist
                      !
                      ! ADD EMISSION TO TRACER TENDENCY
                      DO jc=1, NCELL
                         jgxlg = NINT(iplon(jc))
                         jgylg = NINT(iplat(jc))
                         jklg  = NINT(iplev(jc))
                         IF (jgxlg /= jgx) CYCLE   ! wrong longitude
                         IF (jgylg /= jgy) CYCLE   ! wrong latitude
                         IF (jklg > kmax)  CYCLE   ! don't emit below
                         IF (jklg < kmin)  CYCLE   ! don't emit above
                         SELECT CASE (XPOINT(i)%height_unit)
                         CASE (HU_hPa, HU_Pa)
                            lifc = MIN(pbound(_RI_XYZ__(jgxlg,jgylg,jklg)), &
                               &       pbound(_RI_XYZ__(jgxlg,jgylg,jklg+1)))
                            uifc = MAX(pbound(_RI_XYZ__(jgxlg,jgylg,jklg)), &
                               &       pbound(_RI_XYZ__(jgxlg,jgylg,jklg+1)))
                         CASE DEFAULT
                            CALL error_bi('Emission height in metre not implemented for ATTILA!',substr)
                         END SELECT
                         frac = (MIN(uifc, hmax) - MAX(lifc, hmin)) / (hmax - hmin)
                         SELECT CASE (XPOINT(i)%emission_unit)
                         CASE (EU_kgs)
                            ! (kg/s) / (kg/mol)/(kg/mol) -> mol/mol/s
                            qxtte_a(jc, idt) = qxtte_a(jc, idt)           &
                               + frac * XPOINT(i)%flux(k)                 &
                               * (M_air / XXREF(i)%molarmass(j)) / amcell &
                               / GNCB(jgx,jklg,jgy)
                            ! (kg/s) / (kg(mol) * Na / m^3 -> molec/m^3/s
                            IF (XPOINT(i)%io%type == 4) THEN
                               XXREF(i)%emis(j)%ptr(jc,1,1) =          &
                                    XXREF(i)%emis(j)%ptr(jc,1,1)       &
                                    + frac * XPOINT(i)%flux(k)           &
                                    * (1000._DP / XXREF(i)%molarmass(j)) &
                                    * N_A / grvol_gl(jgxlg,jklg,jgylg)
                            END IF
                         CASE (EU_mols)
                            ! (mol/s) / (kg/(kg/mol)) -> mol/mol/s
                            qxtte_a(jc, idt) = qxtte_a(jc, idt)   &
                                 + frac * XPOINT(i)%flux(k)       &
                                 / ( amcell * ( 1.E3_DP/M_air ) ) &
                                 / GNCB(jgx,jklg,jgy)
                            ! (mol/s) * Na / m^3 -> molec/m^3/s
                            IF (XPOINT(i)%io%type == 4) THEN
                               XXREF(i)%emis(j)%ptr(jc,1,1) =    &
                                    XXREF(i)%emis(j)%ptr(jc,1,1) &
                                    + frac * XPOINT(i)%flux(k)     &
                                    * N_A / grvol_gl(jgxlg,jklg,jgylg)
                            END IF
                         CASE (EU_molecs)
                            ! (molec/s) / Na / (kg/(kg/mol)) -> mol/mol/s
                            qxtte_a(jc, idt) = qxtte_a(jc, idt)      &
                                 + frac * XPOINT(i)%flux(k) / N_A    &
                                 / ( amcell * ( 1.E3_DP/M_air ) )    &
                                 / GNCB(jgx,jklg,jgy)
                            ! (molec/s) / m^3 -> molec/m^3/s
                            IF (XPOINT(i)%io%type == 4) THEN
                               XXREF(i)%emis(j)%ptr(jc,1,1) =    &
                                    XXREF(i)%emis(j)%ptr(jc,1,1) &
                                    + frac * XPOINT(i)%flux(k)     &
                                    / grvol_gl(jgxlg,jklg,jgylg)
                            END IF
                         END SELECT
                      END DO
                   END DO trac_loop_ts2
                END DO height_loop4

             END SELECT
             !
          END DO point_loop2
          !
       END IF force

       ! CLEAN UP
       IF (ASSOCIATED(press_gl)) THEN
          DEALLOCATE(press_gl)
          NULLIFY(press_gl)
       ENDIF

       IF (ASSOCIATED(grvol_gl)) THEN
          DEALLOCATE(grvol_gl)
          NULLIFY(grvol_gl)
       END IF

       DEALLOCATE(efld_lg)
       NULLIFY(efld_lg)

       DEALLOCATE(k0_lg)
       NULLIFY(k0_lg)

       DEALLOCATE(zxt)
       NULLIFY(zxt)
       DEALLOCATE(zxtte)
       NULLIFY(zxtte)
    END IF atilla
!!#D attila -

!!#D clams +
    clams: IF (L_CL) then
       ! INIT
       ! set rest of last row to 1 for division below ...
       if (npromz .lt. nproma) then
          temp_3d(npromz+1:,:,ngpblks) = 1.0_DP
       endif
       ! globalize pressure (for position information) ...
       CALL trp_gpdc_gpgl(1, pbound, press_gl)
       CALL trp_gpdc_gpgl(1, grvol, grvol_gl)
       CALL trp_gpdc_gpgl(1, spc_g, pc_g)
       !
       ALLOCATE(efld_cl(dnparts_max, NTRAC))
       efld_cl(:,:) = 0.0_DP
       !
       ALLOCATE(k0_cl(dnparts_max))
       k0_cl(:) = 0.0_DP
       !
       ALLOCATE(zxt_c(dnparts_max))
       ALLOCATE(zxtte_c(dnparts_max))

       ! INITIALISE EMISSION DIAGNOSTIC
       DO i=1, NPOINT
          IF (TRIM(XPOINT(i)%io%trset) /= TRIM(CLTRSTR)) CYCLE
          DO jt=1, XXREF(i)%ntrac
             IF (ASSOCIATED(XXREF(i)%emis(jt)%ptr)) THEN
                XXREF(i)%emis(jt)%ptr(:,:,:) = 0.0_dp
             END IF
          END DO
       END DO

       ! 1st: FILL EMISSION FIELD FOR ALL TREXP-TRACERS
       tracer_loop_c: DO jt=1, NTRAC
          !
          IF (TRIM(XTR(jt)%io%trset) /= TRIM(CLTRSTR)) CYCLE
          !
          ! SET TRACER ID
          idt = XTRAC(jt)%idt
          !
          point_loop_c: DO i=1, XTRAC(jt)%np
             !
             ii = XTRAC(jt)%ixp(i)
             IF (XPOINT(ii)%ierr /= 0)    CYCLE
             !
             ! GET INDICES IN GLOBAL FIELD
             ! ... lon, lat
             jgx = XPOINT(ii)%jgx  ! lon index
             jgy = XPOINT(ii)%jgy  ! lat index

             SELECT CASE(XPOINT(ii)%io%type)
             CASE (1, 3)
                IF (.NOT. XPOINT(ii)%lnow)   CYCLE   ! ... IS NOT ACTIVE NOW
                ! ... lev
                SELECT CASE (XPOINT(ii)%height_unit)
                CASE(HU_hPa, HU_Pa)
                   height = XPOINT(ii)%io%height * XPOINT(ii)%height_scal
                CASE DEFAULT
                      CALL error_bi('Emission height in metre not implemented for CLAMS!',substr)
                END SELECT
                CALL nl_index(press_gl(jgx,1:nlev,jgy), height, jk)

                trac_loop_cl_c: DO j=1, XXREF(ii)%ntrac
                   !
                   ! RELEASED TRACER IS NON-TREXP
                   IF (.NOT. XXREF(ii)%ltr(j))  CYCLE
                   ! ... IS NOT CURRENT TRACER
                   IF (XXREF(ii)%idt(j) /= idt) CYCLE
                   !
                   DO jc=1, dnparts
                      jgxlg = NINT(pos(jc,1))
                      jgylg = NINT(pos(jc,3))
                      jklg  = NINT(pos(jc,2))
                      IF (jgxlg /= jgx) CYCLE   ! wrong longitude
                      IF (jgylg /= jgy) CYCLE   ! wrong latitude
                      IF (jklg > jk)    CYCLE   ! don't emit below
                      ! emission at this level; if empty then closest above ...
                      ! ... jklg == jk -> nc = pc_g(jgx,jk,jgy)
                      ! ... jklg <  jk -> nc >= pc_g(jgx,jk,jgy)
                      nc = NINT(SUM(pc_g(jgx,jklg:jk,jgy)))
                      ! ... there are still boxes below ...
                      IF (nc > NINT(pc_g(jgx,jklg,jgy))) CYCLE
                      !
                      ! (mol/s) / (kg/(kg/mol)) -> mol/mol/s
                      efld_cl(jc, jt) = efld_cl(jc, jt) + XXREF(ii)%flx(j)  &
                           / ( CMCELL * ( 1.E3_DP/M_air ) ) &
                           / pc_g(jgx,jklg,jgy)
                      ! (mol/s) * Na / m^3 -> molec/m^3/s
                      IF (XPOINT(ii)%io%type == 3) THEN
                         XXREF(ii)%emis(j)%ptr(jc,1,1) = &
                              XXREF(ii)%emis(j)%ptr(jc,1,1) + &
                              XXREF(ii)%flx(j) &
                              * N_A / grvol_gl(jgxlg,jklg,jgylg)
                      ENDIF
                   END DO
                END DO trac_loop_cl_c

             CASE (2, 4)
                height_loop5: DO k=1, SIZE(XPOINT(ii)%height, 1)-1
                   ! ... levels
                   SELECT CASE (XPOINT(ii)%height_unit)
                   CASE(HU_hPa, HU_Pa)
                      hmin = MIN(XPOINT(ii)%height(k), XPOINT(ii)%height(k+1)) &
                         &   * XPOINT(ii)%height_scal
                      hmax = MAX(XPOINT(ii)%height(k), XPOINT(ii)%height(k+1)) &
                         &   * XPOINT(ii)%height_scal
                      CALL nl_index(press_gl(jgx,1:nlev,jgy), hmin, jk1)
                      CALL nl_index(press_gl(jgx,1:nlev,jgy), hmax, jk2)
                   CASE DEFAULT
                      CALL error_bi('Emission height in metre not implemented for CLAMS!',substr)
                   END SELECT

                   kmin = MIN(jk1, jk2)
                   kmax = MAX(jk1, jk2)

                   trac_loop_ts_c: DO j=1, XXREF(ii)%ntrac
                      !
                      ! RELEASED TRACER IS NON-TREXP
                      IF (.NOT. XXREF(ii)%ltr(j))  CYCLE
                      ! ... IS NOT CURRENT TRACER
                      IF (XXREF(ii)%idt(j) /= idt) CYCLE
                      !
                      DO jc=1, dnparts
                         jgxlg = NINT(pos(jc,1))
                         jgylg = NINT(pos(jc,3))
                         jklg  = NINT(pos(jc,2))
                         IF (jgxlg /= jgx) CYCLE   ! wrong longitude
                         IF (jgylg /= jgy) CYCLE   ! wrong latitude
                         IF (jklg > kmax)    CYCLE   ! don't emit below
                         IF (jklg < kmin)    CYCLE   ! don't emit above
                         !
                         SELECT CASE (XPOINT(ii)%height_unit)
                         CASE (HU_hPa, HU_Pa)
                            lifc = MIN(pbound(_RI_XYZ__(jgxlg,jgylg,jklg)), &
                               &       pbound(_RI_XYZ__(jgxlg,jgylg,jklg+1)))
                            uifc = MAX(pbound(_RI_XYZ__(jgxlg,jgylg,jklg)), &
                               &       pbound(_RI_XYZ__(jgxlg,jgylg,jklg+1)))
                         CASE DEFAULT
                            CALL error_bi('Emission height in metre not implemented for ATTILA!',substr)
                         END SELECT
                         frac = (MIN(uifc, hmax) - MAX(lifc, hmin)) / (hmax - hmin)
                         SELECT CASE (XPOINT(ii)%emission_unit)
                         CASE (EU_kgs)
                            ! (kg/s) / (kg/mol)/(kg/mol)) -> mol/mol/s
                            efld_cl(jc, jt) = efld_cl(jc, jt)                &
                                 + frac * XPOINT(ii)%flux(k)                 &
                                 * (M_air / XXREF(ii)%molarmass(j)) / CMCELL &
                                 / pc_g(jgx,jklg,jgy)
                            ! (kg/s) / (kg/mol) * Na / m^3 -> molec/m^3/s
                            IF (XPOINT(ii)%io%type == 3) THEN
                               XXREF(ii)%emis(j)%ptr(jc,1,1) =          &
                                    XXREF(ii)%emis(j)%ptr(jc,1,1) +     &
                                    frac * XPOINT(ii)%flux(k)             &
                                    * (1000._DP / XXREF(ii)%molarmass(j)) &
                                    * N_A / grvol_gl(jgxlg,jklg,jgylg)
                            ENDIF
                         CASE (EU_mols)
                            ! (mol/s) / (kg/(kg/mol)) -> mol/mol/s
                            efld_cl(jc, jt) = efld_cl(jc, jt)     &
                                 + frac * XPOINT(ii)%flux(k)      &
                                 / ( CMCELL * ( 1.E3_DP/M_air ) ) &
                                 / pc_g(jgx,jklg,jgy)
                            ! (mol/s) * Na / m^3 -> molec/m^3/s
                            IF (XPOINT(ii)%io%type == 3) THEN
                               XXREF(ii)%emis(j)%ptr(jc,1,1) =    &
                                    XXREF(ii)%emis(j)%ptr(jc,1,1) &
                                    + frac * XPOINT(ii)%flux(k)     &
                                    * N_A / grvol_gl(jgxlg,jklg,jgylg)
                            ENDIF
                         CASE (EU_molecs)
                            ! (molec/s) / Na / (kg/(kg/mol)) -> mol/mol/s
                            efld_cl(jc, jt) = efld_cl(jc, jt)      &
                                 + frac * XPOINT(ii)%flux(k) / N_A &
                                 / ( CMCELL * ( 1.E3_DP/M_air ) )  &
                                 / pc_g(jgx,jklg,jgy)
                            ! (molec/s) / m^3 -> molec/m^3/s
                            IF (XPOINT(ii)%io%type == 3) THEN
                               XXREF(ii)%emis(j)%ptr(jc,1,1) =    &
                                    XXREF(ii)%emis(j)%ptr(jc,1,1) &
                                    + frac * XPOINT(ii)%flux(k)     &
                                    / grvol_gl(jgxlg,jklg,jgylg)
                            ENDIF
                         END SELECT
                      END DO
                   END DO trac_loop_ts_c
                END DO height_loop5

             END SELECT
             !
          END DO point_loop_c
          !
       END DO tracer_loop_c

       ! 2nd: TRACER TENDENCY: EMISSION + DECAY / REACTION
       tracer_loop2_c: DO jt=1, NTRAC
          !
          IF (TRIM(XTR(jt)%io%trset) /= TRIM(CLTRSTR)) CYCLE
          ! SET TRACER ID
          idt = XTRAC(jt)%idt
          !
          ! REACTION PARTNER
          ! ... start value
          SELECT CASE(XTR(jt)%io%order)
          CASE(0)
             ! DECAY
             k0_cl(:) = XTR(jt)%io%ka
             !
          CASE(1,-1,-2,-3)
             ! REACTION PARTNER
             ! NOTE: in the current implementation, the reaction partner can
             !       only be in GP representation ...
             ALLOCATE(k0_gp(nproma, nlev, ngpblks))
             !
             k0_gp(:,:,:) = XTR(jt)%rp(:,:,:)
             IF (ASSOCIATED(XTR(jt)%rpte)) THEN
                ! TRACER + TENDENCY
                k0_gp(:,:,:) = k0_gp(:,:,:) + XTR(jt)%rpte(:,:,:)*ztmst
             END IF
             ! ... unit conversion -> cm^(-3)
             IF (XTR(jt)%lv2c) THEN
                CALL vmr2conc(k0_gp(:,:,:), pmid(:,:,:), temp_3d(:,:,:))
             ENDIF
             !
             SELECT CASE(XTR(jt)%io%order)
             CASE(1)
                ! -> 1/s
                k0_gp(:,:,:) = k0_gp(:,:,:) * XTR(jt)%io%ka &
                     * exp(- XTR(jt)%io%Ta / temp_3d(:,:,:))
             CASE(-1)
                ! -> 1/s
                k0_gp(:,:,:) = k0_gp(:,:,:) * (XTR(jt)%io%ka + &
                     cair_3d(:,:,:) * XTR(jt)%io%Ta)
             CASE(-2)
                ! special case CH4 + OH -> 1/s
                k0_gp(:,:,:) = k0_gp(:,:,:) * ( &
                     1.85E-20_dp * exp(2.82_dp*log(temp_3d(:,:,:)) &
                     - 987._dp/temp_3d(:,:,:)) )
             CASE(-3)
                ! special case SO2 + OH
                ! reaction rate from mecca.eqn
                k0_gp(:,:,:) = k0_gp(:,:,:) * k_3rd(temp_3d(:,:,:), &
                     cair_3d(:,:,:),3.3E-31_dp, &
                     4.3_dp,1.6E-12_dp,0._dp,0.6_dp)
             END SELECT
             !
             CALL gp2cl_e5(k0_gp, k0_cl)
             !
             DEALLOCATE(k0_gp)
             NULLIFY(k0_gp)
             !
          CASE DEFAULT
             !
             ! SHOULD NEVER BE REACHED
             !
          END SELECT
          !
          ! START VALUE
          zxt_c(:) = qxt_c(:,idt) + qxtte_c(:,idt) * ztmst
          !
          ! INTEGRATE
          CALL solve(zxtte_c(:), zxt_c(:), efld_cl(:,jt), k0_cl(:), ztmst)
          !
          ! ADD TENDENCY
          qxtte_c(:,idt) = qxtte_c(:,idt) + zxtte_c(:)
          !
          !
       END DO tracer_loop2_c

       ! 3rd: IF REQUIRED: ADD EMISSION TENDENCY FOR ALL NON-TREXP-TRACERS
       force_c: IF (l_force_emis) THEN
          !
          point_loop2_c: DO i=1, NPOINT
             !
             IF (TRIM(XPOINT(i)%io%trset) /= TRIM(CLTRSTR)) CYCLE
             !
             IF (XPOINT(i)%ierr /= 0)     CYCLE
             !
             ! GET GLOBAL INDICES
             jgx  = XPOINT(i)%jgx
             jgy  = XPOINT(i)%jgy

             SELECT CASE(XPOINT(i)%io%type)

             CASE (1, 3)

                IF (.NOT. XPOINT(i)%lnow)    CYCLE   ! EMISSION NOT AT THIS TIME

                ! GET LEVEL INDEX
                SELECT CASE (XPOINT(i)%height_unit)
                CASE (HU_hPa, HU_Pa)
                   height = XPOINT(i)%io%height * XPOINT(i)%height_scal
                CASE DEFAULT
                      CALL error_bi('Emission height in metre not implemented for CLAMS!',substr)
                END SELECT
                CALL nl_index(press_gl(jgx,1:nlev,jgy), height, jk)
                !
                trac_loop_cl2_c: DO j=1, XXREF(i)%ntrac
                   !
                   IF (XXREF(i)%ltr(j)) CYCLE  ! TRACER IS TREXP-TRACER
                   !
                   ! GET TRACER ID
                   idt   = XXREF(i)%idt(j)
                   IF (idt == 0) CYCLE        ! tracer does not exist
                   !
                   ! ADD EMISSION TO TRACER TENDENCY
                   ! (mol/s) / (kg/(kg/mol)) -> mol/mol/s
                   DO jc=1, dnparts
                      jgxlg = NINT(pos(jc,1))
                      jgylg = NINT(pos(jc,3))
                      jklg  = NINT(pos(jc,2))
                      IF (jgxlg /= jgx) CYCLE   ! wrong longitude
                      IF (jgylg /= jgy) CYCLE   ! wrong latitude
                      IF (jklg > jk)    CYCLE   ! don't emit below
                      ! emission at this level; if empty then closest above ...
                      ! ... iplev(jc) == jk -> nc = GNCB(jgx,jk,jgy)
                      ! ... iplev(jc) <  jk -> nc >= GNCB(jgx,jk,jgy)
                      nc = NINT(SUM(pc_g(jgx,jklg:jk,jgy)))
                      ! ... there are still boxes below ...
                      IF (nc > NINT(pc_g(jgx,jklg,jgy))) CYCLE
                      !
                      ! (mol/s) / (kg/(kg/mol)) -> mol/mol/s
                      qxtte_c(jc, idt) = qxtte_c(jc, idt)   &
                           + XXREF(i)%flx(j)                &
                           / ( cmcell * ( 1.E3_DP/M_air ) ) &
                           / pc_g(jgx,jklg,jgy)
                      ! (mol/s) * Na / m^3 -> molec/m^3/s
                      IF (XPOINT(i)%io%type == 3) THEN
                         XXREF(i)%emis(j)%ptr(jc,1,1) =    &
                              XXREF(i)%emis(j)%ptr(jc,1,1) &
                              + XXREF(i)%flx(j)              &
                              * N_A / grvol_gl(jgxlg,jklg,jgylg)
                      ENDIF
                      !
                   END DO
                END DO trac_loop_cl2_c

             CASE (2, 4)
                height_loop6: DO k=1, SIZE(XPOINT(i)%height, 1)-1
                   ! GET LEVEL INDICES
                   SELECT CASE (XPOINT(i)%height_unit)
                   CASE (HU_hPa, HU_Pa)
                      hmin = MIN(XPOINT(i)%height(k), XPOINT(i)%height(k+1)) &
                         &   * XPOINT(i)%height_scal
                      hmax = MAX(XPOINT(i)%height(k), XPOINT(i)%height(k+1)) &
                         &   * XPOINT(i)%height_scal
                      CALL nl_index(press_gl(jgx,1:nlev,jgy), hmin, jk1)
                      CALL nl_index(press_gl(jgx,1:nlev,jgy), hmax, jk2)
                   CASE DEFAULT
                      CALL error_bi('Emission height in metre not implemented for CLAMS!',substr)
                   END SELECT

                   kmin = MIN(jk1, jk2)
                   kmax = MAX(jk1, jk2)

                   trac_loop_ts2_c: DO j=1, XXREF(i)%ntrac
                      !
                      IF (XXREF(i)%ltr(j)) CYCLE  ! TRACER IS TREXP-TRACER
                      !
                      ! GET TRACER ID
                      idt   = XXREF(i)%idt(j)
                      IF (idt == 0) CYCLE        ! tracer does not exist
                      !
                      ! ADD EMISSION TO TRACER TENDENCY
                      ! (molec/s) / Na / (kg/(kg/mol)) -> mol/mol/s
                      DO jc=1, dnparts
                         jgxlg = NINT(pos(jc,1))
                         jgylg = NINT(pos(jc,3))
                         jklg  = NINT(pos(jc,2))
                         IF (jgxlg /= jgx) CYCLE   ! wrong longitude
                         IF (jgylg /= jgy) CYCLE   ! wrong latitude
                         IF (jklg > kmax)  CYCLE   ! don't emit below
                         IF (jklg < kmin)  CYCLE   ! don't emit above
                         SELECT CASE (XPOINT(i)%height_unit)
                         CASE (HU_hPa, HU_Pa)
                            lifc = MIN(pbound(_RI_XYZ__(jgxlg,jgylg,jklg)), &
                               &       pbound(_RI_XYZ__(jgxlg,jgylg,jklg+1)))
                            uifc = MAX(pbound(_RI_XYZ__(jgxlg,jgylg,jklg)), &
                               &       pbound(_RI_XYZ__(jgxlg,jgylg,jklg+1)))
                         CASE DEFAULT
                            CALL error_bi('Emission height in metre not implemented for ATTILA!',substr)
                         END SELECT
                         frac = (MIN(uifc, hmax) - MAX(lifc, hmin)) / (hmax - hmin)
                         SELECT CASE (XPOINT(i)%emission_unit)
                         CASE (EU_kgs)
                            ! (kg/s) / (kg/mol)/(kg/mol) -> mol/mol/s
                            qxtte_c(jc, idt) = qxtte_c(jc, idt)           &
                               + frac * XPOINT(i)%flux(k)                 &
                               * (M_air / XXREF(i)%molarmass(j)) / cmcell &
                               / pc_g(jgx,jklg,jgy)
                            ! (kg/s) / (kg(mol) * Na / m^3 -> molec/m^3/s
                            IF (XPOINT(i)%io%type == 4) THEN
                               XXREF(i)%emis(j)%ptr(jc,1,1) =          &
                                    XXREF(i)%emis(j)%ptr(jc,1,1)       &
                                    + frac * XPOINT(i)%flux(k)           &
                                    * (1000._DP / XXREF(i)%molarmass(j)) &
                                    * N_A / grvol_gl(jgxlg,jklg,jgylg)
                            END IF
                         CASE (EU_mols)
                            ! (mol/s) / (kg/(kg/mol)) -> mol/mol/s
                            qxtte_c(jc, idt) = qxtte_c(jc, idt)   &
                                 + frac * XPOINT(i)%flux(k)       &
                                 / ( cmcell * ( 1.E3_DP/M_air ) ) &
                                 / pc_g(jgx,jklg,jgy)
                            ! (mol/s) * Na / m^3 -> molec/m^3/s
                            IF (XPOINT(i)%io%type == 4) THEN
                               XXREF(i)%emis(j)%ptr(jc,1,1) =    &
                                    XXREF(i)%emis(j)%ptr(jc,1,1) &
                                    + frac * XPOINT(i)%flux(k)     &
                                    * N_A / grvol_gl(jgxlg,jklg,jgylg)
                            END IF
                         CASE (EU_molecs)
                            ! (molec/s) / Na / (kg/(kg/mol)) -> mol/mol/s
                            qxtte_c(jc, idt) = qxtte_c(jc, idt)      &
                                 + frac * XPOINT(i)%flux(k) / N_A    &
                                 / ( cmcell * ( 2.E3_DP/M_air ) )    &
                                 / pc_g(jgx,jklg,jgy)
                            ! (molec/s) / m^3 -> molec/m^3/s
                            IF (XPOINT(i)%io%type == 4) THEN
                               XXREF(i)%emis(j)%ptr(jc,1,1) =    &
                                    XXREF(i)%emis(j)%ptr(jc,1,1) &
                                    + frac * XPOINT(i)%flux(k)     &
                                    / grvol_gl(jgxlg,jklg,jgylg)
                            END IF
                         END SELECT

                      END DO
                   END DO trac_loop_ts2_c
                END DO height_loop6

             END SELECT
             !
          END DO point_loop2_c
          !
       END IF force_c

       ! CLEAN UP
       IF (ASSOCIATED(press_gl)) THEN
          DEALLOCATE(press_gl)
          NULLIFY(press_gl)
       ENDIF

       IF (ASSOCIATED(grvol_gl)) THEN
          DEALLOCATE(grvol_gl)
          NULLIFY(grvol_gl)
       END IF

       DEALLOCATE(efld_cl)
       NULLIFY(efld_cl)

       DEALLOCATE(k0_cl)
       NULLIFY(k0_cl)

       DEALLOCATE(zxt_c)
       NULLIFY(zxt_c)
       DEALLOCATE(zxtte_c)
       NULLIFY(zxtte_c)
    ENDIF clams
!!#D clams -

#endif

  END SUBROUTINE trexp_global_end
  ! -------------------------------------------------------------------

  ! -------------------------------------------------------------------
  SUBROUTINE trexp_free_memory

    IMPLICIT NONE

    INTRINSIC :: ASSOCIATED

    ! LOCAL
    INTEGER :: i, jt

    DO jt=1, NTRAC
       IF (ASSOCIATED(XTRAC(jt)%ixp)) THEN
          DEALLOCATE(XTRAC(jt)%ixp)
          NULLIFY(XTRAC(jt)%ixp)
       END IF
    END DO

    DO i=1, NPOINT
       IF (ASSOCIATED(XXREF(i)%idt)) THEN
          DEALLOCATE(XXREF(i)%idt) ; NULLIFY(XXREF(i)%idt)
       END IF
       IF (ASSOCIATED(XXREF(i)%flx)) THEN
          DEALLOCATE(XXREF(i)%flx) ; NULLIFY(XXREF(i)%flx)
       END IF
       IF (ASSOCIATED(XXREF(i)%ltr)) THEN
          DEALLOCATE(XXREF(i)%ltr) ; NULLIFY(XXREF(i)%ltr)
       END IF
        IF (ASSOCIATED(XXREF(i)%emis)) THEN
           DEALLOCATE(XXREF(i)%emis) ; NULLIFY(XXREF(i)%emis)
        END IF
     END DO

    IF (ASSOCIATED(temp_3d)) THEN
       DEALLOCATE(temp_3d) ; NULLIFY(temp_3d)
    END IF
    IF (ASSOCIATED(cair_3d)) THEN
       DEALLOCATE(cair_3d) ; NULLIFY(cair_3d)
    END IF

  END SUBROUTINE trexp_free_memory
  ! -------------------------------------------------------------------

  ! -------------------------------------------------------------------
  SUBROUTINE trexp_read_nml_cpl(status, iou)

    ! trexp MODULE ROUTINE (ECHAM-5 INTERFACE, PRIVATE)
    !
    ! read namelist for 'coupling' to ECHAM5
    !
    ! Author: Patrick Joeckel, MPICH, Dec 2004

    USE messy_main_tools, ONLY: read_nml_open, read_nml_check, read_nml_close
    USE messy_main_tracer_mem_bi, ONLY: NGCELL
    USE messy_main_channel_bi,    ONLY: REPR_LG_CLAMS
    USE messy_main_channel_bi,    ONLY: REPR_UNDEF

    IMPLICIT NONE

    ! I/O
    INTEGER, INTENT(OUT) :: status     ! error status
    INTEGER, INTENT(IN)  :: iou        ! I/O unit

    ! (LOCAL) NAMELIST VARIABLES
    CHARACTER(LEN=*), PARAMETER :: substr = 'trexp_read_nml_cpl'

    NAMELIST /CPL/ L_GP, L_LG, L_CL &
         , l_force_emis, l_allow_ext, l_tf_corr, TR, POINT

    ! LOCAL
    LOGICAL              :: lex      ! file exists ?
    INTEGER              :: fstat    ! file status

    status = 1

    ! INITIALIZE NAMELIST VARIABLES
    ! NOTE: already at definition

    CALL read_nml_open(lex, substr, iou, 'CPL', modstr)
    IF (.not.lex) RETURN    ! <modstr>.nml does not exist

    READ(iou, NML=CPL, IOSTAT=fstat)
    CALL read_nml_check(fstat, substr, iou, 'CPL', modstr)
    IF (fstat /= 0) RETURN  ! error while reading namelist

    CALL read_nml_close(substr, iou, modstr)

    IF ((.NOT.L_GP).AND.(.NOT.L_LG) .AND. (.NOT. L_CL)) THEN
       WRITE(*,*) '*** ERROR: GRIDPOINT AND 2 LAGRANGIAN INTEGRATIONS'
       WRITE(*,*) '           ALL 3 SWITCHED OFF!'
       RETURN ! ERROR
    END IF

    IF (L_LG .AND. NGCELL <= 0) THEN
!!#D attila +
       WRITE(*,*) 'WARNING: LAGRANGIAN INTEGRATION REQUIRES ATTILA'
       WRITE(*,*) ' -> SWITCHING OFF LAGRANGIAN (ATTILA) INTEGRATION FOR TREXP'
!!#D attila -
       L_LG = .FALSE.
    END IF

!!#D clams +
     IF (L_CL .AND. REPR_LG_CLAMS == REPR_UNDEF ) THEN
       WRITE(*,*) 'WARNING: LAGRANGIAN INTEGRATION REQUIRES CLaMS!'
       WRITE(*,*) ' -> SWITCHING OFF LAGRANGIAN (CLaMS) INTEGRATION FOR TREXP'
       L_CL = .FALSE.
    END IF
!!#D clams -

    status = 0  ! no ERROR

  END SUBROUTINE trexp_read_nml_cpl
  ! -------------------------------------------------------------------

! **********************************************************************
END MODULE messy_trexp_si
! **********************************************************************

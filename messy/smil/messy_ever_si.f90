#include "messy_main_ppd_bi.inc"

!
#ifdef ECHAM5
#define _TFCORR
#endif


MODULE messy_ever_si
! **********************************************************************

  ! Explosive Volcanic ERuptions (EVER)
  !
  ! INTERFACE FOR ECHAM5 (MESSy/SMIL)
  !
  ! Author: Matthias Kohl, MPIC, February 2021
  !
  ! - Basic structure taken from TREXP
  ! - allows emissions of all tracers (gas and aerosols) following the
  !   explosive eruptions of volcanoes
  ! - specialized on vertical emission distributions:
  !     - currently possible: linear and Gaussian shape
  !     - emission constantly over a provided time range
  ! - aerosol emissions currently only working combined with GMXe
  !

  USE messy_main_blather_bi,    ONLY: start_message_bi, end_message_bi
#ifdef MESSYTENDENCY
  !tendency budget
  USE messy_main_tendency_bi,   ONLY: mtend_get_handle,       &
                                      mtend_get_start_l,      &
                                      mtend_add_l,            &
                                      mtend_register,           &    
                                      mtend_id_tracer,        &
                                      mtend_id_t, mtend_id_q
#endif
  ! MESSy
  USE messy_main_constants_mem, ONLY: STRLEN_MEDIUM
  USE messy_main_tools,         ONLY: PTR_3D_ARRAY
  USE messy_ever

  IMPLICIT NONE
  PRIVATE

  INTRINSIC :: NULL

  ! NAMELIST LIMITS
  INTEGER, PARAMETER :: NMAXERUPT = 300 ! NUMBER OF VOLCANIC ERUPTIONS
  INTEGER, PARAMETER :: NMAXAERPR = 50  ! NUMBER OF aerosol parameter sets
  INTEGER, PARAMETER :: MAXSTRLEN = 200 ! max. string length

  ! USER INTERFACE
  
  ! ... define aerosol parameters ...
  TYPE AER_IO
     !
     CHARACTER(LEN=MAXSTRLEN) :: par_name = ''     ! name of aerosol parameter set
     REAL(dp)                 :: density  = 0.0_DP ! density of aerosols
     REAL(dp)                 :: diameter = 0.0_DP ! median diameter of aerosols [m]
     REAL(dp)                 :: sigma    = 0.0_DP ! sigma of log-normal aerosols distribution
     CHARACTER(LEN=2)         :: mode     = ''     ! aerosol mode as defined in aer submodel 
  END TYPE AER_IO
  

  ! ... release points in space and time
  TYPE ERUPT_POINT_IO
     !
     REAL(DP) :: lon    =    0.0_DP  ! degrees east
     REAL(DP) :: lat    =    0.0_DP  ! degreas north
     !
     INTEGER  :: vshape     = 1
     REAL(DP) :: minalt     = 0.0_DP  ! altitude lower limit [km]
     REAL(DP) :: maxalt     = 0.0_DP  ! altitude upper limit [km]
     !
     ! only for shape 1 (Gaussian)
     REAL(DP) :: meanalt    = 0.0_DP  ! mean value of vertical Gauss [km]
     REAL(DP) :: vsigma     = 0.0_DP  ! sigma of vertical Gauss [km]
     ! timespan
     INTEGER, DIMENSION(6) :: ds = (/ 0,0,0,0,0,0 /)  ! START: YY, MM, DD, ...
     INTEGER, DIMENSION(6) :: de = (/ 0,0,0,0,0,0 /) ! END:   ... HR, MI, SE
     !
     ! tracers
     CHARACTER(LEN=MAXSTRLEN) :: trlist = '' ! ;-separated list of tracers
     CHARACTER(LEN=MAXSTRLEN) :: trmass = '' ! ;-separated mass list for tracers in [Tg]
     CHARACTER(LEN=MAXSTRLEN) :: traepr = '' ! ;-separated list of aerosol parameters (empty for non-aerosol species)
     !
  END TYPE ERUPT_POINT_IO

  ! WORKSPACE
  ! AEROSOL PARAMETERS
  TYPE AEROSOL_PARS
     TYPE(AER_IO)           :: io
  END TYPE AEROSOL_PARS

  ! ... release points in decomposition; timing
  TYPE ERUPTION_POINT
     TYPE(ERUPT_POINT_IO) :: io
     !
     INTEGER                :: jgx    = 0   ! longitude index (global)
     INTEGER                :: jgy    = 0   ! latitude index (global)
     !
     INTEGER                :: pe     = 0   ! on which CPU [0...p_ncpus]
     INTEGER                :: jp     = 0   ! which column [1...kproma]
     INTEGER                :: jrow   = 0   ! which row    [1...jrow]
     INTEGER                :: ierr   = 0   ! status
     !
     REAL(DP)               :: dt     = 0   ! emission time in seconds
     LOGICAL                :: lnow   = .FALSE. ! EMISSION NOW ?
     !
     REAL(DP), POINTER      :: rnow   => NULL() ! -> channel object
     !
     REAL(DP)               :: estart = 0._DP   ! time of emission start (jd)
     REAL(DP)               :: estop  = 0._DP   ! time when emission ends (jd)
#ifdef _TFCORR
     LOGICAL                :: lfirst = .FALSE.  ! first time step of emission
     LOGICAL                :: l2nd   = .FALSE.  ! 2nd time step of emission
#endif
     !
  END TYPE ERUPTION_POINT

  ! ... cross reference to tracers
  TYPE ERUPT_POINT_XREF
     INTEGER                         :: ntrac = 0       ! NUMBER OF ASSOC. TRAC.
     INTEGER,  DIMENSION(:), POINTER :: idt  => NULL()  ! TRACER IDs
     REAL(DP), DIMENSION(:), POINTER :: flx  => NULL()  ! [mol/s] emission flux
     INTEGER,  DIMENSION(:), POINTER :: nidt => NULL()  ! TRACER IDs of number fluxes
     REAL(DP), DIMENSION(:), POINTER :: nflx => NULL()  ! [1/s] number flux
  END TYPE ERUPT_POINT_XREF

  ! WORKSPACE
  ! ... representation independent
  TYPE(AER_IO),                DIMENSION(NMAXAERPR), SAVE :: AER
  TYPE(AEROSOL_PARS),          DIMENSION(NMAXAERPR), SAVE :: XAER
  TYPE(ERUPT_POINT_IO),        DIMENSION(NMAXERUPT), SAVE :: POINT
  TYPE(ERUPTION_POINT),        DIMENSION(NMAXERUPT), SAVE :: XPOINT
  !
  LOGICAL :: l_tf_corr = .FALSE. ! time filter correction
  LOGICAL :: l_gmxe = .FALSE.    ! flag if GMXe is used
  ! ... representation dependent
  TYPE(ERUPT_POINT_XREF),   DIMENSION(NMAXERUPT), SAVE :: XXREF_GP

  ! ACTUAL NUMBERS
  INTEGER, SAVE :: NPAR   = 0  ! number of aerosol parameter sets
  INTEGER, SAVE :: NPOINT = 0  ! number of release points

  ! REPRESENTATION SWITCHES (currently only implemented for gp)

#ifdef MESSYTENDENCY
  INTEGER :: my_handle
#endif

  ! FOR ADDITIONAL (DIAGNOSTIC) CHANNEL OBJECTS (EMISSIONS)
  TYPE(PTR_3D_ARRAY), DIMENSION(:), POINTER :: emis_gp => NULL()

  PUBLIC :: ever_initialize
  PUBLIC :: ever_init_memory
  PUBLIC :: ever_init_coupling
  !         ! -> init_trac_cpl
  PUBLIC :: ever_global_start
  PUBLIC :: ever_physc
  PUBLIC :: ever_free_memory

  !PRIVATE: ever_read_nml_cpl

CONTAINS

  ! -------------------------------------------------------------------
  SUBROUTINE ever_initialize

    ! ever MODULE ROUTINE (ECHAM-5 INTERFACE)
    !

    ! BML/MESSy
    USE messy_main_mpi_bi,        ONLY: p_parallel_io, p_io, p_bcast
    USE messy_main_blather_bi,    ONLY: error_bi
    USE messy_main_tools,         ONLY: find_next_free_unit
    USE messy_main_timer,         ONLY: time_span_d, gregor2julian

    IMPLICIT NONE

    INTRINSIC :: REAL, TRIM, TINY, ABS

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'ever_initialize'
    INTEGER                     :: iou    ! I/O unit
    INTEGER                     :: status ! error status
    INTEGER                     :: jt, i
    CHARACTER(LEN=2*STRLEN_MEDIUM+1) :: fullname
    
    ! INITIALIZE CPL
    IF (p_parallel_io) THEN
       iou = find_next_free_unit(100,200)
       CALL ever_read_nml_cpl(status, iou)
       IF (status /= 0) CALL error_bi(' ',substr)
    END IF
    CALL p_bcast(l_tf_corr, p_io)
    CALL p_bcast(l_gmxe, p_io)

#ifdef MESSYTENDENCY
    my_handle = mtend_get_handle(modstr)
#endif

    CALL start_message_bi(modstr,'INITIALISATION ',substr)

    IF (p_parallel_io) THEN
       WRITE(*,*) '======================================================'
       WRITE(*,*) '             AEROSOL PARAMETER SETS                   '
       WRITE(*,*) '======================================================'
    END IF




    IF (p_parallel_io) THEN
       ! ############### AEROSOL PARAMETERS #################################
       ! GET NUMBER OF AEROSOL PARAMETER SETS, COPY AND BROADCAST THEM
       NPAR = 1
       ! COPY DATA AND PARSE STR
       parameter_loop: DO jt=1, NMAXAERPR
          IF (TRIM(AER(jt)%par_name) == '') CYCLE
          XAER(NPAR)%io%par_name = AER(jt)%par_name
          WRITE(*,*) 'Aerosol Parameter Set ''',TRIM(XAER(NPAR)%io%par_name),''' ...'
          !
          IF (AER(jt)%density <= 0.0) THEN
             WRITE(*,*) '... WARNING: density <= 0 ...'
          ELSE
             XAER(NPAR)%io%density = AER(jt)%density
             WRITE(*,*) '    ... density : ' &
                  , XAER(NPAR)%io%density
             !
          END IF
          !
          IF (AER(jt)%diameter <= 0.0) THEN
             WRITE(*,*) '... WARNING: diameter <= 0 ...'
          ELSE
             XAER(NPAR)%io%diameter = AER(jt)%diameter
             WRITE(*,*) '    ... diameter : ' &
                  , XAER(NPAR)%io%diameter
             !
          
          ENDIF

          IF (AER(jt)%sigma <= 0.0) THEN
             WRITE(*,*) '... WARNING: sigma <= 0 ...'
          ELSE
             XAER(NPAR)%io%sigma = AER(jt)%sigma
             WRITE(*,*) '    ... sigma : ' &
                  , XAER(NPAR)%io%sigma
             !
          ENDIF

          IF (TRIM(AER(jt)%mode) == '') THEN
             WRITE(*,*) '... WARNING: no aerosol mode provided ...'
          ELSE
             XAER(NPAR)%io%mode = AER(jt)%mode
             WRITE(*,*) '    ... mode : ''' &
                  , TRIM(XAER(NPAR)%io%mode), ''' .'
             !
          
          ENDIF


          ! NEXT PARAMETER
          NPAR = NPAR + 1
          WRITE(*,*) '------------------------------------------------------'
       END DO parameter_loop
       NPAR = NPAR - 1
    END IF
    !
    CALL p_bcast(NPAR, p_io)
    !
    ! BROADCAST RESULTS
    DO jt=1, NPAR
       CALL p_bcast(XAER(jt)%io%par_name,    p_io)
       CALL p_bcast(XAER(jt)%io%density,     p_io)
       CALL p_bcast(XAER(jt)%io%diameter,    p_io)
       CALL p_bcast(XAER(jt)%io%sigma,       p_io)
       CALL p_bcast(XAER(jt)%io%mode,        p_io)
    END DO

    IF (p_parallel_io) THEN
       WRITE(*,*) ' ---> ',NPAR,' AEROSOL PARAMETER(S) SPECIFIED !'
    END IF

    IF (p_parallel_io) THEN
       WRITE(*,*) '======================================================'
       WRITE(*,*) '               VOLCANIC ERUPTIONS                     '
       WRITE(*,*) '======================================================'
    END IF

    IF (p_parallel_io) THEN
       ! ############### ERUPTION POINTS #################################
       ! GET NUMBER OF ERUPTION POINTS, COPY AND BROADCAST THEM
       NPOINT = 1
       ! COPY DATA AND PARSE STR
       point_loop: DO i=1, NMAXERUPT
          !
          ! TRACER LIST
          IF (TRIM(POINT(i)%trlist) == '') THEN
             WRITE(*,*) '    ... empty tracer list ... skipping ...'
             CYCLE
          END IF
          !
          XPOINT(NPOINT)%io%trlist  = TRIM(POINT(i)%trlist)
          WRITE(*,*) '    ... TRACERS  : ',TRIM(XPOINT(NPOINT)%io%trlist)
          !
          ! MASS LIST
          IF (TRIM(POINT(i)%trmass) == '') THEN
             WRITE(*,*) '    ... empty tracer mass list ... skipping ...'
             CYCLE
          END IF
          !
          XPOINT(NPOINT)%io%trmass  = TRIM(POINT(i)%trmass)
          WRITE(*,*) '    ... TRACER MASSES  : ',TRIM(XPOINT(NPOINT)%io%trmass)
          !
          ! PARAMETER LIST
          IF (TRIM(POINT(i)%traepr) == '') THEN
             WRITE(*,*) '    ... empty tracer aerosol parameter list ... skipping ...'
             CYCLE
          END IF
          !
          XPOINT(NPOINT)%io%traepr  = TRIM(POINT(i)%traepr)
          WRITE(*,*) '    ... TRACER AEROSOL PARAMETERS  : ',TRIM(XPOINT(NPOINT)%io%traepr)
          !
          !
          IF ((POINT(i)%lon < -180.0_DP) .OR. (POINT(i)%lon > 180.0_DP)) THEN
             WRITE(*,*) '    ... LONGITUDE out of range [-180, 180] : ' &
                  ,POINT(i)%lon,''' ... skipping'
             CYCLE
          ELSE
             XPOINT(NPOINT)%io%lon = POINT(i)%lon
             WRITE(*,*) '    ... LONGITUDE: ', XPOINT(NPOINT)%io%lon
          END IF
          !
          IF ((POINT(i)%lat < -90.0_DP) .OR. (POINT(i)%lat > 90.0_DP)) THEN
             WRITE(*,*) '    ... LATITUDE out of range [-90, 90] : ' &
                  ,POINT(i)%lat,''' ... skipping'
             CYCLE
          ELSE
             XPOINT(NPOINT)%io%lat = POINT(i)%lat
             WRITE(*,*) '    ... LATITUDE : ', XPOINT(NPOINT)%io%lat
          END IF
          !
          !
          XPOINT(NPOINT)%io%vshape = POINT(i)%vshape
          WRITE(*,*) '    ... SHAPE : ' &
               , XPOINT(NPOINT)%io%vshape
          !
          XPOINT(NPOINT)%io%minalt = POINT(i)%minalt
          WRITE(*,*) '    ... MIN ALTITUDE : ' &
               , XPOINT(NPOINT)%io%minalt,' [km]'
          !
          XPOINT(NPOINT)%io%maxalt = POINT(i)%maxalt
          WRITE(*,*) '    ... MAX ALTITUDE : ' &
               , XPOINT(NPOINT)%io%maxalt,' [km]'
          !
          ! shape 2 (Gaussian)
          IF (XPOINT(NPOINT)%io%vshape == 2) THEN
             !
             XPOINT(NPOINT)%io%meanalt = POINT(i)%meanalt
             WRITE(*,*) '    ... MEAN ALTITUDE : ' &
                  , XPOINT(NPOINT)%io%meanalt,' [km]'
             !
             XPOINT(NPOINT)%io%vsigma = POINT(i)%vsigma
             WRITE(*,*) '    ... sigma of Gaussian vertical distribution : ' &
                  , XPOINT(NPOINT)%io%vsigma,' [km]'
             !
          ENDIF
          !
          XPOINT(NPOINT)%io%ds(:) = POINT(i)%ds(:)
          WRITE(*,*) '    ... START    : ', &
               XPOINT(NPOINT)%io%ds(1)   ,'-' &
               , XPOINT(NPOINT)%io%ds(2) ,'-' &
               , XPOINT(NPOINT)%io%ds(3) ,' ' &
               , XPOINT(NPOINT)%io%ds(4) ,':' &
               , XPOINT(NPOINT)%io%ds(5) ,':' &
               , XPOINT(NPOINT)%io%ds(6)
          XPOINT(NPOINT)%estart = &
               gregor2julian(XPOINT(NPOINT)%io%ds(1) &
               , XPOINT(NPOINT)%io%ds(2) &
               , XPOINT(NPOINT)%io%ds(3), XPOINT(NPOINT)%io%ds(4) &
               , XPOINT(NPOINT)%io%ds(5), XPOINT(NPOINT)%io%ds(6))
          !
          XPOINT(NPOINT)%io%de(:) = POINT(i)%de(:)
          WRITE(*,*) '    ... END      : ', &
               XPOINT(NPOINT)%io%de(1)   ,'-' &
               , XPOINT(NPOINT)%io%de(2) ,'-' &
               , XPOINT(NPOINT)%io%de(3) ,' ' &
               , XPOINT(NPOINT)%io%de(4) ,':' &
               , XPOINT(NPOINT)%io%de(5) ,':' &
               , XPOINT(NPOINT)%io%de(6)
          XPOINT(NPOINT)%estop = &
               gregor2julian(XPOINT(NPOINT)%io%de(1) &
               , XPOINT(NPOINT)%io%de(2) &
               , XPOINT(NPOINT)%io%de(3), XPOINT(NPOINT)%io%de(4) &
               , XPOINT(NPOINT)%io%de(5), XPOINT(NPOINT)%io%de(6))
          
          ! CALCULATE TIME INTERVAL IN SECONDS
          CALL time_span_d(XPOINT(NPOINT)%dt, &
               XPOINT(NPOINT)%io%ds(1), XPOINT(NPOINT)%io%ds(2) &
               , XPOINT(NPOINT)%io%ds(3) &
               , XPOINT(NPOINT)%io%ds(4), XPOINT(NPOINT)%io%ds(5) &
               , XPOINT(NPOINT)%io%ds(6) &
               , XPOINT(NPOINT)%io%de(1), XPOINT(NPOINT)%io%de(2) &
               , XPOINT(NPOINT)%io%de(3) &
               , XPOINT(NPOINT)%io%de(4), XPOINT(NPOINT)%io%de(5) &
               , XPOINT(NPOINT)%io%de(6) &
               )
          ! days -> seconds
          XPOINT(NPOINT)%dt = XPOINT(NPOINT)%dt * 86400.0_dp
          IF (XPOINT(NPOINT)%dt < 0.0_DP) THEN
             WRITE(*,*) '    ... TIME INTERVAL < 0 : ' &
                  ,XPOINT(NPOINT)%dt,' ... skipping'
             CYCLE 
          END IF
          IF (ABS(XPOINT(NPOINT)%dt) < TINY(0.0_DP)) THEN
             WRITE(*,*) '    ... TIME INTERVAL ~ 0 : ' &
                  ,XPOINT(NPOINT)%dt,' ... skipping'
             CYCLE 
          END IF
          WRITE(*,*) '    ... dT       : ',XPOINT(NPOINT)%dt,'[s]'
          !
          ! NEXT POINT
          NPOINT = NPOINT + 1
          WRITE(*,*) '------------------------------------------------------'
       END DO point_loop
       NPOINT = NPOINT - 1
    END IF
    CALL p_bcast(NPOINT, p_io)

    ! BROADCAST RESULTS
    DO i=1, NPOINT
       !

       CALL p_bcast(XPOINT(i)%io%lon,     p_io)
       CALL p_bcast(XPOINT(i)%io%lat,     p_io)
       CALL p_bcast(XPOINT(i)%io%trlist,  p_io)
       CALL p_bcast(XPOINT(i)%io%trmass,  p_io)
       CALL p_bcast(XPOINT(i)%io%traepr, p_io)
       ! jgx, jgy                          -> set in ever_init_coupling
       ! pe, jp, jrow, ierr                -> set in ever_init_coupling
       !
       CALL p_bcast(XPOINT(i)%io%vshape,  p_io)
       CALL p_bcast(XPOINT(i)%io%minalt,    p_io)
       CALL p_bcast(XPOINT(i)%io%maxalt,    p_io)
       CALL p_bcast(XPOINT(i)%io%ds(:),   p_io)
       CALL p_bcast(XPOINT(i)%io%de(:),   p_io)
       CALL p_bcast(XPOINT(i)%dt,     p_io)
       CALL p_bcast(XPOINT(i)%estop,  p_io)
       CALL p_bcast(XPOINT(i)%estart, p_io)
       ! XXREF_GP: ntrac, idt, flx, nidt, nflx -> set in ever_init_coupling
       ! lnow                              -> set in ever_global_start
       !
       ! shape 2 (Gaussian)
       CALL p_bcast(XPOINT(i)%io%meanalt,   p_io)
       CALL p_bcast(XPOINT(i)%io%vsigma,  p_io)       
       !
    END DO

    IF (p_parallel_io) THEN
       WRITE(*,*) ' ---> ',NPOINT,' RELEASE POINT(S) SPECIFIED !'
       WRITE(*,*) '------------------------------------------------------'
    END IF

    CALL end_message_bi(modstr,'INITIALISATION ',substr)

  END SUBROUTINE ever_initialize
  ! -------------------------------------------------------------------

  ! -------------------------------------------------------------------
  SUBROUTINE ever_init_memory

    IMPLICIT NONE

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'ever_init_memory'

#ifdef MESSYTENDENCY
    CALL mtend_register(my_handle, mtend_id_tracer)
#endif

  END SUBROUTINE ever_init_memory
  ! -------------------------------------------------------------------

  ! -------------------------------------------------------------------
  SUBROUTINE ever_init_coupling

    ! ECHAM5/MESSy
    USE messy_main_mpi_bi,           ONLY: p_parallel_io
    USE messy_main_blather_bi,       ONLY: error_bi, info_bi
    USE messy_main_tracer_mem_bi,    ONLY: GPTRSTR, gp_channel, xt, xtte &
                                        , ntrac_gp
    USE messy_main_tracer_tools_bi,  ONLY: tracer_halt
    USE messy_main_transform_bi,     ONLY: locate_in_decomp
    USE messy_main_channel_error_bi,       ONLY: channel_halt
    ! MESSy
    USE messy_main_tracer,           ONLY: get_tracer, tracer_error_str &
                                        , t_trinfo, full2base_sub &
                                        , STRLEN_FNAME
    USE messy_main_tools,            ONLY: strcrack, str2num
    USE messy_main_constants_mem,    ONLY: STRLEN_VLONG, STRLEN_ULONG
    USE messy_main_channel,          ONLY: get_channel_object, get_attribute & 
                                         , new_channel_object, new_channel &
                                         , new_attribute, STRLEN_OBJECT &
                                         , STRLEN_CHANNEL, get_channel_info
    USE messy_main_channel_error_bi, ONLY: channel_halt                                  
    USE messy_main_channel_bi,       ONLY: GP_3D_MID, SCALAR
    USE messy_main_tools,            ONLY: int2str

    IMPLICIT NONE

    INTRINSIC :: REAL, TRIM

    ! LOCAL
    INTEGER                          :: status
    CHARACTER(LEN=*), PARAMETER      :: substr = 'ever_init_coupling'
    INTEGER                          :: jt, i
    CHARACTER(LEN=STRLEN_VLONG)      :: errstr
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
    CHARACTER(LEN=3)                 :: istr = ''
    INTEGER                          :: j, idt
    CHARACTER(LEN=STRLEN_FNAME)      :: trfname = ''
    CHARACTER(LEN=STRLEN_OBJECT)     :: objname = ''
    CHARACTER(LEN=STRLEN_CHANNEL)    :: chaname = ''

    CALL start_message_bi(modstr, 'COUPLING INITIALIZATION', substr)

    ! INITIALIZE RELEASE POINTS
    IF (p_parallel_io) THEN
       WRITE(*,*) '------------------------------------------------------'
       WRITE(*,*) '                VOLCANIC ERUPTIONS                        '
       WRITE(*,*) '------------------------------------------------------'
    END IF

    point_loop: DO i=1, NPOINT
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
             IF (XPOINT(i)%ierr == 1) THEN
                write (*,*) ' longitude out of range: ', XPOINT(i)%io%lon, &
                     '; skipping POINT ', i
             ELSE IF (XPOINT(i)%ierr == 2) THEN
                write (*,*) ' latitude out of range: ',  XPOINT(i)%io%lat, &
                     '; skipping POINT ', i
             ELSE IF (XPOINT(i)%ierr == 555) THEN
                write (*,*) ' POINT not located in model domain: ' &
                     , XPOINT(i)%io%lon, XPOINT(i)%io%lat,         &
                     '; skipping POINT ', i
             ENDIF
          ENDIF
          CYCLE
       ENDIF

       ! DIAGNOSTIC OUTPUT
       IF (p_parallel_io) THEN
          WRITE(*,'(1a,i3.3,1x,a5,f9.4,1x,a5,f9.4,1x,a5,i3,1x,a5,i3,'//&
               &'1x,a5,i3,1x,a5,i3,1x,a5,i3)')  &
               '#',i                                                  &
               , ' LAT=',XPOINT(i)%io%lat, ' LON=',XPOINT(i)%io%lon   &
               , '  PE=',XPOINT(i)%pe,     '  JP=', XPOINT(i)%jp      &
               , 'JROW=',XPOINT(i)%jrow                               &
               , ' JGX=',XPOINT(i)%jgx,    ' JGY=',XPOINT(i)%jgy
       ENDIF
       !
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
       
    END DO point_loop


    ! DEFINE NEW CHANNEL
    CALL new_channel(status, modstr)
    CALL channel_halt(substr, status)

    point_loop2: DO i=1, NPOINT

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

    IF (p_parallel_io) THEN
       WRITE(*,*) '------------------------------------------------------'
       WRITE(*,*) '                GP TRACER COUPLING                    '
       WRITE(*,*) '------------------------------------------------------'
    END IF
    CALL init_trac_cpl(XXREF_GP, GPTRSTR)
    
    IF (p_parallel_io) THEN
       WRITE(*,*) '------------------------------------------------------'
       WRITE(*,*) '             GP EMISSION DIAGNOSTICS                  '
       WRITE(*,*) '------------------------------------------------------'
    END IF

    ALLOCATE(emis_gp(ntrac_gp))

    chaname = modstr//'_gp'

    point_loop3: DO i=1, NPOINT
    
       IF (XPOINT(i)%ierr /= 0) CYCLE

       trac_loop3: DO j=1, XXREF_GP(i)%ntrac
          ! GET TRCER ID
          idt   = XXREF_GP(i)%idt(j)
          IF (idt == 0) CYCLE        ! tracer does not exist

          ! SET TRACER ID, FIND NAME AND CREATE CHANNEL OBJECT NAME
          CALL get_tracer(status, GPTRSTR, idt, fullname=trfname)
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
               p3=emis_gp(idt)%ptr)
          IF (status == 0) THEN
             CALL info_bi('object '//TRIM(objname)//' in channel '//&
                  &TRIM(chaname)//' exists already', substr)
          ELSE
             CALL info_bi('creating object '//TRIM(objname)//&
                  &' in channel '//TRIM(chaname), substr)
             CALL new_channel_object(status, TRIM(chaname) &
                  , TRIM(objname) &
                  , p3=emis_gp(idt)%ptr, reprid=GP_3D_MID)
             CALL channel_halt(substr, status)
             CALL new_attribute(status, TRIM(chaname) &
                  , TRIM(objname) &
                  , 'long_name' &
                  , c='sum of point emissions'//&
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
      USE messy_main_constants_mem, ONLY: I8, pi
      USE messy_main_timer,         ONLY: time_step_len
      
      USE messy_gmxe,               ONLY: modes_gmxe=>cmodes
      USE messy_gmxe_mem,           ONLY: sigma_gmxe=>sigma, nmod
      USE messy_main_switch,        ONLY: USE_GMXE

      IMPLICIT NONE
      INTRINSIC :: MOD, REAL

      ! I/O
      ! x referene table
      TYPE(ERUPT_POINT_XREF), DIMENSION(:), INTENT(INOUT) :: XXREF
      CHARACTER(LEN=*),         INTENT(IN)    :: TRSTR ! tracer set

      ! LOCAL
      INTEGER :: i
      CHARACTER(LEN=2*STRLEN_MEDIUM+1), DIMENSION(:), POINTER :: &
           tracs => NULL(), pars => NULL(), masses => NULL()
      INTEGER :: zntrac, znmass, znaepr
      INTEGER :: j
      INTEGER :: counter 
      ! correction factor for flux to achieve exact total mass
      REAL(DP)    :: fcorr
      INTEGER(I8) :: iets  ! emission time in seconds
      INTEGER(I8) :: idts  ! time step length in seconds
      INTEGER(I8) :: irest
      INTEGER(I8) :: nsteps
      INTEGER(I8) :: jc, jt, mode
      
      REAL(DP)    :: convMtoN, diameter3, sigma_exp_ln, mass_r, sigma
      

      idts = INT(time_step_len, I8)

      point_loop: DO i=1, NPOINT

         IF (XPOINT(i)%ierr /= 0) CYCLE
         ! PARSE TRACER LIST
         CALL strcrack(XPOINT(i)%io%trlist, ';', tracs, zntrac,lempty=.TRUE.)
         CALL strcrack(XPOINT(i)%io%trmass, ';', masses, znmass,lempty=.TRUE.)
         CALL strcrack(XPOINT(i)%io%traepr, ';', pars, znaepr,lempty=.TRUE.)
         
         IF ((zntrac /= znmass) .OR. (zntrac /= znaepr)) THEN
            IF (p_parallel_io) THEN
               WRITE(*,*) 'WARNING: Tracer lists have differing size ... skipping ...'
            END IF
            !
            CYCLE
         END IF
      

         ! INITIALIZE
         XXREF(i)%ntrac = zntrac
         !
         ALLOCATE(XXREF(i)%idt(zntrac))
         XXREF(i)%idt(:) = 0
         ALLOCATE(XXREF(i)%flx(zntrac))
         XXREF(i)%flx(:) = 0.0_DP
         ALLOCATE(XXREF(i)%nidt(zntrac))
         XXREF(i)%nidt(:) = 0
         ALLOCATE(XXREF(i)%nflx(zntrac))
         XXREF(i)%nflx(:) = 0.0_DP         
         !
                 
         trac_loop: DO j=1, zntrac
            !
            ! GET EMISSION IN mol/s FOT TRACER
            CALL str2num(TRIM(masses(j)),mass_r)
            !
            ! CHECK IF TRACER EXISTS
            fullname = TRIM(tracs(j))
            CALL full2base_sub(status, TRIM(fullname) &
                 , basename, subname)
            CALL tracer_halt(substr, status)
            !
            CALL get_tracer(status, TRSTR, basename       &
                 , subname=subname, idx=XXREF(i)%idt(j), trinfo=trinfo)
            !

            ! CHECK STATUS (AND SET FLUX)
            status_how: IF (status /= 0) THEN
               errstr = tracer_error_str(status)
               IF (p_parallel_io) THEN
                  WRITE(*,*) '     -> *** WARNING FOR TRACER ',TRIM(fullname) &
                       ,': ',TRIM(errstr)
               END IF
               !
            ELSE
               !
               molarmass = trinfo%meta%cask_r(R_molarmass)
               !
               ! CALCULATE FLUX IN mol/s
               IF (molarmass <= ZERO_EPS) THEN   ! ZERO
                  XXREF(i)%flx(j) = 0.0_DP
               ELSE
                  XXREF(i)%flx(j) = (mass_r * 1.e9_dp    &
                       * (1.e3_dp/molarmass) / XPOINT(i)%dt)
                       
                  ! FLUX CORRECTION FOR TOTAL MASS
                  iets = INT(XPOINT(i)%dt, I8)
                  irest = MOD(iets, idts)
                  IF (irest == 0_i8 ) THEN
                     nsteps = (iets / idts)
                     fcorr = 1.0_dp
                  ELSE
                     nsteps = (iets / idts) + 1_i8
                     fcorr = XPOINT(i)%dt / &
                          REAL(nsteps * idts, dp)
                  ENDIF
                  XXREF(i)%flx(j) = XXREF(i)%flx(j) * fcorr
                  
               END IF
               !
               IF (p_parallel_io) THEN
                  WRITE(*,*) '     -> ',TRIM(fullname),           &
                       ' (id = ', XXREF(i)%idt(j),';',            &
                       ' molar mass = ',molarmass,' [g/mol];'
                  WRITE(*,*) '           flux = ',XXREF(i)%flx(j) &
                       ,'[mol/s])'
                  WRITE(*,*) '           (corrected with fctor ',fcorr,';'
                  WRITE(*,*) '            DT_emis = ',iets,' s ;' &
                       , ' DT = ',idts,' s ;'
                  WRITE(*,*) '            #steps = ',nsteps,')'
               END IF
               !
            END IF status_how
            !    
            !
            CALL get_tracer(status, TRSTR, 'N'       &
                , subname, idx=XXREF(i)%nidt(j), trinfo=trinfo)
                
            ! CHECK STATUS (AND SET NUMBER FLUX)
            status_how2: IF (status /= 0) THEN
               errstr = tracer_error_str(status)
               IF (p_parallel_io) THEN
                  WRITE(*,*) '     -> *** WARNING FOR TRACER ',TRIM(fullname) &
                       ,': ',TRIM(errstr)
               END IF
            ELSE
               aerpar_loop: DO jc=1, NPAR
                  IF (TRIM(pars(j)) /= TRIM(XAER(jc)%io%par_name)) CYCLE
                  !
                  ! CALCULATE FLUX IN N/s
                  IF (l_gmxe .AND. USE_GMXE) THEN
                     DO jt = 1,nmod
                        IF (TRIM(modes_gmxe(jt)) == TRIM(XAER(jc)%io%mode)) THEN
                           mode = jt
                           EXIT
                        END IF
                     END DO
                     sigma = sigma_gmxe(mode)
                  ELSE
                     sigma = XAER(jc)%io%sigma 
                  END IF
                  IF (sigma <= 0.0) CYCLE    
                  sigma_exp_ln = EXP(4.5*LOG(sigma) * LOG(sigma))
                  !
                  diameter3 = XAER(jc)%io%diameter * XAER(jc)%io%diameter &
                        * XAER(jc)%io%diameter
                  !
                  convMtoN = 6._dp / pi &
                        / (XAER(jc)%io%density * 1.e3_dp) &
                        / (diameter3 * sigma_exp_ln)
                  !
                  ! Set the number flux      
                  XXREF(i)%nflx(j) = (mass_r * 1.e9_dp    &
                        * convMtoN / XPOINT(i)%dt) * fcorr
                               
               END DO aerpar_loop
         
            END IF status_how2
            !
         END DO trac_loop
         !
      END DO point_loop

    END SUBROUTINE init_trac_cpl
    ! -------------------------------------------------------------------

  END SUBROUTINE ever_init_coupling
  ! -------------------------------------------------------------------

  ! -------------------------------------------------------------------
  SUBROUTINE ever_global_start

    USE messy_main_timer,     ONLY: YEAR, MONTH, DAY, HOUR, MINUTE, SECOND &
                                  , time_step_len
    USE messy_main_timer,     ONLY: gregor2julian
#ifdef _TFCORR
    USE messy_main_constants_mem, ONLY: OneDay
#endif

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

  END SUBROUTINE ever_global_start
  ! -------------------------------------------------------------------

  ! -------------------------------------------------------------------
  SUBROUTINE ever_physc

    ! ECHAM5/MESSy
    USE messy_main_mpi_bi,         ONLY: p_pe, p_parallel_io
#ifndef MESSYTENDENCY
    USE messy_main_data_bi,        ONLY: tm1_3d   &
                                       , qm1 => qm1_3d  &
                                       , qte => qte_3d
#else
    USE messy_main_tracer_mem_bi,  ONLY: ntrac_gp
#endif
#ifdef _TFCORR
    USE messy_main_tracer_mem_bi,  ONLY: qxtf, qxt  
    USE messy_main_data_bi,        ONLY: eps         
#endif

    ! _TFCORR or .NOT. MESSYTENDENCY
    USE messy_main_tracer_mem_bi,  ONLY: qxtm1, qxtte, ntrac_gp
    USE messy_main_grid_def_mem_bi,ONLY: jrow, kproma, nlev
    USE messy_main_grid_def_bi,    ONLY: grmass=>grmassdry &
                                       , zbound=>altitudei_msl &
                                       , grvol
    USE messy_main_data_bi,        ONLY: pbound=>pressi_3d
    USE messy_main_timer,          ONLY: ztmst=>time_step_len
    USE messy_main_tracer,         ONLY: get_tracer 

    USE messy_main_tools,          ONLY: nl_index

    USE messy_main_constants_mem,  ONLY: N_A, M_air

    IMPLICIT NONE

    INTRINSIC :: EXP, ASSOCIATED, MIN, MAX

    ! LOCAL
    INTEGER                             :: jt, i,j, ii
    INTEGER                             :: jp, jk, kmin, kmax, jc
    INTEGER                             :: idt, nidt
    !
    REAL(DP)                            :: theight! total height from minalt to maxalt
    REAL(dp),                 PARAMETER :: sqrt2=1.414213562373_dp
    REAL(dp)                            :: erf1, erf2, integral, frac &
                                         , high, low, amin, amax, amean, asigma

    REAL(DP), DIMENSION(_RI_X_ZN_(kproma,nlev,ntrac_gp)), TARGET :: zfxtte

    ! INITIALISE EMISSION DIAGNOSTIC
    DO jt=1, ntrac_gp
       IF (ASSOCIATED(emis_gp(jt)%ptr)) THEN
          emis_gp(jt)%ptr(_RI_XYZ__(:,jrow,:)) = 0.0_dp
       END IF
    END DO
    

    ! ADD EMISSION TENDENCY FOR ALL TRACERS
    !
    zfxtte(:,:,:) = 0.0_dp 
    !
    !
    point_loop: DO i=1, NPOINT
       !
       IF (XPOINT(i)%ierr /= 0)     CYCLE
       IF (XPOINT(i)%pe   /= p_pe)  CYCLE   ! LOCATION IS NOT ON THIS pe
       IF (XPOINT(i)%jrow /= jrow)  CYCLE   ! LOCATION IS NOT IN CURRENT ROW
       ! GET VECTOR INDEX
       jp    = XPOINT(i)%jp 

       IF (.NOT. XPOINT(i)%lnow)    CYCLE   ! EMISSION NOT AT THIS TIME
       !
       amin = MIN(XPOINT(i)%io%minalt, XPOINT(i)%io%maxalt) * 1000 ! in meters
       amax = MAX(XPOINT(i)%io%minalt, XPOINT(i)%io%maxalt) * 1000 ! in meters
       CALL nl_index(zbound(_RI_XYZ__(jp, jrow,1:nlev)), amin, kmax)
       CALL nl_index(zbound(_RI_XYZ__(jp, jrow,1:nlev)), amax, kmin)
       !
       amean = XPOINT(i)%io%meanalt * 1000
       asigma = XPOINT(i)%io%vsigma * 1000
       !        
       IF (XPOINT(i)%io%vshape == 2) THEN
          amean = XPOINT(i)%io%meanalt * 1000
          asigma = XPOINT(i)%io%vsigma * 1000
          erf2 = XERF((amax - amean)/asigma/sqrt2)
          erf1 = XERF((amin - amean)/asigma/sqrt2)
          integral = 0.5_DP * (erf2 - erf1)
       END IF
        
       ! GET TOTAL HEIGHT FROM MINALT TO MAXALT
       theight = amax - amin
       !
       trac_loop: DO j=1, XXREF_GP(i)%ntrac
          !
          ! GET TRACER ID AND TRACER NUMBER ID
          idt   = XXREF_GP(i)%idt(j)
          nidt  = XXREF_GP(i)%nidt(j)
          IF (idt == 0) CYCLE        ! tracer does not exist
          !
          DO jk=kmin, kmax
             !
#ifdef _TFCORR
             IF (l_tf_corr) THEN
                IF (XPOINT(i)%l2nd) THEN
                   zfxtte(_RI_X_ZN_(jp,jk,idt)) = &
                        - eps * qxt(_RI_X_ZN_(jp,jk,idt)) / ztmst
                END IF
             END IF
#endif
             ! 
             ! ##### GET FRACTION FOR EMISSION #####
             high = MIN(amax, zbound(_RI_XYZ__(jp, jrow,jk)))
             low = MAX(amin, zbound(_RI_XYZ__(jp, jrow,jk+1)))
             
             frac = 1.0
             IF (XPOINT(i)%io%vshape == 1 .AND. theight > 0) THEN
                frac = (high - low) / theight
             END IF
             IF (XPOINT(i)%io%vshape == 2 .AND. integral > 0) THEN
                erf2 = XERF((high - amean)/asigma/sqrt2)
                erf1 = XERF((low - amean)/asigma/sqrt2)
                frac = 0.5_DP * (erf2 - erf1) / integral    
             END IF             
             !
             !
             !
             ! ##### ADD EMISSION TO TRACER TENDENCY (mass) #####
             ! (mol/s) / (kg/(kg/mol)) * volume_fraction -> mol/mol/s
             zfxtte(_RI_X_ZN_(jp,jk,idt)) = zfxtte(_RI_X_ZN_(jp,jk,idt)) &
                 + XXREF_GP(i)%flx(j)                           &
                 / ( grmass(_RI_XYZ__(jp,jrow,jk)) * ( 1.E3_DP/M_air ) )   &
                 * frac

             ! (mol/s) * Na / m^3 -> molec/m^3/s
             emis_gp(idt)%ptr(_RI_XYZ__(jp,jrow,jk)) = &
                emis_gp(idt)%ptr(_RI_XYZ__(jp,jrow,jk)) + &
                XXREF_GP(i)%flx(j) * N_A / grvol(_RI_XYZ__(jp,jrow,jk)) &
                * frac
             !
#ifdef _TFCORR
             IF (l_tf_corr) THEN
                IF (XPOINT(i)%lfirst) THEN
                   qxtf(_RI_X_ZN_(jp,jk,idt)) = qxtm1(_RI_X_ZN_(jp,jk,idt))  &
                        + (qxtte(_RI_X_ZN_(jp,jk,idt)) + zfxtte(_RI_X_ZN_(jp,jk,idt))) & 
                        * ztmst
                END IF
             END IF
#endif
             !
             !
             !
             ! ##### ADD NUMBER CONCENTRATION FOR AEROSOL SPECIES #####
             ! GET TRACER NUMBER ID (only for aerosols)
             IF (nidt == 0) CYCLE        ! no aerosol tracer
             !
#ifdef _TFCORR
             IF (l_tf_corr) THEN
                IF (XPOINT(i)%l2nd) THEN
                   zfxtte(_RI_X_ZN_(jp,jk,nidt)) = &
                        - eps * qxt(_RI_X_ZN_(jp,jk,nidt)) / ztmst
                END IF
             END IF
#endif
             ! (N/s) / (cm3) * volume_fraction -> N/cm3/s
             zfxtte(_RI_X_ZN_(jp,jk,nidt)) = zfxtte(_RI_X_ZN_(jp,jk,nidt)) &
                 + XXREF_GP(i)%nflx(j) / grvol(_RI_XYZ__(jp,jrow,jk))      &
                 * frac  
             !    
#ifdef _TFCORR
             IF (l_tf_corr) THEN
                IF (XPOINT(i)%lfirst) THEN
                   qxtf(_RI_X_ZN_(jp,jk,nidt)) = qxtm1(_RI_X_ZN_(jp,jk,nidt))  &
                        + (qxtte(_RI_X_ZN_(jp,jk,nidt)) + zfxtte(_RI_X_ZN_(jp,jk,nidt))) & 
                        * ztmst
                END IF
             END IF
#endif
             !
          END DO
        END DO trac_loop
        !
    END DO point_loop

#ifndef MESSYTENDENCY
       qxtte(_RI_X_ZN_(1:kproma,:,:)) = & 
            qxtte(_RI_X_ZN_(1:kproma,:,:)) + & 
            zfxtte(_RI_X_ZN_(1:kproma,:,:))                            
#else
       CALL mtend_add_l(my_handle, mtend_id_tracer, pxt=zfxtte)                                                        
#endif   
    !

  END SUBROUTINE ever_physc
  ! -------------------------------------------------------------------

  ! -------------------------------------------------------------------
  SUBROUTINE ever_free_memory

    IMPLICIT NONE
    
    INTRINSIC :: ASSOCIATED

    ! LOCAL
    INTEGER :: i

    DO i=1, NPOINT
        IF (ASSOCIATED(XXREF_GP(i)%idt))       DEALLOCATE(XXREF_GP(i)%idt)
        IF (ASSOCIATED(XXREF_GP(i)%flx))       DEALLOCATE(XXREF_GP(i)%flx)
    END DO

    IF (ASSOCIATED(emis_gp)) THEN
       DEALLOCATE(emis_gp) ; NULLIFY(emis_gp)
    ENDIF

  END SUBROUTINE ever_free_memory
  ! -------------------------------------------------------------------

  ! -------------------------------------------------------------------
  SUBROUTINE ever_read_nml_cpl(status, iou)

    ! MESSy
    USE messy_main_tools, ONLY: read_nml_open, read_nml_check, read_nml_close
    ! ECHAM5/MESSy
    IMPLICIT NONE

    ! I/O
    INTEGER, INTENT(OUT) :: status     ! error status
    INTEGER, INTENT(IN)  :: iou        ! I/O unit

    ! (LOCAL) NAMELIST VARIABLES
    CHARACTER(LEN=*), PARAMETER :: substr = 'ever_read_nml_cpl'

    NAMELIST /CPL/ l_tf_corr, l_gmxe, AER, POINT

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

    status = 0  ! no ERROR

  END SUBROUTINE ever_read_nml_cpl
  ! -------------------------------------------------------------------

! **********************************************************************
END MODULE messy_ever_si
! **********************************************************************

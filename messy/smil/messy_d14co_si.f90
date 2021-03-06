#include "messy_main_ppd_bi.inc"
! **********************************************************************
MODULE messy_d14co_si
! **********************************************************************

! **********************************************************************
! INTERFACE TO ECHAM-5 FOR
! 14CO module for OH/STE diagnostic
! Version: see modver in messy_d14co.f90
!
! Author : Patrick Joeckel, MPICH, October  2002
!                                  November 2002
!                                  December 2002
!                                  July     2004
!
! References:
! (1) Patrick Joeckel, Carl A. M. Brenninkmeijer, Mark G. Lawrence,
!     Adriaan B. M. Jeuken, and Peter F. J. van Velthoven,
!     Evaluation of stratosphere - troposphere exchange and the hydroxyl
!     radical distribution in 3-dimensional global atmospheric models using
!     observations of cosmogenic 14CO,
!     J. Geophys. Res., 107(D20), 4446, doi:10.1029/2001JD001324, 2002.
! (2) Patrick Joeckel, and Carl A.M. Brenninkmeijer,
!     The seasonal cycle of cosmogenic 14CO at surface level: A solar
!     cycle adjusted, zonal-average climatology based on observations,
!     J. Geophys. Res., 107(D22), 4656, doi:10.1029/2001JD001104, 2002.
! (3) Patrick Joeckel,
!     Cosmogenic 14CO as tracer for atmospheric chemistry and transport,
!     Dissertation, Combined Faculties for the Natural Sciences and for
!     Mathematics  of the Rupertus Carola University of
!     Heidelberg, Germany, 2000.
!     (http://www.ub.uni-heidelberg.de/archiv/1426)
! (4) D. C. McCabe, T. Gierczak, R. K. Talukdar and A. R. Ravishankara,
!     Kinetics of the reaction OH + CO under atmospheric conditions,
!     Geophys. Res. Lett., 28, 3135-3138, 2001.
!
! **********************************************************************

  USE messy_main_blather_bi,            ONLY: start_message_bi, end_message_bi
  USE messy_main_channel,               ONLY: STRLEN_CHANNEL, STRLEN_OBJECT
  USE messy_d14co

  IMPLICIT NONE
  PRIVATE

  INTRINSIC :: TRIM, SIZE, COS, ASSOCIATED, ADJUSTL, MAX, NULL

  ! GLOBAL PARAMETERS
  !
  INTEGER, PARAMETER :: NMAXSETS       =  10 ! max. number of integration sets
  INTEGER, PARAMETER :: I_TP_DIAG      =  1  ! diagnosed tropopause
  INTEGER, PARAMETER :: I_TP_CLIM      =  2  ! climatological tropopause
  INTEGER, PARAMETER :: I_TP_CONST     =  3  ! const. pressure tropopause
  !
  ! TYPE
  INTEGER, PARAMETER :: IS_UNDEFINED = 0
  INTEGER, PARAMETER :: IS_GPOFFLINE = 1
  INTEGER, PARAMETER :: IS_GPTRACER = 2
  INTEGER, PARAMETER :: IS_GPCHANNEL = 4
!!# attila +
#ifdef ECHAM5
  INTEGER, PARAMETER :: IS_LGTRACER = 3
  INTEGER, PARAMETER :: IS_LGCHANNEL = 5
#endif
!!# attila -

  TYPE TYP_SWITCH
     LOGICAL :: gp_lint  = .FALSE.  ! GP integration on/off
     LOGICAL :: gp_ladjt = .FALSE.  ! GP adjust tendency (14CO = 14COs + 14COt)
     LOGICAL :: lg_lint  = .FALSE.  ! LG integration on/off
     LOGICAL :: lg_ladjt = .FALSE.  ! LG adjust tendency (14CO = 14COs + 14COt)
     LOGICAL :: lg_lmix  = .FALSE.  ! LG tracer mixing
  END TYPE TYP_SWITCH

  TYPE TYP_LINK
     CHARACTER(LEN=STRLEN_CHANNEL) :: channel  = ''
     CHARACTER(LEN=STRLEN_OBJECT ) :: object   = ''
  END TYPE TYP_LINK

  TYPE TYP_TROPOPAUSE
     INTEGER                       :: typ = 0
     TYPE(TYP_LINK)                :: link
     REAL(DP), DIMENSION(2)        :: clim = (/30000._DP, 21500._DP/)
     REAL(DP)                      :: const = 10000._DP
  END TYPE TYP_TROPOPAUSE

  ! CPL-NAMELIST VARIABLES
  !
!!#D attila +
#ifdef ECHAM5
  ! COUPLING TO LAGRANGIAN SCHEME
  CHARACTER(LEN=STRLEN_CHANNEL) :: C_LG_CHANNEL = ''
  CHARACTER(LEN=STRLEN_OBJECT ) :: C_LG_ILON    = ''
  CHARACTER(LEN=STRLEN_OBJECT ) :: C_LG_ILAT    = ''
  CHARACTER(LEN=STRLEN_OBJECT ) :: C_LG_PRESS   = ''
#endif
!!#D attila -
  !
  ! INTEGRATION SETS:
  ! - GLOBAL SWITCHES FOR INTEGRATION OF SET (lint) AND
  !   TENDENCY ADJUSTMENT (ladjt)
  TYPE(TYP_SWITCH),     DIMENSION(NMAXSETS), SAVE :: S_SWITCH
  ! - 14CO SOURCE SPECIFICATION
  TYPE(TYP_LINK),       DIMENSION(NMAXSETS), SAVE :: S_14CO
  ! - SPECIFICATION OF TROPOSPHERIC OH
  TYPE(TYP_LINK),       DIMENSION(NMAXSETS), SAVE :: S_OHt
  ! - SPECIFICATION OF STRATOSPHERIC OH
  TYPE(TYP_LINK),       DIMENSION(NMAXSETS), SAVE :: S_OHs
  ! - SPECIFICATION OF STE TROPOPAUSE
  TYPE(TYP_TROPOPAUSE), DIMENSION(NMAXSETS), SAVE :: S_TP_STE
  ! - SPECIFICATION OF TROPOPAUSE TO MERGE STRAT. AND TROP. OH
  TYPE(TYP_TROPOPAUSE), DIMENSION(NMAXSETS), SAVE :: S_TP_OH

  ! WORKSPACE OF INTEGRATION SET
  TYPE TYP_SET
     ! DIAGNOSTIC CHANNEL OBJECTS
     REAL(DP), DIMENSION(:,:),   POINTER :: ptp   => NULL() ! tropopause press.
     REAL(DP), DIMENSION(:,:),   POINTER :: ptpoh => NULL() ! tropopause press.
     ! POINTERS TO COUPLED CHANNEL OBJECTS
     REAL(DP), DIMENSION(:,:),   POINTER :: xtp_oh   => NULL() ! OH tropopause
     REAL(DP), DIMENSION(:,:),   POINTER :: xtp_ste  => NULL() ! STE tropoause
     !
     ! GRIDPOINT
     ! DIAGNOSTIC CHANNEL OBJECTS
     REAL(DP), DIMENSION(:,:,:), POINTER :: p14co_gp => NULL() ! 14CO prod.
     REAL(DP), DIMENSION(:,:,:), POINTER :: oh_gp    => NULL() ! OH
     REAL(DP), DIMENSION(:,:,:), POINTER :: l14co_gp => NULL() ! 14CO loss.
     ! TRACER INDICES
     INTEGER                             :: idt_gp_14CO  = 0
     INTEGER                             :: idt_gp_14COs = 0
     INTEGER                             :: idt_gp_14COt = 0
     !
!!#D attila +
#ifdef ECHAM5
     ! LAGRANGE
     ! DIAGNOSTIC CHANNEL OBJECTS
     REAL(DP), DIMENSION(:),     POINTER :: p14co_lg => NULL() ! 14CO prod.
     REAL(DP), DIMENSION(:),     POINTER :: oh_lg    => NULL() ! OH
     REAL(DP), DIMENSION(:),     POINTER :: l14co_lg => NULL() ! 14CO loss.
     !
     ! TRACER INDICES
     INTEGER                             :: idt_lg_14CO  = 0
     INTEGER                             :: idt_lg_14COs = 0
     INTEGER                             :: idt_lg_14COt = 0
#endif
!!#D attila -
     !
     ! ### RESOURCES ######################################################
     ! 14CO SOURCE --------------------------------------------------------
     ! - TYPE
     INTEGER                                :: q_is_what = IS_UNDEFINED
     ! - POINTER TO COUPLED CHANNEL OBJECT OR (PROCESSED) OFFLINE FIELD
     REAL(DP),  DIMENSION(:,:,:),   POINTER :: q_gp   => NULL()
!!#D attila +
#ifdef ECHAM5
     REAL(DP),  DIMENSION(:),       POINTER :: q_lg   => NULL()
#endif
!!#D attila -
     ! - SWICH FOR UNIT CONVERSION molec/(g s) -> mol/mol/s
     LOGICAL                                :: lconv_q = .FALSE.
     !
     ! STRATOSPHERIC OH --------------------------------------------------
     ! - TYPE
     INTEGER                                :: ohs_is_what = IS_UNDEFINED
     ! - POINTER TO COUPLED CHANNEL OBJECT OR (PROCESSED) OFFLINE FIELD
     REAL(DP),  DIMENSION(:,:,:),   POINTER :: ohs_gp => NULL()
!!#D attila +
#ifdef ECHAM5
     REAL(DP),  DIMENSION(:),       POINTER :: ohs_lg => NULL()
#endif
!!#D attila -
     ! - SWICH FOR UNIT CONVERSION mol/mol -> cm^-3
     LOGICAL                                :: lconv_OHs = .FALSE.
     ! - OFFLINE FIELD: NUMBER OF ASSOCIATED EVENT
     INTEGER                                :: i_ohs = 0   ! event number
     ! - TRACER: ID FOR OH-TRACER
     INTEGER                                :: idt_ohs = 0 ! tracer ID
     !
     ! TROPOSPHERIC OH --------------------------------------------------
     ! - TYPE
     INTEGER                                :: oht_is_what = IS_UNDEFINED
     ! - POINTER TO COUPLED CHANNEL OBJECT OR (PROCESSED) OFFLINE FIELD
     REAL(DP),  DIMENSION(:,:,:),   POINTER :: oht_gp => NULL()
!!#D attila +
#ifdef ECHAM5
     REAL(DP),  DIMENSION(:),       POINTER :: oht_lg => NULL()
#endif
!!#D attila -
     ! - SWICH FOR UNIT CONVERSION mol/mol -> cm^-3
     LOGICAL                                :: lconv_OHt = .FALSE.
     ! - OFFLINE FIELD: NUMBER OF ASSOCIATED EVENT
     INTEGER                                :: i_oht = 0   ! event number
     ! - TRACER: ID FOR OH-TRACER
     INTEGER                                :: idt_oht = 0 ! tracer ID
     !
     ! additional data required for offline -----------------------------
     ! 
     ! RAW DATA (NOT YET VERTICALLY INTERPOLATED !)
     REAL(DP), DIMENSION(:,:,:), POINTER   :: import_q   => NULL()
     REAL(DP), DIMENSION(:,:,:), POINTER   :: import_OHs => NULL()
     REAL(DP), DIMENSION(:,:,:), POINTER   :: import_OHt => NULL()
     ! GRID DATA
     REAL(DP), DIMENSION(:),     POINTER   :: paxis_q   => NULL()
     REAL(DP), DIMENSION(:),     POINTER   :: paxis_OHs => NULL()
     REAL(DP), DIMENSION(:),     POINTER   :: paxis_OHt => NULL()
  END TYPE TYP_SET

  TYPE(TYP_SET), DIMENSION(NMAXSETS), SAVE  :: XSET

!!#D attila +
#ifdef ECHAM5
  ! INFORMATION FROM LAGRANGIAN SCHEME
  LOGICAL                         :: LCALC_LG = .FALSE.
  REAL(DP), DIMENSION(:), POINTER :: ilon_lg  => NULL() ! LONGITUDE INDEX
  REAL(DP), DIMENSION(:), POINTER :: ilat_lg  => NULL() ! LATITUDE INDEX
  REAL(DP), DIMENSION(:), POINTER :: p_lg     => NULL() ! PRESSURE [Pa]
#endif
  !!#D attila -

  ! -------------------------------------------------------------

  ! PUBLIC ECHAM-5 INTERFACE ROUTINES TO 'messy_SUBMODELS'
  PUBLIC :: d14co_initialize     ! global initialisation of module
  PUBLIC :: d14co_new_tracer     ! define new tracers
  PUBLIC :: d14co_init_memory    ! allocate memory and define channel(s)
  PUBLIC :: d14co_init_coupling  ! check and initialize 'coupling'
  !                              ! to tracers and channel objects
  PUBLIC :: d14co_global_start   ! read offline fields on event
  PUBLIC :: d14co_physc          ! integrate one time step
  PUBLIC :: d14co_global_end     ! read offline fields on event
  PUBLIC :: d14co_free_memory    ! deallocate memory

  ! PRIVATE ECHAM-5 INTERFACE ROUTINES
  ! PRIVATE :: d14co_read_nml_cpl   ! initialize 'coupling' to online tracers
                                    ! ( /CPL/-namelist, /CPL_LG/-namelist )
  ! PRIVARE :: d14co_check_unit

CONTAINS

! ************************************************************************
! PUBLIC ECHAM-5 INTERFACE ROUTINES
! ************************************************************************

! ------------------------------------------------------------------------
  SUBROUTINE  d14co_initialize

    ! 14CO MODULE ROUTINE (ECHAM-5 INTERFACE)
    !
    ! INITIALIZATION OF GLOBAL VARIABLES FROM NAMELIST
    ! IN PARALLEL ENVIRONMENT
    ! INITIALIZATION OF 14CO SPECIFIC EVENTS FROM NAMELIST
    ! IN PARALLEL ENVIRONMENT
    !
    ! Author: Patrick Joeckel, MPICH, Feb 2002

    USE messy_main_mpi_bi,               ONLY: p_parallel_io, p_io, p_bcast
    USE messy_main_blather_bi,           ONLY: error_bi
    USE messy_main_tools,                ONLY: find_next_free_unit

    IMPLICIT NONE

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'd14co_initialize'
    INTEGER                     :: iou    ! I/O unit
    INTEGER                     :: status ! error status
    INTEGER                     :: n

    ! INITIALIZE 'COUPLING'
    IF (p_parallel_io) THEN
       iou = find_next_free_unit(100,200)
       CALL d14co_read_nml_cpl(status, iou)
       IF (status /= 0) CALL error_bi(' ', substr)
    END IF
    ! BROADCAST RESULTS
!!#D attila +
#ifdef ECHAM5
    CALL p_bcast(C_LG_CHANNEL, p_io)
    CALL p_bcast(C_LG_ILON, p_io)
    CALL p_bcast(C_LG_ILAT, p_io)
    CALL p_bcast(C_LG_PRESS, p_io)
#endif
!!#D attila -

    set_loop: DO n = 1, NMAXSETS

       ! S_SWITCH
       CALL p_bcast(S_SWITCH(n)%gp_lint, p_io)
       CALL p_bcast(S_SWITCH(n)%gp_ladjt, p_io)
       CALL p_bcast(S_SWITCH(n)%lg_lint, p_io)
       CALL p_bcast(S_SWITCH(n)%lg_ladjt, p_io)
       CALL p_bcast(S_SWITCH(n)%lg_lmix, p_io)

       ! S_14CO
       CALL p_bcast(S_14CO(n)%channel, p_io)
       CALL p_bcast(S_14CO(n)%object, p_io)
       ! S_OHt
       CALL p_bcast(S_OHt(n)%channel, p_io)
       CALL p_bcast(S_OHt(n)%object, p_io)
       ! S_OHs
       CALL p_bcast(S_OHs(n)%channel, p_io)
       CALL p_bcast(S_OHs(n)%object, p_io)

       ! S_TP_STE
       CALL p_bcast(S_TP_STE(n)%typ, p_io)
       CALL p_bcast(S_TP_STE(n)%link%channel, p_io)
       CALL p_bcast(S_TP_STE(n)%link%object, p_io)
       CALL p_bcast(S_TP_STE(n)%clim(1), p_io)
       CALL p_bcast(S_TP_STE(n)%clim(2), p_io)
       CALL p_bcast(S_TP_STE(n)%const, p_io)
       ! S_TP_OH
       CALL p_bcast(S_TP_OH(n)%typ, p_io)
       CALL p_bcast(S_TP_OH(n)%link%channel, p_io)
       CALL p_bcast(S_TP_OH(n)%link%object, p_io)
       CALL p_bcast(S_TP_OH(n)%clim(1), p_io)
       CALL p_bcast(S_TP_OH(n)%clim(2), p_io)
       CALL p_bcast(S_TP_OH(n)%const, p_io)

    END DO set_loop
       
!!#D attila +
#ifdef ECHAM5
    ! LG CALCULATION ?
    DO n = 1, NMAXSETS
       LCALC_LG = LCALC_LG .OR. S_SWITCH(n)%lg_lint
    END DO
#endif
!!#D attila -

  END SUBROUTINE d14co_initialize
! ------------------------------------------------------------------------

! ------------------------------------------------------------------------
  SUBROUTINE d14co_new_tracer

    ! 14CO MODULE ROUTINE (ECHAM-5 INTERFACE, PRIVATE)
    !
    ! defines 14CO specific tracers
    !
    ! Author: Patrick Joeckel, MPICH, Feb 2002

    USE messy_main_mpi_bi,          ONLY: p_parallel_io
    USE messy_main_tracer_mem_bi,   ONLY: GPTRSTR
    USE messy_main_tracer_tools_bi, ONLY: tracer_halt
    USE messy_main_tools,           ONLY: int2str
    USE messy_main_tracer,          ONLY: new_tracer, set_tracer &
                                        , R_molarmass
!!#D attila +
#ifdef ECHAM5
    USE messy_main_tracer,        ONLY: ON, OFF, I_mix
    USE messy_main_tracer_mem_bi, ONLY: LGTRSTR
#endif
!!#D attila -                         
    
    IMPLICIT NONE

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'd14co_new_tracer'
    INTEGER          :: status
    INTEGER          :: n
    CHARACTER(LEN=2) :: chnn
!!#D attila +
#ifdef ECHAM5
    INTEGER          :: imix
#endif
!!#D attila -

    CALL start_message_bi(modstr,'TRACER REQUEST',substr)

    set_loop: DO n=1, NMAXSETS

       CALL int2str(chnn, n, '0', 'X')

       ! GRIDPOINT
       IF (S_SWITCH(n)%gp_lint) THEN

          IF (p_parallel_io) THEN
             WRITE(*,*) ' SET No. ',n,'(GRIDPOINT)'
          END IF

          IF (p_parallel_io) WRITE(*,*) ' ... CO_14C_'//chnn
          CALL new_tracer(status, GPTRSTR, 'CO', modstr        &
               , subname = '14C_'//chnn, unit='mol/mol'        &
               , idx = XSET(n)%idt_gp_14CO)
          CALL tracer_halt(substr, status)
          CALL set_tracer(status, GPTRSTR, XSET(n)%idt_gp_14CO &
               , R_molarmass, r=30._DP)
          CALL tracer_halt(substr, status)
          
          IF (p_parallel_io) WRITE(*,*) ' ... CO_14Cs_'//chnn
          CALL new_tracer(status, GPTRSTR, 'CO', modstr         &
               , subname = '14Cs_'//chnn, unit='mol/mol'        &
               , idx = XSET(n)%idt_gp_14COs)
          CALL tracer_halt(substr, status)
          CALL set_tracer(status, GPTRSTR, XSET(n)%idt_gp_14COs &
               , R_molarmass, r=30._DP)
          CALL tracer_halt(substr, status)
          
          IF (p_parallel_io) WRITE(*,*) ' ... CO_14Ct_'//chnn
          CALL new_tracer(status, GPTRSTR, 'CO', modstr         &
               , subname = '14Ct_'//chnn, unit='mol/mol'        & 
               , idx = XSET(n)%idt_gp_14COt)
          CALL tracer_halt(substr, status)
          CALL set_tracer(status, GPTRSTR, XSET(n)%idt_gp_14COt &
               , R_molarmass, r=30._DP)
          CALL tracer_halt(substr, status)

       END IF

!!#D attila +
#ifdef ECHAM5
       ! LAGRANGE
       IF (S_SWITCH(n)%lg_lint) THEN

          IF (S_SWITCH(n)%lg_lmix) THEN
             imix = ON
             IF (p_parallel_io) THEN
                WRITE(*,*) ' SET No. ',n,'(LAGRANGE (WITH MIXING))'
             END IF
          ELSE
             imix = OFF
             IF (p_parallel_io) THEN
                WRITE(*,*) ' SET No. ',n,'(LAGRANGE (NO MIXING))'
             END IF
          END IF

          IF (p_parallel_io) WRITE(*,*) ' ... CO_14C_'//chnn
          CALL new_tracer(status, LGTRSTR, 'CO', modstr        &
               , subname = '14C_'//chnn, unit='mol/mol'        &
               , idx = XSET(n)%idt_lg_14CO                     &
               )
          CALL tracer_halt(substr, status)
          CALL set_tracer(status, LGTRSTR, XSET(n)%idt_lg_14CO &
               , R_molarmass, r=30._DP)
          CALL tracer_halt(substr, status)
          CALL set_tracer(status, LGTRSTR, XSET(n)%idt_lg_14CO &
               , I_mix, i=imix)
          CALL tracer_halt(substr, status)
          
          IF (p_parallel_io) WRITE(*,*) ' ... CO_14Cs_'//chnn
          CALL new_tracer(status, LGTRSTR, 'CO', modstr       &
               , subname = '14Cs_'//chnn, unit='mol/mol'      &
               , idx = XSET(n)%idt_lg_14COs                   &
               )
          CALL tracer_halt(substr, status)
          CALL set_tracer(status, LGTRSTR, XSET(n)%idt_lg_14COs &
               , R_molarmass, r=30._DP)
          CALL tracer_halt(substr, status)
          CALL set_tracer(status, LGTRSTR, XSET(n)%idt_lg_14COs &
               , I_mix, i=imix)
          CALL tracer_halt(substr, status)
          
          IF (p_parallel_io) WRITE(*,*) ' ... CO_14Ct_'//chnn
          CALL new_tracer(status, LGTRSTR, 'CO', modstr       &
               , subname = '14Ct_'//chnn, unit='mol/mol'      & 
               , idx = XSET(n)%idt_lg_14COt                   &
               )
          CALL tracer_halt(substr, status)
          CALL set_tracer(status, LGTRSTR, XSET(n)%idt_lg_14COt &
               , R_molarmass, r=30._DP)
          CALL tracer_halt(substr, status)
          CALL set_tracer(status, LGTRSTR, XSET(n)%idt_lg_14COt &
               , I_mix, i=imix)
          CALL tracer_halt(substr, status)

       END IF
#endif
!!#D attila -
    END DO set_loop

    CALL end_message_bi(modstr,'TRACER REQUEST',substr)

  END SUBROUTINE d14co_new_tracer
! ------------------------------------------------------------------------

! ------------------------------------------------------------------------
  SUBROUTINE d14co_init_memory

    ! 14CO MODULE ROUTINE (ECHAM-5 INTERFACE, PRIVATE)
    !
    ! define 14CO specific channel(s) and allocate memory for
    ! global fields
    !
    ! Author: Patrick Joeckel, MPICH, Feb 2002

    USE messy_main_mpi_bi,           ONLY: p_parallel_io
    USE messy_main_channel_error_bi, ONLY: channel_halt
    USE messy_main_channel_bi,       ONLY: GP_3D_MID, GP_2D_HORIZONTAL
!!#D attila +
#ifdef ECHAM5
    USE messy_main_channel_bi,       ONLY: LG_ATTILA
#endif
!!#D attila -
    USE messy_main_channel, ONLY: new_channel, new_channel_object &
                                , new_attribute
    USE messy_main_tools,   ONLY: int2str

    IMPLICIT NONE

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'd14co_init_memory'
    INTEGER                     :: status
    INTEGER                     :: n
    CHARACTER(LEN=2)            :: chnn
    LOGICAL                     :: lfirst_gp
!!#D attila +
#ifdef ECHAM5
    LOGICAL                     :: lfirst_lg
#endif
!!#D attila -

    CALL start_message_bi(modstr,'CHANNEL DEFINITION',substr)

    ! INIT
    lfirst_gp   = .TRUE.
!!#D attila +
#ifdef ECHAM5
    lfirst_lg   = .TRUE.
#endif
!!#D attila -
    set_loop: DO n=1, NMAXSETS

       CALL int2str(chnn, n, '0', 'X')

       ! SHARED
       IF (.NOT. (S_SWITCH(n)%gp_lint .OR. S_SWITCH(n)%lg_lint)) CYCLE

       IF (lfirst_gp) THEN
          ! DEFINE NEW CHANNEL
          CALL new_channel(status, modstr//'_gp')
          CALL channel_halt(substr, status)
          lfirst_gp = .FALSE.
       END IF

       IF (p_parallel_io) THEN
          WRITE(*,*) ' SET No. ',n,'(GP AND LG)'
       END IF

       ! used tropopause pressure (OUTPUT)
       IF (p_parallel_io) WRITE(*,*) ' ... PTP_'//chnn
       CALL new_channel_object(status, modstr//'_gp' &
            , 'PTP_'//chnn &
            , p2=XSET(n)%ptp &
            , reprid = GP_2D_HORIZONTAL )
       CALL channel_halt(substr, status)
       CALL new_attribute(status, modstr//'_gp', 'PTP_'//chnn &
            , 'long_name', c='tropopause pressure')
       CALL channel_halt(substr, status)
       CALL new_attribute(status, modstr//'_gp', 'PTP_'//chnn &
            , 'units', c='Pa')
       CALL channel_halt(substr, status)

       IF (p_parallel_io) WRITE(*,*) ' ... PTPOH_'//chnn
       CALL new_channel_object(status, modstr//'_gp' &
            , 'PTPOH_'//chnn   &
            , p2=XSET(n)%ptpoh &
            , reprid = GP_2D_HORIZONTAL )
       CALL channel_halt(substr, status)
       CALL new_attribute(status, modstr//'_gp', 'PTPOH_'//chnn &
            , 'long_name', c='tropopause pressure for OH')
       CALL channel_halt(substr, status)
       CALL new_attribute(status, modstr//'_gp', 'PTPOH_'//chnn &
            , 'units', c='Pa')
       CALL channel_halt(substr, status)

       ! GRIDPOINT
       IF (S_SWITCH(n)%gp_lint) THEN

          IF (p_parallel_io) THEN
             WRITE(*,*) ' SET No. ',n,'(GRIDPOINT)'
          END IF

          ! 14CO production rate (OUTPUT)
          IF (p_parallel_io) WRITE(*,*) ' ... P14CO_'//chnn
          CALL new_channel_object(status, modstr//'_gp' &
               , 'P14CO_'//chnn   &
               , p3=XSET(n)%p14co_gp &
               , reprid = GP_3D_MID )
          CALL channel_halt(substr, status)
          CALL new_attribute(status, modstr//'_gp', 'P14CO_'//chnn &
               , 'long_name', c='14CO production rate')
          CALL channel_halt(substr, status)
          CALL new_attribute(status, modstr//'_gp', 'P14CO_'//chnn &
               , 'units', c='mol/mol/s')
          CALL channel_halt(substr, status)
          
          ! 14CO loss rate (OUTPUT)
          IF (p_parallel_io) WRITE(*,*) ' ... L14CO_'//chnn
          CALL new_channel_object(status, modstr//'_gp' &
               , 'L14CO_'//chnn   &
               , p3=XSET(n)%l14co_gp &
               , reprid = GP_3D_MID )
          CALL channel_halt(substr, status)
          CALL new_attribute(status, modstr//'_gp', 'L14CO_'//chnn &
               , 'long_name', c='14CO loss rate')
          CALL channel_halt(substr, status)
          CALL new_attribute(status, modstr//'_gp', 'L14CO_'//chnn &
               , 'units', c='1/s')
          CALL channel_halt(substr, status)

          ! OH distribution (OUTPUT)
          IF (p_parallel_io) WRITE(*,*) ' ... OH_'//chnn
          CALL new_channel_object(status, modstr//'_gp' &
               , 'OH_'//chnn   &
               , p3=XSET(n)%oh_gp &
               , reprid = GP_3D_MID )
          CALL channel_halt(substr, status)
          CALL new_attribute(status, modstr//'_gp', 'OH_'//chnn &
               , 'long_name', c='OH')
          CALL channel_halt(substr, status)
          CALL new_attribute(status, modstr//'_gp', 'OH_'//chnn &
               , 'units', c='cm^(-3)')
          CALL channel_halt(substr, status)
       END IF

!!#D attila +
#ifdef ECHAM5
       ! LAGRANGE
       IF (S_SWITCH(n)%lg_lint) THEN

          IF (lfirst_lg) THEN
             ! DEFINE NEW CHANNEL
             CALL new_channel(status, modstr//'_lg')
             CALL channel_halt(substr, status)
             lfirst_lg = .FALSE.
          END IF

          IF (p_parallel_io) THEN
             WRITE(*,*) ' SET No. ',n,'(LAGRANGE)'
          END IF

          ! 14CO production rate (OUTPUT)
          IF (p_parallel_io) WRITE(*,*) ' ... P14CO_'//chnn
          CALL new_channel_object(status, modstr//'_lg' &
               , 'P14CO_'//chnn   &
               , p1=XSET(n)%p14co_lg &
               , reprid = LG_ATTILA )
          CALL channel_halt(substr, status)
          CALL new_attribute(status, modstr//'_lg', 'P14CO_'//chnn &
               , 'long_name', c='14CO production rate')
          CALL channel_halt(substr, status)
          CALL new_attribute(status, modstr//'_lg', 'P14CO_'//chnn &
               , 'units', c='mol/mol/s')
          CALL channel_halt(substr, status)
          
          ! 14CO loss rate (OUTPUT)
          IF (p_parallel_io) WRITE(*,*) ' ... L14CO_'//chnn
          CALL new_channel_object(status, modstr//'_lg' &
               , 'L14CO_'//chnn   &
               , p1=XSET(n)%l14co_lg &
               , reprid = LG_ATTILA )
          CALL channel_halt(substr, status)
          CALL new_attribute(status, modstr//'_lg', 'L14CO_'//chnn &
               , 'long_name', c='14CO loss rate')
          CALL channel_halt(substr, status)
          CALL new_attribute(status, modstr//'_lg', 'L14CO_'//chnn &
               , 'units', c='1/s')
          CALL channel_halt(substr, status)

          ! OH distribution (OUTPUT)
          IF (p_parallel_io) WRITE(*,*) ' ... OH_'//chnn
          CALL new_channel_object(status, modstr//'_lg' &
               , 'OH_'//chnn   &
               , p1=XSET(n)%oh_lg &
               , reprid = LG_ATTILA )
          CALL channel_halt(substr, status)
          CALL new_attribute(status, modstr//'_lg', 'OH_'//chnn &
               , 'long_name', c='OH')
          CALL channel_halt(substr, status)
          CALL new_attribute(status, modstr//'_lg', 'OH_'//chnn &
               , 'units', c='cm^(-3)')
          CALL channel_halt(substr, status)

       END IF
#endif
!!#D attila -

    END DO set_loop

    CALL end_message_bi(modstr,'CHANNEL DEFINITION',substr)

  END SUBROUTINE d14co_init_memory
! ------------------------------------------------------------------------

! ----------------------------------------------------------------------
  SUBROUTINE d14co_init_coupling

    ! 14CO MODULE ROUTINE (ECHAM-5 INTERFACE, PUBLIC)
    !
    ! initialize 'coupling' to online tracers/channels
    !
    ! Author: Patrick Joeckel, MPICH, Jan 2003

    USE messy_main_mpi_bi,           ONLY: p_parallel_io
    USE messy_main_tracer_mem_bi,    ONLY: GPTRSTR, gp_channel
    USE messy_main_tracer_tools_bi,  ONLY: tracer_halt
    USE messy_main_grid_def_mem_bi,  ONLY: nproma, ngpblks, nlev
    USE messy_main_channel_error_bi, ONLY: channel_halt
    USE messy_main_channel_bi,       ONLY: GP_3D_MID
!!#D attila +
#ifdef ECHAM5
    USE messy_main_tracer_mem_bi,    ONLY: LGTRSTR, NCELL, lg_channel
    USE messy_main_channel_bi,       ONLY: LG_ATTILA
#endif
!!#D attila -
    USE messy_main_blather_bi,       ONLY: error_bi
    USE messy_main_tracer,           ONLY: get_tracer, full2base_sub
    USE messy_main_channel,          ONLY: get_channel_object          &
                                         , get_channel_object_info     &
                                         , get_channel_object_dimvalue &
                                         , get_attribute               &
                                         , new_channel_object
    USE messy_main_constants_mem,    ONLY: STRLEN_MEDIUM, STRLEN_ULONG

    IMPLICIT NONE

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER       :: substr = 'd14co_init_coupling'
    INTEGER                           :: n
    INTEGER                           :: status ! status flag
    !
    CHARACTER(LEN=STRLEN_MEDIUM)      :: basename = ''
    CHARACTER(LEN=STRLEN_MEDIUM)      :: subname  = ''
!!#D attila +
#ifdef ECHAM5
    LOGICAL                           :: l_lg = .FALSE.  ! LAGRANGE ?
#endif
!!#D attila -
    CHARACTER(LEN=STRLEN_ULONG)       :: att_unit
    INTEGER                           :: reprid

    CALL start_message_bi(modstr,'COUPLING INITIALIZATION',substr)

    ! #######################################################################
    ! ### PART I: INTEGRATION SETS ##########################################
    ! #######################################################################
    IF (p_parallel_io) THEN
       WRITE(*,*) '================= INTEGRATION SETS ===================='
    END IF

    set_loop: DO n=1, NMAXSETS

!!#D attila +
#ifdef ECHAM5
       l_lg = l_lg .OR. S_SWITCH(n)%lg_lint
#endif
!!#D attila -

       IF (p_parallel_io) THEN
          WRITE(*,*) 'NUMBER ',n
          IF (S_SWITCH(n)%gp_lint) THEN
             WRITE(*,*) ' GP --> ON'
             IF (S_SWITCH(n)%gp_ladjt) THEN
                WRITE(*,*) ' ... WITH TENDENCY ADJUSTMENT'
             ELSE
                WRITE(*,*) ' ... WITHOUT TENDENCY ADJUSTMENT'
             END IF
          ELSE
             WRITE(*,*) ' GP --> OFF'
          END IF

!!#D attila +
#ifdef ECHAM5
          IF (S_SWITCH(n)%lg_lint) THEN
             WRITE(*,*) ' LG --> ON'
             IF (S_SWITCH(n)%lg_ladjt) THEN
                WRITE(*,*) ' ... WITH TENDENCY ADJUSTMENT'
             ELSE
                WRITE(*,*) ' ... WITHOUT TENDENCY ADJUSTMENT'
             END IF
          ELSE
             WRITE(*,*) ' LG --> OFF'
          END IF
#endif
!!#D attila -
       END IF

       IF (.NOT. (S_SWITCH(n)%gp_lint .OR. S_SWITCH(n)%lg_lint)) THEN
          IF (p_parallel_io) THEN
          WRITE(*,*) '-------------------------------------------------------'
          END IF
          CYCLE
       END IF
       
       ! SOURCE OF 14CO --------------------------------------------------
       IF (p_parallel_io) THEN
          WRITE(*,*) 'INITIALIZING 14CO SOURCE DISTRIBUTION ...'
       END IF
       !
       !
       ! CHANNEL OBJECT
       !
       IF (p_parallel_io) THEN
          WRITE(*,*) ' TYPE                      : ONLINE (CHANNEL OBJECT)'
          WRITE(*,*) '   channel                 : ',TRIM(S_14CO(n)%channel)
          WRITE(*,*) '   object                  : ',TRIM(S_14CO(n)%object)
       END IF
       CALL get_channel_object_info(status &
            , TRIM(S_14CO(n)%channel), TRIM(S_14CO(n)%object) &
            , reprid=reprid )
       CALL channel_halt(substr, status)
       IF (reprid == GP_3D_MID) THEN
          CALL get_channel_object(status  &
               , TRIM(S_14CO(n)%channel), TRIM(S_14CO(n)%object) &
               , p3=XSET(n)%q_gp )
          CALL channel_halt(substr, status)
          XSET(n)%q_is_what = IS_GPCHANNEL
!!#D attila +
#ifdef ECHAM5
       ELSEIF (reprid == LG_ATTILA) THEN
          CALL get_channel_object(status  &
               , TRIM(S_14CO(n)%channel), TRIM(S_14CO(n)%object) &
               , p1=XSET(n)%q_lg )
          CALL channel_halt(substr, status)
          XSET(n)%q_is_what = IS_LGCHANNEL
#endif
!!#D attila -
       ELSE
          IF (TRIM(ADJUSTL(S_14CO(n)%channel)) == 'import_grid') THEN
             XSET(n)%q_is_what = IS_GPOFFLINE
             CALL get_channel_object(status                         &
                  , TRIM(S_14CO(n)%channel), TRIM(S_14CO(n)%object) &
                  , p3=XSET(n)%import_q)
             CALL channel_halt(substr, status)
             ! get pressure levels from corresponding dimension variable
             CALL get_channel_object_dimvalue(status                &
                  ,TRIM(S_14CO(n)%channel), TRIM(S_14CO(n)%object)  &
                  , data=XSET(n)%paxis_q , axis='Z')
             CALL channel_halt(substr, status)
             !
             CALL new_channel_object(status          &
                  , modstr//'_gp', TRIM(S_14CO(n)%object)   &
                  , reprid=GP_3D_MID, p3=XSET(n)%q_gp)
             CALL channel_halt(substr, status)
          ELSE
             CALL error_bi(' representation not supported !', substr)
          END IF
       END IF

       ! CHECK UNIT:
       ! potentially required conversion
       CALL get_attribute(status  &
            , TRIM(S_14CO(n)%channel), TRIM(S_14CO(n)%object) &
            , 'units', c=att_unit)
       CALL channel_halt(substr, status)
       !
       CALL d14co_check_unit(substr, 'CS', TRIM(att_unit), XSET(n)%lconv_q)
       !
       IF (p_parallel_io) &
            WRITE(*,*) '   molec/(g s) -> mol/mol/s: ',XSET(n)%lconv_q

!!#D attila +
#ifdef ECHAM5
       SELECT CASE (XSET(n)%q_is_what)
       CASE (IS_GPOFFLINE)
          IF (S_SWITCH(n)%lg_lint) &
               ALLOCATE(XSET(n)%q_lg(NCELL))
       CASE (IS_GPTRACER)
          ! NOT IMPLEMENTED
       CASE (IS_GPCHANNEL)
          IF (S_SWITCH(n)%lg_lint) &
               ALLOCATE(XSET(n)%q_lg(NCELL))
       CASE (IS_LGTRACER)
          ! NOT IMPLEMENTED
       CASE (IS_LGCHANNEL)
          IF (S_SWITCH(n)%gp_lint) &
               ALLOCATE(XSET(n)%q_gp(nproma,nlev,ngpblks))
       CASE DEFAULT
          ! ERROR
       END SELECT
#endif
!!#D attila -
       ! ------------------------------------------------------------------

       ! STRATOSPHERIC OH -------------------------------------------------
       IF (p_parallel_io) THEN
          WRITE(*,*) 'INITIALIZING STRATOSPHERIC OH ...'
       END IF
       ! EVENT, GP-TRACER, LG-TRACER, OR CHANNEL OBJECT
       SELECT CASE(TRIM(ADJUSTL(S_OHs(n)%channel)))
         !
       CASE (gp_channel)
          !
          ! GP-TRACER
          !
          XSET(n)%ohs_is_what = IS_GPTRACER
          ! ... CHECK FOR SUBNAME
          CALL full2base_sub(status, TRIM(S_OHs(n)%object) &
               , basename, subname)
          CALL tracer_halt(substr, status)
          IF (p_parallel_io) THEN
             WRITE(*,*) ' TYPE                      : ONLINE (GP-TRACER)'
             WRITE(*,*) '   basename                : ',TRIM(basename)
             WRITE(*,*) '   subname                 : ',TRIM(subname)
          END IF
          CALL get_tracer(status, GPTRSTR, basename &
               , subname=subname                    &
               , idx=XSET(n)%idt_OHs)
          CALL tracer_halt(substr, status)
          !
!!#D attila +
#ifdef ECHAM5
       CASE (lg_channel)
          !
          ! LG-TRACER
          !
          XSET(n)%ohs_is_what = IS_LGTRACER
          ! ... CHECK FOR SUBNAME
          CALL full2base_sub(status, TRIM(S_OHs(n)%object) &
               , basename, subname)
          CALL tracer_halt(substr, status)
          IF (p_parallel_io) THEN
             WRITE(*,*) ' TYPE                      : ONLINE (LG-TRACER)'
             WRITE(*,*) '   basename                : ',TRIM(basename)
             WRITE(*,*) '   subname                 : ',TRIM(subname)
          END IF
          CALL get_tracer(status, LGTRSTR, basename &
               , subname=subname                  &
               , idx=XSET(n)%idt_OHs)
          CALL tracer_halt(substr, status)
          !
          CALL get_channel_object(status  &
               , TRIM(S_OHs(n)%channel), TRIM(S_OHs(n)%object) &
               , p1=XSET(n)%ohs_lg ) ! set LG-pointer 
#endif
!!#D attila -
          !
       CASE DEFAULT
          !
          ! CHANNEL OBJECT
          !
          IF (p_parallel_io) THEN
             WRITE(*,*) ' TYPE                      : ONLINE (CHANNEL OBJECT)'
             WRITE(*,*) '   channel                 : ',TRIM(S_OHs(n)%channel)
             WRITE(*,*) '   object                  : ',TRIM(S_OHs(n)%object)
          END IF
          CALL get_channel_object_info(status                  &
               , TRIM(S_OHs(n)%channel), TRIM(S_OHs(n)%object) &
               , reprid=reprid )
          CALL channel_halt(substr, status)
          IF (reprid == GP_3D_MID) THEN
             CALL get_channel_object(status                       &
                  , TRIM(S_OHs(n)%channel), TRIM(S_OHs(n)%object) &
                  , p3=XSET(n)%ohs_gp )
             CALL channel_halt(substr, status)
             XSET(n)%ohs_is_what = IS_GPCHANNEL
!!#D attila +
#ifdef ECHAM5
          ELSEIF (reprid == LG_ATTILA) THEN
             CALL get_channel_object(status                       &
                  , TRIM(S_OHs(n)%channel), TRIM(S_OHs(n)%object) &
                  , p1=XSET(n)%ohs_lg )
             CALL channel_halt(substr, status)
             XSET(n)%ohs_is_what = IS_LGCHANNEL
#endif
!!#D attila -
          ELSE
             IF (TRIM(ADJUSTL(S_OHs(n)%channel)) == 'import_grid') THEN
                XSET(n)%ohs_is_what = IS_GPOFFLINE
                CALL get_channel_object(status                       &
                     , TRIM(S_OHs(n)%channel), TRIM(S_OHs(n)%object) &
                     , p3=XSET(n)%import_OHs)
                CALL channel_halt(substr, status)
                ! get pressure levels via dimension variable
                CALL get_channel_object_dimvalue(status              &
                     ,TRIM(S_OHs(n)%channel), TRIM(S_OHs(n)%object)  &
                     , data=XSET(n)%paxis_OHs , axis ='Z')
                CALL channel_halt(substr, status)
                !
                CALL new_channel_object(status, modstr, TRIM(S_OHs(n)%object) &
                     , reprid=GP_3D_MID, p3=XSET(n)%OHs_gp)
                CALL channel_halt(substr, status)
             ELSE
                CALL error_bi(' representation not supported !', substr)
             END IF
          END IF
       END SELECT

       ! CHECK UNIT:
       ! potentially required conversion
       CALL get_attribute(status  &
            , TRIM(S_OHs(n)%channel), TRIM(S_OHs(n)%object) &
            , 'units', c=att_unit)
       CALL channel_halt(substr, status)
       !
       CALL d14co_check_unit(substr, 'OH', TRIM(att_unit) &
            , XSET(n)%lconv_OHs)
       !
       IF (p_parallel_io) &
            WRITE(*,*) '   mol/mol -> cm^(-3)      : ',XSET(n)%lconv_OHs


       SELECT CASE (XSET(n)%ohs_is_what)
       CASE (IS_GPOFFLINE)
!!#D attila +
#ifdef ECHAM5
          IF (S_SWITCH(n)%lg_lint) &
               ALLOCATE(XSET(n)%ohs_lg(NCELL))
#endif
!!#D attila -
       CASE (IS_GPTRACER)
!!#D attila +
#ifdef ECHAM5
          IF (S_SWITCH(n)%lg_lint) &
               ALLOCATE(XSET(n)%ohs_lg(NCELL))
#endif
!!#D attila -
          IF (S_SWITCH(n)%gp_lint) &
               ALLOCATE(XSET(n)%ohs_gp(nproma,nlev,ngpblks))
       CASE (IS_GPCHANNEL)
!!#D attila +
#ifdef ECHAM5
          IF (S_SWITCH(n)%lg_lint) &
               ALLOCATE(XSET(n)%ohs_lg(NCELL))
       CASE (IS_LGTRACER)
          IF (S_SWITCH(n)%gp_lint) &
               ALLOCATE(XSET(n)%ohs_gp(nproma,nlev,ngpblks))
          IF (S_SWITCH(n)%lg_lint) &
               ALLOCATE(XSET(n)%ohs_lg(NCELL))
       CASE (IS_LGCHANNEL)
          IF (S_SWITCH(n)%gp_lint) &
               ALLOCATE(XSET(n)%ohs_gp(nproma,nlev,ngpblks))
#endif
!!#D attila -
       CASE DEFAULT
          ! ERROR
       END SELECT
       ! ------------------------------------------------------------------

       ! TROPOSPHERIC OH -------------------------------------------------
       IF (p_parallel_io) THEN
          WRITE(*,*) 'INITIALIZING TROPOSPHERIC OH ...'
       END IF
       ! EVENT, GP-TRACER, LG-TRACER, OR CHANNEL OBJECT
       SELECT CASE(TRIM(ADJUSTL(S_OHt(n)%channel)))
       CASE (gp_channel)
          !
          ! GP-TRACER
          !
          XSET(n)%oht_is_what = IS_GPTRACER
          ! ... CHECK FOR SUBNAME
          CALL full2base_sub(status, TRIM(S_OHt(n)%object) &
               , basename, subname)
          CALL tracer_halt(substr, status)
          IF (p_parallel_io) THEN
             WRITE(*,*) ' TYPE                      : ONLINE (GP-TRACER)'
             WRITE(*,*) '   basename                : ',TRIM(basename)
             WRITE(*,*) '   subname                 : ',TRIM(subname)
          END IF
          CALL get_tracer(status, GPTRSTR, basename &
               , subname=subname                  &
               , idx=XSET(n)%idt_OHt)
          CALL tracer_halt(substr, status)
          !
!!#D attila +
#ifdef ECHAM5
       CASE (lg_channel)
          !
          ! LG-TRACER
          !
          XSET(n)%oht_is_what = IS_LGTRACER
          ! ... CHECK FOR SUBNAME
          CALL full2base_sub(status, TRIM(S_OHt(n)%object) &
               , basename, subname)
          CALL tracer_halt(substr, status)
          IF (p_parallel_io) THEN
             WRITE(*,*) ' TYPE                      : ONLINE (LG-TRACER)'
             WRITE(*,*) '   basename                : ',TRIM(basename)
             WRITE(*,*) '   subname                 : ',TRIM(subname)
          END IF
          CALL get_tracer(status, LGTRSTR, basename &
               , subname=subname                  &
               , idx=XSET(n)%idt_OHt)
          CALL tracer_halt(substr, status)
          !
          CALL get_channel_object(status  &
               , TRIM(S_OHt(n)%channel), TRIM(S_OHt(n)%object) &
               , p1=XSET(n)%oht_lg ) ! set LG-pointer 
#endif
!!#D attila -
          !
       CASE DEFAULT
          !
          ! CHANNEL OBJECT
          !
          IF (p_parallel_io) THEN
             WRITE(*,*) ' TYPE                      : ONLINE (CHANNEL OBJECT)'
             WRITE(*,*) '   channel                 : ',TRIM(S_OHt(n)%channel)
             WRITE(*,*) '   object                  : ',TRIM(S_OHt(n)%object)
          END IF
          CALL get_channel_object_info(status &
               , TRIM(S_OHt(n)%channel), TRIM(S_OHt(n)%object) &
               , reprid=reprid )
          CALL channel_halt(substr, status)
          IF (reprid == GP_3D_MID) THEN
             CALL get_channel_object(status  &
                  , TRIM(S_OHt(n)%channel), TRIM(S_OHt(n)%object) &
                  , p3=XSET(n)%oht_gp )
             XSET(n)%oht_is_what = IS_GPCHANNEL
!!#D attila +
#ifdef ECHAM5
          ELSEIF (reprid == LG_ATTILA) THEN
             CALL get_channel_object(status                       &
                  , TRIM(S_OHt(n)%channel), TRIM(S_OHt(n)%object) &
                  , p1=XSET(n)%oht_lg )
             XSET(n)%oht_is_what = IS_LGCHANNEL
#endif
!!#D attila -
          ELSE
             IF (TRIM(ADJUSTL(S_OHt(n)%channel)) == 'import_grid') THEN
                XSET(n)%OHt_is_what = IS_GPOFFLINE
                CALL get_channel_object(status                         &
                     , TRIM(S_OHt(n)%channel), TRIM(S_OHt(n)%object)   &
                     , p3=XSET(n)%import_OHt)
                CALL channel_halt(substr, status)
                ! get pressure levels from corresponding dimension variable
                CALL get_channel_object_dimvalue(status                &
                     ,TRIM(S_OHt(n)%channel), TRIM(S_OHt(n)%object)    &
                     , data=XSET(n)%paxis_OHt , axis ='Z')
                CALL channel_halt(substr, status)
                !
                CALL new_channel_object(status             &
                     , modstr, TRIM(S_OHt(n)%object)       &
                     , reprid=GP_3D_MID, p3=XSET(n)%OHt_gp)
                CALL channel_halt(substr, status)
             ELSE
                CALL error_bi(' representation not supported !', substr)
             END IF
          END IF
       END SELECT

       ! CHECK UNIT:
       ! potentially required conversion
       CALL get_attribute(status  &
            , TRIM(S_OHt(n)%channel), TRIM(S_OHt(n)%object) &
            , 'units', c=att_unit)
       CALL channel_halt(substr, status)
       !
       CALL d14co_check_unit(substr, 'OH', TRIM(att_unit) &
            , XSET(n)%lconv_OHt)
       !
       IF (p_parallel_io) &
            WRITE(*,*) '   mol/mol -> cm^(-3)      : ',XSET(n)%lconv_OHt

       SELECT CASE (XSET(n)%oht_is_what)
       CASE (IS_GPOFFLINE)
!!#D attila +
#ifdef ECHAM5
          IF (S_SWITCH(n)%lg_lint) &
               ALLOCATE(XSET(n)%oht_lg(NCELL))
#endif
!!#D attila -
       CASE (IS_GPTRACER)
!!#D attila +
#ifdef ECHAM5
          IF (S_SWITCH(n)%lg_lint) &
               ALLOCATE(XSET(n)%oht_lg(NCELL))
#endif
!!#D attila -
          IF (S_SWITCH(n)%gp_lint) &
               ALLOCATE(XSET(n)%oht_gp(nproma,nlev,ngpblks))
       CASE (IS_GPCHANNEL)
!!#D attila +
#ifdef ECHAM5
          IF (S_SWITCH(n)%lg_lint) &
               ALLOCATE(XSET(n)%oht_lg(NCELL))
       CASE (IS_LGTRACER)
          IF (S_SWITCH(n)%gp_lint) &
               ALLOCATE(XSET(n)%oht_gp(nproma,nlev,ngpblks))
          IF (S_SWITCH(n)%lg_lint) &
               ALLOCATE(XSET(n)%oht_lg(NCELL))
       CASE (IS_LGCHANNEL)
          IF (S_SWITCH(n)%gp_lint) &
               ALLOCATE(XSET(n)%oht_gp(nproma,nlev,ngpblks))
#endif
!!#D attila -
       CASE DEFAULT
          ! ERROR
       END SELECT
       ! ------------------------------------------------------------------

       ! TROPOPAUSE FOR OH ------------------------------------------------
       IF (p_parallel_io) THEN
          WRITE(*,*) 'INITIALIZING TROPOPAUSE FOR OH ...'
       END IF
       SELECT CASE(S_TP_OH(n)%typ)
       CASE(I_TP_DIAG)
          ! CHANNEL OBJECT
          IF (p_parallel_io) THEN
             WRITE(*,*) ' TYPE                      : ', &
                  'DIAGNOSED (CHANNEL OBJECT)'
             WRITE(*,*) '   channel                 : ', &
                  TRIM(S_TP_OH(n)%link%channel)
             WRITE(*,*) '   object                  : ', &
                  TRIM(S_TP_OH(n)%link%object)
          END IF
          CALL get_channel_object(status &
               , TRIM(S_TP_OH(n)%link%channel), TRIM(S_TP_OH(n)%link%object) &
               , p2=XSET(n)%xtp_oh )
          CALL channel_halt(substr, status)          
       CASE(I_TP_CLIM)
          IF (p_parallel_io) THEN
             WRITE(*,*) ' TYPE                      : ', &
                  'CLIMATOLOGICAL'
             WRITE(*,*) '   (',S_TP_OH(n)%clim(1)            &
                  ,' - ',S_TP_OH(n)%clim(2),' * cos^2(latitude) ) Pa'
          END IF
       CASE(I_TP_CONST)
          IF (p_parallel_io) THEN
             WRITE(*,*) ' TYPE                      : ', &
                  'CONSTANT PRESSURE'
             WRITE(*,*) '   (',S_TP_OH(n)%const,' Pa )'
          END IF
       CASE DEFAULT
          CALL error_bi(' unknown tropopause type !', substr)
       END SELECT
       ! ------------------------------------------------------------------

       ! TROPOPAUSE FOR STE -----------------------------------------------
       IF (p_parallel_io) THEN
          WRITE(*,*) 'INITIALIZING TROPOPAUSE FOR STE ...'
       END IF
       SELECT CASE(S_TP_STE(n)%typ)
       CASE(I_TP_DIAG)
          ! CHANNEL OBJECT
          IF (p_parallel_io) THEN
             WRITE(*,*) ' TYPE                      : ', &
                  'DIAGNOSED (CHANNEL OBJECT)'
             WRITE(*,*) '   channel                 : ', &
                  TRIM(S_TP_STE(n)%link%channel)
             WRITE(*,*) '   object                  : ', &
                  TRIM(S_TP_STE(n)%link%object)
          END IF
          CALL get_channel_object(status &
               , TRIM(S_TP_STE(n)%link%channel) &
               , TRIM(S_TP_STE(n)%link%object)  &
               , p2=XSET(n)%xtp_ste )
          CALL channel_halt(substr, status)
       CASE(I_TP_CLIM)
          IF (p_parallel_io) THEN
             WRITE(*,*) ' TYPE                      : ', &
                  'CLIMATOLOGICAL'
             WRITE(*,*) '   (',S_TP_STE(n)%clim(1)            &
                  ,' - ',S_TP_STE(n)%clim(2),' * cos^2(latitude) ) Pa'
          END IF
       CASE(I_TP_CONST)
          IF (p_parallel_io) THEN
             WRITE(*,*) ' TYPE                      : ', &
                  'CONSTANT PRESSURE'
             WRITE(*,*) '   (',S_TP_STE(n)%const,' Pa )'
          END IF
       CASE DEFAULT
          CALL error_bi(' unknown tropopause type !', substr)
       END SELECT
       ! ------------------------------------------------------------------

       IF (p_parallel_io) THEN
          WRITE(*,*) '-------------------------------------------------------'
       END IF

    END DO set_loop

!!#D attila +
#ifdef ECHAM5
    ! #######################################################################
    ! ### PART II: COUPLING TO LAGRANGIAN SCHEME
    ! #######################################################################
    IF (l_lg) THEN
       !
       IF (p_parallel_io) THEN
          WRITE(*,*) '================= COUPLING TO LG-SCHEME ==============='
       END IF
       !
       IF (p_parallel_io) THEN
          WRITE(*,*) ' CHANNEL         : ',TRIM(C_LG_CHANNEL)
       END IF
       !       
       IF (p_parallel_io) THEN
          WRITE(*,*) ' LONGITUDE INDEX : ',TRIM(C_LG_ILON)
       END IF
       CALL get_channel_object(status &
            , TRIM(C_LG_CHANNEL), TRIM(C_LG_ILON) &
            , p1=ilon_lg )
       CALL channel_halt(substr, status)
       !
       IF (p_parallel_io) THEN
          WRITE(*,*) ' LATITUDE INDEX  : ',TRIM(C_LG_ILAT)
       END IF
       CALL get_channel_object(status &
            , TRIM(C_LG_CHANNEL), TRIM(C_LG_ILAT) &
            , p1=ilat_lg )
       CALL channel_halt(substr, status)
       !
       IF (p_parallel_io) THEN
          WRITE(*,*) ' PRESSURE [Pa]   : ',TRIM(C_LG_PRESS)
       END IF
       CALL get_channel_object(status &
            , TRIM(C_LG_CHANNEL), TRIM(C_LG_PRESS) &
            , p1=p_lg )
       CALL channel_halt(substr, status)
       !
    END IF
#endif
!!#D attila -

    CALL end_message_bi(modstr,'COUPLING INITIALIZATION',substr)

  END SUBROUTINE d14co_init_coupling
! ----------------------------------------------------------------------

! ------------------------------------------------------------------------
  SUBROUTINE d14co_global_start

    ! 14CO MODULE ROUTINE (ECHAM-5 INTERFACE, PUBLIC)
    !
    ! INPUT OF OFFLINE FIELDS (OHs, OHt, 14CO source)
    ! TRIGGERED BY SPECIAL EVENTs
    ! I/O in PARALLEL ENVIRONMENT, DATA TRANSFER TO CHANNEL OBJECTS
    !
    ! Author: Patrick Joeckel, MPICH, Feb 2002

    USE messy_main_blather_bi,    ONLY: error_bi
!!#D attila +
#ifdef ECHAM5
    USE messy_attila_tools_e5,    ONLY: lg2gp_e5, gp2lg_e5, LG2GP_AVE
    USE messy_main_tracer_mem_bi, ONLY: xtm1_a, xtte_a, NCELL
    USE messy_main_grid_def_mem_bi,ONLY: nproma, nlev, ngpblks
#endif
!!#D attila -
    USE messy_main_tracer_mem_bi, ONLY: xtm1, xtte
    USE messy_main_timer,         ONLY: ztmst=>time_step_len 

    IMPLICIT NONE

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER     :: substr = 'd14co_global_start'
    INTEGER                         :: n    ! loop counter
!!#D attila +
#ifdef ECHAM5
    REAL(DP), DIMENSION(:),     POINTER :: tmp_lg => NULL()
    REAL(DP), DIMENSION(:,:,:), POINTER :: tmp_gp => NULL()
#endif
!!#D attila -

    ! ### SET GP FIELDS #################################################
    set_loop2: DO n=1, NMAXSETS

       IF (.NOT. (S_SWITCH(n)%gp_lint .OR. S_SWITCH(n)%lg_lint)) CYCLE

       SELECT CASE(XSET(n)%q_is_what)
       CASE (IS_GPOFFLINE)
          ! NOTHING TO DO FOR GP: ALREADY SET AFTER IMPORT 
!!#D attila +
          ! LG: UPDATE POSSIBLE ONLY AFTER VERTICAL INTERPOLATION OF GP-FIELD
!!#D attila -
       CASE (IS_GPTRACER)
          ! NOT IMPLEMENTED
          CALL error_bi(' 14CO SOURCE MUST NOT BE GP-TRACER', substr)
!!#D attila +
#ifdef ECHAM5
       CASE (IS_LGTRACER)
          ! NOT IMPLEMENTED
          CALL error_bi(' 14CO SOURCE MUST NOT BE LG-TRACER', substr)
#endif
!!#D attila -
       CASE (IS_GPCHANNEL)
          ! NOTHING TO DO FOR GP: ALREADY SET IN INIT_COUPLING
!!#D attila +
#ifdef ECHAM5
          ! CONVERT FOR LG:
          IF (S_SWITCH(n)%lg_lint) THEN
             CALL gp2lg_e5(XSET(n)%q_gp, XSET(n)%q_lg, lmcons=.false.)
          END IF
       CASE (IS_LGCHANNEL)
          ! CONVERT FOR GP:
          CALL lg2gp_e5(XSET(n)%q_lg, XSET(n)%q_gp, LG2GP_AVE  &
               , .TRUE., .FALSE.)
          ! NOTHING TO DO FOR LG: ALREADY SET IN INIT_COUPLING
#endif
!!#D attila -
       CASE DEFAULT
          ! ERROR
          CALL error_bi(' 14CO SOURCE IS UNDEFINED', substr)
       END SELECT

       SELECT CASE(XSET(n)%ohs_is_what)
       CASE (IS_GPOFFLINE)
          ! NOTHING TO DO FOR GP: ALREADY SET AFTER IMPORT
!!#D attila +
          ! LG: UPDATE POSSIBLE ONLY AFTER VERTICAL INTERPOLATION OF GP-FIELD
!!#D attila -
       CASE (IS_GPTRACER)
          ! GP (SPACE ALLOCATED IN INIT_COUPLING)
          IF (S_SWITCH(n)%gp_lint) THEN
             XSET(n)%ohs_gp(:,:,:) =                    &
                  xtm1(:,:,XSET(n)%idt_OHs,:)           &
                  + xtte(:,:,XSET(n)%idt_OHs,:) * ztmst 
          END IF
!!#D attila +
#ifdef ECHAM5
          ! LG:
          IF (S_SWITCH(n)%lg_lint) THEN
             ALLOCATE(tmp_gp(nproma,nlev,ngpblks))
             tmp_gp(:,:,:) =                            &
                  xtm1(:,:,XSET(n)%idt_OHs,:)           &
                  + xtte(:,:,XSET(n)%idt_OHs,:) * ztmst
             CALL gp2lg_e5(tmp_gp, XSET(n)%ohs_lg, lmcons=.false.)
             DEALLOCATE(tmp_gp)
          END IF
       CASE (IS_LGTRACER)
          ! GP:
          IF (S_SWITCH(n)%gp_lint) THEN         
             ALLOCATE(tmp_lg(NCELL))
             tmp_lg(:) = xtm1_a(:,1,XSET(n)%idt_OHs,1)    &
                  + xtte_a(:,1,XSET(n)%idt_OHs,1) * ztmst
             CALL lg2gp_e5(tmp_lg, XSET(n)%ohs_gp, LG2GP_AVE  &
                  , .TRUE., .FALSE.)
             DEALLOCATE(tmp_lg)
          END IF
          ! LG (SPACE ALLOCATED IN INIT_COUPLING):
          IF (S_SWITCH(n)%lg_lint) THEN  
             XSET(n)%ohs_lg(:) =                       &
                  xtm1_a(:,1,XSET(n)%idt_OHs,1)        &
                  + xtte(:,1,XSET(n)%idt_OHs,1) * ztmst             
          END IF
#endif
!!#D attila -
       CASE (IS_GPCHANNEL)
          ! NOTHING TO DO FOR GP: ALREADY SET IN INIT_COUPLING
!!#D attila +
#ifdef ECHAM5
          ! CONVERT FOR LG:
          IF (S_SWITCH(n)%lg_lint) THEN
             CALL gp2lg_e5(XSET(n)%ohs_gp, XSET(n)%ohs_lg, lmcons=.false.)
          END IF
       CASE (IS_LGCHANNEL)
          ! CONVERT FOR GP:
          CALL lg2gp_e5(XSET(n)%ohs_lg, XSET(n)%ohs_gp, LG2GP_AVE  &
               , .TRUE., .FALSE.)
          ! NOTHING TO DO FOR LG: ALREADY SET IN INIT_COUPLING
#endif
!!#D attila -
       CASE DEFAULT
          ! ERROR
          CALL error_bi(' OHs IS UNDEFINED', substr)
       END SELECT

       SELECT CASE(XSET(n)%oht_is_what)
       CASE (IS_GPOFFLINE)
          ! NOTHING TO DO FOR GP: ALREADY SET AFTER EVENT UPDATE (ABOVE)
!!#D attila +
          ! LG: UPDATE POSSIBLE ONLY AFTER VERTICAL INTERPOLATION OF GP-FIELD
!!#D attila -
       CASE (IS_GPTRACER)
          ! GP (SPACE ALLOCATED IN INIT_COUPLING)
          IF (S_SWITCH(n)%gp_lint) THEN
             XSET(n)%oht_gp(:,:,:) =                    &
                  xtm1(:,:,XSET(n)%idt_OHt,:)           &
                  + xtte(:,:,XSET(n)%idt_OHt,:) * ztmst 
          END IF
!!#D attila +
#ifdef ECHAM5
          ! LG:
          IF (S_SWITCH(n)%lg_lint) THEN
             ALLOCATE(tmp_gp(nproma,nlev,ngpblks))
             tmp_gp(:,:,:) =                            &
                  xtm1(:,:,XSET(n)%idt_OHt,:)           &
                  + xtte(:,:,XSET(n)%idt_OHt,:) * ztmst
             CALL gp2lg_e5(tmp_gp, XSET(n)%oht_lg, lmcons=.false.)
             DEALLOCATE(tmp_gp)
          END IF
       CASE (IS_LGTRACER)
          ! GP:
          IF (S_SWITCH(n)%gp_lint) THEN         
             ALLOCATE(tmp_lg(NCELL))
             tmp_lg(:) = xtm1_a(:,1,XSET(n)%idt_OHt,1)    &
                  + xtte_a(:,1,XSET(n)%idt_OHt,1) * ztmst
             CALL lg2gp_e5(tmp_lg, XSET(n)%oht_gp, LG2GP_AVE  &
                  , .TRUE., .FALSE.)
             DEALLOCATE(tmp_lg)
          END IF
          ! LG (SPACE ALLOCATED IN INIT_COUPLING):
          IF (S_SWITCH(n)%lg_lint) THEN  
             XSET(n)%oht_lg(:) =                       &
                  xtm1_a(:,1,XSET(n)%idt_OHt,1)        &
                  + xtte(:,1,XSET(n)%idt_OHt,1) * ztmst             
          END IF
#endif
!!#D attila -
       CASE (IS_GPCHANNEL)
          ! NOTHING TO DO FOR GP: ALREADY SET IN INIT_COUPLING
!!#D attila +
#ifdef ECHAM5
          ! CONVERT FOR LG:
          IF (S_SWITCH(n)%lg_lint) THEN
             CALL gp2lg_e5(XSET(n)%oht_gp, XSET(n)%oht_lg, lmcons=.false.)
          END IF
       CASE (IS_LGCHANNEL)
          ! CONVERT FOR GP:
          CALL lg2gp_e5(XSET(n)%oht_lg, XSET(n)%oht_gp, LG2GP_AVE  &
               , .TRUE., .FALSE.)
          ! NOTHING TO DO FOR LG: ALREADY SET IN INIT_COUPLING
#endif
!!#D attila -
       CASE DEFAULT
          ! ERROR
          CALL error_bi(' OHt IS UNDEFINED', substr)
       END SELECT

    END DO set_loop2

  END SUBROUTINE d14co_global_start
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
  SUBROUTINE d14co_physc

    USE messy_main_grid_def_mem_bi, ONLY: jrow, kproma  
    USE messy_main_grid_def_bi,     ONLY: philat_2d
    USE messy_main_data_bi,         ONLY: pmid => press_3d
    USE messy_main_timer,         ONLY: ztmst=>time_step_len
    USE messy_main_constants_mem, ONLY: pi, k_b, N_A, R_gas, M_air, M_H2O

    IMPLICIT NONE

    ! LOCAL
    !CHARACTER(LEN=*), PARAMETER :: substr = 'd14co_physc'
    REAL(dp),         PARAMETER :: vtmpc1 = M_air / M_H2O - 1._dp
    INTEGER     :: jp                           ! column counter (1...kproma)
    INTEGER     :: n                            ! set loop counter
    REAL(DP), DIMENSION(:)  , ALLOCATABLE :: zptpoh  ! OH-tropopause pressure
    REAL(DP), DIMENSION(:)  , ALLOCATABLE :: zptp    ! STE-tropopause pressure

    ! PERFORM VERTICAL INTERPOLATION
    ! + THIS NEEDS TO BE DONE HERE WITHIN THE LOCAL LOOP AND AT EACH
    !   TIME STEP, SINCE IT INVOLVES THE TIME DEPENDENT PRESSURE
    !   FIELD
    interpol_loop: DO n=1, NMAXSETS

       IF (XSET(n)%q_is_what == IS_GPOFFLINE) THEN          
          vector_loop_q: DO jp=1, kproma
             ! RESET
             XSET(n)%q_gp(_RI_XYZ__(jp,jrow,:)) = 0.0_DP
             !
             CALL intpol_p( XSET(n)%import_q(_RI_XYZ__(jp,jrow,:))  &
                  , XSET(n)%paxis_q(:)                              &
                  , pmid(_RI_XYZ__(jp,jrow,:))                      &
                  , XSET(n)%q_gp(_RI_XYZ__(jp,jrow,:))              &
                  , .TRUE. )
          END DO vector_loop_q
       END IF
          
       IF (XSET(n)%OHs_is_what == IS_GPOFFLINE) THEN
          vector_loop_ohs: DO jp=1, kproma
             ! RESET
             XSET(n)%OHs_gp(_RI_XYZ__(jp,jrow,:)) = 0.0_DP
             !
             CALL intpol_p( XSET(n)%import_OHs(_RI_XYZ__(jp,jrow,:))  &
                  , XSET(n)%paxis_OHs(:)                              &
                  , pmid(_RI_XYZ__(jp,jrow,:))                        &
                  , XSET(n)%OHs_gp(_RI_XYZ__(jp,jrow,:))              &
                  , .FALSE. )
          END DO vector_loop_ohs
       END IF
          
       IF (XSET(n)%OHt_is_what == IS_GPOFFLINE) THEN
          vector_loop_oht: DO jp=1, kproma
             ! RESET
             XSET(n)%OHt_gp(_RI_XYZ__(jp,jrow,:)) = 0.0_DP
             !
             CALL intpol_p( XSET(n)%import_OHt(_RI_XYZ__(jp,jrow,:))  &
                  , XSET(n)%paxis_OHt(:)                              &
                  , pmid(_RI_XYZ__(jp,jrow,:))                        &
                  , XSET(n)%OHt_gp(_RI_XYZ__(jp,jrow,:))              &
                  , .FALSE. )          
          END DO vector_loop_oht
       END IF
       
    END DO interpol_loop

    ! OPERATIONS REQUIRED FOR GP AND LG
    ALLOCATE(zptp(kproma))
    ALLOCATE(zptpoh(kproma))

    set_loop: DO n=1, NMAXSETS

       IF (.NOT. (S_SWITCH(n)%gp_lint .OR. S_SWITCH(n)%lg_lint)) CYCLE
       
       ! RESET
       zptp(:)   = 0.0_DP
       zptpoh(:) = 0.0_DP       

       ! STE TROPOPAUSE
       SELECT CASE(S_TP_STE(n)%typ)
       CASE(I_TP_DIAG)
          zptp(:) = XSET(n)%xtp_ste(1:kproma,jrow)
       CASE(I_TP_CLIM)
          DO jp = 1, kproma
             zptp(jp) = S_TP_STE(n)%clim(1) &
               - S_TP_STE(n)%clim(2) * (cos((philat_2d(jp,jrow)/180.)*pi))**2
          END DO
       CASE(I_TP_CONST)
          zptp(:) = S_TP_STE(n)%const
       END SELECT
       ! COPY TO CHANNEL OBJECT
       XSET(n)%ptp(1:kproma,jrow) = zptp(:)

       ! OH TROPOPAUSE
       SELECT CASE(S_TP_OH(n)%typ)
       CASE(I_TP_DIAG)
          zptpoh(:) = XSET(n)%xtp_oh(1:kproma,jrow)
       CASE(I_TP_CLIM)
          DO jp = 1, kproma
             zptpoh(jp) =  S_TP_OH(n)%clim(1) &
                  - S_TP_OH(n)%clim(2) * (cos((philat_2d(jp,jrow)/180.)*pi))**2
          END DO
       CASE(I_TP_CONST)
          zptpoh(:) = S_TP_OH(n)%const
       END SELECT
       XSET(n)%ptpoh(1:kproma,jrow) = zptpoh(:)

    END DO set_loop

    ! CLEAN UP
    DEALLOCATE(zptp)
    DEALLOCATE(zptpoh)

    CALL d14co_physc_gp

  CONTAINS

    ! --------------------------------------------------------------------
    SUBROUTINE d14co_physc_gp

      USE messy_main_tracer_mem_bi,   ONLY: xtm1, pxtte=>qxtte
      USE messy_main_grid_def_mem_bi, ONLY: nlev
      USE messy_main_data_bi,       ONLY: pmid => press_3d  &
                                        , pint => pressi_3d &
                                        , tm1 => tm1_3d     &
                                        , tte => tte_3d     &
                                        , qm1 => qm1_3d     &
                                        , qte => qte_3d

      IMPLICIT NONE

      ! LOCAL
      REAL(DP), DIMENSION(:,:,:), POINTER   :: pxtm1 => NULL() ! tracer at t-1
      ! tracer IDs
      INTEGER                               :: idt_14co, idt_14cos, idt_14cot
      REAL(DP), DIMENSION(:,:), POINTER     :: temp_gp => NULL() ! temp. (K)
      REAL(DP), DIMENSION(:,:), ALLOCATABLE :: zohs    ! stratos. OH field
      REAL(DP), DIMENSION(:,:), ALLOCATABLE :: zoht    ! tropos. OH field
      REAL(DP), DIMENSION(:),   ALLOCATABLE :: ztf     ! tropos. frac. of box
      REAL(DP), DIMENSION(:,:), ALLOCATABLE :: z14co   ! 14CO
      REAL(DP), DIMENSION(:,:), ALLOCATABLE :: z14cot  ! tropos. 14CO
      REAL(DP), DIMENSION(:,:), ALLOCATABLE :: z14cos  ! stratos. 14CO
      REAL(DP), DIMENSION(:,:), ALLOCATABLE :: zt14co  ! 14CO tend.
      REAL(DP), DIMENSION(:,:), ALLOCATABLE :: zt14cot ! troposph. 14CO tend.
      REAL(DP), DIMENSION(:,:), ALLOCATABLE :: zt14cos ! stratosph. 14CO tend.
      INTEGER                               :: n       ! set counter
      INTEGER                               :: jp      ! column counter
      REAL(DP), DIMENSION(:,:), POINTER     :: sphum_gp ! specific humidity
      !                                                != m(H2O)/m(air) [kg/kg]
      REAL(DP), DIMENSION(:,:), POINTER     :: cair_gp  ! conc. of air [cm^-3]


      ! INIT
      pxtm1 => xtm1(:,:,:,jrow)

      ! UPDATE ... 
      ! ... TEMPERATURE
      ALLOCATE(temp_gp(kproma,SIZE(tm1,2)))
      temp_gp(1:kproma,:) = tm1(_RI_XYZ__(1:kproma,jrow,:)) &
           + tte(_RI_XYZ__(1:kproma,jrow,:))*ztmst

      ! ... SPECIFIC HUMIDITY
      ALLOCATE(sphum_gp(kproma,SIZE(qm1,2)))
      sphum_gp(1:kproma,:) = MAX( &
           qm1(_RI_XYZ__(1:kproma,jrow,:))  &
           + qte(_RI_XYZ__(1:kproma,jrow,:)) * ztmst,  &
           0._dp)

      ! ... CONCENTRATION OF AIR
      ALLOCATE(cair_gp(kproma,SIZE(sphum_gp,2)))
      cair_gp(1:kproma,:)  = (N_A/1.E6_dp) * pmid(_RI_XYZ__(1:kproma,jrow,:))/ &
           (R_gas*temp_gp(1:kproma,:)*(1.0_dp+vtmpc1*sphum_gp(1:kproma,:)))

      ! ALLOCATE MEMORY
      ALLOCATE(zohs(kproma,nlev))
      ALLOCATE(zoht(kproma,nlev))
      ALLOCATE(ztf(nlev))
      ALLOCATE(z14co(kproma,nlev))
      ALLOCATE(z14cot(kproma,nlev))
      ALLOCATE(z14cos(kproma,nlev))
      ALLOCATE(zt14co(kproma,nlev))
      ALLOCATE(zt14cot(kproma,nlev))
      ALLOCATE(zt14cos(kproma,nlev))

      set_loop: DO n=1, NMAXSETS

         IF (.NOT. S_SWITCH(n)%gp_lint) CYCLE

         ! RESET
         zohs(:,:)    = 0.0_DP
         zoht(:,:)    = 0.0_DP
         ztf(:)       = 0.0_DP
         z14co(:,:)   = 0.0_DP
         z14cot(:,:)  = 0.0_DP
         z14cos(:,:)  = 0.0_DP
         zt14co(:,:)  = 0.0_DP
         zt14cot(:,:) = 0.0_DP
         zt14cos(:,:) = 0.0_DP

         ! SOURCE
         XSET(n)%p14co_gp(_RI_XYZ__(1:kproma,jrow,:)) = &
              XSET(n)%q_gp(_RI_XYZ__(1:kproma,jrow,:))
         ! UNIT CONVERSION
         IF (XSET(n)%lconv_q) THEN
            DO jp = 1, kproma
               CALL mgs2vmrs(XSET(n)%p14co_gp(_RI_XYZ__(jp,jrow,:)))
            END DO
         END IF
         
         ! OHs
         zohs(:,:) = XSET(n)%ohs_gp(_RI_XYZ__(1:kproma,jrow,:))
         ! UNIT CONVERSION
         IF (XSET(n)%lconv_OHs) THEN
            DO jp = 1, kproma
               CALL vmr2conc(zohs(jp,:),pmid(_RI_XYZ__(jp,jrow,:)) &
                    , temp_gp(jp,:),k_b)
            END DO
         END IF

         ! OHt
         zoht(:,:) = XSET(n)%oht_gp(_RI_XYZ__(1:kproma,jrow,:))
         ! UNIT CONVERSION
         IF (XSET(n)%lconv_OHt) THEN
            DO jp = 1, kproma
               CALL vmr2conc(zoht(jp,:),pmid(_RI_XYZ__(jp,jrow,:)) &
                    , temp_gp(jp,:),k_b)
            END DO
         END IF

         ! MERGE OHs AND OHt (put into channel object)
         DO jp = 1, kproma
            CALL merge_p(zohs(jp,:),zoht(jp,:),pmid(_RI_XYZ__(jp,jrow,:)), &
                 XSET(n)%ptpoh(jp,jrow),XSET(n)%oh_gp(_RI_XYZ__(jp,jrow,:)))
         END DO
         
         ! SET TRACER IDs
         idt_14co  = XSET(n)%idt_gp_14CO
         idt_14cos = XSET(n)%idt_gp_14COs
         idt_14cot = XSET(n)%idt_gp_14COt
         
         ! ADJUST TRACER TENDENCIES TO FORCE 14CO = 14COs + 14COt
         IF (S_SWITCH(n)%gp_ladjt) THEN
            DO jp = 1, kproma
               CALL adj_tend(pxtm1(_RI_X_ZN_(jp,:,idt_14co))     &
                    , pxtte(_RI_X_ZN_(jp,:,idt_14co))  &
                    , pxtm1(_RI_X_ZN_(jp,:,idt_14cos)) &
                    , pxtte(_RI_X_ZN_(jp,:,idt_14cos)) &
                    , pxtm1(_RI_X_ZN_(jp,:,idt_14cot)) &
                    , pxtte(_RI_X_ZN_(jp,:,idt_14cot)) &
                    , ztmst )
            END DO
         END IF
         
         ! INTEGRATE CHEM
         ! 1. START VALUES FOR TIME INTEGRATION
         z14co(:,:)  = pxtm1(_RI_X_ZN_(1:kproma,:,idt_14co))  &
              + pxtte(_RI_X_ZN_(1:kproma,:,idt_14co))*ztmst
         z14cos(:,:) = pxtm1(_RI_X_ZN_(1:kproma,:,idt_14cos)) &
              + pxtte(_RI_X_ZN_(1:kproma,:,idt_14cos))*ztmst
         z14cot(:,:) = pxtm1(_RI_X_ZN_(1:kproma,:,idt_14cot)) &
              + pxtte(_RI_X_ZN_(1:kproma,:,idt_14cot))*ztmst

         DO jp = 1, kproma
            ! 2. FRACTION OF GRID BOX IN TROPOSPHERE
            ztf(:) = &
                 tfrac(pint(jp,1:nlev,jrow), pmid(jp,:,jrow) &
                 , pint(jp,2:nlev+1,jrow), XSET(n)%ptp(jp,jrow))

            ! 3. SOLVE DIFF. EQ.
            CALL int_14CO(ztmst, cair_gp(jp,:)                     &
                 , XSET(n)%p14co_gp(_RI_XYZ__(jp,jrow,:))          &
                 , XSET(n)%oh_gp(_RI_XYZ__(jp,jrow,:))             &
                 , z14co(jp,:),  z14cos(jp,:),  z14cot(jp,:)       &
                 , ztf(:)                                          &
                 , zt14co(jp,:), zt14cos(jp,:), zt14cot(jp,:)      &
                 , XSET(n)%l14co_gp(_RI_XYZ__(jp,jrow,:))          &
                 )

            ! 4. UPDATE TENDENCY
            pxtte(_RI_X_ZN_(jp,:,idt_14co))  = &
                 pxtte(_RI_X_ZN_(jp,:,idt_14co))  + zt14co(jp,:)
            pxtte(_RI_X_ZN_(jp,:,idt_14cos)) = &
                 pxtte(_RI_X_ZN_(jp,:,idt_14cos)) + zt14cos(jp,:)
            pxtte(_RI_X_ZN_(jp,:,idt_14cot)) = &
                 pxtte(_RI_X_ZN_(jp,:,idt_14cot)) + zt14cot(jp,:)
         END DO
         
      END DO set_loop
      
      ! DEALLOCATE MEMORY
      DEALLOCATE(temp_gp)
      DEALLOCATE(sphum_gp)
      DEALLOCATE(cair_gp)
      !
      DEALLOCATE(zohs)
      DEALLOCATE(zoht)
      DEALLOCATE(ztf)
      DEALLOCATE(z14co)
      DEALLOCATE(z14cot)
      DEALLOCATE(z14cos)
      DEALLOCATE(zt14co)
      DEALLOCATE(zt14cot)
      DEALLOCATE(zt14cos)

    END SUBROUTINE d14co_physc_gp
    ! --------------------------------------------------------------------

  END SUBROUTINE d14co_physc
! --------------------------------------------------------------------

! --------------------------------------------------------------------
  SUBROUTINE d14co_global_end

    IMPLICIT NONE
!!#D attila +
#ifdef ECHAM5
    CALL d14co_global_end_lg

  CONTAINS

    ! --------------------------------------------------------------------
    SUBROUTINE d14co_global_end_lg

      USE messy_main_tracer_mem_bi, ONLY: xtm1_a, pxtte=>qxtte_a, NCELL
      USE messy_main_data_bi,       ONLY: tm1 => tm1_3d &
                                        , tte => tte_3d &
                                        , qm1 => qm1_3d &
                                        , qte => qte_3d
      USE messy_main_grid_def_mem_bi, ONLY: nproma, nlev, ngpblks 
      USE messy_attila_tools_e5,    ONLY: gp2lg_e5
      USE messy_main_transform_bi,  ONLY: trp_gpdc_gpgl
      USE messy_main_timer,         ONLY: ztmst=>time_step_len
      USE messy_main_constants_mem, ONLY: N_A, R_gas, M_air, M_H2O, k_b

      IMPLICIT NONE

      INTRINSIC :: NINT

      ! LOCAL
      REAL(dp),         PARAMETER :: vtmpc1 = M_air / M_H2O - 1._dp
      REAL(DP), DIMENSION(:,:),   POINTER :: pxtm1 => NULL() ! tracer at t-1
      ! tracer IDs
      INTEGER                             :: idt_14co, idt_14cos, idt_14cot
      ! temperature GP
      REAL(DP), DIMENSION(:,:,:), POINTER :: temp_gp  => NULL()
      ! temperature LG
      REAL(DP), DIMENSION(:),     POINTER :: temp_lg  => NULL()
      ! globalized tropoause (STE)
      REAL(DP), DIMENSION(:,:),   POINTER :: zptp_g   => NULL()
      ! globalized tropoause (OH)
      REAL(DP), DIMENSION(:,:),   POINTER :: zptpoh_g => NULL()
      REAL(DP), DIMENSION(:), ALLOCATABLE :: zohs    ! stratospheric OH field
      REAL(DP), DIMENSION(:), ALLOCATABLE :: zoht    ! tropospheric OH field
      REAL(DP), DIMENSION(:), ALLOCATABLE :: ztf     ! tropos. frac. of box 
      REAL(DP), DIMENSION(:), ALLOCATABLE :: z14co   ! 14CO
      REAL(DP), DIMENSION(:), ALLOCATABLE :: z14cot  ! tropos. 14CO
      REAL(DP), DIMENSION(:), ALLOCATABLE :: z14cos  ! strat. 14CO
      REAL(DP), DIMENSION(:), ALLOCATABLE :: zt14co  ! 14CO tend.
      REAL(DP), DIMENSION(:), ALLOCATABLE :: zt14cot ! troposph. 14CO tend.
      REAL(DP), DIMENSION(:), ALLOCATABLE :: zt14cos ! stratosph. 14CO tend.
      INTEGER                             :: n       ! set counter
      INTEGER                             :: jc      ! cell counter
      INTEGER                             :: jxg, jyg ! indices in global field
      ! specific humidity GP
      REAL(DP), DIMENSION(:,:,:), POINTER :: sphum_gp  => NULL()
      ! specific humidity LG
      REAL(DP), DIMENSION(:),     POINTER :: sphum_lg  => NULL()
      ! concentration of air LG [cm^-3]
      REAL(DP), DIMENSION(:),     POINTER :: cair_lg  => NULL()

      IF (.NOT. LCALC_LG) RETURN

      ! INIT
      pxtm1 => xtm1_a(:,1,:,1)

      ! UPDATE
      ! ... TEMPERATURE
      ALLOCATE(temp_gp(nproma,nlev,ngpblks))
      ALLOCATE(temp_lg(NCELL))
      temp_gp(:,:,:) = tm1(:,:,:) + tte(:,:,:) * ztmst
      CALL gp2lg_e5(temp_gp, temp_lg, lmcons=.false.)

      ! ... SPECIFIC HUMIDITY
      ALLOCATE(sphum_gp(nproma,nlev,ngpblks))
      ALLOCATE(sphum_lg(NCELL))
      sphum_gp(:,:,:) = MAX(qm1(:,:,:) + qte(:,:,:) * ztmst, 0._dp)
      CALL gp2lg_e5(sphum_gp, sphum_lg, lmcons=.false.)

      ! ... CONCENTRATION OF AIR
      ALLOCATE(cair_lg(NCELL))
      cair_lg(:)  = (N_A/1.E6_dp) * p_lg(:) / &
           (R_gas*temp_lg(:)*(1.0_dp+vtmpc1*sphum_lg(:)))

      ! ALLOCATE MEMORY
      ALLOCATE(zohs(NCELL))
      ALLOCATE(zoht(NCELL))
      ALLOCATE(ztf(NCELL))
      ALLOCATE(z14co(NCELL))
      ALLOCATE(z14cot(NCELL))
      ALLOCATE(z14cos(NCELL))
      ALLOCATE(zt14co(NCELL))
      ALLOCATE(zt14cot(NCELL))
      ALLOCATE(zt14cos(NCELL))

      set_loop: DO n=1, NMAXSETS

         IF (.NOT. S_SWITCH(n)%lg_lint) CYCLE

         ! CONVERT q, ohs, oht IN CASE OF OFFLINE FIELDS
         ! (THIS NEEDS TO BE DONE AFTER THE VERTICAL INTERPOLATION IN GP !)
         IF (XSET(n)%q_is_what == IS_GPOFFLINE) THEN
            CALL gp2lg_e5(XSET(n)%q_gp, XSET(n)%q_lg, lmcons=.false.)
         END IF
         !
         IF (XSET(n)%ohs_is_what == IS_GPOFFLINE) THEN
            CALL gp2lg_e5(XSET(n)%ohs_gp, XSET(n)%ohs_lg, lmcons=.false.)
         END IF
         !
         IF (XSET(n)%oht_is_what == IS_GPOFFLINE) THEN
            CALL gp2lg_e5(XSET(n)%oht_gp, XSET(n)%oht_lg, lmcons=.false.)
         END IF

         ! GLOBALIZE TROPOPAUSE(s)
         CALL trp_gpdc_gpgl(1, XSET(n)%ptp,   zptp_g)
         CALL trp_gpdc_gpgl(1, XSET(n)%ptpoh, zptpoh_g)

         ! SOURCE
         XSET(n)%p14co_lg(:) = XSET(n)%q_lg(:)
         ! UNIT CONVERSION
         IF (XSET(n)%lconv_q) THEN
            CALL mgs2vmrs(XSET(n)%p14co_lg(:))
         END IF
         
         ! OHs
         zohs(:) = XSET(n)%ohs_lg(:)
         ! UNIT CONVERSION
         IF (XSET(n)%lconv_OHs) THEN
            CALL vmr2conc(zohs(:), p_lg(:), temp_lg(:),k_b)
         END IF

         ! OHt
         zoht(:) = XSET(n)%oht_lg(:)
         ! UNIT CONVERSION
         IF (XSET(n)%lconv_OHt) THEN
            CALL vmr2conc(zoht(:), p_lg(:), temp_lg(:),k_b)
         END IF

         ! MERGE OHs AND OHt (put into channel object)
         DO jc = 1, NCELL
            jxg = NINT(ilon_lg(jc))
            jyg = NINT(ilat_lg(jc))
            CALL merge_p(zohs(jc),zoht(jc), p_lg(jc), &
                 zptpoh_g(jxg,jyg), XSET(n)%oh_lg(jc))
         END DO

         ! SET TRACER IDs
         idt_14co  = XSET(n)%idt_lg_14CO
         idt_14cos = XSET(n)%idt_lg_14COs
         idt_14cot = XSET(n)%idt_lg_14COt
         
         ! ADJUST TRACER TENDENCIES TO FORCE 14CO = 14COs + 14COt
         IF (S_SWITCH(n)%lg_ladjt) THEN
            CALL adj_tend(pxtm1(:,idt_14co), pxtte(:,idt_14co) &
                 ,pxtm1(:,idt_14cos),pxtte(:,idt_14cos)        &
                 ,pxtm1(:,idt_14cot),pxtte(:,idt_14cot)        &
                 ,ztmst )
         END IF
         
         ! INTEGRATE CHEM
         ! 1. START VALUES FOR TIME INTEGRATION
         z14co(:)  = pxtm1(:,idt_14co)  + pxtte(:,idt_14co)*ztmst
         z14cos(:) = pxtm1(:,idt_14cos) + pxtte(:,idt_14cos)*ztmst
         z14cot(:) = pxtm1(:,idt_14cot) + pxtte(:,idt_14cot)*ztmst

         ! 2. FRACTION OF GRID BOX IN TROPOSPHERE
         DO jc = 1, NCELL
            jxg = NINT(ilon_lg(jc))
            jyg = NINT(ilat_lg(jc))
            ztf(jc) = tfrac(0.0_DP, p_lg(jc), 0.0_DP, zptp_g(jxg,jyg))
         END DO

         ! 3. SOLVE DIFF. EQ.
         CALL int_14CO(ztmst, cair_lg(:)               &
              , XSET(n)%p14co_lg(:), XSET(n)%oh_lg(:)  &
              , z14co(:),  z14cos(:),  z14cot(:)       &
              , ztf(:)                                 &
              , zt14co(:), zt14cos(:), zt14cot(:)      &
              , XSET(n)%l14co_lg(:)                    &
              )
         
         ! 4. UPDATE TENDENCY
         pxtte(:,idt_14co)  = pxtte(:,idt_14co)  + zt14co(:)
         pxtte(:,idt_14cos) = pxtte(:,idt_14cos) + zt14cos(:)
         pxtte(:,idt_14cot) = pxtte(:,idt_14cot) + zt14cot(:)
         
         ! CLEAN UP
         IF (ASSOCIATED(zptp_g)) THEN
            DEALLOCATE(zptp_g)
            NULLIFY(zptp_g)
         END IF
         IF (ASSOCIATED(zptpoh_g)) THEN
            DEALLOCATE(zptpoh_g)
            NULLIFY(zptpoh_g)
         END IF

      END DO set_loop
      
      ! CLEAN UP
      DEALLOCATE(temp_gp)
      DEALLOCATE(temp_lg)
      DEALLOCATE(sphum_gp)
      DEALLOCATE(sphum_lg)
      DEALLOCATE(cair_lg)
      !
      DEALLOCATE(zohs)
      DEALLOCATE(zoht)
      DEALLOCATE(ztf)
      DEALLOCATE(z14co)
      DEALLOCATE(z14cot)
      DEALLOCATE(z14cos)
      DEALLOCATE(zt14co)
      DEALLOCATE(zt14cot)
      DEALLOCATE(zt14cos)

    END SUBROUTINE d14co_global_end_lg
    ! --------------------------------------------------------------------
#endif
!!#D attila -

  END SUBROUTINE d14co_global_end
! ----------------------------------------------------------------------

! -------------------------------------------------------------------------
  SUBROUTINE d14co_free_memory

    ! 14CO MODULE ROUTINE (ECHAM-5 INTERFACE, PRIVATE)
    !
    ! free memory of global fields
    !
    ! Author: Patrick Joeckel, MPICH, Feb 2002

    IMPLICIT NONE

    ! I/O

    ! LOCAL
    !CHARACTER(LEN=*), PARAMETER :: substr = 'd14co_free_memory'
    INTEGER :: n

    set_loop: DO n=1, NMAXSETS

       SELECT CASE(XSET(n)%q_is_what)
       CASE (IS_GPTRACER)
          ! NOT IMPLEMENTED
!!#D attila +
#ifdef ECHAM5 
      CASE (IS_LGTRACER)
          ! NOT IMPLEMENTED
#endif
!!#D attila -
       CASE (IS_GPCHANNEL, IS_GPOFFLINE)
          NULLIFY(XSET(n)%import_q)
          IF (ASSOCIATED(XSET(n)%paxis_q)) DEALLOCATE(XSET(n)%paxis_q)
          NULLIFY(XSET(n)%paxis_q)
!!#D attila +
#ifdef ECHAM5 
          IF (ASSOCIATED(XSET(n)%q_lg)) DEALLOCATE(XSET(n)%q_lg)
          NULLIFY(XSET(n)%q_lg)
       CASE (IS_LGCHANNEL)
          IF (ASSOCIATED(XSET(n)%q_gp)) DEALLOCATE(XSET(n)%q_gp)
          NULLIFY(XSET(n)%q_gp)
#endif
!!#D attila -
       CASE DEFAULT
          ! ERROR
       END SELECT

       SELECT CASE (XSET(n)%ohs_is_what)
       CASE (IS_GPOFFLINE)
          NULLIFY(XSET(n)%import_q)
          IF (ASSOCIATED(XSET(n)%paxis_q)) DEALLOCATE(XSET(n)%paxis_q)
          NULLIFY(XSET(n)%paxis_q)
!!#D attila +
#ifdef ECHAM5 
          IF (ASSOCIATED(XSET(n)%ohs_lg)) DEALLOCATE(XSET(n)%ohs_lg)
          NULLIFY(XSET(n)%ohs_lg)
#endif
!!#D attila -
       CASE (IS_GPTRACER)
!!#D attila +
#ifdef ECHAM5 
          IF (ASSOCIATED(XSET(n)%ohs_lg)) DEALLOCATE(XSET(n)%ohs_lg)
          NULLIFY(XSET(n)%ohs_lg)
#endif
!!#D attila -
          IF (ASSOCIATED(XSET(n)%ohs_gp)) DEALLOCATE(XSET(n)%ohs_gp)
          NULLIFY(XSET(n)%ohs_gp)
       CASE (IS_GPCHANNEL)
!!#D attila +
#ifdef ECHAM5 
          IF (ASSOCIATED(XSET(n)%ohs_lg)) DEALLOCATE(XSET(n)%ohs_lg)
          NULLIFY(XSET(n)%ohs_lg)
       CASE (IS_LGTRACER)
          IF (ASSOCIATED(XSET(n)%ohs_gp)) DEALLOCATE(XSET(n)%ohs_gp)
          NULLIFY(XSET(n)%ohs_gp)
          IF (ASSOCIATED(XSET(n)%ohs_lg)) DEALLOCATE(XSET(n)%ohs_lg)
          NULLIFY(XSET(n)%ohs_lg)
       CASE (IS_LGCHANNEL)
          IF (ASSOCIATED(XSET(n)%ohs_gp)) DEALLOCATE(XSET(n)%ohs_gp)
          NULLIFY(XSET(n)%ohs_gp)
#endif
!!#D attila -
       CASE DEFAULT
          ! ERROR
       END SELECT

       SELECT CASE (XSET(n)%oht_is_what)
       CASE (IS_GPOFFLINE)
          NULLIFY(XSET(n)%import_q)
          IF (ASSOCIATED(XSET(n)%paxis_q)) DEALLOCATE(XSET(n)%paxis_q)
          NULLIFY(XSET(n)%paxis_q)
!!#D attila +
#ifdef ECHAM5 
          IF (ASSOCIATED(XSET(n)%oht_lg)) DEALLOCATE(XSET(n)%oht_lg)
          NULLIFY(XSET(n)%oht_lg)
#endif
!!#D attila -
       CASE (IS_GPTRACER)
!!#D attila +
#ifdef ECHAM5 
          IF (ASSOCIATED(XSET(n)%oht_lg)) DEALLOCATE(XSET(n)%oht_lg)
          NULLIFY(XSET(n)%oht_lg)
#endif
!!#D attila -
          IF (ASSOCIATED(XSET(n)%oht_gp)) DEALLOCATE(XSET(n)%oht_gp)
          NULLIFY(XSET(n)%oht_gp)
       CASE (IS_GPCHANNEL)
!!#D attila +
#ifdef ECHAM5 
          IF (ASSOCIATED(XSET(n)%oht_lg)) DEALLOCATE(XSET(n)%oht_lg)
          NULLIFY(XSET(n)%oht_lg)
       CASE (IS_LGTRACER)
          IF (ASSOCIATED(XSET(n)%oht_gp)) DEALLOCATE(XSET(n)%oht_gp)
          NULLIFY(XSET(n)%oht_gp)
          IF (ASSOCIATED(XSET(n)%oht_lg)) DEALLOCATE(XSET(n)%oht_lg)
          NULLIFY(XSET(n)%oht_lg)
       CASE (IS_LGCHANNEL)
          IF (ASSOCIATED(XSET(n)%oht_gp)) DEALLOCATE(XSET(n)%oht_gp)
          NULLIFY(XSET(n)%oht_gp)
#endif
!!#D attila -
       CASE DEFAULT
          ! ERROR
       END SELECT

    END DO set_loop

  END SUBROUTINE d14co_free_memory
! -------------------------------------------------------------------------

! ************************************************************************
! PRIVATE ECHAM-5 INTERFACE ROUTINES
! ************************************************************************

! ----------------------------------------------------------------------
  SUBROUTINE d14co_read_nml_cpl(status, iou)

    ! 14CO MODULE ROUTINE (ECHAM-5 INTERFACE, PRIVATE)
    !
    ! read namelist for 'coupling' to online tracers
    !
    ! Author: Patrick Joeckel, MPICH, Feb 2002

    USE messy_main_tools, ONLY: read_nml_open, read_nml_check, read_nml_close
    USE messy_main_tracer_mem_bi, ONLY: NGCELL

    IMPLICIT NONE

    ! I/O
    INTEGER, INTENT(OUT) :: status     ! error status
    INTEGER, INTENT(IN)  :: iou        ! I/O unit

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'd14co_read_nml_cpl'
    LOGICAL              :: lex      ! file exists ?
    INTEGER              :: fstat    ! file status
    INTEGER              :: n
!!#D attila +
#ifdef ECHAM5
    LOGICAL              :: l_lg
#endif
!!#D attila -
    
    NAMELIST /CPL/ &
          S_SWITCH, S_14CO, S_OHt, S_OHs, S_TP_STE, S_TP_OH

!!#D attila +
#ifdef ECHAM5
    NAMELIST /CPL_LG/ &
           C_LG_CHANNEL, C_LG_ILON, C_LG_ILAT, C_LG_PRESS
#endif
!!#D attila -
    
    status = 1

    CALL read_nml_open(lex, substr, iou, 'CPL', modstr)
    IF (.not.lex) RETURN    ! <modstr>.nml does not exist

    READ(iou, NML=CPL, IOSTAT=fstat)
    CALL read_nml_check(fstat, substr, iou, 'CPL', modstr)
    IF (fstat /= 0) RETURN  ! error while reading namelist

    ! CONSISTENCY CHECK
    IF (NGCELL <= 0) THEN
!!#D attila +
#ifdef ECHAM5 
       WRITE(*,*) '!!! NO LAGRANGIAN SCHEME ACTIVATED !!!'
       WRITE(*,*) '!!! SWITCHING OFF LG-INTEGRATION   !!!'
#endif
!!#D attila -
       DO n=1, NMAXSETS
          S_SWITCH(n)%lg_lint = .FALSE.
       END DO
    END IF

    CALL read_nml_close(substr, iou, modstr)

!!#D attila +
#ifdef ECHAM5
    l_lg = .FALSE.
    DO n=1, NMAXSETS
       l_lg = l_lg .OR. S_SWITCH(n)%lg_lint
    END DO
    IF (l_lg) THEN
       CALL read_nml_open(lex, substr, iou, 'CPL_LG', modstr)
       IF (.not.lex) RETURN    ! <modstr>.nml does not exist
       !
       READ(iou, NML=CPL_LG, IOSTAT=fstat)
       CALL read_nml_check(fstat, substr, iou, 'CPL_LG', modstr)
       IF (fstat /= 0) RETURN  ! error while reading namelist
       !
       CALL read_nml_close(substr, iou, modstr)
    END IF
#endif
!!#D attila -
    
    status = 0  ! no ERROR

  END SUBROUTINE d14co_read_nml_cpl
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
  SUBROUTINE d14co_check_unit(callstr, mode, unit, lconv)

    USE messy_main_blather_bi,   ONLY: error_bi

    IMPLICIT NONE

    ! I/O
    CHARACTER(LEN=*), INTENT(IN)  :: callstr
    CHARACTER(LEN=2), INTENT(IN)  :: mode
    CHARACTER(LEN=*), INTENT(IN)  :: unit
    LOGICAL,          INTENT(OUT) :: lconv

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr='d14co_check_unit'

    SELECT CASE(mode)
    CASE('OH')
       !
       SELECT CASE(TRIM(unit))
       CASE('mol/mol','vmr')
          lconv = .TRUE.
       CASE('cm^-3','cm-3','cm^(-3)')
          lconv = .FALSE.
       CASE DEFAULT
          CALL error_bi(' unknown unit for OH distribution !', substr)
       END SELECT
       !
    CASE('CS')
       !
       SELECT CASE(TRIM(unit))
       CASE('molec/(g s)')
          lconv = .TRUE.
       CASE DEFAULT
          CALL error_bi(' unknown unit for 14CO source !', substr)
       END SELECT
       !
    CASE DEFAULT
       !
       CALL error_bi(' called with unknown mode from '//TRIM(callstr), substr)
       !
    END SELECT

  END SUBROUTINE d14co_check_unit
! ----------------------------------------------------------------------

! **********************************************************************
END MODULE messy_d14co_si
! **********************************************************************

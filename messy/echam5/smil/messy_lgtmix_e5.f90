! ************************************************************************
MODULE messy_lgtmix_e5
! ************************************************************************

  ! AUTHOR:
  !  Patrick Joeckel, MPICH, January 2004

  ! ECHAM5/MESSy
  USE messy_main_blather_bi,    ONLY: start_message_bi, end_message_bi
  ! MESSy
  USE messy_main_channel,       ONLY: STRLEN_CHANNEL, STRLEN_OBJECT
  USE messy_lgtmix
  USE messy_main_constants_mem, ONLY: STRLEN_MEDIUM
  USE messy_main_tools,         ONLY: PTR_1D_ARRAY, PTR_3D_ARRAY

  IMPLICIT NONE
  PRIVATE

  INTRINSIC :: TRIM, ABS, EXP, MAX, MIN, ADJUSTL, NINT, NULL &
       , MINVAL, MAXVAL, MERGE

  ! ### STRUCTURES FOR MIXING LAYER INFORMATION #############################
  ! MAX. NUMBER OF MIXING LAYERS IN NAMELIST
  INTEGER, PARAMETER :: NMAXMIXL        = 10

  ! MIXING LAYER INFORMATION
  TYPE IO_MIXL
     CHARACTER(LEN=STRLEN_CHANNEL) :: cha_l1 = '' ! name of layer 1: channel
     CHARACTER(LEN=STRLEN_OBJECT)  :: obj_l1 = '' ! name of layer 1: object
     CHARACTER(LEN=STRLEN_CHANNEL) :: cha_l2 = '' ! name of layer 2: channel
     CHARACTER(LEN=STRLEN_OBJECT)  :: obj_l2 = '' ! name of layer 2: object
     CHARACTER(LEN=STRLEN_CHANNEL) :: cha_mp = '' ! name of mix. parameter: ...
     CHARACTER(LEN=STRLEN_OBJECT)  :: obj_mp = '' ! ... channel , object
     REAL(DP)                      :: mmin = 0.0  ! value for min. mixing (0.0)
     REAL(DP)                      :: mmax = 0.0  ! value for max. mixing (1.0)
  END TYPE IO_MIXL

  TYPE MIXL
     TYPE(IO_MIXL) :: io
     REAL(DP), DIMENSION(:,:),   POINTER :: fl1 => NULL() ! index layer 1
     REAL(DP), DIMENSION(:,:),   POINTER :: fl2 => NULL() ! index layer 2
     REAL(DP), DIMENSION(:,:,:), POINTER :: fmp => NULL() ! mixing parameter
     INTEGER                             :: i1 = -1 ! const. index level 1
     INTEGER                             :: i2 = -1 ! const. index level 2
     REAL(DP)                            :: rmp = -1.0 ! const mixing parameter
  END TYPE MIXL

  INTEGER                                   :: NMX
  TYPE (MIXL),    DIMENSION(NMAXMIXL), SAVE :: XMX
  ! ##########################################################################

  ! ### STRUCTURES FOR TRACER SPECIFIC SCALING  ##############################
  ! MAX. NUMBER OF MIXING LAYERS IN NAMELIST
  TYPE IO_TMIX
     CHARACTER(LEN=STRLEN_MEDIUM)     :: basename    = '' ! name of tracer
     CHARACTER(LEN=STRLEN_MEDIUM)     :: subname     = '' ! OPTIONAL subname
     REAL(DP), DIMENSION(NMAXMIXL)    :: scale       = 1.0_dp
  END TYPE IO_TMIX

  ! MAX. NUMBER OF TRACERS IN NAMELIST
  INTEGER, PARAMETER :: NMAXTRAC = 100

  ! TRACER SCALING MATRIX (TRACER ID x MIXING LEVEL)
  REAL(DP), DIMENSION(:,:), ALLOCATABLE, SAVE :: TSCAL
  ! ##########################################################################

  ! CPL NAMELIST ENTRIES
  ! FORCE GP channel object of mixing parameter
  LOGICAL                                   :: l_force
  TYPE (IO_MIXL), DIMENSION(NMAXMIXL), SAVE :: MX
  TYPE (IO_TMIX), DIMENSION(NMAXTRAC), SAVE :: TX

  ! CALCULATION REQUIRED ?
  LOGICAL :: l_lgtmix_switch = .TRUE.

  ! POINTER TO MIXING PARAMETER [0,1]
  REAL(DP), DIMENSION(:,:,:), POINTER :: fmix_gp => NULL() ! GP REPRESENTATION
  REAL(DP), DIMENSION(:),     POINTER :: fmix_lg => NULL() ! LG REPRESENTATION
  ! MIXING LAYER FLAGS
  TYPE(PTR_3D_ARRAY), DIMENSION(:), POINTER :: fflag_gp => NULL()
  TYPE(PTR_1D_ARRAY), DIMENSION(:), POINTER :: fflag_lg => NULL()
  ! TRACER SPECIFIC COMBINED SCALING (LG)
  REAL(DP), DIMENSION(:), POINTER :: fscal_lg => NULL() 

  ! BACKGROUND MIXING RATIO
  REAL(DP), DIMENSION(:,:,:), POINTER :: bg_gp
  REAL(DP), DIMENSION(:),     POINTER :: bg_lg
  ! LAGRANGIAN TRACER START VALUE
  REAL(DP), DIMENSION(:),     POINTER :: trac_lg

  PUBLIC :: lgtmix_initialize
  PUBLIC :: lgtmix_init_memory
  PUBLIC :: lgtmix_init_coupling
  PUBLIC :: lgtmix_global_end
  PUBLIC :: lgtmix_free_memory
  !PRIVATE :: lgtmix_read_nml_cpl
  !PRIVATE :: str2idx
  !PRIVATE :: str2r

CONTAINS

! -----------------------------------------------------------------------
  SUBROUTINE lgtmix_initialize

    ! ECHAM5/MESSy
    USE messy_main_mpi_bi,     ONLY: p_parallel_io, p_io, p_bcast
    USE messy_main_blather_bi, ONLY: error_bi 
    USE messy_main_tools,      ONLY: find_next_free_unit

    IMPLICIT NONE

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER   :: substr = 'lgtmix_initialize'
    INTEGER                       :: iou    ! I/O unit
    INTEGER                       :: status ! error status
    INTEGER                       :: i

    ! READ NAMELIST CPL
    IF (p_parallel_io) THEN
       iou = find_next_free_unit(100,200)
       CALL lgtmix_read_nml_cpl(status, iou)
       IF (status /= 0) CALL error_bi('', substr)
    END IF
    CALL p_bcast(l_force, p_io)

    CALL start_message_bi(modstr,'INITIALISATION ',substr)

    IF (p_parallel_io) THEN
       NMX = 1
       layer_loop: DO i=1, NMAXMIXL
          IF (TRIM(MX(i)%cha_l1) == '') CYCLE
          !
          WRITE(*,*) 'ANALYZING MXING LAYER ------------------------------'
          !
          ! LEVEL 1
          XMX(NMX)%io%cha_l1 = TRIM(ADJUSTL(MX(i)%cha_l1))
          CALL str2idx(XMX(NMX)%io%cha_l1, XMX(NMX)%i1, status)
          IF (status == 0) THEN
             WRITE(*,*) '  FROM LEVEL (INDEX)       : ',XMX(NMX)%i1
          ELSE
             WRITE(*,*) '  FROM LEVEL (CHANNEL)     : ',XMX(NMX)%io%cha_l1
             XMX(NMX)%io%obj_l1 = TRIM(ADJUSTL(MX(i)%obj_l1))
             IF (TRIM(XMX(NMX)%io%obj_l1) == '') THEN
                WRITE(*,*) '  ... empty channel object ... skipping ...'
                CYCLE
             ELSE
                WRITE(*,*) '             (OBJECT)      : ',XMX(NMX)%io%obj_l1
             END IF
          END IF
          !
          ! LEVEL 2
          XMX(NMX)%io%cha_l2 = TRIM(ADJUSTL(MX(i)%cha_l2))
          CALL str2idx(XMX(NMX)%io%cha_l2, XMX(NMX)%i2, status)
          IF (status == 0) THEN
             WRITE(*,*) '    TO LEVEL (INDEX)       : ',XMX(NMX)%i2
          ELSE
             WRITE(*,*) '    TO LEVEL (CHANNEL)     : ',XMX(NMX)%io%cha_l2
             XMX(NMX)%io%obj_l2 = TRIM(ADJUSTL(MX(i)%obj_l2))
             IF (TRIM(XMX(NMX)%io%obj_l2) == '') THEN
                WRITE(*,*) '  ... empty channel object ... skipping ...'
                CYCLE
             ELSE
                WRITE(*,*) '             (OBJECT)      : ',XMX(NMX)%io%obj_l2
             END IF
          END IF
          !
          ! MIXING PARAMETER
          XMX(NMX)%io%cha_mp = TRIM(ADJUSTL(MX(i)%cha_mp))
          CALL str2r(XMX(NMX)%io%cha_mp, XMX(NMX)%rmp, status)
          IF (status == 0) THEN
             IF ((XMX(NMX)%rmp < 0.0) .OR. (XMX(NMX)%rmp > 1.0)) THEN
                WRITE(*,*) '  ... mixing strength not in interval [0.0,1.0]'&
                     &//' ... skipping ...'
                CYCLE
             END IF
             WRITE(*,*) '  MIXING STRENGTH (VALUE)  : ',XMX(NMX)%rmp
          ELSE
             WRITE(*,*) '  MIXING STRENGTH (CHANNEL): ',XMX(NMX)%io%cha_mp
             IF (TRIM(MX(i)%obj_mp) == '') THEN
                WRITE(*,*) '  ... empty channel object ... skipping ...'
                CYCLE
             ELSE
                XMX(NMX)%io%obj_mp = TRIM(ADJUSTL(MX(i)%obj_mp))
                WRITE(*,*) '                  (OBJECT) : ',XMX(NMX)%io%obj_mp
                ! 
                ! MIXING PARAMETER RANGE
                XMX(NMX)%io%mmin = MX(i)%mmin
                WRITE(*,*) '    MIN. MIXING (0.0) AT   : ',XMX(NMX)%io%mmin
                XMX(NMX)%io%mmax = MX(i)%mmax
                WRITE(*,*) '    MAX. MIXING (1.0) AT   : ',XMX(NMX)%io%mmax
             END IF
          END IF
          !
          NMX = NMX + 1
       END DO layer_loop
       NMX = NMX - 1
       WRITE(*,*) '----------------------------------------------------'
    END IF
    CALL p_bcast(NMX, p_io)
    
    ! BROADCAST RESULTS
    DO i=1, NMX
       CALL p_bcast(XMX(i)%io%cha_l1, p_io)
       CALL p_bcast(XMX(i)%io%obj_l1, p_io)
       CALL p_bcast(XMX(i)%io%cha_l2, p_io)
       CALL p_bcast(XMX(i)%io%obj_l2, p_io)
       CALL p_bcast(XMX(i)%io%cha_mp, p_io)
       CALL p_bcast(XMX(i)%io%obj_mp, p_io)
       CALL p_bcast(XMX(i)%io%mmin, p_io)
       CALL p_bcast(XMX(i)%io%mmax, p_io)
       !
       CALL p_bcast(XMX(i)%i1, p_io)
       CALL p_bcast(XMX(i)%i2, p_io)
       CALL p_bcast(XMX(i)%rmp, p_io)
    END DO

    IF (p_parallel_io) THEN
       WRITE(*,*) ' ---> ',NMX,' MIXING LAYERS IDENTIFIED !'
    END IF

    DO i=1, NMAXTRAC
       CALL p_bcast(TX(i)%basename, p_io)
       CALL p_bcast(TX(i)%subname,  p_io)
       CALL p_bcast(TX(i)%scale(:), p_io)
    END DO

    CALL end_message_bi(modstr,'INITIALISATION ',substr)

  END SUBROUTINE lgtmix_initialize
! -----------------------------------------------------------------------

! -----------------------------------------------------------------------
  SUBROUTINE lgtmix_init_memory

    ! LGTMIX MODULE ROUTINE (ECHAM-5 INTERFACE, PRIVATE)
    !
    ! define LGTMIX specific channel(s) and allocate memory for
    ! global fields
    !
    ! Author: Patrick Joeckel, MPICH, Jan 2004

    ! ECHAM5/MESSy
    USE messy_main_mpi_bi,           ONLY: p_parallel_io, message
    USE messy_main_tracer_mem_bi,    ONLY: ntrac_lg, ti_lg, NGCELL
    USE messy_main_channel_error_bi, ONLY: channel_halt
    USE messy_main_channel_bi,       ONLY: GP_3D_MID, LG_ATTILA
    ! MESSy
    USE messy_main_channel,       ONLY: new_channel, new_channel_object &
                                      , new_attribute
    USE messy_main_tracer,        ONLY: ON, I_mix

    IMPLICIT NONE

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'lgtmix_init_memory'
    INTEGER                     :: status, jt
    INTEGER                     :: i
    CHARACTER(LEN=10)           :: name = ''

    CALL start_message_bi(modstr,'CHANNEL DEFINITION',substr)

    ! CHECK IF LAGRANGIAN TRANSPORT IS ACTIVE
    IF (NGCELL <= 0) THEN
       CALL message(substr,' WARNING! ATTILA NOT ACTIVE!')
       l_lgtmix_switch = .FALSE.
       IF (.NOT. l_force) RETURN
    END IF

    ! CHECK EXISTENCE OF LG-TRACERS
    IF (ntrac_lg  <=0) THEN
       CALL message(substr,' WARNING! NO LG-TRACERS!')
       l_lgtmix_switch = .FALSE.
       IF (.NOT. l_force) RETURN
    END IF
    
    ! CHECK IF MIXING IS REQUESTED FOR LG-TRACERS
    l_lgtmix_switch = .FALSE.
    DO jt = 1, ntrac_lg
       IF (ti_lg(jt)%tp%meta%cask_i(I_mix) == ON) l_lgtmix_switch = .TRUE.
    END DO
    !
    IF (.NOT. l_lgtmix_switch) THEN
       CALL message(substr,' WARNING! NO LG-MIXING REQUESTED!')
       IF (.NOT. l_force) RETURN
    END IF

    IF (l_lgtmix_switch .OR. l_force) THEN

       ! define new channel
       CALL new_channel(status, modstr)
       CALL channel_halt(substr, status)

       ! LGTMIX DIAGNOSTIC OUTPUT
       IF (p_parallel_io) THEN
          WRITE(*,*) 'new channel objects ...'
          WRITE(*,*) ' ... FMIX_GP'
       END IF
       CALL new_channel_object(status, modstr, 'FMIX_GP' &
            , p3=fmix_gp, reprid=GP_3D_MID, lrestreq = .TRUE.)
       CALL channel_halt(substr, status)
       CALL new_attribute(status, modstr, 'FMIX_GP' &
            , 'long_name', c='mixing coefficient (GP)' )
       CALL channel_halt(substr, status)
       CALL new_attribute(status, modstr, 'FMIX_GP' &
            , 'units', c=' ' )
       CALL channel_halt(substr, status)

       ! FLAGS FOR DIFFREENT MIXING LAYERS
       ALLOCATE(fflag_gp(NMX))
       name = 'FLAG_GP_'
       DO i=1, NMX
          WRITE(name(9:10),'(i2.2)') i
          IF (p_parallel_io) THEN
             WRITE(*,*) ' ... ',TRIM(name)
          END IF
          CALL new_channel_object(status, modstr, TRIM(name) &
               , p3=fflag_gp(i)%ptr, reprid=GP_3D_MID, lrestreq = .TRUE.)
          CALL channel_halt(substr, status)
          CALL new_attribute(status, modstr, TRIM(name) &
               , 'long_name', c='mixing layer flag (GP)' )
          CALL channel_halt(substr, status)
          CALL new_attribute(status, modstr, TRIM(name) &
               , 'units', c=' ' )
          CALL channel_halt(substr, status)
       END DO
    END IF

    IF (l_lgtmix_switch) THEN

       IF (p_parallel_io) WRITE(*,*) ' ... FMIX_LG'
       CALL new_channel_object(status, modstr, 'FMIX_LG' &
            , p1=fmix_lg, reprid=LG_ATTILA, lrestreq=.TRUE.)
       CALL channel_halt(substr, status)
       CALL new_attribute(status, modstr, 'FMIX_LG' &
            , 'long_name', c='mixing coefficient (LG)' )
       CALL channel_halt(substr, status)
       CALL new_attribute(status, modstr, 'FMIX_LG' &
            , 'units', c=' ' )
       CALL channel_halt(substr, status)

       ! FLAGS FOR DIFFREENT MIXING LAYERS
       ALLOCATE(fflag_lg(NMX))
       name = 'FLAG_LG_'
       DO i=1, NMX
          WRITE(name(9:10),'(i2.2)') i
          IF (p_parallel_io) THEN
             WRITE(*,*) ' ... ',TRIM(name)
          END IF
          CALL new_channel_object(status, modstr, TRIM(name) &
               , p1=fflag_lg(i)%ptr, reprid=LG_ATTILA, lrestreq = .TRUE.)
          CALL channel_halt(substr, status)
          CALL new_attribute(status, modstr, TRIM(name) &
               , 'long_name', c='mixing layer flag (LG)' )
          CALL channel_halt(substr, status)
          CALL new_attribute(status, modstr, TRIM(name) &
               , 'units', c=' ' )
          CALL channel_halt(substr, status)
       END DO

       ! INTERNAL WORK-SPACE; NOT SUITED FOR OUTPUT:

!    IF (p_parallel_io) WRITE(*,*) ' ... TRAC_LG'
       ! (1) LAGRANGIAN TRACER START VALUE
       CALL new_channel_object(status, modstr, 'TRAC_LG' &
            , p1=trac_lg, reprid=LG_ATTILA)
       CALL channel_halt(substr, status)
       CALL new_attribute(status, modstr, 'TRAC_LG' &
            , 'long_name', c='lagrangian tracer' )
       CALL channel_halt(substr, status)
       CALL new_attribute(status, modstr, 'TRAC_LG' &
            , 'units', c='mol/mol' )
       CALL channel_halt(substr, status)

!    IF (p_parallel_io) WRITE(*,*) ' ... BG_GP'
       ! (2) BACKGROUND VALUE IN GRIDPOINT REPRESENTATION
       CALL new_channel_object(status, modstr, 'BG_GP' &
            , p3=bg_gp, reprid=GP_3D_MID)
       CALL channel_halt(substr, status)
       CALL new_attribute(status, modstr, 'BG_GP' &
            , 'long_name', c='gridpoint tracer' )
       CALL channel_halt(substr, status)
       CALL new_attribute(status, modstr, 'BG_GP' &
            , 'units', c='mol/mol' )
       CALL channel_halt(substr, status)

!    IF (p_parallel_io) WRITE(*,*) ' ... TRAC_LG'
       ! (3) BACKGROUND VALUE IN LAGRANGIAN REPRESENTATION
       CALL new_channel_object(status, modstr, 'BG_LG' &
            , p1=bg_lg, reprid=LG_ATTILA)
       CALL channel_halt(substr, status)
       CALL new_attribute(status, modstr, 'BG_LG' &
            , 'long_name', c='lagrangian tracer' )
       CALL channel_halt(substr, status)
       CALL new_attribute(status, modstr, 'BG_LG' &
            , 'units', c='mol/mol' )
       CALL channel_halt(substr, status)

       !    IF (p_parallel_io) WRITE(*,*) ' ... FSCAL_LG'
       ! (4) TRACER SPECIFIC SCALING IN LAGRANGIAN REPRESENTATION
       CALL new_channel_object(status, modstr, 'FSCAL_LG' &
            , p1=fscal_lg, reprid=LG_ATTILA)
       CALL channel_halt(substr, status)
       CALL new_attribute(status, modstr, 'FSCAL_LG' &
            , 'long_name', c='lagrangian tracer specific scaling' )
       CALL channel_halt(substr, status)
       CALL new_attribute(status, modstr, 'FSCAL_LG' &
            , 'units', c=' ' )
       CALL channel_halt(substr, status)

    END IF

    CALL end_message_bi(modstr,'CHANNEL DEFINITION',substr)
    
  END SUBROUTINE lgtmix_init_memory
! -----------------------------------------------------------------------

! ----------------------------------------------------------------------
  SUBROUTINE lgtmix_init_coupling

    ! LGTMIX MODULE ROUTINE (ECHAM-5 INTERFACE, PRIVATE)
    !
    ! initialize 'coupling' to online tracers/channels
    !
    ! Author: Patrick Joeckel, MPICH, Jan 2004

    ! ECHAM5/MESSy
    USE messy_main_mpi_bi,           ONLY: p_parallel_io
    USE messy_main_channel_error_bi, ONLY: channel_halt
    USE messy_main_tracer_mem_bi,    ONLY: ntrac_lg, LGTRSTR
    ! MESSy
    USE messy_main_channel,       ONLY: get_channel_object
    USE messy_main_tracer,        ONLY: get_tracer

    IMPLICIT NONE

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER       :: substr = 'lgtmix_init_coupling'
    INTEGER                           :: status
    INTEGER                           :: i, j
    INTEGER                           :: idt

    IF (.NOT. (l_lgtmix_switch .OR. l_force)) RETURN

    ! CHECK IF TRACER IS AVAILABLE
    CALL start_message_bi(modstr,'COUPLING INITIALIZATION',substr)

    layer_loop: DO i=1, NMX

       IF (p_parallel_io) THEN
          WRITE(*,*) 'MIXING LAYER NO. ',i,'----------------------'
       END IF

       ! CHECK 1st LEVEL
       IF (XMX(i)%i1 < 0) THEN
          IF (p_parallel_io) THEN
             WRITE(*,*) '  checking for level 1 ...'
             WRITE(*,*) '    channel        : ',TRIM(XMX(i)%io%cha_l1)
             WRITE(*,*) '    object         : ',TRIM(XMX(i)%io%obj_l1)
          END IF
          CALL get_channel_object(status &
               , TRIM(XMX(i)%io%cha_l1), TRIM(XMX(i)%io%obj_l1) &
               , p2=XMX(i)%fl1 )
          CALL channel_halt(substr, status)
       ELSE
          IF (p_parallel_io) THEN
             WRITE(*,*) '  checking for level 1 ...'
             WRITE(*,*) '    constant index : ',XMX(i)%i1
          END IF
       END IF
       
       ! CHECK 2nd LEVEL
       IF (XMX(i)%i2 < 0) THEN
          IF (p_parallel_io) THEN
             WRITE(*,*) '  checking for level 2 ...'
             WRITE(*,*) '    channel        : ',TRIM(XMX(i)%io%cha_l2)
             WRITE(*,*) '    object         : ',TRIM(XMX(i)%io%obj_l2)
          END IF
          CALL get_channel_object(status &
               , TRIM(XMX(i)%io%cha_l2), TRIM(XMX(i)%io%obj_l2) &
               , p2=XMX(i)%fl2 )
          CALL channel_halt(substr, status)
       ELSE
          IF (p_parallel_io) THEN
             WRITE(*,*) '  checking for level 2 ...'
             WRITE(*,*) '    constant index : ',XMX(i)%i2
          END IF
       END IF

       ! CHECK PARAMETER
       IF (XMX(i)%rmp < 0.0) THEN
          IF (p_parallel_io) THEN
             WRITE(*,*) '  checking for mixing parameter ...'
             WRITE(*,*) '    channel        : ',TRIM(XMX(i)%io%cha_mp)
             WRITE(*,*) '    object         : ',TRIM(XMX(i)%io%obj_mp)
          END IF
          CALL get_channel_object(status &
               , TRIM(XMX(i)%io%cha_mp), TRIM(XMX(i)%io%obj_mp) &
               , p3=XMX(i)%fmp )
          CALL channel_halt(substr, status)
       ELSE
          IF (p_parallel_io) THEN
             WRITE(*,*) '  constant mixing parameter: ',XMX(i)%rmp
          END IF
       END IF

       IF (p_parallel_io) THEN
          WRITE(*,*) '--------------------------------------------------'
       END IF

    END DO layer_loop

    ALLOCATE(TSCAL(ntrac_lg,NMAXMIXL))
    TSCAL(:,:) = 1.0_dp

    tracer_in_cpl_nml: DO i=1, NMAXTRAC

       IF (TRIM(TX(i)%basename) == '') CYCLE
       CALL get_tracer(status, LGTRSTR, TRIM(TX(i)%basename) &
            , TRIM(TX(I)%subname), idt)
       IF (status == 0) THEN
!!$          DO j=1, NMX
!!$             IF ( (TX(i)%scale(j) < 0.0_dp) .OR. &
!!$                  (TX(i)%scale(j) > 1.0_dp) ) THEN
!!$                IF (p_parallel_io) THEN
!!$                   IF (TRIM(TX(I)%subname) == '') THEN
!!$                      WRITE(*,*) 'skipping scaling for tracer ' &
!!$                           ,TRIM(TX(i)%basename) &
!!$                           ,' since scaling out of range [0,1] '
!!$                   ELSE
!!$                      WRITE(*,*) 'skipping scaling for tracer ' &
!!$                           ,TRIM(TX(i)%basename)//'_'//TRIM(TX(I)%subname) &
!!$                           ,' since scaling out of range [0,1] '
!!$                   END IF
!!$                END IF
!!$                CYCLE
!!$             END IF
!!$          END DO
          !
          TSCAL(idt,1:NMX) = TX(i)%scale(1:NMX)
          !
          IF (p_parallel_io) THEN
             IF (TRIM(TX(I)%subname) == '') THEN
                WRITE(*,*) 'MIXING FOR TRACER ',TRIM(TX(i)%basename) &
                     ,' SCALED BY ', TSCAL(idt,1:NMX)
             ELSE
                WRITE(*,*) 'MIXING FOR TRACER '&
                     ,TRIM(TX(i)%basename)//'_'//TRIM(TX(I)%subname) &
                     ,' SCALED BY ', TSCAL(idt,1:NMX)
             END IF
          END IF
       END IF

    END DO tracer_in_cpl_nml

    IF (p_parallel_io) THEN
       WRITE(*,*) '--------------------------------------------------'
    END IF

    CALL end_message_bi(modstr,'COUPLING INITIALIZATION',substr)

  END SUBROUTINE lgtmix_init_coupling
! ----------------------------------------------------------------------

! -----------------------------------------------------------------------
  SUBROUTINE lgtmix_global_end

    USE messy_attila_tools_e5,    ONLY: gp2lg_e5, lg2gp_e5, LG2GP_AVE
    USE messy_main_timer,         ONLY: ztmst=>time_step_len
    USE messy_main_grid_def_mem_bi,ONLY: nlev, nproma, npromz, ngpblks
    USE messy_main_tracer_mem_bi, ONLY: qxtm1_a, ti_lg, ntrac_lg, qxtte_a 
    USE messy_main_blather_bi,    ONLY: error_bi
    USE messy_main_tracer,        ONLY: ON, I_mix

    IMPLICIT NONE

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER     :: substr = 'lgtmix_global_end'
    INTEGER                         :: zjrow, zkproma, jp, jt
    INTEGER                         :: i
    INTEGER                         :: jk1h, jk2h, jk1, jk2

    ! CHECK GLOBAL SWITCH
    IF (.NOT. (l_lgtmix_switch .OR. l_force)) RETURN

    ! CALULATE MIXING COEFFICIENT ...
    fmix_gp(:,:,:) = 0.0_DP

    layer_loop: DO i=1, NMX

       ! ... AND LAYER FLAG 
       fflag_gp(i)%ptr(:,:,:) = 0.0_dp

       local_loop: DO zjrow=1, ngpblks
          
          IF ( zjrow == ngpblks ) THEN
             zkproma = npromz
          ELSE
             zkproma = nproma
          END IF

          vector_loop: DO jp=1, zkproma

             ! INDEX OF 1st LEVEL
             IF (XMX(i)%i1 < 0) THEN
                jk1h = NINT(XMX(i)%fl1(jp,zjrow))
             ELSE
                jk1h = XMX(i)%i1
             END IF
             
             ! INDEX OF 2nd LEVEL
             IF (XMX(i)%i2 < 0) THEN
                jk2h = NINT(XMX(i)%fl2(jp,zjrow))
             ELSE
                jk2h = XMX(i)%i2
             END IF

             ! SORT LEVELS
             jk1 = MIN(jk1h, jk2h)
             jk2 = MAX(jk1h, jk2h)
       
             ! LIMIT LEVELS TO [1, nlev]
             jk1 = MIN(jk1, nlev)
             jk1 = MAX(jk1, 1)
             jk2 = MIN(jk2, nlev)
             jk2 = MAX(jk2, 1)

             ! FIELD
             IF (XMX(i)%rmp < 0.0) THEN
                fmix_gp(jp,jk1:jk2,zjrow) = XMX(i)%fmp(jp,jk1:jk2,zjrow)
                ! LINEAR SCALING
                fmix_gp(jp,jk1:jk2,zjrow) =        &
                     MIN(fmix_gp(jp,jk1:jk2,zjrow) &
                     , MAX(XMX(i)%io%mmin, XMX(i)%io%mmax))
                fmix_gp(jp,jk1:jk2,zjrow) =        &
                     MAX(fmix_gp(jp,jk1:jk2,zjrow) &
                     , MIN(XMX(i)%io%mmin, XMX(i)%io%mmax))
                fmix_gp(jp,jk1:jk2,zjrow) =                &
                     (fmix_gp(jp,jk1:jk2,zjrow) -          &
                     MIN(XMX(i)%io%mmin, XMX(i)%io%mmax))  &
                     /ABS(XMX(i)%io%mmax - XMX(i)%io%mmin)
                IF (XMX(i)%io%mmin > XMX(i)%io%mmax) &
                     fmix_gp(jp,jk1:jk2,zjrow) = &
                     1.0_DP - fmix_gp(jp,jk1:jk2,zjrow)
             ELSE
                fmix_gp(jp,jk1:jk2,zjrow) = XMX(i)%rmp
             END IF

             ! FLAGS
             fflag_gp(i)%ptr(jp,jk1:jk2,zjrow) = 1.0_dp

          END DO vector_loop

       END DO local_loop

    END DO layer_loop

    IF (.NOT. l_lgtmix_switch) RETURN

    ! TRANSFORM GRIDPOINT -> LAGRANGE
    CALL gp2lg_e5(fmix_gp, fmix_lg)

    ! FLAGS
    DO i=1, NMX
       CALL gp2lg_e5(fflag_gp(i)%ptr, fflag_lg(i)%ptr)
    END DO

    ! APPLY TRACER MIXING
    tracer_loop: DO jt = 1, ntrac_lg
       IF (ti_lg(jt)%tp%meta%cask_i(I_mix) /= ON) CYCLE  ! NO MIXING REQUESTED

       ! INITIAL VALUE FOR INTEGRATION: <T-1> + <TENDENCY> * dt
       trac_lg(:) = qxtm1_a(:,jt) + qxtte_a(:,jt)*ztmst

       ! TRANSFORM LAGRANGE -> GRIDPOINT
       ! BACKGROUND = AVERAGE OVER ALL LG CELLS
       CALL lg2gp_e5(trac_lg, bg_gp, LG2GP_AVE)
       ! TRANSFORM BACKGROUND BACK: GRIDPOINT -> LAGRANGE
       CALL gp2lg_e5(bg_gp, bg_lg, lmcons=.FALSE.)

       ! TRACER SPECIFIC MIXING SCALING
       fscal_lg(:) = 0.0_dp
       ! NOTE: loop backwards to be consistent with loop for fmix_gp above;
       !       the vertical index ranges of different layers might overlap,
       !       and the last layer determines fmix_gp; here first comes first
       !       serves ...
       DO i=NMX, 1, -1
!qqq          fscal_lg(:) = fscal_lg(:) + fflag_lg(i)%ptr(:) * TSCAL(jt,i)
          fscal_lg(:) = MERGE(fscal_lg(:), fflag_lg(i)%ptr(:)*TSCAL(jt,i) &
               , fscal_lg(:) > 0.0_dp)
       END DO

       ! CHECK
       IF ( (MINVAL(fscal_lg*fmix_lg) < 0.0_dp) .OR. &
            (MAXVAL(fscal_lg*fmix_lg) > 1.0_dp) ) THEN
          CALL error_bi('MIXING STRENGTH OUT OF RANGE [0,1]' ,substr)
       END IF

       ! UPDATE TENDENCY
       qxtte_a(:,jt) = qxtte_a(:,jt) + &
!qqq            ( (bg_lg(:)-trac_lg(:)) * fmix_lg(:) ) / ztmst
            ( (bg_lg(:)-trac_lg(:)) * fmix_lg(:)*fscal_lg(:) ) / ztmst

    END DO tracer_loop

  END SUBROUTINE lgtmix_global_end
! -----------------------------------------------------------------------

! -----------------------------------------------------------------------
  SUBROUTINE lgtmix_free_memory

    IMPLICIT NONE
    
    DEALLOCATE(TSCAL)

    IF (ASSOCIATED(fflag_gp)) THEN
       DEALLOCATE(fflag_gp) ; NULLIFY(fflag_gp)
    END IF


  END SUBROUTINE lgtmix_free_memory
! -----------------------------------------------------------------------

! ************************************************************************
! PRIVATE ECHAM-5 INTERFACE ROUTINES
! ************************************************************************

! ----------------------------------------------------------------------
  SUBROUTINE lgtmix_read_nml_cpl(status, iou)

    ! LGTMIX MODULE ROUTINE (ECHAM-5 INTERFACE, PRIVATE)
    !
    ! read namelist for 'coupling' to channels
    !
    ! Author: Patrick Joeckel, MPICH, Jan 2004

    ! MESSy
    USE messy_main_tools, ONLY: read_nml_open, read_nml_check, read_nml_close

    IMPLICIT NONE

    ! I/O
    INTEGER, INTENT(OUT) :: status     ! error status
    INTEGER, INTENT(IN)  :: iou        ! I/O unit

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'lgtmix_read_nml_cpl'
    LOGICAL              :: lex      ! file exists ?
    INTEGER              :: fstat    ! file status

    NAMELIST /CPL/ l_force, MX, TX

    status = 1

    CALL read_nml_open(lex, substr, iou, 'CPL', modstr)
    IF (.not.lex) RETURN    ! <modstr>.nml does not exist

    READ(iou, NML=CPL, IOSTAT=fstat)
    CALL read_nml_check(fstat, substr, iou, 'CPL', modstr)
    IF (fstat /= 0) RETURN  ! error while reading namelist

    CALL read_nml_close(substr, iou, modstr)

    status = 0  ! no ERROR

  END SUBROUTINE lgtmix_read_nml_cpl
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
SUBROUTINE str2idx(str, idx, status)

  ! ECHAM5/MESSy
  USE messy_main_grid_def_mem_bi, ONLY: nlev

  IMPLICIT NONE

  ! I/O
  CHARACTER(LEN=*), INTENT(INOUT) :: str
  INTEGER,          INTENT(INOUT) :: idx
  INTEGER,          INTENT(OUT)   :: status

  ! LOCAL
  INTEGER :: iostat

  status = 0

  str = ADJUSTL(TRIM(str))
  !
  IF (str(1:1) == '#') THEN
     SELECT CASE(TRIM(str(2:)))
     CASE('GND')
        idx = nlev
     CASE('TOA')
        idx = 1
     CASE DEFAULT
        READ(str(2:),*,IOSTAT=iostat) idx
        IF (iostat /= 0) THEN
           status = 2  ! ERROR IN READING INTEGER
        END IF
     END SELECT
  ELSE
     status = 1  ! STRING DOES NOT START WITH '#'
  END IF

END SUBROUTINE str2idx
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
SUBROUTINE str2r(str, r, status)

  IMPLICIT NONE

  ! I/O
  CHARACTER(LEN=*), INTENT(INOUT) :: str
  REAL(DP),         INTENT(INOUT) :: r
  INTEGER,          INTENT(OUT)   :: status

  ! LOCAL
  INTEGER :: iostat

  status = 0

  str = ADJUSTL(TRIM(str))
  !
  IF (str(1:1) == '=') THEN
     SELECT CASE(TRIM(str(2:)))
     CASE('MAX')
        r = 1.0
     CASE('MIN')
        r = 0.0
     CASE DEFAULT
        READ(str(2:),*,IOSTAT=iostat) r
        IF (iostat /= 0) THEN
           status = 2  ! ERROR IN READING REAL
        END IF
     END SELECT
  ELSE
     status = 1  ! STRING DOES NOT START WITH '~'
  END IF

END SUBROUTINE str2r
! ----------------------------------------------------------------------

! ************************************************************************
END MODULE messy_lgtmix_e5
! ************************************************************************

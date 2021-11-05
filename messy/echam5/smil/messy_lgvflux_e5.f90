! ************************************************************************
MODULE messy_lgvflux_e5
! ************************************************************************

  ! AUTHOR:
  !  Patrick Joeckel, MPICH, January 2006

  ! ECHAM5/MESSy
  USE messy_main_blather_bi,    ONLY: start_message_bi, end_message_bi
  ! MESSy
  USE messy_main_channel,       ONLY: STRLEN_CHANNEL, STRLEN_OBJECT
  USE messy_lgvflux
  USE messy_main_tools,         ONLY: PTR_1D_ARRAY, PTR_2D_ARRAY

  IMPLICIT NONE
  PRIVATE

  INTRINSIC :: NULL
  
  TYPE IO_VDYN
     CHARACTER(LEN=STRLEN_OBJECT)  :: name = ''     ! new object name
     CHARACTER(LEN=STRLEN_CHANNEL) :: cha  = ''     ! channel with hor. surf.
     CHARACTER(LEN=STRLEN_OBJECT)  :: obj  = ''     ! obj with hor. surf.
     REAL(DP)                      :: time = 0.0_DP ! min. transition age
  END TYPE IO_VDYN

  TYPE IO_VFLX
     CHARACTER(LEN=STRLEN_OBJECT)  :: name    = ''  ! new object name
     CHARACTER(LEN=STRLEN_OBJECT)  :: dyname  = ''  ! name of DYN
     CHARACTER(LEN=STRLEN_CHANNEL) :: lgc     = ''  ! LG-quantity to be 
     CHARACTER(LEN=STRLEN_OBJECT)  :: lgo     = ''  ! ... measured ...
     CHARACTER(LEN=STRLEN_CHANNEL) :: lgcte   = ''  ! corresp. tendency.
     CHARACTER(LEN=STRLEN_OBJECT)  :: lgote   = ''  !
  END TYPE IO_VFLX

  TYPE T_LGVFLUX
     TYPE(IO_VDYN)                        :: iodyn
     !
     LOGICAL                              :: ok  = .FALSE. ! OVERALL SWITCH
     !
     LOGICAL                              :: l_age_only = .FALSE. ! no fluxes
     !
     ! COUPLING
     ! horizontal surface [Pa] (GP, in decomposition)
     REAL(DP), DIMENSION(:,:), POINTER :: srf    => NULL()
     !
     ! DYN ANALYSIS
     ! counter (in decomposition)
     ! ... now
     REAL(DP), DIMENSION(:,:), POINTER :: ncu    => NULL() ! upward   (GP)
     REAL(DP), DIMENSION(:,:), POINTER :: ncd    => NULL() ! downward (GP)
     ! ... transition
     REAL(DP), DIMENSION(:,:), POINTER :: ntu    => NULL() ! upward   (GP)
     REAL(DP), DIMENSION(:,:), POINTER :: ntd    => NULL() ! downward (GP)
     ! 
     ! timer and position information (position at transition)
     REAL(DP), DIMENSION(:),   POINTER :: clock  => NULL() ! stop-watch    (LG)
     REAL(DP), DIMENSION(:),   POINTER :: ilat   => NULL() ! position ...  (LG)
     REAL(DP), DIMENSION(:),   POINTER :: ilon   => NULL() ! information   (LG)
     REAL(DP), DIMENSION(:),   POINTER :: pflag  => NULL() ! position-flag (LG)
     !
     ! FLUX ANALYSIS
     INTEGER                                   :: nflux ! number of fluxes
     LOGICAL,            DIMENSION(:), POINTER :: fxok  ! switch
     TYPE(IO_VFLX),      DIMENSION(:), POINTER :: ioflx
     !
     TYPE(PTR_2D_ARRAY), DIMENSION(:), POINTER :: fxu => NULL()
     TYPE(PTR_2D_ARRAY), DIMENSION(:), POINTER :: fxd => NULL()
     !
     TYPE(PTR_1D_ARRAY), DIMENSION(:), POINTER :: qua   => NULL()
     TYPE(PTR_1D_ARRAY), DIMENSION(:), POINTER :: quate => NULL()
     !
     TYPE(PTR_1D_ARRAY), DIMENSION(:), POINTER :: squa   => NULL()
     TYPE(PTR_1D_ARRAY), DIMENSION(:), POINTER :: squate => NULL()
     !
     ! scaling fact. for unit conversion
     REAL(DP),          DIMENSION(:),  POINTER :: uscale => NULL()
     !
  END TYPE T_LGVFLUX

  ! maximum number of dynamical analyses
  INTEGER,                           PARAMETER :: NMAXDYN   = 100
  ! maximum number of flux analyses per dynamical analysis
  INTEGER,                           PARAMETER :: NMAXVFLUX = 120
  INTEGER,                                SAVE :: NDYN      = 0
  TYPE(IO_VDYN),    DIMENSION(NMAXDYN),   SAVE :: vdyn
  TYPE(IO_VFLX),    DIMENSION(NMAXVFLUX), SAVE :: vflx
  TYPE(T_LGVFLUX),  DIMENSION(NMAXDYN),   SAVE :: xvflux

  LOGICAL :: LRUNSM = .TRUE.  ! RUN THIS SUBMODEL

  ! DEFAULT ZERO TENDENCY
  REAL(DP), DIMENSION(:),   POINTER :: zerote => NULL() ! zero tendency
  ! GRID BOX AREA (global)
  REAL(DP), DIMENSION(:,:), POINTER :: garea  => NULL() ! m^2

  ! POSITION INFORMATION
  REAL(DP), DIMENSION(:),   POINTER :: iplat  => NULL() ! index
  REAL(DP), DIMENSION(:),   POINTER :: iplon  => NULL() ! index
  REAL(DP), DIMENSION(:),   POINTER :: ppress => NULL() ! Pa
  REAL(DP),                 POINTER :: amcell => NULL() ! kg

  ! SUBROUTINES
  PUBLIC :: lgvflux_initialize
  PUBLIC :: lgvflux_init_coupling
  PUBLIC :: lgvflux_global_start
  PUBLIC :: lgvflux_global_end
  PUBLIC :: lgvflux_free_memory
  !PRIVATE :: lgvflux_read_nml_cpl

CONTAINS

  ! ######################################################################
  ! PUBLIC SUBROUTINES
  ! ######################################################################
  
  ! ---------------------------------------------------------------------
  SUBROUTINE lgvflux_initialize

    ! lgvflux MODULE ROUTINE (ECHAM5 INTERFACE)
    !
    ! INITIALIZATION OF GLOBAL VARIABLES FROM NAMELIST
    ! IN PARALLEL ENVIRONMENT
    !
    ! Author: Patrick Joeckel, MPICH, Jan 2006
    
    ! ECHAM5/MESSy
    USE messy_main_mpi_bi,     ONLY: p_parallel_io, p_io, p_bcast, finish
    USE messy_main_channel_bi, ONLY: LG_ATTILA, REPR_UNDEF
    USE messy_main_tools,      ONLY: find_next_free_unit

    IMPLICIT NONE

    INTRINSIC :: TRIM
    
    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'lgvflux_initialize'
    INTEGER                 :: iou    ! I/O unit
    INTEGER                 :: status ! error status
    INTEGER                 :: i, j, k, n

    CALL start_message_bi(modstr, 'INITIALISATION', substr)

    LRUNSM = (LG_ATTILA /= REPR_UNDEF)
    IF (.NOT. LRUNSM) THEN
       IF (p_parallel_io) THEN
          WRITE(*,*) 'LAGRANGIAN REPRESENTATION IS UNDEFINED.'
          WRITE(*,*) 'SUBMODEL CANNOT BE RUN.'
       END IF
       CALL end_message_bi(modstr, 'INITIALISATION', substr)
       RETURN
    END IF
    
    ! INITIALIZE CPL
    IF (p_parallel_io) THEN
       iou = find_next_free_unit(100,200)
       CALL lgvflux_read_nml_cpl(status, iou)
       IF (status /= 0) CALL finish(substr)
    END IF

    ! IO (1st part)
    IF (p_parallel_io) THEN

       ! GET NUMBER OF ENTRIES
       WRITE(*,*) '------------------------------------------------------'
       WRITE(*,*) 'VERTICAL FLUX DIAGNOSTIC:'
       ! COPY DATA
       NDYN = 1
       DO i=1, NMAXDYN

          IF (TRIM(vdyn(i)%name) == '') CYCLE
          xvflux(NDYN)%iodyn%name = TRIM(vdyn(i)%name)
          
          IF (TRIM(vdyn(i)%cha) == '') CYCLE
          xvflux(NDYN)%iodyn%cha = TRIM(vdyn(i)%cha)

          IF (TRIM(vdyn(i)%obj) == '') CYCLE
          xvflux(NDYN)%iodyn%obj = TRIM(vdyn(i)%obj)

          xvflux(NDYN)%iodyn%time = vdyn(i)%time

          WRITE(*,*) ' # ', TRIM(xvflux(NDYN)%iodyn%name),': '
          WRITE(*,*) '    surface  : ', TRIM(xvflux(NDYN)%iodyn%cha)    &
               , '(',TRIM(xvflux(NDYN)%iodyn%obj),')'

          xvflux(NDYN)%l_age_only = (xvflux(NDYN)%iodyn%time < 0.0_dp)

          age1: IF (.NOT. xvflux(NDYN)%l_age_only) THEN

             WRITE(*,*) '    min. time: ', xvflux(NDYN)%iodyn%time,' s'

             ! COUNT ASSOCIATED FLUXES
             xvflux(NDYN)%nflux = 1
             DO j=1, NMAXVFLUX
                
                IF (TRIM(vflx(j)%name) == '') CYCLE
                IF (TRIM(vflx(j)%dyname) /= TRIM(xvflux(NDYN)%iodyn%name)) &
                     CYCLE

                WRITE(*,*) '   -> ',TRIM(vflx(j)%name) 
                WRITE(*,*) '      ',' quantity : ' &
                     , TRIM(vflx(j)%lgc) &
                     , '(',TRIM(vflx(j)%lgo),')'
                
                IF (TRIM(vflx(j)%lgc) /= '') THEN
                   WRITE(*,*) '      ', ' tendency : ' &
                        , TRIM(vflx(j)%lgcte) &
                        , '(',TRIM(vflx(j)%lgote),')'
                END IF
                
                ! NEXT ENTRY
                xvflux(NDYN)%nflux = xvflux(NDYN)%nflux + 1
             END DO
             xvflux(NDYN)%nflux = xvflux(NDYN)%nflux - 1

             WRITE(*,*) ' - ',xvflux(NDYN)%nflux,' associated fluxes reqested'
          
          ELSE

             xvflux(NDYN)%nflux = 0
             WRITE(*,*) '    min. time: <infinite>'
             WRITE(*,*) ' - ',' Flux analysis is undefined. '//&
                  &'Only age of air analysis is possible.'

          END IF age1

          ! NEXT ENTRY
          NDYN = NDYN + 1
          WRITE(*,*) '------------------------------------------------------'

       END DO
       NDYN = NDYN - 1
       
    END IF

    ! BROADCAST RESULTS (1st part)
    CALL p_bcast(NDYN, p_io)

    DO i=1, NDYN
       CALL p_bcast(xvflux(i)%iodyn%name, p_io)
       CALL p_bcast(xvflux(i)%iodyn%cha,  p_io)
       CALL p_bcast(xvflux(i)%iodyn%obj,  p_io)
       CALL p_bcast(xvflux(i)%iodyn%time, p_io)
       CALL p_bcast(xvflux(i)%nflux,      p_io)
       CALL p_bcast(xvflux(i)%l_age_only, p_io)
    END DO

    ! ALLOCATE SPACE
    DO i=1, NDYN
       IF (xvflux(i)%nflux > 0) THEN
          n = xvflux(i)%nflux
          !
          ALLOCATE(xvflux(i)%fxok(n))
          xvflux(i)%fxok(:) = .FALSE.
          !
          ALLOCATE(xvflux(i)%uscale(n))
          xvflux(i)%uscale(:) = 1.0_dp
          !
          ALLOCATE(xvflux(i)%ioflx(n))
          !
          ALLOCATE(xvflux(i)%fxu(n))
          ALLOCATE(xvflux(i)%fxd(n))
          ALLOCATE(xvflux(i)%qua(n))
          ALLOCATE(xvflux(i)%quate(n))
          ALLOCATE(xvflux(i)%squa(n))
          ALLOCATE(xvflux(i)%squate(n))

          DO j=1, n
             NULLIFY(xvflux(i)%fxu(j)%ptr)
             NULLIFY(xvflux(i)%fxd(j)%ptr)
             NULLIFY(xvflux(i)%qua(j)%ptr)
             NULLIFY(xvflux(i)%quate(j)%ptr)
             NULLIFY(xvflux(i)%squa(j)%ptr)
             NULLIFY(xvflux(i)%squate(j)%ptr)
          END DO
       END IF
    END DO

    ! IO (2nd part)
    IF (p_parallel_io) THEN
       DO i=1, NDYN
          
          age2: IF (.NOT. xvflux(i)%l_age_only) THEN

             xvflux(i)%nflux = 1
             DO j=1, NMAXVFLUX
                
                IF (TRIM(vflx(j)%name) == '') CYCLE
                IF (TRIM(vflx(j)%dyname) /= TRIM(xvflux(i)%iodyn%name)) CYCLE

                k = xvflux(i)%nflux

                xvflux(i)%ioflx(k)%name   = TRIM(vflx(j)%name)
                xvflux(i)%ioflx(k)%dyname = TRIM(vflx(j)%dyname)
                
                xvflux(i)%ioflx(k)%lgc   = TRIM(vflx(j)%lgc)
                xvflux(i)%ioflx(k)%lgo   = TRIM(vflx(j)%lgo)
                xvflux(i)%ioflx(k)%lgcte = TRIM(vflx(j)%lgcte)
                xvflux(i)%ioflx(k)%lgote = TRIM(vflx(j)%lgote)
                
                ! NEXT ENTRY
                xvflux(i)%nflux = xvflux(i)%nflux + 1
             END DO
             xvflux(i)%nflux = xvflux(i)%nflux - 1

          END IF age2

       END DO
    END IF

    ! BROADCAST RESULTS (2nd part)
    DO i=1, NDYN
       DO k=1, xvflux(i)%nflux
          CALL p_bcast(xvflux(i)%ioflx(k)%name,   p_io)
          CALL p_bcast(xvflux(i)%ioflx(k)%dyname, p_io)
          CALL p_bcast(xvflux(i)%ioflx(k)%lgc,    p_io)
          CALL p_bcast(xvflux(i)%ioflx(k)%lgcte,  p_io)
          CALL p_bcast(xvflux(i)%ioflx(k)%lgo,    p_io)
          CALL p_bcast(xvflux(i)%ioflx(k)%lgote,  p_io)
       END DO
    END DO

    IF (p_parallel_io) THEN
       WRITE(*,*) ' ---> ',NDYN,' FLUX DIAGNOSTIC(S) INITIALIZED !'
    END IF

    LRUNSM = (NDYN > 0)

    CALL end_message_bi(modstr, 'INITIALISATION', substr)

  END SUBROUTINE lgvflux_initialize
  ! ---------------------------------------------------------------------

  ! ---------------------------------------------------------------------
  SUBROUTINE lgvflux_init_coupling

    ! ECHAM5/MESSy
    USE messy_main_mpi_bi,           ONLY: p_parallel_io, message
    USE messy_main_channel_error_bi, ONLY: channel_halt
    USE messy_main_channel_bi,       ONLY: GP_2D_HORIZONTAL, LG_ATTILA
    USE messy_main_tracer_mem_bi,    ONLY: NCELL

    ! MESSy
    USE messy_main_channel,         ONLY: new_channel, new_channel_object &
                                        , new_attribute, get_attribute    &
                                        , get_channel_object              &
                                        , get_channel_object_info
    USE messy_main_constants_mem,   ONLY: STRLEN_ULONG, M_air

    IMPLICIT NONE

    INTRINSIC :: TRIM

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'lgvflux_init_coupling'
    INTEGER :: status
    LOGICAL :: lfirst
    INTEGER :: i, j
    INTEGER :: reprid
    LOGICAL :: lok_unit
    CHARACTER(LEN=STRLEN_ULONG) :: unit, unitte

    IF (.NOT.LRUNSM) RETURN

    CALL start_message_bi(modstr, 'CHANNEL DEFINITION', substr)

    ! ### (1) COUPLING TO ATTILA (POSITION INFORMATION)
    CALL get_channel_object(status, 'attila', 'IPLAT', p1=iplat)
    CALL channel_halt(substr, status)
    CALL get_channel_object(status, 'attila', 'IPLON', p1=iplon)
    CALL channel_halt(substr, status)
    CALL get_channel_object(status, 'attila', 'PPRESS', p1=ppress)
    CALL channel_halt(substr, status)
    CALL get_channel_object(status, 'attila', 'AMCELL', p0=amcell)
    CALL channel_halt(substr, status)

    ALLOCATE(zerote(NCELL))
    zerote(:) = 0.0_DP

    ! ### (2) COUPLING FOR FLUX CALCULATION
    IF (p_parallel_io) THEN
       WRITE(*,*) 'VERTICAL FLUX DIAGNOSTIC:'
    END IF
    lfirst = .TRUE.

    dyn_loop: DO i=1, NDYN

       ! CHECK AVAILABILITY OF GP CHANNEL/OBJECT (HORIZONTAL SURFACE)
       CALL message('  ','# '//TRIM(xvflux(i)%iodyn%name))
       CALL get_channel_object(status &
            , TRIM(xvflux(i)%iodyn%cha) & 
            , TRIM(xvflux(i)%iodyn%obj) &
            , p2=xvflux(i)%srf)
       IF (status /= 0) THEN
          CALL message('  ',' ... '//TRIM(xvflux(i)%iodyn%cha)//' - '//&
               &TRIM(xvflux(i)%iodyn%obj)//' not found ... skipping' )
          CYCLE
       ELSE
          CALL message('  ',' ... '//TRIM(xvflux(i)%iodyn%cha)//' - '//&
               &TRIM(xvflux(i)%iodyn%obj)//' found ... ' )
          ! CHECK REPRESENTATION
          CALL get_channel_object_info(status &
               , TRIM(xvflux(i)%iodyn%cha), TRIM(xvflux(i)%iodyn%obj) &
               , reprid=reprid)
          IF (reprid /= GP_2D_HORIZONTAL) THEN
             CALL message('   ',' ... wrong representation ... skipping')
             CYCLE
          ELSE
             CALL message('   ',' ... GP_2D_HORIZONTAL representation ... ')
          END IF
          ! CHECK UNIT
          unit   = ''
          CALL get_attribute(status &
               , TRIM(xvflux(i)%iodyn%cha), TRIM(xvflux(i)%iodyn%obj) &
               , 'units', c=unit)
          IF (status /= 0) THEN
             CALL message('   ',' ... unknown unit ... skipping')
             CYCLE
          ELSE
             IF (TRIM(unit) /= 'Pa') THEN
                CALL message('   ',' ... wrong unit ['//TRIM(unit)//&
                     &'] ... skipping')
                CYCLE
             ELSE
                CALL message('   ',' ... unit is [Pa] ... ')
             END IF
          END IF
       END IF

       ! DYNAMIC PART IS OK HERE
       xvflux(i)%ok = .TRUE.

       ! CREATE CHANNEL
       IF (lfirst) THEN
          lfirst = .FALSE.
          CALL new_channel(status, modstr)
          CALL channel_halt(substr, status)
       END IF

       ! CREATE CHANNEL OBJECTS WHICH DEPEND ONLY ON SELECTED
       ! GRIDPOINT HORIZONTAL SURFACE (AND LAGRANGIAN TRAJECTORIES)

       ! CLOCK
       CALL new_channel_object(status, modstr &
            , TRIM(xvflux(i)%iodyn%name)//'_clock', p1=xvflux(i)%clock &
            , reprid=LG_ATTILA, lrestreq=.TRUE.)
       !
!!$       IF (.NOT. xvflux(i)%l_age_only) &
            xvflux(i)%clock(:) = -1.0_DP  ! INITIALIZE CORRECTLY
       !
       CALL channel_halt(substr, status)
       CALL new_attribute(status, modstr &
            , TRIM(xvflux(i)%iodyn%name)//'_clock' &
            , 'long_name', c='time since last transition')
       CALL channel_halt(substr, status)
       !
       CALL new_attribute(status, modstr &
            , TRIM(xvflux(i)%iodyn%name)//'_clock' &
            , 'units', c='s')
       CALL channel_halt(substr, status)
       !
       CALL new_attribute(status, modstr &
            , TRIM(xvflux(i)%iodyn%name)//'_clock' &
            , 'channel', c=TRIM(xvflux(i)%iodyn%cha))
       CALL channel_halt(substr, status)
       !
       CALL new_attribute(status, modstr &
            , TRIM(xvflux(i)%iodyn%name)//'_clock' &
            , 'object', c=TRIM(xvflux(i)%iodyn%obj))
       CALL channel_halt(substr, status)
       !
       CALL new_attribute(status, modstr &
            , TRIM(xvflux(i)%iodyn%name)//'_clock' &
            , 'minimum_age', r=xvflux(i)%iodyn%time)
       CALL channel_halt(substr, status)

       ! POSITION-FLAG (above, below)
       CALL new_channel_object(status, modstr &
            , TRIM(xvflux(i)%iodyn%name)//'_pflag', p1=xvflux(i)%pflag &
            , reprid=LG_ATTILA, lrestreq=.TRUE.)
       CALL channel_halt(substr, status)
       !
       CALL new_attribute(status, modstr &
            , TRIM(xvflux(i)%iodyn%name)//'_pflag' &
            , 'long_name', c='position flag (<0: below, >0: above)')
       CALL channel_halt(substr, status)
       !
       CALL new_attribute(status, modstr &
            , TRIM(xvflux(i)%iodyn%name)//'_pflag' &
            , 'channel', c=TRIM(xvflux(i)%iodyn%cha))
       CALL channel_halt(substr, status)
       !
       CALL new_attribute(status, modstr &
            , TRIM(xvflux(i)%iodyn%name)//'_pflag' &
            , 'object', c=TRIM(xvflux(i)%iodyn%obj))
       CALL channel_halt(substr, status)
       !
       CALL new_attribute(status, modstr &
            , TRIM(xvflux(i)%iodyn%name)//'_pflag' &
            , 'minimum_age', r=xvflux(i)%iodyn%time)
       CALL channel_halt(substr, status)

       IF (xvflux(i)%l_age_only) CYCLE

       ! COUNTER
       ! ---
       CALL new_channel_object(status, modstr &
            , TRIM(xvflux(i)%iodyn%name)//'_ncu', p2=xvflux(i)%ncu &
            , reprid=GP_2D_HORIZONTAL, lrestreq=.TRUE.)
       CALL channel_halt(substr, status)
       !
       CALL new_attribute(status, modstr &
            , TRIM(xvflux(i)%iodyn%name)//'_ncu' &
            , 'long_name', c='upward flux event counts (current position)')
       CALL channel_halt(substr, status)
       !
       CALL new_attribute(status, modstr &
            , TRIM(xvflux(i)%iodyn%name)//'_ncu' &
            , 'channel', c=TRIM(xvflux(i)%iodyn%cha))
       CALL channel_halt(substr, status)
       !
       CALL new_attribute(status, modstr &
            , TRIM(xvflux(i)%iodyn%name)//'_ncu' &
            , 'object', c=TRIM(xvflux(i)%iodyn%obj))
       CALL channel_halt(substr, status)
       !
       CALL new_attribute(status, modstr &
            , TRIM(xvflux(i)%iodyn%name)//'_ncu' &
            , 'minimum_age', r=xvflux(i)%iodyn%time)
       CALL channel_halt(substr, status)

       ! ---

       CALL new_channel_object(status, modstr &
            , TRIM(xvflux(i)%iodyn%name)//'_ncd', p2=xvflux(i)%ncd &
            , reprid=GP_2D_HORIZONTAL, lrestreq=.TRUE.)
       !
       CALL channel_halt(substr, status)
       CALL new_attribute(status, modstr &
            , TRIM(xvflux(i)%iodyn%name)//'_ncd' &
            , 'long_name', c='downward flux event counts (current position)')
       CALL channel_halt(substr, status)
       !
       CALL new_attribute(status, modstr &
            , TRIM(xvflux(i)%iodyn%name)//'_ncd' &
            , 'channel', c=TRIM(xvflux(i)%iodyn%cha))
       CALL channel_halt(substr, status)
       !
       CALL new_attribute(status, modstr &
            , TRIM(xvflux(i)%iodyn%name)//'_ncd' &
            , 'object', c=TRIM(xvflux(i)%iodyn%obj))
       CALL channel_halt(substr, status)
       !
       CALL new_attribute(status, modstr &
            , TRIM(xvflux(i)%iodyn%name)//'_ncd' &
            , 'minimum_age', r=xvflux(i)%iodyn%time)
       CALL channel_halt(substr, status)

       ! ---

       CALL new_channel_object(status, modstr &
            , TRIM(xvflux(i)%iodyn%name)//'_ntu', p2=xvflux(i)%ntu &
            , reprid=GP_2D_HORIZONTAL, lrestreq=.TRUE.)
       CALL channel_halt(substr, status)
       !
       CALL new_attribute(status, modstr &
            , TRIM(xvflux(i)%iodyn%name)//'_ntu' &
            , 'long_name', c='upward flux event counts (transition position)')
       CALL channel_halt(substr, status)
       !
       CALL new_attribute(status, modstr &
            , TRIM(xvflux(i)%iodyn%name)//'_ntu' &
            , 'channel', c=TRIM(xvflux(i)%iodyn%cha))
       CALL channel_halt(substr, status)
       !
       CALL new_attribute(status, modstr &
            , TRIM(xvflux(i)%iodyn%name)//'_ntu' &
            , 'object', c=TRIM(xvflux(i)%iodyn%obj))
       CALL channel_halt(substr, status)
       !
       CALL new_attribute(status, modstr &
            , TRIM(xvflux(i)%iodyn%name)//'_ntu' &
            , 'minimum_age', r=xvflux(i)%iodyn%time)
       CALL channel_halt(substr, status)

       ! ---

       CALL new_channel_object(status, modstr &
            , TRIM(xvflux(i)%iodyn%name)//'_ntd', p2=xvflux(i)%ntd &
            , reprid=GP_2D_HORIZONTAL, lrestreq=.TRUE.)
       CALL channel_halt(substr, status)
       !
       CALL new_attribute(status, modstr &
            , TRIM(xvflux(i)%iodyn%name)//'_ntd' &
            , 'long_name'  &
            , c='downward flux event counts (transition position)')
       CALL channel_halt(substr, status)
       !
       CALL new_attribute(status, modstr &
            , TRIM(xvflux(i)%iodyn%name)//'_ntd' &
            , 'channel', c=TRIM(xvflux(i)%iodyn%cha))
       CALL channel_halt(substr, status)
       !
       CALL new_attribute(status, modstr &
            , TRIM(xvflux(i)%iodyn%name)//'_ntd' &
            , 'object', c=TRIM(xvflux(i)%iodyn%obj))
       CALL channel_halt(substr, status)
       !
       CALL new_attribute(status, modstr &
            , TRIM(xvflux(i)%iodyn%name)//'_ntd' &
            , 'minimum_age', r=xvflux(i)%iodyn%time)
       CALL channel_halt(substr, status)

       ! ---

       ! latitude (index) of transition
       CALL new_channel_object(status, modstr &
            , TRIM(xvflux(i)%iodyn%name)//'_ilat', p1=xvflux(i)%ilat &
            , reprid=LG_ATTILA, lrestreq=.TRUE.)
       CALL channel_halt(substr, status)
       !
       CALL new_attribute(status, modstr &
            , TRIM(xvflux(i)%iodyn%name)//'_ilat' &
            , 'long_name', c='latitude index (global grid) of transition')
       CALL channel_halt(substr, status)
       !
       CALL new_attribute(status, modstr &
            , TRIM(xvflux(i)%iodyn%name)//'_ilat' &
            , 'channel', c=TRIM(xvflux(i)%iodyn%cha))
       CALL channel_halt(substr, status)
       !
       CALL new_attribute(status, modstr &
            , TRIM(xvflux(i)%iodyn%name)//'_ilat' &
            , 'object', c=TRIM(xvflux(i)%iodyn%obj))
       CALL channel_halt(substr, status)
       !
       CALL new_attribute(status, modstr &
            , TRIM(xvflux(i)%iodyn%name)//'_ilat' &
            , 'minimum_age', r=xvflux(i)%iodyn%time)
       CALL channel_halt(substr, status)

       ! longitude (index) of transition
       CALL new_channel_object(status, modstr &
            , TRIM(xvflux(i)%iodyn%name)//'_ilon', p1=xvflux(i)%ilon &
            , reprid=LG_ATTILA, lrestreq=.TRUE.)
       CALL channel_halt(substr, status)
       !
       CALL new_attribute(status, modstr &
            , TRIM(xvflux(i)%iodyn%name)//'_ilon' &
            , 'long_name', c='longitude index (global grid) of transition')
       CALL channel_halt(substr, status)
       !
       CALL new_attribute(status, modstr &
            , TRIM(xvflux(i)%iodyn%name)//'_ilon' &
            , 'channel', c=TRIM(xvflux(i)%iodyn%cha))
       CALL channel_halt(substr, status)
       !
       CALL new_attribute(status, modstr &
            , TRIM(xvflux(i)%iodyn%name)//'_ilon' &
            , 'object', c=TRIM(xvflux(i)%iodyn%obj))
       CALL channel_halt(substr, status)
       !
       CALL new_attribute(status, modstr &
            , TRIM(xvflux(i)%iodyn%name)//'_ilon' &
            , 'minimum_age', r=xvflux(i)%iodyn%time)
       CALL channel_halt(substr, status)

       ! SETTINGS FOR ASSOCIATED FLUXES
       flux_loop: DO j=1, xvflux(i)%nflux
    
          ! COUPLING

          ! INIT
          unit   = ''
          unitte = ''
          xvflux(i)%quate(j)%ptr => zerote(:)

          ! CHECK AVAILABILITY OF LG CHANNEL/OBJECT
          CALL message('   ','  >  '//TRIM(xvflux(i)%ioflx(j)%name))
          CALL get_channel_object(status &
               , TRIM(xvflux(i)%ioflx(j)%lgc) & 
               , TRIM(xvflux(i)%ioflx(j)%lgo) &
               , p1=xvflux(i)%qua(j)%ptr)
          IF (status /= 0) THEN
             CALL message('   ','     ... '//&
                  &TRIM(xvflux(i)%ioflx(j)%lgc)//' - '//&
                  &TRIM(xvflux(i)%ioflx(j)%lgo)//' not found ... skipping' )
             CYCLE
          ELSE
             CALL message('   ','     ... '//&
                  &TRIM(xvflux(i)%ioflx(j)%lgc)//' - '//&
                  &TRIM(xvflux(i)%ioflx(j)%lgo)//' found ... ' )
             ! CHECK REPRESENTATION
             CALL get_channel_object_info(status &
                  , TRIM(xvflux(i)%ioflx(j)%lgc) &
                  , TRIM(xvflux(i)%ioflx(j)%lgo) &
                  , reprid=reprid)
             IF (reprid /= LG_ATTILA) THEN
                CALL message('   ',&
                     '     ... wrong representation ... skipping')
                CYCLE
             ELSE
                CALL message('   ','     ... LG_ATTILA representation ... ')
             END IF
             ! CHECK UNIT
             CALL get_attribute(status &
                  , TRIM(xvflux(i)%ioflx(j)%lgc) &
                  , TRIM(xvflux(i)%ioflx(j)%lgo) &
                  , 'units', c=unit)
             IF (status == 0) THEN
                CALL message('   ','     ... unit is ['//TRIM(unit)//'] ... ')
                SELECT CASE(TRIM(unit))
                CASE('mol/mol','mol mol-1')
                   CALL get_attribute(status &
                        , TRIM(xvflux(i)%ioflx(j)%lgc) &
                        , TRIM(xvflux(i)%ioflx(j)%lgo) &
                        , 'molarmass', r=xvflux(i)%uscale(j))
                   IF (status == 0) THEN
                      IF (p_parallel_io) THEN
                         WRITE(*,*) '          ... molar mass is ' &
                              ,xvflux(i)%uscale(j)
                      END IF
                      xvflux(i)%uscale(j) = xvflux(i)%uscale(j)/M_air
                   ELSE
                      xvflux(i)%uscale(j) = 1.0_dp
                   END IF
                   lok_unit = .TRUE.
                CASE('kg/kg','kg kg-1')
                   xvflux(i)%uscale(j) = 1.0_DP
                   lok_unit = .TRUE.
                CASE DEFAULT
                   xvflux(i)%uscale(j) = 1.0_DP
                   lok_unit = .FALSE.
                END SELECT
             ELSE
                CALL message('   ','     ... no unit ... ')
                xvflux(i)%uscale(j) = 1.0_DP
                lok_unit = .FALSE.
             END IF
          END IF
          
          ! CHECK AVAILABILITY OF LG CHANNEL/OBJECT (TENDENCY)
          tendency: IF (TRIM(xvflux(i)%ioflx(j)%lgcte) /= '') THEN
             CALL get_channel_object(status, TRIM(xvflux(i)%ioflx(j)%lgcte) & 
                  , TRIM(xvflux(i)%ioflx(j)%lgote), p1=xvflux(i)%quate(j)%ptr)
             IF (status /= 0) THEN
                CALL message('    ',&
                     '  +  '//TRIM(xvflux(i)%ioflx(j)%lgcte)//' - '//&
                     &TRIM(xvflux(i)%ioflx(j)%lgote)//&
                     &' not found ... ignoring' )
                xvflux(i)%quate(j)%ptr => zerote(:)
             ELSE
                CALL message('    ','  +  '//TRIM(xvflux(i)%ioflx(j)%lgcte)//&
                     &' - '//TRIM(xvflux(i)%ioflx(j)%lgote)//' found ... ' )
                ! CHECK REPRESENTATION
                CALL get_channel_object_info(status &
                     , TRIM(xvflux(i)%ioflx(j)%lgcte) &
                     , TRIM(xvflux(i)%ioflx(j)%lgote) &
                     , reprid=reprid)
                IF (reprid /= LG_ATTILA) THEN
                   CALL message('   ',&
                        '     ... wrong representation ... ignoring')
                   xvflux(i)%quate(j)%ptr => zerote(:)
                ELSE
                   CALL message('   ','     ... LG_ATTILA representation ... ')
                END IF
                ! CHECK UNIT
                CALL get_attribute(status &
                     , TRIM(xvflux(i)%ioflx(j)%lgcte) &
                     , TRIM(xvflux(i)%ioflx(j)%lgote) &
                     , 'units', c=unitte)
                IF (status == 0) THEN
                   CALL message('   ' &
                        , '     ... unit is ['//TRIM(unitte)//'] ... ')
                ELSE
                   CALL message('   ','     ... no unit ... ')
                END IF
             END IF
          END IF tendency
          
          ! NOW EVERYTHING OK FOR THIS FLUX
          xvflux(i)%fxok(j) = .TRUE.
          
          ! CREATE CHANNEL OBJECTS WHICH DEPEND ON SELECTED
          ! LAGRANGIAN QUANTITY
          ! UPWARD FLUX
          CALL new_channel_object(status, modstr &
               , TRIM(xvflux(i)%ioflx(j)%name)//'_fxu' &
               , p2=xvflux(i)%fxu(j)%ptr &
               , reprid=GP_2D_HORIZONTAL, lrestreq=.TRUE.)
          CALL channel_halt(substr, status)
          ! ... long_name
          CALL new_attribute(status, modstr &
               , TRIM(xvflux(i)%ioflx(j)%name)//'_fxu' &
               , 'long_name', c='upward flux of '//&
               &TRIM(xvflux(i)%ioflx(j)%lgc)//'@'//&
               &TRIM(xvflux(i)%ioflx(j)%lgo)//&
               &' through '//&
               &TRIM(xvflux(i)%iodyn%name) &
               )
          CALL channel_halt(substr, status)
          !
          ! ... channel of flux quantity
          CALL new_attribute(status, modstr &
               , TRIM(xvflux(i)%ioflx(j)%name)//'_fxu' &
               , 'channel', c=TRIM(xvflux(i)%ioflx(j)%lgc))
          CALL channel_halt(substr, status)
          !
          ! ... object of flux quantity
          CALL new_attribute(status, modstr &
               , TRIM(xvflux(i)%ioflx(j)%name)//'_fxu' &
               , 'object_lg', c=TRIM(xvflux(i)%ioflx(j)%lgo))
          CALL channel_halt(substr, status)
          !
          ! ... channel/object of tendency
          IF ( (TRIM(xvflux(i)%ioflx(j)%lgcte) /= '') .AND. &
               (TRIM(xvflux(i)%ioflx(j)%lgote) /= '') ) THEN
             CALL new_attribute(status, modstr &
                  , TRIM(xvflux(i)%ioflx(j)%name)//'_fxu' &
                  , 'channel_tendency', c=TRIM(xvflux(i)%ioflx(j)%lgcte))
             CALL channel_halt(substr, status)
             !
             CALL new_attribute(status, modstr &
                  , TRIM(xvflux(i)%ioflx(j)%name)//'_fxu' &
                  , 'object_tendency', c=TRIM(xvflux(i)%ioflx(j)%lgote))
             CALL channel_halt(substr, status)
          END IF
          !
          ! ... channel of surface
          CALL new_attribute(status, modstr &
               , TRIM(xvflux(i)%ioflx(j)%name)//'_fxu' &
               , 'surface', c=TRIM(xvflux(i)%iodyn%name))
          CALL channel_halt(substr, status)
          
          ! DONWARD FLUX
          CALL new_channel_object(status, modstr &
               , TRIM(xvflux(i)%ioflx(j)%name)//'_fxd' &
               , p2=xvflux(i)%fxd(j)%ptr &
               , reprid=GP_2D_HORIZONTAL, lrestreq=.TRUE.)
          CALL channel_halt(substr, status)
          ! ... long_name
          CALL new_attribute(status, modstr &
               , TRIM(xvflux(i)%ioflx(j)%name)//'_fxd' &
               , 'long_name', c='downward flux of '//&
               &TRIM(xvflux(i)%ioflx(j)%lgc)//'@'//&
               &TRIM(xvflux(i)%ioflx(j)%lgo)//&
               &' through '//&
               &TRIM(xvflux(i)%iodyn%name) &
               )
          CALL channel_halt(substr, status)
          !
          ! ... channel of flux quantity
          CALL new_attribute(status, modstr &
               , TRIM(xvflux(i)%ioflx(j)%name)//'_fxd' &
               , 'channel', c=TRIM(xvflux(i)%ioflx(j)%lgc))
          CALL channel_halt(substr, status)
          !
          ! ... object of flux quantity
          CALL new_attribute(status, modstr &
               , TRIM(xvflux(i)%ioflx(j)%name)//'_fxd' &
               , 'object_lg', c=TRIM(xvflux(i)%ioflx(j)%lgo))
          CALL channel_halt(substr, status)
          !
          ! ... channel/object of tendency
          IF ( (TRIM(xvflux(i)%ioflx(j)%lgcte) /= '') .AND. &
               (TRIM(xvflux(i)%ioflx(j)%lgote) /= '') ) THEN
             CALL new_attribute(status, modstr &
                  , TRIM(xvflux(i)%ioflx(j)%name)//'_fxd' &
                  , 'channel_tendency', c=TRIM(xvflux(i)%ioflx(j)%lgcte))
             CALL channel_halt(substr, status)
             !
             CALL new_attribute(status, modstr &
                  , TRIM(xvflux(i)%ioflx(j)%name)//'_fxd' &
                  , 'object_tendency', c=TRIM(xvflux(i)%ioflx(j)%lgote))
             CALL channel_halt(substr, status)
          END IF
          !
          ! ... channel of surface
          CALL new_attribute(status, modstr &
               , TRIM(xvflux(i)%ioflx(j)%name)//'_fxd' &
               , 'surface', c=TRIM(xvflux(i)%iodyn%name))
          CALL channel_halt(substr, status)
          
          
          IF (lok_unit) THEN
             CALL new_attribute(status, modstr &
                  , TRIM(xvflux(i)%ioflx(j)%name)//'_fxu' &
                  , 'units', c='kg/(m^2 s)')
             CALL channel_halt(substr, status)
             !
             CALL new_attribute(status, modstr &
                  , TRIM(xvflux(i)%ioflx(j)%name)//'_fxd' &
                  , 'units', c='kg/(m^2 s)')
             CALL channel_halt(substr, status)
             !
          ELSE
             CALL new_attribute(status, modstr &
                  , TRIM(xvflux(i)%ioflx(j)%name)//'_fxu' &
                  , 'units', c='<unknown>')
             CALL channel_halt(substr, status)
             !
             CALL new_attribute(status, modstr &
                  , TRIM(xvflux(i)%ioflx(j)%name)//'_fxd' &
                  , 'units', c='<unknown>')
             CALL channel_halt(substr, status)
          END IF
          
          ! SAVED (AT TRANSITION) QUANTITY
          CALL new_channel_object(status, modstr &
               , TRIM(xvflux(i)%ioflx(j)%name)//'_s' &
               , p1=xvflux(i)%squa(j)%ptr &
               , reprid=LG_ATTILA, lrestreq=.TRUE.)
          CALL channel_halt(substr, status)
          CALL new_attribute(status, modstr &
               , TRIM(xvflux(i)%ioflx(j)%name)//'_s' &
               , 'long_name', c=TRIM(xvflux(i)%ioflx(j)%lgc)//&
               &'@'//TRIM(xvflux(i)%ioflx(j)%lgo))
          CALL channel_halt(substr, status)
          CALL new_attribute(status, modstr &
               , TRIM(xvflux(i)%ioflx(j)%name)//'_s' &
               , 'units', c=TRIM(unit))
          CALL channel_halt(substr, status)
          
          ! SAVED (AT TRANSITION) TENDENCY
          CALL new_channel_object(status, modstr &
               , TRIM(xvflux(i)%ioflx(j)%name)//'_ste' &
               , p1=xvflux(i)%squate(j)%ptr &
               , reprid=LG_ATTILA, lrestreq=.TRUE.)
          CALL channel_halt(substr, status)
          CALL new_attribute(status, modstr &
               , TRIM(xvflux(i)%ioflx(j)%name)//'_ste' &
               , 'long_name', c='tendency: '//TRIM(xvflux(i)%ioflx(j)%lgcte)//&
               &'@'//TRIM(xvflux(i)%ioflx(j)%lgote))
          CALL channel_halt(substr, status)
          CALL new_attribute(status, modstr &
               , TRIM(xvflux(i)%ioflx(j)%name)//'_ste' &
               , 'units', c=TRIM(unitte))
          CALL channel_halt(substr, status)
          
       END DO flux_loop
       
    END DO dyn_loop
    
    CALL end_message_bi(modstr, 'CHANNEL DEFINITION', substr)
    
  END SUBROUTINE lgvflux_init_coupling
  ! ---------------------------------------------------------------------

  ! ---------------------------------------------------------------------
  SUBROUTINE lgvflux_global_start

    ! ECHAM5/MESSy
    USE messy_main_grid_def_bi,   ONLY: gboxarea_2d
    USE messy_main_transform_bi,  ONLY: trp_gpdc_gpgl

    IMPLICIT NONE

    ! LOCAL
    LOGICAL, SAVE  :: lfirst = .TRUE.

    IF (.NOT.LRUNSM) RETURN

    IF (lfirst) THEN
       CALL trp_gpdc_gpgl(1, gboxarea_2d, garea)
       lfirst = .FALSE.
    END IF

  END SUBROUTINE lgvflux_global_start
  ! ---------------------------------------------------------------------

  ! ---------------------------------------------------------------------
  SUBROUTINE lgvflux_global_end

    USE messy_main_tracer_mem_bi, ONLY: NCELL
    USE messy_main_transform_bi,  ONLY: trp_gpdc_gpgl, M_SUM
    USE messy_main_timer,         ONLY: delta_time, time_step_len
    USE messy_main_grid_def_mem_bi, ONLY:  nlon, ngl

    IMPLICIT NONE

    INTRINSIC :: ABS, NINT

    ! LOCAL
    !CHARACTER(LEN=*), PARAMETER :: substr = 'lgvflux_global_end'
    ! horizontal surface [Pa] (GP, global)
    REAL(DP), DIMENSION(:,:), POINTER :: gsrf   => NULL()
    ! counter (GP, global)
    REAL(DP), DIMENSION(:,:), POINTER :: gncu   => NULL()
    REAL(DP), DIMENSION(:,:), POINTER :: gncd   => NULL()
    REAL(DP), DIMENSION(:,:), POINTER :: gntu   => NULL()
    REAL(DP), DIMENSION(:,:), POINTER :: gntd   => NULL()
    ! fluxes (GP, global)
    REAL(DP), DIMENSION(:,:,:), POINTER :: gfxu => NULL()
    REAL(DP), DIMENSION(:,:,:), POINTER :: gfxd => NULL()
    REAL(DP), DIMENSION(:,:),   POINTER :: p2   => NULL()
    !
    INTEGER  :: i, j, jc, n
    LOGICAL  :: lbelow, labove
    LOGICAL  :: lup, ldown
    LOGICAL  :: lclock
    INTEGER  :: glat_t, glon_t
    INTEGER  :: glat_c, glon_c
    REAL(DP), DIMENSION(:), POINTER :: scalf => NULL()

    IF (.NOT.LRUNSM) RETURN

    ! DYN
    ALLOCATE(gncu(nlon, ngl))
    ALLOCATE(gncd(nlon, ngl))
    ALLOCATE(gntu(nlon, ngl))
    ALLOCATE(gntd(nlon, ngl))

    dyn_loop: DO i=1, NDYN

       IF (.NOT. xvflux(i)%ok) CYCLE

       ! globalize surface
       CALL trp_gpdc_gpgl(1, xvflux(i)%srf, gsrf)

       ! INITIALIZE COUNTERS
       ! (summation per time step only over number of cells)
       ! -> time average via channel output
       gncu(:,:) = 0.0_DP
       gncd(:,:) = 0.0_DP
       gntu(:,:) = 0.0_DP
       gntd(:,:) = 0.0_DP

       n = xvflux(i)%nflux

       ! INITIALIZE FLUXES
       ! (summation per time step only over number of cells)
       ! -> time average via channel output
       ALLOCATE(gfxu(nlon, ngl, n))
       gfxu(:,:,:) = 0.0_DP
       ALLOCATE(gfxd(nlon, ngl, n))
       gfxd(:,:,:) = 0.0_DP
       
       ! INITIALIZE SCALING: *(kg/s)
       ALLOCATE(scalf(n))
       IF (n > 0) THEN
          scalf(:) = xvflux(i)%uscale(:) * amcell / delta_time
       END IF

       cell_loop: DO jc=1, NCELL

          ! RESET FROM LAST TIME STEP
          IF (ABS(xvflux(i)%pflag(jc)) > 1.0_dp) THEN
             IF (xvflux(i)%l_age_only) THEN
                xvflux(i)%clock(jc)  =  0.0_dp  ! reset clock
             ELSE
                xvflux(i)%clock(jc)  = -1.0_dp  ! reset clock
                xvflux(i)%ilon(jc)   =  0.0_dp  ! reset position
                xvflux(i)%ilat(jc)   =  0.0_dp  ! reset position
             END IF
             DO j=1, n
                IF (.NOT. xvflux(i)%fxok(j) ) CYCLE
                xvflux(i)%squa(j)%ptr(jc)   =  0.0_dp  ! reset quantity
                xvflux(i)%squate(j)%ptr(jc) =  0.0_dp  ! reset tendency
             END DO
          END IF

          ! check position
          glon_c = NINT(iplon(jc))
          glat_c = NINT(iplat(jc))

          labove = ( ppress(jc) <  gsrf( glon_c, glat_c ) )
          lbelow = ( ppress(jc) >= gsrf( glon_c, glat_c ) )
          
          ! check transition through horizontal surface
          lup   = labove .AND. (xvflux(i)%pflag(jc) < 0.0_dp)
          ldown = lbelow .AND. (xvflux(i)%pflag(jc) > 0.0_dp)
          
          ! set new position flag
          IF (lbelow) xvflux(i)%pflag(jc) = -1.0_dp
          IF (labove) xvflux(i)%pflag(jc) = +1.0_dp

          ! check clock
          lclock = (xvflux(i)%clock(jc) >= 0.0_dp) ! clock is running

          transition: IF (lup .OR. ldown) THEN
             ! TRANSITION
             IF (lclock) THEN
                ! CLOCK IS ALREADY RUNNING -> 2nd TRANSITION -> RESET
                IF (xvflux(i)%l_age_only) THEN
                   xvflux(i)%clock(jc)  =  0.0_dp  ! reset clock
                ELSE
                   xvflux(i)%clock(jc)  = -1.0_dp  ! reset clock
                   xvflux(i)%ilon(jc)   =  0.0_dp  ! reset position
                   xvflux(i)%ilat(jc)   =  0.0_dp  ! reset position
                END IF
                DO j=1, n
                   IF (.NOT. xvflux(i)%fxok(j) ) CYCLE
                   xvflux(i)%squa(j)%ptr(jc)   =  0.0_dp  ! reset quantity
                   xvflux(i)%squate(j)%ptr(jc) =  0.0_dp  ! reset tendency
                END DO
             ELSE
                ! CLOCK IS NOT YET RUNNING -> NEW TRANSITION -> START CLOCK
                xvflux(i)%clock(jc)  = 0.0_dp    ! start clock
                IF (.NOT. xvflux(i)%l_age_only) THEN
                   xvflux(i)%ilon(jc)   = iplon(jc) ! save current position
                   xvflux(i)%ilat(jc)   = iplat(jc) ! save current position
                END IF
                DO j=1, n
                   IF (.NOT. xvflux(i)%fxok(j) ) CYCLE
                   ! save current quantity and tendency
                   xvflux(i)%squa(j)%ptr(jc)   = xvflux(i)%qua(j)%ptr(jc)
                   xvflux(i)%squate(j)%ptr(jc) = xvflux(i)%quate(j)%ptr(jc)
                END DO
             END IF
          ELSE
             ! NO TRANSITION
             IF (lclock) THEN
                ! CLOCK IS ALREADY RUNNING -> ADD TIME STEP
                xvflux(i)%clock(jc) = xvflux(i)%clock(jc) + delta_time
             !ELSE
             !   ! CLOCK IS NOT RUNNING -> DO NOTHING
             END IF
          END IF transition

          IF (xvflux(i)%l_age_only) CYCLE

          ! CHECK TIME
          ! check residence time; count flux; reset clock and position
          time: IF (xvflux(i)%clock(jc) >= xvflux(i)%iodyn%time) THEN
             ! get saved position of transition
             glon_t = NINT(xvflux(i)%ilon(jc))
             glat_t = NINT(xvflux(i)%ilat(jc))
             ! count flux
             direction: IF (xvflux(i)%pflag(jc) < 0.0_dp) THEN
                ! now below: downward flux
                ! count
                ! ... transition position
                gntd(glon_t, glat_t) = gntd(glon_t, glat_t) + 1.0_dp
                ! ... current position
                gncd(glon_c, glat_c) = gncd(glon_c, glat_c) + 1.0_dp
                ! flux ! kg/(m^2 s)
                DO j=1, n
                   IF (.NOT. xvflux(i)%fxok(j) ) CYCLE
                   gfxd(glon_t, glat_t, j) =  &
                        gfxd(glon_t, glat_t, j) &
                        + (xvflux(i)%squa(j)%ptr(jc) + &
                        xvflux(i)%squate(j)%ptr(jc)*time_step_len) &
                        * (scalf(j) / garea(glon_t, glat_t))
                END DO
             ELSE
                ! now above: upward flux
                ! count
                ! ... transition position
                gntu(glon_t, glat_t) = gntu(glon_t, glat_t) + 1.0_dp
                ! ... current position
                gncu(glon_c, glat_c) = gncu(glon_c, glat_c) + 1.0_dp
                ! flux ! kg/(m^2 s)
                DO j=1, n
                   IF (.NOT. xvflux(i)%fxok(j) ) CYCLE
                   gfxu(glon_t, glat_t, j) =  &
                        gfxu(glon_t, glat_t, j) &
                        + (xvflux(i)%squa(j)%ptr(jc) + &
                        xvflux(i)%squate(j)%ptr(jc)*time_step_len) &
                        * (scalf(j) / garea(glon_t, glat_t))
                END DO
             END IF direction
             ! TRIGGER RESET IN NEXT TIME STEP
             ! SET TO -2.0 (below) OR +2.0 (above)
             xvflux(i)%pflag(jc) = xvflux(i)%pflag(jc)*(2.0_dp)
          END IF time
          
       END DO cell_loop

       age: IF (.NOT. xvflux(i)%l_age_only) THEN

          ! decompose fluxes
          DO j=1, n
             IF (.NOT. xvflux(i)%fxok(j) ) CYCLE
             p2 => gfxu(:,:,j)
             CALL trp_gpdc_gpgl(-1, xvflux(i)%fxu(j)%ptr, p2, M_SUM)
             p2 => gfxd(:,:,j)
             CALL trp_gpdc_gpgl(-1, xvflux(i)%fxd(j)%ptr, p2, M_SUM)
          END DO
          
          CALL trp_gpdc_gpgl(-1, xvflux(i)%ncu, gncu, M_SUM)
          CALL trp_gpdc_gpgl(-1, xvflux(i)%ncd, gncd, M_SUM)
          CALL trp_gpdc_gpgl(-1, xvflux(i)%ntu, gntu, M_SUM)
          CALL trp_gpdc_gpgl(-1, xvflux(i)%ntd, gntd, M_SUM)
          
       END IF age

       ! CLEAN UP
       DEALLOCATE(gsrf)
       NULLIFY(gsrf)
       DEALLOCATE(gfxu)
       NULLIFY(gfxu)
       DEALLOCATE(gfxd)
       NULLIFY(gfxd)
       DEALLOCATE(scalf)
       NULLIFY(scalf)

    END DO dyn_loop

    ! CLEAN UP
    DEALLOCATE(gncu)
    NULLIFY(gncu)
    DEALLOCATE(gncd)
    NULLIFY(gncd)
    DEALLOCATE(gntu)
    NULLIFY(gntu)
    DEALLOCATE(gntd)
    NULLIFY(gntd)

  END SUBROUTINE lgvflux_global_end
  ! ---------------------------------------------------------------------

  ! ---------------------------------------------------------------------
  SUBROUTINE lgvflux_free_memory

    IMPLICIT NONE

    INTRINSIC :: ASSOCIATED

    ! LOCAL
    INTEGER :: i

    IF (.NOT.LRUNSM) RETURN

    IF (ASSOCIATED(zerote)) THEN
       DEALLOCATE(zerote)
       NULLIFY(zerote)
    END IF

    IF (ASSOCIATED(garea)) THEN
       DEALLOCATE(garea)
       NULLIFY(garea)
    END IF

    DO i=1, NDYN
       IF (xvflux(i)%nflux > 0) THEN
          DEALLOCATE(xvflux(i)%fxok)
          DEALLOCATE(xvflux(i)%ioflx)
          NULLIFY(xvflux(i)%ioflx)
          DEALLOCATE(xvflux(i)%fxu)
          NULLIFY(xvflux(i)%fxu)
          DEALLOCATE(xvflux(i)%fxd)
          NULLIFY(xvflux(i)%fxd)
          DEALLOCATE(xvflux(i)%qua)
          NULLIFY(xvflux(i)%qua)
          DEALLOCATE(xvflux(i)%quate)
          NULLIFY(xvflux(i)%quate)
          DEALLOCATE(xvflux(i)%squa)
          NULLIFY(xvflux(i)%squa)
          DEALLOCATE(xvflux(i)%squate)
          NULLIFY(xvflux(i)%squate)
       END IF
    END DO
    
  END SUBROUTINE lgvflux_free_memory
  ! ---------------------------------------------------------------------

  ! ######################################################################
  ! PRIVATE SUBROUTINES
  ! ######################################################################

  ! ---------------------------------------------------------------------
  SUBROUTINE lgvflux_read_nml_cpl(status, iou)

    ! lgvflux MODULE ROUTINE (ECHAM5 INTERFACE, PRIVATE)
    !
    ! read namelist for 'coupling' to ECHAM5
    !
    ! Author: Patrick Joeckel, MPICH, Jan 2006

    ! MESSy
    USE messy_main_tools, ONLY: read_nml_open, read_nml_check, read_nml_close

    IMPLICIT NONE

    ! I/O
    INTEGER, INTENT(OUT) :: status     ! error status
    INTEGER, INTENT(IN)  :: iou        ! I/O unit

    ! (LOCAL) NAMELIST VARIABLES
    CHARACTER(LEN=*), PARAMETER :: substr = 'lgvflux_read_nml_cpl'

    NAMELIST /CPL/ VDYN, VFLX

    ! LOCAL
    LOGICAL          :: lex      ! file exists ?
    INTEGER          :: fstat    ! file status

    status = 1

    CALL read_nml_open(lex, substr, iou, 'CPL', modstr)
    IF (.not.lex) RETURN    ! <modstr>.nml does not exist

    READ(iou, NML=CPL, IOSTAT=fstat)
    CALL read_nml_check(fstat, substr, iou, 'CPL', modstr)
    IF (fstat /= 0) RETURN  ! error while reading namelist

    CALL read_nml_close(substr, iou, modstr)

    status = 0  ! no ERROR

  END SUBROUTINE lgvflux_read_nml_cpl
  ! ---------------------------------------------------------------------

! ************************************************************************
END MODULE messy_lgvflux_e5
! ************************************************************************

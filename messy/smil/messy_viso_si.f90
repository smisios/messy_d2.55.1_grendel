#include "messy_main_ppd_bi.inc"

! **************************************************************************
MODULE messy_viso_si
! **************************************************************************

  ! MODULE FOR VALUES ON (HORIZONTAL) ISOSURFACES
  !
  ! Authors: Patrick Joeckel, MPICH, Feb 2004
  !          - orinal code
  !          Michael Traub,   MPICH, Jul 2004
  !          - index limitation for iso-surface search and search order
  !            in namelist (kskip_t, kskip_b, lrev)

  USE messy_main_blather_bi, ONLY: start_message_bi, end_message_bi
  USE messy_main_channel_bi, ONLY: n_dom
  USE messy_main_channel,    ONLY: STRLEN_CHANNEL, STRLEN_OBJECT, REPR_UNDEF
  USE messy_viso

  IMPLICIT NONE
  PRIVATE
  SAVE

  ! -------------- ISO - SURFACES --------------------------------
  ! MAX. NUMBER OF ISO-SURFACES IN NAMELIST
  INTEGER, PARAMETER :: NMAXISO        = 200

  TYPE IO_ISOSF
     CHARACTER(LEN=STRLEN_OBJECT-2)  :: name    = ''
     CHARACTER(LEN=STRLEN_CHANNEL)   :: channel  = ''
     CHARACTER(LEN=STRLEN_OBJECT)    :: object = ''
     REAL(DP)                        :: value   = 0.0
     LOGICAL                         :: lfrac   = .false.
     LOGICAL                         :: lrev    = .false.
     INTEGER                         :: kskip_t = 0
     INTEGER                         :: kskip_b = 0
  END TYPE IO_ISOSF

  TYPE ISOSF
     TYPE(IO_ISOSF)                     :: io
     INTEGER                            :: reprid = REPR_UNDEF
     REAL(DP), DIMENSION(:,:,:),POINTER :: ptr   => NULL() ! channel object
     REAL(DP), DIMENSION(:,:),  POINTER :: index => NULL() ! index of iso-value
     ! fraction of box 'below' iso-surface
     REAL(DP), DIMENSION(:,:),  POINTER :: frac  => NULL()
     INTEGER                            :: domain_idx    = 0
  END TYPE ISOSF

  TYPE(IO_ISOSF), DIMENSION(NMAXISO)        :: ISO
  INTEGER                                   :: NISO = 0
  TYPE(ISOSF),    DIMENSION(:), ALLOCATABLE :: XISO
  ! -------------- END ISO - SURFACES ----------------------------

  ! -------------- MAPS ------------------------------------------
  ! MAX. NUMBER OF MAPS IN NAMELIST
  INTEGER, PARAMETER :: NMAXMAP        = 300

  TYPE IO_MAPSF
     CHARACTER(LEN=STRLEN_OBJECT)   :: name    = ''
     CHARACTER(LEN=STRLEN_CHANNEL)  :: sfcha   = ''
     CHARACTER(LEN=STRLEN_OBJECT)   :: sfobj   = ''
     CHARACTER(LEN=STRLEN_CHANNEL)  :: channel = ''
     CHARACTER(LEN=STRLEN_OBJECT)   :: object  = ''
  END TYPE IO_MAPSF

  TYPE MAPSF
     TYPE(IO_MAPSF)                     :: io
     REAL(DP), DIMENSION(:,:),  POINTER :: map    => NULL() ! map
     REAL(DP), DIMENSION(:,:),  POINTER :: index  => NULL() ! index
     REAL(DP), DIMENSION(:,:),  POINTER :: frac   => NULL() ! fraction below
     REAL(DP), DIMENSION(:,:,:),POINTER :: src    => NULL() ! channel object
     INTEGER                            :: reprid = REPR_UNDEF
     LOGICAL                            :: llev = .FALSE.
     INTEGER                            :: lev  = 0
     INTEGER                            :: domain_idx    = 0
  END TYPE MAPSF

  TYPE(IO_MAPSF), DIMENSION(NMAXMAP)        :: MAP
  INTEGER                                   :: NMAP = 0
  TYPE(MAPSF),    DIMENSION(:), ALLOCATABLE :: XMAP 
  ! -------------- END MAPs --------------------------------------

  PUBLIC :: viso_initialize
  PUBLIC :: viso_init_coupling
  PUBLIC :: viso_vdiff
  PUBLIC :: viso_local_end
#ifdef COSMO
  PUBLIC :: viso_global_end
#endif
  PUBLIC :: viso_free_memory
  !PRIVATE :: viso_read_nml_cpl

CONTAINS

! ---------------------------------------------------------------------
  SUBROUTINE viso_initialize

    ! viso MODULE ROUTINE (ECHAM-5 INTERFACE)
    !
    ! INITIALIZATION OF GLOBAL VARIABLES FROM NAMELIST
    ! IN PARALLEL ENVIRONMENT
    !
    ! Author: Patrick Joeckel, MPICH, Feb 20034

    USE messy_main_mpi_bi,     ONLY: p_parallel_io, p_io, p_bcast
    USE messy_main_blather_bi, ONLY: error_bi
    USE messy_main_grid_def_mem_bi,ONLY: nlev
    USE messy_main_tools,      ONLY: find_next_free_unit &
                                   , domains_from_string

    IMPLICIT NONE

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'viso_initialize'
    INTEGER                     :: iou    ! I/O unit
    INTEGER                     :: status ! error status
    INTEGER                     :: i, j
    CHARACTER(LEN=32)                        :: varname = '  '
    INTEGER,           DIMENSION(:), POINTER :: domnum  => NULL()
    INTEGER                                  :: nd, num

    ! INITIALIZE CPL
    IF (p_parallel_io) THEN
       iou = find_next_free_unit(100,200)
       CALL viso_read_nml_cpl(status, iou)
       IF (status /= 0) CALL error_bi('', substr)
    END IF

    CALL start_message_bi(modstr,'INITIALISATION ',substr)

    ! ----------- ISO - SURFACES -----------------------------
    IF (p_parallel_io) THEN
       ! GET NUMBER OF ISO-SURFACES
       NISO = 1
       DO i=1, NMAXISO
          IF (TRIM(ISO(i)%name)    == '') CYCLE
          IF (TRIM(ISO(i)%channel) == '') CYCLE
          IF (TRIM(ISO(i)%object)  == '') CYCLE
          ! analyze name for domains
          CALL domains_from_string(status,ISO(i)%name, n_dom, num)
          IF (status /= 0) CALL error_bi('error in namelist parsing (ISO)' &
               , substr)
          NISO = NISO + num
       END DO
       NISO = NISO - 1
    END IF
    CALL p_bcast(NISO, p_io)
    ALLOCATE(XISO(NISO))

    IF (p_parallel_io) THEN
       NISO = 1
       ! CHECK AND COPY DATA
       DO i=1, NMAXISO
          IF (TRIM(ISO(i)%name) == '') CYCLE

          CALL domains_from_string(status, ISO(i)%name, n_dom, num &
                ,varname, dnums=domnum)
          IF (status /= 0) CALL error_bi('error in namelist parsing (ISO)' &
               , substr)

          ! XISO(NISO)%io%name = ISO(i)%name
          domain_loop1: DO nd = 1, SIZE(domnum)
             XISO(NISO)%io%name    = TRIM(varname)
             XISO(NISO)%domain_idx = domnum(nd)
             WRITE(*,*) 'ISO-SURFACE: ''',TRIM(XISO(NISO)%io%name),''''
             !
             IF (TRIM(ISO(i)%channel) == '') THEN
                WRITE(*,*) ' ... empty channel ... skipping ...'
                CYCLE
             ELSE
                WRITE(*,*) '    CHANNEL: ',TRIM(ISO(i)%channel)
                XISO(NISO)%io%channel = ISO(i)%channel
             END IF
             !
             IF (TRIM(ISO(i)%object) == '') THEN
                WRITE(*,*) ' ... empty channel object ... skipping ...'
                CYCLE
             ELSE
                WRITE(*,*) '    OBJECT : ',TRIM(ISO(i)%object)
                XISO(NISO)%io%object = ISO(i)%object
             END IF
             !
             XISO(NISO)%io%value = ISO(i)%value
             WRITE(*,*) '    VALUE  : ', XISO(NISO)%io%value
             ! 
             XISO(NISO)%io%lfrac = ISO(i)%lfrac
             IF (XISO(NISO)%io%lfrac) THEN
                WRITE(*,*) '    METHOD : LINEAR INTERPOLATION'
             ELSE
                WRITE(*,*) '    METHOD : NEAREST NEIGHBOUR'
             END IF
             !

             XISO(NISO)%io%lrev = ISO(i)%lrev
             IF (XISO(NISO)%io%lrev) THEN
                WRITE(*,*) 'SEARCH ORD.: REVERSE'
             ELSE
                WRITE(*,*) 'SEARCH ORD.: NORMAL'
             END IF
             
             XISO(NISO)%io%kskip_b = ISO(i)%kskip_b
             XISO(NISO)%io%kskip_t = ISO(i)%kskip_t
             !
             IF ((XISO(NISO)%io%kskip_t < 0) .OR. &
                  (XISO(NISO)%io%kskip_t) > nlev-1) THEN
                WRITE(*,*) 'NUMBER OF LEVELS TO SKIP FROM THE TOP IS OUT OF RANGE'
                WRITE(*,*) ' -> RE-SET TO ZERO (= SKIPPING IGNORED)'
                XISO(NISO)%io%kskip_t = 0
             END IF
             !
             IF ((XISO(NISO)%io%kskip_b < 0) .OR. &
                  (XISO(NISO)%io%kskip_b) > nlev-1) THEN
                WRITE(*,*) 'NUMBER OF LEVELS TO SKIP FROM THE BOTTOM'//&
                     &' IS OUT OF RANGE'
                WRITE(*,*) ' -> RE-SET TO ZERO (= SKIPPING IGNORED)'
                XISO(NISO)%io%kskip_b = 0
             END IF
             !
             IF (XISO(NISO)%io%kskip_b+XISO(NISO)%io%kskip_t > nlev) THEN 
                WRITE(*,*) 'REQUESTED NUMBER OF LEVELS TO BE SKIPPED'//&
                     &' EXCEEDS TOTAL'
                WRITE(*,*) ' -> RE-SET TO ZERO (= SKIPPING IGNORED)'
                XISO(NISO)%io%kskip_t = 0
                XISO(NISO)%io%kskip_b = 0
             END IF
             !
             WRITE(*,*) 'SKIP AT SFC: ', XISO(NISO)%io%kskip_b
             WRITE(*,*) 'SKIP AT TOA: ', XISO(NISO)%io%kskip_t
             
             ! NEXT ISO-SURFACE
             NISO = NISO + 1
             WRITE(*,*) '------------------------------------------------------'
             
          END DO domain_loop1
          DEALLOCATE(domnum)
          NULLIFY(domnum)
       END DO
       NISO = NISO - 1 
    END IF
    CALL p_bcast(NISO, p_io)
    IF (NISO /= SIZE(XISO)) &
         CALL error_bi('error determing number of iso-surfaces' , substr)

    ! BROADCAST RESULTS
    DO i=1, NISO
       CALL p_bcast(XISO(i)%io%name,    p_io)
       CALL p_bcast(XISO(i)%io%channel, p_io)
       CALL p_bcast(XISO(i)%io%object,  p_io)
       CALL p_bcast(XISO(i)%io%value,   p_io)
       CALL p_bcast(XISO(i)%io%lfrac,   p_io)
       CALL p_bcast(XISO(i)%io%lrev,    p_io)
       CALL p_bcast(XISO(i)%io%kskip_t, p_io)
       CALL p_bcast(XISO(i)%io%kskip_b, p_io)
       CALL p_bcast(XISO(i)%domain_idx,  p_io)
    END DO
    ! ----------- END ISO - SURFACES -------------------------

    ! ----------- MAPS ---------------------------------------
    IF (p_parallel_io) THEN
       ! GET NUMBER OF MAPS
       NMAP = 1
       DO i=1, NMAXMAP
          IF (TRIM(MAP(i)%name)    == '') CYCLE
          IF (TRIM(MAP(i)%sfcha)   == '') CYCLE
          IF (TRIM(MAP(i)%sfobj)   == '') CYCLE
          IF (TRIM(MAP(i)%channel) == '') CYCLE
          IF (TRIM(MAP(i)%object)  == '') CYCLE
          CALL domains_from_string(status, TRIM(MAP(i)%name), n_dom, num)
          IF (status /= 0) CALL error_bi('error in namelist parsing (MAP)' &
               , substr)
          NMAP = NMAP + num
       END DO
       NMAP = NMAP - 1
    END IF
    CALL p_bcast(NMAP, p_io)
    ALLOCATE(XMAP(NMAP))

    IF (p_parallel_io) THEN
       ! GET NUMBER OF MAPS
       NMAP = 1
       ! CHECK AND COPY DATA
       DO i=1, NMAXMAP
          IF (TRIM(MAP(i)%name) == '') CYCLE

          CALL domains_from_string(status,MAP(i)%name,n_dom,num &
               ,varname, dnums=domnum) 
          IF (status /= 0) CALL error_bi('error in namelist parsing (MAP)' &
               , substr)

          domain_loop2: DO nd = 1, SIZE(domnum)

             XMAP(NMAP)%io%name     = MAP(i)%name
             XMAP(NMAP)%domain_idx  = domnum(nd)

             WRITE(*,*) 'MAP                  : '''&
                  ,TRIM(XMAP(NMAP)%io%name),''''

             IF (TRIM(MAP(i)%sfcha) == '') THEN
                WRITE(*,*) ' ... empty channel for SURFACE ... skipping ...'
                CYCLE
             ELSE
                WRITE(*,*) '    CHANNEL (SURFACE): ',TRIM(MAP(i)%sfcha)
                XMAP(NMAP)%io%sfcha = MAP(i)%sfcha
             END IF
             !
             IF (TRIM(MAP(i)%sfobj) == '') THEN
                WRITE(*,*) ' ... empty channel object FOR SURFACE'//&
                     &' ... skipping ...'
                CYCLE
             ELSE
                WRITE(*,*) '    OBJECT  (SURFACE): ',TRIM(MAP(i)%sfobj)
                XMAP(NMAP)%io%sfobj = MAP(i)%sfobj
             END IF
             !
             IF (TRIM(MAP(i)%channel) == '') THEN
                WRITE(*,*) ' ... empty channel ... skipping ...'
                CYCLE
             ELSE
                WRITE(*,*) '    CHANNEL          : ',TRIM(MAP(i)%channel)
                XMAP(NMAP)%io%channel = MAP(i)%channel
             END IF
             !
             IF (TRIM(MAP(i)%object) == '') THEN
                WRITE(*,*) ' ... empty channel object ... skipping ...'
                CYCLE
             ELSE
                WRITE(*,*) '    OBJECT           : ',TRIM(MAP(i)%object)
                XMAP(NMAP)%io%object = MAP(i)%object
             END IF
             !
             ! NEXT MAP-SURFACE
             NMAP = NMAP + 1
             WRITE(*,*) '------------------------------------------------------'
          END DO domain_loop2
          DEALLOCATE(domnum)
          NULLIFY(domnum)
       END DO
       NMAP = NMAP - 1 
    END IF
    CALL p_bcast(NMAP, p_io)
    IF (NMAP /= SIZE(XMAP)) &
         CALL error_bi('error in MAP namelist interpretation', substr) 
    
    ! BROADCAST RESULTS
    DO i=1, NMAP
       CALL p_bcast(XMAP(i)%io%name,     p_io)
       CALL p_bcast(XMAP(i)%io%sfcha,    p_io)
       CALL p_bcast(XMAP(i)%io%sfobj,    p_io)
       CALL p_bcast(XMAP(i)%io%channel,  p_io)
       CALL p_bcast(XMAP(i)%io%object,   p_io)
       CALL p_bcast(XMAP(i)%domain_idx,  p_io)
    END DO
    ! ----------- END MAPS -----------------------------------

    ! It is possible now that viso_init_coupling is called several times.
    ! This implies that the uniqueness of channel object names is not longer
    ! tested there ... it must be checked here already!
    !
    ! Note: names of maps and iso-surfaces might be identical, because
    !       '_i' and/or '_f' is appended to the iso-surfaces -> no check
    !       required
    !
    ! names of iso-surfaces ...
    DO i=1, NISO
       DO j=1, i-1
          IF (TRIM(XISO(j)%io%name) == TRIM(XISO(i)%io%name) &
               .AND. XISO(j)%domain_idx == XISO(i)%domain_idx &
               ) THEN
             CALL error_bi('object name '//&
                  &TRIM(XISO(i)%io%name) //' is not unique', substr)
          END IF
       END DO
    END DO
    !
    ! names of maps ...
    DO i=1, NMAP
       DO j=1, i-1
          IF (TRIM(XMAP(j)%io%name) == TRIM(XMAP(i)%io%name) &
               .AND. XMAP(j)%domain_idx == XMAP(i)%domain_idx &
               ) THEN
             CALL error_bi('object name '//&
                  &TRIM(XMAP(i)%io%name) //' is not unique', substr)
          END IF
       END DO
    END DO

    IF (p_parallel_io) THEN
       WRITE(*,*) ' ---> ',NISO,' ISO-SURFACE(S) REQUESTED !'
       WRITE(*,*) ' ---> ',NMAP,' MAP(S) REQUESTED !'
    END IF

    CALL end_message_bi(modstr,'INITIALISATION ',substr)

  END SUBROUTINE viso_initialize
! ---------------------------------------------------------------------

! ---------------------------------------------------------------------
  SUBROUTINE viso_init_coupling

    ! viso MODULE ROUTINE (ECHAM-5 INTERFACE)
    !
    ! define specific channel and objects(s)
    !
    ! Author: Patrick Joeckel, MPICH, Feb 2004

    USE messy_main_mpi_bi,           ONLY: p_parallel_io
    USE messy_main_blather_bi,       ONLY: warning_bi, info_bi, error_bi
    USE messy_main_channel_error_bi, ONLY: channel_halt
    USE messy_main_channel_bi,       ONLY: GP_3D_MID, GP_3D_INT &
                                         , GP_2D_HORIZONTAL
    USE messy_main_grid_def_mem_bi,  ONLY: nlev
    USE messy_main_channel,          ONLY: new_channel, new_channel_object   &
                                         , new_attribute, get_channel_object &
                                         , get_attribute                     &
                                         , get_channel_object_info           &
                                         , get_channel_info
    USE messy_main_channel_mem,      ONLY: dom_curid
    USE messy_main_constants_mem,    ONLY: STRLEN_ULONG, FLAGGED_BAD
    USE messy_main_tools,            ONLY: str2num, int2str

    IMPLICIT NONE

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER       :: substr = 'viso_init_coupling'
    INTEGER                           :: status
    INTEGER                           :: i
    CHARACTER(LEN=50)                 :: vstr
    INTEGER                           :: reprid
    CHARACTER(LEN=STRLEN_ULONG)       :: att_unit
    LOGICAL                           :: lfirst
    CHARACTER(LEN=3)                  :: str = '   '
    CHARACTER(LEN=2)                  :: sf = '  '
    
    CALL start_message_bi(modstr, 'CHANNEL / OBJECT DEFINITIONS', substr)

    ! viso_init_coupling is called several times, make sure that
    ! channel is only created once
    CALL get_channel_info(status, modstr)
    lfirst = (status /= 0)
    IF (lfirst) THEN
       CALL new_channel(status, modstr, reprid=GP_2D_HORIZONTAL)
       CALL channel_halt(substr, status)
    END IF

    ! ------------- ISO-SURFACES ----------------------------------------
    DO i=1, NISO  ! LOOP OVER ISO-SURFACES
       IF (XISO(i)%domain_idx /= dom_curid) CYCLE

       ! ADD CHANNEL OBJECTS FOR ISO-SURFACES
       IF (p_parallel_io) THEN
          WRITE(*,*) 'ISO-SURFACE: ''',TRIM(XISO(i)%io%name),''''
       END IF

       ! GET REQUESTED CHANNEL OBJECT
       CALL get_channel_object(status &
            , TRIM(XISO(i)%io%channel), TRIM(XISO(i)%io%object) &
            , p3=XISO(i)%ptr)
       IF (status /= 0) THEN
          CALL warning_bi(&
               '   '//TRIM(XISO(i)%io%channel)//' - '//&
               &TRIM(XISO(i)%io%object)//&
               &' not found ... skipping', substr)
          NULLIFY(XISO(i)%ptr)
          CYCLE
       END IF

       ! CHECK REPRESENTATION
       CALL info_bi( &
            '  trying '//TRIM(XISO(i)%io%channel)//' - '//&
            &TRIM(XISO(i)%io%object) &
            , substr)
       CALL get_channel_object_info(status &
            , TRIM(XISO(i)%io%channel), TRIM(XISO(i)%io%object) &
            , reprid=reprid)
       CALL channel_halt(substr, status)
       !
       IF ( (reprid /= GP_3D_MID) .AND. (reprid /= GP_3D_INT) ) THEN
          CALL warning_bi('  representation not supported ... skipping' &
               , substr)
          NULLIFY(XISO(i)%ptr)
          CYCLE
       END IF
       XISO(i)%reprid = reprid
       CALL info_bi('  ... OK', substr)

       ! GET UNITS ATTRIBUTE
       att_unit=''
       CALL get_attribute(status  &
            , TRIM(XISO(i)%io%channel), TRIM(XISO(i)%io%object) &
            , 'units', c=att_unit)

       WRITE(vstr,*) XISO(i)%io%value

       ! ADD CHANNEL OBJECT (INDEX)
       CALL get_channel_object_info(status, modstr, TRIM(XISO(i)%io%name)//'_i')
       IF (status == 0) THEN
          CALL info_bi(TRIM(XISO(i)%io%name)//'_i already defined ...', substr)
       ELSE
          CALL new_channel_object(status, modstr, TRIM(XISO(i)%io%name)//'_i' &
               , p2=XISO(i)%index )
          CALL channel_halt(substr, status)
          CALL new_attribute(status, modstr, TRIM(XISO(i)%io%name)//'_i' &
               , 'long_name', c='index of '//TRIM(XISO(i)%io%object)//'='&
               &//TRIM(vstr)//' '//TRIM(att_unit) )
          CALL channel_halt(substr, status)
          CALL new_attribute(status, modstr, TRIM(XISO(i)%io%name)//'_i' &
               , 'units', c=' ' )
          CALL channel_halt(substr, status)
       END IF

       ! ADD CHANNEL OBJECT (FRACTION 'BELOW')
       IF (XISO(i)%io%lfrac) THEN
          CALL get_channel_object_info(status, modstr &
               , TRIM(XISO(i)%io%name)//'_f')
          IF (status == 0) THEN
             CALL info_bi(TRIM(XISO(i)%io%name)//'_f already defined ...' &
                  , substr)
          ELSE
             CALL new_channel_object(status, modstr &
                  , TRIM(XISO(i)%io%name)//'_f' &
                  , p2=XISO(i)%frac )
             CALL channel_halt(substr, status)
             CALL new_attribute(status, modstr, TRIM(XISO(i)%io%name)//'_f' &
                  , 'long_name', c='fraction of box below '&
                  &//TRIM(XISO(i)%io%object)//'='&
                  &//TRIM(vstr)//' '//TRIM(att_unit) )
             CALL channel_halt(substr, status)
             CALL new_attribute(status, modstr, TRIM(XISO(i)%io%name)//'_f' &
                  , 'units', c=' ' )
             CALL channel_halt(substr, status)
          END IF
       END IF

    END DO ! LOOP OVER ISO-SURFACES
    ! ------------- END ISO-SURFACES ------------------------------------


    ! ------------- MAPS ------------------------------------------------
    DO i=1, NMAP  ! LOOP OVER MAPS
       IF (XMAP(i)%domain_idx /= dom_curid) CYCLE

       ! ADD CHANNEL OBJECTS FOR MAP
       IF (p_parallel_io) THEN
          WRITE(*,*) 'MAP        : ''',TRIM(XMAP(i)%io%name),''''
       END IF

       ! GET REQUESTED CHANNEL OBJECT (SURFACE)
       ! SPECIAL CASE: select specific level
       IF (TRIM(ADJUSTL(XMAP(i)%io%sfcha)) == '#level') THEN
          XMAP(i)%llev = .TRUE.
          CALL str2num(XMAP(i)%io%sfobj, XMAP(i)%lev, status)
          IF (status /=0) THEN
             CALL error_bi('READ ERROR IN #level INFORMATION',substr)
          ENDIF
          ! SPECIAL CASES
          SELECT CASE (XMAP(i)%lev)
          CASE(-1)
             XMAP(i)%lev = nlev
          CASE(-2)
             XMAP(i)%lev = nlev+1
          END SELECT
          CALL int2str(str, XMAP(i)%lev)
          CALL info_bi(' level '//str, substr)
          !
          ! NO INDEX FIELD
          NULLIFY(XMAP(i)%index)
          ! NO FRACTION BELOW
          NULLIFY(XMAP(i)%frac)
       ELSE

          CALL get_channel_object(status &
               , TRIM(XMAP(i)%io%sfcha), TRIM(XMAP(i)%io%sfobj)//'_i' &
               , p2=XMAP(i)%index)
          ! allow also index fields without suffix '_i' ... 
          IF (status == 0) THEN
             sf = '_i'
          ELSE
             CALL get_channel_object(status &
                  , TRIM(XMAP(i)%io%sfcha), TRIM(XMAP(i)%io%sfobj) &
                  , p2=XMAP(i)%index)
             sf = '  '
          ENDIF

          IF (status /= 0) THEN
             CALL warning_bi(&
                  '   '//TRIM(XMAP(i)%io%sfcha)//' - '//&
                  &TRIM(XMAP(i)%io%sfobj)//TRIM(sf)//&
                  &' (SURFACE) not found ... skipping', substr)
             NULLIFY(XMAP(i)%map)
             NULLIFY(XMAP(i)%index)
             NULLIFY(XMAP(i)%frac)
             NULLIFY(XMAP(i)%src)
             CYCLE
          ELSE
             CALL info_bi( &
                  '   '//TRIM(XMAP(i)%io%sfcha)//' - '//&
                  &TRIM(XMAP(i)%io%sfobj)//TRIM(sf)//&
                  &' (SURFACE)     found', substr)
          END IF

          ! GET REQUESTED CHANNEL OBJECT (SURFACE; FRACTION)
          CALL get_channel_object(status &
               , TRIM(XMAP(i)%io%sfcha), TRIM(XMAP(i)%io%sfobj)//'_f' &
               , p2=XMAP(i)%frac )
          IF (status /= 0) THEN
             CALL info_bi( &
                  '   '//TRIM(XMAP(i)%io%sfcha)//' - '//&
                  &TRIM(XMAP(i)%io%sfobj)//'_f'//&
                  &' not found  ---> NEAREST NEIGHBOUR', substr)
             NULLIFY(XMAP(i)%frac)
          ELSE
             CALL info_bi( &
                  '   '//TRIM(XMAP(i)%io%sfcha)//' - '//&
                  &TRIM(XMAP(i)%io%sfobj)//'_f'//&
                  &'     found  ---> LINEAR INTERPOLATION', substr)
          END IF
       END IF

       ! GET REQUESTED CHANNEL OBJECT
       CALL get_channel_object(status &
            , TRIM(XMAP(i)%io%channel), TRIM(XMAP(i)%io%object) &
            , p3=XMAP(i)%src )
       IF (status /= 0) THEN
          CALL warning_bi( &
               '   '//TRIM(XMAP(i)%io%channel)//' - '//&
               &TRIM(XMAP(i)%io%channel)//&
               &' not found ... skipping', substr)
          NULLIFY(XMAP(i)%map)
          NULLIFY(XMAP(i)%index)
          NULLIFY(XMAP(i)%frac)
          NULLIFY(XMAP(i)%src)
          CYCLE
       END IF

       ! CHECK REPRESENTATION
       CALL info_bi( &
            '  trying '//TRIM(XMAP(i)%io%channel)//' - '//&
            &TRIM(XMAP(i)%io%object) &
            , substr)
       CALL get_channel_object_info(status &
            , TRIM(XMAP(i)%io%channel), TRIM(XMAP(i)%io%object) &
            , reprid=reprid )
       IF ( (reprid /= GP_3D_MID) .AND. (reprid /= GP_3D_INT) ) THEN
          CALL warning_bi('  representation not supported ... skipping' &
               , substr)
          NULLIFY(XMAP(i)%map)
          NULLIFY(XMAP(i)%index)
          NULLIFY(XMAP(i)%frac)
          NULLIFY(XMAP(i)%src)
          CYCLE
       END IF
       XMAP(i)%reprid = reprid
       ! SPECIAL FOR SELECTED LEVELS
       IF (XMAP(i)%llev) THEN
          CALL int2str(str, XMAP(i)%lev)
          CALL info_bi(' level '//str, substr)
          IF ( ( (XMAP(i)%reprid == GP_3D_MID) .AND. (XMAP(i)%lev > nlev)) &
               .OR. &
               ( (XMAP(i)%reprid == GP_3D_INT) .AND. (XMAP(i)%lev > nlev+1)) &
               .OR. &
               (XMAP(i)%lev < 1) ) THEN
             CALL error_bi('LEVEL OUT OF VALID RANGE',substr)
          ENDIF
       END IF
       CALL info_bi('  ... OK', substr)

       ! ADD CHANNEL OBJECT (MAP)
       CALL get_channel_object_info(status, modstr, TRIM(XMAP(i)%io%name))
       IF (status == 0) THEN
          CALL info_bi(TRIM(XMAP(i)%io%name)//' already defined ...' &
               , substr)
       ELSE
          CALL new_channel_object(status, modstr, TRIM(XMAP(i)%io%name) &
               , p2= XMAP(i)%map)
          CALL channel_halt(substr, status)
          CALL new_attribute(status, modstr, TRIM(XMAP(i)%io%name) &
               , 'long_name', c=TRIM(XMAP(i)%io%object)//' at '//&
               &TRIM(XMAP(i)%io%sfobj) )
          CALL channel_halt(substr, status)
          CALL new_attribute(status, modstr, TRIM(XMAP(i)%io%name) &
               , 'missing_value', r=FLAGGED_BAD)
          CALL channel_halt(substr, status)
          ! TRANSFER UNIT
          CALL get_attribute(status  &
               , TRIM(XMAP(i)%io%channel), TRIM(XMAP(i)%io%object) &
               , 'units', c=att_unit)
          IF (status == 0) THEN
             CALL new_attribute(status, modstr, TRIM(XMAP(i)%io%name) &
                  , 'units', c=TRIM(att_unit))
             CALL channel_halt(substr, status)
          END IF
       END IF

    END DO ! LOOP OVER MAPS
    ! ------------- END MAPS --------------------------------------------

    CALL end_message_bi(modstr, 'CHANNEL / OBJECT DEFINITIONS', substr)

  END SUBROUTINE viso_init_coupling
! ---------------------------------------------------------------------

! ---------------------------------------------------------------------
  SUBROUTINE viso_vdiff

    ! CALCULATE INDEX (AND FRACTION OF BOX BELOW) ISO-SURFACE(S)

    USE messy_main_channel_bi,  ONLY: GP_3D_MID
    USE messy_main_grid_def_mem_bi, ONLY: jrow, kproma, nlev, nproma
    USE messy_main_tools,       ONLY: iso2ind
    USE messy_main_channel_mem, ONLY: dom_curid

    IMPLICIT NONE

    INTEGER  :: i, jp, k
    INTEGER  :: kmin, kmax
    REAL(DP) :: col(nproma,nlev)
    
    ! LOOP OVER ISO-SURFACES
    DO i=1, NISO
       IF (XISO(i)%domain_idx /= dom_curid) CYCLE
       IF (.NOT.ASSOCIATED(XISO(i)%ptr)) CYCLE

       kmin = 1    + XISO(i)%io%kskip_t
       kmax = nlev - XISO(i)%io%kskip_b

       IF (XISO(i)%reprid == GP_3D_MID) THEN
          col(1:kproma,:) = XISO(i)%ptr(_RI_XYZ__(1:kproma,jrow,:)) ! :,jrow
       ELSE
          ! GP_3D_INT (anything else is ruled out, see above)
          col(1:kproma,:) = 0.5_dp * ( &
               XISO(i)%ptr(_RI_XYZ__(1:kproma,jrow,1:nlev)) +  & ! 1:nlev,jrow
               XISO(i)%ptr(_RI_XYZ__(1:kproma,jrow,2:nlev+1)) )  ! 2:nlev+1,jrow
       END IF

       IF (XISO(i)%io%lfrac) THEN

          DO jp=1, kproma
             CALL iso2ind(col(jp, kmin:kmax)       & ! kmin:kmax,jrow
                  , XISO(i)%io%value               &
                  , k, f=XISO(i)%frac(jp, jrow), lrev=XISO(i)%io%lrev )
             ! ADD OFFSET SKIPPED IN iso2ind TO GET CORRECT INDEX
             XISO(i)%index(jp, jrow) = REAL(k, DP) + &
                  REAL(XISO(i)%io%kskip_t)
          END DO

       ELSE

          DO jp=1, kproma
             CALL iso2ind(col(jp, kmin:kmax)       &
                  , XISO(i)%io%value               & 
                  , k, lrev=XISO(i)%io%lrev)
             ! ADD OFFSET SKIPPED IN iso2ind TO GET CORRECT INDEX
             XISO(i)%index(jp, jrow) = REAL(k, DP) + &
                  REAL(XISO(i)%io%kskip_t)
          END DO

       END IF

    END DO

  END SUBROUTINE viso_vdiff
! ---------------------------------------------------------------------

! ---------------------------------------------------------------------
  SUBROUTINE viso_local_end(jjrow)

    ! CALCULATE VALUES ON (ISO-) SURFACE(S), USING INDEX
    ! AND OPTIONALLY THE FRACTION OF THE BOX BELOW

    USE messy_main_channel_bi,  ONLY: GP_3D_MID
    USE messy_main_grid_def_mem_bi, ONLY: jrow_ext=>jrow, kproma, nproma, nlev
    USE messy_main_tools,       ONLY: ind2val
    USE messy_main_channel_mem, ONLY: dom_curid

    IMPLICIT NONE

    INTEGER, INTENT(IN), OPTIONAL :: jjrow

    INTEGER  :: i, jp, jrow
    REAL(DP) :: col(nproma,nlev)

    IF (PRESENT(jjrow)) THEN
       jrow = jjrow
    ELSE
       jrow = jrow_ext
    END IF

#ifdef MESSYIDTC
    CALL viso_vdiff
#endif

    ! LOOP OVER MAPS
    DO i=1, NMAP
       IF (XMAP(i)%domain_idx /= dom_curid) CYCLE
       IF (.NOT.ASSOCIATED(XMAP(i)%src)) CYCLE

       IF (XMAP(i)%reprid == GP_3D_MID) THEN
          col(1:kproma,:) = XMAP(i)%src(_RI_XYZ__(1:kproma,jrow,:))
       ELSE
          ! GP_3D_INT (anything else is ruled out, see above)
          col(1:kproma,:) = 0.5_dp * (         &
               XMAP(i)%src(_RI_XYZ__(1:kproma,jrow,1:nlev)) +  &
               XMAP(i)%src(_RI_XYZ__(1:kproma,jrow,2:nlev+1)) )
       ENDIF

       IF (ASSOCIATED(XMAP(i)%frac)) THEN

          DO jp=1, kproma
             CALL ind2val(XMAP(i)%map(jp, jrow),  col(jp,:) &
                  , INT(XMAP(i)%index(jp, jrow)), XMAP(i)%frac(jp, jrow) )
          END DO

       ELSE

          IF (XMAP(i)%llev) THEN
             DO jp=1, kproma
                CALL ind2val(XMAP(i)%map(jp, jrow) &
                     , XMAP(i)%src(_RI_XYZ__(jp,jrow,:)) &
                     , XMAP(i)%lev)
             END DO
          ELSE
             DO jp=1, kproma
                CALL ind2val(XMAP(i)%map(jp, jrow), col(jp,:) &
                     , INT(XMAP(i)%index(jp, jrow)) )
             END DO
          END IF

       END IF

    END DO

  END SUBROUTINE viso_local_end
! ---------------------------------------------------------------------
#ifdef COSMO
  SUBROUTINE viso_global_end

    USE messy_main_grid_def_mem_bi, ONLY: ngpblks

    IMPLICIT NONE

    INTEGER :: jrow 

    DO jrow = 1, ngpblks
       CALL viso_local_end(jrow)
    END DO

  END SUBROUTINE viso_global_end
#endif
! ---------------------------------------------------------------------
  
! ---------------------------------------------------------------------
  SUBROUTINE viso_free_memory

    IMPLICIT NONE

    INTEGER :: i

    DO i=1, NISO
       NULLIFY(XISO(i)%ptr)
       NULLIFY(XISO(i)%index)
       NULLIFY(XISO(i)%frac)
    END DO
    DEALLOCATE(XISO)

    DO i=1, NMAP
       NULLIFY(XMAP(i)%map)
       NULLIFY(XMAP(i)%index)
       NULLIFY(XMAP(i)%frac)
       NULLIFY(XMAP(i)%src)
    END DO
    DEALLOCATE(XMAP)

  END SUBROUTINE viso_free_memory
! ---------------------------------------------------------------------
  
! ----------------------------------------------------------------------
  SUBROUTINE viso_read_nml_cpl(status, iou)

    ! viso MODULE ROUTINE (ECHAM-5 INTERFACE, PRIVATE)
    !
    ! read namelist for 'coupling' to ECHAM5
    !
    ! Author: Patrick Joeckel, MPICH, Feb 2004

    USE messy_main_tools, ONLY: read_nml_open, read_nml_check, read_nml_close

    IMPLICIT NONE

    ! I/O
    INTEGER, INTENT(OUT) :: status     ! error status
    INTEGER, INTENT(IN)  :: iou        ! I/O unit

    ! (LOCAL) NAMELIST VARIABLES
    CHARACTER(LEN=*), PARAMETER :: substr = 'viso_read_nml_cpl'

    NAMELIST /CPL/ ISO, MAP

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

    ! DIAGNOSE NAMELIST AND SET GLOBAL SWITCHES
    !
    ! CHECK NAMELIST

    CALL read_nml_close(substr, iou, modstr)

    status = 0  ! no ERROR

  END SUBROUTINE viso_read_nml_cpl
! ----------------------------------------------------------------------

! **************************************************************************
END MODULE messy_viso_si
! **************************************************************************

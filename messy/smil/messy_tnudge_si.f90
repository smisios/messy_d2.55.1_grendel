#include "messy_main_ppd_bi.inc"

!***********************************************************************
MODULE messy_tnudge_si
!***********************************************************************

  ! MODULE FOR TRACER NUDGING  (ECHAM5-MESSy INTERFACE (MESSy-SMIL))
  !
  ! Authors:
  !    Patrick Joeckel, MPICH, December 2003
  !    Patrick Joeckel, DLR,   February 2020
  !
  ! TODO:
  !    - IF MESSYTENDENCY: register only tracers which will be altered

  USE messy_main_blather_bi,    ONLY: start_message_bi, end_message_bi &
                                    , error_bi
#ifdef MESSYTENDENCY
 USE messy_main_tendency_bi,   ONLY: mtend_get_handle,       &
                                     mtend_get_start_l,      &
                                     mtend_add_l,            &
                                     mtend_register,         &
                                     mtend_id_tracer
#endif
  USE messy_main_constants_mem, ONLY: DP, STRLEN_MEDIUM
  USE messy_main_channel,       ONLY: STRLEN_CHANNEL, STRLEN_OBJECT
  USE messy_main_tracer,        ONLY: STRLEN_TRSET
  USE messy_tnudge

  IMPLICIT NONE
  PRIVATE
  SAVE

  INTRINSIC :: NULL, TRIM, ASSOCIATED, MIN, MAX, ABS, NINT, ADJUSTL

  TYPE T_IO_NUDGE
     CHARACTER(LEN=10*STRLEN_TRSET+10)  :: trset  = ''  ! name of tracer set
                                                        ! (and domain)
     CHARACTER(LEN=2*STRLEN_MEDIUM+1) :: trname   = ''  ! name of tracer
     CHARACTER(LEN=STRLEN_CHANNEL) :: channel     = ''  ! name of source-channel
     CHARACTER(LEN=STRLEN_OBJECT ) :: object      = ''  ! name of source-object
     REAL(DP)                      :: coeff  =   0.0_DP ! nudging coeff [s]
     REAL(DP)                      :: latmin = -90.0_DP
     REAL(DP)                      :: latmax =  90.0_DP
     INTEGER                       :: kmin   =     0
     INTEGER                       :: kmax   =     0
     REAL(DP)                      :: lonmin =   0.0_DP
     REAL(DP)                      :: lonmax = 360.0_DP
     LOGICAL                       :: lflux  = .FALSE. ! flux diagnostic ?
     !
     ! channel name of isosurface for kmin
     CHARACTER(LEN=STRLEN_CHANNEL) :: cviso1 = ''
     ! object name of isosurface for kmin
     CHARACTER(LEN=STRLEN_OBJECT)  :: oviso1 = ''
     ! channel name of isosurface for kmax
     CHARACTER(LEN=STRLEN_CHANNEL) :: cviso2 = ''
     ! object name of isosurface for kmax
     CHARACTER(LEN=STRLEN_OBJECT)  :: oviso2 = ''
     ! mask: 0 = everywhere, 1 = over land, 2 = over sea
     INTEGER                       :: land   = 0
     !
     ! VALID RANGE
     REAL(DP), DIMENSION(2)        :: vr = (/ -HUGE(0.0_dp), HUGE(0.0_dp) /)
     !
  END TYPE T_IO_NUDGE

  TYPE T_XNUDGE
     TYPE(T_IO_NUDGE)                    :: io
     INTEGER                             :: idt   = 0        ! tracer index
     INTEGER                             :: type  = 0        ! tracer type
     REAL(DP), DIMENSION(:,:,:), POINTER :: field => NULL()  ! field to nudge
     INTEGER                             :: ndim             ! rank of field
     INTEGER                             :: nfdim            ! rank of flux
     REAL(DP), DIMENSION(:,:,:), POINTER :: flux => NULL()   ! diag. flux
     !
     ! index of kmin if isosurfaces are used
     REAL(DP), DIMENSION(:,:), POINTER   :: iviso1 => NULL()
     ! index of kmax if isosurfaces are used
     REAL(DP), DIMENSION(:,:), POINTER   :: iviso2 => NULL()
     ! account for fractions below the isosurface for kmin
     LOGICAL                             :: lfrac1 = .false.
     ! account for fractions below the isosurface for kmax
     LOGICAL                             :: lfrac2 = .false.
     ! fraction below the isosurface for kmin
     REAL(DP), DIMENSION(:,:), POINTER   :: fviso1 => NULL()
     ! fraction below the isosurface for kmax
     REAL(DP), DIMENSION(:,:), POINTER   :: fviso2 => NULL()
     !
  END TYPE T_XNUDGE

  ! NUDGING ELEMENTS
  ! MAX. NUMBER OF NUDGING FIELDS IN NAMELIST
  INTEGER, PARAMETER                        :: NMAXNUDGE = 800
  TYPE(T_IO_NUDGE), DIMENSION(NMAXNUDGE)    :: TNUDGE
  TYPE(T_XNUDGE),   DIMENSION(:), POINTER   :: XNUDGE => NULL()
  INTEGER                                   :: NNUDGE

  ! TROPOPAUSE
  CHARACTER(LEN=STRLEN_OBJECT)           :: tropopause(2)
  LOGICAL                                :: ltp   ! tropopause available ?
  ! tropopause level index
  REAL(DP), DIMENSION(:,:), POINTER      :: itp => NULL()

  ! BOUNDARY LAYER
  CHARACTER(LEN=STRLEN_OBJECT)           :: boundarylayer(2)
  LOGICAL                                :: lbl   ! boundary layer available ?
  ! boundary layer level index
  REAL(DP), DIMENSION(:,:), POINTER      :: ibl => NULL()

#ifdef ECHAM5
  ! LAGRANGIAN REPRESENTATION (ATTILA)
  LOGICAL                           :: L_LG = .TRUE.
  REAL(DP), DIMENSION(:),   POINTER :: iplon  => NULL() ! index
  REAL(DP), DIMENSION(:),   POINTER :: iplat  => NULL() ! index
  REAL(DP), DIMENSION(:),   POINTER :: plon   => NULL() ! degrees east
  REAL(DP), DIMENSION(:),   POINTER :: plat   => NULL() ! degrees north
  REAL(DP), DIMENSION(:),   POINTER :: iplev  => NULL() ! index
  REAL(DP),                 POINTER :: amcell => NULL() ! kg
#endif
  REAL(DP), DIMENSION(:,:), POINTER :: garea  => NULL() ! m^2

#ifdef ECHAM5
  ! LAGRANGIAN REPRESENTATION (CLaMS)
  LOGICAL                           :: L_CL = .TRUE.
  REAL(DP), DIMENSION(:),   POINTER :: LAT     => NULL()
  REAL(DP), DIMENSION(:),   POINTER :: LON     => NULL()
  REAL(DP), DIMENSION(:,:), POINTER :: POS      => NULL()
#endif

#ifdef MESSYTENDENCY
  INTEGER :: my_handle
#endif

  PUBLIC :: tnudge_initialize
  !           ! -> tnudge_input_check
  PUBLIC :: tnudge_init_memory
  !           ! -> tnudge_setup_mem
  PUBLIC :: tnudge_init_coupling
  !           ! -> tnudge_setup_cpl
  PUBLIC :: tnudge_global_start
  PUBLIC :: tnudge_local_end
  PUBLIC :: tnudge_global_end
  PUBLIC :: tnudge_free_memory
  !PRIVATE :: tnudge_read_nml_cpl

CONTAINS

! ---------------------------------------------------------------------
  SUBROUTINE tnudge_initialize

    ! tnudge MODULE ROUTINE (ECHAM-5 INTERFACE)
    !
    ! INITIALIZATION OF GLOBAL VARIABLES FROM NAMELIST
    ! IN PARALLEL ENVIRONMENT
    !
    ! Author: Patrick Joeckel, MPICH, Dec 2003

    ! ECHAM5/MESSy
    USE messy_main_mpi_bi,         ONLY: p_parallel_io, p_io, p_bcast
#ifdef ECHAM5
    USE messy_main_tracer_mem_bi,  ONLY: NGCELL
    USE messy_main_channel_bi,     ONLY: REPR_LG_CLAMS,&
                                         REPR_UNDEF
#endif
    USE messy_main_tools,          ONLY: find_next_free_unit

    IMPLICIT NONE

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'tnudge_initialize'
    INTEGER                     :: iou    ! I/O unit
    INTEGER                     :: status ! error status
    INTEGER                     :: i

    ! INIT

#ifdef ECHAM5
    ! LAGRANGIAN ?
    L_LG = (NGCELL > 0)
    L_CL = (REPR_LG_CLAMS /= REPR_UNDEF)
#endif

#ifdef MESSYTENDENCY
    my_handle = mtend_get_handle(modstr)
#endif

    ! INITIALIZE CPL
    IF (p_parallel_io) THEN
       iou = find_next_free_unit(100,200)
       CALL tnudge_read_nml_cpl(status, iou)
       IF (status /= 0) CALL error_bi(' ',substr)
    END IF

    ! BROADCAST RESULTS
    CALL p_bcast(tropopause(1), p_io)
    CALL p_bcast(tropopause(2), p_io)
    CALL p_bcast(boundarylayer(1), p_io)
    CALL p_bcast(boundarylayer(2), p_io)

    ! --------------------------
    ! GRID-PPOINT REPRESENTATION
    ! --------------------------
    CALL start_message_bi(modstr, 'INITIALISATION', substr)
    !
    IF (p_parallel_io) &
         CALL tnudge_input_check(status, TNUDGE, XNUDGE, NNUDGE)
    CALL p_bcast(status, p_io)
    IF (status /= 0) CALL error_bi('ERROR IN TNUDGE &CPL NAMELIST', substr)

    CALL p_bcast(NNUDGE, p_io)
    IF (.NOT. p_parallel_io) ALLOCATE(XNUDGE(NNUDGE))

    ! BROADCAST RESULTS
    DO i=1, NNUDGE
       CALL p_bcast(XNUDGE(i)%io%trset,   p_io)
       CALL p_bcast(XNUDGE(i)%io%trname,  p_io)
       CALL p_bcast(XNUDGE(i)%io%channel, p_io)
       CALL p_bcast(XNUDGE(i)%io%object,  p_io)
       CALL p_bcast(XNUDGE(i)%io%coeff,   p_io)
       CALL p_bcast(XNUDGE(i)%io%latmin,  p_io)
       CALL p_bcast(XNUDGE(i)%io%latmax,  p_io)
       CALL p_bcast(XNUDGE(i)%io%kmin,    p_io)
       CALL p_bcast(XNUDGE(i)%io%kmax,    p_io)
       CALL p_bcast(XNUDGE(i)%io%lonmin,  p_io)
       CALL p_bcast(XNUDGE(i)%io%lonmax,  p_io)
       CALL p_bcast(XNUDGE(i)%io%lflux,   p_io)
       !
       CALL p_bcast(XNUDGE(i)%io%cviso1, p_io)
       CALL p_bcast(XNUDGE(i)%io%oviso1, p_io)
       CALL p_bcast(XNUDGE(i)%io%cviso2, p_io)
       CALL p_bcast(XNUDGE(i)%io%oviso2, p_io)
       CALL p_bcast(XNUDGE(i)%io%land,   p_io)
       !
       CALL p_bcast(XNUDGE(i)%io%vr,     p_io)
       !
       ! nfdim, flux
       ! -> set in tnudge_init_memory
       ! idt, type, field, ndim
       ! -> set in tnudge_init_coupling
!!$       CALL  p_bcast(XNUDGE(i)%domain_idx, p_io)
    END DO
    !
    IF (p_parallel_io) THEN
       WRITE(*,*) ' ---> ',NNUDGE,' NUDGING ELEMENT(S) INITIALIZED !'
    END IF
    !
    CALL end_message_bi(modstr, 'INITIALISATION', substr)

  END SUBROUTINE tnudge_initialize
! ---------------------------------------------------------------------

! ---------------------------------------------------------------------
  SUBROUTINE tnudge_init_memory

    USE messy_main_mpi_bi,           ONLY: p_parallel_io
    USE messy_main_channel_error_bi, ONLY: channel_halt
    USE messy_main_channel_bi,       ONLY: GP_3D_MID, GP_3D_1LEV
    USE messy_main_constants_mem,    ONLY: STRLEN_MEDIUM
    USE messy_main_channel,          ONLY: new_channel, new_channel_object &
                                         , new_attribute
    USE messy_main_tracer_mem_bi,    ONLY: GPTRSTR
#ifdef ECHAM5
    USE messy_main_tracer_mem_bi,    ONLY: LGTRSTR, CLTRSTR
#endif

    IMPLICIT NONE

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'tnudge_init_memory'

#ifdef MESSYTENDENCY
    CALL mtend_register(my_handle, mtend_id_tracer)
#endif

    CALL start_message_bi(modstr, 'CHANNEL DEFINITION (GP)', substr)
    CALL tnudge_setup_mem(XNUDGE, NNUDGE, GPTRSTR)
    CALL end_message_bi(modstr, 'CHANNEL DEFINITION (GP)', substr)

#ifdef ECHAM5
    IF (L_LG) THEN
       CALL start_message_bi(modstr, 'CHANNEL DEFINITION (LG)', substr)
       CALL tnudge_setup_mem(XNUDGE, NNUDGE, LGTRSTR)
       CALL end_message_bi(modstr, 'CHANNEL DEFINITION (LG)', substr)
    END IF
#endif

#ifdef ECHAM5
    IF (L_CL) THEN
       CALL start_message_bi(modstr, 'CHANNEL DEFINITION (CL)', substr)
       CALL tnudge_setup_mem(XNUDGE, NNUDGE, CLTRSTR)
       CALL end_message_bi(modstr, 'CHANNEL DEFINITION (CL)', substr)
    END IF
#endif

  CONTAINS

    ! ---------------------------------------------------------------------
    SUBROUTINE tnudge_setup_mem(XNUDGE, NNUDGE, ch)

      IMPLICIT NONE

      ! I/O
      TYPE(T_XNUDGE), DIMENSION(:), INTENT(INOUT) :: XNUDGE
      INTEGER,                      INTENT(IN)    :: NNUDGE
      CHARACTER(LEN=*),             INTENT(IN)    :: ch

      ! LOCAL
      CHARACTER(LEN=*), PARAMETER :: substr = 'tnudge_setup_mem'
      INTEGER                     :: status
      INTEGER                     :: i
      LOGICAL                     :: lfirst
      INTEGER                     :: reprid
      CHARACTER(LEN=2*STRLEN_MEDIUM+1) :: name, tname
      CHARACTER(LEN=STRLEN_OBJECT) :: oname ='' ! op_pj_20120125
      CHARACTER(LEN=17)           :: unit

      lfirst = .TRUE.

      DO i=1, NNUDGE
         ! ONLY CURRENT TRACER SET
         IF (TRIM(XNUDGE(i)%io%trset) /= TRIM(ch)) CYCLE

         ! DIAGNOSTIC FLUX REQUESTED ?
         IF (.NOT. XNUDGE(i)%io%lflux) CYCLE
         ! SET DIMENSION OF FLUX
         IF ((XNUDGE(i)%io%kmax - XNUDGE(i)%io%kmin) == 0) THEN
            IF ((XNUDGE(i)%io%kmin == -4) .AND. &
                 (TRIM(XNUDGE(i)%io%oviso1) /= TRIM(XNUDGE(i)%io%oviso2))) THEN
               XNUDGE(i)%nfdim = 3
               reprid = GP_3D_MID
               unit = 'molecules/(m^3 s)'
            ELSE
               XNUDGE(i)%nfdim = 2
               reprid = GP_3D_1LEV
               unit = 'molecules/(m^2 s)'
            END IF
         ELSE
            XNUDGE(i)%nfdim = 3
            reprid = GP_3D_MID
            unit = 'molecules/(m^3 s)'
         END IF
         ! CREATE CHANNEL / OBJECT
         IF (lfirst) THEN
            IF (p_parallel_io) &
                 WRITE(*,*) '... creating new channel '//modstr//'_'//ch
            CALL new_channel(status, modstr//'_'//ch)
            CALL channel_halt(substr, status)
            lfirst = .FALSE.
         END IF
         name = ''
         name = TRIM(XNUDGE(i)%io%channel)//'_'//TRIM(XNUDGE(i)%io%object)
         tname = ''
         tname = TRIM(ADJUSTL(XNUDGE(i)%io%trname))

         ! resolve naming conflict for multi-layer nudging
         oname = TRIM(tname)//'_'//TRIM(name)//'_flx'

         IF (p_parallel_io) &
              WRITE(*,*) '... adding channel object '//TRIM(oname)
         CALL new_channel_object(status, modstr//'_'//ch &
              , TRIM(oname), reprid=reprid, p3=XNUDGE(i)%flux)
         CALL channel_halt(substr, status)
         CALL new_attribute(status, modstr//'_'//ch, TRIM(oname) &
              , 'long_name', c='flux '//TRIM(name)//' into '//TRIM(tname) )
         CALL channel_halt(substr, status)
         CALL new_attribute(status, modstr//'_'//ch, TRIM(oname) &
              , 'units',c ='flux of '//TRIM(unit) )
         CALL channel_halt(substr, status)
      END DO

    END SUBROUTINE tnudge_setup_mem
    ! ---------------------------------------------------------------------

  END SUBROUTINE tnudge_init_memory
! ---------------------------------------------------------------------

! ---------------------------------------------------------------------
  SUBROUTINE tnudge_global_start

#ifdef ECHAM5

    ! ECHAM5/MESSy
    USE messy_main_grid_def_bi,   ONLY: gboxarea_2d
    USE messy_main_transform_bi,  ONLY: trp_gpdc_gpgl

    IMPLICIT NONE

    ! LOCAL
    LOGICAL, SAVE  :: lfirst = .TRUE.

    IF (.NOT. L_LG .and..not. L_CL) RETURN

    IF (lfirst) THEN
       CALL trp_gpdc_gpgl(1, gboxarea_2d, garea)
       lfirst = .FALSE.
    END IF
#endif
  END SUBROUTINE tnudge_global_start
! ---------------------------------------------------------------------

! ---------------------------------------------------------------------
  SUBROUTINE tnudge_init_coupling

    ! TNUDGE MODULE ROUTINE (ECHAM5 INTERFACE)
    !
    ! set pointers to specific channel(s) and tracer(s)
    !
    ! Author: Patrick Joeckel, MPICH, Dec 2003

    ! ECHAM5/MESSy
    USE messy_main_mpi_bi,           ONLY: p_parallel_io
    USE messy_main_blather_bi,       ONLY: info_bi, error_bi, warning_bi
    USE messy_main_tracer_mem_bi,    ONLY: GPTRSTR
#ifdef ECHAM5
    USE messy_main_tracer_mem_bi,    ONLY: LGTRSTR, CLTRSTR
    USE messy_main_channel_error_bi, ONLY: channel_halt
#endif
    USE messy_main_grid_def_mem_bi,  ONLY: nlev
    USE messy_main_channel_bi,       ONLY: GP_3D_MID &
                                         , GP_3D_1LEV, GP_2D_HORIZONTAL
    ! MESSy
    USE messy_main_tracer,           ONLY: get_tracer
    USE messy_main_channel,          ONLY: get_channel_object &
                                         , get_channel_object_info

    IMPLICIT NONE

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'tnudge_init_coupling'
    INTEGER                     :: status
    INTEGER                     :: i

    CALL start_message_bi(modstr, 'COUPLING', substr)

    ! (1a) GET POINTERS TO TROPOPAUSE INDEX
    ltp = .false.
    IF (p_parallel_io) THEN
       WRITE(*,*) '  INITIALIZING TROPOPAUSE LEVEL INDEX ...'
    END IF
    IF (TRIM(tropopause(1)) /= '') THEN
       CALL get_channel_object(status &
            , TRIM(tropopause(1)), TRIM(tropopause(2)), p2=itp)
       IF (status /= 0) THEN
          CALL warning_bi(&
               TRIM(tropopause(1))//' - '//TRIM(tropopause(2))//&
               &' not found', substr)
       ELSE
          CALL info_bi( &
               TRIM(tropopause(1))//' - '//TRIM(tropopause(2))//&
               &' initialized', substr)
          ltp = .TRUE.
       END IF
    END IF

    ! (1b) GET POINTERS TO PBLH-INDEX, ETC.
    lbl = .false.
    IF (p_parallel_io) THEN
       WRITE(*,*) '  INITIALIZING BOUNDARY LAYER LEVEL INDEX ...'
    END IF
    IF (TRIM(boundarylayer(1)) /= '') THEN
       CALL get_channel_object(status &
            , TRIM(boundarylayer(1)), TRIM(boundarylayer(2)), p2=ibl)
       IF (status /= 0) THEN
          CALL warning_bi(&
               TRIM(boundarylayer(1))//' - '//TRIM(boundarylayer(2))//&
               &' not found' , substr)
       ELSE
          CALL info_bi(&
               TRIM(boundarylayer(1))//' - '//TRIM(boundarylayer(2))//&
               &' initialized', substr)
          lbl = .TRUE.
       END IF
    END IF

#ifdef ECHAM5
    IF (L_LG) THEN
       IF (p_parallel_io) THEN
          WRITE(*,*) '  COUPLING ATTILA CHANNEL OBJECTS (LG):'
       END IF
       CALL get_channel_object(status, 'attila', 'IPLAT', p1=iplat)
       CALL channel_halt(substr, status)
       CALL get_channel_object(status, 'attila', 'IPLON', p1=iplon)
       CALL channel_halt(substr, status)
       CALL get_channel_object(status, 'attila', 'PLAT', p1=plat)
       CALL channel_halt(substr, status)
       CALL get_channel_object(status, 'attila', 'PLON', p1=plon)
       CALL channel_halt(substr, status)
       CALL get_channel_object(status, 'attila', 'IPLEV', p1=iplev)
       CALL channel_halt(substr, status)
       CALL get_channel_object(status, 'attila', 'AMCELL', p0=amcell)
       CALL channel_halt(substr, status)
    END IF

    IF (L_CL) THEN
       IF (p_parallel_io) THEN
          WRITE(*,*) '  COUPLING CLaMS CHANNEL OBJECTS (CL):'
       END IF
       CALL get_channel_object(status, 'clams', 'LAT_OLD', p1=LAT)
       CALL channel_halt(substr, status)
       CALL get_channel_object(status, 'clams', 'LON_OLD', p1=LON)
       CALL channel_halt(substr, status)
       CALL get_channel_object(status, 'clams', 'POS', p2=pos)
       CALL channel_halt(substr, status)
    ENDIF
#endif
    IF (p_parallel_io) THEN
       WRITE(*,*) ' '
       WRITE(*,*) '  INITIALIZING TRACER-INDICES '//&
            &'AND POINTER TO CHANNEL OBJECTS (GP):'
    END IF
    !
    CALL tnudge_setup_cpl(XNUDGE, NNUDGE, GPTRSTR)
    CALL tnudge_cpl_iso(XNUDGE, NNUDGE, GPTRSTR)

#ifdef ECHAM5
    IF (L_LG)  THEN
       IF (p_parallel_io) THEN
          WRITE(*,*) ' '
          WRITE(*,*) '  INITIALIZING TRACER-INDICES '//&
               &'AND POINTER TO CHANNEL OBJECTS (LG):'
       END IF
       !
       CALL tnudge_setup_cpl(XNUDGE, NNUDGE, LGTRSTR)
       CALL tnudge_cpl_iso(XNUDGE, NNUDGE, LGTRSTR)
    END IF

    IF (L_CL)  THEN
       IF (p_parallel_io) THEN
          WRITE(*,*) ' '
          WRITE(*,*) '  INITIALIZING TRACER-INDICES '//&
               &'AND POINTER TO CHANNEL OBJECTS (CL):'
       END IF
       !
       CALL tnudge_setup_cpl(XNUDGE, NNUDGE, CLTRSTR)
       CALL tnudge_cpl_iso(XNUDGE, NNUDGE, CLTRSTR)
    END IF

#endif

    ! check, if all tracers have been found
    DO i=1, NNUDGE
       IF (XNUDGE(i)%idt == 0) THEN
          CALL error_bi( &
               'TRACER '//TRIM(XNUDGE(i)%io%trname)//' of  '&
               &'TRACER SET '//TRIM(XNUDGE(i)%io%trset)//' not found!' &
               , substr)
       END IF
    END DO

    CALL end_message_bi(modstr, 'COUPLING', substr)

  CONTAINS

    ! -----------------------------------------
    SUBROUTINE tnudge_setup_cpl(XNUDGE, NNUDGE, TRSTR)

      USE messy_main_tracer,          ONLY: full2base_sub
      USE messy_main_tracer_tools_bi, ONLY: tracer_halt
      USE messy_main_channel,         ONLY: STRLEN_CHANNEL, STRLEN_OBJECT

      IMPLICIT NONE

      ! I/O
      TYPE(T_XNUDGE), DIMENSION(:), INTENT(INOUT) :: XNUDGE
      INTEGER,                      INTENT(IN)    :: NNUDGE
      CHARACTER(LEN=*),             INTENT(IN)    :: TRSTR

      ! LOCAL
      INTEGER                 :: ierr
      INTEGER                 :: i
      INTEGER                 :: kmin, kmax
      CHARACTER(LEN=STRLEN_CHANNEL+STRLEN_OBJECT+2) :: kminstr, kmaxstr
      INTEGER                 :: reprid
      CHARACTER(LEN=STRLEN_MEDIUM) :: basename    = '' ! name of tracer
      CHARACTER(LEN=STRLEN_MEDIUM) :: subname     = '' ! OPTIONAL subname

      ! (2) GET INDEX OF TRACER AND TENDENCY
      !
      DO i=1, NNUDGE
         ! ONYL CURRENT TRACER SET
         IF (TRIM(XNUDGE(i)%io%trset) /= TRIM(TRSTR)) CYCLE

         ! TRACER AND TENDENCY
         CALL full2base_sub(status, TRIM(XNUDGE(i)%io%trname) &
              , basename, subname)
         CALL tracer_halt(substr, status)
         !
         CALL get_tracer(ierr, TRSTR, TRIM(basename) &
              , subname = TRIM(subname)              &
              , idx = XNUDGE(i)%idt                  &
              , type = XNUDGE(i)%type)
         IF (ierr /= 0) THEN
            CALL error_bi('   '//TRIM(XNUDGE(i)%io%trname)//&
                 &' - tracer not found!', substr)
         ELSE
            CALL info_bi('   '//TRIM(XNUDGE(i)%io%trname)//&
                 &' - tracer found!', substr)
         END IF

         ! CHANNEL OBJECT TO NUDGE
         NULLIFY(XNUDGE(i)%field)
         !
         CALL get_channel_object(status    &
              , TRIM(XNUDGE(i)%io%channel) &
              , TRIM(XNUDGE(i)%io%object)  &
              , p3=XNUDGE(i)%field )
         IF (status /= 0) THEN
            CALL error_bi('    '//&
                 &TRIM(XNUDGE(i)%io%channel)//' - '//&
                 &TRIM(XNUDGE(i)%io%object)//&
                 &' not found ' , substr)
            NULLIFY(XNUDGE(i)%field)
         END IF

         ! SAVE 'REAL' RANK OF CHANNEL OBJECT
         ! (QUASI-2D FIELDS WITH 1 LEVEL)
         CALL get_channel_object_info(status &
              , TRIM(XNUDGE(i)%io%channel), TRIM(XNUDGE(i)%io%object) &
              , reprid=reprid )
         IF (reprid == GP_3D_MID) THEN
            XNUDGE(i)%ndim = 3
         ELSEIF (reprid == GP_2D_HORIZONTAL) THEN
            XNUDGE(i)%ndim = 2
         ELSEIF (reprid == GP_3D_1LEV) THEN
            XNUDGE(i)%ndim = 2
         ELSE
            CALL error_bi( '    '//&
                 &'representation not supported', ' ')
         END IF

         ! LEVELS
         !
         ! -4  : iso surface     ! mz_ho_20150409
         ! -3  : boundary layer
         ! -2  : tropopause
         ! -1  : top layer
         !  0  : bottom layer
         !  n>0: level n
         !
         kmin = XNUDGE(i)%io%kmin
         kmax = XNUDGE(i)%io%kmax
         !
         SELECT CASE(kmin)
         CASE(-4)
            kminstr = TRIM(XNUDGE(i)%io%cviso1)//'::'//TRIM(XNUDGE(i)%io%oviso1)
         CASE(-3)
            kminstr = 'boundary layer'
            IF (.NOT.lbl) THEN
               CALL error_bi( &
                    '    '//'- boundary layer not available', ' ')
            END IF
         CASE(-2)
            kminstr = 'tropopause'
            IF (.NOT.ltp) THEN
               CALL error_bi(&
                    '    '//'- tropopause not available', ' ')
            END IF
         CASE(-1)
            kminstr = 'top'
         CASE(0)
            kminstr = 'surface'
         CASE DEFAULT
            kmin = MAX(kmin,1)
            kmin = MIN(kmin,nlev)
            WRITE(kminstr,'(i4)') kmin
            XNUDGE(i)%io%kmin = kmin
         END SELECT

         SELECT CASE(kmax)
         CASE(-4)
            kmaxstr = TRIM(XNUDGE(i)%io%cviso2)//'::'//TRIM(XNUDGE(i)%io%oviso2)
         CASE(-3)
            kmaxstr = 'boundary layer'
            IF (.NOT.lbl) THEN
               CALL error_bi( &
                    '    '//'boundary layer not available', ' ')
            END IF
         CASE(-2)
            kmaxstr = 'tropopause'
            IF (.NOT.ltp) THEN
               CALL error_bi( &
                    '    '//'tropopause not available', ' ')
            END IF
         CASE(-1)
            kmaxstr = 'top'
         CASE(0)
            kmaxstr = 'surface'
         CASE DEFAULT
            kmax = MIN(kmax,nlev)
            kmax = MAX(kmax,1)
            WRITE(kmaxstr,'(i4)') kmax
            XNUDGE(i)%io%kmax = kmax
         END SELECT

         ! DIAG OUTPUT
         IF (p_parallel_io) THEN
            WRITE(*,*) '   '//TRIM(XNUDGE(i)%io%trname), ' <- '&
                 ,TRIM(XNUDGE(i)%io%object)            &
                 ,' (',TRIM(XNUDGE(i)%io%channel),')'    &
                 , ' RANK = ',XNUDGE(i)%ndim
            WRITE(*,*) '    LATITUDE : ', XNUDGE(i)%io%latmin,' -> ' &
                 ,XNUDGE(i)%io%latmax
            WRITE(*,*) '    LONGITUDE: ', XNUDGE(i)%io%lonmin,' -> ' &
                 ,XNUDGE(i)%io%lonmax
            WRITE(*,*) '    LEVELS   : ', TRIM(kminstr),' -> ', TRIM(kmaxstr)
         END IF

      END DO ! i=1, NNUDGE

    END SUBROUTINE tnudge_setup_cpl
    ! ---------------------------------------------------------------------

    ! ---------------------------------------------------------------------
    SUBROUTINE tnudge_cpl_iso(XNUDGE, NNUDGE, TRSTR)

      IMPLICIT NONE

      ! I/O
      TYPE(T_XNUDGE), DIMENSION(:), INTENT(INOUT) :: XNUDGE
      INTEGER,                      INTENT(IN)    :: NNUDGE
      CHARACTER(LEN=*),             INTENT(IN)    :: TRSTR

      ! LOCAL
      CHARACTER(LEN=*), PARAMETER :: substr = 'tnudge_cpl_iso'
      INTEGER :: i

      DO i=1, NNUDGE  ! LOOP OVER NUDGING FIELDS

         ! ONLY ONE TRACER SET
         IF (TRIM(XNUDGE(i)%io%trset) /= TRIM(TRSTR)) CYCLE

         IF (XNUDGE(i)%io%kmin == -4) THEN
            CALL get_channel_object(status, &
                 TRIM(XNUDGE(i)%io%cviso1), &
                 TRIM(XNUDGE(i)%io%oviso1)//'_i' &
                 , p2 = XNUDGE(i)%iviso1)
            IF (status /= 0) THEN
               CALL error_bi( &
                    '   TNUDGE CHANNEL::OBJECT '//&
                    &TRIM(XNUDGE(i)%io%cviso1)//'::'//&
                    &TRIM(XNUDGE(i)%io%oviso1)//'_i'//' for tracer '//&
                    &TRIM(XNUDGE(i)%io%trname)//&
                    &' not found', substr)
            ELSE
               CALL get_channel_object(status, &
                    TRIM(XNUDGE(i)%io%cviso1), &
                    TRIM(XNUDGE(i)%io%oviso1)//'_f' &
                    , p2 = XNUDGE(i)%fviso1)
               IF (status /= 0) THEN
                  CALL info_bi( &
                       '   TNUDGE CHANNEL::OBJECT '//&
                       &TRIM(XNUDGE(i)%io%cviso1)//'::'//&
                       &TRIM(XNUDGE(i)%io%oviso1)//'_i'//&
                       &' is used as kmin for tracer '//&
                       &TRIM(XNUDGE(i)%io%trname) &
                       , substr)
                  NULLIFY(XNUDGE(i)%fviso1)
                  XNUDGE(i)%lfrac1 = .FALSE.
               ELSE
                  CALL info_bi( &
                       '   TNUDGE CHANNEL::OBJECTs '//&
                       &TRIM(XNUDGE(i)%io%cviso1)//'::'//&
                       &TRIM(XNUDGE(i)%io%oviso1)//'_i and _f'//&
                       &' are used as kmin for tracer '//&
                       &TRIM(XNUDGE(i)%io%trname) &
                       , substr)
                   XNUDGE(i)%lfrac1 = .TRUE.
                END IF
             END IF
          END IF

          IF (XNUDGE(i)%io%kmax == -4) THEN
             CALL get_channel_object(status, &
                  TRIM(XNUDGE(i)%io%cviso2), &
                  TRIM(XNUDGE(i)%io%oviso2)//'_i' &
                  , p2 = XNUDGE(i)%iviso2)
             IF (status /= 0) THEN
                CALL error_bi( &
                     '   TNUDGE CHANNEL::OBJECT '//&
                     &TRIM(XNUDGE(i)%io%cviso2)//'::'//&
                     &TRIM(XNUDGE(i)%io%oviso2)//'_i'//' for tracer '//&
                     &TRIM(XNUDGE(i)%io%trname)//&
                     &' not found', substr)
             ELSE
                CALL get_channel_object(status, &
                     TRIM(XNUDGE(i)%io%cviso2), &
                     TRIM(XNUDGE(i)%io%oviso2)//'_f' &
                     , p2 = XNUDGE(i)%fviso2)
                IF (status /= 0) THEN
                   CALL info_bi( &
                        '   TNUDGE CHANNEL::OBJECT '//&
                        &TRIM(XNUDGE(i)%io%cviso2)//'::'//&
                        &TRIM(XNUDGE(i)%io%oviso2)//'_i'//&
                        &' is used as kmax for tracer '//&
                        &TRIM(XNUDGE(i)%io%trname) &
                        , substr)
                   NULLIFY(XNUDGE(i)%fviso2)
                   XNUDGE(i)%lfrac2 = .FALSE.
                ELSE
                   CALL info_bi( &
                       '   TNUDGE CHANNEL::OBJECTs '//&
                       &TRIM(XNUDGE(i)%io%cviso2)//'::'//&
                       &TRIM(XNUDGE(i)%io%oviso2)//'_i and _f'//&
                       &' are used as kmax for tracer '//&
                       &TRIM(XNUDGE(i)%io%trname) &
                       , substr)
                   XNUDGE(i)%lfrac2 = .TRUE.
                END IF
             END IF
          END IF
       END DO

    END SUBROUTINE tnudge_cpl_iso
    ! ---------------------------------------------------------------------

  END SUBROUTINE tnudge_init_coupling
  ! ---------------------------------------------------------------------

  ! ---------------------------------------------------------------------
  SUBROUTINE tnudge_local_end

    USE messy_main_tracer_bi,     ONLY: main_tracer_fconv_loc
    USE messy_main_tracer_mem_bi, ONLY: GPTRSTR
#ifndef MESSYTENDENCY
    USE messy_main_tracer_mem_bi, ONLY: qxtte, qxtm1
#endif
    USE messy_main_blather_bi,      ONLY: info_bi
    USE messy_main_grid_def_mem_bi, ONLY: nlev, jrow, kproma
    USE messy_main_grid_def_bi,     ONLY: philat_2d, philon_2d, deltaz
    USE messy_main_data_bi,         ONLY: slf, rho_air_dry_3d
    USE messy_main_data_bi,         ONLY: l2tls

    ! MESSy
#ifdef ICON
    USE messy_main_channel,       ONLY: get_channel_object
#endif
    USE messy_main_timer,         ONLY: lstart, timer_get_time_step_len
    USE messy_main_tracer,        ONLY: SINGLE, FAMILY
    USE messy_main_constants_mem, ONLY: N_A, M_air, g

    IMPLICIT NONE

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'tnudge_local_end'
    REAL(DP), PARAMETER         :: uconv = N_A * 1.0e3_DP/M_air
    INTEGER  :: i
    INTEGER  :: idt
    INTEGER  :: jp, jk
    INTEGER  :: kmin, kmax, klev
    INTEGER  :: kmin_tmp, kmax_tmp
    LOGICAL  :: llon, llat, llev, lrange
    INTEGER  :: jtype
    INTEGER  :: klev_f            ! level in diagnostic flux
    REAL(dp), DIMENSION(:,:), ALLOCATABLE :: vstart
    REAL(dp), DIMENSION(:,:), ALLOCATABLE :: vtend
    REAL(DP) :: convert           ! unit conversion
    REAL(DP), DIMENSION(:,:), POINTER     :: zairdens => NULL()
    REAL(DP), DIMENSION(:,:), ALLOCATABLE :: zdz
    REAL(DP) :: ztmst_hard        ! actual time step for 'hard' nudging
    REAL(DP) :: zcoeff            ! nudging coeff.
#ifdef ICON
    INTEGER  :: status
#endif
    ! INIT
    CALL timer_get_time_step_len(ztmst_hard)
    IF (lstart .AND. .NOT. l2tls) ztmst_hard = 2._dp * ztmst_hard

    ALLOCATE(vstart(kproma,nlev))
    ALLOCATE(vtend(kproma,nlev))

    ! CALCULATE AIR DENSITY AND LAYER THICKNESS
    ALLOCATE(zdz(kproma, nlev))
    zairdens => rho_air_dry_3d(_RI_XYZ__(1:kproma,jrow,1:nlev))

#ifdef ICON
    ! (1a) GET POINTERS TO TROPOPAUSE INDEX
    ltp = .false.
    IF (TRIM(tropopause(1)) /= '') THEN
       CALL get_channel_object(status &
            , TRIM(tropopause(1)), TRIM(tropopause(2)), p2=itp)
        IF (status == 0) ltp = .TRUE.
     END IF
     lbl = .false.
     IF (TRIM(boundarylayer(1)) /= '') THEN
        CALL get_channel_object(status &
            , TRIM(boundarylayer(1)), TRIM(boundarylayer(2)), p2=ibl)
        IF (status == 0) lbl = .TRUE.
     END IF

#endif

    zdz(:,:) = 0.0_dp
    zdz(:,2:nlev) = deltaz(_RI_XYZ__(1:kproma,jrow,2:nlev))

    ! LOOP TWICE: 1st SINGLE TRACERS   (i.e., skip families)
    !             2nd FAMILIES  (i.e., convert t2f, skip tracers, convert f2t)
    tracer_type_loop: DO jtype = 1,2

       IF (jtype == 2) CALL main_tracer_fconv_loc('t2f',substr, GPTRSTR)

    ! (1) FILL FLAG ARRAY
    nudging_loop: DO i=1, NNUDGE  ! LOOP OVER NUDGING FIELDS

       ! ONLY grid point tracer set here
       IF (TRIM(XNUDGE(i)%io%trset) /= TRIM(GPTRSTR)) CYCLE

       ! Channel/OBJECT NOT AVAILABLE
       IF (.NOT.ASSOCIATED(XNUDGE(i)%field)) CYCLE
       ! TRACER NOT AVAILABLE / TROPOPAUSE NOT AVAILABLE / BOUNDARY LAYER
       ! HEIGHT NOT AVAILABLE
       IF (XNUDGE(i)%idt == 0) CYCLE

       idt = XNUDGE(i)%idt

       IF ((jtype == 1) .AND. (XNUDGE(i)%type == FAMILY)) CYCLE
       IF ((jtype == 2) .AND. (XNUDGE(i)%type == SINGLE)) CYCLE

       ! INIT FLUX
       IF (XNUDGE(i)%io%lflux) &
            XNUDGE(i)%flux(_RI_XYZ__(:,jrow,:)) = 0._dp

       ! SET NUDGING COEFF
       zcoeff = XNUDGE(i)%io%coeff
!qqq
       IF (ABS(XNUDGE(i)%io%coeff - ztmst_hard) < EPS) THEN
          IF (lstart .AND. .NOT. l2tls)  & !um_ak_20110629 added for COSMO
               zcoeff = zcoeff / 2.0_dp
       END IF

       ! SET START VALUE AND INIT TENDENCY
#ifdef MESSYTENDENCY
       CALL mtend_get_start_l(idt, v0=vstart)
#else
       vstart(:,:) = qxtm1(_RI_X_ZN_(1:kproma,1:nlev,idt)) +  &
                     qxtte(_RI_X_ZN_(1:kproma,1:nlev,idt))*ztmst_hard
#endif
       vtend(:,:) = 0.0_dp

          level_loop: DO jk=1, nlev  ! LOOP OVER LEVELS

             vector_loop: DO jp=1, kproma

                llon = (philon_2d(jp, jrow) >= XNUDGE(i)%io%lonmin) .AND. &
                       (philon_2d(jp, jrow) <= XNUDGE(i)%io%lonmax)
                llat = (philat_2d(jp, jrow) >= XNUDGE(i)%io%latmin) .AND. &
                       (philat_2d(jp, jrow) <= XNUDGE(i)%io%latmax)

                SELECT CASE(XNUDGE(i)%io%kmin)
                CASE(-4)
                   kmin = NINT(XNUDGE(i)%iviso1(jp, jrow))
                CASE (-3)
                   kmin = NINT(ibl(jp, jrow))
                CASE (-2)
                   kmin = NINT(itp(jp, jrow))
                CASE(-1)
                   kmin = 1
                CASE(0)
                   kmin = nlev
                CASE DEFAULT
                   kmin = XNUDGE(i)%io%kmin
                END SELECT

                SELECT CASE(XNUDGE(i)%io%kmax)
                CASE(-4)
                   kmax = NINT(XNUDGE(i)%iviso2(jp, jrow))
                CASE (-3)
                   kmax = NINT(ibl(jp, jrow))
                CASE (-2)
                   kmax = NINT(itp(jp, jrow))
                CASE(-1)
                   kmax = 1
                CASE(0)
                   kmax = nlev
                CASE DEFAULT
                   kmax = XNUDGE(i)%io%kmax
                END SELECT

                kmin_tmp = MAX(MIN(kmin, kmax),1)
                kmax_tmp = MIN(MAX(kmin, kmax),nlev)
                kmin = kmin_tmp
                kmax = kmax_tmp

                llev = (jk >= kmin) .AND. (jk <= kmax)

                SELECT CASE(XNUDGE(i)%ndim)
                   CASE(2)
                      klev = 1
                   CASE(3)
                      klev = jk
                   CASE DEFAULT
                      CALL info_bi( &
                           ' WARNING! CHANNEL OBJECT HAS RANK /= 2 OR 3 !' &
                           , substr)
                END SELECT

                lrange = ( &
                     ( XNUDGE(i)%field(_RI_XYZ__(jp,jrow,klev)) > &
                     XNUDGE(i)%io%vr(1) ) &
                     .AND. &
                     ( XNUDGE(i)%field(_RI_XYZ__(jp,jrow,klev)) < &
                     XNUDGE(i)%io%vr(2) ) )

                ! PERFORM NUDGING
                IF (llon.AND.llat.AND.llev.AND.lrange) THEN

                   vtend(jp,jk) = -1._DP * ( vstart(jp,jk) - &
                        XNUDGE(i)%field(_RI_XYZ__(jp,jrow,klev))) / zcoeff

                   ! only affect tendency of grid cell fraction above land
                   IF (XNUDGE(i)%io%land == 1) &
                        vtend(jp,jk) = vtend(jp,jk) * slf(jp, jrow)
                   ! only affect tendency of grid cell fraction above water
                   IF (XNUDGE(i)%io%land == 2) &
                        vtend(jp,jk) = vtend(jp,jk) * (1._DP-slf(jp, jrow))

                   IF ((jk == kmax) .AND. XNUDGE(i)%lfrac2) THEN
                      IF ((jk == kmin) .AND. XNUDGE(i)%lfrac1) THEN
                         vtend(jp,jk) = vtend(jp,jk) * &
                              MAX(XNUDGE(i)%fviso1(jp, jrow) - &
                              XNUDGE(i)%fviso2(jp, jrow), 0._DP)
                         ! only affect tendency of grid cell fraction
                         ! between interfaces
                      ELSE
                         vtend(jp,jk) = vtend(jp,jk) * &
                              (1._DP-XNUDGE(i)%fviso2(jp, jrow))
                         ! only affect tendency of grid cell fraction
                         ! above lower interface
                      END IF
                   ELSE
                      IF ((jk == kmin) .AND. XNUDGE(i)%lfrac1) &
                           vtend(jp,jk) = vtend(jp,jk) * &
                           XNUDGE(i)%fviso1(jp, jrow)
                      ! only affect tendency of grid cell fraction
                      ! below upper interface
                   END IF

                   ! SAVE FLUX
                   IF (XNUDGE(i)%io%lflux) THEN
                      ! SET LEVEL
                      IF (XNUDGE(i)%nfdim == 2) THEN
                         klev_f = 1
                         ! mol/mol/s -> molecules/(m^2 s)
                         convert = zairdens(jp, jk) &
                              * zdz(jp, jk) &
                              * uconv
                      ELSE
                         klev_f = jk
                         ! mol/mol/s -> molecules/(m^3 s)
                         convert = zairdens(jp, jk) * uconv
                      END IF
                      !
                      XNUDGE(i)%flux(_RI_XYZ__(jp,jrow,klev_f)) = &
                           vtend(jp,jk) * convert

                   END IF

                END IF

             END DO vector_loop

          END DO level_loop        ! LOOP OVER LEVELS

#ifdef MESSYTENDENCY
          CALL mtend_add_l(my_handle, idt, px=vtend)
#else
          qxtte(_RI_X_ZN_(1:kproma,1:nlev,idt)) = &
               qxtte(_RI_X_ZN_(1:kproma,1:nlev,idt)) + &
               vtend(1:kproma,1:nlev)
#endif

       END DO nudging_loop         ! LOOP OVER NUDGING FIELDS

       IF (jtype == 2) CALL main_tracer_fconv_loc('f2t',substr, GPTRSTR)

    END DO tracer_type_loop        ! LOOP OVER TRACER TYPES

    ! CLEAN UP
    DEALLOCATE(vstart)
    DEALLOCATE(vtend)
    NULLIFY(zairdens)
    DEALLOCATE(zdz)

  END SUBROUTINE tnudge_local_end
  ! ---------------------------------------------------------------------

  ! ---------------------------------------------------------------------
  SUBROUTINE tnudge_global_end

#ifdef ECHAM5
    USE messy_main_tracer_mem_bi, ONLY: qxtte_a, qxtm1_a, NCELL
    USE messy_main_tracer_mem_bi, ONLY: qxtte_c, qxt_c
    USE messy_main_tracer_mem_bi, ONLY: dnparts
    USE messy_clams_global,       ONLY: cmcell
    USE messy_main_blather_bi,    ONLY: info_bi
    USE messy_main_timer,         ONLY: time_step_len, lstart
    USE messy_main_grid_def_mem_bi, ONLY: nlev, nlon, ngl
    USE messy_main_data_bi,       ONLY: zairdens => rho_air_dry_3d &
                                      , slf
    USE messy_main_transform_bi,  ONLY: trp_gpdc_gpgl, M_SUM
    USE messy_main_constants_mem, ONLY: N_A, M_air
    USE messy_main_tracer_mem_bi, ONLY: LGTRSTR, CLTRSTR

    IMPLICIT NONE

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'tnudge_global_end'
    REAL(DP), PARAMETER         :: uconv = N_A * 1.0e3_DP/M_air
    INTEGER  :: i
    INTEGER  :: jc, jglon, jglat, jk
    INTEGER  :: kmin, kmax, klev
    INTEGER  :: kmin_tmp, kmax_tmp
    LOGICAL  :: llon, llat, llev, lrange
    INTEGER  :: klev_f            ! level in diagnostic flux
    REAL(DP) :: xten              ! additional nudging tendency
    REAL(DP) :: convert           ! unit conversion
    REAL(DP) :: ztmst_hard        ! actual time step for 'hard' nudging
    REAL(DP) :: zcoeff            ! nudging coeff.
    !
    REAL(DP), DIMENSION(:,:),   POINTER :: ibl_gl   => NULL()
    REAL(DP), DIMENSION(:,:),   POINTER :: itp_gl   => NULL()
    REAL(DP), DIMENSION(:,:,:), POINTER :: field_gl => NULL()
    REAL(DP), DIMENSION(:,:,:), POINTER :: flux_gl  => NULL()
    REAL(DP), DIMENSION(:,:,:), POINTER :: field_cgl => NULL()
    REAL(DP), DIMENSION(:,:,:), POINTER :: flux_cgl  => NULL()
    REAL(DP), DIMENSION(:,:,:), POINTER :: zairdens_gl => NULL()
    REAL(DP), DIMENSION(:,:),   POINTER :: slf_gl   => NULL() ! op_pj_20150428
    !
    REAL(DP), DIMENSION(:,:), POINTER :: gh, lh

    IF (.NOT. L_LG .AND. .NOT. L_CL) RETURN
    ! INIT
    ! SET TIME STEP LENGTH
    IF (lstart) THEN
       ztmst_hard = time_step_len * 2.0_dp
    ELSE
       ztmst_hard = time_step_len
    END IF

    ! ALLOCATE
    ALLOCATE(flux_gl(nlon,nlev,ngl))
    ALLOCATE(flux_cgl(nlon,nlev,ngl))

    ! GLOBALIZE LAYER INFO
    IF (lbl) CALL trp_gpdc_gpgl(1, ibl, ibl_gl)
    IF (ltp) CALL trp_gpdc_gpgl(1, itp, itp_gl)
    CALL trp_gpdc_gpgl(1, zairdens, zairdens_gl)
    ! GLOBALIZE SEA LAND FRACTION
    CALL trp_gpdc_gpgl(1, slf, slf_gl)

    FOR_ATTILA: IF (L_LG) THEN
       ! (1) FILL FLAG ARRAY
       nudging_loop_lg: DO i=1, NNUDGE  ! LOOP OVER NUDGING FIELDS

          ! ONLY THIS TRACER SET
          IF (TRIM(XNUDGE(i)%io%trset) /= TRIM(LGTRSTR)) CYCLE

          ! CHANNEL/OBJECT NOT AVAILABLE
          IF (.NOT.ASSOCIATED(XNUDGE(i)%field)) CYCLE

          ! TRACER NOT AVAILABLE / TROPOPAUSE NOT AVAILABLE / BOUNDARY LAYER
          ! HEIGHT NOT AVAILABLE
          IF (XNUDGE(i)%idt == 0) CYCLE

          ! GLOBALIZE
          CALL trp_gpdc_gpgl(1, XNUDGE(i)%field, field_gl)

          ! INIT FLUX
          flux_gl(:,:,:) = 0.0_dp

          ! SET NUDGING COEFF
          zcoeff = XNUDGE(i)%io%coeff
          IF (ABS(XNUDGE(i)%io%coeff - ztmst_hard) < EPS) THEN
             IF (lstart) &
                  zcoeff = zcoeff / 2.0_dp
          END IF

          cell_loop_lg: DO jc=1, NCELL

             jk    = NINT(iplev(jc))  ! level index
             jglat = NINT(iplat(jc))  ! latitude index
             jglon = NINT(iplon(jc))  ! longitude index

             llon = (plon(jc) >= XNUDGE(i)%io%lonmin) .AND. &
                  (plon(jc) <= XNUDGE(i)%io%lonmax)
             llat = (plat(jc) >= XNUDGE(i)%io%latmin) .AND. &
                  (plat(jc) <= XNUDGE(i)%io%latmax)

             SELECT CASE(XNUDGE(i)%io%kmin)
             CASE(-4)
                kmin = NINT(XNUDGE(i)%iviso1(jglon, jglat))
             CASE (-3)
                kmin = NINT(ibl_gl(jglon, jglat))
             CASE (-2)
                kmin = NINT(itp_gl(jglon, jglat))
             CASE(-1)
                kmin = 1
             CASE(0)
                kmin = nlev
             CASE DEFAULT
                kmin = XNUDGE(i)%io%kmin
             END SELECT

             SELECT CASE(XNUDGE(i)%io%kmax)
             CASE(-4)
                kmax = NINT(XNUDGE(i)%iviso2(jglon, jglat))
             CASE (-3)
                kmax = NINT(ibl_gl(jglon, jglat))
             CASE (-2)
                kmax = NINT(itp_gl(jglon, jglat))
             CASE(-1)
                kmax = 1
             CASE(0)
                kmax = nlev
             CASE DEFAULT
                kmax = XNUDGE(i)%io%kmax
             END SELECT

             kmin_tmp = MAX(MIN(kmin, kmax),1)
             kmax_tmp = MIN(MAX(kmin, kmax),nlev)
             kmin = kmin_tmp
             kmax = kmax_tmp

             llev = (jk >= kmin) .AND. (jk <= kmax)

             SELECT CASE(XNUDGE(i)%ndim)
             CASE(2)
                klev = 1
             CASE(3)
                klev = jk
             CASE DEFAULT
                CALL info_bi(' WARNING! CHANNEL OBJECT HAS RANK /= 2 OR 3 !' &
                     , substr)
             END SELECT

             lrange = ( &
                  ( field_gl(jglon, klev, jglat) > XNUDGE(i)%io%vr(1) ) &
                  .AND. &
                  ( field_gl(jglon, klev, jglat) < XNUDGE(i)%io%vr(2) ) )

             ! PERFORM NUDGING
             IF (llon.AND.llat.AND.llev.AND.lrange) THEN

                xten = -1._DP * &
                     (qxtm1_a(jc, XNUDGE(i)%idt) +                &
                     qxtte_a(jc, XNUDGE(i)%idt)*time_step_len  -  &
                     field_gl(jglon, klev, jglat)) / zcoeff

                ! only affect tendency of grid cell fraction above land
                IF (XNUDGE(i)%io%land == 1) &
                     xten = xten * slf_gl(jglon, jglat)
                ! only affect tendency of grid cell fraction above water
                IF (XNUDGE(i)%io%land == 2) &
                     xten = xten * (1._DP-slf_gl(jglon, jglat))
                IF ((jk == kmax) .AND. XNUDGE(i)%lfrac2) THEN
                   IF ((jk == kmin) .AND. XNUDGE(i)%lfrac1) THEN
                      xten = xten * MAX(XNUDGE(i)%fviso1(jglon, jglat) &
                           - XNUDGE(i)%fviso2(jglon, jglat), 0._DP)
                      ! only affect tendency of grid cell fraction
                      ! between interfaces
                   ELSE
                      xten = xten * (1._DP-XNUDGE(i)%fviso2(jglon, jglat))
                      ! only affect tendency of grid cell fraction above
                      ! lower interface
                   END IF
                ELSE
                   IF ((jk == kmin) .AND. XNUDGE(i)%lfrac1) &
                        xten = xten * XNUDGE(i)%fviso1(jglon, jglat)
                   ! only affect tendency of grid cell fraction below
                   ! upper interface
                END IF

                qxtte_a(jc, XNUDGE(i)%idt) = &
                     qxtte_a(jc, XNUDGE(i)%idt) + xten

                ! SAVE FLUX
                IF (XNUDGE(i)%io%lflux) THEN
                   ! SET LEVEL
                   IF (XNUDGE(i)%nfdim == 2) THEN
                      klev_f = 1
                      ! mol/mol/s -> molecules/(m^2 s)
                      convert = ( amcell / garea(jglon,jglat) ) * uconv
                   ELSE
                      klev_f = jk
                      ! mol/mol/s -> molecules/(m^3 s)
                      convert = zairdens_gl(jglon, jk, jglat) * uconv
                   END IF
                   !
                   flux_gl(jglon, klev_f, jglat) = xten * convert
                END IF
             END IF

          END DO cell_loop_lg

          ! DECOMPOSE FLUX ...
          IF (XNUDGE(i)%io%lflux) THEN
             SELECT CASE(XNUDGE(i)%nfdim)
             CASE(2)
                lh => XNUDGE(i)%flux(:,1,:)
                gh => flux_gl(:,1,:)
                CALL trp_gpdc_gpgl(-1, lh, gh, M_SUM)
             CASE(3)
                CALL trp_gpdc_gpgl(-1, XNUDGE(i)%flux, flux_gl, M_SUM)
             END SELECT
          END IF

          IF (ASSOCIATED(field_gl)) THEN
             DEALLOCATE(field_gl) ; NULLIFY(field_gl)
          ENDIF

       END DO nudging_loop_lg
    ENDIF FOR_ATTILA

    FOR_CLaMS: IF (L_CL) THEN
       ! (1) FILL FLAG ARRAY
       nudging_loop_cl: DO i=1, NNUDGE  ! LOOP OVER NUDGING FIELDS

          ! ONLY THIS TRACER SET
          IF (TRIM(XNUDGE(i)%io%trset) /= TRIM(CLTRSTR)) CYCLE

          ! CHANNEL/OBJECT NOT AVAILABLE
          IF (.NOT.ASSOCIATED(XNUDGE(i)%field)) CYCLE

          ! TRACER NOT AVAILABLE / TROPOPAUSE NOT AVAILABLE / BOUNDARY LAYER
          ! HEIGHT NOT AVAILABLE
          IF (XNUDGE(i)%idt == 0) CYCLE

          ! GLOBALIZE  ###############
          CALL trp_gpdc_gpgl(1, XNUDGE(i)%field, field_cgl)
          ! INIT FLUX
          flux_cgl(:,:,:) = 0.0_dp
          ! SET NUDGING COEFF
          zcoeff = XNUDGE(i)%io%coeff
          IF (ABS(XNUDGE(i)%io%coeff - ztmst_hard) < EPS) THEN
             IF (lstart) &
                  zcoeff = zcoeff / 2.0_dp
          END IF
          cell_loop_cl: DO jc=1, DNPARTS

             jk    = NINT(pos(jc,2))  ! level index
             jglat = NINT(pos(jc,3))  ! latitude index
             jglon = NINT(pos(jc,1))  ! longitude index

             llon = (LON(jc) >= XNUDGE(i)%io%lonmin) .AND. &
                  (LON(jc) <= XNUDGE(i)%io%lonmax)
             llat = (LAT(jc) >= XNUDGE(i)%io%latmin) .AND. &
                  (LAT(jc) <= XNUDGE(i)%io%latmax)

             if (jglat .ge. 1 .or. jglon .ge. 1) then

                SELECT CASE(XNUDGE(i)%io%kmin)
                CASE(-4)
                   kmin = NINT(XNUDGE(i)%iviso1(jglon, jglat))
                CASE (-3)
                   kmin = NINT(ibl_gl(jglon, jglat))
                CASE (-2)
                   kmin = NINT(itp_gl(jglon, jglat))
                CASE(-1)
                   kmin = 1
                CASE(0)
                   kmin = nlev
                CASE DEFAULT
                   kmin = XNUDGE(i)%io%kmin
                END SELECT

                SELECT CASE(XNUDGE(i)%io%kmax)

                CASE(-4)
                   kmax = NINT(XNUDGE(i)%iviso2(jglon, jglat))
                CASE (-3)
                   kmax = NINT(ibl_gl(jglon, jglat))
                CASE (-2)
                   kmax = NINT(itp_gl(jglon, jglat))
                CASE(-1)
                   kmax = 1
                CASE(0)
                   kmax = nlev
                CASE DEFAULT
                   kmax = XNUDGE(i)%io%kmax
                END SELECT

                kmin_tmp = MAX(MIN(kmin, kmax),1)
                kmax_tmp = MIN(MAX(kmin, kmax),nlev)
                kmin = kmin_tmp
                kmax = kmax_tmp

                llev = (jk >= kmin) .AND. (jk <= kmax)
                SELECT CASE(XNUDGE(i)%ndim)
                CASE(2)
                   klev = 1
                CASE(3)
                   klev = jk
                CASE DEFAULT
                   CALL info_bi(' WARNING! CHANNEL OBJECT HAS RANK /= 2 OR 3 !'&
                        , substr)
                END SELECT

                lrange = ( &
                     ( field_cgl(jglon, klev, jglat) > XNUDGE(i)%io%vr(1) ) &
                     .AND. &
                     ( field_cgl(jglon, klev, jglat) < XNUDGE(i)%io%vr(2) ) )

                ! PERFORM NUDGING
                IF (llon.AND.llat.AND.llev.AND.lrange) THEN
                   xten = -1._DP * &
                        (qxt_c(jc, XNUDGE(i)%idt) +                &
                        qxtte_c(jc, XNUDGE(i)%idt)*time_step_len  -  &
                        field_cgl(jglon, klev, jglat)) / zcoeff

                   ! only affect tendency of grid cell fraction above land
                   IF (XNUDGE(i)%io%land == 1) &
                        xten = xten * slf_gl(jglon, jglat)
                   ! only affect tendency of grid cell fraction above water
                   IF (XNUDGE(i)%io%land == 2) &
                        xten = xten * (1._DP-slf_gl(jglon, jglat))
                   IF ((jk == kmax) .AND. XNUDGE(i)%lfrac2) THEN
                      IF ((jk == kmin) .AND. XNUDGE(i)%lfrac1) THEN
                         xten = xten * MAX(XNUDGE(i)%fviso1(jglon, jglat) &
                              - XNUDGE(i)%fviso2(jglon, jglat), 0._DP)
                         ! only affect tendency of grid cell fraction
                         ! between interfaces
                      ELSE
                         xten = xten * (1._DP-XNUDGE(i)%fviso2(jglon, jglat))
                         ! only affect tendency of grid cell fraction above
                         ! lower interface
                      END IF

                   ELSE
                      IF ((jk == kmin) .AND. XNUDGE(i)%lfrac1) &
                           xten = xten * XNUDGE(i)%fviso1(jglon, jglat)
                      ! only affect tendency of grid cell fraction below
                      ! upper interface
                   END IF

                   !! qxt_c(jc, XNUDGE(i)%idt) = &
                   !!       qxt_c(jc, XNUDGE(i)%idt) + xten*time_step_len
                   qxtte_c(jc, XNUDGE(i)%idt) = &
                        qxtte_c(jc, XNUDGE(i)%idt) + xten

                   ! SAVE FLUX
                   IF (XNUDGE(i)%io%lflux) THEN
                      ! SET LEVEL
                      IF (XNUDGE(i)%nfdim == 2) THEN
                         klev_f = 1
                         ! mol/mol/s -> molecules/(m^2 s)
                         convert = (cmcell / garea(jglon,jglat) ) * uconv
                      ELSE
                         klev_f = jk
                         ! mol/mol/s -> molecules/(m^3 s)
                         convert = zairdens_gl(jglon, jk, jglat) * uconv
                      END IF
                      !
                      flux_cgl(jglon, klev_f, jglat) = xten * convert
                   END IF

                END IF
             ENDIF ! op_sb_20190802-
          END DO cell_loop_cl

          ! DECOMPOSE FLUX ...
          IF (XNUDGE(i)%io%lflux) THEN
             SELECT CASE(XNUDGE(i)%nfdim)
             CASE(2)
                lh => XNUDGE(i)%flux(:,1,:)
                gh => flux_cgl(:,1,:)
                CALL trp_gpdc_gpgl(-1, lh, gh, M_SUM)

             CASE(3)
                CALL trp_gpdc_gpgl(-1, XNUDGE(i)%flux, flux_cgl, M_SUM)

             END SELECT
          END IF

          IF (ASSOCIATED(field_cgl)) THEN
             DEALLOCATE(field_cgl) ; NULLIFY(field_cgl)
          ENDIF

       END DO nudging_loop_cl
    ENDIF FOR_CLaMS

    ! CLEAN UP
    IF (ASSOCIATED(ibl_gl)) THEN
       DEALLOCATE(ibl_gl) ; NULLIFY(ibl_gl)
    ENDIF
    IF (ASSOCIATED(itp_gl)) THEN
       DEALLOCATE(itp_gl) ; NULLIFY(itp_gl)
    ENDIF
    IF (ASSOCIATED(slf_gl)) THEN
       DEALLOCATE(slf_gl) ; NULLIFY(slf_gl)
    ENDIF

    DEALLOCATE(flux_gl)
    DEALLOCATE(flux_cgl)

    DEALLOCATE(zairdens_gl)
#endif

  END SUBROUTINE tnudge_global_end
  ! ---------------------------------------------------------------------

  ! ---------------------------------------------------------------------
  SUBROUTINE tnudge_free_memory

    IMPLICIT NONE

    INTRINSIC :: ASSOCIATED

    IF (ASSOCIATED(garea)) THEN
       DEALLOCATE(garea)
       NULLIFY(garea)
    END IF

    IF (ASSOCIATED(XNUDGE)) THEN
       DEALLOCATE(XNUDGE)
       NULLIFY(XNUDGE)
    END IF

  END SUBROUTINE tnudge_free_memory
  ! ---------------------------------------------------------------------

  ! ----------------------------------------------------------------------
  SUBROUTINE tnudge_read_nml_cpl(status, iou)

    ! tnudge MODULE ROUTINE (ECHAM-5 INTERFACE, PRIVATE)
    !
    ! read namelist for 'coupling' to ECHAM5
    !
    ! Author: Patrick Joeckel, MPICH, Dec 2003

    ! MESSy
    USE messy_main_tools, ONLY: read_nml_open, read_nml_check, read_nml_close

    IMPLICIT NONE

    ! I/O
    INTEGER, INTENT(OUT) :: status     ! error status
    INTEGER, INTENT(IN)  :: iou        ! I/O unit

    ! (LOCAL) NAMELIST VARIABLES
    CHARACTER(LEN=*), PARAMETER :: substr = 'tnudge_read_nml_cpl'

    NAMELIST /CPL/ tropopause, boundarylayer, TNUDGE

    ! LOCAL
    LOGICAL              :: lex      ! file exists ?
    INTEGER              :: fstat    ! file status

    status = 1

    ! INITIALIZE NAMELIST VARIABLES
    tropopause(1) = ' '
    tropopause(2) = ' '
    boundarylayer(1) = ' '
    boundarylayer(2) = ' '

    CALL read_nml_open(lex, substr, iou, 'CPL', modstr)
    IF (.not.lex) RETURN    ! <modstr>.nml does not exist

    READ(iou, NML=CPL, IOSTAT=fstat)
    CALL read_nml_check(fstat, substr, iou, 'CPL', modstr)
    IF (fstat /= 0) RETURN  ! error while reading namelist

    ! DIAGNOSE NAMELIST AND SET GLOBAL SWITCHES
    !
    ! CHECK NAMELIST
    WRITE(*,*) 'Using tropopause index ''',TRIM(tropopause(2))   &
         ,''' of channel ''',TRIM(tropopause(1)),''''

    WRITE(*,*) 'Using boundary layer index ''',TRIM(boundarylayer(2))   &
         ,''' of channel ''',TRIM(boundarylayer(1)),''''

    CALL read_nml_close(substr, iou, modstr)

    status = 0  ! no ERROR

  END SUBROUTINE tnudge_read_nml_cpl
  ! ----------------------------------------------------------------------

  ! --------------------------------------------------------------------
  SUBROUTINE tnudge_input_check(status, TNUDGE, XNUDGE, NNUDGE)

#ifdef ICON
    USE messy_main_channel_bi, ONLY: n_dom
#endif
    USE messy_main_tools,      ONLY: strcrack

    IMPLICIT NONE

    INTRINSIC :: ABS, TRIM
#ifdef ICON
    INTRINSIC :: INDEX
#endif
    ! I/O
    INTEGER,                                 INTENT(OUT)   :: status
    TYPE(T_IO_NUDGE), DIMENSION(:),          INTENT(INOUT) :: TNUDGE
    TYPE(T_XNUDGE),   DIMENSION(:), POINTER, INTENT(INOUT) :: XNUDGE
    INTEGER,                                 INTENT(OUT)   :: NNUDGE

    ! LOCAL
    CHARACTER(LEN=*),  PARAMETER             :: substr='tnudge_input_check'
    INTEGER                                  :: i, k, num
    CHARACTER(LEN=STRLEN_TRSET), DIMENSION(:), POINTER :: trs
    CHARACTER(LEN=2*STRLEN_TRSET+1)          :: trsname
    INTEGER                                  :: noe
#ifdef ICON
    INTEGER                                  :: idx, n
#endif
    INTEGER                                  :: INUDGE

    status = 0
    NULLIFY(trs)

    ! GET NUMBER OF ENTRIES
    DO i=1, NMAXNUDGE
       IF (TRIM(TNUDGE(i)%trset)   == '') CYCLE
       IF (TRIM(TNUDGE(i)%trname)  == '') CYCLE

       IF (TRIM(TNUDGE(i)%channel) == '') CYCLE
       IF (TRIM(TNUDGE(i)%object)  == '') CYCLE

       WRITE(*,*) '-----------------------------------------------'
       WRITE(*,*) 'TRACER SET(S)  : '//TRIM(TNUDGE(i)%trset)
       WRITE(*,*) 'TRACER         : '//TRIM(TNUDGE(i)%trname)

       WRITE(*,*) 'CHANNEL        : '//TRIM(TNUDGE(i)%channel)
       WRITE(*,*) 'OBJECT         : '//TRIM(TNUDGE(i)%object)

       ! %coeff

       WRITE(*,*) 'MIN. LATITUDE  : ',TNUDGE(i)%latmin
       IF ((TNUDGE(i)%latmin < -90.0_DP) .OR. &
            (TNUDGE(i)%latmin > 90.0_DP)) THEN
          WRITE(*,*) ' *** ERROR *** OUT OF RANGE [-90,90]'
          status = 1
          RETURN
       END IF

       WRITE(*,*) 'MAX. LATITUDE  : ',TNUDGE(i)%latmax
       IF ((TNUDGE(i)%latmax < -90.0_DP) .OR. &
            (TNUDGE(i)%latmax > 90.0_DP)) THEN
          WRITE(*,*) ' *** ERROR *** OUT OF RANGE [-90,90]'
          status = 2
          RETURN
       END IF

       ! %kmin
       ! %kmax

       WRITE(*,*) 'MIN. LONGITUDE : ',TNUDGE(i)%lonmin
       IF ((TNUDGE(i)%lonmin < -360.0_DP) .OR. &
            (TNUDGE(i)%lonmin > 360.0_DP)) THEN
          WRITE(*,*) ' *** ERROR *** OUT OF RANGE [-360,360]'
          status = 3
          RETURN
       END IF

       WRITE(*,*) 'MAX. LONGITUDE : ',TNUDGE(i)%lonmax
       IF ((TNUDGE(i)%lonmax < -360.0_DP) .OR. &
            (TNUDGE(i)%lonmax > 360.0_DP)) THEN
          WRITE(*,*) ' *** ERROR *** OUT OF RANGE [-360,360]'
          status = 4
          RETURN
       END IF

       IF (TNUDGE(i)%lflux) THEN
          WRITE(*,*) 'DIAGNOSTIC FLUX: ON'
       ELSE
          WRITE(*,*) 'DIAGNOSTIC FLUX: OFF'
       ENDIF

       IF (TNUDGE(i)%kmin == -4) THEN
          IF (TRIM(TNUDGE(i)%cviso1)  == '') THEN
             WRITE(*,*) ' *** ERROR *** MIN. LEVEL = -4,'//&
                  &' BUT NO ISO-SURFACE CHANNEL NAME IS PROVIDED'
             status = 5
             RETURN
          END IF
          IF (TRIM(TNUDGE(i)%oviso1)  == '') THEN
             WRITE(*,*) ' *** ERROR *** MIN. LEVEL = -4,'//&
                  &' BUT NO ISO-SURFACE OBJECT NAME IS PROVIDED'
             status = 6
             RETURN
          END IF
          WRITE(*,*) 'MIN ISO-SURFACE: ',&
               TRIM(TNUDGE(i)%oviso1),' of channel ',TRIM(TNUDGE(i)%cviso1)
       END IF

       IF (TNUDGE(i)%kmax == -4) THEN
          IF (TRIM(TNUDGE(i)%cviso2)  == '') THEN
             WRITE(*,*) ' *** ERROR *** MAX. LEVEL = -4,'//&
                  &' BUT NO ISO-SURFACE CHANNEL NAME IS PROVIDED'
             status = 7
             RETURN
          END IF
          IF (TRIM(TNUDGE(i)%oviso2)  == '') THEN
             WRITE(*,*) ' *** ERROR *** MAX. LEVEL = -4,'//&
                  &' BUT NO ISO-SURFACE OBJECT NAME IS PROVIDED'
             status = 8
             RETURN
          END IF
          WRITE(*,*) 'MAX ISO-SURFACE: ',&
               TRIM(TNUDGE(i)%oviso2),' of channel ',TRIM(TNUDGE(i)%cviso2)
       END IF

       WRITE(*,*) 'LAND FLAG      : ',TNUDGE(i)%land
       IF ((TNUDGE(i)%land < 0) .OR. (TNUDGE(i)%land > 2)) THEN
          WRITE(*,*) ' *** ERROR *** LAND FLAG out of range {0, 1, 2}: ',&
               TNUDGE(i)%land
          status = 9
          RETURN
       END IF

       WRITE(*,*) 'VALID RANGE    : ',TNUDGE(i)%vr(:)

       CALL strcrack(TRIM(TNUDGE(i)%trset), ';', trs, num)
#if defined(ICON)
!! qqq do this only for GPTRSTR???
       ! check for specific patches
       noe = 0
       DO k=1,  num
          idx = INDEX(trsname, '_D')
          IF (idx == 0) THEN
             ! if patch is unspecified, use it as wildcard for all patches;
             ! -> add n_dom entries
             noe = noe + n_dom
          ELSE
             ! only one specific patch
             noe = noe + 1
          END IF
       END DO
#else
       noe = num
#endif
       NNUDGE = NNUDGE + noe
       IF (ASSOCIATED(trs)) THEN
          DEALLOCATE(trs); NULLIFY(trs)
       END IF

       WRITE(*,*) 'NO. OF ELEMENTS: ',noe
       WRITE(*,*) '-----------------------------------------------'
    END DO

    WRITE(*,*) '-----------------------------------------------'
    WRITE(*,*) 'TOTAL NUMBER OF ELEMENTS: ',NNUDGE
    WRITE(*,*) '-----------------------------------------------'

    ALLOCATE(XNUDGE(NNUDGE))

    ! COPY DATA
    INUDGE = 0
    DO i=1, NMAXNUDGE

       IF (TRIM(TNUDGE(i)%trset)   == '') CYCLE
       IF (TRIM(TNUDGE(i)%trname)  == '') CYCLE

       IF (TRIM(TNUDGE(i)%channel) == '') CYCLE
       IF (TRIM(TNUDGE(i)%object)  == '') CYCLE

       ! ALL POTENTIAL ERRORS HAVE BEEN ALREADY TESTED ABOVE

       ! get number of different tracer sets
       CALL strcrack(TRIM(TNUDGE(i)%trset), ';', trs, num)

#if defined(ICON)
!! qqq do this only for GPTRSTR???
       ! check for specific patches
       DO k=1,  num
          trsname = trs(k)
          idx = INDEX(trsname, '_D')
          IF (idx == 0) THEN
             ! if patch is unspecified, use it as wildcard for all patches;
             ! -> add n_dom entries
             DO n=1, n_dom
                INUDGE = INUDGE+1
                IF (INUDGE > NNUDGE) &
                     CALL error_bi('nudging counter out of bounds',substr)
                CALL copy2x(status,XNUDGE(INUDGE),TNUDGE(i),trsname,n)
                IF (status /= 0) RETURN
             END DO
          ELSE
             INUDGE = INUDGE+1
             IF (INUDGE > NNUDGE) &
                  CALL error_bi('nudging counter out of bounds',substr)
             CALL copy2x(status, XNUDGE(INUDGE),TNUDGE(i),trsname)
             IF (status /= 0) RETURN
          END IF
       END DO
#else
       DO k=1, num
          trsname = trs(k)
          INUDGE = INUDGE+1
          IF (INUDGE > NNUDGE) &
               CALL error_bi('nudging counter out of bounds',substr)
          CALL copy2x(status, XNUDGE(INUDGE),TNUDGE(i),trsname)
          IF (status /= 0) RETURN
       END DO
#endif
       IF (ASSOCIATED(trs)) THEN
          DEALLOCATE(trs); NULLIFY(trs)
       END IF

    END DO

  CONTAINS

    SUBROUTINE copy2x(status, X, T, trs, dom)

      USE messy_main_data_bi,     ONLY: l2tls
      USE messy_main_timer,       ONLY: timer_get_time_step_len, lstart
      USE messy_main_channel_mem, ONLY: dom_curid
      USE messy_main_tools,       ONLY: int2str

      IMPLICIT NONE
      INTRINSIC :: MINVAL, MAXVAL

      ! I/O
      INTEGER,                         INTENT(OUT) :: status
      TYPE(T_XNUDGE),                  INTENT(OUT) :: X
      TYPE(T_IO_NUDGE),                INTENT(IN)  :: T
      CHARACTER(LEN=2*STRLEN_TRSET+1), INTENT(IN)  :: trs
      INTEGER, OPTIONAL,               INTENT(IN)  :: dom

      ! LOCAL
      INTEGER          :: zdom
      REAL(dp)         :: ztmst_hard
      CHARACTER(LEN=2) :: dstr

      status = 0

      ! copy everything
      X%io = T

      ! special
      X%io%vr(1) = MINVAL(T%vr(:))
      X%io%vr(2) = MAXVAL(T%vr(:))

      ! set domain
      IF (PRESENT(dom)) THEN
         CALL int2str(dstr, dom)
         X%io%trset   = TRIM(trs)//'_D'//dstr
         zdom = dom
      ELSE
         X%io%trset   = TRIM(trs)
         zdom = dom_curid
      END IF

      ! set nudging coefficient
      WRITE(*,*) '-----------------------------------------------'
      WRITE(*,*) 'TRACER SET     : '//TRIM(X%io%trset)
      WRITE(*,*) 'TRACER         : '//TRIM(X%io%trname)

      CALL timer_get_time_step_len(ztmst_hard, zdom)
!qqq
      IF (lstart .AND. .NOT. l2tls) ztmst_hard = 2._dp * ztmst_hard

      IF (T%coeff < 0.0_DP) THEN
         WRITE(*,*) 'NUDGING COEFF. : <0 => hard nudging'
         X%io%coeff = ztmst_hard
      ELSE
         X%io%coeff = T%coeff
      END IF
      WRITE(*,*) 'NUDGING COEFF. : ',X%io%coeff

      IF (ABS(X%io%coeff) <= EPS) THEN
         WRITE(*,*) ' *** ERROR *** NUDGING COEFF. IS ZERO'
         status = 10
         RETURN
      END IF
      IF (ABS(X%io%coeff) < ztmst_hard) THEN
         WRITE(*,*)  &
              ' *** ERROR *** NUDGING COEFF. IS LOWER THAN MODEL TIME STEP'
         status = 20
         RETURN
      END IF
      WRITE(*,*) '-----------------------------------------------'

      status = 0

    END SUBROUTINE copy2x

  END SUBROUTINE tnudge_input_check
  ! ------------------------------------------------------------------------

!***********************************************************************
END MODULE messy_tnudge_si
!***********************************************************************

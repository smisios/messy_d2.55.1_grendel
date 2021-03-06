#include "messy_main_ppd_bi.inc"

!*************************************************************************
MODULE messy_offemis_si
!*************************************************************************

  ! Author: Patrick Joeckel, MPICH, Mar 2004
  !         Astrid Kerkweg, UNI-MZ, Oct 2009 Umstellung auf OFFEMIS
  ! TODO:
  ! - output lagrangian emission flux as lg-channel object ?
  !   (local pointer as alternative ?)

  ! BMIL/MESSy
  USE messy_main_blather_bi,    ONLY: start_message_bi, end_message_bi &
                                    , error_bi, warning_bi
  USE messy_main_channel_bi,    ONLY: n_dom
#ifdef MESSYTENDENCY
 !tendency budget
 USE messy_main_tendency_bi,    ONLY: mtend_get_handle,       &
                                      mtend_add_l,            &
                                      mtend_register,         &    
                                      mtend_id_tracer
#endif
  ! MESSy
  USE messy_main_constants_mem, ONLY: STRLEN_MEDIUM, STRLEN_ULONG, STRLEN_XLONG
  USE messy_main_channel,       ONLY: STRLEN_CHANNEL, STRLEN_OBJECT
  USE messy_offemis

  IMPLICIT NONE
  INTRINSIC :: TRIM, ADJUSTL
  PRIVATE

  TYPE T_EMIS_IN
     CHARACTER(LEN=10*STRLEN_XLONG) :: tr_data = ''   ! tracer name_subname,scaling, domains
     CHARACTER(LEN=STRLEN_CHANNEL):: ch_name     = '' ! channel name
     CHARACTER(LEN=STRLEN_OBJECT) :: ch_object   = '' ! channel object name
     CHARACTER(LEN=STRLEN_ULONG)  :: action_str  = '' ! emission method
  END TYPE T_EMIS_IN

  INTEGER, PARAMETER                           :: N_MAX_EMIS = 1000
  TYPE(T_EMIS_IN), DIMENSION(N_MAX_EMIS), SAVE :: EMIS_IN

  TYPE T_EMIS_SET
     CHARACTER(LEN=STRLEN_MEDIUM), DIMENSION(:),POINTER :: tr_basename => NULL()
     CHARACTER(LEN=STRLEN_MEDIUM), DIMENSION(:),POINTER :: tr_subname => NULL()
     REAL(DP), DIMENSION(:), POINTER                    :: tr_factor =>NULL()  
     INTEGER                                            :: ntr = 0
     INTEGER                                            :: domain_idx = 0

     CHARACTER(LEN=STRLEN_CHANNEL):: ch_name     = '' ! channel name
     CHARACTER(LEN=STRLEN_OBJECT) :: ch_object   = '' ! channel object name
     CHARACTER(LEN=STRLEN_ULONG)  :: action_str  = '' ! emission method
     CHARACTER(LEN=STRLEN_ULONG)  :: heights     = '' ! emission method
     REAL(DP), DIMENSION(:), POINTER :: z => NULL()! emission heights ...
                                               ! ... for multi level emissions
     ! EMISSION TYPE (1: 2D surface, 2: 3D volume, 3: Nx2D multi level)
     INTEGER                     :: etype = 0
     ! POINTER TO MEMORY (INPUT FIELD)
     REAL(DP), DIMENSION(:,:,:), POINTER :: raw => NULL()
     ! POINTER TO CHANNEL OBJECT FOR INPUT DATA
     REAL(DP), DIMENSION(:,:,:), POINTER :: field => NULL()
     ! Vertical Index at which Emission is applied
     REAL(DP), DIMENSION(:,:,:), POINTER :: index => NULL()
     ! LOGICAL identifying whether imported emission is used a second time
     LOGICAL                             :: obj_exists = .FALSE.
     !
     ! GRIDPOINT TRACERS ---------------------------------------------------
     INTEGER  :: ngpem = 0  ! EMISSION METHOD
     ! -> tracer indices
     INTEGER, DIMENSION(:), POINTER  :: idt_gp =>NULL() ! TRACER idt

     ! ---------------------------------------------------------------------
     ! LAGRANGIAN TRACERS (ATTILA) -----------------------------------------
     INTEGER  :: nlgem = 0  ! EMISSION METHOD
     ! -> tracer indices
     INTEGER, DIMENSION(:),      POINTER :: idt_lg  => NULL() ! TRACER idt
     ! -> 'rest' of flux
     REAL(DP), DIMENSION(:,:,:), POINTER :: rest    => NULL()
     ! -> tracer tendency (OPTIONAL)
     REAL(dp), DIMENSION(:),     POINTER :: xtte_lg => NULL()
     ! ---------------------------------------------------------------------

     ! LAGRANGIAN TRACERS (CLaMS) ------------------------------------------
     INTEGER  :: nclem = 0  ! EMISSION METHOD
     INTEGER, DIMENSION(:), POINTER  :: idt_cl =>NULL() ! TRACER idt
     ! -> 'rest' of flux
     REAL(DP), DIMENSION(:,:,:), POINTER :: rest_cl    => NULL()
     ! -> tracer tendency (OPTIONAL)
     REAL(dp), DIMENSION(:), POINTER :: xtte_cl => NULL()
     ! ---------------------------------------------------------------------
  END TYPE T_EMIS_SET
  !
  TYPE(T_EMIS_SET), DIMENSION(:), POINTER :: EMIS_SET  => NULL()
  INTEGER, SAVE                           :: NEMIS_SET = 0

  ! GLOBAL PARAMETERS (CPL-NAMELIST)
  LOGICAL :: L_LG       = .false.  ! emissions for Lagrangian ATTILA tracers
  LOGICAL :: L_GP       = .false.  ! emissions for Gridpoint tracers
  LOGICAL :: L_CL       = .false.  ! emissions for Lagrangian CLaMS tracers
  LOGICAL :: l_lg_tend  = .false.  ! Lagrangian tracer tendency in channel
  LOGICAL :: l_cl_tend  = .false.  ! Lagrangian tracer tendency in channel
#ifdef MESSYTENDENCY
  INTEGER :: my_handle
#endif

  PUBLIC :: offemis_initialize
  PUBLIC :: offemis_new_tracer
  PUBLIC :: offemis_init_coupling
  PUBLIC :: offemis_vdiff
  PUBLIC :: offemis_global_end
  PUBLIC :: offemis_free_memory

  !PRIVATE :: offemis_read_nml_cpl
  !PRIVATE :: offemis_vdiff_gp
  !PRIVATE :: offemis_global_end_lg
  !PRIVATE :: offemis_global_end_cl

CONTAINS

! ------------------------------------------------------------------------
  SUBROUTINE  offemis_initialize

    ! offemis MODULE ROUTINE (ECHAM-5 INTERFACE)
    !
    ! INITIALIZATION OF GLOBAL VARIABLES FROM NAMELIST
    ! IN PARALLEL ENVIRONMENT
    ! INITIALIZATION OF offemis SPECIFIC EVENTS FROM NAMELIST
    ! IN PARALLEL ENVIRONMENT
    !

    ! BMIL/MESSy
    USE messy_main_mpi_bi,       ONLY: p_parallel_io, p_io, p_bcast
    USE messy_main_tools,        ONLY: find_next_free_unit &
                                     , domains_from_string

    IMPLICIT NONE

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'offemis_initialize'
    INTEGER                     :: iou    ! I/O unit
    INTEGER                     :: status ! error status
    INTEGER                     :: i, n
    INTEGER                     :: j, ntr
    REAL(DP), POINTER, DIMENSION(:)                         :: factor=>NULL()
    CHARACTER(LEN=2*STRLEN_MEDIUM+1), DIMENSION(:), POINTER :: name=>NULL()
    CHARACTER(LEN=2*STRLEN_MEDIUM+1), DIMENSION(:), POINTER :: subname=>NULL()
    INTEGER, PARAMETER                                      :: MAXPSTRLEN = 500
    CHARACTER(LEN=STRLEN_ULONG)                             :: tr_data = '  '
    INTEGER,           DIMENSION(:), POINTER :: domnum  => NULL()
    INTEGER                                  :: nd, num

    CALL start_message_bi(modstr, 'GLOBAL SETUP',substr)

    ! READ CTRL-namelist

    ! READ CPL-namelist
    IF (p_parallel_io) THEN
       iou = find_next_free_unit(100,200)
       CALL offemis_read_nml_cpl(status, iou)
       IF (status /= 0) CALL error_bi(' ',substr)
    END IF
    ! BROADCAST RESULTS
    CALL p_bcast(L_GP, p_io)
    CALL p_bcast(L_LG, p_io)
    CALL p_bcast(L_CL, p_io)
    CALL p_bcast(l_lg_tend, p_io)
    CALL p_bcast(l_cl_tend, p_io)

    ! MEMORY SPACE FOR EMIS-HANDLING
    ! a) COUNT NUMBER OF EMISSIONS
    IF (p_parallel_io) THEN
       n = 0
       DO i=1, N_MAX_EMIS

          ! CYCLE FOR EMPTY entry
          IF (TRIM(EMIS_IN(i)%tr_data) == '') CYCLE

          IF (TRIM(EMIS_IN(i)%ch_name) == '' .OR. &
               TRIM(EMIS_IN(i)%ch_object) == '') &
               CALL error_bi( &
               'empty channel or object for '//TRIM(EMIS_IN(i)%tr_data), substr)
          IF (TRIM(EMIS_IN(i)%action_str) == '') &
               CALL error_bi( &
               'empty action string for '//TRIM(EMIS_IN(i)%tr_data), substr)

          CALL domains_from_string(status,TRIM(EMIS_IN(i)%tr_data), n_dom, num)
          IF (status /= 0) CALL error_bi('error in namelist parsing (EMIS_IN)' &
               , substr)

          n=n+num
       ENDDO
       WRITE(*,*) '--------------------------------------------'
       WRITE(*,*) 'NUMBER OF EMISSION-SETS :', n
       WRITE(*,*) '--------------------------------------------'
    ENDIF
    CALL p_bcast(n, p_io)

    NEMIS_SET=n
    ! b) ALLOCATE MEMORY SPACE
    ALLOCATE(EMIS_SET(NEMIS_SET))

    IF (p_parallel_io) THEN
       n = 0
       DO i=1, N_MAX_EMIS
          IF (TRIM(EMIS_IN(i)%tr_data) == '') CYCLE

          IF (TRIM(EMIS_IN(i)%ch_name) == '' .OR. &
               TRIM(EMIS_IN(i)%ch_object) == '') CYCLE
          IF (TRIM(EMIS_IN(i)%action_str) == '') CYCLE

          CALL domains_from_string(status, TRIM(EMIS_IN(i)%tr_data), n_dom &
               , num &
                ,tr_data, dnums=domnum)
          IF (status /= 0) CALL error_bi('error in namelist parsing (EMIS_IN)' &
               , substr)

          ! using routine parse_datastr to crack 
          ! "tracer[_subname][,scaling];..."
          CALL parse_datatstr(status, MAXPSTRLEN, tr_data, name, &
                                      subname, factor, ntr)
          IF (status /= 0) CALL error_bi('parse_datatstr reported an error' &
               ,substr)

          domain_loop1: DO nd = 1, SIZE(domnum)
             n=n+1

             EMIS_SET(n)%ntr = ntr
             EMIS_SET(n)%domain_idx = domnum(nd)

             ALLOCATE(EMIS_SET(n)%tr_basename(EMIS_SET(n)%ntr), &
                  EMIS_SET(n)%tr_subname( EMIS_SET(n)%ntr), &
                  EMIS_SET(n)%tr_factor(EMIS_SET(n)%ntr ),  &
                  EMIS_set(n)%idt_gp(EMIS_SET(n)%ntr) )

             ! initial / default values
             EMIS_SET(n)%tr_basename(:) = ''
             EMIS_SET(n)%tr_subname(:)  = ''
             EMIS_SET(n)%tr_factor(:)   = 1.0_dp
             EMIS_set(n)%idt_gp(:)      = 0
             IF (L_LG) THEN
                ALLOCATE(EMIS_set(n)%idt_lg(EMIS_SET(n)%ntr)) 
                EMIS_set(n)%idt_lg(:)      = 0
             END IF
             IF (L_CL) THEN
                ALLOCATE(EMIS_set(n)%idt_cl(EMIS_SET(n)%ntr)) 
                EMIS_set(n)%idt_cl(:)      = 0
             END IF
             ! loop over tracers
             do j=1,  EMIS_SET(n)%ntr
                EMIS_SET(n)%tr_basename(j) = trim(name(j))
                EMIS_SET(n)%tr_subname(j)  = trim(subname(j))
                EMIS_SET(n)%tr_factor(j)   = factor(j)
             end do

             EMIS_SET(n)%ch_name     = EMIS_IN(i)%ch_name
             EMIS_SET(n)%ch_object   = EMIS_IN(i)%ch_object
             EMIS_SET(n)%action_str  = EMIS_IN(i)%action_str
          END DO domain_loop1

          IF (ASSOCIATED(name)) THEN
             DEALLOCATE(name) ; NULLIFY(name)
          ENDIF
          IF (ASSOCIATED(subname)) THEN
             DEALLOCATE(subname) ; NULLIFY(subname)
          ENDIF
          IF (ASSOCIATED(factor)) THEN
             DEALLOCATE(factor) ; NULLIFY(factor)
          ENDIF

       ENDDO
    ENDIF

    DO n=1,NEMIS_SET

       CALL p_bcast(EMIS_SET(n)%ntr,        p_io)
       CALL p_bcast(EMIS_SET(n)%domain_idx, p_io)
       ! allocate memory on non-p_io PEs
       IF (.NOT. p_parallel_io) THEN
          ALLOCATE(EMIS_SET(n)%tr_basename(EMIS_SET(n)%ntr), &
                   EMIS_SET(n)%tr_subname( EMIS_SET(n)%ntr), &
                   EMIS_SET(n)%tr_factor(EMIS_SET(n)%ntr ),  &
                   EMIS_set(n)%idt_gp(EMIS_SET(n)%ntr))
          IF (L_LG) ALLOCATE(EMIS_set(n)%idt_lg(EMIS_SET(n)%ntr))
          IF (L_CL) ALLOCATE(EMIS_set(n)%idt_cl(EMIS_SET(n)%ntr))     
       ENDIF

       DO j=1, EMIS_SET(n)%ntr 
          CALL p_bcast(EMIS_SET(n)%tr_basename(j) , p_io)
          CALL p_bcast(EMIS_SET(n)%tr_subname(j)  , p_io)
          CALL p_bcast(EMIS_SET(n)%tr_factor(j)   , p_io)
          EMIS_set(n)%idt_gp(:) = 0
          IF (L_LG) EMIS_set(n)%idt_lg(:) = 0
          IF (L_CL) EMIS_set(n)%idt_cl(:) = 0
       END DO
       
       ! op_mm_20140116-
       CALL p_bcast(EMIS_SET(n)%ch_name     , p_io)
       CALL p_bcast(EMIS_SET(n)%ch_object   , p_io)
       CALL p_bcast(EMIS_SET(n)%action_str  , p_io)
    ENDDO

    ! PARSE ACTION STRING
    DO i=1, NEMIS_SET
       CALL parse_str(status, STRLEN_ULONG, EMIS_SET(i)%action_str  &
            , EMIS_SET(i)%nclem, EMIS_SET(i)%nlgem, EMIS_SET(i)%ngpem )
#if defined(COSMO) || defined(ICON)
       IF (EMIS_SET(i)%ngpem >1) THEN
          CALL warning_bi('COSMO does not allow GP=2',substr)
          CALL warning_bi('switch back to GP=1'      ,substr)
          EMIS_SET(i)%ngpem = 1
       ENDIF
#endif
       !
       ! OUTPUT RESULT
       IF (p_parallel_io) THEN
          WRITE(*,*) 'EMIS-SET NUMBER         : ', i
          WRITE(*,*) 'EMIS-SET No. of TRACERS : ', EMIS_SET(i)%ntr
#ifdef ICON
          WRITE(*,*) 'EMIS-SET on DOMAIN      : ', EMIS_SET(i)%domain_idx
#endif
          WRITE(*,*) 'EMIS-SET EM-METHOD (GP) : ', EMIS_SET(i)%ngpem
          WRITE(*,*) 'EMIS-SET EM-METHOD (LG) : ', EMIS_SET(i)%nlgem
          WRITE(*,*) 'EMIS-SET EM-METHOD (CL) : ', EMIS_SET(i)%nclem
      
          ! op_mm_20131212+ (added do loop)
          do j=1, EMIS_SET(i)%ntr
             IF (TRIM(EMIS_SET(i)%tr_subname(j)) == '') THEN
                WRITE(*,*) 'EMIS-SET TRACER         : ', j, &
                     TRIM(EMIS_SET(i)%tr_basename(j))
             ELSE
                WRITE(*,*) 'EMIS-SET TRACER         : ', j, &
                     TRIM(EMIS_SET(i)%tr_basename(j))//'_'&
                     &//TRIM(EMIS_SET(i)%tr_subname(j))
             END IF
             WRITE(*,*) 'EMIS-SET SCALING        : ', j, &
                  EMIS_SET(i)%tr_factor(j)
          end do
          ! op_mm_20131212-

          WRITE(*,*) 'EMIS-SET CHANNEL        : ', TRIM(EMIS_SET(i)%ch_name)
          WRITE(*,*) 'EMIS-SET CH-OBJECT      : ', TRIM(EMIS_SET(i)%ch_object)
          WRITE(*,*) '--------------------------------------------'
       END IF
       !

       SELECT CASE (status)
       CASE(0)
          ! OK ?
       CASE(1)
          CALL error_bi('ERROR IN ACTION-STRING',substr)
       CASE(2)
          CALL error_bi('EMPTY SPECIFICATION IN ACTION STRING',substr)
       CASE(3)
          CALL error_bi('UNKNOWN SPECIFIER IN ACTION-STRING',substr)
       CASE(4)
          CALL error_bi('ERROR IN LG-SPECIFIER IN ACTION-STRING',substr)
       CASE(5)
          CALL error_bi('ERROR IN GP-SPECIFIER IN ACTION-STRING',substr)
       CASE(6)
          CALL error_bi('ERROR IN LG-SPECIFIER IN ACTION-STRING [0,4]',substr)
       CASE(7)
          CALL error_bi('ERROR IN GP-SPECIFIER IN ACTION-STRING [0,2]',substr)
       CASE(44)
          CALL error_bi('ERROR IN CL-SPECIFIER IN ACTION-STRING',substr)
       CASE(66)
          CALL error_bi('ERROR IN CL-SPECIFIER IN ACTION-STRING [0,4]',substr)
       CASE DEFAULT
          CALL error_bi('ERROR IN ACTION-STRING',substr)
       END SELECT
    END DO

    CALL end_message_bi(modstr, 'GLOBAL SETUP',substr)

  END SUBROUTINE offemis_initialize
! ------------------------------------------------------------------------

! ------------------------------------------------------------------------  
  SUBROUTINE offemis_new_tracer
    
    IMPLICIT NONE

#ifdef MESSYTENDENCY
    ! get handle for tendency treatment
    my_handle = mtend_get_handle(modstr)
    CALL mtend_register(my_handle, mtend_id_tracer) ! op_pj_20160721 see below
#endif

  END SUBROUTINE offemis_new_tracer
! ------------------------------------------------------------------------

! ------------------------------------------------------------------------
  SUBROUTINE offemis_init_coupling

    ! BMIL/MESSy
    USE messy_main_mpi_bi,           ONLY: p_parallel_io
    USE messy_main_tracer_mem_bi,    ONLY: GPTRSTR, LGTRSTR, CLTRSTR 
    USE messy_main_channel_error_bi, ONLY: channel_halt
#if defined(ECHAM5)
    USE messy_main_channel_bi,       ONLY: LG_ATTILA, REPR_LG_CLAMS
#endif
    USE messy_main_channel_bi,       ONLY: GP_3D_MID 
    ! MESSy
    USE messy_main_channel,          ONLY: new_channel, new_channel_object   &
                                         , new_attribute, get_channel_object &
                                         , get_channel_object_info           &
                                         , get_attribute                     &
                                         , get_channel_object_dimvalue
    USE messy_main_channel_mem,      ONLY: dom_curid
    USE messy_main_tracer,           ONLY: get_tracer


    IMPLICIT NONE

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER  :: substr='offemis_init_coupling'
    INTEGER                      :: status
    INTEGER                      :: i
    INTEGER                      :: j        ! op_mm_20140116
    INTEGER                      :: mlev     ! number of levels in emission
                                             ! channel object (2D or 3D)

    CHARACTER(LEN=STRLEN_ULONG)  :: unit     = ''
    CHARACTER(LEN=STRLEN_ULONG)  :: unit_ori = '' ! original object units

    INTEGER                      :: reprid   ! REPRESENTATION ID
    CHARACTER(LEN=4)             :: e_axis   ! axis of channel objects
    INTEGER                      :: ind
    CHARACTER(LEN=2 * STRLEN_MEDIUM + 1):: outname   = '' ! channel object name
    CHARACTER(LEN=STRLEN_ULONG)  :: heightaxis = ' '

     CALL start_message_bi(modstr, 'COUPLING SETUP',substr)

    ! define new channel
    CALL new_channel(status, modstr)
    CALL channel_halt(substr, status)

    offemis_set_loop: DO i=1, NEMIS_SET

       IF (EMIS_SET(i)%domain_idx /= dom_curid) CYCLE

       IF (p_parallel_io) THEN
          WRITE(*,*) '========================================================'
          WRITE(*,*) 'EMISSION-SET       : ',i
       END IF

       CALL get_channel_object(status, TRIM(EMIS_SET(i)%ch_name) &
            , TRIM(EMIS_SET(i)%ch_object), p3=EMIS_SET(i)%raw)
       IF (status /= 0 ) THEN
          CALL error_bi( 'ERROR: '//TRIM(EMIS_SET(i)%ch_name)//' '// &
               TRIM(EMIS_SET(i)%ch_object)//' not available!', substr)
       END IF

       CALL get_channel_object_info(status, TRIM(EMIS_SET(i)%ch_name) &
            , TRIM(EMIS_SET(i)%ch_object), reprid=reprid, axis=e_axis)
       
       mlev = -1
       ! -1 indicates surface emissions in contrast to 
       !  1 indicates Nx2D emissions with elements of Z equals 1  (1 Level! )   
       DO ind=1,4
          IF (e_axis(ind:ind) == 'Z') THEN
             mlev = SIZE(EMIS_SET(i)%raw, ind)
             EXIT
          ENDIF
       END DO

       CALL get_attribute(status, TRIM(EMIS_SET(i)%ch_name) &
            , TRIM(EMIS_SET(i)%ch_object), 'units' &
            , c=unit_ori)
       IF (status /= 0) unit_ori = ''

       level_analysis: IF (mlev > 0) THEN 
          IF (reprid /= GP_3D_MID) THEN
             ! Nx2D (multi level) emission
             EMIS_SET(i)%etype = 3  ! Nx2D (MULTI LEVEL EMISSION)
             ! mlev = number of emission levels in channel object
             ! get pressure levels via dimension variable
             CALL get_channel_object_dimvalue(status       &
                  ,  TRIM(EMIS_SET(i)%ch_name)             &
                  , TRIM(EMIS_SET(i)%ch_object)            &
                  , data=EMIS_SET(i)%z, axis ='Z', levtype=heightaxis)

             IF (TRIM(ADJUSTL(heightaxis)) /= 'heights') CALL error_bi ( &
                  ' vertical axis is not a height axis', substr)

             ! CHECK CONSISTANCY (number of emission levels in namelist)
             IF (mlev /= SIZE(EMIS_SET(i)%z)) THEN
                IF (p_parallel_io) THEN
                   WRITE(*,*) substr//' - ERROR:'
                   WRITE(*,*) 'NUMBER OF EMISSION LEVELS IN ATTRIBUTE  (',&
                        SIZE(EMIS_SET(i)%z),')'
                   WRITE(*,*) 'INCONSISTENT WITH DIMENSION LENGTH OF'
                   WRITE(*,*) 'CHANNEL OBJECT (',mlev,') !!!'
                   CALL error_bi(' ',substr)
                END IF
             END IF
          ELSE
             !
             ! IF THE NUMBER OF LEVELS IS > 1, VOLUME EMISSION (3D) IS ASSUMED,
             ! SINCE THE VERTICAL AXIS INFORMATION HAS BEEN PASSED
             ! THROUGH NCREGRID (HERE RAW-DATA IMPORT), I.E., THE VERTICAL
             ! AXIS IS SPECIFIED IN THE RESPECTIVE &REGRID NAMELIST.
             !
             EMIS_SET(i)%etype = 2       ! 3D (VOLUME EMISSION)
          ENDIF
       ELSE ! (mlev == - 1)
          !
          EMIS_SET(i)%etype = 1       ! 2D (SURFACE EMISSION)
       END IF level_analysis

       SELECT CASE(EMIS_SET(i)%etype)
       CASE(1)
          ! 2D (SURFACE EMISSION) -------------------------------------
          !
          ! NUMBER OF LEVELS IN CHANNEL OBJECT

          ! DIAGNOSTIC OUTPUT
          IF (p_parallel_io) THEN
             WRITE(*,*) '  EMISSION TYPE   : SURFACE (2D)'
             WRITE(*,*) '  -> GP METHOD    : ',EMIS_SET(i)%ngpem
             WRITE(*,*) '  -> LG METHOD    : ',EMIS_SET(i)%nlgem
             WRITE(*,*) '  -> CL METHOD    : ',EMIS_SET(i)%nclem
          END IF

          ! UNIT OF CHANNEL OBJECT
          ! get unit from attribute
          SELECT CASE (ADJUSTL(TRIM(unit_ori)))
             CASE('molec. m-2 s-1','molecules/m^2/s','molec./m^2/s')
                unit = 'molec. m-2 s-1'
             CASE('')
                CALL warning_bi( &
                     'unspecified unit of channel / object '//&
                     TRIM(EMIS_SET(i)%ch_name)//' / '//&
                     TRIM(EMIS_SET(i)%ch_object)//&
                     ' - assuming [molec./m^2/s]' &
                     , substr)
                unit = 'molec. m-2 s-1 (assumed by '//modstr//')'
             CASE DEFAULT
                CALL warning_bi( &
                     'unrecognized unit of channel / object '//&
                     TRIM(EMIS_SET(i)%ch_name)//' / '//&
                     TRIM(EMIS_SET(i)%ch_object)//': ['//&
                     ADJUSTL(TRIM(unit_ori))//&
                     '] - assuming [molec./m^2/s]' &
                     , substr)
                unit = 'molec. m-2 s-1 (assumed by '//modstr//')'
          END SELECT
          ! ----------------------------------------------------------
       CASE(2, 3)
          ! 3D (VOLUME EMISSION) -------------------------------------
          ! Nx2D (MULTI LEVEL EMISSION)
          !    (Nx2D (multi level) emissions are internally converted to
          !     3D (volume) emissions)
          !
          ! RESET EMISSION METHOD APPROPRIATE FOR 3D
          IF (EMIS_SET(i)%ngpem > 1) THEN
             IF (p_parallel_io) THEN
                WRITE(*,*) 'WARNING: GP METHOD FOR 3D RESET !'
             END IF
             EMIS_SET(i)%ngpem = 1
          END IF
          IF (EMIS_SET(i)%nlgem > 1) THEN
             IF (p_parallel_io) THEN
                WRITE(*,*) 'WARNING: LG METHOD FOR 3D RESET !'
             END IF
             EMIS_SET(i)%nlgem = 1
          END IF
          IF (EMIS_SET(i)%nclem > 1) THEN
             IF (p_parallel_io) THEN
                WRITE(*,*) 'WARNING: CL METHOD FOR 3D RESET !'
             END IF
             EMIS_SET(i)%nclem = 1
          END IF

          ! DIAGNOSTIC OUTPUT
          IF (p_parallel_io) THEN
             IF (EMIS_SET(i)%etype == 2) THEN
                WRITE(*,*) '  EMISSION TYPE   : VOLUME (3D)'
             ELSE
                WRITE(*,*) '  EMISSION TYPE   : MULTI LEVEL (Nx2D)'
                WRITE(*,*) '  --> ',mlev,' EMISSION LEVELS'
             END IF
             WRITE(*,*) '  -> GP METHOD    : ',EMIS_SET(i)%ngpem
             WRITE(*,*) '  -> LG METHOD    : ',EMIS_SET(i)%nlgem
             WRITE(*,*) '  -> CL METHOD    : ',EMIS_SET(i)%nclem
          END IF

          ! UNIT OF CHANNEL OBJECT get from attribute
          IF (EMIS_SET(i)%etype == 2) THEN
             ! VOLUME
             SELECT CASE (ADJUSTL(TRIM(unit_ori)))
             CASE('molec. m-3 s-1','molecules/m^3/s','molec./m^3/s')
                unit = 'molec. m-3 s-1'
             CASE('')
                CALL warning_bi( &
                     'unspecified unit of channel / object '//&
                     TRIM(EMIS_SET(i)%ch_name)//' / '//&
                     TRIM(EMIS_SET(i)%ch_object)//&
                     ' - assuming [molec./m^3/s]' &
                     , substr)
                unit = 'molec. m-3 s-1 (assumed by '//modstr//')'
             CASE DEFAULT
                CALL warning_bi( &
                     'unrecognized unit of channel / object '//&
                     TRIM(EMIS_SET(i)%ch_name)//' / '//&
                     TRIM(EMIS_SET(i)%ch_object)//': ['//&
                     ADJUSTL(TRIM(unit_ori))//&
                     '] - assuming [molec./m^3/s]' &
                     , substr)
                unit = 'molec. m-3 s-1 (assumed by '//modstr//')'
             END SELECT
          ELSE
             ! Nx2D
             SELECT CASE (ADJUSTL(TRIM(unit_ori)))
             CASE('molec. m-2 s-1','molecules/m^2/s','molec./m^3/s')
                ! see division by zdz [m] below !
                unit = 'molec. m-3 s-1'
             CASE('')
                CALL warning_bi( &
                     'unspecified unit of channel / object '//&
                     TRIM(EMIS_SET(i)%ch_name)//' / '//&
                     TRIM(EMIS_SET(i)%ch_object)//&
                     ' - assuming [molec./m^2/s]' &
                     , substr)
                ! see division by zdz [m] below !
                unit = 'molec. m-3 s-1 (assumed by '//modstr//')'
             CASE DEFAULT
                CALL warning_bi( &
                     'unrecognized unit of channel / object '//&
                     TRIM(EMIS_SET(i)%ch_name)//' / '//&
                     TRIM(EMIS_SET(i)%ch_object)//': ['//&
                     ADJUSTL(TRIM(unit_ori))//&
                     '] - assuming [molec./m^2/s]' &
                     , substr)
                ! see division by zdz [m] below !
                unit = 'molec. m-3 s-1 (assumed by '//modstr//')'
             END SELECT
          END IF
          ! ----------------------------------------------------------
       CASE DEFAULT
          ! ----------------------------------------------------------
          !
          ! ----------------------------------------------------------
       END SELECT

       !
       ! CHECK FOR GP TRACERS
       IF (L_GP .AND. EMIS_SET(i)%ngpem > 0) THEN
          do j=1, EMIS_SET(i)%ntr
             
             CALL get_tracer(status, GPTRSTR,EMIS_SET(i)%tr_basename(j) &
                  , EMIS_SET(i)%tr_subname(j), idx = EMIS_SET(i)%idt_gp(j))
             IF (status /= 0) THEN
                IF (TRIM(EMIS_SET(i)%tr_subname(j)) == '') THEN
                   CALL error_bi('TRACER '//TRIM(EMIS_SET(i)%tr_basename(j))//&
                        ' NOT FOUND',substr)
                ELSE
                   CALL error_bi('TRACER '//TRIM(EMIS_SET(i)%tr_basename(j))//&
                        '_'//TRIM(EMIS_SET(i)%tr_subname(j))//' NOT FOUND' &
                        ,substr)
                ENDIF
             ENDIF
          end do
       ENDIF
!!#D attila +
       ! CHECK FOR LG TRACERS
       IF (L_LG .AND. EMIS_SET(i)%nlgem > 0) THEN
          do j=1, EMIS_SET(i)%ntr
             CALL get_tracer(status, LGTRSTR, EMIS_SET(i)%tr_basename(j) &
                  , EMIS_SET(i)%tr_subname(j), idx = EMIS_SET(i)%idt_lg(j))
             IF (status /= 0) THEN
                IF (TRIM(EMIS_SET(i)%tr_subname(j)) == '') THEN
                   CALL error_bi('TRACER '//TRIM(EMIS_SET(i)%tr_basename(j))//&
                        ' NOT FOUND',substr)
                ELSE
                   CALL error_bi('TRACER '//TRIM(EMIS_SET(i)%tr_basename(j))//&
                        '_'//TRIM(EMIS_SET(i)%tr_subname(j))//' NOT FOUND' &
                        ,substr)
                ENDIF
             ENDIF
          end do
       END IF
!!#D attila -

!!#D clams +
       ! CHECK FOR CL TRACERS
       IF (L_CL .AND. EMIS_SET(i)%nclem > 0) THEN
          do j=1, EMIS_SET(i)%ntr
             CALL get_tracer(status, CLTRSTR, EMIS_SET(i)%tr_basename(j) &
                  , EMIS_SET(i)%tr_subname(j), idx = EMIS_SET(i)%idt_cl(j))
             IF (status /= 0) THEN
                IF (TRIM(EMIS_SET(i)%tr_subname(j)) == '') THEN
                   CALL error_bi('TRACER '//TRIM(EMIS_SET(i)%tr_basename(j))//&
                        ' NOT FOUND',substr)
                ELSE
                   CALL error_bi('TRACER '//TRIM(EMIS_SET(i)%tr_basename(j))//&
                        '_'//TRIM(EMIS_SET(i)%tr_subname(j))//' NOT FOUND' &
                        ,substr)
                ENDIF
             ENDIF
          end do
       END IF
!!#D clams -

       IF_Nx2Demis: IF (EMIS_SET(i)%etype/=3) THEN
#ifndef VERTICO
          EMIS_SET(i)%field => EMIS_SET(i)%raw
#else
       IF (EMIS_SET(i)%ngpem==0) THEN
          CALL new_channel_object(status, modstr &
               , TRIM(EMIS_SET(i)%ch_object) &
               , p3=EMIS_SET(i)%field, reprid=GP_3D_MID )
          CALL channel_halt(substr, status)
          CALL new_attribute(status, modstr &
               , TRIM(EMIS_SET(i)%ch_object), 'units', c=TRIM(unit) )
          CALL channel_halt(substr, status)
       ELSE
          EMIS_SET(i)%field => EMIS_SET(i)%raw
       ENDIF
#endif
       ELSE
          !
          CALL get_channel_object(status, modstr &
               , TRIM(EMIS_SET(i)%ch_object)     &
               , p3=EMIS_SET(i)%field            )
          IF (status/=0) THEN
             ! CREATE CHANNEL OBJECTS (Nx2D)
             CALL new_channel_object(status, modstr &
                  , TRIM(EMIS_SET(i)%ch_object) &
                  , p3=EMIS_SET(i)%field, reprid=reprid )
             CALL channel_halt(substr, status)
             CALL new_attribute(status, modstr &
                  , TRIM(EMIS_SET(i)%ch_object), 'units', c=TRIM(unit) )
             CALL channel_halt(substr, status)
             ! op_pj_20141016+
             CALL new_attribute(status, modstr &
                  , TRIM(EMIS_SET(i)%ch_object), modstr//'_em_method_gp' &
                  , i=EMIS_SET(i)%ngpem)
             CALL channel_halt(substr, status)
             CALL new_attribute(status, modstr &
                  , TRIM(EMIS_SET(i)%ch_object), modstr//'_em_type' &
                  , i=EMIS_SET(i)%etype)
             CALL channel_halt(substr, status)
             ! op_pj_20141016-
             !
             ! CREATE CHANNEL OBJECTS for emission height index
             CALL new_channel_object(status, modstr &
                  , TRIM(EMIS_SET(i)%ch_object)//'_vind' &
                  , p3=EMIS_SET(i)%index, reprid=reprid )
             CALL channel_halt(substr, status)
             CALL new_attribute(status, modstr &
                  , TRIM(EMIS_SET(i)%ch_object)//'_vind' &
                  , 'name', c='vertical index of Nx2d Emission' )
             CALL channel_halt(substr, status)
          ELSE
             EMIS_SET(i)%obj_exists = .TRUE.
             CALL get_channel_object(status, modstr      &
                  , TRIM(EMIS_SET(i)%ch_object)//'_vind' &
                  , p3=EMIS_SET(i)%index                 )
          END IF
            ! op_sb_20191111+
       END IF IF_Nx2Demis

!!#D attila +
#if defined(ECHAM5)
          ! SPECIAL CHANNEL OBJECTS FOR LAGRANGIAN TRACERS
          lagrange: IF (L_LG) THEN
             ! -> 'REST' FLUX
             ! ... ONLY IF LG-TRACERS ARE PRESENT AND METHOD = 1
             IF ( ANY(EMIS_SET(i)%idt_lg(:) > 0)  &
                  .AND. (EMIS_SET(i)%nlgem == 1) ) THEN
                !
                ! reprid:
                ! Nx2D -> converted to 3D -> GP_3D_MID  (from above)
                ! 2D   -> GP_3D_1LEV (from above)
                ! 3D   -> GP_3D_MID  (from above)
                !
                CALL new_channel_object(status, modstr &
                     , 'R_'//TRIM(EMIS_SET(i)%ch_object)&
                     , p3=EMIS_SET(i)%rest              &
                     , reprid=reprid, lrestreq = .TRUE. )
                CALL channel_halt(substr, status)
                CALL new_attribute(status, modstr &
                     , 'R_'//TRIM(EMIS_SET(i)%ch_object)&
                     , 'units', c='kg mol/mol/s' )
                CALL channel_halt(substr, status)
             END IF
             ! -> TRACER TENDENCY
             ! ... ONLY IF REQUESTED AND LG-TRACERS ARE PRESENT
             IF (l_lg_tend .AND. (ANY(EMIS_SET(i)%idt_lg(:) > 0))) THEN
              DO j=1,EMIS_SET(i)%ntr         
                 IF (TRIM(EMIS_SET(i)%tr_subname(j)) /= '') THEN
                    outname = TRIM(EMIS_SET(i)%tr_basename(j))//'_'//TRIM(EMIS_SET(i)%tr_subname(j))
                 ELSE
                    outname = TRIM(EMIS_SET(i)%tr_basename(j))
                 END IF

                 IF (p_parallel_io) THEN
                    WRITE(*,*) '   CHANNEL OBJECT  : '   &
                        ! ,'(',EMIS_SET(i)%nf,') '        &
!!$                , 'lgte_'//TRIM(EMIS_SET(i)%ch_object)
                   , 'lgte_'//TRIM(outname)          ! op_sb_20191111
                 END IF
                CALL new_channel_object(status, modstr      &
!!$                  , 'lgte_'//TRIM(EMIS_SET(i)%ch_object) &
                     , 'lgte_'//TRIM(outname) &      ! op_sb_20191111
                     , p1=EMIS_SET(i)%xtte_lg               &
                     , reprid=LG_ATTILA )
                CALL channel_halt(substr, status)
                CALL new_attribute(status, modstr &
!!$                  , 'lgte_'//TRIM(EMIS_SET(i)%ch_object) &
                     , 'lgte_'//TRIM(outname) &      ! op_sb_20191111
                     , 'units', c='mol/mol/s' )
                CALL channel_halt(substr, status)
             ENDDO  ! j-loop
             END IF
          END IF lagrange
#endif
!!#D attila -

!!#D clams +
#if defined(ECHAM5)
          ! SPECIAL CHANNEL OBJECTS FOR LAGRANGIAN TRACERS
          clams: IF (L_CL) THEN
             ! -> 'REST' FLUX
             ! ... ONLY IF CL-TRACERS ARE PRESENT AND METHOD = 1
             IF ( ANY(EMIS_SET(i)%idt_cl(:) > 0)  &
                  .AND. (EMIS_SET(i)%nclem == 1) ) THEN
                !
                ! reprid:
                ! Nx2D -> converted to 3D -> GP_3D_MID  (from above)
                ! 2D   -> GP_3D_1LEV (from above)
                ! 3D   -> GP_3D_MID  (from above)
                !
                print*,'vor R'
                CALL new_channel_object(status, modstr &
                     , 'R_'//TRIM(EMIS_SET(i)%ch_object)&
                     , p3=EMIS_SET(i)%rest              &
                     , reprid=reprid, lrestreq = .TRUE. )
                CALL channel_halt(substr, status)
                CALL new_attribute(status, modstr &
                     , 'R_'//TRIM(EMIS_SET(i)%ch_object)&
                     , 'units', c='kg mol/mol/s' )
                CALL channel_halt(substr, status)
               END IF
             ! -> TRACER TENDENCY
             ! ... ONLY IF REQUESTED AND CL-TRACERS ARE PRESENT
             IF (l_cl_tend .AND. (ANY(EMIS_SET(i)%idt_cl(:) > 0))) THEN
             DO j=1,EMIS_SET(i)%ntr
                IF (TRIM(EMIS_SET(i)%tr_subname(j)) /= '') THEN
                      outname = TRIM(EMIS_SET(i)%tr_basename(j))//'_'//TRIM(EMIS_SET(i)%tr_subname(j))
                ELSE
                      outname = TRIM(EMIS_SET(i)%tr_basename(j))
                END IF
                print*,'outname=',outname
                IF (p_parallel_io) THEN
                   WRITE(*,*) '   CHANNEL OBJECT  : '   &
                        ! ,'(',EMIS_SET(i)%nf,') '        &
                   , 'clte_'//TRIM(outname)
                END IF
                CALL new_channel_object(status, modstr      &
                     , 'clte_'//TRIM(outname) &
                     , p1=EMIS_SET(i)%xtte_cl               &
                     , reprid=REPR_LG_CLAMS )
                CALL channel_halt(substr, status)
                CALL new_attribute(status, modstr &
                     , 'clte_'//TRIM(outname) &
                     , 'units', c='mol/mol/s' )
                CALL channel_halt(substr, status)
              ENDDO  ! j-loop
             END IF
          END IF clams
#endif
!!#D clams -

       !   
       IF (p_parallel_io) THEN
          WRITE(*,*) '========================================================'
       END IF

    END DO offemis_set_loop

    CALL end_message_bi(modstr, 'COUPLING SETUP',substr)

  END SUBROUTINE offemis_init_coupling
! ------------------------------------------------------------------------


! ------------------------------------------------------------------------
  SUBROUTINE offemis_vdiff

    IMPLICIT NONE

    IF (L_GP) CALL offemis_vdiff_gp
    
  END SUBROUTINE offemis_vdiff
! ------------------------------------------------------------------------

! ------------------------------------------------------------------------
  SUBROUTINE offemis_vdiff_gp

    ! NOTES:
    ! 2D- (i.e., surcface-) fluxes:
    ! * the unit MUST BE molecules m-2 s-1
    ! * the resulting tracer tendency for method 1 (lowest grid layer) is
    !   mol/mol/s
    ! * for method 2 (lower boundary condition of vertical flux),
    !   pxtems in in vdiff multiplied with ZTMST*G*ZQDP (s m s-2 Pa-1),
    !   therefore pxtems MUST BE in
    !   (mol(tracer)/mol(air)) * kg (air) m-2 s-1
    ! 
    ! 3D- (i.e., volume-) emissions
    ! * the unit MUST BE molecules m-3 s-1
    !
    ! Nx2D- (i.e, multi level) emissions
    ! * the unit MUST BE molecules m-2 s-1
    ! * These emissions are converted to 3D emissions
    !   by dividing by the box height.
    !   The index of the vertical layer is calculated from the user
    !   specified emission height (Z) and the current geopotential.
    !   

    ! ECHAM5/MESSy
    USE messy_main_data_bi,       ONLY: rho_air_dry_3d
#if defined(ECHAM5) || defined(CESM1)
    USE messy_main_data_bi,       ONLY: pxtems
#endif
    USE messy_main_grid_def_mem_bi, ONLY:jrow, kproma, nlev
    USE messy_main_grid_def_bi,     ONLY: deltaz, altitudei_gnd
#ifndef MESSYTENDENCY
    USE messy_main_tracer_mem_bi, ONLY: pxtte => qxtte 
#endif
    ! MESSy
    USE messy_main_channel_mem,   ONLY: dom_curid
    USE messy_main_constants_mem, ONLY: g, N_A, M_air

    IMPLICIT NONE

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER           :: substr = 'offemis_vdiff_gp'
    ! air molecules / kg air
    REAL(DP), PARAMETER                   :: uconv = N_A * 1.0e3_DP/M_air
    REAL(DP), DIMENSION(:,:), POINTER     :: zairdens => NULL() ! air density
    REAL(DP), DIMENSION(:,:), ALLOCATABLE :: zdz         ! layer thickness
    REAL(DP), DIMENSION(:,:), ALLOCATABLE :: zscale      !
    INTEGER                               :: i, jt
    INTEGER                               :: j
    INTEGER                               :: ji, jk, jp   ! for Nx2D
    INTEGER                               :: ml1, mlev    ! levels
#ifdef MESSYTENDENCY 
   REAL(DP), DIMENSION(:,:), ALLOCATABLE :: zxtte
#endif

    ! CALCULATE AIR DENSITY AND LAYER THICKNESS
    ALLOCATE(zdz(kproma, nlev))
    ALLOCATE(zscale(kproma, nlev))
#ifdef MESSYTENDENCY
    ALLOCATE(zxtte(kproma, nlev))
#endif
    zairdens => rho_air_dry_3d(_RI_XYZ__(1:kproma,jrow,1:nlev))
    zdz(:,:) = 0.0
    zdz(:,2:nlev) = deltaz(_RI_XYZ__(1:kproma,jrow,2:nlev))
    zscale(:,:) = 0.0

    offemis_set_loop: DO i=1, NEMIS_SET

       IF (EMIS_SET(i)%domain_idx /= dom_curid) CYCLE

       ! NO CHANNEL OBJECT PRESENT -> SKIP
       IF (.NOT. ASSOCIATED(EMIS_SET(i)%raw)) CYCLE

       ! EMISSION TYPE
       ! NOTE: EMISSION METHOD CHECKED (AND CORRECTED) IN INIT_MEMORY
       SELECT CASE(EMIS_SET(i)%etype)
       CASE(1) ! 2D surface
          ml1  = nlev
          mlev = 1
          zscale(:,2:nlev) = zairdens(:,2:nlev)*zdz(:,2:nlev)
       CASE(2) ! 3D volume
          ml1  = 1
          mlev = nlev
          zscale(:,:) = zairdens(:,:)
       CASE(3) ! Nx2D multi level
          !
          ml1  = 1
          mlev = SIZE(EMIS_SET(i)%z)
          zscale(:,:) = zairdens(:,:)                   
          ! 
          IF (.NOT. EMIS_SET(i)%obj_exists) THEN
             EMIS_SET(i)%field(_RI_XYZ__(:,jrow,:)) = 0.0
             EMIS_SET(i)%index(_RI_XYZ__(:,jrow,:)) = 0.0
             emission_levels1: DO ji=1, mlev
                vector_loop1: DO jp=1, kproma
                   ! calculate level index JK from Z(ji)
                   DO jk=nlev,2,-1
                      IF (EMIS_SET(i)%z(ji) <= altitudei_gnd(_RI_XYZ__(jp,jrow,jk))) THEN
                         EMIS_SET(i)%index(_RI_XYZ__(jp,jrow,ji)) = jk
                         EMIS_SET(i)%field(_RI_XYZ__(jp,jrow,ji)) =  &
                              EMIS_SET(i)%raw(_RI_XYZ__(jp,jrow,ji)) &
                              /zdz(jp,jk)  
                         EXIT
                      END IF
                   END DO
                   IF (jk<=2) CALL error_bi('ERROR NO HEIGHT FOUND', substr)
                END DO vector_loop1
             END DO emission_levels1
          END IF
       CASE DEFAULT
       END SELECT

#ifdef VERTICO
       IF (EMIS_SET(i)%etype/=3) THEN
          !WARNING : TO BE IMPROVED -> switching the axis order 
          ! for arosols emissions 3D (ngpem==0)
          IF (EMIS_SET(i)%ngpem==0) THEN
             EMIS_SET(i)%field(:,:,1:nlev) = EMIS_SET(i)%raw(:,:,nlev:1:-1)
          ENDIF
       ENDIF
#endif

       ! ONLY, IF GP TRACERS ARE PRESENT (see definition of channel objects)

       tr_loop: DO j=1, EMIS_SET(i)%ntr
          IF (EMIS_SET(i)%idt_gp(j) == 0) CYCLE

          jt = EMIS_SET(i)%idt_gp(j)

          ! EMISSION METHOD
          select_gp:SELECT CASE(EMIS_SET(i)%ngpem)
          CASE(0)
             ! DO NOTHING
          CASE(1)
             ! SURFACE (2D): TENDENCY OF LOWEST GRID LAYER
             ! VOLUME  (3D): TENDENCY IN GRID BOX
             select_etype: SELECT CASE (EMIS_SET(i)%etype)
             CASE(1,2) ! 2D surface / 3D volume
#ifndef MESSYTENDENCY
                pxtte(_RI_X_ZN_(1:kproma,ml1:nlev,jt)) = &
                     pxtte(_RI_X_ZN_(1:kproma,ml1:nlev,jt))  &
                     + EMIS_SET(i)%field(_RI_XYZ__(1:kproma,jrow,1:mlev))  &
                     * EMIS_SET(i)%tr_factor(j)                            &
                     / (zscale(1:kproma,ml1:nlev) * uconv)
#else
                zxtte(:,:) = 0._dp
                zxtte(1:kproma,ml1:nlev) = &
                     EMIS_SET(i)%field(_RI_XYZ__(1:kproma,jrow,1:mlev)) &
                     * EMIS_SET(i)%tr_factor(j)            &
                     / (zscale(1:kproma,ml1:nlev) * uconv)
                CALL mtend_add_l(my_handle, jt, px=zxtte)
#endif
             CASE(3) ! Nx2D multi level
                emission_levels2: DO ji=1, mlev
                   vector_loop2: DO jp=1, kproma
                      jk= NINT(EMIS_SET(i)%index(_RI_XYZ__(jp,jrow,ji)))
#ifndef MESSYTENDENCY
                      pxtte(_RI_X_ZN_(jp,jk,jt)) = pxtte(_RI_X_ZN_(jp,jk,jt)) &
                           + EMIS_SET(i)%field(_RI_XYZ__(jp,jrow,ji))         &
                           * EMIS_SET(i)%tr_factor(j)                         &
                           / (zscale(jp,jk) * uconv)
#else
                      zxtte(:,:) = 0._dp
                      zxtte(jp,jk) =                                &
                           EMIS_SET(i)%field(_RI_XYZ__(jp,jrow,ji)) &
                           * EMIS_SET(i)%tr_factor(j)               &
                           / (zscale(jp,jk) * uconv)
                      CALL mtend_add_l(my_handle, jt, px=zxtte)
#endif
                   END DO vector_loop2
                END DO emission_levels2
             CASE DEFAULT
             END SELECT select_etype

#if defined(ECHAM5) || defined(CESM1)
          CASE(2)
             ! SURFACE (2D): LOWER BOUNDARY CONDITION OF VERTICAL FLUX
             ! VOLUME  (3D): NOT APPLICABLE
             pxtems(_RI_XYZN_(1:kproma,jrow,1,jt)) = pxtems(_RI_XYZN_(1:kproma,jrow,1,jt)) &
                  + EMIS_SET(i)%field(_RI_XYZ__(1:kproma,jrow,1))    &
                  * EMIS_SET(i)%tr_factor(j)  &
                  / uconv
#endif
          CASE DEFAULT
             ! ERROR
             CALL error_bi('UNKNOWN EMISSION METHOD',substr)
          END SELECT select_gp
       END DO tr_loop

    END DO offemis_set_loop

    ! CLEAN UP
    NULLIFY(zairdens)
    DEALLOCATE(zdz)
    DEALLOCATE(zscale)
#ifdef MESSYTENDENCY
    DEALLOCATE(zxtte)
#endif

  END SUBROUTINE offemis_vdiff_gp
! ------------------------------------------------------------------------

! ------------------------------------------------------------------------
  SUBROUTINE offemis_global_end

    IMPLICIT NONE

!!#D attila +
#if defined(ECHAM5)
    IF (L_LG) CALL offemis_global_end_lg
#endif
!!#D attila -

!!#D clams +
#if defined(ECHAM5)
    IF (L_CL) CALL offemis_global_end_cl
#endif
!!#D clams -

  END SUBROUTINE offemis_global_end
! ------------------------------------------------------------------------

!!#D attila +
#if defined(ECHAM5)
! ------------------------------------------------------------------------
  SUBROUTINE offemis_global_end_lg

    ! ECHAM5/MESSy
    USE messy_main_data_bi,         ONLY: zairdens => rho_air_dry_3d
    USE messy_main_grid_def_mem_bi, ONLY:nproma, npromz, ngpblks, nlev
    USE messy_main_grid_def_bi,     ONLY: deltaz, altitudei_gnd
    USE messy_main_tracer_mem_bi,   ONLY: NCELL, pxtte => qxtte_a
    USE messy_attila_tools_e5,     ONLY: gpsfemis2lgemis_e5, gp2lg_e5
    ! MESSy
    USE messy_main_constants_mem,  ONLY: g, M_air, N_A

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER             :: substr = 'offemis_global_end_lg'
    REAL(DP), PARAMETER                     :: uconv = N_A * 1.0e3_DP/M_air
    INTEGER                                 :: zkproma
    REAL(DP), DIMENSION(:,:,:), ALLOCATABLE :: zdz       ! layer thickness
    REAL(DP), DIMENSION(:,:,:), ALLOCATABLE :: zscale    !
    REAL(DP), DIMENSION(:,:,:), POINTER     :: zxtte_gp    => NULL()
    REAL(DP), DIMENSION(:,:),   POINTER     :: zxtte_gpsf  => NULL()
    REAL(DP), DIMENSION(:,:),   POINTER     :: rest_2d     => NULL()
    REAL(DP), DIMENSION(:),     POINTER     :: zxtte_lg    => NULL()
    INTEGER                                 :: zjrow
    INTEGER                                 :: i, jt
    INTEGER                                 :: j 
    INTEGER                                 :: ji, ip,jk, jp ! for Nx2D
    INTEGER                                 :: ml1, mlev     ! levels


    ! MEMORY FOR GP-TENDENCY
    ALLOCATE(zxtte_gp(nproma, nlev, ngpblks))
    ! MEMORY FOR LG-TENDENCY
    IF (.NOT. l_lg_tend) THEN
       ! LOCAL SPACE FOR ALL EMISSION SETS / TRACERS
       ALLOCATE(zxtte_lg(NCELL))
       !ELSE
       ! USE POINTERS TO CHANNEL OBJECTS
    END IF

    ! CALCULATE AIR DENSITY AND LAYER THICKNESS
    ALLOCATE(zdz(nproma, nlev, ngpblks))
    zdz(:,:,:) = 1.0_DP
    ALLOCATE(zscale(nproma, nlev, ngpblks))
    zscale(:,:,:) = 1.0_DP
    !
    DO zjrow=1, ngpblks
       IF ( zjrow == ngpblks ) THEN
          zkproma = npromz
       ELSE
          zkproma = nproma
       END IF
       !
       zdz(1:zkproma,2:nlev,zjrow) = deltaz(1:zkproma,2:nlev,zjrow)
    END DO

    offemis_set_loop: DO i=1, NEMIS_SET

       ! NO CHANNEL OBJECT PRESENT -> SKIP
       IF (.NOT. ASSOCIATED(EMIS_SET(i)%field)) CYCLE

       ! ONLY, IF LG TRACERS ARE PRESENT (see definition of channel objects)
       IF (ALL(EMIS_SET(i)%idt_lg(:) == 0)) CYCLE

       ! EMISSION TYPE
       ! NOTE: EMISSION METHOD CHECKED (AND CORRECTED IN INIT_MEMORY)
       SELECT CASE(EMIS_SET(i)%etype)
       CASE(1) ! 2D surface
          ml1  = nlev
          mlev = 1
          zscale(:,2:nlev,:) = zairdens(:,2:nlev,:)*zdz(:,2:nlev,:)
       CASE(2) ! 3D (volume)
          ml1  = 1
          mlev = nlev
          zscale(:,:,:) = zairdens(:,:,:)
       CASE(3) ! Nx2D multi level
          ! convert Nx2D to 3D (raw -> field)
          ! by dividing by the box height.
          ! The index of the vertical layer is calculated
          ! from the user specified emission height (Z) and the
          ! current geopotential.
          !
          ! INITIALIZE
          EMIS_SET(i)%field(:,:,:) = 0.0
          !
          ip = SIZE(EMIS_SET(i)%z)  ! number of emission levels
          local_loop: DO zjrow=1, ngpblks
             IF ( zjrow == ngpblks ) THEN
                zkproma = npromz
             ELSE
                zkproma = nproma
             END IF
             emission_levels: DO ji=1, ip
                vector_loop: DO jp=1, zkproma
                   ! calculate level index JK from Z(ji)
                   DO jk=nlev,2,-1
                      IF (EMIS_SET(i)%z(ji) &
                           <= altitudei_gnd(jp,jk,zjrow)) EXIT
                   END DO
                   ! 
                   EMIS_SET(i)%field(jp,jk,zjrow) = &
                        EMIS_SET(i)%field(jp,jk,zjrow) + &
                        EMIS_SET(i)%raw(jp,ji,zjrow)     &
                        /zdz(jp,jk,zjrow)
                END DO vector_loop
             END DO emission_levels
          END DO local_loop
          !
          ml1  = 1
          mlev = nlev
          zscale(:,:,:) = zairdens(:,:,:)                   
       CASE DEFAULT
       END SELECT
       
       zxtte_gp(:,ml1:nlev,:) = &
            EMIS_SET(i)%field(:,1:mlev,:) &
            / (zscale(:,ml1:nlev,:) * uconv)
       zxtte_gpsf => zxtte_gp(:,nlev,:)
       !
       IF (l_lg_tend) THEN
          ! USE POINTER TO CHANNEL OBJECT
          zxtte_lg => EMIS_SET(i)%xtte_lg
          !ELSE
          ! USE LOCAL MEMORY FOR ALL EMSSION SETS / TRACERS
       END IF
          
       ! EMISSION METHOD
       SELECT CASE(EMIS_SET(i)%nlgem)
       CASE(0)
          ! DO NOTHING
       CASE(1)
          SELECT CASE(EMIS_SET(i)%etype)
          CASE(1) ! 2D surfce
             rest_2d => EMIS_SET(i)%rest(:,nlev,:)
             CALL gpsfemis2lgemis_e5(zxtte_gpsf, zxtte_lg   &
                  , EMIS_SET(i)%nlgem                       &  ! = 1
                  , gprl=rest_2d                            &
                  , lmcons = .true. )
          CASE(2, 3) ! 3D (volume) and Nx2D (multi level)
             CALL gp2lg_e5(zxtte_gp, zxtte_lg     &
                  , gprl=EMIS_SET(i)%rest &
                  , lmcons = .true. )
          CASE DEFAULT
          END SELECT
          !
       CASE(2) ! >1: SURFACE (2D) ONLY 
          CALL gpsfemis2lgemis_e5(zxtte_gpsf, zxtte_lg &
               , EMIS_SET(i)%nlgem                     &  ! = 2
               , lmcons = .true. )                
       CASE(3)
          CALL gpsfemis2lgemis_e5(zxtte_gpsf, zxtte_lg &
               , EMIS_SET(i)%nlgem                     &  ! = 3
               , lmcons = .true. )                
       CASE(4)
          CALL gpsfemis2lgemis_e5(zxtte_gpsf, zxtte_lg &
               , EMIS_SET(i)%nlgem                     &  ! = 4
               , lmcons = .true. )                
       CASE DEFAULT
          ! ERROR
          CALL error_bi('UNKNOWN EMISSION METHOD',substr)
       END SELECT
          
       DO j=1, EMIS_SET(i)%ntr
          ! SET TRACER INDEX
          jt = EMIS_SET(i)%idt_lg(j)

          ! ADD TRACER TENDENCY
          IF (EMIS_SET(i)%nlgem > 0) THEN
             pxtte(:,jt) = pxtte(:,jt) + zxtte_lg(:)*EMIS_SET(i)%tr_factor(j)
          END IF
       END DO
    END DO offemis_set_loop
    
    ! CLEAN UP
    DEALLOCATE(zdz)
    DEALLOCATE(zscale)
    !
    DEALLOCATE(zxtte_gp)
    IF (.NOT. l_lg_tend) THEN
       ! FREE LOCAL SPACE
       DEALLOCATE(zxtte_lg)
       !ELSE
       ! POINTERS TO CHANNEL OBJECTS HAVE BEEN USED
    END IF

  END SUBROUTINE offemis_global_end_lg
! ------------------------------------------------------------------------
#endif
!!#D attila -

!!#D clams +
#if defined(ECHAM5)
! ------------------------------------------------------------------------
  SUBROUTINE offemis_global_end_cl

    ! ECHAM5/MESSy
    USE messy_main_data_bi,         ONLY: zairdens => rho_air_dry_3d
    USE messy_main_grid_def_mem_bi, ONLY: nproma, npromz, ngpblks, nlev
    USE messy_main_grid_def_bi,     ONLY: deltaz, altitudei_gnd
    USE messy_main_tracer_mem_bi,   ONLY: dnparts_max, pxtte => qxtte_c

    USE messy_clams_tools_e5,      ONLY: gpsfemis2clemis_e5, gp2cl_e5
    ! MESSy
    USE messy_main_constants_mem,  ONLY: g, M_air, N_A

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER             :: substr = 'offemis_global_end_cl'
    REAL(DP), PARAMETER                     :: uconv = N_A * 1.0e3_DP/M_air
    INTEGER                                 :: zkproma

    REAL(DP), DIMENSION(:,:,:), ALLOCATABLE :: zdz       ! layer thickness
    REAL(DP), DIMENSION(:,:,:), ALLOCATABLE :: zscale    !
    REAL(DP), DIMENSION(:,:,:), POINTER     :: zxtte_gp    => NULL()
    REAL(DP), DIMENSION(:,:),   POINTER     :: zxtte_gpsf  => NULL()
    REAL(DP), DIMENSION(:,:),   POINTER     :: rest_2d     => NULL()
    REAL(DP), DIMENSION(:),     POINTER     :: zxtte_cl    => NULL()
    INTEGER                                 :: zjrow
    INTEGER                                 :: i, jt
    INTEGER                                 :: j 
    INTEGER                                 :: ji, ip,jk, jp ! for Nx2D
    INTEGER                                 :: ml1, mlev ! levels

    ! MEMORY FOR GP-TENDENCY
    ALLOCATE(zxtte_gp(nproma, nlev, ngpblks))
    ! MEMORY FOR CL-TENDENCY
   
    IF (.NOT. l_cl_tend) THEN
       ! LOCAL SPACE FOR ALL EMISSION SETS / TRACERS
       ALLOCATE(zxtte_cl(dnparts_max))
    END IF

    ! CALCULATE AIR DENSITY AND LAYER THICKNESS
    ALLOCATE(zdz(nproma, nlev, ngpblks))
    zdz(:,:,:) = 1.0_DP
    ALLOCATE(zscale(nproma, nlev, ngpblks))
    zscale(:,:,:) = 1.0_DP
    !
    DO zjrow=1, ngpblks
       IF ( zjrow == ngpblks ) THEN
          zkproma = npromz
       ELSE
          zkproma = nproma
       END IF
       !
       zdz(1:zkproma,2:nlev,zjrow) = deltaz(1:zkproma,2:nlev,zjrow)
    END DO

    offemis_set_loop: DO i=1, NEMIS_SET

       ! NO CHANNEL OBJECT PRESENT -> SKIP
       IF (.NOT. ASSOCIATED(EMIS_SET(i)%field)) CYCLE

       ! ONLY, IF CL TRACERS ARE PRESENT (see definition of channel objects)
       IF (ALL(EMIS_SET(i)%idt_cl(:) == 0)) CYCLE
       ! EMISSION TYPE
       ! NOTE: EMISSION METHOD CHECKED (AND CORRECTED IN INIT_MEMORY)
       SELECT CASE(EMIS_SET(i)%etype)
       CASE(1) ! 2D surface
          ml1  = nlev
          mlev = 1
          zscale(:,2:nlev,:) = zairdens(:,2:nlev,:)*zdz(:,2:nlev,:)
       CASE(2) ! 3D (volume)
          ml1  = 1
          mlev = nlev
          zscale(:,:,:) = zairdens(:,:,:)
       CASE(3) ! Nx2D multi level
          ! convert Nx2D to 3D (raw -> field)
          ! by dividing by the box height.
          ! The index of the vertical layer is calculated
          ! from the user specified emission height (Z) and the
          ! current geopotential.
          !
          ! INITIALIZE
          EMIS_SET(i)%field(:,:,:) = 0.0
          !
          ip = SIZE(EMIS_SET(i)%z)  ! number of emission levels
          local_loop: DO zjrow=1, ngpblks
             IF ( zjrow == ngpblks ) THEN
                zkproma = npromz
             ELSE
                zkproma = nproma
             END IF
             emission_levels: DO ji=1, ip
                vector_loop: DO jp=1, zkproma
                   ! calculate level index JK from Z(ji)
                   DO jk=nlev,2,-1
                      IF (EMIS_SET(i)%z(ji) &
                           <= altitudei_gnd(jp,jk,zjrow)) EXIT
                   END DO
                   ! 
                   EMIS_SET(i)%field(jp,jk,zjrow) = &
                        EMIS_SET(i)%field(jp,jk,zjrow) + &
                        EMIS_SET(i)%raw(jp,ji,zjrow)     &
                        /zdz(jp,jk,zjrow)
                END DO vector_loop
             END DO emission_levels
          END DO local_loop
          !
          ml1  = 1
          mlev = nlev
          zscale(:,:,:) = zairdens(:,:,:)                   
       CASE DEFAULT
       END SELECT
       
       zxtte_gp(:,ml1:nlev,:) = &
            EMIS_SET(i)%field(:,1:mlev,:) &
            / (zscale(:,ml1:nlev,:) * uconv)
       zxtte_gpsf => zxtte_gp(:,nlev,:)
       !
       IF (l_cl_tend) THEN
          ! USE POINTER TO CHANNEL OBJECT
          zxtte_cl => EMIS_SET(i)%xtte_cl
          !ELSE
          ! USE LOCAL MEMORY FOR ALL EMSSION SETS / TRACERS
       END IF
          
       ! EMISSION METHOD
       SELECT CASE(EMIS_SET(i)%nclem)
       CASE(0)
          ! DO NOTHING
       CASE(1)
          SELECT CASE(EMIS_SET(i)%etype)
          CASE(1) ! 2D surfce
             rest_2d => EMIS_SET(i)%rest(:,nlev,:)
             CALL gpsfemis2clemis_e5(zxtte_gpsf, zxtte_cl   &
                  , EMIS_SET(i)%nclem                       &  ! = 1
                  , gprl=rest_2d                            &
                  , lmcons = .true. )
          CASE(2, 3) ! 3D (volume) and Nx2D (multi level)
             CALL gp2cl_e5(zxtte_gp, zxtte_cl     &
                  , gprl=EMIS_SET(i)%rest &
                  , lmcons = .true. )
          CASE DEFAULT
          END SELECT
          !
       CASE(2) ! >1: SURFACE (2D) ONLY 
          CALL gpsfemis2clemis_e5(zxtte_gpsf, zxtte_cl &
               , EMIS_SET(i)%nclem                     &  ! = 2
               , lmcons = .true. )                
       CASE(3)
          CALL gpsfemis2clemis_e5(zxtte_gpsf, zxtte_cl &
               , EMIS_SET(i)%nclem                     &  ! = 3
               , lmcons = .true. )                
       CASE(4)
                   CALL gpsfemis2clemis_e5(zxtte_gpsf, zxtte_cl &
               , EMIS_SET(i)%nclem                     &  ! = 4
               , lmcons = .true. )
       CASE DEFAULT
          ! ERROR
          CALL error_bi('UNKNOWN EMISSION METHOD',substr)
       END SELECT

       DO j=1, EMIS_SET(i)%ntr
          ! SET TRACER INDEX
          jt = EMIS_SET(i)%idt_cl(j)

          ! ADD TRACER TENDENCY
          IF (EMIS_SET(i)%nclem > 0) THEN
             pxtte(:,jt) = pxtte(:,jt) + zxtte_cl(:)*EMIS_SET(i)%tr_factor(j)
          END IF
       END DO
    END DO offemis_set_loop
    
    ! CLEAN UP
    DEALLOCATE(zdz)
    DEALLOCATE(zscale)
    !
    DEALLOCATE(zxtte_gp)
    IF (.NOT. l_cl_tend) THEN
       ! FREE LOCAL SPACE
       DEALLOCATE(zxtte_cl)
       !ELSE
       ! POINTERS TO CHANNEL OBJECTS HAVE BEEN USED
    END IF

  END SUBROUTINE offemis_global_end_cl
! ------------------------------------------------------------------------
#endif
!!#D clams -

! ----------------------------------------------------------------------
  SUBROUTINE offemis_free_memory

    IMPLICIT NONE

    ! LOCAL
    INTEGER :: i
    
    DO i=1, NEMIS_SET

       ! Z
       IF (ASSOCIATED(EMIS_SET(i)%z)) THEN
          DEALLOCATE(EMIS_SET(i)%z)
          NULLIFY(EMIS_SET(i)%z)
       END IF

       ! REST (Lagrange)
       IF (ASSOCIATED(EMIS_SET(i)%rest)) THEN
          DEALLOCATE(EMIS_SET(i)%rest)
          NULLIFY(EMIS_SET(i)%rest)
       END IF

       ! TRACER TENDENCY (Lagrange - ATTILA)
       IF (ASSOCIATED(EMIS_SET(i)%xtte_lg)) THEN
          DEALLOCATE(EMIS_SET(i)%xtte_lg)
          NULLIFY(EMIS_SET(i)%xtte_lg)
       END IF

       ! TRACER TENDENCY (Lagrange - CLAMS)
       IF (ASSOCIATED(EMIS_SET(i)%xtte_cl)) THEN
          DEALLOCATE(EMIS_SET(i)%xtte_cl)
          NULLIFY(EMIS_SET(i)%xtte_cl)
       END IF

       IF (ASSOCIATED(EMIS_SET(i)%tr_basename)) THEN
          DEALLOCATE(EMIS_SET(i)%tr_basename)
          NULLIFY(EMIS_SET(i)%tr_basename)
       END IF
       IF (ASSOCIATED(EMIS_SET(i)%tr_subname)) THEN
          DEALLOCATE(EMIS_SET(i)%tr_subname)
          NULLIFY(EMIS_SET(i)%tr_subname)
       END IF
       IF (ASSOCIATED(EMIS_SET(i)%tr_factor)) THEN
          DEALLOCATE(EMIS_SET(i)%tr_factor)
          NULLIFY(EMIS_SET(i)%tr_factor)
       END IF
       IF (ASSOCIATED(EMIS_set(i)%idt_gp)) THEN
          DEALLOCATE(EMIS_set(i)%idt_gp)
          NULLIFY(EMIS_set(i)%idt_gp)
       END IF
       IF (ASSOCIATED(EMIS_set(i)%idt_lg)) THEN
          DEALLOCATE(EMIS_set(i)%idt_lg)
          NULLIFY(EMIS_set(i)%idt_lg)
       END IF

       IF (ASSOCIATED(EMIS_set(i)%idt_cl)) THEN
          DEALLOCATE(EMIS_set(i)%idt_cl)
          NULLIFY(EMIS_set(i)%idt_cl)
       END IF

    END DO
    
  END SUBROUTINE offemis_free_memory
! ----------------------------------------------------------------------
! ************************************************************************
! PRIVATE ECHAM-5 INTERFACE ROUTINES
! ************************************************************************

! ----------------------------------------------------------------------
  SUBROUTINE offemis_read_nml_cpl(status, iou)

    ! emis MODULE ROUTINE (ECHAM-5 INTERFACE, PRIVATE)
    !
    ! read namelist for 'coupling' to ECHAM5
    !
    ! Author: Patrick Joeckel, MPICH, Mar 2004

    ! MESSy
    USE messy_main_tools,        ONLY: read_nml_open, read_nml_check &
                                     , read_nml_close
  
    IMPLICIT NONE

    ! I/O
    INTEGER, INTENT(OUT) :: status     ! error status
    INTEGER, INTENT(IN)  :: iou        ! I/O unit

    ! (LOCAL) NAMELIST VARIABLES
    CHARACTER(LEN=*), PARAMETER :: substr = 'offemis_read_nml_cpl'

    NAMELIST /CPL/ L_GP, L_LG, L_CL, l_lg_tend, l_cl_tend, EMIS_IN

    ! LOCAL
    LOGICAL          :: lex      ! file exists ?
    INTEGER          :: fstat    ! file status

    status = 1

    ! INITIALIZE NAMELIST VARIABLES
    CALL read_nml_open(lex, substr, iou, 'CPL', modstr)
    IF (.not.lex) RETURN    ! <modstr>.nml does not exist

    READ(iou, NML=CPL, IOSTAT=fstat)
    CALL read_nml_check(fstat, substr, iou, 'CPL', modstr)
    IF (fstat /= 0) RETURN  ! error while reading namelist

    ! DIAGNOSE NAMELIST AND SET GLOBAL SWITCHES
    !
    ! CHECK NAMELIST
    IF (L_GP) THEN
       WRITE(*,*) 'EMISSIONS IN GRIDPOINT SPACE : ON'
    ELSE
       WRITE(*,*) 'EMISSIONS IN GRIDPOINT SPACE : OFF'
    END IF

!!#D attila +
    IF (L_LG) THEN
       WRITE(*,*) 'EMISSIONS IN LAGRANGIAN (ATTILA) SPACE: ON'
       IF (l_lg_tend) THEN
          WRITE(*,*) 'LG-TRACER TENDENCY (ATTILA) IN CHANNEL: ON'
       ELSE
          WRITE(*,*) 'LG-TRACER TENDENCY (ATTILA) IN CHANNEL: OFF'
       END IF
    ELSE
       WRITE(*,*) 'EMISSIONS IN LAGRANGIAN (ATTILA) SPACE: OFF'
    END IF
!!#D attila -

!!#D clams +
    IF (L_CL) THEN
       WRITE(*,*) 'EMISSIONS IN LAGRANGIAN (CLaMS) SPACE: ON'
       IF (l_lg_tend) THEN
          WRITE(*,*) 'LG-TRACER TENDENCY (CLaMS) IN CHANNEL: ON'
       ELSE
          WRITE(*,*) 'LG-TRACER TENDENCY (CLaMS) IN CHANNEL: OFF'
       END IF
    ELSE
       WRITE(*,*) 'EMISSIONS IN LAGRANGIAN (CLaMS) SPACE: OFF'
    END IF
!!#D clams -

    CALL read_nml_close(substr, iou, modstr)

    status = 0  ! no ERROR

  END SUBROUTINE offemis_read_nml_cpl
! ----------------------------------------------------------------------

!*************************************************************************
END MODULE messy_offemis_si
!*************************************************************************

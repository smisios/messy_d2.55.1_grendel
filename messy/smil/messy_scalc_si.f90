#include "messy_main_ppd_bi.inc"
! **********************************************************************
MODULE messy_scalc_si
  ! **********************************************************************
  ! Doing simple mathematical operations on channel objects.
  ! Logic borrowed from messy_sorbit_si.f90
  !
  ! Author: Bastian Kern, MPIC, 08.07.2011
  !         Astrid Kerkweg, Uni Bonn, 2018 / 2019
  !            (Expansion to patches, entrypoint flexibility)
  !         Mariano Mertens, DLR, 2020 (additional operators)
  !         Patrick Joeckel, DLR, 2020 additional DELTA operator
  !                                    - operator with parameters
  !                                    - multiple resulting objects

  USE messy_main_blather_bi,            ONLY: start_message_bi, end_message_bi
  USE messy_main_channel_bi,            ONLY: n_dom
  USE messy_main_import_grid_tools_bi,  ONLY: RGTMAXACTSTR
  USE messy_main_tools,                 ONLY: PTR_3D_ARRAY, PTR_2D_ARRAY
  USE messy_main_channel,               ONLY: STRLEN_OBJECT, STRLEN_CHANNEL
  USE messy_main_constants_mem,         ONLY: STRLEN_ULONG, INT_UNDEF &
                                            , STRLEN_MEDIUM
  USE messy_main_channel_repr,          ONLY: REPR_UNDEF
  USE messy_scalc

  IMPLICIT NONE
  PRIVATE

  INTRINSIC :: NULL

  ! Mathematical functions supported
  INTEGER, PARAMETER              :: OSUM    = 1 ! ordinary sum, needs to be distinguished from INTRINSIC :: SUM
  INTEGER, PARAMETER              :: PSUM    = 2
  INTEGER, PARAMETER              :: PSUMDIA = 3
  INTEGER, PARAMETER              :: SUMGE0  = 4
  INTEGER, PARAMETER              :: MULT    = 5
  INTEGER, PARAMETER              :: DIV     = 6
  INTEGER, PARAMETER              :: QSUM    = 7
  INTEGER, PARAMETER              :: FLAG    = 8
  INTEGER, PARAMETER              :: FLAGU   = 9
  INTEGER, PARAMETER              :: FLAGD   = 10
  INTEGER, PARAMETER              :: DELTA   = 11
  INTEGER, PARAMETER              :: WVA     = 12
  INTEGER, PARAMETER              :: SPLITV  = 13

  ! if we have any PSUMDIA calculation, ldiag will be .TRUE. ->
  ! create diagnostic SCALC output and allocate diag_ptr
  LOGICAL, SAVE :: ldiag = .FALSE.

  TYPE T_CALC_IO
     CHARACTER(LEN=STRLEN_OBJECT) :: object_name = ''  ! object to store result
     CHARACTER(LEN=STRLEN)        :: str = ''          ! channel/object string
     CHARACTER(LEN=STRLEN_MEDIUM) :: math_func_str = ''! math. funct. to use
     CHARACTER(LEN=STRLEN_ULONG)  :: when = 'messy_global_end'
  END TYPE T_CALC_IO

  TYPE T_CALC_CPL
     CHARACTER(LEN=STRLEN_OBJECT)   :: object_name = '' ! object to store result
     CHARACTER(LEN=STRLEN)          :: str = ''         ! channel/object string
     CHARACTER(LEN=STRLEN_MEDIUM)   :: math_func_str = ''! math. funct. to use
     INTEGER                        :: nwhen = 0
     INTEGER, DIMENSION(:), POINTER :: when  => NULL()
     INTEGER, DIMENSION(:), POINTER :: where => NULL()
     CHARACTER(LEN=32)              :: domain_name   = ''
     INTEGER                        :: domain_idx = INT_UNDEF
  END TYPE T_CALC_CPL

  TYPE T_CALC
     TYPE(T_CALC_CPL)             :: io
     INTEGER                      :: nobjnml
     INTEGER                      :: nobj
     LOGICAL                      :: ok
     INTEGER                      :: math_func         ! math. funct. to use
     CHARACTER(LEN=STRLEN_CHANNEL), DIMENSION(:), POINTER :: cha      => NULL()
     CHARACTER(LEN=STRLEN_OBJECT),  DIMENSION(:), POINTER :: obj      => NULL()
     REAL(DP),                      DIMENSION(:), POINTER :: tmp_fac  => NULL()
     LOGICAL,                       DIMENSION(:), POINTER :: lex      => NULL()
     ! POINTER TO DATA
     TYPE(PTR_2D_ARRAY),            DIMENSION(:), POINTER :: dat2d    => NULL()
     TYPE(PTR_3D_ARRAY),            DIMENSION(:), POINTER :: dat3d    => NULL()
     ! rank of data
     INTEGER,                       DIMENSION(:), POINTER :: rank     => NULL()
     ! rank of resulting object
     INTEGER                                              :: res_rank = 0
     ! REPRESENTATION ID
     INTEGER,                       DIMENSION(:), POINTER :: rid      => NULL()
     ! of resulting object
     INTEGER                                            :: res_rid = REPR_UNDEF
     REAL(DP),                      DIMENSION(:), POINTER :: fac      => NULL()
     INTEGER,                       DIMENSION(:), POINTER :: tmp_lev  => NULL()
     INTEGER,                       DIMENSION(:), POINTER :: lev      => NULL()
     ! number of resulting channel objects
     INTEGER                                              :: nres      = 1
     ! number of operator parameters
     INTEGER                                              :: nparam    = 0
     REAL(DP),                      DIMENSION(:), POINTER :: xparam    => NULL()
     TYPE(PTR_2D_ARRAY),            DIMENSION(:), POINTER :: result2d  => NULL()
     TYPE(PTR_3D_ARRAY),            DIMENSION(:), POINTER :: result3d  => NULL()
     ! pointer to 2nd data array for diagnostic purposes
     TYPE(PTR_3D_ARRAY),                          POINTER :: diag     => NULL()
     ! unit string
     CHARACTER(LEN=STRLEN_ULONG), DIMENSION(:),   POINTER :: unit     => NULL()
     ! GEOMETRY in SPACE 2D, 3D-N or 3D-Z
     INTEGER                                              :: spgeo = INT_UNDEF
     ! INDICES FOR GLOBAL FIELD ACCESS in scalc_integrate
     INTEGER, DIMENSION(2)                                 :: xind = INT_UNDEF
     INTEGER, DIMENSION(2)                                 :: yind = INT_UNDEF
  END TYPE T_CALC

  INTEGER, PARAMETER                           :: NMAXCALC = 200
  TYPE(T_CALC_IO), DIMENSION(NMAXCALC),   SAVE :: CALC ! CPL namelist
  TYPE(T_CALC),    DIMENSION(:), POINTER, SAVE :: XCALC => NULL()
  INTEGER,                                SAVE :: NCALC

  PUBLIC :: scalc_initialize
  PUBLIC :: scalc_init_memory
  PUBLIC :: scalc_init_coupler
  PUBLIC :: scalc_init_coupling
  PUBLIC :: scalc_global_start
  PUBLIC :: scalc_convec
  PUBLIC :: scalc_physc
  PUBLIC :: scalc_global_end
  !PRIVATE :: scalc_integrate    ! (formerly scalc_global_end)
  !PRIVATE :: scalc_read_nml_cpl
  PUBLIC :: scalc_free_memory
  !PRIVATE :: parse_when
  !PRIVATE :: get_param

CONTAINS

  ! ======================================================================
  ! PUBLIC ROUTINES
  ! ======================================================================

  ! ----------------------------------------------------------------------
  SUBROUTINE scalc_initialize

    ! SCALC MODULE ROUTINE (ECHAM5 INTERFACE)
    !
    ! read namelist, store and broadcast results
    !
    ! Autor: Bastian Kern, MPIC, 08.07.2011
    !
    USE messy_main_mpi_bi,    ONLY: p_parallel_io, p_io, p_bcast
    USE messy_main_blather_bi,ONLY: error_bi
    USE messy_main_tools,     ONLY: find_next_free_unit, ucase, str &
                                  , domains_from_string
    USE messy_main_control,   ONLY: entry_name

    IMPLICIT NONE

    INTRINSIC :: TRIM

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'scalc_initialize'
    INTEGER                     :: iou    ! I/O unit
    INTEGER                     :: status ! error status
    INTEGER                     :: i, j
    CHARACTER(LEN=32)                        :: varname = '  '
    CHARACTER(LEN=3),  DIMENSION(:), POINTER :: domname => NULL()
    INTEGER,           DIMENSION(:), POINTER :: domnum  => NULL()
    INTEGER          :: nd, num
    CHARACTER(LEN=STRLEN_MEDIUM) :: sopt

    ! initialize CPL
    IF (p_parallel_io) THEN
       iou = find_next_free_unit(100,200)
       CALL scalc_read_nml_cpl(status, iou)
       IF (status /= 0) CALL error_bi('', substr)
    END IF

    CALL start_message_bi(modstr, 'INITIALISATION', substr)

    IF (p_parallel_io) THEN
       NCALC = 1
       DO i=1, NMAXCALC
          IF (TRIM(CALC(i)%object_name) == '') CYCLE

          CALL domains_from_string(status,CALC(i)%object_name,n_dom,num)
          NCALC = NCALC + num
       END DO
       NCALC = NCALC - 1
    END IF
    CALL p_bcast(NCALC, p_io)
    ALLOCATE(XCALC(NCALC))

    IF (p_parallel_io) THEN
       NCALC = 1
       DO i=1, NMAXCALC
          IF (TRIM(CALC(i)%object_name) == '') CYCLE

          CALL domains_from_string(status,CALC(i)%object_name,n_dom,num &
               ,varname,domname,domnum)
          IF (status /= 0) CALL error_bi('error in namelist parsing', substr)

          DO nd = 1, SIZE(domname)
             XCALC(NCALC)%io%object_name = TRIM(varname)
             XCALC(NCALC)%io%domain_idx  = domnum(nd)
             XCALC(NCALC)%io%domain_name = TRIM(domname(nd))
             XCALC(NCALC)%io%str = TRIM(CALC(i)%str)
             CALL ucase(CALC(i)%math_func_str)
             XCALC(NCALC)%io%math_func_str = TRIM(CALC(i)%math_func_str)

             SELECT CASE (TRIM(XCALC(NCALC)%io%math_func_str))
             CASE ('SUM')
                XCALC(NCALC)%math_func = OSUM
             CASE ('PSUM')
                XCALC(NCALC)%math_func = PSUM
             CASE ('PSUMDIA')
                XCALC(NCALC)%math_func = PSUMDIA
                ldiag = .TRUE.
             CASE ('SUMGE0')
                XCALC(NCALC)%math_func = SUMGE0
             CASE ('MULT')
                XCALC(NCALC)%math_func = MULT
             CASE ('DIV')
                XCALC(NCALC)%math_func = DIV
             CASE ('QSUM')
                XCALC(NCALC)%math_func = QSUM
             CASE ('FLAG')
                XCALC(NCALC)%math_func = FLAG
             CASE ('FLAGU')
                XCALC(NCALC)%math_func = FLAGU
             CASE ('FLAGD')
                XCALC(NCALC)%math_func = FLAGD
             CASE ('WVA')
                XCALC(NCALC)%math_func = WVA
             CASE ('SPLITV')
                XCALC(NCALC)%math_func = SPLITV
                ! default if no parameter is given
                XCALC(NCALC)%nparam = 1
                ALLOCATE(XCALC(NCALC)%xparam(1))
                XCALC(NCALC)%xparam = 0.0
                XCALC(NCALC)%nres = 2 ! two resulting channel objects

             CASE DEFAULT
                ! parse operators with parameters, e.g. DELTA(a,b,...)
                CALL get_param(status, &
                     TRIM(XCALC(NCALC)%io%math_func_str) &
                     , sopt                              &
                     , XCALC(NCALC)%xparam               &
                     , XCALC(NCALC)%nparam)
                IF (status /= 0) &
                     CALL error_bi(&
                     'error in parsing operator for parameters OR '//&
                     &'unsupported mathematical function in CALC ('&
                     //TRIM(str(i))//') in '//TRIM(modstr)//'.nml' &
                     , substr)
                !
                SELECT CASE(TRIM(sopt))
                CASE ('DELTA')
                   IF (XCALC(NCALC)%nparam /= 2) &
                        CALL error_bi(&
                        'DELTA operator requires 2 parameters', substr)
                   XCALC(NCALC)%math_func = DELTA
                   XCALC(NCALC)%nres = 2 ! two resulting channel objects
                CASE ('SPLITV')
                   IF (XCALC(NCALC)%nparam > 1) &
                        CALL error_bi(&
                        'SPLITV operator accepts only one parameter', substr)
                   XCALC(NCALC)%math_func = SPLITV
                   XCALC(NCALC)%nres = 2 ! two resulting channel objects
                CASE DEFAULT
                   CALL error_bi('Unsupported mathematical function in CALC ('&
                        //TRIM(str(i))//') in '//TRIM(modstr)//'.nml', substr)
                END SELECT
                !
             END SELECT

             CALL parse_when(status,CALC(i)%when, XCALC(NCALC)%io%nwhen &
                  , XCALC(NCALC)%io%when, XCALC(NCALC)%io%where)
             IF (status /= 0) CALL error_bi('ERROR IN STRING PARSING!' , substr)
             WRITE(*,*) 'scalc: calculate ' &
                  ,TRIM(XCALC(NCALC)%io%object_name),' in:'
             DO j = 1, XCALC(NCALC)%io%nwhen
                WRITE (*,*) '                ...' &
                     ,entry_name(XCALC(NCALC)%io%when(j)) &
                     ,XCALC(NCALC)%io%where(j)
             END DO

             ! set preliminary %nobjnml object number in namelist
             ! (incl. wildcards)
             CALL str2chobfac(status, XCALC(NCALC)%io%str&
                  , XCALC(NCALC)%nobjnml  &
                  , XCALC(NCALC)%cha, XCALC(NCALC)%obj, XCALC(NCALC)%tmp_fac &
                  , XCALC(NCALC)%tmp_lev)
             IF (status /= 0) THEN
                CALL error_bi('Error while reading object string of CALC ('&
                     //TRIM(str(i))//') in '//TRIM(modstr)//'.nml', substr)
             END IF

             ! SOME consistency checks
             SELECT CASE(XCALC(NCALC)%math_func)
                !
             CASE(SUMGE0,FLAG,FLAGU,FLAGD,WVA)
                ! SUMGE0 should only be used with exactly
                ! 2 objects to be summed up
                IF (XCALC(NCALC)%nobjnml /= 2) THEN
                   CALL error_bi(&
                        'operator used in CALC ('//&
                        TRIM(str(i))//') of '//TRIM(modstr)//'.nml'//&
                        ' is only defined for exactly 2 objects' &
                        , substr)
                END IF
                !
             CASE(DELTA)
                IF (XCALC(NCALC)%nobjnml /= 1) THEN
                   CALL error_bi(&
                        'operator used in CALC ('//&
                        TRIM(str(i))//') of '//TRIM(modstr)//'.nml'//&
                        ' is only defined for exactly 1 object' &
                        , substr)
                END IF
                !
             CASE DEFAULT
                IF (XCALC(NCALC)%nobjnml == 0) THEN
                   CALL error_bi('No objects found in CALC ('&
                        //TRIM(str(i))//') of '//TRIM(modstr)//'.nml', substr)
                END IF
                !
             END SELECT

             ! SOME more consistency checks
             SELECT CASE(XCALC(NCALC)%math_func)
                !
             CASE(FLAG,FLAGU,FLAGD,DELTA)
                !
                DO j=1, XCALC(NCALC)%nobjnml
                   IF (XCALC(NCALC)%tmp_fac(j) /= 1.0_dp) THEN
                      CALL error_bi(&
                           'objets of operator used in CALC ('//&
                           TRIM(str(i))//') of '//TRIM(modstr)//'.nml'//&
                           ' must not be scaled' &
                           , substr)
                   END IF
                END DO
                !
             CASE DEFAULT
                ! nothing to do
             END SELECT

             WRITE(*,*) 'scalc:', TRIM(XCALC(NCALC)%io%object_name) &
                  ,'  Number of channels : ', XCALC(NCALC)%nobjnml  &
                  ,'  Method : ', TRIM(XCALC(NCALC)%io%math_func_str)

             WRITE(*,'(1x,a10,1x,a26,1x,a14,1x,a14)') 'CHANNEL','OBJECT(S)' &
                  ,'FACTOR','LEVEL'
             WRITE(*,'(1x,a10,1x,a26,1x,a14,1x,a14)') '-------','---------' &
                  ,'------','------'
             DO j=1, XCALC(NCALC)%nobjnml
                WRITE(*,'(1x,a10,1x,a26,1x,e20.10,1x,i5)') &
                     TRIM(XCALC(NCALC)%cha(j))   &
                     , TRIM(XCALC(NCALC)%obj(j)), XCALC(NCALC)%tmp_fac(j) &
                     , XCALC(NCALC)%tmp_lev(j)
             END DO

             ! next CALC
             NCALC = NCALC + 1
             WRITE(*,*) '------------------------------------------------------'

          END DO

          DEALLOCATE(domname, domnum)
          NULLIFY(domname)
          NULLIFY(domnum)

       END DO ! NMAXCALC
       NCALC = NCALC - 1
    END IF

    IF (NCALC /= SIZE(XCALC)) &
         CALL error_bi('error in namelist interpretation', substr)
    CALL p_bcast(ldiag, p_io)

    ! broadcast all results
    DO i=1, NCALC
       ! I/O user interface
       CALL p_bcast(XCALC(i)%io%object_name, p_io)
       CALL p_bcast(XCALC(i)%io%domain_idx,  p_io)
       CALL p_bcast(XCALC(i)%io%domain_name, p_io)
       CALL p_bcast(XCALC(i)%io%str,  p_io)
       CALL p_bcast(XCALC(i)%io%math_func_str,  p_io)
       CALL p_bcast(XCALC(i)%io%nwhen,  p_io)
       IF (.NOT. p_parallel_io) THEN
          ALLOCATE(XCALC(i)%io%when(XCALC(i)%io%nwhen))
          ALLOCATE(XCALC(i)%io%where(XCALC(i)%io%nwhen))
       END IF
       DO j = 1, XCALC(i)%io%nwhen
          CALL p_bcast(XCALC(i)%io%when(j),  p_io)
          CALL p_bcast(XCALC(i)%io%where(j), p_io)
       END DO
       !
       ! parameters
       CALL p_bcast(XCALC(i)%nparam, p_io)
       IF (XCALC(i)%nparam > 0) THEN
          IF (.NOT. p_parallel_io) THEN
             ALLOCATE(XCALC(i)%xparam(XCALC(i)%nparam))
          END IF
          CALL p_bcast(XCALC(i)%xparam(:), p_io)
       END IF
       ! number of results
       CALL p_bcast(XCALC(i)%nres, p_io)
       !
       ! channel/object names
       CALL p_bcast(XCALC(i)%nobjnml, p_io)
       CALL p_bcast(XCALC(i)%ok, p_io)
       CALL p_bcast(XCALC(i)%math_func, p_io)
       IF (.NOT. p_parallel_io) THEN
          ALLOCATE(XCALC(i)%cha(XCALC(i)%nobjnml))
          ALLOCATE(XCALC(i)%obj(XCALC(i)%nobjnml))
          ALLOCATE(XCALC(i)%tmp_fac(XCALC(i)%nobjnml))
          ALLOCATE(XCALC(i)%tmp_lev(XCALC(i)%nobjnml))
       END IF
       DO j=1, XCALC(i)%nobjnml
          CALL p_bcast(XCALC(i)%cha(j), p_io)
          CALL p_bcast(XCALC(i)%obj(j), p_io)
          CALL p_bcast(XCALC(i)%tmp_fac(j), p_io)
          CALL p_bcast(XCALC(i)%tmp_lev(j), p_io)
       END DO
       !
    END DO

    IF (p_parallel_io) THEN
       WRITE(*,*) ' ---> ',NCALC,'SCALC CHANNEL OBJECTS INITIALIZED !'
    END IF

    CALL end_message_bi(modstr,'INITIALISATION ',substr)

  END SUBROUTINE scalc_initialize
  ! ----------------------------------------------------------------------

  ! ----------------------------------------------------------------------
  SUBROUTINE scalc_init_memory(l_cont)

    ! currently only required for ECHAM

    IMPLICIT NONE

    LOGICAL, OPTIONAL, INTENT(in) :: l_cont

    CALL scalc_init_coupling(l_cont)

  END SUBROUTINE scalc_init_memory
  ! ----------------------------------------------------------------------

  ! ----------------------------------------------------------------------
  SUBROUTINE scalc_init_coupler(l_cont)

    ! currently only required for COSMO

    IMPLICIT NONE

    LOGICAL, OPTIONAL, INTENT(in) :: l_cont

    CALL scalc_init_coupling(l_cont)

  END SUBROUTINE scalc_init_coupler
  ! ---------------------------------------------------------------------

  ! ---------------------------------------------------------------------
  SUBROUTINE scalc_init_coupling(l_cont)

    ! SCALC MODULE ROUTINE
    !
    ! channel/object definitions
    !
    ! Author: Bastian Kern, MPIC, 08.07.2011
    !
    USE messy_main_blather_bi,       ONLY: error_bi, info_bi, warning_bi
    USE messy_main_channel_error_bi, ONLY: channel_halt
    USE messy_main_channel_bi,       ONLY: GP_3D_1LEV, GP_2D_HORIZONTAL
    USE messy_main_mpi_bi,           ONLY: p_parallel_io
    USE messy_main_constants_mem,    ONLY: STRLEN_ULONG
    USE messy_main_channel,          ONLY: new_channel, new_channel_object    &
                                         , new_attribute, get_channel_object  &
                                         , get_channel_info, get_attribute    &
                                         , get_channel_object_info            &
                                         , get_channel_object_slice           &
                                         , get_channel_info
    USE messy_main_channel_mem,      ONLY: dom_curid
    USE messy_main_channel_repr,     ONLY: REPR_UNDEF,get_representation_info &
                                         , IRANK,repr_rank_location
    USE messy_main_tools,            ONLY: match_wild, int2str


    IMPLICIT NONE

    INTRINSIC :: ASSOCIATED, SIZE, TRIM

    LOGICAL, INTENT(IN), OPTIONAL :: l_cont

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER  :: substr = 'scalc_init_coupling'
    INTEGER                      :: status
    INTEGER                      :: old_reprid
    INTEGER                      :: new_reprid
    INTEGER                      :: i, j, nptr, jo
    CHARACTER(LEN=STRLEN_OBJECT), DIMENSION(:), POINTER :: ONAMES => NULL()
    CHARACTER(LEN=STRLEN_ULONG)  :: unit_string, new_unit, old_unit
    LOGICAL                      :: lfound
    INTEGER                      :: old_rank, new_rank
    CHARACTER(LEN=1)             :: third_axis
    CHARACTER(LEN=4)             :: old_axis
    LOGICAL                      :: lfirst = .FALSE.
    LOGICAL                      :: lcycle
    ! continue, do not interrupt simulation as scalc_init_coupling
    ! will be called again and objects might not be defined at the
    ! previous calls
    LOGICAL                      :: llcont = .FALSE.
    LOGICAL                      :: lnew  ! found new data objects reallocate
    INTEGER                      :: size1_2d, size2_2d, size1_3d, size2_3d
    INTEGER                      :: rankidx(IRANK)
    INTEGER                      :: nold
    INTEGER                      :: ires
    CHARACTER(LEN=3)             :: suffix
    REAL(dp), DIMENSION(:,:),   POINTER :: ptr2d
    REAL(dp), DIMENSION(:,:,:), POINTER :: ptr3d

    IF (PRESENT(l_cont)) THEN
       llcont = l_cont
    ELSE
       llcont = .FALSE. ! force simulation interrupt at error
    END IF

    CALL start_message_bi(modstr, 'CHANNEL DEFINITION', substr)

    ! create new channel 'scalc'
    IF (NCALC /= 0) THEN
       ! scalc_init_coupling is called several times, make sure that
       ! channel is only created once
       CALL get_channel_info(status, modstr)
       lfirst = (status /= 0)
       IF (lfirst) THEN
          CALL new_channel(status, modstr)
          CALL channel_halt(substr, status)
          ! create new 'scalc_diag' channel, if any diagnostic is requested
          IF (ldiag) THEN
             CALL new_channel(status, modstr//'_diag')
             CALL channel_halt(substr, status)
          END IF
       END IF
    END IF

    calc_loop: DO i=1, NCALC
       IF ( XCALC(i)%io%domain_idx /= dom_curid) CYCLE
       lcycle = .FALSE.
       nptr = 0
       object_loop_1: DO j=1, XCALC(i)%nobjnml

          ! First check if channels/objects are available.
          ! As long as l_cont = true only a warning is issued.
          ! If not present(l_cont) error_bi is called if objects/channels
          ! are not available.
          CALL get_channel_info(status, TRIM(XCALC(i)%cha(j)), ONAMES=ONAMES)
          IF (status /= 3003) THEN            ! channel (name) does not exist
             CALL channel_halt(substr, status)
          ELSE
             IF (llcont) THEN
                CALL warning_bi(' Channel ''' //TRIM(XCALC(i)%cha(j))// &
                     ''' in scalc.nml does not exist! (1a) ', substr)
#ifdef ICON
                CALL warning_bi('ERROR occured for domain: ''' &
                     //TRIM(XCALC(i)%io%domain_name), substr)
#endif
                lcycle = .TRUE.
                CYCLE
             ELSE
                CALL error_bi(' Channel ''' //TRIM(XCALC(i)%cha(j))// &
                     ''' in scalc.nml does not exist! (1b) ', substr)
             END IF
          END IF
          nold = nptr
          DO jo=1, SIZE(ONAMES)
             IF (match_wild(TRIM(XCALC(i)%obj(j)), TRIM(ONAMES(jo))) ) THEN
                nptr = nptr + 1
             END IF
          END DO
          IF (.NOT. llcont .AND. nptr == nold) CALL error_bi(&
               ' No channel object found for '//TRIM(XCALC(i)%obj(j))// &
               &' in channel '//TRIM(XCALC(i)%cha(j))//' !', substr)
       END DO object_loop_1 ! XCALC(i)%nobjnml
       !
       IF (lcycle) CYCLE

       IF (nptr == 0 .AND. .NOT. llcont) &
            CALL error_bi( ' No Channel object found for''' //&
            &TRIM(XCALC(i)%cha(j))// &
            ''' in scalc.nml does not exist! (1c) ', substr)

       IF (ASSOCIATED(ONAMES)) THEN
          DEALLOCATE(ONAMES)
          NULLIFY(ONAMES)
       END IF

       ! Test if field have already been allocated
       lnew = .FALSE.
       IF (.NOT. lfirst) THEN
          IF (ASSOCIATED(XCALC(i)%dat2d)) THEN
             ! is allocated size equal to number of found objects ?
             IF (SIZE(XCALC(i)%dat2d) /= nptr) THEN
                ! deallocate and reallocate data objects
                write (*,*) ' NEW NUMBER OF OBJECT: reallocate for ' &
                     , TRIM(XCALC(i)%io%object_name) &
                     , ' OLD SIZE: ',SIZE(XCALC(i)%dat2d) &
                     , ' NEW SIZE: ', nptr
                lnew = .TRUE.  ! reallocate
                CALL scalc_freemem(i)
             END IF
          ELSE
             lnew = .TRUE.
          END IF
       END IF

       ! Do not analyse namelist again.
       ! Due to the CYCLE below, the lower part (allocate ff)
       ! is not executed if not lnow/lfirst.
       IF (lnew .OR. lfirst) THEN
          IF (p_parallel_io) &
               write (*,*) ' (RE-)ALLOCATE OBJECT: ' &
               , TRIM(XCALC(i)%io%object_name)       &
               , '   to ', nptr
       ELSE
          IF (p_parallel_io) &
               write (*,*) ' NO (RE-)ALLOCATE REQUIRED FOR OBJECT: ' &
               ,TRIM(XCALC(i)%io%object_name) ,'  .... CYCLE'
          CYCLE
       END IF

       ! allocate space for pointer arrays (data, repr-id, object,
       ! result, scaling factors)
       ALLOCATE(XCALC(i)%dat2d(nptr))   ! pointer to channel object
       ALLOCATE(XCALC(i)%dat3d(nptr))   ! pointer to channel object

       ALLOCATE(XCALC(i)%lex(nptr))   ! channel / object exists?
       XCALC(i)%lex(:) = .FALSE.
       ALLOCATE(XCALC(i)%rid(nptr))   ! representation of channel object
       XCALC(i)%rid(:) = REPR_UNDEF
       ALLOCATE(XCALC(i)%result2d(XCALC(i)%nres))    ! array for data output
       ALLOCATE(XCALC(i)%result3d(XCALC(i)%nres))    ! array for data output
       ALLOCATE(XCALC(i)%fac(nptr))   ! factor
       XCALC(i)%fac=1.0_dp            ! init
       ALLOCATE(XCALC(i)%rank(nptr))  ! rank of channel object
       XCALC(i)%rank(:) = INT_UNDEF
       ALLOCATE(XCALC(i)%lev(nptr))   ! level
       XCALC(i)%lev(:)  = INT_UNDEF

       ! in case of using sum with diagnostics (PSUMDIA),
       ! allocate array for diagnostics
       IF (XCALC(i)%math_func == PSUMDIA) &
            & ALLOCATE(XCALC(i)%diag)

       XCALC(i)%nobj = nptr           ! number of total objects
       ! (after wildcard matching)

       ALLOCATE(XCALC(i)%unit(nptr))

       nptr = 0
       third_axis ='-'

       request_loop: DO j=1, XCALC(i)%nobjnml

          CALL get_channel_info(status, TRIM(XCALC(i)%cha(j)), ONAMES = ONAMES)
          IF (status == 3003) THEN    ! channel (name) does not exist
             CALL error_bi(' Channel ''' //TRIM(XCALC(i)%cha(j))// &
                  ''' in scalc.nml does not exist! (2)', substr)
          ELSE
             CALL channel_halt(substr,status)
          END IF

          lfound = .FALSE.

          object_loop_2: DO jo=1, SIZE(ONAMES)
             ! check if potential object is matching
             IF (match_wild(TRIM(XCALC(i)%obj(j)), TRIM(ONAMES(jo)))) THEN
                nptr = nptr + 1
                lfound = .TRUE.
             ELSE
                CYCLE
             END IF

             CALL info_bi('          '//TRIM(XCALC(i)%cha(j))//' - '//&
                  &TRIM(ONAMES(jo))&
#ifdef ICON
                  &//' - '//TRIM(XCALC(i)%io%domain_name) &
#endif
                  , substr)

             NULLIFY(XCALC(i)%dat2d(nptr)%ptr)
             NULLIFY(XCALC(i)%dat3d(nptr)%ptr)

             CALL get_channel_object_info(status                              &
                  , TRIM(XCALC(i)%cha(j)), TRIM(ONAMES(jo))                   &
                  , reprid =old_reprid)
             CALL channel_halt(substr,status)
             CALL get_representation_info(status, ' '                         &
                  , ID =old_reprid, rank = old_rank, axis = old_axis)
             CALL channel_halt(substr,status)

             IF (old_rank == 3) THEN
                IF (XCALC(i)%tmp_lev(j) < 0) THEN
                   CALL get_channel_object(status                             &
                        , TRIM(XCALC(i)%cha(j)), TRIM(ONAMES(jo))             &
                        , p3=XCALC(i)%dat3d(nptr)%ptr)
                   CALL channel_halt(substr,status)
                ELSE
                   CALL get_channel_object_slice(status                       &
                        , TRIM(XCALC(i)%cha(j)), TRIM(ONAMES(jo))             &
                        , p2=XCALC(i)%dat2d(nptr)%ptr                         &
                        , zslice=XCALC(i)%tmp_lev(j))
                   CALL channel_halt(substr,status)
                   old_rank   = old_rank -1
                   old_reprid = GP_2D_HORIZONTAL
                   old_axis   = 'XY--'
                END IF
             ELSE
                CALL get_channel_object(status                             &
                     , TRIM(XCALC(i)%cha(j)), TRIM(ONAMES(jo))             &
                     , p2=XCALC(i)%dat2d(nptr)%ptr)
             END IF

             ! now check representation:
             ! For correct application of Rank Identifiers
             ! we need to store, if the third rank (for rank = 3)
             ! is N or Z.
             ! Further, it is explicitly checked if the other two dimensions
             ! are X and Y (also for rank 2).
             ! So far, only XY? fields are supported.
             ! The function returns a rankidx with dim = 4,
             ! rankidx(3) gives the dimension-nr of the z-axis,
             ! rankidx(4) gives the dimension-nr of the n-axis.
             ! If z/n axis are not available rankidx(3/4) = INIT_UNDEF.
             ! The same applies, if field has no X or Y axis.
             call repr_rank_location(old_axis,rankidx)

             ! first check, if first two dimensions of rank 2 are X and Y
             IF (old_rank == 2) THEN
                IF ((rankidx(1) == INT_UNDEF) .OR.   &
                     (rankidx(2) == INT_UNDEF)) THEN
                   ! rank 2 without X/Y
                   CALL error_bi(' Channel ''' //TRIM(XCALC(i)%cha(j))//&
                        &':'//TRIM(XCALC(i)%obj(j))//&
#ifdef ICON
                        &''' - '//TRIM(XCALC(i)%io%domain_name)//&
#endif

                        &'''  is rank 2 but has no X or Y axis', substr)
                END IF
             END IF

             !third_axis = '-'
             IF (old_rank == 3) THEN
                IF ((rankidx(1) == INT_UNDEF) .OR.   &
                     (rankidx(2) == INT_UNDEF)) THEN
                   ! rank 3 without X/Y
                   CALL error_bi(' Channel ''' //TRIM(XCALC(i)%cha(j))// &
                        &':'//TRIM(XCALC(i)%obj(j))//&
                        &'''  is rank 3 but has no X or Y axis', substr)
                ELSE
                   third_axis ='-'
                   ! Check what the third axis is
                   ! no z-axis
                   IF (rankidx(3) == INT_UNDEF) third_axis = 'N'
                   ! no n-axis
                   IF (rankidx(4) == INT_UNDEF) third_axis = 'Z'

                   IF (third_axis == '') THEN
                      ! rank 3 and thrid axis is neiter N nor Z
                      CALL error_bi(' Channel ''' //&
                           &TRIM(XCALC(i)%cha(j))//&
                           &':'//TRIM(XCALC(i)%obj(j))//&
#ifdef ICON
                           &''' - '//TRIM(XCALC(i)%io%domain_name)//&
#endif
                           &'''  is rank 3 but third axis is neither'//&
                           &' N nor Z', substr)
                   END IF

                END IF
             END IF

             ! set new object's representation
             XCALC(i)%rank(nptr) = old_rank
             XCALC(i)%rid(nptr)  = old_reprid
             XCALC(i)%lex(nptr)  = .TRUE.
             ! associate right scaling factor to (wildcard expanded) object
             XCALC(i)%fac(nptr) = XCALC(i)%tmp_fac(j)
             ! vertical level index (for reduced 3D fields)
             XCALC(i)%lev(nptr) = XCALC(i)%tmp_lev(j)

             ! get unit of channel object
             ! units are checked later on
             CALL get_attribute(status  &
                  , TRIM(XCALC(i)%cha(j)), TRIM(ONAMES(jo)) &
                  , 'units', c=unit_string                  &
                  )
             IF(status /= 0) THEN
                CALL warning_bi('object '//&
                     &TRIM(XCALC(i)%cha(j))//' in channel '//&
                     TRIM(ONAMES(jo))//&
#ifdef ICON
                     &' in domain '//TRIM(XCALC(i)%io%domain_name)//&
#endif
                     &' has no units attribute', substr)
                unit_string = 'unknown'
             END IF

             XCALC(i)%unit(nptr) = unit_string

          END DO object_loop_2

          ! Test if any channel object found matching the user request
          ! in namelist
          ! This part checks only, if the pointers are available.
          ! The check, if channel objects are available is performed above.

          IF (.NOT. lfound) THEN
             IF (llcont) THEN
                CALL warning_bi(' No channel object found matching '//&
                     &TRIM(XCALC(i)%cha(j))// &
#ifdef ICON
                     &'- '//TRIM(XCALC(i)%io%domain_name)//&
#endif
                     &':'//TRIM(XCALC(i)%obj(j))//''' in scalc.nml!', substr)
                lcycle =.TRUE.
             ELSE
                CALL error_bi(' No channel object found matching '//&
                     &TRIM(XCALC(i)%cha(j))// &
#ifdef ICON
                     &'- '//TRIM(XCALC(i)%io%domain_name)//&
#endif
                     &':'//TRIM(XCALC(i)%obj(j))//''' in scalc.nml!', substr)
             END IF
          END IF

       END DO request_loop
       IF (lcycle) CYCLE

       IF (nptr /= XCALC(i)%nobj) THEN
          CALL error_bi('something went wrong with pointer counting', substr)
       END IF

       IF (ASSOCIATED(ONAMES)) DEALLOCATE(ONAMES)

       ! test if representations are matching ...
       new_reprid = REPR_UNDEF
       old_reprid = REPR_UNDEF

       XCALC(i)%ok = .FALSE.

       object_loop_3: DO jo=1, XCALC(i)%nobj
          IF (jo == 1) THEN
             XCALC(i)%ok = .TRUE.
             old_reprid = XCALC(i)%rid(1)
             old_rank   = XCALC(i)%rank(1)
             CYCLE
          END IF
          new_rank = XCALC(i)%rank(jo)

          SELECT CASE (XCALC(i)%math_func)
          ! IF ( (XCALC(i)%math_func == FLAG) .OR. &
          !      (XCALC(i)%math_func == FLAGU) .OR. &
          !      (XCALC(i)%math_func == FLAGD))  THEN
          CASE(FLAG, FLAGU, FLAGD)
             new_reprid = XCALC(i)%rid(jo)

             ! Note: For FLAG only two objects are possible.
             ! For FLAG the first field can either be GP_3D, GP_2D or NX2D.
             ! The second field (the FLAG field) has to be GP_2D_HORIZONTAL.

             ! Check field which should be flagged
             IF ((old_reprid /= GP_2D_HORIZONTAL) .AND. &
                  (old_reprid /= GP_3D_1LEV )) THEN
                ! is it 3D ?
                IF (old_rank /= 3) THEN
                   CALL warning_bi(&
                        'FIELD TO be flagged is neither GP_2D_HOR nor rank 3 ',&
                        substr)
                   XCALC(i)%ok = .FALSE.
                ELSE
                   ! Check, if size of the two dimensions are
                   ! identical between both fields.
                   size1_2D=size(XCALC(i)%dat2d(2)%ptr,_IX_XY___)
                   size2_2D=size(XCALC(i)%dat2d(2)%ptr,_IY_XY___)

                   size1_3D=size(XCALC(i)%dat3d(1)%ptr,_IX_XYZ__)
                   size2_3D=size(XCALC(i)%dat3d(1)%ptr,_IY_XYZ__)

                   IF (size1_2D/= size1_3D) THEN
                      CALL warning_bi('Warning first horizontal dimensions '//&
                           &'of object 1 and object 2 do not match ', substr)
                      XCALC(i)%ok = .FALSE.
                   ELSE IF (size2_2D/= size2_3D) THEN
                      CALL warning_bi('Warning second horizontal dimensions '//&
                           &'of object 1 and object 2 do not match ', substr)
                      XCALC(i)%ok = .FALSE.
                   END IF
                END IF
             END IF

             ! check FLAG field
             ! Flag field must be 2D, either it GP_2D_Horizontal or GP_3D_1LEV
             IF (new_reprid /= GP_2D_HORIZONTAL) THEN
                IF  (new_reprid /= GP_3D_1LEV ) THEN
                   CALL warning_bi('FLAG field must have representation'//&
                        &' GP_2D_HORIZONTAL or GP_3D_1LEV', substr)
                   XCALC(i)%ok = .FALSE.
                END IF
             END IF

          CASE(SPLITV)

             new_reprid = XCALC(i)%rid(jo)
             IF (new_reprid == GP_2D_HORIZONTAL) THEN
                ! everything is ok
                XCALC(i)%ok = .TRUE.
             ELSE
                CALL warning_bi('SPLITV operator only defined when second field'//&
                     &' has representation GP_2D_HORIZONTAL',substr)
                XCALC(i)%ok = .FALSE.
             END IF

          CASE DEFAULT ! all other cases

             IF (new_rank /= old_rank) THEN
                CALL warning_bi('Ranks do not match!', substr)
                XCALC(i)%ok = .FALSE.
                CYCLE
             ELSE
                new_reprid = XCALC(i)%rid(jo)

                IF (new_reprid /= old_reprid) THEN
                   CALL warning_bi('Representations do not match!', substr)
                   IF (new_reprid == GP_3D_1LEV) THEN
                      IF (old_reprid == GP_2D_HORIZONTAL) THEN
                         CALL warning_bi( &
                              '... that is ok (GP_3D_1LEV / GP_2D_HORIZONTAL)' &
                              , substr)
                         CYCLE
                      END IF
                   END IF
                   IF (old_reprid == GP_3D_1LEV) THEN
                      IF (new_reprid == GP_2D_HORIZONTAL) THEN
                         CALL warning_bi( &
                              '... that is ok (GP_2D_HORIZONTAL / GP_3D_1LEV)' &
                              , substr)
                         CYCLE
                      END IF
                   END IF

                   XCALC(i)%ok = .FALSE.
                   CYCLE
                END IF
                old_reprid = new_reprid
             END IF
          END SELECT
       END DO object_loop_3

       ! - set the rank and representation of the resulting object
       SELECT CASE(XCALC(i)%math_func)
       CASE(WVA)
          !
          XCALC(i)%res_rank = 2
          XCALC(i)%res_rid  = GP_2D_HORIZONTAL
          !
       CASE DEFAULT
          !
          XCALC(i)%res_rank = XCALC(i)%rank(1)
          XCALC(i)%res_rid  = old_reprid
       END SELECT

       ! - check the units
       ! - currently, only a warning is output, if the units are not identical
       new_unit = ''
       old_unit = ''

       SELECT CASE(XCALC(i)%math_func)
       CASE(WVA)
          ! NOTHING TO DO, WEIGHT CAN BE ANYTHING
       CASE DEFAULT
          !
          object_loop_4: DO jo=1, XCALC(i)%nobj
             IF (jo == 1) THEN
                old_unit = XCALC(i)%unit(1)
                CYCLE
             END IF
             new_unit = XCALC(i)%unit(jo)
             IF (new_unit /= old_unit) THEN
                CALL warning_bi('units do not match! The units are '//&
                     &trim(new_unit)//' and '//trim(old_unit), substr)
             END IF
             old_unit = new_unit
          END DO object_loop_4
       END SELECT

       ! if representations match, create new (output) object
       ! qqq At the moment the new object has the unit from the first object
       !     This should be changed in the future!

       is_ok: IF (XCALC(i)%ok) THEN
          results: DO ires=1, XCALC(i)%nres
             NULLIFY(XCALC(i)%result2d(ires)%ptr)
             NULLIFY(XCALC(i)%result3d(ires)%ptr)

             IF (XCALC(i)%nres == 1) THEN
                suffix = ''
             ELSE
                CALL int2str(suffix, ires)
                suffix(1:1) = '_'
             ENDIF

             IF (XCALC(i)%res_rank == 3) THEN
                CALL new_channel_object(status, modstr             &
                     , TRIM(XCALC(i)%io%object_name)//TRIM(suffix) &
                     , p3=XCALC(i)%result3d(ires)%ptr              &
                     , reprid=XCALC(i)%res_rid, lrestreq=.TRUE.    &
                     )
             ELSE
                CALL new_channel_object(status, modstr             &
                     , TRIM(XCALC(i)%io%object_name)//TRIM(suffix) &
                     , p2=XCALC(i)%result2d(ires)%ptr              &
                     , reprid=XCALC(i)%res_rid, lrestreq=.TRUE.    &
                     )
             END IF

             CALL channel_halt(substr, status)
             CALL new_attribute(status, modstr                     &
                  , TRIM(XCALC(i)%io%object_name)//TRIM(suffix)    &
                  , 'math_funct'                                   &
                  , c=TRIM(XCALC(i)%io%math_func_str)              &
                  )
             CALL channel_halt(substr, status)

             CALL new_attribute(status, modstr                     &
                  , TRIM(XCALC(i)%io%object_name)//TRIM(suffix)    &
                  , 'units'                                        &
                  , c=XCALC(i)%unit(1)                             &
                  )
             CALL channel_halt(substr, status)

          END DO results

          ! in case of diagnostics for PSUMDIA, create the diagnostic object in
          ! channel 'scalc_diag'
          IF (XCALC(i)%math_func == PSUMDIA) THEN
             NULLIFY(XCALC(i)%diag%ptr)

             CALL new_channel_object(status, modstr//'_diag'                &
                  , TRIM(XCALC(i)%io%object_name)//'_err'                   &
                  , p3=XCALC(i)%diag%ptr                                    &
                  , reprid=old_reprid, lrestreq=.TRUE.                      &
                  )
             CALL channel_halt(substr, status)
             CALL new_attribute(status, modstr//'_diag'                     &
                  , TRIM(XCALC(i)%io%object_name)//'_err', 'math_funct'     &
                  , c=TRIM(XCALC(i)%io%math_func_str)                       &
                  )
             CALL channel_halt(substr, status)
             CALL new_attribute(status, modstr//'_diag'                     &
                  , TRIM(XCALC(i)%io%object_name)//'_err', 'units'          &
                  , c=XCALC(i)%unit(1)                                      &
                  )
             CALL channel_halt(substr, status)
          END IF

          ! Definition of some specifications used in scalc_integrate

          ! 1. RANK
          ! from here on only the rank of the resulting object(s) matter(s)
          ! set XCALC(.)%spgeo  space geometry to distinguish between
          ! 2D, 3D-'N' and 3D-'Z'.
          IF (XCALC(i)%rank(1) == 2) THEN
             XCALC(i)%spgeo = IS_2D
          ELSE
             IF (third_axis == 'N') XCALC(i)%spgeo = IS_3D_N
             IF (third_axis == 'Z') XCALC(i)%spgeo = IS_3D_Z
          END IF
          DEALLOCATE(XCALC(i)%rank); NULLIFY(XCALC(i)%rank)

          SELECT CASE(XCALC(i)%math_func)
          CASE(WVA)
             !
             IF (XCALC(i)%spgeo == IS_2D) THEN
                CALL error_bi('WVA operator only defined for 3D objects',substr)
             END IF
             !
             ptr2d => XCALC(i)%dat2d(1)%ptr
             ptr3d => XCALC(i)%dat3d(1)%ptr
             !
          CASE(SPLITV)
             !
             IF (XCALC(i)%spgeo /= IS_3D_Z) THEN
                CALL error_bi('SPLITV operator only defined for 3D objects with vertical dimension',substr)
             END IF
             !
             ptr2d => XCALC(i)%dat2d(1)%ptr
             ptr3d => XCALC(i)%dat3d(1)%ptr
             !
          CASE DEFAULT
             !
             ptr2d => XCALC(i)%result2d(1)%ptr
             ptr3d => XCALC(i)%result3d(1)%ptr
             !
          END SELECT

          ! 2. indices for global field access:
          IF (XCALC(i)%spgeo == IS_2D) THEN
             XCALC(i)%xind(1) = LBOUND( ptr2d, _IX_XY___ )
             XCALC(i)%xind(2) = UBOUND( ptr2d, _IX_XY___ )
             XCALC(i)%yind(1) = LBOUND( ptr2d, _IY_XY___ )
             XCALC(i)%yind(2) = UBOUND( ptr2d, _IY_XY___ )
          ELSE IF (XCALC(i)%spgeo == IS_3D_N) THEN
             XCALC(i)%xind(1) = LBOUND( ptr3d, _IX_XY_N_ )
             XCALC(i)%xind(2) = UBOUND( ptr3d, _IX_XY_N_ )
             XCALC(i)%yind(1) = LBOUND( ptr3d, _IY_XY_N_ )
             XCALC(i)%yind(2) = UBOUND( ptr3d, _IY_XY_N_ )
          ELSE IF (XCALC(i)%spgeo == IS_3D_Z) THEN
             XCALC(i)%xind(1) = LBOUND( ptr3d, _IX_XYZ__ )
             XCALC(i)%xind(2) = UBOUND( ptr3d, _IX_XYZ__ )
             XCALC(i)%yind(1) = LBOUND( ptr3d, _IY_XYZ__ )
             XCALC(i)%yind(2) = UBOUND( ptr3d, _IY_XYZ__ )
          END IF

       END IF is_ok

    END DO calc_loop ! NCALC

    CALL end_message_bi(modstr, 'CHANNEL DEFINITION', substr)

  END SUBROUTINE scalc_init_coupling
  ! ---------------------------------------------------------------------

  ! ---------------------------------------------------------------------
  SUBROUTINE scalc_global_start

    IMPLICIT NONE

    CALL scalc_integrate

  END SUBROUTINE scalc_global_start
  ! ---------------------------------------------------------------------

  ! ---------------------------------------------------------------------
  SUBROUTINE scalc_convec(lglobal)

    USE messy_main_grid_def_mem_bi,          ONLY:   jrow, kproma

    IMPLICIT NONE

    LOGICAL, INTENT(IN), OPTIONAL :: lglobal
    !  LOCAL
    LOGICAL :: lglob

    IF (PRESENT(lglobal)) THEN
       lglob = lglobal
    ELSE
       lglob = .FALSE.
    END IF

    IF (lglob) THEN
       CALL scalc_integrate
    ELSE
       CALL scalc_integrate(jrow,kproma)
    END IF

  END SUBROUTINE scalc_convec
  ! ---------------------------------------------------------------------

  ! ---------------------------------------------------------------------
  SUBROUTINE scalc_physc

    USE messy_main_grid_def_mem_bi,          ONLY:   jrow, kproma

    IMPLICIT NONE

    CALL scalc_integrate(jrow,kproma)

  END SUBROUTINE scalc_physc
  ! ---------------------------------------------------------------------

  ! ---------------------------------------------------------------------
  SUBROUTINE scalc_global_end

    IMPLICIT NONE

    CALL scalc_integrate

  END SUBROUTINE scalc_global_end
  ! ---------------------------------------------------------------------

  ! ---------------------------------------------------------------------
  SUBROUTINE scalc_integrate(jrow,kproma)

    ! SCALC MODULE ROUTINE
    !
    ! calculations
    !
    ! Author: Bastian Kern, MPIC, 08.07.2011
    !
    USE messy_main_channel_bi,    ONLY: GP_3D_1LEV
    USE messy_main_tools,         ONLY: spechum2mr, mr2spechum
    USE messy_main_blather_bi,    ONLY: error_bi
    USE messy_main_control,       ONLY: main_control_is_context
    USE messy_main_channel_mem,   ONLY: dom_curid

    IMPLICIT NONE
    INTRINSIC :: SUM

    ! I/O
    INTEGER, OPTIONAL    :: jrow
    INTEGER, OPTIONAL    :: kproma

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER  :: substr = 'scalc_integrate'
    INTEGER :: i, jo, ix, ilev
    INTEGER :: calc_reprid
    INTEGER :: i1, i2, j1, j2      ! local object-wise used indices
    INTEGER :: ii1, ii2, jj1, jj2  ! locally set indices
    INTEGER, PARAMETER :: I_UNSET = -99

    LOGICAL :: is_context   ! check call string

    REAL(DP), DIMENSION(:,:), POINTER :: sum_ptr => NULL()
    REAL(DP), DIMENSION(:,:), POINTER :: dia_ptr => NULL()
    REAL(DP), PARAMETER               :: tiny =  1e-20_dp
    REAL(DP), PARAMETER               :: epsilon =  1e-24_dp

    ! INITIALIZE
    IF (PRESENT(jrow)) THEN
       jj1 = 1
       jj2 = jrow
    ELSE
       jj1 = I_UNSET
       jj2 = I_UNSET
    END IF
    IF (PRESENT(kproma)) THEN
       ii1 = 1
       ii2 = kproma
    ELSE
       ii1 = I_UNSET
       ii2 = I_UNSET
    END IF

    calc_loop: DO i=1, NCALC

       IF (XCALC(i)%io%domain_idx /= dom_curid) CYCLE

       is_context = .FALSE.
       DO ix = 1,XCALC(i)%io%nwhen
          is_context = is_context .OR.  &
               main_control_is_context(XCALC(i)%io%when(ix) &
               , XCALC(i)%io%where(ix))
       END DO

       IF (.NOT. is_context) CYCLE

       is_ok: IF (XCALC(i)%ok) THEN

          ! DETERMINE HORIZITONAL SPACE INDICES
          IF (ii2 == I_UNSET) THEN
             i1=XCALC(i)%xind(1)
             i2=XCALC(i)%xind(2)
          ELSE
             i1 = ii1
             i2 = ii2
          END IF

          IF (jj2 == I_UNSET) THEN
             j1=XCALC(i)%yind(1)
             j2=XCALC(i)%yind(2)
          ELSE
             j1 = jj1
             j2 = jj2
          END IF

          SELECT CASE (XCALC(i)%math_func)

          CASE(OSUM, PSUM, PSUMDIA, QSUM) ! -----------------------------------
             ! sum (with scaling factors), positive confined sum,
             ! positive confined sum with error diagnostics
             calc_reprid = XCALC(i)%rid(1)
             IF (XCALC(i)%spgeo == IS_2D) THEN
                IF (XCALC(i)%math_func == PSUMDIA) &
                     & XCALC(i)%diag%ptr(i1:i2,j1:j2,:) = 0.0_dp

                XCALC(i)%result2d(1)%ptr(i1:i2,j1:j2) = 0.0_dp
                sum_ptr => XCALC(i)%result2d(1)%ptr

                IF (XCALC(i)%math_func == PSUMDIA) THEN
                   dia_ptr => XCALC(i)%diag%ptr(:,:,1)
                   IF (calc_reprid == GP_3D_1LEV) THEN
                      dia_ptr => XCALC(i)%diag%ptr(_RI_XYZ__(:,:,1))
                   END IF
                   dia_ptr = 0.0_dp
                END IF ! PSUMDIA

                object_loop_2d_1: DO jo=1, XCALC(i)%nobj
                   IF (XCALC(i)%lex(jo)) THEN
                      ! For the summation of water tracers, we have to convert
                      ! to mixing ratios, if necessary.
                      ! All calculations are in mixing ratios.
                      IF ((XCALC(i)%math_func == QSUM) .AND.   &
                           ((XCALC(i)%unit(jo) == "kg/kg") .OR. &
                           (XCALC(i)%unit(jo) == "kg kg-1"))   &
                           ) THEN
                         sum_ptr(i1:i2,j1:j2) = sum_ptr(i1:i2,j1:j2)  &
                            + spechum2mr(XCALC(i)%dat2d(jo)%ptr(i1:i2,j1:j2)) &
                            * XCALC(i)%fac(jo)
                      ELSE
                         sum_ptr(i1:i2,j1:j2) = sum_ptr(i1:i2,j1:j2)  &
                              + XCALC(i)%dat2d(jo)%ptr(i1:i2,j1:j2)   &
                              * XCALC(i)%fac(jo)
                      END IF
                   END IF
                END DO object_loop_2d_1 ! XCALC(i)%nobj

                IF (XCALC(i)%math_func == PSUM) THEN
                   WHERE (sum_ptr(i1:i2,j1:j2) .le. 0._dp) &
                        sum_ptr(i1:i2,j1:j2) = 0._dp
                END IF
                IF (XCALC(i)%math_func == PSUMDIA) THEN
                   WHERE (sum_ptr(i1:i2,j1:j2) .le. 0._dp)
                      dia_ptr(i1:i2,j1:j2) = sum_ptr(i1:i2,j1:j2)
                      sum_ptr(i1:i2,j1:j2) = 0._dp
                   END WHERE
                END IF

                ! In case of QSUM, the first object determines output units.
                ! If it is kg/kg, the variable has to be converted to
                ! spefic humidity.
                IF ((XCALC(i)%math_func == QSUM) .AND.  &
                     ((XCALC(i)%unit(1) == "kg/kg") .OR. &
                     (XCALC(i)%unit(1) == "kg kg-1")))  &
                     XCALC(i)%result2d(1)%ptr(i1:i2,j1:j2) = &
                     mr2spechum(XCALC(i)%result2d(1)%ptr(i1:i2,j1:j2))

             ELSE ! END 2D case

                IF (XCALC(i)%spgeo == IS_3D_N) THEN

                   XCALC(i)%result3d(1)%ptr(_RI_XY_N_(i1:i2,j1:j2,:)) = 0.0_dp
                   object_loop_1: DO jo=1, XCALC(i)%nobj
                      IF (XCALC(i)%lex(jo)) THEN
                         ! For the summation of water tracers,
                         ! we have to convert
                         ! to mixing ratios, if necessary.
                         ! All calculations are in mixing ratios.
                         IF ((XCALC(i)%math_func == QSUM) .AND.   &
                              ((XCALC(i)%unit(jo) == "kg/kg") .OR. &
                              (XCALC(i)%unit(jo) == "kg kg-1"))   &
                              ) THEN
                            XCALC(i)%result3d(1)%ptr(_RI_XY_N_(i1:i2,j1:j2,:)) = &
                                 XCALC(i)%result3d(1)%ptr(_RI_XY_N_(i1:i2,j1:j2,:)) &
                                 + spechum2mr(XCALC(i)%dat3d(jo)%ptr(_RI_XY_N_(i1:i2,j1:j2,:))) &
                                 * XCALC(i)%fac(jo)
                         ELSE
                            XCALC(i)%result3d(1)%ptr(_RI_XY_N_(i1:i2,j1:j2,:)) =&
                                 XCALC(i)%result3d(1)%ptr(_RI_XY_N_(i1:i2,j1:j2,:))  &
                                 + XCALC(i)%dat3d(jo)%ptr(_RI_XY_N_(i1:i2,j1:j2,:)) * &
                                 XCALC(i)%fac(jo)
                         END IF
                      END IF

                   END DO object_loop_1 ! XCALC(i)%nobj

                   IF (XCALC(i)%math_func == PSUM) THEN
                      WHERE (XCALC(i)%result3d(1)%ptr(_RI_XY_N_(i1:i2,j1:j2,:)) .le. 0._dp) &
                           & XCALC(i)%result3d(1)%ptr(_RI_XY_N_(i1:i2,j1:j2,:)) = 0._dp
                   END IF
                   IF (XCALC(i)%math_func == PSUMDIA) THEN
                      XCALC(i)%diag%ptr(_RI_XY_N_(i1:i2,j1:j2,:)) = 0.0_dp
                      WHERE (XCALC(i)%result3d(1)%ptr(_RI_XY_N_(i1:i2,j1:j2,:)) .le. 0._dp)
                         XCALC(i)%diag%ptr(_RI_XY_N_(i1:i2,j1:j2,:)) = &
                              XCALC(i)%result3d(1)%ptr(_RI_XY_N_(i1:i2,j1:j2,:))
                         XCALC(i)%result3d(1)%ptr(_RI_XY_N_(i1:i2,j1:j2,:)) = 0._dp
                      END WHERE
                   END IF
                   ! In case of QSUM, the first object determines output units.
                   ! If it is kg/kg, the variable has to be converted
                   ! to spefic humidity.
                   IF ((XCALC(i)%math_func == QSUM) .AND.   &
                        ((XCALC(i)%unit(1) == "kg/kg") .OR. &
                        (XCALC(i)%unit(1) == "kg kg-1")))   &
                        XCALC(i)%result3d(1)%ptr(_RI_XY_N_(i1:i2,j1:j2,:)) = &
                        mr2spechum(XCALC(i)%result3d(1)%ptr(_RI_XY_N_(i1:i2,j1:j2,:)))

                ELSE IF (XCALC(i)%spgeo == IS_3D_Z) THEN
                   XCALC(i)%result3d(1)%ptr(_RI_XYZ__(i1:i2,j1:j2,:)) = 0.0_dp
                   object_loop_1_1: DO jo=1, XCALC(i)%nobj
                      IF (XCALC(i)%lex(jo)) THEN
                         ! For the summation of water tracers,
                         ! we have to convert
                         ! to mixing ratios, if necessary.
                         ! All calculations are in mixing ratios.
                         IF ((XCALC(i)%math_func == QSUM) .AND.   &
                              ((XCALC(i)%unit(jo) == "kg/kg") .OR. &
                              (XCALC(i)%unit(jo) == "kg kg-1"))   &
                              ) THEN
                            XCALC(i)%result3d(1)%ptr(_RI_XYZ__(i1:i2,j1:j2,:)) =&
                                 XCALC(i)%result3d(1)%ptr(_RI_XYZ__(i1:i2,j1:j2,:)) &
                                 + spechum2mr(XCALC(i)%dat3d(jo)%ptr(_RI_XYZ__(i1:i2,j1:j2,:))) &
                                 * XCALC(i)%fac(jo)
                         ELSE
                            XCALC(i)%result3d(1)%ptr(_RI_XYZ__(i1:i2,j1:j2,:)) =&
                                 XCALC(i)%result3d(1)%ptr(_RI_XYZ__(i1:i2,j1:j2,:))  &
                                 + XCALC(i)%dat3d(jo)%ptr(_RI_XYZ__(i1:i2,j1:j2,:)) * &
                                 XCALC(i)%fac(jo)
                         END IF
                      END IF

                   END DO object_loop_1_1 ! XCALC(i)%nobj

                   IF (XCALC(i)%math_func == PSUM) THEN
                      WHERE (XCALC(i)%result3d(1)%ptr(_RI_XYZ__(i1:i2,j1:j2,:)) .le. 0._dp) &
                           & XCALC(i)%result3d(1)%ptr(_RI_XYZ__(i1:i2,j1:j2,:)) = 0._dp
                   END IF
                   IF (XCALC(i)%math_func == PSUMDIA) THEN
                      XCALC(i)%diag%ptr(_RI_XYZ__(i1:i2,j1:j2,:)) = 0.0_dp
                      WHERE (XCALC(i)%result3d(1)%ptr(_RI_XYZ__(i1:i2,j1:j2,:)) .le. 0._dp)
                         XCALC(i)%diag%ptr(_RI_XYZ__(i1:i2,j1:j2,:)) = &
                              XCALC(i)%result3d(1)%ptr(_RI_XYZ__(i1:i2,j1:j2,:))
                         XCALC(i)%result3d(1)%ptr(_RI_XYZ__(i1:i2,j1:j2,:)) = 0._dp
                      END WHERE
                   END IF
                   ! In case of QSUM, the first object determines output units.
                   ! If it is kg/kg, the variable has to be converted
                   ! to spefic humidity.
                   IF ((XCALC(i)%math_func == QSUM) .AND.  &
                        ((XCALC(i)%unit(1) == "kg/kg") .OR. &
                        (XCALC(i)%unit(1) == "kg kg-1")))  &
                        XCALC(i)%result3d(1)%ptr(_RI_XYZ__(i1:i2,j1:j2,:)) = &
                        mr2spechum(XCALC(i)%result3d(1)%ptr(_RI_XYZ__(i1:i2,j1:j2,:)))

                ELSE ! third axis string is neither N nor Z
                   !should never occur
                   CALL error_bi('Third axis string is undefined. &
                        &(This Error should NEVER occur!)', substr)

                END IF !N or Z

             END IF ! rank == 2

          CASE(SUMGE0) ! ------------------------------------------------------
             IF (XCALC(i)%nobj /= 2) THEN
                ! This should never be reached!
                CALL error_bi('For SUMGE0 only 2 objects are allowed. &
                     &(This Error should NEVER occur!)', substr)
             END IF
             ! SUM 2 fields, if field 1 is .le. a tiny number (1E-20)
             IF (XCALC(i)%spgeo == IS_2D) THEN
                XCALC(i)%result2d(1)%ptr(i1:i2,j1:j2) = 0.0_dp
                sum_ptr => XCALC(i)%result2d(1)%ptr
                object_loop_2d_2: DO jo=1, XCALC(i)%nobj
                   IF (XCALC(i)%lex(jo)) THEN
                      where (sum_ptr(i1:i2,j1:j2) .le. tiny) &
                           sum_ptr(i1:i2,j1:j2) = &
                           sum_ptr(i1:i2,j1:j2) + &
                           XCALC(i)%dat2d(jo)%ptr(i1:i2,j1:j2) * &
                           XCALC(i)%fac(jo)
                   END IF
                END DO object_loop_2d_2 ! XCALC(i)%nobj
             ELSE !2D case
                IF (XCALC(i)%spgeo == IS_3D_N) THEN

                   XCALC(i)%result3d(1)%ptr(_RI_XY_N_(i1:i2,j1:j2,:)) = 0.0_dp
                   object_loop_2: DO jo=1, XCALC(i)%nobj
                      IF (XCALC(i)%lex(jo)) THEN
                         where (XCALC(i)%result3d(1)%ptr(_RI_XY_N_(i1:i2,j1:j2,:)) .le. tiny)               &
                              XCALC(i)%result3d(1)%ptr(_RI_XY_N_(i1:i2,j1:j2,:)) = &
                              XCALC(i)%result3d(1)%ptr(_RI_XY_N_(i1:i2,j1:j2,:)) + &
                              XCALC(i)%dat3d(jo)%ptr(_RI_XY_N_(i1:i2,j1:j2,:)) *&
                              XCALC(i)%fac(jo)
                      END IF
                   END DO object_loop_2 ! XCALC(i)%nobj

                ELSE IF (XCALC(i)%spgeo == IS_3D_Z) THEN

                   XCALC(i)%result3d(1)%ptr(_RI_XYZ__(i1:i2,j1:j2,:)) = 0.0_dp
                   object_loop_2_1: DO jo=1, XCALC(i)%nobj
                      IF (XCALC(i)%lex(jo)) THEN
                         where (XCALC(i)%result3d(1)%ptr(_RI_XYZ__(i1:i2,j1:j2,:)) .le. tiny)               &
                              XCALC(i)%result3d(1)%ptr(_RI_XYZ__(i1:i2,j1:j2,:)) = &
                              XCALC(i)%result3d(1)%ptr(_RI_XYZ__(i1:i2,j1:j2,:)) + &
                              XCALC(i)%dat3d(jo)%ptr(_RI_XYZ__(i1:i2,j1:j2,:)) *&
                              XCALC(i)%fac(jo)
                      END IF
                   END DO object_loop_2_1 ! XCALC(i)%nobj

                ELSE ! third axis string is neither N nor Z
                   !should never occur

                   CALL error_bi('Third axis string is undefined. &
                        &(This Error should NEVER occur!)', substr)

                END IF !N or Z

             END IF ! mode_2d

          CASE(MULT) !---------------------------------------------------------
             ! MULTIPLY fields with scaling
             IF (XCALC(i)%spgeo == IS_2D  ) THEN
                XCALC(i)%result2d(1)%ptr(i1:i2,j1:j2) = 1.0_dp
                sum_ptr => XCALC(i)%result2d(1)%ptr
                object_loop_2d_3: DO jo=1, XCALC(i)%nobj
                   IF (XCALC(i)%lex(jo)) THEN
                      sum_ptr(i1:i2,j1:j2) = sum_ptr(i1:i2,j1:j2) * &
                           XCALC(i)%dat2d(jo)%ptr(i1:i2,j1:j2) *&
                           XCALC(i)%fac(jo)
                   END IF
                END DO object_loop_2d_3 ! XCALC(i)%nobj
             ELSE !2D
                IF (XCALC(i)%spgeo == IS_3D_N) THEN
                   XCALC(i)%result3d(1)%ptr(_RI_XY_N_(i1:i2,j1:j2,:)) = 1.0_dp
                   object_loop_3: DO jo=1, XCALC(i)%nobj
                      IF (XCALC(i)%lex(jo)) THEN
                         XCALC(i)%result3d(1)%ptr(_RI_XY_N_(i1:i2,j1:j2,:)) = &
                              XCALC(i)%result3d(1)%ptr(_RI_XY_N_(i1:i2,j1:j2,:)) * &
                              XCALC(i)%dat3d(jo)%ptr(_RI_XY_N_(i1:i2,j1:j2,:))* &
                              XCALC(i)%fac(jo)
                      END IF
                   END DO object_loop_3 ! XCALC(i)%nobj

                ELSE IF (XCALC(i)%spgeo == IS_3D_Z) THEN

                   XCALC(i)%result3d(1)%ptr(_RI_XYZ__(i1:i2,j1:j2,:)) = 1.0_dp
                   object_loop_3_1: DO jo=1, XCALC(i)%nobj
                      IF (XCALC(i)%lex(jo)) THEN
                         XCALC(i)%result3d(1)%ptr(_RI_XYZ__(i1:i2,j1:j2,:)) = &
                              XCALC(i)%result3d(1)%ptr(_RI_XYZ__(i1:i2,j1:j2,:)) * &
                              XCALC(i)%dat3d(jo)%ptr(_RI_XYZ__(i1:i2,j1:j2,:))* &
                              XCALC(i)%fac(jo)
                      END IF
                   END DO object_loop_3_1 ! XCALC(i)%nobj

                ELSE ! third axis string is neither N nor Z
                   !should never occur

                   CALL error_bi('Third axis string is undefined. &
                        &(This Error should NEVER occur!)', substr)

                END IF !N or Z

             END IF ! mode_2d

          CASE(DIV) !----------------------------------------------------------
             ! DIVIDE fields with scaling
             IF (XCALC(i)%spgeo == IS_2D  ) THEN
                XCALC(i)%result2d(1)%ptr(i1:i2,j1:j2) = 0.0_dp
                sum_ptr => XCALC(i)%result2d(1)%ptr
                sum_ptr(i1:i2,j1:j2) = XCALC(i)%dat2d(1)%ptr(i1:i2,j1:j2) *&
                     XCALC(i)%fac(1)
                object_loop_2d_4: DO jo=2, XCALC(i)%nobj
                   IF (XCALC(i)%lex(jo)) THEN
                      WHERE ( ABS( XCALC(i)%dat2d(jo)%ptr(i1:i2,j1:j2) *&
                           XCALC(i)%fac(jo)) &
                           .GE. epsilon)
                         sum_ptr(i1:i2,j1:j2) = sum_ptr(i1:i2,j1:j2) / &
                              (XCALC(i)%dat2d(jo)%ptr(i1:i2,j1:j2) *   &
                              XCALC(i)%fac(jo))
                      ELSEWHERE
                         sum_ptr(i1:i2,j1:j2) = 0.0_dp
                      END WHERE
                   END IF
                END DO object_loop_2d_4 ! XCALC(i)%nobj
             ELSE !2D
                IF (XCALC(i)%spgeo == IS_3D_N) THEN
                   XCALC(i)%result3d(1)%ptr(_RI_XY_N_(i1:i2,j1:j2,:)) = 0.0_dp
                   XCALC(i)%result3d(1)%ptr(_RI_XY_N_(i1:i2,j1:j2,:)) = &
                        XCALC(i)%dat3d(1)%ptr(_RI_XY_N_(i1:i2,j1:j2,:)) * &
                        XCALC(i)%fac(1)
                   object_loop_4: DO jo=2, XCALC(i)%nobj
                      ! if the file exists ... perform calculation
                      IF (XCALC(i)%lex(jo)) THEN
                         WHERE ( ABS(XCALC(i)%dat3d(jo)%ptr(_RI_XY_N_(i1:i2,j1:j2,:)) * &
                              XCALC(i)%fac(jo)) .GE. epsilon)
                            XCALC(i)%result3d(1)%ptr(_RI_XY_N_(i1:i2,j1:j2,:)) = &
                                 XCALC(i)%result3d(1)%ptr(_RI_XY_N_(i1:i2,j1:j2,:)) /   &
                                 (XCALC(i)%dat3d(jo)%ptr(_RI_XY_N_(i1:i2,j1:j2,:)) *      &
                                 XCALC(i)%fac(jo))
                         ELSEWHERE
                            XCALC(i)%result3d(1)%ptr(_RI_XY_N_(i1:i2,j1:j2,:)) =&
                                 0.0_dp
                         END WHERE
                      END IF
                   END DO object_loop_4 ! XCALC(i)%nobj

                ELSE IF (XCALC(i)%spgeo == IS_3D_Z) THEN
                   XCALC(i)%result3d(1)%ptr(_RI_XYZ__(i1:i2,j1:j2,:)) = 0.0_dp
                   XCALC(i)%result3d(1)%ptr(_RI_XYZ__(i1:i2,j1:j2,:)) = &
                        XCALC(i)%dat3d(1)%ptr(_RI_XYZ__(i1:i2,j1:j2,:)) * &
                        XCALC(i)%fac(1)
                   object_loop_4_1: DO jo=2, XCALC(i)%nobj
                      ! if the file exists ... perform calculation
                      IF (XCALC(i)%lex(jo)) THEN
                         WHERE ( ABS(XCALC(i)%dat3d(jo)%ptr(_RI_XYZ__(i1:i2,j1:j2,:)) * &
                              XCALC(i)%fac(jo)) .GE. epsilon)
                            XCALC(i)%result3d(1)%ptr(_RI_XYZ__(i1:i2,j1:j2,:)) = &
                                 XCALC(i)%result3d(1)%ptr(_RI_XYZ__(i1:i2,j1:j2,:)) /   &
                                 (XCALC(i)%dat3d(jo)%ptr(_RI_XYZ__(i1:i2,j1:j2,:)) *      &
                                 XCALC(i)%fac(jo))
                         ELSEWHERE
                            XCALC(i)%result3d(1)%ptr(_RI_XYZ__(i1:i2,j1:j2,:)) = 0.0_dp
                         END WHERE
                      END IF
                   END DO object_loop_4_1 ! XCALC(i)%nobj

                ELSE ! third axis string is neither N nor Z
                   !should never occur

                   CALL error_bi('Third axis string is undefined. &
                        &(This Error should NEVER occur!)', substr)

                END IF !N or Z

             END IF ! mode_2d

          CASE(FLAG) !-------------------------------------------------------

             IF (XCALC(i)%spgeo == IS_2D ) THEN
                XCALC(i)%result2d(1)%ptr(i1:i2,j1:j2) = 1.0_dp
                sum_ptr => XCALC(i)%result2d(1)%ptr
                sum_ptr(i1:i2,j1:j2) = XCALC(i)%dat2d(1)%ptr(i1:i2,j1:j2)*&
                     XCALC(i)%dat2d(2)%ptr(i1:i2,j1:j2)
             ELSE !2D
                IF (XCALC(i)%spgeo == IS_3D_N) THEN
                   XCALC(i)%result3d(1)%ptr(_RI_XY_N_(i1:i2,j1:j2,:)) = 1.0_dp
                   DO ilev = 1,size(XCALC(i)%result3d(1)%ptr,(_IN_XY_N_))
                      XCALC(i)%result3d(1)%ptr(_RI_XY_N_(i1:i2,j1:j2,ilev))= &
                           XCALC(i)%dat3d(1)%ptr(_RI_XY_N_(i1:i2,j1:j2,ilev)) * &
                           XCALC(i)%dat2d(2)%ptr(i1:i2,j1:j2)
                   END DO
                ELSE IF (XCALC(i)%spgeo == IS_3D_Z) THEN
                   XCALC(i)%result3d(1)%ptr(_RI_XYZ__(i1:i2,j1:j2,:)) = 1.0_dp
                   DO ilev = 1,size(XCALC(i)%result3d(1)%ptr,(_IZ_XYZ__))
                      XCALC(i)%result3d(1)%ptr(_RI_XYZ__(i1:i2,j1:j2,ilev))= &
                           XCALC(i)%dat3d(1)%ptr(_RI_XYZ__(i1:i2,j1:j2,ilev)) * &
                           XCALC(i)%dat2d(2)%ptr(i1:i2,j1:j2)
                   END DO

                ELSE ! third axis string is neither N nor Z
                   !should never occur

                   CALL error_bi('Third axis string is undefined. &
                        &(This Error should NEVER occur!)', substr)

                END IF !N or Z

             END IF ! mode_2d

          CASE(FLAGU) !------------------------------------------------------
             IF (XCALC(i)%spgeo == IS_2D ) THEN
                sum_ptr => XCALC(i)%result2d(1)%ptr
                ! Copy flag field on sum ptr to modify the sum ptr
                ! and not the original IMPORT field.
                sum_ptr(i1:i2,j1:j2) = XCALC(i)%dat2d(2)%ptr(i1:i2,j1:j2)
                where(sum_ptr(i1:i2,j1:j2) >= 0.01_dp ) sum_ptr(i1:i2,j1:j2)=1.0_dp
                sum_ptr(i1:i2,j1:j2) = XCALC(i)%dat2d(1)%ptr(i1:i2,j1:j2)* &
                     sum_ptr(i1:i2,j1:j2)
             ELSE !2D
                IF (XCALC(i)%spgeo == IS_3D_N) THEN
                   DO ilev = 1,size(XCALC(i)%result3d(1)%ptr,(_IN_XY_N_))
                      ! Copy flag field on sum ptr to modify the sum ptr
                      ! and not the original IMPORT field .
                      XCALC(i)%result3d(1)%ptr(_RI_XY_N_(i1:i2,j1:j2,ilev))=  &
                           XCALC(i)%dat2d(2)%ptr(i1:i2,j1:j2)
                      where( XCALC(i)%result3d(1)%ptr(_RI_XY_N_(i1:i2,j1:j2,ilev))&
                           >= 0.01_dp ) &
                           XCALC(i)%result3d(1)%ptr(_RI_XY_N_(i1:i2,j1:j2,ilev))  =1.0_dp
                      XCALC(i)%result3d(1)%ptr(_RI_XY_N_(i1:i2,j1:j2,ilev))= &
                           XCALC(i)%dat3d(1)%ptr(_RI_XY_N_(i1:i2,j1:j2,ilev)) * &
                           XCALC(i)%result3d(1)%ptr(_RI_XY_N_(i1:i2,j1:j2,ilev))
                   END DO

                ELSE IF (XCALC(i)%spgeo == IS_3D_Z) THEN
                   DO ilev = 1,size(XCALC(i)%result3d(1)%ptr,(_IZ_XYZ__))
                      ! Copy flag field on sum ptr to modify the sum ptr
                      ! and not the original IMPORT field .
                      XCALC(i)%result3d(1)%ptr(_RI_XYZ__(i1:i2,j1:j2,ilev))=  &
                           XCALC(i)%dat2d(2)%ptr(i1:i2,j1:j2)
                      where( XCALC(i)%result3d(1)%ptr(_RI_XYZ__(i1:i2,j1:j2,ilev))&
                           >= 0.01_dp ) &
                           XCALC(i)%result3d(1)%ptr(_RI_XYZ__(i1:i2,j1:j2,ilev))  = 1.0_dp
                      XCALC(i)%result3d(1)%ptr(_RI_XYZ__(i1:i2,j1:j2,ilev))= &
                           XCALC(i)%dat3d(1)%ptr(_RI_XYZ__(i1:i2,j1:j2,ilev)) * &
                           XCALC(i)%result3d(1)%ptr(_RI_XYZ__(i1:i2,j1:j2,ilev))
                   END DO

                ELSE ! third axis string is neither N nor Z
                   !should never occur

                   CALL error_bi('Third axis string is undefined. &
                        &(This Error should NEVER occur!)', substr)

                END IF !N or Z

             END IF ! mode_2d

          CASE(FLAGD) !------------------------------------------------------
             IF (XCALC(i)%spgeo == IS_2D ) THEN
                sum_ptr => XCALC(i)%result2d(1)%ptr
                ! Copy flag field on sum ptr to modify the sum ptr
                ! and not the original IMPORT field.
                sum_ptr(i1:i2,j1:j2) = XCALC(i)%dat2d(2)%ptr(i1:i2,j1:j2)
                where(sum_ptr(i1:i2,j1:j2) <= 0.99_dp ) sum_ptr(i1:i2,j1:j2)=0.0_dp
                sum_ptr(i1:i2,j1:j2) = XCALC(i)%dat2d(1)%ptr(i1:i2,j1:j2) *&
                     sum_ptr (i1:i2,j1:j2)
             ELSE
                IF (XCALC(i)%spgeo == IS_3D_N) THEN

                   DO ilev = 1,size(XCALC(i)%result3d(1)%ptr,(_IN_XY_N_))
                      ! Copy flag field on sum ptr to modify the sum ptr
                      ! and not the original IMPORT field.
                      XCALC(i)%result3d(1)%ptr(_RI_XY_N_(i1:i2,j1:j2,ilev))=  &
                           XCALC(i)%dat2d(2)%ptr(i1:i2,j1:j2)

                      where( XCALC(i)%result3d(1)%ptr(_RI_XY_N_(i1:i2,j1:j2,ilev)) &
                           <= 0.99_dp ) &
                           XCALC(i)%result3d(1)%ptr(_RI_XY_N_(i1:i2,j1:j2,ilev)) = 0.0_dp

                      XCALC(i)%result3d(1)%ptr(_RI_XY_N_(i1:i2,j1:j2,ilev))= &
                           XCALC(i)%dat3d(1)%ptr(_RI_XY_N_(i1:i2,j1:j2,ilev)) * &
                           XCALC(i)%result3d(1)%ptr(_RI_XY_N_(i1:i2,j1:j2,ilev))
                   END DO


                ELSE IF (XCALC(i)%spgeo == IS_3D_Z) THEN

                   DO ilev = 1,size(XCALC(i)%result3d(1)%ptr,(_IZ_XYZ__))
                      ! Copy flag field on sum ptr to modify the sum ptr
                      ! and not the original IMPORT field.
                      XCALC(i)%result3d(1)%ptr(_RI_XYZ__(i1:i2,j1:j2,ilev))=  &
                           XCALC(i)%dat2d(2)%ptr(i1:i2,j1:j2)

                      where( XCALC(i)%result3d(1)%ptr(_RI_XYZ__(i1:i2,j1:j2,ilev)) &
                           <= 0.99_dp ) &
                           XCALC(i)%result3d(1)%ptr(_RI_XYZ__(i1:i2,j1:j2,ilev)) = 0.0_dp

                      XCALC(i)%result3d(1)%ptr(_RI_XYZ__(i1:i2,j1:j2,ilev))= &
                           XCALC(i)%dat3d(1)%ptr(_RI_XYZ__(i1:i2,j1:j2,ilev)) * &
                           XCALC(i)%result3d(1)%ptr(_RI_XYZ__(i1:i2,j1:j2,ilev))
                   END DO


                ELSE ! third axis string is neither N nor Z
                   !should never occur

                   CALL error_bi('Third axis string is undefined. &
                        &(This Error should NEVER occur!)', substr)

                END IF !N or Z


             END IF ! mode_2d

          CASE(DELTA) !------------------------------------------------------

             delta_2d: IF (XCALC(i)%spgeo == IS_2D ) THEN
                XCALC(i)%result2d(1)%ptr(i1:i2,j1:j2) = 0.0_dp
                XCALC(i)%result2d(2)%ptr(i1:i2,j1:j2) = 0.0_dp
                !
                sum_ptr => XCALC(i)%result2d(1)%ptr(i1:i2,j1:j2) ! OUT
                dia_ptr => XCALC(i)%dat2d(1)%ptr(i1:i2,j1:j2)   ! IN
                CALL delta2frac( &
                     sum_ptr(:,:) &
                     , dia_ptr(:,:) &
                     , XCALC(i)%xparam(1), XCALC(i)%xparam(2) )
                XCALC(i)%result2d(2)%ptr(i1:i2,j1:j2) = &
                     1.0_dp - sum_ptr(:,:)
             ELSE !3D
                delta_3d_n: IF (XCALC(i)%spgeo == IS_3D_N) THEN
                   XCALC(i)%result3d(1)%ptr(_RI_XY_N_(i1:i2,j1:j2,:)) = 0.0_dp
                   XCALC(i)%result3d(2)%ptr(_RI_XY_N_(i1:i2,j1:j2,:)) = 0.0_dp
                   DO ilev = 1,size(XCALC(i)%result3d(1)%ptr,(_IN_XY_N_))
                      sum_ptr => &
                           XCALC(i)%result3d(1)%ptr(_RI_XY_N_(i1:i2,j1:j2,ilev))
                      dia_ptr => &
                           XCALC(i)%dat3d(1)%ptr(_RI_XY_N_(i1:i2,j1:j2,ilev))
                      CALL delta2frac( &
                           sum_ptr(:,:) &
                           , dia_ptr(:,:) &
                           , XCALC(i)%xparam(1), XCALC(i)%xparam(2) )
                      XCALC(i)%result3d(2)%ptr(_RI_XY_N_(i1:i2,j1:j2,ilev)) = &
                           1.0_dp - sum_ptr(:,:)
                   END DO
                ELSE IF (XCALC(i)%spgeo == IS_3D_Z) THEN
                   XCALC(i)%result3d(1)%ptr(_RI_XYZ__(i1:i2,j1:j2,:)) = 0.0_dp
                   XCALC(i)%result3d(1)%ptr(_RI_XYZ__(i1:i2,j1:j2,:)) = 0.0_dp
                   DO ilev = 1,size(XCALC(i)%result3d(1)%ptr,(_IZ_XYZ__))
                      sum_ptr => &
                           XCALC(i)%result3d(1)%ptr(_RI_XYZ__(i1:i2,j1:j2,ilev))
                      dia_ptr => &
                           XCALC(i)%dat3d(1)%ptr(_RI_XYZ__(i1:i2,j1:j2,ilev))
                      CALL delta2frac( &
                           sum_ptr(:,:) &
                           , dia_ptr(:,:) &
                           , XCALC(i)%xparam(1), XCALC(i)%xparam(2) )
                      XCALC(i)%result3d(2)%ptr(_RI_XYZ__(i1:i2,j1:j2,ilev)) = &
                           1.0_dp - sum_ptr(:,:)
                   END DO
                ELSE ! third axis string is neither N nor Z
                   !should never occur
                   CALL error_bi('Third axis string is undefined. &
                        &(This Error should NEVER occur!)', substr)
                END IF delta_3d_n
             END IF delta_2d

          CASE(WVA) !---------------------------------------------------------
             ! weighted vertical average, only defined for 3D, result is 2D
             !
             XCALC(i)%result2d(1)%ptr(i1:i2,j1:j2) = 0.0
             !
             IF (XCALC(i)%spgeo == IS_3D_N) THEN
                XCALC(i)%result2d(1)%ptr(i1:i2,j1:j2) = SUM( &
                       ( XCALC(i)%dat3d(1)%ptr(_RI_XY_N_(i1:i2,j1:j2,:)) * &
                         XCALC(i)%dat3d(2)%ptr(_RI_XY_N_(i1:i2,j1:j2,:)) ) &
                       , DIM=_IN_XY_N_) &
                      / SUM( &
                         XCALC(i)%dat3d(2)%ptr(_RI_XY_N_(i1:i2,j1:j2,:)) &
                       , DIM=_IN_XY_N_)
             ELSE IF (XCALC(i)%spgeo == IS_3D_Z) THEN
                XCALC(i)%result2d(1)%ptr(i1:i2,j1:j2) = SUM( &
                       ( XCALC(i)%dat3d(1)%ptr(_RI_XYZ__(i1:i2,j1:j2,:)) * &
                         XCALC(i)%dat3d(2)%ptr(_RI_XYZ__(i1:i2,j1:j2,:)) ) &
                       , DIM=_IZ_XYZ__) &
                      / SUM( &
                         XCALC(i)%dat3d(2)%ptr(_RI_XYZ__(i1:i2,j1:j2,:)) &
                       , DIM=_IZ_XYZ__)
             ELSE ! third axis string is neither N nor Z
                !should never occur
                CALL error_bi('Third axis string is undefined. &
                     &(This Error should NEVER occur!)', substr)
             END IF !N or Z

          CASE(SPLITV) !------------------------------------------------------

             ! use 2d field to split into result3d(1) and result3d(2)
             XCALC(i)%result3d(1)%ptr(_RI_XYZ__(i1:i2,j1:j2,:)) = 1.0_dp
             XCALC(i)%result3d(2)%ptr(_RI_XYZ__(i1:i2,j1:j2,:)) = 1.0_dp

             DO ilev = 1,size(XCALC(i)%result3d(1)%ptr,(_IZ_XYZ__))
                IF (XCALC(i)%xparam(1) == 0) THEN
                   WHERE ( XCALC(i)%dat2d(2)%ptr(i1:i2,j1:j2) <= REAL(ilev,dp) )
                      XCALC(i)%result3d(1)%ptr(_RI_XYZ__(i1:i2,j1:j2,ilev))= &
                           XCALC(i)%dat3d(1)%ptr(_RI_XYZ__(i1:i2,j1:j2,ilev))
                      XCALC(i)%result3d(2)%ptr(_RI_XYZ__(i1:i2,j1:j2,ilev))= &
                           0.0_dp
                   ELSEWHERE
                      XCALC(i)%result3d(1)%ptr(_RI_XYZ__(i1:i2,j1:j2,ilev))= &
                           0.0_dp
                      XCALC(i)%result3d(2)%ptr(_RI_XYZ__(i1:i2,j1:j2,ilev))= &
                           XCALC(i)%dat3d(1)%ptr(_RI_XYZ__(i1:i2,j1:j2,ilev))
                   END WHERE
                ELSE IF (XCALC(i)%xparam(1) == 1) THEN
                   WHERE ( XCALC(i)%dat2d(2)%ptr(i1:i2,j1:j2) < REAL(ilev,dp) )
                      XCALC(i)%result3d(1)%ptr(_RI_XYZ__(i1:i2,j1:j2,ilev))= &
                           XCALC(i)%dat3d(1)%ptr(_RI_XYZ__(i1:i2,j1:j2,ilev))
                      XCALC(i)%result3d(2)%ptr(_RI_XYZ__(i1:i2,j1:j2,ilev))= &
                           0.0_dp
                   ELSEWHERE
                      XCALC(i)%result3d(1)%ptr(_RI_XYZ__(i1:i2,j1:j2,ilev))= &
                           0.0_dp
                      XCALC(i)%result3d(2)%ptr(_RI_XYZ__(i1:i2,j1:j2,ilev))= &
                           XCALC(i)%dat3d(1)%ptr(_RI_XYZ__(i1:i2,j1:j2,ilev))
                   END WHERE

                ELSE

                   ! Error
                   CALL error_bi('The parameter for SPLITV is not defined', substr)

                END IF
             END DO

          CASE DEFAULT !-------------------------------------------------------
             ! This should never be reached!
             CALL error_bi('Something went horribly wrong in SCALC (E0666)' &
                  , substr)

          END SELECT

       END IF is_ok

    END DO calc_loop

  END SUBROUTINE scalc_integrate
  ! ---------------------------------------------------------------------

  ! ---------------------------------------------------------------------
  SUBROUTINE scalc_free_memory

    ! SCALC MODULE ROUTINE
    !
    ! deallocate/clean memory
    !

    IMPLICIT NONE

    INTRINSIC :: ASSOCIATED

    ! LOCAL
    INTEGER :: i

    calc_loop: DO i=1, NCALC
       !
       IF (ASSOCIATED(XCALC(i)%cha)) DEALLOCATE(XCALC(i)%cha)
       NULLIFY(XCALC(i)%cha)
       IF (ASSOCIATED(XCALC(i)%obj)) DEALLOCATE(XCALC(i)%obj)
       NULLIFY(XCALC(i)%obj)
       IF (ASSOCIATED(XCALC(i)%tmp_fac)) DEALLOCATE(XCALC(i)%tmp_fac)
       NULLIFY(XCALC(i)%tmp_fac)

       CALL scalc_freemem(i)
    END DO calc_loop

    IF (ASSOCIATED(XCALC)) DEALLOCATE(XCALC)
    NULLIFY(XCALC)

  END SUBROUTINE scalc_free_memory
  ! ---------------------------------------------------------------------

  ! ======================================================================
  ! PRIVATE ROUTINES
  ! ======================================================================

  ! ----------------------------------------------------------------------
  SUBROUTINE scalc_read_nml_cpl(status, iou)

    ! SCALC MODULE ROUTINE
    !
    ! read namelist for 'coupling'
    !
    ! Author: Patrick Joeckel, MPICH, Jan 2007

    USE messy_main_tools, ONLY: read_nml_open, read_nml_check, read_nml_close

    IMPLICIT NONE

    ! I/O
    INTEGER, INTENT(OUT) :: status     ! error status
    INTEGER, INTENT(IN)  :: iou        ! I/O unit

    ! (LOCAL) NAMELIST VARIABLES
    CHARACTER(LEN=*), PARAMETER :: substr = 'scalc_read_nml_cpl'

    NAMELIST /CPL/ CALC

    LOGICAL :: lex   ! file exists?
    INTEGER :: fstat ! file status

    status = 1

    CALL read_nml_open(lex, substr, iou, 'CPL', modstr)
    IF (.not. lex) RETURN ! namelist file (<modstr>.nml) not available

    READ (iou, NML=CPL, IOSTAT=fstat)
    CALL read_nml_check(fstat, substr, iou, 'CPL', modstr)
    IF (fstat /= 0) RETURN ! error while reading namelist

    CALL read_nml_close(substr, iou, modstr)

    status = 0 ! ...done, without error


  END SUBROUTINE scalc_read_nml_cpl
  ! ----------------------------------------------------------------------

  ! ----------------------------------------------------------------------
  SUBROUTINE str2chobfac(status, str, n, c, o, f, l)

    ! from messy_main_tools.f90 (SUBROUTINE str2chob)
    ! extended to split namelist string "channel:object%factor"
    ! in channel, object and factor string

    USE messy_main_grid_def_mem_bi,   ONLY: nlev
    USE messy_main_tools,     ONLY: strcrack, str2num

    IMPLICIT NONE
    INTRINSIC :: ASSOCIATED, TRIM

    ! I/O
    INTEGER,           INTENT(OUT)               :: status
    CHARACTER(LEN=*),  INTENT(IN)                :: str
    INTEGER,           INTENT(OUT)               :: n
    CHARACTER(LEN=*),  DIMENSION(:), POINTER     :: c
    CHARACTER(LEN=*),  DIMENSION(:), POINTER     :: o
    REAL(DP),          DIMENSION(:), POINTER     :: f
    INTEGER,           DIMENSION(:), POINTER     :: l

    ! LOCAL
    INTEGER                                            :: nc, no, nf, two, j, i
    INTEGER                                            :: nl
    CHARACTER(LEN=LEN(str)), DIMENSION(:), POINTER :: tmp1 !=> NULL()
    CHARACTER(LEN=LEN(str)), DIMENSION(:), POINTER :: tmp2 !=> NULL()
    CHARACTER(LEN=LEN(str)), DIMENSION(:), POINTER :: tmp3 !=> NULL()
    CHARACTER(LEN=LEN(str)), DIMENSION(:), POINTER :: tmp4 !=> NULL()
    CHARACTER(LEN=LEN(str)), DIMENSION(:), POINTER :: tmp5 !=> NULL()

    ! INIT
    IF (ASSOCIATED(c)) DEALLOCATE(c)
    IF (ASSOCIATED(o)) DEALLOCATE(o)
    IF (ASSOCIATED(f)) DEALLOCATE(f)
    IF (ASSOCIATED(l)) DEALLOCATE(l)
    NULLIFY(c)
    NULLIFY(o)
    NULLIFY(f)
    NULLIFY(l)
    n = 0
    NULLIFY(tmp1)
    NULLIFY(tmp2)
    NULLIFY(tmp3)
    NULLIFY(tmp4)
    NULLIFY(tmp5)

    status = 1 ! ERROR

    ! COUNT OBJECTS
    CALL strcrack(TRIM(str), ';', tmp1, nc)           ! CHANNEL BLOCKS
    DO i=1, nc
       CALL strcrack(TRIM(tmp1(i)), ':', tmp2, two)   ! ONE CHANNEL PER BLOCK
       IF (two /= 2) RETURN      ! ERROR: 0 or more than 1 ':' in string
       CALL strcrack(TRIM(tmp2(2)), ',', tmp3, no)    ! OBJECTS PER CHANNEL
       n = n + no
    END DO

    ! ALLOCATE SPACE
    ALLOCATE(c(n))
    ALLOCATE(o(n))
    ALLOCATE(f(n))
    ALLOCATE(l(n))
    ! INIT
    DO i=1,n
       c(i) = ''
       o(i) = ''
       f(i) = 1.0_dp
       l(i) = -99
    END DO

    ! PARSE STRING
    n = 0
    CALL strcrack(TRIM(str), ';', tmp1, nc)           ! CHANNEL BLOCKS
    DO i=1, nc
       CALL strcrack(TRIM(tmp1(i)), ':', tmp2, two)   ! ONE CHANNEL PER BLOCK
       CALL strcrack(TRIM(tmp2(2)), ',', tmp3, no)    ! OBJECTS PER CHANNEL
       DO j=1, no
          n = n + 1
          c(n) = TRIM(tmp2(1))
          CALL strcrack(TRIM(tmp3(j)), '%', tmp4, nf)
          IF (nf == 1) THEN
             CALL strcrack(TRIM(tmp3(j)), '&', tmp5, nl)
             IF (nl == 1) THEN
                o(n) = TRIM(tmp5(1))
             ELSE
                IF (nl == 2) THEN
                   o(n) = TRIM(tmp5(1))
                   IF (TRIM(ADJUSTL(tmp5(2))) == 'LL') THEN
                      l(n) = nlev
                   ELSE
                      CALL str2num(tmp5(2), l(n), status)
                      IF (status /= 0) RETURN
                   END IF
                ELSE
                   RETURN
                END IF
             END IF
          ELSE
             IF (nf == 2) THEN
                o(n) = TRIM(tmp4(1))
                ! search for vertical rank reduction
                CALL strcrack(TRIM(tmp4(2)), '&', tmp5, nl)
                IF (nl == 1) THEN
                   CALL str2num(tmp4(2), f(n), status)
                   IF (status /= 0) RETURN
                ELSE
                   IF (nl == 2) THEN
                      CALL str2num(tmp5(1), f(n), status)
                      IF (status /= 0) RETURN
                      IF (TRIM(ADJUSTL(tmp5(2))) == 'LL') THEN
                         l(n) = nlev
                      ELSE
                         CALL str2num(tmp5(2), l(n), status)
                         IF (status /= 0) RETURN
                      END IF
                   ELSE
                      RETURN
                   END IF
                END IF
             ELSE
                RETURN
             END IF
          END IF
       END DO
    END DO

    IF (ASSOCIATED(tmp1)) THEN
       DEALLOCATE(tmp1) ; NULLIFY(tmp1)
    ENDIF
    IF (ASSOCIATED(tmp2)) THEN
       DEALLOCATE(tmp2) ; NULLIFY(tmp2)
    ENDIF
    IF (ASSOCIATED(tmp3)) THEN
       DEALLOCATE(tmp3) ; NULLIFY(tmp3)
    ENDIF
    IF (ASSOCIATED(tmp4)) THEN
       DEALLOCATE(tmp4) ; NULLIFY(tmp4)
    ENDIF
    IF (ASSOCIATED(tmp5)) THEN
       DEALLOCATE(tmp5) ; NULLIFY(tmp5)
    ENDIF

    status  = 0  ! NO ERROR

  END SUBROUTINE str2chobfac
  ! ----------------------------------------------------------------------

  ! ---------------------------------------------------------------------
  SUBROUTINE scalc_freemem(i)

    ! SCALC MODULE ROUTINE
    !
    ! deallocate/clean memory
    !

    IMPLICIT NONE

    INTRINSIC :: ASSOCIATED

    ! LOCAL
    INTEGER, INTENT(in) :: i

    IF (ASSOCIATED(XCALC(i)%io%when))  DEALLOCATE(XCALC(i)%io%when)
    NULLIFY(XCALC(i)%io%when)
    IF (ASSOCIATED(XCALC(i)%io%where)) DEALLOCATE(XCALC(i)%io%where)
    NULLIFY(XCALC(i)%io%where)

    IF (ASSOCIATED(XCALC(i)%lex)) DEALLOCATE(XCALC(i)%lex)
    NULLIFY(XCALC(i)%lex)
    !
    IF (ASSOCIATED(XCALC(i)%dat2d)) DEALLOCATE(XCALC(i)%dat2d)
    NULLIFY(XCALC(i)%dat2d)
    IF (ASSOCIATED(XCALC(i)%dat3d)) DEALLOCATE(XCALC(i)%dat3d)
    NULLIFY(XCALC(i)%dat3d)
    IF (ASSOCIATED(XCALC(i)%rank))   DEALLOCATE(XCALC(i)%rank)
    NULLIFY(XCALC(i)%rank)
    IF (ASSOCIATED(XCALC(i)%rid)) DEALLOCATE(XCALC(i)%rid)
    NULLIFY(XCALC(i)%rid)
    IF (ASSOCIATED(XCALC(i)%lev))   DEALLOCATE(XCALC(i)%lev)
    NULLIFY(XCALC(i)%lev)
    IF (ASSOCIATED(XCALC(i)%xparam))   DEALLOCATE(XCALC(i)%xparam)
    NULLIFY(XCALC(i)%xparam)
    IF (ASSOCIATED(XCALC(i)%result2d)) DEALLOCATE(XCALC(i)%result2d)
    NULLIFY(XCALC(i)%result2d)
    IF (ASSOCIATED(XCALC(i)%result3d)) DEALLOCATE(XCALC(i)%result3d)
    NULLIFY(XCALC(i)%result3d)
    IF (ASSOCIATED(XCALC(i)%fac))    DEALLOCATE(XCALC(i)%fac)
    NULLIFY(XCALC(i)%fac)
    IF (ASSOCIATED(XCALC(i)%diag)) DEALLOCATE(XCALC(i)%diag)
    NULLIFY(XCALC(i)%diag)
    IF (ASSOCIATED(XCALC(i)%unit)) DEALLOCATE(XCALC(i)%unit)
    NULLIFY(XCALC(i)%unit)

  END SUBROUTINE scalc_freemem
  ! ---------------------------------------------------------------------

  ! ---------------------------------------------------------------------
  SUBROUTINE parse_when (status, str_in, num, int_when, int_where)

    ! AUTHOR: Astrid Kerkweg, March 2019

    USE messy_main_control,    ONLY: control_get_idx_from_name   &
                                   , MEP_GLOBAL_END, MEPS_ALWAYS &
                                   , MEPS_MAX, MEPS_MIN
    USE messy_main_tools,      ONLY: strcrack, str2num

    IMPLICIT NONE

    ! IO
    INTEGER,                        INTENT(OUT) :: status
    CHARACTER(LEN=*),               INTENT(IN)  :: str_in
    INTEGER,                        INTENT(OUT) :: num
    INTEGER, DIMENSION(:), POINTER, INTENT(OUT) :: int_when
    INTEGER, DIMENSION(:), POINTER, INTENT(OUT) :: int_where

    ! LOCAL
    CHARACTER(LEN=strlen), POINTER :: sl1(:) => NULL()
    CHARACTER(LEN=strlen), POINTER :: sl2(:) => NULL()
    INTEGER                        :: m, k
    INTEGER                        :: i

    status = -1

    NULLIFY(sl1)

    CALL STRCRACK(str_in, ';', sl1, num)

    ALLOCATE(int_when(num))
    int_when  = MEP_GLOBAL_END
    ALLOCATE(int_where(num))
    int_where = MEPS_ALWAYS

    DO i=1, num

       CALL STRCRACK(sl1(i), ',', sl2, m)

       IF (SIZE(sl2) == 2) THEN
          IF (TRIM(ADJUSTL(sl2(2))) == '') THEN
             write (*,*) 'PARSE ENTRY POINT:  SUBPOSITION MISSING DO ALWAYS'
          ELSE
             CALL str2num(sl2(2),k, status)
             IF (status /= 0) THEN
                write (*,*) 'ERROR reading subentry of SCALC set',i
                RETURN
             ELSE
                IF (k>=MEPS_MIN .AND. k<=MEPS_MAX) THEN
                   int_where(i) = k
                ELSE
                   write (*,*) 'ERROR reading subentry of SCALC set',i, k
                END IF
             END IF
          END IF
       ELSE
          CALL control_get_idx_from_name(status, int_when(i), sl2(1))
          IF (status /= 0) THEN
             write  (*,*) 'ERROR evaluating entry of SCALC set',i, sl2(1)
             RETURN
          END IF
       END IF

       IF (ASSOCIATED(sl2)) THEN
          DEALLOCATE(sl2) ; NULLIFY(sl2)
       END IF
    END DO

    IF (ASSOCIATED(sl1)) THEN
       DEALLOCATE(sl1) ; NULLIFY(sl1)
    END IF


  END SUBROUTINE parse_when
  ! ---------------------------------------------------------------------

  ! ---------------------------------------------------------------------
  SUBROUTINE get_param(status, str, op, p, np)

    USE messy_main_tools,      ONLY: strcrack, str2num

    IMPLICIT NONE

    INTRINSIC :: TRIM, ADJUSTL

    ! I/O
    INTEGER,                INTENT(OUT) :: status
    CHARACTER(LEN=*),       INTENT(IN)  :: str
    CHARACTER(LEN=*),       INTENT(OUT) :: op
    REAL(DP), DIMENSION(:), POINTER     :: p   ! INTENT(OUT)
    INTEGER,                INTENT(OUT) :: np

    ! LOCAL
    CHARACTER(LEN=STRLEN_MEDIUM), POINTER :: sl1(:) => NULL()
    CHARACTER(LEN=STRLEN_MEDIUM), POINTER :: sl2(:) => NULL()
    CHARACTER(LEN=STRLEN_MEDIUM), POINTER :: sl3(:) => NULL()
    INTEGER :: m, i

    status = 1

    IF (ASSOCIATED(p)) THEN
       DEALLOCATE(p); NULLIFY(p)
    END IF
    np = 0

    CALL strcrack(str, '(', sl1, m)
    IF (m/=2) RETURN
    op = TRIM(ADJUSTL(sl1(1)))

    CALL strcrack(TRIM(ADJUSTL(sl1(2))), ')', sl2, m)
    IF (m/=1) RETURN

    CALL strcrack(TRIM(ADJUSTL(sl2(1))), ',', sl3, m)
    ALLOCATE(p(m))
    np = m
    DO i=1,m
       CALL str2num(TRIM(ADJUSTL(sl3(i))), p(i), status)
       IF (status /=0) RETURN
    END DO

    IF (ASSOCIATED(sl1)) THEN
       DEALLOCATE(sl1) ; NULLIFY(sl1)
    END IF

    IF (ASSOCIATED(sl2)) THEN
       DEALLOCATE(sl2) ; NULLIFY(sl2)
    END IF

    IF (ASSOCIATED(sl3)) THEN
       DEALLOCATE(sl3) ; NULLIFY(sl3)
    END IF

    status = 0

  END SUBROUTINE get_param
  ! ---------------------------------------------------------------------

  ! **********************************************************************
END MODULE messy_scalc_si
! **********************************************************************

! ************************************************************************
MODULE messy_lggp_e5
! ************************************************************************

  ! AUTHOR:
  !  Patrick Joeckel, MPICH, December 2005

  ! ECHAM5/MESSy
  USE messy_main_blather_bi,    ONLY: start_message_bi, end_message_bi
  ! MESSy
  USE messy_main_channel,       ONLY: STRLEN_CHANNEL, STRLEN_OBJECT
  USE messy_lggp
  USE messy_attila_tools,       ONLY: LG2GP_SUM, LG2GP_AVE, LG2GP_STD &
                                    , LG2GP_AVEGT0

  IMPLICIT NONE
  PRIVATE

  INTRINSIC :: NULL

  ! #####################
  ! LAGRANGE -> GRIDPOINT
  ! #####################

  INTEGER, PARAMETER :: FILL_NONE  = 0
  INTEGER, PARAMETER :: FILL_VALUE = 1
  INTEGER, PARAMETER :: FILL_FIELD = 2

  TYPE IO_LG2GP
     CHARACTER(LEN=STRLEN_OBJECT)  :: name = ''
     CHARACTER(LEN=STRLEN_CHANNEL) :: lg_channel = ''
     CHARACTER(LEN=STRLEN_OBJECT)  :: lg_object = ''
     INTEGER                       :: method = LG2GP_AVE
     LOGICAL                       :: lmasscons = .FALSE.
     INTEGER                       :: fill_flag = FILL_NONE
     REAL(DP)                      :: fill_value = 0.0_dp
     CHARACTER(LEN=STRLEN_CHANNEL) :: fill_channel_gp = ''
     CHARACTER(LEN=STRLEN_OBJECT)  :: fill_object_gp = ''
  END TYPE IO_LG2GP

  TYPE T_LG2GP
     TYPE(IO_LG2GP) :: io
     REAL(DP), DIMENSION(:),     POINTER :: lg  => NULL() ! LAGRANGE
     REAL(DP), DIMENSION(:,:,:), POINTER :: gp  => NULL() ! GRIDPOINT
     REAL(DP), DIMENSION(:,:,:), POINTER :: gpf => NULL() ! FILL (optional)
     LOGICAL                             :: ok = .FALSE.
  END TYPE T_LG2GP

  INTEGER, PARAMETER :: NMAXLG2GP = 100
  TYPE(IO_LG2GP), DIMENSION(NMAXLG2GP), SAVE :: LG2GP
  INTEGER,                              SAVE :: NLG2GP = 0
  TYPE(T_LG2GP),  DIMENSION(NMAXLG2GP), SAVE :: XLG2GP

  ! #####################
  ! GRIDPOINT -> LAGRANGE
  ! #####################

  TYPE IO_GP2LG
     CHARACTER(LEN=STRLEN_OBJECT)  :: name = ''
     CHARACTER(LEN=STRLEN_CHANNEL) :: gp_channel = ''
     CHARACTER(LEN=STRLEN_OBJECT)  :: gp_object = ''
     LOGICAL                       :: lmasscons = .FALSE.
     LOGICAL                       :: lrest = .FALSE.
  END TYPE IO_GP2LG

  TYPE T_GP2LG
     TYPE(IO_GP2LG) :: io
     REAL(DP), DIMENSION(:,:,:), POINTER :: gp  => NULL() ! GRIDPOINT
     REAL(DP), DIMENSION(:),     POINTER :: lg  => NULL() ! LAGRANGE
     REAL(DP), DIMENSION(:,:,:), POINTER :: gpr => NULL() ! GRIDPOINT REST
     LOGICAL                             :: ok = .FALSE.
  END TYPE T_GP2LG
  
  INTEGER, PARAMETER :: NMAXGP2LG = 100
  TYPE(IO_GP2LG), DIMENSION(NMAXGP2LG), SAVE :: GP2LG
  INTEGER,                              SAVE :: NGP2LG = 0
  TYPE(T_GP2LG),  DIMENSION(NMAXGP2LG), SAVE :: XGP2LG

  ! #####################

  LOGICAL :: LRUNSM = .TRUE.  ! RUN THIS SUBMODEL

  ! #####################

  PUBLIC :: lggp_initialize
  PUBLIC :: lggp_init_coupling
  PUBLIC :: lggp_global_end
  !PRIVATE :: lggp_read_nml_cpl

CONTAINS

  ! ######################################################################
  ! PUBLIC SUBROUTINES
  ! ######################################################################
  
  ! ---------------------------------------------------------------------
  SUBROUTINE lggp_initialize

    ! lggp MODULE ROUTINE (ECHAM5 INTERFACE)
    !
    ! INITIALIZATION OF GLOBAL VARIABLES FROM NAMELIST
    ! IN PARALLEL ENVIRONMENT
    !
    ! Author: Patrick Joeckel, MPICH, Dec 2005
    
    ! ECHAM5/MESSy
    USE messy_main_mpi_bi,     ONLY: p_parallel_io, p_io, p_bcast, finish
    USE messy_main_channel_bi, ONLY: LG_ATTILA, REPR_UNDEF
    USE messy_main_tools,      ONLY: find_next_free_unit

    IMPLICIT NONE

    INTRINSIC :: TRIM

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'lggp_initialize'
    INTEGER                :: iou    ! I/O unit
    INTEGER                :: status ! error status
    INTEGER                :: i

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
       CALL lggp_read_nml_cpl(status, iou)
       IF (status /= 0) CALL finish(substr)
    END IF

    IF (p_parallel_io) THEN

       ! ### (1) LG2GP ###
       ! GET NUMBER OF ENTRIES
       WRITE(*,*) '------------------------------------------------------'
       WRITE(*,*) 'TRANSFORMATION LAGRANGE -> GRIDPOINT:'
       ! COPY DATA
       NLG2GP = 1
       DO i=1, NMAXLG2GP

          IF (TRIM(LG2GP(i)%name) == '') CYCLE
          XLG2GP(NLG2GP)%io%name       = TRIM(LG2GP(i)%name)

          IF (TRIM(LG2GP(i)%lg_channel) == '') CYCLE
          XLG2GP(NLG2GP)%io%lg_channel = TRIM(LG2GP(i)%lg_channel)

          IF (TRIM(LG2GP(i)%lg_object) == '') CYCLE
          XLG2GP(NLG2GP)%io%lg_object  = TRIM(LG2GP(i)%lg_object)

          WRITE(*,*) '  ', TRIM(XLG2GP(NLG2GP)%io%name),' <- ', &
               TRIM(XLG2GP(NLG2GP)%io%lg_channel), &
               '(',TRIM(XLG2GP(NLG2GP)%io%lg_object),') ...'

          SELECT CASE (LG2GP(i)%method)
          CASE(LG2GP_SUM)
             WRITE(*,*) '   ... SUM'
          CASE(LG2GP_AVE)
             WRITE(*,*) '   ... AVERAGE'
          CASE(LG2GP_STD)
             WRITE(*,*) '   ... STANDARD DEVIATION'
          CASE(LG2GP_AVEGT0)
             WRITE(*,*) '   ... AVERAGE ALL ELEMENTS > 0'
          CASE DEFAULT
             WRITE(*,*) '   ... UNKNOWN TRANSFORMATION METHOD '//&
                  &'... skipping'
             CYCLE
          END SELECT
          XLG2GP(NLG2GP)%io%method = LG2GP(i)%method
          
          XLG2GP(NLG2GP)%io%lmasscons = LG2GP(i)%lmasscons
          IF (XLG2GP(NLG2GP)%io%lmasscons) THEN
             WRITE(*,*) '   ... mass conserving transformtion: YES'
          ELSE
             WRITE(*,*) '   ... mass conserving transformtion: NO'
          END IF

          SELECT CASE (LG2GP(i)%fill_flag)
          CASE(FILL_NONE)
             WRITE(*,*) '   ... no filling'
          CASE(FILL_VALUE)
             WRITE(*,*) '   ... filling with const. value: ' &
                  , LG2GP(i)%fill_value
          CASE(FILL_FIELD)
             WRITE(*,*) '   ... filling with GP-field: ', &
                  TRIM(LG2GP(i)%fill_channel_gp) &
                  ,'(',TRIM(LG2GP(i)%fill_object_gp),')'             
          CASE DEFAULT
             WRITE(*,*) '   ... UNKNOWN FILL FLAG VALUE ... skipping'
             CYCLE
          END SELECT
          XLG2GP(NLG2GP)%io%fill_flag = LG2GP(i)%fill_flag
          XLG2GP(NLG2GP)%io%fill_value = LG2GP(i)%fill_value
          XLG2GP(NLG2GP)%io%fill_channel_gp = LG2GP(i)%fill_channel_gp
          XLG2GP(NLG2GP)%io%fill_object_gp = LG2GP(i)%fill_object_gp

          ! NEXT ENTRY
          NLG2GP = NLG2GP + 1
          WRITE(*,*) '------------------------------------------------------'

       END DO
       NLG2GP = NLG2GP - 1

       ! ### (2) GP2LG ###
       ! GET NUMBER OF ENTRIES
       WRITE(*,*) '------------------------------------------------------'
       WRITE(*,*) 'TRANSFORMATION GRIDPOINT -> LAGRANGE:'
       ! COPY DATA
       NGP2LG = 1
       DO i=1, NMAXGP2LG

          IF (TRIM(GP2LG(i)%name) == '') CYCLE
          XGP2LG(NGP2LG)%io%name = TRIM(GP2LG(i)%name)

          IF (TRIM(GP2LG(i)%gp_channel) == '') CYCLE
          XGP2LG(NGP2LG)%io%gp_channel = TRIM(GP2LG(i)%gp_channel)

          IF (TRIM(GP2LG(i)%gp_object) == '') CYCLE
          XGP2LG(NGP2LG)%io%gp_object = TRIM(GP2LG(i)%gp_object)

          WRITE(*,*) '  ', TRIM(XGP2LG(NGP2LG)%io%name),' <- ', &
               TRIM(XGP2LG(NGP2LG)%io%gp_channel), &
               '(',TRIM(XGP2LG(NGP2LG)%io%gp_object),') ...'

          XGP2LG(NGP2LG)%io%lmasscons = GP2LG(i)%lmasscons
          IF (XGP2LG(NGP2LG)%io%lmasscons) THEN
             WRITE(*,*) '   ... mass conserving transformtion: YES'
             XGP2LG(NGP2LG)%io%lrest = GP2LG(i)%lrest
             IF (XGP2LG(NGP2LG)%io%lrest) THEN
                WRITE(*,*) '       (with accouning rest)'
             ELSE
                WRITE(*,*) '       (without accouning rest)'
             END IF
          ELSE
             WRITE(*,*) '   ... mass conserving transformtion: NO'
             XGP2LG(NGP2LG)%io%lrest = .FALSE.
          END IF

          ! NEXT ENTRY
          NGP2LG = NGP2LG + 1
          WRITE(*,*) '------------------------------------------------------'

       END DO
       NGP2LG = NGP2LG - 1

    END IF

    ! BROADCAST RESULTS
    CALL p_bcast(NLG2GP, p_io)
    CALL p_bcast(NGP2LG, p_io)


    DO i=1, NLG2GP
       CALL p_bcast(XLG2GP(i)%io%name, p_io)
       CALL p_bcast(XLG2GP(i)%io%lg_channel, p_io)
       CALL p_bcast(XLG2GP(i)%io%lg_object, p_io)
       CALL p_bcast(XLG2GP(i)%io%method, p_io)
       CALL p_bcast(XLG2GP(i)%io%lmasscons, p_io)
       CALL p_bcast(XLG2GP(i)%io%fill_flag, p_io)
       CALL p_bcast(XLG2GP(i)%io%fill_value, p_io)
       CALL p_bcast(XLG2GP(i)%io%fill_channel_gp, p_io)
       CALL p_bcast(XLG2GP(i)%io%fill_object_gp, p_io)
    END DO

    DO i=1, NGP2LG 
       CALL p_bcast(XGP2LG(i)%io%name, p_io)
       CALL p_bcast(XGP2LG(i)%io%gp_channel, p_io)
       CALL p_bcast(XGP2LG(i)%io%gp_object, p_io)
       CALL p_bcast(XGP2LG(i)%io%lmasscons, p_io)
       CALL p_bcast(XGP2LG(i)%io%lrest, p_io)
    END DO

    IF (p_parallel_io) THEN
       WRITE(*,*) ' ---> ',NLG2GP+NGP2LG,' TRANSFORMATIONS(S) INITIALIZED !'
    END IF

    LRUNSM = (NLG2GP + NGP2LG > 0)

    CALL end_message_bi(modstr, 'INITIALISATION', substr)

  END SUBROUTINE lggp_initialize
  ! ---------------------------------------------------------------------

  ! ---------------------------------------------------------------------
  SUBROUTINE lggp_init_coupling

    ! ECHAM5/MESSy
    USE messy_main_mpi_bi,           ONLY: p_parallel_io, message
    USE messy_main_channel_error_bi, ONLY: channel_halt
    USE messy_main_channel_bi,        ONLY: GP_3D_MID, LG_ATTILA

    ! MESSy
    USE messy_main_channel,         ONLY: new_channel, new_channel_object &
                                        , new_attribute, get_attribute    &
                                        , get_channel_object              &
                                        , get_channel_object_info         &
                                        , get_channel_info
    USE messy_main_constants_mem,   ONLY: STRLEN_ULONG

    IMPLICIT NONE

    INTRINSIC :: TRIM

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'lggp_init_coupling'
    INTEGER                     :: status
    INTEGER                     :: i
    INTEGER                     :: reprid
    LOGICAL                     :: lfirst
    CHARACTER(LEN=STRLEN_ULONG) :: unit
    REAL(DP)                    :: molarmass
    LOGICAL, SAVE               :: lcalled = .FALSE.

    IF (.NOT.LRUNSM) RETURN

    CALL start_message_bi(modstr, 'CHANNEL DEFINITION', substr)

    ! (1) LG -> GP
    IF (p_parallel_io) THEN
       WRITE(*,*) '------------------------------------------------------'
       WRITE(*,*) 'TRANSFORMATION LAGRANGE -> GRIDPOINT:'
    END IF
    IF (lcalled) THEN
       CALL get_channel_info(status, modstr//'_gp')
       lfirst = (status /= 0)
    ELSE
       lfirst = .TRUE.
    ENDIF
    loop_lg2gp: DO i=1, NLG2GP

       IF (XLG2GP(i)%ok) CYCLE
       
       ! CHECK AVAILABILITY OF LG CHANNEL/OBJECT
       CALL get_channel_object(status &
            , TRIM(XLG2GP(i)%io%lg_channel), TRIM(XLG2GP(i)%io%lg_object) &
            , p1=XLG2GP(i)%lg &
            )
       IF (status /= 0) THEN
          CALL message('  ',TRIM(XLG2GP(i)%io%lg_channel)//' - '//&
               &TRIM(XLG2GP(i)%io%lg_object)//' not found ... skipping' )
          CYCLE
       ELSE
          CALL message('  ',TRIM(XLG2GP(i)%io%lg_channel)//' - '//&
               &TRIM(XLG2GP(i)%io%lg_object)//' found ...' )
          CALL get_channel_object_info(status &
               , TRIM(XLG2GP(i)%io%lg_channel), TRIM(XLG2GP(i)%io%lg_object) &
               , reprid=reprid)
          IF (reprid /= LG_ATTILA) THEN
             CALL message('  ',' ... wrong representation ... skipping')
             CYCLE
          ELSE
             CALL message('  ',' ... LG_ATTILA representation ... ')
          END IF
       END IF

       ! CHECK AVAILABILITY OF GP_3D_MID FIELD FOR FILLING, IF REQUIRED
       IF (XLG2GP(i)%io%fill_flag == FILL_FIELD) THEN
          CALL message('  ',' ... filling with FIELD requested ...')
          CALL get_channel_object(status &
               , TRIM(XLG2GP(i)%io%fill_channel_gp) &
               , TRIM(XLG2GP(i)%io%fill_object_gp) &
               , p3=XLG2GP(i)%gpf &
               )
          IF (status /= 0) THEN
             CALL message('  ', '  ... '//&
                  &TRIM(XLG2GP(i)%io%fill_channel_gp)//' - '//&
                  &TRIM(XLG2GP(i)%io%fill_object_gp)//&
                  &' not found ... skipping')
             CYCLE
          ELSE
             CALL message('  ', '  ... '//&
                  &TRIM(XLG2GP(i)%io%fill_channel_gp)//' - '//&
                  &TRIM(XLG2GP(i)%io%fill_object_gp)//&
                  &' found ... ')
             CALL get_channel_object_info(status &
                  , TRIM(XLG2GP(i)%io%fill_channel_gp) &
                  , TRIM(XLG2GP(i)%io%fill_object_gp) &
                  , reprid=reprid)
             IF (reprid /= GP_3D_MID) THEN
                CALL message('  ',' ... wrong representation ... skipping')
                CYCLE
             ELSE
                CALL message('  ',' ... GP_3D_MID representation ... ')
             END IF
          END IF
       END IF

       ! EVERYTHIG IS OK NOW
       XLG2GP(i)%ok = .TRUE.

       ! CREATE CHANNEL FOR THE FIRST OBJECT
       IF (lfirst) THEN
          lfirst = .FALSE.
          CALL new_channel(status, modstr//'_gp', reprid=GP_3D_MID)
          CALL channel_halt(substr, status)
       END IF

       ! CREATE NEW CHANNEL OBJECT
       CALL new_channel_object(status, modstr//'_gp'    &
            , TRIM(XLG2GP(i)%io%name), p3=XLG2GP(i)%gp)
       CALL channel_halt(substr, status)

       ! ADD ATTRIBUTES
       ! ... UNITS
       CALL get_attribute(status &
            , TRIM(XLG2GP(i)%io%lg_channel), TRIM(XLG2GP(i)%io%lg_object) &
            , 'units', c=unit)
       IF (status == 0) THEN
          CALL new_attribute(status, modstr//'_gp', TRIM(XLG2GP(i)%io%name) &
               , 'units',c = TRIM(unit) )
       END IF
       ! ... MOLARMASS (FOR TRACERS AND LGVFLUX)
       CALL get_attribute(status &
            , TRIM(XLG2GP(i)%io%lg_channel), TRIM(XLG2GP(i)%io%lg_object) &
            , 'molarmass', r=molarmass)
       IF (status == 0) THEN
          CALL new_attribute(status, modstr//'_gp', TRIM(XLG2GP(i)%io%name) &
               , 'molarmass',r = molarmass )
       END IF
       ! ... ORIGIN
       CALL new_attribute(status, modstr//'_gp', TRIM(XLG2GP(i)%io%name) &
            , 'long_name', c='LG2GP transformation of '//&
            &TRIM(XLG2GP(i)%io%lg_channel)//' - '&
            &//TRIM(XLG2GP(i)%io%lg_object) )
       CALL channel_halt(substr, status)
       ! ... METHOD
       SELECT CASE(XLG2GP(i)%io%method)
       CASE(LG2GP_SUM)
          CALL new_attribute(status, modstr//'_gp', TRIM(XLG2GP(i)%io%name) &
               , 'trafo-method', c='LG2GP_SUM')
       CASE(LG2GP_AVE)
          CALL new_attribute(status, modstr//'_gp', TRIM(XLG2GP(i)%io%name) &
               , 'trafo-method', c='LG2GP_AVE')
       CASE(LG2GP_STD)
          CALL new_attribute(status, modstr//'_gp', TRIM(XLG2GP(i)%io%name) &
               , 'trafo-method', c='LG2GP_STD')
       CASE(LG2GP_AVEGT0)
          CALL new_attribute(status, modstr//'_gp', TRIM(XLG2GP(i)%io%name) &
               , 'trafo-method', c='LG2GP_AVEGT0')
       END SELECT
       CALL channel_halt(substr, status)
       ! ... MASS CONSERVATION
       IF (XLG2GP(i)%io%lmasscons) THEN
          CALL new_attribute(status, modstr//'_gp', TRIM(XLG2GP(i)%io%name) &
               , 'mass-conservations', c='YES')
       ELSE
          CALL new_attribute(status, modstr//'_gp', TRIM(XLG2GP(i)%io%name) &
               , 'mass-conservation', c='NO')
       END IF
       CALL channel_halt(substr, status)
       ! ... FILLING
       SELECT CASE(XLG2GP(i)%io%fill_flag)
       CASE(FILL_NONE)
          CALL new_attribute(status, modstr//'_gp', TRIM(XLG2GP(i)%io%name) &
               , 'filling', c='none')
       CASE(FILL_VALUE)
          CALL new_attribute(status, modstr//'_gp', TRIM(XLG2GP(i)%io%name) &
               , 'filling', r=XLG2GP(i)%io%fill_value)
       CASE(FILL_FIELD)
          CALL new_attribute(status, modstr//'_gp', TRIM(XLG2GP(i)%io%name) &
               , 'filling', c=TRIM(XLG2GP(i)%io%fill_channel_gp)//&
               &' - '//TRIM(XLG2GP(i)%io%fill_object_gp) )
       END SELECT
       CALL channel_halt(substr, status)
    END DO loop_lg2gp

    ! (2) GP -> LG
    IF (p_parallel_io) THEN
       WRITE(*,*) '------------------------------------------------------'
       WRITE(*,*) 'TRANSFORMATION GRIDPOINT -> LAGRANGE:'
    END IF
    IF (lcalled) THEN
       CALL get_channel_info(status, modstr//'_lg')
       lfirst = (status /= 0)
    ELSE
       lfirst = .TRUE.
    END IF
    loop_gp2lg: DO i=1, NGP2LG

       IF (XGP2LG(i)%ok) CYCLE

       ! CHECK AVAILABILITY OF GP CHANNEL/OBJECT
       CALL get_channel_object(status &
            , TRIM(XGP2LG(i)%io%gp_channel), TRIM(XGP2LG(i)%io%gp_object) &
            , p3=XGP2LG(i)%gp &
            )
       IF (status /= 0) THEN
          CALL message('  ',TRIM(XGP2LG(i)%io%gp_channel)//' - '//&
               &TRIM(XGP2LG(i)%io%gp_object)//' not found ... skipping' )
          CYCLE
       ELSE
          CALL message('  ',TRIM(XGP2LG(i)%io%gp_channel)//' - '//&
               &TRIM(XGP2LG(i)%io%gp_object)//' found ...' )
          CALL get_channel_object_info(status &
               , TRIM(XGP2LG(i)%io%gp_channel), TRIM(XGP2LG(i)%io%gp_object) &
               , reprid=reprid)
          IF (reprid /= GP_3D_MID) THEN
             CALL message('  ',' ... wrong REPRESENTATON ... skipping')
             CYCLE
          ELSE
             CALL message('  ',' ... GP_3D_MID REPRESENTATON ... ')
          END IF
       END IF

       ! EVERYTHING IS OK NOW
       XGP2LG(i)%ok = .TRUE.

       ! CREATE CHANNEL FOR THE FIRST OBJECT
       IF (lfirst) THEN
          lfirst = .FALSE.
          CALL new_channel(status, modstr//'_lg', reprid=LG_ATTILA)
          CALL channel_halt(substr, status)
       END IF

       ! CREATE NEW CHANNEL OBJECT
       CALL new_channel_object(status, modstr//'_lg'   &
            , TRIM(XGP2LG(i)%io%name), p1=XGP2LG(i)%lg)
       CALL channel_halt(substr, status)

       ! ADD ATTRIBUTES
       ! ... UNITS
       CALL get_attribute(status &
            , TRIM(XGP2LG(i)%io%gp_channel), TRIM(XGP2LG(i)%io%gp_object) &
            , 'units', c=unit)
       IF (status == 0) THEN
          CALL new_attribute(status, modstr//'_lg', TRIM(XGP2LG(i)%io%name) &
               , 'units',c = TRIM(unit) )
       END IF
       ! ... MOLARMASS (FOR TRACER AND LGVFLUX)
       CALL get_attribute(status &
            , TRIM(XGP2LG(i)%io%gp_channel), TRIM(XGP2LG(i)%io%gp_object) &
            , 'molarmass', r=molarmass)
       IF (status == 0) THEN
          CALL new_attribute(status, modstr//'_lg', TRIM(XGP2LG(i)%io%name) &
               , 'molarmass',r = molarmass )
       END IF
       ! ... ORIGIN
       CALL new_attribute(status, modstr//'_lg', TRIM(XGP2LG(i)%io%name) &
            , 'long_name', c='GP2LG transformation of '//&
            &TRIM(XGP2LG(i)%io%gp_channel)//' - '&
            &//TRIM(XGP2LG(i)%io%gp_object) )
       CALL channel_halt(substr, status)
       ! ... MASS CONSERVATION
       IF (XGP2LG(i)%io%lmasscons) THEN
          CALL new_attribute(status, modstr//'_lg', TRIM(XGP2LG(i)%io%name) &
               , 'mass-conservations', c='YES')
          IF (XGP2LG(i)%io%lrest) THEN
             CALL new_channel_object(status, modstr//'_lg' &
                  , TRIM(XGP2LG(i)%io%name)//'_rest', p3=XGP2LG(i)%gpr &
                  , reprid=GP_3D_MID, lrestreq=.TRUE.)
       CALL channel_halt(substr, status)
          END IF
       ELSE
          CALL new_attribute(status, modstr//'_lg', TRIM(XGP2LG(i)%io%name) &
               , 'mass-conservation', c='NO')
       END IF
       CALL channel_halt(substr, status)

    END DO loop_gp2lg

    lcalled = .TRUE.

    CALL end_message_bi(modstr, 'CHANNEL DEFINITION', substr)

  END SUBROUTINE lggp_init_coupling
  ! ---------------------------------------------------------------------

  ! ---------------------------------------------------------------------
  SUBROUTINE lggp_global_end

    USE messy_attila_tools_e5, ONLY: gp2lg_e5, lg2gp_e5

    IMPLICIT NONE

    ! LOCAL
    !CHARACTER(LEN=*), PARAMETER :: substr = 'lggp_global_end'
    INTEGER :: i

    IF (.NOT.LRUNSM) RETURN

    ! PERFORM TRANSFORMATIONS
    ! (1) LG -> GP
    DO i=1, NLG2GP
       IF (.NOT. XLG2GP(i)%ok) CYCLE

       SELECT CASE(XLG2GP(i)%io%fill_flag)
       CASE(FILL_NONE)
          CALL lg2gp_e5(XLG2GP(i)%lg, XLG2GP(i)%gp, XLG2GP(i)%io%method &
               , XLG2GP(i)%io%lmasscons, ltm1=.FALSE.)          
       CASE(FILL_VALUE)
          CALL lg2gp_e5(XLG2GP(i)%lg, XLG2GP(i)%gp, XLG2GP(i)%io%method &
               , XLG2GP(i)%io%lmasscons, ltm1=.FALSE. &
               , fill_value = XLG2GP(i)%io%fill_value)
       CASE(FILL_FIELD)
          CALL lg2gp_e5(XLG2GP(i)%lg, XLG2GP(i)%gp, XLG2GP(i)%io%method &
               , XLG2GP(i)%io%lmasscons, ltm1=.FALSE. &
               , fill_field = XLG2GP(i)%gpf)
       END SELECT
    END DO

    ! (2) GP -> LG
    DO i=1, NGP2LG
       IF (.NOT. XGP2LG(i)%ok) CYCLE

       IF (XGP2LG(i)%io%lmasscons) THEN
          IF (XGP2LG(i)%io%lrest) THEN
             CALL gp2lg_e5(XGP2LG(i)%gp, XGP2LG(i)%lg &
                  , gprl = XGP2LG(i)%gpr              &
                  , lmcons=XGP2LG(i)%io%lmasscons)
          ELSE
             CALL gp2lg_e5(XGP2LG(i)%gp, XGP2LG(i)%lg &
                  , lmcons=XGP2LG(i)%io%lmasscons)
          END IF
       ELSE
          CALL gp2lg_e5(XGP2LG(i)%gp, XGP2LG(i)%lg &
               , lmcons=XGP2LG(i)%io%lmasscons)
       END IF
    END DO

  END SUBROUTINE lggp_global_end
  ! ---------------------------------------------------------------------

  ! ######################################################################
  ! PRIVATE SUBROUTINES
  ! ######################################################################

  ! ---------------------------------------------------------------------
  SUBROUTINE lggp_read_nml_cpl(status, iou)

    ! lggp MODULE ROUTINE (ECHAM5 INTERFACE, PRIVATE)
    !
    ! read namelist for 'coupling' to ECHAM5
    !
    ! Author: Patrick Joeckel, MPICH, Dec 2005

    ! MESSy
    USE messy_main_tools, ONLY: read_nml_open, read_nml_check, read_nml_close

    IMPLICIT NONE

    ! I/O
    INTEGER, INTENT(OUT) :: status     ! error status
    INTEGER, INTENT(IN)  :: iou        ! I/O unit

    ! (LOCAL) NAMELIST VARIABLES
    CHARACTER(LEN=*), PARAMETER :: substr = 'lggp_read_nml_cpl'

    NAMELIST /CPL/ LG2GP, GP2LG

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

  END SUBROUTINE lggp_read_nml_cpl
  ! ---------------------------------------------------------------------

! ************************************************************************
END MODULE messy_lggp_e5
! ************************************************************************

! ============================================================================
MODULE messy_main_import_ts_bi
! ============================================================================

  USE messy_main_blather_bi,    ONLY: start_message_bi, end_message_bi &
                                    , error_bi
  USE messy_main_import_ts

  IMPLICIT NONE
  INTRINSIC :: NULL
  PRIVATE
  SAVE

  PUBLIC :: import_ts_initialize
  PUBLIC :: import_ts_init_memory
  PUBLIC :: import_ts_global_start
  PUBLIC :: import_ts_free_memory

  INTEGER, DIMENSION(:), POINTER :: list_dimid  => NULL()
  INTEGER, DIMENSION(:), POINTER :: list_reprid => NULL()

CONTAINS

  ! ----------------------------------------------------------------------
  SUBROUTINE import_ts_initialize

    USE messy_main_mpi_bi,       ONLY: p_parallel_io, p_io, p_bcast
    USE messy_main_tools,        ONLY: find_next_free_unit

    IMPLICIT NONE
    INTRINSIC :: TRIM, ADJUSTL

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'import_ts_initialize'
    INTEGER                     :: iou    ! I/O unit
    INTEGER                     :: status ! error status
    INTEGER                     :: i, j

    CALL start_message_bi(submodstr, 'INITIALIZE',substr)

    ! INITIALIZE CTRL
    IF (p_parallel_io) THEN
       iou = find_next_free_unit(100,200)
       CALL import_ts_read_nml_ctrl(status, iou)
       IF (status /= 0) CALL error_bi(' ',substr)
    END IF

    IF (p_parallel_io) THEN
       NTS = 1
       DO i=1, NMAXTS

          IF (TRIM(TS(i)%name) == '') CYCLE

          CALL its_copy_io(XTS(NTS)%io, TS(i))

          WRITE(*,*) 'TIME SERIES   : ', TRIM(XTS(NTS)%io%name)
          WRITE(*,*) '   FILE       : ', TRIM(ADJUSTL(XTS(NTS)%io%fname))
          WRITE(*,*) '   VALID RANGE: ', XTS(NTS)%io%vr(:)

          DO j=1,2
             SELECT CASE(XTS(NTS)%io%cnt(j))
             CASE(TS_BD_STOP)
                WRITE(*,*) '   BND POLICY : stop'
             CASE(TS_BD_CONT)
                WRITE(*,*) '   BND POLICY : continue'
             CASE DEFAULT
                CALL error_bi('UNKNOWN BOUNDARY POLICY',substr)
             END SELECT
          END DO

          SELECT CASE(XTS(NTS)%io%im)
          CASE(TS_IM_PREV)
             WRITE(*,*) '   SELECTION  : previous'
          CASE(TS_IM_NEXT)
             WRITE(*,*) '   SELECTION  : next'
          CASE(TS_IM_LINT)
             WRITE(*,*) '   SELECTION  : linear interpolation'
          CASE DEFAULT
             CALL error_bi('UNKNOWN INTERPOLATION METHOD',substr)
          END SELECT

          WRITE(*,*) '   READING DATA ...'
          CALL its_read_ts(status, XTS(NTS), .FALSE.)
          IF (status /= 0) &
               CALL error_bi('its_read_ts reported an error' ,substr)
          WRITE(*,*) '   ... DONE!'

          ! NEXT TIME SERIES
          NTS = NTS + 1
          WRITE(*,*) '------------------------------------------------------'
       END DO
       NTS = NTS - 1
    END IF
    CALL p_bcast(NTS, p_io)

    ! BROADCAST ALL RESULTS
    DO i=1, NTS
       CALL p_bcast(XTS(i)%io%name, p_io)
       CALL p_bcast(XTS(i)%io%fname, p_io)
       CALL p_bcast(XTS(i)%io%vr(:), p_io)
       CALL p_bcast(XTS(i)%io%cnt(:), p_io)
       CALL p_bcast(XTS(i)%io%im, p_io)

       CALL p_bcast(XTS(i)%io%pdt(:), p_io)
       CALL p_bcast(XTS(i)%io%offset, p_io)

       CALL p_bcast(XTS(i)%nt, p_io)
       CALL p_bcast(XTS(i)%np, p_io)
       IF (.NOT. p_parallel_io) THEN
          ALLOCATE(XTS(i)%jd(XTS(i)%nt))
          ALLOCATE(XTS(i)%par(XTS(i)%np))
          ALLOCATE(XTS(i)%data(XTS(i)%nt,XTS(i)%np))
       END IF
       CALL p_bcast(XTS(i)%jd(:), p_io)
       CALL p_bcast(XTS(i)%par(:), p_io)
       DO j=1, XTS(i)%np
          CALL p_bcast(XTS(i)%data(:,j), p_io)
       END DO
       CALL p_bcast(XTS(i)%dimname, p_io)
    END DO

    CALL end_message_bi(submodstr, 'INITIALIZE',substr)

  END SUBROUTINE import_ts_initialize
  ! ----------------------------------------------------------------------

  ! ----------------------------------------------------------------------
  SUBROUTINE import_ts_init_memory

    USE messy_main_channel_bi,     ONLY: DC_BC
    USE messy_main_blather_bi,     ONLY: info_bi, warning_bi
    USE messy_main_channel,        ONLY: new_channel, new_channel_object   &
                                       , new_attribute
    USE messy_main_channel_error_bi,     ONLY: channel_halt
    USE messy_main_channel_dimensions,   ONLY: get_dimension_info       &
                                             , new_dimension            &
                                             , add_dimension_variable   &
                                             , add_dimension_variable_att &
                                             , DIMID_UNDEF
    USE messy_main_channel_repr,   ONLY: new_representation, AUTO &
                                       , set_representation_decomp &
                                       , IRANK, PIOTYPE_COL, REPR_UNDEF
    USE messy_main_channel_mem,    ONLY: dom_current
    USE messy_main_grid_netcdf,    ONLY: NF90_CHAR, string

    IMPLICIT NONE

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'import_ts_init_memory'
    INTEGER :: i, l, r, j, nr
    INTEGER :: status
    !
    INTEGER :: DIMID
    !
    INTEGER :: REPRID
    INTEGER :: nseg
    INTEGER, DIMENSION(:,:), POINTER :: start => NULL()
    INTEGER, DIMENSION(:,:), POINTER :: cnt   => NULL()
    INTEGER, DIMENSION(:,:), POINTER :: meml  => NULL()
    INTEGER, DIMENSION(:,:), POINTER :: memu  => NULL()
    !
    CHARACTER(LEN=10) :: str = ''
    CHARACTER(LEN=4)  :: nrstr = ''
    !
    INTEGER :: len
    LOGICAL :: lnewrepr = .TRUE.
    LOGICAL :: ldim_ok

    ALLOCATE(list_dimid(NTS))
    list_dimid(:) = DIMID_UNDEF
    ALLOCATE(list_reprid(NTS))
    list_reprid(:) = REPR_UNDEF

    CALL start_message_bi(submodstr, 'CHANNEL DEFINITION', substr)

    CALL new_channel(status, submodstr, dom_id=dom_current)
    CALL channel_halt(substr, status)

    nseg = 1
    ALLOCATE(start(nseg,IRANK))
    ALLOCATE(cnt(nseg,IRANK))
    ALLOCATE(meml(nseg,IRANK))
    ALLOCATE(memu(nseg,IRANK))

    start(:,:) = 1
    cnt(:,:) = 1
    meml(:,:) = 1
    memu(:,:) = 1

    DO i=1, NTS

       ! DEFINE DIMENSION
       CALL get_dimension_info(status, name=TRIM(XTS(i)%dimname) &
            , id=DIMID, len=len)

       dim_exists: IF (status == 0) THEN
          ! DIMENSION EXISTS ALREADY
          CALL info_bi('DIMENSION '//TRIM(XTS(i)%dimname)//' exists already' &
               , substr)
          length: IF (len /= XTS(i)%np) THEN
             ! APPEND ARBITRARY NUMBER TO AVOID ERRORS
             CALL warning_bi('DIMENSION '//TRIM(XTS(i)%dimname)//&
                  &' exists with wrong length', substr)
             ldim_ok = .FALSE.
             DO nr=1, 1000
                WRITE(nrstr,'(i4.4)') nr
                CALL get_dimension_info(status &
                     , name=TRIM(XTS(i)%dimname)//nrstr &
                     , id=DIMID, len=len)
                IF (status /= 0) THEN
                   ! dimension does not yet exist
                   XTS(i)%dimname = TRIM(XTS(i)%dimname)//nrstr
                   ldim_ok = .TRUE.
                   EXIT
                ELSE
                   IF (len == XTS(i)%np) THEN
                      ldim_ok = .TRUE.
                      XTS(i)%dimname = TRIM(XTS(i)%dimname)//nrstr
                   ENDIF
                END IF
             END DO
             IF (.NOT. ldim_ok) THEN
                CALL error_bi('DIMENSION '//TRIM(XTS(i)%dimname)//nrstr//&
                     &' could not be created', substr)
             ENDIF
          ELSE
             ldim_ok = .FALSE. ! exists already with correct length
          END IF length
       ELSE
          ! DIMENSION DOES NOT YET EXIST
          ldim_ok = .TRUE.
       END IF dim_exists

       new_dim_ok: IF (ldim_ok) THEN
          ! NEW DIMENSION
          CALL new_dimension(status, DIMID &
               , TRIM(XTS(i)%dimname), XTS(i)%np)
          CALL channel_halt(substr, status)

          ! DEFINE DIMENSION VARIABLE
          CALL add_dimension_variable(status &
               , TRIM(XTS(i)%dimname), TRIM(XTS(i)%dimname) &
               , XTS(i)%par(:) )
          CALL channel_halt(substr, status)

          ! ADD DIMENSION VARIABLE ATTRIBUTES
          IF (ASSOCIATED(XTS(i)%dimvaratt)) THEN
             DO j=1, SIZE(XTS(i)%dimvaratt)
                IF (XTS(i)%dimvaratt(j)%xtype == NF90_CHAR) THEN
                   CALL add_dimension_variable_att(status &
                        , TRIM(XTS(i)%dimname), TRIM(XTS(i)%dimname) &
                        , TRIM(XTS(i)%dimvaratt(j)%name)             &
                        , c = string(XTS(i)%dimvaratt(j)%dat%vc) )
                   CALL channel_halt(substr, status)
                ELSE
                   CALL info_bi('type of dimension variable attribute ('&
                        &//TRIM(XTS(i)%dimname)//' - '&
                        &//TRIM(XTS(i)%dimvaratt(j)%name)//&
                        &') currently not supported', substr)
                END IF
             END DO
          END IF
       END IF new_dim_ok

       ! SAVE DIMENSION ID
       list_dimid(i) = DIMID

       ! CHECK FOR REPRESENTATION
       lnewrepr = .TRUE.
       DO j=1, i-1
          IF (list_dimid(j) == DIMID) THEN
             REPRID = list_reprid(j)
             lnewrepr = .FALSE.
             EXIT
          END IF
       END DO

       IF (lnewrepr) THEN
          ! DEFINE REPRESENTATION
          CALL new_representation(status, REPRID, 'TSR_'//TRIM(XTS(i)%dimname) &
               , rank = 1, link = 'x---', dctype = DC_BC               &
               , dimension_ids = (/ DIMID /) &
               , ldimlen       = (/ AUTO  /) &
               , axis = 'N---'               &
               , dom_id=dom_current      &
               )
          CALL channel_halt(substr, status)

          start(:,1) = 1
          cnt(:,1)   = XTS(i)%np
          meml(:,1)  = 1
          memu(:,1)  = XTS(i)%np

          CALL set_representation_decomp(status, REPRID &
               , start, cnt, memu, meml, .FALSE., PIOTYPE_COL)
       END IF

       ! SAVE REPRESENTATION ID
       list_reprid(i) = REPRID

       ! DEFINE CHANNEL OBJECTS
       ! ... object with values
       CALL new_channel_object(status &
            , submodstr, TRIM(XTS(i)%io%name), p1=XTS(i)%obj, reprid=REPRID &
            , dom_id = dom_current)
       CALL channel_halt(substr, status)

       ! ADD VARIABLE ATTRIBUTES
       IF (ASSOCIATED(XTS(i)%varatt)) THEN
          DO j=1, SIZE(XTS(i)%varatt)
             IF (XTS(i)%varatt(j)%xtype == NF90_CHAR) THEN
                CALL new_attribute(status              &
                     , submodstr, TRIM(XTS(i)%io%name) &
                     , TRIM(XTS(i)%varatt(j)%name)     &
                     , c = string(XTS(i)%varatt(j)%dat%vc) )
                CALL channel_halt(substr, status)
             ELSE
                CALL info_bi('type of variable attribute ('&
                     &//TRIM(XTS(i)%io%name)//' - '&
                     &//TRIM(XTS(i)%varatt(j)%name)//&
                     &') currently not supported', substr)
             END IF
          END DO
       END IF

       CALL new_attribute(status, submodstr, TRIM(XTS(i)%io%name) &
            , 'from_file', c=TRIM(ADJUSTL(XTS(i)%io%fname)) )
       CALL channel_halt(substr, status)

       CALL new_attribute(status, submodstr, TRIM(XTS(i)%io%name) &
            , 'valid_range_l', r=XTS(i)%io%vr(1) )
       CALL channel_halt(substr, status)
       CALL new_attribute(status, submodstr, TRIM(XTS(i)%io%name) &
            , 'valid_range_u', r=XTS(i)%io%vr(2) )
       CALL channel_halt(substr, status)

       DO j=1,2
          l = (j-1)*5 + 1
          r = l + 4
          SELECT CASE(XTS(NTS)%io%cnt(j))
          CASE(TS_BD_STOP)
             str(l:r) = 'STOP '
          CASE(TS_BD_CONT)
             str(l:r) = 'CONT '
          END SELECT
       END DO

       CALL new_attribute(status, submodstr, TRIM(XTS(i)%io%name) &
            , 'boundary_policy', c=str)
       CALL channel_halt(substr, status)

       SELECT CASE(XTS(NTS)%io%im)
       CASE(TS_IM_PREV)
          str = 'PREV'
       CASE(TS_IM_NEXT)
          str = 'NEXT'
       CASE(TS_IM_LINT)
          str = 'LINT'
       END SELECT

       CALL new_attribute(status, submodstr, TRIM(XTS(i)%io%name) &
            , 'interpolation_method', c=str)
       CALL channel_halt(substr, status)

       ! ... object with flag
       CALL new_channel_object(status, &
            submodstr, TRIM(XTS(i)%io%name)//'_flg', p1=XTS(i)%flg &
            , reprid=REPRID &
            , dom_id=dom_current)
       CALL channel_halt(substr, status)

       CALL new_attribute(status, submodstr, TRIM(XTS(i)%io%name)//'_flg' &
            , 'flag_of_variable', c=TRIM(XTS(i)%io%name))
       CALL channel_halt(substr, status)

    END DO

    DEALLOCATE(start)
    DEALLOCATE(cnt)
    DEALLOCATE(meml)
    DEALLOCATE(memu)

    CALL end_message_bi(submodstr, 'CHANNEL DEFINITION', substr)

  END SUBROUTINE import_ts_init_memory
  ! ----------------------------------------------------------------------

  ! ----------------------------------------------------------------------
  SUBROUTINE import_ts_global_start

    USE messy_main_timer,            ONLY: YEAR, MONTH, DAY &
                                         , HOUR, MINUTE, SECOND
    USE messy_main_channel_error_bi, ONLY: channel_halt
    USE messy_main_channel,          ONLY: GET_CHANNEL_OBJECT
    USE messy_main_channel_mem,      ONLY: dom_current

    IMPLICIT NONE

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'import_ts_global_start'
    INTEGER :: i
    INTEGER :: status

    DO i=1, NTS
       CALL GET_CHANNEL_OBJECT(status, submodstr, TRIM(XTS(i)%io%name) &
            , p1=XTS(i)%obj, dom_id=dom_current)
       CALL channel_halt(substr, status)
       CALL GET_CHANNEL_OBJECT(status, submodstr, TRIM(XTS(i)%io%name)//'_flg' &
            , p1=XTS(i)%flg, dom_id=dom_current)
       CALL channel_halt(substr, status)

       CALL its_set_value_ts(status, XTS(i) &
            , YEAR, MONTH, DAY, HOUR, MINUTE, SECOND)
       IF (status /= 0) &
            CALL error_bi('its_set_value_ts reported an error '&
            &'for time series'//TRIM(XTS(i)%io%name),substr)
    END DO

  END SUBROUTINE import_ts_global_start
  ! ----------------------------------------------------------------------

  ! ----------------------------------------------------------------------
  SUBROUTINE import_ts_free_memory

    IMPLICIT NONE

    INTEGER :: i

    DO i=1, NTS
       CALL its_delete_ts(XTS(i))
    END DO

    IF (ASSOCIATED(list_dimid)) DEALLOCATE(list_dimid)
    NULLIFY(list_dimid)
    IF (ASSOCIATED(list_reprid)) DEALLOCATE(list_reprid)
    NULLIFY(list_reprid)

  END SUBROUTINE import_ts_free_memory
  ! ----------------------------------------------------------------------

! ============================================================================
END MODULE messy_main_import_ts_bi
! ============================================================================

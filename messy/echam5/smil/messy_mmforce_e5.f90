! **********************************************************************
MODULE messy_mmforce_e5
! **********************************************************************

  USE messy_mmforce

  IMPLICIT NONE
  PRIVATE

  INTRINSIC :: NULL

!!$  ! cpl namelist switches
!!$  LOGICAL :: l_diag_out_force = .FALSE.
!!$  LOGICAL :: l_diag_out_hw = .FALSE.

  ! WORKSPACE
  ! pointer to channel objects
  REAL(DP), DIMENSION (:,:,:), POINTER :: hw    => NULL()
  REAL(DP), DIMENSION (:,:,:), POINTER :: forc  => NULL()
  REAL(DP), DIMENSION (:,:,:), POINTER :: press  => NULL()

  ! PUBLIC INTERFACE ROUTINES
  PUBLIC :: mmforce_initialize
  PUBLIC :: mmforce_init_memory
  PUBLIC :: mmforce_init_coupling
  PUBLIC :: mmforce_global_start
  PUBLIC :: mmforce_physc
!!$  !PRIVATE :: mmforce_read_nml_cpl

CONTAINS

! ------------------------------------------------------------------------
  SUBROUTINE mmforce_initialize

    ! ECHAM5/MESSy
    USE messy_main_mpi_bi,   ONLY: p_parallel_io, p_io, p_bcast, finish
    USE messy_main_tools,    ONLY: find_next_free_unit

    IMPLICIT NONE

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'mmforce_initialize'
    INTEGER         :: iou    ! I/O unit
    INTEGER         :: status ! error status

    ! INITIALIZE CTRL
    IF (p_parallel_io) THEN
       iou = find_next_free_unit(100,200)
       ! *** CALL SMCL ROUTINE:
       CALL mmforce_read_nml_ctrl(status, iou)
       IF (status /= 0) CALL finish(substr)
    END IF
    CALL p_bcast(i_month, p_io)
    CALL p_bcast(r_dudt, p_io)

!!$    ! INITIALIZE CPL
!!$    IF (p_parallel_io) THEN
!!$       iou = find_next_free_unit(100,200)
!!$       ! *** CALL SMIL ROUTINE:
!!$       CALL mmforce_read_nml_cpl(status, iou)
!!$       IF (status /= 0) CALL finish(substr)
!!$    END IF
!!$    CALL p_bcast(l_diag_out_force, p_io)
!!$    CALL p_bcast(l_diag_out_hw, p_io)

  END SUBROUTINE mmforce_initialize
! ------------------------------------------------------------------------

! ------------------------------------------------------------------------
  SUBROUTINE mmforce_init_memory

    ! ECHAM5/MESSy
    USE messy_main_channel_error_bi, ONLY: channel_halt
    USE messy_main_channel_bi,       ONLY: GP_3D_MID
    ! MESSy
    USE messy_main_channel,    ONLY: new_channel, new_channel_object &
                                   , new_attribute

    IMPLICIT NONE

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'mmforce_init_memory'
    INTEGER                     :: status

    ! DEFINE NEW CHANNEL
    CALL new_channel(status, modstr)
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr, 'forc' &
         , p3 = forc, reprid = GP_3D_MID )
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'forc'   &
         , 'long_name', c='momentum forcing' )
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'forc'   &
         , 'units', c='m/s^2' )
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr, 'hw' &
         , p3 = hw, reprid = GP_3D_MID )
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'hw'   &
         , 'long_name', c='weighting factor' )
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'hw'   &
         , 'units', c='1' )
    CALL channel_halt(substr, status)

  END SUBROUTINE mmforce_init_memory
! ------------------------------------------------------------------------

! ------------------------------------------------------------------------
  SUBROUTINE mmforce_init_coupling

    !  ECHAM5/MESSy
    USE messy_main_channel_error_bi,  ONLY: channel_halt
    ! MESSy
    USE messy_main_channel,           ONLY: get_channel_object

    IMPLICIT NONE

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER       :: substr = 'mmforce_init_coupling'
    INTEGER                           :: status

    CALL get_channel_object(status, 'ECHAM5', 'press', p3=press)
    CALL channel_halt(substr, status)
    
  END SUBROUTINE mmforce_init_coupling
! ------------------------------------------------------------------------

! ------------------------------------------------------------------------
  SUBROUTINE mmforce_global_start

    USE messy_main_timer,     ONLY: YearDay, current_date
    USE messy_mmforce,        ONLY: dayno

    IMPLICIT NONE

    dayno = YearDay(current_date)

  END SUBROUTINE mmforce_global_start
! ------------------------------------------------------------------------

! ------------------------------------------------------------------------
  SUBROUTINE mmforce_physc

    USE messy_main_data_bi,         ONLY: vom_3d
    USE messy_main_grid_def_bi,     ONLY: philat_2d
    USE messy_main_grid_def_mem_bi, ONLY: jrow

    IMPLICIT NONE

    CALL mmforce_physc_int(vom_3d(:,:,jrow), press(:,:,jrow) &
         , philat_2d(:,jrow), forc(:,:,jrow), hw(:,:,jrow))

  END SUBROUTINE mmforce_physc
! ------------------------------------------------------------------------

! ------
! PRIVATE
! ------------------------------------------------------------------------
!!$  SUBROUTINE mmforce_read_nml_cpl(status, iou)
!!$
!!$    USE messy_main_tools, ONLY: read_nml_open, read_nml_check, read_nml_close
!!$
!!$    IMPLICIT NONE
!!$
!!$    ! I/O
!!$    INTEGER, INTENT(IN)  :: iou    ! logical I/O unit
!!$    INTEGER, INTENT(OUT) :: status ! error status
!!$
!!$    ! LOCAL
!!$    CHARACTER(LEN=*), PARAMETER  :: substr = 'mmforce_read_nml_cpl'
!!$    LOGICAL                      :: lex          ! file exists?
!!$    INTEGER                      :: fstat        ! file status
!!$
!!$    NAMELIST /CPL/ l_diag_out_force, l_diag_out_hw
!!$
!!$    ! initialize
!!$    status = 1 ! error
!!$
!!$    CALL read_nml_open(lex, substr, iou, 'CPL', modstr)
!!$    IF (.NOT.lex) RETURN    ! <modstr>.nml does not exist
!!$
!!$    READ(iou, NML = CPL, IOSTAT = fstat)
!!$    CALL read_nml_check(fstat, substr, iou, 'CPL', modstr)
!!$    IF (fstat /= 0) RETURN  ! error while reading namelist
!!$
!!$    CALL read_nml_close(substr, iou, modstr)
!!$
!!$    status = 0 ! no error    
!!$
!!$  END SUBROUTINE mmforce_read_nml_cpl
! ------------------------------------------------------------------------

! **********************************************************************
END MODULE messy_mmforce_e5
! **********************************************************************

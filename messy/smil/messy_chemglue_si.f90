!*****************************************************************************
!                Time-stamp: <2020-11-27 13:48:33 b302010>
!*****************************************************************************

! SUBMODEL INTERFACE LAYER (SMIL) ROUTINES FOR MESSy SUBMODEL CHEMGLUE
! Author: Rolf Sander, 2014-...

!*****************************************************************************

#include "messy_main_ppd_bi.inc"

MODULE messy_chemglue_si

  USE messy_main_constants_mem, ONLY: DP
  USE messy_main_tracer_mem_bi, ONLY: ti_gp
  USE messy_main_blather_bi, ONLY: start_message_bi, end_message_bi, &
    error_bi, info_bi, warning_bi
  USE messy_chemglue ! mecnum etc.

  IMPLICIT NONE

  ! CPL-NAMELIST PARAMETERS
  INTEGER :: cpl_nml_dummy = 0 ! only dummy, currently not used

  REAL(DP), DIMENSION(:,:,:), POINTER, SAVE :: meccanum_gp
  INTEGER, SAVE :: zidt_C5H8, zidt_APINENE, zidt_TOLUENE ! tracer indices (local)

  PRIVATE
  PUBLIC :: chemglue_initialize    ! initialize submodel
  PUBLIC :: chemglue_init_memory   ! request memory
  PUBLIC :: chemglue_physc         ! entry point in time loop (current vector)

CONTAINS

  !***************************************************************************

  SUBROUTINE chemglue_initialize

    USE messy_main_mpi_bi,    ONLY: p_parallel_io, p_io, p_bcast
    USE messy_main_tools,     ONLY: find_next_free_unit

    IMPLICIT NONE

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'chemglue_initialize'
    INTEGER                     :: status ! error status
    INTEGER                     :: iou    ! I/O unit

    CALL start_message_bi(modstr,'INITIALISATION',substr)

    ! READ CTRL namelist:
    IF (p_parallel_io) THEN
      iou = find_next_free_unit(100,200)
      CALL chemglue_read_nml_ctrl(status, iou)
      IF (status /= 0) CALL error_bi('Error in reading CTRL namelist',substr)
    ENDIF
    CALL p_bcast(ctrl_nml_dummy, p_io)

    ! READ CPL namelist:
    IF (p_parallel_io) THEN
      iou = find_next_free_unit(100,200)
      CALL chemglue_read_nml_cpl(status, iou)
      IF (status /= 0) CALL error_bi('Error in reading CPL namelist',substr)
    ENDIF
    CALL p_bcast(cpl_nml_dummy, p_io)

    ! assign names to the numbers of the MECCA mechanisms:
    CALL assign_mecnum_names

    IF (NMAXMECCA == 1) THEN
      CALL error_bi('CHEMGLUE can only be used if more than one chemistry '// &
        'mechanism has been produced with xpolymecca.', substr)
    ENDIF

    CALL end_message_bi(modstr,'INITIALISATION',substr)

  END SUBROUTINE chemglue_initialize

  ! --------------------------------------------------------------------------

  SUBROUTINE chemglue_init_memory

    USE messy_main_channel_error_bi, ONLY: channel_halt
    USE messy_main_channel_bi,       ONLY: GP_3D_MID
    USE messy_main_channel,          ONLY: new_channel, new_channel_object, &
                                           new_attribute
    USE messy_main_tracer,           ONLY: get_tracer
    USE messy_main_tracer_mem_bi,    ONLY: GPTRSTR

    IMPLICIT NONE
    CHARACTER(LEN=*), PARAMETER :: substr = 'chemglue_init_memory'
    INTEGER                     :: status

    ! create channel and channel object:
    CALL start_message_bi(modstr,'CHANNEL DEFINITION',substr)
    CALL new_channel(status, modstr//'_gp', reprid=GP_3D_MID)
    CALL channel_halt(substr, status)
    CALL new_channel_object(status, modstr//'_gp', 'meccanum', p3=meccanum_gp)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr//'_gp', 'meccanum', &
      'long_name', c='MECCA number')
    CALL channel_halt(substr, status)
    CALL info_bi('channel/object '//modstr//'_gp/meccanum was created')
    CALL end_message_bi(modstr,'CHANNEL DEFINITION',substr)

    ! define tracer indices for chemglue_physc:
    CALL get_tracer(status, GPTRSTR, 'C5H8',    idx=zidt_C5H8)
    IF (status/=0) zidt_C5H8    = 0
    CALL get_tracer(status, GPTRSTR, 'APINENE', idx=zidt_APINENE)
    IF (status/=0) zidt_APINENE = 0
    CALL get_tracer(status, GPTRSTR, 'TOLUENE', idx=zidt_TOLUENE)
    IF (status/=0) zidt_TOLUENE = 0

  END SUBROUTINE chemglue_init_memory

  ! --------------------------------------------------------------------------

  SUBROUTINE chemglue_physc

    USE messy_main_grid_def_mem_bi, ONLY: kproma, nlev, jrow
    USE messy_main_data_bi,         ONLY: press_3d
    USE messy_main_tracer_mem_bi,   ONLY: pxtte => qxtte, pxtm1 => qxtm1
    USE messy_main_timer,           ONLY: time_step_len

    IMPLICIT NONE
    REAL(DP) :: y_C5H8, y_APINENE, y_TOLUENE ! mixing ratios
    CHARACTER(LEN=*), PARAMETER :: substr = 'chemglue_physc'
    INTEGER :: jk, jp
    
    meccanum_gp(_RI_XYZ__(1:kproma,jrow,:)) = 0._DP
    DO jk=1,nlev
      DO jp=1,kproma
        ! --------------------------------------------------------------------
        ! 1) pressure: _RI_XYZ__(jp,jrow,jk)
        ! CALL select_mechanism_from_pressure( &
        !        meccanum_gp(_RI_XYZ__(jp,jrow,jk)), &
        !        press_3d(_RI_XYZ__(jp,jrow,jk))
        ! --------------------------------------------------------------------
        ! 2) concentration:
        IF (zidt_C5H8>0) THEN
          y_C5H8 = pxtm1(jp,jk,zidt_C5H8) + pxtte(jp,jk,zidt_C5H8) * time_step_len
        ELSE
          y_C5H8 = 0.
        ENDIF
        IF (zidt_APINENE>0) THEN
          y_APINENE = pxtm1(jp,jk,zidt_APINENE) + pxtte(jp,jk,zidt_APINENE) * time_step_len
        ELSE
          y_APINENE = 0.
        ENDIF
        IF (zidt_TOLUENE>0) THEN
          y_TOLUENE = pxtm1(jp,jk,zidt_TOLUENE) + pxtte(jp,jk,zidt_TOLUENE) * time_step_len
        ELSE
          y_TOLUENE = 0.
        ENDIF
        CALL select_mechanism_from_mixrat(meccanum_gp(_RI_XYZ__(jp,jrow,jk)), &
          y_C5H8, y_APINENE, y_TOLUENE)
        ! --------------------------------------------------------------------
        ! 3) sea-land mask:
        ! CALL select_mechanism_from_slm(meccanum_gp, slm)
        ! --------------------------------------------------------------------
        ! 4) testing only:
        ! CALL select_mechanism_testing(meccanum_gp(_RI_XYZ__(jp,jrow,jk)), jp)
        ! --------------------------------------------------------------------
      ENDDO
    ENDDO
    IF (MINVAL(meccanum_gp(_RI_XYZ__(1:kproma,jrow,:))).LE.0._DP) THEN
      CALL error_bi('Error: mechanism number not assigned',substr)
    ENDIF

  END SUBROUTINE chemglue_physc

  ! --------------------------------------------------------------------------

  SUBROUTINE chemglue_read_nml_cpl(status, iou)

    USE messy_main_tools,  ONLY: read_nml_open, read_nml_check, read_nml_close
    IMPLICIT NONE
    ! I/O
    INTEGER, INTENT(OUT) :: status     ! error status
    INTEGER, INTENT(IN)  :: iou        ! I/O unit
    NAMELIST /CPL/ cpl_nml_dummy
    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr='chemglue_read_nml_cpl'
    LOGICAL                     :: lex      ! file exists ?
    INTEGER                     :: fstat    ! file status

    status = 1
    CALL read_nml_open(lex, substr, iou, 'CPL', modstr)
    IF (.not.lex) RETURN    ! <modstr>.nml does not exist
    READ(iou, NML=CPL, IOSTAT=fstat)
    CALL read_nml_check(fstat, substr, iou, 'CPL', modstr)
    IF (fstat /= 0) RETURN  ! error while reading namelist
    CALL read_nml_close(substr, iou, modstr)
    status = 0 ! NO ERROR

  END SUBROUTINE chemglue_read_nml_cpl

  !***************************************************************************

END MODULE messy_chemglue_si

!*****************************************************************************

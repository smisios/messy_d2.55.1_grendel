#include "messy_main_ppd_bi.inc"

! **********************************************************************+
MODULE messy_main_tracer_tools_bi
! **********************************************************************+

  USE messy_main_tracer_mem_bi
!!$  USE messy_main_blather_bi,    ONLY: start_message_bi, end_message_bi
!!$  USE messy_main_constants_mem, ONLY: FLAGGED_BAD
!!$#ifdef ICON
!!$  USE messy_main_tools,         ONLY: PTR_3D_ARRAY
!!$#endif
  USE messy_main_tracer

  IMPLICIT NONE
  PRIVATE

  ! SUBMODEL SMIL
  PUBLIC :: tracer_halt
  !
#if defined(ICON)
  PUBLIC :: main_tracer_set_domain
#endif
  !
CONTAINS

  ! -------------------------------------------------------------------
  SUBROUTINE tracer_halt(substr, status)

    USE messy_main_blather_bi,    ONLY: error_bi
    USE messy_main_constants_mem, ONLY: STRLEN_VLONG

    IMPLICIT NONE

    ! I/O
    CHARACTER(LEN=*), INTENT(IN)  :: substr
    INTEGER,          INTENT(IN)  :: status
    ! LOCAL
    CHARACTER(LEN=STRLEN_VLONG)   :: errstr

    IF (status /= 0) THEN
       errstr = tracer_error_str(status)
       CALL error_bi(errstr, substr)
    END IF

  END SUBROUTINE tracer_halt
  ! -------------------------------------------------------------------

  ! -------------------------------------------------------------------
#if defined(ICON)
  SUBROUTINE main_tracer_set_domain(dom_id, jrow)

    IMPLICIT NONE

    ! I/O
    INTEGER, INTENT(IN)           :: dom_id
    INTEGER, INTENT(IN), OPTIONAL :: jrow

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'main_tracer_set_domain'
    INTEGER :: status

    GPTRSTR    => L_GPTRSTR(dom_id)
    !gp_channel => l_gp_channel(dom_id)
    ntrac_icon => l_ntrac_icon(dom_id)

    CALL get_tracer_set(status, GPTRSTR    & 
         , trlist_gp, ti_gp, ntrac_gp      &
         , xt=pxt, xtte=pxtte, xtm1=pxtm1)
    CALL tracer_halt(substr, status)

    xt    => pxt(:,:,:,:,1)
    xtm1  => pxtm1(:,:,:,:,1)
    xtte  => pxtte(:,:,:,:,1)

    IF (PRESENT(jrow)) THEN
       qxt    => xt(_RI_XYZN_(:,jrow,:,:))
       qxtte  => xtte(_RI_XYZN_(:,jrow,:,:))
       qxtm1  => xtm1(_RI_XYZN_(:,jrow,:,:))
    END IF

  END SUBROUTINE main_tracer_set_domain
#endif
  ! -------------------------------------------------------------------

! **********************************************************************+
END MODULE messy_main_tracer_tools_bi
! **********************************************************************+

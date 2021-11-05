! **********************************************************************+
MODULE messy_main_mpi_bi
! **********************************************************************+

  ! THIS MODULE PROVIDES DUMMY PARAMETERS AND SUBROUTINES FOR
  ! REPLACING A PARALLEL ENVIRONMENT

  USE messy_main_constants_mem, ONLY: dp

  IMPLICIT NONE
  PUBLIC

  LOGICAL, PARAMETER :: p_parallel_io = .TRUE.  ! only one PE
  LOGICAL, PARAMETER :: p_parallel    = .FALSE. ! only one PE
  INTEGER, PARAMETER :: p_io          = 0       ! only one PE
  INTEGER, PARAMETER :: p_pe          = 0       ! only one PE
  INTEGER, PARAMETER :: p_nprocs      = 1       ! only one PE

  INTERFACE p_bcast
     MODULE PROCEDURE p_bcast_c
     MODULE PROCEDURE p_bcast_i
     MODULE PROCEDURE p_bcast_r
     MODULE PROCEDURE p_bcast_r_2d
     MODULE PROCEDURE p_bcast_l
     MODULE PROCEDURE p_bcast_l_1d
  END INTERFACE

  ! op_pj_20120925+
  INTERFACE p_sum
     MODULE PROCEDURE p_sum_2d
  END INTERFACE  
  ! op_pj_20120925-

CONTAINS

  ! ====================================================================
  SUBROUTINE finish(caller, message)

    IMPLICIT NONE
    INTRINSIC :: PRESENT, TRIM

    ! I/O
    CHARACTER(LEN=*), INTENT(IN)           :: caller
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: message

    IF (PRESENT(message)) THEN
       WRITE(*,*) 'STOPPED BY '//TRIM(caller)//': '//TRIM(message)
    ELSE
       WRITE(*,*) 'STOPPED BY '//TRIM(caller)
    END IF

    STOP

  END SUBROUTINE finish
  ! ====================================================================

  ! ====================================================================
  SUBROUTINE p_bcast_c(c, p)

    ! I/O
    CHARACTER(LEN=*), INTENT(IN) :: c
    INTEGER,          INTENT(IN) :: p

  END SUBROUTINE p_bcast_c
  ! ====================================================================

  ! ====================================================================
  SUBROUTINE p_bcast_i(i, p)

    ! I/O
    INTEGER, INTENT(IN) :: i
    INTEGER, INTENT(IN) :: p

  END SUBROUTINE p_bcast_i
  ! ====================================================================

  ! ====================================================================
  SUBROUTINE p_bcast_r(r, p)

    ! I/O
    REAL(dp), INTENT(IN) :: r
    INTEGER,  INTENT(IN) :: p

  END SUBROUTINE p_bcast_r
  ! ====================================================================

  ! ====================================================================
  SUBROUTINE p_bcast_r_2d(r, p)

    ! I/O
    REAL(dp), DIMENSION(:,:), INTENT(IN) :: r
    INTEGER,                  INTENT(IN) :: p

  END SUBROUTINE p_bcast_r_2d
  ! ====================================================================

  ! ====================================================================
  SUBROUTINE p_bcast_l(l, p)

    ! I/O
    LOGICAL, INTENT(IN) :: l
    INTEGER, INTENT(IN) :: p

  END SUBROUTINE p_bcast_l
  ! ====================================================================

  ! ====================================================================
  SUBROUTINE p_bcast_l_1d(l, p)

    ! I/O
    LOGICAL, DIMENSION(:), INTENT(IN) :: l
    INTEGER,               INTENT(IN) :: p

  END SUBROUTINE p_bcast_l_1d
  ! ====================================================================

  ! ====================================================================
  SUBROUTINE p_abort

  END SUBROUTINE p_abort
  ! ====================================================================

  ! op_pj_20120925+
  ! ====================================================================
  FUNCTION p_sum_2d (zfield, comm) RESULT (p_sum)
    
    REAL(dp),          INTENT(in) :: zfield(:,:)
    INTEGER, OPTIONAL, INTENT(in) :: comm
    REAL(dp)                      :: p_sum (SIZE(zfield,1),SIZE(zfield,2))

    p_sum = zfield

  END FUNCTION p_sum_2d
  ! ====================================================================
  ! op_pj_20120925-

! **********************************************************************+
END MODULE messy_main_mpi_bi
! **********************************************************************+

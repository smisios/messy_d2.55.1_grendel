! **************************************************************************
MODULE messy_submodel
! **************************************************************************

  ! MODULE FOR ILLUSTRATING THE USE OF CHANNEL / OBJECTS
  !
  ! Author: Patrick Joeckel, MPICH, Jan 2009
  !
  ! ----------- >

  USE messy_main_constants_mem, ONLY: DP

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: DP

  ! ----------- <

  CHARACTER(LEN=*),PARAMETER, PUBLIC :: modstr = 'submodel'
  CHARACTER(LEN=*),PARAMETER, PUBLIC :: modver = '1.0'

  INTEGER, DIMENSION(:), ALLOCATABLE, SAVE :: ISTATE ! STATE VECTOR

  PUBLIC :: RANINI
  PUBLIC :: RANCLEAN
  PUBLIC :: populate_bnds

CONTAINS

  ! -------------------------------------------------------------------
  SUBROUTINE RANINI(STATE, LRESUME, HARVEST)
    
    IMPLICIT NONE
    
    INTRINSIC :: RANDOM_SEED, RANDOM_NUMBER, NINT
    
    ! I/O
    REAL(dp), DIMENSION(:), INTENT(INOUT)         :: STATE
    LOGICAL,                INTENT(IN)            :: LRESUME
    REAL(dp), DIMENSION(:), INTENT(OUT), OPTIONAL :: HARVEST
    
    ! LOCAL
    INTEGER :: NSTATE  ! LENGTH OF RANDOM STATE VECTOR
  
    IF (.NOT.PRESENT(HARVEST)) THEN
       ! INITIALISATION
       CALL RANDOM_SEED(SIZE = NSTATE)
       ALLOCATE(ISTATE(NSTATE))
       IF (LRESUME) THEN
          ! RESTART: RE-SET STATE VECTOR
          ISTATE(:) = NINT(STATE(:))
          CALL RANDOM_SEED(PUT = ISTATE)
       ELSE
          ! START: INITIALISE STATE VECTOR
          CALL RANDOM_SEED()
       END IF
    ELSE
       ISTATE(:) = NINT(STATE(:))
       CALL RANDOM_SEED(PUT = ISTATE)
       CALL RANDOM_NUMBER(HARVEST)
    END IF
    
    ! ALWAYS UPDATE STATE VECTOR FOR RESTART
    CALL RANDOM_SEED(GET = ISTATE)
    STATE(:) = REAL(ISTATE(:), DP)
    
  END SUBROUTINE RANINI
  ! -------------------------------------------------------------------

  ! -------------------------------------------------------------------
  SUBROUTINE RANCLEAN

    IMPLICIT NONE

    IF (ALLOCATED(ISTATE)) DEALLOCATE(ISTATE)

  END SUBROUTINE RANCLEAN
  ! -------------------------------------------------------------------

! mz_ab_20090921+
  ! --------------------------------------------------------------------
  SUBROUTINE populate_bnds(U, s, e)

    IMPLICIT NONE
    INTRINSIC :: REAL

    INTEGER, DIMENSION(:), INTENT(IN) :: s, e
    REAL(DP),DIMENSION(s(1):e(1),s(2):e(2),s(3):e(3)),INTENT(INOUT) :: U

    ! LOCAL
    INTEGER :: i, j, k
    
    write(*,*) modstr,' (populate_bnds): ',s,' -> ',e

    DO i=s(1),e(1)
       DO j=s(3),e(3)
          DO k=s(2),e(2)
             U(i,k,j) = REAL(i+j, DP)
          ENDDO
       ENDDO
    ENDDO
     
  END SUBROUTINE populate_bnds
  ! --------------------------------------------------------------------  
! mz_ab_20090921-

! **************************************************************************
END MODULE messy_submodel
! **************************************************************************

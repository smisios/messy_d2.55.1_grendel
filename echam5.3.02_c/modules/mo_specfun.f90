! Due to a misdesign in xlf 7.1 it is necessary to overload MERGE on
! IBMs. This modules provide the functionality of the intrinsic function 
! MERGE. 
!
! L. Kornblueh, MPI, January 2001

#ifdef __ibm__
#define _OVERLOAD_MERGE
#endif

MODULE mo_specfun

  USE mo_kind, ONLY: sp, dp

  PRIVATE

! for testing:
!#undef _OVERLOAD_MERGE 

#ifdef _OVERLOAD_MERGE
  PUBLIC::merge
 
  INTERFACE merge
     MODULE PROCEDURE emerge_r4, emerge_r8, emerge_i
  END INTERFACE

  LOGICAL, PARAMETER, PUBLIC :: have_intrinsic_merge = .FALSE.

#else

  LOGICAL, PARAMETER, PUBLIC :: have_intrinsic_merge = .TRUE.

#endif

CONTAINS

  ELEMENTAL FUNCTION emerge_r4(a, b, m) RESULT(c)
    REAL(sp), INTENT(in) :: a, b
    LOGICAL,  INTENT(in) :: m
    REAL(sp)             :: c
    IF (m) THEN
      c = a
    ELSE
      c = b
    ENDIF
  END FUNCTION emerge_r4

  ELEMENTAL FUNCTION emerge_r8(a, b, m) RESULT(c)
    REAL(dp), INTENT(in) :: a, b
    LOGICAL,  INTENT(in) :: m
    REAL(dp)             :: c
    IF (m) THEN
      c = a
    ELSE
      c = b
    ENDIF
  END FUNCTION emerge_r8

  ELEMENTAL FUNCTION emerge_i(a, b, m) RESULT(c)
    INTEGER, INTENT(in) :: a, b
    LOGICAL, INTENT(in) :: m
    INTEGER             :: c
    IF (m) THEN
      c = a
    ELSE
      c = b
    ENDIF
  END FUNCTION emerge_i

END MODULE mo_specfun













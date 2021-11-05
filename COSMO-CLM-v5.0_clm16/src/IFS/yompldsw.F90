MODULE YOMPLDSW

! Platform dependent switches

USE PARKIND1  ,ONLY : JPIM     ,JPRB

IMPLICIT NONE

SAVE

!#ifdef RS6K
!LOGICAL :: LOPT_SCALAR=.TRUE.
!LOGICAL :: LOPT_RS6K=.TRUE.
!#else
LOGICAL :: LOPT_SCALAR=.FALSE.
LOGICAL :: LOPT_RS6K=.FALSE.
!#endif
!     ------------------------------------------------------------------
END MODULE YOMPLDSW

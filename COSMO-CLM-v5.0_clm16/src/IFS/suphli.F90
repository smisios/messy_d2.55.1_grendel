!OPTIONS XOPT(HSFUN)
SUBROUTINE SUPHLI

!     ------------------------------------------------------------------

!**   *SUPHLI* IS THE SET-UP ROUTINE FOR COMMON BLOCK *YOEPHLI*

!     J.F. MAHFOUF         E.C.M.W.F.      96/06/23


!     PURPOSE
!     -------

!          THIS ROUTINE INITIALIZES THE CONSTANTS IN COMMON BLOCK
!     *YOEPHLI*

!     INTERFACE.
!     ----------

!     CALL *SUPHLI* FROM *SUPHEC*


!     METHOD.
!     -------

!         INITIALIZATION OF THE CONSTANTS USED IN THE LINEARIZED
!         PHYSICS

!     EXTERNALS.
!     ----------

!        NONE

!     REFERENCE.
!     ----------

!     MODIFICATIONS
!     -------------
!----------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB

USE YOETHF   , ONLY : RTWAT    ,RTICE
USE YOEPHLI  , ONLY : LPHYLIN   ,LENOPERT ,LRAISANEN,&
            &RLPTRC   ,RLPAL1   ,RLPAL2   ,RLPBB    ,&
            &RLPCC    ,RLPDD    ,RLPMIXL  ,RLPBETA  ,&
            &RLPDRAG  ,RLPEVAP  ,RLPP00
USE YOERAD   , ONLY : LRRTM     , NOVLP
!USE YOMLUN   , ONLY : NULNAM
!USE YOMTLEVOL, ONLY : LTLEVOL

IMPLICIT NONE

!#include "namtlevol.h"

!     ------------------------------------------------------------------

!*         1.     SET LOGICAL TO SWICH ON LINEARIZED PHYSICS
!                 ------------------------------------------

LPHYLIN = .FALSE.
!LTLEVOL = .FALSE.

!CALL POSNAM(NULNAM,'NAMTLEVOL')
!READ(NULNAM,NAMTLEVOL)

!*         1.1 No perturbation of surface arrays
!          -------------------------------------

LENOPERT = .TRUE.

!*         2.     SET CONSTANTS RELATED TO WATER MIXED PHASE
!                 ------------------------------------------


RLPTRC=RTICE+(RTWAT-RTICE)/SQRT(2._JPRB)
RLPAL1=0.15_JPRB
RLPAL2=20._JPRB


!*         3.     SET CONSTANTS RELATED TO VERTICAL DIFFUSION
!                 -------------------------------------------


!     CONSTANTS OF THE LOUIS FORMULATION 

RLPBB=5._JPRB
RLPCC=5._JPRB
RLPDD=5._JPRB

!     PSEUDO DEPTH OF THE BOUNDARY LAYER

RLPMIXL=4000._JPRB

!     REDUCTION FACTOR OF THE ASYMPTOTIC MIXING LENGTH

RLPBETA=0.2_JPRB


!*         4.     SET CONSTANTS RELATED TO GRAVITY WAVE DRAG
!                 ------------------------------------------

RLPDRAG=0._JPRB


!*         5.     SET CONSTANTS RELATED TO RAINFALL EVAPORATION 
!                 ---------------------------------------------

RLPEVAP=0._JPRB


!*         6.     SET CONSTANTS RELATED TO RADIATION
!                 ----------------------------------

!     Pressure level above which long-wave cooling is not applied

RLPP00=30000._JPRB

!     Using Raisanen overlap scheme

!*UPG Oper
!IF (LRRTM .AND. (NOVLP.EQ.1)) THEN
!  LRAISANEN = .TRUE.
!ELSE
  LRAISANEN = .FALSE.
!ENDIF

!     -------------------------------------------------------------------

RETURN
END SUBROUTINE SUPHLI

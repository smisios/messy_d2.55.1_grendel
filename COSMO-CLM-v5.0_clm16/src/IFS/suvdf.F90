SUBROUTINE SUVDF

!     ------------------------------------------------------------------

!**   *SUVDF* IS THE SET-UP ROUTINE FOR COMMON BLOCK *YOEVDF*

!     A.C.M. BELJAARS         E.C.M.W.F.       2/11/89

!     PURPOSE
!     -------

!          THIS ROUTINE INITIALIZES THE CONSTANTS IN COMMON BLOCK
!     *YOEVDF*

!     INTERFACE.
!     ----------

!     CALL *SUVDF* FROM *SUPHEC*


!     METHOD.
!     -------

!     EXTERNALS.
!     ----------

!     NONE.

!     REFERENCE.
!     ----------

!     MODIFICATIONS
!     -------------
!     J.-J. MORCRETTE         E.C.M.W.F.      91/07/14
!     ------------------------------------------------------------------


USE PARKIND1  ,ONLY : JPIM     ,JPRB

USE YOEVDF   , ONLY : RLAM     ,RKAP     ,RCHAR    ,RVDIFTS  ,&
            &RZ0ICE   ,REPDU2   ,REPUST   ,RSEZ0H   ,&
            &RSEZ0Q   ,RNUM     ,RNUH     ,RNUQ     ,RENTR    ,&
            &RPAR     ,RPAR1    ,RPARSRF  ,RPARZI   ,RLAMSK   ,&
            &LELWDD


IMPLICIT NONE


!     LOCAL REAL SCALARS
REAL(KIND=JPRB) :: CEPZ0O, ZNU


!     ------------------------------------------------------------------

!*         1.     SET FIRST SET OF CONSTANTS
!                 --------------------------


RLAM   =150._JPRB
RKAP   =0.4_JPRB
RCHAR  =0.018_JPRB
RVDIFTS=1.5_JPRB
RZ0ICE =0.001_JPRB

!     ------------------------------------------------------------------

!*         2.      SET OTHER CONSTANTS
!                  -------------------


CEPZ0O=2._JPRB
REPDU2 =(0.1_JPRB)**2
REPUST=0.0001_JPRB
RSEZ0H=1.4E-5_JPRB
RSEZ0Q=1.3E-4_JPRB

!     KINEMATIC VISCOSITY OF AIR

ZNU   =1.5E-5_JPRB
RNUM  =0.11_JPRB*ZNU
RNUH  =0.40_JPRB*ZNU
RNUQ  =0.62_JPRB*ZNU

!     ENTRAINMENT PARAMETRIZATION

RENTR=0.20_JPRB
RPAR=2._JPRB
RPAR1=0.6_JPRB
RPARSRF=0.1_JPRB
RPARZI=1000._JPRB

!     COMPUTATION OF SKIN TEMPERATURE

RLAMSK=15._JPRB
LELWDD=.TRUE.

!     ------------------------------------------------------------------

RETURN
END SUBROUTINE SUVDF

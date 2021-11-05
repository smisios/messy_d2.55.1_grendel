SUBROUTINE SUCLDP


!**** *SUCLDP*   - INITIALIZE COMMON YOECLD CONTROLLING *CLOUDSC*

!     PURPOSE.
!     --------
!           INITIALIZE YOECLDP

!**   INTERFACE.
!     ----------
!        CALL *SUCLDP* FROM *SUPHEC*
!              ------        ------

!        EXPLICIT ARGUMENTS :
!        --------------------
!        NONE

!        IMPLICIT ARGUMENTS :
!        --------------------
!        COMMON YOECLDP

!     METHOD.
!     -------
!        SEE DOCUMENTATION

!     EXTERNALS.
!     ----------
!        NONE

!     REFERENCE.
!     ----------
!        ECMWF RESEARCH DEPARTMENT DOCUMENTATION OF THE
!     "INTEGRATED FORECASTING SYSTEM"

!     AUTHOR.
!     -------
!        C.JAKOB   *ECMWF*

!     MODIFICATIONS.
!     --------------
!        ORIGINAL : 94-02-07
!     ------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB

USE YOMCST   , ONLY : RG
USE YOECLDP  , ONLY : RAMID    ,RCLDIFF  ,RCLCRIT  ,RKCONV   ,&
            &RPRC1    ,RPRC2    ,RCLDMAX  ,RPECONS  ,RTAUMEL  ,&
            &RENTRTU  ,RENTRRA  ,RAMIN    ,RLMIN    ,RASMICE  ,&
            &RBSMICE, RSATQ



!*       1.    SET VALUES
!              ----------



IMPLICIT NONE
RAMID=0.8_JPRB
RCLDIFF=2.E-6_JPRB  !1e-6 pre 25r1 default 
rcldiff=3.e-6_jprb  

RCLCRIT=0.3E-3_JPRB
RKCONV=1.E-4_JPRB
RPRC1=100._JPRB
RPRC2=0.5_JPRB
RCLDMAX=5.E-3_JPRB

RPECONS=5.44E-4_JPRB/RG
RTAUMEL=5._JPRB*3.6E3_JPRB*1.5_JPRB

RENTRTU=0.5_JPRB
RENTRRA=0.5_JPRB

RAMIN=1.E-8_JPRB
RLMIN=1.E-8_JPRB

RASMICE=0.252_JPRB
RBSMICE=0.837_JPRB

RSATQ=1.0E-7_JPRB

!     -----------------------------------------------------------------

RETURN
END SUBROUTINE SUCLDP



MODULE YOEVDFS

USE PARKIND1  ,ONLY : JPIM     ,JPRB

IMPLICIT NONE

SAVE

!     ------------------------------------------------------------------
!*    ** *YOEVDFS* CONTAINS STABILITY FUNCTION TABLES FOR *VDF...*
!     ------------------------------------------------------------------

INTEGER(KIND=JPIM), PARAMETER :: JPRITBL=101
REAL(KIND=JPRB) :: RITBL(JPRITBL)
REAL(KIND=JPRB) :: ARITBL(JPRITBL)
REAL(KIND=JPRB) :: RCHBA
REAL(KIND=JPRB) :: RCHBB
REAL(KIND=JPRB) :: RCHBC
REAL(KIND=JPRB) :: RCHBD
REAL(KIND=JPRB) :: RCHB23A
REAL(KIND=JPRB) :: RCHBBCD
REAL(KIND=JPRB) :: RCHBCD
REAL(KIND=JPRB) :: RCHETA
REAL(KIND=JPRB) :: RCHETB
REAL(KIND=JPRB) :: RCHETC
REAL(KIND=JPRB) :: RCHBHDL
REAL(KIND=JPRB) :: RCDHALF
REAL(KIND=JPRB) :: RCDHPI2
REAL(KIND=JPRB) :: RIMAX
REAL(KIND=JPRB) :: DRITBL
REAL(KIND=JPRB) :: DRI26

!**   ** *YOEDFS* CONTAINS STABILITY FUNCTION TABLES FOR *VDF...*

!     A.C.M. BELJAARS   E.C.M.W.F.       26/03/90.

!      NAME      TYPE        PURPOSE
!      ----      ----        -------

!     *RCHBA*     REAL       *CONSTANT A IN *HOLTSLAG AND *DEBRUIN
!                            FUNCTIONS FOR STABLE SITUATIONS
!     *RCHBB*     REAL       *CONSTANT B IN *HB* FUNCTIONS
!     *RCHBC*     REAL       *CONSTANT C IN *HB* FUNCTIONS
!     *RCHBD*     REAL       *CONSTANT D IN *HB* FUNCTIONS
!     *RCHB23A    REAL       2./3.*A IN *HB* FUNCTIONS
!     *RCHBBCD    REAL       B*C/D IN *HB* FUNCTIONS
!     *RCHBCD     REAL       C/D IN *HB* FUNCTIONS
!     *RCHETA*    REAL       CONSTANT IN THE *HOGSTROM *ELLISON *TURNER
!                            FUNCTIONS FOR STABLY STRATIFIED TURBULENCE 
!     *RCHETB*    REAL       CONSTANT IN THE *HET* FUNCTIONS     
!     *RCHETC*    REAL       CONSTANT IN *HET* FUNCTIONS  
!     *RCHBHDL*   REAL       MAXIM ZNLEV/L FOR STABLE BOUNDARY LAYER  
!     *RCDHALF    REAL       CONSTANT IN *DYER AND *HICKS FORMULAE
!                            FOR UNSTABLE SITUATIONS
!     *RCDHPI2    REAL       PI/2.
!     *RIMAX*     REAL       *MAXIMIM RICHARDSON NUMBER TABULATED
!     *DRITBL*    REAL       *INCREMENT OF THE RICHARDSON NUMBER
!                            BETWEEN TABULATED VALUES.
!     *DRI26*     REAL       DRITBL**2/6.
!     *RITBL*     REAL ARRAY *TABULATED ETA-VALUES (Z/L) AS A FUNCTION
!                            OF THE RICHARDSON NUMBER FOR STABLE CASES.
!     *ARITBL*    REAL ARRAY *SECOND DERIVATIVES OF TABULATED FUNCTION
!                            FOR SPLINE INTERPOLATION.
!     ------------------------------------------------------------------
END MODULE YOEVDFS

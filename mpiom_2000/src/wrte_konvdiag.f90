!****************************************************************
!
!**** *WRTE_KONVDIAG* - save information on convective overturning.
!
!     CH,    *MPI-Met, HH*    10.04.01
!
!     Modified
!     --------
!     S.Legutke,        *MPI-MaD, HH*    01.10.01
!     - separate routine extracted from OLLIE (MAIN)
!
!     Purpose
!     -------
!
!
!     Method
!     -------
!
!
!**   Interface.
!     ----------
!
!     *CALL*       *WRTE_KONVDIAG(KYEAR2,KMONT2,KDAYS)*
!
!     *PARAMETER*  *PARAM1.h*     - grid size parameters for ocean model.
!     *COMMON*     *COMMO1.h*     - ocean/sediment tracer arrays.
!     *COMMON*     *UNITS.h*      - std I/O logical units.
!
!**   Interface to calling routine (parameter list):
!     ----------------------------------------------
!
!     *INTEGER* *KYEAR2*   - actual year.
!     *INTEGER* *KMONT2*   - actual month.
!     *INTEGER* *KDAYS*    - actual day.
!
!
!     Externals
!     ---------
!     none.
!
!**************************************************************************
      SUBROUTINE WRTE_KONVDIAG(KYEAR2,KMONT2,KDAYS)

      USE MO_PARAM1
      USE MO_PARALLEL
      USE MO_COMMO1
      USE MO_UNITS
      USE MO_KIND

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: KYEAR2, KMONT2, KDAYS

      INTEGER(KIND=i4) I4I1,I4I2,I4I3,I4I4
      ! FIXME: wouldn't writing in single precision be sufficient?
      REAL(wp) :: AUX(IE,JE)

      aux = REAL(kcondep, wp)

      IF(p_pe==p_io) THEN

         IO_OU_F090=90
         I4I1=((KYEAR2*10000)+(KMONT2*100)+KDAYS)
         I4I2=69
         I4I3=0
         I4I4=(IE_G*JE_G)
         WRITE(IO_OU_F090)I4I1,I4I2,I4I3,I4I4

      ENDIF

      CALL WRITE_SLICE_SP(IO_OU_F090,AUX)

      END

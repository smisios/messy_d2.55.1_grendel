      SUBROUTINE ADISITJ(TH,SH,PA,J)
! ------------------------------------------------------------------------------
!
!**** *ADISIT*  - TRANSFORMS POTENTIAL TO IN-SITU TEMPERATURE.
!
!     MODIFIED
!     --------
!     O. BOEHRINGER     *DKRZ*                   95
!        - THIS VERSION USES ONLY 44 % OF THE CPU OF THE ORIGINAL HOPC VERSION
!     UWE MIKOLAJEWICZ 2/99
!     ==>ONE-DIMENSIONAL ARRAY, MERGE LOOPS
!
!     METHOD.
!     --------
!     TRANSFORMATION FROM POTENTIAL TO IN SITU TEMPERATURE
!     ACCORDING TO Bryden, 1973, "New polynomials for thermal expansion, adiabatic temperature gradient
!     and potential temperature of sea water". Deep Sea Research and Oceanographic Abstracts. 20, 401-408 (GILL P.602).
!     WHICH GIVES THE INVERSE TRANSFORMATION
!     FOR AN APPROXIMATE VALUE, ALL TERMS LINEAR IN T ARE TAKEN
!     AFTER THAT ONE NEWTON STEP.
!     FOR THE CHECK VALUE 8.4678516 THE ACCURACY IS 0.2 MIKROKELVIN.
!
!**   INTERFACE.
!     ----------
!     *CALL* *ADISIT(TH,SH,PA)*       CALLED FROM *OCTHER*.
!
!     *COMMON*    *"PARAM1*            - OCEAN GRID DIMENSIONS.
!
!
!     INPUT:
!     -----
!     *TH*        POTENTIAL TEMPERATURE [DEG C]
!     *SH*        SALINITY  [PSU.]
!     *PA*        PRESSURE  [PA]
!
!     OUTPUT:
!     ------
!     *TH*        IN-SITU  TEMPERATURE [DEG C]
!
! ------------------------------------------------------------------------------
!
      USE mo_kind, ONLY: wp
      USE MO_PARAM1
      USE MO_PARAM3
      IMPLICIT NONE
      !
      REAL(wp), INTENT(in) :: pa, SH(IE,JE)
      INTEGER, INTENT(in) :: j
      REAL(wp), INTENT(inout) :: TH(IE,JE)
      INTEGER :: i
      REAL(wp) :: pr, dc, dv, dvs, fne, fst, qc, qn3, qnq, qv, qvs, t, tpo
      !
      PR=PA
!
!  CHECK VALUES
!     TH(1)=8.4678516
!     SH(1)= 25.
!     PR=1000.
!
      QC = PR*(A1 + PR*(C1 - E1*PR))
      QV = PR*(B1 - D*PR)
      DC = 1._wp + PR*(-A2 + PR*(C2 -E2*PR))
      DV = B2*PR
      QNQ  = -PR*(-A3 + PR*C3)
      QN3  = -PR*A4
!
      DO  I=1,IE
!
      QVS = QV*(SH(I,J) - 35._wp) + QC
      DVS = DV*(SH(I,J) - 35._wp) + DC
      TPO     = TH(I,J)
      TH(I,J) = (TH(I,J) + QVS)/DVS
      T       = TH(I,J)
      FNE     = - QVS + T*(DVS + T*(QNQ + T*QN3)) - TPO
      FST     = DVS + T*(2._wp*QNQ + 3._wp*QN3*T)
      TH(I,J) =TH(I,J)-FNE/FST
      ENDDO
!
      RETURN
      END

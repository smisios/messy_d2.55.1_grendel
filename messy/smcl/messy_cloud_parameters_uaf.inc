!=======================================================================
!
! *** INCLUDE FILE 'PARAMETR.INC'
! *** THIS FILE CONTAINS THE DECLARATIONS OF THE GLOBAL CONSTANTS
!     AND VARIABLES. 
!
! *** WRITTEN BY ATHANASIOS NENES 
! *** MODIFIED BY ATHANASIOS NENES AND PRASHANT KUMAR
!
!=======================================================================
!
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
      PARAMETER (NSMX = 7, NpGauss=10)
      REAL      XGS, WGS
      LOGICAL   CRIT2, CCNSPST

      COMMON /INPUTS/ TEMPER, PRES &
!
     &       /CCNSPC/ DPG(NSMX),  SIG(NSMX), Dpc(NSMX), &
     &                MODE(NSMX), NMD &
!
     &       /ACTVPR/ AKOH, SSPLT, ALFA, BET1, CRIT2, CCNSPST &
!
     &       /THERMO/ AMW, AMA, GRAV, DENW, CPAIR, DHV, RGAS, AKA, &
     &                PSAT, DAIR, SURT, Dw &
!
     &       /SLNPAR/ EPSC, MAXIT, NITER &
!
     &       /GAUSSL/ XGS(NpGauss), WGS(NpGauss) &
!
     &       /OTHER/  PI, ZERO, GREAT, SQ2PI, &
     &                D11,D12,D13,D14,D15, &
     &                D21,D22,D23,D24,D25, &
     &                D31,D32,D33,D34,D35, &
     &                D41,D42,D43,D44,D45

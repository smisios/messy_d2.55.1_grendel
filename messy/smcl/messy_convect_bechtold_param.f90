MODULE MESSY_CONVECT_BECHTOLD_PARAM

  USE messy_main_constants_mem,    only: dp,                                      &
                                         RPI => PI, RG => g, RATM => atm2Pa,      &
                                         R => R_GAS, RMD => M_AIR, RMV => M_H2O

  IMPLICIT NONE

  SAVE

!     ------------------------------------------------------------------
!     ###################
!        MODULE YOMCST
!     ###################
!*    Common of physical constants
!     You will find the meanings in the annex 1 of the documentation

!!$! A1.0 Fundamental constants
!!$REAL(dp) :: RPI
!!$REAL(dp) :: RCLUM
!!$REAL(dp) :: RHPLA
!!$REAL(dp) :: RKBOL
!!$REAL(dp) :: RNAVO
!!$! A1.1 Astronomical constants
!!$REAL(dp) :: RDAY
!!$REAL(dp) :: REA
!!$REAL(dp) :: REPSM
!!$REAL(dp) :: RSIYEA
!!$REAL(dp) :: RSIDAY
!!$REAL(dp) :: ROMEGA
!!$! A1.2 Geoide
!!$REAL(dp) :: RA
!!$REAL(dp) :: RG
!!$REAL(dp) :: R1SA
!!$! A1.3 Radiation
!!$REAL(dp) :: RSIGMA
!!$REAL(dp) :: RI0

! A1.4 Thermodynamic gas phase

!!$REAL(dp) :: R
!!$REAL(dp) :: RMD
!!$REAL(dp) :: RMV
!!$REAL(dp) :: RMO3

  REAL(dp), PARAMETER :: RD   = 1000._dp * R / RMD
  REAL(dp), PARAMETER :: RV   = 1000._dp * R / RMV
  REAL(dp), PARAMETER :: RCPD = 3.5_dp * RD
  REAL(dp), PARAMETER :: RCPV = 4.0_dp * RV

!!$REAL(dp) :: RCVD
!!$REAL(dp) :: RCVV
!!$REAL(dp) :: RKAPPA
!!$REAL(dp) :: RETV

! A1.5,6 Thermodynamic liquid,solid phases

  REAL(dp), PARAMETER :: RCW  = 4218._dp
  REAL(dp), PARAMETER :: RCS  = 2106._dp

! A1.7 Thermodynamic transition of phase

  REAL(dp), PARAMETER :: RLVTT = 2.5008E+6_dp
  REAL(dp), PARAMETER :: RLSTT = 2.8345E+6_dp

!!$REAL(dp) :: RLVZER
!!$REAL(dp) :: RLSZER
!!$REAL(dp) :: RLMLT

  REAL(dp), PARAMETER :: RTT  = 273.16_dp
!!$REAL(dp) :: RATM
!!$REAL(dp) :: RDT
!!$! A1.8 Curve of saturation
  REAL(dp), PARAMETER :: RESTT = 611.14_dp
  REAL(dp), PARAMETER :: RGAMW = (RCW - RCPV)/RV
  REAL(dp), PARAMETER :: RGAMS = (RCS - RCPV)/RV
  REAL(dp), PARAMETER :: RBETS = RLSTT/RV + RGAMS*RTT
  REAL(dp), PARAMETER :: RBETW = RLVTT/RV + RGAMW*RTT
  REAL(dp) :: RALPW 
  REAL(dp) :: RALPS
!----------------------------------------------------------------------

  INTEGER  :: JCVEXB ! start vertical computations at
                     ! 1 + JCVEXB = 1 + ( KBDIA - 1 )
  INTEGER  :: JCVEXT ! limit vertical computations to
                     ! KLEV - JCVEXT = KLEV - ( KTDIA - 1 )

!-----------------------------------------------------------------------
!     ###################
!     MODULE YOE_CONVPAR
!     ###################

!!****  *YOE_CONVPAR* - Declaration of convection constants
!!
!!    PURPOSE
!!    -------
!      The purpose of this declarative module is to declare  the
!      constants in the deep convection parameterization.

!!
!!**  IMPLICIT ARGUMENTS
!!    ------------------
!!      None
!!
!!    REFERENCE
!!    ---------
!!      Book2 of documentation of Meso-NH (YOE_CONVPAR)
!!
!!    AUTHOR
!!    ------
!!      P. Bechtold   *Laboratoire d'Aerologie*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    26/03/96
!!   Last modified  15/11/96
!-------------------------------------------------------------------------------


!*       0.   DECLARATIONS
!             ------------


  REAL(dp) :: XA25        ! 25 km x 25 km reference grid area
  
  REAL(dp) :: XCRAD       ! cloud radius
  REAL(dp) :: XCDEPTH     ! minimum necessary cloud depth
  REAL(dp) :: XENTR       ! entrainment constant (m/Pa) = 0.2 (m)
  
  REAL(dp) :: XZLCL       ! maximum allowed allowed height
                          ! difference between departure level and surface
  REAL(dp) :: XZPBL       ! minimum mixed layer depth to sustain convection
  REAL(dp) :: XWTRIG      ! constant in vertical velocity trigger
  REAL(dp) :: XDTHPBL     ! temperature perturbation in PBL for trigger
  REAL(dp) :: XDRVPBL     ! moisture perturbation in PBL for trigger

  REAL(dp) :: XNHGAM      ! accounts for non-hydrost. pressure
                          ! in buoyancy term of w equation
                          ! = 2 / (1+gamma)
  REAL(dp) :: XTFRZ1      ! begin of freezing interval
  REAL(dp) :: XTFRZ2      ! end of freezing interval
  
  REAL(dp) :: XRHDBC      ! relative humidity below cloud in downdraft

  REAL(dp) :: XRCONV      ! constant in precipitation conversion
  REAL(dp) :: XSTABT      ! factor to assure stability in  fractional time
                          ! integration, routine CONVECT_CLOSURE
  REAL(dp) :: XSTABC      ! factor to assure stability in CAPE adjustment,
                          !  routine CONVECT_CLOSURE
  REAL(dp) :: XUSRDPTH    ! pressure thickness used to compute updraft
                          ! moisture supply rate for downdraft
  REAL(dp) :: XMELDPTH    ! layer (Pa) through which precipitation melt is
                          ! allowed below  melting level
  REAL(dp) :: XUVDP       ! constant for pressure perturb in momentum transport

!----------------------------------------------------------------------------------
!     ########################
!     MODULE YOE_CONVPAR_SHAL
!     ########################

!!****  *YOE_CONVPAR_SHAL* - Declaration of convection constants
!!
!!    PURPOSE
!!    -------
!!      The purpose of this declarative module is to declare  the
!!      constants in the deep convection parameterization.
!!
!!
!!**  IMPLICIT ARGUMENTS
!!    ------------------
!!      None
!!
!!    REFERENCE
!!    ---------
!!      Book2 of documentation of Meso-NH (YOE_CONVPAR_SHAL)
!!
!!    AUTHOR
!!    ------
!!      P. Bechtold   *Laboratoire d'Aerologie*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    26/03/96
!!   Last modified  04/10/98
!-------------------------------------------------------------------------------

!#include "tsmbkind.h"

!*       0.   DECLARATIONS
!             ------------

!IMPLICIT NONE

  REAL(dp) :: XA25_S        ! 25 km x 25 km reference grid area
  
  REAL(dp) :: XCRAD_S       ! cloud radius
  REAL(dp) :: XCTIME_SHAL   ! convective adjustment time
  REAL(dp) :: XCDEPTH_S     ! minimum necessary cloud depth
  REAL(dp) :: XCDEPTH_D_S   ! maximum allowed cloud thickness
  REAL(dp) :: XDTPERT_S     ! add small Temp perturb. at LCL
  REAL(dp) :: XENTR_S       ! entrainment constant (m/Pa) = 0.2 (m)

  REAL(dp) :: XZLCL_S       ! maximum allowed allowed height
                            ! difference between departure level and surface
  REAL(dp) :: XZPBL_S       ! minimum mixed layer depth to sustain convection
  REAL(dp) :: XWTRIG_S      ! constant in vertical velocity trigger

  REAL(dp) :: XNHGAM_S      ! accounts for non-hydrost. pressure
                            ! in buoyancy term of w equation
                            ! = 2 / (1+gamma)
  REAL(dp) :: XTFRZ1_S      ! begin of freezing interval
  REAL(dp) :: XTFRZ2_S      ! end of freezing interval

  REAL(dp) :: XSTABT_S      ! factor to assure stability in  fractional time
                            ! integration, routine CONVECT_CLOSURE
  REAL(dp) :: XSTABC_S      ! factor to assure stability in CAPE adjustment,
                            !  routine CONVECT_CLOSURE

CONTAINS


!---------------------------------------------------------------------------
  SUBROUTINE INITIALIZE_SATW

    RALPW = LOG(RESTT) + RBETW/RTT + RGAMW*LOG(RTT)
    RALPS = LOG(RESTT) + RBETS/RTT + RGAMS*LOG(RTT)

  END SUBROUTINE INITIALIZE_SATW

!---------------------------------------------------------------------------
!     ###########################
      SUBROUTINE SU_CONVPAR_SHAL
!     ###########################

!!****  *SU_CONVPAR * - routine to initialize the constants modules
!!
!!    PURPOSE
!!    -------
!!       The purpose of this routine is to initialize  the constants
!!     stored in  modules YOE_CONVPAR_SHAL
!!
!!
!!**  METHOD
!!    ------
!!      The shallow convection constants are set to their numerical values
!!
!!
!!    EXTERNAL
!!    --------
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!      Module YOE_CONVPAR_SHAL   : contains deep convection constants
!!
!!    REFERENCE
!!    ---------
!!      Book2 of the documentation (module YOE_CONVPAR_SHAL, routine SU_CONVPAR)
!!
!!
!!    AUTHOR
!!    ------
!!      P. BECHTOLD       * Laboratoire d'Aerologie *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    26/03/96
!!   Last modified  15/04/98 adapted for ARPEGE
!-------------------------------------------------------------------------------

!#include "tsmbkind.h"

!*       0.    DECLARATIONS
!              ------------

!USE YOE_CONVPAR_SHAL

IMPLICIT NONE

!-------------------------------------------------------------------------------

!*       1.    Set the thermodynamical and numerical constants for
!              the deep convection parameterization
!              ---------------------------------------------------


XA25_S      = 625.E6_dp ! 25 km x 25 km reference grid area

XCRAD_S     = 50._dp    ! cloud radius (m)
XCTIME_SHAL = 10800._dp ! convective adjustment time (s)
XCDEPTH_S   = 0.5E3_dp  ! minimum necessary shallow cloud depth
XCDEPTH_D_S = 2.0E3_dp  ! maximum allowed shallow cloud depth (m)
XDTPERT_S   = .2_dp     ! add small Temp perturbation at LCL (K)
XENTR_S     = 0.03_dp   ! entrainment constant (m/Pa) = 0.2 (m)

XZLCL_S     = 0.5E3_dp  ! maximum allowed height (m)
                        ! difference between the DPL and the surface
XZPBL_S     = 40.E2_dp  ! minimum mixed layer depth to sustain convection

XNHGAM_S    = 1.3333_dp ! accounts for non-hydrost. pressure
                        ! in buoyancy term of w equation
                        ! = 2 / (1+gamma)
XTFRZ1_S    = 273.16_dp ! begin of freezing interval (K)
XTFRZ2_S    = 250.16_dp ! end of freezing interval (K)

XSTABT_S    = 0.75_dp   ! factor to assure stability in  fractional time
                        ! integration, routine CONVECT_CLOSURE
XSTABC_S    = 0.95_dp   ! factor to assure stability in CAPE adjustment,
                        !  routine CONVECT_CLOSURE

END SUBROUTINE SU_CONVPAR_SHAL

!=============================================================================

!   ROUTINE in which the constants are defined in original scheme sucst
!!$!---------------------------------------------------------------------------------------
!!$!*       1.    DEFINE FUNDAMENTAL CONSTANTS.
!!$!              -----------------------------
!!$
!!$
!!$RPI=_TWO_*ASIN(_ONE_)
!!$RCLUM=299792458._dp
!!$RHPLA=6.6260755E-34_dp
!!$RKBOL=1.380658E-23_dp
!!$RNAVO=6.0221367E+23_dp
!!$
!!$!     ------------------------------------------------------------------
!!$
!!$!*       2.    DEFINE ASTRONOMICAL CONSTANTS.
!!$!              ------------------------------
!!$
!!$RDAY=86400._dp
!!$REA=149597870000._dp
!!$REPSM=0.409093_dp
!!$
!!$RSIYEA=365.25_dp*RDAY*_TWO_*RPI/6.283076_dp
!!$RSIDAY=RDAY/(_ONE_+RDAY/RSIYEA)
!!$ROMEGA=_TWO_*RPI/RSIDAY
!!$
!!$IDAT=KDAT
!!$ISSS=KSSS
!!$ID=NDD(IDAT)
!!$IM=NMM(IDAT)
!!$IA=NCCAA(IDAT)
!!$ZJU=RJUDAT(IA,IM,ID)
!!$ZTI=RTIME(IA,IM,ID,ISSS)
!!$RTIMST=ZTI
!!$RTIMTR=ZTI
!!$ZTETA=RTETA(ZTI)
!!$ZRS=RRS(ZTETA)
!!$ZDE=RDS(ZTETA)
!!$ZET=RET(ZTETA)
!!$ZRSREL=ZRS/REA
!!$
!!$!     ------------------------------------------------------------------
!!$
!!$!*       3.    DEFINE GEOIDE.
!!$!              --------------
!!$
!!$RG=9.80665_dp
!!$RA=6371229._dp
!!$R1SA=REAL(_ONE_/REAL(RA,KIND(_ONE_)),KIND(R1SA))
!!$
!!$!     ------------------------------------------------------------------
!!$
!!$!*       4.    DEFINE RADIATION CONSTANTS.
!!$!              ---------------------------
!!$
!!$RSIGMA=_TWO_ * RPI**5 * RKBOL**4 /(15._dp* RCLUM**2 * RHPLA**3)
!!$RI0=1370._dp
!!$
!!$!     ------------------------------------------------------------------
!!$
!!$!*       5.    DEFINE THERMODYNAMIC CONSTANTS, GAS PHASE.
!!$!              ------------------------------------------
!!$
!!$R=RNAVO*RKBOL
!!$RMD=28.9644_dp
!!$RMV=18.0153_dp
!!$RMO3=47.9942_dp
!!$RD=1000._dp*R/RMD
!!$RV=1000._dp*R/RMV
!!$RCPD=3.5_dp*RD
!!$RCVD=RCPD-RD
!!$RCPV=4._dp *RV
!!$RCVV=RCPV-RV
!!$RKAPPA=RD/RCPD
!!$RETV=RV/RD-_ONE_
!!$
!!$!     ------------------------------------------------------------------
!!$
!!$!*       6.    DEFINE THERMODYNAMIC CONSTANTS, LIQUID PHASE.
!!$!              ---------------------------------------------
!!$
!!$RCW=4218._dp
!!$
!!$!     ------------------------------------------------------------------
!!$
!!$!*       7.    DEFINE THERMODYNAMIC CONSTANTS, SOLID PHASE.
!!$!              --------------------------------------------
!!$
!!$RCS=2106._dp
!!$
!!$!     ------------------------------------------------------------------
!!$
!!$!*       8.    DEFINE THERMODYNAMIC CONSTANTS, TRANSITION OF PHASE.
!!$!              ----------------------------------------------------
!!$
!!$RTT=273.16_dp
!!$RDT=11.82_dp
!!$RLVTT=2.5008E+6_dp
!!$RLSTT=2.8345E+6_dp
!!$RLVZER=RLVTT+RTT*(RCW-RCPV)
!!$RLSZER=RLSTT+RTT*(RCS-RCPV)
!!$RLMLT=RLSTT-RLVTT
!!$RATM=100000._dp
!!$
!!$!     ------------------------------------------------------------------
!!$
!!$!*       9.    SATURATED VAPOUR PRESSURE.
!!$!              --------------------------
!!$
!!$RESTT=611.14_dp
!!$RGAMW=(RCW-RCPV)/RV
!!$RBETW=RLVTT/RV+RGAMW*RTT
!!$RALPW=LOG(RESTT)+RBETW/RTT+RGAMW*LOG(RTT)
!!$RGAMS=(RCS-RCPV)/RV
!!$RBETS=RLSTT/RV+RGAMS*RTT
!!$RALPS=LOG(RESTT)+RBETS/RTT+RGAMS*LOG(RTT)
!!$RGAMD=RGAMS-RGAMW
!!$RBETD=RBETS-RBETW
!!$RALPD=RALPS-RALPW
!!$
!!$!     ------------------------------------------------------------------

!==============================================================================

!     ######################
      SUBROUTINE SU_CONVPAR
!     ######################

!!****  *SU_CONVPAR * - routine to initialize the constants modules
!!
!!    PURPOSE
!!    -------
!       The purpose of this routine is to initialize  the constants
!     stored in  modules YOE_CONVPAR, YOMCST, YOE_CONVPAREXT.


!!**  METHOD
!!    ------
!!      The deep convection constants are set to their numerical values
!!
!!
!!    EXTERNAL
!!    --------
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!      Module YOE_CONVPAR   : contains deep convection constants
!!
!!    REFERENCE
!!    ---------
!!      Book2 of the documentation (module YOE_CONVPAR, routine SU_CONVPAR)
!!
!!
!!    AUTHOR
!!    ------
!!      P. BECHTOLD       * Laboratoire d'Aerologie *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    26/03/96
!!   Last modified  15/04/98 adapted for ARPEGE
!-------------------------------------------------------------------------------

!#include "tsmbkind.h"

!*       0.    DECLARATIONS
!              ------------

!USE YOE_CONVPAR

        IMPLICIT NONE

!-------------------------------------------------------------------------------

!*       1.    Set the thermodynamical and numerical constants for
!              the deep convection parameterization
!              ---------------------------------------------------


        XA25     = 625.E6_dp    ! 25 km x 25 km reference grid area

        XCRAD    = 1500._dp     ! cloud radius (m)
        XCDEPTH  = 2.E3_dp      ! minimum necessary cloud depth (m)
 ! mz_ht_20070911+
        ! original value
        XENTR    = 0.03_dp      ! entrainment constant (m/Pa) = 0.2 (m)
!        XENTR    = 0.003_dp      ! entrainment constant (m/Pa) = 0.2 (m)
 ! mz_ht_20070911-

        XZLCL    = 3.5E3_dp     ! maximum allowed allowed height
                                ! difference between the surface and the DPL (m)
        XZPBL    = 60.E2_dp     ! minimum mixed layer depth to sustain convection
        XWTRIG   = 6.00_dp      ! constant in vertical velocity trigger
        XDTHPBL  = .2_dp        ! Temp. perturbation in PBL for trigger (K)
        XDRVPBL  = 1.e-4_dp     ! moisture  perturbation in PBL for trigger (kg/kg)

        XNHGAM   = 1.3333_dp    ! accounts for non-hydrost. pressure
                                ! in buoyancy term of w equation
                                ! = 2 / (1+gamma)
        XTFRZ1   = 273.16_dp    ! begin of freezing interval (K)
        XTFRZ2   = 250.16_dp    ! end of freezing interval (K)

        XRHDBC   = 0.9_dp       ! relative humidity below cloud in downdraft

        ! mz_ht_20070911+
        ! original value
        XRCONV   = 0.015_dp     ! constant in precipitation conversion
        ! mz_ht_20070911-
        XSTABT   = 0.75_dp      ! factor to assure stability in  fractional time
                                ! integration, routine CONVECT_CLOSURE
        ! mz_ht_20070912+
        ! original value
        XSTABC   = 0.95_dp      ! factor to assure stability in CAPE adjustment,
                                !  routine CONVECT_CLOSURE
        ! mz_ht_20070912-
        XUSRDPTH = 165.E2_dp    ! pressure thickness used to compute updraft
                                ! moisture supply rate for downdraft
        XMELDPTH = 200.E2_dp    ! layer (Pa) through which precipitation melt is
                                ! allowed below downdraft
        XUVDP    = 0.9_dp       ! constant for pressure perturb in momentum transport


END SUBROUTINE SU_CONVPAR
!==============================================================================
!     #######################
      SUBROUTINE SU_CONVPAR1
!     #######################

!!****  *SU_CONVPAR * - routine to initialize the constants modules
!!
!!    PURPOSE
!!    -------
!       The purpose of this routine is to initialize  the constants
!     stored in  modules YOE_CONVPAR


!!**  METHOD
!!    ------
!!      The deep convection constants are set to their numerical values
!!
!!
!!    EXTERNAL
!!    --------
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!      Module YOE_CONVPAR   : contains deep convection constants
!!
!!    REFERENCE
!!    ---------
!!
!!    AUTHOR
!!    ------
!!      P. BECHTOLD       * Laboratoire d'Aerologie *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    26/03/96
!!   Last modified  15/04/98 adapted for ARPEGE
!-------------------------------------------------------------------------------

!#include "tsmbkind.h"

!*       0.    DECLARATIONS
!              ------------

!USE YOE_CONVPAR

        IMPLICIT NONE

!-------------------------------------------------------------------------------

!*       1.    Set the thermodynamical and numerical constants for
!              the deep convection parameterization
!              ---------------------------------------------------


        XA25     = 625.E6_dp    ! 25 km x 25 km reference grid area

        XCRAD    =  500._dp     ! cloud radius (m)
!       XCRAD    = 1500._dp     ! cloud radius
        XCDEPTH  = 3.E3_dp      ! minimum necessary cloud depth
        XENTR    = 0.03_dp      ! entrainment constant (m/Pa) = 0.2 (m)

        XZLCL    = 3.5E3_dp     ! maximum allowed allowed height
                                ! difference between the surface and the DPL (m)
        XZPBL    = 60.E2_dp     ! minimum mixed layer depth to sustain convection
        XWTRIG   = 6.00_dp      ! constant in vertical velocity trigger
        XDTHPBL  = .2_dp        ! Temp. perturbation in PBL for trigger (K)
        XDRVPBL  = 1.e-4_dp     ! moisture  perturbation in PBL for trigger (kg/kg)

        XNHGAM   = 1.3333_dp    ! accounts for non-hydrost. pressure
                                ! in buoyancy term of w equation
                                ! = 2 / (1+gamma)
        XTFRZ1   = 273.16_dp    ! begin of freezing interval (K)
        XTFRZ2   = 250.16_dp    ! end of freezing interval (K)

        XRHDBC   = 0.9_dp       ! relative humidity below cloud in downdraft

        XRCONV   = 0.015_dp     ! constant in precipitation conversion
        XSTABT   = 0.75_dp      ! factor to assure stability in  fractional time
                                  ! integration, routine CONVECT_CLOSURE
        XSTABC   = 0.95_dp      ! factor to assure stability in CAPE adjustment,
                                  !  routine CONVECT_CLOSURE
        XUSRDPTH = 165.E2_dp    ! pressure thickness used to compute updraft
                                  ! moisture supply rate for downdraft
        XMELDPTH = 150.E2_dp    ! layer (Pa) through which precipitation melt is
                                  ! allowed below downdraft
        XUVDP    = 0.9_dp       ! constant for pressure perturb in momentum transport
!
!
END SUBROUTINE SU_CONVPAR1

!==============================================================================
END MODULE MESSY_CONVECT_BECHTOLD_PARAM

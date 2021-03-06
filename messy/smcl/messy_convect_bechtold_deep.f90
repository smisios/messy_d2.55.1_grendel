MODULE MESSY_CONVECT_BECHTOLD_DEEP

  USE MESSY_CONVECT_BECHTOLD_PARAM

  IMPLICIT NONE
  PRIVATE
  SAVE
  PUBLIC::  CONVECT_TRIGGER_FUNCT
  PUBLIC::  CONVECT_SATMIXRATIO
  PUBLIC::  CONVECT_UPDRAFT
  PUBLIC::  CONVECT_CONDENS
  PUBLIC::  CONVECT_MIXING_FUNCT
  PUBLIC::  CONVECT_TSTEP_PREF
  PUBLIC::  CONVECT_DOWNDRAFT
  PUBLIC::  CONVECT_PRECIP_ADJUST
  PUBLIC::  CONVECT_CLOSURE
  PUBLIC::  CONVECT_CLOSURE_THRVLCL
  PUBLIC::  CONVECT_CLOSURE_ADJUST
  PUBLIC::  CONVECT_UV_TRANSPORT

  PUBLIC::  CONVECT_CHEM_TRANSPORT

CONTAINS

!######################################################################
 SUBROUTINE CONVECT_TRIGGER_FUNCT( KLON, KLEV,                        &
                              & PPRES, PTH, PTHV, PTHES,              &
                              & PRV, PW, PZ, PDXDY, PHSFLX,           &
                              & PTHLCL, PTLCL, PRVLCL, PWLCL, PZLCL,  &
                              & PTHVELCL, KLCL, KDPL, KPBL, OTRIG,    &
                              & PCAPE )
!######################################################################

!!**** Determine convective columns as well as the cloudy values of theta,
!!     and qv at the lifting condensation level (LCL)
!!
!!    PURPOSE
!!    -------
!!      The purpose of this routine is to determine convective columns
!!
!!
!!
!!**  METHOD
!!    ------
!!      Computations are done at every model level starting from bottom.
!!      The use of masks allows to optimise the inner loops (horizontal loops).
!!      What we look for is the undermost unstable level at each grid point.
!!
!!
!!
!!    EXTERNAL
!!    --------
!!     Routine CONVECT_SATMIXRATIO
!!
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!      Module YOMCST
!!          RG                 ! gravity constant
!!          RATM               ! Reference pressure
!!          RD, RV           ! Gaz  constants for dry air and water vapor
!!          RCPD               ! Cpd (dry air)
!!          RTT                ! triple point temperature
!!          RBETW, RGAMW      ! constants for vapor saturation pressure
!!
!!      Module YOE_CONVPAR
!!          XA25               ! reference grid area
!!          XZLCL              ! maximum height difference between
!!                             ! the surface and the DPL
!!          XZPBL              ! minimum mixed layer depth to sustain convection
!!          XWTRIG             ! constant in vertical velocity trigger
!!          XCDEPTH            ! minimum necessary cloud depth
!!          XNHGAM             ! coefficient for buoyancy term in w eq.
!!                             ! accounting for nh-pressure
!!          XDTHPBL            ! theta perturbation in PBL
!!          XDRVPBL            ! moisture perturbation in PBL
!!
!!      Module YOE_CONVPAREXT
!!          JCVEXB, JCVEXT     ! extra levels on the vertical boundaries
!!
!!    REFERENCE
!!    ---------
!!
!!      Book2 of documentation ( routine TRIGGER_FUNCT)
!!      Fritsch and Chappell (1980), J. Atm. Sci., Vol. 37, 1722-1761.
!!
!!    AUTHOR
!!    ------
!!      P. BECHTOLD       * Laboratoire d'Aerologie *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    07/11/95
!!   Last modified  20/03/97  Select first departure level
!!                            that produces a cloud thicker than XCDEPTH
!-------------------------------------------------------------------------------

!#include "tsmbkind.h"

!*       0.    DECLARATIONS
!              ------------

!USE YOMCST
!USE YOE_CONVPAR
!USE YOE_CONVPAREXT


IMPLICIT NONE

!*       0.1   Declarations of dummy arguments :

INTEGER,    INTENT(IN)                   :: KLON      ! horizontal loop index
INTEGER,    INTENT(IN)                   :: KLEV      ! vertical loop index
REAL(dp), DIMENSION(KLON),     INTENT(IN) :: PDXDY     ! grid area
REAL(dp), DIMENSION(KLON,KLEV),INTENT(IN) :: PTH, PTHV ! theta, theta_v
REAL(dp), DIMENSION(KLON,KLEV),INTENT(IN) :: PTHES     ! envir. satur. theta_e
REAL(dp), DIMENSION(KLON,KLEV),INTENT(IN) :: PRV       ! vapor mixing ratio
REAL(dp), DIMENSION(KLON,KLEV),INTENT(IN) :: PPRES     ! pressure
REAL(dp), DIMENSION(KLON,KLEV),INTENT(IN) :: PZ        ! height of grid point (m)
REAL(dp), DIMENSION(KLON,KLEV),INTENT(IN) :: PW        ! vertical velocity
REAL(dp), DIMENSION(KLON,KLEV),INTENT(IN) :: PHSFLX    ! turbulent sensible heat flux (W/m^2)

REAL(dp), DIMENSION(KLON),     INTENT(OUT):: PTHLCL    ! theta at LCL
REAL(dp), DIMENSION(KLON),     INTENT(OUT):: PTLCL     ! temp. at LCL
REAL(dp), DIMENSION(KLON),     INTENT(OUT):: PRVLCL    ! vapor mixing ratio at  LCL
REAL(dp), DIMENSION(KLON),     INTENT(OUT):: PWLCL     ! parcel velocity at  LCL
REAL(dp), DIMENSION(KLON),     INTENT(OUT):: PZLCL     ! height at LCL (m)
REAL(dp), DIMENSION(KLON),     INTENT(OUT):: PTHVELCL  ! environm. theta_v at LCL (K)
LOGICAL,  DIMENSION(KLON),     INTENT(OUT):: OTRIG     ! logical mask for convection
INTEGER,  DIMENSION(KLON),     INTENT(INOUT):: KLCL    ! contains vert. index of LCL
INTEGER,  DIMENSION(KLON),     INTENT(INOUT):: KDPL    ! contains vert. index of DPL
INTEGER,  DIMENSION(KLON),     INTENT(INOUT):: KPBL    ! contains index of source layer top
REAL(dp), DIMENSION(KLON),     INTENT(OUT):: PCAPE     ! CAPE (J/kg) for diagnostics

!*       0.2   Declarations of local variables :

INTEGER  :: JKK, JK, JKP, JKM, JKDL, JL, JKT, JT! vertical loop index
INTEGER  :: JI                                  ! horizontal loop index
INTEGER  :: IIE, IKB, IKE                       ! horizontal + vertical loop bounds
REAL(dp) :: ZEPS, ZEPSA                         ! R_d / R_v, R_v / R_d
REAL(dp) :: ZCPORD, ZRDOCP                      ! C_pd / R_d, R_d / C_pd

REAL(dp), DIMENSION(KLON) :: ZTHLCL, ZTLCL, ZRVLCL, & ! locals for PTHLCL,PTLCL
                        &  ZWLCL,  ZZLCL, ZTHVELCL  ! PRVLCL, ....
INTEGER,  DIMENSION(KLON) :: IDPL, IPBL, ILCL      ! locals for KDPL, ...
REAL(dp), DIMENSION(KLON) :: ZPLCL    ! pressure at LCL
REAL(dp), DIMENSION(KLON) :: ZZDPL    ! height of DPL
REAL(dp), DIMENSION(KLON) :: ZTHVLCL  ! theta_v at LCL = mixed layer value
REAL(dp), DIMENSION(KLON) :: ZTMIX    ! mixed layer temperature
REAL(dp), DIMENSION(KLON) :: ZEVMIX   ! mixed layer water vapor pressure
REAL(dp), DIMENSION(KLON) :: ZDPTHMIX, ZPRESMIX ! mixed layer depth and pressure
REAL(dp), DIMENSION(KLON) :: ZCAPE    ! convective available energy (m^2/s^2/g)
REAL(dp), DIMENSION(KLON) :: ZTHEUL   ! updraft equiv. pot. temperature (K)
REAL(dp), DIMENSION(KLON) :: ZLV, ZCPH! specific heats of vaporisation, dry air
REAL(dp), DIMENSION(KLON) :: ZDP      ! pressure between LCL and model layer
REAL(dp), DIMENSION(KLON) :: ZTOP     ! estimated cloud top (m)
REAL(dp), DIMENSION(KLON,KLEV):: ZCAP ! CAPE at every level for diagnostics
!INTEGER,    DIMENSION(KLON) :: ITOP  ! work array to store highest test layer
REAL(dp), DIMENSION(KLON) :: ZWORK1, ZWORK2, ZWORK3, ZWORK4 ! work arrays
LOGICAL,  DIMENSION(KLON) :: GTRIG, GTRIG2          ! local arrays for OTRIG
LOGICAL,  DIMENSION(KLON) :: GWORK1                 ! work array


!-------------------------------------------------------------------------------

!*       0.3    Compute array bounds
!               --------------------

IIE = KLON
IKB = 1 + JCVEXB
IKE = KLEV - JCVEXT


!*       1.     Initialize local variables
!               --------------------------

ZEPS       = RD  / RV
ZEPSA      = RV  / RD
ZCPORD     = RCPD / RD
ZRDOCP     = RD  / RCPD

OTRIG(:)   = .FALSE.

IDPL(:)    = KDPL(:)
IPBL(:)    = KPBL(:)
ILCL(:)    = KLCL(:)
!ITOP(:)    = IKB

PWLCL(:)   = 0.0_dp
ZWLCL(:)   = 0.0_dp
PTHLCL(:)  = 1.0_dp
PTHVELCL(:)= 1.0_dp
PTLCL(:)   = 1.0_dp
PRVLCL(:)  = 0.0_dp
PWLCL(:)   = 0.0_dp
PZLCL(:)   = PZ(:,IKB)
ZZDPL(:)   = PZ(:,IKB)
GTRIG2(:)  = .TRUE.
ZCAP(:,:)  = 0.0_dp



!       1.     Determine highest necessary loop test layer
!              -------------------------------------------

JT = IKE - 2
DO JK = IKB + 1, IKE - 2
 ! DO JI = 1, IIE
 !    IF ( PZ(JI,JK) - PZ(JI,IKB) <= XZLCL ) ITOP(JI) = JK
 ! ENDDO
   IF ( PZ(1,JK) - PZ(1,IKB) < 12.E3_dp ) JT = JK
ENDDO


!*       2.     Enter loop for convection test
!               ------------------------------

JKP  = MINVAL( IDPL(:) ) + 1
!JKT = MAXVAL( ITOP(:) )
JKT  = JT

DO JKK = JKP, JKT

     GWORK1(:) = ZZDPL(:) - PZ(:,IKB) < XZLCL
          ! we exit the trigger test when the center of the mixed layer is more
          ! than 3500 m  above soil level.
     WHERE ( GWORK1(:) )
        ZDPTHMIX(:) = 0.0_dp
        ZPRESMIX(:) = 0.0_dp
        ZTHLCL(:)   = 0.0_dp
        ZRVLCL(:)   = 0.0_dp
        ZZDPL(:)    = PZ(:,JKK)
        IDPL(:)     = JKK
     END WHERE


!*       3.     Construct a mixed layer of at least 60 hPa (XZPBL)
!               ------------------------------------------

     DO JK = JKK, IKE - 1
       JKM = JK + 1
       DO JI = 1, IIE
         IF ( GWORK1(JI) .AND. ZDPTHMIX(JI) < XZPBL ) THEN
            IPBL(JI)     = JK
            ZWORK1(JI)   = PPRES(JI,JK) - PPRES(JI,JKM)
            ZDPTHMIX(JI) = ZDPTHMIX(JI) + ZWORK1(JI)
            ZPRESMIX(JI) = ZPRESMIX(JI) + PPRES(JI,JK) * ZWORK1(JI)
            ZTHLCL(JI)   = ZTHLCL(JI)   + PTH(JI,JK)   * ZWORK1(JI)
            ZRVLCL(JI)   = ZRVLCL(JI)   + PRV(JI,JK)   * ZWORK1(JI)
         ENDIF
       ENDDO
        IF ( MINVAL ( ZDPTHMIX(:) ) >= XZPBL ) EXIT
     ENDDO


     WHERE ( GWORK1(:) )

        ZPRESMIX(:) = ZPRESMIX(:) / ZDPTHMIX(:)
      ! ZTHLCL(:)   = ZTHLCL(:)   / ZDPTHMIX(:)
      ! ZRVLCL(:)   = ZRVLCL(:)   / ZDPTHMIX(:)
        ZTHLCL(:)   = ZTHLCL(:)   / ZDPTHMIX(:) + XDTHPBL
        ZRVLCL(:)   = ZRVLCL(:)   / ZDPTHMIX(:) + XDRVPBL
        ZTHVLCL(:)  = ZTHLCL(:) * ( 1.0_dp + ZEPSA * ZRVLCL(:) )                 &
                    &           / ( 1.0_dp + ZRVLCL(:) )

!*       4.1    Use an empirical direct solution ( Bolton formula )
!               to determine temperature and pressure at LCL.
!               Nota: the adiabatic saturation temperature is not
!                     equal to the dewpoint temperature
!               ----------------------------------------------------


        ZTMIX(:)  = ZTHLCL(:) * ( ZPRESMIX(:) / RATM ) ** ZRDOCP
        ZEVMIX(:) = ZRVLCL(:) * ZPRESMIX(:) / ( ZRVLCL(:) + ZEPS )
        ZEVMIX(:) = MAX( 1.E-8_dp, ZEVMIX(:) )
        ZWORK1(:) = LOG( ZEVMIX(:) / 613.3_dp )
              ! dewpoint temperature
        ZWORK1(:) = ( 4780.8_dp - 32.19_dp * ZWORK1(:) ) / ( 17.502_dp - ZWORK1(:) )
              ! adiabatic saturation temperature
        ZTLCL(:)  = ZWORK1(:) - ( .212_dp + 1.571E-3_dp * ( ZWORK1(:) - RTT )      &
                  & - 4.36E-4_dp * ( ZTMIX(:) - RTT ) ) * ( ZTMIX(:) - ZWORK1(:) )
        ZTLCL(:)  = MIN( ZTLCL(:), ZTMIX(:) )
        ZPLCL(:)  = RATM * ( ZTLCL(:) / ZTHLCL(:) ) ** ZCPORD

     END WHERE


!*       4.2    Correct ZTLCL in order to be completely consistent
!               with MNH saturation formula
!               ---------------------------------------------

     CALL CONVECT_SATMIXRATIO( KLON, ZPLCL, ZTLCL, ZWORK1, ZLV, ZWORK2, ZCPH )
     WHERE( GWORK1(:) )
        ZWORK2(:) = ZWORK1(:) / ZTLCL(:) * ( RBETW / ZTLCL(:) - RGAMW ) ! dr_sat/dT
        ZWORK2(:) = ( ZWORK1(:) - ZRVLCL(:) ) /                              &
                  &     ( 1.0_dp + ZLV(:) / ZCPH(:) * ZWORK2(:) )
        ZTLCL(:)  = ZTLCL(:) - ZLV(:) / ZCPH(:) * ZWORK2(:)

     END WHERE


!*       4.3    If ZRVLCL = PRVMIX is oversaturated set humidity
!               and temperature to saturation values.
!               ---------------------------------------------

     CALL CONVECT_SATMIXRATIO( KLON, ZPRESMIX, ZTMIX, ZWORK1, ZLV, ZWORK2, ZCPH )
     WHERE( GWORK1(:) .AND. ZRVLCL(:) > ZWORK1(:) )
        ZWORK2(:) = ZWORK1(:) / ZTMIX(:) * ( RBETW / ZTMIX(:) - RGAMW ) ! dr_sat/dT
        ZWORK2(:) = ( ZWORK1(:) - ZRVLCL(:) ) /                              &
                  &    ( 1.0_dp + ZLV(:) / ZCPH(:) * ZWORK2(:) )
        ZTLCL(:)  = ZTMIX(:) - ZLV(:) / ZCPH(:) * ZWORK2(:)
        ZRVLCL(:) = ZRVLCL(:) - ZWORK2(:)
        ZPLCL(:)  = ZPRESMIX(:)
        ZTHLCL(:) = ZTLCL(:) * ( RATM / ZPLCL(:) ) ** ZRDOCP
        ZTHVLCL(:)= ZTHLCL(:) * ( 1.0_dp + ZEPSA * ZRVLCL(:) )                   &
                  &           / ( 1.0_dp + ZRVLCL(:) )
     END WHERE


!*        5.1   Determine  vertical loop index at the LCL and DPL
!               --------------------------------------------------

    DO JK = JKK, IKE - 1
       DO JI = 1, IIE
         IF ( ZPLCL(JI) <= PPRES(JI,JK) .AND. GWORK1(JI) ) ILCL(JI) = JK + 1
       ENDDO
    ENDDO


!*        5.2   Estimate height and environm. theta_v at LCL
!               --------------------------------------------------

    DO JI = 1, IIE
        JK   = ILCL(JI)
        JKM  = JK - 1
        ZDP(JI)    = LOG( ZPLCL(JI) / PPRES(JI,JKM) ) /                     &
                   & LOG( PPRES(JI,JK) / PPRES(JI,JKM) )
        ZWORK1(JI) = PTHV(JI,JKM) + ( PTHV(JI,JK) - PTHV(JI,JKM) ) * ZDP(JI)
           ! we compute the precise value of the LCL
           ! The precise height is between the levels ILCL and ILCL-1.
        ZWORK2(JI) = PZ(JI,JKM) + ( PZ(JI,JK) - PZ(JI,JKM) ) * ZDP(JI)
    ENDDO
    WHERE( GWORK1(:) )
        ZTHVELCL(:) = ZWORK1(:)
        ZZLCL(:)    = ZWORK2(:)
    END WHERE


!*       6.     Check to see if cloud is bouyant
!               --------------------------------

!*      6.1    Compute grid scale vertical velocity perturbation term ZWORK1
!               -------------------------------------------------------------

             !  normalize w grid scale to a 25 km refer. grid
     DO JI = 1, IIE
        JK  = ILCL(JI)
        JKM = JK - 1
        JKDL= IDPL(JI)
      ! ZWORK1(JI) =  ( PW(JI,JKM)  + ( PW(JI,JK) - PW(JI,JKM) ) * ZDP(JI) )  &
        ZWORK1(JI) =  ( PW(JI,JK)  +  PW(JI,JKDL)*ZZLCL(JI)/PZ(JI,JKDL) ) * 0.5_dp  &
                   &       * SQRT( PDXDY(JI) / XA25 )
!                  &      - 0.02 * ZZLCL(JI) / XZLCL ! avoid spurious convection
     ENDDO
             ! compute sign of normalized grid scale w
        ZWORK2(:) = SIGN( 1.0_dp, ZWORK1(:) )
!        zwork2(:) = 0.0_dp    !mz_ht_20040818
        ZWORK1(:) = XWTRIG * ZWORK2(:) * ABS( ZWORK1(:) ) ** 0.333_dp       &
                  &        * ( RATM / ZPLCL(:) ) ** ZRDOCP

!*       6.2    Compute parcel vertical velocity at LCL
!               ---------------------------------------

!mz_ht_20040818+
!     DO JI = 1, IIE
!        JKDL = IDPL(JI)
!        ZWORK3(JI) = RG * ZWORK1(JI) * ( ZZLCL(JI) - PZ(JI,JKDL) )       &
!                   &   / ( PTHV(JI,JKDL) + ZTHVELCL(JI) )
!     ENDDO
! mz_ht_20040818-

   ! DO JI = 1, IIE
   !    JKDL = IDPL(JI)
   !    JK   = ILCL(JI)
   !    ZWORK4(JI) = RG/RCPD * 0.5_dp * ( PHSFLX(JI,JK) + PHSFLX(JI,JKDL) )   &
   !               &   * ( ZZLCL(JI) - PZ(JI,JKDL) ) / ZTHVELCL(JI)  
   !    ZWORK4(JI) = 3._dp * MAX( 1.E-3_dp, ZWORK4(JI) ) ** .3333_dp
   ! ENDDO

     WHERE( GWORK1(:) )
       ZWLCL(:)  = 1.0_dp !+ 0.5_dp * ZWORK2(:) * SQRT( ABS( ZWORK3(:) ) ) !mz_ht_20040818
     ! ZWLCL(:)  = ZWORK4(:) + .25_dp * ZWORK2(:) * SQRT( ABS( ZWORK3(:) ) ) ! UPG PB
       GTRIG(:)  = ZTHVLCL(:) - ZTHVELCL(:) + ZWORK1(:) > 0.0_dp .AND.       &
                 & ZWLCL(:) > 0.0_dp
   !    GTRIG(:)=.TRUE.    ! mz_ht_20040818
       ZWLCL(:)=MAX( 0.5_dp, ZWLCL(:))
     END WHERE


!*       6.3    Look for parcel that produces sufficient cloud depth.
!               The cloud top is estimated as the level where the CAPE
!               is smaller  than a given value (based on vertical velocity eq.)
!               --------------------------------------------------------------

     ZTHEUL(:) = ZTLCL(:) * ( ZTHLCL(:) / ZTLCL(:) ) **                          &
               &            ( 1.0_dp - 0.28_dp * ZRVLCL(:) )                    &
               &          * EXP( ( 3374.6525_dp / ZTLCL(:) - 2.5403_dp ) *   &
               &                 ZRVLCL(:) * ( 1.0_dp + 0.81_dp * ZRVLCL(:) ) )

     ZCAPE(:) = 0.0_dp
     ZTOP(:)  = 0.0_dp
     ZWORK3(:)= 0.0_dp
     JKM = MINVAL( ILCL(:) )
     DO JL = JKM, JT
        JK = JL + 1
        DO JI = 1, IIE
           ZWORK1(JI) = ( 2.0_dp * ZTHEUL(JI) /                                &
           & ( PTHES(JI,JK) + PTHES(JI,JL) ) - 1.0_dp ) * ( PZ(JI,JK) - PZ(JI,JL) )
           IF ( JL < ILCL(JI) ) ZWORK1(JI) = 0.0_dp
           ZCAPE(JI)  = ZCAPE(JI) + ZWORK1(JI)
           ZCAP(JI,JKK) = ZCAP(JI,JKK) + RG * MAX( 0.0_dp, ZWORK1(JI) ) ! actual CAPE
           ZWORK2(JI) = XNHGAM * RG * ZCAPE(JI) + 1.05_dp * ZWLCL(JI) * ZWLCL(JI)
               ! the factor 1.05 takes entrainment into account
           ZWORK2(JI) = SIGN( 1.0_dp, ZWORK2(JI) )
           ZWORK3(JI) = ZWORK3(JI) + MIN(0.0_dp, ZWORK2(JI) )
           ZWORK3(JI) = MAX( -1.0_dp, ZWORK3(JI) )
               ! Nota, the factors ZWORK2 and ZWORK3 are only used to avoid
               ! if and goto statements, the difficulty is to extract only
               ! the level where the criterium is first fullfilled
           ZTOP(JI)   = PZ(JI,JL) * 0.5_dp * ( 1.0_dp + ZWORK2(JI) ) * ( 1.0_dp + ZWORK3(JI) ) + &
                      & ZTOP(JI)  * 0.5_dp * ( 1.0_dp - ZWORK2(JI) )
         ENDDO
     ENDDO


     WHERE( ZTOP(:) - ZZLCL(:)  >=  XCDEPTH  .AND. GTRIG(:) .AND. GTRIG2(:) )
        GTRIG2(:)   = .FALSE.
        OTRIG(:)    = GTRIG(:)     ! we  select the first departure level
        PTHLCL(:)   = ZTHLCL(:)    ! that gives sufficient cloud depth
        PRVLCL(:)   = ZRVLCL(:)
        PTLCL(:)    = ZTLCL(:)
        PWLCL(:)    = ZWLCL(:)
        PZLCL(:)    = ZZLCL(:)
        PTHVELCL(:) = ZTHVELCL(:)
        KDPL(:)     = IDPL(:)
        KPBL(:)     = IPBL(:)
        KLCL(:)     = ILCL(:)
     END WHERE

ENDDO

     DO JI = 1, IIE
       PCAPE(JI) = MAXVAL( ZCAP(JI,:) ) ! maximum CAPE for diagnostics
     ENDDO


END SUBROUTINE CONVECT_TRIGGER_FUNCT

!=============================================================================

!     ################################################################
      SUBROUTINE CONVECT_SATMIXRATIO( KLON,                          &
                                    & PPRES, PT, PEW, PLV, PLS, PCPH )
!     ################################################################

!!**** Compute vapor saturation mixing ratio over liquid water
!!
!!
!!    PDRPOSE
!!    -------
!!     The purpose of this routine is to determine saturation mixing ratio
!!     and to return values for L_v L_s and C_ph
!!
!!
!!**  METHOD
!!    ------
!!
!!    EXTERNAL
!!    --------
!!     None
!!
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!      Module YOMCST
!!          RALPW, RBETW, RGAMW ! constants for water saturation pressure
!!          RD, RV             ! gaz  constants for dry air and water vapor
!!          RCPD, RCPV           ! specific heat for dry air and water vapor
!!          RCW, RCS             ! specific heat for liquid water and ice
!!          RTT                  ! triple point temperature
!!          RLVTT, RLSTT         ! vaporization, sublimation heat constant
!!
!!
!!    REFERENCE
!!    ---------
!!
!!      Book1.0_dp_TWO_of documentation ( routine CONVECT_SATMIXRATIO)
!!
!!    AUTHOR
!!    ------
!!      P. BECHTOLD       * Laboratoire d'Aerologie *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    0.0_dp/1.0_dp/0.5_dp
!!   Last modified  0.0_dp/1.0_dp/97
!------------------------- ------------------------------------------------------

!#include "tsmbkind.h"

!*       0.0_dp    DECLARATIONS
!              ------------

!USE YOMCST


IMPLICIT NONE

!*       1.1  Declarations of dummy arguments :


INTEGER,                    INTENT(IN) :: KLON    ! horizontal loop index
REAL(dp), DIMENSION(KLON),  INTENT(IN) :: PPRES   ! pressure
REAL(dp), DIMENSION(KLON),  INTENT(IN) :: PT      ! temperature

REAL(dp), DIMENSION(KLON),  INTENT(OUT):: PEW     ! vapor saturation mixing ratio
REAL(dp), DIMENSION(KLON),  INTENT(OUT):: PLV     ! latent heat L_v
REAL(dp), DIMENSION(KLON),  INTENT(OUT):: PLS     ! latent heat L_s
REAL(dp), DIMENSION(KLON),  INTENT(OUT):: PCPH    ! specific heat C_ph

!*       1.2  Declarations of local variables :

REAL(dp), DIMENSION(KLON)              :: ZT      ! temperature
REAL(dp)                               :: ZEPS           ! R_d / R_v


!-------------------------------------------------------------------------------

    ZEPS      = RD / RV

    ZT(:)     = MIN( 400._dp, MAX( PT(:), 10._dp ) ) ! overflow bound
    PEW(:)    = EXP( RALPW - RBETW / ZT(:) - RGAMW * LOG( ZT(:) ) )
    PEW(:)    = ZEPS * PEW(:) / ( PPRES(:) - PEW(:) )

    PLV(:)    = RLVTT + ( RCPV - RCW ) * ( ZT(:) - RTT ) ! compute L_v
    PLS(:)    = RLSTT + ( RCPV - RCS ) * ( ZT(:) - RTT ) ! compute L_i

    PCPH(:)   = RCPD + RCPV * PEW(:)                     ! compute C_ph

END SUBROUTINE CONVECT_SATMIXRATIO

!============================================================================

!###########################################################################
  SUBROUTINE CONVECT_UPDRAFT( KLON, KLEV,                                  &
                       &  KICE, PPRES, PDPRES, PZ, PTHL, PTHV, PTHES, PRW, &
                       &  PTHLCL, PTLCL, PRVLCL, PWLCL, PZLCL, PTHVELCL,   &
                       &  PMFLCL, OTRIG, KLCL, KDPL, KPBL,                 &
                       &  PUMF, PUER, PUDR, PUTHL, PUTHV, PURW,            &
                       &  PURC, PURI, PURR, PURS, PUPR,                    &
                       &  PUTPR, PCAPE, KCTL, KETL, &
                          zpah )
!#############################################################################

!!**** Compute updraft properties from DPL to CTL.
!!
!!
!!    PURPOSE
!!    -------
!!      The purpose of this routine is to determine updraft properties
!!      ( mass flux, thermodynamics, precipitation )
!!
!!
!!**  METHOD
!!    ------
!!      Computations are done at every model level starting from bottom.
!!      The use of masks allows to optimise the inner loops (horizontal loops).
!!
!!
!!
!!    EXTERNAL
!!    --------
!!     Routine CONVECT_MIXING_FUNCT
!!     Routine CONVECT_CONDENS
!!
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!      Module YOMCST
!!          RG                 ! gravity constant
!!          RATM               ! reference pressure
!!          RD, RV           ! gaz  constants for dry air and water vapor
!!          RCPD, RCPV, RCW    ! Cp of dry air, water vapor and liquid water
!!          RTT                ! triple point temperature
!!          RLVTT              ! vaporisation heat at RTT
!!
!!
!!      Module YOE_CONVPAR
!!          XA25               ! reference grid area
!!          XCRAD              ! cloud radius
!!          XCDEPTH            ! minimum necessary cloud depth
!!          XENTR              ! entrainment constant
!!          XRCONV             ! constant in precipitation conversion
!!          XNHGAM             ! coefficient for buoyancy term in w eq.
!!                             ! accounting for nh-pressure
!!          XTFRZ1             ! begin of freezing interval
!!          XTFRZ2             ! begin of freezing interval
!!
!!     Module YOE_CONVPAREXT
!!          JCVEXB, JCVEXT     ! extra levels on the vertical boundaries
!!
!!    REFERENCE
!!    ---------
!!
!!      Book1,2 of documentation ( routine CONVECT_UPDRAFT)
!!      Kain and Fritsch, 1990, J. Atmos. Sci., Vol.
!!      Kain and Fritsch, 1993, Meteor. Monographs, Vol.
!!
!!    AUTHOR
!!    ------
!!      P. BECHTOLD       * Laboratoire d'Aerologie *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    07/11/95
!!   Last modified  10/12/97
!-------------------------------------------------------------------------------

!#include "tsmbkind.h"

!*       0.    DECLARATIONS
!              ------------

!USE YOMCST
!USE YOE_CONVPAR
!USE YOE_CONVPAREXT

IMPLICIT NONE

!*       0.1   Declarations of dummy arguments :

INTEGER, INTENT(IN)                    :: KLON  ! horizontal dimension
INTEGER, INTENT(IN)                    :: KLEV  ! vertical dimension
INTEGER, INTENT(IN)                    :: KICE  ! flag for ice ( 1 = yes,
                                                  !                0 = no ice )
REAL(dp), DIMENSION(KLON,KLEV), INTENT(IN) :: PTHL  ! grid scale enthalpy (J/kg)
REAL(dp), DIMENSION(KLON,KLEV), INTENT(IN) :: PTHV  ! grid scale theta_v
REAL(dp), DIMENSION(KLON,KLEV), INTENT(IN) :: PTHES ! grid scale saturated theta_e
REAL(dp), DIMENSION(KLON,KLEV), INTENT(IN) :: PRW   ! grid scale total water
                                                  ! mixing ratio
REAL(dp), DIMENSION(KLON,KLEV), INTENT(IN) :: PPRES ! pressure (P)
REAL(dp), DIMENSION(KLON,KLEV), INTENT(IN) :: PDPRES! pressure difference between
                                                  ! bottom and top of layer (Pa)
REAL(dp), DIMENSION(KLON,KLEV), INTENT(IN) :: PZ    ! height of model layer (m)
REAL(dp), DIMENSION(KLON,0:KLEV), INTENT(IN) :: zpah! half level pressure

REAL(dp), DIMENSION(KLON),     INTENT(IN) :: PTHLCL ! theta at LCL
REAL(dp), DIMENSION(KLON),     INTENT(IN) :: PTLCL  ! temp. at LCL
REAL(dp), DIMENSION(KLON),     INTENT(IN) :: PRVLCL ! vapor mixing ratio at  LCL
REAL(dp), DIMENSION(KLON),     INTENT(IN) :: PWLCL  ! parcel velocity at LCL (m/s)
REAL(dp), DIMENSION(KLON),     INTENT(IN) :: PMFLCL ! cloud  base unit mass flux
                                                  ! (kg/s)
REAL(dp), DIMENSION(KLON),     INTENT(IN) :: PZLCL  ! height at LCL (m)
REAL(dp), DIMENSION(KLON),     INTENT(IN) :: PTHVELCL  ! environm. theta_v at LCL (K)
LOGICAL, DIMENSION(KLON),  INTENT(INOUT):: OTRIG! logical mask for convection
INTEGER, DIMENSION(KLON),  INTENT(IN) :: KLCL   ! contains vert. index of LCL
INTEGER, DIMENSION(KLON),  INTENT(IN) :: KDPL   ! contains vert. index of DPL
INTEGER, DIMENSION(KLON),  INTENT(IN) :: KPBL   !  " vert. index of source layertop


INTEGER, DIMENSION(KLON),  INTENT(OUT):: KCTL   ! contains vert. index of CTL
INTEGER, DIMENSION(KLON),  INTENT(OUT):: KETL   ! contains vert. index of 
                                                  !equilibrium (zero buoyancy) level
REAL(dp), DIMENSION(KLON,KLEV), INTENT(OUT):: PUMF  ! updraft mass flux (kg/s)
REAL(dp), DIMENSION(KLON,KLEV), INTENT(OUT):: PUER  ! updraft entrainment (kg/s)
REAL(dp), DIMENSION(KLON,KLEV), INTENT(OUT):: PUDR  ! updraft detrainment (kg/s)
REAL(dp), DIMENSION(KLON,KLEV), INTENT(OUT):: PUTHL ! updraft enthalpy (J/kg)
REAL(dp), DIMENSION(KLON,KLEV), INTENT(OUT):: PUTHV ! updraft theta_v (K)
REAL(dp), DIMENSION(KLON,KLEV), INTENT(OUT):: PURW  ! updraft total water (kg/kg)
REAL(dp), DIMENSION(KLON,KLEV), INTENT(OUT):: PURC  ! updraft cloud water (kg/kg)
REAL(dp), DIMENSION(KLON,KLEV), INTENT(OUT):: PURI  ! updraft cloud ice   (kg/kg)
REAL(dp), DIMENSION(KLON,KLEV), INTENT(OUT):: PURR  ! liquid precipit. (kg/kg)
                                                  ! produced in  model layer
REAL(dp), DIMENSION(KLON,KLEV),   INTENT(OUT)::PURS ! solid precipit. (kg/kg)
                                                  ! produced in  model layer
REAL(dp), DIMENSION(KLON,KLEV),   INTENT(OUT)::PUPR ! updraft precipitation in
                                                  ! flux units (kg water / s)
REAL(dp), DIMENSION(KLON),     INTENT(OUT):: PUTPR  ! total updraft precipitation
                                                  ! in flux units (kg water / s)
REAL(dp), DIMENSION(KLON),     INTENT(OUT):: PCAPE  ! available potent. energy

!*       0.2   Declarations of local variables :

INTEGER :: IIE, IKB, IKE  ! horizontal and vertical loop bounds
INTEGER :: JI             ! horizontal loop index
INTEGER :: JK, JKP, JKM, JK1, JK2, JKMIN  ! vertical loop index
REAL(dp)    :: ZEPSA, ZCVOCD  ! R_v / R_d, C_pv / C_pd
REAL(dp)    :: ZCPORD, ZRDOCP ! C_pd / R_d, R_d / C_pd

REAL(dp), DIMENSION(KLON)    :: ZUT             ! updraft temperature (K)
REAL(dp), DIMENSION(KLON)    :: ZUW1, ZUW2      ! square of updraft vert.
                                              ! velocity at levels k and k+1
REAL(dp), DIMENSION(KLON)    :: ZE1,ZE2,ZD1,ZD2 ! fractional entrainm./detrain
                                              ! rates at levels k and k+1
REAL(dp), DIMENSION(KLON)    :: ZMIXF           ! critical mixed fraction
REAL(dp), DIMENSION(KLON)    :: ZCPH            ! specific heat C_ph
REAL(dp), DIMENSION(KLON)    :: ZLV, ZLS        ! latent heat of vaporis., sublim.
REAL(dp), DIMENSION(KLON)    :: ZURV            ! updraft water vapor at level k+1
REAL(dp), DIMENSION(KLON)    :: ZPI             ! Pi=(P0/P)**(Rd/Cpd)
REAL(dp), DIMENSION(KLON)    :: ZTHEU1, ZTHEU2  ! theta_e for undilute ascent
REAL(dp), DIMENSION(KLON)    :: ZWORK1, ZWORK2, ZWORK3  ! work arrays
REAL(dp), DIMENSION(KLON)    :: ZWORK4, ZWORK5, ZWORK6  ! work arrays
INTEGER, DIMENSION(KLON) :: IWORK           ! wok array
LOGICAL, DIMENSION(KLON)      :: GWORK1, GWORK2, GWORK4, GWORK5 ! work arrays
LOGICAL, DIMENSION(KLON,KLEV) :: GWORK6       ! work array


!-------------------------------------------------------------------------------

!        0.3   Set loop bounds
!              ---------------

IKB = 1 + JCVEXB
IKE = KLEV - JCVEXT
IIE = KLON


!*       1.     Initialize updraft properties and local variables
!               -------------------------------------------------

ZEPSA      = RV  / RD
ZCVOCD     = RCPV / RCPD
ZCPORD     = RCPD / RD
ZRDOCP     = RD  / RCPD

PUMF(:,:)  = 0.0_dp
PUER(:,:)  = 0.0_dp
PUDR(:,:)  = 0.0_dp
PUTHL(:,:) = 0.0_dp
PUTHV(:,:) = 0.0_dp
PURW(:,:)  = 0.0_dp
PURC(:,:)  = 0.0_dp
PURI(:,:)  = 0.0_dp
PUPR(:,:)  = 0.0_dp
PURR(:,:)  = 0.0_dp
PURS(:,:)  = 0.0_dp
PUTPR(:)   = 0.0_dp
ZUW1(:)    = PWLCL(:) * PWLCL(:)
ZUW2(:)    = 0.0_dp
ZE1(:)     = 1.0_dp
ZD1(:)     = 0.0_dp
PCAPE(:)   = 0.0_dp
KCTL(:)    = IKB
KETL(:)    = KLCL(:)
GWORK2(:)  = .TRUE.
GWORK5(:)  = .TRUE.
ZPI(:)     = 1.0_dp
ZWORK3(:)  = 0.0_dp
ZWORK4(:)  = 0.0_dp
ZWORK5(:)  = 0.0_dp
ZWORK6(:)  = 0.0_dp
GWORK1(:)  = .FALSE.
GWORK4(:)  = .FALSE.


!*       1.1    Compute undilute updraft theta_e for CAPE computations
!               Bolton (1980) formula.
!               Define accurate enthalpy for updraft
!               -----------------------------------------------------

ZTHEU1(:) = PTLCL(:) * ( PTHLCL(:) / PTLCL(:) ) ** &
          &            ( 1._dp - 0.28_dp * PRVLCL(:) ) &
          &          * EXP( ( 3374.6525_dp / PTLCL(:) - 2.5403_dp )      &
          &          * PRVLCL(:) * ( 1.0_dp + 0.81_dp * PRVLCL(:) ) )


ZWORK1(:) = ( RCPD + PRVLCL(:) * RCPV ) *  PTLCL(:)   &
          & + ( 1.0_dp + PRVLCL(:) ) * RG * PZLCL(:)


!*       2.     Set updraft properties between DPL and LCL
!               ------------------------------------------

JKP = MAXVAL( KLCL(:) )
JKM = MINVAL( KDPL(:) )
DO JK = JKM, JKP
   DO JI = 1, IIE
    IF ( JK >= KDPL(JI) .AND. JK < KLCL(JI) ) THEN
        PUMF(JI,JK)  = PMFLCL(JI)
        PUTHL(JI,JK) = ZWORK1(JI)
        PUTHV(JI,JK) = PTHLCL(JI) * ( 1.0_dp + ZEPSA * PRVLCL(JI) ) /        &
                     &            ( 1.0_dp + PRVLCL(JI) )
        PURW(JI,JK)  = PRVLCL(JI)
   ENDIF
   ENDDO
ENDDO


!*       3.     Enter loop for updraft computations
!               ------------------------------------

JKMIN = MINVAL( KLCL(:) - 1 )
DO JK = MAX( IKB + 1, JKMIN ), IKE - 1
  ZWORK6(:) = 1.0_dp
  JKP = JK + 1

  GWORK4(:) = JK >= KLCL(:) - 1
  GWORK1(:) = GWORK4(:) .AND. GWORK2(:) ! this mask is used to confine
                           ! updraft computations between the LCL and the CTL

  WHERE( JK == KLCL(:) - 1 ) ZWORK6(:) = 0.0_dp ! factor that is used in buoyancy
                                        ! computation at first level above LCL


!*       4.     Estimate condensate, L_v L_i, Cph and theta_v at level k+1
!               ----------------------------------------------------------

    ZWORK1(:) = PURC(:,JK) + PURR(:,JK)
    ZWORK2(:) = PURI(:,JK) + PURS(:,JK)
    CALL CONVECT_CONDENS( KLON, KICE, PPRES(:,JKP), PUTHL(:,JK), PURW(:,JK),&
                        & ZWORK1, ZWORK2, PZ(:,JKP), GWORK1, ZUT, ZURV,     &
                        & PURC(:,JKP), PURI(:,JKP), ZLV, ZLS, ZCPH )


  ZPI(:) = ( RATM / PPRES(:,JKP) ) ** ZRDOCP
  WHERE ( GWORK1(:) )

    PUTHV(:,JKP) = ZPI(:) * ZUT(:) * ( 1.0_dp + ZEPSA * ZURV(:) )           &
                 &       / ( 1.0_dp + PURW(:,JK) )


!*       5.     Compute square of vertical velocity using entrainment
!               at level k
!               -----------------------------------------------------

    ZWORK3(:) = PZ(:,JKP) - PZ(:,JK) * ZWORK6(:) -         &
              &      ( 1.0_dp - ZWORK6(:) ) * PZLCL(:)        ! level thickness
    ZWORK4(:) = PTHV(:,JK) * ZWORK6(:) +                   &
              &  ( 1.0_dp - ZWORK6(:) ) * PTHVELCL(:)
    ZWORK5(:) = 2.0_dp * ZUW1(:) * PUER(:,JK) / MAX( 0.1_dp, PUMF(:,JK) )
    ZUW2(:)   = ZUW1(:) + ZWORK3(:) * XNHGAM * RG *        &
              &   ( ( PUTHV(:,JK) + PUTHV(:,JKP) ) /       &
              &   ( ZWORK4(:) + PTHV(:,JKP) ) - 1.0_dp )    & ! buoyancy term
              & - ZWORK5(:)                                  ! entrainment term


!*       6.     Update total precipitation: dr_r=(r_c+r_i)*exp(-rate*dz)
!               --------------------------------------------------------

!                    compute level mean vertical velocity
    ZWORK2(:)   = 0.5_dp *                                                    &
                &      ( SQRT( MAX( 1.E-2_dp, ZUW2(:) ) ) +                 &
                &        SQRT( MAX( 1.E-2_dp, ZUW1(:) ) ) )
    PURR(:,JKP) = 0.5_dp * ( PURC(:,JK) + PURC(:,JKP) + PURI(:,JK) + PURI(:,JKP) )&
                &     * ( 1.0_dp - EXP( - XRCONV  * ZWORK3(:) / ZWORK2(:) ) )
    PUPR(:,JKP) = PURR(:,JKP) * PUMF(:,JK) ! precipitation rate ( kg water / s)
    PUTPR(:)    = PUTPR(:) + PUPR(:,JKP)   ! total precipitation rate
    ZWORK2(:)   = PURR(:,JKP) / MAX( 1.E-8_dp, PURC(:,JKP) + PURI(:,JKP) )
    PURR(:,JKP) = ZWORK2(:) * PURC(:,JKP)          ! liquid precipitation
    PURS(:,JKP) = ZWORK2(:) * PURI(:,JKP)          ! solid precipitation


!*       7.     Update r_c, r_i, enthalpy, r_w  for precipitation
!               -------------------------------------------------------

    PURW(:,JKP)  = PURW(:,JK) - PURR(:,JKP) - PURS(:,JKP)
    PURC(:,JKP)  = PURC(:,JKP) - PURR(:,JKP)
    PURI(:,JKP)  = PURI(:,JKP) - PURS(:,JKP)
    PUTHL(:,JKP) = ( RCPD + PURW(:,JKP) * RCPV ) * ZUT(:)                     &
                 & + ( 1.0_dp + PURW(:,JKP) ) * RG * PZ(:,JKP)                 &
                 & - ZLV(:) * PURC(:,JKP) - ZLS(:) * PURI(:,JKP)

    ZUW1(:)      = ZUW2(:)

  END WHERE


!*       8.     Compute entrainment and detrainment using conservative
!               variables adjusted for precipitation ( not for entrainment)
!               -----------------------------------------------------------

!*       8.1    Compute critical mixed fraction by estimating unknown
!               T^mix r_c^mix and r_i^mix from enthalpy^mix and r_w^mix
!               We determine the zero crossing of the linear curve
!               evaluating the derivative using ZMIXF=0.1.
!               -----------------------------------------------------

    ZMIXF(:)  = 0.1_dp   ! starting value for critical mixed fraction
    ZWORK1(:) = ZMIXF(:) * PTHL(:,JKP)                                     &
              &      + ( 1.0_dp - ZMIXF(:) ) * PUTHL(:,JKP) ! mixed enthalpy
    ZWORK2(:) = ZMIXF(:) * PRW(:,JKP)                                      &
              &      + ( 1.0_dp - ZMIXF(:) ) * PURW(:,JKP)  ! mixed r_w

    CALL CONVECT_CONDENS( KLON, KICE, PPRES(:,JKP), ZWORK1, ZWORK2,        &
                        & PURC(:,JKP), PURI(:,JKP), PZ(:,JKP), GWORK1, ZUT,&
                        & ZWORK3, ZWORK4, ZWORK5, ZLV, ZLS, ZCPH )
!        put in enthalpy and r_w and get T r_c, r_i (ZUT, ZWORK4-5)

     ! compute theta_v of mixture
    ZWORK3(:) = ZUT(:) * ZPI(:) * ( 1.0_dp + ZEPSA * (                         &
              & ZWORK2(:) - ZWORK4(:) - ZWORK5(:) ) ) / ( 1.0_dp + ZWORK2(:) )
     ! compute final value of critical mixed fraction using theta_v
     ! of mixture, grid-scale and updraft
    ZMIXF(:) = MAX( 0.0_dp, PUTHV(:,JKP) - PTHV(:,JKP) ) * ZMIXF(:) /          &
             &                ( PUTHV(:,JKP) - ZWORK3(:) + 1.E-10_dp )
    ZMIXF(:) = MAX( 0.0_dp, MIN( 1.0_dp, ZMIXF(:) ) )


!*       8.2     Compute final midlevel values for entr. and detrainment
!                after call of distribution function
!                -------------------------------------------------------


    CALL CONVECT_MIXING_FUNCT ( KLON, ZMIXF, 1, ZE2, ZD2 )
!       Note: routine MIXING_FUNCT returns fractional entrainm/detrainm. rates

! ZWORK1(:) = XENTR * PMFLCL(:) * PDPRES(:,JKP) / XCRAD ! rate of env. inflow
!*MOD
  ZWORK3(:) = PZ(:,JKP) - PZ(:,JK) * ZWORK6(:) -         &
              &      ( 1.0_dp - ZWORK6(:) ) * PZLCL(:)        ! level thickness
  ZWORK1(:) = XENTR * RG / XCRAD * PUMF(:,JK) * ZWORK3(:)
! ZWORK1(:) = XENTR * pumf(:,jk) * PDPRES(:,JKP) / XCRAD ! rate of env. inflow
!*MOD
  ZWORK2(:) = 0.0_dp
  WHERE ( GWORK1(:) ) ZWORK2(:) = 1.0_dp
  ze2=0.5_dp; zd2=0.5_dp  ! modif entrainment=detrainment, this avoids
  ze2=.6; zd2=.7  ! modif entrainment=detrainment, this avoids
                          ! too large mass flux values at upper levels
  WHERE ( PUTHV(:,JKP) > PTHV(:,JKP) )
    PUER(:,JKP) = 0.5_dp * ZWORK1(:) * ( ZE1(:) + ZE2(:) ) * ZWORK2(:)
    PUDR(:,JKP) = 0.5_dp * ZWORK1(:) * ( ZD1(:) + ZD2(:) ) * ZWORK2(:)
  ELSEWHERE
    PUER(:,JKP) = 0.0_dp
    PUDR(:,JKP) = ZWORK1(:) * ZWORK2(:)
  END WHERE

!*       8.3     Determine equilibrium temperature level
!                --------------------------------------

   WHERE ( PUTHV(:,JKP) > PTHV(:,JKP) .AND. JK > KLCL(:) + 1 &
         & .AND. GWORK1(:) )
         KETL(:) = JKP            ! equilibrium temperature level
   END WHERE

!*       8.4     If the calculated detrained mass flux is greater than
!                the total updraft mass flux, or vertical velocity is
!                negative, all cloud mass detrains at previous model level,
!                exit updraft calculations - CTL is attained
!                -------------------------------------------------------

  WHERE( GWORK1(:) )                                                   &
      & GWORK2(:) = PUMF(:,JK) - PUDR(:,JKP) > 10._dp .AND. ZUW2(:) > 0.0_dp
  WHERE ( GWORK2(:) ) KCTL(:) = JKP   ! cloud top level
  
  GWORK1(:) = GWORK2(:) .AND. GWORK4(:)

  IF ( COUNT( GWORK2(:) ) == 0 ) EXIT


!*       9.   Compute CAPE for undilute ascent using theta_e and
!             theta_es instead of theta_v. This estimation produces
!             a significantly larger value for CAPE than the actual one.
!             ----------------------------------------------------------

  WHERE ( GWORK1(:) )

    ZWORK3(:)   = PZ(:,JKP) - PZ(:,JK) * ZWORK6(:) -                      &
                & ( 1.0_dp - ZWORK6(:) ) *  PZLCL(:)              ! level thickness
    ZWORK2(:)   = PTHES(:,JK) + ( 1.0_dp - ZWORK6(:) ) *                   &
  &  ( PTHES(:,JKP) - PTHES(:,JK) ) / ( PZ(:,JKP) - PZ(:,JK) ) *          &
  &  ( PZLCL(:) - PZ(:,JK) ) ! linear interpolation for theta_es at LCL
                             ! ( this is only done for model level just above LCL

    ZWORK1(:) = MAX( 0.01_dp, PUER(:,JK) / MAX( .1_dp, PUMF(:,JK) ) )
  ! ZTHEU2(:) = ( 1.0_dp - ZWORK1(:) ) * ZTHEU1(:) + ZWORK1(:) * PTHES(:,JK)
    ztheu2(:) = ztheu1(:)
    ZWORK1(:) = ( ZTHEU1(:) + ZTHEU2(:) ) / ( ZWORK2(:) + PTHES(:,JKP) ) - 1.0_dp
  ! ZTHEU1(:) = ZTHEU2(:)
    PCAPE(:)  = PCAPE(:) + RG * ZWORK3(:) * MAX( 0.0_dp, ZWORK1(:) )


!*       10.   Compute final values of updraft mass flux, enthalpy, r_w
!              at level k+1
!              --------------------------------------------------------

    PUMF(:,JKP)  = PUMF(:,JK) - PUDR(:,JKP) + PUER(:,JKP)
    PUMF(:,JKP)  = MAX( PUMF(:,JKP), 0.1_dp )
    PUTHL(:,JKP) = ( PUMF(:,JK) * PUTHL(:,JK) +                              &
                 &   PUER(:,JKP) * PTHL(:,JK) - PUDR(:,JKP) * PUTHL(:,JK) )  &
                 &  / PUMF(:,JKP) + PUTHL(:,JKP) - PUTHL(:,JK)
    PURW(:,JKP)  = ( PUMF(:,JK) * PURW(:,JK) +                               &
                 &   PUER(:,JKP) * PRW(:,JK) - PUDR(:,JKP) * PURW(:,JK) )    &
                 &  / PUMF(:,JKP) - PURR(:,JKP) - PURS(:,JKP)

    ZE1(:) = ZE2(:) ! update fractional entrainment/detrainment
    ZD1(:) = ZD2(:)

  END WHERE

ENDDO

!*       12.1    Set OTRIG to False if cloud thickness < XCDEPTH
!                or CAPE < 1
!                ------------------------------------------------

    DO JI = 1, IIE
          JK  = KCTL(JI)
          OTRIG(JI) = PZ(JI,JK) - PZLCL(JI) >= XCDEPTH               &
                    & .AND. PCAPE(JI) > 1.0_dp

    ENDDO
    WHERE( .NOT. OTRIG(:) )
          KCTL(:) = IKB
    END WHERE

KCTL(:) = MIN(IKE-2, KCTL(:) )
KETL(:) = MAX( KETL(:), KLCL(:) + 2 )
KETL(:) = MIN( KETL(:), KCTL(:) )


!*       12.2    If the ETL and CTL are the same detrain updraft mass
!                flux at this level
!                -------------------------------------------------------

ZWORK1(:) = 0.0_dp
WHERE ( KETL(:) == KCTL(:) ) ZWORK1(:) = 1.0_dp

DO JI = 1, IIE
    JK = KETL(JI)
    PUDR(JI,JK)   = PUDR(JI,JK) +                                    &
                  &       ( PUMF(JI,JK) - PUER(JI,JK) )  * ZWORK1(JI)
    PUER(JI,JK)   = PUER(JI,JK) * ( 1.0_dp - ZWORK1(JI) )
    PUMF(JI,JK)   = PUMF(JI,JK) * ( 1.0_dp - ZWORK1(JI) )
    JKP = KCTL(JI) + 1
    PUER(JI,JKP)  = 0.0_dp ! entrainm/detr rates have been already computed
    PUDR(JI,JKP)  = 0.0_dp ! at level KCTL+1, set them to zero
    PURW(JI,JKP)  = 0.0_dp
    PURC(JI,JKP)  = 0.0_dp
    PURI(JI,JKP)  = 0.0_dp
    PUTHL(JI,JKP) = 0.0_dp
    PURC(JI,JKP+1)= 0.0_dp
    PURI(JI,JKP+1)= 0.0_dp
ENDDO

!*       12.3    Adjust mass flux profiles, detrainment rates, and
!                precipitation fallout rates to reflect linear decrease
!                in mass flux between the ETL and CTL
!                -------------------------------------------------------

ZWORK1(:) = 0.0_dp
JK1 = MINVAL( KETL(:) )
JK2 = MAXVAL( KCTL(:) )
DO JK = JK1, JK2
    DO JI = 1, IIE
    IF( JK > KETL(JI) .AND. JK <= KCTL(JI) ) THEN
        ZWORK1(JI) = ZWORK1(JI) + PDPRES(JI,JK)
    ENDIF
    ENDDO
ENDDO

DO JI = 1, IIE
    JK = KETL(JI)
    ZWORK1(JI) = PUMF(JI,JK) / MAX( 1.0_dp, ZWORK1(JI) )
ENDDO

DO JK = JK1 + 1, JK2
    JKP = JK - 1
    DO JI = 1, IIE
    IF ( JK > KETL(JI) .AND. JK <= KCTL(JI) ) THEN
      ! PUTPR(JI)    = PUTPR(JI) - ( PURR(JI,JK) + PURS(JI,JK) ) * PUMF(JI,JKP)
        PUTPR(JI)    = PUTPR(JI) - PUPR(JI,JK)
        PUDR(JI,JK)  = PDPRES(JI,JK) * ZWORK1(JI)
        PUMF(JI,JK)  = PUMF(JI,JKP) - PUDR(JI,JK)
        PUPR(JI,JK)  = PUMF(JI,JKP) * ( PURR(JI,JK) + PURS(JI,JK) )
        PUTPR(JI)    = PUTPR(JI) + PUPR(JI,JK)
    ENDIF
    ENDDO
ENDDO

do ji = 1, iie
  jk=kctl(ji)+1
  putpr(ji)=putpr(ji)-pupr(ji,jk)
  pupr(ji,jk)=0.0_dp
enddo

!         12.4   Set mass flux and entrainment in the source layer.
!                Linear increase throughout the source layer.
!                -------------------------------------------------------

IWORK=KLCL(:)-1
DO JI = 1, IIE
     JK  = IKB
     JKP = IWORK(JI)
!          mixed layer depth
     ZWORK2(JI) = PPRES(JI,JK) - PPRES(JI,JKP)
ENDDO

JKP = MAXVAL( IWORK(:) )

     ! mz_ht_20050913+
DO JK = IKB+1, JKP
   DO JI = 1, IIE
   IF ( JK >= IKB+1  .AND. JK <= IWORK(JI) ) THEN
      PUER(JI,JK) = PMFLCL(JI) * (PPRES(JI,JK-1)-PPRES(JI,JK)) / ( ZWORK2(JI) + 0.1 )
      PUMF(JI,JK) = PUMF(JI,JK-1) + PUER(JI,JK)

!!$DO JK = IKB, JKP
!!$  DO JI = 1, IIE
!!$    IF ( JK <= IWORK(JI) ) THEN
!!$     PUER(JI,JK) = PMFLCL(JI) * PDPRES(JI,JK) / ( ZWORK2(JI) + 0.1 )
!!$     IF (IKB == 1) THEN
!!$       PUMF(JI,JK) = PUER(JI,JK)
!!$     ELSE
!!$       PUMF(JI,JK) = PUMF(JI,JK-1) + PUER(JI,JK)
!!$     ENDIF

!!$DO JK = IKB, JKP
!!$  DO JI = 1, IIE
!!$    IF ( JK <= IWORK(JI) ) THEN
!!$     IF (IKB == 1) THEN
!!$       PUER(JI,JK) = PMFLCL(JI) * PDPRES(JI,JK) / ( ZWORK2(JI) + 0.1 )
!!$       PUMF(JI,JK) = PUER(JI,JK)
!!$     ELSE
!!$       PUER(JI,JK) = PMFLCL(JI) * (PPRES(JI,JK-1) - PPRES(JI,JK)) / &
!!$                   ( ZWORK2(JI) + 0.1 )
!!$       PUMF(JI,JK) = PUMF(JI,JK-1) + PUER(JI,JK)
!!$     ENDIF
     ! mz_ht_20050913-

   ENDIF
   ENDDO
ENDDO
! mz_ht_20050915+
  DO JI = 1, IIE
    PUER(JI, IKB) = PMFLCL(JI) * (ZPAH(JI,IKB-1)-ZPAH(JI,IKB)) / ( ZWORK2(JI) + 0.1 )
    PUMF(JI, IKB) = PUER(JI,IKB)
  ENDDO
! mz_ht_20050915-
DO JI = 1, IIE
   JK = KLCL(JI)
    PUDR(JI,JK) =- PUMF(JI,JK) + PUMF(JI,JK-1) + PUER(JI,JK)
END DO

!*       13.   If cloud thickness is smaller than  3 km, no
!              convection is allowed
!              Nota: For technical reasons, we stop the convection
!                    computations in this case and do not go back to
!                    TRIGGER_FUNCT to look for the next unstable LCL
!                    which could produce a thicker cloud.
!              ---------------------------------------------------

GWORK6(:,:) = SPREAD( OTRIG(:), DIM=2, NCOPIES=KLEV )
WHERE ( .NOT. OTRIG(:) ) PUTPR(:) = 0.0_dp
WHERE ( .NOT. GWORK6(:,:) )
    PUMF(:,:)  = 0.0_dp
    PUDR(:,:)  = 0.0_dp
    PUER(:,:)  = 0.0_dp
    PUTHL(:,:) = PTHL(:,:)
    PURW(:,:)  = PRW(:,:)
    PUPR(:,:)  = 0.0_dp
    PURC(:,:)  = 0.0_dp
    PURI(:,:)  = 0.0_dp
    PURR(:,:)  = 0.0_dp
    PURS(:,:)  = 0.0_dp
END WHERE
END SUBROUTINE CONVECT_UPDRAFT

!============================================================================

!     ##########################################################################
      SUBROUTINE CONVECT_CONDENS( KLON,                                        &
                            &  KICE, PPRES, PTHL, PRW, PRCO, PRIO, PZ, OWORK1, &
                            &  PT, PEW, PRC, PRI, PLV, PLS, PCPH   )
!     ##########################################################################

!!**** Compute temperature cloud and ice water content from enthalpy and r_w
!!
!!
!!    PURPOSE
!!    -------
!!     The purpose of this routine is to determine cloud condensate
!!     and to return values for L_v, L_s and C_ph
!!
!!
!!**  METHOD
!!    ------
!!     Condensate is extracted iteratively
!!
!!
!!    EXTERNAL
!!    --------
!!     None
!!
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!
!!      Module YOMCST
!!          RG                   ! gravity constant
!!          RALPW, RBETW, RGAMW ! constants for water saturation pressure
!!          RALPS, RBETS, RGAMS ! constants for ice saturation pressure
!!          RATM                 ! reference pressure
!!          RD, RV             ! gaz  constants for dry air and water vapor
!!          RCPD, RCPV           ! specific heat for dry air and water vapor
!!          RCW, RCS             ! specific heat for liquid water and ice
!!          RTT                  ! triple point temperature
!!          RLVTT, RLSTT         ! vaporization, sublimation heat constant
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!      Module YOE_CONVPAR
!!          XTFRZ1               ! begin of freezing interval
!!          XTFRZ2               ! end of freezing interval
!!
!!    REFERENCE
!!    ---------
!!
!!      Book1,2 of documentation ( routine CONVECT_CONDENS)
!!
!!    AUTHOR
!!    ------
!!      P. BECHTOLD       * Laboratoire d'Aerologie *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    07/11/95
!!   Last modified  04/10/97
!-------------------------------------------------------------------------------

!#include "tsmbkind.h"

!*       0.    DECLARATIONS
!              ------------

!USE YOMCST
!USE YOE_CONVPAR


IMPLICIT NONE

!*       0.1   Declarations of dummy arguments :

INTEGER, INTENT(IN)                :: KLON    ! horizontal loop index
INTEGER, INTENT(IN)                :: KICE    ! flag for ice ( 1 = yes,
                                              !                0 = no ice )
REAL(dp), DIMENSION(KLON),   INTENT(IN) :: PPRES  ! pressure
REAL(dp), DIMENSION(KLON),   INTENT(IN) :: PTHL   ! enthalpy (J/kg)
REAL(dp), DIMENSION(KLON),   INTENT(IN) :: PRW    ! total water mixing ratio
REAL(dp), DIMENSION(KLON),   INTENT(IN) :: PRCO   ! cloud water estimate (kg/kg)
REAL(dp), DIMENSION(KLON),   INTENT(IN) :: PRIO   ! cloud ice   estimate (kg/kg)
REAL(dp), DIMENSION(KLON),   INTENT(IN) :: PZ     ! level height (m)
LOGICAL, DIMENSION(KLON),INTENT(IN) :: OWORK1 ! logical mask


REAL(dp), DIMENSION(KLON),   INTENT(OUT):: PT     ! temperature
REAL(dp), DIMENSION(KLON),   INTENT(OUT):: PRC    ! cloud water mixing ratio(kg/kg)
REAL(dp), DIMENSION(KLON),   INTENT(OUT):: PRI    ! cloud ice mixing ratio  (kg/kg)
REAL(dp), DIMENSION(KLON),   INTENT(OUT):: PLV    ! latent heat L_v
REAL(dp), DIMENSION(KLON),   INTENT(OUT):: PLS    ! latent heat L_s
REAL(dp), DIMENSION(KLON),   INTENT(OUT):: PCPH   ! specific heat C_ph
REAL(dp), DIMENSION(KLON),   INTENT(OUT):: PEW    ! water saturation mixing ratio

!*       0.2   Declarations of local variables KLON

INTEGER :: JITER          ! iteration index
REAL(dp)    :: ZEPS, ZEPSA    ! R_d / R_v, 1 / ZEPS
REAL(dp)    :: ZCVOCD         ! RCPV / RCPD
REAL(dp)    :: ZRDOCP         ! R_d / C_pd

REAL(dp), DIMENSION(KLON)    :: ZEI           ! ice saturation mixing ratio
REAL(dp), DIMENSION(KLON)    :: ZWORK1, ZWORK2, ZWORK3, ZT ! work arrays


!-------------------------------------------------------------------------------

!*       1.     Initialize temperature and Exner function
!               -----------------------------------------

ZRDOCP      = RD   / RCPD
ZEPS        = RD   / RV
ZEPSA       = 1.0_dp / ZEPS
ZCVOCD      = RCPV  / RCPD


    ! Make a first temperature estimate, based e.g. on values of
    !  r_c and r_i at lower level

      !! Note that the definition of ZCPH is not the same as used in
      !! routine CONVECT_SATMIXRATIO
     PCPH(:)   = RCPD + RCPV * PRW(:)
     ZWORK1(:) = ( 1.0_dp + PRW(:) ) * RG * PZ(:)
     PT(:)     = ( PTHL(:) + PRCO(:) * RLVTT + PRIO(:) * RLSTT - ZWORK1(:) )   &
               & / PCPH(:)
     PT(:)     = MAX(180._dp, MIN( 330._dp, PT(:) ) ) ! set overflow bounds in
                                                          ! case that PTHL=0


!*       2.     Enter the iteration loop
!               ------------------------

DO JITER = 1,6
     PEW(:) = EXP( RALPW - RBETW / PT(:) - RGAMW * LOG( PT(:) ) )
     ZEI(:) = EXP( RALPS - RBETS / PT(:) - RGAMS * LOG( PT(:) ) )
     PEW(:) = ZEPS * PEW(:) / ( PPRES(:) - PEW(:) )
     ZEI(:) = ZEPS * ZEI(:) / ( PPRES(:) - ZEI(:) )

     PLV(:)    = RLVTT + ( RCPV - RCW ) * ( PT(:) - RTT ) ! compute L_v
     PLS(:)    = RLSTT + ( RCPV - RCS ) * ( PT(:) - RTT ) ! compute L_i

     ZWORK2(:) = ( PT(:) - XTFRZ2 ) / ( XTFRZ1 - XTFRZ2 ) ! freezing interval
     ZWORK2(:) = MAX( 0.0_dp, MIN(1.0_dp, ZWORK2(:) ) ) 
     IF ( KICE==0 ) ZWORK2(:) = 1.0_dp
     ZWORK2(:) = ZWORK2(:) * ZWORK2(:)
     ZWORK3(:) = ( 1.0_dp - ZWORK2(:) ) * ZEI(:) + ZWORK2(:) * PEW(:)
     PRC(:)    = MAX( 0.0_dp, ZWORK2(:) * ( PRW(:) - ZWORK3(:) ) )
     PRI(:)    = MAX( 0.0_dp, ( 1.0_dp - ZWORK2(:) ) * ( PRW(:) - ZWORK3(:) ) )
     ZT(:)     = ( PTHL(:) + PRC(:) * PLV(:) + PRI(:) * PLS(:) - ZWORK1(:) )   &
               & / PCPH(:)
     PT(:) = PT(:) + ( ZT(:) - PT(:) ) * 0.4_dp  ! force convergence
     PT(:) = MAX( 175._dp, MIN( 330._dp, PT(:) ) )
ENDDO


END SUBROUTINE CONVECT_CONDENS

!============================================================================

!     #######################################################
      SUBROUTINE CONVECT_MIXING_FUNCT( KLON,                &
                                     & PMIXC, KMF, PER, PDR )
!     #######################################################

!!**** Determine the area under the distribution function
!!     KMF = 1 : gaussian  KMF = 2 : triangular distribution function
!!
!!    PURPOSE
!!    -------
!!      The purpose of this routine is to determine the entrainment and
!!      detrainment rate by evaluating the are under the distribution
!!      function. The integration interval is limited by the critical
!!      mixed fraction PMIXC
!!
!!
!!
!!**  METHOD
!!    ------
!!      Use handbook of mathemat. functions by Abramowitz and Stegun, 1968
!!
!!
!!
!!    EXTERNAL
!!    --------
!!      None
!!
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!      None
!!
!!
!!    REFERENCE
!!    ---------
!!
!!      Book2 of documentation ( routine MIXING_FUNCT)
!!      Abramovitz and Stegun (1968), handbook of math. functions
!!
!!    AUTHOR
!!    ------
!!      P. BECHTOLD       * Laboratoire d'Aerologie *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    07/11/95
!!   Last modified  04/10/97
!-------------------------------------------------------------------------------

!#include "tsmbkind.h"

!*       0.    DECLARATIONS
!              ------------


IMPLICIT NONE

!*       0.1   Declarations of dummy arguments :

INTEGER,               INTENT(IN) :: KLON   ! horizontal dimension
INTEGER,               INTENT(IN) :: KMF    ! switch for dist. function
REAL(dp), DIMENSION(KLON), INTENT(IN) :: PMIXC  ! critical mixed fraction

REAL(dp), DIMENSION(KLON), INTENT(OUT):: PER    ! normalized entrainment rate
REAL(dp), DIMENSION(KLON), INTENT(OUT):: PDR    ! normalized detrainment rate

!*       0.2   Declarations of local variables :

REAL(dp)    :: ZSIGMA = 0.166666667_dp         ! standard deviation
REAL(dp)    :: ZFE    = 4.931813949_dp         ! integral normalization

REAL(dp)    :: ZSQRTP = 2.506628_dp ,  ZP  = 0.33267_dp   ! constants
REAL(dp)    :: ZA1    = 0.4361836_dp , ZA2 =-0.1201676_dp ! constants
REAL(dp)    :: ZA3    = 0.9372980_dp , ZT1 = 0.500498_dp  ! constants
REAL(dp)    :: ZE45   = 0.01111_dp                          ! constant

REAL(dp), DIMENSION(KLON) :: ZX, ZY, ZW1, ZW2         ! work variables
REAL(dp)    :: ZW11


!-------------------------------------------------------------------------------

!       1.     Use gaussian function for KMF=1
!              -------------------------------

IF( KMF == 1 ) THEN
    ! ZX(:)  = ( PMIXC(:) - 0.5_dp ) / ZSIGMA
      ZX(:)  = 6._dp * PMIXC(:) - 3._dp
      ZW1(:) = 1.0_dp / ( 1.0_dp+ ZP * ABS ( ZX(:) ) )
      ZY(:)  = EXP( -0.5_dp * ZX(:) * ZX(:) )
      ZW2(:) = ZA1 * ZW1(:) + ZA2 * ZW1(:) * ZW1(:) +                   &
             &   ZA3 * ZW1(:) * ZW1(:) * ZW1(:)
      ZW11   = ZA1 * ZT1 + ZA2 * ZT1 * ZT1 + ZA3 * ZT1 * ZT1 * ZT1
ENDIF

WHERE ( KMF == 1 .AND. ZX(:) >= 0.0_dp )
        PER(:) = ZSIGMA * ( 0.5_dp * ( ZSQRTP - ZE45 * ZW11              &
               & - ZY(:) * ZW2(:) ) + ZSIGMA * ( ZE45 - ZY(:) ) )        &
               & - 0.5_dp * ZE45 * PMIXC(:) * PMIXC(:)
        PDR(:) = ZSIGMA*( 0.5_dp * ( ZY(:) * ZW2(:) - ZE45 * ZW11   )    &
               & + ZSIGMA * ( ZE45 - ZY(:) ) )                           &
               & - ZE45 * ( 0.5_dp + 0.5_dp * PMIXC(:) * PMIXC(:) - PMIXC(:) )
END WHERE
WHERE ( KMF == 1 .AND. ZX(:) < 0.0_dp )
        PER(:) = ZSIGMA*( 0.5_dp * ( ZY(:) * ZW2(:) - ZE45 * ZW11   )    &
               & + ZSIGMA * ( ZE45 - ZY(:) ) )                           &
               & - 0.5_dp * ZE45 * PMIXC(:) * PMIXC(:)
        PDR(:) = ZSIGMA * ( 0.5_dp * ( ZSQRTP - ZE45 * ZW11 - ZY(:)      &
               & * ZW2(:) ) + ZSIGMA * ( ZE45 - ZY(:) ) )                &
               & - ZE45 * ( 0.5_dp + 0.5_dp * PMIXC(:) * PMIXC(:) - PMIXC(:) )
END WHERE

      PER(:) = PER(:) * ZFE
      PDR(:) = PDR(:) * ZFE


!       2.     Use triangular function KMF=2
!              -------------------------------

!     not yet released


END SUBROUTINE CONVECT_MIXING_FUNCT

!============================================================================

!     ######################################################################
      SUBROUTINE CONVECT_TSTEP_PREF( KLON, KLEV,                           &
                                   & PU, PV, PPRES, PZ, PDXDY, KLCL, KCTL, &
                                   & PTIMEA, PPREF )
!     ######################################################################

!!**** Routine to compute convective advection time step and precipitation
!!     efficiency
!!
!!
!!    PURPOSE
!!    -------
!!      The purpose of this routine is to determine the convective
!!      advection time step PTIMEC as a function of the mean ambient
!!      wind as well as the precipitation efficiency as a function
!!      of wind shear and cloud base height.
!!
!!
!!**  METHOD
!!    ------
!!
!!
!!    EXTERNAL
!!    --------
!!     None
!!
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!
!!     Module YOE_CONVPAREXT
!!          JCVEXB, JCVEXT     ! extra levels on the vertical boundaries
!!
!!    REFERENCE
!!    ---------
!!
!!      Book1,2 of documentation
!!      Fritsch and Chappell, 1980, J. Atmos. Sci.
!!      Kain and Fritsch, 1993, Meteor. Monographs, Vol.
!!
!!    AUTHOR
!!    ------
!!      P. BECHTOLD       * Laboratoire d'Aerologie *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    07/11/95
!!   Last modified  04/10/97
!-------------------------------------------------------------------------------

!#include "tsmbkind.h"

!*       0.    DECLARATIONS
!              ------------

!USE YOE_CONVPAREXT


IMPLICIT NONE

!*       0.1   Declarations of dummy arguments :

INTEGER, INTENT(IN)                    :: KLON   ! horizontal dimension
INTEGER, INTENT(IN)                    :: KLEV   ! vertical dimension
REAL(dp), DIMENSION(KLON,KLEV), INTENT(IN) :: PPRES  ! pressure (Pa)
REAL(dp), DIMENSION(KLON,KLEV), INTENT(IN) :: PU     ! grid scale horiz. wind u
REAL(dp), DIMENSION(KLON,KLEV), INTENT(IN) :: PV     ! grid scale horiz. wind v
REAL(dp), DIMENSION(KLON,KLEV), INTENT(IN) :: PZ     ! height of model layer (m)
REAL(dp), DIMENSION(KLON),      INTENT(IN) :: PDXDY  ! grid area (m^2)
INTEGER, DIMENSION(KLON),   INTENT(IN) :: KLCL   ! lifting condensation level index
INTEGER, DIMENSION(KLON),   INTENT(IN) :: KCTL   ! cloud top level index

REAL(dp), DIMENSION(KLON),      INTENT(OUT):: PTIMEA ! advective time period
REAL(dp), DIMENSION(KLON),      INTENT(OUT):: PPREF  ! precipitation efficiency


!*       0.2   Declarations of local variables KLON

INTEGER :: IIE, IKB, IKE                      ! horizontal + vertical loop bounds
INTEGER :: JI                                 ! horizontal loop index
INTEGER :: JK, JKLC, JKP5, JKCT               ! vertical loop index

INTEGER, DIMENSION(KLON)  :: IP500       ! index of 500 hPa levels
REAL(dp), DIMENSION(KLON)     :: ZCBH        ! cloud base height
REAL(dp), DIMENSION(KLON)     :: ZWORK1, ZWORK2, ZWORK3  ! work arrays


!-------------------------------------------------------------------------------

!        0.3   Set loop bounds
!              ---------------

IIE = KLON
IKB = 1 + JCVEXB
IKE = KLEV - JCVEXT


!*       1.     Determine vertical index for 500 hPa levels
!               ------------------------------------------


IP500(:) = IKB
DO JK = IKB, IKE
    WHERE ( PPRES(:,JK) >= 500.E2_dp ) IP500(:) = JK
ENDDO


!*       2.     Compute convective time step
!               ----------------------------

            ! compute wind speed at LCL, 500 hPa, CTL

DO JI = 1, IIE
   JKLC = KLCL(JI)
   JKP5 = IP500(JI)
   JKCT = KCTL(JI)
   ZWORK1(JI) = SQRT( PU(JI,JKLC) * PU(JI,JKLC) +           &
              &       PV(JI,JKLC) * PV(JI,JKLC)  )
   ZWORK2(JI) = SQRT( PU(JI,JKP5) * PU(JI,JKP5) +           &
              &       PV(JI,JKP5) * PV(JI,JKP5)  )
   ZWORK3(JI) = SQRT( PU(JI,JKCT) * PU(JI,JKCT) +           &
              &       PV(JI,JKCT) * PV(JI,JKCT)  )
ENDDO

ZWORK2(:) = MAX( 0.1_dp, 0.5_dp * ( ZWORK1(:) + ZWORK2(:) ) )

PTIMEA(:) = SQRT( PDXDY(:) ) / ZWORK2(:)


!*       3.     Compute precipitation efficiency
!               -----------------------------------

!*       3.1    Precipitation efficiency as a function of wind shear
!               ----------------------------------------------------

ZWORK2(:) = SIGN( 1._dp, ZWORK3(:) - ZWORK1(:) )
DO JI = 1, IIE
    JKLC = KLCL(JI)
    JKCT = KCTL(JI)
    ZWORK1(JI) = ( PU(JI,JKCT) - PU(JI,JKLC) )  *          &
               & ( PU(JI,JKCT) - PU(JI,JKLC) )  +          &
               & ( PV(JI,JKCT) - PV(JI,JKLC) )  *          &
               & ( PV(JI,JKCT) - PV(JI,JKLC) )
    ZWORK1(JI) = 1.E3_dp * ZWORK2(JI) * SQRT( ZWORK1(JI) ) /  &
               & MAX( 1.E-2_dp, PZ(JI,JKCT) - PZ(JI,JKLC) )
ENDDO

PPREF(:)  = 1.591_dp + ZWORK1(:) * ( -.639_dp + ZWORK1(:)       &
          &              * (  9.53E-2_dp - ZWORK1(:) * 4.96E-3_dp ) )
PPREF(:)  = MAX( .4_dp, MIN( PPREF(:), .9_dp ) )

!*       3.2    Precipitation efficiency as a function of cloud base height
!               ----------------------------------------------------------

DO JI = 1, IIE
   JKLC = KLCL(JI)
   ZCBH(JI)   = MAX( 3._dp, ( PZ(JI,JKLC) - PZ(JI,IKB) ) * 3.281E-3_dp )
ENDDO
ZWORK1(:) = .9673_dp + ZCBH(:) * ( -.7003_dp + ZCBH(:) * ( .1622_dp + &
          &   ZCBH(:) *  ( -1.2570E-2_dp + ZCBH(:) * ( 4.2772E-4_dp -   &
          &   ZCBH(:) * 5.44E-6_dp ) ) ) )
ZWORK1(:) = MAX( .4_dp, MIN( .9_dp, 1.0_dp/ ( 1.0_dp + ZWORK1(:) ) ) )

!*       3.3    Mean precipitation efficiency is used to compute rainfall
!               ----------------------------------------------------------

PPREF(:) = 0.5_dp * ( PPREF(:) + ZWORK1(:) )


END SUBROUTINE CONVECT_TSTEP_PREF

!============================================================================

!#######################################################################
 SUBROUTINE CONVECT_DOWNDRAFT( KLON, KLEV,                             &
                          & KICE, PPRES, PDPRES, PZ, PTH, PTHES,       &
                          & PRW, PRC, PRI,                             &
                          & PPREF, KLCL, KCTL, KETL,                   &
                          & PUTHL, PURW, PURC, PURI,                   &
                          & PDMF, PDER, PDDR, PDTHL, PDRW,             &
                          & PMIXF, PDTEVR, KLFS, KDBL, KML,            &
                          & PDTEVRF )
!########################################################################

!!**** Compute downdraft properties from LFS to DBL.
!!
!!
!!    PDRPOSE
!!    -------
!!      The purpose of this routine is to determine downdraft properties
!!      ( mass flux, thermodynamics )
!!
!!
!!**  METHOD
!!    ------
!!      Computations are done at every model level starting from top.
!!      The use of masks allows to optimise the inner loops (horizontal loops).
!!
!!
!!
!!    EXTERNAL
!!    --------
!!     Routine CONVECT_SATMIXRATIO
!!
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!
!!      Module YOMCST
!!          RG                 ! gravity constant
!!          RPI                ! Pi
!!          RATM               ! reference pressure
!!          RD, RV           ! gaz  constants for dry air and water vapor
!!          RCPD               ! Cpd (dry air)
!!          RCPV, RCW, RCS     ! Cp of water vapor, liquid water and ice
!!          RTT                ! triple point temperature
!!          RLVTT, RLSTT       ! vaporisation/sublimation heat at RTT
!!
!!      Module YOE_CONVPAR
!!          XCRAD              ! cloud radius
!!          XZPBL              ! thickness of downdraft detrainment layer
!!          XENTR              ! entrainment constant in pressure coordinates
!!          XRHDBC             ! relative humidity in downdraft below cloud
!!
!!     Module YOE_CONVPAREXT
!!          JCVEXB, JCVEXT     ! extra levels on the vertical boundaries
!!
!!    REFERENCE
!!    ---------
!!
!!      Book1,2 of documentation ( routine CONVECT_DOWNDRAFT)
!!      Kain and Fritsch, 1993, Meteor. Monographs, Vol.
!!
!!    AUTHOR
!!    ------
!!      P. BECHTOLD       * Laboratoire d'Aerologie *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    07/11/95
!!   Last modified  04/10/97
!-------------------------------------------------------------------------------

!#include "tsmbkind.h"

!*       0.    DECLARATIONS
!              ------------

!USE YOMCST
!USE YOE_CONVPAR
!USE YOE_CONVPAREXT


IMPLICIT NONE

!*       0.1   Declarations of dummy arguments :


INTEGER,                    INTENT(IN) :: KLON  ! horizontal dimension
INTEGER,                    INTENT(IN) :: KLEV  ! vertical dimension
INTEGER,                    INTENT(IN) :: KICE  ! flag for ice ( 1 = yes,
                                                  !                0 = no ice )
REAL(dp), DIMENSION(KLON,KLEV), INTENT(IN) :: PTH   ! grid scale theta
REAL(dp), DIMENSION(KLON,KLEV), INTENT(IN) :: PTHES ! grid scale saturated theta_e
REAL(dp), DIMENSION(KLON,KLEV), INTENT(IN) :: PRW   ! grid scale total water
                                                  ! mixing ratio
REAL(dp), DIMENSION(KLON,KLEV), INTENT(IN) :: PRC   ! grid scale r_c (cloud water)
REAL(dp), DIMENSION(KLON,KLEV), INTENT(IN) :: PRI   ! grid scale r_i (cloud ice)
REAL(dp), DIMENSION(KLON,KLEV), INTENT(IN) :: PPRES ! pressure (Pa)
REAL(dp), DIMENSION(KLON,KLEV), INTENT(IN) :: PDPRES! pressure difference between
                                                  ! bottom and top of layer (Pa)
REAL(dp), DIMENSION(KLON,KLEV), INTENT(IN) :: PZ    ! level height (m)
INTEGER, DIMENSION(KLON),   INTENT(IN) :: KLCL  ! contains vert. index of LCL
INTEGER, DIMENSION(KLON),   INTENT(IN) :: KCTL  ! contains vert. index of CTL
INTEGER, DIMENSION(KLON),   INTENT(IN) :: KETL  ! contains vert. index of
                                                  ! equilibrium (zero buoyancy) level
INTEGER, DIMENSION(KLON),   INTENT(IN) :: KML   ! " vert. index of melting level
REAL(dp), DIMENSION(KLON,KLEV), INTENT(IN) :: PUTHL ! updraft enthalpy (J/kg)
REAL(dp), DIMENSION(KLON,KLEV), INTENT(IN) :: PURW  ! updraft total water (kg/kg)
REAL(dp), DIMENSION(KLON,KLEV), INTENT(IN) :: PURC  ! updraft r_c (kg/kg)
REAL(dp), DIMENSION(KLON,KLEV), INTENT(IN) :: PURI  ! updraft r_i (kg/kg)
REAL(dp), DIMENSION(KLON),      INTENT(IN) :: PPREF ! precipitation efficiency


REAL(dp), DIMENSION(KLON,KLEV), INTENT(OUT):: PDMF   ! downdraft mass flux (kg/s)
REAL(dp), DIMENSION(KLON,KLEV), INTENT(OUT):: PDER   ! downdraft entrainment (kg/s)
REAL(dp), DIMENSION(KLON,KLEV), INTENT(OUT):: PDDR   ! downdraft detrainment (kg/s)
REAL(dp), DIMENSION(KLON,KLEV), INTENT(OUT):: PDTHL  ! downdraft enthalpy (J/kg)
REAL(dp), DIMENSION(KLON,KLEV), INTENT(OUT):: PDRW   ! downdraft total water (kg/kg)
REAL(dp), DIMENSION(KLON),      INTENT(OUT):: PMIXF  ! mixed fraction at LFS
REAL(dp), DIMENSION(KLON),      INTENT(OUT):: PDTEVR ! total downdraft evaporation
                                                   ! rate at LFS (kg/s)
REAL(dp), DIMENSION(KLON,KLEV), INTENT(OUT):: PDTEVRF! downdraft evaporation rate
INTEGER, DIMENSION(KLON),  INTENT(OUT):: KLFS    ! contains vert. index of LFS
INTEGER, DIMENSION(KLON),  INTENT(OUT):: KDBL    ! contains vert. index of DBL

!*       0.2   Declarations of local variables :

INTEGER :: IIE, IKB, IKE     ! horizontal + vertical loop bounds
INTEGER :: JK, JKP, JKM, JKT ! vertical loop index
INTEGER :: JI, JL            ! horizontal loop index
INTEGER :: JITER             ! iteration loop index
REAL(dp)    :: ZCPORD, ZRDOCP    ! C_pd / R_d, R_d / C_pd
REAL(dp)    :: ZEPS              ! R_d / R_v
REAL(dp)    :: ZEPSA, ZCVOCD     ! R_v / R_d, C_pv / C_pd

INTEGER, DIMENSION(KLON) :: IDDT      ! top level of detrainm. layer
REAL(dp), DIMENSION(KLON)    :: ZTHE      ! environm. theta_e (K)
REAL(dp), DIMENSION(KLON)    :: ZDT, ZDTP ! downdraft temperature (K)
REAL(dp), DIMENSION(KLON)    :: ZCPH      ! specific heat C_ph
REAL(dp), DIMENSION(KLON)    :: ZLV, ZLS  ! latent heat of vaporis., sublim.
REAL(dp), DIMENSION(KLON)    :: ZDDT      ! thickness (hPa) of detrainm. layer
REAL(dp), DIMENSION(KLON)    :: ZPI       ! Pi=(P0/P)**(Rd/Cpd)
REAL(dp), DIMENSION(KLON)    :: ZWORK1, ZWORK2, ZWORK3  ! work arrays
REAL(dp), DIMENSION(KLON)    :: ZWORK4, ZWORK5          ! work arrays
LOGICAL, DIMENSION(KLON)   :: GWORK1                  ! work array


!-------------------------------------------------------------------------------

!        0.3    Set loop bounds
!               ---------------

IIE = KLON
IKB = 1 + JCVEXB
IKE = KLEV - JCVEXT


!*       1.     Initialize downdraft properties
!               -------------------------------

ZCPORD     = RCPD / RD
ZRDOCP     = RD  / RCPD
ZEPS       = RD  / RV
ZEPSA      = RV  / RD
ZCVOCD     = RCPV / RCPD

PDMF(:,:)  = 0.0_dp
PDER(:,:)  = 0.0_dp
PDDR(:,:)  = 0.0_dp
PDRW(:,:)  = 0.0_dp
PDTHL(:,:) = 0.0_dp
PDTEVR(:)  = 0.0_dp
PMIXF(:)   = 0.0_dp
ZTHE(:)    = 0.0_dp
ZDDT(:)    = PDPRES(:,IKB+2)
KDBL(:)    = IKB + 1
KLFS(:)    = IKB + 1
IDDT(:)    = KDBL(:) + 1


!*       2.     Determine the LFS by looking for minimum of environmental
!               saturated theta_e
!               ----------------------------------------------------------

ZWORK1(:) = 900._dp   ! starting value for search of minimum envir. theta_e
DO JK = MINVAL( KLCL(:) ) + 2, MAXVAL( KETL(:) )
   DO JI = 1, IIE
      GWORK1(JI) = JK >= KLCL(JI) + 2 .AND. JK < KETL(JI)
      IF ( GWORK1(JI) .AND. ZWORK1(JI) > PTHES(JI,JK) ) THEN
         KLFS(JI)   = JK
         ZWORK1(JI) = MIN( ZWORK1(JI), PTHES(JI,JK) )
      ENDIF
   ENDDO
ENDDO


!*       3.     Determine the mixed fraction using environmental and updraft
!               values of theta_e at LFS
!               ---------------------------------------------------------

DO JI = 1, IIE
    JK = KLFS(JI)
    ZPI(JI)    = ( RATM / PPRES(JI,JK) ) ** ZRDOCP
      ! compute updraft theta_e
    ZWORK3(JI) = PURW(JI,JK) - PURC(JI,JK) - PURI(JI,JK)
    ZDT(JI)    = PTH(JI,JK) / ZPI(JI)
    ZLV(JI)    = RLVTT + ( RCPV - RCW ) * ( ZDT(JI) - RTT )
    ZLS(JI)    = RLSTT + ( RCPV - RCS ) * ( ZDT(JI) - RTT )
    ZCPH(JI)   = RCPD + RCPV * PURW(JI,JK)
    ZDT(JI)    = ( PUTHL(JI,JK) - ( 1.0_dp + PURW(JI,JK) ) * RG * PZ(JI,JK) &
               & + ZLV(JI) * PURC(JI,JK) + ZLS(JI) * PURI(JI,JK) ) / ZCPH(JI)
    ZWORK1(JI) = ZDT(JI) * ZPI(JI) ** ( 1.0_dp - 0.28_dp * ZWORK3(JI) )   &
               & * EXP( ( 3374.6525_dp / ZDT(JI) - 2.5403_dp )         &
               & * ZWORK3(JI) * ( 1.0_dp + 0.81_dp * ZWORK3(JI) ) )
      ! compute environmental theta_e
    ZDT(JI)    = PTH(JI,JK) / ZPI(JI)
    ZLV(JI)    = RLVTT + ( RCPV - RCW ) * ( ZDT(JI) - RTT )
    ZLS(JI)    = RLSTT + ( RCPV - RCS ) * ( ZDT(JI) - RTT )
    ZWORK3(JI) = PRW(JI,JK) - PRC(JI,JK) - PRI(JI,JK)
    ZCPH(JI)   = RCPD + RCPV * PRW(JI,JK)
    ZWORK2(JI) = ZDT(JI) * ZPI(JI) ** ( 1.0_dp - 0.28_dp * ZWORK3(JI) ) &
               & * EXP( ( 3374.6525_dp / ZDT(JI) - 2.5403_dp )       &
               & * ZWORK3(JI) * ( 1.0_dp + 0.81_dp * ZWORK3(JI) ) )
      ! compute mixed fraction
    PMIXF(JI)  = MAX( 0.0_dp, ( ZWORK1(JI) - PTHES(JI,JK) ) )            &
               &  / ( ZWORK1(JI) - ZWORK2(JI) + 1.E-10_dp )
    PMIXF(JI)  = MAX(0.0_dp, MIN( 1.0_dp, PMIXF(JI) ) )
    ZWORK4(JI) = PPRES(JI,JK)
ENDDO


!*       4.     Estimate the effect of melting on the downdraft
!               ---------------------------------------------

ZWORK1(:) = 0.0_dp
      ! use total solid precipitation
!DO JK = IKB + 1, IKE
!    ZWORK1(:) = ZWORK1(:) + PURS(:,JK) ! total snow/hail content
!END DO

DO JI = 1, IIE
     JK  = KLCL(JI)
     JKP = KCTL(JI)
     ZWORK1(JI) = 0.5_dp * ( PURW(JI,JK) - PURW(JI,JKP) )
ENDDO

      ! temperature perturbation due to melting at LFS
ZWORK3(:) = 0.0_dp
WHERE( KML(:) > IKB + 2 )
          ZWORK3(:) = ZWORK1(:) * ( ZLS(:) - ZLV(:) ) / ZCPH(:)
          ZDT(:)    = ZDT(:) - ZWORK3(:) * REAL(KICE)
END WHERE


!*       5.     Initialize humidity at LFS as a saturated mixture of
!               updraft and environmental air
!               -----------------------------------------------------

DO JI = 1, IIE
     JK = KLFS(JI)
     PDRW(JI,JK)  = PMIXF(JI) * PRW(JI,JK) + ( 1.0_dp - PMIXF(JI) ) * PURW(JI,JK)
     ZWORK2(JI)   = PDRW(JI,JK) - ( 1.0_dp - PMIXF(JI) )                          &
                  &                  * ( PURC(JI,JK) + PURI(JI,JK) )
ENDDO


!*       6.1    Determine the DBL by looking for level where the envir.
!               theta_es at the LFS corrected by melting effects  becomes
!               larger than envir. value
!               ---------------------------------------------------------

! compute satur. mixing ratio for melting corrected temperature
    CALL CONVECT_SATMIXRATIO( KLON, ZWORK4, ZDT, ZWORK3, ZLV, ZLS, ZCPH )

      ! compute envir. saturated theta_e for melting corrected temperature
    ZWORK1(:) = MIN( ZWORK2(:), ZWORK3(:) )
    ZWORK3(:) = ZWORK3(:) * ZWORK4(:) / ( ZWORK3(:) + ZEPS ) ! sat. pressure
    ZWORK3(:) = LOG( ZWORK3(:) / 613.3_dp )

      ! dewp point temperature
    ZWORK3(:) = ( 4780.8_dp - 32.19_dp * ZWORK3(:) ) / ( 17.502_dp - ZWORK3(:) )

      ! adiabatic saturation temperature
    ZWORK3(:) = ZWORK3(:) - ( .212_dp + 1.571E-3_dp * ( ZWORK3(:) - RTT )          &
              &   - 4.36E-4_dp * ( ZDT(:) - RTT ) ) * ( ZDT(:) - ZWORK3(:) )
    ZWORK4(:) = SIGN(0.5_dp, ZWORK2(:) - ZWORK3(:) )
    ZDT(:)    = ZDT(:) * ( 0.5_dp + ZWORK4(:) ) + ( 0.5_dp - ZWORK4(:) ) * ZWORK3(:)
    ZWORK2(:) = ZDT(:) * ZPI(:) ** ( 1.0_dp - 0.28_dp * ZWORK2(:) )                   &
              &                   * EXP( ( 3374.6525_dp / ZDT(:) - 2.5403_dp )     &
              &                   * ZWORK1(:) * ( 1.0_dp + 0.81_dp * ZWORK1(:) ) )

GWORK1(:) = .TRUE.
JKM = MAXVAL( KLFS(:) )
DO JK = JKM - 1, IKB + 1, -1
  DO JI = 1, IIE
     IF ( JK < KLFS(JI) .AND. ZWORK2(JI) > PTHES(JI,JK) .AND. GWORK1(JI) ) THEN
          KDBL(JI) = JK
          GWORK1(JI) = .FALSE.
     ENDIF
  ENDDO
ENDDO


!*       7.     Define mass flux and entr/detr. rates at LFS
!               -------------------------------------------

DO JI = 1, IIE
     JK = KLFS(JI)
     ZWORK1(JI)  = PPRES(JI,JK) /                                            &
                 & ( RD * ZDT(JI) * ( 1.0_dp + ZEPS * ZWORK1(JI) ) ) ! density
     PDMF(JI,JK) = - ( 1.0_dp - PPREF(JI) ) * ZWORK1(JI) * RPI * XCRAD * XCRAD
     PDTHL(JI,JK)= ZWORK2(JI)   ! theta_l is here actually theta_e
     ZWORK2(JI)  = PDMF(JI,JK)
     PDDR(JI,JK) = 0.0_dp
     PDER(JI,JK) = - PMIXF(JI) * PDMF(JI,JK)
ENDDO


!         7.1   Downdraft detrainment is assumed to occur in a layer
!               of 60 hPa, determine top level IDDT of this layer
!               ---------------------------------------------------------

ZWORK1(:) = 0.0_dp
DO JK = IKB + 2, JKM
      ZWORK1(:) = ZWORK1(:) + PDPRES(:,JK)
      WHERE ( JK > KDBL(:) .AND. ZWORK1(:) <= XZPBL )
           ZDDT(:) = ZWORK1(:)
           IDDT(:) = JK
      END WHERE
ENDDO


!*       8.     Enter loop for downdraft computations. Make a first guess
!               of initial downdraft mass flux.
!               In the downdraft computations we use theta_es instead of
!               enthalpy as it allows to better take into account evaporation
!               effects. As the downdraft detrainment rate is zero apart
!               from the detrainment layer, we just compute enthalpy
!               downdraft from theta_es in this layer.
!               ----------------------------------------------------------


ZWORK5(:) = 0.0_dp

DO JK =  JKM - 1, IKB + 1, -1
  JKP = JK + 1
  DO JI = 1, IIE
    IF ( JK < KLFS(JI) .AND. JK >= IDDT(JI) )  THEN
      PDER(JI,JK)  = - ZWORK2(JI) * XENTR * PDPRES(JI,JKP) / XCRAD
                                               ! DER and DPRES are positive
      PDMF(JI,JK)  = PDMF(JI,JKP) - PDER(JI,JK)
      ZPI(JI)      = ( RATM / PPRES(JI,JK) ) ** ZRDOCP
      ZDT(JI)      = PTH(JI,JK) / ZPI(JI)
      ZWORK1(JI)   = PRW(JI,JK) - PRC(JI,JK) - PRI(JI,JK)
      ZTHE(JI)     = ZDT(JI) * ZPI(JI) ** ( 1.0_dp - 0.28_dp * ZWORK1(JI) )       &
                   &           * EXP( ( 3374.6525_dp / ZDT(JI) - 2.5403_dp )   &
                   &           * ZWORK1(JI) * ( 1.0_dp + 0.81_dp * ZWORK1(JI) ) )
         ! PDTHL is here theta_es, later on in this routine this table is
         ! reskipped to enthalpy
      PDTHL(JI,JK) = ( PDTHL(JI,JKP) * PDMF(JI,JKP) - ZTHE(JI) * PDER(JI,JK) ) /   &
                   & ( PDMF(JI,JK) - 1.E-7_dp )
      PDRW(JI,JK)  = ( PDRW(JI,JKP) * PDMF(JI,JKP) - PRW(JI,JK) * PDER(JI,JK) ) /  &
                   & ( PDMF(JI,JK) - 1.E-7_dp )
    ENDIF
    IF ( JK < IDDT(JI) .AND. JK >= KDBL(JI) )   THEN
      JL = IDDT(JI)
      PDDR(JI,JK)  = - PDMF(JI,JL) * PDPRES(JI,JKP) / ZDDT(JI)
      PDMF(JI,JK)  = PDMF(JI,JKP) + PDDR(JI,JK)
      PDTHL(JI,JK) = PDTHL(JI,JKP)
      PDRW(JI,JK)  = PDRW(JI,JKP)
    ENDIF
  ENDDO
ENDDO


!*       9.     Calculate total downdraft evaporation
!               rate for given mass flux (between DBL and IDDT)
!               -----------------------------------------------

PDTEVRF(:,:) = 0.0_dp

JKT = MAXVAL( IDDT(:) )
DO JK = IKB + 1, JKT

       ZPI(:) = ( RATM / PPRES(:,JK) ) ** ZRDOCP
       ZDT(:) = PTH(:,JK) / ZPI(:)

!*       9.1    Determine wet bulb temperature at DBL from theta_e.
!               The iteration algoritm is similar to that used in
!               routine CONVECT_CONDENS
!               --------------------------------------------------

   DO JITER = 1, 4
       CALL CONVECT_SATMIXRATIO( KLON, PPRES(:,JK), ZDT, ZWORK1, ZLV, ZLS, ZCPH )
    DO JI = 1, IIE
     IF ( JK <= IDDT(JI) ) THEN
       ZDTP(JI) = PDTHL(JI,JK) / ( ZPI(JI) ** ( 1. - 0.28 * ZWORK1(JI) )         &
                      * EXP( ( 3374.6525 / ZDT(JI) - 2.5403 )                 &
                             * ZWORK1(JI) * ( 1. + 0.81 * ZWORK1(JI) ) ) )
       ZDT(JI)  = 0.4 * ZDTP(JI) + 0.6 * ZDT(JI) ! force convergence
     END IF
    ENDDO
   ENDDO


!*       9.2    Sum total downdraft evaporation rate. No evaporation
!               if actual humidity is larger than specified one.
!               -----------------------------------------------------

   ZWORK2(:) = ZWORK1(:) / ZDT(:) * ( RBETW / ZDT(:) - RGAMW ) ! dr_sat/dT
   ZWORK2(:) = ZLV(:) / ZCPH(:) * ZWORK1(:) * ( 1.0_dp - XRHDBC ) /              &
             &      ( 1.0_dp + ZLV(:) / ZCPH(:) * ZWORK2(:) ) ! temperature perturb
                                                             ! due to evaporation
   ZDT(:)    = ZDT(:) + ZWORK2(:)

   CALL CONVECT_SATMIXRATIO( KLON, PPRES(:,JK), ZDT, ZWORK3, ZLV, ZLS, ZCPH )
!
   DO JI = 1, IIE
    IF ( JK <= IDDT(JI) ) THEN
     ZWORK3(JI)    = ZWORK3(JI) * XRHDBC
     ZWORK1(JI)    = MAX( 0.0_dp, ZWORK3(JI) - PDRW(JI,JK) )
     PDTEVR(JI)    = PDTEVR(JI) + ZWORK1(JI) * PDDR(JI,JK)
     PDTEVRF(JI,JK)= PDTEVRF(JI,JK+1) + ZWORK1(JI) * PDDR(JI,JK)
      ! compute enthalpie and humidity in the detrainment layer
     PDRW(JI,JK)   = MAX( PDRW(JI,JK), ZWORK3(JI) )
     PDTHL(JI,JK)  = ( ( RCPD + PDRW(JI,JK) * RCPV ) * ZDT(JI)                    &
                    + ( 1. + PDRW(JI,JK) ) * RG * PZ(JI,JK) )
    ENDIF
   ENDDO
!
ENDDO


!*      12.     If downdraft does not evaporate any water for specified
!               relative humidity, no downdraft is allowed
!               ---------------------------------------------------------
DO JK = IKB, JKM
 DO JI = 1, IIE
  IF (PDTEVR(JI) < 1. .OR. KLFS(JI) == IKB + 1 ) THEN
      KDBL(JI)     = IKB
      KLFS(JI)     = IKB
      PDMF(JI,JK)  = 0.0_dp
      PDER(JI,JK)  = 0.0_dp
      PDDR(JI,JK)  = 0.0_dp
      PDTHL(JI,JK) = 0.0_dp
      PDRW(JI,JK)  = 0.0_dp
      PDTEVR(JI)   = 0.0_dp
      PDTEVRF(JI,JK)= 0.0_dp
  ENDIF
 ENDDO

ENDDO

END SUBROUTINE CONVECT_DOWNDRAFT

!============================================================================

!     ######################################################################
      SUBROUTINE CONVECT_PRECIP_ADJUST( KLON, KLEV,                        &
                                      & PPRES, PUMF, PUER, PUDR,           &
                                      & PUPR, PUTPR, PURW,                 &
                                      & PDMF, PDER, PDDR, PDTHL, PDRW,     &
                                      & PPREF, PTPR, PMIXF, PDTEVR,        &
                                      & KLFS, KDBL, KLCL, KCTL, KETL,      &
                                      & PDTEVRF )
!     ######################################################################

!!**** Adjust up- and downdraft mass fluxes to be consistent with the
!!     mass transport at the LFS given by the precipitation efficiency
!!     relation.
!!
!!
!!    PURPOSE
!!    -------
!!      The purpose of this routine is to adjust up- and downdraft mass
!!      fluxes below the LFS to be consistent with the precipitation
!!      efficiency relation
!!
!!
!!
!!**  METHOD
!!    ------
!!
!!
!!    EXTERNAL
!!    --------
!!     None
!!
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!
!!     Module YOE_CONVPAREXT
!!          JCVEXB, JCVEXT     ! extra levels on the vertical boundaries
!!
!!     Module YOE_CONVPAR
!!        XUSRDPTH             ! pressure depth to compute updraft humidity
!!                             ! supply rate for downdraft
!!
!!    REFERENCE
!!    ---------
!!
!!      Book1,2 of documentation ( routine CONVECT_PRECIP_ADJUST)
!!
!!    AUTHOR
!!    ------
!!      P. BECHTOLD       * Laboratoire d'Aerologie *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    07/11/95
!!   Last modified  04/10/97
!-------------------------------------------------------------------------------

!#include "tsmbkind.h"

!*       0.    DECLARATIONS
!              ------------

!USE YOE_CONVPAREXT
!USE YOE_CONVPAR

IMPLICIT NONE

!*       0.1   Declarations of dummy arguments :


INTEGER,                    INTENT(IN) :: KLON  ! horizontal dimension
INTEGER,                    INTENT(IN) :: KLEV  ! vertical dimension
REAL(dp), DIMENSION(KLON,KLEV), INTENT(IN) :: PPRES ! pressure (Pa)
REAL(dp), DIMENSION(KLON,KLEV), INTENT(IN) :: PURW  ! updraft total water (kg/kg)
REAL(dp), DIMENSION(KLON),      INTENT(IN) :: PUTPR ! updraft  total precipit. (kg/s
REAL(dp), DIMENSION(KLON),      INTENT(IN) :: PPREF ! precipitation efficiency
REAL(dp), DIMENSION(KLON),      INTENT(IN) :: PMIXF ! critical mixed fraction at LCL
INTEGER, DIMENSION(KLON),   INTENT(IN) :: KLCL  ! contains vert. index of LCL
INTEGER, DIMENSION(KLON),   INTENT(IN) :: KCTL  ! contains vert. index of CTL
INTEGER, DIMENSION(KLON),   INTENT(IN) :: KETL  ! contains vert. index of equilibrium
                                                  ! (zero buoyancy) level
INTEGER, DIMENSION(KLON),  INTENT(INOUT) :: KLFS ! contains vert. index of LFS
INTEGER, DIMENSION(KLON),  INTENT(INOUT) :: KDBL ! contains vert. index of DBL

REAL(dp), DIMENSION(KLON),      INTENT(INOUT) :: PDTEVR ! total downdraft evaporation
                                                      ! rate at LFS
REAL(dp), DIMENSION(KLON,KLEV), INTENT(INOUT) :: PDTEVRF! downdraft evaporation rate
REAL(dp), DIMENSION(KLON,KLEV), INTENT(INOUT) :: PUMF   ! updraft mass flux (kg/s)
REAL(dp), DIMENSION(KLON,KLEV), INTENT(INOUT) :: PUER   ! updraft entrainment (kg/s)
REAL(dp), DIMENSION(KLON,KLEV), INTENT(INOUT) :: PUDR   ! updraft detrainment (kg/s)
REAL(dp), DIMENSION(KLON,KLEV), INTENT(INOUT) :: PUPR   ! updraft  precipit. (kg/s)
REAL(dp), DIMENSION(KLON,KLEV), INTENT(INOUT) :: PDMF   ! downdraft mass flux (kg/s)
REAL(dp), DIMENSION(KLON,KLEV), INTENT(INOUT) :: PDER   ! downdraft entrainment (kg/s)
REAL(dp), DIMENSION(KLON,KLEV), INTENT(INOUT) :: PDDR   ! downdraft detrainment (kg/s)
REAL(dp), DIMENSION(KLON,KLEV), INTENT(INOUT) :: PDTHL  ! downdraft enthalpy (J/kg)
REAL(dp), DIMENSION(KLON,KLEV), INTENT(INOUT) :: PDRW   ! downdraft total water (kg/kg)

REAL(dp), DIMENSION(KLON),     INTENT(OUT)   :: PTPR    ! total precipitation (kg/s)
                                                      ! = downdraft precipitation

!*       0.2   Declarations of local variables :

INTEGER :: IIE, IKB, IKE        ! horizontal + vertical loop bounds
INTEGER :: JK, JKT1, JKT2, JKT3 ! vertical loop index
INTEGER :: JI                   ! horizontal loop index

INTEGER, DIMENSION(KLON) :: IPRL
REAL(dp), DIMENSION(KLON)    :: ZWORK1, ZWORK2, ZWORK3, &
                           &  ZWORK4, ZWORK5, ZWORK6    ! work arrays


!-------------------------------------------------------------------------------

!        0.3   Set loop bounds
!              ---------------

IKB  = 1 + JCVEXB
IKE  = KLEV - JCVEXT
IIE  = KLON
JKT1 = MAXVAL( KLFS(:) )
JKT2 = MAXVAL( KCTL(:) )
JKT3 = MINVAL( KLCL(:) )


!        1.    Set some output variables for columns where no downdraft
!              exists. Exit if there is no downdraft at all.
!              --------------------------------------------------------

IPRL(:) = IKB
PTPR(:) = 0.0_dp

WHERE ( PDTEVR(:) == 0.0_dp )
     PTPR(:)    = PUTPR(:)  ! no downdraft evaporation => no downdraft, all
                            ! precipitation occurs in updraft
END WHERE
IF ( COUNT( PDTEVR(:) > 0.0_dp ) == 0 ) RETURN ! exit routine if no downdraft exists

!*       2.     The total mass transported from the updraft to the down-
!               draft at the LFS must be consistent with the three water
!               budget terms :
!               ---------------------------------------------------------

!*       2.1    Downdraft evaporation rate at the DBL. The evaporation
!               rate in downdraft must be consistent with precipitation
!               efficiency relation.
!               --------------------------------------------------------


DO JI = 1, IIE
     JK = KLFS(JI)
     ZWORK1(JI) = PDTEVR(JI) / MIN( -1.E-1_dp, PDMF(JI,JK) )
     ZWORK6(JI) = PDMF(JI,JK)
ENDDO

!*       2.2    Some preliminar computations for downdraft = total
!               precipitation rate. The precipitation is evaluated in
!               a layer thickness DP=XUSRDPTH=165 hPa above the LCL.
!               The difference between updraft precipitation and downdraft
!               precipitation (updraft supply rate) is used to drive the
!               downdraft through evaporational cooling.
!               --------------------------------------------------------

DO JI = 1, IIE
     JK = KLCL(JI)
     ZWORK5(JI) = PPRES(JI,JK)
ENDDO

PTPR(:) = 0.0_dp
DO JK = JKT3, JKT2
    WHERE ( JK >= KLCL(:) .AND. PPRES(:,JK) >= ZWORK5(:) - XUSRDPTH )
        PTPR(:) = PTPR(:) + PUPR(:,JK)
        IPRL(:) = JK
    END WHERE
ENDDO
IPRL(:) = MIN( KETL(:), IPRL(:) )

DO JI = 1, IIE
     JK = IPRL(JI)
     PTPR(JI) = PUMF(JI,JK+1) * PURW(JI,JK+1) + PTPR(JI)
ENDDO

PTPR(:)   = PPREF(:) * MIN( PUTPR(:), PTPR(:) )
ZWORK4(:) = PUTPR(:) - PTPR(:)


!*       2.3    Total amount of precipitation that falls out of the up-
!               draft between the LCL and the LFS.
!               Condensate transfer from up to downdraft at LFS
!               ---------------------------------------------------------

ZWORK5(:) = 0.0_dp
DO JK = JKT3, JKT1
     WHERE ( JK >= KLCL(:) .AND. JK <= KLFS(:) )
           ZWORK5(:) = ZWORK5(:) +  PUPR(:,JK)
     END WHERE
ENDDO

DO JI = 1, IIE
     JK = KLFS(JI)
     ZWORK2(JI) = ( 1.0_dp - PPREF(JI) ) * ZWORK5(JI) *                     &
                & ( 1.0_dp - PMIXF(JI) ) / MAX( 1.E-1_dp, PUMF(JI,JK) )
ENDDO


!*       2.4    Increase the first guess downdraft mass flux to satisfy
!               precipitation efficiency relation.
!               If downdraft does not evaporate any water at the DBL for
!               the specified relative humidity, or if the corrected mass
!               flux at the LFS is positive no downdraft is allowed
!               ---------------------------------------------------------


!ZWORK1(:) = ZWORK4(:) / ( ZWORK1(:) + ZWORK2(:) + 1.E-8_dp )
!ZWORK2(:) = ZWORK1(:) / MIN( -1.E-1_dp, ZWORK6(:) ) ! ratio of budget consistent to actual DMF
ZWORK1(:) = -ZWORK4(:) / ( -ZWORK1(:) + ZWORK2(:) + 1.E-8_dp )
ZWORK2(:) = ZWORK1(:) / MIN( -1.E-1_dp, ZWORK6(:) ) ! ratio of budget consistent to actual DMF


ZWORK3(:) = 1.0_dp
ZWORK6(:) = 1.0_dp
WHERE ( ZWORK1(:) > 0.0_dp .OR. PDTEVR(:) < 1.0_dp )
   KDBL(:)   = IKB
   KLFS(:)   = IKB
   PDTEVR(:) = 0.0_dp
   ZWORK2(:) = 0.0_dp
   ZWORK3(:) = 0.0_dp
   ZWORK6(:) = 0.0_dp
END WHERE

DO JK = IKB, JKT1
     PDMF(:,JK)  = PDMF(:,JK)  * ZWORK2(:)
     PDER(:,JK)  = PDER(:,JK)  * ZWORK2(:)
     PDDR(:,JK)  = PDDR(:,JK)  * ZWORK2(:)
   PDTEVRF(:,JK) = PDTEVRF(:,JK)* ZWORK2(:)
     PDRW(:,JK)  = PDRW(:,JK)  * ZWORK3(:)
     PDTHL(:,JK) = PDTHL(:,JK) * ZWORK3(:)
ENDDO
ZWORK4(:) = ZWORK2(:)


!*       3.     Increase updraft mass flux, mass detrainment rate, and water
!               substance detrainment rates to be consistent with the transfer
!               of the estimated mass from the up- to the downdraft at the LFS
!               --------------------------------------------------------------

DO JI = 1, IIE
    JK = KLFS(JI)
    ZWORK2(JI) = ( 1.0_dp - ZWORK6(JI) ) + ZWORK6(JI) *                   &
               &  ( PUMF(JI,JK) - ( 1.0_dp - PMIXF(JI) ) * ZWORK1(JI) ) / &
               &  MAX( 1.E-1_dp, PUMF(JI,JK) )
ENDDO


JKT1  = MAXVAL( KLFS(:) )  ! value of KLFS might have been reset to IKB above
DO JK = IKB, JKT1
    DO JI = 1, IIE
      IF ( JK <= KLFS(JI) ) THEN
        PUMF(JI,JK)  = PUMF(JI,JK)  * ZWORK2(JI)
        PUER(JI,JK)  = PUER(JI,JK)  * ZWORK2(JI)
        PUDR(JI,JK)  = PUDR(JI,JK)  * ZWORK2(JI)
        PUPR(JI,JK)  = PUPR(JI,JK)  * ZWORK2(JI)
      ENDIF
    ENDDO
ENDDO
!DO JI = 1, IIE
!   JK = KLFS(JI)+1
!   PUDR(JI,JK)=PUMF(JI,JK-1)-PUMF(JI,JK)+PUER(JI,JK)
!END DO



!*       4.     Increase total = downdraft precipitation and evaporation rate
!               -------------------------------------------------------------

WHERE ( PDTEVR(:) > 0.0_dp )
    PTPR(:)    = PTPR(:) + PPREF(:) * ZWORK5(:) * ( ZWORK2(:) - 1.0_dp )
    PDTEVR(:)  = PUTPR(:) - PTPR(:)
    PDTEVRF(:,IKB+1)  = PDTEVR(:)
ELSEWHERE
    PTPR(:)    = PUTPR(:)
END WHERE


END SUBROUTINE CONVECT_PRECIP_ADJUST

!============================================================================

!####################################################################
 SUBROUTINE CONVECT_CLOSURE( KLON, KLEV,                            &
                           & PPRES, PDPRES, PZ, PDXDY, PLMASS,      &
                           & PTHL, PTH, PRW, PRC, PRI, OTRIG1,      &
                           & PTHC, PRWC, PRCC, PRIC, PWSUB,         &
                           & KLCL, KDPL, KPBL, KLFS, KCTL, KML,     &
                           & PUMF, PUER, PUDR, PUTHL, PURW,         &
                           & PURC, PURI, PUPR,                      &
                           & PDMF, PDER, PDDR, PDTHL, PDRW,         &
                           & PTPR, PSPR, PDTEVR,                    &
                           & PCAPE, PTIMEC,                         &
                           & KFTSTEPS,                              &
                           & PDTEVRF, PPRLFLX, PPRSFLX              )
!####################################################################

!!**** Uses modified Fritsch-Chappell closure
!!
!!
!!    PURPOSE
!!    -------
!!      The purpose of this routine is to determine the final adjusted
!!     (over a time step PTIMEC) environmental values of THETA_l, R_w, R_c, R_i
!!      The final convective tendencies can then be evaluated in the main
!!      routine DEEP_CONVECT by (PTHC-PTH)/PTIMEC
!!
!!
!!**  METHOD
!!    ------
!!      Computations are done at every model level starting from bottom.
!!      The use of masks allows to optimise the inner loops (horizontal loops).
!!
!!
!!
!!    EXTERNAL
!!    --------
!!
!!    CONVECT_CLOSURE_THRVLCL
!!    CONVECT_CLOSURE_ADJUST
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!      Module YOMCST
!!          RG                 ! gravity constant
!!          RATM               ! reference pressure
!!          RD, RV           ! gaz  constants for dry air and water vapor
!!          RCPD, RCPV         ! specific heat for dry air and water vapor
!!          RCW, RCS           ! specific heat for liquid water and ice
!!          RTT                ! triple point temperature
!!          RLVTT, RLSTT       ! vaporization, sublimation heat constant
!!
!!      Module YOE_CONVPAR
!!          XA25               ! reference grid area
!!          XSTABT             ! stability factor in time integration
!!          XSTABC             ! stability factor in CAPE adjustment
!!          XMELDPTH           ! allow melting over specific pressure depth
!!
!!     Module YOE_CONVPAREXT
!!          JCVEXB, JCVEXT     ! extra levels on the vertical boundaries
!!
!!
!!    REFERENCE
!!    ---------
!!
!!      Book1,2 of documentation ( routine CONVECT_CLOSURE)
!!      Fritsch and Chappell, 1980, J. Atmos. Sci.
!!      Kain and Fritsch, 1993, Meteor. Monographs, Vol.
!!
!!    AUTHOR
!!    ------
!!      P. BECHTOLD       * Laboratoire d'Aerologie *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    26/03/96
!!   Peter Bechtold 04/10/97 change for enthalpie, r_c + r_i tendencies
!-------------------------------------------------------------------------------

!#include "tsmbkind.h"

!*       0.    DECLARATIONS
!              ------------

!USE YOMCST
!USE YOE_CONVPAR
!USE YOE_CONVPAREXT


IMPLICIT NONE

!*       0.1   Declarations of dummy arguments :

INTEGER,                   INTENT(IN) :: KLON   ! horizontal dimension
INTEGER,                   INTENT(IN) :: KLEV   ! vertical dimension
INTEGER, DIMENSION(KLON),  INTENT(IN) :: KLFS   ! index for level of free sink
INTEGER, DIMENSION(KLON),  INTENT(IN) :: KLCL   ! index lifting condens. level
INTEGER, DIMENSION(KLON),  INTENT(IN) :: KCTL   ! index for cloud top level
INTEGER, DIMENSION(KLON),  INTENT(IN) :: KDPL   ! index for departure level
INTEGER, DIMENSION(KLON),  INTENT(IN) :: KPBL   ! index for top of source layer
INTEGER, DIMENSION(KLON),  INTENT(IN) :: KML    ! index for melting level
REAL(dp), DIMENSION(KLON),  INTENT(INOUT) :: PTIMEC ! convection time step
REAL(dp), DIMENSION(KLON),     INTENT(IN) :: PDXDY  ! grid area (m^2)
REAL(dp), DIMENSION(KLON,KLEV),INTENT(IN) :: PTHL   ! grid scale enthalpy (J/kg)
REAL(dp), DIMENSION(KLON,KLEV),INTENT(IN) :: PTH    ! grid scale theta
REAL(dp), DIMENSION(KLON,KLEV),INTENT(IN) :: PRW    ! grid scale total water
                                                  ! mixing ratio
REAL(dp), DIMENSION(KLON,KLEV),INTENT(IN) :: PRC    ! grid scale r_c
REAL(dp), DIMENSION(KLON,KLEV),INTENT(IN) :: PRI    ! grid scale r_i
LOGICAL, DIMENSION(KLON),  INTENT(IN) :: OTRIG1   ! logical to keep trace of
                                                  ! convective arrays modified in UPDRAFT


REAL(dp), DIMENSION(KLON,KLEV), INTENT(IN) :: PPRES  ! pressure (P)
REAL(dp), DIMENSION(KLON,KLEV), INTENT(IN) :: PDPRES ! pressure difference between
                                                   ! bottom and top of layer (Pa)
REAL(dp), DIMENSION(KLON,KLEV), INTENT(IN) :: PLMASS ! mass of model layer (kg)
REAL(dp), DIMENSION(KLON,KLEV), INTENT(IN) :: PZ     ! height of model layer (m)
REAL(dp), DIMENSION(KLON),     INTENT(IN)  :: PCAPE  ! available potent. energy
INTEGER,                INTENT(OUT)   :: KFTSTEPS! maximum of fract time steps
                                                   ! only used for chemical tracers


REAL(dp), DIMENSION(KLON,KLEV), INTENT(INOUT):: PUMF  ! updraft mass flux (kg/s)
REAL(dp), DIMENSION(KLON,KLEV), INTENT(INOUT):: PUER  ! updraft entrainment (kg/s)
REAL(dp), DIMENSION(KLON,KLEV), INTENT(INOUT):: PUDR  ! updraft detrainment (kg/s)
REAL(dp), DIMENSION(KLON,KLEV), INTENT(INOUT):: PUPR  ! updraft precipitation in
                                                    ! flux units (kg water / s)
REAL(dp), DIMENSION(KLON,KLEV), INTENT(IN)  :: PUTHL  ! updraft enthalpy (J/kg)
REAL(dp), DIMENSION(KLON,KLEV), INTENT(IN)  :: PURW   ! updraft total water (kg/kg)
REAL(dp), DIMENSION(KLON,KLEV), INTENT(IN)  :: PURC   ! updraft cloud water (kg/kg)
REAL(dp), DIMENSION(KLON,KLEV), INTENT(IN)  :: PURI   ! updraft cloud ice   (kg/kg)

REAL(dp), DIMENSION(KLON,KLEV), INTENT(INOUT):: PDMF  ! downdraft mass flux (kg/s)
REAL(dp), DIMENSION(KLON,KLEV), INTENT(INOUT):: PDER  ! downdraft entrainment (kg/s)
REAL(dp), DIMENSION(KLON,KLEV), INTENT(INOUT):: PDDR  ! downdraft detrainment (kg/s)
REAL(dp), DIMENSION(KLON,KLEV), INTENT(IN)   :: PDTHL ! downdraft enthalpy (J/kg)
REAL(dp), DIMENSION(KLON,KLEV), INTENT(IN)   :: PDRW  ! downdraft total water (kg/kg)
REAL(dp), DIMENSION(KLON),      INTENT(INOUT):: PTPR  ! total surf precipitation (kg/s)
REAL(dp), DIMENSION(KLON),      INTENT(OUT)  :: PSPR  ! solid surf precipitation (kg/s)
REAL(dp), DIMENSION(KLON),      INTENT(INOUT):: PDTEVR! donwndraft evapor. (kg/s)

REAL(dp), DIMENSION(KLON,KLEV), INTENT(OUT)  :: PTHC  ! conv. adj. grid scale theta
REAL(dp), DIMENSION(KLON,KLEV), INTENT(OUT)  :: PRWC  ! conv. adj. grid scale r_w
REAL(dp), DIMENSION(KLON,KLEV), INTENT(OUT)  :: PRCC  ! conv. adj. grid scale r_c
REAL(dp), DIMENSION(KLON,KLEV), INTENT(OUT)  :: PRIC  ! conv. adj. grid scale r_i
REAL(dp), DIMENSION(KLON,KLEV), INTENT(OUT)  :: PWSUB ! envir. compensating subsidence(Pa/s)

REAL(dp), DIMENSION(KLON,KLEV), INTENT(INOUT):: PDTEVRF! downdraft evaporation rate
REAL(dp), DIMENSION(KLON,KLEV), INTENT(OUT)  :: PPRLFLX! liquid precip flux
REAL(dp), DIMENSION(KLON,KLEV), INTENT(OUT)  :: PPRSFLX! solid  precip flux

!*       0.2   Declarations of local variables :

INTEGER :: IIE, IKB, IKE  ! horizontal + vertical loop bounds
INTEGER :: IKS            ! vertical dimension
INTEGER :: JK, JKP, JKMAX ! vertical loop index
INTEGER :: JI             ! horizontal loop index
INTEGER :: JITER          ! iteration loop index
INTEGER :: JSTEP          ! fractional time loop index
 REAL(dp)    :: ZCPORD, ZRDOCP ! C_pd / R_d, R_d / C_pd
!REAL(dp)    :: ZCVOCD, ZEPSA  ! C_pv / C_pd, R_v / R_d

REAL(dp), DIMENSION(KLON,KLEV) :: ZTHLC       ! convectively adjusted
                                            ! grid scale enthalpy
REAL(dp), DIMENSION(KLON,KLEV) :: ZOMG        ! conv. environm. subsidence (Pa/s)
REAL(dp), DIMENSION(KLON,KLEV) :: ZUMF        ! non-adjusted updraft mass flux
REAL(dp), DIMENSION(KLON,KLEV) :: ZUER        !   "     updraft  entrainm. rate
REAL(dp), DIMENSION(KLON,KLEV) :: ZUDR        !   "     updraft  detrainm. rate
REAL(dp), DIMENSION(KLON,KLEV) :: ZDMF        !   "   downdraft mass flux
REAL(dp), DIMENSION(KLON,KLEV) :: ZDER        !   "   downdraft  entrainm. rate
REAL(dp), DIMENSION(KLON,KLEV) :: ZDDR        !   "   downdraft  detrainm. rate
REAL(dp), DIMENSION(KLON)     :: ZTPR         !   "   total precipitation
REAL(dp), DIMENSION(KLON)     :: ZDTEVR       !   "   total downdraft evapor.
REAL(dp), DIMENSION(KLON,KLEV):: ZPRLFLX      !   "   liquid precip flux
REAL(dp), DIMENSION(KLON,KLEV):: ZPRSFLX      !   "   solid  precip flux
REAL(dp), DIMENSION(KLON)     :: ZPRMELT      ! melting of precipitation
REAL(dp), DIMENSION(KLON)     :: ZPRMELTO     ! non-adjusted  "
REAL(dp), DIMENSION(KLON)     :: ZADJ         ! mass adjustment factor
REAL(dp), DIMENSION(KLON)     :: ZADJMAX      ! limit value for ZADJ
REAL(dp), DIMENSION(KLON)     :: ZCAPE        ! new CAPE after adjustment
REAL(dp), DIMENSION(KLON)     :: ZTIMEC       ! fractional convective time step
REAL(dp), DIMENSION(KLON,KLEV):: ZTIMC        ! 2D work array for ZTIMEC

REAL(dp), DIMENSION(KLON)     :: ZTHLCL       ! new  theta at LCL
REAL(dp), DIMENSION(KLON)     :: ZRVLCL       ! new  r_v at LCL
REAL(dp), DIMENSION(KLON)     :: ZZLCL        ! height of LCL
REAL(dp), DIMENSION(KLON)     :: ZTLCL        ! temperature at LCL
REAL(dp), DIMENSION(KLON)     :: ZTELCL       ! envir. temper. at LCL
REAL(dp), DIMENSION(KLON)     :: ZTHEU1, ZTHEU2! theta_e for undilute ascent
REAL(dp), DIMENSION(KLON)     :: ZTHES1, ZTHES2! saturation environm. theta_e
REAL(dp), DIMENSION(KLON,KLEV) :: ZTHMFIN, ZTHMFOUT, ZRWMFIN, ZRWMFOUT
REAL(dp), DIMENSION(KLON,KLEV) :: ZRCMFIN, ZRCMFOUT, ZRIMFIN, ZRIMFOUT
                                    ! work arrays for environm. compensat. mass flux
REAL(dp), DIMENSION(KLON)     :: ZPI          ! (P/P00)**R_d/C_pd
REAL(dp), DIMENSION(KLON)     :: ZLV          ! latent heat of vaporisation
REAL(dp), DIMENSION(KLON)     :: ZLS          ! latent heat of sublimation
REAL(dp), DIMENSION(KLON)     :: ZLM          ! latent heat of melting
REAL(dp), DIMENSION(KLON)     :: ZCPH         ! specific heat C_ph
REAL(dp), DIMENSION(KLON)     :: ZMELDPTH     ! actual depth of melting layer
INTEGER, DIMENSION(KLON)  :: ITSTEP       ! fractional convective time step
INTEGER, DIMENSION(KLON)  :: ICOUNT       ! timestep counter
INTEGER, DIMENSION(KLON)  :: ILCL         ! index lifting condens. level
INTEGER, DIMENSION(KLON)  :: IWORK1       ! work array
REAL(dp),  DIMENSION(KLON)      :: ZWORK1, ZWORK2, ZWORK3  ! work arrays
REAL(dp),  DIMENSION(KLON)      :: ZWORK4, ZWORK5          ! work arrays
REAL(dp),  DIMENSION(KLON,KLEV) :: ZWORK6                  ! work array
LOGICAL, DIMENSION(KLON)      :: GWORK1, GWORK3          ! work arrays
LOGICAL, DIMENSION(KLON,KLEV) :: GWORK4                  ! work array


!-------------------------------------------------------------------------------

!*       0.2    Initialize  local variables
!               ----------------------------


PSPR(:)     = 0.0_dp

ZTIMC(:,:)  = 0.0_dp
ZTHES2(:)   = 0.0_dp

ZWORK1(:)   = 0.0_dp
ZWORK2(:)   = 0.0_dp
ZWORK3(:)   = 0.0_dp
ZWORK4(:)   = 0.0_dp
ZWORK5(:)   = 0.0_dp
ZWORK6(:,:) = 0.0_dp

GWORK1(:)   = .FALSE.
GWORK3(:)   = .FALSE.
GWORK4(:,:) = .FALSE.

ILCL(:)     = KLCL(:)

 ZCPORD    = RCPD / RD
 ZRDOCP    = RD  / RCPD
!ZCVOCD    = RCPV / RCPD
!ZEPSA     = RV  / RD

ZADJ(:)   = 1.0_dp
ZWORK5(:) = 1.0_dp
WHERE( .NOT. OTRIG1(:) ) ZWORK5(:) = 0.0_dp


!*       0.3   Compute loop bounds
!              -------------------

IIE    = KLON
IKB    = 1 + JCVEXB
IKS    = KLEV
IKE    = KLEV - JCVEXT
JKMAX  = MAXVAL( KCTL(:) )


!*       2.     Save initial mass flux values to be used in adjustment procedure
!               ---------------------------------------------------------------

ZUMF(:,:)  = PUMF(:,:)
ZUER(:,:)  = PUER(:,:)
ZUDR(:,:)  = PUDR(:,:)
ZDMF(:,:)  = PDMF(:,:)
ZDER(:,:)  = PDER(:,:)
ZDDR(:,:)  = PDDR(:,:)
ZTPR(:)    = PTPR(:)
ZDTEVR(:)  = PDTEVR(:)
ZOMG(:,:)  = 0.0_dp
PWSUB(:,:) = 0.0_dp
ZPRMELT(:) = 0.0_dp
PPRLFLX(:,:) = 0.0_dp
ZPRLFLX(:,:) = 0.0_dp
PPRSFLX(:,:) = 0.0_dp
ZPRSFLX(:,:) = 0.0_dp


!*       2.1    Some preliminar computations for melting of precipitation
!               used later in section 9 and computation of precip fluxes
!               Precipitation fluxes are updated for melting and evaporation
!               ---------------------------------------------------------


ZWORK1(:) = 0.0_dp
ZMELDPTH(:) = 0.0_dp
ZWORK6(:,:) = 0.0_dp

DO JK = JKMAX + 1, IKB + 1, -1
   ! Nota: PUPR is total precipitation flux, but the solid, liquid
   !       precipitation is stored in units kg/kg; therefore we compute here
   !       the solid fraction of the total precipitation flux.
  DO JI = 1, IIE
     ZWORK2(JI)    = PUPR(JI,JK) / ( PURC(JI,JK) + PURI(JI,JK) + 1.E-8_dp )
     ZPRMELT(JI)   = ZPRMELT(JI) + PURI(JI,JK) * ZWORK2(JI)
     ZWORK1(JI)    = ZWORK1(JI) + PURC(JI,JK) * ZWORK2(JI) - (PDTEVRF(JI,JK)-PDTEVRF(JI,JK-1))
     ZPRLFLX(JI,JK)= MAX( 0.0_dp, ZWORK1(JI) )
     ZPRMELT(JI)   = ZPRMELT(JI) + MIN( 0.0_dp, ZWORK1(JI) )
     ZPRSFLX(JI,JK)= ZPRMELT(JI)
     IF ( KML(JI) >= JK .AND. ZMELDPTH(JI) <= XMELDPTH ) THEN
          ZPI(JI)    = ( PPRES(JI,JK) / RATM ) ** ZRDOCP
          ZWORK3(JI) = PTH(JI,JK) * ZPI(JI)            ! temperature estimate
          ZLM(JI)    = RLSTT + ( RCPV - RCS ) * ( ZWORK3(JI) - RTT ) -         &
                   & ( RLVTT + ( RCPV - RCW ) * ( ZWORK3(JI) - RTT ) ) ! L_s - L_v
          ZCPH(JI)   = RCPD + RCPV * PRW(JI,JK)
          ZMELDPTH(JI) = ZMELDPTH(JI) + PDPRES(JI,JK)
        ! ZWORK6(JI,JK)= ZLM(JI) * PTIMEC(JI) / PLMASS(JI,JK) * PDPRES(JI,JK)
          ZWORK6(JI,JK)= ZLM(JI) * PDPRES(JI,JK)
          ZOMG(JI,JK)= 1.0_dp ! at this place only used as work variable
     ENDIF
  ENDDO

ENDDO

ZWORK2(:) = 0.0_dp
DO JK = JKMAX, IKB + 1, -1
    ZWORK1(:) = ZPRMELT(:) * PDPRES(:,JK) / MAX( XMELDPTH, ZMELDPTH(:) )
    ZWORK2(:) = ZWORK2(:) + ZWORK1(:) * ZOMG(:,JK)
    ZPRLFLX(:,JK) = ZPRLFLX(:,JK) + ZWORK2(:)
    ZPRSFLX(:,JK) = ZPRSFLX(:,JK) - ZWORK2(:)
ENDDO
WHERE( ZPRSFLX(:,:) < 1.0_dp ) ZPRSFLX(:,:)=0.0_dp
ZPRMELTO(:) = ZPRMELT(:)


!*       3.     Compute limits on the closure adjustment factor so that the
!               inflow in convective drafts from a given layer can't be larger
!               than the mass contained in this layer initially.
!               ---------------------------------------------------------------

ZADJMAX(:) = 1000._dp
IWORK1(:)  = MAX( ILCL(:), KLFS(:) )
JKP        = MINVAL( KDPL(:) )

DO JK = JKP, IKE
  DO JI = 1, IIE
    IF( JK > KDPL(JI) .AND. JK <= IWORK1(JI) ) THEN
        ZWORK1(JI)  = PLMASS(JI,JK) /                                      &
                & ( ( PUER(JI,JK) + PDER(JI,JK) + 1.E-5_dp ) * PTIMEC(JI) )
        ZADJMAX(JI) = MIN( ZADJMAX(JI), ZWORK1(JI) )
    ENDIF
  ENDDO
ENDDO


GWORK1(:) = OTRIG1(:)  ! logical array to limit adjustment to not definitively
                       ! adjusted columns

DO JK = IKB, IKE
  ZTHLC(:,JK) = PTHL(:,JK) ! initialize adjusted envir. values
  PRWC(:,JK)  = PRW(:,JK)
  PRCC(:,JK)  = PRC(:,JK)
  PRIC(:,JK)  = PRI(:,JK)
  PTHC(:,JK)  = PTH(:,JK)
ENDDO



DO JITER = 1, 6  ! Enter adjustment loop to assure that all CAPE is
                 ! removed within the advective time interval TIMEC

     ZTIMEC(:) = PTIMEC(:)
     GWORK4(:,:)   = SPREAD( GWORK1(:), DIM=2, NCOPIES=IKS )
     WHERE( GWORK4(:,:) ) PWSUB(:,:) = 0.0_dp
     ZOMG(:,:)=0.0_dp

     DO JK = IKB + 1, JKMAX
           JKP = MAX( IKB + 1, JK - 1 )
           WHERE ( GWORK1(:) .AND. JK <= KCTL(:) )


!*       4.     Determine vertical velocity at top and bottom of each layer
!               to satisfy mass continuity.
!               ---------------------------------------------------------------
              ! we compute here Domega/Dp = - g rho Dw/Dz = 1/Dt

             ZWORK1(:)   = - ( PUER(:,JKP) + PDER(:,JKP) -                   &
                         & PUDR(:,JKP) - PDDR(:,JKP) ) / PLMASS(:,JKP)

             PWSUB(:,JK) = PWSUB(:,JKP) - PDPRES(:,JK-1) * ZWORK1(:)
              ! we use PDPRES(JK-1) and not JKP in order to have zero subsidence
              ! at the first layer


!*       5.     Compute fractional time step. For stability or
!               mass conservation reasons one must split full time step PTIMEC)
!               ---------------------------------------------------------------

             ZWORK1(:) = XSTABT * PDPRES(:,JKP) / ( ABS( PWSUB(:,JK) ) + 1.E-10_dp )
              ! the factor XSTABT is used for stability reasons
             ZTIMEC(:) = MIN( ZTIMEC(:), ZWORK1(:) )

              ! transform vertical velocity in mass flux units
             ZOMG(:,JK) = PWSUB(:,JK) * PDXDY(:) / RG
         END WHERE
     ENDDO


     WHERE( GWORK4(:,:) )
           ZTHLC(:,:) = PTHL(:,:) ! reinitialize adjusted envir. values
           PRWC(:,:)  = PRW(:,:)  ! when iteration criterium not attained
           PRCC(:,:)  = PRC(:,:)
           PRIC(:,:)  = PRI(:,:)
           PTHC(:,:)  = PTH(:,:)
     END WHERE


!        6. Check for mass conservation, i.e. ZWORK1 > 1.E-2
!           If mass is not conserved, the convective tendencies
!           automatically become zero.
!           ----------------------------------------------------

    DO JI = 1, IIE
       JK=KCTL(JI)
       ZWORK1(JI) = PUDR(JI,JK) * PDPRES(JI,JK) / ( PLMASS(JI,JK) + .1_dp )    &
                  &                                         - PWSUB(JI,JK)
    ENDDO
    ITSTEP(:) = INT( PTIMEC(:) / ZTIMEC(:) ) + 1
    WHERE( ( GWORK1(:) .AND. ABS( ZWORK1(:) ) - .01_dp > 0.0_dp ) .OR. ITSTEP(:) > 40 )
        GWORK1(:) = .FALSE.
        PTIMEC(:) = 1.E-1_dp
        ZTPR(:)   = 0.0_dp
        ZWORK5(:) = 0.0_dp
    END WHERE
    DO JK = IKB, IKE
        PWSUB(:,JK) = PWSUB(:,JK) * ZWORK5(:)
        ZPRLFLX(:,JK) = ZPRLFLX(:,JK) * ZWORK5(:)
        ZPRSFLX(:,JK) = ZPRSFLX(:,JK) * ZWORK5(:)
    ENDDO
    GWORK4(:,1:IKB) = .FALSE.
    GWORK4(:,IKE:IKS) = .FALSE.

    ZTIMEC(:) = PTIMEC(:) / REAL( ITSTEP(:) ) ! adjust  fractional time step
                                              ! to be an integer multiple of PTIMEC
    ZTIMC(:,:)= SPREAD( ZTIMEC(:), DIM=2, NCOPIES=IKS )
    ICOUNT(:) = 0


    ZTIMC(:,IKB+1:JKMAX) = ZTIMC(:,IKB+1:JKMAX) / PLMASS(:,IKB+1:JKMAX)

    KFTSTEPS = MAXVAL( ITSTEP(:) )
    DO JSTEP = 1, KFTSTEPS ! Enter the fractional time step loop here

         ICOUNT(:) = ICOUNT(:) + 1

             GWORK3(:) =  ITSTEP(:) >= ICOUNT(:) .AND. GWORK1(:)


!*       7.     Assign enthalpy and r_w values at the top and bottom of each
!               layer based on the sign of w
!               ------------------------------------------------------------

             ZTHMFIN(:,:)   = 0.0_dp
             ZRWMFIN(:,:)   = 0.0_dp
             ZRCMFIN(:,:)   = 0.0_dp
             ZRIMFIN(:,:)   = 0.0_dp
             ZTHMFOUT(:,:)  = 0.0_dp
             ZRWMFOUT(:,:)  = 0.0_dp
             ZRCMFOUT(:,:)  = 0.0_dp
             ZRIMFOUT(:,:)  = 0.0_dp

         DO JK = IKB + 1, JKMAX
           JKP = MAX( IKB + 1, JK - 1 )
           DO JI = 1, IIE
           GWORK4(JI,JK) = GWORK3(JI) .AND. JK <= KCTL(JI)
           IF ( GWORK3(JI) ) THEN
               ZWORK1(JI)       = SIGN( 1.0_dp, ZOMG(JI,JK) )
               ZWORK2(JI)       = 0.5_dp * ( 1.0_dp + ZWORK1(JI) )
               ZWORK1(JI)       = 0.5_dp * ( 1.0_dp - ZWORK1(JI) )
               ZTHMFIN(JI,JK)   = - ZOMG(JI,JK) * ZTHLC(JI,JKP) * ZWORK1(JI)
               ZTHMFOUT(JI,JK)  =   ZOMG(JI,JK) * ZTHLC(JI,JK)  * ZWORK2(JI)
               ZRWMFIN(JI,JK)   = - ZOMG(JI,JK) * PRWC(JI,JKP) * ZWORK1(JI)
               ZRWMFOUT(JI,JK)  =   ZOMG(JI,JK) * PRWC(JI,JK)  * ZWORK2(JI)
               ZRCMFIN(JI,JK)   = - ZOMG(JI,JK) * PRCC(JI,JKP) * ZWORK1(JI)
               ZRCMFOUT(JI,JK)  =   ZOMG(JI,JK) * PRCC(JI,JK)  * ZWORK2(JI)
               ZRIMFIN(JI,JK)   = - ZOMG(JI,JK) * PRIC(JI,JKP) * ZWORK1(JI)
               ZRIMFOUT(JI,JK)  =   ZOMG(JI,JK) * PRIC(JI,JK)  * ZWORK2(JI)
           ENDIF
           ENDDO
           DO JI = 1, IIE
           IF ( GWORK3(JI) ) THEN
               ZTHMFIN(JI,JKP)  = ZTHMFIN(JI,JKP)  + ZTHMFOUT(JI,JK) * ZWORK2(JI)
               ZTHMFOUT(JI,JKP) = ZTHMFOUT(JI,JKP) + ZTHMFIN(JI,JK)  * ZWORK1(JI)
               ZRWMFIN(JI,JKP)  = ZRWMFIN(JI,JKP)  + ZRWMFOUT(JI,JK) * ZWORK2(JI)
               ZRWMFOUT(JI,JKP) = ZRWMFOUT(JI,JKP) + ZRWMFIN(JI,JK)  * ZWORK1(JI)
               ZRCMFIN(JI,JKP)  = ZRCMFIN(JI,JKP)  + ZRCMFOUT(JI,JK) * ZWORK2(JI)
               ZRCMFOUT(JI,JKP) = ZRCMFOUT(JI,JKP) + ZRCMFIN(JI,JK)  * ZWORK1(JI)
               ZRIMFIN(JI,JKP)  = ZRIMFIN(JI,JKP)  + ZRIMFOUT(JI,JK) * ZWORK2(JI)
               ZRIMFOUT(JI,JKP) = ZRIMFOUT(JI,JKP) + ZRIMFIN(JI,JK)  * ZWORK1(JI)
           ENDIF
           ENDDO
         ENDDO

         WHERE ( GWORK4(:,:) )

!******************************************************************************

!*       8.     Update the environmental values of enthalpy and r_w at each level
!               NOTA: These are the MAIN EQUATIONS of the scheme
!               -----------------------------------------------------------------


           ZTHLC(:,:) = ZTHLC(:,:) + ZTIMC(:,:) * (                  &
                      &   ZTHMFIN(:,:) + PUDR(:,:) * PUTHL(:,:)  +   &
                      &   PDDR(:,:) * PDTHL(:,:) - ZTHMFOUT(:,:) -   &
                      & ( PUER(:,:) + PDER(:,:) ) * PTHL(:,:)   )

           PRWC(:,:)  = PRWC(:,:) + ZTIMC(:,:)  *  (                 &
                      &   ZRWMFIN(:,:) + PUDR(:,:) * PURW(:,:)  +    &
                      &   PDDR(:,:) * PDRW(:,:) - ZRWMFOUT(:,:) -    &
                      & ( PUER(:,:) + PDER(:,:) ) * PRW(:,:)    )

           PRCC(:,:)  = PRCC(:,:) + ZTIMC(:,:)  *  (                 &
                      &   ZRCMFIN(:,:) + PUDR(:,:) * PURC(:,:)       &
                      &                         - ZRCMFOUT(:,:) -    &
                      & ( PUER(:,:) + PDER(:,:) ) * PRC(:,:)    )

           PRIC(:,:)  = PRIC(:,:) + ZTIMC(:,:)  *  (                 &
                      &   ZRIMFIN(:,:) + PUDR(:,:) * PURI(:,:)       &
                      &                         - ZRIMFOUT(:,:) -    &
                      & ( PUER(:,:) + PDER(:,:) ) * PRI(:,:)    )


!******************************************************************************
           
         END WHERE

    ENDDO ! Exit the fractional time step loop

!!$    do jk=ikb,ike
!!$      do ji=1,iie
!!$        if ( jk == KCTL(ji) .and. jk > 12 .and.&
!!$          (PRWC(ji,JK) - PRCC(ji,JK) - PRIC(ji,JK)) < 0.0_dp)&
!!$          print*, "in clsoure:" ,PRCC(ji,jk),  ZRCMFIN(ji,jk),&
!!$          PUDR(ji,jk), &
!!$          PURC(ji,jk),  ZRCMFIN(ji,jk) + PUDR(ji,jk) * PURC(ji,jk), ZRCMFOUT(ji,jk),  &
!!$          (PUER(ji,jk) + PDER(ji,jk) ) * PRC(ji,jk), prc(ji,jk), ji, jk,              & 
!!$          prwc(ji,jk), pric(ji,jk)
!!$
!!$      end do
!!$    end do

!*           9.    Allow frozen precipitation to melt over a 200 mb deep layer
!                  -----------------------------------------------------------

   !  DO JK = JKMAX, IKB + 1, -1
          ! ZTHLC(:,JK) = ZTHLC(:,JK) -                                &
          ! &  ZPRMELT(:) * ZWORK6(:,JK) / MAX( XMELDPTH, ZMELDPTH(:) )
          ! ZTHLC(:,JK) = ZTHLC(:,JK) -                                &
          !    ZPRMELT(:) * ZWORK6(:,JK) / ( pumf(:,jk)*MAX( XMELDPTH, ZMELDPTH(:))+1.E-10_dp )
   !  ENDDO


!*          10.    Compute final linearized value of theta envir.
!                  ----------------------------------------------

      DO JK = IKB + 1, JKMAX
         DO JI = 1, IIE
         IF( GWORK1(JI) .AND. JK <= KCTL(JI) ) THEN
           ZPI(JI)    = ( RATM / PPRES(JI,JK) ) ** ZRDOCP
           ZCPH(JI)   = RCPD + PRWC(JI,JK) * RCPV
           ZWORK2(JI) = PTH(JI,JK) / ZPI(JI)  ! first temperature estimate
           ZLV(JI)    = RLVTT + ( RCPV - RCW ) * ( ZWORK2(JI) - RTT )
           ZLS(JI)    = RLVTT + ( RCPV - RCS ) * ( ZWORK2(JI) - RTT )
             ! final linearized temperature
           ZWORK2(JI) = ( ZTHLC(JI,JK) + ZLV(JI) * PRCC(JI,JK) + ZLS(JI) * PRIC(JI,JK) &
                      & - (1.0_dp + PRWC(JI,JK) ) * RG * PZ(JI,JK) ) / ZCPH(JI)
           ZWORK2(JI) = MAX( 180._dp, MIN( 340._dp, ZWORK2(JI) ) )
           PTHC(JI,JK)= ZWORK2(JI) * ZPI(JI) ! final adjusted envir. theta
         ENDIF
         ENDDO
      ENDDO


!*         11.     Compute new cloud ( properties at new LCL )
!                     NOTA: The computations are very close to
!                           that in routine TRIGGER_FUNCT
!                  ---------------------------------------------

      CALL CONVECT_CLOSURE_THRVLCL(  KLON, KLEV,                           &
                                  &  PPRES, PTHC, PRWC, PZ, GWORK1,        &
                                  &  ZTHLCL, ZRVLCL, ZZLCL, ZTLCL, ZTELCL, &
                                  &  ILCL, KDPL, KPBL )


       ZTLCL(:)  = MAX( 230._dp, MIN( 335._dp, ZTLCL(:)  ) )  ! set some overflow bounds
       ZTELCL(:) = MAX( 230._dp, MIN( 335._dp, ZTELCL(:) ) )
       ZTHLCL(:) = MAX( 230._dp, MIN( 345._dp, ZTHLCL(:) ) )
       ZRVLCL(:) = MAX(   0.0_dp,  MIN(   1.0_dp,   ZRVLCL(:) ) )


!*         12.    Compute adjusted CAPE
!                 ---------------------

       ZCAPE(:)  = 0.0_dp
       ZPI(:)    = ZTHLCL(:) / ZTLCL(:)
       ZPI(:)    = MAX( 0.95_dp, MIN( 1.5_dp, ZPI(:) ) )
       ZWORK1(:) = RATM / ZPI(:) ** ZCPORD ! pressure at LCL

       CALL CONVECT_SATMIXRATIO( KLON, ZWORK1, ZTELCL, ZWORK3, ZLV, ZLS, ZCPH )
       ZWORK3(:) = MIN(   .1_dp, MAX(   0.0_dp, ZWORK3(:) ) )

                ! compute theta_e updraft undilute
       ZTHEU1(:) = ZTLCL(:) * ZPI(:) ** ( 1.0_dp - 0.28_dp * ZRVLCL(:) )             &
                 &                * EXP( ( 3374.6525_dp / ZTLCL(:) - 2.5403_dp )  &
                 &                * ZRVLCL(:) * ( 1.0_dp + 0.81_dp * ZRVLCL(:) ) )

                ! compute theta_e saturated environment at LCL
       ZTHES1(:) = ZTELCL(:) * ZPI(:) ** ( 1.0_dp - 0.28_dp * ZWORK3(:) )            &
                 &                * EXP( ( 3374.6525_dp / ZTELCL(:) - 2.5403_dp ) &
                 &                * ZWORK3(:) * ( 1.0_dp + 0.81_dp * ZWORK3(:) ) )

      DO JK = MINVAL( ILCL(:) ), JKMAX
        JKP = JK - 1
        DO JI = 1, IIE
          ZWORK4(JI) = 1.0_dp
          IF ( JK == ILCL(JI) ) ZWORK4(JI) = 0.0_dp

           ! compute theta_e saturated environment and adjusted values
           ! of theta

          GWORK3(JI)  = JK >= ILCL(JI) .AND. JK <= KCTL(JI) .AND. GWORK1(JI)

          ZPI(JI)     = ( RATM / PPRES(JI,JK) ) ** ZRDOCP
          ZWORK2(JI)  = PTHC(JI,JK) / ZPI(JI)
        ENDDO

        CALL CONVECT_SATMIXRATIO( KLON, PPRES(:,JK), ZWORK2, ZWORK3, ZLV, ZLS, ZCPH )


        DO JI = 1, IIE
          IF ( GWORK3(JI) ) THEN
              ZTHES2(JI)  = ZWORK2(JI) * ZPI(JI) ** ( 1.0_dp - 0.28_dp * ZWORK3(JI) )  &
                          &        * EXP( ( 3374.6525_dp / ZWORK2(JI) - 2.5403_dp ) &
                          &        * ZWORK3(JI) * ( 1.0_dp + 0.81_dp * ZWORK3(JI) ) )

              ZWORK3(JI)  = PZ(JI,JK) - PZ(JI,JKP) * ZWORK4(JI) -                       &
                         & ( 1.0_dp - ZWORK4(JI) ) * ZZLCL(JI)    ! level thickness
              ZWORK1(JI)  = MAX( 0.01_dp, PUER(JI,JK)/ MAX( .1_dp, PUMF(JI,JK) ) )
              ZTHEU2(JI)  = ( 1.0_dp - ZWORK1(JI) ) * ZTHEU1(JI) + ZWORK1(JI) * ZTHES1(JI)
              ZWORK1(JI)  = ( ZTHEU1(JI) + ZTHEU2(JI) ) / ( ZTHES1(JI) + ZTHES2(JI) ) - 1.0_dp
              ZCAPE(JI)   = ZCAPE(JI) + RG * ZWORK3(JI) * MAX( 0.0_dp, ZWORK1(JI) )
              ZTHEU1(JI)  = ZTHEU2(JI)
              ZTHES1(JI)  = ZTHES2(JI)
          ENDIF
        ENDDO
      ENDDO


!*         13.     Determine mass adjustment factor knowing how much
!                  CAPE has been removed.
!                  -------------------------------------------------

       WHERE ( GWORK1(:) )
           ZWORK1(:) = MAX( PCAPE(:) - ZCAPE(:), 0.2_dp * PCAPE(:) )
           ZWORK2(:) = ZCAPE(:) / ( PCAPE(:) + 1.E-8_dp )

           GWORK1(:) = ZWORK2(:) > 0.2_dp .OR. ZCAPE(:) == 0.0_dp ! mask for adjustment
       END WHERE

       WHERE ( ZCAPE(:) == 0.0_dp .AND. GWORK1(:) )  ZADJ(:) = ZADJ(:) * 0.5_dp
       WHERE ( ZCAPE(:) /= 0.0_dp .AND. GWORK1(:) )                                    &
             & ZADJ(:) = ZADJ(:) * XSTABC * PCAPE(:) / ( ZWORK1(:) + 1.E-8_dp )
       ZADJ(:) = MIN( ZADJ(:), ZADJMAX(:) )


!*         13.     Adjust mass flux by the factor ZADJ to converge to
!                  specified degree of stabilization
!                 ----------------------------------------------------

       CALL CONVECT_CLOSURE_ADJUST( KLON, KLEV, ZADJ,                     &
                                  & PUMF, ZUMF, PUER, ZUER, PUDR, ZUDR,   &
                                  & PDMF, ZDMF, PDER, ZDER, PDDR, ZDDR,   &
                                  & ZPRMELT, ZPRMELTO, PDTEVR, ZDTEVR,    &
                                  & PTPR, ZTPR,                           &
                                  & PPRLFLX, ZPRLFLX, PPRSFLX, ZPRSFLX    )


      IF ( COUNT( GWORK1(:) ) == 0 ) EXIT ! exit big adjustment iteration loop
                                          ! when all columns have reached
                                          ! desired degree of stabilization.

ENDDO  ! end of big adjustment iteration loop


!!$do jk=ikb,ike
!!$  do ji=1,iie
!!$    if (jk == KCTL(JI) .and.  (PRWC(ji,JK) - PRCC(ji,JK) - PRIC(ji,JK)) < 0.0_dp) &
!!$      print*, PRWC(ji,JK), PRCC(ji,JK), PRIC(ji,JK), prw(ji,jk), &
!!$      (PRWC(ji,JK) - PRCC(ji,JK) - PRIC(ji,JK)), KCTL(ji), ji
!!$  end do
!!$End do

        ! skip adj. total water array  to water vapor
DO JK = IKB, IKE
  PRWC(:,JK) = MAX( 0.0_dp, PRWC(:,JK) - PRCC(:,JK) - PRIC(:,JK) )
ENDDO

END SUBROUTINE CONVECT_CLOSURE

!============================================================================

!######################################################################
 SUBROUTINE CONVECT_CLOSURE_THRVLCL( KLON, KLEV,                      &
                              &  PPRES, PTH, PRV, PZ, OWORK1,         &
                              &  PTHLCL, PRVLCL, PZLCL, PTLCL, PTELCL,&
                              &  KLCL, KDPL, KPBL )
!######################################################################

!!**** Determine thermodynamic properties at new LCL
!!
!!    PURPOSE
!!    -------
!!      The purpose of this routine is to determine the thermodynamic
!!      properties at the new lifting condensation level LCL
!!
!!
!!
!!**  METHOD
!!    ------
!!    see CONVECT_TRIGGER_FUNCT
!!
!!
!!
!!    EXTERNAL
!!    --------
!!     Routine CONVECT_SATMIXRATIO
!!
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!      Module YOMCST
!!          RG                 ! gravity constant
!!          RATM               ! Reference pressure
!!          RD, RV           ! Gaz  constants for dry air and water vapor
!!          RCPD               ! Cpd (dry air)
!!          RTT                ! triple point temperature
!!          RBETW, RGAMW      ! constants for vapor saturation pressure
!!
!!      Module YOE_CONVPAR
!!          XA25               ! reference grid area
!!          XZLCL              ! lowest allowed pressure difference between
!!                             ! surface and LCL
!!          XZPBL              ! minimum mixed layer depth to sustain convection
!!          XWTRIG             ! constant in vertical velocity trigger
!!
!!      Module YOE_CONVPAREXT
!!          JCVEXB, JCVEXT     ! extra levels on the vertical boundaries
!!
!!    REFERENCE
!!    ---------
!!
!!      Book2 of documentation ( routine TRIGGER_FUNCT)
!!      Fritsch and Chappell (1980), J. Atm. Sci., Vol. 37, 1722-1761.
!!
!!    AUTHOR
!!    ------
!!      P. BECHTOLD       * Laboratoire d'Aerologie *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    07/11/95
!!   Last modified  04/10/97
!-------------------------------------------------------------------------------

!#include "tsmbkind.h"

!*       0.    DECLARATIONS
!              ------------

!USE YOMCST
!USE YOE_CONVPAR
!USE YOE_CONVPAREXT


IMPLICIT NONE

!*       0.1   Declarations of dummy arguments :

INTEGER,                    INTENT(IN) :: KLON  ! horizontal dimension
INTEGER,                    INTENT(IN) :: KLEV  ! vertical dimension
REAL(dp), DIMENSION(KLON,KLEV), INTENT(IN) :: PTH   ! theta
REAL(dp), DIMENSION(KLON,KLEV), INTENT(IN) :: PRV   ! vapor mixing ratio
REAL(dp), DIMENSION(KLON,KLEV), INTENT(IN) :: PPRES ! pressure
REAL(dp), DIMENSION(KLON,KLEV), INTENT(IN) :: PZ    ! height of grid point (m)
INTEGER, DIMENSION(KLON),   INTENT(IN) :: KDPL  ! contains vert. index of DPL
INTEGER, DIMENSION(KLON),   INTENT(IN) :: KPBL  ! " vert. index of source layer top
LOGICAL,   DIMENSION(KLON),   INTENT(IN) :: OWORK1! logical mask

REAL(dp), DIMENSION(KLON),     INTENT(OUT):: PTHLCL ! theta at LCL
REAL(dp), DIMENSION(KLON),     INTENT(OUT):: PRVLCL ! vapor mixing ratio at  LCL
REAL(dp), DIMENSION(KLON),     INTENT(OUT):: PZLCL  ! height at LCL (m)
REAL(dp), DIMENSION(KLON),     INTENT(OUT):: PTLCL  ! temperature at LCL (m)
REAL(dp), DIMENSION(KLON),     INTENT(OUT):: PTELCL ! environm. temp. at LCL (K)
INTEGER, DIMENSION(KLON),  INTENT(OUT):: KLCL   ! contains vert. index of LCL

!*       0.2   Declarations of local variables :

INTEGER :: JK, JKM, JKMIN, JKMAX      ! vertical loop index
INTEGER :: JI                         ! horizontal loop index
INTEGER :: IIE, IKB, IKE              ! horizontal + vertical loop bounds
REAL(dp)    :: ZEPS, ZEPSA    ! R_d / R_v, R_v / R_d
REAL(dp)    :: ZCPORD, ZRDOCP ! C_pd / R_d, R_d / C_pd

REAL(dp), DIMENSION(KLON) :: ZPLCL    ! pressure at LCL
REAL(dp), DIMENSION(KLON) :: ZTMIX    ! mixed layer temperature
REAL(dp), DIMENSION(KLON) :: ZEVMIX   ! mixed layer water vapor pressure
REAL(dp), DIMENSION(KLON) :: ZDPTHMIX, ZPRESMIX ! mixed layer depth and pressure
REAL(dp), DIMENSION(KLON) :: ZLV, ZCPH! specific heats of vaporisation, dry air
REAL(dp), DIMENSION(KLON) :: ZDP      ! pressure between LCL and model layer
REAL(dp), DIMENSION(KLON) :: ZWORK1, ZWORK2     ! work arrays


!-------------------------------------------------------------------------------

!*       0.3    Compute array bounds
!               --------------------

IIE = KLON
IKB = 1 + JCVEXB
IKE = KLEV - JCVEXT


!*       1.     Initialize local variables
!               --------------------------

ZEPS      = RD  / RV
ZEPSA     = RV  / RD
ZCPORD    = RCPD / RD
ZRDOCP    = RD  / RCPD

ZDPTHMIX(:) = 0.0_dp
ZPRESMIX(:) = 0.0_dp
PTHLCL(:)   = 300._dp
PTLCL(:)    = 300._dp
PTELCL(:)   = 300._dp
PRVLCL(:)   = 0.0_dp
PZLCL(:)    = PZ(:,IKB)
ZTMIX(:)    = 230._dp
ZPLCL(:)    = 1.E4_dp
KLCL(:)     = IKB + 1


!*       2.     Construct a mixed layer as in TRIGGER_FUNCT
!               -------------------------------------------

     JKMAX = MAXVAL( KPBL(:) )
     JKMIN = MINVAL( KDPL(:) )
     DO JK = IKB + 1, JKMAX
        JKM = JK + 1
        DO JI = 1, IIE
        IF ( JK >= KDPL(JI) .AND. JK <= KPBL(JI) ) THEN

            ZWORK1(JI)   = PPRES(JI,JK) - PPRES(JI,JKM)
            ZDPTHMIX(JI) = ZDPTHMIX(JI) + ZWORK1(JI)
            ZPRESMIX(JI) = ZPRESMIX(JI) + PPRES(JI,JK) * ZWORK1(JI)
            PTHLCL(JI)   = PTHLCL(JI)   + PTH(JI,JK)   * ZWORK1(JI)
            PRVLCL(JI)   = PRVLCL(JI)   + PRV(JI,JK)   * ZWORK1(JI)

        ENDIF
        ENDDO
     ENDDO


WHERE ( OWORK1(:) )

        ZPRESMIX(:) = ZPRESMIX(:) / ZDPTHMIX(:)
        PTHLCL(:)   = PTHLCL(:)   / ZDPTHMIX(:)
        PRVLCL(:)   = PRVLCL(:)   / ZDPTHMIX(:)

!*       3.1    Use an empirical direct solution ( Bolton formula )
!               to determine temperature and pressure at LCL.
!               Nota: the adiabatic saturation temperature is not
!                     equal to the dewpoint temperature
!               --------------------------------------------------


        ZTMIX(:)  = PTHLCL(:) * ( ZPRESMIX(:) / RATM ) ** ZRDOCP
        ZEVMIX(:) = PRVLCL(:) * ZPRESMIX(:) / ( PRVLCL(:) + ZEPS )
        ZEVMIX(:) = MAX( 1.E-8_dp, ZEVMIX(:) )
        ZWORK1(:) = LOG( ZEVMIX(:) / 613.3_dp )
              ! dewpoint temperature
        ZWORK1(:) = ( 4780.8_dp - 32.19_dp * ZWORK1(:) ) / &
                  & ( 17.502_dp - ZWORK1(:) )
              ! adiabatic saturation temperature
        PTLCL(:)  = ZWORK1(:) - ( .212_dp + 1.571E-3_dp * ( ZWORK1(:) - RTT )   &
                  & - 4.36E-4_dp * ( ZTMIX(:) - RTT ) ) * ( ZTMIX(:) - ZWORK1(:) )
        PTLCL(:)  = MIN( PTLCL(:), ZTMIX(:) )
        ZPLCL(:)  = RATM * ( PTLCL(:) / PTHLCL(:) ) ** ZCPORD

END WHERE

     ZPLCL(:) = MIN( 2.E5_dp, MAX( 10._dp, ZPLCL(:) ) ) ! bound to avoid overflow


!*       3.2    Correct PTLCL in order to be completely consistent
!               with MNH saturation formula
!               --------------------------------------------------

     CALL CONVECT_SATMIXRATIO( KLON, ZPLCL, PTLCL, ZWORK1, ZLV, ZWORK2, ZCPH )
     WHERE( OWORK1(:) )
        ZWORK2(:) = ZWORK1(:) / PTLCL(:) * ( RBETW / PTLCL(:) - RGAMW ) ! dr_sat/dT
        ZWORK2(:) = ( ZWORK1(:) - PRVLCL(:) ) / &
                  &  ( 1.0_dp + ZLV(:) / ZCPH(:) * ZWORK2(:) )
        PTLCL(:)  = PTLCL(:) - ZLV(:) / ZCPH(:) * ZWORK2(:)

     END WHERE


!*       3.3    If PRVLCL is oversaturated set humidity and temperature
!               to saturation values.
!               -------------------------------------------------------

    CALL CONVECT_SATMIXRATIO( KLON, ZPRESMIX, ZTMIX, ZWORK1, ZLV, ZWORK2, ZCPH )
     WHERE( OWORK1(:) .AND. PRVLCL(:) > ZWORK1(:) )
        ZWORK2(:) = ZWORK1(:) / ZTMIX(:) * ( RBETW / ZTMIX(:) - RGAMW ) ! dr_sat/dT
        ZWORK2(:) = ( ZWORK1(:) - PRVLCL(:) ) / &
                  &  ( 1.0_dp + ZLV(:) / ZCPH(:) * ZWORK2(:) )
        PTLCL(:)  = ZTMIX(:) + ZLV(:) / ZCPH(:) * ZWORK2(:)
        PRVLCL(:) = PRVLCL(:) - ZWORK2(:)
        ZPLCL(:)  = ZPRESMIX(:)
        PTHLCL(:) = PTLCL(:) * ( RATM / ZPLCL(:) ) ** ZRDOCP
     END WHERE


!*        4.1   Determine  vertical loop index at the LCL
!               -----------------------------------------

     DO JK = JKMIN, IKE - 1
        DO JI = 1, IIE
        IF ( ZPLCL(JI) <= PPRES(JI,JK) .AND. OWORK1(JI) ) THEN
            KLCL(JI)  = JK + 1
            PZLCL(JI) = PZ(JI,JK+1)
        ENDIF
        ENDDO
     ENDDO


!*        4.2   Estimate height and environmental temperature at LCL
!               ----------------------------------------------------

    DO JI = 1, IIE
        JK   = KLCL(JI)
        JKM  = JK - 1
        ZDP(JI)     = LOG( ZPLCL(JI) / PPRES(JI,JKM) ) /   &
                    & LOG( PPRES(JI,JK) / PPRES(JI,JKM) )
        ZWORK1(JI)  = PTH(JI,JK)  * ( PPRES(JI,JK)  / RATM ) ** ZRDOCP
        ZWORK2(JI)  = PTH(JI,JKM) * ( PPRES(JI,JKM) / RATM ) ** ZRDOCP
        ZWORK1(JI)  = ZWORK2(JI) + ( ZWORK1(JI) - ZWORK2(JI) ) * ZDP(JI)
           ! we compute the precise value of the LCL
           ! The precise height is between the levels KLCL and KLCL-1.
        ZWORK2(JI) = PZ(JI,JKM) + ( PZ(JI,JK) - PZ(JI,JKM) ) * ZDP(JI)
    ENDDO
    WHERE( OWORK1(:) )
       PTELCL(:) = ZWORK1(:)
       PZLCL(:)  = ZWORK2(:)
    END WHERE



END SUBROUTINE CONVECT_CLOSURE_THRVLCL

!============================================================================

!    #########################################################################
     SUBROUTINE CONVECT_CLOSURE_ADJUST( KLON, KLEV, PADJ,                    &
                                  &   PUMF, PZUMF, PUER, PZUER, PUDR, PZUDR, &
                                  &   PDMF, PZDMF, PDER, PZDER, PDDR, PZDDR, &
                                  &   PPRMELT, PZPRMELT, PDTEVR, PZDTEVR,    &
                                  &   PTPR, PZTPR,                           &
                                  &   PPRLFLX, PZPRLFL, PPRSFLX, PZPRSFL     )
!    #########################################################################

!!**** Uses closure adjustment factor to adjust mass flux and to modify
!!     precipitation efficiency  when necessary. The computations are
!!     similar to routine CONVECT_PRECIP_ADJUST.
!!
!!
!!    PURPOSE
!!    -------
!!      The purpose of this routine is to adjust the mass flux using the
!!      factor PADJ computed in CONVECT_CLOSURE
!!
!!
!!**  METHOD
!!    ------
!!      Computations are done at every model level starting from bottom.
!!      The use of masks allows to optimise the inner loops (horizontal loops).
!!
!!
!!    EXTERNAL
!!    --------
!!     Module YOE_CONVPAREXT
!!          JCVEXB, JCVEXT     ! extra levels on the vertical boundaries
!!
!!    None
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!
!!    None
!!
!!    REFERENCE
!!    ---------
!!
!!      Book1,2 of documentation ( routine CONVECT_CLOSURE_ADJUST)
!!
!!    AUTHOR
!!    ------
!!      P. BECHTOLD       * Laboratoire d'Aerologie *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    26/03/96
!!   Last modified  04/10/97
!-------------------------------------------------------------------------------

!#include "tsmbkind.h"

!*       0.    DECLARATIONS
!              ------------

!USE YOE_CONVPAREXT

IMPLICIT NONE

!*       0.1   Declarations of dummy arguments :


INTEGER,                    INTENT(IN) :: KLON     ! horizontal dimension
INTEGER,                    INTENT(IN) :: KLEV     ! vertical dimension
REAL(dp), DIMENSION(KLON),      INTENT(IN) :: PADJ     ! mass adjustment factor


REAL(dp), DIMENSION(KLON,KLEV), INTENT(INOUT) :: PUMF  ! updraft mass flux (kg/s)
REAL(dp), DIMENSION(KLON,KLEV), INTENT(INOUT) :: PZUMF ! initial value of  "
REAL(dp), DIMENSION(KLON,KLEV), INTENT(INOUT) :: PUER  ! updraft entrainment (kg/s)
REAL(dp), DIMENSION(KLON,KLEV), INTENT(INOUT) :: PZUER ! initial value of  "
REAL(dp), DIMENSION(KLON,KLEV), INTENT(INOUT) :: PUDR  ! updraft detrainment (kg/s)
REAL(dp), DIMENSION(KLON,KLEV), INTENT(INOUT) :: PZUDR ! initial value of  "
REAL(dp), DIMENSION(KLON,KLEV), INTENT(INOUT) :: PDMF  ! downdraft mass flux (kg/s)
REAL(dp), DIMENSION(KLON,KLEV), INTENT(INOUT) :: PZDMF ! initial value of  "
REAL(dp), DIMENSION(KLON,KLEV), INTENT(INOUT) :: PDER  ! downdraft entrainment (kg/s)
REAL(dp), DIMENSION(KLON,KLEV), INTENT(INOUT) :: PZDER ! initial value of  "
REAL(dp), DIMENSION(KLON,KLEV), INTENT(INOUT) :: PDDR  ! downdraft detrainment (kg/s)
REAL(dp), DIMENSION(KLON,KLEV), INTENT(INOUT) :: PZDDR ! initial value of  "
REAL(dp), DIMENSION(KLON),   INTENT(INOUT):: PTPR     ! total precipitation (kg/s)
REAL(dp), DIMENSION(KLON),   INTENT(INOUT):: PZTPR    ! initial value of "
REAL(dp), DIMENSION(KLON),   INTENT(INOUT):: PDTEVR   ! donwndraft evapor. (kg/s)
REAL(dp), DIMENSION(KLON),   INTENT(INOUT):: PZDTEVR  ! initial value of "
REAL(dp), DIMENSION(KLON),   INTENT(INOUT):: PPRMELT  ! melting of precipitation
REAL(dp), DIMENSION(KLON),   INTENT(INOUT):: PZPRMELT ! initial value of "
REAL(dp), DIMENSION(KLON,KLEV),INTENT(INOUT)  :: PPRLFLX! liquid precip flux
REAL(dp), DIMENSION(KLON,KLEV),INTENT(INOUT)  :: PZPRLFL! initial value "
REAL(dp), DIMENSION(KLON,KLEV),INTENT(INOUT)  :: PPRSFLX! solid  precip flux
REAL(dp), DIMENSION(KLON,KLEV),INTENT(INOUT)  :: PZPRSFL! initial value "


!*       0.2   Declarations of local variables :

INTEGER :: IIE, IKB, IKE                 ! horiz. + vert. loop bounds
INTEGER :: JK                            ! vertical loop index


!-------------------------------------------------------------------------------

!*       0.3   Compute loop bounds
!              -------------------

IIE  = KLON
IKB  = 1 + JCVEXB
IKE  = KLEV - JCVEXT


!*       1.     Adjust mass flux by the factor PADJ to converge to
!               specified degree of stabilization
!               ----------------------------------------------------

          PPRMELT(:)  = PZPRMELT(:)   * PADJ(:)
          PDTEVR(:)   = PZDTEVR(:)    * PADJ(:)
          PTPR(:)     = PZTPR(:)      * PADJ(:)

     DO JK = IKB + 1, IKE
          PUMF(:,JK)  = PZUMF(:,JK)   * PADJ(:)
          PUER(:,JK)  = PZUER(:,JK)   * PADJ(:)
          PUDR(:,JK)  = PZUDR(:,JK)   * PADJ(:)
          PDMF(:,JK)  = PZDMF(:,JK)   * PADJ(:)
          PDER(:,JK)  = PZDER(:,JK)   * PADJ(:)
          PDDR(:,JK)  = PZDDR(:,JK)   * PADJ(:)
          PPRLFLX(:,JK) = PZPRLFL(:,JK) * PADJ(:)
          PPRSFLX(:,JK) = PZPRSFL(:,JK) * PADJ(:)
     ENDDO

   END SUBROUTINE CONVECT_CLOSURE_ADJUST

!============================================================================

!#######################################################################
  SUBROUTINE CONVECT_UV_TRANSPORT( KLON, KLEV, PU, PV, PUC, PVC,       &
                                 & KDPL, KPBL, KLCL, KCTL, KLFS, KDBL, &
                                 & PUMF, PUER, PUDR, PDMF, PDER, PDDR, &
                                 & PTIMEC, PDXDY, PMIXF, PLMASS, PWSUB,&
                                 & KFTSTEPS )
!#######################################################################

!!**** Compute  modified horizontal wind components due to convective event
!!
!!
!!    PURPOSE
!!    -------
!!      The purpose of this routine is to determine convective adjusted
!!      horizontal wind components u and v
!!      The final convective tendencies can then be evaluated in the main
!!      routine DEEP_CONVECT by (PUC-PU)/PTIMEC
!!
!!
!!**  METHOD
!!    ------
!!      Identical to the computation of the conservative variables in the
!!      tracer routine but includes pressure term
!!
!!    EXTERNAL
!!    --------
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!      Module YOMCST
!!          RG                 ! gravity constant
!!
!!     Module YOE_CONVPAREXT
!!          JCVEXB, JCVEXT     ! extra levels on the vertical boundaries
!!
!!    AUTHOR
!!    ------
!!      P. BECHTOLD       * Laboratoire d'Aerologie *
!!
!!    MODIFICATIONS
!!    -------------
!!
!!      Original    11/02/02
!!
!-------------------------------------------------------------------------------

!#include "tsmbkind.h"

!*       0.    DECLARATIONS
!              ------------

!USE YOMCST
!USE YOE_CONVPAR
!USE YOE_CONVPAREXT

IMPLICIT NONE

!*       0.1   Declarations of dummy arguments :

INTEGER,                INTENT(IN) :: KLON     ! horizontal dimension
INTEGER,                INTENT(IN) :: KLEV     ! vertical dimension

REAL(dp),DIMENSION(KLON,KLEV),INTENT(IN) :: PU     ! horizontal wind in x (m/s)
REAL(dp),DIMENSION(KLON,KLEV),INTENT(IN) :: PV     ! horizontal wind in x (m/s)
REAL(dp),DIMENSION(KLON,KLEV),INTENT(OUT):: PUC    ! convective adjusted value of u (m/s)
REAL(dp),DIMENSION(KLON,KLEV),INTENT(OUT):: PVC    ! convective adjusted value of v (m/s)

INTEGER, DIMENSION(KLON), INTENT(IN) :: KDPL   ! index for departure level
INTEGER, DIMENSION(KLON), INTENT(IN) :: KPBL   ! index for top of source layer
INTEGER, DIMENSION(KLON), INTENT(IN) :: KLCL   ! index lifting condens. level
INTEGER, DIMENSION(KLON), INTENT(IN) :: KCTL   ! index for cloud top level
INTEGER, DIMENSION(KLON), INTENT(IN) :: KLFS   ! index for level of free sink
INTEGER, DIMENSION(KLON), INTENT(IN) :: KDBL   ! index for downdraft base level

REAL(dp), DIMENSION(KLON,KLEV), INTENT(IN) :: PUMF ! updraft mass flux (kg/s)
REAL(dp), DIMENSION(KLON,KLEV), INTENT(IN) :: PUER ! updraft entrainment (kg/s)
REAL(dp), DIMENSION(KLON,KLEV), INTENT(IN) :: PUDR ! updraft detrainment (kg/s)
REAL(dp), DIMENSION(KLON,KLEV), INTENT(IN) :: PDMF ! downdraft mass flux (kg/s)
REAL(dp), DIMENSION(KLON,KLEV), INTENT(IN) :: PDER ! downdraft entrainment (kg/s)
REAL(dp), DIMENSION(KLON,KLEV), INTENT(IN) :: PDDR ! downdraft detrainment (kg/s)

REAL(dp), DIMENSION(KLON),     INTENT(IN) :: PTIMEC! convection time step
REAL(dp), DIMENSION(KLON),     INTENT(IN) :: PDXDY ! grid area (m^2)
REAL(dp), DIMENSION(KLON),     INTENT(IN) :: PMIXF ! mixed fraction at LFS
REAL(dp), DIMENSION(KLON,KLEV),INTENT(IN) :: PLMASS! mass of model layer (kg)
REAL(dp), DIMENSION(KLON,KLEV),INTENT(IN) :: PWSUB ! envir. compensating subsidence(Pa/s)
INTEGER,                   INTENT(IN) :: KFTSTEPS  ! maximum fractional time steps


!*       0.2   Declarations of local variables :

INTEGER :: IIE, IKB, IKE  ! horizontal + vertical loop bounds
INTEGER :: IKS            ! vertical dimension
INTEGER :: JI             ! horizontal loop index
INTEGER :: JK, JKP        ! vertical loop index
INTEGER :: JSTEP          ! fractional time loop index
INTEGER :: JKLD, JKLP, JKMAX ! loop index for levels

INTEGER, PARAMETER             :: IUV = 2    ! for u and v
REAL(dp), DIMENSION(KLON,KLEV)     :: ZOMG       ! compensat. subsidence (Pa/s)
REAL(dp), DIMENSION(KLON,KLEV,IUV) :: ZUUV, ZDUV ! updraft/downdraft values
REAL(dp), DIMENSION(KLON)          :: ZTIMEC     ! fractional convective time step
REAL(dp), DIMENSION(KLON,KLEV)     :: ZTIMC      ! 2D work array for ZTIMEC
REAL(dp), DIMENSION(KLON,KLEV,IUV) :: ZUVMFIN, ZUVMFOUT
                                   ! work arrays for environm. compensat. mass
REAL(dp), DIMENSION(KLON,IUV)      :: ZWORK1, ZWORK2, ZWORK3

!-------------------------------------------------------------------------------

!*       0.3   Compute loop bounds
!              -------------------

IIE    = KLON
IKB    = 1 + JCVEXB
IKS    = KLEV
IKE    = KLEV - JCVEXT
JKMAX  = MAXVAL( KCTL(:) )


!*      2.      Updraft computations
!               --------------------

ZUUV(:,:,:) = 0.0_dp

!*      2.1     Initialization  at LCL
!               ----------------------------------

DO JI = 1, IIE
    JKLD = KDPL(JI)
    JKLP = KPBL(JI)
    ZWORK1(JI,1) = 0.5_dp * ( PU(JI,JKLD) + PU(JI,JKLP) )
    ZWORK1(JI,2) = 0.5_dp * ( PV(JI,JKLD) + PV(JI,JKLP) )
ENDDO

!*      2.2     Final updraft loop
!               ------------------

DO JK = MINVAL( KDPL(:) ), JKMAX
JKP = JK + 1

     DO JI = 1, IIE
       IF ( KDPL(JI) <= JK .AND. KLCL(JI) > JK ) THEN
            ZUUV(JI,JK,1) = ZWORK1(JI,1)
            ZUUV(JI,JK,2) = ZWORK1(JI,2)
       END IF

       IF ( KLCL(JI) - 1 <= JK .AND. KCTL(JI) > JK ) THEN
                            ! instead of passive tracers equations
                            ! wind equations also include pressure term
           ZUUV(JI,JKP,1) = ( PUMF(JI,JK) * ZUUV(JI,JK,1) +                   &
                            &   PUER(JI,JKP) * PU(JI,JK) )  /                 &
                            & ( PUMF(JI,JKP) + PUDR(JI,JKP) + 1.E-7_dp ) +  &
                            &   XUVDP * ( PU(JI,JKP) - PU(JI,JK) ) 
           ZUUV(JI,JKP,2) = ( PUMF(JI,JK) * ZUUV(JI,JK,2) +                   &
                            &   PUER(JI,JKP) * PV(JI,JK) )  /                 &
                            & ( PUMF(JI,JKP) + PUDR(JI,JKP) + 1.E-7_dp ) +  &
                            &   XUVDP * ( PV(JI,JKP) - PV(JI,JK) ) 
       ENDIF
     ENDDO

ENDDO

!*      3.      Downdraft computations
!               ----------------------

ZDUV(:,:,:) = 0.0_dp

!*      3.1     Initialization at the LFS
!               -------------------------

ZWORK1(:,:) = SPREAD( PMIXF(:), DIM=2, NCOPIES=IUV )
DO JI = 1, IIE
     JK = KLFS(JI)
     ZDUV(JI,JK,1) = ZWORK1(JI,1) * PU(JI,JK) +                          &
                    &            ( 1.0_dp - ZWORK1(JI,1) ) * ZUUV(JI,JK,1)
     ZDUV(JI,JK,2) = ZWORK1(JI,2) * PV(JI,JK) +                          &
                    &            ( 1.0_dp - ZWORK1(JI,2) ) * ZUUV(JI,JK,2)
ENDDO

!*      3.2     Final downdraft loop
!               --------------------

DO JK = MAXVAL( KLFS(:) ), IKB + 1, -1
JKP = JK - 1
    DO JI = 1, IIE
      IF ( JK <= KLFS(JI) .AND. JKP >= KDBL(JI) ) THEN
       ZDUV(JI,JKP,1) = ( ZDUV(JI,JK,1) * PDMF(JI,JK) -                  &
                        &     PU(JI,JK) * PDER(JI,JKP) ) /               &
                        & ( PDMF(JI,JKP) - PDDR(JI,JKP) - 1.E-7_dp ) + & 
                        &   XUVDP * ( PU(JI,JKP) - PU(JI,JK) ) 
       ZDUV(JI,JKP,2) = ( ZDUV(JI,JK,2) * PDMF(JI,JK) -                  &
                        &     PV(JI,JK) *  PDER(JI,JKP) ) /              &
                        & ( PDMF(JI,JKP) - PDDR(JI,JKP) - 1.E-7_dp ) + &
                        &   XUVDP * ( PV(JI,JKP) - PV(JI,JK) ) 
      ENDIF
    ENDDO
ENDDO


!*      4.      Final closure (environmental) computations
!               ------------------------------------------

PUC(:,IKB:IKE) = PU(:,IKB:IKE) ! initialize adjusted envir. values
PVC(:,IKB:IKE) = PV(:,IKB:IKE) ! initialize adjusted envir. values

DO JK = IKB, IKE
   ZOMG(:,JK) = PWSUB(:,JK) * PDXDY(:) / RG ! environmental subsidence
ENDDO

ZTIMEC(:) = PTIMEC(:) / REAL( KFTSTEPS ) ! adjust  fractional time step
                                         ! to be an integer multiple of PTIMEC
WHERE ( PTIMEC(:) < 1.0_dp ) ZTIMEC(:) = 0.0_dp
ZTIMC(:,:)= SPREAD( ZTIMEC(:), DIM=2, NCOPIES=IKS )

ZUVMFIN(:,:,:)   = 0.0_dp
ZUVMFOUT(:,:,:)  = 0.0_dp


DO JSTEP = 1, KFTSTEPS ! Enter the fractional time step loop

      DO JK = IKB + 1, JKMAX
      JKP = MAX( IKB + 1, JK - 1 )
        DO JI = 1, IIE
        IF ( JK <= KCTL(JI) )  THEN
          ZWORK3(JI,1) = ZOMG(JI,JK)
          ZWORK1(JI,1) = SIGN( 1.0_dp, ZWORK3(JI,1) )
          ZWORK2(JI,1) = 0.5_dp * ( 1.0_dp + ZWORK1(JI,1) )
          ZWORK1(JI,1) = 0.5_dp * ( 1.0_dp - ZWORK1(JI,1) )
          ZUVMFIN(JI,JK,1)  = - ZWORK3(JI,1) * PUC(JI,JKP) * ZWORK1(JI,1)
          ZUVMFOUT(JI,JK,1) =   ZWORK3(JI,1) * PUC(JI,JK)  * ZWORK2(JI,1)
          ZUVMFIN(JI,JK,2)  = - ZWORK3(JI,1) * PVC(JI,JKP) * ZWORK1(JI,1)
          ZUVMFOUT(JI,JK,2) =   ZWORK3(JI,1) * PVC(JI,JK)  * ZWORK2(JI,1)
          ZUVMFIN(JI,JKP,1) = ZUVMFIN(JI,JKP,1) + ZUVMFOUT(JI,JK,1) * ZWORK2(JI,1)
          ZUVMFIN(JI,JKP,2) = ZUVMFIN(JI,JKP,2) + ZUVMFOUT(JI,JK,2) * ZWORK2(JI,1)
          ZUVMFOUT(JI,JKP,1)= ZUVMFOUT(JI,JKP,1)+ ZUVMFIN(JI,JK,1)  * ZWORK1(JI,1)
          ZUVMFOUT(JI,JKP,2)= ZUVMFOUT(JI,JKP,2)+ ZUVMFIN(JI,JK,2)  * ZWORK1(JI,1)
        END IF
        END DO
      END DO

!!$do ji=1,IIE
!!$   do jk=1, KCTL(JI)
!!$      if ( pumf(ji,jk)*2400. > plmass(ji,jk) ) print*, "too high massflux", &
!!$           ji,jk, pumf(ji,jk)*2400., plmass(ji,jk)
!!$   enddo
!!$enddo

       DO JK = IKB + 1, JKMAX
        DO JI = 1, IIE
        IF ( JK <= KCTL(JI) ) THEN
         PUC(JI,JK) = PUC(JI,JK) + ZTIMC(JI,JK) / PLMASS(JI,JK) *  (       &
                   &   ZUVMFIN(JI,JK,1) + PUDR(JI,JK) * ZUUV(JI,JK,1) +    &
                   &   PDDR(JI,JK) * ZDUV(JI,JK,1) - ZUVMFOUT(JI,JK,1) -   &
                   &   ( PUER(JI,JK) + PDER(JI,JK) ) * PU(JI,JK)    )
!!$         if (abs(PUC(ji,jk)) > 100._dp) then
!!$            print*, "WARNING, too high speed: ", ji, jk, &
!!$                 puc(ji,jk), ztimc(ji,jk), plmass(ji,jk)
!!$            print*, ZUVMFIN(JI,JK,1) + PUDR(JI,JK) * ZUUV(JI,JK,1) +    &
!!$                 &   PDDR(JI,JK) * ZDUV(JI,JK,1) - ZUVMFOUT(JI,JK,1) -   &
!!$                 &   ( PUER(JI,JK) + PDER(JI,JK) )* PU(ji,jk)
!!$            print*, ZUVMFIN(JI,JK,1), ZOMG(ji,jk), PWSUB(ji,jk),SIGN( 1.0_dp, ZOMG(JI,jk) )
!!$            print*, PUDR(JI,JK), ZUUV(JI,JK,1), PUDR(JI,JK) * ZUUV(JI,JK,1)
!!$            print*, PDDR(JI,JK), ZDUV(JI,JK,1),  PDDR(JI,JK) * ZDUV(JI,JK,1)
!!$            print*, ZUVMFOUT(JI,JK,1),  PUER(JI,JK), PDER(JI,JK),  PU(JI,JK), & 
!!$                 ( PUER(JI,JK) + PDER(JI,JK) ) * PU(JI,JK)
!!$            print*, "comp", ji-1, jk, &
!!$                 puc(ji-1,jk), ztimc(ji-1,jk), plmass(ji-1,jk)
!!$            print*, ZUVMFIN(JI-1,JK,1) + &
!!$                 PUDR(JI-1,JK) * ZUUV(JI-1,JK,1) +    &
!!$                 &   PDDR(JI-1,JK) * ZDUV(JI-1,JK,1) - ZUVMFOUT(JI-1,JK,1) -   &
!!$                 &   ( PUER(JI-1,JK) + PDER(JI-1,JK) )* PU(ji-1,jk)
!!$            print*, &
!!$                 ZUVMFIN(JI-1,JK,1), ZOMG(ji-1,jk), PWSUB(ji-1,jk),SIGN( 1.0_dp, ZOMG(JI-1,jk) )
!!$            print*, PUDR(JI-1,JK), ZUUV(JI-1,JK,1), PUDR(JI-1,JK) * ZUUV(JI-1,JK,1)
!!$            print*,  PDDR(JI-1,JK), ZDUV(JI-1,JK,1),  PDDR(JI-1,JK) * ZDUV(JI-1,JK,1)
!!$            print*, ZUVMFOUT(JI-1,JK,1),  PUER(JI-1,JK), PDER(JI-1,JK),  PU(JI-1,JK), & 
!!$                 ( PUER(JI-1,JK) + PDER(JI-1,JK) ) * PU(JI-1,JK)
!!$         END if
      
            PVC(JI,JK) = PVC(JI,JK) + ZTIMC(JI,JK) / PLMASS(JI,JK) *  (       &
                   &   ZUVMFIN(JI,JK,2) + PUDR(JI,JK) * ZUUV(JI,JK,2) +    &
                   &   PDDR(JI,JK) * ZDUV(JI,JK,2) - ZUVMFOUT(JI,JK,2) -   &
                   &   ( PUER(JI,JK) + PDER(JI,JK) ) * PV(JI,JK)    )
        END IF
        END DO
       END DO

ENDDO ! Exit the fractional time step loop


END SUBROUTINE CONVECT_UV_TRANSPORT

!============================================================================

!#######################################################################
  SUBROUTINE CONVECT_CHEM_TRANSPORT( KLON, KLEV, KCH, PCH1, PCH1C,     &
                                 & KDPL, KPBL, KLCL, KCTL, KLFS, KDBL, &
                                 & PUMF, PUER, PUDR, PDMF, PDER, PDDR, &
                                 & PTIMEC, PDXDY, PMIXF, PLMASS, PWSUB,&
                                 & KFTSTEPS )
!#######################################################################

!!**** Compute  modified chemical tracer values due to convective event
!!
!!
!!    PURPOSE
!!    -------
!!      The purpose of this routine is to determine the final adjusted
!!      environmental values of the chemical tracers
!!      The final convective tendencies can then be evaluated in the main
!!      routine DEEP_CONVECT by (PCH1C-PCH1)/PTIMEC
!!
!!
!!**  METHOD
!!    ------
!!      Identical to the computation of the conservative variables in the
!!      main deep convection code
!!
!!    EXTERNAL
!!    --------
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!      Module YOMCST
!!          RG                 ! gravity constant
!!
!!     Module YOE_CONVPAREXT
!!          JCVEXB, JCVEXT     ! extra levels on the vertical boundaries
!!
!!    AUTHOR
!!    ------
!!      P. BECHTOLD       * Laboratoire d'Aerologie *
!!
!!    MODIFICATIONS
!!    -------------
!!
!!      Original    11/12/97
!!
!-------------------------------------------------------------------------------

!#include "tsmbkind.h"

!*       0.    DECLARATIONS
!              ------------

!USE YOMCST
!USE YOE_CONVPAREXT

IMPLICIT NONE

!*       0.1   Declarations of dummy arguments :

INTEGER,                INTENT(IN) :: KLON     ! horizontal dimension
INTEGER,                INTENT(IN) :: KLEV     ! vertical dimension
INTEGER,                INTENT(IN) :: KCH      ! number of passive tracers

REAL(dp),DIMENSION(KLON,KLEV,KCH),INTENT(IN) :: PCH1 ! grid scale tracer concentr.
REAL(dp),DIMENSION(KLON,KLEV,KCH),INTENT(OUT):: PCH1C! conv adjusted tracer concntr.

INTEGER, DIMENSION(KLON), INTENT(IN) :: KDPL   ! index for departure level
INTEGER, DIMENSION(KLON), INTENT(IN) :: KPBL   ! index for top of source layer
INTEGER, DIMENSION(KLON), INTENT(IN) :: KLCL   ! index lifting condens. level
INTEGER, DIMENSION(KLON), INTENT(IN) :: KCTL   ! index for cloud top level
INTEGER, DIMENSION(KLON), INTENT(IN) :: KLFS   ! index for level of free sink
INTEGER, DIMENSION(KLON), INTENT(IN) :: KDBL   ! index for downdraft base level

REAL(dp), DIMENSION(KLON,KLEV), INTENT(IN) :: PUMF ! updraft mass flux (kg/s)
REAL(dp), DIMENSION(KLON,KLEV), INTENT(IN) :: PUER ! updraft entrainment (kg/s)
REAL(dp), DIMENSION(KLON,KLEV), INTENT(IN) :: PUDR ! updraft detrainment (kg/s)
REAL(dp), DIMENSION(KLON,KLEV), INTENT(IN) :: PDMF ! downdraft mass flux (kg/s)
REAL(dp), DIMENSION(KLON,KLEV), INTENT(IN) :: PDER ! downdraft entrainment (kg/s)
REAL(dp), DIMENSION(KLON,KLEV), INTENT(IN) :: PDDR ! downdraft detrainment (kg/s)

REAL(dp), DIMENSION(KLON),     INTENT(IN) :: PTIMEC! convection time step
REAL(dp), DIMENSION(KLON),     INTENT(IN) :: PDXDY ! grid area (m^2)
REAL(dp), DIMENSION(KLON),     INTENT(IN) :: PMIXF ! mixed fraction at LFS
REAL(dp), DIMENSION(KLON,KLEV),INTENT(IN) :: PLMASS! mass of model layer (kg)
REAL(dp), DIMENSION(KLON,KLEV),INTENT(IN) :: PWSUB ! envir. compensating subsidence(Pa/s)
INTEGER,                   INTENT(IN) :: KFTSTEPS  ! maximum fractional time steps


!*       0.2   Declarations of local variables :

INTEGER :: INCH1          ! number of chemical tracers
INTEGER :: IIE, IKB, IKE  ! horizontal + vertical loop bounds
INTEGER :: IKS            ! vertical dimension
INTEGER :: JI             ! horizontal loop index
INTEGER :: JK, JKP        ! vertical loop index
INTEGER :: JN             ! chemical tracer loop index
INTEGER :: JSTEP          ! fractional time loop index
INTEGER :: JKLD, JKLP, JKMAX ! loop index for levels

REAL(dp), DIMENSION(KLON,KLEV)     :: ZOMG ! compensat. subsidence (Pa/s)
REAL(dp), DIMENSION(KLON,KLEV,KCH) :: ZUCH1, ZDCH1 ! updraft/downdraft values
REAL(dp), DIMENSION(KLON)          :: ZTIMEC  ! fractional convective time step
REAL(dp), DIMENSION(KLON,KLEV)     :: ZTIMC! 2D work array for ZTIMEC
REAL(dp), DIMENSION(KLON,KLEV,KCH) :: ZCH1MFIN, ZCH1MFOUT
                                   ! work arrays for environm. compensat. mass
REAL(dp), DIMENSION(KLON,KCH)      :: ZWORK1, ZWORK2, ZWORK3

!-------------------------------------------------------------------------------

!*       0.3   Compute loop bounds
!              -------------------

INCH1  = KCH
IIE    = KLON
IKB    = 1 + JCVEXB
IKS    = KLEV
IKE    = KLEV - JCVEXT
JKMAX  = MAXVAL( KCTL(:) )

!*      2.      Updraft computations
!               --------------------

ZUCH1(:,:,:) = 0.0_dp

!*      2.1     Initialization  at LCL
!               ----------------------------------

DO JI = 1, IIE
    JKLD = KDPL(JI)
    JKLP = KPBL(JI)
    ZWORK1(JI,:) = 0.5_dp * ( PCH1(JI,JKLD,:) + PCH1(JI,JKLP,:) )
ENDDO

!*      2.2     Final updraft loop
!               ------------------

DO JK = MINVAL( KDPL(:) ), JKMAX
JKP = JK + 1

    DO JN = 1, INCH1
     DO JI = 1, IIE
       IF ( KDPL(JI) <= JK .AND. KLCL(JI) > JK )                             &
          & ZUCH1(JI,JK,JN) = ZWORK1(JI,JN)

       IF ( KLCL(JI) - 1 <= JK .AND. KCTL(JI) > JK ) THEN
                       !if you have reactive i.e. non-passive tracers
                       ! add the corresponding sink term in the following equation
           ZUCH1(JI,JKP,JN) = ( PUMF(JI,JK) * ZUCH1(JI,JK,JN) +              &
                            &   PUER(JI,JKP) * PCH1(JI,JK,JN) )  /           &
                            & ( PUMF(JI,JKP) + PUDR(JI,JKP) + 1.E-7_dp )
       ENDIF
     ENDDO
   ENDDO

ENDDO

!*      3.      Downdraft computations
!               ----------------------

ZDCH1(:,:,:) = 0.0_dp

!*      3.1     Initialization at the LFS
!               -------------------------

ZWORK1(:,:) = SPREAD( PMIXF(:), DIM=2, NCOPIES=INCH1 )
DO JI = 1, IIE
     JK = KLFS(JI)
     ZDCH1(JI,JK,:) = ZWORK1(JI,:) * PCH1(JI,JK,:) +                          &
                    &                  ( 1.0_dp - ZWORK1(JI,:) ) * ZUCH1(JI,JK,:)
ENDDO

!*      3.2     Final downdraft loop
!               --------------------

DO JK = MAXVAL( KLFS(:) ), IKB + 1, -1
JKP = JK - 1
    DO JN = 1, INCH1
    DO JI = 1, IIE
      IF ( JK <= KLFS(JI) .AND. JKP >= KDBL(JI) ) THEN
       ZDCH1(JI,JKP,JN) = ( ZDCH1(JI,JK,JN) * PDMF(JI,JK) -              &
                        &   PCH1(JI,JK,JN) *  PDER(JI,JKP) ) /           &
                        & ( PDMF(JI,JKP) - PDDR(JI,JKP) - 1.E-7_dp )
      ENDIF
    ENDDO
    ENDDO
ENDDO


!*      4.      Final closure (environmental) computations
!               ------------------------------------------

PCH1C(:,IKB:IKE,:) = PCH1(:,IKB:IKE,:) ! initialize adjusted envir. values

DO JK = IKB, IKE
   ZOMG(:,JK) = PWSUB(:,JK) * PDXDY(:) / RG ! environmental subsidence
ENDDO

ZTIMEC(:) = PTIMEC(:) / REAL( KFTSTEPS ) ! adjust  fractional time step
                                           ! to be an integer multiple of PTIMEC
WHERE ( PTIMEC(:) < 1.0_dp ) ZTIMEC(:) = 0.0_dp
ZTIMC(:,:)= SPREAD( ZTIMEC(:), DIM=2, NCOPIES=IKS )

ZCH1MFIN(:,:,:)   = 0.0_dp
ZCH1MFOUT(:,:,:)  = 0.0_dp


DO JSTEP = 1, KFTSTEPS ! Enter the fractional time step loop

      DO JK = IKB + 1, JKMAX
         JKP = MAX( IKB + 1, JK - 1 )
         DO JI = 1, IIE
          ZWORK3(JI,1) = ZOMG(JI,JK)
          ZWORK1(JI,1) = SIGN( 1.0_dp, ZWORK3(JI,1) )
          ZWORK2(JI,1) = 0.5_dp * ( 1. + ZWORK1(JI,1) )
          ZWORK1(JI,1) = 0.5_dp * ( 1. - ZWORK1(JI,1) )
          ZCH1MFIN(JI,JK,:)  = - ZWORK3(JI,1) * PCH1C(JI,JKP,:) * ZWORK1(JI,1)
          ZCH1MFOUT(JI,JK,:) =   ZWORK3(JI,1) * PCH1C(JI,JK,:)  * ZWORK2(JI,1)
          ZCH1MFIN(JI,JKP,:) = ZCH1MFIN(JI,JKP,:) + ZCH1MFOUT(JI,JK,:) * ZWORK2(JI,1)
          ZCH1MFOUT(JI,JKP,:)= ZCH1MFOUT(JI,JKP,:) + ZCH1MFIN(JI,JK,:) * ZWORK1(JI,1)
         END DO
      END DO
!
       DO JK = IKB + 1, JKMAX
       DO JN = 1, INCH1
       DO JI = 1, IIE
         PCH1C(JI,JK,JN) = PCH1C(JI,JK,JN) + ZTIMC(JI,JK) / PLMASS(JI,JK) *  (    &
                      ZCH1MFIN(JI,JK,JN) + PUDR(JI,JK) * ZUCH1(JI,JK,JN) +        &
                      PDDR(JI,JK) * ZDCH1(JI,JK,JN) - ZCH1MFOUT(JI,JK,JN) -       &
                      ( PUER(JI,JK) + PDER(JI,JK) ) * PCH1(JI,JK,JN)    )
      !  PCH1C(JI,JK,JN) = MAX( 0., PCH1C(JI,JK,JN) )
       END DO
       END DO
      END DO

ENDDO ! final values

END SUBROUTINE CONVECT_CHEM_TRANSPORT
!============================================================================

!============================================================================
!============================================================================
 END MODULE MESSY_CONVECT_BECHTOLD_DEEP

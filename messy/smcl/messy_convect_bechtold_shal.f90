MODULE MESSY_CONVECT_BECHTOLD_SHAL

USE MESSY_CONVECT_BECHTOLD_PARAM,  ONLY: dp,                     &
                                         RG, RV, RD, RCPD, RCPV, &
                                         JCVEXB, JCVEXT
IMPLICIT NONE
PRIVATE
SAVE
PUBLIC:: CONVECT_TRIGGER_SHAL, CONVECT_UPDRAFT_SHAL, CONVECT_CLOSURE_SHAL


CONTAINS

!=============================================================================
!######################################################################
 SUBROUTINE CONVECT_TRIGGER_SHAL(  KLON, KLEV,                        &
                             &  PPRES, PTH, PTHV, PTHES,              &
                             &  PRV, PW, PZ, PDXDY, KCOUNT,           &
                             &  PTHLCL, PTLCL, PRVLCL, PWLCL, PZLCL,  &
                             &  PTHVELCL, KLCL, KDPL, KPBL, OTRIG     )
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
!!          XCDEPTH            ! minimum necessary cloud depth
!!          XCDEPTH_D          ! maximum allowed cloud depth
!!          XDTPERT            ! add small Temp peturbation
!!          XNHGAM             ! coefficient for buoyancy term in w eq.
!!                             ! accounting for nh-pressure
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
!USE YOE_CONVPAR_SHAL
!USE YOE_CONVPAREXT
USE MESSY_CONVECT_BECHTOLD_PARAM, ONLY: RATM, RTT, RGAMW, RBETW, &
                                        XZLCL_S, XZPBL_S, XDTPERT_S, XNHGAM_S, XCDEPTH_S
USE MESSY_CONVECT_BECHTOLD_DEEP,  ONLY: CONVECT_SATMIXRATIO

IMPLICIT NONE

!*       0.1   Declarations of dummy arguments :

INTEGER, INTENT(IN)                   :: KLON      ! horizontal loop index
INTEGER, INTENT(IN)                   :: KLEV      ! vertical loop index
REAL(dp), DIMENSION(KLON),     INTENT(IN) :: PDXDY     ! grid area
REAL(dp), DIMENSION(KLON,KLEV),INTENT(IN) :: PTH, PTHV ! theta, theta_v
REAL(dp), DIMENSION(KLON,KLEV),INTENT(IN) :: PTHES     ! envir. satur. theta_e
REAL(dp), DIMENSION(KLON,KLEV),INTENT(IN) :: PRV       ! vapor mixing ratio
REAL(dp), DIMENSION(KLON,KLEV),INTENT(IN) :: PPRES     ! pressure
REAL(dp), DIMENSION(KLON,KLEV),INTENT(IN) :: PZ        ! height of grid point (m)
REAL(dp), DIMENSION(KLON,KLEV),INTENT(IN) :: PW        ! vertical velocity
INTEGER, DIMENSION(KLON),  INTENT(IN) :: KCOUNT    ! convective counter for already
                                                     ! active deep convection points

REAL(dp), DIMENSION(KLON),     INTENT(OUT):: PTHLCL    ! theta at LCL
REAL(dp), DIMENSION(KLON),     INTENT(OUT):: PTLCL     ! temp. at LCL
REAL(dp), DIMENSION(KLON),     INTENT(OUT):: PRVLCL    ! vapor mixing ratio at  LCL
REAL(dp), DIMENSION(KLON),     INTENT(OUT):: PWLCL     ! parcel velocity at  LCL
REAL(dp), DIMENSION(KLON),     INTENT(OUT):: PZLCL     ! height at LCL (m)
REAL(dp), DIMENSION(KLON),     INTENT(OUT):: PTHVELCL  ! environm. theta_v at LCL (K)
LOGICAL,   DIMENSION(KLON),  INTENT(OUT):: OTRIG     ! logical mask for convection
INTEGER, DIMENSION(KLON),  INTENT(INOUT):: KLCL    ! contains vert. index of LCL
INTEGER, DIMENSION(KLON),  INTENT(INOUT):: KDPL    ! contains vert. index of DPL
INTEGER, DIMENSION(KLON),  INTENT(INOUT):: KPBL    ! contains index of source layer top

!*       0.2   Declarations of local variables :

INTEGER :: JKK, JK, JKP, JKM, JKDL, JL, JKT, JT! vertical loop index
INTEGER :: JI                                  ! horizontal loop index
INTEGER :: IIE, IKB, IKE                       ! horizontal + vertical loop bounds
REAL(dp)    :: ZEPS, ZEPSA                         ! R_d / R_v, R_v / R_d
REAL(dp)    :: ZCPORD, ZRDOCP                      ! C_pd / R_d, R_d / C_pd

REAL(dp), DIMENSION(KLON) :: ZTHLCL, ZTLCL, ZRVLCL, & ! locals for PTHLCL,PTLCL
                        &  ZWLCL,  ZZLCL, ZTHVELCL  ! PRVLCL, ....
INTEGER, DIMENSION(KLON) :: IDPL, IPBL, ILCL      ! locals for KDPL, ...
REAL(dp), DIMENSION(KLON) :: ZPLCL    ! pressure at LCL
REAL(dp), DIMENSION(KLON) :: ZZDPL    ! height of DPL
REAL(dp), DIMENSION(KLON) :: ZTHVLCL  ! theta_v at LCL = mixed layer value
REAL(dp), DIMENSION(KLON) :: ZTMIX    ! mixed layer temperature
REAL(dp), DIMENSION(KLON) :: ZEVMIX   ! mixed layer water vapor pressure
REAL(dp), DIMENSION(KLON) :: ZDPTHMIX, ZPRESMIX ! mixed layer depth and pressure
REAL(dp), DIMENSION(KLON) :: ZCAPE    ! convective available energy (m^2/s^2/g)
REAL(dp), DIMENSION(KLON) :: ZCAP     ! pseudo for CAPE
REAL(dp), DIMENSION(KLON) :: ZTHEUL   ! updraft equiv. pot. temperature (K)
REAL(dp), DIMENSION(KLON) :: ZLV, ZCPH! specific heats of vaporisation, dry air
REAL(dp), DIMENSION(KLON) :: ZDP      ! pressure between LCL and model layer
REAL(dp), DIMENSION(KLON) :: ZTOP     ! estimated cloud top (m)
!INTEGER, DIMENSION(KLON) :: ITOP  ! work array to store highest test layer
REAL(dp), DIMENSION(KLON) :: ZWORK1, ZWORK2, ZWORK3    ! work arrays
LOGICAL, DIMENSION(KLON) :: GTRIG, GTRIG2          ! local arrays for OTRIG
LOGICAL, DIMENSION(KLON) :: GWORK1                 ! work array


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



!       1.     Determine highest necessary loop test layer
!              -------------------------------------------

JT = IKE - 2
DO JK = IKB + 1, IKE - 2
 ! DO JI = 1, IIE
 !    IF ( PZ(JI,JK) - PZ(JI,IKB) <= XZLCL_S ) ITOP(JI) = JK
 ! ENDDO
   IF ( PZ(1,JK) - PZ(1,IKB) < 5.E3_dp ) JT = JK
ENDDO


!*       2.     Enter loop for convection test
!               ------------------------------

JKP  = MINVAL( IDPL(:) ) + 1
JKT  = JT

JKT=JKP ! do not loop anymore and only allow departure close to surface
DO JKK = JKP, JKT

     GWORK1(:) = ZZDPL(:) - PZ(:,IKB) < XZLCL_S
          ! we exit the trigger test when the center of the mixed layer is more
          ! than 1500 m  above soil level.
     WHERE ( GWORK1(:) )
        ZDPTHMIX(:) = 0.0_dp
        ZPRESMIX(:) = 0.0_dp
        ZTHLCL(:)   = 0.0_dp
        ZRVLCL(:)   = 0.0_dp
        ZZDPL(:)    = PZ(:,JKK)
        IDPL(:)     = JKK
     END WHERE


!*       3.     Construct a mixed layer of at least 50 hPa (XZPBL_S)
!               ------------------------------------------

     DO JK = JKK, IKE - 1
       JKM = JK + 1
       DO JI = 1, IIE
         IF ( GWORK1(JI) .AND. ZDPTHMIX(JI) < XZPBL_S ) THEN
            IPBL(JI)     = JK
            ZWORK1(JI)   = PPRES(JI,JK) - PPRES(JI,JKM)
            ZDPTHMIX(JI) = ZDPTHMIX(JI) + ZWORK1(JI)
            ZPRESMIX(JI) = ZPRESMIX(JI) + PPRES(JI,JK) * ZWORK1(JI)
            ZTHLCL(JI)   = ZTHLCL(JI)   + PTH(JI,JK)   * ZWORK1(JI)
            ZRVLCL(JI)   = ZRVLCL(JI)   + PRV(JI,JK)   * ZWORK1(JI)
         ENDIF
       ENDDO
        IF ( MINVAL ( ZDPTHMIX(:) ) >= XZPBL_S ) EXIT
     ENDDO


     WHERE ( GWORK1(:) )

        ZPRESMIX(:) = ZPRESMIX(:) / ZDPTHMIX(:)
        ZTHLCL(:)   = ZTHLCL(:)   / ZDPTHMIX(:) + XDTPERT_S ! add small Temp Perturb.
        ZRVLCL(:)   = ZRVLCL(:)   / ZDPTHMIX(:)
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

!            !  normalize w grid scale to a 25 km refer. grid
!    DO JI = 1, IIE
!       JK  = ILCL(JI)
!       JKM = JK - 1
!       ZWORK1(JI) =  ( PW(JI,JKM)  + ( PW(JI,JK) - PW(JI,JKM) ) * ZDP(JI) )  &
!                  &       * SQRT( PDXDY(JI) / XA25 )
!                  &      - 0.02_dp  * ZZLCL(JI) / XZLCL_S ! avoid spurious convection
!    END DO
!            ! compute sign of normalized grid scale w
!       ZWORK2(:) = SIGN( 1.0_dp, ZWORK1(:) )
!       ZWORK1(:) = XWTRIG * ZWORK2(:) * ABS( ZWORK1(:) ) ** 0.333_dp       &
!                 &        * ( RATM / ZPLCL(:) ) ** ZRDOCP

!*       6.2    Compute parcel vertical velocity at LCL
!               ---------------------------------------

!    DO JI = 1, IIE
!       JKDL = IDPL(JI)
!       ZWORK3(JI) = RG * ZWORK1(JI) * ( ZZLCL(JI) - PZ(JI,JKDL) )       &
!                  &   / ( PTHV(JI,JKDL) + ZTHVELCL(JI) )
!    END DO
!    WHERE( GWORK1(:) )
!      ZWLCL(:)  = 1.0_dp + .5_dp * ZWORK2(:) * SQRT( ABS( ZWORK3(:) ) )
!      GTRIG(:)  = ZTHVLCL(:) - ZTHVELCL(:) + ZWORK1(:) > 0.0_dp .AND.       &
!                & ZWLCL(:) > 0.0_dp
!    END WHERE
!
     ZWLCL(:) = 1.0_dp


!*       6.3    Look for parcel that produces sufficient cloud depth.
!               The cloud top is estimated as the level where the CAPE
!               is smaller  than a given value (based on vertical velocity eq.)
!               --------------------------------------------------------------

     ZTHEUL(:) = ZTLCL(:) * ( ZTHLCL(:) / ZTLCL(:) ) **                          &
               &            ( 1.0_dp - 0.28_dp * ZRVLCL(:) )                    &
               &          * EXP( ( 3374.6525_dp / ZTLCL(:) - 2.5403_dp )     &
               &          *      ZRVLCL(:) * ( 1.0_dp + 0.81_dp * ZRVLCL(:) ) )

     ZCAPE(:) = 0.0_dp
     ZCAP(:)  = 0.0_dp
     ZTOP(:)  = 0.0_dp
     ZWORK3(:)= 0.0_dp
     JKM = MINVAL( ILCL(:) )
     DO JL = JKM, JT
        JK = JL + 1
        DO JI = 1, IIE
           ZWORK1(JI) = ( 2.0_dp * ZTHEUL(JI) /                                &
           & ( PTHES(JI,JK) + PTHES(JI,JL) ) - 1.0_dp ) * ( PZ(JI,JK) - PZ(JI,JL) )
           IF ( JL < ILCL(JI) ) ZWORK1(JI) = 0.0_dp
           ZCAP(JI)   = ZCAP(JI) + ZWORK1(JI)
           ZCAPE(JI)  = ZCAPE(JI) + RG * MAX( 0.0_dp, ZWORK1(JI) )
           ZWORK2(JI) = XNHGAM_S * RG * ZCAP(JI) + 1.05_dp * ZWLCL(JI) * ZWLCL(JI)
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


     ZWORK2(:) = ZTOP(:) - ZZLCL(:)
   ! WHERE( ZWORK2(:)   >=  XCDEPTH_S  .AND. ZWORK2(:) < XCDEPTH_D_S .AND. GTRIG2(:) &
     WHERE( ZWORK2(:)   >=  XCDEPTH_S  .AND. GTRIG2(:) .AND. KCOUNT(:) == 0 &
          &  .AND. ZCAPE(:) > 10._dp )
        GTRIG2(:)   = .FALSE.
        OTRIG(:)    = .TRUE.
      ! OTRIG(:)    = GTRIG(:)     ! we  select the first departure level
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


END SUBROUTINE CONVECT_TRIGGER_SHAL
!==============================================================================================

!
 SUBROUTINE CONVECT_UPDRAFT_SHAL( KLON, KLEV,                                &
                           & KICE, PPRES, PDPRES, PZ, PTHL, PTHV, PTHES, PRW,&
                           & PTHLCL, PTLCL, PRVLCL, PWLCL, PZLCL, PTHVELCL,  &
                           & PMFLCL, OTRIG, KLCL, KDPL, KPBL,                &
                           & PUMF, PUER, PUDR, PUTHL, PUTHV, PURW,           &
                           & PURC, PURI, PCAPE, KCTL, KETL )
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
!!      Module YOE_CONVPAR_SHAL
!!          XA25               ! reference grid area
!!          XCRAD              ! cloud radius
!!          XCDEPTH            ! minimum necessary cloud depth
!!          XCDEPTH_D          ! maximum allowed   cloud depth
!!          XENTR              ! entrainment constant
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
!USE YOE_CONVPAR_SHAL
!USE YOE_CONVPAREXT
USE MESSY_CONVECT_BECHTOLD_PARAM,  ONLY: RATM,                             &
                                         XNHGAM_S, XCRAD_S, XENTR_S, XCDEPTH_D_S, XCDEPTH_S
USE MESSY_CONVECT_BECHTOLD_DEEP,   ONLY: CONVECT_CONDENS, CONVECT_MIXING_FUNCT
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
REAL(dp), DIMENSION(KLON),     INTENT(IN) :: PTHLCL ! theta at LCL
REAL(dp), DIMENSION(KLON),     INTENT(IN) :: PTLCL  ! temp. at LCL
REAL(dp), DIMENSION(KLON),     INTENT(IN) :: PRVLCL ! vapor mixing ratio at  LCL
REAL(dp), DIMENSION(KLON),     INTENT(IN) :: PWLCL  ! parcel velocity at LCL (m/s)
REAL(dp), DIMENSION(KLON),     INTENT(IN) :: PMFLCL ! cloud  base unit mass flux
                                                  ! (kg/s)
REAL(dp), DIMENSION(KLON),     INTENT(IN) :: PZLCL  ! height at LCL (m)
REAL(dp), DIMENSION(KLON),     INTENT(IN) :: PTHVELCL  ! environm. theta_v at LCL (K)
LOGICAL, DIMENSION(KLON),  INTENT(INOUT):: OTRIG  ! logical mask for convection
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
REAL(dp), DIMENSION(KLON),     INTENT(OUT):: PCAPE  ! available potent. energy

!*       0.2   Declarations of local variables :

INTEGER :: IIE, IKB, IKE  ! horizontal and vertical loop bounds
INTEGER :: JI             ! horizontal loop index
INTEGER :: JK, JKP, JKM, JK1, JK2, JKMIN   ! vertical loop index
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
LOGICAL, DIMENSION(KLON) :: GWORK1, GWORK2, GWORK4, GWORK5
                                              ! work arrays
LOGICAL, DIMENSION(KLON,KLEV) :: GWORK6       ! work array


!-------------------------------------------------------------------------------

!        0.3   Set loop bounds
!              ---------------

IKB = 1 + JCVEXB
IKE = KLEV - JCVEXT
IIE = KLON


!*       1.     Initialize updraft properties and local variables
!               -------------------------------------------------

ZEPSA      = RV / RD
ZCVOCD     = RCPV / RCPD
ZCPORD     = RCPD / RD
ZRDOCP     = RD / RCPD

PUMF(:,:)  = 0.0_dp
PUER(:,:)  = 0.0_dp
PUDR(:,:)  = 0.0_dp
PUTHL(:,:) = 0.0_dp
PUTHV(:,:) = 0.0_dp
PURW(:,:)  = 0.0_dp
PURC(:,:)  = 0.0_dp
PURI(:,:)  = 0.0_dp
ZUW1(:)    = PWLCL(:) * PWLCL(:)
ZUW2(:)    = 0.0_dp
ZE1(:)     = 0.0_dp
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

ZTHEU1(:) = PTLCL(:) * ( PTHLCL(:) / PTLCL(:) ) ** ( 1._dp - 0.28_dp * PRVLCL(:) )  &
          &          * EXP( ( 3374.6525_dp / PTLCL(:) - 2.5403_dp )       &
          &          * PRVLCL(:) * ( 1.0_dp + 0.81_dp * PRVLCL(:) ) )


ZWORK1(:) = ( RCPD + PRVLCL(:) * RCPV ) * PTLCL(:)                            &
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
        PUTHV(JI,JK) = PTHLCL(JI) * ( 1.0_dp + ZEPSA * PRVLCL(JI) ) /         &
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

    ZWORK1(:) = PURC(:,JK)
    ZWORK2(:) = PURI(:,JK)
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
    ZWORK5(:) = 2.0_dp * ZUW1(:) * PUER(:,JK) / MAX( .1_dp, PUMF(:,JK) )
    ZUW2(:)   = ZUW1(:) + ZWORK3(:) * XNHGAM_S * RG *        &
              &   ( ( PUTHV(:,JK) + PUTHV(:,JKP) ) /       &
              &   ( ZWORK4(:) + PTHV(:,JKP) ) - 1.0_dp )    & ! buoyancy term
              & - ZWORK5(:)                                  ! entrainment term


!*       6.     Update total precipitation: dr_r=(r_c+r_i)*exp(-rate*dz)
!               --------------------------------------------------------

!                    compute level mean vertical velocity
    ZWORK2(:)   = 0.5_dp *                                                    &
                &      ( SQRT( MAX( 1.E-2_dp, ZUW2(:) ) ) +                 &
                &        SQRT( MAX( 1.E-2_dp, ZUW1(:) ) ) )


!*       7.     Update r_c, r_i, enthalpy, r_w  for precipitation
!               -------------------------------------------------------

    PURW(:,JKP)  = PURW(:,JK)
    PURC(:,JKP)  = PURC(:,JKP)
    PURI(:,JKP)  = PURI(:,JKP)
    PUTHL(:,JKP) = PUTHL(:,JK)

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
    ZWORK3(:) = ZUT(:) * ZPI(:) * ( 1.0_dp + ZEPSA * (                          &
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

  ze2=min(zd2,max(.3_dp,ze2))

!       Note: routine MIXING_FUNCT returns fractional entrainm/detrainm. rates

! ZWORK1(:) = XENTR_S * PMFLCL(:) * PDPRES(:,JKP) / XCRAD ! rate of env. inflow
!*MOD
  ZWORK3(:) = PZ(:,JKP) - PZ(:,JK) * ZWORK6(:) -         &
              &      ( 1.0_dp - ZWORK6(:) ) * PZLCL(:)        ! level thickness
  ZWORK1(:) = XENTR_S * RG / XCRAD_S * PUMF(:,JK) * ZWORK3(:)
! ZWORK1(:) = XENTR_S * pumf(:,jk) * PDPRES(:,JKP) / XCRAD_S ! rate of env. inflow
!*MOD
  ZWORK2(:) = 0.0_dp
  WHERE ( GWORK1(:) ) ZWORK2(:) = 1.0_dp
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
   & ( PTHES(:,JKP) - PTHES(:,JK) ) / ( PZ(:,JKP) - PZ(:,JK) ) *          &
   & ( PZLCL(:) - PZ(:,JK) ) ! linear interpolation for theta_es at LCL
                            ! ( this is only done for model level just above LCL

    ZWORK1(:) = MAX( 0.01_dp, PUER(:,JK) / MAX( .1_dp, PUMF(:,JK) ) )
    ZTHEU2(:) = ( 1.0_dp - ZWORK1(:) ) * ZTHEU1(:) + ZWORK1(:) * PTHES(:,JK)
    ZWORK1(:) = ( ZTHEU1(:) + ZTHEU2(:) ) / ( ZWORK2(:) + PTHES(:,JKP) ) - 1.0_dp
    ZTHEU1(:) = ZTHEU2(:)
    PCAPE(:)  = PCAPE(:) + RG * ZWORK3(:) * MAX( 0.0_dp, ZWORK1(:) )


!*       10.   Compute final values of updraft mass flux, enthalpy, r_w
!              at level k+1
!              --------------------------------------------------------

    PUMF(:,JKP)  = PUMF(:,JK) - PUDR(:,JKP) + PUER(:,JKP)
    PUMF(:,JKP)  = MAX( PUMF(:,JKP), 0.1_dp )
    PUTHL(:,JKP) = ( PUMF(:,JK)  * PUTHL(:,JK) +                             &
                 &   PUER(:,JKP) * PTHL(:,JK) - PUDR(:,JKP) * PUTHL(:,JK) )  &
                 &  / PUMF(:,JKP)
    PURW(:,JKP)  = ( PUMF(:,JK)  * PURW(:,JK) +                              &
                 &   PUER(:,JKP) * PRW(:,JK) - PUDR(:,JKP) * PURW(:,JK) )    &
                 &  / PUMF(:,JKP)


    ZE1(:) = ZE2(:) ! update fractional entrainment/detrainment
    ZD1(:) = ZD2(:)

  END WHERE

ENDDO

!*       12.1    Set OTRIG to False if cloud thickness < 0.5km
!                or > 3km (deep convection) or CAPE < 1
!                ------------------------------------------------

    DO JI = 1, IIE
          JK  = KCTL(JI)
          ZWORK1(JI) = PZ(JI,JK) - PZLCL(JI)
          OTRIG(JI) = ZWORK1(JI) >= XCDEPTH_S  .AND. ZWORK1(JI) < XCDEPTH_D_S &
                    & .AND. PCAPE(JI) > 1.0_dp
    ENDDO
    WHERE( .NOT. OTRIG(:) )
          KCTL(:) = IKB
    END WHERE

KCTL(:) = MIN( IKE-2, KCTL(:) )
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
        PUDR(JI,JK)  = PDPRES(JI,JK) * ZWORK1(JI)
        PUMF(JI,JK)  = PUMF(JI,JKP) - PUDR(JI,JK)
    ENDIF
    ENDDO
ENDDO

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
DO JK = IKB+1, JKP
   DO JI = 1, IIE
   IF ( JK >= IKB+1  .AND. JK <= IWORK(JI) ) THEN
       PUER(JI,JK) = PMFLCL(JI) * (PPRES(JI,JK-1)-PPRES(JI,JK)) / ( ZWORK2(JI) + 0.1 )
       PUMF(JI,JK) = PUMF(JI,JK-1) + PUER(JI,JK)
   ENDIF
   ENDDO
ENDDO
DO JI = 1, IIE
   JK = KLCL(JI)
    PUDR(JI,JK) =- PUMF(JI,JK) + PUMF(JI,JK-1) + PUER(JI,JK)
END DO

!*       13.   If cloud thickness is smaller than  .5 km or > 3 km
!              no shallow convection is allowed
!              Nota: For technical reasons, we stop the convection
!                    computations in this case and do not go back to
!                    TRIGGER_FUNCT to look for the next unstable LCL
!                    which could produce a thicker cloud.
!              ---------------------------------------------------

GWORK6(:,:) = SPREAD( OTRIG(:), DIM=2, NCOPIES=KLEV )
WHERE ( .NOT. GWORK6(:,:) )
    PUMF(:,:)  = 0.0_dp
    PUDR(:,:)  = 0.0_dp
    PUER(:,:)  = 0.0_dp
    PUTHL(:,:) = PTHL(:,:)
    PURW(:,:)  = PRW(:,:)
    PURC(:,:)  = 0.0_dp
    PURI(:,:)  = 0.0_dp
END WHERE

END SUBROUTINE CONVECT_UPDRAFT_SHAL

!==============================================================================================

!#######################################################################
 SUBROUTINE CONVECT_CLOSURE_SHAL( KLON, KLEV,                          &
                       &   PPRES, PDPRES, PZ, PDXDY, PLMASS,           &
                       &   PTHL, PTH, PRW, PRC, PRI, OTRIG1,           &
                       &   PTHC, PRWC, PRCC, PRIC, PWSUB,              &
                       &   KLCL, KDPL, KPBL, KCTL,                     &
                       &   PUMF, PUER, PUDR, PUTHL, PURW,              &
                       &   PURC, PURI, PCAPE, PTIMEC, KFTSTEPS         )
!#######################################################################

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
!!    CONVECT_CLOSURE_ADJUST_SHAL
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
!!      Module YOE_CONVPAR_SHAL
!!          XA25               ! reference grid area
!!          XSTABT             ! stability factor in time integration
!!          XSTABC             ! stability factor in CAPE adjustment
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
!!   Peter Bechtold 15/11/96 change for enthalpie, r_c + r_i tendencies
!!      Tony Dore   14/10/96 Initialise local variables
!-------------------------------------------------------------------------------

!#include "tsmbkind.h"

!*       0.    DECLARATIONS
!              ------------

!USE YOMCST
!USE YOE_CONVPAR_SHAL
!USE YOE_CONVPAREXT
USE MESSY_CONVECT_BECHTOLD_PARAM,   ONLY: RATM, RTT, RCW, RLVTT, RCS,     &
                                          XSTABT_S, XSTABC_S
USE MESSY_CONVECT_BECHTOLD_DEEP,    ONLY: CONVECT_CLOSURE_THRVLCL, CONVECT_SATMIXRATIO
IMPLICIT NONE

!*       0.1   Declarations of dummy arguments :

INTEGER,                   INTENT(IN) :: KLON   ! horizontal dimension
INTEGER,                   INTENT(IN) :: KLEV   ! vertical dimension
INTEGER, DIMENSION(KLON),  INTENT(IN) :: KLCL   ! index lifting condens. level
INTEGER, DIMENSION(KLON),  INTENT(IN) :: KCTL   ! index for cloud top level
INTEGER, DIMENSION(KLON),  INTENT(IN) :: KDPL   ! index for departure level
INTEGER, DIMENSION(KLON),  INTENT(IN) :: KPBL   ! index for top of source layer
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
REAL(dp), DIMENSION(KLON,KLEV), INTENT(IN)  :: PUTHL  ! updraft enthalpy (J/kg)
REAL(dp), DIMENSION(KLON,KLEV), INTENT(IN)  :: PURW   ! updraft total water (kg/kg)
REAL(dp), DIMENSION(KLON,KLEV), INTENT(IN)  :: PURC   ! updraft cloud water (kg/kg)
REAL(dp), DIMENSION(KLON,KLEV), INTENT(IN)  :: PURI   ! updraft cloud ice   (kg/kg)

REAL(dp), DIMENSION(KLON,KLEV), INTENT(OUT)  :: PTHC  ! conv. adj. grid scale theta
REAL(dp), DIMENSION(KLON,KLEV), INTENT(OUT)  :: PRWC  ! conv. adj. grid scale r_w
REAL(dp), DIMENSION(KLON,KLEV), INTENT(OUT)  :: PRCC  ! conv. adj. grid scale r_c
REAL(dp), DIMENSION(KLON,KLEV), INTENT(OUT)  :: PRIC  ! conv. adj. grid scale r_i
REAL(dp), DIMENSION(KLON,KLEV), INTENT(OUT)  :: PWSUB ! envir. compensating subsidence(Pa/s)

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
!REAL(dp), DIMENSION(KLON)     :: ZLM          ! latent heat of melting
 REAL(dp), DIMENSION(KLON)     :: ZCPH         ! specific heat C_ph
INTEGER, DIMENSION(KLON)  :: ITSTEP       ! fractional convective time step
INTEGER, DIMENSION(KLON)  :: ICOUNT       ! timestep counter
INTEGER, DIMENSION(KLON)  :: ILCL         ! index lifting condens. level
INTEGER, DIMENSION(KLON)  :: IWORK1       ! work array
REAL(dp),  DIMENSION(KLON)      :: ZWORK1, ZWORK2, ZWORK3 ! work arrays
REAL(dp),  DIMENSION(KLON)      :: ZWORK4, ZWORK5         ! work arrays
LOGICAL, DIMENSION(KLON)      :: GWORK1, GWORK3         ! work arrays
LOGICAL, DIMENSION(KLON,KLEV) :: GWORK4                 ! work array


!-------------------------------------------------------------------------------

!*       0.2    Initialize  local variables
!               ----------------------------


ZTIMC(:,:)  = 0.0_dp
ZTHES2(:)   = 0.0_dp
ZWORK1(:)   = 0.0_dp
ZWORK2(:)   = 0.0_dp
ZWORK3(:)   = 0.0_dp
ZWORK4(:)   = 0.0_dp
ZWORK5(:)   = 0.0_dp
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
ZOMG(:,:)  = 0.0_dp
PWSUB(:,:) = 0.0_dp


!*       3.     Compute limits on the closure adjustment factor so that the
!               inflow in convective drafts from a given layer can't be larger
!               than the mass contained in this layer initially.
!               ---------------------------------------------------------------

ZADJMAX(:) = 1000._dp
IWORK1(:) = ILCL(:)
JKP = MINVAL( KDPL(:) )
DO JK = JKP, IKE
  DO JI = 1, IIE
    IF( JK > KDPL(JI) .AND. JK <= IWORK1(JI) ) THEN
        ZWORK1(JI)  = PLMASS(JI,JK) / ( ( PUER(JI,JK) + 1.E-5_dp ) * PTIMEC(JI) )
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



!DO JITER = 1, 4  ! Enter adjustment loop to assure that all CAPE is
 DO JITER = 1, 2  ! Enter adjustment loop to assure that all CAPE is
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

             ZWORK1(:)   = - ( PUER(:,JKP) - PUDR(:,JKP) ) / PLMASS(:,JKP)

             PWSUB(:,JK) = PWSUB(:,JKP) - PDPRES(:,JK-1) * ZWORK1(:)
              ! we use PDPRES(JK-1) and not JKP in order to have zero subsidence
              ! at the first layer


!*       5.     Compute fractional time step. For stability or
!               mass conservation reasons one must split full time step PTIMEC)
!               ---------------------------------------------------------------

             ZWORK1(:) = XSTABT_S * PDPRES(:,JKP) / ( ABS( PWSUB(:,JK) ) + 1.E-10_dp )
              ! the factor XSTABT_S is used for stability reasons
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
       ZWORK1(JI) = PUDR(JI,JK) * PDPRES(JI,JK) / ( PLMASS(JI,JK) + .1_dp ) &
                  &                                         - PWSUB(JI,JK)
    ENDDO
    WHERE( GWORK1(:) .AND. ABS( ZWORK1(:) ) - .01_dp > 0.0_dp )
        GWORK1(:) = .FALSE.
        PTIMEC(:) = 1.E-1_dp
        ZWORK5(:) = 0.0_dp
    END WHERE
    DO JK = IKB, IKE
        PWSUB(:,JK) = PWSUB(:,JK) * ZWORK5(:)
    ENDDO
    GWORK4(:,1:IKB) = .FALSE.
    GWORK4(:,IKE:IKS) = .FALSE.

    ITSTEP(:) = INT( PTIMEC(:) / ZTIMEC(:) ) + 1
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


           ZTHLC(:,:) = ZTHLC(:,:) + ZTIMC(:,:) * (                    &
                      &    ZTHMFIN(:,:) + PUDR(:,:) * PUTHL(:,:)       &
                      & - ZTHMFOUT(:,:) - PUER(:,:) * PTHL(:,:)   )
           PRWC(:,:)  = PRWC(:,:) + ZTIMC(:,:)  *  (                   &
                      &    ZRWMFIN(:,:) + PUDR(:,:) * PURW(:,:)        &
                      & - ZRWMFOUT(:,:) - PUER(:,:) * PRW(:,:)    )
           PRCC(:,:)  = PRCC(:,:) + ZTIMC(:,:)  *  (                   &
                      &    ZRCMFIN(:,:) + PUDR(:,:) * PURC(:,:)        &
                      & - ZRCMFOUT(:,:) - PUER(:,:) * PRC(:,:)    )
           PRIC(:,:)  = PRIC(:,:) + ZTIMC(:,:)  *  (                   &
                      &    ZRIMFIN(:,:) + PUDR(:,:) * PURI(:,:)        &
                      & - ZRIMFOUT(:,:) - PUER(:,:) * PRI(:,:)    )


!******************************************************************************

         END WHERE

    ENDDO ! Exit the fractional time step loop


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
              ZTHES2(JI)  = ZWORK2(JI) * ZPI(JI) ** ( 1.0_dp - 0.28_dp * ZWORK3(JI) )   &
                          &        * EXP( ( 3374.6525_dp / ZWORK2(JI) - 2.5403_dp )  &
                          &        * ZWORK3(JI) * ( 1.0_dp + 0.81_dp * ZWORK3(JI) ) )

              ZWORK3(JI)  = PZ(JI,JK) - PZ(JI,JKP) * ZWORK4(JI) -                        &
                          & ( 1.0_dp - ZWORK4(JI) ) * ZZLCL(JI)    ! level thickness
              ZWORK1(JI)  = MAX( 0.01_dp, PUER(JI,JK)/ MAX( .1_dp, PUMF(JI,JK) ) )
              ZTHEU2(JI)  = ( 1.0_dp - ZWORK1(JI) ) * ZTHEU1(JI) + ZWORK1(JI) * ZTHES1(JI)
              ZWORK1(JI)  = ( ZTHEU1(JI) + ZTHEU2(JI) ) / ( ZTHES1(JI) + ZTHES2(JI) ) - 1.0_dp
              ZCAPE(JI)   = ZCAPE(JI) + RG * ZWORK3(JI) * MAX( 0.0_dp, ZWORK1(JI) )
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
       WHERE ( ZCAPE(:) /= 0.0_dp .AND. GWORK1(:) )                              &
             & ZADJ(:) = ZADJ(:) * XSTABC_S * PCAPE(:) / ( ZWORK1(:) + 1.E-8_dp )
       ZADJ(:) = MIN( ZADJ(:), ZADJMAX(:) )


!*         13.     Adjust mass flux by the factor ZADJ to converge to
!                  specified degree of stabilization
!                 ----------------------------------------------------

       CALL CONVECT_CLOSURE_ADJUST_SHAL( KLON, KLEV, ZADJ,                     &
                                       & PUMF, ZUMF, PUER, ZUER, PUDR, ZUDR    )


      IF ( COUNT( GWORK1(:) ) == 0 ) EXIT ! exit big adjustment iteration loop
                                          ! when all columns have reached
                                          ! desired degree of stabilization.

ENDDO  ! end of big adjustment iteration loop


        ! skip adj. total water array  to water vapor
DO JK = IKB, IKE
   PRWC(:,JK) = MAX( 0.0_dp, PRWC(:,JK) - PRCC(:,JK) - PRIC(:,JK) )
ENDDO


END SUBROUTINE CONVECT_CLOSURE_SHAL

!==============================================================================================

!    #########################################################################
     SUBROUTINE CONVECT_CLOSURE_ADJUST_SHAL( KLON, KLEV, PADJ,               &
                                   &  PUMF, PZUMF, PUER, PZUER, PUDR, PZUDR  )
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
!!   Last modified  15/11/96
!-------------------------------------------------------------------------------

!#include "tsmbkind.h"

!*       0.    DECLARATIONS
!              ------------

!USE YOE_CONVPAREXT

IMPLICIT NONE

!*       0.1   Declarations of dummy arguments :


INTEGER,                    INTENT(IN) :: KLON     ! horizontal dimension
INTEGER,                    INTENT(IN) :: KLEV     ! vertical dimension
REAL(dp), DIMENSION(KLON),  INTENT(IN) :: PADJ     ! mass adjustment factor


REAL(dp), DIMENSION(KLON,KLEV), INTENT(INOUT) :: PUMF  ! updraft mass flux (kg/s)
REAL(dp), DIMENSION(KLON,KLEV), INTENT(INOUT) :: PZUMF ! initial value of  "
REAL(dp), DIMENSION(KLON,KLEV), INTENT(INOUT) :: PUER  ! updraft entrainment (kg/s)
REAL(dp), DIMENSION(KLON,KLEV), INTENT(INOUT) :: PZUER ! initial value of  "
REAL(dp), DIMENSION(KLON,KLEV), INTENT(INOUT) :: PUDR  ! updraft detrainment (kg/s)
REAL(dp), DIMENSION(KLON,KLEV), INTENT(INOUT) :: PZUDR ! initial value of  "


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

     DO JK = IKB + 1, IKE
          PUMF(:,JK)  = PZUMF(:,JK)   * PADJ(:)
          PUER(:,JK)  = PZUER(:,JK)   * PADJ(:)
          PUDR(:,JK)  = PZUDR(:,JK)   * PADJ(:)
     ENDDO

END SUBROUTINE CONVECT_CLOSURE_ADJUST_SHAL

!==============================================================================================
!==============================================================================================

END MODULE MESSY_CONVECT_BECHTOLD_SHAL

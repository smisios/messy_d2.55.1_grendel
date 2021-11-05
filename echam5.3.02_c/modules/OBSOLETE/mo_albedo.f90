MODULE mo_albedo

  !- Description:
  !
  !  Follows ECHAM4, see also comments in ECHAM4 radint.
  !  U. Schlese, DKRZ, July 1993, original source
  !
  !  This module computes the surface albedo depending on the
  !  surface properties. The albedois computed for a longitude circle.
  !
  !  This module contains :
  !
  !  - constants used in the computation, see below
  !  - the subroutine albedo to do the computation
  !
  !  This module replaces parts of ECHAM4 radint
  !
  !- Author:
  !
  !  Marco Giorgetta, MPI, May 2000

  !=====================================================================

  USE mo_kind,      ONLY: dp
  USE mo_constants, ONLY: tmelt

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: albedos
  PUBLIC :: su_albedo

  REAL(dp), PUBLIC :: calbmns              ! minimum (glacier, snow on ice)
  REAL(dp), PUBLIC :: calbmxs              ! maximum (glacier, snow on ice)
  REAL(dp), PUBLIC :: calbmni              ! minimum (bare sea ice)


  ! constants for albedo computation
  REAL(dp), PARAMETER :: calbsea = 0.07_dp ! sea albedo
  REAL(dp), PARAMETER :: calbmn0 = 0.3_dp  ! for snow albedo minimum on land
  REAL(dp), PARAMETER :: calbmx0 = 0.8_dp  ! for snow albedo maximum on land
  REAL(dp), PARAMETER :: calbcan = 0.2_dp  ! albedo of snow covered canopy
  REAL(dp), PARAMETER :: cskyfac = 1.0_dp  ! constant in sky view factor
  !
  REAL(dp), PARAMETER :: calbmxi = 0.75_dp ! maximum (bare sea ice)
  !=====================================================================

CONTAINS

  SUBROUTINE su_albedo

    USE mo_control,       ONLY: nn, lcouple, lipcc
    USE mo_exception,     ONLY: finish

  ! Define parameters depending on coupled/uncoupled runs
  ! Define parameters depending on resolution (in coupled case)

    IF(lcouple .OR. lipcc) THEN
     IF (nn == 31) THEN
      calbmns = 0.65_dp
      calbmxs = 0.8_dp
      calbmni = 0.55_dp
!mz_ap_20100603+
#ifdef MESSY
     ELSEIF (nn == 42) THEN
      calbmns = 0.65_dp
      calbmxs = 0.8_dp
      calbmni = 0.55_dp
#endif
!mz_ap_20100603-
     ELSE IF (nn == 63 .OR. nn == 106) THEN
      calbmns = 0.75_dp
      calbmxs = 0.85_dp
      calbmni = 0.65_dp
     ELSE
      CALL finish ('su_albedo', 'Truncation not supported in coupled runs.')
     ENDIF
    ELSE IF (nn == 319) THEN
      calbmns = 0.75_dp
      calbmxs = 0.85_dp
      calbmni = 0.65_dp
    ELSE
      calbmns = 0.6_dp
      calbmxs = 0.8_dp
      calbmni = 0.5_dp
    ENDIF

  END SUBROUTINE su_albedo

  SUBROUTINE albedos(klon,                 &
       &             loland,loglac,        &
       &             pforest,pseaice,      &
       &             pcvs,psni,            &
       &             pcvsc,pvlt,           &
       &             pfrl,pfrw,pfri,       &
       &             ptslm1,ptsi,palb,     &
       &             palsol,palsow,palsoi, &
       &             palbedo)

    IMPLICIT NONE

    ! INPUT
    ! -----

    INTEGER ,INTENT(in)                  :: klon  ! number of longitudes
    LOGICAL, INTENT(in), DIMENSION(klon) :: loland ! land mask
    LOGICAL, INTENT(in), DIMENSION(klon) :: loglac ! glacier mask
    REAL(dp),    INTENT(in), DIMENSION(klon) :: pforest! forest cover
    REAL(dp),    INTENT(in), DIMENSION(klon) :: pseaice! sea ice cover
    REAL(dp),    INTENT(in), DIMENSION(klon) :: pcvs   ! snow cover (ground)
    REAL(dp),    INTENT(in), DIMENSION(klon) :: pcvsc  ! snow cover (canopy)
    REAL(dp),    INTENT(in), DIMENSION(klon) :: psni   ! snow depth on ice
    REAL(dp),    INTENT(in), DIMENSION(klon) :: pfrl   ! fraction of land
    REAL(dp),    INTENT(in), DIMENSION(klon) :: pfrw   ! fraction of water
    REAL(dp),    INTENT(in), DIMENSION(klon) :: pfri   ! fraction of ice
    REAL(dp),    INTENT(in), DIMENSION(klon) :: ptslm1 ! land surface temp.
    REAL(dp),    INTENT(in), DIMENSION(klon) :: ptsi   ! ice surface temp.
    REAL(dp),    INTENT(in), DIMENSION(klon) :: palb   ! background albedo
    REAL(dp),    INTENT(in), DIMENSION(klon) :: pvlt   ! leaf area index

    ! OUTPUT
    ! ------

    REAL(dp),    INTENT(out),DIMENSION(klon) :: palsol  ! land albedo
    REAL(dp),    INTENT(out),DIMENSION(klon) :: palsow  ! water albedo
    REAL(dp),    INTENT(out),DIMENSION(klon) :: palsoi  ! ice albedo
    REAL(dp),    INTENT(out),DIMENSION(klon) :: palbedo ! grid-mean albedo

    REAL(dp) :: ztalb    ! upper temp. limit for cold snow albedo
    REAL(dp) :: ztsalb   ! upper temp. limit for cold snow albedo on sea ice
    REAL(dp) :: zalbmax  ! maximum snow albedo
    REAL(dp) :: zalbmin  ! minimum snow albedo
    REAL(dp) :: zdalb    ! snow albedo change per deg C
    REAL(dp) :: zalbsn   ! temperature dependent snow albedo
    REAL(dp) :: zsvf     ! 
    REAL(dp) :: zalgrd   ! 
    REAL(dp) :: zalcan   ! 
    REAL(dp) :: zalfor   ! 

    INTEGER :: jl    ! loop index

    ztalb=tmelt-5.0_dp
    ztsalb=tmelt-1.0_dp

    DO jl = 1, klon

       palsol(jl)=palb(jl)               ! set to background albedo
       palsow(jl)=palb(jl)               !           "
       palsoi(jl)=palb(jl)               !           "

       ! land

       IF (loland(jl)) THEN
          ! minimum and maximum snow albedo
          IF (loglac(jl)) THEN
             ! on glacier covered land
             zalbmin=calbmns
             zalbmax=calbmxs
          ELSE
             ! on glacier-free land
             zalbmin=calbmn0
             zalbmax=calbmx0
          END IF
          ! temperature dependent snow albedo
          IF (ptslm1(jl)>=tmelt) THEN
             zalbsn=zalbmin
          ELSE IF (ptslm1(jl)<ztalb) THEN
             zalbsn=zalbmax
          ELSE
             zdalb=(zalbmax-zalbmin)/(tmelt-ztalb)
             zalbsn=zalbmin+zdalb*(tmelt-ptslm1(jl))
          END IF
          ! final land albedo
          IF (loglac(jl)) THEN
             palsol(jl)=zalbsn
          ELSE
             zsvf=EXP(-cskyfac*pvlt(jl))
             zalgrd=pcvs(jl)*zalbsn+(1._dp-pcvs(jl))*palb(jl)
             zalcan=pcvsc(jl)*calbcan+(1._dp-pcvsc(jl))*palb(jl)
             zalfor=zsvf*zalgrd+(1._dp-zsvf)*zalcan
             palsol(jl)=(1._dp-pforest(jl))*zalgrd+pforest(jl)*zalfor
             palsol(jl)=MAX(palsol(jl),palb(jl))
          END IF
       END IF

       ! ice

       IF (pseaice(jl)>0._dp) THEN
          ! minimum and maximum albedo
          IF (psni(jl)>0.01_dp) THEN
             ! on snow covered sea ice
             zalbmin=calbmns
             zalbmax=calbmxs
          ELSE
             ! on bare sea ice
             zalbmin=calbmni
             zalbmax=calbmxi
          END IF
          ! temperature dependent snow albedo
          IF (ptsi(jl)>=tmelt) THEN
             palsoi(jl)=zalbmin
          ELSE IF (ptsi(jl)<ztsalb) THEN
             palsoi(jl)=zalbmax
          ELSE
             zdalb=(zalbmax-zalbmin)/(tmelt-ztsalb)
             palsoi(jl)=zalbmin+zdalb*(tmelt-ptsi(jl))
          END IF
       END IF

       ! water

       palsow(jl)=calbsea

       ! average

       palbedo(jl)= pfrl(jl)*palsol(jl)+ &
            &       pfrw(jl)*palsow(jl)+ &
            &       pfri(jl)*palsoi(jl)

    END DO

  END SUBROUTINE albedos

!=======================================================================

END MODULE mo_albedo

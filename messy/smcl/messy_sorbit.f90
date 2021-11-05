! **********************************************************************
MODULE messy_sorbit
! **********************************************************************

  ! Sample along satellite orbit (constant local time)
  !
  ! CORE MODULE (MESSy/SMCL)
  !
  ! Author: Patrick Joeckel, MPICH, January 2007
  !

  ! ----------- >

  USE messy_main_constants_mem, ONLY: DP

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: DP

  ! ----------- <

  CHARACTER(LEN=*),PARAMETER, PUBLIC :: modstr = 'sorbit'
  CHARACTER(LEN=*),PARAMETER, PUBLIC :: modver = '1.2'
  ! MAXIMUM LENGTH OF STRING FOR CHANNEL(S)/OBJECT(S)
  INTEGER,         PARAMETER, PUBLIC :: STRLEN = 1000

  PUBLIC :: polst ! POLAR ORBITER LOCAL SOLAR TIME
  PUBLIC :: ltadd ! add seconds to local time, second of day (modulo)

CONTAINS

  ! ----------------------------------------------------------------------
  ELEMENTAL SUBROUTINE polst(lst, incl, flag, lst_eq, lat)

    ! POLAR ORBITER LOCAL SOLAR TIME

    ! References:
    !
    ! STEPHEN S. LEROY, The Effects of Orbital Precession on Remote Climate
    ! Monitoring, Journal of Climate, 14, 4330-4337, 2001.
    !
    ! Bronstein, I. N., K. A. Semendjajew, G. Musiol, H. M\"{u}hlig,
    ! Taschenbuch der Mathematik, 6., vollst. überarb. und erg. Auflage, 
    ! Verlag Harri Deutsch, Frankfurt am Main, ISBN-10 3-8171-2006-0, 2005.
    ! (Note: Equation 3.203g)

    IMPLICIT NONE
    INTRINSIC :: ASIN, TAN

    ! I/O
    REAl(DP), INTENT(OUT) :: lst    ! local solar time [hour of day]
    REAL(DP), INTENT(IN)  :: incl   ! inclination [degrees]
    REAL(DP), INTENT(IN)  :: flag   ! ascending (1), descending (-1) orbit
    REAL(DP), INTENT(IN)  :: lst_eq ! local solar time at equator [hour of day]
    REAL(DP), INTENT(IN)  :: lat    ! latitude [degrees]

    ! LOCAL
    REAL(DP) :: pi     ! pi
    REAL(DP) :: d2r    ! degrees -> radian
    REAL(DP) :: arg

    pi  = 4*ATAN(1.0_dp)
    d2r = pi/180.0_dp

    arg = TAN(lat*d2r)/TAN(incl*d2r)
    IF ((arg >= -1.0_dp) .AND. (arg <= 1.0_dp)) THEN
       lst = flag * ASIN(arg) * 24.0_dp/(2.0_dp*pi) + lst_eq
       IF (lst <   0.0_dp) lst = lst + 24.0_dp
       IF (lst >= 24.0_dp) lst = lst - 24.0_dp
    ELSE
       lst = -1.0_dp
    END IF

  END SUBROUTINE polst
  ! ----------------------------------------------------------------------

  ! ----------------------------------------------------------------------
  INTEGER FUNCTION LTADD(LT, DT)

    IMPLICIT NONE
    
    INTRINSIC :: MOD

    ! I/O
    INTEGER, INTENT(IN) :: LT ! local time [second of day]
    INTEGER, INTENT(IN) :: DT ! time interval [seconds]

    ! LOCAL
    INTEGER, PARAMETER :: SPD = 86400 ! seconds per day

    LTADD = LT + DT
    IF (LTADD < 0) LTADD = LTADD + SPD
    LTADD = MOD(LTADD, SPD)

  END FUNCTION LTADD
  ! ----------------------------------------------------------------------

! **********************************************************************
END MODULE messy_sorbit
! **********************************************************************

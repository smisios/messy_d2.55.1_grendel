! **********************************************************************
MODULE messy_ptracini
! **********************************************************************

  ! TRACER INITIALISATION 
  !
  ! CORE MODULE (MESSy/SMCL)
  !
  ! Author: Patrick Joeckel, MPICH, December 2003
  !

  ! ----------- >

  USE messy_main_constants_mem, ONLY: DP

  IMPLICIT NONE
  PRIVATE

  INTRINSIC :: TINY

  PUBLIC :: DP

  ! ----------- <

  CHARACTER(LEN=*),PARAMETER, PUBLIC :: modstr = 'ptracini'
  CHARACTER(LEN=*),PARAMETER, PUBLIC :: modver = '3.1'
  REAL(DP),        PARAMETER, PUBLIC :: ZERO_EPS = TINY(0.0_DP)    ! zero

  PUBLIC :: solve
  PUBLIC :: vmr2conc

CONTAINS

  ! ----------------------------------------------------------------------
  ELEMENTAL SUBROUTINE solve( xte, x, q, k0, dt )

    ! PTRACINI MODULE ROUTINE (CORE)
    ! MAIN INTERFACE ROUTINE OF PTRACINI
    !
    ! PERFORMS FORCING INTEGRATION FOR ONE POINT/COLUMN/VECTOR
    !
    ! OUTPUT:
    !        xte        tracer tendency for t0-:-t0+dt  (?/s)
    !
    ! INPUT:
    !        x          tracer at t=t0  [mol/mol]
    !        q          source strength [mol/mol/s]
    !        k0         0-order reaction coefficient [1/s]
    !                   ( = decay constant = loss rate)
    !        dt         time step length [s]
    !
    ! Note: [q] and [x] must be in compatible units,
    !       e.g. mol/mol/s and mol/mol
    !
    ! Author: Patrick Joeckel, MPICH, Dec 2004

    IMPLICIT NONE

    INTRINSIC :: EXP

    ! I/O
    REAL(DP), INTENT(INOUT) :: xte   ! tracer tendency [?/s]
    REAL(DP), INTENT(IN)    :: x     ! tracer at t=t0
    REAL(DP), INTENT(IN)    :: q     ! tracer emission  rate [?/s]
    REAL(DP), INTENT(IN)    :: k0    ! 0-order reaction coeff. [1/s]
    REAL(DP), INTENT(IN)    :: dt    ! time step in seconds

    ! LOCAL
    REAL(DP)    :: ep

    ! TENDENCIES:
    ! DIFFERENTIAL EQUATION (HERE ANALYTICALLY SOLVED)
    ! dc/dt = -(1/tau)*c + q ; c(0) = c0
    !    with  1/tau := k0
    ! CASE 1: l <> 0 => tau < infinity:
    !         --> c(t) = exp(-t/tau)*(c0-q*tau+q*tau*exp(t/tau))
    !             c(t) = c0*(exp(-t/tau)-1) + q*tau*(1-exp(-t/tau))
    !             --------------------------------------
    !         ==> c(t) - c0 = (c0-q*tau)*(exp(-t/tau)-1)
    !             --------------------------------------
    ! CASE 2: l = 0  => tau -> infinity
    !         --> dc/dt = q ; c(0) = c0
    !         --> c(t) = c0 + q*t
    !         -------------------
    !         ==> c(t) - c0 = q*t
    !         -------------------
    IF (k0 <= ZERO_EPS) THEN      ! CASE 2
       xte  = q
    ELSE                     ! CASE 1
       ep   = (exp(-dt*k0)-1._DP)/dt
       xte  = (x-q/k0) * ep
    END IF

  END SUBROUTINE solve
  ! ----------------------------------------------------------------------

  ! -------------------------------------------------------------------------
  ELEMENTAL SUBROUTINE vmr2conc (f, p, t)

    ! PTRACINI MODULE ROUTINE (CORE)
    !
    ! PERFORMS CONVERSION OF BOX/COLUMN/VECTOR FROM
    ! VOLUME MIXING RATIO (vmr = mol/mol)
    ! TO CONCENTRATION (conc) in 1/cm^3       !!!
    !
    ! conc = vmr * p/((kb*Na)*t) * Na * 1.0E-06
    !              mol/m^3       | conversion to cm^(-3)
    ! => conc = vmr * (p/kb*t) * 1.OE-06      in cm^(-3)
    !
    ! NEEDS PRESSURE (p) and TEMPERATURE (t)
    !
    ! Note: For the chosen unit (cm^(-3)), Na
    !           (Avogadro-constant) cancels out
    !
    ! Author: Patrick Joeckel, MPICH, Dec 2004

    USE messy_main_constants_mem, ONLY: kb=>k_B

    IMPLICIT NONE

    ! I/O
    REAL(DP),           INTENT(INOUT) :: f  ! volume mixing ratio (mol/mol)
    REAL(DP),           INTENT(IN)    :: p  ! pressure (Pa)
    REAL(DP),           INTENT(IN)    :: t  ! temperature (K)

    ! vmr <- conc
    f = f * (p/(kb*t)) * 1.0E-06_DP

  END SUBROUTINE vmr2conc
  ! -------------------------------------------------------------------------

! **********************************************************************
END MODULE messy_ptracini
! **********************************************************************

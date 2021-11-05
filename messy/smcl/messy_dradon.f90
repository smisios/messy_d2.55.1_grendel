! **********************************************************************
! RADON module for BL diagnostic
! Version: see 'modver' below
!
! Author : Patrick Joeckel, MPICH, July  2003
!          Patrick Joeckel, DLR,   October 2018
!
! This module uses 222-Rn as a tracer for diagnosing the BL.
!
! Reference(s):
!
! Jöckel, P., Kerkweg, A., Pozzer, A., Sander, R., Tost, H., Riede, H.,
! Baumgaertner, A., Gromov, S., & Kern, B.: Development cycle 2 of the Modular
! Earth Submodel System (MESSy2), Geoscientific Model Development, 3, 717-752,
! doi: 10.5194/gmd-3-717-2010, URL http://www.geosci-model-dev.net/3/717/2010/
! (2010)
!
! **********************************************************************
!
! 222Rn -a(3.8d)-> 218Po -a(3min)-> 214Pb -b(27min)-> 214Bi
! 214Bi - (20min)-> 214Po -a(180 micro-s)-> 210Pb
! 210Pb -b(22y)-> 210Bi -b(5d)-> 210Po -a(138d)-> 206Pb
!
! 1 Bq = 1/s
! 
! dN/dt = -N/tau ; tau -> t2 ; N = N0 * 2^(-t/t2)
! => N0-N = N0(1-2^(-t/t2)) 
! =>
! 222Rn:  1 atom/(cm2 s)  => 21.112 mBq/(m2 s)
!         1 atom/(m2 s)   => 21.112E-04 mBq/(m2 s)
!      => 1 mBq/(m2 s)    => 473.66426676772 atoms/(m2 s)
! **********************************************************************

MODULE messy_dradon

  ! ----------- >

  USE messy_main_constants_mem, ONLY: DP

  IMPLICIT NONE
  PRIVATE
  SAVE
  
  INTRINSIC :: EXP, LOG, PRESENT

  PUBLIC :: DP

  ! ----------- <

  ! GLOBAL PARAMETERS
  CHARACTER(LEN=*), PARAMETER, PUBLIC :: modstr = 'dradon' ! name of module
  CHARACTER(LEN=*), PARAMETER, PUBLIC :: modver = '3.0'    ! module version

  ! METHOD FOR RADON FLUX CALCULATION
  ! 0 = const; 1 = offline ; 2 = online
  INTEGER, PARAMETER, PUBLIC :: I_Rn_flux_const   = 0
  INTEGER, PARAMETER, PUBLIC :: I_Rn_flux_offline = 1
  INTEGER, PARAMETER, PUBLIC :: I_Rn_flux_online  = 2

  ! GLOBAL SWITCH
  LOGICAL, PUBLIC :: ldradon = .false.   ! default: no Radon diagnostic

  ! NUMBER OF LEMENTS IN DECAY CHAIN
  INTEGER, PARAMETER, PUBLIC :: NI = 5
  ! HALF-LIFE
  ! e.g., 222Rn-half-life 3.8d
  REAL(DP), DIMENSION(NI), PARAMETER :: t2 = (/ &
        3.8_DP*24.0_DP*3600.0_DP &            ! 222Rn -> 218Po
      , 3._DP*60._DP &                        ! 218Po -> 214Pb
      , 27._DP*60._DP &                       ! 214Pb -> 214Bi
      , 20._DP*60._DP &                       ! 214Bi -> 214Po -> 210Pb
      , 22._DP*365._DP*24.0_DP*3600.0_DP /)   ! 210Pb -> ... -> 206Pb

  ! NAMES
  CHARACTER(LEN=5), DIMENSION(NI), PARAMETER, PUBLIC :: iname = (/ &
       '222Rn', '218Po', '214Pb', '214Bi', '210Pb' /)

  ! PRE-CALCULTED FOR CHAIN INTEGRATION
  ! (LF: workaround for Lahey/Fujitsu Compiler 8.10b)
#if !defined(LF) && !defined(__PGI)
  ! DECAY COEFF. = ln(2) / t_1/2
  REAL(DP), DIMENSION(NI),PARAMETER, PUBLIC :: lam = log(2.0_DP) / t2
#else
  REAL(DP), DIMENSION(NI), PUBLIC :: lam  ! DECAY COEFF. = ln(2) / t_1/2
#endif

  REAL(DP), DIMENSION(NI,NI)    :: PR1   ! PRODUCT 1
  REAL(DP), DIMENSION(NI,NI,NI) :: PR2   ! PRODUCT 2

  ! ----------------------------------------------------------------------
  ! &CTRL NAMELIST PARAMETERS
  INTEGER,  PUBLIC :: I_Rn_flux_method  = I_Rn_flux_const
  ! values over land/sea for const. Rn. flux; unit: atoms m^(-2) s^(-1)
  REAL(DP), PUBLIC :: R_Rn_cflux_land  = 10000.0_DP ! default
  REAL(DP), PUBLIC :: R_Rn_cflux_ocean = 0.0_DP     ! default  
  ! ----------------------------------------------------------------------

  ! DRADON SMCL SUBROUTINES (MAIN ENTRY POINTS)
  PUBLIC :: dradon_read_nml_ctrl   ! read CTRL namelist and initialize
  ! only decay:
  PUBLIC :: int_radon              ! integrate one time step (1D)
  ! integrate chain:
  PUBLIC :: init_radon_chain       ! pre-calculation of constants
  PUBLIC :: int_radon_chain        ! integrate one time step of chain

  ! DRADON HELPER ROUTINES
  PUBLIC :: flux2vmrs              ! converson from atoms / (m^2 s) to mol/mol/s

CONTAINS

! ************************************************************************
! RADON CORE ROUTINES
! ************************************************************************
! ----------------------------------------------------------------------
  SUBROUTINE dradon_read_nml_ctrl(status, iou)

    ! Radon MODULE ROUTINE (CORE)
    !
    ! READ Radon NAMELIST, CHECK IT, AND INITIALIZE GLOBAL VARIABLES
    !
    ! Author: Patrick Joeckel, MPICH, Jul 2003

    USE messy_main_tools, ONLY: read_nml_open, read_nml_check, read_nml_close

    IMPLICIT NONE

    ! I/O
    INTEGER, INTENT(OUT) :: status ! error status
    INTEGER, INTENT(IN)  :: iou    ! logical I/O unit

    NAMELIST /CTRL/   I_Rn_flux_method &
                    , R_Rn_cflux_land  &
                    , R_Rn_cflux_ocean

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'dradon_read_nml_ctrl'
    LOGICAL                     :: lex          ! file exists ?
    INTEGER                     :: fstat        ! file status
#if defined(LF) || defined(__PGI)
    INTEGER                     :: i
#endif

    ! INITIALIZE
    status = 1 ! ERROR

    ! INITIALIZE GLOBAL CONTROL VARIABLES
    ! -> DEFAULT VALUES ARE SET AT DECLARATION ABOVE

    CALL read_nml_open(lex, substr, iou, 'CTRL', modstr)
    IF (.not.lex) RETURN    ! <modstr>.nml does not exist

    READ(iou, NML=CTRL, IOSTAT=fstat)
    CALL read_nml_check(fstat, substr, iou, 'CTRL', modstr)
    IF (fstat /= 0) RETURN  ! error while reading namelist

    ! DIAGNOSE NAMELIST AND SET GLOBAL SWITCHES
    !
    ! CHECK CONSTAN/ONLINE/OFFLINE FLUX
    SELECT CASE(I_Rn_flux_method)
    CASE(I_Rn_flux_const)
       WRITE(*,*) '... flux calculation method : CONSTANT'
       WRITE(*,*) '                     land   : ',R_Rn_cflux_land  &
            , 'atoms m^(-2) s^(-1)'
       WRITE(*,*) '                     ocean  : ',R_Rn_cflux_ocean  &
            , 'atoms m^(-2) s^(-1)'
    CASE(I_Rn_flux_offline)
       WRITE(*,*) '... flux calculation method : OFFLINE'
    CASE(I_Rn_flux_online)
       WRITE(*,*) '... flux calculation method : ONLINE'
       WRITE(*,*) '*** ERROR: ONLINE Rn-flux MODEL NOT YET IMPLEMENTED !'
       RETURN  ! status = 1
    CASE DEFAULT
       WRITE(*,*) '*** ERROR: UNKNOWN Rn-flux CALCULATION METHOD !'
       RETURN  ! status = 1
    END SELECT
    !
    ldradon = .true.

    CALL read_nml_close(substr, iou, modstr)

#if defined(LF) || defined(__PGI)
    ! earliest possible place ...
    DO i=1, NI
       lam(i) = log(2.0_DP) /  t2(i)
    END DO
#endif

    status = 0  ! no ERROR

  END SUBROUTINE dradon_read_nml_ctrl
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
  SUBROUTINE init_radon_chain

    IMPLICIT NONE

    ! LOCAL
    INTEGER  :: i, q, m, j, k

    PR1(:,:) = 1.0_DP
    DO i=1, NI

       DO m=1, i-1
          DO q=m, i-1
             PR1(i,m) = PR1(i,m) * lam(q)
          END DO
       END DO

    END DO

    PR2(:,:,:) = 1.0_DP
    DO i=1, NI

       DO m=1, i-1

          DO k=m, i

             prod_loop: DO j=m,i
                IF (j==k) CYCLE
                PR2(i,m,k) = PR2(i,m,k) * (lam(j)-lam(k))
             END DO prod_loop

          END DO

       END DO

    END DO

  END SUBROUTINE init_radon_chain
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
  SUBROUTINE int_radon_chain(y, x, dt)

    !
    ! INTEGRATION OF THE SET OF DGLs:
    !
    ! dN1/dt = -k1*N1
    ! dN2/dt =  k1*N1 - k2*N2
    ! dN3/dt =          k2*N2 - k3*N3
    ! ...
    ! 
    !
    ! n: species (number in decay chain)
    ! N: number of atoms
    ! t: time interval
    ! k: decay coeff.
    !
    ! GENERAL SOLUTION: BAREMAN SOLUTION
    !
    ! Reference:
    !
    ! Dobromir S. Pressyanov, Short solution of the radioactive decay
    !   chain equations, Am. J. Phys., 70(4), 2002.
    !

    IMPLICIT NONE

    ! I/O
    REAL(DP), DIMENSION(NI), INTENT(OUT) :: y   ! tendency
    REAL(DP), DIMENSION(NI), INTENT(IN)  :: x   ! value at t=0
    REAL(DP)               , INTENT(IN)  :: dt  ! time step in seconds

    ! LOCAL    
    INTEGER                 :: i, m, k
    REAL(DP), DIMENSION(NI) :: xp       ! value at t=dt
    REAL(DP)                :: tmp

    ! INIT
    xp(:) = 0.0

    DO i=1, NI

       xp(i) = x(i)*exp(-lam(i)*dt) 

       sum_loop: DO m=1,i-1
          tmp = 0.0_DP
          DO k=m, i
             tmp = tmp + exp(-lam(k)*dt)/PR2(i,m,k)
          END DO
          xp(i) = xp(i) + x(m) * PR1(i,m)*tmp
       END DO sum_loop

    END DO

    ! TENDENCY
    y(:) = (xp(:)-x(:))/dt

  END SUBROUTINE int_radon_chain
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
  ELEMENTAL SUBROUTINE int_radon(y, x, dt, q)

    ! Radon MODULE ROUTINE (CORE)
    ! MAIN INTERFACE ROUTINE OF Radon DIAGNOSTIC
    !
    ! INPUT:
    !        dt         time step length
    !        q          Radon source strength (mol/mol/s)
    !        x          tracer, e.g., 222-Rn at t=t0
    !
    ! OUTPUT:
    !        y          Radon tendency for t0-:-t0+dt  (?/s)
    !
    ! Note:
    !       [q] and [x] must be in compatible units,
    !           e.g. mol/mol/s and mol/mol
    !
    ! Author: Patrick Joeckel, MPICH, Jul 2003

    IMPLICIT NONE

    ! I/O
    REAL(DP), INTENT(OUT)             :: y   ! Radon tendency
    REAL(DP), INTENT(IN)              :: x   ! Radon
    REAL(DP), INTENT(IN)              :: dt  ! time step in seconds
    REAL(DP), INTENT(IN), OPTIONAL    :: q   ! Radon source

    ! LOCAL
    REAL(DP) :: zq  ! source term
    REAL(DP) :: tau ! decay time

    IF (PRESENT(q)) THEN
       zq = q
    ELSE
       zq = 0.0
    END IF

    ! TENDENCY:
    ! DIFFERENTIAL EQUATION (HERE ANALYTICALLY SOLVED)
    ! dc/dt = -(1/tau)*c + q ; c(0) = c0
    !    with  tau = (t2_Rn / ln(2))
    !         --> c(t) = exp(-t/tau)*(c0-q*tau+q*tau*exp(t/tau))
    !             c(t) = c0*(exp(-t/tau)-1) + q*tau*(1-exp(-t/tau))
    !             --------------------------------------
    !         ==> c(t) - c0 = (c0-q*tau)*(exp(-t/tau)-1)
    !             --------------------------------------

    tau = t2(1) / log(2.0_DP)

    y = ( x - zq * tau ) * ( exp(-dt/tau ) - 1._DP ) / dt

  END SUBROUTINE int_radon
! ----------------------------------------------------------------------

! -------------------------------------------------------------------------
  ELEMENTAL SUBROUTINE flux2vmrs(vmrs, flux, delp, mwair, na, g)

    ! Radon MODULE ROUTINE (HELPER)
    !
    ! CONVERTS Radon flux unit
    ! from atoms/(m^2 s) to mol/mol/s
    !
    ! Note: The required constants
    !            mole-weight of air (mwair) in g/mol
    !       and  Avogadro-Constant  (na)    in molec/mol
    !       and  Earth acceleration (g)     in m/s2
    !       can optionally be passed from outside.
    !
    ! Author: Patrick Joeckel, MPICH, July 2003

    IMPLICIT NONE

    ! I/O
    REAL(DP), INTENT(OUT)             :: vmrs  ! mol/mol/s
    REAL(DP), INTENT(IN)              :: flux  ! atoms / (m^2 s)
    REAL(DP), INTENT(IN)              :: delp  ! pressure difference
    REAL(DP), INTENT(IN), OPTIONAL    :: mwair ! molarmass of dry air (g/mol)
    REAL(DP), INTENT(IN), OPTIONAL    :: na    ! Avogadro constant (molec/mol)
    REAL(DP), INTENT(IN), OPTIONAL    :: g     ! Earth acceleration

    ! LOCAL
    REAL(DP) :: zmwair   ! molarmass of dry air [g/mol]
    REAL(DP) :: zna      ! Avogadro constant
    REAL(DP) :: zg       ! Earth acceleration   [m/s^2]

    IF (PRESENT(mwair)) THEN
       zmwair = mwair
    ELSE
       zmwair = 28.970_dp      ! g/mol
    END IF

    IF (PRESENT(na)) THEN
       zna = na
    ELSE
       zna = 6.022045E+23_dp   ! molec/mol
    END IF

    IF (PRESENT(g)) THEN
       zg = g
    ELSE
       zg = 9.81_dp   ! m/s^2
    END IF

    vmrs = (flux * zg * zmwair) / (1000.0_DP * zna * delp)

  END SUBROUTINE flux2vmrs
! -------------------------------------------------------------------------

END MODULE messy_dradon
! **********************************************************************

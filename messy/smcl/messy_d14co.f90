! **********************************************************************
MODULE messy_d14co
! **********************************************************************

! **********************************************************************
! 14CO module for OH/STE diagnostic
! Version: see "modver" below
!
! Author : Patrick Joeckel, MPICH, October  2002
!                                  November 2002
!                                  December 2002
!                                  July     2004
!
! This module uses 14CO as a tracer for diagnosing
! atmospheric chemistry (OH distribution and seasonality)
! and transport (stratosphere - troposphere exchange)
!
! References:
! (1) Patrick Joeckel, Carl A. M. Brenninkmeijer, Mark G. Lawrence,
!     Adriaan B. M. Jeuken, and Peter F. J. van Velthoven,
!     Evaluation of stratosphere - troposphere exchange and the hydroxyl
!     radical distribution in 3-dimensional global atmospheric models using
!     observations of cosmogenic 14CO,
!     J. Geophys. Res., 107(D20), 4446, doi:10.1029/2001JD001324, 2002.
! (2) Patrick Joeckel, and Carl A.M. Brenninkmeijer,
!     The seasonal cycle of cosmogenic 14CO at surface level: A solar
!     cycle adjusted, zonal-average climatology based on observations,
!     J. Geophys. Res., 107(D22), 4656, doi:10.1029/2001JD001104, 2002.
! (3) Patrick Joeckel,
!     Cosmogenic 14CO as tracer for atmospheric chemistry and transport,
!     Dissertation, Combined Faculties for the Natural Sciences and for
!     Mathematics  of the Rupertus Carola University of
!     Heidelberg, Germany, 2000.
!     (http://www.ub.uni-heidelberg.de/archiv/1426)
! (4) D. C. McCabe, T. Gierczak, R. K. Talukdar and A. R. Ravishankara,
!     Kinetics of the reaction OH + CO under atmospheric conditions,
!     Geophys. Res. Lett., 28, 3135-3138, 2001.
!
! **********************************************************************

  ! ----------- >

  USE messy_main_constants_mem, ONLY: DP

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: DP

  ! ----------- <

  ! 14CO CORE ROUTINES (MAIN INTERFACES)
  PUBLIC :: parse_str                ! parse a string
  PUBLIC :: tfrac                    ! tropospheric fraction of box
  PUBLIC :: int_14co                 ! integrate one time step (1D)

  ! 14CO PUBLIC HELPER ROUTINES
  PUBLIC :: adj_tend                 ! adjust tendencies (linearize)
  PUBLIC :: vmr2conc                 ! convert mol/mol to cm^(-3) (1D)
  PUBLIC :: intpol_p                 ! interpolate on pressure levels (1D)
  PUBLIC :: merge_p                  ! merge two fields (e.g. strat/trop) (1D)
  PUBLIC :: mgs2vmrs                 ! convert molec/(g s) to mol/mol/s (1D)

  ! 14CO PRIVATE HELPER ROUTINES
  ! PRIVATE :: intpol_init           ! for online vertical ...  (1D)
  ! PRIVATE :: intpol_step           ! ... interpolation (pressure axis) (1D)
  ! PRIVATE :: h2p                   ! calc. press from hyb. + surf. press (1D)

  INTRINSIC :: SIZE, ABS, EXP, INT, MAX, MIN &
             , REAL, SIGN, PRESENT, LOG10, TINY, ASSOCIATED

  ! GLOBAL PARAMETER
  CHARACTER(LEN=*), PARAMETER, PUBLIC :: modstr = 'd14co' ! name of module
  CHARACTER(LEN=*), PARAMETER, PUBLIC :: modver = '2.1'   ! module version
  REAL(DP),         PARAMETER :: EPSILON = TINY(0.0_DP)    ! zero

CONTAINS

! ************************************************************************
! 14CO CORE ROUTINES
! ************************************************************************

! ----------------------------------------------------------------------
  SUBROUTINE parse_str(status, strlen, str, var, unit, zscale)

    USE messy_main_tools, ONLY: strcrack

    IMPLICIT NONE

    INTRINSIC :: TRIM, ADJUSTL

    ! I/O
    INTEGER,          INTENT(OUT)   :: status   ! status information
    INTEGER,          INTENT(IN)    :: strlen   ! max. length of strings
    CHARACTER(LEN=*), INTENT(IN)    :: str      ! string to parse
    CHARACTER(LEN=*), INTENT(INOUT) :: var      ! variable name
    CHARACTER(LEN=*), INTENT(INOUT) :: unit     ! unit
    REAL(DP),         INTENT(INOUT) :: zscale   ! z-axis scaling

    ! LOCAL
    CHARACTER(LEN=strlen),               POINTER     :: sl1(:)
    CHARACTER(LEN=strlen),               POINTER     :: sl2(:)
    INTEGER :: n, m
    INTEGER :: i
    INTEGER :: iostat

    status = 1 ! ERROR

    NULLIFY(sl1)
    NULLIFY(sl2)

    CALL strcrack(str, ';', sl1, n)
    DO i=1, n

       CALL strcrack(sl1(i), '=', sl2, m)
       IF (SIZE(sl2) == 2) THEN ! mz_ak_20060731
          IF (TRIM(ADJUSTL(sl2(2))) == '')  THEN
             status = 2    ! EMPTY SPECIFICATION
             RETURN
          ENDIF
       END IF                   ! mz_ak_20060731

       SELECT CASE(TRIM(ADJUSTL(sl2(1))))
          CASE('VAR')
             var = TRIM(ADJUSTL(sl2(2)))
          CASE('UNIT')
             unit = TRIM(ADJUSTL(sl2(2)))
          CASE('ZSCALE')
             READ(sl2(2),*,IOSTAT=iostat) zscale
             IF (iostat /= 0) THEN
                status = 3  ! ERROR IN READING REAL
                RETURN
             END IF
          CASE DEFAULT
             status = 4 ! UNKNOWN SPECIFIER
             RETURN
       END SELECT

    END DO

    ! CLEAN UP
    IF (ASSOCIATED(sl1)) DEALLOCATE(sl1)
    IF (ASSOCIATED(sl2)) DEALLOCATE(sl2)

    status = 0 ! NO ERROR

  END SUBROUTINE parse_str
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
  ELEMENTAL FUNCTION tfrac(pi1, pm, pi2, ptp)

    ! I/O
    REAL(DP)             :: tfrac
    REAL(DP), INTENT(IN) :: pi1    ! interface pressure
    REAL(DP), INTENT(IN) :: pm     ! mid pressure
    REAL(DP), INTENT(IN) :: pi2    ! interface pressure
    REAL(DP), INTENT(IN) :: ptp    ! tropopause pressure

    ! LOCAL
    REAL(DP) :: pu             ! 'upper' pressure
    REAL(DP) :: pl             ! 'lower' pressure
    REAL(DP) :: pd             ! delta pressure
    INTEGER  :: ft             ! flag for troposphere
    INTEGER  :: fs             ! flag for stratosphere
    INTEGER  :: fe             ! flag for fs=ft=1
    INTEGER  :: fp             ! flag for tropopause
    
    ! pu < pm < pl
    pu = MIN(pi1, pi2)
    pl = MAX(pi1, pi2)

    ft = int(SIGN(-0.5_DP,(pm-ptp))+1._DP)  ! ft=1 T, =0 S
    fs = int(SIGN(-0.5_DP,(ptp-pm))+1._DP)  ! fs=0 T, =1 S
    fe = ft*fs
    ft = ft*(1-fe)
    fs = fs*(1-fe)

    fp = int(SIGN(-0.5_DP,(pl-ptp))+1._DP)*  &
         int(SIGN(-0.5_DP,(ptp-pu))+1._DP)
    ft = ft*(1-fp)
    fs = fs*(1-fp)

    pd = (pl - pu) * REAL(fp, DP) + 1.0 * REAL(1-fp, DP)

    tfrac = REAL(ft, DP) + REAL(fp, DP) * (pl-ptp)/pd

  END FUNCTION tfrac
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
!!$  ELEMENTAL SUBROUTINE int_14co( dt, p, q, oh, x, xs, xt, tf &
  ELEMENTAL SUBROUTINE int_14co( dt, c, q, oh, x, xs, xt, tf &
                               , y, ys, yt, l)

    ! 14CO MODULE ROUTINE (CORE)
    ! MAIN INTERFACE ROUTINE OF 14CO DIAGNOSTIC
    !
    ! PERFORMS FORCING INTEGRATION FOR ONE COLUMN/VECTOR
    !
    ! INPUT:
    !        dt         time step length [s]
!!$    !        p          pressure [Pa]
    !        c          concentration of air [cm^(-3)]
    !        q          14CO source strength [mol/mol/s]
    !        oh         OH concentration [cm^(-3)]
    !        x          tracer: 14CO               at t=t0
    !        xs         tracer: stratospheric 14CO at t=t0
    !        xt         tracer: tropospheric  14CO at t=t0
    !        tf         fraction of box in troposphere (e [0,1])
    !
    ! OUTPUT:
    !        y          14CO tendency for        t0-:-t0+dt  (?/s)
    !        ys         strat. 14CO tendency for t0-:-t0+dt  (?/s)
    !        yt         trop.  14CO tendency for t0-:-t0+dt  (?/s)
    !        l          14CO loss rate (14CO+OH)             (1/s)
    !        status     status flag
    !
    ! Note: [OH] must be cm^(-3)
    !       [q] and [x*] must be in compatible units,
    !           e.g. mol/mol/s and mol/mol
    !
    ! Author: Patrick Joeckel, MPICH, Jul 2004

    IMPLICIT NONE

    ! I/O
    REAL(DP), INTENT(IN)    :: dt  ! time step in seconds
!!$    REAL(DP), INTENT(IN)    :: p   ! pressure (in Pa)
    REAL(DP), INTENT(IN)    :: c   ! conc. of air (in cm^(-3))
    REAL(DP), INTENT(IN)    :: q   ! 14CO prod. rate
    REAL(DP), INTENT(IN)    :: oh  ! OH field (in cm^(-3))
    REAL(DP), INTENT(IN)    :: x   ! total 14CO
    REAL(DP), INTENT(IN)    :: xs  ! stratospheric 14CO
    REAL(DP), INTENT(IN)    :: xt  ! tropospheric  14CO
    REAL(DP), INTENT(IN)    :: tf  ! 
    REAL(DP), INTENT(INOUT) :: y   ! tendency of total 14CO (?/s)
    REAL(DP), INTENT(INOUT) :: ys  ! tendency of strat.14CO (?/s)
    REAL(DP), INTENT(INOUT) :: yt  ! tendency of trop. 14CO (?/s)
    REAL(DP), INTENT(INOUT) :: l   ! 14CO loss rate         (1/s)

    ! LOCAL
    REAL(DP)    :: rc              ! reaction coeff. (14CO + OH)
    REAL(DP)    :: ep              !

    ! REACTION RATE COEFFICIENT
!!$    rc = 1.5e-13_DP*(1._DP+0.6_DP*(p/101325._DP))
    rc = 1.57e-13_DP + c * 3.54e-33_DP

    ! LOSS RATE
    ! reset OH to O in case OH < 0 !!!
    l = rc * MAX(oh,0.0_DP)             ! [1/s]

    ! TENDENCIES:
    ! DIFFERENTIAL EQUATION (HERE ANALYTICALLY SOLVED)
    ! dc/dt = -(1/tau)*c + q ; c(0) = c0
    !    with  1/tau := rc*[OH] = l
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
    IF (l <= EPSILON) THEN  ! CASE 2
       y  = q
       ys = q * (1.0_DP-tf)
       yt = q * tf
    ELSE                   ! CASE 1
       ep = (exp(-dt*l)-1._DP)/dt
       y  = (x-q/l) * ep
       ys = (xs - (q*(1.0_DP-tf))/ l) * ep
       yt = (xt - (q*tf)/l) * ep
    END IF
    
    !qqq
    ! ADJUST FOR NEGATIVE CONCENTRATIONS

  END SUBROUTINE int_14co
! ----------------------------------------------------------------------

! ************************************************************************
! 14CO HELPER ROUTINES (PUBLIC)
! ************************************************************************
! -------------------------------------------------------------------------
  SUBROUTINE adj_tend(s, t, f1, t1, f2, t2, dt)

    ! 14CO MODULE ROUTINE (CORE)
    !
    ! LINEARIZE TWO-COMPONENT 'TAGGED' TRACER TENDENCIES:
    ! TRACER TENDENCIES ARE ADJUSTET TO FORCE:
    ! s + t*dt = (f1 + t1*dt) + (f2 + t2*dt)
    ! s is the 'sum' tracer, f1, f2 the 'fractional' tracers,
    ! and t, t1, and t2 the respective tendencies.
    !
    ! This can be used to correct for non-linearities resulting from
    ! tracer gradient dependent processes, such as, e.g., the advection.
    !
    ! Author: Patrick Joeckel, MPICH, Feb 2002

    IMPLICIT NONE

    ! I/O
    REAL(DP), DIMENSION(:), INTENT(IN)    :: s    ! sum tracer
    REAL(DP), DIMENSION(:), INTENT(IN)    :: t    ! tendency of sum tracer
    REAL(DP), DIMENSION(:), INTENT(IN)    :: f1   ! fraction 1 tracer
    REAL(DP), DIMENSION(:), INTENT(INOUT) :: t1   ! tend. of fract. 1 tracer
    REAL(DP), DIMENSION(:), INTENT(IN)    :: f2   ! fraction 2 tracer
    REAL(DP), DIMENSION(:), INTENT(INOUT) :: t2   ! tend. of fract. 2 tracer
    REAL(DP),               INTENT(IN)    :: dt   ! time step length

    ! LOCAL
    REAL(DP) :: a   ! correction factor
    REAL(DP) :: sum ! sum of tracers
    INTEGER  :: n   ! length of list
    INTEGER  :: k   ! loop counter

    n = SIZE(s)
    DO k=1, n
       ! s + t*dt =!= (f1 + t1*dt) + (f2 + t2*dt)
       ! a := (s+t*dt) / ((f1+t1*dt)+(f2+t2*dt))
       ! (f1+t1*dt)' := a*(f1+t1*dt)
       ! (f2+t2*dt)' := a*(f2+t2*dt)
       ! => (f1+t1*dt)' + (f2+t2*dt)' = (s+t*dt)   O.K.
       ! do not change tracers, rather adjust tendencies
       ! cond.: (f2' =!= f2) AND (f1' =!= f1)
       ! => f1 + t1'*dt = a*f1 + a*t1*dt ; f2 + t2'*dt = a*f2 + a*t2*dt
       ! => t1' = f1*(a-1)/dt + a*t1  ; t2' = f2*(a-1)/dt + a*t2
       sum = (f1(k)+t1(k)*dt + f2(k)+t2(k)*dt)
       IF (ABS(sum) >= EPSILON) THEN
          a = (s(k)+t(k)*dt)/sum
       ELSE
          a = 1.0_DP
       ENDIF
       t1(k) = t1(k)*a + f1(k)*(a-1.0_DP)/dt
       t2(k) = t2(k)*a + f2(k)*(a-1.0_DP)/dt
    END DO

  END SUBROUTINE adj_tend
! -------------------------------------------------------------------------

! -------------------------------------------------------------------------
  SUBROUTINE vmr2conc (f, p, t, kb)

    ! 14CO MODULE ROUTINE (CORE)
    !
    ! PERFORMS CONVERSION OF COLUMN FROM
    ! VOLUME MIXING RATIO (vmr = mol/mol)
    ! TO CONCENTRATION (conc) in 1/cm^3       !!!
    !
    ! conc = vmr * p/((kb*Na)*t) * Na * 1.0E-06
    !              mol/m^3       | conversion to cm^(-3)
    ! => conc = vmr * (p/kb*t) * 1.OE-06      in cm^(-3)
    !
    ! NEEDS PRESSURE (p) and TEMPERATURE (t)
    ! OPTIONAL: The values of the
    !           Boltzmann-constant (kb)
    !           can be specified by the subroutine call
    !           (for consistency to the calling model)
    ! Note: For the chosen unit (cm^(-3)), Na
    !           (Avogadro-constant) cancels out
    !
    ! Author: Patrick Joeckel, MPICH, Feb 2002

    IMPLICIT NONE

    ! I/O
    REAL(DP), DIMENSION(:), INTENT(INOUT) :: f  ! volume mixing ratio (mol/mol)
    REAL(DP), DIMENSION(:), INTENT(IN)    :: p  ! pressure (Pa)
    REAL(DP), DIMENSION(:), INTENT(IN)    :: t  ! temperature (K)
    REAL(DP), OPTIONAL,     INTENT(IN)    :: kb ! Boltzmann constant (J/K)

    ! LOCAL
    REAL(DP) :: lkb   ! Boltzmann constant (J/K)

    ! INITIALIZE
    IF (PRESENT(kb)) THEN
       lkb = kb
    ELSE
       lkb = 1.380662e-23_DP       ! Boltzmann constant in J/K
    END IF

    ! vmr <- conc
    f(:) = f(:) * (p(:)/(lkb*t(:))) * 1.0E-06_DP

  END SUBROUTINE vmr2conc
! -------------------------------------------------------------------------

! -------------------------------------------------------------------------
  ELEMENTAL SUBROUTINE merge_p(f1, f2, p, ptp, g)

    ! 14CO MODULE ROUTINE (CORE)
    !
    ! 'MERGING' OF 2 FIELDS (f1, f2) ON PRESSURE AXIS (p)
    !  g :=  f1 for p <= ptp
    !        f2 for p >  ptp
    !
    ! Author: Patrick Joeckel, MPICH, Jul 2004

    IMPLICIT NONE

    ! I/O
    REAL(DP), INTENT(IN)  :: f1, f2  ! INPUT COLUMNS
    REAL(DP), INTENT(IN)  :: p       ! PRESSURE COLUMN
    REAL(DP), INTENT(IN)  :: ptp     ! REF. PRESSURE
    REAL(DP), INTENT(OUT) :: g       ! OUTPUT COLUMN

    ! LOCAL
    INTEGER :: f        ! flag

    f = INT(SIGN(-0.5_DP,(ptp-p))+1._DP) ! f=1 ptp >= p(k), f=0 ptp < p(k)
    g = f1 * REAL(f,DP) + f2 * REAL(1-f,DP)
    
  END SUBROUTINE merge_p
! -------------------------------------------------------------------------

!!$! -------------------------------------------------------------------------
!!$  SUBROUTINE merge_p(f1, f2, p, ptp, g)
!!$
!!$    ! 14CO MODULE ROUTINE (CORE)
!!$    !
!!$    ! 'MERGING' OF 2 FIELD-COLUMNS (f1, f2) ON PRESSURE AXIS (p)
!!$    !  g :=  f1 for p <= ptp
!!$    !        f2 for p >  ptp
!!$    !
!!$    ! Author: Patrick Joeckel, MPICH, Feb 2002
!!$
!!$    IMPLICIT NONE
!!$
!!$    ! I/O
!!$    REAL(DP), DIMENSION(:), INTENT(IN)  :: f1, f2  ! INPUT COLUMNS
!!$    REAL(DP), DIMENSION(:), INTENT(IN)  :: p       ! PRESSURE COLUMN
!!$    REAL(DP)              , INTENT(IN)  :: ptp     ! REF. PRESSURE
!!$    REAL(DP), DIMENSION(:), INTENT(OUT) :: g       ! OUTPUT COLUMN
!!$
!!$    ! LOCAL
!!$    INTEGER :: nk       ! number of levels
!!$    INTEGER :: k        ! counter
!!$    INTEGER :: f        ! flag
!!$
!!$    ! INITIALIZE
!!$    nk = SIZE(p)
!!$
!!$    ! LOOP OVER LEVELS
!!$    DO k=1, nk
!!$       f=int(SIGN(-0.5_DP,(ptp-p(k)))+1._DP) ! f=1 ptp >= p(k), f=0 ptp < p(k)
!!$       g(k) = f1(k)*REAL(f,DP) + f2(k)*REAL(1-f,DP)
!!$    END DO
!!$
!!$  END SUBROUTINE merge_p
!!$! -------------------------------------------------------------------------

! -------------------------------------------------------------------------
  SUBROUTINE intpol_p(f, pi, p ,g, lg)

    ! 14CO MODULE ROUTINE (CORE)
    !
    ! INTERPOLATES ONE COLUMN OF AN INPUT FIELD (f) ON GIVEN PRESSURE
    ! COLUMN (pi) TO ANOTER PRESSURE COLUMN (p);
    ! OUTPUT: g
    !
    ! WITH THE OPTIONAL PARAMETER lg(=.TRUE.) THE INTERPOLATION CAN BE
    ! SWITCHED FROM LINEAR TO LOG-LINEAR
    !
    ! Author: Patrick Joeckel, MPICH, Feb 2002

    IMPLICIT NONE

    ! I/O
    REAL(DP), DIMENSION(:), INTENT(IN) :: f    ! field to be interpolated
    REAL(DP), DIMENSION(:), INTENT(IN) :: pi   ! pressure column of f
    REAL(DP), DIMENSION(:), INTENT(IN) :: p    ! press. column to interpol. on
    REAL(DP), DIMENSION(:), INTENT(INOUT) :: g ! interpolated field
    LOGICAL, OPTIONAL ,     INTENT(IN) :: lg   ! flag for log-interpolation

    ! LOCAL
    INTEGER  :: nko                     ! no. of levels in output column
    INTEGER  :: k                       ! level counter
    INTEGER  :: k1, kn, kd, ki          ! level bounds for interpol.
    INTEGER  :: kdd
    REAL(DP) :: delp                    ! delta pressure
    INTEGER  :: fxt                     ! flag for interpolation loop
    LOGICAL  :: llg                     ! flag for log-interpolation

    ! INITIALIZE
    IF (PRESENT(lg)) THEN
       llg = lg
    ELSE
       llg = .FALSE.   ! DEFAULT: linear interpolation
    END IF
    nko = SIZE(p)
    ! INITIAL CONDITION FOR INTERPOLATION LOOP
    CALL intpol_init(pi, k1, kn, kd)
    ki = k1

    ! NOTE: FOR BETTER PERFORMANCE (EXPENSIVE 'IF') THE
    ! LINEAR/LOG QUESTION IS EVALUATED OUTSIDE THE LOOP

    IF (llg) THEN                  ! LOG - INTERPOLATION

       ! LOOP OVER OUTPUT LEVELS
       DO k=1, nko
          ! FIND APPROPRIATE INTERVAL IN INPUT COLUMN
          kdd = kd
          DO
             CALL intpol_step(p(k), pi, k1, kn, kdd, ki, delp, fxt)
             IF (fxt == 1) EXIT
          END DO      ! END FIND APPROPRIATE INTERVAL IN INPUT COLUMN
          ! LOG-INTERPOLATION
          g(k) = 10**(log10(f(ki)) + delp*(log10(f(ki+kdd))-log10(f(ki))))
       END DO ! LOOP OVER OUTPUT LEVELS

    ELSE                          ! LINEAR INTERPOLATION

       ! LOOP OVER OUTPUT LEVELS
       DO k=1, nko
          ! FIND APPROPRIATE INTERVAL IN INPUT COLUMN
          kdd = kd
          DO
             CALL intpol_step(p(k), pi, k1, kn, kdd, ki, delp, fxt)
             IF (fxt == 1) EXIT
          END DO      ! END FIND APPROPRIATE INTERVAL IN INPUT COLUMN
          ! LINEAR INTERPOLATION
          g(k) = f(ki) + delp*(f(ki+kdd)-f(ki))
       END DO ! LOOP OVER OUTPUT LEVELS

    END IF                        ! INTERPOLATION MODE

  END SUBROUTINE intpol_p
! -------------------------------------------------------------------------

! ************************************************************************
! 14CO HELPER ROUTINES (PRIVATE)
! ************************************************************************
! -------------------------------------------------------------------------
  SUBROUTINE intpol_init(p, k1, kn, kd)

    ! 14CO MODULE ROUTINE (HELPER)
    !
    ! THIS SUBROUTINE INITIALIZES THE LOOP VARIABLES
    ! START (k1 = 1 OR n), STOP (kn = n OR 1), and DELTA (kd = +1 OR -1)
    ! FOR STEPPING THROUGH AN ORDERED (INCREASING OR DECREASING) LIST (p)
    !
    ! Author: Patrick Joeckel, MPICH, Feb 2002

    IMPLICIT NONE

    ! I/O
    REAL(DP), DIMENSION(:), INTENT(IN)  :: p ! ordered list
                            ! start, stop, delta for loop
    INTEGER,                INTENT(OUT) :: k1, kn, kd

    ! LOCAL
    INTEGER     :: i1
    INTEGER     :: n           ! length of ordered list

    n = SIZE(p)

    i1 = int(SIGN(-0.5_DP,(p(n)-p(1)))+1._DP)    ! =1 for p(n)  >= p(1)

    k1 = i1   + (1-i1)*n
    kn = i1*n + (1-i1)
    kd = i1   - (1-i1)

  END SUBROUTINE intpol_init
! -------------------------------------------------------------------------

! -------------------------------------------------------------------------
  SUBROUTINE intpol_step(p, pi, k1, kn, kd, ki, dpf, fex)

    ! 14CO MODULE ROUTINE (HELPER)
    !
    ! THIS SUBROUTINE CHECKS IF p IS IN THE INTERVAL
    ! [pi(ks), pi(ks+kd)] OF AN ORDERED LIST (pi) AND
    ! CALCULATES THE INTERVAL FRACTION dpf FOR LINEAR INTERPOLATION
    ! IF p IS NOT IN [pi(ks), pi(ks+kd)] THEN
    ! THE LOOP INDEX ks IS IN-/DECREMENTED
    ! BY kd (i.e., THE INTERVAL IS SHIFTED)
    ! BEYOND THE INTERVAL SERIES BOUNDS (ks+kd<1 OR ks+kd>kn) THE
    ! VLAUES ARE CONTINUED
    ! A FLAG FOR BREAKING/CONTINUING THE LOOP IS RETURNED (fex)
    !
    !  <---- a   a   b   c   d   e   f   g   h   i   j   j---->
    ! *********|---|---|---|---|---|---|---|---|---|---|**********
    !          |               |   |                   |
    !          1               ks  ks+kd               kn
    !
    ! Author: Patrick Joeckel, MPICH, Feb 2002

    IMPLICIT NONE

    ! I/O
    REAL(DP)              , INTENT(IN)    :: p       ! reference value
    REAL(DP), DIMENSION(:), INTENT(IN)    :: pi      ! ordered list
    INTEGER               , INTENT(IN)    :: k1      ! start index
    INTEGER               , INTENT(IN)    :: kn      ! end index
    INTEGER               , INTENT(INOUT) :: kd      ! delta-index
    INTEGER               , INTENT(INOUT) :: ki      ! current index
    REAL(DP)              , INTENT(OUT)   :: dpf     ! interval fraction
    INTEGER               , INTENT(OUT)   :: fex     ! exit flag

    ! LOCAL
    INTEGER :: n
    INTEGER :: ii1, ii2, ii3, ii4, ii5, ii6           ! flags

    n = SIZE(pi)

    ! RANGE CHECK
    ii1 = INT(SIGN(-0.5_DP,REAL((ki+kd-1),DP))+1._DP)   ! =1 for ki+kd >= 1
    ii2 = INT(SIGN(-0.5_DP,REAL((n-(ki+kd)),DP))+1._DP) ! =1 for n >= ki+kd
    IF (ii1*ii2 == 1) THEN
       ii3 = INT(SIGN(-0.5_DP,(p-pi(ki)))+1._DP)     ! =1 for p  >= p(ki)
       ii4 = INT(SIGN(-0.5_DP,(pi(ki+kd)-p))+1._DP)  ! =1 for p(ki+kd) >= p
       ii5 = 1-MIN(1,ABS(ki-k1))                  ! =1 for ki = k1
       dpf = ((p-pi(ki))/(pi(ki+kd)-pi(ki)))*REAL(ii3*ii4,DP)
       fex = ii3*ii4 + ii5*(1-ii3)*ii4
       ki = ki + kd * (1-fex)
    ELSE
       ii3 = INT(SIGN(-0.5_DP,(pi(ki)-p))+1._DP)    ! =1 for p(ki)  >= p
       ii4 = INT(SIGN(-0.5_DP,(p-pi(ki)))+1._DP)    ! =1 for p  >= p(ki)
       ii5 = 1-MIN(1,ABS(ki-k1))                 ! =1 for ki = k1
       ii6 = 1-MIN(1,ABS(ki-kn))                 ! =1 for ki = kn
       dpf = 0.0_DP
       fex = MIN(1, ii4*ii6 + ii3*ii5)
       kd = kd * (1-fex)
    END IF

  END SUBROUTINE intpol_step
! -------------------------------------------------------------------------

!!$! -------------------------------------------------------------------------
!!$  SUBROUTINE h2p(a,b,p0,ps,p)
!!$
!!$    ! 14CO MODULE ROUTINE (HELPER)
!!$    !
!!$    ! CALCULATES PRESSURE FROM HYBRID-COORDINATES AND SURFACE PRESSURE
!!$    ! FOR ONE COLUMN
!!$    !
!!$    ! Author: Patrick Joeckel, MPICH, Feb 2002
!!$
!!$    IMPLICIT NONE
!!$
!!$    ! I/O
!!$    REAL(DP), DIMENSION(:), INTENT(IN)  :: a, b  ! hybrid a and b coefficients
!!$    REAL(DP)              , INTENT(IN)  :: p0    ! reference pressure
!!$    REAL(DP)              , INTENT(IN)  :: ps    ! surface pressure
!!$    REAL(DP), DIMENSION(:), INTENT(OUT) :: p     ! pressure
!!$
!!$    p(:) = (a(:)*p0 + b(:)*ps)
!!$
!!$  END SUBROUTINE h2p
!!$! -------------------------------------------------------------------------

! -------------------------------------------------------------------------
  SUBROUTINE mgs2vmrs(f, mwair, na)

    ! 14CO MODULE ROUTINE (HELPER)
    !
    ! CONVERTS 14CO source unit  (1 column)
    ! from molec/(g s) to mol/mol/s
    !
    ! Note: The required constants
    !            mole-weight of air (mwair) in g/mol
    !       and  Avogadro-Constant  (na)    in molec/mol
    !       can optionally be passed from outside.
    !
    ! Author: Patrick Joeckel, MPICH, Feb 2002

    IMPLICIT NONE

    ! I/O
    REAL(DP), DIMENSION(:), INTENT(INOUT) :: f
    REAL(DP), OPTIONAL    , INTENT(IN)    :: mwair  ! molarmass of dry air
                                                    ! (g/mol)
    REAL(DP), OPTIONAL    , INTENT(IN)    :: na     ! Avogadro constant
                                                    ! (molec/mol)

    ! LOCAL
    REAL(DP) :: zmwair   ! molarmass of dry air
    REAL(DP) :: zna      ! Avogadro constant

    IF (PRESENT(mwair)) THEN
       zmwair = mwair
    ELSE
       zmwair = 28.970_DP      ! g/mol
    END IF

    IF (PRESENT(na)) THEN
       zna = na
    ELSE
       zna = 6.022045E+23_DP   ! molec/mol
    END IF

    f(:) = f(:) * (zmwair/zna)

  END SUBROUTINE mgs2vmrs
! -------------------------------------------------------------------------

! **********************************************************************
END MODULE messy_d14co
! **********************************************************************

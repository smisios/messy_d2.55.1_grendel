! **********************************************************************
MODULE messy_isopcor
! **********************************************************************

#define _EXTERNAL_LT

  ! REFERENCES:
  !
  ! S. Poluianov, G.A. Kovaltsov, A.L. Mishev, and I.G. Usoskin,
  !    Production of cosmogenic isotopes 7Be, 10Be, 14C, 22Na and 36Cl
  !    in the atmosphere: Altitudinal profiles of yield functions,
  !    JGR, 2016
  !

  ! ----------- >

  USE messy_main_constants_mem, ONLY: DP, pi, g

  IMPLICIT NONE
  INTRINSIC :: RESHAPE
  SAVE

  ! ----------- <

  ! GLOBAL PARAMETER
  CHARACTER(LEN=*), PARAMETER :: modstr = 'isopcor' ! name of module
  CHARACTER(LEN=*), PARAMETER :: modver = '0.9'     ! module version

  !INTEGER,  PARAMETER :: dp = SELECTED_REAL_KIND(12,307)
  !REAL(dp), PARAMETER :: pi = 3.14159265358979323846_dp
  !REAL(dp), PARAMETER :: g  = 9.81_dp ! [m s^(-2)] Earth accelaration

#ifndef _EXTERNAL_LT
  INTEGER, PARAMETER :: NE   =  23 ! number of energy bins
  INTEGER, PARAMETER :: NH   = 110 ! number of height bins
#else
  INTEGER :: NE   ! number of energy bins
  INTEGER :: NH   ! number of height bins
#endif
  INTEGER, PARAMETER :: NISO =   5 ! number of isotopes


  INTEGER, PARAMETER :: I_ISO_14C  = 1
  INTEGER, PARAMETER :: I_ISO_10Be = 2
  INTEGER, PARAMETER :: I_ISO_7Be  = 3
  INTEGER, PARAMETER :: I_ISO_22Na = 4
  INTEGER, PARAMETER :: I_ISO_36Cl = 5

  CHARACTER(LEN=4), DIMENSION(NISO), PARAMETER :: isoname = &
       (/'14C ', '10Be', '7Be ', '22Na', '36Cl' /)

#ifndef _EXTERNAL_LT
  INTEGER, PARAMETER :: NCR  =   2 ! protons and alpha
#else
  INTEGER :: NCR ! protons and alpha
#endif
  INTEGER, PARAMETER :: I_CR_PROTON = 1
  INTEGER, PARAMETER :: I_CR_ALPHA  = 2

#ifndef _EXTERNAL_LT
  ! ENERGY BINS [GeV/nuc]
  REAL(dp), DIMENSION(NE) :: energy = REAL( (/ &
       0.02, 0.025, 0.032, 0.04, &
       0.05, 0.063, 0.079, 0.1, 0.126, 0.159, 0.2, 0.251, 0.3, 0.5, &
       1.0, 2.0, 3.0, 5.0, 10.0, 20.0, 30.0, 50.0, 100.0 /), dp)

  ! ATMOSPHERIC DEPTH [g/cm^2]
  REAL(dp), DIMENSION(NH) :: depth = REAL( (/ &
         0.0,  10.0,  15.0,  20.0,  25.0,  30.0,  35.0,  40.0,  45.0,  50.0, &
        55.0,  60.0,  65.0,  70.0,  75.0,  80.0,  85.0,  90.0,  95.0, 100.0, &
       110.0, 120.0, 130.0, 140.0, 150.0, 160.0, 170.0, 180.0, 190.0, 200.0, &
       210.0, 220.0, 230.0, 240.0, 250.0, 260.0, 270.0, 280.0, 290.0, 300.0, &
       310.0, 320.0, 330.0, 340.0, 350.0, 360.0, 370.0, 380.0, 390.0, 400.0, &
       410.0, 420.0, 430.0, 440.0, 450.0, 460.0, 470.0, 480.0, 490.0, 500.0, &
       510.0, 520.0, 530.0, 540.0, 550.0, 560.0, 570.0, 580.0, 590.0, 600.0, &
       610.0, 620.0, 630.0, 640.0, 650.0, 660.0, 670.0, 680.0, 690.0, 700.0, &
       710.0, 720.0, 730.0, 740.0, 750.0, 760.0, 770.0, 780.0, 790.0, 800.0, &
       810.0, 820.0, 830.0, 840.0, 850.0, 860.0, 870.0, 880.0, 890.0, 900.0, &
       910.0, 920.0, 930.0, 940.0, 950.0, 960.0, 970.0, 980.0, 990.0,1000.0 /) &
       , dp)

  INCLUDE 'messy_isopcor_data.inc'

!!$  ! production function (in atoms of the isotope per one incident
!!$  !   nucleon in gram of air at the given atmospheric depth)
!!$
!!$  REAL, DIMENSION(NE,NH,NCR,NISO), PARAMETER :: S = RESHAPE( (/   &
!!$       S_14C_p,  S_14C_a,  S_10Be_p, S_10Be_a, S_7Be_p, S_7Be_a,  &
!!$       S_22Na_p, S_22Na_a, S_36Cl_p, S_36Cl_a /),                 &
!!$       SHAPE=(/ NE, NH, NCR, NISO /) )
!!$

  ! yield function [atoms g^(-1) cm^2 sr)] ; Y = S * pi
  REAL(DP), DIMENSION(NE,NH,NCR,NISO), PARAMETER :: Y = REAL( RESHAPE( (/ &
       S_14C_p,  S_14C_a,  S_10Be_p, S_10Be_a, S_7Be_p, S_7Be_a,  &
       S_22Na_p, S_22Na_a, S_36Cl_p, S_36Cl_a /),                 &
       SHAPE=(/ NE, NH, NCR, NISO /) ), dp) * pi
#else
  REAL(dp), DIMENSION(:),       POINTER :: energy => NULL()
  REAL(dp), DIMENSION(:),       POINTER :: depth  => NULL()
  REAL(DP), DIMENSION(:,:,:,:), POINTER :: Y      => NULL()
#endif

  ! The yield function, Y(E,h), is defined as the production
  ! (the number of atoms per gram of air) of the isotope, at given
  ! atmospheric depth h, by primary particles of type i with the unit
  ! intensity (one primary particle with kinetic energy per nucleon E in the
  ! interplanetary space per steradian and cm^2 ). The units of Y are
  ! [atoms g^(-1) cm^2 sr].
  !
  ! The production rate Q of cosmogenic isotope at time t is then defined
  ! as an integral of the product of the yield function and the energy
  ! spectrum of cosmic rays J_i(E,t) ( [sr sec cm^2 ]^(-1) ), above the
  ! energy E_c corresponding to the local geomagnetic rigidity cutoff P_c:
  !
  !                        infinity
  !                 --    --
  !                 \     /
  ! Q(t, h, P_c ) = /    /  Y_i(E,h) J_i(E,t) dE,      (1)
  !                 -- --
  !                 i   E_c,i
  !
  ! where the summation is over different types of primary cosmic ray
  ! particles (protons, alpha-particles, etc.). The relation between E_c,i
  ! and P_c (defined independently of the yield function computations) is
  !
  !             (      --------------------       )
  !             (     /    ( Z_i * P_c )^2 |      )
  ! E_c,i = E_r (    / 1 + ( --------- )     - 1) ),  (2)
  !             (  \/      ( A_i * E_r )          )
  !
  ! where Z_i and A_i are the charge and mass numbers of particles,
  ! respectively, E_r = 0.938 GeV is the rest mass of a proton.
  ! For computations of the yield function we considered, as primary particles,
  ! only protons and alpha-particles. Species heavier than helium can be
  ! effectively considered as scaled (by the nucleonic number) alpha-particles
  ! [see Webber and Higbie, 2003].
  !
  ! #########################################################################
  ! Recipe for computation of the cosmogenic isotope production
  ! #########################################################################
  !
  ! Here we present a recipe on how to compute the production rate Q(h,P_c,t)
  ! of a cosmogenic isotope at a given location and time. The location is
  ! defined by the atmospheric depth h and the local geomagnetic rigidity
  ! cutoff P_c.
  !
  ! First, the yield function for a cosmic ray specie i (proton or
  ! alpha-particle) should be computed for the given atmospheric depth as
  ! Y_i(h,E) = pi ? S_i(h,E), where S_i(h,E) is taken from an appropriate
  ! table in the Supporting information. Note that the energy of
  ! alpha-particles should be taken as kinetic energy per nucleon.
  !
  ! Next, the production rate of the isotope should be computed using
  ! Equation (1), where the cutoff energy is calculated from the local
  ! geomagnetic rigidity cutoff P_c using formula (2). For numerical
  ! integration, values of Y can be interpolated by a power-law function
  ! between the tabulated points.
  !
  ! The value of P_c as well as spectra J_i of cosmic-ray protons and
  ! alpha-particles should be known independently.
  ! For the spectra we recommend using the force-field approximation,
  ! where spectra are parameterized via a single parameter, the modulation
  ! potential PHI (see details in Usoskin et al. [2005]).
  ! Values for the modulation potential are given, e.g., by
  ! Usoskin et al. [2011] and updated at http://cosmicrays.oulu.fi/phi/phi.html.
  !


CONTAINS

! ==========================================================================
  ELEMENTAL SUBROUTINE geomag(lon, lat, Ecp, Eca, Pc)

    IMPLICIT NONE
    INTRINSIC :: SIN, COS, SQRT

    ! I/O
    REAL(DP), INTENT(IN)  :: lon ! [deg] geographic longitude
    REAL(DP), INTENT(IN)  :: lat ! [deg] geographic latitude
    ! Ec is energy corresponding to the local geomagnetic rigidity cutoff Pc
    REAL(DP), INTENT(OUT) :: Ecp ! [GeV] Ec of protons
    REAL(DP), INTENT(OUT) :: Eca ! [GeV] Ec of alpha particles
    REAL(DP), INTENT(OUT) :: Pc  ! [GeV] geomagnetic rigidity cutoff

    ! LOCAL
    REAL(dp), PARAMETER :: Er     = 0.938_dp ! [GeV] rest mass of a proton
    REAL(DP), PARAMETER :: Xo     = -6.29E-02_dp
    REAL(DP), PARAMETER :: Yo     =  4.98E-02_dp
    REAL(DP), PARAMETER :: Zo     =  3.28E-02_dp
    REAL(DP), PARAMETER :: cosfi  =  3.13E-01_dp
    REAL(DP), PARAMETER :: sinfi  = -9.50E-01_dp
    REAL(DP), PARAMETER :: costet = 9.84E-01_dp
    REAL(DP), PARAMETER :: sintet = 1.78E-01_dp
    REAL(DP), PARAMETER :: MM     = 7.76E+00_dp
    !
    REAL(DP) :: Coe, l, m, n, del2, tett, fit, lt, mt, nt, colg, A, bbb

    Coe=pi/180.0_dp;
    l=sintet*cosfi
    m=sintet*sinfi
    n=costet

    del2=Xo*Xo + Yo*Yo + Zo*Zo
    tett=(90.0_dp-lat)*Coe
    fit=lon*Coe

    lt = SIN(tett)*COS(fit)
    mt = SIN(tett)*SIN(fit)
    nt = COS(tett)
    colg = l*lt + m*mt + n*nt
    A = Xo*l + Yo*m + Zo*n
    bbb = SQRT(1.0_dp + del2 - 2*(Xo*lt +Yo*mt +Zo*nt))
    colg = (colg-A)/bbb

    ! Geom. cutoff Pc in GV
    Pc = 1.9_dp*MM*(1-colg*colg)*(1-colg*colg)/bbb/bbb

    ! see Eqn. (2):
    ! proton: Z_i = 1; A_i = 1
    Ecp = SQRT(Pc*Pc+Er*Er)-Er
    ! alpha: Z_i = 2; A_i = 4
    Eca = SQRT(0.25_dp*Pc*Pc+Er*Er)-Er

  END SUBROUTINE geomag
! ==========================================================================

! ==========================================================================
  SUBROUTINE gcr_spectrum(phi, Jp, Ja)

    IMPLICIT NONE
    INTRINSIC :: SQRT

    ! I/O
    REAL(DP), INTENT(IN)                 :: phi    ! modulation potential [GV]
    REAL(DP), DIMENSION(NE), INTENT(OUT) :: Ja, Jp

    ! INTEGER
    INTEGER  :: j
    REAL(DP) :: phiA, Tp, Pp, Ta, Pa

    phiA=phi*0.5_dp

    DO j = 1, NE
       !=== protons
       Tp = energy(j) + phi
       Pp = SQRT(Tp*(Tp+1.876_dp))
       Jp(j) = 1.9_dp*(Pp**(-2.78_dp))/(1.0_dp+0.4866_dp*(Pp**(-2.51_dp)))
       Jp(j) = Jp(j)*energy(j)*(energy(j)+1.876_dp)/Tp/(Tp+1.876_dp)

       !===  alpha-particles
       Ta = energy(j) + phiA
       Pa = SQRT(Ta*(Ta+1.876_dp))
       Ja(j) = 0.57_dp*(Pa**(-2.78_dp))/(1.0_dp+0.4866_dp*(Pa**(-2.51_dp)))
       Ja(j) = Ja(j)*energy(j)*(energy(j)+1.876_dp)/Ta/(Ta+1.876_dp)
    END DO

  END SUBROUTINE gcr_spectrum
! ==========================================================================

! ==========================================================================
  SUBROUTINE yield(h, Jp, Ja, Fp, Fa)

    IMPLICIT NONE
    INTRINSIC :: MIN, MAX

    ! I/O
    REAL(DP), INTENT(IN)                      :: h ! atmospheric depth [g/cm2]
    REAL(DP), DIMENSION(NE),      INTENT(IN)  :: Jp, Ja
    REAL(DP), DIMENSION(NE,NISO), INTENT(OUT) :: Fp, Fa

    ! LOCAL
    REAL(DP), PARAMETER :: Ymin = 1.0E-20_dp
    INTEGER  :: II, i, jiso
    REAL(DP) :: w1, w2
    REAL(DP), DIMENSION(NE,NISO) :: YFp, YFa
    REAL(DP) :: zh
    REAL(DP), DIMENSION(NE,NH)   :: Yp
    REAL(DP), DIMENSION(NE,NH)   :: Ya

    zh = MIN(h,depth(NH)) ! no extrapolation

    II = NH - 1
    DO i = 1, NH-1
       IF ((zh >= depth(i)) .AND. (zh < depth(i+1))) THEN
          II=i
          exit
       END IF
    END DO

    w1 = (zh-depth(II)) / (depth(II+1)-depth(II))
    w2 = (depth(II+1)-zh) / (depth(II+1)-depth(II))

    DO jiso = 1, NISO

       ! Note: The limitation to 1E-20 as minimum (i.e. non-zero) is chosen
       !       to allow an efficient yield interpolation (SUBROUTINE yield)!
       Yp = MAX(Y(:,:, I_CR_PROTON, jiso), Ymin)
       Ya = MAX(Y(:,:, I_CR_ALPHA,  jiso), Ymin)

       YFp(:,jiso) = w2*Yp(:,II) + w1*Yp(:,II+1)
       YFa(:,jiso) = w2*Ya(:,II) + w1*Ya(:,II+1)

       Fp(:,jiso) = Jp(:) * YFp(:,jiso)
       Fa(:,jiso) = Ja(:) * YFa(:,jiso)

    END DO

  END SUBROUTINE yield
! ==========================================================================

! ==========================================================================
  SUBROUTINE integrate(Fp, Fa, Ecp, Eca, Q)

    IMPLICIT NONE
    INTRINSIC :: LOG, EXP

    ! I/O
    REAL(dp), DIMENSION(NE), INTENT(IN)  :: Fp, Fa
    REAL(dp),                INTENT(IN)  :: Ecp, Eca
    REAL(dp),                INTENT(OUT) :: Q

    ! LOCAL
    REAL(dp) :: Sump, Suma
    INTEGER  :: j
    REAL(dp) :: b, a, sup
    REAL(dp) :: zFpj, zFpjm1
    REAL(dp) :: zFaj, zFajm1

    Sump = 0.0_dp
    Suma = 0.0_dp

    !====== protons ====
    DO j = 2, NE
       b = LOG(Fp(j)/Fp(j-1))/LOG(energy(j)/energy(j-1))
       a = EXP(LOG(Fp(j))-b*LOG(energy(j)))
       IF ( energy(j) < Ecp ) THEN
          sup = 0.0_dp
       ELSEIF (energy(j-1) < Ecp) THEN
          sup = a/(b+1)*(energy(j)**(b+1)-Ecp**(b+1))
       ELSE
          sup = a/(b+1)*(energy(j)**(b+1)-energy(j-1)**(b+1))
       END IF
       Sump = Sump + sup
    END DO

    Sump = Sump - a/(b+1.0_dp) * energy(NE)**(b+1.0_dp)

    !====== alphas ====
    DO j = 2,NE
       b = LOG(Fa(j)/Fa(j-1))/LOG(energy(j)/energy(j-1))
       a = EXP(LOG(Fa(j))-b*LOG(energy(j)))
       IF (energy(j) < Eca) THEN
          sup = 0.0_dp
       ELSEIF (energy(j-1) < Eca) THEN
          sup = a/(b+1)*(energy(j)**(b+1)-Eca**(b+1))
       ELSE
          sup = a/(b+1)*(energy(j)**(b+1)-energy(j-1)**(b+1))
       END IF
       Suma = Suma + sup
    END DO

    Suma = Suma - a/(b+1.0_dp) * energy(NE)**(b+1.0_dp)

    Q = Suma + Sump

  END SUBROUTINE integrate
! ==========================================================================

! -------------------------------------------------------------------------
  SUBROUTINE adj_tend(s, t, f1, t1, f2, t2, dt)

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
    INTRINSIC :: TINY

    ! I/O
    REAL(DP), DIMENSION(:), INTENT(IN)    :: s    ! sum tracer
    REAL(DP), DIMENSION(:), INTENT(IN)    :: t    ! tendency of sum tracer
    REAL(DP), DIMENSION(:), INTENT(IN)    :: f1   ! fraction 1 tracer
    REAL(DP), DIMENSION(:), INTENT(INOUT) :: t1   ! tend. of fract. 1 tracer
    REAL(DP), DIMENSION(:), INTENT(IN)    :: f2   ! fraction 2 tracer
    REAL(DP), DIMENSION(:), INTENT(INOUT) :: t2   ! tend. of fract. 2 tracer
    REAL(DP),               INTENT(IN)    :: dt   ! time step length

    ! LOCAL
    REAL(DP),         PARAMETER :: EPSILON = TINY(0.0_DP)    ! zero
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

! **********************************************************************
END MODULE messy_isopcor
! **********************************************************************

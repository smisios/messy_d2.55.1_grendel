  SUBROUTINE jval_cal_CH3COCH3

    ! CH3COCH3 -> products Gierzcak et al., CP, (1998)

    INTEGER :: j, k
    REAL    :: dj
    REAL, DIMENSION(0:MAXWAV) :: sig_CH3COCH3

    REAL :: c1, c2, f2, f3, f4, f5
    REAL :: drel
    ! factors that correct the Jval depending on  the "new" quantum yield
    REAL :: nqyf(MAXWAV)
    REAL :: bqyf(MAXWAV)
    REAL :: wqyf(MAXWAV)
    REAL :: gqyf(MAXWAV)
    REAL :: A,B,a0,b0,D0,a1,b1,D1,a2,b2,D2,a3,b3,D3,E3,A4,B4,D4
    REAL :: PHI_ch3co, PHI_co, PhiS
    REAL, PARAMETER :: wave(7) = (/ &
      2.0513E-05, 2.8777E-05, 3.0200E-05, 3.0900E-05, &
      3.2000E-05, 3.7000E-05, 5.8000E-05 /)

    ! T-const parameters:
    REAL, PARAMETER :: a2_CH3COCH3(dim55) = (/ &
      1.0640E-22, 7.0800E-23, 3.8400E-23, 1.2800E-23, -4.3996E-24, -1.4800E-23, &
      -2.0400E-23, -2.2000E-23, -2.2400E-23, -2.1600E-23, -2.0000E-23, -1.8400E-23, &
      -1.6800E-23, -1.5199E-23, -1.4000E-23, -1.2400E-23, -1.1601E-23, -1.0400E-23, &
      -9.5998E-24, -8.7998E-24, -8.4004E-24, -7.1997E-24, -7.2003E-24, -6.3997E-24, &
      -6.0003E-24, -5.5996E-24, -5.2003E-24, -4.8002E-24, -4.3996E-24, -4.4002E-24, &
      -4.0002E-24, -3.5995E-24, -3.2002E-24, -3.2002E-24, -3.2002E-24, -2.7995E-24, &
      -2.4001E-24, -2.8001E-24, -2.0001E-24, -2.4001E-24, -1.7600E-24, -1.2801E-24, &
      -9.5992E-25, -6.4003E-25, -4.8002E-25, -4.0002E-25, -2.4001E-25, -1.6001E-25, &
      -1.5988E-25, -8.0004E-26, -8.0004E-26, -0.0000E+00, -8.0004E-26, -0.0000E+00, &
      -0.0000E+00 /)
    REAL, PARAMETER :: b2_CH3COCH3(dim55) = (/ &
      1.9411E-20, 1.9518E-20, 1.9696E-20, 1.9901E-20, 2.0081E-20, 2.0216E-20, &
      2.0303E-20, 2.0332E-20, 2.0340E-20, 2.0322E-20, 2.0281E-20, 2.0236E-20, &
      2.0187E-20, 2.0135E-20, 2.0092E-20, 2.0031E-20, 1.9999E-20, 1.9947E-20, &
      1.9911E-20, 1.9872E-20, 1.9852E-20, 1.9789E-20, 1.9789E-20, 1.9742E-20, &
      1.9718E-20, 1.9693E-20, 1.9667E-20, 1.9639E-20, 1.9611E-20, 1.9611E-20, &
      1.9581E-20, 1.9550E-20, 1.9518E-20, 1.9518E-20, 1.9518E-20, 1.9482E-20, &
      1.9446E-20, 1.9483E-20, 1.9407E-20, 1.9446E-20, 1.9382E-20, 1.9328E-20, &
      1.9287E-20, 1.9243E-20, 1.9219E-20, 1.9206E-20, 1.9178E-20, 1.9163E-20, &
      1.9163E-20, 1.9146E-20, 1.9146E-20, 1.9127E-20, 1.9147E-20, 1.9126E-20, &
      1.9126E-20 /)
    REAL, PARAMETER :: a3_CH3COCH3(dim55) = (/ &
      -4.9569E-23, -4.1283E-23, -3.2563E-23, -2.4847E-23, -1.8694E-23, -1.4097E-23, &
      -1.0771E-23, -8.3866E-24, -6.6773E-24, -5.4393E-24, -4.5146E-24, -3.8179E-24, &
      -3.2668E-24, -2.8355E-24, -2.4890E-24, -2.1965E-24, -1.9569E-24, -1.7572E-24, &
      -1.5861E-24, -1.4337E-24, -1.3099E-24, -1.1981E-24, -1.1027E-24, -1.0144E-24, &
      -9.4249E-25, -8.7061E-25, -8.1501E-25, -7.5879E-25, -7.1086E-25, -6.6694E-25, &
      -6.2724E-25, -5.9105E-25, -5.5512E-25, -5.2715E-25, -4.9940E-25, -4.7125E-25, &
      -4.4728E-25, -4.2732E-25, -4.0751E-25, -3.8738E-25, -3.3786E-25, -2.7476E-25, &
      -2.2764E-25, -1.9089E-25, -1.6307E-25, -1.3978E-25, -1.2061E-25, -1.0463E-25, &
      -9.1054E-26, -7.9873E-26, -7.0287E-26, -6.0702E-26, -5.4357E-26, -4.7924E-26, &
      -4.1533E-26 /)
    REAL, PARAMETER :: b3_CH3COCH3(dim55) = (/ &
      8.4433E-21, 8.2317E-21, 7.7906E-21, 7.2073E-21, 6.5880E-21, 6.0102E-21, &
      5.5088E-21, 5.0897E-21, 4.7465E-21, 4.4669E-21, 4.2349E-21, 4.0427E-21, &
      3.8768E-21, 3.7362E-21, 3.6146E-21, 3.5046E-21, 3.4085E-21, 3.3234E-21, &
      3.2462E-21, 3.1736E-21, 3.1116E-21, 3.0527E-21, 3.0001E-21, 2.9492E-21, &
      2.9060E-21, 2.8610E-21, 2.8248E-21, 2.7867E-21, 2.7531E-21, 2.7212E-21, &
      2.6914E-21, 2.6632E-21, 2.6344E-21, 2.6113E-21, 2.5877E-21, 2.5630E-21, &
      2.5414E-21, 2.5229E-21, 2.5040E-21, 2.4844E-21, 2.4347E-21, 2.3636E-21, &
      2.3046E-21, 2.2540E-21, 2.2122E-21, 2.1742E-21, 2.1406E-21, 2.1106E-21, &
      2.0834E-21, 2.0596E-21, 2.0380E-21, 2.0152E-21, 1.9993E-21, 1.9824E-21, &
      1.9648E-21 /)
    REAL, PARAMETER :: dsa3_1_CH3COCH3(dim55) = (/ &
      -8.7405E-04, -8.6062E-04, -7.8720E-04, -6.7922E-04, -5.6261E-04, -4.5646E-04, &
      -3.6784E-04, -2.9810E-04, -2.4464E-04, -2.0400E-04, -1.7227E-04, -1.4859E-04, &
      -1.2845E-04, -1.1290E-04, -1.0034E-04, -8.9182E-05, -8.0034E-05, -7.2459E-05, &
      -6.5823E-05, -5.9828E-05, -5.4890E-05, -5.0366E-05, -4.6642E-05, -4.2921E-05, &
      -4.0045E-05, -3.6948E-05, -3.4929E-05, -3.2514E-05, -3.0457E-05, -2.8703E-05, &
      -2.6914E-05, -2.5648E-05, -2.3776E-05, -2.2798E-05, -2.1749E-05, -2.0374E-05, &
      -1.9275E-05, -1.8577E-05, -1.7802E-05, -1.6807E-05, -1.4674E-05, -1.1985E-05, &
      -9.9764E-06, -8.3432E-06, -7.1492E-06, -6.1589E-06, -5.3040E-06, -4.5950E-06, &
      -3.9953E-06, -3.5412E-06, -3.1346E-06, -2.6308E-06, -2.4255E-06, -2.1664E-06, &
      -1.8061E-06 /)
    REAL, PARAMETER :: dsa3_2_CH3COCH3(dim55) = (/ &
      9.1117E-05, 8.7876E-05, 7.8667E-05, 6.6481E-05, 5.4159E-05, 4.3326E-05, &
      3.4607E-05, 2.7929E-05, 2.2909E-05, 1.9017E-05, 1.6106E-05, 1.3845E-05, &
      1.2039E-05, 1.0606E-05, 9.3262E-06, 8.4146E-06, 7.5131E-06, 6.7728E-06, &
      6.1316E-06, 5.6864E-06, 5.1595E-06, 4.7200E-06, 4.4147E-06, 4.0797E-06, &
      3.7396E-06, 3.5361E-06, 3.2833E-06, 3.0966E-06, 2.8868E-06, 2.7184E-06, &
      2.5599E-06, 2.4247E-06, 2.2602E-06, 2.1846E-06, 2.0593E-06, 1.9635E-06, &
      1.8519E-06, 1.7552E-06, 1.6601E-06, 1.5794E-06, 1.4238E-06, 1.1506E-06, &
      9.4632E-07, 8.1659E-07, 6.8318E-07, 5.8587E-07, 5.1928E-07, 4.4674E-07, &
      3.9776E-07, 3.3260E-07, 2.9642E-07, 2.7725E-07, 2.2566E-07, 1.9834E-07, &
      1.8531E-07 /)
    REAL, PARAMETER :: dsb3_1_CH3COCH3(dim55) = (/ &
      3.4148E-01, 3.4114E-01, 3.3742E-01, 3.2926E-01, 3.1752E-01, 3.0418E-01, &
      2.9082E-01, 2.7856E-01, 2.6783E-01, 2.5865E-01, 2.5069E-01, 2.4416E-01, &
      2.3810E-01, 2.3303E-01, 2.2862E-01, 2.2442E-01, 2.2075E-01, 2.1752E-01, &
      2.1453E-01, 2.1168E-01, 2.0920E-01, 2.0682E-01, 2.0477E-01, 2.0262E-01, &
      2.0089E-01, 1.9895E-01, 1.9764E-01, 1.9600E-01, 1.9456E-01, 1.9329E-01, &
      1.9194E-01, 1.9096E-01, 1.8946E-01, 1.8865E-01, 1.8776E-01, 1.8655E-01, &
      1.8556E-01, 1.8491E-01, 1.8417E-01, 1.8320E-01, 1.8106E-01, 1.7803E-01, &
      1.7552E-01, 1.7327E-01, 1.7147E-01, 1.6986E-01, 1.6836E-01, 1.6703E-01, &
      1.6583E-01, 1.6486E-01, 1.6395E-01, 1.6275E-01, 1.6223E-01, 1.6155E-01, &
      1.6056E-01 /)
    REAL, PARAMETER :: dsb3_2_CH3COCH3(dim55) = (/ &
      6.6324E-02, 6.6407E-02, 6.6873E-02, 6.7794E-02, 6.9034E-02, 7.0396E-02, &
      7.1710E-02, 7.2884E-02, 7.3892E-02, 7.4771E-02, 7.5501E-02, 7.6125E-02, &
      7.6668E-02, 7.7136E-02, 7.7585E-02, 7.7928E-02, 7.8289E-02, 7.8605E-02, &
      7.8894E-02, 7.9106E-02, 7.9370E-02, 7.9601E-02, 7.9770E-02, 7.9963E-02, &
      8.0167E-02, 8.0295E-02, 8.0459E-02, 8.0586E-02, 8.0733E-02, 8.0855E-02, &
      8.0974E-02, 8.1079E-02, 8.1211E-02, 8.1274E-02, 8.1381E-02, 8.1465E-02, &
      8.1565E-02, 8.1655E-02, 8.1745E-02, 8.1824E-02, 8.1980E-02, 8.2288E-02, &
      8.2544E-02, 8.2723E-02, 8.2923E-02, 8.3081E-02, 8.3198E-02, 8.3334E-02, &
      8.3433E-02, 8.3571E-02, 8.3653E-02, 8.3698E-02, 8.3828E-02, 8.3899E-02, &
      8.3935E-02 /)
    REAL, PARAMETER :: a4_CH3COCH3(3) = (/ 7.9742E-22, -2.1322E-25, 4.0290E-29 /)
    REAL, PARAMETER :: a5_CH3COCH3(4) = (/ 7.7582E-23, -3.6187E-26, 9.1185E-30, -9.5885E-34 /)

    jval_2d(ip_CH3COCH3)%ptr(:,:) = 0.0_dp

    SELECT CASE(qy_CH3COCH3)
    CASE(1)
    DO k = 1,klev
      DO j = 1,kproma_day
        sig_CH3COCH3(2) = &
          p1(b2_CH3COCH3(i2(j,k)), a2_CH3COCH3(i2(j,k)), &
          v3_du1(j,k))
        sig_CH3COCH3(3) = &
          p1(b3_CH3COCH3(i3(j,k)), a3_CH3COCH3(i3(j,k)),v3_du2(j,k))
        sig_CH3COCH3(4) = &
          p2(a4_CH3COCH3(1),a4_CH3COCH3(2),a4_CH3COCH3(3), v3_du2(j,k))
        sig_CH3COCH3(5) = &
          p3(a5_CH3COCH3(1),a5_CH3COCH3(2),a5_CH3COCH3(3), &
          a5_CH3COCH3(4),v3_du2(j,k))
        ! correction for density dependence
        drel = dens(j,k)/dens_ref
        C1 = p1(DSb3_1_CH3COCH3(i3(j,k)), DSa3_1_CH3COCH3(i3(j,k)), &
          v3_du2(j,k))
        C2 = p1(DSb3_2_CH3COCH3(i3(j,k)), DSa3_2_CH3COCH3(i3(j,k)), &
          v3_du2(j,k))
        F2 = 1./p1(4.3154E-01, 5.6559E-02, DREL)
        F3 = 1./p1(C1,C2,DREL)
        F4 = 1./p1(1.4094E-01, 8.5977E-02, DREL)
        F5 = 1./p1(2.0675E-01, 7.9342E-02, DREL)

        dj=&
          sig_CH3COCH3(2) * F2 * fint(j,k,2)  + &
          sig_CH3COCH3(3) * F3 * fint(j,k,3)  + &
          sig_CH3COCH3(4) * F4 * fint(j,k,4)  + &
          sig_CH3COCH3(5) * F5 * fint(j,k,5)
        jval_2d(ip_CH3COCH3)%ptr(iu0(j),k) = REAL(MAX(0.0, dj*fj_corr(j,6)),dp)
      ENDDO
    ENDDO

    CASE(2)
    DO k = 1,klev
      DO j = 1,kproma_day
        sig_CH3COCH3(2) = &
          p1(b2_CH3COCH3(i2(j,k)), a2_CH3COCH3(i2(j,k)),v3_du1(j,k))
        sig_CH3COCH3(3) = &
          p1(b3_CH3COCH3(i3(j,k)), a3_CH3COCH3(i3(j,k)),v3_du2(j,k))
        sig_CH3COCH3(4) = &
          p2(a4_CH3COCH3(1),a4_CH3COCH3(2),a4_CH3COCH3(3),v3_du2(j,k))
        sig_CH3COCH3(5) = &
          p3(a5_CH3COCH3(1),a5_CH3COCH3(2),a5_CH3COCH3(3), &
          a5_CH3COCH3(4),v3_du2(j,k))
        ! correction for density dependence
        drel = dens(j,k)/dens_ref
        C1 = p1(DSb3_1_CH3COCH3(i3(j,k)), DSa3_1_CH3COCH3(i3(j,k)), &
          v3_du2(j,k))
        C2 = p1(DSb3_2_CH3COCH3(i3(j,k)), DSa3_2_CH3COCH3(i3(j,k)), &
          v3_du2(j,k))
        F2 = 1./p1(4.3154E-01, 5.6559E-02, DREL)
        F3 = 1./p1(C1,C2,DREL)
        F4 = 1./p1(1.4094E-01, 8.5977E-02, DREL)
        F5 = 1./p1(2.0675E-01, 7.9342E-02, DREL)

        ! calculate some correction factors for acetone photolysis,
        ! due to new measurements of the quantum yield
        ! we can skip calculation of 1 & 6 & 7 wavelenght
        !  wavelength(7)
        ! 1 : 2.0513E-05, --> 205
        ! 2 : 2.8777E-05, --> 287 ! used
        ! 3 : 3.0200E-05, --> 302 ! used
        ! 4 : 3.0900E-05, --> 309 ! used
        ! 5 : 3.2000E-05, --> 320 ! used
        ! 6 : 3.7000E-05, --> 370
        ! 7:  5.8000E-05, --> 580
        ! the following lines have been made to increase calculation speed.
        ! If new wavelenght used, check the comments down (real calculation on
        ! all the wavelnght range!)
        ! Quantum Yield of Acetone based on the paper of Gierzack-1998
            A = -15.696 + 0.05707 *WAVE(2) *1.e7
            B = EXP(-88.81 + 0.15161 * WAVE(2)*1.e7 )
            GQYF(2) =  MIN(1.,1./(A + B * DENS(J,K)))

            A = -15.696 + 0.05707 *WAVE(3) *1.e7
            B = EXP(-88.81 + 0.15161 * WAVE(3)*1.e7 )
            GQYF(3) = MIN(1.,1./(A + B * DENS(J,K)))

            A = -130.2 + 0.42884 * WAVE(4)*1.e7
            B = EXP(-55.947 + .044913 * WAVE(4)*1.e7 )
            GQYF(4) =  MIN(1.,1./(A + B * DENS(J,K)))

            A = -130.2 + 0.42884 * WAVE(5)*1.e7
            B = EXP(-55.947 + .044913 * WAVE(5)*1.e7 )
            GQYF(5) =  MIN(1.,1./(A + B * DENS(J,K)))

        ! Quantum yield of Acetone based on the paper of Blitz-2004
            a0 = 0.35*((TEMP(J,K)/295.)**(-1.28))
            b0 =(0.068)*((TEMP(J,K)/295.)**(-2.65))
            D0 = (a0/(1.-a0))*EXP(b0*(WAVE(2)*1.e7-248.))
            PHI_co=(1./(1.+D0))
            a1 = 1.6*1.E-19*((TEMP(J,K)/295.)**(-2.38))
            b1 = 0.55*1.E-3*((TEMP(J,K)/295.)**(-3.19))
            D1 = a1*EXP(-b1*((1.E7/(WAVE(2)*1.e7))-33113.))
            PHI_ch3co = (1.-PHI_co)/(1.+DENS(J,K)*D1)
            BQYF(2)=MIN(1.E+0,PHI_co+PHI_ch3co)

            a0 = 0.35*((TEMP(J,K)/295.)**(-1.28))
            b0 =(0.068)*((TEMP(J,K)/295.)**(-2.65))
            D0 = (a0/(1.-a0))*EXP(b0*(WAVE(3)*1.e7-248.))
            PHI_co=(1./(1.+D0))
            a1 = 1.6*1.E-19*((TEMP(J,K)/295.)**(-2.38))
            b1 = 0.55*1.E-3*((TEMP(J,K)/295.)**(-3.19))
            D1 = a1*EXP(-b1*((1.E7/(WAVE(3)*1.e7))-33113.))
            PHI_ch3co = (1.-PHI_co)/(1.+DENS(J,K)*D1)
            BQYF(3)=MIN(1.E+0,PHI_co+PHI_ch3co)

            a0 = 0.35*((TEMP(J,K)/295.)**(-1.28))
            b0 =(0.068)*((TEMP(J,K)/295.)**(-2.65))
            D0 = (a0/(1.-a0))*EXP(b0*(WAVE(4)*1.e7-248.))
            PHI_co=(1./(1.+D0))
            a2 = 1.62*1.E-17*((TEMP(J,K)/295.)**(-10.03))
            b2 = 1.79*1.E-03*((TEMP(J,K)/295.)**(-1.364))
            D2 = a2 *EXP(-b2*((1.E7/(WAVE(4)*1.e7))-30488.))
            a3 = 26.29 * ((TEMP(J,K)/295.)**(-6.59))
            b3 = 5.72*1.E-7 *((TEMP(J,K)/295.)**(-2.93))
            E3 = 30006. * ((TEMP(J,K)/295.)**(-0.064))
            D3 = a3*EXP(-b3*(((1.E7/(WAVE(4)*1.e7))-E3)**2))
            A4 = 1.67*1.E-15*((TEMP(J,K)/295.)**(-7.25))
            B4 = 2.08*1.E-05*((TEMP(J,K)/295.)**(-1.16))
            D4 = A4 *EXP(-B4*((1.E7/(WAVE(4)*1.e7))-30488.))
            PHI_ch3co = (1.-PHI_co) * (1.+D4*DENS(J,K)+D3) &
              / ((1.+D2*DENS(J,K)+D3)*(1.+D4*DENS(J,K)))
            BQYF(4)=MIN(1.E+0,PHI_co+PHI_ch3co)

            a0 = 0.35*((TEMP(J,K)/295.)**(-1.28))
            b0 =(0.068)*((TEMP(J,K)/295.)**(-2.65))
            D0 = (a0/(1.-a0))*EXP(b0*(WAVE(5)*1.e7-248.))
            PHI_co=(1./(1.+D0))
            a2 = 1.62*1.E-17*((TEMP(J,K)/295.)**(-10.03))
            b2 = 1.79*1.E-03*((TEMP(J,K)/295.)**(-1.364))
            D2 = a2 *EXP(-b2*((1.E7/(WAVE(5)*1.e7))-30488.))
            a3 = 26.29 * ((TEMP(J,K)/295.)**(-6.59))
            b3 = 5.72*1.E-7 *((TEMP(J,K)/295.)**(-2.93))
            E3 = 30006. * ((TEMP(J,K)/295.)**(-0.064))
            D3 = a3*EXP(-b3*(((1.E7/(WAVE(5)*1.e7))-E3)**2))
            A4 = 1.67*1.E-15*((TEMP(J,K)/295.)**(-7.25))
            B4 = 2.08*1.E-05*((TEMP(J,K)/295.)**(-1.16))
            D4 = A4 *EXP(-B4*((1.E7/(WAVE(5)*1.e7))-30488.))
            PHI_ch3co = (1.-PHI_co) * (1.+D4*DENS(J,K)+D3) &
              / ((1.+D2*DENS(J,K)+D3)*(1.+D4*DENS(J,K)))
            BQYF(5)=MIN(1.E+0,PHI_co+PHI_ch3co)

        ! quantum yield correction
          nqyf(2) = bqyf(2)/gqyf(2) ! BLITZ 2004
          nqyf(3) = bqyf(3)/gqyf(3)
          nqyf(4) = bqyf(4)/gqyf(4)
          nqyf(5) = bqyf(5)/gqyf(5)

        ! factors that correct the Jval depending on the "new" quantum yield
        dj=&
          sig_CH3COCH3(2) * F2 * fint(j,k,2) * nqyf(2) + &
          sig_CH3COCH3(3) * F3 * fint(j,k,3) * nqyf(3) + &
          sig_CH3COCH3(4) * F4 * fint(j,k,4) * nqyf(4) + &
          sig_CH3COCH3(5) * F5 * fint(j,k,5) * nqyf(5)
        jval_2d(ip_CH3COCH3)%ptr(iu0(j),k) = REAL(MAX(0.0, dj*fj_corr(j,6)),dp)
      ENDDO
    ENDDO

    CASE(3)

    DO k = 1,klev
      DO j = 1,kproma_day
        sig_CH3COCH3(2) = &
          p1(b2_CH3COCH3(i2(j,k)), a2_CH3COCH3(i2(j,k)),v3_du1(j,k))
        sig_CH3COCH3(3) = &
          p1(b3_CH3COCH3(i3(j,k)), a3_CH3COCH3(i3(j,k)),v3_du2(j,k))
        sig_CH3COCH3(4) = &
          p2(a4_CH3COCH3(1),a4_CH3COCH3(2),a4_CH3COCH3(3), v3_du2(j,k))
        sig_CH3COCH3(5) = &
          p3(a5_CH3COCH3(1),a5_CH3COCH3(2),a5_CH3COCH3(3), &
          a5_CH3COCH3(4),v3_du2(j,k))
        ! correction for density dependence
        drel = dens(j,k)/dens_ref
        C1 = p1(DSb3_1_CH3COCH3(i3(j,k)), DSa3_1_CH3COCH3(i3(j,k)), &
          v3_du2(j,k))
        C2 = p1(DSb3_2_CH3COCH3(i3(j,k)), DSa3_2_CH3COCH3(i3(j,k)), &
          v3_du2(j,k))
        F2 = 1./p1(4.3154E-01, 5.6559E-02, DREL)
        F3 = 1./p1(C1,C2,DREL)
        F4 = 1./p1(1.4094E-01, 8.5977E-02, DREL)
        F5 = 1./p1(2.0675E-01, 7.9342E-02, DREL)

        ! calculate some correction factors for acetone photolysis,
        ! due to new measurements of the quantum yield
        ! we can skip calculation of 1 & 6 & 7 wavelenght
        !  wavelength(7)
        ! 1 : 2.0513E-05, --> 205
        ! 2 : 2.8777E-05, --> 287 ! used
        ! 3 : 3.0200E-05, --> 302 ! used
        ! 4 : 3.0900E-05, --> 309 ! used
        ! 5 : 3.2000E-05, --> 320 ! used
        ! 6 : 3.7000E-05, --> 370
        ! 7:  5.8000E-05, --> 580
        ! the following lines have been made to increase calculation speed.
        ! If new wavelenght used, check the comments down (real calculation on
        ! all the wavelnght range!)
        ! Quantum Yield of Acetone based on the paper of Gierzack-1998
            A = -15.696 + 0.05707 *WAVE(2) *1.e7
            B = EXP(-88.81 + 0.15161 * WAVE(2)*1.e7 )
            GQYF(2) =  MIN(1.,1./(A + B * DENS(J,K)))

            A = -15.696 + 0.05707 *WAVE(3) *1.e7
            B = EXP(-88.81 + 0.15161 * WAVE(3)*1.e7 )
            GQYF(3) = MIN(1.,1./(A + B * DENS(J,K)))

            A = -130.2 + 0.42884 * WAVE(4)*1.e7
            B = EXP(-55.947 + .044913 * WAVE(4)*1.e7 )
            GQYF(4) =  MIN(1.,1./(A + B * DENS(J,K)))

            A = -130.2 + 0.42884 * WAVE(5)*1.e7
            B = EXP(-55.947 + .044913 * WAVE(5)*1.e7 )
            GQYF(5) =  MIN(1.,1./(A + B * DENS(J,K)))

        ! quantum yield of acetone based on Warneck
        ! Atmos. Environ. 35, 5773-5777 (2001)

          PhiS=(1.-0.113)*((1.+EXP(((WAVE(2)*1.E7)-307.5)/3.))**(-1.))+0.113
          WQYF(2) = &
            ((1./PhiS)+7.138E-7*DENS(J,K)*EXP((-8780.6/(WAVE(2)*1.E7))))**(-1.)

          PhiS=(1.-0.113)*((1.+EXP(((WAVE(3)*1.E7)-307.5)/3.))**(-1.))+0.113
          WQYF(3) = &
            ((1./PhiS)+7.138E-7*DENS(J,K)*EXP((-8780.6/(WAVE(3)*1.E7))))**(-1.)

          PhiS=(1.-0.113)*((1.+EXP(((WAVE(4)*1.E7)-307.5)/3.))**(-1.))+0.113
          WQYF(4) = &
            ((1./PhiS)+7.138E-7*DENS(J,K)*EXP((-8780.6/(WAVE(4)*1.E7))))**(-1.)

          PhiS=(1.-0.113)*((1.+EXP(((WAVE(5)*1.E7)-307.5)/3.))**(-1.))+0.113
          WQYF(5) = &
            ((1./PhiS)+7.138E-7*DENS(J,K)*EXP((-8780.6/(WAVE(5)*1.E7))))**(-1.)

        ! quantum yield correction
          nqyf(2) = wqyf(2)/gqyf(2) ! IUPAC
          nqyf(3) = wqyf(3)/gqyf(3)
          nqyf(4) = wqyf(4)/gqyf(4)
          nqyf(5) = wqyf(5)/gqyf(5)

        ! factors that correct the Jval depending on the "new" quantum yield
        dj=&
          sig_CH3COCH3(2) * F2 * fint(j,k,2) * nqyf(2) + &
          sig_CH3COCH3(3) * F3 * fint(j,k,3) * nqyf(3) + &
          sig_CH3COCH3(4) * F4 * fint(j,k,4) * nqyf(4) + &
          sig_CH3COCH3(5) * F5 * fint(j,k,5) * nqyf(5)
        jval_2d(ip_CH3COCH3)%ptr(iu0(j),k) = REAL(MAX(0.0, dj*fj_corr(j,6)),dp)
      ENDDO
    ENDDO

    END SELECT

  END SUBROUTINE jval_cal_CH3COCH3

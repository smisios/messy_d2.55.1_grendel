  SUBROUTINE jval_cal_MACR

    ! MACR -> products JPL 2006, total quantum yield fitted; first channel to CH2=CCH3 + CHO most likely, C3H6 not seen.

    INTEGER :: j, k
    REAL    :: dj
    REAL, DIMENSION(0:MAXWAV) :: sig_MACR

    REAL :: f6

    ! T-const parameters:
    REAL, PARAMETER :: a3_MACR(dim55) = (/ &
      2.2318E-25, 1.4224E-25, 8.4374E-26, 5.3568E-26, 3.7145E-26, 2.7383E-26, &
      2.1222E-26, 1.6775E-26, 1.3863E-26, 1.1467E-26, 9.7569E-27, 8.3861E-27, &
      7.3592E-27, 6.3334E-27, 5.6477E-27, 5.1353E-27, 4.4497E-27, 4.1061E-27, &
      3.7671E-27, 3.4247E-27, 3.2478E-27, 2.7398E-27, 2.7396E-27, 2.3973E-27, &
      2.3932E-27, 2.0547E-27, 2.0549E-27, 1.8835E-27, 1.7094E-27, 1.7124E-27, &
      1.3698E-27, 1.5410E-27, 1.3677E-27, 1.1986E-27, 1.1986E-27, 1.1966E-27, &
      1.0272E-27, 1.0275E-27, 8.5607E-28, 1.0257E-27, 7.1894E-28, 5.8200E-28, &
      4.1066E-28, 3.0813E-28, 2.3965E-28, 2.0531E-28, 1.3693E-28, 1.0271E-28, &
      6.8465E-29, 6.8445E-29, 3.4266E-29, 3.4211E-29, 3.4223E-29, 3.4223E-29, &
      -0.0000E+00 /)
    REAL, PARAMETER :: b3_MACR(dim55) = (/ &
      1.5812E-22, 1.6289E-22, 1.6968E-22, 1.7509E-22, 1.7894E-22, 1.8180E-22, &
      1.8396E-22, 1.8578E-22, 1.8714E-22, 1.8840E-22, 1.8940E-22, 1.9029E-22, &
      1.9101E-22, 1.9179E-22, 1.9235E-22, 1.9280E-22, 1.9344E-22, 1.9378E-22, &
      1.9414E-22, 1.9452E-22, 1.9472E-22, 1.9535E-22, 1.9535E-22, 1.9581E-22, &
      1.9581E-22, 1.9631E-22, 1.9631E-22, 1.9658E-22, 1.9686E-22, 1.9686E-22, &
      1.9746E-22, 1.9715E-22, 1.9747E-22, 1.9780E-22, 1.9780E-22, 1.9780E-22, &
      1.9816E-22, 1.9816E-22, 1.9854E-22, 1.9815E-22, 1.9887E-22, 1.9923E-22, &
      1.9973E-22, 2.0006E-22, 2.0030E-22, 2.0043E-22, 2.0071E-22, 2.0086E-22, &
      2.0102E-22, 2.0102E-22, 2.0120E-22, 2.0120E-22, 2.0120E-22, 2.0120E-22, &
      2.0142E-22 /)
    REAL, PARAMETER :: a4_MACR(2)= (/2.3481E-22,  1.9076E-27/)
    REAL, PARAMETER :: a5_MACR(2)= (/2.5435E-22, -1.2548E-27/)
    REAL, PARAMETER :: a6_MACR(2)= (/3.8102E-23, -2.5458E-28/)
    jval_2d(ip_MACR)%ptr(:,:) = 0.0_dp
    DO k = 1,klev
      DO j = 1,kproma_day
        sig_MACR(3) = &
          p1(b3_MACR(i3(j,k)), a3_MACR(i3(j,k)), &
          v3_du2(j,k))
        sig_MACR(4) = &
          p1(a4_MACR(1),a4_MACR(2),v3_du2(j,k))
        sig_MACR(5) = &
          p1(a5_MACR(1),a5_MACR(2),v3_du2(j,k))
        sig_MACR(6) = &
          p1(a6_MACR(1),a6_MACR(2),v3_du2(j,k))
        f6=1025./0.5856/(1000.+0.0075*press(j,k))
        dj = 0. &
          + sig_MACR(3) *F6 * fint(j,k,3) &
          + sig_MACR(4) *F6 * fint(j,k,4) &
          + sig_MACR(5) *F6 * fint(j,k,5) &
          + sig_MACR(6) *F6 * fint(j,k,6)

        jval_2d(ip_MACR)%ptr(iu0(j),k) = REAL(MAX(0.0, dj*fj_corr(j,6)),dp)
      ENDDO
    ENDDO

  END SUBROUTINE jval_cal_MACR

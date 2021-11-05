  SUBROUTINE jval_cal_GLYOX

    ! CHOCHO -> 2 CHO   (CHOCHO = glyox = glyoxal) Tadic et al, 2006, JPPA 177, only total yield

    INTEGER :: j, k
    REAL    :: dj
    REAL, DIMENSION(0:MAXWAV) :: sig_GLYOX

    REAL :: f6
    REAL :: prel

    ! T-const parameters:
    REAL, PARAMETER :: a3_GLYOX(dim55) = (/ &
      -2.1052E-24, -4.6556E-24, -5.3568E-24, -5.1857E-24, -4.7073E-24, -4.1759E-24, &
      -3.6967E-24, -3.2523E-24, -2.8752E-24, -2.5501E-24, -2.2766E-24, -2.0366E-24, &
      -1.8313E-24, -1.6604E-24, -1.4889E-24, -1.3694E-24, -1.2493E-24, -1.1463E-24, &
      -1.0616E-24, -9.7599E-25, -9.0600E-25, -8.3905E-25, -7.8765E-25, -7.3631E-25, &
      -6.8374E-25, -6.3356E-25, -5.9932E-25, -5.6507E-25, -5.1285E-25, -4.9657E-25, &
      -4.6232E-25, -4.4523E-25, -4.1022E-25, -3.9383E-25, -3.5963E-25, -3.4184E-25, &
      -3.2539E-25, -3.0818E-25, -2.9114E-25, -2.7346E-25, -2.2938E-25, -1.7117E-25, &
      -1.3004E-25, -9.5862E-26, -7.5320E-26, -5.4749E-26, -4.4513E-26, -3.0807E-26, &
      -2.3963E-26, -2.0531E-26, -1.3695E-26, -1.0273E-26, -1.0265E-26, -6.8473E-27, &
      -3.4280E-27 /)
    REAL, PARAMETER :: b3_GLYOX(dim55) = (/ &
      3.5916E-20, 3.6066E-20, 3.6149E-20, 3.6119E-20, 3.6007E-20, 3.5851E-20, &
      3.5683E-20, 3.5501E-20, 3.5324E-20, 3.5153E-20, 3.4993E-20, 3.4839E-20, &
      3.4695E-20, 3.4565E-20, 3.4425E-20, 3.4320E-20, 3.4208E-20, 3.4105E-20, &
      3.4016E-20, 3.3921E-20, 3.3839E-20, 3.3757E-20, 3.3691E-20, 3.3622E-20, &
      3.3548E-20, 3.3475E-20, 3.3423E-20, 3.3369E-20, 3.3283E-20, 3.3256E-20, &
      3.3196E-20, 3.3165E-20, 3.3099E-20, 3.3068E-20, 3.3000E-20, 3.2963E-20, &
      3.2929E-20, 3.2891E-20, 3.2854E-20, 3.2813E-20, 3.2710E-20, 3.2557E-20, &
      3.2437E-20, 3.2327E-20, 3.2255E-20, 3.2177E-20, 3.2135E-20, 3.2075E-20, &
      3.2043E-20, 3.2026E-20, 3.1990E-20, 3.1971E-20, 3.1971E-20, 3.1950E-20, &
      3.1928E-20 /)
    REAL, PARAMETER :: a4_GLYOX(2) = (/3.0906E-20, 1.9869E-25/)
    REAL, PARAMETER :: a5_GLYOX(3) = (/ 1.5819E-20,-1.1829E-24, 8.9498E-29/)
    REAL, PARAMETER :: a6_GLYOX(2) = (/1.7007E-21,-6.2818E-27/)
    REAL, PARAMETER :: a7_GLYOX(2) = (/1.9405E-23, 1.3357E-27/)
    jval_2d(ip_GLYOX)%ptr(:,:) = 0.0_dp
    DO k = 1,klev
      DO j  = 1,kproma_day
        PREL = (MIN(1.E5,PRESS(J,K))/PRESS_REF)  !relative pressure
        sig_GLYOX(3) = &
          p1(b3_GLYOX(i3(j,k)), a3_GLYOX(i3(j,k)), &
          v3_du2(j,k))
        sig_GLYOX(4) = &
          p1(a4_GLYOX(1),a4_GLYOX(2),v3_du2(j,k))
        sig_GLYOX(5) = &
          p2(a5_GLYOX(1),a5_GLYOX(2),a5_GLYOX(3),v3_du2(j,k))
        sig_GLYOX(6) = &
          p1(a6_GLYOX(1),a6_GLYOX(2),v3_du2(j,k))
        sig_GLYOX(7) = &
          p1(a7_GLYOX(1),a7_GLYOX(2),v3_du2(j,k))
        f6=25.7/(6.8+ 1.89*prel)
        dj = 0. &
          + sig_GLYOX(3) *     fint(j,k,3) &
          + sig_GLYOX(4) *     fint(j,k,4) &
          + sig_GLYOX(5) *     fint(j,k,5) &
          + sig_GLYOX(6) * F6* fint(j,k,6) &
          + sig_GLYOX(7) * F6* fint(j,k,7)
        jval_2d(ip_GLYOX)%ptr(iu0(j),k) = REAL(MAX(0.0, dj*fj_corr(j,6)),dp)
      ENDDO
    ENDDO

  END SUBROUTINE jval_cal_GLYOX

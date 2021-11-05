  SUBROUTINE jval_cal_H2SO4

    ! vibrational overtones: Vaida et al 2003 Science 299, Hintze et al 2003,
    !  J.Ph.Ch. A107; Lyalpha: Mills et al, 2005, JGR 110, D24208

    ! mz_rs_20111006+
    ! \cite{2649,2651}, Lyman-alpha from \cite{2650}. However, it is
    ! unclear which spectrum was used to derive the parameters
    ! \verb|a7_H2SO4(:)|.
    ! mz_rs_20111006-

    INTEGER :: j, k
    REAL    :: dj
    REAL, DIMENSION(0:MAXWAV) :: sig_H2SO4

    ! T-const parameters:
    REAL, PARAMETER :: a7_H2SO4(2) = (/4.2204E-25, 6.9058E-30/)

    jval_2d(ip_H2SO4)%ptr(:,:) = 0.0_dp
    DO k = 1,klev
      DO j = 1,kproma_day
        sig_H2SO4(7) = p1(a7_H2SO4(1),a7_H2SO4(2),v3_du2(j,k))
        dj = 0. &
          + 6.4e-18 * r_m(j,k) * phi_la &
          + sig_H2SO4(7) * fint(j,k,7)
        jval_2d(ip_H2SO4)%ptr(iu0(j),k) = REAL(MAX(0.0, dj*fj_corr(j,1)),dp)
      ENDDO
    ENDDO
    
  END SUBROUTINE jval_cal_H2SO4

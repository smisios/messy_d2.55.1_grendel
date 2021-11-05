      SUBROUTINE RHO2(T,S,P,RH)

      USE mo_kind, ONLY: wp
      IMPLICIT NONE

!********************************************************
!SJ***MODIFIED BY SJKIM
! ZUSTANDSGLEICHUNG
! UNTERPROGRAMM NACH ADRIAN GILL (ANHANG)
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
      REAL(wp) S,T,P,RH,S3H
      REAL(wp) :: &
           b0, b1, b2, b3, b4, &
           c0, c1, c2, &
           d0, &
           a0, a1, a2, a3, a4, a5, &
           f0, f1, f2, f3, &
           g0, g1, g2, &
           ai0, ai1, ai2, &
           aj0, &
           am0, am1, am2, &
           e0, e1, e2, e3, e4, &
           h0, h1, h2, h3, &
           ak0, ak1, ak2

      ! FIXME: turn this into PARAMETERs
      DATA B0,B1,B2,B3,B4/ 8.24493E-1_wp, -4.0899E-3_wp, 7.6438E-5_wp, &
     &-8.2467E-7_wp, 5.3875E-9_wp/
      DATA C0,C1,C2/-5.72466E-3_wp, 1.0227E-4_wp, -1.6546e-6_wp/
      DATA D0 /4.8314E-4_wp/
      DATA A0,A1,A2,A3,A4,A5 /999.842594_wp, 6.793952E-2_wp,            &
     &-9.095290e-3_wp, 1.001685e-4_wp, -1.120083E-6_wp, 6.536332E-9_wp/
      DATA F0,F1,F2,F3 / 54.6746_wp, -0.603459_wp,                      &
     &1.09987e-2_wp, -6.1670e-5_wp/
      DATA G0,G1,G2 /7.944E-2_wp, 1.6483E-2_wp, -5.3009E-4_wp/
      DATA AI0,AI1,AI2 /2.2838e-3_wp, -1.0981e-5_wp, -1.6078e-6_wp/
      DATA AJ0 /1.91075e-4_wp/
      DATA AM0,AM1,AM2 /-9.9348E-7_wp, 2.0816E-8_wp, 9.1697E-10_wp/
      DATA E0,E1,E2,E3,E4 /19652.21_wp, 148.4206_wp, -2.327105_wp,      &
     &1.360477E-2_wp, -5.155288E-5_wp/
      DATA H0,H1,H2,H3 /3.239908_wp, 1.43713E-3_wp,                     &
     &1.16092E-4_wp, -5.77905E-7_wp/
      DATA AK0,AK1,AK2 /8.50935E-5_wp, -6.12293E-6_wp, 5.2787E-8_wp/
!
      S3H=SQRT(S**3)
!
      RH = 0.001_wp * (A0+T*(A1+T                                       &
     &       *(A2+T*(A3+T*(A4+T*A5))))                                  &
     &       +S*(B0+T*(B1+T                                             &
     &      *(B2+T*(B3+T*B4))))+D0*S**2                                 &
     &+S3H*(C0+T*(C1+C2*T)) )                                           &
     &       / (1._wp - p/(p*(                                          &
     &   H0+T*(H1+T*(H2+T*H3))                                          &
     &  +S*(AI0+T*(AI1+AI2*T))+AJ0*S3H                                  &
     &  +(AK0+T*(AK1+T*AK2)                                             &
     &  +S*(AM0+T*(AM1+T*AM2)))*P)+                                     &
     &    E0+T*(E1+T*(E2+T*(E3+T*E4)))                                  &
     &      +S*(F0+T*(F1+T*(F2+T*F3)))                                  &
     &      +S3H*(G0+T*(G1+G2*T))))
      RETURN
      END

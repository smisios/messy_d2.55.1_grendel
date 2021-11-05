!********************************************************
! Equation of state following Adrian Gill (Appendix)
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
SUBROUTINE rho1(t, s, p, rh)

  USE mo_kind, ONLY: dp, wp
  USE mo_param1

  IMPLICIT NONE

  REAL(dp), INTENT(INOUT) :: s
  REAL(dp), INTENT(IN) :: t, p
  REAL(dp), INTENT(OUT) :: rh


  REAL(dp), PARAMETER :: a0 =  999.8426_dp,     &
                         a1 =  6.793952e-2_dp,  &
                         a2 = -9.095290e-3_dp,  &
                         a3 =  1.001685e-4_dp,  &
                         a4 = -1.120083e-6_dp,  &
                         a5 =  6.536332e-9_dp
  REAL(dp), PARAMETER :: b0 =  8.24493e-1_dp,   &
                         b1 = -4.0899e-3_dp,    &
                         b2 =  7.6438e-5_dp,    &
                         b3 = -8.2467e-7_dp,    &
                         b4 =  5.3875e-9_dp
  REAL(dp), PARAMETER :: c0 = -5.72466e-3_dp,   &
                         c1 =  1.0227e-4_dp,    &
                         c2 = -1.6546e-6_dp
  REAL(dp), PARAMETER :: d0 =  4.8314e-4_dp
  REAL(dp), PARAMETER :: e0 =  19652.21_dp,     &
                         e1 =  148.4206_dp,     &
                         e2 = -2.327105_dp,     &
                         e3 =  1.360477e-2_dp,  &
                         e4 = -5.155288e-5_dp
  REAL(dp), PARAMETER :: f0 =  54.6746_dp,      &
                         f1 = -0.603459_dp,     &
                         f2 =  1.09987e-2_dp,   &
                         f3 = -6.1670e-5_dp
  REAL(dp), PARAMETER :: g0 =  7.944e-2_dp,     &
                         g1 =  1.6483e-2_dp,    &
                         g2 = -5.3009e-4_dp
  REAL(dp), PARAMETER :: h0 =  3.239908_dp,     &
                         h1 =  1.43713e-3_dp,   &
                         h2 =  1.16092e-4_dp,   &
                         h3 = -5.77905e-7_dp
  REAL(dp), PARAMETER :: ai0 =  2.2838e-3_dp,   &
                         ai1 = -1.0981e-5_dp,   &
                         ai2 = -1.6078e-6_dp
  REAL(dp), PARAMETER :: aj0 =  1.91075e-4_dp
  REAL(dp), PARAMETER :: ak0 =  8.50935e-5_dp,  &
                         ak1 = -6.12293e-6_dp,  &
                         ak2 =  5.2787e-8_dp
  REAL(dp), PARAMETER :: am0 = -9.9348e-7_dp,   &
                         am1 =  2.0816e-8_dp,   &
                         am2 =  9.1697e-10_dp


  REAL(dp) :: s3h


  s=MAX(s,28.0_dp)
  s3h=SQRT(s**3)

  rh=(a0+t*(a1+t                                          &
              *(a2+t*(a3+t*(a4+t*a5))))                   &
              +s*(b0+t*(b1+t                              &
              *(b2+t*(b3+t*b4))))+d0*s**2                 &
              +s3h*(c0+t*(c1+c2*t)) )                     &
              /(1._dp-p/(p*(                                 &
              h0+t*(h1+t*(h2+t*h3))                       &
              +s*(ai0+t*(ai1+ai2*t))+aj0*s3h              &
              +(ak0+t*(ak1+t*ak2)                         &
              +s*(am0+t*(am1+t*am2)))*p)+                 &
              e0+t*(e1+t*(e2+t*(e3+t*e4)))                &
              +s*(f0+t*(f1+t*(f2+t*f3)))                  &
              +s3h*(g0+t*(g1+g2*t))))


END SUBROUTINE rho1







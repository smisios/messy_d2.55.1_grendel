!>
!! Defines coefficients for use in SBR ADISIT.
!!
!! @author O. Boehringer, *DKRZ*
!! @date 19.11.95
!!
MODULE MO_PARAM3

  USE mo_kind, ONLY: wp
  IMPLICIT NONE

  REAL(wp), PARAMETER :: &
       A1=3.6504E-4_wp, A2=8.3198E-5_wp, A3=5.4065E-7_wp, A4=4.0274E-9_wp, &
       B1=1.7439E-5_wp, B2=2.9778E-7_wp, &
       C1=8.9309E-7_wp, C2=3.1628E-8_wp, C3=2.1987E-10_wp, &
       D=4.1057E-9_wp, &
       E1=1.6056E-10_wp, E2=5.0484E-12_wp

END MODULE MO_PARAM3

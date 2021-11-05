      MODULE MO_PARAM3
!-----------------------------------------------------------------------
!
!*   PARAMETER *PARAM3*  - coefficients for use in SBR ADISIT.
!
!    Author.
!    -------
!    O. Boehringer,        *DKRZ*         19.11.95
!
!-----------------------------------------------------------------------

      ! mz_rs_20090909+
      !!$In file mo_param3.f90:11
      !!$
      !!$      PARAMETER(A1=3.6504E-4,A2=8.3198E-5,A3=5.4065E-7,A4=4.0274E-9, &
      !!$                1
      !!$Error: Symbol 'a1' at (1) has no IMPLICIT type
      !!$gmake[1]: *** [mo_param3.o] Error 1
      ! mz_rs_20090909-

      PARAMETER(A1=3.6504E-4,A2=8.3198E-5,A3=5.4065E-7,A4=4.0274E-9,    &
     & B1=1.7439E-5,B2=2.9778E-7,C1=8.9309E-7,C2=3.1628E-8,             &
     & C3=2.1987E-10,D=4.1057E-9,E1=1.6056E-10,E2=5.0484E-12)
      END MODULE MO_PARAM3

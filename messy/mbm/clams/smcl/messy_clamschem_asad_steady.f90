Module messy_clamschem_asad_steady

contains

! *****************************COPYRIGHT*******************************
!
! (c) [University of Cambridge] [2008]. All rights reserved.
! This routine has been licensed to the Met Office for use and
! distribution under the UKCA collaboration agreement, subject 
! to the terms and conditions set out therein.
! [Met Office Ref SC138] 
!
! *****************************COPYRIGHT*******************************
!
!  Description:
!    Computes steady-state species for Newton-Raphson integrator.
!    Part of the ASAD chemical solver.
!
!     ASAD: ycn                      Version: steady.f 4.1 20/7/07
!
!     Purpose
!     -------
!
!  Routine to explicitly define steady state expressions - prevents
!  generality, but removes need for sluggish iteration round family
!  stuff.
!
!  Important Notes:
!     1) Ordering of calculations is important - need to avoid feedbacks!
!     2) Needs to be rewritten whenever reaction numbering is changed.
!
!  Additions:
!     To improve generality, reactions involving steady state species
!  are selected in 'setsteady' and loaded into nss* integer arrays,
!  which are then used in this routine. This causes a very slight
!  increase in CPU time, but removes the need to rewrite the routine
!  whenever a new species is added or a reaction changed.
!
!                                            Oliver   (3 Feb 1998)
! We add more general terms for the steady-state species. It is assumed
! that O(1D), O(3P), H and N can be put in steady state.
!
! 08/2015 Jens-Uwe Grooss (jug): 
!     added the if case for o1d_in_ss to allow schemes in which O(1D)
!     is not handled as a steady state species
!
!
!     Method
!     ------
!
!   We sum up production and loss terms for SS species, and divide.
!   Moreover, corresponding terms in the Jacobian are calculated that
!   account for the dependence of steady-state variables on tracer
!   variables.
!
!
!  Code Description:
!    Language:  FORTRAN 90
!
! ######################################################################
!
      SUBROUTINE asad_steady( kl )

      USE messy_clamschem_asad_mod, ONLY: deriv, y, rk, peps,            &
                                          nspi, nssi, nssrt, nssrx,      &
                                          nssri, nsspt, nsspi, nsst,     &
                                          jlst, nspo1d, nspo3, nspoh,    &
                                          nspo3p, nsph, nuni,            &
                                          nspho2, nspno, nspn, nss_o3p,  &
                                          nss_o1d, nss_n, nss_h,         &
                                          o1d_in_ss, o3p_in_ss, n_in_ss, h_in_ss
      USE messy_clams_global,             ONLY: prec
      USE messy_clamschem_asad_mod_clams, ONLY: theta_field_size
      IMPLICIT NONE


! Subroutine interface
      INTEGER, INTENT(IN) :: kl            ! No. of points

! Local variables
      INTEGER, PARAMETER :: n_o3=1           ! indicies for dssden etc
      INTEGER, PARAMETER :: n_oh=2
      INTEGER, PARAMETER :: n_ho2=3
      INTEGER, PARAMETER :: n_no=4

      INTEGER :: jl
      INTEGER :: jl1
      INTEGER :: jl2
      INTEGER :: jr
      INTEGER :: ix
      INTEGER :: i
      INTEGER :: j

      REAL(PREC) :: ssnum(theta_field_size)
      REAL(PREC) :: ssden(theta_field_size)

! Add here derivatives w.r.t. O3, OH, and HO2, and NO of numerator and
! denominator of steady state species
      REAL(PREC) :: dssnum(theta_field_size,4)
      REAL(PREC) :: dssden(theta_field_size,4)

!
! Set up loops correctly
      jl1 = 1
      jl2 = kl
      IF (kl == 1) THEN
        jl1 = jl1+jlst
        jl2 = jl1
      END IF

! Initialise DERIV - now done in ASAD_INIT, to 1.0
!     deriv(:,:,:) = 0.0

! Loop through steady state species

      DO ix = 1,nsst
        ssnum = 0.
        ssden = 0.
        dssden = 0.
        dssnum = 0.

! Production terms
        DO jr = 1,nsspt(ix)
          i = nsspi(ix,jr)
          IF (i <= nuni) THEN
            ssnum(jl1:jl2) = ssnum(jl1:jl2) +                           &
                rk(jl1:jl2,i)*y(jl1:jl2,nspi(i,1))

            IF ((ix < 5) .AND. (nspi(i,1) == nspo3 ))                   &
! add terms to derivative for d(j[O3])/d[O3] = j_o3
              dssnum(jl1:jl2,n_o3) = dssnum(jl1:jl2,n_o3) +             &
                                       rk(jl1:jl2,i)
            IF ((ix < 5) .AND. (nspi(i,1) == nspno ))                   &
! add terms to derivative for d(j[NO])/d[NO] = j_no
              dssnum(jl1:jl2,n_no) = dssnum(jl1:jl2,n_no) +             &
                                       rk(jl1:jl2,i)
          ELSE
            ssnum(jl1:jl2) = ssnum(jl1:jl2) +                           &
                rk(jl1:jl2,i)*y(jl1:jl2,nspi(i,1))*y(jl1:jl2,nspi(i,2))
            IF (ix < 5) THEN

! add terms for derivative w.r.t. ozone.
              IF (nspi(i,1) == nspo1d)                                  &
                dssnum(jl1:jl2,n_o3) = dssnum(jl1:jl2,n_o3) +           &
                     rk(jl1:jl2,i)                                      &
                    *y(jl1:jl2,nspi(i,2))*deriv(jl1:jl2,nss_o1d,n_o3)

              IF (nspi(i,2) == nspo1d)                                  &
                  dssnum(jl1:jl2,n_o3) = dssnum(jl1:jl2,n_o3) +         &
                    rk(jl1:jl2,i)                                       &
                    *y(jl1:jl2,nspi(i,1))*deriv(jl1:jl2,nss_o1d,n_o3)

              IF (nspi(i,1) == nspo3p)                                  &
                  dssnum(jl1:jl2,n_o3) = dssnum(jl1:jl2,n_o3) +         &
                     rk(jl1:jl2,i)                                      &
                    *y(jl1:jl2,nspi(i,2))*deriv(jl1:jl2,nss_o3p,n_o3)

              IF (nspi(i,2) == nspo3p)                                  &
                  dssnum(jl1:jl2,n_o3) = dssnum(jl1:jl2,n_o3) +         &
                    rk(jl1:jl2,i)                                       &
                    *y(jl1:jl2,nspi(i,1))*deriv(jl1:jl2,nss_o3p,n_o3)

              IF (nspi(i,1) == nspo3)                                   &
                dssnum(jl1:jl2,n_o3) = dssnum(jl1:jl2,n_o3) +           &
                    rk(jl1:jl2,i)                                       &
                    *y(jl1:jl2,nspi(i,2))

              IF (nspi(i,2) == nspo3)                                   &
                dssnum(jl1:jl2,n_o3) = dssnum(jl1:jl2,n_o3) +           &
                    rk(jl1:jl2,i)                                       &
                    *y(jl1:jl2,nspi(i,1))

! add terms for derivative w.r.t OH
              IF (nspi(i,1) == nspo3p)                                  &
! add terms to derivates for d(a[A][B])
                dssnum(jl1:jl2,n_oh) = dssnum(jl1:jl2,n_oh) +           &
                    rk(jl1:jl2,i)                                       &
                    *y(jl1:jl2,nspi(i,2))*deriv(jl1:jl2,nss_o3p,n_oh)

              IF (nspi(i,2) == nspo3p)                                  &
! add terms to derivates for d(a[O2][O1D])/d[O3] and b[N2][O1D]
                dssnum(jl1:jl2,n_oh) = dssnum(jl1:jl2,n_oh) +           &
                    rk(jl1:jl2,i)                                       &
                    *y(jl1:jl2,nspi(i,1))*deriv(jl1:jl2,nss_o3p,n_oh)

              IF (nspi(i,1) == nspoh)                                   &
! add terms to derivates for d(a[O2][O1D])/d[O3] and b[N2][O1D]
                dssnum(jl1:jl2,n_oh) = dssnum(jl1:jl2,n_oh) +           &
                    rk(jl1:jl2,i)                                       &
                    *y(jl1:jl2,nspi(i,2))

              IF (nspi(i,2) == nspoh)                                   &
! add terms to derivates for d(a[O2][O1D])/d[O3] and b[N2][O1D]
                dssnum(jl1:jl2,n_oh) = dssnum(jl1:jl2,n_oh) +           &
                    rk(jl1:jl2,i)                                       &
                    *y(jl1:jl2,nspi(i,1))

! add terms for derivative w.r.t HO2
              IF (nspi(i,1) == nspo3p)                                  &
                dssnum(jl1:jl2,n_ho2) = dssnum(jl1:jl2,n_ho2) +         &
                    rk(jl1:jl2,i)                                       &
                    *y(jl1:jl2,nspi(i,2))*deriv(jl1:jl2,nss_o3p,n_ho2)

              IF (nspi(i,2) == nspo3p)                                  &
                dssnum(jl1:jl2,n_ho2) = dssnum(jl1:jl2,n_ho2) +         &
                    rk(jl1:jl2,i)                                       &
                    *y(jl1:jl2,nspi(i,1))*deriv(jl1:jl2,nss_o3p,n_ho2)

              IF (nspi(i,1) == nspho2)                                  &
                dssnum(jl1:jl2,n_ho2) = dssnum(jl1:jl2,n_ho2) +         &
                    rk(jl1:jl2,i)                                       &
                    *y(jl1:jl2,nspi(i,2))

              IF (nspi(i,2) == nspho2)                                  &
                dssnum(jl1:jl2,n_ho2) = dssnum(jl1:jl2,n_ho2) +         &
                    rk(jl1:jl2,i)                                       &
                    *y(jl1:jl2,nspi(i,1))

! add terms for derivative w.r.t NO
              IF (nspi(i,1) == nspno)                                   &
                dssnum(jl1:jl2,n_no) = dssnum(jl1:jl2,n_no) +           &
                    rk(jl1:jl2,i)                                       &
                    *y(jl1:jl2,nspi(i,2))

              IF (nspi(i,2) == nspno)                                   &
                dssnum(jl1:jl2,n_no) = dssnum(jl1:jl2,n_no) +           &
                    rk(jl1:jl2,i)                                       &
                    *y(jl1:jl2,nspi(i,1))

            END IF
          END IF
        END DO       ! jr
!
! Destruction terms
        DO jr = 1,nssrt(ix)
          i = nssri(ix,jr)
          j = nssrx(ix,jr)
          IF (i <= nuni) THEN
            ssden(jl1:jl2) = ssden(jl1:jl2) + rk(jl1:jl2,i)
          ELSE
            ssden(jl1:jl2) = ssden(jl1:jl2) +                           &
              rk(jl1:jl2,i) * y(jl1:jl2,nspi(i,j))
            IF (ix < 5) THEN
              IF (nspi(i,j) == nspo3 )                                  &
                dssden(jl1:jl2,n_o3 ) = dssden(jl1:jl2,n_o3 ) +         &
                  rk(jl1:jl2,i)
              IF (nspi(i,j) == nspoh )                                  &
                dssden(jl1:jl2,n_oh ) = dssden(jl1:jl2,n_oh ) +         &
                  rk(jl1:jl2,i)
              IF (nspi(i,j) == nspho2)                                  &
                dssden(jl1:jl2,n_ho2) = dssden(jl1:jl2,n_ho2) +         &
                  rk(jl1:jl2,i)
              IF (nspi(i,j) == nspno )                                  &
                dssden(jl1:jl2,n_no ) = dssden(jl1:jl2,n_no ) +         &
                  rk(jl1:jl2,i)
            END IF
          END IF
        END DO    ! jr
!
! Steady state and derivatives of steady state
        y(jl1:jl2,nssi(ix)) = ssnum(jl1:jl2)/ssden(jl1:jl2)
        IF (ix < 5) THEN
          DO jr =1,4
            deriv(jl1:jl2,ix,jr) =                                      &
                (ssden(jl1:jl2)*dssnum(jl1:jl2,jr) -                    &
                 ssnum(jl1:jl2)*dssden(jl1:jl2,jr))/                    &
                 (ssden(jl1:jl2) * ssden(jl1:jl2))
          END DO    ! jr
        END IF
      END DO  ! ix

! rescale deriv to mean [O3]/[O] * d[O]/d[O3], where [O] = [O(1D)] or [O(3P)]
! for O(1D), and O(3P), N and H when these are SS species

      !!! these where statements do not work if species NO is not present in    
      !!! the sheme, i.e. nspno==0; e.g. scheme Pommrich et al, 2014.           
      !!! in this case, no steady state species are involved (jug)              
      !!! thus added switch o1d_in_ss here to check whether o1d is ss species   

      IF (o1d_in_ss) THEN
         WHERE (y(jl1:jl2,nspo1d) > peps)
            deriv(jl1:jl2,nss_o1d,n_o3) = deriv(jl1:jl2,nss_o1d,n_o3 )*     &
                                      y(jl1:jl2,nspo3 )/y(jl1:jl2,nspo1d)
            deriv(jl1:jl2,nss_o1d,n_oh) = deriv(jl1:jl2,nss_o1d,n_oh )*     &
                                      y(jl1:jl2,nspoh )/y(jl1:jl2,nspo1d)
            deriv(jl1:jl2,nss_o1d,n_ho2)= deriv(jl1:jl2,nss_o1d,n_ho2)*     &
                                      y(jl1:jl2,nspho2)/y(jl1:jl2,nspo1d)
            deriv(jl1:jl2,nss_o1d,n_no )= deriv(jl1:jl2,nss_o1d,n_no )*     &
                                      y(jl1:jl2,nspno )/y(jl1:jl2,nspo1d)
         ELSEWHERE
            deriv(jl1:jl2,nss_o1d,n_o3)  = 1.
            deriv(jl1:jl2,nss_o1d,n_oh)  = 1.
            deriv(jl1:jl2,nss_o1d,n_ho2) = 1.
            deriv(jl1:jl2,nss_o1d,n_no)  = 1.
         ENDWHERE
      END IF

      IF (o3p_in_ss) THEN
        WHERE (y(jl1:jl2,nspo3p) > peps)
          deriv(jl1:jl2,nss_o3p,n_o3 )= deriv(jl1:jl2,nss_o3p,n_o3 )*   &
                                  y(jl1:jl2,nspo3 )/y(jl1:jl2,nspo3p)
          deriv(jl1:jl2,nss_o3p,n_oh )= deriv(jl1:jl2,nss_o3p,n_oh )*   &
                                  y(jl1:jl2,nspoh )/y(jl1:jl2,nspo3p)
          deriv(jl1:jl2,nss_o3p,n_ho2)= deriv(jl1:jl2,nss_o3p,n_ho2)*   &
                                  y(jl1:jl2,nspho2)/y(jl1:jl2,nspo3p)
          deriv(jl1:jl2,nss_o3p,n_no )= deriv(jl1:jl2,nss_o3p,n_no )*   &
                                  y(jl1:jl2,nspno )/y(jl1:jl2,nspo3p)
        ELSEWHERE
          deriv(jl1:jl2,nss_o3p,n_o3)  = 1.
          deriv(jl1:jl2,nss_o3p,n_oh)  = 1.
          deriv(jl1:jl2,nss_o3p,n_ho2) = 1.
          deriv(jl1:jl2,nss_o3p,n_no)  = 1.
        ENDWHERE
      END IF

      IF (n_in_ss) THEN
        WHERE (y(jl1:jl2,nspn  ) > peps)
          deriv(jl1:jl2,nss_n,n_o3) = deriv(jl1:jl2,nss_n,n_o3 )*       &
                                y(jl1:jl2,nspo3 )/y(jl1:jl2,nspn)
          deriv(jl1:jl2,nss_n,n_oh) = deriv(jl1:jl2,nss_n,n_oh )*       &
                                y(jl1:jl2,nspoh )/y(jl1:jl2,nspn)
          deriv(jl1:jl2,nss_n,n_ho2)= deriv(jl1:jl2,nss_n,n_ho2)*       &
                                y(jl1:jl2,nspho2)/y(jl1:jl2,nspn)
          deriv(jl1:jl2,nss_n,n_no )= deriv(jl1:jl2,nss_n,n_no )*       &
                                y(jl1:jl2,nspno )/y(jl1:jl2,nspn)
        ELSEWHERE
          deriv(jl1:jl2,nss_n,n_o3)  = 1.
          deriv(jl1:jl2,nss_n,n_oh)  = 1.
          deriv(jl1:jl2,nss_n,n_ho2) = 1.
          deriv(jl1:jl2,nss_n,n_no)  = 1.
        ENDWHERE
      END IF

      IF (h_in_ss) THEN
        WHERE (y(jl1:jl2,nsph  ) > peps)
          deriv(jl1:jl2,nss_h,n_o3) = deriv(jl1:jl2,nss_h,n_o3 )*       &
                              y(jl1:jl2,nspo3 )/y(jl1:jl2,nsph  )
          deriv(jl1:jl2,nss_h,n_oh) = deriv(jl1:jl2,nss_h,n_oh )*       &
                              y(jl1:jl2,nspoh )/y(jl1:jl2,nsph  )
          deriv(jl1:jl2,nss_h,n_ho2)= deriv(jl1:jl2,nss_h,n_ho2)*       &
                              y(jl1:jl2,nspho2)/y(jl1:jl2,nsph  )
          deriv(jl1:jl2,nss_h,n_no )= deriv(jl1:jl2,nss_h,n_no )*       &
                              y(jl1:jl2,nspno )/y(jl1:jl2,nsph  )
        ELSEWHERE
          deriv(jl1:jl2,nss_h,n_o3)  = 1.
          deriv(jl1:jl2,nss_h,n_oh)  = 1.
          deriv(jl1:jl2,nss_h,n_ho2) = 1.
          deriv(jl1:jl2,nss_h,n_no)  = 1.
        ENDWHERE
      END IF

      RETURN
      END SUBROUTINE asad_steady


    End Module messy_clamschem_asad_steady

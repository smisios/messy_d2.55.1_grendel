Module messy_clamschem_asad_impact

contains

!**** asad_impact - driver routine for IMPACT integration scheme
!
!     Glenn Carver,                Centre for Atmospheric Science,
!                                  University of Cambridge.
!
!
! Purpose: Controlling routine for the IMPACT integration scheme
!
!
!          Called from ASAD_CDRIVE
!
!
!     Interface
!     ---------
!     Called from chemistry driver routine asad_cdrive.
!
!     Arguments:
!          lphot - Logical array. One entry per gridpoint. If set .T.
!                  then photolysis is occurring and the scheme will
!                  adjust the accuracy of the Jacobian. See method.
!
!     All other variables are obtained via COMMON.
!
!     Method
!     ------
!     The IMPACT scheme used the Euler backward implicit timescheme
!     as its starting point. This is a first order conservative
!     scheme. A predictor-corrector approach is used to solve this
!     by employing a Newton-Raphson type iteration. However, rather
!     than compute and invert a near full Jacobian matrix, this scheme
!     only considers the diagonal terms only. For mildly stiff chemistry
!     this should work fine. However, if the scheme is not converging
!     rapidly enough during daylight, this is a symptom that the
!     off-diagonal production terms, which are usually ignored, are
!     important. In this case, set the switch ljacx .T. and make sure
!     that the lphot array is set .true. where photolysis is on.
!     Since, only approximate forms of the Jacobian are considered, if
!     switched on, the method will add in the terms in the Jacobian
!     due to photolysis. It does this in an approximate way to avoid
!     an matrix inversion but this greatly improves the performance of
!     the scheme under stiffer conditions.
!
!     The IMPACT scheme is fully described in the forthcoming paper;
!     Carver and Stott, 1997, 'IMPACT: an implicit time integration
!     scheme for chemical species and families', to be submitted to
!     Annales Geophysicae. Preprints available.
!
!     Externals
!     ---------
!     ftoy         Converts f concentrations to y (partitions
!                  families).
!     diffun       Computes tendencies due to chemistry and
!                  main diagonal Jacobian elements.
!
!     Local variables
!     ---------------
!     gconv         .true. if convergence in Newton-Rhapson has been
!                   achieved.
!     ifi           Number of ftoy iterations.
!     zf            Values of f at start of chemistry step.
!     zprf          Value of f on previous Newton-Rhapson iteration.
!     binv          Inverse of the approximate (main diagonal only)
!                   of the matrix ( I - dt.J) where J is Jacobian.
!     dely          Represents the change in variable, y, produced by
!                   multiplying binv to the function we are trying to
!                   find the root of.
!
! Code description:
!   Language: FORTRAN 90
!
! ---------------------------------------------------------------------
!
        SUBROUTINE ASAD_IMPACT(gfirst,n_points)

        USE messy_clamschem_asad_mod
        USE messy_clamschem_asad_mod_clams, ONLY: jpctr, theta_field_size, failed
        USE messy_clams_global,             ONLY: prec
        USE messy_clamschem_asad_diffun,    ONLY: ASAD_DIFFUN
        USE messy_clamschem_asad_ftoy,      ONLY: ASAD_FTOY
        USE messy_clamschem_asad_inimpct,   ONLY: ASAD_INIMPCT
        USE messy_clamschem_asad_jac,       ONLY: ASAD_JAC

        IMPLICIT NONE

!!!!! 
! ju_nt_20140407
        logical :: gfirst

        INTEGER, INTENT(IN ) :: n_points   ! No of spatial points

!       Local variables

        INTEGER :: lphot(theta_field_size)
        INTEGER :: inl
        INTEGER :: nl
        INTEGER :: j                        ! Loop variable
        INTEGER :: jit                      ! Loop variable
        INTEGER :: jl                       ! Loop variable
        INTEGER :: jr                       ! Loop variable
        INTEGER :: jtr                      ! Loop variable
        INTEGER :: jt1                      ! Loop variable
        INTEGER :: isp                      ! Index
        INTEGER :: ifi
        INTEGER :: ntr1
        INTEGER :: ntr2
        INTEGER :: ipos0
        INTEGER :: ipos
        INTEGER :: ir
        INTEGER :: irk
        INTEGER :: njr
        INTEGER :: ireac
        INTEGER :: iprod

        REAL(PREC) :: zf(theta_field_size,jpctr)
        REAL(PREC) :: zprf(theta_field_size,jpctr)
        REAL(PREC) :: binv(theta_field_size,jpctr)
        REAL(PREC) :: dely(theta_field_size,jpctr)
        REAL(PREC) :: dd(theta_field_size)
        REAL(PREC) :: corrn(theta_field_size)

        LOGICAL :: gconv
        LOGICAL :: gfam
! ju_nt_20140407
!!!!! Funktioniert so nicht, da gfirst in jedem Messy-Zeitschritt zurueckgesetzt werden muss
!        LOGICAL, SAVE :: gfirst = .true.


!       1.  Predictor step (first guess)
!           --------- ---- ------ ------

!       Crude way of determining if photolysis is on or not.
!       Needs improving.

        inl = 0
        DO jl = 1, n_points
          IF ( rk(jl,nprkx(1)) > peps ) THEN
            inl = inl + 1
            lphot(inl) = jl
          endif
        ENDDO

        IF ( gfirst ) THEN
          gfirst = .false.
! DEPENDS ON: asad_inimpct
          CALL ASAD_INIMPCT
        ENDIF

!       1.1  Do the linearised first guess to give first approx.
!            solution at y(n+1).

        nl = n_points
        DO jtr = 1, jpctr
          isp = majors(jtr)
          DO jl = 1, n_points
            zf(jl,jtr) = f(jl,jtr)

!           1.2  Test Jacobian has not gone positive. Possible since
!                we only use approximate form. If it has, use
!                an explicit step.

            if ( ej(jl,jtr) > 0.0 ) then
              f(jl,jtr) = zf(jl,jtr) + cdt*fdot(jl,jtr)
            else
              f(jl,jtr) = zf(jl,jtr) + ( cdt*fdot(jl,jtr) )            &
                               / ( 1.0 - cdt*ej(jl,jtr) )
            endif
            if( linfam(jl,jtr) ) f(jl,jtr) = y(jl,isp)

          ENDDO
        ENDDO


!       2.  Newton-Rhapson iteration (corrector step).
!           -------------- --------- ---------- ------

        DO jit = 1, nrsteps

!         2.1  Decide if we are recomputing ratios or not.

! ju_nt_20140313: if-statement changed!
!!!!! ACHTUNG: Abfrage wieder wie im "alten" ASAD-Code !!!
!          IF ( jit < 4 ) THEN
          IF ( jit /= 1 ) THEN
            ifi = 0
          ELSE
            ifi = nitnr
          ENDIF

!         2.2  Work out rates of change and main diagonal of J.

! DEPENDS ON: asad_ftoy
          CALL ASAD_FTOY(.false., ifi, n_points)
! DEPENDS ON: asad_diffun
          CALL ASAD_DIFFUN( nl )
! DEPENDS ON: asad_jac
          CALL ASAD_JAC(n_points)

!         2.3  Do the normal Newton-Raphson iteration with
!              Jacobian approximated to main diagonal.

          DO jtr = 1, jpctr
            isp = majors(jtr)
            DO jl = 1, n_points
              zprf(jl,jtr) = f(jl,jtr)
              binv(jl,jtr) = 1.0 / ( 1.0 - cdt*ej(jl,jtr) )
!!!!!
!!$if (f(jl,jtr)<-1.E28 .or. fdot(jl,jtr)<-1.E28) then
!!$   write (*,*) 'missing_value: ', jl, jtr, f(jl,jtr), fdot(jl,jtr), missing_index(jl)
!!$endif
              dely(jl,jtr) = binv(jl,jtr) *                            &
                     ( (zf(jl,jtr) - f(jl,jtr)) + cdt*fdot(jl,jtr) )
              f(jl,jtr) = f(jl,jtr) + dely(jl,jtr)

              if( linfam(jl,jtr) ) f(jl,jtr) = y(jl,isp)
            ENDDO
          ENDDO

!         2.4  If requested for all points where photolysis is
!              occurring, include the terms in the Jacobian
!              arising from just the photolysis terms. Note
!              that we only do this for points where photolysis
!              is actually turned on.

          IF ( ljacx .AND. inl > 0 ) THEN
            ntr1 = nltrim(0,1)

!           2.4.1  For each tracer that needs correcting (ie.
!                  loop over the rows in the Jacobian matrix).

          DO jt1 = 1, ntr1

!             2.4.2  Compute the contribution from non-zero
!                    elements in the Jacobian due to photolysis.
!                    ie. work our way along the columns. We first
!                    compute the contribution from ordinary tracers
!                    'TR' and then from families.

              ntr2  = nltrim(jt1,2)
              ipos0 = nltrim(jt1,3)
              ipos  = ipos0
              DO jl = 1, n_points
                corrn(jl) = 0.0
              ENDDO

              DO

              ir   = nlpdv(ipos,1)
              irk  = nprkx(ir)
              isp  = nspi(irk,1)
              njr  = nlpdv(ipos,2)
              ireac = madvtr(isp)
              gfam = ireac  ==  0
              DO j = 1, inl
                dd(j) = 0.0
              ENDDO

              IF ( gfam ) THEN
                ireac = moffam(isp)
                DO jr = 1, njr
                  ir  = nlpdv(ipos,1)
                  irk = nprkx(ir)
                  isp = nspi(irk,1)
                  DO j = 1, inl
                    jl    = lphot(j)
                    dd(j) = dd(j) + rk(jl,irk)*y(jl,isp)
                  ENDDO
                  ipos = ipos + 1
                ENDDO
                DO j = 1, inl
                  jl    = lphot(j)
                  dd(j) = dd(j) / zprf(jl,ireac)
                ENDDO

              ELSE

!             2.4.3.  Species of type 'FT' will come here.
!                     Check, at each point, whether it's gone
!                     into the family or not.

                DO jr = 1, njr
                  ir  = nlpdv(ipos,1)
                  irk = nprkx(ir)
                  isp = nspi(irk,1)
                  DO j = 1, inl
                    jl = lphot(j)
                    if ( linfam(jl,ireac) ) then
                      ! ju_jug 08/15 check zprf to avoid division by zero 
                      if(zprf(jl,ireac) > 0.) &
                           dd(j) = dd(j)+rk(jl,irk)*y(jl,isp)/zprf(jl,ireac)
                    else
                      dd(j) = dd(j)+rk(jl,irk)
                    endif
                  ENDDO
                  ipos = ipos + 1
                ENDDO
              ENDIF

              DO j = 1, inl
                jl = lphot(j)
                corrn(j) = corrn(j) + dd(j)*dely(jl,ireac)
              ENDDO

              IF ( ipos - ipos0 >= ntr2 ) EXIT
              ENDDO

!             2.4.3  Now add correction to the tracer (rows). If the
!                    tracer is of type 'FT', then for each pt, we must
!                    check whether the tracer has been put into the f
!                    or not. If it has, we add the correction to the
!                    family and not to the tracer.

              iprod = nltrim(jt1,1)
              isp   = majors(iprod)
              IF ( ctype(isp) == jpif ) THEN
                DO j = 1, inl
                  jl    = lphot(j)
                  iprod = nltrim(jt1,1)
                  IF ( linfam(jl,iprod) ) iprod = moffam(isp)
                  f(jl,iprod) = f(jl,iprod)+cdt*binv(jl,iprod)*corrn(j)

                  ! ju_jug_20150911  output warning in case of neg. conc.
                  if (f(jl,iprod) <0.) then
                     write(*,*)'*** IMPACT warning: negative concentration: ', &
                       speci(isp),jl, iprod, f(jl,iprod), zf(jl,iprod), corrn(j)
                  endif
                ENDDO
              ELSE
                DO j = 1, inl
                  jl = lphot(j)
                  f(jl,iprod) = f(jl,iprod)+cdt*binv(jl,iprod)*corrn(j)

                  ! ju_jug_20150911 warning for neg. conc.
                  if (f(jl,iprod) <0.) then
                     write(*,*)'*** IMPACT warning: negative concentration: ', &
                       speci(isp),jl, iprod, f(jl,iprod), zprf(jl,iprod), corrn(j)
                  endif
                ENDDO
              ENDIF

            ENDDO
          ENDIF      ! End of IF ( ljacx .AND. inl > 0 ) statement

!         9. Check for convergence
!            ----- --- -----------

          gconv = .true.
          DO jtr = 1,jpctr
            DO jl = 1, n_points
              IF ( ABS(f(jl,jtr)-zprf(jl,jtr)) >  ptol*f(jl,jtr)       &
              .AND. f(jl,jtr) >  pmintnd(jl) ) gconv=.false.
! ju_nt_20140212
!     Check for NaNs and unlikely large values (jug, 04/2013)
              if (f(jl,jtr) /= f(jl,jtr)) failed(jl)=.true.
              !if (f(jl,jtr) >= tnd(jl)) failed(jl)=.true.
            ENDDO
            IF ( .NOT. gconv ) EXIT
          ENDDO
          IF (gconv) THEN
            RETURN
          END IF
        ENDDO           ! End of NR (jit) loop

!       9.1  Convergence achieved or max. no. of iterations reached

        RETURN
        END SUBROUTINE ASAD_IMPACT


      End Module messy_clamschem_asad_impact

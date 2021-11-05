Module messy_clamschem_asad_ycn

contains

!
!     YCN  - computes dy/dt for stiff integrators.
!
!     Glenn Carver             Centre for Atmospheric Science
!                              University of Cambridge
!
!
! Purpose: This routine is intended for use with stiff integrators which
!     require the evaluation of df/dt at any point in time.
!
!     This routine is intended to be passed as an argument to the stiff
!     integrators such as the NAG and SVODE drivers.
!
!
!          Called from ASAD_DIFFUN
!
!
!     Interface
!     ---------
!        tz     - on entry, specifies the time (unchanged).
!        f1     - on entry, contains species at time TZ (unchanged).
!        dfdt   - on exit, must contain time derivative of the gridpt.
!
!     As all the chemistry arrays are in common, we have to copy the
!     input to the common array and copy the tendency to the output
!     argument at the end.
!
!     Method
!     ------
!     IMPORTANT!!! This subroutine assumes that only a single gridpt
!     is being worked on by the integrator calling this subroutine.
!
!     The first step is to copy the passed values of the
!     species back into the ASAD common blocks so that the production
!     and loss can be computed. Since the integrator is only integrating
!     the variables, f, passed between the model and ASAD we only copy t
!     the species that will change during the timestep. Note that we onl
!     copy to the first element in the species array, y.
!
!     When no families are in use, we can use the index array nlf since
!     this stores a list of all the species of type TR and nf = jpctr.
!
!     The next step is to compute the rates of change. The routine diffu
!     is used but we tell it only to compute the first gridpt. On exit
!     from diffun, the fdot array will have been assigned.
!
!     Externals
!     ---------
!       diffun    - to compute rates of change.
!
! Code description:
!   Language: FORTRAN 90
!
! ---------------------------------------------------------------------
!
        SUBROUTINE ASAD_YCN(tz,f1,dfdt)

        USE messy_clams_global,          ONLY: prec
        USE messy_clamschem_asad_mod,    ONLY: nf, nlf, y, fdot, method, f, nitfg
        USE messy_clamschem_asad_dummy,  ONLY: ereport
        USE messy_clamschem_asad_diffun, ONLY: ASAD_DIFFUN
        USE messy_clamschem_asad_ftoy,   ONLY: ASAD_FTOY
        USE messy_clamschem_asad_mod_clams, ONLY: jpctr, jl_current, theta_field_size

        IMPLICIT NONE


        REAL(PREC), INTENT(IN)  :: tz
        REAL(PREC), INTENT(IN)  :: f1(jpctr)

        REAL(PREC), INTENT(OUT) :: dfdt(jpctr)

!       Local variables

        INTEGER :: j    ! Loop variable
        INTEGER :: js

! ju_nt_20140304: tzold and fits added
        real, save :: tzold = 0.0
        integer    :: fits
        
        LOGICAL, SAVE :: gfirst = .TRUE.

        CHARACTER(LEN=72) :: cmessage  ! Error message


!       1. Initialise species array; ONLY COPY TO FIRST ELEMENT!
!          ---------- ------- ------ ---- ---- -- ----- --------

        IF ( gfirst ) THEN
! ju_nt_20140304
!!$          IF ( method  >=  10 .and. jpctr  /=  nf ) then
!!$            WRITE (6,*) '** INTERNAL ASAD ERROR: jpctr  /=  nf in ycn',&
!!$            ' There should not be any families in use with the stiff', &
!!$            ' integrators.'
!!$            cmessage = 'jpctr  /=  nf in ycn'
!!$
!!$            CALL EREPORT('ASAD_YCN',jpctr,cmessage)
!!$         ENDIF
         gfirst = .false.
        ENDIF

! ju_nt_20140305 : if (nf==jpctr) added:
!       n.b. nf = jpctr when no families are in use
!         if families are in use then nf > jpctr

        if (nf == jpctr) then
           DO j = 1, nf
              js  = nlf(j)
! ju_nt_20140305:  1 -> jl_current 
!$$              y(1,js) = f1(j)
             y(jl_current,js) = f1(j)
           ENDDO

        else ! ju_nt_20140305: else-block added:
           !  put f values into y array but keep ratios constant.
           !!!!!!!!! should check for svode integrator !!!!!!
           ! tzold gets left with the value of tn from the previous timestep
           ! which will force the ratios to be recomputed. This is probably
           ! wrong since ftoy will have already been called in cdrive.f but
           ! it shouldn't make any difference - it's just a repeated
           ! computation.
           do j = 1, jpctr
              ! changed the following line 1 -> jl_current (jug, 15.5.2006)
              f(jl_current,j) = f1(j)
           end do
           if ( tz /= tzold ) then
              fits = nitfg
           else
              fits = 0
           endif
           tzold = tz
           ! this currently does all gridpoints!!!
           call asad_ftoy( .false., fits, theta_field_size )
        endif
        
!       2.   Compute rates of change.
!            ------- ----- -- -------

! DEPENDS ON: asad_diffun
! ju_nt_20140305:  1 -> jl_current 
!$$        CALL ASAD_DIFFUN( 1 )
        CALL ASAD_DIFFUN( jl_current )

        DO j = 1, jpctr
! ju_nt_20140305:  1 -> jl_current 
!$$          dfdt(j) = fdot(1,j)
           dfdt(j) = fdot(jl_current,j)
        ENDDO

        RETURN
        END SUBROUTINE ASAD_YCN
!--------------------------------------------------------------------

      End Module messy_clamschem_asad_ycn

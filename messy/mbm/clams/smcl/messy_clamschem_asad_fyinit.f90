Module messy_clamschem_asad_fyinit

contains

!****  asad_fyinit - sets or initialises species concentrations
!
!     Glenn Carver & Paul D.Brown, Centre for Atmospheric Science,
!                                  University of Cambridge.
!
!
! Purpose: sets or initialises species concentrations
!
!     Sets the individual species concentrations y from those of
!     the model tracers/families f. If the species concentrations
!     are later to be iterated in ftoy, then they are merely
!     initialised in this routine.
!
!
!          Called from ASAD_FTOY
!
!
!     Arguments:
!        ofirst     Should be .true. if this is first call to FTOY
!                   for the chemistry time loop ie. ofirst should be
!                   true only when starting a new dynamical timestep.
!
!     Method
!     ------
!     If this is the first call to ftoy for the chemical time loop
!     ie. we're starting a new model/dynamical timestep, then we need
!     to initialise the species array, y. For subsequent calls, all we
!     need to do is copy the species of type 'TR' and 'FT' from the f
!     (family/tracer) array to their places in the species array, y.
!
!     When ofirst is .true., family are initialised by setting the
!     major family member to the family concentration and the minor
!     family members to a small number. Steady state species are also
!     set to a small number. Constant species (species type 'CT')
!     are also only set when .ofirst. is true, since these species are
!     never integrated (see diffun.f). Constant field species (type
!     'CF') are allowed to change every chemical time step however.
!
!     Local variables
!     ---------------
!     zcnst          VMR of fixed concentration species.
!
! Code description:
!   Language: FORTRAN 90
!
! ---------------------------------------------------------------------
!
        SUBROUTINE ASAD_FYINIT(ofirst,n_points)

        USE messy_clamschem_asad_mod
        USE messy_clamschem_asad_mod_clams, only: jpspec
        USE messy_clamschem_asad_dummy,     ONLY: ereport
        USE messy_clams_global,             only: prec
        USE messy_clamschem_get_const_2d,   ONLY: get_const_2d
        USE messy_clamschem_asad_inicnt,    ONLY: ASAD_INICNT
        
        IMPLICIT NONE

        INTEGER, INTENT(IN) :: n_points   ! No of spatial points

        LOGICAL, INTENT(IN) :: ofirst     ! True on first call

!       Local variables

        INTEGER, SAVE :: iss
        INTEGER, SAVE :: ict
        INTEGER, SAVE :: iftr
        INTEGER, SAVE :: icf
        INTEGER, SAVE :: istep = 0

        INTEGER       :: j           ! Loop variable
        INTEGER       :: jl          ! Loop variable
        INTEGER       :: js          ! Index
        INTEGER       :: istart
        INTEGER       :: iend
        INTEGER       :: ifam
        INTEGER       :: itr
        INTEGER       :: kstep
        INTEGER       :: inl

        INTEGER       :: errcode    ! variable passed to ereport

        REAL(PREC)          :: peps10
        REAL(PREC)          :: zcnst

        LOGICAL, SAVE :: gonce = .true.
        CHARACTER(len=72) :: cmessage

        INTEGER       :: status = 0
        LOGICAL       :: mixrat = .true.
        REAL(PREC), dimension(n_points) :: zcnst_arr


!       1.  Species initialised only at start of dynamical step.
!           ------- ----------- ---- -- ----- -- --------- -----

        IF ( gonce ) THEN

!         Build lists for this routine. This lot should really be a
!         common variable - next version!

          gonce = .false.
          iss = 0
          ict = 0
          icf = 0
          iftr = 0
          DO js = 1, jpspec
            IF ( ctype(js) == jpna ) THEN
              iss       = iss + 1
              ilss(iss) = js
            END IF
            IF ( ctype(js) == jpco ) THEN
              ict       = ict + 1
              ilct(ict) = js
            END IF
            IF ( ctype(js) == jpcf ) THEN
              icf       = icf + 1
              ilcf(icf) = js
            END IF
            IF (ctype(js) == jpif .OR.                                 &
                ctype(js) == jpsp ) THEN
              iftr        = iftr + 1
              ilftr(iftr) = js
            END IF
          END DO
        END IF           ! End of IF (gonce) statement

        IF ( ofirst ) THEN
          peps10 = 10.0 * peps

!         1.1 Ordinary family members (major then minor )

          istart = nlmajmin(1)
          iend   = nlmajmin(2)
          DO j = istart, iend
            js = nlmajmin(j)
            ifam = moffam(js)
            DO jl = 1, n_points
              y(jl,js) = f(jl,ifam)
            END DO
          END DO

          istart = nlmajmin(3)
          iend   = nlmajmin(4)
          DO j = istart, iend
            js = nlmajmin(j)
            DO jl = 1, n_points
              y(jl,js) = peps10
            END DO
          END DO

!         1.2 Non-model species in steady state

          DO j = 1, iss
            js = ilss(j)
            DO jl = 1, n_points
              y(jl,js) = peps10
            END DO
          END DO

!         1.3 Species set to constants (unaffected by remainder of ftoy

!!!!!
! ju_nt_20140210
          mixrat=.true.

          DO j = 1, ict
            js = ilct(j)
            IF ( speci(js) == 'N2        ' ) THEN
              zcnst = fn2
            ELSE IF ( speci(js) == 'O2        ' ) THEN
              zcnst = fo2
            ELSE IF ( speci(js) == 'CO2       ' ) THEN
              zcnst = fco2
            ELSE IF ( speci(js) == 'H2        ' ) THEN
              zcnst = fh2
            ELSE IF ( speci(js) == 'CH4       ' ) THEN
              zcnst = fch4
!!!!!
! ju_nt_20140210
! jug, 09/2006:
! Einbau "Simple Chemie", d.h. OH, Cl und O(1D) können aus zonalen Mittelwerten
! z.B. aus dem 2-D Modell interpoliert werden...
            ELSE IF ( speci(js) == 'OH        ' .or. &
                      speci(js) == 'O(1D)     ' .or. &
                      speci(js) == 'Cl        ' .or. &
                      speci(js) == 'HO2       ' ) THEN
               status = 0
               call get_const_2d (status,speci(js),zcnst_arr)
               if (status /= 0) CALL EREPORT('ASAD_FYINIT',status,'Error in get_const_2d')
               y(1:n_points,js) = zcnst_arr(1:n_points)
               mixrat=.false.

            ELSE
              cmessage=' Value not supplied for '//speci(js)//          &
                       ',setting value to zero'
              errcode=-1

              CALL EREPORT('ASAD_FYINIT',errcode,cmessage)
              zcnst = 0.0
            END IF

!!!!!
! ju_nt_20140210
            if (mixrat)  then
               DO jl = 1, n_points
                  y(jl,js) = zcnst * tnd(jl)
               END DO
            endif
          END DO
        END IF       ! ofirst

!       2.   Initialise rest of species (potentially each call).
!            ---------- ---- -- ------- ------------ ---- ------

!       2.1  Species of type 'TR' and 'FT' are copied from tracers
!            on every call to this routine.

        DO j = 1, iftr
          js = ilftr(j)
          itr = madvtr(js)
          DO jl = 1, n_points
            y(jl,js) = f(jl,itr)
          END DO
        END DO

!       2.2  Constant but set by user, type 'CF' are 
!            initialised at every chemical step.
!            Currently only for H2O.

        kstep = jsubs
        inl = n_points
        DO j = 1, icf
          js = ilcf(j)
! DEPENDS ON: asad_inicnt
          CALL ASAD_INICNT( speci(js), y(1,js), inl )
        END DO

        RETURN
        END SUBROUTINE ASAD_FYINIT


      End Module messy_clamschem_asad_fyinit

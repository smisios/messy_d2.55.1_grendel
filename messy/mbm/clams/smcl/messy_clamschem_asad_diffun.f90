Module messy_clamschem_asad_diffun

contains

!**** asad_diffun - calculates tendencies due to chemistry
!
!     Glenn Carver & Paul D Brown, Centre for Atmospheric Science,
!                                  University of Cambridge.
!
!
! Purpose: Subroutine to calculates tendencies due to chemistry
!          of the model tracers/families, by adding or equivalencing
!          the tendencies of individual chemical species.
!
!          Called from ASAD_CDRIVE, ASAD_IMPACT, and
!                      ASAD_YCN
!
!     Interface
!     ---------
!        kl   - Inner loop limit over spatial points.
!               Usually this is jpnl but some stiff
!               integrators will compute one gridpt at
!        a time in which case kl is set to 1.  Note however
!        that this severely hampers performance on a vector
!        computer since the inner loop is vectorized.
!
!     Method
!     ------
!     For species which are in or out of the family depending on
!     their lifetime, their tendency is only added to the family
!     if it is in the family. Its species tendency is always copied
!     to the fdot array as we always integrate it. If the species
!     is in the family however, the integrated result will be
!     overwritten with the steady state value.
!
!     If dy/dt = P - Ly where P is the production and L is the loss
!     rate, then the array prod holds 'P' and the array slos holds
!     'Ly'. Both P and L are positive.
!
!     Since we only need to integrate those species used in the
!     calling model, we can ignore non-family steady state species
!     and species that are constant. A list of such species is held
!     in nlf; actually twice because it's used for both product and
!     loss calculations (see prls) so we only need to loop over the
!     first half of the array. By ignoring these species we can
!     save computer time.
!
!     Externals
!     ---------
!     prls      Calculates production and loss terms of
!               individual species.
!
! Code description:
!   Language: FORTRAN 90
!
! ---------------------------------------------------------------------
!
        SUBROUTINE asad_diffun( kl )

        USE messy_clamschem_asad_mod, ONLY: fdot, ydot, prod, slos,       &
                                            linfam, nodd, moffam, madvtr, &
                                            nf, nlf
        USE messy_clamschem_asad_mod_clams,  only: jpctr
        USE messy_clamschem_asad_prls, ONLY: ASAD_PRLS

        IMPLICIT NONE

        INTEGER, INTENT(IN) :: kl      ! No of spatial points

!       Local variables

        INTEGER :: ifam                ! Index
        INTEGER :: itr                 ! Index
        INTEGER :: j                   ! Loop variable
        INTEGER :: jl                  ! Loop variable
        INTEGER :: jtr                 ! Loop variable
        INTEGER :: js                  ! Index

        LOGICAL :: gfam
        LOGICAL :: gtr
        LOGICAL, SAVE :: gdepem = .true.


!       1. Initialise tracer/family chemistry tendencies
!          ---------- ------------- --------- ----------

        DO jtr = 1, jpctr
          DO jl = 1, kl
            fdot(jl,jtr) = 0.0
          ENDDO
        ENDDO

!       2. Calculate production & loss terms of individual species
!          --------- ---------- - ---- ----- -- ---------- -------

! DEPENDS ON: asad_prls
        CALL asad_prls( kl, nf, nlf, gdepem )

!       3. Calculate tendencies
!          --------- ----------

        DO j = 1, nf
          js = nlf(j)
          ifam = moffam(js)
          itr  = madvtr(js)
          gfam = ifam /= 0
          gtr  = itr  /= 0
          DO jl = 1, kl

!           3.1 Tendencies of individual species

            ydot(jl,js) = prod(jl,js) - slos(jl,js)

!           3.2 Tendencies of families 
!            (add in/out species if in family).

            if ( gfam .and. (.not. gtr .or. gtr.and.linfam(jl,itr) ) ) &
               fdot(jl,ifam) = fdot(jl,ifam) + nodd(js) * ydot(jl,js)

!           3.3 Tendencies of non-family tracers 
!           (and in/out species).

            IF ( gtr ) fdot(jl,itr) = ydot(jl,js)

          ENDDO
        ENDDO

        RETURN
        END SUBROUTINE ASAD_DIFFUN


      end Module messy_clamschem_asad_diffun

Module messy_clamschem_asad_fyself

contains

!****  asad_fyself - calculates self-reacting terms for individual species
!
!     Paul D.Brown & Glenn Carver, Centre for Atmospheric Science,
!                                  University of Cambridge.
!
!
! Purpose: Calculates self-reacting terms for individual species.
!
!
!          Called from ASAD_FTOY
!
!
!     Method
!     ------
!     The reaction table is scanned, and the self-reacting terms are
!     calculated from the appropriate rate coefficients.
!
!     Local variables
!     ---------------
!     isp            Index of self reacting species.
!
! Code description:
!   Language: FORTRAN 90
!
! ---------------------------------------------------------------------
!
        SUBROUTINE ASAD_FYSELF(n_points)

        USE messy_clamschem_asad_mod_clams, ONLY: jpnr
        USE messy_clamschem_asad_mod,       ONLY: qa, rk, nstst, nlstst, nspi

        IMPLICIT NONE

        INTEGER, INTENT(IN) :: n_points    ! No of spatial points

!       Local variables

        INTEGER :: j                       ! Loop variable
        INTEGER :: jl                      ! Loop variable
        INTEGER :: jr                      ! Loop variable
        INTEGER :: js                      ! Index
        INTEGER :: isp                     ! Index


!       1.  Initialisation
!           --------------

        DO j = 1, nstst
          js = nlstst(j)
          DO jl = 1, n_points
            qa(jl,js) = 0.0
          ENDDO
        ENDDO

!       2.  Calculate self-reacting terms
!           --------- ------------- -----

        DO jr = 1, jpnr
          isp = nspi(jr,1)
! ju_nt_20140212
!jug following line changed:
!          IF ( isp == nspi(jr,2) ) THEN
          IF ( isp == nspi(jr,2) .and. isp > 0) then
            DO jl = 1, n_points
              qa(jl,isp) = qa(jl,isp) + 2.0 * rk(jl,jr)
            ENDDO
          ENDIF
        ENDDO

        RETURN
        END SUBROUTINE ASAD_FYSELF

      End Module messy_clamschem_asad_fyself

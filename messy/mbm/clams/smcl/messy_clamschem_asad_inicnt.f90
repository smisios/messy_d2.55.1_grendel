Module messy_clamschem_asad_inicnt

contains

!
!     inicnt  - Initialises constant species.
!
!     Glenn Carver,    Centre for Atmospheric Science
!                      University of Cambridge.
!
! Purpose: Species labelled 'CF': ASAD will treat the species as a constant
!     but will call this routine so that the user may set the
!     values differently at each gridpoint for example. Currently
!     only used for setting the water vapour concentration.
!
!
!          Called from ASAD_FYINIT
!
!
!     Interface
!     On entry, the following will be set:
!              species - character name of species to set.
!                        Will be the same as listed in chch.d file
!              klen    - length of array, y.
!
!     On exit, the following must be set:
!              y       - Array of points to set for the species.
!
! Code description:
!   Language: FORTRAN 90
!
! ---------------------------------------------------------------------
!
        SUBROUTINE ASAD_INICNT( species, y, klen )

        USE messy_clams_global, ONLY: prec

        USE messy_clamschem_asad_mod,   ONLY: wp, tnd
        USE messy_clamschem_asad_dummy, ONLY: ereport

        IMPLICIT NONE

        INTEGER, INTENT(IN) :: klen      ! No of spatial points

        CHARACTER (LEN=10), INTENT(IN)  :: species  ! Species char strng

        REAL(PREC), INTENT(OUT)   :: y(klen)   ! Species concentration

!       Local variables

        INTEGER :: errcode                ! Variable passed to ereport
        INTEGER :: jl

        CHARACTER (LEN=72) :: cmessage


!       1.  Copy water into ASAD array.

        IF ( species(1:3) /= 'H2O' ) THEN
           errcode=124
           cmessage= 'Expected species H2O but got '//species

           CALL EREPORT('ASAD_INICNT',errcode,cmessage)
        ENDIF

        DO jl = 1, klen
          y(jl) = wp(jl)*tnd(jl)
        ENDDO

        RETURN
        END SUBROUTINE ASAD_INICNT


      End Module messy_clamschem_asad_inicnt

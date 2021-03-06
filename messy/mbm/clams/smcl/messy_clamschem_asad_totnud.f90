Module messy_clamschem_asad_totnud

contains

!
!     asad_totnud - calculates total number density in gridboxes 
!
!     Glenn Carver & Paul Brown, Centre for Atmospheric Science,
!                                  University of Cambridge.
!
!
! Purpose: Calculates the total number density in gridboxes
!
!          Called from ASAD_CDRIVE
!
!
!     Method
!     ------
!     The total number density at a point is given by: p=nkt
!     p-pressure , n-number density, k-boltzmann's constant,
!     t-temperature.
!
!     local variables
!     ---------------
!     zboltz     boltzmann's constant
!
! Code description:
!   Language: FORTRAN 90
!
! ---------------------------------------------------------------------
!
        SUBROUTINE ASAD_TOTNUD(n_points)

        USE messy_clamschem_asad_mod,  ONLY: tnd, p, t, pmintnd, pmin
        USE messy_clams_global,        ONLY: prec

        IMPLICIT NONE


        INTEGER, INTENT(IN) :: n_points

!       Local variables

! ju_nt_20140203
        real, parameter :: zboltz = 1.3806e-23

        INTEGER :: jl

        REAL(PREC) :: zb


        zb = zboltz*1.0e6

!       1. Total number density (1e6 converts numbers to /cm**3).
!          ----- ------ ------- ---- -------- ------- -- --------

        DO jl = 1, n_points
! ju_nt_20140212          
!     the following if statement avoids division by 0 if jpnl > ntraj
!     (jug, 3.6.98)
           IF (t(jl) > 0.0) THEN
              tnd(jl)     = p(jl) / ( zb * t(jl) )
              pmintnd(jl) = pmin * tnd(jl)
           ENDIF
        ENDDO

        RETURN
        END SUBROUTINE ASAD_TOTNUD

      End Module messy_clamschem_asad_totnud

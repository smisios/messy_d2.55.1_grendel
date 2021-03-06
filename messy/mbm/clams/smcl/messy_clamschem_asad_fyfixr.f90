Module messy_clamschem_asad_fyfixr

contains

!****  asad_fyfixr  - calculates family concentrations from fixed ratios
!
!     Paul D.Brown & Glenn Carver, Centre for Atmospheric Science,
!                                  University of Cambridge.
!
!
! Purpose: Calculates family concentrations from fixed ratios, calculated
!          during a previous call to asad_ftoy.
!
!
!          Called from ASAD_FTOY
!
!
!     Interface
!     ---------
!     Called from ftoy if the number of iterations requested equals
!     zero.
!
!     Method
!     ------
!     The family is partitioned amongst the members
!     using ratios calculated on a previous call to ftoy.
!     The concentrations are found from:
!                  Ym = Rmf*Z, Y1 = R1m*Ym, Y2 = R2m*Ym, .....
!
!     Local variables
!     ---------------
!     ifam           Index of family to which species belongs.
!     imaj           Index of major member of family to which
!                    species belongs.
!     itr            Index of model tracer to which species
!                    corresponds.
!
! Code description:
!   Language: FORTRAN 90
!
! ---------------------------------------------------------------------
!
        SUBROUTINE ASAD_FYFIXR(n_points)

        USE messy_clamschem_asad_mod,   ONLY: y, f, ratio, linfam, nlmajmin, jpif,  &
                                              moffam, madvtr, majors, ctype

        IMPLICIT NONE

        INTEGER, INTENT(IN) :: n_points

!       Local variables

        INTEGER :: istart
        INTEGER :: iend
        INTEGER :: j       ! Loop variable
        INTEGER :: jl      ! Loop variable
        INTEGER :: js      ! Index
        INTEGER :: ifam    ! Index
        INTEGER :: imaj    ! Index
        INTEGER :: itr     ! Index


!       1.  Calculate major species of family
!           --------- ----- ------- -- ------

        istart = nlmajmin(1)
        iend   = nlmajmin(2)
        DO j = istart, iend
          js   = nlmajmin(j)
          ifam = moffam(js)
          DO jl = 1, n_points
            y(jl,js) = f(jl,ifam) * ratio(jl,js)
          ENDDO
        ENDDO

!       3.  Calculate minor species of family
!           --------- ----- ------- -- ------

        istart = nlmajmin(3)
        iend   = nlmajmin(4)
        DO j = istart, iend
          js   = nlmajmin(j)
          ifam = moffam(js)
          imaj = majors(ifam)
          IF ( ctype(js) /= jpif ) THEN
            DO jl = 1, n_points
              y(jl,js) = y(jl,imaj) * ratio(jl,js)
            ENDDO
          ELSE
            itr = madvtr(js)
            DO jl = 1, n_points
              IF ( linfam(jl,itr) ) y(jl,js) =                         &
                                    y(jl,imaj) * ratio(jl,js)
            ENDDO
          ENDIF
        ENDDO

        RETURN
        END SUBROUTINE ASAD_FYFIXR


      End Module messy_clamschem_asad_fyfixr

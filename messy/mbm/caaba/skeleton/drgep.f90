!*****************************************************************************

! Calculation of direct interaction coefficients (DIC, also called
! r_AB), as described by Niemeyer et al., Combust. Flame 157, 1760–1770
! (2010), doi: 10.1016/j.combustflame.2009.12.022

! Authors:
! Kyle Niemeyer: Original code (MODULE skeletal_reduction)
! Rolf Sander:   slightly modified to eliminate dependence on chemkin

! This file is available under the MIT License (MIT)
 
! Permission is hereby granted, free of charge, to any person obtaining
! a copy of this software and associated documentation files (the
! "Software"), to deal in the Software without restriction, including
! without limitation the rights to use, copy, modify, merge, publish,
! distribute, sublicense, and/or sell copies of the Software, and to
! permit persons to whom the Software is furnished to do so, subject to
! the following conditions:
 
! The above copyright notice and this permission notice shall be
! included in all copies or substantial portions of the Software.

! THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
! EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
! MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
! IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
! CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
! TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
! SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

!*****************************************************************************

MODULE drgep

  USE messy_main_constants_mem, ONLY: DP
  IMPLICIT NONE
  REAL(DP), PARAMETER :: ZERO = 0.0_DP, ONE = 1.0_DP

CONTAINS

  ! --------------------------------------------------------------------------

  SUBROUTINE CalcDIC(nuki, rates, rkj, neighbor, n_neigh)

    IMPLICIT NONE

    ! subroutine dummy arguments:
    ! stoichiometric numbers nuki = \nu_{k,i} for species k and reaction i:
    REAL(DP), DIMENSION(:,:), INTENT(IN)  :: nuki
    ! reaction rates (qmatd):
    REAL(DP), DIMENSION(:),   INTENT(IN)  :: rates     
    ! DIC rkl = r_{k,j}, how much does production of k depend on j:
    REAL(DP), DIMENSION(:,:), INTENT(OUT) :: rkj       
    ! neighbor and n_neigh, needed for dijkstra_adj:
    INTEGER, DIMENSION(:,:),  INTENT(OUT) :: neighbor
    INTEGER, DIMENSION(:),    INTENT(OUT) :: n_neigh

    ! local automatic arrays:
    ! SIZE(nuki,1) = num_s = number of species
    ! SIZE(nuki,2) = num_r = number of reactions
    LOGICAL, DIMENSION(:,:) :: rkj_flag(SIZE(nuki,1),SIZE(nuki,1))
    REAL(DP), DIMENSION(:)  :: Pa(SIZE(nuki,1)), Ca(SIZE(nuki,1))
    INTEGER, DIMENSION(:,:) :: edges(SIZE(nuki,1)**2,2)
    ! local scalars:
    INTEGER :: num_r, num_s, nedges, i, j, k
    REAL(DP) :: pcmx, wdotki

    num_s = SIZE(nuki,1)
    num_r = SIZE(nuki,2)

    rkj(:,:) = ZERO
    Pa(:) = ZERO
    Ca(:) = ZERO
    nedges = 0
    edges(:,:) = 0
    rkj_flag(:,:) = .FALSE.
    neighbor(:,:) = 0
    n_neigh(:) = 0

    ! loop through reactions:
    DO i = 1, num_r

      ! loop through species:
      spec_loop1: DO k = 1, num_s

        IF (nuki(k,i) == 0) CYCLE spec_loop1

        wdotki = rates(i) * nuki(k,i)

        ! calculate reaction contribution to denominator of r_AB
        ! production of species k due to reaction i
        Pa(k) = Pa(k) + MAX(zero, wdotki)
        ! consumption of species k due to reaction i
        Ca(k) = Ca(k) + MAX(zero, -wdotki)

        pair_loop1: DO j = 1, num_s
          IF (k == j) CYCLE pair_loop1
          IF (nuki(j,i) == 0) CYCLE pair_loop1

          ! store edge if new
          IF ( .NOT.(rkj_flag(k,j)) ) THEN
            rkj_flag(k,j) = .TRUE.
            nedges = nedges + 1
            edges(nedges, 1) = k
            edges(nedges, 2) = j
            n_neigh(k) = n_neigh(k) + 1
            neighbor(k, n_neigh(k)) = j
          END IF

          ! add contribution to numerator of r_AB
          rkj(k,j) = rkj(k,j) + wdotki
        END DO pair_loop1
      END DO spec_loop1
    END DO

    ! now calculate r_AB by looping through number of edges
    ! (linearly proportional to num_r)
    DO i = 1, nedges
      k = edges(i, 1)
      j = edges(i, 2)
      pcmx = MAX(Pa(k), Ca(k))
      IF (pcmx > zero) THEN
        rkj(k,j) = ABS(rkj(k,j)) / pcmx
      ELSE
        rkj(k,j) = zero
      END IF
    END DO

  END SUBROUTINE CalcDIC

  ! --------------------------------------------------------------------------

END MODULE drgep

!*****************************************************************************

! -*- f90 -*-
  USE messy_main_constants_mem,   ONLY: dp
  USE messy_main_grid_def_mem_bi, ONLY: nlev

IMPLICIT NONE

  REAL(DP), DIMENSION(1) :: longitude = 0.0_dp
  REAL(DP), DIMENSION(1) :: latitude  = 0.0_dp
  REAL(DP), ALLOCATABLE, DIMENSION(:) :: level
  REAL(DP), ALLOCATABLE, DIMENSION(:) :: hyai

  ! SPECIFIC FOR QBO +++++++++++++++++++++++++++++++++++++++++++++++++++++++
  REAL(dp) :: ceta(nlev) = &
       (/ 1.0, 3.2, 5.8, 9.0, 12.7, 17.2, 22.3, 28.1, 34.9, 42.6, &
       51.4, 61.3, 72.2, 84.2, 97.4, 112.0, 128.0, 145., 164., &
       186.2, 209.6, 235.0, 262.7, 292.9, 325.8, 361.6, 400.7, &
       443.3, 489.8, 540.4, 595.4, 655.1, 719.7, 789.7, 865.6, &
       947., 1036., 1133., 1238., 1351., 1474., 1607., &
       1751., 1908., 2077., 2260., 2459., 2674., 2908., &
       3161., 3436., 3736., 4064., 4421., 4813., 5243., &
       5714., 6232., 6802., 7428., 8114., 8867., 9698., &
       10617., 11633., 12760., 14011., 15397., 16932., &
       18628., 20502., 22571., 24855., 27374., 30149., &
       33209., 36585., 40309., 44416., 48939., 53920., 59400.,&
       65421., 72035., 79208., 86583., 93068., 97760., 100621., 102004./)
  REAL(dp) :: philat_2d(1,1)   = 0._dp ! see latitude above !

CONTAINS

		! --------------------------------------------------------------------------
  SUBROUTINE main_grid_def_setup

    IMPLICIT NONE

    INTEGER :: i

    ! DEFINE GRID
    ! LEVELS
    ALLOCATE(level(nlev))
    DO i=1, nlev
       level(i) = REAL(i, DP)
    END DO
    ALLOCATE(hyai(nlev+1))
    hyai(nlev+1) = ceta(nlev) + ( ceta(nlev) - ceta(nlev-1)) / 2
    hyai(nlev)   = ( ceta(nlev) + ceta(nlev-1)) / 2
    DO i=nlev-1,2,-1
       hyai(i) = (ceta(i)+ceta(i-1)) /2.
    END DO
    hyai(1) = ceta(1) -(ceta(2)-ceta(1))/2

  END SUBROUTINE main_grid_def_setup
  ! --------------------------------------------------------------------------

  ! --------------------------------------------------------------------------
  SUBROUTINE main_grid_def_initialize

    IMPLICIT NONE

  END SUBROUTINE MAIN_GRID_DEF_INITIALIZE
  ! --------------------------------------------------------------------------
  ! --------------------------------------------------------------------------
  SUBROUTINE main_grid_def_init_memory

    IMPLICIT NONE

  END SUBROUTINE MAIN_GRID_DEF_INIT_MEMORY
  ! --------------------------------------------------------------------------
  ! --------------------------------------------------------------------------
  SUBROUTINE main_grid_def_read_restart

    IMPLICIT NONE

  END SUBROUTINE MAIN_GRID_DEF_READ_RESTART
  ! --------------------------------------------------------------------------
  ! --------------------------------------------------------------------------
  SUBROUTINE main_grid_def_global_start

    IMPLICIT NONE

  END SUBROUTINE MAIN_GRID_DEF_GLOBAL_START
  ! --------------------------------------------------------------------------
  ! --------------------------------------------------------------------------
  SUBROUTINE main_grid_def_local_start

    IMPLICIT NONE

  END SUBROUTINE main_grid_def_local_start
  ! --------------------------------------------------------------------------
  ! --------------------------------------------------------------------------
  SUBROUTINE main_grid_def_free_memory

    IMPLICIT NONE

    DEALLOCATE(level)

  END SUBROUTINE MAIN_GRID_DEF_FREE_MEMORY
  ! --------------------------------------------------------------------------

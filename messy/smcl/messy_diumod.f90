! **********************************************************************
!
! SUBMODEL CORE LAYER (SMCL) ROUTINES FOR MESSy SUBMODEL DIUMOD
!
! THIS SUBMODEL IS USED TO APPROXIMATE THE DIURNAL ABS. SIN-VARIATION 
! OF A TRACER FROM AN AVERAGE DAILY VALUE
!
! Author : Sergey Gromov, MPI-C, 2014-2015
!
! References:
!
! * SO FAR NONE
!
! **********************************************************************

!> \brief DIUMOD core module.

!> \authors Sergey Gromov, MPI-C, 2014-2015
!>   - original DIUMOD submodel code

!> \version 1.0
!>   - initial release
!>   - successful test with MESSy 2.50 (OHVAR experiment) using imported OH field and TREXP-created OH tracers

!> \todo
!>   - So far nothing


! **********************************************************************
MODULE messy_diumod
! **********************************************************************

  ! ----------- >

  USE messy_main_constants_mem, ONLY: DP

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: DP

  ! ----------- <

  CHARACTER(LEN=*), PARAMETER, PUBLIC :: modstr = 'diumod' !< sub-model name
  CHARACTER(LEN=*), PARAMETER, PUBLIC :: modver = '1.0'    !< sub-model version

  ! ### add your own public subroutines here
  ! PUBLIC :: one

  ! PRIVATE SUBROUTINES
  ! ### add your own private subroutines here
  ! PRIVATE :: two

!CONTAINS

! **********************************************************************
END MODULE messy_diumod
! **********************************************************************


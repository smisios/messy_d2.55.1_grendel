! ************************************************************************
MODULE messy_main_tools_bi
! ************************************************************************

  ! MESSy-SMIL tools (ECHAM5)
  !
  ! Author: Patrick Joeckel, MPICH, Mainz, Sep 2003

  IMPLICIT NONE
  PRIVATE

  ! MESSy
  PUBLIC :: main_tools_initialize

CONTAINS

! --------------------------------------------------------------------
  SUBROUTINE main_tools_initialize

    ! MESSy
    USE messy_main_tools, ONLY: init_convect_tables

    IMPLICIT NONE

    CALL init_convect_tables

  END SUBROUTINE main_tools_initialize
! --------------------------------------------------------------------

! ************************************************************************
END MODULE messy_main_tools_bi
! ************************************************************************

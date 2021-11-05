! -*- f90 -*-
!*****************************************************************************
 MODULE messy_main_grid_def_bi
!*****************************************************************************

   USE messy_main_grid_def_mem_bi, ONLY: nlev, kproma, jrow
   USE messy_main_constants_mem, ONLY: dp
   
   IMPLICIT NONE
   SAVE
   
   CHARACTER(LEN=*), PARAMETER :: modstr = 'grid_def'
   CHARACTER(LEN=*), PARAMETER :: modver = '1.0'
   
   REAL(dp) :: ceta(nlev), hyam(nlev), hybm(nlev)
   REAL(dp) :: gl_twomu(1), gl_sqcst(1)
   
   REAL(dp), DIMENSION(1) :: &
        philat = -60., philon = 0. 
   
   INTEGER, DIMENSION(kproma,jrow) :: &
        ilat = 1, ilon = 1

!#include "grid_def_L39.inc"
#include "grid_def_L90.inc"
!#include "grid_def_L50_Medvedev.inc"
!*****************************************************************************
END MODULE messy_main_grid_def_bi
!*****************************************************************************

! ****************************************************************************
!                Time-stamp: <2020-09-15 16:38:21 sander>
!     Author: Rolf Sander, Max-Planck Institute, Mainz, Germany, 2009-2014
!*****************************************************************************

! This program is free software; you can redistribute it and/or
! modify it under the terms of the GNU General Public License
! as published by the Free Software Foundation; either version 2
! of the License, or (at your option) any later version.
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
! You should have received a copy of the GNU General Public License
! along with this program; if not, get it from:
! http://www.gnu.org/copyleft/gpl.html

! ****************************************************************************

MODULE jvpp_mem

  USE messy_main_constants_mem, ONLY: DP
  USE messy_cmn_photol_mem, ONLY: IP_MAX
  IMPLICIT NONE

  CHARACTER(LEN=*), PARAMETER :: HLINE1 = &
    '************************************************************************'
  CHARACTER(LEN=*), PARAMETER :: HLINE2 = &
    '------------------------------------------------------------------------'
  CHARACTER(LEN=*), PARAMETER :: HLINE2b = &
    '------------------------------------------------------------------'

  ! temporary (workdir) directories:
  CHARACTER(LEN=*), PARAMETER :: DIR_176 = 'workdir_176'
  CHARACTER(LEN=*), PARAMETER :: DIR_EFF = 'workdir_eff'
  CHARACTER(LEN=*), PARAMETER :: DIR_F90 = 'workdir_f90'
  CHARACTER(LEN=*), PARAMETER :: DIR_JNL = 'workdir_jnl'

  ! input directory ("dat_m17" or "dat_lit"):
  CHARACTER(LEN=7) :: inputdir

  ! l_* logicals showing which input files are available:
  LOGICAL, DIMENSION(IP_MAX) :: l_sig     = .FALSE.
  LOGICAL, DIMENSION(IP_MAX) :: l_s2t     = .FALSE.
  LOGICAL, DIMENSION(IP_MAX) :: l_s3t     = .FALSE.
  LOGICAL, DIMENSION(IP_MAX) :: l_phi     = .FALSE.
  LOGICAL, DIMENSION(IP_MAX) :: l_tc1     = .FALSE.
  LOGICAL, DIMENSION(IP_MAX) :: l_tc2     = .FALSE.
  LOGICAL, DIMENSION(IP_MAX) :: l_tc3     = .FALSE.
  LOGICAL, DIMENSION(IP_MAX) :: l_species = .FALSE.
  LOGICAL, DIMENSION(IP_MAX) :: l_tdep    = .FALSE.

  ! IO-units:
  INTEGER, PARAMETER :: IO_LOG = 99 ! logfile
  INTEGER, PARAMETER :: IO_TMP = 98 ! temporary files

  ! IO-units for step 1:
  INTEGER, PARAMETER :: IO_SIG = 12 ! *.sig
  INTEGER, PARAMETER :: IO_176 = 13 ! *.176

  ! IO-units for step 2:
  INTEGER, PARAMETER :: IO_IN  = 21 ! input files
  INTEGER, PARAMETER :: IO_OUT = 22 ! output files
  INTEGER, PARAMETER :: IO_TDP = 23 ! temperature-dependent spectra

  ! IO-units for step 3:
  INTEGER, PARAMETER :: IO_JTC = 32 ! jnl, T-const
  INTEGER, PARAMETER :: IO_JTD = 33 ! jnl, T-dep
  INTEGER, PARAMETER :: IO_ERR = 34 ! error in %
  INTEGER, PARAMETER :: IO_REC = 35 ! recalc
  INTEGER, PARAMETER :: IO_F90 = 36 ! jval_*.f90
  INTEGER, PARAMETER :: IO_NML = 37 ! nml file
  INTEGER, PARAMETER :: IO_CAT = 38 ! cat script
  INTEGER, PARAMETER :: IO_CAL = 39 ! jval_cal subroutine
  INTEGER, PARAMETER :: IO_TEX = 40 ! LaTeX file

  INTEGER, PARAMETER :: NWAV = 176
  INTEGER, PARAMETER :: NBIN = 8
  ! wavelengths [nm] in the middle of the intervals (partially defined
  ! by Bruehl & Crutzen, Clim. Dyn. 2, 173-203, 1988, ref2635)
  ! DO jw=1,13
  !   wave_nm(jw) = 1E7/(56250.-500.*jw)
  ! ENDDO
  ! DO jw=14,45
  !   wave_nm(jw) = 1E7/(49750.-(jw-13)*500.)
  ! ENDDO
  ! DO jw=46,68
  !   wave_nm(jw) = (266.+(jw-13))
  ! ENDDO
  ! DO jw=69,71
  !   wave_nm(jw) = (320.5+2.*(jw-68))
  ! ENDDO
  ! DO jw=72,NWAV
  !   wave_nm(jw) = (325.+5.*(jw-71))
  ! ENDDO
  REAL(dp), PARAMETER :: wave_nm(NWAV) = (/ &
    179.37_dp, 181.00_dp, 182.65_dp, 184.33_dp, 186.05_dp, 187.79_dp, 189.57_dp, 191.39_dp, 193.24_dp, &
    195.12_dp, 197.04_dp, 199.00_dp, 201.01_dp, 203.05_dp, 205.13_dp, 207.25_dp, 209.42_dp, 211.64_dp, &
    213.90_dp, 216.22_dp, 218.58_dp, 220.99_dp, 223.46_dp, 225.99_dp, 228.57_dp, 231.21_dp, 233.92_dp, &
    236.69_dp, 239.52_dp, 242.42_dp, 245.40_dp, 248.45_dp, 251.57_dp, 254.78_dp, 258.06_dp, 261.44_dp, &
    264.90_dp, 268.46_dp, 272.11_dp, 275.86_dp, 279.72_dp, 283.69_dp, 287.77_dp, 291.97_dp, 296.30_dp, &
    299.00_dp, 300.00_dp, 301.00_dp, 302.00_dp, 303.00_dp, 304.00_dp, 305.00_dp, 306.00_dp, 307.00_dp, &
    308.00_dp, 309.00_dp, 310.00_dp, 311.00_dp, 312.00_dp, 313.00_dp, 314.00_dp, 315.00_dp, 316.00_dp, &
    317.00_dp, 318.00_dp, 319.00_dp, 320.00_dp, 321.00_dp, 322.50_dp, 324.50_dp, 326.50_dp, 330.00_dp, &
    335.00_dp, 340.00_dp, 345.00_dp, 350.00_dp, 355.00_dp, 360.00_dp, 365.00_dp, 370.00_dp, 375.00_dp, &
    380.00_dp, 385.00_dp, 390.00_dp, 395.00_dp, 400.00_dp, 405.00_dp, 410.00_dp, 415.00_dp, 420.00_dp, &
    425.00_dp, 430.00_dp, 435.00_dp, 440.00_dp, 445.00_dp, 450.00_dp, 455.00_dp, 460.00_dp, 465.00_dp, &
    470.00_dp, 475.00_dp, 480.00_dp, 485.00_dp, 490.00_dp, 495.00_dp, 500.00_dp, 505.00_dp, 510.00_dp, &
    515.00_dp, 520.00_dp, 525.00_dp, 530.00_dp, 535.00_dp, 540.00_dp, 545.00_dp, 550.00_dp, 555.00_dp, &
    560.00_dp, 565.00_dp, 570.00_dp, 575.00_dp, 580.00_dp, 585.00_dp, 590.00_dp, 595.00_dp, 600.00_dp, &
    605.00_dp, 610.00_dp, 615.00_dp, 620.00_dp, 625.00_dp, 630.00_dp, 635.00_dp, 640.00_dp, 645.00_dp, &
    650.00_dp, 655.00_dp, 660.00_dp, 665.00_dp, 670.00_dp, 675.00_dp, 680.00_dp, 685.00_dp, 690.00_dp, &
    695.00_dp, 700.00_dp, 705.00_dp, 710.00_dp, 715.00_dp, 720.00_dp, 725.00_dp, 730.00_dp, 735.00_dp, &
    740.00_dp, 745.00_dp, 750.00_dp, 755.00_dp, 760.00_dp, 765.00_dp, 770.00_dp, 775.00_dp, 780.00_dp, &
    785.00_dp, 790.00_dp, 795.00_dp, 800.00_dp, 805.00_dp, 810.00_dp, 815.00_dp, 820.00_dp, 825.00_dp, &
    830.00_dp, 835.00_dp, 840.00_dp, 845.00_dp, 850.00_dp /)

END MODULE jvpp_mem

! ****************************************************************************

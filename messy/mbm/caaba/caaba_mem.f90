! Time-stamp: <2019-07-25 12:31:42 sander>

! This file declares global variables for caaba

! Authors:
! Hella Riede, MPICH, Mainz, 2007: original code
! Rolf Sander, MPICH, Mainz, 2007:

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

!*****************************************************************************

MODULE caaba_mem

  USE messy_main_constants_mem, ONLY: STRLEN_VLONG, STRLEN_MEDIUM, &
                                      STRLEN_ULONG, DP

  IMPLICIT NONE

  ! PUBLIC is already default
  SAVE

  ! in SUBROUTINE caaba_init, caaba_version will be set to that of MECCA:
  CHARACTER(LEN=STRLEN_MEDIUM) :: caaba_version = ''

  ! CAABA namelist:
  ! nml SUBMODELS
  LOGICAL  :: USE_CLOUDJ          = .FALSE. ! (messy_main_switch)
  LOGICAL  :: USE_DISSOC          = .FALSE. ! (messy_main_switch)
  LOGICAL  :: USE_JVAL            = .FALSE. ! (messy_main_switch)
  LOGICAL  :: USE_MECCA           = .FALSE. ! (messy_main_switch)
  LOGICAL  :: USE_RADJIMT         = .FALSE. ! (messy_main_switch)
  LOGICAL  :: USE_READJ           = .FALSE. ! (messy_main_switch)
  LOGICAL  :: USE_SAPPHO          = .FALSE. ! (messy_main_switch)
  LOGICAL  :: USE_SEMIDEP         = .FALSE. ! (messy_main_switch)
  LOGICAL  :: USE_TRAJECT         = .FALSE. ! (messy_main_switch)
#ifdef E4CHEM
  LOGICAL  :: USE_E4CHEM          = .FALSE. ! (messy_main_switch)
#endif

  ! nml TIME
  LOGICAL  :: l_steady_state_stop = .FALSE. ! stop caaba at steady state?
  !mz_hr_20160418+ 
  LOGICAL  :: l_groundhogday      = .FALSE. ! repeat first day?
  LOGICAL  :: l_freezetime        = .FALSE. ! freeze time for simulation?
  CHARACTER(LEN=STRLEN_VLONG)  :: runtime_str     = '' ! external runtime
  CHARACTER(LEN=STRLEN_VLONG)  :: timesteplen_str = '' ! time step length string
  !mz_hr_20160418- 

  ! nml PHYS / CHEM
  LOGICAL  :: l_ff                = .FALSE. ! frostflower model run?
  LOGICAL  :: l_skipkpp           = .FALSE. ! skip KPP calculations?
  LOGICAL  :: l_skipoutput        = .FALSE. ! skip output (for benchmarking/aMC)
  LOGICAL  :: l_ignore_relhum     = .FALSE. ! if true, don't convert relhum/spechum 
                                            ! to c(H2O), instead initialize c(H2O)
  LOGICAL  :: l_relhum_wmo        = .FALSE. ! water vapor mass mixing ratios
  !mz_hr_20160503+ 
  !LOGICAL  :: l_psat_emac         = .FALSE. ! psat from EMAC look-up tables
  LOGICAL  :: l_hum_emac          = .FALSE. ! psat and cair functions as in EMAC
  LOGICAL  :: l_psat_liquid       = .FALSE. ! sat. vapor press always ov. liquid
  !mz_hr_20160503- 
  LOGICAL  :: l_RRconc            = .FALSE. ! convert RR to [conc/time]?
  LOGICAL  :: l_injectNOx         = .FALSE. ! additional NOx emissions?
  !mz_hr_20160508+ moved here because they are nml parameters
  REAL(DP) :: temp    = 293._DP    ! temperature [K]
  REAL(DP) :: press   = 101325._DP ! pressure [Pa]
  REAL(DP) :: relhum  = 0.81_DP    ! rel. hum. [0-1]
  REAL(DP) :: spechum = -1._DP     ! (kg water vapor (H2O))/(kg moist air)
  REAL(DP) :: zmbl    = 1000._DP   ! boundary layer height [m]
  !mz_hr_20160508- 
  LOGICAL  :: l_skeleton          = .FALSE. ! skeletal mechanism reduction?
  REAL(DP) :: degree_lat          = 45._DP  ! default lat in degrees
  REAL(DP) :: degree_lon          =  0._DP  ! default lon in degrees (0:360)
  REAL(DP) :: model_start_day     = 80._DP
  REAL(DP) :: Ca_precip           =  0._DP  ! relative CaCO3 precipitation (0...1)
  REAL(DP) :: t_NOxon             = -1._DP  ! start of injection (day of run)
  REAL(DP) :: t_NOxoff            = -1._DP  ! end of injection (day of run)
  !mz_hr_20160510+ 
  CHARACTER(LEN=*), PARAMETER :: H2Og = 'H2O' ! water vapor var in init_spec
  !mz_hr_20160510- 
  CHARACTER(LEN=12) ::     init_scenario = ''
  CHARACTER(LEN=12) ::    photo_scenario = ''
  CHARACTER(LEN=12) :: emission_scenario = ''
  CHARACTER(LEN=12) ::   drydep_scenario = ''
  
  ! nml PHOTOLYSIS
  INTEGER  :: input_readj_index = 1
  !mz_hr_20130705+ 
  REAL(DP) :: efact = 1._DP ! photolysis rate enhancement factor in SAPPHO
  !mz_hr_20130705- 
  CHARACTER(LEN=STRLEN_ULONG)  :: init_spec       = ''
  CHARACTER(LEN=STRLEN_ULONG)  :: input_readj     = ''
  CHARACTER(LEN=STRLEN_MEDIUM) :: photrat_channel = ''

  ! nml TRAJECT
  REAL(DP) :: runlast = -1._DP  ! run last x days of trajectory
  CHARACTER(LEN=STRLEN_ULONG)  :: input_physc     = ''
  CHARACTER(LEN=STRLEN_ULONG)  :: input_jval      = ''
  CHARACTER(LEN=STRLEN_VLONG)  :: vlat            = 'LAT'   ! external latitude variable
  CHARACTER(LEN=STRLEN_VLONG)  :: vlon            = 'LON'   ! external longitude variable
  CHARACTER(LEN=STRLEN_VLONG)  :: vpress          = 'PRESS' ! external pressure variable
  CHARACTER(LEN=STRLEN_VLONG)  :: vtemp           = 'TEMP'  ! external temperature variable
  CHARACTER(LEN=STRLEN_VLONG)  :: vrelhum         = ''      ! external relative humidity variable
  CHARACTER(LEN=STRLEN_VLONG)  :: vspechum        = ''      ! external specific humidity variable
  CHARACTER(LEN=STRLEN_VLONG)  :: vtime           = 'TIME'  ! external time variable

  ! op_ff_20170407+
  ! nml OUTPUT
  INTEGER  :: output_step_freq = 1
  INTEGER  :: output_sync_freq = 1
  ! ouput counter
  INTEGER  :: timestep = 0
  ! op_ff_20170407-

  ! CAABA internal:
  INTEGER  :: istatus(1:20), ierrf(1) ! KPP status and error flag
  ! TIME
  LOGICAL  :: l_runtime_str = .FALSE. ! external runtime given?
  INTEGER  :: t0year, t0month, t0day, t0hour, t0min, t0sec ! start time vars
  REAL(DP) :: model_time, model_start, model_end ! in s
  REAL(DP) :: percent_done = 0._DP
  REAL(DP) :: time0_jul = 0._DP ! Julian day of time origin
  REAL(DP) :: firstjan_jul = 0._DP ! Julian date of 01-JAN-<t0year>
  REAL(DP) :: timesteplen ! time step length
  REAL(DP) :: tuf = 1._DP ! time unit factor
  REAL(DP) :: runtime     = -1._DP  ! in days
  ! time origin (LEN=33 is necessary for open_output_file in traject_init):
  CHARACTER(LEN=33) :: time_string = 'seconds since 2000-01-01 00:00:00'
  ! activate for debugging (should be start of Julian Day 2440000)
  !CHARACTER(LEN=STRLEN_VLONG) :: &
  !  time_string = 'seconds since 1968-05-23 12:00:00' 

  ! PHOTOLYSIS
  LOGICAL  :: l_input_jval  = .FALSE. ! external j-values?
  INTEGER  :: lyear, lmonth, lday, lhour, lmin, lsec ! local time variables
  REAL(DP) :: localtime = 0._DP
  REAL(DP) :: cossza = 1._DP ! (initial dummy value)
  REAL(DP) :: degree_sza
  REAL(DP) :: x_j_no2 ! external J_NO2
  INTEGER  :: photol_clev  ! current pressure level in photolysis model
  ! Definition of default atmosphere for photolysis comparison.
  ! Global average values are extracted with ferret from messy and
  ! cloudj_diag streams using e.g.: "list rhum[i=@ave,j=@ave,l=1]"
  INTEGER, PARAMETER :: DEFAULT_NLEV = 19
  ! vertical ozone column [mcl/cm2]
  REAL(DP), PARAMETER, DIMENSION(DEFAULT_NLEV+1) :: DEFAULT_V3 = (/ &
    3.366E+17, 1.437E+18, 4.085E+18, 5.428E+18, 6.157E+18, 6.583E+18, &
    6.860E+18, 7.070E+18, 7.227E+18, 7.343E+18, 7.436E+18, 7.523E+18, &
    7.605E+18, 7.678E+18, 7.740E+18, 7.788E+18, 7.822E+18, 7.844E+18, &
    7.857E+18, 7.862E+18 /)
  ! relative ozone, i.e. ozone mixing ratio [mol/mol]
  ! Note that although relo3 has the dimension 1:nlev+1, the value
  ! relo3(1) is not used at all here. Also, note that relo3 is
  ! _only_ used for the heating rates. For the calculation of the
  ! J-values, only v3 is used.
  REAL(DP), PARAMETER, DIMENSION(DEFAULT_NLEV+1) :: DEFAULT_RELO3 = (/ &
    7.182E-06, 8.319E-06, 4.172E-06, 2.041E-06, 9.525E-07, 4.334E-07, &
    2.571E-07, 1.514E-07, 9.760E-08, 5.775E-08, 5.064E-08, 4.394E-08, &
    3.980E-08, 3.636E-08, 3.209E-08, 2.807E-08, 2.479E-08, 2.242E-08, &
    2.105E-08, 2.065E-08 /)
  ! pressure [Pa]
  REAL(DP), PARAMETER, DIMENSION(DEFAULT_NLEV)   :: DEFAULT_JPRESS = (/ &
    1000., 3000., 5040., 7339., 10248., 14053., 18935., 24966., 32107., &
    40212., 49027., 58204., 67317., 75897., 83472., 89631., 94099.,     &
    96838., 98169. /)
  ! relative humidity [%]
  REAL(DP), PARAMETER, DIMENSION(DEFAULT_NLEV)   :: DEFAULT_JRHUM = (/ &
    0.23, 1.57, 3.52, 11.73, 24.55, 25.31, 27.45, 36.46, 44.52, 46.27,  &
    46.48, 49.18, 51.73, 57.95, 72.82, 80.71, 81.66, 77.65, 76.18 /)
  ! temperature [K]
  REAL(DP), PARAMETER, DIMENSION(DEFAULT_NLEV)   :: DEFAULT_JTEMP = (/ &
    230.6, 218.2, 211.7, 207.0, 205.6, 210.9, 218.1, 225.8, 235.7, 246.6, &
    256.4, 264.2, 270.6, 275.4, 278.2, 280.9, 283.2, 284.9, 285.7 /)
  ! albedo
  REAL(DP), PARAMETER :: DEFAULT_ALBEDO = 0.07
  
  ! PHYS / CHEM
  LOGICAL :: l_spechum = .TRUE. ! specific (default) or relative humidity for traject input?
  REAL(DP), DIMENSION(:), ALLOCATABLE :: c
  REAL(DP) :: cair, cair_old     ! concentration of air [mcl/cc], old value
  REAL(DP) :: zmix    = 25._DP   ! ocean mixing height [m]
                                 ! http://en.wikipedia.org/wiki/Mixed_layer

END MODULE caaba_mem

!*****************************************************************************

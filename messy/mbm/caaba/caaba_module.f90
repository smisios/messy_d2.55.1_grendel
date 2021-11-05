! Time-stamp: <2019-07-25 13:12:21 sander>

! Author:
! Rolf Sander,   MPICH, Mainz, 2003-2019
! Hella Riede,   MPICH, Mainz, 2007
! Sergey Gromov, MPICH, Mainz, 2009-2019

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

MODULE caaba_module

  USE messy_main_constants_mem,   ONLY: DP, HLINE1, HLINE2

  IMPLICIT NONE

  CHARACTER(LEN=*), PARAMETER :: modstr = 'caaba'  ! name of module

  INTEGER :: ncid_messy

CONTAINS

  !***************************************************************************

  !> \brief Read CTRL namelist
  !> \details Read coupling namelist (based on \c dradon_read_nml_cpl by
  !> P. Joeckel)
  SUBROUTINE caaba_read_nml(status, iou)

    USE messy_main_tools, ONLY: read_nml_open, read_nml_check, read_nml_close, &
                                strcrack, ucase, &
                                init_convect_tables, & ! mz_hr_20100704
                                cair_q, spec2relhum, rel2spechum ! mz_hr_20160507
    USE messy_main_constants_mem, ONLY: R_gas, N_A, OneDay, &
                                        STRLEN_SHORT, STRLEN_MEDIUM
    USE caaba_mem, ONLY: &
      cair, cair_old, runtime, timesteplen, &
      l_input_jval, l_runtime_str, l_spechum, &
      ! namelist:
      USE_CLOUDJ, USE_DISSOC, USE_JVAL, USE_MECCA,     & !< MESSy submodels
      USE_READJ, USE_SAPPHO, USE_SEMIDEP, USE_TRAJECT, & !< MESSy submodels
#ifdef E4CHEM
      USE_E4CHEM,                                    & !< MESSy submodels
#endif
      USE_RADJIMT,                                   & !< MESSy submodels
      init_scenario, photo_scenario,                 & !< scenarios
      emission_scenario, drydep_scenario,            & !< ...
      temp, press, relhum, spechum, zmbl,            & !< CAABA (meteorology)
      l_ignore_relhum, l_hum_emac, l_relhum_wmo,     & !< CAABA (meteorology) !mz_hr_20160428
      l_psat_liquid,                                 & !< CAABA (meteorology) !mz_hr_20160428
      degree_lat, degree_lon, l_ff,                  & !< CAABA (location)
      model_start_day, runtime_str, timesteplen_str, & !< CAABA (time)
      l_skipoutput,                                  & !< CAABA (output)
      Ca_precip, init_spec,                          & !< MECCA-specific
      photrat_channel, l_skipkpp,                    & !< MECCA-specific
      l_steady_state_stop,                           & !< MECCA-specific
      l_groundhogday, l_freezetime,                  & !< MECCA-specific !mz_hr_20160418
      l_RRconc, l_skeleton,                          & !< MECCA-specific
      input_readj, input_readj_index,                & !< READJ-specific
      efact,                                         & !< SAPPHO-specific !mz_hr_20130705
      l_injectNOx, t_NOxon, t_NOxoff,                & !< SEMIDEP-specific
      runlast, input_physc, input_jval,              & !< TRAJECT-specific
      vlon, vlat, vpress, vtemp,                     & !< TRAJECT-specific
      vrelhum, vspechum, vtime,                      & !< TRAJECT-specific
      output_step_freq, output_sync_freq               !< output managing

    IMPLICIT NONE
    INTRINSIC :: ADJUSTL, TRIM

    ! I/O:
    INTEGER, INTENT(OUT) :: status     !< error status
    INTEGER, INTENT(IN)  :: iou        !< I/O unit
    ! local:
    CHARACTER(LEN=*), PARAMETER :: substr = 'caaba_read_nml'
    INTEGER, PARAMETER :: MAX_SCENARIOS = 21
    CHARACTER(LEN=12), PARAMETER, DIMENSION(MAX_SCENARIOS) :: &
      list_of_scenarios = (/ &
      '            ', 'FF_ANTARCTIC', 'FF_ARCTIC   ', 'FREE_TROP   ', &
      'HOOVER      ', 'LAB         ', 'LAB_C15     ', 'MBL         ', &
      'MOM         ', 'OOMPH       ', 'STRATO      ', 'MTCHEM      ', &
      'ISO         ', 'TAG         ', 'CUMULUS     ', 'TROPOPAUSE  ', &
      'LOW_STRATO  ', 'MID_STRATO  ', 'HIGH_STRATO ', 'ZEROAIR     ', &
      'VOLCANO     ' /)
    LOGICAL :: l_init_scenario_ok     = .FALSE.
    LOGICAL :: l_photo_scenario_ok    = .FALSE.
    LOGICAL :: l_emission_scenario_ok = .FALSE.
    LOGICAL :: l_drydep_scenario_ok   = .FALSE.
    LOGICAL :: lex   !< file exists?
    INTEGER :: fstat !< file status
    INTEGER :: i
    INTEGER                        :: nosub
    CHARACTER(LEN=STRLEN_SHORT)    :: tsunit    = '' !< time step length unit
    CHARACTER(LEN=STRLEN_SHORT)    :: rtunit    = '' !< runtime_str unit
    CHARACTER(LEN=STRLEN_MEDIUM), DIMENSION(:), POINTER :: field => NULL()

    NAMELIST /CAABA/ &
      USE_CLOUDJ, USE_DISSOC, USE_JVAL, USE_MECCA,     & !< MESSy submodels
      USE_READJ, USE_SAPPHO, USE_SEMIDEP, USE_TRAJECT, & !< MESSy submodels
#ifdef E4CHEM
      USE_E4CHEM,                                    & !< MESSy submodels
#endif
      USE_RADJIMT,                                   & !< MESSy submodels
      init_scenario, photo_scenario,                 & !< scenarios
      emission_scenario, drydep_scenario,            & !< ...
      temp, press, relhum, spechum, zmbl,            & !< CAABA (meteorology)
      l_ignore_relhum, l_hum_emac, l_relhum_wmo,     & !< CAABA (meteorology)
      l_psat_liquid,                                 & !< CAABA (meteorology)
      degree_lat, degree_lon, l_ff,                  & !< CAABA (location)
      model_start_day, runtime_str, timesteplen_str, & !< CAABA (time)
      l_skipoutput,                                  & !< CAABA (output)
      Ca_precip, init_spec,                          & !< MECCA-specific
      photrat_channel, l_skipkpp,                    & !< MECCA-specific
      l_steady_state_stop,                           & !< MECCA-specific
      l_groundhogday, l_freezetime,                  & !< MECCA-specific !mz_hr_20160418
      l_RRconc, l_skeleton,                          & !< MECCA-specific
      input_readj, input_readj_index,                & !< READJ-specific
      efact,                                         & !< SAPPHO-specific
      l_injectNOx, t_NOxon, t_NOxoff,                & !< SEMIDEP-specific
      runlast, input_physc, input_jval,              & !< TRAJECT-specific
      vlon, vlat, vpress, vtemp,                     & !< TRAJECT-specific
      vrelhum, vspechum, vtime,                      & !< TRAJECT-specific
      output_step_freq, output_sync_freq               !< output managing

    status = 1
    CALL read_nml_open(lex, substr, iou, 'CAABA', modstr)
    IF (.NOT.lex) RETURN    !< <modstr>.nml does not exist
    READ(iou, NML=CAABA, IOSTAT=fstat)
    CALL read_nml_check(fstat, substr, iou, 'CAABA', modstr)
    IF (fstat /= 0) RETURN  !< error while reading namelist
    CALL read_nml_close(substr, iou, modstr)
    status = 0  !< no error

    WRITE(*,*)
    WRITE(*,*) HLINE2
    WRITE(*,*) 'Selected MESSy submodels:'
    IF (USE_CLOUDJ)  WRITE(*,*) '  CLOUDJ'
    IF (USE_DISSOC)  WRITE(*,*) '  DISSOC'
    IF (USE_JVAL)    WRITE(*,*) '  JVAL'
    IF (USE_MECCA)   WRITE(*,*) '  MECCA'
    IF (USE_READJ)   WRITE(*,*) '  READJ'
    IF (USE_SAPPHO)  WRITE(*,*) '  SAPPHO'
    IF (USE_SEMIDEP) WRITE(*,*) '  SEMIDEP'
    IF (USE_TRAJECT) WRITE(*,*) '  TRAJECT'
#ifdef E4CHEM
    IF (USE_E4CHEM)  WRITE(*,*) '  E4CHEM'
#endif
    IF (USE_RADJIMT) WRITE(*,*) '  RADJIMT'
    WRITE(*,*) HLINE2

    ! scenarios:
    DO i=1, MAX_SCENARIOS
      IF (TRIM(list_of_scenarios(i))==TRIM(init_scenario)) &
        l_init_scenario_ok = .TRUE.
      IF (TRIM(list_of_scenarios(i))==TRIM(photo_scenario)) &
        l_photo_scenario_ok = .TRUE.
      IF (TRIM(list_of_scenarios(i))==TRIM(emission_scenario)) &
        l_emission_scenario_ok = .TRUE.
      IF (TRIM(list_of_scenarios(i))==TRIM(drydep_scenario)) &
        l_drydep_scenario_ok = .TRUE.
    ENDDO
    WRITE(*,*) 'Selected scenarios:'
    IF (l_init_scenario_ok) THEN
      WRITE(*,*) '  Init:       ', TRIM(init_scenario)
    ELSE
      WRITE(*,*) 'ERROR: unknown init scenario ', TRIM(init_scenario)
      STOP 1
    ENDIF
    IF (l_photo_scenario_ok) THEN
      WRITE(*,*) '  Photo:      ', TRIM(photo_scenario)
    ELSE
      WRITE(*,*) 'ERROR: unknown photo scenario ', TRIM(photo_scenario)
      STOP 1
    ENDIF
    IF (l_emission_scenario_ok) THEN
      WRITE(*,*) '  Emission:   ', TRIM(emission_scenario)
    ELSE
      WRITE(*,*) 'ERROR: unknown emission scenario ', TRIM(emission_scenario)
      STOP 1
    ENDIF
    IF (l_drydep_scenario_ok) THEN
      WRITE(*,*) '  Deposition: ', TRIM(drydep_scenario)
    ELSE
      WRITE(*,*) 'ERROR: unknown deposition scenario ', TRIM(drydep_scenario)
      STOP 1
    ENDIF
    WRITE(*,*) HLINE2

    WRITE(*,*) 'HUMIDITY'
    IF (l_ignore_relhum) THEN
      WRITE(*,*) 'l_ignore_relhum=TRUE'
      WRITE(*,*) '  Relative humidity, specific humidity, and'
      WRITE(*,*) '  cair are calculated from c(H2O).'
      WRITE(*,*) '  Alternative: l_ignore_relhum=FALSE to calculate c(H2O)'
      WRITE(*,*) '    from specific humidity (relative humidity).'
    ELSE
      WRITE(*,*) 'l_ignore_relhum=FALSE'
      WRITE(*,*) '  c(H2O) is calculated from specific humidity (relative humidity)'
      WRITE(*,*) '  Alternative: l_ignore_relhum=TRUE to calculate specific'
      WRITE(*,*) '    humidity, relative humidity, and cair from c(H2O).'
      IF (l_relhum_wmo) THEN
        WRITE(*,*) 'l_relhum_wmo=TRUE'
        WRITE(*,*) '  WMO definition of relative humidity:'
        WRITE(*,*) '  relhum_WMO = mass mixing ratio H2O / saturation'// &
          ' mass mixing ratio H2O (in dry air)'
        WRITE(*,*) '  always over liquid surface'
      ELSE
        WRITE(*,*) 'l_relhum_wmo=FALSE'
        WRITE(*,*) '  standard definition of relative humidity: '
        WRITE(*,*) '  relhum = p(H2O) / psat(H2O) (in humid air)'
      ENDIF
    ENDIF

    IF (l_hum_emac) THEN
      WRITE(*,*) 'l_hum_emac=TRUE'
      WRITE(*,*) '  psat(H2O) and cair calculated consistent with EMAC'
      CALL init_convect_tables
    ELSE
      WRITE(*,*) 'l_hum_emac=FALSE (default)'
      WRITE(*,*) '  psat(H2O) calculated in psat_mk'
      WRITE(*,*) '  Alternative: set l_hum_emac=T in caaba.nml'
      WRITE(*,*) '    to calculate psat and cair consistent with EMAC'
    ENDIF

    ! process request for saturation vapor pressure over liquid:
    IF (l_psat_liquid) THEN
      WRITE(*,*) 'l_psat_liquid=TRUE'
      WRITE(*,*) '  Calculating saturation water vapor pressure always'//&
                 '  over liquid surface'
    ELSE
      IF (l_relhum_wmo) THEN
        WRITE(*,*) 'l_psat_liquid=TRUE because l_relhum_wmo=TRUE'
        WRITE(*,*) '  WMO defines relative humidity over liquid surface'
        WRITE(*,*) '  Calculating saturation water vapor pressure always'//&
                   '  over liquid surface'
        l_psat_liquid=.TRUE.
      ELSE
        WRITE(*,*) 'l_psat_liquid=FALSE'
        WRITE(*,*) '  Calculating saturation water vapor pressure over'
        WRITE(*,*) '  liquid or ice, depending on temperature'
      ENDIF
    ENDIF

    IF (.NOT. USE_TRAJECT) THEN
      WRITE(*,*) HLINE2
      WRITE(*,'(A,F10.1,A)') ' latitude        = ', degree_lat, ' degree'
      WRITE(*,'(A,F10.1,A)') ' longitude       = ', degree_lon, ' degree'
      WRITE(*,'(A,F10.1,A)') ' T               = ', temp,       ' K'
      WRITE(*,'(A,F10.1,A)') ' p               = ', press,      ' Pa'
      WRITE(*,'(A,F10.1,A)') ' zmbl            = ', zmbl,       ' m'
      ! convert relhum <-> spechum as needed:
      IF (spechum >= 0._DP) THEN
        !> convert spechum given in caaba.nml to relhum:
        relhum = spec2relhum(status, spechum, temp, press, &
          l_hum_emac, l_psat_liquid, l_relhum_wmo)
        WRITE(*,'(A,ES10.3,A)') ' spechum         = ', spechum, ' kg/kg'
        WRITE(*,'(A,F10.1,A)')  ' relhum          = ', 100.*relhum, &
          ' %, calculated from spechum given in caaba.nml'
      ELSE
        !> convert relhum given in caaba.nml to spechum:
        l_spechum = .FALSE.
        IF (l_ignore_relhum) THEN
          ! mz_hh_20190502+
          ! do not call rel2spechum() if l_ignore_relhum
          ! because it may result in spechum<0 at low p
          ! mz_hh_20190502-
          PRINT *, 'skipping rel2spechum because l_ignore_relhum=TRUE'
          status = 0
        ELSE
          spechum = rel2spechum(status, relhum, temp, press, &
            l_hum_emac, l_psat_liquid, l_relhum_wmo)
          WRITE(*,'(A,F10.1,A)')  ' relhum          = ', 100.*relhum, ' %'
          WRITE(*,'(A,ES10.3,A)') ' spechum         = ', spechum, &
            ' kg/kg, calculated from relhum given in caaba.nml'
        ENDIF
      ENDIF
    ENDIF

    !> cair = c(air) in [mcl/cc], conc = cair * mixing ratio
    ! Calculate cair as function of spechum. If l_ignore_relhum=TRUE, cair
    ! will be updated in x0, after c(H2O) is initialized.

    ! cair as function of spechum
    cair = cair_q(spechum, temp, press, l_hum_emac)
    cair_old = cair

    IF (.NOT. USE_TRAJECT) THEN
      WRITE(*,'(A,ES10.3,A)') ' c(air)          = ', cair,           ' mcl/cm3'
      WRITE(*,*) HLINE2
    ENDIF

    ! plausibility checks:
    IF (USE_TRAJECT.AND.(TRIM(input_physc)=='')) THEN
      WRITE(*,*) HLINE1
      PRINT *, 'ERROR: Submodel TRAJECT requires input file'
      PRINT *, '  in variable input_physc of CAABA namelist'
      WRITE(*,*) HLINE1
      STOP 1
    ENDIF

    IF (.NOT.USE_TRAJECT.AND.(TRIM(input_physc)/='')) THEN
      WRITE(*,*) HLINE1
      PRINT *, 'ERROR: CAABA namelist parameter input_physc can only be used'
      PRINT *, '  with TRAJECT submodel (USE_TRAJECT = T in CAABA namelist)'
      WRITE(*,*) HLINE1
      STOP 1
    ENDIF

    IF (USE_TRAJECT.AND.USE_SEMIDEP) THEN
      PRINT *, 'ERROR: Surface emissions and depositions from SEMIDEP make no'
      PRINT *, '  sense for trajectories from TRAJECT with no contact to the'
      PRINT *, '  surface. However, if you know what you are doing, adjust the'
      PRINT *, '  emissions and depositions as needed and comment out this error'
      PRINT *, '  message in caaba.f90.'
      STOP 1
    ENDIF

    ! set runtime
    IF (TRIM(runtime_str) /= '') THEN
      ! crack string into value and unit
      CALL strcrack(TRIM(runtime_str), " ", field, nosub)

      rtunit        = TRIM(ADJUSTL(field(2)(1:7)))
      runtime_str   = TRIM(field(1))
      READ(runtime_str, *) runtime

      CALL ucase(rtunit)
      SELECT CASE (rtunit)
        CASE ('SECONDS', 'SECOND')
          runtime = runtime/OneDay
        CASE ('MINUTES', 'MINUTE')
          runtime = runtime/1440._DP
        CASE ('HOURS', 'HOUR')
          runtime = runtime/24._DP
        CASE ('DAYS', 'DAY')
          ! already in days...
        CASE DEFAULT
          WRITE(*,*) 'Error: unknown unit for runtime in caaba_read_nml:'
          WRITE(*,*) rtunit
          STOP 1
      END SELECT
      l_runtime_str = .TRUE.
    ELSE
      runtime = 8._DP ! days
    ENDIF

    ! set time step if specified in namelist
    IF (TRIM(timesteplen_str) /= '') THEN
      ! crack string into value and unit
      CALL strcrack(TRIM(timesteplen_str), " ", field, nosub)

      tsunit          = TRIM(ADJUSTL(field(2)(1:7)))
      timesteplen_str = TRIM(field(1))
      READ(timesteplen_str, *) timesteplen

      CALL ucase(tsunit)
      SELECT CASE (tsunit)
        CASE ('HOURS', 'HOUR')
          timesteplen = timesteplen * 3600.0_DP
        CASE ('MINUTES', 'MINUTE')
          timesteplen = timesteplen * 60.0_DP
        CASE ('SECONDS', 'SECOND')
          ! already in seconds...
        CASE DEFAULT
          WRITE(*,*) 'ERROR: unknown unit for timesteplen in caaba_read_nml:'
          WRITE(*,*) tsunit
          STOP 1
      END SELECT

    ELSE
      timesteplen = 20._DP * 60._DP    ! seconds
    ENDIF

    WRITE(*,'(A,F10.1,A)') ' timesteplen         = ', timesteplen,    ' s'
    IF (USE_TRAJECT) THEN
      WRITE(*,*)             ' model_start_day     = (calculated later)'
    ELSE
      WRITE(*,'(A,F10.1)')   ' model_start_day     = ', model_start_day
    ENDIF
    WRITE(*,*) 'l_steady_state_stop = ', l_steady_state_stop
    IF (l_steady_state_stop) THEN
      WRITE(*,*) 'runtime             =  until steady state'
    ELSE
      IF ((.NOT.USE_TRAJECT).OR.(TRIM(runtime_str) /= '')) THEN
        WRITE(*,'(A,F10.1,A)') ' runtime             = ', runtime,        ' days'
      ENDIF
    ENDIF
    IF (l_groundhogday) THEN
      WRITE(*,*) 'l_groundhogday      =  T, repeating first day'
    ELSE
      WRITE(*,*) 'l_groundhogday      = ', l_groundhogday
    ENDIF
    IF (l_freezetime) THEN
      WRITE(*,'(A,F10.1,A)') ' l_freezetime        =  T, freezing time'
    ELSE
      WRITE(*,*) 'l_freezetime        = ', l_freezetime
    ENDIF
    ! op_ff_20170407+
    IF (output_step_freq <= 0) THEN
       WRITE(*,*) 'ERROR: Output step length must be >= 1'
       STOP 1
    END IF
    IF (output_step_freq > 1) THEN
       WRITE(*,*) 'Output step length   = ', output_step_freq*timesteplen, ' s'
       WRITE(*,*) '(every ', output_step_freq, ' timesteps)'
    ELSE
       WRITE(*,*) 'Output step length   = ', timesteplen, ' s (eq. timesteplen)'
    END IF
    IF (output_sync_freq <= 0) THEN
       WRITE(*,*) 'ERROR: Output sync freq must be >= 1'
       STOP 1
    END IF
    IF (output_sync_freq > 1) THEN
       WRITE(*,*) 'Output file sync. frequency every ', output_sync_freq, ' output steps'
    ELSE
       WRITE(*,*) 'Output file sync. frequency every output step'
    END IF
    ! op_ff_20170407-
    WRITE(*,*) HLINE2

    WRITE(*,*) 'l_skipkpp           = ', l_skipkpp
    WRITE(*,*) 'l_skipoutput        = ', l_skipoutput
    WRITE(*,*) 'l_ff                = ', l_ff
    WRITE(*,*) 'Ca_precip           = ', Ca_precip, '(CaCO3 precipitation)'
    WRITE(*,*) 'l_RRconc            = ', l_RRconc
    WRITE(*,*) 'l_skeleton          = ', l_skeleton
    WRITE(*,*) HLINE2

    WRITE(*,*) 'CHEMICAL TRACER INITIALIZATION:'
    IF (TRIM(init_spec)/='') THEN
      WRITE(*,*) TRIM(init_spec)
      WRITE(*,*) 'Only species belonging to the chosen mechanism are'
      WRITE(*,*) '  initialized, see mecca.spc,'
      WRITE(*,*) '  messy_mecca_kpp_parameters.f90, or your meccanism file.'
      WRITE(*,*) 'Species in the mechanism not given in the external'
      WRITE(*,*) '  file are initialized in the chosen messy_mecca_box.f90'
      WRITE(*,*) '  init_scenario or with the default in x0_simple.'
    ELSE
      WRITE(*,*) 'init_spec = empty'
      WRITE(*,*) 'No external input for chemical tracer initialization.'
      WRITE(*,*) 'Species are initialized in the chosen messy_mecca_box.f90'
      WRITE(*,*) '  init_scenario or with the default in x0_simple.'
    ENDIF
    WRITE(*,*) HLINE2

    !!! PHOTOLYSIS
    IF (USE_READJ) THEN
      IF (TRIM(input_readj)/='') THEN
        WRITE(*,*) 'netcdf input for READJ:'
        WRITE(*,*) TRIM(input_readj)
        WRITE(*,*) 'index =', input_readj_index
      ELSE
        PRINT *, 'ERROR: ', 'input_readj is empty'
        STOP 1
      ENDIF
      WRITE(*,*) HLINE2
    ENDIF

    IF ((TRIM(photrat_channel)=='jval') .OR. &
      (TRIM(photrat_channel)=='dissoc') .OR. &
      (TRIM(photrat_channel)=='cloudj') .OR. &
      (TRIM(photrat_channel)=='radjimt') .OR. &
      (TRIM(photrat_channel)=='readj') .OR. &
      (TRIM(photrat_channel)=='sappho')) THEN
      WRITE(*,*) 'PHOTOLYSIS RATE COEFFICIENTS ARE TAKEN FROM'
      WRITE(*,*) TRIM(photrat_channel)
    ELSE
      PRINT *, 'ERROR: ', TRIM(photrat_channel), &
        'is not a valid photolysis submodel'
      STOP 1
    ENDIF
    IF ((TRIM(photrat_channel)=='cloudj').AND.(.NOT.USE_CLOUDJ)) THEN
      PRINT *, 'ERROR: photrat_channel=cloudj but USE_CLOUDJ=F'
      STOP 1
    ENDIF
    IF ((TRIM(photrat_channel)=='dissoc').AND.(.NOT.USE_DISSOC)) THEN
      PRINT *, 'ERROR: photrat_channel=dissoc but USE_DISSOC=F'
      STOP 1
    ENDIF
    IF ((TRIM(photrat_channel)=='jval').AND.(.NOT.USE_JVAL)) THEN
      PRINT *, 'ERROR: photrat_channel=jval but USE_JVAL=F'
      STOP 1
    ENDIF
    IF ((TRIM(photrat_channel)=='radjimt').AND.(.NOT.USE_RADJIMT)) THEN
      PRINT *, 'ERROR: photrat_channel=radjimt but USE_RADJIMT=F'
      STOP 1
    ENDIF
    IF ((TRIM(photrat_channel)=='readj').AND.(.NOT.USE_READJ)) THEN
      PRINT *, 'ERROR: photrat_channel=readj but USE_READJ=F'
      STOP 1
    ENDIF
    IF ((TRIM(photrat_channel)=='sappho').AND.(.NOT.USE_SAPPHO)) THEN
      PRINT *, 'ERROR: photrat_channel=sappho but USE_SAPPHO=F'
      STOP 1
    ENDIF
    IF (TRIM(photrat_channel)=='sappho') THEN
      WRITE (*,*) '  enhancement factor = ', efact
    ENDIF

    ! TRAJECT mode: note that some quantities are updated in traject_init
    !   as trajectory data not read in at this point
    IF (USE_TRAJECT) THEN
      WRITE(*,*) HLINE2
      WRITE(*,*) 'NETCDF INPUT FOR PHYSICAL DATA:'
      WRITE(*,*) TRIM(input_physc)
      IF (TRIM(input_jval)/='') THEN
        l_input_jval = .TRUE.
        WRITE(*,*) 'NETCDF INPUT FOR PHOTOLYSIS RATE COEFFICIENTS:'
        WRITE(*,*) '(AT THE MOMENT ONLY J_NO2)'
        WRITE(*,*) TRIM(input_jval)
      ENDIF
      IF (TRIM(vrelhum)/='' .AND. TRIM(vspechum)/='') THEN
        WRITE(*,*) 'ERROR: Submodel TRAJECT allows only ONE active'// &
         ' external humidity variable'
        STOP 1
      ELSEIF (TRIM(vspechum) /= '') THEN
        WRITE(*,*) 'WATER VAPOUR CONTENT GIVEN AS SPECIFIC HUMIDITY'
      ELSEIF (TRIM(vrelhum) /= '') THEN
        l_spechum = .FALSE.
        WRITE(*,*) 'WATER VAPOUR CONTENT GIVEN AS RELATIVE HUMIDITY'
      ELSE ! default is specific humidity with name 'SPECHUM':
        WRITE(*,*) "ASSUMED THAT WATER VAPOUR CONTENT GIVEN AS SPECIFIC"// &
                   " HUMIDITY ('SPECHUM')"
        vspechum = "SPECHUM"
      ENDIF
    ELSE
      IF (TRIM(input_jval)/='') THEN
        WRITE(*,*) 'ERROR: CURRENTLY, INPUT_JVAL CAN ONLY BE USED FOR'
        WRITE(*,*) 'TRAJECTORY MODEL RUNS, I.E. IF USE_TRAJECT=T'
        STOP 1
      ENDIF
    ENDIF
    WRITE(*,*) HLINE2

    IF (ASSOCIATED(field)) THEN
      DEALLOCATE(field)
      NULLIFY(field)
    ENDIF

  END SUBROUTINE caaba_read_nml

  !***************************************************************************

  SUBROUTINE time_control

    USE messy_main_constants_mem, ONLY: OneDay
    USE messy_main_timer,         ONLY: utc2lt, julian2gregor, &
                                        calc_sza
    USE caaba_mem,                ONLY: cossza, model_start, model_time, &
                                        degree_lon, degree_lat, localtime, &
                                        time0_jul, lyear, &
                                        lmonth, lday, lhour, lmin, lsec, &
                                        l_groundhogday, &
                                        l_freezetime, time_string

    IMPLICIT NONE
    INTRINSIC :: FLOOR

    INTEGER  :: status
    REAL(DP) :: time_sza, localtime_jul !> time for SZA calculation

    !> default: current model time used for SZA calculation
    time_sza = model_time

    IF (l_groundhogday) THEN
      ! repetition of the first day
      time_sza = (FLOOR(model_start/OneDay) + MODULO(model_time/OneDay,1._DP)) &
        *OneDay
    ENDIF

    ! no automatic freezing in time for l_steady_state_stop
    ! use namelist boolean l_freezetime
    ! l_freezetime after l_groundhog to overrule if necessary
    !> freeze time at model start?
    IF (l_freezetime) THEN
      time_sza = model_start
    ENDIF

    CALL calc_sza(status, cossza, time_sza, time_string, degree_lon,&
                  degree_lat)
    IF (status /= 0) THEN
      WRITE(*,*) 'ERROR in time_control: calc_sza'
      STOP 1
    ENDIF

    ! current local time
    localtime = utc2lt(status, model_time, degree_lon)
    IF (status /= 0) THEN
      WRITE(*,*) 'ERROR in time_control: utc2lt'
      STOP 1
    ENDIF

    ! local time as Julian day
    localtime_jul = time0_jul + localtime/OneDay

    ! determine local time variables year, month, ...
    CALL julian2gregor(localtime_jul, lyear, lmonth, lday, lhour, lmin, lsec)
    
  END SUBROUTINE time_control

  !***************************************************************************

  SUBROUTINE caaba_result

    USE caaba_mem,                ONLY: model_time, model_start, &
                                        model_end, USE_TRAJECT, &
                                        temp, press, cossza, &
                                        percent_done, degree_sza, &
                                        l_steady_state_stop
                                        !mz_hr_20130422 percent_done, degree_sza
    USE caaba_io,                 ONLY: write_output_file
    USE messy_main_constants_mem, ONLY: OneDay, PI, FLAGGED_BAD

    ! LOCAL
    CHARACTER(LEN=10) :: outunit

    degree_sza   = ACOS(cossza)*180./PI
    IF (l_steady_state_stop) THEN
      percent_done = (model_time-model_start)/OneDay !define as days simulated
      outunit = ' days done'
    ELSE
      percent_done = 100. * (model_time-model_start) / (model_end-model_start)
      outunit = '% done    '
    ENDIF
    IF ( (model_time - model_start > 0.) .OR. (USE_TRAJECT) ) THEN
      WRITE(*,'(A,F9.4,A,F6.2,A,F6.2,A,A)') &
      ' day = ', model_time/OneDay, &
      '    sza = ', degree_sza, &
      '     (', percent_done, TRIM(outunit), ')'
      IF (.NOT. USE_TRAJECT) THEN
        CALL write_output_file(ncid_messy, model_time, &
          (/ press, temp, degree_sza /) )
      ENDIF
    ELSE ! 1st output, sza not calc yet
      WRITE(*,'(A,F9.4,A,F6.2,A,A)') &
      ' day = ', model_time/OneDay, &
      '                     (', percent_done, outunit, ')'
      IF (.NOT. USE_TRAJECT) THEN
        CALL write_output_file(ncid_messy, model_time, &
          (/ press, temp, FLAGGED_BAD /) )
      ENDIF
    ENDIF

  END SUBROUTINE caaba_result

  !***************************************************************************

  SUBROUTINE caaba_init

#ifndef E4CHEM
    USE messy_mecca_kpp,            ONLY: NSPEC
#else
    USE messy_mecca_kpp,            ONLY: NSPEC_mecca => NSPEC
    USE messy_e4chem,               ONLY: NSPEC_fchem => NSPEC
#endif
    USE messy_mecca,                ONLY: mecca_version => modver
    USE caaba_mem,                  ONLY: model_time, model_start,     &
                                          model_end, USE_TRAJECT,      &
                                          caaba_version,  &
#ifdef E4CHEM
                                          USE_MECCA, USE_E4CHEM,       &
#endif
                                          c, runtime, model_start_day, &
                                          caaba_version
    USE messy_main_constants_mem,   ONLY: OneDay
    USE caaba_io,                   ONLY: open_output_file

    IMPLICIT NONE

#ifdef E4CHEM
    INTEGER             :: NSPEC
#endif
    INTEGER, PARAMETER  :: iou = 999   ! I/O unit
    INTEGER             :: status ! error status

    ! set the version number of CAABA to that of MECCA:
    caaba_version = mecca_version

    PRINT *, HLINE1
    PRINT *, '*** START OF CAABA BOX MODEL RUN (VERSION '// &
      TRIM(caaba_version)//')'
    PRINT *, HLINE1
    PRINT *

    CALL caaba_read_nml(status, iou)
    IF (status /= 0) STOP 1

    model_start = model_start_day * OneDay
    model_time  = model_start
    model_end   = model_time + runtime * OneDay ! time end

#ifdef E4CHEM
    IF (USE_MECCA)    NSPEC = NSPEC_mecca
    IF (USE_E4CHEM)   NSPEC = NSPEC_fchem
#endif

    ALLOCATE(c(NSPEC))

    ! open output file mecca_messy.nc:
    IF (.NOT. USE_TRAJECT) THEN
      CALL open_output_file(ncid_messy, 'caaba_messy', &
        (/ 'press', 'temp ', 'sza  ' /), &
        (/ 'Pa ', 'K  ', 'deg' /) )
    ENDIF

  END SUBROUTINE caaba_init

  !***************************************************************************

  SUBROUTINE caaba_physc
    ! USE caaba_mem,                ONLY: model_time, temp, press, relhum, cair
    ! USE messy_main_constants_mem, ONLY: R_gas, N_A
    ! REAL(dp) :: temp_old, press_old

    ! calc_sza renamed to time_control as it handles several
    ! aspects related to time; core functionality of calculating
    ! SZA moved to messy_main_timer as subroutine calc_sza
    !CALL calc_sza
    CALL time_control

    ! If you want to modify p,T,rh during the model run, do it here:
    !
    ! temp_old = temp
    ! press_old = press
    ! temp   = myfunction(model_time)
    ! press  = myfunction(model_time)
    ! relhum = myfunction(model_time)
    !
    ! After changing temp and/or press, the concentration of 'air'
    ! molecules [mcl/cc] must be updated:
    !
    ! cair = ...
    !
    ! To keep the mixing ratios [mol/mol] of chemical species constant,
    ! it is necessary to calculate their new concentrations [mcl/cc]:
    !
    ! c(:) = c(:) * (press*temp_old) / (press_old*temp)
    !
    ! Do not use this correction for water, instead calculate it directly:
    !
    ! c(ind_H2O) = cair * relhum * psat(temp) / press

    !*************************************************************************

  END SUBROUTINE caaba_physc

  !***************************************************************************

  SUBROUTINE caaba_finish

    USE caaba_mem, ONLY: c, USE_TRAJECT, caaba_version
    USE caaba_io,  ONLY: close_file
    USE messy_mecca,  ONLY: modver

    DEALLOCATE(c)

    IF (.NOT. USE_TRAJECT) THEN
      CALL close_file(ncid_messy)
    ENDIF

    PRINT *
    PRINT *, HLINE1
    PRINT *, '*** END OF CAABA BOX MODEL RUN (VERSION '// &
      TRIM(caaba_version)//')'
    PRINT *, HLINE1

  END SUBROUTINE caaba_finish

  !***************************************************************************

END MODULE caaba_module

!*****************************************************************************


! Time-stamp: <2018-02-21 14:47:06 sander>

! This file contains routines needed for box-trajectory procedures.
! It can only be used if the netcdf library is available.

! Author: Hella Riede, MPCH Mainz, 2006-2010
! Citation: H. Riede, P. Joeckel, and R. Sander, Geosci. Model Dev., 2,
!   267-280, 2009.

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

MODULE messy_traject_box

  USE caaba_io,         ONLY: nf90_get_att, nf90_get_var, nf90_inq_varid,    &
                              nf90_inq_dimid, nf90_inquire_dimension,        &
                              nf90_noerr, nf, open_output_file_1d,           &
                              write_output_file_1d, close_file
                              !mz_hr_20160411+
                              !nf90_noerr, nf, open_output_file, &
                              !write_output_file, close_file
                              !mz_hr_20160411- 
  !mz_hr_20160421+ 
  USE messy_main_timer, ONLY: gregor2julian, julian2gregor,& 
                              calc_sza, utc2lt
  !mz_hr_20160421- 
  USE messy_main_constants_mem, ONLY: N_A, R_gas, i4, dp, STRLEN_MEDIUM,     &
                                      STRLEN_VLONG, OneDay
  !mz_hr_20160508 USE caaba_mem, ONLY: tuf, time_string, relhum,    &
  USE caaba_mem,          ONLY: tuf, time_string, relhum, spechum,           &
                                degree_lat, degree_lon,                      &
                                model_time, model_start, model_end,          &
                                runtime, runlast, l_spechum,                 &
                                !mz_hr_20160509+
                                l_ignore_relhum, cair_old,                   &
                                !mz_hr_20160509- 
                                timesteplen, model_start_day,                &
                                c, cair, temp, press, x_j_no2, l_input_jval, &
    !mz_hr_20160421+
                                cossza, localtime, time0_jul,                &
                                lyear, lmonth, lday, lmin, lhour, lsec,      &
    !mz_hr_20160421- 
    !mz_hr_20160428+ 
                                l_psat_liquid
    !mz_hr_20160428- 

  IMPLICIT NONE

  ! LOCAL netCDF vars
  INTEGER :: dimid_time, len_time  ! dimension ID time, no. time values
  !mz_hr_20160515+ 
  !INTEGER :: varid_time, varid_press, varid_relhum, varid_temp, &
  INTEGER :: varid_time, varid_press, varid_hum, varid_temp, &
  !mz_hr_20160515- 
             varid_lat, varid_lon, varid_jno2
  INTEGER :: ncid_physc, ncid_traj, ncid_jval

  ! LOCAL variables
  INTEGER  :: trajpct      ! trajectory point counter
  INTEGER  :: intpct  = 1  ! regular integration point count
  REAL(dp) :: timesteplen_orig ! save original time step, time step varies
  REAL(dp) :: ratio1, ratio2
  REAL(dp) :: lontmp

  ! LOCAL interpolation slots
  REAL(dp) :: time1, time2     ! 2 time slots for interpolation
  !mz_hr_20160515+ 
  !REAL(dp) :: relhum1, relhum2 ! for rel./spec. humidity ...
  REAL(dp) :: hum1, hum2       ! for rel./spec. humidity ...
  !mz_hr_20160515- 
  REAL(dp) :: temp1, temp2     ! for temperature   ...
  REAL(dp) :: press1, press2   ! for pressure      ...
  REAL(dp) :: lat1, lat2       ! for latitude      ...
  REAL(dp) :: lon1, lon2       ! for longitude      ...
  REAL(dp) :: jno21, jno22     ! for J_NO2

CONTAINS

!mz_hr_20160516 changed sequence of subroutines to reflect program flow

  !***************************************************************************

  SUBROUTINE traject_init

    USE caaba_mem,        ONLY: input_physc, input_jval, l_relhum_wmo,       &
                                l_hum_emac
    !USE caaba_mem, ONLY: input_physc, input_jval
    !mz_hr_20160508+ use cair_c and cair_q instead
    !USE messy_main_tools, ONLY: cair_wmo, cair_trad
    !mz_hr_20160508-
    USE messy_main_tools, ONLY: cair_q
    USE caaba_io,         ONLY: open_input_file

    IMPLICIT NONE
    SAVE

    INTRINSIC TRIM

    ! LOCAL
    INTEGER :: status
    REAL(dp) :: localtime_jul !> local time as Julian day
    REAL(dp) :: time_sza !> model time used to calculate solar zenith angle


    !mz_hr_20160422+ 
    WRITE(*,*) "INITIALIZING TRAJECTORY MODE"
    !mz_hr_20160422- 

    !> open external input: trajectory netcdf file
    WRITE(*,*) "  opening trajectory input file ", TRIM(input_physc) !mz_hr_20160422
    CALL open_input_file(ncid_physc, input_physc)
    IF (l_input_jval) THEN
      !> open external input: j-value netcdf file
      WRITE(*,*) "  opening j-value    input file ", TRIM(input_jval) !mz_hr_20160422
      CALL open_input_file(ncid_jval, input_jval)
    ENDIF

    !> initialize variables from trajectory
    CALL get_model_start_end

    !mz_hr_20160421+ define model_time right after get_model_start_end
    model_time      = model_start
    model_start_day = model_start / OneDay

    !> calculate SZA for model_time at end of current time step
    time_sza = model_time
    !> calc_sza called before from caaba_physc; call here again after
    !>   update of longitude to get consistent localtime and sza
    CALL calc_sza(status, cossza, time_sza, time_string, degree_lon,&
                  degree_lat)
    IF (status /= 0) THEN
      WRITE(*,*) 'ERROR in traject_init: calc_sza'
      STOP 1
    ENDIF

    !> update local time variables to be consistent with time used
    !> to calculate SZA
    localtime = utc2lt(status, time_sza, degree_lon)
    IF (status /= 0) THEN
      WRITE(*,*) 'ERROR in traject_init: utc2lt'
      STOP 1
    ENDIF

    !> update local time as Julian day
    localtime_jul = time0_jul + localtime/OneDay

    !> determine local time variables year, month, ...
    CALL julian2gregor(localtime_jul, lyear, lmonth, lday, lhour, lmin, lsec)
    !WRITE(*,'(A,I4,2I0.2,1X,3I0.2)') ' local time (YMD hms): ',&
    !  lyear, lmonth, lday, lhour, lmin, lsec
    !mz_hr_20160420- 

    !mz_hr_20160508+ calculate cair from spechum
    !> update cair; conc = cair * mixing ratio
    cair = cair_q(spechum, temp, press, l_hum_emac)
    cair_old = cair

    !IF (l_relhum_wmo) THEN
    !  cair = cair_wmo(status, relhum, temp, press, l_hum_emac)
    !ELSE
    !  !mz_hr_20160501 cair = cair_trad(status, relhum, temp, press, l_hum_emac)
    !  cair = cair_trad(status, relhum, temp, press, l_hum_emac,&
    !                   l_psat_liquid)
    !ENDIF
    ! spechum checked before, no other check IF (status /= 0) STOP 1
    !mz_hr_20160508- 

    !> write variables for netCDF output file
    ! localtime unit: 33 chars 14+10(yyyy-mm-dd)+1(space)+8(hh:mm:ss)
    IF (l_input_jval) THEN
      !mz_hr_20160411+
      ! alternative open_output_file function in caaba_io
      ! so that not automatically created with dummy dimensions lon,
      ! lat, lev to leave variable names available for trajectory
      ! (non-dimensional) variables in output
      !CALL open_output_file(ncid_traj, 'caaba_messy',       &
        !(/ 'lon_tr   ', 'lat_tr   ', 'press    ', 'temp     ', &
      CALL open_output_file_1d(ncid_traj, 'caaba_messy',       &
        (/ 'lon      ', 'lat      ', 'press    ', 'temp     ', &
      !mz_hr_20160411-
           'relhum   ', 'spechum  ', 'sza      ', 'J_NO2_x  ', &
           'localtime', 'year_loc ', 'month_loc', 'day_loc  ', &
           'hour_loc ', 'min_loc  ', 'sec_loc  ' /), &
        (/ 'degrees_east                     ', &
           'degrees_north                    ', &
           'Pa                               ', &
           'K                                ', &
           '[0-1]                            ', &
           'kg/kg                            ', &
           'deg                              ', &
           '1/s                              ', &
           time_string                        , & 
           '                                 ', &
           '                                 ', &
           '                                 ', &
           '                                 ', &
           '                                 ', &
           '                                 ' /) )
    ELSE
      !mz_hr_20160411+ 
      !CALL open_output_file(ncid_traj, 'caaba_messy',       &
        !(/ 'lon_tr   ', 'lat_tr   ', 'press    ', 'temp     ', &
      CALL open_output_file_1d(ncid_traj, 'caaba_messy',       &
        (/ 'lon      ', 'lat      ', 'press    ', 'temp     ', &
      !mz_hr_20160411-
           'relhum   ', 'spechum  ', 'sza      ', 'localtime', &
           'year_loc ', 'month_loc', 'day_loc  ', 'hour_loc ', &
           'min_loc  ', 'sec_loc  ' /), &
        (/ 'degrees_east                     ', &
           'degrees_north                    ', &
           'Pa                               ', &
           'K                                ', &
           '[0-1]                            ', &
           'kg/kg                            ', &
           'deg                              ', &
           time_string                        , & 
           '                                 ', &
           '                                 ', &
           '                                 ', &
           '                                 ', &
           '                                 ', &
           '                                 ' /) )
    ENDIF
    !mz_hr_20160508 always write out spechum as well IF (l_spechum) THEN
    !ELSE
    !  IF (l_input_jval) THEN
    !    CALL open_output_file_1d(ncid_traj, 'caaba_messy',       &
    !      (/ 'lon      ', 'lat      ', 'press    ', 'temp     ', &
    !         'relhum   ', 'sza      ', 'J_NO2_x  ', 'localtime', &
    !         'year_loc ', 'month_loc', 'day_loc  ', 'hour_loc ', &
    !         'min_loc  ', 'sec_loc  ' /), &
    !...
    !  ENDIF
    !ENDIF
    !mz_hr_20160508- 

  END SUBROUTINE traject_init

  !***************************************************************************

  SUBROUTINE get_model_start_end

    USE caaba_mem,                ONLY: l_runtime_str, vlat, vlon, vpress, &
                                        vtemp, vrelhum, vspechum, vtime, &
                                        !mz_hr_20160503+
                                        !l_relhum_wmo, l_hum_emac,&
                                        l_relhum_wmo, l_hum_emac,&
                                        !mz_hr_20160503-
    !mz_hr_20160422+ 
                                        firstjan_jul, t0year, t0month, t0day,&
                                        t0hour, t0min, t0sec
    !mz_hr_20160422- 
    !mz_hr_20160426 USE messy_main_constants_mem, ONLY: TINY_DP
    !mz_hr_20160426+ 
    !USE messy_main_tools,         ONLY: spec2relhum, spec2relhumwmo, &
    !                                    cair_trad, cair_wmo
    USE messy_main_tools,         ONLY: spec2relhum, rel2spechum
    !mz_hr_20160426- 
    USE messy_main_timer,         ONLY: eval_time_str


    IMPLICIT NONE
    SAVE

    INTRINSIC :: ABS, INT, TRIM

    ! LOCAL
    INTEGER :: status

    !> inquire number of time values
    print *, "Initializing trajectory time ", TRIM(vtime)
    CALL nf(nf90_inq_dimid(ncid_physc, TRIM(vtime), dimid_time)) ! dimid of time?
    CALL nf(nf90_inquire_dimension(ncid_physc, dimid_time, len = len_time))

    !> inquire time units
    CALL nf(nf90_inq_varid(ncid_physc, TRIM(vtime), varid_time)) ! varid of time?
    CALL nf(nf90_get_att(ncid_physc, varid_time, "units", time_string))

    !> eval time units string, determine time unit conversion factor (tuf)
    !mz_hr_20160422+
    !> determine start time year, month, day, hour, minute, second
    !CALL eval_time_str(status, time_string, tuf)
    CALL eval_time_str(status, time_string, tuf, &
                       t0year, t0month, t0day, t0hour, t0min, t0sec)
    !mz_hr_20160518+ 
    !IF (status/=0) STOP 1
    IF (status/=0) THEN
      WRITE(*,*) "ERROR in TRAJECT mode: get_model_start_end: eval_time_str"
      STOP 1
    ENDIF
    !mz_hr_20160518- 

    ! determine Julian date of time origin
    time0_jul = gregor2julian(t0year, t0month, t0day, t0hour, t0min, t0sec)
    !print *, 'time0_jul = ', time0_jul

    ! determine Julian date of 01-JAN of start year (t0year)
    firstjan_jul = gregor2julian(t0year, 1, 1, 0, 0, 0)
    !mz_hr_20160422- 

    !mz_hr_20160422+ 
    !WRITE(*,'(A,F10.3)') '  time in seconds = time * ', tuf
    WRITE(*,'(A,I6)') '   time in seconds = trajectory time * ', INT(tuf)
    !mz_hr_20160422- 

    !> get trajectory end time
    CALL nf(nf90_get_var(ncid_physc, varid_time, model_end, &
      start = (/ len_time /)))
    ! convert to seconds
    model_end   = model_end * tuf

    !> get trajectory start time
    CALL nf(nf90_get_var(ncid_physc, varid_time, model_start, &
      start = (/ 1 /)))
    !convert to seconds
    model_start = model_start * tuf

    !> determine model start time
    IF (runlast > 0._dp) THEN
      IF ( (model_end - runlast*OneDay) < model_start ) THEN
        WRITE(*,'(A,F7.3,A)') '  Error get_model_start_end:&
          & runlast exceeds beginning of trajectory by ', &
          (model_start - (model_end - runlast*OneDay))/OneDay, ' days'
        STOP 1
      !mz_hr_20160422+
      !ENDIF
      ! duplicate to above IF (runlast > 0._DP) THEN
      ELSE IF ( (model_end - runlast*OneDay) > model_start ) THEN
        !mz_hr_20160422+ 
        !WRITE(*,*) 'RUNLAST    = ', runlast, ' days'
        !WRITE(*,*) '  beginning of trajectory clipped'
        WRITE(*,'(A,F7.3,A,F7.3,A)') '   runlast = ', runlast,&
          ' days    => start of trajectory clipped by ',&
          model_end/OneDay - runlast - model_start/OneDay, ' days'
      ENDIF
      !mz_hr_20160422- 
      !> update model start to 'runlast' days before model end
      model_start = model_end - runlast*OneDay
    ENDIF

    !> determine model end time
    !> if runtime and runlast specified in namelist,
    !> an inner section of the trajectory is selected
    IF (l_runtime_str) THEN
      IF ( (model_start + runtime*OneDay) > model_end) THEN
        WRITE(*,'(A,F7.3,A)') '  Error get_model_start_end:&
          & runtime exceeds end of trajectory by ', &
          ((model_start + runtime*OneDay) - model_end)/OneDay, ' days'
        STOP 1
      !mz_hr_20160422+ 
      ELSE IF ((model_start + runtime * OneDay) < model_end) THEN
        !mz_hr_20160422+ 
        !WRITE(*,*) 'RUNTIME    = ', runtime, ' days'
        !WRITE(*,*) '  clips end of trajectory runtime'
        WRITE(*,'(A,F7.3,A,F7.3,A)') '   runtime = ', runtime,&
          ' days    =>   end of trajectory clipped by ',&
          model_end/OneDay - model_start/OneDay - runtime, ' days'
        !mz_hr_20160422- 
      ENDIF
      model_end = model_start + runtime * OneDay
    ELSE
      ! runtime in days determined by end of trajectory
      runtime = (model_end - model_start)/OneDay
    ENDIF


    ! find out between which trajpts starting point is
    print *, "Interpolating trajectory information to start point"
    DO trajpct = 2, len_time, 1

      CALL nf(nf90_get_var(ncid_physc, varid_time, time2, &
          start = (/ trajpct /)))
      time2 = time2 * tuf

      IF (time2 > model_start) THEN
        ! initialize 1st and 2nd value slot of each trajectory var
        CALL nf(nf90_get_var(ncid_physc, varid_time, time1, &
          start = (/ trajpct-1 /)))
        time1 = time1 * tuf

        CALL nf(nf90_inq_varid(ncid_physc, TRIM(vtemp), varid_temp))
        CALL nf(nf90_get_var(ncid_physc, varid_temp, temp1, &
          start = (/ trajpct-1 /)))
        CALL nf(nf90_get_var(ncid_physc, varid_temp, temp2, &
          start = (/ trajpct /)))

        CALL nf(nf90_inq_varid(ncid_physc, TRIM(vpress), varid_press))
        CALL nf(nf90_get_var(ncid_physc, varid_press, press1, &
          start = (/ trajpct-1 /)))
        CALL nf(nf90_get_var(ncid_physc, varid_press, press2, &
          start = (/ trajpct /)))

        IF (l_spechum) THEN
          CALL nf(nf90_inq_varid(ncid_physc, TRIM(vspechum), varid_hum))
        ELSE
          CALL nf(nf90_inq_varid(ncid_physc, TRIM(vrelhum), varid_hum))
        ENDIF
        CALL nf(nf90_get_var(ncid_physc, varid_hum, hum1, &
          start = (/ trajpct-1 /)))
        CALL nf(nf90_get_var(ncid_physc, varid_hum, hum2, &
          start = (/ trajpct /)))

        CALL nf(nf90_inq_varid(ncid_physc, TRIM(vlat), varid_lat))
        CALL nf(nf90_get_var(ncid_physc, varid_lat, lat1, &
          start = (/ trajpct-1 /)))
        CALL nf(nf90_get_var(ncid_physc, varid_lat, lat2, &
          start = (/ trajpct /)))

        CALL nf(nf90_inq_varid(ncid_physc, TRIM(vlon), varid_lon))
        CALL nf(nf90_get_var(ncid_physc, varid_lon, lon1, &
          start = (/ trajpct-1 /)))
        CALL nf(nf90_get_var(ncid_physc, varid_lon, lon2, &
          start = (/ trajpct /)))

        ! TODO better concept for initializing many j-values
        IF (l_input_jval) THEN
          CALL nf(nf90_inq_varid(ncid_jval, "J_NO2", varid_jno2))
          CALL nf(nf90_get_var(ncid_jval, varid_jno2, jno21, &
            start = (/ trajpct-1 /)))
          CALL nf(nf90_get_var(ncid_jval, varid_jno2, jno22, &
            start = (/ trajpct /)))
        ENDIF

        !> if data contains lon up to 360deg east, set to 0
        !mz_hr_20160426+ 
        !IF (ABS(360._dp - lon1) .LE. TINY_DP) THEN
          !print *, 'lon1 = 360 deg east, set to 0'
        IF (lon1 >= 360._dp) THEN
          print *, 'WARNING: TRAJECT mode: get_model_start_end:'//&
            ' lon1 >= 360 deg east, set to 0'
          !mz_hr_20160426- 
          lon1 = 0._dp
        ENDIF
        !mz_hr_20160426+ 
        !IF (ABS(360._dp - lon2) .LE. TINY_DP) THEN
          !print *, 'lon2 = 360 deg east, set to 0'
        IF (lon2 >= 360._dp) THEN
          print *, 'WARNING: TRAJECT mode: get_model_start_end:'//&
            ' lon2 >= 360 deg east, set to 0'
          !mz_hr_20160426- 
          lon2 = 0._dp
        ENDIF

        EXIT
      ENDIF
    ENDDO

    ! save original timesteplen for the rest of the CAABA simulation
    timesteplen_orig = timesteplen

    !> interpolate physical parameters to model_start
    ratio2 = (model_start - time1)/(time2 - time1)
    ratio1 = 1._dp - ratio2
    temp       = ratio1 * temp1  + ratio2 * temp2
    press      = ratio1 * press1 + ratio2 * press2
    relhum     = ratio1 * hum1   + ratio2 * hum2
    degree_lat = ratio1 * lat1   + ratio2 * lat2
    IF ((lon2 - lon1) .LT. -180._dp) THEN
      ! crossing 0 meridian from West to East
      lontmp = lon2 + 360._dp
      degree_lon = ratio1 * lon1 + ratio2 * lontmp
      IF (degree_lon .GE. 360) degree_lon = degree_lon - 360._dp
    ELSEIF ((lon2 - lon1) .GT. 180._dp) THEN
    ! crossing 0 meridian from East to West
      lontmp = lon1 + 360._dp
      degree_lon = ratio1 * lontmp + ratio2 * lon2
      IF (degree_lon .GE. 360) degree_lon = degree_lon - 360._dp
    ELSE
      degree_lon = ratio1 * lon1 + ratio2 * lon2
    ENDIF
    IF (l_input_jval) THEN
      x_j_no2 = ratio1 * jno21   + ratio2 * jno22
    ENDIF

    IF (l_spechum) THEN
      spechum = relhum ! spechum was read into var relhum

      !mz_hr_20160508+ 
      ! check incoming spechum value
      IF (spechum >= 1.0_dp) THEN
        WRITE(*,*) 'WARNING: TRAJECT mode, get_model_start_end: '// &
                   'specific humidity >= 1 (', spechum, ')'
      ELSEIF (spechum < 0._dp) THEN
        WRITE(*,*) 'ERROR: TRAJECT mode, get_model_start_end: '// &
                   'specific humidity < 0 (', spechum, ')'
        STOP 1
      ENDIF 
      !mz_hr_20160508- 

      !mz_hr_20160506+ use spec2relhum with new arg l_relhum_wmo
      relhum = spec2relhum(status, spechum, temp, press, &
                           l_hum_emac, l_psat_liquid, l_relhum_wmo)
      !IF (l_relhum_wmo) THEN
      !  relhum = spec2relhumwmo(status, spechum, temp, press, l_hum_emac)
      !ELSE
      !  relhum = spec2relhum(status, spechum, temp, press, l_hum_emac)
      !ENDIF
      !IF (status/=0) STOP 1
      IF (status > 1) THEN
        WRITE(*,*) 'ERROR: TRAJECT mode, get_model_start_end: spec2relhum'
        STOP 1
      ENDIF
      !mz_hr_20160506- 
    ELSE
      ! check incoming relhum value
      IF (relhum >= 1.0_dp) THEN
        WRITE(*,*) 'WARNING: TRAJECT mode, get_model_start_end: '// &
                   'relative humidity >= 1 (', relhum, ')'
      ELSEIF (relhum < 0._dp) THEN
        WRITE(*,*) 'ERROR: TRAJECT mode, get_model_start_end: '// &
                   'relative humidity < 0 (', relhum, ')'
        STOP 1
      ENDIF 
      !mz_hr_20160508+ 

      spechum = rel2spechum(status, relhum, temp, press, &
                            l_hum_emac, l_psat_liquid, l_relhum_wmo)
      IF (status > 1) THEN
        WRITE(*,*) 'ERROR: TRAJECT mode, get_model_start_end: rel2spechum'
        STOP 1
      ENDIF
      !mz_hr_20160508- 
      !mz_hr_20160506+ only warning, no errors
      !IF ((relhum >= 1.0_dp) .AND. (relhum <= 1.1_dp)) THEN
      !ELSEIF (relhum > 1.1_dp) THEN
      !  WRITE(*,*) 'Error get_model_start_end: relative humidity > 1.1 : ', &
      !              relhum
      !  STOP 1
      !mz_hr_20160506-
    ENDIF

  END SUBROUTINE get_model_start_end

  !***************************************************************************

  SUBROUTINE traject_physc

  !> update variables along trajectory, such as longitude, latitude, j-values

    USE caaba_mem,        ONLY: l_relhum_wmo, l_hum_emac,&
    !mz_hr_20160422+
                                l_groundhogday, l_freezetime,                &
    !mz_hr_20160422-
    !mz_hr_20160510+ 
                                C
    !mz_hr_20160510- 
    USE messy_main_timer, ONLY: calc_sza, utc2lt
    !mz_hr_20160508+ use cair_c, cair_q
    !USE messy_main_tools, ONLY: cair_trad, cair_wmo
    USE messy_main_tools, ONLY: cair_c, cair_q, mr2spechum, spec2relhum
    USE messy_mecca_kpp, ONLY: ind_H2O
    !mz_hr_20160508-

    IMPLICIT NONE
    SAVE

    ! LOCAL variables
    INTEGER  :: status
    !mz_hr_20160421+ 
    REAL(DP) :: time_sza      !> time used to calculate SZA
    REAL(dp) :: localtime_jul !> local time as Julian day
    !mz_hr_20160421- 

    !> get data from trajectory for end of current time step
    CALL get_physc_data

    !mz_hr_20160420+ 
    !> calculate SZA for model_time at end of current time step
    time_sza = model_time + timesteplen

    !mz_hr_20160422+ 
    !> repetition of the first day?
    IF (l_groundhogday) THEN
      time_sza = (FLOOR(model_start/OneDay) +&
                  MODULO((model_time + timesteplen)/OneDay,1._DP))*OneDay
      !WRITE(*,'(A,2F16.6)') "traject_physc: modeltime groundhog",&
      !  (model_time+timesteplen)/86400., time_sza/86400.
    ENDIF

    ! l_freezetime after l_groundhog to overrule if necessary
    !> freeze time at model start?
    !>   (position along trajectory will still change)
    IF (l_freezetime) THEN
      time_sza = model_start
    ENDIF
    !mz_hr_20160422- 

    !> calc_sza called before from caaba_physc; call here again after
    !>   update of longitude to get consistent localtime and sza
    CALL calc_sza(status, cossza, time_sza, time_string, degree_lon,&
                  degree_lat)
    IF (status /= 0) THEN
      WRITE(*,*) 'ERROR in traject_physc: calc_sza'
      STOP 1
    ENDIF
    !print *, "traject_physc: time_sza = ", time_sza

    !> update local time variables to be consistent with time used
    !> to calculate SZA
    localtime = utc2lt(status, time_sza, degree_lon)
    IF (status /= 0) THEN
      WRITE(*,*) 'ERROR in traject_physc: utc2lt'
      STOP 1
    ENDIF

    !> update local time as Julian day
    localtime_jul = time0_jul + localtime/OneDay

    ! determine local time variables year, month, ...
    CALL julian2gregor(localtime_jul, lyear, lmonth, lday, lhour, lmin, lsec)
    !WRITE(*,'(A,I4,2I0.2,1X,3I0.2)') ' local time (YMD hms): ',&
    !  lyear, lmonth, lday, lhour, lmin, lsec
    !mz_hr_20160420- 

    !mz_hr_20160510+ 
    IF (l_ignore_relhum) THEN
      !> calculate spechum, relhum, cair from c(H2O)
      ! correction of cair is currently only needed in TRAJECTORY mode; to
      !   avoid introducing it in mecca_physc, instead c(H2O) is already used
      !   here to update cair, even though the calculation of cair
      !   will be done in mecca_physc later as well (necessary if
      !   CAABA is not in TRAJECT mode)
      ! NOTE: temp, press, and spechum or relhum were at this point already
      !   updated to the end of the coming KPP time step; thus, overwriting
      !   the humidities here with a humidity derived from current c(H2O)
      !   introduces a time inconsistency in data, which is unavoidable of
      !   course, because c(H2O) will only be known after the integration.

      ! update cair with current temp, press (c(H2O) still same as after x0)
      ! this cair will be used in kpp_integrate
      cair = cair_c(c(ind_H2O), temp, press, l_hum_emac)

      ! update concentrations with new cair
      C(:) = C(:) * cair/cair_old
      cair_old = cair

      ! update spechum and relhum with new cair
      spechum = mr2spechum(status, c(ind_H2O)/cair)
      IF (status > 1) THEN
        WRITE(*,*) 'ERROR: traject_physc: mr2spechum'
        STOP 1
      ENDIF
      relhum  = spec2relhum(status, spechum, temp, press, &
                            l_hum_emac, l_psat_liquid, l_relhum_wmo)
      IF (status > 1) THEN
        WRITE(*,*) 'ERROR: traject_physc: spec2relhum'
        STOP 1
      ENDIF
    ELSE
      !mz_hr_20160510- 
      !mz_hr_20160508+ always calculate cair from spechum
      ! update cair with current spechum from trajectory
      cair = cair_q(spechum, temp, press, l_hum_emac)

      ! correct concentrations with new cair
      C(:) = C(:) * cair/cair_old
      cair_old = cair
    ENDIF

    !IF (l_relhum_wmo) THEN
    !  cair = cair_wmo(status, relhum, temp, press, l_hum_emac)
    !ELSE
    !  !mz_hr_20160501 cair = cair_trad(status, relhum, temp, press, l_hum_emac)
    !  cair = cair_trad(status, relhum, temp, press, l_hum_emac,&
    !                   l_psat_liquid)
    !ENDIF
    !mz_hr_20160508- 
    !IF (status /= 0) STOP 1
    !mz_hr_20160508- 

    !> adjust concentrations to new cair so that mixing ratios are unaffected
    !>   by the physical changes (temperature, pressure, humidity (except for
    !>   H2O itself)
    !mz_hr_20160511+
    !ccorr = cair/cair_old
    !C(:) = C(:) * ccorr
    !mz_hr_20160511-
    !mz_hr_20160516+ done above
    !C(:) = C(:) * cair/cair_old
    !cair_old = cair
    !mz_hr_20160516-

  END SUBROUTINE traject_physc

  !***************************************************************************

  !> subroutine to retrieve / update data from external trajectory files
  !> (traj, init, jval); data retrieved for end of model time step to be
  !> consistent with the time written out = time at the end of the time step

  SUBROUTINE get_physc_data

    !mz_hr_20160421 not needed anymore !USE messy_main_constants_mem, ONLY: pi, TINY_DP
    !mz_hr_20160506 USE messy_main_tools, ONLY: spec2relhum, spec2relhumwmo
    USE messy_main_tools,         ONLY: spec2relhum, rel2spechum
    !mz_hr_20160506 USE caaba_mem, ONLY: l_relhum_wmo, l_psatf
    USE caaba_mem,                ONLY: l_relhum_wmo, l_hum_emac

    IMPLICIT NONE
    SAVE

    ! LOCAL
    INTEGER :: status
    LOGICAL  :: l_next_trajp !mz_hr_20160515 only used in this subroutine


    !print *, ''
    !print *, 'get_physc_data: model_time = ', model_time/tuf
    ! calculate how long the next time step should be
    IF (intpct*timesteplen_orig+model_start < time2) THEN
      ! do next regular time step, comes before next trajectory point
      l_next_trajp = .FALSE.
      timesteplen = intpct*timesteplen_orig+model_start - model_time
      !mz_hr_20160405+ moved here from below for clarity
      intpct = intpct + 1
      !mz_hr_20160405- 
      !print *, 'get_physc_data: next point kpp, timesteplen = ', timesteplen/tuf
    ELSEIF (intpct*timesteplen_orig+model_start > time2) THEN
      ! next regular time step would jump over next trajectory point
      ! -> limit jump to next trajectory point
      l_next_trajp = .TRUE.
      timesteplen = time2 - model_time
      !print *, 'get_physc_data: next point traj, timesteplen = ', timesteplen/tuf
    ELSE
      ! next regular time step lands exactly on next trajectory point
      l_next_trajp = .TRUE.
      timesteplen = time2 - model_time
      intpct = intpct + 1
      !print *, 'get_physc_data: next point traj + kpp, timesteplen = ', timesteplen/tuf
    ENDIF
    !print *, 'get_physc_data: model_time+timesteplen = ', (model_time+timesteplen)/tuf

    IF (l_next_trajp) THEN
      ! stop at next trajectory point
      temp       = temp2
      press      = press2
      relhum     = hum2
      degree_lat = lat2
      degree_lon = lon2
      IF (l_input_jval) THEN
        x_j_no2    = jno22
      ENDIF

      ! end of trajectory reached?
      IF (trajpct >= len_time) THEN
        ! end of trajectory reached
        !WRITE(*,*) 'Info get_physc_data: end of trajectory reached, '
        !WRITE(*,*) '  no shift in interpolation slots'
      ELSE
        ! trajpt count up, shift in slots
        trajpct  = trajpct + 1
        time1   = time2
        temp1   = temp2
        press1  = press2
        hum1    = hum2
        lat1    = lat2
        lon1    = lon2
        IF (l_input_jval) THEN
          jno21   = jno22
        ENDIF
        CALL nf(nf90_get_var(ncid_physc, varid_time, time2, &
          start = (/ trajpct /)))
        time2 = time2 * tuf
        CALL nf(nf90_get_var(ncid_physc, varid_temp, temp2, &
          start = (/ trajpct /)))
        CALL nf(nf90_get_var(ncid_physc, varid_press, press2, &
          start = (/ trajpct /)))
        CALL nf(nf90_get_var(ncid_physc, varid_hum, hum2, &
          start = (/ trajpct /)))
        CALL nf(nf90_get_var(ncid_physc, varid_lat, lat2, start = (/ trajpct /)))
        CALL nf(nf90_get_var(ncid_physc, varid_lon, lon2, start = (/ trajpct /)))
        
        !> if data contains lon up to 360deg east, set to 0
        !mz_hr_20160426+ 
        !IF (ABS(360._dp - lon2) .LE. TINY_DP) THEN
          !print *, 'lon2 = 360 deg east, set to 0'
        IF (lon2 >= 360._dp) THEN
          print *, 'WARNING: TRAJECT mode: get_model_start_end:'//&
            ' lon2 >= 360 deg east, set to 0'
          !mz_hr_20160426- 
          lon2 = 0._dp
        ENDIF
        IF (l_input_jval) THEN
          CALL nf(nf90_get_var(ncid_jval, varid_jno2, jno22, start = (/ trajpct /)))
        ENDIF
      ENDIF
      !print *, 'get_physc_data: time1 = ', time1, 'time2 = ', time2
    ELSE ! next regular integration point => interpolation needed
      !mz_hr_20160405+ moved above for clarity
      !intpct = intpct + 1
      !mz_hr_20160405- 

      ! calculate contribution of slot1 and slot2
      ratio2 = ((model_time+timesteplen) - time1)/(time2 - time1)
      !ratio2 = (model_time - time1)/(time2 - time1)
      ratio1 = 1._dp - ratio2
      !print *, 'ratio1 = ', ratio1, 'ratio2 = ', ratio2
      temp       = ratio1 * temp1   + ratio2 * temp2
      press      = ratio1 * press1  + ratio2 * press2
      relhum     = ratio1 * hum1    + ratio2 * hum2
      degree_lat = ratio1 * lat1    + ratio2 * lat2

      IF ((lon2 - lon1) .LT. -180._dp) THEN
        ! crossing 0 meridian from West to East
        lontmp = lon2 + 360._dp
        degree_lon = ratio1 * lon1 + ratio2 * lontmp
        IF (degree_lon .GE. 360) degree_lon = degree_lon - 360._dp
      ELSEIF ((lon2 - lon1) .GT. 180._dp) THEN
      ! crossing 0 meridian from East to West
        lontmp = lon1 + 360._dp
        degree_lon = ratio1 * lontmp + ratio2 * lon2
        IF (degree_lon .GE. 360) degree_lon = degree_lon - 360._dp
      ELSE
        degree_lon = ratio1 * lon1 + ratio2 * lon2
      ENDIF

      IF (l_input_jval) THEN
        x_j_no2 = ratio1 * jno21   + ratio2 * jno22
      ENDIF
    ENDIF

    IF (l_spechum) THEN
      spechum = relhum ! spechum was read into var relhum

      !mz_hr_20160508+ 
      ! check incoming spechum value
      IF (spechum >= 1.0_dp) THEN
        WRITE(*,*) 'WARNING: TRAJECT mode, get_physc_data: '// &
                   'specific humidity >= 1 (', spechum, ')'
      ELSEIF (spechum < 0._dp) THEN
        WRITE(*,*) 'ERROR: TRAJECT mode, get_physc_data: '// &
                   'specific humidity < 0 (', spechum, ')'
        STOP 1
      ENDIF 
      !mz_hr_20160508- 

      !mz_hr_20160506+ use spec2relhum with new arg l_relhum_wmo
      relhum = spec2relhum(status, spechum, temp, press, &
                           l_hum_emac, l_psat_liquid, l_relhum_wmo)
      !IF (l_relhum_wmo) THEN
      !  relhum = spec2relhumwmo(status, spechum, temp, press, l_hum_emac)
      !ELSE
      !  relhum = spec2relhum(status, spechum, temp, press, l_hum_emac)
      !ENDIF
      !IF (status/=0) THEN
      IF (status > 1) THEN
      !mz_hr_20160506- 
        WRITE(*,*) 'ERROR: TRAJECT mode, get_physc_data: spec2relhum'
        STOP 1
      ENDIF
    ELSE
      ! check incoming relhum value
      !mz_hr_20160506+ only warning if relhum >= 1
      IF (relhum >= 1.0_dp) THEN
        WRITE(*,*) 'WARNING: TRAJECT mode, get_physc_data: '// &
                   'relative humidity >= 1 (', relhum, ')'
      ELSEIF (relhum < 0._dp) THEN
        WRITE(*,*) 'ERROR: TRAJECT mode, get_physc_data: '// &
                   'relative humidity < 0 (', relhum, ')'
        STOP 1
      ENDIF
      !IF ((relhum > 1.0_dp) .AND. (relhum <= 1.1_dp)) THEN
      !ELSEIF (relhum > 1.1_dp) THEN
      !  WRITE(*,*) 'Error get_physc_data: relative humidity > 1.1 : ', &
      !    relhum
      !  STOP 1
      !mz_hr_20160506- 

      !mz_hr_20160508+ 
      spechum = rel2spechum(status, relhum, temp, press, &
                            l_hum_emac, l_psat_liquid, l_relhum_wmo)
      IF (status > 1) THEN
        WRITE(*,*) 'ERROR: TRAJECT mode, get_physc_data: rel2spechum'
        STOP 1
      ENDIF
      !mz_hr_20160508- 
    ENDIF

  END SUBROUTINE get_physc_data

  !***************************************************************************

  SUBROUTINE traject_result

  !> Write out physical data such as temperature, pressure, ... to file
  !>   caaba_messy.nc.

    USE caaba_mem,  ONLY: localtime, percent_done, degree_sza,               &
                          lyear, lmonth, lday, lmin, lhour, lsec
    USE messy_main_constants_mem, ONLY: FLAGGED_BAD

    IMPLICIT NONE


    !mz_hr_20160508+ always write out spechum AND relhum
    IF (l_input_jval) THEN
      ! write out external J_NO2
      !mz_hr_20160414 use 1d version instead of write_output_file
      CALL write_output_file_1d(ncid_traj, model_time,  &
        (/degree_lon, degree_lat, press, temp, relhum, spechum, &
          degree_sza, x_j_no2, (localtime/tuf), REAL(lyear,dp), &
          REAL(lmonth,dp), REAL(lday,dp), REAL(lhour,dp),       &
          REAL(lmin,dp), REAL(lsec,dp) /))
    ELSE
      !mz_hr_20160414 CALL write_output_file(ncid_traj, model_time,  &
      CALL write_output_file_1d(ncid_traj, model_time,  &
        (/degree_lon, degree_lat, press, temp, relhum, spechum, &
          degree_sza, (localtime/tuf), REAL(lyear,dp), REAL(lmonth,dp), &
          REAL(lday,dp), REAL(lhour,dp), REAL(lmin,dp), REAL(lsec,dp) /))
    ENDIF

    !mz_hr_20160508+ no special case for 1st time point necessary anymore
    !IF (percent_done > 0.) THEN
    !ELSE ! 1st output, localtime not defined yet
    !  IF (l_spechum) THEN
    !    IF (l_input_jval) THEN
    !      !mz_hr_20160414+ 
    !      !CALL write_output_file(ncid_traj, model_time,  &
    !      CALL write_output_file_1d(ncid_traj, model_time,  &
    !        (/degree_lon, degree_lat, press, temp, relhum, spechum, &
    !          FLAGGED_BAD, x_j_no2, FLAGGED_BAD, FLAGGED_BAD, FLAGGED_BAD, &
    !          FLAGGED_BAD, FLAGGED_BAD, FLAGGED_BAD, FLAGGED_BAD /))
    !      !mz_hr_20160414- 
    !    ELSE
    ! ...

  END SUBROUTINE traject_result

  !***************************************************************************

  SUBROUTINE traject_finish

    CALL close_file(ncid_physc)
    CALL close_file(ncid_traj)
    IF (l_input_jval) THEN
      CALL close_file(ncid_jval)
    ENDIF

  END SUBROUTINE traject_finish

  !***************************************************************************

END MODULE messy_traject_box

!*****************************************************************************

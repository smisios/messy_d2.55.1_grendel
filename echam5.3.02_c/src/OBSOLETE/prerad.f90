SUBROUTINE prerad

  ! Description:
  !
  ! This subroutine computes some parameters for the radiation
  ! scheme which remain unchanged during scan over latitude lines.
  !
  ! Method:
  !
  ! *prerad* is called from *scan1*.
  !
  ! Externals:
  !   *orbit*   compute orbital parameters.
  !   *pre_o3*   compute parameters for ozone.
  !
  ! Authors:
  !
  ! M. Jarraud, ECMWF, June 1983, original source
  ! L. Kornblueh, MPI, May 1998, f90 rewrite
  ! U. Schulzweida, MPI, May 1998, f90 rewrite
  ! U. Schlese, DKRZ and M. Esch, MPI, July 1999, ECHAM5-modifications
  ! M.A. Giorgetta, MPI, March 2000, revision for ECHAM5 radiation
  ! I. Kirchner, MPI, December 2000, time management, new orbit
  ! L. Kornblueh, MPI, January 2001, fixes PCMDI orbit to original specs. 
  ! I. Kirchner, MPI, April 2001, revision
  ! U. Schulzweida, MPI, May 2002, blocking (nproma)
  ! L. Kornblueh, MPI, October 2003, parallelized radiation statistics
  ! M.A. Giorgetta, MPI, May 2005, clean up to make "lrad" working
  ! 
  ! for more details see file AUTHORS
  !

  USE mo_kind           ,ONLY: dp
  USE mo_exception      ,ONLY: message, message_text
  USE mo_control        ,ONLY: nlevp1
  USE mo_physc1         ,ONLY: cdisse ,czen1 ,czen2 ,czen3, &
                               cdissem,czen1m,czen2m,czen3m
  USE mo_constants      ,ONLY: api
  USE mo_param_switches ,ONLY: lrad
  USE mo_radiation      ,ONLY: nmonth,io3, iaero
  USE mo_o3_lwb         ,ONLY: pre_o3_lwb
  USE mo_o3clim         ,ONLY: pre_o3clim, pre_o3clim_3
  USE mo_time_control,   ONLY : l_trigrad, l_diagrad,                       &
                                l_orbvsop87, get_year_len, get_clock,       &
                                get_month_len, lstart, lresume,             &
                                get_year_day, radiation_date, current_date, &
                                get_date_components, start_date
  USE mo_time_conversion,ONLY: TC_get, day_len
  USE mo_time_base,      ONLY: get_calendar_type, JULIAN, CYL360,           &
                               julian_date, ly360_date,                     &
                               Set_JulianDay, Set_Ly360Day
  USE mo_orbit,          ONLY: orbit
  USE mo_geoloc,         ONLY: amu0_x, rdayl_x, amu0m_x, rdaylm_x
  USE m_solang,          ONLY: solang_extended_const
  USE mo_diag_radiation, ONLY: init_diag_radiation

  IMPLICIT NONE

  TYPE(julian_date) :: jul_date_now, jul_date_pal, jul_date_perp
  TYPE(ly360_date)  :: ly360_date_now, ly360_date_pal
  REAL(dp):: zclock, zclockm, zdoy, zdoyve, zdoyvem, zdoym, zytime
  REAL(dp):: zra, zdec, zdis, zgha
  REAL(dp):: zvemar0, zyearve
  INTEGER :: yr, mo, dy, hr, mn, se, jday, jsec, imlen, imlenh, itime

 !--  Compute middle of the month for perpetual runs

  IF (nmonth /= 0) THEN   
    imlen  = get_month_len(1987,nmonth)
    imlenh = imlen/2 + 1
    IF (MOD(imlen,2) == 0) THEN
      itime = 0
    ELSE
      itime = 12
    END IF
    CALL Set_JulianDay(1987, nmonth, imlenh, itime, jul_date_perp)
  END IF

  !--  Compute orbital parameters for present time step

  zclock = get_clock(current_date)

  CALL get_date_components(current_date, yr, mo, dy, hr, mn, se)
  CALL TC_get             (current_date, jday, jsec)
    
  IF (l_orbvsop87) THEN
    ! VSOP87 orbit

    CALL Set_JulianDay(yr, mo, dy, jsec, jul_date_now)
    zdoy = jul_date_now%day + jul_date_now%fraction
    CALL orbit (zdoy, zra, zdec, zdis, zgha)
    cdisse = 1.0_dp/zdis**2

  ELSE
    ! PCMDI-Orbit

    CALL get_date_components(start_date,  yr ,mo ,dy ,hr ,mn ,se)

    SELECT CASE (get_calendar_type())
    CASE (JULIAN)
      CALL Set_JulianDay(yr, mo, dy, jsec, jul_date_now)
      CALL Set_JulianDay(1900, 1, 1, 0, jul_date_pal)
    CASE (CYL360)
      CALL Set_Ly360Day(yr, mo, dy, jsec, ly360_date_now)
      CALL Set_Ly360Day(1900, 1, 1, 0, ly360_date_pal)
    END SELECT
    
    IF(nmonth == 0) THEN

      SELECT CASE (get_calendar_type())
      CASE (JULIAN)
        zdoy =  (jul_date_now%day + jul_date_now%fraction) &
              - (jul_date_pal%day + jul_date_pal%fraction)
      CASE (CYL360)
        zdoy =  (ly360_date_now%day + ly360_date_now%fraction) &
              - (ly360_date_pal%day + ly360_date_pal%fraction)
      END SELECT

    ELSE                   ! perpetual month
    
      zdoy =  (jul_date_perp%day + jul_date_perp%fraction) &
            - (jul_date_pal%day  + jul_date_pal%fraction)

    END IF

    ! Determine day of vernal equinox in March 1900

     zvemar0 = 78.41_dp - 0.0078_dp*(1900-1987) + 0.25_dp*MOD(1900,4)
     zyearve = zdoy + get_year_len() - zvemar0
     zdoyve = MOD(zyearve/get_year_len(),1.0_dp)*2._dp*api
     CALL orbit (zdoyve, cdisse, zdec, zra)

  END IF

  czen1 = SIN(zdec)
  czen2 = COS(zdec)*COS(zclock)
  czen3 = COS(zdec)*SIN(zclock)

  CALL solang_extended_const(czen1, czen2, czen3, amu0_x, rdayl_x)

!  WRITE (message_text,'(a,6e28.18)') '1  ', &
!       doy, zclock, cdisse,  czen1,  czen2,  czen3
!  CALL message('prerad',message_text)

  IF (lrad) THEN

  IF (l_trigrad) THEN

     !--- prepare spectral coeffs for London o3 climatology

     IF ( io3 == 2 ) THEN
        IF (nmonth /= 0) THEN

           zytime = MOD(zdoy/get_year_len(),1.0_dp)*2._dp*api
        ELSE
           CALL get_date_components(radiation_date,yr,mo,dy,hr,mn,se)
           zytime = 2._dp*api*get_year_day(radiation_date)/get_year_len(yr)
        END IF
        CALL pre_o3_lwb(zytime)
     ELSE IF ( io3 == 3 ) THEN
        CALL pre_o3clim_3
     ELSE IF ( io3 == 4 ) THEN
        CALL pre_o3clim
     END IF

     !--  Compute orbital parameters for full radiation time step.

     zclockm = get_clock(radiation_date)
     CALL get_date_components(radiation_date, yr, mo, dy, hr, mn, se)
     CALL TC_get             (radiation_date, jday, jsec)

     IF (l_orbvsop87) THEN
       ! VSOP87 orbit

       CALL Set_JulianDay(yr, mo, dy, jsec, jul_date_now)
       zdoym = jul_date_now%day + jul_date_now%fraction
       CALL orbit (zdoym, zra, zdec, zdis, zgha)
       cdissem = 1.0_dp/zdis**2
       
     ELSE
       ! PCMDI-Orbit

       CALL get_date_components(start_date,  yr ,mo ,dy ,hr ,mn ,se)
       SELECT CASE (get_calendar_type())
       CASE (JULIAN)
         CALL Set_JulianDay(yr, mo, dy, jsec, jul_date_now)
         CALL Set_JulianDay(1900, 1, 1, 0, jul_date_pal)
       CASE (CYL360)
         CALL Set_Ly360Day(yr, mo, dy, jsec, ly360_date_now)
         CALL Set_Ly360Day(1900, 1, 1, 0, ly360_date_pal)
       END SELECT

       IF(nmonth == 0) THEN

         SELECT CASE (get_calendar_type())
         CASE (JULIAN)
           zdoym =  (jul_date_now%day + jul_date_now%fraction) &
                  - (jul_date_pal%day + jul_date_pal%fraction)
         CASE (CYL360)
           zdoym =  (ly360_date_now%day + ly360_date_now%fraction) &
                  - (ly360_date_pal%day + ly360_date_pal%fraction)
         END SELECT

       ELSE                   ! perpetual month
    
         zdoym =  (jul_date_perp%day + jul_date_perp%fraction) &
                - (jul_date_pal%day  + jul_date_pal%fraction)

       END IF

       ! Determine day of vernal equinox in March 1900

       zvemar0 = 78.41_dp - 0.0078_dp*(1900-1987) + 0.25_dp*MOD(1900,4)
       zyearve = zdoym + get_year_len() - zvemar0
       zdoyvem = MOD(zyearve/get_year_len(),1.0_dp)*2.*api
       CALL orbit (zdoyvem, cdissem, zdec, zra)

     END IF
     czen1m = SIN(zdec)
     czen2m = COS(zdec)*COS(zclockm)
     czen3m = COS(zdec)*SIN(zclockm)

     CALL solang_extended_const(czen1m, czen2m, czen3m, amu0m_x, rdaylm_x)

  END IF

  END IF

!  WRITE (message_text,'(a,6e28.18)') '1m ', &
!       doym, zclockm, cdissem, czen1m, czen2m, czen3m
!  CALL message('prerad',message_text)
!

     !--- Initialize diagnostics

     IF (l_diagrad .AND. (lstart .OR. lresume)) THEN
       CALL init_diag_radiation
     END IF

END SUBROUTINE prerad

!**************************************************************************
!
! This file is part of the Chemical Lagrangian Model of the
! Stratosphere (CLaMS). CLaMS is a hierarchy of numerical models and
! preprocessors which simulate transport, chemical reactions and
! mixing processes in the stratosphere.  
!
! Copyright (c) 2016
! Felix Ploeger, Paul Konopka, Nicole Thomas, Baerbel Vogel
! Forschungszentrum Juelich GmbH
! Last Modified By: N.Thomas
! Last Modified On: Wed Jul  1 11:25:54 2020
!
! This program is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 2 of the License, or
! (at your option) any later version. This program is distributed in
! the hope that it will be useful, but WITHOUT ANY WARRANTY; without
! even the implied warranty of MERCHANTABILITY or FITNESS FOR A
! PARTICULAR PURPOSE. See the GNU General Public License for more
! details. You should have received a copy of the GNU General Public
! License along with this program; if not, see <https://www.gnu.org/licenses/>.
!
!********************************************************************************!
! PURPOSE: Add artificial tracers for CLaMS
!********************************************************************************!
!
!    Subroutine set_artificial_tracers (status, LAT, LON, LEV, TRACERS, SGM)
!    integer FUNCTION def_overwrite (date_js, pulse_date, pulse_length, &
!                                    pulse_reset_period, pulse_time_chk,&
!                                    step_perp)  
!    Subroutine error_handler (error)
!    Subroutine clamstracer_read_nml (status, iou)
!    Subroutine clamstracer_read_tracerlist (status, iou, TRACERS)
!
!********************************************************************************!

Module messy_clamstracer

  IMPLICIT NONE

  CHARACTER(LEN=*),PARAMETER, PUBLIC :: modstr = 'clamstracer'
  CHARACTER(LEN=*),PARAMETER, PUBLIC :: modver = '1.0'

CONTAINS 

  
  !**************************************************************************
  !
  !**************************************************************************
  !Subroutine set_artificial_tracers (status, LAT, LON, LEV, TRACERS, SGM)
  Subroutine set_artificial_tracers (status, LAT, LON, LEV, THETA_TROP1, TRACERS, SGM)
    
    USE messy_clams_global,         ONLY: YEAR, MONTH, DAY, HOUR, MINUTE, SECOND
    USE messy_clams_global,         ONLY: rank, dp, prec, mdi, eps
    USE messy_clamstracer_global,   ONLY: tracer_type, ntracer, tracertype, &
                                          lev_bound, step_perp, lsm, gph, &
                                          timestep_tracer, set_to_zero
    USE messy_clamstracer_regional_masks, ONLY: define_regional_mask
    USE messy_clams_tools_dateconv, ONLY: ymds2js

    IMPLICIT NONE
    
    REAL(PREC),         DIMENSION(:), POINTER :: LAT, LON, LEV, SGM
    REAL(PREC),         DIMENSION(:), POINTER :: THETA_TROP1
    type(tracer_type),  DIMENSION(:), POINTER :: TRACERS
    INTEGER                                   :: status, jj

!!!!! festes Datum oder einlesen ???
    ! date0 (used for delta_time)
    integer, parameter :: year0=1979, month0=1, day0=1, hour0=12

    LOGICAL, DIMENSION(:), POINTER :: mask(:),mask_set_zero(:)
    LOGICAL, DIMENSION(:), POINTER :: mask_inv(:) !!! FP--170220

    LOGICAL      :: seasonal_pulsing ! to decide whether the pulsing is on the
!same day every year (reset time in years)  or after a given number of days (reset time in days)
    REAL(DP)     :: date_js, start_js
    REAL(PREC)   :: delta_time, pii
    REAL(PREC)   :: new_value
    REAL(DP)     :: seconds  
    INTEGER      :: finish, start, counts_per_sec
    INTEGER      :: itracer

    status = 0 ! no error

    if (rank==0)  then
       CALL system_clock (count=start)
       write (*,*) 
       write (*,*) 'START OF TRACER'
       write (*,*) 
    endif

    if (maxval(theta_trop1) <= 0.) then
       if (rank==0)   write (*,*) 'WARNING: no valid values on THETA_TROP1 !'
    endif
    
    pii = 4*ATAN(1.)
 
    start_js = ymds2js (year0, month0, day0, hour0*3600)
    if (rank==0) write (*,*) 'start date: ',hour0*3600,day0,month0,year0

    date_js = ymds2js (YEAR, MONTH, DAY, HOUR*3600+MINUTE*60+SECOND)
    if (rank==0) write (*,*) 'current date: ', HOUR*3600+MINUTE*60+SECOND,DAY,MONTH,YEAR

    ! offset in days
    delta_time = (date_js - start_js) / 86400.
    if (rank==0) write (*,*) 'offset in days: ', delta_time


    allocate (mask(size(LAT)))
    allocate (mask_inv(size(LAT))) !!! FP--170220
    allocate (mask_set_zero(size(LAT))) ! !ALJP 190318
    

    do itracer = 1, ntracer

       !-------------------------------------------------
       ! set regional mask
       !-------------------------------------------------
 
       mask = .false.
       mask_set_zero = .false.

       select case (tracertype)
       case ('grid') 
          where (TRACERS(itracer)%latmin <= lat .and. lat <= TRACERS(itracer)%latmax .and. &
                 TRACERS(itracer)%lonmin <= lon .and. lon <= TRACERS(itracer)%lonmax .and. &
                 TRACERS(itracer)%levmin <= lev .and. lev <= TRACERS(itracer)%levmax)
             mask = .true. 
          end where
          where (.not.(TRACERS(itracer)%latmin <= lat .and. lat <= TRACERS(itracer)%latmax .and. &
                 TRACERS(itracer)%lonmin <= lon .and. lon <= TRACERS(itracer)%lonmax) .and. &
                 TRACERS(itracer)%levmin <= lev .and. lev <= TRACERS(itracer)%levmax)
             mask_set_zero = .true.
          end where
       case ('map')
          if (TRACERS(itracer)%seg_no == 0) then
             where (abs((SGM-mdi)/mdi) > eps)
                mask = .true.
             endwhere
          else
             where (TRACERS(itracer)%seg_no == nint(sgm))
                mask = .true.
             endwhere
          endif
       case ('region')
          call define_regional_mask (lon, lat, lev, lev_bound, &
                                     TRACERS(itracer)%name, mask, lsm, gph)
       case default
          write(*,*) '   ERROR: Type not defined!'
          stop
       end select
         
       mask_inv = .true. !!! FP--170220
! op_pj_20190925+: wrong syntax:       
!!$    where (mask .eq. .true.)
       where (mask)
! op_pj_20190925- (wouldn't simply  "mask_inv = .NOT. mask" be correct?)
          mask_inv = .false.
       end where


       !-------------------------------------------------
       ! check how to overwrite the particular tracer
       !-------------------------------------------------
       seasonal_pulsing = .true.
       if (TRACERS(itracer)%type_pulsing == 1) seasonal_pulsing = .false.
       if (TRACERS(itracer)%pulse_length > 0) then

          TRACERS(itracer)%overwrite = def_overwrite (date_js, TRACERS(itracer)%pulse_date, &
                                                TRACERS(itracer)%pulse_length, &
                                                TRACERS(itracer)%pulse_reset_period, &
                                                TRACERS(itracer)%pulse_time_chk, &
                                                TRACERS(1)%pulse_date, &
                                                step_perp, seasonal_pulsing) 

       else

          TRACERS(itracer)%overwrite = 1
          
       endif
       
       !-------------------------------------------------
       ! set value in region
       !-------------------------------------------------

!!!!!
       ! reset to 0 below lev_bound
       if (tracertype == 'region') then
          where (lev <= lev_bound)
             TRACERS(itracer)%values= 0
          end where
       endif
!!!!!!!!
       ! reset to 0 in the layer outside of the geographic region
       if (set_to_zero == 1) then
             where (mask_set_zero)  TRACERS(itracer)%values =0.
       end if

!!!!! Fuer 'map' kein Zuruecksetzen noetig ?!


       ! reset to 0  
       if (TRACERS(itracer)%overwrite == 2) then
          TRACERS(itracer)%values = 0
          if (rank==0) write(*,*)  TRACERS(itracer)%name, '> outside region: = 0'        
       endif

       if (TRACERS(itracer)%overwrite == 0) then

          if (rank==0) write(*,*) TRACERS(itracer)%name, '> in region: = 0'

          where (mask) TRACERS(itracer)%values = 0 

       else if (TRACERS(itracer)%overwrite == 1 .or. TRACERS(itracer)%overwrite == 2) then

!!! AP--190318 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
             
          if ((TRACERS(itracer)%const > eps) .or. (TRACERS(itracer)%lin > eps) .or. (TRACERS(itracer)%quad > eps))  then
             if (TRACERS(itracer)%amplitude < eps) TRACERS(itracer)%periode =1.

             new_value = TRACERS(itracer)%const + TRACERS(itracer)%lin *delta_time + TRACERS(itracer)%lin * delta_time**2 + TRACERS(itracer)%amplitude *COS((2*pii/TRACERS(itracer)%periode)*delta_time+TRACERS(itracer)%phase)
             if (rank==0) write(*,*) TRACERS(itracer)%name,'> in region: =',new_value

             where (mask) TRACERS(itracer)%values = new_value
   
          else
             ! bereits in clamstracer_read_tracerlist abgefangen !
             if (rank==0) write (*,'(A,I3,A)') 'WARNING: TRACERS no ',itracer,': const or lin or quad are not specified !'
             !status = -1
             !return
          endif

          if (TRACERS(itracer)%tau > eps) then
            
          do jj=1,ntracer
           !   print*, 'parent species name',
           !   TRIM(species(jj)%parent_species_name)

              if (TRIM(TRACERS(jj)%parent_specie_name) .eq. TRIM(TRACERS(itracer)%name)) then
             !  print*, 'species ', TRIM(species(jj)%name), ' has for parent ',
             !  TRIM(species(i)%name)
               where (mask_inv)  TRACERS(jj)%values = TRACERS(jj)%values + &
TRACERS(itracer)%values *(1.- exp(-1.*TRACERS(itracer)%tau*timestep_tracer))
               if (rank==0) write (*,*) TRACERS(jj)%name, '> in region: +', &
                    TRACERS(itracer)%values*exp(-1.*TRACERS(itracer)%tau*timestep_tracer)
              end if
          end do

 
             if (rank==0) write (*,*) TRACERS(itracer)%name, '> in region: *', &
                                      exp(-1.*TRACERS(itracer)%tau*timestep_tracer)

             where (.not. mask) TRACERS(itracer)%values = &
                     TRACERS(itracer)%values * exp(-1.*TRACERS(itracer)%tau*timestep_tracer) 
             
          endif

       end if
!!! FP--170220 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    end do


    deallocate (mask)

    IF (rank==0)THEN

       WRITE (*,*) 'Normal termination of TRACER'  
       CALL system_clock (COUNT_RATE=counts_per_sec)
       CALL system_clock (count=finish)
       IF (finish > start) THEN
          seconds = float(finish-start) / float(counts_per_sec)
       ELSE
          seconds = 0.0
       ENDIF
       WRITE (*,*)
       WRITE (*,'(A,F10.2,A)') 'This job has taken ',seconds,' seconds to execute.' 
       
    ENDIF

  End Subroutine set_artificial_tracers

 !*************************************************************************
  ! check how to overwrite the particular tracer
  !*************************************************************************
  integer FUNCTION def_overwrite (date_js, pulse_date, pulse_length, &
                                  pulse_reset_period, pulse_time_chk,&
                                  pulse_date_firstyear, step_perp,seasonal_pulsing)

    USE messy_clams_global,         ONLY: DP, PREC, mdi, eps
    USE messy_clams_tools_dateconv, ONLY: ymds2js

    implicit none

    logical, intent(in) :: seasonal_pulsing
    character(*), intent(in) :: pulse_date, pulse_date_firstyear
    integer,      intent(in) :: pulse_reset_period, pulse_length, step_perp
    real(DP),     intent(in) :: date_js
    real(PREC),   intent(in) :: pulse_time_chk

    character(4) :: y_str, hy_str
    character(2) :: m_str, d_str, h_str
    character(10),  dimension(:), allocatable :: reset_dates
    integer :: i, ios
    integer, parameter :: nmax = 50 ! maximum number of resets for each tracer
    integer  :: y, m, d, h, s
    real(DP) :: js_l, js_u, prec_js = 86400. ! prec_js is there to ensure
!pulsing happens even if the subroutine is not called at the precise pulsing time; should be in the end equal to the time step at which add_artificial_tracer is called (read from general namelist?)
    integer :: diff_date_jd, reminder_days
    integer  :: pulse_first_year, y_ref

    def_overwrite = 0 ! 1=> pulsing; 2=> resetting to 0 then pulsing

    if (seasonal_pulsing) then ! below for seasonal pulsing, reset time in years


       y_str = trim(pulse_date(1:4))
       read(y_str,*,iostat=ios) y
       
       allocate (reset_dates(nmax))
       
       do i=1,nmax
          write(hy_str, "(i4)") y+pulse_reset_period*(i-1)
          reset_dates(i) = trim(hy_str)//trim(pulse_date(5:10))
       end do
       
       ! this is for seasonal tracers; repetition is yearly

       !  !!! - FP-160613
       ! offset in #(yr) to first year (--> defines perpetuum cycle for pulse)
       pulse_first_year = 0
       if (step_perp > 0) then ! only for perpetuum runs! 
          hy_str = trim(pulse_date_firstyear(1:4))
          read(hy_str,*,iostat=ios) y_ref
          pulse_first_year = y - y_ref
       endif
       
       do i=1,nmax
          
          y_str = trim(reset_dates(i)(1:4))
          m_str = trim(reset_dates(i)(5:6))
          d_str = trim(reset_dates(i)(7:8))
          h_str = trim(reset_dates(i)(9:10))
          read(y_str,*,iostat=ios) y
          read(m_str,*,iostat=ios) m
          read(d_str,*,iostat=ios) d
          read(h_str,*,iostat=ios) h
          s = h*3600
          
          js_l = ymds2js (y, m, d, s)
          js_u = js_l + (pulse_length-1)*86400d0
          
          if (date_js >= js_l .and. date_js < js_l+prec_js) then
             def_overwrite = 2
          else if (date_js > js_l .and. date_js <= js_u) then
             def_overwrite = 1
          end if
          
       end do

       deallocate (reset_dates)

    else ! not seasonal pulsing, reset time in days

           y_str = trim(pulse_date(1:4))
           m_str = trim(pulse_date(5:6))
           d_str = trim(pulse_date(7:8))
           h_str = trim(pulse_date(9:10))
           read(y_str,*,iostat=ios) y
           read(m_str,*,iostat=ios) m
           read(d_str,*,iostat=ios) d
           read(h_str,*,iostat=ios) h
           s = h*3600


           js_l = ymds2js (y, m, d, s) ! pulse initial date in julian seconds

           diff_date_jd = int((date_js - js_l)/86400d0) ! diff in days
           !for perpetuum runs (if step_perp > 0): add the perpetuum years
           if (step_perp > 0) then
              diff_date_jd = diff_date_jd + step_perp*int((ymds2js (y+1, 1, 1,0) - ymds2js (y, 1, 1, 0))/86400d0)
           end if
           reminder_days = modulo(diff_date_jd, pulse_reset_period)

           if ((reminder_days >= 0) .and. (reminder_days < pulse_length)) then
              if ((reminder_days == 0) .or. abs((pulse_time_chk-mdi)/mdi)<=eps) then
              ! this case occurs if simulation starts e.g. 1979010112 and first
              ! pulse date is also 1979010112 --> first call of age_spectrum for
              ! 1979010212...
                def_overwrite = 2
              else
                def_overwrite = 1
              end if
           end if

        end if

  END FUNCTION def_overwrite

  !*******************************************************************
  ! 
  !*******************************************************************
  Subroutine error_handler (error)

    USE messy_clams_global, ONLY: rank

    implicit none

    integer :: error

    if (rank==0) then
       write (*,*) '***********************************************************'

       if (error>20 .and. error<50) then
          write (*,*) 'ERROR in tracer description file :'
       else
          write (*,*) 'ERROR:'
       endif
       
       SELECT CASE (error)
       CASE (20)
          write (*,*) '    Cannot open bmix-boundfile list !!!'
       CASE (21)
          WRITE (*,*) '    name of species could not be read !!!'
       CASE (22)
          WRITE (*,*) '    seg-no/region could not be read !!!'
       CASE (23)
          WRITE (*,*) '    lower latitude-range could not be read !!!'
       CASE (24)
          WRITE (*,*) '    upper latitude-range could not be read !!!'
       CASE (25)
          WRITE (*,*) '    lower longitude-range could not be read !!!'
       CASE (26)
          WRITE (*,*) '    upper longitude-range could not be read !!!'
       CASE (27)
          WRITE (*,*) '    lower level-range could not be read !!!'
       CASE (28)
          WRITE (*,*) '    upper level-range could not be read !!!'
       CASE (29)
          WRITE (*,*) '    pulse date could not be read !!!'
       CASE (30)
          WRITE (*,*) '    pulse length could not be read !!!'
       CASE (31)
          WRITE (*,*) '    pulse reset period could not be read !!!'
       CASE (32)
          WRITE (*,*) '    "const" could not be read !!!'
       CASE (33)
          WRITE (*,*) '    "lin" could not be read !!!'
       CASE (34)
          WRITE (*,*) '    "quad" could not be read !!!'
       CASE (35)
          WRITE (*,*) '    "tau" could not be read !!!'
       CASE (36)
          !WRITE (*,*) '    "const","lin","quad" and "tau" = 0. !!!'
          WRITE (*,*) '    "const" and "lin" = 0. !!!'
       CASE (41)
          WRITE (*,*) '    year (pulse date) could not be read !!!'
       CASE (42)
          WRITE (*,*) '    month (pulse date) could not be read !!!'
       CASE (43)
          WRITE (*,*) '    day (pulse date) could not be read !!!'
       CASE (44)
          WRITE (*,*) '    hour (pulse date) could not be read !!!'
       CASE (51)
          WRITE (*,*) '    "periode" could not be read !!!'
       CASE (52)
          WRITE (*,*) '    "amplitude" could not be read !!!'
       CASE (53)
          WRITE (*,*) '    "phase" could not be read !!!'
       CASE (54)
          WRITE (*,*) '    "type_pulsing" could not be read !!!'
       case default
          write (*,*) 'error no ',error
       end select
      write (*,*) '***********************************************************'
    endif

  end subroutine error_handler
 
  !**************************************************************************
  !
  !**************************************************************************
  Subroutine clamstracer_read_nml (status, iou)
      
    USE messy_main_tools,    ONLY: read_nml_open, read_nml_check, read_nml_close
    USE messy_clams_global,  ONLY: rank, ldiagout
    USE messy_clamstracer_global
    
    IMPLICIT NONE
    
    !I/O
    INTEGER, INTENT(OUT) ::   status
    INTEGER, INTENT(IN)  ::   iou
    
    !LOCAL
    CHARACTER(LEN=*),PARAMETER :: substr='clamstracer_read_nml'
    LOGICAL              :: lex      ! file exists ?
    INTEGER              :: fstat    ! file status
    LOGICAL              :: l_print  ! write control output
    
    
    NAMELIST /CTRL/ timestep_tracer, tracertype, step_perp,  &
                    tracer_desc_file, &
                    lev_bound, file_lsm, file_gph, set_to_zero

    status = 1 !ERROR
    
    if (rank==0 .and. ldiagout) then
       l_print = .true.
    else
       l_print = .false.
    endif
    
    !-------------------------------------------------------------------
    ! Read namelist variables:
    !-------------------------------------------------------------------
    
    CALL read_nml_open(lex, substr, iou, 'CTRL', modstr, l_print)
    IF (.not.lex) RETURN    ! <modstr>.nml does not exist
    
    READ(iou, NML=CTRL, IOSTAT=fstat)
    CALL read_nml_check(fstat, substr, iou, 'CTRL', modstr, l_print)
    IF (fstat /= 0) RETURN  ! error while reading namelist
      
    CALL read_nml_close(substr, iou, modstr, l_print)

    status = 0 !NO ERROR
 
  End Subroutine clamstracer_read_nml
   
  !**************************************************************************
  !
  !**************************************************************************
  Subroutine clamstracer_read_tracerlist (status, iou, TRACERS)

    USE messy_clamstracer_global
    USE messy_clams_global,         ONLY: rank, ldiagout, maxspec
    USE messy_clams_tools_dateconv, ONLY: ymds2js

    IMPLICIT NONE

    !I/O
    INTEGER, INTENT(OUT) ::   status
    INTEGER, INTENT(IN)  ::   iou
    type(tracer_type),  dimension(:), pointer :: TRACERS 
    

    !LOCAL
    CHARACTER(LEN=*),PARAMETER :: substr='clamstracer_read_tracerlist'
    CHARACTER(1200)       :: line
    character(4)         :: yy_str
    character(2)         :: mm_str, dd_str, ss_str
    INTEGER              :: ios, pos
    INTEGER              :: itracer
    INTEGER              :: yy, mm, dd, ss

    !-------------------------------------------------------------------
    ! Read tracer description
    !-------------------------------------------------------------------

    ntracer = 0

    if (tracer_desc_file/='') then
       
       if (rank==0) write (*,*) 'file with list of tracers:',tracer_desc_file
       ! open file with list of tracers
       open (iou,file=tracer_desc_file,status="OLD",iostat=ios)
       if (ios /= 0) then 
          status = 20
          call error_handler(status)
          RETURN
       endif
       
       ! ignore comment lines 
       line = '!'
       do while (line(1:1)=='!')
          read (iou,'(A)',iostat=ios) line
          if (ios/=0) exit
          line = TRIM(ADJUSTL(line))
       enddo
       !write (*,*) 'line: ',line
       
       ! list of tracers
       itracer = 1
       do while (ios==0 .and. line/='' .and. itracer<=maxspec) 
          
          ! name of tracer
          read (line,*,iostat=ios) TRACERS(itracer)%name
          if (ios/=0) then
             status = 21
             call error_handler(status)
             RETURN
          endif
          if (rank==0) write (*,*) 'tracer: ',tracers(itracer)%name
          pos = index(line,' ')
          line = adjustl(line(pos+1:))
          
          if (tracertype=='map') then  ! read segment number
            
             ! segment number
             read (line,*,iostat=ios) TRACERS(itracer)%seg_no
             if (ios/=0) then
                status = 22
                call error_handler(status)
                RETURN
             endif
             pos = index(line,' ')
             line = adjustl(line(pos+1:))
            
          elseif (tracertype=='grid') then ! read lat-,lon-,level-range
            
             ! lower latitude-range
             read (line,*,iostat=ios) TRACERS(itracer)%latmin
             if (ios/=0) then
                status = 23
                print*, TRACERS(itracer)%name
                call error_handler(status)
                RETURN
             endif
             pos = index(line,' ')
             line = adjustl(line(pos+1:))
             ! upper latitude-range
             read (line,*,iostat=ios) TRACERS(itracer)%latmax
             if (ios/=0) then
                status = 24
                call error_handler(status)
                RETURN
             endif
             pos = index(line,' ')
             line = adjustl(line(pos+1:))
             ! lower longitude-range
             read (line,*,iostat=ios) TRACERS(itracer)%lonmin
             if (ios/=0) then
                status = 25
                call error_handler(status)
                RETURN
             endif
             pos = index(line,' ')
             line = adjustl(line(pos+1:))
             ! upper longitude-range
             read (line,*,iostat=ios) TRACERS(itracer)%lonmax
             if (ios/=0) then
                status = 26
                call error_handler(status)
                RETURN
             endif
             pos = index(line,' ')
             line = adjustl(line(pos+1:))
             ! lower level-range
             read (line,*,iostat=ios) TRACERS(itracer)%levmin
             if (ios/=0) then
                status = 27
                call error_handler(status)
                RETURN
             endif
             pos = index(line,' ')
             line = adjustl(line(pos+1:))
             ! upper level-range
             read (line,*,iostat=ios) TRACERS(itracer)%levmax
             if (ios/=0) then
                status = 28
                call error_handler(status)
                RETURN
             endif
             pos = index(line,' ')
             line = adjustl(line(pos+1:))

             if (rank==0) write (*,'(A,6F9.2)') '     latmin,latmax,lonmin,lonmax,levmin,levmax:', &
                                      TRACERS(itracer)%latmin, TRACERS(itracer)%latmax, &
                                      TRACERS(itracer)%lonmin, TRACERS(itracer)%lonmax, &
                                      TRACERS(itracer)%levmin, TRACERS(itracer)%levmax
             
          endif
         
          ! pulse-date
          read (line,*,iostat=ios) TRACERS(itracer)%pulse_date
          if (ios/=0) then
             status = 29
             call error_handler(status)
             RETURN
          endif
          yy_str = trim(TRACERS(itracer)%pulse_date(1:4))
          read(yy_str,*,iostat=ios) yy
          if (ios/=0) then
             status = 41
             call error_handler(status)
             RETURN
          endif
          mm_str = trim(TRACERS(itracer)%pulse_date(5:6))
          read(mm_str,*,iostat=ios) mm
          if (ios/=0) then
             status = 42
             call error_handler(status)
             RETURN
          endif
          dd_str = trim(TRACERS(itracer)%pulse_date(7:8))
          read(dd_str,*,iostat=ios) dd
          if (ios/=0) then
             status = 43
             call error_handler(status)
             RETURN
          endif
          ss_str = trim(TRACERS(itracer)%pulse_date(9:10))
          read(ss_str,*,iostat=ios) ss
          if (ios/=0) then
             status = 44
             call error_handler(status)
             RETURN
          endif
          ss = ss*3600.
          TRACERS(itracer)%pulse_time = ymds2js (yy,mm,dd,ss)
          pos = index(line,' ')
          line = adjustl(line(pos+1:))
          ! pulse length [number of days]
          read (line,*,iostat=ios) TRACERS(itracer)%pulse_length
          if (ios/=0) then
             status = 30
             call error_handler(status)
             RETURN
          endif
          pos = index(line,' ')
          line = adjustl(line(pos+1:))
          ! pulse reset period [years]: reset pulse after this time
          read (line,*,iostat=ios) TRACERS(itracer)%pulse_reset_period
          if (ios/=0) then
             status = 31
             call error_handler(status)
             RETURN
          endif
          pos = index(line,' ')
          line = adjustl(line(pos+1:))

          if (rank==0) write (*,*) '    pulse date/lenth/reset-period: ',&
                                   TRACERS(itracer)%pulse_date, &
                                   TRACERS(itracer)%pulse_length, &
                                   TRACERS(itracer)%pulse_reset_period
          
          ! const, lin, quad, tau
          read (line,*,iostat=ios) TRACERS(itracer)%const
          if (ios/=0) then
             status = 32
             call error_handler(status)
             RETURN
          endif
          pos = index(line,' ')
          line = adjustl(line(pos+1:))
          read (line,*,iostat=ios) TRACERS(itracer)%lin
          if (ios/=0) then
             status = 33
             call error_handler(status)
             RETURN
          endif
          pos = index(line,' ')
          line = adjustl(line(pos+1:))
          read (line,*,iostat=ios) TRACERS(itracer)%quad
          if (ios/=0) then
             status = 34
             call error_handler(status)
             RETURN
          endif
          pos = index(line,' ')
          line = adjustl(line(pos+1:))
          read (line,*,iostat=ios) TRACERS(itracer)%tau
          if (ios/=0) then
             status = 35
             call error_handler(status)
             RETURN
          endif

          if (rank==0) write (*,*) '    const, lin, quad, tau: ', &
                                   TRACERS(itracer)%const, TRACERS(itracer)%lin, &
                                   TRACERS(itracer)%quad, TRACERS(itracer)%tau
          
          ! period, amplitude, phase, type_pulsing
          pos = index(line,' ')
          line = adjustl(line(pos+1:))
          read (line,*,iostat=ios) TRACERS(itracer)%periode
          if (ios/=0) TRACERS(itracer)%periode = 0.
!!$          if (ios/=0) then
!!$             call error_handler(51)
!!$             RETURN
!!$          endif
          pos = index(line,' ')
          line = adjustl(line(pos+1:))
          read (line,*,iostat=ios) TRACERS(itracer)%amplitude
          if (ios/=0)  TRACERS(itracer)%amplitude = 0.
!!$          if (ios/=0) then
!!$             call error_handler(52)
!!$             RETURN
!!$          endif
          pos = index(line,' ')
          line = adjustl(line(pos+1:))
          read (line,*,iostat=ios) TRACERS(itracer)%phase
          if (ios/=0) TRACERS(itracer)%phase = 0.
!!$          if (ios/=0) then
!!$             call error_handler(53)
!!$             RETURN
!!$          endif
          pos = index(line,' ')
          line = adjustl(line(pos+1:))
          read (line,*,iostat=ios) TRACERS(itracer)%type_pulsing
          if (ios/=0)  TRACERS(itracer)%type_pulsing = 0
!!$          if (ios/=0) then
!!$             call error_handler(54)
!!$             RETURN
!!$          endif

          if (rank==0) write (*,*) '    period, amplitude, phase, type_pulsing:', &
                                   TRACERS(itracer)%periode, TRACERS(itracer)%amplitude, &
                                   TRACERS(itracer)%phase, TRACERS(itracer)%type_pulsing
          
                                   
          !!!!!!!    ALJP 190318 add secondary tracer
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          pos = index(line,' ')
          line = adjustl(line(pos+1:))
          read (line,*,iostat=ios) TRACERS(itracer)%parent_specie_name
          !,iostat=ios
          if (ios==-1) then
             TRACERS(itracer)%parent_specie_name=''
          elseif (ios==0) then
              pos = index(line,' ')
              line = adjustl(line(pos+1:))
          else
             call error_handler(36)
          end if

          pos = index(line,' ')
          line = adjustl(line(pos+1:))

!!$          if (TRACERS(itracer)%const==0. .and. TRACERS(itracer)%lin==0. .and. &
!!$              TRACERS(itracer)%quad==0. .and. TRACERS(itracer)%tau==0.) then
          if (TRACERS(itracer)%const==0. .and. TRACERS(itracer)%lin==0) then
             !status = 36
             write(*,*) 'WARNING no pulsing' 
             !call error_handler (status)
             !RETURN
          endif
         

          ! go to next line
          itracer = itracer + 1
          read (iou,'(A)',iostat=ios) line
          !write (*,*) 'line: ',line
          line = adjustl(line)
          pos = index(line,'!')
          if (pos /= 0)  line = line(1:pos-1)
          
          ntracer = ntracer + 1

       end do

       if (rank==0) then
          write (*,*) 'number of tracers: ',ntracer
       endif

       close (iou)

    endif

  End Subroutine clamstracer_read_tracerlist




End Module messy_clamstracer

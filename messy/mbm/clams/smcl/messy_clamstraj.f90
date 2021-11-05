! **************************************************************************
!
! This file is part of the Chemical Lagrangian Model of the
! Stratosphere (CLaMS). CLaMS is a hierarchy of numerical models and
! preprocessors which simulate transport, chemical reactions and
! mixing processes in the stratosphere.  
!
! Copyright (c) 2006
! Richard Swinbank, Daniel McKenna, Nicole Thomas, Paul Konopka, 
! Daniela Ciupa, Juergen Ankenbrand, Barbara Deutsch
! Forschungszentrum Juelich GmbH
! Last Modified By: N.Thomas
! Last Modified On: Wed Jun 10 11:20:32 2020
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
! MODULE FOR CLaMS TRAJECTORY MODULE 'TRAJ'
!********************************************************************************!
! THREE DIMENSIONAL ADVECTION                                                    !
! Runs advection backwards or forwards on log level coordinates                  !
! Takes as input winds on isentropic surfaces                                    !
!   and initial positions (lons, lats, levs, time_init) of n air parcels         !
! Outputs longs and lats in degrees, levs in K or hPa, PV and temperature        !
! Uses netcdf- files as input and output files                                   !
!********************************************************************************!
!
! Module messy_clamstraj contains the following subroutines:
!
! SUBROUTINE TRAJ (status, plats, plons, plevs, julsec, param) 
! SUBROUTINE clamstraj_read_nml(status, iou)
! SUBROUTINE error_handler (error)
! SUBROUTINE read_configuration_messy (status, plats, plons, plevs, julsec) 
! SUBROUTINE increase_dt (dt3)
!
!********************************************************************************!

MODULE messy_clamstraj

  USE messy_clams_global, ONLY: DP

  IMPLICIT NONE

  PUBLIC :: clamstraj_read_nml
  PUBLIC :: TRAJ
  PRIVATE:: read_configuration_messy

  CHARACTER(LEN=*),PARAMETER, PUBLIC :: modstr = 'clamstraj'
  CHARACTER(LEN=*),PARAMETER, PUBLIC :: modver = '1.0'

  CONTAINS 

!**************************************************************************
!
!**************************************************************************
SUBROUTINE TRAJ (status, plats, plons, plevs, julsec, param, lcoupled) 
  
  USE netcdf

  USE messy_clams_global,       ONLY: YEAR, MONTH, DAY, HOUR, MINUTE, SECOND, &
                                      delta_time

  USE messy_clams_global,       ONLY: rank, prec, dp, mdi, eps, &
                                      dnparts, dnparts_max, &
                                      init_vertcoorname, dates30, &
                                      longrid, nx, ny, &
                                      paramtype
  USE messy_clamstraj_global,   ONLY: pslat
  USE messy_clams_tools_utils,  ONLY: lowercase
  USE messy_clamstraj_data,     ONLY: deallocate_arrays
  USE messy_clamstraj_timestep, ONLY: TIMSTP
  
  IMPLICIT NONE

! TEST1:
  REAL(PREC) :: plats(dnparts_max), plons(dnparts_max), plevs(dnparts_max)
  REAL(DP)   :: julsec(dnparts_max)
! TEST2:
!!$  REAL(PREC) :: plats(:), plons(:), plevs(:)
!!$  REAL(DP)   :: julsec(:)
! TEST3:
!!$  REAL(PREC), DIMENSION(:), POINTER :: plats(:), plons(:), plevs(:)
!!$  REAL(DP),   DIMENSION(:), POINTER :: julsec(:)

  TYPE(paramtype), DIMENSION(:), POINTER :: param

  INTEGER :: status

  LOGICAL :: lcoupled

  ! Local variables:
  REAL(PREC) :: pi
  REAL(DP)   :: seconds  
  
  INTEGER :: i
  INTEGER :: finish, start, counts_per_sec
  
  status = 0

  ! some constants
  pi=4.*atan(1.)                                                        
  pslat=pi*.4       ! latitude to switch to polar stereo graphic grid

  if(rank==0)THEN
     WRITE (*,*) 
     WRITE (*,*) 'START OF TRAJ'
     WRITE (*,*) 
     write(*,'(A,I8,I4,I4,i4,i4,i6)') 'Current date:', YEAR, MONTH,DAY,HOUR, MINUTE, SECOND
  ENDIF

  CALL system_clock (count=start)

  !--------------------------------------------------------------
  ! Read input parameters from configuration file 
  !--------------------------------------------------------------

  CALL read_configuration_messy (status, plats, plons, plevs, julsec)
  if (status/=0) then
     CALL error_handler (status)
     return
  endif
 
  if(rank==0)THEN
     if (dates30) then
        write (*,*)
        write (*,*) 'ACHTUNG: 30-Tage-Monate !!!'
        write (*,*)
     endif
  endif

  IF(rank==0)THEN
     WRITE (*,*) 'Latitude to switch to polar stereo graphic grid: ',(pslat*180.)/pi
     WRITE (*,*) 
  ENDIF
  
 
   !--------------------------------------------------
   ! Run time loop
   !--------------------------------------------------
   CALL TIMSTP (status, plats, plons, plevs, julsec, param, lcoupled)          
   if (status/=0) return
   
   CALL deallocate_arrays
    
   IF(rank==0)THEN

      WRITE (*,*) 'Normal termination of traj'  
                  
      CALL system_clock (COUNT_RATE=counts_per_sec)
      CALL system_clock (count=finish)
      IF (finish > start) THEN
         seconds = float(finish-start) / float(counts_per_sec)
      ELSE
         seconds = 0.0
      ENDIF
      
      !WRITE (*,*)
      !WRITE (*,*) 'System clock runs at ', counts_per_sec, 'ticks per second'
      WRITE (*,*)
      WRITE (*,'(A,F10.2,A)') 'This job has taken ',seconds,' seconds to execute.' 

   ENDIF
  
   ! convert from polar-stereographic to latitudes and longitudes
   DO i = 1, dnparts
    IF (abs((plons(i)-MDI)/MDI)<=eps) THEN
        plons(i) = MDI
        plats(i) = MDI
        plevs(i) = MDI
     ELSE
! op_pj_20160606+
!!$        plons(i) = mod(plons(i)*180./pi,360.) 
        plons(i) = mod(plons(i)*180._prec/pi,360._prec) 
! op_pj_20160606-
        plats(i) = plats(i)*180./pi                           
        plevs(i) = plevs(i)
     ENDIF
  ENDDO

!!$  IF (TRIM(init_vertcoorname) == 'press') THEN 
!!$     pressure = plevs
!!$  END IF

 END SUBROUTINE TRAJ

!--------------------------------------------------------------------

!--------------------------------------------------------------------
    SUBROUTINE clamstraj_read_nml(status, iou)

      USE messy_main_tools, ONLY: read_nml_open, read_nml_check, read_nml_close
      USE messy_clams_global, ONLY: rank, ldiagout
      USE messy_clamstraj_global

      IMPLICIT NONE

      !I/O
      INTEGER, INTENT(OUT) ::   status
      INTEGER, INTENT(IN)  ::   iou
      !LOCAL
      CHARACTER(LEN=*),PARAMETER :: substr='clamstraj_read_nml'
      LOGICAL              :: lex      ! file exists ?
      INTEGER              :: fstat    ! file status
      LOGICAL              :: l_print  ! write control output

      NAMELIST /CTRL/ type_traj, p_ref, timestep_trajout,  &
                      timestep_trajin, dir_trajin

      status = 1 !ERROR

      if (rank==0 .and. ldiagout) then
         l_print = .true.
      else
         l_print = .false.
      endif

      ! Read namelist variables:
      CALL read_nml_open(lex, substr, iou, 'CTRL', modstr, l_print)
      IF (.not.lex) RETURN    ! <modstr>.nml does not exist

      READ(iou, NML=CTRL, IOSTAT=fstat)
      CALL read_nml_check(fstat, substr, iou, 'CTRL', modstr, l_print)
      IF (fstat /= 0) RETURN  ! error while reading namelist

      CALL read_nml_close(substr, iou, modstr, l_print)

      status = 0 !NO ERROR
      
    END SUBROUTINE clamstraj_read_nml

!--------------------------------------------------------------------

!--------------------------------------------------------------------
SUBROUTINE error_handler (error)

  IMPLICIT NONE

  INTEGER :: error

  !WRITE (*,*)
  WRITE (*,*) '***********************************************************'

  SELECT CASE (error)
  CASE (9)
     WRITE (*,*) '    trajtype must be "1","2" or "3" !!!'
  CASE (10)
     WRITE (*,*) '    end date incorrect'             
  CASE (11)
     WRITE (*,*) '    For backward trajectories end date ', &
          'must be earlier than start date !'
  CASE (12)
     WRITE (*,*) '    For forward trajectories end date ', &
          'must be later than start date !'
  CASE (13)
     WRITE (*,*) '    Start date and end date are the same !'             
  CASE (15)
     WRITE (*,*) '    timestep1 is incorrect !!!'
  CASE (16)
     WRITE (*,*) '    timestep2 is incorrect !!!'
  CASE (17)
     WRITE (*,*) '    (timestep1 MOD timestep2) <> 0 !!!'
  CASE (18)
     WRITE (*,*) '    Writing frequency is incorrect !!!'
  CASE (21)
     WRITE (*,*) '    short name of parameter could not be read!!!'
  CASE (28)
     WRITE (*,*) 'For forward trajectories end date ', &
                 'must be later than all starting times'
  CASE (29)
     WRITE (*,*) 'For backward trajectories end date ', &
                 'must be earlier than all starting times'

  CASE (50)
     WRITE (*,*) 'Initial positions could not be read !!!'
  CASE (51)
     WRITE (*,*) 'Level could not be read !!!'
  CASE (52)
     WRITE (*,*) 'Output file could not be created ! '
  CASE (53)
     WRITE (*,*) 'Output file could not be opened ! '
  CASE (54)
     WRITE (*,*) 'latitudes are not in valid range !!!'
  CASE (55)
     WRITE (*,*) 'longitudes are not in valid range !!!'
  CASE (56)
     WRITE (*,*) 'Wind file could not be opened ! '

  END SELECT

END SUBROUTINE error_handler

!--------------------------------------------------------------------

!--------------------------------------------------------------------
SUBROUTINE read_configuration_messy (status, plats, plons, plevs, julsec) 

   USE messy_clams_global,         ONLY: YEAR, MONTH, DAY, HOUR, &
                                         MINUTE, SECOND, delta_time

   USE messy_clams_global,         ONLY: prec, dp, mdi, eps, rank, &
                                         dnparts, dnparts_max,  &
                                         dates30,  &
                                         init_vertcoorname,  &
                                         levelno
   USE messy_clamstraj_global
   USE messy_clams_tools_dateconv, ONLY: ymds2js_interface, js2ymds_interface, &
                                         d1gtd2, check_date, check_date30
   USE messy_clamstraj_data,       ONLY: allocate_arrays

   IMPLICIT NONE

   INTEGER     :: status

   REAL(PREC)  :: plats(dnparts_max), plons(dnparts_max), plevs(dnparts_max)
   REAL(DP)    :: julsec(dnparts_max)

   ! endcalc : endtime of calculation (in julian seconds)
   REAL(DP)    ::  endcalc, endjultime

   REAL(PREC) :: pi

   INTEGER :: idmmin, idmhr, number,  pos, i 
             

   LOGICAL ::  ok

   real(prec),   dimension(:), allocatable :: helparr

   REAL(DP)         :: start_t_messy, end_t_messy

   INTEGER :: starttime_min, endtime_min, starttime_hour, endtime_hour

   pi=4.*atan(1.)                                                        
   
   status = 0
   number = 0

   !--------------------------------------------------------------
   ! Backward or forward trajectories 
   !--------------------------------------------------------------

! currently only forward trajectories in MESSy
   forward = .TRUE.

   !-------------------------------------------------------------
   ! Read initial positions and starttimes
   !-------------------------------------------------------------

   ! -90<=plats<=90  (oder mdi)
   allocate (helparr(dnparts_max))
   helparr = plats
   where (ABS((plats-mdi)/mdi)<=eps)
      helparr = 0.
      plons = mdi
      plevs = mdi
   end where
   if (minval(helparr)<-90 .or. maxval(helparr)>90.) then
      status = 54 
      return
   endif
   deallocate (helparr)

   ! 0<=plons<=360   (or mdi)  
   where (ABS((plons-mdi)/mdi)<=eps)
      plats = mdi
      plevs = mdi
   end where
   do i = 1, dnparts
      if (abs((plons(i)-mdi)/mdi)>eps) then
         if (plons(i)>360.) then                               
! op_pj_20160606+
!!$            plons(i) = mod(plons(i),360.)
            plons(i) = mod(plons(i),360._prec)
! op_pj_20160606-
         else if (plons(i)<0.) then    
            do while (plons(i)<0.)
               plons(i) = plons(i) + 360.      
            end do
         endif
      endif
   enddo

   where (ABS((plevs-mdi)/mdi)<=eps)
      plats = mdi
      plons = mdi
   end where
   where (plevs<=0.)
      plats = mdi
      plons = mdi
   end where

   CALL allocate_arrays

   levelno=1

   !-------------------------------------------------------------
   ! Change units
   !-------------------------------------------------------------

   where (ABS((plats-mdi)/mdi)>eps)
      plons(:)=plons(:)*(pi/180.)
      plats(:)=plats(:)*(pi/180.)
   end where

   !--------------------------------------------------------------
   ! Start date for calculation
   !--------------------------------------------------------------
!!!!! ???
!All trajectories start at the same time:
   loop_julsec: do i=1,dnparts_max
      if (ABS((julsec(i)-mdi)/mdi)>eps) then
         startfirsttraj = julsec(i)
         startlasttraj  = julsec(i)
         exit loop_julsec
      end if
   end do loop_julsec

!    IF (forward) THEN

!!!!!! Messy starttime
   starttime%year = YEAR
   starttime%month = MONTH
   starttime%day = DAY
   starttime_hour = HOUR
   starttime_min = MINUTE
   starttime%sec = SECOND + MINUTE*60 + HOUR*3600
  

   IF(rank==0)THEN
      WRITE (*,'(A,I4,1X,4(I2.2,1X))') 'Start date for calculation (yyyy mm dd hh mm): ', &
           starttime%year, starttime%month, starttime%day, starttime%sec/3600,&
           & starttime_min
   ENDIF

   !--------------------------------------------------------------
   ! End date for calculation
   !--------------------------------------------------------------

   IF (startfirsttraj == startlasttraj) THEN
      trajtype = 1
   ELSE
      number = type_traj
      SELECT CASE (number)
      CASE (1)
         trajtype = 2
      CASE (2)
         trajtype = 3
      CASE DEFAULT
         status = 9
         return
      END SELECT
   ENDIF

   IF ((trajtype == 1) .OR. (trajtype == 2)) THEN
     ! Calculate endtime from MESSy time & MESSy timestep 
     start_t_messy = ymds2js_interface (YEAR, MONTH, DAY, HOUR*3600+MINUTE*60+SECOND,dates30)
     end_t_messy = start_t_messy + delta_time
     CALL js2ymds_interface (end_t_messy, endtime%year, endtime%month, endtime&
          &%day, endtime%sec, dates30)     
     endtime_hour =  endtime%sec / 3600.

     if (dates30) then
        CALL check_date30 (endtime_hour, endtime%day, endtime%month, endtime%year, ok)
     else
        CALL check_date (endtime_hour, endtime%day, endtime%month, endtime%year, ok)
     endif
     IF (ok) THEN
        IF (d1gtd2(endtime%sec,endtime%day,endtime%month,endtime%year,  &
             starttime%sec,starttime%day,starttime%month,starttime%year)) THEN
           IF (.NOT. forward) then
              status = 11
              return
           endif
        ELSEIF (d1gtd2(starttime%sec,starttime%day,starttime%month, &
             starttime%year,  &
             endtime%sec,endtime%day,endtime%month,endtime%year)) THEN    
           IF (forward) then
              status = 12
              return
           ENDIF
        ELSE 
           status = 13
           return
        ENDIF
        ! check, if no trajectory starts after (backward: before) endtime
        endjultime = ymds2js_interface (endtime%year,endtime%month,endtime%day,endtime%sec,dates30)

     ELSE
        status = 10
        return
     ENDIF
  
   ELSE  
      !READ (3,*,IOSTAT=ios) nrhours   !  Enter trajectory lenght (hours)
      !IF (ios /= 0)  call error_handler(14)
   ENDIF

   !---------------------------------------------------
   ! Timestep
   !---------------------------------------------------

   idtsec = delta_time

   IF (idtsec == 0) THEN
      status = 15
      return
   ELSEIF ((idtsec<60) .OR. (idtsec>3600)) THEN
      status = 15
      return
   ELSEIF (forward .AND. idtsec<0) THEN
      idtsec = -1*idtsec
   ELSEIF (.NOT. forward .AND. idtsec>0) THEN
      idtsec = -1*idtsec
   ENDIF

   dt = idtsec

!!$   IF (idtsec2 == 0) THEN
!!$      status = 16
!!$      return
!!$   ELSEIF (MOD(idtsec,idtsec2) /= 0) THEN
!!$      status = 17
!!$      return
!!$   ELSEIF (forward .AND. idtsec2<0) THEN
!!$      idtsec2 = -1*idtsec2
!!$   ELSEIF (.NOT. forward .AND. idtsec2>0) THEN
!!$      idtsec2 = -1*idtsec2
!!$   ENDIF
!!$
!!$   dt2 = idtsec2
   dt2 = 60 ! wird immer Moment nicht genutzt !

   !---------------------------------------------------
   !
   !---------------------------------------------------

   endtime_hour = endtime%sec/3600
   endtime_min = (endtime%sec - endtime_hour*3600) /60

   IF(rank==0)THEN
      WRITE (*,'(A,I4,1X,4(I2.2,1X))') 'End date for calculation (yyyy mm dd hh mm): ', &
           endtime%year, endtime%month, endtime%day, endtime%sec/3600, endtime_min
      WRITE (*,*)
   ENDIF

   !--------------------------------------------------------------
   ! Writing frequency
   !--------------------------------------------------------------
   ! Return results to MESSy after every timestep: 
   idmsec = idtsec
   idmday = idtsec / 86400
   idmhr  = (idtsec - idmday*86400) / 3600
   idmmin = (idtsec - idmday*86400 - idmhr*3600) / 60

   IF (idmmin<0 .OR. idmmin>=60 .OR. idmhr<0 .OR. idmhr>=24 .OR. idmday<0 &
                .OR. (idmmin==0 .AND. idmhr==0 .AND. idmday==0)) THEN
      status = 18
      return
   ENDIF

   idmsec = idmmin*60 + 3600*idmhr                                                   
   IF (.NOT. forward) THEN
      IF (idmsec > 0.) idmsec=-idmsec                                   
      IF (idmday > 0.) idmday=-idmday                                   
   ENDIF

END SUBROUTINE read_configuration_messy



!--------------------------------------------------------------------

!--------------------------------------------------------------------
  SUBROUTINE increase_dt (dt3,lcoupled)

  USE messy_clams_global,         ONLY: delta_time, YEAR, MONTH, DAY, &
                                        HOUR, MINUTE, SECOND

  USE messy_clams_global,         ONLY: prec, dp, mdi, eps,  &
                                        nx, ny, nz, &
                                        leveldt, levelfut, &
                                        fut_year, fut_month, fut_day, fut_sec, &
                                        udt, vdt, wdt, ufut, vfut, wfut
  USE messy_clamstraj_global
  USE messy_clams_tools_dateconv, ONLY: ymds2js_interface, datsec,&
                                        incdat_interface,ymds2js, js2ymds_interface
  USE messy_clams_global,         ONLY: dnparts_max,rank, timetype, dates30

  IMPLICIT NONE

  REAL(PREC)     :: dt3
  LOGICAL        :: lcoupled
  TYPE(timetype) :: itime, ifutime

  ! local variables

  REAL(PREC),dimension(:,:,:),allocatable :: u, v, w, umid, vmid, wmid   
  real(prec),dimension(:,:,:),allocatable :: level, levelmid

  real(prec),dimension(:,:,:),allocatable :: mdi_feld
  real(prec),dimension(:,:,:),allocatable :: hilfsfeld

  logical,dimension(:,:,:),allocatable :: selection 


  REAL(PREC)      :: pi, er 
  REAL(DP)        :: ijulsec, past, future
  INTEGER         :: year_local, month_local, day_local, sec_local, sod

  sod = HOUR*3600 + MINUTE*60 + SECOND
  ijulsec = ymds2js(YEAR, MONTH, DAY, sod)

  CALL js2ymds_interface (ijulsec, year_local, month_local, day_local, sec_local, dates30)
  itime%year = year_local
  itime%month = month_local
  itime%day = day_local
  itime%sec = sec_local

  ifutime%sec   = fut_sec
  ifutime%day   = fut_day
  ifutime%month = fut_month
  ifutime%year  = fut_year

  allocate (u(nx,ny,nz))
  allocate (v(nx,ny,nz))
  allocate (w(nx,ny,nz))
  allocate (level(nx,ny,nz))
  allocate (umid(nx,ny,nz))
  allocate (vmid(nx,ny,nz))
  allocate (wmid(nx,ny,nz))
  allocate (levelmid(nx,ny,nz))
  allocate (mdi_feld(nx,ny,nz))
  allocate (hilfsfeld(nx,ny,nz))
  allocate (selection(nx,ny,nz))
  
  mdi_feld = mdi
  selection = .true.

  er=6.371e6
  pi=4.*atan(1.)                                                        
         
  ! the function datsec returns the date in seconds from 0,1,1,1960               
  past   = datsec(itime%sec,itime%day,itime%month,itime%year)
  if (lcoupled) then
     future = past + delta_time
  else
     future = datsec(ifutime%sec,ifutime%day,ifutime%month,ifutime%year)
  endif

  !-----------------------------------------------------------------
  ! Calculate:
  ! u,v,w            winds at time t       
  ! umid,vmid,wmid   winds at time t+dt3/2       
  ! udt,vdt,wdt      winds at time t+dt3       
  !-----------------------------------------------------------------
!write(*,*) 'increase_dt','past', past
!write(*,*) 'increase_dt','future', future
!write(*,*) 'increase_dt','delta',future-past

  where(ABS((udt-mdi)/mdi)<=eps)
     selection = .false.
  end where
  where(ABS((vdt-mdi)/mdi)<=eps)
     selection = .false.
  end where
  where(ABS((wdt-mdi)/mdi)<=eps)
     selection = .false.
  end where

  where(ABS((ufut-mdi)/mdi)<=eps)
     selection = .false.
  end where
  where(ABS((vfut-mdi)/mdi)<=eps)
     selection = .false.
  end where
  where(ABS((wfut-mdi)/mdi)<=eps)
     selection = .false.
  end where

  u = merge(udt,mdi_feld,selection)
  v = merge(vdt,mdi_feld,selection)
  w = merge(wdt,mdi_feld,selection)
  level = merge(leveldt,mdi_feld,selection)

  hilfsfeld = (0.5*dt3*ufut+(future-past-0.5*dt3)*u)/(future-past)
  umid = merge(hilfsfeld,mdi_feld,selection)
  hilfsfeld = (0.5*dt3*vfut+(future-past-0.5*dt3)*v)/(future-past)
  vmid = merge(hilfsfeld,mdi_feld,selection)
  hilfsfeld = (0.5*dt3*wfut+(future-past-0.5*dt3)*w)/(future-past)
  wmid = merge(hilfsfeld,mdi_feld,selection)
  hilfsfeld = (0.5*dt3*levelfut+(future-past-0.5*dt3)*level)/(future-past)
  levelmid = merge(hilfsfeld,mdi_feld,selection)
  hilfsfeld = (dt3*ufut+(future-past-dt3)*u)/(future-past)
  udt = merge(hilfsfeld,mdi_feld,selection)
  hilfsfeld = (dt3*vfut+(future-past-dt3)*v)/(future-past)
  vdt = merge(hilfsfeld,mdi_feld,selection)
  hilfsfeld = (dt3*wfut+(future-past-dt3)*w)/(future-past)
  wdt = merge(hilfsfeld,mdi_feld,selection)
  hilfsfeld = (dt3*levelfut+(future-past-dt3)*level)/(future-past)
  leveldt = merge(hilfsfeld,mdi_feld,selection)

  deallocate (u, v, w, level)
  deallocate (umid, vmid, wmid, levelmid)
  deallocate (mdi_feld)
  deallocate (hilfsfeld)
  deallocate (selection)

  END SUBROUTINE increase_dt

! **************************************************************************
END MODULE messy_clamstraj
! **************************************************************************

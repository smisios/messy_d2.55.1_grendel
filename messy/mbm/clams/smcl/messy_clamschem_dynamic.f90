MODULE messy_clamschem_dynamic
!*******************************************************************************
!
!*******************************************************************************

  USE messy_clams_global, ONLY: rank, dnparts_max, dnparts

!*******************************************************************************
CONTAINS
!*******************************************************************************

subroutine indynam(status, prdim)

  ! MESSy Main
  use messy_main_constants_mem, only: DP
  USE messy_clams_global,       ONLY: YEAR_NEXT, MONTH_NEXT, DAY_NEXT,  &
                                      HOUR_NEXT, MINUTE_NEXT, SECOND_NEXT

  ! ASAD
  USE messy_clamschem_asad_mod,        only: dtime
  USE messy_clamschem_asad_mod_clams,  only: failed
  
  ! CLaMS
  USE messy_clamschem_global,     ONLY: maxtraj, ntraj_rank, ntraj, &
                                        maxpacks, npacks, jpnl_count, jpnl_offset, &
                                        dsn_twodavg, ncdt, twodavgfile, iodump, &
                                        js
  USE messy_clams_tools_dateconv, ONLY: ymds2js, js2ymds

  implicit none
  
  !-----------------------------------------------
  !   D u m m y   A r g u m e n t s
  !-----------------------------------------------
  integer  :: status
  real(DP) , intent(out) :: prdim 
  !-----------------------------------------------
  !   L o c a l   V a r i a b l e s
  !-----------------------------------------------
  integer   :: i1,i2,i3,i 
  integer , dimension(7) :: it, if 
  integer   :: itrajtype
  real(DP)  :: jsec_i, jsec_f 
  character :: chunit*8
  INTEGER   :: sod_i, sod_f
  REAL(DP)  ::julsec
  integer   :: rcode
  !-----------------------------------------------
  !
  if (iodump .AND. rank == 0) write (6, *) 'CALLING INDYNAM' 

  status = 0 ! no error 
  
  twodavgfile= dsn_twodavg
  

!!!!! ???
!  ntraj_rank = dnparts_max
  ntraj_rank = dnparts


  ! get number of packages of current rank -> npacks
  if (ntraj_rank <= maxtraj) then
     npacks = 1
  else
     npacks = ntraj_rank/maxtraj + 1
     IF  (ntraj_rank/npacks+ modulo(ntraj_rank,npacks)>maxtraj) npacks=npacks+1
     if (npacks > maxpacks) then
        write(*,*) 'too much airparcels: npacks =',npacks,'; maxpacks =',maxpacks
        write(*,*) 'maxpacks must be increased in gmodules.f90 !'
        status = -111
        return
     endif
  endif

  ! get start index and count of all packages -> jpnl_offset(1:npacks), jpnl_count(1:npacks)
  do i = 1, npacks
     jpnl_count(i)  = ntraj_rank/npacks + (i/npacks * modulo(ntraj_rank,npacks))
     jpnl_offset(i) = sum(jpnl_count(1:i-1)) + 1
  enddo
  if (iodump) then
     write (*,*) 'rank, ntraj, npacks:',rank, ntraj_rank, npacks
     write (*,*) 'rank, count:',rank,jpnl_count(1:npacks)
     write (*,*) 'rank, offset:',rank,jpnl_offset(1:npacks)
  endif
  if (minval(jpnl_count(:npacks)) <= 0 .or. minval(jpnl_offset(:npacks)) <= 0) then
     write (*,*) 'jpnl_count<=0 or jpnl_offset<=0 !!!'
     status = -112
     return
  end if


  ! Trajectory type: Forward trajectories, all start and end at the same time
  itrajtype = 1           !! Other trajtypes not yet implemented for MESSy !!


  if=(/YEAR_NEXT,MONTH_NEXT,DAY_NEXT,HOUR_NEXT,MINUTE_NEXT,SECOND_NEXT,0/)

  sod_f = HOUR_NEXT*3600 +MINUTE_NEXT*60 + SECOND_NEXT
  julsec= ymds2js(YEAR_NEXT, MONTH_NEXT, DAY_NEXT, sod_f) 

  jsec_i = julsec - dtime
  jsec_f = julsec

  CALL js2ymds(jsec_i,it(1),it(2),it(3),sod_i)
  
  it(4) = sod_i / 3600
  it(5) = (sod_i - it(4)*3600) / 60
  it(6) = sod_i - it(5)*60 -it(4)*3600
  it(7) = 0

  if (iodump .AND. rank == 0) then 
        write (6, *) 'Initial Trajectory Time', it 
        write (6, *) '  Final Trajectory Time', if 
        write (6, *) 'Initial Trajectory Time [s]', jsec_i 
        write (6, *) '  Final Trajectory Time [s]', jsec_f 
  endif

  js=jsec_i
  
  ! ju_jug_20150827 
  ! convert pressure units to Pa using factor prdim
  ! now assuming that pressure data in clamstraj are given always in hPa
  prdim = 100.

  if (iodump .AND. rank == 0) write (6, *) 'EXITING  INDYNAM' 

  return  

end subroutine indynam  
!*******************************************************************************
 
!*******************************************************************************
subroutine dynamic(status, ichemstep, prdim, &
                   LAT, LON, LEV, TEMPERATURE, PRESSURE, &
                   LAT_OLD, LON_OLD, LEV_OLD, TEMPERATURE_OLD, PRESSURE_OLD) 
  !-----------------------------------------------
  !   M o d u l e s 
  !-----------------------------------------------
 
  ! MESSy Main
  use messy_main_constants_mem, only: DP

  ! ASAD
  use messy_clamschem_asad_mod,       only: p, t, dtime
  use messy_clamschem_asad_mod_clams, only: jpctr, failed
  
  ! CLaMS
  use messy_clamschem_global, only: ipart, missing_value, ntraj, &
                                    jpnl_count, jpnl_offset, iodump, &
                                    zangle, ftr, ftr_ini, missing_index, &
                                    lats, lons, theta, js, temps, press, sza, slt

  use messy_dissoc,              only: calc_zenith
  use messy_clams_tools_utils,   only: lowercase, uppercase
  USE messy_clams_tools_dateconv, ONLY: js2jd

  implicit none

  !-----------------------------------------------
  !   D u m m y   A r g u m e n t s
  !-----------------------------------------------
  integer :: status
  integer , intent(in) :: ichemstep
  real(DP), intent(in) :: prdim
  REAL(DP):: LAT(dnparts_max)
  REAL(DP):: LON(dnparts_max) 
  REAL(DP):: LEV(dnparts_max)
  REAL(DP):: TEMPERATURE(dnparts_max)
  REAL(DP):: PRESSURE(dnparts_max)
  REAL(DP):: LAT_OLD(dnparts_max)
  REAL(DP):: LON_OLD(dnparts_max)
  REAL(DP):: LEV_OLD(dnparts_max)
  REAL(DP):: TEMPERATURE_OLD(dnparts_max)
  REAL(DP):: PRESSURE_OLD(dnparts_max)
  !-----------------------------------------------
  !   L o c a l   V a r i a b l e s
  !-----------------------------------------------
  integer :: start, end
  integer :: i
  real(DP) :: dlon(ntraj),dt,sza0(ntraj)
  real(DP),parameter :: sza_term=95.   ! solar zenith angle of terminator
!  real,parameter :: pressmax=50000.0 ! maximal pressure (Pa) accepted for stratospheric chemistry (SB 29.07.2003)
  real(DP),parameter :: pressmax=110000.0 ! maximal pressure (Pa) accepted for stratospheric chemistry (all normal values JUG 9/2006)
  integer :: timedimid
  real(DP) :: js_next,js_h, jd
  !-----------------------------------------------
  !
  ! set parameter t_offset to calculate with shifted temperatures with respect
  ! to the trajectory temperature
  real(DP),parameter::t_offset=0.0

  status = 0 ! no error

  if (t_offset /=0 .and. rank ==0) write(*,*) 'Warning: t_offset=',t_offset

!$$  jpnl0=ntraj  ! ntraj = dnparts_max

  if (iodump .AND. rank == 0) write (*, *) 'started dynamic',  prdim

  start = jpnl_offset(ipart)
  end = start + jpnl_count(ipart) - 1

  IF (ichemstep == 0) THEN
       
!!!!!
     allocate (failed(ntraj))
     failed = .false.

!!!!!
     allocate (lats(ntraj), lons(ntraj), theta(ntraj))
     allocate (temps(ntraj), press(ntraj))
     allocate (sza(ntraj), slt(ntraj))

!!!!! aus indynam hierher verschoben:
     press   = 1.E4    ! dummy pressure value
  
     temps = 250.   ! dummy temperature value
       
     lats(:ntraj)  = LAT_OLD(start:end)
     lons(:ntraj)  = LON_OLD(start:end)
     press(:ntraj) = PRESSURE_OLD(start:end)
     press(:ntraj) = press(:ntraj)*prdim 
     temps(:ntraj) = TEMPERATURE_OLD(start:end)
     temps(:ntraj) = temps(:ntraj) + t_offset
     theta(:ntraj) = LEV_OLD(start:end)
     sza(:) = 0.
     slt(:) = 0.

     where(abs(lats(:ntraj)-missing_value) <0.1) 
        missing_index=.true.
     elsewhere
        missing_index=.false.
     endwhere
     where(abs(lats-missing_value) <0.1) 
        lats=0.
        lons=0.
        press=1.E4
        temps=250.
        theta=500.
     endwhere

     ! also check pressure and temperature
     where(press(:ntraj) <0.0 .or. press(:ntraj) >pressmax ) missing_index=.true.
     where(temps(:ntraj) <0.0) missing_index=.true.
     where(press <0.0 .or. press >pressmax .or. temps < 0.0) 
        lats=0.
        lons=0.
        press=1.E4
        temps=250.
        theta=500.
     endwhere

  ELSE IF (ichemstep == 1) THEN

     ! get trajectory data from next timestep
 
     lats(:ntraj)  = LAT(start:end)
     lons(:ntraj)  = LON(start:end)
     theta(:ntraj) = LEV(start:end)
     temps(:ntraj) = TEMPERATURE(start:end)
     temps(:ntraj) = temps(:ntraj)+ t_offset 
     press(:ntraj) = PRESSURE(start:end) * prdim 

     where(abs(lats(:ntraj)-missing_value) <0.1) 
        missing_index=.true.
     elsewhere
        missing_index=.false.
     endwhere
     where(abs(lats-missing_value) <0.1) 
        lats=0.
        lons=0.
        press=1.E4
        temps=250.
        theta=500.
     endwhere

     ! also check pressure and temperature
     where(press(:ntraj) <0.0 .or. press(:ntraj) >pressmax ) missing_index=.true.
     where(temps(:ntraj) <0.0) missing_index=.true.
     where(press <0.0 .or. press >pressmax .or. temps < 0.0) 
        lats=0.
        lons=0.
        press=1.E4
        temps=250.
        theta=500.
     endwhere

  ELSE
     write(*,*) 'Dynamic: ichemstep /= 0/1 !!!'
     write(*,*) 'Program will stop.. '
     status = -114
     return
  ENDIF

  call calc_date_time

  ! calculate zenith angle
  ! now use program by Theo Brauers/Klaus Pfeilsticker; formula from Yoshio Kubo :
  ! model prduces coredump for some reason if sza instead of sza0 is used
  jd = js2jd(js)
  call calc_zenith(lats,lons,ntraj,js,jd,sza0)
  sza(:ntraj)=sza0(:ntraj)

 ! calculate the increment of sunlight time SLT
  !<must be corrected for potential missing values...>
  if (ichemstep > 0) then 
      where (sza < sza_term) slt=slt+dtime
  endif
  
  where(failed(1:ntraj)) missing_index=.true.

  p = press 
  t = temps 
  zangle = sza 
  return  
 
end subroutine dynamic 
!*******************************************************************************

!*******************************************************************************
subroutine calc_date_time
  ! converts js to date_time 

  use messy_main_constants_mem, only: DP

  use messy_clams_tools_dateconv
  use messy_clams_global, only: dates30

  use messy_clamschem_global, only: js, date_time

  implicit none

  integer :: iy,im,id,ih,imin,isec,is
  real(DP)    :: rtime
  ! round to seconds
  js=dble(nint(js))
  if (dates30) then
     CALL js30_to_ymds30(js,iy,im,id,is)
  else
     CALL js2ymds(js,iy,im,id,is)
  endif
  rtime=is/86400.
  ih=int(is/3600)
  imin=int((is-ih*3600)/60.)
  isec=is-ih*3600-imin*60
  date_time=(/iy,im,id,ih,imin,isec,0/)

end subroutine calc_date_time
!*******************************************************************************
END MODULE messy_clamschem_dynamic

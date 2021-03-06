Module messy_clamschem_get_const_2d


contains

SUBROUTINE get_const_2d(status,spec,const_arr) 


! Routine that reads zonally+diurnally averaged quantities, currently
! OH, O(1D), and Cl from a climatology (Mainz 2-D model) and
! interpolates it to the latitudes/pressures from the simulation. This
! is thought for longterm simulations in which the diurnal variability
! of the chemical species does not matter. The chemical species are
! handled as constant variables. This routine should be called from 
! fyinit in ftoy.
!
! 21.09.2006 Jens-Uwe Grooss



  USE messy_main_constants_mem, ONLY: DP

  use netcdf

  USE messy_clams_global,      ONLY: rank, buffersize
  USE messy_clamschem_global,  ONLY: ntraj, &
                                     twod_nlats, twod_lats, twod_press, &
                                     twod_avg, twod_avg_all, twod_lnp, &
                                     twod_npress, twodspecs, twodavgfile, &
                                     twod_data_read, nspecs, &
                                     date_time, lats, press
  use messy_clams_tools_utils, only: delete_control_chars

  implicit none


  INTEGER :: status
  REAL(kind=DP),intent(out),dimension(*):: const_arr
  real(kind=DP):: lat0,dlat,rl,pr,dl,dpr,x1,x2,y1,y2,dt
  CHARACTER*10,intent(in):: spec
  CHARACTER*10:: spec0
  INTEGER::ncid,varid,i,i1,i2,jl,k1,k2,month,m1,m2,day,k(1),ispec
  INTEGER,dimension(twod_nlats):: ilat

  status = 0 ! no error

  ! file from which zonally+diurnally  averaged quantities can be read
  ! will be read in from chem.inp

  ! get month and day for time interpolation
  month=date_time(2)
  day=date_time(3)
  
  if (day > 15) then
     dt=min((day-15)/30.,0.5)
     m1 = month
  else
     dt=0.5+day/30.
     m1 = month-1
  endif

  m2=m1+1
  if (m1==0) m1=12
  if (m2==13) m2=1

  spec0=delete_control_chars(spec)

!  write(*,*) 'started get_const_2d ',month,day 

  select case (trim(spec0))
  case ('OH','Cl','O(1D)','HO2')

    if (.not. twod_data_read) then
      ! read in the 2-D climatology for the chosen species, latitude and pressure
      if (rank==0) write(*,*) 'reading '//twodavgfile
      do ispec=1,nspecs
        status = nf90_open(twodavgfile, nf90_nowrite, ncid, buffersize)
        if (status /= 0) then
           write(*,*)'get_const_2d 1:',nf90_strerror(status)
           return
        endif

        status = nf90_inq_varid(ncid,twodspecs(ispec),varid)
        if (status /= 0) then
           write(*,*)'get_const_2d 2:',nf90_strerror(status),twodspecs(ispec)
           return
        endif
        status = nf90_get_var(ncid, varid, twod_avg)
        if (status /= 0) then
           write(*,*)'get_const_2d 3:',nf90_strerror(status)
           return
        endif
        twod_avg_all(:,:,:,ispec)=twod_avg
      enddo
      status = nf90_inq_varid(ncid,'lat',varid)
      if (status /= 0) then
         write(*,*)'get_const_2d 4:',nf90_strerror(status)
         return
      endif
      status = nf90_get_var(ncid, varid, ilat)
      if (status /= 0) then
         write(*,*)'get_const_2d 5:',nf90_strerror(status)
         return
      endif
      twod_lats=ilat*1.d0

      status = nf90_inq_varid(ncid,'press',varid)
      if (status /= 0) then
         write(*,*)'get_const_2d 6:',nf90_strerror(status)
         return
      endif
      status = nf90_get_var(ncid, varid, twod_press)
      if (status /= 0) then
         write(*,*)'get_const_2d 7:',nf90_strerror(status)
         return
      endif
      twod_lnp=log(twod_press)
      twod_data_read=.true.
      status = nf90_close(ncid)
      if (status /= 0) then
         write(*,*)'get_const_2d 8:',nf90_strerror(status)
         return
      endif
    endif


    ! choose the correct species
    ispec=0
    do i=1,nspecs
      if (trim(twodspecs(i)) == trim(spec0))ispec=i
    enddo

    lat0 =  twod_lats(1)
    dlat =  twod_lats(2)-lat0

    ! loop over air parcels
    do jl=1,ntraj

      ! latitude interpolation
      rl = (lats(jl) -lat0)/dlat + 1
      i1 = max(int(rl),1)
      i2 = min(i1+1,twod_nlats)

      dl = rl -i1
      dl = max(dl,0.d0)
      dl = min(dl,1.d0)
      ! result = data(i1) * (1-dl) + data(i2) * dl
    
      ! pressure interpolation
      pr=press(jl)/100.
      pr=max(twod_press(twod_npress),pr)
      pr=min(twod_press(1),pr)
      k = maxloc(twod_press, (twod_press <= pr))
      k2=k(1)
      k2 = max(k2,2)  ! important if press > lowest lev (980 hPa)
      k1 = max(k2-1,1)

      dpr = (log(pr) -twod_lnp(k1))/(twod_lnp(k2)-twod_lnp(k1))
      !write(*,'(f8.3,2i3,f7.4,2f8.3)') pr,k1,k2,dpr,twod_press(k1:k2)
      dpr = max(dpr,0.d0)
      dpr = min(dpr,1.d0)
      x1 = twod_avg_all(i1,k1,m1,ispec)*(1-dpr) + twod_avg_all(i1,k2,m1,ispec)*dpr
      x2 = twod_avg_all(i2,k1,m1,ispec)*(1-dpr) + twod_avg_all(i2,k2,m1,ispec)*dpr
      y1 = x1*(1-dl) + x2 * dl

      x1 = twod_avg_all(i1,k1,m2,ispec)*(1-dpr) + twod_avg_all(i1,k2,m2,ispec)*dpr
      x2 = twod_avg_all(i2,k1,m2,ispec)*(1-dpr) + twod_avg_all(i2,k2,m2,ispec)*dpr
      y2 = x1*(1-dl) + x2 * dl

      const_arr(jl)=y1*(1-dt) + y2*dt
    enddo ! jl loop
  
  case default
    write (*,*) 'currently only OH,Cl,O(1D) and HO2 are implemented in get_const_2d'
    status = -121
    return
  end select


END SUBROUTINE get_const_2d



end Module messy_clamschem_get_const_2d

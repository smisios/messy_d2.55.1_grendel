!******************************************************************************
! File messy_clams_tools_interpolreg.f90
!
! This file is part of the Chemical Lagrangian Model of the
! Stratosphere (CLaMS). CLaMS is a hierarchy of numerical models and
! preprocessors which simulate transport, chemical reactions and
! mixing processes in the stratosphere.  
!
! Copyright (c) 2012
! N. Thomas, Forschungszentrum Juelich GmbH
! Last Modified By: N.Thomas
! Last Modified On: Fri Oct 11 11:41:40 2019
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
!******************************************************************************
!
! Module  messy_clams_tools_interpolreg
!
! This module contains some interpolation tools:
!
! subroutine get_index_height (levelgrid,nlevs,lev,ilev,xlower,loglev,extrapol)
! subroutine get_index_height (level, lev, nlevs, ilon, ilat, &
!                              ilev, xlower, loglev, extrapol)
!
! function interpolate_height (val1,val2,f,lev,logval)
!
! function interpolate_spatial (data, levelgrid, longrid, latgrid, &
!                               nlon, nlat, nlev, lon, lat, lev, &
!                               asc_lat, loglev, logval, extrapol)
! function interpolate_spatial_modlev (data, levelarr, longrid, latgrid, &
!                                      nlon, nlat, nlev, lon, lat, lev, &
!                                      asc_lat, loglev, logval, extrapol)
!
! function interpolate_time (val1, val2, time_interval, timestep)
! SUBROUTINE interpolate_param (lons, lats, levs, calc, itime, ipasttime, &
!                               iparam, paramname, paramvals)          
!
!
!******************************************************************************

Module messy_clams_tools_interpolreg

  INTERFACE get_index_height
     MODULE PROCEDURE get_index_height_lev1d
     MODULE PROCEDURE get_index_height_lev3d
  END INTERFACE get_index_height

CONTAINS

! op_pj_20160829+: workaround for Lahey/Fujitsu Compiler 8.1b and NAG
#if defined (LF) || defined (NAGFOR) || defined (__PGI)
LOGICAL FUNCTION ISNAN(x)
  use messy_clams_global, only: prec
  IMPLICIT NONE
  REAL(prec), INTENT(in) :: x
  ISNAN = (x/=x)
END FUNCTION ISNAN
#endif
! op_pj_20160829-
!***********************************************************************
!
!***********************************************************************
subroutine get_index_height_lev1d (levelgrid,nlevs,lev,ilev,xlower,&
                                   ascending,loglev,extrapol)

  use messy_clams_global, only: prec

  IMPLICIT NONE

  real(prec), dimension(:), pointer :: levelgrid
  integer,    intent(in)   :: nlevs
  real(prec), intent(in)   :: lev
  integer,    intent(out)  :: ilev
  real(prec), intent(out)  :: xlower
  logical,    intent(in)   :: ascending, loglev, extrapol

  logical :: found

  found = .false.
  ilev = 1

  if (ascending) then  ! levels in ascending order

     if (lev < levelgrid(1)) then
        if (extrapol) then
           ilev   = 1
           xlower = 1.
        else
           ilev = -1
        endif
     elseif (lev >= levelgrid(nlevs)) then
        if (extrapol) then
           ilev   = nlevs-1
           xlower = 0.
        else
           ilev = -1
        endif
     else
        do while (ilev<nlevs .and. .not. found)
           if (levelgrid(ilev)<=lev .and. lev<=levelgrid(ilev+1)) then
              found = .true.
              if (.not. loglev) then  
                 xlower = (levelgrid(ilev+1)-lev) / &
                          (levelgrid(ilev+1)-levelgrid(ilev))
              else
                 xlower = (log(levelgrid(ilev+1))-log(lev)) / &
                          (log(levelgrid(ilev+1))-log(levelgrid(ilev)))
              endif
           else
              ilev = ilev+1
           endif
        enddo
        if (.not. found) ilev=-1
     endif
     
  else ! levels in descending order
     
     if (lev > levelgrid(1)) then
        if (extrapol) then
           ilev   = 1
           xlower = 1.
        else
           ilev = -1
        endif
     elseif (lev <= levelgrid(nlevs)) then
        if (extrapol) then
           ilev   = nlevs-1
           xlower = 0.
        else
           ilev = -1
        endif
     else
        do while (ilev<nlevs .and. .not. found)
           if (levelgrid(ilev)>=lev .and. lev>=levelgrid(ilev+1)) then
              found = .true.
              if (.not. loglev) then 
                 xlower = (levelgrid(ilev+1)-lev) / &
                          (levelgrid(ilev+1)-levelgrid(ilev))
              else
                 xlower = (log(levelgrid(ilev+1))-log(lev)) / &
                          (log(levelgrid(ilev+1))-log(levelgrid(ilev)))
              endif
           else
              ilev = ilev+1
           endif
        enddo
        if (.not. found) ilev=-1
     endif
     
  endif

end subroutine get_index_height_lev1d

!***********************************************************************
!
!***********************************************************************
subroutine get_index_height_lev3d (level, lev, nlevs, ilon, ilat, ilev, &
                                   xlower, ascending, loglev, extrapol)

  use messy_clams_global, only: prec

  IMPLICIT NONE

  real(prec), intent(in)  :: level(:,:,:)
  real(prec), intent(in)  :: lev
  integer,    intent(in)  :: nlevs, ilon, ilat
  integer,    intent(out) :: ilev
  real(prec), intent(out) :: xlower
  logical,    intent(in)  :: ascending, loglev, extrapol

  logical :: found

  found = .false.
  ilev = 1


! !!!!! ACHTUNG: in traj wird "rz" zurueckgegeben in pos_dyn "xlower"
!  => rz = 1. - xlower !!!

  if (ascending) then  ! levels in ascending order
     if (lev<level(ilon,ilat,1)) then
        if (extrapol) then
           ilev   = 1
           xlower = 1.
        else
           ilev = -1
        endif
     elseif (lev>level(ilon,ilat,nlevs)) then
        if (extrapol) then
           ilev   = nlevs-1
           xlower = 0.
        else
           ilev = -1
        endif
     else
        do while (ilev<nlevs .and. .not. found)
           if (level(ilon,ilat,ilev)<=lev .and. lev<=level(ilon,ilat,ilev+1)) then
              found = .true.
              if (.not. loglev) then
                 xlower = (level(ilon,ilat,ilev+1)-lev) / &
                          (level(ilon,ilat,ilev+1)-level(ilon,ilat,ilev))
              else
                 xlower = (log(level(ilon,ilat,ilev+1))-log(lev)) / &
                          (log(level(ilon,ilat,ilev+1))-log(level(ilon,ilat,ilev)))
             endif
           else
              ilev = ilev+1
           endif
        end do
        if (.not. found) ilev=-1
     endif
        
  else  ! levels in descending order

     if (lev>level(ilon,ilat,1)) then
        if (extrapol) then
           ilev   = 1
           xlower = 1.
        else
           ilev = -1
        endif
     elseif(lev<level(ilon,ilat,nlevs)) then
        if (extrapol) then
           ilev   = nlevs-1
           xlower = 0.
        else
           ilev = -1
        endif
     else
        do while (ilev<=nlevs .and. .not. found)
           if (level(ilon,ilat,ilev)>=lev .and. lev>=level(ilon,ilat,ilev+1)) then
              found = .true.
              if (.not. loglev) then
                 xlower = (level(ilon,ilat,ilev+1)-lev) / &
                          (level(ilon,ilat,ilev+1)-level(ilon,ilat,ilev))
              else
                 xlower = (log(level(ilon,ilat,ilev+1))-log(lev)) / &
                          (log(level(ilon,ilat,ilev+1))-log(level(ilon,ilat,ilev)))
              endif
           else
              ilev = ilev+1
           endif
        end do
        if (.not. found) ilev=-1
     endif
     
  endif

 
end subroutine get_index_height_lev3d

!***********************************************************************
!
!***********************************************************************
function interpolate_height (val1,val2,f,lev,logval)

  use messy_clams_global, only: prec, mdi, eps

  implicit none

  real(prec)        :: interpolate_height 
  real(prec)        :: val1, val2, f, lev
  integer, optional :: logval
  
  ! local variables
  real(prec)    :: value, value_lin, value_log
  integer       :: log_val


  log_val = 0  ! lin. interpolation (default)
  if (present(logval)) log_val = logval
  
  if (abs((val1-mdi)/mdi)<=eps .or. abs((val2-mdi)/mdi)<=eps) then
     interpolate_height = mdi
     return
  endif

  if (log_val==0) then  ! linear interpolation
     value = f * val1 + (1-f) * val2       
  elseif (log_val==1) then   ! log. lin. interpolation
     if (val1<=0. .or. val2<=0) then
        write (*,*) 'in sub. interpolate_height: missing_value gesetzt !!!'
        value = mdi
     else
        value = val1**f * val2**(1-f)
     endif
  elseif (log_val==2) then   ! lin. + log.lin. interpolation 
     
     if (lev<=500.) then       ! lin. interpolation
        value = f * val1 + (1-f) * val2
        
     elseif (lev>=1000.) then  ! log. lin. interpolation
        if (val1<=0. .or. val2<=0.) then
           write (*,*) 'sub. interpolate_height: missing_value gesetzt !!!'
           value = mdi
        else
           value = val1**f * val2**(1-f)
        endif
        
     else  ! 500 < level < 1000
        if (val1<=0. .or. val2<=0.) then
           write (*,*) 'sub. interpolate_height: missing_value gesetzt !!!'
           value = mdi
        else  
           value_lin = f * val1 + (1-f) * val2
           value_log = val1**f * val2**(1-f)
           value = f * value_lin + (1-f) * value_log
        endif
        
     endif

  end if

  interpolate_height = value

end function interpolate_height

!***********************************************************************
!
!***********************************************************************
function interpolate_spatial (data, levelgrid, longrid, latgrid, &
                              nlon, nlat, nlev, lon, lat, lev, &
                              asc_lat, asclev, loglev, logval, extrapol)

  USE messy_clams_global,       ONLY: prec, mdi, eps
  USE messy_clams_tools_bval3d, ONLY: bval3d
                                                              
  IMPLICIT NONE

  real(prec)                             :: interpolate_spatial
  REAL(PREC), DIMENSION(:),     POINTER  :: levelgrid
  REAL(PREC), DIMENSION(:),     POINTER  :: longrid, latgrid
  REAL(PREC), DIMENSION(:,:,:), POINTER  :: data
  real(prec)           :: lon, lev, lat
  integer              :: nlon, nlat, nlev
  logical              :: asc_lat
  integer,    optional :: logval
  logical,    optional :: asclev, loglev, extrapol

  integer :: log_val = 0         ! interpolate parameter values lin. log.
  logical :: asc_lev = .true.     ! levels in ascending order?
  logical :: log_lev = .false.   ! true: interpolate log. lin
                                 ! false: interpolate linear
  logical :: extrapolate = .true.
  

  real(prec) :: val1, val2
  real(prec) :: xlower
  integer    :: ilev

  if (present(asclev)) asc_lev = asclev
  if (present(logval)) log_val = logval
  if (present(loglev)) log_lev = loglev
  if (present(extrapol)) extrapolate = extrapol

   if (abs((lev-mdi)/mdi) <= eps) then
     interpolate_spatial = mdi
     return
  endif

  call get_index_height (levelgrid,nlev,lev,ilev,xlower,asc_lev,log_lev,extrapolate)

  if (ilev==-1) then
     interpolate_spatial = mdi
     
  else
     
     ! interpolate horizontal (at the given latitude and longitude) 
     ! for lower level
     val1 = bval3d (data,ilev,lat,lon, &
                    longrid,latgrid,nlon,nlat,asc_lat)   

     ! interpolate horizontal (at the given latitude and longitude)
     ! for upper level
     val2 = bval3d (data,ilev+1,lat,lon, &
                    longrid,latgrid,nlon,nlat,asc_lat)   

     ! interpolate between levels
     if (abs((val1-mdi)/mdi)<=eps .or. abs((val2-mdi)/mdi)<=eps) then
        interpolate_spatial = mdi
     else
        interpolate_spatial = interpolate_height (val1,val2,xlower,lev,logval=log_val)
     endif
     
  endif

end function interpolate_spatial

!***********************************************************************
!
!***********************************************************************
function interpolate_spatial_modlev (data, levelarr, longrid, latgrid, &
                                     nlon, nlat, nlev, lon, lat, lev, &
                                     asc_lat, asclev, loglev, logval, extrapol)

  USE messy_clams_global,       ONLY: prec, mdi, eps
  USE messy_clams_tools_bval3d, ONLY: bval3d_new
                                                    
  IMPLICIT NONE

  real(prec)                            :: interpolate_spatial_modlev
  REAL(PREC), DIMENSION(:),     POINTER :: longrid, latgrid
  REAL(PREC), DIMENSION(:,:,:), POINTER :: levelarr, data
  real(prec)           :: lon, lev, lat
  integer              :: nlon, nlat, nlev
  logical              :: asc_lat
  integer,    optional :: logval
  logical,    optional :: asclev, loglev, extrapol

  integer :: log_val = 0         ! interpolate parameter values lin. log.
  logical :: asc_lev = .true.    ! levels in ascending order
  logical :: log_lev = .false.   ! true: interpolate log. linear
                                 ! false: interpolate linear
  logical :: extrapolate = .true.
  

  real(prec) :: val1, val2, val3, val4
  real(prec) :: xlower
  integer    :: ilev, ilat1, ilat2, ilon1, ilon2, i

  if (present(asclev)) asc_lev = asclev
  if (present(logval)) log_val = logval
  if (present(loglev)) log_lev = loglev
  if (present(extrapol)) extrapolate = extrapol

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   if (abs((lev-mdi)/mdi) <= eps) then
     interpolate_spatial_modlev = mdi
     return
  endif

  ! ilon1 = mod(int((gx0+gdx)/gdx),nlon)+1
  ! ilon2 = mod(int((lon-(gx0+gdx))/gdx+1),nlon)+1
  ! ilat1 = (lat-(gy0+gdy))/gdy+1
  ! ilat2 = (lat-(gy0+gdy))/gdy+2

!!!!! unregelmaessiges Gitter moeglich:
! op_pj_20160606+
!!$  if (lon>360.) lon = mod(lon,360.) !???
  if (lon>360._prec) lon = mod(lon,360._prec) !???
! op_pj_20160606-
  if (lon<longrid(1) .or. lon>longrid(nlon)) then
     ilon1 = nlon
     ilon2 = 1
  else
     do i = 1, nlon-1
        if (longrid(i)<=lon .and. lon<=longrid(i+1)) then 
           ilon1 = i
           ilon2 = i+1
        endif
     enddo
  endif
  
  if (asc_lat) then  ! latitudes in ascending order
     if (lat<latgrid(1)) then
        ilat1 = 1
        ilat2 = 1
     elseif (lat>latgrid(nlat)) then
        ilat1 = nlat
        ilat2 = nlat
     else        
        do i = 1, nlat-1
           if (latgrid(i)<=lat .and. lat<=latgrid(i+1)) then
              ilat1 = i
              ilat2 = i+1
           endif
        enddo
     endif
     
  else ! latitudes in descending order
     if (lat>latgrid(1)) then
        ilat1 = 1
        ilat2 = 1
     elseif  (lat<latgrid(nlat)) then
        ilat1 = nlat
        ilat2 = nlat
     else        
        do i = 1, nlat-1 
           if (latgrid(i)>=lat .and. lat>=latgrid(i+1)) then
              ilat1 = i
              ilat2 = i+1
           endif
        enddo
     endif
  endif
  
  call get_index_height (levelarr, lev, nlev, ilon1, ilat1, ilev, &
                         xlower, asc_lev, log_lev, extrapolate)
  if (ilev==-1) then
     val1=mdi
  else
     val1 = interpolate_height (data(ilon1,ilat1,ilev),data(ilon1,ilat1,ilev+1), &
                                xlower,lev,logval=log_val)
  endif
  
  call get_index_height (levelarr, lev, nlev, ilon2, ilat1, ilev, &
                         xlower, asc_lev, log_lev, extrapolate)
  if (ilev==-1) then
     val2 = mdi
  else
     val2 = interpolate_height (data(ilon2,ilat1,ilev),data(ilon2,ilat1,ilev+1), &
          xlower,lev,logval=log_val)
  endif
  
  if (ilat1==ilat2) then  ! near north or south pole

     ilon1 = mod(ilon1+nlon/2,nlon)
     if (ilon1 == 0) ilon1 = 1
     call get_index_height (levelarr, lev, nlev, ilon1, ilat2, ilev, &
                            xlower, asc_lev, log_lev, extrapolate)
     if (ilev==-1) then
        val3 = mdi
     else
        val3 = interpolate_height (data(ilon1,ilat2,ilev),data(ilon1,ilat2,ilev+1), &
                                   xlower,lev,logval=log_val)
     endif
       
     ilon2 = mod(ilon2+nlon/2,nlon)
     if (ilon2 == 0) ilon2 = 1
     call get_index_height (levelarr, lev, nlev, ilon2, ilat2, ilev, &
                            xlower, asc_lev, log_lev, extrapolate)
     if (ilev==-1) then
        val4 = mdi
     else
        val4 = interpolate_height (data(ilon2,ilat2,ilev),data(ilon2,ilat2,ilev+1), &
                                   xlower,lev,logval=log_val)
     endif
     
  else
     
     call get_index_height (levelarr, lev, nlev, ilon2, ilat2, ilev, &
                            xlower, asc_lev, log_lev, extrapolate)
     if (ilev==-1) then
        val3 = mdi
     else
        val3 = interpolate_height (data(ilon2,ilat2,ilev),data(ilon2,ilat2,ilev+1), &
                                   xlower,lev,logval=log_val)
     endif
     
     call get_index_height (levelarr, lev, nlev, ilon1, ilat2, ilev, &
                            xlower, asc_lev, log_lev, extrapolate)
     if (ilev==-1) then
        val4 = mdi
     else
        val4 = interpolate_height (data(ilon1,ilat2,ilev),data(ilon1,ilat2,ilev+1), &
                                   xlower,lev,logval=log_val)
     endif
     
  endif
  
  ! compute the value at the given latitude and longitude           
  interpolate_spatial_modlev = bval3d_new (val1,val2,val3,val4,lat,lon, &
                                    longrid,latgrid, nlon, nlat,asc_lat)   


end function interpolate_spatial_modlev

!***********************************************************************
!
!***********************************************************************
function interpolate_time (val1, val2, time_interval, timestep)

  USE messy_clams_global,  ONLY: prec, dp, mdi, eps

  implicit none

  REAL(PREC)      :: interpolate_time
  REAL(PREC)      :: val1, val2
  REAL(dp)        :: time_interval, timestep

  if (abs((val1-mdi)/mdi)<=eps .or. abs((val2-mdi)/mdi)<=eps) then
     interpolate_time = mdi
  else
     interpolate_time = ((time_interval-timestep)/time_interval)*val1 + (timestep/time_interval)*val2 
  endif

end function interpolate_time

!*********************************************************************
! interpolate parameter values at given positions 
!*********************************************************************
SUBROUTINE interpolate_param (lons, lats, levs, itime, ipasttime, &
                             iparam, paramname, paramvals, lcoupled, calc)          

   USE messy_clams_global,            ONLY: rank, prec, dp, dnparts_max, dnparts, &
                                            mdi, eps, timetype, dates30, &
                                            PREDATA, FUTDATA, PARAM, irdday, irdsec, &
                                            latgrid, longrid, levelgrid, leveldt, &
                                            nx, ny, nz, &
                                            nparams, paramnames, &
                                            thetagrid, ntheta, &
                                            asc_lat, asc_level, loglev, logpress, &
                                            level_is_vertcoor
   USE messy_clams_tools_dateconv,    ONLY: ymds2js_interface
   USE messy_clams_tools_utils,       ONLY: str_pos

   IMPLICIT NONE

   REAL(PREC)        :: lons(dnparts_max), lats(dnparts_max), levs(dnparts_max)
   LOGICAL, optional :: calc(dnparts_max)
   TYPE(timetype)    :: itime, ipasttime
   INTEGER           :: iparam
   CHARACTER(*)      :: paramname
   LOGICAL           :: lcoupled
 
   REAL(PREC) :: paramvals(dnparts_max)

   ! local variables
   LOGICAL        :: calculate(dnparts_max)
   REAL(DP)       :: ijultime, lastjultime
   REAL(DP)       :: time_interval, timestep
   REAL(PREC)     :: value1, value2
   INTEGER        :: ipart, logval, ipos

 
   calculate = .true.
   if (present(calc)) calculate = calc

   time_interval = irdday*24*3600 +irdsec

   if (.not. lcoupled) then
      lastjultime = ymds2js_interface (ipasttime%year, ipasttime%month, &
                          ipasttime%day, ipasttime%sec, dates30)
      ijultime = ymds2js_interface (itime%year, itime%month, itime%day, itime%sec, dates30)
      timestep = ijultime - lastjultime 
   else
      timestep = 0 ! no interpolation in time
   endif

   if (.not. lcoupled) then
      if (paramname=='EQLAT' .or. paramname=='PV') then
         ipos = str_pos (nparams, paramnames, 'THETA')
      endif
   endif

!!!!!
   ! if (rank==0) then
   !    write (*,*) 'asc_lat, asc_level, loglev, logval: ',asc_lat, asc_level, loglev, logval
   !    write (*,*) 'level_is_vertcoor: ',level_is_vertcoor
   ! endif
   
   DO ipart = 1, dnparts

      IF (calculate(ipart)) THEN

         logval = 0
         if (paramname=='PRESS' .and. logpress) logval = 2

         if (paramname=='EQLAT' .or. paramname=='PV') then

            if (.not. lcoupled) then
               value1 = interpolate_spatial (predata(iparam)%values, &
                              thetagrid, longrid, latgrid, &
                              nx, ny, ntheta, lons(ipart), lats(ipart), param(ipos)%values(ipart), &
                              asc_lat, asclev=asc_level, loglev=loglev, logval=logval, &
                              extrapol=.false.)
               value2 = interpolate_spatial (futdata(iparam)%values, &
                              thetagrid, longrid, latgrid, &
                              nx, ny, ntheta, lons(ipart), lats(ipart), param(ipos)%values(ipart), &
                              asc_lat, asclev=asc_level, loglev=loglev, logval=logval, &
                              extrapol=.false.)
            elseif (paramname=='PV') then
               !!!!! => used in CIRRUS !
               value1 = 1.   
            else
               value1 = mdi
            endif

         else
            if (.not. level_is_vertcoor) then
               value1 = interpolate_spatial_modlev (predata(iparam)%values, &
                          leveldt, longrid, latgrid, &
                          nx, ny, nz, lons(ipart), lats(ipart), levs(ipart), &
                          asc_lat, asclev=asc_level, loglev=loglev,  logval=logval, &
                          extrapol=.false.)
               if (.not. lcoupled) then
                  value2 = interpolate_spatial_modlev (futdata(iparam)%values, &
                          leveldt, longrid, latgrid, &
                          nx, ny, nz, lons(ipart), lats(ipart), levs(ipart), &
                          asc_lat, asclev=asc_level, loglev=loglev,  logval=logval, &
                          extrapol=.false.)
               endif
            else
               value1 = interpolate_spatial (predata(iparam)%values, &
                              levelgrid, longrid, latgrid, &
                              nx, ny, nz, lons(ipart), lats(ipart), levs(ipart), &
                              asc_lat, asclev=asc_level, loglev=loglev, logval=logval, &
                              extrapol=.false.)
               if (.not. lcoupled) then
                  value2 = interpolate_spatial (futdata(iparam)%values, &
                              levelgrid, longrid, latgrid, &
                              nx, ny, nz, lons(ipart), lats(ipart), levs(ipart), &
                              asc_lat, asclev=asc_level, loglev=loglev, logval=logval, &
                              extrapol=.false.)
               endif
            endif

         endif

         if (.not. lcoupled) then
            paramvals(ipart) = interpolate_time (value1, value2, time_interval, timestep)
         else
            paramvals(ipart) = value1
         endif

         if (abs((paramvals(ipart)-mdi)/mdi)<=eps .OR. IsNAN(paramvals(ipart))) then
            IF (IsNAN(paramvals(ipart))) THEN
               write(100,*) 'ipart=',ipart
               write(100,*) 'paramname=',trim(paramname)
               write(100,*) 'value1, value2, paramvals=',value1, value2, paramvals(ipart)
               write(100,*) 'lon', lons(ipart), 'lat', lats(ipart), 'lev', levs(ipart)
               stop
            ENd IF
            if (trim(paramname)=='U' .or. trim(paramname)=='V' .or. &
                trim(paramname)=='W' .or. trim(paramname)=='TEMP' .or. trim(paramname)=='PRESS') then
               lons(ipart) = mdi
               lats(ipart) = mdi
               levs(ipart) = mdi
               paramvals(ipart) = mdi
            endif

            ! If there are missing values on EQLAT: replace with LAT
            if (paramname=="EQLAT") then
               paramvals(ipart) = lats(ipart)
            endif


         endif

      ELSE
         paramvals(ipart) = mdi
         lons(ipart) = mdi
         lats(ipart) = mdi
         levs(ipart) = mdi
      ENDIF

   ENDDO

!!$   if (paramname == 'EQLAT') then
!!$      if (rank==0) then
!!$         write (*,*) 'lats=',lats(1:100)
!!$         write (*,*) 'eqlats=',paramvals(1:100)
!!$         write (*,*) 'lons=',lons(1:100)
!!$         write (*,*) 'levs=',levs(1:100)
!!$         write (*,*) 'thetas=',param(ipos)%values(1:100)
!!$         write (*,*) 'thetagrid=',thetagrid
!!$      endif
!!$   endif


END SUBROUTINE interpolate_param


!*********************************************************************
! interpolate parameter values at given position
!*********************************************************************
SUBROUTINE interpolate_param_one_point (status, lon, lat, lev, itime, ipasttime, &
                                        paramname, paramvalue)          

   USE messy_clams_global,            ONLY: prec, dp, dnparts_max, dnparts, &
                                            mdi, eps, timetype,   &
                                            PREDATA, FUTDATA, met_freq, &
                                            latgrid, longrid, levelgrid, leveldt, &
                                            nx, ny, nz, nparams, &
                                            asc_lat, asc_level, loglev, logpress, &
                                            level_is_vertcoor
   USE messy_clams_tools_dateconv,    ONLY: ymds2js

   IMPLICIT NONE

   REAL(PREC)        :: lon, lat, lev
   REAL(PREC)        :: paramvalue
   REAL(DP)          :: itime, ipasttime
   CHARACTER(*)      :: paramname
   integer           :: status 

   ! local variables
   !REAL(DP)       :: ijultime, lastjultime
   REAL(DP)       :: time_interval, timestep
   REAL(PREC)     :: value1, value2
   INTEGER        :: ipart, logval
   INTEGER        :: iparam, i

   status = 0

   time_interval = met_freq*3600  ! in seconds

   timestep = itime - ipasttime 

   logval = 0
   if (paramname=='PRESS' .and. logpress) logval = 2

   iparam = -1
   do i = 1, nparams
      if (predata(i)%name == paramname) iparam = i
   enddo
   if (iparam==-1) then
      status = -1
      return
   endif
      
   if (.not. level_is_vertcoor) then
      value1 = interpolate_spatial_modlev (predata(iparam)%values, &
                          leveldt, longrid, latgrid, &
                          nx, ny, nz, lon, lat, lev, &
                          asc_lat, asclev=asc_level, loglev=loglev,  logval=logval, &
                          extrapol=.false.)
      value2 = interpolate_spatial_modlev (futdata(iparam)%values, &
                          leveldt, longrid, latgrid, &
                          nx, ny, nz, lon, lat, lev, &
                          asc_lat, asclev=asc_level, loglev=loglev,  logval=logval, &
                          extrapol=.false.)
   else
      value1 = interpolate_spatial (predata(iparam)%values, &
                              levelgrid, longrid, latgrid, &
                              nx, ny, nz, lon, lat, lev, &
                              asc_lat, asclev=asc_level, loglev=loglev, logval=logval, &
                              extrapol=.false.)
      value2 = interpolate_spatial (futdata(iparam)%values, &
                              levelgrid, longrid, latgrid, &
                              nx, ny, nz, lon, lat, lev, &
                              asc_lat, asclev=asc_level, loglev=loglev, logval=logval, &
                              extrapol=.false.)
   endif

   paramvalue = interpolate_time (value1, value2, time_interval, timestep)


END SUBROUTINE interpolate_param_one_point

End Module messy_clams_tools_interpolreg

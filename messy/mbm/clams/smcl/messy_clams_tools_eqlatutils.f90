!******************************************************************************
!
! This file is part of the Chemical Lagrangian Model of the
! Stratosphere (CLaMS). CLaMS is a hierarchy of numerical models and
! preprocessors which simulate transport, chemical reactions and
! mixing processes in the stratosphere.  
!
! Copyright (c) 2006
! J.-U. Grooss, P. Konopka Forschungszentrum Juelich GmbH
! Last Modified By: N.Thomas
! Last Modified On: Thu Jun  8 11:20:34 2017
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
! Module that calculates a field of equivalent latitude based on PV values 
! on isentropic levels. 
!
! SUBROUTINE get_pvlim (pv_field,lat_grid,nlats,nlons,lats,limits)
! SUBROUTINE get_eqlat_int (filename,lats,lons,theta,eqlat,npoints,pv)
!
!******************************************************************************
Module messy_clams_tools_eqlatutils

CONTAINS

!*************************************************************************************
! This subroutine calculates pv values for equivalent latitude bins
! given by lat_grid. 
! Developed from DMcK's routines  "areas_utils"
! J.U. Grooss, October 1997
!*************************************************************************************
SUBROUTINE get_pvlim (status,pv_field,lat_grid,nlats,nlons,lats,limits)

  ! ! use nrtype
  ! ! use nrutil
  ! use nr, only: indexx_sp
  use messy_clams_tools_utils, ONLY: quick_sort_index

  use messy_clams_global, ONLY: prec, dp, mdi, eps
  
  IMPLICIT NONE

  INTEGER,    intent(inout)                       :: status
  INTEGER,    intent(in)                          :: nlats, nlons
  REAL(prec), intent(in),  dimension(nlons,nlats) :: pv_field
  REAL(prec), intent(in),  dimension(nlats)       :: lat_grid
  REAL(prec), intent(out), dimension(nlats)       :: lats,limits

  REAL(kind(0.d0)), dimension(nlats) :: areawt
  REAL(kind(0.d0)), dimension(nlats) :: areas, area_north, area_south
  REAL(DP)   :: torad, awt
  REAL(prec) :: minlat, maxlat, lat1, lat2, dlat
  REAL(prec) :: pvmin, pvmax, pv_val, epsilon
  REAL(prec), dimension(nlats*nlons) :: pv_current

  INTEGER :: npoints, ilat, ilon, cnt, i, in, is
  INTEGER, dimension(nlats*nlons) :: ind_lat, ind_pv

  LOGICAL :: shemi, nhemi, global

  !write (*,*) 'begin get_pvlim'

  status = 0

  npoints = nlats*nlons
  torad   = 4.d0*datan(1.d0)/180.d0

  minlat  = minval(lat_grid)
  maxlat  = maxval(lat_grid)

  ! check wether hemispheric or global data are present
  if (minlat < -75.) minlat=-90.
  if (maxlat >  75.) maxlat=90.
  if (minlat > -15.) minlat=0.
  if (maxlat <  15.) maxlat=0.
  global= (minlat == -90.) .and. (maxlat == 90.)
  nhemi = (minlat == 0.)   .and. (maxlat == 90.)
  shemi = (minlat == -90.) .and. (maxlat == 0.)
  
  if (.not.(global .or. shemi .or. nhemi)) then
     write (*,*) 'get_pvlim: unclear latitude data range' 
     status = -1
     return
  endif

  ! calculate area weight of grid points:
  DO ilat = 1,nlats
     IF (ilat /= 1) THEN
        lat1 = (lat_grid(ilat-1)+lat_grid(ilat))/2.
     ELSE
        dlat = lat_grid(1)-lat_grid(2)
        lat1 = lat_grid(1)+dlat/2.
        IF (lat1+dlat/2. > maxlat) lat1 = maxlat
        IF (lat1+dlat/2. < minlat) lat1 = minlat
     ENDIF
     
     IF (ilat /= nlats) THEN
        lat2 = (lat_grid(ilat)+lat_grid(ilat+1))/2.
     ELSE
        dlat = lat_grid(nlats-1)-lat_grid(nlats)
        lat2 = lat_grid(nlats)-dlat/2.
        IF (lat2-dlat/2. > maxlat) lat2 = maxlat
        IF (lat2-dlat/2. < minlat) lat2 = minlat
     ENDIF
     areawt(ilat) = abs(sin(lat1*torad)-sin(lat2*torad))/(2.d0 *nlons)
  ENDDO

  ! set latitude grid for PV / latitude relation
  lats        = lat_grid
  lats(1)     = maxlat
  lats(nlats) = minlat
  areas       = (1.d0-sin(lats*torad))/2.d0

  ! To account for data points with missing PV values, reduce the array
  ! area (area north of lat) for the corresponding latitude 
  !
  ! get valid area 
  do ilat = 1,nlats
     cnt = 0
     do ilon = 1,nlons
        if (abs((pv_field(ilon,ilat)-mdi)/mdi) <= eps) cnt = cnt+1
     enddo
     areas(ilat:nlats) = areas(ilat:nlats) - cnt * areawt(ilat)
  enddo

  pvmin = minval(pv_field,mask=(abs((pv_field-mdi)/mdi) > eps ))
  pvmax = maxval(pv_field,mask=(abs((pv_field-mdi)/mdi) > eps ))
  if (pvmin == HUGE(pv_field)) then
     write (*,*) 'Fehler in nc_eqlat_utils.f90, sub. get_pvlim:'
     write (*,*) 'PV enthaelt nur missing-values !!!'
     write (*,*)
     status = -1
     return
  endif

  ! bring PV data into one-dimensional array
  pv_current = reshape(pv_field(:,:),(/npoints/))
  ind_lat = reshape(spread((/(i,i=1,nlats)/),dim=1,ncopies=nlons),(/npoints/))

  !sort pv values
  !CALL bubble_sort_ind(pv_current,ind_pv,npoints)
  !call indexx_sp (-pv_current, ind_pv)
  call quick_sort_index (-pv_current, ind_pv)

  if (global) then 
     ! global data present
     area_north(:) = 0.0
     limits(:) = mdi
     in = 2
     do i = 1,npoints
        pv_val = pv_current(ind_pv(i))
        awt  = areawt(ind_lat(ind_pv(i)))
        if((abs((pv_val-mdi)/mdi) > eps)) then
           area_north(in) = area_north(in) + awt
           if (area_north(in) > areas(in)) then
              limits(in) = pv_val
              if (in == nlats) exit
              area_north(in+1) = area_north(in)
              in = in+1
           endif
        endif
     enddo
  endif
  
  if (nhemi) then 
     ! North hemispheric data 
     area_north(:) = 0.0
     limits(:) = mdi
     in = 2
     do i=1,npoints
        pv_val = pv_current(ind_pv(i))
        awt  = areawt(ind_lat(ind_pv(i)))
        if (abs((pv_val-mdi)/mdi) > eps) then
           area_north(in) = area_north(in)+awt
           if (area_north(in) > areas(in)) then
              limits(in) = pv_val
              if (in == nlats) exit
              area_north(in+1) = area_north(in)
              in = in+1
           endif
        endif
     enddo
  endif
  
  if (shemi) then 
     area_south(:) = 1.0
     limits(:) = mdi
     is = nlats -1
     do i=npoints,1,-1
        pv_val = pv_current(ind_pv(i))
        awt  = areawt(ind_lat(ind_pv(i)))
        if((abs((pv_val-mdi)/mdi) > eps)) then
           area_south(is) = area_south(is)-awt
           if (area_south(is).lt.areas(is)) then
              limits(is) = pv_val
              if (is == 1) exit
              area_south(is-1) = area_south(is)
              is = is-1
           endif
        endif
     enddo
  endif
  
  epsilon       = (pvmax-pvmin)*0.0001
  limits(1)     = pvmax+epsilon
  limits(nlats) = pvmin-epsilon
  
  if (any (ABS((limits-mdi)/mdi)<=eps)) then
     write (*,*) 'Unvorgesehener Fall in get_pvlim: Winddaten ueberpruefen'
     status = -1
     return
  endif

  !write (*,*) 'end get_pvlim'
 
END subroutine get_pvlim



!*************************************************************************************
! This subroutine calculates the equivalent latitude field from PV of a given UKMO data 
! set stored in <filename>. Equivalent latitudes are calculated for <npoints> number of
! points given by latitudes <lats> and longitude <lons>. PV values are interpolated from
! the UKMO data points to given location. The result is stored into <eqlat>.
! J.U. Grooss, 25.11.97
!*************************************************************************************
SUBROUTINE get_eqlat_int (status,filename,lats,lons,levs,eqlat,startindex,endindex)

  USE messy_clams_global,        ONLY: PREC, mdi, eps,rank
  USE messy_clams_tools_ncutils, ONLY: nc_grid_descr, nc_get_level, nc_get_values3_cf

  IMPLICIT NONE


  REAL(PREC), DIMENSION(:), POINTER :: lats, lons, levs, eqlat
  INTEGER                           :: status, startindex, endindex
  CHARACTER(*)                      :: filename

  REAL(PREC), dimension(:), pointer :: lat_grid, lon_grid, levels
  REAL(prec), allocatable           :: lat_limits(:), llats(:), llons(:)
  
  REAL(prec), dimension(:),     allocatable :: pv_int
  REAL(PREC), dimension(:,:),   allocatable :: pv_limits
  REAL(prec), dimension(:,:,:), allocatable :: pv_field
  INTEGER,    dimension(:),     allocatable :: lev_index
  INTEGER                                   :: ncid, dimid
  INTEGER                                   :: npoints, nlats, nlons, nlevs, n, na  
  INTEGER                                   :: i, j, il(1), ilev, ilat, ilon, ilonp1
  REAL(PREC)                                :: dlat, dlon, dpv, x1, x2

  !write (*,*) rank, 'begin get_eqlat_int'

  status = 0

  npoints = endindex - startindex + 1
  !write (*,*) 'npoints=',npoints

  CALL nc_grid_descr(filename,nlons,nlats,lon_grid,lat_grid,status)
  if (status/=0) then
     write (*,*) 'Error in subroutine nc_grid_descr !'
     return
  endif
  !write (*,*) 'nlats, nlons=',nlats,nlons

  ALLOCATE (lat_limits(nlats))
  ALLOCATE (llats(npoints), llons(npoints), pv_int(npoints))

  llats = lats(startindex:endindex)
  llons = lons(startindex:endindex)
  
  call nc_get_level(filename,nlevs,levels,err=status)
  if (status/=0) then
     write (*,*) 'Error in subroutine nc_get_level !'
     return
  endif
  !write (*,*) 'Levels in isentropic dataset:',levels

  ALLOCATE (pv_field(nlons,nlats,nlevs))

  !CALL nc_get_values3(status,filename,'PV',pv_field)
  CALL nc_get_values3_cf(status,filename,'PV',pv_field)
  if (status/=0) then
     write (*,*) 'Error in subroutine nc_get_values3_cf !'
     return
  endif
 
  ALLOCATE (lev_index(npoints))
  lev_index = -1
  do i = 1, npoints
     il = minloc(abs(levs(i)-levels))
     lev_index(i) = il(1)
  enddo
  
  ALLOCATE (pv_limits(nlats,nlevs))
  pv_limits = mdi
  !write (*,*) rank, 'min und max level index', minval(lev_index), maxval(lev_index)
  do ilev = minval(lev_index), maxval(lev_index)
     CALL get_pvlim (status, pv_field(:,:,ilev),lat_grid,nlats,nlons,lat_limits,pv_limits(:,ilev))
     if (status/=0) then
        write (*,*) 'Error in subroutine get_pvlim !'
        return
     endif
    !  gets pv limits and latitudes in descending order !!!
  enddo
  
  ! interpolate PV field onto given trajectory locations:
  ! latitudes are stored in decreasing order !!!
  ! and longitudes in increasing order
  DO i=1, npoints
     IF (lons(startindex+i-1) < lon_grid(1)) llons(i)=lons(startindex+i-1)+360.
     IF (lats(startindex+i-1) < minval(lat_grid(1:nlats))) llats(i)=minval(lat_grid(1:nlats))
     IF (lats(startindex+i-1) > maxval(lat_grid(1:nlats))) llats(i)=maxval(lat_grid(1:nlats))
     
     DO j=2,nlats
        ilat = j-1
        dlat = (llats(i)-lat_grid(ilat))/(lat_grid(ilat+1)-lat_grid(ilat))
        IF (dlat >= 0.0 .and. dlat <= 1.0) EXIT
     ENDDO
     IF (dlat > 1. .or. dlat < 0.) THEN
        write(*,*) 'Latitude can not be interpolated! program stop'
        write(*,*)'dlat=',dlat
        status = -1
        return
     ENDIF
     ! modulo calculation because of different conventions of longitude range...
     DO j=1,nlons
        ilon = j
        ilonp1 = ilon+1
        ! whether between the first and lasr longitude grid point the data need
        ! to be interpolated
        if (ilon==nlons) ilonp1=1

! op_pj_20160606+
!!$        dlon = mod(720.+ llons(i) - lon_grid(ilon),360.)/ &
!!$               mod(lon_grid(ilonp1)-lon_grid(ilon)+360.,360.)
        dlon = mod(720.+ llons(i) - lon_grid(ilon),360._prec)/ &
               mod(lon_grid(ilonp1)-lon_grid(ilon)+360._prec,360._prec)
! op_pj_20160606-
        IF (dlon >= 0 .and. dlon <= 1.) EXIT
     ENDDO
     IF (dlon > 1. .or. dlon < 0.) THEN
        write(*,*) 'Longitude can not be interpolated! program stop'
        write(*,*)'dlon=',dlon
        status = -1
        return
     ENDIF
     ! write(*,*)i,' lat ',ilat,lat_grid(ilat:ilat+1),llats(i),dlat
     ! write(*,*)i,' lon ',ilon,lon_grid(ilon),lon_grid(ilonp1),llons(i),dlon
     ! Interpolation of data on a  2-dim lat/lon grid not definite
     ! !here interpolate first in latitude
     ! write(*,*)
     ! write(*,*)lon_grid(ilon),lon_grid(ilonp1)
     ! write(*,*)'-------------------------'
     ! write(*,*)PV_field(ilon,ilat+1,1),PV_field(ilonp1,ilat+1,1),'|',lat_grid(ilat+1)
     ! write(*,*)PV_field(ilon,ilat,1),PV_field(ilonp1,ilat,1),'|',lat_grid(ilat)
     ! write(*,*)
     
     pv_int(i) = mdi
     if (abs ((pv_field(ilon,  ilat,  lev_index(i))-mdi)/mdi) <= eps .and. dlat /= 1.) cycle
     if (abs ((pv_field(ilonp1,ilat,  lev_index(i))-mdi)/mdi) <= eps .and. dlat /= 1.) cycle
     if (abs ((pv_field(ilon,  ilat+1,lev_index(i))-mdi)/mdi) <= eps .and. dlat /= 0.) cycle
     if (abs ((pv_field(ilonp1,ilat+1,lev_index(i))-mdi)/mdi) <= eps .and. dlat /= 0.) cycle
     if (dlat == 1.) then 
        x1=pv_field(ilon,  ilat+1,lev_index(i))
        x2=pv_field(ilonp1,ilat+1,lev_index(i))
     else
        x1=pv_field(ilon,ilat,lev_index(i))+&
             dlat*(pv_field(ilon,ilat+1,lev_index(i))-pv_field(ilon,ilat,lev_index(i)))
        x2=pv_field(ilonp1,ilat,lev_index(i))+&
             dlat*(pv_field(ilonp1,ilat+1,lev_index(i))-pv_field(ilonp1,ilat,lev_index(i)))
     endif
     pv_int(i) = x1+dlon*(x2-x1)
     !write(*,*) rank,'i,x1,x2,pv',i,x1,x2,pv_int(i)
  ENDDO
  
  DO i=1,npoints
     ! find closest PV value in grid
     DO na = 2,nlats
        n = na
        IF (pv_int(i) >= pv_limits(na,lev_index(i))) EXIT
     ENDDO
     ! interpolate eq. latitude to given PV value
     dpv      = (pv_int(i)-pv_limits(n,lev_index(i))) / &
                (pv_limits(n-1,lev_index(i))-pv_limits(n,lev_index(i)))
     eqlat(startindex+i-1) = lat_limits(n) + dpv*(lat_limits(n-1)-lat_limits(n))
     if (abs((pv_int(i)-mdi)/mdi) <= eps) eqlat(startindex+i-1) = mdi
  ENDDO
  !write (*,*) 'end get_eqlat_int'

END subroutine get_eqlat_int



End Module messy_clams_tools_eqlatutils

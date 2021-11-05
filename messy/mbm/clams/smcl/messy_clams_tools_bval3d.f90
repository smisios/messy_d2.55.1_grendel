!******************************************************************************
! File utils/src/lib_bval3d.f90
!
! This file is part of the Chemical Lagrangian Model of the
! Stratosphere (CLaMS). CLaMS is a hierarchy of numerical models and
! preprocessors which simulate transport, chemical reactions and
! mixing processes in the stratosphere.  
!
! Copyright (c) 2006
! J. Ankenbrand, N. Thomas, Forschungszentrum Juelich GmbH
! Last Modified By: N.Thomas
! Last Modified On: Wed Jun  7 13:55:43 2017
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
! Module lib_bval3d
!
! This module contains some interpolation tools:
!
!   function distance(p1,p2)
!   function area(p1,p2,p3)
!   function no_area(p1,p2,p,param,level,index_lat,index_lon,lgnpt,mdi,eps)
!   function bval3d (data,level,plat,plon,lons,lats,nlon,nlat,latint,rmdi,eps)
!   function bval3d_new (val1,val2,val3,val4,plat,plon, &
!                        lons,lats,nlon,nlat,latint,rmdi,eps)
!
!********************************************************************
                         
MODULE  messy_clams_tools_bval3d

CONTAINS


!**************************************************************
! function to get the distance of 2 points in sphere coordinates
!**************************************************************
function distance(p1,p2)

  use messy_clams_global, only: dp

  implicit none

  real(DP),dimension(3) :: p1,p2
  real(DP) :: distance
  
  distance=asin(sqrt((p1(1)-p2(1))**2+(p1(2)-p2(2))**2+(p1(3)-p2(3))**2))

end function distance


!**************************************************************
! function to get the area of a sphere triangle
!**************************************************************
function area(p1,p2,p3)

  use messy_clams_global, only: dp

  implicit none

  real(DP),dimension(3) :: p1,p2,p3
  real(DP) :: area

  real(DP) :: d_p1_p2,d_p1_p3,d_p2_p3
  real(DP) :: alpha,beta,gamma

  real(DP),parameter :: pi=3.14159265359

  ! get the distances of the 3 points
  d_p1_p2 =  distance(p1,p2)
  d_p1_p3 =  distance(p1,p3)
  d_p2_p3 =  distance(p2,p3)

  ! Seitencosinussatz
  alpha=acos((cos(d_p1_p2)-cos(d_p1_p3)*cos(d_p2_p3))/(sin(d_p1_p3)&
       &*sin(d_p2_p3)))
  beta=acos((cos(d_p1_p3)-cos(d_p1_p2)*cos(d_p2_p3))/(sin(d_p1_p2)&
       &*sin(d_p2_p3)))
  gamma=acos((cos(d_p2_p3)-cos(d_p1_p2)*cos(d_p1_p3))/(sin(d_p1_p2)&
       &*sin(d_p1_p3)))

  ! get area from angles
  area=alpha+beta+gamma-pi
  
  if (area <= 0) then 
     write(*,*) 'error: square triangle area not >0'
     write(*,*) 'please choose a greater epsilon in function bval3d!'
  end if

end function area

!**************************************************************
! function to interpolate between 2 points
!**************************************************************
function no_area(p1,p2,p,val1,val2)

  use messy_clams_global, only: dp, prec, mdi, eps

  implicit none 

  real(DP),dimension(3) :: p1,p2,p
  real(prec)    :: val1, val2
  
  real(prec) :: no_area

  ! weight function:
  ! p1p2 / p1p  +  p1p2 / p2p

  if (abs((val1-mdi)/mdi)<=eps .or. abs((val2-mdi)/mdi)<=eps) then

     no_area = mdi

  else
     
     no_area = (distance(p1,p2)/distance(p1,p))/ &
       (distance(p1,p2)/distance(p1,p)+distance(p1,p2)/distance(p2,p)) * val1 + &
       (distance(p1,p2)/distance(p2,p))/ &
       (distance(p1,p2)/distance(p1,p)+distance(p1,p2)/distance(p2,p)) * val2 

  endif

end function no_area

!********************************************************************
! function to interpolate a point from its grid neighbours on a sphere
! the function considers, that the latitude grid is not constant
!********************************************************************
function bval3d (data,level,plat,plon,longrid,latgrid,nlon,nlat,lat_asc)

  use messy_clams_global, only: dp, prec, mdi,eps

  IMPLICIT NONE
 
  real(prec) :: bval3d

  real(prec),dimension(:,:,:),pointer :: data
  real(prec),dimension(:),    pointer :: longrid, latgrid
  integer                             :: level       ! index of theta level
  integer                             :: nlon, nlat
  logical                             :: lat_asc     ! true:  latitudes ascending
                                                     ! false: latitudes descending
  real(prec)                          :: plat,plon   ! 

  integer :: index_lon,index_lonp1,index_lat, index_latp1

  real(prec)    :: val1, val2, val3, val4

  integer :: i

  !get indices of the actual point in the longitude and latitude array
!write (*,*) 'in bval3d'
  index_lon= -1
  index_lonp1 = -1
  index_lat = -1 
  index_latp1 = -1

! op_pj_20160606
!!$  if (plon>360.) plon = mod(plon,360.)
  if (plon>360._prec) plon = mod(plon,360._prec)
! op_pj_20160606-
  if (plon<longrid(1) .or. plon>=longrid(nlon)) then
     index_lon   = nlon
     index_lonp1 = 1
  else
     do i = 1, nlon-1
        if (longrid(i)<=plon .and. plon<=longrid(i+1)) then 
           index_lon   = i
           index_lonp1 = i+1
        endif
     enddo
  endif

  if (lat_asc) then  ! latitudes in ascending order
     if (plat<latgrid(1)) then
        index_lat   = 1
        index_latp1 = 1
     elseif (plat>=latgrid(nlat)) then
        index_lat   = nlat
        index_latp1 = nlat
     else        
        do i = 1, nlat-1
           if (latgrid(i)<=plat .and. plat<=latgrid(i+1)) then
              index_lat   = i
              index_latp1 = i+1
           endif
        enddo
     endif
  else ! latitudes in descending order
     if (plat>latgrid(1)) then
        index_lat   = 1
        index_latp1 = 1
     elseif (plat<=latgrid(nlat)) then
        index_lat   = nlat
        index_latp1 = nlat
     else   
        do i = 1, nlat-1
           if (latgrid(i)>=plat .and. plat>=latgrid(i+1)) then
              index_lat   = i
              index_latp1 = i+1
           endif
        enddo
     endif
  endif

  val1 = data(index_lon,   index_lat, level )
  val2 = data(index_lonp1, index_lat, level)
  if (index_lat==index_latp1) then  ! near north or south pole
     val3 = data(mod(index_lon  +nlon/2,nlon), index_lat, level)
     val4 = data(mod(index_lonp1+nlon/2,nlon), index_lat, level)
  else
     val3 = data(index_lonp1, index_latp1, level)
     val4 = data(index_lon,   index_latp1, level )
  endif

  bval3d = bval3d_new (val1,val2,val3,val4,plat,plon,longrid,latgrid,nlon,nlat,lat_asc)

end function bval3d

!********************************************************************
! function to interpolate a point from its grid neighbours on a sphere
! the function considers, that the latitude grid is not constant
!********************************************************************
function bval3d_new (val1,val2,val3,val4,plat,plon,longrid,latgrid,nlon,nlat,lat_asc)

  use messy_clams_global, only: dp, prec, mdi, eps

  IMPLICIT NONE
 
  real(prec) :: bval3d_new

  real(prec),dimension(:),pointer :: longrid, latgrid
  real(prec)                      :: val1,val2,val3,val4
  real(prec)                      :: plat,plon             ! 
  integer                         :: nlon, nlat
  logical                         :: lat_asc    ! true:  latitudes ascending
                                                ! false: latitudes descending

  real(DP),dimension(3) :: p1,p2,p3,p4,p,p12,p23,p34,p41
  real(DP) :: area1,area2,area3,area4,area_all,w
  real(prec) :: lon_start,lat_start,lon_end,lat_end
  real(prec) :: epslat

  real(DP),parameter :: pi=3.14159265359

  integer :: i
  logical :: latbound = .false.

  epslat=5.e-04  ! 0.05 latitude difference means an epsilon of 0.00087

  ! get indices of the actual point in the longitude and latitude array

! op_pj_20160606+
!!$  if (plon>360.) plon = mod(plon,360.)
  if (plon>360._prec) plon = mod(plon,360._prec)
! op_pj_20160606-
  if (plon<longrid(1) .or. plon>longrid(nlon)) then
     lon_start = longrid(nlon)
     lon_end   = longrid(1)
  else
     do i = 1, nlon-1
        if (longrid(i)<=plon .and. plon<=longrid(i+1)) then 
           lon_start = longrid(i)
           lon_end   = longrid(i+1)
        endif
     enddo
  endif

  latbound = .false.
  if (lat_asc) then  ! latitudes in ascending order
     if (plat<latgrid(1)) then
        lat_start = latgrid(1)
        lat_end   = latgrid(1)
        latbound  = .true.
     elseif (plat>=latgrid(nlat)) then
        lat_start = latgrid(nlat)
        lat_end   = latgrid(nlat)
        latbound  = .true. 
     else        
        do i = 1, nlat-1
           if (latgrid(i)<=plat .and. plat<=latgrid(i+1)) then
              lat_start = latgrid(i)
              lat_end   = latgrid(i+1)
           endif
        enddo
     endif
  else ! latitudes in descending order
     if (plat>latgrid(1)) then
        lat_start = latgrid(1)
        lat_end   = latgrid(1)
        latbound  = .true.
     elseif  (plat<=latgrid(nlat)) then
        lat_start = latgrid(nlat)
        lat_end   = latgrid(nlat)
        latbound  = .true.
     else        
        do i = 1, nlat-1
           if (latgrid(i)>=plat .and. plat>=latgrid(i+1)) then
              lat_start = latgrid(i)
              lat_end   = latgrid(i+1)
           endif
        enddo
     endif
  endif

  !   p1      p12            p2
  !  
  !       A1          A2
  ! 
  !  p41       p            p23
  ! 
  ! 
  !      A4           A3
  ! 
  !   p4      p34            p3
  
!!!! Bei nicht globalen DS muss dies nicht bedeuten, dass man in der Naehe
!!!! eines Pols ist !!!!
!  if ((plat>=latgrid(1) .or. (plat<=lat0+(nlat-1)*latint)) then
  if (latbound) then
 
     !!! p3, p4, p34 geaendert !!!
     !write(*,*) ' p3, p4, p34 geaendert !!!'

!!!!! val3 und val4 muessten auch von der gegenueberliegenden Seite genommen werden
!!!!! (latstart,lon_start+180.) !!!
     
     ! get sphere coordinates of p1,p2,p3,p4 
     p1(1) = cos(lon_start*(pi/180.d0))* sin(abs(lat_start*(pi/180.d0)-pi/2.d0))
     p1(2) = sin(lon_start*(pi/180.d0))* sin(abs(lat_start*(pi/180.d0)-pi/2.d0))
     p1(3) = cos(abs(lat_start*(pi/180.d0)-pi/2.d0))
     p2(1) = cos((lon_end)*(pi/180.d0))* sin(abs(lat_start*(pi/180.d0)-pi/2.d0))
     p2(2) = sin((lon_end)*(pi/180.d0))* sin(abs(lat_start*(pi/180.d0)-pi/2.d0))
     p2(3) = cos(abs(lat_start*(pi/180.d0)-pi/2.d0))
! op_pj_20160606+
!!$     p3(1) = cos((MOD(lon_start+180.,360.))*(pi/180.d0))* sin(abs((lat_start)*(pi/180.d0)-pi/2.d0))
!!$     p3(2) = sin((MOD(lon_start+180.,360.))*(pi/180.d0))* sin(abs((lat_start)*(pi/180.d0)-pi/2.d0))
!!$     p3(3) = cos(abs((lat_start)*(pi/180.d0)-pi/2.d0))
!!$     p4(1) = cos(MOD(lon_end+180.,360.)*(pi/180.d0))* sin(abs((lat_start)*(pi/180.d0)-pi/2.d0))
!!$     p4(2) = sin(MOD(lon_end+180.,360.)*(pi/180.d0))* sin(abs((lat_start)*(pi/180.d0)-pi/2.d0))
     p3(1) = cos((MOD(lon_start+180._prec,360._prec))*(pi/180.d0))* sin(abs((lat_start)*(pi/180.d0)-pi/2.d0))
     p3(2) = sin((MOD(lon_start+180._prec,360._prec))*(pi/180.d0))* sin(abs((lat_start)*(pi/180.d0)-pi/2.d0))
     p3(3) = cos(abs((lat_start)*(pi/180.d0)-pi/2.d0))
     p4(1) = cos(MOD(lon_end+180._prec,360._prec)*(pi/180.d0))* sin(abs((lat_start)*(pi/180.d0)-pi/2.d0))
     p4(2) = sin(MOD(lon_end+180._prec,360._prec)*(pi/180.d0))* sin(abs((lat_start)*(pi/180.d0)-pi/2.d0))
! op_pj_20160606-
     p4(3) = cos(abs((lat_start)*(pi/180.d0)-pi/2.d0))
     
     ! get sphere coordinates of p12,p23,p34 and p41
     p12(1) = cos(plon*(pi/180.d0))* sin(abs(lat_start*(pi/180.d0)-pi/2.d0))
     p12(2) = sin(plon*(pi/180.d0))* sin(abs(lat_start*(pi/180.d0)-pi/2.d0))
     p12(3) = cos(abs(lat_start*(pi/180.d0)-pi/2.d0))
     p23(1) = cos((lon_end)*(pi/180.d0))* sin(abs(plat*(pi/180.d0)-pi/2.d0))
     p23(2) = sin((lon_end)*(pi/180.d0))* sin(abs(plat*(pi/180.d0)-pi/2.d0))
     p23(3) = cos(abs(plat*(pi/180.d0)-pi/2.d0))
! op_pj_20160606+
!!$     p34(1) = cos(MOD(plon+180.,360.)*(pi/180.d0))* sin(abs(lat_start*(pi/180.d0)-pi/2.d0))
!!$     p34(2) = sin(MOD(plon+180.,360.)*(pi/180.d0))* sin(abs(lat_start*(pi/180.d0)-pi/2.d0))
     p34(1) = cos(MOD(plon+180._prec,360._prec)*(pi/180.d0))* sin(abs(lat_start*(pi/180.d0)-pi/2.d0))
     p34(2) = sin(MOD(plon+180._prec,360._prec)*(pi/180.d0))* sin(abs(lat_start*(pi/180.d0)-pi/2.d0))
! op_pj_20160606-
     p34(3) = cos(abs(lat_start*(pi/180.d0)-pi/2.d0))
     p41(1) = cos((lon_start)*(pi/180.d0))* sin(abs(plat*(pi/180.d0)-pi/2.d0))
     p41(2) = sin((lon_start)*(pi/180.d0))* sin(abs(plat*(pi/180.d0)-pi/2.d0))
     p41(3) = cos(abs(plat*(pi/180.d0)-pi/2.d0))
     
  else
     ! get sphere coordinates of p1,p2,p3,p4 
     p1(1) = cos(lon_start*(pi/180.d0))* sin(abs(lat_start*(pi/180.d0)-pi/2.d0))
     p1(2) = sin(lon_start*(pi/180.d0))* sin(abs(lat_start*(pi/180.d0)-pi/2.d0))
     p1(3) = cos(abs(lat_start*(pi/180.d0)-pi/2.d0))
     p2(1) = cos((lon_end)*(pi/180.d0))* sin(abs(lat_start*(pi/180.d0)-pi/2.d0))
     p2(2) = sin((lon_end)*(pi/180.d0))* sin(abs(lat_start*(pi/180.d0)-pi/2.d0))
     p2(3) = cos(abs(lat_start*(pi/180.d0)-pi/2.d0))
     p3(1) = cos((lon_end)*(pi/180.d0))* sin(abs((lat_end)*(pi/180.d0)-pi/2.d0))
     p3(2) = sin((lon_end)*(pi/180.d0))* sin(abs((lat_end)*(pi/180.d0)-pi/2.d0))
     p3(3) = cos(abs((lat_end)*(pi/180.d0)-pi/2.d0))
     p4(1) = cos(lon_start*(pi/180.d0))* sin(abs((lat_end)*(pi/180.d0)-pi/2.d0))
     p4(2) = sin(lon_start*(pi/180.d0))* sin(abs((lat_end)*(pi/180.d0)-pi/2.d0))
     p4(3) = cos(abs((lat_end)*(pi/180.d0)-pi/2.d0))
     
     ! get sphere coordinates of p12,p23,p34 and p41
     p12(1) = cos(plon*(pi/180.d0))* sin(abs(lat_start*(pi/180.d0)-pi/2.d0))
     p12(2) = sin(plon*(pi/180.d0))* sin(abs(lat_start*(pi/180.d0)-pi/2.d0))
     p12(3) = cos(abs(lat_start*(pi/180.d0)-pi/2.d0))
     p23(1) = cos((lon_end)*(pi/180.d0))* sin(abs(plat*(pi/180.d0)-pi/2.d0))
     p23(2) = sin((lon_end)*(pi/180.d0))* sin(abs(plat*(pi/180.d0)-pi/2.d0))
     p23(3) = cos(abs(plat*(pi/180.d0)-pi/2.d0))
     p34(1) = cos(plon*(pi/180.d0))* sin(abs((lat_end)*(pi/180.d0)-pi/2.d0))
     p34(2) = sin(plon*(pi/180.d0))* sin(abs((lat_end)*(pi/180.d0)-pi/2.d0))
     p34(3) = cos(abs((lat_end)*(pi/180.d0)-pi/2.d0))
     p41(1) = cos((lon_start)*(pi/180.d0))* sin(abs(plat*(pi/180.d0)-pi/2.d0))
     p41(2) = sin((lon_start)*(pi/180.d0))* sin(abs(plat*(pi/180.d0)-pi/2.d0))
     p41(3) = cos(abs(plat*(pi/180.d0)-pi/2.d0))
  endif
  
  ! get spere coordinates of p
  p(1) = cos(plon*(pi/180.d0))* sin(abs(plat*(pi/180.d0)-pi/2.d0))
  p(2) = sin(plon*(pi/180.d0))* sin(abs(plat*(pi/180.d0)-pi/2.d0))
  p(3) = cos(abs(plat*(pi/180.d0)-pi/2.d0))
  
  if (distance(p,p12) <= epslat) then 
     ! p is between p1 and p2
     if (distance(p12,p1) <= epslat) then 
        ! p is very near p1
        bval3d_new = val1 
     else if (distance(p12,p2) <= epslat) then 
        !p is very near p2
        bval3d_new = val2 
     else 
        ! interpolate between p1 and p2
        bval3d_new = no_area (p1,p2,p,val1,val2)
     endif
  else if (distance(p,p34) <= epslat) then 
     ! p is between p3 and p4
     if (distance(p34,p4) <= epslat) then 
        ! p is very near p4
        bval3d_new = val4 
     else if (distance(p34,p3) <= epslat) then 
        ! p is very near p3
        bval3d_new = val3 
     else 
        ! interpolate between p3 and p4
        bval3d_new = no_area (p3,p4,p,val3,val4)
     endif
  else if (distance(p,p41) <= epslat) then 
     ! p is between p1 and p4 ; interpolate
     bval3d_new = no_area (p1,p4,p,val1,val4)
  else if (distance(p,p23) <= epslat) then 
     ! p is between p2 and p3 ; interpolate
     bval3d_new = no_area (p2,p3,p,val2,val3)
  else
     
     if (abs((val1-mdi)/mdi)<=eps .or. abs((val2-mdi)/mdi)<=eps .or. &
          abs((val3-mdi)/mdi)<=eps .or. abs((val4-mdi)/mdi)<=eps) then
        
        bval3d_new=mdi
        
     else
        
        if (abs(lat_start-90.) <= 0.01) then
           ! at north pole
           !write (*,*) 'north pole'
           area1 = area(p,p1,p41) ! area of triangle p1,p,p41
           area2 = area(p,p2,p23) ! area of triangle p2,p23,p
           area3 = area(p,p3,p23) + area(p,p3,p34) ! area of quadrangle p23,p3,p34,p
           area4 = area(p,p4,p34) + area(p,p4,p41) ! area of quadrangle p34,p4,p41,p
           
        else if (abs(lat_end+90.) <= 0.01) then
           ! at south pole
           !write (*,*) 'south pole'
           area1 = area(p,p1,p41) + area(p,p1,p12) ! area of quadrangle p1,p12,p,p41
           area2 = area(p,p2,p12) + area(p,p2,p23) ! area of quadrangle p12,p2,p23,p
           area3 = area(p,p3,p23)  ! area of triangle p23,p3,p
           area4 = area(p,p4,p41) ! area of triangle p4,p41,p
        else
           
           ! p is anywhere else between p1,p2,p3 and p4
           area1 = area(p,p1,p41) + area(p,p1,p12) ! area of quadrangle p1,p12,p,p41
           area2 = area(p,p2,p12) + area(p,p2,p23) ! area of quadrangle p12,p2,p23,p
           area3 = area(p,p3,p23) + area(p,p3,p34) ! area of quadrangle p23,p3,p34,p
           area4 = area(p,p4,p34) + area(p,p4,p41) ! area of quadrangle p34,p4,p41,p
           
        endif
        
        area_all=area1+area2+area3+area4 ! area of quadrangle p1,p2,p3,p4
        w=area_all/area1+area_all/area2+area_all/area3+area_all/area4 !weight function
        
        ! interpolate p from p1,p2,p3 and p4 in dependence of area1,area2,area3
        ! and area4
        bval3d_new = real( &
             ((area_all/area1)/w)*val1 + &
             ((area_all/area2)/w)*val2 + &
             ((area_all/area3)/w)*val3 + &
             ((area_all/area4)/w)*val4 )
     end if
     
  endif

end function bval3d_new



END MODULE messy_clams_tools_bval3d







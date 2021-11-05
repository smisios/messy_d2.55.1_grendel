!********************************************************************************!
! File clams/traj/source/calc3d.f90
!
! This file is part of the Chemical Lagrangian Model of the
! Stratosphere (CLaMS). CLaMS is a hierarchy of numerical models and
! preprocessors which simulate transport, chemical reactions and
! mixing processes in the stratosphere.  
!
! Copyright (c) 2006
! Richard Swinbank, Daniel McKenna, Nicole Thomas, Paul Konopka, 
! Forschungszentrum Juelich GmbH
! Last Modified By: N.Thomas
! Last Modified On: Wed Jan 15 13:33:26 2020
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
!
! MODULE calc3d
! -------------
!
! This MODULE contains (numerical) subroutines:
!
! SUBROUTINE rungek3d (plons, plats, plevs, calc, &
!                      u0   ,  v0   ,  w0  , &
!                      umid ,  vmid ,  wmid, & 
!                      udt  ,  vdt  ,  wdt , &              
!                      nx   ,  ny   ,  nz  , & 
!                      longrid,latgrid, levelgrid,    &
!                      nparts,dt,pslat,loglev)     
! subroutine get_lev_box (level, plev, nz, ix, iy, iz, rz, ascending)
! SUBROUTINE intruv3d (u    , v     , w     , &
!                      plons, plats , plevs , &
!                      pu   , pv    , pw    , & 
!                      nx   , ny    , nz    , &
!                      levelgrid,    &
!                      nparts,er,indxps,loglev)
! SUBROUTINE ptoll (plons,plats,indxps,nparts,pi)  
! SUBROUTINE fixll (plons,plats,nparts,pi)  
! SUBROUTINE zen2(zlat,zlong,idum,time,cozen,SZA)               
! CONTAINS: FUNCTION DECLINE(iy,Im,ID)   
!********************************************************************************

MODULE messy_clamstraj_calc3d

CONTAINS


!**************************************************************
!  Runge kutta scheme to advect particles in three dimensions
!**************************************************************
subroutine rungek3d (plons, plats, plevs, calc, tsv, &
                            u0,    v0,    w0, level0, dlevdz0, &
                          umid,  vmid,  wmid, levelmid, dlevdzmid,  & 
                           udt,   vdt,   wdt, leveldt, dlevdzdt, &              
                            nx,    ny,    nz, & 
                           longrid, latgrid, levelgrid,    &
                          nparts,loop_start,loop_end, &
                          dt,pslat,loglev,use_dlevdz)    
  
!
!...Translated by Pacific-Sierra Research 77to90  4.3E  15:09:03   7/16/98  
!...Switches:                     

  use messy_clams_global,     only: prec, rank, mdi, eps

      implicit none

!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer     :: nx, ny, nz, nparts, loop_start, loop_end
      real(prec) , intent(in) :: dt 
      real(prec) , intent(in) :: pslat 
      real(prec)  :: plons(nparts), plats(nparts), plevs(nparts) 
      real(prec)  :: tsv(nparts)
      real(prec)  :: u0(nx,ny,nz), v0(nx,ny,nz), w0(nx,ny,nz)
      real(prec)  :: level0(nx,ny,nz), dlevdz0(nx,ny,nz)
      real(prec)  :: umid(nx,ny,nz), vmid(nx,ny,nz), wmid(nx,ny,nz)
      real(prec)  :: levelmid(nx,ny,nz), dlevdzmid(nx,ny,nz)
      real(prec)  :: udt(nx,ny,nz), vdt(nx,ny,nz), wdt(nx,ny,nz)
      real(prec)  :: leveldt(nx,ny,nz), dlevdzdt(nx,ny,nz)
      real(prec)  :: longrid(1:nx+1) !(nx)
      real(prec)  :: latgrid(0:ny+1) !(ny)
      real(prec)  :: levelgrid(nz) 
      logical, intent(in) :: calc(nparts) 
      logical, intent(in) :: loglev
      logical, optional   :: use_dlevdz

!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer , dimension(nparts) :: indxps 
      integer :: i 
      real(prec) , dimension(nparts) :: pu1, pv1, pw1, pu2, pv2, pw2, pu3, pv3, pw3, &
                                        xwk, ywk, zwk
      real(prec) , dimension(nparts) :: pdlevdz1, pdlevdz2, pdlevdz3
      real(prec) :: pi 
      
      logical :: use_dldz = .false.
!-----------------------------------------------
!
!  Inputs
! plons, plats, plevs         particle coordinates
!    u0,    v0,    w0         winds at time t
!  umid,  vmid,  wmid         winds at time t+dt/2
!   udt,   vdt,   wdt         winds at time t+dt
!    nx,    ny,    nz         grid size (lons,lats,levs)
!              nparts         number of particles
! longrid,latgrid,levelgrid   grid positions
!  dlon,  dlat                grid spacings
!               pslat         latitude to switch to polar stereo graphic grid
!                  dt         time step
 
!  Outputs
! plons, plats, plevs         new positions of particles at t+dt
 
!  Workspace
! pu1,pv1,pw1,         intermediate wind values
! pu2,pv2,pw2,         intermediate wind values
! pu3,pv3,pw3          intermediate wind values
! xwk,ywk,zwk          intermediate positions
! rercos   ps grid:    one over the distance to the Earth's axis in metres
! coslon   ps grid:    cosine of longitude
! sinlon   ps grid:    sine of longitude
! fmap     ps grid:    two over one plus the sine of absolute latitude
! indxps   ps grid:    index of which particles to put on ps grid
      
      use_dldz = .false.
      if (present(use_dlevdz)) use_dldz = use_dlevdz

      pi = 4.*atan(1.) 

   !decide which coordinate system to use for each particle

      indxps = 0
!      DO i = 1,nparts
      DO i = loop_start, loop_end
         if (plats(i) > pslat) then 
            indxps(i) = 1 
         else if (plats(i) < (-pslat)) then 
            indxps(i) = -1 
         else 
            indxps(i) = 0 
         endif 
      end do    

   ! interpolate winds onto particles at 'r(t)' and convert to
   ! polar-stereo where required.

   call intruv3d (u0,   v0,   w0, dlevdz0, level0, &
        plons,plats,plevs,  &   
        pu1,  pv1,  pw1, pdlevdz1, &  
        nx,   ny,   nz,   &         
        longrid, latgrid, levelgrid,   &
        nparts, loop_start, loop_end, indxps, loglev)    

   ! calculate 'r(t)+(k1)/2'
   !do i = 1, nparts
   do i = loop_start, loop_end
      if (ABS((pu1(i)-mdi)/mdi)<=eps) then
         xwk(i) = mdi 
         ywk(i) = mdi 
         zwk(i) = mdi 
      else
         xwk(i) = plons(i) + 0.5*dt*pu1(i) 
         ywk(i) = plats(i) + 0.5*dt*pv1(i) 
         if (use_dldz) then
            zwk(i) = plevs(i) + 0.5*dt*pdlevdz1(i)*tsv(i)
         else
            zwk(i) = plevs(i) + 0.5*dt*pw1(i) 
         endif
      endif
   end do
   
!!!!! Feldoperationen nutzen ????? => Zeiten vergleichen !!!
!!$   if (use_dldz) then
!!$      where (ABS((pu1(:nparts)-mdi)/mdi)<=eps)  
!!$         xwk(:nparts) = mdi 
!!$         ywk(:nparts) = mdi 
!!$         zwk(:nparts) = mdi 
!!$      elsewhere 
!!$         xwk(:nparts) = plons(:nparts) + 0.5*dt*pu1(:nparts) 
!!$         ywk(:nparts) = plats(:nparts) + 0.5*dt*pv1(:nparts) 
!!$         zwk(:nparts) = plevs(:nparts) + 0.5*dt*pdthetadz1(:nparts)*tsv(:nparts)
!!$      end where
!!$   else
!!$      where (ABS((pu1(:nparts)-mdi)/mdi)<=eps)  
!!$         xwk(:nparts) = mdi 
!!$         ywk(:nparts) = mdi 
!!$         zwk(:nparts) = mdi 
!!$      elsewhere 
!!$         xwk(:nparts) = plons(:nparts) + 0.5*dt*pu1(:nparts) 
!!$         ywk(:nparts) = plats(:nparts) + 0.5*dt*pv1(:nparts) 
!!$         zwk(:nparts) = plevs(:nparts) + 0.5*dt*pw1(:nparts) 
!!$      end where
!!$   endif


   ! convert back to longitudes and latitudes
      call ptoll (xwk, ywk, indxps, nparts, loop_start, loop_end, pi) 

   ! interpolate winds onto particles at 'r(t)+(k1)/2' and convert to
   ! polar-stereo where required
   call intruv3d (umid,vmid,wmid,dlevdzmid,levelmid, &
        xwk,   ywk,   zwk, &
        pu2,   pv2,   pw2, pdlevdz2, &
        nx,    ny,    nz,  &
        longrid, latgrid, levelgrid,    &
        nparts,loop_start,loop_end,indxps,loglev) 
 
   ! calculate 'r(t)+(k2)/2'
   !do i = 1, nparts
   do i = loop_start, loop_end
      if (ABS((pu2(i)-mdi)/mdi)<=eps) then
         xwk(i) = mdi 
         ywk(i) = mdi 
         zwk(i) = mdi 
      else
         xwk(i) = plons(i) + 0.5*dt*pu2(i) 
         ywk(i) = plats(i) + 0.5*dt*pv2(i) 
         if (use_dldz) then
            zwk(i) = plevs(i) + 0.5*dt*pdlevdz2(i)*tsv(i)
        else
            zwk(i) = plevs(i) + 0.5*dt*pw2(i) 
         endif
      endif
   enddo 

   ! convert back to longitudes and latitudes
      call ptoll (xwk, ywk, indxps, nparts, loop_start, loop_end, pi) 

   ! interpolate winds onto particles at 'r(t)+(k2)/2' and convert to
   ! polar-stereo where required
   call intruv3d (umid,vmid,wmid,dlevdzmid,levelmid, &
        xwk,   ywk,   zwk, &  
        pu3,   pv3,   pw3, pdlevdz3, &  
        nx,    ny,    nz,  &  
        longrid, latgrid, levelgrid,    &
        nparts,loop_start,loop_end,indxps,loglev)    

   ! calculate 'r(t)+k3'
   !do i = 1, nparts
   do i = loop_start, loop_end
     if (ABS((pu3(i)-mdi)/mdi)<=eps) then
         xwk(i) = mdi 
         ywk(i) = mdi 
         zwk(i) = mdi 
      else
         xwk(i) = plons(i) + dt*pu3(i) 
         ywk(i) = plats(i) + dt*pv3(i) 
         if (use_dldz) then
            zwk(i) = plevs(i) + dt*pdlevdz3(i)*tsv(i)
         else
            zwk(i) = plevs(i) + dt*pw3(i) 
         endif
         pu3(i) = pu2(i) + pu3(i) 
         pv3(i) = pv2(i) + pv3(i) 
         pw3(i) = pw2(i) + pw3(i) 
         if (use_dldz) pdlevdz3(i) = pdlevdz2(i) + pdlevdz3(i)
      endif
   enddo 

   ! convert back to longitudes and latitudes
      call ptoll (xwk, ywk, indxps, nparts, loop_start, loop_end, pi) 

   ! interpolate winds onto particles at 'r(t)+k3' and convert t
   ! polar-stereo where required
   call intruv3d (udt,  vdt, wdt, dlevdzdt, leveldt, &       
        xwk,  ywk, zwk, &            
        pu2,  pv2, pw2, pdlevdz2, &          
        nx,   ny,   nz,   &          
        longrid, latgrid, levelgrid,   &
        nparts,loop_start,loop_end,indxps,loglev) 

   ! calculate particle positions at end of step 'r(t+dt)'
   !Do i = 1, nparts
   Do i = loop_start, loop_end
      if (ABS((pu2(i)-mdi)/mdi)<=eps) then 
         plons(i) = mdi 
         plats(i) = mdi 
         plevs(i) = mdi 
      elseif (calc(i)) then 
         plons(i) = plons(i) + (dt/6.)*(pu1(i)+pu2(i)+2.*pu3(i)) 
         plats(i) = plats(i) + (dt/6.)*(pv1(i)+pv2(i)+2.*pv3(i)) 
         if (use_dldz) then
            plevs(i) = plevs(i) + (dt/6.)*(pdlevdz1(i)+pdlevdz2(i)+2.*pdlevdz3(i))*tsv(i) 
         else
            plevs(i) = plevs(i) + (dt/6.)*(pw1(i)+pw2(i)+2.*pw3(i)) 
         endif
      endif
   end do

   ! convert back to longitudes and latitudes and make sure they are
   ! in range  (0,2*pi)  and  (-.5*pi,.5*pi)

   call ptoll (plons, plats, indxps, nparts, loop_start, loop_end, pi) 
   call fixll (plons, plats, nparts, loop_start, loop_end, pi) 

   return  
   
 end subroutine rungek3d
 
 
!**************************************************************
! 
!**************************************************************
subroutine get_lev_box (level, plev, nz, ix, iy, iz, rz, ascending, loglev)

  use messy_clams_global, only: prec

  implicit none

  real(prec),    intent(in)  :: level(:,:,:)
  real(prec),    intent(in)  :: plev
  integer,       intent(in)  :: nz, ix, iy
  integer,       intent(out) :: iz
  real(prec),    intent(out) :: rz
  logical,       intent(in)  :: ascending, loglev

  logical :: found

  integer :: i

  found = .false.
  iz = 1

  if (ascending) then
     do while (iz<nz .and. .not. found)
        if (level(ix,iy,iz)<=plev .and. plev<=level(ix,iy,iz+1)) then
           found = .true.
        else
           iz = iz+1
        endif
     end do
  else
     do while (iz<nz .and. .not. found)
        if (level(ix,iy,iz)>=plev .and. plev>=level(ix,iy,iz+1)) then
           found = .true.
        else
           iz = iz+1
        endif
     end do
  endif

  if (found) then
     if (.not. loglev) then
        rz = (plev-level(ix,iy,iz)) / &
             (level(ix,iy,iz+1)-level(ix,iy,iz))
     else
        rz = (log(plev)-log(level(ix,iy,iz))) / &
             (log(level(ix,iy,iz+1))-log(level(ix,iy,iz)))
     endif
  else
     iz=-1
  endif

end subroutine get_lev_box

!**************************************************************
! Interpolates gridded winds onto particle positions
!**************************************************************
subroutine intruv3d (u,    v,    w,  dlevdz, level, &
                   plons,plats,plevs, &
                   pu,   pv,   pw,  pdlevdz,  & 
                   nx,   ny,   nz,    &
                   longrid, latgrid, levelgrid,    &
                   nparts,loop_start,loop_end,indxps,loglev)

  use messy_clams_global, only: prec, DP, rank, mdi, eps,  &
                                level_is_vertcoor, asc_level, asc_lat, &
                                levelno, &
                                radius_earth

      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer ,    intent(in) :: nx 
      integer ,    intent(in) :: ny 
      integer ,    intent(in) :: nz 
      integer                 :: nparts, loop_start, loop_end
      integer ,    intent(in) :: indxps(nparts) 
      real(prec) , intent(in) :: u(nx,ny,nz) 
      real(prec) , intent(in) :: v(nx,ny,nz) 
      real(prec) , intent(in) :: w(nx,ny,nz) 
      real(prec) , intent(in) :: dlevdz(nx,ny,nz)
      real(prec) , intent(in) :: level(nx,ny,nz)
      real(prec)              :: plons(nparts) 
      real(prec)              :: plats(nparts) 
      real(prec) , intent(in) :: plevs(nparts) 
      real(prec) , intent(out) :: pu(nparts) 
      real(prec) , intent(out) :: pv(nparts) 
      real(prec) , intent(out) :: pw(nparts) 
      real(prec) , intent(out) :: pdlevdz(nparts)
      real(prec) , intent(in) :: longrid(1:nx+1) 
      real(prec) , intent(in) :: latgrid(0:ny+1) 
      real(prec) , intent(in) :: levelgrid(nz) 
      logical, intent(in) :: loglev
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: ix(nparts), iy(nparts), iz(nparts)
      integer :: iz1(nparts), iz2(nparts), iz3(nparts), iz4(nparts)
      real(prec)    :: rx(nparts), ry(nparts), rz(nparts) 
      real(prec)    :: rz1(nparts), rz2(nparts), rz3(nparts), rz4(nparts)
      real(prec)    :: uwk    (nx+1, 0:ny+1, nz+1) 
      real(prec)    :: vwk    (nx+1, 0:ny+1, nz+1) 
      real(prec)    :: wwk    (nx+1, 0:ny+1, nz+1) 
      real(prec)    :: levelwk(nx+1, 0:ny+1, nz+1)
      real(prec)    :: dlevdzwk(nx+1, 0:ny+1, nz+1)
      integer :: i, j, it 
      real(prec) , dimension(0:ny + 1) :: rercos 
      real(prec) , dimension(nx + 1) :: coslon, sinlon 
      real(prec) , dimension(0:ny + 1) :: fmap 
      real(prec) :: pi, r 
      logical :: found 
!-----------------------------------------------

!  Inputs
! u,v,w                 gridded winds
! plons,plats,plevs     positions of particles
! indxps                index of particels on polar stereographic coords
! nparts                number of particles
! nx,ny,nz              size of grid
 
!  Outputs
! pu,pv,pw              Wind values of particles
 
!  Workspace   see Rungek

      pi = 4.*atan(1.) 
  
   ! pad ends of rows and add rows of cross-polar winds
      uwk    (:nx,1:ny,:nz) = u 
      vwk    (:nx,1:ny,:nz) = v 
      wwk    (:nx,1:ny,:nz) = w 
      levelwk(:nx,1:ny,:nz) = level
      dlevdzwk(:nx,1:ny,:nz) = dlevdz
      uwk    (nx+1,1:ny,:nz) = u(1,:,:) 
      vwk    (nx+1,1:ny,:nz) = v(1,:,:) 
      wwk    (nx+1,1:ny,:nz) = w(1,:,:) 
      levelwk(nx+1,1:ny,:nz) = level(1,:,:) 
      dlevdzwk(nx+1,1:ny,:nz)= dlevdz(1,:,:) 

   !  northern boundary
      uwk    (:nx/2,0,:nz) = -u    (1+nx/2:nx,1,:) 
      vwk    (:nx/2,0,:nz) = -v    (1+nx/2:nx,1,:) 
      wwk    (:nx/2,0,:nz) =  w    (1+nx/2:nx,1,:) 
      levelwk(:nx/2,0,:nz) =  level(1+nx/2:nx,1,:) 
      dlevdzwk(:nx/2,0,:nz)=  dlevdz(1+nx/2:nx,1,:) 
 
      uwk    (1+nx/2:nx,0,:nz) = -u    (:nx/2,1,:) 
      vwk    (1+nx/2:nx,0,:nz) = -v    (:nx/2,1,:) 
      wwk    (1+nx/2:nx,0,:nz) =  w    (:nx/2,1,:) 
      levelwk(1+nx/2:nx,0,:nz) =  level(:nx/2,1,:) 
      dlevdzwk(1+nx/2:nx,0,:nz)=  dlevdz(:nx/2,1,:) 

      uwk    (nx+1,0,:nz) = -u    (1,1,:)  
      vwk    (nx+1,0,:nz) = -v    (1,1,:)  
      wwk    (nx+1,0,:nz) =  w    (1,1,:)  
      levelwk(nx+1,0,:nz) =  level(1,1,:)  
      dlevdzwk(nx+1,0,:nz)=  dlevdz(1,1,:)  

    ! southern boundary
      uwk    (:nx/2,ny+1,:nz) = -u    (1+nx/2:nx,ny,:) 
      vwk    (:nx/2,ny+1,:nz) = -v    (1+nx/2:nx,ny,:) 
      wwk    (:nx/2,ny+1,:nz) =  w    (1+nx/2:nx,ny,:) 
      levelwk(:nx/2,ny+1,:nz) =  level(1+nx/2:nx,ny,:) 
      dlevdzwk(:nx/2,ny+1,:nz)=  dlevdz(1+nx/2:nx,ny,:) 

      uwk    (1+nx/2:nx,ny+1,:nz) = -u   (:nx/2,ny,:) 
      vwk    (1+nx/2:nx,ny+1,:nz) = -v   (:nx/2,ny,:) 
      wwk    (1+nx/2:nx,ny+1,:nz) =  w   (:nx/2,ny,:) 
      levelwk(1+nx/2:nx,ny+1,:nz) = level(:nx/2,ny,:) 
      dlevdzwk(1+nx/2:nx,ny+1,:nz)= dlevdz(:nx/2,ny,:) 

      uwk    (nx+1,ny+1,:nz) = -u    (1,ny,:)  
      vwk    (nx+1,ny+1,:nz) = -v    (1,ny,:)  
      wwk    (nx+1,ny+1,:nz) =  w    (1,ny,:)  
      levelwk(nx+1,ny+1,:nz) =  level(1,ny,:)  
      dlevdzwk(nx+1,ny+1,:nz)=  dlevdz(1,ny,:)  

    ! upper boundary
      uwk    (:,:ny+1,nz+1) = uwk    (:,:ny+1,nz) 
      vwk    (:,:ny+1,nz+1) = vwk    (:,:ny+1,nz) 
      wwk    (:,:ny+1,nz+1) = wwk    (:,:ny+1,nz) 
      levelwk(:,:ny+1,nz+1) = levelwk(:,:ny+1,nz) 
      dlevdzwk(:,:ny+1,nz+1)= dlevdzwk(:,:ny+1,nz) 

     call fixll (plons, plats, nparts, loop_start, loop_end, pi) 

   !--------------------------------------------------------------
   !  locate particles to within grid boxes
   !--------------------------------------------------------------
      !DO i = 1,nparts
      DO i = loop_start, loop_end

         if (ABS((plons(i)-mdi)/mdi)<=eps) then 
            rx(i) = mdi 
            ry(i) = mdi 
            rz(i) = mdi 
            ix(i) = -1 
            iy(i) = -1 
            iz(i) = -1 

         else 

            ! rx(i) = (plons(i)-glon0)/dlon 
            ! ix(i) = rx(i) 
            ! rx(i) = rx(i) - ix(i) 
 
            ! ry(i) = (plats(i)-glat0)/dlat 
            ! iy(i) = ry(i) 
            ! ry(i) = ry(i) - iy(i) 

!!!!! 
            if (plons(i)<longrid(1)) then
               ix(i) = nx
               rx(i) = (plons(i)+2*pi-longrid(nx))/(longrid(1)+2*pi-longrid(nx))
            elseif (plons(i)>=longrid(nx)) then
               ix(i) = nx
               rx(i) = (plons(i)-longrid(nx))/(longrid(1)+2*pi-longrid(nx))
            else
               do j = 1, nx-1
                  if (longrid(j)<=plons(i) .and. plons(i)<=longrid(j+1)) then 
                     ix(i) = j
                     rx(i) = (plons(i)-longrid(j))/(longrid(j+1)-longrid(j))
                  endif
               enddo
            endif
            
            if (asc_lat) then  ! latitudes in ascending order
               if (plats(i)<latgrid(1)) then
                  iy(i) = 0
                  if (latgrid(1)==latgrid(0)) then
                     ry(i) = 1.
                  else
                     ry(i) = (plats(i)-latgrid(0))/(latgrid(1)-latgrid(0))
                  endif
               elseif (plats(i)>=latgrid(ny)) then
                  iy(i) = ny
                  if (latgrid(ny+1)==latgrid(ny)) then
                     ry(i) = 0.
                  else
                     ry(i) = (plats(i)-latgrid(ny))/(latgrid(ny+1)-latgrid(ny))
                  endif
               else        
                  do j = 1, ny-1
                     if (latgrid(j)<=plats(i) .and. plats(i)<=latgrid(j+1)) then
                        iy(i) = j
                        ry(i) = (plats(i)-latgrid(j))/(latgrid(j+1)-latgrid(j))
                     endif
                  enddo
               endif
            else ! latitudes in descending order
               if (plats(i)>latgrid(1)) then
                  iy(i) = 0
                  if (latgrid(1)==latgrid(0)) then
                     ry(i) = 1.
                  else
                     ry(i) = (plats(i)-latgrid(0))/(latgrid(1)-latgrid(0))
                  endif
               elseif  (plats(i)<=latgrid(ny)) then
                  iy(i) = ny
                  if (latgrid(ny+1)==latgrid(ny)) then
                     ry(i) = 0.
                  else
                     ry(i) = (plats(i)-latgrid(ny))/(latgrid(ny+1)-latgrid(ny))
                  endif
               else        
                  do j = 1, ny-1
                     if (latgrid(j)>=plats(i) .and. plats(i)>=latgrid(j+1)) then
                        iy(i) = j
                        ry(i) = (plats(i)-latgrid(j))/(latgrid(j+1)-latgrid(j))
                     endif
                  enddo
               endif
            endif
!!$ write (*,*) 'lat=',plats(i)/pi*180.
!!$ write (*,*) 'latgrid:',latgrid(iy(i))/pi*180.,latgrid(iy(i)+1)/pi*180.
!!$ write (*,*) 'iy,ry:',iy(i),ry(i)
!!$ write (*,*) 'lon=',plons(i)/pi*180.
!!$ write (*,*) 'longrid:',longrid(ix(i))/pi*180.,longrid(ix(i)+1)/pi*180.
!!$ write (*,*) 'ix,rx:',ix(i),rx(i)
!!$ write (*,*)


         ! It has to be:  1 <= ix <= nx   and   0 <= iy <= ny
         ! There can be problems with the precision of Real-Values:
         ! If rx(i)=96.9999999, it will be rounded to 97.0 and so
         ! ix(i) is 97 instead of 96.
         !  => ix can be 0 or nx+1,  iy can be ny+1
         !  => This will be checked and corrected.
         ! Attention: iy is allowed to be 0, but perhaps an underflow
         !            can happen ?!
            if ((ix(i) == 0) .and. (ABS(rx(i)-1) <= 1E-4)) then 
               ix(i) = 1 
               rx(i) = 0. 
            else if (ix(i) < 1) then 
               print *,'gone west!   i=',i,'   ix(i)=',ix(i),plons(i) 
               ix(i) = -1
               rx(i) = mdi 
            else if ((ix(i) == nx+1) .and. (rx(i) <= 1E-4)) then
               ix(i) = nx 
               rx(i) = 1. 
            else if (ix(i) > nx) then
               print *,'gone east!   i=',i,'   ix(i)=',ix(i)
               ix(i) = -1 
               rx(i) = mdi 
            endif 
            if ((iy(i) > ny) .or. (iy(i)<0))  then 
               print *,'gone north or south!    i=',i,'   iy(i)=',iy(i),plats(i)
               !iy(i) = ny 
               !ry(i) = 1. 
               rx(i) = mdi 
               ry(i) = mdi 
               rz(i) = mdi 
               ix(i) = -1 
               iy(i) = -1 
               iz(i) = -1 
            endif 

            ! Es muss gelten: 1 <= iy <= ny-1
            ! Fuer iy=0 oder iy=ny wird nur weitergerechnet, wenn man sich
            ! in Polnaehe befindet.
            ! Evtl. Rundungsfehler ????
            if ((iy(i) == 0) .or. (iy(i) == ny)) then
!               if ((pi/2.-plats(i))>ABS(dlat) .and. (plats(i)+pi/2.)>ABS(dlat)) then
               if ((pi/2.-plats(i))>ABS(latgrid(2)-latgrid(1)) .and.  &
                    (plats(i)+pi/2.)>ABS(latgrid(2)-latgrid(1))) then
                  print *,'gone north or south!     i=',i,'   iy(i)=',iy(i), plats(i)
                  rx(i) = mdi 
                  ry(i) = mdi 
                  rz(i) = mdi 
                  ix(i) = -1 
                  iy(i) = -1 
                  iz(i) = -1 
               endif
            endif

            found = .FALSE.

            ! level in clams-file is not vertical coordinate in windfile:
            if (.not. level_is_vertcoor) then

!!!!! iy(i) koennte 0 werden (wenn kein Gitterpunkt am Pol) !?!

               if (ix(i)==nx) then
                  if (iy(i)==0) then
                     iy(i)=1
                     ry(i) =0
                  end if

                  call get_lev_box (levelwk,plevs(i),nz,ix(i),  iy(i),  iz1(i),rz1(i), asc_level, loglev)
                  
                  call get_lev_box (levelwk,plevs(i),nz,1,      iy(i),  iz2(i),rz2(i), asc_level, loglev)
                  
                  call get_lev_box (levelwk,plevs(i),nz,1,      iy(i)+1,iz3(i),rz3(i), asc_level, loglev)
                  
                  call get_lev_box (levelwk,plevs(i),nz,ix(i),  iy(i)+1,iz4(i),rz4(i), asc_level, loglev)

               else
                  if(iy(i)==0) then
                     iy(i) = 1
                     ry(i) = 0
                  end if

                  call get_lev_box (levelwk,plevs(i),nz,ix(i),  iy(i),  iz1(i),rz1(i), asc_level, loglev)
                  
                  call get_lev_box (levelwk,plevs(i),nz,ix(i)+1,iy(i),  iz2(i),rz2(i), asc_level, loglev)
                  
                  call get_lev_box (levelwk,plevs(i),nz,ix(i)+1,iy(i)+1,iz3(i),rz3(i), asc_level, loglev)
                  
                  call get_lev_box (levelwk,plevs(i),nz,ix(i),  iy(i)+1,iz4(i),rz4(i), asc_level, loglev)

               endif
! write (*,*) 'iz:',iz(i),'iz1(i),iz2(i),iz3(i), iz4(i)'
! write (*,*) iz1(i),iz2(i),iz3(i), iz4(i)

               if (iz1(i)==-1 .or. iz2(i)==-1 .or. iz3(i)==-1 .or. iz4(i)==-1) then
                  ix(i) = -1 
                  iy(i) = -1 
                  rx(i) = mdi 
                  ry(i) = mdi 
                  rz(i) = mdi 
               endif
              
            elseif (asc_level) then
               ! level in clams-file is vertical coordinate in windfile
               ! levels ascending (THETA/ZETA)

               ! Liegt das Level im gleichen Bereich wie im vorhergehenden Aufruf ?
               if (levelgrid(levelno(i))<=plevs(i) .and. plevs(i)<=levelgrid(levelno(i)+1)) then 
                  iz(i) = levelno(i) 
                  if (loglev) then ! log. lin.
                     rz(i) = (log(plevs(i))-log(levelgrid(levelno(i))))/ &
                          (log(levelgrid(levelno(i)+1))-log(levelgrid(levelno(i))))
                  else ! lin.
                     rz(i) = (plevs(i)-levelgrid(levelno(i)))/ &
                          (levelgrid(levelno(i)+1)-levelgrid(levelno(i)))
                  endif
                  found = .TRUE.
                  ! Wenn das Level groesser ist als beim vorhergehenden Aufruf:
                  ! Suche von der aktuellen Position nach oben
               elseif (plevs(i)>levelgrid(levelno(i))) then
                  it = levelno(i)+1 
                  do while(it<=nz-1 .and. .not.found) 
                     if (levelgrid(it)<=plevs(i) .and. plevs(i)<=levelgrid(it+1)) then 
                        iz(i) = it 
                        if (loglev) then ! log. lin.
                           rz(i) = (log(plevs(i))-log(levelgrid(it))) / &
                                (log(levelgrid(it+1))-log(levelgrid(it))) 
                        else ! lin. 
                           rz(i) = (plevs(i)-levelgrid(it))/(levelgrid(it+1)-levelgrid(it)) 
                        endif
                        found = .TRUE. 
                        levelno(i) = it
                     endif
                     it = it + 1 
                  end do
                  ! Oberstes Level auf Gleichheit ueberpruefen 
                  if (.not. found) then
                     if (abs((plevs(i)-levelgrid(nz))/plevs(i))<eps) then 
                        iz(i) = nz 
                        rz(i) = 0 
                        found = .TRUE. 
                        levelno(i) = nz-1  !nicht nz wegen der oben stehenden Abfragen!!!
                     endif
                  endif
                  ! Wenn das Level kleiner ist als beim vorhergehenden Aufruf:
                  ! Suche von der aktuellen Position nach unten
               else 
                  it = levelno(i)-1
                  do while(it>=1 .and. .not.found) 
                     if (levelgrid(it)<=plevs(i) .and. plevs(i)<=levelgrid(it+1)) then 
                        iz(i) = it 
                        if (loglev) then ! log. lin
                           rz(i) = (log(plevs(i))-log(levelgrid(it))) / &
                                (log(levelgrid(it+1))-log(levelgrid(it))) 
                        else ! lin.
                           rz(i) = (plevs(i)-levelgrid(it))/(levelgrid(it+1)-levelgrid(it)) 
                        endif
                        found = .TRUE. 
                        levelno(i) = it
                     endif
                     it = it - 1 
                  end do
                  ! Unterstes Level auf Gleichheit ueberpruefen
                  if (.not. found) then
                     if (abs((plevs(i)-levelgrid(1))/plevs(i))<eps) then 
                        iz(i) = 1 
                        rz(i) = 0 
                        found = .TRUE. 
                        levelno(i) = 1  
                     endif
                  endif
               endif
               
               ! if the level is not in the valid range
               if (.not.found) then 
                  rx(i) = mdi 
                  ry(i) = mdi 
                  rz(i) = mdi 
                  ix(i) = -1 
                  iy(i) = -1 
                  iz(i) = -1 
               endif

            else 
               ! level in clams-file is vertical coordinate in windfile
               ! levels descending (PRESS)
               
               ! Liegt das Level im gleichen Bereich wie im vorhergehenden Aufruf ?
               if (levelgrid(levelno(i))>=plevs(i) .and. plevs(i)>=levelgrid(levelno(i)+1)) then 
                  iz(i) = levelno(i) 
                  if (loglev) then ! log. lin.
                     rz(i) = (log(plevs(i))-log(levelgrid(levelno(i))))/ &
                          (log(levelgrid(levelno(i)+1))-log(levelgrid(levelno(i))))
                  else ! lin.
                     rz(i) = (plevs(i)-levelgrid(levelno(i)))/ &
                          (levelgrid(levelno(i)+1)-levelgrid(levelno(i)))
                  endif
                  found = .TRUE.
                  ! Wenn das Level kleiner ist als beim vorhergehenden Aufruf:
                  ! Suche von der aktuellen Position nach oben
               elseif (plevs(i)<levelgrid(levelno(i))) then
                  it = levelno(i)+1 
                  do while(it<=nz-1 .and. .not.found) 
                     if (levelgrid(it)>=plevs(i) .and. plevs(i)>=levelgrid(it+1)) then 
                        iz(i) = it 
                        if (loglev) then ! log. lin.
                           rz(i) = (log(plevs(i))-log(levelgrid(it))) / &
                                (log(levelgrid(it+1))-log(levelgrid(it))) 
                        else ! lin. 
                           rz(i) = (plevs(i)-levelgrid(it))/(levelgrid(it+1)-levelgrid(it)) 
                        endif
                        found = .TRUE. 
                        levelno(i) = it
                     endif
                     it = it + 1 
                  end do
                  ! Oberstes Level auf Gleichheit ueberpruefen 
                  if (.not. found) then
                     if (abs((plevs(i)-levelgrid(nz))/plevs(i))<eps) then 
                        iz(i) = nz 
                        rz(i) = 0 
                        found = .TRUE. 
                        levelno(i) = nz-1  !nicht nz wegen der oben stehenden Abfragen!!!
                     endif
                  endif
                  ! Wenn das Level groesser ist als beim vorhergehenden Aufruf:
                  ! Suche von der aktuellen Position nach unten
               else 
                  it = levelno(i)-1
                  do while(it>=1 .and. .not.found) 
                     if (levelgrid(it)>=plevs(i) .and. plevs(i)>=levelgrid(it+1)) then 
                        iz(i) = it 
                        if (loglev) then ! log. lin
                           rz(i) = (log(plevs(i))-log(levelgrid(it))) / &
                                (log(levelgrid(it+1))-log(levelgrid(it))) 
                        else ! lin.
                           rz(i) = (plevs(i)-levelgrid(it))/(levelgrid(it+1)-levelgrid(it)) 
                        endif
                        found = .TRUE. 
                        levelno(i) = it
                     endif
                     it = it - 1 
                  end do
                  ! Unterstes Level auf Gleichheit ueberpruefen
                  if (.not. found) then
                     if (abs((plevs(i)-levelgrid(1))/plevs(i))<eps) then 
                        iz(i) = 1 
                        rz(i) = 0 
                        found = .TRUE. 
                        levelno(i) = 1  
                     endif
                  endif
               endif
                              
                                 
               ! if the level is not in the valid range
               if (.not.found) then 
                  rx(i) = mdi 
                  ry(i) = mdi 
                  rz(i) = mdi 
                  ix(i) = -1 
                  iy(i) = -1 
                  iz(i) = -1 
               endif

               
            endif
            
         endif

      end do


!!!!! Wichtig, damit Berechnungen korrekt funktionieren !!!
      if (level_is_vertcoor) then
         iz1 = iz
         iz2 = iz
         iz3 = iz
         iz4 = iz
         rz1 = rz
         rz2 = rz
         rz3 = rz
         rz4 = rz
      endif

      ! coordinates to polar-stereo where required
     
      !Do i = 1,nparts
      Do i = loop_start, loop_end
         if (ABS((plats(i)-mdi)/mdi)<=eps .or. ABS((plons(i)-mdi)/mdi)<=eps) cycle  
         if (indxps(i) == 0) cycle  
         r = 2.*radius_earth*cos(plats(i))/(1. + sin(abs(plats(i)))) 
         plats(i) = r*sin(plons(i)) 
         plons(i) = r*cos(plons(i)) 
      end do
      
   ! interpolate and change coordinates

!!!!!
      ! do j = 0, ny + 1 
      !    rercos(j) = 1./(radius_earth*cos(glat0 + j*dlat)) 
      ! end do 
      do j = 0, ny + 1 
         rercos(j) = 1./(radius_earth*cos(latgrid(j))) 
      end do 
 
!!!!!
      ! do i = 1, nx + 1 
      !    coslon(i) = cos(glon0 + i*dlon) 
      !    sinlon(i) = sin(glon0 + i*dlon) 
      ! end do 
      do i = 1, nx + 1  
         coslon(i) = cos(longrid(i)) 
         sinlon(i) = sin(longrid(i)) 
      end do 
 
!!!!!
      ! do j = 1, ny 
      !    fmap(j) = 2./(1. + sin(abs(glat0 + j*dlat))) 
      ! end do 
      do j = 1, ny 
         fmap(j) = 2./(1. + sin(abs(latgrid(j)))) 
      end do 
      fmap(0) = fmap(1) 
      fmap(ny+1) = fmap(ny) 

      !DO i = 1,nparts
      DO i = loop_start, loop_end

         IF (ix(i) == -1) THEN
            pu(i) = MDI
            pv(i) = MDI
            pw(i) = MDI
            pdlevdz(i) = MDI

         ELSEIF (ABS((uwk(ix(i),  iy(i)  ,iz1(i))-MDI)/MDI)<=eps .OR. &
                 ABS((uwk(ix(i)  ,iy(i)  ,iz1(i)+1)-MDI)/MDI)<=eps .OR. &
                 ABS((uwk(ix(i)+1,iy(i)  ,iz2(i))-MDI)/MDI)<=eps .OR. &
                 ABS((uwk(ix(i)+1,iy(i)  ,iz2(i)+1)-MDI)/MDI)<=eps .OR. &
                 ABS((uwk(ix(i)+1,iy(i)+1,iz3(i))-MDI)/MDI)<=eps .OR. &
                 ABS((uwk(ix(i)+1,iy(i)+1,iz3(i)+1)-MDI)/MDI)<=eps .OR. &
                 ABS((uwk(ix(i)  ,iy(i)+1,iz4(i))-MDI)/MDI)<=eps .OR. &
                 ABS((uwk(ix(i)  ,iy(i)+1,iz4(i)+1)-MDI)/MDI)<=eps) THEN
              
            pu(i) = MDI
            pv(i) = MDI
            pw(i) = MDI
            pdlevdz(i) = MDI

         ELSE

            IF (indxps(i)/=0) THEN         

          pu(i) =   &
             ry(i) * (rx(i) *  (rz3(i) * (-vwk(ix(i)+1,iy(i)+1,iz3(i)+1)) &
                                         *coslon(ix(i)+1)*indxps(i)     &     
                              + rz3(i) * (-uwk(ix(i)+1,iy(i)+1,iz3(i)+1)) &
                                         *sinlon(ix(i)+1)               &
                          + (1.-rz3(i))* (-vwk(ix(i)+1,iy(i)+1,iz3(i)))   &
                                         *coslon(ix(i)+1)*indxps(i)     &     
                          + (1.-rz3(i))* (-uwk(ix(i)+1,iy(i)+1,iz3(i)))   &
                                         *sinlon(ix(i)+1)             ) &     
                 + (1.-rx(i))* (rz4(i) * (-vwk(ix(i),  iy(i)+1,iz4(i)+1)) &
                                         *coslon(ix(i))  *indxps(i)     &   
                              + rz4(i) * (-uwk(ix(i),  iy(i)+1,iz4(i)+1)) &
                                         *sinlon(ix(i))                 &   
                          + (1.-rz4(i))* (-vwk(ix(i),  iy(i)+1,iz4(i)))   &
                                         *coslon(ix(i))  *indxps(i)     &
                          + (1.-rz4(i))* (-uwk(ix(i),  iy(i)+1,iz4(i)))   &
                                         *sinlon(ix(i))              )) &
                 * fmap(iy(i)+1)                                        &
       + (1.-ry(i)) * (rx(i) * (rz2(i) * (-vwk(ix(i)+1,iy(i),  iz2(i)+1)) &
                                         *coslon(ix(i)+1)*indxps(i)     &
                              + rz2(i) * (-uwk(ix(i)+1,iy(i),  iz2(i)+1)) &
                                         *sinlon(ix(i)+1)               &
                          + (1.-rz2(i))* (-vwk(ix(i)+1,iy(i),  iz2(i)))   &
                                         *coslon(ix(i)+1)*indxps(i)     &     
                          + (1.-rz2(i))* (-uwk(ix(i)+1,iy(i),  iz2(i)))   &
                                         *sinlon(ix(i)+1)             ) &     
                 + (1.-rx(i))* (rz1(i) * (-vwk(ix(i),  iy(i),  iz1(i)+1)) &
                                         *coslon(ix(i))  *indxps(i)     & 
                              + rz1(i) * (-uwk(ix(i),  iy(i),  iz1(i)+1)) &
                                         *sinlon(ix(i))                 &
                          + (1.-rz1(i))* (-vwk(ix(i),  iy(i),  iz1(i)))   &
                                         *coslon(ix(i))  *indxps(i)     &
                          + (1.-rz1(i))* (-uwk(ix(i),  iy(i),  iz1(i)))   &
                                         *sinlon(ix(i))         ))      &     
                  *fmap(iy(i))


        pv(i) =  &
             ry(i) * (rx(i) *  (rz3(i) * (-vwk(ix(i)+1,iy(i)+1,iz3(i)+1)) &
                                         *sinlon(ix(i)+1)*indxps(i)     &     
                              + rz3(i) *  uwk(ix(i)+1,iy(i)+1,iz3(i)+1) &
                                         *coslon(ix(i)+1)               &
                          + (1.-rz3(i))* (-vwk(ix(i)+1,iy(i)+1,iz3(i)))   &
                                         *sinlon(ix(i)+1)*indxps(i)     &     
                          + (1.-rz3(i))*  uwk(ix(i)+1,iy(i)+1,iz3(i))   &
                                         *coslon(ix(i)+1)             ) &     
                 + (1.-rx(i))* (rz4(i) * (-vwk(ix(i),  iy(i)+1,iz4(i)+1)) &
                                         *sinlon(ix(i))  *indxps(i)     &   
                              + rz4(i) *  uwk(ix(i),  iy(i)+1,iz4(i)+1) &
                                         *coslon(ix(i))                 &   
                          + (1.-rz4(i))* (-vwk(ix(i),  iy(i)+1,iz4(i)))   &
                                         *sinlon(ix(i))  *indxps(i)     &
                          + (1.-rz4(i))*  uwk(ix(i),  iy(i)+1,iz4(i))   &
                                         *coslon(ix(i))              )) &
                 * fmap(iy(i)+1)                                        &
       + (1.-ry(i)) * (rx(i) * (rz2(i) * (-vwk(ix(i)+1,iy(i),  iz2(i)+1)) &
                                         *sinlon(ix(i)+1)*indxps(i)     &
                              + rz2(i) * uwk(ix(i)+1,iy(i),  iz2(i)+1)  &
                                         *coslon(ix(i)+1)               &
                          + (1.-rz2(i))* (-vwk(ix(i)+1,iy(i),  iz2(i)))   &
                                         *sinlon(ix(i)+1)*indxps(i)     &     
                          + (1.-rz2(i))* uwk(ix(i)+1,iy(i),  iz2(i))    &
                                         *coslon(ix(i)+1)             ) &     
                + (1.-rx(i)) * (rz1(i) * (-vwk(ix(i),  iy(i),  iz1(i)+1)) &
                                         *sinlon(ix(i))  *indxps(i)     & 
                              + rz1(i) * uwk(ix(i),  iy(i),  iz1(i)+1)  &
                                         *coslon(ix(i))                 &
                          + (1.-rz1(i))* (-vwk(ix(i),  iy(i),  iz1(i)))   &
                                         *sinlon(ix(i))  *indxps(i)     &
                          + (1.-rz1(i))* uwk(ix(i),  iy(i),  iz1(i))    &
                                         *coslon(ix(i))         ))      &     
                    *fmap(iy(i))

                
        ELSE                                                                 

        pu(i) =  &
                ry(i) *  (rx(i)  * (rz3(i) * uwk(ix(i)+1,iy(i)+1,iz3(i)+1)  &
                              + (1.-rz3(i))* uwk(ix(i)+1,iy(i)+1,iz3(i)))   &
                    + (1.-rx(i)) * (rz4(i) * uwk(ix(i),  iy(i)+1,iz4(i)+1)  &
                              + (1.-rz4(i))* uwk(ix(i),  iy(i)+1,iz4(i))))  &
                     *rercos(iy(i)+1)                                      &      
         + (1.-ry(i)) *  (rx(i) * (rz2(i) * uwk(ix(i)+1,iy(i),  iz2(i)+1)   &
                             + (1.-rz2(i))* uwk(ix(i)+1,iy(i),  iz2(i)))    &
                    + (1.-rx(i))* (rz1(i) * uwk(ix(i),  iy(i),  iz1(i)+1)   &
                             + (1.-rz1(i))* uwk(ix(i),  iy(i),  iz1(i))))   &
                     *rercos(iy(i))                
       


       pv(i) = (1./radius_earth)*(                                                 &
             ry(i)  * (rx(i)  * (rz3(i) * vwk(ix(i)+1,iy(i)+1,iz3(i)+1)  &
                           + (1.-rz3(i))* vwk(ix(i)+1,iy(i)+1,iz3(i)))   &
                 + (1.-rx(i)) * (rz4(i) * vwk(ix(i),  iy(i)+1,iz4(i)+1)  &
                           + (1.-rz4(i))* vwk(ix(i),  iy(i)+1,iz4(i))))  &
       + (1.-ry(i)) * (rx(i)  * (rz2(i) * vwk(ix(i)+1,iy(i),  iz2(i)+1)  &
                           + (1.-rz2(i))* vwk(ix(i)+1,iy(i),  iz2(i)))   &
                 + (1.-rx(i)) * (rz1(i) * vwk(ix(i),  iy(i),  iz1(i)+1)  &
                           + (1.-rz1(i))* vwk(ix(i),  iy(i),  iz1(i)))))



        ENDIF
   
        pw(i) =                                                    &
            ry(i) * (                                              &
                  rx(i)  * (rz3(i) *wwk(ix(i)+1,iy(i)+1,iz3(i)+1)  &
                      + (1.-rz3(i))*wwk(ix(i)+1,iy(i)+1,iz3(i)))   &
             +(1.-rx(i)) * (rz4(i) *wwk(ix(i)  ,iy(i)+1,iz4(i)+1)  &
                      + (1.-rz4(i))*wwk(ix(i)  ,iy(i)+1,iz4(i))))  &
         + (1.-ry(i))*(                                            &
                  rx(i)  * (rz2(i) *wwk(ix(i)+1,iy(i)  ,iz2(i)+1)  &
                      + (1.-rz2(i))*wwk(ix(i)+1,iy(i)  ,iz2(i)))   &
             +(1.-rx(i)) * (rz1(i) *wwk(ix(i)  ,iy(i),  iz1(i)+1)  &
                      + (1.-rz1(i))*wwk(ix(i)  ,iy(i),  iz1(i))))
     

        pdlevdz(i) =                                         &
             ry(i) * (                                       &
                rx(i)  * (rz3(i) *dlevdzwk(ix(i)+1,iy(i)+1,iz3(i)+1)  &
                    + (1.-rz3(i))*dlevdzwk(ix(i)+1,iy(i)+1,iz3(i)))   &
           +(1.-rx(i)) * (rz4(i) *dlevdzwk(ix(i)  ,iy(i)+1,iz4(i)+1)  &
                    + (1.-rz4(i))*dlevdzwk(ix(i)  ,iy(i)+1,iz4(i))))  &
       + (1.-ry(i))*(                                     &
                rx(i)  * (rz2(i) *dlevdzwk(ix(i)+1,iy(i)  ,iz2(i)+1)  &
                    + (1.-rz2(i))*dlevdzwk(ix(i)+1,iy(i)  ,iz2(i)))   &
           +(1.-rx(i)) * (rz1(i) *dlevdzwk(ix(i)  ,iy(i),  iz1(i)+1)  &
                    + (1.-rz1(i))*dlevdzwk(ix(i)  ,iy(i),  iz1(i))))
                
     ENDIF

     ENDDO    

END SUBROUTINE intruv3d

!********************************************************************
! converts particle coordinates from polar-stereographic to latitudes    
! and longitudes                                                         
!********************************************************************
SUBROUTINE ptoll (plons,plats,indxps,nparts,loop_start,loop_end,pi)                    
          
  use messy_clams_global, only: prec, mdi, eps,  &
                                radius_earth

   IMPLICIT NONE
                                                                     
   INTEGER    :: nparts, loop_start, loop_end, i
   REAL(PREC) :: plons(nparts),plats(nparts)
   INTEGER    :: indxps(nparts)        
   REAL(PREC) :: pi, r

   !DO i = 1,nparts
   DO i = loop_start, loop_end
      !IF ((plons(i) /= MDI) .AND. (plats(i) /= MDI)) THEN
      IF (ABS((plons(i)-MDI)/MDI)>eps .AND. ABS((plats(i)-MDI)/MDI)>eps) THEN
         IF (indxps(i)/=0) THEN                                             
            r = sqrt(plons(i)*plons(i)+plats(i)*plats(i))                     
            plons(i) = atan2(plats(i),plons(i))                               
            plats(i) = indxps(i)*(.5*pi - 2.*atan(.5*r/radius_earth))                   
         ENDIF
      ENDIF
   ENDDO 

END SUBROUTINE ptoll

!********************************************************************
!   Adjusts longitudes and latitudes in plons and plats to ensure that      
!   all values are in range  0,2*pi) and  -.5*pi,.5*pi . It is assumed      
!   that the values passed into this routine are in the range               
!    -2.*pi,3.*pi) for longitude and  -1.5*pi,1.5*pi  for latitude.         
!********************************************************************
!SUBROUTINE fixll (plons,plats,MDI,eps,nparts,pi,glon0,dlon)                   
SUBROUTINE fixll (plons,plats,nparts,loop_start,loop_end,pi)                   
   
  use messy_clams_global, only: prec, mdi, eps 

   IMPLICIT NONE

   INTEGER    :: nparts, loop_start, loop_end, i                               
   REAL(PREC) :: plons(nparts),plats(nparts)  
   REAL(PREC) :: pi

   !DO i = 1,nparts
   DO i = loop_start, loop_end

      !IF ((plons(i) /= MDI) .AND. (plats(i) /= MDI)) THEN
      IF (ABS((plons(i)-MDI)/MDI)>eps .AND. ABS((plats(i)-MDI)/MDI)>eps) THEN
         if (plats(i)>.5*pi) then                                          
            print *,'north before',plats(i)*180./pi                           
            plats(i) = pi - plats(i)                                          
            plons(i) = plons(i) + pi                                          
            print *,'north after',plats(i)*180./pi                            
         elseif (plats(i)<-.5*pi) then                                    
            print *,'south before',plats(i)*180./pi                           
            plats(i) = -pi - plats(i)                                         
            plons(i) = plons(i) + pi                                          
            print *,'south after',plats(i)*180./pi                            
         endif
      ENDIF
   ENDDO
  
   !DO i = 1,nparts
   DO i = loop_start, loop_end
      !IF ((plons(i) /= MDI) .AND. (plats(i) /= MDI)) THEN
      IF (ABS((plons(i)-MDI)/MDI)>eps .AND. ABS((plats(i)-MDI)/MDI)>eps) THEN
         ! if (plons(i)>=2.*pi+glon0+dlon) then                               
         !    plons(i) = plons(i) - 2.*pi                                       
         !    if (plons(i)<glon0+dlon) plons(i)=glon0+dlon
         ! else if (plons(i)<glon0+dlon) then    
         !    plons(i) = plons(i) + 2.*pi      
         !    if (plons(i)>=2.*pi+glon0+dlon) plons(i)=1.9999*pi+glon0+dlon 
         ! endif
         if (plons(i)>=2.*pi) then                               
            plons(i) = plons(i) - 2.*pi                                       
            if (plons(i)<0.) plons(i)=0.
         else if (plons(i)<0.) then    
            plons(i) = plons(i) + 2.*pi      
            if (plons(i)>=2.*pi) plons(i)=1.9999*pi 
         endif
      ENDIF
   ENDDO

END SUBROUTINE fixll


!**************************************************************
! calculates the cosine of the solar zenith angle        
!**************************************************************
SUBROUTINE zen2(zlat,zlong,idum,time,cozen,SZA)               

use messy_clams_global, only: prec

implicit none

  REAL(PREC) :: zlat,zlong,cozen,time,H,phi,PI
  REAL(PREC) :: dec,SZA
  INTEGER    :: IDUM(3)

  ! write(6,*)  'ZEN2 Entered'

  PI=4.0*ATAN(1.0)

  ! Calculate the hour angle 0 for local noon
  ! and 1 hour for every 15 deg.
  ! time is supplied as a fractional day 24hrs = 1.0 

  H = 2*PI*(ZLONG/360.0+time-0.5)
  Phi = (PI/180.0)*zlat
  dec = DECLINE(IDUM(1),IDUM(2),IDUM(3))
  cozen = cos(phi)*cos(dec)*cos(H)+sin(phi)*sin(dec)
  if (cozen<-1 .or. cozen>1.)then
     write (*,*) 'arcus cosinus ausserhalb des Wertebereichs !'
  endif
  !if (cozen<-1. .and. abs(cozen+1.)<eps) then
  !   SZA = 180.
  !elseif (cozen>1 .and. abs(cozen-1.)<eps) then
  !   SZA = 0.
  !else
  !   SZA  = ACOS(cozen)*(180.0/PI)
  !endif
  ! Wenn cozen nicht im Defintionsbereich [-1,1] => Arithmetical exception!
  SZA  = ACOS(cozen)*(180.0/PI)
  
  ! write(6,*)  'ZEN2 Exited'
  
  contains

  REAL(PREC) FUNCTION DECLINE(iy,Im,ID)   

  REAL(PREC) :: DMAX,PI   
  INTEGER    :: months(12),YDAYS,DAYS,sdays,iy,im,id

  integer :: iyr, ires, jm

  DATA months /31,28,31,30,31,30,31,31,30,31,30,31/
  DATA YDAYS,DAYS /365,0/
  DATA DMAX /23.5/
    
  PI= 4.0*ATAN(1.0)

  !    Function to calculate the declination of the
  !    sun in the sky being:-
  !
  !      0   at the Spring Equinox      March 21
  !    +23.5 at the Summer Solstice      June 21.
  !      0   at the Autumnal Equinox Septeber 21
  !    -23.5 at the Winter Solstice  December 21
  
  ! 1. Determine if it is a leap year. i.e. the 
  !    year is divisable by 4

  iyr = (iy/4)*4
  ires = iy-iyr
  IF (ires == 0) THEN
     ydays = 366
     months(2) = 29
  ELSE
     ydays = 365
     months(2) = 28
  ENDIF
   
  ! 2.  Calculate the fractional time of Year
       
  Days = 0
  DO jm = 1, im-1
     Days = Days + Months(jm)
  ENDDO
  Days = Days + id

  ! 3.   Now calculate offset for the spring
  !        Equinox  

  SDAYS = months(1)+months(2)+21

  ! 4.   Now calculate the Declination
  
  DECLINE= (PI/180.0)*DMAX*sin(2*PI*(DAYS-SDAYS)/YDAYS)

  END FUNCTION decline

END SUBROUTINE zen2

END MODULE messy_clamstraj_calc3d

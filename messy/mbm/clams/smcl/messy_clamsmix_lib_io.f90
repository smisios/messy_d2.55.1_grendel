!********************************************************************************!
! File clams/mix/source/lib_io.f90
!
! This file is part of the Chemical Lagrangian Model of the
! Stratosphere (CLaMS). CLaMS is a hierarchy of numerical models and
! preprocessors which simulate transport, chemical reactions and
! mixing processes in the stratosphere.  
!
! Copyright (c) 2006
! Paul Konopka, Nicole Thomas
! Forschungszentrum Juelich GmbH
! Last Modified By: N.Thomas
! Last Modified On: Mon May 18 14:11:39 2020
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
!***********************************************************************!
!                                                                   
!   Module lib_io                                                 
! ----------------
!                                                                   
!   Module for input/output operations in mix                                
!                                                                   
!   Public subroutines:    
!                                       
!      subroutine read_init_mix_config
!      subroutine get_ap_s
!      logical function find_index
!      subroutine save_ap_s 
!
!      For further informations, see inline commentars.             
!                                                                   
!----------------------------------------------------------------------

Module messy_clamsmix_lib_io

contains


subroutine read_init_mix_config (status, file_init, &
                               lev_min, lev_max, nlevs, &
                               lat_down, lat_up, lat_min, lat_max, &
                               r_coarse, r_high, &
                               lev_grid,lev_delta,r_grid, &
                               vertcoorname)

  use messy_clams_global,        only: prec, buffersize
  use messy_clamsmix_global,     only: ctrl_out
  use messy_clams_tools_utils,   only: uppercase, lowercase
  use messy_clams_tools_ncutils, only: nc_check_error
  use netcdf

  implicit none 

  character*(*), intent(in)                        :: file_init
  real(prec), intent(inout)                        :: lev_min, lev_max
  integer, intent(inout)                           :: status, nlevs
  real(prec), intent(inout)                        :: lat_down, lat_up, &
                                                      lat_min, lat_max, &
                                                      r_coarse, r_high
  real(prec), dimension(:),pointer                 :: lev_grid 
  real(prec), dimension(:),pointer                 :: lev_delta
  real(prec), dimension(:),pointer                 :: r_grid
  character(*), intent(in)                         :: vertcoorname

  character(80)                                    :: helpstr

  integer                                          :: rcode, ncid, dummy_id

  integer                                          :: varid
  integer                                          :: nlevs_value

  integer, dimension(1)                            :: data_int
  real(prec), dimension(1)                         :: data_real

  status = 0 ! no error

  if (ctrl_out) print *, 'Netcdf-info file: ', file_init

  status = nf90_open(file_init,nf90_nowrite,ncid, buffersize)
  call nc_check_error (status,'Error on open file '//trim(file_init),abort=.false.)
  if (status/=0) return
    
  ! Read global parameters describing initial distribution
  helpstr = 'exp_POS_'//TRIM(lowercase(vertcoorname))//'_min' ! exp_POS_theta/zeta_min
  status = nf90_get_att (ncid,NF90_GLOBAL,TRIM(helpstr), data_real)
  call nc_check_error (status,'Cannot get attribute '//trim(helpstr),abort=.false.)
  if (status/=0) return
  lev_min = data_real(1)
  helpstr = 'exp_POS_'//TRIM(lowercase(vertcoorname))//'_max' ! exp_POS_theta/zeta_max
  status = nf90_get_att (ncid,NF90_GLOBAL,TRIM(helpstr), data_real)
  call nc_check_error (status,'Cannot get attribute '//trim(helpstr),abort=.false.)
  if (status/=0) return
  lev_max = data_real(1)
  helpstr = 'exp_POS_n'//TRIM(lowercase(vertcoorname))//'s'   ! exp_POS_nthetas/nzetas
  status = nf90_get_att (ncid,NF90_GLOBAL,TRIM(helpstr), data_int)
  call nc_check_error (status,'Cannot get attribute '//trim(helpstr),abort=.false.)
  if (status/=0) return
  nlevs = data_int(1)
  status = nf90_get_att (ncid,NF90_GLOBAL,'exp_POS_lat_down', data_real)
  call nc_check_error (status,'Cannot get attribute exp_POS_lat_down',abort=.false.)
  if (status/=0) return
  lat_down = data_real(1)
  status = nf90_get_att (ncid,NF90_GLOBAL,'exp_POS_lat_up', data_real)
  call nc_check_error (status,'Cannot get attribute exp_POS_lat_up',abort=.false.)
  if (status/=0) return
  lat_up = data_real(1)
  status = nf90_get_att (ncid,NF90_GLOBAL,'exp_POS_lat_min', data_real)
  call nc_check_error (status,'Cannot get attribute exp_POS_lat_min',abort=.false.)
  if (status/=0) return
  lat_min = data_real(1)
  status = nf90_get_att (ncid,NF90_GLOBAL,'exp_POS_lat_max', data_real)
  call nc_check_error (status,'Cannot get attribute exp_POS_lat_max',abort=.false.)
  if (status/=0) return
  lat_max = data_real(1)
  status = nf90_get_att (ncid,NF90_GLOBAL,'exp_POS_r_coarse', data_real)
  call nc_check_error (status,'Cannot get attribute exp_POS_r_coarse',abort=.false.)
  if (status/=0) return
  r_coarse = data_real(1)
  status = nf90_get_att (ncid,NF90_GLOBAL,'exp_POS_r_high', data_real)
  call nc_check_error (status,'Cannot get attribute exp_POS_r_high',abort=.false.)
  if (status/=0) return
  r_high = data_real(1)

  ! Read nlevs (# of lev niveaus)
  helpstr = 'N'//TRIM(uppercase(vertcoorname))//'S'  ! NTHETAS/NZETAS
  status = nf90_inq_dimid(ncid,TRIM(helpstr),dummy_id)     
  call nc_check_error (status, &
       'Cannot find dimension N'//TRIM(uppercase(vertcoorname))//'S',abort=.false.)
  if (status/=0) return

  ! lev information is in the file
  status = nf90_inquire_dimension(ncid,dummy_id,len=nlevs_value)
  call nc_check_error (status, &
       'Cannot read dimension N'//TRIM(uppercase(vertcoorname))//'S',abort=.false.)
  if (status/=0) return
  allocate(lev_grid(nlevs_value),lev_delta(nlevs_value))
  
  helpstr = TRIM(uppercase(vertcoorname))//'_GRID'  ! THETA_GRID/ZETA_GRID
  status = nf90_inq_varid (ncid,TRIM(helpstr),varid)     
  call nc_check_error (status,'Cannot find variable '//trim(helpstr),abort=.false.)
  if (status/=0) return
  status = nf90_get_var (ncid,varid,lev_grid)
  call nc_check_error (status,'Cannot read variable '//trim(helpstr),abort=.false.)
  if (status/=0) return
  helpstr = TRIM(uppercase(vertcoorname))//'_DELTA' ! THETA_DELTA/ZETA_DELTA
  status = nf90_inq_varid (ncid,TRIM(helpstr),varid)
  call nc_check_error (status,'Cannot find variable '//trim(helpstr),abort=.false.)
  if (status/=0) return
  status = nf90_get_var (ncid,varid,lev_delta)
  call nc_check_error (status,'Cannot read variable '//trim(helpstr),abort=.false.)
  if (status/=0) return
  
  rcode = nf90_inq_varid(ncid,'R_GRID',varid)
  if (rcode==0) then
     allocate (r_grid(nlevs_value))
     rcode = nf90_get_var (ncid,varid,r_grid)
     call nc_check_error (rcode,'Cannot read variable R_GRID',abort=.false.)
  endif
  
  rcode=nf90_close(ncid)

end subroutine read_init_mix_config

!###########################################################################

subroutine get_ap_s(status, switch_mixing, ap_s, layer_ind, lat, lon, zeta,&
                    lat_old, lon_old, specarr, &
                    theta, bvf, theta_old, bvf_old)

  use messy_clams_global,         only: prec, dp, species_type
  use messy_clamsmix_global,      only: ap, timestep_mix, ctrl_out
  use messy_clams_tools_dateconv, only: ymds2js
  use messy_clams_global,         only: YEAR, MONTH, DAY, HOUR, &
                                        MINUTE, SECOND, delta_time

  implicit none 

  integer                                         :: status, switch_mixing
  type(ap),           dimension(:), intent(inout) :: ap_s
  integer,            dimension(:)                :: layer_ind
  real(prec),         dimension(:)                :: lat, lon, zeta 
  real(prec),         dimension(:)                :: lat_old, lon_old 
  TYPE(species_type), DIMENSION(:), POINTER       :: specarr
  real(prec),         dimension(:)                :: theta, bvf
  real(prec),         dimension(:)                :: theta_old, bvf_old

  real(prec)    :: pi
  integer       :: nsteps,i,j,size_ap_s, dim_len
  integer       :: sod, julsec
  integer       :: ipos


  status = 0 ! no error

  ! Define pi
  pi=4.*atan(1.)
 
  ! Define the size of the APs structure
  size_ap_s=size(ap_s)

  ! Nullify species-pointer
  do i=1,size_ap_s
     nullify(ap_s(i)%c)
  end do

  ! Read tracer (currently not available)
!!!!!
  !print *, 'Tracer is not available (set on 0.0)'
  ap_s(:)%tracer=0.0

  nsteps = timestep_mix*3600/delta_time +1  

  sod = HOUR*3600 +MINUTE*60 + SECOND
  julsec = ymds2js(YEAR, MONTH, DAY, sod) 
  ap_s(:)%time_init=julsec + delta_time

  ! Latitude and longitude before the MIX timestep
  ap_s(:)%lat=lat_old(layer_ind(:))
  ap_s(:)%lon=lon_old(layer_ind(:))

  ! Transform these positions to (old) spherical coordinates on a unit sphere
  ap_s%coor_old(1) = cos(ap_s%lon*(pi/180.))*sin(abs(ap_s%lat*(pi/180.)-pi/2.)) 
  ap_s%coor_old(2) = sin(ap_s%lon*(pi/180.))*sin(abs(ap_s%lat*(pi/180.)-pi/2.))
  ap_s%coor_old(3) = cos(abs(ap_s%lat*(pi/180.)-pi/2.))

  ! Read lon, lat positions at the end of the time step
  ap_s(:)%lon=lon(layer_ind(:))
  ap_s(:)%lat=lat(layer_ind(:))

  ! Transform these positions to (new) spherical coordinates on a unit sphere
  ap_s%coor(1) = cos(ap_s%lon*(pi/180.))*sin(abs(ap_s%lat*(pi/180.)-pi/2.)) 
  ap_s%coor(2) = sin(ap_s%lon*(pi/180.))*sin(abs(ap_s%lat*(pi/180.)-pi/2.))
  ap_s%coor(3) = cos(abs(ap_s%lat*(pi/180.)-pi/2.))
  
  ! Read Lev at the end of the advection step
  ap_s(:)%lev = zeta(layer_ind(:))


  ! if vert_mixing: set theta and bvf
   if (switch_mixing == 2) then

     ap_s(:)%theta      = lat(layer_ind(:)) 
     ap_s(:)%bvf        = bvf(layer_ind(:))
     ap_s(:)%theta_old  = theta_old(layer_ind(:))
     ap_s(:)%bvf_old    = bvf_old(layer_ind(:))

  else

     ap_s(:)%theta      = 0.
     ap_s(:)%theta_old  = 0.
     ap_s(:)%bvf        = 0.
     ap_s(:)%bvf_old    = 0.

  endif


  ! Read chemical species (if present)
  if (associated(specarr)) then

     dim_len=size(specarr)
     do i = 1, size_ap_s
        allocate(ap_s(i)%c(dim_len))
     enddo

     do j = 1, size_ap_s
        do i = 1, dim_len
           ap_s(j)%c(i) = specarr(i)%values(layer_ind(j))
        enddo
     enddo
     
    if (ctrl_out) print *, 'Tracers have been read'       
    
   endif       ! End read of chem. species

end subroutine get_ap_s

!###########################################################################

logical function find_index (status, lev_min, lev_max, ind_set,  &
                             lat, lev, lat_min, lat_max) result(found)

  use messy_clams_global,    only: prec
  use messy_clamsmix_global, only: ctrl_out

  implicit none 

  integer                                    :: status
  real(prec), intent(in)                     :: lev_min, lev_max
  integer, dimension(:), pointer             :: ind_set
  real(prec), dimension(:), intent(in)       :: lat
  real(prec), dimension(:), intent(in)       :: lev
  real(prec), intent(in), optional           :: lat_min, lat_max

  integer, parameter                         :: layer_limit=1
  integer                                    :: nparts_max
  integer, allocatable, dimension(:)         :: ind
  integer                                    :: i, counter
  real(prec)                                 :: latmin, latmax

  status = 0 ! no error

  nparts_max = size(lat)

! Set default value for the result
  found=.false.

! Check that lev_min and lev_max define a layer
  if (abs(lev_max-lev_min) < 0.2) then 
    print *, '!!!! The layer is too thin !!!!  '
    status = -1
    return
  endif

  allocate(ind(nparts_max))

  if (present(lat_min)) then
     latmin=lat_min
  else
     latmin=-90.
  endif
  if (present(lat_max)) then
     latmax=lat_max
  else
     latmax=90.
  endif

  ! Determine the index set of APs with :
  ! lev_min <= lev < lev_max and  latmin <= lat <= latmax 
  counter=0
  ind=0
  do i=1, nparts_max
    if ((lev(i) >= lev_min) .and. (lev(i) <  lev_max) .and. &
        (lat(i)   >= latmin   ) .and. (lat(i)   <= latmax   )) then 
        counter=counter+1
        ind(counter)=i
    endif 
  enddo
  if (ctrl_out) print *, 'Number of APs in considered levrange :', counter

  ! Create ind_set and produce results
  if (counter < layer_limit) then
    print *, '# of valid APs between', lev_min, lev_max ,'smaller than', layer_limit 
  else 
    found=.true.
    allocate(ind_set(counter))
    ind_set=ind(1:counter)
  endif 

  ! Deallocate all arrays
  deallocate(ind) 
 
end function find_index

!*****************************************************************************

  ! subroutine to save the results 
  subroutine save_ap_s (status, ap_s, lats_out, lons_out, levs_out, &
                        state_out, state_vert_out, specarr_out,  &
                        writepos, count, &
                        lev_min, lev_max, &
                        lat_down, lat_up,  &
                        lev_grid, lev_delta, r_grid, &
                        check_subset)

    use messy_clams_global, only: species_type, prec
    use messy_clamsmix_ap_m_access, only: ap_struct, subset, n, &
                                          lev, lat, lon, state, state_vert, c

    implicit none 

    integer                                :: status
    type(ap_struct)                        :: ap_s
    TYPE(species_type), DIMENSION(:), POINTER :: specarr_out
    real(prec)                             :: lats_out(:)
    real(prec)                             :: lons_out(:)
    real(prec)                             :: levs_out(:)
    real(prec)                             :: state_out(:)
    real(prec)                             :: state_vert_out(:)
    integer, intent(in)                    :: writepos
    integer, intent(out)                   :: count
    real(prec), intent(in)                 :: lev_min, lev_max, &
                                              lat_down, lat_up
    logical, optional                      :: check_subset
   
    real(prec), dimension(:),pointer       :: lev_grid 
    real(prec), dimension(:),pointer       :: lev_delta
    real(prec), dimension(:),pointer       :: r_grid


    integer                                :: i,j,k

    integer, dimension(1)                  :: start
    ! real(dp), dimension(:), &
    !      allocatable                       :: time_init_out
    integer                                :: count_recv

    logical                                :: ctrl_subset
 
    status = 0 ! no error

    ctrl_subset = .false.
    if (present(check_subset)) ctrl_subset = check_subset

    count=0

    do i=1, n

       if (ctrl_subset) then

          if ((lev(ap_s,i) >= lev_min) .and. (lev(ap_s,i) <= lev_max) .and. &
               (lat(ap_s,i) >= lat_down) .and. (lat(ap_s,i) <= lat_up)) then 

             if (subset(ap_s,i)) then
                lats_out   (writepos+count) = lat(ap_s,i)
                lons_out   (writepos+count) = lon(ap_s,i)
                levs_out (writepos+count)   = lev(ap_s,i) 
                !tracer_out(writepos+count) = tracer(ap_s,i)
                !time_init_out(writepos+count)=time_init(ap_s,i)
                state_out (writepos+count)      = state(ap_s,i)
                state_vert_out (writepos+count) = state_vert(ap_s,i)
                count = count+1
             endif

          endif

       else

          if ((lev(ap_s,i) >= lev_min) .and. (lev(ap_s,i) <= lev_max) .and. &
               (lat(ap_s,i) >= lat_down) .and. (lat(ap_s,i) <= lat_up)) then 
          
             lats_out   (writepos+count) = lat(ap_s,i)
             lons_out   (writepos+count) = lon(ap_s,i)
             levs_out (writepos+count) = lev(ap_s,i) 
             !tracer_out(writepos+count)=tracer(ap_s,i)
             !time_init_out(writepos+count)=time_init(ap_s,i)
             state_out (writepos+count)      = state(ap_s,i)
             state_vert_out (writepos+count) = state_vert(ap_s,i)
             count = count+1
             
          endif

       endif

    enddo
    
    if (count == 0) then
      print *, 'No APs between lev_min, lev_max and above lat_down and lat_up'
      status = -1
      return
    else
       !if (ctrl_out) then
       !   write(*,fmt= '(a48,2(f6.0,2x))') 'Lev range (min, max): ', lev_min, lev_max
       !   write(*,fmt= '(a48,2(f6.0,2x))') 'Lat   range (up, down): ', lat_down, lat_up
       !   write(*,fmt= '(a48,i7)')         'Number of APs   : ', count
       !    write(*,fmt= '(a48,i7)')         'Total number of APs   : ', writepos-1 + count
       !endif
    endif

    !chemie
    count = 0

    do i=1,n 
       if (ctrl_subset) then
          if ((lev(ap_s,i) >= lev_min) .and. (lev(ap_s,i) <= lev_max) .and. &
               (lat(ap_s,i) >= lat_down) .and. (lat(ap_s,i) <= lat_up)) then 
             if (subset(ap_s,i)) then
                do k = 1, size(specarr_out)
                   specarr_out(k)%values(writepos+count) = c(ap_s,k,i)
                enddo
                count = count+1
             endif
          endif
       else
          if ((lev(ap_s,i) >= lev_min) .and. (lev(ap_s,i) <= lev_max) .and. &
               (lat(ap_s,i) >= lat_down) .and. (lat(ap_s,i) <= lat_up)) then 
             do k = 1, size(specarr_out)
                specarr_out(k)%values(writepos+count) = c(ap_s,k,i)
             enddo
             count = count+1
          endif
       endif
    end do
 
  end subroutine save_ap_s


end Module messy_clamsmix_lib_io

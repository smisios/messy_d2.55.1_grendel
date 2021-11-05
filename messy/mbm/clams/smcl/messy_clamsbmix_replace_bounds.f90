!********************************************************************************!
!
! This file is part of the Chemical Lagrangian Model of the
! Stratosphere (CLaMS). CLaMS is a hierarchy of numerical models and
! preprocessors which simulate transport, chemical reactions and
! mixing processes in the stratosphere.  
!
! Copyright (c) 2015
! Paul Konopka, Nicole Thomas
! Forschungszentrum Juelich GmbH
! Last Modified By: Nicole Thomas
! Last Modified On: Thu Dec 10 09:36:57 2020
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
! Module messy_clamsbix_replace_bounds
! ------------------------------------
!
!   subroutine replace_boundaries (status, lat, lon, lev, param, specarr)
!
!   subroutine interpol_bounds (status, lat, lon, lev, specarr, ispec, ind)
!
!   subroutine interpol_bounds_lat_time (status, lat, lev, specarr, &
!                                        ispec, ind, ncid, dimids, dimname)
!
!   subroutine interpol_bounds_lon_lat_time (status, lat, lon, lev, specarr, &
!                                            ispec, ind, ncid, dimids, dimname)
!
!   subroutine interpol_bounds_lat_lev (status, lat, lev, specarr, &
!                                       ispec, ind, ncid, dimids, dimname)
!
!   subroutine interpol_bounds_from_metdata (status, lat, lon, lev, specarr, &
!                                        ispec, ind, ncid) 
!
!----------------------------------------------------------------------
Module messy_clamsbmix_replace_bounds

Contains

  
  !****************************************************************************
  !
  !****************************************************************************
  Subroutine replace_boundaries (status, lat, lon, lev, param, specarr)

    USE messy_clams_global,        ONLY: prec, rank, dnparts, &
                                         paramtype, species_type, nspec, &
                                         YEAR_NEXT
    USE messy_clamsbmix_global,    ONLY: nclamsbounds, clamsbound

    implicit none

    REAL(PREC),         DIMENSION(:), POINTER :: lat, lon, lev
    TYPE(paramtype),    DIMENSION(:), POINTER :: param
    TYPE(species_type), DIMENSION(:), POINTER :: specarr
    INTEGER :: status


    INTEGER        :: ispec, ind, i

    if (rank==0) write (*,*) 
    if (rank==0) write (*,*) 'In REPLACE_BOUNDARIES'

    status = 0

    
    DO ispec = 1, nclamsbounds

       !----------------------------------------------------------------
       ! Replace only for given period of time
       !----------------------------------------------------------------
 
       IF (clamsbound(ispec)%startyear<=YEAR_NEXT .AND. YEAR_NEXT<=clamsbound(ispec)%endyear) THEN

          if (rank==0) write (*,*) 
          if (rank==0) write (*,*) 'replace ', clamsbound(ispec)%spec

          !----------------------------------------------------------------
          ! get index of species in specarr
          !----------------------------------------------------------------
          ind = -1
          do i = 1, nspec
             if (clamsbound(ispec)%spec==specarr(i)%name) ind = i
          enddo
          if (ind==-1) then
             if (rank==0) write (*,*) 'WARNING: Species ', &
                  trim(clamsbound(ispec)%spec),' not available !!!'
             cycle
          endif
        
          !----------------------------------------------------------------
          ! replace species below lower and/or above upper boundary with
          !    values interpolated from boundfile
          !----------------------------------------------------------------
          if (clamsbound(ispec)%lno==1 .or. clamsbound(ispec)%uno==1) then
             call interpol_bounds (status, lat, lon, lev, param, specarr, ispec, ind)
             if (status /= 0) return
          endif

          !----------------------------------------------------------------
          ! set to zero
          !----------------------------------------------------------------
          if (clamsbound(ispec)%lno==9 .or. clamsbound(ispec)%uno==9) then
             if (rank==0) write (*,*) 'set to zero'
             do i = 1, dnparts
                if ((clamsbound(ispec)%lno==9 .and. lev(i)<=clamsbound(ispec)%lbound) .or. &
                    (clamsbound(ispec)%uno==9 .and. lev(i)>=clamsbound(ispec)%ubound)) then
                   specarr(ind)%values(i)  = 0.
                endif
             enddo
          endif
         
       ENDIF 

    ENDDO

  End Subroutine replace_boundaries


  !****************************************************************************
  ! replace species below lower and/or above upper boundary
  !****************************************************************************
  subroutine interpol_bounds (status, lat, lon, lev, param, specarr, ispec, ind)

    USE messy_clams_global,        ONLY: YEAR_NEXT, MONTH_NEXT, DAY_NEXT, HOUR_NEXT

    USE messy_clams_global,        ONLY: prec, mdi, rank, met_prefix, met_dir, &
                                         species_type, init_vertcoorname, &
                                         theta_dir, theta_prefix, &
                                         paramtype, nparams, paramnames, &
                                         buffersize, filenamelen
    USE messy_clamsbmix_global,    ONLY: dir_boundfiles, nclamsbounds, clamsbound
    USE messy_clams_tools_utils,   ONLY: uppercase, lowercase, get_metfilename, str_pos
    USE messy_clams_tools_ncutils, ONLY: nc_check_error

    use netcdf

    implicit none

    REAL(PREC),         DIMENSION(:), POINTER :: lat, lon, lev
    TYPE(paramtype),    DIMENSION(:), POINTER :: param
    TYPE(species_type), DIMENSION(:), POINTER :: specarr
    INTEGER                                   :: status, ispec, ind
    
    
    REAL(PREC),               DIMENSION(:), POINTER       :: eqlat
    INTEGER,                  DIMENSION(:), ALLOCATABLE   :: dimids
    CHARACTER(nf90_max_name), DIMENSION(:), ALLOCATABLE   :: dimname

    INTEGER :: ncid, varid, ndims, idim, ipos
                    
    CHARACTER(filenamelen) :: filename
    CHARACTER(10)          :: datestr

    status = 0

    !----------------------------------------------------------------
    ! get filename
    !----------------------------------------------------------------

    filename = trim(clamsbound(ispec)%file)

    if (uppercase(filename)=='METDATA') then
          
       ! METEOROLOGICAL DATA is used for replacement of boundary:
       filename = get_metfilename &
            (met_prefix, met_dir, YEAR_NEXT, MONTH_NEXT, DAY_NEXT, HOUR_NEXT)
             
    elseif (filename(1:4)=='AIRS' .OR. filename(1:6)=='mopitt') then
             
       ! AIRS or MOPITT data: filename = filename + _yyyy010112_yyyy123112.nc
       if (filename(len_trim(filename)-2:len_trim(filename)) /= '.nc') then
          write (datestr,'(I4.4,3I2.2)') YEAR_NEXT,1,1,12
          filename = trim(adjustl(filename))//'_'//datestr
          write (datestr,'(I4.4,3I2.2)') YEAR_NEXT,12,31,12
          filename = trim(adjustl(filename))//'_'//datestr//'.nc'
       endif
       if (dir_boundfiles/=' ') filename = trim(dir_boundfiles)//'/'//trim(filename)
       
    else 
       
       if (dir_boundfiles/=' ') filename = trim(dir_boundfiles)//'/'//trim(filename)
       
    endif
    
    if (rank==0) write (*,*) 'Boundfile ',trim(filename)
    
    !----------------------------------------------------------------
    ! open boundfile
    !----------------------------------------------------------------

    status = nf90_open (filename,nf90_nowrite,ncid, buffersize)
    if (rank==0) call nc_check_error (status,'Error on open file '//trim(filename),abort=.false.)
    if (status/=0) return
    
    !----------------------------------------------------------------
    ! get number of dimensions for species in boundfile
    !----------------------------------------------------------------

    status = nf90_inq_varid (ncid, clamsbound(ispec)%spec_bf, varid)
    if (rank==0) call nc_check_error (status, &
         'Variable '//trim(clamsbound(ispec)%spec_bf)//' not found !',abort=.false.)
    if (status/=0) return
    status = nf90_inquire_variable (ncid, varid, ndims=ndims)
    if (rank==0) call nc_check_error (status, &
          'Cannot get dimensions for variable '//trim(clamsbound(ispec)%spec_bf)//' !', abort=.false.)
    if (status/=0) return
    !if (rank==0) write (*,*) 'ndims=',ndims
    allocate (dimids(ndims))
    status = nf90_inquire_variable (ncid, varid, dimids=dimids)
    if (rank==0) call nc_check_error (status, &
         'Cannot get dimensions for variable '//trim(clamsbound(ispec)%spec_bf)//' !',abort=.false.)
    if (status/=0) return
    
    !----------------------------------------------------------------
    ! get dimension for species from boundfile
    !----------------------------------------------------------------

    allocate (dimname(ndims))
    DO idim = 1, ndims
       status = nf90_inquire_dimension (ncid, dimids(idim), name=dimname(idim))
       if (rank==0) call nc_check_error (status, &
            'Cannot read dimensions  '//trim(dimname(idim))//' !',abort=.false.)
       if (status/=0)  return
       !if (rank==0) write (*,*) trim(dimname(idim))
    ENDDO
    
    !----------------------------------------------------------------
    ! call interpolation routines
    !----------------------------------------------------------------

    if (ndims == 2) then
       
       ! (lat, time) => interpol_bounds_lat_time
       
       if (lowercase(dimname(1))=='lat' .and. lowercase(dimname(2))=='time') then

          if (rank==0) write (*,*) 'call interpol_bounds_lat_time'
          call interpol_bounds_lat_time (status, lat, lev, specarr, &
                                         ispec, ind, ncid, dimids, dimname)
          if (status/=0 .and. rank==0) write (*,*) 'Error in subroutine interpol_bounds_lat_time'
          if (status/=0) return
 
       else
          if (rank==0) write (*,*) 'Type of boundfile unknown !!!'
          status = -2
          return
       endif

    elseif (ndims == 3) then
       
       ! (lat, lon, time) => interpol_bounds_lon_lat_time
       
       if (lowercase(dimname(1))=='lat' .and. lowercase(dimname(2))=='lon' &
            .and. lowercase(dimname(3))=='time') then
          
          if (rank==0) write (*,*) 'call interpol_bounds_lon_lat_time'
          call interpol_bounds_lon_lat_time (status, lat, lon, lev, specarr, &
                                         ispec, ind, ncid, dimids, dimname)
          if (status/=0 .and. rank==0) &
                  write (*,*) 'Error in subroutine interpol_bounds_lon_lat_time'
          if (status/=0) return
          
       ! (lat/eqlat, theta/zeta, month)   => interpol_bounds_lat_lev_month
          
       elseif ((lowercase(dimname(1))=='lat' .or. lowercase(dimname(1))=='eqlat') .and. &
               (lowercase(dimname(2))=='theta' .or. lowercase(dimname(2))=='zeta') .and.  &
                lowercase(dimname(3))=='month') then
          
          if (lowercase(dimname(2)) /= lowercase(init_vertcoorname)) then
             if (rank==0) then
                write (*,*) 'WARNING: vertical coordinate in initfile and boundfile are different!!!'
                write (*,*) 'Vertical coordinate in initfile: ',trim(init_vertcoorname)
                write (*,*) 'Vertical coordinate in boundfile: ',trim(dimname(2))
             endif
          endif
          
          if (lowercase(dimname(1))=='lat') then
             if (rank==0) write (*,*) 'call interpol_bounds_lat_lev_month (lat)'
             call interpol_bounds_lat_lev_month (status, lat, lev, specarr, &
                                           ispec, ind, ncid, dimids, dimname)
             if (status/=0 .and. rank==0) &
                  write (*,*) 'Error in subroutine interpol_bounds_lat_lev_month'
             if (status/=0) return

          else ! use EQLAT

             ipos = str_pos (nparams, paramnames, 'EQLAT')

             ! Use EQLAT from parameter list 
             if (ipos /= -1) then
                               
                if (rank==0) write (*,*) 'Use EQLAT from PARAM'
                if (rank==0) write (*,*) 'call interpol_bounds_lat_lev_month (eqlat)'
                call interpol_bounds_lat_lev_month (status, param(ipos)%values, lev, specarr, &
                                           ispec, ind, ncid, dimids, dimname)
                if (status/=0 .and. rank==0) &
                     write (*,*) 'Error in subroutine interpol_bounds_lat_lev_month'
                if (status/=0) return
                              
             else ! get EQLAT from met. file on theta levels
                
                allocate (eqlat(size(lat)))
                eqlat = mdi
                filename = get_metfilename (theta_prefix, theta_dir, &
                     YEAR_NEXT, MONTH_NEXT, DAY_NEXT, HOUR_NEXT)
                if (rank==0) write (*,*) 'Get EQLAT from file ',trim(filename)

                ! Check if THETA is on paramter list too
                ipos = str_pos (nparams, paramnames, 'THETA')
             
                if (ipos /= -1) then  ! use THETA
                   if (rank==0) write (*,*) 'nutze THETA fuer EQLAT-Interpolation'
                   call get_eqlat_from_metdata (status, filename, lat, lon, param(ipos)%values, eqlat)

                else  ! use vertical coordinate (ZETA)
                   if (rank==0) write (*,*) 'ACHTUNG: nutze ZETA fuer EQLAT-Interpolation !!!'
                   call get_eqlat_from_metdata (status, filename, lat, lon, lev, eqlat)
                endif

                if (status/=0 .and. rank==0) &
                        write (*,*) 'Error in subroutine get_eqlat_from_metdata'
                if (status/=0)  return
                
                if (rank==0) write (*,*) 'call interpol_bounds_lat_lev_month (eqlat)'
                call interpol_bounds_lat_lev_month (status, eqlat, lev, specarr, &
                                              ispec, ind, ncid, dimids, dimname)
                if (status/=0 .and. rank==0) &
                     write (*,*) 'Error in subroutine interpol_bounds_lat_lev_month'
                if (status/=0) return

                deallocate (eqlat)

             endif


          endif

       ! (lat/eqlat, theta/zeta, time) => interpol_bounds_lat_lev_time
          
       elseif ((lowercase(dimname(1))=='lat' .or. lowercase(dimname(1))=='eqlat') .and. &
               (lowercase(dimname(2))=='theta' .or. lowercase(dimname(2))=='zeta') .and.  &
                lowercase(dimname(3))=='time') then
          
          if (lowercase(dimname(2)) /= lowercase(init_vertcoorname)) then
             if (rank==0) then
                write (*,*) 'WARNING: vertical coordinate in initfile and boundfile are different!!!'
                write (*,*) 'Vertical coordinate in initfile: ',trim(init_vertcoorname)
                write (*,*) 'Vertical coordinate in boundfile: ',trim(dimname(2))
             endif
          endif
          
          if (lowercase(dimname(1))=='lat') then
             if (rank==0) write (*,*) 'call interpol_bounds_lat_lev_time (lat)'
             call interpol_bounds_lat_lev_time (status, lat, lev, specarr, &
                                                ispec, ind, ncid, dimids, dimname)
             if (status/=0 .and. rank==0) &
                  write (*,*) 'Error in subroutine interpol_bounds_lat_lev_time'
             if (status/=0) return

          else ! use EQLAT

             ipos = str_pos (nparams, paramnames, 'EQLAT')

             ! Use EQLAT from parameter list 
             if (ipos /= -1) then
                               
                if (rank==0) write (*,*) 'Use EQLAT from PARAM'
                if (rank==0) write (*,*) 'call interpol_bounds_lat_lev_time (eqlat)'
                call interpol_bounds_lat_lev_time (status, param(ipos)%values, lev, specarr, &
                                                    ispec, ind, ncid, dimids, dimname)
                if (status/=0 .and. rank==0) &
                     write (*,*) 'Error in subroutine interpol_bounds_lat_lev_time'
                if (status/=0) return
                              
             else ! get EQLAT from met. file on theta levels
                
                allocate (eqlat(size(lat)))
                eqlat = mdi
                filename = get_metfilename (theta_prefix, theta_dir, &
                                            YEAR_NEXT, MONTH_NEXT, DAY_NEXT, HOUR_NEXT)
                if (rank==0) write (*,*) 'Get EQLAT from file ',trim(filename)

                ! Check if THETA is on paramter list too
                ipos = str_pos (nparams, paramnames, 'THETA')
             
                if (ipos /= -1) then  ! use THETA
                   if (rank==0) write (*,*) 'nutze THETA fuer EQLAT-Interpolation'
                   call get_eqlat_from_metdata (status, filename, lat, lon, param(ipos)%values, eqlat)

                else  ! use vertical coordinate 
                   if (rank==0) write (*,*) 'ACHTUNG: nutze ',trim(init_vertcoorname),' fuer EQLAT-Interpolation !!!'
                   call get_eqlat_from_metdata (status, filename, lat, lon, lev, eqlat)
                endif

                if (status/=0 .and. rank==0) &
                        write (*,*) 'Error in subroutine get_eqlat_from_metdata'
                if (status/=0)  return
                
                if (rank==0) write (*,*) 'call interpol_bounds_lat_lev_time (eqlat)'
                call interpol_bounds_lat_lev_time (status, eqlat, lev, specarr, &
                                                    ispec, ind, ncid, dimids, dimname)
                if (status/=0 .and. rank==0) &
                     write (*,*) 'Error in subroutine interpol_bounds_lat_lev_time'
                if (status/=0) return

                deallocate (eqlat)

             endif

          endif
          
       else
          if (rank==0) write (*,*) 'Type of boundfile unknown !!!'
          status = -2
          return
       endif
       
    elseif (ndims == 4) then
       
       ! (lev, lat, lon, time) => interpol_bounds_from_metdata
       
       if (uppercase(clamsbound(ispec)%file)=='METDATA') then

          if (rank==0) write (*,*) 'call interpol_bounds_from_metdata'
          
          call interpol_bounds_from_metdata (status, lat, lon, lev, specarr, &
                                             ispec, ind, ncid) 
          if (status/=0 .and. rank==0) &
                  write (*,*) 'Error in subroutine interpol_bound_from_metdata'
          if (status/=0) return
          
       else
          if (rank==0) write (*,*) 'Type of boundfile unknown !!!'
          status = -2
          return
       endif
       
    else
       
       if (rank==0) write (*,*) 'Number of dimensions in Boundfile must be 2, 3 or 4 !!!'
       status = -2
       return             
       
    endif
    
    !----------------------------------------------------------------
    ! close boundfile
    !----------------------------------------------------------------
    status = nf90_close (ncid)
    
    !----------------------------------------------------------------
    ! clean up
    !----------------------------------------------------------------
    deallocate (dimname)
    deallocate (dimids)
    
  end subroutine interpol_bounds

  !****************************************************************************
  ! 
  !****************************************************************************
  subroutine interpol_bounds_lat_time (status, lat, lev, specarr, &
                                       ispec, ind, ncid, dimids, dimname)

    USE messy_clams_global,         ONLY: YEAR_NEXT, MONTH_NEXT, DAY_NEXT, HOUR_NEXT

    USE messy_clams_global,         ONLY: PREC, DP, rank, dnparts, species_type
    USE messy_clamsbmix_global,     ONLY: clamsbound
    USE messy_clamsbmix_tools,      ONLY: grid_index, time_index
    USE messy_clams_tools_dateconv, ONLY: ymds2js

    use netcdf

    implicit none

    REAL(PREC),         DIMENSION(:), POINTER :: lat, lev
    TYPE(species_type), DIMENSION(:), POINTER :: specarr
    INTEGER      :: status, ispec, ind, ncid
    INTEGER      :: dimids(2)
    CHARACTER(*) :: dimname(2)

    INTEGER                                 :: it
    INTEGER,    DIMENSION(:),   ALLOCATABLE :: ix
    REAL(PREC)                              :: dt
    REAL(PREC), DIMENSION(:),   ALLOCATABLE :: blat, dx
    REAL(DP),   DIMENSION(:),   ALLOCATABLE :: btime
    REAL(PREC), DIMENSION(:,:), ALLOCATABLE :: bdata

    REAL(DP) :: time

    REAL     :: x1, x2, y1, y2

    INTEGER  :: nlat, ntime, varid, i

    status = 0

    time = ymds2js (YEAR_NEXT, MONTH_NEXT, DAY_NEXT, HOUR_NEXT*3600)

    !--------------------------------------------------------------
    ! read dimensions 
    !--------------------------------------------------------------
    status = nf90_inquire_dimension(ncid,dimids(1),len=nlat)
    if (status /= 0) return
    status = nf90_inquire_dimension(ncid,dimids(2),len=ntime)
    if (status /= 0) return

    !--------------------------------------------------------------
    ! allocate arrays
    !--------------------------------------------------------------
    allocate (blat(nlat), btime(ntime))
    allocate (bdata(nlat, ntime))
    allocate (dx(dnparts),ix(dnparts))

    !--------------------------------------------------------------
    ! get coordinate variables
    !--------------------------------------------------------------
    status = nf90_inq_varid (ncid, dimname(1), varid)
    status = nf90_get_var (ncid, varid, blat)
    if (status /= 0) return
    status = nf90_inq_varid (ncid, dimname(2), varid)
    status = nf90_get_var (ncid, varid, btime)
    if (status /= 0) return
    
    ! if (rank==0) write (*,*) 'blat=',blat
    ! if (rank==0) write (*,*) 'btime=',btime

    !--------------------------------------------------------------
    ! get indices in lat and time
    !--------------------------------------------------------------
    ! get lat indices
    call grid_index (blat, lat(1:dnparts), ix, dx)
    
    ! get time index
    call time_index (btime, time, it, dt)

    !--------------------------------------------------------------
    ! read species from boundfile
    !--------------------------------------------------------------
    status = nf90_inq_varid (ncid, trim(clamsbound(ispec)%spec_bf), varid)
    !if (rank==0) write (*,*) 'read variable ',trim(clamsbound(ispec)%spec_bf)
    status = nf90_get_var (ncid, varid, bdata)
    if (status /= 0) then
       if (rank==0) write (*,*) 'Cannot read ',trim(clamsbound(ispec)%spec_bf)
       return
    endif


    !--------------------------------------------------------------
    ! replace species below lower and above upper boundary
    !--------------------------------------------------------------
    do i = 1, dnparts

       if ((clamsbound(ispec)%lno==1 .and. lev(i)<=clamsbound(ispec)%lbound) .or. &
           (clamsbound(ispec)%uno==1 .and. lev(i)>=clamsbound(ispec)%ubound)) then
        
          !--------------------------------------------------------------
          ! interpolate in lat and time
          !--------------------------------------------------------------
          x1 = bdata(ix(i),it) + &
               dx(i) * ( bdata(ix(i)+1,it) - bdata(ix(i),it))
          x2 = bdata(ix(i),it+1) + &
               dx(i) * ( bdata(ix(i)+1,it+1) - bdata(ix(i),it+1))
          
          specarr(ind)%values(i)  = x1 + dt* (x2 - x1) 

       endif

    enddo

    !--------------------------------------------------------------
    ! deallocate arrays
    !--------------------------------------------------------------
    deallocate (blat, btime, dx, ix)
    deallocate (bdata)


  end subroutine interpol_bounds_lat_time

  !****************************************************************************
  ! 
  !****************************************************************************
  subroutine interpol_bounds_lon_lat_time (status, lat, lon, lev, specarr, &
                                           ispec, ind, ncid, dimids, dimname)

    USE messy_clams_global,         ONLY: YEAR_NEXT, MONTH_NEXT, DAY_NEXT, HOUR_NEXT

    USE messy_clams_global,         ONLY: PREC, DP, rank, dnparts, species_type
    USE messy_clamsbmix_global,     ONLY: clamsbound
    USE messy_clamsbmix_tools,      ONLY: grid_index, time_index
    USE messy_clams_tools_dateconv, ONLY: ymds2js

    USE netcdf

    implicit none 

    REAL(PREC),         DIMENSION(:), POINTER :: lat, lon, lev
    TYPE(species_type), DIMENSION(:), POINTER :: specarr
    INTEGER      :: status, ispec, ind, ncid
    INTEGER      :: dimids(3)
    CHARACTER(*) :: dimname(3)

    INTEGER                                   :: it
    INTEGER,    DIMENSION(:),     ALLOCATABLE :: ix, iy
    REAL(PREC)                                :: dt
    REAL(PREC), DIMENSION(:),     ALLOCATABLE :: blat, blon, dx, dy
    REAL(DP),   DIMENSION(:),     ALLOCATABLE :: btime
    REAL(PREC), DIMENSION(:,:,:), ALLOCATABLE :: bdata

    REAL(DP) :: time

    REAL     :: x1, x2, y1, y2

    INTEGER  :: nlat, nlon, ntime, varid, i

    status = 0

    time = ymds2js (YEAR_NEXT, MONTH_NEXT, DAY_NEXT, HOUR_NEXT*3600)

    !--------------------------------------------------------------
    ! read dimensions
    !--------------------------------------------------------------
    status = nf90_inquire_dimension(ncid,dimids(1),len=nlat)
    if (status /= 0) return
    status = nf90_inquire_dimension(ncid,dimids(2),len=nlon)
    if (status /= 0) return
    status = nf90_inquire_dimension(ncid,dimids(3),len=ntime)
    if (status /= 0) return

    !--------------------------------------------------------------
    ! allocate arrays
    !--------------------------------------------------------------
    allocate (blat(nlat), blon(nlon), btime(ntime))
    allocate (bdata(nlat, nlon,ntime))
    allocate (dx(dnparts),ix(dnparts))
    allocate (dy(dnparts),iy(dnparts))

    !--------------------------------------------------------------
    ! get coordinate variables
    !--------------------------------------------------------------
    status = nf90_inq_varid (ncid, dimname(1), varid)
    status = nf90_get_var (ncid, varid, blat)
    if (status /= 0) return
    status = nf90_inq_varid (ncid, dimname(2), varid)
    status = nf90_get_var (ncid, varid, blon)
    if (status /= 0) return
    status = nf90_inq_varid (ncid, dimname(3), varid)
    status = nf90_get_var (ncid, varid, btime)
    if (status /= 0) return
    
    ! if (rank==0) write (*,*) 'blat=',blat
    ! if (rank==0) write (*,*) 'blon=',blon
    ! if (rank==0) write (*,*) 'btime=',btime

    !--------------------------------------------------------------
    ! get indices in lat, lon and time
    !--------------------------------------------------------------
    ! get lat indices
    call grid_index (blat, lat(1:dnparts), ix, dx)
    
    ! get lon indices
    call grid_index (blon, lon(1:dnparts), iy, dy)

    ! get time index
    call time_index (btime, time, it, dt)

    !--------------------------------------------------------------
    ! read species from boundfile
    !--------------------------------------------------------------
    status = nf90_inq_varid (ncid, trim(clamsbound(ispec)%spec_bf), varid)
    !if (rank==0) write (*,*) 'read variable ',trim(clamsbound(ispec)%spec_bf)
    status = nf90_get_var (ncid, varid, bdata)
    if (status /= 0) then
       if (rank==0) write (*,*) 'Cannot read ',trim(clamsbound(ispec)%spec_bf)
       return
    endif


    !--------------------------------------------------------------
    ! replace species below lower and above upper boundary
    !--------------------------------------------------------------
    do i = 1, dnparts

       if ((clamsbound(ispec)%lno==1 .and. lev(i)<=clamsbound(ispec)%lbound) .or. &
           (clamsbound(ispec)%uno==1 .and. lev(i)>=clamsbound(ispec)%ubound)) then

          !--------------------------------------------------------------
          ! Interpolate in lon, lat and time
          !--------------------------------------------------------------
          y1 = bdata(ix(i),iy(i),it) + &
               dy(i) * (bdata(ix(i),iy(i)+1,it) - bdata(ix(i),iy(i),it))
          y2 = bdata(ix(i)+1,iy(i),it) + &
               dy(i) * (bdata(ix(i)+1,iy(i)+1,it) - bdata(ix(i)+1,iy(i),it))
          x1  = y1 + dx(i) * (y2 -y1)
          
          y1 = bdata(ix(i),iy(i),it+1) + &
               dy(i) * (bdata(ix(i),iy(i)+1,it+1) - bdata(ix(i),iy(i),it+1))
          y2 = bdata(ix(i)+1,iy(i),it+1) + &
               dy(i) * (bdata(ix(i)+1,iy(i)+1,it+1) - bdata(ix(i)+1,iy(i),it+1))
          x2  = y1 + dx(i) * (y2 -y1)
          
          specarr(ind)%values(i)  = x1 + dt* (x2 - x1) 

       endif

    enddo

    !--------------------------------------------------------------
    ! deallocate arrays
    !--------------------------------------------------------------
    deallocate (blat, blon, btime, dx, ix, dy, iy)
    deallocate (bdata)

  end subroutine interpol_bounds_lon_lat_time


  !****************************************************************************
  ! 
  !****************************************************************************
  subroutine interpol_bounds_lat_lev_month (status, lat, lev, specarr, &
                                      ispec, ind, ncid, dimids, dimname)

    USE messy_clams_global,         ONLY: MONTH_NEXT
    USE messy_clams_global,         ONLY: PREC, dnparts, species_type, rank, &
                                          mdi, eps
    USE messy_clamsbmix_global,     ONLY: clamsbound

    use netcdf

    implicit none

    REAL(PREC),         DIMENSION(:), POINTER :: lat, lev
    TYPE(species_type), DIMENSION(:), POINTER :: specarr
    INTEGER      :: status, ispec, ind, ncid
    INTEGER      :: dimids(3)
    CHARACTER(*) :: dimname(3)

    REAL(PREC), DIMENSION(:),     ALLOCATABLE :: blat, blev
    REAL(PREC), DIMENSION(:,:,:), ALLOCATABLE :: bdata
    REAL(PREC)   :: lat0, latd, latx
    REAL(PREC)   :: dlat, dlev, x1, x2
    INTEGER      :: nlat, nlev, nmonth, ilat, ilev, i, k
    INTEGER      :: varid
    

    status = 0

    !--------------------------------------------------------------
    ! read dimensions
    !--------------------------------------------------------------
    status = nf90_inquire_dimension(ncid,dimids(1),len=nlat)
    if (status /= 0) return
    status = nf90_inquire_dimension(ncid,dimids(2),len=nlev)
    if (status /= 0) return
    status = nf90_inquire_dimension(ncid,dimids(3),len=nmonth)
    if (status /= 0) return

    !--------------------------------------------------------------
    ! allocate arrays
    !--------------------------------------------------------------
    allocate (blat(nlat), blev(nlev))
    allocate (bdata(nlat,nlev,nmonth))

    !--------------------------------------------------------------
    ! read coordinate variable lat/eqlat
    !--------------------------------------------------------------
    status = nf90_inq_varid (ncid, dimname(1), varid)
    status = nf90_get_var (ncid, varid, blat)
    if (status /= 0) return

    lat0 = blat(1)
    latd = blat(2) - blat(1)

    !--------------------------------------------------------------
    ! read coordinate variable theta/zeta
    !--------------------------------------------------------------
    status = nf90_inq_varid (ncid, dimname(2), varid)
    status = nf90_get_var (ncid, varid, blev)
    if (status /= 0) return

    !--------------------------------------------------------------
    ! read species from boundfile
    !--------------------------------------------------------------
    status = nf90_inq_varid (ncid, trim(clamsbound(ispec)%spec_bf), varid)
    !if (rank==0) write (*,*) 'read variable ',trim(clamsbound(ispec)%spec_bf)
    status = nf90_get_var (ncid, varid, bdata)
    if (status /= 0) then
       if (rank==0) write (*,*) 'Cannot read ',trim(clamsbound(ispec)%spec_bf)
       return
    endif

    !--------------------------------------------------------------
    ! replace species below lower and above upper boundary
    !--------------------------------------------------------------
    do i = 1, dnparts


       if ((clamsbound(ispec)%lno==1 .and. lev(i)<=clamsbound(ispec)%lbound) .or. &
           (clamsbound(ispec)%uno==1 .and. lev(i)>=clamsbound(ispec)%ubound)) then

          if (ABS((lat(i)-mdi)/mdi)<=eps .or. ABS((lev(i)-mdi)/mdi)<=eps) then
             specarr(ind)%values(i) = mdi
             cycle
          endif

          !-------------------------------------------------------------
          ! get lat index 
          !--------------------------------------------------------------
          ilat = int((lat(i)-lat0)/latd)+1
          ilat = min(ilat,nlat-1)
          ilat = max(ilat,1)
          latx = lat(i)
          if (latx > maxval(blat)) latx = maxval(blat)
          if (latx < minval(blat)) latx = minval(blat)
          dlat = (latx-blat(ilat))/(blat(ilat+1)-blat(ilat))

          !--------------------------------------------------------------
          ! get lev index
          !--------------------------------------------------------------
          if (lev(i) <= blev(1)) then
             ilev = 1
             dlev = 0.
          elseif (lev(i) >= blev(nlev)) then
             ilev = nlev-1
             dlev = 1.
          else
             do k = 1, nlev-1
                if (blev(k)<=lev(i) .and. lev(i)<=blev(k+1)) ilev = k
             enddo
             dlev = (lev(i)-blev(ilev))/(blev(ilev+1)-blev(ilev))
          endif

          !--------------------------------------------------------------
          ! interpolate in latitude and level
          !--------------------------------------------------------------
          x1 = bdata(ilat,ilev,MONTH_NEXT) + &
               dlat * ( bdata(ilat+1,ilev,MONTH_NEXT) - bdata(ilat,ilev,MONTH_NEXT))
          x2 = bdata(ilat,ilev+1,MONTH_NEXT) + &
               dlat * ( bdata(ilat+1,ilev+1,MONTH_NEXT) - bdata(ilat,ilev+1,MONTH_NEXT))
        
          specarr(ind)%values(i)  = x1 + dlev * (x2 - x1) 

       endif

    enddo

    !--------------------------------------------------------------
    ! clean up
    !--------------------------------------------------------------
    deallocate (blat, blev, bdata)


  end subroutine interpol_bounds_lat_lev_month
 
  !****************************************************************************
  ! 
  !****************************************************************************
  subroutine interpol_bounds_lat_lev_time (status, lat, lev, specarr, &
                                           ispec, ind, ncid, dimids, dimname)

    USE messy_clams_global,         ONLY: YEAR_NEXT, MONTH_NEXT, DAY_NEXT, HOUR_NEXT
    USE messy_clams_global,         ONLY: PREC, DP, dnparts, species_type, rank, &
                                          mdi, eps
    USE messy_clamsbmix_global,     ONLY: clamsbound
    USE messy_clamsbmix_tools,      ONLY: time_index
    USE messy_clams_tools_dateconv, ONLY: ymds2js

    use netcdf

    implicit none

    REAL(PREC),         DIMENSION(:), POINTER :: lat, lev
    TYPE(species_type), DIMENSION(:), POINTER :: specarr
    INTEGER      :: status, ispec, ind, ncid
    INTEGER      :: dimids(3)
    CHARACTER(*) :: dimname(3)

    REAL(PREC), DIMENSION(:),     ALLOCATABLE :: blat, blev
    REAL(DP),   DIMENSION(:),     ALLOCATABLE :: btime
    REAL(PREC), DIMENSION(:,:,:), ALLOCATABLE :: bdata
    REAL(PREC)   :: lat0, latd, latx
    REAL(PREC)   :: dlat, dlev, dt, x1, x2, y1, y2
    REAL(DP)     :: time
    INTEGER      :: nlat, nlev, ntime, ilat, ilev, itime, i, k
    INTEGER      :: varid
    
    status = 0

    time = ymds2js (YEAR_NEXT, MONTH_NEXT, DAY_NEXT, HOUR_NEXT*3600)

    !--------------------------------------------------------------
    ! read dimensions
    !--------------------------------------------------------------
    status = nf90_inquire_dimension(ncid,dimids(1),len=nlat)
    if (status /= 0) return
    status = nf90_inquire_dimension(ncid,dimids(2),len=nlev)
    if (status /= 0) return
    status = nf90_inquire_dimension(ncid,dimids(3),len=ntime)
    if (status /= 0) return

    !--------------------------------------------------------------
    ! allocate arrays
    !--------------------------------------------------------------
    allocate (blat(nlat), blev(nlev), btime(ntime))
    allocate (bdata(nlat,nlev,ntime))

    !--------------------------------------------------------------
    ! read coordinate variable lat/eqlat
    !--------------------------------------------------------------
    status = nf90_inq_varid (ncid, dimname(1), varid)
    status = nf90_get_var (ncid, varid, blat)
    if (status /= 0) return

    lat0 = blat(1)
    latd = blat(2) - blat(1)

    !--------------------------------------------------------------
    ! read coordinate variable theta/zeta
    !--------------------------------------------------------------
    status = nf90_inq_varid (ncid, dimname(2), varid)
    status = nf90_get_var (ncid, varid, blev)
    if (status /= 0) return

    !--------------------------------------------------------------
    ! read coordinate variable time
    !--------------------------------------------------------------
    status = nf90_inq_varid (ncid, dimname(3), varid)
    status = nf90_get_var (ncid, varid, btime)
    if (status /= 0) return
    
    !--------------------------------------------------------------
    ! read species from boundfile
    !--------------------------------------------------------------
    status = nf90_inq_varid (ncid, trim(clamsbound(ispec)%spec_bf), varid)
    !if (rank==0) write (*,*) 'read variable ',trim(clamsbound(ispec)%spec_bf)
    status = nf90_get_var (ncid, varid, bdata)
    if (status /= 0) then
       if (rank==0) write (*,*) 'Cannot read ',trim(clamsbound(ispec)%spec_bf)
       return
    endif

    !--------------------------------------------------------------
    ! replace species below lower and above upper boundary
    !--------------------------------------------------------------
    do i = 1, dnparts

       if ((clamsbound(ispec)%lno==1 .and. lev(i)<=clamsbound(ispec)%lbound) .or. &
           (clamsbound(ispec)%uno==1 .and. lev(i)>=clamsbound(ispec)%ubound)) then

          if (ABS((lat(i)-mdi)/mdi)<=eps .or. ABS((lev(i)-mdi)/mdi)<=eps) then
             specarr(ind)%values(i) = mdi
             cycle
          endif

          !--------------------------------------------------------------
          ! get lat index 
          !--------------------------------------------------------------
          ilat = int((lat(i)-lat0)/latd)+1
          ilat = min(ilat,nlat-1)
          ilat = max(ilat,1)
          latx = lat(i)
          if (latx > maxval(blat)) latx = maxval(blat)
          if (latx < minval(blat)) latx = minval(blat)
          dlat = (latx-blat(ilat))/(blat(ilat+1)-blat(ilat))

          !--------------------------------------------------------------
          ! get lev index
          !--------------------------------------------------------------
          if (lev(i) <= blev(1)) then
             ilev = 1
             dlev = 0.
          elseif (lev(i) >= blev(nlev)) then
             ilev = nlev-1
             dlev = 1.
          else
             do k = 1, nlev-1
                if (blev(k)<=lev(i) .and. lev(i)<=blev(k+1)) ilev = k
             enddo
             dlev = (lev(i)-blev(ilev))/(blev(ilev+1)-blev(ilev))
          endif

          !--------------------------------------------------------------
          ! get time index
          !--------------------------------------------------------------
          call time_index (btime, time, itime, dt)
          
          !--------------------------------------------------------------
          ! Interpolate in lat, lev and time
          !--------------------------------------------------------------
          x1 = bdata(ilat,ilev,itime) + &
               dlat * (bdata(ilat+1,ilev,itime) - bdata(ilat,ilev,itime))
          x2 = bdata(ilat,ilev+1,itime) + &
               dlat * (bdata(ilat+1,ilev+1,itime) - bdata(ilat,ilev+1,itime))
          y1  = x1 + dlev * (x2 - x1)
          
          x1 = bdata(ilat,ilev,itime+1) + &
               dlat * (bdata(ilat+1,ilev,itime+1) - bdata(ilat,ilev,itime+1))
          x2 = bdata(ilat,ilev+1,itime+1) + &
               dlat * (bdata(ilat+1,ilev+1,itime+1) - bdata(ilat,ilev+1,itime+1))
          y2  = x1 + dlev * (x2 - x1)
          
          specarr(ind)%values(i)  = y1 + dt * (y2 - y1)
          
       endif

    enddo

    !--------------------------------------------------------------
    ! clean up
    !--------------------------------------------------------------
    deallocate (blat, blev, bdata)


  end subroutine interpol_bounds_lat_lev_time
 

  !****************************************************************************
  ! 
  !****************************************************************************
  subroutine interpol_bounds_from_metdata (status, lat, lon, lev, specarr, &
                                           ispec, ind, ncid) 


    USE messy_clams_global,            ONLY: PREC, dnparts, rank, species_type, &
                                             latgrid, longrid, levelgrid, leveldt, &
                                             nx, ny, nz, &
                                             asc_lat, asc_level, loglev, logpress, &
                                             level_is_vertcoor, &
                                             mdi, eps
    USE messy_clamsbmix_global,        ONLY: clamsbound
    USE messy_clams_tools_interpolreg, ONLY: interpolate_spatial, &
                                             interpolate_spatial_modlev
    USE messy_clams_tools_ncutils,     ONLY: nc_check_error, nc_get_var_cf
    

    use netcdf

    implicit none

    REAL(PREC),         DIMENSION(:), POINTER :: lat, lon, lev
    TYPE(species_type), DIMENSION(:), POINTER :: specarr
    INTEGER                                   :: status, ispec, ind, ncid


    REAL(PREC), DIMENSION(:,:,:), POINTER :: data
    REAL(PREC)                            :: factor
    INTEGER                               :: logval
    INTEGER                               :: varid, i

    status = 0
    
    !--------------------------------------------------------------
    ! read species from boundfile
    !--------------------------------------------------------------
    allocate (data(nx,ny,nz))
    status = nf90_inq_varid (ncid, trim(clamsbound(ispec)%spec_bf), varid)
    if (rank==0) call nc_check_error (status, &
         "Cannot find "//trim(clamsbound(ispec)%spec_bf),abort=.false.)
    call nc_get_var_cf (status, ncid, varid, trim(clamsbound(ispec)%spec_bf), data)
    if (status /= 0) return

    DO i = 1, dnparts
   
       !--------------------------------------------------------------
       ! parameter interpolation linear or log. linear ?
       !--------------------------------------------------------------
       logval = 0  ! lin. interpolation
       if (clamsbound(ispec)%spec=='SH' .or. &
           clamsbound(ispec)%spec=='H2O' .or. &
           clamsbound(ispec)%spec=='O3') then
          logval = 1  ! log. lin. interpolation
       elseif (clamsbound(ispec)%spec=='PRESS' .and. logpress) then
          logval = 2  ! lin. + log. lin. interpolation
       endif

       !--------------------------------------------------------------
       ! set conversion factor
       !--------------------------------------------------------------
       factor = 1.
       ! H2O = (28.9644/18.015) * SH_ECMWF
       if (clamsbound(ispec)%spec=='H2O') then
          factor = 28.9644/18.015
       ! O3 =  (28.9644/48.0) * O3_ECMWF
       elseif (clamsbound(ispec)%spec=='O3') then
          factor = 28.9644/48.0
       endif


       !----------------------------------------------------------------
       ! replace species below lower and above upper boundary
       !----------------------------------------------------------------
       if ((clamsbound(ispec)%lno==1 .and. lev(i)<=clamsbound(ispec)%lbound) .or. &
           (clamsbound(ispec)%uno==1 .and. lev(i)>=clamsbound(ispec)%ubound)) then

          if (.not. level_is_vertcoor) then

             ! interpolate from data on model levels
             specarr(ind)%values(i) = interpolate_spatial_modlev (data, &
                                         leveldt, longrid, latgrid, &
                                         nx, ny, nz, lon(i), lat(i), lev(i), &
                                         asc_lat, asclev=asc_level, loglev=loglev, logval=logval, &
                                         extrapol=.false.)

          else

             ! interpolate from data on theta/zeta/press levels
             specarr(ind)%values(i) = interpolate_spatial (data, &
                                  levelgrid, longrid, latgrid, &
                                  nx, ny, nz, lon(i), lat(i), lev(i), &
                                  asc_lat, asclev=asc_level, loglev=loglev, logval=logval, &
                                  extrapol=.false.)

          endif
 
          ! multiply with conversion factor
          if (abs((specarr(ind)%values(i)-mdi)/mdi)>eps) then
             specarr(ind)%values(i) = factor * specarr(ind)%values(i)
          endif

       endif

    ENDDO

    !--------------------------------------------------------------
    ! clean up
    !--------------------------------------------------------------
    deallocate (data)


  end subroutine interpol_bounds_from_metdata


  !*************************************************************************************
  !
  !*************************************************************************************
  subroutine get_eqlat_from_metdata (status, filename, lat, lon, lev, eqlat)

    USE messy_clams_global,            ONLY: prec, rank, dnparts, buffersize, &
                                             asc_lat, asc_level
    USE messy_clams_tools_ncutils,     ONLY: nc_check_error, nc_get_var_cf
    USE messy_clams_tools_interpolreg, ONLY: interpolate_spatial

    USE netcdf
    
    implicit none

    REAL(PREC), DIMENSION(:), POINTER :: lat, lon, lev, eqlat
    INTEGER                           :: status
    CHARACTER(*)                      :: filename

    REAL(PREC), DIMENSION(:,:,:), POINTER :: data
    REAL(PREC), DIMENSION(:),     POINTER :: longrid, latgrid, levgrid
    INTEGER :: nlev, nlat, nlon
    INTEGER :: ncid, varid, dimid, i
    

    ! open netcdf file (with met. data on theta levels)
    status = nf90_open (filename, nf90_nowrite, ncid, buffersize)
    if (rank==0) call nc_check_error (status,"Error on opening file "//trim(filename),abort=.false.)
    if (status/=0) return
    
    ! get dimensions (lat, lon, theta)
    status = nf90_inq_dimid (ncid, "lat", dimid)
    status = nf90_inquire_dimension (ncid, dimid, len=nlat)
    if (rank==0) call nc_check_error (status,"Cannot find dimension lat",abort=.false.)
    if (status/=0) return
    status = nf90_inq_dimid (ncid, "lon", dimid)
    status = nf90_inquire_dimension (ncid, dimid, len=nlon)
    if (rank==0) call nc_check_error (status,"Cannot find dimension lon",abort=.false.)
    if (status/=0) return
    status = nf90_inq_dimid (ncid, "theta", dimid)
    status = nf90_inquire_dimension (ncid, dimid, len=nlev)
    if (rank==0) call nc_check_error (status,"Cannot find dimension theta",abort=.false.)
    if (status/=0) return

    ! allocate arrays
    allocate (longrid(nlon))
    allocate (latgrid(nlat))
    allocate (levgrid(nlev))
    allocate (data(nlon,nlat,nlev))

    ! read coordinate variables (lon, lat, theta)
    status = nf90_inq_varid (ncid, 'lon', varid)
    status = nf90_get_var (ncid, varid, longrid)
    if (rank==0) call nc_check_error (status,'Cannot read variable lon',abort=.false.)
    if (status/=0) return
    status = nf90_inq_varid (ncid, 'lat', varid)
    status = nf90_get_var (ncid, varid, latgrid)
    if (rank==0) call nc_check_error (status,'Cannot read variable lat',abort=.false.)
    if (status/=0) return
    status = nf90_inq_varid (ncid, 'theta', varid)
    status = nf90_get_var (ncid, varid, levgrid)
    if (rank==0) call nc_check_error (status,'Cannot read variable theta',abort=.false.)
    if (status/=0) return

    ! read EQLAT
    status = nf90_inq_varid (ncid, "EQLAT", varid)
    if (rank==0) call nc_check_error (status,"Cannot find EQLAT",abort=.false.)
    if (status/=0) return
    call nc_get_var_cf (status, ncid, varid, "EQLAT", data)
    if (status/=0) return

    ! interpolate EQLAT from data on theta levels
    DO i = 1, dnparts
!!!!!! ACHTUNG auf lev(i) steht evtl. zeta (NICHT theta) !
        eqlat(i) = interpolate_spatial (data, levgrid, longrid, latgrid, &
                                       nlon, nlat, nlev, lon(i), lat(i), lev(i), &
                                       asc_lat, asclev=asc_level, loglev=.false., logval=0, &
                                       extrapol=.false.)
    ENDDO

    ! close netcdf file
    status = nf90_close (ncid)

    ! clean up
    deallocate (longrid, latgrid, levgrid)
    deallocate (data)

  end subroutine get_eqlat_from_metdata

End Module messy_clamsbmix_replace_bounds

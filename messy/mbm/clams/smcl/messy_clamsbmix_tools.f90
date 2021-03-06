!*****************************************************************************!
!
! MODULE messy_clamsbmix_tools
! -------------------------------
!
! This MODULE contains subroutines used by bmix: 
!
!   SUBROUTINE pack_array (array, mask)
!   subroutine get_nparts_layer (lev, indices, lev_min, lev_max)
!   subroutine time_index (btime, time, it, dt)
!   subroutine grid_index (bvec, vec, ix, dx)
!   subroutine interpol_time_grid (status, specarr, ncid, dim1, dim2, &
!                                  nparts_old, nparts_new,  &
!                                  it, dt, ix, dx, replace_or_add)
!   subroutine interpol_time_grid_3d (status, specarr, ncid, dim1, dim2, dim3, &
!                                  nparts_old, nparts_new,  &
!                                  it, dt, ix, dx, iy, dy, replace_or_add)
!   subroutine get_boundfile_indices (status, boundfile, &
!                                   lat_min, lat_max, lev_min, lev_max, &
!                                   indices, use_oro)
!   subroutine read_boundfile (status, boundfile, indices, lat, lon, lev, specarr)
!   subroutine get_bound_species (status, lat, lon, eqlat, specarr, speclist, &
!                                 startindex, endindex, &
!                                 boundfile, time, replace_or_add)
!   subroutine get_vert_species (status, lat, lon, specarr, nparts_old, nparts_new, &
!                                boundfile, time)
!   Subroutine set_coor (lat, lon, vert_array)
!   subroutine interpolate_species (lat, lon, lev, specarr,  &
!                                   lat_low, lon_low, lev_low, specarr_low, &
!                                   startindex, endindex)
!  
!*****************************************************************************!

Module messy_clamsbmix_tools



contains

  !****************************************************************************
  !
  !****************************************************************************
  SUBROUTINE pack_array (array, mask)

    USE messy_clams_global,     ONLY: mdi, prec

    implicit none

    real(prec),  dimension(:), pointer :: array
    logical,     dimension(:), pointer :: mask

    real(prec),  dimension(:), pointer :: helparr
    integer  :: nparts


    nparts = count(mask)

    allocate (helparr(nparts))
    
    helparr = pack (array,mask)
    array = mdi
    array(1:nparts) = helparr

    deallocate (helparr)

  END SUBROUTINE pack_array

  !****************************************************************************
  ! Count number of points between lev_min and lev_max
  !****************************************************************************
  subroutine get_nparts_layer (lev, indices, lev_min, lev_max, counter)

    use messy_clams_global,   ONLY: prec, dnparts
    
    implicit none
    
    REAL(PREC),  DIMENSION(:), POINTER :: lev
    INTEGER,     DIMENSION(:), POINTER :: indices
    REAL(PREC)                         :: lev_min, lev_max
    INTEGER                            :: counter

    integer,     dimension(:), allocatable :: ind
    integer :: i

    allocate (ind(dnparts))
    ind = 0
    counter = 0

    !write (*,*) 'in get_nparts_layer: dnparts=',dnparts

    do i = 1, dnparts
       if (lev_min <= lev(i) .and. lev(i) <= lev_max) then 
          counter = counter+1
          ind(counter) = i
       endif
    enddo
       
    !write (*,*) 'in get_nparts_layer: counter=',counter

    if (counter > 0) then
       allocate (indices(counter))
       indices = ind(1:counter)
    endif
    
    deallocate (ind)

  end subroutine get_nparts_layer

  !****************************************************************************
  !
  !****************************************************************************
  subroutine time_index (btime, time, it, dt)

    USE messy_clams_global,   ONLY: DP, PREC

      implicit none

      real(DP),  dimension(:), intent(in)     :: btime
      real(DP),                intent(in)     :: time
      integer,                 intent(inout)  :: it
      real(PREC),              intent(inout)  :: dt
      
      integer :: ntime, i

      ntime = size(btime)
      if ( btime(1) >= time) then
         it=1
         dt=0.
      elseif ( btime(ntime) <= time) then
         it=ntime-1
         dt=1.
      else 
         do i=1,ntime-1
            it =i
            if (btime(i+1) > time) exit
         enddo
         dt=(time - btime(it))/(btime(it+1) - btime(it))
      endif

  end subroutine time_index

  !****************************************************************************
  !
  !****************************************************************************
  subroutine grid_index (bvec, vec, ix, dx)

    USE messy_clams_global,   ONLY: PREC

    implicit none

    real(PREC), dimension(:), intent(in)     :: bvec, vec
    integer,    dimension(:), intent(inout)  :: ix
    real(PREC), dimension(:), intent(inout)  :: dx
    
    integer   :: i, nbvec, nvec
    
    nbvec=size(bvec)
    nvec=size(vec)
    ix=1
    do i=1,nbvec-1
       where(vec(:) >= bvec(i)) ix=i
    enddo
    do i=1,nvec
       dx(i) =(vec(i) - bvec(ix(i)))/(bvec(ix(i)+1) - bvec(ix(i)))
    enddo
    where(vec(:) > maxval(bvec))
       dx=1.
       ix=nbvec-1
    endwhere
    where(vec(:) < minval(bvec))
       dx=0.
       ix=1
    endwhere
     
  end subroutine grid_index


  !****************************************************************************
  !
  !****************************************************************************
  subroutine interpol_time_grid (status, specarr, speclist, ncid, dim1, dim2, &
                                 startindex, endindex,  &
                                 it, dt, ix, dx, replace_or_add)
 
    use messy_clams_global,     ONLY: prec, species_type, nspec

    use netcdf

    implicit none

    TYPE(species_type), DIMENSION(:), POINTER    :: specarr
    CHARACTER(*),       DIMENSION(:), POINTER    :: speclist
    integer,                       intent(inout) :: status
    integer,                          intent(in) :: ncid, dim1, dim2
    integer,                          intent(in) :: startindex, endindex
    integer,                          intent(in) :: it
    real(PREC),                       intent(in) :: dt
    integer,            dimension(:), intent(in) :: ix
    real(PREC),         dimension(:), intent(in) :: dx
    integer                                      :: replace_or_add

!       character(*)                                :: spec_tracer
!       logical                                     :: tracer_found

    real(PREC),  dimension(:,:),   allocatable :: bdata

    real(PREC)    :: x1, x2
    integer       :: rcode, varid
    integer       :: nspecnames, ispec, i, ind

    status = 0 ! no error

    allocate (bdata(dim1,dim2))

!    tracer_found=.false.

    nspecnames = size(speclist)

    do ispec = 1, nspecnames

       !------------------------------------------------------------
       ! get index of species in specarr => ind
       !------------------------------------------------------------
       ind = -1
       do i = 1, nspec
          if (speclist(ispec)==specarr(i)%name) ind = i
       enddo
       if (ind==-1) then
          write (*,*) 'Cannot redefine ',trim(speclist(ispec)),' !!!'
          write (*,*) trim(speclist(ispec)),' not in SPECARR !!!'
          cycle
       endif

       !------------------------------------------------------------
       ! read species from boundfile
       !------------------------------------------------------------
       rcode = nf90_inq_varid (ncid, trim(speclist(ispec)), varid)
       if (rcode /= 0) then

          ! Fuer Konstanten und nicht in CHEM genutzte Variablen:
          !   kein Abbruch => status=0 !
          if (specarr(ind)%ctype=='CT' .or. specarr(ind)%ctype==' ') then
             write (*,*) 'WARNING: Variable ',trim(specarr(ind)%name), &
                  ' not found on boundfile  !!!'
             CYCLE
          else
             ! Alle in CHEM genutzten Spezies (mit Ausnahme der konstanten Spezies)
             ! muessen im Boundfile enthalten sein 
             write (*,*) 'ERROR: Variable ',trim(specarr(ind)%name), &
                  ' not found on boundfile !!!'
             status = -1
             return
          endif

       else

          write (*,*) 'read variable ',trim(speclist(ispec))
          status = nf90_get_var (ncid, varid, bdata)
          if (status /= 0) then 
             write(*,*) 'Sub. interpol_time_grid: error reading variable '//trim(speclist(ispec))
             return
          endif

       endif

       !------------------------------------------------------------
       ! interpolate on boundaries
       !------------------------------------------------------------
       do i = 1, endindex-startindex+1
          x1=bdata(ix(i),it) + &
               dx(i) * ( bdata(ix(i)+1,it) - bdata(ix(i),it))
          x2=bdata(ix(i),it+1) + &
               dx(i) * ( bdata(ix(i)+1,it+1) - bdata(ix(i),it+1))
          
          if (replace_or_add==2) then
             specarr(ind)%values(startindex+i-1) = &
                  specarr(ind)%values(startindex+i-1) + x1 + dt* (x2 - x1)
          else
             specarr(ind)%values(startindex+i-1) = x1 + dt* (x2 - x1)
          endif
       enddo

       
!        if (spec == spec_tracer) then
!           write (*,*) 'set tracer ',trim(spec_tracer)
!           tracer_found=.true.
!           do i = 1, nparts 
!              ap_s(i)%tracer= ap_s(i)%c(ipos)
!           enddo
!        endif

    enddo

!     if (.not. tracer_found .and. spec_tracer/=' ') then 
!        write (*,*) 'read tracer ',trim(spec_tracer)
!        rcode = nf90_inq_varid(ncid,spec_tracer,varid)
!        if (rcode == 0) then
!           rcode = nf90_get_var(ncid,varid,bdata)
!           if (rcode /= 0) then 
!              write(*,*) 'error reading variable ', trim(spec_tracer)
!              stop
!           endif
!           do i=1,nparts
!              x1=bdata(ix(i),it) + &
!                   dx(i) * ( bdata(ix(i)+1,it) - bdata(ix(i),it))
!              x2=bdata(ix(i),it+1) + &
!                   dx(i) * ( bdata(ix(i)+1,it+1) - bdata(ix(i),it+1))
!              if (replace_or_add==2) then
!                 ap_s(i)%tracer = ap_s(i)%tracer + x1 + dt* (x2 - x1)
!              else
!                 ap_s(i)%tracer = x1 + dt* (x2 - x1)
!              endif
!           enddo
!        endif
!     endif
!     write (*,*) 

     deallocate (bdata)

  end subroutine interpol_time_grid

  !****************************************************************************
  !
  !****************************************************************************
  subroutine interpol_time_grid_3d (status, specarr, speclist, &
                                    ncid, dim1, dim2, dim3, &
                                    startindex, endindex,  &
                                    it, dt, ix, dx, iy, dy, replace_or_add)

    USE messy_clams_global,     ONLY: prec, species_type, nspec

    use netcdf

    implicit none

    TYPE(species_type), DIMENSION(:), POINTER    :: specarr
    CHARACTER(*),       DIMENSION(:), POINTER    :: speclist
    integer,                       intent(inout) :: status
    integer,                          intent(in) :: ncid, dim1, dim2, dim3
    integer,                          intent(in) :: startindex, endindex
    integer,                          intent(in) :: it
    real(PREC),                       intent(in) :: dt
    integer,            dimension(:), intent(in) :: ix, iy
    real(PREC),         dimension(:), intent(in) :: dx, dy
    integer                                      :: replace_or_add


!      character(*)          :: spec_tracer
!      logical               :: tracer_found

    real(PREC),  dimension(:,:,:),   allocatable :: bdata

    real(PREC)   :: x1, x2, y1, y2
    integer      :: rcode, varid
    integer      :: nspecnames, ispec, i, ind

    status = 0 ! no error

    allocate (bdata(dim1,dim2,dim3))

!       tracer_found=.false.

    nspecnames = size(speclist)

    do ispec = 1, nspecnames

       !------------------------------------------------------------
       ! get index of species in specarr => ind
       !------------------------------------------------------------
       ind = -1
       do i = 1, nspec
          if (speclist(ispec)==specarr(i)%name) ind = i
       enddo
       if (ind==-1) then
          write (*,*) 'Cannot redefine ',trim(speclist(ispec)),' !!!'
          write (*,*) trim(speclist(ispec)),' not in SPECARR !!!'
          cycle
       endif

       !------------------------------------------------------------
       ! read species from boundfile
       !------------------------------------------------------------
       rcode = nf90_inq_varid (ncid, trim(speclist(ispec)), varid)
       if (rcode /= 0) then

          ! Fuer Konstanten und nicht in CHEM genutzte Variablen:
          !   kein Abbruch => status=0 !
          if (specarr(ind)%ctype=='CT' .or. specarr(ind)%ctype==' ') then
             write (*,*) 'WARNING: Variable ',trim(specarr(ind)%name), &
                  ' not found on boundfile  !!!'
             CYCLE
          else
             ! Alle in CHEM genutzten Spezies (mit Ausnahme der konstanten Spezies)
             ! muessen im Boundfile enthalten sein 
             write (*,*) 'ERROR: Variable ',trim(specarr(ind)%name), &
                  ' not found on boundfile !!!'
             status = -1
             return
          endif

       else

          write (*,*) 'read variable ',trim(speclist(ispec))
          status = nf90_get_var (ncid, varid, bdata)
          if (status /= 0) then 
             write(*,*) 'Sub. interpol_time_grid_3d: error reading variable '//trim(speclist(ispec))
             return
          endif

       endif

       !------------------------------------------------------------
       ! interpolate on boundaries
       !------------------------------------------------------------
       do i = 1, endindex-startindex+1
             
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
             
          if (replace_or_add == 2) then
             specarr(ind)%values(startindex+i-1)  = &
                  specarr(ind)%values(startindex+i-1) + x1 + dt* (x2 - x1)
          else
             specarr(ind)%values(startindex+i-1) = x1 + dt* (x2 - x1)
          endif

       enddo        
               
!        if (spec == spec_tracer) then
!           tracer_found=.true.
!           do i = 1, nparts 
!              ap_s(i)%tracer= ap_s(i)%c(ipos)
!           enddo
!        endif

    enddo

!     if (.not. tracer_found .and. spec_tracer/=' ') then 
!        rcode = nf90_inq_varid(ncid,spec_tracer,varid)
!        if (rcode == 0) then
!           rcode = nf90_get_var(ncid,varid,bdata)
!           if (rcode /= 0) then 
!              write(*,*) 'error reading variable', spec_tracer
!              stop
!           endif
!           do i=1,nparts
! !              x1=bdata(ix(i),it) + &
! !                   dx(i) * ( bdata(ix(i)+1,it) - bdata(ix(i),it))
! !              x2=bdata(ix(i),it+1) + &
! !                   dx(i) * ( bdata(ix(i)+1,it+1) - bdata(ix(i),it+1))

!              y1 = bdata(ix(i),iy(i),it) + &
!                      dy(i) * (bdata(ix(i),iy(i)+1,it) - bdata(ix(i),iy(i),it))
!              y2 = bdata(ix(i)+1,iy(i),it) + &
!                      dy(i) * (bdata(ix(i)+1,iy(i)+1,it) - bdata(ix(i)+1,iy(i),it))
!              x1  = y1 + dx(i) * (y2 -y1)

!              y1 = bdata(ix(i),iy(i),it+1) + &
!                      dy(i) * (bdata(ix(i),iy(i)+1,it+1) - bdata(ix(i),iy(i),it+1))
!              y2 = bdata(ix(i)+1,iy(i),it+1) + &
!                      dy(i) * (bdata(ix(i)+1,iy(i)+1,it+1) - bdata(ix(i)+1,iy(i),it+1))
!              x2  = y1 + dx(i) * (y2 -y1)

!              if (replace_or_add == 2) then
!                 ap_s(i)%tracer = ap_s(i)%tracer + x1 + dt* (x2 - x1)
!              else
!                 ap_s(i)%tracer = x1 + dt* (x2 - x1)
!              endif
!           enddo
!        endif
!     endif

    deallocate (bdata)

   end subroutine interpol_time_grid_3d


  !****************************************************************************
  !
  !****************************************************************************
  SUBROUTINE get_boundfile_indices (status, boundfile, &
                                    lat_min, lat_max, lev_min, lev_max, &
                                    indices, use_oro)

    USE messy_clams_global,          ONLY: rank, prec, mdi, eps, init_vertcoorname, &
                                           buffersize
    USE messy_clamsmix_global,       ONLY: levelrange, &
                                           l_nlevs, l_min_act, l_max_act
    USE messy_clamsbmix_global,      ONLY: latlon_grid, lev0, delta_lev
    USE messy_clams_tools_utils,     ONLY: uppercase
    USE messy_clams_tools_ncutils,   ONLY: nc_check_error
    
    USE netcdf

    implicit none

    integer                               :: status
    character(*)                          :: boundfile
    real(prec)                            :: lat_min, lat_max, lev_min, lev_max
    integer,     dimension(:),   pointer  :: indices
    logical, optional                     :: use_oro

    real(prec), dimension(:), allocatable :: lats, lons, levs
    integer,    dimension(:), allocatable :: ind
    ! levelind : Index of THETA/ZETA-Values in THETA/ZETA-Grid
    integer,    dimension(:),    pointer  :: levelind ! (1:dnparts)
    integer    :: counter
    integer    :: nparts_bound
    integer    :: i, ilat, ilon
    integer    :: ncid, dimid, varid
    logical    :: orography


    status = 0 ! no error

    orography = .false.
    if (present(use_oro)) orography = use_oro
    if (rank==0) write (*,*) 'Use orography:', orography    


    ! Open boundary file
    if (rank==0) write (*,*) 'Open bound file: ',trim(boundfile)
    status = nf90_open (trim(boundfile), nf90_nowrite, ncid, buffersize)
    call nc_check_error (status,'Error on open file '//trim(boundfile),abort=.false.)
    if (status/=0) return

    if (rank==0) write (*,*) 'Vertical coordinate in boundfile: ',trim (init_vertcoorname)


    ! Read dimension NPARTS
    status = nf90_inq_dimid (ncid, 'NPARTS', dimid)    
    call nc_check_error (status,'Cannot find dimension NPARTS',abort=.false.)
    if (status/=0) return
    status = nf90_inquire_dimension (ncid, dimid, len=nparts_bound)
    call nc_check_error (status,'Cannot read dimension NPARTS',abort=.false.)
    if (status/=0) return
    if (rank==0) write (*,*) 'NPARTS in boundfile=',nparts_bound
    
    allocate (ind(nparts_bound))
    allocate (lats(nparts_bound))
    allocate (lons(nparts_bound))
    allocate (levs(nparts_bound))
    allocate (levelind(nparts_bound))

    ! Read lats, lons and levs from boundfile
    status = nf90_inq_varid (ncid, 'LAT', varid)
    call nc_check_error (status,'Cannot find variable LAT',abort=.false.)
    if (status/=0) return
    status = nf90_get_var (ncid, varid, lats)
    call nc_check_error (status,'Cannot read variable LAT',abort=.false.)
    if (status/=0) return
    status = nf90_inq_varid (ncid, 'LON', varid)
    call nc_check_error (status,'Cannot find variable LON',abort=.false.)
    if (status/=0) return
    status = nf90_get_var (ncid, varid, lons)
    call nc_check_error (status,'Cannot read variable LON',abort=.false.)
    if (status/=0) return
    status = nf90_inq_varid (ncid, trim(uppercase(init_vertcoorname)), varid)
    call nc_check_error (status, &
         'Cannot find variable '//trim(uppercase(init_vertcoorname)),abort=.false.)
    if (status/=0) return
    status = nf90_get_var (ncid, varid, levs)
    call nc_check_error (status,&
         'Cannot read variable '//trim(uppercase(init_vertcoorname)),abort=.false.)
    if (status/=0) return


!!!!! => levelrange aus mix wird benoetigt
!!!!! => l_min_act, l_max_act, l_nlevs aus mix wird benoetigt !?!

    if (rank==0) write (*,*) 'l_nlevs=',l_nlevs


    ! for all levels: get level index
    do i = 1, l_nlevs
       where (l_min_act(i)<=levs .and. levs<l_max_act(i)) 
          levelind = i
       end where
    enddo
    where (levs<l_min_act(1))
       levelind = 1
    endwhere
    where (levs>=l_max_act(l_nlevs))
       levelind = l_nlevs
    endwhere
    ! levelind fuer mdi setzen 
    where (abs((levs-mdi)/mdi)<eps)
       levelind = 9999
    endwhere

 
!!!!! Indizes ermittlen, bei denen  lev_min<=lev<=lev_max und lat_min<=lat<=lat_max
!!!!! und Level innerhalb des vom aktuellen rank zu bearbeitenden Bereiches:
!!!!! levelrange(1,rank)<=levelind<=levelrange(2,rank)

    ind = 0
    counter=0

    if (.not. orography) then    ! without orography
       
       ! Get all points between lev_down and lev_in_down 
       ! and between lat_in_down and lat_in_up from boundary file.
       ! Get only points within the levelrange for the current rank!

        do i = 1, nparts_bound
          if ((lev_min <= levs(i) .and. levs(i) <= lev_max) .and. &
              (lat_min <= lats(i) .and. lats(i) <= lat_max) .and. &
              (levelrange(1,rank)<=levelind(i) .and. levelind(i)<=levelrange(2,rank)) ) then 
             counter=counter+1
             ind(counter)=i
          endif
       enddo  

    else   ! with orography

       ! Get all points between lev0 and lev0+delta_lev
       ! and between lat_in_down and lat_in_up from boundary file.
       ! Get only points within the levelrange for the current rank!

       do i = 1, nparts_bound

          ! get position in grid (lat,lon)

          ilat = int((lats(i)-latlon_grid%lat0)/latlon_grid%dlat) + 1
! op_pj_20160606+
!!$          ilon = int((modulo(lons(i),360.)-latlon_grid%lon0)/latlon_grid%dlon) + 1
          ilon = int((modulo(lons(i),360._prec)-latlon_grid%lon0)/latlon_grid%dlon) + 1
! op_pj_20160606-

          if ((lev_min <= levs(i) .and. levs(i) <= lev_max) .and. &
              (lat_min <= lats(i) .and. lats(i) <= lat_max) .and. &
              (lev0(ilon,ilat)<=levs(i) .and. levs(i)<=lev0(ilon,ilat)+delta_lev) .and. &
              (levelrange(1,rank)<=levelind(i) .and. levelind(i)<=levelrange(2,rank)) ) then 
             counter=counter+1
             ind(counter)=i
          endif

       enddo

    endif

    if (rank==0) &
         write (*,*)' Bounds:   ',lev_min,lev_max,lat_min,lat_max

!!$    write (*,*) 'rank, levelrange:',rank,levelrange(1,rank),levelrange(2,rank)
!!$    write (*,*) 'rank, counter:   ', rank, counter

    allocate (indices(counter))
    indices = ind(1:counter)


    ! clean up
    deallocate (lats, lons, levs)
    deallocate (levelind, ind)

    ! close boundary file
    status = nf90_close (ncid)

  END SUBROUTINE get_boundfile_indices


  !****************************************************************************
  !
  !****************************************************************************
  SUBROUTINE read_boundfile (status, boundfile, indices, lat, lon, lev, specarr)

    USE messy_clams_global,          ONLY: rank, prec, dnparts, mdi, buffersize, &
                                           species_type, nspec, init_vertcoorname
    USE messy_clams_tools_utils,     ONLY: uppercase
    USE messy_clams_tools_ncutils,   ONLY: nc_check_error

    use netcdf

    implicit none

    integer                                   :: status
    character(*)                              :: boundfile
    integer,            dimension(:), pointer :: indices
    REAL(PREC),         DIMENSION(:), POINTER :: lat, lon, lev
    TYPE(species_type), DIMENSION(:), POINTER :: specarr


    real(prec), dimension(:), pointer :: helparr
    integer :: nparts_all, dnparts_new
    integer :: i
    integer :: ncid, dimid, varid, rcode

    status = 0

    ! Open boundary file
    if (rank==0) write (*,*) 'Open bound file: ',trim(boundfile)
    status = nf90_open (trim(boundfile), nf90_nowrite, ncid, buffersize)
    call nc_check_error (status,'Error on open file '//trim(boundfile),abort=.false.)
    if (status/=0) return

    ! Read dimension NPARTS -> nparts_all
    status = nf90_inq_dimid (ncid, 'NPARTS', dimid)    
    call nc_check_error (status,'Cannot find dimension NPARTS',abort=.false.)
    if (status/=0) return
    status = nf90_inquire_dimension (ncid, dimid, len=nparts_all)
    call nc_check_error (status,'Cannot read dimension NPARTS',abort=.false.)
    if (status/=0) return
    if (rank==0) write (*,*) 'NPARTS in boundfile=',nparts_all
 
    ! Set nparts
    dnparts_new = size(indices)

    allocate (helparr(nparts_all))

    ! Haenge eingelesenen Punkte an bisherige Felder an: 
    ! array(dnparts+1:dnparts+dnparts_new)

    ! Read LAT
    status = nf90_inq_varid (ncid, 'LAT', varid)
    call nc_check_error (status,'Cannot find variable LAT',abort=.false.)
    if (status/=0) return
    status = nf90_get_var (ncid, varid, helparr)
    call nc_check_error (status,'Cannot read variable LAT',abort=.false.)
    if (status/=0) return
    lat(dnparts+1:dnparts+dnparts_new) = helparr(indices)

    ! Read LON
    status = nf90_inq_varid (ncid, 'LON', varid)
    call nc_check_error (status,'Cannot find variable LON',abort=.false.)
    if (status/=0) return
    status = nf90_get_var (ncid, varid, helparr)
    call nc_check_error (status,'Cannot read variable LON',abort=.false.)
    if (status/=0) return
    lon(dnparts+1:dnparts+dnparts_new) = helparr(indices)

    ! Read THETA/ZETA
    status = nf90_inq_varid (ncid, trim(uppercase(init_vertcoorname)), varid)
    call nc_check_error (status, &
         'Cannot find variable '//trim(uppercase(init_vertcoorname)),abort=.false.)
    if (status/=0) return
    status = nf90_get_var (ncid, varid, helparr)
    call nc_check_error (status,&
         'Cannot read variable '//trim(uppercase(init_vertcoorname)),abort=.false.)
    if (status/=0) return
    lev(dnparts+1:dnparts+dnparts_new) = helparr(indices)

    ! Read chemical species
    do i = 1, nspec
       rcode = nf90_inq_varid (ncid, trim(specarr(i)%name), varid)
       if (rcode /= nf90_noerr) then

          ! Fuer Konstanten und nicht in CHEM genutzte Variablen:
          !   kein Abbruch => status=0 !
          if (specarr(i)%ctype=='CT' .or. specarr(i)%ctype==' ') then
             write (*,*) 'WARNING: Variable ',trim(specarr(i)%name), &
                  ' not found on initfile => Boundaries are set to MDI !!!'
             specarr(i)%values(dnparts+1:dnparts+dnparts_new) = MDI
             CYCLE
          else
             ! Alle in CHEM genutzten Spezies (mit Ausnahme der konstanten Spezies)
             ! muessen im Boundfile enthalten sein 
             write (*,*) 'ERROR: Variable ',trim(specarr(i)%name), &
                  ' not found on initfile !!!'
             status = -1
             return
          endif
       endif

       status = nf90_get_var (ncid, varid, helparr)
       call nc_check_error (status, &
            'Cannot read variable '//trim(specarr(i)%name),abort=.false.)
       if (status/=0) return
       specarr(i)%values(dnparts+1:dnparts+dnparts_new) = helparr(indices)     
    end do

!!!!! DNPARTS CHANGED !!!
    ! Set new number of points on rank
    dnparts = dnparts + dnparts_new

    ! clean up
    deallocate (helparr)

    ! close boundary file
    status = nf90_close (ncid)

  END SUBROUTINE read_boundfile

  !****************************************************************************
  ! redefine species from the boundary file
  !****************************************************************************
  subroutine get_bound_species (status, lat, lon, eqlat, specarr, speclist,  &
                                startindex, endindex, &
                                boundfile, time, replace_or_add)

    use messy_clams_global,        ONLY: PREC, DP, rank, species_type, buffersize
    USE messy_clams_tools_ncutils, ONLY: nc_check_error

    use netcdf
    
    implicit none
    
    REAL(PREC),         DIMENSION(:), POINTER :: lat, lon, eqlat
    TYPE(species_type), DIMENSION(:), POINTER :: specarr
    CHARACTER(*),       DIMENSION(:), POINTER :: speclist
    integer           :: status, startindex, endindex
    character(*)      :: boundfile
    real(DP)          :: time
    integer, optional :: replace_or_add
    
    integer       :: ncid, rcode, dimid, varid
    integer       :: nlats, nlons, ntimes
    integer       :: it
    real(PREC)    :: dt
    logical       :: three_dim
    
    integer,     dimension(:),     allocatable :: ix, iy
    real(PREC),  dimension(:),     allocatable :: blats, blons, dx, dy
    real(DP),    dimension(:),     allocatable :: btime
 
    integer       :: replace_add = 1

    ! character(*)  :: spec_tracer

    status = 0 ! no error

    replace_add = 1
    if (present(replace_or_add)) replace_add = replace_or_add

    ! open bound file
    write(*,*) 'reading '//trim(boundfile)
    rcode = nf90_open (boundfile,nf90_nowrite,ncid, buffersize)
    if (rcode /= 0) then 
       write(*,*)
       write(*,*) 'Error opening file '//trim(boundfile)
       write(*,*) 'Continue without update from file ', trim(boundfile)
       write(*,*)
       return   ! kein Abbruch => status=0
    endif

    ! read dimensions
    status = nf90_inq_dimid (ncid,'lat',dimid)
    call nc_check_error (status,'Cannot find dimension lat',abort=.false.)
    if (status/=0) return
    status = nf90_inquire_dimension (ncid,dimid,len=nlats)
    call nc_check_error (status,'Cannot read dimension lat',abort=.false.)
    if (status/=0) return
    status = nf90_inq_dimid (ncid,'time',dimid)
    call nc_check_error (status,'Cannot find dimension time',abort=.false.)
    if (status/=0) return
    status = nf90_inquire_dimension (ncid,dimid,len=ntimes)
    call nc_check_error (status,'Cannot read dimension time',abort=.false.)
    if (status/=0) return

    three_dim = .false.
    rcode = nf90_inq_dimid (ncid,'lon',dimid) ! nicht status setzten !
    if (rcode==nf90_noerr) then
       rcode = nf90_inquire_dimension (ncid,dimid,len=nlons)
       if (rcode==nf90_noerr) three_dim = .true.
    endif
       
    if (.not. three_dim) then

       write (*,*) 'Variablen im Bound-File sind 2-dimensional !'

       ! allocate arrays
       allocate(blats(nlats),btime(ntimes))
       allocate(dx(endindex-startindex+1),ix(endindex-startindex+1))

       ! get coordinate variables
       status = nf90_inq_varid (ncid,'lat',varid)
       call nc_check_error (status,'Cannot find variable lat',abort=.false.)
       if (status/=0) return
       status = nf90_get_var (ncid,varid,blats)
       call nc_check_error (status,'Cannot read variable lat',abort=.false.)
       if (status/=0) return
       status = nf90_inq_varid (ncid,'time',varid)
       call nc_check_error (status,'Cannot find variable time',abort=.false.)
       if (status/=0) return
       status = nf90_get_var (ncid,varid,btime)
       call nc_check_error (status,'Cannot read variable time',abort=.false.)
       if (status/=0) return
      
       ! get index of time
       call time_index (btime, time, it, dt)
       
       ! get indices of latitudes (use EQLAT)
       call grid_index (blats, eqlat(startindex:endindex), ix, dx)
     
       ! interpolate all species
       call interpol_time_grid (status, specarr, speclist, &
                                ncid, nlats, ntimes, &
                                startindex, endindex, &
                                it, dt, ix, dx, replace_add)
       if (status /= 0) then
          write (*,*) 'Sub. get_bound_species: error in interpol_time_grid!!!'
          return
       endif
 
       ! deallocate arrays
       deallocate(blats, btime, dx, ix)

    else

       write (*,*) 'Variablen im Bound-File sind 3-dimensional !'

       ! allocate arrays
       allocate(blats(nlats),blons(nlons),btime(ntimes))
       allocate(dx(endindex-startindex+1),ix(endindex-startindex+1))
       allocate(dy(endindex-startindex+1),iy(endindex-startindex+1))

       ! get coordinate variables
       status = nf90_inq_varid (ncid,'lat',varid)
       call nc_check_error (status,'Cannot find variable lat',abort=.false.)
       if (status/=0) return
       status = nf90_get_var (ncid,varid,blats)
       call nc_check_error (status,'Cannot read variable lat',abort=.false.)
       if (status/=0) return
       status = nf90_inq_varid (ncid,'lon',varid)
       call nc_check_error (status,'Cannot find variable lon',abort=.false.)
       if (status/=0) return
       status = nf90_get_var (ncid,varid,blons)
       call nc_check_error (status,'Cannot read variable lon',abort=.false.)
       if (status/=0) return
       status = nf90_inq_varid (ncid,'time',varid)
       call nc_check_error (status,'Cannot find variable time',abort=.false.)
       if (status/=0) return
       status = nf90_get_var (ncid,varid,btime)
       call nc_check_error (status,'Cannot read variable time',abort=.false.)
       if (status/=0) return
      
       ! get index of time
       call time_index (btime, time, it, dt)

       ! get indices of latitudes
       call grid_index (blats, lat(startindex:endindex), ix, dx)
       
       ! get indices of longitudes 
       call grid_index (blons, lon(startindex:endindex), iy, dy)

       call interpol_time_grid_3d (status, specarr, speclist, &
                                 ncid, nlats, nlons, ntimes, &
                                 startindex, endindex,  &
                                 it, dt, ix, dx, iy, dy, replace_add)      
       if (status /= 0) then
          write (*,*) 'Sub. get_bound_species: error in interpol_time_grid_3d!!!'
          return
       endif

       ! deallocate arrays
       deallocate (blats, blons, btime, dx, ix, dy, iy)

    endif

    rcode = nf90_close(ncid)

  end subroutine get_bound_species

  !****************************************************************************
  !
  !****************************************************************************
  subroutine get_vert_species (status, lev, lon, specarr, speclist, &
                               startindex, endindex, &
                               boundfile, time)

    use messy_clams_global,        ONLY: PREC, DP, species_type, init_vertcoorname, &
                                         buffersize
    USE messy_clams_tools_utils,   ONLY: lowercase
    USE messy_clams_tools_ncutils, ONLY: nc_check_error

    use netcdf
    
    implicit none
    
    REAL(PREC),         DIMENSION(:), POINTER :: lev, lon
    TYPE(species_type), DIMENSION(:), POINTER :: specarr
    CHARACTER(*),       DIMENSION(:), POINTER :: speclist
    integer       :: status, startindex, endindex
    character(*)  :: boundfile
    real(DP)      :: time


    integer       :: ncid, rcode, dimid, varid
    integer       :: nlevs, ntimes, nlons
    integer       :: it
    real(PREC)    :: dt
    logical       :: three_dim

    integer,    dimension(:),  allocatable :: ix, iy
    real(PREC), dimension(:),  allocatable :: blons, blevs, dx, dy
    real(DP),   dimension(:),  allocatable :: btime

    status = 0 ! no error

    rcode = nf90_open (boundfile,nf90_nowrite,ncid, buffersize)
    if (rcode /= 0) then 
       write(*,*)
       write(*,*) 'Error opening file '//trim(boundfile)
       write(*,*) 'Continue without update from file ', trim(boundfile)
       write(*,*)
       return ! kein Abbruch => status=0
    endif
    write(*,*) 'reading '//trim(boundfile)

    ! read dimensions
    status = nf90_inq_dimid (ncid,trim(lowercase(init_vertcoorname)),dimid)
    call nc_check_error (status,&
         'Cannot find dimension '//trim(lowercase(init_vertcoorname)),abort=.false.)
    if (status/=0) return
    status = nf90_inquire_dimension( ncid,dimid,len=nlevs)
    call nc_check_error (status,&
         'Cannot read dimension '//trim(lowercase(init_vertcoorname)),abort=.false.)
    if (status/=0) return
    status = nf90_inq_dimid (ncid,'time',dimid)
    call nc_check_error (status,'Cannot find dimension time',abort=.false.)
    if (status/=0) return
    status = nf90_inquire_dimension (ncid,dimid,len=ntimes)
    call nc_check_error (status,'Cannot find dimension time',abort=.false.)
    if (status/=0) return

    three_dim = .false.
    rcode = nf90_inq_dimid (ncid,'lon',dimid) ! nicht status setzten !
    if (rcode==nf90_noerr) then
       rcode = nf90_inquire_dimension (ncid,dimid,len=nlons)
       if (rcode==nf90_noerr) three_dim = .true.
    endif

    if (.not. three_dim) then

       write (*,*) 'Variablen im Bound-File sind 2-dimensional !'

       ! allocate arrays
       allocate(blevs(nlevs), btime(ntimes))
       allocate(dx(endindex-startindex+1), ix(endindex-startindex+1))
       
       ! get coordinate variables
       status = nf90_inq_varid(ncid,trim(lowercase(init_vertcoorname)),varid)
       call nc_check_error (status, &
            'Cannot find variable '//trim(lowercase(init_vertcoorname)),abort=.false.)
       if (status/=0) return
       status = nf90_get_var(ncid,varid,blevs)
       call nc_check_error (status, &
            'Cannot read variable '//trim(lowercase(init_vertcoorname)),abort=.false.)
       if (status/=0) return
       status = nf90_inq_varid(ncid,'time',varid)
       call nc_check_error (status,'Cannot find variable time',abort=.false.)
       if (status/=0) return
       status = nf90_get_var(ncid,varid,btime)
       call nc_check_error (status,'Cannot read variable time',abort=.false.)
       if (status/=0) return
       
       ! get index of time
       call time_index (btime, time, it, dt)
       
       ! get indices of levels
       call grid_index (blevs, lev(startindex:endindex), ix, dx)
       
       ! use next neighbor vertical interpolation
       ! in order to avoid numerical mixing 
       ! induced by interpolation
       ! where (dx .gt. 0.5) 
       !    dx=1.0
       ! elsewhere
       !    dx=0.0
       ! endwhere

       call interpol_time_grid (status, specarr, speclist, &
                                ncid, nlevs, ntimes, &
                                startindex, endindex, &
                                it, dt, ix, dx, 1)
       if (status /= 0) then
          write (*,*) 'Sub. get_vert_species: error in interpol_time_grid!!!'
          return
       endif
    
       ! deallocate arrays
       deallocate (blevs, btime, dx, ix)

    else

       write (*,*) 'Variablen im Bound-File sind 3-dimensional !'

       ! allocate arrays
       allocate (blons(nlons),blevs(nlevs),btime(ntimes))
       allocate (dx(endindex-startindex+1),ix(endindex-startindex+1))
       allocate (dy(endindex-startindex+1),iy(endindex-startindex+1))

       ! get coordinate variables
       status = nf90_inq_varid(ncid,'lon',varid)
       call nc_check_error (status,'Cannot find variable lon',abort=.false.)
       if (status/=0) return
       status = nf90_get_var(ncid,varid,blons)
       call nc_check_error (status,'Cannot read variable lon',abort=.false.)
       if (status/=0) return
       status = nf90_inq_varid(ncid,trim(lowercase(init_vertcoorname)),varid)
       call nc_check_error (status, &
            'Cannot find variable '//trim(lowercase(init_vertcoorname)),abort=.false.)
       if (status/=0) return
       status = nf90_get_var(ncid,varid,blevs)
       call nc_check_error (status, &
            'Cannot read variable '//trim(lowercase(init_vertcoorname)),abort=.false.)
       if (status/=0) return
       status = nf90_inq_varid(ncid,'time',varid)
       call nc_check_error (status,'Cannot find variable time',abort=.false.)
       if (status/=0) return
       status = nf90_get_var(ncid,varid,btime)
       call nc_check_error (status,'Cannot read variable time',abort=.false.)
       if (status/=0) return
       
       ! get index of time
       call time_index (btime, time, it, dt)
       
       ! get indices of longitudes 
       call grid_index (blons, lon(startindex:endindex), ix, dx)
       
       ! get indices of levels
       call grid_index (blevs, lev(startindex:endindex), iy, dy)

!!!!! ???
       ! use next neighbor vertical interpolation
       ! in order to avoid numerical mixing 
       ! induced by interpolation
!        where (dx .gt. 0.5) 
!           dx=1.0
!        elsewhere
!           dx=0.0
!        endwhere

       call interpol_time_grid_3d (status, specarr, speclist, &
                                   ncid, nlons, nlevs, ntimes, &
                                   startindex, endindex,  &
                                   it, dt, ix, dx, iy, dy, 1)
       if (status /= 0) then
          write (*,*) 'Sub. get_vert_species: error in interpol_time_grid_3d!!!'
          return
       endif


       ! deallocate arrays
       deallocate (blons, blevs, btime, dx, ix, dy, iy)

    endif

    rcode = nf90_close(ncid)

  end subroutine get_vert_species

  !****************************************************************************
  ! Transform ap-coordinates (lon, lat) on kart. coord. on a unit sphere
  !****************************************************************************
  Subroutine set_coor (lat, lon, vert_array)

    USE messy_clams_global,      ONLY: PREC
    
    implicit none
    
    REAL(PREC), DIMENSION(:),   INTENT(IN) :: lat, lon
    REAL(PREC), DIMENSION(:,:), POINTER    :: vert_array
 
    real  :: pi
    real, dimension(:), allocatable :: lat_rad, lon_rad

    pi=4.*atan(1.)

    allocate (lon_rad (size(lon)))
    allocate (lat_rad (size(lat)))

    lon_rad = lon * pi / 180.          ! transform from deg to rad
    lat_rad = (90. - lat) * pi / 180.  ! use spherical coordinates
                                       !   0. < lons < 2*PI
                                       !   0. < lats < PI
    vert_array(1,:) = cos(lon_rad)*sin(lat_rad) ! transf from spheric
    vert_array(2,:) = sin(lon_rad)*sin(lat_rad) ! coord. to  kart.
    vert_array(3,:) = cos(lat_rad)  ! coord on a unit sphere (r=1)
    
    deallocate (lon_rad, lat_rad)

  end Subroutine set_coor

  !****************************************************************************
  ! interpolate species from init file
  !****************************************************************************
  subroutine interpolate_species (status, lat, lon, lev, specarr,  &
                                  lat_init, lon_init, lev_init, specarr_init, &
                                  startindex, endindex)

    USE messy_clams_global,         ONLY: PREC, mdi, eps, species_type, nspec
    USE messy_clamsbmix_global,     ONLY: lev_down, lev_in_down, delta_lev, &
                                          max_dist
    USE messy_clamsmix_global,      ONLY: adapt_par
    USE messy_clams_tools_triang,   ONLY: triang_qhull, find_triangle
    USE messy_clams_tools_interpol, ONLY: determ_weight

    implicit none

    integer                                   :: status
    integer                                   :: startindex, endindex
    REAL(PREC),         DIMENSION(:), POINTER :: lat, lon, lev, &
                                                 lat_init, lon_init, lev_init
    TYPE(species_type), DIMENSION(:), POINTER :: specarr, specarr_init


    REAL,          PARAMETER               :: earth_radius=6371.
    INTEGER,       PARAMETER               :: nb_max = 20 

    REAL(PREC),    DIMENSION(:,:), POINTER :: coor
    REAL(PREC),    DIMENSION(:,:), POINTER :: coor_init
    INTEGER,       DIMENSION(:),   POINTER :: nb_init
    INTEGER,       DIMENSION(:,:), POINTER :: ind_init
    INTEGER,       DIMENSION(:,:), POINTER :: triangles
    REAL(PREC),    DIMENSION(3,3)          :: EckCoor
    REAL(PREC),    DIMENSION(3)            :: lev_arr
    REAL(PREC),    DIMENSION(3)            :: weight
    INTEGER,       DIMENSION(3)            :: index_init
    REAL(PREC)                             :: lev_value
    REAL(PREC)                             :: dlev, ratio
    INTEGER                                :: ntriang, nparts, nparts_init, &
                                              t_act, ipart, ispec, i, n_mdi, &
                                              triangle_index
    LOGICAL                                :: tracer_found

    status = 0 ! no error

    nparts = endindex-startindex+1
    nparts_init = size(lat_init)

    ! allocate arrays
    allocate (coor    (3,nparts))
    allocate (coor_init(3,nparts_init))
    allocate (nb_init  (nparts_init))
    allocate (ind_init (nb_max,nparts_init))
    
    ! Transform ap-coordinates (lon, lat) to kart. coord. on a unit sphere
    call set_coor (lat(startindex:endindex), lon(startindex:endindex), coor)
    call set_coor (lat_init, lon_init, coor_init)
    
    ! Carry out the Delaunay triangulation (using qhull) 
    call triang_qhull  (status, coor_init, nb_init, ind_init, triangles, ntriang)
    if (status/=0) then
       write (*,*) 'Error in triang_qhull !!!'
       return
    endif

    t_act = int((ntriang-1)/2)
    tracer_found = .false.

!!!!! nur fuer TEST     
!    write (*,*) 'max_dist=',max_dist
    
    do ipart = 1, nparts

       ! get index of triangle
       triangle_index = find_triangle (t_act, coor_init, nb_init, ind_init, &
                                       triangles, ntriang, coor(:,ipart), &
                                       write_warnings=.false., max_dist=max_dist)

       ! set missing value if triangle not found
       if (triangle_index<0) then
          do ispec = 1, nspec
             specarr(ispec)%values(startindex+ipart-1) = mdi
!              spec=trim(spec_arr(ispec))
!              if (spec == spec_tracer) then
!                 tracer_found=.true.
!                 ap_s_1(ipart)%tracer = mdi
!              endif
          enddo
          cycle
       endif

       ! get vertices of triangle
       index_init = triangles(:,triangle_index)
       do i=1,3
          EckCoor(:,i) = coor_init(:,index_init(i))
       enddo
     
       ! Umspeichern auf lev_arr fuehrt zu drastischer Laufzeitverbesserung !
       lev_arr = lev_init(index_init)

       ! get weights 

       if (lat(startindex+ipart-1)<adapt_par%lat_down .or. &
           lat(startindex+ipart-1)>adapt_par%lat_up) then
          ratio = adapt_par%r_mean_c/earth_radius
       else
          ratio = adapt_par%r_mean_h/earth_radius
       end if

       if (delta_lev > 0) then  
          dlev = delta_lev
       else
          dlev = lev_in_down - lev_down
       endif

       call determ_weight (EckCoor,lev_arr,&
                          coor(:,ipart),lev(startindex+ipart-1),&
                          weight,lev_value,dlev,ratio,2) 

       do ispec=1, nspec

!           spec=trim(spec_arr(ispec))
 
          ! wenn eines der Gewichte 1. ist:
          if (count(abs(1.-weight(:))<eps) > 0) then
             if (abs(weight(1)-1.)<eps) then
                specarr(ispec)%values(startindex+ipart-1) = specarr_init(ispec)%values(index_init(1))
             elseif (abs(weight(2)-1.)<eps) then
                specarr(ispec)%values(startindex+ipart-1) = specarr_init(ispec)%values(index_init(2))
             else
                specarr(ispec)%values(startindex+ipart-1) = specarr_init(ispec)%values(index_init(3))
             endif
             
          else
             
             ! Ueberpruefe, ob einer der Datenwerte MDI ist

!!!!! Probleme mit NAG und g95 !!!
!             SELECT CASE (count(abs(mdi-specarr_init(ispec)%values(index_init(:)))<eps))

             n_mdi = 0
             if (abs(mdi-specarr_init(ispec)%values(index_init(1)))<eps) n_mdi = n_mdi+1
             if (abs(mdi-specarr_init(ispec)%values(index_init(2)))<eps) n_mdi = n_mdi+1
             if (abs(mdi-specarr_init(ispec)%values(index_init(3)))<eps) n_mdi = n_mdi+1
             SELECT CASE (n_mdi)

             CASE(0) ! kein missing_value vorhanden:
!!!!! Probleme mit NAG und g95 !!!
!               specarr(ispec)%values(startindex+ipart-1) = &
!                      sum(weight(:)*specarr_init(ispec)%values(index(:)))
                specarr(ispec)%values(startindex+ipart-1) = &
                                         weight(1)*specarr_init(ispec)%values(index_init(1)) + &
                                         weight(2)*specarr_init(ispec)%values(index_init(2)) + &
                                         weight(3)*specarr_init(ispec)%values(index_init(3)) 

             CASE(1) ! ein missing_value: Gewicht wird auf die anderen beiden verteilt
                if (abs(mdi-specarr_init(ispec)%values(index_init(1)))<eps) then
                   specarr(ispec)%values(startindex+ipart-1) = &
                        specarr_init(ispec)%values(index_init(2)) * &
                        (weight(2)+weight(1)*weight(2)/(weight(2)+weight(3)))+ &
                        specarr_init(ispec)%values(index_init(3)) * &
                        (weight(3)+weight(1)*weight(3)/(weight(2)+weight(3)))
                elseif (abs(mdi-specarr_init(ispec)%values(index_init(2)))<eps) then
                   specarr(ispec)%values(startindex+ipart-1) = &
                        specarr_init(ispec)%values(index_init(1)) * &
                        (weight(1)+weight(2)*weight(1)/(weight(1)+weight(3)))+ &
                        specarr_init(ispec)%values(index_init(3)) * &
                        (weight(3)+weight(2)*weight(3)/(weight(1)+weight(3)))
                else
                   specarr(ispec)%values(startindex+ipart-1) = &
                        specarr_init(ispec)%values(index_init(1)) * &
                        (weight(1)+weight(3)*weight(1)/(weight(1)+weight(2)))+ &
                        specarr_init(ispec)%values(index_init(2)) * &
                        (weight(2)+weight(3)*weight(2)/(weight(1)+weight(2)))
                endif
                
             CASE DEFAULT ! mehr als ein missing_value:
                specarr(ispec)%values(startindex+ipart-1) = mdi
             END SELECT

          endif
         
!           if (spec == spec_tracer) then
!              tracer_found=.true.
!              ap_s_1(ipart)%tracer= ap_s_1(ipart)%c(ispec)
!           endif

       end do  ! do ispec=1, nspec

    end do  ! do ipart = 1, nparts
  
    ! deallocate arrays
    deallocate (coor)
    deallocate (coor_init, nb_init, ind_init)
    deallocate (triangles)

!     if (.not. tracer_found) then
!        do ipart = 1, nparts
!            ap_s_1(ipart)%tracer=0.
!        enddo
!     endif


  end subroutine interpolate_species



End Module messy_clamsbmix_tools

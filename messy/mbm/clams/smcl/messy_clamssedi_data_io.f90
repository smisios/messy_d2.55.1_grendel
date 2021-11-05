!******************************************************************************
!
! Modul clamssedi_data_io
!
! Written by Gebhard Guenther, Thomas Breuer, Nicole Thomas
! IEK-7, Forschungszentrum Juelich
! Last Modified By: N.Thomas
! Last Modified On: Mon Jan 11 09:07:30 2016
! 
! subroutine get_init_grid (status)
! subroutine get_nairparcels_init (status)
! subroutine read_nucleation_table
!
!****************************************************************************** 

Module messy_clamssedi_data_io

CONTAINS

  !*********************************************************************
  ! Read initial particles from an old clams sedi output
  !*********************************************************************
  subroutine read_clams_particles (status)

    USE netcdf

    use messy_clams_global, only: rank, buffersize
    use messy_clamssedi_global, only: part_init_file, particles, &
                                      nparticles, part_id_max
    USE messy_clams_tools_ncutils, ONLY: nc_check_error
    
    implicit none

    integer       :: status
    integer       :: ncid, dimid, varid
    character(5)       :: molten

    ! open init file
    status = nf90_open(part_init_file, nf90_nowrite, ncid, buffersize)
    CALL nc_check_error(status, 'Cannot open file '//TRIM(part_init_file), abort=.false.)
    if (status/=0) return

    ! check if initial particles file is empty
    status = nf90_get_att (ncid, nf90_global, 'all_molten', molten)
    call nc_check_error (status,'Cannot find attribute all_molten',abort=.false.)
    if (status/=0) return
    
    if(molten == 'true') then
       if (rank == 0) write(*,*) 'part_init_file is empty!!'
       if (rank == 0) write(*,*) 'continue with zero initial particles!'
       return
    endif

    ! determine dimension: number of particles (naero)
    status = nf90_inq_dimid (ncid, "naero", dimid)
    call nc_check_error (status, "Cannot find dimension naero", abort=.false.)
    if (status/=0) return
    status = nf90_inquire_dimension(ncid, dimid, len=nparticles)
    call nc_check_error (status, "Cannot read dimension naero", abort=.false.)
    if (status/=0) return

    ! read and set particles channels
    status = nf90_inq_varid(ncid, "LAT", varid)
    call nc_check_error (status, "Cannot find variable LAT", abort=.false.)
    if (status/=0) return
    status = nf90_get_var(ncid, varid, particles%lat(1:nparticles))
    call nc_check_error (status, "Cannot read variable LAT", abort=.false.)
    if (status/=0) return

    status = nf90_inq_varid(ncid, "LON", varid)
    call nc_check_error (status, "Cannot find variable LON", abort=.false.)
    if (status/=0) return
    status = nf90_get_var(ncid, varid, particles%lon(1:nparticles))
    call nc_check_error (status, "Cannot read variable LON", abort=.false.)
    if (status/=0) return

    status = nf90_inq_varid(ncid, "ZETA", varid)
    call nc_check_error (status, "Cannot find variable ZETA", abort=.false.)
    if (status/=0) return
    status = nf90_get_var(ncid, varid, particles%lev(1:nparticles))
    call nc_check_error (status, "Cannot read variable ZETA", abort=.false.)
    if (status/=0) return

    status = nf90_inq_varid(ncid, "PRESS", varid)
    call nc_check_error (status, "Cannot find variable PRESS", abort=.false.)
    if (status/=0) return
    status = nf90_get_var(ncid, varid, particles%pressure(1:nparticles))
    call nc_check_error (status, "Cannot read variable PRESS", abort=.false.)
    if (status/=0) return

    status = nf90_inq_varid(ncid, "TEMP", varid)
    call nc_check_error (status, "Cannot find variable TEMP", abort=.false.)
    if (status/=0) return
    status = nf90_get_var(ncid, varid, particles%temperature(1:nparticles))
    call nc_check_error (status, "Cannot read variable TEMP", abort=.false.)
    if (status/=0) return

    status = nf90_inq_varid(ncid, "LON", varid)
    call nc_check_error (status, "Cannot find variable LON", abort=.false.)
    if (status/=0) return
    status = nf90_get_var(ncid, varid, particles%lon(1:nparticles))
    call nc_check_error (status, "Cannot read variable LON", abort=.false.)
    if (status/=0) return

    status = nf90_inq_varid(ncid, "RADIUS", varid)
    call nc_check_error (status, "Cannot find variable RADIUS", abort=.false.)
    if (status/=0) return
    status = nf90_get_var(ncid, varid, particles%radius(1:nparticles))
    call nc_check_error (status, "Cannot read variable RADIUS", abort=.false.)
    if (status/=0) return

    status = nf90_inq_varid(ncid, "DENSITY", varid)
    call nc_check_error (status, "Cannot find variable DENSITY", abort=.false.)
    if (status/=0) return
    status = nf90_get_var(ncid, varid, particles%density(1:nparticles))
    call nc_check_error (status, "Cannot read variable DENSITY", abort=.false.)
    if (status/=0) return

    status = nf90_inq_varid(ncid, "SEDIMENTATION", varid)
    call nc_check_error (status, "Cannot find variable SEDIMENTATION", abort=.false.)
    if (status/=0) return
    status = nf90_get_var(ncid, varid, particles%sedimentation(1:nparticles))
    call nc_check_error (status, "Cannot read variable SEDIMENTATION", abort=.false.)
    if (status/=0) return

    status = nf90_inq_varid(ncid, "SETTLING_VELOCITY", varid)
    call nc_check_error (status, "Cannot find variable SETTLING_VELOCITY", abort=.false.)
    if (status/=0) return
    status = nf90_get_var(ncid, varid, particles%tsv(1:nparticles))
    call nc_check_error (status, "Cannot read variable SETTLING_VELOCITY", abort=.false.)
    if (status/=0) return

    status = nf90_inq_varid(ncid, "HNO3", varid)
    call nc_check_error (status, "Cannot find variable HNO3", abort=.false.)
    if (status/=0) return
    status = nf90_get_var(ncid, varid, particles%hno3(1:nparticles))
    call nc_check_error (status, "Cannot read variable HNO3", abort=.false.)
    if (status/=0) return

    status = nf90_inq_varid(ncid, "H2O", varid)
    call nc_check_error (status, "Cannot find variable H2O", abort=.false.)
    if (status/=0) return
    status = nf90_get_var(ncid, varid, particles%h2o(1:nparticles))
    call nc_check_error (status, "Cannot read variable H2O", abort=.false.)
    if (status/=0) return

    status = nf90_inq_varid(ncid, "C_BOX_DENSITY_CHANGE", varid)
    call nc_check_error (status, "Cannot find variable C_BOX_DENSITY_CHANGE", abort=.false.)
    if (status/=0) return
    status = nf90_get_var(ncid, varid, particles%airparcel_density_change(1:nparticles))
    call nc_check_error (status, "Cannot read variable C_BOX_DENSITY_CHANGE", abort=.false.)
    if (status/=0) return

    status = nf90_inq_varid(ncid, "CLASS", varid)
    call nc_check_error (status, "Cannot find variable CLASS", abort=.false.)
    if (status/=0) return
    status = nf90_get_var(ncid, varid, particles%class(1:nparticles))
    call nc_check_error (status, "Cannot read variable CLASS", abort=.false.)
    if (status/=0) return

    status = nf90_inq_varid(ncid, "ICEbin", varid)
    call nc_check_error (status, "Cannot find variable ICEbin", abort=.false.)
    if (status/=0) return
    status = nf90_get_var(ncid, varid, particles%icebin(1:nparticles))
    call nc_check_error (status, "Cannot read variable ICEbin", abort=.false.)
    if (status/=0) return

    status = nf90_inq_varid(ncid, "NATbin", varid)
    call nc_check_error (status, "Cannot find variable NATbin", abort=.false.)
    if (status/=0) return
    status = nf90_get_var(ncid, varid, particles%natbin(1:nparticles))
    call nc_check_error (status, "Cannot read variable NATbin", abort=.false.)
    if (status/=0) return

    status = nf90_inq_varid(ncid, "SICE", varid)
    call nc_check_error (status, "Cannot find variable SICE", abort=.false.)
    if (status/=0) return
    status = nf90_get_var(ncid, varid, particles%sice(1:nparticles))
    call nc_check_error (status, "Cannot read variable SICE", abort=.false.)
    if (status/=0) return

    status = nf90_inq_varid(ncid, "SNAT", varid)
    call nc_check_error (status, "Cannot find variable SNAT", abort=.false.)
    if (status/=0) return
    status = nf90_get_var(ncid, varid, particles%snat(1:nparticles))
    call nc_check_error (status, "Cannot read variable SNAT", abort=.false.)
    if (status/=0) return

    status = nf90_inq_varid(ncid, "PART_ID", varid)
    call nc_check_error (status, "Cannot find variable PART_ID", abort=.false.)
    if (status/=0) return
    status = nf90_get_var(ncid, varid, particles%particle_id(1:nparticles))
    call nc_check_error (status, "Cannot read variable PART_ID", abort=.false.)
    if (status/=0) return
    
    ! close file
    status = nf90_close (ncid)

    ! set max particle id
    part_id_max = maxval(particles%particle_id)

    if (rank == 0) write(*,*) 'read ', nparticles, ' particles from initial file'

    particles%snatmax(1:nparticles) = 0.
    particles%sicemax(1:nparticles) = 0. 
    particles%tmin(1:nparticles) = 0.
    
  end subroutine read_clams_particles
  
  !*********************************************************************                  
  ! Read grid information from init file
  !*********************************************************************                  
  subroutine get_init_grid (status)

    USE netcdf

    USE messy_clams_global,        ONLY: prec, rank, initfile, init_vertcoorname, &
                                         buffersize
    USE messy_clamssedi_global,    ONLY: nlevs, lev_grid, lev_window, &
                                         lev_start, lev_end, &
                                         lev_min, lev_max, &
                                         r_coarse, r_high, lat_down, lat_up, &
                                         pos_nlevs, pos_lev_grid, &
                                         pos_lev_delta, pos_r_grid, &
                                         nhemi, lat_min, lat_max
    USE messy_clams_tools_ncutils, ONLY: nc_check_error
    USE messy_clams_tools_utils,   ONLY: uppercase, lowercase

    implicit none

    integer       :: status

    real(prec), dimension(:),    pointer :: lev_delta    => NULL()
    real(prec), dimension(:),    pointer :: r_grid       => NULL()
    real(prec), dimension(:),    pointer :: lev_grid_org => NULL()
    logical,    dimension(:),    pointer :: mask         => NULL()
    real(prec)    :: lev_act
    integer       :: ncid, dimid, varid, i
    integer       :: vert_type
    character(30) :: helpstr

    ! open init file
    status = nf90_open (initfile,nf90_nowrite,ncid, buffersize)
    CALL nc_check_error (status, 'Cannot open file '//TRIM(initfile),abort=.false.)
    if (status/=0) return
   
    ! get number of theta/zeta levels (NTHETAS/NZETAS)  => nlevs
    helpstr = 'N'//TRIM(uppercase(init_vertcoorname))//'S'  
    status = nf90_inq_dimid(ncid,TRIM(helpstr),dimid)     
    CALL nc_check_error (status, 'Cannot find dimension '//TRIM(helpstr)//' !',abort=.false.)
    if (status/=0) return
    status = nf90_inquire_dimension(ncid,dimid,len=nlevs)
    CALL nc_check_error (status, 'Cannot read dimension '//TRIM(helpstr)//' !',abort=.false.)
    if (status/=0) return
    
    ! read theta_min and theta_max 
    helpstr = 'exp_POS_'//TRIM(lowercase(init_vertcoorname))//'_min'
    status = nf90_get_att (ncid,nf90_global,trim(helpstr),lev_min)
    call nc_check_error (status,'Cannot find attribute '//trim(helpstr),abort=.false.)
    if (status/=0) return
    helpstr = 'exp_POS_'//TRIM(lowercase(init_vertcoorname))//'_max'  
    status = nf90_get_att (ncid,nf90_global,trim(helpstr),lev_max)
    call nc_check_error (status,'Cannot find attribute '//trim(helpstr),abort=.false.)
    if (status/=0) return

    ! read r_coarse, r_high 
    status = nf90_get_att (ncid,nf90_global,'exp_POS_r_coarse',r_coarse)
    call nc_check_error (status,'Cannot find attribute exp_POS_r_coarse',abort=.false.)
    if (status/=0) return
    status = nf90_get_att (ncid,nf90_global,'exp_POS_r_high',r_high)
    call nc_check_error (status,'Cannot find attribute exp_POS_r_high',abort=.false.)
    if (status/=0) return

    ! read lat_down, lat_up
    status = nf90_get_att (ncid,nf90_global,'exp_POS_lat_down',lat_down)
    call nc_check_error (status,'Cannot find attribute exp_POS_lat_down',abort=.false.)
    if (status/=0) return
    status = nf90_get_att (ncid,nf90_global,'exp_POS_lat_up',lat_up)
    call nc_check_error (status,'Cannot find attribute exp_POS_lat_up',abort=.false.)
    if (status/=0) return
  
    ! get THETA/ZETA_GRID => lev_grid
    helpstr = TRIM(uppercase(init_vertcoorname))//'_GRID'  
    status = nf90_inq_varid (ncid,TRIM(helpstr),varid)     
    CALL nc_check_error (status, 'Cannot find variable '//TRIM(helpstr)//' !',abort=.false.)
    if (status/=0) return
    allocate (lev_grid_org(nlevs))
    status = nf90_get_var (ncid,varid,lev_grid_org)
    call nc_check_error (status,'Cannot read variable '//trim(helpstr),abort=.false.)
    if (status/=0) return
    allocate (lev_grid(nlevs+2))
    lev_grid(1) = 0.
    lev_grid(2:nlevs+1) = lev_grid_org 
    lev_grid(nlevs+2) = 10000.


    ! get THETA/ZETA_DELTA => lev_delta
    helpstr = TRIM(uppercase(init_vertcoorname))//'_DELTA' ! THETA_DELTA
    status = nf90_inq_varid (ncid,TRIM(helpstr),varid)
    call nc_check_error (status,'Cannot find variable '//trim(helpstr),abort=.false.)
    if (status/=0) return
    allocate (lev_delta(nlevs))
    status = nf90_get_var (ncid,varid,lev_delta)
    call nc_check_error (status,'Cannot read variable '//trim(helpstr),abort=.false.)
    if (status/=0) return
  
    ! set lev_window
    allocate (lev_window(2,nlevs))
    lev_act = lev_min
    do i = 1, nlevs  
       lev_window(1,i) = lev_act
       lev_window(2,i) = lev_act+lev_delta(i)
       lev_act = lev_act+lev_delta(i)
    end do
    
    ! read r_grid 
    status = nf90_get_att (ncid,nf90_global,'exp_POS_vert_type',vert_type)
    if (status==nf90_noerr) then
       if (vert_type==3) then
          status = nf90_inq_varid(ncid,'R_GRID',varid)
          if (status==nf90_noerr) then
             allocate (r_grid(nlevs))
             status = nf90_get_var (ncid,varid,r_grid)
             call nc_check_error (status,'Cannot read variable R_GRID',abort=.false.)
             if (status/=0) return
          endif
       endif
    endif


    ! set pos_nlevs, pos_lev_grid, pos_lev_delta, pos_r_grid (used in sub. create_positions)
    ! -> grid between lev_start and lev_end
    allocate (mask(nlevs+2))
    mask = .false.
    where (lev_start<=lev_grid_org .and. lev_grid_org<=lev_end)
       mask = .true.
    end where
    pos_nlevs = count(mask)
    allocate (pos_lev_grid (pos_nlevs))
    allocate (pos_lev_delta(pos_nlevs))
    if (associated(r_grid)) allocate (pos_r_grid   (pos_nlevs))
    pos_lev_grid  = pack(lev_grid_org,mask)
    pos_lev_delta = pack(lev_delta,mask)
    if (associated(r_grid)) pos_r_grid  = pack(r_grid,mask)
    deallocate (mask)

    deallocate (lev_grid_org)

    if (rank==0)  then
       write (*,*)
       write (*,*) 'SEDI: subroutine get_init_grid:'
       write (*,*)
       write (*,*) 'nlevs:',nlevs
       write (*,*)
       write (*,*) 'lev_grid:',lev_grid
       write (*,*)
       write (*,*) 'lev_delta:',lev_delta
       write (*,*)
       write (*,*) 'lev_window:',lev_window
       write (*,*)
       if (associated(r_grid)) write (*,*) 'r_grid:',r_grid
       write (*,*)
       write (*,*) 'pos_nlevs:',pos_nlevs
       write (*,*)
       write (*,*) 'pos_lev_grid:',pos_lev_grid
       write (*,*)
       write (*,*) 'pos_lev_delta:',pos_lev_delta
       write (*,*)
       if (associated(pos_r_grid)) write (*,*) 'pos_r_grid:',pos_r_grid
       write (*,*)
    endif

    ! close init file
    status = nf90_close(ncid)

    ! deallocate local arrays
    deallocate (lev_delta)
    if (associated(r_grid)) deallocate (r_grid)

    ! set nhemi 
    if (abs(lat_max - 90) .le. 0.1) then
       nhemi = .true.
       if (rank == 0)  write(*,*) 'Fine grid found in northern hemisphere'
    else if (abs(lat_min + 90) .le. 0.1) then
       if (rank == 0) write(*,*) 'Fine grid found in southern hemisphere'
       nhemi = .false.
    else
       if (rank == 0)  write(*,*) 'WARNING: Cannot determine hemisphere, assuming northern hemisphere !!!'
       nhemi = .true.
    endif
    

  end subroutine get_init_grid

  !*************************************************************************
  ! Determine number of airparcels per lev in initfile
  !*************************************************************************
  
  subroutine get_nairparcels_init (status)

    use netcdf
    use messy_clams_tools_ncutils, only: nc_check_error
    USE messy_clams_tools_utils,   ONLY: uppercase
    use messy_clamssedi_global,    only: nlevs, nhemi, lev_window,  &
                                         airparcels, &
                                         nairparcels50_lev_init   !, nairparcels_lev

    use messy_clams_global,        only: prec, initfile, init_vertcoorname, buffersize
    
    implicit none
    
    integer, intent(out) :: status
    integer              :: ncid, dimid, varid

    integer                           :: ilev, nparts
    logical,    dimension(:), pointer :: mask

    real(prec), dimension(:), pointer :: lat_init, lev_init

    !allocate (nairparcels_lev(nlevs))
    allocate (nairparcels50_lev_init(nlevs))
    !nairparcels_lev   = 0
    nairparcels50_lev_init = 0

    !------------------------------------------------
    ! read initfile
    !------------------------------------------------ 

    ! open file
    status = nf90_open (initfile, nf90_nowrite, ncid, buffersize)
    call nc_check_error (status, "Error on open file "//initfile, abort=.false.)
    if (status/=0) return

    ! determine dimension (nparts)
    status = nf90_inq_dimid (ncid, "NPARTS", dimid)
    call nc_check_error (status, "Cannot find dimension NPARTS", abort=.false.)
    if (status/=0) return
    status = nf90_inquire_dimension (ncid,dimid,len=nparts)
    call nc_check_error (status, "Cannot read dimension NPARTS", abort=.false.)
    if (status/=0) return
    
    ! read LAT and LEV (THETA/ZETA) 
    allocate (lev_init(nparts))
    allocate (lat_init(nparts))
    status = nf90_inq_varid (ncid, "LAT", varid)
    call nc_check_error (status, "Cannot find variable LAT", abort=.false.)
    if (status/=0) return
    status = nf90_get_var (ncid, varid, lat_init)
    call nc_check_error (status, "Cannot read variable LAT", abort=.false.)
    if (status/=0) return
    status = nf90_inq_varid (ncid, uppercase(init_vertcoorname), varid)
    call nc_check_error (status, "Cannot find variable "//uppercase(init_vertcoorname), abort=.false.)
    if (status/=0) return
    status = nf90_get_var (ncid, varid, lev_init)
    call nc_check_error (status, "Cannot read variable "//uppercase(init_vertcoorname), abort=.false.)
    if (status/=0) return
    
    ! close file
    status = nf90_close (ncid)



    !------------------------------------------------
    ! Determine number of airparcels per lev
    !------------------------------------------------
    ! allocate(mask(nparts))
    ! do ilev = 1, nlevs
    !    mask=.false.
    !    where (lev_window(1,ilev)<=airparcels%lev .and. &
    !           airparcels%lev < lev_window(2,ilev)) mask=.true.
    !    nairparcelss_lev(ilev) = count(mask)
    ! end do
    ! deallocate (mask)


    !---------------------------------------------------------------------------
    ! Determine number of airparcels north of 50 deg respectively south of -50 deg
    !---------------------------------------------------------------------------
    allocate(mask(nparts))
    do ilev = 1, nlevs
       mask=.false.
       if (nhemi) then
          where (lev_window(1,ilev)<=lev_init .and. &
                 lev_init < lev_window(2,ilev) .and. &
                 lat_init > 50.) mask=.true.
       else
          where (lev_window(1,ilev)<=lev_init .and. &
                 lev_init < lev_window(2,ilev) .and. &
                 lat_init < -50.) mask=.true.
       endif
       nairparcels50_lev_init(ilev) = count(mask)
    end do
    deallocate (mask)


    !------------------------------------------------
    ! Clean up
    !------------------------------------------------
    deallocate (lat_init)
    deallocate (lev_init)


  end subroutine get_nairparcels_init

  !*************************************************************************
  ! Einlesen der Chemie-Datei und Berechnung der Nukleationsrate
  !
  ! NAT and ice number densities are calculated following Hoyle et al. (2013)
  ! and Engel et al. (2013), respectively. The implementation into CLaMS was
  ! done with the help of a lookup table, which contains 1000 contact angle
  ! bins (in steps of 0.1 deg) and 31 temperature bins (in steps of 1 K)
  !
  ! Jens-Uwe Grooss,  March 2012 -> NAT
  ! Ines Tritscher, February 2015 -> ICE
  !
  !*************************************************************************
  subroutine read_nucleation_table

    use messy_clamssedi_global,   only: nbins, ntbins, snat_table, & 
                                        xnnat_table, xnice_table, sice_table, &
                                        ice_nuc_table, nat_nuc_table
    
    implicit none

    integer :: itime, i, j, ibin, tbin
    real    :: c

    allocate(snat_table(nbins,ntbins), xnnat_table(nbins))

    open(77,file=nat_nuc_table,form='formatted',status='old')
    do i=1,nbins
        read(77,*) c, xnnat_table(i),(snat_table(i,j),j=1,ntbins)
    enddo
    close(77)

    allocate(sice_table(nbins,ntbins), xnice_table(nbins))

    open(78,file=ice_nuc_table,form='formatted',status='old')
    do i=1,nbins
        read(78,*) c, xnice_table(i),(sice_table(i,j),j=1,ntbins)
    enddo
    close(78)

  end subroutine read_nucleation_table



End Module messy_clamssedi_data_io

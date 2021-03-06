MODULE messy_clamsbmix
!***********************************************************************!
!
! BMIX -- Replace boundaries                                    
!
! This module contains the following subroutines:
!    subroutine bmix (status,lat,lon,lev,specarr)
!    subroutine error_handler (error)
!    subroutine clamsbmix_read_nml (status, iou)
!    subroutine define_boundaries
!    subroutine get_orography
!    subroutine delete_boundaries (lat,lon,lev,specarr)
!    subroutine get_lower_bound (status, lat, lon, lev, specarr)
!    subroutine get_upper_bound (status,lat,lon,lev,specarr)
!    subroutine get_vertical_bound (status,lat,lon,lev,specarr,north)
!
!***********************************************************************!
   
  IMPLICIT NONE

  CHARACTER(LEN=*),PARAMETER, PUBLIC :: modstr = 'clamsbmix'
  CHARACTER(LEN=*),PARAMETER, PUBLIC :: modver = '1.0'



  PUBLIC :: bmix
  PUBLIC :: clamsbmix_read_nml
  

!----------- 
CONTAINS
!----------- 

  SUBROUTINE bmix (status, lat, lon, lev, state, specarr, &
                   lat_low, lon_low, lev_low, specarr_low, &
                   lat_upper, lon_upper, lev_upper, specarr_upper, lcoupled)
    
    USE messy_clams_global,         ONLY: YEAR_NEXT, MONTH_NEXT, DAY_NEXT &
                                        , HOUR_NEXT                       &
                                        , rank, PREC, DP, filenamelen,    &
                                          dnparts, species_type, nspec,   &
                                          dates30, &
                                          met_dir, met_prefix, &
                                          theta_dir, theta_prefix
    USE messy_clamsmix_global,      ONLY: adapt_par
    USE messy_clamsbmix_global,     ONLY: file_bounds, &
                                          lev_in_down, lev_in_up, &
                                          lat_in_down, lat_in_up, &
                                          delta_lev, nlevs, lev0, &
                                          replace_low, replace_up, &
                                          replace_north, replace_south, &
                                          time_init_value
    USE messy_clams_tools_utils,    ONLY: get_metfilename    
    USE messy_clams_tools_ncutils,  ONLY: nc_get_level    
    USE messy_clams_tools_dateconv, ONLY: ymds2js

    IMPLICIT NONE

!!!!! lat,lon,lev,specarr(i)%values sind mit dnparts_max_shuffle dimensioniert !!!
    INTEGER                                   :: status
    REAL(PREC),         DIMENSION(:), POINTER :: lat, lon, lev
    real(prec),         dimension(:), pointer :: state
    TYPE(species_type), DIMENSION(:), POINTER :: specarr
    REAL(PREC),         DIMENSION(:), POINTER :: lat_low, lon_low, lev_low
    REAL(PREC),         DIMENSION(:), POINTER :: lat_upper, lon_upper, lev_upper
    TYPE(species_type), DIMENSION(:), POINTER :: specarr_low
    TYPE(species_type), DIMENSION(:), POINTER :: specarr_upper
    LOGICAL                                   :: lcoupled

    INTEGER,  parameter  :: maxlev = 720

    CHARACTER(8)   :: datestr
    CHARACTER(filenamelen) :: file_data_met, file_theta

    INTEGER        :: i
    INTEGER        :: finish, start, counts_per_sec
    REAL(DP)       :: seconds  
    REAL(PREC)     :: dlev

    status = 0  ! no error

    CALL system_clock (count=start)

    if (rank==0) THEN
       WRITE (*,*)
       WRITE (*,*) 'START OF BMIX'
       WRITE (*,*)
    ENDIF
    
    if (dates30) then
       write (*,*)
       write (*,*) 'ACHTUNG: 30-Tage-Monate !!!'
       write (*,*)
    endif
    
    if (rank==0) then
       write (*,*) 
       write (*,*) 'Boundaries from:    ', TRIM(file_bounds)
       write (*,*) 
       write (*,*) 'boundaries:'
       write (*,*) lev_in_down, lev_in_up, lat_in_up, lat_in_down
       write (*,*) 
       write (*,*) 'delta_lev=',delta_lev
       write (*,*) 
       write (*,*) 'lev_grid:'
       write (*,*) adapt_par%lev_grid
       write (*,*) 'lev_delta:'
       write (*,*) adapt_par%lev_delta
       write (*,*) 'lat_down, lat_up, lat_min, lat_max:'
       write (*,*) adapt_par%lat_down, adapt_par%lat_up, adapt_par%lat_min, adapt_par%lat_max
       write (*,*) 'nlevs=',nlevs
       write (*,*) 
    endif
  
    ! Check if the number of layers is large enough
    if (nlevs > 0 .and. nlevs <= 2) then 
       print *, 'Use more than 2 layers for a 3d run'
       status = 301
       return
    endif

    file_data_met = ""
    file_theta = ""
    if (.not. lcoupled) then    
       ! name of wind file
       file_data_met = get_metfilename (met_prefix, met_dir, &
                           YEAR_NEXT, MONTH_NEXT, DAY_NEXT, HOUR_NEXT)
       if (rank==0) write (*,*) 'file_data_met: ',trim(file_data_met)
       ! name of isentropic file
       file_theta = get_metfilename (theta_prefix, theta_dir, &
                           YEAR_NEXT, MONTH_NEXT, DAY_NEXT, HOUR_NEXT)
       if (rank==0) write (*,*) 'file_theta: ',trim(file_theta)
    endif


!!!!! ???
    ! Get time_init
     time_init_value = ymds2js (YEAR_NEXT,MONTH_NEXT,DAY_NEXT,HOUR_NEXT*3600)
     if (rank==0) write(*,*) 'time_init_value=',time_init_value

    !----------------------------------------------------------------
    ! get orography (get latlon_grid, determine lev0)
    !----------------------------------------------------------------
    if (delta_lev>0) then

       if (rank==0) write (*,*) 'get orography'
       call get_orography (status, file_data_met)
       if (status/=0) then
          write (*,*) 'Error in get_orography!!!'
          return
       endif
    endif

    !----------------------------------------------------------------
    ! delete boundaries 
    !----------------------------------------------------------------
    !if (rank==0) write (*,*) 'call delete_boundaries'
    call delete_boundaries (lat,lon,lev,state,specarr)

    !----------------------------------------------------------------
    ! add points on lower boundary
    !----------------------------------------------------------------
    if (replace_low) then
       !if (rank==0) write (*,*) 'get lower boundary layer'
       call get_lower_bound (status, file_theta, lat, lon, lev, &
                             state, specarr, &
                             lat_low, lon_low, lev_low, specarr_low, lcoupled)
       if (status/=0) then
          write (*,*) 'Error in get_lower_bound!!!'
          return
       endif
    endif

    !----------------------------------------------------------------
    ! add points on upper boundary
    !----------------------------------------------------------------
    if (replace_up) then
       !if (rank==0) write (*,*) 'get upper boundary layer'
       call get_upper_bound (status, file_theta, lat, lon, lev, &
                             state, specarr, &
                             lat_upper, lon_upper, lev_upper, specarr_upper, lcoupled)
       if (status/=0) then
          write (*,*) 'Error in get_upper_bound!!!'
          return
       endif
    endif

    !----------------------------------------------------------------
    ! add points on north and/or south boundaries
    !----------------------------------------------------------------
    if (replace_south) then
       !if (rank==0) write (*,*) 'get south boundary '
       call get_vertical_bound (status, lat, lon, lev, specarr,.false.)
       if (status/=0) then
          write (*,*) 'Error in get_vertical_bound!!!'
          return
       endif
    endif
    if (replace_north) then
       !if (rank==0) write (*,*) 'get north boundary '
       call get_vertical_bound (status, lat, lon, lev, specarr,.true.)
       if (status/=0) then
          write (*,*) 'Error in get_vertical_bound!!!'
          return
       endif
    endif

    ! clean up
    !if (rank==0) write (*,*) 'clean up'
    if (delta_lev>0) then
       deallocate (lev0)
    endif

    call system_clock (COUNT_RATE=counts_per_sec)
    call system_clock (count=finish)
    if (finish > start) then
       seconds = float(finish-start) / float(counts_per_sec)
    else
       seconds = 0.0
    endif
    
    if (rank==0) then
       write (*,*)
       write (*,*) 'Normal termination of bmix'
       !write (*,*)
       !write (*,*) 'System clock runs at ', counts_per_sec, 'ticks per second'
       write (*,*)
       write (*,'(A,F10.2,A)') 'This job has taken ',seconds,' seconds to execute.' 
    endif

  END SUBROUTINE bmix

  !*******************************************************************
  ! 
  !*******************************************************************
  subroutine error_handler (error)

    USE messy_clams_global, ONLY: rank

    implicit none

    integer :: error

    if (rank==0) then
       write (*,*) '***********************************************************'

       if (error>=20 .and. error<30) then
          write (*,*) 'ERROR in boundfile list :'
       else
          write (*,*) 'ERROR:'
       endif
       
       select case (error)
       case (8)
          write (*,*) '    Specify 0,1,2 or 3 for interpol_from_init !!!'
       case (20)
          write (*,*) '    Cannot open bmix-boundfile list !!!'
       case (21)
          write (*,*) '    Bmix bound files: Cannot read species !!!'
       case (22)
          write (*,*) '    Bmix bound files: Enter 1 or 2 for lower or upper boundary !!!'
       case (23)
          write (*,*) '    Bmix bound files: Enter 1 or 2 for replace/add !!!'
       case (24)
          write (*,*) '    Bmix bound files: Cannot read name of bound file !!!'
!!$       case (40)
!!$          write (*,*) '    Bound files: Error while reading file_clamsclim_spec !!!'
       case (40)
          write (*,*) '    Cannot open clams-boundfile list !!!'
       case (41)
          write (*,*) '    Clams bound files: Cannot read speciesname !!!'
       case (42)
          write (*,*) '    Clams bound files: Cannot read speciesname (name in boundfile) !!!'
       case (43)
          write (*,*) '    Clams bound files: Cannot read lower boundary !!!'
       case (44)
          write (*,*) '    Clams bound files: Cannot read coordinate for lower boundary !!!'
       case (45)
          write (*,*) '    Clams bound files: Enter 0-9 for action on lower boundary !!!'
       case (46)
          write (*,*) '    Clams bound files: Cannot read upper boundary !!!'
       case (47)
          write (*,*) '    Clams bound files: Cannot read coordinate for upper boundary !!!'
       case (48)
          write (*,*) '    Clams bound files: Enter 1-5 for action on upper boundary !!!'
       case (49)
          write (*,*) '    Clams bound files: Cannot read name of bound file !!!'
       case (50)
          write (*,*) '    Clams bound files: Cannot read startyear !!!'
       case (51)
          write (*,*) '    Clams bound files: Cannot read endyear !!!'
      case default
          write (*,*) 'error no ',error
       end select
      write (*,*) '***********************************************************'
    endif

  end subroutine error_handler

  !****************************************************************************
  ! Read namelists
  !****************************************************************************
  SUBROUTINE clamsbmix_read_nml(status, iou)

    USE messy_main_tools,       ONLY: read_nml_open, read_nml_check, read_nml_close
    USE messy_clams_global,     ONLY: rank, nchemspec, ldiagout
    USE messy_clamsbmix_global

    IMPLICIT NONE

    !I/O
    INTEGER, INTENT(OUT) ::   status
    INTEGER, INTENT(IN)  ::   iou
    !LOCAL
    CHARACTER(LEN=*),PARAMETER :: substr='clamsbmix_read_nml'
    CHARACTER(150)       :: line
    LOGICAL              :: lex      ! file exists ?
    INTEGER              :: fstat    ! file status
    INTEGER              :: ios, pos, i
    LOGICAL              :: l_print  ! write control output

    NAMELIST /CTRL/ timestep_bmix, &
                    lev_in_down, lev_in_up,    &
                    lat_in_down,   lat_in_up,          &
                    delta_lev,   interpol_from_init, &
                    file_bounds,   dir_boundfiles,     &
                    bmix_boundlist, clams_boundlist,    &
                    switch_EMAC_H2O, EMAC_H2O_z, &
                    max_dist

    status = 1 !ERROR

    if (rank==0 .and. ldiagout) then
       l_print = .true.
    else
       l_print = .false.
    endif

    CALL read_nml_open(lex, substr, iou, 'CTRL', modstr, l_print)
    IF (.not.lex) return    ! <modstr>.nml does not exist

    READ(iou, NML=CTRL, IOSTAT=fstat)
    CALL read_nml_check(fstat, substr, iou, 'CTRL', modstr, l_print)
    IF (fstat /= 0) return  ! error while reading namelist
    
    CALL read_nml_close(substr, iou, modstr, l_print)

    status = 0 !NO ERROR

    ! interpol_from_init = 0|1|2|3
    if (interpol_from_init <= 0) then
       if (rank==0) write (*,*) 'No interpolation from init-File !'
    elseif (interpol_from_init==1) then
       if (rank==0) write (*,*) 'Interpolate lower boundary from init-File !'
    elseif (interpol_from_init==2) then
       if (rank==0) write (*,*) 'Interpolate upper boundary from init-File !'
    elseif (interpol_from_init==3) then
       if (rank==0) write (*,*) 'Interpolate lower and upper boundaries from init-File !'
    else
       status = 8
       call error_handler(status)
    endif
   

    !-------------------------------------------------------------------------
    ! read list of boundfiles (used in BMIX)
    !-------------------------------------------------------------------------
    nbmixbounds = 0
    if (bmix_boundlist/='') then
       open (iou,file=bmix_boundlist,status="OLD",iostat=ios)
       if (ios /= 0) then 
          status = 20
          call error_handler(status)
       else
          read (iou,'(A)',iostat=ios) line
          line = adjustl(line)
          pos = index(line,'!')
          if (pos /= 0)  line = line(1:pos-1)
          do while (ios==0 .and. nbmixbounds<maxboundspec .and. line/='') 
             ! name of species
             read (line,*,iostat=ios) bmixbound(nbmixbounds+1)%spec
             if (ios/=0) then
                status = 21
                call error_handler(status)
                exit
             endif
             pos = index(line,' ')
             line = adjustl(line(pos+1:))
             ! lower or upper boundary (1/2)
             read (line,*,iostat=ios) bmixbound(nbmixbounds+1)%no
             if (ios/=0) then
                status = 22
                call error_handler(status)
                exit
             endif
             if (bmixbound(nbmixbounds+1)%no/=1 .and. bmixbound(nbmixbounds+1)%no/=2) then
                status = 22
                call error_handler(status)
                exit
             endif
             pos = index(line,' ')
             line = adjustl(line(pos+1:))
             ! replace or add (1/2)
             read (line,*,iostat=ios) bmixbound(nbmixbounds+1)%add
             if (ios/=0) then
                status = 23
                call error_handler(status)
                exit
             endif
             if (bmixbound(nbmixbounds+1)%add/=1 .and. bmixbound(nbmixbounds+1)%add/=2) then
                status = 23
                call error_handler(status)
                exit
             endif
             pos = index(line,' ')
             line = adjustl(line(pos+1:))
             ! name of bound file
             read (line,*,iostat=ios) bmixbound(nbmixbounds+1)%file
             if (ios/=0) then
                status = 24
                call error_handler(status)
                exit
             endif
             pos = index(line,' ')
             line = adjustl(line(pos+1:))
             ! startyear and endyear (optional)
             read (line,*,iostat=ios) bmixbound(nbmixbounds+1)%startyear, bmixbound(nbmixbounds+1)%endyear
             ! number of bound species
             nbmixbounds = nbmixbounds + 1
             ! read next line
             read (iou,'(A)',iostat=ios) line
             line = adjustl(line)
             pos = index(line,'!')
             if (pos /= 0)  line = line(1:pos-1)
          end do
          
       endif
       close (iou)
    endif

    if (nbmixbounds==0) then
       if (rank==0) write (*,*) 'No bound files for BMIX specified !'
    else
       if (rank==0) then
          write (*,*) 'nbmixbounds=',nbmixbounds
          do i = 1, nbmixbounds
             write (*,'(2A,I4,A,I4,3A,I4,A,I4)') bmixbound(i)%spec,'  ',bmixbound(i)%no,'  ', &
                         bmixbound(i)%add,'  ',trim(bmixbound(i)%file),'   ', &
                         bmixbound(i)%startyear,'   ',bmixbound(i)%endyear
          enddo
       endif
       if (maxval(bmixbound(1:nbmixbounds)%add)==2) then
          if (interpol_from_init<=0) then
             if (rank==0) then
                write (*,*)
                write (*,*) 'WARNING: "add to species" (2) is choosen, ', &
                  ' but no interpolation from init file !!!'
                write (*,*) '         Values are added to start values from bound file !!!'
                write (*,*)
                !status = 30
             endif
          endif
       endif
    endif
    write (*,*)

    !-------------------------------------------------------------------------
    ! read list of boundfiles (replacement after BMIX)
    !-------------------------------------------------------------------------
    nclamsbounds = 0
    if (clams_boundlist/='') then
       open (iou,file=clams_boundlist,status="OLD",iostat=ios)
       if (ios /= 0) then 
          status = 40
          call error_handler(status)
       else
          read (iou,'(A)',iostat=ios) line
          line = adjustl(line)
          pos = index(line,'!')
          if (pos /= 0)  line = line(1:pos-1)
          do while (ios==0 .and. nclamsbounds<maxboundspec .and. line/='') 

             ! name of species 
             read (line,*,iostat=ios) clamsbound(nclamsbounds+1)%spec
             if (ios/=0) then
                status = 41
                call error_handler(status)
                exit
             endif
             pos = index(line,' ')
             line = adjustl(line(pos+1:))

             ! name of species in boundfile
             read (line,*,iostat=ios) clamsbound(nclamsbounds+1)%spec_bf
             if (ios/=0) then
                status = 42
                call error_handler(status)
                exit
             endif
             pos = index(line,' ')
             line = adjustl(line(pos+1:))

             ! lower boundary 
             read (line,*,iostat=ios) clamsbound(nclamsbounds+1)%lbound
             if (ios/=0) then
                status = 43
                call error_handler(status)
                exit
             endif
             pos = index(line,' ')
             line = adjustl(line(pos+1:))

             ! action for lower boundary
             read (line,*,iostat=ios) clamsbound(nclamsbounds+1)%lno
             if (ios/=0) then
                status = 45
                call error_handler(status)
                exit
             endif
             if (clamsbound(nclamsbounds+1)%lno<0 .or. clamsbound(nclamsbounds+1)%lno>9) then
                status = 45
                call error_handler(status)
                exit
             endif
             pos = index(line,' ')
             line = adjustl(line(pos+1:))

             ! upper boundary 
             read (line,*,iostat=ios) clamsbound(nclamsbounds+1)%ubound
             if (ios/=0) then
                status = 46
                call error_handler(status)
                exit
             endif
             pos = index(line,' ')
             line = adjustl(line(pos+1:))

             ! action for upper boundary
             read (line,*,iostat=ios) clamsbound(nclamsbounds+1)%uno
             if (ios/=0) then
                status = 48
                call error_handler(status)
                exit
             endif
             pos = index(line,' ')
             line = adjustl(line(pos+1:))

             ! start year
             read (line,*,iostat=ios) clamsbound(nclamsbounds+1)%startyear
             if (ios/=0) then
                status = 50
                call error_handler(status)
                exit
             endif
             pos = index(line,' ')
             line = adjustl(line(pos+1:))
            
             ! end year
             read (line,*,iostat=ios) clamsbound(nclamsbounds+1)%endyear
             if (ios/=0) then
                status = 51
                call error_handler(status)
                exit
             endif
             pos = index(line,' ')
             line = adjustl(line(pos+1:))
            
             ! name of bound file
             read (line,'(A)',iostat=ios) clamsbound(nclamsbounds+1)%file
             if (ios/=0) then
                status = 49
                call error_handler(status)
                exit
             endif
             pos = index(clamsbound(nclamsbounds+1)%file,'!')
             if (pos /= 0)  clamsbound(nclamsbounds+1)%file = clamsbound(nclamsbounds)%file(1:pos-1)

             ! number of bound species
             nclamsbounds = nclamsbounds + 1

             ! read next line
             read (iou,'(A)',iostat=ios) line
             line = adjustl(line)
             pos = index(line,'!')
             if (pos /= 0)  line = line(1:pos-1)
          end do
          
       endif
       close (iou)
    endif

    if (rank==0) write (*,*) 'nclamsbounds=',nclamsbounds
    if (nclamsbounds==0) then
       if (rank==0) write (*,*) 'No other bound files (used after BMIX) specified !'
    else
       if (rank==0) then
          do i = 1, nclamsbounds
             write (*,'(2A,2(F7.2,A,I3,A),2(I4,A),A)') clamsbound(i)%spec,'  ',&
                  clamsbound(i)%lbound,'  ',clamsbound(i)%lno,'  ', &
                  clamsbound(i)%ubound,'  ',clamsbound(i)%uno,'  ', &
                  clamsbound(i)%startyear,'  ',clamsbound(i)%endyear,'  ',&
                  trim(clamsbound(i)%file)
          enddo
       endif
    endif


  END SUBROUTINE clamsbmix_read_nml

  !****************************************************************************
  !
  !****************************************************************************
  SUBROUTINE replace_yyyy_in_filename (filename, new_filename)

    USE messy_clams_global, ONLY: filenamelen, YEAR_NEXT   
    
    implicit none

    INTEGER :: pos
    CHARACTER(filenamelen) :: filename, new_filename
    CHARACTER(filenamelen) :: str1, str2
    
    pos = index(filename,'yyyy')
    if (pos==0) pos = index(filename,'YYYY')
    if (pos/=0) then
       str1 = filename(1:pos-1)
       str2 = filename(pos+4:)
       write (new_filename,'(A,I4.4,A)') trim(str1),YEAR_NEXT,trim(str2)
    else
       new_filename = trim(filename)
    endif

  END SUBROUTINE replace_yyyy_in_filename
  
  !****************************************************************************
  ! define lev- and lat-boundaries (called in clamsbmix_global_end)
  !****************************************************************************
  SUBROUTINE define_boundaries

    USE messy_clams_global,        ONLY: rank
    USE messy_clamsmix_global,     ONLY: adapt_par
    USE messy_clamsbmix_global,    ONLY: nlevs, lev_down, lev_up, lat_down, lat_up,  &
                                         lev_in_down, lev_in_up, lat_in_down, lat_in_up, &
                                         replace_low, replace_up, replace_north, replace_south

    implicit none

    nlevs = adapt_par%nlevs

    ! Define lev-boundaries
    lev_down= adapt_par%lev_min
    lev_up=adapt_par%lev_max
    if (nlevs == 0) then  ! isentropic case
       lev_in_down=lev_down-5.
       lev_in_up=lev_up+5.
    else
       if (lev_down < lev_in_down  .and. &
            lev_in_down < lev_down+adapt_par%lev_delta(1)) then
          lev_in_down = lev_down+adapt_par%lev_delta(1)
       elseif (lev_in_down < lev_down) then
          lev_in_down = lev_down
       endif
       if (lev_in_up < lev_up  .and. &
            lev_in_up > lev_up-adapt_par%lev_delta(nlevs)) then
          lev_in_up = lev_up-adapt_par%lev_delta(nlevs)
       elseif (lev_in_up > lev_up) then
          lev_in_up = lev_up
       endif
    endif
    
    ! Define lat-boundaries
    lat_down=adapt_par%lat_down
    if (lat_down <= -88.) lat_in_down=lat_down
    lat_up=adapt_par%lat_up
    if (lat_up >= 88.) lat_in_up=lat_up
    
    replace_low = (lev_in_down > lev_down)
    replace_up = (lev_in_up < lev_up)
    replace_north = (lat_in_up < lat_up)
    replace_south = (lat_in_down > lat_down)

    if (rank==0) then
       write (*,*)
       write (*,*) 'replace_low, replace_up, replace_north, replace_south: ', &
            replace_low, replace_up, replace_north, replace_south
       if (nlevs == 0) then 
          write (*,*) 'Main layer', lev_in_down, '---', lev_in_up
          write (*,*) 'Northern vert. bound.', lat_in_up, '---', lat_up
          write (*,*) 'Southern vert. bound.', lat_down, '---', lat_in_down 
          write (*,*)
       else 
          write (*,*) 'Lower layer  ', lev_down, '---', lev_in_down
          write (*,*) 'Main layer   ', lev_in_down, '---', lev_in_up
          write (*,*) 'Upper layer  ', lev_in_up, '---', lev_up
          write (*,*) 'Northern vert. bound.', lat_in_up, '---', lat_up
          write (*,*) 'Southern vert. bound.', lat_down, '---', lat_in_down 
          write (*,*)
       endif
    endif
        
  END SUBROUTINE define_boundaries

  !****************************************************************************
  ! Get orography from isentropic file
  !****************************************************************************
  SUBROUTINE get_orography (status, file)

    use messy_clams_global,        only: prec, buffersize, eps
    use messy_clamsmix_global,     only: ctrl_out, gridtype
    USE messy_clamsbmix_global,    only: lev0, lev_in_down, delta_lev, &
                                         latlon_grid
    use messy_clams_tools_utils,   only: uppercase, lowercase
    use messy_clams_tools_ncutils, only: nc_get_vertcoorname, nc_check_error, &
                                         nc_get_var_cf
    
    use netcdf
    
    implicit none 
    
    integer                    :: status           
    character*(*), intent(in)  :: file
    
    ! local variables
    real(prec), dimension(:), allocatable       :: lev
    real(prec), dimension(:,:,:,:), allocatable :: press
    real(prec)                                  :: mdi
    real(prec), dimension(1)                    :: lat0, lat1, lon0, lon1
    integer                                     :: ntime, nlat, nlon, nlev, &
                                                   ilon, ilat, ilev, &
                                                   ncid, rcode, dimid, varid
    character(80)                               :: vertcoorname

    status = 0 ! no error
 
    ! file name
    if (ctrl_out) print *
    if (ctrl_out) print *, 'NetCDF file in get_orography: ', file
    
    ! open file
    status= nf90_open(file, nf90_nowrite, ncid, buffersize)
    call nc_check_error (status,'Error on open file '//trim(file),abort=.false.)
    if (status/=0) return
    
    ! get name of vertical coordinate 
    call nc_get_vertcoorname (file,vertcoorname)
    
    ! read dimensions
    status = nf90_inq_dimid (ncid,'time',dimid)       
    call nc_check_error (status,'Cannot find dimension time',abort=.false.)
    if (status/=0) return
    status = nf90_inquire_dimension (ncid,dimid,len=ntime)
    call nc_check_error (status,'Cannot read dimension time',abort=.false.)
    if (status/=0) return
    status = nf90_inq_dimid (ncid,'lat',dimid)    
    call nc_check_error (status,'Cannot find dimension lat',abort=.false.)
    if (status/=0) return
    status = nf90_inquire_dimension (ncid,dimid,len=nlat)
    call nc_check_error (status,'Cannot read dimension lat',abort=.false.)
    if (status/=0) return
    status = nf90_inq_dimid (ncid,'lon',dimid)    
    call nc_check_error (status,'Cannot find dimension lon',abort=.false.)
    if (status/=0) return
    status = nf90_inquire_dimension (ncid,dimid,len=nlon)
    call nc_check_error (status,'Cannot read dimension lon',abort=.false.)
    if (status/=0) return
    status = nf90_inq_dimid (ncid,TRIM(lowercase(vertcoorname)),dimid)    
    call nc_check_error (status,'Cannot find dimension '// &
         TRIM(lowercase(vertcoorname)),abort=.false.)
    if (status/=0) return
    status = nf90_inquire_dimension (ncid,dimid,len=nlev)
    call nc_check_error (status,'Cannot read dimension '// &
         TRIM(lowercase(vertcoorname)),abort=.false.)
    if (status/=0) return
    
    ! allocate arrays
    allocate (lev(nlev))
    allocate (press(nlon,nlat,nlev,ntime))
    allocate (lev0(nlon,nlat))
    
    ! read vertical coordinate
    status = nf90_inq_varid (ncid,TRIM(lowercase(vertcoorname)),varid)
    call nc_check_error (status,'Cannot find variable '// &
         TRIM(lowercase(vertcoorname)),abort=.false.)
    if (status/=0) return
    status = nf90_get_var (ncid,varid,lev)
    call nc_check_error (status,'Cannot read variable '// &
         TRIM(lowercase(vertcoorname)),abort=.false.)
    if (status/=0) return

    ! read pressure
    status = nf90_inq_varid (ncid,'PRESS',varid)
    call nc_check_error (status,'Cannot find PRESS',abort=.false.)
    if (status/=0) return
    !status = nf90_get_var (ncid,varid,press)
    call nc_get_var_cf (status, ncid, varid, 'PRESS', press(:,:,:,1))
    call nc_check_error (status,'Cannot read PRESS',abort=.false.)
    if (status/=0) return

    ! get missing_value
    status = nf90_get_att (ncid,varid,'missing_value',mdi)
    if (status /= nf90_noerr) mdi = -1E30
    
    ! get lowest level /= mdi
    do ilon = 1, nlon
       do ilat = 1, nlat
          ilev = 1
          do while (ABS((press(ilon,ilat,ilev,1)-mdi)/mdi)<=eps .and. ilev<nlev) 
             ilev = ilev + 1 
          enddo
          lev0(ilon,ilat) = lev(ilev)
       enddo
    enddo

    ! get grid information
    status = nf90_inq_varid (ncid,'lat',varid) 
    call nc_check_error (status,'Cannot find variable lat',abort=.false.)
    if (status/=0) return
    status = nf90_get_var (ncid,varid,lat0,start=(/1/),count=(/1/))
    call nc_check_error (status,'Cannot read variable lat',abort=.false.)
    if (status/=0) return
    status = nf90_get_var (ncid,varid,lat1,start=(/2/),count=(/1/))
    call nc_check_error (status,'Cannot read variable lat',abort=.false.)
    if (status/=0) return
    status = nf90_inq_varid (ncid,'lon',varid)
    call nc_check_error (status,'Cannot find variable lon',abort=.false.)
    if (status/=0) return
    status = nf90_get_var (ncid,varid,lon0,start=(/1/),count=(/1/))
    call nc_check_error (status,'Cannot read variable lon',abort=.false.)
    if (status/=0) return
    status = nf90_get_var (ncid,varid,lon1,start=(/2/),count=(/1/))
    call nc_check_error (status,'Cannot read variable lon',abort=.false.)
    if (status/=0) return

    latlon_grid%lat0 = lat0(1)
    latlon_grid%lon0 = lon0(1)
    latlon_grid%dlat = lat1(1) - lat0(1)
    latlon_grid%dlon = lon1(1) - lon0(1)
    latlon_grid%nlat = nlat
    latlon_grid%nlon = nlon
 
    ! ueberpruefen, ob unterste Schicht ueber lev_in_down hinausragt !!
    if (maxval(lev0) > lev_in_down) then
       write (*,*) 'Warning: lev0 contains values greater than lev_in_down!'
    elseif (maxval(lev0)+delta_lev > lev_in_down) then
       write (*,*) 'Warning: lev0+delta_lev contains values greater than lev_in_down!'
    endif
    
    ! deallocate arrays
    deallocate (lev, press)
    
    ! close file
    status = nf90_close (ncid)

  END SUBROUTINE get_orography

  !****************************************************************************
  ! Delete all APs on lower/upper/north/south boundaries
  !****************************************************************************
  SUBROUTINE delete_boundaries (lat,lon,lev,state,specarr)

    USE messy_clams_global,     ONLY: rank, dnparts, dnparts_max_shuffle, &
                                      prec, mdi, species_type, nspec
    USE messy_clamsbmix_global, ONLY: replace_low, replace_up, replace_north, &
                                      replace_south, latlon_grid, &
                                      lev_in_down, lev_in_up, lat_in_down, lat_in_up, &
                                      delta_lev, lev0
    use messy_clamsbmix_tools,  ONLY: pack_array

    implicit none

    REAL(PREC),         DIMENSION(:), POINTER :: lat,lon,lev,state
    TYPE(species_type), DIMENSION(:), POINTER :: specarr

    logical,            dimension(:), pointer :: mask

    integer   :: ipart, ilat, ilon, ispec

    allocate (mask(size(lat)))

    !write (*,*) rank,' vor loeschen: dnparts=',dnparts

    !-------------------------------------------------------------------
    ! Delete lower boundary
    !-------------------------------------------------------------------
    if (replace_low) then

       mask = .true.
       if (size(mask)>dnparts) mask(dnparts+1:) = .false.
    
       if (delta_lev==0) then  ! without orography

          ! mark all points with lev<=lev_in_down
          do ipart = 1, dnparts
             if (lev(ipart)<=lev_in_down) mask(ipart)=.false.
          enddo
          
       else  ! with orography

          ! mark all points with lev<=lev0+delta_lev
          !                 and  lev<=lev_in_down
          do ipart = 1, dnparts
     
             ! get position of point in latlon_grid (lat,lon)
             ilat = int((lat(ipart)-latlon_grid%lat0)/latlon_grid%dlat) + 1
! op_pj_20160606+
!!$             ilon = int((modulo(lon(ipart),360.)-latlon_grid%lon0)/latlon_grid%dlon) + 1
             ilon = int((modulo(lon(ipart),360._prec)-latlon_grid%lon0)/latlon_grid%dlon) + 1
! op_pj_20160606-
             
             if (lev(ipart)<=lev_in_down .and. &
                 lev(ipart)<=lev0(ilon,ilat)+delta_lev) mask(ipart) = .false.
                          
          enddo

       endif

       ! delete all marked points
       call pack_array (lat, mask)
       call pack_array (lon, mask)
       call pack_array (lev, mask)
       call pack_array (state, mask)
       do ispec = 1, nspec
          call pack_array (specarr(ispec)%values, mask)
       enddo

       ! get new number of points
       dnparts = count(mask)

       !write (*,*) rank,' nach loeschen (lower bound): dnparts=',dnparts

    endif

    !-------------------------------------------------------------------
    ! Delete upper boundary
    !-------------------------------------------------------------------
    if (replace_up) then

       mask = .true.
       if (size(mask)>dnparts) mask(dnparts+1:) = .false.

       ! mark all points with lev>=lev_in_up
       do ipart = 1, dnparts
          if (lev(ipart)>=lev_in_up) mask(ipart)=.false.
       enddo

       ! delete all marked points
       call pack_array (lat, mask)
       call pack_array (lon, mask)
       call pack_array (lev, mask)
       call pack_array (state, mask)
       do ispec = 1, nspec
          call pack_array (specarr(ispec)%values, mask)
       enddo

       ! get new number of points
       dnparts = count(mask)

       !write (*,*) rank,' nach loeschen (upper bound): dnparts=',dnparts

    endif

    !-------------------------------------------------------------------
    ! Delete north and south boundaries
    !-------------------------------------------------------------------
    mask = .true.
    if (size(mask)>dnparts) mask(dnparts+1:) = .false.

      ! mark all points with lat>=lat_in_up
    if (replace_north) then
       do ipart = 1, dnparts
          if (lev(ipart)>=lat_in_up) mask(ipart)=.false.
       enddo
    endif

    ! mark all points with lat<=lat_in_down
    if (replace_south) then
       do ipart = 1, dnparts
          if (lev(ipart)<=lat_in_down) mask(ipart)=.false.
       enddo
    endif
    
    if (replace_north .or. replace_south) then

       ! delete all marked points
       call pack_array (lat, mask)
       call pack_array (lon, mask)
       call pack_array (lev, mask)
       call pack_array (state, mask)
       do ispec = 1, nspec
          call pack_array (specarr(ispec)%values, mask)
       enddo
       
       ! get new number of points
       dnparts = count(mask)

       !write (*,*) rank,' nach loeschen (north, south): dnparts=',dnparts

    endif

    !-------------------------------------------------------------------
    ! Set missing values  
    !-------------------------------------------------------------------
    lat(dnparts+1:dnparts_max_shuffle) = mdi
    lon(dnparts+1:dnparts_max_shuffle) = mdi
    lev(dnparts+1:dnparts_max_shuffle) = mdi
    do ispec = 1, nspec
       specarr(ispec)%values(dnparts+1:dnparts_max_shuffle) = mdi
    enddo

    !-------------------------------------------------------------------
    ! clean up
    !-------------------------------------------------------------------
    deallocate (mask)

  END SUBROUTINE delete_boundaries

  !****************************************************************************
  ! replace lower boundary layer
  !****************************************************************************
  SUBROUTINE get_lower_bound (status, file_theta, lat, lon, lev, &
                              state, specarr,  &
                              lat_low, lon_low, lev_low, specarr_low, lcoupled)

    USE messy_clams_global,       ONLY: YEAR_NEXT   
    USE messy_clams_global,       ONLY: rank, dnparts, dnparts_max_shuffle, &
                                        prec, species_type, nspec, &
                                        filenamelen, specnamelen
    USE messy_clamsmix_global,    ONLY: adapt_par
    USE messy_clamsbmix_global,   ONLY: dir_boundfiles, nbmixbounds, &
                                        bmixbound, file_bounds, &
                                        interpol_from_init, delta_lev, nlevs, &
                                        lev_down, lev_in_down, lev_in_up, &
                                        lat_in_down, lat_in_up, time_init_value
    USE messy_clamsbmix_tools,    ONLY: get_boundfile_indices, &
                                        read_boundfile, get_bound_species, &
                                        interpolate_species
    USE messy_clams_tools_eqlatutils, ONLY: get_eqlat_int
    use messy_clams_tools_utils,   only: uppercase
    
    implicit none

    integer                                   :: status
    CHARACTER(*)                              :: file_theta
    REAL(PREC),         DIMENSION(:), POINTER :: lat, lon, lev, state, &
                                                 lat_low, lon_low, lev_low
    TYPE(species_type), DIMENSION(:), POINTER :: specarr, specarr_low
    LOGICAL                                   :: lcoupled

    ! above/below this level eqlat/lat will be used for the horizontal boundaries
    real,    parameter        :: lev_eqlat=340.

    REAL(PREC), dimension(:), pointer :: eqlat, helplev
    integer,    dimension(:), pointer :: indices
    
    integer :: dnparts_old, dnparts_new
    integer :: ispec, i, ibound
    
    character(filenamelen) :: filename
    character(specnamelen), dimension(:), pointer  :: speclist

    status = 0 ! no error

    ! Get indices of points 
    if (delta_lev>0) then
       call get_boundfile_indices (status, file_bounds, lat_in_down, lat_in_up, &
            lev_down, lev_in_down, indices, use_oro=.true.)     
    else
       call get_boundfile_indices (status, file_bounds, lat_in_down, lat_in_up, &
            lev_down, lev_in_down, indices)     
    endif
    if (status/=0) then
       write (*,*) 'Sub. get_lower_bound: Error in get_boundfile_indices!!!'
       return
    endif
    
    dnparts_old = dnparts

    dnparts_new = size(indices)
    
    !write (*,*) rank,'get_lower_bound: dnparts_new=', dnparts_new

    ! Ueberpruefen, ob Anz. Punkte auf rank passt: 
    ! dnparts+dnparts_new<=dnparts_max_shuffle
    if (dnparts+dnparts_new>dnparts_max_shuffle) then
       write (*,*) 
       write (*,*) rank, 'Number of points greater than ',dnparts_max_shuffle
       write (*,*) rank, 'dnparts, dnparts_new = ', dnparts, dnparts_new
       write (*,*) 
       status = 302
       return
    endif

    if (dnparts_new > 0) then
       
       ! Haenge eingelesene Punkte an bisherige Felder an: 
       ! lat/lon/lev(dnparts_old+1:dnparts_old+dnparts_new)
       ! dnparts = dnparts_old+dnparts_new

       ! Read points from boundary file and add to arrays
       ! => dnparts CHANGED !
       call read_boundfile (status, file_bounds, indices,  &
                               lat, lon, lev, specarr)
       if (status/=0) then
          write (*,*) 'Sub. get_lower_bound: Error in read_boundfile!!!'
          return
       endif 

       ! Set eqlat to lat for the lower boundary
       allocate (eqlat(size(lat)))
       eqlat = lat

       if (.not. lcoupled) then
          if (delta_lev==0) then   ! without orography
             ! Get EQLAT 
             if (lev_in_down > lev_eqlat) then
                if (rank==0) write (*,*) 'file_theta: ',trim(file_theta)
                allocate (helplev(size(lat)))
                helplev = adapt_par%lev_grid(1)
                call get_eqlat_int (status, trim(file_theta), lat, lon, helplev, eqlat, &
                                    dnparts_old+1, dnparts_old+dnparts_new)
                !write (79,*) eqlat(dnparts_old+1:dnparts_old+dnparts_new)
                deallocate (helplev)
                if (status/=0) then
                   write (*,*) 'Error in get_eqlat_int!!!'
                   return
                endif
             endif
          endif
       endif

       ! Interpolation nur fuer die neu hinzugefuegten Punkte aus Boundfile !!!
       ! Es werden die vollstaendigen Felder uebergeben und nur die Positionen
       ! von dnparts_old+1 bis dnparts_old+dnparts_new bearbeitet !
       
       ! interpolate species from init file
       if (interpol_from_init==1 .or. interpol_from_init==3) then
          if (rank==0) write (*,*) 'interpolate lower boundary from init file'
          call interpolate_species (status, lat, lon, lev, specarr,  &
               lat_low, lon_low, lev_low, specarr_low, &
               dnparts_old+1, dnparts_old+dnparts_new)
          if (status/=0) then
             write (*,*) 'Sub. get_lower_bound: Error in interpolate_species!!!'
             return
          endif
       endif
       

       ! Species nur fuer die neu hinzugefuegten Punkte aus Boundfiles ersetzen !!!
       ! Es werden die vollstaendigen Felder uebergeben und nur die Positionen
       ! von dnparts_old+1 bis dnparts_old+dnparts_new bearbeitet !

       ! redefine species from boundary files 
       do ibound = 1, nbmixbounds
          if (bmixbound(ibound)%no==1) then  ! redefine lower boundary

             IF (bmixbound(ibound)%startyear<=YEAR_NEXT .AND. YEAR_NEXT<=bmixbound(ibound)%endyear) THEN

                if (uppercase(bmixbound(ibound)%spec)=='ALL') then
                   allocate (speclist(nspec))
                   do ispec = 1, nspec
                      speclist(ispec) = specarr(ispec)%name
                   enddo
                else
                   allocate (speclist(1))
                   speclist(1) =  bmixbound(ibound)%spec
                endif
                
                call replace_yyyy_in_filename (bmixbound(ibound)%file, filename)
                filename = trim(dir_boundfiles)//'/'//trim(filename)
                
                call get_bound_species(status, lat, lon, eqlat, specarr, speclist, &
                                       dnparts_old+1, dnparts_old+dnparts_new, &
                                       filename, time_init_value, bmixbound(ibound)%add)
                if (status/=0) then
                   write (*,*) 'Sub. get_lower_bound: Error in get_bound_species!!!'
                   return
                endif
                
                deallocate (speclist)

             ENDIF
                
          end if
       end do

       ! set state to 10 
       do i = dnparts_old+1, dnparts_old+dnparts_new
          state(i) = 10
       end do
    
       deallocate (eqlat)
  
    endif

    deallocate (indices)

    !write (*,*) rank,'get_lower_bound: dnparts=', dnparts

  END SUBROUTINE get_lower_bound

  !****************************************************************************
  ! replace upper boundary layer 
  !****************************************************************************
  subroutine get_upper_bound (status,file_theta,lat,lon,lev,state,specarr, &
                              lat_upper,lon_upper,lev_upper,specarr_upper,lcoupled)

    USE messy_clams_global,     ONLY: YEAR_NEXT   
    USE messy_clams_global,     ONLY: dnparts, dnparts_max_shuffle, rank, &
                                      prec, species_type, nspec, &
                                      filenamelen, specnamelen
    USE messy_clamsmix_global,  ONLY: adapt_par
    USE messy_clamsbmix_global, ONLY: lev_up, nbmixbounds, dir_boundfiles, &
                                      bmixbound, file_bounds, interpol_from_init, &
                                      lat_in_down, lat_in_up, &
                                      lev_in_down, lev_in_up, nlevs, &
                                      time_init_value
    USE messy_clamsbmix_tools,  ONLY: get_boundfile_indices, &
                                      read_boundfile, get_bound_species, &
                                      interpolate_species
    USE messy_clams_tools_eqlatutils, ONLY: get_eqlat_int
    use messy_clams_tools_utils,   only: uppercase

    implicit none

    integer                                   :: status
    CHARACTER(*)                              :: file_theta
    REAL(PREC),         DIMENSION(:), POINTER :: lat, lon, lev, state
    REAL(PREC),         DIMENSION(:), POINTER :: lat_upper, lon_upper, lev_upper
    TYPE(species_type), DIMENSION(:), POINTER :: specarr, specarr_upper
    LOGICAL                                   :: lcoupled

    REAL(PREC), dimension(:), pointer :: eqlat, helplev
    integer,    dimension(:), pointer :: indices
    
    integer :: dnparts_old, dnparts_new
    integer :: ispec, i, ibound

    character(filenamelen) :: filename
    character(specnamelen), dimension(:), pointer  :: speclist

    status = 0

    ! Get indices of points 
    call get_boundfile_indices (status, file_bounds, lat_in_down, lat_in_up, &
                               lev_in_up, lev_up, indices)     
    if (status/=0) then
       write (*,*) 'Sub. get_upper_bound: Error in get_boundfile_indices!!!'
       return
    endif
    
    dnparts_old = dnparts

    dnparts_new = size(indices)

    !write (*,*) rank, 'get_upper_bound: dnparts_new=', dnparts_new
    
    ! Ueberpruefen, ob Anz. Punkte auf rank passt: 
    ! dnparts+dnparts_new<=dnparts_max_shuffle
    if (dnparts+dnparts_new>dnparts_max_shuffle) then
       write (*,*) 
       write (*,*) rank, 'Number of points greater than ',dnparts_max_shuffle
       write (*,*) rank, 'dnparts, dnparts_new = ', dnparts, dnparts_new
       write (*,*) 
       status = 302
       return
    endif

    if (dnparts_new > 0) then
 
       ! Haenge eingelesenen Punkte an bisherige Felder an: 
       ! lat/lon/lev(dnparts_old+1:dnparts_old+dnparts_new)
       ! dnparts = dnparts_old+dnparts_new

       ! Read points from boundary file and add to arrays
       ! => dnparts CHANGED !
       call read_boundfile (status, file_bounds, indices,  &
                            lat, lon, lev, specarr)
       if (status/=0) then
          write (*,*) 'Sub. get_upper_bound: Error in read_boundfile!!!'
          return
       endif

       ! Set eqlat to lat for the upper boundary
       allocate (eqlat(size(lat)))
       eqlat = lat

       
       if (.not. lcoupled) then       
          ! get EQLAT
          if (rank==0) write (*,*) 'file_theta: ',trim(file_theta)
          allocate (helplev(size(lat)))
          helplev = adapt_par%lev_grid(nlevs)
          call get_eqlat_int (status, trim(file_theta), lat, lon, helplev, eqlat, &
                              dnparts_old+1, dnparts_old+dnparts_new)
          deallocate (helplev)
          if (status/=0) then
             write (*,*) 'Error in get_eqlat_int!!!'
             return
          endif
       endif
       
       ! Interpolation nur fuer die neu hinzugefuegten Punkte aus Boundfile !!!
       ! Es werden die vollstaendigen Felder uebergeben und nur die Positionen
       ! von dnparts_old+1 bis dnparts_old+dnparts_new bearbeitet !
       
       ! interpolate species from init file
       if (interpol_from_init==2 .or. interpol_from_init==3) then
          
          if (rank==0) write (*,*) 'interpolate upper boundary from init file'
          call interpolate_species (status, lat, lon, lev, specarr,  &
               lat_upper, lon_upper, lev_upper, specarr_upper, &
               dnparts_old+1, dnparts_old+dnparts_new)
          if (status/=0) then
             write (*,*) 'Sub. get_upper_bound: Error in interpolate_species!!!'
             return
          endif

       endif

       ! Species nur fuer die neu hinzugefuegten Punkte aus Boundfiles ersetzen !!!
       ! Es werden die vollstaendigen Felder uebergeben und nur die Positionen
       ! von dnparts_old+1 bis dnparts_old+dnparts_new bearbeitet !
       
       ! redefine species from boundary files 
       do ibound = 1, nbmixbounds
          if (bmixbound(ibound)%no==2) then  ! redefine upper boundary

             IF (bmixbound(ibound)%startyear<=YEAR_NEXT .AND. YEAR_NEXT<=bmixbound(ibound)%endyear) THEN

                if (uppercase(bmixbound(ibound)%spec)=='ALL') then
                   allocate (speclist(nspec))
                   do ispec = 1, nspec
                      speclist(ispec) = specarr(ispec)%name
                   enddo
                else
                   allocate (speclist(1))
                   speclist(1) =  bmixbound(ibound)%spec
                endif
                
                call replace_yyyy_in_filename (bmixbound(ibound)%file, filename)
                filename = trim(dir_boundfiles)//'/'//trim(filename)
                
                call get_bound_species(status, lat, lon, eqlat, specarr, speclist, &
                                       dnparts_old+1, dnparts_old+dnparts_new, &
                                       filename, time_init_value, bmixbound(ibound)%add)
                if (status/=0) then
                   write (*,*) 'Sub. get_upper_bound: Error in get_bound_species!!!'
                   return
                endif
                
                deallocate (speclist)

             ENDIF
             
          end if
       end do
       
       ! set state to 11 
       do i = dnparts_old+1, dnparts_old+dnparts_new
          state(i) = 11
       end do
    
       deallocate (eqlat)

    endif

    deallocate (indices)

    !write (*,*) rank,'get_upper_bound: dnparts=', dnparts

  end subroutine get_upper_bound

  !****************************************************************************
  ! replace northern or southern vertical layer
  !****************************************************************************
  subroutine get_vertical_bound (status,lat,lon,lev,specarr,north)

    USE messy_clams_global,     ONLY: YEAR_NEXT   
    USE messy_clams_global,     ONLY: dnparts, dnparts_max_shuffle, rank, &
                                      prec, species_type, nspec, &
                                      filenamelen, specnamelen
    USE messy_clamsbmix_global, ONLY: lat_up, lat_down, lev_up, lev_down, &
                                      nbmixbounds, bmixbound, &
                                      dir_boundfiles, file_bounds, &
                                      lat_in_down, lat_in_up, time_init_value
    USE messy_clamsbmix_tools,  ONLY: get_boundfile_indices, &
                                      read_boundfile, get_vert_species
    use messy_clams_tools_utils,only: uppercase

    implicit none

    integer :: status
    logical :: north
    REAL(PREC),         DIMENSION(:), POINTER :: lat, lon, lev
    TYPE(species_type), DIMENSION(:), POINTER :: specarr

    integer, dimension(:), pointer :: indices
    
    integer :: dnparts_old, dnparts_new
    integer :: ispec, ibound
    real    :: latitude

    character(specnamelen), dimension(:), pointer  :: speclist
    character(filenamelen) :: filename
    character(2)   :: latstr

    status = 0

    ! Get indices of points 
    if (north) then
       call get_boundfile_indices (status, file_bounds, lat_in_up, lat_up, &
                        lev_down, lev_up, indices)     
    else
       call get_boundfile_indices (status, file_bounds, lat_down, lat_in_down, &
                        lev_down, lev_up, indices)     
    endif
    if (status/=0) then
       write (*,*) 'Sub. get_vertical_bound: Error in get_boundfile_indices!!!'
       return
    endif
    
    dnparts_old = dnparts

    dnparts_new = size(indices)

    !write (*,*) rank,'get_vert_bound: dnparts_new=', dnparts_new
    
    ! Ueberpruefen, ob Anz. Punkte auf rank passt: 
    ! dnparts+dnparts_new<=dnparts_max_shuffle
    if (dnparts+dnparts_new>dnparts_max_shuffle) then
       write (*,*) 
       write (*,*) rank, 'Number of points greater than ',dnparts_max_shuffle
       write (*,*) rank, 'dnparts, dnparts_new = ', dnparts, dnparts_new
       write (*,*) 
       status = 302
       return
    endif

    if (dnparts_new > 0) then
       
       ! Haenge eingelesenen Punkte an bisherige Felder an: 
       ! lat/lon/lev(dnparts_old+1:dnparts_old+dnparts_new)
       ! dnparts = dnparts_old+dnparts_new

       ! Read points from boundary file and add to arrays
       call read_boundfile (status, file_bounds, indices,  &
                               lat, lon, lev, specarr)
       if (status/=0) then
          write (*,*) 'Sub. get_vertical_bound: Error in read_boundfile!!!'
          return
       endif

       ! Species nur fuer die neu hinzugefuegten Punkte aus Boundfile ersetzen !!!
       ! Es werden die vollstaendigen Felder uebergeben und nur die Positionen
       ! von dnparts_old+1 bis dnparts_old+dnparts_new bearbeitet !

       ! redefine species from boundary files 
       do ibound = 1, nbmixbounds

          if ((north .and. bmixbound(ibound)%no==3) .or. &
              (.not. north .and. bmixbound(ibound)%no==4)) then  ! redefine vertical boundary

             IF (bmixbound(ibound)%startyear<=YEAR_NEXT .AND. YEAR_NEXT<=bmixbound(ibound)%endyear) THEN

                if (uppercase(bmixbound(ibound)%spec)=='ALL') then
                   allocate (speclist(nspec))
                   do ispec = 1, nspec
                      speclist(ispec) = specarr(ispec)%name
                   enddo
                else
                   allocate (speclist(1))
                   speclist(1) =  bmixbound(ibound)%spec
                endif
                
                call replace_yyyy_in_filename (bmixbound(ibound)%file, filename)
                filename = trim(dir_boundfiles)//'/'//trim(filename)
                
                call get_vert_species (status, lev, lon, specarr, speclist, &
                                       dnparts_old+1, dnparts_old+dnparts_new, &
                                       filename, time_init_value)
                if (status/=0) then
                   write (*,*) 'Sub. get_vertical_bound: Error in get_vert_species!!!'
                   return
                endif
                
                deallocate (speclist)

             ENDIF
                
          end if
       end do


    endif

    deallocate (indices)

    !write (*,*) rank,'get_vert_bound: dnparts=', dnparts

  end subroutine get_vertical_bound

END MODULE messy_clamsbmix

program TEM_Diagnostics
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
! RW, January 2019
! 
! This program calculates:
!      - latitudinal and vertical EP flux                (Fy,Fz)
!      - latitudinal contribution to accel_divf       (accel_divf_lat)
!      - vertical    contribution to accel_divf       (accel_divf_vert)
!      ( Thus, it holds: accel_divf = accel_divf_lat + accel_divf_vert )
!      - and the residual mean meridional circulation (vstar,wstar)
! in latitude-pressure coordinates for every input time step.
!
! The calculations are based on Andrwews, Mahlman, Sinclair, 1983
! ( https://journals.ametsoc.org/doi/abs/10.1175/1520-0469%281983%29040%3C2768%3AETWATM%3E2.0.CO%3B2 )
! with equations (2.1) - (2.3) for EP flux
! and  equation  (2.6)         for residual mean meridional circulation
!
! CALLING the program:
! './TEM_Diagnostics_with_EPFD_split.out input_file.nc output_file.nc'
!
!
! BEFORE USE:
! 1. Assign variable names: lon_name, lat_name, lev_name, time_name,
!                           temp_name, u_name, v_name, w_name
! 2. Choose percentage of data that has to exist for zonal averaging: zonal_mean_threshold
!
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!
USE messy_main_constants_mem, ONLY: DP, FLAGGED_BAD
USE messy_main_tools, ONLY: read_nml_open, read_nml_check, read_nml_close
USE netcdf
#if defined (NAGFOR)
USE f90_unix_env, ONLY: iargc, getarg
#endif
       
implicit none
intrinsic      :: trim, maxval, real


! FIRST OF ALL! Specify names of grid dimensions and variables
        integer,                parameter   :: name_len  = 10     ! maximum number of characters for grid and variable names
        
        
        !************************************************************
        character(len=name_len)  :: lon_name  = 'lon'
        character(len=name_len)  :: lat_name  = 'lat'
        character(len=name_len)  :: lev_name  = 'plev'
        character(len=name_len)  :: time_name = 'time'
        
        character(len=name_len)  :: temp_name = 'tm1'
        character(len=name_len)  :: u_name    = 'um1'
        character(len=name_len)  :: v_name    = 'vm1'
        character(len=name_len)  :: w_name    = 'vervel'
        !************************************************************

        NAMELIST /CTRL/ lon_name, lat_name, lev_name, time_name &
             , temp_name, u_name, v_name, w_name
        
        integer                             :: plev_units_length  ! length of units-attribute of plev
        character(len=:),allocatable        :: plev_units         ! to check pressure level units

! input parameter of TEM_Diagnostics program
!!$        integer*4                           :: iargc
        character(len=90)                   :: infile,outfile

! status of netcdf functions
        integer                             :: status

! netcdf input file variables.
        integer                             :: nx           ! number of longitudes
        integer                             :: ny           ! number of latitudes
        integer                             :: nl           ! number of pressure levels
        integer                             :: nt           ! number of time values in netcdf file
        
        integer                             :: ncid_in      ! netcdf file identity (ID) of input file
        integer                             :: ncid_out     ! ID of output file
        
        integer                             :: lonid        ! ID of longitudes
        integer                             :: latid        ! ID of latitudes
        integer                             :: levid        ! ID of pressure levels
        integer                             :: timeid       ! ID of times
        
        integer                             :: uid, vid, wid   ! IDs for wind components
        integer                             :: tempid          ! ID  for temperature
        
        ! input data
        real(dp)  ,allocatable,dimension(:)       :: lon,lat,plev,time
        real(dp),allocatable,dimension(:)       :: time8
        real(dp)  ,allocatable,dimension(:,:,:,:) :: Temp,Uwind,Vwind,Wwind
        ! for flipping dimension
        real(dp),allocatable,dimension(:)       :: coordinate_flipped
        real(dp),allocatable,dimension(:,:,:,:) :: Tflipped,Uflipped,Vflipped,Wflipped

! variables for copying dimensions
        integer                             :: latvinid     ! ID of lat variable of input file
        integer                             :: latvoutid    ! ID of lat variable of output file
        integer                             :: levvinid     ! ID of lev variable of input file
        integer                             :: levvoutid    ! ID of lev variable of output file
        integer                             :: timevinid    ! ID of time variable of input file
        integer                             :: timevoutid   ! ID of time variable of output file
        
        integer                             :: varnatts     ! number of variable attributes
        integer                             :: j_varnatts   ! loop variable for varnatts
        character(len=nf90_max_name)          :: name         ! name of attribute of variable
        
!netcdf output file variables
        integer                             :: accel_divf_lat_id  ! ID of lat. contribution to accel_divf
        integer                             :: accel_divf_vert_id ! ID of vert. contribution to accel_divf
        integer                             :: fyid          ! ID of latitudinal component of EP flux
        integer                             :: fzid          ! ID of vertical component of EP flux
        integer                             :: vstarid       ! ID of latitudinal component of residual circulation
        integer                             :: wstarid       ! ID of vertical component of residual circulation
        
        ! output data
        real(dp),allocatable,dimension(:,:,:)   :: accel_divf_lat,accel_divf_vert
        real(dp),allocatable,dimension(:,:,:)   :: fy,fz,vstar,wstar
        ! for flipping dimensions
        real(dp),allocatable,dimension(:,:,:)   :: accel_divf_lat_flipped,accel_divf_vert_flipped
        real(dp),allocatable,dimension(:,:,:)   :: fy_flipped,fz_flipped,vstar_flipped,wstar_flipped
        

! loop variables
        integer                             :: iy,il,it
        
! logical variables
        ! do not change these logicals here! they are just initialized here and set to their correct values later!
        logical                             :: plev_in_Pa      = .true.   ! are pressure levels units Pa or hPa
        logical                             :: plev_increasing = .true.   ! are pressure levels increasingly ordered
        logical                             :: lat_increasing  = .true.   ! are latitudes increasingly ordered

! for reading the namelist file
        logical                             :: lex ! namelist file exists?
        integer                             :: iou = 23 ! for namelist read
        integer                             :: fstat ! read status
        
!!$! Introduce a number for missing values which should be < 0.0 and ABS(.) > 1.e10
!!$        real,parameter                  :: MISSING_VALUE = -1.0E+34
        real :: MISSING_VALUE

!==================================================================================================
! Main program
! - open and read input nc file
! - convert pressure levels from hPa to Pa if needed
! - transform pressure and latitude to decreasing order if needed
! - pass data to TEM diagnostics subroutine
! - establish original order of pressure and latitude if needed
! - establish original pressure level units if needed
! - write data to output nc file

MISSING_VALUE=REAL(FLAGGED_BAD)

! get input from command line
if (iargc().ne.2) then
 write(*,*)" "
 write(*,*)"TEM_Diagnostics_with_EPFD_split: number of command line arguments incorrect"
 write(*,*)"usage:"
 write(*,*)"   temdiag.exe input_file output_file"
 write(*,*)"abort"
 stop
end if

  write(*,*)' ' 
  write(*,*)'TEM_Diagnostics_with_EPFD_split program called'

  call getarg(1,infile)
  write(*,*)'reading from input file: ',trim(infile)

  call getarg(2,outfile)
  write(*,*)'writing to output file: ',trim(outfile)

  CALL read_nml_open(lex, 'main', iou, 'CTRL', 'temdiag')
  IF (.not.lex) STOP  
  READ(iou, NML=CTRL, IOSTAT=fstat)
  CALL read_nml_check(fstat, 'main', iou, 'CTRL', 'temdiag')
  IF (fstat /= 0) STOP  ! error while reading namelist
  CALL read_nml_close('main', iou, 'temdiag')
  
  write(*,*)'names for dimensions: ',trim(lon_name),', ',trim(lat_name),', ',trim(lev_name),', ',trim(time_name)
  write(*,*)'names for variables: ',trim(temp_name),', ',trim(u_name),', ',trim(v_name),', ',trim(w_name)

!==================================================================================================
! READ INFILE nc DATA

        ! read in lat, lon, levels, time, u,v,w and temperature
        write(*,*) ' '
        write(*,*) '1. Reading from input file: ',trim(infile)

        ! open input netcdf file (nf90_nowrite means read-only access, ncid_in is the ID of the input file)
        status = nf90_open(trim(infile),nf90_nowrite,ncid_in)
        if(status /= nf90_noerr) call handle_nf90_err(status)
        
! READ DIMENSIONS OF THE GRID AND ITS VALUES
        ! times
        status = nf90_inq_dimid(ncid_in,trim(time_name),timeid)   ! return the ID of the dimension time
        if (status /= nf90_noerr) call handle_nf90_err(status)
        status = nf90_inquire_dimension(ncid_in,timeid,len=nt)               ! return the length of the dimension time
        if (status /= nf90_noerr) call handle_nf90_err(status)
        allocate(time(nt))
        allocate(time8(nt))

        status = nf90_inq_varid(ncid_in,trim(time_name),timeid)   ! return the ID of the variable time
        if (status /= nf90_noerr) call handle_nf90_err(status)
        status = nf90_get_var(ncid_in,timeid,time8)        ! read all variable values
        if (status /= nf90_noerr) call handle_nf90_err(status)
!!$        time(:) = real(time8(:))                                ! convert from double to real
        time(:) = time8(:)                                ! copy
        
        !do it=1,nt
        !        write(*,*)'Time variable real',time(it),' and double ',time8(it)
        !end do
        
        deallocate(time8)
        

        !**************************************************************
        ! pressure levels
        status = nf90_inq_dimid(ncid_in,trim(lev_name),levid)
        if (status /= nf90_noerr) call handle_nf90_err(status)
        status = nf90_inquire_dimension(ncid_in,levid,len=nl)
        if (status /= nf90_noerr) call handle_nf90_err(status)
        allocate(plev(nl))

        status = nf90_inq_varid(ncid_in,trim(lev_name),levid)
        if (status /= nf90_noerr) call handle_nf90_err(status)
        status = nf90_get_var(ncid_in,levid,plev)
        if (status /= nf90_noerr) call handle_nf90_err(status)
        
        ! pay attention to the units of the pressure levels
        ! 1. find out length of units-attribute of plev
        status = nf90_inquire_attribute(ncid_in,levid,'units',len=plev_units_length)
        if (status /= nf90_noerr) call handle_nf90_err(status)
        !write(*,*)'length of attribute units of plev: ',plev_units_length
        
        ! 2. allocate string with length of the pressure units
        allocate(character(plev_units_length) :: plev_units)
        
        ! 3. get units of plev
        status = nf90_get_att(ncid_in,levid,'units',plev_units)
        if (status /= nf90_noerr) call handle_nf90_err(status)
        
        ! 4. Transform pressure levels to Pa
        if( (plev_units .eq. 'hPa') .and. (maxval(plev(:)) < 1.0e4) ) then
                plev_in_Pa = .false.
                write(*,*)'Input pressure level units are hPa. Transform pressure level grid to Pa.'
                plev(:) = plev(:) * 100.0
        else if ( plev_units .eq. 'Pa' ) then
                plev_in_Pa = .true.
                write(*,*)'Input pressure level units are Pa. OK!'
        else
                write(*,*)'Unknown units of pressure levels!'
                write(*,*)'abort'
                stop
        end if
        
        deallocate(plev_units)
        
        !do il=1,nl
        !        print*,'plev(',il,')=',plev(il)
        !enddo
        

        !**************************************************************
        ! latitudes
        status = nf90_inq_dimid(ncid_in,trim(lat_name),latid)
        if (status /= nf90_noerr) call handle_nf90_err(status)
        status = nf90_inquire_dimension(ncid_in,latid,len=ny)
        if (status /= nf90_noerr) call handle_nf90_err(status)
        allocate(lat(ny))

        status = nf90_inq_varid(ncid_in,trim(lat_name),latid)
        if (status /= nf90_noerr) call handle_nf90_err(status)
        status = nf90_get_var(ncid_in,latid,lat)
        if (status /= nf90_noerr) call handle_nf90_err(status)

        !**************************************************************
        ! longitudes
        status = nf90_inq_dimid(ncid_in,trim(lon_name),lonid)
        if (status /= nf90_noerr) call handle_nf90_err(status)
        status = nf90_inquire_dimension(ncid_in,lonid,len=nx)
        if (status /= nf90_noerr) call handle_nf90_err(status)
        allocate(lon(nx))
        
        status = nf90_inq_varid(ncid_in,trim(lon_name),lonid)
        if (status /= nf90_noerr) call handle_nf90_err(status)
        status = nf90_get_var(ncid_in,lonid,lon)
        if (status /= nf90_noerr) call handle_nf90_err(status)


! READ VARIABLE VALUES
        ! temperature [K]
        allocate(Temp(nx,ny,nl,nt))
        status = nf90_inq_varid(ncid_in,trim(temp_name),tempid)
        if (status /= nf90_noerr) call handle_nf90_err(status)
        ! read variable values by specifying a corner and a vector of edge lengths
        status = nf90_get_var(ncid_in,tempid,Temp,start=(/1,1,1,1/),count=(/nx,ny,nl,nt/))
        if (status /= nf90_noerr) call handle_nf90_err(status)


        ! u wind [m/s]
        allocate(Uwind(nx,ny,nl,nt))
        status = nf90_inq_varid(ncid_in,trim(u_name),uid)
        if (status /= nf90_noerr) call handle_nf90_err(status)
        status = nf90_get_var(ncid_in,uid,Uwind,start=(/1,1,1,1/),count=(/nx,ny,nl,nt/))
        if (status /= nf90_noerr) call handle_nf90_err(status)
        

        ! v wind [m/s]
        allocate(Vwind(nx,ny,nl,nt))
        status = nf90_inq_varid(ncid_in,trim(v_name),vid)
        if (status /= nf90_noerr) call handle_nf90_err(status)
        status = nf90_get_var(ncid_in,vid,Vwind,start=(/1,1,1,1/),count=(/nx,ny,nl,nt/))
        if (status /= nf90_noerr) call handle_nf90_err(status)


        ! vertical velocity [Pa/s] 
        allocate(Wwind(nx,ny,nl,nt))
        status = nf90_inq_varid(ncid_in,trim(w_name),wid)
        if (status /= nf90_noerr) call handle_nf90_err(status)
        status = nf90_get_var(ncid_in,wid,Wwind,start=(/1,1,1,1/),count=(/nx,ny,nl,nt/))
        if (status /= nf90_noerr) call handle_nf90_err(status)
        

! PAY ATTENTION TO THE ORDER OF PRESSURE LEVELS
! formulas for derivatives in subroutine calc_TEM_Diagnostics use decreasingly ordered pressure levels
        if ( plev(1) > plev(nl) ) then
                plev_increasing = .false.
                write(*,*)'Pressure levels decreasingly ordered. OK!'
        else if ( plev(1) < plev(nl) ) then
                plev_increasing = .true.
                write(*,*)'Pressure levels increasingly ordered.'
                write(*,*)'Flip pressure levels (and wind and temperature fields) to decreasing order.'
                
                allocate(coordinate_flipped(nl))
                allocate(Uflipped(nx,ny,nl,nt))
                allocate(Vflipped(nx,ny,nl,nt))
                allocate(Wflipped(nx,ny,nl,nt))
                allocate(Tflipped(nx,ny,nl,nt))
                
                do il=1,nl
                        coordinate_flipped(il)      = plev(nl + 1 - il)
                        Uflipped(1:nx,1:ny,il,1:nt) = Uwind(1:nx,1:ny,nl + 1 - il,1:nt)
                        Vflipped(1:nx,1:ny,il,1:nt) = Vwind(1:nx,1:ny,nl + 1 - il,1:nt)
                        Wflipped(1:nx,1:ny,il,1:nt) = Wwind(1:nx,1:ny,nl + 1 - il,1:nt)
                        Tflipped(1:nx,1:ny,il,1:nt) = Temp(1:nx,1:ny,nl + 1 - il,1:nt)
                end do
                
                do il=1,nl
                        plev(il)                 = coordinate_flipped(il)
                        Uwind(1:nx,1:ny,il,1:nt) = Uflipped(1:nx,1:ny,il,1:nt)
                        Vwind(1:nx,1:ny,il,1:nt) = Vflipped(1:nx,1:ny,il,1:nt)
                        Wwind(1:nx,1:ny,il,1:nt) = Wflipped(1:nx,1:ny,il,1:nt)
                        Temp(1:nx,1:ny,il,1:nt)  = Tflipped(1:nx,1:ny,il,1:nt)
                end do
                
                
                deallocate(Tflipped)
                deallocate(Wflipped)
                deallocate(Vflipped)
                deallocate(Uflipped)
                deallocate(coordinate_flipped)
        else
                write(*,*)'Something is wrong with the pressure levels!'
                write(*,*)'abort'
                stop
        end if
        
        !do il=1,nl
        !        print*,'plev(',il,')=',plev(il)
        !enddo
        
        


! PAY ATTENTION TO THE ORDER OF LATITUDES
! formulas for derivatives in subroutine calc_TEM_Diagnostics use decreasingly ordered latitudes
        if ( lat(1) > lat(ny) ) then
                lat_increasing = .false.
                write(*,*)'Latitudes decreasingly ordered. OK!'
        else if ( lat(1) < lat(ny) ) then
                lat_increasing = .true.
                write(*,*)'Latitudes increasingly ordered.'
                write(*,*)'Flip latitudes (and wind and temperature fields) to decreasing order.'
                
                allocate(coordinate_flipped(ny))
                allocate(Uflipped(nx,ny,nl,nt))
                allocate(Vflipped(nx,ny,nl,nt))
                allocate(Wflipped(nx,ny,nl,nt))
                allocate(Tflipped(nx,ny,nl,nt))
                
                do iy=1,ny
                        coordinate_flipped(iy)      = lat(ny + 1 - iy)
                        Uflipped(1:nx,iy,1:nl,1:nt) = Uwind(1:nx,ny + 1 - iy,1:nl,1:nt)
                        Vflipped(1:nx,iy,1:nl,1:nt) = Vwind(1:nx,ny + 1 - iy,1:nl,1:nt)
                        Wflipped(1:nx,iy,1:nl,1:nt) = Wwind(1:nx,ny + 1 - iy,1:nl,1:nt)
                        Tflipped(1:nx,iy,1:nl,1:nt) = Temp(1:nx,ny + 1 - iy,1:nl,1:nt)
                end do
                
                do iy=1,ny
                        lat(iy)                  = coordinate_flipped(iy)
                        Uwind(1:nx,iy,1:nl,1:nt) = Uflipped(1:nx,iy,1:nl,1:nt)
                        Vwind(1:nx,iy,1:nl,1:nt) = Vflipped(1:nx,iy,1:nl,1:nt)
                        Wwind(1:nx,iy,1:nl,1:nt) = Wflipped(1:nx,iy,1:nl,1:nt)
                        Temp(1:nx,iy,1:nl,1:nt)  = Tflipped(1:nx,iy,1:nl,1:nt)
                end do
                
                
                deallocate(Tflipped)
                deallocate(Wflipped)
                deallocate(Vflipped)
                deallocate(Uflipped)
                deallocate(coordinate_flipped)
        else
                write(*,*)'Something is wrong with the latitudes!'
                write(*,*)'abort'
                stop
        end if



! keep input file open , to read attributes into output nc file

        ! allocate output variables
        allocate(accel_divf_lat(ny,nl,nt))
        allocate(accel_divf_vert(ny,nl,nt))
        allocate(fy(ny,nl,nt))
        allocate(fz(ny,nl,nt))
        allocate(vstar(ny,nl,nt))
        allocate(wstar(ny,nl,nt))
        
        ! initialize output variables
        accel_divf_lat(:,:,:)  = 0.0
        accel_divf_vert(:,:,:) = 0.0
        fy(:,:,:)              = 0.0
        fz(:,:,:)              = 0.0
        vstar(:,:,:)           = 0.0
        wstar(:,:,:)           = 0.0
        
        
        
        
!==================================================================================================
        write(*,*) ' '
        write(*,*)'2. Calculating TEM Diagnostics: accel_divf_lat, accel_divf_vert, Fy, Fz, vstar, wstar'
        call calc_TEM_Diagnostics(nx,ny,nl,nt,lat,plev,Temp,Uwind,Vwind,Wwind &
                                    ,accel_divf_lat,accel_divf_vert,fy,fz,vstar,wstar)
        
!==================================================================================================
! BRING DIMENSIONS AND VARIABLES INTO ORIGINAL ORDER

        write(*,*)' '
        write(*,*)'3. Establish original order of pressure levels and latitudes (and output data).'
        if ( (.not. plev_increasing) .and. (.not. lat_increasing) ) then
                write(*,*)'Nothing has to be done. OK!'
        end if
        
        if ( plev_increasing ) then
                write(*,*)'Flip pressure levels (and output data) to increasing order again.'
                
                allocate(coordinate_flipped(nl))
                allocate(fy_flipped(ny,nl,nt))
                allocate(fz_flipped(ny,nl,nt))
                allocate(accel_divf_lat_flipped(ny,nl,nt))
                allocate(accel_divf_vert_flipped(ny,nl,nt))
                allocate(vstar_flipped(ny,nl,nt))
                allocate(wstar_flipped(ny,nl,nt))
                
                do il=1,nl
                        coordinate_flipped(il)                = plev(nl + 1 - il)
                        fy_flipped(1:ny,il,1:nt)              = fy(1:ny,nl + 1 - il,1:nt)
                        fz_flipped(1:ny,il,1:nt)              = fz(1:ny,nl + 1 - il,1:nt)
                        accel_divf_lat_flipped(1:ny,il,1:nt)  = accel_divf_lat(1:ny,nl + 1 - il,1:nt)
                        accel_divf_vert_flipped(1:ny,il,1:nt) = accel_divf_vert(1:ny,nl + 1 - il,1:nt)
                        vstar_flipped(1:ny,il,1:nt)           = vstar(1:ny,nl + 1 - il,1:nt)
                        wstar_flipped(1:ny,il,1:nt)           = wstar(1:ny,nl + 1 - il,1:nt)
                end do
                
                do il=1,nl
                        plev(il)                      = coordinate_flipped(il)
                        fy(1:ny,il,1:nt)              = fy_flipped(1:ny,il,1:nt)
                        fz(1:ny,il,1:nt)              = fz_flipped(1:ny,il,1:nt)
                        accel_divf_lat(1:ny,il,1:nt)  = accel_divf_lat_flipped(1:ny,il,1:nt)
                        accel_divf_vert(1:ny,il,1:nt) = accel_divf_vert_flipped(1:ny,il,1:nt)
                        vstar(1:ny,il,1:nt)           = vstar_flipped(1:ny,il,1:nt)
                        wstar(1:ny,il,1:nt)           = wstar_flipped(1:ny,il,1:nt)
                end do
                
                deallocate(wstar_flipped)
                deallocate(vstar_flipped)
                deallocate(accel_divf_vert_flipped)
                deallocate(accel_divf_lat_flipped)
                deallocate(fz_flipped)
                deallocate(fy_flipped)
                deallocate(coordinate_flipped)
        end if
        
        if ( lat_increasing ) then
                write(*,*)'Flip latitudes (and output data) to increasing order again.'
                
                allocate(coordinate_flipped(ny))
                allocate(fy_flipped(ny,nl,nt))
                allocate(fz_flipped(ny,nl,nt))
                allocate(accel_divf_lat_flipped(ny,nl,nt))
                allocate(accel_divf_vert_flipped(ny,nl,nt))
                allocate(vstar_flipped(ny,nl,nt))
                allocate(wstar_flipped(ny,nl,nt))
                
                do iy=1,ny
                        coordinate_flipped(iy)                = lat(ny + 1 - iy)
                        fy_flipped(iy,1:nl,1:nt)              = fy(ny + 1 - iy,1:nl,1:nt)
                        fz_flipped(iy,1:nl,1:nt)              = fz(ny + 1 - iy,1:nl,1:nt)
                        accel_divf_lat_flipped(iy,1:nl,1:nt)  = accel_divf_lat(ny + 1 - iy,1:nl,1:nt)
                        accel_divf_vert_flipped(iy,1:nl,1:nt) = accel_divf_vert(ny + 1 - iy,1:nl,1:nt)
                        vstar_flipped(iy,1:nl,1:nt)           = vstar(ny + 1 - iy,1:nl,1:nt)
                        wstar_flipped(iy,1:nl,1:nt)           = wstar(ny + 1 - iy,1:nl,1:nt)
                end do
                
                do iy=1,ny
                        lat(iy)                       = coordinate_flipped(iy)
                        fy(iy,1:nl,1:nt)              = fy_flipped(iy,1:nl,1:nt)
                        fz(iy,1:nl,1:nt)              = fz_flipped(iy,1:nl,1:nt)
                        accel_divf_lat(iy,1:nl,1:nt)  = accel_divf_lat_flipped(iy,1:nl,1:nt)
                        accel_divf_vert(iy,1:nl,1:nt) = accel_divf_vert_flipped(iy,1:nl,1:nt)
                        vstar(iy,1:nl,1:nt)           = vstar_flipped(iy,1:nl,1:nt)
                        wstar(iy,1:nl,1:nt)           = wstar_flipped(iy,1:nl,1:nt)
                end do
                
                deallocate(wstar_flipped)
                deallocate(vstar_flipped)
                deallocate(accel_divf_vert_flipped)
                deallocate(accel_divf_lat_flipped)
                deallocate(fz_flipped)
                deallocate(fy_flipped)
                deallocate(coordinate_flipped)
        end if

        
        
        ! Establish original pressure level units
        if( .not. plev_in_Pa ) then
                write(*,*)'Transform pressure level grid back to hPa.'
                plev(:) = plev(:) / 100.0
        end if
        
        
!==================================================================================================
! WRITE TEM DIAGNOSTICS CALCULATIONS TO OUTPUT FILE

        !output is for ny,nl,nt dimensions
        write(*,*) ' '
        write(*,*)'4. Writing data to output file: ',trim(outfile)

        ! open output netcdf file
        ! create new netcdf data file with creation mode nf90_clobber...
        ! ...specifying the default behavior: overwriting any existing dataset with the same file name
        status = nf90_create(trim(outfile),nf90_clobber,ncid_out) ! ncid_out is the ID of the output file
        if (status /= nf90_noerr) call handle_nf90_err(status)

! ADD DIMENSIONS TO OUTPUT FILE
        status = nf90_def_dim(ncid_out,lat_name,ny,latid)          ! add a new dimension to an open netCDF dataset
        if (status /= nf90_noerr) call handle_nf90_err(status)
        status = nf90_def_dim(ncid_out,lev_name,nl,levid)
        if (status /= nf90_noerr) call handle_nf90_err(status)
        status = nf90_def_dim(ncid_out,time_name,nf90_unlimited,timeid) 
        if (status /= nf90_noerr) call handle_nf90_err(status)

! COPY ATTRIBUTES OF GRID DIMENSIONS
        ! input file still open to copy attributes of dimensions
        
        ! latitude
        status = nf90_inq_varid(ncid_in,lat_name,latvinid)         ! give variable ID of lat of input file to latvinid
        if (status /= nf90_noerr) call handle_nf90_err(status)
        status = nf90_inquire_variable(ncid_in,latvinid,nAtts=varnatts)   ! give number of variable attributes of lat to varnatts
        if (status /= nf90_noerr) call handle_nf90_err(status)
        ! add the new variable lat to output file with variable type nf90_float and dimension 1
        ! latid is the variable ID of lat of the input file
        ! latvoutid is the returned variable ID of the output file
        status = nf90_def_var(ncid_out,lat_name,nf90_float,dimids = (/latid/), varid=latvoutid)
        if (status /= nf90_noerr) call handle_nf90_err(status)
        do j_varnatts=1,varnatts                              ! for every attribute of lat of input file...
          status = nf90_inq_attname(ncid_in,latvinid,j_varnatts,name)    ! ...give the attribute name...
          if (status /= nf90_noerr) call handle_nf90_err(status)
          status = nf90_copy_att(ncid_in,latvinid,trim(name),ncid_out,latvoutid)   ! ...and copy it to the output file
          if (status /= nf90_noerr) call handle_nf90_err(status)
        end do

       ! pressure level
        status = nf90_inq_varid(ncid_in,lev_name,levvinid)
        if (status /= nf90_noerr) call handle_nf90_err(status)
        status = nf90_inquire_variable(ncid_in,levvinid,nAtts=varnatts)
        if (status /= nf90_noerr) call handle_nf90_err(status)
        status = nf90_def_var(ncid_out,lev_name,nf90_float,dimids=(/levid/),varid=levvoutid)
        if (status /= nf90_noerr) call handle_nf90_err(status)
        do j_varnatts=1,varnatts
          status = nf90_inq_attname(ncid_in,levvinid,j_varnatts,name)
          if (status /= nf90_noerr) call handle_nf90_err(status)
          status = nf90_copy_att(ncid_in,levvinid,trim(name),ncid_out,levvoutid)
          if (status /= nf90_noerr) call handle_nf90_err(status)
        end do
        

! Here the pressure level unit in the output file can be changed
!**************************************************************************************************        
!        ! If pressure level unit was hPa in input file,
!        ! it has been copied to the output file.
!        ! Thus, change pressure level unit to Pa
!        if ( .not. plev_in_Pa ) then
!                ! 1. get plev ID
!                status = nf90_inq_varid(ncid_out,trim(lev_name),levvoutid)
!                if (status /= nf90_noerr) call handle_nf90_err(status)
!                ! 2. put new unit
!                status = nf90_put_att(ncid_out,levvoutid,'units','Pa')
!                if (status /= nf90_noerr) call handle_nf90_err(status)
!        end if
!**************************************************************************************************
        
        
        ! time
        status = nf90_inq_varid(ncid_in,time_name,timevinid)
        if (status /= nf90_noerr) call handle_nf90_err(status)
        status = nf90_inquire_variable(ncid_in,timevinid,nAtts=varnatts)
        if (status /= nf90_noerr) call handle_nf90_err(status)
        status = nf90_def_var(ncid_out,time_name,nf90_float,dimids=(/timeid/),varid=timevoutid)
        if (status /= nf90_noerr) call handle_nf90_err(status)
        do j_varnatts=1,varnatts
          status = nf90_inq_attname(ncid_in,timevinid,j_varnatts,name)
          if (status /= nf90_noerr) call handle_nf90_err(status)
          status = nf90_copy_att(ncid_in,timevinid,trim(name),ncid_out,timevoutid)
          if (status /= nf90_noerr) call handle_nf90_err(status)
        end do
        
! ADD DATA VARIABLES OF EP FLUX DIAGNOSTIC
        ! set up variable dimensions and attributes definitions
        
        ! accel_divf_lat definition
        name = 'accel_divf_lat'
        status = nf90_def_var(ncid_out,trim(name),nf90_float,dimids=(/latid,levid,timeid/),varid=accel_divf_lat_id)
        if (status /= nf90_noerr) call handle_nf90_err(status)
        status = nf90_put_att(ncid_out,accel_divf_lat_id,'long_name','Acceleration by latitudinal contribution to EP flux divergence')
        if (status /= nf90_noerr) call handle_nf90_err(status)
        status = nf90_put_att(ncid_out,accel_divf_lat_id,'units','m/s^2')
        if (status /= nf90_noerr) call handle_nf90_err(status)
        status = nf90_put_att(ncid_out,accel_divf_lat_id,'missing_value',MISSING_VALUE)
        if (status /= nf90_noerr) call handle_nf90_err(status)
        status = nf90_put_att(ncid_out,accel_divf_lat_id,'_FillValue',MISSING_VALUE)
        if (status /= nf90_noerr) call handle_nf90_err(status)
        
        ! accel_divf_vert definition
        name = 'accel_divf_vert'
        status = nf90_def_var(ncid_out,trim(name),nf90_float,dimids=(/latid,levid,timeid/),varid=accel_divf_vert_id)
        if (status /= nf90_noerr) call handle_nf90_err(status)
        status = nf90_put_att(ncid_out,accel_divf_vert_id,'long_name','Acceleration by vertical contribution to EP flux divergence')
        if (status /= nf90_noerr) call handle_nf90_err(status)
        status = nf90_put_att(ncid_out,accel_divf_vert_id,'units','m/s^2')
        if (status /= nf90_noerr) call handle_nf90_err(status)
        status = nf90_put_att(ncid_out,accel_divf_vert_id,'missing_value',MISSING_VALUE)
        if (status /= nf90_noerr) call handle_nf90_err(status)
        status = nf90_put_att(ncid_out,accel_divf_vert_id,'_FillValue',MISSING_VALUE)
        if (status /= nf90_noerr) call handle_nf90_err(status)

        ! Fy definition
        name = 'Fy'
        status = nf90_def_var(ncid_out,trim(name),nf90_float,dimids=(/latid,levid,timeid/),varid=fyid)
        if (status /= nf90_noerr) call handle_nf90_err(status)
        status = nf90_put_att(ncid_out,fyid,'long_name','Latitudinal component of EP flux')
        if (status /= nf90_noerr) call handle_nf90_err(status)
        status = nf90_put_att(ncid_out,fyid,'units','m^3/s^2')
        if (status /= nf90_noerr) call handle_nf90_err(status)
        status = nf90_put_att(ncid_out,fyid,'missing_value',MISSING_VALUE)
        if (status /= nf90_noerr) call handle_nf90_err(status)
        status = nf90_put_att(ncid_out,fyid,'_FillValue',MISSING_VALUE)
        if (status /= nf90_noerr) call handle_nf90_err(status)

        ! Fz definition
        name = 'Fz'
        status = nf90_def_var(ncid_out,trim(name),nf90_float,dimids=(/latid,levid,timeid/),varid=fzid)
        if (status /= nf90_noerr) call handle_nf90_err(status)
        status = nf90_put_att(ncid_out,fzid,'long_name','Vertical component of EP flux')
        if (status /= nf90_noerr) call handle_nf90_err(status)
        status = nf90_put_att(ncid_out,fzid,'units','Pa*m^2/s^2')
        if (status /= nf90_noerr) call handle_nf90_err(status)
        status = nf90_put_att(ncid_out,fzid,'missing_value',MISSING_VALUE)
        if (status /= nf90_noerr) call handle_nf90_err(status)
        status = nf90_put_att(ncid_out,fzid,'_FillValue',MISSING_VALUE)
        if (status /= nf90_noerr) call handle_nf90_err(status)

        ! vstar definition
        name = 'vstar'
        status = nf90_def_var(ncid_out,trim(name),nf90_float,dimids=(/latid,levid,timeid/),varid=vstarid)
        if (status /= nf90_noerr) call handle_nf90_err(status)
        status = nf90_put_att(ncid_out,vstarid,'long_name','Latitudinal component of residual mean meridional circulation')
        if (status /= nf90_noerr) call handle_nf90_err(status)
        status = nf90_put_att(ncid_out,vstarid,'units','m/s')
        if (status /= nf90_noerr) call handle_nf90_err(status)
        status = nf90_put_att(ncid_out,vstarid,'missing_value',MISSING_VALUE)
        if (status /= nf90_noerr) call handle_nf90_err(status)
        status = nf90_put_att(ncid_out,vstarid,'_FillValue',MISSING_VALUE)
        if (status /= nf90_noerr) call handle_nf90_err(status)

       ! wstar definition
        name = 'wstar'
        status = nf90_def_var(ncid_out,trim(name),nf90_float,dimids=(/latid,levid,timeid/),varid=wstarid)
        if (status /= nf90_noerr) call handle_nf90_err(status)
        status = nf90_put_att(ncid_out,wstarid,'long_name','Vertical component of residual mean meridional circulation')
        if (status /= nf90_noerr) call handle_nf90_err(status)
        status = nf90_put_att(ncid_out,wstarid,'units','Pa/s')
        if (status /= nf90_noerr) call handle_nf90_err(status)
        status = nf90_put_att(ncid_out,wstarid,'missing_value',MISSING_VALUE)
        if (status /= nf90_noerr) call handle_nf90_err(status)
        status = nf90_put_att(ncid_out,wstarid,'_FillValue',MISSING_VALUE)
        if (status /= nf90_noerr) call handle_nf90_err(status)



        ! end of definitions, leave define mode (NF90_ENDDEF)
        status = nf90_enddef(ncid_out)
        if (status /= nf90_noerr) call handle_nf90_err(status)

! FILL DIMENSIONS AND VARIABLES WITH DATA

        ! latitude
        status = nf90_put_var(ncid_out,latvoutid,lat,start=(/1/),count=(/ny/))
        if (status /= nf90_noerr) call handle_nf90_err(status)
        
        ! pressure level
        status = nf90_put_var(ncid_out,levvoutid,plev,start=(/1/),count=(/nl/))
        if (status /= nf90_noerr) call handle_nf90_err(status)
        
        ! time
        status = nf90_put_var(ncid_out,timevoutid,time,start=(/1/),count=(/nt/))
        if (status /= nf90_noerr) call handle_nf90_err(status)
        
        ! accel_divf_lat
        status = nf90_put_var(ncid_out,accel_divf_lat_id,accel_divf_lat(:,:,:),start=(/1,1,1/),count=(/ny,nl,nt/))
        if (status /= nf90_noerr) call handle_nf90_err(status)
        
        ! accel_divf_vert
        status = nf90_put_var(ncid_out,accel_divf_vert_id,accel_divf_vert,start=(/1,1,1/),count=(/ny,nl,nt/))
        if (status /= nf90_noerr) call handle_nf90_err(status)

        ! Fy
        status = nf90_put_var(ncid_out,fyid,fy(:,:,:),start=(/1,1,1/),count=(/ny,nl,nt/))
        if (status /= nf90_noerr) call handle_nf90_err(status)

        ! Fz
        status = nf90_put_var(ncid_out,fzid,fz(:,:,:),start=(/1,1,1/),count=(/ny,nl,nt/))
        if (status /= nf90_noerr) call handle_nf90_err(status)
        
        ! vstar
        status = nf90_put_var(ncid_out,vstarid,vstar(:,:,:),start=(/1,1,1/),count=(/ny,nl,nt/))
        if (status /= nf90_noerr) call handle_nf90_err(status)

        ! wstar
        status = nf90_put_var(ncid_out,wstarid,wstar(:,:,:),start=(/1,1,1/),count=(/ny,nl,nt/))
        if (status /= nf90_noerr) call handle_nf90_err(status)

        

! CLOSE FILES
        ! close output file
        status = nf90_close(ncid_out)
        if (status /= nf90_noerr) call handle_nf90_err(status)

        ! close input file
        status = nf90_close(ncid_in)
        if (status /= nf90_noerr) call handle_nf90_err(status)
        
! DEALLOCATE MEMORY
        deallocate(wstar)
        deallocate(vstar)
        deallocate(fz)
        deallocate(fy)
        deallocate(accel_divf_vert)
        deallocate(accel_divf_lat)
        deallocate(Wwind)
        deallocate(Vwind)
        deallocate(Uwind)
        deallocate(Temp)
        deallocate(lon)
        deallocate(lat)
        deallocate(plev)
        deallocate(time)


        ! finished
        write(*,*)'Successfully completed!'
        write(*,*)' '



end program TEM_Diagnostics

! -----------------------------------------------------------------------------
! -----------------------------------------------------------------------------
! -----------------------------------------------------------------------------
! -----------------------------------------------------------------------------
! -----------------------------------------------------------------------------
subroutine calc_TEM_Diagnostics(maxdimx,maxdimy,maxlev,maxtime,lat,plev,temp_in,u_wind_in,v_wind_in,w_wind_in &
                                       ,accel_divf_lat,accel_divf_vert,fy,fz,vstar,wstar)

USE netcdf
USE messy_main_constants_mem, ONLY: dp, er => radius_earth &
                                  , omega, pi &
                                  , gradtorad => DTR &
                                  , MISSING_VALUE => FLAGGED_BAD

implicit none
intrinsic      :: sin, cos, abs


!==================================================================================================


! input arguments of function calc_TEM_Diagnostics that must not be changed
      integer,                                        intent(in)  :: maxdimx,maxdimy,maxlev,maxtime  ! grid dimensions
      real(dp),dimension(maxdimy),                        intent(in)  :: lat                             ! latitude grid
      real(dp),dimension(maxlev),                         intent(in)  :: plev                            ! pressure grid
      real(dp),dimension(maxdimx,maxdimy,maxlev,maxtime), intent(in)  :: temp_in, u_wind_in, v_wind_in, w_wind_in
      
! output variables
      real(dp),dimension(maxdimy,maxlev,maxtime),         intent(out) :: accel_divf_lat, accel_divf_vert
      real(dp),dimension(maxdimy,maxlev,maxtime),         intent(out) :: fy, fz, vstar, wstar

! local parameters
!!$      real(dp), parameter   :: er        = 6371.0e3_dp        ! radius of the earth
!!$      real(dp), parameter   :: omega     = 7.292e-5_dp        ! angular velocity of earth
      real(dp), parameter   :: kappa     = 0.2857_dp          ! = R / c_p = 2/7
      real(dp), parameter   :: pref      = 101325.0_dp        ! reference pressure
!!$      real(dp), parameter   :: pi        = 3.141592654_dp
!!$      real(dp), parameter   :: gradtorad = pi/180.0_dp
      
            
! percentage of data points that must exist for calculation of a zonal mean
      real(dp)              :: zonal_mean_threshold = 0.9_dp
      
! loop variables
      integer                                                          :: ix, iy, ilev, itime

! variables for TEM calculations
      
      ! Coriolis parameter
      real(dp), dimension(maxdimy)                                         :: lat_in_rad, cosphi, rcosphi, f_coriolis
      
      ! denominator of du/dphi and dpsi/dphi for central difference and backward/forward difference
      real(dp)                                                             :: delta_lat_in_rad
      ! denominator of du/dp for central difference and backward/forward difference
      real(dp)                                                             :: delta_plev
      
      ! potential temperature
      real(dp), dimension(maxdimx,maxdimy,maxlev,maxtime)                  :: theta
      
      ! derivatives w.r.t. pressure and latitude
      real(dp), dimension(maxdimy,maxlev,maxtime)                          :: dtdp, dthetadp, dudp, dudphi
      
      ! variable needed for EP flux
      real(dp), dimension(maxdimy,maxlev,maxtime)                          :: psi, dpsidp, dpsidphi

      ! zonal means
      real(dp), dimension(maxdimy,maxlev,maxtime)                          :: thetazon, tzon, uzon, vzon, wzon
      
      ! zonal means of products of departures from zonal means
      real(dp), dimension(maxdimy,maxlev,maxtime)                          :: vsthetas, vsts, vsus, wsus
      



      
! for treatment of missing values
      
      integer                                            :: ntmp       ! number of longitudes with non-missing values
      
      ! contains either 1 or MISSING_VALUE when value exists or not, respectively
      real(dp),    dimension(maxdimx,maxdimy,maxlev,maxtime) :: mv         ! missing values of entire grid
      real(dp),    dimension(maxdimy,maxlev,maxtime)         :: mvzon      ! missing values for zonally averaged quantities
      real(dp),    dimension(maxdimy,maxlev,maxtime)         :: mvzondp    ! missing values for dtdp, dthetadp, dudp
      real(dp),    dimension(maxdimy,maxlev,maxtime)         :: mvzondphi  ! missing values for dudphi
      real(dp),    dimension(maxdimy,maxlev,maxtime)         :: mvpsify    ! missing values for psi, fy
      real(dp),    dimension(maxdimy,maxlev,maxtime)         :: mvfz       ! missing values for fz
      real(dp),    dimension(maxdimy,maxlev,maxtime)         :: mvdivf     ! missing values for accel_divf_lat, accel_divf_vert
      real(dp),    dimension(maxdimy,maxlev,maxtime)         :: mvdpsidp   ! missing values for dpsidp
      real(dp),    dimension(maxdimy,maxlev,maxtime)         :: mvdpsidphi ! missing values for dpsidphi

      

     
     ! Introduce a number for missing values which should be < 0.0 and ABS(.) > 1.e10
!!$     real(dp),parameter                      :: MISSING_VALUE = -1.0E+34
     
!==================================================================================================
! some checks
!do iy = 1,maxdimy
!      print*,'lat(',iy,') = ',lat(iy)
!end do

!do ilev = 1,maxlev
!      print*,'plev(',ilev,') = ',plev(ilev)
!end do

! check initial values
!do iy = 1,maxdimy
!       print*,'uzon(',iy,',1,1) = ',uzon(iy,1,1)
!end do

!write(*,*)' pi = ',pi



! Coriolis parameter
      do iy=1,maxdimy
            lat_in_rad(iy) = lat(iy)*gradtorad
            f_coriolis(iy) = 2.0*omega*sin( lat_in_rad(iy) )
            cosphi(iy)     = cos( lat_in_rad(iy) )
            rcosphi(iy)    = er*cosphi(iy)
      end do
      
! check for zeros in denominator of derivatives
      ! denominator of du/dphi
      do iy=2,maxdimy
            delta_lat_in_rad = lat_in_rad(iy-1) - lat_in_rad(iy)
            
            ! Avoid division by zero. Also, lat is decreasingly ordered
            ! ensuring that delta_lat_in_rad > 0
            if ( delta_lat_in_rad .le. 0.0 ) then
                  write(*,*)'Error in denominator of du/dphi: delta_lat_in_rad(',iy,')=',delta_lat_in_rad
            end if  
      end do
      
     ! denominator of du/dp
      do ilev=2,maxlev
            delta_plev = plev(ilev-1) - plev(ilev)
            
            ! Avoid division by zero. Also, plev is decreasingly ordered
            ! ensuring that delta_plev > 0
            if ( delta_plev .le. 0.0 ) then
                  write(*,*)'Error in denominator of du/dp: delta_plev(',ilev,')=',delta_plev
            end if  
      end do
      
      
! for the treatment of missing values
      ! the correct method would be to read the missing_value for ...
      ! ... each variable from the netCDF variable attribute ...
      ! ... but not every variable has such an attribute ...
      ! ... here the Q&D method is used ...
      do itime=1,maxtime
            do ilev=1,maxlev
                  do iy=1,maxdimy
                        do ix=1,maxdimx
                              if ( &
                                    (abs(u_wind_in(ix,iy,ilev,itime)) > 1.e10 ) .or. &
                                    (abs(v_wind_in(ix,iy,ilev,itime)) > 1.e10 ) .or. &
                                    (abs(w_wind_in(ix,iy,ilev,itime)) > 1.e10 ) .or. &
                                    (temp_in(ix,iy,ilev,itime) < 0.0) ) then
                                    mv(ix,iy,ilev,itime) = MISSING_VALUE
                              else
                                    mv(ix,iy,ilev,itime) = 1.0
                              end if
                        end do
                  end do
            end do
      end do

!==================================================================================================
! Compute TEM diagnostics for entire input data set and for every time step

      !****************************************************************
      ! 1. potential temperature
      theta(:,:,:,:) = MISSING_VALUE
      
      do itime=1,maxtime
            do ilev=1,maxlev
                  do iy=1,maxdimy
                        do ix=1,maxdimx
                              if (mv(ix,iy,ilev,itime) < 0.0) cycle
                              theta(ix,iy,ilev,itime) = temp_in(ix,iy,ilev,itime)*(pref/plev(ilev))**kappa
                              ! print*, 'theta= ', theta(ix,iy,ilev,itime)
                        end do
                  end do
            end do
      end do

      !****************************************************************
      ! 2. zonal means of u, v, w, T, Theta...
      uzon(:,:,:)     = 0.0
      vzon(:,:,:)     = 0.0
      wzon(:,:,:)     = 0.0
      tzon(:,:,:)     = 0.0
      thetazon(:,:,:) = 0.0
      
      ! ...and zonal means of meridional waves: vsus     = zon(v'u')
      !                                         wsus     = zon(w'u')
      !                                         vsts     = zon(v't')
      !                                         vsthetas = zon(v'theta')
      ! = meridional heat and momentum transport
      vsus(:,:,:)     = 0.0
      wsus(:,:,:)     = 0.0
      vsts(:,:,:)     = 0.0
      vsthetas(:,:,:) = 0.0
            
      ! initialize missing values for zonally averaged quantities
      mvzon(:,:,:)    = 1.0
      
      do itime=1,maxtime  
            do ilev=1,maxlev
                  do iy=1,maxdimy
                        !*******************
                        ! 1. compute zonal means of u, v, w, T, Theta
                        ntmp = 0    ! count the number of longitudes with non-missing values
                        do ix=1,maxdimx
                              if (mv(ix,iy,ilev,itime) < 0.0) cycle
                              
                              ! ELSE: no missing value
                              ntmp = ntmp + 1
                              
                              uzon(iy,ilev,itime)     = uzon(iy,ilev,itime)     + u_wind_in(ix,iy,ilev,itime)
                              vzon(iy,ilev,itime)     = vzon(iy,ilev,itime)     + v_wind_in(ix,iy,ilev,itime)
                              wzon(iy,ilev,itime)     = wzon(iy,ilev,itime)     + w_wind_in(ix,iy,ilev,itime)
                              tzon(iy,ilev,itime)     = tzon(iy,ilev,itime)     + temp_in(ix,iy,ilev,itime)
                              thetazon(iy,ilev,itime) = thetazon(iy,ilev,itime) + theta(ix,iy,ilev,itime)
                        end do
                  
                        if (ntmp .ge. zonal_mean_threshold*maxdimx) then
                              uzon(iy,ilev,itime)     =     uzon(iy,ilev,itime)/float(ntmp)
                              vzon(iy,ilev,itime)     =     vzon(iy,ilev,itime)/float(ntmp)
                              wzon(iy,ilev,itime)     =     wzon(iy,ilev,itime)/float(ntmp)
                              tzon(iy,ilev,itime)     =     tzon(iy,ilev,itime)/float(ntmp)
                              thetazon(iy,ilev,itime) = thetazon(iy,ilev,itime)/float(ntmp)
                        else
                              mvzon(iy,ilev,itime)    = MISSING_VALUE
                        
                              uzon(iy,ilev,itime)     = MISSING_VALUE
                              vzon(iy,ilev,itime)     = MISSING_VALUE
                              wzon(iy,ilev,itime)     = MISSING_VALUE
                              tzon(iy,ilev,itime)     = MISSING_VALUE
                              thetazon(iy,ilev,itime) = MISSING_VALUE
                        end if
                        
                        
                        ! 2. compute zonal means of meridional waves
                        if ( ntmp .lt. zonal_mean_threshold*maxdimx) then
                              vsus(iy,ilev,itime)     = MISSING_VALUE
                              wsus(iy,ilev,itime)     = MISSING_VALUE
                              vsts(iy,ilev,itime)     = MISSING_VALUE
                              vsthetas(iy,ilev,itime) = MISSING_VALUE
                        else
                              do ix=1,maxdimx
                                    if (mv(ix,iy,ilev,itime) < 0.0) cycle
                              
                                    vsus(iy,ilev,itime)     = vsus(iy,ilev,itime)+&
                                                              (v_wind_in(ix,iy,ilev,itime) - vzon(iy,ilev,itime))*&
                                                              (u_wind_in(ix,iy,ilev,itime) - uzon(iy,ilev,itime))
              
                                    wsus(iy,ilev,itime)     = wsus(iy,ilev,itime)+&
                                                              (w_wind_in(ix,iy,ilev,itime) - wzon(iy,ilev,itime))*&
                                                              (u_wind_in(ix,iy,ilev,itime) - uzon(iy,ilev,itime))
                  
                                    vsts(iy,ilev,itime)     = vsts(iy,ilev,itime)+&
                                                              (v_wind_in(ix,iy,ilev,itime) - vzon(iy,ilev,itime))*&
                                                              (  temp_in(ix,iy,ilev,itime) - tzon(iy,ilev,itime))
                  
                                    vsthetas(iy,ilev,itime) = vsthetas(iy,ilev,itime)+&
                                                              (v_wind_in(ix,iy,ilev,itime) - vzon(iy,ilev,itime))*&
                                                              (    theta(ix,iy,ilev,itime) - thetazon(iy,ilev,itime))
                              end do
                  
                              if (ntmp > 0) then                              
                                    vsus(iy,ilev,itime)     =     vsus(iy,ilev,itime)/float(ntmp)
                                    wsus(iy,ilev,itime)     =     wsus(iy,ilev,itime)/float(ntmp)
                                    vsts(iy,ilev,itime)     =     vsts(iy,ilev,itime)/float(ntmp)
                                    vsthetas(iy,ilev,itime) = vsthetas(iy,ilev,itime)/float(ntmp)
                              else
                                    write(*,*)'ERROR: ntmp <= 0'
                              end if
                        end if
                        !*******************
                  end do
            end do
      end do
            
      
      
      

      !****************************************************************
      ! 3. derivatives w.r.t. pressure and latitude: dt/dp, dtheta/dp, du/dp, du/dphi
      ! derivatives are calculated via the central difference where possible
      ! boundaries afford forward and backward differences
      !
      ! PRESSURE: The following formulas assume decreasingly ordered levels, i.e. plev(i-1) > plev(i) > plev(i+1)
      !
      !                   --> dt/dp[i] = 0.5*( (t[i-1] - t[i])/(p[i-1] - p[i]) + (t[i] - t[i+1])/(p[i] - p[i+1]) )
      !
      !           (For increasingly ordered pressure levels, plev and u,v,w,t are flipped upside down
      !            in the pressure level dimension above)
      !
      ! LATITUDE: The following formulas assume decreasingly ordered latitudes
      !           beginning at phi = 90°, ending at phi=-90°, i.e. lat(i-1) > lat(i) > lat(i+1)
      !           (The formula for du/dphi is analog to the one for dt/dp.)
      !
         
      dtdp(:,:,:)     = MISSING_VALUE
      dthetadp(:,:,:) = MISSING_VALUE
      dudp(:,:,:)     = MISSING_VALUE
      dudphi(:,:,:)   = MISSING_VALUE
      
      ! initialize missing values for derivatives of zonally averaged quantities
      mvzondp(:,:,:)    = 1.0
      mvzondphi(:,:,:)  = 1.0
         
      ! 3.1 dt/dp, dtheta/dp and du/dp
      do itime=1,maxtime
            ! first, exlude pressure boundary layers, i.e. ilev=1 and ilev=maxlev
            do ilev=2,maxlev-1
                  do iy=1,maxdimy
                        ! Missing values for zonally averaged variables
                        if( (mvzon(iy,ilev-1,itime) < 0.0) .or. (mvzon(iy,ilev,itime) < 0.0) .or. (mvzon(iy,ilev+1,itime) < 0.0) ) then
                        
                              mvzondp(iy,ilev,itime)   = MISSING_VALUE
                        else
                              dtdp(iy,ilev,itime)     = 0.5*( (tzon(iy,ilev-1,itime)-tzon(iy,ilev,itime))/ &
                                                              (   plev(ilev-1)      -   plev(ilev)      )  &
                                                              +                                            &
                                                              (tzon(iy,ilev,itime)-tzon(iy,ilev+1,itime))/ &
                                                              (   plev(ilev)      -   plev(ilev+1)      )  &
                                                            )
                                                      
                              dthetadp(iy,ilev,itime) = 0.5*( (thetazon(iy,ilev-1,itime)-thetazon(iy,ilev,itime))/ &
                                                              (       plev(ilev-1)      -       plev(ilev)      )  &
                                                              +                                                    &
                                                              (thetazon(iy,ilev,itime)-thetazon(iy,ilev+1,itime))/ &
                                                              (       plev(ilev)      -       plev(ilev+1)      )  &
                                                            )

!
!                              dthetadp(iy,ilev,itime) = ( thetazon(iy,ilev,itime) / tzon(iy,ilev,itime)              )* &
!                                                        ( dtdp(iy,ilev,itime) - tzon(iy,ilev,itime)*kappa/plev(ilev) )
!
                             
                              dudp(iy,ilev,itime)     = 0.5*( (uzon(iy,ilev-1,itime)-uzon(iy,ilev,itime))/ &
                                                              (   plev(ilev-1)      -   plev(ilev)      )  &
                                                              +                                            &
                                                              (uzon(iy,ilev,itime)-uzon(iy,ilev+1,itime))/ &
                                                              (   plev(ilev)      -   plev(ilev+1)      )  &
                                                            )
                        end if
                  end do
            end do
         
            ! second, special cases for boundaries, ilev=1...
            do iy=1,maxdimy
                  if( (mvzon(iy,1,itime) < 0.0) .or. (mvzon(iy,2,itime) < 0.0) ) then
                  
                        mvzondp(iy,1,itime)   = MISSING_VALUE
                  else
                        dtdp(iy,1,itime)     = (    tzon(iy,1,itime)-    tzon(iy,2,itime))/(plev(1)-plev(2))
                        
                        dthetadp(iy,1,itime) = (thetazon(iy,1,itime)-thetazon(iy,2,itime))/(plev(1)-plev(2))
                               
                        dudp(iy,1,itime)     = (    uzon(iy,1,itime)-    uzon(iy,2,itime))/(plev(1)-plev(2))
                  end if
            end do
            
            ! ...and ilev=maxlev
            do iy=1,maxdimy
                  if( (mvzon(iy,maxlev-1,itime) < 0.0) .or. (mvzon(iy,maxlev,itime) < 0.0) ) then
                        
                        mvzondp(iy,maxlev,itime)   = MISSING_VALUE
                  else
                        dtdp(iy,maxlev,itime)     = (tzon(iy,maxlev-1,itime)-tzon(iy,maxlev,itime))/&
                                                    (plev(maxlev-1)-plev(maxlev))
                                                  
                        dthetadp(iy,maxlev,itime) = (thetazon(iy,maxlev-1,itime)-thetazon(iy,maxlev,itime))/&
                                                    (plev(maxlev-1)-plev(maxlev))
                        
                        dudp(iy,maxlev,itime)     = (uzon(iy,maxlev-1,itime)-uzon(iy,maxlev,itime))/&
                                                    (plev(maxlev-1)-plev(maxlev))
                  end if
            end do
      end do
      !write(*,*)'dtdp = ', dtdp(32,maxlev)
         
         
      ! 3.2 du/dphi
      do itime=1,maxtime
            ! Again, first exclude the boundaries, i.e. North pole (iy=1) and South pole (iy=mxdimy)
            do ilev=1,maxlev
                  do iy=2,maxdimy-1
                        if( (mvzon(iy-1,ilev,itime) < 0.0) .or. (mvzon(iy,ilev,itime) < 0.0) .or. (mvzon(iy+1,ilev,itime) < 0.0) ) then
                        
                              mvzondphi(iy,ilev,itime) = MISSING_VALUE
                        else

                              dudphi(iy,ilev,itime)= ( 1.0/rcosphi(iy) )*&
                                                     0.5*( ( uzon(iy-1,ilev,itime)*cosphi(iy-1) -   &
                                                             uzon(iy,ilev,itime)*cosphi(iy)      )/ &
                                                           ( lat_in_rad(iy-1) - lat_in_rad(iy)   )  &
                                                           +                                        &
                                                           ( uzon(iy,ilev,itime)*cosphi(iy) -       &
                                                             uzon(iy+1,ilev,itime)*cosphi(iy+1)  )/ &
                                                           ( lat_in_rad(iy) - lat_in_rad(iy+1)   )  &
                                                         )
                                                              
                        end if
                  end do
            end do
            
            ! Then treat North pole...
            do ilev=1,maxlev
                  if( (mvzon(1,ilev,itime) < 0.0) .or. (mvzon(2,ilev,itime) < 0.0) ) then
                  
                        mvzondphi(1,ilev,itime) = MISSING_VALUE
                  else
                        dudphi(1,ilev,itime) = ( 1.0/rcosphi(1) )*&
                                               ( uzon(1,ilev,itime)*cosphi(1) - uzon(2,ilev,itime)*cosphi(2) )/&
                                               (                lat_in_rad(1) -                lat_in_rad(2) )
                  end if
            end do
                  
            ! ...and South pole
            do ilev=1,maxlev
                  if( (mvzon(maxdimy-1,ilev,itime) < 0.0) .or. (mvzon(maxdimy,ilev,itime) < 0.0) ) then 
                  
                        mvzondphi(maxdimy,ilev,itime) = MISSING_VALUE
                  else
                        dudphi(maxdimy,ilev,itime) = ( 1.0/rcosphi(maxdimy) )*&
                                                     ( uzon(maxdimy-1,ilev,itime)*cosphi(maxdimy-1) -&
                                                       uzon(maxdimy,ilev,itime)*cosphi(maxdimy)        )/&
                                                     ( lat_in_rad(maxdimy-1) - lat_in_rad(maxdimy)     )
                  end if
            end do
      end do
      
      !****************************************************************
      ! 4. psi and latitudinal and pressure component of EP flux
      
      psi(:,:,:) = MISSING_VALUE
      fy(:,:,:)  = MISSING_VALUE
      fz(:,:,:)  = MISSING_VALUE
      
      ! initialize missing values for psi, fy and fz
      mvpsify(:,:,:) = 1.0
      mvfz(:,:,:)    = 1.0
      
      do itime=1,maxtime
            do ilev=1,maxlev
                  do iy=1,maxdimy
                        ! at first treat Fy and psi which only need zonal means and derivatives w.r.t. pressure
                        if( (mvzon(iy,ilev,itime) < 0.0) .or. (mvzondp(iy,ilev,itime) < 0.0) ) then
                              
                              mvpsify(iy,ilev,itime) = MISSING_VALUE
                        else
                              ! Two possibilities to compute psi of EP flux
                              ! psi(iy,ilev,itime) = -vsts(iy,ilev,itime)/( kappa*tzon(iy,ilev,itime)/plev(ilev) - dtdp(iy,ilev,itime) )
                              psi(iy,ilev,itime) = vsthetas(iy,ilev,itime)/dthetadp(iy,ilev,itime)
                  
                              fy(iy,ilev,itime) = rcosphi(iy)*&
                                                  ( -vsus(iy,ilev,itime) + psi(iy,ilev,itime)*dudp(iy,ilev,itime) )
                        end if
                        
                        ! second, treat Fz which needs zonal means and both derivatives w.r.t. pressure and latitude
                        if( (mvzon(iy,ilev,itime) < 0.0) .or. (mvzondp(iy,ilev,itime) < 0.0) .or. (mvzondphi(iy,ilev,itime) < 0.0) ) then
                              
                              mvfz(iy,ilev,itime) = MISSING_VALUE
                        else                                            
                              fz(iy,ilev,itime) = rcosphi(iy)*&
                                                  ( -wsus(iy,ilev,itime) - psi(iy,ilev,itime)*(dudphi(iy,ilev,itime) - f_coriolis(iy)) )
                        end if
                  end do
            end do
      end do


      !****************************************************************
      ! 5. EP flux diveregence
      ! The EP flux divergence at (iy,ilev,itime) only exists if
      !          fy exists at (iy-1, ilev ,itime), (iy,ilev,itime) and (iy+1, ilev ,itime)
      !      and fz exists at ( iy ,ilev-1,itime), (iy,ilev,itime) and ( iy ,ilev+1,itime)
      !
      ! Although the latitudinal contribution of the EP flux diveregence can exist independently
      ! of the vertical contribution of the EP flux diveregence and the other way round,
      ! both latitudinal and vertical contributions are decided to be missing, if one of them is missing!
      !
      ! At the same time missing values for dpsi/dp and dpsi/dphi can be determined
      
      ! initialize missing values for accel_divf_lat and accel_divf_vert
      mvdivf(:,:,:) = 1.0
      ! initialize missing values for dpsi/dp and dpsi/dphi
      mvdpsidp(:,:,:)   = 1.0
      mvdpsidphi(:,:,:) = 1.0
      
      ! determine missing values of accel_divf_lat and accel_divf_vert
      do itime=1,maxtime
            ! 1. check for missing values of dfy/dphi in EP flux divergence
            !    which also causes a missing value of dpsi/dphi
            do ilev=1,maxlev
                  ! 1.1 exclude poles
                  do iy=2,maxdimy-1
                        if( ( mvpsify(iy-1,ilev,itime) < 0.0 ) .or. ( mvpsify(iy,ilev,itime) < 0.0 ) .or. ( mvpsify(iy+1,ilev,itime) < 0.0 ) ) then
                        
                              mvdivf(iy,ilev,itime)     = MISSING_VALUE
                              mvdpsidphi(iy,ilev,itime) = MISSING_VALUE
                        end if
                  end do
                  
                  ! 1.2 North pole
                  if( ( mvpsify(1,ilev,itime) < 0.0 ) .or. ( mvpsify(2,ilev,itime) < 0.0 ) ) then
                  
                        mvdivf(1,ilev,itime)     = MISSING_VALUE
                        mvdpsidphi(1,ilev,itime) = MISSING_VALUE
                  end if
                  
                  ! 1.3 South pole
                  if( ( mvpsify(maxdimy-1,ilev,itime) < 0.0 ) .or. ( mvpsify(maxdimy,ilev,itime) < 0.0 ) ) then
                  
                        mvdivf(maxdimy,ilev,itime)     = MISSING_VALUE
                        mvdpsidphi(maxdimy,ilev,itime) = MISSING_VALUE
                  end if
            end do
            
            ! 2. check for missing values of dfz/dp in EP flux divergence
            do iy=1,maxdimy
                  ! 2.1 exclude bottom and top of atmosphere
                  do ilev=2,maxlev-1
                        if( ( mvfz(iy,ilev-1,itime) < 0.0 ) .or. ( mvfz(iy,ilev,itime) < 0.0 ) .or. ( mvfz(iy,ilev+1,itime) < 0.0 ) ) then
                        
                              mvdivf(iy,ilev,itime)   = MISSING_VALUE
                        end if
                  end do
                  
                  ! 2.2 bottom
                  if( ( mvfz(iy,1,itime) < 0.0 ) .or. ( mvfz(iy,2,itime) < 0.0 ) ) then
                        
                        mvdivf(iy,1,itime)   = MISSING_VALUE
                  end if
                  
                  ! 2.3 top
                  if( ( mvfz(iy,maxlev-1,itime) < 0.0 ) .or. ( mvfz(iy,maxlev,itime) < 0.0 ) ) then
                        
                        mvdivf(iy,maxlev,itime)   = MISSING_VALUE
                  end if
            end do
            
            ! 3. check for missing values of dpsi/dp
            do iy=1,maxdimy
                  ! 3.1 exclude bottom and top of atmosphere
                  do ilev=2,maxlev-1
                        if( ( mvpsify(iy,ilev-1,itime) < 0.0 ) .or. ( mvpsify(iy,ilev,itime) < 0.0 ) .or. ( mvpsify(iy,ilev+1,itime) < 0.0 ) ) then
                        
                              mvdpsidp(iy,ilev,itime) = MISSING_VALUE
                        end if
                  end do
                  
                  ! 3.2 bottom
                  if( ( mvpsify(iy,1,itime) < 0.0 ) .or. ( mvpsify(iy,2,itime) < 0.0 ) ) then
                        
                        mvdpsidp(iy,1,itime) = MISSING_VALUE
                  end if
                  
                  ! 3.3 top
                  if( ( mvpsify(iy,maxlev-1,itime) < 0.0 ) .or. ( mvpsify(iy,maxlev,itime) < 0.0 ) ) then
                        
                        mvdpsidp(iy,maxlev,itime) = MISSING_VALUE
                  end if
            end do
      end do
      
      
      ! initialize EP flux diveregence
      accel_divf_lat(:,:,:)  = MISSING_VALUE
      accel_divf_vert(:,:,:) = MISSING_VALUE
      
      ! determine latitudinal contribution to EP flux divergence
      do itime=1,maxtime
            do ilev=1,maxlev
                  ! Exclude poles
                  do iy=2,maxdimy-1
                        if( mvdivf(iy,ilev,itime) < 0.0 ) cycle
            
                        accel_divf_lat(iy,ilev,itime) = ( 1.0/rcosphi(iy) )*&
                                                    0.5 * ( ( fy(iy-1,ilev,itime)*cosphi(iy-1) -    &
                                                              fy(iy,ilev,itime)*cosphi(iy)       )/ &
                                                            ( lat_in_rad(iy-1) - lat_in_rad(iy)  )  &
                                                            +                                       &
                                                            ( fy(iy,ilev,itime)*cosphi(iy) -        &
                                                              fy(iy+1,ilev,itime)*cosphi(iy+1)   )/ &
                                                            ( lat_in_rad(iy) - lat_in_rad(iy+1)  )  &
                                                          )
                  end do
            
                  ! North pole
                  if( .not. mvdivf(1,ilev,itime) < 0.0 ) then
            
                        accel_divf_lat(1,ilev,itime) = ( 1.0/rcosphi(1) )*&
                                                   ( fy(1,ilev,itime)*cosphi(1) -&
                                                     fy(2,ilev,itime)*cosphi(2)    )/&
                                                   ( lat_in_rad(1) - lat_in_rad(2) )
                  end if

                  ! South pole
                  if( .not. mvdivf(maxdimy,ilev,itime) < 0.0 ) then
            
                        accel_divf_lat(maxdimy,ilev,itime) = ( 1.0/rcosphi(maxdimy) )*&
                                                         ( fy(maxdimy-1,ilev,itime)*cosphi(maxdimy-1) -&
                                                             fy(maxdimy,ilev,itime)*cosphi(maxdimy)      )/&
                                                         (  lat_in_rad(maxdimy-1) - lat_in_rad(maxdimy)  )
                  end if
            end do
      end do
      
      
      
      ! determine vertical contribution to EP flux divergence
      do itime=1,maxtime
            do iy=1,maxdimy
                  ! Exclude bottom and top of atmosphere
                  do ilev=2,maxlev-1
                        if( mvdivf(iy,ilev,itime) < 0.0 ) cycle
                        
                        accel_divf_vert(iy,ilev,itime) = 0.5 * ( ( fz(iy,ilev-1,itime) - fz(iy,ilev,itime) )/ &
                                                                 (  plev(ilev-1)       -  plev(ilev)       )  &
                                                                +                                            &
                                                                 ( fz(iy,ilev,itime) - fz(iy,ilev+1,itime) )/ &
                                                                 (  plev(ilev)       -  plev(ilev+1)       )  &
                                                               )
                  end do
      
                  ! bottom
                  if( .not. mvdivf(iy,1,itime) < 0.0 ) then
                  
                        accel_divf_vert(iy,1,itime) = ( fz(iy,1,itime) - fz(iy,2,itime) )/&
                                                      (  plev(1)       -  plev(2)       )
                  end if

                  ! top
                  if( .not. mvdivf(iy,maxlev,itime) < 0.0 ) then
                  
                        accel_divf_vert(iy,maxlev,itime) = ( fz(iy,maxlev-1,itime) - fz(iy,maxlev,itime) )/&
                                                           (  plev(maxlev-1)       -  plev(maxlev)       )
                  end if
            end do
      end do
      
      ! check for unexpected missing values
      ntmp = 0
      
      ! divide latitudinal and vertical contributions to EP flux divergence
      ! by rcosphi to obtain the acceleration exerted by the EP flux divergence
      ! and finally add the two contributions
      do itime=1,maxtime
            do iy=1,maxdimy
                  do ilev=1,maxlev
                        if( mvdivf(iy,ilev,itime) < 0.0 ) cycle
                        
                        accel_divf_lat(iy,ilev,itime)  = accel_divf_lat(iy,ilev,itime)  / rcosphi(iy)
                        accel_divf_vert(iy,ilev,itime) = accel_divf_vert(iy,ilev,itime) / rcosphi(iy)

                        if( (abs(accel_divf_lat(iy,ilev,itime) ) > 1.0e20) .or. &
                            (abs(accel_divf_vert(iy,ilev,itime)) > 1.0e20) )      then
                              ntmp = ntmp + 1
                              write(*,*)'unexpected missing value of accel_divf at:'
                              write(*,*)' accel_divf_lat (iy=',iy,',ilev=',ilev,',itime=',itime,')=',accel_divf_lat(iy,ilev,itime)
                              write(*,*)' accel_divf_vert(iy=',iy,',ilev=',ilev,',itime=',itime,')=',accel_divf_vert(iy,ilev,itime)
                              write(*,*)' mvdivf =',mvdivf(iy,ilev,itime)
                        end if
                  end do
            end do
      end do
      
      
      write(*,*)'Number of unexpected missing values in accel_divf: ',ntmp
      
      
      
      
      
      !****************************************************************
      ! 6. dpsi/dp and dpsi/dphi
      
      ! dpsi/dphi
      dpsidphi(:,:,:) = MISSING_VALUE
      
      do itime=1,maxtime
            do ilev=1,maxlev
            
                  ! Exclude poles
                  do iy=2,maxdimy-1
                        if( mvdpsidphi(iy,ilev,itime) < 0.0 ) cycle
            
                        dpsidphi(iy,ilev,itime) = ( 1.0/rcosphi(iy) )*                             &
                                                  0.5*(   ( psi(iy-1,ilev,itime)*cosphi(iy-1) -    &
                                                            psi(iy,ilev,itime)*cosphi(iy)    )/    &
                                                            ( lat_in_rad(iy-1) - lat_in_rad(iy))   &
                                                         +( psi(iy,ilev,itime)*cosphi(iy) -        &
                                                            psi(iy+1,ilev,itime)*cosphi(iy+1)    )/&
                                                            ( lat_in_rad(iy) - lat_in_rad(iy+1))   &
                                                       )
                  end do
            
                  ! North pole
                  if( .not. mvdpsidphi(1,ilev,itime) < 0.0 ) then
            
                        dpsidphi(1,ilev,itime) = ( 1.0/rcosphi(1) )*&
                                                   ( psi(1,ilev,itime)*cosphi(1) -&
                                                     psi(2,ilev,itime)*cosphi(2)    )/&
                                                   ( lat_in_rad(1) - lat_in_rad(2)  )
                  end if

                  ! South pole
                  if( .not. mvdpsidphi(maxdimy,ilev,itime) < 0.0 ) then
            
                        dpsidphi(maxdimy,ilev,itime) = ( 1.0/rcosphi(maxdimy) )*&
                                                         ( psi(maxdimy-1,ilev,itime)*cosphi(maxdimy-1) -&
                                                             psi(maxdimy,ilev,itime)*cosphi(maxdimy)      )/&
                                                         ( lat_in_rad(maxdimy-1) - lat_in_rad(maxdimy)    )
                  end if
            end do
      end do
      
      
      ! dpsi/dp
      dpsidp(:,:,:)   = MISSING_VALUE
      
      do itime=1,maxtime
            do iy=1,maxdimy
            
                  ! Exclude bottom and top of atmosphere
                  do ilev=2,maxlev-1
                        if( mvdpsidp(iy,ilev,itime) < 0.0 ) cycle
                        
                        dpsidp(iy,ilev,itime)   = 0.5*(  ( psi(iy,ilev-1,itime) - psi(iy,ilev,itime) )/ &
                                                         (   plev(ilev-1)       -   plev(ilev)       )  &
                                                         +                                              &
                                                         ( psi(iy,ilev,itime) - psi(iy,ilev+1,itime) )/ &
                                                         (   plev(ilev)       -   plev(ilev+1)       )  &
                                                      )
                  end do
      
                  ! bottom
                  if( .not. mvdpsidp(iy,1,itime) < 0.0 ) then
                  
                        dpsidp(iy,1,itime)      = ( psi(iy,1,itime) - psi(iy,2,itime) )/&
                                                          ( plev(1) - plev(2) )
                  end if

                  ! top
                  if( .not. mvdpsidp(iy,maxlev,itime) < 0.0 ) then
                  
                        dpsidp(iy,maxlev,itime) = ( psi(iy,maxlev-1,itime) - psi(iy,maxlev,itime) )/&
                                                          ( plev(maxlev-1) - plev(maxlev) )
                  end if
            end do
      end do
      
      
      !****************************************************************
      ! 7. vstar and wstar
      
      vstar(:,:,:) = MISSING_VALUE
      wstar(:,:,:) = MISSING_VALUE
      
      do itime=1,maxtime
            do ilev=1,maxlev
                  do iy=1,maxdimy
                  
                        ! vstar
                        if( (.not. mvzon(iy,ilev,itime) < 0.0) .and. (.not. mvdpsidp(iy,ilev,itime) < 0.0) ) then
                        
                              vstar(iy,ilev,itime) = vzon(iy,ilev,itime) - dpsidp(iy,ilev,itime)
                        
                        end if
                        
                        ! wstar
                        if( (.not. mvzon(iy,ilev,itime) < 0.0) .and. (.not. mvdpsidphi(iy,ilev,itime) < 0.0) ) then
                        
                              wstar(iy,ilev,itime) = wzon(iy,ilev,itime) + dpsidphi(iy,ilev,itime)
                        
                        end if
                  
                  end do
            end do
      end do
      
      

      ! check for unexpected missing values in vstar
      ntmp = 0
      
      do itime=1,maxtime
            do iy=1,maxdimy
                  do ilev=1,maxlev
                        ! vstar
                        if( (.not. mvzon(iy,ilev,itime) < 0.0) .and. (.not. mvdpsidp(iy,ilev,itime) < 0.0) ) then
                              
                              if( abs(vstar(iy,ilev,itime)) > 1.0e10 ) then
                                    ntmp = ntmp + 1
                                    write(*,*)'unexpected missing value of vstar at:'
                                    write(*,*)' vstar(iy=',iy,',ilev=',ilev,',itime=',itime,')=',vstar(iy,ilev,itime)
                              end if
                        end if
                  end do
            end do
      end do
      
      write(*,*)'Number of unexpected missing values in vstar: ',ntmp
      
      
      ! check for unexpected missing values in wstar
      ntmp = 0
      
      do itime=1,maxtime
            do iy=1,maxdimy
                  do ilev=1,maxlev
                        ! wstar
                        if( (.not. mvzon(iy,ilev,itime) < 0.0) .and. (.not. mvdpsidphi(iy,ilev,itime) < 0.0) ) then
                              
                              if( abs(wstar(iy,ilev,itime)) > 1.0e10 ) then
                                    ntmp = ntmp + 1
                                    write(*,*)'unexpected missing value of wstar at:'
                                    write(*,*)' wstar(iy=',iy,',ilev=',ilev,',itime=',itime,')=',wstar(iy,ilev,itime)
                              end if
                        end if
                  end do
            end do
      end do
      
      write(*,*)'Number of unexpected missing values in wstar: ',ntmp

   
end subroutine calc_TEM_Diagnostics




!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
! Deal with a netcdf error message
subroutine handle_nf90_err(status)

        use netcdf
  
        implicit none

        integer,intent(in)      ::      status

        if(status /= nf90_noerr)then
          write(*,*)'netCDF status error: ',trim(nf90_strerror(status))
          stop "Execution Stopped"
        end if

end subroutine handle_nf90_err
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------

!******************************************************************************
!
!
!   function is_netcdf_file(dateiname) result(erg)
!   subroutine nc_check_error (error,str,abort)
!   subroutine nc_get_vertcoorname (filename|ncid, vertcoorname, pref)
!   subroutine nc_grid_descr(file,nlon,nlat,longrid,latgrid,rcode)
!   subroutine nc_get_level  (infile, nlevs, levels, dir, err)
!   subroutine nc_get_values3 (status,filename,varname,values,dir)   
!   subroutine nc_get_var_cf (status, ncid, varid, varname, values)
!   subroutine nc_get_values3_cf (status,filename,varname,values,dir)                             
!   subroutine nc_get_param_2dint (status, ncid, varname, iarr)
!   
!******************************************************************************

MODULE messy_clams_tools_ncutils

  INTERFACE nc_get_vertcoorname
     MODULE PROCEDURE nc_get_vertcoorname_ncid
     MODULE PROCEDURE nc_get_vertcoorname_filename  
  END INTERFACE

CONTAINS

  !*************************************************************************
  !  Function to verify if a file is a netcdf-file
  !*************************************************************************
  function is_netcdf_file(dateiname) result(erg)

    implicit none

    character(len=*)     :: dateiname
    logical              :: erg
    character(len=3)     :: txt

    open(unit=13, file=dateiname, status='old', position='rewind', action='read')
    read(unit=13,fmt=*) txt
    if ( txt .eq. 'CDF' ) then
       erg = .true.
    else
       erg = .false.
    endif
    close(13)
    
  end function is_netcdf_file

  !*************************************************************************
  ! handels NetCDF errors 
  !*************************************************************************
  Subroutine nc_check_error (error,str,abort)

    use netcdf

    implicit none
    
    integer             :: error
    character(*)        :: str
    logical, optional   :: abort
    
    logical :: stop_prg

    if (error/=nf90_noerr) then

       stop_prg = .false.
       if (present(abort)) stop_prg = abort

       write (*,*)
       write (*,*) '***********************************************************'
       write (*,*) 'NetCDF Error:'
       write (*,*) nf90_strerror(error)
       write (*,*) 
       write (*,*) trim(str)
       write (*,*) '***********************************************************'
       write (*,*)
       if (stop_prg) then
          write (*,*) 'THE PROGRAM WILL STOP !!!'   
          write (*,*)
          write (*,*)
          stop
       endif
    endif
    
  end Subroutine nc_check_error

  !*************************************************************                  
  ! get name of vertical coordinate
  !*************************************************************                  
  subroutine nc_get_vertcoorname_filename (filename, vertcoorname, pref)

    use netcdf
    use messy_clams_global, only: PREC
    use messy_clams_tools_utils, only:lowercase

    implicit none

    character(*)         :: filename, vertcoorname
    real(PREC), optional :: pref

!    character(64) :: valid_char='abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789-_'
    character(40) :: helpstr
    integer :: ncid, status, i, k, dummy_id

    helpstr = ""
    vertcoorname = ""

    status = nf90_open(filename,NF90_NOWRITE,ncid)
    if (status /= nf90_noerr) then
       write (*,*) 'Can not open file ',TRIM(filename),' !'
!!!!! besser: error code setzen ?
       stop
    else
       status = nf90_get_att (ncid,NF90_GLOBAL,'exp_VERTCOOR_name',helpstr)
       if (status == nf90_noerr) then
          helpstr = lowercase(ADJUSTL(helpstr))

          ! delete control characters (\0 etc.)
          k = 1
          do i = 1, LEN_TRIM(helpstr)
             if ((IACHAR(helpstr(i:i))>31) .AND. (IACHAR(helpstr(i:i))<127)) then
                vertcoorname(k:k) = helpstr(i:i)
                k = k+1
!               else
!                  print *,'invalid char in vertcoorname deleted !'
             endif
          enddo

!           k = 1
!           do i = 1, LEN_TRIM(helpstr)
!              if (INDEX(valid_char,helpstr(i:i))/=0) then
!                 vertcoorname(k:k) = helpstr(i:i)
!                 k = k+1
!              else
!                 print *,'invalid char in vertcoorname!'
!              endif
!           enddo

       else

          status= nf90_inq_dimid (ncid,'press',dummy_id)     
          if (status==nf90_noerr) then
             vertcoorname = 'press'
          else
             status= nf90_inq_dimid (ncid,'hybrid',dummy_id)     
             if (status==nf90_noerr) then
                vertcoorname = 'hybrid'
             else
                status= nf90_inq_dimid (ncid,'zeta',dummy_id)     
                if (status==nf90_noerr) then
                   vertcoorname = 'zeta'
                else
                   vertcoorname = 'theta'
                endif
             endif
          endif

       endif
    endif

    if (present(pref)) then
       status = nf90_get_att (ncid,NF90_GLOBAL,'exp_VERTCOOR_ref_level',pref)
       if (status/=nf90_noerr) pref=-1
    endif

    status = nf90_close(ncid)

  end subroutine nc_get_vertcoorname_filename

  !*************************************************************                  
  ! get name of vertical coordinate
  !*************************************************************                  
  subroutine nc_get_vertcoorname_ncid (ncid, vertcoorname, pref)

    use netcdf
    use messy_clams_global, only: PREC
    use messy_clams_tools_utils, only: lowercase

    implicit none

    integer              :: ncid
    character(*)         :: vertcoorname
    real(PREC), optional :: pref

    character(40) :: helpstr
    integer      ::  status, i, k, dummy_id

    vertcoorname = ""
    helpstr = ""

    status = nf90_get_att (ncid,NF90_GLOBAL,'exp_VERTCOOR_name',helpstr)
    if (status == nf90_noerr) then

       ! delete control characters (\0 etc.)
       helpstr = lowercase(ADJUSTL(helpstr))

       k = 1
       do i = 1, LEN_TRIM(helpstr)
          if ((IACHAR(helpstr(i:i))>31) .AND. (IACHAR(helpstr(i:i))<127)) then
             vertcoorname(k:k) = helpstr(i:i)
             k = k+1
!           else
!              print *,'invalid char in vertcoorname deleted !'
          endif
       enddo
       
    else

       status= nf90_inq_dimid (ncid,'press',dummy_id)     
       if (status==nf90_noerr) then
          vertcoorname = 'press'
       else
          status= nf90_inq_dimid (ncid,'hybrid',dummy_id)     
          if (status==nf90_noerr) then
             vertcoorname = 'hybrid'
          else
             status= nf90_inq_dimid (ncid,'zeta',dummy_id)     
             if (status==nf90_noerr) then
                vertcoorname = 'zeta'
             else
                vertcoorname = 'theta'
             endif
          endif
       endif
        
    endif

    if (present(pref)) then
       status = nf90_get_att (ncid,NF90_GLOBAL,'exp_VERTCOOR_ref_level',pref)
       if (status/=nf90_noerr) pref=-1
    endif

  end subroutine nc_get_vertcoorname_ncid

  !******************************************************************             
  !
  !******************************************************************             
  subroutine nc_get_var_atts (rcode, filename, varname, longname, units)

    use netcdf

    implicit none

    character(*), intent(in)              :: filename, varname
    character(*), intent(inout), optional :: longname, units
    integer,      intent(out)             :: rcode

    integer      :: ncid, varid
    
    rcode = nf90_open (filename,nf90_nowrite,ncid)
    call nc_check_error (rcode,'Error on open file '//trim(filename),abort=.false.)
    if (rcode/=0) return
   
    rcode = nf90_inq_varid (ncid,varname,varid)
    call nc_check_error (rcode,'Cannot find variable '//trim(varname),abort=.false.)
    if (rcode/=0) return

    if (present(longname)) then
       rcode = nf90_get_att (ncid,varid,"long_name",longname)
       call nc_check_error (rcode,'Cannot find attribute long_name',abort=.false.)
       if (rcode/=0) longname='???'
    endif
    if (present(units)) then
       rcode = nf90_get_att (ncid,varid,"units",units)
       call nc_check_error (rcode,'Cannot find attribute units',abort=.false.)
       if (rcode/=0) units='???'
    endif

!!!!! weitere Attribute !?!

    rcode = nf90_close (ncid)

  end subroutine nc_get_var_atts

  !******************************************************************             
  !
  !******************************************************************             
  SUBROUTINE nc_grid_descr(file,nlon,nlat,longrid,latgrid,rcode)

    use netcdf
    use messy_clams_global, only: PREC, buffersize

    IMPLICIT NONE

    !Input parameters
    character*(*) :: file

    !Output parameters
    REAL(PREC), dimension(:), pointer :: longrid, latgrid 
    INTEGER     :: nlon, nlat, rcode

    ! local variables
    integer :: nin, dimid, varid

    rcode = nf90_open(file,nf90_nowrite,nin, buffersize)
    call nc_check_error (rcode,'Error on open file '//trim(file),abort=.false.)
    if (rcode/=0) return

    rcode = nf90_inq_dimid (nin,"lon",dimid)
    call nc_check_error (rcode,'Cannot find dimension lon',abort=.false.)
    if (rcode/=0) return
    rcode = nf90_inquire_dimension (nin,dimid,len=nlon)
    call nc_check_error (rcode,'Cannot read dimension lon',abort=.false.)
    if (rcode/=0) return

    rcode = nf90_inq_dimid (nin,"lat",dimid)
    call nc_check_error (rcode,'Cannot find dimension lat',abort=.false.)
    if (rcode/=0) return
    rcode = nf90_inquire_dimension (nin,dimid,len=nlat)
    call nc_check_error (rcode,'Cannot read dimension lat',abort=.false.)
    if (rcode/=0) return

    allocate (longrid(nlon))
    allocate (latgrid(nlat))

    rcode = nf90_inq_varid (nin,'lon',varid)
    call nc_check_error (rcode,'Cannot find variable lon',abort=.false.)
    if (rcode/=0) return
    rcode = nf90_get_var (nin,varid,longrid)
    call nc_check_error (rcode,'Cannot read variable lon',abort=.false.)
    if (rcode/=0) return

    rcode = nf90_inq_varid (nin,'lat',varid)
    call nc_check_error (rcode,'Cannot find variable lat',abort=.false.)
    if (rcode/=0) return
    rcode = nf90_get_var (nin,varid,latgrid)
    call nc_check_error (rcode,'Cannot read variable lat',abort=.false.)
    if (rcode/=0) return

    rcode = nf90_close(nin)
    call nc_check_error (rcode,'Error on closing file '//trim(file),abort=.false.)
    if (rcode/=0) return

  END SUBROUTINE nc_grid_descr



  !*************************************************************                  
  !
  !*************************************************************                  
  SUBROUTINE nc_get_level  (infile, nlevs, levels, dir, err)                               

    use netcdf
    use messy_clams_global, only: PREC, buffersize, filenamelen
    use messy_clams_tools_utils, only: lowercase

    IMPLICIT NONE                                                                 

    ! Inputparameter                                                              
    CHARACTER*(*) :: infile                                                       

    ! Outputparameter                                                              
    INTEGER                           :: nlevs                                                
    REAL(PREC), dimension(:), pointer :: levels

    !Optional Parameters                                                          
    INTEGER,INTENT(OUT),OPTIONAL :: err                                           
    CHARACTER*(*),INTENT(IN),OPTIONAL :: dir                                      

    ! local variables                                                             
    INTEGER                :: ncid, rcode, varid, dimid
    character(80)          :: vertcoorname=""
    CHARACTER(filenamelen) :: inputfile                                                     

    ! Check Optional Parameters

    IF (PRESENT(err)) err=0                                                       

    inputfile=infile                                                              
    IF (PRESENT(dir)) THEN                                                        
       IF (dir /= '') THEN                                                         
          inputfile=TRIM(ADJUSTL(dir))//'/'//TRIM(ADJUSTL(infile))                  
       ELSE                                                                        
          inputfile=TRIM(ADJUSTL(infile))                                           
       ENDIF
    ENDIF

    ! open ncdf-inputfile 
    rcode=nf90_open(inputfile,NF90_NOWRITE,ncid, buffersize)

    if (rcode /= 0) then
       PRINT *,'************************************************************'
       PRINT *,'  ERROR !!!! '                                               
       PRINT *,'  CANNOT OPEN NCDF- INPUT-FILE'
       PRINT *,'************************************************************'
       IF (PRESENT(err)) err=-1                                                    
       RETURN                                                                      
    endif

    ! get name of vertical coordinate
    call nc_get_vertcoorname (ncid,vertcoorname)

    ! get level- dimension- id 
    rcode = nf90_inq_dimid(ncid,TRIM(vertcoorname),dimid)
    if (rcode /= nf90_noerr) then
       PRINT *,'************************************************************'
       PRINT *,'  ERROR !!!! '                                               
       PRINT *,'  dimension of vertical coordinate not found ! '
       PRINT *
       PRINT *,'  vertcoorname=',trim(vertcoorname)
       PRINT *
       PRINT *,'  in sub. nc_get_level nach Aufruf von nc_get_vertcoorname'
       PRINT *,'  (geaendert am 27.01.06)'
       PRINT *,'************************************************************'
       IF (PRESENT(err)) err=-2                                                    
       RETURN                                                                      
    endif

    !get level- number
    rcode=nf90_inquire_dimension(ncid,dimid,len=nlevs)


    !get id of theta-lev or press-lev 
    rcode = nf90_inq_varid(ncid,TRIM(vertcoorname),varid)
    if (rcode /= nf90_noerr) then
       PRINT *,'************************************************************'
       PRINT *,'  ERROR !!!! '                                               
       PRINT *,'  vertical coordinate not found ! '
       PRINT *,'************************************************************'       
       IF (PRESENT(err)) err=-4                                                    
       RETURN
    endif

    allocate (levels(nlevs))

    !get theta- or press- values
    rcode = nf90_get_var(ncid,varid,levels)
    if (rcode /= nf90_noerr) then
       if (present(err)) err=-5
    endif

    rcode=nf90_close(ncid)

  end subroutine nc_get_level

   !***********************************************************************       
   ! get all values for a NetCDF variable and return these on a
   ! 3-dimensional field
   ! 
   !***********************************************************************        
   SUBROUTINE nc_get_values3 (status,filename,varname,values,dir)  
   
     use netcdf
     USE messy_clams_global,          ONLY: prec, buffersize, filenamelen
     USE messy_clams_tools_packutils, ONLY: get_packed_var4

     implicit none

     integer                           :: status          
     real(prec)                        :: values(:,:,:)                                                
     character(*),intent(in)           :: filename,varname                                              
     character(*),intent(in),optional  :: dir                                      

     ! local variables
     real(prec), dimension(:,:,:,:), allocatable :: dumarr4
     character(filenamelen)   :: ncdf_file   
     integer                  :: ncid, varid

     status = 0

     ncdf_file=TRIM(ADJUSTL(filename))
     IF (PRESENT(dir)) ncdf_file=TRIM(ADJUSTL(dir))//'/'//TRIM(ADJUSTL(filename))    
     
     ! open netCDF file
     status = nf90_open (ncdf_file,nf90_nowrite,ncid, buffersize)
     call nc_check_error (status,'Error on open netCDF input file',abort=.false.)
     if (status /= nf90_noerr) return
    
     ! get variable ID
     status = nf90_inq_varid (ncid,varname,varid)
     call nc_check_error (status,'Error on reading variable-id',abort=.false.)
     if (status /= nf90_noerr) return

     ! variables in netcdf files are four-dimensional (nlev,nlat,nlon,ntime)
     allocate(dumarr4(size(values,1),size(values,2),size(values,3),1))
     
     ! read variable (and unpack if it is packed)
     status = get_packed_var4 (ncid, varid, dumarr4)
     if (status /= nf90_noerr) then
        write(*,*) 'error on reading variable ',trim(varname)
        write(*,*) nf90_strerror(status)
        return
     endif

     values = RESHAPE(dumarr4,(/size(values,1),size(values,2),size(values,3)/))

     deallocate (dumarr4)
     
     status=nf90_close(ncid)
     call nc_check_error (status,'Error on closing input file',abort=.false.)

   end subroutine nc_get_values3

   
  !***********************************************************************       
  ! get all values for a NetCDF variable and return these on a
  ! 3-dimensional field
  ! If the dataset is not following CF conventions, switch dimensions 
  !  -> (lon,lat,lev)
  !***********************************************************************       
  subroutine nc_get_var_cf (status, ncid, varid, varname, values)

    USE netcdf
    USE messy_clams_tools_packutils, ONLY: get_packed_var4
    USE messy_clams_global,          ONLY: prec

    implicit none
    
    integer                 :: status, ncid, varid
    character(*)            :: varname
    real(prec)              :: values(:,:,:)
    
    real(prec),dimension(:,:,:,:),allocatable :: dumarr4
    character(40)                             :: conventions=''
    integer                                   :: ix, iy, iz
    
    status = 0

    ! get global attribute "Conventions"
    status = nf90_get_att (ncid, NF90_GLOBAL, 'Conventions', conventions)
    if (status /= 0)   conventions = 'noname'

    ! if datasest is following CF conventions: read variable
    if (conventions(1:2) == 'CF') then

       status = nf90_get_var(ncid,varid,values)
       call nc_check_error (status, &
            'Variable '//trim(varname)//' could not be read !!!', abort=.false.)
       if (status /= nf90_noerr) return

    ! if dataset is not following CF conventions: 
    ! read variable, unpack variable and switch dimensions
    else

       allocate(dumarr4(size(values,3),size(values,2), size(values,1),1))

       status = get_packed_var4 (ncid, varid, dumarr4)
       call nc_check_error (status,'Error in get_packed_var4: ' &
            //trim(varname)//' could not be read !!!',abort=.false.)
       if (status /= nf90_noerr) return
              
       do iz = 1, size(values,3)
          do iy = 1, size(values,2)
             do ix= 1, size(values,1)
                values(ix,iy,iz) = dumarr4(iz,iy,ix,1)
             enddo
          enddo
       enddo
      
       deallocate(dumarr4)

    endif

  end subroutine nc_get_var_cf

 
  !***********************************************************************       
  ! get all values for a NetCDF variable and return these on a
  ! 3-dimensional field
  ! If the dataset is not following CF conventions, switch dimensions 
  !  -> (lon,lat,lev)
  !***********************************************************************        
  SUBROUTINE nc_get_values3_cf (status,filename,varname,values,dir)                             

     USE netcdf
     USE messy_clams_global,          ONLY: prec, buffersize, filenamelen

     integer                           :: status          
     real(prec)                        :: values(:,:,:)                                     
     character(*),intent(in)           :: filename,varname                                  
     character(*),intent(in),optional  :: dir                                      

     ! local variables
     character(filenamelen)   :: ncdf_file   
     integer                  :: ncid, varid, rcode

     status = 0

     ncdf_file=TRIM(ADJUSTL(filename))
     IF (PRESENT(dir)) ncdf_file=TRIM(ADJUSTL(dir))//'/'//TRIM(ADJUSTL(filename))    
     
     ! open netCDF file
     status = nf90_open (ncdf_file,nf90_nowrite,ncid, buffersize)
     call nc_check_error (status,'error on open netCDF input file',abort=.false.)
     if (status /= nf90_noerr) return

     ! get variable ID
     status = nf90_inq_varid (ncid,varname,varid)
     call nc_check_error (status,'error on reading variable-id',abort=.false.)
     if (status /= nf90_noerr) return

     ! read variable to values(lon,lat,lev)
     call nc_get_var_cf (status, ncid, varid, varname, values)
     if (status /= nf90_noerr) return

     ! close netCDF file
     rcode=nf90_close(ncid)

  end subroutine nc_get_values3_cf
  
  !*******************************************************************
  !
  !*******************************************************************
  subroutine nc_get_param_2dint (status, ncid, varname, iarr)

    use netcdf
    
    implicit none

    integer                 :: status, ncid
    character(*)            :: varname
    integer, dimension(:,:) :: iarr
    
    integer         :: varid, ix, iy
    character(40)   :: conventions=''
    integer, dimension(:,:), allocatable :: ihelparr

    status = 0
    
    ! get global attribute "Conventions"
    status = nf90_get_att (ncid, NF90_GLOBAL, 'Conventions', conventions)
    if (status /= 0)   conventions = 'noname'
    
    status = nf90_inq_varid (ncid,varname,varid)
    call nc_check_error (status,"Cannot find "//trim(varname),abort=.false.)
    if (status /= nf90_noerr) return
    
    ! if datasest is following CF conventions: read variable
    if (conventions(1:2) == 'CF') then
       
       status = nf90_get_var(ncid,varid,iarr)
       call nc_check_error (status,"Cannot read "//trim(varname),abort=.false.)
       if (status /= nf90_noerr) return
       
    else
       
       allocate (ihelparr(size(iarr,2),size(iarr,1)))
       
       status = nf90_get_var(ncid,varid,ihelparr)
       call nc_check_error (status,"Cannot read "//trim(varname),abort=.false.)
       if (status /= nf90_noerr) return
       
       do iy = 1, size(iarr,2)
          do ix = 1, size(iarr,1)
             iarr(ix,iy) = ihelparr(iy,ix)
          enddo
       enddo
      
       deallocate (ihelparr)
       
    endif

  end subroutine nc_get_param_2dint

END MODULE messy_clams_tools_ncutils

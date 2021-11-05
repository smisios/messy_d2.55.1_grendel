!!! this file includes subroutines to be used to read NetCDF files, in particular EMAC output
!
! getDimLen(infile,nlon,nlat,nlev,nt)
! returns the dimensionlength in infile if time, level, lat, lon
!
! readData(infile,varname,data,nlon,nlat,nlev,nt)
! reads the variable 'varname' from infile into data, where data is of size(nlon,nlat,nlev,nt)
!
! subroutine readData3D(infile,varname,data,nx,ny,nz)
! as above but for 3-dim data
!
! subroutine readData1D(infile,varname,data,nlen)
! as above but fpr 1-dimensional data of length nlen
!
! Hella Garny, 14.3.2011

!------------------------------------------------------------------------------------------------
subroutine getDimLen(infile,nlon,nlat,nlev,nt)

use netcdf

implicit none

character(len=80), intent(in)  :: infile
character(len=nf90_max_name) :: dname

integer :: ncid_in, status, timeid, latid, lonid, levid, varid, nt, nlat, nlon, nlev

! read in lat, lon, levels, time, variables

write(*,*) ' '
write(*,*) 'getting dimension lengths'

! open input netcdf file
status = nf90_open(trim(infile),nf90_nowrite,ncid_in)
if(status /= nf90_noerr) call handle_nf_err(status)

! find the number of times in the input file
status = nf90_inq_dimid(ncid_in,'time',timeid)
if (status /= nf90_noerr) call handle_nf_err(status)
status = nf90_inquire_dimension(ncid_in,timeid,dname,nt)
if (status /= nf90_noerr) call handle_nf_err(status)

! find the number of lats in the input file
status = nf90_inq_dimid(ncid_in,'lat',latid)
if (status /= nf90_noerr) call handle_nf_err(status)
status = nf90_inquire_dimension(ncid_in,latid,dname,nlat)
if (status /= nf90_noerr) call handle_nf_err(status)

! find the number of lons in the input file
status = nf90_inq_dimid(ncid_in,'lon',lonid)
if (status /= nf90_noerr) call handle_nf_err(status)
status = nf90_inquire_dimension(ncid_in,lonid,dname,nlon)
if (status /= nf90_noerr) call handle_nf_err(status)

! find the number of levs in the input file
! op_pj_20110505+
!!$status = nf90_inq_dimid(ncid_in,'mlev',levid)
status = nf90_inq_dimid(ncid_in,'lev',levid)
! op_pj_20110505-
if (status /= nf90_noerr) call handle_nf_err(status)
status = nf90_inquire_dimension(ncid_in,levid,dname,nlev)
if (status /= nf90_noerr) call handle_nf_err(status)

status = nf90_close(ncid_in)
if (status /= nf90_noerr) call handle_nf_err(status)

end subroutine getDimLen

!--------------------------------------------------------------------

subroutine readData(infile,varname,data,nlon,nlat,nlev,nt)

use netcdf

implicit none

character(len=80), intent(in)  :: infile
character(len=*), intent(in) :: varname

integer :: ncid_in, status,  varid

integer, intent (in) ::  nt, nlat, nlon, nlev


real,intent(inout) :: data(nlon,nlat,nlev,nt)

! read in lat, lon, levels, time, variables

write(*,*) ' '
write(*,*) 'reading input file'
write(*,*) varname

! open input netcdf file
status = nf90_open(trim(infile),nf90_nowrite,ncid_in)
if(status /= nf90_noerr) call handle_nf_err(status)

status = nf90_inq_varid(ncid_in,varname,varid)
if (status /= nf90_noerr) call handle_nf_err(status)
status = nf90_get_var(ncid_in,varid,data,(/1,1,1,1/),(/nlon,nlat,nlev,nt/))
if (status /= nf90_noerr) call handle_nf_err(status)

print*, 'data read'

status = nf90_close(ncid_in)
if (status /= nf90_noerr) call handle_nf_err(status)


end subroutine readData

!--------------------------------------------------------------------

subroutine readData3D(infile,varname,data,nx,ny,nz)

use netcdf

implicit none

character(len=80), intent(in)  :: infile
character(len=*), intent(in) :: varname

integer :: ncid_in, status,  varid

integer, intent (in) ::  nx, ny, nz


real,intent(inout) :: data(nx, ny, nz)

! read in lat, lon, levels, time, variables

write(*,*) ' '
write(*,*) 'reading input file'
write(*,*) varname

! open input netcdf file
status = nf90_open(trim(infile),nf90_nowrite,ncid_in)
if(status /= nf90_noerr) call handle_nf_err(status)

status = nf90_inq_varid(ncid_in,varname,varid)
if (status /= nf90_noerr) call handle_nf_err(status)
status = nf90_get_var(ncid_in,varid,data,(/1,1,1/),(/nx, ny, nz/))
if (status /= nf90_noerr) call handle_nf_err(status)

print*, 'data read'

status = nf90_close(ncid_in)
if (status /= nf90_noerr) call handle_nf_err(status)


end subroutine readData3D

!--------------------------------------------------------------------

subroutine readData1D(infile,varname,data,nlen)

use netcdf

implicit none

character(len=80), intent(in)  :: infile
character(len=*), intent(in) :: varname

integer :: ncid_in, status,  varid

integer, intent (in) ::  nlen


real,intent(inout) :: data(nlen)

! read in lat, lon, levels, time, variables

write(*,*) ' '
write(*,*) 'reading input file'
write(*,*) varname

! open input netcdf file
status = nf90_open(trim(infile),nf90_nowrite,ncid_in)
if(status /= nf90_noerr) call handle_nf_err(status)

status = nf90_inq_varid(ncid_in,varname,varid)
if (status /= nf90_noerr) call handle_nf_err(status)
status = nf90_get_var(ncid_in,varid,data,(/1/),(/nlen/))
if (status /= nf90_noerr) call handle_nf_err(status)

print*, 'data read'

status = nf90_close(ncid_in)
if (status /= nf90_noerr) call handle_nf_err(status)


end subroutine readData1D


! -----------------------------------------------------------------------------
! Deal with a netCDF error message
        subroutine handle_nf_err(status)
	use netcdf
        implicit none

        integer,intent(in)      ::      status

        if(status /= nf90_noerr)then
          write(*,*)'netCDF status error: ',trim(nf90_strerror(status))
          stop "Execution Stopped"
        endif

        end subroutine handle_nf_err
! -----------------------------------------------------------------------------

! op_pj_20160511+
! -----------------------------------------------------------------------------
subroutine check_variable(infile,varname,status)

  use netcdf
  implicit none

  character(len=80), intent(in)  :: infile
  character(len=*), intent(in) :: varname
  integer, intent(out) :: status

  integer :: ncid_in, varid

  ! open input netcdf file
  status = nf90_open(trim(infile),nf90_nowrite,ncid_in)
  if(status /= nf90_noerr) call handle_nf_err(status)

  status = nf90_inq_varid(ncid_in,varname,varid)
  if (status /= nf90_noerr) then 
     status = 1
  else
     status = 0
  endif

end subroutine check_variable
! -----------------------------------------------------------------------------
! op_pj_20160511-

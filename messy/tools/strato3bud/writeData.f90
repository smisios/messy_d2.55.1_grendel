!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! programs to create and write netcdf files
!
!
! subroutine createFile4D(file,nlon,nlat,nlev,nt)
! create netcdf file with name 'file' with 4 dimensions: 'lon', 'lat', 'lev', 'time'
! with length nlon, nlat, nlev, nt
!
! subroutine addVar1D(file,VarName,DimName, data, nx)
! add a new variable with name 'VarName' to 'file'. The variable is 1-dimensional alond dimension 'DimName'
! put the data into the file: data of length nx. nx must be equal to the dimension length of DimName.
!
!


! ---------------------------------------------------------------------------
subroutine createFile4D(file,nlon,nlat,nlev,nt)

use netcdf

implicit none

character(len=80), intent(in) :: file
integer, intent(in) :: nt,nlev,nlon,nlat
integer :: i, ncid, status, tID, levID, lonID, latID

write(*,*) 
write(*,*) 'creating new file:'
write(*,*) trim(file) 

status = nf90_create(trim(file),nf90_clobber,ncid)
if (status /= nf90_NoErr) call handle_err(status)

status = nf90_def_dim(ncid,'time',nt,tID)
if (status /= nf90_NoErr) call handle_err(status)
status = nf90_def_dim(ncid,'lev',nlev,levID)
if (status /= nf90_NoErr) call handle_err(status)
status = nf90_def_dim(ncid,'lon',nlon,lonID)
if (status /= nf90_NoErr) call handle_err(status)
status = nf90_def_dim(ncid,'lat',nlat,latID)
if (status /= nf90_NoErr) call handle_err(status)

status = nf90_enddef(ncid)
if (status /= nf90_NoErr) call handle_err(status)

status = nf90_close(ncid)
if (status /= nf90_NoErr) call handle_err(status)

end subroutine createFile4D

! ---------------------------------------------------------------------------
subroutine createFile2D(file,nlat,nlev)

use netcdf

implicit none

character(len=80), intent(in) :: file
integer, intent(in) :: nlev,nlat
integer :: i, ncid, status,  levID, latID

write(*,*)
write(*,*) 'creating new file:'
write(*,*) trim(file)

status = nf90_create(trim(file),nf90_clobber,ncid)
if (status /= nf90_NoErr) call handle_err(status)

status = nf90_def_dim(ncid,'lev',nlev,levID)
if (status /= nf90_NoErr) call handle_err(status)
status = nf90_def_dim(ncid,'lat',nlat,latID)
if (status /= nf90_NoErr) call handle_err(status)

status = nf90_enddef(ncid)
if (status /= nf90_NoErr) call handle_err(status)

status = nf90_close(ncid)
if (status /= nf90_NoErr) call handle_err(status)

end subroutine createFile2D

! ----------------------------------------------------------------------------

subroutine addVar1D(file,VarName,DimName, data, nx)

use netcdf

implicit none

character(len=80), intent(in) :: file
character(len=*), intent(in) :: VarName, DimName

integer, intent(in) :: nx
real, dimension(nx), intent(in) :: data

integer :: i, ncid, status, dimID, varID

write(*,*)
write(*,*)  'writing data for var:'
write(*,*)  VarName

! open file
status = nf90_open(trim(file),nf90_write,ncid)
if(status /= nf90_NoErr) call handle_err(status)

! get ID of dimesion
status = nf90_inq_dimid(ncid,DimName,dimID)
if(status /= nf90_NoErr) call handle_err(status)

! redefine: add variable definition
status = nf90_redef(ncid)
if (status /= nf90_NoErr) call handle_err(status)

status = nf90_def_var(ncid,VarName,nf90_float,(/ dimID /),varID)
if (status /= nf90_NoErr) call handle_err(status)

status = nf90_enddef(ncid)
if (status /= nf90_NoErr) call handle_err(status)

! put data into variable
status = nf90_put_var(ncid,varID,data)
if (status /= nf90_NoErr) call handle_err(status)

status = nf90_close(ncid)
if (status /= nf90_NoErr) call handle_err(status)

end subroutine addVar1D

! ----------------------------------------------------------------------------

subroutine addVar2D(file,VarName,data,nlat,nlev)

use netcdf

implicit none

character(len=80), intent(in) :: file
character(len=*), intent(in) :: VarName

integer, intent(in) :: nlat,nlev

real, intent(in) :: data(nlat,nlev)

integer :: i, ncid, status, varID, latid, levid

write(*,*)
write(*,*)  'writing data for var:'
write(*,*)  VarName


! open file
status = nf90_open(trim(file),nf90_write,ncid)
if(status /= nf90_NoErr) call handle_err(status)

! get ID of dimesions
status = nf90_inq_dimid(ncid,'lat',latid)
if(status /= nf90_NoErr) call handle_err(status)

status = nf90_inq_dimid(ncid,'lev',levid)
if(status /= nf90_NoErr) call handle_err(status)

! redefine: add variable definition
status = nf90_redef(ncid)
if (status /= nf90_NoErr) call handle_err(status)

status = nf90_def_var(ncid,VarName,nf90_float,(/ latid, levid /),varID)
if (status /= nf90_NoErr) call handle_err(status)

status = nf90_enddef(ncid)
if (status /= nf90_NoErr) call handle_err(status)

! put data into variable

status = nf90_put_var(ncid,varID,data)
if (status /= nf90_NoErr) call handle_err(status)

status = nf90_close(ncid)
if (status /= nf90_NoErr) call handle_err(status)

end subroutine addVar2D


! ----------------------------------------------------------------------------

subroutine addVar4D(file,VarName,data,nlon,nlat,nlev,nt )

use netcdf

implicit none

character(len=80), intent(in) :: file
character(len=*), intent(in) :: VarName 
character(len=10) :: dname

integer, intent(in) :: nlon,nlat,nlev,nt

real, intent(in) :: data(nlon,nlat,nlev,nt)

integer :: i, ncid, status, varID, timeid, latid, lonid, levid

write(*,*)
write(*,*)  'writing data for var:'
write(*,*)  VarName


! open file
status = nf90_open(trim(file),nf90_write,ncid)
if(status /= nf90_NoErr) call handle_err(status)

! get ID of dimesions
status = nf90_inq_dimid(ncid,'time',timeid)
if(status /= nf90_NoErr) call handle_err(status)

status = nf90_inq_dimid(ncid,'lat',latid)
if(status /= nf90_NoErr) call handle_err(status)

status = nf90_inq_dimid(ncid,'lon',lonid)
if(status /= nf90_NoErr) call handle_err(status)

status = nf90_inq_dimid(ncid,'lev',levid)
if(status /= nf90_NoErr) call handle_err(status)

! redefine: add variable definition
status = nf90_redef(ncid)
if (status /= nf90_NoErr) call handle_err(status)

status = nf90_def_var(ncid,VarName,nf90_float,(/ lonid, latid, levid, timeid /),varID)
if (status /= nf90_NoErr) call handle_err(status)

status = nf90_enddef(ncid)
if (status /= nf90_NoErr) call handle_err(status)


! put data into variable
status = nf90_inquire_variable(ncid,varID,dname)
if (status /= nf90_NoErr) call handle_err(status)


status = nf90_put_var(ncid,varID,data)
if (status /= nf90_NoErr) call handle_err(status)

status = nf90_close(ncid)
if (status /= nf90_NoErr) call handle_err(status)

end subroutine addVar4D


! -----------------------------------------------------------------------------

subroutine handle_err(status)

use netcdf

implicit none

integer status

print *, "NetCDF error: ", trim(NF90_STRERROR(status))
print *, "Abort."
stop

end ! handle_err


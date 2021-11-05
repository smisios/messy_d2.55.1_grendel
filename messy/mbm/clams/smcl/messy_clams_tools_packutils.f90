!******************************************************************************
! File utils/src/nc_pack_utils.f90
!
! This file is part of the Chemical Lagrangian Model of the
! Stratosphere (CLaMS). CLaMS is a hierarchy of numerical models and
! preprocessors which simulate transport, chemical reactions and
! mixing processes in the stratosphere.  
!
! Copyright (c) 2006
! N. Thomas, Forschungszentrum Juelich GmbH
! Last Modified By: j.-u.grooss
! Last Modified On: Mon Apr  8 12:58:29 2019
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
!******************************************************************************
!
! Module nc_pack_utils
!
! This Module contains subroutines and functions for writing packed netcdf data.
!
!  function write_packed_attributes (ncid, varid, scale, offset, &
!                                    packfv, fv, packmdi, mdi) result(status)
!
!  function delete_packed_attributes (ncid, varid) result(status)
!
!  subroutine get_packing_accuracy (help_rvals, notmdi, relerr, scale, offset)
!
!  function get_packed_var4 (ncid, varid, values, start, count) result(status)
!     
!  function put_packed_var4 (values,nc_out,varid,dimlen) result(status)
!
!  function nc_pack_variable (nc_in,nc_out,varid_in,xtype,ndims,dimids,maxerr) &
!                             result(status) 
!
!  function nc_unpack_variable (nc_in,nc_out,varid_in,varid_out, &
!                                xtype,ndims, dimids) result(status) 
!
! At the moment only arrays of type REAL will be packed/unpacked !!!
! Arrays of type INTEGER or DOUBLE PRECISION will only be copied !!!
!
!***********************************************************************
Module messy_clams_tools_packutils


  ! short_int already defined in module "packing" 
  !integer, parameter :: short_int = SELECTED_INT_KIND(4)

contains

  !***************************************************************************
  ! write variable attributes used for packing
  !***************************************************************************
  function write_packed_attributes (ncid, varid, scale, offset, &
                                    packfv, fv, packmdi, mdi) result(status)

    use netcdf
    use messy_clams_tools_packing   
    use messy_clams_global, only: prec, sp, dp
    
    implicit none

    integer                      :: ncid, varid
!!!!!
    real                         :: fv
    real(prec)                   :: scale, offset
    integer(short_int)           :: packfv
    integer(short_int), optional :: packmdi
    real(prec),         optional :: mdi
    integer                      :: status

    real(prec)               :: valid_min, valid_max
!    integer(short_int) :: valid_min_pk, valid_max_pk

    !status = nf90_redef(ncid)
    !if (status /= nf90_noerr) return

    ! write global attribute "definition_name"
    !status = nf90_put_att (ncid, NF90_GLOBAL, "definition_name",&
    !     "def_icg_unpack")

    ! write attributes "scale_factor" and "add_offset"
    status = nf90_put_att (ncid, varid, "scale_factor", scale)
    if (status /= nf90_noerr) return
    status = nf90_put_att (ncid, varid, "add_offset", offset)
    if (status /= nf90_noerr) return

    ! pack values for attributes "valid_min" and "valid_max"
    status = nf90_get_att (ncid, varid, "valid_min", valid_min)
    if (status == nf90_noerr) then
! If valid_min is smaller than the minimum of data, it can not be converted correctly !!!
!       valid_min_pk = nint((valid_min - offset) / scale)
!       status = nf90_put_att (ncid, varid, "valid_min", valid_min_pk)
!       if (status /= nf90_noerr) return
       status = nf90_del_att (ncid, varid, "valid_min")
       if (status /= nf90_noerr) return
       status = nf90_put_att (ncid, varid, "UNPACK_valid_min", valid_min)
       if (status /= nf90_noerr) return
    endif
    status = nf90_get_att (ncid, varid, "valid_max", valid_max)
    if (status == nf90_noerr) then
! If valid_max is greater than the maximum of data, it can not be converted correctly !!!
!       valid_max_pk = nint((valid_max - offset) / scale)
!       status = nf90_put_att (ncid, varid, "valid_max", valid_max_pk)
!       if (status /= nf90_noerr) return
       status = nf90_del_att (ncid, varid, "valid_max")
       if (status /= nf90_noerr) return
       status = nf90_put_att (ncid, varid, "UNPACK_valid_max", valid_max)
       if (status /= nf90_noerr) return
    endif

    ! write missing values for packed and unpacked data
    if (present(packmdi)) then
       status = nf90_put_att (ncid, varid, "missing_value", packmdi)
       if (status /= nf90_noerr) return
    endif
    if (present(mdi)) then
       status = nf90_put_att (ncid, varid, "UNPACK_missing_value", mdi)
       if (status /= nf90_noerr) return
    endif

    ! write fill values for packed and unpacked data
    if (fv /= NF90_FILL_FLOAT) then
       status = nf90_put_att (ncid, varid, "_FillValue", packfv)
       if (status /= nf90_noerr) return
       status = nf90_put_att (ncid, varid, "UNPACK_FillValue", fv)
       if (status /= nf90_noerr) return
    endif
    
    ! set attribute "PACKED_STATUS"
    status = nf90_put_att (ncid, varid, "PACKED_STATUS", "PACKED")
    if (status /= nf90_noerr) return

    !status = nf90_enddef(ncid)
    
  end function write_packed_attributes

  !***************************************************************************
  ! deletes attributes "scale_factor", "add_offset", "UNPACK_missing_value"
  ! and "UNPACK_FillValue"
  !***************************************************************************
  function delete_packed_attributes (ncid, varid) result(status)

    use messy_clams_global, only: prec
    use netcdf
    use messy_clams_tools_packing   

    implicit none

    integer        :: ncid, varid
    integer        :: status 

    real(prec)    :: valid_min, valid_max
    real(prec)    :: mdi
!!!!!
    real          :: fillvalue    

    !status = nf90_redef(ncid)
    !if (status /= nf90_noerr) return

    ! delete global attribute "definition_name"
    status = nf90_del_att (ncid, NF90_GLOBAL, "definition_name")

    ! delete attributes "scale_factor" and "add_offset"
    status = nf90_del_att (ncid, varid, "scale_factor")
    !if (status /= nf90_noerr) return
    status = nf90_del_att (ncid, varid, "add_offset")
    !if (status /= nf90_noerr) return

    ! set original "valid_min" and "valid_max"
    status = nf90_get_att (ncid, varid, "UNPACK_valid_min", valid_min)
    if (status == nf90_noerr) then
       status = nf90_put_att (ncid, varid, "valid_min", valid_min)
       if (status /= nf90_noerr) return
       status = nf90_del_att (ncid, varid, "UNPACK_valid_min")
       if (status /= nf90_noerr) return
    endif
    status = nf90_get_att (ncid, varid, "UNPACK_valid_max", valid_max)
    if (status == nf90_noerr) then
       status = nf90_put_att (ncid, varid, "valid_max", valid_max)
       if (status /= nf90_noerr) return
       status = nf90_del_att (ncid, varid, "UNPACK_valid_max")
       if (status /= nf90_noerr) return
    endif
    
    ! write missing value and fill value
    status = nf90_get_att (ncid, varid, "UNPACK_missing_value", mdi)
    if (status == nf90_noerr) then
       status = nf90_put_att (ncid, varid, "missing_value", mdi)
       if (status /= nf90_noerr) return
       status = nf90_del_att (ncid, varid, "UNPACK_missing_value")
       if (status /= nf90_noerr) return
    endif
    status = nf90_get_att (ncid, varid, "UNPACK_FillValue", fillvalue)
    if (status == nf90_noerr) then
       status = nf90_put_att (ncid, varid, "_FillValue", fillvalue)
       if (status /= nf90_noerr) return
       status = nf90_del_att (ncid, varid, "UNPACK_FillValue")
       if (status /= nf90_noerr) return
    endif

    ! set "PACKED_STATUS"
    status = nf90_put_att (ncid, varid, "PACKED_STATUS", "UNPACKED")
    if (status /= nf90_noerr) return

    !status = nf90_enddef(ncid)

  end function delete_packed_attributes

  !***************************************************************************
  ! Estimate the accuracy of packing ?!?
  !
  !
  !
  !***************************************************************************
  subroutine get_packing_accuracy (help_rvals, notmdi, relerr, scale, offset)

    use messy_clams_global, only: prec
    use netcdf
    use messy_clams_tools_packing   

    implicit none

    real(prec), dimension(:), pointer :: help_rvals
    logical,    dimension(:), pointer :: notmdi
    real(prec), intent(out)           :: relerr
    real(prec), intent(in)            :: scale, offset 

    real(prec)                        :: min_x, max_x, mean
   
    min_x  = minval(help_rvals,mask=notmdi)
    max_x  = maxval(help_rvals,mask=notmdi)
    
    ! relerr ist nur dann der groesstmoegliche relative Fehler,
    ! wenn alle Werte positiv oder alle Werte negativ sind.
    ! Ansonsten ist der rel. Fehler in der Naehe von 0 groesser!
    if (abs(min_x)>epsilon(min_x) .and. abs(max_x)>epsilon(max_x)) then
       relerr = max(scale/abs(max_x),scale/abs(min_x))
    elseif (abs(max_x)>epsilon(max_x)) then
       relerr = scale/abs(max_x)
    elseif (abs(min_x)>epsilon(min_x)) then
       relerr = scale/abs(min_x)
    else
       relerr = 0.
    endif
    
!     !write (*,*) 'eps=',epsilon(min_x)
!     write (*,*) 'scale, offset= ', scale, offset
!     ! wg. Rundung (mit nint) kann der abs. Fehler max. halb so groß sein
!     ! wie der Skalierungsfaktor
!     write (*,*) 'max. absoluter Fehler:', scale*0.5
    
!     write (*,*) 'minimum, maximum= ', min_x, max_x
!     write (*,*) 'max. rel. Fehler bei min/max= ',relerr*0.5
    
    if (count(notmdi) == 0) then
       relerr = 0.
    else
       mean = sum(help_rvals,mask=notmdi)/count(notmdi)
       !relerr = scale / mean *0.5
!        write (*,*) 'Mittelwert:',mean
!        write (*,*) 'max. rel. Fehler fuer Mittelwert:',scale / mean * 0.5
       !write (*,*) 'test',sum(scale/help_rvals, mask=notmdi)/count(notmdi)
    endif
!    write (*,*) 

  end subroutine get_packing_accuracy

  !**********************************************************************
  ! A four-dimensional variable is read and - if the attribut 
  ! "PACKED_STATUS" is "PACKED" - it is unpacked.
  !**********************************************************************
  Function get_packed_var4 (ncid, varid, values, start, count) result(status)

    use messy_clams_global, only: prec
    use netcdf
    use messy_clams_tools_packing   

    implicit none

    integer                        :: ncid, varid, status
    real(prec), dimension(:,:,:,:) :: values
    integer, optional              :: start(4), count(4)

    integer(short_int) :: mdi_packed, fillvalue_packed
    real(prec)         :: scale, offset, mdi
!!!!!
    real               :: fillvalue
    integer            :: vartype, error
    logical            :: mdi_exist, fv_exist

    integer(short_int), dimension(:,:,:,:), pointer :: packvals 
    real(prec)        , dimension(:),   pointer     :: help_rvals
    integer(short_int), dimension(:),   pointer     :: help_packvals
    character(nf90_max_name)                        :: pkstatus=''

    ! If the attribut "PACKED_STATUS" does not exist,  
    ! the variable is read (without unpacking)
    status = nf90_get_att(ncid,varid,"PACKED_STATUS",pkstatus)
    if (status /= NF90_NOERR) then
       if (present(start)) then
          if (present(count)) then
             status = nf90_get_var(ncid,varid,values,start,count)
             if (status /= nf90_noerr) print*,nf90_strerror(status)
          else
             status = nf90_get_var(ncid,varid,values,start)
             if (status /= nf90_noerr) print*,nf90_strerror(status)
          endif
       else
          if (present(count)) then
             status = nf90_get_var(ncid,varid,values,count)
             if (status /= nf90_noerr) print*,nf90_strerror(status)
          else
             status = nf90_get_var(ncid,varid,values)
             if (status /= nf90_noerr) print*,nf90_strerror(status)
          endif
       endif

    else ! attribut "PACKED_STATUS" exists

       ! If attribut "PACKED_STATUS" is not "PACKED"
       ! the variable is read (without unpacking)
       if (pkstatus /= "PACKED") then
          if (present(start)) then
             if (present(count)) then
                status = nf90_get_var(ncid,varid,values,start,count)
                if (status /= nf90_noerr) print*,nf90_strerror(status)
             else
                status = nf90_get_var(ncid,varid,values,start)
                if (status /= nf90_noerr) print*,nf90_strerror(status)
             endif
          else
             if (present(count)) then
                status = nf90_get_var(ncid,varid,values,count)
                if (status /= nf90_noerr) print*,nf90_strerror(status)
             else
                status = nf90_get_var(ncid,varid,values)
                if (status /= nf90_noerr) print*,nf90_strerror(status)
             endif
          endif

       else ! attribut "PACKED_STATUS" is "PACKED"

          status = nf90_inquire_variable (ncid, varid, xtype=vartype)
          if (status /= nf90_noerr) return
          
          ! If variable is not of type SHORT, it cannot be unpacked !
          if (vartype/=NF90_SHORT) then
             status = -99
             write (*,*) 'Die Variable ist nicht vom Typ SHORT !!!'
             return

          else  ! read and unpack variable 

             ! read scale factor and offset
             status = nf90_get_att (ncid,varid,"scale_factor", scale)
             if (status /= nf90_noerr) return
             status = nf90_get_att (ncid,varid,"add_offset", offset)
             if (status /= nf90_noerr) return

             ! read "missing value" and/or "fill value"
             mdi_exist = .false.
             fv_exist = .false.
             status = nf90_get_att (ncid,varid,"missing_value",mdi_packed)
             if (status == nf90_noerr) then
                status = nf90_get_att (ncid,varid,"UNPACK_missing_value",mdi)
                if (status == nf90_noerr) then
                   mdi_exist = .true.
                else
                   status = -99
                   write (*,*) 'Das Attribut missing_value ist vorhanden '
                   write (*,*) 'aber UNPACK_missing_value nicht !!!'
                   return
                endif
             endif
             status = nf90_get_att (ncid,varid,"_FillValue",fillvalue_packed)
             if (status == nf90_noerr) then
                status = nf90_get_att (ncid,varid,"UNPACK_FillValue",fillvalue)
                if (status == nf90_noerr) then
                   fv_exist = .true.
                else
                   status = -99
                   write (*,*) 'Das Attribut _FillValue ist vorhanden '
                   write (*,*) 'aber UNPACK_FillValue nicht !!!'
                   return
                endif
             endif
             if (.not. mdi_exist .and. .not. fv_exist) then
                status = -99
                write (*,*) 'Es muss mindestens eines der Attribute '
                write (*,*) '"missing_value" oder "_FillValue" vorhanden sein !!!'
                return
             endif
             
             !write (*,*)  'Die Variable wird entpackt'
             !write (*,*)  'scale, offset: ', scale, offset

             allocate(packvals(size(values,1),size(values,2),&
                  & size(values,3),size(values,4)))
             allocate(help_packvals(size(values,1)* &
                  size(values,2)*size(values,3)*size(values,4)))
             allocate(help_rvals(size(values,1)*  &
                  size(values,2)*size(values,3)*size(values,4)))
             
             ! read packed values
             if (present(start)) then
                if (present(count)) then
                   status = nf90_get_var(ncid,varid,packvals,start,count)
                   if (status /= nf90_noerr) return
                else
                   status = nf90_get_var(ncid,varid,packvals,start)
                   if (status /= nf90_noerr) return
                endif
             else
                 if (present(count)) then
                   status = nf90_get_var(ncid,varid,packvals,count)
                   if (status /= nf90_noerr) return
                else
                   status = nf90_get_var(ncid,varid,packvals)
                   if (status /= nf90_noerr) return
                endif
             endif
             
             ! convert to vector
             help_packvals=reshape(packvals, &
                  (/size(values,1)*size(values,2)*size(values,3)*size(values,4)/))
             
             ! unpack values
             if (mdi_exist .and. fv_exist) then
                call unpack_array (help_rvals, help_packvals, scale, offset, error, &
                     mdi=mdi, mdi_packed=mdi_packed, &
                     fillvalue=fillvalue, fillvalue_packed=fillvalue_packed)
             elseif (mdi_exist) then
                call unpack_array (help_rvals, help_packvals, scale, offset, error, &
                     mdi=mdi, mdi_packed=mdi_packed)
             else
                call unpack_array (help_rvals, help_packvals, scale, offset, error, &
                     fillvalue=fillvalue, fillvalue_packed=fillvalue_packed)
             endif
             
             if (error/=0) then
                status = -99
                return
             endif
             
             ! convert to four-dimensional array
             values = reshape(help_rvals, &
                  (/size(values,1),size(values,2),size(values,3),size(values,4)/))
             
             deallocate (packvals)
             deallocate (help_packvals)
             deallocate (help_rvals)
             
           endif
          
       endif
    endif
    
  end Function get_packed_var4

  !**********************************************************************
  ! A four-dimensional variable is packed and written to netcdf file
  !**********************************************************************
  function put_packed_var4 (values,nc_out,varid) result(status)

    use messy_clams_global, only: prec, dp, sp
    use netcdf
    use messy_clams_tools_packing   

    implicit none

    real(prec),dimension(:,:,:,:) :: values
    integer                      :: nc_out, varid, status

    real(prec)         :: missing_value = -1.e30
    integer            :: error
    integer(short_int) :: mdi_packed, fillvalue_packed
    real(prec)         :: scale, offset
    integer,            dimension(4)                :: dimlen
    integer(short_int), dimension(:,:,:,:), pointer :: packvals
    real(prec)        , dimension(:),       pointer :: help_vals
    integer(short_int), dimension(:),       pointer :: help_packvals

    ! get dimensions of 'values'
    dimlen(1) = size(values,1)
    dimlen(2) = size(values,2)
    dimlen(3) = size(values,3)
    dimlen(4) = size(values,4)

    ! allocate arrays
    allocate(packvals (dimlen(1),dimlen(2),dimlen(3),dimlen(4)))
    allocate(help_vals(dimlen(1)*dimlen(2)*dimlen(3)*dimlen(4)))
    allocate(help_packvals(dimlen(1)*dimlen(2)*dimlen(3)*dimlen(4)))
       
    ! transform to vector
    help_vals=reshape(values,(/dimlen(1)*dimlen(2)*dimlen(3)*dimlen(4)/))
       
    ! pack data
    call pack_array (help_vals, help_packvals, scale, offset, &
         mdi_packed, fillvalue_packed, error, &
         mdi=missing_value, fillvalue=NF90_FILL_FLOAT)
    
    ! transform vector with packed data to four-dimensional array
    packvals = reshape(help_packvals,(/dimlen(1),dimlen(2),dimlen(3),dimlen(4)/))
    status = nf90_put_var (nc_out,varid,packvals)
    if (status /= nf90_noerr) return
       
    ! write attributes (scale_factor, add_offset etc)
    status = nf90_redef(nc_out)
    if (status /= nf90_noerr) return
    status = write_packed_attributes (nc_out, varid, &
         scale, offset, fillvalue_packed, NF90_FILL_FLOAT,  &
         packmdi=mdi_packed, mdi=missing_value) 
    if (status /= nf90_noerr) return
    status = nf90_enddef(nc_out)
    if (status /= nf90_noerr) return
       
    deallocate(packvals)
    deallocate(help_vals)
    deallocate(help_packvals)

  end function put_packed_var4
 
  !***************************************************************************
  ! Read variable from input file and write packed variable to
  ! output file (data type "REAL"/"FLOAT").
  !
  ! Variables of data type "INTEGER" or "DOUBLE" are only copied !!!
  !
  ! Input:
  !   nc_in     : netCDF ID of input file
  !   nc_out    : netCDF ID of output file
  !   varid_in  : variable ID in input file
  !   xtype     : data type of unpacked variable
  !   ndims     : number of variable dimensions
  !   dimids    : dimension IDs
  !   maxerr    : ???
  !
  ! Output:
  !   status    : error code
  ! 
  !***************************************************************************
  function nc_pack_variable (nc_in,nc_out,varid_in,xtype,ndims,dimids,maxerr) &
       result(status) 

    use messy_clams_global, only: prec, dp, sp
    use netcdf
    use messy_clams_tools_packing   

    implicit none

    integer            :: nc_in, nc_out       
    integer            :: varid_in            
    integer            :: xtype, ndims        
    integer            :: dimids(ndims)       
    real(prec)         :: maxerr              
    integer            :: status
    
    ! local variables
    integer            :: dimlen(ndims)
    integer            :: error,  i, varid_out, natts, iatt
    integer(short_int) :: mdi_packed, fillvalue_packed
    real(prec)         :: scale, offset, mdi,  relerr
!!!!!
    real               :: fillvalue
    logical            :: mdi_exist
    CHARACTER(nf90_max_name) :: varname='', attname=''

    integer ,allocatable, dimension(:)          :: ivals1 !integer
    integer ,allocatable, dimension(:,:)        :: ivals2 
    integer ,allocatable, dimension(:,:,:)      :: ivals3 
    integer ,allocatable, dimension(:,:,:,:)    :: ivals4 
    integer ,allocatable, dimension(:,:,:,:,:)  :: ivals5 
    integer ,allocatable, dimension(:,:,:,:,:,:):: ivals6 

    real(dp), allocatable, dimension(:)           :: dvals1 !double
    real(dp), allocatable, dimension(:,:)         :: dvals2 
    real(dp), allocatable, dimension(:,:,:)       :: dvals3 
    real(dp), allocatable, dimension(:,:,:,:)     :: dvals4 
    real(dp), allocatable, dimension(:,:,:,:,:)   :: dvals5 
    real(dp), allocatable, dimension(:,:,:,:,:,:) :: dvals6 

    real(prec), dimension(:), pointer           :: rvals1 !real
    real(prec), dimension(:,:), pointer         :: rvals2 
    real(prec), dimension(:,:,:), pointer       :: rvals3 
    real(prec), dimension(:,:,:,:), pointer     :: rvals4 
    real(prec), dimension(:,:,:,:,:), pointer   :: rvals5 
    real(prec), dimension(:,:,:,:,:,:), pointer :: rvals6 

    integer(short_int), dimension(:), pointer           :: packvals1 
    integer(short_int), dimension(:,:), pointer         :: packvals2 
    integer(short_int), dimension(:,:,:), pointer       :: packvals3 
    integer(short_int), dimension(:,:,:,:), pointer     :: packvals4 
    integer(short_int), dimension(:,:,:,:,:), pointer   :: packvals5 
    integer(short_int), dimension(:,:,:,:,:,:), pointer :: packvals6 

    real(prec)              , dimension(:),   pointer :: help_rvals
    integer(short_int), dimension(:),   pointer :: help_packvals

    logical, dimension(:), pointer :: notmdi
   
    status = 0

    if (ndims > 6) then
       write (*,*) 'Die Variable hat mehr als 6 Dimensionen !'
       status = 99
       return
    endif

    ! get dimensions
    do i = 1, ndims
       status = nf90_inquire_dimension (nc_in,dimids(i),len=dimlen(i))
       if (status /= nf90_noerr) return 
    enddo

    status = nf90_inquire_variable (nc_in, varid_in, name=varname, natts=natts)
    IF (status /= nf90_noerr)  write (*,*) nf90_strerror(status)

    select case (xtype)

    ! copy integer variable
    case(nf90_int)

       ! define mode
       status = nf90_redef (nc_out)
       IF (status /= nf90_noerr)  write (*,*) nf90_strerror(status)

       ! define variable
       status = nf90_def_var (nc_out, varname, nf90_int, dimids(1:ndims), varid_out)
       if (status /= nf90_noerr) return 

       ! copy variable attributes
       DO iatt = 1, natts
          status = nf90_inq_attname (nc_in, varid_in, iatt, attname) 
          IF (status /= nf90_noerr)  write (*,*) nf90_strerror(status)
          status = nf90_copy_att (nc_in, varid_in, attname, nc_out, varid_out)
          IF (status /= nf90_noerr)  write (*,*) nf90_strerror(status)
       ENDDO

       ! end define mode
       status = nf90_enddef (nc_out)
       IF (status /= nf90_noerr)  write (*,*) nf90_strerror(status)
       
       ! copy values
       select case (ndims)
       case (1)
          allocate(ivals1(dimlen(1)))
          status = nf90_get_var (nc_in,varid_in,ivals1)
          if (status /= nf90_noerr) return 
          status = nf90_put_var (nc_out,varid_out,ivals1)
          if (status /= nf90_noerr) return 
          deallocate(ivals1)
       case (2)
          allocate(ivals2(dimlen(1),dimlen(2)))
          status = nf90_get_var (nc_in,varid_in,ivals2)
          if (status /= nf90_noerr) return
          status = nf90_put_var (nc_out,varid_out,ivals2)
          if (status /= nf90_noerr) return
          deallocate(ivals2)
       case (3)
          allocate(ivals3(dimlen(1),dimlen(2),dimlen(3)))
          status = nf90_get_var (nc_in,varid_in,ivals3)
          if (status /= nf90_noerr) return
          status = nf90_put_var (nc_out,varid_out,ivals3)
          if (status /= nf90_noerr) return
          deallocate(ivals3)
       case (4)
          allocate(ivals4(dimlen(1),dimlen(2),dimlen(3),dimlen(4)))
          status = nf90_get_var (nc_in,varid_in,ivals4)
          if (status /= nf90_noerr) return
          status = nf90_put_var (nc_out,varid_out,ivals4)
          if (status /= nf90_noerr) return
          deallocate(ivals4)
       case (5)
          allocate(ivals5(dimlen(1),dimlen(2),dimlen(3),dimlen(4),dimlen(5)))
          status = nf90_get_var (nc_in,varid_in,ivals5)
          if (status /= nf90_noerr) return
          status = nf90_put_var (nc_out,varid_out,ivals5)
          if (status /= nf90_noerr) return
          deallocate(ivals5)
       case (6)
          allocate(ivals6(dimlen(1),dimlen(2),dimlen(3),dimlen(4),dimlen(5),dimlen(6)))
          status = nf90_get_var (nc_in,varid_in,ivals6)
          if (status /= nf90_noerr) return
          status = nf90_put_var (nc_out,varid_out,ivals6)
          if (status /= nf90_noerr) return
          deallocate(ivals6)
       end select

    ! pack float variable
    case(nf90_float)

       ! allocate arrays
       ! read values from input file
       ! convert to vector 
       select case (ndims)
       case (1)
          allocate(rvals1(dimlen(1)))
          allocate(packvals1(dimlen(1)))
          allocate(help_packvals(dimlen(1)))
          allocate(help_rvals(dimlen(1)))
          status = nf90_get_var (nc_in,varid_in,rvals1)
          if (status /= nf90_noerr) return
          help_rvals=rvals1
       case (2)
          allocate(rvals2(dimlen(1),dimlen(2)))
          allocate(packvals2(dimlen(1),dimlen(2)))
          allocate(help_packvals(dimlen(1)*dimlen(2)))
          allocate(help_rvals(dimlen(1)*dimlen(2)))
          status = nf90_get_var (nc_in,varid_in,rvals2)
          if (status /= nf90_noerr) return
          help_rvals=reshape(rvals2,(/dimlen(1)*dimlen(2)/))
       case (3)
          allocate(rvals3(dimlen(1),dimlen(2),dimlen(3)))
          allocate(packvals3(dimlen(1),dimlen(2),dimlen(3)))
          allocate(help_packvals(dimlen(1)*dimlen(2)*dimlen(3)))
          allocate(help_rvals(dimlen(1)*dimlen(2)*dimlen(3)))
          status = nf90_get_var (nc_in,varid_in,rvals3)
          if (status /= nf90_noerr) return
          help_rvals=reshape(rvals3,(/dimlen(1)*dimlen(2)*dimlen(3)/))
       case (4)
          allocate(rvals4(dimlen(1),dimlen(2),dimlen(3),dimlen(4)))
          allocate(packvals4(dimlen(1),dimlen(2),dimlen(3),dimlen(4)))
          allocate(help_packvals(dimlen(1)*dimlen(2)*dimlen(3)*dimlen(4)))
          allocate(help_rvals(dimlen(1)*dimlen(2)*dimlen(3)*dimlen(4)))
          status = nf90_get_var (nc_in,varid_in,rvals4)
          if (status /= nf90_noerr) return
          help_rvals=reshape(rvals4,(/dimlen(1)*dimlen(2)*dimlen(3)*dimlen(4)/))
       case (5)
          allocate(rvals5(dimlen(1),dimlen(2),dimlen(3),dimlen(4),dimlen(5)))
          allocate(packvals5(dimlen(1),dimlen(2),dimlen(3),dimlen(4),dimlen(5)))
          allocate(help_packvals(dimlen(1)*dimlen(2)*dimlen(3)*dimlen(4)*dimlen(5)))
          allocate(help_rvals(dimlen(1)*dimlen(2)*dimlen(3)*dimlen(4)*dimlen(5)))
          status = nf90_get_var (nc_in,varid_in,rvals5)
          if (status /= nf90_noerr) return
          help_rvals=reshape(rvals5,(/dimlen(1)*dimlen(2)*dimlen(3)*dimlen(4)*dimlen(5)/))
       case (6)
          allocate(rvals6(dimlen(1),dimlen(2),dimlen(3),dimlen(4),dimlen(5),dimlen(6)))
          allocate(packvals6(dimlen(1),dimlen(2),dimlen(3),dimlen(4),dimlen(5),dimlen(6)))
          allocate(help_packvals(dimlen(1)*dimlen(2)*dimlen(3)*dimlen(4)&
               &*dimlen(5)*dimlen(6)))
          allocate(help_rvals(dimlen(1)*dimlen(2)*dimlen(3)*dimlen(4)*dimlen(5)*dimlen(6)))
          status = nf90_get_var (nc_in,varid_in,rvals6)
          if (status /= nf90_noerr) return
          help_rvals=reshape(rvals6,(/dimlen(1)*dimlen(2)*dimlen(3)*dimlen(4)*dimlen(5)*dimlen(6)/))
       end select
       
       ! get variable attributes "missing_value" and/or "_FillValue"
       mdi_exist = .false.
       status = nf90_get_att (nc_in,varid_in,"missing_value",mdi)
       if (status == nf90_noerr)  mdi_exist = .true.
       status = nf90_get_att (nc_in,varid_in,"_FillValue",fillvalue)
       if (status /= nf90_noerr)  fillvalue = NF90_FILL_FLOAT
       
       ! set valid array "notmdi"
       allocate (notmdi(size(help_rvals)))
       notmdi = .true.
       if (mdi_exist) then
          where (abs(help_rvals-mdi)<epsilon(mdi))
             notmdi = .false.
          end where
       endif
       where (abs(help_rvals-fillvalue)<epsilon(fillvalue))
          notmdi = .false.
       end where

       ! multiply by "scale_factor" if present 
       status = nf90_get_att (nc_in,varid_in,"scale_factor",scale)
       if (status == nf90_noerr) then
          where (notmdi)
             help_rvals = help_rvals*scale
          end where
       endif
       ! add "add_offset" if present
       status = nf90_get_att (nc_in,varid_in,"add_offset",offset)
       if (status == nf90_noerr) then
          where (notmdi)
             help_rvals = help_rvals+offset
          end where
       endif
       
       ! pack values
       if (mdi_exist) then
          call pack_array (help_rvals, help_packvals, scale, offset, &
               mdi_packed, fillvalue_packed, error, &
               mdi=mdi, fillvalue=fillvalue)
       else
          call pack_array (help_rvals, help_packvals, scale, offset, &
               mdi_packed, fillvalue_packed, error, fillvalue=fillvalue)
       endif
       
       !-----------------------------------------------
       ! FEHLERABSCHAETZUNG ???!!!
       !-----------------------------------------------
       call get_packing_accuracy (help_rvals, notmdi, relerr, scale, offset)

       deallocate (notmdi)

       ! Im Moment wird immer gepackt !!!!
       relerr = 0.
       if (relerr > maxerr) then ! copy values
          
          write (*,*) 'Rel. Fehler:',relerr
          write (*,*) 'Der Fehler ist groesser als ',maxerr, '!!!'
          write (*,*) 'Die Variable wird nicht gepackt, sondern kopiert !!!'
          
          ! define mode
          status = nf90_redef (nc_out)
          IF (status /= nf90_noerr)  write (*,*) nf90_strerror(status)
          
          ! define variable
          status = nf90_def_var (nc_out, varname, nf90_float, dimids(1:ndims), varid_out)
          if (status /= nf90_noerr) return 
          
          ! copy variable attributes
          DO iatt = 1, natts
             status = nf90_inq_attname (nc_in, varid_in, iatt, attname) 
             IF (status /= nf90_noerr)  write (*,*) nf90_strerror(status)
             status = nf90_copy_att (nc_in, varid_in, attname, nc_out, varid_out)
             IF (status /= nf90_noerr)  write (*,*) nf90_strerror(status)
          ENDDO
          
          ! end define mode
          status = nf90_enddef (nc_out)
          IF (status /= nf90_noerr)  write (*,*) nf90_strerror(status)
          
          ! write array (unpacked) to output file
          select case (ndims)
          case (1)
             status = nf90_put_var (nc_out,varid_out,rvals1)
          case (2)
             status = nf90_put_var (nc_out,varid_out,rvals2)
          case (3)
             status = nf90_put_var (nc_out,varid_out,rvals3)
          case (4)
             status = nf90_put_var (nc_out,varid_out,rvals4)
          case (5)
             status = nf90_put_var (nc_out,varid_out,rvals5)
          case (6)
             status = nf90_put_var (nc_out,varid_out,rvals6)
          end select
          if (status /= nf90_noerr) return
          
       else  ! pack values

          ! define mode
          status = nf90_redef (nc_out)
          IF (status /= nf90_noerr)  write (*,*) nf90_strerror(status)
          
          ! define variable
          status = nf90_def_var (nc_out, varname, nf90_short, dimids(1:ndims), varid_out)
          if (status /= nf90_noerr) return 
          
          ! copy variable attributes
          DO iatt = 1, natts
             status = nf90_inq_attname (nc_in, varid_in, iatt, attname) 
             IF (status /= nf90_noerr)  write (*,*) nf90_strerror(status)
             if (attname /= '_FillValue') then
                status = nf90_copy_att (nc_in, varid_in, attname, nc_out, varid_out)
                IF (status /= nf90_noerr)  write (*,*) nf90_strerror(status)
             endif
          ENDDO
          
          ! end define mode
          status = nf90_enddef (nc_out)
          IF (status /= nf90_noerr)  write (*,*) nf90_strerror(status)
          
          ! write packed array to output file
          select case (ndims)
          case (1)
             packvals1 = help_packvals
             status = nf90_put_var (nc_out,varid_out,packvals1)
          case (2)
             packvals2 = reshape(help_packvals,(/dimlen(1),dimlen(2)/))
             status = nf90_put_var (nc_out,varid_out,packvals2)
          case (3)
             packvals3 = reshape(help_packvals,(/dimlen(1),dimlen(2),dimlen(3)/))
             status = nf90_put_var (nc_out,varid_out,packvals3)
          case (4)
             packvals4 = reshape(help_packvals,(/dimlen(1),dimlen(2),dimlen(3)&
                  &,dimlen(4)/))
             status = nf90_put_var (nc_out,varid_out,packvals4)
          case (5)
             packvals5 = reshape(help_packvals,(/dimlen(1),dimlen(2),dimlen(3)&
                  &,dimlen(4),dimlen(5)/))
             status = nf90_put_var (nc_out,varid_out,packvals5)
          case (6)
             packvals6 = reshape(help_packvals,(/dimlen(1),dimlen(2),dimlen(3)&
                  &,dimlen(4),dimlen(5),dimlen(6)/))
             status = nf90_put_var (nc_out,varid_out,packvals6)
          end select
          
          if (status /= nf90_noerr) return
          
          ! write attributes (scale_factor,add_offset etc.)
          status = nf90_redef(nc_out)
          if (mdi_exist) then
             status = write_packed_attributes (nc_out, varid_out, scale, offset, &
                  fillvalue_packed, fillvalue, packmdi=mdi_packed, mdi=mdi) 
          else
             status = write_packed_attributes (nc_out, varid_out, scale, offset, &
                  fillvalue_packed, fillvalue) 
          endif
          status = nf90_enddef(nc_out)
          
       endif

       ! deallocate arrays
       select case (ndims)
       case (1)
          deallocate(rvals1, packvals1)
       case (2)
          deallocate(rvals2, packvals2)
       case (3)
          deallocate(rvals3, packvals3)
       case (4)
          deallocate(rvals4, packvals4)
       case (5)
          deallocate(rvals5, packvals5)
       case (6)
          deallocate(rvals6, packvals6)
       end select
       deallocate(help_packvals)
       deallocate(help_rvals)


    ! copy double variable
    case(nf90_double)

       ! define mode
       status = nf90_redef (nc_out)
       IF (status /= nf90_noerr)  write (*,*) nf90_strerror(status)

       ! define variable
       ! im Moment werden DOUBLE-Werte nur kopiert ! (spaeter: => INT*4)
       status = nf90_def_var (nc_out, varname, nf90_double, dimids(1:ndims), varid_out)
       !status = nf90_def_var (nc_out, varname, nf90_int, dimids(1:ndims),
       ! varid_out)
       if (status /= nf90_noerr) return 
       
       ! copy variable attributes
       DO iatt = 1, natts
          status = nf90_inq_attname (nc_in, varid_in, iatt, attname) 
          IF (status /= nf90_noerr)  write (*,*) nf90_strerror(status)
          status = nf90_copy_att (nc_in, varid_in, attname, nc_out, varid_out)
          IF (status /= nf90_noerr)  write (*,*) nf90_strerror(status)
       ENDDO
         
       ! end define mod
       status = nf90_enddef (nc_out)
       IF (status /= nf90_noerr)  write (*,*) nf90_strerror(status)

       ! copy values to output file
       select case (ndims)
       case (1)
          allocate(dvals1(dimlen(1)))
          status = nf90_get_var (nc_in,varid_in,dvals1)
          if (status /= nf90_noerr) return
          status = nf90_put_var (nc_out,varid_out,dvals1)
          if (status /= nf90_noerr) return
          deallocate(dvals1)
       case (2)
          allocate(dvals2(dimlen(1),dimlen(2)))
          status = nf90_get_var (nc_in,varid_in,dvals2)
          if (status /= nf90_noerr) return
          status = nf90_put_var (nc_out,varid_out,dvals2)
          if (status /= nf90_noerr) return
          deallocate(dvals2)
       case (3)
          allocate(dvals3(dimlen(1),dimlen(2),dimlen(3)))
          status = nf90_get_var (nc_in,varid_in,dvals3)
          if (status /= nf90_noerr) return
          status = nf90_put_var (nc_out,varid_out,dvals3)
          if (status /= nf90_noerr) return
          deallocate(dvals3)
       case (4)
          allocate(dvals4(dimlen(1),dimlen(2),dimlen(3),dimlen(4)))
          status = nf90_get_var (nc_in,varid_in,dvals4)
          if (status /= nf90_noerr) return
          status = nf90_put_var (nc_out,varid_out,dvals4)
          if (status /= nf90_noerr) return
          deallocate(dvals4)
       case (5)
          allocate(dvals5(dimlen(1),dimlen(2),dimlen(3),dimlen(4),dimlen(5)))
          status = nf90_get_var (nc_in,varid_in,dvals5)
          if (status /= nf90_noerr) return
          status = nf90_put_var (nc_out,varid_out,dvals5)
          if (status /= nf90_noerr) return
          deallocate(dvals5)
       case (6)
          allocate(dvals6(dimlen(1),dimlen(2),dimlen(3),dimlen(4),dimlen(5),dimlen(6)))
          status = nf90_get_var (nc_in,varid_in,dvals6)
          if (status /= nf90_noerr) return
          status = nf90_put_var (nc_out,varid_out,dvals6)
          if (status /= nf90_noerr) return
          deallocate(dvals6)
       end select
   
    ! invalid data type
    case default
       status = 99
    end select
    
  end function nc_pack_variable

  !***************************************************************************
  ! Read packed variable from input file and write unpacked variable to
  ! output file (data type "REAL"/"FLOAT")
  ! Variables of data type "INTEGER" or "DOUBLE" are only copied !!!
  !
  ! Input:
  !   nc_in     : netCDF ID of input file
  !   nc_out    : netCDF ID of output file
  !   varid_in  : variable ID in input file
  !   varid_out : variable ID in output file
  !   xtype     : data type of unpacked variable
  !   ndims     : number of variable dimensions
  !   dimids    : dimension IDs 
 !
  ! Output:
  !   status    : error code
  ! 
  !***************************************************************************
  function nc_unpack_variable (nc_in,nc_out,varid_in,varid_out, &
                                xtype,ndims, dimids) result(status) 

    use messy_clams_global, only: prec, dp
    use netcdf
    use messy_clams_tools_packing   

    implicit none

    integer            :: nc_in, nc_out       
    integer            :: varid_in,varid_out  
    integer            :: xtype, ndims        
    integer            :: dimids(ndims)       
    integer            :: status

    ! local variables
    integer            :: dimlen(ndims)
    integer            :: error, i
    integer(short_int) :: mdi_packed, fillvalue_packed
    real(prec)         :: scale, offset, mdi
!!!!!
    real               :: fillvalue
    logical            :: mdi_exist, fv_exist

    integer        ,allocatable,dimension(:)           :: ivals1 
    integer        ,allocatable,dimension(:,:)         :: ivals2 
    integer        ,allocatable,dimension(:,:,:)       :: ivals3 
    integer        ,allocatable,dimension(:,:,:,:)     :: ivals4 
    integer        ,allocatable,dimension(:,:,:,:,:)   :: ivals5 
    integer        ,allocatable,dimension(:,:,:,:,:,:) :: ivals6

    real(dp), allocatable,dimension(:)           :: dvals1 
    real(dp), allocatable,dimension(:,:)         :: dvals2 
    real(dp), allocatable,dimension(:,:,:)       :: dvals3 
    real(dp), allocatable,dimension(:,:,:,:)     :: dvals4 
    real(dp), allocatable,dimension(:,:,:,:,:)   :: dvals5 
    real(dp), allocatable,dimension(:,:,:,:,:,:) :: dvals6 

    real(prec), dimension(:), pointer           :: rvals1 
    real(prec), dimension(:,:), pointer         :: rvals2 
    real(prec), dimension(:,:,:), pointer       :: rvals3 
    real(prec), dimension(:,:,:,:), pointer     :: rvals4 
    real(prec), dimension(:,:,:,:,:), pointer   :: rvals5 
    real(prec), dimension(:,:,:,:,:,:), pointer :: rvals6 

    integer(short_int), dimension(:), pointer           :: packvals1 
    integer(short_int), dimension(:,:), pointer         :: packvals2 
    integer(short_int), dimension(:,:,:), pointer       :: packvals3 
    integer(short_int), dimension(:,:,:,:), pointer     :: packvals4 
    integer(short_int), dimension(:,:,:,:,:), pointer   :: packvals5 
    integer(short_int), dimension(:,:,:,:,:,:), pointer :: packvals6 

    real(prec)              , dimension(:),   pointer :: help_rvals
    integer(short_int), dimension(:),   pointer :: help_packvals

    if (ndims > 6) then
       write (*,*) 'Die Variable hat mehr als 6 Dimensionen !'
       status = 99
       return
    endif

    ! get dimensions
    do i = 1, ndims
       status = nf90_inquire_dimension (nc_in,dimids(i),len=dimlen(i))
    end do
   
    select case (xtype)

    ! copy integer variable 
    case(nf90_int)

       select case (ndims)
       case (1)
          allocate(ivals1(dimlen(1)))
          status = nf90_get_var (nc_in,varid_in,ivals1)
          if (status /= nf90_noerr) return
          status = nf90_put_var (nc_out,varid_out,ivals1)
          if (status /= nf90_noerr) return
          deallocate(ivals1)
       case (2)
          allocate(ivals2(dimlen(1),dimlen(2)))
          status = nf90_get_var (nc_in,varid_in,ivals2)
          if (status /= nf90_noerr) return
          status = nf90_put_var (nc_out,varid_out,ivals2)
          if (status /= nf90_noerr) return
          deallocate(ivals2)
       case (3)
          allocate(ivals3(dimlen(1),dimlen(2),dimlen(3)))
          status = nf90_get_var (nc_in,varid_in,ivals3)
          if (status /= nf90_noerr) return
          status = nf90_put_var (nc_out,varid_out,ivals3)
          if (status /= nf90_noerr) return
          deallocate(ivals3)
       case (4)
          allocate(ivals4(dimlen(1),dimlen(2),dimlen(3),dimlen(4)))
          status = nf90_get_var (nc_in,varid_in,ivals4)
          if (status /= nf90_noerr) return
          status = nf90_put_var (nc_out,varid_out,ivals4)
          if (status /= nf90_noerr) return
          deallocate(ivals4)
       case (5)
          allocate(ivals5(dimlen(1),dimlen(2),dimlen(3),dimlen(4),dimlen(5)))
          status = nf90_get_var (nc_in,varid_in,ivals5)
          if (status /= nf90_noerr) return
          status = nf90_put_var (nc_out,varid_out,ivals5)
          if (status /= nf90_noerr) return
          deallocate(ivals5)
       case (6)
          allocate(ivals6(dimlen(1),dimlen(2),dimlen(3),dimlen(4),dimlen(5),dimlen(6)))
          status = nf90_get_var (nc_in,varid_in,ivals6)
          if (status /= nf90_noerr) return
          status = nf90_put_var (nc_out,varid_out,ivals6)
          if (status /= nf90_noerr) return
          deallocate(ivals6)
       end select
       
    ! unpack float variable
    case(nf90_float)
       
       ! allocate arrays
       ! read values from input file
       ! convert to one-dimensional array (help_packvals)
       select case (ndims)
       case (1)
          allocate(rvals1(dimlen(1)))
          allocate(packvals1(dimlen(1)))
          allocate(help_packvals(dimlen(1)))
          allocate(help_rvals(dimlen(1)))
          status = nf90_get_var (nc_in,varid_in,packvals1)
          if (status /= nf90_noerr) return
          help_packvals=packvals1
       case (2)
          allocate(rvals2(dimlen(1),dimlen(2)))
          allocate(packvals2(dimlen(1),dimlen(2)))
          allocate(help_packvals(dimlen(1)*dimlen(2)))
          allocate(help_rvals(dimlen(1)*dimlen(2)))
          status = nf90_get_var (nc_in,varid_in,packvals2)
          if (status /= nf90_noerr) return
          help_packvals=reshape(packvals2,(/dimlen(1)*dimlen(2)/))
       case (3)
          allocate(rvals3(dimlen(1),dimlen(2),dimlen(3)))
          allocate(packvals3(dimlen(1),dimlen(2),dimlen(3)))
          allocate(help_packvals(dimlen(1)*dimlen(2)*dimlen(3)))
          allocate(help_rvals(dimlen(1)*dimlen(2)*dimlen(3)))
          status = nf90_get_var (nc_in,varid_in,packvals3)
          if (status /= nf90_noerr) return
          help_packvals=reshape(packvals3,(/dimlen(1)*dimlen(2)*dimlen(3)/))
       case (4)
          allocate(rvals4(dimlen(1),dimlen(2),dimlen(3),dimlen(4)))
          allocate(packvals4(dimlen(1),dimlen(2),dimlen(3),dimlen(4)))
          allocate(help_packvals(dimlen(1)*dimlen(2)*dimlen(3)*dimlen(4)))
          allocate(help_rvals(dimlen(1)*dimlen(2)*dimlen(3)*dimlen(4)))
          status = nf90_get_var (nc_in,varid_in,packvals4)
          if (status /= nf90_noerr) return
          help_packvals=reshape(packvals4,(/dimlen(1)*dimlen(2)*dimlen(3)*dimlen(4)/))
       case (5)
          allocate(rvals5(dimlen(1),dimlen(2),dimlen(3),dimlen(4),dimlen(5)))
          allocate(packvals5(dimlen(1),dimlen(2),dimlen(3),dimlen(4),dimlen(5)))
          allocate(help_packvals(dimlen(1)*dimlen(2)*dimlen(3)*dimlen(4)*dimlen(5)))
          allocate(help_rvals(dimlen(1)*dimlen(2)*dimlen(3)*dimlen(4)*dimlen(5)))
          status = nf90_get_var (nc_in,varid_in,packvals5)
          if (status /= nf90_noerr) return
          help_packvals=reshape(packvals5,(/dimlen(1)*dimlen(2)*dimlen(3)*dimlen(4)*dimlen(5)/))
       case (6)
          allocate(rvals6(dimlen(1),dimlen(2),dimlen(3),dimlen(4),dimlen(5),dimlen(6)))
          allocate(packvals6(dimlen(1),dimlen(2),dimlen(3),dimlen(4),dimlen(5),dimlen(6)))
          allocate(help_packvals(dimlen(1)*dimlen(2)*dimlen(3)*dimlen(4)*dimlen(5)*dimlen(6)))
          allocate(help_rvals(dimlen(1)*dimlen(2)*dimlen(3)*dimlen(4)*dimlen(5)*dimlen(6)))
          status = nf90_get_var (nc_in,varid_in,packvals6)
          if (status /= nf90_noerr) return
          help_packvals=reshape(packvals6,(/dimlen(1)*dimlen(2)*dimlen(3)*dimlen(4)*dimlen(5)*dimlen(6)/))
       end select
       
       ! get attributes "scale_factor" and "add_offset"
       status = nf90_get_att (nc_in,varid_in,"scale_factor", scale)
       if (status /= nf90_noerr) return
       status = nf90_get_att (nc_in,varid_in,"add_offset", offset)
       if (status /= nf90_noerr) return
       
       ! get attributes "missing_value" and "UNPACK_missing_value" 
       mdi_exist = .false.
       status = nf90_get_att (nc_in,varid_in,"missing_value",mdi_packed)
       if (status == nf90_noerr) then
          status = nf90_get_att (nc_in,varid_in,"UNPACK_missing_value",mdi)
          if (status == nf90_noerr) then
             mdi_exist = .true.
          else
             write (*,*) 'Das Attribut missing_value ist vorhanden '
             write (*,*) 'aber UNPACK_missing_value nicht !!!'
             return
          endif
       endif

       ! get attributes "_FillValue" and "UNPACK_FillValue" 
       fv_exist = .false.
       status = nf90_get_att (nc_in,varid_in,"_FillValue",fillvalue_packed)
       if (status == nf90_noerr) then
          status = nf90_get_att (nc_in,varid_in,"UNPACK_FillValue",fillvalue)
          if (status == nf90_noerr) then
             fv_exist = .true.
          else
             write (*,*) 'Das Attribut _FillValue ist vorhanden '
             write (*,*) 'aber UNPACK_FillValue nicht !!!'
             return
          endif
       endif

       if (.not. mdi_exist .and. .not. fv_exist) then
          write (*,*) 'Es muss mindestens eines der Attribute '
          write (*,*) '"missing_value" oder "_FillValue" vorhanden sein !!!'
          return
       endif
       
       ! unpack array
       if (mdi_exist .and. fv_exist) then
          call unpack_array (help_rvals, help_packvals, scale, offset, error, &
               mdi=mdi, mdi_packed=mdi_packed, &
               fillvalue=fillvalue, fillvalue_packed=fillvalue_packed)
       elseif (mdi_exist) then
          call unpack_array (help_rvals, help_packvals, scale, offset, error, &
               mdi=mdi, mdi_packed=mdi_packed)
       else
          call unpack_array (help_rvals, help_packvals, scale, offset, error, &
               fillvalue=fillvalue, fillvalue_packed=fillvalue_packed)
       endif
       
       ! write unpacked array to output file
       select case (ndims)
       case (1)
          rvals1 = help_rvals
          status = nf90_put_var (nc_out,varid_out,rvals1)
       case (2)
          rvals2 = reshape(help_rvals,(/dimlen(1),dimlen(2)/))
          status = nf90_put_var (nc_out,varid_out,rvals2)
       case (3)
          rvals3 = reshape(help_rvals,(/dimlen(1),dimlen(2),dimlen(3)/))
          status = nf90_put_var (nc_out,varid_out,rvals3)
       case (4)
          rvals4 = reshape(help_rvals,(/dimlen(1),dimlen(2),dimlen(3),dimlen(4)/))
          status = nf90_put_var (nc_out,varid_out,rvals4)
       case (5)
          rvals5 = reshape(help_rvals,(/dimlen(1),dimlen(2),dimlen(3),dimlen(4),dimlen(5)/))
          status = nf90_put_var (nc_out,varid_out,rvals5)
       case (6)
          rvals6 = reshape(help_rvals,(/dimlen(1),dimlen(2),dimlen(3),dimlen(4),dimlen(5),dimlen(6)/))
          status = nf90_put_var (nc_out,varid_out,rvals6)
       end select
       if (status /= nf90_noerr) return
       
       ! delete attributes: scale_factor, add_offset, etc.
       status = nf90_redef(nc_out)
       status = delete_packed_attributes (nc_out,varid_out) 
       status = nf90_enddef(nc_out)
       
       ! deallocate arrays
       select case (ndims)
       case (1)
          deallocate(rvals1, packvals1)
       case (2)
          deallocate(rvals2, packvals2)
       case (3)
          deallocate(rvals3, packvals3)
       case (4)
          deallocate(rvals4, packvals4)
       case (5)
          deallocate(rvals5, packvals5)
       case (6)
          deallocate(rvals6, packvals6)
       end select
       deallocate(help_packvals)
       deallocate(help_rvals)

    ! copy variable of type "double"
    case(nf90_double)

       select case (ndims)
       case (1)
          allocate(dvals1(dimlen(1)))
          status = nf90_get_var (nc_in,varid_in,dvals1)
          if (status /= nf90_noerr) return
          status = nf90_put_var (nc_out,varid_out,dvals1)
          if (status /= nf90_noerr) return
          deallocate(dvals1)
       case (2)
          allocate(dvals2(dimlen(1),dimlen(2)))
          status = nf90_get_var (nc_in,varid_in,dvals2)
          if (status /= nf90_noerr) return
          status = nf90_put_var (nc_out,varid_out,dvals2)
          if (status /= nf90_noerr) return
          deallocate(dvals2)
       case (3)
          allocate(dvals3(dimlen(1),dimlen(2),dimlen(3)))
          status = nf90_get_var (nc_in,varid_in,dvals3)
          if (status /= nf90_noerr) return
          status = nf90_put_var (nc_out,varid_out,dvals3)
          if (status /= nf90_noerr) return
          deallocate(dvals3)
       case (4)
          allocate(dvals4(dimlen(1),dimlen(2),dimlen(3),dimlen(4)))
          status = nf90_get_var (nc_in,varid_in,dvals4)
          if (status /= nf90_noerr) return
          status = nf90_put_var (nc_out,varid_out,dvals4)
          if (status /= nf90_noerr) return
          deallocate(dvals4)
       case (5)
          allocate(dvals5(dimlen(1),dimlen(2),dimlen(3),dimlen(4),dimlen(5)))
          status = nf90_get_var (nc_in,varid_in,dvals5)
          if (status /= nf90_noerr) return
          status = nf90_put_var (nc_out,varid_out,dvals5)
          if (status /= nf90_noerr) return
          deallocate(dvals5)
       case (6)
          allocate(dvals6(dimlen(1),dimlen(2),dimlen(3),dimlen(4),dimlen(5),dimlen(6)))
          status = nf90_get_var (nc_in,varid_in,dvals6)
          if (status /= nf90_noerr) return
          status = nf90_put_var (nc_out,varid_out,dvals6)
          if (status /= nf90_noerr) return
          deallocate(dvals6)
       end select
       
    case default
       status = 99
    end select
    
  end function nc_unpack_variable


End Module messy_clams_tools_packutils

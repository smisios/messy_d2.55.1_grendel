!******************************************************************************
! File utils/src/packing.f90
!
! This file is part of the Chemical Lagrangian Model of the
! Stratosphere (CLaMS). CLaMS is a hierarchy of numerical models and
! preprocessors which simulate transport, chemical reactions and
! mixing processes in the stratosphere.  
!
! Copyright (c) 2006
! N. Thomas, Forschungszentrum Juelich GmbH
! Last Modified By: N.Thomas
! Last Modified On: Tue Jul 26 10:31:04 2011
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
! Module packing
!
! contains the following subroutines:                                      
!
!     subroutine pack_array (x, x_packed, scale, offset, mdi_packed, &
!                         fillvalue_packed, error, mdi, fillvalue)
!     subroutine unpack_array (x, x_packed, scale, offset, error, &
!                       mdi, mdi_packed, fillvalue, fillvalue_packed)
!     
!
!***********************************************************************
Module messy_clams_tools_packing

  integer, parameter :: short_int = SELECTED_INT_KIND(4)

  integer(short_int), parameter :: max_int=32767 
  integer(short_int), parameter :: min_int=-32767

  integer, parameter :: range=65532

contains

  !*******************************************************************************
  ! Pack one-dimensional REAL array to SHORT INT
  !
  ! Input:
  !   x         : one-dimensional REAL array
  !   mdi       : missing value (optional)
  !   fillvalue : fillvalue (optional)
  !
  !   mdi or fillvalue must be present!
  !
  ! Output:
  !   x_packed         : one-dimensional SHORT INT array
  !   scale            : scale factor
  !   offset           : offset
  !   mdi_packed       : missing value for packed data
  !   fillvalue_packed : fillvalue for packed data
  !   error            : error code
  ! 
  !*******************************************************************************
  subroutine pack_array (x, x_packed, scale, offset, mdi_packed, &
                         fillvalue_packed, error, mdi, fillvalue)
    
    use messy_clams_global, only: prec

    implicit none
    
    real(prec),         dimension(:), pointer   :: x          ! input   
    integer(short_int), dimension(:), pointer   :: x_packed   ! output
    real(prec),         intent(out)             :: scale, offset 
    integer(short_int), intent(out)             :: mdi_packed, fillvalue_packed
    integer,            intent(out)             :: error
    real(prec), optional, intent(in)            :: mdi
!!!!!
    real,       optional, intent(in)            :: fillvalue
    
    real(prec) :: min_x, max_x
    logical, allocatable, dimension(:) :: notmdi
    
    error = 0
    mdi_packed = max_int
    fillvalue_packed = min_int
       
    ! mdi or fillvalue must be present
    if (.not. present(mdi) .and. .not. present(fillvalue)) then
       error = 1
       
    else

       ! filter missing values 
       allocate (notmdi(size(x)))
       notmdi = .true.
       if (present(mdi)) then
          where (abs(x-mdi)<epsilon(mdi))
             notmdi = .false.
             x_packed = mdi_packed
          end where
       endif
       if (present(fillvalue)) then
          where (abs(x-fillvalue)<epsilon(fillvalue))
             notmdi = .false.
             x_packed = fillvalue_packed
          end where
       endif
       
       ! get minimum and maximum
       min_x  = minval(x,mask=notmdi)
       max_x  = maxval(x,mask=notmdi)
       
       ! calculate offset and scale
       offset = (max_x + min_x) * 0.5
       scale = (max_x - min_x) / range

       ! if all values are equal: set x_packed=0. (offset=value)
       if (scale == 0.) then
          x_packed = 0.

       else
          where (notmdi)
             x_packed = nint((x - offset) / scale)
          end where
       endif
       
       deallocate (notmdi)
     
    endif
  
  end subroutine pack_array
     
  
  !*******************************************************************************
  ! Unack one-dimensional SHORT INT array to REAL
  !
  ! Input:
  !   x_packed         : one-dimensional SHORT INT array
  !   scale            : scale factor
  !   offset           : offset
  !   mdi              : missing value for unpacked data (optional)
  !   fillvalue        : fillvalue for unpacked data (optional)
  !   mdi_packed       : missing value for packed data (optional)
  !   fillvalue_packed : fillvalue for packed data (optional)
  !
  !   mdi or fillvalue must be present!
  !   if mdi is present, mdi_packed must be present!
  !   if fillvalue is present, fillvalue_packed must be present!
  !
  ! Output:
  !   x                : one-dimensional REAL array
  !   error            : error code
  ! 
  !*******************************************************************************
  subroutine unpack_array (x, x_packed, scale, offset, error, mdi, mdi_packed, &
                           fillvalue, fillvalue_packed)
    
    use messy_clams_global, only: prec

    implicit none
    
    real(prec),         dimension(:), pointer :: x         ! output
    integer(short_int), dimension(:), pointer :: x_packed  ! input
    real(prec),         intent(in)            :: scale, offset
    integer,            intent(out)           :: error
    real(prec), optional, intent(in)          :: mdi
!!!!!
    real,       optional, intent(in)          :: fillvalue
    integer(short_int), optional, intent(in)  :: mdi_packed, fillvalue_packed
    
    logical, allocatable, dimension(:) :: notmdi

    error = 0

    ! mdi or fillvalue must be present
    if (.not. present(mdi) .and. .not. present(fillvalue)) then
       error = 1

    !if mdi is present, mdi_packed must be present
    elseif (present(mdi) .and. .not. present(mdi_packed)) then
       error = 2

    ! if fillvalue is present, fillvalue_packed must be present
    elseif (present(fillvalue) .and. .not. present(fillvalue_packed)) then
       error = 3

    else

       allocate (notmdi(size(x)))
       notmdi = .true.

       ! filter missing values
       if (present(mdi_packed)) then
          where (x_packed == mdi_packed)
             notmdi = .false.
             x = mdi
          end where
       endif
       if (present(fillvalue_packed)) then
          where (x_packed == fillvalue_packed)
             notmdi = .false.
             x = fillvalue
          end where
       endif
       
       ! calculate unpacked values
       where (notmdi)
          x = x_packed * scale + offset
       end where

       deallocate (notmdi)

    endif
    
  end subroutine unpack_array
  
  !*******************************************************************************
  ! 
  !*******************************************************************************
 ! subroutine pack_array_dbl (x, x_packed, scale, offset, mdi_packed, &
 !                            fillvalue_packed, error, mdi, fillvalue)

 !   implicit none
    
 !   real(kind(0d0)), dimension(:), pointer   :: x          ! input   
 !   integer,         dimension(:), pointer   :: x_packed   ! output
 !   real(kind(0d0)), intent(out)             :: scale, offset 
 !   integer, intent(out)                     :: mdi_packed, fillvalue_packed
 !   integer, intent(out)                     :: error
 !   real(kind(0d0)), optional, intent(in)    :: mdi, fillvalue


 ! end subroutine pack_array_dbl

  !*******************************************************************************
  ! 
  !*******************************************************************************
   !subroutine unpack_array_dbl

  ! end subroutine unpack_array_dbl


end Module messy_clams_tools_packing
  

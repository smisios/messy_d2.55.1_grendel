!********************************************************************************!
! File clams/dynmod/source/ap_m_access.f90
!
! This file is part of the Chemical Lagrangian Model of the
! Stratosphere (CLaMS). CLaMS is a hierarchy of numerical models and
! preprocessors which simulate transport, chemical reactions and
! mixing processes in the stratosphere.  
!
! Copyright (c) 2006
! Juergen Ankenbrand, Paul Konopka, Nicole Thomas
! Forschungszentrum Juelich GmbH
! Last Modified By: N.Thomas
! Last Modified On: Mon Sep 28 12:54:18 2015
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
!***********************************************************************
!
!   Module  ap_m_access   
!
!   Data-Modul for access on the grid adaption data structure.
!
!   Main-data-structur: ap_struct 
!
!   Public functions on this data-structur
!
!      function object(feld,dimension)                                    
!      function index(object,i)
!      function lon(object,index)                                        
!      function lat(object,index)                                        
!      function lev(object,index)                                      
!      subroutine define_lev(object,index,lev)
!      function time_init(object,index)                                  
!      function coor(object,index)                                       
!      function coor_old(object,index)                                   
!      function state(object,index)                                      
!      function state_vert(object,index)                                      
!      function subset(object,index)
!      subroutine define_state(object,index, state)                      
!      subroutine define_state_vert(object,index, state)                       
!      subroutine reset_state(object)
!      subroutine reset_state_vert(object)
!      function tracer(object,index)                                     
!      subroutine define_tracer(object,index,tracer)                     
!      function is_c(object)
!      function size_of_c(object)
!      function c(object,index,tracer)                                   
!      subroutine c_all(object,index,c)
!      subroutine define_c_all(object,index,c)
!      function get_nb(object,index)                                     
!      subroutine nb(object,index,neighbours)                            
!      subroutine get_ind(nbfeld,object,index)                           
!      subroutine ind(object,index,nbfeld)                               
!      subroutine free(object)                                           
!                                                                        
!      For further informations, see inline commentars.                  
!                                                                      
!   Operators on this data-structure:                                  
!                                                                      
!      .index.                                                         
!                                                                      
!   used modules:                                                      
!                                                                      
!      types_m   :  data-structures: ap, adapt_set                     
!                                                                      
!***********************************************************************

module messy_clamsmix_ap_m_access

  use messy_clams_global,    only: prec, dp
  use messy_clamsmix_global, only: ap, adapt_set, ctrl_out

  implicit none

  ! non-visible functions of the module
  private :: index
  private :: time_init
  private :: define_state                   
  private :: tracer                             
  private :: define_tracer      
  private :: get_nb                       
  private :: get_ind            


  ! visible functions of the module
  public :: lat
  public :: lon
  public :: lev
  public :: state                      
  public :: state_vert                      
  public :: object
  public :: coor                             
  public :: coor_old                      
  public :: free                                        
  public :: nb                      
  public :: ind  
  public :: reset_state 
  public :: c                 

  ! global variables of the module
  integer, public  :: n                  ! actual number of points(visible for user of the modul)

  !main structure(structure to realize an array of variable size)
  !contains 
  !  a pointer ('konstant')  to the array, which contains the start points of the adaption
  !           the array can be set with the funtion 'object'
  !  the pointer 'dynamic' points to an array, which contains the new points(if the first array is full) 
  !  some integer components to administrate the structure
  type,public :: ap_struct
     type(ap),dimension(:),pointer :: konstant      ! pointer to array with defined dimension 
     type(ap),dimension(:),pointer :: dynamic       ! pointer to array which often is resized
     integer                       :: konstant_size ! size of konstant array
     integer                       :: n             ! actual number of points
     integer                       :: max_points    ! actual maximal array size
     integer                       :: dimension     ! dimension of the used qhull
     logical                       :: old_coor      ! if old coordinates than 1
     real(prec)                    :: r_min         ! minimal initial length of APs on a unit sphere 
     real(prec)                    :: r_max         ! maximal initial length of APs on a unit sphere 
  end type ap_struct


  !operator to get one element of the array(realized by the structure)
  !with var_ap_struct.index.345 you get point 345 (if existing)
  interface operator(.index.) 
     module procedure index
  end interface

contains

  ! function to create an object of this data module
  ! 'feld' has to be an one-dimensional array of type ap
  function object(feld)
    type(ap),dimension(:),target,intent(in) :: feld
    type(ap_struct)                         :: object


    !initialize the structure
    object%konstant => feld
    nullify(object%dynamic)
    n=size(feld)
    object%konstant_size=n
    object%n=n
    object%max_points=n
  end function object



  ! private function for .index. operator
  function index(object,i)
    type(ap_struct),intent(in) :: object
    integer,intent(in)         :: i
    type(ap)                   :: index

    if ((i<=0) .or. (i>object%n)) then
       write(*,*) 'not in index section'
       stop
    else if (i<=object%konstant_size) then
       index = object%konstant(i)
    else
       index = object%dynamic(i-object%konstant_size)
    end if

  end function index

  !function to get the 'lon' component of the point 'index'
  function lon(object,index)
    type(ap_struct)    :: object
    integer,intent(in) :: index
    real(prec)         :: lon

    if ((index<=0) .or. (index>object%n)) then
       write(*,*) 'not in index section'
       stop
    else if (index<=object%konstant_size) then
       lon = object%konstant(index)%lon
    else
       lon = object%dynamic(index-object%konstant_size)%lon
    end if

  end function lon

  !function to get the 'lat' component of the point 'index'
  function lat(object,index)
    type(ap_struct)    :: object
    integer,intent(in) :: index
    real(prec)         :: lat

    if ((index<=0) .or. (index>object%n)) then
       write(*,*) 'not in index section'
       stop
    else if (index<=object%konstant_size) then
       lat = object%konstant(index)%lat
    else
       lat = object%dynamic(index-object%konstant_size)%lat
    end if

  end function lat

  !function to get the 'lev' component of the point 'index'
  function lev(object,index)
    type(ap_struct)    :: object
    integer,intent(in) :: index
    real(prec)         :: lev

    if ((index<=0) .or. (index>object%n)) then
       write(*,*) 'not in index section'
       stop
    else if (index<=object%konstant_size) then
       lev = object%konstant(index)%lev
    else
       lev = object%dynamic(index-object%konstant_size)%lev
    end if

  end function lev

  !function to define the 'lev' component of the point 'index'
  subroutine define_lev(object,index,lev)
    type(ap_struct)     :: object
    integer,intent(in)  :: index
    real(prec)          :: lev

    if ((index<=0) .or. (index>object%n)) then
       write(*,*) 'not in index section'
       stop
    else if (index<=object%konstant_size) then
       object%konstant(index)%lev = lev
    else
       object%dynamic(index-object%konstant_size)%lev = lev
    end if

  end subroutine define_lev

  !function to get the 'time_init' component of the point 'index'
  function time_init(object,index)
    type(ap_struct)    :: object
    integer,intent(in) :: index
    double precision   :: time_init

    if ((index<=0) .or. (index>object%n)) then
       write(*,*) 'not in index section'
       stop
    else if (index<=object%konstant_size) then
       time_init = object%konstant(index)%time_init
    else
       time_init = object%dynamic(index-object%konstant_size)%time_init
    end if

  end function time_init

  !function to get the 'coor' component of the point 'index'
  function coor(object,index)
    type(ap_struct)         :: object
    integer,intent(in)      :: index
    real(prec),dimension(3) :: coor

    if ((index<=0) .or. (index>object%n)) then
       write(*,*) 'not in index section'
       stop
    else if (index<=object%konstant_size) then
       coor(:) = object%konstant(index)%coor(:)
    else
       coor(:) = object%dynamic(index-object%konstant_size)%coor(:)
    end if

  end function coor

  !function to get the 'coor_old' component of the point 'index'
  function coor_old(object,index)
    type(ap_struct)         :: object
    integer,intent(in)      :: index
    real(prec),dimension(3) :: coor_old

    if ((index<=0) .or. (index>object%n)) then
       write(*,*) 'not in index section'
       stop
    else if (index<=object%konstant_size) then
       coor_old(:) = object%konstant(index)%coor_old(:)
    else
       coor_old(:) = object%dynamic(index-object%konstant_size)%coor_old(:)
    end if

  end function coor_old

  !function to get the state of the point 'index'
  function state(object,index)
    type(ap_struct)    :: object
    integer,intent(in) :: index
    integer            :: state

    if ((index<=0) .or. (index>object%n)) then
       write(*,*) 'not in index section'
       stop
    else if (index<=object%konstant_size) then
       state = object%konstant(index)%state
    else
       state = object%dynamic(index-object%konstant_size)%state
    end if

  end function state

  !function to get the state of the point 'index'
  function state_vert(object,index)
    type(ap_struct)    :: object
    integer,intent(in) :: index
    integer            :: state_vert

    if ((index<=0) .or. (index>object%n)) then
       write(*,*) 'not in index section'
       stop
    else if (index<=object%konstant_size) then
       state_vert = object%konstant(index)%state_vert
    else
       state_vert = object%dynamic(index-object%konstant_size)%state_vert
    end if

  end function state_vert

  !function to get the subset of the point 'index'
  function subset(object,index)
    type(ap_struct)    :: object
    integer,intent(in) :: index
    logical            :: subset

    if ((index<=0) .or. (index>object%n)) then
       write(*,*) 'not in index section'
       stop
    else if (index<=object%konstant_size) then
       subset = object%konstant(index)%subset
    else
       subset = object%dynamic(index-object%konstant_size)%subset
    end if

  end function subset

  subroutine define_state(object,index,state)
    type(ap_struct)     :: object
    integer,intent(in)  :: index
    integer             :: state

    if ((index<=0) .or. (index>object%n)) then
       write(*,*) 'not in index section'
       stop
    else if (index<=object%konstant_size) then
       object%konstant(index)%state = state
    else
       object%dynamic(index-object%konstant_size)%state = state
    end if

  end subroutine define_state

  subroutine define_state_vert(object,index,state)
    type(ap_struct)     :: object
    integer,intent(in)  :: index
    integer             :: state

    if ((index<=0) .or. (index>object%n)) then
       write(*,*) 'not in index section'
       stop
    else if (index<=object%konstant_size) then
       object%konstant(index)%state_vert = state
    else
       object%dynamic(index-object%konstant_size)%state_vert = state
    end if

  end subroutine define_state_vert

  subroutine reset_state(object)
    type(ap_struct)     :: object

    object%konstant%state = 0

    if (associated(object%dynamic)) object%dynamic%state=0

  end subroutine reset_state

  subroutine reset_state_vert(object)
    type(ap_struct)     :: object

    object%konstant%state_vert = 0

    if (associated(object%dynamic)) object%dynamic%state_vert=0

  end subroutine reset_state_vert


  !function to get the 'tracer' component of the point 'index'
  function tracer(object,index)
    type(ap_struct)    :: object
    integer,intent(in) :: index
    real(prec)         :: tracer

    if ((index<=0) .or. (index>object%n)) then
       write(*,*) 'not in index section'
       stop
    else if (index<=object%konstant_size) then
       tracer = object%konstant(index)%tracer
    else
       tracer = object%dynamic(index-object%konstant_size)%tracer
    end if

  end function tracer

  !subroutine to define the 'tracer' component of the point 'index'
  subroutine define_tracer(object,index,tracer)
    type(ap_struct)     :: object
    integer,intent(in)  :: index
    real(prec)          :: tracer

    if ((index<=0) .or. (index>object%n)) then
       write(*,*) 'not in index section'
       stop
    else if (index<=object%konstant_size) then
       object%konstant(index)%tracer = tracer
    else
       object%dynamic(index-object%konstant_size)%tracer = tracer
    end if

  end subroutine define_tracer


  !function which says, if c is associated
  function is_c(object)
    type(ap_struct)    :: object
    logical            :: is_c

    is_c = associated(object%konstant(1)%c)

  end function is_c

  !function to get the size of the 'chemie'- component
  function size_of_c(object)
    type(ap_struct)    :: object
    integer            :: size_of_c
    
    if (associated(object%konstant(1)%c)) then
       size_of_c= size(object%konstant(1)%c)
    else
       size_of_c = -1
    endif

  end function size_of_c

  !function to get the 'chemie' component i of the point 'index'
  function c(object,i,index)
    type(ap_struct)    :: object
    integer,intent(in) :: i,index
    real(prec)         :: c

    if ((index<=0) .or. (index>object%n)) then
       write(*,*) 'not in index section'
       stop
    else if (index<=object%konstant_size) then
       c = object%konstant(index)%c(i)
    else
       c = object%dynamic(index-object%konstant_size)%c(i)
    end if

  end function c

  !function to get the 'chemie' component of the point 'index'
  subroutine c_all(object,index,c)
    type(ap_struct)         :: object
    integer,intent(in)      :: index
    real(prec),dimension(:) :: c

    if ((index<=0) .or. (index>object%n)) then
       write(*,*) 'not in index section'
       stop
    else if (index<=object%konstant_size) then
       c=object%konstant(index)%c
    else
       c=object%dynamic(index-object%konstant_size)%c
    end if

  end subroutine c_all

  !function to define the 'chemie' component of the point 'index'
  subroutine define_c_all(object,index,c)
    type(ap_struct)          :: object
    integer,intent(in)       :: index
    real(prec),dimension(:)  :: c

    if ((index<=0) .or. (index>object%n)) then
       write(*,*) 'not in index section'
       stop
    else if (index<=object%konstant_size) then
       object%konstant(index)%c=c
    else
       object%dynamic(index-object%konstant_size)%c=c
    end if

  end subroutine define_c_all

  !function to get the 'nb' component of the point 'index', which contains
  !the number of neighbours of this point
  function get_nb(object,index)
    type(ap_struct)    :: object
    integer,intent(in) :: index
    integer            :: get_nb

    if ((index<=0) .or. (index>object%n)) then
       write(*,*) 'not in index section'
       stop
    else if (index<=object%konstant_size) then
       get_nb=object%konstant(index)%nb
    else
       get_nb=object%dynamic(index-object%konstant_size)%nb
    endif

  end function get_nb

  !subroutine to define the 'nb' component of the point 'index', which contains
  !the number of neighbours of this point
  subroutine nb(object,index,neighbours)
    type(ap_struct)    :: object
    integer,intent(in) :: index
    integer            :: neighbours

    if ((index<=0) .or. (index>object%n)) then
       write(*,*) 'not in index section'
       stop
    else if (index<=object%konstant_size) then
       object%konstant(index)%nb=neighbours
    else
       object%dynamic(index-object%konstant_size)%nb=neighbours       
    endif

  end subroutine nb

  !get object(index)%ind(:)
  !get the neigbours of the point 'index'
  !nbfeld must have the size object(index)%nb
  subroutine get_ind(nbfeld,object,index)
    type(ap_struct),intent(in)       :: object
    integer,intent(in)               :: index
    integer,dimension(:),intent(out) :: nbfeld

    if ((index<=0) .or. (index>object%n)) then
       write(*,*) 'not in index section'
       stop
    else if (index<=object%konstant_size) then
       nbfeld(:) = object%konstant(index)%ind(1:object%konstant(index)%nb)
    else
       nbfeld(:) = object%dynamic(index-object%konstant_size)%ind(1:object%dynamic(index-object%konstant_size)%nb)
    end if

  end subroutine get_ind

  !define object(index)%ind(:)
  !define the neigbours of the point 'index'
  !the size of object(index)%ind(:) is object(index)%nb(so this must be set before)
  subroutine ind(object,index,nbfeld)
    type(ap_struct)       :: object
    integer,intent(in)    :: index
    integer,dimension(:)  :: nbfeld

    if (index<=object%konstant_size) then
       if (associated(object%konstant(index)%ind)) &
            deallocate(object%konstant(index)%ind)
       allocate(object%konstant(index)%ind(object%konstant(index)%nb))
       object%konstant(index)%ind(1:object%konstant(index)%nb)=nbfeld(1:object%konstant(index)%nb)
    else
       if (associated(object%dynamic(index-object%konstant_size)%ind)) &
            deallocate(object%dynamic(index-object%konstant_size)%ind)
       allocate(object%dynamic(index-object%konstant_size)%ind(object%dynamic(index-object%konstant_size)%nb))
       object%dynamic(index-object%konstant_size)%ind(1:object%dynamic(index-object%konstant_size)%nb)= &
            nbfeld(1:object%dynamic(index-object%konstant_size)%nb)
    endif

  end subroutine ind

  !subroutine to get free the memory associated with the object
  subroutine free(object)
    type(ap_struct) :: object 

    integer         :: i

    !give free all neighbourlists of the konstant array
    do i=1,object%konstant_size
       if (associated(object%konstant(i)%ind)) then
          deallocate(object%konstant(i)%ind)
          nullify(object%konstant(i)%ind)
       endif
       if (associated(object%konstant(i)%c)) then
          deallocate(object%konstant(i)%c)
          nullify(object%konstant(i)%c)
       endif
    end do
    nullify(object%konstant)

    !give free all neighbourlists of the dynamic array and delete dynamic array
    if (associated(object%dynamic)) then 
       do i=1,object%max_points-object%konstant_size
          if (associated(object%dynamic(i)%ind)) then
             deallocate(object%dynamic(i)%ind)
             nullify(object%dynamic(i)%ind)
          endif
          if (associated(object%dynamic(i)%c)) then
             deallocate(object%dynamic(i)%c)
             nullify(object%dynamic(i)%c)
          endif
       end do
       deallocate(object%dynamic)
       nullify(object%dynamic)
    endif

  end subroutine free


end module messy_clamsmix_ap_m_access



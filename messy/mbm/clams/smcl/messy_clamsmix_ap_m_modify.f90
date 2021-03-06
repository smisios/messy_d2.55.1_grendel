!********************************************************************************!
! File clams/dynmod/source/ap_m_modify.f90
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
! Last Modified On: Mon May 18 10:20:23 2020
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
!   Module ap_m_modify                                                     
!                                                                       
!   Function-Module for grid adaption.                                   
!                                                                       
!   Main-data-structur: ap_struct 
!                                                                       
!   Public functions on this data-structur:                             
!                                                                       
!      function search_min(object,i,adapt_par,index,optim)              
!      function search_max(object,i,adapt_par,index,optim)              
!      function search_vertical_max(object,i,adapt_par,index,optim)
!      function search_bvf_min(object,i,adapt_par)            
!      function horiz_distance(object,i,j)                              
!      function lev_distance(object,i,j) 
!      function lev_distance_new(object,i,j) 
!      logical function lev_gt_min (object,index1,index2,l_min_act)
!      subroutine eliminate(object,i,index)     
!      subroutine new_point(object,i,index)                             
!      subroutine put_on_level(object,a_p,new)
!      subroutine interpolate(object)                                   
!      subroutine avg_concentration(object,i,index)
!      subroutine avg_concentration_all(object,i)
!      subroutine change_last_point(object,index)                       
!      subroutine change_neighbours_of_new_point(object,optim)          
!      subroutine change_neighbours(object,index)                       
!                                                                       
!   Private functions on this data-structur:  
!
!      function f_limit(object, i, j, adapt_par)
!      function f_min(object, i, j, adapt_par)
!      function f_max(object, i, j, adapt_par)
!      subroutine lose_concentration(object,i,index)
!      subroutine memory_management(object)
!      subroutine index_management(object,point)
!      subroutine get_neighbours (status,object,old1,index1,old2,index2,new)
!      function kart2ukmo(x)
! 
!                                                      
!   used modules:                                                       
!                                                                       
!      types_m     :  data-structures: ap, adapt_set                    
!      ap_m_access :  data-structures: ap_struct,n                      
!                                                                       
!***********************************************************************

module messy_clamsmix_ap_m_modify       

  use messy_clams_global,         only: prec
  use messy_clamsmix_global,      only: ap, adapt_set
  use messy_clamsmix_ap_m_access, only: n,ap_struct

  implicit none

  private:: f_limit, f_max, f_min, lose_concentration,index_management, &
            memory_management,get_neighbours, kart2ukmo

  integer,private :: num_new_points=100 ! number of new points, which will be added,
                                        ! if new points are needed
  integer,private :: num_new_ind=10     ! number of new neighbours, which will be added,
                                        ! if new neighbours are needed


contains

  ! function to set the maximal horizontal distance of two points
  ! APs for the adaptation procedure
  function f_limit(object, i, j, adapt_par, ilev)
    real(prec)                  :: f_limit
    type(ap_struct),intent(in)  :: object
    integer,intent(in)          :: i,j
    type(adapt_set), intent(in) :: adapt_par
    integer, intent(in)         :: ilev

    real(prec)                  :: lat_i, lat_j, lev_i, lev_j

    if (i<=object%konstant_size) then
       if (j<=object%konstant_size) then
          lat_i   = object%konstant(i)%lat
          lat_j   = object%konstant(j)%lat
          lev_i = object%konstant(i)%lev
          lev_j = object%konstant(j)%lev
       else
          lat_i   = object%konstant(i)%lat
          lat_j   = object%dynamic(j-object%konstant_size)%lat
          lev_i = object%konstant(i)%lev
          lev_j = object%dynamic(j-object%konstant_size)%lev
       endif
    else
       lat_i   = object%dynamic(i-object%konstant_size)%lat
       lat_j   = object%dynamic(j-object%konstant_size)%lat
       lev_i = object%dynamic(i-object%konstant_size)%lev
       lev_j = object%dynamic(j-object%konstant_size)%lev
    endif

    if ((lat_i >= adapt_par%lat_min) .and. (lat_i <= adapt_par%lat_max) .and. &         
        (lat_j >= adapt_par%lat_min) .and. (lat_j <= adapt_par%lat_max)) then 
       if ((lev_i<adapt_par%fac_limit_lev_down .and. &
            lev_j<adapt_par%fac_limit_lev_down) .or. &
           (lev_i>adapt_par%fac_limit_lev_up .and.  &
            lev_j>adapt_par%fac_limit_lev_up)) then
          f_limit=adapt_par%r_lim_h_outside(ilev)
       else
          f_limit=adapt_par%r_lim_h_inside(ilev)
       endif
    else  
       if ((lev_i<adapt_par%fac_limit_lev_down .and. &
            lev_j<adapt_par%fac_limit_lev_down) .or. &
           (lev_i>adapt_par%fac_limit_lev_up .and. & 
            lev_j>adapt_par%fac_limit_lev_up)) then
          f_limit=adapt_par%r_lim_c_outside(ilev)
       else
          f_limit=adapt_par%r_lim_c_inside(ilev)
       endif
    endif

  end function f_limit

  !function to set the minimal horiz_distance of two points
  !a point which is nearer to another point than this value, must be
  !eliminated
  function f_min(object, i, j, adapt_par, ilev)
    real(prec)                  :: f_min
    type(ap_struct),intent(in)  :: object
    integer,intent(in)          :: i,j
    type(adapt_set), intent(in) :: adapt_par
    integer, intent(in)         :: ilev

    real(prec)                  :: lat_i, lat_j

    if (i<=object%konstant_size) then
       if (j<=object%konstant_size) then
          lat_i=object%konstant(i)%lat
          lat_j=object%konstant(j)%lat
       else
          lat_i=object%konstant(i)%lat
          lat_j=object%dynamic(j-object%konstant_size)%lat
       endif
    else
       lat_i=object%dynamic(i-object%konstant_size)%lat
       lat_j=object%dynamic(j-object%konstant_size)%lat
    endif

    if ((lat_i >= adapt_par%lat_min) .and. (lat_i <= adapt_par%lat_max) .and. &         
        (lat_j >= adapt_par%lat_min) .and. (lat_j <= adapt_par%lat_max)) then 
       f_min=adapt_par%r_min_h(ilev)
    else
       f_min=adapt_par%r_min_c(ilev)
    endif

  end function f_min

  !function to set the maximal horiz_distance of two points
  !if the horiz_distance between two points i larger, a new point has to be created
  function f_max(object, i, j, adapt_par, ilev)
    real(prec)                  :: f_max
    type(ap_struct),intent(in)  :: object
    integer,intent(in)          :: i,j
    type(adapt_set), intent(in) :: adapt_par
    integer, intent(in)         :: ilev

    real(prec)                  :: lat_i, lat_j

    if (i<=object%konstant_size) then
       if (j<=object%konstant_size) then
          lat_i=object%konstant(i)%lat
          lat_j=object%konstant(j)%lat
       else
          lat_i=object%konstant(i)%lat
          lat_j=object%dynamic(j-object%konstant_size)%lat
       endif
    else
       lat_i=object%dynamic(i-object%konstant_size)%lat
       lat_j=object%dynamic(j-object%konstant_size)%lat
    endif

    if ((lat_i >= adapt_par%lat_min) .and. (lat_i <= adapt_par%lat_max) .and. &         
        (lat_j >= adapt_par%lat_min) .and. (lat_j <= adapt_par%lat_max)) then 
       f_max=adapt_par%r_max_h(ilev)
    else
       f_max=adapt_par%r_max_c(ilev)
    endif

  end function f_max


  !function to check if a point has neighbours which are to near
  !these neigbours have to be eliminated
  !index is the point, which is to near
  !optim is the position of index in the neighbourlist of point i 
  function search_min(object,i,adapt_par,index,optim,ilev)
    logical                     :: search_min
    type(ap_struct),intent(inout)  :: object
    integer,intent(in)          :: i
    type(adapt_set), intent(in) :: adapt_par
    integer,intent(inout)       :: index,optim
    integer, intent(in)         :: ilev
    
    integer                     :: j
    real(prec)                  :: r, delta_lev

    search_min=.false.

    if (optim<1 .or. optim >object%n) optim=n !unsinnige eingabe von optim 
                                              !=> kein punkt wird betrachtet 

    j=optim !ab dieser Stelle in der Nachbarliste von Punkt i

    if (i<=object%konstant_size) then
       !fuer alle nachbarn und solange noch kein punkt zu nahe
       do while ((j<=object%konstant(i)%nb) .and. (.not. search_min))
          if (object%konstant(i)%ind(j) >i) then       !punkt wurde schon vorher bearbeitet
               r=horiz_distance(object,i,object%konstant(i)%ind(j))
               delta_lev=lev_distance(object,i,object%konstant(i)%ind(j))
               if ((r < f_min(object,i,object%konstant(i)%ind(j),adapt_par,ilev)) .and. &
                   (delta_lev < adapt_par%delta_lev)                   ) then
                   search_min=.true.                    !nachbar zu nahe
                   index=object%konstant(i)%ind(j)
                   optim=j                              !ab dieser Stelle muss das naechste Mal
                                                        !betrachtet werden
               endif
          endif
          j=j+1
       end do
    else
       !fuer alle nachbarn und solange noch kein punkt zu nahe
       do while ((j<=object%dynamic(i-object%konstant_size)%nb) .and. (.not. search_min))
          if (object%dynamic(i-object%konstant_size)%ind(j) >i) then  !punkt wurde schon vorher bearbeitet
               r=horiz_distance(object,i,object%dynamic(i-object%konstant_size)%ind(j))
               delta_lev=lev_distance(object,i,object%dynamic(i-object%konstant_size)%ind(j))
               if ((r < f_min(object,i,object%dynamic(i-object%konstant_size)&
                    &%ind(j),adapt_par,ilev)) .and. &
                   (delta_lev < adapt_par%delta_lev)                                       ) then
                   search_min=.true.                  !nachbar zu nahe
                   index=object%dynamic(i-object%konstant_size)%ind(j)
                   optim=j                            !ab dieser Stelle muss das naechste Mal
                                                      !betrachtet werden
               endif
          endif
          j=j+1
       end do
    endif
    
  end function search_min


  !function to check if a point has neighbours which are to far
  !between these points, a new point has to be created
  !index is the point, which is to far
  !optim is the position of index in the neighbourlist of point i 
  function search_max(object,i,adapt_par,index,optim,ilev)
    logical                     :: search_max
    type(ap_struct),intent(inout)  :: object
    integer,intent(in)          :: i
    type(adapt_set), intent(in) :: adapt_par
    integer,intent(inout)       :: index,optim
    integer, intent(in)         :: ilev

    integer                     :: j
    real(prec)                  :: r, delta_lev

    search_max=.false.
 
    if (optim<1 .or. optim >object%n) optim=n

    j=optim !ab dieser Stelle in der Nachbarliste von Punkt i

    if (i<=object%konstant_size) then
       !fuer alle nachbarn und solange noch kein punkt zu weit entfernt
       do while ((j<=object%konstant(i)%nb) .and. (.not. search_max))
          if (object%konstant(i)%ind(j) >i) then     !punkt wurde schon vorher bearbeitet
               r=horiz_distance(object,i,object%konstant(i)%ind(j))
               delta_lev=lev_distance(object,i,object%konstant(i)%ind(j))
               if ((r > f_max(object,i,object%konstant(i)%ind(j),adapt_par,ilev)) .and. &
                   (delta_lev < adapt_par%delta_lev)                     .and. &
                   (r < f_limit(object,i,object%konstant(i)%ind(j),adapt_par,ilev))) then
                  search_max=.true.                  !nachbar zu weit entfernt
                  index=object%konstant(i)%ind(j)
                  optim=j                            !ab dieser Stelle muss das naechste Mal
                                                     !betrachtet werden
               endif
          endif
          j=j+1
       end do
    else
       !fuer alle nachbarn und solange noch kein punkt zu weit entfernt
       do while ((j<=object%dynamic(i-object%konstant_size)%nb) .and. (.not. search_max))
          if (object%dynamic(i-object%konstant_size)%ind(j) >i) then  !punkt wurde schon vorher bearbeitet
               r=horiz_distance(object,i,object%dynamic(i-object%konstant_size)%ind(j))
               delta_lev=lev_distance(object,i,object%dynamic(i-object%konstant_size)%ind(j))
               if ((r > f_max(object,i,object%dynamic(i-object%konstant_size)&
                    &%ind(j),adapt_par,ilev)) .and. &
                   (delta_lev < adapt_par%delta_lev)                                         .and. &
                   (r < f_limit(object,i,object%dynamic(i-object%konstant_size)&
                   &%ind(j),adapt_par,ilev))) then
                  search_max=.true.                 !nachbar zu weit entfernt
                  index=object%dynamic(i-object%konstant_size)%ind(j)
                  optim=j                           !ab dieser Stelle muss das naechste Mal
                                                    !betrachtet werden
               endif
          endif
          j=j+1
       end do
    endif

  end function search_max


  !function to check if a point has neighbours which are to far
  !in vertical coordinate between these points
  !index is the point, which is to far
  !optim is the position of index in the neighbourlist of point i  
  function search_vertical_max(object,i,adapt_par,index,optim)

    logical                        :: search_vertical_max
    type(ap_struct),intent(inout)  :: object
    integer,intent(in)             :: i
    type(adapt_set), intent(in)    :: adapt_par
    integer,intent(inout)          :: index,optim 
    integer                        :: j
    real                           :: delta_lev

    search_vertical_max=.false.
       
    if (optim<1 .or. optim >object%n) optim=n 
    j=optim
    if (i<=object%konstant_size) then
       !fuer alle nachbarn und solange noch kein punkt zu weit entfernt
       do while ((j<=object%konstant(i)%nb) .and. (.not. search_vertical_max))
          if (object%konstant(i)%ind(j) >i) then
             !punkt wurde schon vorher bearbeitet
             delta_lev=lev_distance_new(object,i,object%konstant(i)%ind(j))
             search_vertical_max=.true.
             index=object%konstant(i)%ind(j)
             optim=j
          endif
          j=j+1  
        end do
     else
        do while ((j<=object%dynamic(i-object%konstant_size)%nb) .and. (.not. search_vertical_max))
           if (object%dynamic(i-object%konstant_size)%ind(j) >i) then  
              !punkt wurde schon vorher bearbeitet
              delta_lev=lev_distance_new(object,i,object%dynamic(i-object%konstant_size)%ind(j))
              search_vertical_max=.true.
              index=object%dynamic(i-object%konstant_size)%ind(j)
              optim=j
          endif
          j=j+1
       end do
    endif
           
  end function search_vertical_max


  ! function to check if the bvf or bvfold of point i is smaller than fac_bvf_min
  function search_bvf_min(object,i,adapt_par)

    logical                        :: search_bvf_min
    type(ap_struct),intent(inout)  :: object
    integer,intent(in)             :: i
    type(adapt_set), intent(in)    :: adapt_par
    real                           :: bvf,f_bvf_min
    search_bvf_min=.false.
    f_bvf_min=adapt_par%fac_bvf_min   

     if (i<=object%konstant_size) then

        !bvf =object%konstant(i)%bvf
        bvf =object%konstant(i)%bvf_old
        ! control the missing values, modified 2016.09.08 by M.Tao
        if (bvf < f_bvf_min .and. bvf > -100.)then 
               search_bvf_min=.true.
        endif

     else

        !bvf =object%dynamic(i-object%konstant_size)%bvf
        bvf =object%dynamic(i-object%konstant_size)%bvf_old
        ! control the missing values, modified 2016.09.08 by M.Tao
        if (bvf < f_bvf_min .and. bvf > -100.)then
               search_bvf_min=.true.
        endif
        
     endif
   end function search_bvf_min
             

  !function to get the horiz_distance between two points (cartesian horiz_distance)
  function horiz_distance(object,i,j)
    real(prec)                     :: horiz_distance
    type(ap_struct),intent(in)     :: object
    integer,intent(in)             :: i,j

    real(prec)                     :: norm
    real(prec),dimension(3)        :: coor_i, coor_j


       if (i<=object%konstant_size) then
          if (j<=object%konstant_size) then
             coor_i=object%konstant(i)%coor
             coor_j=object%konstant(j)%coor
          else
             coor_i=object%konstant(i)%coor
             coor_j=object%dynamic(j-object%konstant_size)%coor
          endif
       else
          ! j is always larger than i, so j is in dynamic, too
          coor_i=object%dynamic(i-object%konstant_size)%coor
          coor_j=object%dynamic(j-object%konstant_size)%coor
       endif
       norm=sqrt(coor_i(1)**2+coor_i(2)**2+coor_i(3)**2)
       coor_i=coor_i/norm
       norm=sqrt(coor_j(1)**2+coor_j(2)**2+coor_j(3)**2)
       coor_j=coor_j/norm

       horiz_distance = &
            sqrt((coor_i(1)-coor_j(1))**2+(coor_i(2)-coor_j(2))**2+(coor_i(3)-coor_j(3))**2)

  end function horiz_distance

  !function to get the vertical distance between two points
  ! (abs value of the difference between the theta/zeta levels)
  function lev_distance(object,i,j)
    real(prec)                  :: lev_distance
    type(ap_struct),intent(in)  :: object
    integer,intent(in)          :: i,j

       if (i<=object%konstant_size) then
          if (j<=object%konstant_size) then
             lev_distance = &
                  abs(object%konstant(i)%lev-object%konstant(j)%lev)
          else
             lev_distance = &
                  abs(object%konstant(i)%lev &
                  -object%dynamic(j-object%konstant_size)%lev)
          endif
       else
          ! j is always larger than i, so j is in dynamic, too
             lev_distance = &
                  abs(object%dynamic(i-object%konstant_size)%lev &
                  -object%dynamic(j-object%konstant_size)%lev)
       endif

  end function lev_distance

  !function to calculate the vertical distance between
  !actual vertical coordinate of point i and old vertical coordinate of point j
  function lev_distance_new(object,i,j)

    real                        :: lev_distance_new
    type(ap_struct),intent(in)  :: object
    integer,intent(in)          :: i,j
    
    if (i<=object%konstant_size) then
       if (j<=object%konstant_size) then
          lev_distance_new = &
               abs(object%konstant(i)%theta-object%konstant(j)%theta_old)
       else
          lev_distance_new = abs(object%konstant(i)%theta&
               -object%dynamic(j-object%konstant_size)%theta_old)
       endif
    else 
       ! j is always larger than i, so j is in dynamic, too
       lev_distance_new = &
            abs(object%dynamic(i-object%konstant_size)%theta&
            -object%dynamic(i-object%konstant_size)%theta_old)
    endif

  end function lev_distance_new


  logical function lev_gt_min (object,index1,index2,l_min_act)

    type(ap_struct),intent(inout)  :: object
    integer,        intent(in)     :: index1, index2
    real(prec),     intent(in)     :: l_min_act

    lev_gt_min = .false.

    if (index1<=object%konstant_size) then
       if (index2<=object%konstant_size) then
          if (object%konstant(index1)%lev>=l_min_act .and. &
               object%konstant(index2)%lev>=l_min_act) then
             lev_gt_min = .true.
          endif
       else
          if (object%konstant(index1)%lev>=l_min_act .and. &
               object%dynamic(index2-object%konstant_size)%lev>=l_min_act) then
             lev_gt_min = .true.
          endif
       endif
    else
       ! index2 is always larger than index1, so index2 is in dynamic, too
       if (object%dynamic(index1-object%konstant_size)%lev>=l_min_act .and. &
            object%dynamic(index2-object%konstant_size)%lev>=l_min_act) then
          lev_gt_min = .true.
       endif
    endif

  end function lev_gt_min


  !subroutine to eliminate the point 'index', because the distance to point 'i' is to small
  !tracer is interpolated between the two points
  subroutine eliminate(object,i,index,a_p)
    type(ap_struct),intent(inout)  :: object
    integer,intent(in)             :: i,index
    type(adapt_set), intent(in)    :: a_p

    type(ap), pointer              :: new    
    logical                        :: chem_state

    chem_state=associated(object%konstant(1)%c)
    !interpolate tracer from the two points

    if (i<=object%konstant_size) then
       if (index<=object%konstant_size) then
          !mark the point changed by elimination
          if (object%konstant(i)%state==2 .or. object%konstant(index)%state==2 .or. &
               object%konstant(i)%state==3 .or. object%konstant(index)%state==3) then
             object%konstant(i)%state=3 ! elimination induced by adaptation
          else
             object%konstant(i)%state=1 ! pure elimination
          endif
          object%konstant(i)%lev = &
               (object%konstant(i)%lev+object%konstant(index)%lev)/2
          object%konstant(i)%coor = &
               (object%konstant(i)%coor+object%konstant(index)%coor)/2
          object%konstant(i)%coor_old = &
               (object%konstant(i)%coor_old+object%konstant(index)%coor_old)/2

          if (object%konstant(i)%subset .or. object%konstant(index)%subset) then
             object%konstant(i)%subset = .true.
          else
             object%konstant(i)%subset = .false.
          endif

       else
          object%konstant(i)%state=3 !mark the point changed by elimination induced by adaptation
          object%konstant(i)%lev = &
               (object%konstant(i)%lev+object%dynamic(index-object%konstant_size)%lev)/2
          object%konstant(i)%coor = &
               (object%konstant(i)%coor+object%dynamic(index-object%konstant_size)%coor)/2
          object%konstant(i)%coor_old = &
               (object%konstant(i)%coor_old+ &
               object%dynamic(index-object%konstant_size)%coor_old)/2

          if (object%konstant(i)%subset .or. object%dynamic(index-object%konstant_size)%subset) then
             object%konstant(i)%subset = .true.
          else
             object%konstant(i)%subset = .false.
          endif

       endif
       
       new=>object%konstant(i)
       
     else
       ! index is always larger than i, so index is in dynamic, too
       object%dynamic(i-object%konstant_size)%state=3    !mark the point changed by elimination induced by adaptation
       object%dynamic(i-object%konstant_size)%lev = &
            (object%dynamic(i-object%konstant_size)%lev  &
            +object%dynamic(index-object%konstant_size)%lev)/2       

       object%dynamic(i-object%konstant_size)%coor = &
            (object%dynamic(i-object%konstant_size)%coor+ &
            object%dynamic(index-object%konstant_size)%coor)/2
       object%dynamic(i-object%konstant_size)%coor_old = &
            (object%dynamic(i-object%konstant_size)%coor_old+ &
            object%dynamic(index-object%konstant_size)%coor_old)/2

       if (object%dynamic(i-object%konstant_size)%subset .or. &
            object%dynamic(index-object%konstant_size)%subset) then
          object%dynamic(i-object%konstant_size)%subset = .true.
       else
          object%dynamic(i-object%konstant_size)%subset = .false.
       endif
       
       new=>object%dynamic(i-object%konstant_size)

     endif

     call put_on_level(object,a_p,new)

     call lose_concentration(object,i,index)    

     !delete object(index)
     if (index<=object%konstant_size) then
        deallocate(object%konstant(index)%ind) !give the neighbourlist free
        if (associated(object%konstant(index)%c)) deallocate(object%konstant(index)%c)
        if (index<object%n) then    !else : last point (nothing has to be done) 
           !copy last point to the point 'index'
           if (object%n<=object%konstant_size) then             
              object%konstant(index)=object%konstant(object%n) !hoffentlich wird pointer richtig kopiert
              nullify(object%konstant(object%n)%ind)           !this point doesn't exist any more             
              nullify(object%konstant(object%n)%c) 
           else
              object%konstant(index)=object%dynamic(object%n-object%konstant_size)
              nullify(object%dynamic(object%n-object%konstant_size)%ind)
              nullify(object%dynamic(object%n-object%konstant_size)%c)
           endif
        endif
     else
        deallocate(object%dynamic(index-object%konstant_size)%ind)
        if (associated(object%dynamic(index-object%konstant_size)%c)) &
             deallocate(object%dynamic(index-object%konstant_size)%c)
        if (index<object%n) then    !else : last point (nothing has to be done) 
           !copy last point to the point 'index'
           object%dynamic(index-object%konstant_size) = &
                object%dynamic(object%n-object%konstant_size)
           nullify(object%dynamic(object%n-object%konstant_size)%ind)
           nullify(object%dynamic(object%n-object%konstant_size)%c)
        endif
     endif

     object%n=object%n-1  !one point is deleted
     n=n-1
     
  end subroutine eliminate

  subroutine lose_concentration(object,i,index)
    type(ap_struct),intent(inout)  :: object
    integer,intent(in)             :: i,index

    logical                        :: chem_state 
   
    chem_state=associated(object%konstant(1)%c)

    if (i<=object%konstant_size) then
       if (index<=object%konstant_size) then
          object%konstant(i)%tracer = &
               (object%konstant(i)%tracer+object%konstant(index)%tracer)/2
          if (chem_state) &
               object%konstant(i)%c=(object%konstant(i)%c+object%konstant(index)%c)/2
       else
          object%konstant(i)%tracer = &
               (object%konstant(i)%tracer+ &
               object%dynamic(index-object%konstant_size)%tracer)/2
          if (chem_state) &
               object%konstant(i)%c = &
               (object%konstant(i)%c+object%dynamic(index-object%konstant_size)%c)/2
       end if
    else
       object%dynamic(i-object%konstant_size)%tracer = &
            (object%dynamic(i-object%konstant_size)%tracer+ &
            object%dynamic(index-object%konstant_size)%tracer)/2
       if (chem_state) &
            object%dynamic(i-object%konstant_size)%c = &
            (object%dynamic(i-object%konstant_size)%c+ &
            object%dynamic(index-object%konstant_size)%c)/2             
    end if
  end subroutine lose_concentration

  !subroutine to create a new point between two points, which distance is too far
  subroutine new_point (status,object,i,a_p,index)
    integer                        :: status
    type(ap_struct),intent(inout)  :: object
    integer,intent(in)             :: i,index
    type(adapt_set), intent(in)    :: a_p

    type(ap), pointer              :: old1,old2,new    

    status = 0 ! no error

    !lege neuen punkt an
    if (object%n >= object%konstant_size) then
       !look, if new memory is needed
       call memory_management(object)
    endif

    object%n=object%n+1
    n=n+1

    !ermittle koordinaten
    !ermittle lat,lon,level,time_init
    !ermittle nachbarn
    if (object%n<=object%konstant_size) then
       !new point('n') and the old points are in the konstant array

       !Mark the new point by index 2
       object%konstant(object%n)%state=2
       object%konstant(object%n)%state_vert=0

       ! Define and normalize the coordinates of new point
       object%konstant(object%n)%coor = &
            (object%konstant(i)%coor+object%konstant(index)%coor)/2
       object%konstant(object%n)%coor_old = &
            (object%konstant(i)%coor_old+object%konstant(index)%coor_old)/2

       object%konstant(object%n)%lev = &
            (object%konstant(i)%lev+object%konstant(index)%lev)/2

       ! Define remaining componenents 
       object%konstant(object%n)%time_init = object%konstant(i)%time_init
       object%konstant(object%n)%theta = (object%konstant(i)%theta+object%konstant(index)%theta)/2
       object%konstant(object%n)%theta_old = (object%konstant(i)%theta_old+object%konstant(index)%theta_old)/2
       object%konstant(object%n)%bvf = (object%konstant(i)%bvf+object%konstant(index)%bvf)/2
       object%konstant(object%n)%bvf_old = (object%konstant(i)%bvf_old+object%konstant(index)%bvf_old)/2

       if (object%konstant(i)%subset .or. object%konstant(index)%subset) then
          object%konstant(object%n)%subset = .true.
       else
           object%konstant(object%n)%subset = .false.
       endif

       old1=>object%konstant(i)
       old2=>object%konstant(index)
       new=>object%konstant(object%n)
       !get all neighbours of the new point
       call get_neighbours (status,object,old1,i,old2,index,new)
       if (status/=0) then
          write (*,*) 'Error in get_neighbours !!!'
          return
       endif
                     
    else
       !new point in dynamic array
       if (i<=object%konstant_size) then
          if (index<=object%konstant_size) then
             !i and index in konstant array
             object%dynamic(object%n-object%konstant_size)%coor = &
                  (object%konstant(i)%coor+object%konstant(index)%coor)/2
             object%dynamic(object%n-object%konstant_size)%coor_old = &
                  (object%konstant(i)%coor_old+object%konstant(index)%coor_old)/2
             object%dynamic(object%n-object%konstant_size)%lev = &
                  (object%konstant(i)%lev+object%konstant(index)%lev)/2
             object%dynamic(object%n-object%konstant_size)%time_init = &
                  object%konstant(i)%time_init
             object%dynamic(object%n-object%konstant_size)%theta = &
                  (object%konstant(i)%theta+object%konstant(index)%theta)/2
             object%dynamic(object%n-object%konstant_size)%theta_old = &
                  (object%konstant(i)%theta_old+object%konstant(index)%theta_old)/2
             object%dynamic(object%n-object%konstant_size)%bvf = &
                 (object%konstant(i)%bvf+object%konstant(index)%bvf)/2
             object%dynamic(object%n-object%konstant_size)%bvf_old = &
                 (object%konstant(i)%bvf_old+object%konstant(index)%bvf_old)/2

             if (object%konstant(i)%subset .or. object%konstant(index)%subset) then
                object%dynamic(object%n-object%konstant_size)%subset = .true.
             else
                object%dynamic(object%n-object%konstant_size)%subset = .false.
             endif

             old1=>object%konstant(i)
             old2=>object%konstant(index)
             new=>object%dynamic(object%n-object%konstant_size)
             call get_neighbours (status,object,old1,i,old2,index,new)
             if (status/=0) then
                write (*,*) 'Error in get_neighbours !!!'
                return
             endif
          else
             !i in konstant, index in dynamic array
             object%dynamic(object%n-object%konstant_size)%coor = &
                  (object%konstant(i)%coor+ &
                  object%dynamic(index-object%konstant_size)%coor)/2
             object%dynamic(object%n-object%konstant_size)%coor_old = &
                  (object%konstant(i)%coor_old+ &
                  object%dynamic(index-object%konstant_size)%coor_old)/2
             object%dynamic(object%n-object%konstant_size)%lev = &
                  (object%konstant(i)%lev &
                  +object%dynamic(index-object%konstant_size)%lev)/2
             object%dynamic(object%n-object%konstant_size)%time_init = &
                  object%konstant(i)%time_init
             object%dynamic(object%n-object%konstant_size)%theta = &
                  (object%konstant(i)%theta &
                  +object%dynamic(index-object%konstant_size)%theta)/2
             object%dynamic(object%n-object%konstant_size)%theta_old = &
                  (object%konstant(i)%theta_old &
                  +object%dynamic(index-object%konstant_size)%theta_old)/2
             object%dynamic(object%n-object%konstant_size)%bvf = &
                  (object%konstant(i)%bvf & 
                  +object%dynamic(index-object%konstant_size)%bvf)/2
             object%dynamic(object%n-object%konstant_size)%bvf_old = &
                  (object%konstant(i)%bvf_old & 
                  +object%dynamic(index-object%konstant_size)%bvf_old)/2

             if (object%konstant(i)%subset .or. &
                  object%dynamic(index-object%konstant_size)%subset) then
                object%dynamic(object%n-object%konstant_size)%subset = .true.
             else
                object%dynamic(object%n-object%konstant_size)%subset = .false.
             endif

             old1=>object%konstant(i)
             old2=>object%dynamic(index-object%konstant_size)
             new=>object%dynamic(object%n-object%konstant_size)
             call get_neighbours (status,object,old1,i,old2,index,new)
             if (status/=0) then
                write (*,*) 'Error in get_neighbours !!!'
                return
             endif
          endif
       else
          ! i in dynamic array
          ! index is always larger than i, so index is in dynamic, too          
          object%dynamic(object%n-object%konstant_size)%coor = &
               (object%dynamic(i-object%konstant_size)%coor &
               +object%dynamic(index-object%konstant_size)%coor)/2
          object%dynamic(object%n-object%konstant_size)%coor_old = &
               (object%dynamic(i-object%konstant_size)%coor_old &
               +object%dynamic(index-object%konstant_size)%coor_old)/2
          object%dynamic(object%n-object%konstant_size)%lev = &
               (object%dynamic(i-object%konstant_size)%lev &
               +object%dynamic(index-object%konstant_size)%lev)/2          
          object%dynamic(object%n-object%konstant_size)%time_init = &
               object%dynamic(i-object%konstant_size)%time_init
          object%dynamic(object%n-object%konstant_size)%theta = &
               (object%dynamic(i-object%konstant_size)%theta &
               +object%dynamic(index-object%konstant_size)%theta)/2
          object%dynamic(object%n-object%konstant_size)%theta_old = &
               (object%dynamic(i-object%konstant_size)%theta_old &
               +object%dynamic(index-object%konstant_size)%theta_old)/2
          object%dynamic(object%n-object%konstant_size)%bvf = &
               (object%dynamic(i-object%konstant_size)%bvf & 
               +object%dynamic(index-object%konstant_size)%bvf)/2
          object%dynamic(object%n-object%konstant_size)%bvf_old = &
               (object%dynamic(i-object%konstant_size)%bvf_old & 
               +object%dynamic(index-object%konstant_size)%bvf_old)/2

          if (object%dynamic(i-object%konstant_size)%subset .or. &
               object%dynamic(index-object%konstant_size)%subset) then
             object%dynamic(object%n-object%konstant_size)%subset = .true.
          else
             object%dynamic(object%n-object%konstant_size)%subset = .false.
          endif

          old1=>object%dynamic(i-object%konstant_size)
          old2=>object%dynamic(index-object%konstant_size)
          new=>object%dynamic(object%n-object%konstant_size)
          call get_neighbours (status,object,old1,i,old2,index,new)
          if (status/=0) then
             write (*,*) 'Error in get_neighbours !!!'
             return
          endif
       endif
       !Mark the new point by index 2
       object%dynamic(object%n-object%konstant_size)%state=2
       object%dynamic(object%n-object%konstant_size)%state_vert=0
       
    endif

    !normalize point or put point to his theta/zeta-level
    call put_on_level(object,a_p,new)

  end subroutine new_point

  subroutine put_on_level(object,a_p,new)
    type(ap_struct),intent(inout)  :: object
    type(adapt_set), intent(in)    :: a_p
    type(ap), pointer              :: new    
    
    real(prec)                     :: pi,r
    real(prec)                     :: length
    real(prec),dimension(3)        :: coor_h, coor_h_old
    real(prec),dimension(2)        :: lon_lat

    pi=4.*atan(1.) 
       
    ! Define and normalize the coordinates of new point
    coor_h=new%coor       
    coor_h_old=new%coor_old
    ! Normalize for 3d ghull
    length=sqrt(coor_h(1)**2+coor_h(2)**2+coor_h(3)**2) 
    new%coor=coor_h/length
       
    length=sqrt(coor_h_old(1)**2+coor_h_old(2)**2+coor_h_old(3)**2) 
    new%coor_old=coor_h_old/length
    
    lon_lat=kart2ukmo(new%coor)
    new%lon=lon_lat(1)
    new%lat=lon_lat(2)
    
!!!!!
!     if (object%dimension ==4) then 
          
!        !put point on his theta-level
!        r=1.+((new%theta-a_p%theta_min)/(a_p%theta_max-a_p%theta_min))* &
!             a_p%nthetas*(sqrt(a_p%area_h)/a_p%aspect_ratio)
!        new%coor(1) = r*cos(new%lon*(pi/180.))* sin(abs(new%lat*(pi/180.)-pi/2.)) 
!        new%coor(2) = r*sin(new%lon*(pi/180.))* sin(abs(new%lat*(pi/180.)-pi/2.))
!        new%coor(3) = r*cos(abs(new%lat*(pi/180.)-pi/2.))
       
!        lon_lat=kart2ukmo(new%coor_old)
!        new%coor_old(1) = r*cos(lon_lat(1)*(pi/180.))* sin(abs(lon_lat(2)*(pi/180.)-pi/2.)) 
!        new%coor_old(2) = r*sin(lon_lat(1)*(pi/180.))* sin(abs(lon_lat(2)*(pi/180.)-pi/2.))
!        new%coor_old(3) = r*cos(abs(lon_lat(2)*(pi/180.)-pi/2.))
       
!     end if
    
  end subroutine put_on_level


  !subroutine to interpolate tracer concentration
  !this routine uses only 2 points for interpolation 
  subroutine interpolate(object)
    type(ap_struct),intent(inout)  :: object
    integer                        :: j,size_c
    logical                        :: chem_state

    chem_state=associated(object%konstant(1)%c)
    !look, if chemie information is available
    if (chem_state) then
       size_c=size(object%konstant(1)%c)
    else
       size_c=0
    end if
    
    if (object%n<=object%konstant_size) then
       !last point(=new point) is in konstant array => all points are in konstant array 
       object%konstant(object%n)%tracer=0.
       if (chem_state) then 
          allocate(object%konstant(object%n)%c(size_c))
          object%konstant(object%n)%c=0.
       end if
       do j=1, 2
          object%konstant(object%n)%tracer=object%konstant(object%n)%tracer &
               +object%konstant(object%konstant(object%n)%ind(j))%tracer
       if (chem_state) &
            object%konstant(object%n)%c=object%konstant(object%n)%c &
            +object%konstant(object%konstant(object%n)%ind(j))%c
       end do
       object%konstant(object%n)%tracer=object%konstant(object%n)%tracer/2.
       if (chem_state) &
            object%konstant(object%n)%c=object%konstant(object%n)%c/2.
    else
       object%dynamic(object%n-object%konstant_size)%tracer=0.
       if (chem_state) then 
          allocate(object%dynamic(object%n-object%konstant_size)%c(size_c))
          object%dynamic(object%n-object%konstant_size)%c=0.
       end if
       do j=1, 2
          !for all neighbours
          if (object%dynamic(object%n-object%konstant_size)%ind(j)<=object%konstant_size) then
             !neighbour is in konstant array
             object%dynamic(object%n-object%konstant_size)%tracer = &
                  object%dynamic(object%n-object%konstant_size)%tracer  &
                  +object%konstant(object%dynamic(object%n-object%konstant_size)%ind(j))%tracer
             if (chem_state) &
                  object%dynamic(object%n-object%konstant_size)%c = &
                  object%dynamic(object%n-object%konstant_size)%c  &
                  +object%konstant(object%dynamic(object%n-object%konstant_size)%ind(j))%c
          else 
             !neighbour is in dynamic array
             object%dynamic(object%n-object%konstant_size)%tracer = &
                  object%dynamic(object%n-object%konstant_size)%tracer  &
                  +object%dynamic(object%dynamic(object%n-object%konstant_size)%ind(j) &
                  -object%konstant_size)%tracer
             if (chem_state) &
                  object%dynamic(object%n-object%konstant_size)%c = &
                  object%dynamic(object%n-object%konstant_size)%c  &
                  +object%dynamic(object%dynamic(object%n-object%konstant_size)%ind(j) &
                  -object%konstant_size)%c
          endif
       end do
       object%dynamic(object%n-object%konstant_size)%tracer=&
            object%dynamic(object%n-object%konstant_size)%tracer/2.
       if (chem_state) &
            object%dynamic(object%n-object%konstant_size)%c=&
            object%dynamic(object%n-object%konstant_size)%c/2.
    endif

  end subroutine interpolate

 ! average the concentration of all species and change the state_vert to 1 of point i and
 !  its neighbour index    
  subroutine avg_concentration(object,i,index)
    type(ap_struct),intent(inout)  :: object
    integer,intent(in)             :: i,index
    logical                        :: chem_state 
   
    chem_state=associated(object%konstant(1)%c)


    if (i<=object%konstant_size) then
        object%konstant(i)%state_vert=1
       if (index<=object%konstant_size) then   
          object%konstant(index)%state_vert=1 
          object%konstant(i)%tracer = &
               (object%konstant(i)%tracer+object%konstant(index)%tracer)/2
          if (chem_state) &
               object%konstant(i)%c=(object%konstant(i)%c+object%konstant(index)%c)/2
       else 
          object%dynamic(index-object%konstant_size)%state_vert=1 
          object%konstant(i)%tracer = &
               (object%konstant(i)%tracer+ &
               object%dynamic(index-object%konstant_size)%tracer)/2
          if (chem_state) &
               object%konstant(i)%c = &
               (object%konstant(i)%c+object%dynamic(index-object%konstant_size)%c)/2
       end if
    else
    ! j is always larger than i, so j is also in dynamic
       object%dynamic(i-object%konstant_size)%state_vert=1
       object%dynamic(index-object%konstant_size)%state_vert=1
       object%dynamic(i-object%konstant_size)%tracer = &
            (object%dynamic(i-object%konstant_size)%tracer+ &
            object%dynamic(index-object%konstant_size)%tracer)/2
       if (chem_state) &
            object%dynamic(i-object%konstant_size)%c = &
            (object%dynamic(i-object%konstant_size)%c+ &
            object%dynamic(index-object%konstant_size)%c)/2             
    end if
  end subroutine avg_concentration

 ! average the concentration of all species of point i and
 !  its all neighbours; change their state_vert to 2  
  subroutine avg_concentration_all(object,i)
    type(ap_struct),intent(inout)  :: object
    integer,intent(in)             :: i
    integer                        :: j,index
    logical                        :: chem_state 
    real,dimension(:),allocatable  :: mean_c
    real                           :: mean_t
    chem_state=associated(object%konstant(1)%c)
    allocate(mean_c(size(object%konstant(1)%c)))
    mean_c=0.
    if (i<=object%konstant_size) then
       j=1
       mean_t=object%konstant(i)%tracer
       if (chem_state) mean_c=object%konstant(i)%c
       do while (j<=object%konstant(i)%nb)
         index=object%konstant(i)%ind(j)
         if (index<=object%konstant_size) then
             mean_t =mean_t+(object%konstant(index)%tracer-mean_t)/(j+1)
             if (chem_state) &
             mean_c=mean_c+(object%konstant(index)%c-mean_c)/(j+1)
         else
             mean_t =mean_t+(object%dynamic(index-object%konstant_size)%tracer-mean_t)/(j+1)
             if (chem_state) &
             mean_c=mean_c+(object%dynamic(index-object%konstant_size)%c-mean_c)/(j+1)
         end if 
         j=j+1
       end do
    else
       j=1
       mean_t=object%dynamic(i-object%konstant_size)%tracer
       if (chem_state) mean_c=object%dynamic(i-object%konstant_size)%c
       do while (j<=object%dynamic(i-object%konstant_size)%nb)
         index=object%dynamic(i-object%konstant_size)%ind(j)
         !index could be larger or smaller than i
         if (index<=object%konstant_size)then
             mean_t =mean_t+(object%konstant(index)%tracer-mean_t)/(j+1)
             if (chem_state) &
             mean_c=mean_c+(object%konstant(index)%c-mean_c)/(j+1)
         else      
             mean_t=mean_t+(object%dynamic(index-object%konstant_size)%tracer-mean_t)/(j+1) 
             if (chem_state) &
             mean_c=mean_c+(object%dynamic(index-object%konstant_size)%c-mean_c)/(j+1)
         end if
         j=j+1
       end do         
    end if
    
 
  ! change concentration of point i and all its neighbours
  if (i<=object%konstant_size) then
     object%konstant(i)%state_vert =2
     object%konstant(i)%tracer=mean_t 
     if (chem_state) object%konstant(i)%c=mean_c
     j=1
     do while (j<=object%konstant(i)%nb)
     index=object%konstant(i)%ind(j)
       if (index<=object%konstant_size) then
         object%konstant(index)%state_vert =2
         object%konstant(index)%tracer=mean_t 
         if (chem_state) object%konstant(index)%c=mean_c
       else
         object%dynamic(index-object%konstant_size)%state_vert =2
         object%dynamic(index-object%konstant_size)%tracer=mean_t 
         if (chem_state) object%dynamic(index-object%konstant_size)%c=mean_c
       end if
       j=j+1
     end do
  else
     object%dynamic(i-object%konstant_size)%state_vert =2
     object%dynamic(i-object%konstant_size)%tracer=mean_t 
     if (chem_state) object%dynamic(i-object%konstant_size)%c=mean_c
     j=1
     do while (j<=object%dynamic(i-object%konstant_size)%nb)
       index=object%dynamic(i-object%konstant_size)%ind(j)
       if (index<=object%konstant_size) then
         object%konstant(index)%state_vert =2
         object%konstant(index)%tracer=mean_t 
         if (chem_state) object%konstant(index)%c=mean_c
       else
         object%dynamic(index-object%konstant_size)%state_vert =2
         object%dynamic(index-object%konstant_size)%tracer=mean_t 
         if (chem_state) object%dynamic(index-object%konstant_size)%c=mean_c
       end if 
       j=j+1
     end do
  end if
  deallocate(mean_c)
  end subroutine avg_concentration_all


  !the subroutine is called, when a new point is added
  !it looks if new memory is needed to save the new point
  !or if there is some place left in the dynamic array
  !if new memory is needed, the old dynamic array is copied
  !into a larger array
  subroutine memory_management(object)
    type(ap_struct),intent(inout)            :: object
    type(ap),dimension(:), pointer           :: help_array

    integer                                  :: i

    if (object%n==object%konstant_size) then
       !the konstant array is full => a dynamic array is used to store the next points
       nullify(object%dynamic)
       !num_new_points is the global konstant, which defines the size of the initial 
       !dynamic array
       allocate(object%dynamic(num_new_points))

       do i=1,num_new_points
          nullify(object%dynamic(i)%c)          
          nullify(object%dynamic(i)%ind) !all neighbourlist-pointers are undefined
       end do
       object%max_points=object%max_points+num_new_points !now there is more memory       
    else
       if (object%n==object%max_points) then
          !the dynamic array is full => the dynamic array is resized
          nullify(help_array)
          allocate(help_array(object%max_points-object%konstant_size)) !the size of the dynamic array 
          !copy the whole dynamic array
          help_array=object%dynamic 

          do i=1,object%max_points-object%konstant_size
             if (associated(object%dynamic(i)%ind)) &
                nullify(object%dynamic(i)%ind)  !set all neighbour pointer to NULL
                ! (do not give free the memory!) 
             if (associated(object%dynamic(i)%c)) &
                  nullify(object%dynamic(i)%c)  !set all c pointer to NULL            
          end do

          !resize the dynamic array; num_new_points is the global constant, which defines
          !how many points are added
          deallocate(object%dynamic) 
          nullify(object%dynamic)
          allocate(object%dynamic(object%max_points-object%konstant_size+num_new_points))

          do i=object%max_points-object%konstant_size+1, object%max_points-object%konstant_size+num_new_points
             nullify(object%dynamic(i)%c)
             nullify(object%dynamic(i)%ind)
          enddo
          !copy the whole former dynamic array to the new dynamic array
          object%dynamic(1:object%max_points-object%konstant_size)=help_array

          object%max_points=object%max_points+num_new_points !now there is more memory

          deallocate(help_array)
          nullify(help_array)

       endif
    endif

  end subroutine memory_management

  !subroutine to get all neigbours of a new point
  !this are the two points, between the new point is created, and
  !all same neighbours of these points
  subroutine get_neighbours (status,object,old1,index1,old2,index2,new)
    integer                         :: status
    type(ap_struct),intent(inout)   :: object
    type(ap), pointer               :: old1,old2,new
    integer,intent(in)              :: index1,index2

    integer             :: i,j

    status = 0 ! no error

    nullify(new%ind)
    allocate(new%ind(min(old1%nb,old2%nb)+2)) !maximal number of neigbours

    !the original 2 points
    new%ind(1)=index1
    new%ind(2)=index2
    new%nb=2

    !all other eqal neighbours
    do i=1,old1%nb
       do j=1,old2%nb
          if (old1%ind(i)==old2%ind(j)) then 
             !this point is neighbour of point index1 and point index2 => neigbour of the new point
             new%nb=new%nb+1
             new%ind(new%nb)=old1%ind(i)
          endif
       end do
    end do

    !check 
    if (new%nb > min(old1%nb,old2%nb)+2) then
       print *, 'Err #nb',new%nb
       status = -1
       return
    endif

  end subroutine get_neighbours

  !the soubroutine changes the number 'n+1' to 'index' in all neighbourlists
  !this is, because point 'index' has been eliminated and is replaced by the last point(now 'n+1')
  !only the neighbours of point 'index' have to be changed
  subroutine change_last_point(object,index)
    type(ap_struct),intent(inout)  :: object
    integer,intent(in)             :: index

    integer                        :: i,j,z
    logical                        :: changed

    if (index<=object%konstant_size) then
       do i=1,object%konstant(index)%nb
          !for all neighbours of point 'index'
          changed=.false.

          if (object%konstant(index)%ind(i)<=object%konstant_size) then
             ! this neighbour is in the konstant array
             z=object%konstant(index)%ind(i)

             j=1
             !search in the neighbourlist of this point for point 'index'
             do while ((j<=object%konstant(z)%nb) .and. (.not. changed))
                if (object%konstant(z)%ind(j) == object%n+1) then
                   object%konstant(z)%ind(j)=index  !  n+1 is the old number, index the new number 
                   changed=.true.
                endif
                j=j+1
             end do
          else             
             ! this neighbour is in the dynamic array
             z=object%konstant(index)%ind(i)-object%konstant_size

             j=1
             !search in the neighbourlist of this point for point 'index'
             do while ((j<=object%dynamic(z)%nb) .and. (.not. changed))
                if (object%dynamic(z)%ind(j) == object%n+1) then
                   object%dynamic(z)%ind(j)=index  !  n+1 is the old number, index the new number
                   changed=.true.
                endif
                j=j+1
             end do
          endif
       end do
    else
       !point 'index' is in dynamic array
       do i=1,object%dynamic(index-object%konstant_size)%nb
          !for all neighbours of point 'index'
          changed=.false.

          if (object%dynamic(index-object%konstant_size)%ind(i)<=object%konstant_size) then
             ! this neighbour is in the konstant array
             z=object%dynamic(index-object%konstant_size)%ind(i)

             j=1
             !search in the neighbourlist of this point for point 'index'
             do while ((j<=object%konstant(z)%nb) .and. (.not. changed))
                if (object%konstant(z)%ind(j) == object%n+1) then
                   object%konstant(z)%ind(j)=index !  n+1 is the old number, index the new number 
                   changed=.true.
                endif
                j=j+1
             end do
          else             
             ! this neighbour is in the dynamic array
             z=object%dynamic(index-object%konstant_size)%ind(i)-object%konstant_size

             j=1
             !search in the neighbourlist of this point for point 'index'
             do while ((j<=object%dynamic(z)%nb) .and. (.not. changed))
                if (object%dynamic(z)%ind(j) == object%n+1) then
                   object%dynamic(z)%ind(j)=index   !  n+1 is the old number, index the new number 
                   changed=.true.
                endif
                j=j+1
             end do
          endif
       end do
    endif

  end subroutine change_last_point

  !put the new point in all neighbourlists
  !only the neighbours of the new point have to be changed
  subroutine change_neighbours_of_new_point(object,optim)
    type(ap_struct),intent(inout)  :: object
    integer,intent(in)             :: optim

    integer                        :: j,z,z1
    type(ap), pointer              :: point


    !for the 2 points which distance is too large, the new point replace the old neighbour

    !point 'z1' is the first point
    if (object%n<=object%konstant_size) then
       z1=object%konstant(object%n)%ind(1)
    else
       z1=object%dynamic(object%n-object%konstant_size)%ind(1)
    endif

    !optim is the position of the second point in the neighbourlist of point 2
    if (z1<=object%konstant_size) then
       object%konstant(z1)%ind(optim)=object%n !change it; new point is new neighbour
    else
       object%dynamic(z1-object%konstant_size)%ind(optim)=object%n !change it; new point is new neighbour
    endif

    !point 'z' is the second point
    if (object%n<=object%konstant_size) then
       z=object%konstant(object%n)%ind(2)
    else
       z=object%dynamic(object%n-object%konstant_size)%ind(2)
    endif

    if (z<=object%konstant_size) then
       do j=1,object%konstant(z)%nb !search for point 1 in neighbourlist of point 2
          if (object%konstant(z)%ind(j)==z1) then 
             object%konstant(z)%ind(j)=object%n !change it; new point is new neighbour
             exit !leave loop
          endif
       end do
    else
       do j=1,object%dynamic(z-object%konstant_size)%nb !search for point 1 in neighbourlist of point 2
          if (object%dynamic(z-object%konstant_size)%ind(j)==z1) then
             object%dynamic(z-object%konstant_size)%ind(j)=object%n !change it; new point is new neighbour        
             exit !leave loop
          endif
       end do
    endif

    !all other neighbours

    if (object%n<=object%konstant_size) then
       do j=3,object%konstant(object%n)%nb 
          !for all other neighbours of the new point
          z=object%konstant(object%n)%ind(j)
          if (object%konstant(z)%nb == size(object%konstant(z)%ind)) then !neighbourarray is to small
             point=>object%konstant(z)
             call index_management(object,point)                 !resize neighbourarray
          endif
          object%konstant(z)%nb=object%konstant(z)%nb+1
          object%konstant(z)%ind(object%konstant(z)%nb)=object%n !put new point in neighbourlist         
       end do
    else
       !new point is in dynamic
       do j=3,object%dynamic(object%n-object%konstant_size)%nb
          !for all other neighbours of the new point
          z=object%dynamic(object%n-object%konstant_size)%ind(j)

          if (z<=object%konstant_size) then
             !neighbour is in konstant array
             if (object%konstant(z)%nb == size(object%konstant(z)%ind)) then !neighbourarray is to small
                point=>object%konstant(z)
                call index_management(object,point)               !resize neighbourarray          
             endif
             object%konstant(z)%nb =object%konstant(z)%nb+1
             object%konstant(z)%ind(object%konstant(z)%nb)=object%n !put new point in neighbourlist 
          else 
             !neighbour is in dynamic array
             z=z-object%konstant_size
             if (object%dynamic(z)%nb == size(object%dynamic(z)%ind)) then !neighbourarray is to small
                point=>object%dynamic(z)
                call index_management(object,point)               !resize neighbourarray  
             endif
             object%dynamic(z)%nb =object%dynamic(z)%nb+1
             object%dynamic(z)%ind(object%dynamic(z)%nb)=object%n !put new point in neighbourlist
          endif
       end do
    endif

  end subroutine change_neighbours_of_new_point

  !the subroutine deletes the point 'index' from all neighbourlists
  !only the neighbours of the point 'index' have to be changed
  subroutine change_neighbours(object,index)
    type(ap_struct),intent(inout)  :: object
    integer,intent(in)             :: index

    integer                        :: i,j,z
    logical                        :: deleted

    if (index<=object%konstant_size) then
       do i=1,object%konstant(index)%nb
          !for all neighbours of the point 'index'
          deleted=.false.

          if (object%konstant(index)%ind(i)<=object%konstant_size) then
             ! this neighbour is in the konstant array
             z=object%konstant(index)%ind(i)

             j=1
             !search in the neighbourlist of this point for point 'index'
             do while ((j<=object%konstant(z)%nb) .and. (.not. deleted))
                if (object%konstant(z)%ind(j) == index) then
                   !point 'index' is found; can be deleted in neighbourlist, because the point doesn't further exist
                   !copy last neigbour to the place, which will be free
                   object%konstant(z)%ind(j)=object%konstant(z)%ind(object%konstant(z)%nb)
                   object%konstant(z)%nb=object%konstant(z)%nb-1
                   deleted=.true.       
                endif
                j=j+1
             end do
          else             
             !this neighbour is in the dynamic array
             z=object%konstant(index)%ind(i)-object%konstant_size

             j=1
             !search in the neighbourlist of this point for point 'index'
             do while ((j<=object%dynamic(z)%nb) .and. (.not. deleted))
                if (object%dynamic(z)%ind(j) == index) then
                   !point 'index' is found; can be deleted in neighbourlist, because the point doesn't further exist
                   !copy last neigbour to the place, which will be free
                   object%dynamic(z)%ind(j)=object%dynamic(z)%ind(object%dynamic(z)%nb)
                   object%dynamic(z)%nb=object%dynamic(z)%nb-1
                   deleted=.true.
                endif
                j=j+1
             end do
          endif
       end do
    else
       !point 'index' is in the dynamic array
       do i=1,object%dynamic(index-object%konstant_size)%nb
          !for all neighbours of the point 'index'

          deleted=.false.         

          if (object%dynamic(index-object%konstant_size)%ind(i)<=object%konstant_size) then
             ! this neighbour is in the konstant array
             z=object%dynamic(index-object%konstant_size)%ind(i)

             j=1
             !search in the neighbourlist of this point for point 'index'
             do while ((j<=object%konstant(z)%nb) .and. (.not. deleted))
                if (object%konstant(z)%ind(j) == index) then
                   !point 'index' is found; can be deleted in neighbourlist, because the point doesn't further exist
                   !copy last neigbour to the place, which will be free
                   object%konstant(z)%ind(j)=object%konstant(z)%ind(object%konstant(z)%nb)
                   object%konstant(z)%nb=object%konstant(z)%nb-1
                   deleted=.true.
                endif
                j=j+1
             end do
          else             
             !this neighbour is in the dynamic array
             z=object%dynamic(index-object%konstant_size)%ind(i)-object%konstant_size
             
             j=1
             !search in the neighbourlist of this point for point 'index'
             do while ((j<=object%dynamic(z)%nb) .and. (.not. deleted))
                if (object%dynamic(z)%ind(j) == index) then
                   !point 'index' is found; can be deleted in neighbourlist, because the point doesn't further exist
                   !copy last neigbour to the place, which will be free
                   object%dynamic(z)%ind(j)=object%dynamic(z)%ind(object%dynamic(z)%nb)
                   object%dynamic(z)%nb=object%dynamic(z)%nb-1
                   deleted=.true.
                endif
                j=j+1
             end do
          endif
       end do
    endif

  end subroutine change_neighbours

  !subroutine to resize the array which contains the neighbours of one point
  !the old neighbours are copied into the new neighbourarray
  subroutine index_management(object,point)
    type(ap_struct),intent(inout)            :: object
    type(ap)                                :: point

    integer,dimension(:),allocatable,target :: help_array

    allocate(help_array(point%nb))

    !save the neighbours
    help_array(1:point%nb)=point%ind(1:point%nb)

    !resize the neighbourarray
    deallocate(point%ind)
    nullify(point%ind)
    allocate(point%ind(point%nb+num_new_ind))

    !copy the neighbours to the new neighbourarray
    point%ind(1:point%nb) = help_array(1:point%nb)

    point%ind(point%nb+1:point%nb+num_new_ind)=0 !initialisize the rest

    deallocate(help_array)

  end subroutine index_management


  function kart2ukmo(x)

    ! PURPOSE: Transform kartesian coordinates to a unit sphere (via projection)
    !          and then into spherical coordinates (UKMO notation)
    !                   0 < lon < 360,  -90 < lat < 90 

    implicit none
    real(prec), dimension(3), intent(in)  :: x
    real(prec), dimension(2)              :: kart2ukmo

    real(prec)                            :: pi, r, lat, lon
    real(prec), dimension(3)              :: x_unit

    pi=4.*atan(1.)
    r = sqrt(dot_product(x,x))
    x_unit =x/r

    lon = atan2(x_unit(2),x_unit(1))*180./pi     !   -pi < atan(x,y) < pi
    if (lon < 0) lon = 360.+lon
! op_pj_20160606+
!!$    if (lon>=360.) lon = modulo(lon,360.)
    if (lon>=360.) lon = modulo(lon,360._prec)
! op_pj_20160606-
    lat = -1.*(acos(x_unit(3))*180./pi-90.)      !   0 < acos(x) < pi

    kart2ukmo=(/lon,lat/)

  end function kart2ukmo


end module messy_clamsmix_ap_m_modify

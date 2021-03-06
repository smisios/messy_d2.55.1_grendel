!********************************************************************************!
! File clams/mix/source/lib_mix.f90
!
! This file is part of the Chemical Lagrangian Model of the
! Stratosphere (CLaMS). CLaMS is a hierarchy of numerical models and
! preprocessors which simulate transport, chemical reactions and
! mixing processes in the stratosphere.  
!
! Copyright (c) 2006
! Paul Konopka, Nicole Thomas
! Forschungszentrum Juelich GmbH
! Last Modified By: N.Thomas
! Last Modified On: Mon May  4 12:03:06 2020
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
!***********************************************************************!
!
! MODULE lib_mix
! ---------------
!
! This MODULE contains subroutines used by mix (and bmix)
!
!    function lev2levdelta
!    subroutine adapt_grid(ap_s,adapt_par,switch)
!
!----------------------------------------------------------------------

Module messy_clamsmix_lib_mix

contains

!****************************************************************************
!
!****************************************************************************
function lev2levdelta (lev_grid, lev_delta, lev)
  
  use messy_clams_global, only:prec

  implicit none

  real(prec), dimension(:),pointer  :: lev_grid  ! vertical grid (lev values)
  real(prec), dimension(:),pointer  :: lev_delta ! array of the layer-thickness
  real(prec)                        :: lev
  real(prec)                        :: lev2levdelta       ! middle of interval

  integer :: nlevs, ilev
  logical :: found

  nlevs = size(lev_grid)
  lev2levdelta = -1
 
  if (lev<lev_grid(1)) then
     lev2levdelta = lev_delta(1)
  elseif (lev>lev_grid(nlevs)) then
     lev2levdelta = lev_delta(nlevs)
  else
     ilev = 1
     found = .false.
     do while (.not. found .and. ilev<nlevs)
        if (lev_grid(ilev)<=lev .and. lev<=lev_grid(ilev+1)) then
           found = .true.
           lev2levdelta = lev_delta(ilev) + (lev-lev_grid(ilev))* &
                              (lev_delta(ilev+1)-lev_delta(ilev)) / &
                              (lev_grid(ilev+1)-lev_grid(ilev))
        else
           ilev = ilev + 1
        endif
     enddo
  endif
       
end function lev2levdelta



!****************************************************************************
! subroutine to realise an adaptive grid method
! first all neighbours of a point, which are to near, are eliminated
! then new points are created between two points, which distance is to large
! adapt_par is a structur, which contains the criterias
! there can be memory problems, if there will be created to much points (far about 30000 points)
!****************************************************************************

subroutine adapt_grid (status, ap_s,adapt_par,switch,l_min_act,ilev)
  use messy_clams_global,    only: prec, rank
  use messy_clamsmix_global, only: ctrl_out
  use messy_clamsmix_ap_m_access
  use messy_clamsmix_ap_m_modify
  implicit none

  integer                       :: status
  type(ap_struct),intent(inout) :: ap_s
  type(adapt_set), intent(in)   :: adapt_par
  integer                       :: switch
  real(prec)                    :: l_min_act
  integer, intent(in)           :: ilev

  integer                       :: i,index,optim,n_vert_mix

  status = 0 ! no error

!!$if (rank ==0 ) then
!!$   write(*,*)'in adapt_grid'
!!$   write(*,*) 'adapt_par%no_steps', adapt_par%no_steps
!!$   write(*,*) 'adapt_par%dim', adapt_par%dim
!!$   write(*,*) 'adapt_par%nlevs', adapt_par%nlevs
!!$   write(*,*) 'adapt_par%lat_up', adapt_par%lat_up
!!$   write(*,*) 'adapt_par%lat_down', adapt_par%lat_down
!!$   write(*,*) 'adapt_par%lat_min', adapt_par%lat_min
!!$   write(*,*) 'adapt_par%lat_max', adapt_par%lat_max
!!$   write(*,*) 'adapt_par%lev_min', adapt_par%lev_min
!!$   write(*,*) 'adapt_par%lev_max', adapt_par%lev_max
!!$   write(*,*) 'adapt_par%lexp', adapt_par%lexp
!!$   write(*,*) 'adapt_par%timestep', adapt_par%timestep
!!$   write(*,*) 'adapt_par%fac_min', adapt_par%fac_min(ilev)
!!$   write(*,*) 'adapt_par%fac_max', adapt_par%fac_max(ilev)
!!$   write(*,*) 'adapt_par%fac_limit_outside', adapt_par%fac_limit_outside
!!$   write(*,*) 'adapt_par%fac_limit_inside', adapt_par%fac_limit_inside
!!$   write(*,*) 'adapt_par%fac_limit_lev_down', adapt_par%fac_limit_lev_down
!!$   write(*,*) 'adapt_par%fac_limit_lev_up', adapt_par%fac_limit_lev_up
!!$   write(*,*) 'adapt_par%fac_eliminate', adapt_par%fac_eliminate  
!!$   write(*,*) 'adapt_par%r_mean_c', adapt_par%r_mean_c
!!$   write(*,*) 'adapt_par%r_min_c', adapt_par%r_min_c(ilev)
!!$   write(*,*) 'adapt_par%r_max_c', adapt_par%r_max_c(ilev)
!!$   write(*,*) 'adapt_par%r_lim_c_inside', adapt_par%r_lim_c_inside(ilev)
!!$   write(*,*) 'adapt_par%r_lim_c_outside', adapt_par%r_lim_c_outside(ilev)
!!$   write(*,*) 'adapt_par%r_mean_h', adapt_par%r_mean_h
!!$   write(*,*) 'adapt_par%r_min_h', adapt_par%r_min_h(ilev)
!!$   write(*,*) 'adapt_par%r_max_h', adapt_par%r_max_h(ilev)
!!$   write(*,*) 'adapt_par%r_lim_h_inside', adapt_par%r_lim_h_inside(ilev)
!!$   write(*,*) 'adapt_par%r_lim_h_outside', adapt_par%r_lim_h_outside(ilev)
!!$   write(*,*) 'adapt_par%delta_lev', adapt_par%delta_lev
!!$   write(*,*) 'adapt_par%interp_lev', adapt_par%interp_lev
!!$   write(*,*) 'adapt_par%r_dev', adapt_par%r_dev
!!$   write(*,*) 'adapt_par%lev_grid', adapt_par%lev_grid
!!$   write(*,*) 'adapt_par%lev_delta', adapt_par%lev_delta
!!$!   write(*,*) 'adapt_par%r_grid', adapt_par%r_grid
!!$end if

  if (switch<1 .or. switch>7) write(*,*) 'no adaption; wrong switch code(only 1 to 7)'

  if (ctrl_out) print *, 'Number of APs: ', n, rank
!!!!!
!  print *, rank, 'Number of APs: ', n

  if (switch==1 .or. switch==3 .or. switch==4 .or. switch==6) then

     ! eliminating

     i=1
     do while (i<=n)      ! n is NOT ( I repeat NOT) a constant!!!

        optim=1 !to save the index of the neighbour in the neighbourlist of point i

        do while (search_min(ap_s,i,adapt_par,index,optim,ilev))  !check, if point is to near 

           if (lev_gt_min(ap_s,i,index,l_min_act)) then

              if (subset(ap_s,i) .or. subset(ap_s,index)) then  

                 ! eliminate the point, which will be eliminated in all neighbourlists
                 call change_neighbours(ap_s,index)     
                 
                 ! eliminate the point in the structur;  it is replaced by last point
                 call eliminate(ap_s,i,index,adapt_par) 
                 
                 ! change index of last point in all neighbourlists
                 if (index<=n) call change_last_point(ap_s,index)  

              else
                 optim = optim+1

              endif

           else
              optim = optim+1
           endif

           ! else: last point has been eliminated
        end do
        i=i+1

     end do

     if (ctrl_out) print *, 'Number of APs after eliminating: ', n, rank

  endif


  if (switch==2 .or. switch==3 .or. switch==4 .or. switch==6) then

     ! adapting

     i=1
     do while (i<=n)         ! n is NOT ( I repeat NOT) a constant!!!

        optim=1  !to save the index of the neighbour in the neighbourlist of point i

        do while (search_max(ap_s,i,adapt_par,index,optim, ilev))  !check, if point is to far distanced

           if (lev_gt_min(ap_s,i,index,l_min_act)) then  ! if lev_level of point i and index 
                                                           ! greater than l_min_act

              if (subset(ap_s,i) .or. subset(ap_s,index)) then  

                 !create a new point
                 call new_point (status,ap_s,i,adapt_par,index)           
                 if (status/=0) then
                    write (*,*) 'Error in new_point !!!'
                    return
                 endif
                                  
                 !revert the neighbourlists
                 call change_neighbours_of_new_point(ap_s,optim)   
                 
                 call interpolate(ap_s)

              else
                 optim= optim+1
              endif

           else
              optim=optim+1
           endif

        end do

!         if (ctrl_out) then
!            if (mod(i,2000) == 0) print *, 'Adaptation step', i
!         endif
        i=i+1

     end do
     if (ctrl_out) print *, 'Number of APs after adaptation: ', n, rank

  endif


 if (switch==4 .or. switch==5) then
     !vertical unresolved mixing
     i=1
     n_vert_mix=0  
      do while (i<=n) 

        if(search_bvf_min(ap_s,i,adapt_par)) then 
        
        !check if the bvf is too small  
           if(subset(ap_s,i)) then
                 
            ! average the compositions of point i and neighbour index
            call avg_concentration_all(ap_s,i)
            n_vert_mix=n_vert_mix+1 
           endif
        endif
      i=i+1    
      end do
   print *, 'Number of APs with unresolved vert_mixing: ', n_vert_mix,'out of',n  
 endif  

  if (switch==6 .or. switch==7) then
     !vertical resolved mixing
     n_vert_mix=0
     i=1
      do while (i<=n)
        optim=1 ! to save index of the neighbour in neighbourlist of point i
        do while (search_vertical_max(ap_s,i,adapt_par,index,optim))  !check if
           ! the deltatheta is too large
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!           print *,'aps', i,'index', index,'optim', optim
!           if(i<=ap_s%konstant_size) then
!           print *,'theta', ap_s%konstant(i)%theta
!           else
!           print *,'theta', ap_s%dynamic(i-ap_s%konstant_size)%theta
!           endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
           if(subset(ap_s,i) .or. subset(ap_s,index)) then         
            ! average the compositions of point i and neighbour index
            call avg_concentration(ap_s,i,index)
            n_vert_mix=n_vert_mix+1
           endif  
           optim=optim+1
         
        end do
     i=i+1
      end do
  print *, 'Number of APs with resolved vert_mixing: ', n_vert_mix,'out of',n
  endif

end subroutine adapt_grid


end Module messy_clamsmix_lib_mix

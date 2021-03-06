!******************************************************************************
! File utils/src/lib_interpol.f90
!
! This file is part of the Chemical Lagrangian Model of the
! Stratosphere (CLaMS). CLaMS is a hierarchy of numerical models and
! preprocessors which simulate transport, chemical reactions and
! mixing processes in the stratosphere.  
!
! Copyright (c) 2006
! P. Konopka, N. Thomas, Forschungszentrum Juelich GmbH
! Last Modified By: N.Thomas
! Last Modified On: Mon Sep 21 10:14:25 2015
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
! Module lib_interpol
!
! This module contains some interpolation tools:
!
!  Subroutine determ_weight (EckCoor, theta, coor_exp, theta_exp, weight, &
!                           theta_nearest, delta_theta, ratio, wt)
!
!  subroutine determ_weight_old (EckCoor, theta, coor_exp, theta_exp, weight, &
!                          theta_nearest, wt)
!
!  subroutine vert_interpol (theta_down, theta_up, theta, &
!                            value_down, value_up, value, vert_up, vert_down, wt)
!
!  subroutine log_vert_interpol (theta_down, theta_up, theta, &
!                                value_down, value_up, value, &
!                                vert_up, vert_down, wt)
!
!******************************************************************************

module messy_clams_tools_interpol    
  
contains


  !*****************************************************************************
  ! PURPOSE: determine the weights for the interpolation 
  !          wt = 0 (nearest)
  !          wt = 1 (weighted, method 1)
  !          wt = 2 (weighted, method 2)
  !
  !   Input: EckCoor/theta - cart. coordinates/thetas of the 3 APs (triangle)
  !          surrounding the given AP with coordinates given by 
  !          coor_exp/theta_exp     
  !
  !  Output: weight - vector with 3 weights of APs surrounding the given AP
  !          theta_nearest - weighted thetas of these 3 APs
  !
  ! Remarks: for the option wt=0, weight=[0,1,0] with 1 for the
  !          nearest AP
  !          the weights vector is normalized to 1
  !*****************************************************************************
  Subroutine determ_weight (EckCoor, theta, coor_exp, theta_exp, weight, &
                          theta_nearest, delta_theta, ratio, wt)
    
    USE messy_clams_global,         ONLY: prec

    implicit none
    
    real(prec), dimension(3,3), intent(in)  :: EckCoor
    real(prec), dimension(3)                :: theta
    real(prec), dimension(3),   intent(in)  :: coor_exp
    real(prec),                 intent(in)  :: theta_exp
    real(prec), dimension(3),   intent(out) :: weight
    real(prec),                 intent(out) :: theta_nearest
    real(prec),                 intent(in)  :: delta_theta, ratio
    integer,                    intent(in)  :: wt

    ! Set the tolerance in the case if the g-point has approx the same position
    ! as one of the sourronding APs    
    real(prec),parameter    :: tol = 1e-12
    real(prec),dimension(3) :: r                           
    real(prec),dimension(3) :: x_aps,y_aps,z_aps
    integer                 :: min_index(1)
    
    x_aps = EckCoor(1,:) 
    y_aps = EckCoor(2,:)    
    z_aps = EckCoor(3,:)
    
    ! weighted interpolation: contribution proportional to 1/r
    ! i.e. strongest contribution of the nearest APs   

    if (wt==1) then      !weighted (method 1)
       r = abs(theta_exp-theta)*sqrt((coor_exp(1)-x_aps)**2+&
            &(coor_exp(2)-y_aps)**2+(coor_exp(3)-z_aps)**2)
       where(abs(r)<=tol) r = tol
       weight = (1./r)/(1./r(1)+1./r(2)+1./r(3))
    elseif (wt==2) then  !weighted (method 2)
       r = sqrt((coor_exp(1)-x_aps)**2 + (coor_exp(2)-y_aps)**2 + &
                (coor_exp(3)-z_aps)**2 + ((theta_exp-theta)/delta_theta*ratio)**2)
       where(abs(r)<=tol) r = tol
       weight = (1./r)/(1./r(1)+1./r(2)+1./r(3))
    else                 !nearest
       r = sqrt((coor_exp(1)-x_aps)**2 + (coor_exp(2)-y_aps)**2 + &
                (coor_exp(3)-z_aps)**2 + ((theta_exp-theta)/delta_theta*ratio)**2)
       where(abs(r)<=tol) r = tol
       min_index=minloc(r)
       weight = 0.
       weight(min_index(1)) = 1.
    end if

    theta_nearest = weight(1)*theta(1) + weight(2)*theta(2) + weight(3)*theta(3)
    
  end Subroutine determ_weight

  
  !*****************************************************************************
  ! PURPOSE: determine the weights for the interpolation 
  !          (wt = weighted, nearest)
  !
  !   Input: EckCoor/theta - cart. coordinates/thetas of the 3 APs (triangle)
  !          surrounding the given AP with coordinates given by 
  !          coor_exp/theta_exp     
  !
  !  Output: weight - vector with 3 weights of APs surrounding the given AP
  !          theta_nearest - weighted thetas of these 3 APs
  !
  ! Remarks: for the option wt=false, weight=[0,1,0] with 1 for the
  !          nearest AP
  !          the weights vector is normalized to 1
  !*****************************************************************************
  Subroutine determ_weight_old (EckCoor, theta, coor_exp, theta_exp, weight, &
                          theta_nearest, wt)
    
    USE messy_clams_global,         ONLY: prec

    implicit none
    
    real(prec), dimension(3,3), intent(in)  :: EckCoor
    real(prec), dimension(3)                :: theta
    real(prec), dimension(3),   intent(in)  :: coor_exp
    real(prec),                 intent(in)  :: theta_exp
    real(prec), dimension(3),   intent(out) :: weight
    real(prec),                 intent(out) :: theta_nearest
    logical,                    intent(in)  :: wt

    ! Set the tolerance in the case if the g-point has approx the same position
    ! as one of the sourronding APs    
    real(prec),parameter    :: tol = 1e-12
    real(prec),dimension(3) :: r                           
    real(prec),dimension(3) :: x_aps,y_aps,z_aps
    integer                 :: min_index(1)
    
    x_aps = EckCoor(1,:) 
    y_aps = EckCoor(2,:)    
    z_aps = EckCoor(3,:)
    
    r = abs(theta_exp-theta)*sqrt((coor_exp(1)-x_aps)**2+&
         &(coor_exp(2)-y_aps)**2+(coor_exp(3)-z_aps)**2)
    
    where(abs(r)<=tol) r = tol

    ! weighted interpolation: contribution proportional to 1/r
    ! i.e. strongest contribution of the nearest APs   
    if (wt) then                              !weighted
       weight = (1./r)/(1./r(1)+1./r(2)+1./r(3))
    else                                      !nearest
       min_index=minloc(r)
       weight = 0.
       weight(min_index(1)) = 1.
    end if

    theta_nearest = weight(1)*theta(1) + weight(2)*theta(2) + weight(3)*theta(3)
    
  end Subroutine determ_weight_old
  
  !****************************************************************************
  ! PURPOSE: linear interpolation in the theta coordinates
  !   Input: theta_down and value_down at theta_down,
  !          theta_up and value_up at theta_up
  !          theta - theta-level where the interpolation should be determined
  !  Output: value: interpolated value
  !          vert_up, vert_down
  !****************************************************************************
  subroutine vert_interpol (theta_down, theta_up, theta, &
                            value_down, value_up, value, vert_up, vert_down, wt)
    
    USE messy_clams_global,         ONLY: prec

    implicit none
    
    real(prec),intent(in)    :: theta, theta_down, theta_up, value_down, value_up
    real(prec),intent(out)   :: value, vert_up, vert_down
    logical,   intent(in)    :: wt
    
    if (.not. wt) then
       if (abs(theta_up - theta) <= abs(theta - theta_down)) then
          value = value_up
          vert_up = 1.
          vert_down = 0.
       else
          value = value_down
          vert_up = 0.
          vert_down = 1.
       end if
       
    elseif (abs(theta_up - theta_down) <= 0.1) then
       value = value_up
       vert_up = 1.
       vert_down = 0.
       
    elseif (theta <= theta_down) then
       value = value_down
       vert_up = 0.
       vert_down = 1.
       
    elseif (theta >= theta_up) then
       value = value_up
       vert_up = 1.
       vert_down = 0.
       
    else
       value = value_down + &
            (theta-theta_down)*(value_up-value_down)/(theta_up-theta_down)

       vert_up = (theta-theta_down)/(theta_up-theta_down)
       vert_down = 1. - vert_up

    end if
    
  end subroutine vert_interpol

  !****************************************************************************
  ! PURPOSE: log. linear interpolation in the theta coordinates
  !   Input: theta_down and value_down at theta_down,
  !          theta_up and value_up at theta_up
  !          theta - theta-level where the interpolation should be determined
  !  Output: value: interpolated value
  !          vert_up, vert_down
  !****************************************************************************
  subroutine log_vert_interpol (theta_down, theta_up, theta, &
                                value_down, value_up, value, &
                                vert_up, vert_down, wt)

    USE messy_clams_global,         ONLY: prec

    implicit none
    
    real(prec),intent(in)    :: theta, theta_down, theta_up, value_down, value_up
    real(prec),intent(out)   :: value, vert_up, vert_down
    logical,   intent(in)    :: wt


    if (.not. wt) then
       if (abs(theta_up - theta) <= abs(theta - theta_down)) then
          value = value_up
          vert_up = 1.
          vert_down = 0.
       else
          value = value_down
          vert_up = 0.
          vert_down = 1.
       end if
       
    elseif (abs(theta_up - theta_down) <= 0.1) then
       value = value_up
       vert_up = 1.
       vert_down = 0.
       
    elseif (theta <= theta_down) then
       value = value_down
       vert_up = 0.
       vert_down = 1.
       
    elseif (theta >= theta_up) then
       value = value_up
       vert_up = 1.
       vert_down = 0.
       
    else

       value = value_down + &
            (log(theta)-log(theta_down))*(value_up-value_down)/(log(theta_up)-log(theta_down))

       vert_up = (log(theta)-log(theta_down))/(log(theta_up)-log(theta_down))
       vert_down = 1. - vert_up

    end if


  end subroutine log_vert_interpol

!*************************************************************
! Compute linear interpolation
!*************************************************************
function interpol_lin (x1, y1, x2, y2, x)

  use messy_clams_global, only: prec
  
  implicit none

  real(prec) :: interpol_lin, x1, y1, x2, y2, x
  
  if (x1==x2) then
     interpol_lin = y1
  else     
     interpol_lin = y1 + (y2-y1)/(x2-x1)*(x-x1)
  endif
  
end function interpol_lin

!!!!! ???

function interpol_lin_double (x1, y1, x2, y2, x)
  
  use messy_clams_global, only: prec, dp

  implicit none

  real(prec)  :: interpol_lin_double, y1, y2
  real(dp)    :: x1, x2, x
  
  if (x1==x2) then
     interpol_lin_double = y1
  else
     interpol_lin_double = y1 + (y2-y1)/(x2-x1)*(x-x1)
  endif
  
end function interpol_lin_double




end module messy_clams_tools_interpol


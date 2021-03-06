!**************************************************************************
!
! This file is part of the Chemical Lagrangian Model of the
! Stratosphere (CLaMS). CLaMS is a hierarchy of numerical models and
! preprocessors which simulate transport, chemical reactions and
! mixing processes in the stratosphere.  
!
! Copyright (c) 2020
! Mengchu Tao, Paul Konopka, Nicole Thomas
! Forschungszentrum Juelich GmbH
! Last Modified By: N.Thomas
! Last Modified On: Tue Apr  7 11:51:11 2020
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
!********************************************************************************!
!
!    Subroutine deep_conv (status, lat, lon, zeta, theta,  &
!                          temp, press, bvf_wet, h2o, dc)
!
!    subroutine calc_delta_theta (lat, lon, zeta, theta,  &
!                                 temp, bvf_wet, h2o, dc, ntraj)
!
!    subroutine calc_temp (theta, press, temp, ntraj)
!
!    Subroutine clamsdeepconv_read_nml (status, iou)
!
!********************************************************************************!

Module messy_clamsdeepconv

  IMPLICIT NONE

  CHARACTER(LEN=*),PARAMETER, PUBLIC :: modstr = 'clamsdeepconv'
  CHARACTER(LEN=*),PARAMETER, PUBLIC :: modver = '1.0'

CONTAINS 

  
  !**************************************************************************
  !
  !**************************************************************************
  Subroutine deep_conv (status, lat, lon, zeta, theta,  &
                        temp, press, bvf_wet, h2o, dc)


    USE messy_clams_global,       ONLY: prec, mdi, dnparts, rank
    
    implicit none

    ! Dimension (1:dnparts_max) 
    REAL(PREC),  DIMENSION(:), POINTER :: lat, lon, zeta
    REAL(PREC),  DIMENSION(:), POINTER :: theta
    REAL(PREC),  DIMENSION(:), POINTER :: temp, press, bvf_wet, h2o
    REAL(PREC),  DIMENSION(:), POINTER :: dc
    
    integer :: status

    if (rank==0) THEN
       WRITE (*,*) 
       WRITE (*,*) 'START OF DEEPCONV'
       WRITE (*,*)
    endif

    ! calculate delta theta
    call calc_delta_theta (lat, lon, zeta, theta,  &
                           temp, bvf_wet, h2o, dc, dnparts)

    ! interpolate ZETA and PRESS for new THETA
    call interpol_params (lat, lon, theta, zeta, press, dnparts)

    ! calc TEMP from PRESS and new THETA 
    call calc_temp (theta, press, temp, dnparts)

    
  End Subroutine deep_conv


  
  !*******************************************************************
  ! calculate delta theta
  !*******************************************************************
  subroutine calc_delta_theta (lat, lon, zeta, theta,  &
                               temp, bvf_wet, h2o, dc, ntraj)

    USE messy_clams_global,         ONLY: prec, rank, mdi, eps
    USE messy_clamsdeepconv_global, ONLY: dc_kind, cbvf, crate, min_dtheta, &
                                          lv, cp, zeta_min
    
    implicit none

    real(prec), dimension(:), pointer :: lat, lon, zeta
    real(prec), dimension(:), pointer :: theta, temp, bvf_wet, h2o
    real(prec), dimension(:), pointer :: dc
    
    integer                           :: ntraj
    
    integer       :: counter, itraj
    real(prec)    :: dtheta
    
    if (rank==0) write (*,*) 'Calculate delta-theta'

    if (DC_kind == 0.) then
       if (rank==0) write (*,*) 'DC tracer set to zero'
       DC          = mdi
       DC(1:ntraj) = 0.
    endif

    ! if lat or lon missing, then those theta are valued missing
    where (ABS((lat-mdi)/mdi)<=eps)
       lon   = mdi
       theta = mdi
       zeta  = mdi
    end where
    where (ABS((lon-mdi)/mdi)<=eps)
       lat   = mdi
       theta = mdi
       zeta  = mdi
    end where

    counter = 0
    
    do itraj = 1, ntraj
          
       if (bvf_wet(itraj) < cbvf  .and. &
           zeta(itraj) < zeta_min .and. &
           ABS((h2o(itraj)-mdi)/mdi)>eps .and. &
           ABS((bvf_wet(itraj)-mdi)/mdi)>eps .and. &
           ABS((theta(itraj)-mdi)/mdi)>eps .and. &
           ABS((zeta(itraj)-mdi)/mdi)>eps .and. &
           ABS((temp(itraj)-mdi)/mdi)>eps) then
             
          dtheta = (18. * lv * h2o(itraj) * theta(itraj)) / &
                   (29. * cp * temp(itraj))
             
          if(dtheta>min_dtheta) then
             theta(itraj) = theta(itraj) + crate*dtheta
             DC   (itraj) = 1.
             counter = counter + 1
          endif
          
       endif
       
    enddo

    write (*,*) rank, 'Number of traj with bvf_wet<cbvf: ', counter
    write (*,*) rank, 'Number of traj with DC>0: ', count(DC > 0)


  end subroutine calc_delta_theta


  !*******************************************************************
  ! Get ZETA and PRESS for new THETA values
  !*******************************************************************
  subroutine interpol_params (lat, lon, theta, zeta, press, ntraj)

    USE messy_clams_global,            ONLY: prec, rank, ldiagout, &
                                             level_is_vertcoor, &
                                             longrid, latgrid, nx, ny, nz, &
                                             asc_lat, &
                                             predata, leveldt, &
                                             nparams, paramnames
    
    USE messy_clams_tools_interpolreg, ONLY: interpolate_spatial, &
                                             interpolate_spatial_modlev
    USE messy_clams_tools_utils,       ONLY: str_pos
   
    implicit none

    real(prec), dimension(:), pointer :: lat, lon, theta, zeta, press
    integer                           :: ntraj

    integer      :: i, ind_theta, ind_press
    logical      :: asclev, loglev, extrapol
   
    if (rank==0) write (*,*) 'Interpolate ZETA and PRESS with new THETA'


    loglev = .false.   ! level: theta -> not logarithmic
    asclev = .true.    !              -> ascending levels
    extrapol = .false. !              -> no extrapolation
  
    ! thetadt : einlesen (leveldt enthaelt ZETA-Werte!) -> predata(theta) ?
    ! zeta_met, press_met:  auf PREDATA/FUTDATA ?!?

    ind_theta = str_pos(nparams,paramnames,"THETA")
    ind_press = str_pos(nparams,paramnames,"PRESS")
    if (rank==0 .and. ldiagout) write (*,*) 'ind_theta, ind_press = ', &
                                            ind_theta, ind_press

!!!!! only for data on model levels !?!
    ! ZETA is vertical coordinate in initfile
    ! -> ZETA from windfile read in to leveldt
    do i = 1, ntraj
       zeta(i) = interpolate_spatial_modlev (leveldt, &
                        predata(ind_theta)%values, longrid, latgrid, &
                        nx, ny, nz, lon(i), lat(i), theta(i), &
                        asc_lat, asclev=asclev, loglev=loglev, logval=0, &
                        extrapol=extrapol)
       press(i) = interpolate_spatial_modlev (predata(ind_press)%values, &
                        predata(ind_theta)%values, longrid, latgrid, &
                        nx, ny, nz, lon(i), lat(i), theta(i), &
                        asc_lat, asclev=asclev, loglev=loglev, logval=2, &
                        extrapol=extrapol)
    enddo


  end subroutine interpol_params

  
  !******************************************************************************
  !   2018.05 (M.Tao)
  !   calc temp from press and new theta
  !******************************************************************************
  subroutine calc_temp (theta, press, temp, ntraj)

    USE messy_clams_global,   ONLY: prec, mdi, eps
    
    implicit none
    
    real(prec), dimension(:), pointer :: theta, press, temp
    integer                           :: ntraj
    
    integer                   :: itraj
    real(prec)                :: p_s, r0
    
    p_s = 1000.
    r0  = 0.286
    
    do itraj = 1, ntraj

       if (abs((theta(itraj)-mdi)/mdi) <= eps .or. abs((press(itraj)-mdi)/mdi) <= eps) then
          temp(itraj) = mdi
       else
          temp(itraj) = theta(itraj)/((p_s/press(itraj))**r0)
       endif
       
    enddo
    
  end subroutine calc_temp
  
  !**************************************************************************
  !
  !**************************************************************************
  Subroutine clamsdeepconv_read_nml (status, iou)
      
    USE messy_main_tools,    ONLY: read_nml_open, read_nml_check, read_nml_close
    USE messy_clams_global,  ONLY: rank, ldiagout
    USE messy_clamsdeepconv_global
    
    IMPLICIT NONE
    
    !I/O
    INTEGER, INTENT(OUT) ::   status
    INTEGER, INTENT(IN)  ::   iou
    
    !LOCAL
    CHARACTER(LEN=*),PARAMETER :: substr='clamsdeepconv_read_nml'
    LOGICAL              :: lex      ! file exists ?
    INTEGER              :: fstat    ! file status
    LOGICAL              :: l_print  ! write control output
    
    
    NAMELIST /CTRL/ cbvf, crate, min_dtheta, DC_kind, timestep_deepconv

    status = 1 !ERROR
    
    if (rank==0 .and. ldiagout) then
       l_print = .true.
    else
       l_print = .false.
    endif
    
    !-------------------------------------------------------------------
    ! Read namelist variables:
    !-------------------------------------------------------------------
    
    CALL read_nml_open(lex, substr, iou, 'CTRL', modstr, l_print)
    IF (.not.lex) RETURN    ! <modstr>.nml does not exist
    
    READ(iou, NML=CTRL, IOSTAT=fstat)
    CALL read_nml_check(fstat, substr, iou, 'CTRL', modstr, l_print)
    IF (fstat /= 0) RETURN  ! error while reading namelist
      
    CALL read_nml_close(substr, iou, modstr, l_print)

    status = 0 !NO ERROR
 
  End Subroutine clamsdeepconv_read_nml

End Module messy_clamsdeepconv
   

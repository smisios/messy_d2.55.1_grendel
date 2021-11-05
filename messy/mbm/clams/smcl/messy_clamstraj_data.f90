!********************************************************************************!
!
! This file is part of the Chemical Lagrangian Model of the
! Stratosphere (CLaMS). CLaMS is a hierarchy of numerical models and
! preprocessors which simulate transport, chemical reactions and
! mixing processes in the stratosphere.  
!
! Copyright (c) 2006
! Nicole Thomas
! Forschungszentrum Juelich GmbH
! Last Modified By: N.Thomas
! Last Modified On: Wed Jun 29 10:49:01 2016
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
! MODULE messy_clamstraj_data
! ----------------------------          
!
! This MODULE contains subroutines for handling of data and memory
!
! SUBROUTINE allocate_arrays
! SUBROUTINE deallocate_arrays
! SUBROUTINE get_sza (lons, lats, calc, itime, paramname)
! SUBROUTINE get_positions_and_params (calc, itime, ifutime,  &
!                                      plats, plons, plevs, param)          
!
!********************************************************************
                         
MODULE  messy_clamstraj_data

CONTAINS

!*******************************************************************
!
!*******************************************************************
SUBROUTINE allocate_arrays

  USE messy_clams_global,     ONLY: dnparts_max, levelno

  implicit none

  allocate (levelno(dnparts_max))

END SUBROUTINE allocate_arrays

!*******************************************************************
!
!*******************************************************************
SUBROUTINE deallocate_arrays

  USE messy_clams_global,  ONLY: levelno

  implicit none

  DEALLOCATE (levelno)

END SUBROUTINE deallocate_arrays


!*********************************************************************
!  get current solar zenith angle for a given  date and time <itime>
!  and positions given by an array of latitudes <lats> and 
!  longitudes <lons>.
!*********************************************************************
SUBROUTINE get_sza (lons, lats, calc, itime)

   USE messy_clams_global,     ONLY: dnparts_max, dnparts, timetype, &
                                     prec, mdi, eps
   USE messy_clamstraj_calc3d, ONLY: ZEN2

   IMPLICIT NONE

   REAL(PREC)     :: lons(dnparts_max), lats(dnparts_max)
   LOGICAL        :: calc(dnparts_max)
   TYPE(timetype) :: itime

   ! local variables
   REAL(PREC)    :: sza(dnparts_max)
   REAL(PREC)    :: time1,cozen1
   INTEGER       :: iparts
   INTEGER       :: date1(3)

   date1(1) = itime%year
   date1(2) = itime%month
   date1(3) = itime%day
   time1    = itime%sec/86400.

   DO iparts = 1, dnparts
      IF ((abs((LATS(iparts)-MDI)/MDI)>eps) .AND. &
           (ABS((LONS(iparts)-MDI)/MDI)>eps).and. calc(iparts)) THEN 
         CALL ZEN2(Lats(iparts),lons(iparts),date1,time1,COZEN1,sza(iparts))
      ELSE
         sza(iparts) = MDI        
      ENDIF    
   ENDDO

END SUBROUTINE get_sza

!*********************************************************************
! Convert from polar-stereographic to latitudes and longitudes
! and interpolate parameters at these positions
!*********************************************************************
SUBROUTINE get_positions_and_params (calc, itime, ifutime,  &
                                     plats, plons, plevs, param, lcoupled)          
 
   USE messy_clams_global,            ONLY: prec, mdi, eps, &
                                            dnparts_max, dnparts, &
                                            timetype, dates30, irdday, irdsec, &
                                            nparams, &
                                            paramtype
   USE messy_clamstraj_global,        ONLY: trajtype
   USE messy_clams_tools_interpolreg, ONLY: interpolate_param
   USE messy_clams_tools_dateconv,    ONLY: ymds2js_interface, incdat_interface
   USE messy_clams_tools_utils,       ONLY: uppercase
   USE messy_clamstraj_calc3d,        ONLY: fixll
   
   IMPLICIT NONE

   TYPE(timetype) :: itime, ifutime
   LOGICAL        :: calc(dnparts_max)
   LOGICAL        :: lcoupled

   REAL(PREC)     :: plats(dnparts_max), plons(dnparts_max), plevs(dnparts_max)
   TYPE(paramtype), DIMENSION(:), POINTER :: param

   ! local variables
   TYPE(timetype)  :: ipasttime
   REAL(PREC)      :: lons(dnparts_max), lats(dnparts_max), levs(dnparts_max)
   REAL(PREC)      :: pi
   INTEGER         :: i, ipart

   pi=4.*atan(1.)                                                        
  
   ! Adjusts longitudes and latitudes in plons and plats to ensure that      
   ! all values are in range  (0,2*pi) and  (-.5*pi,.5*pi)      
   call fixll (plons, plats, size(plons), 1, size(plons), pi) 

   ! convert from polar-stereographic to latitudes and longitudes
   ! and from ln(lev) to lev
   DO i = 1, dnparts
      IF ((abs((plons(i)-MDI)/MDI)<=eps) .OR. (.NOT. calc(i) .AND. trajtype/=1)) THEN
         lons(i) = MDI
         lats(i) = MDI
         levs(i) = MDI
      ELSE
! op_pj_20160606+
!!$         lons(i) = mod(plons(i)*180./pi,360.) 
         lons(i) = mod(plons(i)*180._prec/pi,360._prec) 
! op_pj_20160606-
         lats(i) = plats(i)*180./pi                           
         levs(i) = plevs(i)
      ENDIF
   ENDDO
   
   ! get last data time (for predata)
   ipasttime%sec = ifutime%sec
   ipasttime%day = ifutime%day
   ipasttime%month = ifutime%month
   ipasttime%year = ifutime%year
   CALL incdat_interface( ipasttime%sec,ipasttime%day,ipasttime%month,ipasttime%year, &
        -irdsec ,-irdday,0,0, dates30)     

   ! Interpolate parameter values at given positions and time
   DO i=1, nparams
      IF (param(i)%name == 'SZA') THEN
         CALL get_sza (lons,lats,calc,itime)
      ELSE
         CALL interpolate_param (lons,lats,levs,itime,ipasttime,i, &
              param(i)%name, param(i)%values, lcoupled, calc)          
      ENDIF

!!!!! Nur fuer TEMP und PRESS werden bei missing-values die Positionen 
!!!!! auch auf MDI gesetzt !
!!!!! => Sonst kann die Anzahl der Trajektorien je nach Parameterliste unterschiedlich sein!
      IF (param(i)%name=='TEMP' .or. param(i)%name=='PRESS') THEN
         DO ipart = 1, dnparts
            IF (abs((param(i)%values(ipart)-MDI)/MDI)<=eps) THEN
               plons(ipart) = MDI
               plats(ipart) = MDI
               plevs(ipart) = MDI
            END IF
         END DO
      ENDIF

   ENDDO
    
 END SUBROUTINE get_positions_and_params

END MODULE messy_clamstraj_data





!********************************************************************************!
! File clams/traj/source/global.f90
!
! This file is part of the Chemical Lagrangian Model of the
! Stratosphere (CLaMS). CLaMS is a hierarchy of numerical models and
! preprocessors which simulate transport, chemical reactions and
! mixing processes in the stratosphere.  
!
! Copyright (c) 2006
! Nicole Thomas, Forschungszentrum Juelich GmbH
! Last Modified By: jicg1108
! Last Modified On: Thu Nov 15 14:03:12 2018
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
! MODULE global
! -------------
! Global definitions for program traj
!
!********************************************************************************!

MODULE messy_clamstraj_global

  USE messy_clams_global, only: sp, dp, prec, timetype, filenamelen
                                

!---------------------------------------------------------------------------
! Namelist variables
!---------------------------------------------------------------------------

  INTEGER, PUBLIC         :: type_traj
  INTEGER, PUBLIC         :: idtsec2   ! fine resolution timestep in seconds

!---------------------------------------------------------------------------
!
!---------------------------------------------------------------------------

  INTEGER, PUBLIC              :: time_init_julsec
  REAL(DP), PUBLIC             :: p_ref=300    ! ju_ec 20190130 gave default
                                               ! value, units 300 hPa



!!!!!
   ! julsec    : starttime of trajectories
   ! endjulsec : endtime of trajectories 
!   REAL(DP), POINTER :: julsec(:)!, endjulsec(:)

   ! startvalues: positions (lon,lat,lev) and parameter values on starttime
   ! endvalues  : positions and parameter values on endtime
!   TYPE(valuestype) :: startvalues, endvalues

   ! startfirsttraj :: starttime of earliest starting trajectory (in julian seconds)
   ! startlasttraj  :: starttime of latest starting trajectory (in julian seconds)
   ! endfirsttraj   :: endtime of earliest ending trajectory (in julian seconds)
   ! endlasttraj    :: endtime of latest ending trajectory (in julian seconds)
   REAL(DP) :: startfirsttraj, startlasttraj, endfirsttraj, endlasttraj

   ! dt       : main time step
   ! dt2      : fine resolution time step
   ! pslat    : latitude to switch to polar stereo graphic grid
   REAL(PREC)    :: dt, dt2, pslat 

   ! idmday       : writing frequency in days
   ! idmsec       : writing frequency in seconds
   ! idtsec       : main timestep in seconds
   INTEGER :: idmday, idmsec, idtsec
              
   ! trajtype : 1 = the same starttime and endtime for all trajectories
   !            2 = different starttimes but the same endtime for all trajectories 
   !            3 = different starttimes and different endtimes for trajectories 
   !                but the same length (isochronic traj.)
   ! nrhours  : length of trajectories in hours (for isochronic traj.)
   INTEGER :: trajtype = 1
   INTEGER :: nrhours

   ! starttime : starttime of calculation
   ! endtime   : endtime of calculation
   TYPE(timetype) :: starttime, endtime

   ! forward : true, if forward trajectories
   LOGICAL :: forward

   ! Logical to switch interpolation on/off:
   LOGICAL :: linterpolate

   ! timestep for trajectory output [hours]
   INTEGER :: timestep_trajout = 0

   ! timestep for trajectory input [hours]
   INTEGER :: timestep_trajin = 0

   ! directory with trajectory input files
   CHARACTER(filenamelen) :: dir_trajin=''

END MODULE messy_clamstraj_global



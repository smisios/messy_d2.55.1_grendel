!********************************************************************************!
!
! This file is part of the Chemical Lagrangian Model of the
! Stratosphere (CLaMS). CLaMS is a hierarchy of numerical models and
! preprocessors which simulate transport, chemical reactions and
! mixing processes in the stratosphere.  
!
! Copyright (c) 2006
! Richard Swinbank, Daniel McKenna, Nicole Thomas, Paul Konopka, 
! Daniela Ciupa, Juergen Ankenbrand, Barbara Deutsch
! Forschungszentrum Juelich GmbH
! Last Modified By: N.Thomas
! Last Modified On: Wed Jun 29 10:50:17 2016
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
! Module messy_clamstraj_timestep contains the following subroutines:
!
! SUBROUTINE TIMSTP
! SUBROUTINE NEXTPOS
!*******************************************************************************

MODULE messy_clamstraj_timestep


CONTAINS

!*******************************************************************
! TIMESTP is the main program calculating the trajectories
!*******************************************************************
SUBROUTINE TIMSTP (status, plats, plons, plevs, julsec, param, lcoupled)          

  USE messy_clams_global,         ONLY: prec, dp, mdi, eps, &
                                        dnparts_max, timetype, dates30, &
                                        nx, ny, nz, nparams, &
                                        fut_sec, fut_day, fut_month, fut_year, &
                                        paramtype
   USE messy_clamstraj_global,     ONLY: starttime, endtime, &
                                        idmsec, idmday, forward, dt, dt2, &
                                        idtsec,  &
                                        trajtype, startfirsttraj, startlasttraj, &
                                        endfirsttraj, endlasttraj, linterpolate
  USE messy_clams_tools_dateconv, ONLY: ymds2js_interface, incdat_interface
  USE messy_clamstraj_data,       ONLY: get_positions_and_params
  USE messy_clams_global,         ONLY: lstart

   IMPLICIT NONE

   REAL(PREC) :: plats(dnparts_max), plons(dnparts_max), plevs(dnparts_max)
   REAL(DP)   :: julsec(dnparts_max)

   TYPE(paramtype), DIMENSION(:), POINTER :: param

   INTEGER    :: status

   REAL(PREC) :: dt3
   INTEGER    :: idtsec3 
   LOGICAL    :: lcoupled

   ! itime   : current time
   ! ifutime : next data time
   ! iwrtime : next writing time
   TYPE(timetype) :: itime, ifutime, iwrtime

   REAL(DP)   :: indjultime, ijultime, wrjultime,   &
                      ijulnext, ijulold
   REAL(PREC) :: pi, er

   INTEGER    :: nsteps, i
   LOGICAL    :: calc(dnparts_max)

   er=6.371e6
   pi=4.*atan(1.)                                                        

   calc = .FALSE.

   status = 0 ! no error 

   ! Initialize 

   itime = starttime
   ifutime = starttime
   iwrtime = starttime
   
   indjultime = ymds2js_interface (endtime%year, endtime%month, endtime%day, endtime%sec, dates30)
   ijultime = ymds2js_interface (itime%year, itime%month, itime%day, itime%sec, dates30)

   ! nsteps = Number of written timesteps
   nsteps = 0

!!$ if (.not. lcoupled) then
   ifutime%sec   = fut_sec
   ifutime%day   = fut_day
   ifutime%month = fut_month
   ifutime%year  = fut_year
!!$  endif

   !----------------------------------------------------------------------------------
   ! If all trajectories start at the same time, write the initial particle positions
   !----------------------------------------------------------------------------------
   
   ! If trajectories start on calculation starttime, write out 
   ! start positions 
   ! If trajectories don't start on calculation starttime, the
   ! start positions are written out when the trajectory starttime
   ! is reached (in subroutine NEXTPOS)

   IF ((trajtype==1) .AND. (ijultime==startfirsttraj)) THEN

      nsteps = nsteps + 1
      calc = .TRUE.

!!!!!
!!$      IF (lstart) THEN
!!$         CALL get_positions_and_params (calc, itime, ifutime,  &
!!$              plats, plons, plevs, param)          
!!$
!!$        DO i = 1, nparams
!!$            param_old(i)%values = param(i)%values
!!$         ENDDO
!!$      ENDIF

   ENDIF

   ! increment the writing date and time by the writing frequency             
   call incdat_interface (iwrtime%sec,iwrtime%day,iwrtime%month,iwrtime%year, &
        idmsec,idmday,0,0, dates30)          

   !--------------------------------------------------------------
   ! forward trajectories 
   !--------------------------------------------------------------

   IF (forward) THEN   

      ! DO ! integrate until first starting time is reached (up to dt)
      !    IF (ijultime+dt > startfirsttraj) EXIT
      !    CALL NEXTPOS (status, nsteps, itime, ifutime, iwrtime, &
      !                  dt, idtsec, calc, plats, plons, plevs, julsec, param)
      !    ijultime = ymds2js_interface (itime%year, itime%month, itime%day, itime%sec, dates30)
      ! ENDDO

      ! IF (trajtype /= 1) THEN !Trajectories start at different times

      !    ijulold = ijultime
      !    ijulnext = ijultime
      !    ! integrate until first traj is reached (up to dt2, dt2 << dt)
      !    ! and use for this integration only one step dt3
      !    DO 
      !       IF (ijulnext+dt2 > startfirsttraj) EXIT
      !       ijulnext = ijulnext + dt2
      !    ENDDO
      !    dt3 = ijulnext - ijultime
      !    idtsec3 = dt3
      !    IF (idtsec3 /= 0) THEN 
      !       CALL NEXTPOS (status, nsteps, itime, ifutime, iwrtime, &
      !                  dt3, idtsec3, calc, plats, plons, plevs, julsec, param) 
      !       ijultime = ymds2js_interface (itime%year, itime%month, itime%day, itime%sec, dates30)
      !    ENDIF
      !    ! integrate until last traj is reached (up to dt2)
      !    ! using fine integration steps dt2
      !    DO 
      !        IF (ijultime >= startlasttraj) EXIT
      !        CALL NEXTPOS (status, nsteps, itime, ifutime, iwrtime, &
      !                  dt2, idtsec2, calc, plats, plons, plevs, julsec, param)
      !       ijultime = ymds2js_interface (itime%year, itime%month, itime%day, itime%sec, dates30)
      !    ENDDO
 
      !    ! adjust the last time step in order to get the last integration step
      !    ! as t_i=t_start+n*dt, t_i > startlast
      !    DO 
      !       IF (ijulold > ijultime) EXIT 
      !        ijulold = ijulold + dt
      !    ENDDO
         
      !    dt3 = ijulold - ijultime
      !    idtsec3 = dt3

      !    IF (ijultime < indjultime) THEN
      !       CALL NEXTPOS (status, nsteps, itime, ifutime, iwrtime, &
      !            dt3, idtsec3, calc, plats, plons, plevs, julsec, param)
      !       ijultime = ymds2js_interface (itime%year, itime%month, itime%day, itime%sec, dates30)

      !    ENDIF

      ! ELSE      !All Trajectories start at the same time

      !    !adjust the integration step t_i > startlast, t_i=t_start+n*dt
      !    dt3 = startfirsttraj -ijultime
      !    idtsec3 = dt3
      !    IF (idtsec3 /= 0) THEN
      !       CALL NEXTPOS (status, nsteps, itime, ifutime, iwrtime, &
      !                  dt3, idtsec3, calc, plats, plons, plevs, julsec, param)
      !       ijultime = ymds2js_interface (itime%year, itime%month, itime%day, itime%sec, dates30)
      !    ENDIF

      !    dt3 = dt - dt3
      !    idtsec3 = dt3
      !    IF (idtsec3 /= 0) THEN
      !       CALL NEXTPOS (status, nsteps, itime, ifutime, iwrtime, &
      !                  dt3, idtsec3, calc, plats, plons, plevs, julsec, param)
      !       ijultime = ymds2js_interface (itime%year, itime%month, itime%day, itime%sec, dates30)
      !    ENDIF

      ! ENDIF
     
      IF ((trajtype==1) .OR. (trajtype==2)) THEN
         DO 
            IF (ijultime >= indjultime) EXIT ! main integration
            CALL NEXTPOS (status, nsteps, itime, ifutime, iwrtime, &
                          dt, idtsec, calc, plats, plons, plevs, julsec, param, lcoupled)
            ijultime = ymds2js_interface (itime%year, itime%month, itime%day, itime%sec, dates30) 
         ENDDO

      ! ELSE
      !    DO
      !       IF (ijultime+dt > endfirsttraj) EXIT ! main integration
      !       CALL NEXTPOS (status, nsteps, itime, ifutime, iwrtime, &
      !                     dt, idtsec, calc, plats, plons, plevs, julsec, param)
      !       ijultime = ymds2js_interface (itime%year, itime%month, itime%day, itime%sec, dates30)
      !    ENDDO

      !    ! adjust the endpoints similar as the startpoints
      !    ijulold = ijultime
      !    ijulnext = ijultime
      !    DO
      !       IF (ijulnext+dt2 > endfirsttraj) EXIT
      !       ijulnext = ijulnext + dt2
      !    ENDDO
      !    dt3 = ijulnext - ijultime
      !    idtsec3 = dt3
      !    IF (idtsec3 /= 0) THEN
      !       CALL NEXTPOS (status, nsteps, itime, ifutime, iwrtime, &
      !                  dt3, idtsec3, calc, plats, plons, plevs, julsec, param)
      !       ijultime = ymds2js_interface (itime%year, itime%month, itime%day, itime%sec, dates30)
      !    ENDIF

      !    DO 
      !       IF (ijultime >= endlasttraj) EXIT
      !       CALL NEXTPOS (status, nsteps, itime, ifutime, iwrtime, &
      !                     dt2, idtsec2, calc, plats, plons, plevs, julsec, param)
      !       ijultime = ymds2js_interface (itime%year, itime%month, itime%day, itime%sec, dates30)
      !    ENDDO

      ENDIF

   !--------------------------------------------------------------
   ! backward trajectories
   !--------------------------------------------------------------
!    ELSE               

! !      WRITE (*,*) 'Backwards'

!       DO 
!          IF (ijultime+dt < startlasttraj) EXIT
!          CALL NEXTPOS (status, nsteps, itime, ifutime, iwrtime, &
!                        dt, idtsec, calc, plats, plons, plevs, julsec, param)
!          ijultime = ymds2js_interface (itime%year, itime%month, itime%day, itime%sec, dates30)
!       ENDDO

!       IF (trajtype /= 1) THEN ! trajectories start at different times
         
!          ijulold = ijultime
!          ijulnext = ijultime
!          DO 
!             IF (ijulnext+dt2 < startlasttraj) EXIT
!             ijulnext = ijulnext + dt2
!          ENDDO
!          dt3 = ijulnext - ijultime
!          idtsec3 = dt3
!          IF (idtsec3 /= 0) THEN
!             CALL NEXTPOS (status, nsteps, itime, ifutime, iwrtime, &
!                        dt3, idtsec3, calc, plats, plons, plevs, julsec, param)

!             ijultime = ymds2js_interface (itime%year, itime%month, itime%day, itime%sec, dates30)
!          ENDIF
         
!          DO
!             IF (ijultime <= startfirsttraj) EXIT
!             CALL NEXTPOS (status, nsteps, itime, ifutime, iwrtime, &
!                        dt2, idtsec2, calc, plats, plons, plevs, julsec, param)

!              ijultime = ymds2js_interface (itime%year, itime%month, itime%day, itime%sec, dates30)
!          ENDDO
         
!          DO 
!              IF (ijulold < ijultime) EXIT
!             ijulold = ijulold + dt
!          ENDDO
      
!          dt3 = ijulold - ijultime
!          idtsec3 = dt3
         
!          IF (ijultime > indjultime) THEN
!             CALL NEXTPOS (status, nsteps, itime, ifutime, iwrtime, &
!                  dt3, idtsec3, calc, plats, plons, plevs, julsec, param)
!             ijultime = ymds2js_interface (itime%year, itime%month, itime%day, itime%sec, dates30)
!          ENDIF

!       ELSE ! all trajectories start at the same time

!          dt3 = startlasttraj-ijultime
!          idtsec3 = dt3
!          IF (idtsec3 /= 0) THEN
!             CALL NEXTPOS (status, nsteps, itime, ifutime, iwrtime, &
!                        dt3, idtsec3, calc, plats, plons, plevs, julsec, param)
!             ijultime = ymds2js_interface (itime%year, itime%month, itime%day, itime%sec, dates30)
!          ENDIF

!          dt3 = dt - dt3
!          idtsec3 = dt3
!          IF (idtsec3 /= 0) THEN
!             CALL NEXTPOS (status, nsteps, itime, ifutime, iwrtime, &
!                        dt3, idtsec3, calc, plats, plons, plevs, julsec, param)
!             ijultime = ymds2js_interface (itime%year, itime%month, itime%day, itime%sec, dates30)
!          ENDIF
!       ENDIF
         
!       IF ((trajtype==1) .OR. (trajtype==2)) THEN

!          DO 
!             IF (ijultime <= indjultime) EXIT
!             CALL NEXTPOS (status, nsteps, itime, ifutime, iwrtime, &
!                           dt, idtsec, calc, plats, plons, plevs, julsec, param)
!             ijultime = ymds2js_interface (itime%year, itime%month, itime%day, itime%sec, dates30)
!          ENDDO

!       ELSE
!          DO
!             IF (ijultime+dt < endlasttraj) EXIT
!             CALL NEXTPOS (status, nsteps, itime, ifutime, iwrtime, &
!                           dt, idtsec, calc, plats, plons, plevs, julsec, param)
!             ijultime = ymds2js_interface (itime%year, itime%month, itime%day, itime%sec, dates30)
!          ENDDO

!          ijulold = ijultime
!          ijulnext = ijultime
!          DO 
!             IF (ijulnext+dt2 < endlasttraj) EXIT
!             ijulnext = ijulnext + dt2
!          ENDDO
!          dt3 = ijulnext - ijultime
!          idtsec3 = dt3
!          IF (idtsec3 /= 0) THEN
!             CALL NEXTPOS (status, nsteps, itime, ifutime, iwrtime, &
!                        dt3, idtsec3, calc, plats, plons, plevs, julsec, param)
!             ijultime = ymds2js_interface (itime%year, itime%month, itime%day, itime%sec, dates30)
!          ENDIF
          
!          DO 
!             IF (ijultime <= endfirsttraj) EXIT
!             CALL NEXTPOS (status, nsteps, itime, ifutime, iwrtime, &
!                           dt2, idtsec2, calc, plats, plons, plevs, julsec, param)
!             ijultime = ymds2js_interface (itime%year, itime%month, itime%day, itime%sec, dates30)
!          ENDDO

!       ENDIF

   ENDIF

   !---------------------------------------------------------------------
   ! If trajectories start at different times, write starting times
   !---------------------------------------------------------------------
!!$   IF ((trajtype==2) .OR. (trajtype==3)) THEN
!!$      CALL write_startpos
!!$   ENDIF
!!$
   !--------------------------------------------------------------
   ! If all trajectories end on the same time and endtime is not a 
   ! writing time: write positions on endtime
   !--------------------------------------------------------------
   IF ((trajtype==1) .OR. (trajtype==2)) THEN

      CALL incdat_interface (iwrtime%sec,iwrtime%day,iwrtime%month,iwrtime%year, &
           -idmsec,-idmday,0,0, dates30)         
      wrjultime = ymds2js_interface (iwrtime%year, iwrtime%month, iwrtime%day, iwrtime%sec, dates30)

      IF (ijultime /= wrjultime) THEN
         
         ! Number of written timesteps
         nsteps = nsteps + 1
         IF (forward) THEN  
            DO i = 1,dnparts_max
               IF (julsec(i) <= ijultime) THEN
                  calc(i) = .TRUE.
               ELSE
                  calc(i) = .FALSE.
               ENDIF
            ENDDO
         ELSE
            DO i = 1,dnparts_max
               IF (julsec(i) >= ijultime) THEN
                  calc(i) = .TRUE.
               ELSE
                  calc(i) = .FALSE.
               ENDIF
            ENDDO
         ENDIF
         
         IF (linterpolate) THEN
            CALL get_positions_and_params (calc, itime, ifutime,  &
                 plats, plons, plevs, param, lcoupled)          
         END IF
      ENDIF
     
   !--------------------------------------------------------------------
   ! If trajectories end on different times: write positions on endtime
   !--------------------------------------------------------------------
   ELSEIF (trajtype==3) THEN
!!$      CALL write_endpos
      write(*,*) 'trajtype==3 not yet implemented for MESSY'
   ENDIF

END SUBROUTINE TIMSTP
!*******************************************************************
!
!*******************************************************************
! NEXTPOS integrates the trajectory-equation (by use of Runge-Kutta procedure)
! The results are stored in plons, plats and plevs
! udt, vdt, wdt     : wind field for a considered data time point itime
! ufut, vfut, wfut  : wind field for the next data time point ifutime
! nsteps            : time step counter for which results are stored at time iwrtime
! dt3               : time step (as real number) 
! idtsec3           : time step (as integer number)
! If itime=ifutime than a new data set is loaded 
!*******************************************************************
SUBROUTINE NEXTPOS (status, nsteps, itime, ifutime, iwrtime, dt3, idtsec3, calc, &
                    plats, plons, plevs, julsec, param, lcoupled)

  USE messy_clams_global,         ONLY: prec, dp, mdi, eps, timetype, &
                                        nx, ny, nz, &
                                        udt, vdt, wdt, ufut, vfut, wfut, &
                                        leveldt, levelfut, &
                                        dlevdzdt, & !dlevdzfut, &
                                        longrid_rad, latgrid_rad, &
                                        levelgrid, loglev, &
                                        paramtype
  USE messy_clamstraj_global,     ONLY: dt, idmsec, idmday, &
                                        trajtype, forward, linterpolate, &
                                        pslat, startlasttraj, startfirsttraj
  USE messy_clams_tools_dateconv, ONLY: ymds2js_interface, datsec, incdat_interface
  USE messy_clamstraj_calc3d,     ONLY: rungek3d
  USE messy_clamstraj_data,       ONLY: get_positions_and_params
  USE messy_clams_global,         ONLY: dnparts_max, rank, dates30
  USE messy_clams_global,         ONLY: delta_time


  IMPLICIT NONE

  REAL(PREC) :: plats(dnparts_max), plons(dnparts_max), plevs(dnparts_max)
  REAL(DP)   :: julsec(dnparts_max)

  TYPE(paramtype), DIMENSION(:), POINTER :: param

  REAL(PREC)     :: dt3
  INTEGER        :: status, nsteps, idtsec3
  TYPE(timetype) :: itime, ifutime, iwrtime
  LOGICAL        :: calc(dnparts_max), writepos(dnparts_max)
  REAL(PREC)     :: tsv(dnparts_max) !dummy variable
  LOGICAL        :: lcoupled

  ! local variables

  REAL(PREC),dimension(:,:,:),pointer :: u, v, w, umid, vmid, wmid   
  real(prec),dimension(:,:,:),pointer :: level, levelmid, dlevdz, dlevdzmid

  real(prec),dimension(:,:,:),pointer :: mdi_feld
  real(prec),dimension(:,:,:),pointer :: hilfsfeld

  logical,dimension(:,:,:),pointer :: selection 


  REAL(PREC)      :: pi, er 
  REAL(DP)        :: ijulsec, iwrjulsec, ijulsecold, &
                     past, future, calcjulsec

  INTEGER :: i, j, k

  status = 0 ! no error 

  if(rank==0) write (*,'(A,I5,2I3,I6)') 'Current timestep (yyyy mm dd ss): ', itime

  if (abs(dt3)>abs(dt)) then
     write (*,*)
     write (*,*) 'timestep greater than ',dt,' sec !!!'
     write (*,*) 'dt=',abs(dt3)
     write (*,*)
     status = -11
     return
  endif

  allocate (u(nx,ny,nz))
  allocate (v(nx,ny,nz))
  allocate (w(nx,ny,nz))
  allocate (level(nx,ny,nz))
  allocate (dlevdz(nx,ny,nz))
  allocate (umid(nx,ny,nz))
  allocate (vmid(nx,ny,nz))
  allocate (wmid(nx,ny,nz))
  allocate (levelmid(nx,ny,nz))
  allocate (dlevdzmid(nx,ny,nz))
  allocate (mdi_feld(nx,ny,nz))
  allocate (hilfsfeld(nx,ny,nz))
  allocate (selection(nx,ny,nz))
  
  mdi_feld = mdi
  selection = .true.

  er=6.371e6
  pi=4.*atan(1.)                                                        


  if (.not. lcoupled) then
     IF (datsec(itime%sec,itime%day,itime%month,itime%year) ==  &
          datsec(ifutime%sec,ifutime%day,ifutime%month,ifutime%year)) THEN  
        WRITE(*,*) 'MESSy timestep must be lower than the interval between two &
             &windfiles! '
        status = -12
        return
     ENDIF
  endif
         
  ! the function datsec returns the date in seconds from 0,1,1,1960               
  past   = datsec(itime%sec,itime%day,itime%month,itime%year)
  if (.not. lcoupled) then
     future = datsec(ifutime%sec,ifutime%day,ifutime%month,ifutime%year)
  else
     future = past + delta_time
  endif

!!!!!
  if (rank==0) then
     write (*,*) 'in NEXTPOS:'
     write (*,*) 'current time:',itime%year,itime%month,itime%day,itime%sec
     write (*,*) 'future time:',ifutime%year,ifutime%month,ifutime%day,ifutime%sec
  endif

  !-----------------------------------------------------------------
  ! Calculate:
  ! u,v,w            winds at time t       
  ! umid,vmid,wmid   winds at time t+dt3/2       
  ! udt,vdt,wdt      winds at time t+dt3       
  !-----------------------------------------------------------------

!   where(ABS((udt-mdi)/mdi)<=eps)
!      selection = .false.
!   end where
!   where(ABS((vdt-mdi)/mdi)<=eps)
!      selection = .false.
!   end where
!   where(ABS((wdt-mdi)/mdi)<=eps)
!      selection = .false.
!   end where

!   where(ABS((ufut-mdi)/mdi)<=eps)
!      selection = .false.
!   end where
!   where(ABS((vfut-mdi)/mdi)<=eps)
!      selection = .false.
!   end where
!   where(ABS((wfut-mdi)/mdi)<=eps)
!      selection = .false.
!   end where

!   u = merge(udt,mdi_feld,selection)
!   v = merge(vdt,mdi_feld,selection)
!   w = merge(wdt,mdi_feld,selection)
!   level = merge(leveldt,mdi_feld,selection)

! #ifdef MBM_CLAMS
!   hilfsfeld = (0.5*dt3*ufut+(future-past-0.5*dt3)*u)/(future-past)
!   umid = merge(hilfsfeld,mdi_feld,selection)
!   hilfsfeld = (0.5*dt3*vfut+(future-past-0.5*dt3)*v)/(future-past)
!   vmid = merge(hilfsfeld,mdi_feld,selection)
!   hilfsfeld = (0.5*dt3*wfut+(future-past-0.5*dt3)*w)/(future-past)
!   wmid = merge(hilfsfeld,mdi_feld,selection)
!   hilfsfeld = (0.5*dt3*levelfut+(future-past-0.5*dt3)*level)/(future-past)
!   levelmid = merge(hilfsfeld,mdi_feld,selection)
!   hilfsfeld = (dt3*ufut+(future-past-dt3)*u)/(future-past)
!   udt = merge(hilfsfeld,mdi_feld,selection)
!   hilfsfeld = (dt3*vfut+(future-past-dt3)*v)/(future-past)
!   vdt = merge(hilfsfeld,mdi_feld,selection)
!   hilfsfeld = (dt3*wfut+(future-past-dt3)*w)/(future-past)
!   wdt = merge(hilfsfeld,mdi_feld,selection)
!   hilfsfeld = (dt3*levelfut+(future-past-dt3)*level)/(future-past)
!   leveldt = merge(hilfsfeld,mdi_feld,selection)
! #endif
! #ifdef ECHAM5
!  umid = u
!  vmid = v
!  wmid = w
!  levelmid = level
!  udt = u
!  vdt = v
!  wdt = w
!  leveldt = level
! #endif

  DO k=1,nz
     DO j=1,ny      
        DO i=1,nx 
           IF (ABS((udt(i,j,k)-mdi)/mdi)<=eps .OR. ABS((vdt(i,j,k)-mdi)/mdi)<=eps &
                .OR. ABS((wdt(i,j,k)-mdi)/mdi)<=eps) THEN
              u(i,j,k) = MDI
              v(i,j,k) = MDI
              w(i,j,k) = MDI
              level(i,j,k) = MDI
              dlevdz(i,j,k) = MDI
              umid(i,j,k) = MDI 
              vmid(i,j,k) = MDI
              wmid(i,j,k) = MDI
              levelmid(i,j,k) = MDI
              dlevdzmid(i,j,k) = MDI
              udt(i,j,k)  = MDI 
              vdt(i,j,k)  = MDI 
              wdt(i,j,k)  = MDI 
              leveldt(i,j,k) = MDI
              dlevdzdt(i,j,k) = MDI
           ELSEIF (ABS((ufut(i,j,k)-mdi)/mdi)<=eps .OR. ABS((vfut(i,j,k)-mdi)/mdi)<=eps &
                .OR. ABS((wfut(i,j,k)-mdi)/mdi)<=eps) THEN
              u(i,j,k) = MDI
              v(i,j,k) = MDI
              w(i,j,k) = MDI
              level(i,j,k) = MDI
              dlevdz(i,j,k) = MDI
              umid(i,j,k) = MDI 
              vmid(i,j,k) = MDI
              wmid(i,j,k) = MDI
              levelmid(i,j,k) = MDI
              dlevdzmid(i,j,k) = MDI
              udt(i,j,k)  = MDI 
              vdt(i,j,k)  = MDI 
              wdt(i,j,k)  = MDI 
              leveldt(i,j,k) = MDI
              dlevdzdt(i,j,k) = MDI
           ELSE
              u(i,j,k) = udt(i,j,k)
              v(i,j,k) = vdt(i,j,k)
              w(i,j,k) = wdt(i,j,k)
              level(i,j,k) = leveldt(i,j,k)
!!$              dlevdz(i,j,k) = dlevdzdt(i,j,k)
              dlevdz(i,j,k) = -1.
              if (.not. lcoupled) then
                 umid(i,j,k) = &
                   (0.5*dt3*ufut(i,j,k)+(future-past-0.5*dt3)*u(i,j,k))/(future-past)
                 vmid(i,j,k) = & 
                   (0.5*dt3*vfut(i,j,k)+(future-past-0.5*dt3)*v(i,j,k))/(future-past)
                 wmid(i,j,k) = &
                   (0.5*dt3*wfut(i,j,k)+(future-past-0.5*dt3)*w(i,j,k))/(future-past)
                 levelmid(i,j,k) = &
                   (0.5*dt3*levelfut(i,j,k)+(future-past-0.5*dt3)*level(i,j,k))/(future-past)
!!$              dlevdzmid(i,j,k) = &
!!$                   (0.5*dt3*dlevdzfut(i,j,k)+(future-past-0.5*dt3)*dlevdz(i,j,k))/(future-past)
                 dlevdzmid(i,j,k) = -1.
                 udt(i,j,k) = &
                   (    dt3*ufut(i,j,k)+(future-past-    dt3)*u(i,j,k))/(future-past)
                 vdt(i,j,k) = &
                   (    dt3*vfut(i,j,k)+(future-past-    dt3)*v(i,j,k))/(future-past)
                 wdt(i,j,k) = &
                   (    dt3*wfut(i,j,k)+(future-past-    dt3)*w(i,j,k))/(future-past)
                 leveldt(i,j,k) = &
                   (    dt3*levelfut(i,j,k)+(future-past-dt3)*level(i,j,k))/(future-past)
!!$              dlevdzdt(i,j,k) = &
!!$                   (    dt3*dlevdzfut(i,j,k)+(future-past-dt3)*dlevdz(i,j,k))/(future-past)
                 dlevdzdt(i,j,k) = -1.
              endif
           ENDIF
        ENDDO
     ENDDO
  ENDDO

  if (lcoupled) then
     umid = u
     vmid = v
     wmid = w
     levelmid = level
     ! dlevdzmid = dlevdz
     dlevdzmid = -1.
     udt = u
     vdt = v
     wdt = w
     leveldt = level
     ! dlevdzdt = dlevdz
     dlevdzdt = -1.
  endif


  ijulsecold = ymds2js_interface (itime%year, itime%month, itime%day, itime%sec, dates30)

  ! increment the current date by the time step
  CALL incdat_interface (itime%sec,itime%day,itime%month,itime%year,idtsec3,0,0,0, dates30)     

  ijulsec = ymds2js_interface (itime%year, itime%month, itime%day, itime%sec, dates30)

 
  IF (trajtype == 1) THEN  ! all trajectories start on the same time 
     IF (forward) THEN    
      IF (julsec(1) < ijulsec) THEN
           calc = .TRUE.
        ELSE
           calc = .FALSE.
        ENDIF
     
     ELSE
        IF (julsec(1) > ijulsec) THEN
           calc = .TRUE.
        ELSE
           calc = .FALSE.
        ENDIF
     ENDIF
    
  ELSE  ! trajectories start on different times

!!$     IF (forward) THEN     
!!$!        DO i = 1,nparts
!!$        DO i = loop_start, loop_end             
!!$           IF (julsec(i) < ijulsec) THEN 
!!$              IF ((.NOT. calc(i)) .AND. (ijulsec <= endjulsec(i))) THEN
!!$                  CALL get_start_end_values (startvalues, julsec(i), ifutime, i&
!!$                       &, plats, plons, plevs)
!!$               ENDIF
!!$              calc(i) = .TRUE.
!!$           ELSE
!!$              calc(i) = .FALSE.
!!$           ENDIF
!!$
!!$           IF (ijulsec==startlasttraj .and. julsec(i)==ijulsec) THEN
!!$              CALL get_start_end_values (startvalues, julsec(i), ifutime, i,&
!!$                   & plats, plons, plevs)
!!$           ENDIF
!!$        ENDDO
!!$
!!$     ELSE
!!$!        DO i = 1,nparts
!!$        DO i = loop_start, loop_end
!!$           IF (julsec(i) > ijulsec) THEN
!!$              IF ((.NOT. calc(i)) .AND. (ijulsec >= endjulsec(i))) THEN 
!!$                 CALL get_start_end_values (startvalues, julsec(i), ifutime, i,&
!!$                      & plats, plons, plevs)
!!$              ENDIF
!!$              calc(i) = .TRUE.
!!$           ELSE
!!$              calc(i) = .FALSE.
!!$           ENDIF
!!$
!!$           IF (ijulsec==startfirsttraj .and. julsec(i)==ijulsec) THEN
!!$              CALL get_start_end_values (startvalues, julsec(i), ifutime, i,&
!!$                   & plats, plons, plevs)
!!$           ENDIF
!!$        ENDDO
!!$     ENDIF
!!$

  ENDIF

 

!!$  IF (trajtype==3) THEN
!!$
!!$     IF (forward) THEN
!!$        DO i = 1,nparts_max
!!$           IF (endjulsec(i)<ijulsec) calc(i) = .FALSE.
!!$        ENDDO
!!$     ELSE
!!$        DO i = 1,nparts_max
!!$           IF (endjulsec(i)>ijulsec) calc(i) = .FALSE.
!!$        ENDDO
!!$     ENDIF
!!$
!!$  ENDIF

  !-------------------------------------------------------------
  ! Runge kutta scheme to advect particles in three dimensions
  !-------------------------------------------------------------

  CALL rungek3d(plons,plats,plevs,calc,tsv, &
                u,   v,   w,  level, dlevdz,  &
                umid,vmid,wmid,levelmid, dlevdzmid,   &
                udt, vdt, wdt,leveldt, dlevdzdt,   &    
                nx,  ny,  nz,    &
                longrid_rad, latgrid_rad, levelgrid,   &
                dnparts_max,1,dnparts_max,dt3,pslat,loglev)                              

  !-------------------------------------------------------------
  ! Write particle positions
  !-------------------------------------------------------------
  ijulsec = ymds2js_interface (itime%year, itime%month, itime%day, itime%sec, dates30)
  iwrjulsec = ymds2js_interface (iwrtime%year, iwrtime%month, iwrtime%day, iwrtime%sec, dates30)

  IF (forward) THEN

     IF (ijulsec==startfirsttraj .AND. trajtype==1) THEN

        nsteps = nsteps + 1  ! number of written timesteps
        calc=.TRUE.        
        IF (linterpolate) THEN
           CALL get_positions_and_params (calc, itime, ifutime,  &
                plats, plons, plevs, param, lcoupled)          
        END IF
     ELSEIF (ijulsec==iwrjulsec .AND. ijulsec>startfirsttraj) THEN

        writepos = .FALSE.
        WHERE (julsec <= ijulsec) writepos=.TRUE.
        nsteps = nsteps + 1  ! number of written timesteps
        IF (linterpolate) THEN
           CALL get_positions_and_params (writepos, itime, ifutime,  &
                plats, plons, plevs, param, lcoupled)          
        END IF
     ENDIF
     
  ELSE
     
     IF(ijulsec==startlasttraj .AND. trajtype==1) THEN

        nsteps = nsteps + 1  ! number of written timesteps
        calc = .TRUE.
        IF (linterpolate) THEN
           CALL get_positions_and_params (calc, itime, ifutime,  &
                plats, plons, plevs, param, lcoupled)          
        END IF
     ELSEIF (ijulsec==iwrjulsec .AND. ijulsec<startlasttraj) THEN

        writepos = .FALSE.
        WHERE (julsec >= ijulsec) writepos=.TRUE.
        nsteps = nsteps + 1  ! number of written timesteps
        IF (linterpolate) THEN
           CALL get_positions_and_params (writepos, itime, ifutime,  &
                plats, plons, plevs, param, lcoupled)          
        END IF
     ENDIF

  ENDIF

!!$  ! Write positions on endtime
!!$  IF (trajtype==3) THEN
!!$
!!$     IF (forward) THEN
!!$!        DO i = 1,nparts
!!$        DO i = loop_start, loop_end        
!!$           IF ((endjulsec(i)>=ijulsecold) .AND. (endjulsec(i)<=ijulsec)) THEN
!!$              calcjulsec = endjulsec(i)
!!$              CALL get_start_end_values (endvalues, calcjulsec, ifutime, i, plats, plons, plevs)
!!$           ENDIF
!!$        ENDDO
!!$     ELSE
!!$!        DO i = 1,nparts
!!$        DO i = loop_start, loop_end
!!$           IF ((endjulsec(i)<=ijulsecold) .AND. (endjulsec(i)>=ijulsec)) THEN
!!$              calcjulsec = endjulsec(i)
!!$              CALL get_start_end_values (endvalues, calcjulsec, ifutime, i, plats, plons, plevs)
!!$           ENDIF
!!$        ENDDO
!!$     ENDIF
!!$    
!!$  ENDIF
  
  IF (ijulsec==iwrjulsec) THEN
     ! increment the writing date and time by the writing frequency             
     IF (linterpolate) THEN
        CALL incdat_interface (iwrtime%sec,iwrtime%day,iwrtime%month,iwrtime%year, &
             idmsec,idmday,0,0, dates30)          
     END IF
  ENDIF

  deallocate (u, v, w, level,dlevdz)
  deallocate (umid, vmid, wmid, levelmid,dlevdzmid)
  deallocate (mdi_feld)
  deallocate (hilfsfeld)
  deallocate (selection)

END SUBROUTINE NEXTPOS

END MODULE messy_clamstraj_timestep

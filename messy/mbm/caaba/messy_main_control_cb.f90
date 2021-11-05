! Time-stamp: <2018-05-15 21:08:32 sander>

! This module controls the MESSy submodel CALLs

! Authors:
! Rolf Sander, MPICH, Mainz, 2003-2007
! Hella Riede, MPICH, Mainz, 2007

! This program is free software; you can redistribute it and/or
! modify it under the terms of the GNU General Public License
! as published by the Free Software Foundation; either version 2
! of the License, or (at your option) any later version.
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
! You should have received a copy of the GNU General Public License
! along with this program; if not, get it from:
! http://www.gnu.org/copyleft/gpl.html

!*****************************************************************************

MODULE messy_main_control_cb

  ! This is the base model interface layer (BMIL), i.e. the interface between
  ! the main box model program and the MESSy submodels

  USE caaba_mem, ONLY: USE_CLOUDJ, USE_DISSOC, USE_JVAL, USE_MECCA, USE_RADJIMT, USE_READJ, &
#ifdef E4CHEM
                       USE_E4CHEM,                                 &
#endif
                       USE_SAPPHO, USE_SEMIDEP, USE_TRAJECT

  IMPLICIT NONE

CONTAINS

  !***************************************************************************

  SUBROUTINE messy_init

    USE caaba_mem,                ONLY: model_time, model_end, runtime, &
                                        time_string, time0_jul, &
                                        t0year, t0month, t0day, t0hour, &
                                        t0min, t0sec, firstjan_jul, &
                                        l_steady_state_stop
    USE messy_main_constants_mem, ONLY: HLINE1, OneDay, DP
    USE messy_main_timer,         ONLY: gregor2julian, eval_time_str
    USE messy_cloudj_box,         ONLY:  cloudj_init
    USE messy_dissoc_box,         ONLY:  dissoc_init
    USE messy_jval_box,           ONLY:    jval_init
    USE messy_mecca_box,          ONLY:   mecca_init
!!#D radjimt +
    USE messy_radjimt_box,        ONLY: radjimt_init
!!#D radjimt -
    USE messy_readj_box,          ONLY:   readj_init
    USE messy_sappho_box,         ONLY:  sappho_init
    USE messy_traject_box,        ONLY: traject_init
#ifdef E4CHEM
    USE messy_e4chem_box,         ONLY:  e4chem_init
#endif

    IMPLICIT NONE
    INTRINSIC :: ABS, TRIM
    INTEGER :: status

    !-------------------------------------------------------------------------

    ! special requirements:
    ! - traject_init must be called first to have physical boundary conditions
    !   and thus cair consistent when initializing specs in routine x0
    WRITE(*,*) HLINE1
    IF (USE_TRAJECT) CALL traject_init ! must be called before mecca_init
    IF (USE_MECCA)   CALL   mecca_init
!!#D radjimt +
    IF (USE_RADJIMT) CALL radjimt_init
!!#D radjimt -
    IF (USE_READJ)   CALL   readj_init
    IF (USE_SAPPHO)  CALL  sappho_init
    IF (USE_CLOUDJ)  CALL  cloudj_init
    IF (USE_DISSOC)  CALL  dissoc_init
    IF (USE_JVAL)    CALL    jval_init
#ifdef E4CHEM
    IF (USE_E4CHEM)  CALL  e4chem_init
#endif

    !-------------------------------------------------------------------------

    WRITE(*,*) HLINE1
    WRITE(*,*) 'Input/Output time unit and origin: ', TRIM(time_string)
    WRITE(*,'(1X, "Start day  = ", F13.8, "  (=", F18.8, " s)")') &
      model_time/OneDay, model_time
    IF (l_steady_state_stop) THEN       
      WRITE(*,*) 'runtime    = until steady state reached'
    ELSE                                
      WRITE(*,'(1X, "End day    = ", F13.8, "  (=", F18.8, " s)")') &
        model_end/OneDay, model_end
      WRITE(*,'(1X, "runtime    = ", F13.8, "  (=", F18.8, " s)")') &
        runtime, runtime * OneDay 
    ENDIF
    
    IF (.NOT. USE_TRAJECT) THEN
      !> determine start time year, month, day, hour, minute, second
      ! time0_jul used as dummy argument, time factor output not needed here
      ! if in trajectory mode, the following statements have already been done
      !   messy_traject_box:get_model_start_end
      CALL eval_time_str(status, time_string, time0_jul, &
                         t0year, t0month, t0day, t0hour, t0min, t0sec) 
      IF (status/=0) STOP 1
      !print *, 't0year  = ', t0year
      !print *, 't0month = ', t0month
      !print *, 't0day   = ', t0day
      !print *, 't0hour  = ', t0hour
      !print *, 't0min   = ', t0min
      !print *, 't0sec   = ', t0sec

      ! determine Julian date of time origin
      time0_jul = gregor2julian(t0year, t0month, t0day, t0hour, t0min, t0sec)
      !print *, 'time0_jul = ', time0_jul

      ! determine Julian date of 01-JAN of start year (t0year)
      firstjan_jul = gregor2julian(t0year, 1, 1, 0, 0, 0)
    ENDIF

  END SUBROUTINE messy_init

  !***************************************************************************

  SUBROUTINE messy_physc

    USE messy_cloudj_box,  ONLY:  cloudj_physc
    USE messy_dissoc_box,  ONLY:  dissoc_physc
    USE messy_jval_box,    ONLY:    jval_physc
    USE messy_mecca_box,   ONLY:   mecca_physc
!!#D radjimt +
    USE messy_radjimt_box, ONLY: radjimt_physc
!!#D radjimt -
    USE messy_sappho_box,  ONLY:  sappho_physc
    USE messy_semidep_box, ONLY: semidep_physc
    USE messy_traject_box, ONLY: traject_physc
#ifdef E4CHEM
    USE messy_e4chem_box,  ONLY:  e4chem_physc
#endif

    ! Do not change the order of *_physc subroutines!
    IF (USE_TRAJECT) CALL traject_physc
    IF (USE_SEMIDEP) CALL semidep_physc
    IF (USE_CLOUDJ)  CALL  cloudj_physc
    IF (USE_DISSOC)  CALL  dissoc_physc
    IF (USE_JVAL)    CALL    jval_physc
!!#D radjimt +
    IF (USE_RADJIMT) CALL radjimt_physc
!!#D radjimt -
    IF (USE_SAPPHO)  CALL  sappho_physc
    IF (USE_MECCA)   CALL   mecca_physc
#ifdef E4CHEM
    IF (USE_E4CHEM)  CALL  e4chem_physc
#endif

  END SUBROUTINE messy_physc

  !***************************************************************************

  SUBROUTINE messy_result

    USE messy_cloudj_box,  ONLY:  cloudj_result
    USE messy_dissoc_box,  ONLY:  dissoc_result
    USE messy_jval_box,    ONLY:    jval_result
    USE messy_mecca_box,   ONLY:   mecca_result
!!#D radjimt +
    USE messy_radjimt_box, ONLY: radjimt_result
!!#D radjimt -
    USE messy_sappho_box,  ONLY:  sappho_result
    USE messy_traject_box, ONLY: traject_result
#ifdef E4CHEM
    USE messy_e4chem_box,  ONLY:  e4chem_result
#endif

    IF (USE_MECCA)   CALL   mecca_result
    IF (USE_CLOUDJ)  CALL  cloudj_result
    IF (USE_DISSOC)  CALL  dissoc_result
    IF (USE_JVAL)    CALL    jval_result
!!#D radjimt +
    IF (USE_RADJIMT) CALL radjimt_result
!!#D radjimt -
    IF (USE_SAPPHO)  CALL  sappho_result
    IF (USE_TRAJECT) CALL traject_result
#ifdef E4CHEM
    IF (USE_E4CHEM)  CALL  e4chem_result
#endif

  END SUBROUTINE messy_result

  !***************************************************************************

  SUBROUTINE messy_finish

    USE messy_cloudj_box,  ONLY:  cloudj_finish
    USE messy_dissoc_box,  ONLY:  dissoc_finish
    USE messy_jval_box,    ONLY:    jval_finish
    USE messy_mecca_box,   ONLY:   mecca_finish
!!#D radjimt +
    USE messy_radjimt_box, ONLY: radjimt_finish
!!#D radjimt -
    USE messy_sappho_box,  ONLY:  sappho_finish
    USE messy_traject_box, ONLY: traject_finish
#ifdef E4CHEM
    USE messy_e4chem_box,  ONLY:  e4chem_finish
#endif

    IF (USE_CLOUDJ)  CALL  cloudj_finish
    IF (USE_DISSOC)  CALL  dissoc_finish
    IF (USE_JVAL)    CALL    jval_finish
    IF (USE_MECCA)   CALL   mecca_finish
!!#D radjimt +
    IF (USE_RADJIMT) CALL radjimt_finish
!!#D radjimt -
    IF (USE_SAPPHO)  CALL  sappho_finish
    IF (USE_TRAJECT) CALL traject_finish
#ifdef E4CHEM
    IF (USE_E4CHEM)  CALL  e4chem_finish
#endif

  END SUBROUTINE messy_finish

  !***************************************************************************

END MODULE messy_main_control_cb

!*****************************************************************************

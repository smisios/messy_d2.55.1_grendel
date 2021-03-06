#include "messy_main_ppd_bi.inc"
!
!*****************************************************************************
! Authors:
! Benedikt Steil, MPICH,      2015
! Patrick Joeckel, DLR, April 2015 : mtskip flag for tracer
!
! submodel maintainer: Benedikt Steil
!
! Currently only skipping for MECCA is implemented. 
! Please, contact maintainer for skipping of further submodels.
!*****************************************************************************

MODULE messy_mtskip_si

!call tree 
! scan1.i90 
!            407: CALL messy_global_start    <--- messy/echam5/bmil/messy_main_control_e5.f90
!                  .. advection..
!            617: CALL messy_local_start     <--- messy/echam5/bmil/messy_main_control_e5.f90
!            625: CALL gpc     
!                    |------> CALL physc     <--- 79: gpc.i90 
!                                |------> CALL 1018 messy_physc.i90  <--- messy/echam5/bmil/messy_main_control_e5.f90
!            633: CALL messy_local_end       <--- /messy/echam5/bmil/messy_main_control_e5.f90 
!            666: CALL messy_global_end      <--- /messy/echam5/bmil/messy_main_control_e5.f90 

  USE messy_main_blather_bi,    ONLY: start_message_bi, end_message_bi,    &
                                      error_bi

  USE messy_mtskip

#ifdef MESSYTENDENCY
  !tendency budget
  USE messy_main_tendency_bi,    ONLY: mtend_get_handle,       &    
                                       mtend_get_start_l,      &    
                                       mtend_add_l,            &    
                                       mtend_register,         &    
                                       mtend_id_tracer
#endif

  IMPLICIT NONE
  PRIVATE
  SAVE
  INTRINSIC :: NULL

  !input-namelist cpl
  REAL(DP), PUBLIC                  :: dummy

  !logical switches
  INTEGER                           :: nstep_always_mtskip = 10
  LOGICAL                           :: lmodules_mtskip = .FALSE. ! 

  !channel-objects
  !logical switches
  REAL(DP), POINTER                 :: rmode  => NULL() !0 never skip, always newly calculate the process
                                                        !1 skipping active, newly calculate the process
                                                        !2 skipping active, only add tendency
  REAL(DP), POINTER                 :: rmecca => NULL()
  REAL(DP), POINTER                 :: rgmxe  => NULL()
  !REAL(DP), POINTER                 :: rmsbm => NULL()
  
  !time-control
  REAL(DP), POINTER                 :: tsl       => NULL()  &
                                      ,tslphys   => NULL()
  REAL(DP)                          :: tsl_init

! op_pj_20150707+
!!$  REAL(DP), POINTER                 :: cdissephtl => NULL()  ! distance sun - earth [AU]
!!$  REAL(DP), POINTER                 :: decphtl    => NULL()  ! declination of sun
!!$  REAL(DP), POINTER                 :: raphtl     => NULL()  ! right ascension of sun
!!$  ! cos(solar zenith angle)
!!$  REAL(DP), DIMENSION(:,:), POINTER :: cosszaphtl => NULL()
!!$  ! relative day length
!!$  REAL(DP), DIMENSION(:,:), POINTER :: rdaylphtl => NULL()
!!$  !
!!$  ! cos(solar zenith angle); cut-off at zero
!!$  REAL(DP), DIMENSION(:,:), POINTER :: cosszacphtl => NULL()
! op_pj_20150707-

  !Additional memory for prognostic variables, which a subject to skipping
  TYPE, PUBLIC :: mtskip3d
     REAL(DP), POINTER :: ptr(:,:,:) => NULL()
  END TYPE mtskip3d
  TYPE(mtskip3d), PUBLIC, POINTER           :: mtskiptend(:)
  REAL(DP), DIMENSION(:,:,:), POINTER       :: mtskiptend_q,    &
                                               mtskiptend_xl,   &
                                               mtskiptend_xi
  TYPE(mtskip3d), PUBLIC, POINTER           :: mtskiptend_budneg(:)
  REAL(DP), DIMENSION(:,:,:), POINTER       :: mtskiptend_q_budneg,    &
                                               mtskiptend_xl_budneg,   &
                                               mtskiptend_xi_budneg
  !and it's administration
  LOGICAL, DIMENSION(:), POINTER            :: mapmem_tend => NULL()

  ! OFFSET [s] for orbit calculations      ! op_pj_20150707
  REAL(DP), POINTER :: dt_offset => NULL() ! op_pj_20150707

#ifdef MESSYTENDENCY
  INTEGER :: my_handle
#endif

  NAMELIST /CPL/ dummy

  PUBLIC :: mtskip_initialize
  PUBLIC :: mtskip_init_memory 
  PUBLIC :: mtskip_init_coupling
  PUBLIC :: mtskip_local_start
  PUBLIC :: mtskip_global_start
  PUBLIC :: mtskip_free_memory
  
  CONTAINS

    ! **************************************************************************
    ! Public routines
    ! **************************************************************************

    ! **************************************************************************

    SUBROUTINE mtskip_initialize

      ! ECHAM5/MESSy
! op_pj_20150707: finish only available in ECHAM5, replaced by error_bi ...
      USE messy_main_mpi_bi,    ONLY: p_parallel_io, p_io, p_bcast!!$, finish
      ! MESSY
      USE messy_main_tools,     ONLY: find_next_free_unit
      USE messy_main_timer,     ONLY: time_step_len, delta_time,         &
                                      nstep => current_time_step

      IMPLICIT NONE

      ! LOCAL
      CHARACTER(LEN=*), PARAMETER :: substr = 'mtskip_initialize'
      INTEGER                     :: iou            ! I/O unit
      INTEGER                     :: status         ! error status
!!$   LOGICAL                     :: lex = .FALSE.  ! file exists ?
!!$   INTEGER                     :: rst_fsize      ! size of mtskip restart file

      CALL start_message_bi(modstr,'INITIALISATION',modstr)

      IF (p_parallel_io) THEN 
        iou = find_next_free_unit(100,200)
        CALL mtskip_read_nml_ctrl(status, iou) 
        IF (status /= 0) CALL error_bi('error in mtskip_read_nml_ctrl',substr)
      ENDIF
      ! BROADCAST RESULTS
      CALL p_bcast(phtl_offset, p_io)
      CALL p_bcast(nstep_always, p_io)
      CALL p_bcast(trigmtskip, p_io)
      CALL p_bcast(addon_skip, p_io)
      CALL p_bcast(lfixneg, p_io)
      CALL p_bcast(lbudneg, p_io)
      CALL p_bcast(lmecca_mtskip, p_io)
      CALL p_bcast(lgmxe_mtskip, p_io)
      !CALL p_bcast(lmsbm_mtskip, p_io)

      !IF (p_parallel_io) THEN
      !  iou = find_next_free_unit(100,200)
      !  CALL mtskip_read_nml_cpl(status, iou)
      !  IF (status /= 0) CALL error_bi('error in mtskip_read_nml_cpl',substr)
      !ENDIF
      ! BROADCAST RESULTS

      !Go thru basic logic

      !(a) without or with correct skipping 
      IF (trigmtskip <= 1 ) THEN
        IF (p_parallel_io) THEN
          WRITE(*,*) '* ',substr, ': trig_mtskip<=1 stops the model.                       *' 
          WRITE(*,*) '* ',substr, ': The internal default of trig_mtskip is 1.             *'
          WRITE(*,*) '* ',substr, ': trig_mtskip=1 means you intend to calculate processes *'
          WRITE(*,*) '* ',substr, ': every time step, i.e. you do not need mtskip.         *'
          WRITE(*,*) '* ',substr, ': Please switch off mtskip in switch.nml and change     *'
! op_pj_20150707+
!!$          WRITE(*,*) '* ',substr, ': channel for jval_cossza and jval_cdisse in jval.nml   *'
!!$          WRITE(*,*) '* ',substr, ': back to orbit and restart the model.                  *'
          WRITE(*,*) '* ',substr, ': c_offset in orbit.nml and restart the model           *'
! op_pj_20150707-
        ENDIF
!!$        CALL finish(substr, 'cpl namelist: steps >= 1 and steps odd is skipping; steps = 1 means no skipping;')
        CALL error_bi('ctrl namelist: steps >= 1 and steps odd is skipping; steps = 1 means no skipping;',substr)
      ENDIF
      IF ( MOD(trigmtskip,2) == 0 ) &
!!$        CALL finish(substr, 'cpl namelist: only odd multiples of steps are permitted')
        !mz_bs_20151118+
        CALL error_bi('ctrl namelist: only odd multiples of steps are permitted', substr)
        !mz_bs_20151118-

! op_pj_20150708: this check should be always performed here (only once)
      IF ( nstep_always < 0 ) &
           CALL error_bi('ctrl namelist: negative nstep_always are not permitted',substr)

      !(b) consider addon_skip
      IF (addon_skip .LT. 0._dp) THEN
        IF ( trigmtskip .EQ. 3 ) addon_skip=2._dp
        IF ( trigmtskip .GT. 3 ) addon_skip=3._dp
      ENDIF
      !(c)
      !calculate the new timestep-length in case of skipping.
      tsl_init = (addon_skip+REAL(trigmtskip,dp)) * &
                 delta_time

! op_pj_20150708+
! see comment in lresume-block below
!!$      !(d) consider nstep_always 
!!$      IF ( lstart ) THEN
!!$        !new start
!!$        !consider nstep_always
!!$        IF ( nstep_always .LT. 0 ) &
!!$!!!$          CALL finish(substr, 'cpl namelist: negative nstep_always are not permitted')
!!$          CALL error_bi('cpl namelist: negative nstep_always are not permitted',substr)
!!$        nstep_always_mtskip = (nstep_always/trigmtskip)*trigmtskip  &
!!$                              + trigmtskip-1
!!$      ENDIF
! op_pj_20150708-

! op_pj_20150708+: The restart file handling is not required here, since 
!                  all can be controlled via the namelists in channel.nml:
!                  (1) If the "ignore"-flag is set for all MTSKIP channel 
!                      objects (or set as default for the entire MTSKIP channel)
!                      in channel.nml, MTSKIP can be switched on after
!                      restart, although lrestreq=.TRUE.. In that case,
!                      all channel objcts contain zero after initialisation.
!                  (2) Since all MTSKIP objects (including the tendencies)
!                      have lrestreq=.TRUE. they are always dumped for restart,
!                      even though they might still be zero (from 
!                      initialisation) and have not yet been calculated,
!                      because nstep_always > nstep. The file size check is
!                      therefore obsolete (and anyway not very robust!).
!   -> Thus, you only need to consider two cases:
!      a) MTSKIP objects have been read from restart file
!      b) MTSKIP objects have NOT been read from restart file
!         (because MTSKIP has been switched on AFTER restart)
!      This check (of initialisiation status) can be done either at the entry
!      point init_tracer, or at the entry point global_start.
!      Since init_tracer is not used here, I moved it to global_start ...
!   -> For consistency (and readability), I also moved the lstart-check ...
!
!!$      IF ( lresume ) THEN
!!$        !restart from new submit
!!$        lex=.FALSE.
!!$        INQUIRE(file=TRIM('restart_'//modstr//'.nc'),EXIST=lex), &
!!$                SIZE=rst_fsize)
!!$        IF ( (.NOT. lex) .OR. (rst_fsize .LT. 2000000) ) THEN
!!$          !Consider nstep_always.
!!$          !This covers the cases that either no restart file exists or
!!$          !that the restart file contains no tendency fields
!!$          !for prognostic variables. The latter may occure when nstep_always
!!$          !is larger than restart step.
!!$          IF ( nstep_always .LT. 0 ) &
!!$!!!$            CALL finish(substr, 'cpl namelist: negative nstep_always are not permitted')
!!$            CALL error_bi('cpl namelist: negative nstep_always are not permitted', substr)
!!$          nstep_always_mtskip = ((nstep+nstep_always)/trigmtskip)*trigmtskip  &
!!$                                + trigmtskip-1
!!$        ELSE
!!$          !don't consider nstep_always
!!$          nstep_always_mtskip = 0
!!$        ENDIF
!!$     ENDIF
! op_pj_20150708-

      !(e) Check, if there are any submodels to skip?
      lmodules_mtskip = .FALSE.
      IF (lmecca_mtskip               &
          .OR. lgmxe_mtskip           &
      !   .OR. lmsbm_mtskip           &
         ) lmodules_mtskip = .TRUE.

#ifdef MESSYTENDENCY
    my_handle = mtend_get_handle(modstr)
#endif
      
      CALL end_message_bi(modstr,'INITIALISATION',modstr)

    END SUBROUTINE mtskip_initialize

    ! **************************************************************************

    SUBROUTINE mtskip_init_memory 

      ! MESSy/BMIL
      USE messy_main_mpi_bi,           ONLY: p_parallel_io
      USE messy_main_channel_error_bi, ONLY: channel_halt
      USE messy_main_channel_bi,       ONLY: SCALAR,            &
                                             GP_2D_HORIZONTAL, GP_3D_MID
      USE messy_main_tracer_mem_bi,    ONLY: ntrac_gp, ti_gp, GPTRSTR

      ! MESSy
      USE messy_main_channel,          ONLY: new_channel, new_channel_object  &
                                           , new_attribute &
                                           , new_channel_object_reference
      USE messy_main_tracer,           ONLY: get_tracer

      IMPLICIT NONE

      ! LOCAL
      CHARACTER(LEN=*), PARAMETER :: substr = 'mtskip_init_memory'
      INTEGER :: status
      INTEGER :: jt

#ifdef MESSYTENDENCY
    CALL mtend_register(my_handle, mtend_id_tracer)
#endif

      CALL start_message_bi(modstr, 'CHANNEL DEFINITION', substr)

      ! DEFINE NEW CHANNEL
      CALL new_channel(status, modstr, lrestreq=.TRUE.)
      CALL channel_halt(substr, status)

      ! OBJECTS

      CALL new_channel_object(status, modstr, 'tsl', p0=tsl, reprid=SCALAR)
      CALL channel_halt(substr, status)
      CALL new_attribute(status, modstr, 'tsl'   &   
           , 'long_name', c='mtskip time step length' )
      CALL channel_halt(substr, status)
      CALL new_attribute(status, modstr, 'tsl'   &   
           , 'units', c='sec' )
      CALL channel_halt(substr, status)
! op_pj_20150708+: ATTENTION: This will be overwritten after restart!
!                 Thus, in case MTSKIP has been switched on AFTER restart,
!                 it will be overwritten by zero!
!                 -> Moved to global_start ...
      tsl = tsl_init
! op_pj_20150708-

      CALL new_channel_object(status, modstr, 'tslphys', p0=tslphys, reprid=SCALAR)
      CALL channel_halt(substr, status)
      CALL new_attribute(status, modstr, 'tslphys'   &   
           , 'long_name', c='mtskip time step length used for integration, set in global start' )
      CALL channel_halt(substr, status)
      CALL new_attribute(status, modstr, 'tslphys'   &   
           , 'units', c='sec' )
      CALL channel_halt(substr, status)

      CALL new_channel_object(status, modstr, 'rmecca', p0=rmecca, reprid=SCALAR)
      CALL channel_halt(substr, status)
      CALL new_attribute(status, modstr, 'rmecca'   &   
           , 'long_name', c='basic switch: 1/0 skipping mecca on or off' )
      CALL channel_halt(substr, status)
! op_pj_20150708+: ATTENTION: This will be overwritten after restart!
!                 Thus, in case MTSKIP has been switched on AFTER restart,
!                 it will be overwritten by zero!
!                 -> Moved to global_start ...
      IF (lmecca_mtskip) THEN
        rmecca=1._dp
      ELSE
        rmecca=0._dp
      ENDIF
! op_pj_20150708-

      ! op_pj_20150707+
      ! Create a channel object reference named 'rjval' sharing the
      ! memory (and therefore content) with object 'rmecca'.
      ! This has been introduced to have independent switches for the
      ! submodels (JVAL needs not always to run with MECCA ...).
      ! At the moment, this guarnatees that MECCA and JVAL are 
      ! MTSKIPped synchronously. Maybe this will change later, e.g.,
      ! if JVAL is used in a different way ...
      CALL new_channel_object_reference(status, &
           modstr, 'rmecca', modstr, 'rjval', .FALSE.)
      CALL channel_halt(substr, status)
      CALL new_attribute(status, modstr, 'rjval'   &   
           , 'long_name', c='basic switch: 1/0 skipping jval on or off' )
      CALL channel_halt(substr, status)
      ! op_pj_20150707-

      CALL new_channel_object(status, modstr, 'rgmxe', p0=rgmxe, reprid=SCALAR)
      CALL channel_halt(substr, status)
      CALL new_attribute(status, modstr, 'rgmxe'   &   
           , 'long_name', c='basic switch: 1/0 skipping gmxe on or off' )
      CALL channel_halt(substr, status)
! op_pj_20150708+: ATTENTION: This will be overwritten after restart!
!                 Thus, in case MTSKIP has been switched on AFTER restart,
!                 it will be overwritten by zero!
!                 -> Moved to global_start ...
      IF (lgmxe_mtskip) THEN
        rgmxe=1._dp
      ELSE
        rgmxe=0._dp
      ENDIF
! op_pj_20150708-

      !CALL new_channel_object(status, modstr, 'rmsbm', p0=rmsbm, reprid=SCALAR)
      !CALL channel_halt(substr, status)
      !CALL new_attribute(status, modstr, 'rmsbm'   &   
      !     , 'long_name', c='basic switch: 1/0 skipping msbm on or off' )
      !CALL channel_halt(substr, status)
      !IF (lmsbm_mtskip) THEN
      !  rmsbm=1._dp
      !ELSE
      !  rmsbm=0._dp
      !ENDIF

      !.....

      CALL new_channel_object(status, modstr, 'rmode', p0=rmode, reprid=SCALAR)
      CALL channel_halt(substr, status)
      CALL new_attribute(status, modstr, 'rmode'   &   
           , 'long_name', c=                       &
'0 never skip, always newly calculate the process;1 skipping active, newly calculate the process;2 skipping active, only add tendency' )
      CALL channel_halt(substr, status)

      CALL new_channel_object(status, modstr, 'dt_offset' &
           , p0 = dt_offset, reprid=SCALAR)
      CALL channel_halt(substr, status)
      CALL new_attribute(status, modstr, 'dt_offset' &
           , 'long_name', c='offset for orbit calculation')
      CALL channel_halt(substr, status)
      CALL new_attribute(status, modstr, 'dt_offset', 'units', c='s')
      CALL channel_halt(substr, status)

      CALL end_message_bi(modstr, 'CHANNEL DEFINITION', substr)

    END SUBROUTINE mtskip_init_memory 

    ! **************************************************************************

    SUBROUTINE mtskip_init_coupling

      ! MESSy/BMIL
      USE messy_main_mpi_bi,           ONLY: p_parallel_io
      USE messy_main_channel_error_bi, ONLY: channel_halt
      USE messy_main_channel_bi,       ONLY: SCALAR,           &
                                             GP_2D_HORIZONTAL, GP_3D_MID
      USE messy_main_tracer_mem_bi,    ONLY: ntrac_gp, ti_gp, GPTRSTR
      USE messy_main_tracer_tools_bi,  ONLY: tracer_halt
      USE messy_main_tracer,           ONLY: get_tracer, ON, I_MTSKIP
      ! MESSy
      USE messy_main_channel,          ONLY: new_channel, new_channel_object &
                                           , new_attribute
      IMPLICIT NONE


      ! LOCAL
      CHARACTER(LEN=*), PARAMETER :: substr = 'mtskip_init_coupling'
      INTEGER :: status
      INTEGER :: jt, jm, id_mtskip, jt1
      INTEGER :: on_check

#ifdef MESSYTENDENCY
      CALL mtend_register(my_handle, mtend_id_tracer)
#endif

      CALL start_message_bi(modstr, 'TENDENCY-CHANNEL-OBJECT DEFINITION', substr)

      IF (lmodules_mtskip) THEN
     
        ALLOCATE(mapmem_tend(ntrac_gp))
        mapmem_tend(1:ntrac_gp) = .FALSE.
        DO jt = 1, ntrac_gp
          CALL get_tracer(status, GPTRSTR, jt, I_MTSKIP, i=on_check)
          CALL tracer_halt(substr, status)
          IF (on_check == ON ) THEN
            mapmem_tend(jt) = .TRUE.
          ENDIF
          !For H2O-Tracer, I_MTSKIP is set in Mecca
        ENDDO

! op_pj_20160511+
!!$        IF (lgmxe_mtskip) THEN
        IF (lgmxe_mtskip .OR. lmecca_mtskip) THEN
! op_pj_20160511-
          IF (p_parallel_io) WRITE(*,*) ' ... msktq'
          CALL new_channel_object(status, modstr, &
            'msktq', p3=mtskiptend_q, reprid=GP_3D_MID, lrestreq=.TRUE. )
            !'msktq', p3=mtskiptend_q%ptr, reprid=GP_3D_MID, lrestreq=.TRUE. )
          CALL channel_halt(substr, status)
          CALL new_attribute(status, modstr, &
            'msktq', 'units', c='kg/(kg*sec)')
          CALL channel_halt(substr, status)
! op_pj_20160511+
        END IF

        IF (lgmxe_mtskip) THEN
! op_pj_20160511-
          IF (p_parallel_io) WRITE(*,*) ' ... msktxl'
          CALL new_channel_object(status, modstr, &
            'msktxl', p3=mtskiptend_xl, reprid=GP_3D_MID, lrestreq=.TRUE. )
            !'msktxl', p3=mtskiptend_xl%ptr, reprid=GP_3D_MID, lrestreq=.TRUE. )
          CALL channel_halt(substr, status)
          CALL new_attribute(status, modstr, &
            'msktxl', 'units', c='kg/(kg*sec)')
          CALL channel_halt(substr, status)

          IF (p_parallel_io) WRITE(*,*) ' ... msktxi'
          CALL new_channel_object(status, modstr, &
            'msktxi', p3=mtskiptend_xi, reprid=GP_3D_MID, lrestreq=.TRUE. )
            !'msktxi', p3=mtskiptend_xi%ptr, reprid=GP_3D_MID, lrestreq=.TRUE. )
          CALL channel_halt(substr, status)
          CALL new_attribute(status, modstr, &
            'msktxi', 'units', c='kg/(kg*sec)')
          CALL channel_halt(substr, status)
        ENDIF

        !Now we know how much memory is needed.
        !Create channel_objects for tendencies of prognostic variables.
        IF (p_parallel_io) WRITE(*,*) ' ... mskt'
        ALLOCATE(mtskiptend(ntrac_gp))
        DO jt = 1, ntrac_gp
          IF (mapmem_tend(jt)) THEN
            CALL new_channel_object(status, modstr, &
              TRIM(TRIM('mskt')//TRIM(ti_gp(jt)%tp%ident%fullname)), &
              p3=mtskiptend(jt)%ptr, reprid=GP_3D_MID, lrestreq=.TRUE. )
            CALL channel_halt(substr, status)
            CALL new_attribute(status, modstr, &
              TRIM(TRIM('mskt')//TRIM(ti_gp(jt)%tp%ident%fullname)), &
              'units', c=TRIM(ti_gp(jt)%tp%ident%unit)//'/sec')
            CALL channel_halt(substr, status)
          ENDIF
        ENDDO 

        !budget of negatives
        IF (lbudneg) THEN
! op_pj_20160511+
!!$          IF (lgmxe_mtskip) THEN
          IF (lgmxe_mtskip .OR. lmecca_mtskip) THEN
! op_pj_20160511-

            IF (p_parallel_io) WRITE(*,*) ' ... budnegq'
            CALL new_channel_object(status, modstr, &
              'budnegq', p3=mtskiptend_q_budneg, reprid=GP_3D_MID, lrestreq=.TRUE. )
            CALL channel_halt(substr, status)
            CALL new_attribute(status, modstr, &
              'budnegq', 'units', c='kg/(kg*sec)')
            CALL channel_halt(substr, status)
! op_pj_20160511+
          END IF

          IF (lgmxe_mtskip) THEN
! op_pj_20160511-

            IF (p_parallel_io) WRITE(*,*) ' ... budnegxl'
            CALL new_channel_object(status, modstr, &
              'budnegxl', p3=mtskiptend_xl_budneg, reprid=GP_3D_MID, lrestreq=.TRUE. )
            CALL channel_halt(substr, status)
            CALL new_attribute(status, modstr, &
              'budnegxl', 'units', c='kg/(kg*sec)')
            CALL channel_halt(substr, status)

            IF (p_parallel_io) WRITE(*,*) ' ... budnegxi'
            CALL new_channel_object(status, modstr, &
              'budnegxi', p3=mtskiptend_xi_budneg, reprid=GP_3D_MID, lrestreq=.TRUE. )
            CALL channel_halt(substr, status)
            CALL new_attribute(status, modstr, &
              'budnegxi', 'units', c='kg/(kg*sec)')
            CALL channel_halt(substr, status)
          ENDIF

          ALLOCATE(mtskiptend_budneg(ntrac_gp))
          DO jt = 1, ntrac_gp
            IF (mapmem_tend(jt)) THEN
              CALL new_channel_object(status, modstr, &
                TRIM(TRIM('budneg')//TRIM(ti_gp(jt)%tp%ident%fullname)), &
                p3=mtskiptend_budneg(jt)%ptr, reprid=GP_3D_MID, lrestreq=.TRUE. )
              CALL channel_halt(substr, status)
              CALL new_attribute(status, modstr, &
                TRIM(TRIM('budneg')//TRIM(ti_gp(jt)%tp%ident%fullname)), &
                'units', c=TRIM(ti_gp(jt)%tp%ident%unit)//'/sec')
              CALL channel_halt(substr, status)
            ENDIF
          ENDDO 
        ENDIF
 
      ENDIF

      CALL end_message_bi(modstr, 'TENDENCY-CHANNEL-OBJECT DEFINITION', substr)

    END SUBROUTINE mtskip_init_coupling

    ! **************************************************************************

    SUBROUTINE mtskip_global_start

      !mtskip timing calculation have to be done here

      ! ECHAM5/MESSy
      USE messy_main_mpi_bi,    ONLY: p_parallel_io, p_io, p_bcast
      ! op_pj_20160616+
      USE messy_main_blather_bi, ONLY: error_bi
      ! op_pj_20160616-

      ! MESSY
      USE messy_main_timer,     ONLY: time_step_len, delta_time,         &
                                      nstep => current_time_step,        &
                                      current_date, add_date,            &
                                      time_days,                         &
                                      lstart, lresume ! op_pj_20150708

      USE messy_main_grid_def_mem_bi,   ONLY: nproma, npromz, ngpblks
! op_pj_20150707+: non-MESSy-conform usage of other submodel
!!$      USE messy_orbit_si,       ONLY: calc_orbit_param
! This has been modified: MTSKIP calculates here the offset
! dt_offset in sconds. For coupling, select this offset for orbit calculations
! in orbit.nml (c_offset = 'mtskip','dt_offset',).
! op_pj_20150707-
      USE messy_main_channel_error_bi, ONLY: channel_halt
      USE messy_main_channel,          ONLY: get_channel_object_info

      ! mz_bs_20151119+ 
      USE messy_main_constants_mem, ONLY: I8
      ! mz_bs_20151119- 

      IMPLICIT NONE

      ! mz_bs_20151119+ 
      INTRINSIC NINT, ABS
      ! mz_bs_20151119- 

      ! LOCAL
      CHARACTER(LEN=*), PARAMETER :: substr = 'mtskip_global_start'
 
      INTEGER                     :: isw_mtskip = 0  !0 never skip, always newly calculate the process
                                                     !1 skipping active, newly calculate the process
                                                     !2 skipping active, only add tendency
      INTEGER                     :: status
!!$   TYPE (time_days)            :: mtskip_date
      LOGICAL                     :: l_trig_mtskip 
!!$   INTEGER                     :: i, zproma
!!$   REAL(DP)                    :: zcdisse, zdec, zra
!!$   REAL(DP), DIMENSION(nproma,ngpblks) :: zcossza, zcosszac, zrdayl
      LOGICAL                     :: lex    ! op_pj_20150708
      ! mz_bs_20151119+ 
      LOGICAL                     :: lexb, lexc, lexdummy
      INTEGER(I8)                 :: rst_fsize ! size of mtskip restart file
      INTEGER                     :: ilmecca_mtskip, ilgmxe_mtskip !,ilmsbm_mtskip ...
      ! mz_bs_20151119- 

      ! op_pj_20150708+
      IF ( lstart ) THEN
         !new start
         !consider nstep_always
         nstep_always_mtskip = (nstep_always/trigmtskip)*trigmtskip  &
              + trigmtskip-1
      ENDIF
      !
      ! check the restart case here (see comment at original position)
      IF ( lresume ) THEN
         ! check, if data have been read from restart file
         ! (make sure to set the "ignore"-flag for all mstskip objects in
         ! channel.nml !!!)
         CALL get_channel_object_info(status, modstr, 'tsl' &
              , lrestart_read=lex)
         CALL channel_halt(substr, status)
         !mz_bs_20151119+ 
         !The following check has been reintroduced into the code
         !to avoid the error-prone dependence of correctly setting 
         !the "ignore"-flag in channel.nml. It covers the case, where
         !the mtskip restart file exists, but with no tendencies.  
         lexb = .TRUE.
         rst_fsize=4000000
         lexdummy=.FALSE.
#ifndef LF
         INQUIRE(file=TRIM('restart_'//modstr//'.nc'),EXIST=lexdummy, &
                 SIZE=rst_fsize)
#else
         ! op_pj_20160616+: SIZE is not available in Fortran95 standard!
         CALL error_bi( &
          'LF compiler does not support SIZE specifier in INQUIRE statement' &
          , substr)
         ! op_pj_20160616-
#endif
         IF ( ( (rst_fsize .LT. 2000000) .AND. (rst_fsize .GE. 0) ) &
              .OR. (.NOT. lexdummy) ) lexb = .FALSE. 

!        Additional logic is introduced here, to make the usage of 
!        mtskip less errorprone. For the case that after restart due
!        to changes in the namelist, either the number of submodules
!        undergoing skipping or the skipping timestep-length (trigmtskip or
!        delta_time) has changed, a resynchronization due to nstep_always_mtskip
!        is enforced. Note that this also can be done, by simply deleting the
!        restart-file or by setting the "ignore" flag on true for mtskip in
!        channel.nml .
     
         lexc = .TRUE.
         ilmecca_mtskip = 0
         IF ( lmecca_mtskip ) ilmecca_mtskip = 1
         IF ( ilmecca_mtskip .NE. NINT(rmecca) ) lexc = .FALSE.
         ilgmxe_mtskip = 0
         IF ( lgmxe_mtskip ) ilgmxe_mtskip = 1
         IF ( ilgmxe_mtskip .NE. NINT(rgmxe) )   lexc = .FALSE.  
         !ilmsbm_mtskip = 0
         !IF ( lmsbm_mtskip ) imsbm_mtskip = 1
         !IF ( lmsbm_mtskip .NE. NINT(rmsbm) )   lexc = .FALSE.  
         ! ...
         IF ( ABS(tsl_init - tsl) .GT. 1.E-5_dp ) lexc = .FALSE.

         IF ( (.NOT. lex)  .OR.  & 
              (.NOT. lexb) .OR.  & 
              (.NOT. lexc) ) THEN 
            nstep_always_mtskip = ((nstep+nstep_always)/trigmtskip)*trigmtskip &
                 + trigmtskip-1
         ELSE
            !don't consider nstep_always
            nstep_always_mtskip = 0
         END IF
         !mz_bs_20151119-
      END IF
      !
      ! SET channel object contents
      IF (lstart .OR. lresume) THEN
         ! op_pj_20150708: moved from initialize to here, to avoid
         ! that contents is overwritten from restart files
         tsl = tsl_init
         !
         IF (lmecca_mtskip) THEN
            rmecca=1._dp
         ELSE
            rmecca=0._dp
         ENDIF
         !
         IF (lgmxe_mtskip) THEN
            rgmxe=1._dp
         ELSE
            rgmxe=0._dp
         ENDIF
      END IF
      ! op_pj_20150708-

      l_trig_mtskip = .FALSE.
      isw_mtskip = 0
      rmode = 0._dp
      IF ( (nstep .GT. nstep_always_mtskip) .AND. lmodules_mtskip ) THEN
        l_trig_mtskip = mtskip_event_trig(trigmtskip, nstep)
        tslphys = tsl
        IF ( l_trig_mtskip ) THEN 
          isw_mtskip = 1
          rmode = 1._dp
!!$       !calculate cdisse and cossza for jval
!!$       mtskip_date = current_date
          IF (lmecca_mtskip) THEN
!!$         CALL add_date(0,INT(0.5_dp*delta_time)*trigmtskip,mtskip_date)
            dt_offset = 0.5_dp*delta_time*REAL(trigmtskip,dp)
          ELSE
!!$         CALL add_date(0,INT(phtl_offset*delta_time),mtskip_date)
            dt_offset = phtl_offset*delta_time
          ENDIF
!!$       CALL calc_orbit_param(status, mtskip_date, cdissephtl, decphtl, raphtl,       &
!!$                             cosszaphtl, cosszacphtl, rdaylphtl)
!!$
!!$       IF (status /=0) &
!!$         CALL error_bi('calc_orbit_param reported an error',substr)
        ELSE
          isw_mtskip = 2
          rmode = 2._dp
          IF ( .NOT. lmecca_mtskip) THEN
!!$         !calculate cdisse and cossza for jval
!!$         mtskip_date = current_date
!!$         CALL add_date(0,INT(phtl_offset*delta_time),mtskip_date)
            dt_offset = phtl_offset*delta_time
!!$         CALL calc_orbit_param(status, mtskip_date, cdissephtl, decphtl, raphtl,       &
!!$                               cosszaphtl, cosszacphtl, rdaylphtl)
!!$         IF (status /=0) &
!!$           CALL error_bi('calc_orbit_param reported an error',substr)
          ENDIF
        ENDIF
      ELSE
        isw_mtskip = 0
        rmode = 0._dp
        tslphys = time_step_len
!!$     !calculate cdisse and cossza for jval
!!$     mtskip_date = current_date
!!$     CALL add_date(0,INT(phtl_offset*delta_time),mtskip_date)
        dt_offset = phtl_offset*delta_time
!!$     CALL calc_orbit_param(status, mtskip_date, cdissephtl, decphtl, raphtl,       &
!!$                           cosszaphtl, cosszacphtl, rdaylphtl)
!!$     IF (status /=0) &
!!$       CALL error_bi('calc_orbit_param reported an error',substr)
     ENDIF

    END SUBROUTINE mtskip_global_start

    ! **************************************************************************

    SUBROUTINE mtskip_local_start

      CALL mtskip_add_tend

    END SUBROUTINE mtskip_local_start


    ! **************************************************************************

    SUBROUTINE mtskip_free_memory

    IMPLICIT NONE 
    INTRINSIC ASSOCIATED

    IF (ASSOCIATED(mapmem_tend))       DEALLOCATE(mapmem_tend) 
    IF (ASSOCIATED(mtskiptend))        DEALLOCATE(mtskiptend) 
    IF (ASSOCIATED(mtskiptend_budneg)) DEALLOCATE(mtskiptend_budneg) 

    END SUBROUTINE mtskip_free_memory

    ! **************************************************************************

    ! **************************************************************************
    ! PRIVATE ROUTINES
    ! **************************************************************************

    ! **************************************************************************

    SUBROUTINE mtskip_add_tend

      !Add collected tendencies of submodels that are part of skipping on
      !general tendency field.  

      ! ECHAM5/MESSy
      USE messy_main_mpi_bi,         ONLY: p_parallel_io, p_io, p_bcast

      ! MESSY
      USE messy_main_timer,          ONLY: time_step_len, delta_time,      &
                                          nstep => current_time_step
      USE messy_main_grid_def_mem_bi,ONLY: nlev, jrow, kproma

!#ifndef MESSYTENDENCY
      USE messy_main_tracer_mem_bi,  ONLY: pxtte => qxtte, pxtm1 => qxtm1, &
                                             ntrac_gp, ti_gp
!#else
!      USE messy_main_tracer_mem_bi,  ONLY: ntrac_gp, ti_gp
!#endif

       USE messy_main_data_bi,       ONLY:  qm1       & ! specific humidity                      [kg kg-1]
                                          , qte_3d    & ! specific humidity tendency             [kg kg-1 s-1]
                                          , xlm1      & ! cloud water                            [kg kg-1]
                                          , xlte_3d   & ! cloud water tendency                   [kg kg-1 s-1]
                                          , xim1      & ! cloud ice                              [kg kg-1]
                                          , xite_3d   & ! cloud ice tendency                     [kg kg-1 s-1]
                                          , tm1          ! dry air temperature (at time step -1)  [K]

      IMPLICIT NONE

      ! LOCAL
      CHARACTER(LEN=*), PARAMETER :: substr = 'mtskip_add_tend'

      REAL(DP), DIMENSION(:,:), ALLOCATABLE :: ztend
      REAL(DP)                              :: ztend_help
                                               
      REAL(DP)                              :: rtest

      INTEGER                     :: isw = 0  !0 never skip, always newly calculate the process
                                              !1 skipping active, newly calculate the process
                                              !2 skipping active, only add tendency
      INTEGER                     :: jt, jk, jl


      isw = NINT(rmode)

      SELECT CASE (isw)

        CASE(0) 
          RETURN

        CASE(1)
          !set mtskip tendency-fields to zero
          DO jt = 1, ntrac_gp
            IF (mapmem_tend(jt)) THEN
              DO jk=1,nlev
                DO jl=1,kproma
                  mtskiptend(jt)%ptr(_RI_XYZ__(jl,jrow,jk)) = 0._dp
                ENDDO
              ENDDO
              IF (lbudneg) THEN
                DO jk=1,nlev
                  DO jl=1,kproma
                    mtskiptend_budneg(jt)%ptr(_RI_XYZ__(jl,jrow,jk)) = 0._dp
                  ENDDO
                ENDDO
              ENDIF
            ENDIF
          ENDDO
          ! op_pj_20160511+
          IF (lgmxe_mtskip .OR. lmecca_mtskip) THEN
             mtskiptend_q(_RI_XYZ__(1:kproma,jrow,1:nlev)) = 0._dp
             IF (lbudneg) THEN
                mtskiptend_q_budneg(_RI_XYZ__(1:kproma,jrow,1:nlev)) = 0._dp
             END IF
          END IF
          ! op_pj_20160511-
          IF (lgmxe_mtskip) THEN
            mtskiptend_xl(_RI_XYZ__(1:kproma,jrow,1:nlev)) = 0._dp
            mtskiptend_xi(_RI_XYZ__(1:kproma,jrow,1:nlev)) = 0._dp
            IF (lbudneg) THEN
              mtskiptend_xl_budneg(_RI_XYZ__(1:kproma,jrow,1:nlev)) = 0._dp
              mtskiptend_xi_budneg(_RI_XYZ__(1:kproma,jrow,1:nlev)) = 0._dp
            ENDIF
          ENDIF
          RETURN
 
        CASE(2)

          !add tendencies, when the process is not calculated

          IF (lfixneg) THEN
            ALLOCATE(ztend(kproma,nlev))
          ENDIF

! op_pj_20160511+ qm1/qte extracted (for MECCA and/or GMXE)
          IF ( lgmxe_mtskip .OR. lmecca_mtskip ) THEN
            !water vapour, cloud liquid water, cloud ice water
   
            IF (lfixneg) THEN
              IF (lbudneg) THEN
                !lfixneg + lbudneg
                DO jk = 1,nlev
                  DO jl = 1,kproma
                    rtest = qm1(_RI_XYZ__(jl,jrow,jk)) + time_step_len * &
                            ( qte_3d(_RI_XYZ__(jl,jrow,jk)) + &
                            mtskiptend_q(_RI_XYZ__(jl,jrow,jk)) )
                    IF ( rtest .LT. 0.0_dp ) THEN
                      ztend_help = &
                           -1._dp*(qm1(_RI_XYZ__(jl,jrow,jk))/time_step_len+&
                                   qte_3d(_RI_XYZ__(jl,jrow,jk)))
                      mtskiptend_q_budneg(_RI_XYZ__(jl,jrow,jk)) = &
                           ztend_help - mtskiptend_q(_RI_XYZ__(jl,jrow,jk))
                      qte_3d(_RI_XYZ__(jl,jrow,jk)) =   &
                           qte_3d(_RI_XYZ__(jl,jrow,jk)) + ztend_help
                    ELSE
                      qte_3d(_RI_XYZ__(jl,jrow,jk)) = &
                           qte_3d(_RI_XYZ__(jl,jrow,jk)) &
                           + mtskiptend_q(_RI_XYZ__(jl,jrow,jk))
                      mtskiptend_q_budneg(_RI_XYZ__(jl,jrow,jk)) = 0.0_dp
                    ENDIF
                  ENDDO
                ENDDO

              ELSE
                !lfixneg only
                DO jk = 1,nlev
                  DO jl = 1,kproma
                    rtest = qm1(_RI_XYZ__(jl,jrow,jk)) + time_step_len * &
                            ( qte_3d(_RI_XYZ__(jl,jrow,jk)) &
                            + mtskiptend_q(_RI_XYZ__(jl,jrow,jk)))
                    IF ( rtest .LT. 0.0_dp ) THEN
                      qte_3d(_RI_XYZ__(jl,jrow,jk)) =  &
                           -1._dp*(qm1(_RI_XYZ__(jl,jrow,jk))/time_step_len)
                    ELSE
                      qte_3d(_RI_XYZ__(jl,jrow,jk)) =  &
                           qte_3d(_RI_XYZ__(jl,jrow,jk)) &
                                  + mtskiptend_q(_RI_XYZ__(jl,jrow,jk))
                    ENDIF
                  ENDDO
                ENDDO

              ENDIF !lbudneg
            ELSE
              !no fixing of negatives
              IF (lbudneg) THEN
                DO jk = 1,nlev
                  DO jl = 1,kproma
                    rtest = qm1(_RI_XYZ__(jl,jrow,jk)) + time_step_len * &
                            ( qte_3d(_RI_XYZ__(jl,jrow,jk)) + &
                            mtskiptend_q(_RI_XYZ__(jl,jrow,jk)) )
                    IF ( rtest .LT. 0.0_dp ) THEN
                      ztend_help = &
                           -1._dp*(qm1(_RI_XYZ__(jl,jrow,jk))/time_step_len+&
                           qte_3d(_RI_XYZ__(jl,jrow,jk)))
                      mtskiptend_q_budneg(_RI_XYZ__(jl,jrow,jk)) = &
                           ztend_help - mtskiptend_q(_RI_XYZ__(jl,jrow,jk))
                    ELSE
                      mtskiptend_q_budneg(_RI_XYZ__(jl,jrow,jk)) = 0.0_dp
                    ENDIF
                  ENDDO
                ENDDO
              END IF
              qte_3d(_RI_XYZ__(1:kproma,jrow,1:nlev)) =&
                   qte_3d(_RI_XYZ__(1:kproma,jrow,1:nlev)) &
                               + mtskiptend_q(_RI_XYZ__(1:kproma,jrow,1:nlev))
   
            ENDIF !lfixneg
          ENDIF !lgmxe
! op_pj_20160511+

          IF ( lgmxe_mtskip ) THEN
            !water vapour, cloud liquid water, cloud ice water
   
            IF (lfixneg) THEN
              IF (lbudneg) THEN
                !lfixneg + lbudneg
! op_pj_20160511+ extracted and moved above
!!$                DO jk = 1,nlev
!!$                  DO jl = 1,kproma
!!$                    rtest = qm1(_RI_XYZ__(jl,jrow,jk)) + time_step_len * &
!!$                            ( qte_3d(_RI_XYZ__(jl,jrow,jk)) + mtskiptend_q(_RI_XYZ__(jl,jrow,jk)) )
!!$                    IF ( rtest .LT. 0.0_dp ) THEN
!!$                      ztend_help = -1._dp*(qm1(_RI_XYZ__(jl,jrow,jk))/time_step_len+qte_3d(_RI_XYZ__(jl,jrow,jk)))
!!$                      mtskiptend_q_budneg(_RI_XYZ__(jl,jrow,jk)) = ztend_help - mtskiptend_q(_RI_XYZ__(jl,jrow,jk))
!!$                      qte_3d(_RI_XYZ__(jl,jrow,jk)) =  qte_3d(_RI_XYZ__(jl,jrow,jk)) &
!!$                                  + ztend_help
!!$                    ELSE
!!$                      qte_3d(_RI_XYZ__(jl,jrow,jk)) =  qte_3d(_RI_XYZ__(jl,jrow,jk)) &
!!$                                  + mtskiptend_q(_RI_XYZ__(jl,jrow,jk))
!!$                      mtskiptend_q_budneg(_RI_XYZ__(jl,jrow,jk)) = 0.0_dp
!!$                    ENDIF
!!$                  ENDDO
!!$                ENDDO
! op_pj_20160511-
                
                DO jk = 1,nlev
                  DO jl = 1,kproma
                    rtest = xlm1(_RI_XYZ__(jl,jrow,jk)) + time_step_len * &
                            ( xlte_3d(_RI_XYZ__(jl,jrow,jk)) +&
                            mtskiptend_xl(_RI_XYZ__(jl,jrow,jk)) )
                    IF ( rtest .LT. 0.0_dp ) THEN
                      ztend_help = &
                           -1._dp*(xlm1(_RI_XYZ__(jl,jrow,jk))/time_step_len+&
                           xlte_3d(_RI_XYZ__(jl,jrow,jk)))
                      mtskiptend_xl_budneg(_RI_XYZ__(jl,jrow,jk)) = &
                           ztend_help - mtskiptend_xl(_RI_XYZ__(jl,jrow,jk))
                      xlte_3d(_RI_XYZ__(jl,jrow,jk)) =    &
                           xlte_3d(_RI_XYZ__(jl,jrow,jk)) &
                                   + ztend_help
                    ELSE
                      xlte_3d(_RI_XYZ__(jl,jrow,jk)) =  &
                           xlte_3d(_RI_XYZ__(jl,jrow,jk)) &
                                  + mtskiptend_xl(_RI_XYZ__(jl,jrow,jk))
                      mtskiptend_xl_budneg(_RI_XYZ__(jl,jrow,jk)) = 0.0_dp
                    ENDIF
                  ENDDO
                ENDDO
                
                DO jk = 1,nlev
                  DO jl = 1,kproma
                    rtest = xim1(_RI_XYZ__(jl,jrow,jk)) + time_step_len * &
                            ( xite_3d(_RI_XYZ__(jl,jrow,jk)) + &
                            mtskiptend_xi(_RI_XYZ__(jl,jrow,jk)) )
                    IF ( rtest .LT. 0.0_dp ) THEN
                      ztend_help = &
                           -1._dp*(xim1(_RI_XYZ__(jl,jrow,jk))/time_step_len+&
                           xite_3d(_RI_XYZ__(jl,jrow,jk)))
                      mtskiptend_xi_budneg(_RI_XYZ__(jl,jrow,jk)) = &
                           ztend_help - mtskiptend_xi(_RI_XYZ__(jl,jrow,jk))
                      xite_3d(_RI_XYZ__(jl,jrow,jk)) =  &
                           xite_3d(_RI_XYZ__(jl,jrow,jk)) &
                                   + ztend_help
                    ELSE
                      xite_3d(_RI_XYZ__(jl,jrow,jk)) =  &
                           xite_3d(_RI_XYZ__(jl,jrow,jk)) &
                                  + mtskiptend_xi(_RI_XYZ__(jl,jrow,jk))
                      mtskiptend_xi_budneg(_RI_XYZ__(jl,jrow,jk)) = 0.0_dp
                    ENDIF
                  ENDDO
                ENDDO


              ELSE
                !lfixneg only
! op_pj_20160511+ extracted and moved above
!!$                DO jk = 1,nlev
!!$                  DO jl = 1,kproma
!!$                    rtest = qm1(_RI_XYZ__(jl,jrow,jk)) + time_step_len * &
!!$                            ( qte_3d(_RI_XYZ__(jl,jrow,jk)) + mtskiptend_q(_RI_XYZ__(jl,jrow,jk)) )
!!$                    IF ( rtest .LT. 0.0_dp ) THEN
!!$                      qte_3d(_RI_XYZ__(jl,jrow,jk)) = -1._dp*(qm1(_RI_XYZ__(jl,jrow,jk))/time_step_len)
!!$                    ELSE
!!$                      qte_3d(_RI_XYZ__(jl,jrow,jk)) =  qte_3d(_RI_XYZ__(jl,jrow,jk)) &
!!$                                  + mtskiptend_q(_RI_XYZ__(jl,jrow,jk))
!!$                    ENDIF
!!$                  ENDDO
!!$                ENDDO
! op_pj_20160511-                

                DO jk = 1,nlev
                  DO jl = 1,kproma
                    rtest = xlm1(_RI_XYZ__(jl,jrow,jk)) + time_step_len * &
                            ( xlte_3d(_RI_XYZ__(jl,jrow,jk)) + &
                            mtskiptend_xl(_RI_XYZ__(jl,jrow,jk)) )
                    IF ( rtest .LT. 0.0_dp ) THEN
                      xlte_3d(_RI_XYZ__(jl,jrow,jk)) =&
                           -1._dp*(xlm1(_RI_XYZ__(jl,jrow,jk))/time_step_len)
                    ELSE
                      xlte_3d(_RI_XYZ__(jl,jrow,jk)) = &
                           xlte_3d(_RI_XYZ__(jl,jrow,jk)) &
                                  + mtskiptend_xl(_RI_XYZ__(jl,jrow,jk))
                    ENDIF
                  ENDDO
                ENDDO
                
                DO jk = 1,nlev
                  DO jl = 1,kproma
                    rtest = xim1(_RI_XYZ__(jl,jrow,jk)) + time_step_len * &
                            ( xite_3d(_RI_XYZ__(jl,jrow,jk)) + &
                            mtskiptend_xi(_RI_XYZ__(jl,jrow,jk)) )
                    IF ( rtest .LT. 0.0_dp ) THEN
                      xite_3d(_RI_XYZ__(jl,jrow,jk)) = &
                           -1._dp*(xim1(_RI_XYZ__(jl,jrow,jk))/time_step_len)
                    ELSE
                      xite_3d(_RI_XYZ__(jl,jrow,jk)) =  &
                           xite_3d(_RI_XYZ__(jl,jrow,jk)) &
                                  + mtskiptend_xi(_RI_XYZ__(jl,jrow,jk))
                    ENDIF
                  ENDDO
                ENDDO
              ENDIF !lbudneg
            ELSE
              !no fixing of negatives
              IF (lbudneg) THEN
! op_pj_20160511+ extracted and moved above
!!$                DO jk = 1,nlev
!!$                  DO jl = 1,kproma
!!$                    rtest = qm1(_RI_XYZ__(jl,jrow,jk)) + time_step_len * &
!!$                            ( qte_3d(_RI_XYZ__(jl,jrow,jk)) + mtskiptend_q(_RI_XYZ__(jl,jrow,jk)) )
!!$                    IF ( rtest .LT. 0.0_dp ) THEN
!!$                      ztend_help = -1._dp*(qm1(_RI_XYZ__(jl,jrow,jk))/time_step_len+qte_3d(_RI_XYZ__(jl,jrow,jk)))
!!$                      mtskiptend_q_budneg(_RI_XYZ__(jl,jrow,jk)) = ztend_help - mtskiptend_q(_RI_XYZ__(jl,jrow,jk))
!!$                    ELSE
!!$                      mtskiptend_q_budneg(_RI_XYZ__(jl,jrow,jk)) = 0.0_dp
!!$                    ENDIF
!!$                  ENDDO
!!$                ENDDO
! op_pj_20160511-
                
                DO jk = 1,nlev
                  DO jl = 1,kproma
                    rtest = xlm1(_RI_XYZ__(jl,jrow,jk)) + time_step_len * &
                            ( xlte_3d(_RI_XYZ__(jl,jrow,jk)) &
                            + mtskiptend_xl(_RI_XYZ__(jl,jrow,jk)) )
                    IF ( rtest .LT. 0.0_dp ) THEN
                      ztend_help = &
                           -1._dp*(xlm1(_RI_XYZ__(jl,jrow,jk))/time_step_len+&
                           xlte_3d(_RI_XYZ__(jl,jrow,jk)))
                      mtskiptend_xl_budneg(_RI_XYZ__(jl,jrow,jk)) = &
                           ztend_help - mtskiptend_xl(_RI_XYZ__(jl,jrow,jk))
                    ELSE
                      mtskiptend_xl_budneg(_RI_XYZ__(jl,jrow,jk)) = 0.0_dp
                    ENDIF
                  ENDDO
                ENDDO
                
                DO jk = 1,nlev
                  DO jl = 1,kproma
                    rtest = xim1(_RI_XYZ__(jl,jrow,jk)) + time_step_len * &
                            ( xite_3d(_RI_XYZ__(jl,jrow,jk)) + &
                            mtskiptend_xi(_RI_XYZ__(jl,jrow,jk)) )
                    IF ( rtest .LT. 0.0_dp ) THEN
                      ztend_help =&
                           -1._dp*(xim1(_RI_XYZ__(jl,jrow,jk))/time_step_len+&
                           xite_3d(_RI_XYZ__(jl,jrow,jk)))
                      mtskiptend_xi_budneg(_RI_XYZ__(jl,jrow,jk)) = &
                           ztend_help - mtskiptend_xi(_RI_XYZ__(jl,jrow,jk))
                    ELSE
                      mtskiptend_xi_budneg(_RI_XYZ__(jl,jrow,jk)) = 0.0_dp
                    ENDIF
                  ENDDO
                ENDDO
              ENDIF

              xlte_3d(_RI_XYZ__(1:kproma,jrow,1:nlev)) =&
                   xlte_3d(_RI_XYZ__(1:kproma,jrow,1:nlev)) &
                   + mtskiptend_xl(_RI_XYZ__(1:kproma,jrow,1:nlev))
  
              xite_3d(_RI_XYZ__(1:kproma,jrow,1:nlev)) =  &
                   xite_3d(_RI_XYZ__(1:kproma,jrow,1:nlev))&
                              + mtskiptend_xi(_RI_XYZ__(1:kproma,jrow,1:nlev))
            ENDIF !lfixneg
          ENDIF !lgmxe

          !tracer

          DO jt=1,ntrac_gp
            IF (mapmem_tend(jt)) THEN
              ztend(1:kproma,1:nlev) = &
                   mtskiptend(jt)%ptr(_RI_XYZ__(1:kproma,jrow,1:nlev))
              IF (lfixneg) THEN
                IF (lbudneg) THEN
                  !lfixneg + lbudneg
                  DO jk = 1,nlev
                    DO jl = 1,kproma
                      rtest = pxtm1(_RI_X_ZN_(jl,jk,jt)) + time_step_len * &
                              ( pxtte(_RI_X_ZN_(jl,jk,jt)) + ztend(jl,jk) )
                      IF ( rtest .LT. 0.0_dp ) THEN
                        ztend_help = &
                             -1._dp*( pxtm1(_RI_X_ZN_(jl,jk,jt))/time_step_len &
                                            + pxtte(_RI_X_ZN_(jl,jk,jt)) )
                        mtskiptend_budneg(jt)%ptr(_RI_XYZ__(jl,jrow,jk)) = &
                                             ztend_help - ztend(jl,jk)
                        ztend(jl,jk) = ztend_help 
                      ELSE
                        mtskiptend_budneg(jt)%ptr(_RI_XYZ__(jl,jrow,jk)) = 0.0_dp
                      ENDIF
                    ENDDO
                  ENDDO
                ELSE
                  !lfixneg
                  DO jk = 1,nlev
                    DO jl = 1,kproma
                      rtest = pxtm1(_RI_X_ZN_(jl,jk,jt)) + time_step_len * &
                              ( pxtte(_RI_X_ZN_(jl,jk,jt)) + ztend(jl,jk) )
                      IF ( rtest .LT. 0.0_dp ) &
                        ztend(jl,jk) = &
                        -1._dp*( pxtm1(_RI_X_ZN_(jl,jk,jt))/time_step_len & 
                                              + pxtte(_RI_X_ZN_(jl,jk,jt)) )
                    ENDDO
                  ENDDO
                ENDIF
              ELSE
                !no fixneg
                IF (lbudneg) THEN
                  !lbudneg
                  DO jk = 1,nlev
                    DO jl = 1,kproma
                      rtest = pxtm1(_RI_X_ZN_(jl,jk,jt))  + time_step_len * &
                              ( pxtte(_RI_X_ZN_(jl,jk,jt))  + ztend(jl,jk) )
                      IF ( rtest .LT. 0.0_dp ) THEN
                        ztend_help = &
                            -1._dp*( pxtm1(_RI_X_ZN_(jl,jk,jt)) /time_step_len &
                                            + pxtte(_RI_X_ZN_(jl,jk,jt))  )
                        mtskiptend_budneg(jt)%ptr(_RI_XYZ__(jl,jrow,jk)) = &
                                             ztend_help - ztend(jl,jk)
                      ELSE
                        mtskiptend_budneg(jt)%ptr(_RI_XYZ__(jl,jrow,jk)) = 0.0_dp
                      ENDIF
                    ENDDO
                  ENDDO
                ENDIF
              ENDIF      
              !add tendency
#ifdef MESSYTENDENCY
!!$              CALL mtend_add_l(my_handle, mtend_id_tracer, &
!!$                               px=ztend(1:kproma,1:nlev), idt=jt)
              CALL mtend_add_l(my_handle, jt, px=ztend(1:kproma,1:nlev))
#else
              pxtte(_RI_X_ZN_(1:kproma,1:nlev,jt)) = &
                   pxtte(_RI_X_ZN_(1:kproma,1:nlev,jt)) &
                                          + ztend(1:kproma,1:nlev)
#endif
            ENDIF
          ENDDO

          IF (lfixneg) THEN
            DEALLOCATE(ztend)
          ENDIF

      END SELECT

      RETURN

    END SUBROUTINE mtskip_add_tend

    ! **************************************************************************

    SUBROUTINE mtskip_read_nml_cpl(status, iou)

      ! MESSy
      USE messy_main_tools, ONLY: read_nml_open, read_nml_check, read_nml_close

      IMPLICIT NONE

      ! I/O
      INTEGER, INTENT(OUT) :: status     ! error status
      INTEGER, INTENT(IN)  :: iou        ! I/O unit

      ! LOCAL
      CHARACTER(LEN=*), PARAMETER :: substr='mtskip_read_nml_cpl'
      LOGICAL              :: lex      ! file exists ?
      INTEGER              :: fstat    ! file status

      status = 1

      CALL read_nml_open(lex, substr, iou, 'CPL', modstr)
      IF (.NOT.lex) RETURN    ! <modstr>.nml does not exist

      READ(iou, NML=CPL, IOSTAT=fstat)
      CALL read_nml_check(fstat, substr, iou, 'CPL', modstr)
      IF (fstat /= 0) RETURN  ! error while reading namelist

      CALL read_nml_close(substr, iou, modstr)
      status = 0 ! NO ERROR

    END SUBROUTINE mtskip_read_nml_cpl

END MODULE messy_mtskip_si

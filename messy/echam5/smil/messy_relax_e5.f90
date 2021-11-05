! **********************************************************************
!
! SUBMODEL INTERFACE LAYER (SMIL) ROUTINES FOR MESSy SUBMODEL RELAX
!
! Author : Hella Garny, DLR-IPA, June 2016
!
! References: see messy_relax.f90
!
! **********************************************************************

! **********************************************************************
MODULE messy_relax_e5
! **********************************************************************

  ! MESSy BASEMODEL INTERFACE LAYER (BMIL)
  USE messy_main_blather_bi,    ONLY: start_message_bi, end_message_bi, &
                                      error_bi
!                                     warning_bi ! not used in this module
!                                     info_bi    ! not used in this module
#ifdef MESSYTENDENCY
  USE messy_main_tendency_bi,   ONLY:  mtend_get_handle,             &
                                       mtend_add_l,                  &
                                       mtend_register,               &
                                       mtend_id_t,                   &
                                       mtend_id_u, mtend_id_v
#endif

  ! SMCL
  USE messy_main_constants_mem, ONLY: STRLEN_LONG ! op_pj_20180718
!!$  USE messy_main_channel,       ONLY:  t_chaobj_cpl ! op_pj_20180718
  USE messy_main_channel,       ONLY: STRLEN_CHANNEL
  USE messy_relax

  IMPLICIT NONE
  INTRINSIC :: NULL
  PRIVATE
  SAVE

  ! op_pj_20100827+
  TYPE t_chaobj_cpl_rel
     CHARACTER(LEN=STRLEN_CHANNEL)  :: cha  = ''
     CHARACTER(LEN=2*STRLEN_LONG+5) :: obj  = ''
  END TYPE t_chaobj_cpl_rel
  ! op_pj_20100827-

  ! CPL-NAMELIST PARAMETERS
  ! op_pj_20100827: type changed from t_chaobj_cpl to t_chaobj_cpl_rel
  LOGICAL                             :: lrayfr            = .TRUE.
  LOGICAL                             :: lnewco            = .TRUE.
  LOGICAL                             :: liheat_cc_tropics = .FALSE. !op_rw_20190212
  LOGICAL                             :: liheat_waves      = .FALSE. !op_rw_20181102
  LOGICAL                             :: liheat_mons       = .FALSE. !op_mn_20180620
  LOGICAL                             :: l_no_polar_vortex = .FALSE. !op_rw_20190710
  LOGICAL                             :: l_Butler_heat     = .TRUE.  !op_rw_20200213
  TYPE(t_chaobj_cpl_rel)              :: rayfr_k_inp
  TYPE(t_chaobj_cpl_rel)              :: newco_t_inp
  TYPE(t_chaobj_cpl_rel)              :: newco_k_inp
  TYPE(t_chaobj_cpl_rel)              :: cct_h_inp   !op_rw_20190212
  TYPE(t_chaobj_cpl_rel)              :: waves_h_inp !op_rw_20181102
  TYPE(t_chaobj_cpl_rel)              :: reght_h_inp !op_mn_20180620
  TYPE(t_chaobj_cpl_rel)              :: tmpht_h_inp !op_mn_20180620


  ! WORKSPACE
  ! pointer to channel objects
  REAL(DP), DIMENSION (:,:,:), POINTER :: kdamp           => NULL()
  REAL(DP), DIMENSION (:,:,:), POINTER :: tequ            => NULL()
  REAL(DP), DIMENSION (:,:,:), POINTER :: kappa           => NULL()
  REAL(DP), DIMENSION (:,:,:), POINTER :: tteh_cc_tropics => NULL() !op_rw_20190212
  REAL(DP), DIMENSION (:,:,:), POINTER :: tteh_waves      => NULL() !op_rw_20181102
  REAL(DP), DIMENSION (:,:,:), POINTER :: tteh_mons       => NULL() !op_mn_20180620


#ifdef MESSYTENDENCY
  ! variable for tendency budget
  integer                         :: my_handle
#endif

  ! PUBLIC SUBROUTINES (called from messy_main_control_e5.f90)
  ! NOTE: in case you activate further entry points, make sure to call them
  !       in messy_main_control_e5.f90
  PUBLIC :: relax_initialize    ! initialize submodel
  PUBLIC :: relax_init_memory   ! request memory
!!$  PUBLIC :: bufly_new_tracer    ! define new tracers
  PUBLIC :: relax_init_coupling ! set pointers for coupling to BM and other SMs
!!$  PUBLIC :: bufly_init_tracer   ! initilize tracers
!!$  PUBLIC :: bufly_global_start  ! entry point in time loop (all vectors)
!!$  PUBLIC :: bufly_physc         ! entry point in time loop (current vector)
  PUBLIC :: relax_physc
!!$  PUBLIC :: bufly_free_memory   ! free allocated memory

  ! PRIVATE SUBROTINES
  !PRIVATE :: relax_read_nml_cpl

CONTAINS

  ! ####################################################################
  ! PUBLIC SUBROUTINES
  ! ####################################################################

  ! ====================================================================
  SUBROUTINE relax_initialize

    ! ------------------------------------------------------------------
    ! This subroutine is used to
    ! - read (and broadcast) the CTRL-namelist,
    ! - read (and broadcast) the CPL-namelist,
    ! - perform the basic setup of the submodel.
    ! ------------------------------------------------------------------

    ! MESSy BASEMODEL INTERFACE LAYER (BMIL)
    USE messy_main_mpi_bi,    ONLY: p_parallel_io, p_io, p_bcast
    USE messy_main_tools,     ONLY: find_next_free_unit

    IMPLICIT NONE

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'relax_initialize'
    INTEGER                     :: status ! error status
    INTEGER                     :: iou    ! I/O unit

    CALL start_message_bi(modstr,'INITIALISATION',substr)  ! log-output

    ! READ CPL namelist
    IF (p_parallel_io) THEN                  ! read only on I/O-PE
       iou = find_next_free_unit(100,200)    ! find next free I/O unit
       CALL relax_read_nml_cpl(status, iou)  ! read CPL-namelist
       ! terminate if error
       IF (status /= 0) CALL error_bi('Error in reading CPL namelist',substr)
    END IF
    ! BROADCAST CPL namelist entries from I/O-PE to ALL OTHER PEs
    CALL p_bcast(lrayfr, p_io)
    CALL p_bcast(lnewco, p_io)
    CALL p_bcast(liheat_cc_tropics, p_io) !op_rw_20190212
    CALL p_bcast(liheat_waves, p_io)      !op_rw_20181102
    CALL p_bcast(liheat_mons, p_io)       !op_mn_20180620
    CALL p_bcast(l_no_polar_vortex, p_io) !op_rw_20190710
    CALL p_bcast(l_Butler_heat, p_io)     !op_rw_20200213
    CALL p_bcast(rayfr_k_inp%cha, p_io)
    CALL p_bcast(rayfr_k_inp%obj, p_io)
    CALL p_bcast(newco_t_inp%cha, p_io)
    CALL p_bcast(newco_t_inp%obj, p_io)
    CALL p_bcast(newco_k_inp%cha, p_io)
    CALL p_bcast(newco_k_inp%obj, p_io)

    ! +op_rw_20190212
    CALL p_bcast(cct_h_inp%cha, p_io)
    CALL p_bcast(cct_h_inp%obj, p_io)
    ! -op_rw_20190212

    ! +op_rw_20181102
    CALL p_bcast(waves_h_inp%cha, p_io)
    CALL p_bcast(waves_h_inp%obj, p_io)
    ! -op_rw_20181102

    ! +op_mn_20180620
    CALL p_bcast(reght_h_inp%cha, p_io)
    CALL p_bcast(reght_h_inp%obj, p_io)
    CALL p_bcast(tmpht_h_inp%cha, p_io)
    CALL p_bcast(tmpht_h_inp%obj, p_io)
    ! -op_mn_20180620




    ! ### PERFORM INITIAL SETUP (CALL RESPECTIVE SMCL ROUTINE(S)) HERE

    ! ub_ak_20181026+
    ! moved to init_memory
!!$#ifdef MESSYTENDENCY
!!$    my_handle = mtend_get_handle(modstr)
!!$    CALL mtend_register (my_handle, mtend_id_t)
!!$#endif
    ! ub_ak_20181026-


    CALL end_message_bi(modstr,'INITIALISATION',substr)  ! log-output

  END SUBROUTINE relax_initialize
  ! ====================================================================

  ! ====================================================================
  SUBROUTINE relax_init_memory

    ! ------------------------------------------------------------------
    ! This subroutine is used to request memory for the submodel.
    ! The preferable method is to use "channel objects".
    ! Allocate your own memory, only if absolutely required.
    ! ------------------------------------------------------------------

    ! BMIL
    USE messy_main_channel_error_bi, ONLY: channel_halt
    USE messy_main_channel_bi,       ONLY: GP_3D_MID
    USE messy_main_channel,          ONLY: new_channel, new_attribute
!                                             new_channel_object  ! not used in this subroutine

    IMPLICIT NONE

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'relax_init_memory'
    INTEGER                     :: status

    ! ub_ak_20181026+
    ! moved from initialize
#ifdef MESSYTENDENCY
    my_handle = mtend_get_handle(modstr)
    CALL mtend_register (my_handle, mtend_id_t)
    CALL mtend_register (my_handle, mtend_id_u)
    CALL mtend_register (my_handle, mtend_id_v)
#endif
    ! ub_ak_20181026-
    !!CALL start_message_bi(modstr,'MEMORY ALLOCATION',substr)  ! log-output
    !! ### ALLOCATE OWN MEMORY HERE, BUT ONLY IF ABSOLUTELY REQUIRED!
    !!CALL end_message_bi(modstr,'MEMORY ALLOCATION',substr)  ! log-output

    ! CHANNEL AND CHANNEL OBJECTS
    CALL start_message_bi(modstr,'CHANNEL DEFINITION',substr)  ! log-output

    ! new channel
    CALL new_channel(status, modstr, reprid=GP_3D_MID)
    CALL channel_halt(substr, status)

!!$    ! object with attributes
!!$    CALL new_channel_object(status, modstr//'_gp', 'objname', &
!!$         p3=objptr)
!!$    CALL channel_halt(substr, status)
!!$    CALL new_attribute(status, modstr//'_gp', 'objname', &
!!$         'long_name', c='explanation')
!!$    CALL channel_halt(substr, status)
!!$    CALL new_attribute(status, modstr//'_gp', 'objname', &
!!$         'units', c=' ')
!!$    CALL channel_halt(substr, status)
!!$    CALL info_bi('channel/object '//modstr//'_gp/'// &
!!$         'objname'//' was created')
!!$
!!$    ! ### ADD MORE CHANNEL OBJECTS HERE
!!$
    CALL end_message_bi(modstr,'CHANNEL DEFINITION',substr)  ! log-output


  END SUBROUTINE relax_init_memory
  ! ====================================================================

  ! ====================================================================
  SUBROUTINE relax_init_coupling

    ! ------------------------------------------------------------------
    ! This soubroutine is used to set pointers
    ! (channel objects and/or tracers) for coupling to the
    ! basemodel and to other submodels.
    ! ------------------------------------------------------------------

    ! MESSy BASEMODEL INTERFACE LAYER (BMIL)
    USE messy_main_channel_error_bi, ONLY: channel_halt
    USE messy_main_channel_bi,       ONLY: GP_3D_MID
    USE messy_main_channel,          ONLY: get_channel_object, new_attribute &
                                         , new_channel_object
!                                       new_channel ! not used in this subroutine


    IMPLICIT NONE
    INTRINSIC :: TRIM

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'relax_init_coupling'
    INTEGER                     :: status


    CALL start_message_bi(modstr,'COUPLING',substr)  ! log-output


    IF (lrayfr) THEN
    ! ### check if equ temperture input is from channel or fct/const, and define new channel / get channel
        IF ( (TRIM(rayfr_k_inp%cha) .EQ. '#const') .OR. &
               (TRIM(rayfr_k_inp%cha) .EQ. '#fct') ) THEN

             CALL new_channel_object(status, modstr, 'kdamp' &
                  ,   p3=kdamp , reprid=GP_3D_MID)
             CALL channel_halt(substr, status)
             CALL new_attribute(status, modstr, 'kdamp'   &
                                , 'long_name', c='wind damping coefficient' )
             CALL channel_halt(substr, status)
             CALL new_attribute(status, modstr, 'kdamp'   &
                               , 'units', c='s-1' )
             CALL channel_halt(substr, status)

          ELSE

             CALL get_channel_object(status,TRIM(rayfr_k_inp%cha) &
                  , TRIM(rayfr_k_inp%obj), p3 = kdamp)
             CALL  channel_halt(substr, status)

        END IF
     END IF

   IF (lnewco) THEN
   ! ### check if equ temperture input is from channel or fct/const, and define new channel / get channel
        IF ( (TRIM(newco_t_inp%cha) .EQ. '#const') .OR. &
               (TRIM(newco_t_inp%cha) .EQ. '#fct') ) THEN

             CALL new_channel_object(status, modstr, 'tequ' &
                  ,   p3=tequ , reprid=GP_3D_MID)
             CALL channel_halt(substr, status)
             CALL new_attribute(status, modstr, 'tequ'   &
                                , 'long_name', c='equilibrium temperature' )
             CALL channel_halt(substr, status)
             CALL new_attribute(status, modstr, 'tequ'   &
                               , 'units', c='K' )
             CALL channel_halt(substr, status)

          ELSE

             CALL get_channel_object(status,TRIM(newco_t_inp%cha) &
                  , TRIM(newco_t_inp%obj), p3 = tequ)
             CALL  channel_halt(substr, status)

        END IF


    ! ### check if kappa input is from channel or fct/const, and define new channel / get channel
        IF ( (TRIM(newco_k_inp%cha) .EQ. '#const') .OR. &
               (TRIM(newco_k_inp%cha) .EQ. '#fct') ) THEN


             CALL new_channel_object(status, modstr, 'kappa' &
                  ,   p3=kappa , reprid=GP_3D_MID)
             CALL channel_halt(substr, status)
             CALL new_attribute(status, modstr, 'kappa'   &
                                , 'long_name', c='inverse relaxation time scale' )
             CALL channel_halt(substr, status)
             CALL new_attribute(status, modstr, 'kappa'   &
                               , 'units', c='1/s' )
             CALL channel_halt(substr, status)

          ELSE

             CALL get_channel_object(status,TRIM(newco_k_inp%cha) &
                  , TRIM(newco_k_inp%obj), p3 = kappa)
             CALL  channel_halt(substr, status)

        END IF
    END IF


    ! +op_rw_20190212
    IF (liheat_cc_tropics) THEN
    ! ### check if tteh_cc_tropics input is from channel or fct/const, and define new channel / get channel
        IF ( (TRIM(cct_h_inp%cha) .EQ. '#const') .OR. &
               (TRIM(cct_h_inp%cha) .EQ. '#fct') ) THEN

             CALL new_channel_object(status, modstr, 'tteh_cc_tropics' &
                  ,   p3=tteh_cc_tropics , reprid=GP_3D_MID)
             CALL channel_halt(substr, status)
             CALL new_attribute(status, modstr, 'tteh_cc_tropics'   &
                                , 'long_name', c='temperature tendency due to climate change-like idealized heating' )
             CALL channel_halt(substr, status)
             CALL new_attribute(status, modstr, 'tteh_cc_tropics'   &
                               , 'units', c='K/s' )
             CALL channel_halt(substr, status)

          ELSE

             CALL get_channel_object(status,TRIM(cct_h_inp%cha) &
                  , TRIM(cct_h_inp%obj), p3 = tteh_cc_tropics)
             CALL channel_halt(substr, status)

        END IF
     END IF
    ! -op_rw_20190212


    ! +op_rw_20181102
    IF (liheat_waves) THEN
    ! ### check if tteh_waves  input is from channel or fct/const, and define new channel / get channel
        IF ( (TRIM(waves_h_inp%cha) .EQ. '#const') .OR. &
               (TRIM(waves_h_inp%cha) .EQ. '#fct') ) THEN

             CALL new_channel_object(status, modstr, 'tteh_waves' &
                  ,   p3=tteh_waves , reprid=GP_3D_MID)
             CALL channel_halt(substr, status)
             CALL new_attribute(status, modstr, 'tteh_waves'   &
                                , 'long_name', c='temperature tendency due to idealized heating for planetary wave generation' )
             CALL channel_halt(substr, status)
             CALL new_attribute(status, modstr, 'tteh_waves'   &
                               , 'units', c='K/s' )
             CALL channel_halt(substr, status)

          ELSE

             CALL get_channel_object(status,TRIM(waves_h_inp%cha) &
                  , TRIM(waves_h_inp%obj), p3 = tteh_waves)
             CALL channel_halt(substr, status)

        END IF
     END IF
    ! -op_rw_20181102


    ! +op_mn_20180620
    IF (liheat_mons) THEN
       ! ### check if tteh_mons input is from channel or fct/const, and define new channel / get channel
        IF ( (TRIM(reght_h_inp%cha) .EQ. '#const') .OR. &
               (TRIM(reght_h_inp%cha) .EQ. '#fct') ) THEN

             CALL new_channel_object(status, modstr, 'tteh_mons' &
                  ,   p3=tteh_mons , reprid=GP_3D_MID)
             CALL channel_halt(substr, status)
             CALL new_attribute(status, modstr, 'tteh_mons'   &
                                , 'long_name', c='temperature tendency due to idealized heating for monsoon' )
             CALL channel_halt(substr, status)
             CALL new_attribute(status, modstr, 'tteh_mons'   &
                               , 'units', c='K/s' )
             CALL channel_halt(substr, status)

          ELSE

             CALL get_channel_object(status,TRIM(reght_h_inp%cha) &
                  , TRIM(reght_h_inp%obj), p3 = tteh_mons)
             CALL channel_halt(substr, status)

         END IF

    END IF
    ! -op_mn_20180620



    CALL end_message_bi(modstr,'COUPLING',substr)  ! log-output


  END SUBROUTINE relax_init_coupling
  ! ====================================================================


  ! ====================================================================
  SUBROUTINE relax_physc

    ! ------------------------------------------------------------------
    ! This subroutine is called within the time loop.
    ! It constitutes the main entry point for additional processes
    ! or diagnostics.
    ! Here, only the current vector of the grid-point-fields is
    ! accessible.
    !
    ! calculates wind tendency due to rayleigh friction (instead of vdiff and/or upper sponge, thus replacing vdiff)
    ! and / or temperature tendency due to newtonian cooling (instead of radiative heating, thus replacing radheat)
    ! calculates additionally (if adjusted in relax-namelist):
    !        - idealized heating for climate change-like tropical upper-tropospheric warming
    !        - idealized heating for planetary wave generation
    !        - idealized heating for monsoon
    !
    ! called from messy_physc
    ! ------------------------------------------------------------------

    ! MESSy BASEMODEL INTERFACE LAYER (BMIL)
    USE messy_main_timer,         ONLY: delta_time, current_time_step
    USE messy_main_grid_def_mem_bi, ONLY: kproma, jrow, nlev
    USE messy_main_grid_def_bi,     ONLY: philon_2d, philat_2d &
                                        , hyam, hybm
    USE messy_main_data_bi,       ONLY: um1, vm1, apm1 &
#ifndef MESSYTENDENCY
                                      , vol_3d, vom_3d, tte_3d     &
#endif
                                      , um1, vm1, apm1             &
                                      , tm1, aps



    USE messy_main_tools,         ONLY: str2num, strcrack
    USE messy_main_constants_mem, ONLY: STRLEN_MEDIUM !op_rw_20180620 test if STRLEN_ULONG is still needed

    IMPLICIT NONE
!    INTRINSIC :: SIZE
    INTRINSIC :: REAL, TRIM

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER         :: substr = 'relax_physc'
    INTEGER                             :: status
    REAL(DP)                            :: const_val
    REAL(DP), DIMENSION(kproma,nlev)    :: my_vol, my_vom, my_tte      ! perturbation tendency
    REAL(DP), DIMENSION(kproma,nlev)    :: my_kdamp, my_tequ           ! local array for damping coeff. and modified tequ
!    CHARACTER(LEN=STRLEN_ULONG), POINTER :: outstring(:)       => NULL() ! op_rw_20180620 test if STRLEN_ULONG is still needed
    CHARACTER(LEN=STRLEN_MEDIUM), POINTER :: outstring(:)       => NULL()
! +op_rw_20190212: one outstring is sufficient
!    CHARACTER(LEN=STRLEN_MEDIUM), POINTER :: outstring1(:)      => NULL() ! op_mn_20180620
!    CHARACTER(LEN=STRLEN_MEDIUM), POINTER :: outstring2(:)      => NULL() ! op_mn_20180620
!    CHARACTER(LEN=STRLEN_MEDIUM), POINTER :: outstring_waves(:) => NULL() ! op_rw_20181102
! -op_rw_20190212

    INTEGER                             :: m, im, i_in

    ! DEFAULT VALUES
    ! ---Newtonian cooling-------------------------------------------------------------------------
    ! 1.) HS set-up for tequ
    ! name the parameters of tequ
    CHARACTER(LEN=*), DIMENSION(7), PARAMETER :: HS_tequ_parameter_name = (/'hfac   ', 'p0     ', 'T0     ' &
                                                                 , 'T1     ', 'Ty     ', 'Tz     ', 'eps_abs'/)
    ! define default values for parameters of tequ
    REAL(DP), DIMENSION(7), PARAMETER :: HS_tequ_parameter_default = (/0.0_dp, 101325.0_dp, 200.0_dp &
                                                                 , 315.0_dp, 60.0_dp, 10.0_dp, 0.0_dp/)
    ! in this variable the further on used parameter values are written
    REAL(DP), DIMENSION(7)            :: HS_tequ_parameter_value
    ! -------------------------------------------------------------------------
    ! 2.) PK set-up for tequ
    ! name the parameters of tequ
    CHARACTER(LEN=*), DIMENSION(11), PARAMETER :: PK_tequ_parameter_name = (/'gammaPK' &
                                                                 , 'hfac   ', 'p0     ', 'T1     ', 'Ty     ' &
                                                                 , 'Tz     ', 'eps_abs', 'l0_abs ', 'dl     ' &
                                                                 , 'pT_SH  ', 'pT_WH  '/)
    ! define default values for parameters of tequ
    REAL(DP), DIMENSION(11), PARAMETER :: PK_tequ_parameter_default = (/4.0_dp &
                                                                 , 1.0_dp, 101325.0_dp, 315.0_dp, 60.0_dp &
                                                                 , 10.0_dp, 10.0_dp, 50.0_dp, 10.0_dp &
                                                                 , 10000.0_dp, 10000.0_dp/)
    ! in this variable the further on used parameter values are written
    REAL(DP), DIMENSION(11)            :: PK_tequ_parameter_value
    ! -------------------------------------------------------------------------
    ! 3.) HS set-up for kappa which is the inverse relax time for both tequ set-ups
    ! name the parameters of kappa
    CHARACTER(LEN=*), DIMENSION(3), PARAMETER :: HS_kappa_parameter_name = (/'ta     ', 'ts     ', 'sigb   '/)
    ! define default values for parameters of kappa
    REAL(DP), DIMENSION(3), PARAMETER :: HS_kappa_parameter_default = (/40.0_dp , 4.0_dp, 0.7_dp/)
    ! in this variable the further on used parameter values are written
    REAL(DP), DIMENSION(3)            :: HS_kappa_parameter_value
    ! ---------------------------------------------------------------------------------------------

    ! ---Rayleigh friction-------------------------------------------------------------------------
    ! 1.) HS set-up for wind damping close to surface
    ! name the parameters of kdamp
    CHARACTER(LEN=*), DIMENSION(2), PARAMETER :: HS_kdamp_parameter_name = (/'kmaxHS ', 'sig0   '/)
    ! define default values for parameters of kdamp
    REAL(DP), DIMENSION(2), PARAMETER :: HS_kdamp_parameter_default = (/1.1574e-05_dp, 0.7_dp/)
    ! in this variable the further on used parameter values are written
    REAL(DP), DIMENSION(2)            :: HS_kdamp_parameter_value
    ! -------------------------------------------------------------------------
    ! 2.) PK set-up for wind damping at model top
    ! name the parameters of kdamp
    CHARACTER(LEN=*), DIMENSION(2), PARAMETER :: PK_kdamp_parameter_name = (/'kmaxPK ', 'psp    '/)
    ! define default values for parameters of kdamp
    REAL(DP), DIMENSION(2), PARAMETER :: PK_kdamp_parameter_default = (/2.3148e-05_dp, 50.0_dp/)
    ! in this variable the further on used parameter values are written
    REAL(DP), DIMENSION(2)            :: PK_kdamp_parameter_value
    ! -------------------------------------------------------------------------
    ! 3.) EH set-up for wind damping at model top
    ! name the parameters of kdamp
    CHARACTER(LEN=*), DIMENSION(2), PARAMETER :: EH_kdamp_parameter_name = (/'spdrag ', 'enfac  '/)
    !CHARACTER(LEN=*), DIMENSION(3), PARAMETER :: EH_kdamp_parameter_name = (/'spdrag ', 'enfac  ', 'nlevs  '/)
    ! define default values for parameters of kdamp
    REAL(DP), DIMENSION(2), PARAMETER :: EH_kdamp_parameter_default = (/5.0200e-07_dp, 1.5238_dp/)
    !INTEGER,                PARAMETER :: EH_kdamp_nlevs_default     = 10

    ! in this variable the further on used parameter values are written
    REAL(DP), DIMENSION(2)            :: EH_kdamp_parameter_value
    !INTEGER                           :: EH_kdamp_nlevs
    ! ---------------------------------------------------------------------------------------------


    ! ---Climate change-like idealized heating for tropical upper-tropospheric warming-------------
    ! name of the parameters
    CHARACTER(LEN=*), DIMENSION(5), PARAMETER :: tteh_cc_tropics_parameter_name = (/'q0_cct   ', 'lat0     ' &
                                                                                  , 'sigma_lat', 'z0       ' &
                                                                                  , 'sigma_z  '/)
    ! default values
    REAL(DP), DIMENSION(5), PARAMETER :: tteh_cc_tropics_parameter_default = (/0.5_dp, 0.0_dp &
                                                                             , 0.4_dp, 0.3_dp &
                                                                             , 0.11_dp/)
    ! further used parameter values
    REAL(DP), DIMENSION(5)            :: tteh_cc_tropics_parameter_value
    ! ---------------------------------------------------------------------------------------------


    ! ---Idealized heating planetary wave generation-----------------------------------------------
    ! name of the parameters
    CHARACTER(LEN=*), DIMENSION(6), PARAMETER :: tteh_waves_parameter_name = (/'q0       ', 'm_WN     ' &
                                                                             , 'phi0     ', 'sigma_phi' &
                                                                             , 'p_bot    ', 'p_top    '/)
    ! default values
    REAL(DP), DIMENSION(6), PARAMETER :: tteh_waves_parameter_default = (/6.0_dp, 2.0_dp &
                                                                       , 45.0_dp, 0.175_dp &
                                                                       , 80000.0_dp, 20000.0_dp/)
    ! further used parameter values
    REAL(DP), DIMENSION(6)            :: tteh_waves_parameter_value
    ! ---------------------------------------------------------------------------------------------


    ! ---Idealized heating monsoon-----------------------------------------------------------------
    REAL(DP)                                :: pres0, lat0, lon0, presd, latd, lond ! op_mn_20180620
    REAL(DP)                                :: ofht, amht, peht, suht               ! op_mn_20180620
    !----------------------------------------------------------------------------------------------


    !INTRINSIC                           :: minval, maxval ! op_mn_20180620 op_rw_20180620 seems not to be needed anymore

    ! init my_vol/m
    my_vol(:,:) = 0.0_dp
    my_vom(:,:) = 0.0_dp
    my_tte(:,:) = 0.0_dp

    ! ---------------------------------------------------------------------------------------------
    IF (lrayfr) THEN
        ! init kdamp
        kdamp(1:kproma,:,jrow) = 0.0_dp
        ! then add damping coefficients

        ! set kdamp (if const or fct)
        IF (TRIM(rayfr_k_inp%cha) .EQ. '#const') THEN

            CALL str2num(rayfr_k_inp%obj, const_val, status)
            kdamp(:,:,:) = const_val

        ELSEIF (TRIM(rayfr_k_inp%cha) .EQ. '#fct') THEN

            CALL strcrack(TRIM(rayfr_k_inp%obj),',',outstring,m,.TRUE.)
            !print *, '#*# File string with input parameters: ', outstring
            !print *, '#*# Number of parameters: ', m

            IF ( (m /= 0) .AND. (m /= 3) .AND. (m /= 6) .AND. (m /= 9) ) THEN

                CALL error_bi('syntax error in #obj declaration: input has to be composed of &
                                     &"HS,kmaxHS,sig0", "PK,kmaxPK,psp", "EH,spdrag,enfac" ',substr)

            END IF

            ! -------------------------------------------------------------------------------------

            im = 1

            DO WHILE ( im .LE. m )

                SELECT CASE(TRIM(outstring(im)))

                    CASE('HS')

                        ! set parameter values for rayleigh friction at bottom layers
                        DO i_in = 1, 2

                            IF ( TRIM(outstring(im + i_in)) == '') THEN
                                !print *, '#*# case of empty string: use default value'
                                HS_kdamp_parameter_value(i_in) = HS_kdamp_parameter_default(i_in)
                            ELSE
                                !print *, '#*# case of non empty string: use value of rayfr_k_inp in relax.nml'
                                CALL str2num(outstring(im + i_in), HS_kdamp_parameter_value(i_in), status)
                            END IF

                            !print *, '#*# use ', HS_kdamp_parameter_name(i_in), ' = ', HS_kdamp_parameter_value(i_in)

                        END DO

                        CALL relax_kdamphs(&
                            my_kdamp(1:kproma,:) &
                          , apm1(1:kproma,:) &
                          , aps(1:kproma,jrow) &
                          , HS_kdamp_parameter_value(1) &
                          , HS_kdamp_parameter_value(2) &
                          )

                        kdamp(1:kproma,:,jrow) = kdamp(1:kproma,:,jrow) + my_kdamp

                    CASE('PK')

                        ! set parameter values for rayleigh friction at top layers
                        DO i_in = 1, 2

                            IF ( TRIM(outstring(im + i_in)) == '') THEN
                                !print *, '#*# case of empty string: use default value'
                                PK_kdamp_parameter_value(i_in) = PK_kdamp_parameter_default(i_in)
                            ELSE
                                !print *, '#*# case of non empty string: use value of rayfr_k_inp in relax.nml'
                                CALL str2num(outstring(im + i_in), PK_kdamp_parameter_value(i_in), status)
                            END IF

                            !print *, '#*# use ', PK_kdamp_parameter_name(i_in), ' = ', PK_kdamp_parameter_value(i_in)

                        END DO

                        CALL relax_kdamppk(&
                            my_kdamp(1:kproma,:) &
                          , apm1(1:kproma,:) &
                          , PK_kdamp_parameter_value(1) &
                          , PK_kdamp_parameter_value(2) &
                          )

                        kdamp(1:kproma,:,jrow) = kdamp(1:kproma,:,jrow) + my_kdamp



                    CASE('EH')

                        ! set parameter values for rayleigh friction at top layers
                        DO i_in = 1, 2

                            IF ( TRIM(outstring(im + i_in)) == '') THEN
                                !print *, '#*# case of empty string: use default value'
                                EH_kdamp_parameter_value(i_in) = EH_kdamp_parameter_default(i_in)
                            ELSE
                                !print *, '#*# case of non empty string: use value of rayfr_k_inp in relax.nml'
                                CALL str2num(outstring(im + i_in), EH_kdamp_parameter_value(i_in), status)
                            END IF

                            !print *, '#*# use ', EH_kdamp_parameter_name(i_in), ' = ', EH_kdamp_parameter_value(i_in)

                        END DO

                        ! +op_rw_20190927
                        ! set nlev
                        !IF ( TRIM(outstring(7)) == '') THEN
                        !        !print *, '#*# case of empty string: use default value'
                        !        EH_kdamp_nlevs = EH_kdamp_nlevs_default
                        !    ELSE
                        !        !print *, '#*# case of non empty string: use value of rayfr_k_inp in relax.nml'
                        !        CALL str2num(outstring(7), EH_kdamp_nlevs, status)
                        !    END IF
                        !
                        !    !print *, '#*# use ', EH_kdamp_parameter_name(3), ' = ', EH_kdamp_nlevs
                        ! -op_rw_20190927

                        CALL relax_kdampeh(&
                            my_kdamp(1:kproma,:) &
                          , apm1(1:kproma,:) &
                          , EH_kdamp_parameter_value(1) &
                          , EH_kdamp_parameter_value(2) &
                          !, EH_kdamp_nlevs &
                          )

                        kdamp(1:kproma,:,jrow) = kdamp(1:kproma,:,jrow) + my_kdamp


                    CASE DEFAULT

                        CALL error_bi('choosen function name '//&
                            &TRIM(rayfr_k_inp%obj)//&
                            ' not supported!', &
                            substr)

                END SELECT

                ! every sponge label (HS,PK,EH) is followed by two parameters
                im = im + 3

            END DO

            DEALLOCATE(outstring) ; NULLIFY(outstring)

        END IF ! (if const or fct)

        ! call SMCL routine to calculate horizontal wind tendendy
        CALL relax_rayfr_smcl(um1(1:kproma,:,jrow), vm1(1:kproma,:,jrow), kdamp(1:kproma,:,jrow)  &
        , my_vom(:,:), my_vol(:,:))

#ifdef MESSYTENDENCY
        !tendency budget
        call mtend_add_l (my_handle, mtend_id_u &
        , px  = my_vom)

        call mtend_add_l (my_handle, mtend_id_v &
        , px  = my_vol)

#else
        ! add tendency
        vom_3d(1:kproma,1:nlev,jrow) = vom_3d(1:kproma,1:nlev,jrow) &
                               + my_vom(:,:)  ! -> tendency [m/s/s]

        vol_3d(1:kproma,1:nlev,jrow) = vol_3d(1:kproma,1:nlev,jrow) &
                               + my_vol(:,:)  ! -> tendency [m/s/s]

#endif
    END IF ! (lrayfr)
    ! ---------------------------------------------------------------------------------------------

    ! ---------------------------------------------------------------------------------------------
    IF (lnewco) THEN
        ! init tequ and kappa
        tequ(1:kproma,:,jrow)  = 0.0_dp
        kappa(1:kproma,:,jrow) = 0.0_dp

        ! set Tequ and kappa (if const or fct)
        IF (TRIM(newco_t_inp%cha) .EQ. '#const') THEN

            CALL str2num(newco_t_inp%obj, const_val, status)
            tequ(:,:,:) = const_val

        ELSEIF (TRIM(newco_t_inp%cha) .EQ. '#fct') THEN

            ! crack input file string with (m-1) commas in m parts
            CALL strcrack(TRIM(newco_t_inp%obj),',',outstring,m,.TRUE.)
            !print *, '#*# File string with input parameters: ', outstring
            !print *, '#*# Number of parameters: ', m

            ! length of outstring(i) is 24 due to STRLEN_MEDIUM
            ! print *, '#*# ', LEN(outstring(2))

            SELECT CASE(TRIM(outstring(1)))

                CASE('HS')

                    IF (m /= 8) THEN
                        CALL error_bi('syntax error in #obj declaration: HS need input of form &
                                     &"HS,hfac,p0,T0,T1,Ty,Tz,eps_abs" ',substr)
                    END IF

                    ! set parameter value for newtonian cooling
                    DO i_in = 1, 7

                        IF ( TRIM(outstring(i_in + 1)) == '') THEN
                            !print *, '#*# case of empty string: use default value'
                            HS_tequ_parameter_value(i_in) = HS_tequ_parameter_default(i_in)
                        ELSE
                            !print *, '#*# case of non empty string: use value of newco_t_inp in relax.nml'
                            CALL str2num(outstring(i_in + 1),  HS_tequ_parameter_value(i_in), status)
                        END IF

                        !print *, '#*# use ', HS_tequ_parameter_name(i_in), ' = ', HS_tequ_parameter_value(i_in)

                    END DO

                    CALL relax_tequhs(&
                           tequ(1:kproma,:,jrow) &
                         , philat_2d(1:kproma,jrow) &
                         , apm1(1:kproma,:) &
                         , HS_tequ_parameter_value(1) &
                         , HS_tequ_parameter_value(2) &
                         , HS_tequ_parameter_value(3) &
                         , HS_tequ_parameter_value(4) &
                         , HS_tequ_parameter_value(5) &
                         , HS_tequ_parameter_value(6) &
                         , HS_tequ_parameter_value(7) &
                         )

                CASE('PK')

                    IF (m /=12) THEN
                        CALL error_bi('syntax error in #obj declaration: PK need input of form &
                                     &"PK,gamma,hfac,p0,T1,Ty,Tz,eps_abs,l0_abs,dl,pT_SH,pT_WH" ',substr)
                    END IF

                    ! set parameter value for newtonian cooling
                    DO i_in = 1, 11

                        IF ( TRIM(outstring(i_in + 1)) == '') THEN
                            !print *, '#*# case of empty string: use default value'
                            PK_tequ_parameter_value(i_in) = PK_tequ_parameter_default(i_in)
                        ELSE
                            !print *, '#*# case of non empty string: use value of newco_t_inp in relax.nml'
                            CALL str2num(outstring(i_in + 1),  PK_tequ_parameter_value(i_in), status)
                        END IF

                        !print *, '#*# use ', PK_tequ_parameter_name(i_in), ' = ', PK_tequ_parameter_value(i_in)

                    END DO


                    CALL relax_tequpk(&
                           tequ(1:kproma,:,jrow) &
                         , philat_2d(1:kproma,jrow) &
                         , apm1(1:kproma,:) &
                         , PK_tequ_parameter_value(1) &
                         , PK_tequ_parameter_value(2) &
                         , PK_tequ_parameter_value(3) &
                         , PK_tequ_parameter_value(4) &
                         , PK_tequ_parameter_value(5) &
                         , PK_tequ_parameter_value(6) &
                         , PK_tequ_parameter_value(7) &
                         , PK_tequ_parameter_value(8) &
                         , PK_tequ_parameter_value(9) &
                         , PK_tequ_parameter_value(10) &
                         , PK_tequ_parameter_value(11) &
                         , l_no_polar_vortex &
                         )

                CASE DEFAULT

                    CALL error_bi('choosen function name '//&
                    &TRIM(newco_t_inp%obj)//&
                    ' not supported!', &
                    substr)
            END SELECT

            DEALLOCATE(outstring) ; NULLIFY(outstring)

        ELSE !if imported channel object

            ! interpolate imported tequ to current pressure profile
            CALL relax_intpol_p(tequ(1:kproma,:,jrow),apm1(1:kproma,:),hyam,hybm,my_tequ(1:kproma,:))
        END IF
        ! -----------------------------------------------------------------------------------------

        ! -----------------------------------------------------------------------------------------
        IF (TRIM(newco_k_inp%cha) .EQ. '#const') THEN

            CALL str2num(newco_k_inp%obj, const_val, status)
            kappa(:,:,:) = const_val

        ELSEIF (TRIM(newco_k_inp%cha) .EQ. '#fct') THEN

            ! crack input file string with (m-1) commas in m parts
            CALL strcrack(TRIM(newco_k_inp%obj),',',outstring,m,.TRUE.)
            !print *, '#*# File string with input parameters: ', outstring
            !print *, '#*# Number of parameters: ', m, ' need 4'

            SELECT CASE(TRIM(outstring(1)))

                ! At the moment, there is only one function (relax_kappahs) to calculate kappa!; op_rw_20181105
                CASE('HS')

                    IF (m /= 4) THEN
                        CALL error_bi('syntax error in #obj declaration: HS need input of form &
                                     &"HS,ta,ts,sigb" ',substr)
                    END IF

                    ! set parameter value for newtonian cooling
                    DO i_in = 1, 3

                        IF ( TRIM(outstring(i_in + 1)) == '') THEN
                            !print *, '#*# case of empty string: use default value'
                            HS_kappa_parameter_value(i_in) = HS_kappa_parameter_default(i_in)
                        ELSE
                            !print *, '#*# case of non empty string: use value of newco_k_inp in relax.nml'
                            CALL str2num(outstring(i_in + 1),  HS_kappa_parameter_value(i_in), status)
                        END IF

                        !print *, '#*# use ', HS_kappa_parameter_name(i_in), ' = ', HS_kappa_parameter_value(i_in)

                    END DO

                    CALL relax_kappahs(&
                           kappa(1:kproma,:,jrow) &
                         , philat_2d(1:kproma,jrow) &
                         , apm1(1:kproma,:) &
                         , aps(1:kproma,jrow) &
                         , HS_kappa_parameter_value(1) &
                         , HS_kappa_parameter_value(2) &
                         , HS_kappa_parameter_value(3) &
                         )

                CASE('PK')

                    IF (m /= 4) THEN
                        CALL error_bi('syntax error in #obj declaration: PK need input of form &
                                     &"PK,ta,ts,sigb" ',substr)
                    END IF

                    ! set parameter value for newtonian cooling
                    DO i_in = 1, 3

                        IF ( TRIM(outstring(i_in + 1)) == '') THEN
                            !print *, '#*# case of empty string: use default value'
                            HS_kappa_parameter_value(i_in) = HS_kappa_parameter_default(i_in)
                        ELSE
                            !print *, '#*# case of non empty string: use value of newco_k_inp in relax.nml'
                            CALL str2num(outstring(i_in + 1),  HS_kappa_parameter_value(i_in), status)
                        END IF

                        !print *, '#*# use ', HS_kappa_parameter_name(i_in), ' = ', HS_kappa_parameter_value(i_in)

                    END DO

                    CALL relax_kappahs(&
                           kappa(1:kproma,:,jrow) &
                         , philat_2d(1:kproma,jrow) &
                         , apm1(1:kproma,:) &
                         , aps(1:kproma,jrow) &
                         , HS_kappa_parameter_value(1) &
                         , HS_kappa_parameter_value(2) &
                         , HS_kappa_parameter_value(3) &
                         )

                CASE DEFAULT
                    CALL error_bi('choosen function name '//&
                    &TRIM(newco_k_inp%obj)//&
                    ' not supported!', &
                    substr)

            END SELECT

            DEALLOCATE(outstring) ; NULLIFY(outstring)

        END IF

        ! call SMCL routine to calculate temperature tendency
        IF (TRIM(newco_t_inp%cha) .EQ. '#const' .OR. TRIM(newco_t_inp%cha) .EQ. '#fct' ) THEN

            CALL relax_newco_smcl(tm1(1:kproma,:,jrow), tequ(1:kproma,:,jrow), kappa(1:kproma,:,jrow), my_tte(:,:))

        ELSE

            CALL relax_newco_smcl(tm1(1:kproma,:,jrow), my_tequ(1:kproma,:),kappa(1:kproma,:,jrow), my_tte(:,:))

        END IF

! +op_mn_20180620
! This part is shifted to the end of the subroutine
!#ifdef MESSYTENDENCY
!        !tendency budget
!        call mtend_add_l (my_handle, mtend_id_t &
!                          , px  = my_tte)
!#else
!        ! add tendency
!        tte(1:kproma,1:nlev) = tte(1:kproma,1:nlev) &
!                                + my_tte(:,:)  ! -> tendency [K/s]
!#endif
! -op_mn_20180620
    END IF !endif IF (lnewco)
    ! ---------------------------------------------------------------------------------------------


    ! ---------------------------------------------------------------------------------------------
    ! +op_rw_20190212
    IF (liheat_cc_tropics) THEN
        ! init tteh_cc_tropics
        tteh_cc_tropics(1:kproma,:,jrow) = 0.0_dp

        IF (TRIM(cct_h_inp%cha) .EQ. '#const') THEN

            CALL error_bi('error in #obj declaration: #const not defined for cct_h_inp%cha',substr)

        ELSEIF (TRIM(cct_h_inp%cha) .EQ. '#fct') THEN

            CALL strcrack(TRIM(cct_h_inp%obj),',',outstring,m,.TRUE.)
            !print *, '#*# File string with input parameters: ', outstring
            !print *, '#*# Number of parameters: ', m, ' need 5'

            IF (m /= 5) THEN

                CALL error_bi('syntax error in #obj declaration: tteh_cc_tropics need input of form &
                               &"q0_cct,lat0,sigma_lat,z0,sigma_z" ',substr)

            END IF

            ! set parameter value for climate change-like idealized heating for tropical upper-tropospheric warming
            DO i_in = 1, 5

                IF ( TRIM(outstring(i_in)) == '' ) THEN
                    !print *, '#*# case of empty string: use default value'
                    tteh_cc_tropics_parameter_value(i_in) = tteh_cc_tropics_parameter_default(i_in)
                ELSE
                    !print *, '#*# case of non empty string: use value of cct_h_inp in relax.nml'
                    CALL str2num(outstring(i_in), tteh_cc_tropics_parameter_value(i_in), status)
                END IF

                !print *, '#*# use ', tteh_cc_tropics_parameter_name(i_in), ' = ', tteh_cc_tropics_parameter_value(i_in)

            END DO

            DEALLOCATE(outstring) ; NULLIFY(outstring)

            ! call SMCL to calculate tte due to climate change-like idealized heating for tropical upper-tropospheric warming
            CALL relax_tteh_cc_tropics_smcl(apm1(1:kproma,:) &
                                          , aps(1:kproma,jrow) &
                                          , philat_2d(1:kproma,jrow) &
                                          , tteh_cc_tropics_parameter_value(1) &
                                          , tteh_cc_tropics_parameter_value(2) &
                                          , tteh_cc_tropics_parameter_value(3) &
                                          , tteh_cc_tropics_parameter_value(4) &
                                          , tteh_cc_tropics_parameter_value(5) &
                                          , l_Butler_heat &
                                          , tteh_cc_tropics(1:kproma,:,jrow))

        END IF

    END IF ! IF (liheat_cc_tropics)
    ! -op_rw_20190212
    ! ---------------------------------------------------------------------------------------------


    ! ---------------------------------------------------------------------------------------------
    ! +op_rw_20181102
    IF (liheat_waves) THEN
        ! init tteh_waves
        tteh_waves(1:kproma,:,jrow) = 0.0_dp

        IF (TRIM(waves_h_inp%cha) .EQ. '#const') THEN

            CALL error_bi('error in #obj declaration: #const not defined for waves_h_inp%cha',substr)

        ELSEIF (TRIM(waves_h_inp%cha) .EQ. '#fct') THEN

            CALL strcrack(TRIM(waves_h_inp%obj),',',outstring,m,.TRUE.)
            !print *, '#*# File string with input parameters: ', outstring
            !print *, '#*# Number of parameters: ', m, ' need 6'

            IF (m /= 6) THEN

                CALL error_bi('syntax error in #obj declaration: tteh_waves need input of form &
                               &"q0,m_WN,phi0,sigma_phi,p_bot,p_top" ',substr)

            END IF

            ! set parameter value for idealized heating for planetary wave generation
            DO i_in = 1, 6

                IF ( TRIM(outstring(i_in)) == '' ) THEN
                    !print *, '#*# case of empty string: use default value'
                    tteh_waves_parameter_value(i_in) = tteh_waves_parameter_default(i_in)
                ELSE
                    !print *, '#*# case of non empty string: use value of waves_h_inp in relax.nml'
                    CALL str2num(outstring(i_in), tteh_waves_parameter_value(i_in), status)
                END IF

                !print *, '#*# use ', tteh_waves_parameter_name(i_in), ' = ', tteh_waves_parameter_value(i_in)

            END DO

            DEALLOCATE(outstring) ; NULLIFY(outstring)

            ! call SMCL to calculate tte due to heating for planetary wave generation
            CALL relax_tteh_waves_smcl(apm1(1:kproma,:), philat_2d(1:kproma,jrow), philon_2d(1:kproma,jrow) &
                                       , tteh_waves_parameter_value(1) &
                                       , tteh_waves_parameter_value(2) &
                                       , tteh_waves_parameter_value(3) &
                                       , tteh_waves_parameter_value(4) &
                                       , tteh_waves_parameter_value(5) &
                                       , tteh_waves_parameter_value(6) &
                                       , tteh_waves(1:kproma,:,jrow))

        END IF

    END IF ! IF (liheat_waves)
    ! -op_rw_20181102
    ! ---------------------------------------------------------------------------------------------


! +op_mn_20180620
    ! example from NML
    !reght_h_inp='#fct','500.00,20.0,130.0,200.0,20.0,30.0'  ! press (hPa), lat(deg) , lon, dec press (hPa), declat, declon
    !tmpht_h_inp='#fct','5.0,2.0,15.0,30.0'     ! offset (K/day), amplitude (K/day), heating period (days), spin up (days)


    IF (liheat_mons) THEN
        ! init tteh_mons
        tteh_mons(1:kproma,:,jrow) = 0.0_dp

        IF (TRIM(reght_h_inp%cha) .EQ. '#fct') THEN

            CALL strcrack(TRIM(reght_h_inp%obj),',',outstring, m)

            IF (m /=6) THEN
                CALL error_bi('syntax error in #obj declaration: heating needs input of form "pres,lat,lon,decpres,declat,declon" ' &
                     ,substr)
            END IF

            CALL str2num(outstring(1), pres0, status)
            CALL str2num(outstring(2), lat0 , status)
            CALL str2num(outstring(3), lon0 , status)
            CALL str2num(outstring(4), presd, status)
            CALL str2num(outstring(5), latd , status)
            CALL str2num(outstring(6), lond , status)

            DEALLOCATE(outstring) ; NULLIFY(outstring)
        END IF

        IF (TRIM(tmpht_h_inp%cha) .EQ. '#fct') THEN

            CALL strcrack(TRIM(tmpht_h_inp%obj),',',outstring, m)

            IF (m /=4) THEN
                CALL error_bi('syntax error in #obj declaration: heating needs input of form "ofht,amht,hperiod,spin-up" ' &
                     ,substr)
            END IF

            CALL str2num(outstring(1), ofht, status)
            CALL str2num(outstring(2), amht, status)
            CALL str2num(outstring(3), peht, status)
            CALL str2num(outstring(4), suht, status)

            DEALLOCATE(outstring) ; NULLIFY(outstring)
        END IF

        ! op_mn_20180620 call SMCL to calculate tte due to heating
        ! WHY SOME WITH JROW and some without
        CALL relax_tteh_mons_smcl(apm1(1:kproma,:), philat_2d(1:kproma,jrow), philon_2d(1:kproma,jrow) &
                               , delta_time, REAL(current_time_step,dp), pres0, lat0, lon0, presd &
                               , latd, lond, ofht, peht, amht, suht, tteh_mons(1:kproma,:,jrow))

    END IF ! endif IF liheat_mons


    ! ---Add temperature tendencies (diabatic heating for climate change and wave generation and for monsoon)----------
    ! ---to my_tte in case of active newtonian cooling-----------------------------------------------------------------
    IF (lnewco) THEN

        ! +op_rw_20190212
        IF (liheat_cc_tropics) THEN

            my_tte(:,:) = my_tte(:,:) + tteh_cc_tropics(1:kproma,1:nlev,jrow)

        END IF
        ! -op_rw_20190212

        ! +op_rw_20181102
        IF (liheat_waves) THEN

            my_tte(:,:) = my_tte(:,:) + tteh_waves(1:kproma,1:nlev,jrow)

        END IF
        ! -op_rw_20181102

        IF (liheat_mons) THEN

            my_tte(:,:) = my_tte(:,:) + tteh_mons(1:kproma,1:nlev,jrow)

        END IF

#ifdef MESSYTENDENCY
        !tendency budget
        call mtend_add_l(my_handle, mtend_id_t &
                         , px  = my_tte)
#else
        ! add tendency
        tte_3d(1:kproma,1:nlev,jrow) = tte_3d(1:kproma,1:nlev,jrow) &
                               + my_tte(:,:)  ! -> tendency [K/s]
#endif

    ! ---Add temperature tendencies directly to tte--------------------------------------------------------------------
    ELSE
        ! +op_rw_20190212
        IF (liheat_cc_tropics) THEN

#ifdef MESSYTENDENCY
        !tendency budget
        call mtend_add_l(my_handle, mtend_id_t &
                         , px  = tteh_cc_tropics(1:kproma,1:nlev,jrow))
#else
        ! add tendency
        tte_3d(1:kproma,1:nlev,jrow) = &
            tte_3d(1:kproma,1:nlev,jrow) + tteh_cc_tropics(1:kproma,1:nlev,jrow)

#endif

        END IF
        ! -op_rw_20190212

        ! +op_rw_20181102
        IF (liheat_waves) THEN

#ifdef MESSYTENDENCY
        !tendency budget
        call mtend_add_l(my_handle, mtend_id_t &
                         , px  = tteh_waves(1:kproma,1:nlev,jrow))
#else
        ! add tendency
        tte_3d(1:kproma,1:nlev,jrow) =&
             tte_3d(1:kproma,1:nlev,jrow) + tteh_waves(1:kproma,1:nlev,jrow)

#endif

        END IF
        ! -op_rw_20181102

        IF (liheat_mons) THEN

#ifdef MESSYTENDENCY
        !tendency budget
        call mtend_add_l(my_handle, mtend_id_t &
                         , px  = tteh_mons(1:kproma,1:nlev,jrow))
#else
        ! add tendency
        tte_3d(1:kproma,1:nlev,jrow) = &
             tte_3d(1:kproma,1:nlev,jrow) + tteh_mons(1:kproma,1:nlev,jrow)

#endif

        END IF

    END IF
! -op_mn_20180620





  END SUBROUTINE relax_physc
  ! ====================================================================

!!$  ! ====================================================================
!!$  SUBROUTINE bufly_free_memory
!!$
!!$    ! ------------------------------------------------------------------
!!$    ! This subroutine is used to deallocate the memory, which has
!!$    ! been "manually" allocated in bufly_init_memory.
!!$    ! Note: channel object memory must not be deallocated! This is
!!$    !       performed centrally.
!!$    ! ------------------------------------------------------------------
!!$
!!$    IMPLICIT NONE
!!$
!!$    ! LOCAL
!!$    CHARACTER(LEN=*), PARAMETER :: substr = 'bufly_free_memory'
!!$    INTEGER                     :: status
!!$
!!$  END SUBROUTINE bufly_free_memory
!!$  ! ====================================================================

  ! ####################################################################
  ! PRIVATE SUBROUTINES
  ! ####################################################################

  ! ====================================================================
  SUBROUTINE relax_read_nml_cpl(status, iou)

    ! ------------------------------------------------------------------
    ! This subroutine is used to read the CPL-namelist of the submodel.
    ! ------------------------------------------------------------------

    ! MESSy
    USE messy_main_tools,  ONLY: read_nml_open, read_nml_check, read_nml_close

    IMPLICIT NONE

    ! I/O
    INTEGER, INTENT(OUT) :: status     ! error status
    INTEGER, INTENT(IN)  :: iou        ! I/O unit


    NAMELIST /CPL/ lrayfr, lnewco, liheat_cc_tropics, liheat_waves, liheat_mons &
                  , rayfr_k_inp, newco_t_inp, l_no_polar_vortex, newco_k_inp &
                  , cct_h_inp, l_Butler_heat, waves_h_inp, reght_h_inp, tmpht_h_inp
    !op_mn_20180620 added liheat_mons and reght_h_inp, tmpht_h_inp as nml objects for monsoon-like idealized heating
    !op_rw_20181102 added liheat_waves and waves_h_inp        as nml objects for planetary wave generation
    !op_rw_20190212 added liheat_cc_tropics and cct_h_inp     as nml objects for tropical climate change
    !op_rw_20190710 added l_no_polar_vortex                   as nml object  for turn on/off polar vortex
    !op_rw_20200213 added l_Butler_heat                       as nml object  for Butler heat or log-p gaussian heat



    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr='relax_read_nml_cpl'
    LOGICAL                     :: lex      ! file exists ?
    INTEGER                     :: fstat    ! file status

    status = 1

    ! DEFAULTS
    !rayfr_k_inp%cha = ''
    !rayfr_k_inp%obj = ''

    CALL read_nml_open(lex, substr, iou, 'CPL', modstr)
    IF (.not.lex) RETURN    ! <modstr>.nml does not exist

    READ(iou, NML=CPL, IOSTAT=fstat)
    CALL read_nml_check(fstat, substr, iou, 'CPL', modstr)
    IF (fstat /= 0) RETURN  ! error while reading namelist

    ! ### ADD HERE DIAGNOSTIC OUTPUT FOR LOG-FILE

    CALL read_nml_close(substr, iou, modstr)
    status = 0 ! NO ERROR

  END SUBROUTINE relax_read_nml_cpl
  ! ====================================================================

! **********************************************************************
END MODULE messy_relax_e5
! **********************************************************************

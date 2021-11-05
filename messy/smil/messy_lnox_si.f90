#include "messy_main_ppd_bi.inc"

#ifdef COSMO
!#if defined(COSMOv4s8)
! COSMO > v4s8 :  Adjustment to humidity tracers required    
!#define DAHL2000
!#endif
#endif

! ***********************************************************************
MODULE messy_lnox_si
  ! ***********************************************************************

  ! MESSy-SMIL FOR (SUB-)SUBMODEL LNOX
  !
  ! LNOX = LIGHTNING NOx EMISSION PARAMETERIZATION (LNOx)
  !
  ! Authors: 
  ! Patrick Joeckel, MPICH,  Aug 2003
  ! Pozzer Andrea,   MPICH,  Giu 2005, Lagrangian 
  ! Holger Tost,     MPICH,  Jan 2007, added AaP parameterisations
  ! Astrid Kerkweg,  UNI-MZ, Jun 2012, added Dahl parameterisation
  ! Patrick Joeckel, DLR,    Nov 2012, revision, 2 modes of operation
  ! Patrick Joeckel, DLR,    Sep 2014, - added Finney parameterisation
  !                                    - entire submodel completely revised;
  !                                      Note: results differ numerically,
  !                                      due to usage of grid-box area:
  !                                      - old: grvol/grheight
  !                                      - new: gboxarea_2d
  ! Francisco Javier Perez-Invernon, DLR, Aug 2020: extended Finney param.
  !
  ! TODO:
  !  1) calculate scaling factors for various resolutions
  !  2) set pointers to channel objects, only if required
  !  3) clean up scheme of Dahl, 2010 (COSMO only)
  !     - update to COSMO 5.0 or above
  !     - channel output
  !     - avoid trigger, if possible
  !     - avoid MPI_ALLREDUCE (SUM) (decomposition dependent result!)
  !     - ...

  USE messy_main_blather_bi,    ONLY: start_message_bi, end_message_bi
  USE messy_main_channel,       ONLY: STRLEN_OBJECT
#ifdef MESSYTENDENCY
 USE messy_main_tendency_bi,    ONLY: mtend_get_handle,       &
                                      mtend_get_start_l,      &
                                      mtend_add_l,            &
                                      mtend_register,           &    
                                      mtend_id_tracer
#endif

  USE messy_lnox
#ifdef DAHL2000
  USE messy_main_timer_event,   ONLY: time_event, io_time_event
  USE messy_lnox_dahl2000
#endif

  IMPLICIT NONE
  INTRINSIC :: NULL
  PRIVATE
  SAVE

#ifdef MESSYTENDENCY
  INTEGER :: my_handle
#endif

  ! PARAMETERISATION SPECIFIC CHANNEL OBJECTS:
  ! Grewe: 
  ! - mean updraft velocity [m/s]
  REAL(DP),  DIMENSION(:,:), POINTER    :: muv => NULL()
  !
  ! Allen & Pickering, massflux:
  ! - (convective) updraft mass flux [kg/m^2/s]
  !   at sigma = p/ps = 0.44 (or at p=440 hPa ?)
  REAL(DP),  DIMENSION(:,:), POINTER    :: umf => NULL()
  !
  ! Allen & Pickering, precipitation:
  ! - convective precipitation at surface [kg m^-2 s^-1]
  REAL(DP),  DIMENSION(:,:), POINTER    :: precon => NULL()
  !
  ! Finney et al., cloud ice flux
  ! - updraft cloud ice flux at 440 hPa [kg m^-2 s^-1]
  REAL(DP),  DIMENSION(:,:), POINTER    :: phiice => NULL()
  !
  ! Finney et al., cloud ice flux, extended (through external iso-surface)
  ! - updraft cloud ice flux at 440 hPa [kg m^-2 s^-1]
  REAL(DP),  DIMENSION(:,:), POINTER    :: phiice_ext => NULL()

  ! POINTERS TO CHANNEL OBJECTS FOR COUPLING
  ! conv. cloud top index
  REAL(DP), DIMENSION(:,:),   POINTER  :: cu_jkt  => NULL()
  ! conv. cloud base index
  REAL(DP), DIMENSION(:,:),   POINTER  :: cu_jkb  => NULL()
  !
!!$  ! updraft velocity [m/s]
!!$  REAL(DP), DIMENSION(:,:,:), POINTER  :: cu_xupdr
!!$  ! 0 degree level index
!!$  REAL(DP), DIMENSION(:,:),   POINTER  :: cu_jkfreeze
!!$  ! SPECIAL FOR MIDLEVEL CONVECTION
!!$  ! conv. cloud top index
!!$  REAL(DP), DIMENSION(:,:),   POINTER  :: cu_jkt_mid
!!$  ! conv. cloud base index
!!$  REAL(DP), DIMENSION(:,:),   POINTER  :: cu_jkb_mid
!!$  ! 0 deg. level index
!!$  REAL(DP), DIMENSION(:,:),   POINTER  :: cu_jkfreeze_mid
  !
  ! conv. updraft massflux [kg m-^2 s^-1]
  REAL(DP), DIMENSION(:,:,:), POINTER  :: umassf 
  ! conv. precipitation [kg m-^2 s^-1]
  REAL(DP), DIMENSION(:,:,:), POINTER  :: precflx_cv 

  ! ice-flux through external iso-surface
  REAL(DP), DIMENSION(:,:), POINTER  :: fice_i
!!$  REAL(DP), DIMENSION(:,:), POINTER  :: fice_f

  ! GLOBAL NAMELIST PARAMETERS ('COUPLING' PARAMETERS)
  ! channel,object for cloud bottom level
  CHARACTER(LEN=STRLEN_OBJECT) :: c_bot(2)    = ''
  ! channel,object for cloud top level
  CHARACTER(LEN=STRLEN_OBJECT) :: c_top(2)    = ''
  ! channel,object for external iso-surface (for extended ice-flux param.)
  CHARACTER(LEN=STRLEN_OBJECT) :: c_iif(2)    = ''
  !
!!$  ! channel,object for updraft velocity
!!$  CHARACTER(LEN=STRLEN_OBJECT) :: c_updr(2)   = ''
!!$  ! channel,object for cloud freezing level
!!$  CHARACTER(LEN=STRLEN_OBJECT) :: c_freeze(2) = ''
  !
!!$  ! SPECIAL FOR MID-LEVEL CONVECTION
!!$  LOGICAL            :: l_midlevel = .FALSE.
!!$  ! channel,object for cloud bottom level
!!$  CHARACTER(LEN=STRLEN_OBJECT)  :: c_bot_mid(2)    = ''
!!$  ! channel,object for cloud top level
!!$  CHARACTER(LEN=STRLEN_OBJECT)  :: c_top_mid(2)    = ''
!!$  ! channel,object for cloud freezing level
!!$  CHARACTER(LEN=STRLEN_OBJECT)  :: c_freeze_mid(2) = ''

  ! channel,object for conv. updraft massflux
  CHARACTER(LEN=STRLEN_OBJECT)  :: c_massfu(2) = ''
  ! channel,object for conv. precipitation
  CHARACTER(LEN=STRLEN_OBJECT)  :: c_precflx(2) = ''

  ! CALCLUATE lAGRANGIAN RATE COEFFICIENTS?
  LOGICAL :: l_calc_lg = .FALSE.
  !
  ! lightning emiss [kgN/s/M3] LG
  REAL(DP), DIMENSION(:),     POINTER  :: xnox_lg
  ! l.NOx tend. [mol/mol/s] LG
  REAL(DP), DIMENSION(:),     POINTER  :: telnox_lg

  ! GLOBAL PARMETERS
  ! WHICH TRACERS SHOULD SEE LNOx ?
  INTEGER                            :: nlntrac_gp = 0 ! no of lightning NOx GP
  INTEGER                            :: nlntrac_lg = 0 ! no of lightning NOx LG
  INTEGER, DIMENSION(:), POINTER     :: idt_list_gp => NULL()
  INTEGER, DIMENSION(:), POINTER     :: idt_list_lg => NULL()

  ! pointer for channel objects:
  TYPE lnox_set
     ! cloud properties
     REAL(DP), DIMENSION(:,:),   POINTER  :: cth       => NULL()
     REAL(DP), DIMENSION(:,:),   POINTER  :: cbh       => NULL()
     REAL(DP), DIMENSION(:,:),   POINTER  :: czh       => NULL()
     ! flash frequencies, densities, ...
     REAL(DP), DIMENSION(:,:),   POINTER  :: ff        => NULL()
     REAL(DP), DIMENSION(:,:),   POINTER  :: pg        => NULL()
     REAL(DP), DIMENSION(:,:),   POINTER  :: fpscg     => NULL()
     REAL(DP), DIMENSION(:,:),   POINTER  :: fpsic     => NULL()
     REAL(DP), DIMENSION(:,:),   POINTER  :: fpsm2cg   => NULL()
     REAL(DP), DIMENSION(:,:),   POINTER  :: fpsm2ic   => NULL()
     REAL(DP), DIMENSION(:,:),   POINTER  :: npcanz    => NULL()
     ! (GP) NOx production ...
     REAL(DP), DIMENSION(:,:),   POINTER  :: NOxcg     => NULL()
     REAL(DP), DIMENSION(:,:),   POINTER  :: NOxic     => NULL()
     REAL(DP), DIMENSION(:,:,:), POINTER  :: xnox_gp   => NULL()
     REAL(DP), DIMENSION(:,:,:), POINTER  :: telnox_gp => NULL()
!!$     ! (LG)
!!$     REAL(DP), DIMENSION(:),     POINTER  :: xnox_lg   => NULL()
!!$     REAL(DP), DIMENSION(:),     POINTER  :: telnox_lg => NULL()
  END TYPE lnox_set
  TYPE(lnox_set), DIMENSION(:), POINTER  :: lnox      => NULL()
  INTEGER                                :: numbers
  INTEGER                                :: i_ff_cpl = 1 ! default

#ifdef DAHL2000
  TYPE(io_time_event) :: &
       LIGHT_IOEVENT = io_time_event(1, 'steps','first',0)
  TYPE(time_event) :: LIGHT_EVENT 
#endif

  ! SUBROUTINES/FUNCTIONS
  PUBLIC :: lnox_initialize
  PUBLIC :: lnox_init_memory
  PUBLIC :: lnox_init_coupling
  PUBLIC :: lnox_physc
  PUBLIC :: lnox_global_end
  PUBLIC :: lnox_free_memory
  !PRIVATE :: lnox_read_nml_cpl

CONTAINS

  ! ========================================================================
  SUBROUTINE lnox_initialize

    ! LIGHTNING NOx EMISSION PARAMETERIZATION INITIALIZATION
    !
    ! Author: Patrick Joeckel, MPICH, Aug 2003

    USE messy_main_blather_bi, ONLY: error_bi
    USE messy_main_mpi_bi,     ONLY: p_parallel_io, p_io, p_bcast
#ifdef DAHL2000
    USE messy_main_timer_bi,   ONLY: p_bcast_event, timer_event_init 
#endif
    USE messy_main_tools,      ONLY: find_next_free_unit

    IMPLICIT NONE

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'lnox_initialize'
    INTEGER                     :: iou    ! I/O unit
    INTEGER                     :: status ! error status

    ! INITIALIZE MAIN-CTRL
    IF (p_parallel_io) THEN
       iou = find_next_free_unit(100,200)
       CALL lnox_read_nml_ctrl(status, iou)
       IF (status /= 0) CALL error_bi(' ',substr)
    END IF
    CALL p_bcast(l_mode_scal, p_io)
    CALL p_bcast(i_ffcalc, p_io)
    CALL p_bcast(i_iccg, p_io)
    CALL p_bcast(i_shape, p_io)
    CALL p_bcast(r_scal_ff(:), p_io)
    CALL p_bcast(r_noxpf(:), p_io)
    CALL p_bcast(r_eff(:), p_io)
    CALL p_bcast(r_Grewe_A, p_io)
    CALL p_bcast(r_Grewe_B, p_io)

    ! INITIALIZE COUPLING-CONTROL
    IF (p_parallel_io) THEN
       iou = find_next_free_unit(100,200)
       CALL lnox_read_nml_cpl(status, iou)
       IF (status /= 0) CALL error_bi(' ',substr)
    END IF
!!$    CALL p_bcast(c_updr(1), p_io)
!!$    CALL p_bcast(c_updr(2), p_io)
    CALL p_bcast(c_top(1), p_io)
    CALL p_bcast(c_top(2), p_io)
    CALL p_bcast(c_bot(1), p_io)
    CALL p_bcast(c_bot(2), p_io)
    CALL p_bcast(c_iif(1), p_io)
    CALL p_bcast(c_iif(2), p_io)
!!$    CALL p_bcast(c_freeze(1), p_io)
!!$    CALL p_bcast(c_freeze(2), p_io)
    !
!!$    CALL p_bcast(l_midlevel, p_io)
!!$    CALL p_bcast(c_top_mid(1), p_io)
!!$    CALL p_bcast(c_top_mid(2), p_io)
!!$    CALL p_bcast(c_bot_mid(1), p_io)
!!$    CALL p_bcast(c_bot_mid(2), p_io)
!!$    CALL p_bcast(c_freeze_mid(1), p_io)
!!$    CALL p_bcast(c_freeze_mid(2), p_io)
    CALL p_bcast(l_calc_lg, p_io)
    CALL p_bcast(i_ff_cpl, p_io)
    CALL p_bcast(c_massfu(1),  p_io)
    CALL p_bcast(c_massfu(2),  p_io)
    CALL p_bcast(c_precflx(1), p_io)
    CALL p_bcast(c_precflx(2), p_io)

#ifdef DAHL2000
    CALL p_bcast_event(LIGHT_IOEVENT, p_io) 
    CALL timer_event_init(LIGHT_EVENT, LIGHT_IOEVENT, 'LIGHT_EVENT', 'present')
    ! Dahl, 2010
    ! qqq what about i_ffcalc == IPARAM_ALL ???
    IF (i_ffcalc == IPARAM_DahlC) CALL dahl_2010_initialize
#endif

#ifdef MESSYTENDENCY
    my_handle = mtend_get_handle(modstr)
#endif

  END SUBROUTINE lnox_initialize
  ! ========================================================================

  ! ========================================================================
  SUBROUTINE lnox_init_memory

    ! LNOx MODULE ROUTINE (ECHAM-5 INTERFACE)
    !
    ! - define LNOx specific channel(s) and allocate memory for
    !   global fields
    ! - look for NO tracer(s)
    !
    ! Author: Patrick Joeckel, MPICH, Aug 2003

    USE messy_main_channel_error_bi, ONLY: channel_halt
    USE messy_main_channel_bi,       ONLY: GP_2D_HORIZONTAL &
#if defined(ECHAM5)
    , GP_3D_MID, LG_ATTILA
#else
    , GP_3D_MID
#endif
    USE messy_main_tracer_mem_bi,   ONLY: GPTRSTR, LGTRSTR
    USE messy_main_tracer_tools_bi, ONLY: tracer_halt
    USE messy_main_mpi_bi,          ONLY: p_parallel_io
    USE messy_main_tracer,          ONLY: get_tracer_list
    USE messy_main_constants_mem,   ONLY: STRLEN_MEDIUM
    USE messy_main_channel,         ONLY: new_channel, new_channel_object &
                                        , new_attribute

    IMPLICIT NONE

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER             :: substr = 'lnox_init_memory'
    INTEGER                                 :: status
    INTEGER                                 :: ji
    INTEGER                                 :: ind
    INTEGER                                 :: jt
    CHARACTER(LEN=STRLEN_MEDIUM), DIMENSION(:), POINTER :: subnames => NULL()

    IF (i_ffcalc /= IPARAM_ALL) THEN 
       numbers   = 1
    ELSE
!qqq not yet tested
!#ifdef DAHL2000
!      numbers = 7
!#else
       numbers = 6
!#endif
    ENDIF
    ALLOCATE(lnox(numbers))

    CALL start_message_bi(modstr,'CHANNEL DEFINITION GRID POINT',substr)

    parameterisation_loop: DO ji = 1, numbers

       IF (numbers == 1) THEN
          ind = i_ffcalc
       ELSE
          ind = ji
       END IF

       ! define new channel
       CALL new_channel(status, modstr//paramnames(ind)//'_gp',&
            reprid=GP_2D_HORIZONTAL)
       CALL channel_halt(substr, status)

       IF (l_mode_scal) THEN
          CALL new_attribute(status, modstr//paramnames(ind)//'_gp', &
               TRIM(modstr)//'_mode', c='test simulation to readjust scaling')
          CALL channel_halt(substr, status)
       ELSE
          CALL new_attribute(status, modstr//paramnames(ind)//'_gp', &
               TRIM(modstr)//'_mode', c='production simulation')
          CALL channel_halt(substr, status)
       ENDIF

       CALL new_attribute(status, modstr//paramnames(ind)//'_gp', &
            TRIM(modstr)//'_r_fpm2s_min', r=fpm2s_min )
       CALL channel_halt(substr, status)
       CALL new_attribute(status, modstr//paramnames(ind)//'_gp', &
            TRIM(modstr)//'_r_fpm2s_min_unit', c='1/(m^2 s)' )
       CALL channel_halt(substr, status)

       CALL new_attribute(status, modstr//paramnames(ind)//'_gp', &
            TRIM(modstr)//'_r_scal_ff', r=r_scal_ff(ind))
       CALL channel_halt(substr, status)

       CALL new_attribute(status, modstr//paramnames(ind)//'_gp', &
            TRIM(modstr)//'_r_noxpf', r=r_noxpf(ind))
       CALL channel_halt(substr, status)

       CALL new_attribute(status, modstr//paramnames(ind)//'_gp', &
            TRIM(modstr)//'_r_eff', r=r_eff(ind))
       CALL channel_halt(substr, status)

       ! define channel elements  
       IF (l_mode_scal) THEN
          CALL new_channel_object(status, modstr//paramnames(ind)//'_gp',     &
               'ff', p2=lnox(ji)%ff)
          CALL channel_halt(substr, status)
          CALL new_attribute(status, modstr//paramnames(ind)//'_gp', 'ff', &
               'long_name', c='unscaled flash frequency')
          CALL channel_halt(substr, status)
          CALL new_attribute(status, modstr//paramnames(ind)//'_gp', 'ff', &
               'units', c='1/s')
          CALL channel_halt(substr, status)

          CALL new_channel_object(status, modstr//paramnames(ind)//'_gp',     &
               'pg', p2=lnox(ji)%pg)
          CALL channel_halt(substr, status)
          CALL new_attribute(status, modstr//paramnames(ind)//'_gp', 'pg', &
               'long_name', c='fraction of CG')
          CALL channel_halt(substr, status)
          CALL new_attribute(status, modstr//paramnames(ind)//'_gp', 'pg', &
               'units', c='[0,1]')
          CALL channel_halt(substr, status)
       END IF

       CALL new_channel_object(status, modstr//paramnames(ind)//'_gp',     &
            'fpscg', p2=lnox(ji)%fpscg)
       CALL channel_halt(substr, status)
       CALL new_attribute(status, modstr//paramnames(ind)//'_gp', 'fpscg', &
            'long_name', c='CG flash frequency')
       CALL channel_halt(substr, status)
       CALL new_attribute(status, modstr//paramnames(ind)//'_gp', 'fpscg', &
            'units', c='1/s')
       CALL channel_halt(substr, status)

       CALL new_channel_object(status, modstr//paramnames(ind)//'_gp',     &
            'fpsic', p2=lnox(ji)%fpsic)
       CALL channel_halt(substr, status)
       CALL new_attribute(status, modstr//paramnames(ind)//'_gp', 'fpsic', &
            'long_name', c='IC flash frequency')
       CALL channel_halt(substr, status)
       CALL new_attribute(status, modstr//paramnames(ind)//'_gp', 'fpsic', &
            'units', c='1/s')
       CALL channel_halt(substr, status)

       CALL new_channel_object(status, modstr//paramnames(ind)//'_gp',       &
            'fpsm2cg', p2=lnox(ji)%fpsm2cg)
       CALL channel_halt(substr, status)
       CALL new_attribute(status, modstr//paramnames(ind)//'_gp', 'fpsm2cg', &
            'long_name', c='CG flash density')
       CALL channel_halt(substr, status)
       CALL new_attribute(status, modstr//paramnames(ind)//'_gp', 'fpsm2cg', &
            'units', c='1/s/m2')
       CALL channel_halt(substr, status)

       CALL new_channel_object(status, modstr//paramnames(ind)//'_gp',       &
            'fpsm2ic', p2=lnox(ji)%fpsm2ic)
       CALL channel_halt(substr, status)
       CALL new_attribute(status, modstr//paramnames(ind)//'_gp', 'fpsm2ic', &
            'long_name', c='IC flash density')
       CALL channel_halt(substr, status)
       CALL new_attribute(status, modstr//paramnames(ind)//'_gp', 'fpsm2ic', &
            'units', c='1/s/m2')
       CALL channel_halt(substr, status)

       CALL new_channel_object(status, modstr//paramnames(ind)//'_gp',       &
            'npcanz', p2=lnox(ji)%npcanz)
       CALL channel_halt(substr, status)
       CALL new_attribute(status, modstr//paramnames(ind)//'_gp', 'npcanz',  &
            'long_name', c='no. of lightning events')
       CALL channel_halt(substr, status)
       CALL new_attribute(status, modstr//paramnames(ind)//'_gp', 'npcanz',  &
            'units', c=' ')
       CALL channel_halt(substr, status)

       CALL new_channel_object(status, modstr//paramnames(ind)//'_gp',       &
            'NOxcg', p2=lnox(ji)%NOxcg)
       CALL channel_halt(substr, status)
       CALL new_attribute(status, modstr//paramnames(ind)//'_gp', 'NOxcg',   &
            'long_name', c='CG NOx lightning emission')
       CALL channel_halt(substr, status)
       CALL new_attribute(status, modstr//paramnames(ind)//'_gp', 'NOxcg',   &
            'units', c='kg(N)')
       CALL channel_halt(substr, status)

       CALL new_channel_object(status, modstr//paramnames(ind)//'_gp',       &
            'NOxic', p2=lnox(ji)%NOxic)
       CALL channel_halt(substr, status)
       CALL new_attribute(status, modstr//paramnames(ind)//'_gp', 'NOxic',   &
            'long_name', c='IC NOx lightning emission')
       CALL channel_halt(substr, status)
       CALL new_attribute(status, modstr//paramnames(ind)//'_gp', 'NOxic',   &
            'units', c='kg(N)')
       CALL channel_halt(substr, status)

       CALL new_channel_object(status, modstr//paramnames(ind)//'_gp',       &
            'cth', p2=lnox(ji)%cth)
       CALL channel_halt(substr, status)
       CALL new_attribute(status, modstr//paramnames(ind)//'_gp', 'cth',     &
            'long_name', c='cloud top height')
       CALL channel_halt(substr, status)
       CALL new_attribute(status, modstr//paramnames(ind)//'_gp', 'cth',     &
            'units', c='m')
       CALL channel_halt(substr, status)

       CALL new_channel_object(status, modstr//paramnames(ind)//'_gp',       &
            'cbh', p2=lnox(ji)%cbh)
       CALL channel_halt(substr, status)
       CALL new_attribute(status, modstr//paramnames(ind)//'_gp', 'cbh',     &
            'long_name', c='cloud bottom height')
       CALL channel_halt(substr, status)
       CALL new_attribute(status, modstr//paramnames(ind)//'_gp', 'cbh',     &
            'units', c='m')
       CALL channel_halt(substr, status)

       CALL new_channel_object(status, modstr//paramnames(ind)//'_gp',       &
            'czh', p2=lnox(ji)%czh)
       CALL channel_halt(substr, status)
       CALL new_attribute(status, modstr//paramnames(ind)//'_gp', 'czh',     &
            'long_name', c='depth of cloud below 0 degC')
       CALL channel_halt(substr, status)
       CALL new_attribute(status, modstr//paramnames(ind)//'_gp', 'czh',     &
            'units', c='m')
       CALL channel_halt(substr, status)

       ! lrestreq is needed for MECO(n), in case case xnox is scaled down 
       ! using mmdclnt; otherwise the results are restart dependend 
       CALL new_channel_object(status, modstr//paramnames(ind)//'_gp',       &
            'xnox', p3=lnox(ji)%xnox_gp, reprid = GP_3D_MID, lrestreq=.TRUE.)
       CALL channel_halt(substr, status)
       CALL new_attribute(status, modstr//paramnames(ind)//'_gp', 'xnox',    &
            'long_name', c='lightning NOx emission')
       CALL channel_halt(substr, status)
       CALL new_attribute(status, modstr//paramnames(ind)//'_gp', 'xnox',    &
            'units', c='kg(N)/s/m3')
       CALL channel_halt(substr, status)

       CALL new_channel_object(status, modstr//paramnames(ind)//'_gp',       &
            'telnox', p3=lnox(ji)%telnox_gp, reprid = GP_3D_MID &
            ,lrestreq=.TRUE.)
       CALL channel_halt(substr, status)
       CALL new_attribute(status, modstr//paramnames(ind)//'_gp', 'telnox',  &
            'long_name', c='lightning NOx emission tendency')
       CALL channel_halt(substr, status)
       CALL new_attribute(status, modstr//paramnames(ind)//'_gp', 'telnox',  &
            'units', c='mol/mol/s')
       CALL channel_halt(substr, status)

       ! ADD OBJECTY WHICH ARE SPECIFIC FOR THE ACTUAL PARAMETERISATION
       SELECT CASE(ind)
       CASE(IPARAM_PaR_T)
          !
          ! NO ADDITIONAL PARAMETER (cth is output for all parameterisations)
          !
       CASE(IPARAM_Grewe)
          !
          CALL new_channel_object(status, modstr//paramnames(ind)//'_gp', &
               'muv', p2=muv, reprid = GP_2D_HORIZONTAL )
          CALL channel_halt(substr, status)
          CALL new_attribute(status, modstr//paramnames(ind)//'_gp', 'muv', &
               'long_name', c='mean updraft velocity')
          CALL channel_halt(substr, status)
          CALL new_attribute(status, modstr//paramnames(ind)//'_gp', 'muv',  &
               'units', c='m/s')
          CALL channel_halt(substr, status)
          !
       CASE(IPARAM_AaP_M)
          !
          CALL new_channel_object(status, modstr//paramnames(ind)//'_gp', &
               'umf', p2=umf, reprid = GP_2D_HORIZONTAL )
          CALL channel_halt(substr, status)
          CALL new_attribute(status, modstr//paramnames(ind)//'_gp', 'umf', &
               'long_name', c='updraft mass flux at sigma = p/ps = 0.44')
          CALL channel_halt(substr, status)
          CALL new_attribute(status, modstr//paramnames(ind)//'_gp', 'umf',  &
               'units', c='kg/m^2/s')
          CALL channel_halt(substr, status)
          !
       CASE(IPARAM_AaP_P)
          !
          CALL new_channel_object(status, modstr//paramnames(ind)//'_gp', &
               'precon', p2=precon, reprid = GP_2D_HORIZONTAL )
          CALL channel_halt(substr, status)
          CALL new_attribute(status, modstr//paramnames(ind)//'_gp', 'precon', &
               'long_name', c='convective precipitation at ground')
          CALL channel_halt(substr, status)
          CALL new_attribute(status, modstr//paramnames(ind)//'_gp', 'precon', &
               'units', c='mm/s')
          CALL channel_halt(substr, status)
          !
       CASE(IPARAM_FinIF)
          !
          CALL new_channel_object(status, modstr//paramnames(ind)//'_gp', &
               'phiice', p2=phiice, reprid = GP_2D_HORIZONTAL )
          CALL channel_halt(substr, status)
          CALL new_attribute(status, modstr//paramnames(ind)//'_gp', 'phiice', &
               'long_name', c='convective cloud ice flux at 440 hPa')
          CALL channel_halt(substr, status)
          CALL new_attribute(status, modstr//paramnames(ind)//'_gp', 'phiice', &
               'units', c='kg/m^2/s')
          CALL channel_halt(substr, status)
          !
       CASE(IPARAM_extIF)
          !
          CALL new_channel_object(status, modstr//paramnames(ind)//'_gp', &
               'phiice', p2=phiice_ext, reprid = GP_2D_HORIZONTAL )
          CALL channel_halt(substr, status)
          CALL new_attribute(status, modstr//paramnames(ind)//'_gp', 'phiice', &
               'long_name', c='convective cloud ice flux through '//&
               &TRIM(c_iif(1))//'::'//TRIM(c_iif(2)))
          CALL channel_halt(substr, status)
          CALL new_attribute(status, modstr//paramnames(ind)//'_gp', 'phiice', &
               'units', c='kg/m^2/s')
          CALL channel_halt(substr, status)
       END SELECT

    END DO parameterisation_loop

    CALL end_message_bi(modstr,'CHANNEL DEFINITION GRID POINT',substr)

!!#D attila +
#if defined(ECHAM5)
    IF (l_calc_lg) THEN

       CALL start_message_bi(modstr, &
            'CHANNEL DEFINITION FOR LAGRANGIAN STUDIES',substr)

       ! define new channel
       CALL new_channel(status, modstr//'_lg', reprid=LG_ATTILA)
       CALL channel_halt(substr, status)

       ! output: only final result 
       CALL new_channel_object(status, modstr//'_lg', 'xnox', p1=xnox_lg &
            , reprid = LG_ATTILA )
       CALL channel_halt(substr, status)
       CALL new_attribute(status, modstr//'_lg', 'xnox' &
            , 'long_name', c='lightning NOx emission')
       CALL channel_halt(substr, status)
       CALL new_attribute(status, modstr//'_lg', 'xnox', 'units', c='kgN/s/m3')
       CALL channel_halt(substr, status)

       CALL new_channel_object(status, modstr//'_lg', 'telnox', p1=telnox_lg &
            , reprid = LG_ATTILA )
       CALL channel_halt(substr, status)
       CALL new_attribute(status, modstr//'_lg', 'telnox' &
            , 'long_name', c='lightning NOx emission tendency')
       CALL channel_halt(substr, status)
       CALL new_attribute(status, modstr//'_lg', 'telnox', 'units' &
            , c='mol/mol/s')
       CALL channel_halt(substr, status)

       CALL end_message_bi(modstr, &
            'CHANNEL DEFINITION FOR LAGRANGIAN STUDIES',substr)

    END IF
#endif
!!#D attila -

    CALL start_message_bi(modstr, 'LOOKING FOR TRACERS', substr)

    CALL get_tracer_list(status, GPTRSTR, 'NO', idt_list_gp, subnames)
    CALL tracer_halt(substr, status)
    nlntrac_gp = SIZE(idt_list_gp)
    IF (p_parallel_io) WRITE(*,*) 'GP-TRACERS:'
    DO jt=1, nlntrac_gp
       IF (p_parallel_io) THEN
          IF (TRIM(subnames(jt)) == '') THEN
             WRITE(*,*) ' ... NO'
          ELSE
             WRITE(*,*) ' ... NO_'//TRIM(subnames(jt))
          END IF
       END IF
#ifdef MESSYTENDENCY
!       CALL mtend_register(my_handle, mtend_id_tracer, idt=idt_list_gp(jt))
       CALL mtend_register(my_handle, idt_list_gp(jt))
#endif
    END DO

!!#D attila +
#if defined(ECHAM5)
    IF (l_calc_lg) THEN
       CALL get_tracer_list(status, LGTRSTR, 'NO', idt_list_lg, subnames)
       CALL tracer_halt(substr, status)
       nlntrac_lg = SIZE(idt_list_lg)
       IF (p_parallel_io) WRITE(*,*) 'LG-TRACERS:'
       DO jt=1, nlntrac_lg
          IF (p_parallel_io) THEN
             IF (TRIM(subnames(jt)) == '') THEN
                WRITE(*,*) ' ... NO'
             ELSE
                WRITE(*,*) ' ... NO_'//TRIM(subnames(jt))
             END IF
          END IF
!#ifdef MESSYTENDENCY ! not yet implemented for LG tracers
!       CALL mtend_register(my_handle, mtend_id_tracer, idt=idt_list_lg(jt))
!#endif          
       END DO
    END IF
#endif
!!#D attila -

    IF (ASSOCIATED(subnames)) DEALLOCATE(subnames)
    NULLIFY(subnames)

    CALL end_message_bi(modstr, 'LOOKING FOR TRACERS', substr)

  END SUBROUTINE lnox_init_memory
  ! ========================================================================

  ! ========================================================================
  SUBROUTINE lnox_init_coupling

    ! LNOX MODULE ROUTINE (ECHAM-5 INTERFACE, PRIVATE)
    !
    ! initialize 'coupling' to online channels
    !
    ! Author: Patrick Joeckel, MPICH, Jan 2004

    USE messy_main_mpi_bi,           ONLY: p_parallel_io
!!$    USE messy_main_tracer_mem_bi, ONLY: GPTRSTR, LGTRSTR
!!$    USE messy_main_tracer_bi,     ONLY: tracer_halt
    USE messy_main_channel_error_bi, ONLY: channel_halt
#ifdef DAHL2000
    USE messy_main_blather_bi,       ONLY: error_bi
    USE messy_main_data_bi,          ONLY: qg_3d
#endif
    USE messy_main_channel,          ONLY: get_channel_object
!!$    USE messy_main_tracer,        ONLY: get_tracer_list
!!$    USE messy_main_constants_mem, ONLY: STRLEN_MEDIUM

    IMPLICIT NONE
    INTRINSIC :: TRIM, ASSOCIATED, SIZE

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER       :: substr = 'lnox_init_coupling'
    INTEGER                           :: status

    CALL start_message_bi(modstr,'COUPLING INITIALIZATION',substr)

!!$    IF (p_parallel_io) THEN
!!$       WRITE(*,*) 'Checking for updraft velocity ...'
!!$       WRITE(*,*) '    channel: ',TRIM(c_updr(1))
!!$       WRITE(*,*) '    object : ',TRIM(c_updr(2))
!!$    END IF
!!$    CALL get_channel_object(status, TRIM(c_updr(1)), TRIM(c_updr(2)) &
!!$         , p3=cu_xupdr)
!!$    CALL channel_halt(substr, status)

    IF (p_parallel_io) THEN
       WRITE(*,*) 'Checking for cloud bottom level ...'
       WRITE(*,*) '    channel: ',TRIM(c_bot(1))
       WRITE(*,*) '    object : ',TRIM(c_bot(2))
    END IF
    CALL get_channel_object(status, TRIM(c_bot(1)), TRIM(c_bot(2))  &
         , p2=cu_jkb)
    CALL channel_halt(substr, status)

    IF (p_parallel_io) THEN
       WRITE(*,*) 'Checking for cloud top level ...'
       WRITE(*,*) '    channel: ',TRIM(c_top(1))
       WRITE(*,*) '    object : ',TRIM(c_top(2))
    END IF
    CALL get_channel_object(status, TRIM(c_top(1)), TRIM(c_top(2))  &
         , p2=cu_jkt)
    CALL channel_halt(substr, status)

!!$    IF (p_parallel_io) THEN
!!$       WRITE(*,*) 'Checking for cloud freezing level ...'
!!$       WRITE(*,*) '    channel: ',TRIM(c_freeze(1))
!!$       WRITE(*,*) '    object : ',TRIM(c_freeze(2))
!!$    END IF
!!$    CALL get_channel_object(status, TRIM(c_freeze(1)), TRIM(c_freeze(2))  &
!!$         , p2=cu_jkfreeze)
!!$    CALL channel_halt(substr, status)

!!$    IF (l_midlevel) THEN
!!$       IF (p_parallel_io) THEN
!!$          WRITE(*,*) 'Checking for cloud bottom level (midlevel) ...'
!!$          WRITE(*,*) '    channel: ',TRIM(c_bot_mid(1))
!!$          WRITE(*,*) '    object : ',TRIM(c_bot_mid(2))
!!$       END IF
!!$       CALL get_channel_object(status, TRIM(c_bot_mid(1)) &
!!$            , TRIM(c_bot_mid(2))  &
!!$            , p2=cu_jkb_mid)
!!$       CALL channel_halt(substr, status)
!!$       
!!$       IF (p_parallel_io) THEN
!!$          WRITE(*,*) 'Checking for cloud top level (midlevel) ...'
!!$          WRITE(*,*) '    channel: ',TRIM(c_top_mid(1))
!!$          WRITE(*,*) '    object : ',TRIM(c_top_mid(2))
!!$       END IF
!!$       CALL get_channel_object(status, TRIM(c_top_mid(1)) &
!!$            , TRIM(c_top_mid(2))  &
!!$            , p2=cu_jkt_mid)
!!$       CALL channel_halt(substr, status)
!!$
!!$       IF (p_parallel_io) THEN
!!$          WRITE(*,*) 'Checking for cloud freezing level (midlevel) ...'
!!$          WRITE(*,*) '    channel: ',TRIM(c_freeze_mid(1))
!!$          WRITE(*,*) '    object : ',TRIM(c_freeze_mid(2))
!!$       END IF
!!$       CALL get_channel_object(status, TRIM(c_freeze_mid(1)) &
!!$            , TRIM(c_freeze_mid(2))  &
!!$            , p2=cu_jkfreeze_mid)
!!$       CALL channel_halt(substr, status)
!!$    END IF

    IF (p_parallel_io) THEN
       WRITE(*,*) 'Checking for convective massflux ...'
       WRITE(*,*) '    channel: ', TRIM(c_massfu(1))
       WRITE(*,*) '    object : ', TRIM(c_massfu(2))
    END IF
    CALL get_channel_object(status, TRIM(c_massfu(1)), TRIM(c_massfu(2)) &
         , p3=umassf)
    CALL channel_halt(substr, status)

    IF (p_parallel_io) THEN
       WRITE(*,*) 'Checking for convective precipitation flux ...'
       WRITE(*,*) '    channel: ', TRIM(c_precflx(1))
       WRITE(*,*) '    object : ', TRIM(c_precflx(2))
    END IF
    CALL get_channel_object(status, TRIM(c_precflx(1)), TRIM(c_precflx(2)) &
         , p3=precflx_cv)
    CALL channel_halt(substr, status)

    IF ((i_ffcalc == IPARAM_extIF) .OR. (i_ffcalc == IPARAM_ALL)) THEN
       IF (p_parallel_io) THEN
          WRITE(*,*) 'Checking for external iso-surface for extended'//&
               &' ice-flux param. ...'
          WRITE(*,*) '    channel: ',TRIM(c_iif(1))
          WRITE(*,*) '    object : ',TRIM(c_iif(2))
       END IF

       CALL get_channel_object(status, TRIM(c_iif(1)), TRIM(c_iif(2))//'_i' &
            , p2=fice_i)
       CALL channel_halt(substr, status)
!!$    CALL get_channel_object(status, TRIM(c_iif(1)), TRIM(c_iif(2))//'_f' &
!!$         , p2=fice_f)
!!$    CALL channel_halt(substr, status)
    END IF
       
#ifdef DAHL2000
    IF ( (i_ffcalc == IPARAM_DahlC) .OR. (i_ffcalc == IPARAM_ALL) &
         .AND. .NOT. ASSOCIATED(qg_3d) ) THEN
       CALL error_bi('Dahl-2010 scheme only applicable with graupel scheme' &
            , substr)
    ENDIF
#endif

    CALL end_message_bi(modstr,'COUPLING INITIALIZATION',substr)

  END SUBROUTINE lnox_init_coupling
  ! ========================================================================

  ! ========================================================================
  SUBROUTINE lnox_physc

    ! LNOx MODULE ROUTINE (ECHAM-5 INTERFACE)
    !
    ! - transfer required fields from/to ECHAM5 for the
    !   calculation of NOx emission by lightning
    ! - add lightning NOx to tendency of NO - tracer(s)
    !
    ! Author: Patrick Joeckel, DLR, Sep 2014

    USE messy_main_constants_mem, ONLY: N_A, M_air, R_gas, g, MN
    USE messy_main_data_bi,       ONLY: &
           slf                & ! land - sea fraction [0 (sea) ... 1 (land)]
         , press_3d           & ! atmospheric pressure [Pa]
         , aps                & ! surface pressure [Pa]
         , tm1_3d             & ! temperature at t-1 [K]
         , tte_3d             & ! temperature tendency [K/s]
         , qm1_3d             & ! specific humidity at t-1 [kg/kg]
         , qte_3d             & ! specific humidity tendency [kg/kg/s]
         , xim1_3d            & ! ice water content at t-1 [kg/kg]
         , xite_3d            & ! ice water content tendency [kg/kg/s]
         , aclc               !& ! (large scale) cloud cover             
    USE messy_main_grid_def_mem_bi, ONLY:jrow, kproma, nlev
    USE messy_main_grid_def_bi,     ONLY: &
           grmass             & ! grid mass [kg]
         , grvol              & ! grid vol. [m3]
         , coslat_2d          & ! cos(latitude)
         , gboxarea_2d        & ! grid box area [m^2]
         , altitude_gnd       & ! height above ground [m]
         , deltaz             !& ! vertical layer thickness [m]
    USE messy_main_timer,         ONLY: time_step_len, delta_time 
#ifndef MESSYTENDENCY
    USE messy_main_tracer_mem_bi, ONLY: pxtte=>qxtte
#endif
    USE messy_main_blather_bi,    ONLY: error_bi

    IMPLICIT NONE
    INTRINSIC :: INT, LOG, MAX, SIZE, REAL

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'lnox_physc'
    INTEGER                     :: jk
    !
    ! POINTERS TO DIAGNOSTIC CHANNEL OBJECTS
    ! flash frequency (unscaled)
    REAL(DP), DIMENSION(:),   POINTER  :: ff => NULL()
    ! CG fraction
    REAL(DP), DIMENSION(:),   POINTER  :: pg => NULL()
    ! flashes per second CG
    REAL(DP), DIMENSION(:),   POINTER  :: fpscg => NULL()
    ! flashes per second IC
    REAL(DP), DIMENSION(:),   POINTER  :: fpsic => NULL()
    ! flashes / (s m^2) CG
    REAL(DP), DIMENSION(:),   POINTER  :: fpsm2cg => NULL()
    ! flashes / (s m^2) IC
    REAL(DP), DIMENSION(:),   POINTER  :: fpsm2ic => NULL()
    ! number of conv. events
    REAL(DP), DIMENSION(:),   POINTER  :: npcanz => NULL()
    ! CG NOx lightn. em. [kg(N)]
    REAL(DP), DIMENSION(:),   POINTER  :: NOxcg => NULL()
    ! IC NOx lightn. em. [kg(N)]
    REAL(DP), DIMENSION(:),   POINTER  :: NOxic => NULL()
    ! cloud top height [m]
    REAL(DP), DIMENSION(:),   POINTER  :: cth => NULL()
    ! cloud bottom height [m]
    REAL(DP), DIMENSION(:),   POINTER  :: cbh => NULL()
    ! cloud depth < 0 degC [m]
    REAL(DP), DIMENSION(:),   POINTER  :: czh => NULL()
    ! lightning NOx emiss [kgN/s/M3] GP
    REAL(DP), DIMENSION(:,:), POINTER  :: xnox_gp   => NULL()
    ! l.NOx tend. [mol/mol/s] GP
    REAL(DP), DIMENSION(:,:), POINTER  :: telnox_gp => NULL()
    !
    ! ADDITIONAL LOCAL FIELDS
    REAL(DP), DIMENSION(kproma,nlev) :: temp   ! temperature [K]
    REAL(DP), DIMENSION(kproma,nlev) :: q      ! specific humidity [kg/kg]
    REAL(DP), DIMENSION(kproma,nlev) :: tvirt  ! virtual temperature
    REAL(DP), DIMENSION(kproma,nlev) :: prhoa  ! air density [kg m-3]
    INTEGER,  DIMENSION(kproma)      :: npcbot ! cloud base [level]
    INTEGER,  DIMENSION(kproma)      :: npctop ! cloud top  [level]
    !
    REAL(DP), DIMENSION(kproma)      :: zcth    ! cloud top height [m]
    REAL(DP), DIMENSION(kproma)      :: zcbh    ! cloud bottom height [m]
    ! depth of cloud above 0 degC level (i.e., T < 0 degC) [m]
    REAL(DP), DIMENSION(kproma)      :: zczh
    INTEGER,  DIMENSION(kproma)      :: npcdh    ! level (index) of 0 degC
    LOGICAL,  DIMENSION(kproma)      :: lconvect ! deep cloud?/active conv.?
    LOGICAL,  DIMENSION(kproma)      :: lcutoff  ! low flash density?
    !
    REAL(DP), DIMENSION(kproma)      :: zpg  ! fraction of CG flashes
    REAL(DP), DIMENSION(kproma)      :: zzpg ! clipped fraction of CG flashes
    !
    REAL(DP), DIMENSION(kproma)      :: zff  ! flash frequency [1/s]
    !
    INTEGER,  DIMENSION(kproma)      :: kref ! reference level (index)
    !
    ! VERTICAL SHAPE
    ! w is the "weight" in [mol(N)/mol(air) / kg(N)], i.e.
    ! multiplied with kg(N)/dtime it yields the tendency in [mol/mol/s]
    ! _l: land, _s: sea/ice, cg: cloud-to-ground, ic: intra-cloud
    REAL(DP), DIMENSION(kproma,nlev) :: wcg_l, wcg_s, wic_l, wic_s
    !
    INTEGER                          :: ji, ind, jp, jt, idt

    IF (i_ffcalc == IPARAM_DahlC) RETURN

    ! UPDATE TEMPERATURE AND HUMIDITY
    temp(:,:) = tm1_3d(_RI_XYZ__(1:kproma,jrow,:)) + tte_3d(_RI_XYZ__(1:kproma,jrow,:)) &
         * time_step_len
    q(:,:) = qm1_3d(_RI_XYZ__(1:kproma,jrow,:)) + qte_3d(_RI_XYZ__(1:kproma,jrow,:)) &
         * time_step_len
    ! VIRTUAL TEMPERATURE [K]
    tvirt(:,:) = temp(:,:) * ( 1._dp + 0.607717_dp * q(:,:) )
    ! AIR DENSITY [kg/m^3]
    prhoa(:,:) = press_3d(_RI_XYZ__(1:kproma,jrow,:)) * M_air * 1.0E-3_dp / &
         (tvirt(:,:) * R_gas)
    
    ! cloud bottom level (index)
    npcbot(:) = INT(cu_jkb(1:kproma,jrow))
    ! cloud top level (index)
    npctop(:) = INT(cu_jkt(1:kproma,jrow))

    ! calculate cloud properties 
    CALL cloud_heights(kproma, nlev             &  ! IN
         , deltaz(_RI_XYZ__(1:kproma,jrow,:))         &  ! IN
         , altitude_gnd(_RI_XYZ__(1:kproma,jrow,:))   &  ! IN
         , temp, npcbot, npctop                 &  ! IN
         , zcth, zcbh, zczh, npcdh, lconvect)      ! OUT

    ! calculate fraction of cloud-to-ground flashes
    CALL cg_fraction(zpg, zczh)

    ! CALCULATE WEIGHTS FOR VERTICAL DISTRIBUTION OF NOx
    SELECT CASE(i_shape)
    CASE(1)
       ! FLAT
       ! LAND
       CALL shape_flat(kproma, nlev, grmass(_RI_XYZ__(1:kproma,jrow,:)) & ! IN
            , npctop(:), npcbot(:), npcdh(:)                  & ! IN
            , .TRUE.                                          & ! IN
            , wcg_l(:,:), wic_l(:,:))
   
       ! SEA/ICE
       CALL shape_flat(kproma, nlev, grmass(_RI_XYZ__(1:kproma,jrow,:)) & ! IN
            , npctop(:), npcbot(:), npcdh(:)                  & ! IN
            , .FALSE.                                         & ! IN
            , wcg_s(:,:), wic_s(:,:))

    CASE(2)
       ! C-shape
       ! LAND
       CALL shape_C(kproma, nlev, grmass(_RI_XYZ__(1:kproma,jrow,:)) & ! IN
            , grvol(_RI_XYZ__(1:kproma,jrow,:))                      & ! IN
            , npctop(:), npcbot(:)                         & ! IN
            , .TRUE.                                       & ! IN
            , wcg_l(:,:))
       wic_l(:,:) = wcg_l(:,:)

       ! SEA/ICE
       CALL shape_C(kproma, nlev, grmass(_RI_XYZ__(1:kproma,jrow,:)) & ! IN
            , grvol(_RI_XYZ__(1:kproma,jrow,:))                      & ! IN
            , npctop(:), npcbot(:)                         & ! IN
            , .FALSE.                                      & ! IN
            , wcg_s(:,:))
       wic_s(:,:) = wcg_s(:,:)            
         
    CASE DEFAULT
       ! NOT REACHED (see lnox_read_nml_ctrl)
    END SELECT

    parameterisation_loop: DO ji = 1, numbers
       IF (numbers == 1) THEN
          ind = i_ffcalc
       ELSE
          ind = ji
       ENDIF

       ! SET POINTER TO ACTUAL CHANNEL OBJECTS
       IF (l_mode_scal) THEN
          ff           => lnox(ji)%ff(1:kproma, jrow)
          pg           => lnox(ji)%pg(1:kproma, jrow)
       ELSE
          NULLIFY(ff)
          NULLIFY(pg)
       END IF
       !
       cth       => lnox(ji)%cth(1:kproma, jrow)
       cbh       => lnox(ji)%cbh(1:kproma, jrow)
       czh       => lnox(ji)%czh(1:kproma, jrow)
       !
       fpscg     => lnox(ji)%fpscg(1:kproma, jrow)
       fpsic     => lnox(ji)%fpsic(1:kproma, jrow)
       fpsm2cg   => lnox(ji)%fpsm2cg(1:kproma, jrow)
       fpsm2ic   => lnox(ji)%fpsm2ic(1:kproma, jrow)
       npcanz    => lnox(ji)%npcanz(1:kproma, jrow)
       !
       NOxcg     => lnox(ji)%NOxcg(1:kproma, jrow)
       NOxic     => lnox(ji)%NOxic(1:kproma, jrow)
       xnox_gp   => lnox(ji)%xnox_gp(_RI_XYZ__(1:kproma,jrow,:))
       telnox_gp => lnox(ji)%telnox_gp(_RI_XYZ__(1:kproma,jrow,:))

       zff(:) = 0.0_dp

       SELECT CASE (ind)

       CASE(IPARAM_PaR_T)

          CALL ff_PaR_cth(zff, zcth, slf(1:kproma,jrow))

       CASE(IPARAM_Grewe)

         CALL mean_updraft_vel(kproma, nlev  & ! IN
               , deltaz(_RI_XYZ__(1:kproma,jrow,:)) & ! IN
               , umassf(_RI_XYZ__(1:kproma,jrow,:))     & ! IN
               , prhoa(:,:)                   & ! IN
               , npcbot(:), npctop(:)         & ! IN
               , muv(1:kproma,jrow) )           ! OUT

          CALL ff_Grewe_muv(zff, zcth, zcbh, muv(1:kproma,jrow) )

       CASE(IPARAM_AaP_M)

          CALL search_level_index(kref(:)   & ! OUT
               , kproma, nlev, 0.44_dp      & ! IN   ! sigma = 0.44
               , press_3d(_RI_XYZ__(1:kproma,jrow,:)) & ! IN
               , aps(1:kproma,jrow) )         ! IN

          DO jp=1, kproma
             jk = kref(jp)
             IF (jk > 0) THEN
                umf(jp,jrow) = umassf(_RI_XYZ__(jp,jrow,jk))
             ELSE
                umf(jp,jrow) = 0.0_dp
             END IF
          END DO

          CALL ffcg_AaP_umf(fpscg(:)          & ! OUT
               , umf(1:kproma,jrow)           & ! IN
               , gboxarea_2d(1:kproma,jrow) )   ! IN  

       CASE(IPARAM_AaP_P)

          precon(1:kproma,jrow) = precflx_cv(_RI_XYZ__(1:kproma,jrow,nlev))      

          CALL ffcg_AaP_precip(fpscg(:)       & ! OUT
               , precon(1:kproma,jrow)        & ! IN
               , slf(1:kproma,jrow)           & ! IN
               , gboxarea_2d(1:kproma,jrow) )   ! IN

       CASE(IPARAM_FinIF)

          CALL search_level_index(kref(:)     & ! OUT
               , kproma, nlev, 44000.0_dp     & ! IN   ! p = 440 hPa
               , press_3d(_RI_XYZ__(1:kproma,jrow,:)) )   ! IN

          DO jp=1, kproma
             jk = kref(jp)
             IF ( (jk > 0) .AND. &
                  (aclc(_RI_XYZ__(jp,jrow,jk)) >= 0.01_dp) ) THEN
                phiice(jp,jrow) = &
                     ( xim1_3d(_RI_XYZ__(jp,jrow,jk)) + &
                     xite_3d(_RI_XYZ__(jp,jrow,jk)) * time_step_len ) &
                     * umassf(_RI_XYZ__(jp,jrow,jk)) &
                     / aclc(_RI_XYZ__(jp,jrow,jk))
             ELSE
                phiice(jp,jrow) = 0.0_dp
             END IF
          END DO

          CALL ff_Finney_cif(zff(:)           & ! OUT
               , phiice(1:kproma,jrow)        & ! IN
               , slf(1:kproma,jrow)           & ! IN
               , gboxarea_2d(1:kproma,jrow)   & ! IN
               , ind)                           ! IN

       CASE(IPARAM_extIF)

          DO jp=1, kproma
             jk = NINT(fice_i(jp,jrow))
             IF ( (jk > 0) .AND. &
                  (aclc(_RI_XYZ__(jp,jrow,jk)) >= 0.01_dp) ) THEN
                phiice_ext(jp,jrow) = &
                     ( xim1_3d(_RI_XYZ__(jp,jrow,jk)) + &
                     xite_3d(_RI_XYZ__(jp,jrow,jk)) * time_step_len ) &
                     * umassf(_RI_XYZ__(jp,jrow,jk)) &
                     / aclc(_RI_XYZ__(jp,jrow,jk))
             ELSE
                phiice_ext(jp,jrow) = 0.0_dp
             END IF
          END DO

          CALL ff_Finney_cif(zff(:)           & ! OUT
               , phiice_ext(1:kproma,jrow)    & ! IN
               , slf(1:kproma,jrow)           & ! IN
               , gboxarea_2d(1:kproma,jrow)   & ! IN
               , ind)                           ! IN

       END SELECT

       ! calculate CG and IC flash frequency <-> total flash frequency
       CALL ffs(ind, zpg(:), zff(:), fpscg(:), fpsic(:))

       ! APPLY CUT-OFF CRITERIA
       ! 1st criterion: lconvect -> cloud at least 3000m thick (see above)
       ! 2nd criterion: lflash   -> cut-off low flash densities,
       !                            if not in scaling-mode
       lcutoff(:) = lflash(zff           &
            , gboxarea_2d(1:kproma,jrow) &
            , coslat_2d(1:kproma,jrow) )
       !
       DO jp=1, kproma
          IF (lconvect(jp) .AND. lcutoff(jp) ) THEN
             npcanz(jp) = 1.0_dp
          ELSE
             npcanz(jp) = 0.0_dp
          ENDIF
       END DO
       !
       zff(:)   = zff(:)   * npcanz(:)
       zzpg(:)  = zpg(:)   * npcanz(:)
       fpscg(:) = fpscg(:) * npcanz(:)
       fpsic(:) = fpsic(:) * npcanz(:)

       ! SAVE CLOUD PROPERTIES, WHERE CONDITIONS ABOVE ARE FULFILLED
       cth(:) = zcth(:) * npcanz(:)
       cbh(:) = zcbh(:) * npcanz(:)
       czh(:) = zczh(:) * npcanz(:)

       ! SAVE UNSCALED (!) FLASH FREQUENCY AND CG FRACTION
       ! IN SCALING MODE; NOTE: lcutoff = .TRUE. for l_mode_scal = .TRUE.
       IF (ASSOCIATED(ff)) ff(:) = zff(:)  ! copy result from above
       IF (ASSOCIATED(pg)) pg(:) = zzpg(:) ! copy result from above

       ! flash densities [1/m^2/s]
       fpsm2cg(:) = fpscg(:)/gboxarea_2d(1:kproma,jrow)
       fpsm2ic(:) = fpsic(:)/gboxarea_2d(1:kproma,jrow)

       SELECT CASE(ind)
          !
       CASE(IPARAM_PaR_T)
          ! NOTING TO DO
       CASE(IPARAM_Grewe)
          muv(1:kproma,jrow) = muv(1:kproma,jrow) * npcanz(:)
       CASE(IPARAM_AaP_M)
          umf(1:kproma,jrow) = umf(1:kproma,jrow) * npcanz(:)
       CASE(IPARAM_AaP_P)
          precon(1:kproma,jrow) = precon(1:kproma,jrow) * npcanz(:)
       CASE(IPARAM_FinIF)
          phiice(1:kproma,jrow) = phiice(1:kproma,jrow) * npcanz(:)
       CASE(IPARAM_extIF)
          phiice_ext(1:kproma,jrow) = phiice_ext(1:kproma,jrow) * npcanz(:)
          !
       END SELECT

       ! CALCULATE total, vertically integrated NOx production
       ! [kg(N)] for current time step in actual column
       CALL nox_prod(ind, delta_time & ! IN 
            , fpscg(:), fpsic(:)     & ! IN
            , NOxcg(:), NOxic(:) )     ! OUT

       ! CALCULATE LIGHTNING NOX TENDENCY [mol/mol/s]
       DO jp=1,kproma
          telnox_gp(jp,:) = ( &
               NOxcg(jp) * wcg_l(jp,:) * slf(jp,jrow)            +   &
               NOxcg(jp) * wcg_s(jp,:) * (1.0_dp - slf(jp,jrow)) +   &
               NOxic(jp) * wic_l(jp,:) * slf(jp,jrow)            +   &
               NOxic(jp) * wic_s(jp,:) * (1.0_dp - slf(jp,jrow))   ) &
               / delta_time
               
       END DO

       ! CALCULATE LIGHTNING NOX EMISSION [kg(N)/m^3/s]
       ! mol(N)/mol(air)/s:  * (k)g(N)/mol(N)      -> (k)g(N)/mol(air)/s
       !                     / (k)g(air)/mol(air)  -> (k)g(N)/(k)g(air)/s
       !                     * kg(air)/m^3         -> kg(N)/m^3/s
       xnox_gp(:,:) = telnox_gp(:,:) * ( MN / M_air) * prhoa(:,:)

    END DO parameterisation_loop

    ind = 1
    if (i_ffcalc == IPARAM_ALL) ind = i_ff_cpl

    IF (nlntrac_gp > 0) THEN
       ! LOOP OVER ALL NO-TRACERS AND ADD LIGHTNING NOx TO TENDENCY
       DO jt=1, nlntrac_gp
          idt = idt_list_gp(jt)
#ifndef MESSYTENDENCY
          pxtte(_RI_X_ZN_(1:kproma,:,idt)) = &
               pxtte(_RI_X_ZN_(1:kproma,:,idt)) + &
               lnox(ind)%telnox_gp(_RI_XYZ__(1:kproma,jrow,:))
#else
          CALL mtend_add_l(my_handle, idt, &
               px=lnox(ind)%telnox_gp(_RI_XYZ__(1:kproma,jrow,:)))
#endif
       END DO
    END IF

  END SUBROUTINE lnox_physc
  ! ========================================================================

  ! ========================================================================
  SUBROUTINE lnox_global_end
!!#D attila +
#ifdef ECHAM5
    USE messy_main_tracer_mem_bi, ONLY: pxtte_a=>qxtte_a, NCELL 
    USE messy_attila_tools_e5,    ONLY: gp2lg_e5

    IMPLICIT NONE

    INTEGER :: jt, jn, ind

    ! Lagrangian calculations
    IF (l_calc_lg) THEN
       ind = 1
       if (i_ffcalc == IPARAM_ALL) ind = i_ff_cpl
       CALL gp2lg_e5(lnox(ind)%telnox_gp, telnox_lg, lmcons=.FALSE.)
       CALL gp2lg_e5(lnox(ind)%xnox_gp, xnox_lg, lmcons=.FALSE.)
       ! LOOP OVER ALL NO-TRACERS AND ADD LIGHTNING NOx TO TENDENCY
       IF (nlntrac_lg > 0) THEN
          DO jt=1, nlntrac_lg
!#ifndef MESSYTENDENCY ! not yet implemented for LG tracers
             DO jn=1,NCELL
                pxtte_a(jn,idt_list_lg(jt)) = pxtte_a(jn,idt_list_lg(jt)) + &
                     telnox_lg(jn)
             END DO
!#else
!             CALL mtend_add_?(my_handle, mtend_id_tracer &
!               px=lnox(ind)%telnox_lg(:), idt=idt_list_lg(jt))
!#endif
          END DO
       END IF
    END IF
#endif
!!#D attila -

#ifdef DAHL2000
    CALL lnox_dahl2000_global_end

#endif

  END SUBROUTINE lnox_global_end
  ! ========================================================================

  ! ========================================================================
  SUBROUTINE lnox_free_memory

    ! LNOx MODULE ROUTINE (ECHAM-5 INTERFACE)
    !
    ! deallocate memory
    !
    ! Author: Patrick Joeckel, MPICH, Aug 2003

    IMPLICIT NONE
    INTRINSIC :: ASSOCIATED

    IF (ASSOCIATED(idt_list_gp)) DEALLOCATE(idt_list_gp)
    IF (ASSOCIATED(idt_list_lg)) DEALLOCATE(idt_list_lg)
    IF (ASSOCIATED(lnox))        DEALLOCATE(lnox)

  END SUBROUTINE lnox_free_memory
  ! ========================================================================

  ! ========================================================================
  SUBROUTINE lnox_read_nml_cpl(status, iou)

    ! LNOx MODULE ROUTINE (ECHAM-5 INTERFACE, PRIVATE)
    !
    ! read namelist for 'coupling' to ECHAM5
    !
    ! Author: Patrick Joeckel, MPICH, Aug 2003

    USE messy_main_tools, ONLY: read_nml_open, read_nml_check, read_nml_close
    USE messy_main_tracer_mem_bi, ONLY: NGCELL

    IMPLICIT NONE
    INTRINSIC :: TRIM

    ! I/O
    INTEGER, INTENT(OUT) :: status     ! error status
    INTEGER, INTENT(IN)  :: iou        ! I/O unit
    ! Switch for skipping calculation of Lagrangian rate coefficients.
    ! It is local,not broadcasted!
    LOGICAL                      :: l_skip_lg = .FALSE.

!!$    NAMELIST /CPL/ c_updr, c_top, c_bot, c_freeze &
!!$         , l_midlevel, c_top_mid, c_bot_mid, c_freeze_mid, l_skip_lg, i_ff_cpl &
!!$         , c_massfu, c_precflx, c_iif
    NAMELIST /CPL/ c_top, c_bot &
         , l_skip_lg, i_ff_cpl &
         , c_massfu, c_precflx &
         , c_iif

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr='lnox_read_nml_cpl'
    LOGICAL              :: lex      ! file exists ?
    INTEGER              :: fstat    ! file status

    status = 1

    CALL read_nml_open(lex, substr, iou, 'CPL', modstr)
    IF (.NOT.lex) RETURN    ! <modstr>.nml does not exist

    READ(iou, NML=CPL, IOSTAT=fstat)
    CALL read_nml_check(fstat, substr, iou, 'CPL', modstr)
    IF (fstat /= 0) RETURN  ! error while reading namelist

    ! DIAGNOSE NAMELIST AND SET GLOBAL SWITCHES
!!$    WRITE(*,*) 'ONLINE UPDRAFT VELOCITY:'
!!$    WRITE(*,*) '    channel: ',TRIM(c_updr(1))
!!$    WRITE(*,*) '    object : ',TRIM(c_updr(2))

    WRITE(*,*) 'ONLINE CLOUD BOTTOM LEVEL:'
    WRITE(*,*) '    channel: ',TRIM(c_bot(1))
    WRITE(*,*) '    object : ',TRIM(c_bot(2))

    WRITE(*,*) 'ONLINE CLOUD TOP LEVEL:'
    WRITE(*,*) '    channel: ',TRIM(c_top(1))
    WRITE(*,*) '    object : ',TRIM(c_top(2))

!!$    WRITE(*,*) 'ONLINE CLOUD FREEZING LEVEL:'
!!$    WRITE(*,*) '    channel: ',TRIM(c_freeze(1))
!!$    WRITE(*,*) '    object : ',TRIM(c_freeze(2))

    IF (i_ffcalc /= IPARAM_ALL) THEN
       IF (i_ff_cpl /= i_ffcalc) WRITE(*,*) &
            'Coupling adjusted according to available LNOX parameterisation!'
       i_ff_cpl = i_ffcalc
    ENDIF

    WRITE(*,*) 'Coupling flash frequency of parameterisation...'
    SELECT CASE (i_ff_cpl)
    CASE(IPARAM_PaR_T)
       WRITE(*,*) '... Price and Rind, 1992'
    CASE(IPARAM_Grewe)
       WRITE(*,*) '... Grewe et al., 2001'
    CASE(IPARAM_AaP_M)
       WRITE(*,*) '... Allen and Pickering (Massflux), 2002'
    CASE(IPARAM_AaP_P)
       WRITE(*,*) '... Allen and Pickering (Precipitation), 2002'
    CASE(IPARAM_FinIF)
       WRITE(*,*) '... Finney et al. (Cloud ice flux), 2014'
    CASE(IPARAM_extIF)
       WRITE(*,*) '... Finney et al. (Cloud ice flux), 2014, extended'
    CASE(IPARAM_DahlC)
#ifdef DAHL2000
       WRITE(*,*) '... Dahl (2010)'
!#else 
!       WRITE(*,*) 'ERROR: Dahl (2010) not applicable in EMAC'
!       RETURN
#endif
    CASE DEFAULT
       WRITE (*,*) 'ERROR: Parameterisation unknown! ', i_ff_cpl
       RETURN
    END SELECT
    WRITE(*,*) '... to NOx emission.'

#if defined(ECHAM5)
    IF ((.NOT. l_skip_lg) .AND. (NGCELL > 0)) THEN
       l_calc_lg = .TRUE.
!!#D attila +
       WRITE(*,*) 'Lagrangian: ON'
!!#D attila -
    ELSE
       IF (.NOT. l_skip_lg) THEN
!!#D attila +
          WRITE(*,*) 'l_skip_lg = T in namelist'
          WRITE(*,*) 'However no Lagrangian scheme activated ...'
          WRITE(*,*) ' ... setting l_calc_lg = F'
!!#D attila -
       ENDIF
       l_calc_lg = .FALSE.
!!#D attila +
       WRITE(*,*) 'Lagrangian: OFF'
!!#D attila -
    ENDIF
#endif

    CALL read_nml_close(substr, iou, modstr)
    status = 0 ! NO ERROR

  END SUBROUTINE lnox_read_nml_cpl
  ! ========================================================================


#ifdef DAHL2000
  ! ========================================================================
  SUBROUTINE lnox_dahl2000_global_end

!#if defined(COSMOv4s8)
! COSMO > v4s8 :  Adjustment to humidity tracers required     
                                ! Dahl, 2010
                                ! qqq this needs to be cleaned up
    USE messy_main_blather_bi,    ONLY: error_bi
    USE messy_main_mpi_bi,        ONLY: p_parallel_io, p_nprocs, p_pe &
         , ij_local, exchange_boundaries          &
         , global_fields
    USE messy_main_grid_def_mem_bi, ONLY: ie, je, ke, ie_tot, je_tot, ke_tot & 
         , jstartpar, jendpar, istartpar, iendpar, dlon, dlat         &
         , startlat_tot, startlon_tot
    USE messy_main_grid_def_mem_bi, ONLY:hhl 
    USE messy_main_data_bi,       ONLY: &
         , qi_3d, qs_3d, qg_3d, qc_3d, rho        &
         , tm1_3d, ttens
    USE messy_main_timer_bi,      ONLY: event_state
    USE messy_main_timer,         ONLY: current_time_step, time_step_len &
         , current_date
    USE messy_main_timer_event,   ONLY: event_delta
    USE messy_main_constants_mem, ONLY: pi, r_earth => radius_earth
    USE messy_main_tools,         ONLY: find_next_free_unit, int2str

    IMPLICIT NONE

    INTRINSIC :: EXP, MAXVAL, NINT, REAL, SQRT

    ! LOCAL 
    ! Structure containing relevant/final cell info
    ! ---------------------------------------------

    TYPE storm
       INTEGER  ::  centroid_x                          ! [1  ]
       INTEGER  ::  centroid_y                          ! [1  ]
       REAL(DP) ::  flash_rate                          ! [1/s]
    END TYPE storm

    TYPE(storm), DIMENSION(:), POINTER :: cell
    ! Number of vertically overlapping clusters
    INTEGER                            :: number_overlaps 
    ! number of CBs (length of final csize) this is just the number of 
    ! nontrivial elements (either length of CSIZE or number of overlapping 
    ! capacitor plates if itype_lightning == 1; assigned in LABEL 
    ! and re-assigned in DISTRIBUTE_FLASHES routines) 
    INTEGER, SAVE                      :: cs_count
    INTEGER                            :: accum_flashes
    INTEGER                            :: vol_index
    INTEGER                            :: mg_index

    INTEGER, PARAMETER :: csize_length = 1000000 ! Length of global csize array 
    INTEGER, PARAMETER :: n_overlapping_elements = 100000! max allowed number of 
    ! overlapping clusters (x2)
    REAL(dp), PARAMETER ::  &
         qsi_binmar  = 0.1_dp,       & 
                                !grar_binmar = 0.1_dp,       &
         mg_binmar   = 0.1_dp!,       &
    !w_binmar    = 2.0_dp
    REAL(dp) :: sinclight 

    REAL(dp), DIMENSION(ie,je,ke) :: temp

    REAL(dp), DIMENSION(ie,je,ke) :: qci
    REAL(dp), DIMENSION(ie,je,ke) :: qsi
    REAL(dp), DIMENSION(ie,je,ke) :: mg3d

    INTEGER,  DIMENSION(ie,je,ke) :: qsi_bin
    INTEGER,  DIMENSION(ie,je,ke) :: mg_bin
    INTEGER,  DIMENSION(ie,je,ke) :: mg_bin_lab

    REAL(dp), ALLOCATABLE  ::  &
         fl_lon(:),                         &  ! lon/lat/t of discharges
         fl_lat(:),                         &
         fl_time(:)

    REAL(dp)  :: &
         plate_depth,         &  ! Mean of upper and lower plate depths
         plate_dist,          &
         geometric_term_num,  &
         geometric_term_den,  &
         geometric_term,      &
         sigma_crit,          &
         area,                &  ! area of the capacitor plates
         plate_strength,      &  ! GRAR or QG, determining charge density        
         sgr_diam,            &
         sv_gr,               &
         srho_charge,         &
         charge_volume,       &
         sdelta_q,            &
         e_crit,              & ! e_crit 
         alt_be,              & ! Breakeven altitude in m
         q_crit,              & ! amount of charge required for E_crit
         rho_crit,            & ! respective 3D chage density  
         j_gen,               & ! generator current density
         charge_after_flash,  & ! 
         sigma_after_flash,   & ! 2D charge density after the flash
         efac,                & ! factor for E-field expression
         e_after_flash,       & ! E-field after flash
         efne,                & ! E-Field Neutralization Efficiency
         crme,                & ! Charge Removal Efficiency 
         sigma_efficiency,    & ! same but w.r.t. sigma
         f_rate_term1,        & ! factors of flash-rate expression
         f_rate_term2,        & !
         lf_dl10,             & ! flash rate 
         diameter,            & ! graupel-region's equivalent circular diameter      
         mcs_diameter,        & ! Diameter beyond which storm is MCS (def 15 km)
         gauss_width,         & ! Describes decay rate towards 0.4 from 1.0
         mcs_correction,      & ! Correction factor [0.4, 1.0]
         lf_dl10_uc,          & ! non-MCS-corrected flash rate (too large)
         gauss_arg              ! Exponent of Gauss function

    LOGICAL  :: lverbose = .TRUE. ! qqq AK   namelist ?
    LOGICAL  :: levent = .FALSE.
    CHARACTER (LEN=25) :: yerrmsg  ! for MPI error message
    INTEGER  :: status
    INTEGER  :: ii, jj
    CHARACTER(LEN=*), PARAMETER :: substr = 'lnox_global_end'

    IF (i_ffcalc /= IPARAM_DahlC) RETURN

    levent    = event_state(LIGHT_EVENT, current_date)
    sinclight = event_delta(LIGHT_EVENT)

    IF (.NOT. levent) RETURN

    ! preparation of physical variables
    temp(1:ie,1:je,1:ke) = &
         tm1_3d(1:ie,1:je,1:ke) + ttens(1:ie,1:je,1:ke) * time_step_len

    ! Prepare upper region (positive charge region)
    ! ---------------------------------------------
    qci(:,:,:) = 1.0E3_dp * ( qc_3d(:,:,:) + qi_3d(:,:,:) )
    qsi(:,:,:) = 1.0E3_dp * ( qs_3d(:,:,:) + qi_3d(:,:,:) )
    qsi_bin(:,:,:) = BIN_FIELD(qsi, qsi_binmar)

    ! Prepare lower region (negative-charge region)
    ! ----------------------------------------------
    mg3d(:,:,:) =  1.0E3_dp * rho(:,:,:) * qg_3d(:,:,:)  ! in g/kg

    ! Only consider that part of the field where T < -13 C (charge-reversal)
    WHERE (temp(:,:,:) > tcr) mg3d(:,:,:) = 0.0_dp
    mg_bin(:,:,:) = BIN_FIELD(mg3d, mg_binmar)

    !---------------------------------------------------------------------------
    ! Section 2.1: Search storm cells
    !   In the subroutine cluster_overlaps, the capacitor_details routine is 
    !   called which updates global variables that are used farther below (like
    !   number_overlaps).
    !---------------------------------------------------------------------------

    CALL cluster_overlaps (qsi_bin, mg_bin, mg_bin_lab)

    !----------------------------------------------------------------------------
    ! Section 3: Calculate the flash frequency
    !----------------------------------------------------------------------------

    ! Some exception handling

    IF (p_pe == 0 .AND. number_overlaps > max_number_cap) THEN
       WRITE (*,*) '     *** SRC_LIGHTNING.DAHL_2010 ERROR ***  '
       WRITE (*,*) 'SET max_number_cap VARIABLE TO AT LEAST', max_number_cap
    ENDIF

    ALLOCATE (cell(number_overlaps))

    cell % flash_rate = 0.0_dp
    cell % centroid_x = 0
    cell % centroid_y = 0

    IF (number_overlaps > 0) THEN

       DO ii = 1, number_overlaps

          plate_depth = &
               0.5_dp * (info_cap(ii) % top_depth + info_cap(ii) % bot_depth)

          plate_dist = info_cap(ii) % plate_distance 

          area = info_cap(ii) % bot_area

          !----------------------------------------------------------------------
          ! Section 3.1:
          !   Determine variables that are related MG:
          !   graupel diameter (sgr_diam), sedimentation velocity (sv_gr), and
          !   charge density in generator current (srho_charge)
          !----------------------------------------------------------------------

          plate_strength = 1.0E3_dp * info_cap(ii) % graupel

          IF (plate_strength >= MAXVAL(mg_range)) THEN
             mg_index = 100
          ELSE
             mg_index = 1
             mg_loop:  DO jj = 1, 100
                IF (mg_range(jj) >= plate_strength) THEN
                   mg_index = jj
                   EXIT mg_loop
                ENDIF
             ENDDO mg_loop
          ENDIF

          sgr_diam    = gr_diam   (mg_index)
          sv_gr       = v_gr      (mg_index)  
          srho_charge = rho_charge(mg_index)

          !----------------------------------------------------------------------
          ! Section 3.3: Set up variables related to space-charge volume
          !----------------------------------------------------------------------

          charge_volume = plate_depth * area

          IF (charge_volume >= MAXVAL(v_range)) THEN
             vol_index = 100
          ELSE
             vol_index = 1
             vol_loop:  DO jj = 1, 100
                IF (v_range(jj) >= charge_volume) THEN
                   vol_index = jj
                   EXIT vol_loop
                ENDIF
             ENDDO vol_loop
          ENDIF

          sdelta_q = delta_q(vol_index)  

          !----------------------------------------------------------------------
          ! Section 3.4: Calculate lightning frequency
          !----------------------------------------------------------------------

          ! Breakeven field strength 
          ! ------------------------ 

          alt_be = info_cap(ii) % breakdown_alt  
          e_crit = 1.0E3_dp * efac1 * EXP(-alt_be / efac2)   

          ! Geometric term
          ! --------------

          geometric_term_num = plate_dist
          geometric_term_den = &
               SQRT(area / pi + (0.5_dp * plate_dist) * (0.5_dp * plate_dist)) 

          geometric_term = geometric_term_num / geometric_term_den

          ! Critical charge (densities)
          ! -------------------------------------

          sigma_crit = 2.0_dp * eps / (geometric_term - 2.0_dp) * e_crit

          q_crit     = sigma_crit * area 
          rho_crit   = q_crit / (area * plate_depth) 

          ! Generator current
          ! -----------------

          j_gen = srho_charge * sv_gr         

          ! Make sure no more charge is transferred than is present initially

          IF (sdelta_q > q_crit) THEN
             sdelta_q = q_crit
          ENDIF

          charge_after_flash = q_crit - sdelta_q
          sigma_after_flash  = charge_after_flash / area

          ! E-field after the flash
          ! -----------------------  

          efac = 0.5_dp * sigma_after_flash / eps
          e_after_flash = -sigma_after_flash / eps + efac * geometric_term

          ! E-field-neutralization and charge-removal efficiencies (efne, crme)
          ! -------------------------------------------------------------------

          efne = (e_crit - e_after_flash) / e_crit
          crme = (q_crit - charge_after_flash) / q_crit

          sigma_efficiency = (sigma_crit - sigma_after_flash) / sigma_crit

          ! Calculate flash rate
          ! --------------------

          f_rate_term1 = (geometric_term / (2.0_dp * eps) - &
               1.0_dp/eps) / e_crit
          f_rate_term2 = j_gen / efne
          lf_dl10_uc   = gamma * f_rate_term1 * f_rate_term2  ! MCS-uncorrected

          ! Apply MCS-correction (flash rate reduced by up to 60 %)
          ! -------------------------------------------------------

          diameter     = 2.0E-3_dp * SQRT(area / pi)  ! in km
          mcs_diameter = 15.0_dp                      ! in km
          gauss_width  = 0.085_dp

          IF (diameter    <= mcs_diameter) THEN 
             mcs_correction = 1.0_dp
          ELSEIF (diameter > mcs_diameter) THEN 
             gauss_arg      = gauss_width * (diameter - mcs_diameter)
             mcs_correction = 0.4_dp + 0.6_dp * EXP(-gauss_arg * gauss_arg)
          ENDIF

          lf_dl10 = mcs_correction * lf_dl10_uc 

          cell(ii) % flash_rate = lf_dl10

          !----------------------------------------------------------------------
          ! Section 3.5: Print thunderstorm cell information
          !----------------------------------------------------------------------

          IF (lverbose .AND. p_pe == 0) THEN     
             print *, '****************************************************'
             print *, '     SRC_LIGHTNING.DAHL_2010: Verbose output'
             print *, '****************************************************'
             print *, ''     
             print '(a35, i7, a25, i7, f7.2)', 'Max number of cells:',          &
                  number_overlaps, 'at time step/hour:',  current_time_step &
                  , REAL(current_time_step,dp) * time_step_len / 3600.0_dp 
             print '(a35, i5)', 'Details of storm (label):', info_cap(ii) % bot_label 
             print *, ''
             print *, '             GEOMETRIC PARAMETERS:'
             print *, '             ---------------------'
             print '(a35, i6, a1, i4)', 'Position (X,Y):',    &
                  info_cap(ii) % x_pos, ',',  info_cap(ii) % y_pos 
             print '(a35, f6.1, a6)', 'Equivalent circular diameter:',    &
                  2.0E-3_dp * SQRT(area / pi),                ' km'
             print '(a35, f6.1, a6)', 'Centroid distance:',    &
                  1.0E-3_dp * plate_dist, ' km'
             print '(a35, f6.1, a6)', 'Plate SFC distance:',    &
                  1.0E-3_dp * info_cap(ii) % separation,         ' km' 
             print '(a35, f6.1, a6)', 'Plate depths:',    & 
                  1.0E-3_dp * plate_depth,                   ' km' 
             print '(a35, f6.1, a6)', 'Storm-top height:',    &
                  1.0E-3_dp * info_cap(ii) % total_height,       ' km' 
             print '(a35, f6.1, a6)', 'Total storm depth:',    &
                  1.0E-3_dp*(info_cap(ii)%separation+2.0*plate_depth), ' km'
             print '(a35, f6.1, a6)', 'Charge volume:',    &
                  1.0E-9_dp * charge_volume,                   ' km**3' 
             print '(a35, f6.1, a8)', 'Max graupel content (QG)', &
                  plate_strength, 'g/kg'
             print '(a35, f8.4)', 'Geometric term', geometric_term
             print *, '' 
             print *, '              ELECTRIC PARAMETERS:'
             print *, '             ---------------------'
             print '(a35, f6.1, a6)', 'Critical field:',    & 
                  1.0E-3_dp * e_crit,                            ' kV/m'
             print '(a35, f6.1, a6)', 'Total charge per plate:', q_crit,  ' C'
             print '(a35, e12.3, a8)', 'Charge per volume: ', rho_crit,    'C/m**3' 
             print '(a35, e12.3, a8)', 'Charge per area:', sigma_crit,     'C/m**2'
             print '(a35, f6.3, a6)', 'Graupel diameter:', 100.0 * sgr_diam, ' cm'
             print '(a35, f6.2, a6)', 'Sedimentation velocity:', sv_gr,    ' m/s'
             print '(a35, e12.3, a8)', 'Charging current density:', srho_charge, 'C/m**3'
             print '(a35, f6.1, a6)', 'Removed charge:',          &
                  q_crit - charge_after_flash,                       ' C'
             print '(a35, f6.1, a6)', 'Removed E-field:',          &
                  1.0E-3_dp * (e_crit - e_after_flash),          ' kV/m' 
             print '(a35, f6.3)', 'Charge removal efficiency:', crme             
             print '(a35, f6.3)', 'Field neutralization efficiency:', efne
             print '(a35, f8.4, a10)', 'Uncorrected flash rate:', &
                  lf_dl10_uc * 60.0_dp,                          ' 1/min'
             print '(a35, f8.4, a10)', 'MCS-correction:', mcs_correction
             print '(a35, f8.4, a10)', 'Flash rate:', &
                  lf_dl10 * 60.0_dp,                             ' 1/min'
             print *, ''
          ENDIF

          cell(ii) % centroid_x = info_cap(ii) % x_pos
          cell(ii) % centroid_y = info_cap(ii) % y_pos

       ENDDO
    ELSE IF (number_overlaps == 0) THEN
       lf_dl10 = 0.0_dp
    ENDIF

    !----------------------------------------------------------------------------
    ! Section 4: Distribute flashes around the centroid (space and time)
    !----------------------------------------------------------------------------

    accum_flashes = 0

    DO jj = 1, number_overlaps

       accum_flashes = accum_flashes + &
            !            NINT( cell(jj) % flash_rate * 3600.0_dp * hinclight )
            NINT( cell(jj) % flash_rate * sinclight )

    ENDDO

    !WRITE (*,'(f4.1, a25, i6)') hinclight * 60., &
    !            '-min accumulated flashes:', accum_flashes                         

    IF (accum_flashes > 0) THEN
       ALLOCATE(fl_lon (accum_flashes))
       ALLOCATE(fl_lat (accum_flashes))
       ALLOCATE(fl_time(accum_flashes))
    ELSEIF (accum_flashes == 0) THEN
       ALLOCATE(fl_lon (2))
       ALLOCATE(fl_lat (2))
       ALLOCATE(fl_time(2))
    ENDIF

    CALL distribute_flashes (cell, accum_flashes, fl_lon, fl_lat, &
         fl_time)

    !----------------------------------------------------------------------------
    ! Deallocate dynamic arrays
    !----------------------------------------------------------------------------

    DEALLOCATE(cell)
    DEALLOCATE(fl_lon)
    DEALLOCATE(fl_lat)
    DEALLOCATE(fl_time)

  CONTAINS

    !=========================================================================
    ! Subroutine cluster_overlaps
    !=========================================================================

    !-------------------------------------------------------------------------
    ! Description:
    !   This internal subroutine looks for overlapping clusters of two classes.
    !   The idea is to find the positively and negatively charged plates of
    !   the capacitor.  The positive plate is assumed to be the anvil 
    !   (or part of it), and the negative plate is the region where graupel
    !   exists.    
    !   Hence, regions where an ice cloud overlaps with a region of  
    !   graupel are identified as potential thundercloud.  
    !   The information about the top and bottom clusters are stored in the
    !   derived-type variable ovl_clusters (for overlapping clusters).
    !
    ! Method:
    !   Calls of LABEL and CLUSTER_ANALYSIS routines separately for both 
    !   cluster types (e.g., MG and QI+QS).  Overlaps are assumed if at least 
    !   half of the lower area is covered by the cluster aloft.  
    !   This condition is checked 
    !   by testing whether occupied pixels are found above the centroid of the
    !   lower cluster.  If true, the information about the lower and upper 
    !   clusters are stored in the ovl_clusters structure.
    !
    ! Input:
    ! 
    !   As input, this subroutine accepts the binary fields
    !   
    !   * top   : upper cluster
    !
    !   * bottom: lower cluster  
    ! 
    !--------------------------------------------------------------------------

    SUBROUTINE cluster_overlaps (top, bottom, bottom_lab)

      !------------------------------------------------------------------------
      ! Subroutine arguments
      !------------------------------------------------------------------------

      IMPLICIT NONE

      INTRINSIC :: MAX, INT

      ! Input arguments
      ! ---------------

      INTEGER, INTENT(IN) :: top (ie,je,ke)   ! binary field of upper cluster
      INTEGER, INTENT(IN) :: bottom(ie,je,ke) ! binary field of lower cluster

      ! Output arguments
      ! ----------------

      INTEGER  :: bottom_lab (ie,je,ke)  ! Labeled field

      !------------------------------------------------------------------------
      ! Local variables
      !------------------------------------------------------------------------

      ! Derived data types

      TYPE(cluster_info), ALLOCATABLE  ::     &
           info_top   (:),                    &   
           info_bottom(:)

      ! 3D-arrays
      ! --------

      REAL (dp)  :: &
                                !!$not usedrho_glob_temp (ie_tot, je_tot, ke_tot),    &
                                !!$not usedrho_glob      (ie_tot, je_tot, ke_tot),    &
           rhfl_glob_red (ie_tot, je_tot, ke_tot),    &
           rhfl_glob     (ie_tot, je_tot, ke_tot),    &
           qg_glob_temp  (ie_tot, je_tot, ke_tot),    &
           qg_glob       (ie_tot, je_tot, ke_tot),    &
           rhfl          (ie,je,ke),                  &  ! full levels
           graupel_max_arr(ie_tot, je_tot, ke_tot)

      INTEGER  ::                   &
           top_lab             (ie,je,ke),           &
                                ! non-reduced labeled upper field
           top_lab_glob        (ie_tot, je_tot, ke), & 
                                ! non-reduced labeled upper field
           bottom_lab_glob     (ie_tot, je_tot, ke), & 
                                ! all-reduced labeled upper field
           bottom_lab_glob_red (ie_tot, je_tot, ke), & 
                                ! all-reduced labeled upper field
           top_lab_glob_red    (ie_tot, je_tot, ke)    

      ! 2D-arrays
      ! ---------

      ! None

      ! 1D-arrays
      ! ---------

      INTEGER ::  &
           y_glob(je), &   ! result of i/j_glob function
           tcsize(csize_length)  , &
           bcsize(csize_length)

      INTEGER,  ALLOCATABLE ::                 &
           bottom_centroid_x(:),               &
           bottom_centroid_y(:),               &
           bottom_centroid_z(:),               &   
           overlap_top_ut(:),                  &  ! Contains overlapping top
           overlap_bot_ut(:)                      ! Contains overlapping bottom

      ! Local scalars
      ! -------------

      INTEGER ::   &
           ierror,                    & ! error flags, logical unit
           tcs_len,                   & ! length of top csize
           bcs_len,                   & ! length of mottom csize 
           i, j, k, ii,               & ! loop variables
           xind, yind, zind,          & ! indices for centroid positions
           size_out,                  & ! size overlaps utility 
           icount,                    & ! counter
           vals,                      & ! buffer size (elements) for allreduce
           iglob, jglob,              & ! global indizes
           izx, izy,                  & ! utility position indices
           xix, xiy, xiz,             & ! utility position indices
           search_range                 ! for area where max is searched


      REAL (dp)  :: grau_max      ! Maximum of 1D slice thru QG centroid

      !========================================================================
      ! Start routine
      !========================================================================

      !------------------------------------------------------------------------
      ! Section 1: Label clusters and analyze their properties
      !------------------------------------------------------------------------

      ! Label bottom clusters and find centroid positions

      bottom_lab = 0
      cs_count   = 0
      bcsize     = 0

      CALL label (bottom, bottom_lab, bcsize)

      bcs_len = cs_count

      ALLOCATE(info_bottom       (bcs_len))
      ALLOCATE(bottom_centroid_x (bcs_len))
      ALLOCATE(bottom_centroid_y (bcs_len))
      ALLOCATE(bottom_centroid_z (bcs_len))

      IF (bcs_len > 0) THEN

         CALL cluster_analysis(info_bottom, bottom_lab, bcsize)

         DO ii = 1, bcs_len
            bottom_centroid_x(ii) = info_bottom(ii) % cent_pos_x
            bottom_centroid_y(ii) = info_bottom(ii) % cent_pos_y
            bottom_centroid_z(ii) = info_bottom(ii) % cent_pos_z
         ENDDO

         ! Label top clusters

         top_lab = 0

         CALL label (top, top_lab, tcsize)
         tcs_len = cs_count

         ALLOCATE(info_top(tcs_len))
         IF (tcs_len > 0) THEN
            CALL cluster_analysis(info_top, top_lab, tcsize)
         ENDIF

         size_out = MAX(bcs_len, tcs_len)    ! overlap, utility

         ALLOCATE(overlap_top_ut(size_out))
         ALLOCATE(overlap_bot_ut(size_out))

         overlap_top_ut = 0
         overlap_bot_ut = 0

         ! The labeled field (upper plate) is all-reduced to all processors,
         ! and hence of dimension (ie_tot, je_tot, ke).  This way, it can
         ! be accessed with global indices.  Doing it the other way round,
         ! i.e., finding the local (i,j) and the respective local PE, of
         ! the bottom centroids would be utterly inefficient.

         rhfl(:,:,1:ke) = 0.5_dp * (hhl(:,:,2:ke+1) + hhl(:,:,1:ke))

!!$not usedrho_glob_temp   = 0.0_dp
         rhfl_glob       = 0.0_dp
         top_lab_glob    = 0
         bottom_lab_glob = 0
         qg_glob_temp    = 0.0_dp

         y_glob = global_j

         DO k = 1, ke
            DO j = jstartpar, jendpar
               DO i = istartpar, iendpar
                  iglob = i
                  jglob = y_glob(j)
                  top_lab_glob   (iglob, jglob, k) = top_lab   (i,j,k)
                  bottom_lab_glob(iglob, jglob, k) = bottom_lab(i,j,k)
                  rhfl_glob      (iglob, jglob, k) = rhfl      (i,j,k)
                  qg_glob_temp   (iglob, jglob, k) = qg_3d     (i,j,k) 
!!$not usedrho_glob_temp  (iglob, jglob, k) = rho       (i,j,k)
               ENDDO
            ENDDO
         ENDDO

!!$not used           rho_glob            = 0.0_dp
         top_lab_glob_red    = 0
         bottom_lab_glob_red = 0
         rhfl_glob_red       = 0.0_dp  
         qg_glob             = 0.0_dp

         vals    = ie_tot * je_tot * ke_tot

         CALL global_fields(top_lab_glob, top_lab_glob_red, vals, 'SUM'  &
              , yerrmsg, ierror)

!!$            CALL MPI_ALLREDUCE                                       &
!!$                (top_lab_glob, top_lab_glob_red, vals,                &
!!$                imp_integers, MPI_SUM, icomm_cart, ierror)

         CALL global_fields(bottom_lab_glob, bottom_lab_glob_red, vals  &
              , 'SUM', yerrmsg, ierror)
!!$           CALL MPI_ALLREDUCE                                       &
!!$                (bottom_lab_glob, bottom_lab_glob_red, vals,          &
!!$                imp_integers, MPI_SUM, icomm_cart, ierror)

         CALL global_fields(rhfl_glob, rhfl_glob_red, vals  &
              , 'MAX', yerrmsg, ierror)
!!$          CALL MPI_ALLREDUCE                                       &
!!$                (rhfl_glob, rhfl_glob_red, vals,                        &
!!$                imp_reals, MPI_MAX, icomm_cart, ierror)

!!$not used           CALL global_fields(rho_glob_temp, rho_glob, vals  &
!!$not used                , 'MAX', yerrmsg, ierror)
!!$           CALL MPI_ALLREDUCE                                       &
!!$                (rho_glob_temp, rho_glob, vals,                         &
!!$                imp_reals, MPI_MAX, icomm_cart, ierror)

         CALL global_fields(qg_glob_temp,qg_glob , vals  &
              , 'MAX', yerrmsg, ierror)
!!$           CALL MPI_ALLREDUCE                                       &
!!$                (qg_glob_temp, qg_glob, vals,                           &
!!$                imp_reals, MPI_MAX, icomm_cart, ierror)

         ! Find label of upper pair member

         icount = 0
         cs_loop: DO ii = 1, bcs_len
            xind = bottom_centroid_x(ii)
            yind = bottom_centroid_y(ii)
            zind = bottom_centroid_z(ii)
            k_loop: DO k = zind, 1, -1  ! going upward
               IF (top_lab_glob_red(xind, yind, k) /= 0) THEN 
                  icount = icount + 1
                  overlap_top_ut(icount) = top_lab_glob_red(xind, yind, k)
                  overlap_bot_ut(icount) = info_bottom(ii) % label
                  EXIT k_loop
               ENDIF
            ENDDO k_loop
         ENDDO cs_loop

         number_overlaps = icount 

         ! Check if size of info_cap structure is sufficient

         IF (number_overlaps > max_number_cap) THEN
            IF (p_parallel_io) THEN
               print *, ' *** SRC_LIGHTNING.CLUSTER_OVERLAPS ERROR ***'
               print *, &
                    'CLUSTER_OVERLAPS: Capacitor structure assigned too few elements.'
               print *, 'Increase MAX_NUMBER_CAP to at least', number_overlaps
            ENDIF
            CALL error_bi(' ', ' ')
            !              CALL MPI_ABORT(MPI_COMM_WORLD, 100, ierror)
         ENDIF

         ! Obtain structure contents and store them in arrays for better access 
         ! (vectorization)

         ! Initialize and assign capacitor structure

         info_cap % cap_label      = 0
         info_cap % top_label      = 0
         info_cap % top_pixels     = 0
         info_cap % top_area       = 0.0_dp
         info_cap % top_depth      = 0.0_dp
         info_cap % separation     = 0.0_dp
         info_cap % plate_distance = 0.0_dp
         info_cap % breakdown_alt  = 0.0_dp
         info_cap % x_pos          = 0
         info_cap % y_pos          = 0
         info_cap % bot_label      = 0
         info_cap % bot_pixels     = 0
         info_cap % bot_area       = 0.0_dp
         info_cap % bot_depth      = 0.0_dp
         info_cap % graupel        = 0.0_dp

         DO ii = 1, number_overlaps
            IF (overlap_top_ut(ii) /= 0) THEN
               info_cap(ii) % cap_label  = ii
               info_cap(ii) % top_label  = overlap_top_ut(ii) 
               info_cap(ii) % top_pixels = &
                    info_top(overlap_top_ut(ii)) % pixels  
               info_cap(ii) % x_pos      = &
                    info_bottom(overlap_bot_ut(ii)) % cent_pos_x
               info_cap(ii) % y_pos      = &
                    info_bottom(overlap_bot_ut(ii)) % cent_pos_y 
               info_cap(ii) % bot_label  = &
                    overlap_bot_ut(ii)
               info_cap(ii) % bot_pixels = &
                    info_bottom(overlap_bot_ut(ii)) % pixels
            ENDIF
         ENDDO

         ! Calculate "strength" of lower plate (GRAR, QG,...)

         DO ii = 1, number_overlaps

            izx = info_cap(ii) % x_pos
            izy = info_cap(ii) % y_pos 

            search_range    = 2
            graupel_max_arr = 0.0_dp
            grau_max        = 0.0_dp

            IF (izx > search_range .AND. xind < ie - search_range .AND. &
                 izy > search_range .AND. yind < je_tot - search_range) THEN
               DO k = 1, ke
                  DO j = izy-search_range, izy+search_range
                     DO i = izx-search_range, izx+search_range  
                        ! Adjust here for density
                        graupel_max_arr(i,j,k) = rho(i,j,k) * qg_glob(i,j,k)  
                     ENDDO
                  ENDDO
               ENDDO
               grau_max = MAXVAL(graupel_max_arr(:,:,:))
            ENDIF

            info_cap(ii) % graupel = graupel_mass_correction * grau_max

         ENDDO

         ! Determine the breakdown altitude, i.e., the altitude between the 
         ! plates. This is just the average of the upper and lower z-centroid 
         ! positions. Though the centroids may be strongly offset laterally, 
         ! this estimate should be OK.

         DO ii = 1, number_overlaps
            xiz = INT(0.5_dp * REAL(info_top(ii)%cent_pos_z + info_bottom(ii)%cent_pos_z,dp) )
            xiy = info_cap(ii) % y_pos
            xix = info_cap(ii) % x_pos
            info_cap(ii) % breakdown_alt = rhfl_glob_red(xix, xiy, xiz)   
         ENDDO

         ! Complete entries of info_cap structure

         CALL capacitor_details  (info_top, info_bottom, &
              top_lab_glob_red, bottom_lab_glob_red,     &
              number_overlaps, rhfl_glob_red, ke         &
              , ie_tot, je_tot, ke_tot, dlon &
              , dlat, startlat_tot, p_parallel_io)

      ELSEIF (bcs_len == 0) THEN

         info_cap % cap_label     = 0
         info_cap % top_label     = 0
         info_cap % top_pixels    = 0
         info_cap % top_area      = 0.0_dp
         info_cap % top_depth     = 0.0_dp
         info_cap % separation    = 0.0_dp
         info_cap % x_pos         = 0
         info_cap % y_pos         = 0
         info_cap % bot_label     = 0
         info_cap % bot_pixels    = 0
         info_cap % bot_area      = 0.0_dp
         info_cap % bot_depth     = 0.0_dp
         info_cap % graupel       = 0.0_dp
         info_cap % breakdown_alt = 0.0_dp

      ENDIF

      !-----------------------------------------------------------------------
      ! Deallocate dynamic arrays
      !------------------------------------------------------------------------

      DEALLOCATE(bottom_centroid_x)
      DEALLOCATE(bottom_centroid_y)
      DEALLOCATE(bottom_centroid_z)
      DEALLOCATE(info_top         )
      DEALLOCATE(info_bottom      )
      DEALLOCATE(overlap_top_ut   )
      DEALLOCATE(overlap_bot_ut   )

      !=======================================================================
      ! End subroutine cluster_overlaps
      !=======================================================================

    END SUBROUTINE cluster_overlaps

    !=========================================================================
    ! Subroutine cluster_analysis
    !=========================================================================

    !--------------------------------------------------------------------------
    ! Description:
    !   This module procedure calculates the vertical and horizontal extent of
    !   the different clusters, their volume, and their centroid positions.
    !
    ! Method:
    !   The centroid is simply the arithmetic mean of the positions of the
    !   occupied gridopints.  The extents (and volume) are proportional to
    !   the number of occupied gridpoints.
    !
    ! Input:
    !
    !   * Properly and consecutively labeled field
    !   * final (consecutive) CSIZE array
    !   * Optionally, the field upon which BIN_FIELD function acted
    !
    ! Output:
    !   Structure containing arrays with cluster parameters
    !
    !--------------------------------------------------------------------------

    SUBROUTINE cluster_analysis (cluster, labeled, lcsize) 

      !------------------------------------------------------------------------
      ! Subroutine arguments
      !------------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER, INTENT(IN)   ::   &
           labeled(ie,je,ke),    &   ! properly/consecutively labeled field
           lcsize(:)                 ! csize


      ! Derived types
      !--------------

      TYPE(cluster_info), INTENT(OUT)  :: cluster(cs_count)

      ! Local parameters
      !-----------------

      ! None
      INTEGER ::  lcs_len    ! number of clusters
      INTEGER ::  ierror 
      ! Local scalars
      !--------------

      INTEGER   ::   &
           i, j, k, ii,                    &   ! loop variables
           iz1, iz2, iz3,                  &   ! loop variables
           lerror,                  &
           top_x_loc, top_y_loc,           &   ! ij_local output  
           bottom_x_loc, bottom_y_loc,     &   ! ij_local output
           lprocess                            !     - " -  

      ! Local arrays
      !-------------

      INTEGER  ::    &
           top_xx(cs_count),         & ! x-position of uppermost pixel
           top_yy(cs_count),         &
           top_zz(cs_count),         &
           top_x_glob(cs_count),     &
           top_y_glob(cs_count),     & 
           top_z_glob(cs_count),     & ! Globally reduced (MIN (= max height)) 
           bottom_xx(cs_count),      &   ! x-position of lowermost pixel
           bottom_yy(cs_count),      &   !
           bottom_zz(cs_count),      &   !
           bottom_x_glob(cs_count),  &   ! Globally-reduced
           bottom_y_glob(cs_count),  &
           bottom_z_glob(cs_count),  &
           jglobind(je)

      REAL (dp)           ::    &
           cent_pos_x_util(cs_count),       & 
           cent_pos_y_util(cs_count),       &
           cent_pos_z_util(cs_count),       &
           cent_pos_x_arr_glob(cs_count),   & ! Globally reduced (sum)
           cent_pos_y_arr_glob(cs_count),   &
           cent_pos_z_arr_glob(cs_count),   &
           height_top(cs_count),            & ! Utility array (= top_height)
           height_top_glob(cs_count),       &
           height_bottom(cs_count),         & ! Utility array (= bottom_height)
           height_bottom_glob(cs_count)

      REAL (dp) :: hfl(ie,je,ke)           ! Full levels

      ! For performance debugging

!!$        REAL (dp) ::                            &
!!$             ts_reduce_pos1, ts_reduce_pos2, ts_reduce_pos, &
!!$             ts_bottom1, ts_bottom2, ts_bottom
!!$
      ! Dynamic arrays
      !---------------

      ! None

      LOGICAL ::  verbose = .FALSE. 

      !========================================================================
      ! Begin procedure
      !========================================================================

      !------------------------------------------------------------------------
      ! Section 1.1: Find cluster centroids
      !------------------------------------------------------------------------

      ! Initializations (necessary as not all processors contain all clusters)

      lcs_len = cs_count

      hfl(:,:,1:ke) = 0.5 * (hhl(:,:,2:ke+1) + hhl(:,:,1:ke))

      top_xx = 0
      top_yy = 0

      jglobind = global_j

      ! Array elements are first determined by each processes separately
      ! and subsequently summed up in a global reduction operation.

      DO ii = 1, lcs_len
         iz1 = 0
         iz2 = 0
         iz3 = 0
         DO k = 1, ke
            DO j = jstartpar, jendpar
               DO i = istartpar, iendpar
                  IF (labeled(i,j,k) == ii) THEN
                     iz1 = iz1 + i
                     iz2 = iz2 + jglobind(j)   
                     iz3 = iz3 + k
                  ENDIF
               ENDDO
            ENDDO
         ENDDO

         cent_pos_x_util(ii) =  REAL(iz1,dp) / REAL(lcsize(ii),dp)
         cent_pos_y_util(ii) =  REAL(iz2,dp) / REAL(lcsize(ii),dp)
         cent_pos_z_util(ii) =  REAL(iz3,dp) / REAL(lcsize(ii),dp)

      ENDDO

      ! CALL MPI_BARRIER (icomm_cart, lerror)

      !       ts_reduce_pos1 = MPI_WTIME()

      CALL global_fields(cent_pos_x_util,cent_pos_x_arr_glob , lcs_len  &
           , 'SUM', yerrmsg, ierror)

!!$        CALL MPI_ALLREDUCE                            &
!!$             (cent_pos_x_util, cent_pos_x_arr_glob, lcs_len, &
!!$             imp_reals, MPI_SUM, icomm_cart, lerror)

      CALL global_fields(cent_pos_y_util, cent_pos_y_arr_glob, lcs_len  &
           , 'SUM', yerrmsg, ierror)

!!$        CALL MPI_ALLREDUCE                            &
!!$             (cent_pos_y_util, cent_pos_y_arr_glob, lcs_len, &
!!$             imp_reals, MPI_SUM, icomm_cart, lerror)


      CALL global_fields(cent_pos_z_util, cent_pos_z_arr_glob, lcs_len  &
           , 'SUM', yerrmsg, ierror)
!!$        CALL MPI_ALLREDUCE                            &
!!$             (cent_pos_z_util, cent_pos_z_arr_glob, lcs_len, &
!!$             imp_reals, MPI_SUM, icomm_cart, lerror)

      !       ts_reduce_pos2 = MPI_WTIME()

      !        ts_reduce_pos = ts_reduce_pos2 - ts_reduce_pos1

      !------------------------------------------------------------------------
      ! Section 1.2: Find bottom height of clusters
      !------------------------------------------------------------------------

      !ts_bottom1 = MPI_WTIME()

      bottom_zz     = 0
      bottom_x_glob = 0
      bottom_y_glob = 0

      DO ii = 1, lcs_len
         vertical_loop: DO k = ke, 1, -1    ! From bottom to top
            DO j = jstartpar, jendpar
               DO i = istartpar, iendpar
                  IF (labeled(i,j,k) == ii) THEN
                     bottom_zz(ii) = k
                     EXIT vertical_loop 
                  ENDIF
               ENDDO
            ENDDO
         ENDDO vertical_loop
      ENDDO

      CALL global_fields(bottom_zz, bottom_z_glob, lcs_len  &
           , 'MAX', yerrmsg, ierror)
!!$     CALL MPI_ALLREDUCE (bottom_zz, bottom_z_glob, lcs_len, imp_integers,    &
!!$             MPI_MAX, icomm_cart, lerror)

      ! Find horizontal coordinates of cluster bottom

      bottom_xx = 0
      bottom_yy = 0

      DO ii = 1, lcs_len
         k_loop: DO k = ke, 1, -1   ! From bottom to top
            DO j = jstartpar, jendpar
               DO i = istartpar, iendpar
                  IF (labeled(i,j,k) == ii .AND. k == bottom_z_glob(ii)) THEN
                     bottom_xx(ii) = i
                     bottom_yy(ii) = global_j(j)
                     EXIT k_loop
                  ENDIF
               ENDDO
            ENDDO
         ENDDO k_loop
      ENDDO

      ! Let other processes know

      CALL global_fields(bottom_xx, bottom_x_glob, lcs_len  &
           , 'MAX', yerrmsg, ierror)
!!$    CALL MPI_ALLREDUCE (bottom_xx, bottom_x_glob, lcs_len, imp_integers,   &
!!$             MPI_MAX, icomm_cart, lerror)

      CALL global_fields(bottom_yy ,bottom_y_glob, lcs_len  &
           , 'MAX', yerrmsg, ierror)
!!$     CALL MPI_ALLREDUCE (bottom_yy, bottom_y_glob, lcs_len, imp_integers,  &
!!$             MPI_MAX, icomm_cart, lerror)

!!$        ts_bottom2 = MPI_WTIME()
!!$        ts_bottom = ts_bottom2 - ts_bottom1

      !------------------------------------------------------------------------
      ! Section 1.2.1: Calculate the geometric heights of the clusters' bottoms
      !------------------------------------------------------------------------

      height_bottom = 0._dp

      DO ii = 1, lcs_len
         CALL ij_local                                    & 
              (bottom_x_glob(ii), bottom_y_glob(ii), bottom_x_loc,   &
              bottom_y_loc, lprocess, lerror)
         IF (p_pe == lprocess) THEN
            height_bottom(ii) = &
                 hfl(bottom_x_loc, bottom_y_loc, bottom_z_glob(ii))
         ENDIF
      ENDDO

      CALL global_fields(height_bottom ,  height_bottom_glob, lcs_len  &
           , 'MAX', yerrmsg, ierror)
      ! CALL MPI_ALLREDUCE (height_bottom, height_bottom_glob, lcs_len, imp_reals, &
      !           MPI_MAX, icomm_cart, lerror)

      !------------------------------------------------------------------------
      ! Section 1.3: Find highest point occupied by individual clusters
      !------------------------------------------------------------------------

      ! Initialize top_z with a large number (GT ke)

      top_zz     = 300
      top_x_glob = 0
      top_y_glob = 0

      DO ii = 1, lcs_len
         k_loop2: DO k = 1, ke
            DO j = jstartpar, jendpar
               DO i = istartpar, iendpar
                  IF (labeled(i,j,k) == ii) THEN
                     top_zz(ii) = k
                     EXIT k_loop2
                  ENDIF
               ENDDO
            ENDDO
         ENDDO k_loop2
      ENDDO

      CALL global_fields(top_zz, top_z_glob, lcs_len  &
           , 'MIN', yerrmsg, ierror)

      !$        CALL MPI_ALLREDUCE (top_zz, top_z_glob, lcs_len, imp_integers,       &
!!$             MPI_MIN, icomm_cart, lerror)                   

      ! Find horizontal coordinates of cluster top

      DO ii = 1, lcs_len
         k_loop3:  DO k = 1, ke
            DO j = jstartpar, jendpar
               DO i = istartpar, iendpar
                  IF (labeled(i,j,k) == ii .AND. k == top_z_glob(ii)) THEN
                     top_xx(ii) = i
                     top_yy(ii) = global_j(j)
                     EXIT k_loop3
                  ENDIF
               ENDDO
            ENDDO
         ENDDO k_loop3
      ENDDO

      ! Let other processes know

      CALL global_fields(top_xx, top_x_glob, lcs_len  &
           , 'MAX', yerrmsg, ierror)
!!$     CALL MPI_ALLREDUCE (top_xx, top_x_glob, lcs_len, imp_integers,         &
!!$             MPI_MAX, icomm_cart, lerror)

      CALL global_fields(top_yy, top_y_glob, lcs_len  &
           , 'MAX', yerrmsg, ierror)
!!$     CALL MPI_ALLREDUCE (top_yy, top_y_glob, lcs_len, imp_integers,         &
!!$           MPI_MAX, icomm_cart, lerror)

      !-----------------------------------------------------------------------
      ! Section 1.3.1: Calculate the geometric heights of the clusters' tops
      !-----------------------------------------------------------------------

      height_top = 0._dp

      DO ii = 1, lcs_len
         CALL ij_local                                    & 
              (top_x_glob(ii), top_y_glob(ii), top_x_loc,   &
              top_y_loc, lprocess, lerror)
         IF (p_pe == lprocess) THEN
            height_top(ii) = hfl(top_x_loc, top_y_loc, top_z_glob(ii))
         ENDIF
      ENDDO

      CALL global_fields(height_top, height_top_glob, lcs_len  &
           , 'MAX', yerrmsg, ierror)
!!$     CALL MPI_ALLREDUCE (height_top, height_top_glob, lcs_len, imp_reals,   &
!!$             MPI_MAX, icomm_cart, lerror)

      !------------------------------------------------------------------------
      ! Section 1.4: Calculate average value of cluster property (like mean 
      !              upward velocity in w-cluster)
      !------------------------------------------------------------------------

      cluster(ii) % average = 0.0_dp

      !  in case of MESSY, only itype_light==1      IF (itype_light == 5) THEN 

!!$        cluster_sum_loc = 0.0_dp
!!$
!!$        DO ii = 1, lcs_len
!!$           DO k = 1, ke
!!$              DO j =  jstartpar, jendpar
!!$                 DO i = istartpar, iendpar
!!$                    IF (labeled(i,j,k) == ii) THEN
!!$                       cluster_sum_loc(ii) = &
!!$                            cluster_sum_loc(ii) + element_values(i,j,k) 
!!$                    ENDIF
!!$                 ENDDO
!!$              ENDDO
!!$           ENDDO
!!$        ENDDO
!!$
!!$        cluster_sum_glob = 0.0_dp
!!$
!!$           CALL global_fields(cluster_sum_loc, cluster_sum_glob, lcs_len  &
!!$                , 'SUM', yerrmsg, ierror)
!!$! CALL MPI_ALLREDUCE (cluster_sum_loc, cluster_sum_glob, lcs_len, imp_reals,&
!!$!                MPI_SUM, icomm_cart, lerror)
!!$
!!$        DO ii = 1, lcs_len
!!$           cluster(ii) % average = cluster_sum_glob(ii) / REAL(lcsize(ii),dp)
!!$        ENDDO

      ! ENDIF

      !---------------------------------------------------------------------
      ! Section 2: Pass contents to cluster structure
      !---------------------------------------------------------------------

      DO ii = 1, lcs_len

         cluster(ii) % label         = ii
         cluster(ii) % pixels        = lcsize(ii)
         cluster(ii) % cent_pos_x    = NINT(cent_pos_x_arr_glob(ii))
         cluster(ii) % cent_pos_y    = NINT(cent_pos_y_arr_glob(ii))
         cluster(ii) % cent_pos_z    = NINT(cent_pos_z_arr_glob(ii))
         cluster(ii) % top_x         = top_x_glob(ii)
         cluster(ii) % top_y         = top_y_glob(ii)
         cluster(ii) % top_z         = top_z_glob(ii)
         cluster(ii) % bottom_height = height_bottom_glob(ii)
         cluster(ii) % top_height    = height_top_glob(ii)

      ENDDO

      IF (verbose) THEN
         print *, p_pe, 'Time step is: ', current_time_step
         DO ii = 1, lcs_len
            print *, p_pe, 'label: ',    &
                 cluster(ii) % label,      &
                 cluster(ii) % pixels,     &
                 cluster(ii) % cent_pos_x, &
                 cluster(ii) % cent_pos_y, &
                 cluster(ii) % cent_pos_z, & 
                 cluster(ii) % top_x,      &
                 cluster(ii) % top_y,      &
                 cluster(ii) % top_z,      & 
                 cluster(ii) % bottom_height, & 
                 cluster(ii) % top_height  
         ENDDO
      ENDIF

      !=======================================================================
      ! End module procedure cluster_analysis
      !=======================================================================

    END SUBROUTINE cluster_analysis

    !==========================================================================
    ! Subroutine LABEL
    !==========================================================================

    !--------------------------------------------------------------------------
    !
    ! Description:
    !
    ! Identify and label contiguous regions of properties defined in bin_field
    ! function.  The algorithm is based on Hoshen and Kopelman (1976).  This
    ! procedure accepts input of an integer array containing 0's and -1's,
    ! where the clusters of adjacent -1's are labeled.  The NEWS-neighborhood
    ! rule is used to determine pixel adjacency.  Parallelization has been
    ! realized as in Constantin et al. (supercomp. appl., 1997)
    !
    ! Method:
    !
    ! This program identifies and labels clusters using the NEWS neighborhood
    ! rule with the Hoshen-Kopelman algorithm.  Its fast processing time is
    ! owed to the maintainence of a separate array, CSIZE, which contains
    ! information about the number of pixels in one cluster, as well as cluster
    ! coalescence.  Initially, the pixels are given temporary labels, which
    ! are adjusted according to CSIZE information in a second pass through
    ! the array.
    !
    ! The array is traversed row by row, and previously-labeled pixels are
    ! checked every time an occupied site is encountered.
    !
    ! CSIZE is a 1D array, initialized to contain all clusters (checkerboard
    ! distribution).  Every cluster label is used as index to a CSIZE element.
    ! This element either contains the number of pixels in the cluster, or
    ! a negative number.  The absolute value of this number is the index
    ! to the cluster which the current cluster belongs to.  It may point
    ! to another CSIZE element which is negative, but eventually a positive
    ! number is reached, representing the number of elements in this cluster
    ! (and of all the clusters pointing to this cluster).
    ! The pointer path is found iteratively by the function PROPLAB (for
    ! "proper label").
    ! CSIZE is updated every time an occupied site is encountered.
    !
    ! In the parallel implementation, CSIZE is declared as array that
    ! has 1E6 elements.
    ! CSIZE is defined for all processes such as to be able to contain all
    ! clusters of the TOTAL domain;  however, each process uses only a 
    ! well-defined part of CSIZE.
    !--------------------------------------------------------------------------

    SUBROUTINE label (field_in, field_out, csize_consec) 

      IMPLICIT NONE        

      INTRINSIC :: COUNT, INT, MAX, MIN
      !------------------------------------------------------------------------
      ! Subroutine arguments
      !------------------------------------------------------------------------

      INTEGER, INTENT(IN)  ::  field_in    (ie,je,ke)

      INTEGER, INTENT(OUT) ::  field_out   (ie,je,ke)                  

      INTEGER, INTENT(OUT) ::  csize_consec(csize_length)  

      ! Local constants and parameters; utility variables
      !---------------------------------------------------

      INTEGER, PARAMETER   ::  dcount = 1

      INTEGER  :: iejeh, i, j, k, ii,  &         !  Loop variables
           incr

      CHARACTER (LEN=25)        ::  yerrmsg   ! for MPI error message

      ! Local scalar variables
      !-----------------------

      INTEGER   :: cs_size, num_clust,  &    ! size of csize;
           proper!,                     &    ! number of clusters
      ! lun                               ! unit number for write

      INTEGER   ::  ccount, counter   ! COUNT function output,
      ! counter variable

      ! Naming convention: "nb" stands for neighbor; 
      !                    "t" for target (of pointer path)

      INTEGER ::        &
           s,                             &
           min_n, max_n,                    &
           target_min, target_max,          &
           target_n1, target_n2, target_n3, &  ! targets of neighbors
           n1, n2, n3,                      &  ! neighbors
           position, eval, oval,            &  ! output from
                                ! two_neighbors
           target_eval, target_oval,        &  ! their targets
           min_2nb, max_2nb,                &  ! min, max neighbors
           positiont,                       &  ! output from two_n.
           evalt, ovalt,                    &  !     -- " --
           min_neighbor,                    &  ! min neighbor
           min3_index,                      &  ! index of -"-
           onb1, onb2,                      &  ! non-minumum neighb.
           mint3nb, onb1tar, onb2tar,       &  ! targets
           position2nb, eval2nb, oval2nb,   &  ! output two_elements
           target1nb,                       &  ! one neighb. target
           nb2eq, target_nb2eq,             &  ! neighbor if both
                                ! identical; its target
           min_2nbt, max_2nbt,              &  ! min, 2 neighbors; targ.
           ierror                              ! general and MPI error code 

      INTEGER  ::        &
           s_ovlap   , n_ovlap   ,          &  ! Overlapping clusters (north and south)
           loclen    ,                      &  ! "local" length of csize
           count_fe  , count_se  ,          &  ! Ordered-pair elements of merge_array
           merlen    ,                      &  ! length of merge_array 
           maxol     ,                      &  ! Maximum of overlapping clusters
                                ! per row
           ope1    , ope2,                  &  ! ordered-pair elements 1 and 2
           count_feu, count_seu,            &  !
           target_north, target_south

      !------------------------------------------------------------------------
      ! Local arrays
      !------------------------------------------------------------------------

      ! 1D-arrays
      ! ---------

      INTEGER  :: &   
           nb_targets(3),             & 
           neighbors (3),             &
           kzdims   (24),             &            ! Vertical dimension of sendbuf variables
           op        (2)                           ! ordered pair  

      ! 3D-arrays
      ! ---------

      REAL (DP) :: field_prel(ie,je,ke)   ! prelim. field; field_out has 
      ! wrong data type (INT) for MPI
      INTEGER ::    matrix(ie+1,je+1,ke+1)

      !------------------------------------------------------------------------
      ! Dynamic arrays
      !------------------------------------------------------------------------

      INTEGER, ALLOCATABLE :: &
           csize     (:),                        &
           csize_glob(:),                        &          
           csize_work(:),                        &
           merge_array(:),                       &  ! Local merge_array        
           utility_ma (:),                       &  ! Utility for merge array
           merge_array_glob(:)                   ! Globally-reduced version

      !------------------------------------------------------------------------
      ! Logical variables
      !------------------------------------------------------------------------

      LOGICAL                               :: &
           printdb = .FALSE.,                  & ! if true, print cluster info 
           flag, flagt, flag2nb

      !-----------------------------------------------------------------------
      ! Allocate memory 
      !-----------------------------------------------------------------------

      iejeh  = csize_length 
      maxol  = n_overlapping_elements
      loclen = NINT(REAL(iejeh,dp) / REAL(p_nprocs,dp))  
      merlen = p_nprocs * maxol

      IF (p_parallel_io .AND. INT(p_nprocs) * INT(loclen) < INT(iejeh)) THEN
         print *, '***********************************************************'
         print *, 'SRC_LIGHTNING.LABEL WARNING:' 
         print *, 'Sum of local CSIZE segments bigger than entire CSIZE array!'
         print *, 'nproc, loclen, nproc * loclen, iejeh:', &
              p_nprocs, loclen, p_nprocs * loclen, iejeh
         print *, '***********************************************************'
      ENDIF

      ALLOCATE (csize           (iejeh))
      ALLOCATE (csize_glob      (iejeh) )  ! csize after global reduction
      ALLOCATE (utility_ma      (merlen))
      ALLOCATE (merge_array     (merlen))  ! Merge_Array
      ALLOCATE (merge_array_glob(merlen))

      !========================================================================
      ! Begin program
      !========================================================================

      !------------------------------------------------------------------------
      ! Section 1:
      !
      ! Assign working array "matrix" (m+1)x(n+1)x(k+1) matrix
      ! such that a semi-halo around the upper-left border as well as at the 
      ! top of the cube, is created.
      !------------------------------------------------------------------------

      matrix = 0
      DO k = 1, ke
         DO j = jstartpar, jendpar
            DO i = istartpar, iendpar  
               matrix(i+1,j+1,k+1) = field_in(i,j,k) 
            ENDDO
         ENDDO
      ENDDO

      !------------------------------------------------------------------------
      ! Section 1.1:
      !
      ! Define processor-dependent variables used for parallel computations
      ! the appropriate segmentation of CSIZE would imply that incr starts from
      !------------------------------------------------------------------------

      incr = p_pe * loclen    

      csize = 0

      !------------------------------------------------------------------------
      ! Section 2:
      ! Traverse field and label it.  No cluster-fragment linking is
      ! performed during the first pass and has to be done in a second pass 
      ! with the aid of csize.
      !------------------------------------------------------------------------

      DO k = 2, ke+1                ! First the uppermost horizontal slice
         DO j =  2, je+1            ! row-wise traversing
            DO i = 2, ie+1

               IF (matrix(i,j,k) == -1) THEN

                  ! Some self-explanatory utility variables

                  n1        = matrix(i-1,j,k)         ! western neighbor
                  n2        = matrix(i,j-1,k)         ! northern neighbor
                  n3        = matrix(i,j,k-1)         ! upper neighbor
                  neighbors = (/n1, n2, n3/)
                  min_n     = MIN(n1, n2, n3)
                  max_n     = MAX(n1, n2, n3)

                  !------------------------------------------------------------
                  ! Case 1: No neighbors.
                  !         Simply assign new label given by the counter
                  !------------------------------------------------------------

                  IF (n1 == 0 .AND. n2 == 0 .AND. n3 == 0) THEN
                     s = incr + dcount     

                     matrix(i,j,k) = s                
                     csize(s)      = 1            
                     incr          = incr + dcount

                     !--------------------------------------------------------
                     ! Case 2: All three neighbors labeled previously.
                     !         Altogether, 6 cases to be considered
                     !--------------------------------------------------------

                  ELSEIF (n1 > 0 .AND. n2 > 0 .AND. n3 > 0) THEN

                     ! Some utility assignments

                     matrix(i,j,k) = min_n         ! Assign temporary label

                     CALL PROPLAB(csize, n1, target_n1, status)
                     IF (status /= 0) CALL error_bi( &
                          'error in PROPLAB call for target_n1', substr)
                     CALL PROPLAB(csize, n2, target_n2, status)
                     IF (status /= 0) CALL error_bi( &
                          'error in PROPLAB call for target_n2', substr)
                     CALL PROPLAB(csize, n3, target_n3, status)
                     IF (status /= 0) CALL error_bi( &
                          'error in PROPLAB call for target_n3', substr)
                     CALL PROPLAB(csize, min_n, target_min, status)
                     IF (status /= 0) CALL error_bi( &
                          'error in PROPLAB call for target_min', substr)

                     nb_targets    = (/target_n1, target_n2, target_n3/)

                     ! Attempt, to call later - otherwise it's called for every
                     ! site, and that will slow down the program ...

                     CALL TWO_ELEMENTS(neighbors, flag, position, eval, oval)
                     CALL TWO_ELEMENTS(nb_targets, flagt, positiont, evalt  &
                          , ovalt)

                     !---------------------------------------------------------
                     ! Case 2.1: All neighbors are different
                     !---------------------------------------------------------

                     IF (n1 /= n2 .AND. n2 /= n3 .AND. n1 /= n3) THEN

                        !------------------------------------------------------
                        ! Case 2.1.1: All targets equal
                        !------------------------------------------------------

                        IF (target_n1 == target_n2 .AND. &
                             target_n2 == target_n3) THEN

                           csize(target_n1) = 1 + csize(target_n1)

                           !---------------------------------------------------
                           ! Case 2.1.2: Two of the three targets equal
                           !--------------------------_------------------------

                        ELSEIF (flagt) THEN

                           ! Find position of minimum

                           min_neighbor = MIN(n1, n2, n3)

                           ccount = 0
                           DO ii = 1, 3
                              ccount = ccount + 1
                              IF (neighbors(ii) == min_neighbor) &
                                   min3_index = ccount
                           ENDDO

                           ! Find position of other elements

                           IF (min3_index == 1) THEN
                              onb1 = 2                   ! other neighbor 1
                              onb2 = 3                   ! other neighbor 2
                           ELSEIF (min3_index == 2) THEN
                              onb1 = 1
                              onb2 = 3
                           ELSEIF (min3_index == 3) THEN
                              onb1 = 1
                              onb2 = 2
                           ENDIF

                           ! target of minimum neighbor

                           CALL PROPLAB(csize, min_neighbor, mint3nb, status)
                           IF (status /= 0) CALL error_bi(&
                                'error in PROPLAB call for mint3nb', substr)

                           ! targets of other beighbors
                           CALL PROPLAB(csize, neighbors(onb1),onb1tar, status)
                           IF (status /= 0) CALL error_bi(&
                                'error in PROPLAB call for onb1tar', substr)
                           CALL PROPLAB(csize, neighbors(onb2),onb2tar, status)
                           IF (status /= 0) CALL error_bi(&
                                'error in PROPLAB call for onb2tar', substr)

                           IF (mint3nb == onb1tar .AND. &
                                onb1tar /= onb2tar) THEN
                              csize(mint3nb) = &
                                   1 + csize(mint3nb) + csize(onb2tar)
                              csize(onb2tar) = -mint3nb
                           ELSEIF (mint3nb == onb2tar .AND. &
                                onb2tar /= onb1tar) THEN
                              csize(mint3nb) = &
                                   1 + csize(mint3nb) + csize(onb1tar)
                              csize(onb1tar) = -mint3nb
                           ELSEIF (onb1tar == onb2tar .AND. &
                                onb1tar /= mint3nb) THEN
                              csize(mint3nb) = &
                                   1 + csize(mint3nb) + csize(onb1tar)
                              csize(onb1tar) = -mint3nb
                              csize(onb2tar) = -mint3nb
                           ENDIF

                           !--------------------------------------------------
                           ! Case 2.1.3: All three targets are different
                           !--------------------------------------------------

                        ELSEIF (target_n1 /= target_n2 .AND. &
                             target_n2 /= target_n3    .AND. &
                             target_n1 /= target_n3)   THEN

                           ! Find position of minimum

                           min_neighbor = MIN(n1, n2, n3)

                           ccount = 0
                           DO ii = 1, 3
                              ccount = ccount + 1
                              IF (neighbors(ii) == min_neighbor) &
                                   min3_index = ccount
                           ENDDO

                           ! Find position of other elements

                           IF (min3_index == 1) THEN
                              onb1 = 2                   ! other neighbor 1
                              onb2 = 3                   ! other neighbor 2
                           ELSEIF (min3_index == 2) THEN
                              onb1 = 1
                              onb2 = 3
                           ELSEIF (min3_index == 3) THEN
                              onb1 = 1
                              onb2 = 2
                           ENDIF

                           ! Do cluster-fragment linking

                           CALL PROPLAB(csize, min_neighbor, mint3nb, status)
                           IF (status /= 0) CALL error_bi(&
                                'error in PROPLAB call for mint3nb (b)', substr)

                           CALL PROPLAB(csize, neighbors(onb1),onb1tar, status)
                           IF (status /= 0) CALL error_bi(&
                                'error in PROPLAB call for onb1tar (b)', substr)
                           CALL PROPLAB(csize, neighbors(onb2),onb2tar, status)
                           IF (status /= 0) CALL error_bi(&
                                'error in PROPLAB call for onb2tar (b)', substr)

                           csize(mint3nb) =  1 + csize(mint3nb) &
                                + csize(onb1tar) +  csize(onb2tar)
                           csize(onb1tar) = - mint3nb
                           csize(onb2tar) = - mint3nb

                        ENDIF       ! Check of number of different targets

                        !-----------------------------------------------------
                        ! Case 2.2: Two of three neighbors are different
                        !-----------------------------------------------------

                     ELSEIF (flag) THEN

                        CALL PROPLAB(csize, oval,target_oval , status)
                        IF (status /= 0) CALL error_bi(&
                             'error in PROPLAB call for target_oval(b)', substr)
                        CALL PROPLAB(csize, eval, target_eval, status)
                        IF (status /= 0) CALL error_bi(&
                             'error in PROPLAB call for target_eval(b)', substr)

                        !------------------------------------------------------
                        ! Case 2.2.1: Both targets are equal
                        !---------------------_--------------------------------

                        IF (target_oval == target_eval) THEN
                           csize(target_oval) = 1 + csize(target_oval)

                           !---------------------------------------------------
                           ! Case 2.2.2: Both targets are different
                           !---------------------------------------------------

                        ELSEIF (target_oval /= target_eval) THEN

                           CALL PROPLAB(csize,min_n ,min_2nbt , status)
                           IF (status /= 0) CALL error_bi(&
                                'error in PROPLAB call for min_2nbt (b)', substr)
                           CALL PROPLAB(csize,max_n ,max_2nbt , status)
                           IF (status /= 0) CALL error_bi(&
                                'error in PROPLAB call for max_2nbt (b)', substr)

                           csize(min_2nbt) = &
                                1 + csize(min_2nbt) + csize(max_2nbt)
                           csize(max_2nbt) = -min_2nbt

                        ENDIF

                        !-----------------------------------------------------
                        ! Case 2.3: All neighbors are equal 
                        !           (and hence, all targets)
                        !-----------------------------------------------------

                     ELSEIF (n1 == n2 .AND. n2 == n3) THEN
                        csize(target_n1) = 1 + csize(target_n1)

                     ENDIF   ! How many neighbors of the three are equal

                     !--------------------------------------------------------
                     ! Case 3: Two neighbors out of three are occupied 
                     !         (and labeled)
                     !
                     ! Two cases are checked: Two neighbors having been 
                     !                        labeled previously and one 
                     !                        neighbor having been labeled
                     !                        previously. The latter case 
                     !                        implies that
                     !                        there is only one neighbor.
                     !--------------------------------------------------------

                  ELSEIF ( n1 > 0 .AND. n2 > 0 .AND. n3 == 0 .OR.         &
                       n1 > 0 .AND. n3 > 0 .AND. n2 == 0 .OR.         &
                       n2 > 0 .AND. n3 > 0 .AND. n1 == 0 ) THEN

                     !--------------------------------------------------------
                     ! Case 3.1: Both neighbors are different
                     !--------------------------------------------------------

                     ! Call of two_elements to obtain logical flag:
                     ! If .FALSE., all three values are different; since one of
                     ! them
                     ! is zero, it directs the program to the desired branch.

                     CALL TWO_ELEMENTS(neighbors, flag2nb, position2nb &
                          , eval2nb, oval2nb)

                     IF (.NOT. flag2nb) THEN

                        max_2nb    = MAX(n1, n2, n3)

                        ! Temporarily overwrite the zero to find non-trivial
                        ! minimum

                        WHERE (neighbors == 0) neighbors = 1000000

                        n1 = neighbors(1)
                        n2 = neighbors(2)
                        n3 = neighbors(3)

                        min_2nb = MIN(n1, n2, n3)

                        matrix(i,j,k) = min_2nb

                        ! Recover neighbors

                        WHERE (neighbors == 1000000) neighbors = 0

                        n1 = neighbors(1)
                        n2 = neighbors(2)
                        n3 = neighbors(3)

                        CALL PROPLAB(csize,min_2nb , target_min , status)
                        IF (status /= 0) CALL error_bi(&
                             'error in PROPLAB call for target_min (c)', substr)
                        CALL PROPLAB(csize,max_2nb , target_max , status)
                        IF (status /= 0) CALL error_bi(&
                             'error in PROPLAB call for target_max (c)', substr)

                        !-----------------------------------------------------
                        ! Case 3.1.1: Both targets are different
                        !-----------------------------------------------------

                        IF (target_max /= target_min) THEN
                           csize(target_min) = &
                                1 + csize(target_min) + csize(target_max)
                           csize(target_max) = - target_min 

                           !--------------------------------------------------
                           ! Case 3.1.2: Both targets are equal
                           !--------------------------------------------------

                        ELSEIF (target_max == target_min) THEN
                           csize(target_min) = 1 + csize(target_min)
                        ENDIF

                        !-----------------------------------------------------
                        ! Case 3.2: Both neighbors (and hence, targets)
                        !           are identical
                        !-----------------------------------------------------

                     ELSE    ! flag2nb true, i.e., both neighbors are identical
                        nb2eq            = MAX(n1, n2, n3)
                        ! output from two_neighbors is (0, n, n)
                        matrix(i,j,k)    = nb2eq  

                        CALL PROPLAB(csize, nb2eq , target_nb2eq , status)
                        IF (status /= 0) CALL error_bi(&
                             'error in PROPLAB call for target_nb2eq', substr)
                        csize(target_nb2eq) = 1 + csize(target_nb2eq)
                     ENDIF

                     !---------------------------------------------------------
                     ! Case 4: Only one neighbor previously labeled
                     !---------------------------------------------------------

                  ELSEIF (n1 == 0 .AND. n2 == 0 .AND. n3 > 0 .OR.          &
                       n1 == 0 .AND. n3 == 0 .AND. n2 > 0 .OR.          &
                       n2 == 0 .AND. n3 == 0 .AND. n1 > 0) THEN

                     matrix(i,j,k) = max_n

                     ! Find proper label
                     CALL PROPLAB(csize, max_n , target1nb , status)
                     IF (status /= 0) CALL error_bi(&
                          'error in PROPLAB call for target1nb', substr)
                     csize(target1nb) = 1 + csize(target1nb)

                  ENDIF     ! Number of neighbor
               ENDIF       ! site occupied

            ENDDO            ! i-loop
         ENDDO
      ENDDO

      !------------------------------------------------------------------------
      ! Section 3.1:  
      !  Do the proper labeling. If CSIZE has no non-zero elements, no further
      !  statistics is attempted and the number of clusters is set to zero.
      !------------------------------------------------------------------------

      cs_size   = COUNT(csize /= 0)
      num_clust = COUNT(csize > 0)

      IF (cs_size == 0 .AND. printdb) &
           print *, 'Process ', p_pe, 'LABEL: NO CLUSTERS FOUND'

      IF (cs_size > 0) THEN
         ALLOCATE(csize_work(cs_size))
         csize_work = 0
         csize_work = &
              csize(1+(p_pe*loclen) : p_pe*loclen+cs_size)

         ! Proper labeling

         DO i = 1+(p_pe*loclen), p_pe*loclen+cs_size   
            IF (csize(i) < 0) THEN

               CALL PROPLAB(csize, i , proper , status)
               IF (status /= 0) CALL error_bi(&
                    'error in PROPLAB call for proper', substr)
               WHERE (matrix == i) matrix = proper
            ENDIF
         ENDDO

      ENDIF  ! csize containing non-trivial elements

      ! Assign 3D field (ie,je,ke) with contents of matrix

      field_prel(1:ie,1:je,1:ke) = REAL(matrix(2:ie+1,2:je+1,2:ke+1),dp)

      !-----------------------------------------------------------------------
      ! Section 4: Exchange boundaries and set up merge_array ordered-pair 
      !            array 
      !            Update global CSIZE
      !  
      ! Method: The upper and lower three rows are used as halos so that they 
      !         are occupied by the values of the respective neighbors. Then, 
      !         the row at jendpar+1 is compared with the row at j = jendpar.
      !         Where
      !         clusters overlap (i.e. "touch" each other), csize is updated.
      !
      !------------------------------------------------------------------------

      ! Assign kzdims

      kzdims(:) = 0
      kzdims(1) = ke 

      ! Exchange values; halo depth is 3 (nboundlines) grid points.

      CALL exchange_boundaries(kzdims, field_prel,  ierror, yerrmsg)
!!$        CALL exchg_boundaries                                               &
!!$             (nnow+39, sendbuf, isendbuflen, imp_reals, icomm_cart, ie, je, &
!!$         kzdims, jstartpar, jendpar,nboundlines, nboundlines, my_cart_neigh,&
!!$         21000+ntstep, .FALSE., ncomm_type, izerror, yerrmsg,               &
!!$             field_prel(:,:,:))

      IF (ierror /= 0) THEN
         yerrmsg = 'LABEL: MPI EXCHANGE ERROR' 
         RETURN
      ENDIF

      !-----------------------------------------------------------------------
      ! Section 4.1:
      !
      ! * Spread local CSIZE contents to all processes via MPI global reduction
      !   routine.
      ! 
      ! * Find processor that covers northern edge of domain to limit the loop
      !   to those processors that have a northern neighbor.  
      !
      ! * Scan through the edge rows, look for overlapping clusters, 
      !   and set up merge_array (ordered neighbor pairs)
      !-----------------------------------------------------------------------

      ! Distribute entire CSIZE array to all processes

      !csize_glob = 0

      CALL global_fields(csize, csize_glob, iejeh  &
           , 'SUM', yerrmsg, ierror)
      !        CALL MPI_ALLREDUCE (csize, csize_glob, iejeh, imp_integers &
      !             , MPI_SUM, icomm_cart, izerror)       

      ! Determine processor at northern boundary edge (nedgeproc)


      !um_ak_20110606 +
      ! -3 includes an assumption about nboundlines, why not simply
      !    use ie_tot, je_tot
      ! anyway: this does not work, as this only yields the processor
      !  at the edge, but for the following you need to exclude all
      !  processors at the northern boundary
      !CALL ij_local (ie_tot-3, je_tot-3, i_loc, j_loc, nedgeproc, ierror)
      !        CALL ij_local (ie_tot, je_tot, i_loc, j_loc, nedgeproc, ierror)

      !        use meaningful criteria: see below
      !um_ak_20110606 -

      ! Set up Merge_Array

      merge_array = 0
      count_fe    = p_pe * maxol + 1 
      count_se    = p_pe * maxol + 2

      !um_ak_20110606 IF (p_pe /= nedgeproc) THEN
      IF (je /= jendpar) THEN
         DO  k = 1, ke
            DO i = 1, ie
               n_ovlap = INT(field_prel(i,jendpar+1,k))    
               s_ovlap = INT(field_prel(i,jendpar  ,k)) 
               IF ( s_ovlap  /= 0 .AND. n_ovlap /= 0 ) THEN
                  merge_array(count_fe:count_se) = (/ n_ovlap, s_ovlap /)
                  count_fe = count_fe + 2
                  count_se = count_se + 2
               ENDIF
            ENDDO
         ENDDO
      ENDIF

      ! Scan through merge_array in parallel and eliminate redundant entries
      ! within the local processor scopes 

      utility_ma = 0 
      count_feu  = p_pe * maxol + 1
      count_seu  = p_pe * maxol + 2

      DO i = p_pe * maxol + 1, (p_pe+1) * maxol, 2 
         ope1 = merge_array(i)
         ope2 = merge_array(i+1)
         op = (/ ope1, ope2 /) 
         IF (ope1 /= 0 .AND. ope2 /= 0) THEN
            DO ii = p_pe * maxol + 1, (p_pe+1) * maxol, 2
               IF (merge_array(ii)   == ope1 .AND. &
                    merge_array(ii+1) == ope2) THEN
                  merge_array(ii)                  = 0
                  merge_array(ii+1)                = 0
                  utility_ma (count_feu:count_seu) = op
               ENDIF
            ENDDO
            count_feu = count_feu + 2
            count_seu = count_seu + 2
         ENDIF
      ENDDO

      ! Broadcast result to all processes

      !        merge_array_glob = 0

      CALL global_fields(utility_ma, merge_array_glob, merlen  &
           , 'SUM', yerrmsg, ierror)
!!$        CALL MPI_ALLREDUCE                             &
!!$             (utility_ma, merge_array_glob, merlen,      &
!!$             imp_integers, MPI_SUM, icomm_cart, izerror)

      counter = 0
      DO ii = 1, merlen
         IF (merge_array_glob(ii) /= 0) THEN
            counter = counter + 1
            IF (printdb) print *, p_pe, 'global merge array: ',  &
                 counter, merge_array_glob(ii)
         ENDIF
      ENDDO

      !------------------------------------------------------------------------
      ! Section 4.2: Update global csize
      !  The usual cluster-fragment linking is performed, but the global matrix
      !  (field_prelim) and merge_array fields are updated right away.  This
      !  avoids PROPLAB-Calls and the problem of linking multiply connected
      !  regions, some of which are properly labeled while others are not.
      !------------------------------------------------------------------------

      DO i = 1, merlen-1, 2
         n_ovlap = merge_array_glob(i)
         s_ovlap = merge_array_glob(i+1)

         IF (n_ovlap /= 0 .AND. s_ovlap /= 0 .AND. s_ovlap == n_ovlap) THEN 
            IF (printdb) THEN
               print *, ''
               print *, 'Previous merging has occurred', &
                    s_ovlap, n_ovlap
            ENDIF
         ELSEIF (n_ovlap /= 0 .AND. s_ovlap /= 0 .AND. n_ovlap /= s_ovlap) THEN

            IF (printdb) THEN
               print *, p_pe, 'MERGE_ARRAY ORIGINAL: ' &
                    , merge_array_glob(i), merge_array_glob(i+1)
               print *, ''
               print *, p_pe, 'Overlaps', n_ovlap, s_ovlap 
            ENDIF

            target_north = n_ovlap 
            target_south = s_ovlap

            IF (printdb) &
                 print *, p_pe, 'before: ', csize_glob(target_north) &
                 , csize_glob(target_south)

            csize_glob(target_north) = &
                 csize_glob(target_north) + csize_glob(target_south)
            csize_glob(target_south) = -target_north

            IF (printdb) &
                 print *, p_pe, 'after: ', csize_glob(target_north) &
                 , csize_glob(target_south) 

            WHERE (field_prel == REAL(target_south,dp)) &
                 field_prel = REAL(target_north)
            WHERE (merge_array_glob == target_south) &
                 merge_array_glob = target_north     

            IF (printdb) THEN
               print *, p_pe, 'Global field relabeling: ', &
                    target_south, ' becomes ', target_north
               print *, p_pe, 'MERGE_ARRAY MODIFIED: ', &
                    merge_array_glob(i), merge_array_glob(i+1)
               print *, ''
            ENDIF

         ENDIF
      ENDDO

      !------------------------------------------------------------------------
      ! Section 4.3: 
      !  Two final passes for proper and consecutive labeling with the aid of 
      !  the global CSIZE array
      !------------------------------------------------------------------------

      ! Do consecutive labeling for entire domain ...

      cs_count     = 0
      csize_consec = 0

      DO i = 1, iejeh
         IF (csize_glob(i) > 0) THEN
            cs_count = cs_count + 1
            csize_consec(cs_count) = csize_glob(i)
            WHERE (field_prel == REAL(i,dp)) field_prel = REAL(cs_count,dp)
         ENDIF
      ENDDO

      field_out = INT(field_prel)

      !-----------------------------------------------------------------------
      ! Section 5: Some cluster statistics
      !---------------------------------_-------------------------------------

      num_clust = 0

      IF (cs_size > 0) THEN

         IF (printdb) THEN
            print *, 'Local statistics for process number ', p_pe
            print *, ''
         ENDIF

         num_clust = COUNT(csize > 0)

         IF (printdb) THEN
            DO i = 1+p_pe * loclen, p_pe* loclen+cs_size
               print *, p_pe, i, csize(i), csize_glob(i)
            ENDDO

            print *, ''
            print *, p_pe, ' Process  Cluster  number of elements'
            print *, '-------------------------------------------------'

            DO i = 1, cs_size  
               IF (csize_work(i) > 0) print '(3i8)', i, csize_work(i)
            ENDDO
            print *, p_pe,  'Number of clusters:', num_clust
            print *, ''
         ENDIF              ! verbose output
      ENDIF              ! csize containing nontrivial elements  

      IF (printdb) THEN
         print *, p_pe, 'Global statistics at timestep:', current_time_step
         print *, ''

         print *, p_pe, ' Cluster       number of elements'
         print *,            '----------------------------------'

         DO i = 1, cs_count  
            print '(3i8)', p_pe, i, csize_consec(i)
         ENDDO
      ENDIF

      !------------------------------------------------------------------------
      ! Section 6: Deallocate dynamic arrays
      !-----------------------_------------------------------------------------

      IF (cs_size > 0) THEN
         DEALLOCATE (csize_work  )
      ENDIF

      DEALLOCATE (csize           )
      DEALLOCATE (csize_glob      )
      DEALLOCATE (utility_ma      )
      DEALLOCATE (merge_array     )
      DEALLOCATE (merge_array_glob)

      !========================================================================
      ! End subroutine label
      !========================================================================

    END SUBROUTINE LABEL

    !=========================================================================
    ! Subroutine distribute_flashes
    !=========================================================================

    !-------------------------------------------------------------------------
    ! Description:
    !   Based on the location of the cell centroid (i,j) as well as their
    !   instantaneous lightning frequencies, three arrays are created, 
    !   containing flash locations (lon, lat in rotated coordinates) as well 
    !   as the times of their occurrence.   The spatial distribution is 
    !   circularly-symmetric about the cell centroid
    !   position.  Upstream distribution has been attempted but yielded no
    !   satisfactory results.  The optimum seems to be a circular distribution
    !   where the radius corresponds to the width of the cell.
    !
    ! Method:
    !   A random-number generator is used to distribute the flashes. For the 
    !   spatial distribution, a Gauss-weighing is applied, while the temporal
    !   distribution is purely random.
    !  
    !   The time step as well as the interval when the lightning package are 
    !   called, determine the time interval across which the flashes have to
    !   be distributed.  The instantaneous flash rate is given in s^(-1), so
    !   the number of flashes to be distributed is simply the calling interval
    !   (in steps) times the step size times the flash frequency.
    !
    !    Those cells that exhibit too weak a lightning rate to result in
    !    at least one lightning per period are included but their lon/lat
    !    couplet is set to -999.  These values must be filtered by the
    !    post-processing software. 
    !--------------------------------------------------------------------------

    SUBROUTINE distribute_flashes (dcell, n_tot, lon, lat, time)

      !------------------------------------------------------------------------
      ! Subroutine arguments
      !------------------------------------------------------------------------

      IMPLICIT NONE

      INTRINSIC :: COS, EXP, INT, NINT, REAL, SIN, SQRT

      ! INTENT(IN)

      ! (cs_count) Assumed size; either cs_count or number_overlaps
      TYPE(storm), INTENT(IN)    :: dcell(:) 

      INTEGER                    :: n_tot  ! total number of strikes per
      ! hinclight

      ! INTENT(OUT)

      REAL (dp), INTENT(OUT)      ::  &   ! Assumed shape from call
           lon (:),                               &
           lat (:),                               &
           time(:)

      !------------------------------------------------------------------------
      ! Local variables
      !------------------------------------------------------------------------

      ! Tempo variables of various kinds for debugging purposes
      ! -------------------------------------------------------

      ! Local parameters
      ! ----------------

      ! 3 GPs width of lightning path: 0.075 deg
      REAL (dp), PARAMETER :: min_width  = 0.17_dp 
      !  width_min  = 0.17 looks good for some reason   

      ! Local scalars
      ! -------------

      ! DECL  (mark)

      REAL (dp)        ::  &
           !             sinclight,           & ! hinclight in seconds
           flash_freq,          &
           deg2rad, rad2deg,    &
           cent_x_angl,         & ! x-position of centroid in COSMO coordinates
           cent_y_angl,         & ! y-position of centroid in COSMO coordinates
           llsigma, weighfac,   & ! coefficients of weighing function
           equiv_diam_km          ! Equivalent circular cell diameter

      INTEGER ::  &
           i, ii,              &      ! loop variables
           lun, lc_count

      INTEGER  ::  &
                                ! current time and time at beginning of interval
           t_current, t_start,          & 
           n_cell,                     & ! flashes per cell and interval
                                ! start/end indices for output array concatenation 
           arr_start, arr_end

      ! Local arrays
      ! ------------

      ! 1D arrays
      ! ---------

      INTEGER  ::  &
           n_cell_arr(cs_count+1)   ! array containing the number of flashes per 
      ! cell per time interval (contains 1 as 
      ! first element)

      ! Local dynamic arrays
      ! --------------------

      INTEGER, ALLOCATABLE  :: &
           argument(:),               &     ! independent variables for Gauss
           li_cell_arr(:),            &     ! same as n_cell_arr but w/o zeros
           xpos_cell_arr(:),          &     ! x-coordinates of flashing cells
           ypos_cell_arr(:)                 ! y-coordinates of flashing cells

      REAL (dp), ALLOCATABLE  :: &
           weighing(:),                &  ! Gauss function
           dist(:),                    &  ! radial distance 
           angl(:),                    &  ! azimuth angle
           weighted(:),                &  ! weighted distance  
           x_loc(:),                   &  ! x-coord of flash rel to centroid
           y_loc(:),                   &  ! y-coord of flash
           x_loc_final(:),             &  ! x-coord in COSMO coordinates
           y_loc_final(:),             &  ! y-coord in COSMO coordinates
           yy_ut(:),                   &  ! Utility field for random numbers
           pwidth(:),                        &
           pwidth_flashing(:),               & 
           steering_alt(:)

      ! Characters
      ! ----------
      CHARACTER (LEN=5) :: time_arr

      LOGICAL           :: printdb = .FALSE.   


      INTRINSIC :: RANDOM_SEED, RANDOM_number

      !======================================================================
      ! Start procedure
      !======================================================================

      ! itype == 1 DAHL
      !        IF (itype_light == 1 .AND. cs_count < number_overlaps .AND. my_cart_id == 0) THEN
      IF (cs_count < number_overlaps .AND. p_pe == 0) THEN
         print *, 'SRC_LIGHTNING.DISTRIBUTE_FLASHES ERROR'
         print *,'Inconsistent number of cells:', cs_count, number_overlaps
         print *, 'Accum. flashes, time step', n_tot, current_time_step
      ENDIF

      !IF (itype_light == 1) THEN
      cs_count = number_overlaps      ! Adjust trip counts for 2-class cell
      !ENDIF                             ! global var from LABEL overwritten

      deg2rad = pi/180.0_dp
      rad2deg = 180.0_dp/pi

      !sinclight = lightning_step * time_step_len ! calling interval in seconds
      ! seconds after initialization
      t_current = current_time_step * INT(time_step_len)
      t_start   = t_current - INT(sinclight)

      !-----------------------------------------------------------------------
      ! Section 1.1: Determine the width of cells
      !-----------------------------------------------------------------------

      ALLOCATE (pwidth(cs_count)); pwidth = 0.0_dp

      !-----------------------------------------------------------------------
      ! Section 1.2: Prepare fields needed for determining storm motion 
      !-----------------------------------------------------------------------

      !notused hfl(:,:,1:ke) = 0.5_dp * (hhl(:,:,1:ke) + hhl(:,:,2:ke+1))

      !-----------------------------------------------------------------------
      ! Section 2: Continue serially on processor #0 - distribute flashes
      !-----------------------------------------------------------------------

      IF (p_pe == 0) THEN

         IF (n_tot == 0) THEN
            ! Arrays have two elements (as allocated above)
            time = -999.0_dp       
            lon  = -999.0_dp
            lat  = -999.0_dp
         ELSEIF (n_tot > 0) THEN

            !------------------------------------------------------------------
            ! Section 2.1:
            !  Create times of occurrence in seconds after initialization
            !------------------------------------------------------------------

            CALL RANDOM_SEED
            CALL RANDOM_NUMBER (time)

            time = REAL(t_start,dp) + (sinclight * time)

            !-----------------------------------------------------------------
            ! Section 2.2:
            ! Distribute flashes spatially
            !-----------------------------------------------------------------

            ! Array containing the number of flashes per call interval and 
            ! cell are calculated in a separate loop 

            n_cell_arr    = 0
            DO ii = 1, cs_count
               flash_freq = dcell(ii) % flash_rate
               !                 n_cell =   NINT( flash_freq * 3600.0_dp * hinclight )
               n_cell =   NINT( flash_freq * sinclight )
               n_cell_arr(ii+1) = n_cell
            ENDDO

            ! Initialize with dummy values

            lon = -999.0_dp
            lat = -999.0_dp

            ! Create new n_cell array and cs_count, containing only those cells
            ! that actually produce lightning.  First loop is to find length
            ! of lightning-cell array, second loop is for assignment.

            lc_count = 0
            DO ii = 2, cs_count+1
               IF (n_cell_arr(ii) /= 0) THEN
                  lc_count = lc_count + 1
               ENDIF
            ENDDO

            ALLOCATE (li_cell_arr(lc_count+1)  )
            ALLOCATE (xpos_cell_arr(lc_count)  )
            ALLOCATE (ypos_cell_arr(lc_count)  )
            ALLOCATE (pwidth_flashing(lc_count))
            ALLOCATE(steering_alt(lc_count)    )

            li_cell_arr   = 0
            xpos_cell_arr = 0
            ypos_cell_arr = 0

            lc_count    = 0

            DO ii = 2, cs_count+1
               IF (n_cell_arr(ii) /= 0) THEN
                  lc_count = lc_count + 1
                  li_cell_arr(lc_count+1) = n_cell_arr(ii)
                  xpos_cell_arr(lc_count) = dcell(ii-1) % centroid_x        
                  ypos_cell_arr(lc_count) = dcell(ii-1) % centroid_y
                  steering_alt(lc_count)  = REAL(zrs(ii-1),dp)
               ENDIF
            ENDDO

            arr_start = 1
            arr_end   = 0

            ! Loop over flashing cells
            ! ------------------------

            DO ii = 1, lc_count

               n_cell = li_cell_arr(ii+1) 

               ALLOCATE(argument   (n_cell))
               ALLOCATE(weighing   (n_cell))
               ALLOCATE(dist       (n_cell)) 
               ALLOCATE(angl       (n_cell))
               ALLOCATE(weighted   (n_cell))
               ALLOCATE(x_loc      (n_cell))
               ALLOCATE(y_loc      (n_cell))
               ALLOCATE(x_loc_final(n_cell))
               ALLOCATE(y_loc_final(n_cell))
               ALLOCATE(yy_ut(3 * n_cell))

               !--------------------------------------------------------------
               ! Section 2.4:
               ! Distribute flashes in a circle around the cell      
               !---------------------------------------------------------------

               ! Diameter of the cell

               equiv_diam_km = 2.0_dp * SQRT(1.E-6 * info_cap(ii) % bot_area / pi) 

               ! Form-preserving weighing (within [0, n_cell] interval)
               ! 0 <= weighing(x) < 1 for x eps [0, ncell]

               llsigma    = 0.4_dp * REAL(n_cell,dp)
               weighfac = - 1.0_dp / (2.0_dp * (llsigma * llsigma))

               ! Create data in polar coordinates and thereafter transform into
               ! Cartesian grid

               ! AK here itype_light always = 1
!!$                 IF (itype_light == 1) THEN
               pwidth_flashing(ii) = 0.5_dp * equiv_diam_km / &    
                                ! Radius 
                    (r_earth * 1.0E-3_dp) * rad2deg 
!!$                 ELSEIF (itype_light > 1) THEN
!!$                    pwidth_flashing(ii) = 0.05_dp + 2.5E-4_dp * n_cell 
!!$                 ENDIF

               IF (pwidth_flashing(ii) <= min_width) &
                    pwidth_flashing(ii) = min_width

               IF (printdb) WRITE (*,*) 'TRACK_WIDTH', ii, pwidth_flashing(ii)

               CALL RANDOM_SEED
               CALL RANDOM_NUMBER (dist)

               dist = pwidth_flashing(ii) * dist 

               ! Create uncorrelated random-number series

               CALL RANDOM_SEED
               CALL RANDOM_NUMBER (yy_ut)

               angl(1:n_cell) = yy_ut(3:n_cell+2)

               angl = 2.0_dp * pi * angl      ! in radians

               ! Set up Gauss function; starting from x = 0

               argument = 0
               DO i = 2, n_cell
                  argument(i) =  argument(i) + (i-1)
               ENDDO
               argument(1) =  0

               weighing = EXP(weighfac * (REAL(argument) * REAL(argument)))

               weighted = dist * weighing

               ! Cartesian coordinates

               x_loc = weighted * COS(angl)
               y_loc = weighted * SIN(angl)

               cent_x_angl = startlon_tot + dlon *  REAL(xpos_cell_arr(ii) - 1 , dp) 
               cent_y_angl = startlat_tot + dlat *  REAL(ypos_cell_arr(ii) - 1 , dp)

               x_loc_final = x_loc + cent_x_angl
               y_loc_final = y_loc + cent_y_angl

               arr_start = arr_start + li_cell_arr(ii)
               arr_end   = arr_end   + li_cell_arr(ii+1)

               ! Concatenate into final array

               lon(arr_start:arr_end) = x_loc_final
               lat(arr_start:arr_end) = y_loc_final 

               DEALLOCATE(argument   )
               DEALLOCATE(weighing   )
               DEALLOCATE(dist       )
               DEALLOCATE(angl       )
               DEALLOCATE(weighted   )
               DEALLOCATE(x_loc      )
               DEALLOCATE(y_loc      )
               DEALLOCATE(x_loc_final)
               DEALLOCATE(y_loc_final)
               DEALLOCATE(yy_ut      )

            ENDDO    ! loop over flashing cells
         ENDIF    ! ntot > 0
      ENDIF    ! my_cart_id == 0  

      !------------------------------------------------------------------------
      ! Section 3:
      !   Write output as ASCII
      !   files:
      !
      !   The target directory is: 
      !     /e/uhome/dlrdahl/lightning_output/
      !
      !   File-name format:
      !     lightning_sinclight.txt
      ! 
      !------------------------------------------------------------------------

      IF (p_pe == 0) THEN

!!$           IF (t_current >= 100 .AND. t_current < 1000) THEN
!!$              WRITE(time_arr,'(i3)') INT(t_current)
!!$              time_arr = '00'//time_arr
!!$           ELSEIF (t_current >= 1000 .AND. t_current < 10000) THEN
!!$              WRITE(time_arr,'(i4)') INT(t_current)
!!$              time_arr = '0'//time_arr
!!$           ELSEIF(t_current >= 10000 .AND. t_current < 100000) THEN
!!$              WRITE(time_arr,'(i5)') INT(t_current)
!!$              !time_arr = time_ut_arr3
!!$           ENDIF
         CALL int2str(time_arr, t_current)

         lun = find_next_free_unit(100,200)
         !CALL get_free_unit(lun)

         OPEN &
              (unit = lun, file = &
              'lightning_'// time_arr //'.txt')
         WRITE (lun, *)  '       COSMO Lightning Data' 
         WRITE (lun, *)  'TIME           LON          LAT' 
         WRITE (lun, *)  '-----------------------------------'
         WRITE (lun, *)

         IF (n_tot > 0) THEN 
            DO i = 1, n_tot ! cs_count
               WRITE(lun, '(3f12.5)') time(i), lon(i), lat(i)
            ENDDO
         ELSE IF (n_tot == 0) THEN
            DO i = 1, 2
               WRITE(lun, '(3f12.5)') time(i), lon(i), lat(i)
            ENDDO
         ENDIF
         CLOSE (lun)

      ENDIF

      !------------------------------------------------------------------------
      ! Deallocate dynamic arrays
      !------------------------------------------------------------------------

      DEALLOCATE(steering_alt    )
      DEALLOCATE (pwidth         )
      DEALLOCATE (li_cell_arr    )
      DEALLOCATE (xpos_cell_arr  ) 
      DEALLOCATE (ypos_cell_arr  )
      DEALLOCATE (pwidth_flashing)

      !========================================================================
      ! End subroutine distribute_flashes
      !========================================================================

    END SUBROUTINE distribute_flashes

!#endif

  END SUBROUTINE lnox_dahl2000_global_end
  ! ========================================================================
!#endif
! CSOMOv...
#endif
! DAHL2000

  ! ***********************************************************************
END MODULE messy_lnox_si
! ***********************************************************************


#include "messy_main_ppd_bi.inc"

! **********************************************************************
! INTERFACE TO ECHAM5 FOR
! Radon module for BL/convection/scavenging diagnostic
! Version: see 'modver' in messy_dradon.f90
!
! Author : Patrick Joeckel, MPICH, July  2002
!
! References:
!
! **********************************************************************

! **********************************************************************
MODULE messy_dradon_si
! **********************************************************************

  ! ECHAM5/MESSy
  USE messy_main_blather_bi,    ONLY: start_message_bi, end_message_bi
  ! MESSy
  USE messy_main_constants_mem, ONLY: STRLEN_MEDIUM
  USE messy_main_channel,       ONLY: t_chaobj_cpl
  USE messy_dradon
#ifdef MESSYTENDENCY
 !tendency budget
 USE messy_main_tendency_bi,   ONLY: mtend_get_handle,       &
                                     mtend_get_start_l,      &
                                     mtend_add_l,            &
                                     mtend_register,           &    
                                     mtend_id_tracer
#endif

  IMPLICIT NONE
  PRIVATE
  SAVE
  INTRINSIC :: NULL

  ! ... CHANNEL OBJECTS: (DIAGNOSED FIELDS)
  ! Rn_222 flux density
  REAL(DP), DIMENSION(:,:),   POINTER :: rnflux       => NULL()
  ! delta-p at surface
  REAL(DP), DIMENSION(:,:),   POINTER :: dp_sfc       => NULL()
  ! Rn_222 tendency
  REAL(DP), DIMENSION(:,:),   POINTER :: rnflux_vmrs  => NULL()

#ifdef ECHAM5
  ! LAGRANGIAN
  ! Rn_222 flux
  REAL(DP), DIMENSION(:),     POINTER :: lg_rnflux_vmrs   => NULL()
  ! Rn_222 flux
  REAL(DP), DIMENSION(:,:),   POINTER :: rnflux_vmrs_rest => NULL()
#endif
  ! GLOBAL COUPLING SWITCHES
  ! GRIDPOINT
  LOGICAL :: L_GP               = .TRUE.  ! GRIDPOINT CALCULATION
  INTEGER :: I_GP_emis_method   = 1
  ! ... chain calculation
  LOGICAL :: L_GP_chain         = .FALSE. ! 222Rn -> ... -> 210Pb
  CHARACTER(LEN=STRLEN_MEDIUM) :: C_GP_210Pb_aermod  = ''  ! aerosol model
  INTEGER :: I_GP_210Pb_mode    = 0             ! mode
  ! LAGRANGE
  LOGICAL :: L_LG               = .FALSE. ! LAGRANGE CALCULATION
  INTEGER :: I_LG_emis_method   = 1       ! emission method
  LOGICAL :: L_LG_emis_mcons    = .TRUE.  ! emission mass conserving ?
  INTEGER :: I_LG_emis_rest     = -1      ! how to handle emission 'rest'
  LOGICAL :: L_LG_emis_rest_int = .TRUE.  ! integrate 'rest'-flux
  ! ... chain calculation
  LOGICAL :: L_LG_chain         = .FALSE. ! 222Rn -> ... -> 210Pb
  CHARACTER(LEN=STRLEN_MEDIUM) :: C_LG_210Pb_aermod  = ''  ! aerosol model
  INTEGER :: I_LG_210Pb_mode    = 0             ! mode

  ! EXTERNAL FLUX
  TYPE(t_chaobj_cpl)                  :: C_Rn_flux ! CPL-namelist
  REAL(DP), DIMENSION(:,:),   POINTER :: ext_flux => NULL()

  ! POINTERS TO EXTERNAL OBJECTS
  REAL(DP), DIMENSION(:,:),   POINTER :: cvs => NULL()  ! snow cover fraction
  REAL(DP), DIMENSION(:,:),   POINTER :: slf => NULL()  ! sea-land fraction
  REAL(DP), DIMENSION(:,:,:), POINTER :: pint => NULL() ! pressure at interfaces

#ifdef MESSYTENDENCY
  INTEGER :: my_handle
#endif
  ! -------------------------------------------------------------

  ! PUBLIC ECHAM-5 INTERFACE ROUTINES TO 'messy_SUBMODELS'
  PUBLIC :: dradon_initialize       ! global initialisation of module
  PUBLIC :: dradon_new_tracer       ! define new tracers
  PUBLIC :: dradon_init_memory      ! allocate memory and define channels
  PUBLIC :: dradon_init_coupling    ! coupling to external fields
  PUBLIC :: dradon_global_start     ! read offline fields on event
  PUBLIC :: dradon_vdiff                ! integrate one time step
  !         -> dradon_vdiff_gp          ! GP
  !            -> dradon_vdiff_gp_tend  ! GP 1: tendency of lowest box)
  !                                     !       source + sink (decay)
  !            -> dradon_vdiff_gp_flux  ! GP 2: lower flux boundary cond.)
  !                                     !       source
  PUBLIC :: dradon_physc                ! integrate one time step
  !         -> dradon_physc_gp          ! GP
  !            -> dradon_physc_gp_decay ! GP 2: sink (decay)
  PUBLIC :: dradon_global_end
  !         -> dradon_global_end_lg     ! LG
  !                                     ! LG 1...4: source and sink

  ! PRIVATE ECHAM-5 INTERFACE ROUTINES
  ! PRIVATE :: dradon_read_nml_cpl    ! initialize 'coupling' to E5/ATTILA
                                     ! ( /CPL/-namelist )

  ! TRACER INDICES
  INTEGER, DIMENSION(NI) :: idx
#ifdef ECHAM5
  INTEGER, DIMENSION(NI) :: idx_lg
#endif

CONTAINS

! ************************************************************************
! PUBLIC ECHAM-5 INTERFACE ROUTINES
! ************************************************************************

! ------------------------------------------------------------------------
  SUBROUTINE  dradon_initialize

    ! Radon MODULE ROUTINE (ECHAM-5 INTERFACE)
    !
    ! INITIALIZATION OF GLOBAL VARIABLES FROM NAMELIST
    ! IN PARALLEL ENVIRONMENT
    ! INITIALIZATION OF Radon SPECIFIC EVENTS FROM NAMELIST
    ! IN PARALLEL ENVIRONMENT
    !
    ! Author: Patrick Joeckel, MPICH, July 2002

    ! MESSy/BMIL
    USE messy_main_mpi_bi,       ONLY: p_parallel_io, p_io, p_bcast
    USE messy_main_blather_bi,   ONLY: error_bi
    USE messy_main_tools,        ONLY: find_next_free_unit

    IMPLICIT NONE

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'dradon_initialize'
    INTEGER         :: iou    ! I/O unit
    INTEGER         :: status ! error status

    ! INITIALIZE MAIN-CTRL
    IF (p_parallel_io) THEN
       iou = find_next_free_unit(100,200)
       CALL dradon_read_nml_ctrl(status, iou)
       IF (status /= 0) CALL error_bi(' ',substr)
    END IF
    CALL p_bcast(ldradon, p_io)
    IF (.NOT.ldradon) RETURN
    !
    CALL p_bcast(I_Rn_flux_method, p_io)
    CALL p_bcast(R_Rn_cflux_land, p_io)
    CALL p_bcast(R_Rn_cflux_ocean, p_io)

    ! INITIALIZE CPL
    IF (p_parallel_io) THEN
       iou = find_next_free_unit(100,200)
       CALL dradon_read_nml_cpl(status, iou)
       IF (status /= 0) CALL error_bi(' ',substr)
    END IF
    CALL p_bcast(ldradon, p_io)
    IF (.NOT.ldradon) RETURN
    !
    !
    CALL p_bcast(L_GP, p_io)
    CALL p_bcast(I_GP_emis_method,   p_io)
    CALL p_bcast(L_GP_chain, p_io)
    CALL p_bcast(C_GP_210Pb_aermod, p_io) 
    CALL p_bcast(I_GP_210Pb_mode, p_io) 
    !
    CALL p_bcast(L_LG, p_io)
    CALL p_bcast(I_LG_emis_method,   p_io)
    CALL p_bcast(L_LG_emis_mcons,    p_io)
    CALL p_bcast(I_LG_emis_rest,     p_io)
    CALL p_bcast(L_LG_emis_rest_int, p_io)
    CALL p_bcast(L_LG_chain, p_io)
    CALL p_bcast(C_LG_210Pb_aermod, p_io) 
    CALL p_bcast(I_LG_210Pb_mode, p_io) 

    CALL p_bcast(C_Rn_flux%cha, p_io)
    CALL p_bcast(C_Rn_flux%obj, p_io)

    ! PRE-CALCULATION FOR CHAIN INTEGRATION
    IF (L_GP_chain .OR. L_LG_chain) CALL init_radon_chain

  END SUBROUTINE dradon_initialize
! ------------------------------------------------------------------------

! ------------------------------------------------------------------------
  SUBROUTINE dradon_new_tracer

    ! Radon MODULE ROUTINE (ECHAM-5 INTERFACE, PRIVATE)
    !
    ! defines Radon specific tracers
    !
    ! Author: Patrick Joeckel, MPICH, July 2003

    ! MESSy/BMIL
    USE messy_main_mpi_bi,        ONLY: p_parallel_io
    USE messy_main_tracer_mem_bi, ONLY: GPTRSTR
#ifdef ECHAM5
    USE messy_main_tracer_mem_bi, ONLY: LGTRSTR
#endif
    ! MESSy
    USE messy_main_tracer,        ONLY: new_tracer, set_tracer &
                                      , AEROSOL, ON, R_molarmass, I_scav &
                                      , I_drydep, I_sedi, S_aerosol_model &
                                      , I_aerosol_mode &
#ifdef ECHAM5
                                      , I_mix             &
#endif
                                      , R_aerosol_density

    IMPLICIT NONE
    INTRINSIC :: TRIM

    ! LOCAL
    INTEGER :: status
    CHARACTER(LEN=*), PARAMETER :: substr = 'dradon_new_tracer'

    IF (.NOT.ldradon) RETURN

    CALL start_message_bi(modstr, 'TRACER REQUEST', substr)

    idx(:)    = 0
#ifdef ECHAM5
    idx_lg(:) = 0
#endif

#ifdef MESSYTENDENCY
    my_handle = mtend_get_handle(modstr)
#endif
    ! REQUEST TRACER FIELDS
    IF (L_GP) THEN
       IF (p_parallel_io) THEN
          WRITE(*,*) ' ... Rn222 (GRIDPOINT)'
       END IF
       CALL new_tracer(status, GPTRSTR, 'Rn222', modstr      &
            , idx=idx(1) , unit='mol/mol')
       CALL set_tracer(status, GPTRSTR, idx(1), R_molarmass, r=222._dp)

       ! DECAY CHAIN
       IF (L_GP_chain) THEN
          IF (p_parallel_io) WRITE(*,*) ' ... Po218 (GRIDPOINT)'
          CALL new_tracer(status, GPTRSTR, 'Po218', modstr      &
               , idx=idx(2) , unit='mol/mol')
          CALL set_tracer(status, GPTRSTR, idx(2), R_molarmass, r=218._dp)
          !
          IF (p_parallel_io) WRITE(*,*) ' ... Pb214 (GRIDPOINT)'
          CALL new_tracer(status, GPTRSTR, 'Pb214', modstr      &
               , idx=idx(3) , unit='mol/mol')
          CALL set_tracer(status, GPTRSTR, idx(3), R_molarmass, r=214._dp)
          !
          IF (p_parallel_io) WRITE(*,*) ' ... Bi214 (GRIDPOINT)'
          CALL new_tracer(status, GPTRSTR, 'Bi214', modstr      &
               , idx=idx(4) , unit='mol/mol')
          CALL set_tracer(status, GPTRSTR, idx(4), R_molarmass, r=214._dp)
          !
          IF (p_parallel_io) WRITE(*,*) ' ... Pb210 (GRIDPOINT)'
          CALL new_tracer(status, GPTRSTR, 'Pb210', modstr      &
               , idx=idx(5) , unit='mol/mol', medium=AEROSOL)
          CALL set_tracer(status, GPTRSTR, idx(5), R_molarmass, r=210._dp)
          CALL set_tracer(status, GPTRSTR, idx(5), I_scav,   i=ON)
          CALL set_tracer(status, GPTRSTR, idx(5), I_drydep, i=ON)
          CALL set_tracer(status, GPTRSTR, idx(5), I_sedi,   i=ON)
          CALL set_tracer(status, GPTRSTR, idx(5), S_aerosol_model &
               , s=TRIM(C_GP_210Pb_aermod))
          CALL set_tracer(status, GPTRSTR, idx(5), I_aerosol_mode &
               , i=I_GP_210Pb_mode)
          CALL set_tracer(status, GPTRSTR, idx(5), R_aerosol_density &
               , r=2000._dp) 
          !
       END IF
       !
#ifdef MESSYTENDENCY
       IF (L_GP_chain) THEN
          CALL mtend_register(my_handle, mtend_id_vec=idx)
       ELSE
          CALL mtend_register(my_handle, idx(1))
       END IF
#endif
    END IF

#ifdef ECHAM5
    IF (L_LG) THEN
       IF (p_parallel_io) THEN
          WRITE(*,*) ' ... Rn222 (LAGRANGE)'
       END IF
       CALL new_tracer(status, LGTRSTR, 'Rn222', modstr  &
            , unit='mol/mol', idx=idx_lg(1))
       CALL set_tracer(status, LGTRSTR, idx_lg(1), R_molarmass, r=222._dp)

       ! DECAY CHAIN
       IF (L_LG_chain) THEN
          IF (p_parallel_io) WRITE(*,*) ' ... Po218 (LAGRANGE)'
          CALL new_tracer(status, LGTRSTR, 'Po218', modstr      &
               , idx=idx_lg(2) , unit='mol/mol')
          CALL set_tracer(status, LGTRSTR, idx_lg(2), R_molarmass, r=218._dp)
          !
          IF (p_parallel_io) WRITE(*,*) ' ... Pb214 (LAGRANGE)'
          CALL new_tracer(status, LGTRSTR, 'Pb214', modstr      &
               , idx=idx_lg(3) , unit='mol/mol')
          CALL set_tracer(status, LGTRSTR, idx_lg(3), R_molarmass, r=214._dp)
          !
          IF (p_parallel_io) WRITE(*,*) ' ... Bi214 (LAGRANGE)'
          CALL new_tracer(status, LGTRSTR, 'Bi214', modstr      &
               , idx=idx_lg(4) , unit='mol/mol')
          CALL set_tracer(status, LGTRSTR, idx_lg(4), R_molarmass, r=214._dp)
          !
          IF (p_parallel_io) WRITE(*,*) ' ... Pb210 (LAGRANGE)'
          CALL new_tracer(status, LGTRSTR, 'Pb210', modstr      &
               , idx=idx_lg(5) , unit='mol/mol', medium=AEROSOL)
          CALL set_tracer(status, LGTRSTR, idx_lg(5), R_molarmass, r=210._dp)
          CALL set_tracer(status, LGTRSTR, idx_lg(5), I_scav,   i=ON)
          CALL set_tracer(status, LGTRSTR, idx_lg(5), I_drydep, i=ON)
          CALL set_tracer(status, LGTRSTR, idx_lg(5), I_sedi,   i=ON)
          CALL set_tracer(status, LGTRSTR, idx_lg(5), S_aerosol_model &
               ,   s=TRIM(C_LG_210Pb_aermod))
          CALL set_tracer(status, LGTRSTR, idx_lg(5), I_aerosol_mode &
               , i=I_LG_210Pb_mode)
          CALL set_tracer(status, LGTRSTR, idx_lg(5), I_mix, i=ON)
          CALL set_tracer(status, LGTRSTR, idx_lg(5), R_aerosol_density &
               , r=2000._dp) 
          !
       END IF

    END IF
#endif
    CALL end_message_bi(modstr, 'TRACER REQUEST', substr)

  END SUBROUTINE dradon_new_tracer
! ------------------------------------------------------------------------

! ------------------------------------------------------------------------
  SUBROUTINE dradon_init_memory

    ! Radon MODULE ROUTINE (BASEMODEL INTERFACE)
    !
    ! define Radon specific channel(s) and allocate memory for
    ! global fields
    !
    ! Author: Patrick Joeckel, MPICH, July 2003

    ! MESSy/BMIL
    USE messy_main_mpi_bi,           ONLY: p_parallel_io
    USE messy_main_channel_error_bi, ONLY: channel_halt
    USE messy_main_channel_bi,       ONLY: GP_2D_HORIZONTAL, SCALAR
#ifdef ECHAM5
    USE messy_main_channel_bi,       ONLY: LG_ATTILA
#endif
    ! MESSy
    USE messy_main_channel,          ONLY: new_channel, new_channel_object &
                                         , new_attribute

    IMPLICIT NONE

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'dradon_init_memory'
    INTEGER                     :: status
    INTEGER                     :: nspec, jt
    REAL(DP), POINTER           :: lambda

    IF (.NOT.ldradon) RETURN

    CALL start_message_bi(modstr, 'CHANNEL DEFINITION', substr)

    ! DEFINE NEW CHANNEL
    CALL new_channel(status, modstr)
    CALL channel_halt(substr, status)

    IF (L_GP_chain .OR. L_LG_chain) THEN
       nspec = NI
    ELSE
       nspec = 1
    END IF
    DO jt=1, nspec
       NULLIFY(lambda)
       CALL new_channel_object(status, modstr, 'lambda_'//TRIM(iname(jt)) &
            , p0 = lambda, reprid = SCALAR, lstatic=.TRUE.)
       CALL channel_halt(substr, status)
       lambda = lam(jt)
       CALL new_attribute(status, modstr, 'lambda_'//TRIM(iname(jt))   &
            , 'long_name', c='decay constant of '//TRIM(iname(jt)) )
       CALL channel_halt(substr, status)
       CALL new_attribute(status, modstr, 'lambda_'//TRIM(iname(jt))   &
            , 'units', c='1/s' )
       CALL channel_halt(substr, status)
    END DO

    ! OBJECTS USED FOR GP AND LG
    IF (p_parallel_io) WRITE(*,*) ' ... Rn_fd'
    CALL new_channel_object(status, modstr, 'Rn_fd' &
         , p2 = rnflux &
         , reprid = GP_2D_HORIZONTAL )
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'Rn_fd'   &
         , 'long_name', c='222Radon flux density' )
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'Rn_fd'   &
         , 'units', c='atoms/m2/s' )
    CALL channel_halt(substr, status)

    IF (p_parallel_io) WRITE(*,*) ' ... Rn_fte'
    CALL new_channel_object(status, modstr, 'Rn_fte' &
         , p2 = rnflux_vmrs &
         , reprid = GP_2D_HORIZONTAL )
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'Rn_fte'   &
         , 'long_name', c='222Radon source tendency (GP)' )
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'Rn_fte'   &
         , 'units', c='mol/mol/s' )
    CALL channel_halt(substr, status)

    IF (p_parallel_io) WRITE(*,*) ' ... dp_sfc'
    CALL new_channel_object(status, modstr, 'dp_sfc' &
         , p2 = dp_sfc &
         , reprid = GP_2D_HORIZONTAL )
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'dp_sfc'   &
         , 'long_name', c='pressure difference at surface layer' )
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'dp_sfc'   &
         , 'units', c='Pa' )
    CALL channel_halt(substr, status)
    
#ifdef ECHAM5
    IF (L_LG) THEN

       IF (p_parallel_io) WRITE(*,*) ' ... Rn_fte_lg'
       CALL new_channel_object(status, modstr, 'Rn_fte_lg' &
            , p1 = lg_rnflux_vmrs &
            , reprid = LG_ATTILA )
       CALL channel_halt(substr, status)
       CALL new_attribute(status, modstr, 'Rn_fte_lg'   &
            , 'long_name', c='222Radon source tendency (LG)' )
       CALL channel_halt(substr, status)
       CALL new_attribute(status, modstr, 'Rn_fte_lg'   &
            , 'units', c='mol/mol/s' )
       CALL channel_halt(substr, status)

       ! ADDITIONAL DIAGNOSTIC OUTPUT FOR LG
       IF (I_LG_emis_rest == 1) THEN
          IF (p_parallel_io) WRITE(*,*) ' ... Rn_fte_rest'
          CALL new_channel_object(status, modstr, 'Rn_fte_rest' &
               , p2 = rnflux_vmrs_rest &
               , reprid = GP_2D_HORIZONTAL )
          CALL channel_halt(substr, status)
          CALL new_attribute(status, modstr, 'Rn_fte_rest'   &
               , 'long_name', c='non-captured 222Radon source tendency' )
          CALL channel_halt(substr, status)
          CALL new_attribute(status, modstr, 'Rn_fte_rest'   &
               , 'units', c='kg * mol/mol/s' )
          CALL channel_halt(substr, status)
       ELSE
          NULLIFY(rnflux_vmrs_rest)
       END IF

    END IF
#endif
    CALL end_message_bi(modstr, 'CHANNEL DEFINITION', substr)

  END SUBROUTINE dradon_init_memory
! ------------------------------------------------------------------------

! ------------------------------------------------------------------------
  SUBROUTINE dradon_init_coupling

    USE messy_main_blather_bi,       ONLY: info_bi
    USE messy_main_data_bi,          ONLY: basemodel => modstr
    USE messy_main_channel_error_bi, ONLY: channel_halt
    USE messy_main_channel,          ONLY: get_channel_object

    IMPLICIT NONE

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'dradon_init_coupling'
    INTEGER :: status
    
    CALL info_bi('checking for snow cover ('//TRIM(basemodel)//', cvs)')
    CALL get_channel_object(status, TRIM(basemodel), 'cvs', p2=cvs)
    CALL channel_halt(substr, status)

    CALL info_bi('checking for sea-land fractrion ('//TRIM(basemodel)//', slf)')
    CALL get_channel_object(status, TRIM(basemodel), 'slf', p2=slf)
    CALL channel_halt(substr, status)

    CALL info_bi('checking for pressure at layer interfaces ('//&
         &TRIM(basemodel)//', pressi)')
    CALL get_channel_object(status, TRIM(basemodel), 'pressi', p3=pint)
    CALL channel_halt(substr, status)

    IF (I_Rn_flux_method == I_Rn_flux_offline) THEN
       CALL info_bi('checking for external flux ('//&
            &TRIM(C_Rn_flux%cha)//', '//TRIM(C_Rn_flux%obj)//')')
       CALL get_channel_object(status                  &
            , TRIM(C_Rn_flux%cha), TRIM(C_Rn_flux%obj) &
            , p2 = ext_flux)
       CALL channel_halt(substr, status)
    END IF

  END SUBROUTINE dradon_init_coupling
! ------------------------------------------------------------------------

! ------------------------------------------------------------------------
  SUBROUTINE dradon_global_start

    ! Radon MODULE ROUTINE (BASEMODEL INTERFACE, PUBLIC)
    !
    ! CALCULATION OF 222Rn/FLUX in atoms/(m2 s):
    ! 0) constant flux over land/ocean
    ! 1) INPUT OF OFFLINE FIELDS (Rn - flux)
    !    TRIGGERED BY SPECIAL EVENTs
    !    I/O in PARALLEL ENVIRONMENT, DATA TRANSFER TO CHANNEL OBJECTS
    ! 2) Rn-flux model
    !
    ! Author: Patrick Joeckel, MPICH, July 2003

    ! MESSy/BMIL
    USE messy_main_mpi_bi,       ONLY: p_parallel_io
    USE messy_main_blather_bi,   ONLY: error_bi

    IMPLICIT NONE

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'dradon_global_start'

    IF (.NOT.ldradon) RETURN

    ! INITIALIZE

    SELECT CASE(I_Rn_flux_method)
    CASE(I_Rn_flux_const)    ! CONSTANT ------------------------------------
       ! LAND: no flux when covered by snow
       rnflux(:,:) =  slf(:,:)        * (1._DP - cvs(:,:)) * R_Rn_cflux_land  &
                   +  (1._DP - slf(:,:)) * R_Rn_cflux_ocean

    CASE(I_Rn_flux_offline)  ! OFFLINE -------------------------------------

       ! no flux when covered by snow
       rnflux(:,:) = ext_flux(:,:) * (1._DP - cvs(:,:))

    CASE(I_Rn_flux_online)   ! ONLINE --------------------------------------

       IF (p_parallel_io) THEN
          WRITE(*,*) '******************************************************'
          WRITE(*,*) '*** ERROR: ONLINE Rn-flux MODEL NOT YET IMPLEMENTED !'
          WRITE(*,*) '******************************************************'
       END IF
       CALL error_bi('ONLINE Rn-flux MODEL NOT YET IMPLEMENTED !',substr)

    CASE DEFAULT

       IF (p_parallel_io) THEN
          WRITE(*,*) '******************************************************'
          WRITE(*,*) '*** ERROR: UNKNOWN Rn-flux CALCULATION METHOD !'
          WRITE(*,*) '******************************************************'
       END IF
       CALL error_bi('UNKNOWN Rn-flux CALCULATION METHOD !',substr)

    END SELECT

  END SUBROUTINE dradon_global_start
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
  SUBROUTINE dradon_vdiff

    ! MESSy/BMIL
    USE messy_main_grid_def_mem_bi, ONLY: jrow, kproma, nlev
    USE messy_main_timer,           ONLY: ztmst => time_step_len
    USE messy_main_constants_mem,   ONLY: N_A, g, M_air

    IMPLICIT NONE

    INTRINSIC :: ABS

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'dradon_vdiff'

    ! CHECK SWITCH
    IF (.NOT.ldradon) RETURN

    ! CONVERT FLUX FROM atoms/(m2 s) TO mol/mol/s
    dp_sfc(1:kproma, jrow) = ABS( pint(_RI_XYZ__(1:kproma,jrow,nlev+1)) &
         - pint(_RI_XYZ__(1:kproma,jrow,nlev)) )
    CALL flux2vmrs(rnflux_vmrs(1:kproma,jrow), rnflux(1:kproma,jrow), &
         dp_sfc(1:kproma, jrow), mwair=M_air, na=N_A, g=g)

    IF (L_GP) CALL dradon_vdiff_gp

  CONTAINS

    ! -----------------------------------------------------------------
    SUBROUTINE dradon_vdiff_gp

      ! MESSy/BMIL
      USE messy_main_blather_bi,      ONLY: error_bi

      IMPLICIT NONE

      ! LOCAL
      CHARACTER(LEN=*), PARAMETER :: substr = 'dradon_vdiff_gp'
      
      SELECT CASE(I_GP_emis_method)
      CASE(0)
         ! DO NOTHING
      CASE(1)
         CALL dradon_vdiff_gp_tend
      CASE(2)
#ifdef ECHAM5
         CALL dradon_vdiff_gp_flux
#else
         CALL error_bi('FLUX OPTION ONLY IMPLEMENTED FOR ECHAM5',substr)
#endif

      CASE DEFAULT
         CALL error_bi(' *** UNKNOWN GP EMISSION METHOD ***',substr)
      END SELECT

    END SUBROUTINE dradon_vdiff_gp
    ! -----------------------------------------------------------------

    ! -----------------------------------------------------------------
    SUBROUTINE dradon_vdiff_gp_tend

      ! MESSy/BMIL
#ifndef MESSYTENDENCY
      USE messy_main_tracer_mem_bi, ONLY: pxtte => qxtte, pxtm1 => qxtm1
#endif

      IMPLICIT NONE

      ! LOCAL
      REAL(DP), DIMENSION(:,:,:), ALLOCATABLE :: zrn   ! radon field
      REAL(DP), DIMENSION(:,:,:), ALLOCATABLE :: zrnte ! radon field tendency
      INTEGER :: nspec      ! number of species
      INTEGER :: n
      INTEGER :: jp         ! column counter (1...kproma)
      INTEGER :: jk         ! level counter
#ifdef MESSYTENDENCY
      REAL(DP), DIMENSION(:,:), ALLOCATABLE :: ztmp   ! radon flux (h,v)
#endif
      
      ! ALLOCATE MEMORY
      IF (L_GP_chain) THEN
         nspec = NI
      ELSE
         nspec = 1
      END IF
      !
      ALLOCATE(zrn(kproma, nlev, nspec))
      ALLOCATE(zrnte(kproma, nlev, nspec))

#ifndef MESSYTENDENCY
      DO n=1, nspec
         DO jk=1,nlev
            zrn(1:kproma,jk,n) =  &
                 pxtm1(_RI_X_ZN_(1:kproma,jk,idx(n)))  &
                 + pxtte(_RI_X_ZN_(1:kproma,jk,idx(n)))*ztmst
            zrnte(1:kproma,jk,n) = 0.0_DP
         ENDDO
      END DO
#else
      DO n=1, nspec
         CALL mtend_get_start_l(idx(n), v0=zrn(:,:,n))
         zrnte(:,:,n) = 0.0_DP
      END DO
#endif

      chain: IF (L_GP_chain) THEN

#ifndef MESSYTENDENCY
         ! UPDATE TENDENCY AT SURFACE LEVELS (OPERATOR SPLITTING)
         n=1
         pxtte(_RI_X_ZN_(:,nlev,idx(n)))  = &
              pxtte(_RI_X_ZN_(:,nlev,idx(n))) + rnflux_vmrs(:, jrow)
         ! RECALCULATE START VALUE
         DO jk = 1, nlev
            zrn(1:kproma,jk,1) =  &
                 pxtm1(_RI_X_ZN_(1:kproma,jk,idx(n))) &
                 + pxtte(_RI_X_ZN_(1:kproma,jk,idx(n)))*ztmst
         ENDDO
#else
         n=1
         ALLOCATE(ztmp(kproma,nlev))
         ztmp(:,:) = 0.0_dp
         ztmp(:,nlev) = rnflux_vmrs(1:kproma, jrow)
         CALL mtend_add_l(my_handle, idx(n), px=ztmp)
         ! RECALCULATE START VALUE
         CALL mtend_get_start_l(idx(n), v0=zrn(:,:,n))
         DEALLOCATE(ztmp)
#endif

         DO jk = 1, nlev
            DO jp=1, kproma
               CALL int_radon_chain(zrnte(jp,jk,:), zrn(jp,jk,:), ztmst)
            END DO
         END DO

      ELSE

         ! TENDENCY AT SURFACE LEVEL: DECAY + SOURCE
         DO jp=1, kproma
            CALL int_radon(zrnte(jp, nlev,1), zrn(jp, nlev,1), ztmst &
                 ,rnflux_vmrs(jp, jrow))
         END DO

         ! TENDENCY AT NON-SURFACE LEVELS: DECAY
         DO jk = 1, nlev-1
            DO jp=1, kproma
               CALL int_radon(zrnte(jp,jk,1), zrn(jp,jk,1), ztmst)
            END DO
         END DO

      ENDIF chain

      ! UPDATE TENDENCIES
#ifndef MESSYTENDENCY
      DO n=1, nspec
         pxtte(_RI_X_ZN_(1:kproma,:,idx(n))) =  &
              pxtte(_RI_X_ZN_(1:kproma,:,idx(n))) + zrnte(:,:,n)
      END DO
#else
      DO n=1, nspec
         CALL mtend_add_l(my_handle, idx(n), px=zrnte(:,:,n))
      END DO
#endif

      ! DEALLOCATE MEMORY
      DEALLOCATE(zrn)
      DEALLOCATE(zrnte)

    END SUBROUTINE dradon_vdiff_gp_tend
    ! -----------------------------------------------------------------

    ! -----------------------------------------------------------------
    SUBROUTINE dradon_vdiff_gp_flux

#ifdef ECHAM5
      ! ECHAM5/MESSy
      USE messy_main_data_bi,       ONLY: pxtems
      USE messy_main_constants_mem, ONLY: N_A, M_air

      IMPLICIT NONE

      ! LOCAL
      ! air molecules / kg air
      REAL(DP), PARAMETER                :: uconv = N_A * 1.0e3_DP/M_air
      REAL(DP), DIMENSION(:,:), POINTER  :: zxtems => NULL()

      ! INIT
      zxtems => pxtems(:,1,:,jrow)

      ! [rnflux] = atoms / (m^2 s)
      !
      ! required unit:   (mol(tracer)/mol(air)) * kg (air) m-2 s-1
      zxtems(1:kproma,idx(1)) = zxtems(1:kproma,idx(1)) + &
           rnflux(1:kproma,jrow) / uconv
#endif

    END SUBROUTINE dradon_vdiff_gp_flux
    ! -----------------------------------------------------------------

  END SUBROUTINE dradon_vdiff
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
  SUBROUTINE dradon_physc

    ! ECHAM5/MESSy
    USE messy_main_grid_def_mem_bi,ONLY: nlev, kproma
    USE messy_main_timer,          ONLY: ztmst => time_step_len

    IMPLICIT NONE

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'dradon_physc'

    ! CHECK SWITCH
    IF (.NOT.ldradon) RETURN

    IF (L_GP) CALL dradon_physc_gp

    ! NOTHING FOR L_LG

  CONTAINS

    ! -----------------------------------------------------------------
    SUBROUTINE dradon_physc_gp

      ! ECHAM5/MESSy
      USE messy_main_blather_bi, ONLY: error_bi

      IMPLICIT NONE

      SELECT CASE(I_GP_emis_method)
      CASE(0)
         ! DO NOTHING
      CASE(1)
         ! DO NOTHING
      CASE(2)
         CALL dradon_physc_gp_decay
      CASE DEFAULT
         CALL error_bi(' *** UNKNOWN GP EMISSION METHOD ***',substr)
      END SELECT

    END SUBROUTINE dradon_physc_gp
    ! -----------------------------------------------------------------

    ! -----------------------------------------------------------------
    SUBROUTINE dradon_physc_gp_decay

      ! MESSy/BMIL
#ifndef MESSYTENDENCY
      USE messy_main_tracer_mem_bi,  ONLY: pxtte => qxtte, pxtm1 => qxtm1
#endif

      IMPLICIT NONE

      ! LOCAL
      REAL(DP), DIMENSION(:,:,:), ALLOCATABLE :: zrn   ! radon field
      REAL(DP), DIMENSION(:,:,:), ALLOCATABLE :: zrnte ! radon field tendency
      INTEGER :: nspec, n
      INTEGER :: jp         ! column counter (1...kproma)
      INTEGER :: jk         ! level counter

      ! ALLOCATE MEMORY
      IF (L_GP_chain) THEN
         nspec = NI
      ELSE
         nspec = 1
      END IF
      !
      ALLOCATE(zrn(kproma, nlev, nspec))
      ALLOCATE(zrnte(kproma, nlev, nspec))

#ifndef MESSYTENDENCY
      DO n=1, nspec
         zrn(:,:,n) = pxtm1(_RI_X_ZN_(1:kproma,:,idx(n))) &
              + pxtte(_RI_X_ZN_(1:kproma,:,idx(n)))*ztmst
         zrnte(:,:,n) = 0.0_DP
      END DO
#else
      DO n=1, nspec
         CALL mtend_get_start_l(idx(n), v0=zrn(:,:,n))
         zrnte(:,:,n) = 0.0_DP
      END DO
#endif

      chain: IF (L_GP_chain) THEN

         ! TENDENCY AT ALL LEVELS: DECAY
         DO jk = 1, nlev
            DO jp=1, kproma
               CALL int_radon_chain(zrnte(jp,jk,:), zrn(jp,jk,:), ztmst)
            END DO
         END DO

      ELSE

         ! TENDENCY AT ALL LEVELS: DECAY
         DO jk = 1, nlev
            DO jp=1, kproma
               CALL int_radon(zrnte(jp,jk,1), zrn(jp,jk,1), ztmst)
            END DO
         END DO

      ENDIF chain

      ! UPDATE TENDENCIES
#ifndef MESSYTENDENCY
      DO n=1, nspec
         pxtte(_RI_X_ZN_(1:kproma,:,idx(n))) = &
              pxtte(_RI_X_ZN_(1:kproma,:,idx(n))) + zrnte(:,:,n)
      END DO
#else
      DO n=1, nspec
         CALL mtend_add_l(my_handle, idx(n), px=zrnte(:,:,n))
      END DO
#endif

      ! DEALLOCATE MEMORY
      DEALLOCATE(zrn)
      DEALLOCATE(zrnte)

    END SUBROUTINE dradon_physc_gp_decay
    ! -----------------------------------------------------------------

  END SUBROUTINE dradon_physc
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
  SUBROUTINE dradon_global_end
!!#D attila +
    IMPLICIT NONE

    ! CHECK SWITCH
    IF (.NOT.ldradon) RETURN

    IF (L_LG) CALL dradon_global_end_lg

  CONTAINS

    ! -----------------------------------------------------------------
    SUBROUTINE dradon_global_end_lg

#ifdef ECHAM5
      ! ECHAM5/MESSy
      USE messy_main_tracer_mem_bi,  ONLY: pxtte_a => qxtte_a &
                                         , pxtm1_a => qxtm1_a &
                                         , NCELL
      USE messy_attila_tools_e5,      ONLY: gpsfemis2lgemis_e5
      USE messy_main_grid_def_mem_bi, ONLY: nproma, npromz, ngpblks 
      USE messy_main_timer,           ONLY: ztmst => time_step_len

      IMPLICIT NONE

      ! LOCAL
      REAL(DP), DIMENSION(:,:), ALLOCATABLE :: zrn_lg     ! local radon field
      REAL(DP), DIMENSION(:,:), ALLOCATABLE :: zrnte_lg   ! tendency
      INTEGER   :: jn ! CELL counter (1...NCELL)
      ! SPECIAL FOR GRIDPOINT REST
      REAL(DP), DIMENSION(:,:), ALLOCATABLE :: zrn_sf   ! radon field
      REAL(DP), DIMENSION(:,:), ALLOCATABLE :: zrnte_sf ! radon field tendency
      INTEGER   :: jp, zjrow, zkproma, n
      INTEGER   :: nspec

      ! ALLOCATE MEMORY
      IF (L_LG_chain) THEN
         nspec = NI
      ELSE
         nspec = 1
      END IF
      !
      ALLOCATE(zrn_lg(NCELL,nspec))
      ALLOCATE(zrnte_lg(NCELL,nspec))

      ! (MASS CONSERVING) TRANSFORMATION FROM GP-SPACE TO LAGRANGIAN-SPACE
      ! - INCLUDING DECOMPOSITION (CHANNEL OBJECT -> CHANNEL OBJECT)
      IF (I_LG_emis_rest == 1) THEN
         CALL gpsfemis2lgemis_e5(rnflux_vmrs, lg_rnflux_vmrs    &
              , I_LG_emis_method                                &
              , gprl=rnflux_vmrs_rest, lmcons = L_LG_emis_mcons)
      ELSE
         CALL gpsfemis2lgemis_e5(rnflux_vmrs, lg_rnflux_vmrs    &
              , I_LG_emis_method                                &
              , lmcons = L_LG_emis_mcons)
      END IF
      
      ! UPDATE START VALUES ; RESET LOCAL TENDENCIES
      DO n=1, nspec
         zrn_lg(:,n)   = pxtm1_a(:,idx_lg(n)) + pxtte_a(:,idx_lg(n))*ztmst
         zrnte_lg(:,n) = 0.0_DP
      END DO

      chain: IF (L_LG_chain) THEN
      
         ! UPDATE TENDENCY (EMISSION AT SURFACE LEVELS ; OPERATOR SPLITTING)
         pxtte_a(:,idx_lg(1))  = pxtte_a(:,idx_lg(1)) + lg_rnflux_vmrs(:)
         ! RECALCULATE START VALUE
         zrn_lg(:,1) = pxtm1_a(:,idx_lg(1)) + pxtte_a(:,idx_lg(1))*ztmst

         DO jn=1, NCELL
            CALL int_radon_chain(zrnte_lg(jn,:), zrn_lg(jn,:), ztmst)
         END DO

      ELSE

         DO jn=1, NCELL
            CALL int_radon(zrnte_lg(jn,1), zrn_lg(jn,1), ztmst &
                 ,lg_rnflux_vmrs(jn) )
         END DO
         
      END IF chain

      ! UPDATE TENDENCIES
      DO n=1, nspec
         pxtte_a(:,idx_lg(n)) = pxtte_a(:,idx_lg(n)) + zrnte_lg(:,n)
      END DO

      ! CLEAN MEMORY
      DEALLOCATE(zrn_lg)
      DEALLOCATE(zrnte_lg)

      IF (L_LG_emis_rest_int) THEN
         ! REST (MASS) IN GRIDPOINT SPACE NEEDS ALSO TO BE INTEGRATED:
         ! THE ACCUMULATED 'REST-(MASS)-FLUX' WHERE NO LAGRANGIAN CELLS
         ! WERE LOCATED, DECAYS WITH TIME ...
         ! [rnflux_vmrs_rest] = kg * (mol/mol/s)
         !
         ALLOCATE(zrnte_sf(nproma, ngpblks))
         zrnte_sf(:,:) = 0.0                  ! INITIALIZE TENDENCY
         ALLOCATE(zrn_sf(nproma, ngpblks))
         zrn_sf(:,:) = 0.0                    ! INITIAL
         !
         DO zjrow = 1, ngpblks
            IF ( zjrow == ngpblks ) THEN
               zkproma = npromz
            ELSE
               zkproma = nproma
            END IF
            DO jp=1, zkproma   ! INTEGRATE IN GRID POINT SPACE
                               ! (SOURCE + DECAY FROM ZERO INITIAL)
                               ! ONE TIME STEP
               CALL int_radon(zrnte_sf(jp, zjrow), zrn_sf(jp, zjrow)  &
                    , ztmst, rnflux_vmrs_rest(jp, zjrow))
            END DO
         END DO
         ! UPDATE FLUX WITH INTEGRATED TENDENCY
         rnflux_vmrs_rest(:,:) = zrnte_sf(:,:)
         !
         DEALLOCATE(zrn_sf)
         DEALLOCATE(zrnte_sf)
      END IF
#endif
    END SUBROUTINE dradon_global_end_lg
    ! -----------------------------------------------------------------
!!#D attila -
  END SUBROUTINE dradon_global_end
! ----------------------------------------------------------------------

! ************************************************************************
! PRIVATE ECHAM-5 INTERFACE ROUTINES
! ************************************************************************

! ----------------------------------------------------------------------
  SUBROUTINE dradon_read_nml_cpl(status, iou)

    ! DRADON MODULE ROUTINE (ECHAM-5 INTERFACE, PRIVATE)
    !
    ! read namelist for 'coupling' to ECHAM5/ATTILA
    !
    ! Author: Patrick Joeckel, MPICH, Oct 2003

    ! MESSy
    USE messy_main_tools,         ONLY: read_nml_open, read_nml_check &
                                      , read_nml_close
#ifdef ECHAM5
    ! MESSy/BMIL
    USE messy_main_tracer_mem_bi, ONLY: NGCELL
#endif

    IMPLICIT NONE

    ! I/O
    INTEGER, INTENT(OUT) :: status     ! error status
    INTEGER, INTENT(IN)  :: iou        ! I/O unit

    ! (LOCAL) NAMELIST VARIABLES
    CHARACTER(LEN=*), PARAMETER :: substr = 'dradon_read_nml_cpl'

    NAMELIST /CPL/ L_GP, I_GP_emis_method      &
         , L_GP_chain, C_GP_210Pb_aermod, I_GP_210Pb_mode &
         , L_LG                               &
         , I_LG_emis_method, L_LG_emis_mcons  &
         , I_LG_emis_rest, L_LG_emis_rest_int &
         , L_LG_chain, C_LG_210Pb_aermod, I_LG_210Pb_mode &
         , C_Rn_flux

    ! LOCAL
    LOGICAL :: lex      ! file exists ?
    INTEGER :: fstat    ! file status

    status = 1

    CALL read_nml_open(lex, substr, iou, 'CPL', modstr)
    IF (.not.lex) RETURN    ! <modstr>.nml does not exist

    READ(iou, NML=CPL, IOSTAT=fstat)
    CALL read_nml_check(fstat, substr, iou, 'CPL', modstr)
    IF (fstat /= 0) RETURN  ! error while reading namelist

    ! DIAGNOSE NAMELIST AND SET GLOBAL SWITCHES
    !
    ! CHECK NAMELIST
    IF (L_GP) THEN
       WRITE(*,*) 'GRIDPOINT  INTEGRATION: ON'
       !
       SELECT CASE(I_GP_emis_method)
       CASE(1)
          WRITE(*,*) 'GP EMISSION METHOD    : SRUFACE LAYER TENDENCY'
       CASE(2)
          WRITE(*,*) 'GP EMISSION METHOD    : '&
               &//'LOWER BOUNDARY CONDITION OF FLUX'
       CASE DEFAULT
          ! ERROR: UNKNOWN EMISSION METHOD
          WRITE(*,*) '*** ERROR *** UNKNOWN EMISSION METHOD'
          RETURN ! ERROR
       END SELECT
       !
    ELSE
       WRITE(*,*) 'GRIDPOINT  INTEGRATION: OFF'
    END IF

    IF (L_GP_chain) THEN
       WRITE(*,*) 'GP DECAY CHAIN        : ON'
    ELSE
       WRITE(*,*) 'GP DECAY CHAIN        : OFF'
    END IF

    IF (L_LG) THEN
       WRITE(*,*) 'LAGRANGIAN INTEGRATION: ON'
    ELSE
       WRITE(*,*) 'LAGRANGIAN INTEGRATION: OFF'
    END IF

    IF ((.NOT.L_GP).AND.(.NOT.L_LG)) THEN
       WRITE(*,*) 'WARNING: GRIDPOINT AND LAGRANGIAN INTEGRATION'
       WRITE(*,*) '         BOTH SWITCHED OFF!'
       WRITE(*,*) ' -> SWITCHING OFF DRADON ENTIRELY'
       ldradon = .FALSE.
    END IF

#ifdef ECHAM5
    IF (L_LG) THEN
       IF (NGCELL <= 0) THEN
          WRITE(*,*) 'WARNING: LAGRANGIAN INTEGRATION REQUIRES'
          WRITE(*,*) '         ATTILA!'
          WRITE(*,*) ' -> SWITCHING OFF LAGRANGIAN INTEGRATION FOR DRADON'
          L_LG = .FALSE.
       ELSE
          !
          SELECT CASE(I_LG_emis_method)
          CASE(1)
             WRITE(*,*) 'LG EMISSION METHOD    : SRUFACE LAYER ONLY'
          CASE(2)
             WRITE(*,*) 'LG EMISSION METHOD    : '&
                  &//'LOWEST CELL(S) IN BOUNDARY LAYER'
          CASE(3)
             WRITE(*,*) 'LG EMISSION METHOD    : ALL CELL(S) IN BOUNDARY LAYER'
             WRITE(*,*) '                        WEIGHTED'
          CASE(4)
             WRITE(*,*) 'LG EMISSION METHOD    : ALL CELL(S) IN BOUNDARY LAYER'
             WRITE(*,*) '                        UNWEIGHTED'
          CASE DEFAULT
             ! ERROR: UNKNOWN EMISSION METHOD
             WRITE(*,*) '*** ERROR *** UNKNOWN EMISSION METHOD'
             RETURN ! ERROR
          END SELECT
          !
          IF (L_LG_emis_mcons) THEN
             WRITE(*,*) 'LG EMISSION           : MASS CONSERVING ( = GP)'
          ELSE
             WRITE(*,*) 'LG EMISSION           : NON-MASS CONSERVING ( != GP)'
          END IF
          !
          SELECT CASE(I_LG_emis_rest)
          CASE(-1)
             WRITE(*,*) 'LG EMISSION REST      : AUTOMATIC ...'
             SELECT CASE(I_LG_emis_method)
             CASE(1)
                WRITE(*,*) '                  ... : ON'
                I_LG_emis_rest = 1
             CASE(2,3,4)
                WRITE(*,*) '                  ... : OFF'
                I_LG_emis_rest = 0
             END SELECT
          CASE(0)
             WRITE(*,*) 'LG EMISSION REST      : OFF'
          CASE(1)
             WRITE(*,*) 'LG EMISSION REST      : ON'
          CASE DEFAULT
             ! ERROR: UNKNOWN EMISSION REST HANDLING
             WRITE(*,*) '*** ERROR *** UNKNOWN EMISSION REST HANDLING'
             RETURN ! ERROR
          END SELECT
          !
          IF (L_LG_emis_mcons .AND. (I_LG_emis_rest == 1)) THEN
             IF (L_LG_emis_rest_int) THEN
                WRITE(*,*) 'LG EMISSION REST INT. : ON'
             ELSE
                WRITE(*,*) 'LG EMISSION REST INT. : OFF'
             END IF
          ELSE
             WRITE(*,*) 'LG EMISSION REST INT. : OBSOLETE'
             WRITE(*,*) &
                  '                        WITHOUT L_LG_emis_mcons = .true.'  
             WRITE(*,*) &
                  '                AND/OR  WITHOUT I_LG_emis_rest  = 1' 
             L_LG_emis_rest_int = .false.
          END IF
          !
          IF (L_LG_chain) THEN
             WRITE(*,*) 'LG DECAY CHAIN        : ON'
          ELSE
             WRITE(*,*) 'LG DECAY CHAIN        : OFF'
          END IF
          !
       END IF
    ELSE
       WRITE(*,*) ' *_LG_* SWITCHES OBSOLETE WITHOUT L_LG = T !'
    END IF
#endif    

#ifdef ECHAM5
    ! AVOID ERRORS
    IF (L_LG_emis_rest_int .AND. L_LG_chain) THEN
       WRITE(*,*) 'ERROR: L_LG_emis_rest_int .AND. L_LG_chain'
       WRITE(*,*) '       is in conflict with mass-conservation'
       status = 1
       RETURN
    END IF
#endif
    CALL read_nml_close(substr, iou, modstr)

    status = 0  ! no ERROR

  END SUBROUTINE dradon_read_nml_cpl
! ----------------------------------------------------------------------

! **********************************************************************
END MODULE messy_dradon_si
! **********************************************************************

MODULE messy_h2o_e5

  !   H2O CORRECTIVE MODULE FOR ECHAM5 (MESSy-INTERFACE)
  !
  !   THIS MODULE
  !      - DEFINES A TRACER H2O
  !      - INITIALIZES THE H2O TRACER CORRECTLY (STRATOSPHERE, I.E
  !        'ABOVE' 100hPa) FROM SATELLITE DATA
  !        (SEE h2o_t.nml)
  !      - CHANGES THE SPECIFIC HUMIDITY IN CONSISTENCE WITH THE
  !        H2O TRACER
  !      - CALCULATES A CLIMATOLOGICAL CHEMICAL TENDENCY OF H2O, DUE TO
  !        REACTIONS CH4 + OH,O1D,Cl
  !      - DECIDES, WHICH CHEMICAL WATER VAPOUR TENDECIES ARE FEEDBACK
  !        INTO THE SPECIFIC HUMIDITY
  !        (SEE NAMELIST 'CPL' IN h2o.nml)
  !
  !   Authors: 
  !   Christoph Bruehl, Oktober 2003:
  !     - original code
  !
  !   Patrick Joeckel, November 2003:
  !     - restructured according to MESSy-standard
  !     - time filter consistency ('synchronize' (q,qte) and (H2O, H2Ote))
  !     - simplified
  !     - made independent of vertical resolution
  !
!!#D attila +
  !   Sabine Brinkop, December 2008:
  !     - inclusion of LG water variables/tendencies plus transformation to GP
!!#D attila -
  !
  !   TODO:
  !         - check for ECHAM5 time filter consistency !!!
  !         - forcheck analysis
  !         - testing
  !         - check if there is a difference between the results
  !           for I_H2O_TENDENCY = 1 and 2, even if no other SM changes the
  !           H2O tendency
  
  ! =========================================================================
  !            MODULE H2O  (C) Patrick Joeckel, MPICH, Nov 2003
  !   ----------------------------------------------------------------------
  !                                 H2O,te=0 <----- h2o_new_tracer
  !                                    |             (DEFINE TRACER)
  !                                    |
  !                                    V
!!$  !                                 H2O,te=0 <-----  h2o_tracer_init
!!$  !                                    |           (INITIALIZE TRACER)
!!$  !                                    |         (STRAT.: SATELLITE DATA)
!!$  !                                    |         (TROP.:  q)
  !                                    |
  !   --------------------------->  H2O,te   <------ h2o_global_start    
  !   ^   IF (lstart)  q   <- H2O      |        (UPDATE OFFLINE FILEDS ...
  !   |   ELSE         q   -> H2O      |         ... FOR CLIM. CHEM. TEND.)
  !   |   ALWAYS       qte -> H2Ote    |        (TRANSFER (q, qte) 
  !   |                                |             <--> (H2O, H2Ote)
  !  =|====ECHAM5 CORE =======         |
  !  ||             ^        |         |
  !  ||       I     |        |         |             h2o_physc
  !  ||   C   N     |        |         |      (CALC CLIM. CHEM. TEND.)        
  !  || A O V T     |        |         |               |  
  !  || D N D E     |        |         |               |  
  !  || V V I G     |        |         |               |    other SUBMODELS
  !  || E E F R     |        |         |               |    (add chem. tend.
  !  || C C F A     |        |         |               |      to H2Ote)
  !  || T T   T     |        |         |               |              |
  !  ||       E     |        |         |               ------|        |
  !  || + + + + [q, qte]     |         |                     |        |
  !  || - - - + [H2O, H2OTE] |         |                     |        |
  !  =|=======================         |                     |        |
  !   |                                |                     |        |
  !   |                                V                     |        |
  !   |<--------------------------  H2O,te <-- h2o_local_end |        |
  !            qte <- qte + H2Ote             -2: H2Ote = 0  |        |
  !          H2Ote <- qte                                    |        |
  !                                           -1: H2Ote = 0  |        |
  !                                                          |        |
  !                                            0: H2Ote = <-----------|
  !                                                          |        |
  !                                            1: H2Ote = <--- [+] <--|
  !                                                          |
  !                                            2: H2Ote = <--|
  ! =========================================================================

  ! ECHAM5/MESSy
  USE messy_main_blather_bi,             ONLY: start_message_bi &
                                             , end_message_bi   &
                                             , error_bi
  ! MESSy
  USE messy_main_constants_mem,  ONLY: M_air, M_H2O
  USE messy_h2o
#ifdef MESSYTENDENCY
  USE messy_main_tendency_bi,    ONLY: mtend_get_handle,      &
                                       mtend_add_g,           &
                                       mtend_register,        &
                                       mtend_id_q,            &
                                       mtend_id_tracer
#endif

  IMPLICIT NONE
  PRIVATE

  ! PUBLIC INTERFACE ROUTINES
  PUBLIC :: h2o_initialize
  PUBLIC :: h2o_init_coupling  ! mim_sb_20091009
  PUBLIC :: h2o_new_tracer     ! define tracers
  PUBLIC :: h2o_init_memory
  PUBLIC :: h2o_init_tracer    ! initialize tracers
  PUBLIC :: h2o_global_start
  PUBLIC :: h2o_global_end     ! LG, feedback
  PUBLIC :: h2o_physc
  PUBLIC :: h2o_local_end
  PUBLIC :: h2o_free_memory
 
  ! PRIVATE HELPER ROUTINES
  !PRIVATE :: h2o_read_nml_cpl

  ! GLOABL CPL-NAMELIST PARAMETERS
  ! DEFAULT: USE EXCLUSIVELY TEND. OF THIS MODULE
  LOGICAL, PUBLIC :: L_GP = .TRUE.
  LOGICAL, PUBLIC :: L_LG = .FALSE.  ! mz_pj_20080201
  INTEGER, PUBLIC :: I_FEEDBACK = 1  ! 0=no sync between LG and GP, 1=GP, 2=LG
  INTEGER, PUBLIC :: I_H2O_TENDENCY = 2
  LOGICAL, PUBLIC :: L_MODINI_QZAV = .FALSE.

  ! GLOBAL PARAMETERS
  REAL(DP), PARAMETER :: scvmr = M_air/M_H2O
  REAL(DP), PARAMETER :: pmerge = 10000.0_DP       ! 100 hPa [Pa]

  ! Lagrangian moist convection
  LOGICAL, PUBLIC :: L_LGMC = .FALSE.    ! mim_sb_20091009
  ! GLOBAL VARIABLES / PARAMETERS / POINTERS
  INTEGER                         :: idt_h2o_gp
  INTEGER                         :: idt_h2o_lg     ! mim_sb_20080201
  INTEGER                         :: idt_h2o_liq_lg ! mim_sb_20091009
  INTEGER                         :: idt_h2o_ice_lg ! mim_sb_20091009
  ! DIAGNOSTIC CHANNEL OBJECTS
  REAL(DP), DIMENSION(:,:,:), POINTER :: ch4s_3d => NULL()
  REAL(DP), DIMENSION(:,:,:), POINTER :: CH4_H   => NULL()
  REAL(DP), DIMENSION(:,:,:), POINTER :: CH4OH   => NULL()
  REAL(DP), DIMENSION(:,:,:), POINTER :: CH4O1D  => NULL()
  REAL(DP), DIMENSION(:,:,:), POINTER :: CH4CL   => NULL()
  ! for LG
  REAL(DP), DIMENSION(:),     POINTER :: ch4s_1d => NULL()
  REAL(DP), DIMENSION(:,:,:), POINTER :: qte_adv => NULL()

  REAL(DP),                   POINTER :: rLG_H2O => NULL()
  ! Channel OBJECTS
  REAL(DP), DIMENSION(:,:,:), POINTER :: conv_qte  => NULL()
#ifdef MESSYTENDENCY
  ! variable for tendency budget
  INTEGER                                     :: my_handle
  REAL(kind=dp),DIMENSION(:,:,:),ALLOCATABLE  :: lo_h2ote
#endif

  INTRINSIC MAX, SIZE

CONTAINS

!************************************************************************
  SUBROUTINE h2o_initialize

    ! ECHAM5/MESSy
    USE messy_main_mpi_bi,               ONLY: p_parallel_io, p_io, p_bcast
    USE messy_main_tools,                ONLY: find_next_free_unit 

    IMPLICIT NONE

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'h2o_initialize'
    INTEGER                     :: iou    ! I/O unit
    INTEGER                     :: status ! error status

    ! INITIALIZE CPL-NAMELIST
    IF (p_parallel_io) THEN
       iou = find_next_free_unit(100,200)
       CALL h2o_read_nml_cpl(status, iou)
       IF (status /= 0) CALL error_bi(' ', substr)
    END IF
    CALL p_bcast(L_GP, p_io)
    CALL p_bcast(L_LG, p_io) ! mz_pj_20080201
    CALL p_bcast(I_FEEDBACK, p_io)
    CALL p_bcast(I_H2O_TENDENCY, p_io)
    CALL p_bcast(L_MODINI_QZAV, p_io)

  END SUBROUTINE h2o_initialize
!************************************************************************

!************************************************************************
  SUBROUTINE h2o_init_coupling

    ! ECHAM5/MESSy
    USE messy_main_channel_error_bi,       ONLY: channel_halt

    ! MESSy
    USE messy_main_channel,          ONLY: get_channel_object, get_channel_info

    IMPLICIT NONE

    CHARACTER(LEN=*), PARAMETER :: substr='h2o_init_coupling'
    INTEGER                     :: status

    IF (I_H2O_TENDENCY > 0) THEN
       CALL get_channel_object(status, 'import_grid', 'A_CH4_H', p3=CH4_H)
       CALL channel_halt(substr, status)
       
       CALL get_channel_object(status, 'import_grid', 'B_CH4OH', p3=CH4OH)
       CALL channel_halt(substr, status)
       
       CALL get_channel_object(status, 'import_grid', 'C_CH4O1D', p3=CH4O1D)
       CALL channel_halt(substr, status)
       
       CALL get_channel_object(status, 'import_grid', 'D_CH4CL', p3=CH4Cl)
       CALL channel_halt(substr, status)
    END IF

!!#D attila +
    IF (.NOT. L_LG) RETURN

    CALL start_message_bi(modstr,'COUPLING INITIALIZATION',substr)
    status = 0

    CALL get_channel_info(status, 'lgmc')
    L_LGMC = (status == 0)

    IF (L_LGMC .and. I_FEEDBACK == 2) THEN

       CALL get_channel_object(status, 'convect', 'conv_qte', p3=conv_qte)
       CALL channel_halt(substr, status)

     ENDIF

    CALL end_message_bi(modstr,'COUPLING INITIALIZATION',substr)    

!!#D attila -
  END SUBROUTINE h2o_init_coupling
!************************************************************************

!************************************************************************
  SUBROUTINE h2o_new_tracer
    
    ! ECHAM5/MESSy
    USE messy_main_tracer_mem_bi,   ONLY: GPTRSTR, LGTRSTR
    ! MESSy
    USE messy_main_tracer,          ONLY: new_tracer &
                                        , AIR, OFF, ON, AMOUNTFRACTION &
                                        , I_ADVECT, I_CONVECT, I_INTEGRATE &
                                        , I_VDIFF, R_molarmass, set_tracer
    USE messy_main_tracer_tools_bi, ONLY: tracer_halt
    USE messy_main_constants_mem,   ONLY: M_H2O

    IMPLICIT NONE

    ! LOCAL
    INTEGER :: status
    CHARACTER(LEN=*), PARAMETER :: substr = 'h2o_new_tracer'

    ! moved from initialise
#ifdef MESSYTENDENCY
     my_handle = mtend_get_handle(modstr)
     CALL mtend_register (my_handle,mtend_id_q)
     CALL mtend_register (my_handle,mtend_id_tracer)
#endif
    CALL start_message_bi(modstr,'TRACER REQUEST',substr)

    IF (L_GP) THEN
      CALL new_tracer(status, GPTRSTR, 'H2O', modstr     &
            ,unit = 'mol/mol'      &
            ,quantity = AMOUNTFRACTION  &
            ,medium = AIR  &
            ,idx = idt_H2O_gp)
      CALL tracer_halt(substr, status)

      CALL set_tracer(status, GPTRSTR, idt_H2O_gp, R_molarmass, M_H2O)
      CALL tracer_halt(substr, status)
      CALL set_tracer(status, GPTRSTR, idt_H2O_gp, I_ADVECT,    OFF)
      CALL tracer_halt(substr, status)
      CALL set_tracer(status, GPTRSTR, idt_H2O_gp, I_VDIFF,     OFF)
      CALL tracer_halt(substr, status)
      CALL set_tracer(status, GPTRSTR, idt_H2O_gp, I_CONVECT,   OFF)
      CALL tracer_halt(substr, status)
      CALL set_tracer(status, GPTRSTR, idt_H2O_gp, I_INTEGRATE, ON)
      CALL tracer_halt(substr, status)
    END IF

!!#D attila +
    IF (L_LG) THEN 
       CALL new_tracer(status, LGTRSTR, 'H2O', modstr     &
            ,unit = 'mol/mol'      &
            ,quantity = AMOUNTFRACTION  &
            ,medium = AIR  &
            ,idx = idt_H2O_lg)
      CALL tracer_halt(substr, status)

      CALL set_tracer(status, GPTRSTR, idt_H2O_lg, R_molarmass, M_H2O)
      CALL tracer_halt(substr, status)
      CALL set_tracer(status, GPTRSTR, idt_H2O_lg, I_ADVECT,    OFF)
      CALL tracer_halt(substr, status)
      CALL set_tracer(status, GPTRSTR, idt_H2O_lg, I_VDIFF,     OFF)
      CALL tracer_halt(substr, status)
      CALL set_tracer(status, GPTRSTR, idt_H2O_lg, I_CONVECT,   OFF)
      CALL tracer_halt(substr, status)
      CALL set_tracer(status, GPTRSTR, idt_H2O_lg, I_INTEGRATE, ON)
      CALL tracer_halt(substr, status)

      CALL new_tracer(status, LGTRSTR, 'LIQ', modstr     &
            ,unit = 'mol/mol'      &
            ,quantity = AMOUNTFRACTION  &
            ,medium = AIR &
            ,idx = idt_H2O_liq_lg)

      CALL tracer_halt(substr, status)

      CALL set_tracer(status, GPTRSTR, idt_H2O_liq_lg, R_molarmass, M_H2O)
      CALL tracer_halt(substr, status)
      CALL set_tracer(status, GPTRSTR, idt_H2O_liq_lg, I_ADVECT,    OFF)
      CALL tracer_halt(substr, status)
      CALL set_tracer(status, GPTRSTR, idt_H2O_liq_lg, I_VDIFF,     OFF)
      CALL tracer_halt(substr, status)
      CALL set_tracer(status, GPTRSTR, idt_H2O_liq_lg, I_CONVECT,   OFF)
      CALL tracer_halt(substr, status)
      CALL set_tracer(status, GPTRSTR, idt_H2O_liq_lg, I_INTEGRATE, ON)
      CALL tracer_halt(substr, status)
      CALL new_tracer(status, LGTRSTR, 'ICE', modstr      &
            ,unit = 'mol/mol'      &
            ,quantity = AMOUNTFRACTION  &
            ,medium = AIR &
            ,idx = idt_H2O_ice_lg)
      CALL tracer_halt(substr, status)

      CALL set_tracer(status, GPTRSTR, idt_H2O_ice_lg, R_molarmass, M_H2O)
      CALL tracer_halt(substr, status)
      CALL set_tracer(status, GPTRSTR, idt_H2O_ice_lg, I_ADVECT,    OFF)
      CALL tracer_halt(substr, status)
      CALL set_tracer(status, GPTRSTR, idt_H2O_ice_lg, I_VDIFF,     OFF)
      CALL tracer_halt(substr, status)
      CALL set_tracer(status, GPTRSTR, idt_H2O_ice_lg, I_CONVECT,   OFF)
      CALL tracer_halt(substr, status)
      CALL set_tracer(status, GPTRSTR, idt_H2O_ice_lg, I_INTEGRATE, ON)
      CALL tracer_halt(substr, status)
    END IF
!!#D attila -

    CALL end_message_bi(modstr,'TRACER REQUEST',substr)

  END SUBROUTINE h2o_new_tracer
!************************************************************************

!************************************************************************
  SUBROUTINE h2o_init_memory

    ! ECHAM5/MESSy
    USE messy_main_channel_error_bi, ONLY: channel_halt
    USE messy_main_channel_bi,       ONLY: GP_3D_MID, LG_ATTILA, SCALAR
    ! MESSy
    USE messy_main_channel,    ONLY: new_channel, new_channel_object &
                                   , new_attribute

    IMPLICIT NONE
    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'h2o_init_memory'
    INTEGER                     :: status

    CALL start_message_bi(modstr,'CHANNEL DEFINITION',substr)

!!#D attila +
    IF (L_LG) THEN
       CALL new_channel(status, modstr//'_lg', reprid=GP_3D_MID)
       CALL channel_halt(substr, status)

       CALL new_channel_object(status, modstr//'_lg', 'qte_adv' &
            , p3=qte_adv, reprid=GP_3D_MID)
       CALL channel_halt(substr, status)
       CALL new_attribute(status, modstr//'_lg', 'qte_adv' &
            , 'long_name', c='advection tendency')
       CALL channel_halt(substr, status)
       CALL new_attribute(status, modstr//'_lg', 'qte_adv' &
            , 'units', c='kg/kg/s')
       CALL channel_halt(substr, status)

       CALL new_channel_object(status, modstr//'_lg', 'rLG_H2O' &
            , p0=rLG_H2O, reprid=SCALAR)
       CALL channel_halt(substr, status)
       CALL new_attribute(status, modstr//'_lg', 'rLG_H2O' &
            , 'long_name', c='H2O flag')
       CALL channel_halt(substr, status)
       CALL new_attribute(status, modstr//'_lg', 'rLG_H2O' &
            , 'units', c='(1,0)')
       CALL channel_halt(substr, status)  
       rLG_H2O = 1._dp
    END IF
!!#D attila -

    ! required for both, LG and GP
    IF (I_H2O_TENDENCY > 0) THEN
       CALL new_channel(status, modstr//'_gp', reprid=GP_3D_MID)
       CALL channel_halt(substr, status)

       CALL new_channel_object(status, modstr//'_gp', 'CH4S', p3=CH4S_3d)
       CALL channel_halt(substr, status)
       CALL new_attribute(status, modstr//'_gp', 'CH4S' &
            , 'long_name', c='CH4 sink')
       CALL channel_halt(substr, status)
       CALL new_attribute(status, modstr//'_gp', 'CH4S' &
            , 'units', c='mol/mol/s')
       CALL channel_halt(substr, status)

!!$       CALL new_channel_object(status, modstr//'_gp', 'CH4_H', p3=CH4_H)
!!$       CALL channel_halt(substr, status)
!!$       CALL new_attribute(status, modstr//'_gp', 'CH4_H' &
!!$            , 'long_name', c='CH4 HALOE')
!!$       CALL channel_halt(substr, status)
!!$       CALL new_attribute(status, modstr//'_gp', 'CH4_H' &
!!$            , 'units', c='mol/mol')
!!$       CALL channel_halt(substr, status)
!!$
!!$       CALL new_channel_object(status, modstr//'_gp', 'CH4OH', p3=CH4OH)
!!$       CALL channel_halt(substr, status)
!!$       CALL new_attribute(status, modstr//'_gp', 'CH4OH' &
!!$            , 'long_name', c='CH4+OH loss rate')
!!$       CALL channel_halt(substr, status)
!!$       CALL new_attribute(status, modstr//'_gp', 'CH4OH', 'units', c='1/s')
!!$       CALL channel_halt(substr, status)
!!$
!!$       CALL new_channel_object(status, modstr//'_gp', 'CH4O1D', p3=CH4O1D)
!!$       CALL channel_halt(substr, status)
!!$       CALL new_attribute(status, modstr//'_gp', 'CH4O1D' &
!!$            , 'long_name', c='CH4+O1D loss rate')
!!$       CALL channel_halt(substr, status)
!!$       CALL new_attribute(status, modstr//'_gp', 'CH4O1D', 'units', c='1/s')
!!$       CALL channel_halt(substr, status)
!!$
!!$       CALL new_channel_object(status, modstr//'_gp', 'CH4Cl', p3=CH4Cl)
!!$       CALL channel_halt(substr, status)
!!$       CALL new_attribute(status, modstr//'_gp', 'CH4Cl' &
!!$            , 'long_name', c='CH4+Cl loss rate')
!!$       CALL channel_halt(substr, status)
!!$       CALL new_attribute(status, modstr//'_gp', 'CH4Cl', 'units', c='1/s')
!!$       CALL channel_halt(substr, status)

!!#D attila +
       IF (L_LG) THEN 
          CALL new_channel_object(status, modstr//'_lg' &
               , 'CH4S_LG', p1=CH4S_1d &
               , reprid=LG_ATTILA)
          CALL channel_halt(substr, status)
          CALL new_attribute(status, modstr//'_lg', 'CH4S_LG' &
               , 'long_name', c='CH4 sink')
          CALL channel_halt(substr, status)
          CALL new_attribute(status, modstr//'_lg' &
               , 'CH4S_LG', 'units', c='mol/mol/s')
          CALL channel_halt(substr, status)
       END IF
!!#D attila -
    END IF

    CALL end_message_bi(modstr,'CHANNEL DEFINITION',substr)

  END SUBROUTINE h2o_init_memory
!************************************************************************

!************************************************************************
  SUBROUTINE h2o_init_tracer

    ! ECHAM5/MESSy
    USE messy_main_transform_bi,    ONLY: zonal_average
    USE messy_main_blather_bi,      ONLY: error_bi
    USE messy_main_tracer_tools_bi, ONLY: tracer_halt
    USE messy_main_timer,           ONLY: lresume
    USE messy_main_grid_def_mem_bi, ONLY: nvclev, vct, nlev &
                                        , ngpblks, nproma, npromz

    USE messy_main_grid_def_bi,     ONLY: apzero
    USE messy_main_data_bi,         ONLY: q
    ! MESSy
    USE messy_main_tracer,          ONLY: get_tracer, tracer_iniflag

    IMPLICIT NONE

    ! LOCAL
    INTEGER :: status
    ! H2O TRACER
    ! zonal average q in decomp
    REAL(DP), DIMENSION(:,:,:), POINTER :: zaq => NULL()
    ! modulation around zaq
    REAL(DP), DIMENSION(:,:,:), POINTER :: modulation => NULL()
    INTEGER                             :: kproma, jrow, jk, kmerge

    ! TRACER INITIALISATION FROM netCDF-FILE
    ! (SKIP IN CASE I_H2O_TENDENCY == -2)

    ! RERUN: H2O HAS BEEN INITIALIZED FROM RERUN-FILE
    !        -> NO FURTHER CHANGES REQUIRED
    IF (lresume) RETURN

    ! CALCULATE LONGITUDINAL VARIATION (ZONAL AVERAGE +/- MODULATION)
    IF ( (I_H2O_TENDENCY > -2) .AND. L_MODINI_QZAV ) THEN
       !
       CALL zonal_average(q, zaq)
       !
       ALLOCATE(modulation(SIZE(q,1),SIZE(q,2),SIZE(q,3)))
       modulation(:,:,:) = 1.0_DP
       !
       DO jrow=1, ngpblks
          IF ( jrow == ngpblks ) THEN
             kproma = npromz
          ELSE
             kproma = nproma
          END IF
          
          modulation(1:kproma,:,jrow) = &
               q(1:kproma,:,jrow)/zaq(1:kproma,:,jrow)
       END DO
    END IF

    ! (A) OVERWRITE INITIALIZATION
    IF (I_H2O_TENDENCY == -2) THEN
       ! EVERYWHERE (NO INITIALIZATION)
       kmerge = 0
    ELSE
       ! SET EQUIVALENT TO q IN TROPOSPHERE (PRESSURE > 100 hPa)
       DO jk=1, nvclev 
          IF ( (vct(jk) + vct(nvclev+jk)*apzero) >= pmerge ) EXIT
       END DO
       kmerge = jk
    END IF

    IF (L_GP) CALL h2o_init_tracer_gp
!!#D attila +
    IF (L_LG) CALL h2o_init_tracer_lg
!!#D attila -

    ! CLEAN UP MEMORY
    IF (ASSOCIATED(modulation)) THEN
       DEALLOCATE(modulation)
       NULLIFY(modulation)
    END IF
    IF (ASSOCIATED(zaq)) THEN
       DEALLOCATE(zaq)
       NULLIFY(zaq)
    END IF

CONTAINS

    SUBROUTINE h2o_init_tracer_gp

      USE messy_main_tracer_mem_bi, ONLY: GPTRSTR

      IMPLICIT NONE

      ! H2O TRACER
      CHARACTER(LEN=*), PARAMETER         :: substr = 'h2o_init_tracer_gp'
      REAL(DP), DIMENSION(:,:,:), POINTER :: h2o         ! TRACER
      LOGICAL                             :: linit
      INTEGER                             :: jk

      ! SET POINTER TO TRACER
      CALL get_tracer(status, GPTRSTR, 'H2O', pxt=h2o, idx=idt_h2o_gp)
      IF (I_H2O_TENDENCY > -2) THEN
         CALL tracer_iniflag(status, GPTRSTR, idt_h2o_gp, lget=linit)
         CALL tracer_halt(substr, status)
         IF (.NOT. linit) THEN
            CALL error_bi('H2O TRACER WAS NOT INITIALISED!', substr)
         END IF
      END IF

      ! (A) OVERWRITE INITIALIZATION
      DO jk=kmerge+1, nlev
         h2o(:,jk,:) = scvmr * (q(:,jk,:)/(1.0_DP-q(:,jk,:)))
      END DO
      
      IF ((I_H2O_TENDENCY > -2) .AND. L_MODINI_QZAV) THEN
         ! (B) KEEP LONGITUDINAL VARIATION (ZONAL AVERAGE +/- MODULATION)
         !     IN STRATOSPHERE (PRESSURE < 100 hPa)
         DO jk=1, kmerge
            h2o(:,jk,:) = h2o(:,jk,:)*modulation(:,jk,:)
         END DO
      END IF
      
      ! (C) RE-SET NEGATIVES TO ZERO
      h2o(:,:,:) = MAX(h2o(:,:,:), 1.0E-40_DP)
      
    END SUBROUTINE h2o_init_tracer_gp

!!#D attila +
    ! -------------------------------------------------------------

    SUBROUTINE h2o_init_tracer_lg

      USE messy_main_tracer_mem_bi, ONLY: LGTRSTR
      USE messy_attila_tools_e5,    ONLY: gp2lg_e5, lg2gp_e5, LG2GP_AVE

      IMPLICIT NONE

      ! H2O TRACER
       CHARACTER(LEN=*), PARAMETER        :: substr = 'h2o_init_tracer_lg'
      REAL(DP), DIMENSION(:,:,:), POINTER :: h2o_tmp
      REAL(DP), DIMENSION(:,:,:), POINTER :: h2o_help    ! TRACER
      REAL(DP), DIMENSION(:),     POINTER :: h2o         ! TRACER
      LOGICAL                             :: linit
      INTEGER                             :: jk

      ! SET POINTER TO TRACER
      CALL get_tracer(status, LGTRSTR, 'H2O',pxt=h2o_help)
      h2o => h2o_help(:,1,1)

      IF (I_H2O_TENDENCY > -2) THEN
         CALL tracer_iniflag(status, LGTRSTR, idt_h2o_lg, lget=linit)
         CALL tracer_halt(substr, status)
         IF (.NOT. linit) THEN
            CALL error_bi('H2O TRACER WAS NOT INITIALISED!', substr)
         END IF
      END IF

      ! MEMORY
      ALLOCATE(h2o_tmp(nproma,nlev,ngpblks))

      ! (A) OVERWRITE INITIALIZATION
      CALL lg2gp_e5(h2o, h2o_tmp, LG2GP_AVE)
      DO jk=kmerge+1, nlev
         h2o_tmp(:,jk,:) = scvmr * (q(:,jk,:)/(1.0_DP-q(:,jk,:)))
      END DO
      CALL gp2lg_e5(h2o_tmp, h2o)
      
      IF ((I_H2O_TENDENCY > -2) .AND. L_MODINI_QZAV) THEN
         ! (B) KEEP LONGITUDINAL VARIATION (ZONAL AVERAGE +/- MODULATION)
         !     IN STRATOSPHERE (PRESSURE < 100 hPa)
         CALL lg2gp_e5(h2o, h2o_tmp, LG2GP_AVE)
         DO jk=1, kmerge
            h2o_tmp(:,jk,:) = h2o_tmp(:,jk,:)*modulation(:,jk,:)
         END DO
         CALL gp2lg_e5(h2o_tmp, h2o)
      END IF
      
      ! (C) RE-SET NEGATIVES TO ZERO
      h2o(:) = MAX(h2o(:), 1.0E-40_DP)

      ! FREE MEMORY
      IF (ASSOCIATED(h2o_tmp)) THEN
         DEALLOCATE(h2o_tmp)
         NULLIFY(h2o_tmp)
      END IF

    END SUBROUTINE h2o_init_tracer_lg
!!#D attila -

  END SUBROUTINE h2o_init_tracer
!************************************************************************

!************************************************************************
  SUBROUTINE h2o_global_start

    ! ECHAM5/MESSy
    USE messy_main_timer,         ONLY: lstart 
    USE messy_main_data_bi,       ONLY: qte_scb, q, qm1

    IMPLICIT NONE

!!$    IF (I_H2O_TENDENCY > 0) THEN
!!$       ! UPDATE FIELDS FOR CLIM. CHEM TEND. (IF REQUIRED)
!!$       CALL RGTEVENT_READ(RGT, modstr, 'A', RGREAD_NCVAR   &
!!$            ,CH4_H, lstop=.TRUE.)
!!$       CALL RGTEVENT_READ(RGT, modstr, 'B', RGREAD_NCVAR   &
!!$            ,CH4OH, lstop=.TRUE.)
!!$       CALL RGTEVENT_READ(RGT, modstr, 'C', RGREAD_NCVAR   &
!!$            ,CH4O1D, lstop=.TRUE.)
!!$       CALL RGTEVENT_READ(RGT, modstr, 'D', RGREAD_NCVAR   &
!!$            ,CH4CL, lstop=.TRUE.)
!!$    END IF

    IF (L_GP) CALL h2o_global_start_gp

!!#D attila +
    IF (L_LG) CALL h2o_global_start_lg ! mz_pj_20080201
!!#D attila -

  CONTAINS

    SUBROUTINE h2o_global_start_gp

      USE messy_main_tracer_mem_bi, ONLY: xt, xtte, xtm1

      IMPLICIT NONE

      ! H2O TRACER
      REAL(DP), DIMENSION(:,:,:), POINTER :: h2o         ! TRACER
      REAL(DP), DIMENSION(:,:,:), POINTER :: h2om1       ! TRACER at t-1
      REAL(DP), DIMENSION(:,:,:), POINTER :: h2ote       ! TRACER-TENDENCY
      
      ! SET POINTER TO TRACER AND TENDENCY
      h2o   => xt(:,:,idt_h2o_gp,:)
      h2om1 => xtm1(:,:,idt_h2o_gp,:)
      h2ote => xtte(:,:,idt_h2o_gp,:)
      
      ! CONVERT SPECIFIC HUMIDITY TO H2O TRACER (INCL. TENDENCY),
      ! EXCEPT AT FIRST MODEL TIME STEP (lstart),
      ! WHERE H2O TRACER HAS BEEN INITIALIZED CORRECTLY, BUT q IS WRONG
      IF (lstart) THEN
         ! ALWAYS CONVERT CORRECTLY INITIALIZED TRACER TO SPECIFIC HUMIDITY
         q(:,:,:) = h2o(:,:,:)/(scvmr+h2o(:,:,:))
         qm1(:,:,:) = h2om1(:,:,:)/(scvmr+h2om1(:,:,:))
         ! 
      ELSE
         IF (I_FEEDBACK /= 0) THEN  ! op_sb_20130301
            ! CONVERT SPECIFIC HUMIDITY TO TRACER
            h2o(:,:,:)   = scvmr * (q(:,:,:)/(1.0_DP-q(:,:,:)))
            h2om1(:,:,:) = scvmr * (qm1(:,:,:)/(1.0_DP-qm1(:,:,:)))
         END IF                     ! op_sb_20130301
      END IF
      
      ! CONSISTENT TENDENCIES (BEFORE CHEMISTRY)
      !
      !         d(H2O)   d             q                   dq/dt
      ! H2Ote = ------ = --  (scvmr * ----- ) = scvmr * -----------
      !           dt     dt           1 - q               (1-q)^2
      !
      h2ote(:,:,:) = scvmr * qte_scb(:,:,:)/(1.0_DP - q(:,:,:))**2
      
#ifdef MESSYTENDENCY
      if (.not. allocated(lo_h2ote)) then
         allocate(lo_h2ote(size(h2ote,1),size(h2ote,2),size(h2ote,3)))
      endif

      lo_h2ote = 0.0_dp
      lo_h2ote(:,:,:) = h2ote(:,:,:)
#endif

    END SUBROUTINE h2o_global_start_gp

    ! ------------------------------------------------
!!#D attila +
    SUBROUTINE h2o_global_start_lg

      USE messy_main_timer,         ONLY: lstart, ztmst=>time_step_len
      USE messy_main_grid_def_mem_bi, ONLY: nproma, nlev, ngpblks
      USE messy_main_tracer_mem_bi, ONLY: xt_a, xtm1_a
      USE messy_attila_tools_e5,    ONLY: lg2gp_e5, LG2GP_AVE

      IMPLICIT NONE

      ! LOCAL
      REAL(DP), DIMENSION(:),     POINTER :: h2o_lg    => NULL()  ! TR
      REAL(DP), DIMENSION(:),     POINTER :: h2om1_lg  => NULL()  ! TR at t-1
      REAL(DP), DIMENSION(:,:,:), POINTER :: h2o_tmp   => NULL()  ! TR
      REAL(DP), DIMENSION(:,:,:), POINTER :: h2om1_tmp => NULL()  ! TR at t-1
      REAL(DP), DIMENSION(:,:,:), POINTER :: fillf     => NULL()
      REAL(DP), DIMENSION(:,:,:), POINTER :: h2om1lt_tmp => NULL() ! TR at t-1
      REAL(DP), DIMENSION(:,:,:), POINTER :: qm1vlt => NULL()

      REAL(DP), DIMENSION(:),     POINTER :: h2o_liq_lg   => NULL()
      REAL(DP), DIMENSION(:),     POINTER :: h2om1_liq_lg => NULL()

      ! SET POINTERS
      h2o_lg   => xt_a(:,1,idt_h2o_lg,1)
      h2om1_lg => xtm1_a(:,1,idt_h2o_lg,1)

      h2o_liq_lg   => xt_a(:,1,idt_h2o_liq_lg,1)
      h2om1_liq_lg => xtm1_a(:,1,idt_h2o_liq_lg,1)

      ! MEMORY
      ALLOCATE(h2o_tmp(nproma,nlev,ngpblks))
      ALLOCATE(h2om1_tmp(nproma,nlev,ngpblks))
      ALLOCATE(fillf(nproma, nlev, ngpblks))
      ALLOCATE(qm1vlt(nproma, nlev, ngpblks))

      SELECT CASE(I_FEEDBACK)

      CASE(0)
         !print*,'in h2o i-feedback =0 in global start'
         IF (lstart) THEN
            ! CONVERT CORRECTLY INITIALIZED TRACER TO SPECIFIC HUMIDITY
            fillf(:,:,:) = scvmr * (q(:,:,:)/(1.0_DP-q(:,:,:)))
            CALL lg2gp_e5(h2o_lg, h2o_tmp, LG2GP_AVE, fill_field=fillf)
            q(:,:,:) = h2o_tmp(:,:,:)/(scvmr+h2o_tmp(:,:,:))

            fillf(:,:,:) = scvmr * (qm1(:,:,:)/(1.0_DP-qm1(:,:,:)))
            CALL lg2gp_e5(h2om1_lg, h2om1_tmp, LG2GP_AVE, fill_field=fillf)
            qm1(:,:,:) = h2om1_tmp(:,:,:)/(scvmr+h2om1_tmp(:,:,:))
         END IF

      !CASE(1)
         ! nothing to do

      CASE(2)
         qte_adv(:,:,:) = 0._dp

         IF (lstart) THEN
            ! CONVERT CORRECTLY INITIALIZED TRACER TO SPECIFIC HUMIDITY
            fillf(:,:,:) = scvmr * (q(:,:,:)/(1.0_DP-q(:,:,:)))
            CALL lg2gp_e5(h2o_lg, h2o_tmp, LG2GP_AVE, fill_field=fillf)
            q(:,:,:) = h2o_tmp(:,:,:)/(scvmr+h2o_tmp(:,:,:))
            
            fillf(:,:,:) = scvmr * (qm1(:,:,:)/(1.0_DP-qm1(:,:,:)))
            CALL lg2gp_e5(h2om1_lg, h2om1_tmp, LG2GP_AVE, fill_field=fillf)
            qm1(:,:,:) = h2om1_tmp(:,:,:)/(scvmr+h2om1_tmp(:,:,:))
            ! 
         ELSE
            IF (L_LGMC) THEN

               ALLOCATE(h2om1lt_tmp(nproma,nlev,ngpblks))

               ! at the beginning of the time loop syncronize  lg ==> gp
               !
               !---------------
               ! water vapour -
               !---------------
               !
               ! SAVE ADVECTIVE ATTILA TENDENCY (WHICH INCLUDES THE MC-DIFFUSION
               ! AS EQUIVALENT TO VDIFF) FOR LATER USE TO CALULCATE
               ! HUMIDITY CONVERGENCE (USED IN CONVECTION)
               !
               ! 1. old humidity field at old positions
               fillf(:,:,:) = scvmr * (qm1(:,:,:)/(1.0_DP-qm1(:,:,:)))
               CALL lg2gp_e5(h2om1_lg,h2om1lt_tmp,LG2GP_AVE,ltm1=(.not.lstart) &
                    ,fill_field=fillf)

               ! 2. old humidity field at new positions
               !    ( this is also used to synchronize q with H2O_lg)
               fillf(:,:,:) = scvmr * (qm1(:,:,:)/(1.0_DP-qm1(:,:,:)))
               CALL lg2gp_e5(h2om1_lg, h2om1_tmp, LG2GP_AVE, fill_field=fillf)

               ! 3. new humidity field at new positions
               !   ( this is only used to synchronize q with H2O_lg)
               fillf(:,:,:) = scvmr * (q(:,:,:)/(1.0_DP-q(:,:,:)))
               CALL lg2gp_e5(h2o_lg, h2o_tmp, LG2GP_AVE, fill_field=fillf)

               ! CONVERT CORRECTLY INITIALIZED TRACER TO SPECIFIC HUMIDITY
               qm1(:,:,:) = h2om1lt_tmp(:,:,:)/(scvmr+h2om1lt_tmp(:,:,:))
               qm1vlt(:,:,:) = h2om1_tmp(:,:,:)/(scvmr+h2om1_tmp(:,:,:))
               q(:,:,:) = h2o_tmp(:,:,:)/(scvmr+h2o_tmp(:,:,:))

               ! last step: save advective (attila) tendency
               !    (including mc-diff == vdiff)
               qte_adv(:,:,:)  = (qm1vlt(:,:,:) - qm1(:,:,:))/ztmst 
            ENDIF 

            DEALLOCATE(h2om1lt_tmp) ; NULLIFY(h2om1lt_tmp)

         END IF

         ! NOTE: In LG it is not desired to convert from GP back to LG
         !       (diffusion !). Since a lagrangian q does not exist, no
         !       update of h2ote_lg is required.

         qte_scb(:,:,:)  = qte_adv(:,:,:)

      END SELECT

      DEALLOCATE(qm1vlt) ;      NULLIFY(qm1vlt)
      DEALLOCATE(fillf)  ;      NULLIFY(fillf)
      DEALLOCATE(h2o_tmp) ;     NULLIFY(h2o_tmp)
      DEALLOCATE(h2om1_tmp) ;   NULLIFY(h2om1_tmp)

    END SUBROUTINE h2o_global_start_lg
!!#D attila -

  END SUBROUTINE h2o_global_start
!************************************************************************

!************************************************************************
  SUBROUTINE h2o_physc

    ! ECHAM5
    ! ECHAM5/MESSy
    USE messy_main_grid_def_mem_bi,      ONLY: jrow, kproma

    IMPLICIT NONE

    IF (I_H2O_TENDENCY <= 0) RETURN

    ! UPDATE METHANE LOSS RATE [mol/mol/s]
    ch4s_3d(1:kproma,:,jrow) = ch4_H(1:kproma,:,jrow) *       &
         (ch4oh(1:kproma,:,jrow) + ch4o1d(1:kproma,:,jrow) +  &
         ch4cl(1:kproma,:,jrow))

    ! NOTE: OTHER SUBMODELS CAN POTENTIALLY CHANGE THE
    !       H2O-TRACER TENDENCY AT THIS LEVEL; THE DECISION
    !       WHICH TENDENCY TO FEEDBACK INTO q MUST BE TAKEN
    !       LATER --> h2o_local_end or h2o_global_end ???

  END SUBROUTINE h2o_physc
!************************************************************************

!************************************************************************
  SUBROUTINE h2o_local_end

    IMPLICIT NONE

    IF (L_GP) CALL h2o_local_end_gp
!!#D attila +
    ! FOR LG THIS NEEDS TO BE PERFORMED IN global_end ...
!!#D attila -

  CONTAINS

    SUBROUTINE h2o_local_end_gp

      ! ECHAM5/MESSy
      USE messy_main_grid_def_mem_bi,       ONLY: jrow, kproma
      USE messy_main_tracer_mem_bi, ONLY: xtte=>qxtte

      IMPLICIT NONE

      ! LOCAL
      ! H2O TRACER TENDENCY
      REAL(DP), DIMENSION(:,:), POINTER :: h2ote

      ! SET POINTER TO TRACER AND TENDENCY
      h2ote => xtte(:,:,idt_h2o_gp)
      
      ! 'RE-SYNCHRONIZE' H2O-TENDENCY AND q-TENDENCY
      !
      ! NOTES:
      ! (1) q AND H2O MUST NOT HAVE BEEN CHANGED BY ANY SUBMODEL
      !     AND MUST THEREFORE BE EQUIVALENT TO EACH OTHER !
      !
      ! (2) qte IS NOT CHANGED BY ANY OTHER SUBMODEL (ONLY ECHAM5)
      ! 
      ! (3) THE TRACER TENDENCY H2OTE IS PROBABLY CHANGED BY OTHER
      !     SUBMODELS
      !
      ! =>  DECIDE HERE (VIA NAMELIST SWITCH 'I_H2O_TENDENCY') WHICH
      !     TENDENCY SHOULD BE FEEDBACK TO THE SPECIFIC HUMIDITY
      
      ! STEP 1: UPDATE H2O TRACER TENDENCY
      !         (NOTE: H2O PRODUCTION APPROX 2*CH4 LOSS)
      SELECT CASE(I_H2O_TENDENCY)
      CASE(-2, -1)
         ! NO FEEDBACK AT ALL: BLANK ALL TENDENCIES
         h2ote(1:kproma,:) = 0.0
      CASE(0) 
         ! NO CLIMATOLOGICAL TENDENCY; KEEP THOSE FROM OTHER SUBMODELS
      CASE(1)
         ! ADD CLIMATOLOGICAL TENDENCY TO THOSE FROM OTHER SUBMODELS
         h2ote(1:kproma,:) = h2ote(1:kproma,:) &
              + 2. * ch4s_3d(1:kproma,:, jrow)
      CASE(2)
         ! USE CLIMATOLOGICAL TENDENCY EXCLUSIVELY AS CHEM. TEND.!
         ! !!! OVERWRITE: IGNORE ALL TENDENCIES, WHICH HAVE PROBABLY
         !                BEEN CALCULATED BY OTHER SUBMODELS
         h2ote(1:kproma,:) = 2. * ch4s_3d(1:kproma,:,jrow)
      CASE DEFAULT
         ! ERROR
      END SELECT

    END SUBROUTINE h2o_local_end_gp

    ! ------------------------------------------------------

  END SUBROUTINE h2o_local_end
!************************************************************************

!************************************************************************
  SUBROUTINE h2o_global_end

    ! ECHAM5/MESSy
    USE messy_main_grid_def_mem_bi, ONLY: nproma, ngpblks, nlev
    USE messy_main_data_bi,       ONLY: qte_scb, q, qm1       &
                                      , xl, xlte_scb, xlm1
    USE messy_main_timer,         ONLY: ztmst => time_step_len
    USE messy_main_tracer_mem_bi, ONLY: xt, xtte, xtm1, xt_a, xtte_a, xtm1_a &
                                      , NCELL
!!#D attila +
    USE messy_attila_tools_e5,    ONLY: gp2lg_e5, lg2gp_e5, LG2GP_AVE &
                                      , lggpte2lgte_e5
!!#D attila -

    IMPLICIT NONE

    ! LOCAL
    REAL(DP), DIMENSION(:,:,:), POINTER :: zq
    ! H2O TRACER TENDENCY
    REAL(DP), DIMENSION(:,:,:), POINTER :: h2ote_gp => NULL()
    REAL(DP), DIMENSION(:,:,:), POINTER :: h2o_gp   => NULL()
    REAL(DP), DIMENSION(:,:,:), POINTER :: h2om1_gp => NULL()
    !
    REAL(DP), DIMENSION(:)    , POINTER :: h2ote_lg => NULL()
    REAL(DP), DIMENSION(:)    , POINTER :: h2o_lg => NULL()
    REAL(DP), DIMENSION(:)    , POINTER :: h2om1_lg => NULL()
    !
    REAL(DP), DIMENSION(:,:,:), POINTER :: fillf     => NULL()
    !
    REAL(DP), DIMENSION(:,:,:), POINTER :: h2ote_tmp => NULL()
    REAL(DP), DIMENSION(:,:,:), POINTER :: h2ovn_gp => NULL()
    REAL(DP), DIMENSION(:,:,:), POINTER :: h2ote_scb_gp => NULL()
    !
    REAL(DP), DIMENSION(:)    , POINTER :: h2ovn_lg => NULL()
    REAL(DP), DIMENSION(:)    , POINTER :: h2ote_scb_lg => NULL()
    ! LIQ TRACER TENDENCY
    !
    REAL(DP), DIMENSION(:)    , POINTER :: liqte_lg => NULL()
    REAL(DP), DIMENSION(:)    , POINTER :: liq_lg   => NULL()
    REAL(DP), DIMENSION(:)    , POINTER :: liqm1_lg => NULL()
    !
    REAL(DP), DIMENSION(:,:,:), POINTER :: liqte_tmp    => NULL()
    REAL(DP), DIMENSION(:,:,:), POINTER :: liqvn_gp     => NULL()
    REAL(DP), DIMENSION(:,:,:), POINTER :: liqte_scb_gp => NULL()
    !
    REAL(DP), DIMENSION(:)    , POINTER :: liqvn_lg     => NULL()
    REAL(DP), DIMENSION(:)    , POINTER :: liqte_scb_lg => NULL()

#ifdef MESSYTENDENCY
    REAL(kind=dp),DIMENSION(size(qte_scb,1),size(qte_scb,2),size(qte_scb,3)) :: lo_qte

    lo_qte           =  0.0_dp
    lo_qte(:,:,:)    =  qte_scb(:,:,:)
#endif

    ! SET LOCAL POINTER
    zq    => q(:,:,:)

    SELECT CASE(I_FEEDBACK)
       !CASE(0) ! no feedback

       CASE(1) ! GP
          !
          ! SET GP POINTER TO TRACER AND TENDENCY
          h2ote_gp => xtte(:,:,idt_h2o_gp,:)
          h2o_gp   => xt(:,:,idt_h2o_gp,:)
          h2om1_gp => xtm1(:,:,idt_h2o_gp,:)
          !
          ! STEP 2: ADD SELECTED TRACER TENDENCY TO qte FOR FEEDBACK
          !
          !        dq    d      H2O           scvmr * d(H2O)/dt
          ! qte = ---- = -- (-----------) = ---------------------
          !        dt    dt  scvmr + H2O      (scvmr + H2O)^2
          !
          qte_scb(:,:,:) = qte_scb(:,:,:) + &
               scvmr * h2ote_gp(:,:,:) / (scvmr + h2o_gp(:,:,:))**2
          !
          ! STEP 3: RE-SYNCHRONIZE H2O AND q TENDENCIES FOR TIME INTEGRATION
          !         (CONSISTENT TENDENCIES (AFTER CHEMISTRY))
          h2ote_gp(:,:,:) = scvmr * qte_scb(:,:,:)/(1.0_DP - zq(:,:,:))**2
          !
          !
!!#D attila +
          ! STEP 4: SYNCHRONIZE LG TRACER IF IT IS PRESENT
          IF (L_LG) THEN
             h2o_lg   => xt_a(:,1,idt_h2o_lg,1)
             h2om1_lg => xtm1_a(:,1,idt_h2o_lg,1)
             h2ote_lg => xtte_a(:,1,idt_h2o_lg,1)
             ! mim_sb_20080605+
             liq_lg   => xt_a(:,1,idt_h2o_liq_lg,1)
             liqm1_lg => xtm1_a(:,1,idt_h2o_liq_lg,1)
             liqte_lg => xtte_a(:,1,idt_h2o_liq_lg,1)
             ! mim_sb_20080605-
             !
             CALL gp2lg_e5(h2ote_gp, h2ote_lg)
             CALL gp2lg_e5(h2om1_gp, h2om1_lg)
             CALL gp2lg_e5(h2o_gp, h2o_lg)
             ! mim_sb_20080605+
             CALL gp2lg_e5(xlte_scb, liqte_lg)
             CALL gp2lg_e5(xlm1, liqm1_lg)
             CALL gp2lg_e5(xl, liq_lg)
             ! mim_sb_20080605-
          END IF
          !
       CASE(2) ! LG
          !
          CALL h2o_global_end_lg
          ! 
          ! SET LG POINTER TO TRACER AND TENDENCY
          h2ote_lg => xtte_a(:,1,idt_h2o_lg,1)
          h2om1_lg => xtm1_a(:,1,idt_h2o_lg,1)
          h2o_lg   => xt_a(:,1,idt_h2o_lg,1)

          liqte_lg => xtte_a(:,1,idt_h2o_liq_lg,1)
          liqm1_lg => xtm1_a(:,1,idt_h2o_liq_lg,1)
          liq_lg   => xt_a(:,1,idt_h2o_liq_lg,1)
          !
          ! MEMORY
          ALLOCATE(h2ote_tmp(nproma,nlev,ngpblks))
          ALLOCATE(h2ovn_gp(nproma,nlev,ngpblks))
          ALLOCATE(h2ote_scb_gp(nproma,nlev,ngpblks))
          ALLOCATE(h2ovn_lg(NCELL))
          ALLOCATE(h2ote_scb_lg(NCELL))

          ALLOCATE(liqte_tmp(nproma,nlev,ngpblks))
          ALLOCATE(liqvn_gp(nproma,nlev,ngpblks))
          ALLOCATE(liqte_scb_gp(nproma,nlev,ngpblks))
          ALLOCATE(liqvn_lg(NCELL))
          ALLOCATE(liqte_scb_lg(NCELL))

          ! ADD GP-tendency of specific humidity to LG-tendency

          ! Tendency of humidity mixing ratio in gp
          ! substract in case of moist convection grid point convection 
          ! tendencies
          if (L_LGMC) then
             h2ote_scb_gp(:,:,:) = scvmr * (qte_scb(:,:,:) - conv_qte(:,:,:)) &
                  /(1.0_DP - q(:,:,:))**2
          else
             h2ote_scb_gp(:,:,:) = scvmr * qte_scb(:,:,:) &
                  /(1.0_DP - q(:,:,:))**2
          endif

          ! updated values for h2o variables in lg
          ! (for tendency interpolation/scaling (GP->LG) necessary)
          h2ovn_lg(:) = h2om1_lg(:) + h2ote_lg(:) * ztmst
          liqvn_lg(:) = liqm1_lg(:) + liqte_lg(:) * ztmst ! mim_sb_20091009
          !
          ALLOCATE(fillf(nproma, nlev, ngpblks))
          fillf(:,:,:) = scvmr * (qm1(:,:,:)/(1.0_DP-qm1(:,:,:)))
          CALL lg2gp_e5(h2ovn_lg,h2ovn_gp,LG2GP_AVE, &
                       fill_field=fillf)

          fillf(:,:,:) = 0.
          CALL lg2gp_e5(liqvn_lg,liqvn_gp,LG2GP_AVE, &
                       fill_field=fillf)

          DEALLOCATE(fillf)
          NULLIFY(fillf)

          ! Interpolation of gp-tendency to lg
          CALL lggpte2lgte_e5(h2ovn_gp, h2ote_scb_gp, h2ovn_lg, h2ote_scb_lg)
          CALL lggpte2lgte_e5(liqvn_gp, liqte_scb_gp, liqvn_lg, liqte_scb_lg)

          ! Interpolation of lg-tendeny to gp
          CALL lg2gp_e5(h2ote_lg, h2ote_tmp, LG2GP_AVE)
          ! mim_sb_20091009+
          !!! Vorsicht !!!
          ! nach Konvektion dürfen die Tendenzen so nicht transformiert werden!!!
!          CALL lg2gp_e5(liqte_lg, liqte_tmp, LG2GP_AVE)
          ! mim_sb_20091009-

          ! Total tendency in lg
          h2ote_lg(:) = h2ote_lg(:) + h2ote_scb_lg(:)

          ! Total tendency in gp
          qte_scb(:,:,:) = qte_scb(:,:,:) +      &
               scvmr * h2ote_tmp(:,:,:)          &
               / (scvmr + h2ovn_gp(:,:,:))**2

          ! STEP 4: SYNCHRONIZE GP TRACER IF IT IS PRESENT 
          IF (L_GP) THEN
             h2ote_gp => xtte(:,:,idt_h2o_gp,:)
             h2o_gp   => xt(:,:,idt_h2o_gp,:)
             h2om1_gp => xtm1(:,:,idt_h2o_gp,:)
             !
             h2ote_gp(:,:,:) = scvmr * qte_scb(:,:,:)/(1.0_DP - q(:,:,:))**2
             h2o_gp(:,:,:)   = scvmr * (q(:,:,:)/(1.0_DP-q(:,:,:)))
             h2om1_gp(:,:,:) = scvmr * (qm1(:,:,:)/(1.0_DP-qm1(:,:,:)))
          END IF

          !
          DEALLOCATE(h2ote_tmp)  ; NULLIFY(h2ote_tmp)
          DEALLOCATE(h2ote_scb_gp); NULLIFY(h2ote_scb_gp)
          DEALLOCATE(h2ovn_gp)  ; NULLIFY(h2ovn_gp)
          DEALLOCATE(h2ovn_lg)  ; NULLIFY(h2ovn_lg)
          DEALLOCATE(h2ote_scb_lg) ; NULLIFY(h2ote_scb_lg)

          DEALLOCATE(liqte_tmp)  ; NULLIFY(liqte_tmp)
          DEALLOCATE(liqte_scb_gp); NULLIFY(liqte_scb_gp)
          DEALLOCATE(liqvn_gp)  ; NULLIFY(liqvn_gp)
          DEALLOCATE(liqvn_lg)  ; NULLIFY(liqvn_lg)
          DEALLOCATE(liqte_scb_lg) ; NULLIFY(liqte_scb_lg)

!!#D attila -
          !
       END SELECT

#ifdef MESSYTENDENCY
        lo_qte(:,:,:)   = qte_scb(:,:,:) - lo_qte(:,:,:)
        qte_scb(:,:,:)  = qte_scb(:,:,:) - lo_qte(:,:,:)

        lo_h2ote(:,:,:) = h2ote_gp(:,:,:) - lo_h2ote(:,:,:)
        h2ote_gp(:,:,:) = h2ote_gp(:,:,:) - lo_h2ote(:,:,:)

        call mtend_add_g (my_handle, mtend_id_q, px = lo_qte)
        call mtend_add_g (my_handle, idt_h2o_gp, px = lo_h2ote)
#endif

!!#D attila +
CONTAINS
    ! -----------------------------------------------

    SUBROUTINE h2o_global_end_lg

      ! ECHAM5/MESSy
      USE messy_attila_tools_e5,    ONLY: gp2lg_e5
      USE messy_main_tracer_mem_bi, ONLY: xt_a, xtte_a, xtm1_a

      IMPLICIT NONE

      REAL(DP), DIMENSION(:),     POINTER :: h2ote   => NULL()
      REAL(DP), DIMENSION(:),     POINTER :: h2om1   => NULL()
      REAL(DP), DIMENSION(:),     POINTER :: h2o     => NULL()

      ! PHYSC:
      ! TRANSFORM CLIMATOLOGICAL CHEMICAL TENDENCY
      IF (I_H2O_TENDENCY > 0) THEN
         CALL gp2lg_e5(ch4s_3d,ch4s_1d)
      END IF

      ! LOCAL_END:
      ! - SET POINTER TO TRACER AND TENDENCY
      h2ote => xtte_a(:,1,idt_h2o_lg,1)
      h2om1 => xtm1_a(:,1,idt_h2o_lg,1)
      h2o   => xt_a(:,1,idt_h2o_lg,1)

      ! STEP 1: UPDATE H2O TRACER TENDENCY
      !         (NOTE: H2O PRODUCTION APPROX 2*CH4 LOSS)
      SELECT CASE(I_H2O_TENDENCY)
      CASE(-2, -1)
         ! NO FEEDBACK AT ALL: BLANK ALL TENDENCIES
         h2ote(:) = 0.0
      CASE(0) 
         ! NO CLIMATOLOGICAL TENDENCY; KEEP THOSE FROM OTHER SUBMODELS
      CASE(1)
         ! ADD CLIMATOLOGICAL TENDENCY TO THOSE FROM OTHER SUBMODELS
         h2ote(:) = h2ote(:) + 2. * ch4s_1d(:)
      CASE(2)
         ! USE CLIMATOLOGICAL TENDENCY EXCLUSIVELY AS CHEM. TEND.!
         ! !!! OVERWRITE: IGNORE ALL TENDENCIES, WHICH HAVE PROBABLY
         !                BEEN CALCULATED BY OTHER SUBMODELS
         h2ote(:) = 2. * ch4s_1d(:)
      CASE DEFAULT
         ! ERROR
      END SELECT

    END SUBROUTINE h2o_global_end_lg
!!#D attila -

  END SUBROUTINE h2o_global_end
!************************************************************************

!************************************************************************
  SUBROUTINE h2o_free_memory

    IMPLICIT NONE

#ifdef MESSYTENDENCY
      IF (ALLOCATED(lo_h2ote)) THEN
         DEALLOCATE(lo_h2ote)
      ENDIF
#endif

  END SUBROUTINE h2o_free_memory
!************************************************************************

! ************************************************************************
! PRIVATE ECHAM-5 INTERFACE ROUTINES
! ************************************************************************

! ************************************************************************
  SUBROUTINE h2o_read_nml_cpl(status, iou)

    ! H2O MODULE ROUTINE (ECHAM-5 INTERFACE, PRIVATE)
    !
    ! read namelist for 'coupling' to ECHAM5
    !
    ! Author: Patrick Joeckel, MPICH, Nov 2003

    ! MESSy
    USE messy_main_tools, ONLY: read_nml_open, read_nml_check, read_nml_close
    USE messy_main_tracer_mem_bi, ONLY: NGCELL

    IMPLICIT NONE

    ! I/O
    INTEGER, INTENT(OUT) :: status     ! error status
    INTEGER, INTENT(IN)  :: iou        ! I/O unit

    ! (LOCAL) NAMELIST VARIABLES
    CHARACTER(LEN=*), PARAMETER :: substr = 'h2o_read_nml_cpl'

    NAMELIST /CPL/ I_H2O_TENDENCY, L_MODINI_QZAV &
         , L_GP, L_LG, I_FEEDBACK

    ! LOCAL
    LOGICAL          :: lex      ! file exists ?
    INTEGER          :: fstat    ! file status

    status = 1

    CALL read_nml_open(lex, substr, iou, 'CPL', modstr)
    IF (.not.lex) RETURN    ! <modstr>.nml does not exist

    READ(iou, NML=CPL, IOSTAT=fstat)
    CALL read_nml_check(fstat, substr, iou, 'CPL', modstr)
    IF (fstat /= 0) RETURN  ! error while reading namelist

    ! DIAGNOSE NAMELIST AND SET GLOBAL SWITCHES
    !
    ! CHECK NAMELIST
    IF (L_LG) THEN
       IF (NGCELL > 0) THEN
!!#D attila +
          WRITE(*,*) 'LAGRANGIAN CALCULATION OF H2O: ON'
!!#D attila -
       ELSE
          IF (L_LG) THEN
!!#D attila +
             WRITE(*,*) 'L_LG = T in namelist'
             WRITE(*,*) 'However no Lagrangian scheme activated ...'
             WRITE(*,*) ' ... setting L_LG = F'
!!#D attila -
          END IF
          L_LG = .FALSE.
       END IF
    ELSE
!!#D attila +
       WRITE(*,*) 'LAGRANGIAN CALCULATION OF H2O: OFF'
!!#D attila -
    END IF

    IF (L_GP) THEN
       WRITE(*,*) 'GRIDPOINT  CALCULATION OF H2O: ON'
    ELSE
       WRITE(*,*) 'GRIDPOINT  CALCULATION OF H2O: OFF'
    END IF

    IF ((.NOT.L_LG) .AND. (.NOT.L_GP)) THEN
      CALL error_bi( &
           ' ERROR: NEITHER GP NOR LG CALCULATION OF H2O ACTIVATED ... ' &
           , substr) 
    END IF

    SELECT CASE(I_H2O_TENDENCY)
    CASE(-2)
       WRITE(*,*) 'NO FEEDBACK: ADD ZERO H2O-TRACER-TENDENCY'
       WRITE(*,*) '             (NOTE: TRANSPORT TENDENCIES IN Q ARE KEPT)'
       WRITE(*,*) 'NO INITIALIZATION IN STRATOSPHERE !'
    CASE(-1)
       WRITE(*,*) 'NO FEEDBACK: ADD ZERO H2O-TRACER-TENDENCY'
       WRITE(*,*) '             (NOTE: TRANSPORT TENDENCIES IN Q ARE KEPT)'
    CASE(0)
       WRITE(*,*) 'FEEDBACK: ADD H2O-TRACER TENDENCIES'
       WRITE(*,*) '          POTENTIALLY CALCULATED BY OTHER SUBMODELS'
    CASE(1)
       WRITE(*,*) 'FEEDBACK: ADD CLIMATOLOGICAL CHEMICAL H2O-TRACER TENDENCY'
       WRITE(*,*) '          (CH4 + OH,O1D,Cl) AND H2O-TRACER TENDENCIES'
       WRITE(*,*) '          POTENTIALLY CALCULATED BY OTHER SUBMODELS'
    CASE(2)
       WRITE(*,*) 'FEEDBACK: ADD EXCLUSIVELY CLIMATOLOGICAL CHEMICAL'
       WRITE(*,*) '          H2O-TRACER TENDENCY; IGNORE ALL OTHERS'
       WRITE(*,*) '          POTENTIALLY CALCULATED BY OTHER SUBMODELS'
    CASE DEFAULT
       WRITE(*,*) 'ERROR: UNRECOGNIZED VALUE FOR I_H2O_TENDENCY'
       RETURN
    END SELECT

    IF (L_MODINI_QZAV) THEN
       WRITE(*,*) 'MODULATION OF INITIAL H2O CLIMATOLOGY WITH Q/Q_zave: ON'
    ELSE
       WRITE(*,*) 'MODULATION OF INITIAL H2O CLIMATOLOGY WITH Q/Q_zave: OFF'
    END IF

    SELECT CASE(I_FEEDBACK)
    ! op_sb_20130301+
    CASE(0)
        IF (.NOT. L_LG) THEN
          CALL error_bi('LG OFFLINE HYDROLOGICAL CYCLE IMPOSSIBLE FOR L_LG=F' &
               , substr)
        ELSE
          WRITE(*,*) 'LAGRANGIAN HYDROLOGICAL CYCLE WITHOUT FEEDBACK IS SELECTED'
        END IF
    ! op_sb_20130301-
    CASE(1)
       IF (.NOT. L_GP) THEN
          CALL error_bi('GP FEEDBACK IMPOSSIBLE FOR L_GP=F', substr)
       ELSE
          WRITE(*,*) 'FEEDBACK FROM GRIDPOINT H2O'
       END IF
    CASE(2)
       IF (.NOT. L_LG) THEN
          CALL error_bi('LG FEEDBACK IMPOSSIBLE FOR L_LG=F', substr)
       ELSE
          WRITE(*,*) 'FEEDBACK FROM LAGRANGIAN H2O'
       END IF
    CASE DEFAULT
          CALL error_bi('UNKNOWN FEEDBACK OPTION', substr)
    END SELECT

    CALL read_nml_close(substr, iou, modstr)

    status = 0  ! no ERROR

  END SUBROUTINE h2o_read_nml_cpl

! ************************************************************************
END MODULE messy_h2o_e5
! ************************************************************************ 

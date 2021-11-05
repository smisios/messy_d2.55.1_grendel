! ====================================================================
!
! TITLE    : H2OEMIS
! 
! PROJECT  : STRATOFLY
!
! Module   : smil/messy_h2oemis_si
!
! DATE     : 01.10.2019
!
! NOTE     : Comment format '!>' is used for generation of code 
!              documentation from source file via 'doxygen'.
! 
!> \brief
!> This is  documentation for the  Submodel Interface Layer (SMIL) routines. 
!> Here, the submodel adds H2O emissions to specific humidity, to keep them in the
!> 'active' cycle. Development is done mostly here. \n
!
!> \authors Johannes Emmerig, DLR Oberpfaffenhofen, 2019,
!>          (Johannes.Emmerig@dlr.de)
!>          - Development of 'H2OEMIS'
!
!> \authors Patrick Jöckel, DLR Oberpfaffenhofen, 2019,
!>          (Patrick.Joeckel@dlr.de)
!>          - Mentoring the development
!
!> \version 1.0
!>     - First working SMIL version of 'H2OEMIS'
!
!> \todo
!>     - Implement h2oemis_init_tracer
!
! ====================================================================


! ********************************************************************

MODULE messy_h2oemis_si

#ifdef ECHAM5

  ! BMIL
  USE messy_main_blather_bi,    ONLY: start_message_bi, end_message_bi, &
                                      error_bi, info_bi, warning_bi

  ! SMCL
  USE messy_h2oemis

#ifdef MESSYTENDENCY
  USE messy_main_tendency_bi,   ONLY: mtend_get_handle,                 &
                                      mtend_add_l,                      &
                                      mtend_register,                   &
                                      mtend_id_q
#endif

  ! STANDARDS
  IMPLICIT NONE
  INTRINSIC :: NULL
  PRIVATE

  !> @var Initialize pointer for external emissions from import grid
  REAL(DP), DIMENSION(:,:,:), POINTER :: emis_ptr => NULL()

#ifdef MESSYTENDENCY
  INTEGER :: my_handle
#endif

  ! from bmil/messy_main_control_bi.f90
  ! --------------------------------------------------------------------
  PUBLIC :: h2oemis_initialize    ! Initialize submodel
  PUBLIC :: h2oemis_init_coupling ! Set pointers for coupling to BM and other SMs
  PUBLIC :: h2oemis_physc         ! Entry point in time loop (current vector)
  ! --------------------------------------------------------------------

CONTAINS

  !> \brief Initialize submodel 'H2OEMIS'
  !> \details Call handle and register handle for specific humidity.
  !> This is a MESSYTENDENCY procedure.
  ! ====================================================================
  SUBROUTINE h2oemis_initialize
    
    ! BMIL
    USE messy_main_mpi_bi,    ONLY: p_parallel_io, p_io, p_bcast
    USE messy_main_tools,     ONLY: find_next_free_unit

    IMPLICIT NONE

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'h2oemis_initialize'
    INTEGER                     :: status ! error status
    INTEGER                     :: iou    ! I/O unit

#ifdef MESSYTENDENCY
    my_handle = mtend_get_handle(modstr)
    CALL mtend_register(my_handle,mtend_id_q)
#endif

  END SUBROUTINE h2oemis_initialize
  ! ====================================================================

  !> \brief Initialize coupling of 'H2OEMIS'
  !> \details Define external emissions (channel object from e.g. 
  !>  import grid)
  !> \param EMIS_H2O Pointer "emis_ptr" refers to channel.nml object "AIR_H2O"
  ! ====================================================================
  SUBROUTINE h2oemis_init_coupling

    ! BMIL
    USE messy_main_channel_error_bi, ONLY: channel_halt
    USE messy_main_channel,          ONLY: get_channel_object
    
    ! STANDARD
    IMPLICIT NONE

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr   = 'h2oemis_init_coupling'
    INTEGER                     :: status

    ! PARAMETER 
    CHARACTER(LEN=*), PARAMETER :: EMIS_H2O = 'AIR_H2O'

    CALL start_message_bi(modstr,'COUPLING INITIALIZATION',substr)  ! log-output

    ! SET POINTERS TO CHANNEL OBJECT
    CALL get_channel_object(status, 'import_grid', EMIS_H2O, p3=emis_ptr)
    CALL channel_halt(substr, status)

    CALL end_message_bi(modstr,'COUPLING INITIALIZATIO INITIALIZATION',substr)  ! log-output

  END SUBROUTINE h2oemis_init_coupling
  ! ====================================================================

  !> \brief Do calculations for 'H2OEMIS' 
  !> \details The subroutine 'physc' is called within the time loop. It 
  !> constitutes the main entry point for additional processes or
  !> diagnostics. Only the current vector of the grid-point-fields
  !> is accessible. \n
  !> The emission input is converted from # molecules/m3 to kg/m3. With
  !> recalculation of mixing ratio to specific humidity, the emissions 
  !> are added.
  !> @todo Decide on one calculation
  !> @test Test both possible calculations on orders of magnitude
  !> \param mconv Constants for conversion of number of molecules to kg
  !> @return qte_emis(:,:)  Tendency of spec. humidity with additional part of emission input
  ! ====================================================================
  SUBROUTINE h2oemis_physc

    ! MESSy BASEMODEL INTERFACE LAYER (BMIL)
    USE messy_main_constants_mem,    ONLY: N_A, R_gas, M_air, M_H2O
    USE messy_main_timer,            ONLY: time_step_len
    USE messy_main_data_bi,          ONLY: press_3d,              &
                                           qte_3d,                &
                                           qm1_3d,                &
                                           tm1_3d
    USE messy_main_grid_def_mem_bi, ONLY:  kproma,                &
                                           jrow,                  &
                                           nlev
    ! STANDARDS
    IMPLICIT NONE

    ! LOCAL PARAMETERS
    CHARACTER(LEN=*), PARAMETER :: substr = 'h2oemis_physc' 
    INTEGER                     :: status 

    ! CONSTANTS FOR MOLAR CONVERSION
    REAL(DP), PARAMETER :: mconv   = M_H2O/N_A     ! Conversion of # molec to kg 
    REAL(DP), PARAMETER :: molrat  = M_air/M_H2O   ! Molar ratio of air to water

    ! LOCAL VARIABLES FOR CALCULATING SPECIFIC HUMIDITY
    REAL(DP),DIMENSION(kproma,nlev) :: dryairdens  ! Density of dry air
    REAL(DP),DIMENSION(Kproma,nlev) :: mmr_emis    ! Molar mixing ratio of input emission
    REAL(DP),DIMENSION(kproma,nlev) :: qte_emis    ! Specific humidity tendency
    REAL(DP),DIMENSION(kproma,nlev) :: h2o_in      ! Converted emis_ptr in [kg.m-3]
    REAL(DP),DIMENSION(kproma,nlev) :: h2ote       ! H2O mixing ratio tendency calculated from spec. hum.
    REAL(DP),DIMENSION(kproma,nlev) :: h2o         ! H2O mixing ratio calculated from global spec. hum.

    ! CONVERSION OF INPUT TO MIX. RATIO
    ! -------------------------------------------------- 

    ! MOLEC/VOL TO KG/VOL AND 1/S TO 1
    h2o_in(:,:) = emis_ptr(1:kproma,:,jrow) &     ! Function of emis_ptr conversion 
                   * mconv * time_step_len   

    ! CALC: DENSITY OF DRY AIR IN KG/VOL 
    dryairdens(:,:) = press_3d(1:kproma,:,jrow) & 
                      * M_air &     
                      / ( R_gas * tm1_3d(1:kproma,:,jrow) )

    ! CALC: MOLAR MIX. RATIO OF H2O-EMISSION
    mmr_emis(:,:) = ( h2o_in(:,:) &               ! is h2ote(t=dt)
                    / dryairdens(:,:) ) &           
                    * molrat                      ! [n/n] = mol(H2O)/mol(dryair)
    ! -------------------------------------------------- 

    ! ADDITION OF H2O EMISSION TO SPEC. HUMIDITY
    ! -------------------------------------------------- 

    ! SAVE LOCAL COPY OF QTE
    !qte_local(:,:) = qte_3d(1:kproma,:,jrow)

    ! -------------------------------------------------- 

  !!$  ! CALC. q(t) = qm1_3d(t=0) + qte(t=dt)*dt
  !!$  q_3d(1:kproma,:,jrow) = qm1_3d(1:kproma,:,jrow) &
  !!$                          + qte_3d(1:kproma,:,jrow) &
  !!$                          * time_step_len
  !!$
  !!$  ! PREPARE H2O MIX. RATIO
  !!$  h2o(:,:) = ( qm1_3d(1:kproma,:,jrow) &
  !!$             / ( 1.0_dp - qm1_3d(1:kproma,:,jrow) ) ) &
  !!$             * molrat
  !!$
  !!$  h2o(:,:) = h2o(:,:) + mmr_emis(:,:)
  !!$
  !!$  ! PREPARE H2Ote(t=0) FROM H2Ote = H2Ote(t=0) + H2Ote(t=dt)
  !!$  h2ote(:,:) = ( qte_3d(1:kproma,:,jrow) &                     
  !!$               / ( ( 1.0_dp - q_3d(1:kproma,:,jrow) ) & 
  !!$               * ( 1.0_dp - qm1_3d(1:kproma,:,jrow) ) ) ) &
  !!$               * molrat
  !!$
  !!$  ! ADD H2Ote(t=dt) TO H2Ote(t=0)
  !!$  h2ote(:,:) = h2ote(:,:) &
  !!$               + ( mmr_emis(:,:) / time_step_len )
  !!$
  !!$  ! FINALIZE: MOLAR MIX RATIO TO SPEC. HUMIDITY 
  !!$  qte_emis(:,:) = ( h2ote(:,:) &                     
  !!$                  / ( ( molrat + h2o(:,:) ) &
  !!$                  * ( molrat + h2o(:,:) + h2ote(:,:) * time_step_len ) ) ) &
  !!$                  * molrat
  !!$
    ! -------------------------------------------------- 
 
    ! PREPARE H2O MIX. RATIO
    h2o(:,:) = ( qm1_3d(1:kproma,:,jrow) &
               / ( 1.0_dp - qm1_3d(1:kproma,:,jrow) ) ) &
               * molrat
 
    h2o(:,:) = h2o(:,:) + mmr_emis(:,:)
 
    ! PREPARE H2Ote(t=0) FROM H2Ote = H2Ote(t=0) + H2Ote(t=dt)
    h2ote(:,:) = ( qte_3d(1:kproma,:,jrow) &     
                 / ( 1.0_dp - qm1_3d(1:kproma,:,jrow) )**2 ) &
                 * molrat
 
    ! ADD H2Ote(t=dt) TO H2Ote(t=0)
    h2ote(:,:) = h2ote(:,:) &
                 + ( mmr_emis(:,:) / time_step_len )
 
    ! FINALIZE: MOLAR MIX RATIO TO SPEC. HUMIDITY 
    qte_emis(:,:) = ( h2ote(:,:) &                     
                    / ( molrat + h2o(:,:) )**2 ) &
                    * molrat
 
    ! -------------------------------------------------- 

    ! OUTPUT: QTE OF H2O EMISSION VIA TENDENCY OR QTE_3D
    ! -------------------------------------------------- 
#ifdef MESSYTENDENCY
    qte_emis(:,:) = qte_emis(:,:) - qte_3d(1:kproma,:,jrow)
    CALL mtend_add_l(my_handle, mtend_id_q, px=qte_emis)
#else
    qte_3d(1:kproma,:,jrow) = qte_emis(:,:)
#endif
    ! -------------------------------------------------- 


  END SUBROUTINE h2oemis_physc
  ! ====================================================================

#endif

END MODULE messy_h2oemis_si
! **********************************************************************

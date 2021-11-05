MODULE messy_hamocc_e5

  !
  !  AUTHOR:  Pozzer Andrea, MPICH, May 2007
  !
  !  mz_bk_20120628 code re-formated
  !
  !  MESSy Interface for the Submodel kernel
  !
  !  UNIT SYSTEM : 
  !
  !   [M]  = [mol/L] = [kmol/m^3]
  !   [nM] = [ 10^-9 mol/L] = [10^-9 mol/dm^3]
  !        = [10^-6 mol/m^3] = [10^-9 kmol/m^3] 
  !   [M]  = [mol/dm^3] = [mol/L] =  [kmol/m^3]

  ! kernel
  USE messy_hamocc

  !MESSY
  USE messy_main_timer_bi,      ONLY: p_bcast_event, timer_event_init
  USE messy_main_timer_event,   ONLY: io_time_event, TRIG_FIRST               &
       &                            , time_event
  USE messy_main_mpi_bi,        ONLY: message, finish
  USE messy_main_timer,         ONLY: YEAR, MONTH, DAY, HOUR, MINUTE, SECOND  &
       &                            , lstart

  IMPLICIT NONE 
  INTRINSIC TRIM
  PRIVATE

  ! mz_bk_20120613 make all module variables SAVE!!
  SAVE

  LOGICAL :: LRUNSM = .TRUE.  ! RUN THIS SUBMODEL

  ! NEW REPRESENTATION(S)
  INTEGER :: GP_3D_SED
  INTEGER :: DIMID_DEPTH_SED

  ! NEEDED FOR TIME LOOP
  LOGICAL ::  l_trig_hamocc = .TRUE.  
  TYPE(io_time_event), PUBLIC :: trig_hamocc 
  TYPE(time_event),    PUBLIC :: ev_trig_hamocc
  LOGICAL ::  l_trig_hamocc_day = .TRUE.  
  TYPE(io_time_event), PUBLIC :: trig_hamocc_day 
  TYPE(time_event),    PUBLIC :: ev_trig_hamocc_day

  ! CPL-NAMELIST PARAMETER
  LOGICAL, PUBLIC :: L_DUST
  ! information which channel and channel object to use for "online" dust input
  CHARACTER(LEN=STRLEN_MEDIUM),PUBLIC :: dust_input(2) = ''

  ! mz_bk_20100519+
  ! skip integration for passive tracer tests in the ocean
  LOGICAL, PUBLIC :: L_SKIP_INT
  ! mz_bk_20100519-

  ! CLIMATOLOGICAL DUST INPUT : dustin
  ! in kg/m2/year (0 if ddpo(i,j,1).gt.0.5)
  REAL(DP), DIMENSION (:,:,:), POINTER :: dustin

  ! TIME STEP
  INTEGER :: DT_steps
  
  REAL(dp), DIMENSION(:,:,:), POINTER :: sao   => NULL() 
  REAL(dp), DIMENSION(:,:,:), POINTER :: tho   => NULL() 
  REAL(dp), DIMENSION(:,:),   POINTER :: fswr_in  => NULL() 

  ! Indices for oceanic tracers
  INTEGER :: idx_sco212   &
       &   , idx_alkali   &
       &   , idx_phosph   &
       &   , idx_oxygen   &
       &   , idx_gasnit   &
       &   , idx_ano3     &
       &   , idx_silica   &
       &   , idx_doc      &
       &   , idx_phy      &
       &   , idx_zoo      &
       &   , idx_an2o     &
       &   , idx_dms      &
       &   , idx_iron     &
       &   , idx_fdust    
  INTEGER :: idx_hp       & 
       &   , idx_co3mm
  INTEGER :: idx_co
  INTEGER :: idx_isop
  INTEGER :: idx_ch3i
  INTEGER :: idx_det      & 
       &   , idx_calc     &   
       &   , idx_opal   
  ! In case of L_AGG
  INTEGER :: idx_nos      &
       &   , idx_adust
   

  ! mz_bk_20120613+ not used
  ! INTEGER :: ii,jj,kk
  ! mz_bk_20120613-
  LOGICAL :: l_init_bgc
  
  ! PUBLIC ECHAM-5 INTERFACE ROUTINES TO 'messy_SUBMODELS'
  PUBLIC :: hamocc_initialize
  PUBLIC :: hamocc_new_tracer
  PUBLIC :: hamocc_init_memory
  PUBLIC :: hamocc_init_coupling
  PUBLIC :: hamocc_global_start
  PUBLIC :: hamocc_global_end
  PUBLIC :: hamocc_free_memory

CONTAINS

  ! ---------------------------------------------------------------------------
  SUBROUTINE hamocc_initialize
    !
    ! HAMOCC MODULE ROUTINE (ECHAM-5 INTERFACE)
    !
    ! INITIALIZATION OF GLOBAL VARIABLES FROM NAMELIST
    ! IN PARALLEL ENVIRONMENT
    !
    ! Author: Pozzer Andrea, MPICH, May 2007
    ! ECHAM5/MESSy

    USE messy_main_blather_bi,       ONLY: start_message_bi, end_message_bi    &
         &                            , info_bi, warning_bi, error_bi
    USE messy_main_mpi_bi,           ONLY: p_parallel_io, p_io, p_pe, p_bcast
    USE messy_main_timer,            ONLY: delta_time        
    USE messy_main_tools,            ONLY: find_next_free_unit
    USE messy_main_channel_error_bi, ONLY: channel_halt
    USE messy_main_channel_bi,       ONLY: DC_GP_MPIOM                        &
            &                            , DIMID_LON_MPIOM, DIMID_LAT_MPIOM   &
            &                            , GP_3D_MPIOM, REPR_UNDEF
    USE messy_mpiom_mem_e5,          ONLY: DT
    USE messy_main_channel_repr,     ONLY: new_representation, AUTO
    USE messy_main_channel_dimensions, ONLY: new_dimension                    &
         &                                 , add_dimension_variable_att       &
         &                                 , add_dimension_variable
    ! TODO : to be eliminated --> through interface!
    USE messy_mpiom_mem_e5,          ONLY: ie ,je

    IMPLICIT NONE

    INTEGER i,j,k,l
    CHARACTER(LEN=*), PARAMETER :: substr='hamocc_initialize'
    INTEGER                     :: status
    INTEGER                     :: iou    ! I/O unit

    ! helper array for sediment depths
    REAL(dp), DIMENSION(:), ALLOCATABLE :: array
   
    CALL start_message_bi(modstr, 'INITIALIZE', substr)

    !-----------------------------------------------------------------------
    ! CTRL / CPL NAMELIST
    !-----------------------------------------------------------------------

    ! INITIALIZE CTRL
    IF (p_parallel_io) THEN
       iou = find_next_free_unit(100,200)
       ! *** CALL CORE ROUTINE:
       CALL hamocc_read_nml_ctrl(status, iou)
       IF (status /= 0) CALL finish(substr)
    END IF

    ! BROADCAST RESULTS
    CALL p_bcast(L_CHECK, p_io)   
    CALL p_bcast(L_AGG,   p_io)   
    ! mz_bk_20120606+
    CALL p_bcast(L_BGC_DIAG,   p_io)   
    ! mz_bk_20120606-
    CALL p_bcast(isac,    p_io)
    CALL p_bcast(rmasks,  p_io)
    CALL p_bcast(hamocc_ini_files_path, p_io)
    ! mz_bk_20120724+
    CALL p_bcast(ro2ut,   p_io)
    CALL p_bcast(rcar,    p_io)
    CALL p_bcast(rnit,    p_io)
    CALL p_bcast(nitdem,  p_io)
    CALL p_bcast(n2prod,  p_io)
    CALL p_bcast(rcalc,   p_io)
    CALL p_bcast(ropal,   p_io)
    CALL p_bcast(n2_fixation, p_io)
    CALL p_bcast(rno3,    p_io)
    CALL p_bcast(perc_diron, p_io)
    CALL p_bcast(riron,   p_io)
    CALL p_bcast(fesoly,  p_io)
    CALL p_bcast(relaxfe, p_io)
    CALL p_bcast(pi_alpha, p_io)
    CALL p_bcast(fPAR,    p_io)
#if defined MPIOM_13B
    CALL p_bcast(ctochl,  p_io)
    CALL p_bcast(atten_w, p_io)
    CALL p_bcast(atten_f, p_io)
#endif
    CALL p_bcast(phytomi, p_io)
    CALL p_bcast(bkphy,   p_io)
    CALL p_bcast(bkopal,  p_io)
    CALL p_bcast(remido,  p_io)
    CALL p_bcast(dyphy,   p_io)
    CALL p_bcast(gammap,  p_io)
    CALL p_bcast(bkzoo,   p_io)
    CALL p_bcast(grami,   p_io)
    CALL p_bcast(zinges,  p_io)
    CALL p_bcast(epsher,  p_io)
    CALL p_bcast(grazra,  p_io)
    CALL p_bcast(spemor,  p_io)
    CALL p_bcast(gammaz,  p_io)
    CALL p_bcast(ecan,    p_io)
    CALL p_bcast(sinkspeed_poc, p_io)
    CALL p_bcast(sinkspeed_opal, p_io)
    CALL p_bcast(sinkspeed_cal, p_io)
    CALL p_bcast(drempoc, p_io)
    CALL p_bcast(dremdoc, p_io)
    CALL p_bcast(dremn2o, p_io)
    CALL p_bcast(denitrification, p_io)
    CALL p_bcast(sulfate_reduction, p_io)
    CALL p_bcast(dremopal, p_io)
    CALL p_bcast(dremcalc, p_io)
    CALL p_bcast(dphymor, p_io)
    CALL p_bcast(dzoomor, p_io)
    DO i=1,6
       CALL p_bcast(dmspar(i),  p_io)
    ENDDO
    CALL p_bcast(calmax, p_io)

    ! mz_bk_20120804+
    ! moved to hamocc_init_memory, because dtb is not set here.
    ! dtb will be set up in hamocc_init_memory
    ! remido_dtb   = remido   * dtb
    ! dyphy_dtb    = dyphy    * dtb
    ! grazra_dtb   = grazra   * dtb
    ! spemor_dtb   = spemor   * dtb
    ! gammap_dtb   = gammap   * dtb
    ! gammaz_dtb   = gammaz   * dtb
    ! drempoc_dtb  = drempoc  * dtb
    ! dremdoc_dtb  = dremdoc  * dtb
    ! dremopal_dtb = dremopal * dtb
    ! dremcalc_dtb = dremcalc * dtb
    ! dremn2o_dtb  = dremn2o  * dtb
    ! dphymor_dtb  = dphymor  * dtb
    ! dzoomor_dtb  = dzoomor  * dtb
    ! relaxfe_dtb  = relaxfe  * dtb
    ! mz_bk_20120804-
    ! mz_bk_20120724-


    ! INITIALIZE CPL
    IF (p_parallel_io) THEN
       iou = find_next_free_unit(100,200)
       ! *** CALL CORE ROUTINE:
       CALL hamocc_read_nml_cpl(status, iou)
       IF (status /= 0) CALL finish(substr)
    END IF

    ! BROADCAST RESULTS
    CALL p_bcast(L_DUST,     p_io)   
    CALL p_bcast(dust_input(1), p_io)   
    CALL p_bcast(dust_input(2), p_io)
    CALL p_bcast(L_SKIP_INT, p_io) ! mz_bk_20100519

    !-----------------------------------------------------------------------
    ! CTRL NAMELIST
    !-----------------------------------------------------------------------
   
    ! mz_bk_20120610 some more output...
    IF (p_parallel_io) THEN
       WRITE(*,*) " --- HAMOCC CTRL namelist --------------------"
       WRITE(*,*) " rmasks = ", rmasks
       WRITE(*,*) " isac (acceleration factor for sediment) = ", isac
       WRITE(*,*) " hamocc checks = ", L_CHECK
       WRITE(*,*) " BGC diagnostics output = ", L_BGC_DIAG
       WRITE(*,*) " --- HAMOCC CPL namelist ---------------------"
       WRITE(*,*) " climatological dust input = ", L_DUST
       IF (.NOT. L_DUST) THEN 
          ! mz_bk_20120614+
          ! Using the following two lines, there is a "deadlock" with some
          ! setup(s?) (E5MPIOMTROP/AO_SPINUP) on blizzard and an
          ! executable compiled with compiler xlf v.13.1.0.8, optimized -O3,
          ! without error checks.
          ! WRITE(*,*) "  dust input from channel object ",TRIM(dust_input(1)) &
          !      &    ," ", TRIM(dust_input(2))
          ! mz_bk_20120614-
          WRITE(*,*) "  dust input from channel object "                      &
               &     //TRIM(dust_input(1))//"%"//TRIM(dust_input(2))
       ENDIF
       WRITE(*,*) " hamocc init files in: "//TRIM(hamocc_ini_files_path)
       WRITE(*,*) " ---------------------------------------------"
    ENDIF

    !-----------------------------------------------------------------------
    !  CHECK MPIOM MODEL PRESENCE           
    !-----------------------------------------------------------------------

    LRUNSM = (GP_3D_MPIOM /= REPR_UNDEF)
    IF (.NOT. LRUNSM) THEN
       IF (p_parallel_io) THEN
          WRITE(*,*) 'MPIOM REPRESENTATION GP_3D_MPIOM IS UNDEFINED.'
          WRITE(*,*) 'SUBMODEL CANNOT BE RUN.'
       END IF
       CALL end_message_bi(modstr, 'INITIALISATION', substr)
       CALL finish(substr)
    END IF

    !-----------------------------------------------------------------------
    !  INITIALIZATION           
    !-----------------------------------------------------------------------

    DT_steps=INT(DT/delta_time)

    write(*,*) "DT_steps = ", DT_steps

    trig_hamocc = io_time_event (DT_steps, 'steps',TRIG_FIRST, 0)
    trig_hamocc_day = io_time_event (1, 'days',TRIG_FIRST, 0)

    CALL p_bcast_event (trig_hamocc, p_io)
    CALL p_bcast_event (trig_hamocc_day, p_io)

    CALL timer_event_init(ev_trig_hamocc, trig_hamocc,                        &
         &                'hamocc computation', 'present')
    CALL timer_event_init(ev_trig_hamocc_day, trig_hamocc_day,                &
         &                'hamocc computation', 'present')

    WRITE (*,*) 'trig_hamocc: '      ,trig_hamocc
    WRITE (*,*) 'trig_hamocc_day: '  ,trig_hamocc_day

    !-----------------------------------------------------------------------
    !  NEW REPRESENTATION(S)
    !-----------------------------------------------------------------------
    ! NEW DIMENSIONS
    !-----------------LEVELS DEPTHS SEDIMENT MODEL----------------------
    ALLOCATE(dzs(ksp))
    ALLOCATE(array(ks))

    ! define sediment layer thickness [m]
    dzs(1)  = 0.001_dp
    dzs(2)  = 0.003_dp
    dzs(3)  = 0.005_dp
    dzs(4)  = 0.007_dp
    dzs(5)  = 0.009_dp
    dzs(6)  = 0.011_dp
    dzs(7)  = 0.013_dp
    dzs(8)  = 0.015_dp
    dzs(9)  = 0.017_dp
    dzs(10) = 0.019_dp
    dzs(11) = 0.021_dp
    dzs(12) = 0.023_dp
    dzs(13) = 0.025_dp

    IF (p_pe == p_io) THEN
       WRITE(*,*)  ' '
       WRITE(*,*)  'Sediment layer thickness [m] : '
       WRITE(*,'(5F9.3)') dzs
       WRITE(*,*)  ' '
    ENDIF

    array(1) = (dzs(1) / 2._dp)
    DO i=2,ks
       array(i) = array(i-1) + (dzs(i) + dzs(i+1))/ 2._dp
    ENDDO

    CALL new_dimension(status, DIMID_DEPTH_SED, 'sediment_depth', ks)
    CALL channel_halt(substr, status)

    CALL add_dimension_variable(status, 'sediment_depth'                      &
         &                    , 'sediment_depth',array )
    CALL channel_halt(substr, status)

    CALL add_dimension_variable_att(status, 'sediment_depth'                  &
         &                        , 'sediment_depth', 'long_name'             &
         &                        , c='ocean sediment_depth ')
    CALL channel_halt(substr, status)
    !
    CALL add_dimension_variable_att(status, 'sediment_depth'                  &
         , 'sediment_depth', 'units', c='m')
    CALL channel_halt(substr, status)

    DEALLOCATE(array)

    ! NEW REPRESENTATION
    !-----------------SEDIMENT MODEL-----------------------------------

    CALL new_representation(status, GP_3D_SED, 'GP_3D_SED'                    &
         &                , rank = 3, link = 'xxx-', dctype = DC_GP_MPIOM     &
         &                , dimension_ids = (/ DIMID_LON_MPIOM                &
         &                              , DIMID_LAT_MPIOM, DIMID_DEPTH_SED /) &
         &                , ldimlen       = (/ ie , je, AUTO /)               &
         &                , output_order  = (/ 1,2,3 /)                       &
         &                , axis = 'XYZ-')


  END SUBROUTINE hamocc_initialize
  ! ---------------------------------------------------------------------------

  ! ---------------------------------------------------------------------------
  SUBROUTINE hamocc_new_tracer

    USE messy_main_tracer, ONLY : new_tracer, OCEAN, CONCENTRATION,ON,OFF     &
         &                      , MAX_CASK_I, MAX_CASK_R,I_ADVECT,I_CONVECT   &
         &                      , I_VDIFF,I_DRYDEP,I_SEDI,I_SCAV              &
         &                      , I_MIX,I_FORCE_COL,I_INTEGRATE,I_TIMEFILTER  &
         &                      , I_FORCE_INIT, R_MOLARMASS
    USE messy_main_tracer_tools_bi, ONLY : tracer_halt

    INTEGER   :: status
    INTEGER, DIMENSION(MAX_CASK_I) :: HAMOCC_CASK_I 
    REAL(DP),    DIMENSION(MAX_CASK_R) :: HAMOCC_CASK_R 

    CHARACTER(LEN=*), PARAMETER :: substr='hamocc_new_tracer'


    HAMOCC_CASK_I(:) = OFF

    HAMOCC_CASK_I(I_ADVECT    ) = ON   ! ADVECTION
    HAMOCC_CASK_I(I_CONVECT   ) = ON   ! CONVECTION
    HAMOCC_CASK_I(I_VDIFF     ) = ON   ! VERTICAL DIFFUSION
    HAMOCC_CASK_I(I_DRYDEP    ) = OFF  ! DRY DEPOSITION
    HAMOCC_CASK_I(I_SEDI      ) = OFF  ! SEDIMENTATION
    HAMOCC_CASK_I(I_SCAV      ) = OFF  ! SCAVENGING
    HAMOCC_CASK_I(I_MIX       ) = ON   ! TURBULENT MIXING
    HAMOCC_CASK_I(I_FORCE_COL ) = OFF  ! FORCING IN COLUMN MODE
    HAMOCC_CASK_I(I_INTEGRATE ) = ON   ! TIME INTEGRATION
    HAMOCC_CASK_I(I_TIMEFILTER) = ON   ! TIME FILTER
    HAMOCC_CASK_I(I_FORCE_INIT) = OFF  ! FORCE INIT AFTER RESTART

    HAMOCC_CASK_R(:) = 0.0_dp


    !TODO: WE NEED FLUXES OF: CO2, O2, N2,N2O
    !--------------------------------------!
    !              TRACERS                 !
    !--------------------------------------!
    !-----------------------------------------------------------!
    ! N.B. For the time being the molar mass is                 !
    !      NOT always of the tracer.                            ! 
    ! As example for N2O we have the molar mass of N; hence     !
    ! the calculation of pdef are in [Kg N] and not in [Kg]     !
    !-----------------------------------------------------------!

    ! ----------- DMS --------------- !
    !  HAMOCC_CASK_R(R_PSS  )      = 0.0 ! we do not need it!
    HAMOCC_CASK_R(R_MOLARMASS)  = 62.13
    HAMOCC_CASK_I(I_ADVECT    ) = ON  ! ADVECTION
    
    CALL new_tracer(status, setname ='om'                                     &
         &        , basename='DMS'                                            &
         &        , idx = idx_dms                                             &
         &        , submodel=modstr                                           &
         &        , longname='dimethyl sulfide'                               & 
         &        , unit='kmol/m^3'                                           &
         &        , medium = OCEAN, quantity = CONCENTRATION                  &
         &        , cask_i=HAMOCC_CASK_I, cask_r=HAMOCC_CASK_R ) 
    CALL tracer_halt(substr, status)

    ! ----------- DUST --------------- !
    HAMOCC_CASK_R(R_MOLARMASS)  = 1 !??????? 
    HAMOCC_CASK_I(I_ADVECT    ) = ON  ! ADVECTION

    CALL new_tracer(status, setname ='om'                                     &
         &        , basename='DUST'                                           &
         &        , idx = idx_fdust                                           &
         &        , submodel=modstr                                           &
         &        , longname='free (non aggregate) dust'                      &
         &        , unit='kmol/m^3'                                           &
         &        , medium = OCEAN, quantity = CONCENTRATION                  &
         &        , cask_i=HAMOCC_CASK_I, cask_r=HAMOCC_CASK_R ) 
    CALL tracer_halt(substr, status)

    ! ----------- FE ---------------- !
    HAMOCC_CASK_R(R_MOLARMASS)  = 55.845
    HAMOCC_CASK_I(I_ADVECT    ) = ON  ! ADVECTION

    CALL new_tracer(status, setname ='om'                                     &
         &        , basename='Fe'                                             &
         &        , idx = idx_iron                                            &
         &        , submodel=modstr                                           &
         &        , longname='Dissolved Iron'                                 &
         &        , unit='kmol/m^3'                                           &
         &        , medium = OCEAN, quantity = CONCENTRATION                  &
         &        , cask_i=HAMOCC_CASK_I, cask_r=HAMOCC_CASK_R )
    CALL tracer_halt(substr, status)

    ! ----------- PO4 --------------- !
    HAMOCC_CASK_R(R_MOLARMASS)  = 30.9737 !(P MW)
    HAMOCC_CASK_I(I_ADVECT    ) = ON  ! ADVECTION

    CALL new_tracer(status, setname ='om'                                     &
         &        , basename='PO4'                                            &
         &        , idx = idx_phosph                                          &
         &        , submodel=modstr                                           &
         &        , longname='Phosphate'                                      &
         &        , unit='kmol/m^3'                                           &
         &        , medium = OCEAN, quantity = CONCENTRATION                  &
         &        , cask_i=HAMOCC_CASK_I, cask_r=HAMOCC_CASK_R )
    CALL tracer_halt(substr, status)

    ! ----------- N2O --------------- !
    HAMOCC_CASK_R(R_MOLARMASS)  = 44.013
    HAMOCC_CASK_I(I_ADVECT    ) = ON  ! ADVECTION

    CALL new_tracer(status, setname ='om'                                     &
         &        , basename='N2O'                                            &
         &        , idx = idx_an2o                                            &
         &        , submodel=modstr                                           &
         &        , longname='nitrous oxide'                                  &
         &        , unit='kmol/m^3'                                           &
         &        , medium = OCEAN, quantity = CONCENTRATION                  &
         &        , cask_i=HAMOCC_CASK_I, cask_r=HAMOCC_CASK_R )
    CALL tracer_halt(substr, status)

    ! ----------- NO3 --------------- !
    HAMOCC_CASK_R(R_MOLARMASS)  = 14.0067
    HAMOCC_CASK_I(I_ADVECT    ) = ON  ! ADVECTION

    CALL new_tracer(status, setname ='om'                                     &
         &        , basename='NO3'                                            &
         &        , idx = idx_ano3                                            &
         &        , submodel=modstr                                           &
         &        , longname='Nitrate'                                        &
         &        , unit='kmol/m^3'                                           &
         &        , medium = OCEAN, quantity = CONCENTRATION                  &
         &        , cask_i=HAMOCC_CASK_I, cask_r=HAMOCC_CASK_R )
    CALL tracer_halt(substr, status)

    ! ----------- N2 --------------- !
    HAMOCC_CASK_R(R_MOLARMASS)  = 28.0134
    HAMOCC_CASK_I(I_ADVECT    ) = ON  ! ADVECTION

    CALL new_tracer(status, setname ='om'                                     &
         &        , basename='N2'                                             &
         &        , idx = idx_gasnit                                          &
         &        , submodel=modstr                                           &
         &        , longname='Dinitrogen'                                     &
         &        , unit='kmol/m^3'                                           &
         &        , medium = OCEAN, quantity = CONCENTRATION                  &
         &        , cask_i=HAMOCC_CASK_I, cask_r=HAMOCC_CASK_R )
    CALL tracer_halt(substr, status)

    ! ------------ ZOO ---------------- !             
    HAMOCC_CASK_R(R_MOLARMASS)  = 30.973762
    HAMOCC_CASK_I(I_ADVECT    ) = ON  ! ADVECTION

    CALL new_tracer(status, setname ='om'                                     &
         &        , basename='ZOO'                                            &
         &        , idx = idx_zoo                                             &
         &        , submodel=modstr                                           &
         &        , longname='zooplankton'                                    &
         &        , unit='kmol/m^3'                                           &
         &        , medium = OCEAN, quantity = CONCENTRATION                  &
         &        , cask_i=HAMOCC_CASK_I, cask_r=HAMOCC_CASK_R )
    CALL tracer_halt(substr, status)

    ! ------------- PHY --------------- !             
    HAMOCC_CASK_R(R_MOLARMASS)  = 30.973762
    HAMOCC_CASK_I(I_ADVECT    ) = ON  ! ADVECTION

    CALL new_tracer(status, setname ='om'                                     &
         &        , basename='PHY'                                            &
         &        , idx = idx_phy                                             &
         &        , submodel=modstr                                           &
         &        , longname='phytoplankton'                                  &
         &        , unit='kmol/m^3'                                           &
         &        , medium = OCEAN, quantity = CONCENTRATION                  &
         &        , cask_i=HAMOCC_CASK_I, cask_r=HAMOCC_CASK_R )
    CALL tracer_halt(substr, status)

    ! ---------- DOM ---------------- !             
    HAMOCC_CASK_R(R_MOLARMASS)  = 30.973762
    HAMOCC_CASK_I(I_ADVECT    ) = ON  ! ADVECTION

    CALL new_tracer(status, setname ='om'                                     &
         &        , basename='DOM'                                            &
         &        , idx = idx_doc                                             &
         &        , submodel=modstr                                           &
         &        , longname='Dissolved organic matter '                      & 
         &        , unit='kmol/m^3'                                           &
         &        , medium = OCEAN, quantity = CONCENTRATION                  &
         &        , cask_i=HAMOCC_CASK_I, cask_r=HAMOCC_CASK_R )
    CALL tracer_halt(substr, status)

    ! ---------- DET ---------------- !             
    HAMOCC_CASK_R(R_MOLARMASS)  = 30.973762
    ! non-advected (fast sinking) tracer
    ! mz_bk_20120628: in the newer hamocc version of mpiesm-1.0.00, detritus
    ! is advected. In case of an update, also change it here!!
    ! mz_bk_20120706+
    ! HAMOCC_CASK_I(I_ADVECT    ) = OFF  ! ADVECTION
    HAMOCC_CASK_I(I_ADVECT    ) = ON  ! ADVECTION
    ! mz_bk_20120706-

    CALL new_tracer(status, setname ='om'                                     &
         &        , basename='DET'                                            &
         &        , idx = idx_det                                             &
         &        , submodel=modstr                                           &
         &        , longname='detritus'                                       &
         &        , unit='kmol/m^3'                                           &
         &        , medium = OCEAN, quantity = CONCENTRATION                  &
         &        , cask_i=HAMOCC_CASK_I, cask_r=HAMOCC_CASK_R )
    CALL tracer_halt(substr, status)

    ! ----------- CaCO3 -------------- !             
    HAMOCC_CASK_R(R_MOLARMASS)  = 12.0107
    ! non-advected (fast sinking) tracer
    HAMOCC_CASK_I(I_ADVECT    ) = OFF  ! ADVECTION

    CALL new_tracer(status, setname ='om'                                     &
         &        , basename='CaCO3'                                          &
         &        , idx = idx_calc                                            &
         &        , submodel=modstr                                           &
         &        , longname='calcium carbonate shells'                       &
         &        , unit='kmol/m^3'                                           &
         &        , medium = OCEAN, quantity = CONCENTRATION                  &
         &        , cask_i=HAMOCC_CASK_I, cask_r=HAMOCC_CASK_R )
    CALL tracer_halt(substr, status)

    ! ------------ C12 --------------- !             
    HAMOCC_CASK_R(R_MOLARMASS)  = 44.0095
    HAMOCC_CASK_I(I_ADVECT    ) = ON  ! ADVECTION

    CALL new_tracer(status, setname ='om'                                     &
         &        , basename='C12'                                            &
         &        , idx = idx_sco212                                          &
         &        , submodel=modstr                                           &
         &        , longname='total dissolved inorganic'                      &
         &        , unit='kmol/m^3'                                           &
         &        , medium = OCEAN, quantity = CONCENTRATION                  &
         &        , cask_i=HAMOCC_CASK_I, cask_r=HAMOCC_CASK_R )
    CALL tracer_halt(substr, status)

    ! ----------- O2 ---------------- !
    HAMOCC_CASK_R(R_MOLARMASS)  = 31.9988
    HAMOCC_CASK_I(I_ADVECT    ) = ON  ! ADVECTION

    CALL new_tracer(status, setname ='om'                                     &
         &        , basename='O2'                                             &
         &        , idx = idx_oxygen                                          &
         &        , submodel=modstr                                           &
         &        , longname='Oxygen'                                         &
         &        , unit='kmol/m^3'                                           &
         &        , medium = OCEAN, quantity = CONCENTRATION                  &
         &        , cask_i=HAMOCC_CASK_I, cask_r=HAMOCC_CASK_R )
    CALL tracer_halt(substr, status)

    ! ----------- Si ---------------- !
    HAMOCC_CASK_R(R_MOLARMASS)  = 28.0855
    HAMOCC_CASK_I(I_ADVECT    ) = ON  ! ADVECTION

    CALL new_tracer(status, setname ='om'                                     &
         &        , basename='Si'                                             &
         &        , idx = idx_silica                                          &
         &        , submodel=modstr                                           &
         &        , longname='Silicate'                                       &
         &        , unit='kmol/m^3'                                           &
         &        , medium = OCEAN, quantity = CONCENTRATION                  &
         &        , cask_i=HAMOCC_CASK_I, cask_r=HAMOCC_CASK_R )
    CALL tracer_halt(substr, status)

    ! ----------- Opal ---------------- !
    HAMOCC_CASK_R(R_MOLARMASS)  = 28.0855
    ! non-advected (fast sinking) tracer
    HAMOCC_CASK_I(I_ADVECT    ) = OFF  ! ADVECTION

    CALL new_tracer(status, setname ='om'                                     &
         &        , basename='OPAL'                                           &
         &        , idx = idx_opal                                            &
         &        , submodel=modstr                                           &
         &        , longname='Opal Shell'                                     &
         &        , unit='kmol/m^3'                                           &
         &        , medium = OCEAN, quantity = CONCENTRATION                  &
         &        , cask_i=HAMOCC_CASK_I, cask_r=HAMOCC_CASK_R ) 
    CALL tracer_halt(substr, status)

    ! ----------- Total alkalinity ---------------- !
    HAMOCC_CASK_R(R_MOLARMASS)  = 1.0
    HAMOCC_CASK_I(I_ADVECT    ) = ON  ! ADVECTION

    CALL new_tracer(status, setname ='om'                                     &
         &        , basename='AT'                                             &
         &        , idx = idx_alkali                                          &
         &        , submodel=modstr                                           &
         &        , longname='Total alkalinity'                               &
         &        , unit='kmol/m^3'                                           &
         &        , medium = OCEAN, quantity = CONCENTRATION                  &
         &        , cask_i=HAMOCC_CASK_I, cask_r=HAMOCC_CASK_R ) 
    CALL tracer_halt(substr, status)

    IF (L_AGG) THEN

       ! ----------- Number of snow aggregates ---------------- !
       HAMOCC_CASK_R(R_MOLARMASS)  = 1
       HAMOCC_CASK_I(I_ADVECT    ) = ON  ! ADVECTION

       CALL new_tracer(status, setname ='om'                                  &
            &        , basename='NOS'                                         &
            &        , idx = idx_nos                                          &
            &        , submodel=modstr                                        &
            &        , longname='number of snow aggregates'                   &
            &        , unit='particles cm-3'                                  &
            &        , medium = OCEAN, quantity = CONCENTRATION               &
            &        , cask_i=HAMOCC_CASK_I, cask_r=HAMOCC_CASK_R )
       CALL tracer_halt(substr, status)

       ! ----------- Aggregate dust ---------------- !
       HAMOCC_CASK_R(R_MOLARMASS)  = 1
       HAMOCC_CASK_I(I_ADVECT    ) = ON  ! ADVECTION

       CALL new_tracer(status, setname ='om'                                  &
            &        , basename='ADUST'                                       &
            &        , idx = idx_adust                                        &
            &        , submodel=modstr                                        &
            &        , longname='aggregate dust'                              &
            &        , unit='Kg m-3'                                          &
            &        , medium = OCEAN, quantity = CONCENTRATION               &
            &        , cask_i=HAMOCC_CASK_I, cask_r=HAMOCC_CASK_R )
       CALL tracer_halt(substr, status)

    ENDIF ! L_AGG

    !-------DIAGNOSTIC TRACERS -------!

    ! ----------- H+ ---------------- !
    HAMOCC_CASK_R(R_MOLARMASS)  = 1.007
    HAMOCC_CASK_I(I_ADVECT    ) = OFF  ! ADVECTION
    HAMOCC_CASK_I(I_MIX       ) = OFF  ! MIXING
    HAMOCC_CASK_I(I_VDIFF     ) = OFF  ! DIFFUSION

    CALL new_tracer(status, setname ='om'                                     &
         &        , basename='Hp'                                             &
         &        , idx = idx_hp                                              &
         &        , submodel=modstr                                           &
         &        , longname='proton concentration'                           &
         &        , unit='kmol/m^3'                                           &
         &        , medium = OCEAN, quantity = CONCENTRATION                  &
         &        , cask_i=HAMOCC_CASK_I, cask_r=HAMOCC_CASK_R )
    CALL tracer_halt(substr, status)

    ! ----------- CO3-- ---------------- !
    HAMOCC_CASK_R(R_MOLARMASS)  = 58.00
    HAMOCC_CASK_I(I_ADVECT    ) = OFF  ! ADVECTION
    HAMOCC_CASK_I(I_MIX       ) = OFF  ! MIXING
    HAMOCC_CASK_I(I_VDIFF     ) = OFF  ! DIFFUSION

    CALL new_tracer(status, setname ='om'                                     &
         &        , basename='CO3mm'                                          &
         &        , idx = idx_co3mm                                           &
         &        , submodel=modstr                                           &
         &        , longname='CO3--'                                          &
         &        , unit='kmol/m^3'                                           &
         &        , medium = OCEAN, quantity = CONCENTRATION                  &
         &        , cask_i=HAMOCC_CASK_I, cask_r=HAMOCC_CASK_R )
    CALL tracer_halt(substr, status)

    ! ----------- CO ---------------- !
    HAMOCC_CASK_R(R_MOLARMASS)  = 28.0101
    HAMOCC_CASK_I(I_ADVECT    ) = ON  ! ADVECTION
    HAMOCC_CASK_I(I_MIX       ) = ON  ! MIXING
    HAMOCC_CASK_I(I_VDIFF     ) = ON  ! DIFFUSION

    CALL new_tracer(status, setname ='om'                                     &
         &        , basename='CO'                                             &
         &        , idx = idx_co                                              &
         &        , submodel=modstr                                           &
         &        , longname='Carbon monoxide'                                & 
         &        , unit='kmol/m^3'                                           &
         &        , medium = OCEAN, quantity = CONCENTRATION                  &
         &        , cask_i=HAMOCC_CASK_I, cask_r=HAMOCC_CASK_R )
    CALL tracer_halt(substr, status)

    ! ----------- ISOPRENE ---------------- !
    HAMOCC_CASK_R(R_MOLARMASS)  = 68.11
    HAMOCC_CASK_I(I_ADVECT    ) = OFF  ! ADVECTION
    HAMOCC_CASK_I(I_MIX       ) = OFF  ! MIXING
    HAMOCC_CASK_I(I_VDIFF     ) = OFF  ! DIFFUSION

    CALL new_tracer(status, setname ='om'                                     &
         &        , basename='ISOP'                                           &
         &        , idx = idx_isop                                            &
         &        , submodel=modstr                                           &
         &        , longname='Isoprene'                                       &
         &        , unit='kmol/m^3'                                           &
         &        , medium = OCEAN, quantity = CONCENTRATION                  &
         &        , cask_i=HAMOCC_CASK_I, cask_r=HAMOCC_CASK_R )
    CALL tracer_halt(substr, status)

    ! ----------- CH3I (Methil iodide) ---- !
    HAMOCC_CASK_R(R_MOLARMASS)  = 141.94
    HAMOCC_CASK_I(I_ADVECT    ) = ON  ! ADVECTION
    HAMOCC_CASK_I(I_MIX       ) = ON  ! MIXING
    HAMOCC_CASK_I(I_VDIFF     ) = ON  ! DIFFUSION

    CALL new_tracer(status, setname ='om'                                     &
         &        , basename='CH3I'                                           &
         &        , idx = idx_ch3i                                            &
         &        , submodel=modstr                                           &
         &        , longname='Methyliodide'                                   &
         &        , unit='kmol/m^3'                                           &
         &        , medium = OCEAN, quantity = CONCENTRATION                  &
         &        , cask_i=HAMOCC_CASK_I, cask_r=HAMOCC_CASK_R )
    CALL tracer_halt(substr, status)

  END SUBROUTINE hamocc_new_tracer
  ! ---------------------------------------------------------------------------

  ! ---------------------------------------------------------------------------
  SUBROUTINE hamocc_init_memory

    USE messy_main_channel_error_bi, ONLY: channel_halt
    USE messy_main_channel_bi,       ONLY: REPR_UNDEF                           &
         &                               , DC_GP_MPIOM, GP_3D_MPIOM, GP_2D_MPIOM
    USE messy_main_channel,          ONLY: new_channel, new_channel_object    &
         &                               , new_attribute                      &
         &                               , new_channel_object_reference
    USE messy_main_blather_bi,       ONLY: start_message_bi, end_message_bi
    USE messy_main_timer,            ONLY: delta_time, lstart, lresume 
    USE messy_main_tracer_mem_bi,    ONLY: xt_om
    ! TODO: eliminate!
    USE messy_mpiom_mem_e5,          ONLY: ddpo, ie, je, ke, ie_g, je_g, dt   &
         &                               , global_min, tiestw

    INTRINSIC NINT, TRIM

    ! CHANNEL MANAGEMENT
    INTEGER                      :: status
    CHARACTER(LEN=STRLEN_MEDIUM) :: reprname
    ! LOCAL
    CHARACTER(len=*), PARAMETER  :: substr = 'hamocc_init_memory'
    INTEGER i,j,k,l

    ! LOCAL
    REAL(DP) :: zmini


    CALL start_message_bi(modstr, 'MEMORY INITIALIZATION', substr)

    IF (GP_3D_MPIOM == REPR_UNDEF) THEN
       CALL message(substr, ' module MPIOM not running')  
       CALL message(substr, ' not possible to run HAMOCC')
       CALL finish(substr,'Run terminated')
    END IF

    !FOR OUTPUT
    ! (a) ALLOCATE POINTERS

    ! (b) new channel

    CALL new_channel(status, modstr, reprid=GP_3D_MPIOM,lrestreq=.TRUE.)
    CALL channel_halt(substr,status)

    ! (c) new channel objects

    CALL new_channel_object_reference(status, 'mpiom', 'longitude', modstr    &
         &                          , 'longitude', .TRUE.)
    CALL channel_halt(substr, status)

    CALL new_channel_object_reference(status, 'mpiom', 'latitude', modstr     &
         &                          , 'latitude', .TRUE.)
    CALL channel_halt(substr, status)

    !!----------------------------------------------------------------------
    !! THIS IS EQUAL TO INI_BGC.f90
    !!----------------------------------------------------------------------

    ! dtbgc = delta_time                       ! time step length [sec]
    dtbgc = dt                               ! time step length [sec]
    ! ndtdaybgc=NINT(86400./dtbgc)             ! time steps per day [no.]
    ! dtb=1./ndtdaybgc                         ! time step length [days]
    ! mz_bk_20120725+
    ! dtb = 1._dp / NINT(86400._dp / dtbgc)    ! time step length [days]
    dtb = 1._dp / NINT(86400._dp / dtbgc)  ! time step length [days]
    ! mz_bk_20120725-
    
    ! get first level, deeper than 90 metres (level 1 = surface level)
    n90depth = 1
    DO k = 1, ke
       if (tiestw(k) .lt. 90) n90depth = k
    ENDDO

    ! mz_bk_20120804+
    ! set up all the parameters geiven in [d]^-1 to [step]^-1
    remido_dtb   = remido   * dtb
    dyphy_dtb    = dyphy    * dtb
    grazra_dtb   = grazra   * dtb
    spemor_dtb   = spemor   * dtb
    gammap_dtb   = gammap   * dtb
    gammaz_dtb   = gammaz   * dtb
    drempoc_dtb  = drempoc  * dtb
    dremdoc_dtb  = dremdoc  * dtb
    dremopal_dtb = dremopal * dtb
    dremcalc_dtb = dremcalc * dtb
    dremn2o_dtb  = dremn2o  * dtb
    dphymor_dtb  = dphymor  * dtb
    dzoomor_dtb  = dzoomor  * dtb
    relaxfe_dtb  = relaxfe  * dtb
    ! mz_bk_20120804-
    !__________________________________________________________________________!
    !
    ! Allocate memory : biology
    !__________________________________________________________________________!

    ALLOCATE (kbo(ie,je))
    ALLOCATE (strahl(ie,je))
    ALLOCATE (bolay(ie,je))
    ALLOCATE (alar1max(ie,je))
    ALLOCATE (TSFmax(ie,je))
    ALLOCATE (TMFmax(ie,je))
    !__________________________________________________________________________!
    !
    ! Allocate memory : sediment
    !__________________________________________________________________________!

    ALLOCATE (sedlay(nsedtra))
    ALLOCATE (burial(nsedtra))
    ALLOCATE (powtra(npowtra))

    ALLOCATE (silpro(ie,je))
    ALLOCATE (prorca(ie,je))
    ALLOCATE (pror13(ie,je))
    ALLOCATE (pror14(ie,je))
    ALLOCATE (prcaca(ie,je))
    ALLOCATE (prca13(ie,je))
    ALLOCATE (prca14(ie,je))
    ALLOCATE (produs(ie,je))
    ALLOCATE (seddzi(ksp))
    ALLOCATE (seddw(ks))
    ALLOCATE (porsol(ks))
    ALLOCATE (porwah(ks))
    ALLOCATE (porwat(ks))
    !__________________________________________________________________________!
    !
    ! Allocate memory : inorganic carbon cycle
    !__________________________________________________________________________!

    ALLOCATE (ocetra(nocetra))

    ALLOCATE (c14pool(ie,je,ke))
    ALLOCATE (chemcm(ie,je,8))
!!$!!        ALLOCATE (sedfluxo(ie,je,npowtra))
    ALLOCATE (satn2o(ie,je))
    ALLOCATE (aksp(ie,je,ke))
    ALLOCATE (satoxy(ie,je,ke))
    ALLOCATE (ak23(ie,je,ke))
    ALLOCATE (ak13(ie,je,ke))
    ALLOCATE (akb3(ie,je,ke))
    ALLOCATE (akw3(ie,je,ke))

    IF (L_DUST) THEN
       ALLOCATE (dusty(ie,je))
       ALLOCATE (dustin(ie,je,12))
    ENDIF

    ! mz_bk_20120606+
    ! allocate field for diagnostics
    IF (L_BGC_DIAG) THEN
       ALLOCATE (bgcprod(ibgcprod))
    ENDIF
    ! mz_bk_20120606-

    !__________________________________________________________________________!
    !
    ! Initialize sediment layering
    !__________________________________________________________________________!

    CALL HAMOCC_BODENSED(ie,je,ke,ddpo)
    !__________________________________________________________________________!
    !
    ! Initialize dust from external file.
    !__________________________________________________________________________!

    IF (L_DUST) CALL HAMOCC_GET_DUST(ie,je,ke,ie_g,je_g,ddpo)
    !__________________________________________________________________________!
    !
    ! Initialize constants.
    !__________________________________________________________________________!

    IF (L_AGG) THEN
       zmini = 8000._dp
       DO j=1,je
          DO i=1,ie
             DO k=1,kbo(i,j)-1
                if (ddpo(i,j,k) .gt. 0.5_dp) then
                   zmini = min(ddpo(i,j,k),zmini)
                endif
             ENDDO
          ENDDO
       ENDDO
       CALL global_min(zmini)
    ELSE
       zmini = 0._dp
    ENDIF

    CALL HAMOCC_BELEG(ie,je,ke,ddpo,zmini)
    !__________________________________________________________________________!
    !
    ! Initialize variables.
    !__________________________________________________________________________!

    !-- abs_oce - light absorption --------------------------------------------
    CALL new_channel_object(status, modstr                                    &
         &                , 'abs_oce'                                         &
         &                , p3=abs_oce, reprid = GP_3D_MPIOM)
    CALL new_attribute(status, modstr , 'abs_oce'                             &
         &           , 'long_name', c='ocean light absorption' )
    CALL channel_halt(substr, status)
    !--------------------------------------------------------------------------

    !-- co2 - surface concentration -------------------------------------------
    CALL new_channel_object(status, modstr                                    &
         &                , 'co2'                                             &
         &                , p2=co2, reprid = GP_2D_MPIOM)
    CALL new_attribute(status, modstr , 'co2'                                 &
         &           , 'long_name', c='surface CO2 concentration' )
    CALL channel_halt(substr, status)
    !--------------------------------------------------------------------------

    !-- sediment pH -----------------------------------------------------------
    CALL new_channel_object(status, modstr                                    &
         &                , 'sedhpl'                                          &
         &                , p3=sedhpl, reprid = GP_3D_SED)
    CALL new_attribute(status, modstr , 'sedhpl'                              &
         &           , 'long_name', c='sediment PH' )
    CALL channel_halt(substr, status)
    !--------------------------------------------------------------------------

    !-- detritus in sediment --------------------------------------------------
    CALL new_channel_object(status, modstr                                    &
         &                , 'DETs'                                            &
         &                , p3=sedlay(issso12)%ptr, reprid = GP_3D_SED)
    CALL new_attribute(status, modstr , 'DETs'                                &
         &           , 'long_name', c='detritus sediment' )
    CALL channel_halt(substr, status)
    !--------------------------------------------------------------------------

    !-- CaCO3 in sediment -----------------------------------------------------
    CALL new_channel_object(status, modstr                                    &
         &                , 'CACO3s'                                          &
         &                , p3=sedlay(isssc12)%ptr , reprid = GP_3D_SED)
    CALL new_attribute(status, modstr , 'CACO3s'                              &
         &           , 'long_name', c='calcium carbonate sediment' )
    CALL channel_halt(substr, status)
    !--------------------------------------------------------------------------

    !-- opal in sediment ------------------------------------------------------
    CALL new_channel_object(status, modstr                                    &
         &                , 'OPALs'                                           &
         &                , p3=sedlay(issssil)%ptr , reprid = GP_3D_SED)
    CALL new_attribute(status, modstr , 'OPALs'                               &
         &           , 'long_name', c='opal sediment' )
    CALL channel_halt(substr, status)
    !--------------------------------------------------------------------------

    !-- dust in sediment ------------------------------------------------------
    CALL new_channel_object(status, modstr                                    &
         &                , 'DUSTs' &
         &                , p3=sedlay(issster)%ptr , reprid = GP_3D_SED)
    CALL new_attribute(status, modstr , 'DUSTs'                               &
         &           , 'long_name'                                            &
         &           , c='clay(from free and aggregate dust) sediment' )
    CALL channel_halt(substr, status)
    !--------------------------------------------------------------------------

    !-- C12 in sediment -------------------------------------------------------
    CALL new_channel_object(status, modstr                                    &
         &                , 'C12s'                                            &
         &                , p3=powtra(ipowaic)%ptr, reprid = GP_3D_SED)
    CALL new_attribute(status, modstr , 'C12s'                                &
         &           , 'long_name'                                            &
         &           , c='dissolved inorganic carbon (pore water)' )
    CALL channel_halt(substr, status)
    !--------------------------------------------------------------------------

    !-- alkalinity of pore water ----------------------------------------------
    CALL new_channel_object(status, modstr                                    &
         &                , 'alks'                                            &
         &                , p3=powtra(ipowaal)%ptr, reprid = GP_3D_SED)
    CALL new_attribute(status, modstr , 'alks'                                &
         &           , 'long_name', c='alkalinity (pore water)' )
    CALL channel_halt(substr, status)
    !--------------------------------------------------------------------------

    !-- phosphate in pore water -----------------------------------------------
    CALL new_channel_object(status, modstr                                    &
         &                , 'PO4s'                                            &
         &                , p3=powtra(ipowaph)%ptr, reprid = GP_3D_SED)
    CALL new_attribute(status, modstr , 'PO4s'                                &
         &           , 'long_name', c='phospate (pore water)' )
    CALL channel_halt(substr, status)
    !--------------------------------------------------------------------------

    !-- oxygen in pore water --------------------------------------------------
    CALL new_channel_object(status, modstr                                    &
         &                , 'O2s'                                             &
         &                , p3=powtra(ipowaox)%ptr, reprid = GP_3D_SED)
    CALL new_attribute(status, modstr , 'O2s'                                 &
         &           , 'long_name', c='oxygen (pore water)' )
    CALL channel_halt(substr, status)
    !--------------------------------------------------------------------------

    !-- nitrogen in pore water ------------------------------------------------
    CALL new_channel_object(status, modstr                                    &
         &                , 'N2s'                                             &
         &                , p3=powtra(ipown2)%ptr, reprid = GP_3D_SED)
    CALL new_attribute(status, modstr , 'N2s'                                 &
         &           , 'long_name', c='N2 (pore water)' )
    CALL channel_halt(substr, status)
    !--------------------------------------------------------------------------

    !-- nitrate in pore water -------------------------------------------------
    CALL new_channel_object(status, modstr                                    &
         &                , 'NO3s'                                            &
         &                , p3=powtra(ipowno3)%ptr, reprid = GP_3D_SED)
    CALL new_attribute(status, modstr , 'NO3s'                                &
         &           , 'long_name', c='nitrate (pore water)' )
    CALL channel_halt(substr, status)
    !--------------------------------------------------------------------------

    !-- silicate in pore water ------------------------------------------------
    CALL new_channel_object(status, modstr                                    &
         &                , 'SiOHs' &
         &                , p3=powtra(ipowasi)%ptr, reprid = GP_3D_SED)
    CALL new_attribute(status, modstr , 'SiOHs'                               &
         &           , 'long_name', c='silicate (pore water)' )
    CALL channel_halt(substr, status)
    !--------------------------------------------------------------------------

    !-- burried detritus in sediment ------------------------------------------
    CALL new_channel_object(status, modstr                                    &
         &                , 'DETb'                                            &
         &                , p2=burial(issso12)%ptr, reprid = GP_2D_MPIOM)
    CALL new_attribute(status, modstr , 'DETb'                                &
         &           , 'long_name', c='detritus sediment (burial)' )
    CALL channel_halt(substr, status)
    !--------------------------------------------------------------------------

    !-- burried CaCO3 in sediment ---------------------------------------------
    CALL new_channel_object(status, modstr                                    &
         &                , 'CACO3b'                                          &
         &                , p2=burial(isssc12)%ptr , reprid = GP_2D_MPIOM)
    CALL new_attribute(status, modstr , 'CACO3b'                              &
         &           , 'long_name', c='calcium carbonate sediment (burial)' )
    CALL channel_halt(substr, status)
    !--------------------------------------------------------------------------

    !-- burried opal in sediment ----------------------------------------------
    CALL new_channel_object(status, modstr                                    &
         &                , 'OPALb'                                           &
         &                , p2=burial(issssil)%ptr , reprid = GP_2D_MPIOM)
    CALL new_attribute(status, modstr , 'OPALb'                               &
         &           , 'long_name', c='opal sediment (burial)' )
    CALL channel_halt(substr, status)
    !--------------------------------------------------------------------------

    !-- burried dust in sediment ----------------------------------------------
    CALL new_channel_object(status, modstr                                    &
         &                , 'DUSTb'                                           &
         &                , p2=burial(issster)%ptr , reprid = GP_2D_MPIOM)
    CALL new_attribute(status, modstr , 'DUSTb'                               &
         &         , 'long_name'                                              &
         &         , c='clay(from free and aggregate dust) sediment (burial)' )
    CALL channel_halt(substr, status)
    !--------------------------------------------------------------------------

    ! set pointers from internal fields to ocean tracer fields
    ocetra(isco212)%ptr => xt_om(:,:,idx_sco212,:)
    ocetra(ialkali)%ptr => xt_om(:,:,idx_alkali,:)
    ocetra(iphosph)%ptr => xt_om(:,:,idx_phosph,:)
    ocetra(ioxygen)%ptr => xt_om(:,:,idx_oxygen,:)
    ocetra(igasnit)%ptr => xt_om(:,:,idx_gasnit,:)
    ocetra(iano3  )%ptr => xt_om(:,:,idx_ano3  ,:)
    ocetra(isilica)%ptr => xt_om(:,:,idx_silica,:)
    ocetra(idoc   )%ptr => xt_om(:,:,idx_doc   ,:)
    ocetra(iphy   )%ptr => xt_om(:,:,idx_phy   ,:)
    ocetra(izoo   )%ptr => xt_om(:,:,idx_zoo   ,:)
    ocetra(ian2o  )%ptr => xt_om(:,:,idx_an2o  ,:)
    ocetra(idms   )%ptr => xt_om(:,:,idx_dms   ,:)
    ocetra(iiron  )%ptr => xt_om(:,:,idx_iron  ,:)
    ocetra(ifdust )%ptr => xt_om(:,:,idx_fdust ,:)
    ocetra(idet   )%ptr => xt_om(:,:,idx_det   ,:)
    ocetra(icalc  )%ptr => xt_om(:,:,idx_calc  ,:)
    ocetra(iopal  )%ptr => xt_om(:,:,idx_opal  ,:)
    hi                  => xt_om(:,:,idx_hp    ,:)
    co3                 => xt_om(:,:,idx_co3mm ,:)

    IF (L_AGG) THEN
       ocetra(inos)%ptr => xt_om(:,:,idx_nos  ,:)
       ocetra(iadust)%ptr => xt_om(:,:,idx_adust  ,:)
    ENDIF

    ! mz_bk_20120606+
    !-- diagnostic fields -----------------------------------------------------
    IF (L_BGC_DIAG) THEN
       !-- DMS production -----------------------------------------------------
       CALL new_channel_object(status, modstr                                 &
            &                , 'dmsprod'                                      &
            &                , p3=bgcprod(kdmsprod)%ptr, reprid = GP_3D_MPIOM &
            &                , lrestreq=.FALSE.)
       CALL new_attribute(status, modstr , 'dmsprod'                          &
            &           , 'long_name', c='DMS production' )
       CALL new_attribute(status, modstr , 'dmsprod'                          &
            &           , 'units', c='kmol(S)/m^3/s' )
       CALL channel_halt(substr, status)
       !-----------------------------------------------------------------------

       !-- DMS bacterial consumption ------------------------------------------
       CALL new_channel_object(status, modstr                                 &
            &                , 'dms_bac'                                      &
            &                , p3=bgcprod(kdms_bac)%ptr, reprid = GP_3D_MPIOM &
            &                , lrestreq=.FALSE.)
       CALL new_attribute(status, modstr , 'dms_bac'                          &
            &           , 'long_name', c='DMS bacterial consumption' )
       CALL new_attribute(status, modstr , 'dms_bac'                          &
            &           , 'units', c='kmol(S)/m^3/s' )
       CALL channel_halt(substr, status)
       !-----------------------------------------------------------------------

       !-- DMS destruction by UV radiation ------------------------------------
       CALL new_channel_object(status, modstr                                 &
            &                , 'dms_uv'                                       &
            &                , p3=bgcprod(kdms_uv)%ptr, reprid = GP_3D_MPIOM  &
            &                , lrestreq=.FALSE.)
       CALL new_attribute(status, modstr , 'dms_uv'                           &
            &           , 'long_name', c='DMS decay due to UV radiation' )
       CALL new_attribute(status, modstr , 'dms_uv'                           &
            &           , 'units', c='kmol(S)/m^3/s' )
       CALL channel_halt(substr, status)
       !-----------------------------------------------------------------------

       !-- export production --------------------------------------------------
       CALL new_channel_object(status, modstr                                 &
            &                , 'export'                                       &
            &                , p3=bgcprod(kexport)%ptr, reprid = GP_3D_MPIOM  &
            &                , lrestreq=.FALSE.)
       CALL new_attribute(status, modstr , 'export'                           &
            &           , 'long_name', c='carbon export' )
       CALL new_attribute(status, modstr , 'export'                           &
            &           , 'units', c='kmol(C)/m^3/s' )
       CALL channel_halt(substr, status)
       !-----------------------------------------------------------------------

       !-- production of calcium carbonate ------------------------------------
       CALL new_channel_object(status, modstr                                 &
            &                , 'delcar'                                       &
            &                , p3=bgcprod(kdelcar)%ptr, reprid = GP_3D_MPIOM  &
            &                , lrestreq=.FALSE.)
       CALL new_attribute(status, modstr , 'delcar'                           &
            &           , 'long_name', c='calcium carbonate production' )
       CALL new_attribute(status, modstr , 'delcar'                           &
            &           , 'units', c='kmol(C)/m^3/s' )
       CALL channel_halt(substr, status)
       !-----------------------------------------------------------------------

       !-- production of opal -------------------------------------------------
       CALL new_channel_object(status, modstr                                 &
            &                , 'delsil'                                       &
            &                , p3=bgcprod(kdelsil)%ptr, reprid = GP_3D_MPIOM  &
            &                , lrestreq=.FALSE.)
       CALL new_attribute(status, modstr , 'delsil'                           &
            &           , 'long_name', c='opal production' )
       CALL new_attribute(status, modstr , 'delsil'                           &
            &           , 'units', c='kmol(Si)/m^3/s' )
       CALL channel_halt(substr, status)
       !-----------------------------------------------------------------------

       !-- photosenthesis -----------------------------------------------------
       CALL new_channel_object(status, modstr                                 &
            &                , 'phosy'                                        &
            &                , p3=bgcprod(kphosy)%ptr, reprid = GP_3D_MPIOM   &
            &                , lrestreq=.FALSE.)
       CALL new_attribute(status, modstr , 'phosy'                            &
            &           , 'long_name', c='photosynthesis' )
       CALL new_attribute(status, modstr , 'phosy'                            &
            &           , 'units', c='kmol(P)/m^3/s' )
       CALL channel_halt(substr, status)
       !-----------------------------------------------------------------------

       !-- grazing ------------------------------------------------------------
       CALL new_channel_object(status, modstr                                 &
            &                , 'grazing'                                      &
            &                , p3=bgcprod(kgraz)%ptr, reprid = GP_3D_MPIOM    &
            &                , lrestreq=.FALSE.)
       CALL new_attribute(status, modstr , 'grazing'                          &
            &           , 'long_name', c='grazing by zooplankton' )
       CALL new_attribute(status, modstr , 'grazing'                          &
            &           , 'units', c='kmol(P)/m^3/s' )
       CALL channel_halt(substr, status)
       !-----------------------------------------------------------------------

       !-- iron limited zones -------------------------------------------------
       CALL new_channel_object(status, modstr                                 &
            &                , 'flim'                                         &
            &                , p3=bgcprod(kflim)%ptr, reprid = GP_3D_MPIOM    &
            &                , lrestreq=.FALSE.)
       CALL new_attribute(status, modstr , 'flim'                             &
            &           , 'long_name', c='iron limitation' )
       CALL new_attribute(status, modstr , 'flim'                             &
            &           , 'units', c='' )
       CALL channel_halt(substr, status)
       !-----------------------------------------------------------------------

       !-- phosphorous limited zones ------------------------------------------
       CALL new_channel_object(status, modstr                                 &
            &                , 'plim'                                         &
            &                , p3=bgcprod(kplim)%ptr, reprid = GP_3D_MPIOM    &
            &                , lrestreq=.FALSE.)
       CALL new_attribute(status, modstr , 'plim'                             &
            &           , 'long_name', c='phosphorous limitation' )
       CALL new_attribute(status, modstr , 'plim'                             &
            &           , 'units', c='' )
       CALL channel_halt(substr, status)
       !-----------------------------------------------------------------------

       !-- nitrogen limited zones ---------------------------------------------
       CALL new_channel_object(status, modstr                                 &
            &                , 'nlim'                                         &
            &                , p3=bgcprod(knlim)%ptr, reprid = GP_3D_MPIOM    &
            &                , lrestreq=.FALSE.)
       CALL new_attribute(status, modstr , 'nlim'                             &
            &           , 'long_name', c='nitrogen limitation' )
       CALL new_attribute(status, modstr , 'nlim'                             &
            &           , 'units', c='' )
       CALL channel_halt(substr, status)
       !-----------------------------------------------------------------------
    ENDIF  ! L_BGC_DIAG
    ! mz_bk_20120606-

    CALL end_message_bi(modstr, 'MEMORY INITIALIZATION', substr)

  END SUBROUTINE hamocc_init_memory
  ! ---------------------------------------------------------------------------

  ! ---------------------------------------------------------------------------
  SUBROUTINE hamocc_init_coupling

    ! MESSy
    USE messy_main_mpi_bi,        ONLY: p_parallel_io
    USE messy_main_tracer_mem_bi, ONLY: xt_om
    USE messy_main_channel,       ONLY: get_channel_object, get_channel_info  &
         &                            , get_channel_object_info

    IMPLICIT NONE

    CHARACTER(len=*), PARAMETER  :: substr = 'hamocc_init_coupling'

    INTEGER                 :: status
    INTEGER :: reprid

    ! CALL get_channel_info(status, 'otphysc')
    ! IF (status /= 0) THEN
    !    CALL finish(substr,'No Ocean Trace Physics (OTPHYSC) running... STOP')
    CALL message(substr, ' IS OTPHYSC RUNNING??????????')  
    ! END IF

    ! Test if channel objects for coupling present
    CALL get_channel_object(status                                            &
         &                , oname='tho', cname='mpiom'                        & 
         &                , p3=tho )
    IF (status /= 0) THEN
       CALL finish(substr,'ERROR in the coupling procedure with MPIOM(THO)')
    END IF

    CALL get_channel_object(status                                            &
         &                , oname='sao', cname='mpiom'                        &
         &                , p3=sao )
    IF (status /= 0) THEN
       CALL finish(substr,'ERROR in the coupling procedure with MPIOM(SAO)')
    END IF

    CALL get_channel_object(status                                            &
         &                , oname='AOFLSHWO', cname='a2o'                     &
         &                , p2=fswr_in )
    IF (status /= 0) THEN
       ! N.B. Sadly fswr is not a channel! call error if it does not exist 
       ! (i.e.) no coupling ECHAM5-->MPIOM
       CALL finish(substr,'ERROR in the coupling procedure with a2o(AOFLSHWO)')
    END IF

    IF (.not.L_DUST) THEN
       CALL message(substr, 'IMPORT DUST ... ')
       CALL get_channel_info(status, dust_input(1))
       IF (status /= 0) THEN
          IF (p_parallel_io) THEN
             WRITE(*,*) ' No DUST available from other submodel :  '
             WRITE(*,*) ' channel not present ! (check namelist)'
             CALL finish(substr, 'dust_input channel (1) not present!')
          END IF
       ELSE
          CALL get_channel_object_info(status                                 &
               &                     , cname=dust_input(1)                    &
               &                     , oname=dust_input(2)                    & 
               &                     , reprid=reprid )
          IF (status /= 0) THEN
             IF (p_parallel_io) THEN
                WRITE(*,*) ' No DUST available from other submodels '
                WRITE(*,*) ' object not present ! (check namelist)'
                CALL finish(substr, 'dust_input channel (2) not present!')
             END IF
          ELSE
             IF (p_parallel_io) THEN
                WRITE(*,*) ' Oceanic DUST input available from channel object'&
                     &     , dust_input(1), '%', dust_input(2) 
             END IF
             CALL get_channel_object(status, cname=dust_input(1)              &
                  &                , oname=dust_input(2)                      &
                  &                , p2=dusty )
          END IF
       END IF
    END IF

  END SUBROUTINE hamocc_init_coupling
  ! ---------------------------------------------------------------------------

  ! ---------------------------------------------------------------------------
  SUBROUTINE hamocc_global_start

    USE messy_main_timer_bi,      ONLY: event_state
    USE messy_main_timer,         ONLY: current_date, previous_date           &
         &                            , time_days, MONTH, lresume,lstart 
    USE messy_main_mpi_bi,    ONLY: p_parallel_io, p_io, p_bcast, finish, p_pe

    ! TODO : to be eliminated --> through interface!
    USE messy_mpiom_mem_e5,   ONLY: ddpo, dpio, tiestu, sicomo, ie, je, ke    &
         &                        , ie_g,je_g, gather_arr
    USE messy_main_tracer_mem_bi, ONLY: xt_om

    IMPLICIT NONE

    INTEGER :: i,j,k,l
    REAL(dp), DIMENSION(:,:), POINTER  :: check_field => NULL()

    !TODO: for rerun:
    ! * all OCETRA (tracer) [x]
    ! * hi  [H+]            [x]
    ! * co3 [CO3--]         [x]
    ! * all SEDLAY          [x] 
    ! * all BURIAL          [x]
    ! * sedhpl              [x]
    ! * all POWTRA          [x]

    l_trig_hamocc       = event_state(ev_trig_hamocc, previous_date) 
    ! l_trig_hamocc_day   = event_state(ev_trig_hamocc_day, current_date)
    l_trig_hamocc_day   = event_state(ev_trig_hamocc_day, previous_date)


    ! check for tracer initialisation in 3 steps:
    ! 1) if rerun -> no initialisation
    ! 2) if start -> overwrite value : new start! we need initialisation
    ! 3) if tracer = 0.0_dp -> new start (see later) or rerun without hamocc
    !                          restart files! we need initialisation!

    l_init_bgc = .FALSE.

    IF (lresume) l_init_bgc = .FALSE.
    IF (lstart)  l_init_bgc = .TRUE.

    ! mz_bk_20120721+
    IF (lresume) THEN
    ! mz_bk_20120721-
       ALLOCATE(check_field(ie_g,je_g))
       CALL gather_arr(ocetra(isco212)%ptr(:,:,1),check_field,p_io)
       IF (p_pe==p_io) THEN
          IF (MAXVAL(check_field) .eq. 0.0_dp) THEN
             ! mz_bk_20120609+
             l_init_bgc = .TRUE.
             write(*,*) " TRACER MISSING FROM HAMOCC RESTART!"
             ! mz_bk_20120609-
          ENDIF
       ENDIF
       CALL p_bcast(l_init_bgc,   p_io)
       DEALLOCATE(check_field)
    ! mz_bk_20120721+
    ENDIF
    ! mz_bk_20120721-

    IF (l_init_bgc) THEN
       ! here we have a new start or no restart file present
       ! initialise ocean absorption with 1
       abs_oce(:,:,:) = 1._dp
       CALL HAMOCC_INIT
       ! initialize with constants if no files present!
       IF (l_init_bgc) THEN
          IF (p_pe==p_io) THEN
             write(*,*) " WARNING: Missing initialisation file for HAMOCC..."
             write(*,*) "   ... initialisation with constants!     "
          END IF
          CALL HAMOCC_INIT_VAR(ie,je,ke,ddpo)
          l_init_bgc = .FALSE.
       ENDIF
    ENDIF

    IF (L_SKIP_INT) return   ! mz_bk_20100519

    IF (l_trig_hamocc .and. .not. lstart) THEN  ! HAMOCC TIME LOOP 

       IF (p_pe==p_io) THEN
          WRITE(*,*) " HAMOCC TIME STEP"
       ENDIF

       !!--------------------------------------------------------------------
       !!                     BGC routines
       !!--------------------------------------------------------------------
       !--------------------------------------------------------------------
       ! Net solar radiation: multiply  with sea ice concentration

       DO  j=1,je
          DO  i=1,ie
             ! IT IS NEEDED IN OCPROD
             ! NB: it is not updated before a2o...
             strahl(i,j) = fswr_in(i,j) * (1._dp - sicomo(i,j))
          ENDDO
       ENDDO
       !---------------------------------------------------------------
       ! CALCULATION OF REACTION CONSTANTS
       !
       CALL HAMOCC_CHEMCON(ie,je,ke,sao,tho,ddpo,tiestu)  

       IF (L_CHECK) CALL hamocc_check("CHEMCON")
 
       !---------------------------------------------------------------
       ! BIOGEOCHEMISTRY: OCEAN PRODUCTION 
       !
       IF (L_DUST) dusty(:,:) = dustin(:,:,MONTH)
       CALL HAMOCC_OCPROD(IE,JE,KE,THO,ddpo,dpio)

       IF (L_CHECK) CALL hamocc_check("OCPROD")

       !---------------------------------------------------------------
       ! NITRATE REDUCTION by cyano bacteria 
       ! (2NO3 + O2 => N2O + O2).
       !
       CALL HAMOCC_CYANO(IE,JE,KE,ddpo)

       IF (L_CHECK) CALL hamocc_check("CYANO")

       !---------------------------------------------------------------
       ! INORGANIC CARBON CYCLE.
       !
       CALL HAMOCC_CARCHM(ie,je,ke,sao,tho,ddpo)

       IF (L_CHECK) CALL hamocc_check("CARCHM")

       !---------------------------------------------------------------
       ! SEDIMENT MODULE
       !
       CALL HAMOCC_POWACH(ie,je,ke,sao)

       IF (L_CHECK) CALL hamocc_check("POWAC")

       !---------------------------------------------------------------
       ! SEDIMENT SHIFTING
       !
       IF (l_trig_hamocc_day) THEN
          CALL HAMOCC_SEDSHI(ie,je)
       ENDIF
 
       IF (L_CHECK) CALL hamocc_check("SEDSHI")

       !---------------------------------------------------------------
       ! TRACER PHOTOPRODUCTION
       ! OR DIAGNOSTIC TRACERS

       CALL HAMOCC_PHOTOPROD(ie, je, ke                             & ! indices
            &              , xt_om(:,:,idx_co,:)                    & ! CO
            &              , xt_om(:,:,idx_isop,:)                  & ! isop
            &              , xt_om(:,:,idx_ch3i,:)                  & ! CH3I
            &              , sao(:,:,:), tho(:,:,:)+273.15_dp)

       IF (L_CHECK) CALL hamocc_check("PHOTO")
       !---------------------------------------------------------------

    ENDIF  ! HAMOCC TIME LOOP 

  END SUBROUTINE hamocc_global_start
  ! ---------------------------------------------------------------------------

  ! ---------------------------------------------------------------------------
  SUBROUTINE hamocc_global_end

    USE messy_main_timer_bi,      ONLY: event_state
    USE messy_main_timer,         ONLY: current_date, previous_date

    ! l_trig_hamocc_day   = event_state(ev_trig_hamocc_day, current_date)     &
    !      &                .OR. lstart

    ! INTEGER :: I,J,K
     
  END SUBROUTINE hamocc_global_end
  ! ---------------------------------------------------------------------------

  ! ---------------------------------------------------------------------------
  SUBROUTINE hamocc_free_memory

    USE messy_main_tracer_mem_bi,   ONLY: ntrac_om

    IMPLICIT NONE
    INTRINSIC ASSOCIATED

    INTEGER :: I

    IF (ASSOCIATED(dzs))        DEALLOCATE(dzs)
    NULLIFY(dzs)

    IF (ASSOCIATED(kbo))        DEALLOCATE(kbo)
    NULLIFY(kbo)
    IF (ASSOCIATED(strahl))     DEALLOCATE(strahl)
    NULLIFY(strahl)
    IF (ASSOCIATED(bolay))      DEALLOCATE(bolay)
    NULLIFY(bolay)
    IF (ASSOCIATED(alar1max))   DEALLOCATE(alar1max)
    NULLIFY(alar1max)
    IF (ASSOCIATED(TSFmax))     DEALLOCATE(TSFmax)
    NULLIFY(TSFmax)
    IF (ASSOCIATED(TMFmax))     DEALLOCATE(TMFmax)
    NULLIFY(TMFmax)

    IF (ASSOCIATED(silpro))     DEALLOCATE(silpro)
    NULLIFY(silpro)
    IF (ASSOCIATED(prorca))     DEALLOCATE(prorca)
    NULLIFY(prorca)
    IF (ASSOCIATED(pror13))     DEALLOCATE(pror13)
    NULLIFY(pror13)
    IF (ASSOCIATED(pror14))     DEALLOCATE(pror14)
    NULLIFY(pror14)
    IF (ASSOCIATED(prcaca))     DEALLOCATE(prcaca)
    NULLIFY(prcaca)
    IF (ASSOCIATED(prca13))     DEALLOCATE(prca13)
    NULLIFY(prca13)
    IF (ASSOCIATED(prca14))     DEALLOCATE(prca14)
    NULLIFY(prca14)
    IF (ASSOCIATED(produs))     DEALLOCATE(produs)
    NULLIFY(produs)
    IF (ASSOCIATED(seddzi))     DEALLOCATE(seddzi)
    NULLIFY(seddzi)
    IF (ASSOCIATED(seddw))      DEALLOCATE(seddw)
    NULLIFY(seddw)
    IF (ASSOCIATED(porsol))     DEALLOCATE(porsol)
    NULLIFY(porsol)
    IF (ASSOCIATED(porwah))     DEALLOCATE(porwah)
    NULLIFY(porwah)
    IF (ASSOCIATED(porwat))     DEALLOCATE(porwat)
    NULLIFY(porwat)

    !special
    IF (ASSOCIATED(sedlay))     DEALLOCATE(sedlay)
    NULLIFY(sedlay)
    IF (ASSOCIATED(burial))     DEALLOCATE(burial)
    NULLIFY(burial)
    IF (ASSOCIATED(powtra))     DEALLOCATE(powtra)
    NULLIFY(powtra)
    IF (ASSOCIATED(ocetra))     DEALLOCATE(ocetra)
    NULLIFY(ocetra)

    IF (ASSOCIATED(c14pool))    DEALLOCATE(c14pool)
    NULLIFY(c14pool)
    IF (ASSOCIATED(chemcm))     DEALLOCATE(chemcm)
    NULLIFY(chemcm)
!!$!!        ALLOCATE(sedfluxo)
    IF (ASSOCIATED(satn2o))     DEALLOCATE(satn2o)
    NULLIFY(satn2o)
    IF (ASSOCIATED(aksp))       DEALLOCATE(aksp)
    NULLIFY(aksp)
    IF (ASSOCIATED(satoxy))     DEALLOCATE(satoxy)
    NULLIFY(satoxy)
    IF (ASSOCIATED(ak23))       DEALLOCATE(ak23)
    NULLIFY(ak23)
    IF (ASSOCIATED(ak13))       DEALLOCATE(ak13)
    NULLIFY(ak13)
    IF (ASSOCIATED(akb3))       DEALLOCATE(akb3)
    NULLIFY(akb3)
    IF (ASSOCIATED(akw3))       DEALLOCATE(akw3)
    NULLIFY(akw3)

    IF (L_DUST) THEN 
       IF (ASSOCIATED(dusty))   DEALLOCATE(dusty)
       NULLIFY(dusty)
       IF (ASSOCIATED(dustin))  DEALLOCATE(dustin)
       NULLIFY(dustin)
    ENDIF

    ! mz_bk_20120606+
    IF (L_BGC_DIAG) THEN
       IF (ASSOCIATED(bgcprod)) DEALLOCATE(bgcprod)
       NULLIFY(bgcprod)
    ENDIF
    ! mz_bk_20120606-

  END SUBROUTINE hamocc_free_memory
  ! ---------------------------------------------------------------------------


  ! ===========================================================================
  ! PRIVATE HAMOCC INTERFACE ROUTINES
  ! ===========================================================================
  !
  ! ---------------------------------------------------------------------------

  SUBROUTINE hamocc_check(text)

    USE messy_main_mpi_bi,    ONLY: p_parallel_io, p_io, p_bcast, finish, p_pe

    IMPLICIT NONE

    CHARACTER(LEN=*), INTENT(IN) :: text
    IF (p_pe==p_io) THEN
       WRITE(*,*) "HAMOCC at: "
       WRITE(*,*) text
       WRITE(*,*) "------------------------------------------"
       WRITE(*,*) "------------------------------------------"
    ENDIF

  END SUBROUTINE hamocc_check
  ! ---------------------------------------------------------------------------
  
  ! ---------------------------------------------------------------------------
  SUBROUTINE hamocc_init

    USE messy_main_mpi_bi,        ONLY: p_parallel_io, p_io, p_bcast          &
         &                            , finish, p_pe
    USE messy_mpiom_mem_e5,       ONLY: ke
    USE messy_main_tracer_mem_bi, ONLY: xt_om

    implicit none

    INCLUDE 'netcdf.inc'
    INTEGER :: ncid,ncstat,ncvarid
    CHARACTER(LEN=400) :: input_file
    CHARACTER(LEN=8) :: date
    !
    ! Open netCDF data file
    !      
    IF(p_pe==p_io) THEN
       write (date(1:4),'(I4.4)') YEAR
       write (date(5:6),'(I2.2)') MONTH
       write (date(7:8),'(I2.2)') DAY
       input_file = TRIM(hamocc_ini_files_path)//'/hamocc_ini_'               &
            &       //TRIM(date)//'.nc'  
       ncstat = NF_OPEN(TRIM(input_file),NF_NOWRITE, ncid)
       IF (ncstat.NE.NF_NOERR ) THEN
          ! If we do not have initialization files initialize with constants!
          CALL message('  '                                                   &
               &     , 'WARNING: INIT OCEAN BGC variables: Problem with netCDF')
          l_init_bgc=.TRUE.
       ELSE
          write(*,*) " Obtain hamocc fields from input files..."
          l_init_bgc=.FALSE.
       ENDIF
    END IF

    CALL p_bcast(l_init_bgc,   p_io)   
    IF (l_init_bgc) RETURN

    !
    ! Read netCDF data file
    !
    !tracers
    call read_netcdf_var('DMS',  ocetra(idms)%ptr   ,ncid, ke)
    call read_netcdf_var('DUST', ocetra(ifdust)%ptr ,ncid, ke)
    call read_netcdf_var('FE',   ocetra(iiron)%ptr  ,ncid, ke)
    call read_netcdf_var('PO4',  ocetra(iphosph)%ptr,ncid, ke)
    call read_netcdf_var('N2O',  ocetra(ian2o)%ptr  ,ncid, ke)
    call read_netcdf_var('NO3',  ocetra(iano3)%ptr  ,ncid, ke)
    call read_netcdf_var('N2',   ocetra(igasnit)%ptr,ncid, ke)
    call read_netcdf_var('ZOO',  ocetra(izoo)%ptr   ,ncid, ke)
    call read_netcdf_var('PHY',  ocetra(iphy)%ptr   ,ncid, ke)
    call read_netcdf_var('DOM',  ocetra(idoc)%ptr   ,ncid, ke)
    call read_netcdf_var('DET',  ocetra(idet)%ptr   ,ncid, ke)
    call read_netcdf_var('CACO3',ocetra(icalc)%ptr  ,ncid, ke)
    call read_netcdf_var('C12',  ocetra(isco212)%ptr,ncid, ke)
    call read_netcdf_var('O2',   ocetra(ioxygen)%ptr,ncid, ke)
    call read_netcdf_var('SI',   ocetra(isilica)%ptr,ncid, ke)
    call read_netcdf_var('OPAL', ocetra(iopal)%ptr  ,ncid, ke)
    call read_netcdf_var('AT',   ocetra(ialkali)%ptr,ncid, ke)
    IF (L_AGG) THEN
       call read_netcdf_var('NOS',   ocetra(inos)%ptr,   ncid, ke)
       call read_netcdf_var('ADUST', ocetra(iadust)%ptr, ncid, ke)
    ENDIF
    call read_netcdf_var('CO',   xt_om(:,:,idx_co,:), ncid, ke)
    call read_netcdf_var('HP',   hi,  ncid, ke)
    call read_netcdf_var('CO3MM',co3, ncid, ke)

    !all sedlay
    call read_netcdf_var('DETS'  ,sedlay(issso12)%ptr, ncid, ks)
    call read_netcdf_var('CACO3S',sedlay(isssc12)%ptr, ncid, ks)
    call read_netcdf_var('OPALS' ,sedlay(issssil)%ptr, ncid, ks)
    call read_netcdf_var('DUSTS' ,sedlay(issster)%ptr, ncid, ks)

    !all powtra
    call read_netcdf_var('C12S'  ,powtra(ipowaic)%ptr,  ncid, ks)
    call read_netcdf_var('ALKS'  ,powtra(ipowaal)%ptr,  ncid, ks)
    call read_netcdf_var('PO4S'  ,powtra(ipowaph)%ptr,  ncid, ks)
    call read_netcdf_var('O2S'   ,powtra(ipowaox)%ptr,  ncid, ks)
    call read_netcdf_var('N2S'   ,powtra(ipown2)%ptr,  ncid, ks)
    call read_netcdf_var('NO3S'  ,powtra(ipowno3)%ptr, ncid, ks)
    call read_netcdf_var('SIOHS' ,powtra(ipowasi)%ptr, ncid, ks)

    !all burial
    call read_netcdf_var('DETB'   ,burial(issso12)%ptr, ncid, 1)
    call read_netcdf_var('CACO3B' ,burial(isssc12)%ptr, ncid, 1)
    call read_netcdf_var('OPALB'  ,burial(issssil)%ptr, ncid, 1)
    call read_netcdf_var('DUSTB'  ,burial(issster)%ptr, ncid, 1)

    !sedphl
    call read_netcdf_var('SEDHPL' ,sedhpl, ncid, ks)

    ! Close file
    IF(p_pe==p_io) THEN
       ncstat = NF_CLOSE(ncid)
       IF ( ncstat .NE. NF_NOERR )                                            &
            &         CALL finish('Fields initialization: Problem with netCDF')
    END IF

  END SUBROUTINE hamocc_init
  ! ---------------------------------------------------------------------------

  ! ---------------------------------------------------------------------------
  SUBROUTINE read_netcdf_var(var,arr, ncid, levels)

    USE messy_main_mpi_bi,    ONLY: p_parallel_io, p_io, p_bcast              &
         &                        , finish, p_pe
    USE messy_mpiom_mem_e5,   ONLY: scatter_arr, ddpo, ie, je                 & 
         &                        , ie_g, je_g

    implicit none

    INCLUDE 'netcdf.inc'

    CHARACTER(*), INTENT(IN) :: var
    INTEGER,INTENT(IN)       :: ncid    
    INTEGER,INTENT(IN)       :: levels    
    REAL(DP), INTENT(OUT)    :: arr(ie,je,levels)
    !LOCAL
    INTEGER :: k
    REAL(DP) :: arr_g(ie_g,je_g,levels)
    INTEGER :: ncstat,ncvarid
    CHARACTER (LEN=80) err_text

    ! Read NETCDF data

    IF(p_pe==p_io) THEN
       err_text = 'Problem reading '//var
       ncstat = NF_INQ_VARID(ncid,TRIM(var),ncvarid)
       IF ( ncstat .NE. NF_NOERR ) CALL finish(err_text)
       ncstat = NF_GET_VAR_DOUBLE(ncid,ncvarid,arr_g)
       IF ( ncstat .NE. NF_NOERR ) CALL finish(err_text)
    ENDIF

    DO k=1,levels
       CALL scatter_arr(arr_g(:,:,k),arr(:,:,k),p_io)
    ENDDO

  END SUBROUTINE read_netcdf_var
  ! ---------------------------------------------------------------------------

  ! ---------------------------------------------------------------------------
  SUBROUTINE hamocc_get_dust(ie,je,ke,ie_g,je_g,ddpo)

    USE messy_main_mpi_bi,    ONLY: p_parallel_io, p_io, p_bcast              &
         &                        , finish, p_pe
    USE messy_mpiom_mem_e5,   ONLY: scatter_arr

    implicit none
    INTEGER,  INTENT(IN):: ie, je, ke
    INTEGER,  INTENT(IN):: ie_g, je_g
    REAL(DP), INTENT(IN):: ddpo(ie,je,ke)

    INTEGER :: i, j , k, l

    ! define the fields

    ! annual iron input
    REAL(DP) :: iron_sum

    INCLUDE 'netcdf.inc'
    INTEGER :: ncid,ncstat,ncvarid
    REAL(DP) :: arr(ie,je,12)
    REAL(DP) :: arr_g(ie_g,je_g,12)
    CHARACTER (LEN=80) err_text

    !
    ! Open netCDF data file
    !      
    IF(p_pe==p_io) THEN
       write(*,*) TRIM(hamocc_ini_files_path)//'/INPDUST.nc'
       ncstat = NF_OPEN(TRIM(hamocc_ini_files_path)//'/INPDUST.nc',           &
            &           NF_NOWRITE, ncid)
       IF (ncstat.NE.NF_NOERR ) CALL finish('get_dust: Problem with netCDF')
    END IF
    !
    ! call read_netcdf_var(ncid,'DUST',dustin(1,1,1),12)

    ! Read NETCDF data

    IF(p_pe==p_io) THEN
       err_text = 'Problem reading '//'DUST'
       ncstat = NF_INQ_VARID(ncid,'DUST',ncvarid)
       IF ( ncstat .NE. NF_NOERR ) CALL finish(err_text)
       ncstat = NF_GET_VAR_DOUBLE(ncid,ncvarid,arr_g)
       IF ( ncstat .NE. NF_NOERR ) CALL finish(err_text)
    ENDIF

    ! write (*,*) 'x, y, z,  arr_g'                                           &
    !      &      , SIZE(arr_g,DIM=1), SIZE(arr_g,DIM=2), SIZE(arr_g,DIM=3)
    ! write (*,*) 'x, y, z, dustin'                                           &
    !      &      , SIZE(dustin,DIM=1), SIZE(dustin,DIM=2), SIZE(dustin,DIM=3)
    DO k=1,12
       CALL scatter_arr(arr_g(:,:,k),dustin(:,:,k),p_io)
    ENDDO

    !
    ! Close file
    IF(p_pe==p_io) THEN
       ncstat = NF_CLOSE(ncid)
       IF ( ncstat .NE. NF_NOERR )                                            &
            &                   CALL finish('get_dust: Problem with netCDF200')
    END IF


    ! done before hamocc_ocprod   !mz_ap_20090525
    !! set to missing value (0.) over land
    !      do l=1,12
    !        do j=1,je
    !          do i=1,ie
    !            if(ddpo(i,j,1).gt.0.5) then 
    !              dusty(i,j,l) = dustin(i,j,l)
    !            else 
    !              dusty(i,j,l) = 0.
    !            endif    
    !          enddo
    !         enddo
    !      enddo

    !  can be done offline !mz_ap_20090525
    !! sum up annual iron input from dust (reactivated 16.02.06 js)
    !      iron_sum=0.
    ! 
    !      do l=1,12
    !        do j=1,je
    !          do i=1,ie
    !            if(ddpo(i,j,1).gt.0.5) then
    !              iron_sum = iron_sum + dustin(i,j,l) * 0.00035 / 55.85
    !! js 070706 removed /12 since input is per month
    !            endif
    !          enddo
    !         enddo
    !      enddo
    !
    !      if (p_pe == p_io) write(*,*) 'total iron input ',iron_sum

  END SUBROUTINE hamocc_get_dust
  ! ---------------------------------------------------------------------------

  ! ===========================================================================
  ! PRIVATE ROUTINES
  ! ===========================================================================

  ! ---------------------------------------------------------------------------
  SUBROUTINE hamocc_read_nml_cpl(status, iou)
    !
    ! read namelist for 'coupling' to ECHAM5
    !
    ! Author: Andrea Pozzer, MPICH, Aug 2003

    ! MESSy
    USE messy_main_tools,         ONLY: read_nml_open, read_nml_check         &
         &                            , read_nml_close
    USE messy_main_mpi_bi,        ONLY: finish

    IMPLICIT NONE

    ! I/O
    INTEGER, INTENT(OUT) :: status     ! error status
    INTEGER, INTENT(IN)  :: iou        ! I/O unit

    NAMELIST /CPL/ L_DUST, dust_input                                         &
         &       , L_SKIP_INT ! mz_bk_20100519+ skip integration              
    !                           passive tracer for testing mass balance

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr='hamocc_read_nml_cpl'
    LOGICAL                     :: lex      ! file exists ?
    INTEGER                     :: fstat    ! file status

    status = 1

    !----------------------------------------------------------------------
    ! DEFAULT PARAMETER SETTINGS - CAN BE OVERWRITEN BY THE NAMELIST
    L_DUST  = .TRUE.
    L_SKIP_INT = .FALSE.  ! mz_bk_20100519
    !----------------------------------------------------------------------

    CALL read_nml_open(lex, substr, iou, 'CPL', modstr)
    IF (.NOT.lex) RETURN    ! <modstr>.nml does not exist

    READ(iou, NML=CPL, IOSTAT=fstat)
    CALL read_nml_check(fstat, substr, iou, 'CPL', modstr)
    IF (fstat /= 0) RETURN  ! error while reading namelist

    CALL read_nml_close(substr, iou, modstr)
    status = 0 ! NO ERROR

  END SUBROUTINE hamocc_read_nml_cpl
  ! ---------------------------------------------------------------------------
  ! ===========================================================================

END MODULE messy_hamocc_e5

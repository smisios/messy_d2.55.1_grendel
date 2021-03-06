#include "messy_main_ppd_bi.inc"

MODULE messy_megan_si


  !-------------------------------------------------------------------------------------------------------
  !  MEGAN : Model of Emissions of Gases and Aerosols from Nature 
  !    
  !  AUTHOR:  Pozzer Andrea, MPICH, Feb 2011
  !
  !  MESSy Interface  for the Submodel kernel
  !-------------------------------------------------------------------------------------------------------
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
! 
!     Scientific algorithm 
! 
!             Emission = [EF][GAMMA][RHO] 
!           where [EF]    = emission factor (ug/m2h) 
!                 [GAMMA] = emission activity factor (non-dimension) 
!                 [RHO]   = production and loss within plant canopies 
!                           (non-dimensino) 
!                 Assumption: [RHO] = 1  (11/27/06)   (See PDT_LOT_CP.EXT) 
! 
!             GAMMA  = [GAMMA_CE][GAMMA_age][GAMMA_SM] 
!           where [GAMMA_CE]  = canopy correction factor 
!                 [GAMMA_age] = leaf age correction factor 
!                 [GAMMA_SM]  = soil moisture correction factor 
!                 Assumption: [GAMMA_SM]  = 1  (11/27/06) 
!             GAMMA_CE = [GAMMA_LAI][GAMMA_P][GAMMA_T] 
!           where [GAMMA_LAI] = leaf area index factor 
!                 [GAMMA_P]   = PPFD emission activity factor 
!                 [GAMMA_T]   = temperature response factor 
! 
!             Emission = [EF][GAMMA_LAI][GAMMA_P][GAMMA_T][GAMMA_age] 
!        Derivation: 
!             Emission = [EF][GAMMA_etc](1-LDF) + [EF][GAMMA_etc][LDF][GAMMA_P] 
!             Emission = [EF][GAMMA_etc]{ (1-LDF) + [LDF][GAMMA_P] } 
!             Emission = [EF][GAMMA_etc]{ (1-LDF) + [LDF][GAMMA_P] } 
!           where LDF = light dependent function (non-dimension) 
!                               (See LD_FCT.EXT) 
! 
!        Final Equation 
!             Emission = [EF][GAMMA_LAI][GAMMA_T][GAMMA_age]* 
!                        { (1-LDF) + [LDF][GAMMA_P] } 
! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

  USE messy_main_blather_bi,    ONLY: start_message_bi, end_message_bi
  USE messy_main_tools,         ONLY: PTR_2D_ARRAY,PTR_1D_ARRAY
  USE messy_megan
  USE messy_main_timer_bi,      ONLY: timer_event_init
  USE messy_main_timer_event,   ONLY: io_time_event, TRIG_LAST   &
                                    , time_event

  IMPLICIT NONE
  PRIVATE

  INTRINSIC NULL

  TYPE(PTR_2D_ARRAY), DIMENSION(:), POINTER, SAVE  :: flux_gp => NULL()  !final emissions
  TYPE(PTR_2D_ARRAY), DIMENSION(:), POINTER, SAVE  :: gam_smt => NULL()  ! soil moisture correction factor 
  TYPE(PTR_2D_ARRAY), DIMENSION(:), POINTER, SAVE  :: gam_tmp => NULL()  ! temperature response factor
  TYPE(PTR_2D_ARRAY), DIMENSION(:), POINTER, SAVE  :: gam_lht => NULL()  ! leaf area index factor 
  TYPE(PTR_2D_ARRAY), DIMENSION(:), POINTER, SAVE  :: gam_pho => NULL()  ! PPFD emission activity factor
  TYPE(PTR_2D_ARRAY), DIMENSION(:), POINTER, SAVE  :: gam_age => NULL()  ! leaf age correction factor
  REAL(dp)                                 , SAVE  :: rho                ! production and loss within plant canopies
  REAL(dp)                                 , SAVE  :: ldf                ! light dependent fraction  
  TYPE(PTR_2D_ARRAY), DIMENSION(:), POINTER, SAVE  :: gam_total => NULL()! total before EF multiplication
  ! COUPLING TO OTHER SUBMODEL
  ! op_pj_20120116+
  REAL(dp), POINTER, DIMENSION(:,:), PUBLIC :: cossza_2d => NULL()
  ! op_pj_20120116-
  REAL(dp), DIMENSION(:,:),POINTER   :: input_srfl => NULL() ! surface raditiation information 
  REAL(dp), DIMENSION(:,:),POINTER   :: input_lai  => NULL() ! leaf area index information
  REAL(dp), DIMENSION(:,:),POINTER   :: input_laip => NULL() ! leaf area index information (of previous month/timestep)
  REAL(dp), DIMENSION(:,:),POINTER   :: input_btr  => NULL() ! broadleaf coverage
  REAL(dp), DIMENSION(:,:),POINTER   :: input_ntr  => NULL() ! needleleaf coverage
  REAL(dp), DIMENSION(:,:),POINTER   :: input_shr  => NULL() ! shrub coverage
  REAL(dp), DIMENSION(:,:),POINTER   :: input_hrb  => NULL() ! herb/grass/crop coverage
!
  CHARACTER(LEN=STRLEN_MEDIUM),PUBLIC,SAVE :: sw_flux(2) = ''      !short wave flux at surface
  CHARACTER(LEN=STRLEN_MEDIUM),PUBLIC,SAVE :: lai(2) = ''          !leaf area index
  CHARACTER(LEN=STRLEN_MEDIUM),PUBLIC,SAVE :: laip(2) = ''         !leaf area index
  CHARACTER(LEN=STRLEN_MEDIUM),PUBLIC,SAVE :: btr_frac(2) = ''     ! broadleaf coverage
  CHARACTER(LEN=STRLEN_MEDIUM),PUBLIC,SAVE :: ntr_frac(2) = ''     ! needleleaf coverage
  CHARACTER(LEN=STRLEN_MEDIUM),PUBLIC,SAVE :: shr_frac(2) = ''     ! shrub coverage
  CHARACTER(LEN=STRLEN_MEDIUM),PUBLIC,SAVE :: hrb_frac(2) = ''     ! herb/grass/crop coverage
  ! op_pj_20120116+
  CHARACTER(LEN=STRLEN_MEDIUM),PUBLIC,SAVE :: cossza(2) = ''
  ! op_pj_20120116-

  INTEGER, PUBLIC, SAVE  :: nmegan = 0                  !actual number of basic tracer  
  INTEGER, PUBLIC, PARAMETER :: NMAXNTRAC_MEGAN =100        !total number of basic tracer interaction available 

  CHARACTER(LEN=STRLEN_MEDIUM), PUBLIC, SAVE       :: MGN_SPC(NMAXNTRAC_MEGAN)      !name of tracer
  REAL(dp),DIMENSION(NMAXNTRAC_MEGAN), PUBLIC,SAVE :: LDF_FCT     ! "light dependent" factor
  REAL(dp),DIMENSION(NMAXNTRAC_MEGAN), PUBLIC,SAVE :: RHO_FCT     ! production and loss within canopy factor
  REAL(dp),DIMENSION(NMAXNTRAC_MEGAN), PUBLIC,SAVE :: TDF_PRM     ! "temperature dependent" parameter for light-independent emissions
  REAL(dp),DIMENSION(NMAXNTRAC_MEGAN), PUBLIC,SAVE :: EF_BT       ! emission factor for broadleaf   
  REAL(dp),DIMENSION(NMAXNTRAC_MEGAN), PUBLIC,SAVE :: EF_NT       ! emission factor for needleleaf  
  REAL(dp),DIMENSION(NMAXNTRAC_MEGAN), PUBLIC,SAVE :: EF_SB       ! emission factor for shrub  
  REAL(dp),DIMENSION(NMAXNTRAC_MEGAN), PUBLIC,SAVE :: EF_HB       ! emission factor for herb/grass/crop


  ! MEGAN --> CHEM MECHANISM VARIABLE!!!
  INTEGER, PUBLIC, SAVE  :: nmech = 0                  !actual number of emissions  
  INTEGER, PUBLIC, PARAMETER :: NMAXNTRAC_MECH =200    !total number of tracer interaction available 

  CHARACTER(LEN=STRLEN_MEDIUM), PUBLIC, SAVE      :: MECH_SPC(NMAXNTRAC_MECH)      !name of tracer
  CHARACTER(LEN=STRLEN_MEDIUM), PUBLIC, SAVE      :: MGN_SRC(NMAXNTRAC_MECH)       ! MEGAN species "source" to be used
  REAL(dp),DIMENSION(NMAXNTRAC_MECH), PUBLIC,SAVE :: MECH_MWT    ! Molar weight of tracer
  REAL(dp),DIMENSION(NMAXNTRAC_MECH), PUBLIC,SAVE :: EFFS_BT     ! speciacion for broadleaf   
  REAL(dp),DIMENSION(NMAXNTRAC_MECH), PUBLIC,SAVE :: EFFS_NT     ! speciacion for needleleaf  
  REAL(dp),DIMENSION(NMAXNTRAC_MECH), PUBLIC,SAVE :: EFFS_SB     ! speciacion for shrub  
  REAL(dp),DIMENSION(NMAXNTRAC_MECH), PUBLIC,SAVE :: EFFS_HB     ! speciacion for herb/grass/crop
! mz_dc_20131118+
   REAL(dp),DIMENSION(NMAXNTRAC_MECH), PUBLIC,SAVE :: GS         ! Global scale factor for each tracer in the mechanism
! mz_dc_20131118-

  INTEGER, PUBLIC, SAVE  :: nout = 0                                               !actual number of emitted tracers (i.e. output)
  CHARACTER(LEN=STRLEN_MEDIUM), PUBLIC, SAVE       :: MGN_OUT(NMAXNTRAC_MECH)      !name of effective emitted tracers/ i.e. output
  
  ! TRACER INDEX (MECH -> MESSy TRACER)
  INTEGER, PUBLIC,  SAVE :: trac_id_gp(NMAXNTRAC_MECH)     !total number of tracer present in the program
  LOGICAL, PUBLIC,  SAVE :: l_trac_exist(NMAXNTRAC_MECH)
  

  REAL(dp), POINTER, SAVE :: PPFD(:,:)     ! calculated PPFD (Photosyntethic photon flux density) [umol/m2.s] (used in GAMMA_P)
                                             ! Note: should be similar to PAR (Photosyntethic active radiance) 
  REAL(dp), POINTER, SAVE :: D_PPFD(:,:)   ! DAILY AVERAGE PPFD (Photosyntethic photon flux density) [umol/m2.s] (used in GAMMA_P)
  REAL(dp), POINTER, SAVE :: TEMP(:,:)     ! Temperature (used in GAMMA_T)
  REAL(dp), POINTER, SAVE :: D_TEMP(:,:)   ! DAILY AVERAGE TEMPERATURE  (used in GAMMA_T)
  REAL(dp), POINTER, SAVE :: EF(:,:)       ! EMISSION FACTOR

  !ACCUMULATED FIELD!
  REAL(dp), POINTER       :: dt_megan
  REAL(dp), POINTER, SAVE :: acc_PPFD(:,:) ! accumulated PPFD (Photosyntethic photon flux density) [umol/m2.s] (used in GAMMA_P)
  REAL(dp), POINTER, SAVE :: acc_TEMP(:,:) ! accumulated TEMPERATURE  (used in GAMMA_T)



  !TIME MANAGMENT!
  TYPE(io_time_event), PUBLIC, SAVE :: trig_megan_daily = &
       io_time_event (1,'days',TRIG_LAST,0)
  TYPE(time_event),    PUBLIC, SAVE :: ev_trig_megan



  ! PUBLIC ECHAM-5 INTERFACE ROUTINES TO 'messy_SUBMODELS'
  PUBLIC :: megan_initialize
  PUBLIC :: megan_init_memory
  PUBLIC :: megan_init_coupling
  PUBLIC :: megan_global_end
  PUBLIC :: megan_vdiff
  PUBLIC :: megan_free_memory

CONTAINS

  ! ---------------------------------------------------------------------------

  SUBROUTINE megan_initialize

    ! MEGAN MODULE ROUTINE (ECHAM-5 INTERFACE)
    !
    ! INITIALIZATION OF GLOBAL VARIABLES FROM NAMELIST
    ! IN PARALLEL ENVIRONMENT
    !
 
    ! ECHAM5/MESSy
    USE messy_main_mpi_bi,     ONLY: p_parallel_io, p_io, p_bcast
    USE messy_main_blather_bi, ONLY: error_bi
    USE messy_main_tools,      ONLY: find_next_free_unit

    IMPLICIT NONE

    ! DEFINE NCREGRID EVENT TRIGGERs
    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr='megan_initialize'
    INTEGER                 :: status
    INTEGER                 :: iou    ! I/O unit
    INTEGER                 :: jt,jl
    LOGICAL                 :: l_exist

    CALL start_message_bi(modstr, 'MEMORY INITIALIZATION', substr)

    ! INITIALIZE CTRL
    IF (p_parallel_io) THEN
       iou = find_next_free_unit(100,200)
       ! *** CALL CORE ROUTINE:
       CALL megan_read_nml_ctrl(status, iou)
       IF (status /= 0) CALL error_bi(' ',substr)
    END IF

    ! BROADCAST RESULTS
    CALL p_bcast( l_tendency,      p_io)
    ! INITIALIZE CPL
    IF (p_parallel_io) THEN
       iou = find_next_free_unit(100,200)
       CALL megan_read_nml_cpl(status, iou)
       IF (status /= 0) CALL error_bi(' ',substr)
    END IF

    ! BROADCAST RESULTS
    ! CPL NAMELIST
    CALL p_bcast(sw_flux(1),  p_io)
    CALL p_bcast(sw_flux(2),  p_io)
    CALL p_bcast(lai(1),      p_io)
    CALL p_bcast(lai(2),      p_io)
    CALL p_bcast(laip(1),     p_io)
    CALL p_bcast(laip(2),     p_io)
    CALL p_bcast(btr_frac(1), p_io)
    CALL p_bcast(btr_frac(2), p_io)
    CALL p_bcast(ntr_frac(1), p_io)
    CALL p_bcast(ntr_frac(2), p_io)
    CALL p_bcast(shr_frac(1), p_io)
    CALL p_bcast(shr_frac(2), p_io)
    CALL p_bcast(hrb_frac(1), p_io)
    CALL p_bcast(hrb_frac(2), p_io)
    ! op_pj_20120116+
    CALL p_bcast(cossza(1),   p_io)
    CALL p_bcast(cossza(2),   p_io)
    ! op_pj_20120116-


    ! INITIALIZE MGN

    ! INIT
    DO jt=1, NMAXNTRAC_MEGAN
          MGN_SPC(jt)  = ''
          LDF_FCT(jt)  = 1.0_dp
          RHO_FCT(jt)  = 0.0_dp
          TDF_PRM(jt)  = 1.0_dp
          EF_BT(jt)    = 0.0_dp
          EF_NT(jt)    = 0.0_dp
          EF_SB(jt)    = 0.0_dp
          EF_HB(jt)    = 0.0_dp
    END DO

    ! MGN namelist
    IF (p_parallel_io) THEN
       iou = find_next_free_unit(100,200)
       CALL megan_read_nml_mgn(status, iou)
       IF (status /= 0) CALL error_bi(' ',substr)
       ! COUNT TRACERS
       nmegan = 0
       DO jt=1, NMAXNTRAC_MEGAN
          IF (TRIM(MGN_SPC(jt)) == '') CYCLE
          nmegan = nmegan + 1
          !
          MGN_SPC(nmegan) = TRIM(MGN_SPC(jt))
          LDF_FCT(nmegan) = LDF_FCT(jt)   
          RHO_FCT(nmegan) = RHO_FCT(jt)
          TDF_PRM(nmegan) = TDF_PRM(jt)
          EF_BT(nmegan)   = EF_BT(jt)
          EF_NT(nmegan)   = EF_NT(jt)
          EF_SB(nmegan)   = EF_SB(jt)
          EF_HB(nmegan)   = EF_HB(jt)
       END DO
    write (*,*) " Total number of tracer in megan_nml (20 originally)= " , nmegan 
    END IF

    CALL p_bcast(nmegan, p_io)
    DO jt=1, nmegan
       CALL p_bcast(MGN_SPC(jt), p_io)
       CALL p_bcast(LDF_FCT(jt), p_io)
       CALL p_bcast(RHO_FCT(jt), p_io)
       CALL p_bcast(TDF_PRM(jt), p_io)
       CALL p_bcast(EF_BT(jt),   p_io)
       CALL p_bcast(EF_NT(jt),   p_io)
       CALL p_bcast(EF_SB(jt),   p_io)
       CALL p_bcast(EF_HB(jt),   p_io)
    END DO


    ! INITIALIZE MGN2MECH

    !  INIT
    DO jt=1, NMAXNTRAC_MECH
          MECH_SPC(jt)  = ''
          MGN_SRC(jt)   = ''
          MECH_MWT(jt)  = 1.0_dp
          EFFS_BT(jt)   = 0.0_dp
          EFFS_NT(jt)   = 0.0_dp
          EFFS_SB(jt)   = 0.0_dp
          EFFS_HB(jt)   = 0.0_dp
          GS(jt)        = 1.0_dp
    END DO

    ! MGN2MECH namelist
    IF (p_parallel_io) THEN
       iou = find_next_free_unit(100,200)
       CALL megan_read_nml_mgn2mech(status, iou)
       IF (status /= 0) CALL error_bi(' ',substr)
       ! COUNT TRACERS
       nmech = 0
       DO jt=1, NMAXNTRAC_MECH
          IF (TRIM(MECH_SPC(jt)) == '') CYCLE
          nmech = nmech + 1
          !
          MECH_SPC(nmech) = TRIM(MECH_SPC(jt))
          MGN_SRC(nmech) = TRIM(MGN_SRC(jt))   
          MECH_MWT(nmech)  = MECH_MWT(jt)
          EFFS_BT(nmech)   = EFFS_BT(jt)
          EFFS_NT(nmech)   = EFFS_NT(jt)
          EFFS_SB(nmech)   = EFFS_SB(jt)
          EFFS_HB(nmech)   = EFFS_HB(jt)
          GS(nmech)        = GS(jt)
       END DO
    write (*,*) " Total number of emissions via MEGAN= " , nmech 
    END IF


    CALL p_bcast(nmech,  p_io)
    DO jt=1, nmech
       CALL p_bcast(MECH_SPC(jt), p_io)
       CALL p_bcast(MGN_SRC(jt),  p_io)   
       CALL p_bcast(MECH_MWT(jt), p_io)
       CALL p_bcast(EFFS_BT(jt),  p_io)
       CALL p_bcast(EFFS_NT(jt),  p_io)
       CALL p_bcast(EFFS_SB(jt),  p_io)
       CALL p_bcast(EFFS_HB(jt),  p_io)
       CALL p_bcast(GS(jt),       p_io) 
    END DO

    ! calculate single emitted species
    IF (p_parallel_io) THEN
       ! COUNT TRACERS
       nout = 0
       DO jt=1, nmech
          l_exist = .FALSE.
          DO jl=1,nout
             IF (TRIM(MECH_SPC(jt)) == TRIM(MGN_OUT(jl))) l_exist = .TRUE. 
          ENDDO
          IF (.NOT.l_exist) THEN 
             nout = nout + 1
             MGN_OUT(jl) = TRIM(MECH_SPC(jt))
          ENDIF
       END DO
    write (*,*) " Number of tracer emissions estimated via MEGAN= " , nout
    END IF

    CALL p_bcast(nout,  p_io)
    DO jt=1, nout
       CALL p_bcast(MGN_OUT(jt),  p_io)   
    END DO

    ! TIME EVENT INITIALISATION FOR DAILY AVERAGES
    CALL timer_event_init(ev_trig_megan, trig_megan_daily, &
           'megan computation', 'present') 

  END SUBROUTINE megan_initialize 

  ! ---------------------------------------------------------------------------

  SUBROUTINE megan_init_memory
   
    ! ECHAM5/MESSy
    USE messy_main_tracer_mem_bi   , ONLY: ntrac_gp, ti_gp     
    USE messy_main_grid_def_mem_bi,  ONLY: ngpblks, nproma
    ! MESSy
    USE messy_main_channel,          ONLY: new_channel, new_channel_object &
                                         , new_attribute
    USE messy_main_channel_error_bi, ONLY: channel_halt
    USE messy_main_channel_bi,       ONLY: GP_2D_HORIZONTAL, SCALAR
    USE messy_main_mpi_bi,           ONLY: p_parallel_io
    USE messy_main_blather_bi,       ONLY: error_bi

    IMPLICIT NONE
    INTRINSIC TRIM


    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'megan_init_memory'
    INTEGER :: status
    INTEGER :: jt,jl
    LOGICAl :: l_exist_src
    CHARACTER(LEN=4) :: echar = '    '

    CALL start_message_bi(modstr, 'MEMORY INITIALIZATION ', substr)

    ALLOCATE(gam_smt(nmegan))
    ALLOCATE(gam_tmp(nmegan))
    ALLOCATE(gam_lht(nmegan))
    ALLOCATE(gam_pho(nmegan))
    ALLOCATE(gam_age(nmegan))
    ALLOCATE(gam_total(nmegan))

! output fluxes
    ALLOCATE(flux_gp(nout))

    ! LOCAL variables
    ALLOCATE(PPFD(nproma,ngpblks))
    ALLOCATE(TEMP(nproma,ngpblks))
    ALLOCATE(EF(nproma,ngpblks))

    CALL new_channel(status,modstr//'_gp', reprid=GP_2D_HORIZONTAL)
    CALL channel_halt(substr, status)  

     DO jt=1,nout 
        CALL new_channel_object(status, modstr//'_gp'     &
             , TRIM(MGN_OUT(jt))//'_flux'                &
             , p2=flux_gp(jt)%ptr)                            
        CALL channel_halt(substr, status)
        CALL new_attribute(status, modstr//'_gp'                &
             , TRIM(MGN_OUT(jt))//'_flux'                      &
             , 'longname', c='flux for '//TRIM(MGN_OUT(jt)) )
        CALL channel_halt(substr, status)
        CALL new_attribute(status, modstr//'_gp'                &
             , TRIM(MGN_OUT(jt))//'_flux'                      &
             , 'units', c='mol / (m^2 s)')
        CALL channel_halt(substr, status)
     END DO

     DO jt=1,nmegan 
        CALL new_channel_object(status, modstr//'_gp'      &
             , TRIM(MGN_SPC(jt))//'_gam_smt'               &
             , p2=gam_smt(jt)%ptr)                      
        CALL channel_halt(substr, status) 
        CALL new_attribute(status, modstr//'_gp'           &
             , TRIM(MGN_SPC(jt))//'_gam_smt'               &
             , 'longname', c='soil mosture factor for '//TRIM(MGN_SPC(jt))) 
        CALL channel_halt(substr, status) 

        CALL new_channel_object(status, modstr//'_gp'      &
             , TRIM(MGN_SPC(jt))//'_gam_tmp'               &
             , p2=gam_tmp(jt)%ptr)
        CALL channel_halt(substr, status) 
        CALL new_attribute(status, modstr//'_gp'           &
             , TRIM(MGN_SPC(jt))//'_gam_tmp'               &
             , 'longname', c='temperature factor for '//TRIM(MGN_SPC(jt)))
        CALL channel_halt(substr, status) 

        CALL new_channel_object(status, modstr//'_gp'      &
             , TRIM(MGN_SPC(jt))//'_gam_lht'               &
             , p2=gam_lht(jt)%ptr)
        CALL channel_halt(substr, status) 
        CALL new_attribute(status, modstr//'_gp'           &
             , TRIM(MGN_SPC(jt))//'_gam_lht'               &
             , 'longname', c='LAI factor for '//TRIM(MGN_SPC(jt)))
        CALL channel_halt(substr, status) 

        CALL new_channel_object(status, modstr//'_gp'      &
             , TRIM(MGN_SPC(jt))//'_gam_pho'               &
             , p2=gam_pho(jt)%ptr)
        CALL channel_halt(substr, status) 
        CALL new_attribute(status, modstr//'_gp'           &
             , TRIM(MGN_SPC(jt))//'_gam_pho'               &
             , 'longname', c='PPFD factor for '//TRIM(MGN_SPC(jt)))
        CALL channel_halt(substr, status) 

        CALL new_channel_object(status, modstr//'_gp'      &
             , TRIM(MGN_SPC(jt))//'_gam_age'               &
             , p2=gam_age(jt)%ptr)
        CALL channel_halt(substr, status) 
        CALL new_attribute(status, modstr//'_gp'           &
             , TRIM(MGN_SPC(jt))//'_gam_age'               &
             , 'longname', c='Leaf age factor for '//TRIM(MGN_SPC(jt)))
        CALL channel_halt(substr, status) 

        CALL new_channel_object(status, modstr//'_gp'      &
             , TRIM(MGN_SPC(jt))//'_gam_total'               &
             , p2=gam_total(jt)%ptr)
        CALL channel_halt(substr, status) 
        CALL new_attribute(status, modstr//'_gp'           &
             , TRIM(MGN_SPC(jt))//'_gam_total'               &
             , 'longname', c='total gamma factor for '//TRIM(MGN_SPC(jt)))
        CALL channel_halt(substr, status) 
     END DO

     !ALLOCATE FIELDS FOR TEMPORARY STORING
     CALL new_channel_object(status, modstr//'_gp'      &
          , 'acc_PPFD',  p2=acc_PPFD, lrestreq=.TRUE.)          
     CALL channel_halt(substr, status) 
     CALL new_attribute(status, modstr//'_gp'           &
          , 'acc_PPFD', 'longname', c='accumulated PPFD ')
     CALL channel_halt(substr, status) 

     !ALLOCATE FIELDS FOR TEMPORARY STORING
     CALL new_channel_object(status, modstr//'_gp'      &
          , 'acc_TEMP',  p2=acc_TEMP, lrestreq=.TRUE.)
     CALL channel_halt(substr, status) 
     CALL new_attribute(status, modstr//'_gp'           &
          , 'acc_TEMP', 'longname', c='accumulated TEMP ')
     CALL channel_halt(substr, status) 

     !ALLOCATE FIELDS FOR TEMPORARY STORING
     CALL new_channel_object(status, modstr//'_gp'      &
          , 'dt_megan',  p0=dt_megan, REPRID=SCALAR, lrestreq=.TRUE.)                                   
     CALL channel_halt(substr, status) 
     CALL new_attribute(status, modstr//'_gp'           &
          , 'dt_megan', 'longname', c='dt megan ')
     CALL channel_halt(substr, status) 

     !ALLOCATE AVERAGED FIELDS 
     CALL new_channel_object(status, modstr//'_gp'      &
          , 'D_PPFD',  p2=D_PPFD, lrestreq=.TRUE.)          
     CALL channel_halt(substr, status) 
     CALL new_attribute(status, modstr//'_gp'           &
          , 'D_PPFD', 'longname', c='previous day average PPFD ')
     CALL channel_halt(substr, status) 

     !ALLOCATE AVERAGED FIELDS 
     CALL new_channel_object(status, modstr//'_gp'      &
          , 'D_TEMP',  p2=D_TEMP, lrestreq=.TRUE.)          
     CALL channel_halt(substr, status) 
     CALL new_attribute(status, modstr//'_gp'           &
          , 'D_TEMP', 'longname', c='previous day average temperature ')
     CALL channel_halt(substr, status) 
     CALL new_attribute(status, modstr//'_gp'           &
          , 'D_TEMP', 'units', c='K')
     CALL channel_halt(substr, status) 

!CHECK FOR TRACER EXISTANCE
      DO jl = 1,nout
        l_trac_exist(jl) = .FALSE.
        DO jt = 1,ntrac_gp
         IF(ti_gp(jt)%tp%ident%basename == MGN_OUT(jl)) THEN
           trac_id_gp(jl) = jt
           l_trac_exist(jl) = .TRUE.
           IF (p_parallel_io) THEN
               WRITE(*,*) ' TRACER ',MGN_OUT(jl), ' exist! '
               WRITE(*,*) ' TRACER UPDATED WITH MEGAN EMISSIONS'
           END IF
         ENDIF
        END DO
        IF (.not.l_trac_exist(jl)) THEN 
          IF (p_parallel_io) THEN
               WRITE(*,*) ' !!! -> TRACER ',MGN_OUT(jl), ' DO NOT exist! '
               WRITE(*,*) ' ONLY DIAGNOSTIC IN MEGAN CHANNEL'
          END IF
         ENDIF
      END DO
   
      !  loop over chemical mechanism
      DO jl= 1, nmech
        l_exist_src = .FALSE.
        !loop over megan species
        DO jt= 1, nmegan
           IF (TRIM(MGN_SRC(jl))==TRIM(MGN_SPC(jt))) l_exist_src=.TRUE.
        END DO
        IF (.not.l_exist_src) CALL error_bi(TRIM(MGN_SRC(jl))//' as source tracer do not exist',substr)
      END DO      

 
    CALL end_message_bi(modstr, 'MEMORY INITIALIZATION ', substr)

  END SUBROUTINE megan_init_memory

  ! ---------------------------------------------------------------------------

  SUBROUTINE megan_init_coupling

    ! MESSy
    USE messy_main_channel,          ONLY: get_channel_object, get_channel_info &
                                         , get_channel_object_info
    USE messy_main_channel_error_bi, ONLY: channel_halt 
    USE messy_main_mpi_bi,        ONLY: p_parallel_io
    USE messy_main_blather_bi,    ONLY: error_bi, info_bi

    USE messy_main_tracer_mem_bi, ONLY: ti_gp
    USE messy_main_tracer,        ONLY: R_molarmass

    IMPLICIT NONE

    CHARACTER(LEN=*), PARAMETER :: substr = 'megan_init_coupling'
    ! I/O
    INTEGER :: status
    INTEGER :: jl
    INTEGER :: reprid

    CALL start_message_bi(modstr, 'COUPLING', substr)

    !CALL message(substr, 'IMPORT DATA: ')
    CALL info_bi('IMPORT DATA: ', substr)
    call get_channel_info(status, sw_flux(1))
    IF (status /= 0) THEN
       IF (p_parallel_io) THEN
           WRITE(*,*) ' no radiation at surface avaliable  '
           WRITE(*,*) ' channel not present ! (check namelist)'
!           CALL FINISH(substr,'surface radiation not present!')       
           CALL error_bi('surface radiation not present!', substr)       
       END IF
    ELSE
       call get_channel_object_info(status, cname=sw_flux(1), oname=sw_flux(2), reprid=reprid )        
       IF (status /= 0) THEN
          IF (p_parallel_io) THEN
              WRITE(*,*) ' no radiation at surface available '
              WRITE(*,*) ' object not present ! (check namelist)'
!              CALL FINISH(substr,'surface radiation not present!')       
              CALL error_bi('surface radiation not present!', substr)       
          END IF
       ELSE
          IF (p_parallel_io) THEN
              WRITE(*,*) ' Radiation at surface available from channel ', sw_flux(2) 
          END IF
          CALL get_channel_object(status, cname=sw_flux(1), oname=sw_flux(2) &
               , p2=input_srfl )
          CALL channel_halt(substr, status)
       END IF
    END IF

    call get_channel_info(status, lai(1))
    IF (status /= 0) THEN
       IF (p_parallel_io) THEN
           WRITE(*,*) ' no Leaf Area Index (LAI) avaliable  '
           WRITE(*,*) ' channel not present ! (check namelist)'
!           CALL FINISH(substr,'Leaf Area Index (LAI) not present!')       
           CALL error_bi('Leaf Area Index (LAI) not present!', substr)
       END IF
    ELSE
       call get_channel_object_info(status, cname=lai(1), oname=lai(2), reprid=reprid)
       IF (status /= 0) THEN
          IF (p_parallel_io) THEN
              WRITE(*,*) ' no Leaf Area Index (LAI) avaliable  '
              WRITE(*,*) ' object not present ! (check namelist)'
!              CALL FINISH(substr,'Leaf Area Index (LAI) not present!')       
              CALL error_bi('Leaf Area Index (LAI) not present!', substr)
          END IF
       ELSE
          IF (p_parallel_io) THEN
              WRITE(*,*) ' Leaf Area Index (LAI) avaliable from channel ', lai(2) 
          END IF
          CALL get_channel_object(status, cname=lai(1), oname=lai(2) &
               , p2=input_lai )
          CALL channel_halt(substr, status)
       END IF
    END IF

    call get_channel_info(status, laip(1))
    IF (status /= 0) THEN
       IF (p_parallel_io) THEN
           WRITE(*,*) ' no Leaf Area Index (LAIp) avaliable  '
           WRITE(*,*) ' channel not present ! (check namelist)'
!           CALL FINISH(substr,'Leaf Area Index (LAIp) not present!')       
           CALL error_bi('Leaf Area Index (LAIp) not present!', substr)
       END IF
    ELSE
       call get_channel_object_info(status, cname=laip(1), oname=laip(2), reprid=reprid)
       IF (status /= 0) THEN
          IF (p_parallel_io) THEN
              WRITE(*,*) ' no Leaf Area Index (LAIp) avaliable  '
              WRITE(*,*) ' object not present ! (check namelist)'
!              CALL FINISH(substr,'Leaf Area Index (LAIp) not present!')       
              CALL error_bi('Leaf Area Index (LAIp) not present!', substr)
          END IF
       ELSE
          IF (p_parallel_io) THEN
              WRITE(*,*) ' Previous  Leaf Area Index (LAIp) avaliable from channel ', laip(2) 
          END IF
          CALL get_channel_object(status, cname=laip(1), oname=laip(2) &
               , p2=input_laip )
          CALL channel_halt(substr, status)
       END IF
    END IF

    call get_channel_info(status, btr_frac(1))
    IF (status /= 0) THEN
       IF (p_parallel_io) THEN
           WRITE(*,*) ' no broadleaf coverage avaliable  '
           WRITE(*,*) ' channel not present ! (check namelist)'
!           CALL FINISH(substr,'broadleaf coverage not present!')       
           CALL error_bi('broadleaf coverage not present!',substr)
       END IF
    ELSE
       call get_channel_object_info(status, cname=btr_frac(1), oname=btr_frac(2), reprid=reprid)
       IF (status /= 0) THEN
          IF (p_parallel_io) THEN
              WRITE(*,*) ' no broadleaf coverage avaliable  '
              WRITE(*,*) ' object not present ! (check namelist)'
!              CALL FINISH(substr,'broadleaf coverage not present!')       
              CALL error_bi('broadleaf coverage not present!',substr)
          END IF
       ELSE
          IF (p_parallel_io) THEN
              WRITE(*,*) ' broadleaf coverage avaliable from channel ', btr_frac(2) 
          END IF
          CALL get_channel_object(status, cname=btr_frac(1), oname=btr_frac(2) &
               , p2=input_btr )
          CALL channel_halt(substr, status)
       END IF
    END IF

    call get_channel_info(status, ntr_frac(1))
    IF (status /= 0) THEN
       IF (p_parallel_io) THEN
           WRITE(*,*) ' no needleleaf coverage avaliable  '
           WRITE(*,*) ' channel not present ! (check namelist)'
!           CALL FINISH(substr,'needleleaf coverage not present!')       
           CALL error_bi('needleleaf coverage not present!',substr)
       END IF
    ELSE
       call get_channel_object_info(status, cname=ntr_frac(1), oname=ntr_frac(2), reprid=reprid)
       IF (status /= 0) THEN
          IF (p_parallel_io) THEN
              WRITE(*,*) ' no needleleaf coverage avaliable  '
              WRITE(*,*) ' object not present ! (check namelist)'
!              CALL FINISH(substr,'needleleaf coverage not present!')       
              CALL error_bi('needleleaf coverage not present!',substr)
          END IF
       ELSE
          IF (p_parallel_io) THEN
              WRITE(*,*) ' needleleaf coverage avaliable from channel ', ntr_frac(2) 
          END IF
          CALL get_channel_object(status, cname=ntr_frac(1), oname=ntr_frac(2) &
               , p2=input_ntr )
          CALL channel_halt(substr, status)
       END IF
    END IF

    call get_channel_info(status, shr_frac(1))
    IF (status /= 0) THEN
       IF (p_parallel_io) THEN
           WRITE(*,*) ' no shrub coverage avaliable  '
           WRITE(*,*) ' channel not present ! (check namelist)'
!           CALL FINISH(substr,'shrub coverage not present!')       
           CALL error_bi('shrub coverage not present!',substr)
       END IF
    ELSE
       call get_channel_object_info(status, cname=shr_frac(1), oname=shr_frac(2), reprid=reprid)
       IF (status /= 0) THEN
          IF (p_parallel_io) THEN
              WRITE(*,*) ' no shrub coverage avaliable  '
              WRITE(*,*) ' object not present ! (check namelist)'
!              CALL FINISH(substr,'shrub coverage not present!')       
              CALL error_bi('shrub coverage not present!',substr)
          END IF
       ELSE
          IF (p_parallel_io) THEN
              WRITE(*,*) ' shrub coverage avaliable from channel ', shr_frac(2) 
          END IF
          CALL get_channel_object(status, cname=shr_frac(1), oname=shr_frac(2) &
               , p2=input_shr )
          CALL channel_halt(substr, status)
       END IF
    END IF

    call get_channel_info(status, hrb_frac(1))
    IF (status /= 0) THEN
       IF (p_parallel_io) THEN
           WRITE(*,*) ' no herb/grass/crop coverage avaliable  '
           WRITE(*,*) ' channel not present ! (check namelist)'
!           CALL FINISH(substr,'herb/grass/crop coverage not present!')       
           CALL error_bi('herb/grass/crop coverage not present!',substr)
       END IF
    ELSE
       call get_channel_object_info(status, cname=hrb_frac(1), oname=hrb_frac(2), reprid=reprid)
       IF (status /= 0) THEN
          IF (p_parallel_io) THEN
              WRITE(*,*) ' no herb/grass/crop coverage avaliable  '
              WRITE(*,*) ' object not present ! (check namelist)'
!              CALL FINISH(substr,'herb/grass/crop coverage not present!')       
              CALL error_bi('herb/grass/crop coverage not present!',substr)
          END IF
       ELSE
          IF (p_parallel_io) THEN
              WRITE(*,*) ' herb/grass/crop coverage avaliable from channel ', hrb_frac(2) 
          END IF
          CALL get_channel_object(status, cname=hrb_frac(1), oname=hrb_frac(2) &
               , p2=input_hrb )
          CALL channel_halt(substr, status)
       END IF
    END IF

    ! op_pj_20120116+
    CALL info_bi('Looking for cos(solar zenith angle) (cossza)')
    CALL info_bi('       channel: '//cossza(1))
    CALL info_bi('       object : '//cossza(2))

    CALL get_channel_object(status  &
         ,  TRIM(cossza(1)), TRIM(cossza(2)), p2=cossza_2d)
    CALL channel_halt(substr, status)    
    ! op_pj_20120116-

    CALL end_message_bi(modstr, 'COUPLING', substr)

  END SUBROUTINE megan_init_coupling

  ! ---------------------------------------------------------------------------

  SUBROUTINE megan_vdiff
 
    USE messy_main_blather_bi,      ONLY: error_bi
    USE messy_main_grid_def_mem_bi, ONLY: jrow, kproma, nlev
    USE messy_main_grid_def_bi,     ONLY: grmass, grvol, deltaz
    USE messy_main_data_bi,         ONLY: &
! mz_rj_20140801+
#if defined (ECHAM5)
                                        pxtems,                       &
#endif
! mz_rj_20140801-
                                        tslm1, alake
    USE messy_main_tracer_mem_bi,   ONLY: pxtte => qxtte 
    USE messy_main_timer,         ONLY: JulianMonthLength
    USE messy_main_constants_mem, ONLY: g, M_air
    USE messy_main_timer,         ONLY: lstart, DAY, MONTH, YEAR

    IMPLICIT NONE

    CHARACTER(LEN=*), PARAMETER :: substr='megan_vdiff'
    REAL(dp), PARAMETER :: hr2sec = 3600     ! convert hr to second
    REAL(dp), PARAMETER :: ug2g = 1E-6       ! convert microgram to gram
    INTEGER :: jt,jc,jl,idt ! tracer loop
    INTEGER :: nday_month ! number of days in this month

    ! air mol / kg air
    REAL(DP), PARAMETER                   :: uconv = 1.0e3_DP/M_air
    REAL(DP), DIMENSION(:), ALLOCATABLE :: zairdens  ! air density
    REAL(DP), DIMENSION(:), ALLOCATABLE :: zdz       ! layer thickness
    REAL(DP), DIMENSION(:), ALLOCATABLE :: zscale    !
! to be cleand.... possibly not used
!    REAL(DP), DIMENSION(:), ALLOCATABLE :: land_mask !

    IF (lstart) THEN
      D_PPFD(1:kproma,jrow) = 0.0_dp
      D_TEMP(1:kproma,jrow) = 273.0_dp
    ENDIF

    ! calculation of PPFD (Photosyntethic photon flux density)
    !        PPFD: input_srfl - short wave from sun (W/m2)  !TODO: CHECK IF CORRECT INPUT OBJECT
    !        assuming 4.766 (umol m-2 s-1) per (W m-2)
    !        assume 1/2 of SRAD is in 400-700nm band
    PPFD(1:kproma,jrow)=input_srfl(1:kproma,jrow) * 4.766 * 0.5
 
    ! Temperature at the surface
    TEMP(1:kproma,jrow)=tslm1(1:kproma,jrow) !surface temperature 
 
    ! CALCULATE AIR DENSITY, LAYER THICKNESS AT SURFACE AND LAND MASK
    ALLOCATE(zairdens(kproma))
    ALLOCATE(zdz(kproma))
    ALLOCATE(zscale(kproma))

    zairdens(1:kproma) = grmass(_RI_XYZ__(1:kproma,jrow,nlev))/grvol(_RI_XYZ__(1:kproma,jrow,nlev))
    zdz(1:kproma) = deltaz(_RI_XYZ__(1:kproma,jrow,nlev))
    zscale(1:kproma) = zairdens(1:kproma)*zdz(1:kproma)

    DO jt=1,nmegan
       !----------------------------------------------------
       !  TEMPERATURE DEPENDENCE FACTOR (GAMMA_T)--> GAM_TMP
       !----------------------------------------------------
       SELECT CASE (MGN_SPC(jt))
       CASE ('ISOP')
          CALL GAMMA_TISOP( TEMP(1:kproma,jrow), D_TEMP(1:kproma,jrow), gam_tmp(jt)%ptr(1:kproma,jrow) )
       CASE ( 'MBO','MYRC','SABI','LIMO','CAR3','OCIM','BPIN',        &
             'APIN','FARN','BCAR','MEOH','ACTO','ACTA','FORM',        &
              'CH4',  'NO','OMTP','OSQT',  'CO'                )
          CALL GAMMA_TNISP( TDF_PRM(jt), TEMP(1:kproma,jrow), gam_tmp(jt)%ptr(1:kproma,jrow))
       CASE DEFAULT
!          CALL FINISH(substr,'Error: GAMMA_T, invalid variable:'//MGN_SPC(jt))       
          CALL error_bi('Error: GAMMA_T, invalid variable:'//MGN_SPC(jt),substr)       
       ENDSELECT
       !----------------------------------------------------
       !  LAI CORRECTION FACTOR (GAMMA_LAI) --> GAM_LHT
       !----------------------------------------------------
       CALL GAMMA_LAI(input_lai(1:kproma,jrow),gam_lht(jt)%ptr(1:kproma,jrow))
       !----------------------------------------------------
       !  LIGHT CORRECTION FACTOR (GAMMA_P) --> GAM_PHO
       !----------------------------------------------------
       ! NB: cossza_2d : COS(zenith angle)
       CALL GAMMA_P(DAY, cossza_2d(1:kproma,jrow), PPFD(1:kproma,jrow), D_PPFD(1:kproma,jrow), &
                    gam_pho(jt)%ptr(1:kproma,jrow) )
       !----------------------------------------------------
       !  AGE CORRECTION FACTOR (GAMMA_A) --> GAM_AGE
       !----------------------------------------------------
       nday_month = JulianMonthLength(YEAR, MONTH)
       CALL GAMMA_A(MGN_SPC(jt),input_laip(1:kproma,jrow),input_lai(1:kproma,jrow),nday_month, &
                    D_TEMP(1:kproma,jrow), gam_age(jt)%ptr(1:kproma,jrow) )
       !----------------------------------------------------
       !  SOIL MOSTURE DEPENDENCE FACTOR (GAMMA_S) --> GAM_SMT
       !----------------------------------------------------
        CALL GAMMA_S(gam_smt(jt)%ptr(1:kproma,jrow))
       !----------------------------------------------------
       ! ADDITIONAL TERMS IN THE EQUATIONS (FROM NAMELIST)
       !----------------------------------------------------
       !----------------------------------------------------
       !  PRODUCTION AND LOSSES WITHIN THE CANOPY (RHO_FCT)--> RHO
       !----------------------------------------------------
       ! in this case is always equal to 1 (see namelist) 
       RHO = RHO_FCT(jt)
       !----------------------------------------------------
       !  light dependent fraction (LDF_FCT)--> LDF
       !----------------------------------------------------
       LDF = LDF_FCT(jt)
       ! final emissions ER = ( EF * GAM_TMP * GAM_AGE * GAM_LHT * GAM_SMT * RHO ) *  &
       !                       ( (1-LDF) + (GAM_PHO*LDF) ) 
       ! which is equivalent to : ER =  EF * (GAM_TMP * GAM_AGE * GAM_LHT * GAM_SMT * RHO ) *  &
       !                          ( (1-LDF) + (GAM_PHO*LDF) )
       ! which is equivalent to : ER =  EF * (GAM_TOTAL)   
       ! with GAM_TOTAL = (GAM_TMP * GAM_AGE * GAM_LHT * GAM_SMT * RHO ) * ( (1-LDF) + (GAM_PHO*LDF) )
       gam_total(jt)%ptr(1:kproma,jrow) = ( gam_tmp(jt)%ptr(1:kproma,jrow)                       & 
                                          * gam_age(jt)%ptr(1:kproma,jrow)                       &  
                                          * gam_lht(jt)%ptr(1:kproma,jrow)                       & 
                                          * gam_smt(jt)%ptr(1:kproma,jrow) * RHO )               & 
                                          * ((1-LDF)+(gam_pho(jt)%ptr(1:kproma,jrow) *LDF)) 
       WHERE (gam_total(jt)%ptr(1:kproma,jrow) < 0.0_dp) gam_total(jt)%ptr(1:kproma,jrow) = 0.0_dp
    END DO !jt loop (nmegan)

  ! SPECIATION TO ACTUAL CHEMICAL MECHANISM  !

   !loop over total emissions outputted
   DO jl=1,nout
      ! nullify any total flux calculated before
      flux_gp(jl)%ptr(1:kproma,jrow) = 0.0_dp
      !  loop over chemical mechanism
      DO jc=1,nmech
         IF (TRIM(MECH_SPC(jc)) == TRIM(MGN_OUT(jl))) THEN
            !loop over megan species
            DO jt=1,nmegan
               IF (TRIM(MGN_SRC(jc))==TRIM(MGN_SPC(jt))) THEN
               !----------------------------------------------------
               !  Estimation of emissions factor (EF) based on the PFT (Plant Functional Type)
               !----------------------------------------------------
         !     ORIGINAL PART OF MEGAN: WE DO NOT HAVE % OF COVERAGE BUT FRACTIONS (SEE IMPORT!).
         !     EF(1:kproma,jrow) =(( EF_BT(jt)*EFFS_BT(jc)*input_btr(1:kproma,jrow)/100)  + &
         !                         ( EF_NT(jt)*EFFS_NT(jc)*input_ntr(1:kproma,jrow)/100)  + &
         !                         ( EF_SB(jt)*EFFS_SB(jc)*input_shr(1:kproma,jrow)/100)  + &
         !                         ( EF_HB(jt)*EFFS_HB(jc)*input_hrb(1:kproma,jrow)/100))
         !mz_dc_20131118+ 
               EF(1:kproma,jrow) =(( EF_BT(jt)*EFFS_BT(jc)*input_btr(1:kproma,jrow))  + &
                                   ( EF_NT(jt)*EFFS_NT(jc)*input_ntr(1:kproma,jrow))  + &
                                   ( EF_SB(jt)*EFFS_SB(jc)*input_shr(1:kproma,jrow))  + &
                                   ( EF_HB(jt)*EFFS_HB(jc)*input_hrb(1:kproma,jrow)))

        !mz_dc_20131118- 
               !----------------------------------------------------
               !  FLUX CALCULATION in mol/m^2s
               !----------------------------------------------------
               ! global scale factor (GS) for final global tuning of emissions
               ! 
                    flux_gp(jl)%ptr(1:kproma,jrow) = flux_gp(jl)%ptr(1:kproma,jrow)+ &
                                                     ( EF(1:kproma,jrow) * gam_total(jt)%ptr(1:kproma,jrow)) &
                                                     *((1/hr2sec)*ug2g / MECH_MWT(jc) * GS(jc))
               END IF
            END DO
         ENDIF
      END DO
      !----------------------------------------------------
      !             UPGRADE CONCENTRATION.... ONLY IN GRID POINT!
      !----------------------------------------------------
      IF (l_trac_exist(jl)) THEN 
#ifdef ECHAM5
         IF (.NOT. l_tendency ) THEN
            ! mol/m^2/s -> mol(X)Kg(air)/mol(air)m^2s
           pxtems(1:kproma,1,trac_id_gp(jl),jrow) = pxtems(1:kproma,1,trac_id_gp(jl),jrow) + &
                                             flux_gp(jl)%ptr(1:kproma,jrow)/uconv
         ELSE  
#endif
           ! mol/m^2/s -> mol/mol s
            ! um_ak_20110725+
 !          pxtte(1:kproma,nlev,trac_id_gp(jc)) = pxtte(1:kproma,nlev,trac_id_gp(jc)) +       & 
            idt = trac_id_gp(jl)
            pxtte(_RI_X_ZN_(1:kproma,nlev,idt)) = pxtte(_RI_X_ZN_(1:kproma,nlev,idt)) +       & 
            ! um_ak_20110725-
! op_pj_20180713+
!!$              flux_gp(jc)%ptr(1:kproma,jrow)/(zscale(1:kproma)*uconv)
                 flux_gp(jl)%ptr(1:kproma,jrow)/(zscale(1:kproma)*uconv)
! op_pj_20180713-
#ifdef ECHAM5
         END IF
#endif
     END IF
   END DO

   ! CLEAN UP
   DEALLOCATE(zairdens)
   DEALLOCATE(zdz)
   DEALLOCATE(zscale)

  END SUBROUTINE megan_vdiff

  ! ---------------------------------------------------------------------------

  SUBROUTINE megan_global_end

    USE messy_main_timer_bi,      ONLY: event_state
    USE messy_main_timer,         ONLY: current_date,lstart,delta_time

    IMPLICIT NONE

    CHARACTER(LEN=*), PARAMETER :: substr='megan_global_end'

   IF (lstart) THEN
     acc_PPFD(:,:) = 0.0_dp
     acc_TEMP(:,:) = 0.0_dp
     dt_megan = 0.0_dp
   ENDIF
! We need daily average temperature and radiation
   IF (event_state(ev_trig_megan,current_date)) THEN
     D_PPFD(:,:) = acc_PPFD(:,:)/dt_megan
     D_TEMP(:,:) = acc_TEMP(:,:)/dt_megan
     ! reset to zero !
     acc_PPFD(:,:) = 0.0_dp
     acc_TEMP(:,:) = 0.0_dp
     dt_megan = 0.0_dp
   ELSE
     acc_PPFD(:,:) = acc_PPFD(:,:) + PPFD(:,:)*delta_time
     acc_TEMP(:,:) = acc_TEMP(:,:) + TEMP(:,:)*delta_time
     dt_megan = dt_megan + delta_time
   ENDIF

  END SUBROUTINE

  ! ---------------------------------------------------------------------------

  SUBROUTINE megan_free_memory

    IMPLICIT NONE
    INTRINSIC ASSOCIATED

    INTEGER :: jt

    IF (ASSOCIATED(gam_smt))   DEALLOCATE(gam_smt)
    IF (ASSOCIATED(gam_tmp))   DEALLOCATE(gam_tmp)
    IF (ASSOCIATED(gam_lht))   DEALLOCATE(gam_lht)
    IF (ASSOCIATED(gam_pho))   DEALLOCATE(gam_pho)
    IF (ASSOCIATED(gam_age))   DEALLOCATE(gam_age)
    IF (ASSOCIATED(gam_total)) DEALLOCATE(gam_total)

    IF (ASSOCIATED(flux_gp))   DEALLOCATE(flux_gp)

    IF (ASSOCIATED(PPFD))      DEALLOCATE(PPFD)
    IF (ASSOCIATED(TEMP))      DEALLOCATE(TEMP)
    IF (ASSOCIATED(EF))        DEALLOCATE(EF)

  END SUBROUTINE megan_free_memory

! ===========================================================================
! PRIVATE ROUTINES
! ===========================================================================

  SUBROUTINE megan_read_nml_cpl(status, iou)
   
    ! read namelist for 'coupling' to ECHAM5
    !
    ! Author: Andrea Pozzer, MPICH, Feb 2011

    ! MESSy
    USE messy_main_tools,         ONLY: read_nml_open, read_nml_check, read_nml_close
    !USE messy_main_tracer_mem_bi, ONLY: NCELL
    !USE messy_main_mpi_bi,        ONLY: finish

    IMPLICIT NONE

    ! I/O
    INTEGER, INTENT(OUT) :: status     ! error status
    INTEGER, INTENT(IN)  :: iou        ! I/O unit

    NAMELIST /CPL/ sw_flux, lai, laip,              &
            btr_frac, ntr_frac, shr_frac, hrb_frac, cossza
         

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr='megan_read_nml_cpl'
    LOGICAL                     :: lex      ! file exists ?
    INTEGER                     :: fstat    ! file status
    INTEGER                     :: jt

    status = 1

    CALL read_nml_open(lex, substr, iou, 'CPL', modstr)
    IF (.NOT.lex) RETURN    ! <modstr>.nml does not exist

    READ(iou, NML=CPL, IOSTAT=fstat)
    CALL read_nml_check(fstat, substr, iou, 'CPL', modstr)
    IF (fstat /= 0) RETURN  ! error while reading namelist

    CALL read_nml_close(substr, iou, modstr)
    status = 0 ! NO ERROR

  END SUBROUTINE megan_read_nml_cpl

! ===========================================================================

  SUBROUTINE megan_read_nml_mgn(status, iou)
   
    ! read namelist for 'coupling' to ECHAM5
    !
    ! Author: Andrea Pozzer, MPICH, Feb 2011

    ! MESSy
    USE messy_main_tools,         ONLY: read_nml_open, read_nml_check, read_nml_close
    !USE messy_main_tracer_mem_bi, ONLY: NCELL
    !USE messy_main_mpi_bi,        ONLY: finish

    IMPLICIT NONE

    ! I/O
    INTEGER, INTENT(OUT) :: status     ! error status
    INTEGER, INTENT(IN)  :: iou        ! I/O unit

    NAMELIST /MGN/  MGN_SPC, LDF_FCT,RHO_FCT,TDF_PRM,       &
            EF_BT, EF_NT, EF_SB, EF_HB                    
         

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr='megan_read_nml_mgn'
    LOGICAL                     :: lex      ! file exists ?
    INTEGER                     :: fstat    ! file status
    INTEGER                     :: jt

    status = 1

    CALL read_nml_open(lex, substr, iou, 'MGN', modstr)
    IF (.NOT.lex) RETURN    ! <modstr>.nml does not exist

    READ(iou, NML=MGN, IOSTAT=fstat)
    CALL read_nml_check(fstat, substr, iou, 'MGN', modstr)
    IF (fstat /= 0) RETURN  ! error while reading namelist

    CALL read_nml_close(substr, iou, modstr)
    status = 0 ! NO ERROR

   WRITE(*,*) '.........................................................'
   WRITE(*,*) '             MEGAN  TRACER DEFINED                       '
   WRITE(*,*) '.........................................................'

    DO jt=1, NMAXNTRAC_MEGAN
       IF (TRIM(MGN_SPC(jt)) == '') CYCLE

       WRITE(*,*) '  TRACER NO.          ',jt
       WRITE(*,*) '  NAME              = ', TRIM(MGN_SPC(jt))
       WRITE(*,*) '  LDF_FCT           = ', LDF_FCT(jt)
       WRITE(*,*) '  RHO_FCT           = ', RHO_FCT(jt)
       WRITE(*,*) '  TDF_PRM           = ', TDF_PRM(jt)
       WRITE(*,*) '  EF_BT             = ', EF_BT(jt)
       WRITE(*,*) '  EF_NT             = ', EF_NT(jt)
       WRITE(*,*) '  EF_SB             = ', EF_SB(jt)
       WRITE(*,*) '  EF_HB             = ', EF_HB(jt)
       WRITE(*,*) '.........................................................'

    END DO 

  END SUBROUTINE megan_read_nml_mgn
! ===========================================================================

  SUBROUTINE megan_read_nml_mgn2mech(status, iou)
   
    ! read namelist for 'coupling' to ECHAM5
    !
    ! Author: Andrea Pozzer, MPICH, Feb 2011

    ! MESSy
    USE messy_main_tools,         ONLY: read_nml_open, read_nml_check, read_nml_close
    !USE messy_main_tracer_mem_bi, ONLY: NCELL
    !USE messy_main_mpi_bi,        ONLY: finish

    IMPLICIT NONE

    ! I/O
    INTEGER, INTENT(OUT) :: status     ! error status
    INTEGER, INTENT(IN)  :: iou        ! I/O unit

! mz_dc_20131118+
    NAMELIST /MGN2MECH/  MECH_SPC, MGN_SRC,MECH_MWT, &
                    EFFS_BT,EFFS_NT,EFFS_SB,EFFS_HB, GS           

! mz_dc_20131118-
    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr='megan_read_nml_mgn2mech'
    LOGICAL                     :: lex      ! file exists ?
    INTEGER                     :: fstat    ! file status
    INTEGER                     :: jt

    status = 1

    CALL read_nml_open(lex, substr, iou, 'MGN2MECH', modstr)
    IF (.NOT.lex) RETURN    ! <modstr>.nml does not exist

    READ(iou, NML=MGN2MECH, IOSTAT=fstat)
    CALL read_nml_check(fstat, substr, iou, 'MGN2MECH', modstr)
    IF (fstat /= 0) RETURN  ! error while reading namelist

    CALL read_nml_close(substr, iou, modstr)
    status = 0 ! NO ERROR

   WRITE(*,*) '.........................................................'
   WRITE(*,*) '             MECHANISM  TRACER DEFINED                   '
   WRITE(*,*) '.........................................................'

    DO jt=1, NMAXNTRAC_MECH
       IF (TRIM(MECH_SPC(jt)) == '') CYCLE

       WRITE(*,*) '  TRACER NO.          ',jt
       WRITE(*,*) '  NAME              = ', TRIM(MECH_SPC(jt))
       WRITE(*,*) '  SOURCE (MEGAN)    = ', TRIM(MGN_SRC(jt))
       WRITE(*,*) '  MOLAR WEIGHT      = ', MECH_MWT(jt)
       WRITE(*,*) '  EFFS_BT           = ', EFFS_BT(jt)
       WRITE(*,*) '  EFFS_NT           = ', EFFS_NT(jt)
       WRITE(*,*) '  EFFS_SB           = ', EFFS_SB(jt)
       WRITE(*,*) '  EFFS_HB           = ', EFFS_HB(jt)
! mz_dc_20131118+
       WRITE(*,*) '  GS                = ', GS(jt)
! mz_dc_20131118-
       WRITE(*,*) '.........................................................'

    END DO 

  END SUBROUTINE megan_read_nml_mgn2mech
! ===========================================================================
 
END MODULE messy_megan_si

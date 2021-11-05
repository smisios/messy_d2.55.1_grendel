#include "messy_main_ppd_bi.inc"

MODULE messy_m7_si

  ! MODULE FOR M7-ECHAM5 INTERFACE
  !
  ! Astrid Kerkweg (akerkweg@mpch-mainz.mpg.de), MPI-CHEM, Dec 2003
  ! Swen Metzger    (metzger@mpch-mainz.mpg.de), MPI-CHEM, Dec 2003
  ! M7 adopted to the structure of the Mainz Earth Submodel System (MESSy)
  ! M7 was originally implemented by Philip Stier, MPI-Met, Hamburg, 2001-2003
  ! Original M7 source code (box model) by 
  ! J. Wilson, E. Vignati, JRC/EI, 09/2000

  ! BMIL/MESSy
  USE messy_main_tracer_mem_bi, ONLY: ti_gp
  USE messy_main_blather_bi,    ONLY: start_message_bi, end_message_bi &
                                    , error_bi, info_bi
#ifdef MESSYTENDENCY
 !tendency budget
 USE messy_main_tendency_bi,    ONLY: mtend_get_handle,       &
                                      mtend_get_start_l,      &
                                      mtend_add_l,            &
                                      mtend_register,         &    
                                      mtend_id_tracer,        &    
                                      mtend_id_q,             &    
                                      mtend_id_t
#endif
  ! MESSy
  USE messy_main_constants_mem, ONLY: STRLEN_MEDIUM

  ! SUBMODEL M7
  USE messy_m7

  IMPLICIT NONE
  SAVE
  PRIVATE

  INTRINSIC  TRIM, ABS, ANY, MAX, MIN, NULL

  ! define radius and density 5D pointer
  REAL(dp), DIMENSION(:,:,:,:,:), POINTER :: wetradius => NULL()
  REAL(dp), DIMENSION(:,:,:,:,:), POINTER :: dryradius => NULL()
  REAL(dp), DIMENSION(:,:,:,:,:), POINTER :: densaer => NULL()

  ! define radius and density 4D pointer
  REAL(dp), DIMENSION(:,:,:,:),   POINTER :: wetrad_4d    => NULL()
  REAL(dp), DIMENSION(:,:,:,:),   POINTER :: dryrad_4d    => NULL()
  REAL(dp), DIMENSION(:,:,:,:),   POINTER :: densaer_4d   => NULL()

  ! seasalt emission calculations 
  REAL(dp), DIMENSION(:,:), POINTER :: Mss_as => NULL()
  REAL(dp), DIMENSION(:,:), POINTER :: Mss_cs => NULL()
  REAL(dp), DIMENSION(:,:), POINTER :: Nss_as => NULL()
  REAL(dp), DIMENSION(:,:), POINTER :: Nss_cs => NULL()
  
  !carbon emission mass fluxes insoluble part in kg/kg
  REAL(dp) ,DIMENSION(:,:), POINTER :: OC_sum_insol => NULL()
  !carbon emission mass fluxes soluble part in kg/kg
  REAL(dp) ,DIMENSION(:,:), POINTER :: OC_sum_sol => NULL()
  REAL(dp) ,DIMENSION(:,:), POINTER :: BC_sum_insol => NULL()
  !number emission fluxes resulting from carbon emissions in 1 /kg
  REAL(dp) ,DIMENSION(:,:), POINTER :: Nemis_ks => NULL()
  REAL(dp) ,DIMENSION(:,:), POINTER :: Nemis_ki => NULL()
  ! wildfire emission channel 3d
  REAL(dp) ,DIMENSION(:,:,:), POINTER :: BC_wf => NULL()
  REAL(dp) ,DIMENSION(:,:,:), POINTER :: OC_wf => NULL()

  ! dust emission mass flux
  REAL(dp), DIMENSION(:,:), POINTER :: dust_emis_ai => NULL()
  REAL(dp), DIMENSION(:,:), POINTER :: dust_emis_ci => NULL()

  ! sulfate emission fluxes in kg(SO2) m-2 s-1
  REAL(dp), DIMENSION(:,:), POINTER :: SO2_emis_low  => NULL()
  REAL(dp), DIMENSION(:,:), POINTER :: SO2_emis_high => NULL()

  ! aerosol column mass
  REAL(dp), DIMENSION(:,:), POINTER :: colmass_SO4   => NULL()
  REAL(dp), DIMENSION(:,:), POINTER :: colmass_BC    => NULL()
  REAL(dp), DIMENSION(:,:), POINTER :: colmass_OC    => NULL()
  REAL(dp), DIMENSION(:,:), POINTER :: colmass_SS    => NULL()
  REAL(dp), DIMENSION(:,:), POINTER :: colmass_DU    => NULL()
  REAL(dp), DIMENSION(:,:), POINTER :: colmass_DU_as => NULL()
  REAL(dp), DIMENSION(:,:), POINTER :: colmass_DU_cs => NULL()
  REAL(dp), DIMENSION(:,:), POINTER :: colmass_DU_ai => NULL()
  REAL(dp), DIMENSION(:,:), POINTER :: colmass_DU_ci => NULL()

  ! aerosol optical depth
  REAL(dp), DIMENSION(:,:), POINTER :: AOD_SO4_COLUMN   => NULL()
  REAL(dp), DIMENSION(:,:), POINTER :: AOD_BC_COLUMN    => NULL()
  REAL(dp), DIMENSION(:,:), POINTER :: AOD_OC_COLUMN    => NULL()
  REAL(dp), DIMENSION(:,:), POINTER :: AOD_SS_COLUMN    => NULL()
  REAL(dp), DIMENSION(:,:), POINTER :: AOD_DU_COLUMN    => NULL()
  REAL(dp), DIMENSION(:,:), POINTER :: AOD_TOTAL_COLUMN => NULL()

  ! CPL - NAMELIST 
  LOGICAL :: lcpl_gasphase      = .false.
  CHARACTER (LEN=STRLEN_MEDIUM) :: chemmodule    = ''
  CHARACTER (LEN=STRLEN_MEDIUM) :: H2SO4_gas(2)  = ''
  CHARACTER (LEN=STRLEN_MEDIUM) :: H2SO4_as(2)   = ''
  CHARACTER (LEN=STRLEN_MEDIUM) :: H2SO4_cs(2)   = ''
  CHARACTER (LEN=STRLEN_MEDIUM) :: HSO4m_cs(2)   = ''
  CHARACTER (LEN=STRLEN_MEDIUM) :: SO4mm_cs(2)   = ''

  LOGICAL :: l_notrac_SO4       = .false.
  LOGICAL :: l_notrac_OCBC      = .false.
  LOGICAL :: l_notrac_SS        = .false.
  LOGICAL :: l_notrac_DU        = .false.

  LOGICAL :: l_calc_emis        = .false.
  LOGICAL :: l_tendency         = .false.
  ! sea salt emissions
  LOGICAL :: l_ss               = .false.
  CHARACTER (LEN=STRLEN_MEDIUM) :: SSemis_channel = ''
  CHARACTER (LEN=STRLEN_MEDIUM) :: SS_mass_as    = ''
  CHARACTER (LEN=STRLEN_MEDIUM) :: SS_num_as     = ''
  CHARACTER (LEN=STRLEN_MEDIUM) :: SS_mass_cs    = ''
  CHARACTER (LEN=STRLEN_MEDIUM) :: SS_num_cs     = ''
  ! carbon emissions
  LOGICAL :: l_carbon           = .false.
  CHARACTER (LEN=STRLEN_MEDIUM) :: Cemis_channel  = ''
  ! - organic carbon mass
  CHARACTER (LEN=STRLEN_MEDIUM) :: emis_OC_sol   = ''
  CHARACTER (LEN=STRLEN_MEDIUM) :: emis_OC_insol = ''
  ! - black carbon mass
  CHARACTER (LEN=STRLEN_MEDIUM) :: emis_BC_insol = ''
  ! - number carbon
  CHARACTER (LEN=STRLEN_MEDIUM) :: emis_N_sol    = ''
  CHARACTER (LEN=STRLEN_MEDIUM) :: emis_N_insol  = ''

  ! dust emissions
  ! Logical l_dust, l_dust_tegen replaced by l_dust_ai, l_dust_ci
  !          Schulze Scheme: only l_dust_ci
  !          Tegen Scheme:   l_dust_ci, l_dust_ai
  !          Astitha Scheme: l_dust_ci, l_dust_ai
  LOGICAL :: l_dust_ci   = .false.
  LOGICAL :: l_dust_ai   = .false.
  CHARACTER (LEN=STRLEN_MEDIUM) :: Duemis_channel = ''
  CHARACTER (LEN=STRLEN_MEDIUM) :: emis_dust_ci  = ''
  CHARACTER (LEN=STRLEN_MEDIUM) :: emis_dust_ai  = ''
  ! sulfate emissions
  LOGICAL :: l_so4               = .false.
  CHARACTER (LEN=STRLEN_MEDIUM) :: SO4emis_channel = ''
  CHARACTER (LEN=STRLEN_MEDIUM) :: emis_so2_low   = ''
  CHARACTER (LEN=STRLEN_MEDIUM) :: emis_so2_high  = ''
  ! calculate aerosol column masses
  LOGICAL :: l_colmass           = .false.

  ! calculate aerosol optical depth
  LOGICAL :: l_aod               = .false.

#ifdef MESSYTENDENCY
  INTEGER :: my_handle
#endif

    NAMELIST /CPL/ lcpl_gasphase, chemmodule, H2SO4_gas,        &
         H2SO4_as, H2SO4_cs,   HSO4m_cs,   SO4mm_cs,            &
         l_calc_emis,  l_tendency,                              &
         l_notrac_SO4, l_notrac_OCBC, l_notrac_SS, l_notrac_DU, &
         l_ss,         SSemis_channel,                          &
         SS_mass_as,   SS_num_as,     SS_mass_cs, SS_num_cs,    &
         Cemis_channel, l_carbon,                               &
         emis_OC_sol,  emis_OC_insol, emis_BC_insol,            &
         emis_N_sol,   emis_N_insol,                            &
         Duemis_channel, emis_dust_ci, emis_dust_ai,            &
         l_so4,        SO4emis_channel,                         &
         emis_so2_low, emis_so2_high,                           &
         l_colmass,  l_aod

  ! TRACER INDICES
  INTEGER :: idt_so4_m7 = 0     ! mass mixing ratio so4
  INTEGER :: idt_h2so4_as = 0   ! mass mixing ratio h2so4
  INTEGER :: idt_h2so4_cs = 0   ! mass mixing ratio h2so4
  INTEGER :: idt_hso4m_as = 0   ! mass mixing ratio hso4m
  INTEGER :: idt_hso4m_cs = 0   ! mass mixing ratio hso4m
  INTEGER :: idt_so4mm_as = 0   ! mass mixing ratio so4mm
  INTEGER :: idt_so4mm_cs = 0   ! mass mixing ratio so4mm
  !
  INTEGER :: idt_ms4ns  = 0 ! mass mixing ratio sulfate        nuclea. soluble
  INTEGER :: idt_ms4ks  = 0 ! mass mixing ratio sulfate        aitken  soluble
  INTEGER :: idt_ms4as  = 0 ! mass mixing ratio sulfate        accum.  soluble
  INTEGER :: idt_ms4cs  = 0 ! mass mixing ratio sulfate        coarse  soluble
  INTEGER :: idt_mbcki  = 0 ! mass mixing ratio black carbon   aitken  insol.
  INTEGER :: idt_mbcks  = 0 ! mass mixing ratio black carbon   aitken  soluble
  INTEGER :: idt_mbcas  = 0 ! mass mixing ratio black carbon   accum.  soluble
  INTEGER :: idt_mbccs  = 0 ! mass mixing ratio black carbon   coarse  soluble
  INTEGER :: idt_mocki  = 0 ! mass mixing ratio organic carbon aitken  insol.
  INTEGER :: idt_mocks  = 0 ! mass mixing ratio organic carbon aitken  soluble
  INTEGER :: idt_mocas  = 0 ! mass mixing ratio organic carbon accum.  soluble
  INTEGER :: idt_moccs  = 0 ! mass mixing ratio organic carbon coarse  soluble
  INTEGER :: idt_mssas  = 0 ! mass mixing ratio seasalt        accum.  soluble
  INTEGER :: idt_msscs  = 0 ! mass mixing ratio seasalt        coarse  soluble
  INTEGER :: idt_mduai  = 0 ! mass mixing ratio dust           accum.  insol.
  INTEGER :: idt_mduas  = 0 ! mass mixing ratio dust           accum.  soluble
  INTEGER :: idt_mduci  = 0 ! mass mixing ratio dust           coarse  insol.
  INTEGER :: idt_mducs  = 0 ! mass mixing ratio dust           coarse  soluble
  !
  INTEGER :: idt_nns    = 0 ! number mixing ratio              nuclea. soluble
  INTEGER :: idt_nki    = 0 ! number mixing ratio              aitken  insol.
  INTEGER :: idt_nks    = 0 ! number mixing ratio              aitken  soluble
  INTEGER :: idt_nai    = 0 ! number mixing ratio              accum.  insol.
  INTEGER :: idt_nas    = 0 ! number mixing ratio              accum.  soluble
  INTEGER :: idt_nci    = 0 ! number mixing ratio              coarse  insol.
  INTEGER :: idt_ncs    = 0 ! number mixing ratio              coarse  soluble
  !
  REAL(dp), DIMENSION(:), POINTER :: sigma_cha

  ! SUBROUTINES
  PUBLIC :: m7_initialize      ! initialize m7
  !PRIVATE :: m7_read_nml_cpl
  PUBLIC :: m7_new_tracer      ! define tracers
  PUBLIC :: m7_init_memory     ! allocate memory
  PUBLIC :: m7_init_coupling   ! initialize coupling
  PUBLIC :: m7_vdiff           ! distribute online emissions
  PUBLIC :: m7_physc           ! calculate m7-"chemistry"
  PUBLIC :: m7_local_end       ! diagnose column mass and aod ! um_gg_20091030
  PUBLIC :: m7_free_memory     ! deallocate radius field

CONTAINS

!-----------------------------------------------------------------------------
  SUBROUTINE m7_initialize

    ! ECHAM5/MESSy
    USE messy_main_mpi_bi,             ONLY: p_parallel_io, p_io, p_bcast
    USE messy_main_tools,              ONLY: find_next_free_unit
    
    IMPLICIT NONE
  
    ! LOCAL
    CHARACTER(LEN=*), PARAMETER  :: substr = 'm7_initialize'
    INTEGER                      :: iou    ! I/O unit
    INTEGER(i4)                  :: status ! error status

    !--- 1) Read M7 namelist and control variables:

    ! INITIALIZE MAIN-CTRL
    IF (p_parallel_io) THEN
       iou = find_next_free_unit(100,200)
       ! *** CALL CORE ROUTINE:
       CALL m7_read_nml_ctrl(status, iou)
       IF (status /= 0)  CALL error_bi('error in m7_read_nml_ctrl', substr)
    END IF
    ! BROADCAST RESULTS
    CALL p_bcast (lm7,        p_io)
    CALL p_bcast (lmass_diag, p_io)
    CALL p_bcast (lcdnc,      p_io)
    CALL p_bcast (licnc,      p_io)
    CALL p_bcast (lscoag,     p_io)
    CALL p_bcast (lscond,     p_io)
    CALL p_bcast (lsnucl,     p_io)
    CALL p_bcast (nnucl,      p_io)

    !--- 2)  Read M7 CPL namelist
    IF (p_parallel_io) THEN
       iou = find_next_free_unit(100,200)
       CALL m7_read_nml_cpl(status, iou)
       IF (status /= 0) CALL error_bi(' ',substr)
    END IF
    ! BROADCAST RESULTS
    CALL p_bcast(lcpl_gasphase,     p_io)
    !
    CALL p_bcast(chemmodule,    p_io)
    CALL p_bcast(H2SO4_gas(1),  p_io)
    CALL p_bcast(H2SO4_gas(2),  p_io)
    !
    CALL p_bcast(H2SO4_as(1),  p_io)
    CALL p_bcast(H2SO4_as(2),  p_io)
    CALL p_bcast(H2SO4_cs(1),  p_io)
    CALL p_bcast(H2SO4_cs(2),  p_io)
    CALL p_bcast(HSO4m_cs(1),  p_io)
    CALL p_bcast(HSO4m_cs(2),  p_io)
    CALL p_bcast(SO4mm_cs(1),  p_io)
    CALL p_bcast(SO4mm_cs(2),  p_io)
    !
    CALL p_bcast(l_calc_emis,   p_io)

    CALL p_bcast(l_notrac_SO4,  p_io)
    CALL p_bcast(l_notrac_OCBC, p_io)
    CALL p_bcast(l_notrac_SS,   p_io)
    CALL p_bcast(l_notrac_DU,   p_io)

    CALL p_bcast(l_tendency,    p_io)
    !
    CALL p_bcast(l_ss,          p_io)
    CALL p_bcast(SSemis_channel, p_io)
    !
    CALL p_bcast(SS_mass_as,    p_io)
    CALL p_bcast(SS_num_as,     p_io)
     !
    CALL p_bcast(SS_mass_cs,    p_io)
    CALL p_bcast(SS_num_cs,     p_io)
    !
    CALL p_bcast(l_carbon,      p_io)
    CALL p_bcast(Cemis_channel,  p_io)
    !
    CALL p_bcast(emis_OC_sol,   p_io)
    CALL p_bcast(emis_OC_insol, p_io)
    CALL p_bcast(emis_BC_insol, p_io)
    CALL p_bcast(emis_N_sol,    p_io)
    CALL p_bcast(emis_N_insol,  p_io)
    !
    CALL p_bcast(Duemis_channel, p_io)
    CALL p_bcast(emis_dust_ci,   p_io)
    CALL p_bcast(emis_dust_ai,   p_io)
    ! 
    l_dust_ai = (TRIM(emis_dust_ai) /= "" )
    l_dust_ci = (TRIM(emis_dust_ci) /= "" )
    !
    !
    CALL p_bcast(l_so4,           p_io)
    CALL p_bcast(SO4emis_channel, p_io)
    !
    CALL p_bcast(emis_so2_low,   p_io)
    CALL p_bcast(emis_so2_high,  p_io)
    !
    CALL p_bcast(l_colmass,      p_io)
    CALL p_bcast(l_aod,          p_io)

    !--- 4) Initialize M7:
    CALL m7_initialize_core

#ifdef MESSYTENDENCY
    my_handle = mtend_get_handle(modstr)
#endif

  END SUBROUTINE m7_initialize
!-----------------------------------------------------------------------------

!-----------------------------------------------------------------------------
  SUBROUTINE m7_new_tracer

    ! ECHAM5/MESSy
    USE messy_main_blather_bi,      ONLY: warning_bi
    USE messy_main_tracer_mem_bi,   ONLY: GPTRSTR
    USE messy_main_tracer_tools_bi, ONLY: tracer_halt
    ! MESSy
    USE messy_main_tracer,        ONLY: new_tracer, set_tracer, get_tracer &
                                      , AIR, AEROSOL  &
                                      , AMOUNTFRACTION, NUMBERDENSITY      &
                                      , OFF, ON, MODAL                     &
                                      , R_MOLARMASS, R_AEROSOL_DENSITY     &
                                      , R_dryreac_sf                       &
                                      , I_AEROSOL_SOL, S_AEROSOL_MODEL     &
                                      , I_AEROSOL_METHOD, I_AEROSOL_MODE   &
                                      , I_DRYDEP, I_SEDI, I_SCAV

    !--- Module variables:
    !
    !    Tracer indices:
    !
    !    Legend: idt_ABBCD
    !
    !      A:  m  = particle mass mixing ratio, n number mixing ratio
    !      BB: s4 = sulfate, bc/oc = black/organic carbon, 
    !               du = dust, ss = seasalt
    !      C:  n  = nucleation , k = Aitken, a = accumulation, c = coarse mode
    !      D:  i  = insoluble,  s = soluble
    !
    ! Parameters:
    ! -----------
    ! User defined flags: density   density                    [kg m-3]
    !                     osm       osmotic coefficient        [???]
    !                     nion      number of ions the tracer 
    !                               dissolves into             [1]

    ! LOCAL
    INTEGER :: status
    CHARACTER(LEN=*), PARAMETER :: substr = 'm7_new_tracer'
    INTEGER :: idt_dummy = 0

    CALL start_message_bi(modstr,'REQUEST M7 TRACER', substr)

    ! create dummy tracer to avoid if causes in time loop
    CALL new_tracer(status, GPTRSTR, 'dummy', modstr      &
         , idx=idt_dummy , unit='mol/mol',medium=AEROSOL  &
         , quantity=AMOUNTFRACTION )
    CALL tracer_halt(substr, status)
    CALL set_tracer(status, GPTRSTR, idt_dummy, R_molarmass,       r=1._dp)
    CALL tracer_halt(substr, status)
    CALL set_tracer(status, GPTRSTR, idt_dummy, R_AEROSOL_DENSITY, r=1._dp)
    CALL tracer_halt(substr, status)
    CALL set_tracer(status, GPTRSTR, idt_dummy, I_AEROSOL_MODE,     1)
    CALL tracer_halt(substr, status)
    CALL set_tracer(status, GPTRSTR, idt_dummy, I_AEROSOL_METHOD,   MODAL)
    CALL tracer_halt(substr, status)
    CALL set_tracer(status, GPTRSTR, idt_dummy, S_AEROSOL_MODEL,   modstr)
    CALL tracer_halt(substr, status)
#ifdef MESSYTENDENCY
    CALL mtend_register(my_handle, idt_dummy)
#endif

    if (lcpl_gasphase) then
       if (l_notrac_SO4) then
          l_notrac_SO4 = .false.
          CALL warning_bi('because lcpl_gasphase = .true.',substr)
          CALL warning_bi('setting l_notrac_SO4 = .false.',substr)
       endif
    else
       idt_so4_m7 = idt_dummy
    endif

    CALL get_tracer(status, GPTRSTR, H2SO4_gas(1) &
         ,subname= H2SO4_gas(2), idx=idt_so4_m7)
    if (status.ne.0) then 
       CALL new_tracer(status, GPTRSTR, 'H2SO4', modstr      &
            , idx=idt_so4_m7 , unit='mol/mol', medium=AIR    &
            , quantity=AMOUNTFRACTION )
       CALL tracer_halt(substr, status)
       CALL set_tracer(status, GPTRSTR, idt_so4_m7, R_molarmass,  r=98.076_dp)
       CALL tracer_halt(substr, status)
       CALL set_tracer(status, GPTRSTR, idt_so4_m7, I_DRYDEP, ON)
       CALL tracer_halt(substr, status)
       CALL set_tracer(status, GPTRSTR, idt_so4_m7, I_SEDI,   ON)
       CALL tracer_halt(substr, status)
       CALL set_tracer(status, GPTRSTR, idt_so4_m7, I_SCAV,   ON)
       CALL tracer_halt(substr, status)
       CALL set_tracer(status, GPTRSTR, idt_so4_m7, R_dryreac_sf, r=0._dp)
       CALL tracer_halt(substr, status)
    endif
#ifdef MESSYTENDENCY
    CALL mtend_register(my_handle, idt_so4_m7)
#endif

    !--- 2) Allocate aerosol masses according to the modes to obtain 
    !       succeding tracer identifiers:

    !--- Mode 1 - Nucleation Soluble:

    if (l_notrac_SO4) then
       idt_ms4ns = idt_dummy
    else
       CALL new_tracer(status, GPTRSTR, 'SO4', modstr , subname='ns'     &
            , idx=idt_ms4ns , unit='mol/mol', medium=AEROSOL  &
            , quantity=AMOUNTFRACTION )
       CALL tracer_halt(substr, status)
       CALL set_tracer(status, GPTRSTR, idt_ms4ns, R_molarmass,      r=96.076_dp)
       CALL tracer_halt(substr, status)
       CALL set_tracer(status, GPTRSTR, idt_ms4ns, R_AEROSOL_DENSITY, r=1841._dp)
       CALL tracer_halt(substr, status)
       CALL set_tracer(status, GPTRSTR, idt_ms4ns, I_AEROSOL_MODE,     1)
       CALL tracer_halt(substr, status)
       CALL set_tracer(status, GPTRSTR, idt_ms4ns, I_AEROSOL_METHOD,   MODAL)
       CALL tracer_halt(substr, status)
       CALL set_tracer(status, GPTRSTR, idt_ms4ns, S_AEROSOL_MODEL,   modstr)
       CALL tracer_halt(substr, status)
       CALL set_tracer(status, GPTRSTR, idt_ms4ns, I_DRYDEP, ON)
       CALL tracer_halt(substr, status)
       CALL set_tracer(status, GPTRSTR, idt_ms4ns, I_SEDI,   ON)
       CALL tracer_halt(substr, status)
       CALL set_tracer(status, GPTRSTR, idt_ms4ns, I_SCAV,   ON)
       CALL tracer_halt(substr, status)
#ifdef MESSYTENDENCY
    CALL mtend_register(my_handle, idt_ms4ns)
#endif
    endif
    !--- Mode 2 - Aitken Soluble:

    if (l_notrac_SO4) then
       idt_ms4ks = idt_dummy
    else
       CALL new_tracer(status, GPTRSTR, 'SO4', modstr , subname='ks'     &
            , idx=idt_ms4ks , unit='mol/mol', medium=AEROSOL  &
            , quantity=AMOUNTFRACTION )
       CALL tracer_halt(substr, status)
       CALL set_tracer(status, GPTRSTR, idt_ms4ks, R_molarmass,    r=96.076_dp)
       CALL tracer_halt(substr, status)
       CALL set_tracer(status, GPTRSTR, idt_ms4ks, R_AEROSOL_DENSITY, r=1841._dp)
       CALL tracer_halt(substr, status)
       CALL set_tracer(status, GPTRSTR, idt_ms4ks, I_AEROSOL_MODE,     2)
       CALL tracer_halt(substr, status)
       CALL set_tracer(status, GPTRSTR, idt_ms4ks, I_AEROSOL_METHOD,   MODAL)
       CALL tracer_halt(substr, status)
       CALL set_tracer(status, GPTRSTR, idt_ms4ks, S_AEROSOL_MODEL,   modstr)
       CALL tracer_halt(substr, status)
       CALL set_tracer(status, GPTRSTR, idt_ms4ks, I_DRYDEP, ON)
       CALL tracer_halt(substr, status)
       CALL set_tracer(status, GPTRSTR, idt_ms4ks, I_SEDI,   ON)
       CALL tracer_halt(substr, status)
       CALL set_tracer(status, GPTRSTR, idt_ms4ks, I_SCAV,   ON)
       CALL tracer_halt(substr, status)
#ifdef MESSYTENDENCY
    CALL mtend_register(my_handle, idt_ms4ks)
#endif
    endif

    if (l_notrac_OCBC) then
       idt_mbcks = idt_dummy
    else
       CALL new_tracer(status, GPTRSTR, 'BC', modstr , subname='ks'     &
            , idx=idt_mbcks , unit='mol/mol', medium=AEROSOL  &
            , quantity=AMOUNTFRACTION )
       CALL tracer_halt(substr, status)
       CALL set_tracer(status, GPTRSTR, idt_mbcks, R_molarmass,      r=12.011_dp)
       CALL tracer_halt(substr, status)
       CALL set_tracer(status, GPTRSTR, idt_mbcks, R_AEROSOL_DENSITY, r=2000._dp)
       CALL tracer_halt(substr, status)
       CALL set_tracer(status, GPTRSTR, idt_mbcks, I_AEROSOL_MODE,     2)
       CALL tracer_halt(substr, status)
       CALL set_tracer(status, GPTRSTR, idt_mbcks, I_AEROSOL_METHOD,   MODAL)
       CALL tracer_halt(substr, status)
       CALL set_tracer(status, GPTRSTR, idt_mbcks, S_AEROSOL_MODEL,   modstr)
       CALL tracer_halt(substr, status)
       CALL set_tracer(status, GPTRSTR, idt_mbcks, I_DRYDEP, ON)
       CALL tracer_halt(substr, status)
       CALL set_tracer(status, GPTRSTR, idt_mbcks, I_SEDI,   ON)
       CALL tracer_halt(substr, status)
       CALL set_tracer(status, GPTRSTR, idt_mbcks, I_SCAV,   ON)
       CALL tracer_halt(substr, status)
#ifdef MESSYTENDENCY
    CALL mtend_register(my_handle, idt_mbcks)
#endif
    endif

    if (l_notrac_OCBC) then
       idt_mocks = idt_dummy
    else
       CALL new_tracer(status, GPTRSTR, 'OC', modstr , subname='ks'     &
            , idx=idt_mocks , unit='mol/mol', medium=AEROSOL  &
            , quantity=AMOUNTFRACTION )
       CALL tracer_halt(substr, status)
       CALL set_tracer(status, GPTRSTR, idt_mocks, R_molarmass,      r=12.011_dp)
       CALL tracer_halt(substr, status)
       CALL set_tracer(status, GPTRSTR, idt_mocks, R_AEROSOL_DENSITY, r=2000._dp)
       CALL tracer_halt(substr, status)
       CALL set_tracer(status, GPTRSTR, idt_mocks, I_AEROSOL_MODE,     2)
       CALL tracer_halt(substr, status)
       CALL set_tracer(status, GPTRSTR, idt_mocks, I_AEROSOL_METHOD,   MODAL)
       CALL tracer_halt(substr, status)
       CALL set_tracer(status, GPTRSTR, idt_mocks, S_AEROSOL_MODEL,   modstr)
       CALL tracer_halt(substr, status)
       CALL set_tracer(status, GPTRSTR, idt_mocks, I_DRYDEP, ON)
       CALL tracer_halt(substr, status)
       CALL set_tracer(status, GPTRSTR, idt_mocks, I_SEDI,   ON)
       CALL tracer_halt(substr, status)
       CALL set_tracer(status, GPTRSTR, idt_mocks, I_SCAV,   ON)
       CALL tracer_halt(substr, status)
#ifdef MESSYTENDENCY
    CALL mtend_register(my_handle, idt_mocks)
#endif
    endif

    !--- Mode 3 - Accumulation Soluble:
    if (l_notrac_SO4) then
       idt_ms4as = idt_dummy
    else
       CALL new_tracer(status, GPTRSTR, 'SO4', modstr , subname='as'     &
            , idx=idt_ms4as , unit='mol/mol', medium=AEROSOL  &
            , quantity=AMOUNTFRACTION )
       CALL tracer_halt(substr, status)
       CALL set_tracer(status, GPTRSTR, idt_ms4as, R_molarmass, r=96.076_dp)
       CALL tracer_halt(substr, status)
       CALL set_tracer(status, GPTRSTR, idt_ms4as, R_AEROSOL_DENSITY, r=1841._dp)
       CALL tracer_halt(substr, status)
       CALL set_tracer(status, GPTRSTR, idt_ms4as, I_AEROSOL_MODE,    3)
       CALL tracer_halt(substr, status)
       CALL set_tracer(status, GPTRSTR, idt_ms4as, I_AEROSOL_METHOD,  MODAL)
       CALL tracer_halt(substr, status)
       CALL set_tracer(status, GPTRSTR, idt_ms4as, S_AEROSOL_MODEL,   modstr)
       CALL tracer_halt(substr, status)
       CALL set_tracer(status, GPTRSTR, idt_ms4as, I_DRYDEP, ON)
       CALL tracer_halt(substr, status)
       CALL set_tracer(status, GPTRSTR, idt_ms4as, I_SEDI,   ON)
       CALL tracer_halt(substr, status)
       CALL set_tracer(status, GPTRSTR, idt_ms4as, I_SCAV,   ON)
       CALL tracer_halt(substr, status)
#ifdef MESSYTENDENCY
    CALL mtend_register(my_handle, idt_ms4as)
#endif
    endif
    
    if (l_notrac_OCBC) then
       idt_mbcas = idt_dummy
    else
       CALL new_tracer(status, GPTRSTR, 'BC', modstr , subname='as'     &
            , idx=idt_mbcas , unit='mol/mol', medium=AEROSOL  &
         , quantity=AMOUNTFRACTION )
       CALL tracer_halt(substr, status)
       CALL set_tracer(status, GPTRSTR, idt_mbcas, R_molarmass,      r=12.011_dp)
       CALL tracer_halt(substr, status)
       CALL set_tracer(status, GPTRSTR, idt_mbcas, R_AEROSOL_DENSITY, r=2000._dp)
       CALL tracer_halt(substr, status)
       CALL set_tracer(status, GPTRSTR, idt_mbcas, I_AEROSOL_MODE,    3)
       CALL tracer_halt(substr, status)
       CALL set_tracer(status, GPTRSTR, idt_mbcas, I_AEROSOL_METHOD,  MODAL)
       CALL tracer_halt(substr, status)
       CALL set_tracer(status, GPTRSTR, idt_mbcas, S_AEROSOL_MODEL,   modstr)
       CALL tracer_halt(substr, status)
       CALL set_tracer(status, GPTRSTR, idt_mbcas, I_DRYDEP, ON)
       CALL tracer_halt(substr, status)
       CALL set_tracer(status, GPTRSTR, idt_mbcas, I_SEDI,   ON)
       CALL tracer_halt(substr, status)
       CALL set_tracer(status, GPTRSTR, idt_mbcas, I_SCAV,   ON)
       CALL tracer_halt(substr, status)
#ifdef MESSYTENDENCY
    CALL mtend_register(my_handle, idt_mbcas)
#endif
    endif

    if (l_notrac_OCBC) then
       idt_mocas = idt_dummy
    else
       CALL new_tracer(status, GPTRSTR, 'OC', modstr , subname='as'     &
            , idx=idt_mocas , unit='mol/mol', medium=AEROSOL  &
            , quantity=AMOUNTFRACTION )
       CALL tracer_halt(substr, status)
       CALL set_tracer(status, GPTRSTR, idt_mocas, R_molarmass,      r=12.011_dp)
       CALL tracer_halt(substr, status)
       CALL set_tracer(status, GPTRSTR, idt_mocas, R_AEROSOL_DENSITY, r=2000._dp)
       CALL tracer_halt(substr, status)
       CALL set_tracer(status, GPTRSTR, idt_mocas, I_AEROSOL_MODE,    3)
       CALL tracer_halt(substr, status)
       CALL set_tracer(status, GPTRSTR, idt_mocas, I_AEROSOL_METHOD,  MODAL)
       CALL tracer_halt(substr, status)
       CALL set_tracer(status, GPTRSTR, idt_mocas, S_AEROSOL_MODEL,   modstr)
       CALL tracer_halt(substr, status)
       CALL set_tracer(status, GPTRSTR, idt_mocas, I_DRYDEP, ON)
       CALL tracer_halt(substr, status)
       CALL set_tracer(status, GPTRSTR, idt_mocas, I_SEDI,   ON)
       CALL tracer_halt(substr, status)
       CALL set_tracer(status, GPTRSTR, idt_mocas, I_SCAV,   ON)
       CALL tracer_halt(substr, status)
#ifdef MESSYTENDENCY
    CALL mtend_register(my_handle, idt_mocas)
#endif
    endif

    if (l_notrac_SS) then
       idt_mssas = idt_dummy
    else
       CALL new_tracer(status, GPTRSTR, 'SS', modstr , subname='as'     &
            , idx=idt_mssas , unit='mol/mol', medium=AEROSOL  &
            , quantity=AMOUNTFRACTION )
       CALL tracer_halt(substr, status)
       CALL set_tracer(status, GPTRSTR, idt_mssas, R_molarmass,       r=58.44_dp)
       CALL tracer_halt(substr, status)
       CALL set_tracer(status, GPTRSTR, idt_mssas, R_AEROSOL_DENSITY, r=2165._dp)
       CALL tracer_halt(substr, status)
       CALL set_tracer(status, GPTRSTR, idt_mssas, I_AEROSOL_MODE,    3)
       CALL tracer_halt(substr, status)
       CALL set_tracer(status, GPTRSTR, idt_mssas, I_AEROSOL_METHOD,  MODAL)
       CALL tracer_halt(substr, status)
       CALL set_tracer(status, GPTRSTR, idt_mssas, S_AEROSOL_MODEL,   modstr)
       CALL tracer_halt(substr, status)
       CALL set_tracer(status, GPTRSTR, idt_mssas, I_DRYDEP, ON)
       CALL tracer_halt(substr, status)
       CALL set_tracer(status, GPTRSTR, idt_mssas, I_SEDI,   ON)
       CALL tracer_halt(substr, status)
       CALL set_tracer(status, GPTRSTR, idt_mssas, I_SCAV,   ON)
       CALL tracer_halt(substr, status)
#ifdef MESSYTENDENCY
    CALL mtend_register(my_handle, idt_mssas)
#endif
    endif

    if (l_notrac_DU .or. l_notrac_SO4) then
       idt_mduas = idt_dummy
    else
       CALL new_tracer(status, GPTRSTR, 'DU', modstr , subname='as'     &
            , idx=idt_mduas , unit='mol/mol', medium=AEROSOL  &
            , quantity=AMOUNTFRACTION )
       CALL tracer_halt(substr, status)
       CALL set_tracer(status, GPTRSTR, idt_mduas, R_molarmass,       r=40.08_dp)
       CALL tracer_halt(substr, status)
       CALL set_tracer(status, GPTRSTR, idt_mduas, R_AEROSOL_DENSITY, r=2650._dp)
       CALL tracer_halt(substr, status)
       CALL set_tracer(status, GPTRSTR, idt_mduas, I_AEROSOL_MODE,    3)
       CALL tracer_halt(substr, status)
       CALL set_tracer(status, GPTRSTR, idt_mduas, I_AEROSOL_METHOD,  MODAL)
       CALL tracer_halt(substr, status)
       CALL set_tracer(status, GPTRSTR, idt_mduas, S_AEROSOL_MODEL,   modstr)
       CALL tracer_halt(substr, status)
       CALL set_tracer(status, GPTRSTR, idt_mduas, I_DRYDEP, ON)
       CALL tracer_halt(substr, status)
       CALL set_tracer(status, GPTRSTR, idt_mduas, I_SEDI,   ON)
       CALL tracer_halt(substr, status)
       CALL set_tracer(status, GPTRSTR, idt_mduas, I_SCAV,   ON)
       CALL tracer_halt(substr, status)
#ifdef MESSYTENDENCY
    CALL mtend_register(my_handle, idt_mduas)
#endif
    endif

    !--- Mode 4 - Coarse Soluble:
  
    if (l_notrac_SO4) then
       idt_ms4cs = idt_dummy
    else
       CALL new_tracer(status, GPTRSTR, 'SO4', modstr , subname='cs'     &
            , idx=idt_ms4cs , unit='mol/mol', medium=AEROSOL  &
            , quantity=AMOUNTFRACTION )
       CALL tracer_halt(substr, status)
       CALL set_tracer(status, GPTRSTR, idt_ms4cs, R_molarmass,      r=96.076_dp)
       CALL tracer_halt(substr, status)
       CALL set_tracer(status, GPTRSTR, idt_ms4cs, R_AEROSOL_DENSITY, r=1841._dp)
       CALL tracer_halt(substr, status)
       CALL set_tracer(status, GPTRSTR, idt_ms4cs, I_AEROSOL_MODE,    4)
       CALL tracer_halt(substr, status)
       CALL set_tracer(status, GPTRSTR, idt_ms4cs, I_AEROSOL_METHOD,  MODAL)
       CALL tracer_halt(substr, status)
       CALL set_tracer(status, GPTRSTR, idt_ms4cs, S_AEROSOL_MODEL,   modstr)
       CALL tracer_halt(substr, status)
       CALL set_tracer(status, GPTRSTR, idt_ms4cs, I_DRYDEP, ON)
       CALL tracer_halt(substr, status)
       CALL set_tracer(status, GPTRSTR, idt_ms4cs, I_SEDI,   ON)
       CALL tracer_halt(substr, status)
       CALL set_tracer(status, GPTRSTR, idt_ms4cs, I_SCAV,   ON)
       CALL tracer_halt(substr, status)
#ifdef MESSYTENDENCY
    CALL mtend_register(my_handle, idt_ms4cs)
#endif
    endif

    if (l_notrac_OCBC) then
       idt_mbccs = idt_dummy
    else
    CALL new_tracer(status, GPTRSTR, 'BC', modstr , subname='cs'     &
         , idx=idt_mbccs , unit='mol/mol', medium=AEROSOL  &
         , quantity=AMOUNTFRACTION )
    CALL tracer_halt(substr, status)
    CALL set_tracer(status, GPTRSTR, idt_mbccs, R_molarmass,       r=12.011_dp)
    CALL tracer_halt(substr, status)
    CALL set_tracer(status, GPTRSTR, idt_mbccs, R_AEROSOL_DENSITY, r=2000._dp)
    CALL tracer_halt(substr, status)
    CALL set_tracer(status, GPTRSTR, idt_mbccs, I_AEROSOL_MODE,    4)
    CALL tracer_halt(substr, status)
    CALL set_tracer(status, GPTRSTR, idt_mbccs, I_AEROSOL_METHOD,  MODAL)
    CALL tracer_halt(substr, status)
    CALL set_tracer(status, GPTRSTR, idt_mbccs, S_AEROSOL_MODEL,   modstr)
    CALL tracer_halt(substr, status)
    CALL set_tracer(status, GPTRSTR, idt_mbccs, I_DRYDEP, ON)
    CALL tracer_halt(substr, status)
    CALL set_tracer(status, GPTRSTR, idt_mbccs, I_SEDI,   ON)
    CALL tracer_halt(substr, status)
    CALL set_tracer(status, GPTRSTR, idt_mbccs, I_SCAV,   ON)
    CALL tracer_halt(substr, status)
#ifdef MESSYTENDENCY
    CALL mtend_register(my_handle, idt_mbccs)
#endif
    endif

    if (l_notrac_OCBC) then
       idt_moccs = idt_dummy
    else
    CALL new_tracer(status, GPTRSTR, 'OC', modstr , subname='cs'     &
         , idx=idt_moccs , unit='mol/mol', medium=AEROSOL  &
         , quantity=AMOUNTFRACTION )
    CALL tracer_halt(substr, status)
    CALL set_tracer(status, GPTRSTR, idt_moccs, R_molarmass,       r=12.011_dp)
    CALL tracer_halt(substr, status)
    CALL set_tracer(status, GPTRSTR, idt_moccs, R_AEROSOL_DENSITY, r=2000._dp)
    CALL tracer_halt(substr, status)
    CALL set_tracer(status, GPTRSTR, idt_moccs, I_AEROSOL_MODE,    4)
    CALL tracer_halt(substr, status)
    CALL set_tracer(status, GPTRSTR, idt_moccs, I_AEROSOL_METHOD,  MODAL)
    CALL tracer_halt(substr, status)
    CALL set_tracer(status, GPTRSTR, idt_moccs, S_AEROSOL_MODEL,   modstr)
    CALL tracer_halt(substr, status)
    CALL set_tracer(status, GPTRSTR, idt_moccs, I_DRYDEP, ON)
    CALL tracer_halt(substr, status)
    CALL set_tracer(status, GPTRSTR, idt_moccs, I_SEDI,   ON)
    CALL tracer_halt(substr, status)
    CALL set_tracer(status, GPTRSTR, idt_moccs, I_SCAV,   ON)
    CALL tracer_halt(substr, status)
#ifdef MESSYTENDENCY
    CALL mtend_register(my_handle, idt_moccs)
#endif
    endif

    if (l_notrac_SS) then
       idt_msscs = idt_dummy
    else
    CALL new_tracer(status, GPTRSTR, 'SS', modstr , subname='cs'     &
         , idx=idt_msscs , unit='mol/mol', medium=AEROSOL  &
         , quantity=AMOUNTFRACTION )
    CALL tracer_halt(substr, status)
    CALL set_tracer(status, GPTRSTR, idt_msscs, R_molarmass,       r=58.44_dp)
    CALL tracer_halt(substr, status)
    CALL set_tracer(status, GPTRSTR, idt_msscs, R_AEROSOL_DENSITY, r=2165._dp)
    CALL tracer_halt(substr, status)
    CALL set_tracer(status, GPTRSTR, idt_msscs, I_AEROSOL_MODE,    4)
    CALL tracer_halt(substr, status)
    CALL set_tracer(status, GPTRSTR, idt_msscs, I_AEROSOL_METHOD,  MODAL)
    CALL tracer_halt(substr, status)
    CALL set_tracer(status, GPTRSTR, idt_msscs, S_AEROSOL_MODEL,   modstr)
    CALL tracer_halt(substr, status)
    CALL set_tracer(status, GPTRSTR, idt_msscs, I_DRYDEP, ON)
    CALL tracer_halt(substr, status)
    CALL set_tracer(status, GPTRSTR, idt_msscs, I_SEDI,   ON)
    CALL tracer_halt(substr, status)
    CALL set_tracer(status, GPTRSTR, idt_msscs, I_SCAV,   ON)
    CALL tracer_halt(substr, status)
#ifdef MESSYTENDENCY
    CALL mtend_register(my_handle, idt_msscs)
#endif
    endif

    if (l_notrac_DU .or. l_notrac_SO4) then
       idt_mducs = idt_dummy
    else
    CALL new_tracer(status, GPTRSTR, 'DU', modstr , subname='cs'     &
         , idx=idt_mducs , unit='mol/mol', medium=AEROSOL  &
         , quantity=AMOUNTFRACTION )
    CALL tracer_halt(substr, status)
    CALL set_tracer(status, GPTRSTR, idt_mducs, R_molarmass,       r=40.08_dp)
    CALL tracer_halt(substr, status)
    CALL set_tracer(status, GPTRSTR, idt_mducs, R_AEROSOL_DENSITY, r=2650._dp)
    CALL tracer_halt(substr, status)
    CALL set_tracer(status, GPTRSTR, idt_mducs, I_AEROSOL_MODE,    4)
    CALL tracer_halt(substr, status)
    CALL set_tracer(status, GPTRSTR, idt_mducs, I_AEROSOL_METHOD,  MODAL)
    CALL tracer_halt(substr, status)
    CALL set_tracer(status, GPTRSTR, idt_mducs, S_AEROSOL_MODEL,   modstr)
    CALL tracer_halt(substr, status)
    CALL set_tracer(status, GPTRSTR, idt_mducs, I_DRYDEP, ON)
    CALL tracer_halt(substr, status)
    CALL set_tracer(status, GPTRSTR, idt_mducs, I_SEDI,   ON)
    CALL tracer_halt(substr, status)
    CALL set_tracer(status, GPTRSTR, idt_mducs, I_SCAV,   ON)
    CALL tracer_halt(substr, status)
#ifdef MESSYTENDENCY
    CALL mtend_register(my_handle, idt_mducs)
#endif
    endif
    
    !--- Mode 5 - Aitken Insoluble:

    if (l_notrac_OCBC) then
       idt_mbcki = idt_dummy
    else
    CALL new_tracer(status, GPTRSTR, 'BC', modstr , subname='ki'     &
         , idx=idt_mbcki , unit='mol/mol', medium=AEROSOL  &
         , quantity=AMOUNTFRACTION )
    CALL tracer_halt(substr, status)
    CALL set_tracer(status, GPTRSTR, idt_mbcki, R_molarmass,       r=12.011_dp)
    CALL tracer_halt(substr, status)
    CALL set_tracer(status, GPTRSTR, idt_mbcki, R_AEROSOL_DENSITY, r=2000._dp)
    CALL tracer_halt(substr, status)
    CALL set_tracer(status, GPTRSTR, idt_mbcki, I_AEROSOL_MODE,    5)
    CALL tracer_halt(substr, status)
    CALL set_tracer(status, GPTRSTR, idt_mbcki, I_AEROSOL_METHOD,  MODAL)
    CALL tracer_halt(substr, status)
    CALL set_tracer(status, GPTRSTR, idt_mbcki, S_AEROSOL_MODEL,   modstr)
    CALL tracer_halt(substr, status)
    CALL set_tracer(status, GPTRSTR, idt_mbcki, I_DRYDEP, ON)
    CALL tracer_halt(substr, status)
    CALL set_tracer(status, GPTRSTR, idt_mbcki, I_SEDI,   ON)
    CALL tracer_halt(substr, status)
    CALL set_tracer(status, GPTRSTR, idt_mbcki, I_SCAV,   ON)
    CALL tracer_halt(substr, status)
    CALL set_tracer(status, GPTRSTR, idt_mbcki, I_AEROSOL_SOL,   OFF)
    CALL tracer_halt(substr, status)
#ifdef MESSYTENDENCY
    CALL mtend_register(my_handle, idt_mbcki)
#endif
    endif

    if (l_notrac_OCBC) then
       idt_mocki = idt_dummy
    else
    CALL new_tracer(status, GPTRSTR, 'OC', modstr , subname='ki'     &
         , idx=idt_mocki , unit='mol/mol', medium=AEROSOL  &
         , quantity=AMOUNTFRACTION )
    CALL tracer_halt(substr, status)
    CALL set_tracer(status, GPTRSTR, idt_mocki, R_molarmass,       r=12.011_dp)
    CALL tracer_halt(substr, status)
    CALL set_tracer(status, GPTRSTR, idt_mocki, R_AEROSOL_DENSITY, r=2000._dp)
    CALL tracer_halt(substr, status)
    CALL set_tracer(status, GPTRSTR, idt_mocki, I_AEROSOL_MODE,    5)
    CALL tracer_halt(substr, status)
    CALL set_tracer(status, GPTRSTR, idt_mocki, I_AEROSOL_METHOD,  MODAL)
    CALL tracer_halt(substr, status)
    CALL set_tracer(status, GPTRSTR, idt_mocki, S_AEROSOL_MODEL,   modstr)
    CALL tracer_halt(substr, status)
    CALL set_tracer(status, GPTRSTR, idt_mocki, I_DRYDEP, ON)
    CALL tracer_halt(substr, status)
    CALL set_tracer(status, GPTRSTR, idt_mocki, I_SEDI,   ON)
    CALL tracer_halt(substr, status)
    CALL set_tracer(status, GPTRSTR, idt_mocki, I_SCAV,   ON)
    CALL tracer_halt(substr, status)
    CALL set_tracer(status, GPTRSTR, idt_mocki, I_AEROSOL_SOL,   OFF)
    CALL tracer_halt(substr, status)
#ifdef MESSYTENDENCY
    CALL mtend_register(my_handle, idt_mocki)
#endif
    endif

    !--- Mode 6 - Accumulation Insoluble:
    if (l_notrac_DU) then
       idt_mduai = idt_dummy
    else
    CALL new_tracer(status, GPTRSTR, 'DU', modstr , subname='ai'     &
         , idx=idt_mduai , unit='mol/mol', medium=AEROSOL  &
         , quantity=AMOUNTFRACTION )
    CALL tracer_halt(substr, status)
    CALL set_tracer(status, GPTRSTR, idt_mduai, R_molarmass,       r=40.08_dp)
    CALL tracer_halt(substr, status)
    CALL set_tracer(status, GPTRSTR, idt_mduai, R_AEROSOL_DENSITY, r=2650._dp)
    CALL tracer_halt(substr, status)
    CALL set_tracer(status, GPTRSTR, idt_mduai, I_AEROSOL_MODE,    6)
    CALL tracer_halt(substr, status)
    CALL set_tracer(status, GPTRSTR, idt_mduai, I_AEROSOL_METHOD,  MODAL)
    CALL tracer_halt(substr, status)
    CALL set_tracer(status, GPTRSTR, idt_mduai, S_AEROSOL_MODEL,   modstr)
    CALL tracer_halt(substr, status)
    CALL set_tracer(status, GPTRSTR, idt_mduai, I_DRYDEP, ON)
    CALL tracer_halt(substr, status)
    CALL set_tracer(status, GPTRSTR, idt_mduai, I_SEDI,   ON)
    CALL tracer_halt(substr, status)
    CALL set_tracer(status, GPTRSTR, idt_mduai, I_SCAV,   ON)
    CALL tracer_halt(substr, status)
    CALL set_tracer(status, GPTRSTR, idt_mduai, I_AEROSOL_SOL,   OFF)
    CALL tracer_halt(substr, status)
#ifdef MESSYTENDENCY
    CALL mtend_register(my_handle, idt_mduai)
#endif
    endif

    !--- Mode 7 - Coarse Insoluble:
    if (l_notrac_DU) then
       idt_mduci = idt_dummy
    else
    CALL new_tracer(status, GPTRSTR, 'DU', modstr , subname='ci'     &
         , idx=idt_mduci , unit='mol/mol', medium=AEROSOL  &
         , quantity=AMOUNTFRACTION )
    CALL tracer_halt(substr, status)
    CALL set_tracer(status, GPTRSTR, idt_mduci, R_molarmass,       r=40.08_dp)
    CALL tracer_halt(substr, status)
    CALL set_tracer(status, GPTRSTR, idt_mduci, R_AEROSOL_DENSITY, r=2650._dp)
    CALL tracer_halt(substr, status)
    CALL set_tracer(status, GPTRSTR, idt_mduci, I_AEROSOL_MODE,    7)
    CALL tracer_halt(substr, status)
    CALL set_tracer(status, GPTRSTR, idt_mduci, I_AEROSOL_METHOD,  MODAL)
    CALL tracer_halt(substr, status)
    CALL set_tracer(status, GPTRSTR, idt_mduci, S_AEROSOL_MODEL,   modstr)
    CALL tracer_halt(substr, status)
    CALL set_tracer(status, GPTRSTR, idt_mduci, I_DRYDEP, ON)
    CALL tracer_halt(substr, status)
    CALL set_tracer(status, GPTRSTR, idt_mduci, I_SEDI,   ON)
    CALL tracer_halt(substr, status)
    CALL set_tracer(status, GPTRSTR, idt_mduci, I_SCAV,   ON)
    CALL tracer_halt(substr, status)
    CALL set_tracer(status, GPTRSTR, idt_mduci, I_AEROSOL_SOL,   OFF)
    CALL tracer_halt(substr, status)
#ifdef MESSYTENDENCY
    CALL mtend_register(my_handle, idt_mduci)
#endif
    endif

    !--- 3) Aerosol Numbers:
    ! NUM renamed to N, because this name is hardcoded in EMDEP 
    ! otherwise no emissions
    ! unit has to be written 1/mol not 1 cm-3 also caused by EMDEP
    ! choose unit as 1/kg(air) or 1/mol(air)

    if (l_notrac_SO4) then
       idt_nns = idt_dummy
    else
    CALL new_tracer(status, GPTRSTR, 'N', modstr , subname='ns'     &
         , idx=idt_nns , unit='1/mol', medium=AEROSOL  &
         , quantity=NUMBERDENSITY )
    CALL tracer_halt(substr, status)
    CALL set_tracer(status, GPTRSTR, idt_nns, R_AEROSOL_DENSITY, r=1._dp)
    CALL tracer_halt(substr, status)
    CALL set_tracer(status, GPTRSTR, idt_nns, I_AEROSOL_MODE,    1)
    CALL tracer_halt(substr, status)
    CALL set_tracer(status, GPTRSTR, idt_nns, I_AEROSOL_METHOD,  MODAL)
    CALL tracer_halt(substr, status)
    CALL set_tracer(status, GPTRSTR, idt_nns, S_AEROSOL_MODEL,   modstr)
    CALL tracer_halt(substr, status)
    CALL set_tracer(status, GPTRSTR, idt_nns, I_DRYDEP, ON)
    CALL tracer_halt(substr, status)
    CALL set_tracer(status, GPTRSTR, idt_nns, I_SEDI,   ON)
    CALL tracer_halt(substr, status)
    CALL set_tracer(status, GPTRSTR, idt_nns, I_SCAV,   ON)
    CALL tracer_halt(substr, status)
#ifdef MESSYTENDENCY
    CALL mtend_register(my_handle, idt_nns)
#endif
    endif

    if (l_notrac_SO4 .and. l_notrac_OCBC) then
       idt_nks = idt_dummy
    else
    CALL new_tracer(status, GPTRSTR, 'N', modstr , subname='ks'     &
         , idx=idt_nks , unit='1/mol', medium=AEROSOL             &
         , quantity=NUMBERDENSITY )
    CALL tracer_halt(substr, status)
    CALL set_tracer(status, GPTRSTR, idt_nks, R_AEROSOL_DENSITY, r=1._dp)
    CALL tracer_halt(substr, status)
    CALL set_tracer(status, GPTRSTR, idt_nks, I_AEROSOL_MODE,    2)
    CALL tracer_halt(substr, status)
    CALL set_tracer(status, GPTRSTR, idt_nks, I_AEROSOL_METHOD,  MODAL)
    CALL tracer_halt(substr, status)
    CALL set_tracer(status, GPTRSTR, idt_nks, S_AEROSOL_MODEL,   modstr)
    CALL tracer_halt(substr, status)
    CALL set_tracer(status, GPTRSTR, idt_nks, I_DRYDEP, ON)
    CALL tracer_halt(substr, status)
    CALL set_tracer(status, GPTRSTR, idt_nks, I_SEDI,   ON)
    CALL tracer_halt(substr, status)
    CALL set_tracer(status, GPTRSTR, idt_nks, I_SCAV,   ON)
    CALL tracer_halt(substr, status)
#ifdef MESSYTENDENCY
    CALL mtend_register(my_handle, idt_nks)
#endif
    endif
    if (l_notrac_SO4 .and. l_notrac_OCBC .and. l_notrac_SS) then
       idt_nas = idt_dummy
    else
    CALL new_tracer(status, GPTRSTR, 'N', modstr , subname='as'  &
         , idx=idt_nas , unit='1/mol', medium=AEROSOL            &
         , quantity=NUMBERDENSITY )
    CALL tracer_halt(substr, status)
    CALL set_tracer(status, GPTRSTR, idt_nas, R_AEROSOL_DENSITY, r=1._dp)
    CALL tracer_halt(substr, status)
    CALL set_tracer(status, GPTRSTR, idt_nas, I_AEROSOL_MODE,    3)
    CALL tracer_halt(substr, status)
    CALL set_tracer(status, GPTRSTR, idt_nas, I_AEROSOL_METHOD,  MODAL)
    CALL tracer_halt(substr, status)
    CALL set_tracer(status, GPTRSTR, idt_nas, S_AEROSOL_MODEL,   modstr)
    CALL tracer_halt(substr, status)
    CALL set_tracer(status, GPTRSTR, idt_nas, I_DRYDEP, ON)
    CALL tracer_halt(substr, status)
    CALL set_tracer(status, GPTRSTR, idt_nas, I_SEDI,   ON)
    CALL tracer_halt(substr, status)
    CALL set_tracer(status, GPTRSTR, idt_nas, I_SCAV,   ON)
    CALL tracer_halt(substr, status)
#ifdef MESSYTENDENCY
    CALL mtend_register(my_handle, idt_nas)
#endif
    endif

    if (l_notrac_SO4 .and. l_notrac_OCBC .and. l_notrac_SS) then
       idt_ncs = idt_dummy
    else
    CALL new_tracer(status, GPTRSTR, 'N', modstr , subname='cs'     &
         , idx=idt_ncs , unit='1/mol', medium=AEROSOL  &
         , quantity=NUMBERDENSITY )
    CALL tracer_halt(substr, status)
    CALL set_tracer(status, GPTRSTR, idt_ncs, R_AEROSOL_DENSITY, r=1._dp)
    CALL tracer_halt(substr, status)
    CALL set_tracer(status, GPTRSTR, idt_ncs, I_AEROSOL_MODE,    4)
    CALL tracer_halt(substr, status)
    CALL set_tracer(status, GPTRSTR, idt_ncs, I_AEROSOL_METHOD,  MODAL)
    CALL tracer_halt(substr, status)
    CALL set_tracer(status, GPTRSTR, idt_ncs, S_AEROSOL_MODEL,   modstr)
    CALL tracer_halt(substr, status)
    CALL set_tracer(status, GPTRSTR, idt_ncs, I_DRYDEP, ON)
    CALL tracer_halt(substr, status)
    CALL set_tracer(status, GPTRSTR, idt_ncs, I_SEDI,   ON)
    CALL tracer_halt(substr, status)
    CALL set_tracer(status, GPTRSTR, idt_ncs, I_SCAV,   ON)
    CALL tracer_halt(substr, status)
#ifdef MESSYTENDENCY
    CALL mtend_register(my_handle, idt_ncs)
#endif
    endif

    if (l_notrac_OCBC) then
       idt_nki = idt_dummy
    else
    CALL new_tracer(status, GPTRSTR, 'N', modstr , subname='ki'     &
         , idx=idt_nki , unit='1/mol', medium=AEROSOL  &
         , quantity=NUMBERDENSITY )
    CALL tracer_halt(substr, status)
    CALL set_tracer(status, GPTRSTR, idt_nki, R_AEROSOL_DENSITY, r=1._dp)
    CALL tracer_halt(substr, status)
    CALL set_tracer(status, GPTRSTR, idt_nki, I_AEROSOL_MODE,    5)
    CALL tracer_halt(substr, status)
    CALL set_tracer(status, GPTRSTR, idt_nki, I_AEROSOL_METHOD,  MODAL)
    CALL tracer_halt(substr, status)
    CALL set_tracer(status, GPTRSTR, idt_nki, S_AEROSOL_MODEL,   modstr)
    CALL tracer_halt(substr, status)
    CALL set_tracer(status, GPTRSTR, idt_nki, I_DRYDEP, ON)
    CALL tracer_halt(substr, status)
    CALL set_tracer(status, GPTRSTR, idt_nki, I_SEDI,   ON)
    CALL tracer_halt(substr, status)
    CALL set_tracer(status, GPTRSTR, idt_nki, I_SCAV,   ON)
    CALL tracer_halt(substr, status)
    CALL set_tracer(status, GPTRSTR, idt_nki, I_AEROSOL_SOL,   OFF)
    CALL tracer_halt(substr, status)
#ifdef MESSYTENDENCY
    CALL mtend_register(my_handle, idt_nki)
#endif
    endif

    if (l_notrac_OCBC .and. l_notrac_DU) then
       idt_nai = idt_dummy
    else
    CALL new_tracer(status, GPTRSTR, 'N', modstr , subname='ai'     &
         , idx=idt_nai , unit='1/mol', medium=AEROSOL  &
         , quantity=NUMBERDENSITY )
    CALL tracer_halt(substr, status)
    CALL set_tracer(status, GPTRSTR, idt_nai, R_AEROSOL_DENSITY, r=1._dp)
    CALL tracer_halt(substr, status)
    CALL set_tracer(status, GPTRSTR, idt_nai, I_AEROSOL_MODE,    6)
    CALL tracer_halt(substr, status)
    CALL set_tracer(status, GPTRSTR, idt_nai, I_AEROSOL_METHOD,  MODAL)
    CALL tracer_halt(substr, status)
    CALL set_tracer(status, GPTRSTR, idt_nai, S_AEROSOL_MODEL,   modstr)
    CALL tracer_halt(substr, status)
    CALL set_tracer(status, GPTRSTR, idt_nai, I_DRYDEP, ON)
    CALL tracer_halt(substr, status)
    CALL set_tracer(status, GPTRSTR, idt_nai, I_SEDI,   ON)
    CALL tracer_halt(substr, status)
    CALL set_tracer(status, GPTRSTR, idt_nai, I_SCAV,   ON)
    CALL tracer_halt(substr, status)
    CALL set_tracer(status, GPTRSTR, idt_nai, I_AEROSOL_SOL,   OFF)
    CALL tracer_halt(substr, status)
    endif

    if (l_notrac_OCBC .and. l_notrac_DU) then
       idt_nci = idt_dummy
    else
    CALL new_tracer(status, GPTRSTR, 'N', modstr , subname='ci'     &
         , idx=idt_nci , unit='1/mol', medium=AEROSOL  &
         , quantity=NUMBERDENSITY )
    CALL tracer_halt(substr, status)
    CALL set_tracer(status, GPTRSTR, idt_nci, R_AEROSOL_DENSITY, r=1._dp)
    CALL tracer_halt(substr, status)
    CALL set_tracer(status, GPTRSTR, idt_nci, I_AEROSOL_MODE,    7)
    CALL tracer_halt(substr, status)
    CALL set_tracer(status, GPTRSTR, idt_nci, I_AEROSOL_METHOD,  MODAL)
    CALL tracer_halt(substr, status)
    CALL set_tracer(status, GPTRSTR, idt_nci, S_AEROSOL_MODEL,   modstr)
    CALL tracer_halt(substr, status)
    CALL set_tracer(status, GPTRSTR, idt_nci, I_DRYDEP, ON)
    CALL tracer_halt(substr, status)
    CALL set_tracer(status, GPTRSTR, idt_nci, I_SEDI,   ON)
    CALL tracer_halt(substr, status)
    CALL set_tracer(status, GPTRSTR, idt_nci, I_SCAV,   ON)
    CALL tracer_halt(substr, status)
    CALL set_tracer(status, GPTRSTR, idt_nci, I_AEROSOL_SOL,   OFF)
    CALL tracer_halt(substr, status)
#ifdef MESSYTENDENCY
    CALL mtend_register(my_handle, idt_nci)
#endif
    endif

    CALL end_message_bi(modstr,'REQUEST M7 TRACER', substr)

  END SUBROUTINE m7_new_tracer
!-----------------------------------------------------------------------------

!-----------------------------------------------------------------------------
  SUBROUTINE m7_init_memory

    ! ECHAM5/MESSy
    USE messy_main_grid_def_mem_bi,  ONLY: nproma, ngpblks, nlev
    USE messy_main_channel_error_bi, ONLY: channel_halt
    USE messy_main_channel_bi,       ONLY: GP_3D_MID, DC_GP                &
                                         , DIMID_LON, DIMID_LAT, DIMID_LEV &
                                         , DC_BC                           &
                                         , gp_nseg, gp_start, gp_cnt       &
                                         , gp_meml, gp_memu                &
                                         , GP_2D_HORIZONTAL
    ! MESSy
    USE messy_main_channel,          ONLY: new_channel, new_channel_object &
                                         , new_attribute
    USE messy_main_channel_dimensions, ONLY: new_dimension
    USE messy_main_channel_repr,       ONLY: new_representation, AUTO  &
                                           , set_representation_decomp &
                                           , IRANK, PIOTYPE_COL        &
                                           , repr_def_axes
    USE messy_main_tracer_mem_bi,      ONLY: GPTRSTR

    IMPLICIT NONE

    ! LOCAL
    CHARACTER(len=*), PARAMETER :: substr='m7_init_memory'
    ! auxiliary pointer for data managment
    INTEGER                               :: jmod, jsol, status
    REAL(dp), DIMENSION(:,:,:,:),   POINTER ::  mem => NULL()
    INTEGER                               :: DIMID_NMODE
    INTEGER                               :: REPR_M7_4D_NMOD
    INTEGER                               :: REPR_M7_1D
    CHARACTER(LEN=1)                      :: ichar

    ! PARALLEL DECOMPOSITION
    INTEGER                          :: nseg = 0
    INTEGER, DIMENSION(:,:), POINTER :: start => NULL()
    INTEGER, DIMENSION(:,:), POINTER :: cnt   => NULL()
    INTEGER, DIMENSION(:,:), POINTER :: meml  => NULL()
    INTEGER, DIMENSION(:,:), POINTER :: memu  => NULL()

    CHARACTER(LEN=*), PARAMETER :: modstr_gp = modstr//'_'//GPTRSTR

    CALL start_message_bi(modstr,'INITIALIZE MEMORY', substr)

    ALLOCATE(wetradius(_RI_XYZN_(nproma,ngpblks,nlev,nmod),1))
    wetradius(:,:,:,:,:) = 0._dp
    ALLOCATE(dryradius(_RI_XYZN_(nproma,ngpblks,nlev,nmod),1))
    dryradius(:,:,:,:,:) = 0._dp
    ALLOCATE(densaer(_RI_XYZN_(nproma,ngpblks,nlev,nmod),1))
    densaer(:,:,:,:,:)   = 0._dp

    !--- 1) Construct the m7 channel: --------------------------------------!
    CALL new_dimension(status, DIMID_NMODE, 'M7_NMODE', nmod)
    CALL channel_halt(substr, status)

    ! NEW REPRESENTATIONS
    CALL new_representation(status, REPR_M7_4D_NMOD, &
         'REPR_M7_4D_NMOD', rank = 4, link = 'xxxx', dctype = DC_GP   &
         , dimension_ids = (/ &
           _RI_XYZN_(DIMID_LON, DIMID_LAT, DIMID_LEV, DIMID_NMODE) /) &
           , ldimlen       = (/ &
           _RI_XYZN_(nproma, ngpblks, AUTO, AUTO) /)       &
           , output_order  = (/ _IN_XYZN_, _IX_XYZN_        & ! E: 3,1,4,2
                              , _IY_XYZN_, _IZ_XYZN_ /)     & ! C: 3,1,2,4
         , axis = repr_def_axes(_RI_XYZN_('X','Y','Z','N')) &
         )    
    CALL channel_halt(substr, status)

    nseg = gp_nseg
    ALLOCATE(start(nseg,IRANK))
    ALLOCATE(cnt(nseg,IRANK))
    ALLOCATE(meml(nseg,IRANK))
    ALLOCATE(memu(nseg,IRANK))
    
    start(:,:) = gp_start(:,:)
    cnt(:,:) = gp_cnt(:,:)
    meml(:,:) = gp_meml(:,:)
    memu(:,:) = gp_memu(:,:)
    
    cnt(:,_IN_XYZN_) = nmod
    memu(:,_IN_XYZN_) = nmod
    
    CALL set_representation_decomp(status, REPR_M7_4D_NMOD &
         , start, cnt, memu, meml, .FALSE., PIOTYPE_COL)
    CALL channel_halt(substr, status)
    
    DEALLOCATE(start) ; NULLIFY(start)
    DEALLOCATE(cnt)   ; NULLIFY(cnt)
    DEALLOCATE(meml)  ; NULLIFY(meml)
    DEALLOCATE(memu)  ; NULLIFY(memu)

    CALL new_representation(status, REPR_M7_1D,     &
         'REPR_M7_1D'                               &
         , rank = 1, link = 'x---', dctype = DC_BC  &
         , dimension_ids = (/ DIMID_NMODE /)        &
         , ldimlen       = (/ AUTO /)               &
         , axis = 'N---'                            &
         )
    CALL channel_halt(substr, status)

    nseg = 1
    ALLOCATE(start(nseg,IRANK))
    ALLOCATE(cnt(nseg,IRANK))
    ALLOCATE(meml(nseg,IRANK))
    ALLOCATE(memu(nseg,IRANK))
    
    start(:,:) = 1
    cnt(:,:) = 1
    meml(:,:) = 1
    memu(:,:) = 1
    
    start(:,1) = 1
    cnt(:,1)   = nmod
    meml(:,1)  = 1
    memu(:,1)  = nmod
    
    CALL set_representation_decomp(status, REPR_M7_1D &
         , start, cnt, memu, meml, .FALSE., PIOTYPE_COL)
    CALL channel_halt(substr, status)
    
    DEALLOCATE(start) ; NULLIFY(start)
    DEALLOCATE(cnt)   ; NULLIFY(cnt)
    DEALLOCATE(meml)  ; NULLIFY(meml)
    DEALLOCATE(memu)  ; NULLIFY(memu)

    CALL new_channel(status, modstr_gp, reprid=GP_3D_MID)
    CALL channel_halt(substr, status)

    !--- 2) Aerosol Properties:
    ! WET RADII
    mem => wetradius(:,:,:,:,1)
    CALL new_channel_object(status, modstr_gp, 'wetradius' &
         , p4=wetrad_4d, lrestreq=.FALSE., reprid=REPR_M7_4D_NMOD, mem=mem)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr_gp, 'wetradius' &
         , 'units', c='m')
    CALL channel_halt(substr, status)

    DO jmod=1, nmod
       mem => wetradius(_RI_XYZN_(:,:,:,jmod),:)
       WRITE(ichar,'(i1)') jmod
       CALL new_channel_object(status, modstr_gp, 'wetrad_M'//ichar &
            , mem=mem, lrestreq=.TRUE.)
       CALL channel_halt(substr, status)
       CALL new_attribute(status, modstr_gp, 'wetrad_M'//ichar &
            , 'units', c='m')
       CALL channel_halt(substr, status)
    END DO

    ! DRY RADII
    mem => dryradius(:,:,:,:,1)
    CALL new_channel_object(status, modstr_gp, 'dryradius' &
         , p4=dryrad_4d, lrestreq=.FALSE., reprid=REPR_M7_4D_NMOD, mem=mem)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr_gp, 'dryradius' &
         , 'units', c='m')
    CALL channel_halt(substr, status)

    DO jsol=1, nsol
       mem => dryradius(_RI_XYZN_(:,:,:,jsol),:)
       WRITE(ichar,'(i1)') jsol
       CALL new_channel_object(status, modstr_gp, 'dryrad_M'//ichar &
            , mem=mem, lrestreq=.TRUE.)
       CALL channel_halt(substr, status)
       CALL new_attribute(status, modstr_gp, 'dryrad_M'//ichar &
            , 'units', c='m')
       CALL channel_halt(substr, status)
    END DO

    ! DENSITY
    mem => densaer(:,:,:,:,1)
    CALL new_channel_object(status, modstr_gp, 'densaer' &
         , p4=densaer_4d, lrestreq=.FALSE., reprid=REPR_M7_4D_NMOD, mem=mem)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr_gp, 'densaer' &
         , 'units', c='kg m-3')
    CALL channel_halt(substr, status)

    DO jmod=1, nmod
       mem => densaer(_RI_XYZN_(:,:,:,jmod),:)
       WRITE(ichar,'(i1)') jmod
       CALL new_channel_object(status, modstr_gp, 'densaer_M'//ichar &
            , mem=mem, lrestreq=.TRUE.)
       CALL channel_halt(substr, status)
       CALL new_attribute(status, modstr_gp, 'densaer_M'//ichar &
            , 'units', c='kg m-3')
       CALL channel_halt(substr, status)
    END DO

    ! AEROSOL COLUMN MASSES
    IF (l_colmass) THEN
       CALL new_channel_object(status, modstr_gp, 'colmass_SO4', &
            p2=colmass_SO4, reprid=GP_2D_HORIZONTAL )
       CALL channel_halt(substr, status)
       CALL new_attribute(status, modstr_gp, &
            'colmass_SO4', 'long_name', c='column mass of all SO4 modes' )
       CALL channel_halt(substr, status)
       CALL new_attribute(status, modstr_gp, 'colmass_SO4', &
            'units', c='mg m-2' )
       CALL channel_halt(substr, status)

       CALL new_channel_object(status, modstr_gp, 'colmass_BC', &
            p2=colmass_BC, reprid=GP_2D_HORIZONTAL )
       CALL channel_halt(substr, status)
       CALL new_attribute(status, modstr_gp, &
            'colmass_BC', 'long_name', c='column mass of all BC modes' )
       CALL channel_halt(substr, status)
       CALL new_attribute(status, modstr_gp, 'colmass_BC', &
            'units', c='mg m-2' )
       CALL channel_halt(substr, status)

       CALL new_channel_object(status, modstr_gp, 'colmass_OC', &
            p2=colmass_OC, reprid=GP_2D_HORIZONTAL )
       CALL channel_halt(substr, status)
       CALL new_attribute(status, modstr_gp, &
            'colmass_OC', 'long_name', c='column mass of all OC modes' )
       CALL channel_halt(substr, status)
       CALL new_attribute(status, modstr_gp, 'colmass_OC', &
            'units', c='mg m-2' )
       CALL channel_halt(substr, status)

       CALL new_channel_object(status, modstr_gp, 'colmass_SS', &
            p2=colmass_SS, reprid=GP_2D_HORIZONTAL )
       CALL channel_halt(substr, status)
       CALL new_attribute(status, modstr_gp, &
            'colmass_SS', 'long_name', c='column mass of all SS modes' )
       CALL channel_halt(substr, status)
       CALL new_attribute(status, modstr_gp, 'colmass_SS', &
            'units', c='mg m-2' )
       CALL channel_halt(substr, status)

       CALL new_channel_object(status, modstr_gp, 'colmass_DU', &
            p2=colmass_DU, reprid=GP_2D_HORIZONTAL )
       CALL channel_halt(substr, status)
       CALL new_attribute(status, modstr_gp, &
            'colmass_DU', 'long_name', c='column mass of all DU modes' )
       CALL channel_halt(substr, status)
       CALL new_attribute(status, modstr_gp, 'colmass_DU', &
            'units', c='mg m-2' )
       CALL channel_halt(substr, status)

       CALL new_channel_object(status, modstr_gp, 'colmass_DU_as', &
            p2=colmass_DU_as, reprid=GP_2D_HORIZONTAL )
       CALL channel_halt(substr, status)
       CALL new_attribute(status, modstr_gp, &
            'colmass_DU_as', 'long_name', c='column mass of DU_as' )
       CALL channel_halt(substr, status)
       CALL new_attribute(status, modstr_gp, 'colmass_DU_as', &
            'units', c='mg m-2' )
       CALL channel_halt(substr, status)

       CALL new_channel_object(status, modstr_gp, 'colmass_DU_cs', &
            p2=colmass_DU_cs, reprid=GP_2D_HORIZONTAL )
       CALL channel_halt(substr, status)
       CALL new_attribute(status, modstr_gp, &
            'colmass_DU_cs', 'long_name', c='column mass of DU_cs' )
       CALL channel_halt(substr, status)
       CALL new_attribute(status, modstr_gp, 'colmass_DU_cs', &
            'units', c='mg m-2' )
       CALL channel_halt(substr, status)

       CALL new_channel_object(status, modstr_gp, 'colmass_DU_ai', &
            p2=colmass_DU_ai, reprid=GP_2D_HORIZONTAL )
       CALL channel_halt(substr, status)
       CALL new_attribute(status, modstr_gp, &
            'colmass_DU_ai', 'long_name', c='column mass of DU_ai' )
       CALL channel_halt(substr, status)
       CALL new_attribute(status, modstr_gp, 'colmass_DU_ai', &
            'units', c='mg m-2' )
       CALL channel_halt(substr, status)

       CALL new_channel_object(status, modstr_gp, 'colmass_DU_ci', &
            p2=colmass_DU_ci, reprid=GP_2D_HORIZONTAL )
       CALL channel_halt(substr, status)
       CALL new_attribute(status, modstr_gp, &
            'colmass_DU_ci', 'long_name', c='column mass of DU_ci' )
       CALL channel_halt(substr, status)
       CALL new_attribute(status, modstr_gp, 'colmass_DU_ci', &
            'units', c='mg m-2' )
       CALL channel_halt(substr, status)
    ENDIF ! l_colmass

    ! AEROSOL OPTICAL DEPTH
    IF (l_aod) THEN
       CALL new_channel_object(status, modstr_gp, 'AOD_SO4_COLUMN', &
            p2=AOD_SO4_COLUMN, reprid=GP_2D_HORIZONTAL )
       CALL channel_halt(substr, status)
       CALL new_attribute(status, modstr_gp, &
            'AOD_SO4_COLUMN', 'long_name', c='aod (all SO4 modes)' )
       CALL channel_halt(substr, status)
       CALL new_attribute(status, modstr_gp, 'AOD_SO4_COLUMN', &
            'units', c='mg m-2' )
       CALL channel_halt(substr, status)

       CALL new_channel_object(status, modstr_gp, 'AOD_BC_COLUMN', &
            p2=AOD_BC_COLUMN, reprid=GP_2D_HORIZONTAL )
       CALL channel_halt(substr, status)
       CALL new_attribute(status, modstr_gp, &
            'AOD_BC_COLUMN', 'long_name', c='aod (all BC modes)' )
       CALL channel_halt(substr, status)
       CALL new_attribute(status, modstr_gp, 'AOD_BC_COLUMN', &
            'units', c='mg m-2' )
       CALL channel_halt(substr, status)

       CALL new_channel_object(status, modstr_gp, 'AOD_OC_COLUMN', &
            p2=AOD_OC_COLUMN, reprid=GP_2D_HORIZONTAL )
       CALL channel_halt(substr, status)
       CALL new_attribute(status, modstr_gp, &
            'AOD_OC_COLUMN', 'long_name', c='aod (all OC modes)' )
       CALL channel_halt(substr, status)
       CALL new_attribute(status, modstr_gp, 'AOD_OC_COLUMN', &
            'units', c='mg m-2' )
       CALL channel_halt(substr, status)

       CALL new_channel_object(status, modstr_gp, 'AOD_SS_COLUMN', &
            p2=AOD_SS_COLUMN, reprid=GP_2D_HORIZONTAL )
       CALL channel_halt(substr, status)
       CALL new_attribute(status, modstr_gp, &
            'AOD_SS_COLUMN', 'long_name', c='aod (all SS modes)' )
       CALL channel_halt(substr, status)
       CALL new_attribute(status, modstr_gp, 'AOD_SS_COLUMN', &
            'units', c='mg m-2' )
       CALL channel_halt(substr, status)

       CALL new_channel_object(status, modstr_gp, 'AOD_DU_COLUMN', &
            p2=AOD_DU_COLUMN, reprid=GP_2D_HORIZONTAL )
       CALL channel_halt(substr, status)
       CALL new_attribute(status, modstr_gp, &
            'AOD_DU_COLUMN', 'long_name', c='aod (all DU modes)' )
       CALL channel_halt(substr, status)
       CALL new_attribute(status, modstr_gp, 'AOD_DU_COLUMN', &
            'units', c='mg m-2' )
       CALL channel_halt(substr, status)

       CALL new_channel_object(status, modstr_gp, 'AOD_TOTAL_COLUMN', &
            p2=AOD_TOTAL_COLUMN, reprid=GP_2D_HORIZONTAL )
       CALL channel_halt(substr, status)
       CALL new_attribute(status, modstr_gp, &
            'AOD_TOTAL_COLUMN', 'long_name', c='aod (total aerosol)' )
       CALL channel_halt(substr, status)
       CALL new_attribute(status, modstr_gp, 'AOD_TOTAL_COLUMN', &
            'units', c='mg m-2' )
       CALL channel_halt(substr, status)
    ENDIF ! l_aod

    CALL new_channel_object(status, modstr_gp, &
         'sigma', p1=sigma_cha, reprid=REPR_M7_1D)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr_gp, 'sigma', &
         'long_name', c='standard deviation of M7 modes' )
    CALL channel_halt(substr, status)
    
    !--- Initialization of channel objects containing M7 parameter

    sigma_cha    = sigma

    CALL end_message_bi(modstr,'INITIALIZE MEMORY', substr)

  END SUBROUTINE m7_init_memory
!-----------------------------------------------------------------------------

!-----------------------------------------------------------------------------
  SUBROUTINE m7_init_coupling

    ! SUBROUTINE to check for existence of required coupling tracer

    ! Author: Astrid Kerkweg, MPICH, Jan. 2004

    ! ECHAM5/MESSy
    USE messy_main_tracer_mem_bi,   ONLY: GPTRSTR
    USE messy_main_tracer_tools_bi, ONLY: tracer_halt
    ! MESSy
    USE messy_main_tracer,        ONLY: get_tracer
    USE messy_main_channel,       ONLY: get_channel_object

    IMPLICIT NONE

    CHARACTER(LEN=*), PARAMETER :: substr = 'm7_init_coupling'
    LOGICAL                     :: lok      ! error status   
    INTEGER                     :: status

    CALL start_message_bi(modstr,'CHECK NAMELIST SETTINGS', substr)

    IF (lcpl_gasphase) THEN
       IF  (H2SO4_gas(1) == '' .AND. H2SO4_gas(2) == '') THEN
          CALL info_bi('no gas-phase H2SO4 tracer required', substr)
          CALL error_bi('require H2SO4 tracer for gas phase coupling', substr)
       ELSE
          CALL get_tracer(status, GPTRSTR, H2SO4_gas(1) &
               ,subname= H2SO4_gas(2), idx=idt_so4_m7)
          lok = (status == 0)
          IF (lok) THEN
             IF (TRIM(ti_gp(idt_so4_m7)%tp%ident%submodel) == &
                  TRIM(chemmodule)) THEN
                CALL info_bi( &
                     'fetching chemistry module H2SO4 tracer idt', substr)
#ifdef MESSYTENDENCY
    CALL mtend_register(my_handle, idt_so4_m7)
#endif
             ELSE
                CALL error_bi( &
                     'H2SO4 tracer not initialized by chosen chemistry module'&
                     , substr)
             ENDIF
          ELSE
             CALL error_bi(&
                  'no chemistry module gas-phase H2SO4 tracer available' &
                  , substr)
             CALL tracer_halt(substr, status)
          ENDIF
       ENDIF
    ENDIF

    ! couple to MECCA-AERO H2SO4_as
    IF  (H2SO4_as(1) == '' .AND. H2SO4_as(2) == '') THEN
       CALL info_bi('no aerosol-phase H2SO4 tracer required', substr)
    ELSE
       CALL get_tracer(status, GPTRSTR, H2SO4_as(1) &
            ,subname= H2SO4_as(2), idx=idt_h2so4_as)
       lok = (status == 0)
       IF (lok) THEN
          IF (TRIM(ti_gp(idt_h2so4_as)%tp%ident%submodel) == &
               TRIM(chemmodule)) THEN
             CALL info_bi('fetching chemistry module H2SO4_as tracer idt' &
                  , substr)
#ifdef MESSYTENDENCY
    CALL mtend_register(my_handle, idt_h2so4_as)
#endif
          ELSE
             idt_h2so4_as =0 
             CALL info_bi('no MECCA-AERO H2SO4_as tracer initialized by'//&
                  'chosen chemistry module', substr)
          ENDIF
       ELSE
          idt_h2so4_as =0 
          CALL info_bi('no MECCA-AERO H2SO4_as tracer initialized by'//&
               'chosen chemistry module', substr)
       ENDIF
    ENDIF
    ! find idts for the other MECCA-AERO AS tracers   
    CALL get_tracer(status, GPTRSTR, 'HSO4m','as', idx=idt_hso4m_as)
    lok = (status == 0)
    IF (lok) THEN
       IF (TRIM(ti_gp(idt_hso4m_as)%tp%ident%submodel) == &
            TRIM(chemmodule)) THEN
          CALL info_bi('fetching chemistry module HSO4m_as tracer idt', substr)
#ifdef MESSYTENDENCY
    CALL mtend_register(my_handle, idt_hso4m_as)
#endif
       ELSE
          idt_hso4m_as =0 
          CALL info_bi( &
            'no MECCA-AERO HSO4m_as tracer initialized by chosen chem. module'&
          , substr)
       ENDIF
    ELSE
       idt_hso4m_as =0 
       CALL info_bi(&
            'no MECCA-AERO HSO4m_as tracer initialized by chosen chem. module' &
            , substr)
    ENDIF
!-----------------
    CALL get_tracer(status, GPTRSTR, 'SO4mm','as', idx=idt_so4mm_as)
    lok = (status == 0)
    IF (lok) THEN
       IF (TRIM(ti_gp(idt_so4mm_as)%tp%ident%submodel) == &
            TRIM(chemmodule)) THEN
          CALL info_bi('fetching chemistry module SO4mm_as tracer idt', substr)
#ifdef MESSYTENDENCY
    CALL mtend_register(my_handle, idt_so4mm_as)
#endif
       ELSE
          idt_so4mm_as =0 
          CALL info_bi('no MECCA-AERO SO4mm_as tracer initialized by'//& 
               'chosen chemistry module', substr)
       ENDIF
    ELSE
       idt_so4mm_as =0 
       CALL info_bi('no MECCA-AERO SO4mm_as tracer initialized by'//&
            'chosen chemistry module', substr)
    ENDIF

    ! couple to MECCA-AERO H2SO4_cs
    IF  (H2SO4_cs(1) == '' .AND. H2SO4_cs(2) == '') THEN
       CALL INFO_BI('no aerosol-phase H2SO4 tracer required', substr)
    ELSE
       CALL get_tracer(status, GPTRSTR, H2SO4_cs(1) &
            ,subname= H2SO4_cs(2), idx=idt_h2so4_cs)
       lok = (status == 0)
       IF (lok) THEN
          IF (TRIM(ti_gp(idt_h2so4_cs)%tp%ident%submodel) == &
               TRIM(chemmodule)) THEN
             CALL info_bi('fetching chemistry module H2SO4_cs tracer idt' &
                  , substr)
#ifdef MESSYTENDENCY
    CALL mtend_register(my_handle, idt_h2so4_cs)
#endif
          ELSE
             idt_h2so4_cs =0 
             CALL info_bi('no MECCA-AERO H2SO4_cs tracer initialized by'//&
                  'chosen chemistry module', substr)
          ENDIF
       ELSE   ! lok
          idt_h2so4_cs =0 
          CALL info_bi('no MECCA-AERO H2SO4_cs tracer initialized by'//&
               'chosen chemistry module', substr)
       ENDIF  !lok
    ENDIF
       
    ! find idts for the other MECCA-AERO CS tracers   
    IF (HSO4m_cs(1) == '' .AND. HSO4m_cs(2) == '') THEN
       CALL info_bi('no aerosol-phase HSO4m tracer required', substr)
    ELSE
       CALL get_tracer(status, GPTRSTR, HSO4m_cs(1), subname=HSO4m_cs(2)&
            , idx=idt_hso4m_cs)
       lok = (status == 0)
       IF (lok) THEN
          IF (TRIM(ti_gp(idt_hso4m_cs)%tp%ident%submodel) == &
               TRIM(chemmodule)) THEN
             CALL info_bi('fetching chemistry module HSO4m coarse mode tracer idt' &
                  , substr)
#ifdef MESSYTENDENCY
    CALL mtend_register(my_handle, idt_hso4m_cs)
#endif
          ELSE
             idt_hso4m_cs =0 
             CALL info_bi('no MECCA-AERO HSO4m tracer initialized by'//&
                  'chosen chemistry module', substr)
          ENDIF
       ELSE
          idt_hso4m_cs =0 
          CALL info_bi('no MECCA-AERO HSO4m_cs tracer initialized by'//&
               'chosen chemistry module', substr)
       ENDIF
    ENDIF
    !-----------------
    IF (SO4mm_cs(1) == '' .AND. SO4mm_cs(2) == '') THEN
       CALL INFO_BI('no aerosol-phase SO4mm tracer required', substr)
       idt_so4mm_cs =0 
    ELSE
       CALL get_tracer(status, GPTRSTR,SO4mm_cs(1),subname=SO4mm_cs(2)&
            , idx=idt_so4mm_cs)
       lok = (status == 0)
       IF (lok) THEN
          IF (TRIM(ti_gp(idt_so4mm_cs)%tp%ident%submodel) == &
               TRIM(chemmodule)) THEN
             CALL info_bi( &
                  'fetching chemistry module SO4mm coarse mode tracer idt' &
                  , substr)
#ifdef MESSYTENDENCY
!    CALL mtend_register(my_handle, mtend_id_tracer, idt=idt_so4mm_cs)
    CALL mtend_register(my_handle, idt_so4mm_cs)
#endif
          ELSE
             idt_so4mm_cs =0 
             CALL info_bi('no MECCA-AERO SO4mm_cs tracer initialized by'//&
                  'chosen chemistry module', substr)
          ENDIF
       ELSE !lok
          CALL info_bi('no MECCA-AERO SO4mm _cs tracer initialized by'//&
               'chosen chemistry module', substr)
       ENDIF
    ENDIF
    

    emis: IF (l_calc_emis) THEN
       if_ss: IF (l_SS) THEN
          ! get seasalt emission 
          CALL get_channel_object(status, &
               TRIM(SSemis_channel), TRIM(SS_mass_as), p2=Mss_as)
          IF (status /= 0) THEN
             CALL error_bi(&
                  'channel object for accumulation mode seasalt'//&
                  'mass emission flux not available', substr)               
          ELSE
             CALL info_bi('fetching  seasalt mass(as) emission flux', substr)
          END IF

          ! get accumulation mode number emission flux
          CALL get_channel_object(status, &
               TRIM(SSemis_channel), TRIM(SS_num_as), p2=Nss_as)
          IF (status /= 0) THEN
             CALL error_bi ('channel object for accumulation mode '//&
                  'seasalt number emission flux not available', substr)
          ELSE
             CALL info_bi('fetching  seasalt number(as) emission flux', substr)
          ENDIF
          
          ! get coarse mode mass emission flux
          CALL get_channel_object(status, &
               TRIM(SSemis_channel), TRIM(SS_mass_cs), p2=Mss_cs)
          IF (status /= 0) THEN
             CALL error_bi ('channel object for coarse mode seasalt'//&
                  'mass emission fluxnot available', substr)
          ELSE
             CALL info_bi('fetching  seasalt mass(cs) emission flux', substr)
          ENDIF

          ! get coarse mode number emission flux     
          CALL get_channel_object(status, &
               TRIM(SSemis_channel), TRIM(SS_num_cs), p2=Nss_cs)
          IF (status /= 0) THEN
             CALL error_bi ('channel object for coarse mode seasalt'//&
                  'number emission flux not available', substr)
          ELSE
             CALL info_bi('fetching  seasalt number(cs) emission flux', substr)
          ENDIF
       ENDIF if_ss

       if_carbon: IF(l_carbon) THEN

          ! get emission channel object for soluble OC emissions
          CALL get_channel_object(status, &
               TRIM(Cemis_channel), TRIM(emis_OC_sol), p2=OC_sum_sol)
          IF (status /= 0) THEN 
             CALL error_bi ('channel object for OC_sum_sol'//&
                  'emission flux not available', substr)
          ELSE
             CALL info_bi (' got channel object for OC_sum_sol', substr)
          ENDIF

          ! get emission channel object for insoluble OC emissions
          CALL get_channel_object(status, &
               TRIM(Cemis_channel), TRIM(emis_OC_insol), p2=OC_sum_insol)
          IF (status /= 0) THEN
               CALL error_bi (&
               'channel object for OC_sum_insol emission flux not available'&
               , substr) 
          ELSE
             CALL info_bi (' got channel object for OC_sum_insol', substr)
          ENDIF

          ! get emission channel object for insoluble BC emissions
          CALL get_channel_object(status, &
               TRIM(Cemis_channel), TRIM(emis_BC_insol), p2=BC_sum_insol)
          IF (status /= 0) THEN           
               CALL error_bi (&
               'channel object for BC_sum_insol emission flux not available'&
               , substr) 
          ELSE
             CALL info_bi (' got channel object for BC_sum_insol',substr)
          ENDIF

          ! get emission channel object for number soluble carbon emissions
          CALL get_channel_object(status, &
               TRIM(Cemis_channel), TRIM(emis_N_sol), p2=Nemis_ks)
          IF (status /= 0) THEN       
               CALL error_bi ('channel object for N_sol emission flux'//&
               ' not available', substr) 
          ELSE
             CALL info_bi (' got channel object for N_sol', substr)
          ENDIF

          ! get emission channel object for number insoluble carbon emissions
          CALL get_channel_object(status, &
               TRIM(Cemis_channel), TRIM(emis_N_insol), p2=Nemis_ki)
          IF (status /= 0) THEN
               CALL error_bi ('channel object for N_insol emission'//&
               'flux not available', substr) 
          ELSE
             CALL info_bi (' got channel object for N_insol', substr)
          ENDIF

          ! get 3D emission channel object BC wildfire emission
          CALL get_channel_object(status, 'offemis', 'aero1_emis_bc_wf', p3=BC_wf)
          IF (status /= 0) THEN
               CALL info_bi ('BC wildfire emission channel object'//&
               ' not available', substr) 
          ELSE
             CALL info_bi (substr,' got channel object for BC wildfire')
          ENDIF

          ! get 3D emission channel object OC wildfire emission
          CALL get_channel_object(status, 'offemis', 'aero2_emis_oc_wf', p3=OC_wf)
          IF (status /= 0) THEN         
               CALL info_bi ('OC wildfire emission channel object'//&
               ' not available', substr) 
          ELSE
             CALL info_bi (' got channel object for OC wildfire', substr)
          ENDIF

       ENDIF if_carbon

       IF (l_dust_ci) THEN ! dust emission required
          ! get emission channel object for coarse mode insoluble dust emissions
          CALL get_channel_object(status, &
               TRIM(Duemis_channel), TRIM(emis_dust_ci), p2=dust_emis_ci)
          IF (status /= 0) THEN               
             CALL error_bi ('channel object for ci dust '//&
                  'emission flux not available', substr)
          ELSE
             CALL info_bi (' got channel object for ci dust emissions', substr)
          ENDIF
       ENDIF
       IF (l_dust_ai) THEN
          ! get emission channel object for dust emissions
          CALL get_channel_object(status, &
               TRIM(Duemis_channel), TRIM(emis_dust_ai), p2=dust_emis_ai)
          IF (status /= 0) THEN               
             CALL error_bi ('channel object for ai dust '//&
                  'emission flux not available', substr)
          ELSE
             CALL info_bi (' got channel object for ai DUST emissions', substr)
          ENDIF
       ENDIF
          
       if_so4: IF (l_so4) THEN

          ! get emission channel object for lowest layer sufate emissions
          CALL get_channel_object(status, &
               TRIM(SO4emis_channel), TRIM(emis_so2_low), p2=so2_emis_low)
          IF (status /= 0) THEN
             CALL error_bi ('channel object for lowest layer sulfate'//&
                  'emission flux not available', substr)
          ELSE
             CALL info_bi (' got channel object for SO2 emission low', substr)
          ENDIF
         
          ! get emission channel object for 2nd lowest layer sufate emissions
          CALL get_channel_object(status, &
               TRIM(SO4emis_channel), TRIM(emis_so2_high), p2=so2_emis_high)
          IF (status /= 0) THEN
               CALL error_bi ('channel object for 2nd lowest'//&
               &' layer sulfate emission flux not available', substr)
          ELSE
             CALL info_bi (' got channel object for SO2 emission high', substr)
          ENDIF
          
       ENDIF if_so4

    ENDIF emis

    CALL end_message_bi(modstr,'CHECK NAMELIST SETTINGS', substr)
    
  END SUBROUTINE m7_init_coupling
  
!=============================================================================!

SUBROUTINE m7_vdiff

    ! ECHAM5/MESSy
    USE messy_main_grid_def_mem_bi, ONLY: nlev, jrow, kproma, nproma
    USE messy_main_data_bi,       ONLY: pressi_3d
#ifndef MESSYTENDENCY
    USE messy_main_tracer_mem_bi, ONLY: pxtte => qxtte
#endif
#if defined (ECHAM5)
    USE messy_main_data_bi,       ONLY: pxtems 
#endif
    ! MESSy
    USE messy_main_constants_mem, ONLY: M_air, g, pi
    USE messy_main_tracer,        ONLY: R_molarmass, R_aerosol_density

    INTRINSIC :: EXP, LOG

    ! LOCAL
    CHARACTER(len=*), PARAMETER :: substr='m7_vdiff'
    REAL(dp), POINTER :: zxtems(:,:)
    REAL(dp) :: zdp(nproma)
    REAL(dp) :: zdp2d(nproma, nlev)
    REAL(dp) :: fac

    !--- Dust source mass median diameter [m]:
    REAL(dp), PARAMETER :: srcmmd= 2.5e-6_dp
    !--- Dust source standard deviation:
    REAL(dp), PARAMETER :: srcsigma= 2._dp

    REAL(dp) :: dustconvMtoN_ci    = 0._dp
    REAL(dp) :: dustconvMtoN_ai    = 0._dp

    ! parameters for Tegen emission scheme
    !--- Dust source mass median diameter [m]:
    REAL(dp), PARAMETER :: srcmmd_ai= 0.7e-6_dp
    REAL(dp), PARAMETER :: srcmmd_ci= 3.5e-6_dp
    !--- Dust source standard deviation:
    REAL(dp), PARAMETER :: srcsigma_ai= 1.59_dp
    REAL(dp), PARAMETER :: srcsigma_ci= 2._dp
    !

    ! sulfate emission source: assumed median number radius
    REAL(dp), PARAMETER :: cmr_sk   = 0.03e-6_dp   ! aitken mode
    REAL(dp), PARAMETER :: cmr_sa   = 0.075e-6_dp  ! accumulation mode
    REAL(dp), PARAMETER :: cmr_sc   = 0.75e-6_dp   ! coarse mode 
    REAL(dp), PARAMETER :: cmr_bb   = 0.075E-6_dp  ! biomass burning

    ! (Seinfeld and Pandis, 1998, p709;  Ferek et al., JGR, 1998) 
    REAL(dp), PARAMETER :: zom2oc         = 1.4_dp
    ! Biom. Burn. Percentage of  Water Soluble OC (WSOC) [1]
    ! (M.O. Andreae; Talk: Smoke and Climate)
    REAL(dp), PARAMETER :: zbb_wsoc_perc  = 0.65_dp
    ! Assume same Percentage of WSOC for biogenic OC

    ! carbon number emission aerosol mode scaling factors
    REAL(dp) :: zm2n_C_sol
    REAL(dp) :: zm2n_C_insol
    REAL(dp) :: zm2n_bb
#ifdef MESSYTENDENCY
    REAL(dp) :: zxtte(nproma,nlev)
#endif
    INTEGER :: idt

#if defined (ECHAM5)
   zxtems => pxtems(:,1,:,jrow)
#endif
    IF (.not. l_calc_emis) RETURN

    ! determine the factor for number emissions
    ! the assumption that all number densities have the same unit is made
    fac = M_air * 1.E-3

    ! precalculate conversion factor for number emission in aitken mode
    ! attention: the assumption is made that OC and BC have the same density 
    IF (l_carbon) THEN
       zm2n_C_sol   = 1./( ti_gp(idt_mocks)%tp%meta%cask_r(R_aerosol_density) &
            * (cmr2ram(iaits))**3.)
       zm2n_C_insol = 1./(ti_gp(idt_mbcki)%tp%meta%cask_r(R_aerosol_density) &
            * (cmr2ram(iaiti))**3.)
       zm2n_bb= 3./(4.*pi*(cmr_bb)**3.)* &
            1./(ti_gp(idt_mbcki)%tp%meta%cask_r(R_aerosol_density) &
            * (cmr2ram(iaiti))**3.)
    ENDIF

    ! calculate mass to number conversion factor for dust emissions
    IF (l_dust_ci) THEN
       dustconvMtoN_ci =1./pi*6./&
            ti_gp(idt_mduci)%tp%meta%cask_r(R_aerosol_density)/(srcmmd_ci**3 &
            *EXP(4.5*LOG(srcsigma_ci)**2)) ! convert mass to number flux
    ENDIF
    IF (l_dust_ai) THEN
       dustconvMtoN_ai =1./pi*6./&
            ti_gp(idt_mduai)%tp%meta%cask_r(R_aerosol_density)/(srcmmd_ai**3 &
            *EXP(4.5*LOG(srcsigma_ai)**2)) ! convert mass to number flux
    ENDIF

#if defined (ECHAM5)
    tendency: If (.not. l_tendency) THEN
       
       l_ss1: IF (l_ss) THEN
          ! xtems for accumulation mass flux
          zxtems(1:kproma,idt_mssas) =  zxtems(1:kproma,idt_mssas) + &
               M_air / ti_gp(idt_mssas)%tp%meta%cask_r(R_molarmass)    * &
               Mss_as(1:kproma,jrow)
          
          ! xtems for coarse mass flux
          zxtems(1:kproma,idt_msscs) =  zxtems(1:kproma,idt_msscs) + &
               M_air / ti_gp(idt_msscs)%tp%meta%cask_r(R_molarmass)    * &
               Mss_cs(1:kproma,jrow) 
          
          ! xtems for accumulation number flux
          zxtems(1:kproma,idt_nas) =  zxtems(1:kproma,idt_nas) + &
               Nss_as(1:kproma, jrow) * fac
          ! xtems for coarse mass flux
          zxtems(1:kproma,idt_ncs) =  zxtems(1:kproma,idt_ncs) + &
               Nss_cs(1:kproma, jrow) * fac
     
       END IF l_ss1

       l_carbon1:IF (l_carbon) THEN
          ! OC soluble
          zxtems(1:kproma,idt_mocks) =  zxtems(1:kproma,idt_mocks) + &
               M_air / ti_gp(idt_mocks)%tp%meta%cask_r(R_molarmass) &
               * OC_sum_sol(1:kproma,jrow)
          ! BC soluble emission does not exist
          ! OC insoluble
          zxtems(1:kproma,idt_mocki) =  zxtems(1:kproma,idt_mocki) + &
               M_air / ti_gp(idt_mocki)%tp%meta%cask_r(R_molarmass) &
               * OC_sum_insol(1:kproma,jrow)
          ! BC insoluble
          zxtems(1:kproma,idt_mbcki) =  zxtems(1:kproma,idt_mbcki) + &
               M_air / ti_gp(idt_mbcki)%tp%meta%cask_r(R_molarmass)  &
               * BC_sum_insol(1:kproma,jrow)
          
          ! Number emissions aitken soluble
          zxtems(1:kproma,idt_nks) =  zxtems(1:kproma,idt_nks)  &
               + Nemis_ks(1:kproma,jrow) * zm2n_C_sol   * fac 
          ! Number emissions aitken insoluble
          zxtems(1:kproma,idt_nki) =  zxtems(1:kproma,idt_nki)  &
               + Nemis_ki(1:kproma,jrow) * zm2n_C_insol * fac 

       ENDIF l_carbon1

       ll_dust_ci:IF (l_dust_ci) THEN
          ! dust mass (coarse insoluble)
          zxtems(1:kproma,idt_mduci) =  zxtems(1:kproma,idt_mduci) + &
               M_air / ti_gp(idt_mduci)%tp%meta%cask_r(R_molarmass)  &
               * dust_emis_ci(1:kproma,jrow)
          ! Number emissions (coarse insoluble)
          zxtems(1:kproma,idt_nci) =  zxtems(1:kproma,idt_nci) +     &
                dust_emis_ci(1:kproma,jrow) * fac * dustconvMtoN_ci 
       ENDIF ll_dust_ci

       ! um_gg_20090923+
       ll_dust_ai:IF (l_dust_ai) THEN
          ! dust mass (aitken insoluble)
          zxtems(1:kproma,idt_mduai) =  zxtems(1:kproma,idt_mduai) + &
               M_air / ti_gp(idt_mduai)%tp%meta%cask_r(R_molarmass)  &
               * dust_emis_ai(1:kproma,jrow)
          ! Number emissions (aitken insoluble)
          zxtems(1:kproma,idt_nai) =  zxtems(1:kproma,idt_nai) +     &
                dust_emis_ai(1:kproma,jrow) * fac * dustconvMtoN_ai 
       ENDIF ll_dust_ai
       ! um_gg_20090923-

        l_so4_1:IF (l_so4) THEN
          ! I) 50 % of primary SO4 emissions go into aitken mode particles
          !    and 50 % in accumulation mode particles

          ! II) 2.5 % of SO2 emissions are assumed to be primary sulfate 

          ! so4 mass (aitken soluble)
          zxtems(1:kproma,idt_ms4ks) =  zxtems(1:kproma,idt_ms4ks) +  &
               M_air / ti_gp(idt_ms4ks)%tp%meta%cask_r(R_molarmass)       &
               * SO2_emis_low(1:kproma,jrow) * 0.5_dp * 0.025_dp
          ! Number emissions ( aitken soluble)
          zxtems(1:kproma,idt_nks) =  zxtems(1:kproma,idt_nks) +        & 
                SO2_emis_low(1:kproma,jrow)  *0.5_dp * 0.025_dp *fac    &
                * 3._dp/(4.*pi* ti_gp(idt_ms4ks)%tp%meta%cask_r(R_aerosol_density) &
                * (cmr_sk*cmr2ram(iaits))**3.)

          ! so4 mass (accumulation soluble)
          zxtems(1:kproma,idt_ms4as) =  zxtems(1:kproma,idt_ms4as) +  &
               M_air / ti_gp(idt_ms4as)%tp%meta%cask_r(R_molarmass)       &
               * SO2_emis_low(1:kproma,jrow) *0.5_dp *0.025_dp
          ! Number emissions ( accumulation soluble)
          zxtems(1:kproma,idt_nci) =  zxtems(1:kproma,idt_nci) +        &
               SO2_emis_low(1:kproma,jrow)  *0.5_dp *0.025_dp   *fac    &
                * 3._dp/(4.*pi* ti_gp(idt_ms4as)%tp%meta%cask_r(R_aerosol_density) &
                * (cmr_sa*cmr2ram(iaccs))**3.)
       ENDIF l_so4_1

    ELSE !tendency
#endif
       ! pressure difference per box
       zdp(1:kproma) = pressi_3d(_RI_XYZ__(1:kproma,jrow,nlev+1)) &
            - pressi_3d(_RI_XYZ__(1:kproma,jrow,nlev)) 

#ifdef MESSYTENDENCY
       zxtte(:,:) = 0._dp
#endif

       l_ss2:IF(l_ss) THEN
          ! tendency seasalt mass accumulation mode
#ifndef MESSYTENDENCY
          idt = idt_mssas
          pxtte(_RI_X_ZN_(1:kproma,nlev,idt)) =  &
               pxtte(_RI_X_ZN_(1:kproma,nlev,idt)) + &
               Mss_as(1:kproma,jrow) / zdp(1:kproma) * g /                   &
               ti_gp(idt_mssas)%tp%meta%cask_r(R_molarmass) * M_air 
#else
          zxtte(1:kproma,nlev) = Mss_as(1:kproma,jrow) / zdp(1:kproma) * g / &
               ti_gp(idt_mssas)%tp%meta%cask_r(R_molarmass) * M_air 
          CALL mtend_add_l(my_handle, idt_mssas, px=zxtte)
       
#endif          

          ! tendency seasalt mass coarse mode
#ifndef MESSYTENDENCY
          idt = idt_msscs
          pxtte(_RI_X_ZN_(1:kproma,nlev,idt)) =  &
               pxtte(_RI_X_ZN_(1:kproma,nlev,idt)) + &
               Mss_cs(1:kproma,jrow) / zdp(1:kproma) * g /                   &
               ti_gp(idt_msscs)%tp%meta%cask_r(R_molarmass) * M_air 
#else
          zxtte(1:kproma,nlev) = Mss_cs(1:kproma,jrow) / zdp(1:kproma) * g / &
               ti_gp(idt_msscs)%tp%meta%cask_r(R_molarmass) * M_air 
          CALL mtend_add_l(my_handle, idt_msscs, px=zxtte)
#endif          
          ! tendency number accumulation mode
#ifndef MESSYTENDENCY
          idt = idt_nas
          pxtte(_RI_X_ZN_(1:kproma,nlev,idt)) =  &
               pxtte(_RI_X_ZN_(1:kproma,nlev,idt)) + &
               Nss_as(1:kproma,jrow) *fac / zdp(1:kproma) * g         
#else
          zxtte(1:kproma,nlev) = Nss_as(1:kproma,jrow) *fac / zdp(1:kproma) * g
          CALL mtend_add_l(my_handle, idt_nas, px=zxtte)
#endif          
          ! tendency number coarse mode
#ifndef MESSYTENDENCY
          idt = idt_ncs
          pxtte(_RI_X_ZN_(1:kproma,nlev,idt)) =  &
               pxtte(_RI_X_ZN_(1:kproma,nlev,idt)) + &
               Nss_cs(1:kproma,jrow) *fac / zdp(1:kproma) * g  
#else
          zxtte(1:kproma,nlev) = Nss_cs(1:kproma,jrow) *fac / zdp(1:kproma) * g
          CALL mtend_add_l(my_handle, idt_ncs, px=zxtte)
#endif          
       ENDIF l_ss2

    ! calculate OC/BC emissions
    l_carbon2:IF (l_carbon) THEN

       ! aitken soluble OC mass
#ifndef MESSYTENDENCY
       idt=idt_mocks
       pxtte(_RI_X_ZN_(1:kproma,nlev,idt)) =  &
            pxtte(_RI_X_ZN_(1:kproma,nlev,idt)) +      &
            OC_sum_sol(1:kproma,jrow) / zdp(1:kproma) *g /         &
            ti_gp(idt_mocks)%tp%meta%cask_r(R_molarmass) * M_air 
#else
       zxtte(1:kproma,nlev) = OC_sum_sol(1:kproma,jrow) / zdp(1:kproma) *g / &
            ti_gp(idt_mocks)%tp%meta%cask_r(R_molarmass) * M_air 
          CALL mtend_add_l(my_handle, idt_mocks, px=zxtte)
#endif
       ! aitken insoluble OC mass 
#ifndef MESSYTENDENCY
       idt=idt_mocki
       pxtte(_RI_X_ZN_(1:kproma,nlev,idt)) =  &
            pxtte(_RI_X_ZN_(1:kproma,nlev,idt)) +    &
            OC_sum_insol(1:kproma,jrow) / zdp(1:kproma) *g /     &
            ti_gp(idt_mocki)%tp%meta%cask_r(R_molarmass) * M_air 
#else
       zxtte(1:kproma,nlev) = OC_sum_insol(1:kproma,jrow) / zdp(1:kproma) *g / &
            ti_gp(idt_mocki)%tp%meta%cask_r(R_molarmass) * M_air 
          CALL mtend_add_l(my_handle, idt_mocki, px=zxtte)
#endif

       ! aitken insoluble BC mass 
#ifndef MESSYTENDENCY
       idt = idt_mbcki
       pxtte(_RI_X_ZN_(1:kproma,nlev,idt)) =  &
            pxtte(_RI_X_ZN_(1:kproma,nlev,idt)) +     &
            BC_sum_insol(1:kproma,jrow) / zdp(1:kproma) *g /      &
            ti_gp(idt_mbcki)%tp%meta%cask_r(R_molarmass) * M_air 
#else
       zxtte(1:kproma,nlev) = BC_sum_insol(1:kproma,jrow) / zdp(1:kproma) *g / &
            ti_gp(idt_mocki)%tp%meta%cask_r(R_molarmass) * M_air 
          CALL mtend_add_l(my_handle, idt_mbcki, px=zxtte)
#endif

       ! aitken soluble number 
#ifndef MESSYTENDENCY
       idt=idt_nks
       pxtte(_RI_X_ZN_(1:kproma,nlev,idt)) =  &
            pxtte(_RI_X_ZN_(1:kproma,nlev,idt)) +     &
            Nemis_ks(1:kproma,jrow) / zdp(1:kproma) *g            &
            * fac * zm2n_C_sol 
#else
       zxtte(1:kproma,nlev) = Nemis_ks(1:kproma,jrow) / zdp(1:kproma) *g   &
            * fac * zm2n_C_sol 
       CALL mtend_add_l(my_handle, idt_nks, px=zxtte)
#endif

       ! aitken insoluble number 
#ifndef MESSYTENDENCY
       idt = idt_nki
       pxtte(_RI_X_ZN_(1:kproma,nlev,idt)) =  &
            pxtte(_RI_X_ZN_(1:kproma,nlev,idt)) + &
            Nemis_ki(1:kproma,jrow) / zdp(1:kproma) *g * fac * zm2n_C_insol 
#else
       zxtte(1:kproma,nlev) = Nemis_ki(1:kproma,jrow) / zdp(1:kproma) *g   &
            * fac * zm2n_C_sol 
       CALL mtend_add_l(my_handle, idt_nki, px=zxtte)
#endif
    ENDIF l_carbon2

    ! calculate dust emissions 
    l_dust2_ci: IF (l_dust_ci) THEN
#ifndef MESSYTENDENCY
       idt=idt_mduci
       pxtte(_RI_X_ZN_(1:kproma,nlev,idt)) =  &
            pxtte(_RI_X_ZN_(1:kproma,nlev,idt)) +     &
            dust_emis_ci(1:kproma,jrow) / zdp(1:kproma) *g /         &
            ti_gp(idt_mduci)%tp%meta%cask_r(R_molarmass) * M_air 
#else
       zxtte(1:kproma,nlev) = dust_emis_ci(1:kproma,jrow) / zdp(1:kproma) *g /  &
            ti_gp(idt_mduci)%tp%meta%cask_r(R_molarmass) * M_air 
       CALL mtend_add_l(my_handle, idt_mduci, px=zxtte)
#endif

#ifndef MESSYTENDENCY
       idt=idt_nci
       pxtte(_RI_X_ZN_(1:kproma,nlev,idt)) =  &
            pxtte(_RI_X_ZN_(1:kproma,nlev,idt)) +     &
            dust_emis_ci(1:kproma,jrow) / zdp(1:kproma) *g           &
            *fac * dustconvMtoN_ci 
#else
       zxtte(1:kproma,nlev) = dust_emis_ci(1:kproma,jrow) / zdp(1:kproma) *g  &
            *fac * dustconvMtoN_ci 
       CALL mtend_add_l(my_handle, idt_nci, px=zxtte)
#endif
    ENDIF l_dust2_ci

    l_dust2_ai: IF (l_dust_ai) THEN
#ifndef MESSYTENDENCY
       idt = idt_mduai
       pxtte(_RI_X_ZN_(1:kproma,nlev,idt)) =  &
            pxtte(_RI_X_ZN_(1:kproma,nlev,idt)) +     &
            dust_emis_ai(1:kproma,jrow) / zdp(1:kproma) *g /      &
            ti_gp(idt_mduai)%tp%meta%cask_r(R_molarmass) * M_air 
#else
       zxtte(1:kproma,nlev) = dust_emis_ai(1:kproma,jrow) / zdp(1:kproma) *g / &
            ti_gp(idt_mduai)%tp%meta%cask_r(R_molarmass) * M_air 
       CALL mtend_add_l(my_handle, idt_mduai, px=zxtte)
#endif

#ifndef MESSYTENDENCY
       idt=idt_nai
       pxtte(_RI_X_ZN_(1:kproma,nlev,idt)) =  &
            pxtte(_RI_X_ZN_(1:kproma,nlev,idt)) +     &
            dust_emis_ai(1:kproma,jrow) / zdp(1:kproma) *g        &
            *fac * dustconvMtoN_ai 
#else
       zxtte(1:kproma,nlev) = dust_emis_ai(1:kproma,jrow) / zdp(1:kproma) *g  &
            *fac * dustconvMtoN_ai 
       CALL mtend_add_l(my_handle, idt_nai, px=zxtte)
#endif
    ENDIF l_dust2_ai

    l_so4_2:IF (l_so4) THEN
       ! I) 50 % of primary SO4 emissions go into aitken mode particles
       !    and 50 % in accumulation mode particles
       
       ! II) 2.5 % of SO2 emissions are assumed to be primary sulfate 
       
       ! so4 mass (aitken soluble)
#ifndef MESSYTENDENCY
       idt = idt_ms4ks
       pxtte(_RI_X_ZN_(1:kproma,nlev,idt)) =  &
            pxtte(_RI_X_ZN_(1:kproma,nlev,idt)) +     &
            M_air / ti_gp(idt_ms4ks)%tp%meta%cask_r(R_molarmass)  &  
            / zdp(1:kproma) *g                                    &
            * so2_emis_low(1:kproma,jrow) * 0.5_dp * 0.025_dp
#else
       zxtte(1:kproma,nlev) = so2_emis_low(1:kproma,jrow) * 0.5_dp * 0.025_dp &
            * M_air / ti_gp(idt_ms4ks)%tp%meta%cask_r(R_molarmass)  &  
            / zdp(1:kproma) *g                                      
       CALL mtend_add_l(my_handle, idt_ms4ks, px=zxtte)
#endif

       ! Number emissions ( aitken soluble)
#ifndef MESSYTENDENCY
       idt = idt_nks
       pxtte(_RI_X_ZN_(1:kproma,nlev,idt)) =  &
            pxtte(_RI_X_ZN_(1:kproma,nlev,idt)) +      &
            g / zdp(1:kproma)                                       &
            * so2_emis_low(1:kproma,jrow) * 0.5_dp * 0.025_dp *fac    &
            * 3._dp/(4.*pi* ti_gp(idt_ms4ks)%tp%meta%cask_r(R_aerosol_density)  &
            * (cmr_sk*cmr2ram(iaits))**3.)
#else
       zxtte(1:kproma,nlev) = g / zdp(1:kproma)                       &
            * so2_emis_low(1:kproma,jrow) * 0.5_dp * 0.025_dp *fac    &
            * 3._dp/(4.*pi* ti_gp(idt_ms4ks)%tp%meta%cask_r(R_aerosol_density)  &
            * (cmr_sk*cmr2ram(iaits))**3.)
       CALL mtend_add_l(my_handle, idt_nks, px=zxtte)
#endif
           
       ! so4 mass (accumulation soluble)
#ifndef MESSYTENDENCY
       idt = idt_ms4as
       pxtte(_RI_X_ZN_(1:kproma,nlev,idt)) =  &
            pxtte(_RI_X_ZN_(1:kproma,nlev,idt)) +     &
            M_air / ti_gp(idt_ms4as)%tp%meta%cask_r(R_molarmass)  &  
            / zdp(1:kproma) *g                                    &
            * so2_emis_low(1:kproma,jrow) * 0.5_dp * 0.025_dp
#else
       zxtte(1:kproma,nlev) = &
            M_air / ti_gp(idt_ms4as)%tp%meta%cask_r(R_molarmass)  &  
            / zdp(1:kproma) *g                                    &
            * so2_emis_low(1:kproma,jrow) * 0.5_dp * 0.025_dp
       CALL mtend_add_l(my_handle, idt_ms4as, px=zxtte)
#endif
       ! Number emissions ( accumulation soluble)
#ifndef MESSYTENDENCY
       idt = idt_nas
       pxtte(_RI_X_ZN_(1:kproma,nlev,idt)) =  &
            pxtte(_RI_X_ZN_(1:kproma,nlev,idt)) +      &
              g  / zdp(1:kproma)                                     &
            * so2_emis_low(1:kproma,jrow) * 0.5_dp * 0.025_dp *fac    &
            * 3._dp/(4.*pi* ti_gp(idt_ms4as)%tp%meta%cask_r(R_aerosol_density) &
            * (cmr_sa*cmr2ram(iaccs))**3.)
#else
       zxtte(1:kproma,nlev) = g  / zdp(1:kproma)                      &
            * so2_emis_low(1:kproma,jrow) * 0.5_dp * 0.025_dp *fac    &
            * 3._dp/(4.*pi* ti_gp(idt_ms4as)%tp%meta%cask_r(R_aerosol_density) &
            * (cmr_sa*cmr2ram(iaccs))**3.)
       CALL mtend_add_l(my_handle, idt_nas, px=zxtte)
#endif
        ENDIF l_so4_2
#if defined (ECHAM5)
    ENDIF tendency
#endif

    ! pressure difference per box
    zdp2d(1:kproma, 1:nlev) = pressi_3d(_RI_XYZ__(1:kproma,jrow,2:nlev+1)) &
         - pressi_3d(_RI_XYZ__(1:kproma,jrow,1:nlev)) 
    
    ! 3D wildfire OC/BC emissions provided by AEROCOM
    l_carbon3: if (l_carbon) then
       IF (ASSOCIATED(BC_wf)) THEN
         ! mass
#ifndef MESSYTENDENCY
          idt=idt_mbcki
          pxtte(_RI_X_ZN_(1:kproma,1:nlev,idt)) =      &
               pxtte(_RI_X_ZN_(1:kproma,1:nlev,idt))   &
             + BC_wf(_RI_XYZ__(1:kproma,jrow,1:nlev))  &
             / zdp2d(1:kproma,1:nlev) *g   &
             / ti_gp(idt_mbcki)%tp%meta%cask_r(R_molarmass) * M_air 
          ! number
#else
          zxtte(1:kproma,1:nlev) = &
               BC_wf(_RI_XYZ__(1:kproma,jrow,1:nlev))  / zdp2d(1:kproma,1:nlev) *g   &
               / ti_gp(idt_mbcki)%tp%meta%cask_r(R_molarmass) * M_air 
          CALL mtend_add_l(my_handle, idt_mbcki, px=zxtte)
#endif
#ifndef MESSYTENDENCY
          idt = idt_nki
          pxtte(_RI_X_ZN_(1:kproma,1:nlev,idt)) =  &
               pxtte(_RI_X_ZN_(1:kproma,1:nlev,idt))           &
               + BC_wf(_RI_XYZ__(1:kproma,jrow,1:nlev))  /     &
               zdp2d(1:kproma,1:nlev) *g  &
               * zm2n_bb * fac
#else
          zxtte(1:kproma,1:nlev) = &
               + BC_wf(_RI_XYZ__(1:kproma,jrow,1:nlev))  / zdp2d(1:kproma,1:nlev) *g  &
               * zm2n_bb * fac
          CALL mtend_add_l(my_handle, idt_nki, px=zxtte)
#endif
        ENDIF
        
        IF (ASSOCIATED(OC_wf)) THEN
          ! aitken insoluble
          ! mass
#ifndef MESSYTENDENCY
          idt = idt_mocki
          pxtte(_RI_X_ZN_(1:kproma,1:nlev,idt)) = &
               pxtte(_RI_X_ZN_(1:kproma,1:nlev,idt))          &
             + OC_wf(_RI_XYZ__(1:kproma,jrow,1:nlev))  /      &
             zdp2d(1:kproma,1:nlev) *g   &
             / ti_gp(idt_mocki)%tp%meta%cask_r(R_molarmass) * M_air    &
             * zom2oc * (1._dp-zbb_wsoc_perc)
#else
          zxtte(1:kproma,1:nlev) = &
               OC_wf(_RI_XYZ__(1:kproma,jrow,1:nlev))  / zdp2d(1:kproma,1:nlev) *g   &
             / ti_gp(idt_mocki)%tp%meta%cask_r(R_molarmass) * M_air    &
             * zom2oc * (1._dp-zbb_wsoc_perc)          
          CALL mtend_add_l(my_handle, idt_mocki, px=zxtte)
#endif
           ! number
#ifndef MESSYTENDENCY
          idt = idt_nki
          pxtte(_RI_X_ZN_(1:kproma,1:nlev,idt)) =          &
               pxtte(_RI_X_ZN_(1:kproma,1:nlev,idt))       &
               + OC_wf(_RI_XYZ__(1:kproma,jrow,1:nlev))  / &
               zdp2d(1:kproma,1:nlev) *g                   &
               * zm2n_bb * fac                             &
               * zom2oc * (1._dp-zbb_wsoc_perc)
#else
          zxtte(1:kproma,1:nlev) =                       &
               OC_wf(_RI_XYZ__(1:kproma,jrow,1:nlev))  / &
               zdp2d(1:kproma,1:nlev) *g                 &
               * zm2n_bb * fac                           &
               * zom2oc * (1._dp-zbb_wsoc_perc)   
          CALL mtend_add_l(my_handle, idt_nki, px=zxtte)
#endif
          ! aitken soluble
          ! mass
#ifndef MESSYTENDENCY
          idt =idt_mocks
          pxtte(_RI_X_ZN_(1:kproma,1:nlev,idt)) =             &
               pxtte(_RI_X_ZN_(1:kproma,1:nlev,idt))          &
             + OC_wf(_RI_XYZ__(1:kproma,jrow,1:nlev))  /      &
             zdp2d(1:kproma,1:nlev) *g                        &
             / ti_gp(idt_mocki)%tp%meta%cask_r(R_molarmass) * M_air    &
             * zom2oc * zbb_wsoc_perc
#else
          zxtte(1:kproma,1:nlev) =                                     &
               OC_wf(_RI_XYZ__(1:kproma,jrow,1:nlev))  /               &
               zdp2d(1:kproma,1:nlev) *g                               &
             / ti_gp(idt_mocki)%tp%meta%cask_r(R_molarmass) * M_air    &
             * zom2oc * zbb_wsoc_perc        
          CALL mtend_add_l(my_handle, idt_mocks, px=zxtte)
#endif

           ! number
#ifndef MESSYTENDENCY
          idt = idt_nks
          pxtte(_RI_X_ZN_(1:kproma,1:nlev,idt)) =              &
               pxtte(_RI_X_ZN_(1:kproma,1:nlev,idt))           &
               + OC_wf(_RI_XYZ__(1:kproma,jrow,1:nlev))  /     &
               zdp2d(1:kproma,1:nlev) *g                       &
               * zm2n_bb * fac                                 &
               * zom2oc * zbb_wsoc_perc
#else
          zxtte(1:kproma,1:nlev) =                       &
               OC_wf(_RI_XYZ__(1:kproma,jrow,1:nlev))  / &
               zdp2d(1:kproma,1:nlev) *g                 &
               * zm2n_bb * fac                           &
               * zom2oc * zbb_wsoc_perc
          CALL mtend_add_l(my_handle, idt_nks, px=zxtte)
#endif
      ENDIF
    ENDIF l_carbon3
    ! sulfate emissions in second lowest layer
    ! only emission into tendency possible

    l_so4_3:IF (l_so4) THEN
       ! pressure difference per box
       zdp(1:kproma) = pressi_3d(_RI_XYZ__(1:kproma,jrow,nlev+1)) &
            - pressi_3d(_RI_XYZ__(1:kproma,jrow,nlev)) 

       ! I) 50 % of primary SO4 emissions go into accumulation  mode particles
       !    and 50 % in coarse mode particles
       
       ! II) 2.5 % of SO2 emissions are assumed to be primary sulfate 
       
       ! so4 mass (accumulation soluble)
#ifndef MESSYTENDENCY
       idt=idt_ms4as
       pxtte(_RI_X_ZN_(1:kproma,nlev-1,idt)) =                    &
            pxtte(_RI_X_ZN_(1:kproma,nlev-1,idt)) +               &
            M_air / ti_gp(idt_ms4as)%tp%meta%cask_r(R_molarmass)  &  
            / zdp(1:kproma) *g                                    &
            * so2_emis_high(1:kproma,jrow) * 0.5_dp * 0.025_dp
#else
       zxtte(:,:) = 0._dp
       zxtte(1:kproma,nlev-1) =                                   &
            M_air / ti_gp(idt_ms4as)%tp%meta%cask_r(R_molarmass)  &  
            / zdp(1:kproma) *g                                    &
            * so2_emis_high(1:kproma,jrow) * 0.5_dp * 0.025_dp
       CALL mtend_add_l(my_handle, idt_ms4as, px=zxtte)
#endif
       ! Number emissions ( accumulation soluble)
#ifndef MESSYTENDENCY
       idt = idt_nas
       pxtte(_RI_X_ZN_(1:kproma,nlev-1,idt)) =                        &
            pxtte(_RI_X_ZN_(1:kproma,nlev-1,idt))                     &
            + g / zdp(1:kproma)                                       &
            * so2_emis_high(1:kproma,jrow) * 0.5_dp * 0.025_dp *fac   &
            * 3._dp/(4.*pi* ti_gp(idt_ms4as)%tp%meta%cask_r(R_aerosol_density) &
            * (cmr_sa*cmr2ram(iaccs))**3.)
#else
       zxtte(:,:) = 0._dp
       zxtte(1:kproma,nlev-1) =                                       &
            g / zdp(1:kproma)                                         &
            * so2_emis_high(1:kproma,jrow) * 0.5_dp * 0.025_dp *fac   &
            * 3._dp/(4.*pi* ti_gp(idt_ms4as)%tp%meta%cask_r(R_aerosol_density) &
            * (cmr_sa*cmr2ram(iaccs))**3.)
       CALL mtend_add_l(my_handle, idt_nas, px=zxtte)
#endif

       ! so4 mass (coarse soluble)
#ifndef MESSYTENDENCY
       idt = idt_ms4cs
       pxtte(_RI_X_ZN_(1:kproma,nlev-1,idt )) =                    &
            pxtte(_RI_X_ZN_(1:kproma,nlev-1,idt)) +                &
            M_air / ti_gp(idt_ms4cs)%tp%meta%cask_r(R_molarmass)   &  
            / zdp(1:kproma) *g                                     &
            * so2_emis_high(1:kproma,jrow) * 0.5_dp * 0.025_dp
#else
       zxtte(:,:) = 0._dp
       zxtte(1:kproma,nlev-1) =                                    &
            M_air / ti_gp(idt_ms4cs)%tp%meta%cask_r(R_molarmass)   &  
            / zdp(1:kproma) *g                                     &
            * so2_emis_high(1:kproma,jrow) * 0.5_dp * 0.025_dp
       CALL mtend_add_l(my_handle, idt_ms4cs, px=zxtte)
#endif
       ! Number emissions ( coarse soluble)
#ifndef MESSYTENDENCY
       idt = idt_ncs
       pxtte(_RI_X_ZN_(1:kproma,nlev-1,idt)) =                             &
            pxtte(_RI_X_ZN_(1:kproma,nlev-1,idt))                          & 
            + g / zdp(1:kproma)                                            &
            * so2_emis_high(1:kproma,jrow) * 0.5_dp * 0.025_dp *fac        &
            * 3._dp/(4.*pi* ti_gp(idt_ms4cs)%tp%meta%cask_r(R_aerosol_density) &
            * (cmr_sc*cmr2ram(icoas))**3.)
#else
       zxtte(:,:) = 0._dp
       zxtte(1:kproma,nlev-1) =                                           &
            g / zdp(1:kproma)                                              &
            * so2_emis_high(1:kproma,jrow) * 0.5_dp * 0.025_dp *fac        &
            * 3._dp/(4.*pi* ti_gp(idt_ms4cs)%tp%meta%cask_r(R_aerosol_density) &
            * (cmr_sc*cmr2ram(icoas))**3.)
       CALL mtend_add_l(my_handle, idt_ncs, px=zxtte)
#endif

    ENDIF l_so4_3


END SUBROUTINE m7_vdiff

!=============================================================================!

SUBROUTINE m7_physc

#ifndef MESSYTENDENCY
  USE messy_main_tracer_mem_bi,  ONLY: pxtte => qxtte, pxtm1 => qxtm1
  USE messy_main_data_bi,        ONLY: tm1_3d, tte_3d, qm1_3d, qte_3d
#endif
  USE messy_main_grid_def_bi,    ONLY: gboxarea_2d
  USE messy_main_grid_def_mem_bi, ONLY: nlev, jrow, kproma, nproma
  USE messy_main_data_bi,        ONLY: press_3d, pressi_3d, rhum_3d
  USE messy_main_timer,          ONLY: time_step_len
  USE messy_main_blather_bi,     ONLY: error_bi
  USE messy_main_constants_mem,  ONLY: M_air, M_H2O, g, N_A, R_gas
  USE messy_main_tracer,         ONLY: R_molarmass

  !--- Parameter list:
  !
  !  zww = aerosol water content for each mode [kg(water) m-3(air)]

  ! LOCAL
  !CHARACTER(LEN=*), PARAMETER :: substr = 'm7_physc'
  REAL(dp),         PARAMETER :: vtmpc1 = M_air / M_H2O -1._dp
  REAL(dp) :: papp1(nproma,nlev), paphp1(nproma,nlev+1)
  INTEGER(i4) :: jl,jk,jc
  REAL(dp)     :: ztmst
  REAL(dp) :: zqtmst
  REAL(dp) :: zmass_pre(5),               zmass_post(5)
  REAL(dp) :: zgboxarea(nproma)
  REAL(dp) :: zgso4(nproma,nlev),      &
              zrelhum(nproma,nlev),        zdpg(nproma,nlev)
  REAL(dp) :: zaerml(nproma,nlev,naermod), zaernl(nproma,nlev,nmod), &
              zm6rp(nproma,nlev,nmod),     zm6dry(nproma,nlev,nsol), &
              zrhop(nproma,nlev,nmod),     zww(nproma,nlev,nmod)

  LOGICAL  :: loabort(5)=.FALSE.
  REAL(dp) :: cair(nproma,nlev), crhoair(nproma,nlev)
  REAl(dp) :: temp(nproma,nlev), sphum(nproma,nlev)
  REAl(dp) :: sulfmass_cs(nproma,nlev)
  REAl(dp) :: sulfmass_as(nproma,nlev)
#ifdef MESSYTENDENCY
  REAl(dp) :: startval(nproma,nlev)
  REAL(dp) :: zgso4_bef(nproma,nlev)
  REAL(dp) :: zaerml_bef(nproma,nlev,naermod), zaernl_bef(nproma,nlev,nmod)
  REAL(dp) :: zxtte(nproma,nlev)
#endif
  INTEGER :: idt, jm

  !--- 0. Initialisations: -----------------------------------------------

  ztmst  = time_step_len
  zqtmst = 1./time_step_len
  ! actual temperature
#ifndef MESSYTENDENCY
  temp(1:kproma,1:nlev) =   tm1_3d(_RI_XYZ__(1:kproma,jrow,1:nlev))+ &
                            tte_3d(_RI_XYZ__(1:kproma,jrow,1:nlev))*ztmst
#else
  CALL mtend_get_start_l(mtend_id_t, v0=temp)
#endif
  !actual specific humidity
#ifndef MESSYTENDENCY
  sphum(1:kproma, 1:nlev) = qm1_3d(_RI_XYZ__(1:kproma,jrow,1:nlev))+ &
                            qte_3d(_RI_XYZ__(1:kproma,jrow,1:nlev))*ztmst
#else
  CALL mtend_get_start_l(mtend_id_q, v0=sphum)
#endif

  crhoair(1:kproma,1:nlev) = &
       press_3d(_RI_XYZ__(1:kproma, jrow,1:nlev)) / &
       (R_gas*temp(1:kproma, 1:nlev)*(1.0_dp+vtmpc1*sphum(1:kproma, 1:nlev)))

  cair(1:kproma,1:nlev) = &
       (N_A/1.E6) * crhoair(1:kproma,1:nlev)

  zmass_pre(:)  = 0._dp
  zmass_post(:) = 0._dp


  !--- Calculate ambient properties: ----------------------------------

  papp1(1:kproma,1:nlev)   = press_3d(_RI_XYZ__(1:kproma,jrow,1:nlev))
  paphp1(1:kproma,1:nlev+1)= pressi_3d(_RI_XYZ__(1:kproma,jrow,1:nlev+1))
  DO jk=1,nlev
     DO  jl=1,kproma
        
        !--- Air mass auxiliary variable:
        
        zdpg(jl,jk)=(paphp1(jl,jk+1)-paphp1(jl,jk))/g
        
     ENDDO
  ENDDO
  zrelhum (1:kproma,1:nlev)   = &
       MAX(0._dp, MIN(0.95_dp, rhum_3d(_RI_XYZ__(1:kproma,jrow,1:nlev)) / 100._dp))
  
  !--- 3) Convert units and add the tendencies of preceeding proceses: ----
  !       Transform:
  !         - mass of sulfur species from [mass(S)/mass(air)]
  !           into [molecules/cm+3]
  !         - mass of all other species into micro-gram/cubic-meter
  !         - particle numbers are transformed from [N/kg(air)]
  !           into [N/cm+3]
  !
  !--- 3.1 Gases:
#ifndef MESSYTENDENCY

  DO jk=1,nlev
     DO jl=1,kproma

        !--- 3.1 Gases:
        idt = idt_so4_m7

        zgso4(jl,jk)=cair(jl,jk)* &
             (pxtm1(_RI_X_ZN_(jl,jk,idt))+ pxtte(_RI_X_ZN_(jl,jk,idt))*ztmst)

        !--- Discard negative values:

        zgso4(jl,jk)= MAX(zgso4(jl,jk),0._dp)

        !--- 3.2 Particle mass:

        !--- 3.2.1) Sulfate mass:

        idt = idt_ms4ns
        zaerml(jl,jk,iso4ns)=cair(jl,jk)*(pxtm1(_RI_X_ZN_(jl,jk,idt))+ &
                                 pxtte(_RI_X_ZN_(jl,jk,idt))*ztmst)
        idt = idt_ms4ks
        zaerml(jl,jk,iso4ks)=cair(jl,jk)*(pxtm1(_RI_X_ZN_(jl,jk,idt))+ &
                                 pxtte(_RI_X_ZN_(jl,jk,idt))*ztmst)
        IF ((idt_h2so4_as/=0) .and. (idt_hso4m_as/=0).and. &
             (idt_so4mm_as/=0)) THEN
           idt = idt_h2so4_as
           sulfmass_as(jl,jk) = &
                pxtm1(_RI_X_ZN_(jl,jk,idt)) + pxtte(_RI_X_ZN_(jl,jk,idt))*ztmst
           idt = idt_hso4m_as
           sulfmass_as(jl,jk) = sulfmass_as(jl,jk) + &
                pxtm1(_RI_X_ZN_(jl,jk,idt)) + pxtte(_RI_X_ZN_(jl,jk,idt))*ztmst
           idt = idt_so4mm_as
           sulfmass_as(jl,jk) = sulfmass_as(jl,jk) + &
                pxtm1(_RI_X_ZN_(jl,jk,idt))+ pxtte(_RI_X_ZN_(jl,jk,idt))*ztmst

           zaerml(jl,jk,iso4as)=cair(jl,jk)*sulfmass_as(jl,jk)
        ELSE
           idt = idt_ms4as
           zaerml(jl,jk,iso4as)=cair(jl,jk)*(pxtm1(_RI_X_ZN_(jl,jk,idt))+ &
                                 pxtte(_RI_X_ZN_(jl,jk,idt))*ztmst)
        ENDIF

        ! IMPROVED FOR MECCA-AERO COUPLING
        IF ((idt_h2so4_cs/=0) .or. (idt_hso4m_cs/=0).or. &
             (idt_so4mm_cs/=0)) THEN
           sulfmass_cs(jl,jk) = 0

           idt = idt_h2so4_cs
           IF (idt_h2so4_cs/=0) &
               sulfmass_cs(jl,jk) = pxtm1(_RI_X_ZN_(jl,jk,idt)) &
                + pxtte(_RI_X_ZN_(jl,jk,idt)) *ztmst

           idt = idt_hso4m_cs
           IF (idt_hso4m_cs/=0) &
                sulfmass_cs(jl,jk) = sulfmass_cs(jl,jk) &
                + pxtm1(_RI_X_ZN_(jl,jk,idt)) + pxtte(_RI_X_ZN_(jl,jk,idt)) *ztmst

           idt = idt_so4mm_cs
           IF (idt_so4mm_cs/=0) &
                sulfmass_cs(jl,jk) = sulfmass_cs(jl,jk) &
                + pxtm1(_RI_X_ZN_(jl,jk,idt)) + pxtte(_RI_X_ZN_(jl,jk,idt)) *ztmst
           zaerml(jl,jk,iso4cs)=cair(jl,jk)*sulfmass_cs(jl,jk)
        ELSE
           idt = idt_ms4cs
           zaerml(jl,jk,iso4cs)=cair(jl,jk)*(pxtm1(_RI_X_ZN_(jl,jk,idt))+ &
                pxtte(_RI_X_ZN_(jl,jk,idt))*ztmst)
        ENDIF
        !--- 3.2.3) Black Carbon:
        idt = idt_mbcks
        zaerml(jl,jk,ibcks)=crhoair(jl,jk)*             &
             ti_gp(idt_mbcks)%tp%meta%cask_r(R_molarmass) * &
             1.E6*(pxtm1(_RI_X_ZN_(jl,jk,idt))+ pxtte(_RI_X_ZN_(jl,jk,idt))*ztmst)

        idt = idt_mbcas
        zaerml(jl,jk,ibcas)=crhoair(jl,jk)*             &
             ti_gp(idt_mbcas)%tp%meta%cask_r(R_molarmass) * &
             1.E6*(pxtm1(_RI_X_ZN_(jl,jk,idt))+ pxtte(_RI_X_ZN_(jl,jk,idt))*ztmst)

        idt = idt_mbccs
        zaerml(jl,jk,ibccs)=crhoair(jl,jk)*             &
             ti_gp(idt_mbccs)%tp%meta%cask_r(R_molarmass) * &
             1.E6*(pxtm1(_RI_X_ZN_(jl,jk,idt))+ pxtte(_RI_X_ZN_(jl,jk,idt))*ztmst)

        idt = idt_mbcki
        zaerml(jl,jk,ibcki)=crhoair(jl,jk)*             &
             ti_gp(idt_mbcki)%tp%meta%cask_r(R_molarmass) * &
             1.E6*(pxtm1(_RI_X_ZN_(jl,jk,idt))+ pxtte(_RI_X_ZN_(jl,jk,idt))*ztmst)

        !--- 3.2.3) Organic Carbon:
        
        idt = idt_mocks
        zaerml(jl,jk,iocks)=crhoair(jl,jk)*             &
             ti_gp(idt_mocks)%tp%meta%cask_r(R_molarmass) * &
             1.E6*(pxtm1(_RI_X_ZN_(jl,jk,idt))+ pxtte(_RI_X_ZN_(jl,jk,idt))*ztmst)

        idt = idt_mocas
        zaerml(jl,jk,iocas)=crhoair(jl,jk)*             &
             ti_gp(idt_mocas)%tp%meta%cask_r(R_molarmass) * &
             1.E6*(pxtm1(_RI_X_ZN_(jl,jk,idt))+ pxtte(_RI_X_ZN_(jl,jk,idt))*ztmst)

        idt = idt_moccs
        zaerml(jl,jk,ioccs)=crhoair(jl,jk)*             &
             ti_gp(idt_moccs)%tp%meta%cask_r(R_molarmass) * &
             1.E6*(pxtm1(_RI_X_ZN_(jl,jk,idt))+ pxtte(_RI_X_ZN_(jl,jk,idt))*ztmst)

        idt = idt_mocki
        zaerml(jl,jk,iocki)=crhoair(jl,jk)*             &
             ti_gp(idt_mocki)%tp%meta%cask_r(R_molarmass) * &
             1.E6*(pxtm1(_RI_X_ZN_(jl,jk,idt))+ pxtte(_RI_X_ZN_(jl,jk,idt))*ztmst)

        !--- 3.2.4) Sea Salt:

        idt = idt_mssas
        zaerml(jl,jk,issas)=crhoair(jl,jk)*             &
             ti_gp(idt_mssas)%tp%meta%cask_r(R_molarmass) * &
             1.E6*(pxtm1(_RI_X_ZN_(jl,jk,idt))+ pxtte(_RI_X_ZN_(jl,jk,idt))*ztmst)

        idt = idt_msscs
        zaerml(jl,jk,isscs)=crhoair(jl,jk)*             &
             ti_gp(idt_msscs)%tp%meta%cask_r(R_molarmass) * &
             1.E6*(pxtm1(_RI_X_ZN_(jl,jk,idt))+ pxtte(_RI_X_ZN_(jl,jk,idt))*ztmst)
 
        !--- 3.2.5) Dust:
        idt =idt_mduas
        zaerml(jl,jk,iduas)=crhoair(jl,jk)*             &
             ti_gp(idt_mduas)%tp%meta%cask_r(R_molarmass) * &
             1.E6*(pxtm1(_RI_X_ZN_(jl,jk,idt))+ pxtte(_RI_X_ZN_(jl,jk,idt))*ztmst)
        
        idt =idt_mducs
        zaerml(jl,jk,iducs)=crhoair(jl,jk)*             &
             ti_gp(idt_mducs)%tp%meta%cask_r(R_molarmass) * &
             1.E6*(pxtm1(_RI_X_ZN_(jl,jk,idt))+  pxtte(_RI_X_ZN_(jl,jk,idt))*ztmst)
        
        idt =idt_mduai
        zaerml(jl,jk,iduai)=crhoair(jl,jk)*             &
             ti_gp(idt_mduai)%tp%meta%cask_r(R_molarmass) * &
             1.E6*(pxtm1(_RI_X_ZN_(jl,jk,idt))+  pxtte(_RI_X_ZN_(jl,jk,idt))*ztmst)
        
        idt =idt_mduci
        zaerml(jl,jk,iduci)=crhoair(jl,jk)*             &
             ti_gp(idt_mduci)%tp%meta%cask_r(R_molarmass) * &
             1.E6*(pxtm1(_RI_X_ZN_(jl,jk,idt))+ pxtte(_RI_X_ZN_(jl,jk,idt))*ztmst)

        !--- Discard negative values:
        zaerml(jl,jk,:)=MAX(zaerml(jl,jk,:),0._dp)

        !--- 3.3) Particle numbers:
        SELECT CASE (TRIM(ti_gp(idt_ncs)%tp%ident%unit))
           CASE('mol(part)/mol')
              idt = idt_nns
              zaernl(jl,jk,inucs)= crhoair(jl,jk)/1.e6 &
                   * (pxtm1(_RI_X_ZN_(jl,jk,idt))+ pxtte(_RI_X_ZN_(jl,jk,idt))*ztmst) * N_A
              idt = idt_nks
              zaernl(jl,jk,iaits)=crhoair(jl,jk)/1.e6  &
                   * (pxtm1(_RI_X_ZN_(jl,jk,idt))+ pxtte(_RI_X_ZN_(jl,jk,idt))*ztmst) * N_A
              idt = idt_nas
              zaernl(jl,jk,iaccs)=crhoair(jl,jk)/1.e6 &
                   * (pxtm1(_RI_X_ZN_(jl,jk,idt))+ pxtte(_RI_X_ZN_(jl,jk,idt))*ztmst) * N_A
              idt = idt_ncs
              zaernl(jl,jk,icoas)=crhoair(jl,jk)/1.e6 &
                   * (pxtm1(_RI_X_ZN_(jl,jk,idt))+ pxtte(_RI_X_ZN_(jl,jk,idt))*ztmst) * N_A
              idt = idt_nki
              zaernl(jl,jk,iaiti)=crhoair(jl,jk)/1.e6 &
                   * (pxtm1(_RI_X_ZN_(jl,jk,idt))+ pxtte(_RI_X_ZN_(jl,jk,idt))*ztmst) * N_A
              idt = idt_nai
              zaernl(jl,jk,iacci)=crhoair(jl,jk)/1.e6 &
                   * (pxtm1(_RI_X_ZN_(jl,jk,idt))+ pxtte(_RI_X_ZN_(jl,jk,idt))*ztmst) * N_A
              idt = idt_nci
              zaernl(jl,jk,icoai)=crhoair(jl,jk)/1.e6 &
                   * (pxtm1(_RI_X_ZN_(jl,jk,idt))+ pxtte(_RI_X_ZN_(jl,jk,idt))*ztmst) * N_A
           CASE('1/mol')
              idt =idt_nns
              zaernl(jl,jk,inucs)= crhoair(jl,jk)/1.e6 &
                   * (pxtm1(_RI_X_ZN_(jl,jk,idt))+ pxtte(_RI_X_ZN_(jl,jk,idt))*ztmst)
              idt =idt_nks
              zaernl(jl,jk,iaits)=crhoair(jl,jk)/1.e6 &
                   * (pxtm1(_RI_X_ZN_(jl,jk,idt))+ pxtte(_RI_X_ZN_(jl,jk,idt))*ztmst)
              idt =idt_nas
              zaernl(jl,jk,iaccs)=crhoair(jl,jk)/1.e6 &
                   * (pxtm1(_RI_X_ZN_(jl,jk,idt))+ pxtte(_RI_X_ZN_(jl,jk,idt))*ztmst)
              idt =idt_ncs
              zaernl(jl,jk,icoas)=crhoair(jl,jk)/1.e6 &
                   * (pxtm1(_RI_X_ZN_(jl,jk,idt))+ pxtte(_RI_X_ZN_(jl,jk,idt))*ztmst)
              idt =idt_nki
              zaernl(jl,jk,iaiti)=crhoair(jl,jk)/1.e6 &
                   * (pxtm1(_RI_X_ZN_(jl,jk,idt))+ pxtte(_RI_X_ZN_(jl,jk,idt))*ztmst)
              idt =idt_nai
              zaernl(jl,jk,iacci)=crhoair(jl,jk)/1.e6 &
                   * (pxtm1(_RI_X_ZN_(jl,jk,idt))+ pxtte(_RI_X_ZN_(jl,jk,idt))*ztmst)
              idt =idt_nci
              zaernl(jl,jk,icoai)=crhoair(jl,jk)/1.e6 &
                   * (pxtm1(_RI_X_ZN_(jl,jk,idt))+ pxtte(_RI_X_ZN_(jl,jk,idt))*ztmst)
           CASE('1/kg')
              idt =idt_nns
              zaernl(jl,jk,inucs)= crhoair(jl,jk)*M_air/1.e9* &
                    (pxtm1(_RI_X_ZN_(jl,jk,idt))+ pxtte(_RI_X_ZN_(jl,jk,idt))*ztmst)
              idt =idt_nks
              zaernl(jl,jk,iaits)=crhoair(jl,jk)*M_air/1.e9*  &
                    (pxtm1(_RI_X_ZN_(jl,jk,idt))+ pxtte(_RI_X_ZN_(jl,jk,idt))*ztmst)
              idt =idt_nas
              zaernl(jl,jk,iaccs)=crhoair(jl,jk)*M_air/1.e9* &
                    (pxtm1(_RI_X_ZN_(jl,jk,idt))+ pxtte(_RI_X_ZN_(jl,jk,idt))*ztmst)
              idt =idt_ncs
              zaernl(jl,jk,icoas)=crhoair(jl,jk)*M_air/1.e9* &
                    (pxtm1(_RI_X_ZN_(jl,jk,idt))+ pxtte(_RI_X_ZN_(jl,jk,idt))*ztmst)
              idt =idt_nki
              zaernl(jl,jk,iaiti)=crhoair(jl,jk)*M_air/1.e9* &
                    (pxtm1(_RI_X_ZN_(jl,jk,idt))+ pxtte(_RI_X_ZN_(jl,jk,idt))*ztmst)
              idt =idt_nai
              zaernl(jl,jk,iacci)=crhoair(jl,jk)*M_air/1.e9* &
                    (pxtm1(_RI_X_ZN_(jl,jk,idt))+ pxtte(_RI_X_ZN_(jl,jk,idt))*ztmst)
              idt =idt_nci
              zaernl(jl,jk,icoai)=crhoair(jl,jk)*M_air/1.e9* &
                    (pxtm1(_RI_X_ZN_(jl,jk,idt))+ pxtte(_RI_X_ZN_(jl,jk,idt))*ztmst)
           END SELECT
                
        !--- Discard negative values:
        zaernl(jl,jk,:)=MAX(zaernl(jl,jk,:),0._dp)

     ENDDO
  ENDDO

#else
  ! TRACER BUDGET
        !--- 3.1 Gases:

        CALL mtend_get_start_l(idt_so4_m7, v0=zgso4)
        zgso4(:,:) = MAX(zgso4(:,:),0._dp)
        zgso4_bef(:,:) = zgso4(:,:)
        !--- 3.2 Particle mass:

        !--- 3.2.1) Sulfate mass:

        CALL mtend_get_start_l(idt_ms4ns, v0=zaerml(:,:,iso4ns))        
        CALL mtend_get_start_l(idt_ms4ks, v0=zaerml(:,:,iso4ks))

        IF ((idt_h2so4_as/=0) .and. (idt_hso4m_as/=0).and. &
             (idt_so4mm_as/=0)) THEN
           CALL mtend_get_start_l(idt_h2so4_as, v0=sulfmass_as(:,:))

          CALL mtend_get_start_l(idt_hso4m_as, v0=startval)
          sulfmass_as(:,:) = sulfmass_as(:,:) + startval

          CALL mtend_get_start_l(idt_so4mm_as, v0=startval)
          sulfmass_as(:,:) = sulfmass_as(:,:) + startval

          zaerml(:,:,iso4as)=cair(:,:)*sulfmass_as(:,:)

       ELSE
          CALL mtend_get_start_l(idt_ms4as, v0=zaerml(:,:,iso4as))
        ENDIF

        ! IMPROVED FOR MECCA-AERO COUPLING
        IF ((idt_h2so4_cs/=0) .or. (idt_hso4m_cs/=0).or. &
             (idt_so4mm_cs/=0)) THEN

           IF (idt_h2so4_cs/=0) THEN
              CALL mtend_get_start_l(idt_h2so4_cs, v0=sulfmass_cs(:,:))
           ENDIF

           IF (idt_hso4m_cs/=0) THEN
              CALL mtend_get_start_l(idt_hso4m_cs, v0=startval)
              sulfmass_cs(:,:) = sulfmass_cs(:,:) + startval
           ENDIF
           IF (idt_so4mm_cs/=0) THEN
              CALL mtend_get_start_l(idt_so4mm_cs, v0=startval)
              sulfmass_cs(:,:) = sulfmass_cs(:,:) + startval
           ENDIF
           zaerml(:,:,iso4cs)=cair(:,:)*sulfmass_cs(:,:)
        ELSE
              CALL mtend_get_start_l(idt_ms4cs, v0=zaerml(:,:,iso4cs))
        ENDIF

        !--- 3.2.3) Black Carbon:

        CALL mtend_get_start_l(idt_mbcks, v0=zaerml(:,:,ibcks))
        CALL mtend_get_start_l(idt_mbcas, v0=zaerml(:,:,ibcas))
        CALL mtend_get_start_l(idt_mbccs, v0=zaerml(:,:,ibccs))
        CALL mtend_get_start_l(idt_mbcki, v0=zaerml(:,:,ibcki))

        !--- 3.2.3) Organic Carbon:
        
        CALL mtend_get_start_l(idt_mocks, v0=zaerml(:,:,iocks))
        CALL mtend_get_start_l(idt_mocas, v0=zaerml(:,:,iocas))
        CALL mtend_get_start_l(idt_moccs, v0=zaerml(:,:,ioccs))
        CALL mtend_get_start_l(idt_mocki, v0=zaerml(:,:,iocki))

        !--- 3.2.4) Sea Salt:

        CALL mtend_get_start_l(idt_mssas, v0=zaerml(:,:,issas))
        CALL mtend_get_start_l(idt_msscs, v0=zaerml(:,:,isscs))
 
        !--- 3.2.5) Dust:
        
        CALL mtend_get_start_l(idt_mduas, v0=zaerml(:,:,iduas))
        CALL mtend_get_start_l(idt_mducs, v0=zaerml(:,:,iducs))
        CALL mtend_get_start_l(idt_mduai, v0=zaerml(:,:,iduai))
        CALL mtend_get_start_l(idt_mduci, v0=zaerml(:,:,iduci))

        !--- Discard negative values:
         
        zaerml(:,:,:)=MAX(zaerml(:,:,:),0._dp)
        zaerml_bef(1:kproma,:,:) = zaerml(1:kproma,:,:)

        !--- 3.3) Particle numbers:
        
        SELECT CASE (TRIM(ti_gp(idt_ncs)%tp%ident%unit))
           CASE('mol(part)/mol')
              CALL mtend_get_start_l(idt_nns, v0=startval)
              zaernl(:,:,inucs) = crhoair(:,:)/1.e6 * startval(:,:) * N_A
 
              CALL mtend_get_start_l(idt_nks, v0=startval)
              zaernl(:,:,iaits) = crhoair(:,:)/1.e6 * startval(:,:) * N_A
 
              CALL mtend_get_start_l(idt_nas, v0=startval)
              zaernl(:,:,iaccs) = crhoair(:,:)/1.e6 * startval(:,:) * N_A
 
              CALL mtend_get_start_l(idt_ncs, v0=startval )
              zaernl(:,:,icoas) = crhoair(:,:)/1.e6 * startval(:,:) * N_A
 
              CALL mtend_get_start_l(idt_nki, v0=startval )
              zaernl(:,:,iaiti) = crhoair(:,:)/1.e6 * startval(:,:) * N_A
 
              CALL mtend_get_start_l(idt_nai, v0=startval )
              zaernl(:,:,iacci) = crhoair(:,:)/1.e6 * startval(:,:) * N_A
 
              CALL mtend_get_start_l(idt_nci, v0=startval )
              zaernl(:,:,icoai) = crhoair(:,:)/1.e6 * startval(:,:) * N_A
 
           CASE('1/mol')

              CALL mtend_get_start_l(idt_nns, v0=startval )
              zaernl(:,:,inucs) = crhoair(:,:)/1.e6 * startval(:,:) 
 
              CALL mtend_get_start_l(idt_nks, v0=startval )
              zaernl(:,:,iaits) = crhoair(:,:)/1.e6 * startval(:,:) 
 
              CALL mtend_get_start_l(idt_nas, v0=startval)
              zaernl(:,:,iaccs) = crhoair(:,:)/1.e6 * startval(:,:) 
 
              CALL mtend_get_start_l(idt_ncs, v0=startval)
              zaernl(:,:,icoas) = crhoair(:,:)/1.e6 * startval(:,:) 
 
              CALL mtend_get_start_l(idt_nki, v0=startval)
              zaernl(:,:,iaiti) = crhoair(:,:)/1.e6 * startval(:,:) 

              CALL mtend_get_start_l(idt_nai, v0=startval)
              zaernl(:,:,iacci) = crhoair(:,:)/1.e6 * startval(:,:) 
 
              CALL mtend_get_start_l(idt_nci, v0=startval)
              zaernl(:,:,icoai) = crhoair(:,:)/1.e6 * startval(:,:) 
 
           CASE('1/kg')

              CALL mtend_get_start_l(idt_nns, v0=startval)
              zaernl(:,:,inucs) = crhoair(:,:)*M_air/1.e9 * startval(:,:) 
 
              CALL mtend_get_start_l(idt_nks, v0=startval)
              zaernl(:,:,iaits) = crhoair(:,:)*M_air/1.e9 * startval(:,:) 

              CALL mtend_get_start_l(idt_nas, v0=startval)
              zaernl(:,:,iaccs) = crhoair(:,:)*M_air/1.e9 * startval(:,:) 

              CALL mtend_get_start_l(idt_ncs, v0=startval)
              zaernl(:,:,icoas) = crhoair(:,:)*M_air/1.e9 * startval(:,:) 

              CALL mtend_get_start_l(idt_nki, v0=startval)
              zaernl(:,:,iaiti) = crhoair(:,:)*M_air/1.e9 * startval(:,:) 

              CALL mtend_get_start_l(idt_nai, v0=startval)
              zaernl(:,:,iacci) = crhoair(:,:)*M_air/1.e9 * startval(:,:) 

              CALL mtend_get_start_l(idt_nci, v0=startval)
              zaernl(:,:,icoai) = crhoair(:,:)*M_air/1.e9 * startval(:,:) 

           END SELECT
                
        !--- Discard negative values:

        zaernl(:,:,:)=MAX(zaernl(:,:,:),0._dp)
        zaernl_bef(1:kproma,:,:) = zaernl(1:kproma,:,:)
#endif

  !--- Sum total mass of all compounds for mass diagnostics:
  IF(lmass_diag) CALL sum_mass(zmass_pre)


  !--- 4) Call of m7: -----------------------------------------------------

  If (lm7) &
  CALL m7_main(kproma,nproma,   nlev, ztmst,     &  ! ECHAM indices
          papp1, zrelhum, temp,           &  !   "   thermodynamics
          zgso4,          zaerml, zaernl, &  !  M7   tracers
          zm6rp, zm6dry,  zrhop,  zww      ) !   "   aerosol properties 

  !--- 5) Reconvert masses and numbers into mixing ratios, other ----------
  !       quantities to SI units and calculate the tendencies (xte):
  !
#ifndef MESSYTENDENCY
  DO jk=1,nlev
     DO jl=1,kproma

        !--- 5.1 Gases:
        idt =idt_so4_m7
        pxtte(_RI_X_ZN_(jl,jk,idt))=(zgso4(jl,jk)/cair(jl,jk)-&
                                              pxtm1(_RI_X_ZN_(jl,jk,idt)))*zqtmst

        !--- 5.2) Particle mass:
        
        !--- 5.2.1) Sulfate mass:
        
        idt =idt_ms4ns
        pxtte(_RI_X_ZN_(jl,jk,idt))=(zaerml(jl,jk,iso4ns)/cair(jl,jk)-&
             pxtm1(_RI_X_ZN_(jl,jk,idt)))*zqtmst
        idt =idt_ms4ks
        pxtte(_RI_X_ZN_(jl,jk,idt))=(zaerml(jl,jk,iso4ks)/cair(jl,jk)-&
             pxtm1(_RI_X_ZN_(jl,jk,idt)))*zqtmst

        ! included for MECCA-AERO coupling
        IF (idt_h2so4_as /= 0 .AND. idt_hso4m_as/= 0 .AND. &
             idt_so4mm_as /=0) THEN
           idt =idt_h2so4_as
           pxtte(_RI_X_ZN_(jl,jk,idt))=  pxtte(_RI_X_ZN_(jl,jk,idt))      &
                +(zaerml(jl,jk,iso4as)/cair(jl,jk)- &
                sulfmass_as(jl,jk))*zqtmst
        ENDIF
        idt =idt_ms4as
        pxtte(_RI_X_ZN_(jl,jk,idt))=(zaerml(jl,jk,iso4as)/cair(jl,jk)-&
                                   pxtm1(_RI_X_ZN_(jl,jk,idt)))*zqtmst

        ! included for MECCA-AERO coupling
        IF (idt_h2so4_cs /= 0) THEN
             idt = idt_h2so4_cs
             pxtte(_RI_X_ZN_(jl,jk,idt))= pxtte(_RI_X_ZN_(jl,jk,idt)) &
             +(zaerml(jl,jk,iso4cs)/cair(jl,jk)- &
             sulfmass_cs(jl,jk))*zqtmst
          ENDIF
        idt =idt_ms4cs
        pxtte(_RI_X_ZN_(jl,jk,idt))=(zaerml(jl,jk,iso4cs)/cair(jl,jk)-&
             pxtm1(_RI_X_ZN_(jl,jk,idt)))*zqtmst

        !--- 5.2.2) Black Carbon:

        idt =idt_mbcks
        pxtte(_RI_X_ZN_(jl,jk,idt))=(zaerml(jl,jk,ibcks)* 1.E-6/ crhoair(jl,jk)&
             /ti_gp(idt_mbcks)%tp%meta%cask_r(R_molarmass)-&
             pxtm1(_RI_X_ZN_(jl,jk,idt)))*zqtmst
        idt =idt_mbcas
        pxtte(_RI_X_ZN_(jl,jk,idt))=(zaerml(jl,jk,ibcas)* 1.E-6/ &
             crhoair(jl,jk)/ti_gp(idt_mbcas)%tp%meta%cask_r(R_molarmass)-&
             pxtm1(_RI_X_ZN_(jl,jk,idt)))*zqtmst
        idt =idt_mbccs
        pxtte(_RI_X_ZN_(jl,jk,idt))=(zaerml(jl,jk,ibccs)* 1.E-6/&
             crhoair(jl,jk)/ti_gp(idt_mbccs)%tp%meta%cask_r(R_molarmass)-&
             pxtm1(_RI_X_ZN_(jl,jk,idt)))*zqtmst
        idt =idt_mbcki
        pxtte(_RI_X_ZN_(jl,jk,idt))=(zaerml(jl,jk,ibcki)* 1.E-6/&
             crhoair(jl,jk)/ti_gp(idt_mbcki)%tp%meta%cask_r(R_molarmass)-&
             pxtm1(_RI_X_ZN_(jl,jk,idt)))*zqtmst

        !--- 5.2.3) Organic Carbon:

        idt = idt_mocks
        pxtte(_RI_X_ZN_(jl,jk,idt))=(zaerml(jl,jk,iocks)* 1.E-6/ &
             crhoair(jl,jk)/ti_gp(idt_mocks)%tp%meta%cask_r(R_molarmass)-&
             pxtm1(_RI_X_ZN_(jl,jk,idt)))*zqtmst
        idt = idt_mocas
        pxtte(_RI_X_ZN_(jl,jk,idt))=(zaerml(jl,jk,iocas)* 1.E-6/&
             crhoair(jl,jk)/ti_gp(idt_mocas)%tp%meta%cask_r(R_molarmass)-&
             pxtm1(_RI_X_ZN_(jl,jk,idt)))*zqtmst
        idt = idt_moccs
        pxtte(_RI_X_ZN_(jl,jk,idt))=(zaerml(jl,jk,ioccs)* 1.E-6/&
             crhoair(jl,jk)/ti_gp(idt_moccs)%tp%meta%cask_r(R_molarmass)-&
             pxtm1(_RI_X_ZN_(jl,jk,idt)))*zqtmst
        idt = idt_mocki
        pxtte(_RI_X_ZN_(jl,jk,idt))=(zaerml(jl,jk,iocki)* 1.E-6/&
             crhoair(jl,jk)/ti_gp(idt_mocki)%tp%meta%cask_r(R_molarmass)-&
             pxtm1(_RI_X_ZN_(jl,jk,idt)))*zqtmst
        
        !--- 5.2.4) Sea Salt:
        
        idt = idt_mssas
        pxtte(_RI_X_ZN_(jl,jk,idt))=(zaerml(jl,jk,issas)* 1.E-6/&
             crhoair(jl,jk)/ti_gp(idt_mssas)%tp%meta%cask_r(R_molarmass)-&
             pxtm1(_RI_X_ZN_(jl,jk,idt)))*zqtmst

        idt = idt_msscs
        pxtte(_RI_X_ZN_(jl,jk,idt))=(zaerml(jl,jk,isscs)* 1.E-6/&
             crhoair(jl,jk)/ti_gp(idt_msscs)%tp%meta%cask_r(R_molarmass)-&
             pxtm1(_RI_X_ZN_(jl,jk,idt)))*zqtmst
        
        !--- 5.2.5) Dust:
        
        idt = idt_mduas
        pxtte(_RI_X_ZN_(jl,jk,idt))=(zaerml(jl,jk,iduas)* 1.E-6/&
             crhoair(jl,jk)/ti_gp(idt_mduas)%tp%meta%cask_r(R_molarmass)-&
             pxtm1(_RI_X_ZN_(jl,jk,idt)))*zqtmst
        idt = idt_mducs
        pxtte(_RI_X_ZN_(jl,jk,idt))=(zaerml(jl,jk,iducs)* 1.E-6/&
             crhoair(jl,jk)/ti_gp(idt_mducs)%tp%meta%cask_r(R_molarmass)-&
             pxtm1(_RI_X_ZN_(jl,jk,idt)))*zqtmst

        idt = idt_mduai
        pxtte(_RI_X_ZN_(jl,jk,idt))=(zaerml(jl,jk,iduai)* 1.E-6/&
             crhoair(jl,jk)/ti_gp(idt_mduai)%tp%meta%cask_r(R_molarmass)-&
             pxtm1(_RI_X_ZN_(jl,jk,idt)))*zqtmst
        idt = idt_mduci
        pxtte(_RI_X_ZN_(jl,jk,idt))=(zaerml(jl,jk,iduci)* 1.E-6/&
             crhoair(jl,jk)/ti_gp(idt_mduci)%tp%meta%cask_r(R_molarmass)-&
             pxtm1(_RI_X_ZN_(jl,jk,idt)))*zqtmst
        
        !--- 5.3 Particle numbers:
        ! factor for 1/mol:  1.e6/crhoair(jl,jk)*
        
        SELECT CASE (TRIM(ti_gp(idt_ncs)%tp%ident%unit))
        CASE('mol(part)/mol')
           idt = idt_nns
           pxtte(_RI_X_ZN_(jl,jk,idt))=( 1.e6/crhoair(jl,jk)*zaernl(jl,jk,inucs)-&
                pxtm1(_RI_X_ZN_(jl,jk,idt)))*zqtmst / N_A
           idt = idt_nks
           pxtte(_RI_X_ZN_(jl,jk,idt))=( 1.e6/crhoair(jl,jk)*zaernl(jl,jk,iaits)-&
                pxtm1(_RI_X_ZN_(jl,jk,idt)))*zqtmst / N_A
           idt = idt_nas
           pxtte(_RI_X_ZN_(jl,jk,idt))=( 1.e6/crhoair(jl,jk)*zaernl(jl,jk,iaccs)-&
                pxtm1(_RI_X_ZN_(jl,jk,idt)))*zqtmst / N_A
           idt = idt_ncs
           pxtte(_RI_X_ZN_(jl,jk,idt))=( 1.e6/crhoair(jl,jk)*zaernl(jl,jk,icoas)-&
                pxtm1(_RI_X_ZN_(jl,jk,idt)))*zqtmst / N_A
           idt = idt_nki
           pxtte(_RI_X_ZN_(jl,jk,idt))=( 1.e6/crhoair(jl,jk)*zaernl(jl,jk,iaiti)-&
                pxtm1(_RI_X_ZN_(jl,jk,idt)))*zqtmst / N_A
           idt = idt_nai
           pxtte(_RI_X_ZN_(jl,jk,idt))=( 1.e6/crhoair(jl,jk)*zaernl(jl,jk,iacci)-&
                pxtm1(_RI_X_ZN_(jl,jk,idt)))*zqtmst / N_A
           idt = idt_nci
           pxtte(_RI_X_ZN_(jl,jk,idt))=( 1.e6/crhoair(jl,jk)*zaernl(jl,jk,icoai)-&
                pxtm1(_RI_X_ZN_(jl,jk,idt)))*zqtmst / N_A
        CASE('1/mol')
           idt = idt_nns
           pxtte(_RI_X_ZN_(jl,jk,idt))=( 1.e6/crhoair(jl,jk)*zaernl(jl,jk,inucs)-&
                pxtm1(_RI_X_ZN_(jl,jk,idt)))*zqtmst
           idt = idt_nks
           pxtte(_RI_X_ZN_(jl,jk,idt))=( 1.e6/crhoair(jl,jk)*zaernl(jl,jk,iaits)-&
                pxtm1(_RI_X_ZN_(jl,jk,idt)))*zqtmst
           idt = idt_nas
           pxtte(_RI_X_ZN_(jl,jk,idt))=( 1.e6/crhoair(jl,jk)*zaernl(jl,jk,iaccs)-&
                pxtm1(_RI_X_ZN_(jl,jk,idt)))*zqtmst
           idt = idt_ncs
           pxtte(_RI_X_ZN_(jl,jk,idt))=( 1.e6/crhoair(jl,jk)*zaernl(jl,jk,icoas)-&
                pxtm1(_RI_X_ZN_(jl,jk,idt)))*zqtmst
           idt = idt_nki
           pxtte(_RI_X_ZN_(jl,jk,idt))=( 1.e6/crhoair(jl,jk)*zaernl(jl,jk,iaiti)-&
                pxtm1(_RI_X_ZN_(jl,jk,idt)))*zqtmst
           idt = idt_nai
           pxtte(_RI_X_ZN_(jl,jk,idt))=( 1.e6/crhoair(jl,jk)*zaernl(jl,jk,iacci)-&
                pxtm1(_RI_X_ZN_(jl,jk,idt)))*zqtmst
           idt = idt_nci
           pxtte(_RI_X_ZN_(jl,jk,idt))=( 1.e6/crhoair(jl,jk)*zaernl(jl,jk,icoai)-&
                pxtm1(_RI_X_ZN_(jl,jk,idt)))*zqtmst
        CASE('1/kg')
           idt = idt_nns
           pxtte(_RI_X_ZN_(jl,jk,idt))=( 1.e9/crhoair(jl,jk)/M_air* &
                zaernl(jl,jk,inucs)-pxtm1(_RI_X_ZN_(jl,jk,idt)))*zqtmst
           idt = idt_nks
           pxtte(_RI_X_ZN_(jl,jk,idt))=( 1.e9/crhoair(jl,jk)/M_air* &
                zaernl(jl,jk,iaits)-pxtm1(_RI_X_ZN_(jl,jk,idt)))*zqtmst
           idt = idt_nas
           pxtte(_RI_X_ZN_(jl,jk,idt))=( 1.e9/crhoair(jl,jk)/M_air &
                *zaernl(jl,jk,iaccs)-pxtm1(_RI_X_ZN_(jl,jk,idt)))*zqtmst
           idt = idt_ncs
           pxtte(_RI_X_ZN_(jl,jk,idt))=( 1.e9/crhoair(jl,jk)/M_air* &
                zaernl(jl,jk,icoas)-pxtm1(_RI_X_ZN_(jl,jk,idt)))*zqtmst
           idt = idt_nki
           pxtte(_RI_X_ZN_(jl,jk,idt))=( 1.e9/crhoair(jl,jk)/M_air* &
                zaernl(jl,jk,iaiti)- pxtm1(_RI_X_ZN_(jl,jk,idt)))*zqtmst
           idt = idt_nai
           pxtte(_RI_X_ZN_(jl,jk,idt))=( 1.e9/crhoair(jl,jk)/M_air* &
                zaernl(jl,jk,iacci)-pxtm1(_RI_X_ZN_(jl,jk,idt)))*zqtmst
           idt = idt_nci
           pxtte(_RI_X_ZN_(jl,jk,idt))=( 1.e9/crhoair(jl,jk)/M_air* &
                zaernl(jl,jk,icoai)-pxtm1(_RI_X_ZN_(jl,jk,idt)))*zqtmst
           END SELECT
     ENDDO
  ENDDO
#else
        !--- 5.1 Gases:
        zxtte(:,:) = 0._dp  
        zxtte(1:kproma,1:nlev) = zqtmst * &
             (zgso4(1:kproma,1:nlev) - zgso4_bef(1:kproma,1:nlev))/cair(:,:)
        CALL mtend_add_l(my_handle, idt_so4_m7, px=zxtte)

        !--- 5.2) Particle mass:
        
        !--- 5.2.1) Sulfate mass:
        
        zxtte(1:kproma,1:nlev) = zqtmst / cair(:,:) * &
             (zaerml(1:kproma,1:nlev,iso4ns)-zaerml_bef(1:kproma,1:nlev,iso4ns))
        CALL mtend_add_l(my_handle, idt_ms4ns, px=zxtte)

        zxtte(1:kproma,1:nlev) = zqtmst / cair(:,:) * &
             (zaerml(1:kproma,1:nlev,iso4ks)-zaerml_bef(1:kproma,1:nlev,iso4ks))
        CALL mtend_add_l(my_handle, idt_ms4ks, px=zxtte)

        ! included for MECCA-AERO coupling
        IF (idt_h2so4_as /= 0 .AND. idt_hso4m_as/= 0 .AND. &
             idt_so4mm_as /=0) THEN
           idt =idt_h2so4_as

           zxtte(1:kproma,1:nlev) = zqtmst *                &
                (zaerml(1:kproma,1:nlev,iso4as) / cair(:,:) &
                - sulfmass_as(1:kproma,1:nlev))
           CALL mtend_add_l(my_handle, idt_h2so4_as, px=zxtte)
        ENDIF
        zxtte(1:kproma,1:nlev) = zqtmst / cair(:,:) * &
             (zaerml(1:kproma,1:nlev,iso4as)-zaerml_bef(1:kproma,1:nlev,iso4as))
        CALL mtend_add_l(my_handle, idt_ms4as, px=zxtte)

        ! included for MECCA-AERO coupling
        IF (idt_h2so4_cs /= 0) THEN
           zxtte(1:kproma,1:nlev) = zqtmst * &
                (zaerml(1:kproma,1:nlev,iso4cs) / cair(:,:) &
                - sulfmass_cs(1:kproma,1:nlev))
           CALL mtend_add_l(my_handle, idt_h2so4_cs, px=zxtte)
        ENDIF
        zxtte(1:kproma,1:nlev) = zqtmst / cair(:,:) * &
             (zaerml(1:kproma,1:nlev,iso4cs)-zaerml_bef(1:kproma,1:nlev,iso4cs))
        CALL mtend_add_l(my_handle, idt_ms4cs, px=zxtte)

        !--- 5.2.2) Black Carbon:

        zxtte(1:kproma,1:nlev) = zqtmst * 1.e-6/ crhoair(:,:) &
             /ti_gp(idt_mbcks)%tp%meta%cask_r(R_molarmass) * &
             (zaerml(1:kproma,1:nlev,ibcks)-zaerml_bef(1:kproma,1:nlev,ibcks))
        CALL mtend_add_l(my_handle, idt_mbcks, px=zxtte)

        zxtte(1:kproma,1:nlev) = zqtmst * 1.e-6/ crhoair(:,:)&
             /ti_gp(idt_mbcas)%tp%meta%cask_r(R_molarmass) * &
             (zaerml(1:kproma,1:nlev,ibcas)-zaerml_bef(1:kproma,1:nlev,ibcas))
        CALL mtend_add_l(my_handle, idt_mbcas, px=zxtte)

        zxtte(1:kproma,1:nlev) = zqtmst * 1.e-6/ crhoair(:,:)&
             /ti_gp(idt_mbccs)%tp%meta%cask_r(R_molarmass) * &
             (zaerml(1:kproma,1:nlev,ibccs)-zaerml_bef(1:kproma,1:nlev,ibccs))
        CALL mtend_add_l(my_handle, idt_mbccs, px=zxtte)

        zxtte(1:kproma,1:nlev) = zqtmst * 1.e-6/ crhoair(:,:)  &
             /ti_gp(idt_mbcki)%tp%meta%cask_r(R_molarmass) * &
             (zaerml(1:kproma,1:nlev,ibcki)-zaerml_bef(1:kproma,1:nlev,ibcki))
        CALL mtend_add_l(my_handle, idt_mbcki, px=zxtte)

        !--- 5.2.3) Organic Carbon:

        zxtte(1:kproma,1:nlev) = zqtmst * 1.e-6/ crhoair(:,:)&
             /ti_gp(idt_mocks)%tp%meta%cask_r(R_molarmass) * &
             (zaerml(1:kproma,1:nlev,iocks)-zaerml_bef(1:kproma,1:nlev,iocks))
        CALL mtend_add_l(my_handle, idt_mocks, px=zxtte)

        zxtte(1:kproma,1:nlev) = zqtmst * 1.e-6/ crhoair(:,:)&
             /ti_gp(idt_mocas)%tp%meta%cask_r(R_molarmass) * &
             (zaerml(1:kproma,1:nlev,iocas)-zaerml_bef(1:kproma,1:nlev,iocas))
        CALL mtend_add_l(my_handle, idt_mocas, px=zxtte)

        zxtte(1:kproma,1:nlev) = zqtmst * 1.e-6/ crhoair(:,:)&
             /ti_gp(idt_moccs)%tp%meta%cask_r(R_molarmass) * &
             (zaerml(1:kproma,1:nlev,ioccs)-zaerml_bef(1:kproma,1:nlev,ioccs))
        CALL mtend_add_l(my_handle, idt_moccs , px=zxtte)

        zxtte(1:kproma,1:nlev) = zqtmst * 1.e-6/ crhoair(:,:)&
             /ti_gp(idt_mocki)%tp%meta%cask_r(R_molarmass) * &
             (zaerml(1:kproma,1:nlev,iocki)-zaerml_bef(1:kproma,1:nlev,iocki))
        CALL mtend_add_l(my_handle, idt_mocki, px=zxtte)

        !--- 5.2.4) Sea Salt:
        
        zxtte(1:kproma,1:nlev) = zqtmst * 1.e-6/ crhoair(:,:)  &
             /ti_gp(idt_mssas)%tp%meta%cask_r(R_molarmass) * &
             (zaerml(1:kproma,1:nlev,issas)-zaerml_bef(1:kproma,1:nlev,issas))
        CALL mtend_add_l(my_handle, idt_mssas, px=zxtte)

        zxtte(1:kproma,1:nlev) = zqtmst * 1.e-6/ crhoair(:,:)  &
             /ti_gp(idt_msscs)%tp%meta%cask_r(R_molarmass) * &
             (zaerml(1:kproma,1:nlev,isscs)-zaerml_bef(1:kproma,1:nlev,isscs))
        CALL mtend_add_l(my_handle, idt_msscs, px=zxtte)

        !--- 5.2.5) Dust:
        
        zxtte(1:kproma,1:nlev) = zqtmst * 1.e-6/ crhoair(:,:)&
             /ti_gp(idt_mduas)%tp%meta%cask_r(R_molarmass) * &
             (zaerml(1:kproma,1:nlev,iduas)-zaerml_bef(1:kproma,1:nlev,iduas))
        CALL mtend_add_l(my_handle, idt_mduas, px=zxtte)

        zxtte(1:kproma,1:nlev) = zqtmst * 1.e-6/ crhoair(:,:)&
             /ti_gp(idt_mducs)%tp%meta%cask_r(R_molarmass) * &
             (zaerml(1:kproma,1:nlev,iducs)-zaerml_bef(1:kproma,1:nlev,iducs))
        CALL mtend_add_l(my_handle, idt_mducs, px=zxtte)

        zxtte(1:kproma,1:nlev) = zqtmst * 1.e-6/ crhoair(:,:)&
             /ti_gp(idt_mduai)%tp%meta%cask_r(R_molarmass) * &
             (zaerml(1:kproma,1:nlev,iduai)-zaerml_bef(1:kproma,1:nlev,iduai))
        CALL mtend_add_l(my_handle, idt_mduai, px=zxtte)

        zxtte(1:kproma,1:nlev) = zqtmst * 1.e-6/ crhoair(:,:)&
             /ti_gp(idt_mduci)%tp%meta%cask_r(R_molarmass) * &
             (zaerml(1:kproma,1:nlev,iduci)-zaerml_bef(1:kproma,1:nlev,iduci))
        CALL mtend_add_l(my_handle, idt_mduci, px=zxtte)

        !--- 5.3 Particle numbers:
        ! factor for 1/mol:  1.e6/crhoair(jl,jk)*
        
        SELECT CASE (TRIM(ti_gp(idt_ncs)%tp%ident%unit))
        CASE('mol(part)/mol')

           zxtte(1:kproma,1:nlev) = zqtmst *1.e6 / crhoair(:,:) / N_A * &
             (zaernl(1:kproma,1:nlev,inucs)-zaernl_bef(1:kproma,1:nlev,inucs))
           CALL mtend_add_l(my_handle, idt_nns, px=zxtte)

           zxtte(1:kproma,1:nlev) = zqtmst *1.e6 / crhoair(:,:) / N_A * &
             (zaernl(1:kproma,1:nlev,iaits)-zaernl_bef(1:kproma,1:nlev,iaits))
           CALL mtend_add_l(my_handle, idt_nks, px=zxtte)

           zxtte(1:kproma,1:nlev) = zqtmst *1.e6 / crhoair(:,:) / N_A * &
             (zaernl(1:kproma,1:nlev,iaccs)-zaernl_bef(1:kproma,1:nlev,iaccs))
           CALL mtend_add_l(my_handle, idt_nas , px=zxtte)

           zxtte(1:kproma,1:nlev) = zqtmst *1.e6 / crhoair(:,:) / N_A * &
             (zaernl(1:kproma,1:nlev,icoas)-zaernl_bef(1:kproma,1:nlev,icoas))
           CALL mtend_add_l(my_handle,idt_ncs, px=zxtte)

           zxtte(1:kproma,1:nlev) = zqtmst *1.e6 / crhoair(:,:) / N_A * &
             (zaernl(1:kproma,1:nlev,iaiti)-zaernl_bef(1:kproma,1:nlev,iaiti))
           CALL mtend_add_l(my_handle, idt_nki, px=zxtte )

           zxtte(1:kproma,1:nlev) = zqtmst *1.e6 / crhoair(:,:) / N_A * &
             (zaernl(1:kproma,1:nlev,iacci)-zaernl_bef(1:kproma,1:nlev,iacci))
           CALL mtend_add_l(my_handle, idt_nai, px=zxtte )

           zxtte(1:kproma,1:nlev) = zqtmst *1.e6 / crhoair(:,:) / N_A * &
             (zaernl(1:kproma,1:nlev,icoai)-zaernl_bef(1:kproma,1:nlev,icoai))
           CALL mtend_add_l(my_handle, idt_nci, px=zxtte )

        CASE('1/mol')

           zxtte(1:kproma,1:nlev) = zqtmst *1.e6 / crhoair(:,:) * &
             (zaernl(1:kproma,1:nlev,inucs)-zaernl_bef(1:kproma,1:nlev,inucs))
           CALL mtend_add_l(my_handle, idt_nns, px=zxtte )

           zxtte(1:kproma,1:nlev) = zqtmst *1.e6 / crhoair(:,:) * &
              (zaernl(1:kproma,1:nlev,iaits)-zaernl_bef(1:kproma,1:nlev,iaits))
           CALL mtend_add_l(my_handle, idt_nks, px=zxtte )

          zxtte(1:kproma,1:nlev) = zqtmst *1.e6 / crhoair(:,:) * &
              (zaernl(1:kproma,1:nlev,iaccs)-zaernl_bef(1:kproma,1:nlev,iaccs))
           CALL mtend_add_l(my_handle, idt_nas, px=zxtte )

           zxtte(1:kproma,1:nlev) = zqtmst *1.e6 / crhoair(:,:) * &
              (zaernl(1:kproma,1:nlev,icoas)-zaernl_bef(1:kproma,1:nlev,icoas))
           CALL mtend_add_l(my_handle, idt_ncs, px=zxtte )

           zxtte(1:kproma,1:nlev) = zqtmst *1.e6 / crhoair(:,:) * &
               (zaernl(1:kproma,1:nlev,iaiti)-zaernl_bef(1:kproma,1:nlev,iaiti))
           CALL mtend_add_l(my_handle, idt_nki, px=zxtte )

           zxtte(1:kproma,1:nlev) = zqtmst *1.e6 / crhoair(:,:) * &
             (zaernl(1:kproma,1:nlev,iacci)-zaernl_bef(1:kproma,1:nlev,iacci))
           CALL mtend_add_l(my_handle, idt_nai, px=zxtte )

           zxtte(1:kproma,1:nlev) = zqtmst *1.e6 / crhoair(:,:) * &
               (zaernl(1:kproma,1:nlev,icoai)-zaernl_bef(1:kproma,1:nlev,icoai))
           CALL mtend_add_l(my_handle, idt_nci, px=zxtte )

        CASE('1/kg')

           zxtte(1:kproma,1:nlev) = zqtmst *1.e9 / crhoair(:,:)/M_air * &
               (zaernl(1:kproma,1:nlev,inucs)-zaernl_bef(1:kproma,1:nlev,inucs))
           CALL mtend_add_l(my_handle, idt_nns, px=zxtte)

           zxtte(1:kproma,1:nlev) = zqtmst *1.e9 / crhoair(:,:)/M_air * &
               (zaernl(1:kproma,1:nlev,iaits)-zaernl_bef(1:kproma,1:nlev,iaits))
           CALL mtend_add_l(my_handle, idt_nks, px=zxtte)

           zxtte(1:kproma,1:nlev) = zqtmst *1.e9 / crhoair(:,:)/M_air * &
               (zaernl(1:kproma,1:nlev,iaccs)-zaernl_bef(1:kproma,1:nlev,iaccs))
           CALL mtend_add_l(my_handle, idt_nas, px=zxtte)

           zxtte(1:kproma,1:nlev) = zqtmst *1.e9 / crhoair(:,:)/M_air * &
               (zaernl(1:kproma,1:nlev,icoas)-zaernl_bef(1:kproma,1:nlev,icoas))
           CALL mtend_add_l(my_handle, idt_ncs, px=zxtte)

           zxtte(1:kproma,1:nlev) = zqtmst *1.e9 / crhoair(:,:)/M_air * &
               (zaernl(1:kproma,1:nlev,iaiti)-zaernl_bef(1:kproma,1:nlev,iaiti))
           CALL mtend_add_l(my_handle,idt_nki, px=zxtte)

           zxtte(1:kproma,1:nlev) = zqtmst *1.e9 / crhoair(:,:)/M_air * &
               (zaernl(1:kproma,1:nlev,iacci)-zaernl_bef(1:kproma,1:nlev,iacci))
           CALL mtend_add_l(my_handle, idt_nai, px=zxtte)

           zxtte(1:kproma,1:nlev) = zqtmst *1.e9 / crhoair(:,:)/M_air * &
               (zaernl(1:kproma,1:nlev,icoai)-zaernl_bef(1:kproma,1:nlev,icoai))
           CALL mtend_add_l(my_handle, idt_nci, px=zxtte)

           END SELECT
#endif

  !--- 6)) Convert M7 quantities to SI units and store in channels:

  DO jk=1, nlev
     DO jl=1, kproma
        
        
        !--- Ambient Count Mean Radius from [cm] to [m]:
        !--- Mean mode density from [g/cm3] to [kg/m3]:
        
        jm = 1
        wetrad_4d(_RI_XYZN_(jl,jrow,jk,jm)) = zm6rp(jl,jk,1)/100.
        dryrad_4d(_RI_XYZN_(jl,jrow,jk,jm)) = zm6dry(jl,jk,1)/100.
        densaer_4d(_RI_XYZN_(jl,jrow,jk,jm)) = zrhop(jl,jk,1)*1.E3

        jm = 2
        wetrad_4d(_RI_XYZN_(jl,jrow,jk,jm)) = zm6rp(jl,jk,2)/100.
        dryrad_4d(_RI_XYZN_(jl,jrow,jk,jm)) = zm6dry(jl,jk,2)/100.
        densaer_4d(_RI_XYZN_(jl,jrow,jk,jm)) = zrhop(jl,jk,2)*1.E3

        jm = 3
        wetrad_4d(_RI_XYZN_(jl,jrow,jk,jm)) = zm6rp(jl,jk,3)/100.
        dryrad_4d(_RI_XYZN_(jl,jrow,jk,jm)) = zm6dry(jl,jk,3)/100.
        densaer_4d(_RI_XYZN_(jl,jrow,jk,jm)) = zrhop(jl,jk,3)*1.E3

        jm = 4    
        wetrad_4d(_RI_XYZN_(jl,jrow,jk,jm)) = zm6rp(jl,jk,4)/100.
        dryrad_4d(_RI_XYZN_(jl,jrow,jk,jm)) = zm6dry(jl,jk,4)/100.
        densaer_4d(_RI_XYZN_(jl,jrow,jk,jm)) = zrhop(jl,jk,4)*1.E3

        jm = 5    
        wetrad_4d(_RI_XYZN_(jl,jrow,jk,jm)) = zm6rp(jl,jk,5)/100.
        dryrad_4d(_RI_XYZN_(jl,jrow,jk,jm)) = zm6rp(jl,jk,5)/100.
        densaer_4d(_RI_XYZN_(jl,jrow,jk,jm)) = zrhop(jl,jk,5)*1.E3

        jm = 6    
        wetrad_4d(_RI_XYZN_(jl,jrow,jk,jm)) = zm6rp(jl,jk,6)/100.
        dryrad_4d(_RI_XYZN_(jl,jrow,jk,jm)) = zm6rp(jl,jk,6)/100.
        densaer_4d(_RI_XYZN_(jl,jrow,jk,jm)) = zrhop(jl,jk,6)*1.E3

        jm = 7    
        wetrad_4d(_RI_XYZN_(jl,jrow,jk,jm)) = zm6rp(jl,jk,7)/100.
        dryrad_4d(_RI_XYZN_(jl,jrow,jk,jm)) = zm6rp(jl,jk,7)/100.
        densaer_4d(_RI_XYZN_(jl,jrow,jk,jm)) = zrhop(jl,jk,7)*1.E3

     END DO
  END DO

  !--- 7) Perform mass conservation check: 

  IF(lmass_diag) THEN

     !--- Sum total mass of all compounds for mass diagnostics:

     CALL sum_mass(zmass_post)

     !--- Perform mass conservation check:

     DO jc=1, 5

        IF( ABS(zmass_pre(jc)-zmass_post(jc)) > &
             0.100*ABS(MAX(zmass_pre(jc),zmass_post(jc)))) THEN

           WRITE(*,*) 'm7_interface: ',&
          'm7 violates the mass conservation for compound number = ', jc, ">10%"

           IF (MAX(zmass_pre(jc),zmass_post(jc)) >0.) THEN
              WRITE(*,*) '              change in mass : ', &
                   (zmass_post(jc)-zmass_pre(jc))/&
                   MAX(zmass_pre(jc),zmass_post(jc))*100.,"%"
           ELSE
              WRITE(*,*) '              change in mass : ', &
                            zmass_post(jc)-zmass_pre(jc)
           END IF

           loabort(jc)=.TRUE.        
        
        ELSE IF( ABS(zmass_pre(jc)-zmass_post(jc)) > &
             0.010*ABS(MAX(zmass_pre(jc),zmass_post(jc)))) THEN

           WRITE(*,*) 'm7_interface: m7 violates the mass conservation for compound number = ', jc, ">1%"

           IF (MAX(zmass_pre(jc),zmass_post(jc)) >0.) THEN
              WRITE(*,*) '              change in mass : ', &
                            (zmass_post(jc)-zmass_pre(jc))/&
                            MAX(zmass_pre(jc),zmass_post(jc))*100.,"%"
           ELSE
              WRITE(*,*) '              change in mass : ', &
                            zmass_post(jc)-zmass_pre(jc)
           END IF

        ELSE IF( ABS(zmass_pre(jc)-zmass_post(jc)) > &
             0.001*ABS(MAX(zmass_pre(jc),zmass_post(jc)))) THEN

           WRITE(*,*) 'm7_interface: m7 violates the mass conservation for compound number = ', jc, ">0.1%"

           IF (MAX(zmass_pre(jc),zmass_post(jc)) >0.) THEN
              WRITE(*,*) '              change in mass : ', &
                            (zmass_post(jc)-zmass_pre(jc))/&
                            MAX(zmass_pre(jc),zmass_post(jc))*100.,"%"
           ELSE
              WRITE(*,*) '              change in mass : ', &
                            zmass_post(jc)-zmass_pre(jc)
           END IF

        END IF

     END DO

     IF (ANY(loabort))THEN
        WRITE(*,*) 'Mass conservation failure for compounds:', loabort
        CALL error_bi('m7 violates the mass conservation','m7_interface')
     END IF

  END IF

  !--- Auxiliary routines:

  CONTAINS

    SUBROUTINE sum_mass(pmasssum)

      REAL(dp) :: pmasssum(5)
      INTEGER  :: jk, jl

      pmasssum=0.

      !--- Store channel in local array for optimization:

      zgboxarea(1:kproma)=gboxarea_2d(1:kproma,jrow)

      DO jk=1,nlev
         DO jl=1,kproma
            pmasssum(1)=pmasssum(1) + (zgso4(jl,jk)         + &
                                       zaerml(jl,jk,iso4ns) + &
                                       zaerml(jl,jk,iso4ks) + &
                                       zaerml(jl,jk,iso4as) + &
                                       zaerml(jl,jk,iso4cs))* &
                                       zdpg(jl,jk)*zgboxarea(jl)

            pmasssum(2)=pmasssum(2) + (zaerml(jl,jk,ibcks)  + &
                                       zaerml(jl,jk,ibcas)  + &
                                       zaerml(jl,jk,ibccs)  + &
                                       zaerml(jl,jk,ibcki))*  &
                                       zdpg(jl,jk)*zgboxarea(jl)

            pmasssum(3)=pmasssum(3) + (zaerml(jl,jk,iocks)  + &
                                       zaerml(jl,jk,iocas)  + &
                                       zaerml(jl,jk,ioccs)  + &
                                       zaerml(jl,jk,iocki))*  &
                                       zdpg(jl,jk)*zgboxarea(jl)

            pmasssum(4)=pmasssum(4) + (zaerml(jl,jk,issas)  + &
                                       zaerml(jl,jk,isscs))*  &
                                       zdpg(jl,jk)*zgboxarea(jl)

            pmasssum(5)=pmasssum(5) + (zaerml(jl,jk,iduas)  + &
                                       zaerml(jl,jk,iducs)  + &
                                       zaerml(jl,jk,iduai)  + &
                                       zaerml(jl,jk,iduci))*  &
                                       zdpg(jl,jk)*zgboxarea(jl)
         ENDDO
      ENDDO

    END SUBROUTINE sum_mass

END SUBROUTINE m7_physc
!-----------------------------------------------------------------------------

!-----------------------------------------------------------------------------

SUBROUTINE m7_local_end
  !-----
  ! 
  ! 
  !-----

  ! ECHAM5/MESSy
  USE messy_main_grid_def_mem_bi, ONLY: jrow, kproma, nproma, nlev
  USE messy_main_data_bi,       ONLY: pressi_3d,rhum_3d
  USE messy_main_timer,         ONLY: time_step_len
  USE messy_main_tracer_mem_bi, ONLY: pxtte => qxtte, pxtm1 => qxtm1
  USE messy_main_tracer,        ONLY: R_molarmass
  USE messy_main_constants_mem, ONLY: M_air,g

  IMPLICIT NONE

  ! LOCAL
  REAL(dp)      :: ztmst, g_M_air
  REAL(dp)      :: moltomg_SO4, moltomg_BC, moltomg_OC, moltomg_SS, moltomg_DU
  REAL(dp)      :: trac_mass(nproma,nlev,5)
  REAL(dp)      :: zdp(nproma,nlev), aod_f_rh(nproma,nlev)
  ! exponents used to calculate aod_f_rh empirically
  REAL(dp), PARAMETER :: aod_ap0       =  0.196581266_dp
  REAL(dp), PARAMETER :: aod_ap1       =  0.100650483_dp
  REAL(dp), PARAMETER :: aod_ap2       = -0.00389645252_dp
  REAL(dp), PARAMETER :: aod_ap3       =  5.63596671e-5_dp
  REAL(dp), PARAMETER :: aod_ap4       =  6.4988432e-8_dp
  REAL(dp), PARAMETER :: aod_ap5       = -7.95507834e-9_dp
  REAL(dp), PARAMETER :: aod_ap6       =  4.95445298e-11_dp
  ! Alpha is the assumed mass extinction effciency of each species (m2 mg-1)
  ! Alpha values from AeroCom median 
  ! (Kinne ACP - Table 4, An AeroCom initial assesment)
  REAL(dp), PARAMETER :: aod_alpha_SO4 =  0.0085_dp
  REAL(dp), PARAMETER :: aod_alpha_BC  =  0.0089_dp
  REAL(dp), PARAMETER :: aod_alpha_OC  =  0.0057_dp
  REAL(dp), PARAMETER :: aod_alpha_SS  =  0.003_dp
  REAL(dp), PARAMETER :: aod_alpha_DU  =  0.00095_dp
  INTEGER   :: jk
  INTEGER   :: idt


  IF (l_colmass .or. l_aod) THEN

     ! pressure difference per box
     zdp(1:kproma,1:nlev) = pressi_3d(_RI_XYZ__(1:kproma,jrow,2:nlev+1)) &
          - pressi_3d(_RI_XYZ__(1:kproma,jrow,1:nlev)) 
     
     ztmst = time_step_len
     
     g_M_air = 1.E6/(g*M_air)
     
     moltomg_SO4 = ti_gp(idt_ms4as)%tp%meta%cask_r(R_molarmass)*g_M_air
     moltomg_BC  = ti_gp(idt_mbcas)%tp%meta%cask_r(R_molarmass)*g_M_air
     moltomg_OC  = ti_gp(idt_mocas)%tp%meta%cask_r(R_molarmass)*g_M_air
     moltomg_SS  = ti_gp(idt_mssas)%tp%meta%cask_r(R_molarmass)*g_M_air
     moltomg_DU  = ti_gp(idt_mduas)%tp%meta%cask_r(R_molarmass)*g_M_air
     
     
     ! 3D tracer mass in mg m-2
     ! SO4
     idt = idt_ms4ns
     trac_mass(1:kproma,1:nlev,1) = zdp(1:kproma,:) * moltomg_SO4 *           &
          (pxtm1(_RI_X_ZN_(1:kproma,1:nlev,idt)) + pxtte(_RI_X_ZN_(1:kproma,1:nlev,idt))*ztmst)
     idt = idt_ms4ks
     trac_mass(1:kproma,1:nlev,1) = trac_mass(1:kproma,1:nlev,1) +            &
          zdp(1:kproma,:) * moltomg_SO4 *                                     &
          (pxtm1(_RI_X_ZN_(1:kproma,1:nlev,idt)) +  pxtte(_RI_X_ZN_(1:kproma,1:nlev,idt)) *ztmst) 
     idt = idt_ms4as
     trac_mass(1:kproma,1:nlev,1) = trac_mass(1:kproma,1:nlev,1) +            &
          zdp(1:kproma,:) * moltomg_SO4 *                                     &
          (pxtm1(_RI_X_ZN_(1:kproma,1:nlev,idt)) + pxtte(_RI_X_ZN_(1:kproma,1:nlev,idt))*ztmst)
     idt = idt_ms4cs
     trac_mass(1:kproma,1:nlev,1) = trac_mass(1:kproma,1:nlev,1) +            &
          zdp(1:kproma,:) * moltomg_SO4 *                                     &
          ( pxtm1(_RI_X_ZN_(1:kproma,1:nlev,idt)) + pxtte(_RI_X_ZN_(1:kproma,1:nlev,idt))*ztmst)
     ! BC
     idt = idt_mbcks
     trac_mass(1:kproma,1:nlev,2) = zdp(1:kproma,:) * moltomg_BC *            &
          (pxtm1(_RI_X_ZN_(1:kproma,1:nlev,idt))+ pxtte(_RI_X_ZN_(1:kproma,1:nlev,idt))*ztmst)
     idt = idt_mbcas
     trac_mass(1:kproma,1:nlev,2) =      trac_mass(1:kproma,1:nlev,2)  +      &
          zdp(1:kproma,:) * moltomg_BC *                                      &
          (pxtm1(_RI_X_ZN_(1:kproma,1:nlev,idt))+ pxtte(_RI_X_ZN_(1:kproma,1:nlev,idt))*ztmst)
     idt = idt_mbccs
     trac_mass(1:kproma,1:nlev,2) =      trac_mass(1:kproma,1:nlev,2)  +      &
          zdp(1:kproma,:) * moltomg_BC *                                      &
          ( pxtm1(_RI_X_ZN_(1:kproma,1:nlev,idt)) + pxtte(_RI_X_ZN_(1:kproma,1:nlev,idt)) *ztmst)
     idt = idt_mbcki
     trac_mass(1:kproma,1:nlev,2) =      trac_mass(1:kproma,1:nlev,2)  +      & 
          zdp(1:kproma,:) * moltomg_BC *                                      & 
          ( pxtm1(_RI_X_ZN_(1:kproma,1:nlev,idt)) + pxtte(_RI_X_ZN_(1:kproma,1:nlev,idt))*ztmst)
     ! OC
     idt =idt_mocks 
     trac_mass(1:kproma,1:nlev,3) = zdp(1:kproma,:) * moltomg_OC *                 &
          (pxtm1(_RI_X_ZN_(1:kproma,1:nlev,idt)) + pxtte(_RI_X_ZN_(1:kproma,1:nlev,idt))*ztmst)
     idt = idt_mocas
     trac_mass(1:kproma,1:nlev,3) = trac_mass(1:kproma,1:nlev,3) +       &
          zdp(1:kproma,:) * moltomg_OC *                       &
          (pxtm1(_RI_X_ZN_(1:kproma,1:nlev,idt)) + pxtte(_RI_X_ZN_(1:kproma,1:nlev,idt))*ztmst)
     idt = idt_moccs
     trac_mass(1:kproma,1:nlev,3) = trac_mass(1:kproma,1:nlev,3) + &
          zdp(1:kproma,:) * moltomg_OC *                 &
          (pxtm1(_RI_X_ZN_(1:kproma,1:nlev,idt)) +  pxtte(_RI_X_ZN_(1:kproma,1:nlev,idt)) *ztmst)
     idt = idt_mocki
     trac_mass(1:kproma,1:nlev,3) = trac_mass(1:kproma,1:nlev,3) + &
          zdp(1:kproma,:) * moltomg_OC *                 &
          (pxtm1(_RI_X_ZN_(1:kproma,1:nlev,idt)) + pxtte(_RI_X_ZN_(1:kproma,1:nlev,idt)) *ztmst)
     ! SS
     idt = idt_mssas
     trac_mass(1:kproma,1:nlev,4) = zdp(1:kproma,:) * moltomg_SS *                 &
          (pxtm1(_RI_X_ZN_(1:kproma,1:nlev,idt)) + pxtte(_RI_X_ZN_(1:kproma,1:nlev,idt))*ztmst)
     idt = idt_msscs
     trac_mass(1:kproma,1:nlev,4) = trac_mass(1:kproma,1:nlev,4) + &
          zdp(1:kproma,:) * moltomg_SS *                 &
          (pxtm1(_RI_X_ZN_(1:kproma,1:nlev,idt)) + pxtte(_RI_X_ZN_(1:kproma,1:nlev,idt))*ztmst)

     ! DU
     idt = idt_mduas
     trac_mass(1:kproma,1:nlev,5) = zdp(1:kproma,:) * moltomg_DU *  &     
          (pxtm1(_RI_X_ZN_(1:kproma,1:nlev,idt)) + pxtte(_RI_X_ZN_(1:kproma,1:nlev,idt))*ztmst)
     idt = idt_mducs
     trac_mass(1:kproma,1:nlev,5) = trac_mass(1:kproma,1:nlev,5) +  &
          zdp(1:kproma,:) * moltomg_DU *                 &
          (pxtm1(_RI_X_ZN_(1:kproma,1:nlev,idt)) + pxtte(_RI_X_ZN_(1:kproma,1:nlev,idt))*ztmst) 
     idt = idt_mduai
     trac_mass(1:kproma,1:nlev,5) = trac_mass(1:kproma,1:nlev,5) +  &
          zdp(1:kproma,:) * moltomg_DU *                 &
          (pxtm1(_RI_X_ZN_(1:kproma,1:nlev,idt)) + pxtte(_RI_X_ZN_(1:kproma,1:nlev,idt))*ztmst)
     idt = idt_mduci
     trac_mass(1:kproma,1:nlev,5) = trac_mass(1:kproma,1:nlev,5) +  & 
          zdp(1:kproma,:) * moltomg_DU *                 &
          ( pxtm1(_RI_X_ZN_(1:kproma,1:nlev,idt)) + pxtte(_RI_X_ZN_(1:kproma,1:nlev,idt))*ztmst)


     IF (l_colmass) THEN 

        colmass_SO4(1:kproma,jrow)   = 0._dp
        colmass_BC(1:kproma,jrow)    = 0._dp
        colmass_OC(1:kproma,jrow)    = 0._dp
        colmass_SS(1:kproma,jrow)    = 0._dp
        colmass_DU(1:kproma,jrow)    = 0._dp
        colmass_DU_as(1:kproma,jrow) = 0._dp
        colmass_DU_cs(1:kproma,jrow) = 0._dp
        colmass_DU_ai(1:kproma,jrow) = 0._dp
        colmass_DU_ci(1:kproma,jrow) = 0._dp

        DO jk=1,nlev
           ! SO4
           colmass_SO4(1:kproma,jrow) = colmass_SO4(1:kproma,jrow) +          &
                trac_mass(1:kproma,jk,1)
           ! BC
           colmass_BC(1:kproma,jrow) = colmass_BC(1:kproma,jrow) +            &
                trac_mass(1:kproma,jk,2)
           ! OC
           colmass_OC(1:kproma,jrow) = colmass_OC(1:kproma,jrow) +            &
                trac_mass(1:kproma,jk,3)
           ! SS
           colmass_SS(1:kproma,jrow) = colmass_SS(1:kproma,jrow) +            &
                trac_mass(1:kproma,jk,4)
           ! DU
           colmass_DU(1:kproma,jrow) = colmass_DU(1:kproma,jrow) +            &
                trac_mass(1:kproma,jk,5)

           idt = idt_mduas
           colmass_DU_as(1:kproma,jrow) = colmass_DU_as(1:kproma,jrow) +      &
             zdp(1:kproma,jk) * moltomg_DU *                                  &
             (pxtm1(_RI_X_ZN_(1:kproma,jk,idt)) + pxtte(_RI_X_ZN_(1:kproma,jk,idt))*ztmst)
           idt = idt_mducs
           colmass_DU_cs(1:kproma,jrow) = colmass_DU_cs(1:kproma,jrow) +      &
             zdp(1:kproma,jk) * moltomg_DU *  &
             (pxtm1(_RI_X_ZN_(1:kproma,jk,idt)) + pxtte(_RI_X_ZN_(1:kproma,jk,idt))*ztmst)
           idt = idt_mduai
           colmass_DU_ai(1:kproma,jrow) = colmass_DU_ai(1:kproma,jrow) +      &
             zdp(1:kproma,jk) * moltomg_DU * &
             (pxtm1(_RI_X_ZN_(1:kproma,jk,idt)) +  pxtte(_RI_X_ZN_(1:kproma,jk,idt))*ztmst)
           idt = idt_mduci 
           colmass_DU_ci(1:kproma,jrow) = colmass_DU_ci(1:kproma,jrow) +      &
             zdp(1:kproma,jk) * moltomg_DU * &
             (pxtm1(_RI_X_ZN_(1:kproma,jk,idt)) +  pxtte(_RI_X_ZN_(1:kproma,jk,idt))*ztmst)
        ENDDO

     ENDIF ! l_colmass

     IF (l_aod) THEN

        AOD_SO4_COLUMN(1:kproma,jrow) = 0._dp
        AOD_BC_COLUMN(1:kproma,jrow)  = 0._dp
        AOD_OC_COLUMN(1:kproma,jrow)  = 0._dp
        AOD_SS_COLUMN(1:kproma,jrow)  = 0._dp
        AOD_DU_COLUMN(1:kproma,jrow)  = 0._dp

        aod_f_rh(1:kproma,1:nlev) =                        &
             aod_ap0*rhum_3d(_RI_XYZ__(1:kproma,jrow,1:nlev))**0. +   &
             aod_ap1*rhum_3d(_RI_XYZ__(1:kproma,jrow,1:nlev))     +   &
             aod_ap2*rhum_3d(_RI_XYZ__(1:kproma,jrow,1:nlev))**2. +   & 
             aod_ap3*rhum_3d(_RI_XYZ__(1:kproma,jrow,1:nlev))**3. +   & 
             aod_ap4*rhum_3d(_RI_XYZ__(1:kproma,jrow,1:nlev))**4. +   & 
             aod_ap5*rhum_3d(_RI_XYZ__(1:kproma,jrow,1:nlev))**5. +   & 
             aod_ap6*rhum_3d(_RI_XYZ__(1:kproma,jrow,1:nlev))**6. 

        ! AOD_species = aod_f_rh * mass extinction efficiency * mass_species in mg m-2
        ! summing over all levels gives the AOD for the whole column
        DO jk=1,nlev
           ! SO4
           AOD_SO4_COLUMN(1:kproma,jrow) = AOD_SO4_COLUMN(1:kproma,jrow) +    &
                aod_f_rh(1:kproma,jk)*aod_alpha_SO4*trac_mass(1:kproma,jk,1)
           ! BC
           AOD_BC_COLUMN(1:kproma,jrow) = AOD_BC_COLUMN(1:kproma,jrow) +      &
                aod_f_rh(1:kproma,jk)*aod_alpha_BC*trac_mass(1:kproma,jk,2)
           ! OC
           AOD_OC_COLUMN(1:kproma,jrow) = AOD_OC_COLUMN(1:kproma,jrow) +      &
                aod_f_rh(1:kproma,jk)*aod_alpha_OC*trac_mass(1:kproma,jk,3)
           ! SS
           AOD_SS_COLUMN(1:kproma,jrow) = AOD_SS_COLUMN(1:kproma,jrow) +      &
                aod_f_rh(1:kproma,jk)*aod_alpha_SS*trac_mass(1:kproma,jk,4)
           ! DU
           AOD_DU_COLUMN(1:kproma,jrow) = AOD_DU_COLUMN(1:kproma,jrow) +      &
                aod_f_rh(1:kproma,jk)*aod_alpha_DU*trac_mass(1:kproma,jk,5)
        ENDDO
           
        ! TOTAL
        AOD_TOTAL_COLUMN(1:kproma,jrow) = AOD_SO4_COLUMN(1:kproma,jrow) +     &
             AOD_BC_COLUMN(1:kproma,jrow) + AOD_OC_COLUMN(1:kproma,jrow) +    &
             AOD_SS_COLUMN(1:kproma,jrow) + AOD_DU_COLUMN(1:kproma,jrow)

     ENDIF ! l_aod

  ENDIF ! l_colmass .or. l_aod

END SUBROUTINE m7_local_end

!-----------------------------------------------------------------------------

!-----------------------------------------------------------------------------
SUBROUTINE m7_free_memory

INTRINSIC ASSOCIATED

IF (ASSOCIATED(wetradius))  DEALLOCATE(wetradius)
IF (ASSOCIATED(dryradius))  DEALLOCATE(dryradius)
IF (ASSOCIATED(densaer))    DEALLOCATE(densaer)

END SUBROUTINE m7_free_memory
!-----------------------------------------------------------------------------

!-----------------------------------------------------------------------------
  SUBROUTINE m7_read_nml_cpl(status, iou)

    ! M7 MODULE ROUTINE (INTERFACE, PRIVATE)
    ! read namelist for 'coupling' to online tracers
    ! Author: Astrid Kerkweg, MPICH, JAN 2004

    ! MESSy
    USE messy_main_tools, ONLY: read_nml_open, read_nml_check, read_nml_close

    IMPLICIT NONE

    INTRINSIC TRIM

    ! I/O
    INTEGER(I4), INTENT(OUT) :: status     ! error status
    INTEGER,     INTENT(IN)  :: iou        ! I/O unit

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'm7_read_nml_cpl'
    LOGICAL                     :: lex      ! file exists ?
    INTEGER                     :: fstat    ! file status

    ! INITIALIZE
    status = 1

    CALL read_nml_open(lex, substr, iou, 'CPL', modstr)
    IF (.NOT.lex) RETURN    ! <modstr>.nml does not exist

    READ(iou, NML=CPL, IOSTAT=fstat)
    CALL read_nml_check(fstat, substr, iou, 'CPL', modstr)
    IF (fstat /= 0) RETURN  ! error while reading namelist

    ! DIAGNOSE NAMELIST AND SET GLOBAL SWITCHES
    !
    ! CHECK NAMELIST
    IF ( (TRIM(chemmodule) == 'mecca1') .OR. &
         (TRIM(chemmodule) == 'mecca') ) THEN
       IF (lcpl_gasphase) THEN
          WRITE(*,*) 'using MECCA gas phase chemistry'
       ELSE
          WRITE(*,*) 'using M7 gas phase chemistry'
       END IF
    ELSE
       write(*,*) 'no other chemistry module than MECCA available'
       write(*,*) 'using M7 chemistry'
    ENDIF
    IF (l_calc_emis) THEN
       WRITE(*,*) ' calculation of emission required'
    ELSE
       WRITE(*,*) ' expecting emission from somewhere else ...'
    END IF

    CALL read_nml_close(substr, iou, modstr)

    status = 0  ! no ERROR

  END SUBROUTINE m7_read_nml_cpl
!-----------------------------------------------------------------------------

END MODULE messy_m7_si

! -*- f90 -*-
!==============================================================================
!
! This file controls the MESSy submodel CALLs
! Authors: Patrick Joeckel and Rolf Sander, MPICH, 2002-2007-2009
!          Patrick Joeckel, DLR, 2009-2014
!
!==============================================================================

SUBROUTINE messy_setup(flag)

  USE messy_main_switch   ! ONLY: USE_*
  USE messy_main_control  ! ONLY: entrypoint, subentry
  USE messy_main_blather_bi, ONLY: main_blather_setup
  USE messy_main_switch_bi,  ONLY: main_switch_setup !initialize
  USE messy_main_qtimer_bi,  ONLY: main_qtimer_setup !initialize
  USE messy_main_timer_bi,   ONLY: main_timer_setup  !initialize
  USE messy_main_channel_bi, ONLY: main_channel_setup, n_dom
  USE messy_main_tracer_bi,  ONLY: main_tracer_setup
  USE messy_main_decomp_bi,  ONLY: main_decomp_setup

  USE messy_clams_si,        ONLY: clams_setup
  USE messy_clamsmix_si,     ONLY: clamsmix_setup
  USE messy_clamssedi_si,    ONLY: clamssedi_setup
  USE messy_clamstraj_si,    ONLY: clamstraj_setup
  USE messy_oasis3mct_si,    ONLY: oasis3mct_setup

  IMPLICIT NONE

  INTEGER, INTENT(IN) :: flag
  ! LOCAL
  INTEGER             :: nclock

  CALL main_control_set_context(MEP_SETUP, flag)

  SELECT CASE(flag)
  CASE(0)
     CALL main_blather_setup
     CALL main_switch_setup  ! initialize ! define USE_* logicals
     CALL oasis3mct_setup(USE_OASIS3MCT)

  CASE(1)
     ! TIMER SPLIT NEEDED HERE, TO ENABLE OTHER MODELS TO OVERWRITE NAMELIST
     ! ENTRY OF START DATE AND STOP DATE (e.g. MMD_CLIENT)
     ! n_dom added to avoid cicurlar dependency
     ! additionally, the parameter needs to be INTENT(INOUT) for ICON,
     ! therefore we change to local parameter nclock here.
     nclock = n_dom
     CALL main_timer_setup(1, nclock)
     CALL main_channel_setup     ! -> channel_init_restart_bi
     CALL main_timer_setup(2 ,nclock)
     ! do NOT change the following order of main_* subroutines
     CALL main_qtimer_setup

     IF (USE_CLAMS)     CALL clams_setup
     IF (USE_CLAMSTRAJ) CALL clamstraj_setup
     IF (USE_CLAMSSEDI) CALL clamssedi_setup
     IF (USE_CLAMSMIX)  CALL clamsmix_setup

  CASE(2)
     CALL main_tracer_setup  ! -> read namelists

  CASE(3)
     CALL main_decomp_setup ! -> called separately, ngl and nlon are not
                            !    known before

  END SELECT

  CALL main_control_unset_context

END SUBROUTINE messy_setup

!==============================================================================

SUBROUTINE messy_initialize

  ! initialize submodels, read CTRL and CPL namelists (called by initialize)

  USE messy_main_switch   ! ONLY: USE_*
  USE messy_main_control  ! ONLY: entrypoint, subentry
  USE messy_main_channel_bi,ONLY: main_channel_initialize
  USE messy_main_timer_bi,  ONLY: main_timer_initialize
  USE messy_main_data_bi,   ONLY: main_data_initialize
  USE messy_main_grid_bi,   ONLY: main_grid_initialize
  USE messy_main_tracer_bi, ONLY: main_tracer_initialize
  USE messy_main_tools_bi,  ONLY: main_tools_initialize
  USE messy_main_import_bi, ONLY: main_import_initialize
#ifdef MESSYTENDENCY
  USE messy_main_tendency_bi, ONLY: main_tendency_initialize
#endif

  USE messy_a2o_e5,        ONLY:         a2o_initialize
  USE messy_accf_e5,       ONLY:        accf_initialize
  USE messy_aeropt_si,     ONLY:      aeropt_initialize
  USE messy_airsea_si,     ONLY:      airsea_initialize
  USE messy_airtrac_e5,    ONLY:     airtrac_initialize
  USE messy_airtraf_e5,    ONLY:     airtraf_initialize
  USE messy_attila_e5,     ONLY:      attila_initialize
  USE messy_aveout_si,     ONLY:      aveout_initialize
  USE messy_bioburn_si,    ONLY:     bioburn_initialize
  USE messy_bufly_si,      ONLY:       bufly_initialize
  USE messy_cat_e5,        ONLY:         cat_initialize
  USE messy_ch4_si,        ONLY:         ch4_initialize
  USE messy_clams_si,      ONLY:       clams_initialize
  USE messy_clamstraj_si,  ONLY:   clamstraj_initialize
  USE messy_clamschem_si,  ONLY:   clamschem_initialize
  USE messy_clamscheme5_si,ONLY: clamscheme5_initialize
  USE messy_clamsmix_si,   ONLY:    clamsmix_initialize
  USE messy_clamsbmix_si,  ONLY:   clamsbmix_initialize
  USE messy_clamsrdfrc_si, ONLY:  clamsrdfrc_initialize
  USE messy_clamssedi_si,  ONLY:   clamssedi_initialize
  USE messy_clamscirrus_si,ONLY: clamscirrus_initialize
  USE messy_clamstracer_si,ONLY: clamstracer_initialize
  USE messy_chemglue_si,   ONLY:    chemglue_initialize
  USE messy_cloud_si,      ONLY:       cloud_initialize
  USE messy_cloudj_si,     ONLY:      cloudj_initialize
  USE messy_cloudopt_si,   ONLY:    cloudopt_initialize
  USE messy_contrail_si,   ONLY:    contrail_initialize
  USE messy_convect_si,    ONLY:     convect_initialize
  USE messy_crm_e5,        ONLY:         crm_initialize
  USE messy_cvtrans_si,    ONLY:     cvtrans_initialize
  USE messy_d14co_si,      ONLY:       d14co_initialize
  USE messy_ddep_si,       ONLY:        ddep_initialize
  USE messy_dissoc_si,     ONLY:      dissoc_initialize
  USE messy_diumod_si,     ONLY:      diumod_initialize
  USE messy_domint_si,     ONLY:      domint_initialize
  USE messy_dradon_si,     ONLY:      dradon_initialize
  USE messy_e4chem_si,     ONLY:      e4chem_initialize
  USE messy_e5vdiff_si,    ONLY:     e5vdiff_initialize
  USE messy_edith_e5,      ONLY:       edith_initialize
  USE messy_ever_si,       ONLY:        ever_initialize
  USE messy_gec_si,        ONLY:         gec_initialize
  USE messy_gmxe_si,       ONLY:        gmxe_initialize
  USE messy_gwave_si,      ONLY:       gwave_initialize
  USE messy_h2o_e5,        ONLY:         h2o_initialize
  USE messy_h2oemis_si,    ONLY:     h2oemis_initialize
  USE messy_h2oiso_e5,     ONLY:      h2oiso_initialize
  USE messy_hd_e5,         ONLY:          hd_initialize
  USE messy_hamocc_e5,     ONLY:      hamocc_initialize
  USE messy_ions_si,       ONLY:        ions_initialize
  USE messy_isopcor_si,    ONLY:     isopcor_initialize
  USE messy_jval_si,       ONLY:        jval_initialize
  USE messy_jvst_si,       ONLY:        jvst_initialize
  USE messy_lggp_e5,       ONLY:        lggp_initialize
  USE messy_lgtmix_e5,     ONLY:      lgtmix_initialize
  USE messy_lgvflux_e5,    ONLY:     lgvflux_initialize
  USE messy_lnox_si,       ONLY:        lnox_initialize
  USE messy_m7_si,         ONLY:          m7_initialize
  USE messy_made_si,       ONLY:        made_initialize
  USE messy_made3_si,      ONLY:       made3_initialize
  USE messy_mecca_si,      ONLY:       mecca_initialize
  USE messy_megan_si,      ONLY:       megan_initialize
  USE messy_mesoenergy_e5, ONLY:  mesoenergy_initialize
  USE messy_mlocean_e5,    ONLY:     mlocean_initialize
  USE messy_mmforce_e5,    ONLY:     mmforce_initialize
  USE messy_mpiom_e5,      ONLY:       mpiom_initialize
  USE messy_msbm_si,       ONLY:        msbm_initialize
  USE messy_mtskip_si,     ONLY:      mtskip_initialize
  USE messy_nan_si,        ONLY:         nan_initialize
  USE messy_o3orig_si,     ONLY:      o3orig_initialize
  USE messy_oasis3mct_si,  ONLY:   oasis3mct_initialize
  USE messy_offemis_si,    ONLY:     offemis_initialize
  USE messy_onemis_si,     ONLY:      onemis_initialize
  USE messy_oracle_si,     ONLY:      oracle_initialize
  USE messy_orogw_si,      ONLY:       orogw_initialize
  USE messy_orbit_si,      ONLY:       orbit_initialize
  USE messy_otphysc_e5,    ONLY:     otphysc_initialize
  USE messy_plumegas_si,   ONLY:    plumegas_initialize
  USE messy_ptrac_si,      ONLY:       ptrac_initialize
  USE messy_ptracini_si,   ONLY:    ptracini_initialize
  USE messy_qbo_si,        ONLY:         qbo_initialize
  USE messy_rad_si,        ONLY:         rad_initialize
  USE messy_relax_e5,      ONLY:       relax_initialize
  USE messy_rndtest_si,    ONLY:     rndtest_initialize
  USE messy_scalc_si,      ONLY:       scalc_initialize
  USE messy_scav_si,       ONLY:        scav_initialize
  USE messy_satsims_e5,    ONLY:     satsims_initialize
  USE messy_scout_si,      ONLY:       scout_initialize
  USE messy_sedi_si,       ONLY:        sedi_initialize
  USE messy_s4d_si,        ONLY:         s4d_initialize
  USE messy_sf6_e5,        ONLY:         sf6_initialize
  USE messy_sorbit_si,     ONLY:      sorbit_initialize
  USE messy_spe_e5,        ONLY:         spe_initialize
  USE messy_spacenox_e5,   ONLY:    spacenox_initialize
  USE messy_svoc_si,       ONLY:        svoc_initialize
  USE messy_surface_si,    ONLY:     surface_initialize
  USE messy_tagging_si,    ONLY:     tagging_initialize
  USE messy_tbudget_si,    ONLY:     tbudget_initialize
  USE messy_tpulse_si,     ONLY:      tpulse_initialize
  USE messy_timepos_e5,    ONLY:     timepos_initialize
  USE messy_tnudge_si,     ONLY:      tnudge_initialize
  USE messy_trexp_si,      ONLY:       trexp_initialize
  USE messy_tropop_si,     ONLY:      tropop_initialize
  USE messy_trsync_si,     ONLY:      trsync_initialize
  USE messy_ubcnox_e5,     ONLY:      ubcnox_initialize
  USE messy_vahr_e5,       ONLY:        vahr_initialize
  USE messy_vaxtra_si,     ONLY:      vaxtra_initialize
  USE messy_vertex_si,     ONLY:      vertex_initialize
  USE messy_viso_si,       ONLY:        viso_initialize
  USE messy_mmd2way_si,    ONLY:     mmd2way_initialize

  IMPLICIT NONE

  CALL main_control_set_context(MEP_INITIALIZE)
  ! do NOT change the following order of main_* subroutines
  CALL main_timer_initialize
  CALL main_grid_initialize
  CALL main_data_initialize
  CALL main_channel_initialize
  CALL main_tracer_initialize ! namelist(s)
  CALL main_tools_initialize
  CALL main_import_initialize
#ifdef MESSYTENDENCY
  CALL main_tendency_initialize(1)
  CALL main_tendency_initialize(2)
#endif

  ! special requirements:
  ! - attila_initialize must be called first
  !   (define number of lagrangian cells)
  ! - ions must be called before gec
  ! - (the rest is in alphabetical order)
  IF (USE_ATTILA)      CALL      attila_initialize
  IF (USE_ACCF)        CALL        accf_initialize
  IF (USE_AEROPT)      CALL      aeropt_initialize
  IF (USE_AIRSEA)      CALL      airsea_initialize
  IF (USE_AIRTRAC)     CALL     airtrac_initialize
  IF (USE_AIRTRAF)     CALL     airtraf_initialize
  IF (USE_AVEOUT)      CALL      aveout_initialize
  IF (USE_BIOBURN)     CALL     bioburn_initialize
  IF (USE_BUFLY)       CALL       bufly_initialize
  IF (USE_CAT)         CALL         cat_initialize
  IF (USE_CH4)         CALL         ch4_initialize
  IF (USE_CLAMS)       CALL       clams_initialize
  IF (USE_CLAMSTRAJ)   CALL   clamstraj_initialize
  IF (USE_CLAMSCHEM)   CALL   clamschem_initialize
  IF (USE_CLAMSCHEME5) CALL clamscheme5_initialize
  IF (USE_CLAMSMIX)    CALL    clamsmix_initialize
  IF (USE_CLAMSBMIX)   CALL   clamsbmix_initialize
  IF (USE_CLAMSSEDI)   CALL   clamssedi_initialize
  IF (USE_CLAMSRDFRC)  CALL  clamsrdfrc_initialize
  IF (USE_CLAMSCIRRUS) CALL clamscirrus_initialize
  IF (USE_CLAMSTRACER) CALL clamstracer_initialize
  IF (USE_CHEMGLUE)    CALL    chemglue_initialize
  IF (USE_CLOUD)       CALL       cloud_initialize
  IF (USE_CLOUDOPT)    CALL    cloudopt_initialize
  IF (USE_CONTRAIL)    CALL    contrail_initialize
  IF (USE_CONVECT)     CALL     convect_initialize
  IF (USE_CRM)         CALL         crm_initialize
  IF (USE_CVTRANS)     CALL     cvtrans_initialize
  IF (USE_D14CO)       CALL       d14co_initialize
  IF (USE_DDEP)        CALL        ddep_initialize
  IF (USE_DISSOC)      CALL      dissoc_initialize
  IF (USE_DIUMOD)      CALL      diumod_initialize
  IF (USE_DOMINT)      CALL      domint_initialize
  IF (USE_DRADON)      CALL      dradon_initialize
  IF (USE_E4CHEM)      CALL      e4chem_initialize
  IF (USE_E5VDIFF)     CALL     e5vdiff_initialize
  IF (USE_EDITH)       CALL       edith_initialize
  IF (USE_EVER)        CALL        ever_initialize
  IF (USE_IONS)        CALL        ions_initialize
  IF (USE_GEC)         CALL         gec_initialize
  IF (USE_GMXE)        CALL        gmxe_initialize
  IF (USE_GWAVE)       CALL       gwave_initialize
  IF (USE_H2O)         CALL         h2o_initialize
  IF (USE_H2OEMIS)     CALL     h2oemis_initialize
  IF (USE_H2OISO)      CALL      h2oiso_initialize
  IF (USE_HD)          CALL          hd_initialize
  IF (USE_ISOPCOR)     CALL     isopcor_initialize
  IF (USE_CLOUDJ)      CALL      cloudj_initialize
  IF (USE_JVAL)        CALL        jval_initialize
  IF (USE_JVST)        CALL        jvst_initialize
  IF (USE_LGGP)        CALL        lggp_initialize
  IF (USE_LGTMIX)      CALL      lgtmix_initialize
  IF (USE_LGVFLUX)     CALL     lgvflux_initialize
  IF (USE_LNOX)        CALL        lnox_initialize
  IF (USE_M7)          CALL          m7_initialize
  IF (USE_MADE)        CALL        made_initialize
  IF (USE_MADE3)       CALL       made3_initialize
  IF (USE_MECCA)       CALL       mecca_initialize
  IF (USE_MEGAN)       CALL       megan_initialize
  IF (USE_MESOENERGY)  CALL  mesoenergy_initialize
  IF (USE_MLOCEAN)     CALL     mlocean_initialize
  IF (USE_MMFORCE)     CALL     mmforce_initialize
  IF (USE_MPIOM)       CALL       mpiom_initialize
  IF (USE_NAN)         CALL         nan_initialize
  IF (USE_A2O)         CALL         a2o_initialize
  IF (USE_OASIS3MCT)   CALL   oasis3mct_initialize
  IF (USE_OTPHYSC)     CALL     otphysc_initialize
  IF (USE_HAMOCC)      CALL      hamocc_initialize
  IF (USE_MSBM)        CALL        msbm_initialize
  IF (USE_MTSKIP)      CALL      mtskip_initialize
  IF (USE_O3ORIG)      CALL      o3orig_initialize
  IF (USE_OFFEMIS)     CALL     offemis_initialize
  IF (USE_ONEMIS)      CALL      onemis_initialize
  IF (USE_ORACLE)      CALL      oracle_initialize
  IF (USE_ORBIT)       CALL       orbit_initialize
  IF (USE_OROGW)       CALL       orogw_initialize
  IF (USE_PLUMEGAS)    CALL    plumegas_initialize
  IF (USE_PTRAC)       CALL       ptrac_initialize
  IF (USE_PTRACINI)    CALL    ptracini_initialize
  IF (USE_QBO)         CALL         qbo_initialize
  IF (USE_RAD)         CALL         rad_initialize
  IF (USE_RELAX)       CALL       relax_initialize
  IF (USE_RNDTEST)     CALL     rndtest_initialize
  IF (USE_SATSIMS)     CALL     satsims_initialize
  IF (USE_SCALC)       CALL       scalc_initialize
  IF (USE_SCAV)        CALL        scav_initialize
  IF (USE_SCOUT)       CALL       scout_initialize
  IF (USE_SEDI)        CALL        sedi_initialize
  IF (USE_S4D)         CALL         s4d_initialize
  IF (USE_SF6)         CALL         sf6_initialize
  IF (USE_SORBIT)      CALL      sorbit_initialize
  IF (USE_SPE)         CALL         spe_initialize
  IF (USE_SPACENOX)    CALL    spacenox_initialize
  IF (USE_SVOC)        CALL        svoc_initialize
  IF (USE_SURFACE)     CALL     surface_initialize
  IF (USE_TAGGING)     CALL     tagging_initialize
  IF (USE_TBUDGET)     CALL     tbudget_initialize
  IF (USE_TPULSE)      CALL      tpulse_initialize
  IF (USE_TIMEPOS)     CALL     timepos_initialize
  IF (USE_TNUDGE)      CALL      tnudge_initialize
  IF (USE_TREXP)       CALL       trexp_initialize
  IF (USE_TROPOP)      CALL      tropop_initialize
  IF (USE_TRSYNC)      CALL      trsync_initialize
  IF (USE_UBCNOX)      CALL      ubcnox_initialize
  IF (USE_VAHR)        CALL        vahr_initialize
  IF (USE_VAXTRA)      CALL      vaxtra_initialize
  IF (USE_VERTEX)      CALL      vertex_initialize
  IF (USE_VISO)        CALL        viso_initialize
  IF (USE_MMD2WAY)     CALL     mmd2way_initialize

  CALL main_control_unset_context

END SUBROUTINE messy_initialize

!==============================================================================

SUBROUTINE messy_new_tracer

  ! request tracers (called by initialize)

  USE messy_main_switch    ! ONLY: USE_*
  USE messy_main_control   ! ONLY: entrypoint, subentry
  USE messy_main_tracer_bi,  ONLY: main_tracer_new_tracer

  USE messy_airtrac_e5,     ONLY:     airtrac_new_tracer
  USE messy_ch4_si,         ONLY:         ch4_new_tracer
  USE messy_clams_si,       ONLY:       clams_new_tracer
  USE messy_clamscheme5_si, ONLY: clamscheme5_new_tracer
  USE messy_cloud_si,       ONLY:       cloud_new_tracer
  USE messy_d14co_si,       ONLY:       d14co_new_tracer
  USE messy_dradon_si,      ONLY:      dradon_new_tracer
  USE messy_e4chem_si,      ONLY:      e4chem_new_tracer
  USE messy_gmxe_si,        ONLY:        gmxe_new_tracer
  USE messy_hamocc_e5,      ONLY:      hamocc_new_tracer
  USE messy_h2o_e5,         ONLY:         h2o_new_tracer
  USE messy_h2oiso_e5,      ONLY:      h2oiso_new_tracer
  USE messy_m7_si,          ONLY:          m7_new_tracer
  USE messy_made_si,        ONLY:        made_new_tracer
  USE messy_made3_si,       ONLY:       made3_new_tracer
  USE messy_mecca_si,       ONLY:       mecca_new_tracer
  USE messy_msbm_si,        ONLY:        msbm_new_tracer
  USE messy_o3orig_si,      ONLY:      o3orig_new_tracer
  USE messy_offemis_si,     ONLY:     offemis_new_tracer
  USE messy_oracle_si,      ONLY:      oracle_new_tracer
  USE messy_plumegas_si,    ONLY:    plumegas_new_tracer
  USE messy_ptrac_si,       ONLY:       ptrac_new_tracer
  USE messy_scav_si,        ONLY:        scav_new_tracer
  USE messy_svoc_si,        ONLY:        svoc_new_tracer
  USE messy_sf6_e5,         ONLY:         sf6_new_tracer
  USE messy_tagging_si,     ONLY:     tagging_new_tracer
  USE messy_tbudget_si,     ONLY:     tbudget_new_tracer
  USE messy_tpulse_si,      ONLY:      tpulse_new_tracer
  USE messy_trexp_si,       ONLY:       trexp_new_tracer
  USE messy_ubcnox_e5,      ONLY:      ubcnox_new_tracer

  IMPLICIT NONE

  CALL main_control_set_context(MEP_NEW_TRACER)

  CALL main_tracer_new_tracer(1)  ! define tracer set

  ! special requirements:
  ! - mecca_new_tracer must be called first!
  ! - h2o_new_tracer should be called before  gmXe_new_tracer
  ! - oracle_new_tracer should be called after gmxe_new_tracer
  ! - cloud_new_tracer must be called after h2oiso_new_tracer
  ! - (the rest is in alphabetical order)
  IF (USE_MECCA)       CALL       mecca_new_tracer
  IF (USE_AIRTRAC)     CALL     airtrac_new_tracer
  IF (USE_CH4)         CALL         ch4_new_tracer
  IF (USE_CLAMS)       CALL       clams_new_tracer
  IF (USE_CLAMSCHEME5) CALL clamscheme5_new_tracer
  IF (USE_D14CO)       CALL       d14co_new_tracer
  IF (USE_DRADON)      CALL      dradon_new_tracer
  IF (USE_E4CHEM)      CALL      e4chem_new_tracer
  IF (USE_HAMOCC)      CALL      hamocc_new_tracer
  IF (USE_H2O)         CALL         h2o_new_tracer
  IF (USE_H2OISO)      CALL      h2oiso_new_tracer
  IF (USE_CLOUD)       CALL       cloud_new_tracer
  IF (USE_GMXE)        CALL        gmxe_new_tracer
  IF (USE_M7)          CALL          m7_new_tracer
  IF (USE_MADE)        CALL        made_new_tracer
  IF (USE_MADE3)       CALL       made3_new_tracer
  IF (USE_MSBM)        CALL        msbm_new_tracer
  IF (USE_O3ORIG)      CALL      o3orig_new_tracer
  IF (USE_OFFEMIS)     CALL     offemis_new_tracer
  IF (USE_ORACLE)      CALL      oracle_new_tracer
  IF (USE_PLUMEGAS)    CALL    plumegas_new_tracer
  IF (USE_PTRAC)       CALL       ptrac_new_tracer
  IF (USE_SCAV)        CALL        scav_new_tracer
  IF (USE_SVOC)        CALL        svoc_new_tracer
  IF (USE_SF6)         CALL         sf6_new_tracer
  IF (USE_TAGGING)     CALL     tagging_new_tracer
  IF (USE_TBUDGET)     CALL     tbudget_new_tracer
  IF (USE_TPULSE)      CALL      tpulse_new_tracer
  IF (USE_TREXP)       CALL       trexp_new_tracer
  IF (USE_UBCNOX)      CALL      ubcnox_new_tracer

  CALL main_tracer_new_tracer(2) ! define tracer families
  CALL main_tracer_new_tracer(3) ! diagnostic output

  CALL main_control_unset_context

END SUBROUTINE messy_new_tracer

!==============================================================================

SUBROUTINE messy_init_memory

  ! allocate memory and define output channels
  ! (called by init_memory in mo_memory_streams)

  USE messy_main_switch    ! ONLY: USE_*
  USE messy_main_control   ! ONLY: entrypoint, subentry
  USE messy_main_qtimer_bi,  ONLY: main_qtimer_init_memory
  USE messy_main_data_bi,    ONLY: main_data_init_memory
  USE messy_main_tracer_bi,  ONLY: main_tracer_init_memory
  USE messy_main_channel_bi, ONLY: main_channel_init_memory
  USE messy_main_grid_bi,    ONLY: main_grid_init_memory
  USE messy_main_import_bi,  ONLY: main_import_init_memory

  USE messy_a2o_e5,         ONLY:         a2o_init_memory
  USE messy_accf_e5,        ONLY:        accf_init_memory
  USE messy_aeropt_si,      ONLY:      aeropt_init_memory
  USE messy_airsea_si,      ONLY:      airsea_init_memory
  USE messy_airtrac_e5,     ONLY:     airtrac_init_memory
  USE messy_airtraf_e5,     ONLY:     airtraf_init_memory
  USE messy_attila_e5,      ONLY:      attila_init_memory
  USE messy_bioburn_si,     ONLY:     bioburn_init_memory
  USE messy_cat_e5,         ONLY:         cat_init_memory
  USE messy_ch4_si,         ONLY:         ch4_init_memory
  USE messy_clams_si,       ONLY:       clams_init_memory
  USE messy_clamstraj_si,   ONLY:   clamstraj_init_memory
  USE messy_clamschem_si,   ONLY:   clamschem_init_memory
  USE messy_clamscheme5_si, ONLY: clamscheme5_init_memory
  USE messy_clamsmix_si,    ONLY:    clamsmix_init_memory
  USE messy_clamsbmix_si,   ONLY:   clamsbmix_init_memory
  USE messy_clamssedi_si,   ONLY:   clamssedi_init_memory
  USE messy_clamsrdfrc_si,  ONLY:  clamsrdfrc_init_memory
  USE messy_clamscirrus_si, ONLY: clamscirrus_init_memory
  USE messy_clamstracer_si, ONLY: clamstracer_init_memory
  USE messy_chemglue_si,    ONLY:    chemglue_init_memory
  USE messy_cloud_si,       ONLY:       cloud_init_memory
  USE messy_cloudj_si,      ONLY:      cloudj_init_memory
  USE messy_cloudopt_si,    ONLY:    cloudopt_init_memory
  USE messy_contrail_si,    ONLY:    contrail_init_memory
  USE messy_convect_si,     ONLY:     convect_init_memory
  USE messy_crm_e5,         ONLY:         crm_init_memory
  USE messy_cvtrans_si,     ONLY:     cvtrans_init_memory
  USE messy_d14co_si,       ONLY:       d14co_init_memory
  USE messy_ddep_si,        ONLY:        ddep_init_memory
  USE messy_dissoc_si,      ONLY:      dissoc_init_memory
  USE messy_diumod_si,      ONLY:      diumod_init_memory
  USE messy_domint_si,      ONLY:      domint_init_memory
  USE messy_dradon_si,      ONLY:      dradon_init_memory
  USE messy_e4chem_si,      ONLY:      e4chem_init_memory
  USE messy_e5vdiff_si,     ONLY:     e5vdiff_init_memory
  USE messy_ec2cosmo_e5,    ONLY:    ec2cosmo_init_memory
  USE messy_edith_e5,       ONLY:       edith_init_memory
  USE messy_ever_si,        ONLY:        ever_init_memory
  USE messy_gec_si,         ONLY:         gec_init_memory
  USE messy_gmxe_si,        ONLY:        gmxe_init_memory
  USE messy_gwave_si,       ONLY:       gwave_init_memory
  USE messy_h2o_e5,         ONLY:         h2o_init_memory
  USE messy_h2oiso_e5,      ONLY:      h2oiso_init_memory
  USE messy_hd_e5,          ONLY:          hd_init_memory
  USE messy_hamocc_e5,      ONLY:      hamocc_init_memory
  USE messy_isopcor_si,     ONLY:     isopcor_init_memory
  USE messy_ions_si,        ONLY:        ions_init_memory
  USE messy_jval_si,        ONLY:        jval_init_memory
  USE messy_jvst_si,        ONLY:        jvst_init_memory
  USE messy_lgtmix_e5,      ONLY:      lgtmix_init_memory
  USE messy_lnox_si,        ONLY:        lnox_init_memory
  USE messy_m7_si,          ONLY:          m7_init_memory
  USE messy_made_si,        ONLY:        made_init_memory
  USE messy_made3_si,       ONLY:       made3_init_memory
  USE messy_mecca_si,       ONLY:       mecca_init_memory
  USE messy_megan_si,       ONLY:       megan_init_memory
  USE messy_mesoenergy_e5,  ONLY:  mesoenergy_init_memory
  USE messy_mlocean_e5,     ONLY:     mlocean_init_memory
  USE messy_mmforce_e5,     ONLY:     mmforce_init_memory
  USE messy_mmd2way_si,     ONLY:     mmd2way_init_memory
  USE messy_mpiom_e5,       ONLY:       mpiom_init_memory
  USE messy_msbm_si,        ONLY:        msbm_init_memory
  USE messy_mtskip_si,      ONLY:      mtskip_init_memory
  USE messy_nan_si,         ONLY:         nan_init_memory
  USE messy_o3orig_si,      ONLY:      o3orig_init_memory
  USE messy_onemis_si,      ONLY:      onemis_init_memory
  USE messy_oracle_si,      ONLY:      oracle_init_memory
  USE messy_orbit_si,       ONLY:       orbit_init_memory
  USE messy_orogw_si,       ONLY:       orogw_init_memory
  USE messy_otphysc_e5,     ONLY:     otphysc_init_memory
  USE messy_plumegas_si,    ONLY:    plumegas_init_memory
  USE messy_ptrac_si,       ONLY:       ptrac_init_memory
  USE messy_ptracini_si,    ONLY:    ptracini_init_memory
  USE messy_qbo_si,         ONLY:         qbo_init_memory
  USE messy_rad_si,         ONLY:         rad_init_memory
  USE messy_relax_e5,       ONLY:       relax_init_memory
  USE messy_rndtest_si,     ONLY:     rndtest_init_memory
  USE messy_s4d_si,         ONLY:         s4d_init_memory
  USE messy_satsims_e5,     ONLY:     satsims_init_memory
  USE messy_scav_si,        ONLY:        scav_init_memory
  USE messy_scout_si,       ONLY:       scout_init_memory
  USE messy_sedi_si,        ONLY:        sedi_init_memory
  USE messy_sf6_e5,         ONLY:         sf6_init_memory
  USE messy_spe_e5,         ONLY:         spe_init_memory
  USE messy_spacenox_e5,    ONLY:    spacenox_init_memory
  USE messy_surface_si,     ONLY:     surface_init_memory
  USE messy_svoc_si,        ONLY:        svoc_init_memory
  USE messy_tagging_si,     ONLY:     tagging_init_memory
  USE messy_tbudget_si,     ONLY:     tbudget_init_memory
  USE messy_tpulse_si,      ONLY:      tpulse_init_memory
  USE messy_tnudge_si,      ONLY:      tnudge_init_memory
  USE messy_trexp_si,       ONLY:       trexp_init_memory
  USE messy_tropop_si,      ONLY:      tropop_init_memory
  USE messy_trsync_si,      ONLY:      trsync_init_memory
  USE messy_ubcnox_e5,      ONLY:      ubcnox_init_memory
  USE messy_vaxtra_si,      ONLY:      vaxtra_init_memory
  USE messy_vertex_si,      ONLY:      vertex_init_memory
  USE messy_oasis3mct_si,   ONLY: oasis3mct_init_memory, oasis3mct_init_coupling
  USE messy_scalc_si,       ONLY:       scalc_init_memory

  IMPLICIT NONE

  CALL main_control_set_context(MEP_INIT_MEMORY)
  CALL main_qtimer_init_memory
  ! associate ECHAM5 streams to MESSy data channels
  CALL main_channel_init_memory   ! setup tracer memory
  CALL main_tracer_init_memory(1) ! create data channel BML <-> SMIL
  CALL main_grid_init_memory(1)   ! CALL GRID_DEF before DATA
  CALL main_data_init_memory
  CALL main_grid_init_memory(2)   !
  CALL main_import_init_memory

  ! special requirements:
  ! - HAMOCC after MPIOM
  ! - IONS must be called befoer GEC
  ! - (in alphabetical order)
  IF (USE_A2O)         CALL        a2o_init_memory
  IF (USE_ACCF)        CALL       accf_init_memory
  IF (USE_AIRTRAC)     CALL    airtrac_init_memory
  IF (USE_ATTILA)      CALL     attila_init_memory
  IF (USE_AEROPT)      CALL     aeropt_init_memory
  IF (USE_AIRSEA)      CALL     airsea_init_memory
  IF (USE_AIRTRAF)     CALL    airtraf_init_memory
  IF (USE_BIOBURN)     CALL    bioburn_init_memory
  IF (USE_CAT)         CALL        cat_init_memory
  IF (USE_CH4)         CALL        ch4_init_memory
  IF (USE_CLAMS)       CALL       clams_init_memory
  IF (USE_CLAMSTRAJ)   CALL   clamstraj_init_memory
  IF (USE_CLAMSCHEM)   CALL   clamschem_init_memory
  IF (USE_CLAMSCHEME5) CALL clamscheme5_init_memory
  IF (USE_CLAMSSEDI)   CALL   clamssedi_init_memory
  IF (USE_CLAMSMIX)    CALL    clamsmix_init_memory
  IF (USE_CLAMSBMIX)   CALL   clamsbmix_init_memory
  IF (USE_CLAMSCIRRUS) CALL clamscirrus_init_memory
  IF (USE_CLAMSRDFRC)  CALL  clamsrdfrc_init_memory
  IF (USE_CLAMSTRACER) CALL clamstracer_init_memory
  IF (USE_CHEMGLUE)    CALL   chemglue_init_memory
  IF (USE_CLOUD)       CALL      cloud_init_memory
  IF (USE_CLOUDOPT)    CALL   cloudopt_init_memory
  IF (USE_CONTRAIL)    CALL   contrail_init_memory
  IF (USE_CONVECT)     CALL    convect_init_memory
  IF (USE_CRM)         CALL        crm_init_memory
  IF (USE_CVTRANS)     CALL    cvtrans_init_memory
  IF (USE_D14CO)       CALL      d14co_init_memory
  IF (USE_DDEP)        CALL       ddep_init_memory
  IF (USE_DISSOC)      CALL     dissoc_init_memory
  IF (USE_DIUMOD)      CALL     diumod_init_memory
  IF (USE_DOMINT)      CALL     domint_init_memory
  IF (USE_DRADON)      CALL     dradon_init_memory
  IF (USE_E4CHEM)      CALL     e4chem_init_memory
  IF (USE_E5VDIFF)     CALL    e5vdiff_init_memory
  IF (USE_EC2COSMO)    CALL   ec2cosmo_init_memory
  IF (USE_EDITH)       CALL      edith_init_memory
  IF (USE_EVER)        CALL       ever_init_memory
  IF (USE_IONS)        CALL       ions_init_memory
  IF (USE_GEC)         CALL        gec_init_memory
  IF (USE_GMXE)        CALL       gmxe_init_memory
  IF (USE_GWAVE)       CALL      gwave_init_memory
  IF (USE_H2O)         CALL        h2o_init_memory
  IF (USE_H2OISO)      CALL     h2oiso_init_memory
  IF (USE_HD)          CALL         hd_init_memory
  IF (USE_ISOPCOR)     CALL    isopcor_init_memory
  IF (USE_CLOUDJ)      CALL     cloudj_init_memory
  IF (USE_JVAL)        CALL       jval_init_memory
  IF (USE_JVST)        CALL       jvst_init_memory
  IF (USE_LGTMIX)      CALL     lgtmix_init_memory
  IF (USE_LNOX)        CALL       lnox_init_memory
  IF (USE_M7)          CALL         m7_init_memory
  IF (USE_MADE)        CALL       made_init_memory
  IF (USE_MADE3)       CALL      made3_init_memory
  IF (USE_MECCA)       CALL      mecca_init_memory
  IF (USE_MEGAN)       CALL      megan_init_memory
  IF (USE_MESOENERGY)  CALL mesoenergy_init_memory
  IF (USE_MLOCEAN)     CALL    mlocean_init_memory
  IF (USE_MMD2WAY)     CALL    mmd2way_init_memory(1)
  IF (USE_MMFORCE)     CALL    mmforce_init_memory
  IF (USE_MPIOM)       CALL      mpiom_init_memory
  IF (USE_MTSKIP)      CALL     mtskip_init_memory
  IF (USE_HAMOCC)      CALL     hamocc_init_memory
  IF (USE_MSBM)        CALL       msbm_init_memory
  IF (USE_NAN)         CALL        nan_init_memory
  IF (USE_O3ORIG)      CALL     o3orig_init_memory
  IF (USE_ONEMIS)      CALL     onemis_init_memory
  IF (USE_ORACLE)      CALL     oracle_init_memory
  IF (USE_ORBIT)       CALL      orbit_init_memory
  IF (USE_OROGW)       CALL      orogw_init_memory
  IF (USE_OTPHYSC)     CALL    otphysc_init_memory
  IF (USE_PLUMEGAS)    CALL   plumegas_init_memory
  IF (USE_PTRAC)       CALL      ptrac_init_memory
  IF (USE_PTRACINI)    CALL   ptracini_init_memory
  IF (USE_QBO)         CALL        qbo_init_memory
  IF (USE_RAD)         CALL        rad_init_memory
  IF (USE_RELAX)       CALL      relax_init_memory
  IF (USE_RNDTEST)     CALL    rndtest_init_memory
  IF (USE_S4D)         CALL        s4d_init_memory
  IF (USE_SATSIMS)     CALL    satsims_init_memory
  IF (USE_SCAV)        CALL       scav_init_memory
  IF (USE_SCOUT)       CALL      scout_init_memory
  IF (USE_SEDI)        CALL       sedi_init_memory
  IF (USE_SF6)         CALL        sf6_init_memory
  IF (USE_SPE)         CALL        spe_init_memory
  IF (USE_SPACENOX)    CALL   spacenox_init_memory
  IF (USE_SURFACE)     CALL    surface_init_memory
  IF (USE_SVOC)        CALL       svoc_init_memory
  IF (USE_TAGGING)     CALL    tagging_init_memory
  IF (USE_TBUDGET)     CALL    tbudget_init_memory
  IF (USE_TPULSE)      CALL     tpulse_init_memory
  IF (USE_TNUDGE)      CALL     tnudge_init_memory
  IF (USE_TREXP)       CALL      trexp_init_memory
  IF (USE_TROPOP)      CALL     tropop_init_memory
  IF (USE_TRSYNC)      CALL     trsync_init_memory
  IF (USE_UBCNOX)      CALL     ubcnox_init_memory
  IF (USE_VAXTRA)      CALL     vaxtra_init_memory
  IF (USE_VERTEX)      CALL     vertex_init_memory

  ! - associate tracer memory to MESSy channel(s)
  ! - setting meta information of family-members to fraction
  !   (for advection initialization)
  CALL main_tracer_init_memory(2)

!!#D oasis3mct +
  IF (USE_OASIS3MCT) THEN
     CALL     oasis3mct_init_memory
     IF (USE_SCALC)   CALL  scalc_init_memory(l_cont=.TRUE.)
    CALL   oasis3mct_init_coupling
  END IF
!!#D oasis3mct -

  CALL main_control_unset_context

END SUBROUTINE messy_init_memory

!==============================================================================

SUBROUTINE messy_init_coupling

  ! (will be called by control.f90)

  USE messy_main_switch   ! ONLY: USE_*
  USE messy_main_control  ! ONLY: entrypoint, subentry
  USE messy_main_channel_bi, ONLY: main_channel_init_coupling
  USE messy_main_data_bi,    ONLY: main_data_init_coupling
  USE messy_main_tracer_bi,  ONLY: main_tracer_init_coupling
#ifdef MESSYTENDENCY
  USE messy_main_tendency_bi, ONLY: main_tendency_init_coupling
#endif
  USE messy_main_rnd_bi,     ONLY: main_rnd_init_coupling

  USE messy_a2o_e5,          ONLY:         a2o_init_coupling
  USE messy_accf_e5,         ONLY:        accf_init_coupling
  USE messy_aeropt_si,       ONLY:      aeropt_init_coupling
  USE messy_airsea_si,       ONLY:      airsea_init_coupling
  USE messy_airtrac_e5,      ONLY:     airtrac_init_coupling
  USE messy_airtraf_e5,      ONLY:     airtraf_init_coupling
  USE messy_attila_e5,       ONLY:      attila_init_coupling
  USE messy_aveout_si,       ONLY:      aveout_init_coupling
  USE messy_bioburn_si,      ONLY:     bioburn_init_coupling
  USE messy_bufly_si,        ONLY:       bufly_init_coupling
  USE messy_cat_e5,          ONLY:         cat_init_coupling
  USE messy_ch4_si,          ONLY:         ch4_init_coupling
  !USE messy_chemglue_si,    ONLY:    chemglue_init_coupling
  USE messy_clams_si,        ONLY:       clams_init_coupling
  USE messy_clamstraj_si,    ONLY:   clamstraj_init_coupling
  USE messy_clamschem_si,    ONLY:   clamschem_init_coupling
  USE messy_clamscheme5_si,  ONLY: clamscheme5_init_coupling
  USE messy_clamssedi_si,    ONLY:   clamssedi_init_coupling
  USE messy_clamsmix_si,     ONLY:    clamsmix_init_coupling
  USE messy_clamsbmix_si,    ONLY:   clamsbmix_init_coupling
  USE messy_clamscirrus_si,  ONLY: clamscirrus_init_coupling
  USE messy_clamsrdfrc_si,   ONLY:  clamsrdfrc_init_coupling
  USE messy_clamstracer_si,  ONLY: clamstracer_init_coupling
  USE messy_cloud_si,        ONLY:       cloud_init_coupling
  USE messy_cloudj_si,       ONLY:      cloudj_init_coupling
  USE messy_cloudopt_si,     ONLY:    cloudopt_init_coupling
  USE messy_convect_si,      ONLY:     convect_init_coupling
  USE messy_crm_e5,          ONLY:         crm_init_coupling
  USE messy_contrail_si,     ONLY:    contrail_init_coupling
  USE messy_cvtrans_si,      ONLY:     cvtrans_init_coupling
  USE messy_d14co_si,        ONLY:       d14co_init_coupling
  USE messy_ddep_si,         ONLY:        ddep_init_coupling
  USE messy_dissoc_si,       ONLY:      dissoc_init_coupling
  USE messy_diumod_si,       ONLY:      diumod_init_coupling
  USE messy_domint_si,       ONLY:      domint_init_coupling
  USE messy_dradon_si,       ONLY:      dradon_init_coupling
  USE messy_e4chem_si,       ONLY:      e4chem_init_coupling
  USE messy_e5vdiff_si,      ONLY:     e5vdiff_init_coupling
  USE messy_ec2cosmo_e5,     ONLY:    ec2cosmo_init_coupling
  USE messy_edith_e5,        ONLY:       edith_init_coupling
  USE messy_ever_si,         ONLY:        ever_init_coupling
  USE messy_gec_si,          ONLY:         gec_init_coupling
  USE messy_gmxe_si,         ONLY:        gmxe_init_coupling
  USE messy_hd_e5,           ONLY:          hd_init_coupling
  USE messy_h2o_e5,          ONLY:         h2o_init_coupling
  USE messy_h2oemis_si,      ONLY:     h2oemis_init_coupling
  USE messy_h2oiso_e5,       ONLY:      h2oiso_init_coupling
  USE messy_hamocc_e5,       ONLY:      hamocc_init_coupling
  USE messy_isopcor_si,      ONLY:     isopcor_init_coupling
  USE messy_ions_si,         ONLY:        ions_init_coupling
  USE messy_jval_si,         ONLY:        jval_init_coupling
  USE messy_jvst_si,         ONLY:        jvst_init_coupling
  USE messy_lggp_e5,         ONLY:        lggp_init_coupling
  USE messy_lgtmix_e5,       ONLY:      lgtmix_init_coupling
  USE messy_lgvflux_e5,      ONLY:     lgvflux_init_coupling
  USE messy_lnox_si,         ONLY:        lnox_init_coupling
  USE messy_m7_si,           ONLY:          m7_init_coupling
  USE messy_made_si,         ONLY:        made_init_coupling
  USE messy_made3_si,        ONLY:       made3_init_coupling
  USE messy_mecca_si,        ONLY:       mecca_init_coupling
  USE messy_megan_si,        ONLY:       megan_init_coupling
  USE messy_mesoenergy_e5,   ONLY:  mesoenergy_init_coupling
  USE messy_mlocean_e5,      ONLY:     mlocean_init_coupling
  USE messy_mmforce_e5,      ONLY:     mmforce_init_coupling
  USE messy_mpiom_e5,        ONLY:       mpiom_init_coupling
  USE messy_msbm_si,         ONLY:        msbm_init_coupling
  USE messy_mtskip_si,       ONLY:      mtskip_init_coupling
  USE messy_nan_si,          ONLY:         nan_init_coupling
  USE messy_o3orig_si,       ONLY:      o3orig_init_coupling
  USE messy_offemis_si,      ONLY:     offemis_init_coupling
  USE messy_onemis_si,       ONLY:      onemis_init_coupling
  USE messy_oracle_si,       ONLY:      oracle_init_coupling
  USE messy_orbit_si,        ONLY:       orbit_init_coupling
  USE messy_orogw_si,        ONLY:       orogw_init_coupling
  USE messy_otphysc_e5,      ONLY:     otphysc_init_coupling
  USE messy_plumegas_si,     ONLY:    plumegas_init_coupling
  USE messy_ptracini_si,     ONLY:    ptracini_init_coupling
  USE messy_qbo_si,          ONLY:         qbo_init_coupling
  USE messy_rad_si,          ONLY:         rad_init_coupling
  USE messy_relax_e5,        ONLY:       relax_init_coupling
  USE messy_s4d_si,          ONLY:         s4d_init_coupling
  USE messy_sorbit_si,       ONLY:      sorbit_init_coupling
  USE messy_satsims_e5,      ONLY:     satsims_init_coupling
  USE messy_scalc_si,        ONLY:       scalc_init_coupling
  USE messy_scav_si,         ONLY:        scav_init_coupling
  USE messy_scout_si,        ONLY:       scout_init_coupling
  USE messy_sedi_si,         ONLY:        sedi_init_coupling
  USE messy_sf6_e5,          ONLY:         sf6_init_coupling
  USE messy_spe_e5,          ONLY:         spe_init_coupling
  USE messy_spacenox_e5,     ONLY:    spacenox_init_coupling
  USE messy_surface_si,      ONLY:     surface_init_coupling
  USE messy_svoc_si,         ONLY:        svoc_init_coupling
  USE messy_tagging_si,      ONLY:     tagging_init_coupling
  USE messy_tbudget_si,      ONLY:     tbudget_init_coupling
  USE messy_timepos_e5,      ONLY:     timepos_init_coupling
  USE messy_tnudge_si,       ONLY:      tnudge_init_coupling
  USE messy_trexp_si,        ONLY:       trexp_init_coupling
  USE messy_trsync_si,       ONLY:      trsync_init_coupling
  USE messy_ubcnox_e5,       ONLY:      ubcnox_init_coupling
  USE messy_vahr_e5,         ONLY:        vahr_init_coupling
  USE messy_vaxtra_si,       ONLY:      vaxtra_init_coupling
  USE messy_vertex_si,       ONLY:      vertex_init_coupling
  USE messy_viso_si,         ONLY:        viso_init_coupling
  USE messy_mmd2way_si,      ONLY:     mmd2way_init_coupling

  IMPLICIT NONE

  CALL main_control_set_context(MEP_INIT_COUPLING)
  ! - resetting meta information of family-members (to tracers)
  !   (after advection initialisation)
  CALL main_tracer_init_coupling(1)

  CALL main_data_init_coupling

  ! - modify attributes according to ADD_ATT in CTRL namelist
  CALL main_channel_init_coupling(1)

  ! special requirements:
  ! - Some SMs (as indicated) need to define new channel objects in
  !   init_coupling. Therefore, the diagnostic SMs S4D, SORBIT, SCOUT, SCALC,
  !   VISO, LGGP, and LGVFLUX
  !   must be called at the end, since all other channel
  !   objects need to be present.
  ! - VISO must be called before LGVFLUX, because LGVFLUX requires the
  !   iso-surfaces defined in VISO
  ! - VISO (2nd call!) must be called after MMD2WAY
  ! - LGGP must be called before LGVFLUX
  ! - LGGP must be called twice (2nd time after  LGVFLUX) in order to be
  !   able to convert Lagrangian information from LGVFLUX
  ! - CH4 must be called after LGGP to enable usage of transformed fields,
  !   e.g., from IMPORT
  ! - GMXE must be called after offemis and onemis, since it uses objects
  !   from both submodels
  ! - ORBIT must be called before RAD, since it creates new objects for RAD
  ! - CONTRAIL must be called between calls to LGGP to enable GP-2-LG
  !   conversion of CLOUD parameters for CONTRAIL and LG-CONTRAIL results
  !   for LGGP output
  ! - CONTRAIL must be called after TREXP, because it uses objects from TREXP
  ! - TBUDGET needs to be called after SCALC, because summed elements are
  !   required
  ! - TNUDGE must be called after VISO for usage of iso-surface information
  ! - MTSKIP must be called in between two calls (with flag 1 and 2,
  !   respectively) of submodels equipped with mtskip (currently: MECCA)
  ! - (the rest is in alphabetical order)
  ! - DIUMOD needs to be called after ORBIT, since it might use the solar
  !   zenith angle shifted in time (created in ORBIT in init_coupling)
  ! - SCALC must be called before RAD to have sums of tracers for radiation
  !   calculation
  ! - TNUDGE must be called last since all channel objects, which are used by
  !   tnudge must be present
  ! - VAXTRA must be called last (except SCOUT), since all other channel objects must be
  !   present
  ! - VEG must be called before onemis for coupling of vegetation to emissions
  ! - CRM must be calld after e5vdiff/vertex
  ! - DOMINT must be called after all channels for variables for domain
  !   integration are created (like for SCOUT)
  ! - LNOX must be called after VISO (for ice flux iso-surface)

  IF (USE_A2O)         CALL         a2o_init_coupling
  IF (USE_ACCF)        CALL        accf_init_coupling
  IF (USE_ATTILA)      CALL      attila_init_coupling(1)
  IF (USE_AEROPT)      CALL      aeropt_init_coupling
  IF (USE_AIRSEA)      CALL      airsea_init_coupling
  IF (USE_AIRTRAC)     CALL     airtrac_init_coupling
  IF (USE_AIRTRAF)     CALL     airtraf_init_coupling
  IF (USE_AVEOUT)      CALL      aveout_init_coupling
  IF (USE_BIOBURN)     CALL     bioburn_init_coupling
  IF (USE_BUFLY)       CALL       bufly_init_coupling
  IF (USE_CAT)         CALL         cat_init_coupling
  !IF (USE_CHEMGLUE)   CALL    chemglue_init_coupling
  IF (USE_CLAMS)       CALL       clams_init_coupling(1)
  IF (USE_CLAMSTRAJ)   CALL   clamstraj_init_coupling
  IF (USE_CLAMSCHEM)   CALL   clamschem_init_coupling
  IF (USE_CLAMSCHEME5) CALL clamscheme5_init_coupling
  IF (USE_CLAMSSEDI)   CALL   clamssedi_init_coupling
  IF (USE_CLAMSMIX)    CALL    clamsmix_init_coupling
  IF (USE_CLAMSBMIX)   CALL   clamsbmix_init_coupling
  IF (USE_CLAMSCIRRUS) CALL clamscirrus_init_coupling
  IF (USE_CLAMSRDFRC)  CALL  clamsrdfrc_init_coupling
  IF (USE_CLAMSTRACER) CALL clamstracer_init_coupling
  IF (USE_CLOUD)       CALL       cloud_init_coupling
  IF (USE_CLOUDOPT)    CALL    cloudopt_init_coupling
  IF (USE_CONVECT)     CALL     convect_init_coupling
  IF (USE_CVTRANS)     CALL     cvtrans_init_coupling
  IF (USE_D14CO)       CALL       d14co_init_coupling
  IF (USE_DDEP)        CALL        ddep_init_coupling ! new channel objects
  IF (USE_DISSOC)      CALL      dissoc_init_coupling
  IF (USE_DRADON)      CALL      dradon_init_coupling
  IF (USE_E4CHEM)      CALL      e4chem_init_coupling
  IF (USE_E5VDIFF)     CALL     e5vdiff_init_coupling
  IF (USE_EC2COSMO)    CALL    ec2cosmo_init_coupling
  IF (USE_EDITH)       CALL       edith_init_coupling
  IF (USE_IONS)        CALL        ions_init_coupling
  IF (USE_GEC)         CALL         gec_init_coupling
  IF (USE_HD)          CALL          hd_init_coupling
  IF (USE_H2O)         CALL         h2o_init_coupling
  IF (USE_H2OEMIS)     CALL     h2oemis_init_coupling
  IF (USE_H2OISO)      CALL      h2oiso_init_coupling(1)
  IF (USE_ISOPCOR)     CALL     isopcor_init_coupling
  IF (USE_CLOUDJ)      CALL      cloudj_init_coupling
  IF (USE_JVAL)        CALL        jval_init_coupling
  IF (USE_JVST)        CALL        jvst_init_coupling
  IF (USE_LGTMIX)      CALL      lgtmix_init_coupling
  IF (USE_M7)          CALL          m7_init_coupling
  IF (USE_MADE)        CALL        made_init_coupling
  IF (USE_MADE3)       CALL       made3_init_coupling
  IF (USE_MECCA)       CALL       mecca_init_coupling(1)
  IF (USE_MEGAN)       CALL       megan_init_coupling
  IF (USE_MESOENERGY)  CALL  mesoenergy_init_coupling
  IF (USE_MLOCEAN)     CALL     mlocean_init_coupling
  IF (USE_MMFORCE)     CALL     mmforce_init_coupling
  IF (USE_OTPHYSC)     CALL     otphysc_init_coupling
  IF (USE_MPIOM)       CALL       mpiom_init_coupling
  IF (USE_HAMOCC)      CALL      hamocc_init_coupling
  IF (USE_MSBM)        CALL        msbm_init_coupling
  IF (USE_NAN)         CALL         nan_init_coupling
  IF (USE_O3ORIG)      CALL      o3orig_init_coupling
  IF (USE_SCALC)       CALL       scalc_init_coupling(l_cont=.TRUE.)
                                                ! qqq to be done BEFORE OFFEMIS
!                                               !     for FLAGGING of emissions
  IF (USE_OFFEMIS)     CALL     offemis_init_coupling ! new channel objects
  IF (USE_ONEMIS)      CALL      onemis_init_coupling
  IF (USE_ORACLE)      CALL      oracle_init_coupling
  IF (USE_ORBIT)       CALL       orbit_init_coupling ! new channel objects
  IF (USE_OROGW)       CALL       orogw_init_coupling
  IF (USE_DIUMOD)      CALL      diumod_init_coupling
  IF (USE_GMXE)        CALL        gmxe_init_coupling
  IF (USE_PLUMEGAS)    CALL    plumegas_init_coupling
  IF (USE_PTRACINI)    CALL    ptracini_init_coupling
  IF (USE_QBO)         CALL         qbo_init_coupling
!!$  IF (USE_RAD)      CALL         rad_init_coupling
  IF (USE_SATSIMS)     CALL     satsims_init_coupling
  IF (USE_SCAV)        CALL        scav_init_coupling
  IF (USE_SEDI)        CALL        sedi_init_coupling ! new channel objects
  IF (USE_SF6)         CALL         sf6_init_coupling
  IF (USE_SVOC)        CALL        svoc_init_coupling
  IF (USE_TIMEPOS)     CALL     timepos_init_coupling
!+  IF (USE_TREXP)       CALL       trexp_init_coupling ! new channel objects
  IF (USE_TRSYNC)      CALL      trsync_init_coupling
  IF (USE_SPE)         CALL         spe_init_coupling
  IF (USE_SPACENOX)    CALL    spacenox_init_coupling
  IF (USE_VAHR)        CALL        vahr_init_coupling
  IF (USE_TAGGING)     CALL     tagging_init_coupling
  IF (USE_VISO)        CALL        viso_init_coupling ! new channel objects
! IF (USE_TNUDGE)      CALL      tnudge_init_coupling
  IF (USE_LGGP)        CALL        lggp_init_coupling ! new channel objects (1)
  IF (USE_CH4)         CALL         ch4_init_coupling
  IF (USE_SCALC)       CALL       scalc_init_coupling ! new channel objects
  IF (USE_TREXP)       CALL       trexp_init_coupling ! new channel objects
  IF (USE_EVER)        CALL        ever_init_coupling ! new channel objects
  IF (USE_CONTRAIL)    CALL    contrail_init_coupling
  IF (USE_MMD2WAY)     CALL     mmd2way_init_coupling ! new channel objects
!+  IF (USE_SCALC)       CALL       scalc_init_coupling ! new channel objects
  IF (USE_RAD)         CALL         rad_init_coupling
  IF (USE_RELAX)       CALL       relax_init_coupling ! new channel objects
  IF (USE_SURFACE)     CALL     surface_init_coupling
  IF (USE_TBUDGET)     CALL     tbudget_init_coupling
  IF (USE_VERTEX)      CALL      vertex_init_coupling
  IF (USE_CRM)         CALL         crm_init_coupling
  IF (USE_LGVFLUX)     CALL     lgvflux_init_coupling ! new channel objects
  IF (USE_LGGP)        CALL        lggp_init_coupling ! new channel objects (2)
  IF (USE_S4D)         CALL         s4d_init_coupling ! new channel objects
  IF (USE_SORBIT)      CALL      sorbit_init_coupling ! new channel objects
  IF (USE_VISO)        CALL        viso_init_coupling ! new channel objects (2)
  IF (USE_LNOX)        CALL        lnox_init_coupling
  IF (USE_UBCNOX)      CALL      ubcnox_init_coupling
  IF (USE_MTSKIP)      CALL      mtskip_init_coupling
  IF (USE_MECCA)       CALL      mecca_init_coupling(2)
  IF (USE_TNUDGE)      CALL      tnudge_init_coupling
  IF (USE_VAXTRA)      CALL      vaxtra_init_coupling
  IF (USE_SCOUT)       CALL       scout_init_coupling ! new channel objects
  IF (USE_DOMINT)      CALL      domint_init_coupling

  ! ALLOCATE TENDENCY MEMORY
  ! => CALL mtend_register BEFORE
  ! => CALL mtend_request  AFTER
#ifdef MESSYTENDENCY
  ! - must be called after all submodels have placed their requests
  !   to alter prognostic variables (and tracers!) and after the
  !   number of tracers is fixated (i.e., after main_tracer_init_memory)
  ! - must be called after main_tracer_init_coupling to have correct
  !   meta-information of family member tracers
  ! - BUT before the coupling of the submodels, as channels are generated
  !   here which might be needed in other submodels (i.e. through mtend_request)
  CALL main_tendency_init_coupling
#endif
  IF (USE_ATTILA)      CALL      attila_init_coupling(2)
  IF (USE_CLAMS)       CALL       clams_init_coupling(2)
  IF (USE_H2OISO)      CALL      h2oiso_init_coupling(2)

  CALL main_rnd_init_coupling      ! create channel objects for state vectors
  CALL main_channel_init_coupling(2)

  CALL main_control_unset_context

END SUBROUTINE messy_init_coupling

!==============================================================================

SUBROUTINE messy_read_restart

  USE messy_main_control       ! ONLY: entrypoint, subentry
  USE messy_main_channel_bi,   ONLY: main_channel_read_restart
  USE messy_main_grid_bi,      ONLY: main_grid_read_restart

  CALL main_control_set_context(MEP_READ_RESTART)
  CALL main_channel_read_restart
  CALL main_grid_read_restart
  CALL main_control_unset_context

END SUBROUTINE messy_read_restart

!==============================================================================

SUBROUTINE messy_init_tracer

  ! tracer initialization (called by xtini in mo_tracer.f90)

  USE messy_main_switch   ! ONLY: USE_*
  USE messy_main_control  ! ONLY: entrypoint, subentry
  USE messy_main_tracer_bi, ONLY: main_tracer_init_tracer
  USE messy_main_rnd_bi,    ONLY: main_rnd_init_tracer
  USE messy_main_import_bi, ONLY: main_import_init_tracer

  USE messy_attila_e5,      ONLY:      attila_init_tracer
  USE messy_clams_si,       ONLY:       clams_init_tracer
!!$  USE messy_d14co_si,    ONLY:       d14co_init_tracer
!!$  USE messy_dradon_si,   ONLY:      dradon_init_tracer
  USE messy_e4chem_si,      ONLY:      e4chem_init_tracer
  USE messy_h2o_e5,         ONLY:         h2o_init_tracer
  USE messy_h2oiso_e5,      ONLY:      h2oiso_init_tracer
  USE messy_mecca_si,       ONLY:       mecca_init_tracer
  USE messy_tagging_si,     ONLY:     tagging_init_tracer
  USE messy_trsync_si,      ONLY:      trsync_init_tracer
!!$  USE messy_msbm_si,     ONLY:        msbm_init_tracer
!!$  USE messy_ptrac_si,    ONLY:       ptrac_init_tracer
!  USE messy_ptracini_si,   ONLY:    ptracini_init_tracer
  USE messy_clamscheme5_si, ONLY: clamscheme5_init_tracer
  USE messy_ubcnox_e5,      ONLY:      ubcnox_init_tracer
!!$  USE messy_clams_si,       ONLY:       clams_init_tracer

  IMPLICIT NONE

  CALL main_control_set_context(MEP_INIT_TRACER)

  CALL main_rnd_init_tracer       ! reset state vectors after restart
  CALL main_tracer_init_tracer(1) ! check tracer init from rerun
  CALL main_import_init_tracer

  ! special requirements:
  ! - attila_init_tracer must be called first (before
  !   main_tracer_init_tracer(4), to ensure that all fields for
  !   tracer initialisation (mass conserving transformation between gridpoint
  !   and Lagrangian space) are present
  IF (USE_ATTILA)  CALL  attila_init_tracer    ! lagrange !
  IF (USE_CLAMS)   CALL   clams_init_tracer    ! lagrange !
  CALL main_tracer_init_tracer(4) ! initialize via tracer_init (tracer.nml)

  ! CHECK / MODIFY INDIVIDUALLY
  ! special requirements:
  ! - (the rest is in alphabetical order)
!!$  IF (USE_D14CO)    CALL       d14co_init_tracer    ! lagrange !
!!$  IF (USE_DRADON)   CALL      dradon_init_tracer    ! lagrange !
  IF (USE_E4CHEM)      CALL      e4chem_init_tracer
  IF (USE_H2O)         CALL         h2o_init_tracer    ! lagrange !
  IF (USE_H2OISO)      CALL      h2oiso_init_tracer
  IF (USE_MECCA)       CALL       mecca_init_tracer    ! lagrange !
  IF (USE_TAGGING)     CALL     tagging_init_tracer
  IF (USE_TRSYNC)      CALL      trsync_init_tracer
!!$  IF (USE_MSBM)     CALL        msbm_init_tracer
!!$  IF (USE_PTRAC)    CALL       ptrac_init_tracer    ! lagrange !
! IF (USE_PTRACINI)    CALL    ptracini_init_tracer    ! lagrange !
  IF (USE_CLAMSCHEME5) CALL clamscheme5_init_tracer
  IF (USE_UBCNOX)      CALL      ubcnox_init_tracer
!!$  IF (USE_CLAMS)       CALL       clams_init_tracer

  CALL main_tracer_init_tracer(2) ! check tracer init by init_tracer
                                  ! set to constant if not yet initialized
  CALL main_tracer_init_tracer(3) ! diagnose tracer initialization

  CALL main_control_unset_context

END SUBROUTINE messy_init_tracer

!==============================================================================

SUBROUTINE messy_time(flag)

  USE messy_main_switch
  USE messy_main_control
  USE messy_oasis3mct_si,  ONLY: oasis3mct_time

  USE messy_main_timer_bi, ONLY: main_timer_time

  INTEGER, INTENT(IN) :: flag

  CALL main_control_set_context(MEP_TIME,flag)

  CALL main_timer_time(flag)      ! 1: set, 2: reset

  IF (USE_OASIS3MCT .AND. flag == 1) CALL oasis3mct_time(MEP_TIME,MEP_ONE)

  CALL main_control_unset_context

END SUBROUTINE messy_time

!==============================================================================

#ifdef MESSYTENDENCY
SUBROUTINE messy_tendency_reset

  !USE messy_main_switch    ! ONLY: USE_*
  USE messy_main_control  ! ONLY: entrypoint, subentry
   USE messy_main_tendency_bi,  ONLY: main_tendency_reset

   CALL main_control_set_context(MEP_TENDENCY_RESET)
   CALL main_tendency_reset
   CALL main_control_unset_context

END SUBROUTINE messy_tendency_reset
#endif

!==============================================================================

SUBROUTINE messy_global_start

  ! (called by stepon)

  USE messy_main_switch    ! ONLY: USE_*
  USE messy_main_control   ! ONLY: entrypoint, subentry
  USE messy_main_channel_bi, ONLY: main_channel_global_start
  USE messy_main_data_bi,    ONLY: main_data_global_start
  USE messy_main_grid_bi,    ONLY: main_grid_global_start
  USE messy_main_tracer_bi,  ONLY: main_tracer_global_start
  USE messy_main_import_bi,  ONLY: main_import_global_start

  USE messy_a2o_e5,        ONLY:       a2o_global_start
!!$  USE messy_attila_e5,  ONLY:    attila_global_start
  USE messy_bufly_si,      ONLY:     bufly_global_start
  USE messy_cat_e5,        ONLY:       cat_global_start
  USE messy_ch4_si,        ONLY:       ch4_global_start
  !USE messy_chemglue_si,  ONLY:  chemglue_global_start
  USE messy_clams_si,      ONLY:     clams_global_start
  USE messy_clamstraj_si,  ONLY: clamstraj_global_start
  USE messy_clamssedi_si,  ONLY: clamssedi_global_start
  USE messy_clamsmix_si,   ONLY:  clamsmix_global_start
  USE messy_cloud_si,      ONLY:     cloud_global_start
  USE messy_cloudj_si,     ONLY:    cloudj_global_start
  USE messy_clamstracer_si,ONLY: clamstracer_global_start
  USE messy_d14co_si,      ONLY:     d14co_global_start
  USE messy_ddep_si,       ONLY:      ddep_global_start
  USE messy_dissoc_si,     ONLY:    dissoc_global_start
  USE messy_domint_si,     ONLY:    domint_global_start
  USE messy_dradon_si,     ONLY:    dradon_global_start
  USE messy_e5vdiff_si,    ONLY:   e5vdiff_global_start
  USE messy_ec2cosmo_e5,   ONLY:  ec2cosmo_global_start
  USE messy_ever_si,       ONLY:      ever_global_start
  USE messy_gmxe_si,       ONLY:      gmxe_global_start
  USE messy_h2o_e5,        ONLY:       h2o_global_start
!!$  USE messy_h2oemis_si,    ONLY:   h2oemis_global_start
  USE messy_h2oiso_e5,     ONLY:    h2oiso_global_start
  USE messy_hamocc_e5,     ONLY:    hamocc_global_start
  USE messy_isopcor_si,    ONLY:   isopcor_global_start
  USE messy_jval_si,       ONLY:      jval_global_start
  USE messy_jvst_si,       ONLY:      jvst_global_start
  USE messy_lgvflux_e5,    ONLY:   lgvflux_global_start
  USE messy_mlocean_e5,    ONLY:   mlocean_global_start
  USE messy_mmforce_e5,    ONLY:   mmforce_global_start
  USE messy_mpiom_e5,      ONLY:     mpiom_global_start
  USE messy_mtskip_si,     ONLY:    mtskip_global_start
  USE messy_oasis3mct_si,  ONLY: oasis3mct_global_start
  USE messy_onemis_si,     ONLY:    onemis_global_start
  USE messy_orbit_si,      ONLY:     orbit_global_start
  USE messy_otphysc_e5,    ONLY:   otphysc_global_start
  USE messy_ptracini_si,   ONLY:  ptracini_global_start
  USE messy_rad_si,        ONLY:       rad_global_start
  USE messy_s4d_si,        ONLY:       s4d_global_start
  USE messy_scalc_si,      ONLY:     scalc_global_start
  USE messy_sorbit_si,     ONLY:    sorbit_global_start
  USE messy_spe_e5,        ONLY:       spe_global_start
  USE messy_tagging_si,    ONLY:   tagging_global_start
  USE messy_timepos_e5,    ONLY:   timepos_global_start
  USE messy_tpulse_si,     ONLY:    tpulse_global_start
  USE messy_tnudge_si,     ONLY:    tnudge_global_start
  USE messy_trexp_si,      ONLY:     trexp_global_start
  USE messy_ubcnox_e5,     ONLY:    ubcnox_global_start
  USE messy_vertex_si,     ONLY:    vertex_global_start
  USE messy_mmd2way_si,    ONLY:   mmd2way_global_start

  IMPLICIT NONE

  CALL main_control_set_context(MEP_GLOBAL_START)
  CALL main_channel_global_start
  CALL main_data_global_start
  CALL main_grid_global_start  ! GRID requires vars from DATA
  CALL main_tracer_global_start
  CALL main_import_global_start

  ! special requirements:
  ! - bufly_global_start must be called very first, because it modifies the
  !   initial condition
!!$  ! - attila_global_start must be called first to ensure that all fields for
!!$  !   mass conserving transformation between gridpoint and
!!$  !   Lagrangian space are present
  !   (cdisse must be up to date)
  ! - RAD must be called before and after ORBIT (with flag = 1,2 respectively)
  !       1. RAD calculates trigger and time offset
  !       2. ORBIT calculates orbit parameters with time offset
  !       3. RAD calculates radiation
  ! - MTSKIP must be called before orbit to update the time offset correctly
  ! - PTRACINI should be called last before the coupling routines to ensure that
  !   the necessary fields are calculated
  ! - scalc must be called first to have updated objects from last time step
  ! - (the rest is in alphabetical order)
  IF (USE_SCALC)    CALL     scalc_global_start
  IF (USE_OASIS3MCT) CALL oasis3mct_global_start(MEP_GLOBAL_START,MEPS_BEGIN)
  IF (USE_BUFLY)    CALL     bufly_global_start
!!$  IF (USE_ATTILA) CALL   attila_global_start    ! lagrange !
  IF (USE_CAT)      CALL       cat_global_start
  IF (USE_CH4)      CALL       ch4_global_start
  !IF (USE_CHEMGLUE) CALL chemglue_global_start
  IF (USE_CLAMS)    CALL     clams_global_start
  IF (USE_CLAMSSEDI)CALL clamssedi_global_start
  IF (USE_CLAMSTRAJ)CALL clamstraj_global_start
  IF (USE_CLAMSMIX) CALL  clamsmix_global_start
  IF (USE_CLAMSTRACER) CALL clamstracer_global_start
  IF (USE_CLOUD)    CALL     cloud_global_start
  IF (USE_D14CO)    CALL     d14co_global_start    ! lagrange !
  IF (USE_DDEP)     CALL      ddep_global_start    ! lagrange !
  IF (USE_DISSOC)   CALL    dissoc_global_start
  IF (USE_DOMINT)   CALL    domint_global_start
  IF (USE_DRADON)   CALL    dradon_global_start    ! lagrange !
  IF (USE_E5VDIFF)  CALL   e5vdiff_global_start
  IF (USE_EVER)     CALL      ever_global_start
  IF (USE_GMXE)     CALL      gmxe_global_start
  IF (USE_H2O)      CALL       h2o_global_start    ! lagrange !
!!$  IF (USE_H2OEMIS)  CALL   h2oemis_global_start
  IF (USE_H2OISO)   CALL    h2oiso_global_start
  IF (USE_ISOPCOR)  CALL   isopcor_global_start
  IF (USE_LGVFLUX)  CALL   lgvflux_global_start    ! lagrange !
  IF (USE_MLOCEAN)  CALL   mlocean_global_start
  IF (USE_MMFORCE)  CALL   mmforce_global_start
  IF (USE_MPIOM)    CALL     mpiom_global_start
  IF (USE_HAMOCC)   CALL    hamocc_global_start
  IF (USE_OTPHYSC)  CALL   otphysc_global_start
  IF (USE_A2O)      CALL       a2o_global_start
  IF (USE_ONEMIS)   CALL    onemis_global_start    ! lagrange !
  IF (USE_RAD)      CALL       rad_global_start(1)
  IF (USE_MTSKIP)   CALL    mtskip_global_start
  IF (USE_ORBIT)    CALL     orbit_global_start
  IF (USE_RAD)      CALL       rad_global_start(2)
  IF (USE_CLOUDJ)   CALL    cloudj_global_start    ! lagrange !
  IF (USE_JVAL)     CALL      jval_global_start    ! lagrange !
  IF (USE_JVST)     CALL      jvst_global_start    ! lagrange !
  IF (USE_S4D)      CALL       s4d_global_start
  IF (USE_SORBIT)   CALL    sorbit_global_start
  IF (USE_SPE)      CALL       spe_global_start
  IF (USE_TAGGING)  CALL   tagging_global_start
  IF (USE_TIMEPOS)  CALL   timepos_global_start
  IF (USE_TNUDGE)   CALL    tnudge_global_start
  IF (USE_TREXP)    CALL     trexp_global_start
  IF (USE_TPULSE)   CALL    tpulse_global_start
  IF (USE_UBCNOX)   CALL    ubcnox_global_start
  IF (USE_VERTEX)   CALL    vertex_global_start
  IF (USE_PTRACINI) CALL  ptracini_global_start
  IF (USE_EC2COSMO) CALL  ec2cosmo_global_start
  IF (USE_MMD2WAY)  CALL   mmd2way_global_start
  IF (USE_OASIS3MCT) CALL oasis3mct_global_start(MEP_GLOBAL_START,MEPS_END)

  CALL main_control_unset_context

END SUBROUTINE messy_global_start

!==============================================================================

SUBROUTINE messy_beforeadv

  !USE messy_main_switch   ! ONLY: USE_*
  USE messy_main_control  ! ONLY: entrypoint, subentry
  USE messy_main_tracer_bi, ONLY: main_tracer_beforeadv

  CALL main_control_set_context(MEP_BEFOREADV)
  CALL main_tracer_beforeadv
  CALL main_control_unset_context

END SUBROUTINE messy_beforeadv

!==============================================================================

SUBROUTINE messy_afteradv

  !USE messy_main_switch   ! ONLY: USE_*
  USE messy_main_control  ! ONLY: entrypoint, subentry
  USE messy_main_tracer_bi, ONLY: main_tracer_afteradv

  CALL main_control_set_context(MEP_AFTERADV)
  CALL main_tracer_afteradv
  CALL main_control_unset_context

END SUBROUTINE messy_afteradv

!==============================================================================

SUBROUTINE messy_local_start

  ! (will be called by scan1)

  USE messy_main_switch   ! ONLY: USE_*
  USE messy_main_control  ! ONLY: entrypoint, subentry
  USE messy_main_data_bi,   ONLY: main_data_local_start
  USE messy_main_grid_bi,   ONLY: main_grid_local_start
  USE messy_main_tracer_bi, ONLY: main_tracer_local_start

  USE messy_attila_e5,      ONLY: attila_local_start
  USE messy_clams_si,       ONLY: clams_local_start
  USE messy_diumod_si,      ONLY: diumod_local_start
  USE messy_e4chem_si,      ONLY: e4chem_local_start
  USE messy_h2oiso_e5,      ONLY: h2oiso_local_start
  USE messy_msbm_si,        ONLY: msbm_local_start
  USE messy_mtskip_si,      ONLY: mtskip_local_start
  USE messy_tagging_si,     ONLY: tagging_local_start

  IMPLICIT NONE

  CALL main_control_set_context(MEP_LOCAL_START)
  CALL main_grid_local_start ! GRID needs to be called first !
  CALL main_data_local_start
  CALL main_tracer_local_start

  ! special requirements:
  ! - call msbm_local_start before modules which need flt_pscreg
  IF (USE_ATTILA)  CALL  attila_local_start
  IF (USE_CLAMS)   CALL   clams_local_start
  IF (USE_MSBM)    CALL    msbm_local_start
  IF (USE_DIUMOD)  CALL  diumod_local_start
  IF (USE_E4CHEM)  CALL  e4chem_local_start
  IF (USE_H2OISO)  CALL  h2oiso_local_start
  IF (USE_TAGGING) CALL tagging_local_start
  IF (USE_MTSKIP)  CALL  mtskip_local_start

  CALL main_control_unset_context

END SUBROUTINE messy_local_start

!==============================================================================

SUBROUTINE messy_radiation

  USE messy_main_switch   ! ONLY: USE_*
  USE messy_main_control  ! ONLY: entrypoint, subentry

  USE messy_main_grid_bi,  ONLY: main_grid_radiation

  USE messy_aeropt_si,     ONLY:   aeropt_radiation
  USE messy_cloud_si,      ONLY:    cloud_radiation
  USE messy_cloudopt_si,   ONLY: cloudopt_radiation
  USE messy_crm_e5,        ONLY:      crm_radiation
  USE messy_gmxe_si,       ONLY:     gmxe_radiation
  USE messy_nan_si,        ONLY:      nan_radiation
  USE messy_rad_si,        ONLY:      rad_radiation
  USE messy_surface_si,    ONLY:  surface_radiation

  IMPLICIT NONE

  CALL main_control_set_context(MEP_RADIATION)
  CALL main_grid_radiation

  ! special requirements:
  ! - surface needs to be called first
  !   (to keep original ECHAM5 sequence of operations)
  ! - aeropt should be called directly before radiation
  !   such that no other process can influence the aerosol
  !   distribution after the optical properties have been determined
  ! - nan_radiation(1) and (2) must be called before and after gmxe,
  !   respectively
  IF (USE_SURFACE)  CALL  surface_radiation
  IF (USE_CLOUD)    CALL    cloud_radiation
  IF (USE_NAN)      CALL      nan_radiation(1)
  IF (USE_GMXE)     CALL     gmxe_radiation
  IF (USE_NAN)      CALL      nan_radiation(2)
  IF (USE_CRM)      CALL      crm_radiation
  IF (USE_AEROPT)   CALL   aeropt_radiation
  IF (USE_CLOUDOPT) CALL cloudopt_radiation
  IF (USE_RAD)      CALL      rad_radiation

  CALL main_control_unset_context

END SUBROUTINE messy_radiation

!==============================================================================

SUBROUTINE messy_vdiff

  ! (called by physc)   ! mz_ho_20160412

  USE messy_main_switch   ! ONLY: USE_*
  USE messy_main_control  ! ONLY: entrypoint, subentry

  USE messy_airsea_si,         ONLY:     airsea_vdiff
  USE messy_bioburn_si,        ONLY:    bioburn_vdiff
  USE messy_ddep_si,           ONLY:       ddep_vdiff
  USE messy_dradon_si,         ONLY:     dradon_vdiff
  USE messy_e5vdiff_si,        ONLY:    e5vdiff_vdiff
  USE messy_gmxe_si,           ONLY:       gmxe_vdiff
  USE messy_ch4_si,            ONLY:        ch4_vdiff
  USE messy_m7_si,             ONLY:         m7_vdiff
  USE messy_made_si,           ONLY:       made_vdiff
  USE messy_made3_si,          ONLY:      made3_vdiff
  USE messy_mecca_si,          ONLY:      mecca_vdiff
  USE messy_megan_si,          ONLY:      megan_vdiff
  USE messy_offemis_si,        ONLY:    offemis_vdiff
  USE messy_onemis_si,         ONLY:     onemis_vdiff
  USE messy_oracle_si,         ONLY:     oracle_vdiff
  USE messy_ptracini_si,       ONLY:   ptracini_vdiff
  USE messy_svoc_si,           ONLY:       svoc_vdiff
  USE messy_tagging_si,        ONLY:    tagging_vdiff
  USE messy_tropop_si,         ONLY:     tropop_vdiff
  USE messy_vertex_si,         ONLY:     vertex_vdiff
  USE messy_viso_si,           ONLY:       viso_vdiff

  IMPLICIT NONE

  CALL main_control_set_context(MEP_VDIFF)
  ! special requirements:

  ! - call mecca_vdiff(1) after emissions (ONLEM, OFFLEM, OFFEMIS, AIRSEA,
  !   MEGAN)
  !   and before deposition (DRYDEP) to get the "real" new emitted aerosol
  ! - call mecca_vdiff(2) after DRYDEP
  ! - call m7_vdiff and gmxe_vdiff after emissions
  !   (ONLEM, OFFLEM, OFFEMIS, AIRSEA, MEGAN) to get the
  !   current emissions if m7/gmxe calculates its own emissions
  ! - tropop_vdiff must be called before viso_vdiff
  ! - viso_vdiff and tropop_vdiff
  !   must be called first, since they calulate diagnostic quantities
  !   needed in other submodels
  ! - e5vdiff_vdiff(1) or vertex_vdiff(1) must be called first
  ! - e5vdiff_vdiff(2) or vertex_vdiff(2) must be called last
  ! - call svoc_vdiff(1) before AIRSEA and after ONEMIS and OFFEMIS
  !   and call svoc_vdiff(2) called after DDEP and M7/GMXE
  ! - bioburn must be called before offemis for up-to-date fluxes

  IF (USE_E5VDIFF) CALL e5vdiff_vdiff(1)
  IF (USE_VERTEX)  CALL  vertex_vdiff(1)
  IF (USE_CH4)     CALL     ch4_vdiff
  IF (USE_TROPOP)  CALL  tropop_vdiff
  IF (USE_VISO)    CALL    viso_vdiff
  IF (USE_BIOBURN) CALL bioburn_vdiff
  IF (USE_DRADON)  CALL  dradon_vdiff    ! lagrange !
  IF (USE_OFFEMIS) CALL offemis_vdiff
  IF (USE_TAGGING) CALL tagging_vdiff(1)
  IF (USE_ONEMIS)  CALL  onemis_vdiff    ! lagrange !
  IF (USE_SVOC)    CALL    svoc_vdiff(1)
  IF (USE_AIRSEA)  CALL  airsea_vdiff    ! lagrange !
  IF (USE_TAGGING) CALL tagging_vdiff(2)
  IF (USE_MEGAN)   CALL   megan_vdiff
  IF (USE_M7)      CALL      m7_vdiff
  IF (USE_MADE)    CALL    made_vdiff
  IF (USE_MADE3)   CALL   made3_vdiff
  IF (USE_GMXE)    CALL    gmxe_vdiff
  IF (USE_MECCA)   CALL   mecca_vdiff(1)
  IF (USE_DDEP)    CALL    ddep_vdiff    ! lagrange !
  IF (USE_SVOC)    CALL    svoc_vdiff(2)
  IF (USE_MECCA)   CALL   mecca_vdiff(2)
  IF (USE_ORACLE)  CALL  oracle_vdiff
  IF (USE_PTRACINI) CALL ptracini_vdiff
  IF (USE_E5VDIFF) CALL e5vdiff_vdiff(2)
  IF (USE_VERTEX)  CALL  vertex_vdiff(2)

  CALL main_control_unset_context

END SUBROUTINE messy_vdiff

!==============================================================================

SUBROUTINE messy_radheat

  USE messy_main_switch   ! ONLY: USE_*
  USE messy_main_control  ! ONLY: entrypoint, subentry

  USE messy_h2oiso_e5,     ONLY:  h2oiso_radheat
  USE messy_rad_si,        ONLY:     rad_radheat

  IMPLICIT NONE

  CALL main_control_set_context(MEP_RADHEAT)
  !op_re_20140919+
  ! h2oiso_radheat must be called here (right after vdiff),
  ! since a call within messy_vdiff would come too early
  IF (USE_H2OISO)  CALL  h2oiso_radheat
  !op_re_20140919-

  ! special requirements:
  IF (USE_RAD)     CALL     rad_radheat

  CALL main_control_unset_context

END SUBROUTINE messy_radheat

!==============================================================================

SUBROUTINE messy_gwdrag

  USE messy_main_switch   ! ONLY: USE_*
  USE messy_main_control  ! ONLY: entrypoint, subentry

  USE messy_gwave_si,       ONLY: gwave_physc
  USE messy_orogw_si,       ONLY: orogw_physc

  IMPLICIT NONE

  CALL main_control_set_context(MEP_GWDRAG)

  ! special requirements:
  IF (USE_GWAVE) CALL gwave_physc
  IF (USE_OROGW) CALL orogw_physc

  CALL main_control_unset_context

END SUBROUTINE messy_gwdrag

!==============================================================================

SUBROUTINE messy_convec

  ! (called by physc before convection subroutine cucall)

  USE messy_main_switch   ! ONLY: USE_*
  USE messy_main_control  ! ONLY: entrypoint, subentry

  USE messy_contrail_si,    ONLY: contrail_convec
  USE messy_convect_si,     ONLY:  convect_convec
  USE messy_cvtrans_si,     ONLY:  cvtrans_convec
  USE messy_h2oiso_e5,      ONLY:   h2oiso_convec
  USE messy_cloud_si,       ONLY:    cloud_convec
  USE messy_scalc_si,       ONLY:    scalc_convec
  USE messy_scav_si,        ONLY:     scav_convec
  USE messy_svoc_si,        ONLY:     svoc_convec

  IMPLICIT NONE

  CALL main_control_set_context(MEP_CONVEC)

  ! special requirements:
  ! - convect_convec must be called first, because it calculates
  !   the convection and all other submodels use convection data
  ! - cvtrans_convec(use_scav,1), scav_convec, and cvtrans_convec(use_scav,2)
  !   must be in this order and together
  ! - h2oiso_convec must be called before convect_convec, because
  !   it pre-calculates various intermediate steps of convect_convec for
  !   isotopologues
  IF (USE_H2OISO)   CALL   h2oiso_convec
  IF (USE_CONVECT)  CALL  convect_convec
  IF (USE_CVTRANS)  CALL  cvtrans_convec(1)
  IF (USE_SCAV)     CALL     scav_convec
  IF (USE_SVOC)     CALL     svoc_convec
  IF (USE_CVTRANS)  CALL  cvtrans_convec(2)
  IF (USE_CONTRAIL) CALL contrail_convec
  IF (USE_CLOUD)    CALL    cloud_convec

  IF (USE_SCALC)    CALL    scalc_convec

  CALL main_control_unset_context

END SUBROUTINE messy_convec

!==============================================================================

SUBROUTINE messy_mixlo

  USE messy_main_switch   ! ONLY: USE_*
  USE messy_main_control  ! ONLY: entrypoint, subentry
  USE messy_h2oiso_e5,     ONLY: h2oiso_mixlo
  USE messy_mlocean_e5,    ONLY: mlocean_mixlo
  USE messy_surface_si,    ONLY: surface_mixlo

  CALL main_control_set_context(MEP_MIXLO)
  ! special requirements:
  !  h2oiso_mixlo must be called right after surface_mixlo(1)
  IF (USE_SURFACE) CALL  surface_mixlo(1) ! surf, lake, licetemp
  IF (USE_H2OISO)  CALL  h2oiso_mixlo
  IF (USE_MLOCEAN) CALL  mlocean_mixlo
  IF (USE_SURFACE) CALL  surface_mixlo(2) ! sicetemp, tsurf

  CALL main_control_unset_context

END SUBROUTINE messy_mixlo

!==============================================================================

SUBROUTINE messy_physc

  ! main submodel call (called by xtdriver in mo_tracer)

  ! ONLY FOR DEBUGGING +
  USE messy_main_grid_def_mem_bi, ONLY: jrow, ngpblks, kproma
  USE messy_main_timer,     ONLY: YEAR, MONTH, DAY, HOUR, MINUTE, SECOND
  USE messy_main_mpi_bi,    ONLY: p_pe
  ! ONLY FOR DEBUGGING -

  USE messy_main_switch   ! ONLY: USE_*
  USE messy_main_control  ! ONLY: entrypoint, subentry

  USE messy_accf_e5,        ONLY:       accf_physc
  USE messy_airtraf_e5,     ONLY:    airtraf_physc
  USE messy_bufly_si,       ONLY:      bufly_physc
  USE messy_cat_e5,         ONLY:        cat_physc
  USE messy_ch4_si,         ONLY:        ch4_physc
  USE messy_chemglue_si,    ONLY:   chemglue_physc
  USE messy_cloudj_si,      ONLY:     cloudj_physc
  USE messy_contrail_si,    ONLY:   contrail_physc
  USE messy_d14co_si,       ONLY:      d14co_physc
  USE messy_dissoc_si,      ONLY:     dissoc_physc
  USE messy_dradon_si,      ONLY:     dradon_physc
  USE messy_e4chem_si,      ONLY:     e4chem_physc
  USE messy_edith_e5,       ONLY:      edith_physc
  USE messy_ever_si,        ONLY:       ever_physc
  USE messy_gec_si,         ONLY:        gec_physc
  USE messy_gmxe_si,        ONLY:       gmxe_physc
  USE messy_h2o_e5,         ONLY:        h2o_physc
  USE messy_h2oemis_si,     ONLY:    h2oemis_physc
  USE messy_h2oiso_e5,      ONLY:     h2oiso_physc
  USE messy_isopcor_si,     ONLY:    isopcor_physc
  USE messy_ions_si,        ONLY:       ions_physc
  USE messy_jval_si,        ONLY:       jval_physc
  USE messy_jvst_si,        ONLY:       jvst_physc
  USE messy_lnox_si,        ONLY:       lnox_physc
  USE messy_m7_si,          ONLY:         m7_physc
  USE messy_made_si,        ONLY:       made_physc
  USE messy_made3_si,       ONLY:      made3_physc
  USE messy_mecca_si,       ONLY:      mecca_physc
  USE messy_mesoenergy_e5,  ONLY: mesoenergy_physc
  USE messy_mmforce_e5,     ONLY:    mmforce_physc
  USE messy_msbm_si,        ONLY:       msbm_physc
  USE messy_nan_si,         ONLY:        nan_physc
  USE messy_o3orig_si,      ONLY:     o3orig_physc
  USE messy_oracle_si,      ONLY:     oracle_physc
  USE messy_plumegas_si,    ONLY:   plumegas_physc
  USE messy_ptracini_si,    ONLY:   ptracini_physc
  USE messy_qbo_si,         ONLY:        qbo_physc
  USE messy_relax_e5,       ONLY:      relax_physc
  USE messy_satsims_e5,     ONLY:    satsims_physc
  USE messy_scav_si,        ONLY:       scav_physc
  USE messy_scalc_si,       ONLY:      scalc_physc
  USE messy_sedi_si,        ONLY:       sedi_physc
  USE messy_sf6_e5,         ONLY:        sf6_physc
  USE messy_spe_e5,         ONLY:        spe_physc
  USE messy_spacenox_e5,    ONLY:   spacenox_physc
  USE messy_svoc_si,        ONLY:       svoc_physc
  USE messy_tagging_si,     ONLY:    tagging_physc
  USE messy_trexp_si,       ONLY:      trexp_physc
  USE messy_trsync_si,      ONLY:     trsync_physc
  USE messy_ubcnox_e5,      ONLY:     ubcnox_physc
  USE messy_vahr_e5,        ONLY:       vahr_physc

  IMPLICIT NONE

  ! ONLY FOR DEBUGGING +
  IF (L_TIME_INFO) &
  WRITE(*, '(1x,i4,a1,i2.2,a1,i2.2,a1,i2.2,a1,i2.2,a1,i2.2,a3'//&
       &',a4,i4,a8,i5,a3,i5,a10,i5)') &
       YEAR,'-',MONTH,'-',DAY,' ',HOUR,':', MINUTE,':',SECOND,' # ' &
       ,' PE=',p_pe,' # JROW=',jrow,' / ',ngpblks,' # KPROMA=',kproma
  ! ONLY FOR DEBUGGING -

  ! special requirements:
  ! - bufly_physc must be called first, because it modifies the initial cond.
  ! - cloudj_physc must be called before mecca_physc
  ! - jval_physc must be called before mecca_physc
  ! - msbm_physc must be called before mecca_physc
  ! - chemglue_physc must be called before mecca_physc
  ! - scav_physc(1) must be called before mecca_physc, m7_physc
  ! - scav_physc(2) must be called last, since species in the cloud phases are
  !   released back into the gas phase
  ! - sedi_physc must be called before mecca_physc
  ! - m7_physc and gmxe_physc should be called before mecca_physc
  !   - after all aerosol modules (PSC, MSBM, M7, PTRAC)
  !   - directly before mecca_physc, to avoid problems due to the conversion
  !     of first order reaction rates in second order.
  ! - d14co_physc and trexp_physc must be called at the end, i.e., after all
  !   other chemistry
  ! - made_physc and made3_physc must be called after mecca_physc
  !   because of H2SO4 mixing ratio
  ! - trsync_physc(1) must be called before ch4 and mecca
  ! - trsync_physc(2) must be called after ch4 and mecca
  ! - h2oiso_physc must be called after trsync_physc(2)
  ! - trsync must be called before (1) and after (2) chemistry
  !   (mecca, ch4, e4chem)
  ! - contrail_physc called at the end to minimize difference between potcov(gp)
  !   and potcov(lg)
  ! - svoc_physc(1) must be called after scav_physc(1), svoc_physc(2) must be
  !   called before mecca_physc and svoc_physc(3) after scav_physc(2)

  CALL main_control_set_context(MEP_PHYSC)

  IF (USE_SCALC)      CALL    scalc_physc ! requrired for ECHAM-FESOM coupling
  IF (USE_ACCF)       CALL     accf_physc
  IF (USE_AIRTRAF)    CALL  airtraf_physc
  IF (USE_BUFLY)      CALL    bufly_physc
  IF (USE_CAT)        CALL      cat_physc
!  IF (USE_CONTRAIL)  CALL contrail_physc
  IF (USE_EDITH)      CALL    edith_physc
  IF (USE_QBO)        CALL      qbo_physc
  IF (USE_SPE)        CALL      spe_physc
  IF (USE_SPACENOX)   CALL spacenox_physc
  IF (USE_DRADON)     CALL   dradon_physc    ! lagrange !
  IF (USE_IONS)       CALL     ions_physc
  IF (USE_GEC)        CALL      gec_physc
  IF (USE_H2O)        CALL      h2o_physc    ! lagrange !
  IF (USE_H2OEMIS)    CALL  h2oemis_physc
  IF (USE_LNOX)       CALL     lnox_physc    ! lagrange !
  IF (USE_CLOUDJ)     CALL   cloudj_physc    ! lagrange !
  IF (USE_ISOPCOR)    CALL  isopcor_physc
  IF (USE_JVAL)       CALL     jval_physc    ! lagrange !
  IF (USE_JVST)       CALL     jvst_physc    ! lagrange !
  IF (USE_DISSOC)     CALL   dissoc_physc
  IF (USE_MSBM)       CALL     msbm_physc
  IF (USE_H2OISO)     CALL   h2oiso_physc(1)
  IF (USE_SCAV)       CALL     scav_physc(1)
  IF (USE_SVOC)       CALL     svoc_physc(1)
  IF (USE_SEDI)       CALL     sedi_physc
!!$  IF (USE_CH4)        CALL      ch4_physc
!!$  IF (USE_H2OISO)     CALL   h2oiso_physc
  IF (USE_M7)         CALL       m7_physc
  IF (USE_GMXE)       CALL     gmxe_physc
  IF (USE_NAN)        CALL      nan_physc
  IF (USE_O3ORIG)     CALL   o3orig_physc(1)
  IF (USE_E4CHEM)     CALL   e4chem_physc
  IF (USE_TAGGING)    CALL  tagging_physc(1)
  IF (USE_TRSYNC)     CALL   trsync_physc(1)
  IF (USE_CH4)        CALL      ch4_physc
  IF (USE_CHEMGLUE)   CALL chemglue_physc
  IF (USE_SVOC)       CALL     svoc_physc(2)
  IF (USE_MECCA)      CALL    mecca_physc    ! lagrange !
  IF (USE_TRSYNC)     CALL   trsync_physc(2)
  IF (USE_TAGGING)    CALL  tagging_physc(2)
  IF (USE_H2OISO)     CALL   h2oiso_physc(2)
  IF (USE_MADE)       CALL     made_physc
  IF (USE_MADE3)      CALL    made3_physc
  IF (USE_MESOENERGY) CALL mesoenergy_physc
  IF (USE_O3ORIG)     CALL   o3orig_physc(2)
  IF (USE_ORACLE)     CALL   oracle_physc
  IF (USE_MMFORCE)    CALL  mmforce_physc
  IF (USE_TREXP)      CALL    trexp_physc
  IF (USE_EVER)       CALL     ever_physc
  IF (USE_PTRACINI)   CALL ptracini_physc
  IF (USE_D14CO)      CALL    d14co_physc    ! lagrange !
  IF (USE_RELAX)      CALL    relax_physc
  IF (USE_SATSIMS)    CALL  satsims_physc
  IF (USE_SCAV)       CALL     scav_physc(2)
  IF (USE_SVOC)       CALL     svoc_physc(3)
  IF (USE_SF6)        CALL      sf6_physc
  IF (USE_UBCNOX)     CALL   ubcnox_physc
  IF (USE_VAHR)       CALL     vahr_physc
  IF (USE_PLUMEGAS)   CALL plumegas_physc
  IF (USE_CONTRAIL)   CALL contrail_physc

  CALL main_control_unset_context

END SUBROUTINE messy_physc

!==============================================================================

SUBROUTINE messy_local_end

  ! diagnostics at the end of a time step (called by xtdiagn in mo_tracer)

  USE messy_main_switch   ! ONLY: USE_*
  USE messy_main_control  ! ONLY: entrypoint, subentry

  USE messy_attila_e5,       ONLY:  attila_local_end
  USE messy_clams_si,        ONLY:   clams_local_end
  USE messy_h2o_e5,          ONLY:     h2o_local_end
  USE messy_m7_si,           ONLY:      m7_local_end
  USE messy_msbm_si,         ONLY:    msbm_local_end
  USE messy_tbudget_si,      ONLY: tbudget_local_end
  USE messy_tnudge_si,       ONLY:  tnudge_local_end
  USE messy_viso_si,         ONLY:    viso_local_end

  IMPLICIT NONE

  CALL main_control_set_context(MEP_LOCAL_END)
  ! special requirements:
  ! - tnudge_local_end and viso_local_end must be called last,
  !   since all channel objects must be up to date
  ! - (the rest is in alphabetical order)
  IF (USE_ATTILA)  CALL  attila_local_end
  IF (USE_CLAMS)   CALL   clams_local_end
  IF (USE_H2O)     CALL     h2o_local_end    ! lagrange !
  IF (USE_MSBM)    CALL    msbm_local_end
  IF (USE_M7)      CALL      m7_local_end
  IF (USE_TBUDGET) CALL tbudget_local_end
  IF (USE_TNUDGE)  CALL  tnudge_local_end
  IF (USE_VISO)    CALL    viso_local_end

  CALL main_control_unset_context

END SUBROUTINE messy_local_end

!==============================================================================

SUBROUTINE messy_global_end

  ! diagnostics at the end of a time step; outside the region loop
  ! (called by stepon)

  USE messy_main_switch    ! ONLY: USE_*
  USE messy_main_control   ! ONLY: entrypoint, subentry
  USE messy_main_qtimer_bi,  ONLY: main_qtimer_global_end
  USE messy_main_tracer_bi,  ONLY: main_tracer_global_end
!!$  USE messy_main_channel_bi, ONLY: main_channel_global_end
#ifdef MESSYTENDENCY
  USE messy_main_tendency_bi,ONLY: main_tendency_global_end
#endif

  USE messy_a2o_e5,         ONLY:         a2o_global_end
  USE messy_attila_e5,      ONLY:      attila_global_end
  USE messy_airsea_si,      ONLY:      airsea_global_end
  USE messy_airtrac_e5,     ONLY:     airtrac_global_end
  USE messy_airtraf_e5,     ONLY:     airtraf_global_end
  USE messy_aveout_si,      ONLY:      aveout_global_end
  USE messy_bioburn_si,     ONLY:     bioburn_global_end
  USE messy_ch4_si,         ONLY:         ch4_global_end
  USE messy_clams_si,       ONLY:       clams_global_end
  USE messy_clamstraj_si,   ONLY:   clamstraj_global_end
  USE messy_clamschem_si,   ONLY:   clamschem_global_end
  USE messy_clamscheme5_si, ONLY: clamscheme5_global_end
  USE messy_clamssedi_si,   ONLY:   clamssedi_global_end
  USE messy_clamsmix_si,    ONLY:    clamsmix_global_end
  USE messy_clamsbmix_si,   ONLY:   clamsbmix_global_end
  USE messy_clamscirrus_si, ONLY: clamscirrus_global_end
  USE messy_clamsrdfrc_si,  ONLY:  clamsrdfrc_global_end
  USE messy_clamstracer_si, ONLY: clamstracer_global_end
  USE messy_cloudj_si,      ONLY:      cloudj_global_end
  USE messy_contrail_si,    ONLY:    contrail_global_end
  USE messy_cvtrans_si,     ONLY:     cvtrans_global_end
  USE messy_dissoc_si,      ONLY:      dissoc_global_end
  USE messy_ddep_si,        ONLY:        ddep_global_end
  USE messy_domint_si,      ONLY:      domint_global_end
  USE messy_dradon_si,      ONLY:      dradon_global_end
  USE messy_d14co_si,       ONLY:       d14co_global_end
  USE messy_e5vdiff_si,     ONLY:     e5vdiff_global_end
  USE messy_gwave_si,       ONLY:       gwave_global_end
  USE messy_h2o_e5,         ONLY:         h2o_global_end
  USE messy_hd_e5,          ONLY:          hd_global_end
  USE messy_hamocc_e5,      ONLY:      hamocc_global_end
  USE messy_jval_si,        ONLY:        jval_global_end
  USE messy_jvst_si,        ONLY:        jvst_global_end
  USE messy_lggp_e5,        ONLY:        lggp_global_end
  USE messy_lgtmix_e5,      ONLY:      lgtmix_global_end
  USE messy_lgvflux_e5,     ONLY:     lgvflux_global_end
  USE messy_lnox_si,        ONLY:        lnox_global_end
  USE messy_mlocean_e5,     ONLY:     mlocean_global_end
  USE messy_mecca_si,       ONLY:       mecca_global_end
  USE messy_megan_si,       ONLY:       megan_global_end
  USE messy_mpiom_e5,       ONLY:       mpiom_global_end
  USE messy_mmd2way_si,     ONLY:     mmd2way_global_end
  USE messy_offemis_si,     ONLY:     offemis_global_end
  USE messy_onemis_si,      ONLY:      onemis_global_end
  USE messy_otphysc_e5,     ONLY:     otphysc_global_end
  USE messy_rndtest_si,     ONLY:     rndtest_global_end
  USE messy_scalc_si,       ONLY:       scalc_global_end
!  USE messy_scav_si,        ONLY:        scav_global_end
  USE messy_sedi_si,        ONLY:        sedi_global_end
  USE messy_surface_si,     ONLY:     surface_global_end
  USE messy_tnudge_si,      ONLY:      tnudge_global_end
  USE messy_trexp_si,       ONLY:       trexp_global_end
  USE messy_vaxtra_si,      ONLY:      vaxtra_global_end
  USE messy_vertex_si,      ONLY:      vertex_global_end

  IMPLICIT NONE

  CALL main_control_set_context(MEP_GLOBAL_END)

  ! special requirements:
  ! - attila_global_end must be called before all other lg-submodels
  !   since positions need to be updated
  ! - d14co_global_end and trexp_global_end must be called after
  !   all other lg-chemistry
  ! - cloudj_global_end must be called before mecca_global_end
  ! - jval_global_end must be called before mecca_global_end
  ! - lnox_global_end must be called before mecca_global_end
  ! - h2o_global_end must be called after mecca_global_end
  ! - lggp_global_end and lgvflux_global end must be the last (in this
  !   order) of all lg-submodels, since fields need to be up to date
  ! - ch4_global_end must be called after lggp_global_end to enable up-to-date
  !   fields from gp-to-lg transformation
  ! - mecca_global_end must be called after emissions (OFFLEM, OFFEMIS, ONLEM,
  !   AIRSEA, MEGAN) and before depositions (DRYDEP, DDEP) ; see vdiff for GP
  ! - scav_global_end(1) should be called before mecca_global_end to be
  !   consistent with grid point order
  ! - scav_global_end(2) should be called after all tracer modifying submodels
  ! - hd_global_end and a2o_global_end must be called before mpiom_global_end
  ! - scalc_global_end must be called after depositions (DRYDEP, DDEP) and
  !   before a2o_global_end
  ! - call contrail after lggp to have up to date input data and after trexp to
  !   have up to date lagrangian emission data
  ! - call airtrac_global_end and contrail_global_end after
  !   trexp_global_end, because new tendencies from trexp are used
  ! - vaxtra must be called last to have all channel objects up to date
  !   for transformation (for output)
  ! - aveout must be called last to have all other channel objects up to date
  ! - DOMINT must be called last to have all other channel objects up to date
  ! - clams_global_end must be called before OFFEMIS
  IF (USE_RNDTEST)     CALL     rndtest_global_end
  IF (USE_AIRTRAF)     CALL     airtraf_global_end
  IF (USE_ATTILA)      CALL      attila_global_end    ! lagrange !
  IF (USE_CLAMS)       CALL       clams_global_end(1)
  IF (USE_CLAMSSEDI)   CALL   clamssedi_global_end
  IF (USE_CLAMSTRAJ)   CALL   clamstraj_global_end
  IF (USE_CLAMSCIRRUS) CALL clamscirrus_global_end
  IF (USE_CLAMSCHEM)   CALL   clamschem_global_end
  IF (USE_CLAMSCHEME5) CALL clamscheme5_global_end
  IF (USE_CLAMSMIX)    CALL    clamsmix_global_end
  IF (USE_CLAMSBMIX)   CALL   clamsbmix_global_end
  IF (USE_CLAMSRDFRC)  CALL  clamsrdfrc_global_end
  IF (USE_CLAMSTRACER) CALL clamstracer_global_end
  IF (USE_CLAMS)       CALL       clams_global_end(2)
  IF (USE_DRADON)      CALL      dradon_global_end    ! lagrange !
  IF (USE_E5VDIFF)     CALL     e5vdiff_global_end
  IF (USE_OFFEMIS)     CALL     offemis_global_end    ! lagrange !
  IF (USE_ONEMIS)      CALL      onemis_global_end    ! lagrange !
  IF (USE_AIRSEA)      CALL      airsea_global_end    ! lagrange !
  IF (USE_BIOBURN)     CALL     bioburn_global_end
  IF (USE_MEGAN)       CALL       megan_global_end
  IF (USE_GWAVE)       CALL       gwave_global_end
  IF (USE_CVTRANS)     CALL     cvtrans_global_end    ! lagrange !
!  IF (USE_SCAV)        CALL        scav_global_end(1) ! lagrange !
  IF (USE_SEDI)        CALL        sedi_global_end    ! lagrange !
  IF (USE_CLOUDJ)      CALL      cloudj_global_end    ! lagrange !
  IF (USE_JVAL)        CALL        jval_global_end    ! lagrange !
  IF (USE_JVST)        CALL        jvst_global_end    ! lagrange !
  IF (USE_LNOX)        CALL        lnox_global_end    ! lagrange !
  IF (USE_MECCA)       CALL       mecca_global_end    ! lagrange !
  IF (USE_DDEP)        CALL        ddep_global_end    ! lagrange !
  IF (USE_DISSOC)      CALL      dissoc_global_end
  IF (USE_H2O)         CALL         h2o_global_end    ! lagrange !
  IF (USE_D14CO)       CALL       d14co_global_end    ! lagrange !
!  IF (USE_SCAV)        CALL        scav_global_end(2) ! lagrange !
! IF (USE_LGGP)        CALL        lggp_global_end    ! lagrange !
! IF (USE_LGTMIX)      CALL      lgtmix_global_end    ! lagrange !
! IF (USE_CH4)         CALL         ch4_global_end    ! lagrange !
! IF (USE_LGVFLUX)     CALL     lgvflux_global_end    ! lagrange !
  IF (USE_SURFACE)     CALL     surface_global_end(1) ! surf + lake + licetemp
  IF (USE_MLOCEAN)     CALL     mlocean_global_end
  IF (USE_SURFACE)     CALL     surface_global_end(2) ! sicetemp
  IF (USE_HD)          CALL          hd_global_end
  IF (USE_A2O)         CALL         a2o_global_end
  IF (USE_MPIOM)       CALL       mpiom_global_end
  IF (USE_HAMOCC)      CALL      hamocc_global_end
  IF (USE_OTPHYSC)     CALL     otphysc_global_end
  IF (USE_TNUDGE)      CALL      tnudge_global_end    ! lagrange !
  IF (USE_TREXP)       CALL       trexp_global_end    ! lagrange !
  IF (USE_VERTEX)      CALL      vertex_global_end
  IF (USE_AIRTRAC)     CALL     airtrac_global_end    ! lagrange !
  IF (USE_CONTRAIL)    CALL    contrail_global_end    ! lagrange !
  IF (USE_LGGP)        CALL        lggp_global_end    ! lagrange !
  IF (USE_LGTMIX)      CALL      lgtmix_global_end    ! lagrange !
  IF (USE_CH4)         CALL         ch4_global_end    ! lagrange !
  IF (USE_LGVFLUX)     CALL     lgvflux_global_end    ! lagrange !
  IF (USE_MMD2WAY)     CALL     mmd2way_global_end(1)
  IF (USE_MMD2WAY)     CALL     mmd2way_global_end(3)
  IF (USE_MMD2WAY)     CALL     mmd2way_global_end(2)
  IF (USE_SCALC)       CALL       scalc_global_end
!!$  IF (USE_CLAMS)       CALL       clams_global_end
!!$  IF (USE_CLAMSSEDI)   CALL   clamssedi_global_end
!!$  IF (USE_CLAMSTRAJ)   CALL   clamstraj_global_end
!!$  IF (USE_CLAMSCIRRUS) CALL clamscirrus_global_end
!!$  IF (USE_CLAMSCHEM)   CALL   clamschem_global_end
!!$  IF (USE_CLAMSCHEME5) CALL clamscheme5_global_end
!!$  IF (USE_CLAMSMIX)    CALL    clamsmix_global_end
!!$  IF (USE_CLAMSBMIX)   CALL   clamsbmix_global_end
!!$  IF (USE_CLAMSRDFRC)  CALL  clamsrdfrc_global_end
!!$  IF (USE_CLAMSTRACER) CALL clamstracer_global_end
  IF (USE_VAXTRA)      CALL      vaxtra_global_end
  IF (USE_AVEOUT)      CALL      aveout_global_end
  IF (USE_DOMINT)      CALL      domint_global_end

  CALL main_tracer_global_end
#ifdef MESSYTENDENCY
  CALL main_tendency_global_end
#endif
  CALL main_qtimer_global_end

  CALL main_control_unset_context

END SUBROUTINE messy_global_end

!==============================================================================

SUBROUTINE messy_write_output

  USE messy_main_switch       ! ONLY: USE_*
  USE messy_main_control      ! ONLY: entrypoint, subentry
  USE messy_main_tracer_bi,   ONLY: main_tracer_write_output
  USE messy_main_channel_bi,  ONLY: main_channel_write_output
  USE messy_main_rnd_bi,      ONLY: main_rnd_write_output

  USE messy_h2oiso_e5,        ONLY:  h2oiso_write_output
  USE messy_s4d_si,           ONLY:     s4d_write_output
  USE messy_scout_si,         ONLY:   scout_write_output
  USE messy_sorbit_si,        ONLY:  sorbit_write_output
  USE messy_timepos_e5,       ONLY: timepos_write_output
  USE messy_oasis3mct_si,     ONLY: oasis3mct_write_output

  CALL main_control_set_context(MEP_WRITE_OUTPUT)
  CALL main_tracer_write_output(1)

  IF (USE_H2OISO)  CALL  h2oiso_write_output
  IF (USE_SCOUT)   CALL   scout_write_output
  IF (USE_S4D)     CALL     s4d_write_output    ! force output & new file
  IF (USE_SORBIT)  CALL  sorbit_write_output    ! force output
  IF (USE_TIMEPOS) CALL timepos_write_output

  CALL main_rnd_write_output     ! save integer state vectors in channel objs
  CALL main_channel_write_output
  CALL main_control_unset_context

  ! added for ECHAM-FESOM2 coupling
  IF (USE_OASIS3MCT) CALL oasis3mct_write_output(I_WRITE_OUTPUT, I_END)

END SUBROUTINE messy_write_output

!==============================================================================

SUBROUTINE messy_write_restart

  USE messy_main_switch
  USE messy_main_control       ! ONLY: entrypoint, subentry
  USE messy_main_channel_bi,   ONLY: main_channel_write_restart

  CALL main_control_set_context(MEP_WRITE_RESTART)
  CALL main_channel_write_restart
  CALL main_control_unset_context

END SUBROUTINE messy_write_restart

!==============================================================================

SUBROUTINE messy_free_memory

  ! free memory (called by free_memory in mo_memory_streams)

  USE messy_main_switch   ! ONLY: USE_*
  USE messy_main_control  ! ONLY: entrypoint, subentry
  USE messy_main_qtimer_bi,    ONLY: main_qtimer_free_memory
  USE messy_main_tracer_bi,    ONLY: main_tracer_free_memory
  USE messy_main_channel_bi,   ONLY: main_channel_free_memory
  USE messy_main_data_bi,      ONLY: main_data_free_memory
  USE messy_main_import_bi,    ONLY: main_import_free_memory
#ifdef MESSYTENDENCY
  USE messy_main_tendency_bi,  ONLY: main_tendency_free_memory
#endif
  USE messy_main_grid_bi,      ONLY: main_grid_free_memory
  USE messy_main_timer_bi,     ONLY: main_timer_free_memory

  USE messy_a2o_e5,         ONLY:         a2o_free_memory
  USE messy_aeropt_si,      ONLY:      aeropt_free_memory
  USE messy_airsea_si,      ONLY:      airsea_free_memory
  USE messy_airtrac_e5,     ONLY:     airtrac_free_memory
  USE messy_airtraf_e5,     ONLY:     airtraf_free_memory
  USE messy_attila_e5,      ONLY:      attila_free_memory
  USE messy_bioburn_si,     ONLY:     bioburn_free_memory
  USE messy_cat_e5,         ONLY:         cat_free_memory
  USE messy_ch4_si,         ONLY:         ch4_free_memory
  USE messy_clams_si,       ONLY:       clams_free_memory
  USE messy_clamstraj_si,   ONLY:   clamstraj_free_memory
  USE messy_clamschem_si,   ONLY:   clamschem_free_memory
  USE messy_clamscheme5_si, ONLY: clamscheme5_free_memory
  USE messy_clamssedi_si,   ONLY:   clamssedi_free_memory
  USE messy_clamsmix_si,    ONLY:    clamsmix_free_memory
  USE messy_clamsbmix_si,   ONLY:   clamsbmix_free_memory
  USE messy_clamscirrus_si, ONLY: clamscirrus_free_memory
  USE messy_clamstracer_si, ONLY: clamstracer_free_memory
  USE messy_cloud_si,       ONLY:       cloud_free_memory
  USE messy_cloudj_si,      ONLY:      cloudj_free_memory
  USE messy_cloudopt_si,    ONLY:    cloudopt_free_memory
  USE messy_convect_si,     ONLY:     convect_free_memory
  USE messy_crm_e5,         ONLY:         crm_free_memory
  USE messy_cvtrans_si,     ONLY:     cvtrans_free_memory
  USE messy_d14co_si,       ONLY:       d14co_free_memory
  USE messy_ddep_si,        ONLY:        ddep_free_memory
  USE messy_dissoc_si,      ONLY:      dissoc_free_memory
  USE messy_diumod_si,      ONLY:      diumod_free_memory
  USE messy_domint_si,      ONLY:      domint_free_memory
  USE messy_e5vdiff_si,     ONLY:     e5vdiff_free_memory
  USE messy_e4chem_si,      ONLY:      e4chem_free_memory
  USE messy_ever_si,        ONLY:        ever_free_memory
  USE messy_gec_si,         ONLY:         gec_free_memory
  USE messy_gmxe_si,        ONLY:        gmxe_free_memory
  USE messy_gwave_si,       ONLY:       gwave_free_memory
  USE messy_h2o_e5,         ONLY:         h2o_free_memory
  USE messy_h2oiso_e5,      ONLY:      h2oiso_free_memory
  USE messy_hamocc_e5,      ONLY:      hamocc_free_memory
  USE messy_hd_e5,          ONLY:          hd_free_memory
  USE messy_isopcor_si,     ONLY:     isopcor_free_memory
  USE messy_jval_si,        ONLY:        jval_free_memory
  USE messy_jvst_si,        ONLY:        jvst_free_memory
  USE messy_lnox_si,        ONLY:        lnox_free_memory
  USE messy_lgtmix_e5,      ONLY:      lgtmix_free_memory
  USE messy_lgvflux_e5,     ONLY:     lgvflux_free_memory
  USE messy_m7_si,          ONLY:          m7_free_memory
  USE messy_made_si,        ONLY:        made_free_memory
  USE messy_made3_si,       ONLY:       made3_free_memory
  USE messy_mecca_si,       ONLY:       mecca_free_memory
  USE messy_megan_si,       ONLY:       megan_free_memory
  USE messy_mmd2way_si,     ONLY:     mmd2way_free_memory
  USE messy_mtskip_si,      ONLY:      mtskip_free_memory
  USE messy_nan_si,         ONLY:         nan_free_memory
  USE messy_o3orig_si,      ONLY:      o3orig_free_memory
  USE messy_offemis_si,     ONLY:     offemis_free_memory
  USE messy_onemis_si,      ONLY:      onemis_free_memory
  USE messy_oracle_si,      ONLY:      oracle_free_memory
  USE messy_orbit_si,       ONLY:       orbit_free_memory
  USE messy_otphysc_e5,     ONLY:     otphysc_free_memory
  USE messy_mpiom_e5,       ONLY:       mpiom_free_memory
  USE messy_plumegas_si,    ONLY:    plumegas_free_memory
  USE messy_qbo_si,         ONLY:         qbo_free_memory
  USE messy_rndtest_si,     ONLY:     rndtest_free_memory
  USE messy_satsims_e5,     ONLY:     satsims_free_memory
  USE messy_scalc_si,       ONLY:       scalc_free_memory
  USE messy_scav_si,        ONLY:        scav_free_memory
  USE messy_scout_si,       ONLY:       scout_free_memory
  USE messy_sedi_si,        ONLY:        sedi_free_memory
  USE messy_s4d_si,         ONLY:         s4d_free_memory
  USE messy_sorbit_si,      ONLY:      sorbit_free_memory
  USE messy_spe_e5,         ONLY:         spe_free_memory
  USE messy_spacenox_e5,    ONLY:    spacenox_free_memory
  USE messy_svoc_si,        ONLY:        svoc_free_memory
  USE messy_tagging_si,     ONLY:     tagging_free_memory
  USE messy_timepos_e5,     ONLY:     timepos_free_memory
  USE messy_tnudge_si,      ONLY:      tnudge_free_memory
  USE messy_tpulse_si,      ONLY:      tpulse_free_memory
  USE messy_trexp_si,       ONLY:       trexp_free_memory
  USE messy_trsync_si,      ONLY:      trsync_free_memory
  USE messy_ubcnox_e5,      ONLY:      ubcnox_free_memory
  USE messy_vaxtra_si,      ONLY:      vaxtra_free_memory
  USE messy_vertex_si,      ONLY:      vertex_free_memory
  USE messy_viso_si,        ONLY:        viso_free_memory

  IMPLICIT NONE

  CALL main_control_set_context(MEP_FREE_MEMORY)
  ! special requirements:
  ! - (none, all CALLs are in alphabetical order)
  IF (USE_A2O)         CALL         a2o_free_memory
  IF (USE_AEROPT)      CALL      aeropt_free_memory
  IF (USE_AIRSEA)      CALL      airsea_free_memory
  IF (USE_BIOBURN)     CALL     bioburn_free_memory
  IF (USE_AIRTRAC)     CALL     airtrac_free_memory
  IF (USE_AIRTRAF)     CALL     airtraf_free_memory
  IF (USE_ATTILA)      CALL      attila_free_memory
  IF (USE_CAT)         CALL         cat_free_memory
  IF (USE_CH4)         CALL         ch4_free_memory
  IF (USE_CLAMSMIX)    CALL    clamsmix_free_memory
  IF (USE_CLAMSBMIX)   CALL   clamsbmix_free_memory
  IF (USE_CLAMSCIRRUS) CALL clamscirrus_free_memory
  IF (USE_CLAMSCHEM)   CALL   clamschem_free_memory
  IF (USE_CLAMSSEDI)   CALL   clamssedi_free_memory
  IF (USE_CLAMSTRAJ)   CALL   clamstraj_free_memory
  IF (USE_CLAMSTRACER) CALL clamstracer_free_memory
  IF (USE_CLAMS)       CALL       clams_free_memory
  IF (USE_CLAMSCHEME5) CALL clamscheme5_free_memory
  IF (USE_CLOUD)       CALL       cloud_free_memory
  IF (USE_CLOUDOPT)    CALL    cloudopt_free_memory
  IF (USE_CONVECT)     CALL     convect_free_memory
  IF (USE_CRM)         CALL         crm_free_memory
  IF (USE_CVTRANS)     CALL     cvtrans_free_memory
  IF (USE_D14CO)       CALL       d14co_free_memory
  IF (USE_DDEP)        CALL        ddep_free_memory
  IF (USE_DISSOC)      CALL      dissoc_free_memory
  IF (USE_DIUMOD)      CALL      diumod_free_memory
  IF (USE_DOMINT)      CALL      domint_free_memory
  IF (USE_E4CHEM)      CALL      e4chem_free_memory
  IF (USE_E5VDIFF)     CALL     e5vdiff_free_memory
  IF (USE_EVER)        CALL        ever_free_memory
  IF (USE_GEC)         CALL         gec_free_memory
  IF (USE_GMXE)        CALL        gmxe_free_memory
  IF (USE_GWAVE)       CALL       gwave_free_memory
  IF (USE_H2O)         CALL         h2o_free_memory
  IF (USE_H2OISO)      CALL      h2oiso_free_memory
  IF (USE_ISOPCOR)     CALL     isopcor_free_memory
  IF (USE_CLOUDJ)      CALL      cloudj_free_memory
  IF (USE_JVAL)        CALL        jval_free_memory
  IF (USE_JVST)        CALL        jvst_free_memory
  IF (USE_LGTMIX)      CALL      lgtmix_free_memory
  IF (USE_LGVFLUX)     CALL     lgvflux_free_memory
  IF (USE_LNOX)        CALL        lnox_free_memory
  IF (USE_M7)          CALL          m7_free_memory
  IF (USE_MADE)        CALL        made_free_memory
  IF (USE_MADE3)       CALL       made3_free_memory
  IF (USE_MECCA)       CALL       mecca_free_memory
  IF (USE_MEGAN)       CALL       megan_free_memory
  IF (USE_MMD2WAY)     CALL     mmd2way_free_memory
  IF (USE_MPIOM)       CALL       mpiom_free_memory
  IF (USE_MTSKIP)      CALL      mtskip_free_memory
  IF (USE_NAN)         CALL         nan_free_memory
  IF (USE_HD)          CALL          hd_free_memory
  IF (USE_HAMOCC)      CALL      hamocc_free_memory
  IF (USE_O3ORIG)      CALL      o3orig_free_memory
  IF (USE_OFFEMIS)     CALL     offemis_free_memory
  IF (USE_ONEMIS)      CALL      onemis_free_memory
  IF (USE_ORACLE)      CALL      oracle_free_memory
  IF (USE_ORBIT)       CALL       orbit_free_memory
  IF (USE_OTPHYSC)     CALL     otphysc_free_memory
  IF (USE_PLUMEGAS)    CALL    plumegas_free_memory
  IF (USE_QBO)         CALL         qbo_free_memory
  IF (USE_RNDTEST)     CALL     rndtest_free_memory
  IF (USE_SATSIMS)     CALL     satsims_free_memory
  IF (USE_SCALC)       CALL       scalc_free_memory
  IF (USE_SCAV)        CALL        scav_free_memory
  IF (USE_SCOUT)       CALL       scout_free_memory
  IF (USE_SEDI)        CALL        sedi_free_memory
  IF (USE_S4D)         CALL         s4d_free_memory
  IF (USE_SORBIT)      CALL      sorbit_free_memory
  IF (USE_SPE)         CALL         spe_free_memory
  IF (USE_SPACENOX)    CALL    spacenox_free_memory
  IF (USE_SVOC)        CALL        svoc_free_memory
  IF (USE_TAGGING)     CALL     tagging_free_memory
  IF (USE_TIMEPOS)     CALL     timepos_free_memory
  IF (USE_TNUDGE)      CALL      tnudge_free_memory
  IF (USE_TPULSE)      CALL      tpulse_free_memory
  IF (USE_TREXP)       CALL       trexp_free_memory
  IF (USE_TRSYNC)      CALL      trsync_free_memory
  IF (USE_UBCNOX)      CALL      ubcnox_free_memory
  IF (USE_VAXTRA)      CALL      vaxtra_free_memory
  IF (USE_VERTEX)      CALL      vertex_free_memory
  IF (USE_VISO)        CALL        viso_free_memory

  CALL main_import_free_memory
  CALL main_tracer_free_memory
  CALL main_channel_free_memory
  CALL main_data_free_memory
  CALL main_qtimer_free_memory
#ifdef MESSYTENDENCY
  CALL main_tendency_free_memory
#endif
  CALL main_grid_free_memory
  CALL main_timer_free_memory
  CALL main_control_unset_context

END SUBROUTINE messy_free_memory

!==============================================================================

 SUBROUTINE messy_finalize

   USE messy_main_switch
   USE messy_main_control  ! ONLY: entrypoint, subentry
   USE messy_main_blather_bi, ONLY: main_blather_finalize

   ! finalize OASIS coupled simulations
   USE messy_oasis3mct_si,  ONLY: oasis3mct_finalize

   CALL main_control_set_context(MEP_FINALIZE)
   ! LOCAL
   CALL oasis3mct_finalize(USE_OASIS3MCT)
   CALL main_blather_finalize

   CALL main_control_unset_context

 END SUBROUTINE messy_finalize

!==============================================================================

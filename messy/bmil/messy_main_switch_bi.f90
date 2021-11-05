! ***************************************************************************
MODULE messy_main_switch_bi
! ***************************************************************************

  ! MODULE FOR SWITCHING SUBMODEL ON/OFF
  ! Authors:
  !   Patrick Joeckel, MPICH, September 2002

  ! NOTE: TO ADD A NEW SWITCH LOOK FOR
  !       '### ADD HERE'
  !       AND ADD IT ALSO TO messy_main_switch.f90

  ! MESSy
  USE messy_main_switch

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: main_switch_setup

CONTAINS

  ! ------------------------------------------------------------------------
  SUBROUTINE main_switch_setup

    ! ECHAM5/MESSy
    USE messy_main_blather_bi,       ONLY: error_bi
    USE messy_main_mpi_bi,           ONLY: p_parallel_io, p_io, p_bcast
    USE messy_main_channel_error_bi, ONLY: channel_halt
    USE messy_main_tools,            ONLY: find_next_free_unit

    IMPLICIT NONE

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'main_switch_setup'
    INTEGER     :: iou    ! I/O unit
    INTEGER     :: status

    IF (p_parallel_io) THEN
       iou = find_next_free_unit(100,200)
       CALL messy_main_read_nml_ctrl(status, iou)
       IF (status /= 0) CALL error_bi(' ',substr)
    END IF

    ! BROADCAST RESULTS
    !
    CALL p_bcast(L_TIME_INFO, p_io)
    !
    CALL p_bcast(USE_A2O,     p_io)
    CALL p_bcast(USE_ACCF,    p_io)
    CALL p_bcast(USE_AEROPT,  p_io)
    CALL p_bcast(USE_AIRSEA,  p_io)
    CALL p_bcast(USE_AIRTRAC, p_io)
    CALL p_bcast(USE_AIRTRAF, p_io)
    CALL p_bcast(USE_ATTILA,  p_io)
    CALL p_bcast(USE_AVEOUT,  p_io)
    CALL p_bcast(USE_BIOBURN, p_io)
    CALL p_bcast(USE_BUFLY,   p_io)
    CALL p_bcast(USE_CAT,     p_io)
    CALL p_bcast(USE_CH4,     p_io)
    CALL p_bcast(USE_CHEMGLUE,p_io)
    CALL p_bcast(USE_CLOUD,   p_io)
    CALL p_bcast(USE_CLOUDJ,  p_io)
    CALL p_bcast(USE_CLOUDOPT,p_io)
    CALL p_bcast(USE_CONTRAIL, p_io)
    CALL p_bcast(USE_CONVECT, p_io)
    CALL p_bcast(USE_CRM,     p_io)
    CALL p_bcast(USE_COSMOTOY,p_io)
    CALL p_bcast(USE_CVTRANS, p_io)
    CALL p_bcast(USE_D14CO,   p_io)
    CALL p_bcast(USE_DDEP,    p_io)
    CALL p_bcast(USE_DISSOC,  p_io)
    CALL p_bcast(USE_DIUMOD,  p_io)
    CALL p_bcast(USE_DOMINT,  p_io)
    CALL p_bcast(USE_DRADON,  p_io)
    CALL p_bcast(USE_E4CHEM,  p_io)
    CALL p_bcast(USE_E5VDIFF, p_io)
    CALL p_bcast(USE_EC2COSMO,p_io)
    CALL p_bcast(USE_EDITH,   p_io)
    CALL p_bcast(USE_EVER,    p_io)
    CALL p_bcast(USE_GMXE,    p_io)
    CALL p_bcast(USE_GWAVE,   p_io)
    CALL p_bcast(USE_H2O,     p_io)
    CALL p_bcast(USE_H2OEMIS, p_io)
    CALL p_bcast(USE_H2OISO,  p_io)
    CALL p_bcast(USE_HD,      p_io)
    CALL p_bcast(USE_HAMOCC,  p_io)
    CALL p_bcast(USE_IONS,    p_io)
    CALL p_bcast(USE_ISOPCOR, p_io)
    CALL p_bcast(USE_JVAL,    p_io)
    CALL p_bcast(USE_JVST,    p_io)
    CALL p_bcast(USE_LGGP,    p_io)
    CALL p_bcast(USE_LGTMIX,  p_io)
    CALL p_bcast(USE_LGVFLUX, p_io)
    CALL p_bcast(USE_LNOX,    p_io)
    CALL p_bcast(USE_M7,      p_io)
    CALL p_bcast(USE_MADE,    p_io)
    CALL p_bcast(USE_MADE3,   p_io)
    CALL p_bcast(USE_MECCA,   p_io)
    CALL p_bcast(USE_MEGAN,   p_io)
    CALL p_bcast(USE_MESOENERGY, p_io)
    CALL p_bcast(USE_MLOCEAN, p_io)
    CALL p_bcast(USE_MMFORCE, p_io)
    CALL p_bcast(USE_MPIOM,   p_io)
    CALL p_bcast(USE_MSBM,    p_io)
    CALL p_bcast(USE_MTSKIP,  p_io)
    CALL p_bcast(USE_NAN,     p_io)
    CALL p_bcast(USE_O3ORIG,  p_io)
    CALL p_bcast(USE_OASIS3MCT, p_io)
    CALL p_bcast(USE_OFFEMIS, p_io)
    CALL p_bcast(USE_ONEMIS,  p_io)
    CALL p_bcast(USE_ORACLE,  p_io)
    CALL p_bcast(USE_ORBIT,   p_io)
    CALL p_bcast(USE_OROGW,   p_io)
    CALL p_bcast(USE_OTPHYSC, p_io)
    CALL p_bcast(USE_PLUMEGAS,p_io)
    CALL p_bcast(USE_PTRAC,   p_io)
    CALL p_bcast(USE_PTRACINI,p_io)
    CALL p_bcast(USE_QBO,     p_io)
    CALL p_bcast(USE_RAD    , p_io)
    CALL p_bcast(USE_RELAX,   p_io)
    CALL p_bcast(USE_RNDTEST, p_io)
    CALL p_bcast(USE_SATSIMS, p_io)
    CALL p_bcast(USE_SCALC,   p_io)
    CALL p_bcast(USE_SCAV,    p_io)
    CALL p_bcast(USE_SCOUT,   p_io)
    CALL p_bcast(USE_SEDI,    p_io)
    CALL p_bcast(USE_S4D,     p_io)
    CALL p_bcast(USE_SF6,     p_io)
    CALL p_bcast(USE_SORBIT,  p_io)
    CALL p_bcast(USE_SPE,     p_io)
    CALL p_bcast(USE_SPACENOX,p_io)
    CALL p_bcast(USE_SVOC,    p_io)
    CALL p_bcast(USE_GEC,     p_io)
    CALL p_bcast(USE_SURFACE, p_io)
    CALL p_bcast(USE_TAGGING, p_io)
    CALL p_bcast(USE_TBUDGET, p_io)
    CALL p_bcast(USE_TIMEPOS, p_io)
    CALL p_bcast(USE_TNUDGE,  p_io)
    CALL p_bcast(USE_TPULSE,  p_io)
    CALL p_bcast(USE_TREXP,   p_io)
    CALL p_bcast(USE_TROPOP,  p_io)
    CALL p_bcast(USE_TRSYNC,  p_io)
    CALL p_bcast(USE_UBCNOX,  p_io)
    CALL p_bcast(USE_VAHR,    p_io)
    CALL p_bcast(USE_VAXTRA,  p_io)
    CALL p_bcast(USE_VERTDIFF,p_io)
    CALL p_bcast(USE_VERTEX,  p_io)
    CALL p_bcast(USE_VISO,    p_io)
    CALL p_bcast(USE_VISOP,   p_io)
    CALL p_bcast(USE_MMD2WAY, p_io)
    ! ### example submodels
    CALL p_bcast(USE_SUBMOD1, p_io)
    CALL p_bcast(USE_SUBMOD2, p_io)
    CALL p_bcast(USE_SUBMOD3, p_io)
    CALL p_bcast(USE_TESTEVENT,p_io)

    CALL p_bcast(USE_CLAMS,       p_io)
    CALL p_bcast(USE_CLAMSTRAJ,   p_io)
    CALL p_bcast(USE_CLAMSCHEM,   p_io)
    CALL p_bcast(USE_CLAMSCHEME5, p_io)
    CALL p_bcast(USE_CLAMSSEDI,   p_io)
    CALL p_bcast(USE_CLAMSMIX ,   p_io)
    CALL p_bcast(USE_CLAMSBMIX,   p_io)
    CALL p_bcast(USE_CLAMSCIRRUS, p_io)
    CALL p_bcast(USE_CLAMSRDFRC,  p_io)
    CALL p_bcast(USE_CLAMSTRACER, p_io)
    CALL p_bcast(USE_CLAMSDEEPCONV, p_io)
! ju_ch_20110429-

    CALL p_bcast(USE_MXL,     p_io)

    ! ### ADD HERE

    CALL switch_init(status)
    CALL channel_halt(substr, status)

  END SUBROUTINE main_switch_setup
  ! ------------------------------------------------------------------------

  ! -------------------------------------------------------------
  SUBROUTINE switch_init(status)

    ! MESSy
    USE messy_main_channel,       ONLY: new_attribute
    USE messy_main_constants_mem, ONLY: modstr_MESSy => modstr &
                                      , modver_MESSy => modver
    USE messy_main_tracer,        ONLY: modstr_tracer => modstr &
                                      , modver_tracer => modver
    USE messy_main_channel,       ONLY: modstr_channel => modstr &
                                      , modver_channel => modver
    USE messy_main_timer,         ONLY: modstr_timer => modstr &
                                      , modver_timer => modver
    USE messy_main_qtimer,        ONLY: modstr_qtimer => modstr &
                                      , modver_qtimer => modver
    USE messy_main_import,        ONLY: modstr_import => modstr &
                                      , modver_import => modver
    USE messy_main_rnd,           ONLY: modstr_rnd => modstr &
                                      , modver_rnd => modver
#ifdef MESSYTENDENCY
    USE messy_main_tendency,      ONLY: modstr_tendency => modstr &
                                      , modver_tendency => modver
#endif

#if defined(ECHAM5) || defined(COSMO) || defined(CESM1)
    ! MESSY SUBMODELS
    USE messy_a2o,        ONLY: modstr_a2o=>modstr,     modver_a2o=>modver
    USE messy_accf,       ONLY: modstr_accf=>modstr,    modver_accf=>modver
    USE messy_aeropt,     ONLY: modstr_aeropt=>modstr,  modver_aeropt=>modver
    USE messy_airsea,     ONLY: modstr_airsea=>modstr,  modver_airsea=>modver
    USE messy_airtrac,    ONLY: modstr_airtrac=>modstr, modver_airtrac=>modver
    USE messy_airtraf,    ONLY: modstr_airtraf=>modstr, modver_airtraf=>modver
    USE messy_attila,     ONLY: modstr_attila=>modstr,  modver_attila=>modver
    USE messy_aveout,     ONLY: modstr_aveout=>modstr,  modver_aveout=>modver
    USE messy_bioburn,    ONLY: modstr_bioburn=>modstr, modver_bioburn=>modver
    USE messy_bufly,      ONLY: modstr_bufly=>modstr,   modver_bufly=>modver
    USE messy_cat,        ONLY: modstr_cat=>modstr,     modver_cat=>modver
    USE messy_ch4,        ONLY: modstr_ch4=>modstr,     modver_ch4=>modver
    USE messy_chemglue,   ONLY: modstr_chemglue=>modstr,modver_chemglue=>modver
    USE messy_cloud,      ONLY: modstr_cloud=>modstr,   modver_cloud=>modver
    USE messy_cloudj,     ONLY: modstr_cloudj=>modstr,  modver_cloudj=>modver
    USE messy_cloudopt,   ONLY: modstr_cloudopt=>modstr,modver_cloudopt=>modver
    USE messy_contrail,   ONLY: modstr_contrail=>modstr,modver_contrail=>modver
    USE messy_convect,    ONLY: modstr_convect=>modstr, modver_convect=>modver
    USE messy_crm,        ONLY: modstr_crm=>modstr,     modver_crm=>modver
    USE messy_cosmotoy,   ONLY: modstr_cosmotoy=>modstr,modver_cosmotoy=>modver
    USE messy_cvtrans,    ONLY: modstr_cvtrans=>modstr, modver_cvtrans=>modver
    USE messy_d14co,      ONLY: modstr_d14co=>modstr,   modver_d14co=>modver
    USE messy_ddep,       ONLY: modstr_ddep=>modstr,    modver_ddep=>modver
    USE messy_dissoc,     ONLY: modstr_dissoc=>modstr,  modver_dissoc=>modver
    USE messy_diumod,     ONLY: modstr_diumod=>modstr,  modver_diumod=>modver
    USE messy_domint,     ONLY: modstr_domint=>modstr,  modver_domint=>modver
    USE messy_dradon,     ONLY: modstr_dradon=>modstr,  modver_dradon=>modver
    USE messy_e4chem,     ONLY: modstr_e4chem=>modstr,  modver_e4chem=>modver
    USE messy_e5vdiff,    ONLY: modstr_e5vdiff=>modstr, modver_e5vdiff=>modver
    USE messy_ec2cosmo,   ONLY: modstr_ec2cosmo=>modstr,modver_ec2cosmo=>modver
    USE messy_edith,      ONLY: modstr_edith=>modstr,   modver_edith=>modver
    USE messy_ever,       ONLY: modstr_ever=>modstr,    modver_ever=>modver
    USE messy_gmxe,       ONLY: modstr_gmxe=>modstr,    modver_gmxe=>modver
    USE messy_gwave,      ONLY: modstr_gwave=>modstr,   modver_gwave=>modver
    USE messy_h2o,        ONLY: modstr_h2o=>modstr,     modver_h2o=>modver
    USE messy_h2oemis,    ONLY: modstr_h2oemis=>modstr, modver_h2oemis=>modver
    USE messy_h2oiso,     ONLY: modstr_h2oiso=>modstr,  modver_h2oiso=>modver
    USE messy_hamocc,     ONLY: modstr_hamocc=>modstr,  modver_hamocc=>modver
    USE messy_hd,         ONLY: modstr_hd=>modstr,      modver_hd=>modver
    USE messy_ions,       ONLY: modstr_ions=>modstr,    modver_ions=>modver
    USE messy_isopcor,    ONLY: modstr_isopcor=>modstr, modver_isopcor=>modver
    USE messy_jval,       ONLY: modstr_jval=>modstr,    modver_jval=>modver
    USE messy_jvst,       ONLY: modstr_jvst=>modstr,    modver_jvst=>modver
    USE messy_lggp,       ONLY: modstr_lggp=>modstr,    modver_lggp=>modver
    USE messy_lgtmix,     ONLY: modstr_lgtmix=>modstr,  modver_lgtmix=>modver
    USE messy_lgvflux,    ONLY: modstr_lgvflux=>modstr, modver_lgvflux=>modver
    USE messy_lnox,       ONLY: modstr_lnox=>modstr,    modver_lnox=>modver
    USE messy_m7,         ONLY: modstr_m7=>modstr,      modver_m7=>modver
    USE messy_made,       ONLY: modstr_made=>modstr,    modver_made=>modver
    USE messy_made3,      ONLY: modstr_made3=>modstr,   modver_made3=>modver
    USE messy_mecca,      ONLY: modstr_mecca=>modstr,   modver_mecca=>modver
    USE messy_megan,      ONLY: modstr_megan=>modstr,   modver_megan=>modver
    USE messy_mesoenergy, ONLY: modstr_mesoenergy=>modstr, modver_mesoenergy=>modver
    USE messy_mlocean,    ONLY: modstr_mlocean=>modstr, modver_mlocean=>modver
    USE messy_mmforce,    ONLY: modstr_mmforce=>modstr, modver_mmforce=>modver
    USE messy_mpiom,      ONLY: modstr_mpiom=>modstr,   modver_mpiom=>modver
    USE messy_msbm,       ONLY: modstr_msbm=>modstr,    modver_msbm=>modver
    USE messy_mtskip,     ONLY: modstr_mtskip=>modstr,  modver_mtskip=>modver
    USE messy_nan,        ONLY: modstr_nan=>modstr,     modver_nan=>modver
    USE messy_o3orig,     ONLY: modstr_o3orig=>modstr,  modver_o3orig=>modver
    USE messy_oasis3mct,  ONLY: modstr_oasis3mct=>modstr,modver_oasis3mct=>modver
    USE messy_offemis,    ONLY: modstr_offemis=>modstr, modver_offemis=>modver
    USE messy_onemis,     ONLY: modstr_onemis=>modstr,  modver_onemis=>modver
    USE messy_oracle,     ONLY: modstr_oracle=>modstr,  modver_oracle=>modver
    USE messy_orbit,      ONLY: modstr_orbit=>modstr,   modver_orbit=>modver
    USE messy_orogw,      ONLY: modstr_orogw=>modstr,   modver_orogw=>modver
    USE messy_otphysc,    ONLY: modstr_otphysc=>modstr, modver_otphysc=>modver
    USE messy_plumegas,   ONLY: modstr_plumegas=>modstr,modver_plumegas=>modver
    USE messy_ptrac,      ONLY: modstr_ptrac=>modstr,   modver_ptrac=>modver
    USE messy_ptracini,   ONLY: modstr_ptracini=>modstr,modver_ptracini=>modver
    USE messy_qbo,        ONLY: modstr_qbo=>modstr,     modver_qbo=>modver
    USE messy_rad,        ONLY: modstr_rad=>modstr,     modver_rad=>modver
    USE messy_relax,      ONLY: modstr_relax=>modstr,   modver_relax=>modver
    USE messy_rndtest,    ONLY: modstr_rndtest=>modstr, modver_rndtest=>modver
    USE messy_satsims,    ONLY: modstr_satsims=>modstr, modver_satsims=>modver
    USE messy_scalc,      ONLY: modstr_scalc=>modstr,   modver_scalc=>modver
    USE messy_scav,       ONLY: modstr_scav=>modstr,    modver_scav=>modver
    USE messy_scout,      ONLY: modstr_scout=>modstr,   modver_scout=>modver
    USE messy_sedi,       ONLY: modstr_sedi=>modstr,    modver_sedi=>modver
    USE messy_s4d,        ONLY: modstr_s4d=>modstr,     modver_s4d=>modver
    USE messy_sf6,        ONLY: modstr_sf6=>modstr,     modver_sf6=>modver
    USE messy_sorbit,     ONLY: modstr_sorbit=>modstr,  modver_sorbit=>modver
    USE messy_spacenox,   ONLY: modstr_spacenox=>modstr,modver_spacenox=>modver
    USE messy_svoc,       ONLY: modstr_svoc=>modstr,    modver_svoc=>modver
    USE messy_gec,        ONLY: modstr_gec=>modstr,     modver_gec=>modver
    USE messy_spe,        ONLY: modstr_spe=>modstr,     modver_spe=>modver
    USE messy_surface,    ONLY: modstr_surface=>modstr, modver_surface=>modver
    USE messy_tagging,    ONLY: modstr_tagging=>modstr, modver_tagging=>modver
    USE messy_tbudget,    ONLY: modstr_tbudget=>modstr, modver_tbudget=>modver
    USE messy_timepos,    ONLY: modstr_timepos=>modstr, modver_timepos=>modver
    USE messy_tnudge,     ONLY: modstr_tnudge=>modstr,  modver_tnudge=>modver
    USE messy_trexp,      ONLY: modstr_trexp=>modstr,   modver_trexp=>modver
    USE messy_tpulse,     ONLY: modstr_tpulse=>modstr,  modver_tpulse=>modver
    USE messy_tropop,     ONLY: modstr_tropop=>modstr,  modver_tropop=>modver
    USE messy_trsync,     ONLY: modstr_trsync=>modstr,  modver_trsync=>modver
    USE messy_ubcnox,     ONLY: modstr_ubcnox=>modstr,  modver_ubcnox=>modver
    USE messy_vahr,       ONLY: modstr_vahr=>modstr,    modver_vahr=>modver
    USE messy_vaxtra,     ONLY: modstr_vaxtra=>modstr,  modver_vaxtra=>modver
    USE messy_vertdiff,   ONLY: modstr_vertdiff=>modstr,modver_vertdiff=>modver
    USE messy_vertex,     ONLY: modstr_vertex=>modstr,  modver_vertex=>modver
    USE messy_viso,       ONLY: modstr_viso=>modstr,    modver_viso=>modver
    USE messy_visop,      ONLY: modstr_visop=>modstr,   modver_visop=>modver
    USE messy_mmd2way,    ONLY: modstr_mmd2way=>modstr, modver_mmd2way=>modver

    USE messy_clams,       ONLY: modstr_clams=>modstr,       modver_clams=>modver
    USE messy_clamstraj,   ONLY: modstr_clamstraj=>modstr,   modver_clamstraj=>modver
    USE messy_clamschem,   ONLY: modstr_clamschem=>modstr,   modver_clamschem=>modver
    USE messy_clamscheme5, ONLY: modstr_clamscheme5=>modstr, modver_clamscheme5=>modver
    USE messy_clamssedi,   ONLY: modstr_clamssedi=>modstr,   modver_clamssedi=>modver
    USE messy_clamsmix,    ONLY: modstr_clamsmix=>modstr,    modver_clamsmix=>modver
    USE messy_clamsbmix,   ONLY: modstr_clamsbmix=>modstr,   modver_clamsbmix=>modver
    USE messy_clamscirrus, ONLY: modstr_clamscirrus=>modstr, modver_clamscirrus=>modver
    USE messy_clamsrdfrc,  ONLY: modstr_clamsrdfrc=>modstr,  modver_clamsrdfrc=>modver
    USE messy_clamstracer, ONLY: modstr_clamstracer=>modstr, modver_clamstracer=>modver
    USE messy_clamsdeepconv, ONLY: modstr_clamsdeepconv=>modstr, modver_clamsdeepconv=>modver
! ju_ch_20110429-
#endif

#ifdef ICON
!!$    USE messy_tropop,      ONLY: modstr_tropop=>modstr,      modver_tropop=>modver
    USE messy_ptrac,      ONLY: modstr_ptrac=>modstr,   modver_ptrac=>modver
#endif

#ifdef MBM_MPIOM
    USE messy_mpiom,      ONLY: modstr_mpiom=>modstr,   modver_mpiom=>modver
#endif

#ifdef BLANK
#if defined (MBM_BLANK)
    USE messy_ptrac,     ONLY: modstr_ptrac=>modstr,   modver_ptrac=>modver
    USE messy_orbit,     ONLY: modstr_orbit=>modstr,   modver_orbit=>modver
    USE messy_submod1,   ONLY: modstr_submod1=>modstr, modver_submod1=>modver
    USE messy_submod2,   ONLY: modstr_submod2=>modstr, modver_submod2=>modver
    USE messy_submod3,   ONLY: modstr_submod3=>modstr, modver_submod3=>modver
    USE messy_testevent, ONLY: modstr_testevent=>modstr, modver_testevent=>modver
#endif
#if defined(MBM_QBO)
    USE messy_qbo,        ONLY: modstr_qbo=>modstr,     modver_qbo=>modver
#endif
#if defined(MBM_DISSOC)
    USE messy_dissoc,     ONLY: modstr_dissoc=>modstr,  modver_dissoc=>modver
#endif
#if defined(MBM_RAD)
    USE messy_aeropt,     ONLY: modstr_aeropt=>modstr,  modver_aeropt=>modver
    USE messy_cloudopt,   ONLY: modstr_cloudopt=>modstr,modver_cloudopt=>modver
    USE messy_orbit,      ONLY: modstr_orbit=>modstr,   modver_orbit=>modver
    USE messy_rad,        ONLY: modstr_rad=>modstr,     modver_rad=>modver
#endif
#endif

#ifdef MBM_CLAMS
    USE messy_clams,       ONLY: modstr_clams=>modstr,       modver_clams=>modver
    USE messy_clamstraj,   ONLY: modstr_clamstraj=>modstr,   modver_clamstraj=>modver
    USE messy_clamschem,   ONLY: modstr_clamschem=>modstr,   modver_clamschem=>modver
!    USE messy_clamscheme5, ONLY: modstr_clamscheme5=>modstr, modver_clamscheme5=>modver
    USE messy_clamssedi,   ONLY: modstr_clamssedi=>modstr,   modver_clamssedi=>modver
    USE messy_clamsmix,    ONLY: modstr_clamsmix=>modstr,    modver_clamsmix=>modver
    USE messy_clamsbmix,   ONLY: modstr_clamsbmix=>modstr,   modver_clamsbmix=>modver
    USE messy_clamscirrus, ONLY: modstr_clamscirrus=>modstr, modver_clamscirrus=>modver
!    USE messy_clamsrdfrc,  ONLY: modstr_clamsrdfrc=>modstr,  modver_clamsrdfrc=>modver
    USE messy_clamstracer, ONLY: modstr_clamstracer=>modstr, modver_clamstracer=>modver
    USE messy_clamsdeepconv, ONLY: modstr_clamsdeepconv=>modstr, modver_clamsdeepconv=>modver
#endif

#if defined(VERTICO)
    USE messy_ptrac,      ONLY: modstr_ptrac=>modstr,   modver_ptrac=>modver
    USE messy_ddep,       ONLY: modstr_ddep=>modstr,    modver_ddep=>modver
    USE messy_mxl,        ONLY: modstr_mxl=>modstr,     modver_mxl=>modver
    USE messy_mecca,      ONLY: modstr_mecca=>modstr,   modver_mecca=>modver
    USE messy_megan,      ONLY: modstr_megan=>modstr,   modver_megan=>modver
    USE messy_offemis,    ONLY: modstr_offemis=>modstr, modver_offemis=>modver
    USE messy_onemis,     ONLY: modstr_onemis=>modstr,  modver_onemis=>modver
    USE messy_jval,       ONLY: modstr_jval=>modstr,    modver_jval=>modver
#endif

    ! ### ADD HERE

    IMPLICIT NONE

    ! I/O
    INTEGER, INTENT(OUT) :: status

    ! put MESSy information to output file
    CALL new_attribute(status, 'MESSy' &
         , c = modstr_MESSy//' version '//modver_MESSy//&
         &', http://www.messy-interface.org' )
    IF (status /= 0) RETURN
    !
    CALL put_submodel_att (status, .TRUE., modstr, modver) ! switch
    IF (status /= 0) RETURN
    CALL put_submodel_att (status, .TRUE., modstr_channel, modver_channel)
    IF (status /= 0) RETURN
    CALL put_submodel_att (status, .TRUE., modstr_tracer,  modver_tracer)
    IF (status /= 0) RETURN
    CALL put_submodel_att (status, .TRUE., modstr_timer,   modver_timer)
    IF (status /= 0) RETURN
    CALL put_submodel_att (status, .TRUE., modstr_qtimer,  modver_qtimer)
    IF (status /= 0) RETURN
    CALL put_submodel_att (status, .TRUE., modstr_import,  modver_import)
    IF (status /= 0) RETURN
    CALL put_submodel_att (status, .TRUE., modstr_rnd,     modver_rnd)
    IF (status /= 0) RETURN
#ifdef MESSYTENDENCY
    CALL put_submodel_att (status, .TRUE., modstr_tendency, modver_tendency)
    IF (status /= 0) RETURN
#endif

#if defined(ECHAM5) || defined(COSMO) || defined(CESM1)
    ! put submodel version into output file
    ! MESSy SUBMODELS
    CALL put_submodel_att (status, USE_A2O,      modstr_a2o,      modver_a2o)      ; IF (status/=0) RETURN
    CALL put_submodel_att (status, USE_ACCF,     modstr_accf,     modver_accf)     ; IF (status/=0) RETURN
    CALL put_submodel_att (status, USE_AEROPT,   modstr_aeropt,   modver_aeropt)   ; IF (status/=0) RETURN
    CALL put_submodel_att (status, USE_AIRSEA,   modstr_airsea,   modver_airsea)   ; IF (status/=0) RETURN
    CALL put_submodel_att (status, USE_AIRTRAC,  modstr_airtrac,  modver_airtrac)  ; IF (status/=0) RETURN
    CALL put_submodel_att (status, USE_AIRTRAF,  modstr_airtraf,  modver_airtraf)   ; IF (status/=0) RETURN
    CALL put_submodel_att (status, USE_ATTILA,   modstr_attila,   modver_attila)   ; IF (status/=0) RETURN
    CALL put_submodel_att (status, USE_AVEOUT,   modstr_aveout,   modver_aveout)   ; IF (status/=0) RETURN
    CALL put_submodel_att (status, USE_BIOBURN,  modstr_bioburn,  modver_bioburn)  ; IF (status/=0) RETURN
    CALL put_submodel_att (status, USE_BUFLY,    modstr_bufly,    modver_bufly)    ; IF (status/=0) RETURN
    CALL put_submodel_att (status, USE_CAT,      modstr_cat,      modver_cat)      ; IF (status/=0) RETURN
    CALL put_submodel_att (status, USE_CH4,      modstr_ch4,      modver_ch4)      ; IF (status/=0) RETURN
    CALL put_submodel_att (status, USE_CHEMGLUE, modstr_chemglue, modver_chemglue) ; IF (status/=0) RETURN
    CALL put_submodel_att (status, USE_CLOUD,    modstr_cloud,    modver_cloud)    ; IF (status/=0) RETURN
    CALL put_submodel_att (status, USE_CLOUDJ,   modstr_cloudj,   modver_cloudj)   ; IF (status/=0) RETURN
    CALL put_submodel_att (status, USE_CLOUDOPT, modstr_cloudopt, modver_cloudopt) ; IF (status/=0) RETURN
    CALL put_submodel_att (status, USE_CONTRAIL, modstr_contrail, modver_contrail) ; IF (status/=0) RETURN
    CALL put_submodel_att (status, USE_CONVECT,  modstr_convect,  modver_convect)  ; IF (status/=0) RETURN
    CALL put_submodel_att (status, USE_CRM,      modstr_crm,      modver_crm)      ; IF (status/=0) RETURN
    CALL put_submodel_att (status, USE_COSMOTOY, modstr_cosmotoy, modver_cosmotoy) ; IF (status/=0) RETURN
    CALL put_submodel_att (status, USE_CVTRANS,  modstr_cvtrans,  modver_cvtrans)  ; IF (status/=0) RETURN
    CALL put_submodel_att (status, USE_D14CO,    modstr_d14co,    modver_d14co)    ; IF (status/=0) RETURN
    CALL put_submodel_att (status, USE_DDEP,     modstr_ddep,     modver_ddep)     ; IF (status/=0) RETURN
    CALL put_submodel_att (status, USE_DISSOC,   modstr_dissoc,   modver_dissoc)     ; IF (status/=0) RETURN
    CALL put_submodel_att (status, USE_DIUMOD,   modstr_diumod,   modver_diumod)     ; IF (status/=0) RETURN
    CALL put_submodel_att (status, USE_DOMINT,   modstr_domint,   modver_domint)   ; IF (status/=0) RETURN
    CALL put_submodel_att (status, USE_DRADON,   modstr_dradon,   modver_dradon)   ; IF (status/=0) RETURN
    CALL put_submodel_att (status, USE_EC2COSMO, modstr_ec2cosmo, modver_ec2cosmo) ; IF (status/=0) RETURN
    CALL put_submodel_att (status, USE_EDITH,    modstr_edith,    modver_edith)    ; IF (status/=0) RETURN
    CALL put_submodel_att (status, USE_EVER,     modstr_ever,     modver_ever)     ; IF (status/=0) RETURN
    CALL put_submodel_att (status, USE_GMXE,     modstr_gmxe,     modver_gmxe)     ; IF (status/=0) RETURN
    CALL put_submodel_att (status, USE_GWAVE,    modstr_gwave,    modver_gwave)    ; IF (status/=0) RETURN
    CALL put_submodel_att (status, USE_H2O,      modstr_h2o,      modver_h2o)      ; IF (status/=0) RETURN
    CALL put_submodel_att (status, USE_H2OEMIS,  modstr_h2oemis,  modver_h2oemis)  ; IF (status/=0) RETURN
    CALL put_submodel_att (status, USE_H2OISO,   modstr_h2oiso,   modver_h2oiso)   ; IF (status/=0) RETURN
    CALL put_submodel_att (status, USE_HAMOCC,   modstr_hamocc,   modver_hamocc)   ; IF (status/=0) RETURN
    CALL put_submodel_att (status, USE_HD,       modstr_hd,       modver_hd)       ; IF (status/=0) RETURN
    CALL put_submodel_att (status, USE_IONS,     modstr_ions,     modver_ions)  ; IF (status/=0) RETURN
    CALL put_submodel_att (status, USE_ISOPCOR,  modstr_isopcor,  modver_isopcor)  ; IF (status/=0) RETURN
    CALL put_submodel_att (status, USE_JVAL,     modstr_jval,     modver_jval)     ; IF (status/=0) RETURN
    CALL put_submodel_att (status, USE_jvst,     modstr_jvst,     modver_jvst)     ; IF (status/=0) RETURN
    CALL put_submodel_att (status, USE_LGGP,     modstr_lggp,     modver_lggp)     ; IF (status/=0) RETURN
    CALL put_submodel_att (status, USE_LGTMIX,   modstr_lgtmix,   modver_lgtmix)   ; IF (status/=0) RETURN
    CALL put_submodel_att (status, USE_LGVFLUX,  modstr_lgvflux,  modver_lgvflux)  ; IF (status/=0) RETURN
    CALL put_submodel_att (status, USE_LNOX,     modstr_lnox,     modver_lnox)     ; IF (status/=0) RETURN
    CALL put_submodel_att (status, USE_M7,       modstr_m7,       modver_m7)       ; IF (status/=0) RETURN
    CALL put_submodel_att (status, USE_MADE,     modstr_made,     modver_made)     ; IF (status/=0) RETURN
    CALL put_submodel_att (status, USE_MADE3,    modstr_made3,    modver_made3)    ; IF (status/=0) RETURN
    CALL put_submodel_att (status, USE_MECCA,    modstr_mecca,    modver_mecca)    ; IF (status/=0) RETURN
    CALL put_submodel_att (status, USE_MEGAN,    modstr_megan,    modver_megan)    ; IF (status/=0) RETURN
    CALL put_submodel_att (status, USE_MESOENERGY, modstr_mesoenergy, modver_mesoenergy) ; IF (status/=0) RETURN
    CALL put_submodel_att (status, USE_MLOCEAN,  modstr_mlocean,  modver_mlocean)  ; IF (status/=0) RETURN
    CALL put_submodel_att (status, USE_MMFORCE,  modstr_mmforce,  modver_mmforce)  ; IF (status/=0) RETURN
    CALL put_submodel_att (status, USE_MPIOM,    modstr_mpiom,    modver_mpiom)    ; IF (status/=0) RETURN
    CALL put_submodel_att (status, USE_MSBM,     modstr_msbm,     modver_msbm)     ; IF (status/=0) RETURN
    CALL put_submodel_att (status, USE_MTSKIP,   modstr_mtskip,   modver_mtskip)   ; IF (status/=0) RETURN
    CALL put_submodel_att (status, USE_NAN,      modstr_nan,      modver_nan)  ; IF (status/=0) RETURN
    CALL put_submodel_att (status, USE_O3ORIG,   modstr_o3orig,   modver_o3orig)   ; IF (status/=0) RETURN
    CALL put_submodel_att (status, USE_OASIS3MCT,modstr_oasis3mct,modver_oasis3mct)  ; IF (status/=0) RETURN
    CALL put_submodel_att (status, USE_OFFEMIS,  modstr_offemis,  modver_offemis)  ; IF (status/=0) RETURN
    CALL put_submodel_att (status, USE_ONEMIS,   modstr_onemis,   modver_onemis)   ; IF (status/=0) RETURN
    CALL put_submodel_att (status, USE_ORACLE,   modstr_oracle,   modver_oracle)   ; IF (status/=0) RETURN
    CALL put_submodel_att (status, USE_ORBIT,    modstr_orbit,    modver_orbit)    ; IF (status/=0) RETURN
    CALL put_submodel_att (status, USE_OROGW,    modstr_orogw,    modver_orogw)    ; IF (status/=0) RETURN
    CALL put_submodel_att (status, USE_OTPHYSC,  modstr_otphysc,  modver_otphysc)  ; IF (status/=0) RETURN
    CALL put_submodel_att (status, USE_PLUMEGAS, modstr_plumegas, modver_plumegas) ; IF (status/=0) RETURN
    CALL put_submodel_att (status, USE_PTRAC,    modstr_ptrac,    modver_ptrac)    ; IF (status/=0) RETURN
    CALL put_submodel_att (status, USE_PTRACINI, modstr_ptracini, modver_ptracini) ; IF (status/=0) RETURN
    CALL put_submodel_att (status, USE_QBO,      modstr_qbo,      modver_qbo)      ; IF (status/=0) RETURN
    CALL put_submodel_att (status, USE_RAD,      modstr_rad,      modver_rad)      ; IF (status/=0) RETURN
    CALL put_submodel_att (status, USE_RELAX,    modstr_relax,    modver_relax)    ; IF (status/=0) RETURN
    CALL put_submodel_att (status, USE_RNDTEST,  modstr_rndtest,  modver_rndtest)    ; IF (status/=0) RETURN
    CALL put_submodel_att (status, USE_SATSIMS,  modstr_satsims,  modver_satsims)  ; IF (status/=0) RETURN
    CALL put_submodel_att (status, USE_SCALC,    modstr_scalc,    modver_scalc)    ; IF (status/=0) RETURN
    CALL put_submodel_att (status, USE_SCAV,     modstr_scav,     modver_scav)     ; IF (status/=0) RETURN
    CALL put_submodel_att (status, USE_SCOUT,    modstr_scout,    modver_scout)    ; IF (status/=0) RETURN
    CALL put_submodel_att (status, USE_SEDI,     modstr_sedi,     modver_sedi)     ; IF (status/=0) RETURN
    CALL put_submodel_att (status, USE_S4D,      modstr_s4d,      modver_s4d)      ; IF (status/=0) RETURN
    CALL put_submodel_att (status, USE_SF6,      modstr_sf6,      modver_sf6)      ; IF (status/=0) RETURN
    CALL put_submodel_att (status, USE_SORBIT,   modstr_sorbit,   modver_sorbit)   ; IF (status/=0) RETURN
    CALL put_submodel_att (status, USE_SPACENOX, modstr_spacenox, modver_spacenox) ; IF (status/=0) RETURN
    CALL put_submodel_att (status, USE_SVOC,     modstr_svoc,     modver_svoc)     ; IF (status/=0) RETURN
    CALL put_submodel_att (status, USE_GEC,      modstr_gec,      modver_gec)      ; IF (status/=0) RETURN
    CALL put_submodel_att (status, USE_E4CHEM,   modstr_e4chem,   modver_e4chem)   ; IF (status/=0) RETURN
    CALL put_submodel_att (status, USE_E5VDIFF,  modstr_e5vdiff,  modver_e5vdiff)    ; IF (status/=0) RETURN
    CALL put_submodel_att (status, USE_SPE,      modstr_spe,      modver_spe)      ; IF (status/=0) RETURN
    CALL put_submodel_att (status, USE_SURFACE,  modstr_surface,  modver_surface)  ; IF (status/=0) RETURN
    CALL put_submodel_att (status, USE_TBUDGET,  modstr_tbudget,  modver_tbudget)  ; IF (status/=0) RETURN
    CALL put_submodel_att (status, USE_TAGGING,  modstr_tagging,  modver_tagging)  ; IF (status/=0) RETURN
    CALL put_submodel_att (status, USE_TIMEPOS,  modstr_timepos,  modver_timepos)  ; IF (status/=0) RETURN
    CALL put_submodel_att (status, USE_TNUDGE,   modstr_tnudge,   modver_tnudge)   ; IF (status/=0) RETURN
    CALL put_submodel_att (status, USE_TPULSE,   modstr_tpulse,   modver_tpulse)   ; IF (status/=0) RETURN
    CALL put_submodel_att (status, USE_TREXP,    modstr_trexp,    modver_trexp)    ; IF (status/=0) RETURN
    CALL put_submodel_att (status, USE_TROPOP,   modstr_tropop,   modver_tropop)   ; IF (status/=0) RETURN
    CALL put_submodel_att (status, USE_TRSYNC,   modstr_trsync,   modver_trsync)   ; IF (status/=0) RETURN
    CALL put_submodel_att (status, USE_UBCNOX,   modstr_ubcnox,   modver_ubcnox)     ; IF (status/=0) RETURN
    CALL put_submodel_att (status, USE_VAHR,     modstr_vahr,     modver_vahr)     ; IF (status/=0) RETURN
    CALL put_submodel_att (status, USE_VAXTRA,   modstr_vaxtra,   modver_vaxtra)     ; IF (status/=0) RETURN
    CALL put_submodel_att (status, USE_VERTDIFF, modstr_vertdiff, modver_vertdiff) ; IF (status/=0) RETURN
    CALL put_submodel_att (status, USE_VERTEX,   modstr_vertex,   modver_vertex)   ; IF (status/=0) RETURN
    CALL put_submodel_att (status, USE_VISO,     modstr_viso,     modver_viso)     ; IF (status/=0) RETURN
    CALL put_submodel_att (status, USE_VISOP,    modstr_visop,    modver_visop) ; IF (status /= 0) RETURN
    CALL put_submodel_att (status, USE_MMD2WAY,  modstr_mmd2way,  modver_mmd2way)  ; IF (status/=0) RETURN

    CALL put_submodel_att (status, USE_CLAMS,       modstr_clams,       modver_clams); IF (status/=0) RETURN
    CALL put_submodel_att (status, USE_CLAMSTRAJ,   modstr_clamstraj,   modver_clamstraj); IF (status/=0) RETURN
    CALL put_submodel_att (status, USE_CLAMSCHEM,   modstr_clamschem,   modver_clamschem); IF (status/=0) RETURN
    CALL put_submodel_att (status, USE_CLAMSCHEME5, modstr_clamscheme5, modver_clamscheme5); IF (status/=0) RETURN
    CALL put_submodel_att (status, USE_CLAMSSEDI,   modstr_clamssedi,   modver_clamssedi); IF (status/=0) RETURN
    CALL put_submodel_att (status, USE_CLAMSMIX,    modstr_clamsmix,    modver_clamsmix); IF (status/=0) RETURN
    CALL put_submodel_att (status, USE_CLAMSBMIX,   modstr_clamsbmix,   modver_clamsbmix); IF (status/=0) RETURN
    CALL put_submodel_att (status, USE_CLAMSRDFRC,  modstr_clamsrdfrc,  modver_clamsrdfrc); IF (status/=0) RETURN
    CALL put_submodel_att (status, USE_CLAMSCIRRUS, modstr_clamscirrus, modver_clamscirrus); IF (status/=0) RETURN
    CALL put_submodel_att (status, USE_CLAMSTRACER, modstr_clamstracer, modver_clamstracer); IF (status/=0) RETURN
    CALL put_submodel_att (status, USE_CLAMSDEEPCONV, modstr_clamsdeepconv, modver_clamsdeepconv); IF (status/=0) RETURN
! ju_ch_20110429-
#endif

#ifdef ICON
!!$    CALL put_submodel_att (status, USE_TROPOP,      modstr_tropop,      modver_tropop); IF (status /= 0) RETURN
    CALL put_submodel_att (status, USE_PTRAC,    modstr_ptrac,    modver_ptrac)    ; IF (status/=0) RETURN
#endif

#ifdef MBM_MPIOM
    CALL put_submodel_att (status, USE_MPIOM,    modstr_mpiom,    modver_mpiom)    ; IF (status/=0) RETURN
#endif

#ifdef BLANK
#if defined (MBM_BLANK)
    CALL put_submodel_att (status, USE_PTRAC,    modstr_ptrac,    modver_ptrac)    ; IF (status/=0) RETURN
    CALL put_submodel_att (status, USE_ORBIT,    modstr_orbit,    modver_orbit)    ; IF (status/=0) RETURN
    CALL put_submodel_att (status, USE_SUBMOD1,  modstr_submod1,  modver_submod1)  ; IF (status/=0) RETURN
    CALL put_submodel_att (status, USE_SUBMOD2,  modstr_submod2,  modver_submod2)  ; IF (status/=0) RETURN
    CALL put_submodel_att (status, USE_SUBMOD3,  modstr_submod3,  modver_submod3)  ; IF (status/=0) RETURN
    CALL put_submodel_att (status, USE_TESTEVENT,  modstr_testevent,  modver_testevent)  ; IF (status/=0) RETURN
#endif
#if defined(MBM_QBO)
    CALL put_submodel_att (status, USE_QBO,      modstr_qbo,      modver_qbo)    ; IF (status/=0) RETURN
#endif
#if defined(MBM_DISSOC)
    CALL put_submodel_att (status, USE_DISSOC,   modstr_dissoc,   modver_dissoc)    ; IF (status/=0) RETURN
#endif
#if defined(MBM_RAD)
    CALL put_submodel_att (status, USE_ORBIT,    modstr_orbit,    modver_orbit)    ; IF (status/=0) RETURN
    CALL put_submodel_att (status, USE_CLOUDOPT, modstr_cloudopt, modver_cloudopt) ; IF (status/=0) RETURN
    CALL put_submodel_att (status, USE_RAD,      modstr_rad,      modver_rad)      ; IF (status/=0) RETURN
    CALL put_submodel_att (status, USE_AEROPT,   modstr_aeropt,   modver_aeropt)   ; IF (status/=0) RETURN
#endif
#endif

#ifdef MBM_CLAMS
    CALL put_submodel_att (status, USE_CLAMS,       modstr_clams,       modver_clams); IF (status/=0) RETURN
    CALL put_submodel_att (status, USE_CLAMSTRAJ,   modstr_clamstraj,   modver_clamstraj); IF (status/=0) RETURN
    CALL put_submodel_att (status, USE_CLAMSCHEM,   modstr_clamschem,   modver_clamschem); IF (status/=0) RETURN
!    CALL put_submodel_att (status, USE_CLAMSCHEME5, modstr_clamscheme5, modver_clamscheme5); IF (status/=0) RETURN
    CALL put_submodel_att (status, USE_CLAMSSEDI,   modstr_clamssedi,   modver_clamssedi); IF (status/=0) RETURN
    CALL put_submodel_att (status, USE_CLAMSMIX,    modstr_clamsmix,    modver_clamsmix); IF (status/=0) RETURN
    CALL put_submodel_att (status, USE_CLAMSBMIX,   modstr_clamsbmix,   modver_clamsbmix); IF (status/=0) RETURN
!    CALL put_submodel_att (status, USE_CLAMSRDFRC,  modstr_clamsrdfrc,  modver_clamsrdfrc); IF (status/=0) RETURN
    CALL put_submodel_att (status, USE_CLAMSCIRRUS, modstr_clamscirrus, modver_clamscirrus); IF (status/=0) RETURN
    CALL put_submodel_att (status, USE_CLAMSTRACER, modstr_clamstracer, modver_clamstracer); IF (status/=0) RETURN
    CALL put_submodel_att (status, USE_CLAMSDEEPCONV, modstr_clamsdeepconv, modver_clamsdeepconv); IF (status/=0) RETURN
! ju_ch_20110429-
#endif

#if defined(VERTICO)
    CALL put_submodel_att (status, USE_PTRAC,    modstr_ptrac,    modver_ptrac)    ; IF (status/=0) RETURN
    CALL put_submodel_att (status, USE_DDEP,     modstr_ddep,     modver_ddep)     ; IF (status/=0) RETURN
    CALL put_submodel_att (status, USE_MXL,      modstr_mxl,      modver_mxl)      ; IF (status/=0) RETURN
    CALL put_submodel_att (status, USE_MEGAN,    modstr_megan,    modver_megan)    ; IF (status/=0) RETURN
    CALL put_submodel_att (status, USE_MECCA,    modstr_mecca,    modver_mecca)    ; IF (status/=0) RETURN
    CALL put_submodel_att (status, USE_OFFEMIS,  modstr_offemis,  modver_offemis)  ; IF (status/=0) RETURN
    CALL put_submodel_att (status, USE_ONEMIS,   modstr_onemis,   modver_onemis)   ; IF (status/=0) RETURN
    CALL put_submodel_att (status, USE_JVAL,     modstr_jval,     modver_jval)     ; IF (status/=0) RETURN
#endif
    ! ### ADD HERE

    status = 0

  END SUBROUTINE switch_init
  ! -------------------------------------------------------------

  ! -------------------------------------------------------------
  SUBROUTINE put_submodel_att(status, LUSE, modstr, modver)

    ! MESSy
    USE messy_main_channel,       ONLY: new_attribute

    ! I/O
    INTEGER,      INTENT(OUT) :: status
    LOGICAL,      INTENT(in)  :: LUSE
    CHARACTER(*), INTENT(in)  :: modstr
    CHARACTER(*), INTENT(in)  :: modver

    status = 0

    IF (LUSE) THEN
       CALL new_attribute(status,  &
            'MESSy_'//TRIM(modstr), c= 'version '//TRIM(modver))
    END IF

  END SUBROUTINE put_submodel_att
  ! -------------------------------------------------------------

! ***************************************************************************
END MODULE messy_main_switch_bi
! ***************************************************************************

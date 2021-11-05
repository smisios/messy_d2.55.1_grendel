
MODULE messy_mpiom_mem_e5

!
#if defined MPIOM_13B
  !MESSy-MPIOM
  USE mo_param1            
  USE mo_parallel
  USE mo_commo1
  USE mo_levitus
  USE mo_commoau1
  USE mo_commoau2
  USE mo_commoau3
  USE mo_ocectlcom
  USE mo_octdiff
  USE mo_elicom
  USE mo_para2
  USE mo_commobbl
  USE mo_tro
  USE mo_adpo
  USE mo_units   
  USE mo_octher
  USE mo_ocice  
!tidal mode
  USE MO_TIDAL
!#ifdef CORE
! USE MO_NCAR_OCEAN_FLUXES
!#else
  USE MO_OMIP
!#endif
#elif defined MPIOM_2000  
  !MESSy-MPIOM
  USE mo_contro,      ONLY: contro
  USE mo_param1,      ONLY: ie, je, ie_g, je_g, ke, icycli, stabn, jto, kbb   &
                          , ielimi, kep, set_param1
  USE mo_parallel,    ONLY: p_deco, gather_arr => gather                      &
                          , scatter_arr => scatter, nprocxy, nprocx, nprocy   &
                          , nprocy, have_g_is, have_g_ie, have_g_js           &
                          , have_g_je, global_min
!mz_bk_20110222           , global_array_sum,                                 &
!mz_bk_20110222           , stop_all
  USE mo_commo1,      ONLY: lgmdiag, lforcediag, lhfldiag, lconvdiag          &
                          , ldiffdiag, lcalcdifi, lisopyc, iocad, ladpo       &
                          , ladfs, imean, ibolk, ibbl_transport, iocaduv      &
                          , aus, nfixYearLen, CAH00, CSTABEPS, CRELSAL        &
                          , CRELTEM, IAUFR, IAUFW, ISTART, I3DREST            &
                          , IOASISFLUX, CDVOCON, CAVOCON, rleadclose, IOCONV  &
                          , ICONTRO, lmpitype, lnonblock, ltidal              &
                          , ltidal_diag, lswr_jerlov, SAF, TAF, DZW, dz, DTI  &
                          , NDTDAY, AULAPUV, TIESTW, conn, cono, fclou, fprec &
                          , KM, KBM, gila, giph, dlxp, dlyp, zo, ddpo         &
                          , dt, lyears, lmonts, ldays, lyear1, lmont1         &
                          , istart_new_topo_update, istart_new_run            &
                          , istart_new_topo, THO, SICTHO, SICOMO, SAO, SICSNO &
                          , uoo, voe, dvo, avo, wo, zo, z1o, sicuo, sicve     &
                          , hibzeto, hibeto, hibzete, hibete, tice, fu10      &
                          , txo, tye, lwith_barotropic_stokes_drift           &
                          , stabio, uaccel, vaccel, lbounds_exch_tp, po       &
                          , eminpo, rivrun, ltstranspose, ltswrite            &
                          , dlxp_g, dlyp_g, dlxu_g, dlyu_g, dlxv_g, dlyv_g    &
                          , amsuo, amsue, dpio, tiestu                        &
                          , almzer, weto, rhoo, kbot                          &
                          , alloc_mem_forcing, alloc_mem_stokes_drift         &
                          , alloc_mem_commo1, alloc_mem_gmbolus
  USE mo_levitus,     ONLY: spongezone, spzndamp_time                         &
                          , init_levitus_2d, init_levitus_3d, relax_ts        &
                          , levitus_read_3d_restore                           &
                          , levitus_read_3d_stratification                    &
                          , levitus_horizontal_stratification                 &
                          , levitus_per_month_setup
  USE mo_commoau1,    ONLY: h0, hmin, armin, armax, hsntoice, sicthmin, sice  &
                          , isnflg, tmelt
  USE mo_commoau2,    ONLY: alloc_mem_commoau2
  USE mo_commoau3,    ONLY: alloc_mem_commoau3
  USE mo_constants,   ONLY: api  
!mz_bk_20110222  USE mo_nudge_ts
  USE mo_ocvisc,      ONLY: bofric, rayfric, ocvisc
  USE mo_ocean_vertical_mixing, ONLY: av0, dv0                                &
                          , setup_ocean_vertical_mixing, calc_rinum
  USE mo_grid,        ONLY: boden, coriol, setup_grid, lwith_one_layer_shelfs &
!mz_bk_20110222           , wrte_gridinfo                                     &
                          , cell_thickness, gila_g, giph_g, weto_g 
  USE mo_runoff,      ONLY: river_runoff_omip_ini, river_runoff_stations_ini  &
                          , glac_calv_ini, luse_river_runoff_stations, numriv &
                          , luse_glac_calv, numglac
  USE mo_swr_absorption, ONLY : alloc_mem_swr_absorption                      &
!mz_bk_20110222           , jerlov_swr_absorption, old_swr_absorption         &
                          , jerlov_atten, jerlov_bluefrac,lfb_bgc_oce         &
                          , swrab, swsum
  USE mo_diffusion,   ONLY: alloc_mem_octdiff, octdiff_base, octdiff_trf      &
                          , aulapts, ah00
  USE mo_planetary_constants, ONLY:rhoref_water, rhoref_snow, rhoicwa, rhosnwa
!mz_bk_20110222           , rhoref_ice
  USE mo_elicom,      ONLY: alloc_mem_elicom
  USE mo_para2,       ONLY: iter_sor, iter_sor_hack, rtsorpar, rtsorpar_hack  &
                          , alloc_mem_para2
  USE mo_commobbl,    ONLY: alloc_mem_commobbl, findalfa, slopetrans
  USE mo_tro,         ONLY: itprep, trotest2, troneu, ocvtro, update_zo       &
                          , alloc_mem_dilcor, correct_zo
  USE mo_adpo,        ONLY: ocadpo_trf, ocadpo_base, alloc_mem_adpo
  USE mo_units,       ONLY: io_stdout, setunits
  USE mo_octher,      ONLY: octher, calc_dens 
  USE MO_OCICE,       ONLY: lsaoclose, ocice
  USE mo_tidal,       ONLY: alloc_mem_tidal,foreph_ini,foreph,tipouv, mmccdt
  USE mo_boundsexch,  ONLY: bounds_exch, alloc_mem_bounds_exch_halo
  USE mo_diagnosis,   ONLY: calc_icecutoff                                    &
       ! mz_bk_20110315+
                          , calc_global_mean, global_sum_salinity             &
                          , global_sum_temperature, global_mass               &
                          , global_volume, global_salt_content
       ! mz_bk_20110315-
  USE mo_convection,  ONLY: convection 
!  USE mo_mpi,         ONLY: p_io, p_pe, p_bcast, p_abort, init_mpi_datatypes
  

!#ifdef CORE
! USE MO_NCAR_OCEAN_FLUXES
!#else
!mz_bk_20110222  USE MO_OMIP
!#endif
#endif

END MODULE messy_mpiom_mem_e5

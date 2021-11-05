! **************************************************************************
MODULE messy_main_bmluse_bi
! **************************************************************************

! BMLUSE should be empty for models which do not change much anymore
#if defined(ECHAM5) || defined(COSMO) || defined(CESM1)

  IMPLICIT NONE

#endif 

! BMLUSE should be the only(!) use point of ICON variables for access in MESSy
#ifdef ICON
    USE mo_cdi_constants,               ONLY: GRID_REGULAR_LONLAT
    USE mo_communication,               ONLY: exchange_data, idx_1d, blk_no, idx_no         &
       &                                    , t_comm_gather_pattern                         &
       &                                    , get_np_send, get_np_recv, get_pelist_recv     &
       &                                    , t_ScatterPattern                              &
       &                                    , t_comm_pattern, delete_comm_pattern
    USE mo_decomposition_tools,         ONLY: get_local_index
    USE mo_dynamics_config,             ONLY: nnow,nnew, nnow_rcf, nnew_rcf
    USE mo_dist_dir,                    ONLY: dist_dir_get_owners
    USE mo_exception,                   ONLY: finish
    USE mo_grid_config,                 ONLY: n_phys_dom, grid_sphere_radius
    USE mo_grid_subset,                 ONLY: get_index_range
    USE mo_impl_constants,              ONLY: max_dom, min_rlcell_int, min_rlcell
    USE mo_impl_constants_grf,          ONLY: grf_bdywidth_c, grf_bdywidth_e, grf_bdyintp_start_c &
       &                                    , grf_bdyintp_end_c
    USE mo_io_config,                   ONLY: default_read_method                           &
       &                                    , read_netcdf_broadcast_method                  &
       &                                    , read_netcdf_distribute_method
    USE mo_kind,                        ONLY: wp, sp
    USE mo_loopindices,                 ONLY: get_indices_c
    USE mo_lonlat_grid,                 ONLY: compute_lonlat_blocking, compute_lonlat_specs &
       &                                    , t_lon_lat_grid
    USE mo_math_types,                  ONLY: t_geographical_coordinates
    USE mo_model_domain,                ONLY: t_patch, p_patch, t_phys_patch, p_phys_patch  &
       &                                    , t_grid_cells, p_patch_local_parent
    USE mo_mpi,                         ONLY: my_process_is_io, my_process_is_mpi_workroot  &
       &                                    , my_process_is_mpi_seq
    USE mo_name_list_output_gridinfo,   ONLY: collect_all_grid_info                         &
       &                                    , GRID_INFO_NONE, GRID_INFO_FILE                &
       &                                    , GRID_INFO_BCAST
    USE mo_name_list_output_init,       ONLY: patch_info
    USE mo_name_list_output_types,      ONLY: l_output_phys_patch                           &
       &                                    , t_reorder_info                                &
       &                                    , icell, ivert, iedge
    USE mo_nonhydro_state,              ONLY: p_nh_state, p_nh_state_lists
    USE mo_ext_data_state,              ONLY: ext_data ! ub_ak_20190729
    USE mo_parallel_config,             ONLY: nproma
    USE mo_read_interface,              ONLY: openInputFile, closeFile                      &
       &                                    , t_stream_id                                   &
       &                                    , on_cells, on_edges, on_vertices               &
       &                                    , read_3D_extdim, read_2d_extdim &
       &                                    , read_3D, read_2d                              &
       &                                    , read_1D, read_bcast_REAL_2D                   &
       &                                    , read_0D_real &
                                            , read_inq_varexists
    USE mo_run_config,                  ONLY: output_mode, number_of_grid_used, iqv, iqm_max, ntracer
    USE mo_sync,                        ONLY: sync_patch_array, sync_patch_array_mult, SYNC_C
    USE mo_time_config,                 ONLY: time_config
    USE mo_util_string,                 ONLY: int2string
    USE mo_vertical_coord_table,        ONLY: vct
    USE mtime,                          ONLY: MAX_DATETIME_STR_LEN, datetime, timedelta      &
       &                                    , newDatetime, newTimedelta, deallocateTimedelta &
       &                                    , getTotalMilliSecondsTimeDelta                  &
       &                                    , OPERATOR(>=), OPERATOR(<=), OPERATOR(-)

    USE mo_lnd_nwp_config,              ONLY: nlev_snow, nlev_soil
    USE mo_echam_phy_memory,            ONLY: prm_field

    ! FOR CHANNEL
    USE mo_vertical_coord_table,        ONLY: vct_a, vct_b
    USE mo_util_uuid,                   ONLY: uuid_unparse
    
    ! FOR TRACER
    USE mo_nwp_phy_state,               ONLY: prm_nwp_tend
    USE mo_linked_list,                 ONLY: t_list_element, t_var_list
    USE mo_var_metadata_types,          ONLY: t_var_metadata, t_var_metadata_dynamic
    USE mo_var_list,                    ONLY: get_var_name
    USE mo_grid_config,                 ONLY: n_dom
    USE mo_atm_phy_nwp_config,          ONLY: atm_phy_nwp_config
    USE mo_advection_config,            ONLY: advection_config
    
#ifdef ICON_2_1_01
    USE mo_cdi_constants,               ONLY: ZA_HYBRID, ZA_SURFACE
    USE mo_communication,               ONLY: setup_comm_pattern
#else
    USE mo_communication_factory,       ONLY: setup_comm_pattern
    USE mo_zaxis_type,                  ONLY: ZA_HYBRID, ZA_SURFACE
#endif
    USE src_turbdiff,                   ONLY: modvar

    IMPLICIT NONE
    PUBLIC
#endif

! **************************************************************************
END MODULE messy_main_bmluse_bi
! **************************************************************************

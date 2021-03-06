MODULE MESSY_AEROPT_SETS

  USE MESSY_AEROPT_MEM

  IMPLICIT NONE


  CONTAINS

!===============================================================================

    SUBROUTINE MAP_LUT_STRUCT(idx)

      ! This subroutine maps the pointers of the struct for the lookuptable for 
      ! the aerosol optical properties with the index number "idx" to the 
      ! pointers used in the core for the calculations

      INTEGER, INTENT(IN) :: idx

      num_sw_intervals => tab_set(idx)%num_sw_intervals
      num_lw_intervals => tab_set(idx)%num_lw_intervals
      num_aerosols     => tab_set(idx)%num_aerosols
      nmodrad          => tab_set(idx)%nmodrad

      idim_nr          => tab_set(idx)%idim_nr         
      idim_ni          => tab_set(idx)%idim_ni         
      idim_msp         => tab_set(idx)%idim_msp        
                                                         
      sigma_g          => tab_set(idx)%sigma_g

      sw_intervals     => tab_set(idx)%sw_intervals     
      lw_intervals     => tab_set(idx)%lw_intervals     
      nr_min           => tab_set(idx)%nr_min           
      nr_max           => tab_set(idx)%nr_max           
      ni_min           => tab_set(idx)%ni_min           
      ni_max           => tab_set(idx)%ni_max           
      msp_min          => tab_set(idx)%msp_min          
      msp_max          => tab_set(idx)%msp_max          
      nr_step          => tab_set(idx)%nr_step          
      ni_step          => tab_set(idx)%ni_step          
      msp_step         => tab_set(idx)%msp_step         
      log_ni_min       => tab_set(idx)%log_ni_min       
      log_msp_min      => tab_set(idx)%log_msp_min     

      lambda_sw        => tab_set(idx)%lambda_sw        
      lambda_lw        => tab_set(idx)%lambda_lw        
      weight_sw        => tab_set(idx)%weight_sw        
      weight_lw        => tab_set(idx)%weight_lw             
      sumweight_sw     => tab_set(idx)%sumweight_sw    
      sumweight_lw     => tab_set(idx)%sumweight_lw    
      ref_re_sw        => tab_set(idx)%ref_re_sw        
      ref_re_lw        => tab_set(idx)%ref_re_lw        
      ref_im_sw        => tab_set(idx)%ref_im_sw        
      ref_im_lw        => tab_set(idx)%ref_im_lw        
      lut_sw_sigma     => tab_set(idx)%lut_sw_sigma    
      lut_lw_sigma     => tab_set(idx)%lut_lw_sigma    
      lut_sw_omega     => tab_set(idx)%lut_sw_omega    
      lut_sw_gamma     => tab_set(idx)%lut_sw_gamma    
      rad_sw_filename  => tab_set(idx)%rad_sw_filename
      rad_lw_filename  => tab_set(idx)%rad_lw_filename

    END SUBROUTINE MAP_LUT_STRUCT

!===============================================================================

    SUBROUTINE ALLOCATE_LUT_STRUCT(IDX)

      INTEGER, INTENT(IN) :: IDX

      INTEGER :: nmod

      nmod = tab_set(idx)%nmodrad

      ALLOCATE(tab_set(idx)%idim_nr(2))
      ALLOCATE(tab_set(idx)%idim_ni(2))
      ALLOCATE(tab_set(idx)%idim_msp(2))

      ALLOCATE(tab_set(idx)%sw_intervals(nsw))
      ALLOCATE(tab_set(idx)%lw_intervals(jpband))
      ALLOCATE(tab_set(idx)%sumweight_sw(nsw))
      ALLOCATE(tab_set(idx)%sumweight_lw(jpband))


      ALLOCATE(tab_set(idx)%sigma_g(2,nmod))

      ALLOCATE(tab_set(idx)%nr_min(2,nmod))
      ALLOCATE(tab_set(idx)%nr_max(2,nmod))
      ALLOCATE(tab_set(idx)%ni_min(2,nmod))
      ALLOCATE(tab_set(idx)%ni_max(2,nmod))
      ALLOCATE(tab_set(idx)%msp_min(2,nmod))
      ALLOCATE(tab_set(idx)%msp_max(2,nmod))
      ALLOCATE(tab_set(idx)%nr_step(2,nmod))
      ALLOCATE(tab_set(idx)%ni_step(2,nmod))
      ALLOCATE(tab_set(idx)%msp_step(2,nmod))       
      ALLOCATE(tab_set(idx)%log_ni_min(2,nmod))
      ALLOCATE(tab_set(idx)%log_msp_min(2,nmod))

      ALLOCATE(tab_set(idx)%lambda_sw(max_sw_intervals))
      ALLOCATE(tab_set(idx)%lambda_lw(max_lw_intervals))
      ALLOCATE(tab_set(idx)%weight_sw(max_sw_intervals))
      ALLOCATE(tab_set(idx)%weight_lw(max_lw_intervals))

    END SUBROUTINE ALLOCATE_LUT_STRUCT

!===============================================================================

    SUBROUTINE MAP_OPTSET_STRUCT(optset)


      ! This subroutine maps the pointers for the STRUCT information
      ! for all aspects of the optset call to the pointers used for
      ! the calculations

      TYPE(t_aero_set), INTENT(IN), POINTER :: optset

      aero_set_name         => optset%aero_set_name
      aot_sw                => optset%aot_sw
      aot_lw                => optset%aot_lw
      omega_sw              => optset%omega_sw
      gamma_sw              => optset%gamma_sw
      l_tracer              => optset%l_tracer
      exclude_string_full   => optset%exclude_string_full
      exclude_string_spec   => optset%exclude_string_spec
      aermodelname          => optset%aermodelname
      lcalc_jval            => optset%lcalc_jval
      n_jv_bands            => optset%n_jv_bands
      IF (lcalc_jval) THEN
        jv_asca             => optset%jv_asca
        jv_aabs             => optset%jv_aabs
        jv_ga               => optset%jv_ga
      END IF
      IF (.NOT. optset%TANRE) THEN
        tracerset           => optset%tracerset
        num_opt_wavelens    => optset%num_opt_wavelens
        aot_opt             => optset%aot_opt
        extcoeff_opt        => optset%extcoeff_opt
        wetradius           => optset%wetradius
        dryradius           => optset%dryradius
        aernumber           => optset%aernumber
        lcalc_seasalt       => optset%lcalc_seasalt
        l_extmixt           => optset%l_extmixt
        aerspec             => optset%aerspec_e5
        sigma               => optset%sigma
        naertot             => optset%naertot
        nmod                => optset%nmod
        nsol                => optset%nsol
        lut_number          => optset%lut_number
        ns                  => optset%ns
        ks                  => optset%ks
        as                  => optset%as
        cs                  => optset%cs
        radmod              => optset%radmod
        rad_diag_wavelen    => optset%diag_wavelen
        
        lambda              => optset%lambda
        lambda_squared      => optset%lambda_squared
        weight              => optset%weight
        ref_re              => optset%ref_re
        ref_im              => optset%ref_im
        int2band            => optset%int2band
        znwave              => optset%znwave

      END IF

    END SUBROUTINE MAP_OPTSET_STRUCT
    

END MODULE MESSY_AEROPT_SETS

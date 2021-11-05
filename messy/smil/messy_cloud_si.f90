
MODULE messy_cloud_si

#if defined (ECHAM5) || defined (CESM1)

 USE messy_main_grid_def_mem_bi, ONLY: nproma, ngpblks
 USE messy_main_data_bi,         ONLY: aclc, acdnc
 USE messy_main_mpi_bi,          ONLY: p_pe,p_bcast,p_io
 USE messy_main_rnd_bi,          ONLY: RND_MP_PIN, rnd_init_bi, rnd_finish_bi  &
                                     , rnd_number_bi
 USE messy_main_tools,         ONLY: PTR_3D_ARRAY
 USE messy_main_constants_mem, ONLY: dp, STRLEN_MEDIUM
 USE messy_cloud
#ifdef MESSYTENDENCY
 !tendency budget
 USE messy_main_tendency_bi,   ONLY: mtend_get_handle,       &
                                     mtend_get_start_l,      &
                                     mtend_add_l,            &
                                     mtend_register,         &    
                                     mtend_id_t,             &
                                     mtend_id_q,             &
                                     mtend_id_xl,            &
                                     mtend_id_xi,            &
                                     mtend_id_tracer
#endif

 !WISO++
 USE messy_main_data_wiso_bi,  ONLY: l_wiso
 !WISO--

 IMPLICIT NONE

 SAVE
 INTRINSIC :: NULL

 INTEGER,  ALLOCATABLE, DIMENSION(:)   :: invb
 REAL(dp), ALLOCATABLE, DIMENSION(:,:) :: zbetaa
 REAL(dp), ALLOCATABLE, DIMENSION(:,:) :: zbetab
 REAL(dp), ALLOCATABLE, DIMENSION(:,:) :: zbetass

! pointers for channel objects of this module

 REAL(dp), POINTER, DIMENSION(:,:,:) :: frain    => NULL()      ! rainflux
 REAL(dp), POINTER, DIMENSION(:,:,:) :: fsnow    => NULL()      ! snowflux
 ! rainflux without production of 
 ! new rain from cloud processes
 REAL(dp), POINTER, DIMENSION(:,:,:) :: frain_no => NULL()     
 ! snowflux without production of 
 ! new snow from cloud processes
 REAL(dp), POINTER, DIMENSION(:,:,:) :: fsnow_no => NULL()    
 ! liquid water content kg/kg, 
 ! only inside cloud, not averaged over box!!!!
 REAL(dp), POINTER, DIMENSION(:,:,:) :: mlwc     => NULL()    
 ! snow/ice content kg/kg, 
 ! only inside cloud, not averaged over box!!!!
 REAL(dp), POINTER, DIMENSION(:,:,:) :: msic     => NULL()    
 ! precipitation formation rate kg/kg, averaged over box!!!!
 REAL(dp), POINTER, DIMENSION(:,:,:) :: mratep   => NULL()    
 ! snow/ice formation rate kg/kg, averaged over box!!!!
 REAL(dp), POINTER, DIMENSION(:,:,:) :: mratesi  => NULL()    
 ! precipitation evaporation rate kg/kg, averaged over box!!!!
 REAL(dp), POINTER, DIMENSION(:,:,:) :: mrevap   => NULL()    
 ! snow/ice sublimation rate kg/kg, averaged over box!!!!
 REAL(dp), POINTER, DIMENSION(:,:,:) :: mssubl   => NULL()   
 ! precipitating cloud cover
 REAL(dp), POINTER, DIMENSION(:,:,:) :: preccover => NULL()
 ! melting of frozen precip into the liquid phase
 REAL(dp), POINTER, DIMENSION(:,:,:) :: mimelt    => NULL()  
 ! ice sedimentation 
 REAL(dp), POINTER, DIMENSION(:,:,:) :: misedi    => NULL()   
 ! mz_ak_20051221+
 ! condensation of liquid water in cloud covered area (absolute value, no rate)
 REAL(dp), POINTER, DIMENSION(:,:,:) :: condens   => NULL()     
 ! mz_ak_20051221-

 !mz_sb_20170703+
 REAL(dp), POINTER, DIMENSION(:,:,:) :: w_sub            => NULL()
 REAL(dp), POINTER, DIMENSION(:,:,:) :: w_ave            => NULL()
 REAL(dp), POINTER, DIMENSION(:,:,:) :: w_grid           => NULL()
 REAL(dp), POINTER, DIMENSION(:,:,:) :: sigwBN           => NULL()
 REAL(dp), POINTER, DIMENSION(:,:,:) :: ndropBN          => NULL()
 REAL(dp), POINTER, DIMENSION(:,:,:) :: dsulfBN          => NULL()
 REAL(dp), POINTER, DIMENSION(:,:,:) :: ndust_aiBN       => NULL()
 REAL(dp), POINTER, DIMENSION(:,:,:) :: ddust_aiBN       => NULL()
 REAL(dp), POINTER, DIMENSION(:,:,:) :: ndust_ciBN       => NULL()
 REAL(dp), POINTER, DIMENSION(:,:,:) :: ddust_ciBN       => NULL()
 REAL(dp), POINTER, DIMENSION(:,:,:) :: norgBN           => NULL()
 REAL(dp), POINTER, DIMENSION(:,:,:) :: dorgBN           => NULL()
 REAL(dp), POINTER, DIMENSION(:,:,:) :: nsootBN          => NULL()
 REAL(dp), POINTER, DIMENSION(:,:,:) :: dsootBN          => NULL()
 REAL(dp), POINTER, DIMENSION(:,:,:) :: nsolo_ksBN       => NULL()
 REAL(dp), POINTER, DIMENSION(:,:,:) :: dsolo_ksBN       => NULL()
 REAL(dp), POINTER, DIMENSION(:,:,:) :: nsolo_asBN       => NULL()
 REAL(dp), POINTER, DIMENSION(:,:,:) :: dsolo_asBN       => NULL()
 REAL(dp), POINTER, DIMENSION(:,:,:) :: nsolo_csBN       => NULL()
 REAL(dp), POINTER, DIMENSION(:,:,:) :: dsolo_csBN       => NULL()
 REAL(dp), POINTER, DIMENSION(:,:,:) :: smaxice_cirrusBN => NULL()
 REAL(dp), POINTER, DIMENSION(:,:,:) :: sc_ice_cirrusBN  => NULL()
 REAL(dp), POINTER, DIMENSION(:,:,:) :: nlim_cirrusBN    => NULL()
 REAL(dp), POINTER, DIMENSION(:,:,:) :: nhet_cirrusBN    => NULL()
 REAL(dp), POINTER, DIMENSION(:,:,:) :: nice_cirrusBN    => NULL()
 REAL(dp), POINTER, DIMENSION(:,:,:) :: dice_cirrusBN    => NULL()
 REAL(dp), POINTER, DIMENSION(:,:,:) :: sigwpre_cirrusBN => NULL()
 REAL(dp), POINTER, DIMENSION(:,:,:) :: smaxice_immBN    => NULL()
 REAL(dp), POINTER, DIMENSION(:,:,:) :: sc_ice_immBN     => NULL()
 REAL(dp), POINTER, DIMENSION(:,:,:) :: nice_immBN       => NULL()
 REAL(dp), POINTER, DIMENSION(:,:,:) :: dice_immBN       => NULL()
 REAL(dp), POINTER, DIMENSION(:,:,:) :: sigwpre_immBN    => NULL()
 REAL(dp), POINTER, DIMENSION(:,:,:) :: newIC_cirri      => NULL()
 REAL(dp), POINTER, DIMENSION(:,:,:) :: newICR_cirri     => NULL()
 REAL(dp), POINTER, DIMENSION(:,:,:) :: newIC_cnt_therm  => NULL()
 REAL(dp), POINTER, DIMENSION(:,:,:) :: newIC_imm        => NULL()
 REAL(dp), POINTER, DIMENSION(:,:,:) :: newIC_mix        => NULL()
 REAL(dp), POINTER, DIMENSION(:,:,:) :: newIC_mix_fin    => NULL()
 !mz_sb_20170703-

! cloud droplet number concentration (E5 parameterisation, changing with time)
 REAL(dp), POINTER, DIMENSION(:,:,:) :: acdnc_1      => NULL()  
! cloud droplet number concentration after Rotstayn
 REAL(dp), POINTER, DIMENSION(:,:,:) :: acdnc_2      => NULL()  
! cloud droplet number concentration after Jones
 REAL(dp), POINTER, DIMENSION(:,:,:) :: acdnc_3      => NULL()   
! cloud droplet number concentration after Menon
 REAL(dp), POINTER, DIMENSION(:,:,:) :: acdnc_4      => NULL()
! cloud droplet number concentration (activation) after Abdul-Razzak - Ghan
 REAL(dp), POINTER, DIMENSION(:,:,:)     :: acdnc_5     => NULL()   
 REAL(dp), POINTER, DIMENSION(:,:,:,:,:) :: p_arg_frac  => NULL() 
 REAL(dp), POINTER, DIMENSION(:,:,:,:,:) :: p_arg_frac2 => NULL() 
 REAL(dp), POINTER, DIMENSION(:,:,:,:,:) :: p_arg_mfrac => NULL() 
 REAL(dp), POINTER, DIMENSION(:,:,:,:,:) :: p_arg_num   => NULL() 
 REAL(dp), POINTER, DIMENSION(:,:,:,:)   :: arg_frac    => NULL() 
 REAL(dp), POINTER, DIMENSION(:,:,:,:)   :: arg_frac2   => NULL() 
 REAL(dp), POINTER, DIMENSION(:,:,:,:)   :: arg_mfrac   => NULL() 
 REAL(dp), POINTER, DIMENSION(:,:,:,:)   :: arg_num     => NULL() 
 REAL(dp), POINTER, DIMENSION(:,:,:,:,:) :: p_arg_rad   => NULL() 
 REAL(dp), POINTER, DIMENSION(:,:,:,:)   :: arg_rad     => NULL() 
 REAL(dp), POINTER, DIMENSION(:,:,:,:,:) :: p_arg_scrit => NULL() 
 REAL(dp), POINTER, DIMENSION(:,:,:,:)   :: arg_scrit => NULL() 
! cloud droplet number concentration (activation) after Lin & Leaitch
 REAL(dp), POINTER, DIMENSION(:,:,:)     :: acdnc_6      => NULL()   
 REAL(dp), POINTER, DIMENSION(:,:,:)     :: np_pot_activ => NULL()
 REAL(dp), POINTER, DIMENSION(:,:,:,:,:) :: p_lin_frac   => NULL() 
 REAL(dp), POINTER, DIMENSION(:,:,:,:,:) :: p_lin_mfrac  => NULL() 
 REAL(dp), POINTER, DIMENSION(:,:,:,:,:) :: p_lin_num    => NULL() 
 REAL(dp), POINTER, DIMENSION(:,:,:,:)   :: lin_frac     => NULL() 
 REAL(dp), POINTER, DIMENSION(:,:,:,:)   :: lin_mfrac    => NULL() 
 REAL(dp), POINTER, DIMENSION(:,:,:,:)   :: lin_num      => NULL() 
! cloud droplet number concentration (activation) after UAF (Karydis et al., 2016)
 REAL(dp), POINTER, DIMENSION(:,:,:) :: acdnc_7      => NULL()
 REAL(dp), POINTER, DIMENSION(:,:,:,:,:):: p_uaf_num     => NULL()
 REAL(dp), POINTER, DIMENSION(:,:,:,:,:):: p_uaf_frac    => NULL()
 REAL(dp), POINTER, DIMENSION(:,:,:,:,:):: p_uaf_mfrac   => NULL()
 REAL(dp), POINTER, DIMENSION(:,:,:,:,:):: p_uaf_AKK     => NULL()
 REAL(dp), POINTER, DIMENSION(:,:,:,:,:):: p_uaf_ei      => NULL()
 REAL(dp), POINTER, DIMENSION(:,:,:,:,:):: p_uaf_TPI     => NULL()
 REAL(dp), POINTER, DIMENSION(:,:,:,:,:):: p_uaf_DPG     => NULL()
 REAL(dp), POINTER, DIMENSION(:,:,:,:,:):: p_uaf_SIG     => NULL()
 REAL(dp), POINTER, DIMENSION(:,:,:,:):: uaf_num     => NULL()
 REAL(dp), POINTER, DIMENSION(:,:,:,:):: uaf_frac    => NULL()
 REAL(dp), POINTER, DIMENSION(:,:,:,:):: uaf_mfrac   => NULL()
 REAL(dp), POINTER, DIMENSION(:,:,:,:):: uaf_AKK     => NULL()
 REAL(dp), POINTER, DIMENSION(:,:,:,:):: uaf_TPI     => NULL()
 REAL(dp), POINTER, DIMENSION(:,:,:,:):: uaf_DPG     => NULL()
 REAL(dp), POINTER, DIMENSION(:,:,:,:):: uaf_SIG     => NULL()
 REAL(dp), POINTER, DIMENSION(:,:,:,:):: uaf_ei      => NULL()
! aerosol sulphate mass
 REAL(dp), POINTER, DIMENSION(:,:,:) :: aersulf_0    => NULL()   
! aerosol number pointer (in case of gm7 or eqsam) 
! with no aerosol number tracers
 REAL(dp), POINTER, DIMENSION(:,:,:,:) :: anumber    => NULL()
! switches
 LOGICAL     :: l_cdnc_calc     ! calculation of CDNCs
 INTEGER     :: i_cdnc_calc = 0 ! number of parameterisation for CDNCs
 INTEGER     :: i_cdnc_cpl  = 0 ! feedback of CDNC parameterisation
 INTEGER     :: sup_sat_scheme = 0 ! scheme for critical supersaturation

 CHARACTER(LEN=STRLEN_MEDIUM) :: aer_stream

! pointers for the coupling ro other submodels, routines
 CHARACTER(LEN=STRLEN_MEDIUM), ALLOCATABLE :: base_name(:)    !mz_vk_20170703
 LOGICAL :: USE_PSC=.FALSE.
 REAL(dp), POINTER, DIMENSION(:,:,:)   :: flt_pscreg => NULL()
 REAL(dp), POINTER, DIMENSION(:,:,:)   :: pqtec      => NULL()
 REAL(dp), POINTER, DIMENSION(:,:,:)   :: pxtecl     => NULL()
 REAL(dp), POINTER, DIMENSION(:,:,:)   :: pxteci     => NULL()
 REAL(dp), POINTER, DIMENSION(:,:,:)   :: pxtecnl    => NULL()
 REAL(dp), POINTER, DIMENSION(:,:,:)   :: pxtecni    => NULL()
 REAL(dp), POINTER, DIMENSION(:,:,:)   :: pvmixtau   => NULL()
 REAL(dp), POINTER, DIMENSION(:,:,:)   :: pvdiffp    => NULL()
#ifndef MESSYTENDENCY
 REAL(dp), POINTER, DIMENSION(:,:,:)   :: sdtdt_cloud  => NULL()
#endif

 REAL(dp), POINTER, DIMENSION(:,:)     :: conv_type  => NULL()
 REAL(dp), POINTER, DIMENSION(:,:)     :: conv_bot   => NULL()
 REAL(dp), POINTER, DIMENSION(:,:)     :: PCAPE      => NULL()

 REAL(dp), POINTER, DIMENSION(:,:,:)   :: pdiga    => NULL()

 REAL(dp)              :: fac_sulf, fac_ss, fac_om
 INTEGER               :: sulf_count, ss_count, ss2_count, om_count
 INTEGER, ALLOCATABLE  :: sulf_idt(:), om_idt(:)
 INTEGER, ALLOCATABLE  :: ss_idt(:), ss2_idt(:,:)

 REAL(dp), POINTER, DIMENSION(:,:,:)   :: w_gwd_kpr  => NULL()
 REAL(dp), POINTER, DIMENSION(:,:,:)   :: ampl_gwd   => NULL()
 REAL(dp), POINTER, DIMENSION(:,:)     :: mask_orogw => NULL()
 REAL(dp), POINTER, DIMENSION(:,:)     :: l_z        => NULL()

! related to ARG and UAF activation
 INTEGER               :: ksol, kmod, nspec, ktot
 REAL(dp), POINTER, DIMENSION(:)       :: aerosol_num => NULL()
 REAL(dp), POINTER, DIMENSION(:)       :: nu          => NULL()
 REAL(dp), POINTER, DIMENSION(:)       :: PHI         => NULL()
 REAL(dp), POINTER, DIMENSION(:)       :: EPS         => NULL()
 REAL(dp), POINTER, DIMENSION(:)       :: KAPPA       => NULL()
 REAL(dp), POINTER, DIMENSION(:)       :: molar       => NULL()
 REAL(dp), POINTER, DIMENSION(:)       :: aerdensity  => NULL()
 REAL(dp), POINTER, DIMENSION(:,:)     :: mass        => NULL()
 INTEGER,  POINTER, DIMENSION(:)       :: num_idx     => NULL()
 INTEGER,  POINTER, DIMENSION(:,:)     :: aer_idx     => NULL()

 REAL(dp), POINTER, DIMENSION(:)       :: sigma       => NULL()
 REAL(dp), POINTER, DIMENSION(:)       :: crdiv       => NULL()
 REAL(dp), POINTER, DIMENSION(:,:,:,:) :: wetradius   => NULL()
 REAL(dp), POINTER, DIMENSION(:,:,:,:) :: dryradius   => NULL()

 REAL(dp), POINTER, DIMENSION(:)       :: cmr2mmr     => NULL()
 LOGICAL,  POINTER, DIMENSION(:)       :: FHH
 TYPE (PTR_3D_ARRAY),POINTER, DIMENSION(:):: kappa_uaf 

 ! critical humidity for natural clouds, necessary for contrail calculation
 REAL(dp), POINTER, DIMENSION(:,:,:) :: rhc   => NULL()

! Name and indexes of the soluble modes
    CHARACTER(LEN=STRLEN_MEDIUM), DIMENSION(:), ALLOCATABLE :: sol_name
    INTEGER, DIMENSION(:), ALLOCATABLE :: sol_idx

! no tracer for aerosol number concentrations found
 LOGICAL :: no_trac = .false. 
!!$!----------------------------------------------------------------
!!$ LOGICAL :: l_tf           ! switch for correct use of time filter
!!$!----------------------------------------------------------------
! for cloud acc 
 INTEGER,  POINTER, DIMENSION(:)            :: IDT_H2O => NULL()

 REAL(DP), POINTER, DIMENSION(:,:,:)        :: FPR     => NULL()

 TYPE (PTR_3D_ARRAY), POINTER, DIMENSION(:) :: S_crit_in
 REAL(dp), POINTER, DIMENSION(:,:,:) :: xwat => NULL() 
 REAL(dp), POINTER, DIMENSION(:,:,:) :: xsat => NULL() 

 REAL(dp), POINTER, DIMENSION(:,:,:) :: Scrit   => NULL() 
 REAL(dp), POINTER, DIMENSION(:,:,:) :: WPARC   => NULL() !mz_vk_20170703

#ifdef MESSYTENDENCY
 integer                        :: my_handle
#endif

 INTEGER, PARAMETER :: ntrac_tot = 25 !ntrac_cl+2
 INTEGER, PARAMETER :: ntrac_cl = 23
 CHARACTER(LEN=8), DIMENSION(ntrac_cl) :: trac_cl = (/ &
      'COVER1  ', 'COVER2  ', &
      'CCCOV1DT', 'CCCOV2DT', 'CCCOV3DT', 'CCCOV4DT', 'CCCOV5DT', 'CCCOV   ', &
      'CCVOL1DT', 'CCVOL2DT', 'CCVOL3DT', 'CCVOL4DT', 'CCVOL5DT', 'CCVOL   ', &
      'CCLEN1DT', 'CCLEN2DT', 'CCLEN3DT', 'CCLEN4DT', 'CCLEN5DT', 'CCLEN   ', &
      'CCLENOU ',                                                             &
      'CCIWC   ', 'CCICNC  ' /)
 INTEGER, DIMENSION(ntrac_cl)  :: idt_cl = 0

 ! Handling of random numbers
 INTEGER, PARAMETER :: START_SEED = 260115  ! start seed
 INTEGER  :: RNDID = 0  ! pseudo random number generator
 INTEGER  :: RND_METHOD = 2  ! Mersenne Twister (machine independent)
 REAL(dp), DIMENSION(:), ALLOCATABLE :: HARVEST

!----------------------------------------------------------------
 CONTAINS
!===============================================================================
  SUBROUTINE  cloud_initialize

    ! CLOUD MODULE ROUTINE (ECHAM-5 INTERFACE)
    !
    ! INITIALIZATION OF GLOBAL VARIABLES FROM NAMELIST
    ! IN PARALLEL ENVIRONMENT
    ! INITIALIZATION OF XTSURF SPECIFIC EVENTS FROM NAMELIST
    ! IN PARALLEL ENVIRONMENT
    ! 
    ! Author: H. Tost, MPICH, 22-10-2004

    USE messy_main_mpi_bi,     ONLY: p_parallel_io, p_io, p_bcast
    USE messy_main_blather_bi, ONLY: error_bi, info_bi
    USE messy_main_tools,      ONLY: find_next_free_unit
    USE messy_cloud_mem,       ONLY: ncdnc, nicnc, du_nfrac_thr
    USE messy_cloud_lohmann10,   ONLY: cloud_read_nml_ctrl_l10 &
                                     , cdncmin_L10 => cdncmin, limm_BN09
    USE messy_cloud_kuebbeler14, ONLY: cloud_read_nml_ctrl_k14 &
                                     , cdncmin_K14 => cdncmin, nfrzaer, nexti &
                                     , ninp, max_ninp, inp_properties &
                                     , vervel_scheme &
                                     , scale_v_ls_cirrus &
                                     , scale_v_tke_cirrus &
                                     , scale_v_orogw_cirrus &
                                     , l_ic2snow

    IMPLICIT NONE

    ! LOCAL
    INTEGER     :: iou         ! I/O unit
    INTEGER     :: status      ! status
    INTEGER     :: i
    CHARACTER(LEN=70)           :: str=''
    CHARACTER(LEN=*), PARAMETER :: substr='cloud_initialize'

    ! INITIALIZE MAIN-CTRL
    status = 1
    IF (p_parallel_io) THEN
       iou = find_next_free_unit(100,200)
       ! *** CALL CLOUD CORE ROUTINE:
       CALL cloud_read_nml_ctrl(iou, status)
       if (status/=0) CALL error_bi('CLOUD INIT', substr)
    END IF
   
    CALL info_bi ('CLOUD active', substr)
    
    CALL P_BCAST(cloud_param, p_io)
    IF (CLOUD_PARAM == 1)   &
      CALL INFO_BI ( 'Standard ECHAM5 CLOUD scheme selected', substr)
    IF (CLOUD_PARAM == 2)   &
      CALL INFO_BI ( &
      'ECHAM5 - CDNC CLOUD scheme (Lohmann et al. 2007) selected', substr)
    IF (CLOUD_PARAM == 3)   &
      CALL INFO_BI ( &
      'ECHAM5 - CDNC/ICNC CLOUD scheme (Lohmann et al. 2007) selected', substr)
    IF (CLOUD_PARAM == 4)   &
      CALL INFO_BI ( &
      'ECHAM5 - CDNC/ICNC CLOUD scheme (Lohmann et al. 2010) selected', substr)
    IF (CLOUD_PARAM == 5)   &
      CALL info_bi ( &
      'ECHAM5 - CDNC/ICNC CLOUD scheme (Kuebbeler et al. 2014) selected' &
      , substr)
    IF (CLOUD_PARAM == 6)   &
      CALL info_bi ( &
      'ECHAM5 - CDNC/ICNC CLOUD scheme (Lohmann et al. 2010) and CCMod selected', substr)

    CALL P_BCAST(ncdnc, p_io)
    IF ( ncdnc == 1 )       &
      CALL INFO_BI ( &
      'CDNC activation scheme (Lin & Leaitch) selected', substr)
    IF ( ncdnc == 2 )       &
      CALL INFO_BI ( &
      'CDNC activation scheme (Abdul-Razzak & Ghan) selected', substr)
    IF ( ncdnc == 3 )       &                          
      CALL INFO_BI ( &                            
      'CDNC activation scheme (UAF) selected', substr)         

    CALL P_BCAST(nicnc, p_io)
    IF ( nicnc == 0 )       &
      CALL INFO_BI ( &
      'ECHAM5 - no ICNC calculations selected', substr)
    IF ( nicnc == 1 )       &
      CALL INFO_BI ( &
      'ECHAM5 - ICNC according to Lohmann selected', substr)
    IF ( nicnc == 2 )       &
      CALL INFO_BI ( &
      'ECHAM5 - ICNC according to Kaercher and Lohmann selected', substr)
    IF ( nicnc == 3 )       &
      CALL INFO_BI ( &
      'ECHAM5 - ICNC according to Barahona and Nenes selected', substr)
    CALL P_BCAST(lcover, p_io)

    CALL p_bcast(rset_ccraut%l, p_io)
    CALL p_bcast(rset_ccraut%v, p_io)
    CALL p_bcast(rset_ccsaut%l, p_io)
    CALL p_bcast(rset_ccsaut%v, p_io)
    CALL p_bcast(rset_cauloc%l, p_io)
    CALL p_bcast(rset_cauloc%v, p_io)
    CALL p_bcast(rset_csatsc%l, p_io)
    CALL p_bcast(rset_csatsc%v, p_io)
    CALL p_bcast(rset_crhsc%l, p_io)
    CALL p_bcast(rset_crhsc%v, p_io)

    IF ((CLOUD_PARAM == 4) .OR. (CLOUD_PARAM == 6)) THEN ! Lohmann (2010)
       IF (p_parallel_io) THEN
          iou = find_next_free_unit(100,200)
          CALL cloud_read_nml_ctrl_l10(iou, status, modstr)
          if (status/=0) CALL ERROR_BI('CTRL_L10 namelist error',substr)
       END IF
       CALL P_BCAST(cdncmin_L10, p_io)
       CALL P_BCAST(limm_BN09,   p_io)
    END IF

    IF (CLOUD_PARAM == 5) THEN ! Kuebbeler (2014)
       IF (p_parallel_io) THEN
          iou = find_next_free_unit(100,200)
          CALL cloud_read_nml_ctrl_k14(iou, status, modstr)
          if (status/=0) CALL ERROR_BI('CTRL_K14 namelist error', substr)
       END IF
       CALL P_BCAST(cdncmin_K14,   p_io)
       CALL P_BCAST(nfrzaer,       p_io)
       DO i = 1, max_ninp
          CALL P_BCAST(inp_properties(i)%name,     p_io)
          CALL P_BCAST(inp_properties(i)%Scrit,    p_io)
          CALL P_BCAST(inp_properties(i)%f_active, p_io)
       END DO
       CALL P_BCAST(nexti,         p_io)
       CALL P_BCAST(vervel_scheme, p_io)
       CALL P_BCAST(scale_v_ls_cirrus, p_io)
       CALL P_BCAST(scale_v_tke_cirrus, p_io)
       CALL P_BCAST(scale_v_orogw_cirrus, p_io)
       CALL P_BCAST(l_ic2snow,     p_io)
       CALL P_BCAST(du_nfrac_thr,  p_io)

       IF (nicnc == 2) THEN

          ! Check number of input INP properties
          IF (nfrzaer.gt.1) THEN
             DO i = 1, max_ninp
                IF (TRIM(inp_properties(i)%name) == '') EXIT
             END DO
             ninp = i - 1
             IF (ninp.ne.nfrzaer - 1) &
                  CALL ERROR_BI('CTRL_K14: number of inp_properties ' // &
                       'must be equal nfrzaer -1',substr)
          ELSE
             ninp = 0
          END IF

          ! Check unique and monotonically increasing Scrit values
          DO i = 2, ninp - 1
             IF (inp_properties(i)%Scrit <= inp_properties(i-1)%Scrit .OR. &
                 inp_properties(i)%Scrit >= inp_properties(i+1)%Scrit) &
                CALL ERROR_BI('CTRL_K14: input Scrit values must be ' // &
                            'strictly monotonically increasing', substr)
             IF (inp_properties(i)%Scrit < 1.2_dp) &
                CALL ERROR_BI('CTRL_K14: Scrit < 1.2 not supported ' // &
                            '(conflicting with 1.1<Scrit<1.2 range for DUdep' &
                            , substr)
          ENDDO

          ! Write summary
          WRITE(str,'(A,1X,F4.1,1X,A)') 'Minimum CDNC',cdncmin_K14/1.e6,'cm-3'
          CALL INFO_BI(str, substr)
          WRITE(str,'(I1,1X,A)') nfrzaer, 'freezing mode(s) in cirrus scheme:'
          CALL INFO_BI( str, substr)
          DO i = 1, ninp
             WRITE(str,'(A,A,1X,A,F5.2,1X,A,F6.3)') &
                  '  - ', TRIM(inp_properties(i)%name), &
                  'Scrit=',inp_properties(i)%Scrit,     &
                  'f_act=',inp_properties(i)%f_active
             CALL INFO_BI(str, substr)
          END DO
          CALL INFO_BI('  - homogeneous freezing', substr)
          IF (nexti == 1) &
               CALL INFO_BI( &
                    'Pre-existing ice crystals are taken into account', substr)
          IF (nexti == 0) &
               CALL INFO_BI(&
                    'Pre-existing ice crystals are not taken into account'&
                    , substr)
          WRITE(str, '(A,1X,F4.2)') &
               'Assuming dust dominance where N_DU/N_cm >=',du_nfrac_thr
          CALL INFO_BI( str, substr)
          CALL INFO_BI('Vertical velocity scheme:', substr)
          SELECT CASE (vervel_scheme)
          CASE(1)
             CALL INFO_BI('  Standard ECHAM5: large-scale + TKE', substr)
          CASE(2)
             CALL INFO_BI( &
                  '  Kuebbeler et al. (2014): large-scale + TKE/orogw', substr)
          CASE(3)
             CALL INFO_BI( &
                  '  Joos et al. (2008): large-scale + TKE + orogw', substr)
          CASE(4)
             CALL INFO_BI( &
                  '  Penner et al. (2018): large-scale + Laplace distribution'&
                  , substr)
          END SELECT
          WRITE(str, '(A,1X,F4.2)') &
               'Scaling factor for large-scale vertical velocity in cirrus =', &
               scale_v_ls_cirrus
          CALL INFO_BI( str, substr)
          WRITE(str, '(A,1X,F4.2)') &
               'Scaling factor for TKE vertical velocity in cirrus =', &
               scale_v_tke_cirrus
          CALL INFO_BI(str, substr)
          WRITE(str, '(A,1X,F4.2)') &
               'Scaling factor for orographic vertical velocity in cirrus =', &
               scale_v_orogw_cirrus
          CALL INFO_BI( str, substr)
          IF (l_ic2snow) &
               CALL INFO_BI( &
                    'Ice crystals larger than 100um are transformed to snow'&
                    , substr)
       ENDIF
    END IF

    ! INITIALIZE COUPLING-CONTROL
    IF (p_parallel_io) THEN
       iou = find_next_free_unit(100,200)
       CALL cloud_read_nml_cpl(status, iou)
       IF (status /= 0) CALL error_bi('cloud COUPLING INIT', substr)
    END IF  
    
    CALL P_BCAST(l_cdnc_calc, p_io)
    CALL P_BCAST(i_cdnc_calc, p_io)
    CALL P_BCAST(i_cdnc_cpl, p_io)
    CALL P_BCAST(aer_stream, p_io)
!!$    CALL P_BCAST(l_tf, p_io)
    CALL P_BCAST(sup_sat_scheme, p_io)

    !WISO++
    IF ( l_wiso .AND. (CLOUD_PARAM /= 1) ) THEN
       CALL error_bi('H2OISO only implemented for CLOUD_PARAM=1' ,substr)
    END IF
    !WISO--
    
  END SUBROUTINE cloud_initialize
! ==============================================================================

  SUBROUTINE cloud_read_nml_cpl(status, iou)

    ! cloud MODULE ROUTINE (ECHAM-5 INTERFACE, PRIVATE)
    !
    ! read namelist for 'coupling' to ECHAM5
    !
    ! Author: H. Tost, MPICH, March 2004

    USE messy_main_tools, ONLY: read_nml_open, read_nml_check, read_nml_close

    IMPLICIT NONE

    NAMELIST /CPL/ l_cdnc_calc, i_cdnc_calc, i_cdnc_cpl, aer_stream, sup_sat_scheme !!$, l_tf

    ! I/O
    INTEGER, INTENT(OUT) :: status     ! error status
    INTEGER, INTENT(IN)  :: iou        ! I/O unit

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr='cloud_read_nml_cpl'
    LOGICAL              :: lex      ! file exists ?
    INTEGER              :: fstat    ! file status

    status = 1
    l_cdnc_calc = .FALSE.

    ! INITIALIZE NAMELIST VARIABLES
 
    CALL read_nml_open(lex, substr, iou, 'CPL', modstr)
    IF (.not.lex) RETURN    ! <modstr>.nml does not exist

    READ(iou, NML=CPL, IOSTAT=fstat)
    CALL read_nml_check(fstat, substr, iou, 'CPL', modstr)
    IF (fstat /= 0) RETURN  ! error while reading namelist
  
    CALL read_nml_close(substr, iou, modstr)
    status = 0 ! NO ERROR

  END SUBROUTINE cloud_read_nml_cpl

! ==============================================================================

  SUBROUTINE cloud_new_tracer

    ! ECHAM5/MESSy
    USE messy_main_tracer_mem_bi,   ONLY: GPTRSTR, LGTRSTR
    ! MESSy
    USE messy_main_tracer,          ONLY: new_tracer, get_tracer,        &
                                          set_tracer,                    &
                                          CLOUD, ON, OFF, NUMBERDENSITY, &
                                          CONCENTRATION,                 &
                                          I_ADVECT, I_CONVECT,           &
                                          I_VDIFF,                       &
                                          I_DRYDEP, I_SEDI,              &
                                          I_SCAV, I_MIX,                 &
                                          R_MOLARMASS
    USE messy_main_tracer_tools_bi, ONLY: tracer_halt
    USE messy_cloud_mem,            ONLY: idt_cdnc, idt_icnc, ncdnc
    !WISO++
    USE messy_main_data_wiso_bi,    ONLY: kphase, mwiso, idiso
    !WISO--
 
    INTEGER :: js, status, idt
    INTEGER :: ji
    CHARACTER(LEN=*), DIMENSION(2), PARAMETER :: setname = (/GPTRSTR, LGTRSTR/)
    CHARACTER(LEN=*), PARAMETER :: substr = 'cloud_new_tracer'
#ifdef MESSYTENDENCY
    INTEGER :: jphase, jwiso
#endif
    
#ifdef MESSYTENDENCY
    my_handle = mtend_get_handle(modstr)
    CALL mtend_register (my_handle, mtend_id_t)
    CALL mtend_register (my_handle, mtend_id_q)
    CALL mtend_register (my_handle, mtend_id_xl)
    CALL mtend_register (my_handle, mtend_id_xi)

    !WISO++
    IF (l_wiso) THEN
       ! cloud_new_tracer must be called after h2oiso_new_tracer
       DO jphase = 1, kphase
          DO jwiso = 1, mwiso
             CALL mtend_register(my_handle, idiso(jphase, jwiso))
          END DO
       END DO
    END IF
    !WISO--    
#endif

    if ((cloud_param > 1) .and. (ncdnc > 0) ) THEN

      DO js=1,2

        IF (js == 1) THEN
          IF (.NOT.lcloud_gp) CYCLE
        END IF

        IF (js == 2) THEN
          IF (.NOT.lcloud_lg) CYCLE
        END IF

        CALL get_tracer(status, setname(js), 'CDNC',idx=idt_CDNC)
        IF (status /=0)  THEN
          CALL new_tracer(status, setname(js), 'CDNC','cloud',    &
                          quantity = numberdensity, unit='1/mol', &
                          medium = CLOUD, idx = idt_CDNC)
          CALL tracer_halt(substr,status)
          idt = idt_CDNC
#ifdef MESSYTENDENCY
          CALL mtend_register (my_handle, idt_CDNC) 
#endif
          CALL set_tracer(status, setname(js), idt, I_advect          , ON)
          CALL tracer_halt(substr, status)
          CALL set_tracer(status, setname(js), idt, I_convect         , ON)
          CALL tracer_halt(substr, status)
          CALL set_tracer(status, setname(js), idt, I_vdiff           , ON)
          CALL tracer_halt(substr, status)
          CALL set_tracer(status, setname(js), idt, I_scav            , OFF)
          CALL tracer_halt(substr, status)
          CALL set_tracer(status, setname(js), idt, I_drydep          , OFF)
          CALL tracer_halt(substr, status)
          CALL set_tracer(status, setname(js), idt, I_sedi            , OFF)
          CALL tracer_halt(substr, status)
          CALL set_tracer(status, setname(js), idt, I_mix             , OFF)
          CALL tracer_halt(substr, status)
          CALL set_tracer(status, setname(js), idt, R_molarmass       , 1._dp) 
          CALL tracer_halt(substr, status)    
        ENDIF

      END DO

    ENDIF
    if ((cloud_param >= 3) .and. (ncdnc > 0) ) THEN

      DO js=1,2

        IF (js == 1) THEN
          IF (.NOT.lcloud_gp) CYCLE
        END IF

        IF (js == 2) THEN
          IF (.NOT.lcloud_lg) CYCLE
        END IF

        CALL get_tracer(status, setname(js), 'ICNC',idx=idt_ICNC)
        IF (status /=0)  THEN
          CALL new_tracer(status, setname(js), 'ICNC','cloud',    &
                          quantity = numberdensity, unit='1/mol', &
                          medium = CLOUD, idx = idt_ICNC)
          CALL tracer_halt(substr,status)
          idt = idt_ICNC
#ifdef MESSYTENDENCY
          CALL mtend_register (my_handle, idt_ICNC) 
#endif
          CALL set_tracer(status, setname(js), idt, I_advect          , ON)
          CALL tracer_halt(substr, status)
          CALL set_tracer(status, setname(js), idt, I_convect         , ON)
          CALL tracer_halt(substr, status)
          CALL set_tracer(status, setname(js), idt, I_vdiff           , ON)
          CALL tracer_halt(substr, status)
          CALL set_tracer(status, setname(js), idt, I_scav            , OFF)
          CALL tracer_halt(substr, status)
          CALL set_tracer(status, setname(js), idt, I_drydep          , OFF)
          CALL tracer_halt(substr, status)
          CALL set_tracer(status, setname(js), idt, I_sedi            , OFF)
          CALL tracer_halt(substr, status)
          CALL set_tracer(status, setname(js), idt, I_mix             , OFF)
          CALL tracer_halt(substr, status)
          CALL set_tracer(status, setname(js), idt, R_molarmass       , 1._dp) 
          CALL tracer_halt(substr, status)   
        END IF
      END DO

    ENDIF

    if ((cloud_param == 6) .and. (ncdnc > 0) ) THEN

      DO js=1,2

        IF (js == 1) THEN
          IF (.NOT.lcloud_gp) CYCLE
        END IF

        IF (js == 2) THEN
          IF (.NOT.lcloud_lg) CYCLE
        END IF

        DO ji=1,ntrac_cl-2 
          CALL get_tracer(status, setname(js), trac_cl(ji), idx=idt_cl(ji))
          IF (status /=0)  THEN
            CALL new_tracer(status, setname(js), trac_cl(ji),'cloud',    &
                            quantity = 0, unit='', &
                            medium = CLOUD, idx = idt_cl(ji))
            CALL tracer_halt(substr,status)
            idt = idt_cl(ji)
#ifdef MESSYTENDENCY
            CALL mtend_register (my_handle, idt_cl(ji)) 
#endif
            CALL set_tracer(status, setname(js), idt, I_advect          , ON)
            CALL tracer_halt(substr, status)
            CALL set_tracer(status, setname(js), idt, I_convect         , ON)
            CALL tracer_halt(substr, status)
            CALL set_tracer(status, setname(js), idt, I_vdiff           , ON)
            CALL tracer_halt(substr, status)
            CALL set_tracer(status, setname(js), idt, I_scav            , OFF)
            CALL tracer_halt(substr, status)
            CALL set_tracer(status, setname(js), idt, I_drydep          , OFF)
            CALL tracer_halt(substr, status)
            CALL set_tracer(status, setname(js), idt, I_sedi            , OFF)
            CALL tracer_halt(substr, status)
            CALL set_tracer(status, setname(js), idt, I_mix             , OFF)
            CALL tracer_halt(substr, status)
            CALL set_tracer(status, setname(js), idt, R_molarmass       , 1._dp)
            CALL tracer_halt(substr, status)
          END IF
        END DO


        CALL get_tracer(status, setname(js), 'CCIWC',idx=idt_cl(22))
        IF (status /=0)  THEN
          CALL new_tracer(status, setname(js), 'CCIWC','cloud',    &
                          quantity = concentration, unit='kg/kg', &
                          medium = CLOUD, idx = idt_cl(22))
          CALL tracer_halt(substr,status)
          idt = idt_cl(22)
#ifdef MESSYTENDENCY
          CALL mtend_register (my_handle, idt_cl(22)) 
#endif
          CALL set_tracer(status, setname(js), idt, I_advect          , ON)
          CALL tracer_halt(substr, status)
          CALL set_tracer(status, setname(js), idt, I_convect         , ON)
          CALL tracer_halt(substr, status)
          CALL set_tracer(status, setname(js), idt, I_vdiff           , ON)
          CALL tracer_halt(substr, status)
          CALL set_tracer(status, setname(js), idt, I_scav            , OFF)
          CALL tracer_halt(substr, status)
          CALL set_tracer(status, setname(js), idt, I_drydep          , OFF)
          CALL tracer_halt(substr, status)
          CALL set_tracer(status, setname(js), idt, I_sedi            , OFF)
          CALL tracer_halt(substr, status)
          CALL set_tracer(status, setname(js), idt, I_mix             , OFF)
          CALL tracer_halt(substr, status)
          CALL set_tracer(status, setname(js), idt, R_molarmass       , 1._dp)
          CALL tracer_halt(substr, status)
        END IF

        CALL get_tracer(status, setname(js), 'CCNI',idx=idt_cl(23))
        IF (status /=0)  THEN
          CALL new_tracer(status, setname(js), 'CCNI','cloud',    &
                          quantity = numberdensity, unit='1/mol', &
                          medium = CLOUD, idx = idt_cl(23))
          CALL tracer_halt(substr,status)
          idt = idt_cl(23)
#ifdef MESSYTENDENCY
          CALL mtend_register (my_handle, idt_cl(23)) 
#endif
          CALL set_tracer(status, setname(js), idt, I_advect          , ON)
          CALL tracer_halt(substr, status)
          CALL set_tracer(status, setname(js), idt, I_convect         , ON)
          CALL tracer_halt(substr, status)
          CALL set_tracer(status, setname(js), idt, I_vdiff           , ON)
          CALL tracer_halt(substr, status)
          CALL set_tracer(status, setname(js), idt, I_scav            , OFF)
          CALL tracer_halt(substr, status)
          CALL set_tracer(status, setname(js), idt, I_drydep          , OFF)
          CALL tracer_halt(substr, status)
          CALL set_tracer(status, setname(js), idt, I_sedi            , OFF)
          CALL tracer_halt(substr, status)
          CALL set_tracer(status, setname(js), idt, I_mix             , OFF)
          CALL tracer_halt(substr, status)
          CALL set_tracer(status, setname(js), idt, R_molarmass       , 1._dp)
          CALL tracer_halt(substr, status)
        END IF

     END DO

    ENDIF

  END SUBROUTINE cloud_new_tracer

!===============================================================================

  SUBROUTINE cloud_init_memory


! channel objects for cloud parameters are constructed

    USE messy_main_grid_def_mem_bi,  ONLY: nlon
    USE messy_main_blather_bi,       ONLY: start_message_bi, end_message_bi 
    USE messy_main_channel_error_bi, ONLY: channel_halt
    USE messy_main_channel_bi,       ONLY: GP_3D_MID, GP_2D_HORIZONTAL
    USE messy_main_channel,    ONLY: new_channel, new_channel_object &
                                   , new_attribute &
                                   , new_channel_object_reference
    USE messy_cloud_mem,       ONLY: ncdnc, reffl, swat,                    &
                                     qeva,  qnuc,  qfre, qmel, qacc, qaut,  &
                                     cdnc,  cdnc_burden,                    &
                                     cdnc_acc, cdnc_burden_acc, cloud_tm1,  &
                                     sice, reffi, icnc, icnc_burden,        &
                                     icnc_acc, icnc_burden_acc, random_2d,  &
                                     ccfh2o, ccfkme,                        &
                                     cccov_tot, cccov_tot2, cccov,          &
                                     ccvol, ccicnc, cciwc, ccreffi, cctau
    USE messy_cloud_kuebbeler14,  ONLY: vervel_scheme
    USE messy_main_constants_mem, ONLY: ceffmax

    IMPLICIT NONE


    ! LOCAL    
    CHARACTER(LEN=*), PARAMETER::substr='cloud_init_memory'   
    INTEGER :: status

    CALL start_message_bi(modstr,'MEMORY INITIALIZATION', substr)

    CALL new_channel(status, modstr, reprid=GP_3D_MID)
    CALL channel_halt(substr, status)

#ifdef ECHAM5
    CALL new_channel_object_reference(status, 'g3b', 'acdnc', modstr, 'acdnc')
    CALL channel_halt(substr, status)
#endif
#ifdef CESM1
    CALL new_channel_object_reference(status, 'CESM1', 'acdnc', modstr, 'acdnc')
#endif

#ifdef ECHAM5
    CALL new_channel_object_reference(status, &
         'g3b', 'aclc', modstr, 'aclc')
#endif
#ifdef CESM1
    CALL new_channel_object_reference(status, &
         'CESM1', 'aclc', modstr, 'aclc')
#endif
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'aclc', &
         'long_name', c='large scale cloud cover')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'aclc', 'units', c='-' )
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr, 'prec_cover', p3=preccover)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'prec_cover', &
         'long_name', c='large scale precipitation cloud cover' )
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'prec_cover', 'units', c='-' )
    CALL channel_halt(substr, status)


    CALL new_channel_object(status, modstr, 'rainflux', p3=frain)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'rainflux', &
         'long_name', c='large scale rain precipitation flux' )
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'rainflux', 'units', c='kg/(m^2 * s)' )
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr, 'snowflux', p3=fsnow)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'snowflux', &
         'long_name', c='large scale snow precipitation flux' )
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'snowflux', 'units', c='kg/(m^2 * s)' )
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr, 'rainflux_no', p3=frain_no)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'rainflux_no', &
         'long_name', &
         c='large scale rain precipitation flux without'//&
         &' cloud production of new rain' )
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'rainflux_no', &
         'units', c='kg/(m^2 * s)' )
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr, 'snowflux_no', p3=fsnow_no)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'snowflux_no', &
         'long_name', &
         c='large scale snow precipitation flux without'//&
         &' cloud production of new snow' )
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'snowflux_no', &
         'units', c='kg/(m^2 * s)' )
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr, 'rain_form', p3=mratep)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'rain_form', &
         'long_name', c='large scale rain formation inside cloud' )
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'rain_form', 'units', c='kg/kg' )
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr, 'snow_form', p3=mratesi)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'snow_form', &
         'long_name', c='large scale snow formation inside cloud' )
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'snow_form', 'units', c='kg/kg' )
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr, 'rain_evap', p3=mrevap)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'rain_evap', &
         'long_name', c='large scale rain evaporation' )
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'rain_evap', 'units', c='kg/kg' )
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr, 'snow_subl', p3=mssubl)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'snow_subl', &
         'long_name', c='large scale snow sublimation' )
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'snow_subl', 'units', c='kg/kg' )
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr, 'lwc', p3=mlwc)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'lwc', &
         'long_name', c='large scale cloud liquid water content' )
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'lwc', 'units', c='kg/kg' )
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr, 'iwc', p3=msic)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'iwc', &
         'long_name', c='large scale cloud snow/ice content' )
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'iwc', 'units', c='kg/kg' )
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr, 'condensation', p3=condens)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'condensation', &
         'long_name', c='condensate in cloud covered part of gridbox' )
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'condensation', 'units', c='-' )
    CALL channel_halt(substr, status)
 
    CALL new_channel_object(status, modstr, 'mimelt', p3=mimelt)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'mimelt', &
         'long_name', c='large scale frozen precipitation melting' )
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'mimelt', 'units', c='kg/(m^2*s)' )
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr, 'misedi', p3=misedi)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'misedi', &
         'long_name', c='large scale ice sedimentation' )
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'misedi', 'units', c='kg/kg' )
    CALL channel_halt(substr, status)

#ifndef MESSYTENDENCY
    CALL new_channel_object(status, modstr, 'sdtdt_cloud', p3=sdtdt_cloud)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'sdtdt_cloud', &
         'long_name', c='cloud temp tendency' )
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'sdtdt_cloud', 'units', c='K/s' )
    CALL channel_halt(substr, status)
#endif

    if (l_cdnc_calc) then

      CALL new_channel_object(status, modstr, 'sulfaer', p3=aersulf_0)
      CALL channel_halt(substr, status)
      CALL new_attribute(status, modstr, 'sulfaer', &
        'long_name', c='aerosol sulfate' )
      CALL channel_halt(substr, status)
      CALL new_attribute(status, modstr, 'sulfaer', 'units', c='mug m^-3' )
      CALL channel_halt(substr, status)

      if ( (i_cdnc_calc == 1) .or. (i_cdnc_calc == 9) ) then
        CALL new_channel_object(status, modstr, 'CDN1', p3=acdnc_1)
        CALL channel_halt(substr, status)
        CALL new_attribute(status, modstr, 'CDN1', 'long_name', &
          c='cloud droplet number concentration (E5 changing with time) ' )
        CALL channel_halt(substr, status)
        CALL new_attribute(status, modstr, 'CDN1', 'units', c='m^-3' )
        CALL channel_halt(substr, status)
      endif

      if ( (i_cdnc_calc == 2) .or. (i_cdnc_calc == 9) ) then
        CALL new_channel_object(status, modstr, 'CDN2', p3=acdnc_2)
        CALL channel_halt(substr, status)
        CALL new_attribute(status, modstr, 'CDN2', &
          'long_name', c='cloud droplet number concentration (Rotstayn) ' )
        CALL channel_halt(substr, status)
        CALL new_attribute(status, modstr, 'CDN2', 'units', c='m^-3' )
        CALL channel_halt(substr, status)
      endif
      
      if ( (i_cdnc_calc == 3) .or. (i_cdnc_calc == 9) ) then
        CALL new_channel_object(status, modstr, 'CDN3', p3=acdnc_3)
        CALL channel_halt(substr, status)
        CALL new_attribute(status, modstr, 'CDN3', &
          'long_name', c='cloud droplet number concentration (Jones) ' )
        CALL channel_halt(substr, status)
        CALL new_attribute(status, modstr, 'CDN3', 'units', c='m^-3' )
        CALL channel_halt(substr, status)
      endif

      if ( (i_cdnc_calc == 4) .or. (i_cdnc_calc == 9) ) then
        CALL new_channel_object(status, modstr, 'CDN4', p3=acdnc_4)
        CALL channel_halt(substr, status)
        CALL new_attribute(status, modstr, 'CDN4', &
          'long_name', c='cloud droplet number concentration (Menon) ' )
        CALL channel_halt(substr, status)
        CALL new_attribute(status, modstr, 'CDN4', 'units', c='m^-3' )
        CALL channel_halt(substr, status)
      end if
    END if

    if ( (i_cdnc_calc == 5) .or. &
         (i_cdnc_calc == 9) .or. &
         (ncdnc       == 2) ) then

      CALL new_channel_object(status, modstr, 'CDN5', p3=acdnc_5)
      CALL channel_halt(substr, status)
      CALL new_attribute(status, modstr, 'CDN5', &
        'long_name', c='cloud droplet number activation (ARG) ' )
      CALL channel_halt(substr, status)
      CALL new_attribute(status, modstr, 'CDN5', 'units', c='m^-3' )
      CALL channel_halt(substr, status)   
     
    end if

    if ( (i_cdnc_calc == 6) .or. &
         (i_cdnc_calc == 9) .or. &
         (ncdnc       == 1) ) then
      CALL new_channel_object(status, modstr, 'CDN6', p3=acdnc_6)
      CALL channel_halt(substr, status)
      CALL new_attribute(status, modstr, 'CDN6', &
        'long_name', c='cloud droplet number activation (Lin) ' )
      CALL channel_halt(substr, status)
      CALL new_attribute(status, modstr, 'CDN6', 'units', c='m^-3' )
      CALL channel_halt(substr, status)
      CALL new_channel_object(status, modstr, 'NP_POTACT', p3=np_pot_activ)
      CALL channel_halt(substr, status)
      CALL new_attribute(status, modstr, 'NP_POTACT', &
        'long_name', c='potentially activated aerosol particle number(Lin) ' )
      CALL channel_halt(substr, status)
      CALL new_attribute(status, modstr, 'NP_POTACT', 'units', c='m^-3' )
      CALL channel_halt(substr, status)
      
    end if
    
    if ( (i_cdnc_calc == 7) .or. &
         (i_cdnc_calc == 9) .or. &
         (ncdnc       == 3) ) then

      CALL new_channel_object(status, modstr, 'CDN7', p3=acdnc_7)
      CALL channel_halt(substr, status)
      CALL new_attribute(status, modstr, 'CDN7', &
       'long_name', c='cloud droplet number activation (UAF) ' )
      CALL channel_halt(substr, status)
      CALL new_attribute(status, modstr, 'CDN7', 'units', c='m^-3' )
      CALL channel_halt(substr, status)

    end if
    
    if (ncdnc > 0) THEN
      CALL new_channel_object(status, modstr, 'SWAT', p3=swat)
      CALL channel_halt(substr, status)
      CALL new_attribute(status, modstr, 'SWAT', &
        'long_name', c='super saturation over water' )
      CALL channel_halt(substr, status)
      CALL new_attribute(status, modstr, 'SWAT', 'units', c='-' )
      CALL channel_halt(substr, status)

      CALL new_channel_object(status, modstr, 'REFFL', p3=REFFL &
           , lrestreq=.TRUE.)
      CALL channel_halt(substr, status)
      CALL new_attribute(status, modstr, 'REFFL', &
        'long_name', c='cloud droplet effective radius' )
      CALL channel_halt(substr, status)
      CALL new_attribute(status, modstr, 'REFFL', 'units', c='um' )
      CALL channel_halt(substr, status)
      REFFL = 40._dp

      CALL new_channel_object(status, modstr, 'QNUC', p3=QNUC)
      CALL channel_halt(substr, status)
      CALL new_attribute(status, modstr, 'QNUC', &
        'long_name', c='cloud droplet nucleation rate' )
      CALL channel_halt(substr, status)
      CALL new_attribute(status, modstr, 'QNUC', 'units', c='m-3 s-1' )
      CALL channel_halt(substr, status)

      CALL new_channel_object(status, modstr, 'QAUT', p3=QAUT)
      CALL channel_halt(substr, status)
      CALL new_attribute(status, modstr, 'QAUT', &
        'long_name', c='cloud droplet autoconversion rate' )
      CALL channel_halt(substr, status)
      CALL new_attribute(status, modstr, 'QAUT', 'units', c='m-3 s-1' )
      CALL channel_halt(substr, status)

      CALL new_channel_object(status, modstr, 'QACC', p3=QACC)
      CALL channel_halt(substr, status)
      CALL new_attribute(status, modstr, 'QACC', &
        'long_name', c='cloud droplet accretion rate' )
      CALL channel_halt(substr, status)
      CALL new_attribute(status, modstr, 'QACC', 'units', c='m-3 s-1' )
      CALL channel_halt(substr, status)

      CALL new_channel_object(status, modstr, 'QFRE', p3=QFRE)
      CALL channel_halt(substr, status)
      CALL new_attribute(status, modstr, 'QFRE', &
        'long_name', c='cloud droplet freezing rate' )
      CALL channel_halt(substr, status)
      CALL new_attribute(status, modstr, 'QFRE', 'units', c='m-3 s-1' )
      CALL channel_halt(substr, status)

      CALL new_channel_object(status, modstr, 'QEVA', p3=QEVA)
      CALL channel_halt(substr, status)
      CALL new_attribute(status, modstr, 'QEVA', &
        'long_name', c='cloud droplet evaporation rate' )
      CALL channel_halt(substr, status)
      CALL new_attribute(status, modstr, 'QEVA', 'units', c='m-3 s-1' )
      CALL channel_halt(substr, status)

      CALL new_channel_object(status, modstr, 'QMEL', p3=QMEL)
      CALL channel_halt(substr, status)
      CALL new_attribute(status, modstr, 'QMEL', &
        'long_name', c='cloud droplet melting (of ice) rate' )
      CALL channel_halt(substr, status)
      CALL new_attribute(status, modstr, 'QMEL', 'units', c='m-3 s-1' )
      CALL channel_halt(substr, status)

      CALL new_channel_object(status, modstr, 'CDNC', p3=CDNC)
      CALL channel_halt(substr, status)
      CALL new_attribute(status, modstr, 'CDNC', &
        'long_name', c='cloud droplet number concentration (whole grid box)' )
      CALL channel_halt(substr, status)
      CALL new_attribute(status, modstr, 'CDNC', 'units', c='m-3' )
      CALL channel_halt(substr, status)
      
      CALL new_channel_object(status, modstr, 'CDNC_acc', p3=CDNC_acc)
      CALL channel_halt(substr, status)
      CALL new_attribute(status, modstr, 'CDNC_acc', &
        'long_name', c='cloud droplet number concentration (in cloud)')
      CALL channel_halt(substr, status)
      CALL new_attribute(status, modstr, 'CDNC_acc', 'units', c='m-3' )
      CALL channel_halt(substr, status)
      
      CALL new_channel_object(status, modstr, 'CDNC_burden', p2=CDNC_burden, &
        reprid=GP_2d_horizontal)
      CALL channel_halt(substr, status)
      CALL new_attribute(status, modstr, 'CDNC_burden', &
        'long_name', c='cloud droplet number concentration column burden (whole grid box)' )
      CALL channel_halt(substr, status)
      CALL new_attribute(status, modstr, 'CDNC_burden', 'units', c='m-2' )
      CALL channel_halt(substr, status)

      CALL new_channel_object(status, modstr, 'CDNC_burden_acc', &
        p2=CDNC_burden_acc, reprid=GP_2d_horizontal)
      CALL channel_halt(substr, status)
      CALL new_attribute(status, modstr, 'CDNC_burden_acc', &
        'long_name', &
        c='cloud droplet number concentration column burden (in cloud only)')
      CALL channel_halt(substr, status)
      CALL new_attribute(status, modstr, 'CDNC_burden_acc', 'units', c='m-2' )
      CALL channel_halt(substr, status)

      CALL new_channel_object(status, modstr, 'cloud_tm1', p3=cloud_tm1 &
           , lrestreq = .TRUE.) ! op_mr_20161202
      CALL channel_halt(substr, status)
      CALL new_attribute(status, modstr, 'cloud_tm1', &
        'long_name', c='cloud cover (previous time step)')
      CALL channel_halt(substr, status)
      CALL new_attribute(status, modstr, 'cloud_tm1', 'units', c='-' )
      CALL channel_halt(substr, status)
    endif
    
    IF ((cloud_param >= 3).AND.(cloud_param <= 6)) THEN

      CALL new_channel_object(status, modstr, 'sice', p3=sice)
      CALL channel_halt(substr, status)
      CALL new_attribute(status, modstr, 'sice', &
        'long_name', c='ice supersaturation')
      CALL channel_halt(substr, status)
      CALL new_attribute(status, modstr, 'sice', 'units', c='-' )
      CALL channel_halt(substr, status)
      
      CALL new_channel_object(status, modstr, 'REFFI', p3=REFFI &
           , lrestreq=.TRUE.)
      CALL channel_halt(substr, status)
      CALL new_attribute(status, modstr, 'REFFI', &
        'long_name', c='ice crystal effective radius' )
      CALL channel_halt(substr, status)
      CALL new_attribute(status, modstr, 'REFFI', 'units', c='um' )
      CALL channel_halt(substr, status)
      REFFI = ceffmax

      CALL new_channel_object(status, modstr, 'ICNC', p3=ICNC)
      CALL channel_halt(substr, status)
      CALL new_attribute(status, modstr, 'ICNC', &
        'long_name', c='ice crystal number concentration (whole grid box)' )
      CALL channel_halt(substr, status)
      CALL new_attribute(status, modstr, 'ICNC', 'units', c='m-3' )
      CALL channel_halt(substr, status)
      
      CALL new_channel_object(status, modstr, 'ICNC_acc', p3=ICNC_acc)
      CALL channel_halt(substr, status)
      CALL new_attribute(status, modstr, 'ICNC_acc', &
        'long_name', c='ice crystal number concentration (in cloud)')
      CALL channel_halt(substr, status)
      CALL new_attribute(status, modstr, 'ICNC_acc', 'units', c='m-3' )
      CALL channel_halt(substr, status)
      
      CALL new_channel_object(status, modstr, 'ICNC_burden', p2=ICNC_burden, &
        reprid=GP_2d_horizontal)
      CALL channel_halt(substr, status)
      CALL new_attribute(status, modstr, 'ICNC_burden', &
        'long_name', c='ice crystal number concentration column burden (whole grid box)' )
      CALL channel_halt(substr, status)
      CALL new_attribute(status, modstr, 'ICNC_burden', 'units', c='m-2' )
      CALL channel_halt(substr, status)

      CALL new_channel_object(status, modstr, 'ICNC_burden_acc', &
        p2=ICNC_burden_acc, reprid=GP_2d_horizontal)
      CALL channel_halt(substr, status)
      CALL new_attribute(status, modstr, 'ICNC_burden_acc', &
        'long_name', &
        c='ice crystal number concentration column burden (in cloud only)')
      CALL channel_halt(substr, status)
      CALL new_attribute(status, modstr, 'ICNC_burden_acc', 'units', c='m-2' )
      CALL channel_halt(substr, status)
    ENDIF

    CALL new_channel_object(status, modstr, 'rhc', p3=rhc)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'rhc', &
         'long_name', c='critical relative humidity for natural clouds' )
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'rhc', 'units', c=' ' )
    CALL channel_halt(substr, status)

    IF (cloud_param == 6) THEN

      CALL new_channel_object(status, modstr, 'ccfh2o', p3=ccfh2o)
      CALL channel_halt(substr, status)
      CALL new_attribute(status, modstr, 'ccfh2o', &
            'long_name', c='water vapor emission of aviation' )
      CALL channel_halt(substr, status)
      CALL new_attribute(status, modstr, 'ccfh2o', 'units', c='kg/kg/s' )
      CALL channel_halt(substr, status)

      CALL new_channel_object(status, modstr, 'ccfkme', p3=ccfkme)
      CALL channel_halt(substr, status)
      CALL new_attribute(status, modstr, 'ccfkme', &
            'long_name', c='flight (slant) distance of aviation' )
      CALL channel_halt(substr, status)
      CALL new_attribute(status, modstr, 'ccfkme', 'units', c='m/s' )
      CALL channel_halt(substr, status)

      CALL new_channel_object(status, modstr, 'ccvol', p3=ccvol)
      CALL channel_halt(substr, status)
      CALL new_attribute(status, modstr, 'ccvol', &
            'long_name', c='contrail cirrus volume' )
      CALL channel_halt(substr, status)
      CALL new_attribute(status, modstr, 'ccvol', 'units', c='' )
      CALL channel_halt(substr, status)

      CALL new_channel_object(status, modstr, 'cccov', p3=cccov)
      CALL channel_halt(substr, status)
      CALL new_attribute(status, modstr, 'cccov', &
            'long_name', c='contrail cirrus coverage' )
      CALL channel_halt(substr, status)
      CALL new_attribute(status, modstr, 'cccov', 'units', c='' )
      CALL channel_halt(substr, status)

      CALL new_channel_object(status, modstr, 'cccov_tot', p2=cccov_tot, &
        reprid=GP_2d_horizontal)
      CALL channel_halt(substr, status)
      CALL new_attribute(status, modstr, 'cccov_tot', &
            'long_name', c='total contrail cirrus coverage overlapped' )
      CALL channel_halt(substr, status)
      CALL new_attribute(status, modstr, 'cccov_tot', 'units', c='' )
      CALL channel_halt(substr, status)

      CALL new_channel_object(status, modstr, 'cccov_tot2', p2=cccov_tot2, &
        reprid=GP_2d_horizontal)
      CALL channel_halt(substr, status)
      CALL new_attribute(status, modstr, 'cccov_tot2', &
            'long_name', c='total contrail cirrus coverage overlapped, tau>0.05' )
      CALL channel_halt(substr, status)
      CALL new_attribute(status, modstr, 'cccov_tot2', 'units', c='' )
      CALL channel_halt(substr, status)

      CALL new_channel_object(status, modstr, 'ccicnc', p3=ccicnc)
      CALL channel_halt(substr, status)
      CALL new_attribute(status, modstr, 'ccicnc', &
            'long_name', &
            c='ice crystal number concentration in contrail cirrus (in-cloud)' )
      CALL channel_halt(substr, status)
      CALL new_attribute(status, modstr, 'ccicnc', 'units', c='m-3' )
      CALL channel_halt(substr, status)

      CALL new_channel_object(status, modstr, 'cciwc', p3=cciwc)
      CALL channel_halt(substr, status)
      CALL new_attribute(status, modstr, 'cciwc', &
            'long_name', c='ice water content in contrail cirrus (in-cloud)' )
      CALL channel_halt(substr, status)
      CALL new_attribute(status, modstr, 'cciwc', 'units', c='kg/m-3' )
      CALL channel_halt(substr, status)

      CALL new_channel_object(status, modstr, 'ccreffi', p3=ccreffi)
      CALL channel_halt(substr, status)
      CALL new_attribute(status, modstr, 'ccreffi', &
            'long_name', c='ice crystal radius in contrail cirrus (in-cloud)' )
      CALL channel_halt(substr, status)
      CALL new_attribute(status, modstr, 'ccreffi', 'units', c='um' )
      CALL channel_halt(substr, status)

      CALL new_channel_object(status, modstr, 'cctau', p3=cctau)
      CALL channel_halt(substr, status)
      CALL new_attribute(status, modstr, 'cctau', &
            'long_name', c='optical depth of contrail cirrus (in-cloud)' )
      CALL channel_halt(substr, status)
      CALL new_attribute(status, modstr, 'cctau', 'units', c='' )
      CALL channel_halt(substr, status)

    END IF

    IF (cloud_param == 5 .AND. vervel_scheme == 4) THEN
       ! Initialize pseudo number streams
       CALL rnd_init_bi(RNDID, RND_METHOD, RND_MP_PIN, START_SEED)
       ! Allocate array for random numbers
       ALLOCATE(HARVEST(nlon))

       CALL new_channel_object(status, modstr, 'random_2d', p2=random_2d, &
            reprid=GP_2d_horizontal)
       CALL channel_halt(substr, status)
       CALL new_attribute(status, modstr, 'random_2d', &
            'long_name', c='two-dimensional array of random numbers')
       CALL channel_halt(substr, status)
       CALL new_attribute(status, modstr, 'random_2d', 'units', c='-')
       CALL channel_halt(substr, status)
    END IF

    CALL end_message_bi(modstr,'MEMORY INITIALIZATION', substr)

  END SUBROUTINE cloud_init_memory

!===============================================================================

  SUBROUTINE cloud_global_start

    USE messy_cloud_kuebbeler14, ONLY: vervel_scheme
    
    IMPLICIT NONE

    IF (cloud_param == 5 .AND. vervel_scheme == 4) THEN
       ! Generate array of random numbers
       CALL rnd_number_bi(RNDID, HARVEST)

       ! Distribute between -0.5 and 0.5
       HARVEST = HARVEST - 0.5_dp
    END IF
       
  END SUBROUTINE cloud_global_start

!===============================================================================
 
  SUBROUTINE cloud_free_memory

    USE messy_cloud_kuebbeler14, ONLY: vervel_scheme
    
    IMPLICIT NONE

    if (ALLOCATED(sulf_idt)) DEALLOCATE(sulf_idt)
    if (ALLOCATED(ss_idt))   DEALLOCATE(ss_idt)
    if (ALLOCATED(om_idt))   DEALLOCATE(om_idt)

    IF (cloud_param == 5 .AND. vervel_scheme == 4) THEN
       ! Handling of random numbers
       IF (ALLOCATED(HARVEST)) DEALLOCATE(HARVEST)
       CALL rnd_finish_bi(RNDID)
    END IF

  END SUBROUTINE cloud_free_memory

!============================================================================== 
   SUBROUTINE cloud_init_coupling
  
    USE messy_main_mpi_bi,           ONLY: p_parallel_io
    USE messy_main_blather_bi,       ONLY: start_message_bi, end_message_bi &
                                         , error_bi, info_bi
    USE messy_main_channel_error_bi, ONLY: channel_halt
    USE messy_main_channel_bi,       ONLY: GP_3D_MID, DC_GP                &
                                         , DIMID_LON, DIMID_LAT, DIMID_LEV &
                                         , gp_nseg, gp_start, gp_cnt &
                                         , gp_meml, gp_memu
    USE messy_main_tracer_mem_bi,    ONLY: GPTRSTR, ntrac => ntrac_gp, ti_gp
    USE messy_main_data_bi,          ONLY: tke, ltdiag

#if defined(CESM1)
    USE messy_cloud_ori,       ONLY: sucloud, tbetai0, ncctop, ccraut, ccsaut &
                                   , cauloc, csatsc, crhsc
    USE messy_main_grid_def_mem_bi, ONLY: nlev, nlevp1, nvclev, nn, lmidatm, vct
    USE messy_main_data_bi,         ONLY: modstr_base=>modstr

#endif
#if defined(ECHAM5)
    USE messy_cloud_ori,       ONLY: sucloud_1, tbetai0, ncctop, ccraut, ccsaut &
                                   , cauloc, csatsc, crhsc
    USE messy_main_grid_def_mem_bi, ONLY: nlev, nlevp1, nvclev, nn, lmidatm, vct
    USE messy_main_data_bi,         ONLY: modstr_base=>modstr, lcouple
#endif
    USE messy_cloud_mem,          ONLY: ncdnc,                    &
                                        cdncact_cv, nbcsol_strat, &
                                        ndusol_strat, nbcsol_cirrus, nbctagsol_cirrus, &
                                        ndusol_cirrus, nbcinsol, nbctaginsol,  &
                                        nduinsolai, nduinsolci,   &
                                        nocsolks,  nocsolas, nocsolcs, nocinsol,  &
                                        naerinsol, naersol,         &
                                        inucs, iaits, iaccs, icoas, &
                                        iaitm, iaccm, icoam, &
                                        iaiti, iacci, icoai, &
                                        B_cc, B_co,                   &
                                        gboxarea_2d,                  &
                                        ccfh2o_invent, ccfkme_invent, &
                                        twc_conv, conv_time, nfrzmod, &
                                        vervel_ls, vervel_tke, vervel_gw,     &
                                        vervel_p18, cdnc_insitu,              &
                                        Nice_preex, Nice_DUdep, Nice_DUimm,   &
                                        Nice_BC,    Nice_BCtag, Nice_homog,   &
                                        N_CIRRUS_TEMP, N_CIRRUS_BIN_IWC,      &
                                        N_CIRRUS_BIN_Nice, N_CIRRUS_BIN_Rice, &
                                        N_CIRRUS_BIN_RHi, N_CIRRUS_BIN_VEL,   &
                                        CIRRUS_TEMP,                          &
                                        CIRRUS_BIN_IWC,  CIRRUS_IBIN_IWC,     &
                                        CIRRUS_BIN_Nice, CIRRUS_IBIN_Nice,    &
                                        CIRRUS_BIN_Rice, CIRRUS_IBIN_Rice,    &
                                        CIRRUS_BIN_RHi,  CIRRUS_IBIN_RHi,     &
                                        CIRRUS_BIN_VEL,  CIRRUS_IBIN_VEL,     &
                                        CIRRUS_IWC, CIRRUS_Nice, CIRRUS_Rice, &
                                        CIRRUS_RHi_cloud, CIRRUS_RHi_clear,   &
                                        CIRRUS_vervel, CIRRUS_Nice_ML
    USE messy_cloud_cover,        ONLY: init_cover_ori
    USE messy_cloud_lohmann07,    ONLY: idt_ncs, idt_nas, idt_nks,          &
                                        idt_nci, idt_nai, idt_nki,          &
                                        idt_moccs, idt_mocas, idt_mocks,    &
                                        idt_mocki, idt_mbccs, idt_mbcas,    &
                                        idt_mbcks, idt_mbcki, idt_ms4ks,    &
                                        idt_ms4as, idt_ms4cs, idt_mduas,    &
                                        idt_mducs, idt_mssas, idt_msscs
    USE messy_main_channel,       ONLY: get_channel_object, new_channel,    &
                                        new_channel_object, new_attribute 
    USE messy_main_channel_dimensions, ONLY: new_dimension,                 &
                                             add_dimension_variable,        &
                                             add_dimension_variable_att

    USE messy_main_channel_repr,       ONLY: new_representation, AUTO &
                                           , set_representation_decomp &
                                           , IRANK, PIOTYPE_COL
    USE messy_main_constants_mem, ONLY: M_air
    USE messy_main_tools,         ONLY: int2str, match_wild
    USE messy_main_tracer,        ONLY: get_tracer,                         &
                                        AMOUNTFRACTION, NUMBERDENSITY,      &
                                        AEROSOL, OFF,                       &
                                        R_molarmass, R_aerosol_density,     &
                                        I_aerosol_mode, I_aerosol_sol


    IMPLICIT NONE

    INTEGER  :: status, jt, i, ji, jn, jm, ibin
    LOGICAL  :: laerosol(ntrac)
    REAL(dp) :: mmass(ntrac), aerdens(ntrac)
    CHARACTER(LEN=*), PARAMETER::substr='cloud_init_cpl'

    CHARACTER(LEN=2)             :: mod_number, str_num
    CHARACTER(LEN=15)            :: name, modname
    CHARACTER(LEN=1)             :: char1
    
    REAL(dp), DIMENSION(:,:,:,:),   POINTER ::  mem => NULL()
    INTEGER                               :: DIMID_NSOLMODE
    INTEGER                               :: DIMID_NTOTMODE
    INTEGER                               :: REPR_CLOUD_4D_KSOL
    INTEGER                               :: REPR_CLOUD_4D_KTOT

    ! For cirrus diagnostics
    INTEGER :: DIMID_CIRRUS_TEMP
    INTEGER :: DIMID_CIRRUS_BIN_IWC
    INTEGER :: DIMID_CIRRUS_BIN_Nice
    INTEGER :: DIMID_CIRRUS_BIN_Rice
    INTEGER :: DIMID_CIRRUS_BIN_RHi
    INTEGER :: DIMID_CIRRUS_BIN_VEL
    INTEGER :: REPR_CIRRUS_4D_IWC
    INTEGER :: REPR_CIRRUS_4D_Nice
    INTEGER :: REPR_CIRRUS_4D_Rice
    INTEGER :: REPR_CIRRUS_4D_RHi
    INTEGER :: REPR_CIRRUS_4D_VEL

    ! PARALLEL DECOMPOSITION
    INTEGER                          :: nseg = 0
    INTEGER, DIMENSION(:,:), POINTER :: start => NULL()
    INTEGER, DIMENSION(:,:), POINTER :: cnt   => NULL()
    INTEGER, DIMENSION(:,:), POINTER :: meml  => NULL()
    INTEGER, DIMENSION(:,:), POINTER :: memu  => NULL()
    LOGICAL :: l_use_psc = .FALSE.

    INTRINSIC ::  MAXVAL
    status      = 0
    laerosol(:) = .FALSE.
    mmass(:)    = 0._dp
    aerdens(:)  = 0._dp
    
    CALL start_message_bi(modstr,'COUPLING INITIALIZATION',substr)

    IF ( (cloud_param == 1) .or. (cloud_param == 2) .or. (cloud_param == 3) &
                 .or. ((cloud_param  >= 4) .AND. (cloud_param <= 6)) ) &
      THEN
#if defined(CESM1)
       CALL sucloud(status, nlev, nlevp1, nvclev, vct, nn, lmidatm)
#endif
#if defined(ECHAM5)
       CALL sucloud_1(status, nlev, nlevp1, nvclev, vct, &
                                       nn, lmidatm, lcouple, .FALSE.)
#endif
    endif
#if defined(ECHAM5)
    IF ( (cloud_param > 2)  ) CALL set_cloud_parameters(nlev, .FALSE.)
    ! Overwrite ccraut and ccsaut with the values provided via namelist
    if (rset_ccraut%l) ccraut = rset_ccraut%v
    if (rset_ccsaut%l) ccsaut = rset_ccsaut%v
    if (p_parallel_io) write(*,*) "Autoconversion rate = ",ccraut
    if (p_parallel_io) write(*,*) "Aggregation rate = ",ccsaut
    ! Overwrite cauloc, csatsc, and crhsc with the values provided via namelist
    IF (rset_cauloc%l) cauloc = rset_cauloc%v
    IF (rset_csatsc%l) csatsc = rset_csatsc%v
    IF (rset_crhsc%l)  crhsc = rset_crhsc%v
    IF (p_parallel_io) THEN
       WRITE(*,*) "rset_cauloc = ",rset_cauloc," cauloc = ",cauloc
       WRITE(*,*) "rset_csatsc = ",rset_csatsc," csatsc = ",csatsc
       WRITE(*,*) "rset_crhsc = ",rset_crhsc," crhsc = ",crhsc
    END IF
#endif

    if (status.eq.1) call error_bi(  &
        'problem with invalid vertical resolution' &
        ,'messy_cloud_ori.f90 : sucloud')
    if (status.eq.2) call error_bi(  &
      'problem with invalid combination of horizontal and vertical resolution'&
       , 'messy_cloud_ori.f90 : sucloud')
    if (status.eq.3) call error_bi( &
      'problem with invalid truncation','messy_cloud_ori.f90 : sucloud')   
    IF (p_parallel_io) &
      WRITE (*,*) 'highest level for condensation: ncctop= ',ncctop


    IF  (p_parallel_io) &
      WRITE(*,*) 'Initializing COVER lookup tables...'
    call init_cover_ori
    if (maxval(tbetai0(:,:)).le.0._dp) &
      CALL ERROR_BI ('init lookup tables for COVER', modstr)

    IF (p_parallel_io) &
         WRITE(*,*) 'Checking for PSC module ...'

!!#D psc +
    CALL get_channel_object(status, 'psc', 'PSC_region', p3=flt_pscreg)
    l_use_psc = (status == 0)
    USE_PSC = USE_PSC .OR. l_use_psc

    IF (l_use_psc)  THEN
       IF (p_parallel_io) WRITE(*,*) ' ... PSC active'
    ENDIF
!!#D psc -

!!#D msbm +
    CALL get_channel_object(status, 'msbm', 'STRAT_region', p3=flt_pscreg)
    l_use_psc = (status == 0)
    USE_PSC = USE_PSC .OR. l_use_psc
    IF (l_use_psc)  THEN
       IF (p_parallel_io) WRITE(*,*) ' ... MSBM active'
    END IF
!!#D msbm -

!!#D e4chem +
    CALL get_channel_object(status, 'e4chem', 'PSC_region', p3=flt_pscreg)
    l_use_psc = (status == 0)
    USE_PSC = USE_PSC .OR. l_use_psc
    IF (l_use_psc)  THEN
       IF (p_parallel_io) WRITE(*,*) ' ... E4CHEM active'
    END IF
!!#D e4chem -

    IF (.NOT. USE_PSC)  THEN
       IF (p_parallel_io) WRITE(*,*) ' ... none active'
    END IF

    ! NOTE: This should ultimately be replaced by local pointers set
    !       with get_channel_object, however, this needs to be done
    !       for ALL basemodels (and not all have channel objects already).
    IF (p_parallel_io) &
      WRITE(*,*) 'Checking for grid box parameters from vertical diffusion ...'
    IF (.NOT. ASSOCIATED(tke)) &
         call error_bi('tke is not associated', substr)

    IF (p_parallel_io) &
      WRITE(*,*) 'Checking for grid box parameters from base model ...'
    CALL get_channel_object(status, modstr_base, 'qtec', p3=pqtec)
    IF (status /= 0) &
      call error_bi('channel object for qtec not found', substr)
    CALL get_channel_object(status, modstr_base, 'xtecl', p3=pxtecl)
    IF (status /= 0) &
      call error_bi('channel object for xtecl not found', substr)
    CALL get_channel_object(status, modstr_base, 'xteci', p3=pxteci)
    IF (status /= 0) &
      call error_bi('channel object for xteci not found', substr)
    CALL get_channel_object(status, modstr_base, 'vmixtau', p3=pvmixtau)
    IF (status /= 0) &
      call error_bi('channel object for vmixtau not found', substr)    
    CALL get_channel_object(status, modstr_base, 'vdiffp', p3=pvdiffp)
    IF (status /= 0) &
      call error_bi('channel object for vdiffp not found', substr)
    

    CALL get_channel_object(status, modstr_base, 'xtecnl', p3=pxtecnl)
    IF (status /= 0) &
      call error_bi('channel object for xtecnl not found', substr)
    CALL get_channel_object(status, modstr_base, 'xtecni', p3=pxtecni)
    IF (status /= 0) &
      call error_bi('channel object for xtecni not found', substr)

    IF (p_parallel_io) &
      WRITE(*,*) 'Checking for convective parameters ...'

    CALL get_channel_object(status, 'convect', 'conv_type', p2=conv_type)
    IF (status /= 0) THEN
      call info_bi('channel object for conv_type not found',substr)
      call info_bi( &
        'assuming a convection type of zero for each grid column',substr)
    END IF

    CALL get_channel_object(status, 'convect', 'conv_bot', p2=conv_bot)
    IF (status /= 0) THEN
      call info_bi( &
        'channel object for bottom level of convection not found',substr)
      call info_bi( &
        'assuming a bottom level of convection of zero for each grid column'&
        ,substr)
    END IF

    CALL get_channel_object(status, 'convect', 'CAPE', p2=pcape)
    IF (status /= 0) THEN
      call info_bi('channel object for CAPE not found',substr)
      call info_bi('assuming a CAPE of zero for each grid column',substr)
    END IF

    IF (ltdiag) THEN
       CALL info_bi('looking for channel / object tdiag / PDIGA19','  ')
       CALL get_channel_object(status, 'tdiag', 'PDIGA19', p3=pdiga)
       CALL channel_halt(substr, status)
    END IF
    
!------------------------------------------------------------------------
! for cdnc coupling
    if ( (i_cdnc_calc < 1)   .or. &
         (i_cdnc_calc >= 10) .and.&
         (cloud_param == 1) ) l_cdnc_calc = .FALSE.

    do jt=1,ntrac
       IF ((ti_gp(jt)%tp%ident%medium==AEROSOL)) laerosol(jt)  = .true.
    enddo

    if (l_cdnc_calc .OR. ncdnc > 0) then

      CALL INFO_BI('Checking for aerosol compounds for coupling of aerosol species to cloud droplet numbers!', substr)
! sulphate
      sulf_count = 0
      do jt=1,ntrac
        IF ((ti_gp(jt)%tp%ident%medium==AEROSOL)) laerosol(jt)  = .true.
        if ( ( ti_gp(jt)%tp%ident%basename=='SO4' .or. &
               ti_gp(jt)%tp%ident%basename=='SO4mm' ) .and. &
              (laerosol(jt)) ) sulf_count = sulf_count + 1
      enddo

      if (sulf_count > 0) then
        allocate(sulf_idt(sulf_count))
        sulf_idt(:) = 0
        i=0
        do jt=1,ntrac
          if ( ( ti_gp(jt)%tp%ident%basename =='SO4'     .or.  &
                 ti_gp(jt)%tp%ident%basename =='SO4mm' ) .and. &
               ( laerosol(jt) ) ) then
            i=i+1
            CALL get_tracer(status, GPTRSTR, trim(ti_gp(jt)%tp%ident%basename),&
                            subname=trim(ti_gp(jt)%tp%ident%subname),          &
                            idx=sulf_idt(i) )
            fac_sulf =  ti_gp(jt)%tp%meta%cask_r(R_molarmass) / M_air
!!$        print*, sulf_idt, ti_gp(jt)%tp%ident%fullname, fac_sulf,  &
!!$                ti_gp(jt)%tp%meta%cask_r(R_molarmass), M_air
          endif
        enddo
      else
        if (i_cdnc_calc > 1) &
        CALL error_bi( &
        'No aerosol sulphate found, but required for Menon CDNC calculation!'&
        , substr)
      endif
!-------
! seasalt
      ss_count  = 0
      ss2_count = 0
      do jt=1,ntrac
        if ( ( ti_gp(jt)%tp%ident%basename=='SS' ) .and. &
             (laerosol(jt)) ) ss_count = ss_count + 1
        if ( ( ti_gp(jt)%tp%ident%basename=='Nap' ) .and. &
             (laerosol(jt)) ) ss2_count = ss2_count + 1
      enddo
      
      if (ss_count > 0) then
         allocate(ss_idt(ss_count))
         ss_idt(:) = 0
         i=0
         do jt=1,ntrac
            if ( ( ti_gp(jt)%tp%ident%basename =='SS' ) .and. &
                 ( laerosol(jt) ) ) then
               i=i+1
               CALL get_tracer(status, GPTRSTR, trim(ti_gp(jt)%tp%ident%basename),&
                    subname=trim(ti_gp(jt)%tp%ident%subname),          &
                    idx=ss_idt(i) )
               fac_ss =  ti_gp(jt)%tp%meta%cask_r(R_molarmass) / M_air
!!$        print*, ss_idt, ti_gp(jt)%tp%ident%fullname, fac_ss,  &
!!$                ti_gp(jt)%tp%meta%cask_r(R_molarmass), M_air
            endif
         enddo
      ENDIF
      i = 0
      If (ss2_count > 0) THEN
        allocate (ss2_idt(ss2_count,2))
        ss2_idt(:,:) = 0
        do jt=1,ntrac
          if ( ( trim(ti_gp(jt)%tp%ident%basename) =='Nap' ) .and. &
               ( laerosol(jt) ) ) then
            i = i+1
            CALL get_tracer(status, GPTRSTR, trim(ti_gp(jt)%tp%ident%basename),&
                            subname=trim(ti_gp(jt)%tp%ident%subname),          &
                            idx=ss2_idt(i,1) )
            kmod = ti_gp(jt)%tp%meta%cask_i(I_AEROSOL_MODE)
            DO ji=1,ntrac
              IF (ti_gp(ji)%tp%ident%medium /=AEROSOL) CYCLE
              IF (ti_gp(ji)%tp%meta%cask_i(I_AEROSOL_MODE) /= kmod ) CYCLE
              if (trim(ti_gp(ji)%tp%ident%basename) =='Clm' ) &
                CALL get_tracer(status, GPTRSTR, trim(ti_gp(ji)%tp%ident%basename),&
                                subname=trim(ti_gp(ji)%tp%ident%subname),          &
                                idx=ss2_idt(i,2) )
            END DO
          END if
        END do
     END If

      IF (ss_count + ss2_count == 0) THEN
        if (i_cdnc_calc == 4) then
          CALL INFO_BI( &
            'No aerosol seasalt found, but required for Menon CDNC calculation!', substr)
          CALL INFO_BI( &
            'Calculation may be possible, but is likely to be wrong (incomplete)!', substr)
        endif
      endif
      
!---------
! organic matter
      om_count = 0
      do jt=1,ntrac
        IF (.NOT. laerosol(jt)) CYCLE
        IF (ti_gp(jt)%tp%meta%cask_i(I_AEROSOL_SOL) == OFF) CYCLE
        if (ti_gp(jt)%tp%ident%basename=='OC' ) &
          om_count = om_count + 1
        IF (MATCH_WILD( 'WSOC*', ti_gp(jt)%tp%ident%basename) ) &
          om_count = om_count + 1
      enddo
      
      if (om_count > 0) then
        allocate(om_idt(om_count))
        om_idt(:) = 0
        i=0
        do jt=1,ntrac
          IF (.NOT. laerosol(jt)) CYCLE
          IF (ti_gp(jt)%tp%meta%cask_i(I_AEROSOL_SOL) == OFF) CYCLE
          if (ti_gp(jt)%tp%ident%basename=='OC') then
            i=i+1
            CALL get_tracer(status, GPTRSTR, trim(ti_gp(jt)%tp%ident%basename),&
                            subname=trim(ti_gp(jt)%tp%ident%subname),          &
                            idx=om_idt(i) )
        ! fac_om is multiplied by 1.3 to take into account not only the
        ! carbon but other atoms according to the original work of 
        ! Menon et al. JAS, 2002
            fac_om = ti_gp(jt)%tp%meta%cask_r(R_molarmass) / M_air
          ENDIF
          IF (MATCH_WILD( 'WSOC*', ti_gp(jt)%tp%ident%basename) ) THEN
            i = i+1
            CALL get_tracer(status, GPTRSTR, trim(ti_gp(jt)%tp%ident%basename),&
                            subname=trim(ti_gp(jt)%tp%ident%subname),          &
                            idx=om_idt(i) )
!!$        print*, om_idt, ti_gp(jt)%tp%ident%fullname, fac_om,  &
!!$                ti_gp(jt)%tp%meta%cask_r(R_molarmass), M_air
          endif
        enddo
      else
        if (i_cdnc_calc == 4) then
          CALL INFO_BI( &
            'No organic aerosol found, but required for Menon CDNC calculation!', substr)
          CALL INFO_BI( &
            'Calculation may be possible, but is likely to be wrong (incomplete)!', substr)
        endif
      endif

      IF ( (i_cdnc_calc == 5) .or. &
           (i_cdnc_calc == 6) .or. &
           (i_cdnc_calc == 7) .or. &
           (i_cdnc_calc == 9) .or. &
           (ncdnc       >  0) ) THEN
         
         CALL get_channel_object(status, TRIM(aer_stream), 'sigma', p1=sigma)
         IF (status /= 0) &
              call error_bi('channel object for sigma not found', substr)
         IF (status == 0) THEN
            ALLOCATE(cmr2mmr(SIZE(sigma)))
            DO jt=1,SIZE(sigma)
               cmr2mmr(jt) = EXP(3.5*(LOG(sigma(jt)))**2)
            END DO
         ENDIF
         
         CALL get_channel_object(status, TRIM(aer_stream),&
              'wetradius', p4=wetradius)
         IF (status /= 0) &
              CALL info_bi('wet radius channel object not found !', substr)
         CALL get_channel_object(status, TRIM(aer_stream),&
              'dryradius', p4=dryradius)
         IF (status /= 0) &
              CALL info_bi('dry radius channel object not found !', substr)
         
! Explicit assignment of ksol and of the indexes for freezing schemes
         SELECT CASE (TRIM(aer_stream))
         CASE ('made3_gp')
            ksol = 6
            ktot = 6                  !mz_vk_20170703: needed because of the changes done for the gmxe case to use UAF
            ALLOCATE(sol_name(ksol))
            ALLOCATE(sol_idx(ksol))
            sol_name(:) = (/"ks", "km", "as", "am", "cs", "cm"/)
            sol_idx(:) = (/1, 2, 4, 5, 7, 8/)
            inucs = 1 ! not avail., set to ks (not used anyway)
            iaits = 1
            iaitm = 2
            iaiti = 3
            iaccs = 4
            iaccm = 5
            iacci = 6
            icoas = 7
            icoam = 8
            icoai = 9
         CASE ('gmxe_gp')
            ksol = 4

            !Since UAF considers 7 modes while the variable ksol
            !includs just the soluble modes, we differetiate here the two cases
            IF ( i_cdnc_calc == 7 ) THEN  !Only if UAF is used, ktot =7
               ktot = 7 
               ALLOCATE(sol_name(ktot))
               ALLOCATE(sol_idx(ktot))
               sol_name(:) = (/"ns", "ks", "as", "cs", "ki", "ai", "ci"/)
               sol_idx(:) = (/1, 2, 3, 4, 5, 6, 7/)
            ELSE
               ktot = ksol
               ALLOCATE(sol_name(ksol))
               ALLOCATE(sol_idx(ksol))
               sol_name(:) = (/"ns", "ks", "as", "cs"/)
               sol_idx(:) = (/1, 2, 3, 4/)
            END IF

            inucs = 1
            iaits = 2
            iaccs = 3
            icoas = 4
            iaitm = 2
            iaccm = 3
            icoam = 4
            iaiti = 5
            iacci = 6
            icoai = 7
         CASE DEFAULT
            CALL error_bi('coupling for ' // aer_stream // ' not available'&
                 , substr)
         END SELECT
      
         ! create new channel objects for activated fraction(s) to be used by SCAV for nucleation scavenging
         CALL new_dimension(status, DIMID_NSOLMODE, 'CLOUD_SOLMODE', ksol)
         CALL channel_halt(substr, status)
         
         ! NEW REPRESENTATIONS
         CALL new_representation(status, REPR_CLOUD_4D_KSOL, &
              'REPR_CLOUD_4D_KSOL'    &
              , rank = 4, link = 'xxxx', dctype = DC_GP               &
              , dimension_ids = (/ DIMID_LON, DIMID_LEV               &
              ,                    DIMID_NSOLMODE, DIMID_LAT /)       &
              , ldimlen       = (/ nproma, AUTO, AUTO, ngpblks   /)   &
              , output_order  = (/ 3,1,4,2 /)                         &
              , axis = 'XZNY'                                         &
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
         
         cnt(:,3) = ksol
         memu(:,3) = ksol
         
         CALL set_representation_decomp(status, REPR_CLOUD_4D_KSOL &
              , start, cnt, memu, meml, .FALSE., PIOTYPE_COL)
         CALL channel_halt(substr, status)
         
         DEALLOCATE(start) ; NULLIFY(start)
         DEALLOCATE(cnt)   ; NULLIFY(cnt)
         DEALLOCATE(meml)  ; NULLIFY(meml)
         DEALLOCATE(memu)  ; NULLIFY(memu)
         
         IF ( (i_cdnc_calc == 5) .or. &
              (i_cdnc_calc == 9) .or. &
              (ncdnc       == 2) ) THEN
            ALLOCATE(p_arg_frac(nproma,nlev,ksol,ngpblks,1))
            mem => p_arg_frac(:,:,:,:,1)
            CALL new_channel_object(status, modstr, 'ARG_ACT_FRAC' &
                 , p4 = arg_frac, reprid=REPR_CLOUD_4D_KSOL, mem=mem)
            CALL channel_halt(substr, status)
            CALL new_attribute(status, modstr, 'ARG_ACT_FRAC'  &
                 , 'long_name', c='activated number fraction (ARG scheme)')
            CALL channel_halt(substr, status)
            CALL new_attribute(status, modstr, 'ARG_ACT_FRAC'   &
                 , 'units', c='-')
            CALL channel_halt(substr, status)
            CALL new_attribute(status, modstr, 'ARG_ACT_FRAC'   &
                 , 'aermod', c=TRIM(aer_stream) )
            CALL channel_halt(substr, status)

            ALLOCATE(p_arg_frac2(nproma,nlev,ksol,ngpblks,1))
            mem => p_arg_frac2(:,:,:,:,1)
            CALL new_channel_object(status, modstr, 'ARG_ACT_FRAC2' &
                 , p4 = arg_frac2, reprid=REPR_CLOUD_4D_KSOL, mem=mem)
            CALL channel_halt(substr, status)
            CALL new_attribute(status, modstr, 'ARG_ACT_FRAC2'  &
                 , 'long_name', c='activated number fraction2 (ARG scheme)')
            CALL channel_halt(substr, status)
            CALL new_attribute(status, modstr, 'ARG_ACT_FRAC2'   &
                 , 'units', c='-')
            CALL channel_halt(substr, status)
            CALL new_attribute(status, modstr, 'ARG_ACT_FRAC2'   &
                 , 'aermod', c=TRIM(aer_stream) )
            CALL channel_halt(substr, status)
            
            ALLOCATE(p_arg_mfrac(nproma,nlev,ksol,ngpblks,1))
            mem => p_arg_mfrac(:,:,:,:,1)
            CALL new_channel_object(status, modstr, 'ARG_ACT_MFRAC' &
                 , p4 = arg_mfrac, reprid=REPR_CLOUD_4D_KSOL, mem=mem)
            CALL channel_halt(substr, status)
            CALL new_attribute(status, modstr, 'ARG_ACT_MFRAC'  &
                 , 'long_name', c='activated mass fraction (ARG scheme)')
            CALL channel_halt(substr, status)
            CALL new_attribute(status, modstr, 'ARG_ACT_MFRAC'   &
                 , 'units', c='-')
            CALL channel_halt(substr, status)
            CALL new_attribute(status, modstr, 'ARG_ACT_MFRAC'   &
                 , 'aermod', c=TRIM(aer_stream) )
            CALL channel_halt(substr, status)
            
            ALLOCATE(p_arg_num(nproma,nlev,ksol,ngpblks,1))
            mem => p_arg_num(:,:,:,:,1)
            CALL new_channel_object(status, modstr, 'ARG_ACT_NUM' &
                 , p4 = arg_num, reprid=REPR_CLOUD_4D_KSOL, mem=mem)
            CALL channel_halt(substr, status)
            CALL new_attribute(status, modstr, 'ARG_ACT_NUM'  &
                , 'long_name', c='activated number (ARG scheme)')
            CALL channel_halt(substr, status)
            CALL new_attribute(status, modstr, 'ARG_ACT_NUM'   &
                 , 'units', c='1/m^3')
            CALL channel_halt(substr, status)

! New channel output for critical radius
            ALLOCATE(p_arg_rad(nproma,nlev,ksol,ngpblks,1))
            mem => p_arg_rad(:,:,:,:,1)
            CALL new_channel_object(status, modstr, 'ARG_ACT_RAD' &
                 , p4 = arg_rad, reprid=REPR_CLOUD_4D_KSOL, mem=mem)
            CALL channel_halt(substr, status)
            CALL new_attribute(status, modstr, 'ARG_ACT_RAD'  &
                , 'long_name', c='critical radius (ARG scheme)')
            CALL channel_halt(substr, status)
            CALL new_attribute(status, modstr, 'ARG_ACT_RAD'   &
                 , 'units', c='m')
            CALL channel_halt(substr, status)

            ALLOCATE(p_arg_scrit(nproma,nlev,ksol,ngpblks,1))
            mem => p_arg_scrit(:,:,:,:,1)
            CALL new_channel_object(status, modstr, 'ARG_S_CRIT' &
                 , p4 = arg_scrit, reprid=REPR_CLOUD_4D_KSOL, mem=mem)
            CALL channel_halt(substr, status)
            CALL new_attribute(status, modstr, 'ARG_S_CRIT'  &
                , 'long_name', c='critical supersaturation (ARG scheme)')
            CALL channel_halt(substr, status)
            CALL new_attribute(status, modstr, 'ARG_S_CRIT'   &
                 , 'units', c='1')
            CALL channel_halt(substr, status)

            DO jm=1,ksol
               CALL int2str(char1, jm)
               ! activated fraction per mode
               WRITE(name,'(A13,A1)') 'ARG_ACT_FRAC_',char1
               mem => p_arg_frac(:,:,jm,:,:)
               CALL new_channel_object(status, modstr, TRIM(name), mem=mem)
               CALL channel_halt(substr, status)
               CALL new_attribute(status, modstr, TRIM(name) &
                    , 'long_name', &
                    c='activated number fraction in Mode '//char1//' (ARG scheme)' )
               CALL channel_halt(substr, status)
               CALL new_attribute(status, modstr, TRIM(name) &
                    , 'units', c='-' )
               CALL channel_halt(substr, status)

               WRITE(name,'(A14,A1)') 'ARG_ACT_FRAC2_',char1
               mem => p_arg_frac2(:,:,jm,:,:)
               CALL new_channel_object(status, modstr, TRIM(name), mem=mem)
               CALL channel_halt(substr, status)
               CALL new_attribute(status, modstr, TRIM(name) &
                    , 'long_name', &
                    c='activated number2 fraction in Mode '//char1//' (ARG scheme)' )
               CALL channel_halt(substr, status)
               CALL new_attribute(status, modstr, TRIM(name) &
                    , 'units', c='-' )
               CALL channel_halt(substr, status)
               
               WRITE(name,'(A14,A1)') 'ARG_ACT_MFRAC_',char1
               mem => p_arg_mfrac(:,:,jm,:,:)
               CALL new_channel_object(status, modstr, TRIM(name), mem=mem)
               CALL channel_halt(substr, status)
               CALL new_attribute(status, modstr, TRIM(name) &
                    , 'long_name', &
                    c='activated mass fraction in Mode '//char1//' (ARG scheme)' )
               CALL channel_halt(substr, status)
               CALL new_attribute(status, modstr, TRIM(name) &
                    , 'units', c='-' )
               CALL channel_halt(substr, status)
               
               WRITE(name,'(A12,A1)') 'ARG_ACT_NUM_',char1
               mem => p_arg_num(:,:,jm,:,:)
               CALL new_channel_object(status, modstr, TRIM(name), mem=mem)
               CALL channel_halt(substr, status)
               CALL new_attribute(status, modstr, TRIM(name) &
                    , 'long_name', &
                    c='activated numbers in Mode '//char1//' (ARG scheme)' )
               CALL channel_halt(substr, status)
               CALL new_attribute(status, modstr, TRIM(name) &
                    , 'units', c='1/m^3' )
               CALL channel_halt(substr, status)

               WRITE(name,'(A12,A1)') 'ARG_ACT_RAD_',char1
               mem => p_arg_rad(:,:,jm,:,:)
               CALL new_channel_object(status, modstr, TRIM(name), mem=mem)
               CALL channel_halt(substr, status)
               CALL new_attribute(status, modstr, TRIM(name) &
                    , 'long_name', &
                    c='critical radius in Mode '//char1//' (ARG scheme)' )
               CALL channel_halt(substr, status)
               CALL new_attribute(status, modstr, TRIM(name), 'units', c='m')
               CALL channel_halt(substr, status)

               WRITE(name,'(A11,A1)') 'ARG_S_CRIT_',char1
               mem => p_arg_scrit(:,:,jm,:,:)
               CALL new_channel_object(status, modstr, TRIM(name), mem=mem)
               CALL channel_halt(substr, status)
               CALL new_attribute(status, modstr, TRIM(name) &
                    , 'long_name', &
                    c='critical supersaturation in Mode '//char1//' (ARG scheme)' )
               CALL channel_halt(substr, status)
               CALL new_attribute(status, modstr, TRIM(name), 'units', c='1')
               CALL channel_halt(substr, status)
            END DO
            p_arg_frac(:,:,:,:,:)  = 0._dp
            p_arg_frac2(:,:,:,:,:) = 0._dp
            p_arg_mfrac(:,:,:,:,:) = 0._dp
            p_arg_num(:,:,:,:,:)   = 0._dp
            p_arg_rad(:,:,:,:,:)   = 0._dp
            p_arg_scrit(:,:,:,:,:) = 0._dp
             
            CALL new_channel_object(status, modstr, 'Scrit', p3=scrit)
            CALL channel_halt(substr, status)
            CALL new_attribute(status, modstr, 'Scrit' &
                 , 'long_name', c='maximum supersaturation')
            CALL channel_halt(substr, status)
            CALL new_attribute(status, modstr, 'Scrit', &
                 'units', c='%' )
            CALL channel_halt(substr, status)
         END IF

         IF ( (i_cdnc_calc == 6) .or. &
              (i_cdnc_calc == 9) .or. &
              (ncdnc       == 1) ) THEN
            
            ALLOCATE(p_lin_frac(nproma,nlev,ksol,ngpblks,1))
            mem => p_lin_frac(:,:,:,:,1)
            CALL new_channel_object(status, modstr, 'LIN_ACT_FRAC' &
                 , p4 = lin_frac, reprid=REPR_CLOUD_4D_KSOL, mem=mem)
            CALL channel_halt(substr, status)
            CALL new_attribute(status, modstr, 'LIN_ACT_FRAC'  &
                 , 'long_name', c='activated number fraction (LIN scheme)')
            CALL channel_halt(substr, status)
            CALL new_attribute(status, modstr, 'LIN_ACT_FRAC'   &
                 , 'units', c='-')
            CALL channel_halt(substr, status)
            CALL new_attribute(status, modstr, 'LIN_ACT_FRAC'   &
                 , 'aermod', c=TRIM(aer_stream) )
            CALL channel_halt(substr, status)
            
            ALLOCATE(p_lin_mfrac(nproma,nlev,ksol,ngpblks,1))
            mem => p_lin_mfrac(:,:,:,:,1)
            CALL new_channel_object(status, modstr, 'LIN_ACT_MFRAC' &
                 , p4 = lin_mfrac, reprid=REPR_CLOUD_4D_KSOL, mem=mem)
            CALL channel_halt(substr, status)
            CALL new_attribute(status, modstr, 'LIN_ACT_MFRAC'  &
                 , 'long_name', c='activated mass fraction (LIN scheme)')
            CALL channel_halt(substr, status)
            CALL new_attribute(status, modstr, 'LIN_ACT_MFRAC'   &
                 , 'units', c='-')
            CALL channel_halt(substr, status)
            CALL new_attribute(status, modstr, 'LIN_ACT_MFRAC'   &
                 , 'aermod', c=TRIM(aer_stream) )
            CALL channel_halt(substr, status)
            
            ALLOCATE(p_lin_num(nproma,nlev,ksol,ngpblks,1))
            mem => p_lin_num(:,:,:,:,1)
            CALL new_channel_object(status, modstr, 'LIN_ACT_NUM' &
                 , p4 = lin_num, reprid=REPR_CLOUD_4D_KSOL, mem=mem)
            CALL channel_halt(substr, status)
            CALL new_attribute(status, modstr, 'LIN_ACT_NUM'  &
                 , 'long_name', c='activated number (LIN scheme)')
            CALL channel_halt(substr, status)
            CALL new_attribute(status, modstr, 'LIN_ACT_NUM'   &
                 , 'units', c='1/cm^3')
            CALL channel_halt(substr, status)
            DO jm=1,ksol
               CALL int2str(char1, jm)
               ! activated fraction per mode
               WRITE(name,'(A13,A1)') 'LIN_ACT_FRAC_',char1
               mem => p_lin_frac(:,:,jm,:,:)
               CALL new_channel_object(status, modstr, TRIM(name), mem=mem)
               CALL channel_halt(substr, status)
               CALL new_attribute(status, modstr, TRIM(name) &
                    , 'long_name', &
                    c='activated number fraction in Mode '//char1//' (LIN scheme)' )
               CALL channel_halt(substr, status)
               CALL new_attribute(status, modstr, TRIM(name) &
                    , 'units', c='-' )
               CALL channel_halt(substr, status)
               
               WRITE(name,'(A14,A1)') 'LIN_ACT_MFRAC_',char1
               mem => p_lin_mfrac(:,:,jm,:,:)
               CALL new_channel_object(status, modstr, TRIM(name), mem=mem)
               CALL channel_halt(substr, status)
               CALL new_attribute(status, modstr, TRIM(name) &
                    , 'long_name', &
                    c='activated mass fraction in Mode '//char1//' (LIN scheme)' )
               CALL channel_halt(substr, status)
               CALL new_attribute(status, modstr, TRIM(name) &
                    , 'units', c='-' )
               CALL channel_halt(substr, status)
               
               WRITE(name,'(A12,A1)') 'LIN_ACT_NUM_',char1
               mem => p_lin_num(:,:,jm,:,:)
               CALL new_channel_object(status, modstr, TRIM(name), mem=mem)
               CALL channel_halt(substr, status)
               CALL new_attribute(status, modstr, TRIM(name) &
                    , 'long_name', &
                    c='activated numbers in Mode '//char1//' (LIN scheme)' )
               CALL channel_halt(substr, status)
               CALL new_attribute(status, modstr, TRIM(name) &
                    , 'units', c='1/cm^3' )
               CALL channel_halt(substr, status)
            END DO
            p_lin_frac(:,:,:,:,:)  = 0._dp
            p_lin_mfrac(:,:,:,:,:) = 0._dp
            p_lin_num(:,:,:,:,:)   = 0._dp
         END IF
          
!For the UAF parameterization 

         IF ( (i_cdnc_calc == 7) .or. &
              (i_cdnc_calc == 9) .or. &
              (ncdnc       == 3) ) THEN 

            ! Get the effective hygroscopicity of each soluble mode for the UAF Parameterization --- VAK
            ALLOCATE (kappa_uaf(ksol))
            DO ji=1,ksol
               write (str_num,'(I1)') ji
               CALL get_channel_object(status, aer_stream, &
                    'DIAGAER_M'//TRIM(ADJUSTL(str_num))//'_25',&
                    p3=kappa_uaf(ji)%ptr )
               IF (status /= 0) THEN 
                  call error_bi('channel object kappa value not found', substr)
               ENDIF
            END DO

         ! create new channel objects for activated fraction(s) to be used by SCAV for nucleation scavenging --- VAK
         CALL new_dimension(status, DIMID_NTOTMODE, 'CLOUD_TOTMODE', ktot)
         CALL channel_halt(substr, status)
     
         ! NEW REPRESENTATIONS
         CALL new_representation(status, REPR_CLOUD_4D_KTOT, &
              'REPR_CLOUD_4D_KTOT'    &    
              , rank = 4, link = 'xxxx', dctype = DC_GP               &    
              , dimension_ids = (/ DIMID_LON, DIMID_LEV               &    
              ,                    DIMID_NTOTMODE, DIMID_LAT /)       &    
              , ldimlen       = (/ nproma, AUTO, AUTO, ngpblks   /)   &
              , output_order  = (/ 3,1,4,2 /)                         &    
              , axis = 'XZNY'                                         &    
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

         cnt(:,3) = ktot
         memu(:,3) = ktot

         CALL set_representation_decomp(status, REPR_CLOUD_4D_KTOT &
              , start, cnt, memu, meml, .FALSE., PIOTYPE_COL)
         CALL channel_halt(substr, status)

         DEALLOCATE(start) ; NULLIFY(start)
         DEALLOCATE(cnt)   ; NULLIFY(cnt)
         DEALLOCATE(meml)  ; NULLIFY(meml)
         DEALLOCATE(memu)  ; NULLIFY(memu)

            CALL new_channel_object(status, modstr, 'Scrit', p3=scrit)
            CALL channel_halt(substr, status)
            CALL new_attribute(status, modstr, 'Scrit' &
                 , 'long_name', c='critical supersaturation')
            CALL channel_halt(substr, status)
            CALL new_attribute(status, modstr, 'Scrit', &
                 'units', c='%' )
            CALL channel_halt(substr, status)


            CALL new_channel_object(status, modstr, 'WPARC', p3=WPARC)
            CALL channel_halt(substr, status)
            CALL new_attribute(status, modstr, 'WPARC' &
                 , 'long_name', c='updraft velocity')
            CALL channel_halt(substr, status)
            CALL new_attribute(status, modstr, 'WPARC', &
                 'units', c='m s-1' )
            CALL channel_halt(substr, status)

            ALLOCATE(p_uaf_num(nproma,nlev,ktot,ngpblks,1))
            mem => p_uaf_num(:,:,:,:,1)
            CALL new_channel_object(status, modstr, 'UAF_ACT_NUM' &
                 , p4 = uaf_num, reprid=REPR_CLOUD_4D_KTOT, mem=mem)
            CALL channel_halt(substr, status)
            CALL new_attribute(status, modstr, 'UAF_ACT_NUM'  &
                 , 'long_name', c='activated number (UAF scheme)')
            CALL channel_halt(substr, status)
            CALL new_attribute(status, modstr, 'UAF_ACT_NUM'   &
                 , 'units', c='1/m^3')
            CALL channel_halt(substr, status)

            ALLOCATE(p_uaf_frac(nproma,nlev,ktot,ngpblks,1))
            mem => p_uaf_frac(:,:,:,:,1)
            CALL new_channel_object(status, modstr, 'UAF_ACT_FRAC' &
                 , p4 = uaf_frac, reprid=REPR_CLOUD_4D_KTOT, mem=mem)
            CALL channel_halt(substr, status)
            CALL new_attribute(status, modstr, 'UAF_ACT_FRAC'  &
                 , 'long_name', c='activated number fraction (UAF scheme)')
            CALL channel_halt(substr, status)
            CALL new_attribute(status, modstr, 'UAF_ACT_FRAC'   &
                 , 'units', c='-')
            CALL channel_halt(substr, status)
            CALL new_attribute(status, modstr, 'UAF_ACT_FRAC'   &
                 , 'aermod', c=TRIM(aer_stream) )
            CALL channel_halt(substr, status)

            ALLOCATE(p_uaf_mfrac(nproma,nlev,ktot,ngpblks,1))
            mem => p_uaf_mfrac(:,:,:,:,1)
            CALL new_channel_object(status, modstr, 'UAF_ACT_MFRAC' &
                 , p4 = uaf_mfrac, reprid=REPR_CLOUD_4D_KTOT, mem=mem)
            CALL channel_halt(substr, status)
            CALL new_attribute(status, modstr, 'UAF_ACT_MFRAC'  &
                 , 'long_name', c='activated mass fraction (UAF scheme)')
            CALL channel_halt(substr, status)
            CALL new_attribute(status, modstr, 'UAF_ACT_MFRAC'   &
                 , 'units', c='-')
            CALL new_attribute(status, modstr, 'UAF_ACT_MFRAC'   &
                 , 'aermod', c=TRIM(aer_stream) )
            CALL channel_halt(substr, status)

            ALLOCATE(p_uaf_AKK(nproma,nlev,ktot,ngpblks,1))
            mem => p_uaf_AKK(:,:,:,:,1)
            CALL new_channel_object(status, modstr, 'UAF_SOL_AKK' &
                 , p4 = uaf_AKK, reprid=REPR_CLOUD_4D_KTOT, mem=mem)
            CALL channel_halt(substr, status)
            CALL new_attribute(status, modstr, 'UAF_SOL_AKK'  &
                 , 'long_name', c='hygroscopivity of soluble fraction (UAF scheme)')
            CALL channel_halt(substr, status)
            CALL new_attribute(status, modstr, 'UAF_SOL_AKK'   &
                 , 'units', c='-')
            CALL channel_halt(substr, status)

            ALLOCATE(p_uaf_TPI(nproma,nlev,ktot,ngpblks,1))
            mem => p_uaf_TPI(:,:,:,:,1)
            CALL new_channel_object(status, modstr, 'UAF_MOD_TPI' &
                 , p4 = uaf_TPI, reprid=REPR_CLOUD_4D_KTOT, mem=mem)
            CALL channel_halt(substr, status)
            CALL new_attribute(status, modstr, 'UAF_MOD_TPI'  &
                 , 'long_name', c='aerosol number (UAF scheme)')
            CALL channel_halt(substr, status)
            CALL new_attribute(status, modstr, 'UAF_MOD_TPI'   &
                 , 'units', c='1/m^3')
            CALL channel_halt(substr, status)

            ALLOCATE(p_uaf_DPG(nproma,nlev,ktot,ngpblks,1))
            mem => p_uaf_DPG(:,:,:,:,1)
            CALL new_channel_object(status, modstr, 'UAF_MOD_DPG' &
                 , p4 = uaf_DPG, reprid=REPR_CLOUD_4D_KTOT, mem=mem)
            CALL channel_halt(substr, status)
            CALL new_attribute(status, modstr, 'UAF_MOD_DPG'  &
                 , 'long_name', c='dry diameterty(UAF scheme)')
            CALL channel_halt(substr, status)
            CALL new_attribute(status, modstr, 'UAF_MOD_DPG'   &
                 , 'units', c='m')
            CALL channel_halt(substr, status)

            ALLOCATE(p_uaf_SIG(nproma,nlev,ktot,ngpblks,1))
            mem => p_uaf_SIG(:,:,:,:,1)
            CALL new_channel_object(status, modstr, 'UAF_MOD_SIG' &
                 , p4 = uaf_SIG, reprid=REPR_CLOUD_4D_KTOT, mem=mem)
            CALL channel_halt(substr, status)
            CALL new_attribute(status, modstr, 'UAF_MOD_SIG'  &
                 , 'long_name', c='sigma (UAF scheme)')
            CALL channel_halt(substr, status)
            CALL new_attribute(status, modstr, 'UAF_MOD_SIG'   &
                 , 'units', c='-')
            CALL channel_halt(substr, status)

            ALLOCATE(p_uaf_ei(nproma,nlev,ktot,ngpblks,1))
            mem => p_uaf_ei(:,:,:,:,1)
            CALL new_channel_object(status, modstr, 'UAF_INSOL_FRAC' &
                 , p4 = uaf_ei, reprid=REPR_CLOUD_4D_KTOT, mem=mem)
            CALL channel_halt(substr, status)
            CALL new_attribute(status, modstr, 'UAF_INSOL_FRAC'  &
                 , 'long_name', c='insoluble fraction (UAF scheme)')
            CALL channel_halt(substr, status)
            CALL new_attribute(status, modstr, 'UAF_INSOL_FRAC'   &
                 , 'units', c='-')
            CALL channel_halt(substr, status)


            DO jm=1,ktot
               CALL int2str(char1, jm)
               ! activated fraction per mode
               WRITE(name,'(A12,A1)') 'UAF_ACT_NUM_',char1
               mem => p_uaf_num(:,:,jm,:,:)
               CALL new_channel_object(status, modstr, TRIM(name), mem=mem)
               CALL channel_halt(substr, status)
               CALL new_attribute(status, modstr, TRIM(name) &
                    , 'long_name', &
                    c='activated number in Mode '//char1//' (UAF scheme)' )
               CALL channel_halt(substr, status)
               CALL new_attribute(status, modstr, TRIM(name) &
                    , 'units', c='1/m^3' )
               CALL channel_halt(substr, status)

               ! activated fraction per mode
               WRITE(name,'(A13,A1)') 'UAF_ACT_FRAC_',char1
               mem => p_uaf_frac(:,:,jm,:,:)
               CALL new_channel_object(status, modstr, TRIM(name), mem=mem)
               CALL channel_halt(substr, status)
               CALL new_attribute(status, modstr, TRIM(name) &
                    , 'long_name', &
                    c='activated number fraction in Mode '//char1//' (UAF scheme)' )
               CALL channel_halt(substr, status)
               CALL new_attribute(status, modstr, TRIM(name) &
                    , 'units', c='-' )
               CALL channel_halt(substr, status)

               ! activated mass fraction per mode
               WRITE(name,'(A14,A1)') 'UAF_ACT_MFRAC_',char1
               mem => p_uaf_mfrac(:,:,jm,:,:)
               CALL new_channel_object(status, modstr, TRIM(name), mem=mem)
               CALL channel_halt(substr, status)
               CALL new_attribute(status, modstr, TRIM(name) &
                    , 'long_name', &
                    c='activated mass fraction in Mode '//char1//' (UAF scheme)' )
               CALL channel_halt(substr, status)
               CALL new_attribute(status, modstr, TRIM(name) &
                    , 'units', c='-' )
               CALL channel_halt(substr, status)

               ! hygroscopicity per mode
               WRITE(name,'(A12,A1)') 'UAF_SOL_AKK_',char1
               mem => p_uaf_AKK(:,:,jm,:,:)
               CALL new_channel_object(status, modstr, TRIM(name), mem=mem)
               CALL channel_halt(substr, status)
               CALL new_attribute(status, modstr, TRIM(name) &
                    , 'long_name', &
                    c='hygroscopicity of soluble fraction in Mode '//char1//' (UAF scheme)' )
               CALL channel_halt(substr, status)
               CALL new_attribute(status, modstr, TRIM(name) &
                    , 'units', c='-' )
               CALL channel_halt(substr, status)

              ! aerosol number per mode
               WRITE(name,'(A12,A1)') 'UAF_MOD_TPI_',char1
               mem => p_uaf_TPI(:,:,jm,:,:)
               CALL new_channel_object(status, modstr, TRIM(name), mem=mem)
               CALL channel_halt(substr, status)
               CALL new_attribute(status, modstr, TRIM(name) &
                    , 'long_name', &
                    c='Aerosol number in Mode '//char1//' (UAF scheme)' )
               CALL channel_halt(substr, status)
               CALL new_attribute(status, modstr, TRIM(name) &
                    , 'units', c='1/m^3' )
               CALL channel_halt(substr, status)

               ! dry diameter per mode
               WRITE(name,'(A12,A1)') 'UAF_MOD_DPG_',char1
               mem => p_uaf_DPG(:,:,jm,:,:)
               CALL new_channel_object(status, modstr, TRIM(name), mem=mem)
               CALL channel_halt(substr, status)
               CALL new_attribute(status, modstr, TRIM(name) &
                    , 'long_name', &
                    c='Dry diameter of Mode '//char1//' (UAF scheme)' )
               CALL channel_halt(substr, status)
               CALL new_attribute(status, modstr, TRIM(name) &
                    , 'units', c='m' )
               CALL channel_halt(substr, status)

               ! sigma per mode
               WRITE(name,'(A12,A1)') 'UAF_MOD_SIG_',char1
               mem => p_uaf_SIG(:,:,jm,:,:)
               CALL new_channel_object(status, modstr, TRIM(name), mem=mem)
               CALL channel_halt(substr, status)
               CALL new_attribute(status, modstr, TRIM(name) &
                    , 'long_name', &
                    c='Sigma of Mode '//char1//' (UAF scheme)' )
               CALL channel_halt(substr, status)
               CALL new_attribute(status, modstr, TRIM(name) &
                    , 'units', c='-' )
               CALL channel_halt(substr, status)

               ! activated fraction per mode
               WRITE(name,'(A14,A1)') 'UAF_INSOL_FRC_',char1
               mem => p_uaf_ei(:,:,jm,:,:)
               CALL new_channel_object(status, modstr, TRIM(name), mem=mem)
               CALL channel_halt(substr, status)
               CALL new_attribute(status, modstr, TRIM(name) &
                    , 'long_name', &
                    c='insoluble fraction in Mode '//char1//' (UAF scheme)' )
               CALL channel_halt(substr, status)
               CALL new_attribute(status, modstr, TRIM(name) &
                    , 'units', c='-' )
               CALL channel_halt(substr, status)
            END DO
            p_uaf_frac(:,:,:,:,:)   = 0._dp
            p_uaf_mfrac(:,:,:,:,:)   = 0._dp
            p_uaf_num(:,:,:,:,:)   = 0._dp
            p_uaf_AKK(:,:,:,:,:)   = 0._dp
            p_uaf_TPI(:,:,:,:,:)   = 0._dp
            p_uaf_DPG(:,:,:,:,:)   = 0._dp
            p_uaf_SIG(:,:,:,:,:)   = 0._dp
            p_uaf_ei(:,:,:,:,:)   = 0._dp

         END IF     !End case UAF

         ALLOCATE(aerosol_num(ktot))
         ALLOCATE(num_idx(ktot))

         num_idx(:) = 0

         SELECT CASE (TRIM(aer_stream))
         CASE ('made3_gp')
            DO jm = 1,ksol
               DO jt=1,ntrac
                  IF ((laerosol(jt)) .AND. &
                       (ti_gp(jt)%tp%ident%quantity == NUMBERDENSITY)) THEN
                     IF (ti_gp(jt)%tp%ident%subname == sol_name(jm)) THEN
                        num_idx(jm) = jt
                     ENDIF
                  ENDIF
               ENDDO
            ENDDO
         CASE ('gmxe_gp')
            DO jm = 1,ktot
               DO jt=1,ntrac
                  IF ((laerosol(jt)) .AND. &
                       (ti_gp(jt)%tp%ident%quantity == NUMBERDENSITY)) THEN
                     IF (ti_gp(jt)%tp%ident%subname == sol_name(jm)) THEN
                        num_idx(jm) = jt
                     ENDIF
                  ENDIF
               ENDDO
            ENDDO
            ! No aerosol number tracers, but stream elements
            do jm = 1,ktot
               if (num_idx(jm) == 0) THEN
                  no_trac=.TRUE.
                  write(mod_number,"(i2)") jm
                  CALL INFO_BI(&
                      'Mode '//mod_number//': No tracers with number concentration found for this mode', substr)
                  CALL INFO_BI('Looking for stream elements', substr)
               endif
            enddo
            IF (no_trac) THEN
               CALL get_channel_object(status, TRIM(aer_stream), 'anumber', &
                    p4=anumber)
               IF (status /= 0) &
                    CALL info_bi('aerosol number stream element not found !' &
                    , substr)
            ENDIF

         CASE DEFAULT
            CALL error_bi('coupling for ' // aer_stream // ' not available' &
                 , substr)
         END SELECT


         allocate(base_name(ntrac))
         nspec      = 0
         base_name(:) = ""

         DO jt = 1,ntrac
            IF ( (laerosol(jt)) .AND.                                  &
                 (ti_gp(jt)%tp%ident%quantity == AMOUNTFRACTION) ) THEN

! Aerosol water and Hpres shall be excluded from the total mass calculation
! (used in the denominator of Eq. 4 in Abdul-Razzak and Ghan (2000))
               IF ( ANY(base_name(:) == TRIM(ti_gp(jt)%tp%ident%basename)).OR. &
                    "H2O" == TRIM(ti_gp(jt)%tp%ident%basename).OR. &
                    "Hpres" == TRIM(ti_gp(jt)%tp%ident%basename)) THEN
                  CYCLE
               ELSE
                  nspec = nspec + 1 
                  base_name(nspec) = TRIM(ti_gp(jt)%tp%ident%basename)
                  mmass(nspec)     = ti_gp(jt)%tp%meta%cask_r(R_molarmass)
                  aerdens(nspec)   = ti_gp(jt)%tp%meta%cask_r(R_aerosol_density)
               ENDIF
            ENDIF
         ENDDO
         ALLOCATE(mass(ktot,nspec))
         ALLOCATE(aer_idx(ktot,nspec))
         ALLOCATE(nu(nspec))
         ALLOCATE(eps(nspec))
         ALLOCATE(phi(nspec))
         ALLOCATE(kappa(nspec))
         ALLOCATE(molar(nspec))
         ALLOCATE(aerdensity(nspec))
         ALLOCATE(FHH(nspec))

! Added MADE and MADE3 cases
         aerdensity(:) = 0._dp
         aer_idx(:,:)  = 0 
         molar(1:nspec)      = mmass(1:nspec)
         DO jt=1,nspec
            IF (aerdens(jt) >= 1._dp) THEN
               aerdensity(jt) = aerdens(jt)
            ELSE
               aerdensity(jt) = 1000._dp
            ENDIF
         ENDDO

         SELECT CASE (TRIM(aer_stream))
         CASE ('made3_gp')
            NU(:) = 1.
            EPS(:) = 1._dp
            PHI(:) = 1._dp

            DO jt = 1,nspec
               SELECT CASE (TRIM(base_name(JT)))
               CASE('BC', 'BCtag', 'POM', 'DU')
                  NU(jt)  = 0._dp
                  EPS(jt) = 0._dp
                  PHI(jt) = 0._dp
                  KAPPA(jt) = 0._dp
               CASE('SS')
                  NU(jt)  = 2._dp
                  KAPPA(jt) = 1.12_dp
               CASE('Na', 'Cl')
                  KAPPA(jt) = 1.12_dp
               CASE('SO4')
                  KAPPA(jt) = 1.19_dp
               CASE('NO3', 'NH4')
                  KAPPA(jt) = 0.67_dp
               END SELECT
            ENDDO
         CASE ('gmxe_gp')
            NU(:)         = 1._dp
            EPS(:)        = 1._dp
            PHI(:)        = 0.7_dp
            FHH(:)        = .false.

            DO jt = 1,nspec
               SELECT CASE (TRIM(base_name(JT)))
               CASE('SS')
                  NU(jt)  = 2._dp
               CASE('OC')
                  EPS(jt) = 0.5_dp
                  PHI(jt) = 0.2_dp
               CASE('BC')
                  EPS(jt) = 0.001_dp
                  PHI(jt) = 0.1_dp
               CASE('DU')
                  EPS(jt) = 0.001_dp
                  PHI(jt) = 0.1_dp
                  FHH(jt) = .true.
               END SELECT
               IF (MATCH_WILD( 'WSOC*', base_name(jt)) ) THEN
                  EPS(jt) = 0.5_dp
                  PHI(jt) = 0.2_dp
               ENDIF
            ENDDO
         CASE DEFAULT
            CALL error_bi('coupling for ' // aer_stream // ' not available' &
                 , substr)
         END SELECT

! Modified to account for MADE3 case (skip insoluble modes)
         DO ji = 1,nspec
            DO jn = 1,ktot
               DO jt = 1,ntrac
                  IF ((TRIM(ti_gp(jt)%tp%ident%basename) == TRIM(base_name(ji))) .AND. &
                       (TRIM(ti_gp(jt)%tp%ident%subname) == TRIM(sol_name(jn)))) THEN
                     aer_idx(jn,ji) = jt
                     EXIT
                  ENDIF
               ENDDO
            ENDDO
         ENDDO
         DEALLOCATE(sol_name)

! Print a summary table
         if (p_parallel_io) then
            write(*,*) "======= ABDUL-RAZZAK AND GHAN COEFFICIENTS ======="
            write(*,*) "-Species----------eps-----nu------phi-----kappa---"
            do ji=1,ksol
               do jt=1,nspec
                  IF (aer_idx(ji,jt) /= 0) &
                     write(*, '(1X,A15,4(2X,F6.3))') ti_gp(aer_idx(ji,jt))%tp%ident%fullname // "_" // &
                                                     ti_gp(aer_idx(ji,jt))%tp%ident%basename, &
                                                     EPS(jt),NU(jt),PHI(jt),KAPPA(jt)
               enddo
            enddo
            write(*,*) "------------------------------------"
            do ji=1,ksol
               IF (num_idx(ji) /= 0) &
                    print *,"Found number tracer ",ti_gp(num_idx(ji))%tp%ident%fullname
            enddo
            write(*,*) "------------------------------------"
         endif

         ALLOCATE (S_crit_in(ksol))
         SELECT CASE(sup_sat_scheme)
         CASE(1)
            CALL info_bi(&
                 'critical supersaturation calculated according to ARG' &
                 , substr)
            DO ji=1,ksol
               S_crit_in(ji)%ptr => NULL()
            ENDDO
         CASE(2)
            CALL info_bi(&
                 'critical supersaturation from gmxe_gp channel object', substr)
            if (TRIM(aer_stream) /= 'gmxe_gp') &
                 CALL error_bi (&
                 'sup_sat_scheme = 2 only available with gmxe_gp', substr)
            DO ji=1,ksol
               write (str_num,'(I1)') ji
               CALL get_channel_object(status, aer_stream, 'DIAGAER_M'//TRIM(ADJUSTL(str_num))//'_22', &
                                       p3=S_crit_in(ji)%ptr)
               IF (status /= 0) &
                    CALL error_bi (&
                    'critical supersaturation not found in the aerosol channel!', substr)
            END DO
         CASE(3)
            CALL info_bi( 'critical supersaturation calculated according to Petters and Kreidenweis', substr)
            DO ji=1,ksol
               S_crit_in(ji)%ptr => NULL()
            ENDDO
         END SELECT
      ENDIF    ! end ARG
   END IF

   IF (cloud_param == 1 ) THEN
      IF ( (i_cdnc_cpl /= i_cdnc_calc) .AND. &
           (i_cdnc_cpl /= 0)           .AND. &
           (i_cdnc_calc /= 9)         ) THEN
         CALL error_bi( &
              'Coupling between CDNC calculation and feedback to basemodel cloud routine inconsistent!', substr)
      endif
   ELSE
      IF ( i_cdnc_cpl /= 0) &
           CALL INFO_BI( &
           'Cloud droplet number concentration is determined independent of the simplified parameterisations!',substr)
   END IF

   IF (cloud_param == 3) THEN
      ! check for M7 stream / tracers
      CALL get_channel_object(status, aer_stream, 'sigma', p1=sigma)
      IF ( (status /= 0) .or. (SIZE(sigma) /=7) ) &
           CALL error_bi(&
           'This cloud routine can only work with M7 or '//       &
           'equivalently structured aerosol submodels!', substr)
      CALL get_channel_object(status, aer_stream,&
           'wetradius', p4=wetradius)
      IF (status /= 0) &
           CALL error_bi('wet radius channel object not found !', substr)
      ! collect tracer idts from M7 or equivalent aerosol module
      CALL GET_TRACER(status, GPTRSTR, 'SS',  subname='cs', idx=idt_msscs)
      CALL GET_TRACER(status, GPTRSTR, 'DU',  subname='cs', idx=idt_mducs)
      CALL GET_TRACER(status, GPTRSTR, 'SO4', subname='cs', idx=idt_ms4cs)
      CALL GET_TRACER(status, GPTRSTR, 'BC',  subname='cs', idx=idt_mbccs)
      CALL GET_TRACER(status, GPTRSTR, 'OC',  subname='cs', idx=idt_moccs)
      CALL GET_TRACER(status, GPTRSTR, 'SS',  subname='as', idx=idt_mssas)
      CALL GET_TRACER(status, GPTRSTR, 'DU',  subname='as', idx=idt_mduas)
      CALL GET_TRACER(status, GPTRSTR, 'SO4', subname='as', idx=idt_ms4as)
      CALL GET_TRACER(status, GPTRSTR, 'BC',  subname='as', idx=idt_mbcas)
      CALL GET_TRACER(status, GPTRSTR, 'OC',  subname='as', idx=idt_mocas)
      CALL GET_TRACER(status, GPTRSTR, 'SO4', subname='ks', idx=idt_ms4ks)
      CALL GET_TRACER(status, GPTRSTR, 'BC',  subname='ks', idx=idt_mbcks)
      CALL GET_TRACER(status, GPTRSTR, 'OC',  subname='ks', idx=idt_mocks)
      CALL GET_TRACER(status, GPTRSTR, 'BC',  subname='ki', idx=idt_mbcki)
      CALL GET_TRACER(status, GPTRSTR, 'OC',  subname='ki', idx=idt_mocki)
      
      CALL GET_TRACER(status, GPTRSTR, 'N',   subname='ki', idx=idt_nki)
      CALL GET_TRACER(status, GPTRSTR, 'N',   subname='ai', idx=idt_nai)
      CALL GET_TRACER(status, GPTRSTR, 'N',   subname='ci', idx=idt_nci)
      CALL GET_TRACER(status, GPTRSTR, 'N',   subname='ks', idx=idt_nks)
      CALL GET_TRACER(status, GPTRSTR, 'N',   subname='as', idx=idt_nas)
      CALL GET_TRACER(status, GPTRSTR, 'N',   subname='cs', idx=idt_ncs)

      IF (idt_ms4cs == 0)&
           CALL GET_TRACER(status, GPTRSTR, 'SO4mm', subname='cs', idx=idt_ms4cs)
      IF (idt_ms4as == 0)&
           CALL GET_TRACER(status, GPTRSTR, 'SO4mm', subname='as', idx=idt_ms4as)
      IF (idt_ms4ks == 0)&
           CALL GET_TRACER(status, GPTRSTR, 'SO4mm', subname='ks', idx=idt_ms4ks)
      nfrzmod = 3
   ENDIF
!-------------------------------------

!=====================================
   IF ((cloud_param >= 4) .AND. (cloud_param <= 6)) THEN

    CALL INFO_BI("Searching for coupling parameters for CLOUDPARAM = 4,5,6 !"&
         ,substr)
      CALL get_channel_object(status, aer_stream, 'sigma', p1=sigma)
      CALL channel_halt(substr, status)           
      CALL get_channel_object(status, aer_stream, 'crdiv_mid', p1=crdiv)

      IF (status /= 0.AND.(MATCH_WILD( 'made*', aer_stream))) THEN
         CALL INFO_BI('crdiv_mid channel not found for MADE/MADE3 but not required', substr)
      ELSE
         CALL channel_halt(substr, status)
      ENDIF
      CALL get_channel_object(status, aer_stream, 'wetradius', p4=wetradius)
      CALL channel_halt(substr, status)
      CALL get_channel_object(status, aer_stream, 'dryradius', p4=dryradius)
      CALL channel_halt(substr, status)


      ! perpare coupling of aerosols for usage in clouds
      CALL CPL_SPECIES_SI
      ! create a new channel with objects that are required for the 
      ! exchange of information between the aerosol and the cloud scheme of
      ! Lohmann et al., ACP, 2010 or Kuebbeler et al., ACP, 2014.

      modname = 'cloud_aer'
      if ((cloud_param == 4 .OR. cloud_param == 6)) &
           modname = 'cloud_aer_lohmann'
      if (cloud_param == 5) modname = 'cloud_aer_kuebbel'

      CALL new_channel(status, modname, reprid=GP_3D_MID)
      CALL channel_halt(substr, status)

      CALL new_channel_object(status, modname, 'cdncact_cv', &
        p3 = cdncact_cv)
      CALL channel_halt(substr, status)
      CALL new_attribute(status, modname, 'cdncact_cv'    , &
        'long_name', c='convective activated cloud droplets')
      CALL channel_halt(substr, status)
      CALL new_attribute(status, modname, 'cdncact_cv'    , &
        'units', c='1/cm^3' )
      CALL channel_halt(substr, status)
      
      CALL new_channel_object(status, modname, 'nbcsol_strat', &
        p3 = nbcsol_strat)
      CALL channel_halt(substr, status)
      CALL new_attribute(status, modname, 'nbcsol_strat'    , &
        'long_name', &
        c='number of BC particles for stratiform clouds in soluble modes')
      CALL channel_halt(substr, status)
      CALL new_attribute(status, modname, 'nbcsol_strat'    , &
        'units', c='1/m^3' )
      CALL channel_halt(substr, status)
      
      CALL new_channel_object(status, modname, 'ndusol_strat', &
        p3 = ndusol_strat)
      CALL channel_halt(substr, status)
      CALL new_attribute(status, modname, 'ndusol_strat'    , &
        'long_name', &
        c='number of dust particles for stratiform clouds in soluble modes')
      CALL channel_halt(substr, status)
      CALL new_attribute(status, modname, 'ndusol_strat'    , &
        'units', c='1/m^3' )
      CALL channel_halt(substr, status)

      CALL new_channel_object(status, modname, 'nbcsol_cirrus', &
        p3 = nbcsol_cirrus)
      CALL channel_halt(substr, status)
      CALL new_attribute(status, modname, 'nbcsol_cirrus'    , &
        'long_name', &
        c='number of BC particles for cirrus clouds in soluble modes')
      CALL channel_halt(substr, status)
      CALL new_attribute(status, modname, 'nbcsol_cirrus'    , &
        'units', c='1/m^3' )
      CALL channel_halt(substr, status)

      CALL new_channel_object(status, modname, 'nbctagsol_cirrus', &
        p3 = nbctagsol_cirrus)
      CALL channel_halt(substr, status)
      CALL new_attribute(status, modname, 'nbctagsol_cirrus'    , &
        'long_name', &
        c='number of tagged BC particles for cirrus clouds in soluble modes')
      CALL channel_halt(substr, status)
      CALL new_attribute(status, modname, 'nbctagsol_cirrus'    , &
        'units', c='1/m^3' )
      CALL channel_halt(substr, status)

      CALL new_channel_object(status, modname, 'ndusol_cirrus', &
        p3 = ndusol_cirrus)
      CALL channel_halt(substr, status)
      CALL new_attribute(status, modname, 'ndusol_cirrus'    , &
        'long_name', &
        c='number of dust particles for cirrus clouds in soluble modes')
      CALL channel_halt(substr, status)
      CALL new_attribute(status, modname, 'ndusol_cirrus'    , &
        'units', c='1/m^3' )
      CALL channel_halt(substr, status)

      !Soluble OC are needed for ice nucleation when BN09 + PDA13 spectrum is used 
      CALL new_channel_object(status, modname, 'nocsolks', p3 = nocsolks)
      CALL channel_halt(substr, status)
      CALL new_attribute(status, modname, 'nocsolks'    , &
        'long_name', c='number concentration of OC in Aitken soluble mode')
      CALL channel_halt(substr, status)
      CALL new_attribute(status, modname, 'nocsolks' , 'units', c='1/m^3' )
      CALL channel_halt(substr, status)

      CALL new_channel_object(status, modname, 'nocsolas', p3 = nocsolas)
      CALL channel_halt(substr, status)
      CALL new_attribute(status, modname, 'nocsolas'    , &
        'long_name', c='number concentration of OC in accumulation soluble mode')
      CALL channel_halt(substr, status)
      CALL new_attribute(status, modname, 'nocsolas' , 'units', c='1/m^3' )
      CALL channel_halt(substr, status)

      CALL new_channel_object(status, modname, 'nocsolcs', p3 = nocsolcs)
      CALL channel_halt(substr, status)
      CALL new_attribute(status, modname, 'nocsolcs'    , &
        'long_name', c='number concentration of OC in coarse soluble mode')
      CALL channel_halt(substr, status)
      CALL new_attribute(status, modname, 'nocsolcs' , 'units', c='1/m^3' )
      CALL channel_halt(substr, status)

      !Insoluble OC are needed for ice nucleation when BN09 + PDA08 spectrum is used 
      CALL new_channel_object(status, modname, 'nocinsol', p3 = nocinsol)
      CALL channel_halt(substr, status)
      CALL new_attribute(status, modname, 'nocinsol'    , &
        'long_name', c='number concentration of OC in Aitken insoluble mode')
      CALL channel_halt(substr, status)
      CALL new_attribute(status, modname, 'nocinsol' , 'units', c='1/m^3' )
      CALL channel_halt(substr, status)
      
      CALL new_channel_object(status, modname, 'nbcinsol', &
        p3 = nbcinsol)
      CALL channel_halt(substr, status)
      CALL new_attribute(status, modname, 'nbcinsol'    , &
        'long_name', &
        c='number of BC particles for stratiform clouds in insoluble modes')
      CALL channel_halt(substr, status)
      CALL new_attribute(status, modname, 'nbcinsol'    , &
        'units', c='1/m^3' )
      CALL channel_halt(substr, status)

      CALL new_channel_object(status, modname, 'nbctaginsol', &
        p3 = nbctaginsol)
      CALL channel_halt(substr, status)
      CALL new_attribute(status, modname, 'nbctaginsol'    , &
        'long_name', &
        c='number of tagged BC particles for stratiform clouds in insoluble modes')
      CALL channel_halt(substr, status)
      CALL new_attribute(status, modname, 'nbctaginsol'    , &
        'units', c='1/m^3' )
      CALL channel_halt(substr, status)

      CALL new_channel_object(status, modname, 'nduinsolai', &
        p3 = nduinsolai)
      CALL channel_halt(substr, status)
      CALL new_attribute(status, modname, 'nduinsolai'    , &
        'long_name', &
        c='number of dust particles for stratiform clouds in insoluble acc mode')
      CALL channel_halt(substr, status)
      CALL new_attribute(status, modname, 'nduinsolai'    , &
        'units', c='1/m^3' )
      CALL channel_halt(substr, status)
      
      CALL new_channel_object(status, modname, 'nduinsolci', &
        p3 = nduinsolci)
      CALL channel_halt(substr, status)
      CALL new_attribute(status, modname, 'nduinsolci'    , &
        'long_name', &
        c='number of dust particles for stratiform clouds in insoluble coarse mode')
      CALL channel_halt(substr, status)
      CALL new_attribute(status, modname, 'nduinsolci'    , &
        'units', c='1/m^3' )
      CALL channel_halt(substr, status)

      CALL new_channel_object(status, modname, 'naerinsol', &
        p3 = naerinsol)
      CALL channel_halt(substr, status)
      CALL new_attribute(status, modname, 'naerinsol'    , &
        'long_name', &
        c='total number of particles for stratiform clouds in insoluble modes')
      CALL channel_halt(substr, status)
      CALL new_attribute(status, modname, 'naerinsol'    , &
        'units', c='1/m^3' )
      CALL channel_halt(substr, status)

      CALL new_channel_object(status, modname, 'naersol', &
        p3 = naersol)
      CALL channel_halt(substr, status)
      CALL new_attribute(status, modname, 'naersol'    , &
        'long_name', &
        c='total number of particles for stratiform clouds in soluble modes')
      CALL channel_halt(substr, status)
      CALL new_attribute(status, modname, 'naersol'    , &
        'units', c='1/m^3' )
      CALL channel_halt(substr, status)

      CALL new_channel_object(status, modname, 'twc_conv', &
        p3 = twc_conv, lrestreq=.TRUE.)
      CALL channel_halt(substr, status)
      CALL new_attribute(status, modname, 'twc_conv'    , &
        'long_name', &
        c='total convective cloud water')
      CALL channel_halt(substr, status)
      CALL new_attribute(status, modname, 'twc_conv'    , &
        'units', c='????' )
      CALL channel_halt(substr, status)
      
      CALL new_channel_object(status, modname, 'conv_time', &
        p3 = conv_time, lrestreq=.TRUE.)
      CALL channel_halt(substr, status)
      CALL new_attribute(status, modname, 'conv_time'    , &
        'long_name', &
        c='total convective time')
      CALL channel_halt(substr, status)
      CALL new_attribute(status, modname, 'conv_time'    , &
        'units', c='????' )
      CALL channel_halt(substr, status)

   IF (cloud_param == 4) THEN

       !create new channel for some variables related to ice nucleation
       CALL new_channel(status, 'cloud_ice', reprid=GP_3D_MID)
       CALL channel_halt(substr, status)
    
       CALL new_channel_object(status, 'cloud_ice', 'newIC_cirri', p3=newIC_cirri)
       CALL channel_halt(substr, status)
       CALL new_attribute(status, 'cloud_ice', 'newIC_cirri', &
         'long_name', c='new ice crystals formed in the cirrus regime' )
       CALL channel_halt(substr, status)
       CALL new_attribute(status, 'cloud_ice', 'newIC_cirri', 'units', c='m-3' )
       CALL channel_halt(substr, status)

       CALL new_channel_object(status, 'cloud_ice', 'newICR_cirri', p3=newICR_cirri)
       CALL channel_halt(substr, status)
       CALL new_attribute(status, 'cloud_ice', 'newICR_cirri', &
         'long_name', c='radius of new ice crystals formed in the cirrus regime' )
       CALL channel_halt(substr, status)
       CALL new_attribute(status, 'cloud_ice', 'newICR_cirri', 'units', c='m' )
       CALL channel_halt(substr, status)

       CALL new_channel_object(status, 'cloud_ice', 'newIC_cnt_therm', p3=newIC_cnt_therm)
       CALL channel_halt(substr, status)
       CALL new_attribute(status, 'cloud_ice', 'newIC_cnt_therm', &
         'long_name', c='new ice crystals formed in the mixed-phase regime via contact and thermophoresis nucl.' )
       CALL channel_halt(substr, status)
       CALL new_attribute(status, 'cloud_ice', 'newIC_cnt_therm', 'units', c='m-3' )
       CALL channel_halt(substr, status)

       CALL new_channel_object(status, 'cloud_ice', 'newIC_imm', p3=newIC_imm)
       CALL channel_halt(substr, status)
       CALL new_attribute(status, 'cloud_ice', 'newIC_imm', &
         'long_name', c='new ice crystals formed in the mixed-phase regime via immersion nucleation' )
       CALL channel_halt(substr, status)
       CALL new_attribute(status, 'cloud_ice', 'newIC_imm', 'units', c='m-3' )
       CALL channel_halt(substr, status)

       CALL new_channel_object(status, 'cloud_ice', 'newIC_mix', p3=newIC_mix)
       CALL channel_halt(substr, status)
       CALL new_attribute(status, 'cloud_ice', 'newIC_mix', &
         'long_name', c='new ice crystals formed in mixed-phase regime (cnt+imm+therm)' )
       CALL channel_halt(substr, status)
       CALL new_attribute(status, 'cloud_ice', 'newIC_mix', 'units', c='m-3' )
       CALL channel_halt(substr, status)

       CALL new_channel_object(status, 'cloud_ice', 'newIC_mix_fin', p3=newIC_mix_fin)
       CALL channel_halt(substr, status)
       CALL new_attribute(status, 'cloud_ice', 'newIC_mix_fin', &
         'long_name', c='FINAL value of new ice crystals formed in the mixed-phase regime (cnt+imm+therm' )
       CALL channel_halt(substr, status)
       CALL new_attribute(status, 'cloud_ice', 'newIC_mix_fin', 'units', c='m-3' )
       CALL channel_halt(substr, status)

       CALL new_channel_object(status, 'cloud_ice', 'w_sub', p3=w_sub)
       CALL channel_halt(substr, status)
       CALL new_attribute(status, 'cloud_ice', 'w_sub', &
         'long_name', c='subgrid vertical velocity (turbulent contribution)' )
       CALL channel_halt(substr, status)
       CALL new_attribute(status, 'cloud_ice', 'w_sub', 'units', c='cm/s' )
       CALL channel_halt(substr, status)

       CALL new_channel_object(status, 'cloud_ice', 'w_ave', p3=w_ave)
       CALL channel_halt(substr, status)
       CALL new_attribute(status, 'cloud_ice', 'w_ave', &
         'long_name', c='average (resolved) vertical velocity' )
       CALL channel_halt(substr, status)
       CALL new_attribute(status, 'cloud_ice', 'w_ave', 'units', c='cm/s' )
       CALL channel_halt(substr, status)

       CALL new_channel_object(status, 'cloud_ice', 'w_grid', p3=w_grid)
       CALL channel_halt(substr, status)
       CALL new_attribute(status, 'cloud_ice', 'w_grid', &
         'long_name', c='vertical velocity in the grid cell (= w_sub + w_ave)' )
       CALL channel_halt(substr, status)
       CALL new_attribute(status, 'cloud_ice', 'w_grid', 'units', c='cm/s' )
       CALL channel_halt(substr, status)

          !Input of BN09
          CALL new_channel_object(status, 'cloud_ice', 'sigwBN', p3=sigwBN)
          CALL channel_halt(substr, status)
          CALL new_attribute(status, 'cloud_ice', 'sigwBN', &
            'long_name', c='width of the distribution of updraft velocity (input for BN09)')
          CALL channel_halt(substr, status)
          CALL new_attribute(status, 'cloud_ice', 'sigwBN', 'units', c='m/s' )
          CALL channel_halt(substr, status)

          CALL new_channel_object(status, 'cloud_ice', 'ndropBN', p3=ndropBN)
          CALL channel_halt(substr, status)
          CALL new_attribute(status, 'cloud_ice', 'ndropBN', &
            'long_name', c='total cloud droplet number concentration (input for BN09)')
          CALL channel_halt(substr, status)
          CALL new_attribute(status, 'cloud_ice', 'ndropBN', 'units', c='m-3' )
          CALL channel_halt(substr, status)

          CALL new_channel_object(status, 'cloud_ice', 'dsulfBN', p3=dsulfBN)
          CALL channel_halt(substr, status)
          CALL new_attribute(status, 'cloud_ice', 'dsulfBN', &
            'long_name', c='geometric mean diameter of sulfate (input for BN09)')
          CALL channel_halt(substr, status)
          CALL new_attribute(status, 'cloud_ice', 'dsulfBN', 'units', c='m' )
          CALL channel_halt(substr, status)

          CALL new_channel_object(status, 'cloud_ice', 'ndust_aiBN', p3=ndust_aiBN)
          CALL channel_halt(substr, status)
          CALL new_attribute(status, 'cloud_ice', 'ndust_aiBN', &
            'long_name', c='number concentration of dust DU acc. insoluble mode (input BN09 + CNT,PDA08,PDA13)')
          CALL channel_halt(substr, status)
          CALL new_attribute(status, 'cloud_ice', 'ndust_aiBN', 'units', c='m-3' )
          CALL channel_halt(substr, status)

          CALL new_channel_object(status, 'cloud_ice', 'ddust_aiBN', p3=ddust_aiBN)
          CALL channel_halt(substr, status)
          CALL new_attribute(status, 'cloud_ice', 'ddust_aiBN', &
            'long_name', c='geometric mean diameter of DU acc. insoluble mode (input BN09 + CNT,PDA08,PDA13)')
          CALL channel_halt(substr, status)
          CALL new_attribute(status, 'cloud_ice', 'ddust_aiBN', 'units', c='m' )
          CALL channel_halt(substr, status)

          CALL new_channel_object(status, 'cloud_ice', 'ndust_ciBN', p3=ndust_ciBN)
          CALL channel_halt(substr, status)
          CALL new_attribute(status, 'cloud_ice', 'ndust_ciBN', &
            'long_name', c='number concentration of dust DU coarse insoluble mode (input BN09 + CNT,PDA08,PDA13)')
          CALL channel_halt(substr, status)
          CALL new_attribute(status, 'cloud_ice', 'ndust_ciBN', 'units', c='m-3' )
          CALL channel_halt(substr, status)

          CALL new_channel_object(status, 'cloud_ice', 'ddust_ciBN', p3=ddust_ciBN)
          CALL channel_halt(substr, status)
          CALL new_attribute(status, 'cloud_ice', 'ddust_ciBN', &
            'long_name', c='geometric mean diameter of DU coarse insoluble mode (input BN09 + CNT,PDA08,PDA13)')
          CALL channel_halt(substr, status)
          CALL new_attribute(status, 'cloud_ice', 'ddust_ciBN', 'units', c='m' )
          CALL channel_halt(substr, status)

          CALL new_channel_object(status, 'cloud_ice', 'norgBN', p3=norgBN)
          CALL channel_halt(substr, status)
          CALL new_attribute(status, 'cloud_ice', 'norgBN', &
            'long_name', c='number concentration of insoluble organics OC (input BN09 + PDA08)')
          CALL channel_halt(substr, status)
          CALL new_attribute(status, 'cloud_ice', 'norgBN', 'units', c='m-3' )
          CALL channel_halt(substr, status)

          CALL new_channel_object(status, 'cloud_ice', 'dorgBN', p3=dorgBN)
          CALL channel_halt(substr, status)
          CALL new_attribute(status, 'cloud_ice', 'dorgBN', &
            'long_name', c='geometric mean diameter of insoluble OC (input BN09 + PDA08)')
          CALL channel_halt(substr, status)
          CALL new_attribute(status, 'cloud_ice', 'dorgBN', 'units', c='m' )
          CALL channel_halt(substr, status)

          CALL new_channel_object(status, 'cloud_ice', 'nsootBN', p3=nsootBN)
          CALL channel_halt(substr, status)
          CALL new_attribute(status, 'cloud_ice', 'nsootBN', &
            'long_name', c='number concentration of insoluble soot BC (input BN09 + CNT,PDA08,PDA13)')
          CALL channel_halt(substr, status)
          CALL new_attribute(status, 'cloud_ice', 'nsootBN', 'units', c='m-3' )
          CALL channel_halt(substr, status)

          CALL new_channel_object(status, 'cloud_ice', 'dsootBN', p3=dsootBN)
          CALL channel_halt(substr, status)
          CALL new_attribute(status, 'cloud_ice', 'dsootBN', &
            'long_name', c='geometric mean diameter of insoluble BC (input BN09 + CNT,PDA08,PDA13)')
          CALL channel_halt(substr, status)
          CALL new_attribute(status, 'cloud_ice', 'dsootBN', 'units', c='m' )
          CALL channel_halt(substr, status)

          CALL new_channel_object(status, 'cloud_ice', 'nsolo_ksBN', p3=nsolo_ksBN)
          CALL channel_halt(substr, status)
          CALL new_attribute(status, 'cloud_ice', 'nsolo_ksBN', &
            'long_name', c='number concentration of OC Aitken soluble mode (input BN09 + PDA13)')
          CALL channel_halt(substr, status)
          CALL new_attribute(status, 'cloud_ice', 'nsolo_ksBN', 'units', c='m-3' )
          CALL channel_halt(substr, status)

          CALL new_channel_object(status, 'cloud_ice', 'dsolo_ksBN', p3=dsolo_ksBN)
          CALL channel_halt(substr, status)
          CALL new_attribute(status, 'cloud_ice', 'dsolo_ksBN', &
            'long_name', c='geometric mean diameter of OC ait. soluble mode (input BN09 + PDA13)')
          CALL channel_halt(substr, status)
          CALL new_attribute(status, 'cloud_ice', 'dsolo_ksBN', 'units', c='m' )
          CALL channel_halt(substr, status)

          CALL new_channel_object(status, 'cloud_ice', 'nsolo_asBN', p3=nsolo_asBN)
          CALL channel_halt(substr, status)
          CALL new_attribute(status, 'cloud_ice', 'nsolo_asBN', &
            'long_name', c='number concentration of OC acc. soluble mode (input BN09 + PDA13)')
          CALL channel_halt(substr, status)
          CALL new_attribute(status, 'cloud_ice', 'nsolo_asBN', 'units', c='m-3' )
          CALL channel_halt(substr, status)

          CALL new_channel_object(status, 'cloud_ice', 'dsolo_asBN', p3=dsolo_asBN)
          CALL channel_halt(substr, status)
          CALL new_attribute(status, 'cloud_ice', 'dsolo_asBN', &
            'long_name', c='geometric mean diameter of OC acc. soluble mode (input BN09 + PDA13)')
          CALL channel_halt(substr, status)
          CALL new_attribute(status, 'cloud_ice', 'dsolo_asBN', 'units', c='m' )
          CALL channel_halt(substr, status)

          CALL new_channel_object(status, 'cloud_ice', 'nsolo_csBN', p3=nsolo_csBN)
          CALL channel_halt(substr, status)
          CALL new_attribute(status, 'cloud_ice', 'nsolo_csBN', &
            'long_name', c='number concentration of OC coarse soluble mode (input BN09 + PDA13)')
          CALL channel_halt(substr, status)
          CALL new_attribute(status, 'cloud_ice', 'nsolo_csBN', 'units', c='m-3' )
          CALL channel_halt(substr, status)
 
          CALL new_channel_object(status, 'cloud_ice', 'dsolo_csBN', p3=dsolo_csBN)
          CALL channel_halt(substr, status)
          CALL new_attribute(status, 'cloud_ice', 'dsolo_csBN', &
            'long_name', c='geometric mean diameter of OC coarse soluble mode (input BN09 + PDA13)')
          CALL channel_halt(substr, status)
          CALL new_attribute(status, 'cloud_ice', 'dsolo_csBN', 'units', c='m' )
          CALL channel_halt(substr, status)

          !Output of BN09
          CALL new_channel_object(status, 'cloud_ice', 'smaxice_cirrusBN', p3=smaxice_cirrusBN)
          CALL channel_halt(substr, status)
          CALL new_attribute(status, 'cloud_ice', 'smaxice_cirrusBN', &
            'long_name', c='maximum supersaturation with respect to ice (output of BN09)')
          CALL channel_halt(substr, status)
          CALL new_attribute(status, 'cloud_ice', 'smaxice_cirrusBN', 'units', c='-' )
          CALL channel_halt(substr, status)

          CALL new_channel_object(status, 'cloud_ice', 'sc_ice_cirrusBN', p3=sc_ice_cirrusBN)
          CALL channel_halt(substr, status)
          CALL new_attribute(status, 'cloud_ice', 'sc_ice_cirrusBN', &
            'long_name', c='characteristic freezing point of the aerosol population (output of BN09)')
          CALL channel_halt(substr, status)
          CALL new_attribute(status, 'cloud_ice', 'sc_ice_cirrusBN', 'units', c='-' )
          CALL channel_halt(substr, status)

          CALL new_channel_object(status, 'cloud_ice', 'nlim_cirrusBN', p3=nlim_cirrusBN)
          CALL channel_halt(substr, status)
          CALL new_attribute(status, 'cloud_ice', 'nlim_cirrusBN', &
            'long_name', c='limiting IN concentration (output of BN09)')
          CALL channel_halt(substr, status)
          CALL new_attribute(status, 'cloud_ice', 'nlim_cirrusBN', 'units', c='m-3' )
          CALL channel_halt(substr, status)

          CALL new_channel_object(status, 'cloud_ice', 'nhet_cirrusBN', p3=nhet_cirrusBN)
          CALL channel_halt(substr, status)
          CALL new_attribute(status, 'cloud_ice', 'nhet_cirrusBN', &
            'long_name', c='ice crystal concentration by het freezing (output of BN09)')
          CALL channel_halt(substr, status)
          CALL new_attribute(status, 'cloud_ice', 'nhet_cirrusBN', 'units', c='m-3' )
          CALL channel_halt(substr, status)

          CALL new_channel_object(status, 'cloud_ice', 'nice_cirrusBN', p3=nice_cirrusBN)
          CALL channel_halt(substr, status)
          CALL new_attribute(status, 'cloud_ice', 'nice_cirrusBN', &
            'long_name', c='nucleated ice crystal number concentration (output of BN09)')
          CALL channel_halt(substr, status)
          CALL new_attribute(status, 'cloud_ice', 'nice_cirrusBN', 'units', c='m-3' )
          CALL channel_halt(substr, status)

          CALL new_channel_object(status, 'cloud_ice', 'dice_cirrusBN', p3=dice_cirrusBN)
          CALL channel_halt(substr, status)
          CALL new_attribute(status, 'cloud_ice', 'dice_cirrusBN', &
            'long_name', c='diameter of ice crystals (output of BN09)')
          CALL channel_halt(substr, status)
          CALL new_attribute(status, 'cloud_ice', 'dice_cirrusBN', 'units', c='m' )
          CALL channel_halt(substr, status)

          CALL new_channel_object(status, 'cloud_ice', 'wsub_cirrusBN', p3=sigwpre_cirrusBN)
          CALL channel_halt(substr, status)
          CALL new_attribute(status, 'cloud_ice', 'wsub_cirrusBN', &
            'long_name', c='Subgrid vertical velocity used in BN09' )
          CALL channel_halt(substr, status)
          CALL new_attribute(status, 'cloud_ice', 'wsub_cirrusBN', 'units', c='m/s' )
          CALL channel_halt(substr, status)

          !Output of BN09
          CALL new_channel_object(status, 'cloud_ice', 'smaxice_immBN', p3=smaxice_immBN)
          CALL channel_halt(substr, status)
          CALL new_attribute(status, 'cloud_ice', 'smaxice_immBN', &
            'long_name', c='maximum supersaturation with respect to ice (output of BN09)')
          CALL channel_halt(substr, status)
          CALL new_attribute(status, 'cloud_ice', 'smaxice_immBN', 'units', c='-' )
          CALL channel_halt(substr, status)

          CALL new_channel_object(status, 'cloud_ice', 'sc_ice_immBN', p3=sc_ice_immBN)
          CALL channel_halt(substr, status)
          CALL new_attribute(status, 'cloud_ice', 'sc_ice_immBN', &
            'long_name', c='characteristic freezing point of the aerosol population (output of BN09)')
          CALL channel_halt(substr, status)
          CALL new_attribute(status, 'cloud_ice', 'sc_ice_immBN', 'units', c='-' )
          CALL channel_halt(substr, status)

          CALL new_channel_object(status, 'cloud_ice', 'nice_immBN', p3=nice_immBN)
          CALL channel_halt(substr, status)
          CALL new_attribute(status, 'cloud_ice', 'nice_immBN', &
            'long_name', c='nucleated ice crystal number concentration (output of BN09)')
          CALL channel_halt(substr, status)
          CALL new_attribute(status, 'cloud_ice', 'nice_immBN', 'units', c='m-3' )
          CALL channel_halt(substr, status)

          CALL new_channel_object(status, 'cloud_ice', 'dice_immBN', p3=dice_immBN)
          CALL channel_halt(substr, status)
          CALL new_attribute(status, 'cloud_ice', 'dice_immBN', &
            'long_name', c='diameter of ice crystals (output of BN09)')
          CALL channel_halt(substr, status)
          CALL new_attribute(status, 'cloud_ice', 'dice_immBN', 'units', c='m' )
          CALL channel_halt(substr, status)

          CALL new_channel_object(status, 'cloud_ice', 'wsub_immBN', p3=sigwpre_immBN)
          CALL channel_halt(substr, status)
          CALL new_attribute(status, 'cloud_ice', 'wsub_immBN', &
            'long_name', c='Subgrid vertical velocity used in BN09' )
          CALL channel_halt(substr, status)
          CALL new_attribute(status, 'cloud_ice', 'wsub_immBN', 'units', c='m/s' )
          CALL channel_halt(substr, status)
        !END IF

   END IF   !cloud_param = 4

      IF (cloud_param == 5) then

         ! Coupling parameters from OROGW
         CALL INFO_BI("Searching for coupling parameters for CLOUDPARAM = 5 !" &
              , substr)
         CALL get_channel_object(status, 'orogw', 'w_gwd_kpr', p3=w_gwd_kpr)
         CALL channel_halt(substr, status)
         CALL get_channel_object(status, 'orogw', 'ampl_gwd', p3=ampl_gwd)
         CALL channel_halt(substr, status)
         CALL get_channel_object(status, 'orogw', 'mask_orogw', p2=mask_orogw)
         CALL channel_halt(substr, status)
         CALL get_channel_object(status, 'orogw', 'l_z', p2=l_z)
         CALL channel_halt(substr, status)

         ! Vertical velocities
         CALL new_channel_object(status, modname, 'vervel_ls', p3=vervel_ls)
         CALL channel_halt(substr, status)
         CALL new_attribute(status, modname, 'vervel_ls', &
              'long_name', c='Vertical velocity (large scale component)')
         CALL channel_halt(substr, status)
         CALL new_attribute(status, modname, 'vervel_ls', 'units', c='cm/s' )
         CALL channel_halt(substr, status)

         CALL new_channel_object(status, modname, 'vervel_tke', p3=vervel_tke)
         CALL channel_halt(substr, status)
         CALL new_attribute(status, modname, 'vervel_tke', &
              'long_name', c='Vertical velocity (TKE component)')
         CALL channel_halt(substr, status)
         CALL new_attribute(status, modname, 'vervel_tke', 'units', c='cm/s' )
         CALL channel_halt(substr, status)

         CALL new_channel_object(status, modname, 'vervel_gw', p3=vervel_gw)
         CALL channel_halt(substr, status)
         CALL new_attribute(status, modname, 'vervel_gw', &
              'long_name', c='Vertical velocity (orogr. waves component)')
         CALL channel_halt(substr, status)
         CALL new_attribute(status, modname, 'vervel_gw', 'units', c='cm/s' )
         CALL channel_halt(substr, status)

         CALL new_channel_object(status, modname, 'vervel_p18', p3=vervel_p18)
         CALL channel_halt(substr, status)
         CALL new_attribute(status, modname, 'vervel_p18', &
              'long_name', c='Vertical velocity (from Penner et al. 2018)')
         CALL channel_halt(substr, status)
         CALL new_attribute(status, modname, 'vervel_p18', 'units', c='cm/s' )
         CALL channel_halt(substr, status)

         ! CDNC for comparison with observational data
         CALL new_channel_object(status, modname, 'CDNC_insitu', p3=cdnc_insitu)
         CALL channel_halt(substr, status)
         CALL new_attribute(status, modname, 'CDNC_insitu', &
              'long_name', c='In-cloud CDNC (for LWC>0.01 g/m3)')
         CALL channel_halt(substr, status)
         CALL new_attribute(status, modname, 'CDNC_insitu', 'units', c='m-3' )
         CALL channel_halt(substr, status)

         ! --------------------------------------------------------------------

         ! Cirrus-specific diagnostics for aircraft data by M. Kraemer (FZJ)

         ! Temperature dimension (1-K bins from 182. to 243 K)
         DO ibin = 1, N_CIRRUS_TEMP
            CIRRUS_TEMP(ibin) = 182._dp + ibin - 1
         END DO
         CALL new_dimension(status, DIMID_CIRRUS_TEMP, &
              'CIRRUS_TEMP', N_CIRRUS_TEMP)
         CALL channel_halt(substr, status)
         CALL add_dimension_variable(status, DIMID_CIRRUS_TEMP, &
              'CIRRUS_TEMP', CIRRUS_TEMP)
         CALL channel_halt(substr, status)
         CALL add_dimension_variable_att(status, DIMID_CIRRUS_TEMP, &
              'CIRRUS_TEMP', 'units', c='K')
         CALL add_dimension_variable_att(status, DIMID_CIRRUS_TEMP, &
              'CIRRUS_TEMP', 'ibin', c='bin +/- 0.5')
         CALL channel_halt(substr, status)

         ! Bin dimension for IWC (0.2 log bins from -2.9 to 2.9)
         DO ibin = 1, N_CIRRUS_BIN_IWC
            CIRRUS_BIN_IWC(ibin) = -2.9_dp + 0.2_dp * (ibin - 1)
            CIRRUS_IBIN_IWC(ibin) = CIRRUS_BIN_IWC(ibin) - 0.1_dp
         END DO
         CIRRUS_IBIN_IWC(N_CIRRUS_BIN_IWC + 1) = &
              CIRRUS_BIN_IWC(N_CIRRUS_BIN_IWC) + 0.1_dp
         CIRRUS_BIN_IWC(:) = 10**CIRRUS_BIN_IWC(:)
         CIRRUS_IBIN_IWC(:) = 10**CIRRUS_IBIN_IWC(:)
         CALL new_dimension(status, DIMID_CIRRUS_BIN_IWC, &
              'CIRRUS_BIN_IWC', N_CIRRUS_BIN_IWC)
         CALL channel_halt(substr, status)
         CALL add_dimension_variable(status, DIMID_CIRRUS_BIN_IWC, &
              'CIRRUS_BIN_IWC', CIRRUS_BIN_IWC)
         CALL channel_halt(substr, status)
         CALL add_dimension_variable_att(status, DIMID_CIRRUS_BIN_IWC, &
              'CIRRUS_BIN_IWC', 'units', c='ppmv') 
         CALL channel_halt(substr, status)
         CALL add_dimension_variable_att(status, DIMID_CIRRUS_BIN_IWC, &
              'CIRRUS_BIN_IWC', 'ibin', c='10^(log10(bin) +/- 0.1)') 
         CALL channel_halt(substr, status)
     
         CALL new_representation(status, REPR_CIRRUS_4D_IWC, &
              'REPR_CIRRUS_4D_IWC'    &    
              , rank = 4, link = 'xxxx', dctype = DC_GP               &    
              , dimension_ids = (/ DIMID_LON, DIMID_CIRRUS_TEMP       &    
              ,                    DIMID_CIRRUS_BIN_IWC, DIMID_LAT /) &    
              , ldimlen       = (/ nproma, AUTO, AUTO, ngpblks   /)   &
              , output_order  = (/ 3,1,4,2 /)                         &    
              , axis = 'XNNY'                                         &    
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
         cnt(:,2) = N_CIRRUS_TEMP
         memu(:,2) = N_CIRRUS_TEMP
         cnt(:,3) = N_CIRRUS_BIN_IWC
         memu(:,3) = N_CIRRUS_BIN_IWC

         CALL set_representation_decomp(status, REPR_CIRRUS_4D_IWC &
              , start, cnt, memu, meml, .FALSE., PIOTYPE_COL)
         CALL channel_halt(substr, status)

         DEALLOCATE(start) ; NULLIFY(start)
         DEALLOCATE(cnt)   ; NULLIFY(cnt)
         DEALLOCATE(meml)  ; NULLIFY(meml)
         DEALLOCATE(memu)  ; NULLIFY(memu)

         ! Bin dimension for Nice (0.2 log bins from -3.9 to 2.9)
         DO ibin = 1, N_CIRRUS_BIN_Nice
            CIRRUS_BIN_Nice(ibin) = -3.9_dp + 0.2_dp * (ibin - 1)
            CIRRUS_IBIN_Nice(ibin) = CIRRUS_BIN_Nice(ibin) - 0.1_dp
         END DO
         CIRRUS_IBIN_Nice(N_CIRRUS_BIN_Nice + 1) = &
              CIRRUS_BIN_Nice(N_CIRRUS_BIN_Nice) + 0.1_dp
         CIRRUS_BIN_Nice(:) = 10**CIRRUS_BIN_Nice(:)
         CIRRUS_IBIN_Nice(:) = 10**CIRRUS_IBIN_Nice(:)
         CALL new_dimension(status, DIMID_CIRRUS_BIN_Nice, &
              'CIRRUS_BIN_Nice', N_CIRRUS_BIN_Nice)
         CALL channel_halt(substr, status)
         CALL add_dimension_variable(status, DIMID_CIRRUS_BIN_Nice, &
              'CIRRUS_BIN_Nice', CIRRUS_BIN_Nice)
         CALL channel_halt(substr, status)
         CALL add_dimension_variable_att(status, DIMID_CIRRUS_BIN_Nice, &
              'CIRRUS_BIN_Nice', 'units', c='cm-3')
         CALL channel_halt(substr, status)
         CALL add_dimension_variable_att(status, DIMID_CIRRUS_BIN_Nice, &
              'CIRRUS_BIN_Nice', 'ibin', c='10^(log10(bin) +/- 0.1)') 
         CALL channel_halt(substr, status)
     
         CALL new_representation(status, REPR_CIRRUS_4D_Nice, &
              'REPR_CIRRUS_4D_Nice'    &    
              , rank = 4, link = 'xxxx', dctype = DC_GP               &    
              , dimension_ids = (/ DIMID_LON, DIMID_CIRRUS_TEMP       &    
              ,                    DIMID_CIRRUS_BIN_Nice, DIMID_LAT /) &    
              , ldimlen       = (/ nproma, AUTO, AUTO, ngpblks   /)   &
              , output_order  = (/ 3,1,4,2 /)                         &    
              , axis = 'XNNY'                                         &    
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
         cnt(:,2) = N_CIRRUS_TEMP
         memu(:,2) = N_CIRRUS_TEMP
         cnt(:,3) = N_CIRRUS_BIN_Nice
         memu(:,3) = N_CIRRUS_BIN_Nice

         CALL set_representation_decomp(status, REPR_CIRRUS_4D_Nice &
              , start, cnt, memu, meml, .FALSE., PIOTYPE_COL)
         CALL channel_halt(substr, status)

         DEALLOCATE(start) ; NULLIFY(start)
         DEALLOCATE(cnt)   ; NULLIFY(cnt)
         DEALLOCATE(meml)  ; NULLIFY(meml)
         DEALLOCATE(memu)  ; NULLIFY(memu)

         ! Bin dimension for Rice (0.12 log bins from 0.06 to 2.94)
         DO ibin = 1, N_CIRRUS_BIN_Rice
            CIRRUS_BIN_Rice(ibin) = 0.06_dp + 0.12_dp * (ibin - 1)
            CIRRUS_IBIN_Rice(ibin) = CIRRUS_BIN_Rice(ibin) - 0.06_dp
         END DO
         CIRRUS_IBIN_Rice(N_CIRRUS_BIN_Rice + 1) = &
              CIRRUS_BIN_Rice(N_CIRRUS_BIN_Rice) + 0.06_dp
         CIRRUS_BIN_Rice(:) = 10**CIRRUS_BIN_Rice(:)
         CIRRUS_IBIN_Rice(:) = 10**CIRRUS_IBIN_Rice(:)
         CALL new_dimension(status, DIMID_CIRRUS_BIN_Rice, &
              'CIRRUS_BIN_Rice', N_CIRRUS_BIN_Rice)
         CALL channel_halt(substr, status)
         CALL add_dimension_variable(status, DIMID_CIRRUS_BIN_Rice, &
              'CIRRUS_BIN_Rice', CIRRUS_BIN_Rice)
         CALL channel_halt(substr, status)
         CALL add_dimension_variable_att(status, DIMID_CIRRUS_BIN_Rice, &
              'CIRRUS_BIN_Rice', 'units', c='um')
         CALL channel_halt(substr, status)
         CALL add_dimension_variable_att(status, DIMID_CIRRUS_BIN_Rice, &
              'CIRRUS_BIN_Rice', 'ibin', c='10^(log10(bin) +/- 0.06)') 
         CALL channel_halt(substr, status)
     
         CALL new_representation(status, REPR_CIRRUS_4D_Rice, &
              'REPR_CIRRUS_4D_Rice'    &    
              , rank = 4, link = 'xxxx', dctype = DC_GP               &    
              , dimension_ids = (/ DIMID_LON, DIMID_CIRRUS_TEMP       &    
              ,                    DIMID_CIRRUS_BIN_Rice, DIMID_LAT /) &    
              , ldimlen       = (/ nproma, AUTO, AUTO, ngpblks   /)   &
              , output_order  = (/ 3,1,4,2 /)                         &    
              , axis = 'XNNY'                                         &    
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
         cnt(:,2) = N_CIRRUS_TEMP
         memu(:,2) = N_CIRRUS_TEMP
         cnt(:,3) = N_CIRRUS_BIN_Rice
         memu(:,3) = N_CIRRUS_BIN_Rice

         CALL set_representation_decomp(status, REPR_CIRRUS_4D_Rice &
              , start, cnt, memu, meml, .FALSE., PIOTYPE_COL)
         CALL channel_halt(substr, status)

         DEALLOCATE(start) ; NULLIFY(start)
         DEALLOCATE(cnt)   ; NULLIFY(cnt)
         DEALLOCATE(meml)  ; NULLIFY(meml)
         DEALLOCATE(memu)  ; NULLIFY(memu)

         ! Bin dimension for RHi (10-% bins from 10 to 240)
         DO ibin = 1, N_CIRRUS_BIN_RHi
            CIRRUS_BIN_RHi(ibin) = 10._dp + 10._dp * (ibin - 1)
            CIRRUS_IBIN_RHi(ibin) = CIRRUS_BIN_RHi(ibin) - 5._dp
         END DO
         CIRRUS_IBIN_RHi(N_CIRRUS_BIN_RHi + 1) = &
              CIRRUS_BIN_RHi(N_CIRRUS_BIN_RHi) + 5._dp
         CALL new_dimension(status, DIMID_CIRRUS_BIN_RHi, &
              'CIRRUS_BIN_RHi', N_CIRRUS_BIN_RHi)
         CALL channel_halt(substr, status)
         CALL add_dimension_variable(status, DIMID_CIRRUS_BIN_RHi, &
              'CIRRUS_BIN_RHi', CIRRUS_BIN_RHi)
         CALL channel_halt(substr, status)
         CALL add_dimension_variable_att(status, DIMID_CIRRUS_BIN_RHi, &
              'CIRRUS_BIN_RHi', 'units', c='%')
         CALL channel_halt(substr, status)
         CALL add_dimension_variable_att(status, DIMID_CIRRUS_BIN_RHi, &
              'CIRRUS_BIN_RHi', 'ibin', c='bin +/- 5') 
         CALL channel_halt(substr, status)
     
         CALL new_representation(status, REPR_CIRRUS_4D_RHi, &
              'REPR_CIRRUS_4D_RHi'    &    
              , rank = 4, link = 'xxxx', dctype = DC_GP               &    
              , dimension_ids = (/ DIMID_LON, DIMID_CIRRUS_TEMP       &    
              ,                    DIMID_CIRRUS_BIN_RHi, DIMID_LAT /) &    
              , ldimlen       = (/ nproma, AUTO, AUTO, ngpblks   /)   &
              , output_order  = (/ 3,1,4,2 /)                         &    
              , axis = 'XNNY'                                         &    
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
         cnt(:,2) = N_CIRRUS_TEMP
         memu(:,2) = N_CIRRUS_TEMP
         cnt(:,3) = N_CIRRUS_BIN_RHi
         memu(:,3) = N_CIRRUS_BIN_RHi

         CALL set_representation_decomp(status, REPR_CIRRUS_4D_RHi &
              , start, cnt, memu, meml, .FALSE., PIOTYPE_COL)
         CALL channel_halt(substr, status)

         DEALLOCATE(start) ; NULLIFY(start)
         DEALLOCATE(cnt)   ; NULLIFY(cnt)
         DEALLOCATE(meml)  ; NULLIFY(meml)
         DEALLOCATE(memu)  ; NULLIFY(memu)

         ! Bin dimension for VEL (0.1 log bins from 0.05 to 2.95)
         DO ibin = 1, N_CIRRUS_BIN_VEL
            CIRRUS_BIN_VEL(ibin) = 0.05_dp + 0.1_dp * (ibin - 1)
            CIRRUS_IBIN_VEL(ibin) = CIRRUS_BIN_VEL(ibin) - 0.05_dp
         END DO
         CIRRUS_IBIN_VEL(N_CIRRUS_BIN_VEL + 1) = &
              CIRRUS_BIN_VEL(N_CIRRUS_BIN_VEL) + 0.05_dp
         CIRRUS_BIN_VEL(:) = 10**CIRRUS_BIN_VEL(:)
         CIRRUS_IBIN_VEL(:) = 10**CIRRUS_IBIN_VEL(:)
         CALL new_dimension(status, DIMID_CIRRUS_BIN_VEL, &
              'CIRRUS_BIN_VEL', N_CIRRUS_BIN_VEL)
         CALL channel_halt(substr, status)
         CALL add_dimension_variable(status, DIMID_CIRRUS_BIN_VEL, &
              'CIRRUS_BIN_VEL', CIRRUS_BIN_VEL)
         CALL channel_halt(substr, status)
         CALL add_dimension_variable_att(status, DIMID_CIRRUS_BIN_VEL, &
              'CIRRUS_BIN_VEL', 'units', c='cm/s') 
         CALL channel_halt(substr, status)
         CALL add_dimension_variable_att(status, DIMID_CIRRUS_BIN_VEL, &
              'CIRRUS_BIN_VEL', 'ibin', c='10^(log10(bin) +/- 0.05)') 
         CALL channel_halt(substr, status)
     
         CALL new_representation(status, REPR_CIRRUS_4D_VEL, &
              'REPR_CIRRUS_4D_VEL'    &    
              , rank = 4, link = 'xxxx', dctype = DC_GP               &    
              , dimension_ids = (/ DIMID_LON, DIMID_CIRRUS_TEMP       &    
              ,                    DIMID_CIRRUS_BIN_VEL, DIMID_LAT /) &    
              , ldimlen       = (/ nproma, AUTO, AUTO, ngpblks   /)   &
              , output_order  = (/ 3,1,4,2 /)                         &    
              , axis = 'XNNY'                                         &    
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
         cnt(:,2) = N_CIRRUS_TEMP
         memu(:,2) = N_CIRRUS_TEMP
         cnt(:,3) = N_CIRRUS_BIN_VEL
         memu(:,3) = N_CIRRUS_BIN_VEL

         CALL set_representation_decomp(status, REPR_CIRRUS_4D_VEL &
              , start, cnt, memu, meml, .FALSE., PIOTYPE_COL)
         CALL channel_halt(substr, status)

         DEALLOCATE(start) ; NULLIFY(start)
         DEALLOCATE(cnt)   ; NULLIFY(cnt)
         DEALLOCATE(meml)  ; NULLIFY(meml)
         DEALLOCATE(memu)  ; NULLIFY(memu)

         ! New channel for cirrus diagnostics
         ! These variables should be saved as monthly mean averages
         ! To get the total number of points in each bin, multiply then by
         ! the number of timesteps in the given month (ndays * 24 * 3600 / dt)
         CALL new_channel(status, 'cloud_cirrus')
         CALL channel_halt(substr, status)

         CALL new_channel_object(status, 'cloud_cirrus', 'CIRRUS_IWC', &
              p4=CIRRUS_IWC, reprid=REPR_CIRRUS_4D_IWC)
         CALL channel_halt(substr, status)
         CALL new_attribute(status, 'cloud_cirrus', 'CIRRUS_IWC', &
              'long_name', c='PDF of ice water content (in-cloud)')
         CALL channel_halt(substr, status)
         CALL new_attribute(status, 'cloud_cirrus', &
              'CIRRUS_IWC', 'units', c='-' )
         CALL channel_halt(substr, status)

         CALL new_channel_object(status, 'cloud_cirrus', 'CIRRUS_Nice', &
              p4=CIRRUS_Nice, reprid=REPR_CIRRUS_4D_Nice)
         CALL channel_halt(substr, status)
         CALL new_attribute(status, 'cloud_cirrus', 'CIRRUS_Nice', &
              'long_name', c='PDF of ICNC (in-cloud)')
         CALL channel_halt(substr, status)
         CALL new_attribute(status, 'cloud_cirrus', &
              'CIRRUS_Nice', 'units', c='-' )
         CALL channel_halt(substr, status)

         CALL new_channel_object(status, 'cloud_cirrus', 'CIRRUS_Nice_ML', &
              p4=CIRRUS_Nice_ML, reprid=REPR_CIRRUS_4D_Nice)
         CALL channel_halt(substr, status)
         CALL new_attribute(status, 'cloud_cirrus', 'CIRRUS_Nice_ML', &
              'long_name', c='PDF of ICNC (in-cloud) for ML-CIRRUS')
         CALL channel_halt(substr, status)
         CALL new_attribute(status, 'cloud_cirrus', &
              'CIRRUS_Nice_ML', 'units', c='-' )
         CALL channel_halt(substr, status)

         CALL new_channel_object(status, 'cloud_cirrus', 'CIRRUS_Rice', &
              p4=CIRRUS_Rice, reprid=REPR_CIRRUS_4D_Rice)
         CALL channel_halt(substr, status)
         CALL new_attribute(status, 'cloud_cirrus', 'CIRRUS_Rice', &
              'long_name', c='PDF of IC volume-mean Rice (in-cloud)')
         CALL channel_halt(substr, status)
         CALL new_attribute(status, 'cloud_cirrus', &
              'CIRRUS_Rice', 'units', c='-' )
         CALL channel_halt(substr, status)

         CALL new_channel_object(status, 'cloud_cirrus', 'CIRRUS_RHi_cloud', &
              p4=CIRRUS_RHi_cloud, reprid=REPR_CIRRUS_4D_RHi)
         CALL channel_halt(substr, status)
         CALL new_attribute(status, 'cloud_cirrus', 'CIRRUS_RHi_cloud', &
              'long_name', c='PDF of RH over ice (in-cloud)')
         CALL channel_halt(substr, status)
         CALL new_attribute(status, 'cloud_cirrus', &
              'CIRRUS_RHi_cloud', 'units', c='-' )
         CALL channel_halt(substr, status)

         CALL new_channel_object(status, 'cloud_cirrus', 'CIRRUS_RHi_clear', &
              p4=CIRRUS_RHi_clear, reprid=REPR_CIRRUS_4D_RHi)
         CALL channel_halt(substr, status)
         CALL new_attribute(status, 'cloud_cirrus', 'CIRRUS_RHi_clear', &
              'long_name', c='PDF of RH over ice (clear)')
         CALL channel_halt(substr, status)
         CALL new_attribute(status, 'cloud_cirrus', &
              'CIRRUS_RHi_clear', 'units', c='-' )
         CALL channel_halt(substr, status)

         CALL new_channel_object(status, 'cloud_cirrus', 'CIRRUS_vervel', &
              p4=CIRRUS_vervel, reprid=REPR_CIRRUS_4D_VEL)
         CALL channel_halt(substr, status)
         CALL new_attribute(status, 'cloud_cirrus', 'CIRRUS_vervel', &
              'long_name', c='PDF of vertical velocity')
         CALL channel_halt(substr, status)
         CALL new_attribute(status, 'cloud_cirrus', &
              'CIRRUS_vervel', 'units', c='-' )
         CALL channel_halt(substr, status)

         ! --------------------------------------------------------------------

         CALL new_channel_object(status, modname, 'Nice_preex', p3=Nice_preex)
         CALL channel_halt(substr, status)
         CALL new_attribute(status, modname, 'Nice_preex', &
              'long_name', c='Number concentration of IC (pre-existing)')
         CALL channel_halt(substr, status)
         CALL new_attribute(status, modname, 'Nice_preex', 'units', c='cm-3' )
         CALL channel_halt(substr, status)

         CALL new_channel_object(status, modname, 'Nice_DUdep', p3=Nice_DUdep)
         CALL channel_halt(substr, status)
         CALL new_attribute(status, modname, 'Nice_DUdep', &
              'long_name', c='Number concentration of IC (dust deposition)')
         CALL channel_halt(substr, status)
         CALL new_attribute(status, modname, 'Nice_DUdep', 'units', c='cm-3' )
         CALL channel_halt(substr, status)

         CALL new_channel_object(status, modname, 'Nice_DUimm', p3=Nice_DUimm)
         CALL channel_halt(substr, status)
         CALL new_attribute(status, modname, 'Nice_DUimm', &
              'long_name', c='Number concentration of IC (dust immersion)')
         CALL channel_halt(substr, status)
         CALL new_attribute(status, modname, 'Nice_DUimm', 'units', c='cm-3' )
         CALL channel_halt(substr, status)

         CALL new_channel_object(status, modname, 'Nice_BC', p3=Nice_BC)
         CALL channel_halt(substr, status)
         CALL new_attribute(status, modname, 'Nice_BC', &
              'long_name', c='Number concentration of IC (BC deposition)')
         CALL channel_halt(substr, status)
         CALL new_attribute(status, modname, 'Nice_BC', 'units', c='cm-3' )
         CALL channel_halt(substr, status)

         CALL new_channel_object(status, modname, 'Nice_BCtag', p3=Nice_BCtag)
         CALL channel_halt(substr, status)
         CALL new_attribute(status, modname, 'Nice_BCtag', &
              'long_name', c='Number concentration of IC (BCtag deposition)')
         CALL channel_halt(substr, status)
         CALL new_attribute(status, modname, 'Nice_BCtag', 'units', c='cm-3' )
         CALL channel_halt(substr, status)

         CALL new_channel_object(status, modname, 'Nice_homog', p3=Nice_homog)
         CALL channel_halt(substr, status)
         CALL new_attribute(status, modname, 'Nice_homog', &
              'long_name', c='Number concentration of IC (homog. freezing)')
         CALL channel_halt(substr, status)
         CALL new_attribute(status, modname, 'Nice_homog', 'units', c='cm-3' )
         CALL channel_halt(substr, status)

      ENDIF

   ENDIF

   IF (cloud_param == 6) THEN
      
      CALL INFO_BI("Searching for coupling parameters for CLOUDPARAM = 6!" &
           ,substr)
      CALL get_channel_object(status, 'contrail_gp', 'b_cc', p3=B_cc)
      CALL channel_halt(substr, status)           
      CALL get_channel_object(status, 'contrail_gp', 'potcov', p3=B_co)
      CALL channel_halt(substr, status)           
   
      CALL get_channel_object(status, 'geoloc', 'gboxarea', p2=gboxarea_2d)
      CALL channel_halt(substr, status)           
 
      !aviation inventory
      CALL get_channel_object(status, 'import_grid', 'AEDT1_H2O', p3=ccfh2o_invent)
      CALL channel_halt(substr, status)           
      CALL get_channel_object(status, 'import_grid', 'AEDT2_slantdist', p3=ccfkme_invent)
      CALL channel_halt(substr, status)           
 
   END IF

   CALL end_message_bi(modstr,'COUPLING INITIALIZATION',substr)

 END SUBROUTINE cloud_init_coupling

!===============================================================================

     SUBROUTINE cloud_radiation

       USE messy_main_grid_def_mem_bi, ONLY: jrow, kproma, nlev, nlevp1
       USE messy_main_grid_def_bi,     ONLY: grmass, grvol
       USE messy_main_data_bi,         ONLY: aphm1, apm1, loland_2d, loglac_2d,&
                                             acdnc, slm, seaice,           &
                                             rintop, xvar, xskew, aclc,    &
                                             xlm1, xim1, qm1,              &
                                             vervel => vervel_3d,          &
                                             slf,                          &
                                             tm1, tte_3d,                  &
                                             lcouple
       USE messy_main_data_bi,        ONLY: qte_3d, xlte_3d, xite_3d
       USE messy_main_timer,          ONLY: time_step_len, lstart
       USE messy_main_blather_bi,     ONLY: error_bi
       USE messy_main_tracer_mem_bi,  ONLY: pxtte => qxtte, pxtm1 => qxtm1
       USE messy_cloud_ori,           ONLY: cloud_droplet_nc
       USE messy_cloud_cover,         ONLY: cover_ori
       USE messy_cloud_droplet,       ONLY: cloud_droplet_rotstayn, &
                                            cloud_droplet_menon,    &
                                            cloud_droplet_jones
       USE messy_main_tools,          ONLY: jptlucu1, jptlucu2
       USE messy_cloud_mem,           ONLY: ncdnc &
            , qnuc, qaut, qacc, qfre, qeva, qmel
       USE messy_main_data_bi,        ONLY: betaa_tmp, betab_tmp &
                                          , betass_tmp, nvb_tmp

       INTRINSIC :: ASSOCIATED, NINT

       INTEGER  :: kbdim, ktdia, status
       INTEGER  :: c_type(nproma)
       REAL(dp) :: zfrw(nproma)

       REAL(dp) :: zrho(nproma,nlev)
       REAL(dp) :: aersulf(nproma,nlev), aer_ss(nproma, nlev), &
                   aer_om(nproma,nlev)
       INTEGER  :: i, jk, jl, idt, idt1, idt2

       REAL(dp) :: help
       REAL(dp), POINTER :: temp(:,:) => NULL()

       REAL(DP),DIMENSION(:,:),POINTER:: zbetaa2, zbetab2, zbetass2
       REAL(DP),DIMENSION(:),POINTER:: invb2
       zbetaa2  => betaa_tmp(1:kproma,:,jrow)
       zbetab2  => betab_tmp(1:kproma,:,jrow)
       zbetass2 => betass_tmp(1:kproma,:,jrow)
       invb2    => nvb_tmp(1:kproma,jrow)

       ! reset diagnsotic fields to avoid restart frequency dependent results
       ! (note that the fields are not recalculated in each box in each time
       !  step!); objects are only present for this condition ...
       if (ncdnc > 0) THEN
          !swat
          QNUC(:,:,jrow) = 0.0_dp
          QAUT(:,:,jrow) = 0.0_dp
          QACC(:,:,jrow) = 0.0_dp
          QFRE(:,:,jrow) = 0.0_dp
          QEVA(:,:,jrow) = 0.0_dp
          QMEL(:,:,jrow) = 0.0_dp
          !CDNC
          !CDNC_acc
          !CDNC_burden
          !CDNC_burden_acc
          !cloud_tm1
          !
       end if

       kbdim =  nproma
       ktdia = 1

! standard ECHAM5 cloud droplet number concentration
! calculated only at initialisation

       IF (LSTART) CALL cloud_droplet_nc(kproma, kbdim, nlev, apm1, &
                                         acdnc(:,:,jrow),           &
                                         loland_2d(:,jrow), loglac_2d(:,jrow))

!------------------------------------------------------------------------
! additional for cdnc
       if (l_cdnc_calc) then

!-------------------
! calculate rho = airdensity
         do jk=1,nlev
           zrho(1:kproma,jk) = grmass(1:kproma,jk,jrow) / &
                               grvol (1:kproma,jk,jrow)
         enddo
!-----------------
!   sum up all sulphate aerosol components
!   to determine total sulphate mass mug / m^3
         aersulf(:,:) = 0._dp
         do i=1,sulf_count
           idt = sulf_idt(i)
           do jk=1,nlev
             do jl=1,kproma
! sulphate in mol/mol             
               aersulf(jl,jk) = aersulf(jl,jk)   + &
                                pxtm1(jl,jk,idt) + &
                                pxtte(jl,jk,idt) * time_step_len
             enddo
           enddo
         enddo
! convert mixing ratio to mass concentration
! * fac_sulf : mol/mol -> kg/kg
! / zrho     : kg/kg / kg/m^3 -> kg/m^3
! * 1e9      : kg/m^3 -> mug/m^3
         do jk=1,nlev
           do jl=1,kproma
             aersulf(jl,jk) = MAX(aersulf(jl,jk) * 1.0e9_dp * &
                              fac_sulf / zrho (jl,jk), 0._dp) 
           enddo
         enddo
!--------------------------------
!--------------------------
         if (i_cdnc_calc >= 3) then

!   sum up all seasalt aerosol components
!   to determine total seasalt mass mug / m^3
           aer_ss(:,:) = 0._dp
           do i=1,ss_count
             idt = ss_idt(i)
             do jk=1,nlev
               do jl=1,kproma
! seasalt in mol/mol             
                 aer_ss(jl,jk) = aer_ss(jl,jk)    + &
                                 pxtm1(jl,jk,idt) + &
                                 pxtte(jl,jk,idt) * time_step_len
               enddo
             enddo
           enddo
           ! consider the minimum of Na+ and Cl-
           ! to take NaCl into account
           do i=1,ss2_count
             idt1 = ss2_idt(i,1)
             idt2 = ss2_idt(i,2)
             do jk=1,nlev
               do jl=1,kproma
                 help = MIN(pxtm1(jl,jk,idt1) +            &
                   pxtte(jl,jk,idt1) * time_step_len,      &
                   pxtm1(jl,jk,idt2) + pxtte(jl,jk,idt2) * time_step_len)
! seasalt in mol/mol             
                 aer_ss(jl,jk) = aer_ss(jl,jk) + help
               enddo
             enddo
           enddo
           
! convert mixing ratio to mass concentration
! * fac_ss   : mol/mol -> kg/kg
! / zrho     : kg/kg / kg/m^3 -> kg/m^3
! * 1e9      : kg/m^3 -> mug/m^3
           do jk=1,nlev
             do jl=1,kproma
               aer_ss(jl,jk) = aer_ss(jl,jk) * 1.0e9_dp * &
                               fac_ss / zrho (jl,jk) 
             enddo
           enddo
!-------------
!   sum up all organic matter aerosol components
!   to determine total organic matter mass mug / m^3
           aer_om(:,:) = 0._dp
           do i=1,om_count
             idt = om_idt(i)
             do jk=1,nlev
               do jl=1,kproma
! sulphate in mol/mol             
                 aer_om(jl,jk) = aer_om(jl,jk)    + &
                                 pxtm1(jl,jk,idt) + &
                                 pxtte(jl,jk,idt) * time_step_len
               enddo
             enddo
           enddo
           
! convert mixing ratio to mass concentration
! * fac_om   : mol/mol -> kg/kg
! / zrho     : kg/kg / kg/m^3 -> kg/m^3
! * 1e9      : kg/m^3 -> mug/m^3
           do jk=1,nlev
             do jl=1,kproma
               aer_om(jl,jk) = aer_om(jl,jk) * 1.0e9_dp * &
                               fac_om / zrho (jl,jk) 
             enddo
           enddo
         endif
      END IF
!------------------------------
      IF (L_CDNC_CALC) THEN
         SELECT CASE(i_cdnc_calc)

         CASE(1)
! standard ECHAM5 cloud droplet number concentration
! calculated every time step
            CALL cloud_droplet_nc(kproma, kbdim, nlev, &
                                  apm1, acdnc_1(:,:,jrow),  &
                                  loland_2d(:,jrow), loglac_2d(:,jrow))
         CASE(2)
            CALL cloud_droplet_rotstayn(nlev, kproma, slf(1:kproma,jrow), &
                                        aersulf(1:kproma,:),              &
                                        acdnc_2(1:kproma,:,jrow))

         CASE(3)
            CALL cloud_droplet_jones(nlev, kproma, fac_sulf,                &
                                     aersulf(1:kproma,:),                   &
                                     acdnc_3(1:kproma,:,jrow))

         CASE(4)
            CALL cloud_droplet_menon(nlev, kproma, slf(1:kproma,jrow),     &
                                  aersulf(1:kproma,:), aer_ss(1:kproma,:), &
                                  aer_om(1:kproma,:), acdnc_4(1:kproma,:,jrow))

         CASE(5)
            CALL cloud_activate_ARG
           
         CASE(6)
            CALL cloud_activate_LIN
     
         CASE(7)                    
            CALL cloud_activate_UAF

         CASE(9)
            CALL cloud_droplet_nc(kproma, kbdim, nlev,     &
                                  apm1, acdnc_1(:,:,jrow), &
                                  loland_2d(:,jrow), loglac_2d(:,jrow))

            CALL cloud_droplet_rotstayn(nlev, kproma, slf(1:kproma,jrow), &
                                        aersulf(1:kproma,:),              &
                                        acdnc_2(1:kproma,:,jrow))
        
            CALL cloud_droplet_jones(nlev, kproma, fac_sulf,              &
                                     aersulf(1:kproma,:),                 &
                                     acdnc_3(1:kproma,:,jrow))

            CALL cloud_droplet_menon(nlev, kproma, slf(1:kproma,jrow),        &
                                     aersulf(1:kproma,:), aer_ss(1:kproma,:), &
                                     aer_om(1:kproma,:), acdnc_4(1:kproma,:,jrow))

            CALL cloud_activate_ARG
            
            CALL cloud_activate_LIN
           
            CALL cloud_activate_UAF

         END SELECT
        !-----------------------------------
         aersulf_0(1:kproma,1:nlev,jrow) = aersulf(1:kproma,1:nlev)
      ENDIF

!======================================================
! coupling of alternative CDNCs to ECHAM CDNC field

      IF (ncdnc == 0) THEN
         SELECT CASE (i_cdnc_cpl)

         CASE(1)
            acdnc(1:kproma,1:nlev,jrow) = acdnc_1(1:kproma,1:nlev,jrow)
         CASE(2)
            acdnc(1:kproma,1:nlev,jrow) = acdnc_2(1:kproma,1:nlev,jrow)
         CASE(3)
            acdnc(1:kproma,1:nlev,jrow) = acdnc_3(1:kproma,1:nlev,jrow)
         CASE(4)
            acdnc(1:kproma,1:nlev,jrow) = acdnc_4(1:kproma,1:nlev,jrow)
         END SELECT
      ENDIF
       
!------------------------------------------------------------------

      IF (LSTART) RETURN
      
       ! compiler workaround for xlf95_r
      IF (ALLOCATED(invb))    DEALLOCATE(invb)
      IF (ALLOCATED(zbetaa))  DEALLOCATE(zbetaa)
      IF (ALLOCATED(zbetab))  DEALLOCATE(zbetab)
      IF (ALLOCATED(zbetass)) DEALLOCATE(zbetass)

      ALLOCATE(invb(kbdim))
      ALLOCATE(zbetaa(kbdim,nlev))
      ALLOCATE(zbetab(kbdim,nlev))
      ALLOCATE(zbetass(kbdim,nlev))
      !  bugfix initialise allocated variables
      invb    = 0._dp
      zbetaa  = 0._dp
      zbetab  = 0._dp
      zbetass = 0._dp
      ! 

      c_type(:) = 0
      if (ASSOCIATED(conv_type)) c_type(1:kproma) = &
           NINT(conv_type(1:kproma,jrow))
      zfrw(1:kproma)=(1.0_dp-slm(1:kproma,jrow))*(1.0_dp-seaice(1:kproma,jrow))
      status = 0

      temp => tm1(:,:,jrow)
      CALL COVER_ori(kproma, kbdim, ktdia,   nlev,             nlevp1,        &
                     c_type,                 zfrw,             invb,          &
                     rintop(:,jrow),         aphm1,            apm1,          &
                     qm1(:,:,jrow),          tm1(:,:,jrow),    xlm1(:,:,jrow),&
                     xim1(:,:,jrow),         vervel(:,:,jrow), xvar(:,:,jrow),&
                     xskew(:,:,jrow),        aclc(:,:,jrow),   zbetaa,        &
                     zbetab,                 zbetass,          lcover,        &
                     lcouple,                status,                          &
                     rhc(:,:,jrow) )

!!$#endif

      if (status == 1) THEN
         do jk=1,nlev
            do jl=1,kproma
               if ( ( INT(temp(jl,jk)*1000.) <jptlucu1 .OR.            &
                    INT(temp(jl,jk)*1000.) >jptlucu2) .AND.            &
                    (apm1(jk,jl) >= 1.0_dp) )                          &
                    print*, jk, jl,temp(jl,jk)*1000., tm1(jl,jk,jrow), &
                    tte_3d(jl,jk,jrow)*time_step_len
            enddo
         enddo
         CALL error_bi('lookuperror in cover', modstr)
      ENDIF

      DO jk = 1, nlev
         zbetaa2(1:kproma,jk)  = zbetaa(1:kproma,jk)
         zbetab2(1:kproma,jk)  = zbetab(1:kproma,jk)
         zbetass2(1:kproma,jk) = zbetass(1:kproma,jk)
      ENDDO
      invb2(1:kproma) = REAL(invb(1:kproma),DP)

     END SUBROUTINE cloud_radiation

!===============================================================================

     SUBROUTINE cloud_convec

       USE messy_main_data_bi,         ONLY: ltdiag
       USE messy_main_grid_def_mem_bi, ONLY: jrow, kproma, nlev, nlevp1
       USE messy_main_grid_def_bi,     ONLY: coslat_2d
       USE messy_main_data_bi,   ONLY: xlm1, xim1, tm1, qm1,                 &
                                       um1, vm1,                             &
                                       xtec, xvar, xskew, aclc,              &
                                       aclcac, aclcov, relhum,               &
                                       aprl, aprs, qvi, xlvi, xivi,          &
                                       qte_3d, tte_3d, vervel=> vervel_3d,   &
                                       xlte_3d, xite_3d, aphm1, apm1, app1,  &
                                       geopot_3d, rsfl_2d, ssfl_2d, vo_scb,  &
                                       tke,   slm, glac, acdnc
#ifdef CESM1
       USE messy_main_data_bi,   ONLY: aprflux
#endif

       USE messy_main_timer,      ONLY: time_step_len, lstart
       USE messy_main_blather_bi, ONLY: error_bi, info_bi
       USE messy_cloud_ori,       ONLY: vtmpc1, cloud_ori, lookupoverflow, &
                                        ctaus, ctaul, ctauk
       USE messy_cloud_lohmann07,ONLY: cloud_cdnc, cloud_cdnc_icnc
       USE messy_cloud_lohmann10,ONLY: cloud_cdnc_icnc2 &
                                     , cloud_cdnc_icnc_cc
       USE messy_cloud_kuebbeler14,ONLY: cloud_cdnc_icnc3, vervel_scheme
       USE messy_cloud_mem,        ONLY: ncdnc, idt_cdnc, idt_icnc, nmod, &
                                         vervel_p18
       USE messy_main_constants_mem,  ONLY: ccwmin
       USE messy_main_tracer_mem_bi,  ONLY: pxtte => qxtte, pxtm1 => qxtm1, &
                                            ntrac => ntrac_gp
       USE messy_main_tools,          ONLY: jptlucu1, jptlucu2
!WISO++
       USE messy_main_data_wiso_bi,   ONLY: l_wiso, mwiso, l_wiso_nocloud_dd &
                                          , idiso &
                                          , i_vap, i_liq, i_ice &
                                          , pwisoqm1, pwisoxlm1, pwisoxim1 &
                                          , pwisoqte, pwisoxlte, pwisoxite &
                                          , pwiso3, pwiso2 &
                                          , i_xtec, i_qtec & ! 3D
                                          , i_aprl, i_qvi, i_xlvi, i_xivi &
                                          , i_ssfl, i_rsfl, i_aprs
       USe messy_main_tools_wiso,     ONLY: cwisomin, tnat, i_HHO, i_HDO
!WISO--
       
       IMPLICIT NONE

       INTRINSIC :: ABS, MAX, MIN, ALLOCATED, NINT

       INTEGER  :: ktdia, kbdim, jl, jk, jm
       REAL(dp) :: ztmst

       REAL(dp) :: ztvm1(nproma, nlev), ztp1(nproma, nlev)
       REAL(dp) :: zhmixtau(nproma, nlev)

       REAL(dp) :: zxtec0(nproma, nlev)
       REAL(dp) :: zaclc0(nproma, nlev)
       REAL(dp) :: zaclcac0(nproma, nlev)
       REAL(dp) :: zqte0(nproma, nlev)
       REAL(dp) :: ztte0(nproma, nlev)
       REAL(dp) :: zxlte0(nproma, nlev)
       REAL(dp) :: zxite0(nproma, nlev)
       REAL(dp) :: c_bottom(nproma), zwcape(nproma)
       LOGICAL  :: lo, lo1
       REAL(dp) :: zna(nproma,nlev), nucl(nproma,nlev), nuclm(nproma,nlev,nmod)
!WISO++
       REAL(DP) :: zpwisoqte0(nproma,nlev,mwiso)
       REAL(DP) :: zpwisoxlte0(nproma,nlev,mwiso)
       REAL(DP) :: zpwisoxite0(nproma,nlev,mwiso)
       REAL(DP) :: zdtime
       INTEGER  :: jwiso
       REAL(DP) :: HHOVAP_cloud(nproma,nlev), HDOVAP_cloud(nproma,nlev)
!WISO--
#ifdef MESSYTENDENCY
       REAL(dp) :: save_tte(nproma, nlev), save_qte(nproma, nlev),   &
            save_xlte(nproma, nlev), save_xite(nproma, nlev)
       REAL(dp), DIMENSION(:,:,:), POINTER :: save_xtte => NULL()  
!WISO++
       REAL(DP) :: save_wiso_qte(nproma,nlev,mwiso)
       REAL(DP) :: save_wiso_xlte(nproma,nlev,mwiso)
       REAL(DP) :: save_wiso_xite(nproma,nlev,mwiso)
!WISO--
#endif
       REAL(dp) :: zxtp1(nproma, nlev, ntrac)
       REAL(dp) :: zxtm1(nproma,nlev,ntrac_tot), &
                   zxtte(nproma,nlev,ntrac_tot)
       
       LOGICAL, DIMENSION(nproma,nlev) :: val_psc
       CHARACTER(LEN=32)               :: status_string
       CHARACTER(LEN=70)               :: cirrus_string
       INTEGER                         :: cirrus_status
       CHARACTER(LEN=*), PARAMETER     :: substr='cloud_convec'

       ! Note: These are only set below for cloud_param = 5,
       !       but they are used for all cloud_param. Thus, they
       !       need to be initialised to avoid the model termination.
       cirrus_status = 0
       cirrus_string = ''

       ztmst=time_step_len
       zxtp1 = 0._dp

       IF (ltdiag) THEN
         ! prepare next fields for COND
         pdiga(1:kproma,:,jrow) = pdiga(1:kproma,:,jrow) &
              - tte_3d(1:kproma,:,jrow)
       ENDIF

       kbdim =  nproma
       ktdia = 1
       zhmixtau = 0._dp
       DO jk=1,nlev
         DO jl=1,kproma
           zhmixtau(jl,jk)=ctauk*ABS(vo_scb(jl,jk,jrow)) 
           zhmixtau(jl,jk)=MIN(ctaus,MAX(ctaul,zhmixtau(jl,jk)))
         enddo
       enddo
       ztvm1 = 0._dp     ! um_ak_20091105 initialise
       ztvm1(1:kproma,:) = tm1(1:kproma,:,jrow) *                        &
                          (1.0_dp+vtmpc1*qm1(1:kproma,:,jrow) -          &
                          (xlm1(1:kproma,:,jrow)+xim1(1:kproma,:,jrow)))

       ! Store old tendency values
       IF (USE_PSC) then
         zxtec0(1:kproma,:)    = xtec(1:kproma,:,jrow)
         zaclc0(1:kproma,:)    = aclc(1:kproma,:,jrow)
         ! ADJUST LARGE SCALE CLOUD COVER, BECAUSE xl and xi might have been
         ! modified by other processes
         DO jk=1,nlev
            DO jl=1,kproma
               lo  = (xlm1(jl,jk,jrow) + xlte_3d(jl,jk,jrow)*time_step_len) < ccwmin
               lo1 = (xim1(jl,jk,jrow) + xite_3d(jl,jk,jrow)*time_step_len) < ccwmin
               zaclc0(jl,jk) = MERGE(0.0_dp, zaclc0(jl,jk), lo.AND.lo1)
            END DO
         END DO
         zaclcac0(1:kproma,:)  = aclcac(1:kproma,:,jrow)
         zqte0(1:kproma,:)     = qte_3d(1:kproma,:,jrow)
         ztte0(1:kproma,:)     = tte_3d(1:kproma,:,jrow)
         zxlte0(1:kproma,:)    = xlte_3d(1:kproma,:,jrow)
         zxite0(1:kproma,:)    = xite_3d(1:kproma,:,jrow)
         IF (l_wiso) THEN
            ! rank 3 is nwiso:
            zpwisoqte0(1:kproma,:,:)  = pwisoqte(1:kproma,:,:)
            zpwisoxlte0(1:kproma,:,:) = pwisoxlte(1:kproma,:,:)
            zpwisoxite0(1:kproma,:,:) = pwisoxite(1:kproma,:,:)
         END IF
       ENDIF

       mlwc(:,:,jrow)      = 0.0_dp
       msic(:,:,jrow)      = 0.0_dp
       frain(:,:,jrow)     = 0.0_dp
       fsnow(:,:,jrow)     = 0.0_dp
       frain_no(:,:,jrow)  = 0.0_dp
       fsnow_no(:,:,jrow)  = 0.0_dp
       mratep(:,:,jrow)    = 0.0_dp
       mratesi(:,:,jrow)   = 0.0_dp
       mrevap(:,:,jrow)    = 0.0_dp
       mssubl(:,:,jrow)    = 0.0_dp
       preccover(:,:,jrow) = 0.0_dp
       condens(:,:,jrow)   = 0.0_dp
       mimelt(:,:,jrow)    = 0.0_dp
       misedi(:,:,jrow)    = 0.0_dp
       
       no_lstart: IF (.not. lstart) THEN

#ifdef MESSYTENDENCY
          save_tte(1:kproma,:) = tte_3d(1:kproma,:,jrow)
          save_qte(1:kproma,:) = qte_3d(1:kproma,:,jrow)
          save_xlte(1:kproma,:) = xlte_3d(1:kproma,:,jrow)
          save_xite(1:kproma,:) = xite_3d(1:kproma,:,jrow)
          IF (l_wiso) THEN
             ! rank 3 is 1:mwiso, pwiso* already 1:kprom !
             save_wiso_qte(1:kproma,:,:)  = pwisoqte(1:kproma,:,:)
             save_wiso_xlte(1:kproma,:,:) = pwisoxlte(1:kproma,:,:)
             save_wiso_xite(1:kproma,:,:) = pwisoxite(1:kproma,:,:)
          END IF
#else
          sdtdt_cloud(1:kproma,:,jrow) = tte_3d(1:kproma,:,jrow) 
#endif
    
          select_cloud_param: SELECT CASE (cloud_param)

          CASE(1)

             CALL CLOUD_ORI(kproma, kbdim, ktdia, nlev, nlevp1, ztmst,    &
                  aphm1,              vervel(:,:,jrow),   apm1,   app1,   &
                  acdnc(:,:,jrow),    qm1(:,:,jrow),      tm1(:,:,jrow),  &
                  ztvm1,              xlm1(:,:,jrow),     xim1(:,:,jrow), &
                  xtec(:,:,jrow),     xvar(:,:,jrow),     xskew(:,:,jrow),&
                  pqtec(:,:,jrow),    zbetaa,             zbetab,         &
                  pvdiffp(:,:,jrow),  zhmixtau(:,:),                      &
                  pvmixtau(:,:,jrow), geopot_3d(:,:,jrow),zbetass,        &
                  invb,               aclc(:,:,jrow),                     &
                  aclcac(:,:,jrow),   relhum(:,:,jrow),                   &
                  aclcov(:,jrow),     aprl(:,jrow),       qvi(:,jrow),    &
                  xlvi(:,jrow),       xivi(:,jrow),                       &
                  ssfl_2d(:,jrow),    rsfl_2d(:,jrow),                    &
                  qte_3d(:,:,jrow),   tte_3d(:,:,jrow),                   &
                  xlte_3d(:,:,jrow),  xite_3d(:,:,jrow),                  &
                  aprs(:,jrow),       lcover,                             &
                  !                     added for channel output
                  mlwc(:,:,jrow),     msic(:,:,jrow),                     &
                  frain(:,:,jrow),    fsnow(:,:,jrow),                    &
                  frain_no(:,:,jrow), fsnow_no(:,:,jrow),                 &
                  mratep(:,:,jrow),   mratesi(:,:,jrow),                  &
                  mrevap(:,:,jrow),   mssubl(:,:,jrow),                   &
                  preccover(:,:,jrow),condens(:,:,jrow),                  &
                  mimelt(:,:,jrow),   misedi(:,:,jrow),                   &
                  !WISO++
                  l_wiso, mwiso,                                          &
                  pwisoqm1, pwisoxlm1, pwisoxim1,                         &
                  pwisoqte, pwisoxlte, pwisoxite,                         &
                  pwiso3(i_xtec)%ptr, pwiso3(i_qtec)%ptr,                 &
                  pwiso2(i_aprl)%ptr, pwiso2(i_qvi)%ptr,                  &
                  pwiso2(i_xlvi)%ptr, pwiso2(i_xivi)%ptr,                 &
                  pwiso2(i_ssfl)%ptr, pwiso2(i_rsfl)%ptr,                 &
                  pwiso2(i_aprs)%ptr                                      &
                  !WISO++
                  )

#ifdef MESSYTENDENCY
! NOTE THAT H2OISO / WISO WORK ONLY WITH MESSYTENDENCY ENABLED
             IF (l_wiso) THEN
                IF (l_wiso_nocloud_dd) THEN
                   ! for sensitivity study without the effect of cloud on deltaD
                   HHOVAP_cloud(1:kproma,:) = pwisoqm1(1:kproma,:,i_HHO) + &
                        pwisoqte(1:kproma,:,i_HHO) * ztmst
                   HDOVAP_cloud(1:kproma,:) = pwisoqm1(1:kproma,:,i_HDO) + &
                        pwisoqte(1:kproma,:,i_HDO) * ztmst
                   ! keep the delta value of vapour as the delta value before
                   ! cloud
                   DO jk = 1, nlev
                      DO jl = 1, kproma
                         IF(HHOVAP_cloud(jl,jk) > cwisomin) THEN
                            pwisoqte(jl,jk,i_HDO) = &
                                 save_wiso_qte(jl,jk,i_HHO) * &
                                 ( HDOVAP_cloud(jl,jk)/HHOVAP_cloud(jl,jk) )
                         ELSE
                            pwisoqte(jl,jk,i_HDO) = &
                                 save_wiso_qte(jl,jk,i_HHO) * tnat(i_HDO)
                         ENDIF
                      ENDDO
                   ENDDO
                ENDIF
             END IF
#endif             

          CASE(2)
           
             c_bottom(:) = 0._dp
             IF (ASSOCIATED(conv_bot)) c_bottom(:) = conv_bot(:,jrow)
             zwcape(:)   = 0._dp
             IF (ASSOCIATED(PCAPE)) THEN
                DO jl=1,kproma
                   IF (PCAPE(jl,jrow) > 0._dp) &
                        zwcape(jl) = 0.5_dp*SQRT(PCAPE(jl,jrow)/0.05_dp)
                END Do
             ENDIF
             zna(1:kproma,1:nlev) = 0._dp
             SELECT CASE(ncdnc)
             CASE(1) 
                nucl(1:kproma,1:nlev) = acdnc_6(1:kproma,1:nlev,jrow)
                zna(1:kproma,1:nlev)  = np_pot_activ(1:kproma,1:nlev,jrow)
             CASE(2) 
                nucl(1:kproma,1:nlev) = acdnc_5(1:kproma,1:nlev,jrow)
             END SELECT

#ifdef MESSYTENDENCY
             ALLOCATE(save_xtte(kproma,nlev,2))
             save_xtte(:,:,1) = pxtte(1:kproma,:,idt_cdnc)
             save_xtte(:,:,2) = pxtte(1:kproma,:,idt_icnc)
#endif
             
             CALL cloud_cdnc (kproma, kbdim, ktdia, nlev, nlevp1, ztmst,      &
!---Included for in-cloud scavenging (Philip Stier, 28/03/01):----------
                  ntrac,  jrow,                                               &
!---End Included for scavenging-----------------------------------------
!-----------------------------------------------------------------------
! - INPUT  2D .
                  aphm1,              vervel(:,:,jrow),   apm1,  app1,    &
                  acdnc(:,:,jrow),    qm1(:,:,jrow),      tm1(:,:,jrow),  &
                  ztvm1,              xlm1(:,:,jrow),     xim1(:,:,jrow), &
                  xtec(:,:,jrow),     xvar(:,:,jrow),     xskew(:,:,jrow),&
                  pqtec(:,:,jrow),    zbetaa,             zbetab,         &
                  pvdiffp(:,:,jrow),  zhmixtau,                           &
                  pvmixtau(:,:,jrow), geopot_3d(:,:,jrow),zbetass,        &
                  !---Included for in-cloud scavenging (Philip Stier, 28/03/01):----------
                  pxtm1,                                                  &
!---End Included for scavenging-----------------------------------------
!--- Included for prognostic CDNC/IC scheme ----------------------------
                  tke(:,:,jrow),       c_bottom(:),       zwcape(:),      &
!--- End included for CDNC/IC scheme -----------------------------------
! - INPUT  1D .
                  invb,                aclc(:,:,jrow),                    &
                  aclcac(:,:,jrow),    relhum(:,:,jrow),                  &
                  aclcov(:,jrow),      aprl(:,jrow),      qvi(:,jrow),    & 
                  xlvi(:,jrow),        xivi(:,jrow),                      &
                  ssfl_2d(:,jrow),     rsfl_2d(:,jrow),                   &
                  qte_3d(:,:,jrow),    tte_3d(:,:,jrow),                  &
                  xlte_3d(:,:,jrow),   xite_3d(:,:,jrow),                 &
!---Included for in-cloud scavenging (Philip Stier, 28/03/01):----------
                  pxtte,                                                  &
!---End Included for scavenging-----------------------------------------
! - INPUT/OUTPUT 1D .
                  aprs,                                                   &
                  ! mz_ht_20070629+
                  status_string,     lcover,                              &
                  slm,               glac,                nucl,           &
                  zna,                                                    &
                  ! added for channel output
                  mlwc(:,:,jrow),     msic(:,:,jrow),                     &
                  frain(:,:,jrow),    fsnow(:,:,jrow),                    &
                  frain_no(:,:,jrow), fsnow_no(:,:,jrow),                 &
                  mratep(:,:,jrow),   mratesi(:,:,jrow),                  &
                  mrevap(:,:,jrow),   mssubl(:,:,jrow),                   &
                  ! mz_ak_20051221 added condens
                  preccover(:,:,jrow),condens(:,:,jrow),                  &
                  mimelt(:,:,jrow),   misedi(:,:,jrow)      )
                      
             if (lookupoverflow) CALL INFO_BI(status_string, modstr)

#ifdef MESSYTENDENCY
             save_xtte(1:kproma,:,1) = &
                  pxtte(1:kproma,:,idt_cdnc) - save_xtte(1:kproma,:,1)
             save_xtte(1:kproma,:,2) = &
                  pxtte(1:kproma,:,idt_icnc) - save_xtte(1:kproma,:,2)
             !
             CALL mtend_add_l(my_handle, idt_cdnc,  save_xtte(:,:,1) &
                  , l_add = .FALSE.)
             CALL mtend_add_l(my_handle, idt_icnc,  save_xtte(:,:,2) &
                  , l_add = .FALSE.)
             !
             DEALLOCATE(save_xtte); NULLIFY(save_xtte)
#endif
             
        CASE(3)

           c_bottom(:) = 0._dp
           IF (ASSOCIATED(conv_bot)) c_bottom(:) = conv_bot(:,jrow)
           zwcape(:)   = 0._dp
           IF (ASSOCIATED(PCAPE)) THEN
              DO jl=1,kproma
                 IF (PCAPE(jl,jrow) > 0._dp) &
                      zwcape(jl) = 0.5_dp*SQRT(PCAPE(jl,jrow)/0.05_dp)
              END Do
           ENDIF
           zna(1:kproma,1:nlev) = 0._dp
           SELECT CASE(ncdnc)
           CASE(1) 
              nucl(1:kproma,1:nlev) = acdnc_6(1:kproma,1:nlev,jrow)
              zna(1:kproma,1:nlev)  = np_pot_activ(1:kproma,1:nlev,jrow)
           CASE(2) 
              nucl(1:kproma,1:nlev) = acdnc_5(1:kproma,1:nlev,jrow)
           END SELECT

#ifdef MESSYTENDENCY
             ALLOCATE(save_xtte(kproma,nlev,2))
             save_xtte(:,:,1) = pxtte(1:kproma,:,idt_cdnc)
             save_xtte(:,:,2) = pxtte(1:kproma,:,idt_icnc)
#endif

           CALL cloud_cdnc_icnc(kproma, kbdim, ktdia, nlev, nlevp1, ztmst,    &
!---Included for in-cloud scavenging (Philip Stier, 28/03/01):----------
                ntrac,  jrow,                                           &
!---End Included for scavenging-----------------------------------------
!-----------------------------------------------------------------------
! - INPUT  2D .
                aphm1,              vervel(:,:,jrow),   apm1,  app1,    &
                acdnc(:,:,jrow),    qm1(:,:,jrow),      tm1(:,:,jrow),  &
                ztvm1,              xlm1(:,:,jrow),     xim1(:,:,jrow), &
                xtec(:,:,jrow),     xvar(:,:,jrow),     xskew(:,:,jrow),&
                pqtec(:,:,jrow),    zbetaa,             zbetab,         &
                pvdiffp(:,:,jrow),  zhmixtau,                           &
                pvmixtau(:,:,jrow), geopot_3d(:,:,jrow),zbetass,        &
!---Included for in-cloud scavenging (Philip Stier, 28/03/01):----------
                pxtm1,                                                  &
!---End Included for scavenging-----------------------------------------
!--- Included for prognostic CDNC/IC scheme ----------------------------
                tke(:,:,jrow),       c_bottom(:),       zwcape(:),      &
                pxtecl(:,:,jrow),    pxteci(:,:,jrow),                  &
!--- End included for CDNC/IC scheme -----------------------------------
! - INPUT  1D .
                invb,                aclc(:,:,jrow),                    &
                aclcac(:,:,jrow),    relhum(:,:,jrow),                  &
                aclcov(:,jrow),      aprl(:,jrow),      qvi(:,jrow),    & 
                xlvi(:,jrow),        xivi(:,jrow),                      &
                ssfl_2d(:,jrow),     rsfl_2d(:,jrow),                   &
                qte_3d(:,:,jrow),    tte_3d(:,:,jrow),                  &
                xlte_3d(:,:,jrow),   xite_3d(:,:,jrow),                 &
!---Included for in-cloud scavenging (Philip Stier, 28/03/01):----------
                pxtte,                                                  &
!---End Included for scavenging-----------------------------------------
! - INPUT/OUTPUT 1D .
                aprs,                                                   &
                status_string,     lcover,                              &
                slm,               glac,                nucl,           &
                zna,                                                    &
                ! added for channel output
                mlwc(:,:,jrow),     msic(:,:,jrow),                     &
                frain(:,:,jrow),    fsnow(:,:,jrow),                    &
                frain_no(:,:,jrow), fsnow_no(:,:,jrow),                 &
                mratep(:,:,jrow),   mratesi(:,:,jrow),                  &
                mrevap(:,:,jrow),   mssubl(:,:,jrow),                   &
                preccover(:,:,jrow),condens(:,:,jrow),                  &
                mimelt(:,:,jrow),   misedi(:,:,jrow),                   &
                sigma(:), wetradius(:,:,:,jrow) )

           if (lookupoverflow) CALL INFO_BI( status_string, modstr)

#ifdef MESSYTENDENCY
             save_xtte(1:kproma,:,1) = &
                  pxtte(1:kproma,:,idt_cdnc) - save_xtte(1:kproma,:,1)
             save_xtte(1:kproma,:,2) = &
                  pxtte(1:kproma,:,idt_icnc) - save_xtte(1:kproma,:,2)
             !
             CALL mtend_add_l(my_handle, idt_cdnc,  save_xtte(:,:,1) &
                  , l_add = .FALSE.)
             CALL mtend_add_l(my_handle, idt_icnc,  save_xtte(:,:,2) &
                  , l_add = .FALSE.)
             !
             DEALLOCATE(save_xtte); NULLIFY(save_xtte)
#endif

        CASE(4, 5)
           ! Handle cases 4 and 5 together
#ifndef MESSYTENDENCY
           zxtp1(1:kproma,:,:) = &
                pxtm1(1:kproma,:,:) + pxtte(1:kproma,:,:) * ztmst
#else
           call mtend_get_start_l(mtend_id_tracer, v0t = zxtp1)
#endif
           zxtm1(1:kproma,:,1) = pxtm1(1:kproma,:,idt_cdnc)
           zxtm1(1:kproma,:,2) = pxtm1(1:kproma,:,idt_icnc)
           zxtte(1:kproma,:,1) = pxtte(1:kproma,:,idt_cdnc)
           zxtte(1:kproma,:,2) = pxtte(1:kproma,:,idt_icnc)

#ifdef MESSYTENDENCY
             ALLOCATE(save_xtte(kproma,nlev,2))
             save_xtte(:,:,1) = pxtte(1:kproma,:,idt_cdnc)
             save_xtte(:,:,2) = pxtte(1:kproma,:,idt_icnc)
#endif

           ! update activation calculation
           c_bottom(:) = 0._dp
           IF (ASSOCIATED(conv_bot)) c_bottom(:) = conv_bot(:,jrow)
           zwcape(:)   = 0._dp
           IF (ASSOCIATED(PCAPE)) THEN
             DO jl=1,kproma
               IF (PCAPE(jl,jrow) > 0._dp) &
                 zwcape(jl) = 0.5_dp*SQRT(PCAPE(jl,jrow)/0.05_dp)
             END Do
           ENDIF

           nuclm(1:kproma,1:nlev,:) = 0._dp
           SELECT CASE(ncdnc)
           CASE(1) 
              CALL cloud_activate_LIN
              nucl(1:kproma,1:nlev) = acdnc_6(1:kproma,1:nlev,jrow)
              zna(1:kproma,1:nlev)  = np_pot_activ(1:kproma,1:nlev,jrow)
              do jm=1,ksol
                 nuclm(1:kproma,1:nlev,sol_idx(jm)) = &
                      lin_num(1:kproma,1:nlev,jm,jrow)
              end do
           CASE(2) 
              CALL cloud_activate_ARG
              nucl(1:kproma,1:nlev) = acdnc_5(1:kproma,1:nlev,jrow)
              do jm=1,ksol
                 nuclm(1:kproma,1:nlev,sol_idx(jm)) = &
                      arg_num(1:kproma,1:nlev,jm,jrow)
              end do
           CASE(3)
              CALL cloud_activate_UAF
              nucl(1:kproma,1:nlev) = acdnc_7(1:kproma,1:nlev,jrow)
              do jm=1,ktot
                 nuclm(1:kproma,1:nlev,sol_idx(jm)) = &
                      uaf_num(1:kproma,1:nlev,jm,jrow)
              end do
           END SELECT

           CALL prepare_freezing(kproma, nlev, jrow, zxtp1(1:kproma,:,:) &
                , nuclm(1:kproma,:,:))

           IF (cloud_param == 5 .AND. vervel_scheme == 4) &
                CALL calc_vervel_penner18(kproma, nlev, jrow &
                , vervel_p18(1:kproma,:,jrow))

           if (cloud_param == 4) &
                CALL cloud_cdnc_icnc2(kproma, kproma, ktdia, nlev, nlevp1 &
                , ztmst  &
!---Included for in-cloud scavenging (Philip Stier, 28/03/01):----------
!               , ntrac,  jrow                                                &
                , 2,  jrow                                                    &
!---End Included for scavenging-----------------------------------------
                , aphm1(1:kproma,:), vervel(1:kproma,:,jrow)                  &
                , apm1(1:kproma,:),  app1(1:kproma,:), acdnc(1:kproma,:,jrow) &
                , qm1(1:kproma,:,jrow), tm1(1:kproma,:,jrow), ztvm1(1:kproma,:)      &
                , xlm1(1:kproma,:,jrow), xim1(1:kproma,:,jrow), xtec(1:kproma,:,jrow)  &
                , xvar(1:kproma,:,jrow), xskew(1:kproma,:,jrow)               &
                , pqtec(1:kproma,:,jrow)                                      &
                , zbetaa(1:kproma,:),       zbetab(1:kproma,:)                &
                , pvdiffp(1:kproma,:,jrow), zhmixtau(1:kproma,:)              &
                , pvmixtau(1:kproma,:,jrow)                                   &
                , geopot_3d(1:kproma,:,jrow), zbetass(1:kproma,:)             &
!---Included for in-cloud scavenging (Philip Stier, 28/03/01):----------
                , zxtm1(1:kproma,:,1:2)                                       &
!---End Included for scavenging-----------------------------------------
!---Included for prognostic CDNC/IC scheme (Philip Stier, 30/11/2002)---
                , tke(1:kproma,:,jrow), c_bottom(1:kproma)                    &
                , zwcape(1:kproma)                                            &
                , pxtecl(1:kproma,:,jrow),  pxteci(1:kproma,:,jrow)           & 
!---Included for prognostic CDNC/IC scheme (Ulrike Lohmann, 11/02/2007)---
               , pxtecnl(1:kproma,:,jrow), pxtecni(1:kproma,:,jrow)           &
!--- End included for CDNC/IC scheme -----------------------------------
               , invb                                                     &
               , aclc(1:kproma,:,jrow),    aclcac(1:kproma,:,jrow)        &
               , relhum(1:kproma,:,jrow)                                  &
               , w_sub(1:kproma,:,jrow), w_ave(1:kproma,:,jrow)           &
               , w_grid(1:kproma,:,jrow)                                  &
               , sigwBN(1:kproma,:,jrow),ndropBN(1:kproma,:,jrow)         &
               , dsulfBN(1:kproma,:,jrow)                                 &
               , ndust_aiBN(1:kproma,:,jrow), ddust_aiBN(1:kproma,:,jrow) &
               , ndust_ciBN(1:kproma,:,jrow), ddust_ciBN(1:kproma,:,jrow) &
               , norgBN(1:kproma,:,jrow),     dorgBN(1:kproma,:,jrow)     &
               , nsootBN(1:kproma,:,jrow),    dsootBN(1:kproma,:,jrow)    &   
               , nsolo_ksBN(1:kproma,:,jrow), dsolo_ksBN(1:kproma,:,jrow) &
               , nsolo_asBN(1:kproma,:,jrow), dsolo_asBN(1:kproma,:,jrow) &
               , nsolo_csBN(1:kproma,:,jrow), dsolo_csBN(1:kproma,:,jrow) &
               , smaxice_cirrusBN(1:kproma,:,jrow), sc_ice_cirrusBN(1:kproma,:,jrow)  &
               , nlim_cirrusBN(1:kproma,:,jrow), nhet_cirrusBN(1:kproma,:,jrow)       &
               , nice_cirrusBN(1:kproma,:,jrow)                                       &
               , dice_cirrusBN(1:kproma,:,jrow), sigwpre_cirrusBN(1:kproma,:,jrow)    &
               , smaxice_immBN(1:kproma,:,jrow), sc_ice_immBN(1:kproma,:,jrow)        &
               , nice_immBN(1:kproma,:,jrow)                                          &
               , dice_immBN(1:kproma,:,jrow), sigwpre_immBN(1:kproma,:,jrow)          &
               , newIC_cirri(1:kproma,:,jrow), newICR_cirri(1:kproma,:,jrow)          &
               , newIC_cnt_therm(1:kproma,:,jrow), newIC_imm(1:kproma,:,jrow)         &
               , newIC_mix(1:kproma,:,jrow), newIC_mix_fin(1:kproma,:,jrow)           &
               , aclcov(1:kproma,jrow),    aprl(1:kproma,jrow)            &
               , qvi(1:kproma,jrow)                                       &
               , xlvi(1:kproma,jrow),      xivi(1:kproma,jrow)            &
               , ssfl_2d(1:kproma,jrow),   rsfl_2d(1:kproma,jrow)         &
               , qte_3d(1:kproma,:,jrow),  tte_3d(1:kproma,:,jrow)        &
               , xlte_3d(1:kproma,:,jrow), xite_3d(1:kproma,:,jrow)       &
!---Included for in-cloud scavenging (Philip Stier, 28/03/01):----------
               , zxtte(1:kproma,:,1:2)                                    &
!---End Included for scavenging-----------------------------------------
               , aprs(1:kproma,jrow)                                      &
               , status_string,     lcover                                &
               , slm(1:kproma,jrow),        glac(1:kproma,jrow)           &
               , nucl(1:kproma,:)                                         &
               ! added for channel output
               , mlwc(1:kproma,:,jrow),     msic(1:kproma,:,jrow)         &
               , frain(1:kproma,:,jrow),    fsnow(1:kproma,:,jrow)        &
               , frain_no(1:kproma,:,jrow), fsnow_no(1:kproma,:,jrow)     &
               , mratep(1:kproma,:,jrow),   mratesi(1:kproma,:,jrow)      &
               , mrevap(1:kproma,:,jrow),   mssubl(1:kproma,:,jrow)       &
               , preccover(1:kproma,:,jrow),condens(1:kproma,:,jrow)      &
               , mimelt(1:kproma,:,jrow),   misedi(1:kproma,:,jrow)       &
               ) 

           IF (cloud_param == 5) &
                                  ! indices
              CALL cloud_cdnc_icnc3(kproma, kproma, ktdia, nlev, nlevp1,    &
                 ! timestep, # of CDNC/ICNC tracers, row, cos(lat)
                 ztmst, 2, jrow, coslat_2d(1:kproma,jrow),                  &
                 ! variables for orographic waves
                 w_gwd_kpr(1:kproma,:,jrow), mask_orogw(1:kproma,jrow),     &
                 l_z(1:kproma,jrow), ampl_gwd(1:kproma,:,jrow),             &
                 ! pressure, vertical velocity
                 aphm1(1:kproma,:), vervel(1:kproma,:,jrow),                &
                 ! vertical velocity from Penner et al. 2018, pressure (full)
                 vervel_p18(1:kproma,:,jrow), apm1(1:kproma,:),             &
                 ! pressure (mid), ECHAM5 CDNC, specific humidity
                 app1(1:kproma,:), acdnc(1:kproma,:,jrow), qm1(1:kproma,:,jrow), &
                 ! temperature, virtual temperature, LWC
                 tm1(1:kproma,:,jrow), ztvm1(1:kproma,:), xlm1(1:kproma,:,jrow),     &
                 ! IWC, detrained convective water/ice
                 xim1(1:kproma,:,jrow), xtec(1:kproma,:,jrow),              &
                 ! variance of total water, skewness of total water
                 xvar(1:kproma,:,jrow), xskew(1:kproma,:,jrow),             &
                 ! convective detrained humidity, beta distrib minimum a
                 pqtec(1:kproma,:,jrow), zbetaa(1:kproma,:),                &
                 ! beta distrib. maximum b, rate of change of q due to vdiff
                 zbetab(1:kproma,:),pvdiffp(1:kproma,:,jrow),               &
                 ! inverse timescale for horizontal and vertical turbulence
                 zhmixtau(1:kproma,:), pvmixtau(1:kproma,:,jrow),           &
                 ! geopotential, ???
                 geopot_3d(1:kproma,:,jrow), zbetass(1:kproma,:),           &
                 ! CDNC/ICNC tracers, turbulent kinetic energy
                 zxtm1(1:kproma,:,1:2), tke(1:kproma,:,jrow),              &
                 ! convective cloud base, convective detrained water
                 c_bottom(1:kproma), pxtecl(1:kproma,:,jrow),               &
                 ! convective detrained ice, convective detrained water number
                 pxteci(1:kproma,:,jrow), pxtecnl(1:kproma,:,jrow),         &
                 ! convective detrained ice number, ???, LS cloud cover
                 pxtecni(1:kproma,:,jrow), invb, aclc(1:kproma,:,jrow),     &
                 ! LS cloud cover accumulated, relative humidity (diag)
                 aclcac(1:kproma,:,jrow), relhum(1:kproma,:,jrow),          &
                 ! total cloud cover, accumulated LS rain at the surface
                 aclcov(1:kproma,jrow), aprl(1:kproma,jrow),                &
                 ! vertically integrated water vapour, LWP
                 qvi(1:kproma,jrow), xlvi(1:kproma,jrow),                   &
                 ! IWP, LS surface-level snow and rain rate
                 xivi(1:kproma,jrow), ssfl_2d(1:kproma,jrow),               &
                 rsfl_2d(1:kproma,jrow),                                    &
                 ! specific humidity and temperature tendency
                 qte_3d(1:kproma,:,jrow), tte_3d(1:kproma,:,jrow),          &
                 ! cloud water and cloud ice tendency
                 xlte_3d(1:kproma,:,jrow), xite_3d(1:kproma,:,jrow),        &
                 ! CDNC/ICNC tendency, accumulated LS snow at the surface
                 zxtte(1:kproma,:,1:2), aprs(1:kproma,jrow),              &
                 ! lookupoverflow, cirrus status, cirrus error,
                 status_string, cirrus_status, cirrus_string,               &
                 ! switch for cover calc, land-sea mask
                 lcover, slm(1:kproma,jrow),  &
                 ! fraction of glacier cover, CDNC from activation scheme,
                 glac(1:kproma,jrow), nucl(1:kproma,:),                     &
                 ! large scale LWC and IWC
                 mlwc(1:kproma,:,jrow), msic(1:kproma,:,jrow),              &
                 ! large scale rain and snow flux
                 frain(1:kproma,:,jrow), fsnow(1:kproma,:,jrow),            &
                 ! large scale rain and snow flux without cloud production
                 frain_no(1:kproma,:,jrow), fsnow_no(1:kproma,:,jrow),      &
                 ! large scale rain and snow formation inside cloud
                 mratep(1:kproma,:,jrow), mratesi(1:kproma,:,jrow),         &
                 ! rain evaporation and snow sublimation rate
                 mrevap(1:kproma,:,jrow), mssubl(1:kproma,:,jrow),          &
                 ! LS precipitation cloud cover, in-cloud condensate
                 preccover(1:kproma,:,jrow), condens(1:kproma,:,jrow),      &
                 ! LS frozen precipitation melting, LS ice sedimentation
                 mimelt(1:kproma,:,jrow), misedi(1:kproma,:,jrow))

           if (lookupoverflow) CALL INFO_BI(status_string, modstr)
           if (cirrus_status /= 0) CALL ERROR_BI(cirrus_string, substr)

           pxtte(1:kproma,:,idt_cdnc) = zxtte(1:kproma,:,1)
           pxtte(1:kproma,:,idt_icnc) = zxtte(1:kproma,:,2)
           
#ifdef MESSYTENDENCY
           save_xtte(:,:,1) = pxtte(1:kproma,:,idt_cdnc) - &
                save_xtte(:,:,1)
           save_xtte(:,:,2) = pxtte(1:kproma,:,idt_icnc) - &
                save_xtte(:,:,2)

           CALL mtend_add_l(my_handle, idt_cdnc, save_xtte(1:kproma,:,1) &
                , l_add = .FALSE.)
           CALL mtend_add_l(my_handle, idt_icnc, save_xtte(1:kproma,:,2) &
                , l_add = .FALSE.)

           DEALLOCATE(save_xtte); NULLIFY(save_xtte)
#endif

        CASE(6)
#ifndef MESSYTENDENCY
           zxtp1(1:kproma,:,:) = &
                pxtm1(1:kproma,:,:) + pxtte(1:kproma,:,:) * ztmst
#else
           call mtend_get_start_l(mtend_id_tracer, v0t = zxtp1)
#endif
           zxtm1(1:kproma,:,1) = pxtm1(1:kproma,:,idt_cdnc)
           zxtm1(1:kproma,:,2) = pxtm1(1:kproma,:,idt_icnc)
           zxtte(1:kproma,:,1) = pxtte(1:kproma,:,idt_cdnc)
           zxtte(1:kproma,:,2) = pxtte(1:kproma,:,idt_icnc)

           DO jm=1, ntrac_cl
              zxtm1(1:kproma,:,jm+2) = pxtm1(1:kproma,:,idt_cl(jm))
              zxtte(1:kproma,:,jm+2) = pxtte(1:kproma,:,idt_cl(jm))
           END DO

#ifdef MESSYTENDENCY
             ALLOCATE(save_xtte(kproma,nlev,ntrac_tot))
             save_xtte(:,:,1) = pxtte(1:kproma,:,idt_cdnc)
             save_xtte(:,:,2) = pxtte(1:kproma,:,idt_icnc)
             DO jm=1, ntrac_cl
                save_xtte(:,:,jm+2) = pxtte(1:kproma,:,idt_cl(jm))
             END DO
#endif
             
           ! update activation calculation
           c_bottom(:) = 0._dp
           IF (ASSOCIATED(conv_bot)) c_bottom(:) = conv_bot(:,jrow)
           zwcape(:)   = 0._dp
           IF (ASSOCIATED(PCAPE)) THEN
              DO jl=1,kproma
                 IF (PCAPE(jl,jrow) > 0._dp) &
                      zwcape(jl) = 0.5_dp*SQRT(PCAPE(jl,jrow)/0.05_dp)
              END DO
           ENDIF

           SELECT CASE(ncdnc)
           CASE(1) 
              CALL cloud_activate_LIN
              nucl(1:kproma,1:nlev) = acdnc_6(1:kproma,1:nlev,jrow)
              zna(1:kproma,1:nlev)  = np_pot_activ(1:kproma,1:nlev,jrow)
           CASE(2) 
              CALL cloud_activate_ARG
              nucl(1:kproma,1:nlev) = acdnc_5(1:kproma,1:nlev,jrow)
           CASE(3)
              CALL cloud_activate_UAF
              nucl(1:kproma,1:nlev) = acdnc_7(1:kproma,1:nlev,jrow)
           END SELECT
           CALL prepare_freezing(kproma, nlev, jrow, zxtp1(1:kproma,:,:) &
                , nucl(1:kproma,:))

           CALL cloud_cdnc_icnc_cc(kproma, kproma, ktdia, nlev, nlevp1, ztmst  &
!---Included for in-cloud scavenging (Philip Stier, 28/03/01):----------
               , ntrac_tot,  jrow                                            &
               , um1(1:kproma,:,jrow), vm1(1:kproma,:,jrow)                  &
!---End Included for scavenging-----------------------------------------
               , aphm1(1:kproma,:), vervel(1:kproma,:,jrow)                  &
               , apm1(1:kproma,:),  app1(1:kproma,:), acdnc(1:kproma,:,jrow) &
               , qm1(1:kproma,:,jrow), tm1(1:kproma,:,jrow), ztvm1(1:kproma,:)      &
               , xlm1(1:kproma,:,jrow),  xim1(1:kproma,:,jrow), xtec(1:kproma,:,jrow)  &
               , xvar(1:kproma,:,jrow), xskew(1:kproma,:,jrow)               &
               , pqtec(1:kproma,:,jrow)                                      &
               , zbetaa(1:kproma,:),       zbetab(1:kproma,:)                &
               , pvdiffp(1:kproma,:,jrow), zhmixtau(1:kproma,:)              &
               , pvmixtau(1:kproma,:,jrow)  &
               , geopot_3d(1:kproma,:,jrow), zbetass(1:kproma,:)             &
!---Included for in-cloud scavenging (Philip Stier, 28/03/01):----------
               , zxtm1(1:kproma,:,:)                                         &
!---End Included for scavenging-----------------------------------------
!---Included for prognostic CDNC/IC scheme (Philip Stier, 30/11/2002)---
               , tke(1:kproma,:,jrow),     c_bottom(1:kproma),   zwcape(1:kproma)    &
               , pxtecl(1:kproma,:,jrow),  pxteci(1:kproma,:,jrow)                   & 
!---Included for prognostic CDNC/IC scheme (Ulrike Lohmann, 11/02/2007)---
               , pxtecnl(1:kproma,:,jrow), pxtecni(1:kproma,:,jrow)                  &
!--- End included for CDNC/IC scheme -----------------------------------

               , invb                                                     &
               , aclc(1:kproma,:,jrow),    aclcac(1:kproma,:,jrow)        &
               , relhum(1:kproma,:,jrow)                                  &
               , aclcov(1:kproma,jrow),    aprl(1:kproma,jrow)            &
               , qvi(1:kproma,jrow)                                       &
               , xlvi(1:kproma,jrow),      xivi(1:kproma,jrow)            &
               , ssfl_2d(1:kproma,jrow),   rsfl_2d(1:kproma,jrow)         &
               , qte_3d(1:kproma,:,jrow),  tte_3d(1:kproma,:,jrow)        &
               , xlte_3d(1:kproma,:,jrow), xite_3d(1:kproma,:,jrow)       &
!---Included for in-cloud scavenging (Philip Stier, 28/03/01):----------
               , zxtte(1:kproma,:,:)                                      &
!---End Included for scavenging-----------------------------------------
               , aprs(1:kproma,jrow)                                      &
               , status_string,     lcover                                &
               , slm(1:kproma,jrow),        glac(1:kproma,jrow)           &
               , nucl(1:kproma,:)                                         &
               ! added for channel output
               , mlwc(1:kproma,:,jrow),     msic(1:kproma,:,jrow)         &
               , frain(1:kproma,:,jrow),    fsnow(1:kproma,:,jrow)        &
               , frain_no(1:kproma,:,jrow), fsnow_no(1:kproma,:,jrow)     &
               , mratep(1:kproma,:,jrow),   mratesi(1:kproma,:,jrow)      &
               , mrevap(1:kproma,:,jrow),   mssubl(1:kproma,:,jrow)       &
               , preccover(1:kproma,:,jrow),condens(1:kproma,:,jrow)      &
               , mimelt(1:kproma,:,jrow),   misedi(1:kproma,:,jrow)       &
               ) 

           if (lookupoverflow) CALL INFO_BI( status_string, modstr)

! XXX WARNING, M1 values are not supposed to be changed
           pxtm1(1:kproma,:,idt_cl(1)) = zxtm1(1:kproma,:,3)
           pxtm1(1:kproma,:,idt_cl(2)) = zxtm1(1:kproma,:,4)
           
           ! COPY ABSOLUTE TENDENCIES BACK TO TOTAL TENDENCIES
           pxtte(1:kproma,:,idt_cdnc) = zxtte(1:kproma,:,1)
           pxtte(1:kproma,:,idt_icnc) = zxtte(1:kproma,:,2)
           DO jm=1, ntrac_cl
              pxtte(1:kproma,:,idt_cl(jm)) = zxtte(1:kproma,:,jm+2)
           END DO

#ifdef MESSYTENDENCY
           save_xtte(:,:,1) = pxtte(1:kproma,:,idt_cdnc) - &
                save_xtte(:,:,1)
           save_xtte(:,:,2) = pxtte(1:kproma,:,idt_icnc) - &
                save_xtte(:,:,2)

           CALL mtend_add_l(my_handle, idt_cdnc, save_xtte(1:kproma,:,1) &
                , l_add = .FALSE.)
           CALL mtend_add_l(my_handle, idt_icnc, save_xtte(1:kproma,:,2) &
                , l_add = .FALSE.)

           DO jm=1, ntrac_cl
              save_xtte(:,:,jm+2) = pxtte(1:kproma,:,idt_cl(jm)) &
                   - save_xtte(:,:,jm+2)
              CALL mtend_add_l(my_handle, idt_cl(jm) &
                   , save_xtte(:,:,jm+2), l_add = .FALSE.)
           END DO

           DEALLOCATE(save_xtte); NULLIFY(save_xtte)
#endif

        END SELECT select_cloud_param

        psc: IF (USE_PSC) THEN
           ! store information about psc relevant region locally in val_psc
           WHERE (NINT(flt_pscreg(1:kproma,:,jrow))==1) 
              val_psc(1:kproma,:)=.true.
           ELSEWHERE
              val_psc(1:kproma,:)=.false.
           END WHERE

           ! values of output variables stay unchanged for psc relevant region
           WHERE(val_psc(1:kproma,:))
              xtec(1:kproma,:,jrow)    = zxtec0(1:kproma,:)
              aclc(1:kproma,:,jrow)    = zaclc0(1:kproma,:)
              aclcac(1:kproma,:,jrow)  = zaclcac0(1:kproma,:)

              qte_3d(1:kproma,:,jrow)  = zqte0(1:kproma,:)
              tte_3d(1:kproma,:,jrow)  = ztte0(1:kproma,:)
              xlte_3d(1:kproma,:,jrow) = zxlte0(1:kproma,:)
              xite_3d(1:kproma,:,jrow) = zxite0(1:kproma,:)
           END WHERE

!qqq+ POSSIBLE INCONSISTENCY FOR SCHEMES > 2:
           !     --> TRACER TENDENCIES MUST BE RESET AS WELL, BUT BEFORE !!!
           !     MESSYTENDENCY DIAGNOSTIC IS APPLIED !!!
!qqq-           
           IF (l_wiso) THEN
              DO jwiso=1,mwiso
                 WHERE(val_psc(1:kproma,:))
                    pwisoqte(1:kproma,:,jwiso)  = zpwisoqte0(1:kproma,:,jwiso)
                    pwisoxlte(1:kproma,:,jwiso) = zpwisoxlte0(1:kproma,:,jwiso)
                    pwisoxite(1:kproma,:,jwiso) = zpwisoxite0(1:kproma,:,jwiso)
                 END WHERE
              END DO
           END IF
        END IF psc
        
#ifdef MESSYTENDENCY
        save_tte(1:kproma,:)  = &
             tte_3d(1:kproma,:,jrow) - save_tte(1:kproma,:)
        save_qte(1:kproma,:)  = &
             qte_3d(1:kproma,:,jrow) - save_qte(1:kproma,:)
        save_xlte(1:kproma,:) = &
             xlte_3d(1:kproma,:,jrow) - save_xlte(1:kproma,:)
        save_xite(1:kproma,:) = &
             xite_3d(1:kproma,:,jrow) - save_xite(1:kproma,:)
        
        CALL mtend_add_l(my_handle, mtend_id_t,  save_tte, l_add = .FALSE.)
        CALL mtend_add_l(my_handle, mtend_id_q,  save_qte, l_add = .FALSE.)
        CALL mtend_add_l(my_handle, mtend_id_xl, save_xlte, l_add = .FALSE.)
        CALL mtend_add_l(my_handle, mtend_id_xi, save_xite, l_add = .FALSE.)

        IF (l_wiso) THEN
           DO jwiso=1,mwiso
              save_wiso_qte(1:kproma,:,jwiso)  = &
                   pwisoqte(1:kproma,:,jwiso) - &
                   save_wiso_qte(1:kproma,:,jwiso)
              save_wiso_xlte(1:kproma,:,jwiso) = &
                   pwisoxlte(1:kproma,:,jwiso) - &
                   save_wiso_xlte(1:kproma,:,jwiso)
              save_wiso_xite(1:kproma,:,jwiso) = &
                   pwisoxite(1:kproma,:,jwiso) - &
                   save_wiso_xite(1:kproma,:,jwiso)
              CALL mtend_add_l(my_handle, idiso(i_vap,jwiso) &
                   ,  save_wiso_qte(1:kproma,:,jwiso), l_add = .FALSE.)
              CALL mtend_add_l(my_handle, idiso(i_liq,jwiso) &
                   ,  save_wiso_xlte(1:kproma,:,jwiso), l_add = .FALSE.)
              CALL mtend_add_l(my_handle, idiso(i_ice,jwiso) &
                   ,  save_wiso_xite(1:kproma,:,jwiso), l_add = .FALSE.)
           END DO
        END IF
#else
        sdtdt_cloud(1:kproma,:,jrow) = tte_3d(1:kproma,:,jrow) - &
             sdtdt_cloud(1:kproma,:,jrow)
#endif
        
     ENDIF no_lstart

     IF (lookupoverflow) THEN 
        do jk=1,nlev
           do jl=1,kproma
              ztp1(jl,jk) = tm1(jl,jk,jrow) + tte_3d(jl,jk,jrow)*ztmst 
              if ( INT(ztp1(jl,jk)*1000.) <jptlucu1 .OR. &
                   INT(ztp1(jl,jk)*1000.) >jptlucu2)     &
                   print*, jk, jl,ztp1(jl,jk)*1000., tm1(jl,jk,jrow), &
                   tte_3d(jl,jk,jrow)*ztmst, ztvm1(jl,jk), qm1(jl,jk,jrow),   &
                   xlm1(jl,jk,jrow), xim1(jl,jk,jrow)
           enddo
        enddo
        CALL error_bi ('cloud_convec - lookuperror', substr)
     ENDIF

       IF (ltdiag) THEN
          ! store fields for COND
          pdiga(1:kproma,:,jrow) = pdiga(1:kproma,:,jrow) + &
               tte_3d(1:kproma,:,jrow)
       ENDIF

#ifdef CESM1
       ! for CESM, averaging is done through CHANNEL 
       ! instead of the echam laccu stream objects routines.
       ! Therefore, accumulation in smcl routines is wrong and corrected here:
       aprl(1:kproma,jrow) = rsfl_2d(1:kproma,jrow)+ssfl_2d(1:kproma,jrow)
       ! convect_convec is before cloud_convec
       aprs(1:kproma,jrow) = aprs(1:kproma,jrow) + ssfl_2d(1:kproma,jrow)
       aclcac(1:kproma,:,jrow) = aclcac(1:kproma,:,jrow)/ztmst*2._dp
       aclcov(1:kproma,jrow)   = aclcov(1:kproma,jrow)/ztmst*2._dp
       qvi(1:kproma,jrow)      = qvi(1:kproma,jrow)/ztmst*2._dp
       xlvi(1:kproma,jrow)     = xlvi(1:kproma,jrow)/ztmst*2._dp
       xivi(1:kproma,jrow)     = xivi(1:kproma,jrow)/ztmst*2._dp

       aprflux(1:kproma,jrow)  = aprflux(1:kproma,jrow) &
            + rsfl_2d(1:kproma,jrow)+ssfl_2d(1:kproma,jrow)
#endif

       IF (ALLOCATED(invb))    DEALLOCATE(invb)
       IF (ALLOCATED(zbetaa))  DEALLOCATE(zbetaa)
       IF (ALLOCATED(zbetab))  DEALLOCATE(zbetab)
       IF (ALLOCATED(zbetass)) DEALLOCATE(zbetass)

     END SUBROUTINE cloud_convec

!===============================================================================

     SUBROUTINE CPL_SPECIES_SI

       ! This subroutine should determine all the required index information
       ! for aerosol cloud coupling

       USE messy_main_mpi_bi,           ONLY: p_parallel_io
       USE messy_main_blather_bi,       ONLY: error_bi
       USE messy_cloud_mem,             ONLY: nmod, mode, nsol, nfrzmod, &
                                              freez_spec, comp_str
       USE messy_main_tracer_mem_bi,    ONLY: ntrac => ntrac_gp, ti_gp
       USE messy_main_tracer,           ONLY: ON, NUMBERDENSITY,          &
                                              I_AEROSOL_MODE, AEROSOL,    &
                                              I_AEROSOL_SOL, R_MOLARMASS, &
                                              R_AEROSOL_DENSITY
       USE messy_main_tools,            ONLY: match_wild

       INTEGER :: jm, jt, idx, ji, counter
       CHARACTER(LEN=*), PARAMETER :: substr='CPL_SPECIES_E5'

       nmod = SIZE(SIGMA)
       nsol = 0

       ALLOCATE (mode(nmod))

       DO jm=1,nmod
          mode(jm)%sigma = sigma(jm)
       ENDDO

       SELECT CASE (TRIM(aer_stream))
       CASE ('made3_gp')
          DO jm=1,nmod
             mode(jm)%crdiv = 1.01e-7_dp  ! [cm] set to minimum radius in the cirrus parameterization
          ENDDO
       CASE ('gmxe_gp')
          DO jm=1,nmod
             mode(jm)%crdiv = crdiv(jm)
          ENDDO
       CASE DEFAULT
          CALL error_bi('coupling for ' // aer_stream // ' not available' &
               ,substr)
       END SELECT

       ! determine the number tracer indices
       DO jt=1,ntrac
         IF (ti_gp(jt)%tp%ident%medium /=AEROSOL) CYCLE
         IF (MATCH_WILD( 'PASSAER*', ti_gp(jt)%tp%ident%basename) ) CYCLE
         IF (ti_gp(jt)%tp%ident%quantity == NUMBERDENSITY) THEN
           idx = ti_gp(jt)%tp%meta%cask_i(I_AEROSOL_MODE)
           mode(idx)%nr = jt
         ENDIF
       END DO

       ! check if a mode is hydrophobic or hydrophilic
       DO jm=1,nmod
         idx = mode(jm)%nr 
         IF (ti_gp(idx)%tp%meta%cask_i(I_AEROSOL_SOL) == ON ) THEN
           mode(jm)%sol=.TRUE.
           nsol = nsol + 1
         ELSE
           mode(jm)%sol=.FALSE.
         ENDIF
       END DO

       ! check for nucleation mode(s)
       ! and determine the number of modes used for freezing calculations
       DO jm=1,nmod
          SELECT CASE (TRIM(aer_stream))
          CASE ('made3_gp')
             mode(jm)%nucl = .FALSE.
          CASE ('gmxe_gp')
             IF (mode(jm)%crdiv > 5.e-9_dp ) THEN
                mode(jm)%nucl=.FALSE.
                NFRZMOD = NFRZMOD + 1
             ELSE
                mode(jm)%nucl=.TRUE.
             END IF
          CASE DEFAULT
            CALL error_bi('coupling for ' // aer_stream // ' not available' &
                 , substr)
          END SELECT
       END DO

! Added MADE and MADE3 cases
       SELECT CASE (TRIM(aer_stream))
       CASE ('made3_gp')
          NFRZMOD = 3
       CASE ('gmxe_gp')

         ! NFRZMOD should be only half of the modes to account for
         !         same size classes in hydrophobic and hydrophlic
         !         categories
         ! Therefore, devide the value by 2 (for the standard case this 
         ! defaults back to a value of 3)
         !     
         NFRZMOD = NFRZMOD / 2
       CASE DEFAULT
          CALL error_bi('coupling for ' // aer_stream // ' not available' &
               , substr)
       END SELECT

       ! For the Kuebbeler et al. (2014) scheme, always NFRZMOD = 1
       if (cloud_param == 5) NFRZMOD = 1

       DO jm=1,nmod
         mode(jm)%no_freezspec = 0
         DO jt=1,ntrac
           IF (ti_gp(jt)%tp%ident%medium /=AEROSOL) CYCLE
           IF (ti_gp(jt)%tp%meta%cask_i(I_AEROSOL_MODE) /= jm) CYCLE
           IF (MATCH_WILD( 'PASSAER*', ti_gp(jt)%tp%ident%basename) ) CYCLE
           DO ji=1,freez_spec
             IF (ti_gp(jt)%tp%ident%basename == (TRIM(COMP_STR(ji))) ) &
               mode(jm)%no_freezspec =  mode(jm)%no_freezspec + 1
           END DO
         END DO
         ALLOCATE(mode(jm)%freezspec(mode(jm)%no_freezspec))
       END DO

       DO jm=1,nmod
         DO ji=1,mode(jm)%no_freezspec
           mode(jm)%freezspec(ji)%idt = 0
         END DO
         
         counter = 0
         DO jt=1,ntrac
           IF (ti_gp(jt)%tp%ident%medium /=AEROSOL) CYCLE
           IF (ti_gp(jt)%tp%meta%cask_i(I_AEROSOL_MODE) /= jm) CYCLE
           IF (MATCH_WILD( 'PASSAER*', ti_gp(jt)%tp%ident%basename) ) CYCLE
           DO ji=1,freez_spec
             IF (ti_gp(jt)%tp%ident%basename==(TRIM(COMP_STR(ji))) ) THEN
               counter = counter + 1
               mode(jm)%freezspec(counter)%idt = jt
               mode(jm)%freezspec(counter)%molarmass = &
                 ti_gp(jt)%tp%meta%cask_R(R_MOLARMASS)
               mode(jm)%freezspec(counter)%density = &
                 ti_gp(jt)%tp%meta%cask_R(R_AEROSOL_DENSITY)
               mode(jm)%freezspec(counter)%name = &
                 ti_gp(jt)%tp%ident%basename
             END IF
           END DO
         END DO
       END DO

       ! span up the struct with information about the aerosol compounds
       ! in each mode
       DO jm=1,nmod
         mode(jm)%no_aerspec = 0
         DO jt=1,ntrac
           IF (ti_gp(jt)%tp%ident%medium /=AEROSOL) CYCLE
           IF (ti_gp(jt)%tp%ident%quantity == NUMBERDENSITY) CYCLE
           IF (ti_gp(jt)%tp%meta%cask_i(I_AEROSOL_MODE) /= jm) CYCLE
           IF (MATCH_WILD( 'PASSAER*', ti_gp(jt)%tp%ident%basename) ) CYCLE
           mode(jm)%no_aerspec =  mode(jm)%no_aerspec + 1
         END DO
         ALLOCATE(mode(jm)%aerspec(mode(jm)%no_aerspec))
         counter = 0
         DO jt=1,ntrac
           IF (ti_gp(jt)%tp%ident%medium /=AEROSOL) CYCLE
           IF (ti_gp(jt)%tp%ident%quantity == NUMBERDENSITY) CYCLE
           IF (ti_gp(jt)%tp%meta%cask_i(I_AEROSOL_MODE) /= jm) CYCLE
           IF (MATCH_WILD( 'PASSAER*', ti_gp(jt)%tp%ident%basename) ) CYCLE
           counter = counter + 1
           mode(jm)%aerspec(counter)%idt = jt
           mode(jm)%aerspec(counter)%density = &
             ti_gp(jt)%tp%meta%cask_R(R_AEROSOL_DENSITY)
           mode(jm)%aerspec(counter)%molarmass = &
             ti_gp(jt)%tp%meta%cask_R(R_MOLARMASS)
         ENDDO
       ENDDO

! Print a summary table
      if (p_parallel_io) then
         write(*,*) "====== COUPLING IN CPL_SPECIES_E5 ======"
          do jm=1,nmod
             write(*,*) "-----------------------------------------"
             write(*,*) "Mode sol nucl #freez #aero sigma  #idx"
             write(*,'(3X,I1,3X,L1,4X,L1,5X,I1,5X,I2,4X,F4.1,1X,I3)') &
                  jm,mode(jm)%sol, mode(jm)%nucl, &
                  mode(jm)%no_freezspec, mode(jm)%no_aerspec, &
                  mode(jm)%sigma, mode(jm)%nr
             do ji=1,mode(jm)%no_aerspec
                write(*,'(3X,A20,1X,I3,1X,A10)') &
                      "--> Aerosol species", mode(jm)%aerspec(ji)%idt, &
                      ti_gp(mode(jm)%aerspec(ji)%idt)%tp%ident%basename
             end do
             write (*,*) " "
             do ji=1,mode(jm)%no_freezspec
                write(*,'(3X,A20,1X,I3,1X,A10)') &
                      "--> Freezing species", mode(jm)%freezspec(ji)%idt, &
                      mode(jm)%freezspec(ji)%name
             end do
          enddo
          write(*,*) "----------------------------------------"
          write(*,*) "NFRZMOD = ",NFRZMOD
          write(*,*) "----------------------------------------"
       endif
    
     END SUBROUTINE CPL_SPECIES_SI

!==================================================================================

     SUBROUTINE prepare_freezing(kproma, nlev, jrow, pxtp1, N_act)

       ! This subroutine calculates the values required for the freezing scheme

       USE messy_main_blather_bi,       ONLY: error_bi
       USE messy_main_constants_mem,    ONLY: M_air
       USE messy_main_tracer_mem_bi,    ONLY: ntrac => ntrac_gp
       USE messy_cloud_mem,             ONLY: nmod, mode,   &
                                              du_nfrac_thr, &
                                              nbcsol_strat, &
                                              ndusol_strat, nbcsol_cirrus, &
                                              nbctagsol_cirrus, &
                                              ndusol_cirrus, nbcinsol,   &
                                              nbctaginsol,      &
                                              nduinsolai, nduinsolci,   &
                                              nocinsol, nocsolks,       &
                                              nocsolas, nocsolcs,       &
                                              naerinsol, naersol,       &
                                              iaitm, iaiti, iaits,      &
                                              iaccm, iacci, iaccs,      &
                                              icoam, icoai, icoas
       USE messy_main_grid_def_bi,      ONLY: grmass, grvol

       INTEGER,  INTENT(IN) :: kproma, nlev, jrow
       REAL(dp), INTENT(IN) :: pxtp1(kproma,nlev,ntrac)
       REAL(dp), INTENT(IN) :: N_act(kproma,nlev,nmod)

       REAL(dp), DIMENSION(kproma,nlev) :: mass, vol, val, zrho
       REAL(dp), DIMENSION(kproma,nlev) :: mass2, zratio, numb, numb2, help
       REAL(dp)                         :: molar, dens, delta

       INTEGER :: jl, jk, ji, jt, idt, idx, jm
       CHARACTER(LEN=*), PARAMETER :: substr='prepare_freezing'      

       ndusol_strat(1:kproma,:,jrow) = 0._dp
       nbcsol_strat(1:kproma,:,jrow) = 0._dp
       ndusol_cirrus(1:kproma,:,jrow)    = 0._dp
       nbcsol_cirrus(1:kproma,:,jrow)    = 0._dp
       nbctagsol_cirrus(1:kproma,:,jrow) = 0._dp
       nocsolks(1:kproma,:,jrow)     = 0._dp
       nocsolas(1:kproma,:,jrow)     = 0._dp
       nocsolcs(1:kproma,:,jrow)     = 0._dp
       naerinsol(1:kproma,:,jrow)    = 0._dp
       naersol(1:kproma,:,jrow)      = 0._dp
       nbcinsol(1:kproma,:,jrow)     = 0._dp
       nbctaginsol(1:kproma,:,jrow)  = 0._dp
       nocinsol(1:kproma,:,jrow)     = 0._dp
       nduinsolai(1:kproma,:,jrow)   = 0._dp
       nduinsolci(1:kproma,:,jrow)   = 0._dp
       do jk=1,nlev
          zrho(1:kproma,jk) = grmass(1:kproma,jk,jrow) / &
                              grvol (1:kproma,jk,jrow)
       enddo

       ! Add MADE3 case
       SELECT CASE(TRIM(aer_stream))
       CASE ('made3_gp')

          DO jm=1,nmod

             mode(jm)%wetrad => wetradius(1:kproma,:,jm,jrow)
             mode(jm)%dryrad => dryradius(1:kproma,:,jm,jrow)

             idt = mode(jm)%nr

             ! Number concentration in the mode
             IF (ASSOCIATED(mode(jm)%aernum)) DEALLOCATE(mode(jm)%aernum)
             NULLIFY(mode(jm)%aernum)
             ALLOCATE(mode(jm)%aernum(kproma,nlev))
             ! [1/mol] --> [1/m3]
             mode(jm)%aernum(1:kproma,:) = pxtp1(1:kproma,:,idt)/m_air * &
                  1000._dp * zrho(1:kproma,:)

             ! Total number of particles in purely soluble modes [1/m3]
             ! (for homogeneous freezing)
             IF (jm == iaits .OR. jm == iaccs .OR. jm == icoas) THEN
                naersol(1:kproma,:,jrow) = &
                  naersol(1:kproma,:,jrow) + mode(jm)%aernum(1:kproma,:)
             ENDIF

             ! Total number of particles in insolube modes [1/m3]
             IF (.NOT.mode(jm)%sol) THEN
                naerinsol(1:kproma,:,jrow) = &
                  naerinsol(1:kproma,:,jrow) + mode(jm)%aernum(1:kproma,:)
             ENDIF

             ! Dust mass in the mode
             mass(1:kproma,:) = 0._dp
             DO jt=1,mode(jm)%no_freezspec
                IF (TRIM(mode(jm)%freezspec(jt)%name) == 'DU' ) THEN
                   idx = mode(jm)%freezspec(jt)%idt
                   mass(1:kproma,1:nlev) = pxtp1(1:kproma,1:nlev,idx) * &
                        mode(jm)%freezspec(jt)%molarmass / M_air * &
                        zrho(1:kproma,:)  ! [mol/mol] --> [kg/m3]
                ENDIF
             END DO

             ! BCtag mass in the mode
             mass2(1:kproma,:) = 0._dp
             DO jt=1,mode(jm)%no_freezspec
                IF (TRIM(mode(jm)%freezspec(jt)%name) == 'BCtag' ) THEN
                   idx = mode(jm)%freezspec(jt)%idt
                   mass2(1:kproma,1:nlev) = pxtp1(1:kproma,1:nlev,idx) * &
                        mode(jm)%freezspec(jt)%molarmass / M_air * &
                        zrho(1:kproma,:)  ! [mol/mol] --> [kg/m3]
                ENDIF
             END DO

             ! Mass-to-number conversion for BCtag in the aitken modes
             ! assuming D=0.025 micron, sigma=1.55 and rho=1500 kg/m3
             ! 3.433428e+19 1/kg

             IF (jm == iaitm) THEN
                numb2(1:kproma,:) = 3.433428E19_dp * mass2(1:kproma,:)
                numb2(1:kproma,:) = &
                     MIN(numb2(1:kproma,:), mode(jm)%aernum(1:kproma,:))

                ! Cirrus clouds
                nbctagsol_cirrus(1:kproma,:,jrow) = &
                     nbctagsol_cirrus(1:kproma,:,jrow) + &
                     numb2(1:kproma,:)

                nbcsol_cirrus(1:kproma,:,jrow) = &
                     nbcsol_cirrus(1:kproma,:,jrow) + &
                     MAX(mode(jm)%aernum(1:kproma,:) - numb2(1:kproma,:),0._dp)

                ! Mixed-phase clouds
                nbcsol_strat(1:kproma,:,jrow) = &
                     nbcsol_strat(1:kproma,:,jrow) + &
                     N_act(1:kproma,:,jm)
             ENDIF

             IF (jm == iaiti) THEN
                numb2(1:kproma,:) = 3.433428E19_dp * mass2(1:kproma,:)
                numb2(1:kproma,:) = &
                     MIN(numb2(1:kproma,:), mode(jm)%aernum(1:kproma,:))

                ! Cirrus clouds
                nbctaginsol(1:kproma,:,jrow) = &
                     nbctaginsol(1:kproma,:,jrow) + &
                     numb2(1:kproma,:)

                nbcinsol(1:kproma,:,jrow) = &
                     nbcinsol(1:kproma,:,jrow) + &
                     MAX(mode(jm)%aernum(1:kproma,:) - numb2(1:kproma,:),0._dp)
             ENDIF

             ! Mass-to-number conversion for dust in the accumulation modes
             ! assuming D=0.42 micron, sigma=1.59 and rho=2500 kg/m3
             ! 3.917756e+15 1/kg

             ! Mass-to-number conversion for BCtag in the accumulation modes
             ! assuming D=0.15 micron, sigma=1.65 and rho=1500 kg/m3
             ! 1.220503e+17 1/kg

             IF (jm == iaccm) THEN
                numb(1:kproma,:) = 3.917756E15_dp * mass(1:kproma,:)
                numb(1:kproma,:) = &
                     MIN(numb(1:kproma,:), mode(jm)%aernum(1:kproma,:))  ! DU

                numb2(1:kproma,:) = 1.220503E17_dp * mass2(1:kproma,:)
                numb2(1:kproma,:) = &
                     MIN(numb2(1:kproma,:), mode(jm)%aernum(1:kproma,:)) ! BCtag

                ! Cirrus clouds
                ndusol_cirrus(1:kproma,:,jrow) = &
                     ndusol_cirrus(1:kproma,:,jrow) + numb(1:kproma,:)
                nbctagsol_cirrus(1:kproma,:,jrow) = &
                     nbctagsol_cirrus(1:kproma,:,jrow) + numb2(1:kproma,:)
                nbcsol_cirrus(1:kproma,:,jrow) = &
                     nbcsol_cirrus(1:kproma,:,jrow) + &
                     MAX(mode(jm)%aernum(1:kproma,:) - numb(1:kproma,:) - numb2(1:kproma,:),0._dp)

                ! Mixed-phase clouds
                help(1:kproma,:) = MIN(1._dp, N_act(1:kproma,:,jm) / &
                     (mode(jm)%aernum(1:kproma,:) + EPSILON(1._dp)))
                numb(1:kproma,:) = numb(1:kproma,:) * help(1:kproma,:)

                ndusol_strat(1:kproma,:,jrow) = &
                     ndusol_strat(1:kproma,:,jrow) + numb(1:kproma,:)
                nbcsol_strat(1:kproma,:,jrow) = &
                     nbcsol_strat(1:kproma,:,jrow) + &
                     MAX(N_act(1:kproma,:,jm) - numb(1:kproma,:),0._dp)
             ENDIF

             IF (jm == iacci) THEN
                numb(1:kproma,:) = 3.917756E15_dp * mass(1:kproma,:)
                numb(1:kproma,:) = &
                     MIN(numb(1:kproma,:), mode(jm)%aernum(1:kproma,:)) ! DU

                numb2(1:kproma,:) = 1.220503E17_dp * mass2(1:kproma,:)
                numb2(1:kproma,:) = &
                     MIN(numb2(1:kproma,:), mode(jm)%aernum(1:kproma,:)) ! BCtag

                nduinsolai(1:kproma,:,jrow) = &
                     nduinsolai(1:kproma,:,jrow) + numb(1:kproma,:)
                nbctaginsol(1:kproma,:,jrow) = &
                     nbctaginsol(1:kproma,:,jrow) + numb2(1:kproma,:)
                nbcinsol(1:kproma,:,jrow) = &
                     nbcinsol(1:kproma,:,jrow) + &
                     MAX(mode(jm)%aernum(1:kproma,:) - numb(1:kproma,:) - numb2(1:kproma,:),0._dp)
             ENDIF

             ! Mass-to-number conversion for dust in the coarse modes
             ! assuming D=1.30 micron, sigma=2.00 and rho=2500 kg/m3
             ! 4.001934e+13 1/kg

             IF (jm == icoam) THEN

                numb(1:kproma,:) = 4.001934E13_dp * mass(1:kproma,:)
                numb(1:kproma,:) = &
                     MIN(numb(1:kproma,:), mode(jm)%aernum(1:kproma,:))

                ! Ratio N_DU/N_aer
                zratio(1:kproma,:) = MAX(0._dp, MIN(1._dp, numb(1:kproma,:) / &
                     (mode(jm)%aernum(1:kproma,:) + EPSILON(1._dp))))

                val(1:kproma,:) = &
                     MERGE(mode(jm)%aernum(1:kproma,:), &
                           numb(1:kproma,:), &
                           zratio(1:kproma,:).ge.du_nfrac_thr)

                ! Cirrus clouds
                ndusol_cirrus(1:kproma,:,jrow) = &
                     ndusol_cirrus(1:kproma,:,jrow) + val(1:kproma,:)
                nbcsol_cirrus(1:kproma,:,jrow) = &
                     nbcsol_cirrus(1:kproma,:,jrow) + & 
                     MERGE(0._dp, &
                           MAX(mode(jm)%aernum(1:kproma,:) - &
                               val(1:kproma,:),0._dp), &
                           zratio(1:kproma,:).ge.du_nfrac_thr)

                ! Mixed-phase clouds
                help(1:kproma,:) = MIN(1._dp, N_act(1:kproma,:,jm) / &
                     (mode(jm)%aernum(1:kproma,:) + EPSILON(1._dp)))
                val(1:kproma,:) = val(1:kproma,:) * help(1:kproma,:)

                ndusol_strat(1:kproma,:,jrow) = &
                     ndusol_strat(1:kproma,:,jrow) + val(1:kproma,:)
                nbcsol_strat(1:kproma,:,jrow) = &
                     nbcsol_strat(1:kproma,:,jrow) + & 
                     MERGE(0._dp, &
                           MAX(N_act(1:kproma,:,jm) - &
                               val(1:kproma,:),0._dp), &
                           zratio(1:kproma,:).ge.du_nfrac_thr)
             ENDIF

             IF (jm == icoai) THEN
                nduinsolci(1:kproma,:,jrow) = &
                     nduinsolci(1:kproma,:,jrow) + &
                     mode(jm)%aernum(1:kproma,:)
             ENDIF

          END DO

       CASE ('gmxe_gp')          

          DO jm=1,nmod

             mode(jm)%wetrad => wetradius(1:kproma,:,jm,jrow)
             mode(jm)%dryrad => dryradius(1:kproma,:,jm,jrow)

             idt = mode(jm)%nr

             IF (ASSOCIATED(mode(jm)%aernum)) DEALLOCATE(mode(jm)%aernum)
             NULLIFY(mode(jm)%aernum)
             ALLOCATE(mode(jm)%aernum(kproma,nlev))
             ! conversion from 1/mol -> 1/kg
             mode(jm)%aernum(1:kproma,:) = pxtp1(1:kproma,:,idt)/m_air * 1000._dp
             ! conversion from 1/kg -> 1/m^3
             mode(jm)%aernum(1:kproma,:) = mode(jm)%aernum(1:kproma,:) * zrho(1:kproma,:) 

             IF (mode(jm)%sol) THEN
                IF (.NOT. mode(jm)%nucl) THEN
                  naersol(1:kproma,:,jrow) = &
                  naersol(1:kproma,:,jrow) + mode(jm)%aernum(1:kproma,:)
                ENDIF
             ELSE
                  naerinsol(1:kproma,:,jrow) = &
                  naerinsol(1:kproma,:,jrow) + mode(jm)%aernum(1:kproma,:)
             ENDIF

             DO ji=1,mode(jm)%no_freezspec
                ALLOCATE(mode(jm)%freezspec(ji)%ptr(kproma,nlev))
             ENDDO
             ! calculate total aerosol properties

             IF (mode(jm)%no_freezspec == 0) CYCLE
         
             mass(:,:) = 0._dp
             vol(:,:)  = 0._dp

             DO jt = 1,mode(jm)%no_aerspec
                idx   = mode(jm)%aerspec(jt)%idt
                molar = mode(jm)%aerspec(jt)%molarmass
                dens  = mode(jm)%aerspec(jt)%density
                DO jk=1,nlev
                   DO jl=1,kproma
                      ! total aerosol mass per mode [in kg/kg]
                      delta = pxtp1(jl,jk,idx) * molar/m_air
                      mass(jl,jk) = mass(jl,jk) + delta
                       
                      ! total volume per mode in [m^3/kg]
                      vol(jl,jk)  = vol(jl,jk) + delta / dens 
                 
                   ENDDO
                END DO
             END DO
             MASS(1:kproma,1:nlev) = MAX(1.e-20_dp, MASS(1:kproma,1:nlev) )
             VOL(1:kproma,1:nlev)  = MAX(1.e-20_dp, VOL(1:kproma,1:nlev) )
        

             ! determine ratios of individual compounds      

             DO jt=1,mode(jm)%no_freezspec
                idx = mode(jm)%freezspec(jt)%idt
           
                IF (idx/=0 ) THEN
                   val(1:kproma,1:nlev) = pxtp1(1:kproma,1:nlev,idx) * &
                        mode(jm)%freezspec(jt)%molarmass / M_air

                   IF (mode(jm)%sol) THEN
                      ! used for immersion freezing
                      ! use volume and CCNs here
                      val(1:kproma,1:nlev) = val(1:kproma,1:nlev) / &
                           mode(jm)%freezspec(jt)%density
                      zratio(1:kproma,1:nlev) = MAX(0._dp, MIN(1._dp,&
                           val(1:kproma,1:nlev) / vol(1:kproma,1:nlev)))
                      numb(1:kproma,1:nlev)  = N_act(1:kproma,1:nlev,jm)
                   ELSE
                      ! used for heterogeneous freezing
                      ! use mass here
                      zratio(1:kproma,1:nlev) = MAX(0._dp, MIN(1._dp, &
                           val(1:kproma,1:nlev) / mass(1:kproma,1:nlev)))
                      numb(1:kproma,1:nlev)  = mode(jm)%aernum(1:kproma,1:nlev)
                                        ! pxtp1(1:kproma,1:nlev,idt)
                   ENDIF

                   mode(jm)%freezspec(jt)%ptr(1:kproma,:) = &
                        zratio(1:kproma,:)**(2._dp/3._dp) * numb(1:kproma,:)
                ENDIF
             END DO

             !copy into corresponding arrays
             ! not nice but working !!!!

             DO ji=1,mode(jm)%no_freezspec
                IF (mode(jm)%sol) THEN
                   IF (TRIM(mode(jm)%freezspec(ji)%name) == 'DU' ) THEN
                      ndusol_strat(1:kproma,:,jrow) =    &
                           ndusol_strat(1:kproma,:,jrow) +  &
                           mode(jm)%freezspec(ji)%ptr(1:kproma,:)
                   ENDIF
                   IF (TRIM(mode(jm)%freezspec(ji)%name) == 'BC' ) THEN
                      nbcsol_strat(1:kproma,:,jrow) =    &
                           nbcsol_strat(1:kproma,:,jrow) +  &
                           mode(jm)%freezspec(ji)%ptr(1:kproma,:)
                   ENDIF
                   IF (TRIM(mode(jm)%freezspec(ji)%name) == 'OC' ) THEN
                      IF (jm==2) &                                      
                         nocsolks(1:kproma,:,jrow) = &               !OC in ks
                             mode(jm)%freezspec(ji)%ptr(1:kproma,:)
                      IF (jm==3) &
                         nocsolas(1:kproma,:,jrow) = &               !OC in as
                             mode(jm)%freezspec(ji)%ptr(1:kproma,:)
                      IF (jm==4) &     
                         nocsolcs(1:kproma,:,jrow) = &               !OC in cs
                             mode(jm)%freezspec(ji)%ptr(1:kproma,:)
                   END IF
                ELSE
                   IF (TRIM(mode(jm)%freezspec(ji)%name) == 'DU' ) THEN
                      IF (jm==6) &
                           nduinsolai(1:kproma,:,jrow) = &
                           mode(jm)%freezspec(ji)%ptr(1:kproma,:)
                      IF (jm==7) &
                           nduinsolci(1:kproma,:,jrow) = &
                           mode(jm)%freezspec(ji)%ptr(1:kproma,:)
                   ENDIF
                   IF (TRIM(mode(jm)%freezspec(ji)%name) == 'BC' ) THEN
                      nbcinsol(1:kproma,:,jrow) = &
                           mode(jm)%freezspec(ji)%ptr(1:kproma,:)
                   ENDIF
                   IF (TRIM(mode(jm)%freezspec(ji)%name) == 'OC' ) THEN
                      nocinsol(1:kproma,:,jrow) = &
                          mode(jm)%freezspec(ji)%ptr(1:kproma,:)
                   END IF
                ENDIF
             END DO

             nbcsol_cirrus(1:kproma,:,jrow) = nbcsol_strat(1:kproma,:,jrow)
             ndusol_cirrus(1:kproma,:,jrow) = ndusol_strat(1:kproma,:,jrow)

             DO ji=1,mode(jm)%no_freezspec
                IF ( ASSOCIATED(mode(jm)%freezspec(ji)%ptr) )&
                     DEALLOCATE(mode(jm)%freezspec(ji)%ptr)
             ENDDO
          END DO

       CASE DEFAULT
          CALL error_bi('coupling for ' // aer_stream // ' not available' &
               , substr)
       END SELECT

     END SUBROUTINE prepare_freezing

!=======================================================================
     
     SUBROUTINE cloud_activate_ARG

       USE messy_main_grid_def_mem_bi, ONLY: jrow, kproma, nlev
       USE messy_main_grid_def_bi,     ONLY: grmass, grvol
       USE messy_main_data_bi,         ONLY: apm1, vervel_3d, &
                                             tke, tm1, tte_3d
       USE messy_main_timer,          ONLY: time_step_len
       USE messy_main_tracer_mem_bi,  ONLY: pxtte => qxtte, pxtm1 => qxtm1, &
                                            ntrac => ntrac_gp
       USE messy_cloud_droplet,       ONLY: cloud_droplet_ARG
       USE messy_main_constants_mem,  ONLY: g, m_air

       INTRINSIC :: ASSOCIATED

       INTEGER  :: jk, jl, ji, jt

       CHARACTER(LEN=*), PARAMETER::substr='cloud_activate_ARG'

       REAL(dp) :: zrho(nproma,nlev), temp(nproma,nlev)
       REAL(dp) :: aer_crit(nproma,nlev,ksol), &
                   S_max(nproma,nlev), velo(nproma,nlev)
       REAL(dp) :: pxtp1(nproma,nlev,0:ntrac)
       REAL(dp) :: help
       REAL(dp) :: S_crit(ksol), inp_dryrad(ksol), inp_sigma(ksol)
       REAL(dp) :: nfrac(ksol), mfrac(ksol), kfrac(ksol), cmr_to_mmr(ksol)

       pxtp1(:,:,0) = 0._dp
       do jk=1,nlev
         zrho(1:kproma,jk) = grmass(1:kproma,jk,jrow) / &
                             grvol (1:kproma,jk,jrow)
       enddo
       DO ji=1,ksol
         IF (ASSOCIATED(cmr2mmr)) THEN
           cmr_to_mmr(ji) = cmr2mmr(ji)
         ELSE
           cmr_to_mmr(ji) = 1._dp
         END IF
       ENDDO
       DO jk=1,nlev
         DO jl=1,kproma
           temp(jl,jk) = tm1(jl,jk,jrow) + tte_3d(jl,jk,jrow) * time_step_len
           pxtp1(jl,jk,1:ntrac) = pxtm1(jl,jk,1:ntrac) + &
                                  pxtte(jl,jk,1:ntrac) * time_step_len
               
           velo(jl,jk) = - vervel_3d(jl,jk,jrow) / (zrho(jl,jk) * g)
           velo(jl,jk) = velo(jl,jk) + 1.33_dp * sqrt(tke(jl,jk,jrow))
           velo(jl,jk) = MAX(velo(jl,jk),0._dp)
           DO JI = 1,ksol
             DO jt = 1,nspec
               mass(JI,JT) = MAX(pxtp1(jl,jk,aer_idx(JI,JT)),1.e-28_dp)

               ! Convert to mass mixing ratio [kg/kg]
               mass(ji,jt) = mass(ji,jt) * molar(jt) / M_air
             END DO
           END DO
           DO JI = 1,ksol
             ! conversion from 1/mol to 1 /m^3
             ! 1/mol -> 1/kg : 1/mol / M_air[g/mol] * 1000
             ! 1/kg  -> 1/m^3: 1/kg * rho [kg/m^3]
             IF (num_idx(JI) /= 0) THEN
               help = pxtp1(jl,jk,num_idx(JI)) / M_air * 1000._dp
               help = help * zrho(jl,jk)
             ELSE
               IF (ASSOCIATED(anumber)) &
                 help = anumber(jl,jk,JI,jrow) * 1.e6_dp
             ENDIF
             AEROSOL_NUM(JI) = MAX(help,1.e-6_dp)
             inp_dryrad(ji) = dryradius(jl,jk,sol_idx(ji),jrow)
             inp_sigma(ji) = sigma(sol_idx(ji))
             S_crit(JI) = 0._dp

             IF (ASSOCIATED(S_crit_in(ji)%ptr)) THEN
               S_crit(JI) = S_crit_in(ji)%ptr(jl,jk,jrow)/100.
             ENDIF
          ENDDO

           CALL cloud_droplet_ARG(temp(jl,jk),     apm1(jl,jk),                 &
                                  velo(jl,jk),     inp_dryrad(:),               &
                                  aerdensity(:),   inp_sigma(:),                &
                                  aerosol_num(:),  mass(:,:),                   &
                                  nu(:),           PHI(:),                      &
                                  EPS(:), KAPPA(:), molar(:),                   &
                                  S_crit(:), sup_sat_scheme,                    &

                                  nfrac(:),        S_max(jl,jk),                &
                                  aer_crit(jl,jk,:),                            &
                                  acdnc_5(jl,jk,jrow),                          &
                                  jl,          nlev,                            &
                                  ksol,        nspec,    kfrac(:),              &
                                  mfrac(:),    cmr_to_mmr(:) )
           
          DO ji=1,ksol
             arg_frac(jl,jk,ji,jrow)  = nfrac(ji)
             arg_mfrac(jl,jk,ji,jrow) = mfrac(ji)
             arg_frac2(jl,jk,ji,jrow) = kfrac(ji)
             arg_num(jl,jk,ji,jrow)  = aerosol_num(ji) * nfrac(ji)
             arg_rad(jl,jk,ji,jrow)  =  aer_crit(jl,jk,ji)
             arg_scrit(jl,jk,ji,jrow) = s_crit(ji)
           END DO
           Scrit(jl,jk,jrow) = S_max(jl,jk)
         enddo
       enddo

     END SUBROUTINE cloud_activate_ARG

!=======================================================================


     SUBROUTINE cloud_activate_LIN

       USE messy_main_grid_def_mem_bi, ONLY: jrow, kproma, nlev
       USE messy_main_grid_def_bi,     ONLY: grmass, grvol
       USE messy_main_data_bi,         ONLY: vervel_3d, &
                                             tke, tm1, tte_3d
       USE messy_main_timer,          ONLY: time_step_len
       USE messy_main_tracer_mem_bi,  ONLY: pxtte => qxtte, pxtm1 => qxtm1, &
                                            ntrac => ntrac_gp
       USE messy_cloud_droplet,       ONLY: cloud_droplet_LIN
       USE messy_cloud_mem,           ONLY: idt_CDNC

       USE messy_main_constants_mem,  ONLY: g, m_air


       INTRINSIC :: ASSOCIATED

       INTEGER  :: jk, jl, ji

       REAL(dp) :: zrho(nproma,nlev), temp(nproma,nlev)
       REAL(dp) :: velo(nproma,nlev)
       REAL(dp) :: pxtp1(nproma,nlev,0:ntrac)
       REAL(dp) :: help, N_old
       REAL(dp) :: lfrac(ksol), mfrac(ksol), cmr_to_mmr(ksol)
       REAL(dp) :: inp_sigma(ksol), inp_wetrad(ksol)

       pxtp1(:,:,0) = 0._dp
       do jk=1,nlev
         zrho(1:kproma,jk) = grmass(1:kproma,jk,jrow) / &
                             grvol (1:kproma,jk,jrow)
       enddo
       DO ji=1,ksol
         IF (ASSOCIATED(cmr2mmr)) THEN
           cmr_to_mmr(ji) = cmr2mmr(ji)
         ELSE
           cmr_to_mmr(ji) = 1._dp
         END IF
       ENDDO
       DO jk=1,nlev
         DO jl=1,kproma
           temp(jl,jk) = tm1(jl,jk,jrow) + tte_3d(jl,jk,jrow) * time_step_len
           pxtp1(jl,jk,1:ntrac) = pxtm1(jl,jk,1:ntrac) + &
                                  pxtte(jl,jk,1:ntrac) * time_step_len
               
           velo(jl,jk) = - vervel_3d(jl,jk,jrow) / (zrho(jl,jk) * g)
           velo(jl,jk) = velo(jl,jk) + 1.33_dp * sqrt(tke(jl,jk,jrow))
           velo(jl,jk) = MAX(velo(jl,jk),0._dp)
           DO JI = 1,ksol
             ! conversion from 1/mol to 1 /m^3
             ! 1/mol -> 1/kg : 1/mol / M_air[g/mol] * 1000
             ! 1/kg  -> 1/m^3: 1/kg * rho [kg/m^3]
             IF (num_idx(JI) /= 0) THEN
               help = pxtp1(jl,jk,num_idx(JI)) / M_air * 1000._dp
               help = help * zrho(jl,jk)
             ELSE
               IF (ASSOCIATED(anumber)) &
                 help = anumber(jl,jk,JI,jrow) * 1.e6_dp
             ENDIF
             AEROSOL_NUM(JI) = MAX(help,1.e-6_dp)

             inp_wetrad(ji) = wetradius(jl,jk,sol_idx(ji),jrow)
             inp_sigma(ji) = sigma(sol_idx(ji))

           ENDDO
           IF (idt_cdnc /= 0) THEN
             N_old = pxtp1(jl,jk,idt_cdnc) / M_air * 1000._dp * zrho(jl,jk)
             N_old = MAX(N_old,0._dp)
           ELSE
             N_old = MAX(acdnc_6(jl,jk,jrow),0._dp)
           ENDIF

           CALL cloud_droplet_Lin(ksol,  AEROSOL_NUM(:), velo(jl,jk),      &
                                        N_old, inp_sigma(:), inp_wetrad(:),&

                                  acdnc_6(jl,jk,jrow),                     &
                                  np_pot_activ(jl,jk,jrow), temp(jl,jk),   &
                                  lfrac(:), mfrac(:), cmr_to_mmr(:) )
           DO ji=1,ksol
             lin_frac(jl,jk,ji,jrow)  = lfrac(ji)
             lin_mfrac(jl,jk,ji,jrow) = mfrac(ji)
             lin_num(jl,jk,ji,jrow)   = aerosol_num(ji) * lfrac(ji)
           END DO
           
         ENDDO
       ENDDO

     END SUBROUTINE cloud_activate_LIN

!=======================================================================

     SUBROUTINE cloud_activate_UAF

       USE messy_main_grid_def_mem_bi, ONLY: jrow, kproma, nlev
       USE messy_main_grid_def_bi,     ONLY: grmass, grvol
       USE messy_main_data_bi,         ONLY: vervel_3d,   &
                                             tke, tm1, tte_3d, press_3d
       USE messy_main_timer,          ONLY: time_step_len
       USE messy_main_tracer_mem_bi,  ONLY: pxtte => qxtte, pxtm1 => qxtm1, &
                                            ntrac => ntrac_gp
       USE messy_main_constants_mem,  ONLY: g, m_air

       include 'messy_cloud_parameters_uaf.inc'

       INTRINSIC :: ASSOCIATED

       INTEGER  :: jk, jl, ji, jt, jn

       REAL(dp) :: zrho(nproma,nlev), temp(nproma,nlev),zpress(nproma,nlev)
       REAL(dp) :: S_max(nproma,nlev), velo(nproma,nlev)
       REAL(dp) :: pxtp1(nproma,nlev,0:ntrac)
       REAL(dp) :: help
       REAL(dp) :: nfrac(ktot)
       REAL(dp) :: mfrac(ktot)
       REAL(dp) :: modmass(ktot),dumass(ktot),DPGI(ktot),AKKI(ktot),ei(ktot),A,B,SG(ktot),ACCOM,SIGW,NDM(ktot) !VAK
       REAL(dp) :: modvol(ktot),duvol(ktot)
       LOGICAL  :: FHH2

       mfrac(:)     = 0._dp
       pxtp1(:,:,0) = 0._dp
        do jk=1,nlev
          zrho(1:kproma,jk) = grmass(1:kproma,jk,jrow) / &
                              grvol (1:kproma,jk,jrow)
        enddo
           DO jk=1,nlev
             DO jl=1,kproma
               temp(jl,jk) = tm1(jl,jk,jrow) + tte_3d(jl,jk,jrow) * time_step_len
               pxtp1(jl,jk,1:ntrac) = pxtm1(jl,jk,1:ntrac) + &
                                      pxtte(jl,jk,1:ntrac) * time_step_len

               zpress  (jl,jk)      = press_3d(jl,jk,jrow)     ! [Pa]

               velo(jl,jk) = - vervel_3d(jl,jk,jrow) / (zrho(jl,jk) * g)
               velo(jl,jk) = velo(jl,jk) + 0.7_dp * sqrt(tke(jl,jk,jrow)) !VAK: instead of : 1.33_dp * sqrt(tke(jl,jk,jrow))
               velo(jl,jk) = MAX(velo(jl,jk),0._dp)
               DO JI = 1,ktot
                 modmass(JI)=0._dp
                 modvol(JI)=0._dp
                 !bulkvol(JI)=0._dp
                 dumass(JI)  =0._dp
                 DO jt = 1,nspec
                   mass(JI,JT) = MAX(pxtp1(jl,jk,aer_idx(JI,JT)),1.e-28_dp)
                   modmass(JI) = modmass(JI)+mass(JI,JT)
                   modvol(JI)  = modvol(JI)+mass(JI,JT)*molar(JT)/aerdensity(JT)
                   SELECT CASE (TRIM(base_name(JT)))
                    CASE('DU')
                     dumass(JI)  = mass(JI,JT)
                     duvol(JI)  = mass(JI,JT)*molar(JT)/aerdensity(JT)
                   END SELECT
                 END DO
               END DO

               DO JI = 1,ktot
                 ! conversion from 1/mol to 1 /m^3
                 ! 1/mol -> 1/kg : 1/mol / M_air[g/mol] * 1000
                 ! 1/kg  -> 1/m^3: 1/kg * rho [kg/m^3]
                 IF (num_idx(JI) /= 0) THEN
                   help = pxtp1(jl,jk,num_idx(JI)) / M_air * 1000._dp
                   help = help * zrho(jl,jk)
                 ELSE
                   IF (ASSOCIATED(anumber)) &
                     help = anumber(jl,jk,JI,jrow) * 1.e6_dp
                 ENDIF
                 AEROSOL_NUM(JI) = MAX(help,1.e-6_dp)
               ENDDO


!C ** INITIALIZE: CALCULATE GAUSS QUADRATURE POINTS *********************
!C
!      CALL PROPS
!      CALL GAULEG (XGS, WGS, Npgauss)
!C
!C ** SPECIFY INPUT FOR CALCULATING CCN SPECTRUM ************************
!C
      ACCOM  = 0.06d0           ! Accommodation coefficient
      SIGW  = 0.0           ! STDEV of updraft distribution
                            ! 0 = calculate for single updraft (no PDF)
      A = 2.25d0
      B = 1.2d0

      AKKI(:)=0.d0
      ei(:)=0.d0
!      modmass(:)=0.d0

      do jn=1,ksol
          ei(jn)=dumass(jn)/modmass(jn)
        if (dumass(jn).lt.1.e-27_dp.or.ei(jn).lt.1.e-02_dp) then
          MODE (jn) = 1             !     MODE = 1 DENOTES KOHLER MODE (soluble aerosol)
        else
          MODE (jn) = 2             !     MODE = 2 DENOTES unified MODE (insoluble+soluble aerosol)
        end if
      end do
      do jn=ksol+1,ktot
          ei(jn)=dumass(jn)/modmass(jn)
        if (dumass(jn).lt.1.e-27_dp.or.ei(jn).lt.1.e-02_dp) then
          MODE (jn) = 1             !     MODE = 1 DENOTES KOHLER MODE (soluble aerosol)
          aerosol_num(jn)=1.e-6_dp
        else
          MODE (ksol+2:ktot) = 2        !     MODE = 2 and ei = 1 DENOTEs FHH MODE (insoluble aerosol)
!!!!          aerosol_num(jn)=ei(jn)*aerol_num(j)
          ei (jn) = 1._dp
        end if
      end do
!C

      do jn = 1,ktot
        if (jn.le.ksol) then
        AKKI(jn) = ( kappa_uaf(jn)%ptr(jl,jk,jrow)*modvol(jn) - 0.03*duvol(jn) ) /  &     !0.03 = kappa of dust
        &          (modvol(jn)-duvol(jn))    
        end if
        DPGI(jn)=2._dp*dryradius(jl,jk,jn,jrow)
      end do

               DO JI = ksol+1,ktot !VAK
                 FHH2=.false.
                 DO jt = 1,nspec
                  if(FHH(jt)) FHH2=.true.
                 END DO
                 if(.not.FHH2) AEROSOL_NUM(JI) = 1.e-6_dp
               END DO             !VAK

!C ---------------------------------------------------------
!C ---------------------------------------------------------

!C        Running the parameterization:
!C        1) Call the CCN spectrum     
            CALL CCNSPEC (aerosol_num,DPGI,sigma,temp(jl,jk),zpress(jl,jk),ktot,AKKI,ei,A,B,SG)
!C        2) Call the activation routine that computes SMAX and NACT. 
            CALL PDFACTIV (velo(jl,jk),aerosol_num,AKKI,ei,A,B,ACCOM,SG,SIGW,temp(jl,jk),zpress(jl,jk),&
     &              acdnc_7(jl,jk,jrow),S_max(jl,jk),NDM)

           DO ji=1,ktot
             nfrac(ji) = NDM(ji)/aerosol_num(ji)
             uaf_frac(jl,jk,ji,jrow)  = nfrac(ji)
             uaf_mfrac(jl,jk,ji,jrow) = mfrac(ji)
             uaf_num(jl,jk,ji,jrow)  = NDM(ji)
             uaf_AKK(jl,jk,ji,jrow) = AKKI(ji)
             uaf_TPI(jl,jk,ji,jrow) = aerosol_num(ji)
             uaf_DPG(jl,jk,ji,jrow) = DPGI(ji)
             uaf_SIG(jl,jk,ji,jrow) = sigma(ji)
             uaf_ei(jl,jk,ji,jrow) = ei(ji)
           END DO
           Scrit(jl,jk,jrow) = S_max(jl,jk)*100._dp
           WPARC(jl,jk,jrow) = velo(jl,jk)
             enddo
           enddo
     END SUBROUTINE cloud_activate_UAF

!=======================================================================

     SUBROUTINE set_cloud_parameters(nlev, lipcc)

       USE messy_cloud_ori,       ONLY: ccsaut, ccraut
       USE messy_cloud_mem,       ONLY: ncdnc, nicnc

       INTEGER, INTENT(IN) :: nlev
       LOGICAL, INTENT(IN) :: lipcc

       !-- overwrite values for coupled CDNC/ICNC cloud schemce

       IF (ncdnc>0) THEN
          ccsaut = 250._dp
          ccraut = 3.7_dp
       ENDIF
       IF (nicnc>0)  THEN
          IF (nlev == 19) THEN
             ccsaut = 350._dp
             ccraut = 5.8_dp
          ENDIF
          IF (nlev == 31) THEN
             ccsaut=350._dp
             ccraut = 6.3_dp
             IF (lipcc) THEN
                ccsaut = 400._dp
                ccraut = 10.0_dp
             ENDIF
          ENDIF
       ENDIF

     END SUBROUTINE set_cloud_parameters

!=======================================================================

     SUBROUTINE calc_vervel_penner18(kproma, nlev, jrow, wvel)
     
       ! Implement the Penner et al. (2018) Ansatz for sub-grid scale vertical
       ! velocity (see their Supporting Information).
       !
       ! References:
       !   Penner et al., J. Geophys. Res. (2018)
       !   Podglajen et al., Geophys. Res. Lett. (2016)
       !   Gary, Atmos. Chem. Phys. (2006)
       !   Gary, Atmos. Chem. Phys. (2008)
       !
       USE messy_main_grid_def_bi,    ONLY: grmass, grvol, philat_2d, ilon
       USE messy_main_data_bi,        ONLY: tpot_3d, geopot_3d, &
                                            geosp
       USE messy_cloud_mem,           ONLY: random_2d
       USE messy_main_constants_mem,  ONLY: pi, g
       USE messy_main_timer,          ONLY: DAYOFYEAR, YEAR, HOUR, MINUTE

       IMPLICIT NONE

       INTRINSIC :: SQRT, EXP, SIGN, LOG, ABS, SIN, MINLOC, MAXLOC

       INTEGER, INTENT(IN) :: kproma, nlev, jrow
       REAL(dp), INTENT(INOUT) :: wvel(kproma, nlev)
       INTEGER :: jp, jk, ndays, id1, id2
       REAL(dp), PARAMETER :: hh0 = 19000._dp ! reference altitude [m]
       REAL(dp), DIMENSION(kproma) :: bvfreq0
       REAL(dp), DIMENSION(kproma, nlev) :: hscale, bvfreq, hh, zrho
       REAL(dp) :: fracday, w0, AA, AA0, wness

       ! Create 2D array of random numbers (for channel output)
       DO jp = 1, kproma
          random_2d(jp,jrow) =  HARVEST(ilon(jp,jrow))
       END DO

       ! Altitude [m]
       DO jp = 1, kproma
          hh(jp,:) = geopot_3d(jp,:,jrow) + geosp(jp, jrow)
       END DO
       hh = hh / g

       ! Density scale height [m]
       zrho(1:kproma,:) = grmass(1:kproma,:,jrow) / grvol (1:kproma,:,jrow)
       hscale(1:kproma,nlev) = 0._dp
       DO jk = 1, nlev-1
          hscale(1:kproma,jk) = - hh(1:kproma,jk) / &
               LOG(zrho(1:kproma,jk)/zrho(1:kproma,nlev))
       END DO
       hscale(:,:) = MAX(hscale, 0._dp)

       ! Wintriness parameter, based on Eq. (3) Gary (2006)
       ! Note that Eq. (3) has an error, in the sin argument pi/2 has to be
       ! replaced by 2pi
       IF ((MOD(YEAR,4)==0 .AND. MOD(YEAR,100)/=0) .OR. MOD(YEAR,400)==0) THEN
          ndays = 366._dp
       ELSE
          ndays = 365._dp
       ENDIF
       fracday = DAYOFYEAR + (HOUR + MINUTE/60._dp) / 24._dp
       wness = 0.5 * (1 + &
            sin(2._dp * pi * (fracday - 295._dp) / ndays))

       ! Brunt-Vaisala frequency at all levels [1/s]
       bvfreq(:,:) = 0._dp
       DO jp = 1, kproma
          DO jk = 2, nlev - 1  ! exclude top and bottom level
             bvfreq(jp,jk) = (g / tpot_3d(jp,jk,jrow)) * &
                  ((tpot_3d(jp,jk-1,jrow) - tpot_3d(jp,jk+1,jrow)) / &
                  ((geopot_3d(jp,jk-1,jrow) - geopot_3d(jp,jk+1,jrow)) / g))
             IF (bvfreq(jp,jk).gt.0._dp) THEN
                bvfreq(jp,jk) = SQRT(bvfreq(jp,jk))
             ELSE
                bvfreq(jp,jk) = 0._dp
             ENDIF
          END DO
       END DO

       ! Interpolate Brunt-Vaisala frequency at reference level (hh0)
       DO jp = 1, kproma
          id1 = MAXLOC(hh(jp,:), DIM=1, MASK=(hh(jp,:).le.hh0))
          id2 = MINLOC(hh(jp,:), DIM=1, MASK=(hh(jp,:).gt.hh0))
          bvfreq0(jp) = (bvfreq(jp,id2) - bvfreq(jp,id1)) / &
               (hh(jp,id2) - hh(jp,id1)) * (hh0 - hh(jp,id1)) + bvfreq(jp,id1)
       END DO
       
       ! Calculate vertical velocity
       DO jp = 1, kproma

          ! Distribution at the equator at 19 km altitude
          w0 = -0.1202_dp * SIGN(1._dp,random_2d(jp,jrow)) * &
               LOG(1._dp - 2 * ABS(random_2d(jp,jrow)))  ! [m/s]

          ! Scale with latitude
          IF (philat_2d(jp,jrow).gt.0.) THEN  ! Northern Hemisphere
             AA0 = 112._dp
             AA = AA0 - 1.21_dp * philat_2d(jp,jrow) + &
                  2.20_dp * wness * philat_2d(jp,jrow)  ! Eq. (5) Gary 2008
             AA = AA / AA0
          ELSE ! Southern Hemisphere
             AA0 = 114._dp
             AA = AA0 - 0.42_dp * philat_2d(jp,jrow) + &
                  0.84 * wness * philat_2d(jp,jrow)  ! Eq. (4) Gary 2008
             AA = AA / AA0
          END IF
          w0 = w0 * AA

          ! Scale with altitude
          wvel(jp,:) = 0._dp
          DO jk = 2, nlev-1
             IF (bvfreq0(jp).eq.0._dp .OR. bvfreq(jp,jk).eq.0._dp) THEN
                wvel(jp,jk) = 0._dp
             ELSE
                wvel(jp,jk) = SQRT(bvfreq0(jp)/bvfreq(jp,jk)) * &
                     EXP((hh(jp,jk) - hh0) / (2. * hscale(jp,jk))) * w0
             END IF
          END DO
       END DO

       ! Exclude negative values
       wvel = MAX(0._dp, wvel)

       wvel = 100. * wvel  ! [cm/s]

     END SUBROUTINE calc_vervel_penner18

#endif

END MODULE messy_cloud_si

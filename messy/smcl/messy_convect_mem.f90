MODULE MESSY_CONVECT_MEM

! contains channel objects for convection channel

  USE MESSY_MAIN_CONSTANTS_MEM,   ONLY: dp

  IMPLICIT NONE
  PRIVATE
  SAVE
  INTRINSIC :: NULL
  
! Pointers of the channel objects of the convection channel
! op_mm_20140131 Added 2d_pointer
! op_mm_20140520 Added 1d_pointer
  REAL(dp), PUBLIC, DIMENSION(:,:),   POINTER :: conv_type      => NULL()
  REAL(dp), PUBLIC, DIMENSION(:,:),   POINTER :: cu_bot         => NULL()
  REAL(dp), PUBLIC, DIMENSION(:),     POINTER :: cu_bot_1d      => NULL()
  REAL(dp), PUBLIC, DIMENSION(:,:),   POINTER :: cu_top         => NULL()
  REAL(dp), PUBLIC, DIMENSION(:),     POINTER :: cu_top_1d      => NULL()
  REAL(dp), PUBLIC, DIMENSION(:,:),   POINTER :: conv_top       => NULL()
  REAL(dp), PUBLIC, DIMENSION(:),     POINTER :: conv_top_1d    => NULL()
  REAL(dp), PUBLIC, DIMENSION(:,:),   POINTER :: conv_bot       => NULL()
  REAL(dp), PUBLIC, DIMENSION(:),     POINTER :: conv_bot_1d    => NULL()
  REAL(dp), PUBLIC, DIMENSION(:,:),   POINTER :: cu_freeze      => NULL()
  REAL(dp), PUBLIC, DIMENSION(:),     POINTER :: cu_freeze_1d   => NULL()
  REAL(dp), PUBLIC, DIMENSION(:,:,:), POINTER :: cu_uvelo       => NULL()
  REAL(dp), PUBLIC, DIMENSION(:,:)  , POINTER :: cu_uvelo_2d    => NULL()
  REAL(dp), PUBLIC, DIMENSION(:,:,:), POINTER :: cv_cover       => NULL()
  REAL(dp), PUBLIC, DIMENSION(:,:,:), POINTER :: cv_cover_sikma => NULL()
  REAL(dp), PUBLIC, DIMENSION(:,:,:), POINTER :: massfu         => NULL()
  REAL(dp), PUBLIC, DIMENSION(:,:,:), POINTER :: massfd         => NULL()
  REAL(dp), PUBLIC, DIMENSION(:,:,:), POINTER :: u_entr         => NULL()
  REAL(dp), PUBLIC, DIMENSION(:,:),   POINTER :: u_entr_2d      => NULL() 
  REAL(dp), PUBLIC, DIMENSION(:,:,:), POINTER :: u_detr         => NULL()
  REAL(dp), PUBLIC, DIMENSION(:,:),   POINTER :: u_detr_2d      => NULL()
  REAL(dp), PUBLIC, DIMENSION(:,:,:), POINTER :: d_entr         => NULL()
  REAL(dp), PUBLIC, DIMENSION(:,:),   POINTER :: d_entr_2d      => NULL()
  REAL(dp), PUBLIC, DIMENSION(:,:,:), POINTER :: d_detr         => NULL()
  REAL(dp), PUBLIC, DIMENSION(:,:),   POINTER :: d_detr_2d      => NULL()
  REAL(dp), PUBLIC, DIMENSION(:,:),   POINTER :: base_f1        => NULL() 
  REAL(dp), PUBLIC, DIMENSION(:,:),   POINTER :: base_f2        => NULL() 
  REAL(dp), PUBLIC, DIMENSION(:,:),   POINTER :: base_f3        => NULL() 
  REAL(dp), PUBLIC, DIMENSION(:,:),   POINTER :: base_f4        => NULL() 
  REAL(dp), PUBLIC, DIMENSION(:),     POINTER :: base_f1_1d     => NULL() 
  REAL(dp), PUBLIC, DIMENSION(:),     POINTER :: base_f2_1d     => NULL() 
  REAL(dp), PUBLIC, DIMENSION(:),     POINTER :: base_f3_1d     => NULL() 
  REAL(dp), PUBLIC, DIMENSION(:),     POINTER :: base_f4_1d     => NULL() 
  REAL(dp), PUBLIC, DIMENSION(:,:,:), POINTER :: massfu_asc     => NULL()
  REAL(dp), PUBLIC, DIMENSION(:,:),   POINTER :: massfu_asc_2d  => NULL()
  REAL(dp), PUBLIC, DIMENSION(:,:,:), POINTER :: massfd_draf    => NULL()
  REAL(dp), PUBLIC, DIMENSION(:,:),   POINTER :: massfd_draf_2d => NULL()
  REAL(dp), PUBLIC, DIMENSION(:,:,:), POINTER :: cv_precflx     => NULL()
  REAL(dp), PUBLIC, DIMENSION(:,:),   POINTER :: cv_precflx_2d  => NULL()
  REAL(dp), PUBLIC, DIMENSION(:,:,:), POINTER :: cv_precnew     => NULL()
  REAL(dp), PUBLIC, DIMENSION(:,:),   POINTER :: cv_precnew_2d  => NULL()
  REAL(dp), PUBLIC, DIMENSION(:,:,:), POINTER :: cv_snowflx     => NULL()
  REAL(dp), PUBLIC, DIMENSION(:,:),   POINTER :: cv_snowflx_2d  => NULL()
  REAL(dp), PUBLIC, DIMENSION(:,:,:), POINTER :: cv_snownew     => NULL() 
  REAL(dp), PUBLIC, DIMENSION(:,:),   POINTER :: cv_snownew_2d  => NULL() 
  REAL(dp), PUBLIC, DIMENSION(:,:,:), POINTER :: conv_tte       => NULL()
  REAL(dp), PUBLIC, DIMENSION(:,:),   POINTER :: conv_tte_2d    => NULL()
  REAL(dp), PUBLIC, DIMENSION(:,:,:), POINTER :: conv_qte       => NULL()
  REAL(dp), PUBLIC, DIMENSION(:,:),   POINTER :: conv_qte_2d    => NULL()
  REAL(dp), PUBLIC, DIMENSION(:,:,:), POINTER :: conv_lte       => NULL() 
  REAL(dp), PUBLIC, DIMENSION(:,:,:), POINTER :: conv_ite       => NULL()
  REAL(dp), PUBLIC, DIMENSION(:,:,:), POINTER :: conv_ute       => NULL()
  REAL(dp), PUBLIC, DIMENSION(:,:,:), POINTER :: conv_vte       => NULL()
  REAL(dp), PUBLIC, DIMENSION(:,:,:), POINTER :: conv_covte     => NULL()
  ! mim_sb_20090202+
  REAL(dp), PUBLIC, DIMENSION(:,:,:), POINTER :: ttp1_gp        => NULL()
  REAL(dp), PUBLIC, DIMENSION(:,:,:), POINTER :: ptu_gp         => NULL()
  REAL(dp), PUBLIC, DIMENSION(:,:,:), POINTER :: ptd_gp         => NULL()
  REAL(dp), PUBLIC, DIMENSION(:,:,:), POINTER :: conv_tte_up    => NULL()
  REAL(dp), PUBLIC, DIMENSION(:,:),   POINTER :: conv_tte_up_2d => NULL()
  REAL(dp), PUBLIC, DIMENSION(:,:,:), POINTER :: conv_tte_do    => NULL()
  REAL(dp), PUBLIC, DIMENSION(:,:),   POINTER :: conv_tte_do_2d => NULL()
  REAL(dp), PUBLIC, DIMENSION(:,:,:), POINTER :: conv_tte_up_cond      => NULL()
  REAL(dp), PUBLIC, DIMENSION(:,:),   POINTER :: conv_tte_up_cond_2d   => NULL()
  REAL(dp), PUBLIC, DIMENSION(:,:,:), POINTER :: conv_tte_up_freeze    => NULL()
  REAL(dp), PUBLIC, DIMENSION(:,:),   POINTER :: conv_tte_up_freeze_2d => NULL()
  REAL(dp), PUBLIC, DIMENSION(:,:,:), POINTER :: conv_tte_do_verd      => NULL()
  REAL(dp), PUBLIC, DIMENSION(:,:),   POINTER :: conv_tte_do_verd_2d   => NULL()
  REAL(dp), PUBLIC, DIMENSION(:,:,:), POINTER :: conv_tte_do_melt      => NULL()
  REAL(dp), PUBLIC, DIMENSION(:,:)  , POINTER :: conv_tte_do_melt_2d   => NULL()
  REAL(dp), PUBLIC, DIMENSION(:,:,:), POINTER :: conv_tte_su    => NULL()
  REAL(dp), PUBLIC, DIMENSION(:,:),   POINTER :: conv_tte_su_2d => NULL()
  REAL(dp), PUBLIC, DIMENSION(:,:,:), POINTER :: conv_tte_ev    => NULL()
  REAL(dp), PUBLIC, DIMENSION(:,:),   POINTER :: conv_tte_ev_2d => NULL()
  REAL(dp), PUBLIC, DIMENSION(:,:,:), POINTER :: conv_qte_up    => NULL()
  REAL(dp), PUBLIC, DIMENSION(:,:,:), POINTER :: conv_qte_do    => NULL()
  REAL(dp), PUBLIC, DIMENSION(:,:,:), POINTER :: conv_qte_su    => NULL() 
  REAL(dp), PUBLIC, DIMENSION(:,:,:), POINTER :: conv_qte_ev    => NULL()
  REAL(dp), PUBLIC, DIMENSION(:,:,:), POINTER :: conv_pqtec     => NULL()
  REAL(dp), PUBLIC, DIMENSION(:,:,:), POINTER :: conv_pxtec     => NULL()
  REAL(dp), PUBLIC, DIMENSION(:,:),   POINTER :: aprsc          => NULL() 
  REAL(dp), PUBLIC, DIMENSION(:,:),   POINTER :: aprss          => NULL() 
  ! op_mm_20140226+
  REAL(dp), PUBLIC, DIMENSION(:,:),   POINTER :: vgustcon       => NULL() 
  ! op_mm_20140226-



  ! ... for CPL namelist
  LOGICAL, PUBLIC :: l_lgmc_diag = .FALSE.  ! switch for additional diagnostic (LG-moist convection)
  ! mim_sb_20090202-

! mz_pj_20050615+
  REAL(dp), PUBLIC, DIMENSION(:,:),   POINTER :: cu_bot_mid    => NULL()
  REAL(dp), PUBLIC, DIMENSION(:),     POINTER :: cu_bot_mid_1d => NULL()
  REAL(dp), PUBLIC, DIMENSION(:,:),   POINTER :: cu_top_mid    => NULL()
  REAL(dp), PUBLIC, DIMENSION(:),     POINTER :: cu_top_mid_1d => NULL()
  REAL(dp), PUBLIC, DIMENSION(:,:),   POINTER :: cu_freeze_mid => NULL()
  REAL(dp), PUBLIC, DIMENSION(:),     POINTER :: cu_freeze_mid_1d=> NULL()
  ! mz_pj_20050615-
  REAL(dp), PUBLIC, DIMENSION(:,:),   POINTER :: WAT_DIAG    => NULL()
  REAL(dp), PUBLIC, DIMENSION(:,:),   POINTER :: CAPE        => NULL()
  REAL(dp), PUBLIC, DIMENSION(:,:,:), POINTER :: udetr_h     => NULL()
  REAL(dp), PUBLIC, DIMENSION(:,:)  , POINTER :: udetr_h_2d  => NULL()
  REAL(dp), PUBLIC, DIMENSION(:,:,:), POINTER :: cv_lwc      => NULL()
  REAL(dp), PUBLIC, DIMENSION(:,:),   POINTER :: cv_lwc_2d   => NULL()
  REAL(dp), PUBLIC, DIMENSION(:,:,:), POINTER :: cv_iwc      => NULL()
  REAL(dp), PUBLIC, DIMENSION(:,:),   POINTER :: cv_iwc_2d   => NULL()

  REAL(dp), PUBLIC, DIMENSION(:,:,:), POINTER :: cv_rform    => NULL()
  REAL(dp), PUBLIC, DIMENSION(:,:)  , POINTER :: cv_rform_2d => NULL()
  REAL(dp), PUBLIC, DIMENSION(:,:,:), POINTER :: cv_sform    => NULL()
  REAL(dp), PUBLIC, DIMENSION(:,:),   POINTER :: cv_sform_2d => NULL()
  REAL(dp), PUBLIC, DIMENSION(:,:,:), POINTER :: tketconv    => NULL() ! op_mm_20140227

  REAL(dp), PUBLIC, DIMENSION(:,:,:), POINTER :: c_cldfrac   => NULL()
  REAL(dp), PUBLIC, DIMENSION(:,:,:), POINTER :: c_liqamt    => NULL()
  REAL(dp), PUBLIC, DIMENSION(:,:,:), POINTER :: c_liqsize   => NULL()
  REAL(dp), PUBLIC, DIMENSION(:,:,:), POINTER :: c_iceamt    => NULL()
  REAL(dp), PUBLIC, DIMENSION(:,:,:), POINTER :: c_icesize   => NULL()
  REAL(dp), PUBLIC, DIMENSION(:,:,:), POINTER :: m_cldfrac   => NULL()
  REAL(dp), PUBLIC, DIMENSION(:,:,:), POINTER :: m_liqamt    => NULL()
  REAL(dp), PUBLIC, DIMENSION(:,:,:), POINTER :: m_liqsize   => NULL()
  REAL(dp), PUBLIC, DIMENSION(:,:,:), POINTER :: m_iceamt    => NULL()
  REAL(dp), PUBLIC, DIMENSION(:,:,:), POINTER :: m_icesize   => NULL()
  REAL(dp), PUBLIC, DIMENSION(:,:,:), POINTER :: d_humarea   => NULL()
  REAL(dp), PUBLIC, DIMENSION(:,:,:), POINTER :: d_humratio  => NULL()

  ! mz_ak_20051221+
  ! liquid water content of convective cloud
  REAL(dp), PUBLIC, DIMENSION(:,:,:), POINTER :: cv_cldwater    => NULL()
  REAL(dp), PUBLIC, DIMENSION(:,:),   POINTER :: cv_cldwater_2d => NULL()
  ! change in liquid water
  REAL(dp), PUBLIC, DIMENSION(:,:,:), POINTER :: del_liqwat     => NULL()
  REAL(dp), PUBLIC, DIMENSION(:,:),   POINTER :: del_liqwat_2d  => NULL()
  ! mz_ak_20051221-

  LOGICAL,  PUBLIC :: ODEEP
  LOGICAL,  PUBLIC :: OSHAL
  LOGICAL,  PUBLIC :: ODOWN
  LOGICAL,  PUBLIC :: OREFRESH_ALL
  LOGICAL,  PUBLIC :: OSETTADJ
  LOGICAL,  PUBLIC :: OUVTRANS
  LOGICAL,  PUBLIC :: OCHTRANS

  INTEGER,  PUBLIC :: KENSM
  INTEGER,  PUBLIC :: KICE
  REAL(dp), PUBLIC :: PTADJD
  REAL(dp), PUBLIC :: PTADJS

! local pointer for information of other submodels
  LOGICAL,  PUBLIC :: cvtrans                           ! information use_cvtrans

END MODULE MESSY_CONVECT_MEM

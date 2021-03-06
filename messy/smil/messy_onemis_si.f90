#include "messy_main_ppd_bi.inc"

!*****************************************************************************
! Authors:
!   of ONLEM
!   Rolf Sander,     MPICH, 2004
!   Astrid Kerkweg,  MPICH, 2004, 2008, 2012
!   Patrick Joeckel, MPICH, 2004
!   Joerg Steinkamp, MPICH, 2008, improved version of Yienger and Lewy, 1995
!                                 (NO_yl95sl10)
!   Susannah Burrows, MPICH, 2009, expanded by bioaerosol emissions
!   Gregor Glaeser,   Uni-Mz, 2010, Dust emission from Tegen et al.
!   Marina Astitha,   CyI-EEWRC, Dust emissions scheme new implementation
!                     (Astitha et al. 2012)
!
!   of ONEMIS:
!    Astrid Kerkweg, Uni-Mz, Mar 2010 ONLEM rewritten to ONEMIS
!                                     INPUT is done by coupling via CHANNEL
!                                     no more 'self-regridding'
!                            2011     adapted to COSMO
!   Stefanie Falk,   KIT,    2017, Bromine emission scheme from snow and sea ice
MODULE messy_onemis_si

  ! ECHAM5/MESSy
  USE messy_main_blather_bi,    ONLY: start_message_bi, end_message_bi &
                                    , error_bi, warning_bi
#ifdef MESSYTENDENCY
 !tendency budget
 USE messy_main_tendency_bi,    ONLY: mtend_get_handle,       &
                                      mtend_get_start_l,      &
                                      mtend_add_l,            &
                                      mtend_register,         &
                                      mtend_id_tracer
#endif
  ! MESSy
  USE messy_main_constants_mem, ONLY: dp, STRLEN_MEDIUM
  USE messy_main_tools,         ONLY: PTR_3D_ARRAY, PTR_2D_ARRAY
  USE messy_main_channel,       ONLY: STRLEN_OBJECT, STRLEN_CHANNEL &
                                    , t_chaobj_cpl
  USE messy_onemis

    ! um_gg_20130527+
#ifdef ECHAM5
  USE messy_main_grid_def_mem_bi,   ONLY: nlon
#endif
#ifdef COSMO
  USE messy_main_grid_def_mem_bi,   ONLY: dlon, dlat
#endif
    ! um_gg_20130527-

  IMPLICIT NONE
  PRIVATE
  SAVE

  INTRINSIC :: ADJUSTL, ASSOCIATED, EXP, MAX, NULL, SIZE, SQRT, TRIM

  ! CONVERT FLUXES TO TRACERS ------------------------------------------
  INTEGER, PARAMETER :: MAXPSTRLEN = 500
  TYPE T_F2T_IO
     CHARACTER(LEN=STRLEN_OBJECT)    :: name      = ''
     CHARACTER(LEN=MAXPSTRLEN)       :: string_gp = ''
     CHARACTER(LEN=MAXPSTRLEN)       :: string_lg = ''
  END TYPE T_F2T_IO
  !
  TYPE T_F2T
     CHARACTER(LEN=STRLEN_OBJECT)    :: name = ''
     ! POINTER TO CHANNEL OBJECT (FLUX)
     REAL(DP),                   DIMENSION(:,:,:), POINTER :: ptr => NULL()
     !
     ! GRIDPOINT
     CHARACTER(LEN=2*STRLEN_MEDIUM+1), DIMENSION(:), POINTER :: &
          tracer_gp => NULL()
     INTEGER                                             :: ngpt = 0
     INTEGER,                      DIMENSION(:), POINTER :: mgp => NULL()
     REAL(DP),                     DIMENSION(:), POINTER :: efact_gp => NULL()
     INTEGER,                      DIMENSION(:), POINTER :: idt_gp => NULL()
     !
     ! LAGRANGIAN
     CHARACTER(LEN=2*STRLEN_MEDIUM+1), DIMENSION(:), POINTER :: &
          tracer_lg => NULL()
     INTEGER                                             :: nlgt = 0
     INTEGER,                      DIMENSION(:), POINTER :: mlg => NULL()
     REAL(DP),                     DIMENSION(:), POINTER :: efact_lg => NULL()
     INTEGER,                      DIMENSION(:), POINTER :: idt_lg => NULL()
     ! REST FLUX FOR LAGRANGIAN (METHOD LG=1 ONLY !)
     TYPE(PTR_3D_ARRAY),           DIMENSION(:), POINTER :: ptr_rest => NULL()
     !
  END TYPE T_F2T
  !
  INTEGER, PARAMETER                          :: NMAXFLUXES = 100
  TYPE(T_F2T_IO), DIMENSION(NMAXFLUXES)       :: F2T   ! INPUT FROM namelist
  TYPE(T_F2T),    DIMENSION(NMAXFLUXES)       :: XF2T  ! WORKSPACE
  INTEGER                                     :: NFLUXES = 0
  ! --------------------------------------------------------------------

  ! pointers to external data
  REAL(dp), DIMENSION(:,:), POINTER :: wind10_2d => NULL()

  ! channel objects:
  ! dms emissions
  REAL(dp), DIMENSION(:,:), POINTER :: emis_dms_sea
  REAL(dp), DIMENSION(:,:), POINTER :: seawater_dms => NULL()
  TYPE(T_CHAOBJ_CPL)                :: imp_seawater_dms

  ! sea salt emissions mass (mss) and number(nss) fluxes by M. Schulz (lsce)
  ! directly distributed to a accumulation (as) and a coarse (cs) mode
  REAL(dp), DIMENSION(:,:), POINTER :: mss_as_lsce, mss_cs_lsce
  REAL(dp), DIMENSION(:,:), POINTER :: nss_as_lsce, nss_cs_lsce
  ! sea salt emissions mass (mss) and number(nss) fluxes by Monahan (monahan)
  ! directly distributed to a accumulation (as) and a coarse (cs) mode
  REAL(dp), DIMENSION(:,:), POINTER :: mss_as_monahan, mss_cs_monahan
  REAL(dp), DIMENSION(:,:), POINTER :: nss_as_monahan, nss_cs_monahan

  ! sea salt emissions mass (mss) and number(nss) fluxes defined by aerocom
  ! (aerocom)
  ! directly distributed to a accumulation (as) and a coarse (cs) mode
  REAL(dp), DIMENSION(:,:), POINTER :: mss_as_aerocom, mss_cs_aerocom
  REAL(dp), DIMENSION(:,:), POINTER :: nss_as_aerocom, nss_cs_aerocom
  ! regridded emissions as defined aerocom (unit: per gridbox and per day)
  REAL(dp), DIMENSION(:,:), POINTER :: numflx_as_aerocom => NULL()
  REAL(dp), DIMENSION(:,:), POINTER :: numflx_cs_aerocom => NULL()
  TYPE(T_CHAOBJ_CPL) :: imp_numflx_as_aerocom, imp_numflx_cs_aerocom
  REAL(dp), DIMENSION(:,:), POINTER :: massflx_as_aerocom => NULL()
  REAL(dp), DIMENSION(:,:), POINTER :: massflx_cs_aerocom => NULL()
  TYPE(T_CHAOBJ_CPL) :: imp_massflx_as_aerocom, imp_massflx_cs_aerocom

  ! organic contribution of seasalt
  ! originally implemented by Susannah Burrows
  ! using the SS flux calculated or used in onemis
  ! Sea WIOC emissions
  ! channel objects
  TYPE SS_OC_EMIS
    REAL(dp), DIMENSION(:,:), POINTER :: emis_flx    => NULL()
    REAL(dp), DIMENSION(:,:), POINTER :: emis_flxsum => NULL()
    REAL(dp), DIMENSION(:,:), POINTER :: ss_flx      => NULL()
    CHARACTER(LEN=STRLEN_CHANNEL)     :: channel
    CHARACTER(LEN=STRLEN_OBJECT)      :: object
  END TYPE SS_OC_EMIS
  INTEGER :: n_ss_poc_aqua, n_ss_poc_swifs, n_ss_wioc_aqua, n_ss_wioc_blend
  TYPE(SS_OC_EMIS), DIMENSION(:), POINTER :: emis_poc_aqua       => NULL()
  TYPE(SS_OC_EMIS), DIMENSION(:), POINTER :: emis_poc_swifs      => NULL()
  TYPE(SS_OC_EMIS), DIMENSION(:), POINTER :: emis_wioc_aqua      => NULL()
  TYPE(SS_OC_EMIS), DIMENSION(:), POINTER :: emis_wioc_blend     => NULL()
  ! working pointers
  REAL(dp), DIMENSION(:), POINTER :: emis_poc       => NULL()
  REAL(dp), DIMENSION(:), POINTER :: emisflxsum_poc => NULL()
  REAL(dp), DIMENSION(:), POINTER :: ss_flux        => NULL()
  ! seawater concentrations
  REAL(dp), DIMENSION(:,:), POINTER :: poc_aqua           => NULL()
  REAL(dp), DIMENSION(:,:), POINTER :: poc_seawifs        => NULL()
  REAL(dp), DIMENSION(:,:), POINTER :: chlor_a_aqua_wioc  => NULL()
  REAL(dp), DIMENSION(:,:), POINTER :: chlor_a_blend_wioc => NULL()
  TYPE(T_CHAOBJ_CPL) :: imp_poc_aqua,  imp_poc_swifs
  TYPE(T_CHAOBJ_CPL) :: imp_wioc_aqua, imp_wioc_blen
  TYPE(T_CHAOBJ_CPL) :: imp_ss_poc_aqua,  imp_ss_poc_swifs
  TYPE(T_CHAOBJ_CPL) :: imp_ss_wioc_aqua, imp_ss_wioc_blen

  ! diagnostic output
  REAL(dp), DIMENSION(:,:), POINTER :: emisflx_ss_cs_sum => NULL()
  REAL(dp), DIMENSION(:,:), POINTER :: emisflx_ss_as_sum => NULL()

  ! carbon emission calculations
  REAL(dp), DIMENSION(:,:), POINTER :: OC_ag  => NULL()
  REAL(dp), DIMENSION(:,:), POINTER :: OC_ant => NULL()
  REAL(dp), DIMENSION(:,:), POINTER :: OC_bge => NULL()
  REAL(dp), DIMENSION(:,:), POINTER :: OC_wf  => NULL()
  REAL(dp), DIMENSION(:,:), POINTER :: BC_ag  => NULL()
  REAL(dp), DIMENSION(:,:), POINTER :: BC_ant => NULL()
  REAL(dp), DIMENSION(:,:), POINTER :: BC_wf  => NULL()
  TYPE(T_CHAOBJ_CPL) :: imp_OC_ag
  TYPE(T_CHAOBJ_CPL) :: imp_OC_ant
  TYPE(T_CHAOBJ_CPL) :: imp_OC_bge
  TYPE(T_CHAOBJ_CPL) :: imp_OC_wf
  TYPE(T_CHAOBJ_CPL) :: imp_BC_ag
  TYPE(T_CHAOBJ_CPL) :: imp_BC_ant
  TYPE(T_CHAOBJ_CPL) :: imp_BC_wf

  ! carbon emission fluxes
  ! orcanic carbon insoluble
  REAL(dp), DIMENSION(:,:), POINTER :: OC_sum_sol   => NULL()
  ! orcanic carbon soluble
  REAL(dp), DIMENSION(:,:), POINTER :: OC_sum_insol => NULL()
  ! black carbon insoluble
  REAL(dp), DIMENSION(:,:), POINTER :: BC_sum_insol => NULL()
  ! number carbon insoluble (still to be multiplied with mode information
  ! in the respective aerosol module (see e.g. M7))
  REAL(dp), DIMENSION(:,:), POINTER :: N_insol      => NULL()
  ! number carbon soluble (still to be multiplied with mode information
  ! in the respective aerosol module (see e.g. M7))
  REAL(dp), DIMENSION(:,:), POINTER :: N_sol        => NULL()

  ! mz_ht_20110301+
  REAL(dp), DIMENSION(:,:), POINTER :: OC_soa_sol   => NULL()
  REAL(dp), DIMENSION(:,:), POINTER :: OC_bb_sol    => NULL()
  REAL(dp), DIMENSION(:,:), POINTER :: OC_ff_sol    => NULL()
  ! orcanic carbon soluble
  REAL(dp), DIMENSION(:,:), POINTER :: OC_soa_insol => NULL()
  REAL(dp), DIMENSION(:,:), POINTER :: OC_bb_insol  => NULL()
  REAL(dp), DIMENSION(:,:), POINTER :: OC_ff_insol  => NULL()
  ! black carbon insoluble
  REAL(dp), DIMENSION(:,:), POINTER :: BC_bb_insol  => NULL()
  REAL(dp), DIMENSION(:,:), POINTER :: BC_ff_insol  => NULL()
  REAL(dp), DIMENSION(:,:), POINTER :: N_insol_bc   => NULL()
  REAL(dp), DIMENSION(:,:), POINTER :: N_insol_oc   => NULL()
  ! mz_ht_20110301-

  ! O3 emisflux over ice
  REAL(dp), DIMENSION(:,:), POINTER :: O3_emflux => NULL()

  ! CH4 emisflux
  REAL(dp), DIMENSION(:,:), POINTER :: CH4_emflux    => NULL()
  ! climatological value
  REAL(dp), DIMENSION(:,:), POINTER :: CH4_conc_clim => NULL()
  TYPE(T_CHAOBJ_CPL)                :: imp_CH4_conc_clim
  ! special: id of associated tracer
  INTEGER :: idt_CH4 = 0

  !VOC emisflux etc.
  ! cos zenith angle
  LOGICAL                           :: l_cossza = .FALSE.
  REAL(dp), DIMENSION(:,:), POINTER :: cossza_2d     => NULL()
  TYPE(T_CHAOBJ_CPL)                :: imp_cossza
  REAL(dp), DIMENSION(:,:), POINTER :: drymatter     => NULL() ! dry matter
  TYPE(T_CHAOBJ_CPL)                :: imp_drymatter           ! dry matter
  REAL(dp), DIMENSION(:,:), POINTER :: isop_emflux   => NULL() ! emis. fluxes
  REAL(dp), DIMENSION(:,:), POINTER :: isop_emisfac  => NULL() ! emis. fluxes
  TYPE(T_CHAOBJ_CPL)                :: imp_isop_emisfac        ! emis. fluxes
  REAL(dp), DIMENSION(:,:), POINTER :: mterp_emflux  => NULL() ! emis. fluxes
  REAL(dp), DIMENSION(:,:), POINTER :: mterp_emisfac => NULL() ! emis. fluxes
  TYPE(T_CHAOBJ_CPL)                :: imp_mterp_emisfac       ! emis. fluxes

  ! channel object pointer; to be read in for NO and VOC emissions
  TYPE(T_CHAOBJ_CPL)   :: imp_lai         ! LAI
  TYPE(T_CHAOBJ_CPL)   :: imp_lad_top     ! LAD profile
  TYPE(T_CHAOBJ_CPL)   :: imp_lad_soil    ! LAD profile
  TYPE(T_CHAOBJ_CPL)   :: imp_lad_lay2    ! LAD profile
  TYPE(T_CHAOBJ_CPL)   :: imp_lad_lay3    ! LAD profile
  TYPE(T_CHAOBJ_CPL)   :: imp_hc          ! canopy height
  TYPE(T_CHAOBJ_CPL)   :: imp_drag        !  drag coeff.
  TYPE(T_CHAOBJ_CPL)   :: imp_disp        ! displacement height
  TYPE(T_CHAOBJ_CPL)   :: imp_forestfr    ! forest fraction
  REAL(dp), DIMENSION(:,:), POINTER :: lai         => NULL()! LAI
  REAL(dp), DIMENSION(:,:), POINTER :: lad_top     => NULL()! LAD profile
  REAL(dp), DIMENSION(:,:), POINTER :: lad_soil    => NULL()! LAD profile
  REAL(dp), DIMENSION(:,:), POINTER :: lad_lay2    => NULL()! LAD profile
  REAL(dp), DIMENSION(:,:), POINTER :: lad_lay3    => NULL()! LAD profile
  REAL(dp), DIMENSION(:,:), POINTER :: hc          => NULL()! canopy height
  REAL(dp), DIMENSION(:,:), POINTER :: drag        => NULL()!  drag coeff.
  REAL(dp), DIMENSION(:,:), POINTER :: disp        => NULL()! displacement height
  REAL(dp), DIMENSION(:,:), POINTER :: forestfr    => NULL()! forest fraction
  REAL(dp), DIMENSION(:,:), POINTER :: z0m         => NULL()! surface roughness
  REAL(dp) ,DIMENSION(:,:,:), POINTER :: lad  => NULL()  ! LAD profile

  ! NO emisfluxes etc.
  ! ratio veg./emis class 1
  REAL(dp), DIMENSION(:,:,:), POINTER  :: NOemisclass1 => NULL()
  TYPE(T_CHAOBJ_CPL)                   :: imp_NOemisclass1
  ! ratio veg./emis class 1
  REAL(dp), DIMENSION(:,:,:), POINTER  :: NOemisclass2 => NULL()
  TYPE(T_CHAOBJ_CPL)                   :: imp_NOemisclass2
  ! cultivation intensity
  REAL(dp), DIMENSION(:,:),   POINTER  :: cultiv  => NULL()
  TYPE(T_CHAOBJ_CPL)                   :: imp_cultiv
  ! fertilizer application
  REAL(dp), DIMENSION(:,:),   POINTER  :: fertil  => NULL()
  TYPE(T_CHAOBJ_CPL)                   :: imp_fertil
  ! NO emission factor for wet soils
  REAL(dp), DIMENSION(:,:),   POINTER  :: noemis_w     => NULL()
  ! NO emission factor for dry soils
  REAL(dp), DIMENSION(:,:),   POINTER  :: noemis_d     => NULL()
  ! monthly mean total precipitation
  REAL(dp), DIMENSION(:,:),   POINTER  :: prectot      => NULL()

  !   channel object pointer for NO pulsing
  ! conv. prec. record previous month
  REAL(dp), DIMENSION(:,:,:), POINTER :: cpold      => NULL()
  ! large scale prec. record previous month
  REAL(dp), DIMENSION(:,:,:), POINTER :: lspold     => NULL()
  ! conv. prec.
  REAL(dp), DIMENSION(:,:),   POINTER :: cp         => NULL()
  ! large scale prec.
  REAL(dp), DIMENSION(:,:),   POINTER :: lsp        => NULL()
  ! status of the pulsing event
  REAL(dp), DIMENSION(:,:),   POINTER :: plsday     => NULL()
  ! pulsing regime
  REAL(dp), DIMENSION(:,:),   POINTER :: pulsing    => NULL()
  ! duration of pulse
  REAL(dp), DIMENSION(:,:),   POINTER :: plsdurat   => NULL()
  ! soil-biogenic NO emission fluxes
  REAL(dp), DIMENSION(:,:),   POINTER :: NO_emflux  => NULL()
  ! pulse of the soil emission flux
  REAL(dp), DIMENSION(:,:),   POINTER :: pls        => NULL()

  REAL(dp), DIMENSION(:,:),   POINTER :: noslflux      => NULL()
  REAL(dp), DIMENSION(:,:,:), POINTER :: noslflux_diag => NULL()

  ! added variables for the orignal Yienger & Levy algorithm
  REAL(dp), DIMENSION(:,:),   POINTER :: pulse            => NULL()
  REAL(dp), DIMENSION(:,:),   POINTER :: pulseday         => NULL()
  REAL(dp), DIMENSION(:,:),   POINTER :: pulseregime      => NULL()
  REAL(dp), DIMENSION(:,:,:), POINTER :: prec_hist        => NULL()
  REAL(dp), DIMENSION(:,:),   POINTER :: crf              => NULL()
  TYPE(T_CHAOBJ_CPL)                  :: imp_lai_yl95sl10
  REAL(dp), DIMENSION(:,:),   POINTER :: lai_yl95sl10     => NULL()
  TYPE(T_CHAOBJ_CPL)                  :: imp_fertil_yl95sl10
  REAL(dp), DIMENSION(:,:),   POINTER :: fertil_yl95sl10  => NULL()
  TYPE(T_CHAOBJ_CPL)                  :: imp_NOemclass_yl95sl10
  REAL(dp), DIMENSION(:,:),   POINTER :: noslflux_yl95sl10  => NULL()
  REAL(dp), DIMENSION(:,:),   POINTER :: noemflux_yl95sl10  => NULL()
  REAL(dp), DIMENSION(:,:,:), POINTER :: NOemclass_yl95sl10 => NULL()
  REAL(dp), DIMENSION(:,:,:), POINTER :: noslflux_diag_yl95sl10 => NULL()
  REAL(dp), DIMENSION(:,:),   POINTER :: tsoil_top          => NULL()
  REAL(dp), DIMENSION(:,:),   POINTER :: vsm                => NULL()
  TYPE(T_CHAOBJ_CPL)                  :: imp_rootdepth_yl95sl10
  REAL(dp), DIMENSION(:,:),   POINTER :: root_depth         => NULL()
  TYPE(T_CHAOBJ_CPL)                  :: imp_rootdepth_mask_yl95sl10
  REAL(dp), DIMENSION(:,:),   POINTER :: root_depth_mask    => NULL()

  ! dust emission for Schulz scheme
  REAL(dp), DIMENSION(:,:),   POINTER :: du_cla           => NULL()
  TYPE(T_CHAOBJ_CPL)                  :: imp_du_cla
  REAL(dp), DIMENSION(:,:),   POINTER :: du_src           => NULL()
  TYPE(T_CHAOBJ_CPL)                  :: imp_du_src
  REAL(dp), DIMENSION(:,:),   POINTER :: du_thr           => NULL()
  TYPE(T_CHAOBJ_CPL)                  :: imp_du_thr
  REAL(dp), DIMENSION(:,:),   POINTER :: du_emflux_B_ci     => NULL()
  REAL(dp), DIMENSION(:,:),   POINTER :: zprecipinsoil_2d => NULL()

  REAL(dp), DIMENSION(:,:),     POINTER :: du_cla2     => NULL()
  TYPE(T_CHAOBJ_CPL)                    :: imp_du_cla2
  REAL(dp), DIMENSION(:,:),     POINTER :: rdepth      => NULL()
  TYPE(T_CHAOBJ_CPL)                    :: imp_rdepth
  REAL(dp), DIMENSION(:,:,:,:), POINTER :: dustsrc     => NULL()
  TYPE(T_CHAOBJ_CPL)                    :: imp_dustsrc
  REAL(dp), DIMENSION(:,:),     POINTER :: lai_in      => NULL()
  TYPE(T_CHAOBJ_CPL)                    :: imp_lai_in
  REAL(dp), DIMENSION(:,:,:,:), POINTER :: soiltexture => NULL()
  TYPE(T_CHAOBJ_CPL)                    :: imp_soiltexture

  REAL(dp), DIMENSION(:,:),     POINTER :: du_nap     => NULL()
  TYPE(T_CHAOBJ_CPL)                    :: imp_du_nap
  REAL(dp), DIMENSION(:,:),     POINTER :: du_kp      => NULL()
  TYPE(T_CHAOBJ_CPL)                    :: imp_du_kp
  REAL(dp), DIMENSION(:,:),     POINTER :: du_capp     => NULL()
  TYPE(T_CHAOBJ_CPL)                    :: imp_du_capp
  REAL(dp), DIMENSION(:,:),     POINTER :: du_mgpp     => NULL()
  TYPE(T_CHAOBJ_CPL)                    :: imp_du_mgpp
  REAL(dp), DIMENSION(:,:),     POINTER :: du_misc     => NULL()
  TYPE(T_CHAOBJ_CPL)                    :: imp_du_misc
  TYPE(T_CHAOBJ_CPL) :: imp_kkdu_clay
  TYPE(T_CHAOBJ_CPL) :: imp_kkdu_mask
  TYPE(T_CHAOBJ_CPL) :: imp_kkdu_lai
  TYPE(T_CHAOBJ_CPL) :: imp_kkdu_topo
  TYPE(T_CHAOBJ_CPL) :: imp_kkdu_nap
  TYPE(T_CHAOBJ_CPL) :: imp_kkdu_kp
  TYPE(T_CHAOBJ_CPL) :: imp_kkdu_capp
  TYPE(T_CHAOBJ_CPL) :: imp_kkdu_mgpp
  TYPE(T_CHAOBJ_CPL) :: imp_kkdu_misc
  REAL(dp), DIMENSION(:,:), POINTER :: kkdu_clay => NULL()
  REAL(dp), DIMENSION(:,:,:), POINTER :: kkdu_mask => NULL()
  REAL(dp), DIMENSION(:,:), POINTER :: kkdu_lai => NULL()
  REAL(dp), DIMENSION(:,:), POINTER :: kkdu_topo => NULL()
  REAL(dp), DIMENSION(:,:), POINTER :: kkdu_nap => NULL()
  REAL(dp), DIMENSION(:,:), POINTER :: kkdu_kp => NULL()
  REAL(dp), DIMENSION(:,:), POINTER :: kkdu_capp => NULL()
  REAL(dp), DIMENSION(:,:), POINTER :: kkdu_mgpp => NULL()
  REAL(dp), DIMENSION(:,:), POINTER :: kkdu_misc => NULL()
  REAL(dp), DIMENSION(:,:), POINTER :: kkdu_nap_emflux_ai => NULL()
  REAL(dp), DIMENSION(:,:), POINTER :: kkdu_nap_emflux_ci => NULL()
  REAL(dp), DIMENSION(:,:), POINTER :: kkdu_kp_emflux_ai => NULL()
  REAL(dp), DIMENSION(:,:), POINTER :: kkdu_kp_emflux_ci => NULL()
  REAL(dp), DIMENSION(:,:), POINTER :: kkdu_capp_emflux_ai => NULL()
  REAL(dp), DIMENSION(:,:), POINTER :: kkdu_capp_emflux_ci => NULL()
  REAL(dp), DIMENSION(:,:), POINTER :: kkdu_mgpp_emflux_ai => NULL()
  REAL(dp), DIMENSION(:,:), POINTER :: kkdu_mgpp_emflux_ci => NULL()
  REAL(dp), DIMENSION(:,:), POINTER :: kkdu_misc_emflux_ai => NULL()
  REAL(dp), DIMENSION(:,:), POINTER :: kkdu_misc_emflux_ci => NULL()

  REAL(dp), DIMENSION(:,:),     POINTER :: du_emflux_A1_ai => NULL()
  REAL(dp), DIMENSION(:,:),     POINTER :: du_emflux_A1_ci => NULL()
  REAL(dp), DIMENSION(:,:),     POINTER :: du_emflux_A2_ai => NULL()
  REAL(dp), DIMENSION(:,:),     POINTER :: du_emflux_A2_ci => NULL()
  REAL(dp), DIMENSION(:,:),     POINTER :: horflux         => NULL()
  REAL(dp), DIMENSION(:,:),     POINTER :: ustarthr1       => NULL()
  REAL(dp), DIMENSION(:,:,:),   POINTER :: ustarthr        => NULL()

  REAL(dp), DIMENSION(:,:),     POINTER :: du_nap_emflux_ai => NULL()
  REAL(dp), DIMENSION(:,:),     POINTER :: du_nap_emflux_ci => NULL()
  REAL(dp), DIMENSION(:,:),     POINTER :: du_kp_emflux_ai => NULL()
  REAL(dp), DIMENSION(:,:),     POINTER :: du_kp_emflux_ci => NULL()
  REAL(dp), DIMENSION(:,:),     POINTER :: du_capp_emflux_ai => NULL()
  REAL(dp), DIMENSION(:,:),     POINTER :: du_capp_emflux_ci => NULL()
  REAL(dp), DIMENSION(:,:),     POINTER :: du_mgpp_emflux_ai => NULL()
  REAL(dp), DIMENSION(:,:),     POINTER :: du_mgpp_emflux_ci => NULL()
  REAL(dp), DIMENSION(:,:),     POINTER :: du_misc_emflux_ai => NULL()
  REAL(dp), DIMENSION(:,:),     POINTER :: du_misc_emflux_ci => NULL()

  REAL(dp), DIMENSION(:,:),  POINTER :: mat_s2         => NULL()
  TYPE(T_CHAOBJ_CPL)                 :: imp_mat_s2
  REAL(dp), DIMENSION(:,:),  POINTER :: mat_s3         => NULL()
  TYPE(T_CHAOBJ_CPL)                 :: imp_mat_s3
  REAL(dp), DIMENSION(:,:),  POINTER :: mat_s4         => NULL()
  TYPE(T_CHAOBJ_CPL)                 :: imp_mat_s4
  REAL(dp), DIMENSION(:,:),  POINTER :: mat_s6         => NULL()
  TYPE(T_CHAOBJ_CPL)                 :: imp_mat_s6
  REAL(dp), DIMENSION(:,:),  POINTER :: mat_psrc       => NULL()
  TYPE(T_CHAOBJ_CPL)                 :: imp_mat_psrc
  REAL(dp), DIMENSION(:,:),  POINTER :: k_fpar_eff     => NULL()
  TYPE(T_CHAOBJ_CPL)                 :: imp_k_fpar_eff
  REAL(dp), DIMENSION(:,:),  POINTER :: du_emflux_T_ai => NULL()
  REAL(dp), DIMENSION(:,:),  POINTER :: du_emflux_T_ci => NULL()
  REAL(dp)                           :: cuscale ! um_gg_20130527

  ! SO2 and/or sulphate aersosol emission
  ! input fields
  REAL(dp), DIMENSION(:,:),  POINTER  :: SO2_ant_high  => NULL()
  TYPE(T_CHAOBJ_CPL)                  :: imp_SO2_ant_high
  REAL(dp), DIMENSION(:,:),  POINTER  :: SO2_ant_low   => NULL()
  TYPE(T_CHAOBJ_CPL)                  :: imp_SO2_ant_low
  REAL(dp), DIMENSION(:,:,:), POINTER :: SO2_ant        => NULL()
  ! output fields
  REAL(dp), DIMENSION(:,:,:), POINTER :: SO2_emflux => NULL()

  ! bacteria -- primary aerosol emission
  ! OLSON LUMPED EMISSIONS
  INTEGER, PARAMETER :: nolclass = 11
  CHARACTER(LEN=21), DIMENSION(nolclass)  :: olson_name = (/    &
       'olson_emis_seas      ','olson_emis_landice   ' &
      ,'olson_emis_deserts   ','olson_emis_forests   ' &
      ,'olson_emis_grasslands','olson_emis_crops     ' &
      ,'olson_emis_wetlands  ','olson_emis_shrubs    ' &
      ,'olson_emis_coastal   ','olson_emis_urban     ' &
      ,'olson_emis_tundra    ' /)

  TYPE(PTR_2D_ARRAY), DIMENSION(nolclass) :: olson_emis
  REAL(dp), DIMENSION(:,:,:), POINTER     :: olson       => NULL()
  TYPE(T_CHAOBJ_CPL)                      :: imp_olson

  ! MODIS EMISSIONS
  INTEGER, PARAMETER :: nmodisclass = 18 ! Number of classes in MODIS dataset.
  CHARACTER(LEN=15), DIMENSION(nmodisclass) :: modis_name = (/  &
       'water          ', 'ever_need      ', 'ever_broad     ', &
       'deci_need      ', 'deci_broad     ', 'mixed_forest   ', &
       'closed_shrubs  ', 'open_shrubs    ', 'woody_savannas ', &
       'savannas       ', 'grasslands     ', 'perm_wetlands  ', &
       'crops          ', 'urban          ', 'crop_nature    ', &
       'snow_ice       ', 'barren         ', 'unclass        '/)
  TYPE(PTR_2D_ARRAY), DIMENSION(nmodisclass) :: modis_emis
  REAL(dp), DIMENSION(:,:,:), POINTER :: modis     => NULL()
  TYPE(T_CHAOBJ_CPL)                  :: imp_modis

!################################################
!#     MODIS LAI EMISSIONS (time-dependent)     #
!################################################
   REAL(dp), DIMENSION(:,:),  POINTER :: heald_emis => NULL()
   REAL(dp), DIMENSION(:,:),  POINTER :: js_emis => NULL()
   REAL(dp), DIMENSION(:,:),  POINTER :: hummel_emis => NULL()
   REAL(dp), DIMENSION(:,:),  POINTER :: modis_lai => NULL()
!################################################

  ! terrestrial isotopic carbon 13C signature
  REAL(dp), DIMENSION(:,:), POINTER :: terr13C_delta => NULL()
  TYPE(T_CHAOBJ_CPL)                :: imp_terr13C_delta

  ! isotopic ISOP emisflux
  REAL(dp), DIMENSION(:,:), POINTER :: I12ISOP_emflux => NULL()
  REAL(dp), DIMENSION(:,:), POINTER :: I13ISOP_emflux => NULL()

  ! AirSnow
  ! Import channels
  REAL(dp), DIMENSION(:,:), POINTER :: ddepflux_HOBr  => NULL()
  REAL(dp), DIMENSION(:,:), POINTER :: ddepflux_BrNO3 => NULL()
  REAL(dp), DIMENSION(:,:), POINTER :: ddepflux_HBr   => NULL()
  REAL(dp), DIMENSION(:,:), POINTER :: ddepflux_O3    => NULL()
  REAL(dp), DIMENSION(:,:), POINTER :: sic_multi_year => NULL()
  TYPE(T_CHAOBJ_CPL) :: imp_ddepflux_HOBr
  TYPE(T_CHAOBJ_CPL) :: imp_ddepflux_BrNO3
  TYPE(T_CHAOBJ_CPL) :: imp_ddepflux_HBr
  TYPE(T_CHAOBJ_CPL) :: imp_ddepflux_O3
  TYPE(T_CHAOBJ_CPL) :: imp_sic_multi_year
  ! Export channels bromine emission flux over snow/ice
  REAL(dp), DIMENSION(:,:), POINTER :: snow_air_flux_br2  => NULL()
  REAL(dp), DIMENSION(:,:), POINTER :: snow_air_flux_brcl => NULL()

  !#########################################################
  !### add channel object pointer for new emissions here ###
  !#########################################################

  ! GLOBAL PARAMETERS (CPL-NAMELIST)
  LOGICAL :: L_LG       = .false.  ! emissions for Lagrangian tracers

#ifdef MESSYTENDENCY
  INTEGER :: my_handle
#endif

  namelist /CPL/ L_LG, F2T

  namelist /CPL_IMPORT/ imp_seawater_dms, imp_CH4_conc_clim, imp_drymatter     &
       , imp_isop_emisfac, imp_mterp_emisfac, imp_lai, imp_lad_top             &
       , imp_cossza                                                            &
       , imp_lad_soil, imp_lad_lay2, imp_lad_lay3, imp_hc, imp_drag            &
       , imp_disp, imp_forestfr, imp_NOemisclass1, imp_NOemisclass2            &
       , imp_cultiv, imp_fertil, imp_fertil_yl95sl10, imp_rootdepth_yl95sl10   &
       , imp_NOemclass_yl95sl10, imp_lai_yl95sl10, imp_rootdepth_mask_yl95sl10 &
       , imp_du_cla, imp_du_src, imp_du_thr                                    &
       , imp_du_cla2, imp_rdepth, imp_dustsrc,imp_lai_in,imp_soiltexture       &
       , imp_du_nap, imp_du_kp, imp_du_capp, imp_du_mgpp, imp_du_misc &
       , imp_kkdu_clay, imp_kkdu_mask, imp_kkdu_lai, imp_kkdu_topo &
       , imp_kkdu_nap, imp_kkdu_kp, imp_kkdu_capp, imp_kkdu_mgpp, imp_kkdu_misc&
       , imp_mat_s2, imp_mat_s3, imp_mat_s4, imp_mat_s6, imp_mat_psrc          &
       , imp_k_fpar_eff, imp_numflx_as_aerocom, imp_numflx_cs_aerocom          &
       , imp_massflx_as_aerocom, imp_massflx_cs_aerocom, imp_OC_ag, imp_OC_ant &
       , imp_OC_bge, imp_OC_wf, imp_BC_ag, imp_BC_ant, imp_BC_wf               &
       , imp_SO2_ant_high, imp_SO2_ant_low                                     &
       , imp_olson, imp_modis                                                  &
       , imp_terr13C_delta                                                     &
       , imp_poc_aqua,  imp_poc_swifs, imp_wioc_aqua, imp_wioc_blen            &
       , imp_ss_poc_aqua,imp_ss_poc_swifs,imp_ss_wioc_aqua,imp_ss_wioc_blen    &
       , imp_ddepflux_HOBR, imp_ddepflux_BrNO3                                 &
       , imp_ddepflux_HBr, imp_ddepflux_O3                                     &
       , imp_sic_multi_year

  PUBLIC :: onemis_initialize
  PUBLIC :: onemis_init_memory
  PUBLIC :: onemis_init_coupling
  PUBLIC :: onemis_global_start
  PUBLIC :: onemis_vdiff
  PUBLIC :: onemis_global_end
  PUBLIC :: onemis_free_memory
  !PRIVATE :: onemis_read_nml_cpl
  !PRIVATE :: onemis_global_end_lg

CONTAINS

  ! --------------------------------------------------------------------------
  SUBROUTINE onemis_initialize

    ! ECHAM5/MESSy
    USE messy_main_mpi_bi,       ONLY: p_parallel_io, p_bcast, p_io
    USE messy_main_tools,        ONLY: find_next_free_unit
    USE messy_onemis,            ONLY: onemis_read_nml_ctrl, EMIS_TYPE

    IMPLICIT NONE

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'onemis_initialize'
    INTEGER                     :: iou    ! I/O unit
    INTEGER                     :: status ! error status
    INTEGER                     :: i
    INTEGER                     :: jf, jt ! loop indices

    CALL start_message_bi(modstr, 'GLOBAL setup',substr)

    ! read CTRL namelist
    IF (p_parallel_io) THEN
       iou = find_next_free_unit(100,200)
       CALL onemis_read_nml_ctrl(status, iou)
       IF (status /= 0) CALL error_bi( &
            'error in onemis_read_nml_ctrl', substr)
    END IF

    ! BROADCAST CTRL
    DO i=1, max_emis
       IF (p_parallel_io .AND. TRIM(EMIS_TYPE(i)) /= '') &
            write(*,*) 'EMIS_TYPE: ', EMIS_TYPE(i),' switched on!'
       CALL p_bcast(EMIS_TYPE(i), p_io)

       ! READ IN EXTRA CONTROL NAMELIST
       IF (TRIM(EMIS_TYPE(i)) == 'DU_tegen') THEN
          ! read CTRL_DU namelist
          IF (p_parallel_io) THEN
             iou = find_next_free_unit(100,200)
             CALL onemis_read_nml_ctrl_DU(status, iou)
             IF (status /= 0) CALL error_bi( &
                  'error in onemis_read_nml_ctrl_DU', substr)
          END IF
          ! Broadcast emission type specific switches
          CALL p_bcast(cuscale_in, p_io)
          CALL p_bcast(l_nudging,  p_io)
       END IF

       ! READ IN EXTRA CONTROL NAMELIST
       IF (TRIM(EMIS_TYPE(i)) == 'NO_yl95sl10') THEN
          ! read CTRL NOsl10 namelist
          IF (p_parallel_io) THEN
             iou = find_next_free_unit(100,200)
             CALL onemis_read_nml_ctrl_NOsl10(status, iou)
             IF (status /= 0) CALL error_bi( &
                  'error in onemis_read_nml_ctrl_NOsl10', substr)
          END IF
          ! Broadcast emission type specific switches
          CALL p_bcast(noemfact_wet_yl95sl10, p_io)
          CALL p_bcast(noemfact_dry_yl95sl10, p_io)
          CALL p_bcast(smoist_method, p_io)
       END IF

       ! READ IN EXTRA CONTROL NAMELIST
       IF (TRIM(EMIS_TYPE(i)) == 'AirSnow' ) THEN
          ! read CTRL AirSnow namelist
          IF (p_parallel_io) THEN
             iou = find_next_free_unit(100,200)
             CALL onemis_read_nml_ctrl_AirSnow(status, iou)
             IF (status /= 0) CALL error_bi ( &
                  'error in onemis_read_nml_ctrl_AirSnow', substr)
          END IF
       END IF
    END DO

   ! read CPL namelist
    IF (p_parallel_io) THEN
       iou = find_next_free_unit(100,200)
       CALL onemis_read_nml_cpl(status, iou)
       IF (status /= 0) CALL error_bi( &
            'error in onemis_read_nml_cpl', substr)
    END IF

    ! read CPL namelist
    IF (p_parallel_io) THEN
       iou = find_next_free_unit(100,200)
       CALL onemis_read_nml_cpl_import(status, iou)
       IF (status /= 0) CALL error_bi( &
            'error in onemis_read_nml_cpl_import' &
            , substr)
    END IF

    call p_bcast(imp_seawater_dms%CHA , p_io)
    call p_bcast(imp_seawater_dms%OBJ , p_io)
    call p_bcast(imp_CH4_conc_clim%CHA, p_io)
    call p_bcast(imp_CH4_conc_clim%OBJ, p_io)
    call p_bcast(imp_cossza%CHA       , p_io)
    call p_bcast(imp_cossza%OBJ       , p_io)
    call p_bcast(imp_drymatter%CHA    , p_io)
    call p_bcast(imp_drymatter%OBJ    , p_io)
    call p_bcast(imp_isop_emisfac%CHA , p_io)
    call p_bcast(imp_isop_emisfac%OBJ , p_io)
    call p_bcast(imp_mterp_emisfac%CHA, p_io)
    call p_bcast(imp_mterp_emisfac%OBJ, p_io)
    call p_bcast(imp_lai%CHA          , p_io)
    call p_bcast(imp_lai%OBJ          , p_io)
    call p_bcast(imp_lad_top%CHA      , p_io)
    call p_bcast(imp_lad_top%OBJ      , p_io)
    call p_bcast(imp_lad_soil%CHA     , p_io)
    call p_bcast(imp_lad_soil%OBJ     , p_io)
    call p_bcast(imp_lad_lay2%CHA     , p_io)
    call p_bcast(imp_lad_lay2%OBJ     , p_io)
    call p_bcast(imp_lad_lay3%CHA     , p_io)
    call p_bcast(imp_lad_lay3%OBJ     , p_io)
    call p_bcast(imp_hc%CHA           , p_io)
    call p_bcast(imp_hc%OBJ           , p_io)
    call p_bcast(imp_drag%CHA         , p_io)
    call p_bcast(imp_drag%OBJ         , p_io)
    call p_bcast(imp_disp%CHA         , p_io)
    call p_bcast(imp_disp%OBJ         , p_io)
    call p_bcast(imp_forestfr%CHA     , p_io)
    call p_bcast(imp_forestfr%OBJ     , p_io)
    call p_bcast(imp_NOemisclass1%CHA , p_io)
    call p_bcast(imp_NOemisclass1%OBJ , p_io)
    call p_bcast(imp_NOemisclass2%CHA , p_io)
    call p_bcast(imp_NOemisclass2%OBJ , p_io)
    call p_bcast(imp_cultiv%CHA       , p_io)
    call p_bcast(imp_cultiv%OBJ       , p_io)
    call p_bcast(imp_fertil%CHA       , p_io)
    call p_bcast(imp_fertil%OBJ       , p_io)
    call p_bcast(imp_lai_yl95sl10%CHA,           p_io)
    call p_bcast(imp_lai_yl95sl10%OBJ,           p_io)
    call p_bcast(imp_fertil_yl95sl10%CHA,        p_io)
    call p_bcast(imp_fertil_yl95sl10%OBJ,        p_io)
    call p_bcast(imp_NOemclass_yl95sl10%CHA,     p_io)
    call p_bcast(imp_NOemclass_yl95sl10%OBJ,     p_io)
    call p_bcast(imp_rootdepth_yl95sl10%CHA,     p_io)
    call p_bcast(imp_rootdepth_yl95sl10%OBJ,     p_io)
    call p_bcast(imp_rootdepth_mask_yl95sl10%CHA,p_io)
    call p_bcast(imp_rootdepth_mask_yl95sl10%OBJ,p_io)
    call p_bcast(imp_du_cla%CHA,                 p_io)
    call p_bcast(imp_du_cla%OBJ,                 p_io)
    call p_bcast(imp_du_src%CHA,                 p_io)
    call p_bcast(imp_du_src%OBJ,                 p_io)
    call p_bcast(imp_du_thr%CHA,                 p_io)
    call p_bcast(imp_du_thr%OBJ,                 p_io)

    call p_bcast(imp_du_cla2%CHA,                p_io)
    call p_bcast(imp_du_cla2%OBJ,                p_io)
    call p_bcast(imp_rdepth%CHA,                 p_io)
    call p_bcast(imp_rdepth%OBJ,                 p_io)
    call p_bcast(imp_dustsrc%CHA,                p_io)
    call p_bcast(imp_dustsrc%OBJ,                p_io)
    call p_bcast(imp_lai_in%CHA,                 p_io)
    call p_bcast(imp_lai_in%OBJ,                 p_io)
    call p_bcast(imp_soiltexture%CHA,            p_io)
    call p_bcast(imp_soiltexture%OBJ,            p_io)

    call p_bcast(imp_kkdu_clay%CHA, p_io)
    call p_bcast(imp_kkdu_clay%OBJ, p_io)
    call p_bcast(imp_kkdu_mask%CHA, p_io)
    call p_bcast(imp_kkdu_mask%OBJ, p_io)
    call p_bcast(imp_kkdu_lai%CHA, p_io)
    call p_bcast(imp_kkdu_lai%OBJ, p_io)
    call p_bcast(imp_kkdu_topo%CHA, p_io)
    call p_bcast(imp_kkdu_topo%OBJ, p_io)
    call p_bcast(imp_kkdu_nap%CHA, p_io)
    call p_bcast(imp_kkdu_nap%OBJ, p_io)
    call p_bcast(imp_kkdu_kp%CHA, p_io)
    call p_bcast(imp_kkdu_kp%OBJ, p_io)
    call p_bcast(imp_kkdu_capp%CHA, p_io)
    call p_bcast(imp_kkdu_capp%OBJ, p_io)
    call p_bcast(imp_kkdu_mgpp%CHA, p_io)
    call p_bcast(imp_kkdu_mgpp%OBJ, p_io)
    call p_bcast(imp_kkdu_misc%CHA, p_io)
    call p_bcast(imp_kkdu_misc%OBJ, p_io)

    IF (imp_du_nap%CHA  == '') l_ducomp=.FALSE.
    IF (imp_du_nap%OBJ  == '') l_ducomp=.FALSE.
    IF (imp_du_kp%CHA   == '') l_ducomp=.FALSE.
    IF (imp_du_kp%OBJ   == '') l_ducomp=.FALSE.
    IF (imp_du_capp%CHA == '') l_ducomp=.FALSE.
    IF (imp_du_capp%OBJ == '') l_ducomp=.FALSE.
    IF (imp_du_mgpp%CHA == '') l_ducomp=.FALSE.
    IF (imp_du_mgpp%OBJ == '') l_ducomp=.FALSE.
    IF (imp_du_misc%CHA == '') l_ducomp=.FALSE.
    IF (imp_du_misc%OBJ == '') l_ducomp=.FALSE.

    call p_bcast(l_ducomp,                      p_io)

    IF (l_ducomp) THEN
       call p_bcast(imp_du_nap%CHA,                p_io)
       call p_bcast(imp_du_nap%OBJ,                p_io)
       call p_bcast(imp_du_kp%CHA,                 p_io)
       call p_bcast(imp_du_kp%OBJ,                 p_io)
       call p_bcast(imp_du_capp%CHA,               p_io)
       call p_bcast(imp_du_capp%OBJ,               p_io)
       call p_bcast(imp_du_mgpp%CHA,               p_io)
       call p_bcast(imp_du_mgpp%OBJ,               p_io)
       call p_bcast(imp_du_misc%CHA,               p_io)
       call p_bcast(imp_du_misc%OBJ,               p_io)
    ENDIF

    call p_bcast(imp_mat_s2%CHA,                 p_io)
    call p_bcast(imp_mat_s2%OBJ,                 p_io)
    call p_bcast(imp_mat_s3%CHA,                 p_io)
    call p_bcast(imp_mat_s3%OBJ,                 p_io)
    call p_bcast(imp_mat_s4%CHA,                 p_io)
    call p_bcast(imp_mat_s4%OBJ,                 p_io)
    call p_bcast(imp_mat_s6%CHA,                 p_io)
    call p_bcast(imp_mat_s6%OBJ,                 p_io)
    call p_bcast(imp_mat_psrc%CHA,               p_io)
    call p_bcast(imp_mat_psrc%OBJ,               p_io)
    call p_bcast(imp_k_fpar_eff%CHA,             p_io)
    call p_bcast(imp_k_fpar_eff%OBJ,             p_io)
    call p_bcast(imp_numflx_as_aerocom%CHA,      p_io)
    call p_bcast(imp_numflx_as_aerocom%OBJ,      p_io)
    call p_bcast(imp_numflx_cs_aerocom%CHA,      p_io)
    call p_bcast(imp_numflx_cs_aerocom%OBJ,      p_io)
    call p_bcast(imp_massflx_as_aerocom%CHA,     p_io)
    call p_bcast(imp_massflx_as_aerocom%OBJ,     p_io)
    call p_bcast(imp_massflx_cs_aerocom%CHA,     p_io)
    call p_bcast(imp_massflx_cs_aerocom%OBJ,     p_io)
    call p_bcast(imp_OC_ag%CHA,         p_io)
    call p_bcast(imp_OC_ag%OBJ,         p_io)
    call p_bcast(imp_OC_ant%CHA,        p_io)
    call p_bcast(imp_OC_ant%OBJ,        p_io)
    call p_bcast(imp_OC_bge%CHA,        p_io)
    call p_bcast(imp_OC_bge%OBJ,        p_io)
    call p_bcast(imp_OC_wf%CHA,         p_io)
    call p_bcast(imp_OC_wf%OBJ,         p_io)
    call p_bcast(imp_BC_ag%CHA,         p_io)
    call p_bcast(imp_BC_ag%OBJ,         p_io)
    call p_bcast(imp_BC_ant%CHA,        p_io)
    call p_bcast(imp_BC_ant%OBJ,        p_io)
    call p_bcast(imp_BC_wf%CHA,         p_io)
    call p_bcast(imp_BC_wf%OBJ,         p_io)
    call p_bcast(imp_SO2_ant_high%CHA,  p_io)
    call p_bcast(imp_SO2_ant_high%OBJ,  p_io)
    call p_bcast(imp_SO2_ant_low%CHA,   p_io)
    call p_bcast(imp_SO2_ant_low%OBJ,   p_io)
    call p_bcast(imp_olson%CHA,         p_io)
    call p_bcast(imp_olson%OBJ,         p_io)
    call p_bcast(imp_modis%CHA,         p_io)
    call p_bcast(imp_modis%OBJ,         p_io)
    call p_bcast(imp_terr13C_delta%CHA, p_io)
    call p_bcast(imp_terr13C_delta%OBJ, p_io)

    call p_bcast(imp_poc_aqua%CHA,      p_io)
    call p_bcast(imp_poc_aqua%OBJ,      p_io)
    call p_bcast(imp_poc_swifs%CHA,     p_io)
    call p_bcast(imp_poc_swifs%OBJ,     p_io)
    call p_bcast(imp_wioc_aqua%CHA,     p_io)
    call p_bcast(imp_wioc_aqua%OBJ,     p_io)
    call p_bcast(imp_wioc_blen%CHA,     p_io)
    call p_bcast(imp_wioc_blen%OBJ,     p_io)
    call p_bcast(imp_ss_poc_aqua%CHA,   p_io)
    call p_bcast(imp_ss_poc_aqua%OBJ,   p_io)
    call p_bcast(imp_ss_poc_swifs%CHA,  p_io)
    call p_bcast(imp_ss_poc_swifs%OBJ,  p_io)
    call p_bcast(imp_ss_wioc_aqua%CHA,  p_io)
    call p_bcast(imp_ss_wioc_aqua%OBJ,  p_io)
    call p_bcast(imp_ss_wioc_blen%CHA,  p_io)
    call p_bcast(imp_ss_wioc_blen%OBJ,  p_io)

    call p_bcast(imp_ddepflux_HOBr%CHA, p_io)
    call p_bcast(imp_ddepflux_HOBr%OBJ, p_io)
    call p_bcast(imp_ddepflux_BrNO3%CHA,p_io)
    call p_bcast(imp_ddepflux_BrNO3%OBJ,p_io)
    call p_bcast(imp_ddepflux_HBr%CHA,  p_io)
    call p_bcast(imp_ddepflux_HBr%OBJ,  p_io)
    call p_bcast(imp_ddepflux_O3%CHA,   p_io)
    call p_bcast(imp_ddepflux_O3%OBJ,   p_io)
    call p_bcast(imp_sic_multi_year%CHA,p_io)
    call p_bcast(imp_sic_multi_year%OBJ,p_io)

   ! BROADCAST RESULTS
    CALL p_bcast(L_LG, p_io)

    ! ANALYZE F2T
    IF (p_parallel_io) THEN
       ! GET NUMBER OF FLUXES
       NFLUXES = 1
       DO jf=1, NMAXFLUXES
          ! CHECK name
          IF (TRIM(F2T(jf)%name) == '') CYCLE
          XF2T(NFLUXES)%name = TRIM(ADJUSTL(F2T(jf)%name))
          WRITE(*,*) ' FLUX<->TRACER SET : ',TRIM(XF2T(NFLUXES)%name)
          !
          ! PARSE STRING
          CALL parse_f2tstr(status, MAXPSTRLEN, F2T(jf)%string_gp &
               , XF2T(NFLUXES)%tracer_gp                          &
               , XF2T(NFLUXES)%mgp, XF2T(NFLUXES)%efact_gp        &
               )
          CALL check_status(status, substr, '(F2T(GP))')
          XF2T(NFLUXES)%ngpt = SIZE(XF2T(NFLUXES)%tracer_gp)
          ALLOCATE(XF2T(NFLUXES)%idt_gp(XF2T(NFLUXES)%ngpt))
          XF2T(NFLUXES)%idt_gp(:) = 0
          !
          WRITE(*,*) '   NO. OF GRIDPOINT TRACER(S): ',XF2T(NFLUXES)%ngpt
          DO jt=1, XF2T(NFLUXES)%ngpt
             WRITE(*,*) '    ',TRIM(XF2T(NFLUXES)%tracer_gp(jt)),' (',&
                  'METHOD=',XF2T(NFLUXES)%mgp(jt), &
                  ', SC=',XF2T(NFLUXES)%efact_gp(jt),')'
#ifdef COSMO
             IF (XF2T(NFLUXES)%mgp(jt) ==2 ) THEN
                write(*,*) '**** METHOD = 2 not possible for COSMO ****'
                write(*,*) '****        GP METHOD set to 1         ****'
                XF2T(NFLUXES)%mgp(jt)=1
             ENDIF
#endif
          END DO

!!#D attila +
#ifdef ECHAM5
          IF (L_LG) THEN
             CALL parse_f2tstr(status, MAXPSTRLEN, F2T(jf)%string_lg &
                  , XF2T(NFLUXES)%tracer_lg &
                  , XF2T(NFLUXES)%mlg, XF2T(NFLUXES)%efact_lg &
                  )
             CALL check_status(status, substr, '(F2T(LG))')
             XF2T(NFLUXES)%nlgt = SIZE(XF2T(NFLUXES)%tracer_lg)
             ALLOCATE(XF2T(NFLUXES)%idt_lg(XF2T(NFLUXES)%nlgt))
             XF2T(NFLUXES)%idt_lg(:) = 0
             !
             WRITE(*,*) '   NO. OF LAGRANGIAN TRACER(S): ',XF2T(NFLUXES)%nlgt
             DO jt=1, XF2T(NFLUXES)%nlgt
                WRITE(*,*) '    ',TRIM(XF2T(NFLUXES)%tracer_lg(jt)),' (',&
                     'METHOD=',XF2T(NFLUXES)%mlg(jt), &
                     ', SC=',XF2T(NFLUXES)%efact_lg(jt),')'
             END DO
          END IF
#endif
!!#D attila -
          !
          ! NEXT SET
          NFLUXES = NFLUXES + 1
          WRITE(*,*) '------------------------------------------------------'
       END DO
       NFLUXES = NFLUXES - 1
    END IF
    !
    ! BROADCAST RESULTS
    CALL p_bcast(NFLUXES, p_io)
    !
    DO jf=1, NFLUXES
       ! GRIDPOINT
       CALL p_bcast(XF2T(jf)%name, p_io)
       !
       CALL p_bcast(XF2T(jf)%ngpt, p_io)
       IF (.NOT. p_parallel_io) THEN
          ALLOCATE(XF2T(jf)%mgp(XF2T(jf)%ngpt))
          ALLOCATE(XF2T(jf)%efact_gp(XF2T(jf)%ngpt))
          ALLOCATE(XF2T(jf)%tracer_gp(XF2T(jf)%ngpt))
          ALLOCATE(XF2T(jf)%idt_gp(XF2T(jf)%ngpt))
       END IF
       CALL p_bcast(XF2T(jf)%mgp, p_io)
       CALL p_bcast(XF2T(jf)%efact_gp, p_io)
       CALL p_bcast(XF2T(jf)%idt_gp, p_io)
       DO jt=1, XF2T(jf)%ngpt
          CALL p_bcast(XF2T(jf)%tracer_gp(jt), p_io)
       END DO

!!#D attila +
#ifdef ECHAM5
       ! LAGRANGE
       IF (L_LG) THEN
          CALL p_bcast(XF2T(jf)%nlgt, p_io)
          !
          CALL p_bcast(XF2T(jf)%nlgt, p_io)
          IF (.NOT. p_parallel_io) THEN
             ALLOCATE(XF2T(jf)%mlg(XF2T(jf)%nlgt))
             ALLOCATE(XF2T(jf)%efact_lg(XF2T(jf)%nlgt))
             ALLOCATE(XF2T(jf)%tracer_lg(XF2T(jf)%nlgt))
             ALLOCATE(XF2T(jf)%idt_lg(XF2T(jf)%nlgt))
          END IF
          CALL p_bcast(XF2T(jf)%mlg, p_io)
          CALL p_bcast(XF2T(jf)%efact_lg, p_io)
          CALL p_bcast(XF2T(jf)%idt_lg, p_io)
          DO jt=1, XF2T(jf)%nlgt
             CALL p_bcast(XF2T(jf)%tracer_lg(jt), p_io)
          END DO
       END IF
#endif
!!#D attila -

    END DO

    CALL end_message_bi(modstr, 'GLOBAL setup',substr)

CONTAINS

  SUBROUTINE check_status(status, substr, str)

    INTEGER,          INTENT(IN) :: status
    CHARACTER(LEN=*), INTENT(IN) :: substr
    CHARACTER(LEN=*), INTENT(IN) :: str

    SELECT CASE (status)
    CASE(0)
       RETURN
    CASE(1)
       ! F2T
       CALL error_bi( 'EMPTY TRACER-NAME IN TRACER SPECIFICATION '&
            &//str,substr)
    CASE(2)
       CALL error_bi( 'SYNTAX ERROR IN TRACER SPECIFICATION '&
            &//str,substr)
    CASE(3)
       CALL error_bi( 'TOO MANY SWITCH SPECIFICATIONS '&
            &//str,substr)
    CASE(4)
       CALL error_bi( 'MISSING VALUE IN SWITCH SPECIFICATION '&
            &//str,substr)
    CASE(5)
       CALL error_bi( 'SYNTAX ERROR IN SWITCH SPECIFICATION '&
            &//str,substr)
    CASE(6)
       CALL error_bi( 'READ ERROR FOR EMISSION METHOD '&
            &//str,substr)
    CASE(7)
       CALL error_bi( 'READ ERROR FOR SCALING FACTOR '&
            &//str,substr)
    CASE(8)
       CALL error_bi( 'UNKNOWN SWITCH SPECIFICATION '&
            &//str,substr)
       ! RGT-EVENT action
    CASE(50)
       CALL error_bi( 'EMPTY RGT-EVENT ACTION'&
            &//str,substr)
       !
    CASE DEFAULT
       CALL error_bi( 'SEVERE ERROR !',substr)
       !
    END SELECT

  END SUBROUTINE check_status

END SUBROUTINE onemis_initialize
! --------------------------------------------------------------------------

  ! --------------------------------------------------------------------------
  SUBROUTINE onemis_init_memory

    USE messy_main_timer,            ONLY: lstart
    USE messy_main_grid_def_mem_bi,  ONLY: nproma, ngpblks
    USE messy_main_channel_error_bi, ONLY: channel_halt
    USE messy_main_channel_bi,       ONLY: DC_GP             &
                                         , DIMID_LON, DIMID_LAT            &
                                         , GP_3D_MID, GP_2D_HORIZONTAL     &
                                         , GP_3D_1LEV                      &
                                         , gp_nseg, gp_start, gp_cnt       &
                                         , gp_meml, gp_memu
    USE messy_main_channel_dimensions,  ONLY: new_dimension, get_dimension_info
    USE messy_main_channel_repr,        ONLY: new_representation, AUTO  &
                                            , set_representation_decomp &
                                            , get_representation_id     &
                                            , IRANK, PIOTYPE_COL        &
                                            , repr_def_axes
    USE messy_main_channel,             ONLY: new_channel, new_channel_object &
                                            , new_attribute
    USE messy_main_tools,               ONLY: strcrack ! mz_ht_20110513
    IMPLICIT NONE

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'onemis_init_memory'
    INTEGER                     :: status
    INTEGER                     :: je
    INTEGER                     :: DIMID_ncl_yl95
    INTEGER                     :: DIMID_ncl_yl95sl10
    INTEGER                     :: DIMID_numveglay_hr
    INTEGER                     :: DIMID_ndrydays
    INTEGER                     :: REPRID_ncl_yl95
    INTEGER                     :: REPRID_ncl_yl95sl10
    INTEGER                     :: REPRID_numveglay_hr
    INTEGER                     :: REPRID_ndrydays
    INTEGER                     :: i

    INTEGER                     :: DIMID_du3_maxsizes
    INTEGER                     :: DIMID_du3_npt
    INTEGER                     :: REPRID_du3_maxsizes
    INTEGER                     :: DU3_npt

    ! PARALLEL DECOMPOSITION
    INTEGER                          :: nseg = 0
    INTEGER, DIMENSION(:,:), POINTER :: start => NULL()
    INTEGER, DIMENSION(:,:), POINTER :: cnt   => NULL()
    INTEGER, DIMENSION(:,:), POINTER :: meml  => NULL()
    INTEGER, DIMENSION(:,:), POINTER :: memu  => NULL()
    INTEGER                          :: rep_id

    CHARACTER(LEN=STRLEN_CHANNEL), POINTER  :: strname_cha(:) => NULL()
    CHARACTER(LEN=STRLEN_OBJECT),  POINTER  :: strname_obj(:) => NULL()
    CHARACTER(LEN=STRLEN_OBJECT)            :: name
    INTEGER                                 :: ni

    ! moved  to init_memory
#ifdef MESSYTENDENCY
    my_handle = mtend_get_handle(modstr)
    CALL mtend_register(my_handle, mtend_id_tracer)
#endif

    CALL start_message_bi(modstr, 'CHANNEL DEFINITION', substr)

    ! define new channel
    CALL new_channel(status, modstr, reprid=GP_3D_1LEV)
    CALL channel_halt(substr, status)

    emistype_loop1: DO je = 1, max_emis
       IF (TRIM(EMIS_TYPE(je)) == '' ) CYCLE

       SELECT CASE(TRIM(EMIS_TYPE(je)))

       CASE('DMS')
          CALL new_channel_object(status, modstr, &
               'emis_dms_sea', p2=emis_dms_sea )
          CALL channel_halt(substr, status)
          CALL new_attribute(status, modstr, &
               'emis_dms_sea', 'long_name', c='DMS emission' )
          CALL channel_halt(substr, status)
          CALL new_attribute(status, modstr, &
               'emis_dms_sea', 'units', c='molec m-2 s-1' )
          CALL channel_halt(substr, status)

       CASE('O3ice')
          CALL new_channel_object(status, modstr, 'O3_emflux', p2=O3_emflux )
          CALL channel_halt(substr, status)
          CALL new_attribute(status, modstr, &
               'O3_emflux', 'long_name', c='O3 over ice emission' )
          CALL channel_halt(substr, status)
          CALL new_attribute(status, modstr, &
               'O3_emflux', 'units', c='molec m-2 s-1' )
          CALL channel_halt(substr, status)

       CASE ('SS_lsce')
          CALL new_channel_object(status, modstr, &
               'mss_as_lsce', p2=mss_as_lsce )
          CALL channel_halt(substr, status)
          CALL new_attribute(status, modstr, &
               'mss_as_lsce', 'long_name',   &
               c='accumulation sea-salt mass flux (lsce)' )
          CALL channel_halt(substr, status)
          CALL new_attribute(status, modstr, &
               'mss_as_lsce', 'units', c='kg m-2 s-1' )
          CALL channel_halt(substr, status)

          CALL new_channel_object(status, modstr, 'mss_cs_lsce', p2=mss_cs_lsce )
          CALL channel_halt(substr, status)
          CALL new_attribute(status, modstr, &
               'mss_cs_lsce', 'long_name',   &
               c='coarse sea-salt mass flux (lsce)' )
          CALL channel_halt(substr, status)
          CALL new_attribute(status, modstr, &
               'mss_cs_lsce', 'units', c='kg m-2 s-1' )
          CALL channel_halt(substr, status)

          CALL new_channel_object(status, modstr, &
               'nss_as_lsce', p2=nss_as_lsce )
          CALL channel_halt(substr, status)
          CALL new_attribute(status, modstr, &
               'nss_as_lsce', 'long_name',   &
               c='accumulation sea-salt number flux (lsce)' )
          CALL channel_halt(substr, status)
          CALL new_attribute(status, modstr, &
               'nss_as_lsce', 'units', c='m-2 s-1' )
          CALL channel_halt(substr, status)

          CALL new_channel_object(status, modstr, &
               'nss_cs_lsce', p2=nss_cs_lsce )
          CALL channel_halt(substr, status)
          CALL new_attribute(status, modstr, &
               'nss_cs_lsce', 'long_name',   &
               c='coarse sea-salt number flux (lsce)' )
          CALL channel_halt(substr, status)
          CALL new_attribute(status, modstr, &
               'nss_cs_lsce', 'units', c='m-2 s-1' )
          CALL channel_halt(substr, status)

          CALL new_channel_object(status, modstr, &
               'emisflx_ss_as_sum', p2=emisflx_ss_as_sum, lrestreq=.TRUE. )
          CALL channel_halt(substr, status)
          CALL new_attribute(status, modstr,     &
               'emisflx_ss_as_sum', 'long_name', &
               c='sum of accumulation mode seasalt emissions' )
          CALL channel_halt(substr, status)
          CALL new_attribute(status, modstr, &
               'emisflx_ss_as_sum', 'units', &
               c='mol(NaCl) / mol(air) * kg(air) m-2 s-1' )
          CALL channel_halt(substr, status)

          CALL new_channel_object(status, modstr, &
               'emisflx_ss_cs_sum', p2=emisflx_ss_cs_sum, lrestreq=.TRUE. )
          CALL channel_halt(substr, status)
          CALL new_attribute(status, modstr, &
               'emisflx_ss_cs_sum', 'long_name', &
               c='sum of coarse mode seasalt emissions' )
          CALL channel_halt(substr, status)
          CALL new_attribute(status, modstr, &
               'emisflx_ss_cs_sum', 'units', &
               c='mol(NaCl) / mol(air) * kg(air) m-2 s-1' )
          CALL channel_halt(substr, status)

       CASE('SS_monahan')
          CALL new_channel_object(status, modstr, &
               'mss_as_monahan', p2=mss_as_monahan )
          CALL channel_halt(substr, status)
          CALL new_attribute(status, modstr,  &
               'mss_as_monahan', 'long_name', &
               c='accumulation sea-salt mass flux (Monaham)' )
          CALL channel_halt(substr, status)
          CALL new_attribute(status, modstr, &
               'mss_as_monahan', 'units', c='kg m-2 s-1' )
          CALL channel_halt(substr, status)

          CALL new_channel_object(status, modstr, &
               'mss_cs_monahan', p2= mss_cs_monahan)
          CALL channel_halt(substr, status)
          CALL new_attribute(status, modstr,  &
               'mss_cs_monahan', 'long_name', &
               c='coarse sea-salt mass flux (Monaham)' )
          CALL channel_halt(substr, status)
          CALL new_attribute(status, modstr, &
               'mss_cs_monahan', 'units', c='kg m-2 s-1' )
          CALL channel_halt(substr, status)

          CALL new_channel_object(status, modstr, &
               'nss_as_monahan', p2=nss_as_monahan )
          CALL channel_halt(substr, status)
          CALL new_attribute(status, modstr,  &
               'nss_as_monahan', 'long_name', &
               c='accumulation sea-salt number flux (Monaham)' )
          CALL channel_halt(substr, status)
          CALL new_attribute(status, modstr, &
               'nss_as_monahan', 'units', c='m-2 s-1' )
          CALL channel_halt(substr, status)

          CALL new_channel_object(status, modstr, &
               'nss_cs_monahan', p2=nss_cs_monahan )
          CALL channel_halt(substr, status)
          CALL new_attribute(status, modstr,  &
               'nss_cs_monahan', 'long_name', &
               c='coarse sea-salt number flux (Monaham)' )
          CALL channel_halt(substr, status)
          CALL new_attribute(status, modstr, &
               'nss_cs_monahan', 'units', c='m-2 s-1' )
          CALL channel_halt(substr, status)

       CASE ('SS_aerocom')
          ! define output channel objects
          CALL new_channel_object(status, modstr, &
               'mss_as_aerocom', p2=mss_as_aerocom )
          CALL channel_halt(substr, status)
          CALL new_attribute(status, modstr,  &
               'mss_as_aerocom', 'long_name', &
               c='accumulation sea-salt mass flux (aerocom)' )
          CALL channel_halt(substr, status)
          CALL new_attribute(status, modstr, &
               'mss_as_aerocom', 'units', c='kg m-2 s-1' )
          CALL channel_halt(substr, status)

          CALL new_channel_object(status, modstr, &
               'mss_cs_aerocom', p2=mss_cs_aerocom )
          CALL channel_halt(substr, status)
          CALL new_attribute(status, modstr,  &
               'mss_cs_aerocom', 'long_name', &
               c='coarse sea-salt mass flux (aerosom)' )
          CALL channel_halt(substr, status)
          CALL new_attribute(status, modstr, &
               'mss_cs_aerocom', 'units', c='kg m-2 s-1' )
          CALL channel_halt(substr, status)

          CALL new_channel_object(status, modstr, &
               'nss_as_aerocom', p2=nss_as_aerocom )
          CALL channel_halt(substr, status)
          CALL new_attribute(status, modstr,  &
               'nss_as_aerocom', 'long_name', &
               c='accumulation sea-salt number flux (aerocom)' )
          CALL channel_halt(substr, status)
          CALL new_attribute(status, modstr, &
               'nss_as_aerocom', 'units', c='m-2 s-1' )
          CALL channel_halt(substr, status)

          CALL new_channel_object(status, modstr, &
               'nss_cs_aerocom', p2=nss_cs_aerocom )
          CALL channel_halt(substr, status)
          CALL new_attribute(status, modstr,  &
               'nss_cs_aerocom', 'long_name', &
               c='coarse sea-salt number flux (aerocom)' )
          CALL channel_halt(substr, status)
          CALL new_attribute(status, modstr, &
               'nss_cs_aerocom', 'units', c='m-2 s-1' )
          CALL channel_halt(substr, status)

          ! define diagnostic channel objects
          IF(.not. associated(emisflx_ss_as_sum)) THEN
             CALL new_channel_object(status, modstr, &
                  'emisflx_ss_as_sum', p2=emisflx_ss_as_sum, lrestreq=.TRUE. )
             CALL channel_halt(substr, status)
             CALL new_attribute(status, modstr,  &
                  'emisflx_ss_as_sum', 'long_name', &
                  c='sum of accumulation mode seasalt emissions' )
             CALL channel_halt(substr, status)
             CALL new_attribute(status, modstr, &
                  'emisflx_ss_as_sum', 'units', &
                  c='mol(NaCl) / mol(air) * kg(air) m-2 s-1' )
             CALL channel_halt(substr, status)
          END IF

          IF(.not. associated(emisflx_ss_cs_sum)) THEN
             CALL new_channel_object(status, modstr, &
                  'emisflx_ss_cs_sum', p2=emisflx_ss_cs_sum, lrestreq=.TRUE. )
             CALL channel_halt(substr, status)
             CALL new_attribute(status, modstr,  &
                  'emisflx_ss_cs_sum', 'long_name', &
                  c='sum of coarse mode seasalt emissions' )
             CALL channel_halt(substr, status)
             CALL new_attribute(status, modstr, &
                  'emisflx_ss_cs_sum', 'units', &
                  c='mol(NaCl) / mol(air) * kg(air) m-2 s-1' )
             CALL channel_halt(substr, status)
          END IF

       CASE('SS_POC_AQUA')
         if (associated(strname_cha)) DEALLOCATE (strname_cha)
         NULLIFY(strname_cha)
         if (associated(strname_obj)) DEALLOCATE (strname_obj)
         NULLIFY(strname_obj)
         CALL strcrack(imp_ss_poc_aqua%cha,';', strname_cha, n_ss_poc_aqua)
         CALL strcrack(imp_ss_poc_aqua%obj,';', strname_obj, n_ss_poc_aqua)

         ALLOCATE(emis_poc_aqua(n_ss_poc_aqua))

         DO ni=1,n_ss_poc_aqua
           emis_poc_aqua(ni)%channel= strname_cha(ni)
           emis_poc_aqua(ni)%object = strname_obj(ni)

           name='emis_poc_aqua_'//TRIM(emis_poc_aqua(ni)%object)
           CALL new_channel_object(status, modstr, &
                TRIM(name), p2=emis_poc_aqua(ni)%emis_flx, lrestreq=.TRUE. )
           CALL channel_halt(substr, status)
           CALL new_attribute(status, modstr,  &
                TRIM(name), 'long_name', &
                c='Ocean POC emissions (based on AQUA)' )
           CALL channel_halt(substr, status)
           CALL new_attribute(status, modstr, &
                TRIM(name), 'units', &
                c='kg (POC)  m-2 s-1' )
           CALL channel_halt(substr, status)

           name='emissum_poc_aqua_'//TRIM(emis_poc_aqua(ni)%object)
           CALL new_channel_object(status, modstr, &
                TRIM(name), p2=emis_poc_aqua(ni)%emis_flxsum, lrestreq=.TRUE. )
           CALL channel_halt(substr, status)
           CALL new_attribute(status, modstr,  &
                TRIM(name), 'long_name', &
                c='accumulated Ocean POC emissions (based on AQUA)' )
           CALL channel_halt(substr, status)
           CALL new_attribute(status, modstr, &
                TRIM(name), 'units', &
                c='kg (POC)  m-2' )
           CALL channel_halt(substr, status)
         ENDDO

      CASE('SS_POC_SWIFS')
         if (associated(strname_cha)) DEALLOCATE (strname_cha)
         NULLIFY(strname_cha)
         if (associated(strname_obj)) DEALLOCATE (strname_obj)
         NULLIFY(strname_obj)
         CALL strcrack(imp_ss_poc_swifs%cha,';', strname_cha, n_ss_poc_swifs)
         CALL strcrack(imp_ss_poc_swifs%obj,';', strname_obj, n_ss_poc_swifs)

         ALLOCATE(emis_poc_swifs(n_ss_poc_swifs))

         DO ni=1,n_ss_poc_swifs
           emis_poc_swifs(ni)%channel= strname_cha(ni)
           emis_poc_swifs(ni)%object = strname_obj(ni)

           name='emis_poc_swifs_'//TRIM(emis_poc_swifs(ni)%object)
           CALL new_channel_object(status, modstr, &
                TRIM(name), p2=emis_poc_swifs(ni)%emis_flx, lrestreq=.TRUE. )
           CALL channel_halt(substr, status)
           CALL new_attribute(status, modstr,  &
                TRIM(name), 'long_name', &
                c='Ocean POC emissions (based on SEAWIFS)' )
           CALL channel_halt(substr, status)
           CALL new_attribute(status, modstr, &
                TRIM(name), 'units', &
                  c='kg (POC)  m-2 s-1' )
           CALL channel_halt(substr, status)

           name='emissum_poc_swifs_'//TRIM(emis_poc_swifs(ni)%object)
           CALL new_channel_object(status, modstr, &
                TRIM(name), p2=emis_poc_swifs(ni)%emis_flxsum, &
                lrestreq=.TRUE. )
           CALL channel_halt(substr, status)
           CALL new_attribute(status, modstr,  &
                TRIM(name), 'long_name', &
                c='accumulated Ocean POC emissions (based on SEAWIFS)' )
           CALL channel_halt(substr, status)
           CALL new_attribute(status, modstr, &
                TRIM(name), 'units', &
                c='kg (POC)  m-2' )
           CALL channel_halt(substr, status)
         ENDDO

       CASE('SS_WIOC_AQUA')

         if (associated(strname_cha)) DEALLOCATE (strname_cha)
         NULLIFY(strname_cha)
         if (associated(strname_obj)) DEALLOCATE (strname_obj)
         NULLIFY(strname_obj)
         CALL strcrack(imp_ss_wioc_aqua%cha,';', strname_cha, n_ss_wioc_aqua)
         CALL strcrack(imp_ss_wioc_aqua%obj,';', strname_obj, n_ss_wioc_aqua)

         ALLOCATE(emis_wioc_aqua(n_ss_wioc_aqua))

         DO ni=1,n_ss_wioc_aqua
           emis_wioc_aqua(ni)%channel= strname_cha(ni)
           emis_wioc_aqua(ni)%object = strname_obj(ni)

           name='emis_wioc_aqua_'//TRIM(emis_wioc_aqua(ni)%object)
           CALL new_channel_object(status, modstr, &
                TRIM(name), p2=emis_wioc_aqua(ni)%emis_flx, lrestreq=.TRUE. )
           CALL channel_halt(substr, status)
           CALL new_attribute(status, modstr,  &
                TRIM(name), 'long_name', &
                c='Ocean WIOC emissions (based on AQUA)' )
           CALL channel_halt(substr, status)
           CALL new_attribute(status, modstr, &
                TRIM(name), 'units', &
                  c='kg (WIOC)  m-2 s-1' )
           CALL channel_halt(substr, status)

           name='emissum_wioc_aqua_'//TRIM(emis_wioc_aqua(ni)%object)
           CALL new_channel_object(status, modstr, &
                TRIM(name), p2=emis_wioc_aqua(ni)%emis_flxsum, &
                lrestreq=.TRUE. )
           CALL channel_halt(substr, status)
           CALL new_attribute(status, modstr,  &
                TRIM(name), 'long_name', &
                c='accumulated Ocean WIOC emissions (based on AQUA)' )
           CALL channel_halt(substr, status)
           CALL new_attribute(status, modstr, &
                TRIM(name), 'units', &
                c='kg (WIOC)  m-2' )
           CALL channel_halt(substr, status)
         ENDDO

       CASE('SS_WIOC_BLEN')
         if (associated(strname_cha)) DEALLOCATE (strname_cha)
         NULLIFY(strname_cha)
         if (associated(strname_obj)) DEALLOCATE (strname_obj)
         NULLIFY(strname_obj)
         CALL strcrack(imp_ss_wioc_blen%cha,';', strname_cha, n_ss_wioc_blend)
         CALL strcrack(imp_ss_wioc_blen%obj,';', strname_obj, n_ss_wioc_blend)

         ALLOCATE(emis_wioc_blend(n_ss_wioc_blend))

         DO ni=1,n_ss_wioc_blend
           emis_wioc_blend(ni)%channel= strname_cha(ni)
           emis_wioc_blend(ni)%object = strname_obj(ni)

           name='emis_wioc_blend_'//TRIM(emis_wioc_blend(ni)%object)
           CALL new_channel_object(status, modstr, &
                TRIM(name), p2=emis_wioc_blend(ni)%emis_flx, lrestreq=.TRUE. )
           CALL channel_halt(substr, status)
           CALL new_attribute(status, modstr,  &
                TRIM(name), 'long_name', &
                c='Ocean WIOC emissions (based on CZCS-blend)' )
           CALL channel_halt(substr, status)
           CALL new_attribute(status, modstr, &
                TRIM(name), 'units', &
                  c='kg (WIOC)  m-2 s-1' )
           CALL channel_halt(substr, status)

           name='emissum_wioc_blend_'//TRIM(emis_wioc_blend(ni)%object)
           CALL new_channel_object(status, modstr, &
                TRIM(name), p2=emis_wioc_blend(ni)%emis_flxsum, &
                lrestreq=.TRUE. )
           CALL channel_halt(substr, status)
           CALL new_attribute(status, modstr,  &
                TRIM(name), 'long_name', &
                c='accumulated Ocean WIOC emissions (based on CZCS-blend)' )
           CALL channel_halt(substr, status)
           CALL new_attribute(status, modstr, &
                TRIM(name), 'units', &
                c='kg (WIOC)  m-2' )
           CALL channel_halt(substr, status)
         ENDDO

       CASE('OC/BC')
          CALL new_channel_object(status, modstr, &
               'BC_sum_insol', p2=BC_sum_insol )
          CALL channel_halt(substr, status)
          CALL new_attribute(status, modstr, &
               'BC_sum_insol', 'long_name',  &
               c='black carbon mass insoluble emission flux' )
          CALL channel_halt(substr, status)
          CALL new_attribute(status, modstr, &
               'BC_sum_insol', 'units', c='kg (BC) m-2 s-1' )
          CALL channel_halt(substr, status)

          CALL new_channel_object(status, modstr, &
               'OC_sum_insol', p2=OC_sum_insol )
          CALL channel_halt(substr, status)
          CALL new_attribute(status, modstr, &
               'OC_sum_insol', 'long_name',  &
               c='organic carbon mass insoluble emission flux' )
          CALL channel_halt(substr, status)
          CALL new_attribute(status, modstr, &
               'OC_sum_insol', 'units', c='kg (OC) m-2 s-1' )
          CALL channel_halt(substr, status)

          CALL new_channel_object(status, modstr, &
               'OC_sum_sol', p2=OC_sum_sol )
          CALL channel_halt(substr, status)
          CALL new_attribute(status, modstr, &
               'OC_sum_sol', 'long_name',    &
               c='organic carbon mass soluble emission flux' )
          CALL channel_halt(substr, status)
          CALL new_attribute(status, modstr, &
               'OC_sum_sol', 'units', c='kg (OC) m-2 s-1' )
          CALL channel_halt(substr, status)

          CALL new_channel_object(status, modstr, 'Num_sol', p2=N_sol )
          CALL channel_halt(substr, status)
          CALL new_attribute(status, modstr, &
               'Num_sol', 'units', c='m-2 s-1 *kg(aerosol)/m3(aerosol)' )
          CALL channel_halt(substr, status)

          CALL new_channel_object(status, modstr, &
               'Num_insol', p2=N_insol )
          CALL channel_halt(substr, status)
          CALL new_attribute(status, modstr, &
               'Num_insol', 'units', c='m-2 s-1 *kg(aerosol)/m3(aerosol)' )
          CALL channel_halt(substr, status)

          CALL new_channel_object(status, modstr, &
               'OC_soa_sol', p2=OC_soa_sol )
          CALL channel_halt(substr, status)
          CALL new_attribute(status, modstr, &
               'OC_soa_sol', 'long_name',    &
               c='organic carbon mass SOA soluble emission flux' )
          CALL channel_halt(substr, status)
          CALL new_attribute(status, modstr, &
               'OC_soa_sol', 'units', c='kg (OC) m-2 s-1' )
          CALL channel_halt(substr, status)

          CALL new_channel_object(status, modstr, &
               'OC_bb_sol', p2=OC_bb_sol )
          CALL channel_halt(substr, status)
          CALL new_attribute(status, modstr, &
               'OC_bb_sol', 'long_name',    &
               c='organic carbon mass BB soluble emission flux' )
          CALL channel_halt(substr, status)
          CALL new_attribute(status, modstr, &
               'OC_bb_sol', 'units', c='kg (OC) m-2 s-1' )
          CALL channel_halt(substr, status)

          CALL new_channel_object(status, modstr, &
               'OC_ff_sol', p2=OC_ff_sol )
          CALL channel_halt(substr, status)
          CALL new_attribute(status, modstr, &
               'OC_ff_sol', 'long_name',    &
               c='organic carbon mass FF soluble emission flux' )
          CALL channel_halt(substr, status)
          CALL new_attribute(status, modstr, &
               'OC_ff_sol', 'units', c='kg (OC) m-2 s-1' )
          CALL channel_halt(substr, status)

          CALL new_channel_object(status, modstr, &
               'OC_soa_insol', p2=OC_soa_insol )
          CALL channel_halt(substr, status)
          CALL new_attribute(status, modstr, &
               'OC_soa_insol', 'long_name',    &
               c='organic carbon mass SOA insoluble emission flux' )
          CALL channel_halt(substr, status)
          CALL new_attribute(status, modstr, &
               'OC_soa_insol', 'units', c='kg (OC) m-2 s-1' )
          CALL channel_halt(substr, status)

          CALL new_channel_object(status, modstr, &
               'OC_bb_insol', p2=OC_bb_insol )
          CALL channel_halt(substr, status)
          CALL new_attribute(status, modstr, &
               'OC_bb_insol', 'long_name',    &
               c='organic carbon mass BB insoluble emission flux' )
          CALL channel_halt(substr, status)
          CALL new_attribute(status, modstr, &
               'OC_bb_insol', 'units', c='kg (OC) m-2 s-1' )
          CALL channel_halt(substr, status)

          CALL new_channel_object(status, modstr, &
               'OC_ff_insol', p2=OC_ff_insol )
          CALL channel_halt(substr, status)
          CALL new_attribute(status, modstr, &
               'OC_ff_insol', 'long_name',    &
               c='organic carbon mass FF insoluble emission flux' )
          CALL channel_halt(substr, status)
          CALL new_attribute(status, modstr, &
               'OC_ff_insol', 'units', c='kg (OC) m-2 s-1' )
          CALL channel_halt(substr, status)

          CALL new_channel_object(status, modstr, &
               'BC_bb_insol', p2=BC_bb_insol )
          CALL channel_halt(substr, status)
          CALL new_attribute(status, modstr, &
               'BC_bb_insol', 'long_name',    &
               c='Black carbon mass BB insoluble emission flux' )
          CALL channel_halt(substr, status)
          CALL new_attribute(status, modstr, &
               'BC_bb_insol', 'units', c='kg (OC) m-2 s-1' )
          CALL channel_halt(substr, status)

          CALL new_channel_object(status, modstr, &
               'BC_ff_insol', p2=BC_ff_insol )
          CALL channel_halt(substr, status)
          CALL new_attribute(status, modstr, &
               'BC_ff_insol', 'long_name',    &
               c='Black carbon mass FF insoluble emission flux' )
          CALL channel_halt(substr, status)
          CALL new_attribute(status, modstr, &
               'BC_ff_insol', 'units', c='kg (OC) m-2 s-1' )
          CALL channel_halt(substr, status)

          CALL new_channel_object(status, modstr, &
               'Num_insol_bc', p2=N_insol_bc )
          CALL channel_halt(substr, status)
          CALL new_attribute(status, modstr, &
               'Num_insol_bc', 'long_name',    &
               c='Black carbon insoluble number flux' )
          CALL channel_halt(substr, status)
          CALL new_attribute(status, modstr, &
               'Num_insol_bc', 'units', c='m-2 s-1 *kg(aerosol)/m3(aerosol)' )
          CALL channel_halt(substr, status)

          CALL new_channel_object(status, modstr, &
               'Num_insol_oc', p2=N_insol_oc )
          CALL channel_halt(substr, status)
          CALL new_attribute(status, modstr, &
               'Num_insol_oc', 'long_name',    &
               c='organic carbon insoluble number flux' )
          CALL channel_halt(substr, status)
          CALL new_attribute(status, modstr, &
               'Num_insol_oc', 'units', c='m-2 s-1 *kg(aerosol)/m3(aerosol)' )
          CALL channel_halt(substr, status)

       CASE ('CH4')
          CALL new_channel_object(status, modstr, &
               'CH4_emflux', p2=CH4_emflux )
          CALL channel_halt(substr, status)
          CALL new_attribute(status, modstr, &
               'CH4_emflux', 'long_name', c='CH4 emission flux' )
          CALL channel_halt(substr, status)
          CALL new_attribute(status, modstr, &
               'CH4_emflux', 'units', c='molec. m-2 s-1' )
          CALL channel_halt(substr, status)

       CASE ('VOC')
          CALL new_channel_object(status, modstr, &
               'ISOP_emflux', p2=ISOP_emflux )
          CALL channel_halt(substr, status)
          CALL new_attribute(status, modstr, &
               'ISOP_emflux', 'long_name', c='isoprene emission flux' )
          CALL channel_halt(substr, status)
          CALL new_attribute(status, modstr, &
               'ISOP_emflux', 'units', c='molec. m-2 s-1' )
          CALL channel_halt(substr, status)

          CALL new_channel_object(status, modstr, &
               'MTERP_emflux', p2=MTERP_emflux )
          CALL channel_halt(substr, status)
          CALL new_attribute(status, modstr, &
               'MTERP_emflux', 'long_name', c='monoterpene emission flux' )
          CALL channel_halt(substr, status)
          CALL new_attribute(status, modstr, &
               'MTERP_emflux', 'units', c='molec. m-2 s-1' )
          CALL channel_halt(substr, status)

          CALL new_channel_object(status, modstr, 'z0m', p2=z0m )
          CALL channel_halt(substr, status)
          CALL new_attribute(status, modstr, &
               'z0m', 'long_name', c='surface roughness' )
          CALL channel_halt(substr, status)
          CALL new_attribute(status, modstr, 'z0m', 'units', c='m' )
          CALL channel_halt(substr, status)

          ! NEW DIMENSION
          CALL new_dimension(status, DIMID_numveglay_hr, &
               'numveglay_hr', nveglay_hr)
          CALL channel_halt(substr, status)
          ! NEW REPRESENTATION
          CALL new_representation(status, REPRID_numveglay_hr,             &
               'REPRID_numveglay_hr'                                       &
               , rank = 3, link = 'xxx-', dctype = DC_GP                   &
               , dimension_ids = (/&
                 _RI_XY_N_(DIMID_LON, DIMID_LAT, DIMID_numveglay_hr) /)    &
               , ldimlen       = (/ _RI_XY_N_(nproma, ngpblks, AUTO)  /)   &
               , output_order  = (/ _IX_XY_N_ , _IY_XY_N_ , _IN_XY_N_ /)   &
               , axis = repr_def_axes(_RI_XY_N_('X','Y','N'),'-')          &
               )
          CALL channel_halt(substr, status)

          nseg = gp_nseg
          ALLOCATE(start(nseg,IRANK))
          ALLOCATE(cnt(nseg,IRANK))
          ALLOCATE(meml(nseg,IRANK))
          ALLOCATE(memu(nseg,IRANK))

          start(:,1) = gp_start(:,1)
          cnt(:,1) = gp_cnt(:,1)
          meml(:,1) = gp_meml(:,1)
          memu(:,1) = gp_memu(:,1)

          start(:,_IY_XY_N_) = gp_start(:,_IY_XYZN_)
          cnt(:,_IY_XY_N_)   = gp_cnt(:,_IY_XYZN_)
          meml(:,_IY_XY_N_)  = gp_meml(:,_IY_XYZN_)
          memu(:,_IY_XY_N_)  = gp_memu(:,_IY_XYZN_)

          start(:,_IN_XY_N_) = 1
          cnt(:,_IN_XY_N_) = nveglay_hr
          meml(:,_IN_XY_N_) = 1
          memu(:,_IN_XY_N_) = nveglay_hr
          start(:,4) = 1
          cnt(:,4) = 1
          meml(:,4) = 1
          memu(:,4) = 1

          CALL set_representation_decomp(status, REPRID_numveglay_hr &
               , start, cnt, memu, meml, .FALSE., PIOTYPE_COL)
          CALL channel_halt(substr, status)

          DEALLOCATE(start) ; NULLIFY(start)
          DEALLOCATE(cnt)   ; NULLIFY(cnt)
          DEALLOCATE(meml)  ; NULLIFY(meml)
          DEALLOCATE(memu)  ; NULLIFY(memu)

          CALL new_channel_object(status, modstr, 'lad', p3=lad &
               , reprid=REPRID_numveglay_hr )
          CALL channel_halt(substr, status)
          CALL new_attribute(status, modstr, &
               'lad', 'long_name', c='Leaf Area Density' )
          CALL channel_halt(substr, status)
          CALL new_attribute(status, modstr, &
               'lad', 'units', c='m2 m-3' )
          CALL channel_halt(substr, status)

       CASE('NO')
          CALL new_channel_object(status, modstr, 'NO_emflux', p2=NO_emflux )
          CALL channel_halt(substr, status)
          CALL new_attribute(status, modstr, &
               'NO_emflux', 'long_name',     &
               c='soil-biogenic NO emission fluxes' )
          CALL channel_halt(substr, status)
          CALL new_attribute(status, modstr, &
               'NO_emflux', 'units', c='molec. m-2 s-1' )
          CALL channel_halt(substr, status)

          ! NEW DIMENSION
          CALL new_dimension(status, DIMID_ncl_yl95, 'ncl_yl95', ncl_yl95)
          CALL channel_halt(substr, status)
          ! NEW REPRESENTATION
          CALL new_representation(status, REPRID_ncl_yl95, 'REPRID_ncl_yl95' &
               , rank = 3, link = 'xxx-', dctype = DC_GP                     &
               , dimension_ids = (/ &
                       _RI_XY_N_(DIMID_LON, DIMID_LAT,DIMID_ncl_yl95) /)     &
               , ldimlen       = (/ _RI_XY_N_(nproma, ngpblks, AUTO)  /)     &
               , output_order  = (/ _IX_XY_N_ , _IY_XY_N_ , _IN_XY_N_ /)     &
               , axis = repr_def_axes(_RI_XY_N_('X','Y','N'),'-')            &
               )
          CALL channel_halt(substr, status)

          nseg = gp_nseg
          ALLOCATE(start(nseg,IRANK))
          ALLOCATE(cnt(nseg,IRANK))
          ALLOCATE(meml(nseg,IRANK))
          ALLOCATE(memu(nseg,IRANK))

#ifndef VERTICO
          start(:,1) = gp_start(:,1)
          cnt(:,1)   = gp_cnt(:,1)
          meml(:,1)  = gp_meml(:,1)
          memu(:,1)  = gp_memu(:,1)

          start(:,_IN_XY_N_) = 1
          cnt(:,_IN_XY_N_)   = ncl_yl95
          meml(:,_IN_XY_N_)  = 1
          memu(:,_IN_XY_N_)  = ncl_yl95

          start(:,_IY_XY_N_) = gp_start(:,_IY_XYZN_)
          cnt(:,_IY_XY_N_)   = gp_cnt(:,_IY_XYZN_)
          meml(:,_IY_XY_N_)  = gp_meml(:,_IY_XYZN_)
          memu(:,_IY_XY_N_)  = gp_memu(:,_IY_XYZN_)
#endif
! not VERTICO
#ifdef VERTICO
          start(:,1) = 1
          cnt(:,1)   = 1
          meml(:,1)  = 1
          memu(:,1)  = 1

          start(:,2) = 1
          cnt(:,2)   = 1
          meml(:,2)  = 1
          memu(:,2)  = 1

          start(:,3) = 1
          cnt(:,3)   = ncl_yl95
          meml(:,3)  = 1
          memu(:,3)  = ncl_yl95
#endif
          start(:,4) = 1
          cnt(:,4)   = 1
          meml(:,4)  = 1
          memu(:,4)  = 1

          CALL set_representation_decomp(status, REPRID_ncl_yl95 &
               , start, cnt, memu, meml, .FALSE., PIOTYPE_COL)
          CALL channel_halt(substr, status)

          DEALLOCATE(start) ; NULLIFY(start)
          DEALLOCATE(cnt)   ; NULLIFY(cnt)
          DEALLOCATE(meml)  ; NULLIFY(meml)
          DEALLOCATE(memu)  ; NULLIFY(memu)

          CALL new_channel_object(status, modstr, &
               'noslflux', p2=noslflux, reprid=GP_2D_HORIZONTAL)
          CALL channel_halt(substr, status)
          CALL new_attribute(status, modstr, &
               'noslflux', 'long_name'       &
               , c='NO soil emission flux' )
          CALL channel_halt(substr, status)
          CALL new_attribute(status, modstr, 'noslflux' &
               , 'units', c='ng(N)/(m^2*s)' )
          CALL channel_halt(substr, status)

          CALL new_channel_object(status, modstr, &
               'noslflux_diag', p3=noslflux_diag, reprid=REPRID_ncl_yl95 )
          CALL channel_halt(substr, status)
          CALL new_attribute(status, modstr, &
               'noslflux_diag', 'long_name'  &
               , c='NO soil emission flux per ecosystem' )
          CALL channel_halt(substr, status)
          CALL new_attribute(status, modstr, 'noslflux_diag' &
               , 'units', c='ng(N)/(m^2*s)' )
          CALL channel_halt(substr, status)

          CALL new_channel_object(status, modstr, 'NOemis_w', p2=NOemis_w )
          CALL channel_halt(substr, status)
          CALL new_attribute(status, modstr, &
               'NOemis_w', 'long_name', c='NO emission factor, wet soils' )
          CALL channel_halt(substr, status)
          CALL new_attribute(status, modstr, &
               'NOemis_w', 'units', c='ng(N) m-2 s-1' )
          CALL channel_halt(substr, status)

          CALL new_channel_object(status, modstr, 'NOemis_d', p2=NOemis_d )
          CALL channel_halt(substr, status)
          CALL new_attribute(status, modstr, &
               'NOemis_d', 'long_name', c='NO emission factor, dry soils' )
          CALL channel_halt(substr, status)
          CALL new_attribute(status, modstr, &
               'NOemis_d', 'units', c='ng(N) m-2 s-1' )
          CALL channel_halt(substr, status)

       CASE('NOpls')
          ! Store the information of the status of the pulse in channel,
          ! so that it is known even for a change in month/rerun simulation

          CALL get_representation_id(status, 'ndrydays', rep_id)
          IF ( (status /=0) .OR. (rep_id /= REPRID_ndrydays)) THEN
             CALL make_repr_ndrydays
          ENDIF

          CALL new_channel_object(status, modstr, 'cpold', p3=cpold &
               , lrestreq=.TRUE., reprid=REPRID_ndrydays )
          CALL channel_halt(substr, status)
          CALL new_attribute(status, modstr, &
               'cpold', 'long_name', c='conv. prec. record' )
          CALL channel_halt(substr, status)
          CALL new_attribute(status, modstr, 'cpold', 'units', c='m' )
          CALL channel_halt(substr, status)

          CALL new_channel_object(status, modstr, 'lspold', p3=lspold &
               , lrestreq=.TRUE., reprid=REPRID_ndrydays )
          CALL channel_halt(substr, status)
          CALL new_attribute(status, modstr, &
               'lspold', 'long_name', c='large scale. prec. record' )
          CALL channel_halt(substr, status)
          CALL new_attribute(status, modstr, 'lspold', 'units', c='m' )
          CALL channel_halt(substr, status)

          CALL new_channel_object(status, modstr, 'cp', p2=cp &
               , lrestreq = .TRUE.)
          CALL channel_halt(substr, status)
          CALL new_attribute(status, modstr, &
               'cp', 'long_name', c='daily accumulated conv. precip.' )
          CALL channel_halt(substr, status)
          CALL new_attribute(status, modstr, 'cp', 'units', c='m' )
          CALL channel_halt(substr, status)

          CALL new_channel_object(status, modstr, 'lsp', p2=lsp &
               , lrestreq = .TRUE.)
          CALL channel_halt(substr, status)
          CALL new_attribute(status, modstr, &
               'lsp', 'long_name', c='daily accumulated large scale precip.' )
          CALL channel_halt(substr, status)
          CALL new_attribute(status, modstr, 'lsp', 'units', c='m' )
          CALL channel_halt(substr, status)

          CALL new_channel_object(status, modstr, 'pulsing', p2=pulsing &
               , lrestreq = .TRUE.)
          CALL channel_halt(substr, status)
          CALL new_attribute(status, modstr, &
               'pulsing', 'long_name', c='pulsing regime' )
          CALL channel_halt(substr, status)
          CALL new_attribute(status, modstr, &
               'pulsing', 'units', c='Index [1-3]' )
          CALL channel_halt(substr, status)

          CALL new_channel_object(status, modstr, 'plsday', p2=plsday &
               , lrestreq = .TRUE.)
          CALL channel_halt(substr, status)
          CALL new_attribute(status, modstr, &
               'plsday', 'long_name', c='timing of pulse' )
          CALL channel_halt(substr, status)
          CALL new_attribute(status, modstr, 'plsday', 'units', c='-' )
          CALL channel_halt(substr, status)

          CALL new_channel_object(status, modstr, 'plsdurat', p2=plsdurat &
               , lrestreq = .TRUE.)
          CALL channel_halt(substr, status)
          CALL new_attribute(status, modstr, &
               'plsdurat', 'long_name', c='duration of pulse' )
          CALL channel_halt(substr, status)
          CALL new_attribute(status, modstr, &
               'plsdurat', 'units', c='days' )
          CALL channel_halt(substr, status)

          CALL new_channel_object(status, modstr, 'pls', p2=pls &
               , lrestreq = .TRUE.)
          CALL channel_halt(substr, status)
          CALL new_attribute(status, modstr, 'pls', 'long_name', c='pulse' )
          CALL channel_halt(substr, status)
          CALL new_attribute(status, modstr, 'pls', 'units', c='1 -  ' )
          CALL channel_halt(substr, status)

          CALL new_channel_object(status, modstr, 'prectot', p2=prectot, &
               lrestreq = .TRUE.)
          CALL channel_halt(substr, status)
          CALL new_attribute(status, modstr, &
               'prectot', 'long_name', c='monthly mean total precipitation' )
          CALL channel_halt(substr, status)
          CALL new_attribute(status, modstr, &
               'prectot', 'units', c='m' )
          CALL channel_halt(substr, status)

          IF (lstart) THEN
             cp(:,:)=0._dp
             lsp(:,:)=0._dp

             plsday(:,:)=0._dp
             pulsing(:,:)=0._dp
             plsdurat(:,:)=0._dp
             cpold(:,:,:)=4.e-4_dp
             lspold(:,:,:)=4.e-4_dp
             ! for initialization a small value has been
             ! selected to avoid the calculation of pulse for first timesteps.
             ! The value 4e-4 results to an accumulated amount of precipitation
             ! of ndrydays*(4e-3+4e-3)*1000.=11.2 mm of rainfall
          ENDIF

       CASE('NO_yl95sl10')

          CALL get_representation_id(status, 'ndrydays', rep_id)
          IF ( (status /=0) .OR. (rep_id /= REPRID_ndrydays)) THEN
             CALL make_repr_ndrydays
          ENDIF

          ! NEW DIMENSION
          CALL new_dimension(status, DIMID_ncl_yl95sl10, &
               'nyl95sl10', ncl_yl95sl10)
          CALL channel_halt(substr, status)
          ! NEW REPRESENTATION
          CALL new_representation(status, REPRID_ncl_yl95sl10                  &
               , 'REPRID_nyl95sl10', rank = 3, link = 'xxx-', dctype = DC_GP   &
               , dimension_ids = (/ &
                       _RI_XY_N_(DIMID_LON, DIMID_LAT, DIMID_ncl_yl95sl10) /)  &
               , ldimlen       = (/ _RI_XY_N_(nproma, ngpblks, AUTO)   /)      &
               , output_order  = (/ _IX_XY_N_ , _IY_XY_N_ , _IN_XY_N_ /)   &
               , axis = repr_def_axes(_RI_XY_N_('X','Y','N'),'-')          &
               )
          CALL channel_halt(substr, status)
          nseg = gp_nseg
          ALLOCATE(start(nseg,IRANK))
          ALLOCATE(cnt(nseg,IRANK))
          ALLOCATE(meml(nseg,IRANK))
          ALLOCATE(memu(nseg,IRANK))

          start(:,1) = gp_start(:,1)
          cnt(:,1)   = gp_cnt(:,1)
          meml(:,1)  = gp_meml(:,1)
          memu(:,1)  = gp_memu(:,1)

          start(:,_IY_XY_N_) = gp_start(:,_IY_XYZN_)
          cnt(:,_IY_XY_N_)   = gp_cnt(:,_IY_XYZN_)
          meml(:,_IY_XY_N_)  = gp_meml(:,_IY_XYZN_)
          memu(:,_IY_XY_N_)  = gp_memu(:,_IY_XYZN_)

          start(:,_IN_XY_N_) = 1
          cnt(:,_IN_XY_N_)   = ncl_yl95sl10
          meml(:,_IN_XY_N_)  = 1
          memu(:,_IN_XY_N_)  = ncl_yl95sl10

          start(:,4) = 1
          cnt(:,4)   = 1
          meml(:,4)  = 1
          memu(:,4)  = 1

          CALL set_representation_decomp(status, REPRID_ncl_yl95sl10 &
               , start, cnt, memu, meml, .FALSE., PIOTYPE_COL)
          CALL channel_halt(substr, status)

          DEALLOCATE(start) ; NULLIFY(start)
          DEALLOCATE(cnt)   ; NULLIFY(cnt)
          DEALLOCATE(meml)  ; NULLIFY(meml)
          DEALLOCATE(memu)  ; NULLIFY(memu)

          ! added variables for the orignal Yienger & Levy algorithm
          ! NOemclass_yl95sl10
          CALL new_channel_object(status, modstr, 'pulse', p2=pulse &
               , lrestreq = .TRUE.)
          CALL channel_halt(substr, status)
          CALL new_attribute(status, modstr, 'pulse' &
               , 'long_name', c='pulsing factor')
          CALL channel_halt(substr, status)
          CALL new_attribute(status, modstr,'pulse' , 'units', c='1 -  ' )
          CALL channel_halt(substr, status)

          CALL new_channel_object(status, modstr &
               , 'pulseregime', p2=pulseregime &
               , lrestreq = .TRUE.)
          CALL channel_halt(substr, status)
          CALL new_attribute(status, modstr, 'pulseregime' &
               , 'long_name', c='pulsing regime' )
          CALL channel_halt(substr, status)
          CALL new_attribute(status, modstr, 'pulseregime' &
               , 'units', c='Index [0-3]' )
          CALL channel_halt(substr, status)

          CALL new_channel_object(status, modstr,'pulseday' , p2=pulseday &
               , lrestreq = .TRUE.)
          CALL channel_halt(substr, status)
          CALL new_attribute(status, modstr,'pulseday' &
               , 'long_name', c='duration of pulsing' )
          CALL channel_halt(substr, status)
          CALL new_attribute(status, modstr,'pulseday' , 'units', c='days' )
          CALL channel_halt(substr, status)

          CALL new_channel_object(status, modstr,'prec_hist' , p3=prec_hist &
               , lrestreq = .TRUE., reprid=REPRID_ndrydays)
          CALL channel_halt(substr, status)
          CALL new_attribute(status, modstr, 'prec_hist'&
               , 'long_name', c='daily precipitation history' )
          CALL channel_halt(substr, status)
          CALL new_attribute(status, modstr, 'prec_hist', 'units', c='m' )
          CALL channel_halt(substr, status)

          CALL new_channel_object(status, modstr,'crf' , p2=crf)
          CALL channel_halt(substr, status)
          CALL new_attribute(status, modstr, 'crf'&
               , 'long_name', c='canopy reduction factor' )
          CALL channel_halt(substr, status)
          CALL new_attribute(status, modstr, 'crf', 'units', c='0 - 1' )
          CALL channel_halt(substr, status)

          CALL new_channel_object(status, modstr,'noslflux_yl95sl10' &
               , p2=noslflux_yl95sl10)
          CALL channel_halt(substr, status)
          CALL new_attribute(status, modstr, 'noslflux_yl95sl10'&
               , 'long_name', c='NO soil emission flux' )
          CALL channel_halt(substr, status)
          CALL new_attribute(status, modstr,'noslflux_yl95sl10' &
               , 'units', c='ng(N)/(m^2*s)' )
          CALL channel_halt(substr, status)

          CALL new_channel_object(status, modstr, 'noemflux_yl95sl10' &
               , p2=noemflux_yl95sl10)
          CALL channel_halt(substr, status)
          CALL new_attribute(status, modstr, 'noemflux_yl95sl10'&
               , 'long_name', c='above canopy NO emission flux from soils' )
          CALL channel_halt(substr, status)
          CALL new_attribute(status, modstr,'noemflux_yl95sl10' &
               , 'units', c='molec/(m^2*s)' )
          CALL channel_halt(substr, status)

          CALL new_channel_object(status, modstr,'noslflux_diag_yl95sl10' &
               , p3=noslflux_diag_yl95sl10, reprid=REPRID_ncl_yl95sl10)
          CALL channel_halt(substr, status)
          CALL new_attribute(status, modstr, 'noslflux_diag_yl95sl10'&
               , 'long_name', c='NO soil emission flux per landcover' )
          CALL channel_halt(substr, status)
          CALL new_attribute(status, modstr,'noslflux_diag_yl95sl10' &
               , 'units', c='ng(N)/(m^2*s)' )
          CALL channel_halt(substr, status)

          CALL new_channel_object(status, modstr,'tsoil_top' &
               , p2=tsoil_top, reprid=GP_2D_HORIZONTAL)
          CALL channel_halt(substr, status)
          CALL new_attribute(status, modstr, 'tsoil_top'&
               , 'long_name', c='soil temperature (top layer)' )
          CALL channel_halt(substr, status)
          CALL new_attribute(status, modstr,'tsoil_top' &
               , 'units', c='-' )
          CALL channel_halt(substr, status)

          CALL new_channel_object(status, modstr,'vsm' &
               , p2=vsm, reprid=GP_2D_HORIZONTAL)
          CALL channel_halt(substr, status)
          CALL new_attribute(status, modstr, 'vsm'&
               , 'long_name', c='volumetric soil moisture' )
          CALL channel_halt(substr, status)
          CALL new_attribute(status, modstr,'vsm' &
               , 'units', c='-' )
          CALL channel_halt(substr, status)

       IF (lstart) THEN
          prec_hist(:,:,:) = 0.011_dp
          pulse(:,:)       = 1.
          pulseregime(:,:) = 0.
       ENDIF

       CASE('DU')
          CALL new_channel_object(status, modstr, &
               'du_emflux_B_ci', p2=du_emflux_B_ci )
          CALL channel_halt(substr, status)
          CALL new_attribute(status, modstr, &
               'du_emflux_B_ci', 'long_name', c='Dust emission flux' )
          CALL channel_halt(substr, status)
          CALL new_attribute(status, modstr, &
               'du_emflux_B_ci', 'units', c='kg m-2 s-1' )
          CALL channel_halt(substr, status)

          CALL new_channel_object(status, modstr, &
               'zprecipinsoil', p2=zprecipinsoil_2d, lrestreq = .TRUE.)
          CALL channel_halt(substr, status)
          CALL new_attribute(status, modstr, &
               'zprecipinsoil', 'long_name', &
               c='soil moisture (below surface)' )
          CALL channel_halt(substr, status)
          CALL new_attribute(status, modstr, &
               'zprecipinsoil', 'units', c='m' )
          CALL channel_halt(substr, status)

       CASE('DU_tegen')
          CALL new_channel_object(status, modstr, &
               'du_emflux_T_ci', p2=du_emflux_T_ci )
          CALL channel_halt(substr, status)
          CALL new_attribute(status, modstr, &
               'du_emflux_T_ci', 'long_name', c='Dust emission flux ci' )
          CALL channel_halt(substr, status)
          CALL new_attribute(status, modstr, &
               'du_emflux_T_ci', 'units', c='kg m-2 s-1' )
          CALL channel_halt(substr, status)

          CALL new_channel_object(status, modstr, &
               'du_emflux_T_ai', p2=du_emflux_T_ai )
          CALL channel_halt(substr, status)
          CALL new_attribute(status, modstr, &
               'du_emflux_T_ai', 'long_name', c='Dust emission flux ai' )
          CALL channel_halt(substr, status)
          CALL new_attribute(status, modstr, &
               'du_emflux_T_ai', 'units', c='kg m-2 s-1' )
          CALL channel_halt(substr, status)

       CASE('DU_Astitha1')
          CALL new_channel_object(status, modstr, &
               'du_emflux_A1_ci', p2=du_emflux_A1_ci )
          CALL channel_halt(substr, status)
          CALL new_attribute(status, modstr,   &
               'du_emflux_A1_ci', 'long_name', &
               c='Dust emission vertical flux- ci' )
          CALL channel_halt(substr, status)
          CALL new_attribute(status, modstr, &
               'du_emflux_A1_ci', 'units', c='kg m-2 s-1' )
          CALL channel_halt(substr, status)

          CALL new_channel_object(status, modstr, &
               'du_emflux_A1_ai', p2=du_emflux_A1_ai )
          CALL channel_halt(substr, status)
          CALL new_attribute(status, modstr,   &
               'du_emflux_A1_ai', 'long_name', &
               c='Dust emission vertical flux-ai' )
          CALL channel_halt(substr, status)
          CALL new_attribute(status, modstr, &
               'du_emflux_A1_ai', 'units', c='kg m-2 s-1' )
          CALL channel_halt(substr, status)

          CALL new_channel_object(status, modstr, &
               'horflux', p2=horflux )
          CALL channel_halt(substr, status)
          CALL new_attribute(status, modstr, &
               'horflux', 'long_name',       &
               c='Dust emission horiz. flux-DU_Astitha2' )
          CALL channel_halt(substr, status)
          CALL new_attribute(status, modstr, &
               'horflux', 'units', c='kg m-1 s' )
          CALL channel_halt(substr, status)


          CALL new_channel_object(status, modstr, &
               'ustarthr1', p2=ustarthr1, lrestreq = .TRUE. )
          CALL channel_halt(substr, status)
          CALL new_attribute(status, modstr, 'ustarthr1'&
               , 'long_name', c='threshold wind friction velocity-DU_Astitha2' )
          CALL channel_halt(substr, status)
          CALL new_attribute(status, modstr, &
               'ustarthr1', 'units', c='m s-1' )
          CALL channel_halt(substr, status)

          if (l_ducomp) then
             CALL new_channel_object(status, modstr, 'du_nap_emflux_ai', &
                                     p2 = du_nap_emflux_ai)
             CALL channel_halt(substr, status)
             CALL new_attribute(status, modstr, 'du_nap_emflux_ai' &
                  , 'long_name',  c = 'Dust Na+ emission flux, acc. mode')
             CALL channel_halt(substr, status)
             CALL new_attribute(status, modstr, 'du_nap_emflux_ai', 'units', &
                                c = 'kg m-2 s-1' )
             CALL channel_halt(substr, status)

             CALL new_channel_object(status, modstr, 'du_nap_emflux_ci', &
                                     p2 = du_nap_emflux_ci)
             CALL channel_halt(substr, status)
             CALL new_attribute(status, modstr, 'du_nap_emflux_ci', 'long_name',&
                                c = 'Dust Na+ emission flux, coarse mode')
             CALL channel_halt(substr, status)
             CALL new_attribute(status, modstr, 'du_nap_emflux_ci', 'units', &
                                c = 'kg m-2 s-1' )
             CALL channel_halt(substr, status)

             CALL new_channel_object(status, modstr, 'du_kp_emflux_ai', &
                                     p2 = du_kp_emflux_ai)
             CALL channel_halt(substr, status)
             CALL new_attribute(status, modstr, 'du_kp_emflux_ai', 'long_name', &
                                c = 'Dust K+ emission flux, acc. mode')
             CALL channel_halt(substr, status)
             CALL new_attribute(status, modstr, 'du_kp_emflux_ai', 'units', &
                                c = 'kg m-2 s-1' )
             CALL channel_halt(substr, status)

             CALL new_channel_object(status, modstr, 'du_kp_emflux_ci', &
                                     p2 = du_kp_emflux_ci)
             CALL channel_halt(substr, status)
             CALL new_attribute(status, modstr, 'du_kp_emflux_ci', 'long_name', &
                                c = 'Dust K+ emission flux, coarse mode')
             CALL channel_halt(substr, status)
             CALL new_attribute(status, modstr, 'du_kp_emflux_ci', 'units', &
                                c = 'kg m-2 s-1' )
             CALL channel_halt(substr, status)

             CALL new_channel_object(status, modstr, 'du_capp_emflux_ai', &
                                     p2 = du_capp_emflux_ai)
             CALL channel_halt(substr, status)
             CALL new_attribute(status, modstr, 'du_capp_emflux_ai','long_name',&
                                c = 'Dust Ca++ emission flux, acc. mode')
             CALL channel_halt(substr, status)
             CALL new_attribute(status, modstr, 'du_capp_emflux_ai', 'units', &
                                c = 'kg m-2 s-1' )
             CALL channel_halt(substr, status)

             CALL new_channel_object(status, modstr, 'du_capp_emflux_ci', &
                                     p2 = du_capp_emflux_ci)
             CALL channel_halt(substr, status)
             CALL new_attribute(status, modstr, 'du_capp_emflux_ci', 'long_name', &
                                c = 'Dust Ca++ emission flux, coarse mode')
             CALL channel_halt(substr, status)
             CALL new_attribute(status, modstr, 'du_capp_emflux_ci', 'units', &
                                c = 'kg m-2 s-1' )
             CALL channel_halt(substr, status)

             CALL new_channel_object(status, modstr, 'du_mgpp_emflux_ai', &
                                     p2 = du_mgpp_emflux_ai)
             CALL channel_halt(substr, status)
             CALL new_attribute(status, modstr, 'du_mgpp_emflux_ai', &
                  'long_name', c = 'Dust Mg++ emission flux, acc. mode')
             CALL channel_halt(substr, status)
             CALL new_attribute(status, modstr, 'du_mgpp_emflux_ai', 'units', &
                                c = 'kg m-2 s-1' )
             CALL channel_halt(substr, status)

             CALL new_channel_object(status, modstr, 'du_mgpp_emflux_ci', &
                                     p2 = du_mgpp_emflux_ci)
             CALL channel_halt(substr, status)
             CALL new_attribute(status, modstr, 'du_mgpp_emflux_ci', &
                  'long_name',  c = 'Dust Mg++ emission flux, coarse mode')
             CALL channel_halt(substr, status)
             CALL new_attribute(status, modstr, 'du_mgpp_emflux_ci', 'units', &
                                c = 'kg m-2 s-1' )
             CALL channel_halt(substr, status)

             CALL new_channel_object(status, modstr, 'du_misc_emflux_ai', &
                                     p2 = du_misc_emflux_ai)
             CALL channel_halt(substr, status)
             CALL new_attribute(status, modstr, 'du_misc_emflux_ai', &
                  'long_name', c = 'Misc. dust emission flux, acc. mode')
             CALL channel_halt(substr, status)
             CALL new_attribute(status, modstr, 'du_misc_emflux_ai', 'units', &
                                c = 'kg m-2 s-1' )
             CALL channel_halt(substr, status)

             CALL new_channel_object(status, modstr, 'du_misc_emflux_ci', &
                                     p2 = du_misc_emflux_ci)
             CALL channel_halt(substr, status)
             CALL new_attribute(status, modstr, 'du_misc_emflux_ci', &
                  'long_name',  c = 'Misc. dust emission flux, coarse mode')
             CALL channel_halt(substr, status)
             CALL new_attribute(status, modstr, 'du_misc_emflux_ci', 'units', &
                                c = 'kg m-2 s-1' )
             CALL channel_halt(substr, status)
          endif

       CASE('DU_Astitha2')
          CALL new_channel_object(status, modstr, &
               'du_emflux_A2_ci', p2=du_emflux_A2_ci )
          CALL channel_halt(substr, status)
          CALL new_attribute(status, modstr, &
               'du_emflux_A2_ci', 'long_name', &
               c='Dust emission vertical flux- ci' )
          CALL channel_halt(substr, status)
          CALL new_attribute(status, modstr, &
               'du_emflux_A2_ci', 'units', c='kg m-2 s-1' )
          CALL channel_halt(substr, status)

          CALL new_channel_object(status, modstr, &
               'du_emflux_A2_ai', p2=du_emflux_A2_ai )
          CALL channel_halt(substr, status)
          CALL new_attribute(status, modstr, &
               'du_emflux_A2_ai', 'long_name', &
               c='Dust emission vertical flux-ai' )
          CALL channel_halt(substr, status)
          CALL new_attribute(status, modstr, &
               'du_emflux_A2_ai', 'units', c='kg m-2 s-1' )
          CALL channel_halt(substr, status)

          CALL new_dimension(status, DIMID_du3_maxsizes, 'maxsizes', maxsizes)
          CALL channel_halt(substr, status)
          CALL new_dimension(status, DIMID_du3_npt, 'DU3_npt', DU3_npt)
          CALL channel_halt(substr, status)
          ! NEW REPRESENTATION
          CALL new_representation(status, REPRID_du3_maxsizes                  &
              , 'REPRID_du3_maxsizes', rank = 3, link = 'xxx-', dctype = DC_GP &
              , dimension_ids = (/&
                     _RI_XY_N_( DIMID_LON, DIMID_LAT,DIMID_du3_maxsizes) /)    &
               , ldimlen       = (/ _RI_XY_N_(nproma, ngpblks,  AUTO) /)       &
               , output_order  = (/ _IX_XY_N_ , _IY_XY_N_ , _IN_XY_N_ /)   &
               , axis = repr_def_axes(_RI_XY_N_('X','Y','N'),'-')          &
               )
          CALL channel_halt(substr, status)
          nseg = gp_nseg
          ALLOCATE(start(nseg,IRANK))
          ALLOCATE(cnt(nseg,IRANK))
          ALLOCATE(meml(nseg,IRANK))
          ALLOCATE(memu(nseg,IRANK))

          start(:,1) = gp_start(:,1)
          cnt(:,1)   = gp_cnt(:,1)
          meml(:,1)  = gp_meml(:,1)
          memu(:,1)  = gp_memu(:,1)

          start(:,_IY_XY_N_) = gp_start(:,_IY_XYZN_)
          cnt(:,_IY_XY_N_)   = gp_cnt(:,_IY_XYZN_)
          meml(:,_IY_XY_N_)  = gp_meml(:,_IY_XYZN_)
          memu(:,_IY_XY_N_)  = gp_memu(:,_IY_XYZN_)

          start(:,_IN_XY_N_) = 1
          cnt(:,_IN_XY_N_)   = maxsizes
          meml(:,_IN_XY_N_)  = 1
          memu(:,_IN_XY_N_)  = maxsizes

          start(:,4) = 1
          cnt(:,4)   = 1
          meml(:,4)  = 1
          memu(:,4)  = 1

          CALL set_representation_decomp(status, REPRID_du3_maxsizes &
               , start, cnt, memu, meml, .FALSE., PIOTYPE_COL)
          CALL channel_halt(substr, status)

          DEALLOCATE(start) ; NULLIFY(start)
          DEALLOCATE(cnt)   ; NULLIFY(cnt)
          DEALLOCATE(meml)  ; NULLIFY(meml)
          DEALLOCATE(memu)  ; NULLIFY(memu)

          CALL new_channel_object(status, modstr, 'ustarthr' &
               , p3=ustarthr, lrestreq = .TRUE., reprid=REPRID_du3_maxsizes)
          CALL channel_halt(substr, status)
          CALL new_attribute(status, modstr, 'ustarthr', 'long_name' &
               , c='threshold wind friction velocity per particle size-DU3' )
          CALL channel_halt(substr, status)
          CALL new_attribute(status, modstr, 'ustarthr', 'units', c='m s-1' )
          CALL channel_halt(substr, status)

          if (l_ducomp) then
             CALL new_channel_object(status, modstr, 'du_nap_emflux_ai', &
                                     p2 = du_nap_emflux_ai)
             CALL channel_halt(substr, status)
             CALL new_attribute(status, modstr, 'du_nap_emflux_ai', &
                  'long_name', c = 'Dust Na+ emission flux, acc. mode')
             CALL channel_halt(substr, status)
             CALL new_attribute(status, modstr, 'du_nap_emflux_ai', 'units', &
                                c = 'kg m-2 s-1' )
             CALL channel_halt(substr, status)

             CALL new_channel_object(status, modstr, 'du_nap_emflux_ci', &
                                     p2 = du_nap_emflux_ci)
             CALL channel_halt(substr, status)
             CALL new_attribute(status, modstr, 'du_nap_emflux_ci', &
                  'long_name', c = 'Dust Na+ emission flux, coarse mode')
             CALL channel_halt(substr, status)
             CALL new_attribute(status, modstr, 'du_nap_emflux_ci', 'units', &
                                c = 'kg m-2 s-1' )
             CALL channel_halt(substr, status)

             CALL new_channel_object(status, modstr, 'du_kp_emflux_ai', &
                                     p2 = du_kp_emflux_ai)
             CALL channel_halt(substr, status)
             CALL new_attribute(status, modstr, 'du_kp_emflux_ai', 'long_name', &
                                c = 'Dust K+ emission flux, acc. mode')
             CALL channel_halt(substr, status)
             CALL new_attribute(status, modstr, 'du_kp_emflux_ai', 'units', &
                                c = 'kg m-2 s-1' )
             CALL channel_halt(substr, status)

             CALL new_channel_object(status, modstr, 'du_kp_emflux_ci', &
                                     p2 = du_kp_emflux_ci)
             CALL channel_halt(substr, status)
             CALL new_attribute(status, modstr, 'du_kp_emflux_ci', 'long_name', &
                                c = 'Dust K+ emission flux, coarse mode')
             CALL channel_halt(substr, status)
             CALL new_attribute(status, modstr, 'du_kp_emflux_ci', 'units', &
                                c = 'kg m-2 s-1' )
             CALL channel_halt(substr, status)

             CALL new_channel_object(status, modstr, 'du_capp_emflux_ai', &
                                     p2 = du_capp_emflux_ai)
             CALL channel_halt(substr, status)
             CALL new_attribute(status, modstr, 'du_capp_emflux_ai', &
                  'long_name', c = 'Dust Ca++ emission flux, acc. mode')
             CALL channel_halt(substr, status)
             CALL new_attribute(status, modstr, 'du_capp_emflux_ai', 'units', &
                                c = 'kg m-2 s-1' )
             CALL channel_halt(substr, status)

             CALL new_channel_object(status, modstr, 'du_capp_emflux_ci', &
                                     p2 = du_capp_emflux_ci)
             CALL channel_halt(substr, status)
             CALL new_attribute(status, modstr, 'du_capp_emflux_ci',&
                  'long_name', c = 'Dust Ca++ emission flux, coarse mode')
             CALL channel_halt(substr, status)
             CALL new_attribute(status, modstr, 'du_capp_emflux_ci', 'units', &
                                c = 'kg m-2 s-1' )
             CALL channel_halt(substr, status)

             CALL new_channel_object(status, modstr, 'du_mgpp_emflux_ai', &
                                     p2 = du_mgpp_emflux_ai)
             CALL channel_halt(substr, status)
             CALL new_attribute(status, modstr, 'du_mgpp_emflux_ai', &
                  'long_name', c = 'Dust Mg++ emission flux, acc. mode')
             CALL channel_halt(substr, status)
             CALL new_attribute(status, modstr, 'du_mgpp_emflux_ai', 'units', &
                                c = 'kg m-2 s-1' )
             CALL channel_halt(substr, status)

             CALL new_channel_object(status, modstr, 'du_mgpp_emflux_ci', &
                                     p2 = du_mgpp_emflux_ci)
             CALL channel_halt(substr, status)
             CALL new_attribute(status, modstr, 'du_mgpp_emflux_ci', &
                  'long_name', c = 'Dust Mg++ emission flux, coarse mode')
             CALL channel_halt(substr, status)
             CALL new_attribute(status, modstr, 'du_mgpp_emflux_ci', 'units', &
                                c = 'kg m-2 s-1' )
             CALL channel_halt(substr, status)

             CALL new_channel_object(status, modstr, 'du_misc_emflux_ai', &
                                     p2 = du_misc_emflux_ai)
             CALL channel_halt(substr, status)
             CALL new_attribute(status, modstr, 'du_misc_emflux_ai',  &
                  'long_name', c = 'Misc. dust emission flux, acc. mode')
             CALL channel_halt(substr, status)
             CALL new_attribute(status, modstr, 'du_misc_emflux_ai', 'units', &
                                c = 'kg m-2 s-1' )
             CALL channel_halt(substr, status)

             CALL new_channel_object(status, modstr, 'du_misc_emflux_ci', &
                                     p2 = du_misc_emflux_ci)
             CALL channel_halt(substr, status)
             CALL new_attribute(status, modstr, 'du_misc_emflux_ci', &
                  'long_name', c = 'Misc. dust emission flux, coarse mode')
             CALL channel_halt(substr, status)
             CALL new_attribute(status, modstr, 'du_misc_emflux_ci', 'units', &
                                c = 'kg m-2 s-1' )
             CALL channel_halt(substr, status)
          endif

       CASE('KKDU_Astitha')
          CALL new_channel_object(status, modstr, 'du_emflux_A1_ci', &
                p2 = du_emflux_A1_ci)
          CALL channel_halt(substr, status)
          CALL new_attribute(status, modstr, 'du_emflux_A1_ci', 'long_name', &
                c = 'Dust emission vertical flux- ci' )
          CALL channel_halt(substr, status)
          CALL new_attribute(status, modstr, 'du_emflux_A1_ci', 'units', &
                c = 'kg m-2 s-1')
          CALL channel_halt(substr, status)

          CALL new_channel_object(status, modstr, 'du_emflux_A1_ai', &
                p2 = du_emflux_A1_ai)
          CALL channel_halt(substr, status)
          CALL new_attribute(status, modstr, 'du_emflux_A1_ai', 'long_name', &
                c = 'Dust emission vertical flux-ai')
          CALL channel_halt(substr, status)
          CALL new_attribute(status, modstr, 'du_emflux_A1_ai', 'units', &
                c = 'kg m-2 s-1')
          CALL channel_halt(substr, status)

          CALL new_channel_object(status, modstr, 'horflux', p2 = horflux)
          CALL channel_halt(substr, status)
          CALL new_attribute(status, modstr, 'horflux', 'long_name', &
                c = 'Dust emission horiz. flux-DU_Astitha2')
          CALL channel_halt(substr, status)
          CALL new_attribute(status, modstr, 'horflux', 'units', c = 'kg m-1 s')
          CALL channel_halt(substr, status)

          CALL new_channel_object(status, modstr, 'ustarthr1', p2 = ustarthr1, &
                lrestreq = .TRUE.)
          CALL channel_halt(substr, status)
          CALL new_attribute(status, modstr, 'ustarthr1', 'long_name', &
                c = 'threshold wind friction velocity-DU_Astitha2')
          CALL channel_halt(substr, status)
          CALL new_attribute(status, modstr, 'ustarthr1', 'units', c = 'm s-1')
          CALL channel_halt(substr, status)

          CALL new_channel_object(status, modstr, 'kkdu_nap_emflux_ai', &
                p2 = kkdu_nap_emflux_ai)
          CALL channel_halt(substr, status)
          CALL new_attribute(status, modstr, 'kkdu_nap_emflux_ai', &
                'long_name', c = 'Dust Na+ emission flux, acc. mode')
          CALL channel_halt(substr, status)
          CALL new_attribute(status, modstr, 'kkdu_nap_emflux_ai', 'units', &
                c = 'kg m-2 s-1' )
          CALL channel_halt(substr, status)

          CALL new_channel_object(status, modstr, 'kkdu_nap_emflux_ci', &
                p2 = kkdu_nap_emflux_ci)
          CALL channel_halt(substr, status)
          CALL new_attribute(status, modstr, 'kkdu_nap_emflux_ci', &
                'long_name', c = 'Dust Na+ emission flux, coarse mode')
          CALL channel_halt(substr, status)
          CALL new_attribute(status, modstr, 'kkdu_nap_emflux_ci', 'units', &
                c = 'kg m-2 s-1' )
          CALL channel_halt(substr, status)

          CALL new_channel_object(status, modstr, 'kkdu_kp_emflux_ai', &
                p2 = kkdu_kp_emflux_ai)
          CALL channel_halt(substr, status)
          CALL new_attribute(status, modstr, 'kkdu_kp_emflux_ai', 'long_name', &
                c = 'Dust K+ emission flux, acc. mode')
          CALL channel_halt(substr, status)
          CALL new_attribute(status, modstr, 'kkdu_kp_emflux_ai', 'units', &
                c = 'kg m-2 s-1' )
          CALL channel_halt(substr, status)

          CALL new_channel_object(status, modstr, 'kkdu_kp_emflux_ci', &
                p2 = kkdu_kp_emflux_ci)
          CALL channel_halt(substr, status)
          CALL new_attribute(status, modstr, 'kkdu_kp_emflux_ci', 'long_name', &
                c = 'Dust K+ emission flux, coarse mode')
          CALL channel_halt(substr, status)
          CALL new_attribute(status, modstr, 'kkdu_kp_emflux_ci', 'units', &
                c = 'kg m-2 s-1' )
          CALL channel_halt(substr, status)

          CALL new_channel_object(status, modstr, 'kkdu_capp_emflux_ai', &
                                  p2 = kkdu_capp_emflux_ai)
          CALL channel_halt(substr, status)
          CALL new_attribute(status, modstr, 'kkdu_capp_emflux_ai', &
                'long_name', c = 'Dust Ca++ emission flux, acc. mode')
          CALL channel_halt(substr, status)
          CALL new_attribute(status, modstr, 'kkdu_capp_emflux_ai', 'units', &
                c = 'kg m-2 s-1' )
          CALL channel_halt(substr, status)

          CALL new_channel_object(status, modstr, 'kkdu_capp_emflux_ci', &
                p2 = kkdu_capp_emflux_ci)
          CALL channel_halt(substr, status)
          CALL new_attribute(status, modstr, 'kkdu_capp_emflux_ci', &
                'long_name', c = 'Dust Ca++ emission flux, coarse mode')
          CALL channel_halt(substr, status)
          CALL new_attribute(status, modstr, 'kkdu_capp_emflux_ci', 'units', &
                c = 'kg m-2 s-1' )
          CALL channel_halt(substr, status)

          CALL new_channel_object(status, modstr, 'kkdu_mgpp_emflux_ai', &
                p2 = kkdu_mgpp_emflux_ai)
          CALL channel_halt(substr, status)
          CALL new_attribute(status, modstr, 'kkdu_mgpp_emflux_ai', &
                'long_name', c = 'Dust Mg++ emission flux, acc. mode')
          CALL channel_halt(substr, status)
          CALL new_attribute(status, modstr, 'kkdu_mgpp_emflux_ai', 'units', &
                c = 'kg m-2 s-1' )
          CALL channel_halt(substr, status)

          CALL new_channel_object(status, modstr, 'kkdu_mgpp_emflux_ci', &
                p2 = kkdu_mgpp_emflux_ci)
          CALL channel_halt(substr, status)
          CALL new_attribute(status, modstr, 'kkdu_mgpp_emflux_ci', &
                'long_name', c = 'Dust Mg++ emission flux, coarse mode')
          CALL channel_halt(substr, status)
          CALL new_attribute(status, modstr, 'kkdu_mgpp_emflux_ci', 'units', &
                c = 'kg m-2 s-1' )
          CALL channel_halt(substr, status)

          CALL new_channel_object(status, modstr, 'kkdu_misc_emflux_ai', &
                p2 = kkdu_misc_emflux_ai)
          CALL channel_halt(substr, status)
          CALL new_attribute(status, modstr, 'kkdu_misc_emflux_ai', &
                'long_name', c = 'Misc. dust emission flux, acc. mode')
          CALL channel_halt(substr, status)
          CALL new_attribute(status, modstr, 'kkdu_misc_emflux_ai', 'units', &
                c = 'kg m-2 s-1' )
          CALL channel_halt(substr, status)

          CALL new_channel_object(status, modstr, 'kkdu_misc_emflux_ci', &
                p2 = kkdu_misc_emflux_ci)
          CALL channel_halt(substr, status)
          CALL new_attribute(status, modstr, 'kkdu_misc_emflux_ci', &
                'long_name', c = 'Misc. dust emission flux, coarse mode')
          CALL channel_halt(substr, status)
          CALL new_attribute(status, modstr, 'kkdu_misc_emflux_ci', 'units', &
                c = 'kg m-2 s-1' )
          CALL channel_halt(substr, status)

    CASE('SO2_ant')
       ! output channel objects for SO2 tracer emissions
       CALL new_channel_object(status, modstr, &
            'SO2_emflux', p3=so2_emflux, reprid=GP_3D_MID )
       CALL channel_halt(substr, status)
       CALL new_attribute(status, modstr,  &
            'SO2_emflux', 'long_name', &
            c='SO2 emissions' )
       CALL channel_halt(substr, status)
       CALL new_attribute(status, modstr, &
            'SO2_emflux', 'units', c='molecules m-3 s-1' )
       CALL channel_halt(substr, status)

       ! input channel objects
       CALL new_channel_object(status, modstr, &
            'SO2_ant', p3=so2_ant, reprid=GP_3D_MID )
       CALL channel_halt(substr, status)
       CALL new_attribute(status, modstr,  &
            'SO2_ant', 'long_name', &
            c='SO2 emissions' )
       CALL channel_halt(substr, status)
       CALL new_attribute(status, modstr, &
            'SO2_ant', 'units', c='kg m-2 s-1' )
       CALL channel_halt(substr, status)

       CASE('BIOO')
          DO i=1, nolclass
             CALL new_channel_object(status, modstr, TRIM(olson_name(i)), &
                  p2=olson_emis(i)%ptr)
             CALL channel_halt(substr, status)
             CALL new_attribute(status, modstr, TRIM(olson_name(i)), &
                  'long_name', &
                  c='Bacterial emission flux('//TRIM(olson_name(i))//')' )
             CALL channel_halt(substr, status)
             CALL new_attribute(status, modstr, TRIM(olson_name(i)), &
                  'units', c='num. bacteria m-2 s-1' )
             CALL channel_halt(substr, status)
          END DO

       CASE('BIOM')
          ! MODIS emissions
          DO i=1, nmodisclass
             CALL new_channel_object(status, modstr, TRIM(modis_name(i)), &
                  p2=modis_emis(i)%ptr)
             CALL channel_halt(substr, status)
             CALL new_attribute(status, modstr, TRIM(modis_name(i)), &
                  'long_name', &
                  c='Bacterial emission flux('//TRIM(modis_name(i))//')' )
             CALL channel_halt(substr, status)
             CALL new_attribute(status, modstr, TRIM(modis_name(i)), &
                  'units', c='num. bacteria m-2 s-1' )
             CALL channel_halt(substr, status)
          END DO

!##########################################
!# LAI-based bioaerosol emissions         #
!##########################################
    CASE('BIOL')
          CALL new_channel_object(status, modstr, &
               'modis_lai', p2=modis_lai)
          CALL channel_halt(substr, status)
          CALL new_attribute(status, modstr, &
               'modis_lai', 'long_name', c='Leaf Area Index (monthly)' )
          CALL channel_halt(substr, status)
          CALL new_attribute(status, modstr, &
               'modis_lai', 'units', c='m2 m-2' )
          CALL channel_halt(substr, status)

          CALL new_channel_object(status, modstr, &
               'heald_emis', p2=heald_emis)
          CALL channel_halt(substr, status)
          CALL new_attribute(status, modstr, &
               'heald_emis', 'long_name',    &
               c='Fungal spore emissions (Heald and Spracklen)' )
          CALL channel_halt(substr, status)
          CALL new_attribute(status, modstr, &
               'heald_emis', 'units', c='m-2 s-1' )
          CALL channel_halt(substr, status)

          CALL new_channel_object(status, modstr, &
               'hummel_emis', p2=hummel_emis)
          CALL channel_halt(substr, status)
          CALL new_attribute(status, modstr, &
               'hummel_emis', 'long_name',    &
               c='Fungal spore emissions (Hummel)' )
          CALL channel_halt(substr, status)
          CALL new_attribute(status, modstr, &
               'hummel_emis', 'units', c='m-2 s-1' )
          CALL channel_halt(substr, status)

          CALL new_channel_object(status, modstr, &
               'js_emis', p2=js_emis)
          CALL channel_halt(substr, status)
          CALL new_attribute(status, modstr, &
               'js_emis', 'long_name',    &
               c='Pollen emission (Jacobson and Streets)' )
          CALL channel_halt(substr, status)
          CALL new_attribute(status, modstr, &
               'js_emis', 'units', c='m-2 s-1' )
          CALL channel_halt(substr, status)
!##########################################

       CASE ('terr13C')

          CALL new_channel_object(status, modstr, &
               'I12ISOP_emflux', p2=I12ISOP_emflux )
          CALL channel_halt(substr, status)
          CALL new_attribute(status, modstr, &
               'I12ISOP_emflux', 'long_name', c='12C isoprene emission flux' )
          CALL channel_halt(substr, status)
          CALL new_attribute(status, modstr, &
               'I12ISOP_emflux', 'units', c='molec. m-2 s-1' )
          CALL channel_halt(substr, status)

          CALL new_channel_object(status, modstr, &
               'I13ISOP_emflux', p2=I13ISOP_emflux )
          CALL channel_halt(substr, status)
          CALL new_attribute(status, modstr, &
               'I13ISOP_emflux', 'long_name', c='13C isoprene emission flux' )
          CALL channel_halt(substr, status)
          CALL new_attribute(status, modstr, &
               'I13ISOP_emflux', 'units', c='molec. m-2 s-1' )
          CALL channel_halt(substr, status)

       CASE('AirSnow')
          CALL new_channel_object(status, modstr, &
               'snow_air_flux_br2', p2=snow_air_flux_br2 )
          CALL channel_halt(substr, status)
          CALL new_attribute(status, modstr, &
               'snow_air_flux_br2', 'long_name', &
               c='Br2 over snow/ice emission' )
          CALL channel_halt(substr, status)
          CALL new_attribute(status, modstr, &
               'snow_air_flux_br2', 'units', c='molec m-2 s-1' )
          CALL channel_halt(substr, status)

          CALL new_channel_object(status, modstr, &
               'snow_air_flux_brcl', p2=snow_air_flux_brcl )
          CALL channel_halt(substr, status)
          CALL new_attribute(status, modstr, &
               'snow_air_flux_brcl', 'long_name', &
               c='BrCl over snow/ice emission' )
          CALL channel_halt(substr, status)
          CALL new_attribute(status, modstr, &
               'snow_air_flux_brcl', 'units', c='molec m-2 s-1' )
          CALL channel_halt(substr, status)

       !CASE()
          !##############################
          !### add new emissions here ###
          !##############################

       CASE DEFAULT
          CALL error_bi( 'unknown emission type',substr)
       END SELECT
    ENDDO emistype_loop1

    CALL end_message_bi(modstr, 'CHANNEL DEFINITION', substr)

CONTAINS

  SUBROUTINE make_repr_ndrydays

    IMPLICIT NONE

    CALL get_dimension_info(status,'ndrydays')
    IF (status == 0) RETURN

    CALL new_dimension(status, DIMID_ndrydays, &
         'numdrydays', ndrydays)
    CALL channel_halt(substr, status)
    ! NEW REPRESENTATION
    CALL new_representation(status, REPRID_ndrydays, 'REPRID_numdrydays' &
         , rank = 3, link = 'xxx-', dctype = DC_GP                       &
         , dimension_ids = (/ &
                _RI_XY_N_(DIMID_LON, DIMID_LAT,DIMID_ndrydays) /)        &
         , ldimlen       = (/ _RI_XY_N_(nproma, ngpblks, AUTO)   /)      &
         , output_order  = (/ _IX_XY_N_ , _IY_XY_N_ , _IN_XY_N_ /)   &
         , axis = repr_def_axes(_RI_XY_N_('X','Y','N'),'-')          &
         )
    CALL channel_halt(substr, status)
    nseg = gp_nseg
    ALLOCATE(start(nseg,IRANK))
    ALLOCATE(cnt(nseg,IRANK))
    ALLOCATE(meml(nseg,IRANK))
    ALLOCATE(memu(nseg,IRANK))

    start(:,1) = gp_start(:,1)
    cnt(:,1)   = gp_cnt(:,1)
    meml(:,1)  = gp_meml(:,1)
    memu(:,1)  = gp_memu(:,1)

    start(:,_IY_XY_N_) = gp_start(:,_IY_XYZN_)
    cnt(:,_IY_XY_N_)   = gp_cnt(:,_IY_XYZN_)
    meml(:,_IY_XY_N_)  = gp_meml(:,_IY_XYZN_)
    memu(:,_IY_XY_N_)  = gp_memu(:,_IY_XYZN_)

    start(:,_IN_XY_N_) = 1
    cnt(:,_IN_XY_N_)   = ndrydays
    meml(:,_IN_XY_N_)  = 1
    memu(:,_IN_XY_N_)  = ndrydays

    start(:,4) = 1
    cnt(:,4)   = 1
    meml(:,4)  = 1
    memu(:,4)  = 1

    CALL set_representation_decomp(status, REPRID_ndrydays &
         , start, cnt, memu, meml, .FALSE., PIOTYPE_COL)
    CALL channel_halt(substr, status)

    DEALLOCATE(start) ; NULLIFY(start)
    DEALLOCATE(cnt)   ; NULLIFY(cnt)
    DEALLOCATE(meml)  ; NULLIFY(meml)
    DEALLOCATE(memu)  ; NULLIFY(memu)

  END SUBROUTINE make_repr_ndrydays
  !------------------------------------------------------------------------

!------------------------------------------------------------------------

END SUBROUTINE onemis_init_memory
! --------------------------------------------------------------------------

  ! --------------------------------------------------------------------------
  SUBROUTINE onemis_init_coupling

    ! BML/MESSy
    USE messy_main_mpi_bi,           ONLY: p_parallel_io
    USE messy_main_tracer_mem_bi,    ONLY: GPTRSTR
    USE messy_main_tracer_tools_bi,  ONLY: tracer_halt
    USE messy_main_channel_error_bi, ONLY: channel_halt
    USE messy_main_data_bi,          ONLY: basemodstr=>modstr
    ! MESSy
    USE messy_main_tracer,           ONLY: get_tracer, tracer_error_str &
                                         , full2base_sub
    USE messy_main_channel,          ONLY: new_channel_object, new_attribute &
                                         , get_channel_object, get_channel_info

    IMPLICIT NONE

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER  :: substr = 'onemis_init_coupling'
    INTEGER                      :: status
    INTEGER                      :: jf, jt
    CHARACTER(LEN=STRLEN_MEDIUM) :: basename = ''
    CHARACTER(LEN=STRLEN_MEDIUM) :: subname  = ''
    CHARACTER(LEN=2*STRLEN_MEDIUM+1) :: fullname = ''
    INTEGER                      :: je
    INTEGER                      :: i

    ! AVOID EXPENSIVE SIMULATIONS WITH MISLEADING OUTPUT ...
    CALL get_channel_info(status, 'tagging')
    IF (status == 0) THEN
       fl: DO jf=1, NFLUXES
          gptl: DO jt=1, XF2T(jf)%ngpt
             IF ( (XF2T(jf)%mgp(jt) == 2) .AND.               &
                ( (TRIM(XF2T(jf)%tracer_gp(jt)) == 'NO') .OR. &
                  (TRIM(XF2T(jf)%tracer_gp(jt)) == 'C5H8')    &
                ) ) THEN
                CALL error_bi( &
                     'TAGGING GIVES WRONG RESULTS, IF TRACER '// &
                     &TRIM(XF2T(jf)%tracer_gp(jt))//' IS EMITTED WITH '//&
                     &'METHOD M=2. USE M=1 INSTEAD.', substr)
             END IF
          END DO gptl
       END DO fl
    ENDIF

    emistype_loop2: DO je = 1, max_emis

       IF (TRIM(EMIS_TYPE(je)) == '') CYCLE

       SELECT CASE(TRIM(EMIS_TYPE(je)))
       CASE('DMS')
          CALL get_channel_object(status, imp_seawater_dms%CHA &
               , imp_seawater_dms%OBJ, p2=seawater_dms )
          CALL channel_halt(substr, status)
       CASE('O3ice')
          ! nothing to do
       CASE ('SS_lsce')
          ! nothing to do
       CASE('SS_monahan')
          ! nothing to do

       CASE ('SS_aerocom')
          ! define input channel objects
          CALL get_channel_object(status, imp_numflx_as_aerocom%CHA  &
               , imp_numflx_as_aerocom%OBJ, p2=numflx_as_aerocom )
          CALL channel_halt(substr, status)

          CALL get_channel_object(status, imp_numflx_cs_aerocom%CHA &
               , imp_numflx_cs_aerocom%OBJ, p2=numflx_cs_aerocom )

          CALL get_channel_object(status, imp_massflx_as_aerocom%CHA &
               , imp_massflx_as_aerocom%OBJ, p2=massflx_as_aerocom )
          CALL channel_halt(substr, status)

          CALL get_channel_object(status, imp_massflx_cs_aerocom%CHA &
               , imp_massflx_cs_aerocom%OBJ, p2=massflx_cs_aerocom )
          CALL channel_halt(substr, status)

       CASE('SS_POC_AQUA')
         ! get OC ocean concentration
         CALL get_channel_object(status, imp_poc_aqua%cha  &
               , imp_poc_aqua%obj, p2=poc_aqua )
         CALL channel_halt(substr, status)
         ! get seasalt input flux(es)
         DO i=1,n_ss_poc_aqua
           CALL get_channel_object(status, emis_poc_aqua(i)%channel  &
               , emis_poc_aqua(i)%object, p2=emis_poc_aqua(i)%ss_flx )
           CALL channel_halt(substr, status)
         END DO
       CASE('SS_POC_SWIFS')
         ! get OC ocean concentration
         CALL get_channel_object(status, imp_poc_swifs%cha  &
               , imp_poc_swifs%obj, p2=poc_seawifs )
         CALL channel_halt(substr, status)
         ! get seasalt input flux(es)
         DO i=1,n_ss_poc_swifs
           CALL get_channel_object(status, emis_poc_swifs(i)%channel  &
               , emis_poc_swifs(i)%object, p2=emis_poc_swifs(i)%ss_flx )
           CALL channel_halt(substr, status)
         END DO
       CASE('SS_WIOC_AQUA')
         ! get OC ocean concentration
         CALL get_channel_object(status, imp_wioc_aqua%cha  &
               , imp_wioc_aqua%obj, p2=chlor_a_aqua_wioc )
         CALL channel_halt(substr, status)
         ! get seasalt input flux(es)
         DO i=1,n_ss_wioc_aqua
           CALL get_channel_object(status, emis_wioc_aqua(i)%channel  &
               , emis_wioc_aqua(i)%object, p2=emis_wioc_aqua(i)%ss_flx )
           CALL channel_halt(substr, status)
         END DO
       CASE('SS_WIOC_BLEN')
         ! get OC ocean concentration
         CALL get_channel_object(status, imp_wioc_blen%cha  &
               , imp_wioc_blen%obj, p2=chlor_a_blend_wioc )
         CALL channel_halt(substr, status)
         ! get seasalt input flux(es)
         DO i=1,n_ss_wioc_blend
           CALL get_channel_object(status, emis_wioc_blend(i)%channel  &
               , emis_wioc_blend(i)%object, p2=emis_wioc_blend(i)%ss_flx )
           CALL channel_halt(substr, status)
         END DO

       CASE('OC/BC')
          CALL get_channel_object(status, imp_BC_ag%CHA &
               , imp_BC_ag%OBJ, p2=BC_ag )
          CALL channel_halt(substr, status)

          CALL get_channel_object(status, imp_BC_ant%CHA &
               , imp_BC_ant%OBJ, p2=BC_ant )
          CALL channel_halt(substr, status)

          IF (TRIM(imp_BC_wf%OBJ) /= '') THEN
             CALL get_channel_object(status, imp_BC_wf%CHA  &
                  , imp_BC_wf%OBJ, p2=BC_wf )
             CALL channel_halt(substr, status)
          ELSE
             CALL new_channel_object(status, modstr, 'BC_wf', p2=BC_wf )
             CALL channel_halt(substr, status)
             CALL new_attribute(status, modstr, &
                  'BC_wf', 'long_name', c='BC wilfire emission dummy field' )
             CALL channel_halt(substr, status)
             CALL new_attribute(status, modstr, 'BC_wf', 'units', c='-' )
             CALL channel_halt(substr, status)
          ENDIF
          CALL get_channel_object(status, imp_OC_ag%CHA &
               , imp_OC_ag%OBJ, p2=OC_ag )
          CALL channel_halt(substr, status)

          CALL get_channel_object(status, imp_OC_ant%CHA &
               , imp_OC_ant%OBJ, p2=OC_ant )
          CALL channel_halt(substr, status)

          CALL get_channel_object(status, imp_OC_bge%CHA &
               , imp_OC_bge%OBJ, p2=OC_bge )
          CALL channel_halt(substr, status)

           IF (TRIM(imp_OC_wf%OBJ) /= '') THEN
              CALL get_channel_object(status, imp_OC_wf%CHA &
                   , imp_OC_wf%OBJ, p2=OC_wf )
              CALL channel_halt(substr, status)
          ELSE
             CALL new_channel_object(status, modstr, 'OC_wf', p2=OC_wf )
             CALL channel_halt(substr, status)
             CALL new_attribute(status, modstr, &
                  'OC_wf', 'long_name', c='OC wilfire emission dummy field' )
             CALL channel_halt(substr, status)
             CALL new_attribute(status, modstr, 'OC_wf', 'units', c='-' )
             CALL channel_halt(substr, status)
          ENDIF

       CASE ('CH4')
          CALL get_channel_object(status, imp_ch4_conc_clim%CHA, &
               imp_ch4_conc_clim%OBJ, p2=ch4_conc_clim )
          CALL channel_halt(substr, status)

       CASE ('VOC')
          CALL get_channel_object(status, imp_ISOP_emisfac%CHA , &
               imp_ISOP_emisfac%OBJ, p2=ISOP_emisfac )
          CALL channel_halt(substr, status)

          CALL get_channel_object(status, imp_MTERP_emisfac%CHA , &
               imp_MTERP_emisfac%OBJ, p2=MTERP_emisfac )
          CALL channel_halt(substr, status)

          CALL get_channel_object(status, imp_drymatter%CHA  , &
               imp_drymatter%OBJ, p2=drymatter )
          CALL channel_halt(substr, status)

          IF (.NOT. ASSOCIATED(lai)) THEN
             CALL get_channel_object(status, imp_lai%CHA &
                  , imp_lai%OBJ, p2=lai )
             CALL channel_halt(substr, status)
          END IF

          CALL get_channel_object(status, imp_hc%CHA, &
               imp_hc%OBJ, p2=hc )
          CALL channel_halt(substr, status)

          CALL get_channel_object(status, imp_drag%CHA &
               , imp_drag%OBJ, p2=drag )
          CALL channel_halt(substr, status)

          CALL get_channel_object(status, imp_disp%CHA &
               , imp_disp%OBJ, p2=disp )
          CALL channel_halt(substr, status)

          CALL get_channel_object(status, imp_forestfr%CHA &
               , imp_forestfr%OBJ, p2=forestfr )
          CALL channel_halt(substr, status)

          CALL get_channel_object(status, imp_lad_top%CHA &
               , imp_lad_top%OBJ, p2=lad_top)
          CALL channel_halt(substr, status)

          CALL get_channel_object(status, imp_lad_lay2%CHA &
               , imp_lad_lay2%OBJ, p2=lad_lay2)
          CALL channel_halt(substr, status)

          CALL get_channel_object(status, imp_lad_lay3%CHA &
               , imp_lad_lay3%OBJ, p2=lad_lay3)
          CALL channel_halt(substr, status)

          CALL get_channel_object(status, imp_lad_soil%CHA &
               , imp_lad_soil%OBJ, p2=lad_soil)
          CALL channel_halt(substr, status)

          l_cossza = .TRUE.

       CASE('NO')

          CALL get_channel_object(status, imp_NOemisclass1%CHA &
               , imp_NOemisclass1%OBJ, p3=NOemisclass1)
          CALL channel_halt(substr, status)

          CALL get_channel_object(status, imp_NOemisclass2%CHA &
               , imp_NOemisclass2%OBJ, p3=NOemisclass2)
          CALL channel_halt(substr, status)

          CALL get_channel_object(status, imp_cultiv%CHA &
               , imp_cultiv%OBJ, p2=cultiv )
          CALL channel_halt(substr, status)

          CALL get_channel_object(status, imp_fertil%CHA  &
               , imp_fertil%OBJ, p2=fertil )
          CALL channel_halt(substr, status)

          IF (.NOT. ASSOCIATED(lai)) THEN
             CALL get_channel_object(status, imp_lai%CHA  &
                  , imp_lai%OBJ , p2=lai)
             CALL channel_halt(substr, status)
          END IF

          IF (SIZE(noemisclass1,_IN_XY_N_) /= ncl_yl95) &
               CALL error_bi( &
               'number of NO emission classes (1) do not'&
               &//' match emission algorithm !',substr)
          IF (SIZE(noemisclass2,_IN_XY_N_) /= ncl_yl95) &
               CALL error_bi( &
               'number of NO emission classes (2) do not'&
               &//' match emission algorithm !',substr)

       CASE('NOpls')
          ! nothing to do

       CASE('NO_yl95sl10')
           CALL get_channel_object(status, imp_lai_yl95sl10%CHA          &
                ,imp_lai_yl95sl10%OBJ , p2=lai_yl95sl10)
          CALL channel_halt(substr, status)

          CALL get_channel_object(status,imp_fertil_yl95sl10%CHA         &
               , imp_fertil_yl95sl10%OBJ, p2=fertil_yl95sl10)
          CALL channel_halt(substr, status)

          CALL get_channel_object(status, imp_NOemclass_yl95sl10%CHA     &
               , imp_NOemclass_yl95sl10%OBJ, p3=NOemclass_yl95sl10)
          CALL channel_halt(substr, status)

          CALL get_channel_object(status, imp_rootdepth_yl95sl10%CHA     &
               , imp_rootdepth_yl95sl10%OBJ, p2=root_depth)
          CALL channel_halt(substr, status)

          CALL get_channel_object(status, imp_rootdepth_mask_yl95sl10%CHA&
               , imp_rootdepth_mask_yl95sl10%OBJ, p2=root_depth_mask)
          CALL channel_halt(substr, status)

         IF (SIZE(NOemclass_yl95sl10,_IN_XY_N_) /= ncl_yl95sl10) &
               CALL error_bi( &
               'number of NO emission classes (1) do not'&
               &//' match emission algorithm !',substr)

       CASE('DU')

          CALL get_channel_object(status, imp_du_thr%CHA &
               , imp_du_thr%OBJ, p2=du_thr )
          CALL channel_halt(substr, status)

          CALL get_channel_object(status, imp_du_src%CHA &
               , imp_du_src%OBJ, p2=du_src )
          CALL channel_halt(substr, status)

          CALL get_channel_object(status, imp_du_cla%CHA &
               , imp_du_cla%OBJ, p2=du_cla )
          CALL channel_halt(substr, status)

       CASE('DU_tegen')
          CALL get_channel_object(status, imp_mat_s2%CHA &
               , imp_mat_s2%OBJ, p2=mat_s2 )
          CALL channel_halt(substr, status)

          CALL get_channel_object(status, imp_mat_s3%CHA &
               , imp_mat_s3%OBJ, p2=mat_s3 )
          CALL channel_halt(substr, status)

          CALL get_channel_object(status, imp_mat_s4%CHA &
               , imp_mat_s4%OBJ, p2=mat_s4 )
          CALL channel_halt(substr, status)

          CALL get_channel_object(status, imp_mat_s6%CHA &
               , imp_mat_s6%OBJ, p2=mat_s6 )
          CALL channel_halt(substr, status)

          CALL get_channel_object(status, imp_mat_psrc%CHA &
               , imp_mat_psrc%OBJ, p2=mat_psrc )
          CALL channel_halt(substr, status)

          CALL get_channel_object(status, imp_k_fpar_eff%CHA &
               , imp_k_fpar_eff%OBJ, p2=k_fpar_eff )
          CALL channel_halt(substr, status)

          ! INITIALISE DUST CONSTANTS

          ! Calculation of scale factor for wind stress threshold to adjust
          ! dust emission occurence (Tegen et al., 2004) -> cuscale
          IF(cuscale_in == -1._dp) THEN
#ifdef ECHAM5
             ! nlon(T42) = 128, nlon(T159) = 480
             IF(nlon < 128) THEN ! i.e. HRES < T42
                cuscale = 0.6_dp
             ELSEIF (nlon > 480) THEN ! i.e. HRES > T159
                cuscale = 1._dp
             ELSE ! i.e. T42 <= HRES <= T159
                cuscale = 0.0894_dp * (360._dp/REAL(nlon))**2._dp &
                     - 0.4787_dp * (360._dp/REAL(nlon)) + 1.2962_dp
             ENDIF
             IF (.NOT. l_nudging) cuscale = cuscale + 0.06_dp
#endif
#ifdef COSMO
             IF(dlon /= dlat) dlon = MIN(dlon,dlat)
             cuscale = -0.142_dp * LOG(dlon) + 0.4872_dp
#endif
             cuscale = MIN(1._dp,cuscale)
             cuscale = MAX(0.6_dp,cuscale)
             IF (p_parallel_io) THEN
                WRITE(*,*) 'Correction factor for Tegen dust emission scheme:'
                WRITE(*,'(a10,f9.4)') ' cuscale =', cuscale
             END IF
          ELSEIF(cuscale_in >= 0.6_dp .AND. cuscale_in <= 1._dp) THEN
             cuscale = cuscale_in
             IF (p_parallel_io) THEN
                WRITE(*,*) 'Correction factor for Tegen dust emission scheme (from nml):'
                WRITE(*,'(a10,f9.4)') ' cuscale =', cuscale
             END IF
          ELSE
#ifdef ECHAM5
             ! nlon(T42) = 128, nlon(T159) = 480
             IF(nlon < 128) THEN ! i.e. HRES < T42
                cuscale = 0.6_dp
             ELSEIF (nlon > 480) THEN ! i.e. HRES > T159
                cuscale = 1._dp
             ELSE ! i.e. T42 <= HRES <= T159
                cuscale = 0.0894_dp * (360._dp/REAL(nlon))**2._dp &
                     - 0.4787_dp * (360._dp/REAL(nlon)) + 1.2962_dp
             ENDIF
             IF (.NOT. l_nudging) cuscale = cuscale + 0.06_dp
#endif
#ifdef COSMO
             IF(dlon /= dlat) dlon = MIN(dlon,dlat)
             cuscale = -0.142_dp * LOG(dlon) + 0.4872_dp
#endif
             cuscale = MIN(1._dp,cuscale)
             cuscale = MAX(0.6_dp,cuscale)
             IF (p_parallel_io) THEN
                WRITE(*,*) 'WARNING: Correction factor for Tegen dust emission scheme!'
                WRITE(*,*) 'given value for cuscale outside reasonalbe range = [0.6,1.0]'
                WRITE(*,*) 'Correction factor for Tegen dust emission scheme recalculated:'
                WRITE(*,'(a10,f9.4)') ' cuscale =', cuscale
             END IF
          ENDIF
          CALL dust_tegen_init(cuscale)

       CASE('DU_Astitha1')
          CALL get_channel_object(status, imp_du_cla2%CHA &
               , imp_du_cla2%OBJ, p2=du_cla2 )
          CALL channel_halt(substr, status)

          CALL get_channel_object(status, imp_rdepth%CHA &
               , imp_rdepth%OBJ, p2=rdepth )
          CALL channel_halt(substr, status)

          CALL get_channel_object(status, imp_dustsrc%CHA &
               , imp_dustsrc%OBJ, p4=dustsrc )
          CALL channel_halt(substr, status)

          IF (SIZE(dustsrc,_IN_XY_N_) /= nbiomes) &
               CALL error_bi('number of dustsrc emission do not'&
               &//' match emission algorithm !',substr)

          CALL get_channel_object(status, imp_lai_in%CHA &
               , imp_lai_in%OBJ, p2=lai_in )
          CALL channel_halt(substr, status)
          if (l_ducomp) then
             CALL get_channel_object(status &
                  , imp_du_nap%CHA , imp_du_nap%OBJ, p2=du_nap)
             CALL channel_halt(substr, status)
             CALL get_channel_object(status&
                  , imp_du_kp%CHA , imp_du_kp%OBJ, p2=du_kp)
             CALL channel_halt(substr, status)
             CALL get_channel_object(status&
                  , imp_du_capp%CHA , imp_du_capp%OBJ, p2=du_capp)
             CALL channel_halt(substr, status)
             CALL get_channel_object(status&
                  , imp_du_mgpp%CHA , imp_du_mgpp%OBJ, p2=du_mgpp)
             CALL channel_halt(substr, status)
             CALL get_channel_object(status&
                  , imp_du_misc%CHA , imp_du_misc%OBJ, p2=du_misc)
             CALL channel_halt(substr, status)
          endif

        CASE('DU_Astitha2')
           CALL get_channel_object(status, imp_du_cla2%CHA &
               , imp_du_cla2%OBJ, p2=du_cla2 )
          CALL channel_halt(substr, status)

          CALL get_channel_object(status, imp_rdepth%CHA &
               , imp_rdepth%OBJ, p2=rdepth )
          CALL channel_halt(substr, status)

          CALL get_channel_object(status, imp_dustsrc%CHA &
               , imp_dustsrc%OBJ, p4=dustsrc )
          CALL channel_halt(substr, status)

          IF (SIZE(dustsrc,_IN_XY_N_) /= nbiomes) &
               CALL error_bi('number of dustsrc emission do not'&
               &//' match emission algorithm !',substr)

          CALL get_channel_object(status, imp_lai_in%CHA &
               , imp_lai_in%OBJ, p2=lai_in )
          CALL channel_halt(substr, status)

          CALL get_channel_object(status, imp_soiltexture%CHA &
               , imp_soiltexture%OBJ, p4=soiltexture )
          CALL channel_halt(substr, status)

          IF (SIZE(soiltexture,_IN_XY_N_) /= nclasses) &
               CALL error_bi( &
               'number of soiltexture do not'&
               &//' match emission algorithm !',substr)

          if (l_ducomp) then
             CALL get_channel_object(status &
                  , imp_du_nap%CHA , imp_du_nap%OBJ, p2=du_nap)
             CALL channel_halt(substr, status)
             CALL get_channel_object(status &
                  , imp_du_kp%CHA , imp_du_kp%OBJ, p2=du_kp)
             CALL channel_halt(substr, status)
             CALL get_channel_object(status &
                  , imp_du_capp%CHA , imp_du_capp%OBJ, p2=du_capp)
             CALL channel_halt(substr, status)
             CALL get_channel_object(status &
                  , imp_du_mgpp%CHA , imp_du_mgpp%OBJ, p2=du_mgpp)
             CALL channel_halt(substr, status)
             CALL get_channel_object(status &
                  , imp_du_misc%CHA , imp_du_misc%OBJ, p2=du_misc)
             CALL channel_halt(substr, status)
          endif

       CASE('KKDU_Astitha')
          CALL get_channel_object(status, imp_kkdu_clay%CHA, &
                imp_kkdu_clay%OBJ, p2 = kkdu_clay)
          CALL channel_halt(substr, status)
          CALL get_channel_object(status, imp_kkdu_mask%CHA, &
                imp_kkdu_mask%OBJ, p3 = kkdu_mask)
          CALL channel_halt(substr, status)
          CALL get_channel_object(status, imp_kkdu_lai%CHA, &
                imp_kkdu_lai%OBJ, p2 = kkdu_lai)
          CALL channel_halt(substr, status)
          CALL get_channel_object(status, imp_kkdu_topo%CHA, &
                imp_kkdu_topo%OBJ, p2 = kkdu_topo)
          CALL channel_halt(substr, status)
          CALL get_channel_object(status, imp_kkdu_nap%CHA, imp_kkdu_nap%OBJ, &
                p2 = kkdu_nap)
          CALL channel_halt(substr, status)
          CALL get_channel_object(status, imp_kkdu_kp%CHA, imp_kkdu_kp%OBJ, &
                p2 = kkdu_kp)
          CALL channel_halt(substr, status)
          CALL get_channel_object(status, imp_kkdu_capp%CHA, &
                imp_kkdu_capp%OBJ, p2 = kkdu_capp)
          CALL channel_halt(substr, status)
          CALL get_channel_object(status, imp_kkdu_mgpp%CHA, &
                imp_kkdu_mgpp%OBJ, p2 = kkdu_mgpp)
          CALL channel_halt(substr, status)
          CALL get_channel_object(status, imp_kkdu_misc%CHA, &
                imp_kkdu_misc%OBJ, p2 = kkdu_misc)
          CALL channel_halt(substr, status)

    CASE('SO2_ant')
       CALL get_channel_object(status, imp_SO2_ant_high%CHA  &
            , imp_SO2_ant_high%OBJ, p2=so2_ant_high)
       CALL channel_halt(substr, status)
       CALL get_channel_object(status, imp_SO2_ant_low%CHA  &
            , imp_SO2_ant_low%OBJ, p2=so2_ant_low)
       CALL channel_halt(substr, status)

       CASE('BIOO')
          ! OLSON LUMPED EMISSIONS
          CALL get_channel_object(status, imp_olson%CHA &
               , imp_olson%OBJ, p3=olson)

          IF (SIZE(olson,_IN_XY_N_) /= nolclass) &
               CALL error_bi( &
               'Number of lumped ecosystem emission classes does not'&
               &//' match emissions algorithm !',substr)

       CASE('BIOM')
          ! MODIS emissions

          CALL get_channel_object(status, imp_modis%CHA &
               , imp_modis%OBJ,p3=modis)
          CALL channel_halt(substr, status)

          IF (SIZE(modis,_IN_XY_N_) /= nmodisclass) &
               CALL error_bi( &
               'Number of ecosystem emission classes does not'&
               &//' match emissions algorithm !',substr)

       CASE ('terr13C')

          CALL get_channel_object(status, imp_terr13C_delta%CHA  &
               , imp_terr13C_delta%OBJ, p2=terr13C_delta )
          CALL channel_halt(substr, status)

          l_cossza = .TRUE.

       CASE ('AirSnow')
          ! Get ddepflux for HOBr, BrNO3 (BrONO2), HBr, and O3
          CALL get_channel_object(status, imp_ddepflux_HOBr%CHA &
               , imp_ddepflux_HOBr%OBJ, p2=ddepflux_HOBr )
          CALL channel_halt(substr, status)

          CALL get_channel_object(status, imp_ddepflux_BrNO3%CHA &
               , imp_ddepflux_BrNO3%OBJ, p2=ddepflux_BrNO3 )
          CALL channel_halt(substr, status)

          CALL get_channel_object(status, imp_ddepflux_HBr%CHA &
               , imp_ddepflux_HBr%OBJ, p2=ddepflux_HBr )
          CALL channel_halt(substr, status)

          CALL get_channel_object(status, imp_ddepflux_O3%CHA &
               , imp_ddepflux_O3%OBJ, p2=ddepflux_O3 )
          CALL channel_halt(substr, status)

          CALL get_channel_object(status, imp_sic_multi_year%CHA &
               , imp_sic_multi_year%OBJ, p2=sic_multi_year )
          CALL channel_halt(substr, status)

      !CASE()
          !##############################
          !### add new emissions here ###
          !##############################

       CASE DEFAULT
          CALL error_bi( 'unknown emission type',substr)
       END SELECT
    ENDDO emistype_loop2

    IF (l_cossza) THEN
       CALL get_channel_object(status, imp_cossza%CHA , &
            imp_cossza%OBJ, p2=cossza_2d )
       CALL channel_halt(substr, status)
    END IF

    CALL get_channel_object(status, TRIM(basemodstr), 'wind10', p2=wind10_2d)
    CALL channel_halt(substr, status)

    CALL end_message_bi(modstr, 'END CHANNEL COUPLING', substr)


    CALL start_message_bi(modstr, 'COUPLING INITIALIZATION', substr)

    ! SET POINTER TO CHANNEL OBJECTS
    IF (p_parallel_io) THEN
       WRITE(*,*) 'LOOKING IN CHANNEL ... ''',TRIM(modstr),''''
    END IF

    fluxes_loop: DO jf=1, NFLUXES
       ! ASSOCIATE F2T TO CHANNEL OBJECT POINTERS
       IF (p_parallel_io) THEN
          WRITE(*,*) '  > ',TRIM(XF2T(jf)%name)
       END IF
       CALL get_channel_object(status, modstr &
            , TRIM(XF2T(jf)%name), p3=XF2T(jf)%ptr )

       IF (status /= 0) THEN
          CALL warning_bi('WARNING: channel object '''//&
               &TRIM(XF2T(jf)%name)//'''not found !',substr)
          CYCLE  ! SKIP TRACER ASSOCIATION, IF CHANNEL OBJECT NOT PRESENT
       END IF
       ! mz_pj_20060328+

       IF (SIZE(XF2T(jf)%ptr,_IN_XY_N_) > 1) THEN
          DO jt=1, XF2T(jf)%ngpt
             IF (XF2T(jf)%mgp(jt) > 1) THEN
                IF (p_parallel_io) THEN
                   WRITE(*,*) 'WARNING: GP METHOD FOR 3D RESET !'
                END IF
                XF2T(jf)%mgp(jt) = 1
             END IF
          END DO
       END IF
!!#D attila +
       IF (L_LG) THEN
          IF (SIZE(XF2T(jf)%ptr,_IN_XY_N_) > 1) THEN
             DO jt=1, XF2T(jf)%nlgt
                IF (XF2T(jf)%mlg(jt) > 1) THEN
                   IF (p_parallel_io) THEN
                      WRITE(*,*) 'WARNING: LG METHOD FOR 3D RESET !'
                   END IF
                   XF2T(jf)%mlg(jt) = 1
                END IF
             END DO
          END IF
       END IF
!!#D attila -

       IF (p_parallel_io) THEN
          WRITE(*,*) '  ... ASSOCIATED GRIDPOINT TRACERS ...'
       END IF
       gp_tracer_loop: DO jt=1, XF2T(jf)%ngpt

          fullname = TRIM(XF2T(jf)%tracer_gp(jt))
          IF (p_parallel_io) THEN
             WRITE(*,*) '  ... ... ',TRIM(fullname)
          END IF

          CALL full2base_sub(status, TRIM(fullname) &
               , basename, subname)
          CALL tracer_halt(substr, status)
          CALL get_tracer(status, GPTRSTR, basename &
               , subname=subname, idx=XF2T(jf)%idt_gp(jt))

          IF (status /= 0) THEN
               CALL error_bi('tracer '''&
               &//TRIM(fullname)//''' not found: '&
               &//tracer_error_str(status),substr)

              XF2T(jf)%idt_gp(jt) = 0
            ENDIF

          IF ((XF2T(jf)%mgp(jt) < 0) .OR. (XF2T(jf)%mgp(jt) > 2)) THEN
             CALL error_bi( 'unknown gp-emission method !',substr)
          ELSE
             IF (p_parallel_io) THEN
                WRITE(*,*) '  ... ... ... method  ',XF2T(jf)%mgp(jt)
                WRITE(*,*) '  ... ... ... scaling ',XF2T(jf)%efact_gp(jt)
             END IF
          END IF

       END DO gp_tracer_loop

       ! SPECIAL CASE: CH4-EMISSION MUST BE ASSOCIATED TO EXACTLY 1 TRACER
       IF ( (TRIM(XF2T(jf)%name)) == 'CH4_emflux') THEN
          IF (XF2T(jf)%ngpt /= 1) &
               CALL error_bi( &
               'CH4_emflux MUST BE ASSOCIATED TO 1 TRACER',substr)
          IF (XF2T(jf)%idt_gp(1) == 0) &
               CALL error_bi( &
               'CH4_emflux MUST BE ASSOCIATED TO EXISTING TRACER',substr)
          idt_CH4 = XF2T(jf)%idt_gp(1)
       END IF

    END DO fluxes_loop

    CALL end_message_bi(modstr, 'COUPLING INITIALIZATION', substr)

  END SUBROUTINE onemis_init_coupling

  ! --------------------------------------------------------------------------

  SUBROUTINE onemis_global_start

    ! BML/MESSy
    USE messy_main_grid_def_mem_bi,      ONLY: nlev

    IMPLICIT NONE

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'onemis_global_start'

    INTEGER :: je, i

    emistype_loop3: DO je=1, max_emis

       ! INCLUDE EMISSION SPECIFIC CALCULATIONS
       SELECT CASE(TRIM(EMIS_TYPE(je)))
       CASE('VOC')
          z0m(:,:)=100./EXP(1./SQRT(MAX(0.001_dp,drag(:,:))))
          lad(_RI_XY_N_(:,:,1))=lad_top(:,:)
          lad(_RI_XY_N_(:,:,2))=lad_lay3(:,:)
          lad(_RI_XY_N_(:,:,3))=lad_lay2(:,:)
          lad(_RI_XY_N_(:,:,4))=lad_soil(:,:)
       CASE('NO')

          NOemis_w(:,:) = 0._dp
          NOemis_d(:,:) = 0._dp

          trop_class = 11
          DO i=1, ncl_yl95
             IF (i == trop_class)  CYCLE
             NOemis_w(:,:) = NOemis_w(:,:) + REAL(noemfact_wet(i),dp)/2.   &
                    * (noemisclass1(_RI_XY_N_(:,:,i)) + &
                          noemisclass2(_RI_XY_N_(:,:,i)))
             NOemis_d(:,:) = NOemis_d(:,:) + REAL(noemfact_dry(i),dp)/2.   &
                  * (noemisclass1(_RI_XY_N_(:,:,i)) +   &
                         noemisclass2(_RI_XY_N_(:,:,i)))
          ENDDO

       CASE('NO_yl95sl10')
          root_depth = root_depth * root_depth_mask
       CASE('SO2_ant')
          so2_ant(_RI_X_ZN_(:,nlev-1,:)) =  so2_ant_high(:,:)
          so2_ant(_RI_X_ZN_(:,nlev,:))   =  so2_ant_low(:,:)
       CASE DEFAULT
          ! NOTHING TO DO HERE
       END SELECT
    ENDDO emistype_loop3

  END SUBROUTINE onemis_global_start

  ! --------------------------------------------------------------------------

  SUBROUTINE onemis_vdiff

    ! Notes from PJ:
    ! * the unit MUST BE molecules m-2 s-1
    ! * the resulting tracer tendency for method 1 (lowest grid layer) is
    !   mol/mol/s
    ! * for method 2 (lower boundary condition of vertical flux),
    !   pxtems in in vdiff multiplied with ZTMST*G*ZQDP (s m s-2 Pa-1),
    !   therefore pxtems MUST BE in
    !   (mol(tracer)/mol(air)) * kg (air) m-2 s-1

    ! ECHAM5/MESSy
#ifndef MESSYTENDENCY
    USE messy_main_tracer_mem_bi,   ONLY: pxtte => qxtte, pxtm1 => qxtm1
#endif
    USE messy_main_grid_def_mem_bi, ONLY: nlev, jrow, kproma, ngpblks, nproma
    USE messy_main_grid_def_bi,     ONLY: philat_2d, philon_2d, gboxarea_2d &
                                        , deltaz
    USE messy_main_data_bi,         ONLY: tsw, rho_air_dry_3d              &
#if defined(ECHAM5) || defined(CESM1)
                                       , pxtems                            &
#endif
                                       , slm, seaice, slf, cvs             &
                                       , srfl,tslm1, prc, prl              &
                                       , tsurf_2d                          &
#if defined(VERTICO)
                                       , ws => wsoil                       &
#else
                                       , ws                                &
#endif
                                       , tsoil, alake                      &
                                       , glac, wsmx                        &
                                       , zust_2d, press_3d, tm1_3d, qm1_3d
    USE messy_main_timer,          ONLY: time_step_len, delta_time   &
                                       , nstep=>current_time_step    &
                                       , imonth=>MONTH               &
                                       , INIT_STEP
    USE messy_main_constants_mem,  ONLY: g, M_air, N_A, OneDay

    IMPLICIT NONE

    ! LOCAL
    !CHARACTER(LEN=*), PARAMETER         :: substr = 'onemis_vdiff'
    REAL(DP), PARAMETER                 :: uconv = N_A * 1.E3_dp/M_air

    REAL(DP), DIMENSION(:,:), POINTER     :: zairdens => NULL() ! air density
    REAL(DP), DIMENSION(:,:), ALLOCATABLE :: zdz       ! layer thickness
    REAL(DP), DIMENSION(:),   ALLOCATABLE :: zscale    !
#ifdef MESSYTENDENCY
    REAL(DP), DIMENSION(:,:), ALLOCATABLE :: zxtte     ! local tracer tendency
    REAL(dp) :: ch4start(nproma,nlev)
#endif
    INTEGER                             :: jt, jf, je, jl, ni
    INTEGER                             :: idt !um_ak_20120625 , idx, idx2
    INTEGER                             :: jp ! op_pj_20140129
    !
    REAL(DP) :: voc_emisfac(nproma,3,ngpblks)
    REAL(DP) :: voc_emflux(nproma,3,ngpblks)

    ! definition for CH4 emissions:
    REAL(dp)   :: concch4(nproma)

     ! definitions for carbon emissions:
    REAL(dp) :: zCseason(nproma)
    REAL(dp), PARAMETER :: zCseason_nh(1:12)  =  &
         (/1.146_dp, 1.139_dp, 1.081_dp, 0.995_dp, 0.916_dp, 0.920_dp, &
         0.910_dp, 0.907_dp, 0.934_dp, 0.962_dp, 1.019_dp, 1.072_dp  /) ! INIT

    REAL(dp) :: prect(nproma)

    ! CALCULATE AIR DENSITY AND LAYER THICKNESS
    ALLOCATE(zdz(kproma, nlev))
    ALLOCATE(zscale(kproma))
#ifdef MESSYTENDENCY
    ALLOCATE(zxtte(kproma,nlev))
#endif
    zairdens => rho_air_dry_3d(_RI_XYZ__(1:kproma,jrow,1:nlev))

    zdz(:,1:nlev) = deltaz(_RI_XYZ__(1:kproma,jrow,1:nlev))
    zscale(1:kproma) = zairdens(1:kproma,nlev)*zdz(1:kproma,nlev)

    emistype_loop4: DO je=1, max_emis

       ! STEP 1: CALCULATE EMISSION FLUXES

       SELECT CASE(TRIM(EMIS_TYPE(je)))
       CASE('DMS')
          CALL dms_emissions(               &
               emis_dms_sea(1:kproma,jrow), &
               seawater_dms(1:kproma,jrow), &
               wind10_2d(1:kproma,jrow),    &
               tsw(1:kproma,jrow),          &
               slm(1:kproma,jrow),          &
               seaice(1:kproma,jrow)        )
       CASE('SS_lsce')
          CALL seasalt_emissions_lsce(     wind10_2d(1:kproma,jrow),   &
               slf(1:kproma,jrow),         seaice(1:kproma,jrow),      &
               alake(1:kproma,jrow),                                   &
               mss_as_lsce(1:kproma,jrow), mss_cs_lsce(1:kproma,jrow), &
               nss_as_lsce(1:kproma,jrow), nss_cs_lsce(1:kproma,jrow))

          emisflx_ss_cs_sum(1:kproma,jrow) =  emisflx_ss_cs_sum(1:kproma,jrow) &
               + mss_cs_lsce(1:kproma,jrow) * M_air / 58.44 * delta_time
          emisflx_ss_as_sum(1:kproma,jrow) =  emisflx_ss_as_sum(1:kproma,jrow) &
               + mss_as_lsce(1:kproma,jrow) * M_air / 58.44 * delta_time


       CASE('SS_monahan')
          CALL seasalt_emissions_monahan(     wind10_2d(1:kproma,jrow),      &
               slf(1:kproma,jrow),            seaice(1:kproma,jrow),         &
               alake(1:kproma,jrow),                                         &
               mss_as_monahan(1:kproma,jrow), mss_cs_monahan(1:kproma,jrow), &
               nss_as_monahan(1:kproma,jrow), nss_cs_monahan(1:kproma,jrow))

       CASE('SS_aerocom')
          CALL seasalt_emissions_aerocom( &
               mss_as_aerocom(1:kproma,jrow), mss_cs_aerocom(1:kproma,jrow), &
               nss_as_aerocom(1:kproma,jrow), nss_cs_aerocom(1:kproma,jrow), &
               numflx_as_aerocom(1:kproma,jrow),                             &
               numflx_cs_aerocom(1:kproma,jrow),                             &
               massflx_as_aerocom(1:kproma,jrow),                            &
               massflx_cs_aerocom(1:kproma,jrow),                            &
               gboxarea_2d(1:kproma,jrow)   )

       CASE('OC/BC')
          zCseason(1:kproma) = 1._dp
          WHERE (philat_2d(1:kproma,jrow) > 0._dp)
             zCseason(1:kproma) = zCseason_nh(imonth)
          ENDWHERE
          CALL carbon_emissions(kproma,     OC_ag(1:kproma,jrow),           &
               OC_ant(1:kproma,jrow),       OC_bge(1:kproma,jrow),          &
               OC_wf(1:kproma,jrow),        BC_ag(1:kproma,jrow),           &
               BC_ant(1:kproma,jrow),       BC_wf(1:kproma,jrow),           &
               BC_sum_insol(1:kproma,jrow), OC_sum_insol(1:kproma,jrow),    &
               OC_sum_sol(1:kproma,jrow),   N_sol(1:kproma,jrow),           &
               N_insol(1:kproma,jrow),      zCseason(1:kproma),             &
               N_insol_bc(1:kproma,jrow),   N_insol_oc(1:kproma,jrow),      &
               OC_soa_sol(1:kproma,jrow),   OC_bb_sol(1:kproma,jrow),       &
               OC_ff_sol(1:kproma,jrow),    OC_ff_insol(1:kproma,jrow),     &
               OC_soa_insol(1:kproma,jrow), OC_bb_insol(1:kproma,jrow),     &
               BC_ff_insol(1:kproma,jrow),  BC_bb_insol(1:kproma,jrow) )
          CASE('O3ice')
             CALL O3ice_emissions(kproma, slm(1:kproma,jrow)      &
                                        , seaice(1:kproma,jrow)   &
                                        , cvs(1:kproma,jrow)      &
                                        , O3_emflux(1:kproma,jrow))

          CASE('CH4')
             idt = idt_CH4
#ifndef MESSYTENDENCY
             concch4(1:kproma) =  pxtm1(_RI_X_ZN_(1:kproma,nlev,idt))    &
                  + pxtte(_RI_X_ZN_(1:kproma,nlev,idt)) * time_step_len
#else
          CALL mtend_get_start_l(idt, v0=ch4start)
          concch4(1:kproma) = ch4start(1:kproma,nlev)
#endif
             CALL CH4_emissions(kproma,  zairdens(1:kproma,nlev)          &
                  , zdz(1:kproma,nlev),  time_step_len                    &
                  , concch4(1:kproma),   CH4_conc_clim(1:kproma,jrow) &
                  , CH4_emflux(1:kproma,jrow))
          CASE('VOC')
             voc_emisfac(:,:,jrow) = 0.0_DP
             voc_emflux(:,:,jrow) = 0.0_DP
             voc_emisfac(1:kproma,1,jrow) = isop_emisfac(1:kproma,jrow)
             voc_emisfac(1:kproma,2,jrow) = mterp_emisfac(1:kproma,jrow)
            CALL VOC_emissions(kproma,          srfl(1:kproma,jrow)        &
                  ,cossza_2d(1:kproma,jrow),    lai(1:kproma,jrow)         &
                  ,lad(_RI_XYZ__(1:kproma,jrow,:)),      drymatter(1:kproma,jrow)   &
                  ,voc_emisfac(1:kproma,:,jrow),tslm1(1:kproma,jrow)       &
                  ,voc_emflux(1:kproma,:,jrow))
             isop_emflux(1:kproma,jrow)  = voc_emflux(1:kproma,1,jrow)
             mterp_emflux(1:kproma,jrow) = voc_emflux(1:kproma,2,jrow)

          CASE('terr13C')
             ! calculating regular isoprene flux
             voc_emisfac(:,:,jrow) = 0.0_DP
             voc_emflux(:,:,jrow) = 0.0_DP
             voc_emisfac(1:kproma,1,jrow) = isop_emisfac(1:kproma,jrow)
             voc_emisfac(1:kproma,2,jrow) = mterp_emisfac(1:kproma,jrow)
             CALL VOC_emissions(kproma,         srfl(1:kproma,jrow)          &
                  ,cossza_2d(1:kproma,jrow),    lai(1:kproma,jrow)           &
                  ,lad(_RI_XYZ__(1:kproma,jrow,:)),      drymatter(1:kproma,jrow)     &
                  ,voc_emisfac(1:kproma,:,jrow),tslm1(1:kproma,jrow)         &
                  ,voc_emflux(1:kproma,:,jrow))
             ! calling fracturizer
             CALL terr13C_frac(kproma, &
                  terr13C_delta(1:kproma,jrow), 5,     &  ! delta & atoms count
                  voc_emflux(1:kproma,:,jrow),         &  ! total flux
                  I12ISOP_emflux(1:kproma,jrow),       &  ! 12C fractured flux
                  I13ISOP_emflux(1:kproma,jrow)        )  ! 13C fractured flux

          CASE('NO_yl95sl10')
             tsoil_top(1:kproma,jrow) = tsoil(1:kproma, 1, jrow)
             CALL NOemis_yl95sl10_pulsing(kproma, delta_time, nstep     &
                  , prl(1:kproma,jrow), prc(1:kproma,jrow), INT(OneDay) &
                  , pulseregime(1:kproma,jrow), pulseday(1:kproma,jrow) &
                  , pulse(1:kproma,jrow)                                &
                  , prec_hist(_RI_XY_N_(1:kproma,jrow,1:ndrydays)) )

              CALL NOemis_yl95sl10_crf(kproma, imonth                           &
                   , philat_2d(1:kproma,jrow)                                   &
                   , NOemclass_yl95sl10(_RI_XY_N_(1:kproma,jrow,1:ncl_yl95sl10))&
                   , lai_yl95sl10(1:kproma,jrow), crf(1:kproma,jrow))

              vsm(:,jrow) = 0.0_dp
              DO jp=1, kproma
                 vsm(jp,jrow) = root_depth(jp,jrow) * &
                      MAX(ws(jp,jrow)/wsmx(jp,jrow), 1.0_dp)
              END DO

              CALL NOemis_yl95sl10(kproma, imonth, philat_2d(1:kproma,jrow) &
                   , philon_2d(1:kproma, jrow)                              &
                   , prec_hist(_RI_XY_N_(1:kproma,jrow,1:ndrydays) )        &
                   , NOemclass_yl95sl10(_RI_XY_N_(1:kproma,jrow,1:ncl_yl95sl10))&
                   , fertil_yl95sl10(1:kproma, jrow), ws(1:kproma, jrow)    &
                   , vsm(1:kproma,jrow), tsoil(_RI_XYZ__(1:kproma,jrow,1))  &
                   , noslflux_yl95sl10(1:kproma, jrow)                      &
                   , noslflux_diag_yl95sl10 &
                                      (_RI_XY_N_(1:kproma,jrow,1:ncl_yl95sl10)))

              ! multiply soil flux with pulsing
              noslflux_yl95sl10(1:kproma,jrow) = &
                   noslflux_yl95sl10(1:kproma,jrow) * pulse(1:kproma, jrow)
              DO jl=1,ncl_yl95sl10
                 noslflux_diag_yl95sl10(_RI_XY_N_(1:kproma,jrow,jl)) = &
                      noslflux_diag_yl95sl10(_RI_XY_N_(1:kproma,jrow,jl)) &
                      * pulse(1:kproma, jrow)
              END DO

              ! change units: from ng(N)/(m^2*s) to molec/(m^2*s)
              ! for above canopy flux
              noemflux_yl95sl10(1:kproma, jrow) = &
                   ((noslflux_yl95sl10(1:kproma, jrow) * 1.e-9_dp) &
                   / 14.007_dp) * N_A
              ! multiply soil flux with CRF factors -> above canopy flux
              noemflux_yl95sl10(1:kproma, jrow) = &
                   noemflux_yl95sl10(1:kproma, jrow) * crf(1:kproma, jrow)

          CASE('NOpls')
             CALL NOpls_emissions(kproma,    INT(OneDay)           &
                  , delta_time, init_step,   nstep                 &
                  , prc(1:kproma,jrow),      prl(1:kproma,jrow)    &
                  , cpold(_RI_XY_N_(1:kproma,jrow,1:ndrydays))    &
                  , lspold(_RI_XY_N_(1:kproma,jrow,1:ndrydays))   &
                  , pulsing(1:kproma,jrow),  plsday(1:kproma,jrow) &
                  , plsdurat(1:kproma,jrow), cp(1:kproma,jrow)     &
                  , prectot(1:kproma,jrow)                         &
                  , lsp(1:kproma,jrow),      pls(1:kproma,jrow) )
           CASE('NO')
              prect(:)        = 0._dp
              IF (.not. ASSOCIATED(pls)) THEN
                 ! setting the total precipitation prectot to a value
                 ! indicating
                 ! the wet and dry season in the tropics, used to determine
                 ! the soil NO emission flux for tropical forest.
                 ! Normally this prectot is calculated in the pulsing routine
                 ! from the rainfall history
                 ! default setting to value of 151 mm, above threshold
                 DO jl =1,kproma
                    prect(jl)=151._dp
                    ! Northern hemisphere
                    IF (philat_2d(jl,jrow) >= 0._dp) THEN
                       IF (imonth > 4 .AND. imonth < 10) & ! 5 driest months
                            prect(jl)=0.
                    ELSE                             ! Southern hemisphere
                       IF (imonth < 4 .OR. imonth > 10) & ! 5 driest months
                            prect(jl)=0.
                    ENDIF
                 ENDDO
              ELSE
                 prect(1:kproma) = prectot(1:kproma, jrow)
              ENDIF
              CALL NO_emissions( kproma,     imonth                      &
                  , cultiv(1:kproma,jrow), fertil(1:kproma,jrow)         &
                  , tsoil(_RI_XYZ__(1:kproma,jrow,1)), ws(1:kproma,jrow) &
                  , prect(1:kproma)                                      &
                  , NOemis_w(1:kproma,jrow),   NOemis_d(1:kproma,jrow)   &
                  , noemisclass1(_RI_XYZ__(1:kproma,jrow,:))             &
                  , lai(1:kproma,jrow),    philat_2d(1:kproma,jrow)      &
                  , NO_emflux(1:kproma,jrow), noslflux(1:kproma,jrow)    &
                  , noslflux_diag(_RI_XY_N_(1:kproma,jrow,1:ncl_yl95))   )
              IF (ASSOCIATED(pls))  NO_emflux(1:kproma,jrow) =       &
                   NO_emflux(1:kproma,jrow) *  pls(1:kproma,jrow)
           CASE('DU')
              CALL dust_emissions(kproma, time_step_len,          &
                   wind10_2d(1:kproma,jrow), slm(1:kproma,jrow),  &
                   tslm1(1:kproma,jrow), prc(1:kproma,jrow),      &
                   prl(1:kproma,jrow), du_cla(1:kproma,jrow),     &
                   du_thr(1:kproma,jrow), du_src(1:kproma,jrow),  &
                   zprecipinsoil_2d(1:kproma,jrow),               &
                   du_emflux_B_ci(1:kproma,jrow))
           CASE('DU_tegen')
              CALL dust_emissions_tegen(kproma,                              &
                slf(1:kproma,jrow), glac(1:kproma,jrow), cvs(1:kproma,jrow), &
                ws(1:kproma,jrow), wsmx(1:kproma,jrow),                      &
                mat_s2(1:kproma,jrow), mat_s3(1:kproma,jrow),                &
                mat_s4(1:kproma,jrow), mat_s6(1:kproma,jrow),                &
                mat_psrc(1:kproma,jrow), k_fpar_eff(1:kproma,jrow),          &
                wind10_2d(1:kproma,jrow), cuscale,                           &
                du_emflux_T_ci(1:kproma,jrow), du_emflux_T_ai(1:kproma,jrow))
           CASE('DU_Astitha1')

               if (l_ducomp) then

                  CALL dust_emissions_DU_Astitha1(kproma &
                       , zairdens(1:kproma,nlev), zust_2d(1:kproma,jrow)       &
                       , slf(1:kproma,jrow), ws(1:kproma,jrow)                 &
                       , cvs(1:kproma,jrow), du_cla2(1:kproma,jrow)            &
                       , rdepth(1:kproma,jrow)                                 &
                       , dustsrc(_RI_XY_N_(1:kproma,jrow,1:nbiomes),1)         &
                       , lai_in(1:kproma,jrow), du_emflux_A1_ai(1:kproma,jrow) &
                       , du_emflux_A1_ci(1:kproma,jrow)                        &
                       , horflux(1:kproma,jrow)                                &
                       , ustarthr1(1:kproma,jrow)                              &
                       , du_nap(1:kproma,jrow), du_kp(1:kproma,jrow)           &
                       , du_capp(1:kproma,jrow), du_mgpp(1:kproma,jrow)        &
                       , du_misc(1:kproma,jrow)                                &
                       , du_nap_emflux_ai(1:kproma,  jrow)                     &
                       , du_nap_emflux_ci(1:kproma,  jrow)                     &
                       , du_kp_emflux_ai(1:kproma,   jrow)                     &
                       , du_kp_emflux_ci(1:kproma,   jrow)                     &
                       , du_capp_emflux_ai(1:kproma, jrow)                     &
                       , du_capp_emflux_ci(1:kproma, jrow)                     &
                       , du_mgpp_emflux_ai(1:kproma, jrow)                     &
                       , du_mgpp_emflux_ci(1:kproma, jrow)                     &
                       , du_misc_emflux_ai(1:kproma, jrow)                     &
                       , du_misc_emflux_ci(1:kproma, jrow))
               else
                  CALL dust_emissions_DU_Astitha1(kproma                       &
                       , zairdens(1:kproma,nlev), zust_2d(1:kproma,jrow)       &
                       , slf(1:kproma,jrow), ws(1:kproma,jrow)                 &
                       , cvs(1:kproma,jrow), du_cla2(1:kproma,jrow)            &
                       , rdepth(1:kproma,jrow)                                 &
                       , dustsrc(_RI_XY_N_(1:kproma,jrow,1:nbiomes),1)         &
                       , lai_in(1:kproma,jrow), du_emflux_A1_ai(1:kproma,jrow) &
                       , du_emflux_A1_ci(1:kproma,jrow), horflux(1:kproma,jrow)&
                       , ustarthr1(1:kproma,jrow))
               endif
           CASE('DU_Astitha2')
               if (l_ducomp) then
                  CALL dust_emissions_DU_Astitha2(kproma                       &
                       , zairdens(1:kproma,nlev), zust_2d(1:kproma,jrow)       &
                       , slf(1:kproma,jrow)                                    &
                       , ws(1:kproma,jrow), cvs(1:kproma,jrow)                 &
                       , du_cla2(1:kproma,jrow)                                &
                       , du_nap(1:kproma,jrow), du_kp(1:kproma,jrow)           &
                       , du_capp(1:kproma,jrow), du_mgpp(1:kproma,jrow)        &
                       , du_misc(1:kproma,jrow)                                &
                       , rdepth(1:kproma,jrow)                                 &
                       , lai_in(1:kproma,jrow)                                 &
                       , dustsrc(_RI_XY_N_(1:kproma,jrow,1:nbiomes),1)         &
                       , soiltexture(_RI_XY_N_(1:kproma,jrow,1:nclasses),1)    &
                       , du_emflux_A2_ai(1:kproma,jrow)                        &
                       , du_emflux_A2_ci(1:kproma,jrow)                        &
                       , ustarthr(_RI_XY_N_(1:kproma,jrow,1:maxsizes) )        &
                       , du_nap_emflux_ai(1:kproma, jrow)                      &
                       , du_nap_emflux_ci(1:kproma, jrow)                      &
                       , du_kp_emflux_ai(1:kproma, jrow)                       &
                       , du_kp_emflux_ci(1:kproma, jrow)                       &
                       , du_capp_emflux_ai(1:kproma, jrow)                     &
                       , du_capp_emflux_ci(1:kproma, jrow)                     &
                       , du_mgpp_emflux_ai(1:kproma, jrow)                     &
                       , du_mgpp_emflux_ci(1:kproma, jrow)                     &
                       , du_misc_emflux_ai(1:kproma, jrow)                     &
                       , du_misc_emflux_ci(1:kproma, jrow))
               else
                  CALL dust_emissions_DU_Astitha2(kproma                       &
                       , zairdens(1:kproma,nlev), zust_2d(1:kproma,jrow)       &
                       , slf(1:kproma,jrow), ws(1:kproma,jrow)                 &
                       , cvs(1:kproma,jrow)                                    &
                       , du_cla2(1:kproma,jrow)                                &
                       , rdepth(1:kproma,jrow)                                 &
                       , lai_in(1:kproma,jrow)                                 &
                       , dustsrc(_RI_XY_N_(1:kproma,jrow,1:nbiomes) ,1)        &
                       , soiltexture(_RI_XY_N_(1:kproma,jrow,1:nclasses) ,1)   &
                       , du_emflux_A2_ai(1:kproma,jrow)                        &
                       , du_emflux_A2_ci(1:kproma,jrow)                        &
                       , ustarthr(_RI_XY_N_(1:kproma,jrow,1:maxsizes) ))
               endif

       CASE('KKDU_Astitha')
          CALL dust_emissions_KKDU_Astitha1(kproma = kproma, &
                airdens = zairdens(1:kproma, nlev), &
                zust_2d = zust_2d(1:kproma, jrow), &
                slf = slf(1:kproma, jrow), &
                cvs = cvs(1:kproma, jrow), &
                kkdu_clay = kkdu_clay(1:kproma, jrow), &
                kkdu_topo = kkdu_topo(1:kproma, jrow), &
                kkdu_nap = kkdu_nap(1:kproma, jrow), &
                kkdu_kp = kkdu_kp(1:kproma, jrow), &
                kkdu_capp = kkdu_capp(1:kproma, jrow), &
                kkdu_mgpp = kkdu_mgpp(1:kproma, jrow), &
                kkdu_misc = kkdu_misc(1:kproma, jrow), &
                kkdu_mask = kkdu_mask(_RI_XY_N_(1:kproma,jrow,1:2)), &
                kkdu_lai = kkdu_lai(1:kproma, jrow), &
                duflux_ai = du_emflux_A1_ai(1:kproma, jrow), &
                duflux_ci = du_emflux_A1_ci(1:kproma, jrow), &
                kkdu_nap_emflux_ai = kkdu_nap_emflux_ai(1:kproma, jrow), &
                kkdu_nap_emflux_ci = kkdu_nap_emflux_ci(1:kproma, jrow), &
                kkdu_kp_emflux_ai = kkdu_kp_emflux_ai(1:kproma, jrow), &
                kkdu_kp_emflux_ci = kkdu_kp_emflux_ci(1:kproma, jrow), &
                kkdu_capp_emflux_ai = kkdu_capp_emflux_ai(1:kproma, jrow), &
                kkdu_capp_emflux_ci = kkdu_capp_emflux_ci(1:kproma, jrow), &
                kkdu_mgpp_emflux_ai = kkdu_mgpp_emflux_ai(1:kproma, jrow), &
                kkdu_mgpp_emflux_ci = kkdu_mgpp_emflux_ci(1:kproma, jrow), &
                kkdu_misc_emflux_ai = kkdu_misc_emflux_ai(1:kproma, jrow), &
                kkdu_misc_emflux_ci = kkdu_misc_emflux_ci(1:kproma, jrow), &
                horflux = horflux(1:kproma, jrow), &
                ustarthr = ustarthr1(1:kproma, jrow))

             CASE('SO2_ant')
                CALL so2_emissions(zdz(1:kproma,:) &
                     , so2_ant(_RI_XYZ__(1:kproma,jrow,:))  &
                     , so2_emflux(_RI_XYZ__(1:kproma, jrow,:)))

             CASE('BIOO')
                ! OLSON LUMPED EMISSIONS
                CALL bioaer_emissions_olson(kproma        &
                     , olson_emis(1)%ptr(1:kproma,jrow)   & ! _seas
                     , olson_emis(2)%ptr(1:kproma,jrow)   & ! _landice
                     , olson_emis(3)%ptr(1:kproma,jrow)   & ! _deserts
                     , olson_emis(4)%ptr(1:kproma,jrow)   & ! _forests
                     , olson_emis(5)%ptr(1:kproma,jrow)   & ! _grasslands
                     , olson_emis(6)%ptr(1:kproma,jrow)   & ! _crops
                     , olson_emis(7)%ptr(1:kproma,jrow)   & ! _wetlands
                     , olson_emis(8)%ptr(1:kproma,jrow)   & ! _shrubs
                     , olson_emis(9)%ptr(1:kproma,jrow)   & ! _coastal
                     , olson_emis(10)%ptr(1:kproma,jrow)  & ! _urban
                     , olson_emis(11)%ptr(1:kproma,jrow)  & ! _tundra
                     , olson(_RI_XY_N_(1:kproma,jrow,:))  &
                     , nolclass                           &
                     )

             CASE('BIOM')
                ! MODIS EMISSIONS
                CALL bioaer_emissions_modis(kproma       &
                     , modis_emis(1)%ptr(1:kproma,jrow)  & ! water
                     , modis_emis(2)%ptr(1:kproma,jrow)  & ! ever_need
                     , modis_emis(3)%ptr(1:kproma,jrow)  & ! ever_broad
                     , modis_emis(4)%ptr(1:kproma,jrow)  & ! deci_need
                     , modis_emis(5)%ptr(1:kproma,jrow)  & ! deci_broad
                     , modis_emis(6)%ptr(1:kproma,jrow)  & ! mixed_forest
                     , modis_emis(7)%ptr(1:kproma,jrow)  & ! closed_shrubs
                     , modis_emis(8)%ptr(1:kproma,jrow)  & ! open_shrubs
                     , modis_emis(9)%ptr(1:kproma,jrow)  & ! woody_savannas
                     , modis_emis(10)%ptr(1:kproma,jrow) & ! savannas
                     , modis_emis(11)%ptr(1:kproma,jrow) & ! grasslands
                     , modis_emis(12)%ptr(1:kproma,jrow) & ! perm_wetlands
                     , modis_emis(13)%ptr(1:kproma,jrow) & ! crops
                     , modis_emis(14)%ptr(1:kproma,jrow) & ! urban
                     , modis_emis(15)%ptr(1:kproma,jrow) & ! crop_nature
                     , modis_emis(16)%ptr(1:kproma,jrow) & ! snow_ice
                     , modis_emis(17)%ptr(1:kproma,jrow) & ! barren
                     , modis_emis(18)%ptr(1:kproma,jrow) & ! unclass
                     , modis(_RI_XY_N_(1:kproma,jrow,:)) &
                     , nmodisclass                       &
                     )

!##########################################
!# LAI-BASED BIOAEROSOL EMISSIONS         #
!##########################################
          CASE('BIOL')
                CALL bio_emissions_lai(kproma                        &
                     , heald_emis(1:kproma,jrow)                     &
                     , js_emis(1:kproma,jrow)                        &
                     , hummel_emis(1:kproma,jrow)                    &
                     , modis_lai(1:kproma,jrow)                      &
                     , qm1_3d(_RI_XYZ__(1:kproma,jrow,nlev))         &
                     , tm1_3d(_RI_XYZ__(1:kproma,jrow,nlev))         &
                     , philat_2d(1:kproma,jrow)                      &
                     , imonth                                        &
                     )
!##########################################
! Organic seasalt emissions for the accumulation mode only
             CASE('SS_POC_AQUA')
               DO ni=1,n_ss_poc_aqua
                 ss_flux        => emis_poc_aqua(ni)%ss_flx(1:kproma,jrow)
                 emis_poc       => emis_poc_aqua(ni)%emis_flx(1:kproma,jrow)
                 emisflxsum_poc => emis_poc_aqua(ni)%emis_flxsum(1:kproma,jrow)

                 CALL POC_EMIS_SS(ss_flux(1:kproma), slf(1:kproma,jrow),  &
                      poc_aqua(1:kproma,jrow), emis_poc(1:kproma) )

                 emisflxsum_poc(1:kproma) = emisflxsum_poc(1:kproma)      &
                                          + emis_poc(1:kproma) * delta_time
               ENDDO
             CASE('SS_POC_SWIFS')
               DO ni=1,n_ss_poc_swifs
                 ss_flux        => emis_poc_swifs(ni)%ss_flx(1:kproma,jrow)
                 emis_poc       => emis_poc_swifs(ni)%emis_flx(1:kproma,jrow)
                 emisflxsum_poc => emis_poc_swifs(ni)%emis_flxsum(1:kproma,jrow)

                 CALL POC_EMIS_SS(ss_flux(1:kproma), slf(1:kproma,jrow),  &
                      poc_seawifs(1:kproma,jrow), emis_poc(1:kproma) )

                 emisflxsum_poc(1:kproma) = emisflxsum_poc(1:kproma)      &
                                          + emis_poc(1:kproma) * delta_time
               ENDDO
             CASE('SS_WIOC_AQUA')
               DO ni=1,n_ss_wioc_aqua
                 ss_flux        => emis_wioc_aqua(ni)%ss_flx(1:kproma,jrow)
                 emis_poc       => emis_wioc_aqua(ni)%emis_flx(1:kproma,jrow)
                 emisflxsum_poc => emis_wioc_aqua(ni)%emis_flxsum(1:kproma,jrow)

                 CALL WIOC_EMIS_SS(ss_flux(1:kproma), slf(1:kproma,jrow),  &
                      chlor_a_aqua_wioc(1:kproma,jrow), emis_poc(1:kproma) )

                 emisflxsum_poc(1:kproma) = emisflxsum_poc(1:kproma)       &
                                          + emis_poc(1:kproma) * delta_time
               ENDDO
             CASE('SS_WIOC_BLEN')
               DO ni=1,n_ss_wioc_blend
                 ss_flux        => emis_wioc_blend(ni)%ss_flx(1:kproma,jrow)
                 emis_poc       => emis_wioc_blend(ni)%emis_flx(1:kproma,jrow)
                 emisflxsum_poc => emis_wioc_blend(ni)%emis_flxsum(1:kproma,jrow)

                 CALL WIOC_EMIS_SS(ss_flux(1:kproma), slf(1:kproma,jrow),  &
                      chlor_a_blend_wioc(1:kproma,jrow), emis_poc(1:kproma) )

                 emisflxsum_poc(1:kproma) = emisflxsum_poc(1:kproma)       &
                                          + emis_poc(1:kproma) * delta_time
               ENDDO

       CASE('AirSnow')
          ! Initaliation
          snow_air_flux_br2(1:kproma,jrow) = 0._dp
          snow_air_flux_brcl(1:kproma,jrow) = 0._dp
          CALL airsnow_emissions( kproma, &
               snow_air_flux_br2(1:kproma,jrow),  &
               snow_air_flux_brcl(1:kproma,jrow), &
               tsurf_2d(1:kproma,jrow),           &
               cvs(1:kproma,jrow),                &
               seaice(1:kproma,jrow),             &
               cossza_2d(1:kproma,jrow),          &
               sic_multi_year(1:kproma,jrow),     &
               ddepflux_HOBr(1:kproma,jrow),      &
               ddepflux_BrNO3(1:kproma, jrow),    &
               ddepflux_HBr(1:kproma,jrow),       &
               ddepflux_O3(1:kproma,jrow) )

             !CASE()
             !###############################
             ! ### add new emissions here ###
             !###############################
          CASE DEFAULT
       END SELECT

    END DO emistype_loop4
    ! STEP 2: PUT EMISSION FLUXES INTO TRACERS

    fluxes_loop: DO jf=1, NFLUXES

       ! GRIDPOINT TRACERS
       gp_tracer_loop: DO jt=1, XF2T(jf)%ngpt

          IF (XF2T(jf)%idt_gp(jt) == 0) CYCLE ! TRACER NOT AVAILABLE
          idt = XF2T(jf)%idt_gp(jt) ! um_ak_20100802 (RI)

         ! SELECT EMISSION METHOD
          SELECT CASE(XF2T(jf)%mgp(jt))
          CASE(0)
             ! DO NOTHING
          CASE(1)
             IF (SIZE(XF2T(jf)%ptr,_IZ_XYZ__) > 1) THEN
                ! 3D
#ifndef MESSYTENDENCY
                pxtte(_RI_X_ZN_(1:kproma,:,idt)) =    &
                     pxtte(_RI_X_ZN_(1:kproma,:,idt)) &
                     + ( XF2T(jf)%ptr(_RI_XYZ__(1:kproma,jrow,:))   &
                     / (zairdens(1:kproma,:) * uconv) )   &
                     * XF2T(jf)%efact_gp(jt)
#else
                zxtte(:,:) = 0._dp
                zxtte(1:kproma,:) = XF2T(jf)%ptr(_RI_XYZ__(1:kproma,jrow,:))   &
                     / (zairdens(1:kproma,:) * uconv)                &
                     * XF2T(jf)%efact_gp(jt)
                CALL mtend_add_l(my_handle, jt, px=zxtte)
#endif
             ELSE
                ! 2D
#ifndef MESSYTENDENCY
                pxtte(_RI_X_ZN_(1:kproma,nlev,idt)) =    &
                     pxtte(_RI_X_ZN_(1:kproma,nlev,idt)) &
                     + ( XF2T(jf)%ptr(_RI_XYZ__(1:kproma,jrow,1))       &
                     / (zscale(1:kproma) * uconv) )           &
                     * XF2T(jf)%efact_gp(jt)
#else
                zxtte(:,:) = 0._dp
                zxtte(1:kproma,nlev) = XF2T(jf)%ptr(_RI_XYZ__(1:kproma,jrow,1))  &
                     / (zscale(1:kproma) * uconv)                      &
                     * XF2T(jf)%efact_gp(jt)
                CALL mtend_add_l(my_handle, jt, px=zxtte)
#endif
             END IF
#if defined(ECHAM5) || defined(CESM1)
          CASE(2)
             pxtems(_RI_XYZN_(1:kproma,jrow,1,idt)) =    &
                  pxtems(_RI_XYZN_(1:kproma,jrow,1,idt)) &
                  + ( XF2T(jf)%ptr(_RI_XYZ__(1:kproma,jrow,1))      &
                  /  uconv )                              &
                  * XF2T(jf)%efact_gp(jt)
#endif
          CASE DEFAULT
             ! DO NOTHING (CHECKED ABOVE)
          END SELECT

       END DO gp_tracer_loop

    END DO fluxes_loop

    NULLIFY(zairdens)
    DEALLOCATE(zdz)
    DEALLOCATE(zscale)
#ifdef MESSYTENDENCY
    DEALLOCATE(zxtte)
#endif
  END SUBROUTINE onemis_vdiff
! ---------------------------------------------------------------------------

! ---------------------------------------------------------------------------
  SUBROUTINE onemis_global_end

    IMPLICIT NONE
!!#D attila +
#ifdef ECHAM5
    IF (L_LG) CALL onemis_global_end_lg
#endif
!!#D attila -

  END SUBROUTINE onemis_global_end
! ---------------------------------------------------------------------------

!!#D attila +
#ifdef ECHAM5
! ----------------------------------------------------------------
  SUBROUTINE onemis_global_end_lg

    ! ECHAM5/MESSy
    USE messy_main_tracer_mem_bi,   ONLY: pxtte => qxtte_a, NCELL
    USE messy_main_grid_def_mem_bi, ONLY: nlev, ngpblks, kproma &
                                        , nproma, npromz
    USE messy_main_grid_def_bi,     ONLY: deltaz
    USE messy_main_data_bi,         ONLY: zairdens => rho_air_dry_3d

    USE messy_attila_tools_e5,    ONLY: gpsfemis2lgemis_e5, gp2lg_e5
    ! MESSy
    USE messy_main_constants_mem, ONLY: g, M_air, N_A

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER             :: substr = 'onemis_global_end_lg'
    REAL(DP), PARAMETER                     :: uconv = N_A * 1.0e3_DP/M_air
    REAL(DP), DIMENSION(:,:,:), ALLOCATABLE :: zdz       ! layer thickness
    REAL(DP), DIMENSION(:,:,:), ALLOCATABLE :: zscale    !
    REAL(DP), DIMENSION(:,:,:), POINTER     :: zxtte_gp    => NULL()
    REAL(DP), DIMENSION(:,:),   POINTER     :: zxtte_gpsf  => NULL()
    REAL(DP), DIMENSION(:,:),   POINTER     :: rest_2d     => NULL()
    REAL(DP), DIMENSION(:),     POINTER     :: zxtte_lg    => NULL()
    INTEGER                                 :: zjrow
    INTEGER                                 :: jt, jf
    INTEGER                                 :: nm1count
    INTEGER                                 :: idt
    INTEGER                                 :: ml1, mlev

    ! MEMORY FOR GP-TENDENCY
    ALLOCATE(zxtte_gp(nproma, nlev, ngpblks))
    ! MEMORY FOR LG-TENDENCY
    ALLOCATE(zxtte_lg(NCELL))

    ! CALCULATE AIR DENSITY AND LAYER THICKNESS
    ALLOCATE(zdz(nproma, nlev, ngpblks))
    zdz(:,:,:) = 1.0_DP
    ALLOCATE(zscale(nproma, nlev, ngpblks))
    zscale(:,:,:) = 1.0_DP
    !
    DO zjrow=1, ngpblks
       IF ( zjrow == ngpblks ) THEN
          kproma = npromz
       ELSE
          kproma = nproma
       END IF
       !
       zdz(1:kproma,2:nlev,zjrow) = deltaz(1:kproma,2:nlev,zjrow)
    END DO
    ! ONLY 2D SURFACE EMISSIONS (SO FAR ...!)
    zscale(:,2:nlev,:) = zairdens(:,2:nlev,:)*zdz(:,2:nlev,:)

    fluxes_loop: DO jf=1, NFLUXES

       IF (SIZE(XF2T(jf)%ptr,_IN_XY_N_) > 1) THEN
          ! 3D
          ml1  = 1
          mlev = nlev
       ELSE
          ! 2D
          ml1  = nlev
          mlev = 1
       END IF

       nm1count = 0
       lg_tracer_loop: DO jt=1, XF2T(jf)%nlgt

          IF (XF2T(jf)%idt_lg(jt) == 0) CYCLE ! TRACER DOES NOT EXIST

          IF (XF2T(jf)%mlg(jt) == 0) CYCLE    ! LG-METHOD = 0

          zxtte_gp(:,:,:) = 0._dp
          zxtte_gp(:,ml1:nlev,:) =              &
               XF2T(jf)%ptr(:,1:mlev,:)         &
               / (zscale(:,ml1:nlev,:) * uconv) &
               * XF2T(jf)%efact_lg(jt)
          zxtte_gpsf => zxtte_gp(:,nlev,:)

          SELECT CASE(XF2T(jf)%mlg(jt))
          CASE(0)
             ! CANNOT BE REACHED (SEE ABOVE)
          CASE(1)
             !
             nm1count = nm1count + 1 ! count pointer-index to rest-flux
             !
             IF (SIZE(XF2T(jf)%ptr,_IZ_XYZ__) > 1) THEN
                ! 3D
                CALL gp2lg_e5(zxtte_gp, zxtte_lg             &
                     , gprl=XF2T(jf)%ptr_rest(nm1count)%ptr  &
                     , lmcons = .true. )
             ELSE
                ! 2D
                rest_2d => XF2T(jf)%ptr_rest(nm1count)%ptr(:,nlev,:)
                CALL gpsfemis2lgemis_e5(zxtte_gpsf, zxtte_lg   &
                     , XF2T(jf)%mlg(jt)                        &  ! = 1
                     , gprl=rest_2d                            &
                     , lmcons = .true. )
             END IF
             !
          CASE(2, 3, 4)
             !
             CALL gpsfemis2lgemis_e5(zxtte_gpsf, zxtte_lg &
                  , XF2T(jf)%mlg(jt)                      &  ! = 2,3,4
                  , lmcons = .true. )
             !
          CASE DEFAULT
             !
             ! ERROR
             CALL error_bi( 'UNKNOWN EMISSION METHOD',substr)
             !
          END SELECT

          ! SET TRACER INDEX
          idt = XF2T(jf)%idt_lg(jt)

          ! ADD TRACER TENDENCY
          pxtte(:,idt) = pxtte(:,idt) + zxtte_lg(:)

       END DO lg_tracer_loop

    END DO fluxes_loop

    ! CLEAN UP
    DEALLOCATE(zdz)
    DEALLOCATE(zscale)
    !
    DEALLOCATE(zxtte_gp)
    DEALLOCATE(zxtte_lg)

  END SUBROUTINE onemis_global_end_lg
! ----------------------------------------------------------------
#endif
!!#D attila -

! ---------------------------------------------------------------------------
SUBROUTINE onemis_free_memory

  INTEGER :: jf

  CHARACTER(len=*), PARAMETER :: substr='onemis_free_memory'

  CALL start_message_bi(modstr,'FREE EMISSET MEMORY',substr)

  DO jf=1, NFLUXES
     ! CHANNEL OBJECT; DO NOT DEALLOCATE
     ! GP
     IF (ASSOCIATED(XF2T(jf)%tracer_gp)) DEALLOCATE(XF2T(jf)%tracer_gp)
     IF (ASSOCIATED(XF2T(jf)%mgp))       DEALLOCATE(XF2T(jf)%mgp)
     IF (ASSOCIATED(XF2T(jf)%efact_gp))  DEALLOCATE(XF2T(jf)%efact_gp)
     IF (ASSOCIATED(XF2T(jf)%idt_gp))    DEALLOCATE(XF2T(jf)%idt_gp)
!!#D attila +
#ifdef ECHAM5
     ! LG
     IF (ASSOCIATED(XF2T(jf)%ptr_rest)) THEN
        ! CHANNEL OBJECT:
        DEALLOCATE(XF2T(jf)%ptr_rest)
     END IF
     IF (ASSOCIATED(XF2T(jf)%tracer_lg)) DEALLOCATE(XF2T(jf)%tracer_lg)
     IF (ASSOCIATED(XF2T(jf)%mlg))       DEALLOCATE(XF2T(jf)%mlg)
     IF (ASSOCIATED(XF2T(jf)%efact_lg))  DEALLOCATE(XF2T(jf)%efact_lg)
     IF (ASSOCIATED(XF2T(jf)%idt_lg))    DEALLOCATE(XF2T(jf)%idt_lg)
#endif
!!#D attila -
  END DO
  !
  CALL end_message_bi(modstr,'FREE EMISSET MEMORY',substr)

END SUBROUTINE onemis_free_memory

! ---------------------------------------------------------------------------
! ************************************************************************
! PRIVATE ECHAM-5 INTERFACE ROUTINES
! ************************************************************************
! ---------------------------------------------------------------------------

! ----------------------------------------------------------------------
  SUBROUTINE onemis_read_nml_cpl(status, iou)

    ! onemis MODULE ROUTINE (ECHAM-5 INTERFACE, PRIVATE)
    !
    ! read namelist for 'coupling' to ECHAM5
    !
    ! Author: Patrick Joeckel, MPICH, Mar 2004

    ! MESSy
    USE messy_main_tools,         ONLY: read_nml_open, read_nml_check &
                                      , read_nml_close
    ! ECHAM5/MESSy
    USE messy_main_tracer_mem_bi, ONLY: NGCELL

    IMPLICIT NONE

    ! I/O
    INTEGER, INTENT(OUT) :: status     ! error status
    INTEGER, INTENT(IN)  :: iou        ! I/O unit

    ! (LOCAL) NAMELIST VARIABLES
    CHARACTER(LEN=*), PARAMETER :: substr = 'onemis_read_nml_cpl'

    ! LOCAL
    LOGICAL              :: lex      ! file exists ?
    INTEGER              :: fstat    ! file status

    status = 1

    CALL read_nml_open(lex, substr, iou, 'CPL', modstr)
    IF (.not.lex) RETURN    ! <modstr>.nml does not exist

    READ(iou, NML=CPL, IOSTAT=fstat)
    CALL read_nml_check(fstat, substr, iou, 'CPL', modstr)
    IF (fstat /= 0) RETURN  ! error while reading namelist

    ! DIAGNOSE NAMELIST AND SET GLOBAL SWITCHES
    !
    ! CHECK NAMELIST
    IF (L_LG .AND. (NGCELL <= 0)) THEN
!!#D attila +
       WRITE(*,*) 'L_LG = T in namelist'
       WRITE(*,*) 'However no Lagrangian scheme activated ...'
       WRITE(*,*) ' ... setting L_LG = F'
!!#D attila -
       L_LG = .false.
    END IF

!!#D attila +
    IF (L_LG) THEN
       WRITE(*,*) 'EMISSIONS IN LAGRANGIAN SPACE: ON'
    ELSE
       WRITE(*,*) 'EMISSIONS IN LAGRANGIAN SPACE: OFF'
    END IF
!!#D attila -

    CALL read_nml_close(substr, iou, modstr)

    status = 0  ! no ERROR

  END SUBROUTINE onemis_read_nml_cpl
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
  SUBROUTINE onemis_read_nml_cpl_import(status, iou)

    ! onemis MODULE ROUTINE (ECHAM-5 INTERFACE, PRIVATE)
    !
    ! read namelist for 'coupling' to ECHAM5
    !
    ! Author: Patrick Joeckel, MPICH, Mar 2004

    ! MESSy
    USE messy_main_tools,         ONLY: read_nml_open, read_nml_check &
                                      , read_nml_close

    IMPLICIT NONE

    ! I/O
    INTEGER, INTENT(OUT) :: status     ! error status
    INTEGER, INTENT(IN)  :: iou        ! I/O unit

    ! (LOCAL) NAMELIST VARIABLES
    CHARACTER(LEN=*), PARAMETER :: substr = 'onemis_read_nml_cpl'

    ! LOCAL
    LOGICAL              :: lex      ! file exists ?
    INTEGER              :: fstat    ! file status

    status = 1

    CALL read_nml_open(lex, substr, iou, 'CPL_IMPORT', modstr)
    IF (.not.lex) RETURN    ! <modstr>.nml does not exist

    READ(iou, NML=CPL_IMPORT, IOSTAT=fstat)
    CALL read_nml_check(fstat, substr, iou, 'CPL_IMPORT', modstr)
    IF (fstat /= 0) RETURN  ! error while reading namelist

    ! DIAGNOSE NAMELIST AND SET GLOBAL SWITCHES
    !
    CALL read_nml_close(substr, iou, modstr)

    status = 0  ! no ERROR

  END SUBROUTINE onemis_read_nml_cpl_import
! ----------------------------------------------------------------------

!*****************************************************************************
END MODULE messy_onemis_si
!*****************************************************************************

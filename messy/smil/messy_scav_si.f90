#include "messy_main_ppd_bi.inc"

MODULE messy_scav_si

!  AUTHOR:   Holger Tost, MPI Chemie, Mainz
!            last modified 13.04.2006
!  CONTRIBUTORS:
!            Christopher Kaiser, DLR, Oberpfaffenhofen
!              - Adaptations for MADE and MADE3 (2013-2016)
!              - Equipped with TENDENCY (2016)
!            Mattia Righi, DLR, Oberpfaffenhofen
!              - Adaptations for MADE

  ! ECHAM5/MESSy
  USE messy_main_grid_def_mem_bi, ONLY: nproma, ngpblks
  USE messy_main_mpi_bi,        ONLY: p_pe,p_bcast,p_io
  USE messy_main_tracer_mem_bi, ONLY: GPTRSTR, ntrac_gp, ti_gp, &
                                      LGTRSTR, ntrac_lg, ti_lg, &
                                      t_trinfo_tp, t_trinfo

  ! MESSy
  USE messy_main_tracer,        ONLY: I_convect, I_aerosol_mode, &
                                      I_aerosol_method, S_aerosol_model, &
                                      I_aerosol_sol, I_scav, &
                                      I_aerosol_hetice 

  USE messy_main_constants_mem, ONLY: I4, dp, STRLEN_MEDIUM, STRLEN_ULONG
  USE messy_scav
  USE messy_scav_mem
  USE messy_main_channel,       ONLY: STRLEN_OBJECT, STRLEN_CHANNEL
  USE messy_scav_inter,         ONLY: in_flux, out_flux, in_cloud, &
                                      out_cloud, evapo, evapo_aer, &
                                      thres_val, gboxarea_2d,      &
                                      melting, sedi_cloud,         &
                                      aer_scav_drv, cloud2trac,    &
                                      cloudtrac2aer, aer2cloudtrac
  USE messy_scav_aer,           ONLY: aermodel, max_mode
  USE messy_main_tools,         ONLY: PTR_3D_ARRAY

#ifdef MECCA_TAG
  ! this code is included only when tagging is used
  ! mecca_tag routine to process scav. for related tagging tracers
  USE messy_mecca_tag_si,       ONLY: mecca_tag_calc_xtte4scav
#endif

  IMPLICIT NONE

  SAVE

  INTRINSIC ::  EXP, MIN, INT, REAL, ANY, MAXVAL, ABS, TRIM, MAX, EPSILON
  LOGICAL :: cvtrans_use

  CHARACTER(LEN=STRLEN_OBJECT), PUBLIC  :: &
! char for convective precipitation
                                rain_cv(2)              ,&  
!      for            snow flux
                                snow_cv(2)              ,&  
!      for            cloud cover
                                ccover(2)               ,&  
!      for            freshly formed precipitation
                                cvprec(2)               ,&  
!      for            freshly formed snow
                                cvsnow(2)               ,&  
!      for large scale cloud cover
                                lcover(2)               ,&  
!      for             precipitating cover  
                                rcover(2)               ,&  
!      for             precipitation formation rate  
                                ratep(2)                ,&  
!      for             precipitation flux       
                                prec(2)                 ,&  
!      for             snow/ice formation rate 
                                ratesi(2)               ,&  
!      for             snow/ice precipitation flux  
                                fsi(2)                  ,&  
!      for             liquid water content
                                lwc(2)                  ,&  
!      for             snow / ice water content
                                iwc(2)                  ,&  
!      for             melting of frozen precip
                                imelt(2)                ,&  
!      for             ice sedimentation flux  
                                isedi(2)                ,&  
!      for grid mass
                                mass(2)                 ,&  
!      for grid volume
                                vol(2)                  ,&  
!      for pressure
                                press(2), pressi(2)     ,&  
!      for grid box area
                                boxarea(2)              ,&  
!      for photolysis rate
                                photol(2)               ,&    
!      for convective lwc
                                cvlwc(2)                ,&
!      for convective precip formation rate
                                cvrform(2)              ,& 
!      for convective lwc
                                cviwc(2)                ,&
!      for convective precip formation rate
                                cvsform(2)              ,& 
!      for nucleation scavenging number fraction from the 
!          cloud module(= activated fraction)
                                nfrac_nuc(2)            ,& 
!      for nucleation scavenging mass fraction from the 
!          cloud module(= activated fraction)
                                mfrac_nuc(2)              

! STRING FOR AEROSOL MODULE NAME
  CHARACTER(LEN=32), PUBLIC :: AERMOD_GP, AERMOD_LG

! MADE3 parameters (mainly used for mode assignment of cloud residuals)
  TYPE t_made3_nml
    INTEGER          :: ks = 0
    INTEGER          :: km = 0
    INTEGER          :: ki = 0
    INTEGER          :: as = 0
    INTEGER          :: am = 0
    INTEGER          :: ai = 0
    INTEGER          :: cs = 0
    INTEGER          :: cm = 0
    INTEGER          :: ci = 0
    INTEGER          :: nmod = 0
    INTEGER          :: i_mode = 0
    CHARACTER(LEN=2) :: csubname = ''
    INTEGER          :: evap_mode = 0
    INTEGER          :: n_resmod = 0
  END TYPE t_made3_nml
  TYPE(t_made3_nml), PUBLIC :: made3params

  INTEGER, PUBLIC :: idt_scav_NH3        = 0
  INTEGER, PUBLIC :: idt_scav_NO3mres_cs = 0
  INTEGER, PUBLIC :: idt_scav_SO4res_cs  = 0
  INTEGER, PUBLIC :: idt_scav_nh4pres_cs = 0
  INTEGER, PUBLIC :: idt_scav_Clmres_cs  = 0 
  INTEGER, PUBLIC :: idt_scav_Hpres_cs   = 0
  
  REAL(dp), DIMENSION(:,:,:), POINTER :: bclwc_cv => NULL() 
  REAL(dp), DIMENSION(:,:,:), POINTER :: bclwc_ls => NULL() 
  REAL(dp), DIMENSION(:,:,:), POINTER :: bc_cov   => NULL()

! for output of HNO3 partitioning parameters
  REAL(dp), DIMENSION(:,:,:), POINTER :: phi_i_cv => NULL()
  REAL(dp), DIMENSION(:,:,:), POINTER :: phi_i_ls => NULL()
  REAL(dp), DIMENSION(:,:,:), POINTER :: mju_i_cv => NULL()
  REAL(dp), DIMENSION(:,:,:), POINTER :: mju_i_ls => NULL()
  REAL(dp), DIMENSION(:,:,:), POINTER :: iwc_cv   => NULL()
  REAL(dp), DIMENSION(:,:,:), POINTER :: iwc_T_cv => NULL()
  REAL(dp), DIMENSION(:,:,:), POINTER :: iwc_ls   => NULL()
  REAL(dp), DIMENSION(:,:,:), POINTER :: iwc_T_ls => NULL()

  LOGICAL                                    :: L_TEND = .FALSE.
  CHARACTER(LEN=3000)                        :: TE_STRING = ''
  INTEGER, DIMENSION(:), POINTER             :: te_idt => NULL()
  TYPE(PTR_3D_ARRAY), DIMENSION(:), POINTER  :: tte_scav => NULL()
  INTEGER                                    :: n_te

#ifdef MESSYTENDENCY
  INTEGER :: tend_cv, tend_ls, tend_ev, tend_now
#endif

  LOGICAL :: column = .FALSE.

  INTEGER, ALLOCATABLE :: kpp_out_nr(:)

  INTEGER              :: n_out, n_out_aer

  LOGICAL :: l_trac_loc = .FALSE.   ! aerosol tracer locally defined

  ! CPL NAMELIST SWITCH
  REAL(dp), DIMENSION(2) :: LIQFRAC

  INTEGER  :: nsets = 0   ! number of tracer sets
  INTEGER  :: ntrac
  INTEGER  :: max_lev_scav
  REAL(dp) :: altitude

  TYPE(t_trinfo_tp), DIMENSION(:), POINTER :: ti => NULL()

! grid point approach

  TYPE w_flux_set
    REAL(dp), DIMENSION(:,:), POINTER :: flux_ls_sum  => NULL()
    REAL(dp), DIMENSION(:,:), POINTER :: flux_cv_sum  => NULL()
  END TYPE w_flux_set
  TYPE(w_flux_set), DIMENSION(:,:), POINTER :: wet_flx  => NULL()

  TYPE w_flux_set_aer
    REAL(dp), DIMENSION(:,:), POINTER :: flux_ls_aer_sum  => NULL()
    REAL(dp), DIMENSION(:,:), POINTER :: flux_cv_aer_sum  => NULL()
    INTEGER                           :: aer_out_nr
  END TYPE w_flux_set_aer
  TYPE(w_flux_set_aer), DIMENSION(:,:), POINTER :: wet_flx_aer  => NULL()

  TYPE spec_flux_set
    REAL(dp), POINTER, DIMENSION(:,:) :: nitrate_flx     => NULL()
    REAL(dp), POINTER, DIMENSION(:,:) :: nitrate_flx_sum => NULL()
    REAL(dp), POINTER, DIMENSION(:,:) :: sulfate_flx     => NULL()
    REAL(dp), POINTER, DIMENSION(:,:) :: sulfate_flx_sum => NULL()
    REAL(dp), POINTER, DIMENSION(:,:) :: ammoni_flx      => NULL()
    REAL(dp), POINTER, DIMENSION(:,:) :: ammoni_flx_sum  => NULL()
  END TYPE spec_flux_set
  TYPE(spec_flux_set), DIMENSION(:), POINTER :: special_flux => NULL()
  
 ! p5 -ptr for ls scav gas species
  REAL(dp), DIMENSION(:,:,:,:,:), POINTER :: spec_field_ls => NULL()
 ! p5 -ptr for cv scav gas species
  REAL(dp), DIMENSION(:,:,:,:,:), POINTER :: spec_field_cv => NULL()  
  
  REAL(dp), DIMENSION(:,:,:),   POINTER :: xtte_scav => NULL()

 ! KPP_field
  TYPE kpp_spec_arrays
    REAL(dp), POINTER, DIMENSION(:,:,:,:) :: cloud_field  => NULL()
    REAL(dp), POINTER, DIMENSION(:,:,:,:) :: precip_field => NULL()
  END TYPE kpp_spec_arrays
  TYPE(kpp_spec_arrays), DIMENSION(:,:), POINTER :: ls_field => NULL()
  TYPE(kpp_spec_arrays), DIMENSION(:,:), POINTER :: cv_field => NULL()

 ! set for process related coefficients
  TYPE all_proc
    REAL(dp), DIMENSION(:,:,:,:,:), POINTER :: evafrac => NULL()
  END TYPE all_proc
  ! set of parameterisation type (LS or CV)
  TYPE param
    TYPE(all_proc), DIMENSION(:), POINTER :: proc => NULL()
  END TYPE param
  TYPE trac_set
    TYPE(param), DIMENSION(:), POINTER    :: para => NULL()
  END TYPE trac_set
  TYPE(trac_set), DIMENSION(:), POINTER   :: set => NULL()


 ! for kpp field allocation (MAX(LSPEC,ISPEC))
  INTEGER :: KSPEC   

! l_trac_aer: make aerosol tracers for kpp species;  using this switch, 
!             but not l_trac_aer_all will create tracers only for charged
!             compounds
#ifdef COSMO
  LOGICAL :: l_trac_cloud = .TRUE.
#else
  LOGICAL :: l_trac_cloud = .FALSE.
#endif
  LOGICAL :: l_trac_aer = .FALSE.
  LOGICAL :: l_trac_aer_all = .FALSE. ! make tracer tracers for all kpp species

CONTAINS
!==============================================================================
!==============================================================================

  SUBROUTINE  scav_initialize

    ! SCAVENGING MODULE ROUTINE (ECHAM-5 INTERFACE)
    !
    ! INITIALIZATION OF GLOBAL VARIABLES FROM NAMELIST
    ! IN PARALLEL ENVIRONMENT
    ! INITIALIZATION OF XTSURF SPECIFIC EVENTS FROM NAMELIST
    ! IN PARALLEL ENVIRONMENT
    ! 

    ! Author, H. Tost, MPICH, 31-10-2003

    ! ECHAM5/MESSy
    USE messy_main_mpi_bi,       ONLY: p_parallel_io
    USE messy_main_blather_bi,   ONLY: error_bi, info_bi, warning_bi
    ! MESSy
    USE messy_main_tools,        ONLY: find_next_free_unit
#ifdef MESSYTENDENCY
    USE messy_main_tendency_bi,  ONLY: mtend_get_handle
#endif
    USE messy_scav_liq,          ONLY: use_schwartz
    USE messy_scav_aer,          ONLY: coeff_para

    USE MESSY_SCAV_L_KPP,        ONLY : LICNTRL => ICNTRL
    USE MESSY_SCAV_I_KPP,        ONLY : IICNTRL => ICNTRL
    IMPLICIT NONE

    ! LOCAL
    INTEGER     :: iou         ! I/O unit
    INTEGER     :: status      ! status
    CHARACTER(LEN=*), PARAMETER::substr='scav_initialize'


    ! INITIALIZE MAIN-CTRL
    IF (p_parallel_io) THEN
       iou = find_next_free_unit(100,200)
       ! *** CALL SCAV CORE ROUTINE:
       CALL init_scav(iou, status)
       IF (status/=0) CALL error_bi(' ','SCAV INIT')
       IF (coeff_para.EQ.0) THEN
         lscav_aer=.FALSE.
         WRITE(*,*) 'No scav_parameter chosen for aerosol scavenging !!!!' 
         WRITE(*,*) 'Aerosol scavenging switched off'
       ENDIF
    END IF

    ! Solver Method for KPP Liquid phase
    LICNTRL(3) = 2               ! =ROS3
!    LICNTRL(3) = 4               ! =RODAS3
    ! Solver Method for KPP Ice phase
    IICNTRL(3) = 2               ! =ROS3
    


    CALL p_bcast(lscav, p_io)
    IF (lscav) CALL info_bi( 'Scavenging active', substr)

    CALL p_bcast(lscav_ls, p_io)
    IF (lscav_ls) &
      CALL info_bi( 'Scavenging active for large scale precipitation', substr)
    CALL p_bcast(lscav_cv, p_io)
    IF (lscav_cv) &
      CALL info_bi( 'Scavenging active for convective precipitation', substr)
    CALL p_bcast(lscav_nuc, p_io)
    IF (lscav_nuc) CALL info_bi( 'Nucleation - Scavenging active', substr)
    CALL p_bcast(lscav_imp, p_io)
    IF (lscav_imp) CALL info_bi( 'Impaction - Scavenging active', substr)
    CALL p_bcast(lscav_l, p_io)
    IF (lscav_l) &
      CALL info_bi( 'Scavenging active for liquid water precipitation', substr)
    CALL p_bcast(lscav_i, p_io)
    IF (lscav_i) &
      CALL info_bi( 'Scavenging active for snow / ice precipitation', substr)

    CALL p_bcast(lscav_gas, p_io)
    IF (lscav_gas) &
      CALL info_bi( 'Scavenging active for gasphase species', substr)
    CALL p_bcast(lscav_aer, p_io)
    IF (lscav_aer) &
      CALL info_bi( 'Scavenging active for aerosol species', substr)
    CALL p_bcast(lscav_easy, p_io)
    IF (lscav_easy) &
      CALL info_bi( 'Easy - Scavenging for gas / liquid  phase active', substr)
    CALL p_bcast(l_scav_easy, p_io)
    IF (lscav_easy) THEN
      IF (l_scav_easy == 1) &
        CALL info_bi('Easy-Scavenging with fixed coefficients active', substr)
      IF (l_scav_easy == 2) &
        CALL info_bi('Easy-Scavenging with absolute Henry`s law', substr)
      IF (l_scav_easy == 3) &
        CALL info_bi('Easy-Scavenging with effective Henry`s law', substr)
    ENDIF
    CALL p_bcast(iscav_easy, p_io)
    IF (iscav_easy) &
      CALL info_bi( 'Easy - Scavenging for gas / ice  phase active', substr)
    CALL p_bcast(i_scav_easy, p_io)
    IF (iscav_easy) THEN
      IF (i_scav_easy == 1) &
        CALL info_bi(&
        'Easy-Scavenging with fixed coefficients for ice active', substr)
      IF (i_scav_easy == 2)  &
        CALL info_bi(&
        'Easy-Scavenging with pseudo Henry coefficients for ice active', substr)
      IF (i_scav_easy == 3)  &
        CALL info_bi(&
        'Easy-Scavenging with iterative Langmuir uptake for ice active', substr)
      if (i_scav_easy == 4)  &
        CALL info_bi (&
        'Easy-Scavenging with iterative HNO3 Trapping for ice active', substr)
    ENDIF
    CALL p_bcast(loverwrite_henry, p_io)
    IF (loverwrite_henry) CALL warning_bi(&
         'Henry coefficients from CHEMPROP will be overwritten!!', substr)
    CALL p_bcast(loverwrite_alpha, p_io)
    IF (loverwrite_alpha) THEN
       CALL warning_bi(&
            'alpha coefficients from CHEMPROP will be overwritten!!', substr)
       luse_empiric_alpha = .TRUE.
    END IF
    CALL p_bcast(luse_empiric_alpha, p_io)
    IF (luse_empiric_alpha) THEN
       CALL info_bi(&
         'use more complex alpha calculation for some species!', substr)
    ENDIF
       
    CALL p_bcast(use_schwartz, p_io)

    CALL p_bcast(coeff_para, p_io)
    IF (lscav_aer) THEN
      IF (i_evap == 1) CALL info_bi( &
        'After evaporation aerosols will be redistributed into the mode where they originate from.', substr)
      IF (i_evap == 2) CALL info_bi( &
        'After evaporation aerosols will be redistributed into the largest available soluble mode of the aerosol module.', substr)
    ENDIF
    CALL p_bcast(i_evap, p_io)
    CALL p_bcast(cpl_aerosol, p_io)
    CALL p_bcast(frac_resnum, p_io)

    CALL p_bcast(iscav_rate_hiTmp, p_io)
    CALL p_bcast(iscav_rate_loTmp_het, p_io)
    CALL p_bcast(iscav_rate_loTmp_hom, p_io)

    IF (.NOT.lscav) RETURN

    ! INITIALIZE COUPLING-CONTROL
    IF (p_parallel_io) THEN
       iou = find_next_free_unit(100,200)
       CALL scav_read_nml_cpl(status, iou)
       IF (status /= 0) CALL error_bi(' ', 'SCAV COUPLING INIT')
    END IF
   
    CALL p_bcast(rain_cv(1), p_io)
    CALL p_bcast(rain_cv(2), p_io)
    CALL p_bcast(snow_cv(1), p_io)
    CALL p_bcast(snow_cv(2), p_io)
    CALL p_bcast(ccover(1), p_io)  
    CALL p_bcast(ccover(2), p_io)  
    CALL p_bcast(cvprec(1), p_io) 
    CALL p_bcast(cvprec(2), p_io) 
    CALL p_bcast(cvsnow(1), p_io) 
    CALL p_bcast(cvsnow(2), p_io) 
    CALL p_bcast(lcover(1), p_io)    
    CALL p_bcast(lcover(2), p_io)    
    CALL p_bcast(mass(1), p_io)  
    CALL p_bcast(mass(2), p_io)  
    CALL p_bcast(vol(1), p_io)  
    CALL p_bcast(vol(2), p_io)  

    CALL p_bcast(press(1), p_io)  
    CALL p_bcast(press(2), p_io)  
    CALL p_bcast(pressi(1), p_io)  
    CALL p_bcast(pressi(2), p_io) 
    CALL p_bcast(boxarea(1), p_io)  
    CALL p_bcast(boxarea(2), p_io) 
    CALL p_bcast(ratep(1), p_io)
    CALL p_bcast(ratep(2), p_io)
    CALL p_bcast(prec(1), p_io)  
    CALL p_bcast(prec(2), p_io)
    CALL p_bcast(ratesi(1), p_io)   
    CALL p_bcast(ratesi(2), p_io)   
    CALL p_bcast(fsi(1), p_io)    
    CALL p_bcast(fsi(2), p_io)      
    CALL p_bcast(rcover(1), p_io)  
    CALL p_bcast(rcover(2), p_io)        
    CALL p_bcast(lwc(1), p_io) 
    CALL p_bcast(lwc(2), p_io) 
    CALL p_bcast(iwc(1), p_io) 
    CALL p_bcast(iwc(2), p_io) 
    CALL p_bcast(imelt(1), p_io) 
    CALL p_bcast(imelt(2), p_io) 
    CALL p_bcast(isedi(1), p_io) 
    CALL p_bcast(isedi(2), p_io) 
    CALL p_bcast(photol(1), p_io)
    CALL p_bcast(photol(2), p_io)
    CALL p_bcast(out_string, p_io)  
    CALL p_bcast(out_string_aer, p_io)
    CALL p_bcast(aermod_gp, p_io)
    IF (TRIM(aermod_gp) == 'made' .OR. TRIM(aermod_gp) == 'MADE') THEN
      n_resmod = 1
    END IF
    CALL p_bcast(aermod_lg, p_io)

    CALL p_bcast(made3params%ks, p_io)
    CALL p_bcast(made3params%km, p_io)
    CALL p_bcast(made3params%ki, p_io)
    CALL p_bcast(made3params%as, p_io)
    CALL p_bcast(made3params%am, p_io)
    CALL p_bcast(made3params%ai, p_io)
    CALL p_bcast(made3params%cs, p_io)
    CALL p_bcast(made3params%cm, p_io)
    CALL p_bcast(made3params%ci, p_io)
    CALL p_bcast(made3params%nmod, p_io)
    CALL p_bcast(made3params%i_mode, p_io)
    CALL p_bcast(made3params%csubname, p_io)
    CALL p_bcast(made3params%evap_mode, p_io)
    CALL p_bcast(made3params%n_resmod, p_io)
    IF (TRIM(aermod_gp) == 'made3' .OR. TRIM(aermod_gp) == 'MADE3') &
         n_resmod = made3params%n_resmod

    CALL p_bcast(liqfrac(:), p_io)
    CALL p_bcast(altitude, p_io)
    CALL p_bcast(thres_val, p_io)

    CALL p_bcast(cvlwc(1),   p_io)
    CALL p_bcast(cvlwc(2),   p_io)
    CALL p_bcast(cvrform(1), p_io)
    CALL p_bcast(cvrform(2), p_io)
    CALL p_bcast(cviwc(1),   p_io)
    CALL p_bcast(cviwc(2),   p_io)
    CALL p_bcast(cvsform(1), p_io)
    CALL p_bcast(cvsform(2), p_io)
    CALL p_bcast(nfrac_nuc(1), p_io)
    CALL p_bcast(nfrac_nuc(2), p_io)
    CALL p_bcast(mfrac_nuc(1), p_io)
    CALL p_bcast(mfrac_nuc(2), p_io)

    CALL p_bcast(lscav_gp, p_io)
    IF (lscav_gp) &
      CALL info_bi( 'Scavenging active for gridpoint tracers', substr)

    CALL p_bcast(te_string, p_io)

    CALL p_bcast(lscav_lg, p_io)
!!#D attila +
#ifdef ECHAM5
    IF (lscav_lg) &
      CALL info_bi( 'Scavenging active for Lagrangian tracers', substr)
#endif
!!#D attila -
    CALL info_bi('Setting up internal species structure of SCAV', substr)
    CALL INIT_SCAV_SPEC_STRUCT

#ifdef MESSYTENDENCY
    tend_cv = mtend_get_handle(modstr//'_cv')
    tend_ls = mtend_get_handle(modstr//'_ls')
    tend_ev = mtend_get_handle(modstr//'_ev')
#endif

  END SUBROUTINE scav_initialize
! ==============================================================================

  SUBROUTINE scav_read_nml_cpl(status, iou)

    ! SCAV MODULE ROUTINE (ECHAM-5 INTERFACE, PRIVATE)
    !
    ! read namelist for 'coupling' to ECHAM5
    !
    ! Author: H. Tost, MPICH, March 2004

    ! MESSy
    USE messy_main_tools, ONLY: read_nml_open, read_nml_check, read_nml_close
#ifdef ECHAM5
    USE messy_main_tracer_mem_bi, ONLY: NGCELL
#endif
    IMPLICIT NONE

    ! I/O
    INTEGER, INTENT(OUT) :: status     ! error status
    INTEGER, INTENT(IN)  :: iou        ! I/O unit

    NAMELIST /CPL/ rain_cv, ccover, cvprec, cvsnow, snow_cv,              &
                   lcover, rcover, ratep, prec, ratesi, fsi, lwc, iwc,    &
                   cvlwc, cviwc, cvrform, cvsform, imelt, isedi,          &
                   mass, vol, press, pressi, photol,                      &
                   boxarea, &
                   nfrac_nuc, mfrac_nuc,                                  &
                   out_string, out_string_aer,                            &
                   LIQFRAC, thres_val, altitude,                          &
                   lscav_lg, lscav_gp,                                    &
                   aermod_gp, aermod_lg,                                  &
                   made3params,                                           &
                   te_string
    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr='scav_read_nml_cpl'
    LOGICAL              :: lex      ! file exists ?
    INTEGER              :: fstat    ! file status

    status = 1

    ! INITIALIZE NAMELIST VARIABLES
 
    rain_cv(:) = ''
    snow_cv(:) = ''
    ccover(:)  = ''
    cvprec(:)  = ''
    cvsnow(:)  = ''
    lcover(:)  = ''    
    rcover(:)  = ''
    ratep(:)   = ''
    prec(:)    = ''          
    ratesi(:)  = '' 
    fsi(:)     = ''
    lwc(:)     = ''
    iwc(:)     = ''
    imelt(:)   = ''
    isedi(:)   = ''
    mass(:)    = ''
    vol(:)     = ''
    press(:)   = ''
    pressi(:)  = ''
    boxarea(:) = ''
    photol(:)  = ''
    out_string = ''
    out_string_aer = ''
    aermod_gp  = ''
    aermod_lg  = ''
    mfrac_nuc(:) = ''
    nfrac_nuc(:) = ''

    liqfrac(:) = (/ 0.0_dp, 273.15_dp /)
    thres_val  = 1.e-7_dp
    altitude   = 1._dp

    CALL read_nml_open(lex, substr, iou, 'CPL', modstr)
    IF (.NOT.lex) RETURN    ! <modstr>.nml does not exist
    READ(iou, NML=CPL, IOSTAT=fstat)
    CALL read_nml_check(fstat, substr, iou, 'CPL', modstr)
    IF (fstat /= 0) RETURN  ! error while reading namelist
  
    CALL read_nml_close(substr, iou, modstr)

#ifdef ECHAM5
    ! CHECK NAMELIST
    IF (lscav_lg) THEN
       IF (NGCELL > 0) THEN
!!#D attila +
          WRITE(*,*) 'LAGRANGIAN CALCULATION: ON'
!!#D attila -
       ELSE
          IF (lscav_lg) THEN
!!#D attila +
             WRITE(*,*) 'lscav_lg = T in namelist'
             WRITE(*,*) 'However no Lagrangian scheme activated ...'
             WRITE(*,*) ' ... setting lscav_lg = F'
!!#D attila -
          END IF
          lscav_lg = .FALSE.
       END IF
    ELSE
!!#D attila +
       WRITE(*,*) 'LAGRANGIAN CALCULATION: OFF'
!!#D attila -
    END IF
#endif
    ! LIMIT LIQFRAC
    LIQFRAC(1) = MAX(0.0_dp, LIQFRAC(1))
    LIQFRAC(1) = MIN(1.0_dp, LIQFRAC(1))
    thres_val  = REAL(THRES_VAL,DP)

    IF ((TRIM(aermod_gp) == 'made3' .OR. TRIM(aermod_gp) == 'MADE3') .AND. &
         (made3params%ks .EQ. 0 .OR. &
         made3params%km .EQ. 0 .OR. &
         made3params%ki .EQ. 0 .OR. &
         made3params%as .EQ. 0 .OR. &
         made3params%am .EQ. 0 .OR. &
         made3params%ai .EQ. 0 .OR. &
         made3params%cs .EQ. 0 .OR. &
         made3params%cm .EQ. 0 .OR. &
         made3params%ci .EQ. 0 .OR. &
         made3params%nmod .EQ. 0 .OR. &
         made3params%i_mode .EQ. 0 .OR. &
         made3params%csubname == '' .OR. &
         made3params%evap_mode .EQ. 0 .OR. &
         made3params%n_resmod .EQ. 0)) THEN
      status = 3
      RETURN
    END IF

    status = 0 ! NO ERROR

  END SUBROUTINE scav_read_nml_cpl
! ==============================================================================
  SUBROUTINE scav_new_tracer

    ! MESSy
    USE messy_main_tracer,          ONLY: AIR, ON, AEROSOL,                  &
                                          AMOUNTFRACTION, OFF, CLOUD, MODAL, &
                                          I_ADVECT, I_CONVECT,               &
                                          I_VDIFF,                           &
                                          I_DRYDEP, I_SEDI,                  &
                                          I_SCAV, I_MIX, I_CHARGE,           &
                                          I_AEROSOL_METHOD, I_AEROSOL_MODE,  &
                                          I_AEROSOL_SOL, S_AEROSOL_MODEL,    &
                                          R_MOLARMASS, R_AEROSOL_DENSITY,    &
                                          R_pss,      R_dryreac_sf,          &
                                          new_tracer, set_tracer, get_tracer,&
                                          get_chemprop
    USE messy_main_tracer_tools_bi, ONLY: tracer_halt 
    USE messy_main_blather_bi,      ONLY: error_bi
    USE messy_main_tools,           ONLY: strcrack
    USE messy_main_constants_mem,   ONLY: MH, MN, MO, MS, MCl
    USE messy_scav_l_kpp,           ONLY: str_field_kpp_l => spc_names
    USE messy_scav_i_kpp,           ONLY: str_field_kpp_i => spc_names

    IMPLICIT NONE
    INTEGER :: status


    INTEGER :: i, idt_dummy, dummy, jt, idx, idt
    INTEGER :: i_mode, status_l, status_i
    LOGICAL :: tracer_found

    ! number of aerosol species to scan through when looking for a
    ! specific component

    INTEGER, PARAMETER :: num_NO3 = 6
    INTEGER, PARAMETER :: num_NH4 = 6
    INTEGER, PARAMETER :: num_SO4 = 8
    INTEGER, PARAMETER :: num_Cl  = 6
    INTEGER, PARAMETER :: num_Hp  = 4
    INTEGER            :: js

    CHARACTER(LEN=*), DIMENSION(2), PARAMETER :: setname = (/GPTRSTR, LGTRSTR/)
    CHARACTER(LEN=STRLEN_MEDIUM)           :: aermod =''
    CHARACTER(LEN=STRLEN_MEDIUM), POINTER  :: strname(:) => NULL()
    CHARACTER(LEN=STRLEN_MEDIUM)           :: c_name, f_name
    CHARACTER(LEN=2)            :: csubname   = ''
    TYPE(t_trinfo)              :: ti_init
    CHARACTER(LEN=*), PARAMETER :: substr='scav_new_tracer'

    CHARACTER (LEN=8) NO3str1(num_NO3), NO3str2(num_NO3)
    CHARACTER (LEN=8) NH4str1(num_NH4), NH4str2(num_NH4)
    CHARACTER (LEN=8) SO4str1(num_SO4), SO4str2(num_SO4)
    CHARACTER (LEN=8) Clstr1(num_Cl),   Clstr2(num_Cl)
    CHARACTER (LEN=8) Hpstr1(num_Hp),   Hpstr2(num_Hp)

    INTEGER :: iicharge
    INTEGER :: istat

    ! names of aerosol species to scan through when looking for a
    ! specific component

    ! nitrate, NO3-
    DATA NO3str1 / 'NO3m', 'NO3m', 'NO3', 'NO3', 'NO3m', 'NO3m' /
    DATA NO3str2 / 'cs',   'as',   'cm',  'am',  'a01',  'a02'  /

    ! ammonium, NH4+
    DATA NH4str1 / 'NH4p', 'NH4p', 'NH4p', 'NH4p', 'NH4', 'NH4' /
    DATA NH4str2 / 'cs',   'as',   'a01',  'a02',  'cm',  'am'  /

    ! sulfate, SO4--
    DATA SO4str1 / 'SO4', 'SO4', 'SO4', 'SO4', 'SO4', 'SO4', 'SO4mm', 'SO4mm' /
    DATA SO4str2 / 'cs' , 'as' , 'a01', 'a02', 'cm' , 'am',  'cs'   , 'as'    /

    ! chloride, Cl-
    DATA Clstr1  / 'Clm', 'Clm', 'Clm', 'Clm', 'Cl', 'Cl' /
    DATA Clstr2  / 'cs',  'as',  'a01', 'a02', 'cm', 'am' /

    ! Hplus, H+
    DATA Hpstr1 / 'Hp', 'Hp',  'Hp', 'Hp'/
    DATA Hpstr2 / 'cs', 'as', 'a01', 'a02' /

    !
    ! request tracers from ECHAM
    !
 
    IF (.NOT.lscav) RETURN
    IF (lscav_easy) RETURN
    IF (.NOT.lscav_gas) RETURN
    
    DO js=1,2

       IF (js == 1) THEN
          IF (.NOT.lscav_gp) THEN
            CYCLE
          ELSE
            aermod = aermod_gp
            ntrac  = ntrac_gp
          ENDIF       
       END IF

       IF (js == 2) THEN
          IF (.NOT.lscav_lg) THEN
            CYCLE
          ELSE
            aermod = aermod_lg
            ntrac  = ntrac_lg
          ENDIF       
       END IF
       
       ! build gas phase tracers if they are required by the SCAV mechanism
       ! either they have emission source, or they are produced in the 
       ! aqueous/ice phase and exchanged with the gas phase
       ! liquid
       DO jt = 1,lspec_gas
         idx = kpp0_l_idx%gas_spec(jt,gas_idx)
         CALL get_tracer(status, setname(js), TRIM(str_field_kpp_l(idx)), &
                         idx = idt)
         IF (status /=0) THEN
            CALL new_tracer(status, setname(js), TRIM(str_field_kpp_l(idx)),  &
                 modstr, quantity = AMOUNTFRACTION, unit = 'mol/mol',         &
                 medium = AIR, idx = idt)
            CALL tracer_halt(substr,status)
            CALL set_tracer(status, setname(js), idt, I_advect          , ON)
            CALL tracer_halt(substr, status)
            CALL set_tracer(status, setname(js), idt, I_convect         , ON)
            CALL tracer_halt(substr, status)
            CALL set_tracer(status, setname(js), idt, I_vdiff           , ON)
            CALL tracer_halt(substr, status)
            CALL set_tracer(status, setname(js), idt, I_scav            , ON)
            CALL tracer_halt(substr, status)
            CALL set_tracer(status, setname(js), idt, I_drydep          , ON)
            CALL tracer_halt(substr, status)
            CALL set_tracer(status, setname(js), idt, I_sedi            , ON)
            CALL tracer_halt(substr, status)
            CALL set_tracer(status, setname(js), idt, I_mix             , OFF)
            CALL tracer_halt(substr, status)
         END IF
       END DO
      ! ice
      DO jt = 1,ispec_gas
         idx = kpp0_i_idx%gas_spec(jt,gas_idx)
         CALL get_tracer(status, setname(js), TRIM(str_field_kpp_i(idx)), &
              idx=dummy)
         IF (status /=0) THEN
           CALL new_tracer(status, setname(js), TRIM(str_field_kpp_i(idx)),  &
            modstr, quantity = AMOUNTFRACTION, unit = 'mol/mol',         &
            medium = AIR, idx = idt)
           CALL tracer_halt(substr,status)
           CALL set_tracer(status, setname(js), idt, I_advect          , ON)
           CALL tracer_halt(substr, status)
           CALL set_tracer(status, setname(js), idt, I_convect         , ON)
           CALL tracer_halt(substr, status)
           CALL set_tracer(status, setname(js), idt, I_vdiff           , ON)
           CALL tracer_halt(substr, status)
           CALL set_tracer(status, setname(js), idt, I_scav            , ON)
           CALL tracer_halt(substr, status)
           CALL set_tracer(status, setname(js), idt, I_drydep          , ON)
           CALL tracer_halt(substr, status)
           CALL set_tracer(status, setname(js), idt, I_sedi            , ON)
           CALL tracer_halt(substr, status)
           CALL set_tracer(status, setname(js), idt, I_mix             , OFF)
           CALL tracer_halt(substr, status)

         END IF
       END DO

       SELECT CASE (TRIM(aermod))
       CASE('made','MADE')
         i_mode = 2   ! accumulation mode, 
                      ! since this is chemically active
         csubname = 'am'
       CASE('made3','MADE3')
         i_mode   = made3params%i_mode
         csubname = made3params%csubname
       CASE DEFAULT
         i_mode = 4   ! coarse, soluble (cs)
         csubname = 'cs'
       END SELECT


       IF (l_trac_aer) THEN
         ! build aerosol tracers for all liquid compounds enforcing 
         ! a 100% mass conservation after evaporation
         DO jt = 1,lspec_liq
           idx = kpp0_l_idx%liq_spec(jt,liq_idx)
           IF (ASSOCIATED(strname)) DEALLOCATE (strname)
           NULLIFY(strname) 
           CALL strcrack(STR_FIELD_KPP_L(idx),'_', strname, dummy)
           IF (TRIM(strname(1)) == "Prod") CYCLE
           c_name=TRIM(strname(1))
           CALL get_tracer(status, setname(js), TRIM(c_name), &
             subname = csubname, idx = dummy)
           IF (status==0 ) THEN
              CALL get_tracer(istat, setname(js),dummy, I_CHARGE, iicharge)
              IF (istat /= 0) iicharge  = 0
           ELSE
              istat = get_chemprop(TRIM(c_name), I_CHARGE, iicharge)
              IF (istat /=0)  CALL error_bi (&
                   'species '//TRIM(c_name)//' not part of chemprop', substr)
           END IF
           IF ( (.NOT. l_trac_aer_all) .AND. (iicharge == 0) ) CYCLE
           IF (status /= 0) THEN
             CALL new_tracer(status, setname(js), TRIM(c_name),  &
               modstr, subname=csubname, quantity = AMOUNTFRACTION,      &
               unit = 'mol/mol', medium = AEROSOL, idx = idt)
          
             CALL tracer_halt(substr,status)
             CALL set_tracer(status, setname(js), idt, I_advect          , ON)
             CALL tracer_halt(substr, status)
             CALL set_tracer(status, setname(js), idt, I_convect         , ON)
             CALL tracer_halt(substr, status)
             CALL set_tracer(status, setname(js), idt, I_vdiff           , ON)
             CALL tracer_halt(substr, status)
             CALL set_tracer(status, setname(js), idt, I_scav            , ON)
             CALL tracer_halt(substr, status)
             CALL set_tracer(status, setname(js), idt, I_drydep          , ON)
             CALL tracer_halt(substr, status)
             CALL set_tracer(status, setname(js), idt, I_sedi            , ON)
             CALL tracer_halt(substr, status)
             CALL set_tracer(status, setname(js), idt, I_mix             , OFF)
             CALL tracer_halt(substr, status)
             CALL set_tracer(status, setname(js), idt, I_aerosol_mode    , &
               i_mode)
             CALL tracer_halt(substr, status)
             CALL set_tracer(status, setname(js), idt, S_aerosol_model   , &
               TRIM(aermod) )
             CALL tracer_halt(substr, status)
             CALL set_tracer(status, setname(js), idt, I_aerosol_method , MODAL)
             CALL tracer_halt(substr, status)
             CALL set_tracer(status, setname(js), idt, I_aerosol_sol     , 1)
             CALL tracer_halt(substr, status)
           END IF
         END DO
         DO jt = 1,ispec_ice
           idx = kpp0_i_idx%ice_spec(jt,ice_idx)
           IF (ASSOCIATED(strname)) DEALLOCATE (strname)
           NULLIFY(strname) 
           CALL strcrack(STR_FIELD_KPP_I(idx),'_', strname, dummy)
           IF (TRIM(strname(1)) == "Prod") CYCLE
           c_name=TRIM(strname(1))
           CALL get_tracer(status, setname(js), TRIM(c_name), &
             subname = csubname, idx = dummy)
           IF (status == 0) THEN
              CALL get_tracer(istat, setname(js), dummy, I_CHARGE, iicharge)
              IF (istat /= 0) iicharge = 0
           ELSE
              istat = get_chemprop(TRIM(c_name), I_CHARGE, iicharge)
              IF (istat /=0)  CALL error_bi (&
                   'species '//TRIM(c_name)//' not part of chemprop', substr)
           END IF

           IF ( (.NOT. l_trac_aer_all) .AND. (iicharge == 0)) CYCLE
           IF (status /= 0) THEN
             CALL new_tracer(status, setname(js), TRIM(c_name),  &
               modstr, subname=csubname, quantity = AMOUNTFRACTION,      &
               unit = 'mol/mol', medium = AEROSOL, idx = idt)
          
             CALL tracer_halt(substr,status)
             CALL set_tracer(status, setname(js), idt, I_advect         , ON)
             CALL tracer_halt(substr, status)
             CALL set_tracer(status, setname(js), idt, I_convect        , ON)
             CALL tracer_halt(substr, status)
             CALL set_tracer(status, setname(js), idt, I_vdiff          , ON)
             CALL tracer_halt(substr, status)
             CALL set_tracer(status, setname(js), idt, I_scav           , ON)
             CALL tracer_halt(substr, status)
             CALL set_tracer(status, setname(js), idt, I_drydep         , ON)
             CALL tracer_halt(substr, status)
             CALL set_tracer(status, setname(js), idt, I_sedi           , ON)
             CALL tracer_halt(substr, status)
             CALL set_tracer(status, setname(js), idt, I_mix            , OFF)
             CALL tracer_halt(substr, status)
             CALL set_tracer(status, setname(js), idt, I_aerosol_mode   , &
               i_mode)
             CALL tracer_halt(substr, status)
             CALL set_tracer(status, setname(js), idt, S_aerosol_model  , &
               TRIM(aermod) )
             CALL tracer_halt(substr, status)
             CALL set_tracer(status, setname(js), idt, I_aerosol_method , MODAL)
             CALL tracer_halt(substr, status)
             CALL set_tracer(status, setname(js), idt, I_aerosol_sol    , 1)
             CALL tracer_halt(substr, status)
           END IF
         END DO
         
       END IF ! l_trac_aer

       IF (l_trac_cloud) THEN
         ! build cloud tracers for all liquid/ice compounds 
         ! for transport of in-cloud components and no artificial evaporation 
         ! at the end of each time step

         DO jt = 1,lspec_liq
           idx = kpp0_l_idx%liq_spec(jt,liq_idx)
           IF (ASSOCIATED(strname)) DEALLOCATE (strname)
           NULLIFY(strname) 
           CALL strcrack(STR_FIELD_KPP_L(idx),'_', strname, dummy)
           IF (TRIM(strname(1)) == "Prod") CYCLE
           IF (TRIM(strname(1)) == "H2O") CYCLE
           c_name=TRIM(strname(1))
           CALL get_tracer(status, setname(js), TRIM(c_name), &
             subname = 'l', idx = dummy)
           IF (status /= 0) THEN
             CALL new_tracer(status, setname(js), TRIM(c_name),  &
               modstr, subname='l', quantity = AMOUNTFRACTION,      &
               unit = 'mol/mol', medium = CLOUD, idx = idt)
          
             CALL tracer_halt(substr,status)
             CALL set_tracer(status, setname(js), idt, I_advect          , ON)
             CALL tracer_halt(substr, status)
             CALL set_tracer(status, setname(js), idt, I_convect         , ON)
             CALL tracer_halt(substr, status)
             CALL set_tracer(status, setname(js), idt, I_vdiff           , ON)
             CALL tracer_halt(substr, status)
             CALL set_tracer(status, setname(js), idt, I_scav            , ON)
             CALL tracer_halt(substr, status)
             CALL set_tracer(status, setname(js), idt, I_drydep          , OFF)
             CALL tracer_halt(substr, status)
             CALL set_tracer(status, setname(js), idt, I_sedi            , OFF)
             CALL tracer_halt(substr, status)
             CALL set_tracer(status, setname(js), idt, I_mix             , OFF)
             CALL tracer_halt(substr, status)
           END IF
         END DO
         DO jt = 1,ispec_ice
           idx = kpp0_i_idx%ice_spec(jt,ice_idx)
           IF (ASSOCIATED(strname)) DEALLOCATE (strname)
           NULLIFY(strname) 
           CALL strcrack(STR_FIELD_KPP_I(idx),'_', strname, dummy)
           IF (TRIM(strname(1)) == "Prod") CYCLE
           c_name=TRIM(strname(1))
           CALL get_tracer(status, setname(js), TRIM(c_name), &
             subname = 'i', idx = dummy)
           IF (status /= 0) THEN
             CALL new_tracer(status, setname(js), TRIM(c_name),  &
               modstr, subname='i', quantity = AMOUNTFRACTION,      &
               unit = 'mol/mol', medium = CLOUD, idx = idt)
          
             CALL tracer_halt(substr,status)
             CALL set_tracer(status, setname(js), idt, I_advect          , ON)
             CALL tracer_halt(substr, status)
             CALL set_tracer(status, setname(js), idt, I_convect         , ON)
             CALL tracer_halt(substr, status)
             CALL set_tracer(status, setname(js), idt, I_vdiff           , ON)
             CALL tracer_halt(substr, status)
             CALL set_tracer(status, setname(js), idt, I_scav            , ON)
             CALL tracer_halt(substr, status)
             CALL set_tracer(status, setname(js), idt, I_drydep          , OFF)
             CALL tracer_halt(substr, status)
             CALL set_tracer(status, setname(js), idt, I_sedi            , OFF)
             CALL tracer_halt(substr, status)
             CALL set_tracer(status, setname(js), idt, I_mix             , OFF)
             CALL tracer_halt(substr, status)
           END IF
         END DO
         
       END IF ! l_trac_cloud

        CALL get_tracer(status, setname(js), 'NH3',idx=idt_scav_NH3)
        IF (status /=0) THEN
             CALL new_tracer(status, setname(js), 'NH3', 'scav',             &
                          idx=idt_scav_NH3, unit='mol/mol',                  &
                          medium = AIR, quantity=amountfraction   )
             CALL tracer_halt(substr,status)

             CALL set_tracer(status, setname(js), idt_scav_NH3 &
                  , R_molarmass,  MN+3.*MH)
             CALL tracer_halt(substr, status)
             CALL set_tracer(status, setname(js), idt_scav_NH3 &
                  , R_pss,        2.e4_dp)
             CALL tracer_halt(substr, status)
             CALL set_tracer(status, setname(js), idt_scav_NH3 &
                  , R_dryreac_sf, 1.0_dp)
             CALL tracer_halt(substr, status)
             CALL set_tracer(status, setname(js), idt_scav_NH3 &
                  , I_scav      , ON)
             CALL tracer_halt(substr, status)
             CALL set_tracer(status, setname(js), idt_scav_NH3 &
                  , I_drydep    , ON)
             CALL tracer_halt(substr, status)
          END IF
     IF (cpl_aerosol.GT.0) THEN

       ! Check wether dummy tracer need to be created by SCAV

        SELECT CASE (TRIM(aermod))
           CASE('made','MADE')
              i_mode = 2   ! accumulation mode, 
                           ! since this is chemically active
              csubname = 'am'
           CASE('made3','MADE3')
              i_mode   = made3params%i_mode
              csubname = made3params%csubname
           CASE DEFAULT
              i_mode = 4   ! coarse, soluble (cs)
              csubname = 'cs'
        END SELECT

       !nitrate
       tracer_found = .FALSE.
       DO i = 1, num_NO3
          ! Look for aerosol NO3 tracer (coarse or accumulation mode)...
          CALL get_tracer(status, setname(js), TRIM(NO3str1(i)), &
                          subname=TRIM(NO3str2(i)), idx=idt_dummy)
          IF (status == 0) THEN
             tracer_found = .TRUE.
             EXIT
          END IF
       END DO
       IF (.NOT.tracer_found) THEN  ! create dummy tracer for NO3
         CALL new_tracer(status, setname(js), 'NO3mres', 'scav',       &
           subname=csubname, unit='mol/mol',                           &
           medium = AEROSOL, quantity = AMOUNTFRACTION,                &
           idx=idt_scav_NO3mres_cs)
         CALL tracer_halt(substr,status)

         CALL set_tracer(status, setname(js), idt_scav_NO3mres_cs &
              , R_molarmass,  MN+3.*MO)
         CALL tracer_halt(substr, status)
         CALL set_tracer(status, setname(js), idt_scav_NO3mres_cs &
              , I_scav      , ON)
         CALL tracer_halt(substr, status)
         CALL set_tracer(status, setname(js), idt_scav_NO3mres_cs &
              , I_drydep    , ON)
         CALL tracer_halt(substr, status)
         CALL set_tracer(status, setname(js), idt_scav_NO3mres_cs &
              , I_sedi    , ON)
         CALL tracer_halt(substr, status)
         CALL set_tracer(status, setname(js), idt_scav_NO3mres_cs &
              , S_aerosol_model    , TRIM(aermod))
         CALL tracer_halt(substr, status)
         CALL set_tracer(status, setname(js), idt_scav_NO3mres_cs &
              , I_aerosol_mode    , i_mode)
         CALL tracer_halt(substr, status)
         CALL set_tracer(status, setname(js), idt_scav_NO3mres_cs &
              , R_aerosol_density    , 1513._dp)
         CALL tracer_halt(substr, status)
       END IF

       !ammonium
       tracer_found = .FALSE.
       DO i = 1, num_NH4
          ! Look for aerosol NH4 tracer (coarse or accumulation mode)...
          CALL get_tracer(status, setname(js), TRIM(NH4str1(i)),           &
                          subname=TRIM(NH4str2(i)), idx=idt_dummy)
          IF (status == 0) THEN
             tracer_found = .TRUE.
             EXIT
          END IF
       END DO
       IF (.NOT.tracer_found) THEN  ! create dummy tracer for NH4
         CALL new_tracer(status, setname(js), 'NH4pres', 'scav',       &
               subname=csubname,                                               &
               idx=idt_scav_NH4pres_cs, unit='mol/mol',  &
               Quantity = Amountfraction, medium = AEROSOL)
          CALL tracer_halt(substr,status)
          
          CALL set_tracer(status, setname(js), idt_scav_NH4pres_cs &
               , R_molarmass,  MN+4.*MH)
          CALL tracer_halt(substr, status)
          CALL set_tracer(status, setname(js), idt_scav_NH4pres_cs &
               , I_scav      , ON)
          CALL tracer_halt(substr, status)
          CALL set_tracer(status, setname(js), idt_scav_NH4pres_cs &
               , I_drydep    , ON)
          CALL tracer_halt(substr, status)
          CALL set_tracer(status, setname(js), idt_scav_NH4pres_cs &
               , I_sedi    , ON)
          CALL tracer_halt(substr, status)
          CALL set_tracer(status, setname(js), idt_scav_NH4pres_cs &
               , S_aerosol_model    , TRIM(aermod))
          CALL tracer_halt(substr, status)
          CALL set_tracer(status, setname(js), idt_scav_NH4pres_cs &
               , I_aerosol_mode    , i_mode)
          CALL tracer_halt(substr, status)
          CALL set_tracer(status, setname(js), idt_scav_NH4pres_cs &
               , R_aerosol_density    , 696.0_dp)
          CALL tracer_halt(substr, status)
       END IF

       !sulfate
       tracer_found = .FALSE.
       DO i = 1, num_SO4
          ! Look for aerosol SO4 tracer (coarse or accumulation mode)...
          CALL get_tracer(status, setname(js), TRIM(SO4str1(i)),      &
                          subname=TRIM(SO4str2(i)), idx=idt_dummy)
          IF (status == 0) THEN
             tracer_found = .TRUE.
             EXIT
          END IF
       END DO
       IF (.NOT.tracer_found) THEN  ! create dummy tracer for SO4
          CALL new_tracer(status, setname(js), 'SO4res', 'scav',            &
               subname=csubname,                                               &
               idx=idt_scav_SO4res_cs, unit='mol/mol', medium = AEROSOL,       &
               Quantity = Amountfraction)
          CALL tracer_halt(substr,status)
          
          CALL set_tracer(status, setname(js), idt_scav_SO4res_cs &
               , R_molarmass,  MS+4.*MO)
          CALL tracer_halt(substr, status)
          CALL set_tracer(status, setname(js), idt_scav_SO4res_cs &
               , I_scav      , ON)
          CALL tracer_halt(substr, status)
          CALL set_tracer(status, setname(js), idt_scav_SO4res_cs &
               , I_drydep    , ON)
          CALL tracer_halt(substr, status)
          CALL set_tracer(status, setname(js), idt_scav_SO4res_cs &
               , I_sedi    , ON)
          CALL tracer_halt(substr, status)
          CALL set_tracer(status, setname(js), idt_scav_SO4res_cs &
               , S_aerosol_model    , TRIM(aermod))
          CALL tracer_halt(substr, status)
          CALL set_tracer(status, setname(js), idt_scav_SO4res_cs &
               , I_aerosol_mode    , i_mode)
          CALL tracer_halt(substr, status)
          CALL set_tracer(status, setname(js), idt_scav_SO4res_cs &
               , R_aerosol_density    , 1841.0_dp)
          CALL tracer_halt(substr, status)
       END IF

       !chloride
       tracer_found = .FALSE.
       DO i = 1, num_Cl
          ! Look for aerosol Cl tracer (coarse or accumulation mode)...
          CALL get_tracer(status, setname(js), TRIM(Clstr1(i)),              &
                          subname=TRIM(Clstr2(i)), idx=idt_dummy)
          IF (status == 0) THEN
             tracer_found = .TRUE.
             EXIT
          END IF
       END DO

       ! Check if HCl is available. If not, no dummy Cl tracer is needed.
       CALL get_tracer(status, setname(js), 'HCl', idx=idt_dummy)

       IF ((.NOT.tracer_found).AND.(status == 0)) THEN
          ! HCl is present, but no corresponding aerocol Cl tracer.
          ! ---> Create dummy Cl tracer.
         CALL new_tracer(status, setname(js), 'Clmres', 'scav',           &
           subname=csubname,                                              &
           idx=idt_scav_Clmres_cs, unit='mol/mol',                        &
           medium = AEROSOL, Quantity = Amountfraction)
          CALL tracer_halt(substr,status)
          
          CALL set_tracer(status, setname(js), idt_scav_Clmres_cs &
               , R_molarmass,  MCl)
          CALL tracer_halt(substr, status)
          CALL set_tracer(status, setname(js), idt_scav_Clmres_cs &
               , I_scav      , ON)
          CALL tracer_halt(substr, status)
          CALL set_tracer(status, setname(js), idt_scav_Clmres_cs &
               , I_drydep    , ON)
          CALL tracer_halt(substr, status)
          CALL set_tracer(status, setname(js), idt_scav_Clmres_cs &
               , I_sedi    , ON)
          CALL tracer_halt(substr, status)
          CALL set_tracer(status, setname(js), idt_scav_Clmres_cs &
               , S_aerosol_model    , TRIM(aermod))
          CALL tracer_halt(substr, status)
          CALL set_tracer(status, setname(js), idt_scav_Clmres_cs &
               , I_aerosol_mode    , i_mode)
          CALL tracer_halt(substr, status)
          CALL set_tracer(status, setname(js), idt_scav_Clmres_cs &
               , R_aerosol_density    , 1490.0_dp)
          CALL tracer_halt(substr, status)
       END IF

       !Hplus
       tracer_found = .FALSE.
       DO i = 1, num_Hp
          ! Look for aerosol H+ tracer (coarse or accumulation mode)...
          CALL get_tracer(status, setname(js), TRIM(Hpstr1(i)),          &
                          subname=TRIM(Hpstr2(i)), idx=idt_dummy)
          IF (status == 0) THEN
             tracer_found = .TRUE.
             EXIT
          END IF
       END DO
       IF (.NOT.tracer_found) THEN  ! create dummy tracer for H+
          CALL new_tracer(status,setname(js), 'Hpres', 'scav'  &
               , subname=csubname,&
               idx=idt_scav_Hpres_cs, unit='mol/mol', medium = AEROSOL,  &
               Quantity = Amountfraction)
          CALL tracer_halt(substr,status)
          
          CALL set_tracer(status, setname(js), idt_scav_Hpres_cs &
               , R_molarmass,  MH)
          CALL tracer_halt(substr, status)
          CALL set_tracer(status, setname(js), idt_scav_Hpres_cs &
               , I_scav      , ON)
          CALL tracer_halt(substr, status)
          CALL set_tracer(status, setname(js), idt_scav_Hpres_cs &
               , I_drydep    , ON)
          CALL tracer_halt(substr, status)
          CALL set_tracer(status, setname(js), idt_scav_Hpres_cs &
               , I_sedi    , ON)
          CALL tracer_halt(substr, status)
          CALL set_tracer(status, setname(js), idt_scav_Hpres_cs &
               , S_aerosol_model    , TRIM(aermod))
          CALL tracer_halt(substr, status)
          CALL set_tracer(status, setname(js), idt_scav_Hpres_cs &
               , I_aerosol_mode    , i_mode)
          CALL tracer_halt(substr, status)
          CALL set_tracer(status, setname(js), idt_scav_Hpres_cs &
               , R_aerosol_density    , 999.97_dp)
          CALL tracer_halt(substr, status)
       END IF

     ENDIF
     jt = 0
     IF (l_trac_cloud) THEN
        DO 
          jt = jt + 1 
          CALL get_tracer(status, setname(js), jt, fullname=f_name, &
            subname=csubname, basename=c_name, trinfo=ti_init)
          IF (TRIM(f_name) == '') EXIT
          IF (ti_init%ident%medium /= AEROSOL) CYCLE
          IF ( (TRIM(c_name) == 'H2O') .OR. &
               (TRIM(c_name) == 'N'  ) ) CYCLE
          IF (status == 0) THEN

            CALL get_tracer(status_l, setname(js), TRIM(c_name), &
            subname = 'l', idx = dummy)

            IF (status_l /= 0) THEN
              CALL new_tracer(status_l, setname(js), TRIM(c_name),  &
                modstr, subname='l', quantity = AMOUNTFRACTION,      &
                unit = 'mol/mol', medium = CLOUD, idx = idt)
          
              CALL tracer_halt(substr,status_l)
              CALL set_tracer(status_l, setname(js), idt, I_advect        , ON)
              CALL tracer_halt(substr, status_l)
              CALL set_tracer(status_l, setname(js), idt, I_convect       , ON)
              CALL tracer_halt(substr, status_l)
              CALL set_tracer(status_l, setname(js), idt, I_vdiff         , ON)
              CALL tracer_halt(substr, status_L)
              CALL set_tracer(status_l, setname(js), idt, I_scav          , ON)
              CALL tracer_halt(substr, status_l)
              CALL set_tracer(status_l, setname(js), idt, I_drydep        , OFF)
              CALL tracer_halt(substr, status_l)
              CALL set_tracer(status_l, setname(js), idt, I_sedi          , OFF)
              CALL tracer_halt(substr, status_l)
              CALL set_tracer(status_l, setname(js), idt, I_mix           , OFF)
              CALL tracer_halt(substr, status_l)
              CALL set_tracer(status_l, setname(js), idt, R_molarmass     , &
                ti_init%meta%cask_r(R_molarmass))
              CALL tracer_halt(substr, status_l)
              CALL set_tracer(status_l, setname(js), idt, R_aerosol_density , &
                ti_init%meta%cask_r(R_aerosol_density))
              CALL tracer_halt(substr, status_l)
            END IF

            CALL get_tracer(status_i, setname(js), TRIM(c_name), &
            subname = 'i', idx = dummy)

            IF (status_i /= 0) THEN
              CALL new_tracer(status_i, setname(js), TRIM(c_name),  &
                modstr, subname='i', quantity = AMOUNTFRACTION,      &
                unit = 'mol/mol', medium = CLOUD, idx = idt)
          
              CALL tracer_halt(substr,status_i)
              CALL set_tracer(status_i, setname(js), idt, I_advect        , ON)
              CALL tracer_halt(substr, status_i)
              CALL set_tracer(status_i, setname(js), idt, I_convect       , ON)
              CALL tracer_halt(substr, status_i)
              CALL set_tracer(status_i, setname(js), idt, I_vdiff         , ON)
              CALL tracer_halt(substr, status_i)
              CALL set_tracer(status_i, setname(js), idt, I_scav          , ON)
              CALL tracer_halt(substr, status_i)
              CALL set_tracer(status_i, setname(js), idt, I_drydep        , OFF)
              CALL tracer_halt(substr, status_i)
              CALL set_tracer(status_i, setname(js), idt, I_sedi          , OFF)
              CALL tracer_halt(substr, status_i)
              CALL set_tracer(status_i, setname(js), idt, I_mix           , OFF)
              CALL tracer_halt(substr, status_i)
              CALL set_tracer(status_i, setname(js), idt, R_molarmass     , &
                ti_init%meta%cask_r(R_molarmass))
              CALL tracer_halt(substr, status_i)
              CALL set_tracer(status_i, setname(js), idt, R_aerosol_density , &
                ti_init%meta%cask_r(R_aerosol_density))
              CALL tracer_halt(substr, status_i)
            END IF

          ENDIF
        END DO
      ENDIF

    END DO ! tracer_sets

  END SUBROUTINE scav_new_tracer

 
!===============================================================================
  SUBROUTINE scav_init_memory

    !-----------------------------------------------------------------
    ! construct the scavenging tables
    ! all information specific to this tables is set in this subroutine
    !-----------------------------------------------------------------

    ! ECHAM5/MESSy
    USE messy_main_blather_bi,       ONLY: start_message_bi, end_message_bi &
                                         , info_bi, warning_bi
    USE messy_main_channel_error_bi, ONLY: channel_halt
    USE messy_main_channel_bi,       ONLY: GP_3D_MID &
                                         , GP_2D_HORIZONTAL, GP_3D_1LEV
    USE messy_main_grid_def_mem_bi,  ONLY: nlev
    USE messy_main_tracer_mem_bi,    ONLY: GPTRSTR
    ! MESSy
    USE messy_main_tracer,           ONLY: ON, AEROSOL, FAMILY, R_molarmass &
                                         , get_tracer
    USE messy_main_channel,          ONLY: new_channel, new_channel_object, &
                                           new_attribute
    USE messy_main_tools,            ONLY: strcrack, match_wild, int2str
#ifdef MESSYTENDENCY
    USE messy_main_tendency_bi,      ONLY: mtend_register, mtend_id_tracer
#endif

    USE messy_scav_l_kpp,            ONLY: INITVAL_L => INITIALIZE_L
    USE messy_scav_i_kpp,            ONLY: INITVAL_I => INITIALIZE_I
    USE messy_scav_ice,              ONLY: ALLOC_SCAV_VALUES_I,        &
                                           SET_EASY_SCAV_COEFFICIENTS_I
    USE messy_scav_liq,              ONLY: ALLOC_SCAV_VALUES_L,        &
                                           SET_EASY_SCAV_COEFFICIENTS, &
                                           DEFINE_Ks_for_effective_henry
    USE messy_scav_l_kpp,            ONLY: LSPEC=>NSPEC, L_NVAR => NVAR,  &
                                           l_rtol=>rtol, l_atol=>atol,    &
                                           str_field_kpp_l => spc_names
    USE messy_scav_i_kpp,            ONLY: ISPEC=>NSPEC, I_NVAR=> NVAR,   &
                                           i_rtol=>rtol, i_atol=>atol
    USE messy_scav_inter,            ONLY: idx_evap_num, process


    IMPLICIT NONE
    INTEGER :: i, jt, j, js, dummy, kprod
    CHARACTER(LEN=30) :: string
    CHARACTER(LEN=26), POINTER     :: strname(:) => NULL()
    INTEGER :: status
    REAL(dp), POINTER, DIMENSION(:,:,:,:) :: mem => NULL()
    CHARACTER(LEN=*), PARAMETER::substr='scav_init_memory'
    CHARACTER(len=3)  :: strmod
    CHARACTER(LEN=2) :: nm
    CHARACTER(LEN=26), POINTER, DIMENSION(:) :: trname => NULL()
    INTEGER :: n2
    CHARACTER(LEN=2) :: modestr
    CHARACTER(LEN=STRLEN_MEDIUM) :: hlp

    REAL(dp) :: l_rtols, i_rtols

    INTRINSIC :: ASSOCIATED


    ! LOCAL    
    CALL start_message_bi(modstr,'MEMORY INITIALIZATION',substr)

    CALL initval_l
    CALL initval_i

    IF (.NOT.lscav) RETURN

    IF (lscav_gp .AND. ntrac_gp == 0) THEN
      lscav_gp = .FALSE.
      CALL info_bi( &
        'Gridpoint scavenging switched off - no tracers available', substr)
    ENDIF
    IF (lscav_lg .AND. ntrac_lg == 0) THEN
      lscav_lg = .FALSE.
!!#D attila +
#ifdef ECHAM5
      CALL info_bi( &
        'Lagrangian scavenging switched off - no tracers available', substr)
#endif
!!#D attila -
    ENDIF

    nsets = 0
    IF (lscav_gp) nsets = nsets + 1
    IF (lscav_lg) nsets = nsets + 1

    ! solver accuracy for liquid kpp
    l_rtols = 1.e-2_dp
    DO i=1,L_NVAR
       l_rtol(i) = l_rtols      ! must be defined for ros3
       l_atol(i) = 1.e1_dp      ! must be defined for ros3
    END DO
    ! solver accuracy for ice kpp
    i_rtols = 1.e-2_dp
    DO i=1,I_NVAR
       i_rtol(i) = i_rtols      ! must be defined for ros3
       i_atol(i) = 1.e1_dp      ! must be defined for ros3
    END DO
    
    CALL strcrack(out_string, ';', strname, n_out)
    ALLOCATE (kpp_out_nr(n_out))

    ALLOCATE(attr(nsets))

    ALLOCATE(set(nsets))

    DO js=1,nsets
      SELECT CASE (nsets)
      CASE(1)
        IF (lscav_gp) &
          ntrac = ntrac_gp
        IF (lscav_lg) &
          ntrac = ntrac_lg
      CASE(2)
        IF (js == 1)  ntrac = ntrac_gp
        IF (js == 2)  ntrac = ntrac_lg
      END SELECT

      ALLOCATE(attr(js)%log_att(ntrac,log_max))
      ALLOCATE(attr(js)%int_att(ntrac,int_max))
      ALLOCATE(attr(js)%evap_trac(n_resmod,ntrac))

      ALLOCATE(set(js)%para(2))
      DO i=1,2
        ALLOCATE (set(js)%para(i)%proc(4))
        DO j = 1, 4
          ALLOCATE &
            (set(js)%para(i)%proc(j)%evafrac(_RI_XYZN_(nproma,ngpblks,nlev,ntrac),n_resmod))
          set(js)%para(i)%proc(j)%evafrac(:,:,:,:,:) = 0.0_dp ! op_pj_20170112

        END DO
      END DO
    ENDDO
    ALLOCATE(process(4))
    
    DO j  = 1, 4
        ALLOCATE(process(j)%frac_eva(ntrac, n_resmod))
    END DO
    
    ALLOCATE(L_spec(0:LSPEC,nproma))
    ALLOCATE(I_spec(0:ISPEC,nproma))

    KSPEC = MAX(LSPEC,ISPEC)

    IF (LSCAV_LS) THEN
      ALLOCATE(ls_field(nsets,2))
      DO js=1,nsets
        DO i=1,2
          ALLOCATE(ls_field(js,i)%cloud_field(nproma, KSPEC, NLEV, NGPBLKS))
          ls_field(js,i)%cloud_field  = 0._dp
        END DO
      END DO
    ENDIF
    IF (LSCAV_CV) THEN
      ALLOCATE(cv_field(nsets,2))
      DO js=1,nsets
        DO i=1,2
          ALLOCATE(cv_field(js,i)%cloud_field(nproma, KSPEC, NLEV, NGPBLKS))
          cv_field(js,i)%cloud_field  = 0._dp
        END DO
      END DO
    ENDIF

    ALLOCATE(spec_field_ls(_RI_XYZN_(nproma, ngpblks, KSPEC, nsets), 2))
    ALLOCATE(spec_field_cv(_RI_XYZN_(nproma, ngpblks, KSPEC, nsets), 2))

    spec_field_ls = 0._dp
    spec_field_cv = 0._dp   

    ALLOCATE(special_flux(nsets))
    ALLOCATE(ph(nsets))
    ALLOCATE(KPP_INFO(nsets))
    ALLOCATE(wet_flx(n_out,nsets))

!----------------------------------------

    ! DEFAULT: USE Henry from CHEMPROP
    IF (loverwrite_henry) CALL  ALLOC_SCAV_VALUES_L
    ! call this simplifying subroutine, as values over ice need a
    ! special treatment
    CALL  ALLOC_SCAV_VALUES_I

    IF (ISCAV_EASY) CALL SET_EASY_SCAV_COEFFICIENTS_I
    IF (lscav_easy) THEN
      CALL SET_EASY_SCAV_COEFFICIENTS
      CALL DEFINE_Ks_for_effective_henry
    ENDIF

! In this section it is tested if the used gasphase species of the 
! scavenging code already exist in the model and an IDT should 
! be given to those. The species of the scavenging code that do not exist 
! MECCA were given a dummy-idt to prevent a crashing or memory leak or 
! misplaced memory during compiling and runtime
!
! If a tracer has the flag not to be scavenged (nscav=0) it won't be given an 
! IDT_SCAV and will therefore not be treated in the scavenging KPP mechanism
    CALL info_bi('... setting up scavenging behaviour of tracers !', substr)

!    for aerosols the important parameter is also the nscav switch which is locally stored
!    in the lwetdep array

    l_anynconvect  = .FALSE.
    ALLOCATE(aer_count(nsets)) 

    DO js=1,nsets
      attr(js)%log_att(:,lwetdep)       = .FALSE.
      attr(js)%log_att(:,laerosol)      = .FALSE.
      attr(js)%int_att(:,lquantity)     = 1
      attr(js)%int_att(:,tmode)         = 0
      attr(js)%int_att(:,aerosol_index) = 0 

      ALLOCATE(aer_count(js)%numbers(aercount_max))

      aer_count(js)%numbers(:) = 0
    ENDDO

    IF (lscav_gp) THEN
      DO jt = 1, ntrac_gp
        IF (ti_gp(jt)%tp%meta%cask_i(I_convect)==ON) l_anynconvect   = .TRUE.
      ENDDO
    ENDIF

    DO js=1,nsets
      SELECT CASE (nsets)
      CASE(1)
        IF (lscav_gp) THEN
          ntrac = ntrac_gp
          ti    => ti_gp
        ENDIF
        IF (lscav_lg) THEN
          ntrac = ntrac_lg
          ti    => ti_lg
        ENDIF
      CASE(2)
        IF (js == 1)  THEN
          ntrac = ntrac_gp
          ti    => ti_gp
        ENDIF
        IF (js == 2)  THEN
          ntrac = ntrac_lg
          ti    => ti_lg
        ENDIF
      END SELECT

      DO jt=1,ntrac
        IF ((ti(jt)%tp%ident%medium==AEROSOL)) &
          attr(js)%log_att(jt,laerosol) = .TRUE.
        IF ((ti(jt)%tp%meta%cask_i(I_scav)==ON))        &
          attr(js)%log_att(jt,lwetdep)  = .TRUE. 
        IF ((ti(jt)%tp%ident%type==FAMILY))    &
          attr(js)%log_att(jt,lwetdep)  = .FALSE.
        IF (attr(js)%log_att(jt,lwetdep) .AND. &
          attr(js)%log_att(jt,laerosol) )   THEN
          aer_count(js)%numbers(c_all) = aer_count(js)%numbers(c_all) + 1
        ENDIF
        attr(js)%int_att(jt,lquantity) = ti(jt)%tp%ident%quantity
        attr(js)%int_att(jt,tmode)     = ti(jt)%tp%meta%cask_i(I_aerosol_mode)
      END DO
    ENDDO

#ifdef MESSYTENDENCY
    IF (lscav_gp) THEN
      CALL mtend_register(tend_cv, mtend_id_tracer)
      CALL mtend_register(tend_ls, mtend_id_tracer)
      CALL mtend_register(tend_ev, mtend_id_tracer)
    ENDIF
#endif
    
    ALLOCATE(aer_attr(nsets))
    ALLOCATE(kpp_l_idx(nsets))
    ALLOCATE(kpp_i_idx(nsets))
    ALLOCATE(idx_evap_num(nsets))
    ALLOCATE(aermodel(nsets))   
    
    DO js=1,nsets
      SELECT CASE (nsets)
      CASE(1)
        IF (lscav_gp) THEN
          ntrac = ntrac_gp
          ti    => ti_gp
          strmod = '_gp'
        ENDIF
        IF (lscav_lg) THEN
          ntrac = ntrac_lg
          ti    => ti_lg
          strmod = '_lg'
        ENDIF
      CASE(2)
        IF (js == 1)  THEN
          ntrac = ntrac_gp
          ti    => ti_gp
          strmod = '_gp'
        ENDIF
        IF (js == 2)  THEN
          ntrac = ntrac_lg
          ti    => ti_lg
          strmod = '_lg'
        ENDIF
      END SELECT
      
      ALLOCATE(aer_attr(js)%aer_spec(aer_count(js)%numbers(c_all),aer_max))
      ALLOCATE(aer_attr(js)%name(aer_count(js)%numbers(c_all)))
      ALLOCATE(aer_attr(js)%mw(aer_count(js)%numbers(c_all)))

      aer_attr(js)%name(:) = ''
      aer_attr(js)%mw(:)   = 0._dp
      i = 1
      DO jt=1,ntrac
        IF (attr(js)%log_att(jt,lwetdep) .AND. &
          attr(js)%log_att(jt,laerosol) ) THEN 
          aer_attr(js)%name(i) = TRIM(ti(jt)%tp%ident%fullname) 
          attr(js)%int_att(jt,aerosol_index) = i
          aer_attr(js)%aer_spec(i,aer_idx)   = jt
          IF (attr(js)%int_att(jt, lquantity) == 1) THEN
            aer_count(js)%numbers(c_mass) = aer_count(js)%numbers(c_mass) + 1
            IF ( .NOT. (MATCH_WILD( 'PASSAER*', aer_attr(js)%name(i)))) THEN
               aer_attr(js)%mw(i) = ti(jt)%tp%meta%cask_r(R_molarmass)
            ENDIF
          ELSE
            aer_count(js)%numbers(c_num)  = aer_count(js)%numbers(c_num)  + 1
          ENDIF
          i = i +1 
        ENDIF
      ENDDO
!---------------------------
! building output fields

      CALL new_channel(status, modstr//strmod, reprid=GP_2D_HORIZONTAL)
      CALL channel_halt(substr, status)
      
      CALL new_channel_object(status, modstr//strmod, 'wetflx_nitrate', &
        p2=special_flux(js)%nitrate_flx)
      CALL channel_halt(substr, status)
      CALL new_attribute(status, modstr//strmod, 'wetflx_nitrate', &
      'long_name', c='wet dep. flux of nitrate')
      CALL channel_halt(substr, status)
      CALL new_attribute(status, modstr//strmod, 'wetflx_nitrate', &
        'units', c='molecules / (m^2 * s)')
      CALL channel_halt(substr, status)

      CALL new_channel_object(status, modstr//strmod, 'wetflx_sum_nitrate', &
        p2=special_flux(js)%nitrate_flx_sum, lrestreq=.TRUE.)
      CALL channel_halt(substr, status)
      CALL new_attribute(status, modstr//strmod, 'wetflx_sum_nitrate', &
        'long_name', c='time integral of wet dep. flux of nitrate')
      CALL channel_halt(substr, status)
      CALL new_attribute(status, modstr//strmod, 'wetflx_sum_nitrate', &
        'units', c='molecules / m^2')
      CALL channel_halt(substr, status)

      CALL new_channel_object(status, modstr//strmod, 'wetflx_sulfate', &
        p2=special_flux(js)%sulfate_flx)
      CALL channel_halt(substr, status)
      CALL new_attribute(status, modstr//strmod, 'wetflx_sulfate', &
        'long_name', c='wet dep. flux of sulfate')
      CALL channel_halt(substr, status)
      CALL new_attribute(status, modstr//strmod, 'wetflx_sulfate', &
        'units', c='molecules / (m^2 * s)')
      CALL channel_halt(substr, status)
      
      CALL new_channel_object(status, modstr//strmod, 'wetflx_sum_sulfate', &
        p2=special_flux(js)%sulfate_flx_sum, lrestreq=.TRUE.)
      CALL channel_halt(substr, status)
      CALL new_attribute(status, modstr//strmod, 'wetflx_sum_sulfate', &
        'long_name', c='time integral of wet dep. flux of sulfate')
      CALL channel_halt(substr, status)
      CALL new_attribute(status, modstr//strmod, 'wetflx_sum_sulfate', &
        'units', c='molecules / m^2')
      CALL channel_halt(substr, status)

      CALL new_channel_object(status, modstr//strmod, 'wetflx_ammoni', &
        p2=special_flux(js)%ammoni_flx)
      CALL channel_halt(substr, status)
      CALL new_attribute(status, modstr//strmod, 'wetflx_ammoni', &
        'long_name', c='wet dep. flux of ammonia/ammonium')
      CALL channel_halt(substr, status)
      CALL new_attribute(status, modstr//strmod, 'wetflx_ammoni', &
        'units', c='molecules / (m^2 * s)')
      CALL channel_halt(substr, status)
      
      CALL new_channel_object(status, modstr//strmod, 'wetflx_sum_ammoni', &
        p2=special_flux(js)%ammoni_flx_sum, lrestreq=.TRUE.)
      CALL channel_halt(substr, status)
      CALL new_attribute(status, modstr//strmod, 'wetflx_sum_ammoni', &
        'long_name', c='time integral of wet dep. flux of ammonia/ammonium')
      CALL channel_halt(substr, status)
      CALL new_attribute(status, modstr//strmod, 'wetflx_sum_ammoni', &
        'units', c='molecules / m^2')
      CALL channel_halt(substr, status)
      

      IF (lscav_ls) THEN
        kpp_rain_ls => spec_field_ls(:,:,:,:,1)
        kpp_snow_ls => spec_field_ls(:,:,:,:,2)
        DO i=1, n_out
          string = TRIM(strname(i))
          DO j=1, LSPEC
            IF ( TRIM(string) == TRIM(str_field_kpp_l(j))) THEN
              kpp_out_nr(i) = j
              mem => spec_field_ls(_RI_XYZN_(:,:,kpp_out_nr(i),js:js),1:1)
              CALL new_channel_object(status, modstr//strmod, &
                'wetflx_ls_'//TRIM(strname(i)),       &
                reprid=GP_3D_1LEV, mem=mem)
              CALL channel_halt(substr, status)
              CALL new_attribute(status, modstr//strmod,      &
                'wetflx_ls_'//TRIM(strname(i)),       &
                'long_name', c='ls wet dep. flux '//TRIM(strname(i)) )
              CALL channel_halt(substr, status)
              CALL new_attribute(status, modstr//strmod,      &
                'wetflx_ls_'//TRIM(strname(i)),       &
                'units', c='molecules /(m^2 * s)')
              CALL channel_halt(substr, status)
              
              CALL new_channel_object(status, modstr//strmod, &
                'wetflx_ls_sum_'//TRIM(strname(i)),   &
                p2=wet_flx(i,js)%flux_ls_sum, lrestreq=.TRUE.)
              CALL channel_halt(substr, status)
              CALL new_attribute(status, modstr//strmod,      &
                'wetflx_ls_sum_'//TRIM(strname(i)),   &
                'long_name', c='time integral of ls wet dep. flux'//&
                &TRIM(strname(i)) )
              CALL channel_halt(substr, status)
              CALL new_attribute(status, modstr//strmod,      &
                'wetflx_ls_sum_'//TRIM(strname(i)),   &
                'units', c='molecules /(m^2)')
              CALL channel_halt(substr, status)
            ENDIF
          ENDDO
        ENDDO

        loop_modes: DO j = 1, n_resmod

          jt = 0
          loop_tracers: DO

            CALL int2str(modestr, j, '0', 'X')
            IF (n_resmod > 2) THEN
              jt = jt + 1
              hlp = '_'//TRIM(aer_attr(js)%name(jt))//'_m'//modestr
            ELSE
              jt = aer_count(js)%numbers(c_all)
              hlp = '_m'//modestr
            END IF

            i = aer_attr(js)%aer_spec(jt,aer_idx)

            IF (lscav_imp) THEN
              IF (LSCAV_I) THEN
                mem => set(js)%para(1)%proc(2)%evafrac(_RI_XYZN_(:,:,:,i),j:j)
                CALL new_channel_object(status, modstr//strmod, &
                     'frac_evap_snow_ls'//hlp, mem=mem, reprid=GP_3D_MID)
                CALL channel_halt(substr, status)
                CALL new_attribute(status, modstr//strmod,      &
                     'frac_evap_snow_ls'//hlp, 'long_name', &
                     c='fraction of evaporation into mode '//modestr//' (large scale snow)' )
                CALL channel_halt(substr, status)
                CALL new_attribute(status, modstr//strmod,      &
                     'frac_evap_snow_ls'//hlp, 'units', c='-')
                CALL channel_halt(substr, status)
              ENDIF
              IF (LSCAV_L) THEN
                mem => set(js)%para(1)%proc(1)%evafrac(_RI_XYZN_(:,:,:,i),j:j)
                CALL new_channel_object(status, modstr//strmod, &
                     'frac_evap_rain_ls'//hlp, mem=mem, reprid=GP_3D_MID)
                CALL channel_halt(substr, status)
                CALL new_attribute(status, modstr//strmod,      &
                     'frac_evap_rain_ls'//hlp, 'long_name', &
                     c='fraction of evaporation into mode '//modestr//' (large scale rain)' )
                CALL channel_halt(substr, status)
                CALL new_attribute(status, modstr//strmod,      &
                     'frac_evap_rain_ls'//hlp, 'units', c='-')
                CALL channel_halt(substr, status)
              ENDIF
            END IF
            IF (lscav_nuc) THEN
              IF (LSCAV_I) THEN
                mem => set(js)%para(1)%proc(4)%evafrac(_RI_XYZN_(:,:,:,i),j:j)
                CALL new_channel_object(status, modstr//strmod, &
                     'frac_evap_iwc_ls'//hlp, mem=mem, reprid=GP_3D_MID)
                CALL channel_halt(substr, status)
                CALL new_attribute(status, modstr//strmod,      &
                     'frac_evap_iwc_ls'//hlp, 'long_name', &
                     c='fraction of evaporation into mode '//modestr//' (large scale cloud ice)' )
                CALL channel_halt(substr, status)
                CALL new_attribute(status, modstr//strmod,      &
                     'frac_evap_iwc_ls'//hlp, 'units', c='-')
                CALL channel_halt(substr, status)
              ENDIF
              IF (LSCAV_L) THEN
                mem => set(js)%para(1)%proc(3)%evafrac(_RI_XYZN_(:,:,:,i),j:j)
                CALL new_channel_object(status, modstr//strmod, &
                     'frac_evap_lwc_ls'//hlp, mem=mem, reprid=GP_3D_MID)
                CALL channel_halt(substr, status)
                CALL new_attribute(status, modstr//strmod,      &
                     'frac_evap_lwc_ls'//hlp, 'long_name', &
                     c='fraction of evaporation into mode '//modestr//' (large scale cloud water)' )
                CALL channel_halt(substr, status)
                CALL new_attribute(status, modstr//strmod,      &
                     'frac_evap_lwc_ls'//hlp, 'units', c='-')
                CALL channel_halt(substr, status)
              ENDIF
            END IF

            IF (jt .EQ. aer_count(js)%numbers(c_all)) EXIT loop_tracers

          END DO loop_tracers

        END DO loop_modes

      END IF

      IF (lscav_cv) THEN
        kpp_rain_cv => spec_field_cv(:,:,:,:,1)
        kpp_snow_cv => spec_field_cv(:,:,:,:,2)
        DO i=1, n_out
          string = TRIM(strname(i))
          DO j=1, LSPEC
            IF ( TRIM(string) == TRIM(str_field_kpp_l(j))) THEN
              kpp_out_nr(i) = j
              mem => spec_field_cv(_RI_XYZN_(:,:,kpp_out_nr(i),js:js),1:1)
              CALL new_channel_object(status, modstr//strmod,   &
                'wetflx_cv_'//TRIM(strname(i)),         &
                reprid=GP_3D_1LEV, mem=mem)
              CALL channel_halt(substr, status)
              CALL new_attribute(status, modstr//strmod,        &
                'wetflx_cv_'//TRIM(strname(i)),         &
                'long_name', c='cv wet dep. flux '//TRIM(strname(i)) )
              CALL channel_halt(substr, status)
              CALL new_attribute(status, modstr//strmod,        &
                'wetflx_cv_'//TRIM(strname(i)),         &
                'units', c='molecules /(m^2 * s)')
              CALL channel_halt(substr, status)
              
              CALL new_channel_object(status, modstr//strmod,   &
                'wetflx_cv_sum_'//TRIM(strname(i)),     &
                p2=wet_flx(i,js)%flux_cv_sum, lrestreq=.TRUE.)
              CALL channel_halt(substr, status)
              CALL new_attribute(status, modstr//strmod,        &
                'wetflx_cv_sum_'//TRIM(strname(i)),     &
                'long_name',                            &
                c='time integral of cv wet dep. flux '//TRIM(strname(i)))
              CALL channel_halt(substr, status)
              CALL new_attribute(status, modstr//strmod,        &
                'wetflx_cv_sum_'//TRIM(strname(i)),     &
                'units', c='molecules / m^2')
              CALL channel_halt(substr, status)
            ENDIF
          ENDDO
        ENDDO

        loop2_modes: DO j = 1, n_resmod

          jt = 0
          loop2_tracers: DO

            CALL int2str(modestr, j, '0', 'X')
            IF (n_resmod > 2) THEN
              jt = jt + 1
              hlp = '_'//TRIM(aer_attr(js)%name(jt))//'_m'//modestr
            ELSE
              jt = aer_count(js)%numbers(c_all)
              hlp = '_m'//modestr
            END IF

            i = aer_attr(js)%aer_spec(jt,aer_idx)

            IF (lscav_imp) THEN
              IF (LSCAV_I) THEN
                mem => set(js)%para(2)%proc(2)%evafrac(_RI_XYZN_(:,:,:,i),j:j)
                CALL new_channel_object(status, modstr//strmod, &
                     'frac_evap_snow_cv'//hlp, mem=mem, reprid=GP_3D_MID)
                CALL channel_halt(substr, status)
                CALL new_attribute(status, modstr//strmod,      &
                     'frac_evap_snow_cv'//hlp, 'long_name', &
                     c='fraction of evaporation into mode '//modestr//' (convective snow)' )
                CALL channel_halt(substr, status)
                CALL new_attribute(status, modstr//strmod,      &
                     'frac_evap_snow_cv'//hlp, 'units', c='-')
                CALL channel_halt(substr, status)
              ENDIF
              IF (LSCAV_L) THEN
                mem => set(js)%para(2)%proc(1)%evafrac(_RI_XYZN_(:,:,:,i),j:j)
                CALL new_channel_object(status, modstr//strmod, &
                     'frac_evap_rain_cv'//hlp, mem=mem, reprid=GP_3D_MID)
                CALL channel_halt(substr, status)
                CALL new_attribute(status, modstr//strmod,      &
                     'frac_evap_rain_cv'//hlp, 'long_name', &
                     c='fraction of evaporation into mode '//modestr//' (convective rain)' )
                CALL channel_halt(substr, status)
                CALL new_attribute(status, modstr//strmod,      &
                     'frac_evap_rain_cv'//hlp, 'units', c='-')
                CALL channel_halt(substr, status)
              ENDIF
            END IF
            IF (lscav_nuc) THEN
              IF (LSCAV_I) THEN
                mem => set(js)%para(2)%proc(4)%evafrac(_RI_XYZN_(:,:,:,i),j:j)
                CALL new_channel_object(status, modstr//strmod, &
                     'frac_evap_iwc_cv'//hlp, mem=mem, reprid=GP_3D_MID)
                CALL channel_halt(substr, status)
                CALL new_attribute(status, modstr//strmod,      &
                     'frac_evap_iwc_cv'//hlp, 'long_name', &
                     c='fraction of evaporation into mode '//modestr//' (convective cloud ice)' )
                CALL channel_halt(substr, status)
                CALL new_attribute(status, modstr//strmod,      &
                     'frac_evap_iwc_cv'//hlp, 'units', c='-')
                CALL channel_halt(substr, status)
              ENDIF
              IF (LSCAV_L) THEN
                mem => set(js)%para(2)%proc(3)%evafrac(_RI_XYZN_(:,:,:,i),j:j)
                CALL new_channel_object(status, modstr//strmod, &
                     'frac_evap_lwc_cv'//hlp, mem=mem, reprid=GP_3D_MID)
                CALL channel_halt(substr, status)
                CALL new_attribute(status, modstr//strmod,      &
                     'frac_evap_lwc_cv'//hlp, 'long_name', &
                     c='fraction of evaporation into mode '//modestr//' (convective cloud water)' )
                CALL channel_halt(substr, status)
                CALL new_attribute(status, modstr//strmod,      &
                     'frac_evap_lwc_cv'//hlp, 'units', c='-')
                CALL channel_halt(substr, status)
              ENDIF
            END IF

            IF (jt .EQ. aer_count(js)%numbers(c_all)) EXIT loop2_tracers

          END DO loop2_tracers

        END DO loop2_modes
! op_ck_20140522-

      ENDIF

      IF (lscav_ls) THEN
        CALL new_channel_object(status, modstr//strmod, 'cloudpH_ls',  &
          p3=ph(js)%cloudpH_ls, reprid=GP_3D_MID)
        CALL channel_halt(substr, status)
        CALL new_attribute(status, modstr//strmod, 'cloudpH_ls',       &
          'long_name', c='large scale cloud pH')
        CALL channel_halt(substr, status)
        CALL new_attribute(status, modstr//strmod, 'cloudpH_ls',       &
          'units', c='0-14, NA=19.99')
        CALL channel_halt(substr, status)
        
        CALL new_channel_object(status, modstr//strmod, 'Hp_cloud_ls', &
          p3=ph(js)%Hp_cloud_ls, reprid=GP_3D_MID)
        CALL channel_halt(substr, status)
        CALL new_attribute(status, modstr//strmod, 'Hp_cloud_ls',      &
          'long_name', c='large scale cloud H+ concentration')
        CALL channel_halt(substr, status)
        CALL new_attribute(status, modstr//strmod, 'Hp_cloud_ls',      &
          'units', c='molecules/kg(air)')
        CALL channel_halt(substr, status)

        CALL new_channel_object(status, modstr//strmod, &
          'KPP_steps_cloud_ls', &
          p3=kpp_info(js)%ls_cloud_steps, reprid=GP_3D_MID)
        CALL channel_halt(substr, status)
        CALL new_attribute(status, modstr//strmod, 'KPP_steps_cloud_ls',&
          'long_name', c='# of KPP steps for ls cloud integration')
        CALL channel_halt(substr, status)
        CALL new_attribute(status, modstr//strmod, 'KPP_steps_cloud_ls', &
          'units', c=' - ')
        CALL channel_halt(substr, status)

        CALL new_channel_object(status, modstr//strmod, &
          'KPP_rsteps_cloud_ls', &
          p3=kpp_info(js)%ls_cloud_rsteps, reprid=GP_3D_MID)
        CALL channel_halt(substr, status)
        CALL new_attribute(status, modstr//strmod, 'KPP_rsteps_cloud_ls',&
          'long_name', &
          c='# of rejected KPP steps for ls cloud integration')
        CALL channel_halt(substr, status)
        CALL new_attribute(status, modstr//strmod, 'KPP_rsteps_cloud_ls',&
          'units', c=' - ')
        CALL channel_halt(substr, status)
        
        CALL new_channel_object(status, modstr//strmod, 'rainpH_ls',  &
          p3=ph(js)%rainpH_ls, reprid=GP_3D_MID)
        CALL channel_halt(substr, status)
        CALL new_attribute(status, modstr//strmod, 'rainpH_ls',       &
          'long_name', c='large scale rain pH')
        CALL channel_halt(substr, status)
        CALL new_attribute(status, modstr//strmod, 'rainpH_ls',       &
          'units', c='0-14, NA=19.99')
        CALL channel_halt(substr, status)

        CALL new_channel_object(status, modstr//strmod, &
          'KPP_steps_rain_ls', &
          p3=kpp_info(js)%ls_rain_steps, reprid=GP_3D_MID)
        CALL channel_halt(substr, status)
        CALL new_attribute(status, modstr//strmod, 'KPP_steps_rain_ls',&
          'long_name', c='# of KPP steps for ls rain integration')
        CALL channel_halt(substr, status)
        CALL new_attribute(status, modstr//strmod, 'KPP_steps_rain_ls', &
          'units', c=' - ')
        CALL channel_halt(substr, status)

        CALL new_channel_object(status, modstr//strmod, &
          'KPP_rsteps_rain_ls', &
          p3=kpp_info(js)%ls_rain_rsteps, reprid=GP_3D_MID)
        CALL channel_halt(substr, status)
        CALL new_attribute(status, modstr//strmod, 'KPP_rsteps_rain_ls',&
          'long_name', &
          c='# of rejected KPP steps for ls rain integration')
        CALL channel_halt(substr, status)
        CALL new_attribute(status, modstr//strmod, 'KPP_rsteps_rain_ls',&
          'units', c=' - ')
        CALL channel_halt(substr, status)

        ! for output of partitioning parameters
        CALL new_channel_object(status, modstr//strmod, &
          'phi_i_ls', &
          p3=phi_i_ls, reprid=GP_3D_MID)
        CALL channel_halt(substr, status)
        CALL new_attribute(status, modstr//strmod, 'phi_i_ls',&
          'long_name', &
          c='large scale [HNO3_ice]/[HNO3_tot]')
        CALL channel_halt(substr, status)
        CALL new_attribute(status, modstr//strmod, 'phi_i_ls',&
          'units', c=' - ')
        CALL channel_halt(substr, status)
        
        CALL new_channel_object(status, modstr//strmod, &
          'mju_i_ls', &
          p3=mju_i_ls, reprid=GP_3D_MID)
        CALL channel_halt(substr, status)
        CALL new_attribute(status, modstr//strmod, 'mju_i_ls',&
          'long_name', &
          c='large scale [HNO3_ice]/[H2O_ice]')
        CALL channel_halt(substr, status)
        CALL new_attribute(status, modstr//strmod, 'mju_i_ls',&
          'units', c=' - ')
        CALL channel_halt(substr, status)

        CALL new_channel_object(status, modstr//strmod, &
          'iwc_ls', &
          p3=iwc_ls, reprid=GP_3D_MID)
        CALL channel_halt(substr, status)
        CALL new_attribute(status, modstr//strmod, 'iwc_ls',&
          'long_name', &
          c='large scale IWC')
        CALL channel_halt(substr, status)
        CALL new_attribute(status, modstr//strmod, 'iwc_ls',&
          'units', c='kg/m^3')
        CALL channel_halt(substr, status)

        CALL new_channel_object(status, modstr//strmod, &
          'iwc_T_ls', &
          p3=iwc_T_ls, reprid=GP_3D_MID)
        CALL channel_halt(substr, status)
        CALL new_attribute(status, modstr//strmod, 'iwc_T_ls',&
          'long_name', &
          c='T based large scale IWC')
        CALL channel_halt(substr, status)
        CALL new_attribute(status, modstr//strmod, 'iwc_T_ls',&
          'units', c='kg/m^3')
        CALL channel_halt(substr, status)
      ENDIF
      
      IF (lscav_cv) THEN
        CALL new_channel_object(status, modstr//strmod, 'cloudpH_cv', &
          p3=ph(js)%cloudpH_cv, reprid=GP_3D_MID)
        CALL channel_halt(substr, status)
        CALL new_attribute(status, modstr//strmod, 'cloudpH_cv',      &
          'long_name', c='convective cloud pH')
        CALL channel_halt(substr, status)
        CALL new_attribute(status, modstr//strmod, 'cloudpH_cv',      &
          'units', c='0-14, NA=19.99')
        CALL channel_halt(substr, status)
        
        CALL new_channel_object(status, modstr//strmod, 'Hp_cloud_cv',&
          p3=ph(js)%Hp_cloud_cv, reprid=GP_3D_MID)
        CALL channel_halt(substr, status)
        CALL new_attribute(status, modstr//strmod, 'Hp_cloud_cv',     &
          'long_name', c='convective cloud H+ concentration')
        CALL channel_halt(substr, status)
        CALL new_attribute(status, modstr//strmod, 'Hp_cloud_cv',     &
          'units', c='molecules/kg(air)')
        CALL channel_halt(substr, status)
        
        CALL new_channel_object(status, modstr//strmod, 'rainpH_cv',  &
          p3=ph(js)%rainpH_cv, reprid=GP_3D_MID)
        CALL channel_halt(substr, status)
        CALL new_attribute(status, modstr//strmod, 'rainpH_cv',       &
          'long_name', c='convective rain pH')
        CALL channel_halt(substr, status)
        CALL new_attribute(status, modstr//strmod, 'rainpH_cv',       &
          'units', c='0-14, NA=19.99')
        CALL channel_halt(substr, status)

        CALL new_channel_object(status, modstr//strmod, &
          'KPP_steps_cloud_cv', &
          p3=kpp_info(js)%cv_cloud_steps, reprid=GP_3D_MID)
        CALL channel_halt(substr, status)
        CALL new_attribute(status, modstr//strmod, 'KPP_steps_cloud_cv',&
          'long_name', c='# of KPP steps for cv cloud integration')
        CALL channel_halt(substr, status)
        CALL new_attribute(status, modstr//strmod, 'KPP_steps_cloud_cv', &
          'units', c=' - ')
        CALL channel_halt(substr, status)

        CALL new_channel_object(status, modstr//strmod, &
          'KPP_rsteps_cloud_cv', &
          p3=kpp_info(js)%cv_cloud_rsteps, reprid=GP_3D_MID)
        CALL channel_halt(substr, status)
        CALL new_attribute(status, modstr//strmod, 'KPP_rsteps_cloud_cv',&
          'long_name', &
          c='# of rejected KPP steps for cv cloud integration')
        CALL channel_halt(substr, status)
        CALL new_attribute(status, modstr//strmod, 'KPP_rsteps_cloud_cv',&
          'units', c=' - ')
        CALL channel_halt(substr, status)

        CALL new_channel_object(status, modstr//strmod, &
          'KPP_steps_rain_cv', &
          p3=kpp_info(js)%cv_rain_steps, reprid=GP_3D_MID)
        CALL channel_halt(substr, status)
        CALL new_attribute(status, modstr//strmod, 'KPP_steps_rain_cv',&
          'long_name', c='# of KPP steps for cv rain integration')
        CALL channel_halt(substr, status)
        CALL new_attribute(status, modstr//strmod, 'KPP_steps_rain_cv', &
          'units', c=' - ')
        CALL channel_halt(substr, status)

        CALL new_channel_object(status, modstr//strmod, &
          'KPP_rsteps_rain_cv', &
          p3=kpp_info(js)%cv_rain_rsteps, reprid=GP_3D_MID)
        CALL channel_halt(substr, status)
        CALL new_attribute(status, modstr//strmod, 'KPP_rsteps_rain_cv',&
          'long_name', &
          c='# of rejected KPP steps for cv rain integration')
        CALL channel_halt(substr, status)
        CALL new_attribute(status, modstr//strmod, 'KPP_rsteps_rain_cv',&
          'units', c=' - ')
        CALL channel_halt(substr, status)

        ! for output of partitioning parameters
        CALL new_channel_object(status, modstr//strmod, &
          'phi_i_cv', &
          p3=phi_i_cv, reprid=GP_3D_MID)
        CALL channel_halt(substr, status)
        CALL new_attribute(status, modstr//strmod, 'phi_i_cv',&
          'long_name', c='convective [HNO3_ice]/[HNO3_tot]')
        CALL channel_halt(substr, status)
        CALL new_attribute(status, modstr//strmod, 'phi_i_cv',&
          'units', c=' - ')
        CALL channel_halt(substr, status)
        
        CALL new_channel_object(status, modstr//strmod, &
          'mju_i_cv', &
          p3=mju_i_cv, reprid=GP_3D_MID)
        CALL channel_halt(substr, status)
        CALL new_attribute(status, modstr//strmod, 'mju_i_cv',&
          'long_name', c='convective [HNO3_ice]/[H2O_ice]')
        CALL channel_halt(substr, status)
        CALL new_attribute(status, modstr//strmod, 'mju_i_cv',&
          'units', c=' - ')
        CALL channel_halt(substr, status)

        CALL new_channel_object(status, modstr//strmod, &
          'iwc_cv', &
          p3=iwc_cv, reprid=GP_3D_MID)
        CALL channel_halt(substr, status)
        CALL new_attribute(status, modstr//strmod, 'iwc_cv',&
          'long_name', c='convective IWC')
        CALL channel_halt(substr, status)
        CALL new_attribute(status, modstr//strmod, 'iwc_cv',&
          'units', c='kg/m^3')
        CALL channel_halt(substr, status)

        CALL new_channel_object(status, modstr//strmod, &
          'iwc_T_cv', &
          p3=iwc_T_cv, reprid=GP_3D_MID)
        CALL channel_halt(substr, status)
        CALL new_attribute(status, modstr//strmod, 'iwc_T_cv',&
          'long_name', c='T based convective IWC')
        CALL channel_halt(substr, status)
        CALL new_attribute(status, modstr//strmod, 'iwc_T_cv',&
          'units', c='kg/m^3')
        CALL channel_halt(substr, status)
      ENDIF

     L_TEND = (TRIM(te_string) /= '')
      IF (L_TEND .AND. (js==1)) THEN
         ! ... recycle vriable here
         IF (ASSOCIATED(strname)) DEALLOCATE (strname)
         NULLIFY(strname)
         ! get number and name(s) of tracers
         CALL strcrack(te_string, ';', strname, n_te)
         ! allocate space for tracer IDs
         ALLOCATE(te_idt(n_te))
         te_idt(:) = 0
         ! allocate space for pointers to tendency copies
         ALLOCATE(tte_scav(n_te))

         DO jt=1, n_te
            NULLIFY(tte_scav(jt)%ptr)
            CALL strcrack(strname(jt),'_', trname, n2)
            IF (n2==1) THEN 
               CALL get_tracer(status, GPTRSTR, TRIM(trname(1)),    &
                    idx=te_idt(jt))
             ELSE
               CALL get_tracer(status, GPTRSTR, TRIM(trname(1)),    &
                    subname=TRIM(trname(2)),idx=te_idt(jt))
            END IF
            IF (te_idt(jt) /= 0) THEN
               
               ! new channel object for preserving tendency from scav
               CALL new_channel_object(status, modstr//strmod &
                    , TRIM(strname(jt))//'_scte' &
                    , p3=tte_scav(jt)%ptr, reprid=GP_3D_MID)
               CALL channel_halt(substr, status)
               CALL new_attribute(status, modstr//strmod &
                    , TRIM(strname(jt))//'_scte' &
                    , 'long_name' &
                    , c=TRIM(strname(jt))//' tendency from scavenging')
               CALL channel_halt(substr, status)
               CALL new_attribute(status, modstr//strmod &
                    , TRIM(strname(jt))//'_scte' &
                    , 'units', c='mol/mol/s')
               CALL channel_halt(substr, status)

            ELSE
               
               CALL warning_bi('tracer '//TRIM(strname(jt))//' not present' &
                    ,substr)

            END IF
         END DO
      ENDIF
      
!-------------------------------------------------------------------------------
    ENDDO  ! nsets loop for kpp_species

    IF (ASSOCIATED(strname)) DEALLOCATE (strname)
    NULLIFY(strname)

!  building aerosol output fields
       
    CALL strcrack(out_string_aer, ';', strname, n_out_aer)

    ALLOCATE(wet_flx_aer(n_out_aer,nsets))
    ALLOCATE(aero_flx(nsets))
    !    allocate (aer_out_nr(n_out_aer))


    DO js = 1,nsets
      SELECT CASE (nsets)
      CASE(1)
        IF (lscav_gp) THEN
          ntrac  =  ntrac_gp
          ti     => ti_gp
          strmod = '_gp'
        ENDIF
        IF (lscav_lg) THEN
          ntrac  =  ntrac_lg
          ti     => ti_lg
          strmod = '_lg'
        ENDIF
      CASE(2)
        IF (js == 1)  THEN
          ntrac  =  ntrac_gp
          ti     => ti_gp
          strmod = '_gp'
        ENDIF
        IF (js == 2)  THEN
          ntrac  =  ntrac_lg
          ti     => ti_lg
          strmod = '_lg'
        ENDIF
      END SELECT

      IF (lscav_ls) THEN
        ALLOCATE(aero_flx(js)%aero_field_ls( &
             _RI_XYZN_( nproma,ngpblks, aer_count(js)%numbers(c_all), 2), 1))
        ALLOCATE(aero_flx(js)%cloud_field_ls_aer( &
             _RI_XYZN_(nproma, ngpblks,aer_count(js)%numbers(c_all), 2),nlev))
        aero_flx(js)%aero_field_ls      = 0._dp
        aero_flx(js)%cloud_field_ls_aer = 0._dp
        aero_flx(js)%aer_flx_ls => aero_flx(js)%aero_field_ls(:,:,:,:,1)

        DO i=1, n_out_aer
          DO j=1, aer_count(js)%numbers(c_all)
            IF ( TRIM(strname(i)) == TRIM(aer_attr(js)%name(j)) ) THEN
              wet_flx_aer(i,js)%aer_out_nr = j
              mem => &
              aero_flx(js)%aero_field_ls(&
              _RI_XYZN_(:,:,wet_flx_aer(i,js)%aer_out_nr,1:1),:)
              CALL new_channel_object(status, modstr//strmod,&
                'wetflx_aer_ls_'//TRIM(strname(i)), mem=mem, &
                reprid=GP_3D_1LEV)
              CALL channel_halt(substr, status)
              CALL new_attribute(status, modstr//strmod,     &
                'wetflx_aer_ls_'//TRIM(strname(i)),          &
                'long_name',                                 &
                c='ls aerosol wet dep. flux '//TRIM(strname(i)))
              CALL channel_halt(substr, status)
              CALL new_attribute(status, modstr//strmod,     &
                'wetflx_aer_ls_'//TRIM(strname(i)),          &
                'units', c='molecules / (m^2 * s)')
              CALL channel_halt(substr, status)
              
              CALL new_channel_object(status, modstr//strmod,&
                'wetflx_aer_ls_sum_'//TRIM(strname(i)),      &
                p2=wet_flx_aer(i,js)%flux_ls_aer_sum, lrestreq=.TRUE.)
              CALL channel_halt(substr, status)
              CALL new_attribute(status, modstr//strmod,     &
                'wetflx_aer_ls_sum_'//TRIM(strname(i)),      &
                'long_name',                                 &
                c='sum of ls aerosol wet dep. flux '//TRIM(strname(i)))
              CALL channel_halt(substr, status)
              CALL new_attribute(status, modstr//strmod,     &
                'wetflx_aer_ls_sum_'//TRIM(strname(i)),      &
                'units', c='molecules / m^2')
              CALL channel_halt(substr, status)
            ENDIF
          ENDDO
        ENDDO
        
      ENDIF
      
      IF (lscav_cv) THEN
        ALLOCATE(aero_flx(js)%aero_field_cv(&
             _RI_XYZN_(nproma, ngpblks, aer_count(js)%numbers(c_all), 2), 1))   
        ALLOCATE(aero_flx(js)%cloud_field_cv_aer(&
             _RI_XYZN_(nproma, ngpblks,aer_count(js)%numbers(c_all), 2), nlev))
        aero_flx(js)%aero_field_cv      = 0._dp
        aero_flx(js)%cloud_field_cv_aer = 0._dp
        aero_flx(js)%aer_flx_cv => aero_flx(js)%aero_field_cv(:,:,:,:,1)
        
        DO i=1, n_out_aer
          DO j=1, aer_count(js)%numbers(c_all)
            IF ( TRIM(strname(i)) ==  TRIM(aer_attr(js)%name(j)) ) THEN
              wet_flx_aer(i,js)%aer_out_nr = j
              mem => &
              aero_flx(js)%aero_field_cv(_RI_XYZN_(:,:,wet_flx_aer(i,js)%aer_out_nr,1:1),:)
              CALL new_channel_object(status, modstr//strmod, &
                'wetflx_aer_cv_'//TRIM(strname(i)), mem=mem,  &
                reprid=GP_3D_1LEV)
              CALL channel_halt(substr, status)
              CALL new_attribute(status, modstr//strmod,      &
                'wetflx_aer_cv_'//TRIM(strname(i)),           &
                'long_name',                                  &
                c='cv aerosol wet dep. flux '//TRIM(strname(i)))
              CALL channel_halt(substr, status)
              CALL new_attribute(status, modstr//strmod,      &
                'wetflx_aer_cv_'//TRIM(strname(i)),           &
                'units', c='molecules / (m^2 * s)')
              CALL channel_halt(substr, status)
              
              CALL new_channel_object(status, modstr//strmod, &
                'wetflx_aer_cv_sum_'//TRIM(strname(i)),       &
                p2=wet_flx_aer(i,js)%flux_cv_aer_sum, lrestreq=.TRUE.)
              CALL channel_halt(substr, status)
              CALL new_attribute(status, modstr//strmod,      &
                'wetflx_aer_cv_sum_'//TRIM(strname(i)),       &
                'long_name',                                  &
                c='sum of cv aerosol wet dep. flux '//TRIM(strname(i)))
              CALL channel_halt(substr, status)
              CALL new_attribute(status, modstr//strmod,      &
                'wetflx_aer_cv_sum_'//TRIM(strname(i)),       &
                'units', c='molecules / m^2')
              CALL channel_halt(substr, status)
            ENDIF
          ENDDO
        ENDDO
              
      ENDIF !lscav_cv
      
    ENDDO ! nsets
    
    IF (ASSOCIATED(strname))      DEALLOCATE(strname)

    !--------------- 
    ! scavenging parameters for all tracer sets
    !---------------
    CALL new_channel(status, modstr, reprid=GP_2D_HORIZONTAL)
    CALL channel_halt(substr, status)

    CALL new_channel_object(status, modstr, 'kconbot', p2=kconbot)
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'kconbot', &
      'long_name', c='index of layers of scav_convec')
    CALL channel_halt(substr, status)
    CALL new_attribute(status, modstr, 'kconbot', &
      'units', c='-')
    CALL channel_halt(substr, status)

    IF (lscav_ls) THEN
      CALL new_channel_object(status, modstr, 'bc_lwc_ls', &
        p3=bclwc_ls, reprid=GP_3D_MID)
      CALL channel_halt(substr, status)
      CALL new_attribute(status, modstr, 'bc_lwc_ls', &
        'long_name', c='liquid water content below cloud LS')
      CALL channel_halt(substr, status)
      CALL new_attribute(status, modstr, 'bc_lwc_ls',&
        'units', c='kg/m^3')
      CALL channel_halt(substr, status)
      
      CALL new_channel_object(status, modstr, 'oldcov', &
        p3=bc_cov, reprid=GP_3D_MID)
      CALL channel_halt(substr, status)
      CALL new_attribute(status, modstr, 'oldcov', &
        'long_name', c='precipitating column cloud cover')
      CALL channel_halt(substr, status)
      CALL new_attribute(status, modstr, 'oldcov', &
        'units', c='-')
      CALL channel_halt(substr, status)
    END IF

    IF (lscav_cv) THEN
      CALL new_channel_object(status, modstr, 'bc_lwc_cv',  &
        p3=bclwc_cv, reprid=GP_3D_MID)
      CALL channel_halt(substr, status)
      CALL new_attribute(status, modstr, 'bc_lwc_cv', &
        'long_name', c='liquid water content below cloud CV')
      CALL channel_halt(substr, status)
      CALL new_attribute(status, modstr, 'bc_lwc_cv', &
        'units', c='kg/m^3')
      CALL channel_halt(substr, status)

      CALL new_channel_object(status, modstr, 'lwc_cv', &
        p3=lwc_cv, reprid=GP_3D_MID)
      CALL channel_halt(substr, status)
      CALL new_attribute(status, modstr, 'lwc_cv', &
        'long_name', c='approx. convective cloud water content')
      CALL channel_halt(substr, status)
      CALL new_attribute(status, modstr, 'lwc_cv', &
        'units', c='kg(water)/kg(air)')
      CALL channel_halt(substr, status)
    END IF

    nprod = 0
    DO jt=1,lspec
      IF (MATCH_WILD( '*Prod*', str_field_kpp_l(jt) )) THEN
        nprod = nprod + 1
      END IF
    END DO


    ALLOCATE(PR(nsets))
    DO js = 1,nsets
      SELECT CASE (nsets)
      CASE(1)
        IF (lscav_gp) THEN
          ntrac  =  ntrac_gp
          ti     => ti_gp
          strmod = '_gp'
        ENDIF
      END SELECT

      ALLOCATE(PR(js)%PROD(nprod))
      ALLOCATE(PR(js)%PROD_IDX(nprod))
      ALLOCATE(PR(js)%PROD_CHAR(nprod))

      PR(js)%PROD_CHAR(:) = ""
      PR(js)%PROD_IDX(:)  = 0

      kprod = 0
      DO jt=1,lspec
        IF (ASSOCIATED(strname)) DEALLOCATE (strname)
        NULLIFY(strname)
        IF (MATCH_WILD( '*Prod*', str_field_kpp_l(jt) )) THEN
          kprod = kprod + 1
          CALL strcrack(str_field_kpp_l(jt), '_', strname, dummy)
          PR(js)%PROD_CHAR(kprod) = strname(2)
          PR(js)%PROD_IDX(kprod) = jt
        END IF
      END DO
      

      DO j = 1, nprod

        nm = PR(js)%PROD_CHAR(j)

        IF (lscav_ls) THEN
          CALL new_channel_object(status, modstr//strmod, &
            'CL_LS_Prod_'//TRIM(nm),       &
            p3 = PR(js)%PROD(j)%cl_ls, reprid=GP_3D_MID)
          CALL channel_halt(substr, status)
          CALL new_attribute(status, modstr//strmod, &
            'CL_LS_Prod_'//TRIM(nm), &
            'long_name', c=&
            'In-Cloud (LS) Production rate from Reaction '//TRIM(nm))
          CALL channel_halt(substr, status)
          CALL new_attribute(status, modstr//strmod, &
            'CL_LS_Prod_'//TRIM(nm), &
            'units', c='molecules/cm^3/s')
          CALL channel_halt(substr, status)

          CALL new_channel_object(status, modstr//strmod, &
            'RA_LS_Prod_'//TRIM(nm),       &
            p3 = PR(js)%PROD(j)%ra_ls, reprid=GP_3D_MID)
          CALL channel_halt(substr, status)
          CALL new_attribute(status, modstr//strmod, &
            'RA_LS_Prod_'//TRIM(nm), &
            'long_name', c=&
            'In-Rain (LS) Production rate from Reaction '//TRIM(nm))
          CALL channel_halt(substr, status)
          CALL new_attribute(status, modstr//strmod, &
            'RA_LS_Prod_'//TRIM(nm), &
            'units', c='molecules/cm^3/s')
          CALL channel_halt(substr, status)

        ENDIF
        IF (lscav_cv) THEN
          CALL new_channel_object(status, modstr//strmod, &
            'CL_CV_Prod_'//TRIM(nm),       &
            p3 = PR(js)%PROD(j)%cl_cv, reprid=GP_3D_MID)
          CALL channel_halt(substr, status)
          CALL new_attribute(status, modstr//strmod, &
            'CL_CV_Prod_'//TRIM(nm), &
            'long_name', c=&
            'In-Cloud (CV) Production rate from Reaction '//TRIM(nm))
          CALL channel_halt(substr, status)
          CALL new_attribute(status, modstr//strmod, &
            'CL_CV_Prod_'//TRIM(nm), &
            'units', c='molecules/cm^3/s')
          CALL channel_halt(substr, status)

          CALL new_channel_object(status, modstr//strmod, &
            'RA_CV_Prod_'//TRIM(nm),       &
            p3 = PR(js)%PROD(j)%ra_cv, reprid=GP_3D_MID)
          CALL channel_halt(substr, status)
          CALL new_attribute(status, modstr//strmod, &
            'RA_CV_Prod_'//TRIM(nm), &
            'long_name', c=&
            'In-Rain (CV) Production rate from Reaction '//TRIM(nm))
          CALL channel_halt(substr, status)
          CALL new_attribute(status, modstr//strmod, &
            'RA_CV_Prod_'//TRIM(nm), &
            'units', c='molecules/cm^3/s')
          CALL channel_halt(substr, status)

        ENDIF
      ENDDO
    ENDDO

    
    CALL end_message_bi(modstr,'MEMORY INITIALIZATION',substr)

  END SUBROUTINE scav_init_memory

!==============================================================================

  SUBROUTINE scav_free_memory

    USE messy_scav_inter,        ONLY: idx_evap_num, made3 &
                                     , process

    IMPLICIT NONE

    INTEGER   :: js, i, j

    INTRINSIC :: ALLOCATED, ASSOCIATED

    IF (.NOT.lscav) RETURN

    IF (ASSOCIATED(wet_flx))       DEALLOCATE (wet_flx)
    NULLIFY(wet_flx)
    IF (ASSOCIATED(wet_flx_aer))   DEALLOCATE (wet_flx_aer)
    NULLIFY(wet_flx_aer)
    IF (ALLOCATED(kpp_out_nr))     DEALLOCATE (kpp_out_nr)
    IF (ALLOCATED(idx_evap_num))   DEALLOCATE (idx_evap_num)
    IF (ALLOCATED(made3))          DEALLOCATE (made3)
    
    IF (ASSOCIATED(spec_field_ls)) DEALLOCATE(spec_field_ls)
    NULLIFY(spec_field_ls)
    IF (ASSOCIATED(spec_field_cv)) DEALLOCATE(spec_field_cv)
    NULLIFY(spec_field_cv)
    DO js=1,nsets

      DO i = 1, 2
        DO j = 1, 4
          IF (ASSOCIATED(set(js)%para(i)%proc(j)%evafrac)) &
            DEALLOCATE(set(js)%para(i)%proc(j)%evafrac)
          NULLIFY(set(js)%para(i)%proc(j)%evafrac)
        END DO
      END DO

      IF (ASSOCIATED(aero_flx(js)%aero_field_ls)) &
        DEALLOCATE(aero_flx(js)%aero_field_ls)
        NULLIFY(aero_flx(js)%aero_field_ls)
      IF (ASSOCIATED(aero_flx(js)%aero_field_cv)) &
        DEALLOCATE(aero_flx(js)%aero_field_cv)
        NULLIFY(aero_flx(js)%aero_field_cv)
      IF (ASSOCIATED(aero_flx(js)%cloud_field_ls_aer)) &
        DEALLOCATE(aero_flx(js)%cloud_field_ls_aer)
        NULLIFY(aero_flx(js)%cloud_field_ls_aer)
      IF (ASSOCIATED(aero_flx(js)%cloud_field_cv_aer)) &
        DEALLOCATE(aero_flx(js)%cloud_field_cv_aer)
        NULLIFY(aero_flx(js)%cloud_field_cv_aer)

      IF (ASSOCIATED(PR(js)%PROD)) DEALLOCATE(PR(js)%PROD)
      NULLIFY(PR(js)%PROD)
      IF (ASSOCIATED(PR(js)%PROD_IDX)) DEALLOCATE(PR(js)%PROD_IDX)
      NULLIFY(PR(js)%PROD_IDX)
      IF (ASSOCIATED(PR(js)%PROD_CHAR)) DEALLOCATE(PR(js)%PROD_CHAR)
      NULLIFY(PR(js)%PROD_CHAR)
    ENDDO
    
    IF (ASSOCIATED(aero_flx)) DEALLOCATE(aero_flx) ; NULLIFY(aero_flx)
    IF (ASSOCIATED(PR))       DEALLOCATE(PR)       ; NULLIFY(PR)

    IF (ASSOCIATED(te_idt)) THEN
       DEALLOCATE(te_idt) ; NULLIFY(te_idt)
    END IF
    IF (ASSOCIATED(tte_scav)) THEN
       DEALLOCATE(tte_scav) ; NULLIFY(tte_scav)
    END IF

    IF(ASSOCIATED(process)) THEN
       DO j=1,4
          IF (ASSOCIATED(process(j)%frac_eva)) THEN 
             DEALLOCATE(process(j)%frac_eva) ; NULLIFY(process(j)%frac_eva)
          END IF
       ENDDO
       DEALLOCATE(process); NULLIFY(process)
    END IF
       
  END SUBROUTINE scav_free_memory

!=============================================================================  
! coupling with cvtrans for convective scavenging

   SUBROUTINE scav_init_coupling
  
    ! ECHAM5/MESSy
    USE messy_main_mpi_bi,           ONLY: p_parallel_io
    USE messy_main_blather_bi,       ONLY: start_message_bi, end_message_bi &
                                         , error_bi, info_bi
    USE messy_main_channel_error_bi, ONLY: channel_halt
    USE messy_main_tracer_tools_bi,  ONLY: tracer_halt
    USE messy_main_grid_def_mem_bi,  ONLY: nlev, nlevp1

#if defined(ECHAM5) || defined(CESM1)           
    USE messy_main_grid_def_mem_bi,  ONLY:nvclev, vct
#endif
#if defined(COSMO) || defined(MESSYDWARF)
    USE messy_main_grid_def_bi,      ONLY:h_a => hyai, h_b=>hybi
#endif    

     ! MESSy
    USE messy_main_tracer,      ONLY: ON, OFF, CLOUD, get_tracer &
                                    , R_Henry_T0, R_henry_Tdep   &
                                    , R_alpha_T0, R_alpha_Tdep   &
                                    , R_MOLARMASS, get_chemprop
    USE messy_main_tools,       ONLY: strcrack
    USE messy_main_channel,     ONLY: get_channel_object, get_channel_info &
                                    , set_channel_object_restreq, get_attribute
    USE messy_scav_liq,         ONLY: jval_h2o2, jval_hno3, jval_o3
    USE messy_scav_aer,         ONLY: NMAX_AEROMODELS,                     &
                                      aer_int_max, aer_log_max,            &
                                      spec_sol, spec_hetice, spec_idx,     &
                                      spec_mode, spec_idx_comp

    USE messy_scav_inter,       ONLY: get_core_mode_made3, idx_evap_num, made3
    USE messy_scav_l_kpp,       ONLY: lspec => nspec
    USE messy_scav_i_kpp,       ONLY: ispec => nspec
    USE messy_scav_l_kpp,       ONLY: str_field_kpp_l => spc_names
    USE messy_scav_i_kpp,       ONLY: str_field_kpp_i => spc_names

    IMPLICIT NONE
    INTRINSIC :: EPSILON

    INTEGER  :: status, l, jk, js, i, j, k, jt, jn, idx2
    INTEGER  :: jmod
    INTEGER  :: status1, status2

    REAL(dp) :: hypi(nlevp1), sfpress
#if defined(ECHAM5) || defined(CESM1)
    REAL(dp) :: h_a(nvclev), h_b(nvclev)
#endif    
    CHARACTER(LEN=*), PARAMETER::substr='scav_init_coupling'
    CHARACTER(LEN=30) :: string
    CHARACTER(LEN=STRLEN_ULONG) :: dummy_str
    CHARACTER(LEN=26), POINTER, DIMENSION(:)     :: strname => NULL()
    CHARACTER(len=3)  :: strmod
    CHARACTER(len=7)  :: setname
    LOGICAL           :: lnew
    INTEGER           :: naermo, dummy, dummy_idt
    INTEGER           :: evap_max_mode, evap_mode
    INTEGER           :: count_g, count_l, indx
    CHARACTER(LEN=STRLEN_CHANNEL) :: channel_name 
    CHARACTER(LEN=32) :: aermod
    CHARACTER(LEN=52) :: line

    INTRINSIC :: TRIM, SIZE

    CALL start_message_bi(modstr,'COUPLING INITIALIZATION',substr)

    CALL get_channel_info(status, 'cvtrans')
    cvtrans_use = (status == 0)
    IF (.NOT.cvtrans_use) THEN
       CALL info_bi(&
            'required channel of cvtrans module not found', substr)
       CALL info_bi(&
            'cvtrans obviously not used......, ', substr)
       CALL info_bi( '    END - Coupling to CVTRANS          ', substr)
       CALL info_bi('***************************************', substr)
    END IF

    IF (cvtrans_use) THEN
      CALL get_channel_object(status, 'cvtrans', 'cumassf', p3=updraft)
      IF (status == 1) &
        CALL error_bi('channel object for corrected updraft not found', substr)
      CALL get_channel_object(status, 'cvtrans', 'trac_field', p4=trac_field)
      IF (status == 1) &
        CALL error_bi('channel object for tracer field not found', substr)
      
      CALL get_channel_object(status, 'cvtrans', 'kbot', p2=kbot)
      IF (status == 1) &
        CALL error_bi(&
        'channel object for bottom level of convection not found', substr)
      
      CALL get_channel_object(status, 'cvtrans', 'kraintop', p2=kraintop)
      IF (status == 1) CALL error_bi(&
        'channel object for top level of convective precipitation not found' &
        , substr)

      CALL get_channel_object(status, 'cvtrans', 'column', p0=col)
      IF (status == 1) CALL error_bi(&
        'channel object for column logical/real not found', substr)
      IF (col > 0._dp) column = .TRUE.
    ENDIF

    IF (lscav_cv) THEN
       CALL info_bi('Checking for convective parameters ...', substr)

       CALL get_channel_object(status, &
            TRIM(rain_cv(1)), TRIM(rain_cv(2)), p3=precflx_cv)
       IF (status /= 0) CALL error_bi(&
            ' channel object for convective precipitation not found !', substr)
       
       CALL get_channel_object(status, &
            TRIM(ccover(1)), TRIM(ccover(2)), p3=conv_cover)
       IF (status /= 0) CALL error_bi(&
            ' channel object for convective cloud cover not found !', substr)
       
       CALL get_channel_object(status, &
            TRIM(cvprec(1)), TRIM(cvprec(2)), p3=pcvdprec)
       IF (status /= 0) CALL error_bi(&
            ' channel object for freshly formed convective precipitation not found !' &
            , substr)
       
       CALL get_channel_object(status, &
            TRIM(snow_cv(1)), TRIM(snow_cv(2)), p3=snowflx_cv)
       IF (status /= 0) CALL error_bi(&
            ' channel object for convective snowflux not found !', substr)

       CALL get_channel_object(status, &
            TRIM(cvsnow(1)), TRIM(cvsnow(2)), p3=pcvdsnow)
       IF (status /= 0) CALL error_bi(&
            ' channel object for freshly formed convective snow not found !', substr)
    
       CALL get_channel_object(status, &
            TRIM(cvlwc(1)), TRIM(cvlwc(2)), p3=pcvlwc)
       IF (status /= 0) CALL error_bi(&
            ' channel object for cv lwc not found !', substr)

       CALL get_channel_object(status, &
            TRIM(cvrform(1)), TRIM(cvrform(2)), p3=pcvrform)
       IF (status /= 0) CALL error_bi(&
            ' channel object for cv rain formation not found !', substr)

       CALL get_channel_object(status, &
            TRIM(cviwc(1)), TRIM(cviwc(2)), p3=pcviwc)
       IF (status /= 0) CALL error_bi(&
            ' channel object for cv iwc not found !', substr)
       
       CALL get_channel_object(status, &
            TRIM(cvsform(1)), TRIM(cvsform(2)), p3=pcvsform)
       IF (status /= 0) CALL error_bi(&
            ' channel object for cv snow formation not found !', substr)
       
    END IF

    IF (lscav_ls) THEN
       CALL info_bi( 'Checking for large-scale cloud parameters ...', substr)
 
       CALL get_channel_object(status, &
            TRIM(lcover(1)), TRIM(lcover(2)), p3=pclcover)
       IF (status /= 0) CALL error_bi(&
            ' channel object for ls cloud cover not found !', substr)

       CALL get_channel_object(status, &
            TRIM(rcover(1)), TRIM(rcover(2)), p3=prcover)
       IF (status /= 0) CALL error_bi(&
            ' channel object for precipitating ls cloud cover not found !', substr)

       CALL get_channel_object(status, &
            TRIM(ratep(1)), TRIM(ratep(2)), p3=pmratep)
       IF (status /= 0) CALL error_bi(&
            ' channel object for ls rain formation rate not found !', substr)
       
       CALL get_channel_object(status, &
            TRIM(prec(1)), TRIM(prec(2)), p3=pfprec)
       IF (status /= 0) CALL error_bi(&
            ' channel object for ls rain not found !', substr)

       CALL get_channel_object(status, &
            TRIM(ratesi(1)), TRIM(ratesi(2)), p3=pmratesi)
       IF (status /= 0) CALL error_bi(&
            ' channel object for ls snow formation rate not found !', substr)

       CALL get_channel_object(status, &
            TRIM(fsi(1)), TRIM(fsi(2)), p3=pfsi)
       IF (status /= 0) CALL error_bi(&
            ' channel object for ls snow not found !', substr)

       CALL get_channel_object(status, &
            TRIM(lwc(1)), TRIM(lwc(2)), p3=pmlwc)
       IF (status /= 0) CALL error_bi(&
            ' channel object for cloud water content not found !', substr)

       CALL get_channel_object(status, &
            TRIM(iwc(1)), TRIM(iwc(2)), p3=pmiwc)
       IF (status /= 0) CALL error_bi(&
            ' channel object for cloud ice content not found !', substr)

       CALL get_channel_object(status, &
            TRIM(imelt(1)), TRIM(imelt(2)), p3=pimelt)
       IF (status /= 0) CALL error_bi(&
            ' channel object for melting of frozen precip not found !', substr)

       CALL get_channel_object(status, &
            TRIM(isedi(1)), TRIM(isedi(2)), p3=pisedi)
       IF (status /= 0) CALL error_bi(&
            ' channel object for large-scale ice sedimentation not found !', substr)
    END IF
!=================================================      

    CALL info_bi('Checking for aerosol parameters ...', substr)
    DO js=1,nsets
    CALL info_bi("***********************************", substr)
      SELECT CASE (nsets)
      CASE(1)
        IF (lscav_gp) THEN
          ntrac  =  ntrac_gp
          ti     => ti_gp
          setname = GPTRSTR
          strmod = '_gp'
          aermod = aermod_gp
          CALL info_bi( '... treating tracer set '//GPTRSTR, substr)
        ENDIF
        IF (lscav_lg) THEN
          ntrac  =  ntrac_lg
          ti     => ti_lg
          setname = LGTRSTR
          strmod = '_lg'
          aermod = aermod_lg
          CALL info_bi( '... treating tracer set '//LGTRSTR, substr)
        ENDIF
      CASE(2)
        IF (js == 1)  THEN
          ntrac  =  ntrac_gp
          ti     => ti_gp
          setname = GPTRSTR
          strmod = '_gp'
          aermod = aermod_gp
          CALL info_bi( '... treating tracer set '//GPTRSTR, substr)
        ENDIF
        IF (js == 2)  THEN
          ntrac  =  ntrac_lg
          ti     => ti_lg
          setname = LGTRSTR
          strmod = '_lg'
          aermod = aermod_lg
          CALL info_bi( '... treating tracer set '//LGTRSTR, substr)
        ENDIF
      END SELECT
      
      DO jt=1,ntrac
        ! check for running aerosol models
        IF (TRIM(ti(jt)%tp%meta%cask_s(S_aerosol_model)) == '') CYCLE
        IF (.NOT. attr(js)%log_att(jt,laerosol)) CYCLE
        lnew = .TRUE.
        ! CHECK IF AEROSOL MODEL ALREADY IN LIST
        DO naermo = 1, aermodel(js)%aeromodnum
          IF ( TRIM(ti(jt)%tp%meta%cask_s(S_aerosol_model)) &
            == aermodel(js)%aermodname(naermo) ) THEN
            lnew = .FALSE.
            EXIT
          END IF
        END DO
        ! NEW ENTRY IN LIST
        IF (lnew) THEN
          aermodel(js)%aeromodnum = aermodel(js)%aeromodnum + 1
          IF (aermodel(js)%aeromodnum > NMAX_AEROMODELS) &
            CALL error_bi( &
            'RECOMPILATION WITH INCREASED NMAX_AEROMODELS REQUIRED', substr)
          aermodel(js)%aermodname(aermodel(js)%aeromodnum)   = &
            TRIM(ti(jt)%tp%meta%cask_s(S_aerosol_model))
          aermodel(js)%aermodmethod(aermodel(js)%aeromodnum) = &
            ti(jt)%tp%meta%cask_i(I_aerosol_method)
          naermo = aermodel(js)%aeromodnum
        END IF
        attr(js)%int_att(jt,aer_mod_num) = naermo
      ENDDO

      CALL info_bi("... Running aerosol submodels ...", substr)
      IF (p_parallel_io) THEN
        WRITE(*,*) "... Number of aerosol submodels: ",aermodel(js)%aeromodnum
        WRITE(*,*) "... and their respective names: ",&
        aermodel(js)%aermodname(1:aermodel(js)%aeromodnum)
      END IF
      evap_max_mode = 0
      nucl_mode = 0
      SELECT CASE (TRIM(aermod))
      CASE('m7','M7','ptrac','PTRAC','gmxe','GMXE')
        evap_mode = 4
        nucl_mode = 1
      CASE('made','MADE')
        evap_mode = 2
! op_ck_20140514+
      CASE('made3','MADE3')
        ALLOCATE(made3(1))
        made3(1)%ks        = made3params%ks
        made3(1)%km        = made3params%km
        made3(1)%ki        = made3params%ki
        made3(1)%as        = made3params%as
        made3(1)%am        = made3params%am
        made3(1)%ai        = made3params%ai
        made3(1)%cs        = made3params%cs
        made3(1)%cm        = made3params%cm
        made3(1)%ci        = made3params%ci
        made3(1)%nmod      = made3params%nmod
        made3(1)%i_mode    = made3params%i_mode
        made3(1)%csubname  = made3params%csubname
        made3(1)%evap_mode = made3params%evap_mode
        made3(1)%n_resmod  = made3params%n_resmod
        evap_mode = made3(1)%evap_mode
! op_ck_20140514-
      CASE DEFAULT
        evap_mode = 0
      END SELECT
      
      indx = aermodel(js)%aeromodnum
      ALLOCATE(aermodel(js)%aer_mass_attr(indx))
      ALLOCATE(aermodel(js)%aer_num_attr(indx))
      ALLOCATE(aermodel(js)%aer_input(indx))

      DO jn = 1,aermodel(js)%aeromodnum
        aermodel(js)%aer_mass_attr(jn)%number = 0
        aermodel(js)%aer_num_attr(jn)%number = 0
      ENDDO
      
      DO jt=1,ntrac
        IF (attr(js)%log_att(jt,lwetdep) .AND. &
          attr(js)%log_att(jt,laerosol) ) THEN 
          DO jn = 1,aermodel(js)%aeromodnum
            IF (attr(js)%int_att(jt,aer_mod_num) == jn) THEN 
              IF (attr(js)%int_att(jt, lquantity) == 1) THEN
                aermodel(js)%aer_mass_attr(jn)%number = &
                  aermodel(js)%aer_mass_attr(jn)%number + 1
              ENDIF
              IF (attr(js)%int_att(jt, lquantity) == 2) THEN
                aermodel(js)%aer_num_attr(jn)%number = &
                  aermodel(js)%aer_num_attr(jn)%number + 1
              ENDIF
            ENDIF
          ENDDO
        ENDIF
      ENDDO
      
      DO jn = 1,aermodel(js)%aeromodnum 
        indx = aermodel(js)%aer_mass_attr(jn)%number
        ALLOCATE(aermodel(js)%aer_mass_attr(jn)%mass_int_att(indx,aer_int_max))
        ALLOCATE(aermodel(js)%aer_mass_attr(jn)%mass_log_att(indx,aer_log_max))
        indx = aermodel(js)%aer_num_attr(jn)%number 
        ALLOCATE(aermodel(js)%aer_num_attr(jn)%num_int_att(indx,aer_int_max))
        ALLOCATE(aermodel(js)%aer_num_attr(jn)%num_log_att(indx,aer_log_max))
      ENDDO
      
      DO jn = 1,aermodel(js)%aeromodnum 
        aermodel(js)%aer_mass_attr(jn)%mass_log_att(:,spec_sol) = .FALSE.
        aermodel(js)%aer_num_attr(jn)%num_log_att(:,spec_sol)   = .FALSE.
! op_ck_20130114+
        aermodel(js)%aer_mass_attr(jn)%mass_log_att(:,spec_hetice) = .FALSE.
        aermodel(js)%aer_num_attr(jn)%num_log_att(:,spec_hetice)   = .FALSE.
! op_ck_20130114-
        
        i = 1
        j = 1
        
        DO jt=1,ntrac
          IF (attr(js)%log_att(jt,lwetdep) .AND. &
            attr(js)%log_att(jt,laerosol) ) THEN 
            IF (attr(js)%int_att(jt,aer_mod_num) == jn) THEN 
              IF (attr(js)%int_att(jt,lquantity) == 1) THEN
                aermodel(js)%aer_mass_attr(jn)%mass_int_att(i,spec_idx)  = jt
                aermodel(js)%aer_mass_attr(jn)%mass_int_att(i,spec_mode) = &
                  attr(js)%int_att(jt,tmode)
                aermodel(js)%aer_mass_attr(jn)%mass_int_att(i,spec_idx_comp) = &
                  attr(js)%int_att(jt,aerosol_index)
                aermodel(js)%aer_mass_attr(jn)%mass_log_att(i, spec_sol)     = &
                  (ti(jt)%tp%meta%cask_i(I_aerosol_sol) == ON)
! op_ck_20130114+
                aermodel(js)%aer_mass_attr(jn)%mass_log_att(i, spec_hetice)  = &
                  (ti(jt)%tp%meta%cask_i(I_aerosol_hetice) == ON)
! op_ck_20130114-
                i = i + 1 
              ENDIF
              IF (attr(js)%int_att(jt,lquantity) == 2) THEN
                aermodel(js)%aer_num_attr(jn)%num_int_att(j,spec_idx)    = jt
                aermodel(js)%aer_num_attr(jn)%num_int_att(j,spec_mode)   = &
                  attr(js)%int_att(jt,tmode)
                aermodel(js)%aer_num_attr(jn)%num_int_att(j,spec_idx_comp) = & 
                  attr(js)%int_att(jt,aerosol_index)
                aermodel(js)%aer_num_attr(jn)%num_log_att(j, spec_sol)     = &
                  (ti(jt)%tp%meta%cask_i(I_aerosol_sol) == ON)
! op_ck_20130114+
                aermodel(js)%aer_num_attr(jn)%num_log_att(j, spec_hetice)  = &
                  (ti(jt)%tp%meta%cask_i(I_aerosol_hetice) == ON)
! op_ck_20130114-
                j = j + 1
              ENDIF
            ENDIF
          ENDIF
        ENDDO
      ENDDO
      
!!$      if (p_parallel_io) then
!!$        print*, "***** aerosol_species indices * START ******"
!!$        do jn = 1,aermodel(js)%aeromodnum
!!$          print*, "***** MASSES *****"
!!$          do jt=1, aermodel(js)%aer_mass_attr(jn)%number
!!$            print*, jn, jt, &
!!$              aermodel(js)%aer_mass_attr(jn)%mass_int_att(jt,spec_idx), &
!!$              ti(aermodel(js)%aer_mass_attr(jn)%mass_int_att(jt,spec_idx))%tp%ident%fullname, &
!!$              aermodel(js)%aer_mass_attr(jn)%mass_int_att(jt,spec_idx_comp)     
!!$          enddo
!!$          print*, "***** NUMBERS *****"
!!$          do jt=1, aermodel(js)%aer_num_attr(jn)%number
!!$            print*, jn, jt, &
!!$              aermodel(js)%aer_num_attr(jn)%num_int_att(jt,spec_idx), &
!!$              ti(aermodel(js)%aer_num_attr(jn)%num_int_att(jt,spec_idx))%tp%ident%fullname, &
!!$              aermodel(js)%aer_num_attr(jn)%num_int_att(jt,spec_idx_comp)     
!!$          enddo
!!$        enddo
!!$        print*, "***** aerosol_species indices * END ******"
!!$      endif

      DO jt=1, ntrac
        IF ((ti(jt)%tp%meta%cask_i(I_aerosol_sol)==ON)) &
          evap_max_mode = MAX(attr(js)%int_att(jt,tmode), evap_max_mode)
      ENDDO
      
      evap_max_mode = MIN(evap_max_mode, evap_mode)
      
! op_ck_20140529+
!!$      attr(js)%int_att(:,evap_trac_idt)  = 0
!!$      attr(js)%int_att(:,evap_trac_idt2) = 0
      attr(js)%evap_trac(:,:) = 0

      SELECT CASE (TRIM(aermod))

! op_mr_20140902+
      CASE ('made','MADE')

        DO jt = 1, ntrac
           
          ! Aerosol mass (and number)
          IF ( attr(js)%log_att(jt,laerosol) .AND. &
               attr(js)%log_att(jt,lwetdep) .AND. &
               (ti(jt)%tp%meta%cask_i(I_aerosol_mode) <= evap_max_mode) ) THEN
            DO j=1,ntrac
              ! tracer jt will be associated with tracer j after evaporation
              ! if j has the same basename and is in evap_max_mode
              IF ( (TRIM(ti(jt)%tp%ident%basename) ==    &
                   TRIM(ti(j)%tp%ident%basename)).AND.  &
                   (ti(j)%tp%meta%cask_i(I_aerosol_mode) == &
                   evap_max_mode) )  THEN
                ! n_resmod = 1 if MADE is used as `aermod'
                attr(js)%evap_trac(n_resmod,jt) = j
                EXIT
              ENDIF
            ENDDO
          END IF
          ! if no other tracer was found, tracer jt will be associated with
          ! itself after evaporation
          ! n_resmod = 1 if MADE is used as `aermod'
          IF  (attr(js)%evap_trac(n_resmod,jt) == 0) &
               attr(js)%evap_trac(n_resmod,jt) = jt

          ! Aerosol number
          ! Assuming that the aerosol develops from particles of the
          ! largest available soluble mode and smaller particles coalesce
          ! with these, only the number of the evap_max_mode remains after 
          ! evaporation, and the numbers of the smaller modes will vanish
          IF ( (attr(js)%int_att(jt,lquantity) == 2) .AND. &
               ( (ti(jt)%tp%meta%cask_i(I_aerosol_mode) < evap_max_mode) .OR. &
                 (ti(jt)%tp%meta%cask_i(I_aerosol_sol)==OFF) ) ) THEN
            attr(js)%evap_trac(n_resmod,jt) = 0
          ENDIF

        END DO

        ! Determine the index of the number tracer of the largest soluble mode
        idx_evap_num(js) = 0
        DO jt = 1,ntrac
          IF ( (attr(js)%log_att(jt,laerosol))       .AND. &
               (attr(js)%int_att(jt,lquantity) == 2) ) THEN
            IF (ti(jt)%tp%meta%cask_i(I_aerosol_mode) == evap_max_mode) THEN
              idx_evap_num(js) = jt
              EXIT
            END IF
          ENDIF
        ENDDO

! op_mr_20140902-
      CASE ('made3','MADE3')

        DO jt = 1, ntrac
          
          ! Consider only wet deposited aerosol tracers
          IF (attr(js)%log_att(jt,laerosol) .AND. &
               attr(js)%log_att(jt,lwetdep)) THEN
            DO j = 1, ntrac
              ! tracer jt will be associated with tracer j after evaporation if
              ! j has the same basename and is in one of MADE3's residual modes
              IF ((TRIM(ti(jt)%tp%ident%basename) ==                           &
                   TRIM(ti(j)%tp%ident%basename)) .AND.                        &
                   ((ti(j)%tp%meta%cask_i(I_aerosol_mode).EQ.made3(1)%as) .OR. &
                   (ti(j)%tp%meta%cask_i(I_aerosol_mode).EQ.made3(1)%am) .OR.  &
                   (ti(j)%tp%meta%cask_i(I_aerosol_mode).EQ.made3(1)%cs) .OR.  &
                   (ti(j)%tp%meta%cask_i(I_aerosol_mode).EQ.made3(1)%cm))) THEN
                jmod = get_core_mode_made3(ti(j)%tp%meta%cask_i(I_aerosol_mode))
                attr(js)%evap_trac(jmod,jt) = j
              END IF
            END DO
          END IF

        END DO

        idx_evap_num(js) = -1

      CASE DEFAULT
! op_ck_20140529-

      DO jt=1,ntrac
        IF (attr(js)%log_att(jt,laerosol) .AND. &
          attr(js)%log_att(jt,lwetdep) ) THEN
          IF ( (ti(jt)%tp%meta%cask_i(I_aerosol_mode)) <= evap_max_mode) THEN
            
!!$.or. &
!!$             ( (ti(jt)%tp%meta%cask_i(I_aerosol_sol) == OFF) .and.          &
!!$               (TRIM(ti(jt)%tp%meta%cask_s(S_aerosol_model)) == 'made') ) ) then
            DO k=evap_max_mode,1,-1
              DO j=1,ntrac
                ! tracer jt will be associated with tracer j after evaporation
                ! if j is has the same basename and is of in evap_max_mode, 
                ! else in the next lowest mode, until it finds itself
                IF ( (TRIM(ti(jt)%tp%ident%basename) ==    &
                  TRIM(ti(j)%tp%ident%basename)).AND.  &
                  (ti(j)%tp%meta%cask_i(I_aerosol_mode) == k ))  THEN
! op_ck_20140529+
!!$                  attr(js)%int_att(jt,evap_trac_idt) = j
                  attr(js)%evap_trac(n_resmod,jt) = j
! op_ck_20140529-
                  EXIT
                ENDIF
              ENDDO
! op_ck_20140529+
!!$              IF  (attr(js)%int_att(jt,evap_trac_idt) > 0) EXIT
              IF  (attr(js)%evap_trac(n_resmod,jt) > 0) EXIT
! op_ck_20140529-
            ENDDO

            DO k=evap_max_mode-1,1,-1
              DO j=1,ntrac
                ! tracer jt will be associated with tracer j after evaporation
                ! if j is has the same basename and is of in evap_max_mode, 
                ! else in the next lowest mode, until it finds itself
                IF ( (TRIM(ti(jt)%tp%ident%basename) ==              &
                  TRIM(ti(j)%tp%ident%basename))               .AND. &
                  (ti(j)%tp%meta%cask_i(I_aerosol_mode) == k ) .AND. &
! op_ck_20140529+
!!$                  (attr(js)%int_att(jt,evap_trac_idt) /= j) )  THEN
!!$                  attr(js)%int_att(jt,evap_trac_idt2) = j
                  (attr(js)%evap_trac(n_resmod,jt) /= j) )  THEN
                  attr(js)%evap_trac(n_resmod-1,jt) = j
! op_ck_20140529-
                  EXIT
                ENDIF
              ENDDO
! op_ck_20140529+
!!$              IF  (attr(js)%int_att(jt,evap_trac_idt2) > 0) EXIT
              IF  (attr(js)%evap_trac(n_resmod-1,jt) > 0) EXIT
! op_ck_20140529-
            ENDDO

          ENDIF
          ! check for insoluble tracers to be placed in the largest 
          ! available soluble mode
! op_ck_20140529+
!!$          IF (attr(js)%int_att(jt,evap_trac_idt) == 0) THEN
          IF (attr(js)%evap_trac(n_resmod,jt) == 0) THEN
! op_ck_20140529-
            IF ( (ti(jt)%tp%meta%cask_i(I_aerosol_sol) == OFF) .AND.        &
              (TRIM(ti(jt)%tp%meta%cask_s(S_aerosol_model)) /= 'made') ) THEN
              DO k=evap_max_mode,1,-1
                DO j=1,ntrac
                  ! tracer jt will be associated with tracer j after evaporation
                  ! if j is has the same basename and is of in evap_max_mode, 
                  ! else in the next lowest mode, until it finds itself
                  IF ( (TRIM(ti(jt)%tp%ident%basename) ==    &
                    TRIM(ti(j)%tp%ident%basename)).AND.  &
                    (ti(j)%tp%meta%cask_i(I_aerosol_mode) == k ))  THEN
! op_ck_20140529+
!!$                    attr(js)%int_att(jt,evap_trac_idt) = j
                    attr(js)%evap_trac(n_resmod,jt) = j
! op_ck_20140529-
                    EXIT
                  ENDIF
                ENDDO
! op_ck_20140529+
!!$                IF  (attr(js)%int_att(jt,evap_trac_idt) > 0) EXIT
                IF  (attr(js)%evap_trac(n_resmod,jt) > 0) EXIT
! op_ck_20140529-
              ENDDO
            ENDIF
          ENDIF

! op_ck_20140529+
!!$          IF (attr(js)%int_att(jt,evap_trac_idt2) == 0) THEN
          IF (attr(js)%evap_trac(n_resmod-1,jt) == 0) THEN
! op_ck_20140529-
            IF ( (ti(jt)%tp%meta%cask_i(I_aerosol_sol) == OFF) .AND.        &
              (TRIM(ti(jt)%tp%meta%cask_s(S_aerosol_model)) /= 'made') ) THEN
              DO k=evap_max_mode-1,1,-1
                DO j=1,ntrac
                  ! tracer jt will be associated with tracer j after evaporation
                  ! if j is has the same basename and is of in evap_max_mode, 
                  ! else in the next lowest mode, until it finds itself
                  IF ( (TRIM(ti(jt)%tp%ident%basename) ==               &
                    TRIM(ti(j)%tp%ident%basename))                .AND. &
                    (ti(j)%tp%meta%cask_i(I_aerosol_mode) == k )  .AND. &
! op_ck_20140529+
!!$                    (attr(js)%int_att(jt,evap_trac_idt) /= j) )   THEN
!!$                    attr(js)%int_att(jt,evap_trac_idt2) = j
                    (attr(js)%evap_trac(n_resmod,jt) /= j) )   THEN
                    attr(js)%evap_trac(n_resmod-1,jt) = j
! op_ck_20140529-
                    EXIT
                  ENDIF
                ENDDO
! op_ck_20140529+
!!$                IF  (attr(js)%int_att(jt,evap_trac_idt2) > 0) EXIT
                IF  (attr(js)%evap_trac(n_resmod-1,jt) > 0) EXIT
! op_ck_20140529-
              ENDDO
            ENDIF
          ENDIF
! op_ck_20140529+
!!$          IF  (attr(js)%int_att(jt,evap_trac_idt2) == 0) &
!!$            attr(js)%int_att(jt,evap_trac_idt2) = jt
          IF  (attr(js)%evap_trac(n_resmod-1,jt) == 0) &
            attr(js)%evap_trac(n_resmod-1,jt) = jt
! op_ck_20140529-
          
!!! WARNING: treating aerosol numbers differently !!!
          ! assuming that the aerosol develops from particles of the
          ! largest available soluble mode and smaller particles coalesce
          ! with these only the number of the evap_max_mode remain after 
          ! evaporation, and the numbers of the smaller modes will vanish
          IF (attr(js)%int_att(jt,lquantity) == 2) THEN
            DO k=evap_max_mode,1,-1
              IF ( (ti(jt)%tp%meta%cask_i(I_aerosol_mode) < evap_max_mode - 2 )) THEN!.OR.&
!                 ( (ti(jt)%tp%meta%cask_i(I_aerosol_sol)==OFF)            .AND.  &
!                   (ti(jt)%tp%meta%cask_i(I_aerosol_mode) < 7) ) )               THEN
!!$              IF ( (ti(jt)%tp%meta%cask_i(I_aerosol_mode) < evap_max_mode - 1 ) .OR.&
!!$                 ( (ti(jt)%tp%meta%cask_i(I_aerosol_sol)==OFF)            .AND.  &
!!$                   (ti(jt)%tp%meta%cask_i(I_aerosol_mode) < 7) ) )               THEN
! op_ck_20140529+
!!$                attr(js)%int_att(jt,evap_trac_idt) = 0
!!$                attr(js)%int_att(jt,evap_trac_idt2) = 0
                attr(js)%evap_trac(n_resmod-1:n_resmod,jt) = 0
! op_ck_20140529-
              ENDIF
            ENDDO
          ENDIF
        ENDIF
      ENDDO
      
      ! determine the index of the number tracer of the largest soluble mode
      idx_evap_num(js) = 0
      DO jt = 1,ntrac
        IF ( (attr(js)%log_att(jt,laerosol))       .AND. &
          (attr(js)%int_att(jt,lquantity) == 2) ) THEN
          DO k = evap_max_mode,1,-1
            IF (ti(jt)%tp%meta%cask_i(I_aerosol_mode) ==  k) THEN
               idx_evap_num(js) = jt
               EXIT
            END IF
          END DO
        ENDIF
      ENDDO
!!$      if (idx_evap_num(js) /= 0) then
!!$        if (p_parallel_io) &
!!$          print*, "The number is: ", idx_evap_num(js), " for tracer ", &
!!$          ti(idx_evap_num(js))%tp%ident%fullname
!!$      else
!!$        if (p_parallel_io) print*, "No aerosol number tracer available"
!!$      endif
      
!!$      if (p_parallel_io) then
!!$        do jt=1,ntrac
!!$          if ( (attr(js)%log_att(jt,laerosol)) ) then
!!$            if (attr(js)%int_att(jt,evap_trac_idt) > 0) then
!!$              print*, "Tracer ",TRIM(ti(jt)%tp%ident%fullname), &
!!$                " evaporates into Tracer ",                     &
!!$                TRIM(ti(attr(js)%int_att(jt,evap_trac_idt))%tp%ident%fullname)
!!$            else
!!$              print*, "Tracer ",TRIM(ti(jt)%tp%ident%fullname), " vanishes!"
!!$            endif
!!$          endif
!!$        enddo
!!$      endif
      
! op_ck_20140529+
      END SELECT
! op_ck_20140529-

! op_ck_20141112+
      ! Testing residual mode assignment 
      IF (p_parallel_io) THEN
        line = ' RESIDUAL1'
        IF (n_resmod .GT. 1) THEN
          DO j = 2, n_resmod
            WRITE (line,'(A,I1)') TRIM(line)//'    RESIDUAL', j
          END DO
        END IF
        WRITE (*,*)
        WRITE (*,*) ' Aerosol residual assignment upon evaporation '
        WRITE (*,*) '------------------------------------------------------------------'
        WRITE (*,*) '    ORIGIN    --->', line
        WRITE (*,*) '------------------------------------------------------------------'
        DO jt=1,ntrac
          IF (attr(js)%log_att(jt,laerosol) .and. &
               attr(js)%log_att(jt,lwetdep)) THEN
!!$            IF (i_evap == 2) THEN
            line = ''
            DO j = 1, n_resmod
              k = attr(js)%evap_trac(j,jt)
              IF (k .EQ. 0) THEN
                WRITE(line,*) TRIM(line)//'             '
              ELSE
                WRITE(line,'(A,A12)') TRIM(line)//' ' &
                     , TRIM(ti_gp(k)%tp%ident%fullname)
              END IF
            END DO
            WRITE(*,'(A18,A52)') &
                 ' '//TRIM(ti_gp(jt)%tp%ident%fullname)//' --->', line
!!$            ELSE
!!$              j = jt
!!$              k = jt
!!$            END IF
          END IF
        END DO
        WRITE (*,*) '------------------------------------------------------------------'
        k = idx_evap_num(js)
        IF (k .GT. 0) THEN
          WRITE (*,*) ' idx_evap_num ---> ', TRIM(ti_gp(k)%tp%ident%fullname)
        ELSE
          WRITE (*,'(A,I2,A)') ' idx_evap_num ---> not set (= ', k, ')'
        END IF
        WRITE (*,*) '------------------------------------------------------------------'
      END IF
! op_ck_20141112-

      DO naermo = 1,aermodel(js)%aeromodnum
        channel_name = TRIM(aermodel(js)%aermodname(naermo))//strmod
        CALL info_bi(&
          'Checking for aerosol parameters from: '//channel_name//' ', substr) 
        CALL get_channel_info(status, channel_name)
        IF (status /= 0) THEN
           CALL error_bi(&
                'requested aerosol model '//TRIM(channel_name)//' not running!'&
                  , substr)
!!$          CALL info_bi(&
!!$            ' aerosol channel not found for '//channel_name//' !', substr)
!!$          CALL info_bi(' using dummy radius for aerosol scavenging !', substr)
!!$          cpl_aerosol = 1   
!!$          CALL info_bi(&
!!$            ' No aerosol properties available from the chosen channel!', substr)
!!$          CALL info_bi(&
!!$            ' Coupling of scavenging of gases and aerosols reset to level 1' &
!!$            , substr )
!!$          aermodel(js)%aer_input(naermo)%lmode = 7
!!$          CALL info_bi(&
!!$            ' Assuming seven aerosol modes only as in M7 (mode 4 = coarse soluble)', substr)
        END IF
        CALL get_channel_object(status, channel_name,&
          'sigma', p1=aermodel(js)%aer_input(naermo)%aersigma)
        IF (status /= 0) THEN
          CALL info_bi(' sigma channel object not found !', substr)
        ELSE
          aermodel(js)%aer_input(naermo)%lmode = &
            SIZE(aermodel(js)%aer_input(naermo)%aersigma)
        ENDIF

!        IF ((aermodel(js)%aer_input(naermo)%lmode < 4) .AND. l_trac_loc) &
!          CALL error_bi( ' at least 4 modes required', substr)
        
        CALL get_channel_object(status, channel_name,&
          'wetradius', p4=aermodel(js)%aer_input(naermo)%inputrad)
        IF (status /= 0) &
          CALL info_bi(' radius channel object not found !', substr)

        ALLOCATE(aermodel(js)%aer_input(naermo)%cmr2mmr(aermodel(js)%aer_input(naermo)%lmode))
        aermodel(js)%aer_input(naermo)%cmr2mmr(:) = 1._dp
        DO l=1,aermodel(js)%aer_input(naermo)%lmode
! op_pj_20171108: bug-fix to avoid LOG(0) for inactive modes ...
          IF (ASSOCIATED(aermodel(js)%aer_input(naermo)%aersigma)) THEN
             IF (aermodel(js)%aer_input(naermo)%aersigma(l) > EPSILON(0.0)) &
            aermodel(js)%aer_input(naermo)%cmr2mmr(l) =            &
!               EXP(3.5_dp *                                        &
               EXP(3._dp *                                        &
               (LOG(aermodel(js)%aer_input(naermo)%aersigma(l)))**2._dp)
         END IF
         !   aermodel(js)%aer_input(naermo)%cmr2mmr(l) =            &
         !      exp(1.5_dp *                                        &
         !      (log(aermodel(js)%aer_input(naermo)%aersigma(l)))**2._dp)
        ENDDO

! op_ck_20130116+
        CALL get_channel_object(status, channel_name, &
             'philfrac', p4=aermodel(js)%aer_input(naermo)%inputpf)
        IF (status /= 0 .AND. (TRIM(channel_name) == 'made_gp' .OR. &
             TRIM(channel_name) == 'made3_gp')) &
             CALL info_bi(' philfrac channel object not found ! ' &
             // 'Scavenging routines will be applied to total ' &
             // 'number for all modes with number tracer attribute ' &
             // 'I_SCAV=ON' , substr)
! op_ck_20130116-

        CALL info_bi('Checking for activation parameters ...', substr)


        aermodel(js)%aer_input(naermo)%lcalc_nucl_aer = .TRUE.

        CALL get_attribute(status, TRIM(nfrac_nuc(1)), &
          TRIM(nfrac_nuc(2)), 'aermod', c=dummy_str)
        
        IF (TRIM(channel_name) == TRIM(dummy_str) ) THEN
          aermodel(js)%aer_input(naermo)%lcalc_nucl_aer = .FALSE.
          CALL get_channel_object(status1, &
            TRIM(nfrac_nuc(1)), TRIM(nfrac_nuc(2)), &
            p4=aermodel(js)%aer_input(naermo)%nfracnuc)
          IF (status1 /= 0) CALL info_bi(&
            ' channel object for activated number fraction not found !', substr)
          CALL get_channel_object(status2, &
            TRIM(mfrac_nuc(1)), TRIM(mfrac_nuc(2)), &
            p4=aermodel(js)%aer_input(naermo)%mfracnuc)
          IF (status2 /= 0) CALL info_bi(&
            ' channel object for activated mass fraction not found !', substr)
        
          IF ( (status1 /= 0) .or. (status2 /= 0) )THEN
            CALL info_bi(&
              'No channel objects found for activated fractions, use internal parameterisation!',&
              substr)
            aermodel(js)%aer_input(naermo)%lcalc_nucl_aer = .TRUE.
          ENDIF
        ELSE
          CALL info_bi(&
            'Activation routine and aerosol submodel do not match, use internal parameterisation!',&
            substr)
          aermodel(js)%aer_input(naermo)%lcalc_nucl_aer = .TRUE.
        ENDIF

      ENDDO

      max_mode = 0
      DO naermo = 1,aermodel(js)%aeromodnum
        max_mode = MAX(max_mode,aermodel(js)%aer_input(naermo)%lmode)
      ENDDO

!-------------------------------------------------------------------------
!!$! op_cf_20120416+
!!$! get tracer for HNO3-Tendency for AIRTRAC
!!$      IF (LTEND_HNO3) THEN
!!$         CALL get_tracer(status,GPTRSTR,'HNO3',idx=idt_hno3)
!!$         CALL tracer_halt(substr, status)
!!$      ENDIF
!!$! op_cf_20120416-

! new mapping of KPP species and tracers for SCAV

!      gas-liquid interaction

      CALL info_bi(' Setting up: Gas - Liquid - Interactions!', substr)

      ALLOCATE(kpp_l_idx(js)%gas_spec(lspec_gas, gas_max))
      ALLOCATE(kpp_l_idx(js)%gas_attr(lspec_gas, gatt_max))
      ALLOCATE(kpp_l_idx(js)%liq_spec(lspec_liq, liq_max))
      ALLOCATE(kpp_l_idx(js)%liq_attr(lspec_liq, latt_max))
      count_g = 0
      count_l = 0
      DO jt=1,lspec
         dummy = 0
         IF (ASSOCIATED(strname)) DEALLOCATE (strname)
         NULLIFY(strname)
         CALL strcrack(STR_FIELD_KPP_L(jt),'_', strname, dummy)
         IF (TRIM(strname(1)) == 'Prod') CYCLE
         IF (dummy==1) THEN
            CALL get_tracer(status, setname, TRIM(strname(1)),    &
                 idx=dummy_idt)
            IF (status /= 0) dummy_idt = 0
            count_g = count_g + 1
            kpp_l_idx(js)%gas_spec(count_g,gas_idx)  = jt
            kpp_l_idx(js)%gas_spec(count_g,gas2trac) = dummy_idt
            IF (status == 0) THEN
               CALL get_tracer(status, setname, dummy_idt, R_molarmass&
                    , kpp_l_idx(js)%gas_attr(count_g,gas_mw))
               CALL tracer_halt(substr,status)
               IF (.NOT. loverwrite_Henry) THEN
                  CALL get_tracer(status, setname,dummy_idt,R_Henry_T0, HS0(jt))
                  CALL tracer_halt(substr,status)
                  IF (HS0(jt) < 0._dp) CALL error_bi(&
                       'Henry_T0 coefficient missing for '//TRIM(strname(1))&
                       , substr)
                  CALL get_tracer(status, setname,dummy_idt &
                       , R_Henry_Tdep, DHT(jt))
                  CALL tracer_halt(substr,status)
               END IF
               IF (.NOT. loverwrite_alpha) THEN
                  CALL get_tracer(status, setname,dummy_idt &
                       , R_alpha_T0, alpha0(jt))
                  CALL tracer_halt(substr,status)
                  CALL get_tracer(status, setname,dummy_idt &
                       , R_alpha_Tdep, alpha_T(jt))
                  CALL tracer_halt(substr,status)
               END IF
            ELSE
               status = get_chemprop(TRIM(strname(1)),R_molarmass &
                    , kpp_l_idx(js)%gas_attr(count_g,gas_mw) )
               IF (status /= 0) CALL error_bi (&
                    'species '//TRIM(strname(1))//' not part of chemprop'&
                    , substr)
               IF (.NOT. loverwrite_Henry) THEN
                  status = get_chemprop(TRIM(strname(1)), R_Henry_T0, HS0(jt))
                  IF (HS0(jt) < 0._dp .OR. status /= 0) CALL error_bi(&
                       'Henry_T0 coefficient missing for '//TRIM(strname(1))&
                       , substr)
                  status = get_chemprop(TRIM(strname(1)),R_Henry_Tdep, dht(jt))
                  IF (status /=0 )  CALL error_bi (&
                       'species '//TRIM(strname(1))//' read dht chemprop' &
                       , substr)
               END IF
               IF (.NOT. loverwrite_alpha) THEN
                  status = get_chemprop(TRIM(strname(1)),R_alpha_T0, alpha0(jt))
                  IF (status /=0 )  CALL error_bi (&
                       'species '//TRIM(strname(1))//' read alpha0 chemprop'&
                       , substr)
                  status = get_chemprop(TRIM(strname(1)), R_alpha_Tdep &
                       , alpha_T(jt))
                  IF (status /=0 )  CALL error_bi (&
                       'species '//TRIM(strname(1))//' read alpha_T chemprop'&
                       , substr)
               END IF
            END IF
            IF (.NOT. attr(js)%log_att(dummy_idt,lwetdep)) THEN
               SELECT CASE (TRIM(strname(1)))
               CASE ("O2","N2")
                  CYCLE
               CASE DEFAULT
                  CALL error_bi('Species '//TRIM(strname(1))//' is part of the mechanism, '&
                       //'but shall not be scavenged.'&
                       //' This causes mass violations! Change mechanism and rerun!', substr)
               END SELECT
            ENDIF

         ELSE
            CALL get_tracer(status, setname, TRIM(strname(1)),    &
                 subname=TRIM(strname(2)),idx=dummy_idt)
            IF (status /= 0) dummy_idt = 0
            count_l = count_l + 1
            kpp_l_idx(js)%liq_spec(count_l,liq_idx)  = jt
            kpp_l_idx(js)%liq_spec(count_l,liq2trac) = dummy_idt
            IF (status == 0) THEN
               CALL get_tracer(status, setname, dummy_idt, R_molarmass &
                    , kpp_l_idx(js)%liq_attr(count_l,liq_mw))
               CALL tracer_halt(substr,status)
            ELSE
               status = get_chemprop(strname(1),R_molarmass &
                    , kpp_l_idx(js)%liq_attr(count_l,liq_mw) )
               IF (status /= 0) CALL error_bi (&
                    'species '//TRIM(strname(1))//' not part of chemprop'&
                    , substr)
            END IF
         ENDIF
      ENDDO

!--------------------------------

      CALL info_bi(' Setting up: Gas - Ice - Interactions!', substr)
      
      ALLOCATE(kpp_i_idx(js)%gas_spec(ispec_gas, gas_max))
      ALLOCATE(kpp_i_idx(js)%gas_attr(ispec_gas, gatt_max))
      ALLOCATE(kpp_i_idx(js)%ice_spec(ispec_ice, ice_max))
      ALLOCATE(kpp_i_idx(js)%ice_attr(ispec_ice, iatt_max))

      count_g = 0
      count_l = 0
      DO jt=1,ispec
         dummy=0
         IF (ASSOCIATED(strname)) DEALLOCATE (strname)
         NULLIFY(strname)
         CALL strcrack(STR_FIELD_KPP_I(jt),'_', strname, dummy)
         IF (TRIM(strname(1)) == 'Prod') CYCLE
         IF (dummy==1) THEN
            CALL get_tracer(status, setname, TRIM(strname(1)),idx=dummy_idt)
            count_g = count_g + 1
            kpp_i_idx(js)%gas_spec(count_g,gas_idx) = jt
            kpp_i_idx(js)%gas_spec(count_g,gas2trac) = dummy_idt
            IF (status == 0) THEN
               CALL get_tracer(status, setname, dummy_idt, R_molarmass &
                    , kpp_i_idx(js)%gas_attr(count_g,gas_mw))
               CALL tracer_halt(substr, status)
! ju_ak+
!! ice still uses predefined values
! activating this might overwrite the (correct) values for liquied phase
! as in my assessment species associated to one specific index (idx)
!  are not the same for liquid and ice phase (need additional hs0_ice field
! here for correct coupling (which currently is not used!)
!!$          CALL get_tracer(status, setname, dummy_idt, R_molarmass&
!!$           , kpp_i_idx(js)%gas_attr(count_g,gas_mw))
!!$          CALL tracer_halt(substr,status)
!!$          CALL get_tracer(status, setname, TRIM(strname(1)),    &
!!$               dummy_idt, R_Henry_T0, HS0(jt))
!!$          CALL tracer_halt(substr,status)
!!$          CALL get_tracer(status, setname, TRIM(strname(1)),    &
!!$               dummy_idt, R_Henry_Tdep, HSt(jt))
!!$          CALL tracer_halt(substr,status)
!!$          CALL get_tracer(status, setname, TRIM(strname(1)),    &
!!$               dummy_idt, R_alpha_T0, alpha0(jt))
!!$          CALL tracer_halt(substr,status)
!!$          CALL get_tracer(status, setname, TRIM(strname(1)),    &
!!$               dummy_idt, R_alpha_Tdp, alpha_T(jt))
!!$          CALL tracer_halt(substr,status)
! ju_ak-
            ELSE
               status = get_chemprop(TRIM(strname(1)),R_molarmass &
                    , kpp_i_idx(js)%gas_attr(count_g,gas_mw) )
               IF (status /= 0) CALL error_bi (&
                    'species '//TRIM(strname(1))//' not part of chemprop'&
                    , substr)
               ! add henry & Co here, if ice phase no longer hard-coded
            END IF

            IF (.NOT. attr(js)%log_att(dummy_idt,lwetdep)) THEN
               SELECT CASE (TRIM(strname(1)))
               CASE ("O2","N2")
                  CYCLE
               CASE DEFAULT
                  CALL error_bi('Species '//TRIM(strname(1))//' is part of the mechanism, '&
                       //'but shall not be scavenged.'&
                       //' This causes mass violations! Change mechanism and rerun!', substr)
               END SELECT
            ENDIF

         ELSE
            CALL get_tracer(status, setname, TRIM(strname(1)),    &
                 subname=TRIM(strname(2)),idx=dummy_idt)
            IF (status /= 0) dummy_idt = 0
            count_l = count_l + 1
            kpp_i_idx(js)%ice_spec(count_l,ice_idx) = jt
            kpp_i_idx(js)%ice_spec(count_l,ice2trac) = dummy_idt
            IF (status == 0) THEN
               CALL get_tracer(status, setname, dummy_idt,R_molarmass &
                    , kpp_i_idx(js)%ice_attr(count_l,ice_mw))
               CALL tracer_halt(substr, status)
            ELSE
               status = get_chemprop(TRIM(strname(1)),R_molarmass &
                    , kpp_i_idx(js)%ice_attr(count_l,ice_mw) )
               IF (status /= 0) CALL error_bi (&
                    'species '//TRIM(strname(1))//' not part of chemprop'&
                    , substr)
            END IF
         ENDIF
      ENDDO

!----------------------------------

!     aerosol_liquid interaction
       
      CALL info_bi(' Setting up: Aerosol - Liquid - Interactions!', substr)
    
      aer_attr(js)%aer_spec(:,aer2l_idx) = 0
      DO jt=1,aer_count(js)%numbers(c_all)
        string= & 
          TRIM(ti(aer_attr(js)%aer_spec(jt,aer_idx))%tp%ident%basename)//'_l'
        DO j=1,lspec
          IF (TRIM(STR_FIELD_KPP_L(j)) == TRIM(string) ) THEN
            aer_attr(js)%aer_spec(jt,aer2l_idx) = j           
          ENDIF
        ENDDO
      ENDDO
    
!   for specific artificial tracers from the aerosol modules 

      DO jt = 1, aer_count(js)%numbers(c_all)
        IF (aer_attr(js)%aer_spec(jt,aer2l_idx) /= 0) CYCLE
        string= TRIM(ti(aer_attr(js)%aer_spec(jt,aer_idx))%tp%ident%basename)
        SELECT CASE (string)
        CASE('SO4')
          DO j=1,lspec
            IF (TRIM(STR_FIELD_KPP_L(j)) == 'SO4mm_l' ) &
              aer_attr(js)%aer_spec(jt,aer2l_idx) = j
          ENDDO
! op_ck_20150213+
!!$          IF ((TRIM(ti(aer_attr(js)%aer_spec(jt,aer_idx))%tp%ident%fullname)&
!!$            == 'SO4_ns') .OR.                                               &
!!$            (TRIM(ti(aer_attr(js)%aer_spec(jt,aer_idx))%tp%ident%fullname)  &
!!$            == 'SO4_ks') )                                                  &
!!$!qqq       .or. (TRIM(ti(aer_attr(js)%aer_spec(jt,aer_idx))%tp%ident%fullname)  &
!!$!qqq         == 'SO4_cs') ) &
!!$            aer_attr(js)%aer_spec(jt,aer2l_idx) = 0
! op_ck_20150213-
        CASE('NO3')
          DO j=1,lspec
            IF (TRIM(STR_FIELD_KPP_L(j)) == 'NO3m_l' ) &
              aer_attr(js)%aer_spec(jt,aer2l_idx) = j
          ENDDO
        CASE('NH4')
          DO j=1,lspec
            IF (TRIM(STR_FIELD_KPP_L(j)) == 'NH4p_l' ) &
              aer_attr(js)%aer_spec(jt,aer2l_idx) = j
          ENDDO
        CASE('Cl')
          DO j=1,lspec
            IF (TRIM(STR_FIELD_KPP_L(j)) == 'Clm_l' ) &
              aer_attr(js)%aer_spec(jt,aer2l_idx) = j
          ENDDO
        CASE('SO4res')
          DO j=1,lspec
            IF (TRIM(STR_FIELD_KPP_L(j)) == 'SO4mm_l' ) &
              aer_attr(js)%aer_spec(jt,aer2l_idx) = j
          ENDDO
        CASE('NO3mres')
          DO j=1,lspec
            IF (TRIM(STR_FIELD_KPP_L(j)) == 'NO3m_l' ) &
              aer_attr(js)%aer_spec(jt,aer2l_idx) = j
          ENDDO
        CASE('NH4pres')
          DO j=1,lspec
            IF (TRIM(STR_FIELD_KPP_L(j)) == 'NH4p_l' ) &
              aer_attr(js)%aer_spec(jt,aer2l_idx) = j
          ENDDO
        CASE('Hpres')
          DO j=1,lspec
            IF (TRIM(STR_FIELD_KPP_L(j)) == 'Hp_l' ) &
              aer_attr(js)%aer_spec(jt,aer2l_idx) = j
          ENDDO
        CASE('Clmres')
          DO j=1,lspec
            IF (TRIM(STR_FIELD_KPP_L(j)) == 'Clm_l' ) &
              aer_attr(js)%aer_spec(jt,aer2l_idx) = j
          ENDDO
        END SELECT
      ENDDO
     
!!$      if (p_parallel_io) then
!!$        do jt=1,aer_count(js)%numbers(c_all)
!!$          if (aer_attr(js)%aer_spec(jt,aer2l_idx) /= 0 ) then
!!$            print*, "Tracer ", aer_attr(js)%aer_spec(jt,aer_idx),             &
!!$              TRIM(ti(aer_attr(js)%aer_spec(jt,aer_idx))%tp%ident%fullname),  &
!!$              " is used in liquid phase for KPP_SPECIES",                     &
!!$              aer_attr(js)%aer_spec(jt,aer2l_idx),                            &
!!$              (TRIM(STR_FIELD_KPP_L(aer_attr(js)%aer_spec(jt,aer2l_idx))))
!!$          else
!!$            print*, "Tracer ", aer_attr(js)%aer_spec(jt,aer_idx),             &
!!$              TRIM(ti(aer_attr(js)%aer_spec(jt,aer_idx))%tp%ident%fullname),  &
!!$              " is not used in liquid phase chemistry!"
!!$          endif
!!$        enddo
!!$      endif

!--------------------------------------

      CALL info_bi(' Setting up: Aerosol - Ice - Interactions!', substr)
    
      aer_attr(js)%aer_spec(:,aer2i_idx) = 0
      DO jt=1,aer_count(js)%numbers(c_all)
        string= & 
          TRIM(ti(aer_attr(js)%aer_spec(jt,aer_idx))%tp%ident%basename)//'_i'
        DO j=1,ispec
          IF (TRIM(STR_FIELD_KPP_I(j)) == TRIM(string) ) THEN
            aer_attr(js)%aer_spec(jt,aer2i_idx) = j           
          ENDIF
        ENDDO
      ENDDO
    
!   for specific artificial tracers from the aerosol modules 

      DO jt = 1, aer_count(js)%numbers(c_all)
        IF (aer_attr(js)%aer_spec(jt,aer2i_idx) /= 0) CYCLE
        string= TRIM(ti(aer_attr(js)%aer_spec(jt,aer_idx))%tp%ident%basename)
        SELECT CASE (string)
        CASE('SO4')
          DO j=1,ispec
            IF (TRIM(STR_FIELD_KPP_i(j)) == 'SO4mm_i' ) &
              aer_attr(js)%aer_spec(jt,aer2i_idx) = j
          ENDDO
! op_ck_20150213+
!!$          IF ((TRIM(ti(aer_attr(js)%aer_spec(jt,aer_idx))%tp%ident%fullname)&
!!$            == 'SO4_ns') .OR.                                               &
!!$            (TRIM(ti(aer_attr(js)%aer_spec(jt,aer_idx))%tp%ident%fullname)  &
!!$            == 'SO4_ks') )                                                  &
!!$!qqq      .or. (TRIM(ti(aer_attr(js)%aer_spec(jt,aer_idx))%tp%ident%fullname)  &
!!$!qqq        == 'SO4_cs') )                                                  &
!!$            aer_attr(js)%aer_spec(jt,aer2i_idx) = 0
! op_ck_20150213-
        CASE('NO3')
          DO j=1,ispec
            IF (TRIM(STR_FIELD_KPP_I(j)) == 'NO3m_i' ) &
              aer_attr(js)%aer_spec(jt,aer2i_idx) = j
          ENDDO
        CASE('NH4')
          DO j=1,ispec
            IF (TRIM(STR_FIELD_KPP_I(j)) == 'NH4p_i' ) &
              aer_attr(js)%aer_spec(jt,aer2i_idx) = j
          ENDDO
        CASE('Cl')
          DO j=1,ispec
            IF (TRIM(STR_FIELD_KPP_I(j)) == 'Clm_i' ) &
              aer_attr(js)%aer_spec(jt,aer2i_idx) = j
          ENDDO
        CASE('SO4res')
          DO j=1,ispec
            IF (TRIM(STR_FIELD_KPP_I(j)) == 'SO4mm_i' ) &
              aer_attr(js)%aer_spec(jt,aer2i_idx) = j
          ENDDO
        CASE('NO3mres')
          DO j=1,ispec
            IF (TRIM(STR_FIELD_KPP_I(j)) == 'NO3m_i' ) &
              aer_attr(js)%aer_spec(jt,aer2i_idx) = j
          ENDDO
        CASE('NH4pres')
          DO j=1,ispec
            IF (TRIM(STR_FIELD_KPP_I(j)) == 'NH4p_i' ) &
              aer_attr(js)%aer_spec(jt,aer2i_idx) = j
          ENDDO
        CASE('Hpres')
          DO j=1,ispec
            IF (TRIM(STR_FIELD_KPP_I(j)) == 'Hp_i' ) &
              aer_attr(js)%aer_spec(jt,aer2i_idx) = j
          ENDDO
        CASE('Clmres')
          DO j=1,ispec
            IF (TRIM(STR_FIELD_KPP_I(j)) == 'Clm_i' ) &
              aer_attr(js)%aer_spec(jt,aer2i_idx) = j
          ENDDO
        END SELECT
      ENDDO
     
!!$      if (p_parallel_io) then
!!$        do jt=1,aer_count(js)%numbers(c_all)
!!$          if (aer_attr(js)%aer_spec(jt,aer2i_idx) /= 0 ) then
!!$            print*, "Tracer ", aer_attr(js)%aer_spec(jt,aer_idx),           &
!!$              TRIM(ti(aer_attr(js)%aer_spec(jt,aer_idx))%tp%ident%fullname),&
!!$              " is used in ice phase for KPP_SPECIES",                      &
!!$              aer_attr(js)%aer_spec(jt,aer2i_idx),                          &
!!$              (TRIM(STR_FIELD_KPP_I(aer_attr(js)%aer_spec(jt,aer2i_idx))))
!!$          else
!!$            print*, "Tracer ", aer_attr(js)%aer_spec(jt,aer_idx),           &
!!$              TRIM(ti(aer_attr(js)%aer_spec(jt,aer_idx))%tp%ident%fullname),&
!!$              " is not used in ice phase chemistry!"
!!$          endif
!!$        enddo
!!$      endif

!----------------------------------------

!     liquid - evaporation - interaction
      CALL info_bi(' Setting up: Liquid - Evaporation - Interactions!', substr)
      kpp_l_idx(js)%liq_spec(:,liq2evap) = 0
      kpp_l_idx(js)%liq_spec(:,liq2eva2) = 0

      DO j=1,lspec_liq
        if (P_parallel_io) print*, "Searching for: ", &
          TRIM(STR_FIELD_KPP_L(kpp_l_idx(js)%liq_spec(j,liq_idx)))
!!$        DO jt=1,ntrac
!!$          string=TRIM(ti(jt)%tp%ident%basename)//'_l'
!!$          if (P_parallel_io) print*, "Trying to compare: ", TRIM(string), &
!!$            " ", TRIM(ti(jt)%tp%ident%fullname), " ", &
!!$            TRIM(ti(jt)%tp%ident%basename), attr(js)%log_att(jt,laerosol), &
!!$            attr(js)%int_att(jt,tmode)
!!$          DO k=evap_max_mode,1,-1
!!$            IF ( (TRIM(STR_FIELD_KPP_L(kpp_l_idx(js)%liq_spec(j,liq_idx))) &
!!$              == TRIM(string)) .AND. attr(js)%log_att(jt,laerosol) .AND.     &
!!$              (attr(js)%int_att(jt,tmode) == k) ) THEN
!!$              kpp_l_idx(js)%liq_spec(j,liq2evap) = jt
!!$              if (P_parallel_io) print*,"Found ! ", &
!!$                TRIM(ti(jt)%tp%ident%fullname)
!!$              EXIT
!!$            ENDIF
!!$          ENDDO ! mode
!!$          DO k=evap_max_mode-1,1,-1
!!$            IF ( (TRIM(STR_FIELD_KPP_L(kpp_l_idx(js)%liq_spec(j,liq_idx))) &
!!$              == TRIM(string)) .AND. attr(js)%log_att(jt,laerosol) .AND.     &
!!$              (attr(js)%int_att(jt,tmode) == k) ) THEN
!!$              kpp_l_idx(js)%liq_spec(j,liq2eva2) = jt
!!$              EXIT
!!$            ENDIF
!!$          ENDDO ! mode
!!$        ENDDO   ! ntrac


! op_ck_20150128+
        IF (.NOT. ALLOCATED(made3)) THEN
! op_ck_20150128-

        DO k=evap_max_mode,1,-1
          IF (kpp_l_idx(js)%liq_spec(j,liq2evap) /= 0) CYCLE
          DO jt=1,ntrac
            string=TRIM(ti(jt)%tp%ident%basename)//'_l'
            IF ( (TRIM(STR_FIELD_KPP_L(kpp_l_idx(js)%liq_spec(j,liq_idx))) &
              == TRIM(string)) .AND. attr(js)%log_att(jt,laerosol) .AND.     &
              (attr(js)%int_att(jt,tmode) == k) ) THEN
              kpp_l_idx(js)%liq_spec(j,liq2evap) = jt
              EXIT
            END IF
          END DO
        END DO

! op_mr_20140908+
        IF ((TRIM(aermod) /= 'made') .AND. (TRIM(aermod) /= 'MADE')) THEN
! op_mr_20140908-

        DO k=evap_max_mode-1,1,-1
          IF (kpp_l_idx(js)%liq_spec(j,liq2eva2) /= 0) CYCLE
          DO jt=1,ntrac
            string=TRIM(ti(jt)%tp%ident%basename)//'_l'
            IF ( (TRIM(STR_FIELD_KPP_L(kpp_l_idx(js)%liq_spec(j,liq_idx))) &
              == TRIM(string)) .AND. attr(js)%log_att(jt,laerosol) .AND.     &
              (attr(js)%int_att(jt,tmode) == k) ) THEN
              kpp_l_idx(js)%liq_spec(j,liq2eva2) = jt
              EXIT
            END IF
          END DO
        END DO          

! op_mr_20140908+
        END IF
! op_mr_20140908-

! op_ck_20150128+
        END IF
! op_ck_20150128-

        IF (kpp_l_idx(js)%liq_spec(j,liq2evap) == 0 ) THEN
          DO jt=1,ntrac
            string=TRIM(ti(jt)%tp%ident%fullname)//'_l'
            IF ( TRIM(STR_FIELD_KPP_L(kpp_l_idx(js)%liq_spec(j,liq_idx))) &
              == TRIM(string) )   THEN
               kpp_l_idx(js)%liq_spec(j,liq2evap) = jt
               EXIT
            ENDIF
          ENDDO
        ENDIF

! op_ck_20150128+
        IF (.NOT. ALLOCATED(made3)) THEN
! op_ck_20150128-

        ! special cases
        IF (kpp_l_idx(js)%liq_spec(j,liq2evap) /= 0 ) CYCLE
        DO k=evap_max_mode,1,-1
          DO jt=1,ntrac
            IF ( attr(js)%int_att(jt,tmode) == k) THEN
              string=TRIM(ti(jt)%tp%ident%basename)
              IF (kpp_l_idx(js)%liq_spec(j,liq2evap) /= 0 ) CYCLE
              SELECT CASE (TRIM(string))
              CASE('SO4')
                IF (TRIM(STR_FIELD_KPP_L(kpp_l_idx(js)%liq_spec(j,liq_idx))) &
                  == 'SO4mm_l') kpp_l_idx(js)%liq_spec(j,liq2evap) = jt
                IF (TRIM(STR_FIELD_KPP_L(kpp_l_idx(js)%liq_spec(j,liq_idx))) &
                  == 'HSO4m_l') kpp_l_idx(js)%liq_spec(j,liq2evap) = jt
              CASE('NO3')
                IF (TRIM(STR_FIELD_KPP_L(kpp_l_idx(js)%liq_spec(j,liq_idx))) &
                  == 'NO3m_l') kpp_l_idx(js)%liq_spec(j,liq2evap) = jt
              CASE('NH4')
                IF (TRIM(STR_FIELD_KPP_L(kpp_l_idx(js)%liq_spec(j,liq_idx))) &
                  == 'NH4p_l') kpp_l_idx(js)%liq_spec(j,liq2evap) = jt
              CASE('SO4res')
                IF (TRIM(STR_FIELD_KPP_L(kpp_l_idx(js)%liq_spec(j,liq_idx))) &
                  == 'SO4mm_l') kpp_l_idx(js)%liq_spec(j,liq2evap) = jt
                IF (TRIM(STR_FIELD_KPP_L(kpp_l_idx(js)%liq_spec(j,liq_idx))) &
                  == 'HSO4m_l') kpp_l_idx(js)%liq_spec(j,liq2evap) = jt
              CASE('NO3mres')
                IF (TRIM(STR_FIELD_KPP_L(kpp_l_idx(js)%liq_spec(j,liq_idx))) &
                  == 'NO3m_l') kpp_l_idx(js)%liq_spec(j,liq2evap) = jt
              CASE('NH4pres')
                IF (TRIM(STR_FIELD_KPP_L(kpp_l_idx(js)%liq_spec(j,liq_idx))) &
                  == 'NH4p_l') kpp_l_idx(js)%liq_spec(j,liq2evap) = jt
              CASE('Hpres')
                IF (TRIM(STR_FIELD_KPP_L(kpp_l_idx(js)%liq_spec(j,liq_idx))) &
                  == 'Hp_l') kpp_l_idx(js)%liq_spec(j,liq2evap) = jt
              CASE('Clmres')
                IF (TRIM(STR_FIELD_KPP_L(kpp_l_idx(js)%liq_spec(j,liq_idx))) &
                  == 'Clm_l') kpp_l_idx(js)%liq_spec(j,liq2evap) = jt
              END SELECT
            ENDIF
          ENDDO
        ENDDO
      
! op_mr_20140908+
        IF ((TRIM(aermod) /= 'made') .AND. (TRIM(aermod) /= 'MADE')) THEN
! op_mr_20140908-

        IF (kpp_l_idx(js)%liq_spec(j,liq2eva2) /= 0 ) CYCLE
        DO k=evap_max_mode-1,1,-1
          DO jt=1,ntrac
            IF ( attr(js)%int_att(jt,tmode) == k) THEN
              string=TRIM(ti(jt)%tp%ident%basename)
              IF (kpp_l_idx(js)%liq_spec(j,liq2eva2) /= 0 ) CYCLE
              SELECT CASE (TRIM(string))
              CASE('SO4')
                IF (TRIM(STR_FIELD_KPP_L(kpp_l_idx(js)%liq_spec(j,liq_idx))) &
                  == 'SO4mm_l') kpp_l_idx(js)%liq_spec(j,liq2eva2) = jt
                IF (TRIM(STR_FIELD_KPP_L(kpp_l_idx(js)%liq_spec(j,liq_idx))) &
                  == 'HSO4m_l') kpp_l_idx(js)%liq_spec(j,liq2eva2) = jt
              CASE('NO3')
                IF (TRIM(STR_FIELD_KPP_L(kpp_l_idx(js)%liq_spec(j,liq_idx))) &
                  == 'NO3m_l') kpp_l_idx(js)%liq_spec(j,liq2eva2) = jt
              CASE('NH4')
                IF (TRIM(STR_FIELD_KPP_L(kpp_l_idx(js)%liq_spec(j,liq_idx))) &
                  == 'NH4p_l') kpp_l_idx(js)%liq_spec(j,liq2eva2) = jt
              CASE('SO4res')
                IF (TRIM(STR_FIELD_KPP_L(kpp_l_idx(js)%liq_spec(j,liq_idx))) &
                  == 'SO4mm_l') kpp_l_idx(js)%liq_spec(j,liq2eva2) = jt
                IF (TRIM(STR_FIELD_KPP_L(kpp_l_idx(js)%liq_spec(j,liq_idx))) &
                  == 'HSO4m_l') kpp_l_idx(js)%liq_spec(j,liq2eva2) = jt
              CASE('NO3mres')
                IF (TRIM(STR_FIELD_KPP_L(kpp_l_idx(js)%liq_spec(j,liq_idx))) &
                  == 'NO3m_l') kpp_l_idx(js)%liq_spec(j,liq2eva2) = jt
              CASE('NH4pres')
                IF (TRIM(STR_FIELD_KPP_L(kpp_l_idx(js)%liq_spec(j,liq_idx))) &
                  == 'NH4p_l') kpp_l_idx(js)%liq_spec(j,liq2eva2) = jt
              CASE('Hpres')
                IF (TRIM(STR_FIELD_KPP_L(kpp_l_idx(js)%liq_spec(j,liq_idx))) &
                  == 'Hp_l') kpp_l_idx(js)%liq_spec(j,liq2eva2) = jt
              CASE('Clmres')
                IF (TRIM(STR_FIELD_KPP_L(kpp_l_idx(js)%liq_spec(j,liq_idx))) &
                  == 'Clm_l') kpp_l_idx(js)%liq_spec(j,liq2eva2) = jt
              END SELECT
            ENDIF
          ENDDO
        ENDDO

        ! if a second compounds could not be found use the first 
        ! found compound
        IF (kpp_l_idx(js)%liq_spec(j,liq2eva2) /= 0 ) CYCLE
        kpp_l_idx(js)%liq_spec(j,liq2eva2) = kpp_l_idx(js)%liq_spec(j,liq2evap)

! op_mr_20140908+
        END IF
! op_mr_20140908-

! op_ck_20150128+
        END IF
! op_ck_20150128-

      END DO

! op_ck_20150128+
      IF (.NOT. ALLOCATED(made3)) THEN
! op_ck_20150128-
      if (p_parallel_io) then
        do j=1,lspec_liq
          if (kpp_l_idx(js)%liq_spec(j,liq2evap) == 0 ) then
            print*, "This species does not evaporate:",  j,            &
              STR_FIELD_KPP_L(kpp_l_idx(js)%liq_spec(j,liq_idx))
          else
            print*, j, STR_FIELD_KPP_L(kpp_l_idx(js)%liq_spec(j,liq_idx)), &
              kpp_l_idx(js)%liq_spec(j,liq2evap),                        &
              TRIM(ti(kpp_l_idx(js)%liq_spec(j,liq2evap))%tp%ident%fullname)
          end if
        enddo
      end if
! op_ck_20150128+
      END IF
! op_ck_20150128-
!---------------------------------------
!     liquid - (kpp)gas interaction (for easy scav only)
      kpp_l_idx(js)%gas_spec(:,gas_g2p) = 0
      DO jt=1,lspec_gas
        idx2 = KPP_L_IDX(JS)%GAS_SPEC(JT,GAS_IDX)
        DO J = 1,LSPEC
          IF ( TRIM(STR_FIELD_KPP_L(IDX2))//'_l' == &
               TRIM(STR_FIELD_KPP_L(J)) )         THEN
!!$            print*,"liquid: ",STR_FIELD_KPP_L(IDX2), STR_FIELD_KPP_L(J), &
!!$                    IDX2, J
            kpp_l_idx(js)%gas_spec(jt,gas_g2p) = J
          ENDIF
        END DO
      END DO
!---------------------------------------
!     ice - (kpp)gas interaction (for easy scav only)
      kpp_i_idx(js)%gas_spec(:,gas_g2p) = 0
      DO jt=1,ispec_gas
        idx2 = KPP_I_IDX(JS)%GAS_SPEC(JT,GAS_IDX)
        DO J = 1,ISPEC
          IF ( TRIM(STR_FIELD_KPP_I(IDX2))//'_i' == &
               TRIM(STR_FIELD_KPP_I(J)) )         THEN
!!$            print*,"ice: ",STR_FIELD_KPP_I(IDX2), STR_FIELD_KPP_I(J), &
!!$                    IDX2, J
            kpp_i_idx(js)%gas_spec(jt,gas_g2p) = J
          ENDIF
        END DO
      END DO
!---------------------------------------
!     ice - evaporation - interaction
      CALL info_bi(' Setting up: Ice - Evaporation - Interactions!', substr)

      kpp_i_idx(js)%ice_spec(:,ice2evap) = 0

      DO j=1,ispec_ice
! op_ck_20150129+
        IF (.NOT. ALLOCATED(made3)) THEN
! op_ck_20150129-
        DO k=evap_max_mode,1,-1
          IF (kpp_i_idx(js)%ice_spec(j,ice2evap) /= 0) CYCLE
          DO jt=1,ntrac
            string=TRIM(ti(jt)%tp%ident%basename)//'_i'
            IF ( (TRIM(STR_FIELD_KPP_I(kpp_i_idx(js)%ice_spec(j,ice_idx))) &
              == string) .AND. attr(js)%log_att(jt,laerosol) .AND.     &
              (attr(js)%int_att(jt,tmode) == k) ) THEN
               kpp_i_idx(js)%ice_spec(j,ice2evap) = jt
               EXIT
            END IF
          END DO
        END DO
! op_ck_20150129+
        END IF
! op_ck_20150129-
        IF (kpp_i_idx(js)%ice_spec(j,ice2evap) == 0 ) THEN
          DO jt=1,ntrac
            string=TRIM(ti(jt)%tp%ident%fullname)//'_i'
            IF ( TRIM(STR_FIELD_KPP_I(kpp_i_idx(js)%ice_spec(j,ice_idx))) &
              == TRIM(string) )    THEN
               kpp_i_idx(js)%ice_spec(j,ice2evap) = jt
               EXIT
            ENDIF
          ENDDO
        ENDIF

!!$
!!$
!!$        DO jt=1,ntrac
!!$          string=TRIM(ti(jt)%tp%ident%basename)//'_i'
!!$          DO k=evap_max_mode,1,-1
!!$            IF ( (TRIM(STR_FIELD_KPP_I(kpp_i_idx(js)%ice_spec(j,ice_idx))) &
!!$              == string) .AND. attr(js)%log_att(jt,laerosol) .AND.     &
!!$              (attr(js)%int_att(jt,tmode) == k) ) THEN
!!$              kpp_i_idx(js)%ice_spec(j,ice2evap) = jt
!!$              EXIT
!!$            ENDIF
!!$          ENDDO ! mode
!!$        ENDDO   ! ntrac
!!$        IF (kpp_i_idx(js)%ice_spec(j,ice2evap) == 0 ) THEN
!!$          DO jt=1,ntrac
!!$            string=TRIM(ti(jt)%tp%ident%fullname)//'_i'
!!$            IF ( TRIM(STR_FIELD_KPP_I(kpp_i_idx(js)%ice_spec(j,ice_idx))) &
!!$              == TRIM(string) )    THEN
!!$              kpp_i_idx(js)%ice_spec(j,ice2evap) = jt
!!$              EXIT
!!$            ENDIF
!!$          ENDDO
!!$        ENDIF

! op_ck_20150129+
        IF (.NOT. ALLOCATED(made3)) THEN
! op_ck_20150129-
        ! special cases
        IF (kpp_i_idx(js)%ice_spec(j,ice2evap) /= 0 ) CYCLE
        DO k=evap_max_mode,1,-1
          DO jt=1,ntrac
            IF ( attr(js)%int_att(jt,tmode) == k) THEN
              string=TRIM(ti(jt)%tp%ident%basename)
              IF (kpp_i_idx(js)%ice_spec(j,ice2evap) /= 0 ) CYCLE
              SELECT CASE (string)
              CASE('SO4')
                IF (TRIM(STR_FIELD_KPP_I(kpp_i_idx(js)%ice_spec(j,ice_idx))) &
                  == 'SO4mm_i') kpp_i_idx(js)%ice_spec(j,ice2evap) = jt
                IF (TRIM(STR_FIELD_KPP_I(kpp_i_idx(js)%ice_spec(j,ice_idx))) &
                  == 'HSO4m_i') kpp_i_idx(js)%ice_spec(j,ice2evap) = jt
              CASE('NO3')
                IF (TRIM(STR_FIELD_KPP_I(kpp_i_idx(js)%ice_spec(j,ice_idx))) &
                  == 'NO3m_i') kpp_i_idx(js)%ice_spec(j,ice2evap) = jt
              CASE('NH4')
                IF (TRIM(STR_FIELD_KPP_I(kpp_i_idx(js)%ice_spec(j,ice_idx))) &
                  == 'NH4p_i') kpp_i_idx(js)%ice_spec(j,ice2evap) = jt
              CASE('SO4res')
                IF (TRIM(STR_FIELD_KPP_I(kpp_i_idx(js)%ice_spec(j,ice_idx))) &
                  == 'SO4mm_i') kpp_i_idx(js)%ice_spec(j,ice2evap) = jt
                IF (TRIM(STR_FIELD_KPP_I(kpp_i_idx(js)%ice_spec(j,ice_idx))) &
                  == 'HSO4m_i') kpp_i_idx(js)%ice_spec(j,ice2evap) = jt
              CASE('NO3mres')
                IF (TRIM(STR_FIELD_KPP_I(kpp_i_idx(js)%ice_spec(j,ice_idx))) &
                  == 'NO3m_i') kpp_i_idx(js)%ice_spec(j,ice2evap) = jt
              CASE('NH4pres')
                IF (TRIM(STR_FIELD_KPP_I(kpp_i_idx(js)%ice_spec(j,ice_idx))) &
                  == 'NH4p_i') kpp_i_idx(js)%ice_spec(j,ice2evap) = jt
              CASE('Hpres')
                IF (TRIM(STR_FIELD_KPP_I(kpp_i_idx(js)%ice_spec(j,ice_idx))) &
                  == 'Hp_i') kpp_i_idx(js)%ice_spec(j,ice2evap) = jt
              CASE('Clmres')
                IF (TRIM(STR_FIELD_KPP_I(kpp_i_idx(js)%ice_spec(j,ice_idx))) &
                  == 'Clm_i') kpp_i_idx(js)%ice_spec(j,ice2evap) = jt
              END SELECT
            ENDIF
          ENDDO
        ENDDO
! op_ck_20150129+
        END IF
! op_ck_20150129-
      ENDDO

! op_ck_20150129+
      IF (.NOT. ALLOCATED(made3)) THEN
! op_ck_20150129-
      if (p_parallel_io) then
        do j=1,ispec_ice
          if (kpp_i_idx(js)%ice_spec(j,ice2evap) == 0 ) then
            print*, "This species does not evaporate:",  j,            &
              STR_FIELD_KPP_I(kpp_i_idx(js)%ice_spec(j,ice_idx))
          else
            print*, j, STR_FIELD_KPP_I(kpp_i_idx(js)%ice_spec(j,ice_idx)), &
              kpp_i_idx(js)%ice_spec(j,ice2evap),                        &
              TRIM(ti(kpp_i_idx(js)%ice_spec(j,ice2evap))%tp%ident%fullname)
          end if
        enddo
      end if
! op_ck_20150129+
      END IF
! op_ck_20150129-
!--------------------------------------------------------------

      CALL info_bi(' Setting up: Ice - Liquid - Interactions!', substr)
      kpp_i_idx(js)%ice_spec(:,ice2liq) = 0
      DO j=1,ispec_ice
        IF (ASSOCIATED(strname)) DEALLOCATE (strname)
        NULLIFY(strname)     
        CALL strcrack(STR_FIELD_KPP_I(kpp_i_idx(js)%ice_spec(j,ice_idx)), &
                      '_', strname, dummy)
        DO jt=1,lspec
          string = TRIM(strname(1))//'_l'
          IF ( TRIM(STR_FIELD_KPP_L(jt)) == TRIM(string) )    &
            kpp_i_idx(js)%ice_spec(j,ice2liq) = jt
        ENDDO
      ENDDO
      kpp_i_idx(js)%gas_spec(:,gas_i2l) = 0
      DO j=1,ispec_gas
        DO jt=1,lspec
          IF ( TRIM(STR_FIELD_KPP_L(jt)) == &
            TRIM(STR_FIELD_KPP_I(KPP_I_IDX(JS)%GAS_SPEC(J,GAS_IDX))) ) &
            kpp_i_idx(js)%gas_spec(j,gas_i2l) = jt
        ENDDO
      ENDDO
!---------------------------------------------------
 ! determine the cloud tracer index for each aerosol compound
      aer_attr(js)%aer_spec(:,aer2tracl) = 0
      aer_attr(js)%aer_spec(:,aer2traci) = 0
      DO jt = 1, aer_count(js)%numbers(c_all)
        idx2 = aer_attr(js)%aer_spec(jt,AER_IDX)
        DO j=1,ntrac           
          IF ( (TRIM(ti(idx2)%tp%ident%basename) ==  &
            TRIM(ti(j)%tp%ident%basename)).AND.    &
            (ti(j)%tp%ident%medium==CLOUD ) )     THEN
            IF (TRIM(ti(j)%tp%ident%subname) == 'l') THEN
              aer_attr(js)%aer_spec(jt,aer2tracl) = j
            ELSE IF (TRIM(ti(j)%tp%ident%subname) == 'i') THEN
              aer_attr(js)%aer_spec(jt,aer2traci) = j
            END IF
          END IF
        END DO
      END DO

    ENDDO   ! tracer sets
 
!=============================
    CALL info_bi( "***********************************", substr)
    CALL info_bi( 'Checking for grid box parameters ...', substr)
    
    CALL get_channel_object(status, &
      TRIM(mass(1)), TRIM(mass(2)), p3=grmass)
    IF (status /= 0) &
      CALL error_bi('grmass channel object not found !', substr)

    CALL get_channel_object(status, &
      TRIM(vol(1)), TRIM(vol(2)), p3=grvol)
    IF (status /= 0) &
      CALL error_bi('grvol channel object not found !', substr)

! op_pj_20120307+
!!$    CALL get_channel_object(status, &
!!$      TRIM(tm13d(1)), TRIM(tm13d(2)), p3=tm1_3d)
!!$    IF (status /= 0) &
!!$      CALL error_bi('tm1 channel object not found !', substr)
!!$    
!!$    CALL get_channel_object(status, &
!!$      TRIM(tte3d(1)), TRIM(tte3d(2)), p3=tte_3d)
!!$    IF (status /= 0) &
!!$      CALL error_bi('tte channel object not found !', substr)
! op_pj_20120307-
    
    CALL get_channel_object(status, &
      TRIM(press(1)), TRIM(press(2)), p3=press_3d)
    IF (status /= 0) &
      CALL error_bi('pressure channel object not found !', substr)
    
    CALL get_channel_object(status, &
      TRIM(pressi(1)), TRIM(pressi(2)), p3=pressi_3d)
    IF (status /= 0) &
      CALL error_bi('pressure interface channel object not found !', substr)

    CALL get_channel_object(status, &
      TRIM(boxarea(1)), TRIM(boxarea(2)), p2=gboxarea_2d) ! um_ak_20110601 namelist parameter
    IF (status /= 0) &
      CALL error_bi('gridbox area channel object not found !', substr)

    

    CALL info_bi('... getting photolysis rates in liquid phase !', substr)
    
    CALL get_channel_object(status, TRIM(photol(1)), 'J_H2O2', p3=jval_h2o2)
    IF (status /= 0) THEN
      CALL info_bi('J_H2O2 channel object not found !', substr)
      CALL info_bi('no aqueous phase photolysis of H2O2', substr)
      CALL info_bi('There will be no photolysis in aqueous phase!', substr)
    ELSE
      CALL info_bi( &
        'WARNING !!! Rerun activated for photolysis rate of H2O2', substr)
      CALL set_channel_object_restreq(status, TRIM(photol(1)), 'J_H2O2' )
      CALL channel_halt(substr, status)
    END IF

    CALL get_channel_object(status, TRIM(photol(1)), 'J_HNO3', p3=jval_hno3)
    IF (status /= 0) THEN
      CALL info_bi('J_HNO3 channel object not found !', substr)
      CALL info_bi('no aqueous phase photolysis of HNO3', substr)
      CALL info_bi('There will be no photolysis in aqueous phase!', substr)
    ELSE
      CALL info_bi( &
        'WARNING !!! Rerun activated for photolysis rate of HNO3', substr)
      CALL set_channel_object_restreq(status, TRIM(photol(1)), 'J_HNO3' )
      CALL channel_halt(substr, status)
    END IF

    CALL get_channel_object(status, TRIM(photol(1)), 'J_O3P', p3=jval_o3)
    IF (status /= 0) THEN
      CALL info_bi('J_O3P channel object not found !', substr)
      CALL info_bi('no aqueous phase photolysis of O3', substr)
      CALL info_bi('There will be no photolysis in aqueous phase!', substr)
    ELSE
      CALL info_bi( &
        'WARNING !!! Rerun activated for photolysis rate of O3', substr)
      CALL set_channel_object_restreq(status, TRIM(photol(1)), 'J_O3P' )
      CALL channel_halt(substr, status)
    END IF
    
!Defining top level of scavenging for performance improvements
       
!!$       altitude = 5000. ! op_pj_20120310
       IF (p_parallel_io) &
         PRINT*, substr, ': Restrict scavenging to a height below ', &
         altitude/100, ' hPa!'
#ifndef COSMO
       DO jk=1,nvclev
         h_a(jk) = vct(jk)
         h_b(jk) = vct(jk+nvclev)
       ENDDO
#endif
       sfpress     = 1.e5_dp   ! reference pressure of 1000 hPa
       DO jk=1,nlev+1
         hypi(jk)      = h_a(jk) + h_b(jk) * sfpress
       ENDDO
!!$#else
!!$       DO jk=1,nlev+1
!!$         hypi(jk)      = sigmr(jk) * p0sl
!!$       ENDDO
!!$#endif
       max_lev_scav = 1
       DO jk=1,nlev
         IF (hypi(jk) < altitude .AND. hypi(jk+1) >= altitude) THEN
           max_lev_scav = jk
           EXIT
         ENDIF
       END DO
       max_lev_scav = MAX(2, max_lev_scav) ! op_pj_20120310 (jk-1 !!!)
       IF (p_parallel_io) &
         PRINT*,                                                         &
         substr, ': Restrict scavenging to a level number higher than ', &
         max_lev_scav, ' !'

!!#D attila +
!      Coupling to ATTILA for pseudo - lagrangian scavenging
       IF (lscav_lg) THEN
         CALL get_channel_object(status, 'attila','NCB', p3=NCB)
         IF (status /= 0) &
           CALL error_bi(&
           'no ATTILA channel object for number of cells per gid box found!'&
           , substr)
       ENDIF
!!#D attila -

       CALL end_message_bi(modstr,'COUPLING INITIALIZATION',substr)

     END SUBROUTINE scav_init_coupling


! ==============================================================================
!===============================================================================
  SUBROUTINE scav_convec

#ifdef MESSYTENDENCY
    USE messy_main_tendency_bi,    ONLY: mtend_get_start_l, mtend_id_tracer &
                                         , mtend_add_l
#else
    USE messy_main_tracer_mem_bi,  ONLY: pxtte => qxtte,   pxtm1 => qxtm1
#endif
    USE messy_main_grid_def_mem_bi, ONLY: jrow, nlev, kproma
    USE messy_main_timer,          ONLY: time_step_len
    USE messy_scav_liq,            ONLY: jval_h2o2, jval_h2o2_2d, &
                                         jval_o3,   jval_o3_2d,   &
                                         jval_hno3, jval_hno3_2d

    IMPLICIT NONE

    INTEGER  :: js
    REAL(dp), POINTER, DIMENSION(:,:,:) :: pxtp1 => NULL()

    INTEGER  :: jl, jk, jt, i, jm, mode

    IF (.NOT. lscav_gp) RETURN
    IF (.NOT. cvtrans_use) RETURN

    l_lg  = .FALSE.
    ntrac = ntrac_gp
    ti    => ti_gp
    js = 1
    ALLOCATE(pxtp1(nproma,nlev,ntrac))
    pxtp1(:,:,:)     = 0._dp
    ALLOCATE(xtte_scav(nproma,nlev,ntrac))
    xtte_scav(:,:,:) = 0._dp
! um_ak_20110615+ 
!!$ tracer_conv  => trac_field(:,:,:,jrow)
!!$ jval_h2o2_2d => jval_h2o2(:,:,jrow)
    tracer_conv  => trac_field(_RI_XYZN_(:,jrow,:,:))
    jval_h2o2_2d => jval_h2o2(_RI_XYZ__(:,jrow,:))
    jval_hno3_2d => jval_hno3(_RI_XYZ__(:,jrow,:))
    jval_o3_2d   => jval_o3(_RI_XYZ__(:,jrow,:))
! um_ak_20110615-

    ! op_pg_20130225+
    IF (L_TEND) THEN
       DO jt=1, n_te
          IF (te_idt(jt) /= 0) THEN
             tte_scav(jt)%ptr(_RI_XYZ__(:,jrow,:)) = 0._dp
          END IF
       END DO
    END IF
    ! op_pg_20130225-


    IF (column) THEN
      DO jt=1,ntrac
        DO jk=max_lev_scav, nlev
          DO jl=1,kproma
            pxtp1(jl,jk,jt) = tracer_conv(_RI_X_ZN_(jl,jk,jt))
          ENDDO
        ENDDO
      ENDDO
    ELSE
#ifdef MESSYTENDENCY
       ! ub_ak_20190612+
!      CALL mtend_get_start_l(mtend_id_tracer, v0t=pxtp1)
       DO jt=1,ntrac
!          CALL mtend_get_start_l(jt, v0=pxtp1(_RI_X_ZN_(:,1:nlev,jt)))
          CALL mtend_get_start_l(jt, v0=pxtp1(:,:,jt))
       ENDDO
       ! ub_ak_20190612-
#else
      DO jt=1,ntrac
        DO jk=max_lev_scav, nlev
          DO jl=1,kproma
            pxtp1(jl,jk,jt) = pxtm1(_RI_X_ZN_(jl,jk,jt)) + &
                 pxtte(_RI_X_ZN_(jl,jk,jt))*time_step_len
          ENDDO
        ENDDO
      ENDDO
#endif
    ENDIF

!!$    do jt=1,ntrac
!!$      do jk=1,nlev
!!$        do jl=1,kproma
!!$          if (pxtp1(jl,jk,jt) < -1.e-20_dp) &
!!$          print*, "in scav_convec start", js, jl, jk, jt, jrow, &
!!$            pxtp1(jl,jk,jt), pxtm1(jl,jk,jt), pxtte(jl,jk,jt),  &
!!$            tracer_conv(jl,jk,jt)
!!$        enddo
!!$      enddo
!!$    enddo

!   choose correct aerosol wetradius
    DO i = 1, aermodel(js)%aeromodnum
      mode = aermodel(js)%aer_input(i)%lmode
      IF (ASSOCIATED(aermodel(js)%aer_input(i)%inputrad)) THEN
        ALLOCATE(aermodel(js)%aer_input(i)%wetrad(nproma,mode,nlev))
        DO jm=1,mode
          DO jl=1,kproma
            DO jk=1,nlev
              aermodel(js)%aer_input(i)%wetrad(jl,jm,jk) = &
                aermodel(js)%aer_input(i)%inputrad(_RI_XYZN_(jl,jrow,jk,jm))
            END DO
          ENDDO
        ENDDO
      END IF
! op_ck_20130116+
      IF (ASSOCIATED(aermodel(js)%aer_input(i)%inputpf)) THEN
        ALLOCATE(aermodel(js)%aer_input(i)%philfrac(nproma,mode,nlev))
        DO jm=1,mode
          DO jl=1,kproma
            DO jk=1,nlev
              aermodel(js)%aer_input(i)%philfrac(jl,jm,jk) = &
                aermodel(js)%aer_input(i)%inputpf(_RI_XYZN_(jl,jrow,jk,jm))
            END DO
          ENDDO
        ENDDO
      END IF
! op_ck_20130116-
    ENDDO

!      CALL scav_cv(js,jrow,kproma,pxtp1)
    CALL scav_cv_new(js,jrow,kproma,pxtp1)

    DO i = 1, aermodel(js)%aeromodnum
      IF (ASSOCIATED(aermodel(js)%aer_input(i)%inputrad)) THEN
        DEALLOCATE(aermodel(js)%aer_input(i)%wetrad)
        NULLIFY(aermodel(js)%aer_input(i)%wetrad)
      ENDIF
! op_ck_20130116+
      IF (ASSOCIATED(aermodel(js)%aer_input(i)%inputpf)) THEN
        DEALLOCATE(aermodel(js)%aer_input(i)%philfrac)
        NULLIFY(aermodel(js)%aer_input(i)%philfrac)
      ENDIF
! op_ck_20130116-
    ENDDO

! mz_sg_20120425+
#ifdef MECCA_TAG
  ! this code is executed only when tagging is used
  ! calling mecca_tag to process scav. for related tagging tracers
    CALL mecca_tag_calc_xtte4scav(xtte_scav, max_lev_scav, pxtp1, kproma)
#endif
! mz_sg_20120425-

    IF (column) THEN
      DO jt=1,ntrac
        DO jk=max_lev_scav, nlev
          DO jl=1,kproma
            tracer_conv(_RI_X_ZN_(jl,jk,jt)) = pxtp1(jl,jk,jt) + &
              xtte_scav(jl,jk,jt) * time_step_len
          ENDDO
        ENDDO
      ENDDO
    ELSE
       ! um_ak_20100419+
       !pxtte(1:kproma,1:nlev,1:ntrac) = pxtte(1:kproma,1:nlev,1:ntrac) +  &
       !                                xtte_scav(1:kproma,1:nlev,1:ntrac) 
#ifdef MESSYTENDENCY
       ! ub_ak_20190612+
           
!      CALL mtend_add_l(tend_cv, mtend_id_tracer, &
!            pxt=xtte_scav(1:kproma,1:nlev,1:ntrac))
       DO jt = 1, ntrac
          CALL mtend_add_l(tend_cv, jt, px=xtte_scav(1:kproma,1:nlev,jt))
       END DO
       ! ub_ak_20190612-
#else
       ! ub_ak_20190612-
       DO jt=1,ntrac
          DO jk=1, nlev
             pxtte(_RI_X_ZN_(1:kproma,jk,jt)) =  pxtte(_RI_X_ZN_(1:kproma,jk,jt)) +  &
                  xtte_scav(1:kproma,jk,jt) 
        ENDDO
      ENDDO
#endif
      ! um_ak_20100419-
    ENDIF

! op_pg_20130225+
    IF (L_TEND) THEN
       DO jt=1, n_te
          IF (te_idt(jt) /= 0) THEN

             DO jk=1, nlev
                DO jl=1,kproma
                   tte_scav(jt)%ptr(_RI_XYZ__(jl,jrow,jk)) = &
                        tte_scav(jt)%ptr(_RI_XYZ__(jl,jrow,jk)) + &
                        xtte_scav(jl,jk,te_idt(jt))
                ENDDO
             ENDDO

          END IF
       END DO
    END IF
! op_pg_20130225-

!!$    do jt=1,ntrac
!!$      do jk=1,nlev
!!$        do jl =1,kproma
!!$          if (pxtm1(jl,jk,jt) + pxtte(jl,jk,jt) * time_step_len < -1.e-20_dp) &
!!$          print*, "in scav_convec finish", js, jl, jk, jt, jrow, &
!!$            pxtm1(jl,jk,jt) + pxtte(jl,jk,jt) * time_step_len,     &
!!$            xtte_scav(jl,jk,jt), pxtp1(jl,jk,jt), pxtte(jl,jk,jt)
!!$        end do
!!$      enddo
!!$    enddo
    DEALLOCATE(xtte_scav)
    NULLIFY(xtte_scav)
    DEALLOCATE(pxtp1)
    NULLIFY(pxtp1)

  END SUBROUTINE scav_convec

!===============================================================================

  SUBROUTINE scav_physc(call_type)

#ifdef MESSYTENDENCY
    USE messy_main_tendency_bi,    ONLY: mtend_get_start_l, mtend_id_tracer &
                                         , mtend_add_l
#else
    USE messy_main_tracer_mem_bi,  ONLY: pxtte => qxtte,   pxtm1 => qxtm1
#endif
    USE messy_main_grid_def_mem_bi, ONLY: jrow, kproma, nlev
    USE messy_scav_l_kpp,          ONLY: lspec => nspec
    USE messy_scav_i_kpp,          ONLY: ispec => nspec
    USE messy_scav_liq,            ONLY: jval_h2o2, jval_h2o2_2d, &
                                         jval_o3,   jval_o3_2d,   &
                                         jval_hno3, jval_hno3_2d

    USE messy_scav_inter,          ONLY: process
    USE messy_main_timer,          ONLY: time_step_len  

    IMPLICIT NONE

    INTEGER, INTENT(IN) :: call_type

    INTEGER  :: js
    REAL(dp), POINTER, DIMENSION(:,:,:) :: pxtp1           => NULL()
    REAL(dp), POINTER, DIMENSION(:,:,:) :: cloud_field_aer => NULL()
    REAL(dp), POINTER, DIMENSION(:,:,:) :: cloud_field_kpp => NULL()
    INTEGER  :: jl, jk, jt, jm, nspec, phase, i, mode
    INTEGER  :: idt, jjrow, j ! op_pj_20160926

    IF (.NOT. lscav_gp) RETURN

    l_lg  = .FALSE.
    ntrac = ntrac_gp
    ti    => ti_gp
    js = 1
    ALLOCATE(pxtp1(nproma,nlev,ntrac))
    pxtp1(:,:,:) = 0.0_dp
    ALLOCATE(xtte_scav(nproma,nlev,ntrac))
    xtte_scav(:,:,:) = 0._dp
    
#ifdef MESSYTENDENCY
    ! ub_ak_20190612+
    !      CALL mtend_get_start_l(mtend_id_tracer, v0t=pxtp1)
    DO jt=1,ntrac
!       CALL mtend_get_start_l(jt, v0=pxtp1(_RI_X_ZN_(:,1:nlev,jt)))
       CALL mtend_get_start_l(jt, v0=pxtp1(:,:,jt))
    ENDDO
    ! ub_ak_20190612-
#else
  ! um_ak_20100419+
    DO jk=1,nlev
       DO jt = 1, ntrac
          !pxtp1(:,:,:) = pxtm1(:,:,:) + pxtte(:,:,:) * time_step_len
          pxtp1(1:kproma,jk,jt) = pxtm1(_RI_X_ZN_(1:kproma,jk,jt))  &
               + pxtte(_RI_X_ZN_(1:kproma,jk,jt)) * time_step_len
       ENDDO
    ENDDO
  ! um_ak_20100419-
#endif

    SELECT CASE(call_type)

    CASE(1)
! calculate ls_scavenging

       ! um_ak_20100419-
!!$    do jt=1,ntrac
!!$      do jk=1,nlev
!!$        do jl=1,kproma
!!$          if (pxtp1(jl,jk,jt) < -1.e-20_dp) &
!!$          print*, "in scav_physc start", js, jl, jk, jt, jrow, &
!!$            pxtp1(jl,jk,jt), pxtm1(jl,jk,jt), pxtte(jl,jk,jt)
!!$        enddo
!!$      enddo
!!$    enddo

      jval_h2o2_2d => jval_h2o2(_RI_XYZ__(:,jrow,:))
      jval_hno3_2d => jval_hno3(_RI_XYZ__(:,jrow,:))
      jval_o3_2d   => jval_o3(_RI_XYZ__(:,jrow,:))

!   choose correct aerosol wetradius
      DO i = 1, aermodel(js)%aeromodnum
        mode = aermodel(js)%aer_input(i)%lmode
        IF (ASSOCIATED(aermodel(js)%aer_input(i)%inputrad)) THEN
          ALLOCATE(aermodel(js)%aer_input(i)%wetrad(nproma,mode,nlev))
          DO jm=1,mode
            DO jl=1,kproma
              DO jk=1,nlev
                aermodel(js)%aer_input(i)%wetrad(jl,jm,jk) = &
                  aermodel(js)%aer_input(i)%inputrad(_RI_XYZN_(jl,jrow,jk,jm))
              END DO
            ENDDO
          ENDDO
        END IF
! op_ck_20130116+
        IF (ASSOCIATED(aermodel(js)%aer_input(i)%inputpf)) THEN
          ALLOCATE(aermodel(js)%aer_input(i)%philfrac(nproma,mode,nlev))
          DO jm=1,mode
            DO jl=1,kproma
              DO jk=1,nlev
                aermodel(js)%aer_input(i)%philfrac(jl,jm,jk) = &
                  aermodel(js)%aer_input(i)%inputpf(_RI_XYZN_(jl,jrow,jk,jm))
              END DO
            ENDDO
          ENDDO
        END IF
! op_ck_20130116-
        IF (ASSOCIATED(aermodel(js)%aer_input(i)%mfracnuc)) THEN
          mode = SIZE(aermodel(js)%aer_input(i)%mfracnuc,3)
          ALLOCATE(aermodel(js)%aer_input(i)%mfrac(nproma,mode,nlev))
          DO jm=1,mode
            DO jl=1,kproma
              DO jk=1,nlev
                aermodel(js)%aer_input(i)%mfrac(jl,jm,jk) = &
                  aermodel(js)%aer_input(i)%mfracnuc(_RI_XYZN_(jl,jrow,jk,jm))
 !               IF (.NOT. ((aermodel(js)%aer_input(i)%mfrac(jl,jm,jk) > 0._dp) .AND. &
 !                   (aermodel(js)%aer_input(i)%mfrac(jl,jm,jk) < 1._dp) ))&
 !                   print*, "problem mfrac:", i, js,jl, jm, jk, &
 !                   aermodel(js)%aer_input(i)%mfrac(jl,jm,jk), &
 !                   aermodel(js)%aer_input(i)%mfracnuc(_RI_XYZN_(jl,jrow,jk,jm))
              END DO
            ENDDO
          ENDDO
        END IF
        IF (ASSOCIATED(aermodel(js)%aer_input(i)%nfracnuc)) THEN
          mode = SIZE(aermodel(js)%aer_input(i)%nfracnuc,3)
          ALLOCATE(aermodel(js)%aer_input(i)%nfrac(nproma,mode,nlev))
          DO jm=1,mode
            DO jl=1,kproma
              DO jk=1,nlev
                aermodel(js)%aer_input(i)%nfrac(jl,jm,jk) = &
                  aermodel(js)%aer_input(i)%nfracnuc(_RI_XYZN_(jl,jrow,jk,jm))
 !               IF (.NOT. ((aermodel(js)%aer_input(i)%nfrac(jl,jm,jk) > 0._dp) .AND. &
 !                          (aermodel(js)%aer_input(i)%nfrac(jl,jm,jk) < 1._dp) ) )&
 !                   print*, "problem nfrac:", i, js,jl, jm, jk, &
 !                   aermodel(js)%aer_input(i)%nfrac(jl,jm,jk),  &
 !                   aermodel(js)%aer_input(i)%nfracnuc(_RI_XYZN_(jl,jrow,jk,jm))

              END DO
            ENDDO
          ENDDO
        END IF
        
      ENDDO

!      CALL scav_ls(js,jrow,kproma,pxtp1)
      CALL scav_ls_new(js,jrow,kproma,pxtp1) ! changes xtte_scav

      DO i = 1, aermodel(js)%aeromodnum
        IF (ASSOCIATED(aermodel(js)%aer_input(i)%inputrad)) THEN
          DEALLOCATE(aermodel(js)%aer_input(i)%wetrad)
          NULLIFY(aermodel(js)%aer_input(i)%wetrad)
        ENDIF
        IF (ASSOCIATED(aermodel(js)%aer_input(i)%mfrac)) THEN
          DEALLOCATE(aermodel(js)%aer_input(i)%mfrac)
          NULLIFY(aermodel(js)%aer_input(i)%mfrac)
        ENDIF
        IF (ASSOCIATED(aermodel(js)%aer_input(i)%nfrac)) THEN
          DEALLOCATE(aermodel(js)%aer_input(i)%nfrac)
          NULLIFY(aermodel(js)%aer_input(i)%nfrac)
        ENDIF
! op_ck_20130116+
        IF (ASSOCIATED(aermodel(js)%aer_input(i)%inputpf)) THEN
          DEALLOCATE(aermodel(js)%aer_input(i)%philfrac)
          NULLIFY(aermodel(js)%aer_input(i)%philfrac)
        ENDIF
! op_ck_20130116-
      ENDDO

#ifdef MESSYTENDENCY
      tend_now = tend_ls
#endif


    CASE(2)
!!$! mz_sg_20120425+
!!$#ifndef MECCA_TAG
!!$  ! this code is NOT executed when tagging is used
!!$  ! pxtp1() value is not in use below throughout the routine but zeroed; 
!!$  ! it should be preserved for mecca_tag_calc_xtte4scav()
!!$
!!$      pxtp1(:,:,:) = 0._dp
!!$
!!$#endif
!!$! mz_sg_20120425-

! calculate re-evaporation of clouds and release of dissolved / adsorbed species
! from the cloud phases
      DO jm=1,4
        SELECT CASE(jm)
        CASE(1)
          IF (lscav_cv .AND. lscav_i) THEN
            cloud_field_kpp => cv_field(js,2)%cloud_field(:,:,:,jrow)
            cloud_field_aer => aero_flx(js)%cloud_field_cv_aer(_RI_XYZN_(:,jrow,:,2),:)
            phase = 2
            nspec = ispec
            jjrow = jrow
            DO idt = 1, ntrac
               DO j = 1, n_resmod
                  process(3)%frac_eva(idt,j)%ptr => &
                      set(js)%para(2)%proc(4)%evafrac(_RI_XYZN_(:,jjrow,:,idt),j) 
               END DO
            END DO
          ELSE 
            CYCLE
          ENDIF
        CASE(2)
          IF (lscav_cv .AND. lscav_l) THEN
            cloud_field_kpp => cv_field(js,1)%cloud_field(:,:,:,jrow)
            cloud_field_aer => aero_flx(js)%cloud_field_cv_aer(_RI_XYZN_(:,jrow,:,1),:)
            phase = 1
            nspec = lspec
            jjrow = jrow
            DO idt = 1, ntrac
               DO j = 1, n_resmod
                  process(3)%frac_eva(idt,j)%ptr => &
                      set(js)%para(2)%proc(3)%evafrac(_RI_XYZN_(:,jjrow,:,idt),j)
               END DO
            END DO
          ELSE 
            CYCLE
          ENDIF
        CASE(3)
          IF (lscav_ls .AND. lscav_i) THEN
            cloud_field_kpp => ls_field(js,2)%cloud_field(:,:,:,jrow)
            cloud_field_aer => aero_flx(js)%cloud_field_ls_aer(_RI_XYZN_(:,jrow,:,2),:)
            phase = 2
            nspec = ispec
            jjrow = jrow
            DO idt = 1, ntrac
               DO j = 1, n_resmod
                  process(3)%frac_eva(idt,j)%ptr => &
                      set(js)%para(1)%proc(4)%evafrac(_RI_XYZN_(:,jjrow,:,idt),j)
               END DO
            END DO
          ELSE 
            CYCLE
          ENDIF
        CASE(4)
          IF (lscav_ls .AND. lscav_l) THEN
            cloud_field_kpp => ls_field(js,1)%cloud_field(:,:,:,jrow)
            cloud_field_aer => aero_flx(js)%cloud_field_ls_aer(_RI_XYZN_(:,jrow,:,1),:)
            phase = 1
            nspec = lspec
            jjrow = jrow
            DO idt = 1, ntrac
               DO j = 1, n_resmod
                  process(3)%frac_eva(idt,j)%ptr => &
                      set(js)%para(1)%proc(3)%evafrac(_RI_XYZN_(:,jjrow,:,idt),j)
               END DO
            END DO
          ELSE 
            CYCLE
          ENDIF
        END SELECT
!        print*, "dimensions aer",jrow, js, phase, lbound(cloud_field_aer), &
!          ubound(cloud_field_aer), &
!          "kpp",lbound(cloud_field_kpp), ubound(cloud_field_kpp) 

        CALL evapo_cloud(js, jrow, kproma, phase, nspec, ntrac, &
                         cloud_field_aer(1:kproma,:,:),     &
                         cloud_field_kpp(1:kproma,:,:) )

        cloud_field_aer(:,:,:) = 0._dp
        cloud_field_kpp(:,:,:) = 0._dp
      ENDDO

#ifdef MESSYTENDENCY
      tend_now = tend_ev
#endif

    END SELECT

! mz_sg_20120425+
#ifdef MECCA_TAG
  ! this code is executed only when tagging is used
  ! calling mecca_tag to process scav. for related tagging tracers
    CALL mecca_tag_calc_xtte4scav(xtte_scav, max_lev_scav, pxtp1, kproma)
#endif
! mz_sg_20120425-

#ifndef MESSYTENDENCY
        ! um_ak_20100419+
       !pxtte(1:kproma,1:nlev,1:ntrac) = pxtte(1:kproma,1:nlev,1:ntrac) +  &
       !                                xtte_scav(1:kproma,1:nlev,1:ntrac) 
       DO jt=1,ntrac
          DO jk=1, nlev
             pxtte(_RI_X_ZN_(1:kproma,jk,jt)) =  pxtte(_RI_X_ZN_(1:kproma,jk,jt)) +  &
                  xtte_scav(1:kproma,jk,jt) 
        ENDDO
      ENDDO
      ! um_ak_20100419-
#else
      ! ub_ak_20190612+
      !CALL mtend_add_l(tend_now, mtend_id_tracer, &
      !    pxt=xtte_scav(1:kproma,1:nlev,1:ntrac))
      DO jt = 1, ntrac
         CALL mtend_add_l(tend_now, jt, px=xtte_scav(1:kproma,1:nlev,jt))
      END DO
      ! ub_ak_20190612-
#endif

! op_pg_20130225+
    IF (L_TEND) THEN
       DO jt=1, n_te
          IF (te_idt(jt) /= 0) THEN

             DO jk=1, nlev
                DO jl=1,kproma
                   tte_scav(jt)%ptr(_RI_XYZ__(jl,jrow,jk)) = &
                        tte_scav(jt)%ptr(_RI_XYZ__(jl,jrow,jk)) + &
                        xtte_scav(jl,jk,te_idt(jt))
                ENDDO
             ENDDO

          END IF
       END DO
    END IF
! op_pg_20130225-

!!$    do jt=1,ntrac
!!$      do jk=1,nlev
!!$        do jl =1,kproma
!!$          if (pxtm1(jl,jk,jt) + pxtte(jl,jk,jt) * time_step_len < -1.e-15_dp) &
!!$          print*, "in scav_physc finish", js, jl, jk, jt, jrow, &
!!$            pxtm1(jl,jk,jt) + pxtte(jl,jk,jt) * time_step_len,     &
!!$            xtte_scav(jl,jk,jt), pxtp1(jl,jk,jt)
!!$        end do
!!$      enddo
!!$    enddo
    DEALLOCATE(xtte_scav)
    NULLIFY(xtte_scav)
    DEALLOCATE(pxtp1)
    NULLIFY(pxtp1)

  END SUBROUTINE scav_physc

!===============================================================================

!!$  SUBROUTINE scav_global_end(call_type)
!!$
!!$    IMPLICIT NONE
!!$
!!$    INTEGER, INTENT(IN) :: call_type
!!$
!!$  END SUBROUTINE scav_global_end

!==============================================================================
! This subroutine collects all the specific parameters for 
! CONVECTIVE scavenging, calls the scav_main driver, and 
! determines the tracer tendencies from convective scavenging

  SUBROUTINE scav_cv_new(js, jrow, kproma, pxtp1)

    ! ECHAM5/MESSy
    USE messy_main_data_bi,         ONLY: tte_3d, tm1
    USE messy_main_grid_def_mem_bi, ONLY: nlev
    ! MESSy
    USE messy_main_timer,           ONLY: time_step_len
    USE messy_main_constants_mem,   ONLY: avo => N_A, rd, R_gas, tmelt
    USE messy_scav_inter,           ONLY: scav_main, trac2l, trac2i, &
                                          process
    USE messy_scav_l_kpp,           ONLY: lspec => nspec
    USE messy_scav_i_kpp,           ONLY: ispec => nspec

    IMPLICIT NONE

    INTEGER,  INTENT(IN) :: js, jrow, kproma
    REAL(dp), INTENT(IN) :: pxtp1(nproma,nlev,ntrac)
    INTEGER :: idt, j, jjrow, p ! op_pj_20160926

! op_ck_20140514+
    CHARACTER(LEN=*), PARAMETER :: substr = 'scav_cv_new'
    INTEGER  :: status = 0
! op_ck_20140514-
    INTEGER  :: jl, jk, jt, idx1, idx2
    INTEGER  :: mlev ! um_ak_20100802 (RI)
    REAL(dp) :: zt(kproma,nlev), zrhoa(kproma, nlev)
    REAL(dp) :: mc(kproma,nlev), cm(kproma, nlev)
    REAL(dp) :: zxtp1(kproma, nlev, 0:ntrac)
    REAL(dp) :: zxt_vol(kproma, ntrac)
    REAL(dp) :: zcover(kproma,nlev), zcov(kproma)
    REAL(dp) :: zrainprod(kproma,nlev), zsnowprod(kproma,nlev)
    REAL(dp) :: zrlwc(kproma,nlev), zriwc(kproma,nlev)
    REAL(dp) :: rain(kproma,nlev),  snow(kproma,nlev)
    REAL(dp) :: lfrac(kproma,nlev), melt(kproma,nlev)
    REAL(dp) :: xrprod(kproma), xsprod(kproma)
    REAL(dp) :: xlwc(kproma), xiwc(kproma)
    REAL(dp) :: dummy(kproma,nlev)
    REAL(dp) :: prod_array(kproma,nlev,nprod,2)
    LOGICAL  :: lnucl_i(kproma,nlev), lnucl_l(kproma,nlev)
    LOGICAL  :: pseudo_ls = .FALSE.

    IF (.NOT.lscav_cv) RETURN
    IF (.NOT.cvtrans_use) THEN
      IF ( (lscav_gp) .AND. (js == 1) ) RETURN
    ENDIF
    pseudo_ls = .FALSE.
    IF (.NOT. lscav_ls) pseudo_ls = .TRUE.

    DO p=1, 4
       IF ( ASSOCIATED(set(js)%para(2)%proc(p)%evafrac) ) THEN
          jjrow = jrow
          DO idt = 1, ntrac
             DO j = 1, n_resmod
                process(p)%frac_eva(idt,j)%ptr  => &
                     set(js)%para(2)%proc(p)%evafrac(_RI_XYZN_(:,jjrow,:,idt),j)
                process(p)%frac_eva(idt,j)%ptr(:,:) = 0.0_dp
             END DO
          END DO
       END IF
    END DO
    ! op_pj_20160926-


! initialising
    dummy     = 0._dp
    melt      = 0._dp
    zrlwc     = 0._dp
    zriwc     = 0._dp
    zrainprod = 0._dp
    zsnowprod = 0._dp
    ! op_pj_20101209+
    zt        = 0._dp
    zrhoa     = 0._dp
    zxtp1     = 0._dp
    ! op_pj_20101209-
    zcover    = 0._dp     ! op_pj_20120310

    lwc_cv(_RI_XYZ__(:,jrow,:)) = 0._dp

    DO jk=max_lev_scav,nlev
      DO jl=1,kproma
        zt(jl,jk) = tm1(_RI_XYZ__(jl,jrow,jk)) &
             + tte_3d(_RI_XYZ__(jl,jrow,jk)) * time_step_len
        zrhoa(jl,jk) = press_3d(_RI_XYZ__(jl,jrow,jk)) / (rd * zt(jl,jk)) 
        ! convert mixing ratio [mol/mol] to concentration [mcl/cc]
        ! mc = mixing ratio to concentration         
        mc(jl,jk) = (avo/1.e6_dp) * press_3d(_RI_XYZ__(jl,jrow,jk)) / (R_gas * zt(jl,jk))
        ! convert concentration [mcl/cc] to mixing ratio [mol/mol]
        ! cm = concentration to mixing ratio
        cm(jl,jk)= 1._dp/mc(jl,jk)
      
        ! CALCULATE LIQUID FRACTION OF SNOWFLUX (SUPER-COOLED WATER)        
        lfrac(jl,jk) = 0.0_dp
        IF ((tmelt-LIQFRAC(2)) > EPSILON(tmelt)) THEN
          lfrac(jl,jk) = MAX( 0.0_dp, LIQFRAC(1)/(tmelt-LIQFRAC(2)) * &
                              (zt(jl,jk) - LIQFRAC(2)) )
        END IF

        !      rain(jl,jk) = precflx_cv(jl,jk,jrow)
        rain(jl,jk) = precflx_cv(_RI_XYZ__(jl,jrow,jk)) + &
                      lfrac(jl,jk)*snowflx_cv(_RI_XYZ__(jl,jrow,jk))
        !      snow(jl,jk) = snowflx_cv(jl,jk,jrow)
        snow(jl,jk) = (1._dp-lfrac(jl,jk))*snowflx_cv(_RI_XYZ__(jl,jrow,jk))
        IF (lfrac(jl,jk) > 0._dp) melt(jl,jk) = rain(jl,jk) - rain(jl,jk-1)
      ENDDO
    ENDDO

    ph(js)%rainph_cv(_RI_XYZ__(:,jrow,:))   = 19.99_dp
    ph(js)%cloudph_cv(_RI_XYZ__(:,jrow,:))  = 19.99_dp
    ph(js)%Hp_cloud_cv(_RI_XYZ__(:,jrow,:)) = 0._dp

! um_ak_20110615+
    KPP_INFO(js)%cv_rain_steps(_RI_XYZ__(:,jrow,:))  = 0._dp
    KPP_INFO(js)%cv_cloud_steps(_RI_XYZ__(:,jrow,:)) = 0._dp
    KPP_INFO(js)%cv_rain_rsteps(_RI_XYZ__(:,jrow,:))  = 0._dp
    KPP_INFO(js)%cv_cloud_rsteps(_RI_XYZ__(:,jrow,:)) = 0._dp
!!$    KPP_INFO(js)%cv_rain_steps(:,:,jrow)  = 0._dp
!!$    KPP_INFO(js)%cv_cloud_steps(:,:,jrow) = 0._dp
!!$    KPP_INFO(js)%cv_rain_rsteps(:,:,jrow)  = 0._dp
!!$    KPP_INFO(js)%cv_cloud_rsteps(:,:,jrow) = 0._dp
! um_ak_20110615-

    DO jt=1,ntrac
      DO jk=max_lev_scav, nlev
        DO jl=1,kproma
          zxtp1(jl,jk,jt) = pxtp1(jl,jk,jt) * mc(jl,jk)
        ENDDO
      ENDDO
    ENDDO

    aero_flx(js)%aer_flx_cv(_RI_XYZN_(:,jrow,:,:))  = 0.0_dp
    kpp_rain_cv(_RI_XYZN_(:,jrow,:,:)) = 0.0_dp
    kpp_snow_cv(_RI_XYZN_(:,jrow,:,:)) = 0.0_dp

    DO jk=max_lev_scav,nlev

       ! op_pj_20120310+
       IF (jk < nlev) THEN
          mlev = jk+1 ! um_ak_20100802 (RI)
       ELSE
          mlev = jk   ! um_ak_20100802 (RI)
       END IF
       ! op_pj_20120310-

      DO jl=1,kproma
        zcover(jl,jk) = MAXVAL(conv_cover(_RI_XYZ__(jl,jrow,1:mlev)))

        IF (column) zcover(jl,jk) = 1._dp

        IF (zcover(jl,jk) > 1.e-10_dp) THEN
          xlwc(jl)      = pcvlwc(_RI_XYZ__(jl,jrow,jk)) + &
                          lfrac(jl,jk) * pcviwc(_RI_XYZ__(jl,jrow,jk))
          xiwc(jl)      = (1._dp - lfrac(jl,jk)) * pcviwc(_RI_XYZ__(jl,jrow,jk))
          
          xrprod(jl)    = pcvrform(_RI_XYZ__(jl,jrow,jk)) + &
                          lfrac(jl,jk) * pcvsform(_RI_XYZ__(jl,jrow,jk))
          xsprod(jl)    = (1._dp - lfrac(jl,jk)) * pcvsform(_RI_XYZ__(jl,jrow,jk))

          zrlwc(jl,jk)  = xlwc(jl) * zrhoa(jl,jk) 
          zriwc(jl,jk)  = xiwc(jl) * zrhoa(jl,jk) 

          zrainprod(jl,jk) = xrprod(jl) * zrhoa(jl,jk) 
          zsnowprod(jl,jk) = xsprod(jl) * zrhoa(jl,jk) 
          IF (zrlwc(jl,jk) > thres_val) &
            lwc_cv(_RI_XYZ__(jl,jrow,jk)) = zrlwc(jl,jk)/zrhoa(jl,jk)
        END IF
      END DO
    END DO
    
    DO jt=1,nprod
      DO jk=1,nlev
        DO jl=1,kproma
          prod_array(jl,jk,jt,1) = 0._dp
          prod_array(jl,jk,jt,2) = 0._dp
        ENDDO
      ENDDO
    ENDDO

!-------------------------------------------------------------------

    IF (l_trac_cloud) THEN
      DO jk=max_lev_scav,nlev
        DO jt = 1,ntrac 
          zxt_vol(1:kproma,jt) = zxtp1(1:kproma,jk,jt) * &
                                 grvol(_RI_XYZ__(1:kproma,jrow,jk)) * 1.e6_dp
        ENDDO
        lnucl_i(1:kproma,jk) = .FALSE.
        lnucl_l(1:kproma,jk) = .FALSE.

        CALL trac2l(cv_field(js,1)%cloud_field(1:kproma,1:lspec,jk,jrow), &
                    zxt_vol(1:kproma,:), ntrac, js, kproma,&
                    lspec )
        CALL trac2i(cv_field(js,2)%cloud_field(1:kproma,1:ispec,jk,jrow), &
                    zxt_vol(1:kproma,:), ntrac, js, kproma,&
                    ispec )
        CALL cloudtrac2aer( &
                    aero_flx(js)%cloud_field_cv_aer(&
                    _RI_XYZN_(1:kproma,jrow,:,1),jk), &
                    zxt_vol(1:kproma,:), js, kproma, ntrac, 1 )
        CALL cloudtrac2aer( &
                    aero_flx(js)%cloud_field_cv_aer(&
                    _RI_XYZN_(1:kproma,jrow,:,2),jk), &
                    zxt_vol(1:kproma,:), js, kproma, ntrac, 2 )
      END DO
    END IF

!-------------------------------------------------------------------

! call scav_main
! op_ck_20140514+
!!$    CALL SCAV_MAIN(kproma, nlev, ntrac, jrow, js, kspec,                 &
    CALL SCAV_MAIN(status, kproma, nlev, ntrac, jrow, js, kspec,         &
! op_ck_20140514-
                   max_lev_scav, time_step_len,                          &
                   rain(1:kproma,:),          snow(1:kproma,:),          &
                   zrainprod(1:kproma,:),     zsnowprod(1:kproma,:),     &
                   zrlwc(1:kproma,:),         zriwc(1:kproma,:),         &
                   zt(1:kproma,:),            press_3d(_RI_XYZ__(1:kproma,jrow,:)), &
                   zcover(1:kproma,:),        zcover(1:kproma,:),        &
                   grvol(_RI_XYZ__(1:kproma,jrow,:)),    zrhoa(1:kproma,:),         &
                   melt(1:kproma,:),          dummy(1:kproma,:),         &
                   bclwc_cv(_RI_XYZ__(1:kproma,jrow,:)),                 &
                   gboxarea_2d(1:kproma,jrow),zxtp1(1:kproma,:,1:ntrac), &
                   cv_field(js,1)%cloud_field(1:kproma,1:lspec,:,jrow),  &
                   cv_field(js,2)%cloud_field(1:kproma,1:ispec,:,jrow),  &
                   aero_flx(js)%cloud_field_cv_aer(&
                                        _RI_XYZN_(1:kproma,jrow,:,1),:), &
                   aero_flx(js)%cloud_field_cv_aer(&
                                        _RI_XYZN_(1:kproma,jrow,:,2),:), &
                   kpp_rain_cv(_RI_XYZN_(1:kproma,jrow,1:lspec,js)),     &
                   kpp_snow_cv(_RI_XYZN_(1:kproma,jrow,1:ispec,js)),     &
                   aero_flx(js)%aer_flx_cv(_RI_XYZN_(1:kproma,jrow,:,1)),&
                   aero_flx(js)%aer_flx_cv(_RI_XYZN_(1:kproma,jrow,:,2)),&
                   ph(js)%rainpH_cv(_RI_XYZ__(1:kproma,jrow,:)),                   &
                   ph(js)%cloudpH_cv(_RI_XYZ__(1:kproma,jrow,:)),                  &
                   ph(js)%Hp_cloud_cv(_RI_XYZ__(1:kproma,jrow,:)),                 &
                   pseudo_ls, lnucl_i(1:kproma,:), lnucl_l(1:kproma,:),  &
                   KPP_INFO(js)%cv_rain_steps(_RI_XYZ__(1:kproma,jrow,:)),         &
                   KPP_INFO(js)%cv_cloud_steps(_RI_XYZ__(1:kproma,jrow,:)),        &
                   KPP_INFO(js)%cv_rain_rsteps(_RI_XYZ__(1:kproma,jrow,:)),        &
                   KPP_INFO(js)%cv_cloud_rsteps(_RI_XYZ__(1:kproma,jrow,:)),       &
                   prod_array(1:kproma,:,:,:),                           &
!op_kg_20110126+
! um_ak_20110615+
                   phi_i_cv(_RI_XYZ__(1:kproma,jrow,:)),                            &
                   mju_i_cv(_RI_XYZ__(1:kproma,jrow,:)),                            &
                   iwc_cv(_RI_XYZ__(1:kproma,jrow,:)),iwc_T_cv(_RI_XYZ__(1:kproma,jrow,:)) )
!!$                   phi_i_cv(1:kproma,:,jrow),                            &
!!$                   mju_i_cv(1:kproma,:,jrow),                            &
!!$                   iwc_cv(1:kproma,:,jrow),iwc_T_cv(1:kproma,:,jrow) )
! um_ak_20110615-
!op_kg_20110126-
! op_ck_20140514+
    CALL scav_halt(substr, status)
! op_ck_20140514-

!-------------------------------------------------------------------

      IF (l_trac_cloud) THEN
        DO jt=1,lspec_liq
          idx1 = kpp_l_idx(js)%liq_spec(jt,liq2trac)
          idx2 = kpp_l_idx(js)%liq_spec(jt,liq_idx)
          DO jk=1,nlev
            DO jl=1,kproma
              IF (lnucl_l(jl,jk) ) THEN
                zxtp1(jl,jk,idx1) = 0._dp
              ELSE
                IF (pseudo_ls) zxtp1(jl,jk,idx1) = 0._dp
                cv_field(js,1)%cloud_field(jl,idx2,jk,jrow) = 0._dp
              ENDIF
            END DO
          END DO
        END DO
        DO jt=1,ispec_ice
          idx1 = kpp_i_idx(js)%ice_spec(jt,ice2trac)
          idx2 = kpp_i_idx(js)%ice_spec(jt,ice_idx)
          DO jk=1,nlev
            DO jl=1,kproma
              IF (lnucl_i(jl,jk) ) THEN
                zxtp1(jl,jk,idx1) = 0._dp
              ELSE
                IF (pseudo_ls) zxtp1(jl,jk,idx1) = 0._dp
                cv_field(js,2)%cloud_field(jl,idx2,jk,jrow) = 0._dp
              ENDIF
            END DO
          END DO
        END DO
        DO jt=1,aer_count(js)%numbers(c_all)
          idx1 = aer_attr(js)%aer_spec(jt,aer2tracl)
          idx2 = aer_attr(js)%aer_spec(jt,aer2traci)
          DO jk=1,nlev
            DO jl=1,kproma
              IF (lnucl_l(jl,jk) ) THEN
                zxtp1(jl,jk,idx1) = 0._dp
              ELSE 
                IF (pseudo_ls) zxtp1(jl,jk,idx1) = 0._dp
                aero_flx(js)%cloud_field_cv_aer(_RI_XYZN_(jl,jrow,jt,1),jk) = 0._dp
              ENDIF
              IF (lnucl_i(jl,jk) ) THEN
                zxtp1(jl,jk,idx2) = 0._dp
              ELSE 
                IF (pseudo_ls) zxtp1(jl,jk,idx2) = 0._dp
                aero_flx(js)%cloud_field_cv_aer(_RI_XYZN_(jl,jrow,jt,2),jk) = 0._dp
              ENDIF
            END DO
          END DO
        END DO
      ENDIF

!-------------------------------------------------------------------
    
    DO jt=1,nprod
      DO jk=1,nlev
        DO jl=1,kproma
          PR(js)%PROD(JT)%cl_cv(_RI_XYZ__(jl,jrow,jk)) = prod_array(jl,jk,jt,1)
          PR(js)%PROD(JT)%ra_cv(_RI_XYZ__(jl,jrow,jk)) = prod_array(jl,jk,jt,2)
        ENDDO
      ENDDO
    ENDDO



!-------------------------------------------------------
! update tendency

    DO jt=1,ntrac
      IF (.NOT. attr(js)%log_att(jt,lwetdep)) CYCLE
      DO jk=max_lev_scav,nlev
        DO jl=1,kproma

!!$   !          if (zxtp1(jl,jk,jt)*cm(jl,jk) < -1.e-15_dp) &
!!$   !            print*, "WARNING, cv scav negative",js,jl,jk,jt,jrow,&
!!$   !            zxtp1(jl,jk,jt)*cm(jl,jk), pxtp1(jl,jk,jt)

          xtte_scav(jl,jk,jt) = ( zxtp1(jl,jk,jt) * cm(jl,jk) - &
                                  pxtp1(jl,jk,jt) ) / time_step_len
        ENDDO
      ENDDO
    ENDDO

!--------------------------------------------------------
! correct wet deposition fluxes (if necessary)
    IF (column) THEN
      zcov(:) = 0._dp
      DO jk=max_lev_scav,nlev
        DO jl=1,kproma
          zcov(jl) = MAX(zcov(jl),conv_cover(_RI_XYZ__(jl,jrow,jk)))
        ENDDO
      ENDDO
      DO jt=1,lspec
        DO jl=1,kproma
          kpp_rain_cv(_RI_XYZN_(jl,jrow,jt,js)) = &
               kpp_rain_cv(_RI_XYZN_(jl,jrow,jt,js)) * &
                                       zcov(jl)
        ENDDO
      ENDDO
      DO jt=1,ispec
        DO jl=1,kproma
          kpp_snow_cv(_RI_XYZN_(jl,jrow,jt,js)) = &
               kpp_snow_cv(_RI_XYZN_(jl,jrow,jt,js)) * &
                                       zcov(jl)
        ENDDO
      ENDDO
      DO jt=1,aer_count(js)%numbers(c_all)
        DO jl=1,kproma
          aero_flx(js)%aer_flx_cv(_RI_XYZN_(jl,jrow,jt,1)) =         &
            aero_flx(js)%aer_flx_cv(_RI_XYZN_(jl,jrow,jt,1)) * zcov(jl)
          aero_flx(js)%aer_flx_cv(_RI_XYZN_(jl,jrow,jt,2)) =         &
            aero_flx(js)%aer_flx_cv(_RI_XYZN_(jl,jrow,jt,2)) * zcov(jl)
        ENDDO
      ENDDO
    ENDIF

  END SUBROUTINE scav_cv_new
!==============================================================================
! This subroutine collects all the specific parameters for 
! LARGE-SCALE scavenging, calls the scav_main driver, and 
! determines the tracer tendencies from large-scale scavenging

  SUBROUTINE scav_ls_new(js, jrow, kproma, pxtp1)

    ! ECHAM5/MESSy
    USE messy_main_grid_def_mem_bi, ONLY: nlev
    USE messy_main_data_bi,         ONLY: tm1, tte_3d
    ! MESSy
    USE messy_main_timer,          ONLY: time_step_len, delta_time
    USE messy_main_constants_mem,  ONLY: avo => N_A, rd, R_gas
    USE messy_scav_i_kpp,          ONLY: ispec => nspec
    USE messy_scav_l_kpp,          ONLY: lspec => nspec
    USE messy_scav_l_kpp,          ONLY: IND_Hp_l, IND_NH3_l, &
                                         IND_NH4p_l, IND_HNO3_l,     &
                                         IND_NO3m_l, IND_H2SO4_l,    &
                                         IND_HSO4m_l, IND_SO4mm_l
    USE messy_scav_inter,          ONLY: scav_main, trac2l, trac2i,  &
                                         process

    IMPLICIT NONE
    INTEGER,  INTENT(IN) :: js, jrow, kproma
    REAL(dp), INTENT(IN) :: pxtp1(nproma, nlev, ntrac)

    ! local variables
! op_ck_20140514+
    CHARACTER(LEN=*), PARAMETER :: substr = 'scav_ls_new'
    INTEGER  :: status = 0
! op_ck_20140514-
    INTEGER  :: jl, jk, jt, idx1, idx2
    REAL(dp) :: zrlwc(kproma,nlev), zriwc(kproma,nlev)
    REAL(dp) :: zrainprod(kproma,nlev), zsnowprod(kproma,nlev)
    REAL(dp) :: zt(kproma,nlev), zrhoa(kproma,nlev)
    REAL(dp) :: mc(kproma,nlev), cm(kproma,nlev)
    REAL(dp) :: zxtp1(kproma, nlev, 0:ntrac)
    REAL(dp) :: zxt_vol(kproma, ntrac)
    REAL(dp) :: prod_array(kproma,nlev,nprod,2)
    LOGICAL  :: lnucl_i(kproma,nlev), lnucl_l(kproma,nlev)
    INTEGER :: idt, j, jjrow, p ! op_pj_20160926

    ! resetting/initialising values
    zrlwc     = 0._dp
    zriwc     = 0._dp
    zrainprod = 0._dp
    zsnowprod = 0._dp
    ! op_pj_20101209+
    zrhoa     = 0._dp
    zt        = 0._dp
    zxtp1     = 0._dp
    zxt_vol   = 0._dp

    DO p=1, 4
       IF ( ASSOCIATED(set(js)%para(1)%proc(p)%evafrac) ) THEN
          jjrow = jrow
          DO idt = 1, ntrac
             DO j = 1, n_resmod
                process(p)%frac_eva(idt,j)%ptr  => &
                     set(js)%para(1)%proc(p)%evafrac(_RI_XYZN_(:,jjrow,:,idt),j)
                process(p)%frac_eva(idt,j)%ptr(:,:) = 0.0_dp
             END DO
          END DO
       END IF
    END DO

    DO jk=max_lev_scav,nlev
      DO jl=1,kproma
        zt(jl,jk) = tm1(_RI_XYZ__(jl,jrow,jk)) &
             + tte_3d(_RI_XYZ__(jl,jrow,jk)) * time_step_len
        zrhoa(jl,jk) = press_3d(_RI_XYZ__(jl,jrow,jk)) / (rd * zt(jl,jk)) 
        ! convert mixing ratio [mol/mol] to concentration [mcl/cc]
        ! mc = mixing ratio to concentration         
        mc(jl,jk) = (avo/1.e6_dp) * press_3d(_RI_XYZ__(jl,jrow,jk)) / (R_gas * zt(jl,jk))
        ! convert concentration [mcl/cc] to mixing ratio [mol/mol]
        ! cm = concentration to mixing ratio
        cm(jl,jk)= 1._dp/mc(jl,jk)
      ENDDO
    ENDDO

    !=======================================================
    ! mz_ht_20030430+, resetting the cloud and rainwater pH, which is
    !     only being done here whenever the convective scavenging is 
    !     turned-off
    IF (lscav_cv) THEN
      DO jl=1,kproma
        IF (kpp_rain_cv(_RI_XYZN_(jl,jrow,IND_hp_l,js)) <  1.e-5_dp)  &
          ph(js)%rainpH_cv(_RI_XYZ__(jl,jrow,:)) = 19.99_dp
      ENDDO
    ENDIF

    DO jk=max_lev_scav,nlev
      DO jt=1,ntrac
        DO jl=1,kproma
          zxtp1(jl,jk,jt) = pxtp1(jl,jk,jt) * mc(jl,jk)
        ENDDO
      ENDDO
    ENDDO

    IF (lscav_ls) THEN
      ph(js)%cloudpH_ls(_RI_XYZ__(:,jrow,:))  = 19.99_dp
      ph(js)%rainpH_ls(_RI_XYZ__(:,jrow,:))   = 19.99_dp
!      bc_cov(_RI_XYZ__(:,jrow,:))         = 0._dp
      ph(js)%Hp_cloud_ls(_RI_XYZ__(:,jrow,:)) = 0._dp

      kpp_rain_ls(_RI_XYZN_(:,jrow,:,:))  = 0.0_dp
      kpp_snow_ls(_RI_XYZN_(:,jrow,:,:))  = 0.0_dp
      aero_flx(js)%aer_flx_ls(_RI_XYZN_(:,jrow,:,:))   = 0.0_dp

      kpp_info(js)%ls_rain_steps(_RI_XYZ__(:,jrow,:))  = 0._dp
      kpp_info(js)%ls_cloud_steps(_RI_XYZ__(:,jrow,:)) = 0._dp
      kpp_info(js)%ls_rain_rsteps(_RI_XYZ__(:,jrow,:))  = 0._dp
      kpp_info(js)%ls_cloud_rsteps(_RI_XYZ__(:,jrow,:)) = 0._dp
!      bclwc_ls(:,:,jrow)       = 0._dp
    ENDIF

    large_scale_scav: IF (lscav_ls) THEN
      IF (lscav_nuc) THEN
        DO jk=max_lev_scav,nlev 
          DO jl=1,kproma

!         unit conversion of kg / kg(Luft) in kg / m^3(Luft)
            zrlwc(jl,jk) = pmlwc(_RI_XYZ__(jl,jrow,jk)) * zrhoa(jl,jk)
            zriwc(jl,jk) = pmiwc(_RI_XYZ__(jl,jrow,jk)) * zrhoa(jl,jk)
             
!         scaled with cloud cover, since pmratep, pmratesi are grid box 
!         mean values, and not inside the cloud
            IF (pclcover(_RI_XYZ__(jl,jrow,jk)) > 1.e-7_dp ) THEN
              zrainprod(jl,jk) = pmratep(_RI_XYZ__(jl,jrow,jk))  * &
                                 zrhoa(jl,jk) / pclcover(_RI_XYZ__(jl,jrow,jk))
              zsnowprod(jl,jk) = pmratesi(_RI_XYZ__(jl,jrow,jk)) * &
                                 zrhoa(jl,jk) / pclcover(_RI_XYZ__(jl,jrow,jk))
            ELSE
              zrainprod(jl,jk) = 0._dp
              zsnowprod(jl,jk) = 0._dp
            ENDIF
                    
          ENDDO
        END DO
      END IF
      DO jt=1,nprod
        DO jk=1,nlev
          DO jl=1,kproma
            prod_array(jl,jk,jt,1) = 0._dp
            prod_array(jl,jk,jt,2) = 0._dp
          ENDDO
        ENDDO
      ENDDO

!-------------------------------------------------------

      IF (l_trac_cloud) THEN
        DO jk=max_lev_scav,nlev
          DO jt=1,ntrac
            zxt_vol(1:kproma,jt) = MAX(zxtp1(1:kproma,jk,jt),0._dp) * &
                                   grvol(_RI_XYZ__(1:kproma,jrow,jk)) * 1.e6_dp
          ENDDO
          lnucl_i(1:kproma,jk) = .FALSE.
          lnucl_l(1:kproma,jk) = .FALSE.

          CALL trac2l(ls_field(js,1)%cloud_field(1:kproma,1:lspec,jk,jrow), &
                      zxt_vol(1:kproma,:), ntrac, js, kproma,&
                      lspec )

          CALL trac2i(ls_field(js,2)%cloud_field(1:kproma,1:ispec,jk,jrow), &
                      zxt_vol(1:kproma,:), ntrac, js, kproma,&
                      ispec )

          CALL cloudtrac2aer( &
                      aero_flx(js)%cloud_field_ls_aer(_RI_XYZN_(1:kproma,jrow,:,1),jk), &
                      zxt_vol(1:kproma,:), js, kproma, ntrac, 1 )
          CALL cloudtrac2aer( &
                      aero_flx(js)%cloud_field_ls_aer(_RI_XYZN_(1:kproma,jrow,:,2),jk), &
                      zxt_vol(1:kproma,:), js, kproma, ntrac, 2 )
        END DO
      END IF
        
!--------------------------------------
! op_ck_20140514+
!!$      CALL SCAV_MAIN(kproma, nlev, ntrac, jrow, js, kspec,                 &
      CALL SCAV_MAIN(status, kproma, nlev, ntrac, jrow, js, kspec,         &
! op_ck_20140514-
                     max_lev_scav, time_step_len,                          &
                     pfprec(_RI_XYZ__(1:kproma,jrow,:)),   pfsi(_RI_XYZ__(1:kproma,jrow,:)), &
                     zrainprod(1:kproma,:),     zsnowprod(1:kproma,:),     &
                     zrlwc(1:kproma,:),         zriwc(1:kproma,:),         &
                     zt(1:kproma,:),            press_3d(_RI_XYZ__(1:kproma,jrow,:)),  &
                     pclcover(_RI_XYZ__(1:kproma,jrow,:)), prcover(_RI_XYZ__(1:kproma,jrow,:)), &
                     grvol(_RI_XYZ__(1:kproma,jrow,:)),    zrhoa(1:kproma,:),          &
                     pimelt(_RI_XYZ__(1:kproma,jrow,:)),   pisedi(_RI_XYZ__(1:kproma,jrow,:)),  &
                     bclwc_ls(_RI_XYZ__(1:kproma,jrow,:)),                 &
                     gboxarea_2d(1:kproma,jrow),zxtp1(1:kproma,:,1:ntrac), &
                     ls_field(js,1)%cloud_field(1:kproma,1:lspec,:,jrow),  &
                     ls_field(js,2)%cloud_field(1:kproma,1:ispec,:,jrow),  &
                     aero_flx(js)%cloud_field_ls_aer(_RI_XYZN_(1:kproma,jrow,:,1),:), &
                     aero_flx(js)%cloud_field_ls_aer(_RI_XYZN_(1:kproma,jrow,:,2),:), &
                     kpp_rain_ls(_RI_XYZN_(1:kproma,jrow,1:lspec,js)),     &
                     kpp_snow_ls(_RI_XYZN_(1:kproma,jrow,1:ispec,js)),     &
                     aero_flx(js)%aer_flx_ls(_RI_XYZN_(1:kproma,jrow,:,1)),&
                     aero_flx(js)%aer_flx_ls(_RI_XYZN_(1:kproma,jrow,:,2)),&
                     ph(js)%rainpH_ls(_RI_XYZ__(1:kproma,jrow,:)),                   &
                     ph(js)%cloudpH_ls(_RI_XYZ__(1:kproma,jrow,:)),                  &
                     ph(js)%Hp_cloud_ls(_RI_XYZ__(1:kproma,jrow,:)),                 &
                     .TRUE., lnucl_i(1:kproma,:), lnucl_l(1:kproma,:),     &
                     KPP_INFO(js)%ls_rain_steps(_RI_XYZ__(1:kproma,jrow,:)),         &
                     KPP_INFO(js)%ls_cloud_steps(_RI_XYZ__(1:kproma,jrow,:)),        &
                     KPP_INFO(js)%ls_rain_rsteps(_RI_XYZ__(1:kproma,jrow,:)),        &
                     KPP_INFO(js)%ls_cloud_rsteps(_RI_XYZ__(1:kproma,jrow,:)),       &
                     prod_array(1:kproma,:,:,:),                           &
!op_kg_20110126+
! um_ak_20110615+
!!$                     phi_i_ls(1:kproma,:,jrow),                            &
!!$                     mju_i_ls(1:kproma,:,jrow),                            &
!!$                     iwc_ls(1:kproma,:,jrow),iwc_T_ls(1:kproma,:,jrow) )
                     phi_i_ls(_RI_XYZ__(1:kproma,jrow,:)),                            &
                     mju_i_ls(_RI_XYZ__(1:kproma,jrow,:)),                            &
                     iwc_ls(_RI_XYZ__(1:kproma,jrow,:)),iwc_T_ls(_RI_XYZ__(1:kproma,jrow,:)) )
! um_ak_20110615-
!op_kg_20110126- 
! op_ck_20140514+
    CALL scav_halt(substr, status)
! op_ck_20140514-

!-------------------------------------------------------------------

      IF (l_trac_cloud) THEN
        DO jt=1,lspec_liq
          idx1 = kpp_l_idx(js)%liq_spec(jt,liq2trac)
          idx2 = kpp_l_idx(js)%liq_spec(jt,liq_idx)
          DO jk=1,nlev
            DO jl=1,kproma
              zxtp1(jl,jk,idx1) = 0._dp
              IF (.NOT. lnucl_l(jl,jk) ) THEN
                ls_field(js,1)%cloud_field(jl,idx2,jk,jrow) = 0._dp
              ENDIF
            END DO
          END DO
        END DO
        DO jt=1,ispec_ice
          idx1 = kpp_i_idx(js)%ice_spec(jt,ice2trac)
          idx2 = kpp_i_idx(js)%ice_spec(jt,ice_idx)
          DO jk=1,nlev
            DO jl=1,kproma
              zxtp1(jl,jk,idx1) = 0._dp
              IF (.NOT. lnucl_i(jl,jk) ) THEN
                ls_field(js,2)%cloud_field(jl,idx2,jk,jrow) = 0._dp
              ENDIF
            END DO
          END DO
        END DO
        DO jt=1,aer_count(js)%numbers(c_all)
          idx1 = aer_attr(js)%aer_spec(jt,aer2tracl)
          idx2 = aer_attr(js)%aer_spec(jt,aer2traci)
          DO jk=1,nlev
            DO jl=1,kproma
              zxtp1(jl,jk,idx1) = 0._dp
              zxtp1(jl,jk,idx2) = 0._dp
              IF (.NOT. lnucl_l(jl,jk) ) THEN
                aero_flx(js)%cloud_field_ls_aer(_RI_XYZN_(jl,jrow,jt,1),jk) = 0._dp
              ENDIF
              IF (.NOT. lnucl_i(jl,jk) ) THEN
                aero_flx(js)%cloud_field_ls_aer(_RI_XYZN_(jl,jrow,jt,2),jk) = 0._dp
              ENDIF
            END DO
          END DO
        END DO
      ENDIF

!-------------------------------------------------------------------

      DO jt=1,nprod
        DO jk=1,nlev
          DO jl=1,kproma
            PR(js)%PROD(JT)%cl_ls(_RI_XYZ__(jl,jrow,jk)) = prod_array(jl,jk,jt,1)
            PR(js)%PROD(JT)%ra_ls(_RI_XYZ__(jl,jrow,jk)) = prod_array(jl,jk,jt,2)
          ENDDO
        ENDDO
      ENDDO
!------------------------------------------------------------------------
      !   calculating the tracer tendency from the updated concentration, 
      !   recalculated from molecules cm-3 to mixing ratio, 
      !   and the old concentration
        
      DO jt=1,ntrac
        IF (.NOT. attr(js)%log_att(jt,lwetdep)) CYCLE
        DO jk=max_lev_scav,nlev
          DO jl=1,kproma
!!$            if (zxtp1(jl,jk,jt)*cm(jl,jk) < -1.e-15_dp) &
!!$               print*, "WARNING, ls scav negative",js,jl,jk,jt,jrow,&
!!$               zxtp1(jl,jk,jt)*cm(jl,jk), pxtp1(jl,jk,jt)
            xtte_scav(jl,jk,jt) = ( zxtp1(jl,jk,jt) * cm(jl,jk) - &
                                    pxtp1(jl,jk,jt) ) / time_step_len
          ENDDO
        ENDDO
      ENDDO
!--------------------------------------------------------------------
    END IF large_scale_scav
!-------------------------------------------------------------------------------
! sum up wet deposition fluxes (both phases)

      IF (lscav_ls) THEN
        DO jt=1,ispec_ice
          idx1 = kpp_i_idx(js)%ice_spec(jt,ice2liq)
          idx2 = kpp_i_idx(js)%ice_spec(jt,ice_idx)
          kpp_rain_ls(_RI_XYZN_(1:kproma,jrow,idx1,js) ) =    &
            kpp_rain_ls(_RI_XYZN_(1:kproma,jrow,idx1,js)) +   &
            kpp_snow_ls(_RI_XYZN_(1:kproma,jrow,idx2,js))
        ENDDO
        DO jt=1,lspec
          kpp_rain_ls(_RI_XYZN_(1:kproma,jrow,jt,js) ) =    &
            kpp_rain_ls(_RI_XYZN_(1:kproma,jrow,jt,js)) /   &
           (gboxarea_2d(1:kproma,jrow) * time_step_len)
        ENDDO

        IF (lscav_aer) THEN
          DO jt=1, aer_count(js)%numbers(c_all)
            aero_flx(js)%aer_flx_ls(_RI_XYZN_(1:kproma,jrow,jt,1)) =     &
              ( aero_flx(js)%aer_flx_ls(_RI_XYZN_(1:kproma,jrow,jt,1)) + &
                aero_flx(js)%aer_flx_ls(_RI_XYZN_(1:kproma,jrow,jt,2)) ) &
              / (gboxarea_2d(1:kproma,jrow) * time_step_len)
          ENDDO
        ENDIF

        DO jt=1, n_out
          IF (ASSOCIATED(wet_flx(jt,js)%flux_ls_sum) )      &
            wet_flx(jt,js)%flux_ls_sum(1:kproma,jrow) =     &
            wet_flx(jt,js)%flux_ls_sum(1:kproma,jrow) +     &
            kpp_rain_ls(_RI_XYZN_(1:kproma,jrow,kpp_out_nr(jt),js)) * delta_time
        ENDDO
        DO jt=1, n_out_aer
          IF (ASSOCIATED(wet_flx_aer(jt,js)%flux_ls_aer_sum) )          &
            wet_flx_aer(jt,js)%flux_ls_aer_sum(1:kproma,jrow) =         &
            wet_flx_aer(jt,js)%flux_ls_aer_sum(1:kproma,jrow) +         &
            aero_flx(js)%aer_flx_ls(&
            _RI_XYZN_(1:kproma,jrow, wet_flx_aer(jt,js)%aer_out_nr,1))  &
            * delta_time
        ENDDO
      ENDIF
      
      IF (lscav_cv) THEN
        DO jt=1,ispec_ice
          idx1 = kpp_i_idx(js)%ice_spec(jt,ice2liq)
          idx2 = kpp_i_idx(js)%ice_spec(jt,ice_idx)
          kpp_rain_cv(_RI_XYZN_(1:kproma,jrow,idx1,js)) =      &
            kpp_rain_cv(_RI_XYZN_(1:kproma,jrow,idx1,js)) +   &
            kpp_snow_cv(_RI_XYZN_(1:kproma,jrow,idx2,js) )  
        ENDDO
        DO jt=1,lspec
          kpp_rain_cv(_RI_XYZN_(1:kproma,jrow,jt,js) ) =    &
            kpp_rain_cv(_RI_XYZN_(1:kproma,jrow,jt,js)) /   &
           (gboxarea_2d(1:kproma,jrow) * time_step_len)
        ENDDO

        IF (lscav_aer) THEN
          DO jt=1, aer_count(js)%numbers(c_all)
            aero_flx(js)%aer_flx_cv(_RI_XYZN_(1:kproma,jrow,jt,1)) =     &
              ( aero_flx(js)%aer_flx_cv(_RI_XYZN_(1:kproma,jrow,jt,1)) + &
                aero_flx(js)%aer_flx_cv(_RI_XYZN_(1:kproma,jrow,jt,2)) ) &
              / (gboxarea_2d(1:kproma, jrow) * time_step_len)
          ENDDO
        ENDIF
        DO jt=1, n_out
          IF (ASSOCIATED(wet_flx(jt,js)%flux_cv_sum) )                 &
            wet_flx(jt,js)%flux_cv_sum(1:kproma,jrow) =                &
            wet_flx(jt,js)%flux_cv_sum(1:kproma,jrow) +                &
            kpp_rain_cv(_RI_XYZN_(1:kproma,jrow,kpp_out_nr(jt),js)) * delta_time
        ENDDO
        DO jt=1, n_out_aer
          IF (ASSOCIATED(wet_flx_aer(jt,js)%flux_cv_aer_sum) )         &
            wet_flx_aer(jt,js)%flux_cv_aer_sum(1:kproma,jrow) =        &
            wet_flx_aer(jt,js)%flux_cv_aer_sum(1:kproma,jrow) +        &
            aero_flx(js)%aer_flx_cv(&
            _RI_XYZN_(1:kproma,jrow,wet_flx_aer(jt,js)%aer_out_nr,1))  &
            * delta_time
        ENDDO
      ENDIF
      
      special_flux(js)%nitrate_flx(1:kproma,jrow) = 0.0_dp
      special_flux(js)%sulfate_flx(1:kproma,jrow) = 0.0_dp
      special_flux(js)%ammoni_flx(1:kproma,jrow)  = 0.0_dp
      IF (lscav_ls) THEN
        special_flux(js)%nitrate_flx(1:kproma,jrow) =    &
          kpp_rain_ls(_RI_XYZN_(1:kproma,jrow, IND_HNO3_l, js))  +  &
          kpp_rain_ls(_RI_XYZN_(1:kproma,jrow, IND_NO3m_l, js))
        special_flux(js)%sulfate_flx(1:kproma,jrow) =    &
          kpp_rain_ls(_RI_XYZN_(1:kproma,jrow, IND_HSO4m_l, js)) +  &
          kpp_rain_ls(_RI_XYZN_(1:kproma,jrow, IND_SO4mm_l, js)) +  &
          kpp_rain_ls(_RI_XYZN_(1:kproma,jrow, IND_H2SO4_l, js))
        special_flux(js)%ammoni_flx(1:kproma,jrow)  =    &
          kpp_rain_ls(_RI_XYZN_(1:kproma,jrow, IND_NH3_l, js))   +  &
          kpp_rain_ls(_RI_XYZN_(1:kproma,jrow, IND_NH4p_l, js))
      ENDIF
      IF (lscav_cv) THEN
        special_flux(js)%nitrate_flx(1:kproma,jrow) =    &
          special_flux(js)%nitrate_flx(1:kproma,jrow) +  & 
          kpp_rain_cv(_RI_XYZN_(1:kproma,jrow, IND_HNO3_l, js))  +  &
          kpp_rain_cv(_RI_XYZN_(1:kproma,jrow, IND_NO3m_l, js))
        special_flux(js)%sulfate_flx(1:kproma,jrow) =    &
        special_flux(js)%sulfate_flx(1:kproma,jrow)   +  &
        kpp_rain_cv(_RI_XYZN_(1:kproma,jrow, IND_HSO4m_l, js))   +  &
        kpp_rain_cv(_RI_XYZN_(1:kproma,jrow, IND_SO4mm_l, js))   +  &
        kpp_rain_cv(_RI_XYZN_(1:kproma,jrow, IND_H2SO4_l, js))
        special_flux(js)%ammoni_flx(1:kproma,jrow)  =    &
          special_flux(js)%ammoni_flx(1:kproma,jrow)  +  &
          kpp_rain_cv(_RI_XYZN_(1:kproma,jrow, IND_NH3_l, js))   +  &
          kpp_rain_cv(_RI_XYZN_(1:kproma,jrow, IND_NH4p_l, js))
      ENDIF
      special_flux(js)%nitrate_flx_sum(1:kproma,jrow) =    &
        special_flux(js)%nitrate_flx_sum(1:kproma,jrow) +  &
        special_flux(js)%nitrate_flx(1:kproma,jrow) * delta_time
      special_flux(js)%sulfate_flx_sum(1:kproma,jrow) =    &
        special_flux(js)%sulfate_flx_sum(1:kproma,jrow) +  &
        special_flux(js)%sulfate_flx(1:kproma,jrow) * delta_time
      special_flux(js)%ammoni_flx_sum(1:kproma,jrow)  =    &
        special_flux(js)%ammoni_flx_sum(1:kproma,jrow)  +  &
        special_flux(js)%ammoni_flx(1:kproma,jrow)  * delta_time

      RETURN
    END SUBROUTINE scav_ls_new

!==============================================================================
! This subroutine prepares the cloud fields for evaporation after all chemical 
! transformations in the gas phase
! to be called after all other chemical/aerosol processes

  SUBROUTINE EVAPO_CLOUD(js, jrow, kproma, phase, spec, ntrac,   &
                         cloud_field_aer, cloud_field_kpp)

    USE messy_main_grid_def_mem_bi, ONLY: nlev
    USE messy_main_data_bi,         ONLY: tte_3d, tm1
    USE messy_main_timer,           ONLY: time_step_len
    USE messy_main_constants_mem,   ONLY: avo => N_A, R_gas
    USE MESSY_SCAV_INTER,           ONLY: process, made3

    INTEGER,  INTENT(IN)    :: js, jrow, kproma, phase, spec, ntrac
    REAL(dp), INTENT(INOUT) :: cloud_field_kpp(kproma,spec,nlev)
    REAL(dp), INTENT(INOUT) :: &
      cloud_field_aer(kproma,aer_count(js)%numbers(c_all),nlev)

    INTEGER  :: jl, lproma, jt, jk, i
    INTEGER  :: lwork(kproma)
    REAL(dp) :: laerfld(aer_count(js)%numbers(c_all), kproma)
    REAL(dp) :: lkppfld(kspec,kproma)
    REAL(dp) :: xt_vec(ntrac,kproma)
    REAL(dp) :: zxtp1(kproma,nlev, ntrac)
    REAL(dp) :: zt(kproma,nlev), mc(kproma,nlev), cm(kproma, nlev)
    REAL(dp) :: evfrac(kproma)
    INTEGER  :: idt, j
! op_ck_20140523+
    REAL(dp) :: evfrac_aer(kproma,ntrac,n_resmod)
    INTRINSIC :: REAL
! op_ck_20140523-

    DO jk=max_lev_scav,nlev 
      
! packing for clouds
      lproma   = 0
      lwork(:) = 0

      xt_vec(:,:)  = 0._dp
      laerfld(:,:) = 0._dp
      lkppfld(:,:) = 0._dp

      DO jl=1,kproma
        IF ( (MAXVAL(cloud_field_aer(jl,:,jk)) > 1.e-34_dp ) .OR.        &
             (MAXVAL(cloud_field_kpp(jl,1:spec,jk)) > 1.e-34_dp ) ) THEN
          lproma        = lproma + 1
          lwork(lproma) = jl
        ENDIF
      ENDDO
      
!!$   IF (lproma > 0 ) THEN     ! op_pj_20170113

        DO jl=1,kproma
          zt(jl,jk)    = tm1(_RI_XYZ__(jl,jrow,jk)) + tte_3d(_RI_XYZ__(jl,jrow,jk)) * time_step_len
        ! convert mixing ratio [mol/mol] to concentration [mcl/cc]
        ! mc = mixing ratio to concentration
          mc(jl,jk) = (avo/1.e6_dp) * press_3d(_RI_XYZ__(jl,jrow,jk)) / (R_gas * zt(jl,jk))
        ! convert concentration [mcl/cc] to mixing ratio [mol/mol]
        ! cm = concentration to mixing ratio
          cm(jl,jk)= 1._dp/mc(jl,jk)
        ENDDO
        zxtp1(:,jk,:) = 0._dp

      IF (lproma > 0 ) THEN  ! op_pj_20170113
    
        DO jt=1,ntrac
          DO jl=1,lproma
            i = lwork(jl)
            xt_vec(jt,jl) = zxtp1(i,jk,jt)
          ENDDO
        ENDDO
        DO jt=1,aer_count(js)%numbers(c_all)
          DO jl=1,lproma
            i = lwork(jl)
            laerfld(jt,jl) = cloud_field_aer(i,jt,jk) &
                           / (grvol(_RI_XYZ__(i,jrow,jk)) * 1.e6_dp)
          ENDDO
        ENDDO
        DO jt=1,spec
          DO jl=1,lproma
            i = lwork(jl)
            lkppfld(jt,jl) = cloud_field_kpp(i,jt,jk) &
                           / (grvol(_RI_XYZ__(i,jrow,jk)) * 1.e6_dp)
          ENDDO
        ENDDO
        DO jl=1,lproma
          i = lwork(jl)
! op_ck_20140523+
!!$          evfrac(jl) = process(3)%frac_eva(i,jk)
! op_pj_20160926+
!!$          evfrac(jl) = SUM(process(3)%frac_eva(i,jk,:,n_resmod)) &
!!$               / REAL(ntrac, kind=dp)
!!$          evfrac_aer(jl,:,:) = process(3)%frac_eva(i,jk,:,:)
          evfrac(jl) = 0.0_dp
          DO idt = 1, ntrac
             evfrac(jl) = evfrac(jl) + &
                  process(3)%frac_eva(idt,n_resmod)%ptr(i,jk) 
          END DO
          evfrac(jl) = evfrac(jl) / REAL(ntrac, kind=dp)
          DO idt = 1, ntrac
             DO j=1, n_resmod
                evfrac_aer(jl,idt,j) = process(3)%frac_eva(idt,j)%ptr(i,jk)
             END DO
          END DO
! op_pj_20160926-
! op_ck_20140523-
        ENDDO
! op_ck_20140523+
        IF (ALLOCATED(made3)) evfrac(1:lproma) = 1._dp
! op_ck_20140523-

        IF (l_trac_cloud) THEN
          IF (MAXVAL(laerfld(:,:)) > 0._dp)                       &
            CALL aer2cloudtrac( laerfld(:,1:lproma),              &
                          xt_vec(1:ntrac,1:lproma),               &
                          js, lproma, ntrac, phase)
          IF (MAXVAL(lkppfld(1:spec,:)) > 0._dp)                  &
            CALL cloud2trac(lkppfld(1:spec,1:lproma),             &
                       xt_vec(1:ntrac,1:lproma),                  &
                       ntrac, lproma, js, phase, spec)
        END IF
! call evaporation routine
        IF (MAXVAL(laerfld(:,:)) > 0._dp) THEN
          CALL evapo_aer(xt_vec(:,1:lproma), laerfld(:,1:lproma), &
! op_ck_20140523+
!!$                   SPREAD(1.e-6_dp,1,lproma), evfrac(1:lproma),   &
                   SPREAD(1.e-6_dp,1,lproma), evfrac_aer(1:lproma,:,:), &
! op_ck_20140523-
                   ntrac, lproma, lwork(1:lproma), js, jk, jrow)
        END IF

        IF (MAXVAL(lkppfld(1:spec,:)) > 0._dp) THEN
          CALL evapo(lkppfld(1:spec,1:lproma), xt_vec(:,1:lproma), &
                   SPREAD(1.e-6_dp,1,lproma), evfrac(1:lproma),    &
                   ntrac, lproma,       &
                   lwork(1:lproma), js, jk, jrow, phase, spec)
        ENDIF

! unpack
        DO jt=1,ntrac
          DO jl=1,lproma
            i = lwork(jl)
            zxtp1(i,jk,jt) = xt_vec(jt,jl)
          ENDDO
        ENDDO
        
      END IF        ! on lproma
    ENDDO           ! vertical level loop
          
    DO jt=1,ntrac
      IF (.NOT. attr(js)%log_att(jt,lwetdep)) CYCLE
      DO jk=max_lev_scav,nlev 
        DO jl=1,kproma
          xtte_scav(jl,jk,jt) = xtte_scav(jl,jk,jt) + &
                              ( zxtp1(jl,jk,jt) * cm(jl,jk) )/time_step_len
        ENDDO
      ENDDO
    END DO

  END SUBROUTINE EVAPO_CLOUD

! op_ck_20140514+
!==============================================================================
! This subroutine handles potential errors in scav_main.

  SUBROUTINE scav_halt(substr, status)

    ! MESSy
    USE messy_main_blather_bi, ONLY: error_bi

    ! I/O
    CHARACTER(LEN=*), INTENT(IN)  :: substr
    INTEGER,          INTENT(IN)  :: status

    IF (status == 0) RETURN

    SELECT CASE (status)
    CASE (11)
      CALL error_bi('Could not determine mode of aerosol number tracer !', &
           'messy_scav_inter::calc_evapo_frac')
    CASE (12)
      CALL error_bi('Could not determine mode of aerosol tracer !', &
           'messy_scav_inter::calc_evapo_frac')
    CASE DEFAULT
      CALL error_bi('Unspecified error', 'messy_scav_si::'//substr)
    END SELECT

  END SUBROUTINE scav_halt
! op_ck_20140514-
!==============================================================================
END MODULE messy_scav_si
